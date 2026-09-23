# Amino-acid substitution model inference in IQ-TREE 3

**Audience:** a reader who knows phylogenetics and maximum likelihood and needs the exact method
IQ-TREE 3 uses to infer amino-acid substitution matrices, time-reversible (Q) and
non-time-reversible (nQ), stated as mathematics together with the code that implements it.

**Baseline:** commit `63c330d9`, upstream tag `v3.1.4`, read on 2026-09-15. Every formula,
constant, and default below was read out of the source in that tree. Line numbers drift, symbol
names do not, so grep the name if an anchor misses.

**Companions:** `docs/agent/ARCHITECTURE.md` (how the codebase is organised) and
`docs/agent/FILE_INDEX.md` (which file to open when). This document is about the method, not the
plumbing. Section 16 records which claims were verified directly and which were not.

**What this document is not.** It describes how IQ-TREE 3 infers amino-acid models today, from
the mathematical, optimization, and algorithmic side; it is one of the three reference
documents on the current code, with `ARCHITECTURE.md` (programming architecture) and
`FILE_INDEX.md` (navigation). It is not a plan, a design, or a decision record for this project,
and nothing in it is an instruction about what the project should build. Where it mentions a
π-constrained model, it does so only for comparison with what exists. The project's design is in
`docs/agent/design/`, its decisions in `docs/agent/DECISIONS.md`, and its plan in
`docs/agent/PLAN.md`.

---

## 0. Scope, and what "maximum likelihood" leaves unsaid

"IQ-TREE estimates amino-acid models by maximum likelihood" is true and almost content-free.
The content lies in seven separate choices, and the amino-acid model families IQ-TREE ships make
them differently:

1. what is a free parameter and what is fixed by construction,
2. how the rate matrix is parameterised, and therefore what the feasible set is,
3. how overall scale is fixed, since a rate matrix and a branch length are confounded,
4. where the stationary distribution comes from: supplied, counted, derived, or optimised,
5. what likelihood is maximised: one alignment, or a sum over many,
6. what optimiser runs, and in particular whether its gradients are analytic or numerical,
7. what happens to branch lengths, rate heterogeneity, the root, and the topology while the
   matrix is being estimated.

This document answers those seven questions for every amino-acid model family in IQ-TREE 3.

### Boundaries

In scope: everything IQ-TREE applies to `SEQ_PROTEIN` data. That is 28 fixed empirical
reversible matrices (LG, WAG, JTT, Q.pfam and the rest), the fully parameterised reversible
`GTR20`, six fixed empirical non-reversible matrices (`NQ.*`), the fully parameterised
non-reversible `NONREV`, the generic `UNREST`, profile mixtures (C10 to C60, CF4, EX2, EX3, EHO,
UL2, UL3, EX_EHO), full-matrix mixtures fused with a rate model (LG4M, LG4X), the
linked-exchangeability GTRpmix procedure, the posterior-mean site-frequency (PMSF) procedure,
and the site-specific mutation-selection models new in v3.1.4.

Out of scope: DNA, codon, binary, morphological, and PoMo models, except where shared machinery
changes the amino-acid answer. The `cmaple/` directory is a separate vendored engine carrying
its own duplicate `NONREV` and `GTR20`; it is not the model layer IQ-TREE's ML search uses and
nothing in it is described here.

Rate heterogeneity across sites (`+I`, `+G`, `+R`, `+H`) is not a substitution model, but it
cannot be omitted, because the matrix and the rate model are estimated inside the same
alternating loop and each changes the other's estimate. Section 2.4 gives it exactly the depth
that dependence requires.

### Reading order

Sections 1 and 2 define the object and the objective. Section 3 is the catalogue. Section 4 is
the shared numerical machinery, and it is where the surprises are. Sections 5 to 11 are the
seven estimation procedures, each answering the seven questions above. Sections 12 to 14 cover
model selection, parameter counting, and output formats. Section 15 lists defects and dead code,
several of which silently change results.

---

## 1. The object being estimated

### 1.1 The continuous-time Markov chain

Amino-acid substitution at a site is modelled as a homogeneous continuous-time Markov chain on
20 states, indexed in IQ-TREE's order (`alignment/alignment.cpp:34`, `symbols_protein`):

```
index: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
state: A R N D C Q E G H I  L  K  M  F  P  S  T  W  Y  V
```

The chain is specified by an instantaneous rate matrix Q with non-negative off-diagonals and
zero row sums,

```
Q_ij >= 0  for i != j,      Q_ii = - sum_{j != i} Q_ij,
```

and transition probabilities over a branch of length t are given by the matrix exponential

```
P(t) = exp(Q t),     P_ij(t) = Pr(state j at the far end | state i at the near end).
```

A stationary distribution satisfies `pi^T Q = 0` with `sum_i pi_i = 1`, and is unique whenever Q
is irreducible, which holds for any Q with strictly positive off-diagonals.

The chain is **time-reversible** if detailed balance holds:

```
pi_i Q_ij = pi_j Q_ji     for all i, j.
```

Reversibility is not a technical convenience. It is what makes the likelihood of an unrooted
tree well defined. Dropping it changes the tree object, the likelihood kernel, the parameter
count, the eigenproblem, the meaning of pi, and the cost. Every difference between the two paths
in IQ-TREE traces back to this one equation.

### 1.2 The reversible parameterisation

A reversible Q is written as a symmetric exchangeability matrix R times a diagonal frequency
matrix:

```
Q_ij = R_ij * pi_j     (i != j),     R_ij = R_ji >= 0.
```

Detailed balance is then automatic, since `pi_i R_ij pi_j = pi_j R_ji pi_i`. Any strictly
positive pi can be paired with any symmetric R and remain stationary, which is exactly what
makes profile mixture models possible (Section 9).

In code, `rates[]` holds the **upper triangle of R**, row-major, of length `n(n-1)/2 = 190`:

```
k     :  0      1      2    ...  18      19    ...  189
pair  : (0,1)  (0,2)  (0,3)     (0,19)  (1,2)      (18,19)
```

Two code paths build Q from R, and they agree exactly.
`EigenDecomposition::computeRateMatrix` (`utils/eigendecomposition.cpp:478-518`) and the inline
Eigen3 version (`model/modelmarkov.cpp:1502-1514`) both compute

```
Q_ij  <- R_ij * pi_j                       (i != j)
m_i   <- sum_{j != i} Q_ij
mu    <- sum_i pi_i * m_i                  (expected substitution rate)
delta <- total_num_subst / mu
Q     <- delta * (Q - diag(m))
```

The `pi_j` multiplication is guarded by `!ignore_state_freq`, and `ignore_state_freq` is true
exactly for non-reversible models.

### 1.3 The non-reversible parameterisation

A non-reversible Q is stored directly, with no exchangeability matrix and no `diag(pi)` factor.
`rates[]` holds **all off-diagonal entries of Q, row-major, skipping the diagonal**, of length
`n(n-1) = 380`:

```
k     :  0      1    ...  18      19     20   ...  379
entry : (0,1)  (0,2)     (0,19)  (1,0)  (1,2)     (19,18)
```

The diagonal is derived and never stored (`model/modelmarkov.cpp:1256-1263`):

```cpp
for (i = 0, k = 0; i < num_states; i++) {
    double *rate_row = rate_matrix+(i*num_states);
    double row_sum = 0.0;
    for (j = 0; j < num_states; j++)
        if (j != i) row_sum += (rate_row[j] = rates[k++]);
    rate_row[i] = -row_sum;
}
```

**Switching between the two.** `ModelMarkov::setReversible(false, adapt_tree)`
(`model/modelmarkov.cpp:75-152`) expands a reversible half matrix into a full one through
detailed balance, `rate_matrix[i][j] = rates[k] * pi_j` and
`rate_matrix[j][i] = rates[k] * pi_i`, repacks into the 380-entry layout, sets
`ignore_state_freq = true`, allocates the complex eigen arrays, and **converts the tree from
unrooted to rooted**. That last step is why `-m NONREV` yields a rooted tree with no request
from the user.

### 1.4 Identifiability, and how scale is fixed

Q and branch length are confounded, since `exp(Qt)` is invariant under `Q -> cQ, t -> t/c`.
IQ-TREE removes the redundancy twice, once by normalisation and once by pinning a parameter.

**Normalisation.** Every call to `decomposeRateMatrix` rescales Q so that the expected number of
substitutions per unit branch length equals `total_num_subst`, which is 1.0 for every
non-mixture model (`utils/eigendecomposition.cpp:32`):

```
sum_i pi_i * ( -Q_ii ) = total_num_subst = 1.
```

Branch lengths are therefore expected substitutions per site.

**Pinning.** Because normalisation already fixes the scale, one entry of `rates[]` is redundant.
`num_params = getNumRateEntries() - 1` (`model/modelprotein.cpp:1177`, `1233`), and
`ModelMarkov::setVariables` copies only the first `num_params` entries into the optimiser vector
(`model/modelmarkov.cpp:1018-1027`). The pinned entry keeps its seed value:

| model | free rates | pinned entry | identity |
|---|---|---|---|
| `GTR20` | 189 of 190 | `rates[189]` | the Tyr-Val exchangeability |
| `NONREV` | 379 of 380 | `rates[379]` | the directed rate Val to Tyr |

This is an unusual but harmless convention: the more common choice is to pin one rate to 1.0,
whereas IQ-TREE pins it to whatever the seed matrix carried after `rescaleRates` scaled the seed
so its largest entry is 10.0 (`model/modelprotein.cpp:1090-1105`, `AA_SCALE = 10.0`). Since the
normalisation rescales everything afterwards, the seed value only chooses a point on the scale
ray, not the fitted model.

**Frequencies are parameterised projectively too.** When pi is free (`+FO`, `FREQ_ESTIMATE`)
only `state_freq[0..18]` enter the optimiser, `state_freq[19]` is held fixed, the vector may have
any positive sum during the search, and `ModelMarkov::getStateFrequency` normalises on every
read (`model/modelmarkov.cpp:875-893`). After the optimiser returns, `scaleStateFreq(true)`
renormalises in place (`model/modelmarkov.cpp:1226-1229`). The comment in the source, dated
2015-09-07, says so explicitly: "relax the sum of state_freq to be 1, this will be done at the
end of optimization".

### 1.5 Where pi comes from, the single most important asymmetry

**Reversible models: pi is an input.** `ModelMarkov::init_state_freq`
(`model/modelmarkov.cpp:299-351`) fills `state_freq` according to the `StateFreqType`
(`utils/tools.h:469-482`):

| suffix | `StateFreqType` | pi |
|---|---|---|
| none, on a named matrix | `FREQ_USER_DEFINED` | the vector shipped with the matrix in its NEXUS block |
| `+F`, `+FC` | `FREQ_EMPIRICAL` | counted from the alignment, then held fixed |
| `+FO`, `+Fo` | `FREQ_ESTIMATE` | 19 extra free parameters, optimised jointly with the rates |
| `+FQ` | `FREQ_EQUAL` | 1/20 each |
| `+F{p1,...,p20}` | `FREQ_USER_DEFINED` | the literal vector |

The protein default is decided inside `ModelProtein::init`, not in `ModelFactory`
(`model/modelfactory.cpp:439`, `case SEQ_PROTEIN: break; // let ModelProtein decide by itself`):
a named empirical matrix defaults to its own pi, `GTR20` defaults to `FREQ_EMPIRICAL`, and
`NONREV` defaults to `FREQ_ESTIMATE`.

**Non-reversible models: pi is an output.** It is re-solved from Q on *every likelihood
evaluation*, inside `decomposeRateMatrixNonrev` (`model/modelmarkov.cpp:1247-1267`):

```cpp
if (freq_type != FREQ_USER_DEFINED || Params::getInstance().optimize_from_given_params)
    for (i = 0; i < num_states; i++) state_freq[i] = 1.0/num_states;   // reset to uniform
... unpack rates[] into rate_matrix, derive the diagonal ...
if (freq_type != FREQ_USER_DEFINED || Params::getInstance().optimize_from_given_params)
    computeStateFreqFromQMatrix(rate_matrix, state_freq, num_states);  // solve pi^T Q = 0
```

The solver (`model/modelmarkov.cpp:2118-2131`) stacks the normalisation onto the balance
equations and solves the resulting overdetermined system in least squares:

```cpp
int computeStateFreqFromQMatrix (double Q[], double pi[], int n) {
    MatrixXd A(n+1, n);
    A.topRows(1).setOnes();                       // 1^T pi = 1
    A.bottomRows(n) = Map<MatrixXd>(Q, n, n);     // column-major Map transposes: Q^T pi = 0
    VectorXd b(n+1); b.setZero(); b(0) = 1.0;
    Map<VectorXd> freq(pi, n);
    freq = A.colPivHouseholderQr().solve(b);
    ASSERT(fabs(freq.sum()-1.0) < 1e-4);
    return 0;
}
```

In mathematics, solve

```
[ 1^T ]        [ 1 ]
[ Q^T ] pi  =  [ 0 ]
```

by column-pivoted Householder QR: 21 equations in 20 unknowns. The `Map<MatrixXd>` is
column-major while `rate_matrix` is row-major, so the map reads Q transposed, which is what
`pi^T Q = 0` requires. The system has rank 20 because the 20 balance equations have rank 19
(their sum vanishes identically, since `sum_j sum_i pi_i Q_ij = sum_i pi_i * 0 = 0`) and the
normalisation supplies the missing equation.

**Three consequences that regularly surprise people.**

First, `+F` does nothing on a non-reversible amino-acid model. `-m NQ.pfam+F` sets
`freq_type = FREQ_EMPIRICAL`, `init_state_freq` writes the empirical composition into
`state_freq`, and the very next `decomposeRateMatrixNonrev` overwrites it with uniform and
re-solves from Q. The empirical vector never reaches the likelihood.

Second, `-m NONREV+F{...}` is the one case where a supplied pi survives, and it survives
unsoundly. A literal frequency vector sets `freq_type = FREQ_USER_DEFINED`, which makes both
guards false, so the solve is skipped and `state_freq` stays at the user's vector while Q is
estimated with no constraint tying the two together. The root distribution used in the
likelihood is then **not** the stationary distribution of the fitted Q. This path is live and
should not be changed without a decision entry.

Third, for `-m NQ.pfam` the pi printed in the NEXUS block is read, immediately discarded, and
replaced by the solved stationary vector. `readParameters` compares the two and prints a warning
if any entry differs by more than 1e-3 (`model/modelmarkov.cpp:1833-1842`). That warning is the
only consistency check in the codebase between a supplied pi and the Q it is supposed to belong
to.

---

## 2. The objective: what likelihood is maximised

### 2.1 Data reduction and the top-level formula

Identical alignment columns are collapsed into patterns with integer multiplicities. The tree
log-likelihood is a multiplicity-weighted sum over patterns
(`tree/phylotreesse.cpp:535-542`, `tree/phylokernelnew.h:2966-2976`):

```
log L = sum over patterns p of  n_p * [ log( L_p + iota_p ) + s_p * log(2^-256) ]
```

where `n_p` is the pattern count, `iota_p` the invariable-site term (Section 2.4), `L_p` the
mixture sum over rate and model classes, and `s_p` an integer count of rescalings applied during
pruning to avoid underflow. The constant is
`LOG_SCALING_THRESHOLD = -177.4456782233459932741 = -256 log 2` (`tree/phylotree.h:74-80`).
Log-space accumulation happens only at the top of the recursion, once per pattern.

The likelihood is always evaluated **at a branch**, not at a node
(`PhyloTree::computeLikelihood`, `tree/phylotree.cpp:1274-1315`). The branch is picked as the
pendant edge of the farthest leaf and cached in `current_it`.

### 2.2 Reversible: Felsenstein pruning in the eigenvector basis

The mathematical content is standard Felsenstein pruning. The implementation is not, and the
difference matters for anyone reading the kernel: **pi never appears explicitly in the reversible
recursion.** It is absorbed into the eigenvectors.

With `S = Pi^(1/2) Q Pi^(-1/2)` symmetric and `S = V Lambda V^T` its orthonormal
decomposition, IQ-TREE stores (`model/modelmarkov.cpp:1519-1562`)

```
U      = Pi^(-1/2) V
U^(-1) = V^T Pi^(1/2)
Q      = U Lambda U^(-1),      P(t) = U exp(Lambda t) U^(-1)
```

and the identity that both code paths satisfy exactly is `(U^(-1))_{ba} = pi_a U_{ab}`. Hence

```
sum_a pi_a f1_a (P(t) f2)_a  =  sum_b exp(lambda_b t) (U^(-1) f1)_b (U^(-1) f2)_b,
```

and the right-hand side is literally what the code computes. Tip vectors are pre-multiplied by
`U^(-1)` once, for every state including ambiguity codes
(`tree/phylotreesse.cpp:364-372`, `multiplyWithInvEigenvector`), and the recursion for an
internal node v with children L and R, for rate category c, is

```
x_v = U^(-1) [ ( U exp(Lambda r_c t_L) x_L )  .*  ( U exp(Lambda r_c t_R) x_R ) ]
```

which in the state basis (`f = U x`) is exactly

```
f_v(a) = [ sum_b P_ab(r_c t_L) f_L(b) ] * [ sum_b P_ab(r_c t_R) f_R(b) ].
```

Multifurcations multiply over all children. The complete reversible formula is

```
log L = sum_p n_p log [ p_inv * iota'_p
                        + sum_m w_m sum_c rho_c sum_a pi_a f^(m,c)_L(a) sum_b P^(m)_ab(r_c t) f^(m,c)_R(b) ]
```

with `t` the length of the branch carrying the virtual root, `rho_c = site_rate->getProp(c)`,
and `w_m = model->getMixtureWeight(m)`.

**The pulley principle.** Reversibility makes this independent of which branch carries the
virtual root. IQ-TREE exploits it in two visible places: `computeLikelihood` and
`computePatternLhCat` pick different branches without changing the value, and a rooted input
tree is silently unrooted when the model is reversible
(`main/treetesting.cpp:1210-1214`).

### 2.3 Non-reversible: a rooted tree with pi as the root distribution

**The tree is rooted.** `PhyloTree::convertToRooted` (`tree/phylotree.cpp:5906-5960`) adds an
extra degree-one leaf named `ROOT_NAME` with `id == leafNum`, plus an internal node, and places
it either at the outgroup taxon given by `-o` or at the **midpoint of the longest path**. Kernels
test for it with `MTree::isRootLeaf` (`tree/mtree.h:135-137`).

**`state_freq` is literally the root distribution.** In the non-reversible kernel the root leaf's
partial-likelihood vector is filled from `model->getStateFrequency(...)` and is **not** multiplied
by any transition matrix (`tree/phylokernelnonrev.h:1149-1158`,
`tree/phylokernelnew.h:960-965`):

```cpp
if (isRootLeaf(dad)) {
    for (size_t c = 0; c < ncat_mix; c++) {
        double *lh_node = partial_lh_node + (c*nstates);
        size_t m = c/denom;
        model->getStateFrequency(lh_node, m);
        double prop = site_rate->getProp(c%ncat) * model->getMixtureWeight(m);
        for (size_t i = 0; i < nstates; i++) lh_node[i] *= prop;
    }
}
```

Consistently, the root's pendant branch length is ignored: `optimizeOneBranch` returns
immediately for a branch incident on the root (`tree/phylotree.cpp:2635-2637`) and
`getNBranchParameters` subtracts one degree of freedom for it (`tree/phylotree.cpp:2614-2626`).
So the rooted non-reversible site likelihood is

```
L_p^(c) = rho_c w_m sum_a f_u(a) sum_b P^(m)_ab(r_c t) f_v(b),   f_rootleaf(a) = pi_a.
```

Since pi is the solved stationary distribution of Q, the root distribution is stationary by
construction under the default path, and only under it. Under `-m NONREV+F{...}` it is whatever
the user supplied.

**No eigen-space trick.** Non-reversible tip vectors are *not* transformed by `U^(-1)`
(`tree/phylotreesse.cpp:354-355`, `369-372`), so all arithmetic runs in the state basis with
explicit `computeTransMatrix` calls. The transition matrix for a child pointing toward the root
is **transposed** before use (`tree/phylokernelnew.h:936-948`), which is how directionality is
implemented. This is the second genuine difference: reversible means eigen-space with pi baked
into `U^(-1)`, non-reversible means explicit P(t) matrices with pi applied once, at the root.

**Root position is estimated.** `PhyloTree::optimizeRootPosition`
(`tree/phylotree.cpp:3072-3125`) enumerates every branch within `root_move_dist` (default 2,
`utils/tools.cpp:7117`) of the current root, calls `moveRoot` to reinsert the root at the
midpoint of the candidate branch, re-optimises all branch lengths for each candidate with up to
100 sweeps, and keeps the best. It is called during NNI search whenever the likelihood improves
(`tree/iqtree.cpp:3257-3259`) and from `ModelFactory::optimizeParameters` when `root_find` is on
(`model/modelfactory.cpp:1733-1734`). Separately, `--root-test` scans **every** branch and writes
a ranked list for the standard topology tests (`tree/phylotree.cpp:3130-3215`, driven from
`main/phyloanalysis.cpp:3957-3961`), and `--rootstrap` counts how often each branch is the root
across bootstrap replicates. Rooting by a non-reversible model is a first-class feature, not a
side effect.

**The unrooted non-reversible path.** If a non-reversible model is run on an unrooted tree
(only reachable via `kernel_nonrev`), pi is folded into the branch matrix as `pi_a P_ab(r_c t)`,
reproducing the reversible-style formula (`tree/phylokernelnonrev.h:1119-1130`). That is the only
place where the two non-reversible code paths differ mathematically.

### 2.4 Rate heterogeneity, because it is inside the same loop

Everywhere in every kernel the combination is exactly

```cpp
double len  = site_rate->getRate(c) * branch->getLength(c);
double prop = site_rate->getProp(c) * model->getMixtureWeight(m);
```

so a rate category scales the branch length and contributes with weight `rho_c w_m`.

**Discrete Gamma (`+G`), mean of category, the default.** `RateGamma::computeRatesMean`
(`model/rategamma.cpp:158-172`) implements Yang (1994) equations 9 and 10 with
`alpha = beta = gamma_shape` so the mean is 1:

```
q_k = chi2_{k/K}(2 alpha) / (2 alpha)                  cut points, the k/K quantile of Gamma(alpha, alpha)
F_k = I( alpha q_k , alpha + 1 )                       regularised incomplete gamma, shape alpha+1
r_k = K ( F_k - F_{k-1} ) = E[ X | q_{k-1} < X <= q_k ]
```

Weights are uniform, `getProp(c) = 1/K` (`model/rategamma.h:120`). The incomplete gamma is the
AS32 routine with tolerance 1e-8 (`model/rategamma.cpp:325-382`), the quantile is AS91. `+G`
with no number means `+G4`, since `num_rate_cats` defaults to 4 (`utils/tools.cpp:7265`). The
`--gamma-median` option instead uses the midpoint quantile `r_k = chi2_{(2k+1)/2K}(2 alpha)/(2 alpha)`
rescaled to mean 1 (`model/rategamma.cpp:109-127`).

**Invariable sites (`+I`).** `p_invar` is not an extra rate category. It is precomputed per
pattern (`tree/phylotreesse.cpp:547-612`) and added *after* the mixture sum:

```
iota_p =  p_inv * pi_a                    pattern constant at state a
          p_inv * sum_{a in A} pi_a       pattern constant at an ambiguity set A (B, Z, J)
          p_inv                           gap-only pattern
          0                               variable pattern
```

Under `+I+G` the Gamma rates are divided by `1 - p_inv` and the weights multiplied by it
(`model/rategammainvar.cpp:25-43`, `model/rategammainvar.h:69`), so `sum_c rho_c r_c = 1` holds
over all sites including the invariable ones. A change in either `p_inv` or `alpha` triggers a
full recomputation of the category rates (`model/rategammainvar.cpp:107-116`).

**FreeRate (`+R`).** Both rates and weights are free, and both constraints are enforced by
reparameterisation rather than by penalty (`model/ratefree.cpp:435-482`). With free variables
`u_i = w_i / w_K` and `v_i = r_i / r_K` for `i < K`,

```
w_i = u_i / (1 + sum_j u_j),   w_K = 1 / (1 + sum_j u_j)      so sum_i w_i = 1
S   = w_K + sum_{j<K} w_j v_j
r_i = v_i / S,                 r_K = 1 / S                    so sum_i w_i r_i = 1
```

Both identities hold exactly, by construction, at every optimiser step. The dimension is
`2K - 2` (`model/ratefree.cpp:286-296`).

**Fused rate models (`*G`, `*R`, `*H`).** Writing `*` instead of `+` ties the mixture class index
to the rate category index one-to-one (`model/modelfactory.cpp:818-864`). This is how LG4M and
LG4X work, and it changes `ncat_mix` from a product to a diagonal.

**Ascertainment bias (`+ASC`).** Lewis' correction subtracts
`N log(1 - sum_a Pr(constant at a))` with N the number of **sites**
(`tree/phylokernelnew.h:3262-3284`); Holder's `+ASC_MIS` uses a per-pattern denominator
(`tree/phylokernelnew.h:3238-3261`). The non-reversible kernel implements Lewis only, and
mixtures reject `+ASC` outright with `outError` (`model/modelmixture.cpp:4356`).

### 2.5 The general mixture form

The class index runs over the Cartesian product of mixture components and rate categories, with
`ncat_mix = ncat * nmixtures` normally and `ncat_mix = ncat` when fused
(`tree/phylokernelnew.h:2701-2703`, `2726`, `2749-2756`):

```
L_p = iota_p + sum_{m=1..M} sum_{c=1..K} w_m rho_c L_p^(m,c)(r_c t)          unfused
L_p = iota_p + sum_{c=1..K} w_c rho_c L_p^(c,c)                              fused
```

This is a **site-wise mixture**: no site is ever assigned a class. Posterior class
probabilities `P(c | site) = w_c L_c(site) / L(site)` are formed only inside the EM optimiser
(`model/modelmixture.cpp:4202-4206`), never in the likelihood itself. The contrast with
`ModelSet`, which assigns exactly one Q per site with no weights, is drawn in Section 10.

### 2.6 What the model layer actually optimises

Every model-layer objective is the negative full-tree log-likelihood at the current branch
lengths (`model/modelmarkov.cpp:1091-1117`):

```cpp
double ModelMarkov::targetFunk(double x[]) {
    bool changed = getVariables(x);
    if (changed) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
    }
    for (int i = 0; i < num_states; i++)
        if (state_freq[i] < 0 || (state_freq[i] > 0 && state_freq[i] < Params::getInstance().min_state_freq))
            return 1.0e+30;
    return -phylo_tree->computeLikelihood();
}
```

Three things to note. The objective has the side effect of writing the parameter vector into the
model. A full re-decomposition and a full invalidation of every partial likelihood happen on each
evaluation where anything changed, so one objective evaluation costs one complete tree traversal.
The `1.0e+30` return is a barrier, not an exception: infeasible regions are signalled by an
enormous function value, with `min_state_freq = MIN_FREQUENCY = 1e-4` by default
(`alignment/alignment.h:21`, `utils/tools.cpp:7263`).

---

## 3. The catalogue of amino-acid models

Every amino-acid model is constructed by `createModel` (`model/modelmixture.cpp:3122-3266`),
which dispatches on sequence type to `ModelProtein`, except for `UNREST` and the Lie-Markov
names, which are reachable through `ModelMarkov::getModelByName`
(`model/modelmarkov.cpp:1954-1963`). `ModelProtein::init` (`model/modelprotein.cpp:1106-1245`)
is a single `if / else` chain over the model-name string, and `GTR20` and `NONREV` are two of its
branches. There is no `ModelNonrev` class: `model/modelnonrev.{cpp,h}` are zero-byte
placeholders not in the build.

| family | reversible | free rate params | pi | estimation |
|---|---|---|---|---|
| named empirical (LG, WAG, ...) | yes | 0 | matrix's own, or `+F`, `+FO`, `+FQ` | nothing, or 19 frequencies |
| `GTR20` | yes | 189 | empirical by default | BFGS, Section 6 |
| `NQ.*` (6 matrices) | no | 0 | solved from the fixed Q | nothing |
| `NONREV` | no | 379 | solved from Q each evaluation | BFGS, Section 7 |
| `UNREST` | no | 379 on protein | solved from Q | BFGS, same as `NONREV` |
| C10 to C60, CF4 | yes | 0 (Poisson R) | K fixed profiles | weights, by EM or BFGS |
| `LG+C20+F` | yes | 0 (LG R fixed) | 20 fixed profiles plus 1 empirical | 20 free weights |
| EX2, EX3, EHO, UL2, UL3, EX_EHO | yes | 0 per class | each class's own | weights only |
| LG4M, LG4X | yes | 0 per class | each class's own | fused rate model, 1 or 6 params |
| GTRpmix (`--link-exchange`) | yes | 189 shared | K estimated or profile | EM then BFGS, Section 9 |
| PMSF (`-ft`, `-fs`) | yes | 0 | one vector per site | two-pass, Section 10 |
| MUTSEL | see Section 11 | see Section 11 | one vector per site | external, Section 11 |

### 3.1 The 28 built-in reversible matrices

`main/phylotesting.cpp:161-163`:

```cpp
const char* aa_model_names[] = {"LG", "WAG", "JTT", "Q.pfam", "Q.bird", "Q.mammal", "Q.insect",
    "Q.plant", "Q.yeast", "JTTDCMut", "DCMut", "VT", "PMB", "Blosum62", "Dayhoff",
    "mtREV", "mtART", "mtZOA", "mtMet", "mtVer", "mtInv", "mtMAM", "FLAVI",
    "HIVb", "HIVw", "FLU", "rtREV", "cpREV"};
```

The matrices themselves live in `builtin_prot_models`, a raw-string NEXUS block occupying
`model/modelprotein.cpp:31-1105`. A reversible entry is 190 numbers in **lower-triangular**
(PAML) order followed by a 20-entry frequency line. Since `rates[]` is an upper triangle,
`ModelProtein::readRates` (`model/modelprotein.cpp:1265+`) transposes on read with the index map
`id = col*(2*num_states - col - 1)/2 + (row - col - 1)`.

For these models `num_params = 0`, so `getNDim()` returns 0 and the matrix is not optimised at
all. The only free parameters a user can add are the frequencies (`+F` costs 19 degrees of
freedom in `getNDimFreq`, `+FO` costs 19 in `getNDim`) and the rate-heterogeneity parameters.

### 3.2 The six built-in non-reversible matrices

`main/phylotesting.cpp:165`:

```cpp
const char* aa_model_names_nonrev[] = {"NQ.bird","NQ.insect","NQ.mammal","NQ.pfam","NQ.plant","NQ.yeast"};
```

These are stored as a **full 20x20 Q with negative diagonal and rows summing to zero**, followed
by a 20-entry pi line (`model/modelprotein.cpp:915-1052`). `ModelMarkov::readRates` validates
every row sum against a tolerance of 1e-3 and throws otherwise
(`model/modelmarkov.cpp:1665-1698`). **Reversibility is auto-detected from the sign of the first
number in the file**: negative means a full Q, hence non-reversible
(`model/modelmarkov.cpp:1798-1810`). The published nQ matrices come from Dang et al. (2022),
whose method is exactly the procedure in Section 8.

### 3.3 Two fully parameterised models

`GTR20` is the general reversible 20-state model: 189 free exchangeabilities plus a frequency
vector. `ModelProtein::init` warns that it "will estimate 189 substitution rates that might be
overfitting", unless `--link-model` is in force. It seeds from `Params::gtr20_model`, default
`POISSON` (`utils/tools.cpp:7290`), overridable with `--init-exchange` or `--init-model`.

`NONREV` is the general non-reversible 20-state model: 379 free directed rates, pi derived. It
seeds from LG converted by `setReversible(false)`, overridable with `--init-model`, and carries
the same overfitting warning.

Neither name appears in any ModelFinder candidate array, so both are opt-in only, through
`-m` directly, `-mset`, or `-madd`.

### 3.4 `UNREST` on protein data

`ModelUnrest` (`model/modelunrest.{h,cpp}`, 55 + 138 lines) is written generically over
`num_states` and is reachable for any sequence type through `ModelMarkov::validModelName`. On
protein data `-m UNREST` gives 379 free rates, the same dimension as `NONREV`. It differs in
three details: the last rate is pinned to exactly 1.0 by an overridden `setRates`
(`model/modelunrest.cpp:106-110`), `setStateFrequency` is a deliberate no-op, and the seed is
all rates equal to 1.0 rather than LG. It is the shortest complete non-reversible model class in
the tree and the right template for a new one.

### 3.5 Where the profile mixtures live

`builtin_mixmodels_definition` (`model/modelmixture.cpp:23-3115`) is a second embedded NEXUS
block holding 159 `model` declarations and 164 `frequency` declarations, including the CAT
profile families `C10pi1..10`, `C20pi1..20`, `C30pi1..30`, `C40pi1..40`, `C50pi1..50`,
`C60pi1..60`, and `Fclass1..4`. Structure and estimation are in Section 9.

---

## 4. The shared numerical machinery

This section is where the answers depart most sharply from what "maximum likelihood" suggests.

### 4.1 Empirical frequencies are not simple counts

`Alignment::computeStateFreq` does not divide observed counts by their total. It runs an eight
round EM-style fixed point that distributes ambiguous characters in proportion to the current
estimate (`alignment/alignment.cpp`, `convertCountToFreq`):

```
pi^(0)_j = 1/20
repeat 8 times:
    for each observed character code s with count c_s:
        a_j    = 1 if state j is compatible with s, else 0     (all ones for X or gap)
        num_j  = pi_j a_j
        add  c_s * num_j / sum_k num_k  to the accumulator
    pi <- accumulator, normalised
```

The compatibility vectors come from `Alignment::getAppearance`: B resolves to {N, D}, Z to
{Q, E}, J to {I, L}, and X or a gap to all 20 states. The iteration count is a hard-coded
`NUM_TIME = 8`, not a convergence test.

A post-processing step, `Alignment::convfreq`, raises any frequency below
`Params::min_state_freq` (default 1e-4) to that floor and absorbs the deficit into the largest
frequency, unless `--keep-ident`-style `keep_zero_freq` is set. So an amino acid absent from the
alignment gets 1e-4, not 0, which matters because `ZERO_FREQ = 1e-10`
(`utils/eigendecomposition.h:24`) would otherwise delete the state from the matrix entirely.

### 4.2 Eigendecomposition, reversible

The default technique is `MET_EIGEN3LIB_DECOMPOSITION` (`utils/tools.cpp:7533`). The reversible
path builds Q, symmetrises it by conjugation with `diag(sqrt(pi))`, and calls Eigen's
`SelfAdjointEigenSolver` (Householder tridiagonalisation plus implicit-shift QL), then maps back
(`model/modelmarkov.cpp:1450-1601`):

```cpp
Q = pi_sqrt * Q * pi_sqrt_inv;                          // S = Pi^(1/2) Q Pi^(-1/2)
SelfAdjointEigenSolver<MatrixXd> eigensolver(Q);
evec     = pi_sqrt_inv * eigensolver.eigenvectors();    // U      = Pi^(-1/2) V
inv_evec = eigensolver.eigenvectors().transpose() * pi_sqrt;   // U^(-1) = V^T Pi^(1/2)
```

Three tolerance tests can abort this path and fall back to the older hand-written
`eigensystem_sym` (Numerical Recipes `tred2` plus `tqli`):

| test | threshold | action |
|---|---|---|
| asymmetry of S, `max abs(S - S^T)` | `>= 0.01` | diagnostic dump, continue |
| asymmetry of S | `> 0.1` | fall back to `decomposeRateMatrixRev` |
| `eigensolver.info()` | not `Success` | fall back |
| largest eigenvalue | `> 1e-4` | fall back (eigenvalues must be non-positive) |

A NaN eigenvalue is silently replaced by 0. `F81`-style models with `num_params == -1` use a
closed-form decomposition instead, with `mu = total_num_subst / (1 - sum_i pi_i^2)` and all
non-zero eigenvalues equal to `-mu`.

States with `pi_i <= ZERO_FREQ = 1e-10` are **removed** from the matrix before decomposition and
reinserted afterwards as identity rows and columns with eigenvalue 0
(`model/modelmarkov.cpp:1563-1597`, `utils/eigendecomposition.cpp:521-555`). Any procedure that
supplies a pi with near-zero entries must expect those states to drop out silently.

### 4.3 Eigendecomposition, non-reversible

A non-reversible Q is real and generally non-symmetric, so its spectrum is complex and it may be
defective. IQ-TREE calls Eigen's `EigenSolver` (real Schur by Francis double-shift QR, then
eigenvector extraction) and stores complex eigenvalues and eigenvectors in `ceval`, `cevec`,
`cinv_evec` (`model/modelmarkov.cpp:1292-1387`). Those three arrays are
`std::complex<double>` aliases over the same allocation as the real ones, so they must never be
freed or reallocated independently.

The inverse eigenvector matrix is obtained by `FullPivLU<MatrixXcd>::inverse()`, guarded by
`isInvertible()`. If the eigenvector matrix is singular, `nondiagonalizable` is set and a warning
is printed. There is no explicit numeric threshold in IQ-TREE here: Eigen's default full-pivot
rank threshold is used.

`eigensystem_nonrev` (`utils/eigendecomposition.cpp:319-470`) is the alternative, reached only
under `--eigen`. It is the classical EISPACK chain `elmhes` (Gaussian Hessenberg reduction with
partial pivoting), `eltran` (accumulate the transformations), `hqr2` (double-shift implicit QR,
iteration budget `30n`, exceptional shift every ten iterations), and `luinverse` (Crout LU with
implicit scaling). It returns complex eigenvalues as adjacent real and imaginary parts with the
eigenvector matrix in real quasi-Schur form. **This path is broken for non-reversible models in
release builds**: see Section 15.

### 4.4 Matrix exponentiation

**Reversible.** `P(t) = U exp(Lambda * t / total_num_subst) U^(-1)`, computed as a dense triple
product by Eigen (`model/modelmarkov.cpp:507-568`). Since the decomposition already scaled Q by
`total_num_subst`, the two factors cancel and the effective rate is always 1. On x86 the default
path applies **no clamping and no row-sum renormalisation**; only the ARM NEON branch clamps
negative entries to 0, fills the lower triangle by detailed balance
`P_ji = (pi_i/pi_j) P_ij`, and closes each row by setting the diagonal to `1 - sum of the row`.

**Non-reversible.** `P(t) = exp(Q t)` with `t` used **raw**, because `Q` was already scaled by
`delta = total_num_subst / mu` inside `decomposeRateMatrixNonrev`
(`model/modelmarkov.cpp:463-505`). The complex spectrum is exponentiated in full complex
arithmetic and the imaginary part discarded:

```cpp
ceval_exp = (eval*time).exp().matrix();
MatrixXcd res = cevectors * ceval_exp.asDiagonal() * cinv_evectors;
map_trans = res.real();
```

A row-sum sanity test follows. If any row sum leaves `[0.9999, 1.0001]`, the call **recomputes**
that P(t) by scaling and squaring, which is Eigen's `MatrixBase::exp()` implementing Higham's
algorithm with degree-13 diagonal Pade approximants and adaptive scaling. That fallback is also
taken unconditionally when `nondiagonalizable` is set, or when `--scaling-squaring` is passed.
The scaling-and-squaring branch only warns about negative entries, it does not fix them.

**The `total_num_subst` asymmetry.** For reversible models the factor is applied in the
decomposition and divided out at exponentiation time, so it cancels. For non-reversible models it
is applied in the decomposition and *not* divided out. With the default value of 1.0 the two
agree exactly. They would differ for a non-reversible mixture component with a per-class rate,
which is not a configuration IQ-TREE currently reaches.

### 4.5 The optimiser: BFGS with numerical gradients

This is the single most consequential implementation fact in the document.

**The routine is `Optimization::minimizeMultiDimen` wrapping `dfpmin`**
(`utils/optimization.cpp:750-906`). Despite the Numerical Recipes name, `dfpmin` implements
**BFGS**, not DFP, on the inverse Hessian. With `s = x_{k+1} - x_k`, `y = g_{k+1} - g_k`,
`H` the approximate inverse Hessian, and `u = s/(y^T s) - (H y)/(y^T H y)`:

```
H_{k+1} = H_k + (s s^T)/(y^T s) - (H_k y y^T H_k)/(y^T H_k y) + (y^T H_k y) u u^T
```

The first two terms are the DFP update; the third is the BFGS correction. The update is skipped
entirely when `(y^T s)^2 <= EPS * ||y||^2 * ||s||^2` with `EPS = 3.0e-8`. The new direction is
`-H g` with no positive-definiteness check and no reset to steepest descent on failure.

**Constants** (`utils/optimization.cpp:781-785`, `639-640`):

| constant | value | role |
|---|---|---|
| `ITMAX` | 200 | maximum BFGS iterations, exhaustion is silent |
| `EPS` | 3.0e-8 | curvature safeguard |
| `TOLX` (dfpmin) | 1.2e-7 | relative parameter-change convergence test |
| `TOLX` (lnsrch) | 1.0e-7 | minimum relative step length |
| `STPMX` | 100.0 | maximum step length, `100 * max(norm of x0, n)` |
| `ALF` | 1.0e-4 | Armijo sufficient-decrease constant |
| `ERROR_X` | 1.0e-4 | finite-difference step fraction |
| `MAX_ITER` | 3 | maximum random restarts |

**Convergence tests**, in order, each iteration:

```
test 1 (parameters):  max_i |s_i| / max(|x_i|, 1)  <  1.2e-7
test 2 (gradient)  :  max_i |g_i| max(|x_i|, 1) / max(|f|, 1)  <  gtol
```

For amino-acid models `gtol = max(gradient_epsilon, TOL_RATE)` with `TOL_RATE = 1e-4`.

**Gradients are numerical, one-sided, forward differences.**
`Optimization::derivativeFunk` (`utils/optimization.cpp:916-939`):

```cpp
double fx = targetFunk(x);
for (dim = 1; dim <= ndim; dim++) {
    temp = x[dim];
    h[dim] = ERROR_X * fabs(temp);
    if (h[dim] == 0.0) h[dim] = ERROR_X;
    x[dim] = temp + h[dim];
    h[dim] = x[dim] - temp;           // exactly representable step
    dfx[dim] = targetFunk(x);
    x[dim] = temp;
}
for (dim = 1; dim <= ndim; dim++) dfx[dim] = (dfx[dim] - fx) / h[dim];
```

So `d f / d x_i ≈ [ f(x + h_i e_i) - f(x) ] / h_i` with `h_i = 1e-4 |x_i|`, or `1e-4` when
`x_i` is exactly zero. No substitution model anywhere in IQ-TREE overrides `derivativeFunk`; the
only override in the entire tree is `PhyloTreeMixlen::derivativeFunk` for heterotachy branch
lengths (`tree/phylotreemixlen.cpp:476`). **Every amino-acid rate and frequency gradient is a
numerical forward difference.**

Three consequences follow, and they should be stated in any methods section that reports a
`NONREV` or `GTR20` fit.

*Cost.* One gradient costs `ndim + 1` full-tree likelihood evaluations. For `NONREV` that is 380
traversals per gradient; for `GTR20+FO`, 209. The worst case with 200 iterations and 4 BFGS runs
is of order 150,000 likelihood evaluations.

*Accuracy.* A forward difference has truncation error of order `h |f''| = 1e-4 |f''|`. The
gradient convergence tolerance is also 1e-4. The stopping test therefore sits at the noise floor
of the gradient it is testing. The step size is also far from the `sqrt(machine epsilon)` that
would minimise the truncation-plus-roundoff tradeoff for an exactly computed function, though it
is defensible given the summation noise of a log-likelihood.

*Bounds.* The probe `x_i + h_i` is not clamped, so a parameter sitting at the upper bound of 100
is probed at 100.01, outside the feasible box, and the resulting one-sided derivative is taken
from the infeasible side.

**Box constraints are handled by clamping the trial point, not by projection of the direction.**
The only constraint mechanism is `fixBound` (`utils/optimization.cpp:149-156`), called once per
line-search trial, after forming `x = x_old + lambda p` and before evaluating the objective. The
Armijo test then compares the objective at the *clamped* point against a predicted decrease
computed from the *unclamped* direction, so the two are inconsistent whenever clamping is active.
When the current point is on a bound and the direction points outward, the clamped point equals
the old point, the Armijo test can never be satisfied, backtracking runs until the step falls
below `alamin`, and the parameter stays pinned at the boundary. This failure mode is the reason
`restartParameters` exists.

**Bounds for amino-acid rates** are `[MIN_RATE, MAX_RATE] = [1e-4, 100]`
(`model/modelmarkov.h:30-32`, `setBounds` at `model/modelmarkov.cpp:1139-1164`), and for free
frequencies `[min_state_freq, 1.0] = [1e-4, 1.0]`.

**Restarts.** `restartParameters` (`utils/optimization.cpp:726-748`) triggers when any
coordinate flagged in `bound_check[]` lies within 1e-4 of a bound, and then re-randomises **every**
coordinate uniformly on the lower third of its box, `x_i ~ U[l_i, l_i + (u_i - l_i)/3]`, up to
three times, keeping the best of the runs. This uses `random_double()`, so restarts are not
reproducible unless the seed is fixed.

`ModelMarkov::setBounds` sets `bound_check[i] = false` for every coordinate, so **plain `GTR20`
and `NONREV` never restart**. The exception matters: `ModelMixture::setBounds` forces
`bound_check[i] = true` and `upper_bound[i] = 100` for all protein exchangeabilities when
`optimize_linked_gtr` is on (`model/modelmixture.cpp:4648-4652`, `4663-4673`), so the GTRpmix
procedure of Section 9 does restart, stochastically, and is the one amino-acid procedure whose
result depends on the random seed through the optimiser.

**L-BFGS-B is present, vendored, and not used for matrices.** `Optimization::L_BFGS_B`
(`utils/optimization.cpp:1118-1169`) wraps a translation of the Byrd, Lu, Nocedal, and Zhu
algorithm with `m = 10` stored pairs, `factr = 1e7`, `nbd[i] = 2` for every variable, and
`maxit` defaulting to 5. The call in `ModelMarkov::optimizeParameters` is **commented out**, with
the note `2019-09-05: REMOVED due to numerical issue (NAN) with L-BFGS-B`
(`model/modelmarkov.cpp:1198-1224`). It remains live for `RateFree` and for mixture-length
branch optimisation. A genuinely bound-constrained optimiser therefore exists in the tree and is
not used for amino-acid matrices.

**No standard errors.** `dfpmin` allocates its inverse Hessian locally and frees it on every
return path, with no write-back. The converged `H`, which approximates the asymptotic covariance
of the rate parameters, is discarded. IQ-TREE reports no standard errors for amino-acid
exchangeabilities or directed rates, and copying `hessin` out before `FREEALL` is the hook that
would provide them.

### 4.6 The alternating optimisation loop

`ModelFactory::optimizeParameters` (`model/modelfactory.cpp:1570-1700`) is the outer loop for a
single alignment:

```
cur_lh <- computeLikelihood()
for i = 2 .. num_param_iterations:                        # default 100
    optimise branch lengths, min(i,3) sweeps, tolerance logl_epsilon
    optimise model parameters      -> ModelMarkov::optimizeParameters(gradient_epsilon)
    optimise rate parameters       -> site_rate->optimizeParameters(gradient_epsilon)
    if the improvement is <= logl_epsilon: final branch-length pass, break
rescale site rates so the mean rate is 1, absorbing the factor into branch lengths
```

Model and rate parameters are optimised **sequentially**, not jointly, unless
`--opt-model-rate-joint` is passed, in which case `ModelFactory::optimizeAllParameters` runs one
BFGS over the concatenated vector (`model/modelfactory.cpp:1331-1377`). Defaults:
`num_param_iterations = 100`, `modelEps = 0.01`, `loglh_epsilon = 0.001`, and
`modelfinder_eps = 0.1` for candidate evaluation inside ModelFinder.

During tree search the whole loop is re-entered every time NNI improves the likelihood by more
than `modelEps` (`tree/iqtree.cpp:3244-3254`), with the tolerance loosened tenfold, and for a
rooted tree the root position is re-optimised immediately afterwards.

### 4.7 Branch lengths

`optimizeOneBranch` (`tree/phylotree.cpp:2635-2694`) uses a **safeguarded Newton root-finder on
the first derivative** by default (`minimizeNewton`, `utils/optimization.cpp:422-506`), falling
back to bisection whenever the second derivative is non-negative or the Newton step would leave
the bracket. Tolerance is `params->min_branch_length`, the step cap is `maxNRStep = 100`, and the
upper bound is `max_branch_length = 10.0`.

The derivatives are **analytic**, unlike the model-parameter gradients. For reversible models the
eigen-space representation gives all three quantities from one dot product
(`tree/phylokernelnew.h:2326-2364`, `2513-2517`):

```
val0_i = w_c exp(lambda_i r_c t),    val1_i = lambda_i r_c val0_i,    val2_i = (lambda_i r_c)^2 val0_i
L_p = <val0, theta>,   dL_p/dt = <val1, theta>,   d2L_p/dt2 = <val2, theta>
d log L / dt   = sum_p n_p L'_p / L_p
d2 log L / dt2 = sum_p n_p [ L''_p / L_p - (L'_p / L_p)^2 ]
```

where `theta` is the branch-independent product of the two conditional vectors, cached so that
re-evaluating at a new `t` costs one dot product. For non-reversible models the same derivatives
come from explicit matrices, `dP/dt = Q P(t)` and `d2P/dt2 = Q^2 P(t)`
(`model/modelmarkov.cpp:735-768`), and the `theta` cache is disabled
(`computeLikelihoodFromBufferPointer = nullptr`), so every branch-length evaluation recomputes
the branch likelihood in full. That is a large part of the roughly twofold cost of a
non-reversible analysis.

A branch incident on the root of a rooted tree is skipped entirely.

### 4.8 The initial tree

Model selection and model estimation both need a tree before they have a model. The default
start tree is a randomised stepwise-addition parsimony tree (`STT_PLL_PARSIMONY`,
`utils/tools.cpp:7461`), with `-t BIONJ` and `-t <file>` as alternatives
(`tree/iqtree.cpp:573-760`).

---

## 5. Procedure A: a fixed empirical matrix

Invocation: `-m LG`, `-m LG+F`, `-m WAG+FO+G4`, `-m NQ.pfam`.

**What is estimated.** Nothing about the matrix. `ModelProtein::init` reads the NEXUS block,
calls `rescaleRates` so the largest entry is 10.0, and sets `num_params = 0`. `getNDim()`
returns 0, so `ModelMarkov::optimizeParameters` returns immediately. The only free parameters are
the frequencies, if `+F` or `+FO` was requested, and the rate-heterogeneity parameters.

**Reversible matrices.** Q is rebuilt as `Q_ij = R_ij pi_j`, row sums zeroed, rescaled to unit
expected rate. With `+FO` the 19 free frequencies are optimised by the BFGS of Section 4.5, with
each evaluation requiring a fresh eigendecomposition because `U^(-1) = V^T Pi^(1/2)` depends on
pi.

**Non-reversible matrices (`NQ.*`).** The file's Q is used verbatim. `readParametersString`
calls `decomposeRateMatrix` once as a consistency check, which re-solves pi from Q and warns on
any disagreement with the file's pi beyond 1e-3. The final `freq_type` is `FREQ_USER_DEFINED`
carrying the **solved** vector, not the file's. The tree is rooted by `setReversible(false)`, the
root position is then subject to the search in Section 2.3, and the non-reversible kernel runs.
Adding `+F` changes nothing, for the reason in Section 1.5.

**Cost.** This is the cheapest amino-acid analysis, and the non-reversible variants cost roughly
twice a reversible one from the complex arithmetic, the loss of the `theta` cache, and the
doubled kernel memory.

---

## 6. Procedure B: `GTR20`, the general reversible matrix

Invocation: `-m GTR20`, `-m GTR20+FO+G4`.

**Objective.** The full-tree log-likelihood of one alignment, `-targetFunk`.

**Free parameters.** 189 exchangeabilities, the 190th pinned, plus 19 frequencies if `+FO`.
Default `freq_type` is `FREQ_EMPIRICAL`, so by default the frequencies are counted once and held
fixed and `getNDim() = 189`.

**Feasible set.** The box `[1e-4, 100]^189`, which is a relaxation of the true constraint
`R_ij >= 0`. The lower bound of 1e-4 means an exchangeability can never actually reach zero.

**Seed.** `Params::gtr20_model`, default `POISSON` (all exchangeabilities equal), overridable
with `--init-exchange NAME` or `--init-model NAME|FILE`. A non-reversible seed is rejected with
`outError("Cannot initialize from non-reversible model")`.

**Algorithm.** One call to `minimizeMultiDimen` per outer round, inside the alternating loop of
Section 4.6. BFGS, forward-difference gradients at 190 likelihood evaluations each, box
clamping, and no restarts because `bound_check` is false throughout. After the optimiser returns,
`scaleStateFreq(true)` renormalises pi if it was free, then `decomposeRateMatrix` and
`clearAllPartialLH` and one final likelihood evaluation.

**Normalisation.** `sum_i pi_i (-Q_ii) = 1` after every decomposition. The pinned 190th
exchangeability, at its seed value, is the only thing fixing which point of the scale ray the
free parameters live on.

**Practical note.** 189 parameters estimated from one alignment is the overfitting IQ-TREE warns
about at construction. The published `Q.*` matrices were not estimated this way; they came from
Procedure D.

---

## 7. Procedure C: `NONREV`, the general non-reversible matrix

Invocation: `-m NONREV`, `-m NONREV+G4`.

**Objective.** The full-tree log-likelihood of one alignment, on a **rooted** tree.

**Free parameters.** 379 of the 380 directed rates. pi contributes **no** free dimensions, since
`ModelMarkov::getNDim` returns `num_params` early for non-reversible models. Default `freq_type`
is `FREQ_ESTIMATE`, which for a non-reversible model is a misnomer: nothing is estimated, and its
only effect is to leave the solve-from-Q guard true.

**Feasible set.** The box `[1e-4, 100]^379` on the off-diagonals. Zero row sums hold by
construction, and the stationarity of pi holds because pi is derived rather than constrained.
There is **no equality-constraint solver anywhere in IQ-TREE**; box bounds are the only
constraint mechanism.

**Seed.** LG, converted by `setReversible(false)` through `Q_ij = R_ij pi_j` with LG's own pi,
then repacked to 380 entries. `--init-model NAME|FILE` overrides, and a reversible seed is
converted automatically.

**Per-evaluation work.** Each of the roughly 380 likelihood evaluations per gradient does, in
order: unpack 379 free rates plus the pinned one into a 20x20 matrix; derive the diagonal; reset
pi to uniform; solve `pi^T Q = 0` by column-pivoted QR; rescale Q so `sum_i pi_i (-Q_ii) = 1`;
compute a complex eigendecomposition and invert the complex eigenvector matrix; invalidate every
partial likelihood; traverse the tree.

**Root.** `setReversible(false)` roots the tree by midpoint (or at the `-o` outgroup) as a side
effect of model construction. The root position is then optimised by local search during tree
search, and `--root-test` will scan all branches.

**Degrees of freedom.** 380 off-diagonal entries, minus one for scale, gives 379. That is the
free-parameter count IQ-TREE reports. For comparison, constraining Q to have a *supplied*
stationary distribution would remove a further `n - 1 = 19` degrees of freedom, because
`pi^T Q = 0` is a set of linear equations in the entries of Q whose rank is `n - 1` rather than
`n` (the n equations sum to zero identically given zero row sums). A pi-constrained nQ would
therefore have 360 free parameters, nested inside the 379 of `NONREV`, which is the arithmetic
this fork's project rests on. No such constrained model exists in IQ-TREE 3 for 20 states;
`ModelLieMarkov::setBasis` implements exactly that idea for 4-state DNA.

---

## 8. Procedure D: QMaker and nQMaker, one matrix across many alignments

Invocation: `iqtree3 -s <concat> -p <partitions.nex> --model-joint GTR20` for QMaker and
`--model-joint NONREV` for nQMaker. This is how the published `Q.*` and `NQ.*` matrices were
estimated, and it is the procedure to imitate for any new amino-acid matrix.

**The objective.** With partitions `k = 1..K`, data `D_k`, a shared topology `T`, per-partition
rate multipliers `r_k`, branch lengths `b`, and per-partition rate-heterogeneity parameters
`phi_k`, the shared parameter vector `theta` maximises

```
L(theta) = sum over k in S(M) of  log Pr( D_k | Q(theta), T, r_k b, phi_k )
```

where `S(M)` is the set of partitions carrying the linked model name M.
`PartitionModel::targetFunk` (`model/partitionmodel.cpp:299-338`) computes it as an
**unweighted** sum of per-partition log-likelihoods. There is no explicit per-partition weight,
so a partition contributes in proportion to its size, which is the intended QMaker objective.
Partitions with a different model name contribute zero during the theta update and are added back
by a full supertree recomputation afterwards.

The sum is accumulated into a per-partition array and summed **serially in a fixed order** after
an OpenMP parallel loop, rather than by a `reduction`. The source comment says why: this makes
the objective bit-reproducible independent of thread count, which matters because the BFGS
gradient is a finite difference and would otherwise be perturbed by scheduling noise. This is a
detail worth recording in a methods section.

**One object or many.** Every partition has its own `ModelFactory` and its own `ModelSubst`.
`linked_models` maps a model *name* to the **first** partition's model, which acts as the
representative. `setVariables` packs from the representative only; `getVariables` broadcasts to
every partition whose model name matches (`model/partitionmodel.cpp:733-744`). Matching is by
name string, not by pointer. Any new parameterisation must therefore be a pure function of the
variable vector plus per-model constants, with no hidden per-partition state.

**The `fixed_parameters` gate.** Every per-partition model is constructed with
`fixParameters(true)`, which forces `getNDim()` to 0 and makes `optimizeParameters` a no-op, so
the per-partition optimisers never touch the matrix. `optimizeLinkedModels`
(`model/partitionmodel.cpp:842-868`) unfixes all partitions of one linked model, runs the BFGS,
re-fixes them, checkpoints, and moves to the next linked model. Distinct linked models are
optimised sequentially, in `unordered_map` order.

**The full alternation.** Two different outer loops exist and the flag chooses between them.

| flag | edge model | factory | outer rounds | stopping test |
|---|---|---|---|---|
| `-q` / `--edge equal` | all partitions share `b` exactly | `PartitionModelPlen` | 100 | improvement `< logl_epsilon` |
| `-p` / `--edge scale` | partition k uses `r_k b` | `PartitionModelPlen` | 100 | improvement `< logl_epsilon` |
| `-Q` / `-M` / `--edge unlink` | each partition has its own `b^(k)` | `PartitionModel` | 10 | improvement `< 10 logl_epsilon` |
| `-S` | separate topologies too | `PartitionModel` | 10 | as above |

The published invocation uses `-p`, so the default QMaker scheme is edge-linked-proportional with
`PartitionModelPlen` and up to 100 rounds. One round is
(`model/partitionmodelplen.cpp:75-194`):

```
for i = 1 .. 99:
  1. in parallel over partitions: optimise phi_k only (no branch lengths)
     rescale each partition's site rates to mean 1, folding the factor into b and r_k
  2. if --link-alpha: optimise one shared Gamma shape by 1-D Brent
  3. for each linked model M, sequentially:
       theta_M <- argmax of the summed log-likelihood, by BFGS
       broadcast theta_M to every partition in S(M), re-decompose, checkpoint to disk
  4. optimise the partition rate multipliers r_k, then renormalise so
     sum_k r_k n_k / sum_k n_k = 1, folding the residual into b
  5. optimise b on the supertree, min(5, i+1) sweeps
  until |change in log-likelihood| < logl_epsilon
```

The number of rounds for the `-Q` path is controlled by `--loop-model`, default 10
(`utils/tools.cpp:7239`).

**Parameter counting.** The shared matrix is counted **once**, not K times
(`model/partitionmodel.cpp:242-256`), so a `--model-joint NONREV` analysis of K partitions adds
379 to the degrees of freedom, not `379 K`.

**Seeding.** `NONREV` from LG, `GTR20` from POISSON, or `--init-model NAME|FILE`.
`--init-model DIVMAT`, the empirical divergence matrix, is **implemented but disabled**: the code
path opens with `ASSERT(0 && "init_by_div_mat not working")`, and the intended
`Q = log(diag(pi) F)` step is commented out. `ASSERT` compiles to a no-op only under `NDEBUG`,
which the gcc and clang Release configuration does not define (see Section 15, item 9), so in
those builds the path stops at the assertion (not run). In a build with
`NDEBUG`, the remaining code would run, handing `setFullRateMatrix` the raw row-normalised
divergence matrix rather than its logarithm and applying it only to the representative model
rather than broadcasting it. Do not use `--init-model DIVMAT`.

For the record, the divergence matrix itself (`alignment/alignment.cpp:5617-5674`) counts, over
all sites and all unordered taxon pairs, how often each pair of states co-occurs:

```
D_ij = sum_p f_p n_p(i) n_p(j)          i != j
D_ii = sum_p f_p n_p(i) (n_p(i) - 1) / 2
```

then symmetrises and row-normalises. There is no tree and no correction for shared ancestry. It
also has a genuine inconsistency: the state counts that form pi omit the `* it->frequency`
weighting that every pair count has, so pi would be computed over distinct patterns while D is
computed over sites.

**Interleaving with tree search.** There is no separate two-pass driver. The shared-matrix step
sits inside the ordinary IQ-TREE search: every time NNI improves the likelihood by more than
`modelEps`, the entire alternation including `optimizeLinkedModels` is re-run on the new topology
at a tenfold looser tolerance. The training/testing split and the "fix the trees, estimate Q,
iterate" structure of the published QMaker protocol is a shell-level workflow, typically a
fixed-tree run (`-te` or `-n 0`) followed by feeding the estimated matrix back with
`--init-model`.

**Citations printed.** The discriminator is purely whether `model_joint` contains the substring
`NONREV` (`main/phyloanalysis.cpp:164-179`). Anything else, including `GTR20+G`, is attributed to
Minh, Dang, Vinh and Lanfear (2021), *Systematic Biology* 70:1046-1060, doi
10.1093/sysbio/syab010. Anything containing `NONREV` is attributed to Dang, Minh, McShea, Masel,
James, Vinh and Lanfear (2022), *Systematic Biology* 71:1110-1123, doi 10.1093/sysbio/syac007.
`--link-model` without `--model-joint` triggers neither.

---

## 9. Procedure E: profile mixtures, full-matrix mixtures, and GTRpmix

### 9.1 What a profile mixture is

A profile mixture pairs **one** exchangeability matrix with **many** frequency vectors. Because
detailed balance holds for any positive pi paired with a symmetric R, each class
`Q_k = R diag(pi_k)` is a valid reversible generator, and the site likelihood is the weighted sum
of Section 2.5. This construction has no naive non-reversible analogue, which is one motivation
for constrained nQ work.

The CAT families C10 to C60 are defined in `builtin_mixmodels_definition`
(`model/modelmixture.cpp:735-977`) as, for example:

```
model C20 = POISSON+G+FMIX{C20pi1:1:0.0559910600, ... , C20pi20:1:w20};
```

Reading the syntax `name:rate:weight`: the exchangeability matrix is **POISSON** (all `R_ij = 1`,
the F81 case, handled by the closed-form eigendecomposition), there are K frequency vectors, each
class has relative rate 1, and each has a fixed weight from the original CAT analysis. Because
every component carries an explicit weight, `fix_prop` is set and the C-series weights are
**fixed**. The `C*Opt` variants omit the weights and leave them free. Writing `LG+C20` replaces
POISSON with LG.

`LG+C20+F` is a **21-class** mixture, not 20: a trailing `+F` on a `+FMIX` model prepends an
extra empirical-frequency class and sets `optimize_mixmodel_weight`, which frees all the weights
(`model/modelfactory.cpp:518-524`, `model/modelmixture.cpp:3585-3587`). `+FO` does the same with
an estimated rather than counted profile.

### 9.2 Full-matrix mixtures and fusion

`EX2`, `EX3`, `EHO`, `UL2`, `UL3`, `EX_EHO` are mixtures of complete reversible models, each
class carrying its own 190 exchangeabilities and its own pi, with a fixed relative rate written
into the definition and free weights.

`LG4M` and `LG4X` are four complete models **fused** with a rate model:

```
model LG4M = MIX{LG4M1,LG4M2,LG4M3,LG4M4}*G4;
model LG4X = MIX{LG4X1,LG4X2,LG4X3,LG4X4}*R4;
```

The `*` operator ties class index to rate-category index one-to-one. With fusion,
`getMixtureWeight(m)` returns 1.0 for every class and each class's `total_num_subst` is forced to
1.0, so the class probability comes entirely from `site_rate->getProp(c)` and the rate variation
entirely from the rate model. LG4M therefore has **one** free parameter, the Gamma shape; LG4X has
six, the free-rate weights and rates.

### 9.3 Normalisation of a mixture

Each class Q is separately normalised so `sum_i pi^k_i (-Q^k_ii) = mu_k`, where `mu_k` is that
class's `total_num_subst`. The mixture-level invariant is imposed at construction and re-imposed
after every optimisation round (`model/modelmixture.cpp:3572-3580`, `4438-4451`):

```
sum_k w_k mu_k = 1
```

so a branch length is expected substitutions per site averaged over the mixture. Fused models are
explicitly excluded from the re-imposition, keeping `mu_k = 1` and letting the rate model carry
the scale. For the C-series every `mu_k` is 1 before and after, so the relative-rate degree of
freedom is unused and only the profiles differ.

### 9.4 How mixture parameters are estimated

Weights are parameterised on the simplex by ratio to the last class
(`model/modelmixture.cpp:4579-4628`), not by log-ratio: free variables `x_i = w_i / w_{K-1}` map
back as `w_i = x_i / (1 + sum_j x_j)` and `w_{K-1} = 1 / (1 + sum_j x_j)`. Bounds are
`[0.001, 1000]`.

Three optimisation regimes exist (`model/modelmixture.cpp:4348-4453`):

**BFGS**, the default for non-fused Q-mixtures. One joint `minimizeMultiDimen` over the whole
flattened vector, all class parameters plus `K - 1` weight ratios, with
`ModelMixture::targetFunk` re-decomposing every class that has free parameters and then
recomputing the tree likelihood.

**EM**, used when the model is fused, when `-optalg_qmix EM` is passed, or when
`--link-exchange` is on. `optimizeWithEM` (`model/modelmixture.cpp:4140-4314`) does an E step
giving `P(class c | pattern p)`, an M step updating the weights, and then optimises each class
**one at a time** against a private tree whose pattern frequencies have been replaced by that
class's posteriors. This is the standard EM decoupling: maximise
`sum_p P(c|p) log L_c(p)` per class.

**Weights-only EM**, `optimizeWeights` (`model/modelmixture.cpp:4066-4138`), the closed-form
update `w_c <- (1/N) sum_p n_p P(c | p)` from Wang, Li, Susko and Roger (2008), with a floor of
1e-10 and a convergence test of 1e-4 on every weight.

Mixtures reject `+ASC` with an error.

### 9.5 GTRpmix: estimating the exchangeability matrix shared by a profile mixture

Invocation: `--link-exchange` (or `--link-exchange-rates`) with a model string containing
`GTR` or `GTR20`, for example `-m GTR20+C60+G`. The option also forces
`optimize_alg_qmix = "EM"` and `reset_method = "const"` (`utils/tools.cpp:1255-1259`). A second,
implicit trigger rewrites `R+F<k>` for integer k into `MIX{R+FO, ... k times}` and sets the same
flag (`utils/tools.cpp:1014-1079`), so `GTR20+F60+G` means sixty free profiles sharing one
estimated R.

**The objective is the ordinary mixture log-likelihood.** What the option changes is parameter
sharing: one 190-entry exchangeability vector is shared by all K classes while each class keeps
its own pi_k and w_k.

**The procedure, once per optimisation round:**

1. Freeze R by setting `num_params = 0` on every class, and renormalise R so its last entry is 1.
2. Run EM: posteriors, then weights, then each class's pi_k against posterior-weighted pattern
   frequencies.
3. Run `optimizeLinkedSubst` (`model/modelmixture.cpp:4464-4537`): a single BFGS over **189**
   exchangeabilities, packed from class 0 with its `freq_type` temporarily forced to
   `FREQ_USER_DEFINED` so frequencies are excluded, and unpacked into **every** class. The
   objective is the full mixture likelihood.

Step 3 is the one amino-acid procedure where `bound_check` is true, upper bounds are forced to
100, and `restartParameters` therefore fires, re-randomising all 189 coordinates up to three
times when any lands within 1e-4 of a bound. Its result depends on the random seed.

The seed for R is `Params::gtr20_model`, default POISSON, or the reference matrices `EAL` and
`ELM` supplied in the NEXUS block for exactly this purpose
(`model/modelmixture.cpp:984`, `1008`), or any matrix given to `--init-exchange`.

**Output.** `<prefix>.GTRPMIX.nex` holds **only** the estimated exchangeabilities, written
lower-triangular, followed by a placeholder uniform frequency line. The pi are deliberately
omitted because the matrix is meant to be reused with a profile mixture that supplies its own.
The citation printed is Banos et al. (2025), *Molecular Biology and Evolution*
(`main/phyloanalysis.cpp:122-129`).

---

## 10. Procedure F: PMSF, posterior-mean site frequencies

Invocation, two stages:

```
iqtree3 -s aln -m LG+C20+F+G -ft guide.treefile      # stage 1, writes <prefix>.sitefreq
iqtree3 -s aln -m LG+G -fs <prefix>.sitefreq          # stage 2, or continue in the same process
```

**Stage 1** (`main/phyloanalysis.cpp:4839-4894`). The topology is **fixed** to the `-ft` tree. A
profile mixture, which must genuinely be a mixture or the run errors out, is fitted for model
parameters and branch lengths only, at ten times the usual tolerance for speed. Then, for each
pattern, the posterior-mean frequency vector is computed
(`PhyloTree::computePatternStateFreq`, `tree/phylotree.cpp:1439-1495`):

```
pihat_s(x) = sum over classes m of  P(class m | D_s, T, thetahat) * F_m(x)

P(m | D_s) = w_m sum_c rho_c L(D_s | m, c)  /  sum_{m'} w_{m'} sum_c rho_c L(D_s | m', c)
```

The conditioning is on the **entire site column** given the whole tree, not on one sequence or
one node. The class weight is already folded into `_pattern_lh_cat` by the kernel, and
`transformPatternLhCat` sums over rate categories within each class. `F_m` are the fixed class
profiles. For `LG+C20+F`, M is 21.

`--freq-max` selects the posterior **mode** instead: the profile of the single highest-posterior
class is copied verbatim rather than averaged (`tree/phylotree.cpp:1476-1493`). Default with
`-ft` is the mean.

**Stage 2** (`model/modelfactory.cpp:643-672`). A `ModelSet` is built with one `ModelMarkov` per
pattern. The exchangeability matrix is taken **once** from the base model and shared by all
patterns; each pattern's pi is its own posterior-mean vector. The mixture is gone: stage 2 of
`LG+C20+F+G` is effectively `LG+G` with one pi per site. `normalize_matrix` stays true, so
**each site's Q is individually normalised to unit expected rate**.

`ModelSet` is a per-site *assignment*, not a mixture: `getPtnModelID(ptn)` returns `ptn`, there
are no weights, `isMixture()` is false, and the kernels index the eigen arrays by pattern under a
`SITE_MODEL` template flag. Memory rises to `nptn * nstates * leafNum` for tip vectors. The
exchangeabilities remain a single shared, optimisable parameter set: `getNDim()` returns
`front()->getNDim()` and `getVariables` broadcasts the same vector to every site model.

**File format.** `<prefix>.sitefreq` is one line per site: a 1-based site index then 20 numbers.
Reading it back (`alignment/alignment.cpp:6151-6295`) enforces strictly increasing site indices
with no duplicates, frequencies strictly inside (0,1), and a sum within 1e-4 of 1 or a
normalisation warning. Unlisted sites keep the model's own frequencies. If two sites sharing a
pattern get different vectors, the alignment's patterns are **split**.

Citation printed: Wang, Susko, Minh and Roger (2018), *Systematic Biology* 67:216-235.

**A command-line collision worth knowing.** `--site-freq` is parsed twice in the same loop, first
as an AliSim option taking `MEAN`, `SAMPLING` or `MODEL`, and only later as the PMSF file option.
The first branch wins, so `--site-freq <filename>` throws
`"Use --site-freq MEAN/SAMPLING/MODEL"`. **Only `-fs <file>` selects a site-frequency file**,
despite the help text advertising both spellings.

---

## 11. Procedure G: MUTSEL, site-specific mutation-selection models

New in v3.1.4, implemented in Rust, **not compiled by default**.

Invocation:

```
iqtree3 -s aln -m MUTSEL -ft guide.treefile
iqtree3 -s aln -m MUTSEL{0.65/1.05/1.86} -ft guide.treefile
iqtree3 -s aln --site-model-file <prefix>.sitemodel        # reuse a previous fit, no Rust needed
```

The guide tree must carry branch lengths. The model-name string is forwarded verbatim to Rust,
which panics on anything other than `MUTSEL` or `MUTSEL{a/b/c}`, so suffixes such as `+G4` are
invalid.

**The mathematics.** Site s has its own stationary distribution `pi^s`. A single mutation process
is shared: a symmetric exchangeability matrix `M` with its own equilibrium `pi^mut`. Site-specific
scaled fitnesses are **derived**, not estimated separately:

```
f^s_i = log pi^s_i - log pi^mut_i
Q^s_ij = M_ij * g( f^s_j - f^s_i ),    g(x) = x / (1 - exp(-x)),    i != j
```

with `g` replaced by its Taylor expansion `1 + x/2 + x^2/12` for `|x| < 1e-6`
(`mutsel_rust/src/model.rs:7-13`). This is the Halpern and Bruno fixation factor, with no
explicit effective population size: the selection coefficients are already scaled. Because M is
reversible with respect to `pi^mut`, the definition of `f^s` makes `pi^s` the stationary
distribution of `Q^s` exactly, and the Rust unit test asserts the equivalent symmetry of
`diag(sqrt(pi^s)) Q^s diag(sqrt(pi^s))^{-1}` to within 1e-12. **MUTSEL is therefore reversible.**
Only 190 numbers per site cross the interface, the upper triangle of `R^s`, so a non-reversible Q
could not be represented even in principle.

**What is estimated, and how.** A shared 190-entry mutation matrix plus its 20-entry equilibrium,
one free log-frequency vector per **site** (not per pattern), a global scaling, and one
re-scaling factor per branch with the geometric mean pinned. The optimiser is **AdamW with
learning rate 0.05 and reverse-mode automatic differentiation**, not BFGS, not EM, and not MCMC
(`mutsel_rust/src/optimization.rs:453-515`). The Felsenstein likelihood and its gradient come
from the external crate `phylo_grad`. Stopping: at least 100 and at most 500 iterations, halting
when the relative improvement is at or below 1e-5 for five consecutive iterations.

**This is a MAP estimate, not an ML estimate.** Three quadratic penalties in log space act as
Gaussian priors, with default strengths 0.65 on the site frequencies, 1.05 on the mutation
matrix, and 1.86 on the branch-length rescaling. The prior centres are: for the mutation matrix,
`--mutsel-prior-rate-file` in PAML lower-triangle format with 20 equilibrium frequencies
appended, defaulting to a built-in codon-structure-derived table; for the site frequencies,
`--mutsel-prior-freq-file`, or, if absent, an internally computed lightweight PMSF over a
built-in UDM256 profile mixture plus one empirical class at weight 0.1. Neither prior option
appears in the help text.

**Rate heterogeneity.** IQ-TREE always requests a single free-rate category, so there is no
across-site rate mixture inside the fit. All across-site rate variation is carried by the
differing per-site `mu^s`.

**Normalisation.** Enforced on the Rust side, not in C++: Q is divided by the arithmetic mean of
the per-site rates so that

```
(1/L) sum over sites s of mu^s = 1,      mu^s = - sum_j pi^s_j Q^s_jj
```

On the C++ side each site model has `normalize_matrix = false`, so per-site normalisation is
explicitly disabled and only this cross-site average holds. **No C++ code re-checks the
invariant**, so a `.sitemodel` file supplied with `--site-model-file` is trusted as given. Because
`mu^s` can greatly exceed 1 in this time unit, `max_branch_length` is forced to 200.0.

The contrast with PMSF is sharp and both use the same `ModelSet` container: PMSF normalises each
site's Q individually to unit rate, MUTSEL normalises only the average across sites.

**Build gating.** `USE_MUTSEL` is tested with `if (USE_MUTSEL STREQUAL "ON")` but is never
declared with `option(...)`, so it is undefined unless the user passes `-DUSE_MUTSEL=ON`. A
default build compiles stubs that print "Mutsel support not compiled in!" and call `exit(1)`.
Continuous integration always enables it.

**Corrections to the existing documentation.** `docs/agent/ARCHITECTURE.md` describes
`mutsel_rust/` as 13 files and roughly 13k lines of Rust. It is **9 files and 10,839 lines**, of
which `data.rs` alone is 7,854. Two sibling models in `model.rs`, `calc_rate_matrix_relaxpmsf`
and `calc_rate_matrix_mutselapprox`, are unreachable from IQ-TREE, which hard-codes
`SubstitutionModel::MutSel`. There is no user documentation of MUTSEL anywhere in the tree, no
named paper reference, and no test coverage of the C++ path; the citation block prints only the
PMSF citation even for a MUTSEL run.

---

## 12. Model selection: which amino-acid model, by ModelFinder

Choosing among models is estimation too, and the choice is made by exhaustive scoring with two
heuristic filters, not by a greedy walk.

### 12.1 The candidate set

`CandidateModelSet::generate` (`main/phylotesting.cpp:1593-1721`) forms a strict three-factor
cross product: matrix, then frequency, then rate heterogeneity. For protein data there are no
exclusion rules at all; the only restrictions in that code block are codon-specific.

**Matrices**: the 28 names of `aa_model_names` by default. `-mset` selects a preset
(`phyml`, `raxml`, `mrbayes`, `beast1`, `beast2`, `modelomatic`, `non-reversible`), a comma list,
or a comma list prefixed with `+` to append to the 28. `-msub` intersects with `nuclear`,
`mitochondrial`, `chloroplast` or `viral`. All names are upper-cased.

**Frequencies**: `{"", "+F"}` by default. `-mfreq FULL` gives
`{"", "+F", "+FQ", "+FO", "C10", ..., "C60"}` and `-mfreq COMPLETE` omits `+FQ` and `+FO`. The
C-series entries become profile-mixture suffixes, and this is the only default-reachable route to
C10 through C60.

**Rate heterogeneity** (`getRateHet`, `main/phylotesting.cpp:1235-1338`). For `-m MFP` or
`-m MF` on a normal alignment the option set is `{"", "+I", "+G", "+I+G", "+R", "+I+R"}`. When
the alignment has **no** invariant sites, `+I` is dropped entirely and the set is
`{"", "+G", "+R"}`. `+R` and `+I+R` are then expanded over category counts from `--cmin` to
`--cmax`, default 2 to 10. `+G` carries no number and means `+G4`. The default protein list has
22 entries.

**The default `-m MFP` candidate count on protein data is therefore 28 x 2 x 22 = 1232**, or
28 x 2 x 11 = 616 when there are no invariant sites. `-madd` appends raw model strings with a
flag that exempts them from every filter, and it is the route by which `GTR20`, `NONREV`,
`LG4X` or `LG+C60+F+G` enter a ModelFinder run.

**Non-reversible models are opt-in and cannot be mixed.** No `NQ.*` name appears in the default
array. `-mset non-reversible` selects the six, and also sets `contain_nonrev`, which is
**mandatory**: a candidate set containing a non-reversible name without `--nonrev-model` is a
hard error, and a set containing both reversible and non-reversible names is a hard error
(`main/phylotesting.cpp:448-489`, `1157-1159`). `NONREV` itself is in no array, so it escapes
that check by name and must be added with `-madd` or `-mset`.

### 12.2 The tree used, and what is re-optimised

`computeFastMLTree` (`main/phylotesting.cpp:740-904`) builds **one fixed topology** used for
every candidate: a parsimony or BIONJ start tree, then `LG` plus the first rate option
(`+I+G` normally, `+G` with no invariant sites) optimised at fifty times the ModelFinder epsilon,
then a fast NNI hill-climb, skipped if a user tree was given. `modelfinder_ml_tree` is true by
default.

Branch lengths **are** re-optimised for every candidate, from the checkpointed tree, with
`BRLEN_OPTIMIZE` and tolerance `modelfinder_eps = 0.1`. `+R` candidates warm-start from the
`+R(k-1)` fit. `-mtree` instead runs a full tree search per candidate.

### 12.3 Scoring

`computeInformationScores` (`main/phylotesting.cpp:632-636`):

```
AIC  = -2 logL + 2 k
AICc = AIC + 2 k (k+1) / max(n - k - 1, 1)
BIC  = -2 logL + k log n
```

The AICc denominator is clamped at 1, so it never goes negative. The default criterion is
**BIC** (`utils/tools.cpp:7444`). The sample size `n` is the number of **sites**, not patterns
(`Alignment::getNSite()` returns `site_pattern.size()`), overridable with `-ms`.

`k` is `ModelFactory::getNParameters` (`model/modelfactory.cpp:1253-1257`):

```
k = model->getNDim() + model->getNDimFreq() + site_rate->getNDim() + tree->getNBranchParameters(brlen_type)
```

**Branch lengths are counted.** `getNBranchParameters` (`tree/phylotree.cpp:2608-2634`) returns
`branchNum` for an unrooted tree, subtracts one for a rooted tree's root pendant edge, and
subtracts a **second** one if the model is reversible, because only the sum of the two edges at
the root affects a reversible likelihood. Multiply by `getNMixlen()` for heterotachy. Returns 0
under `-blfix`, and 1 under `BRLEN_SCALE`.

### 12.4 Search strategy and its two filters

The loop is exhaustive over the generated vector, with three pruning mechanisms.

**`filterRates` is on by default.** `ratehet_set` defaults to `AUTO`, which makes the first
matrix (LG with no `+F`) be tested against **all** 22 rate options; only those scoring within
`--score-diff` of the best are then applied to the remaining 55 matrix-and-frequency
combinations. `--score-diff` defaults to **10.0**, and `--score-diff ALL` or `-m MF1` disables
the filter (`main/phylotesting.cpp:3055-3078`). This is the single most important thing to state
when reporting a ModelFinder result: the search is not exhaustive by default.

**`filterSubst`** is the mirror image and requires `-mset AUTO`, which is not the default.

**The `+R` ladder stops early.** As soon as `+R(k)` scores worse than `+R(k-1)` for the same
matrix, every `+R(j)` with `j > k` for that matrix is skipped
(`main/phylotesting.cpp:3283-3321`, `3434-3445`). The same treatment applies to `+H`.

All three criteria are tracked simultaneously, a ranked list is written, and Akaike weights with
a 0.05 confidence set are computed in `ModelCheckpoint::getOrderedModels`.

### 12.5 Partitioned data

`-m MFP+MERGE` runs PartitionFinder-style greedy merging, defaulting to **fast relaxed
clustering** (`MERGE_RCLUSTERF`). Partitions are first scored independently, candidate pairs are
the closest `--rcluster` percent (default 10, capped at 10 times the partition count) by absolute
difference in tree length, each pair is scored by concatenating and rescaling the subtree to the
geometric mean of the two tree lengths, and up to half the partitions are merged per round in
disjoint pairs. The merge-phase candidate set is deliberately tiny: `merge_models` and
`merge_rates` both default to `"1"`, which for protein means **only `LG+I+G4` and `LG+F+I+G4`**
are evaluated for each trial merge.

**`--model-joint` suppresses ModelFinder** for partitions that would otherwise have had an empty
model (`main/phylotesting.cpp:1380-1386`), because one joint matrix is being estimated instead.

---

## 13. Parameter counts, for the record

For 20 states, with B branch parameters from Section 12.3:

| model string | `getNDim` | `getNDimFreq` | rate df | total df |
|---|---|---|---|---|
| `LG` | 0 | 0 | 0 | B |
| `LG+F` | 0 | 19 | 0 | B + 19 |
| `LG+FO` | 19 | 0 | 0 | B + 19 |
| `LG+FQ` | 0 | 0 | 0 | B |
| `LG+G4` | 0 | 0 | 1 | B + 1 |
| `LG+I+G4` | 0 | 0 | 2 | B + 2 |
| `LG+F+R k` | 0 | 19 | 2k - 1 | B + 19 + 2k - 1 |
| `GTR20` | 189 | 19 | 0 | B + 208 |
| `GTR20+FO` | 208 | 0 | 0 | B + 208 |
| `NQ.pfam` | 0 | 0 | 0 | B' |
| `NONREV` | 379 | 0 | 0 | B' + 379 |
| `LG+C20+F` | 20 weights | 19 (counted once) | 0 | B + 39 |
| `LG4M` | 0 | 0 | 1 | B + 1 |
| `LG4X` | 0 | 0 | 6 | B + 6 |
| `--model-joint NONREV`, K partitions | 379, counted once | 0 | per partition | see Section 8 |

B is the unrooted branch count `2s - 3` for s taxa; B' is the rooted count, one larger, because a
rooted tree loses only one degree of freedom for the root pendant edge when the model is
non-reversible instead of two.

Note that `+F` and `+FO` cost the same 19 degrees of freedom but are booked in different
functions, `getNDimFreq` and `getNDim` respectively. `getNDimFreq` deliberately counts an
empirical frequency vector **once** across mixture classes rather than once per class.

---

## 14. What comes out, and in what format

**The `.iqtree` report.** `reportModel` (`main/phyloanalysis.cpp:581-660`) prints, for a protein
model with more than 20 free dimensions, a warning about overfitting and then the matrix itself:
reversible models as a **lower-triangular list in PAML format**, non-reversible models as the
**full 20x20 Q including the diagonal**, both followed by the 20 state frequencies. For a
non-reversible model the printed frequencies are the solved stationary vector, so the block is
self-consistent and can be fed back with `-m <file>`.

For `--model-joint`, the linked matrix appears once, under `Linked model of substitution:`,
after a temporary `fixParameters(false)` so that the parameter list is not suppressed
(`main/phyloanalysis.cpp:1388-1402`).

**`<prefix>.GTRPMIX.nex`** is written **only** under `--link-exchange`. It holds one NEXUS
`model GTRPMIX = ...;` entry: lower-triangular exchangeabilities for a reversible model, the full
Q for a non-reversible one, and then a **uniform 1/20 frequency line in both cases**. For the
non-reversible branch that placeholder is correct, because pi is already folded into Q and will
be re-derived on reading. For the reversible branch it silently discards the estimated pi, and
the user must re-supply `+F` or `+FO`.

**`<prefix>.sitefreq`** is the PMSF per-site frequency table, one line per site.

**`<prefix>.sitemodel`** is the MUTSEL binary per-site model: a length-prefixed rate-model
string, the site count, then `L * 20` frequencies and `L * 190` upper-triangular rates as
little-endian doubles.

**Checkpointing.** Every model must round-trip its rates through
`startCheckpoint` / `saveCheckpoint` / `restoreCheckpoint`, with the base class called **last on
save and first on restore**, and `decomposeRateMatrix()` followed by `clearAllPartialLH()` after
restore (`model/modelunrest.cpp:113-138` is the reference). A linked model is checkpointed and
dumped to disk after every matrix update, so a long QMaker run resumes at that granularity.

---

## 15. Defects, dead code, and things that silently change results

Found while reading, listed because each one can change a number or waste a week.

1. **`--eigen` is broken for non-reversible models.** `decomposeRateMatrixNonrev` returns early
   after `eigensystem_nonrev`, but `computeTransMatrixNonrev` has no branch for
   `MET_EIGEN_DECOMPOSITION` and falls into `ASSERT(0 && "this line should not be reached")`.
   `ASSERT` is a no-op only under `NDEBUG` (item 9); in a build with `NDEBUG`, `trans_matrix` is
   left untouched.
   `computeTransMatrixEigen`, the function that would consume that output, has zero callers.

2. **`--matrix-exp` does not exist.** An earlier version of `docs/agent/FILE_INDEX.md` named it;
   the index has since been corrected. The real flags are `--eigenlib`, `--eigen`,
   `--scaling-squaring` and `--lie-markov` (`utils/tools.cpp:5093-5108`), and the default is
   `MET_EIGEN3LIB_DECOMPOSITION`.

3. **`--init-model DIVMAT` is guarded by `ASSERT(0 && "init_by_div_mat not working")`.** In a
   build with `NDEBUG` that assertion vanishes and the remaining code runs, feeding
   `setFullRateMatrix` the raw row-normalised divergence matrix instead of its logarithm, and
   applying it only to the representative partition; without `NDEBUG` (item 9) the assertion
   stops the run. Do not use it.

4. **`--site-freq <file>` never reaches the PMSF parser**, because an AliSim branch expecting
   `MEAN`, `SAMPLING` or `MODEL` is tested first in the same loop. Use `-fs`.

5. **`+F` on a non-reversible amino-acid model is silently ignored**, for the reason in
   Section 1.5.

6. **`-m NONREV+F{...}` pins pi and skips the stationarity solve while leaving Q unconstrained.**
   The root distribution used in the likelihood is then, in general, not stationary for the
   fitted Q: a nonstationary-root model. This is a live path, not a bug to be fixed casually.

7. **`luinverse` is called with `num_state` rather than `new_num`**
   (`utils/eigendecomposition.cpp:387`, `159`), so when any state has `pi <= 1e-10` it inverts an
   n-by-n matrix of which only the leading block was filled. Reachable only under `--eigen`.

8. **No clamping and no row-sum renormalisation on the default x86 reversible path.** A negative
   transition probability propagates into the kernel unchanged. Only the ARM NEON branch clamps
   and closes rows.

9. **Every `ASSERT`-based sanity check vanishes under `NDEBUG`** (`utils/tools.h:60-68`),
   including the row-sum checks on P(t), the eigen-equation residual test, and the stationarity
   check on the solved pi. Whether a given build defines `NDEBUG` is a build-configuration fact.
   Correction, 2026-09-23: an earlier version of this item assumed that Release builds define
   it. For gcc and clang, the root `CMakeLists.txt` replaces `CMAKE_CXX_FLAGS_RELEASE` with flags
   that omit `-DNDEBUG` (lines 400 and 420), and nothing else in the build defines it
   (`ncl/ncl.h:72` does so only for the Metrowerks compiler). So in gcc and clang Release builds
   the `ASSERT`s stay active; confirmed from the configure output of a clang 14 Release build on
   2026-09-23, whose flags contained no `-DNDEBUG`.

10. **`tqli` gives up after 100 iterations with a warning and continues** with an unconverged
    eigenvalue (`utils/eigendecomposition.cpp:682-685`); `hqr2` instead calls `exit(1)`.

11. **`ModelMarkov::setRates()` is `ASSERT(0)`.** Any indirect parameterisation must override it
    and must actually call it from `getVariables`.

12. **Dead code that misleads searches.** `model/modelnonrev.{cpp,h}` are zero bytes and unbuilt.
    `model/modelgtr.cpp` is unbuilt and includes a header that does not exist. `cmaple/` is a
    separate engine with its own `NONREV` and `GTR20`. `phylotesting.cpp:2379-3039` is a
    commented-out legacy `testPartitionModel` whose greedy-merge loop reads almost identically to
    the live one. `PartitionModel::reportLinkedModel` has no callers. `Optimization::brent` and
    `Optimization::dbrent` are declared but undefined, so calling either is a link error.

13. **`mutsel_rust/` is 9 files and 10,839 lines**, not the 13 files and roughly 13k lines stated
    in `docs/agent/ARCHITECTURE.md`. That document's "arrived upstream and not yet reviewed"
    section should now be updated: the subproject has been read, and Section 11 records what it
    does.

---

## 16. Verification status

**Read directly and quoted in this document:** `model/modelmarkov.cpp` (constructor,
`setReversible`, `init_state_freq`, `getRateMatrix` through `getNDimFreq`, `setVariables`,
`getVariables`, `targetFunk`, `setBounds`, `optimizeParameters`, `decomposeRateMatrixNonrev`,
`decomposeRateMatrix`, `decomposeRateMatrixRev`, `readRates`, `readParameters`,
`computeTransMatrix`, `computeStateFreqFromQMatrix`); `model/modelmarkov.h`;
`model/modelprotein.cpp` (`init`, `rescaleRates`, `readRates`, `getNameParams`, the NEXUS block
format); `model/modelunrest.cpp` in full; `model/modelfactory.cpp` (frequency parsing, the
`ModelSet` branches, `getNParameters`, `optimizeParametersOnly`, `optimizeAllParameters`,
`optimizeParameters`, the fused-rate conditions); `model/modelmixture.cpp` (`createModel`,
`setBounds`, `computeTransMatrix`, the mixture rate assignment and rescaling);
`utils/eigendecomposition.cpp` (`computeRateMatrix`, `eliminateZero`, `symmetrizeRateMatrix`);
`utils/optimization.cpp` (`minimizeMultiDimen`, `restartParameters`, `derivativeFunk`);
`alignment/alignment.cpp` (`computeStateFreq`, `countStates`, `convertCountToFreq`, `convfreq`,
`getAppearance`); `tree/phylotree.cpp` (`getNBranchParameters`, `optimizeRootPosition`);
`tree/iqtree.cpp` (`initializeModel`, `optimizeModelParameters`, `computeInitialTree`, the NNI
re-optimisation block); `tree/phylotreesse.cpp` and `tree/phylokernelnew.h` (tip-vector
preparation and the root-leaf frequency fill); `main/phyloanalysis.cpp` (`reportModel`,
`reportNexusFile` and its call sites); `utils/tools.h` and `utils/tools.cpp` (the `StateFreqType`
enum, defaults, and option parsing); `mutsel_rust/src/model.rs` (the fixation function);
`CMakeLists.txt` (the `USE_MUTSEL` gate).

**Established by delegated reading and spot-checked but not read line by line by the author of
this document:** the internals of `dfpmin` and `lnsrch`; the EISPACK chain in
`eigensystem_nonrev`; the SIMD likelihood kernels beyond the passages quoted; `RateGamma`,
`RateInvar` and `RateFree`; `partitionmodel.cpp` and `partitionmodelplen.cpp` beyond
`getNParameters` and the factory selection; `phylotesting.cpp` beyond the model-name arrays;
`modelmixture.cpp` beyond the functions listed above; the Rust optimiser and prior machinery.
Two claims in those reports that contradicted existing project documentation were checked
directly and confirmed: the absence of `--matrix-exp`, and the size of `mutsel_rust/`. The
forward-difference gradient and the protein-specific `bound_check` override in
`ModelMixture::setBounds` were also checked directly and confirmed.

**Derivations, not code.** The statement that `pi^T Q = 0` has rank `n - 1` given zero row sums,
the resulting 360-parameter count for a pi-constrained non-reversible amino-acid model, and the
verification that the Halpern and Bruno fitness definition makes `pi^s` stationary for the MUTSEL
`Q^s`, are derivations presented as such. The first two have not been independently checked; the
third is corroborated by a unit test in the Rust source.

**Not verified.** No claim here rests on running IQ-TREE. Nothing in this document was checked
against numerical output, because the build is not currently working on this machine.

---

## Quick reference

| question | answer | anchor |
|---|---|---|
| Objective for a single alignment | full-tree log-likelihood at current branch lengths | `model/modelmarkov.cpp:1091` |
| Objective for many alignments | unweighted sum of per-partition log-likelihoods | `model/partitionmodel.cpp:299` |
| Optimiser | BFGS (`dfpmin`), 200 iterations, `gtol = 1e-4` | `utils/optimization.cpp:781` |
| Gradients | numerical forward difference, step 1e-4 times the parameter | `utils/optimization.cpp:916` |
| Constraints | box clamping of the trial point only | `utils/optimization.cpp:149` |
| Rate bounds | `[1e-4, 100]` | `model/modelmarkov.h:30` |
| Restarts | up to 3, random, only when `bound_check` is set | `utils/optimization.cpp:726` |
| Reversible free params | 189 exchangeabilities, plus 19 if `+FO` | `model/modelprotein.cpp:1205` |
| Non-reversible free params | 379 directed rates, pi derived | `model/modelprotein.cpp:1233` |
| Q normalisation | `sum_i pi_i (-Q_ii) = 1` | `model/modelmarkov.cpp:1275` |
| pi for reversible | input: matrix, `+F`, `+FO`, `+FQ` | `model/modelmarkov.cpp:299` |
| pi for non-reversible | solved from Q every evaluation | `model/modelmarkov.cpp:1267` |
| The Q-to-pi solver | column-pivoted QR on `[1^T; Q^T] pi = e_1` | `model/modelmarkov.cpp:2118` |
| Reversible eigen path | `SelfAdjointEigenSolver` on `Pi^(1/2) Q Pi^(-1/2)` | `model/modelmarkov.cpp:1519` |
| Non-reversible eigen path | complex `EigenSolver`, `FullPivLU` inverse | `model/modelmarkov.cpp:1324` |
| Matrix exponential fallback | Eigen scaling and squaring, Pade 13 | `model/modelmarkov.cpp:469` |
| Branch lengths | safeguarded Newton on the analytic first derivative | `tree/phylotree.cpp:2635` |
| Root position | local search over branches within distance 2 | `tree/phylotree.cpp:3072` |
| Model selection criterion | BIC, `n` = number of sites | `main/phylotesting.cpp:632` |
| Default candidate count, protein | 1232, filtered by `--score-diff 10` | `main/phylotesting.cpp:1593` |
| Empirical pi | 8-round fixed point over ambiguity codes | `alignment/alignment.cpp` |
| Outer alternation | 100 rounds, branch lengths then model then rates | `model/modelfactory.cpp:1570` |
