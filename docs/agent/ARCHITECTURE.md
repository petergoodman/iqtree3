# IQ-TREE 3 — Architecture Map for Agents

**Audience:** AI coding agents extending IQ-TREE 3's substitution-model layer, specifically for
constrained optimization of non-reversible amino-acid rate matrices (nQ) under a fixed target
stationary distribution π.

**Purpose:** give an agent, in one read, enough of IQ-TREE's structure and conventions to write
code that looks like the code already here — without re-deriving it from a thousand greps.

**Companion:** `AGENT_CONTEXT_FILE_INDEX.md` — the per-file index with line anchors. Read this
document first, then use that one to locate specific code.

## Baseline and staleness

Verified against local commit `8977d31a`. Line numbers drift, names and structure do not.

At the time of writing, this working tree sits 64 commits behind `origin/master`
(`63c330d9`), and those commits rewrite several files described here:

| File | Upstream churn | Status of this document |
|---|---|---|
| `model/modelmarkov.{cpp,h}` | none | anchors valid |
| `model/modelprotein.cpp` | none | anchors valid |
| `model/modelunrest.cpp`, `model/modelliemarkov.cpp`, `model/modelmixture.cpp` | none | anchors valid |
| `main/phyloanalysis.cpp` | 377 lines | anchors likely stale |
| `model/partitionmodel.cpp` | 155 lines | anchors likely stale |
| `main/phylotesting.cpp` | 135 lines | anchors likely stale |
| `model/modelfactory.cpp` | 132 lines | anchors likely stale |
| `tree/phylotree.cpp` | 126 lines | anchors likely stale |
| `utils/tools.{cpp,h}` | 93 lines | anchors likely stale |
| `model/modelfactorymixlen.{cpp,h}` | deleted upstream | remove references when rebasing |

The entire Tier 0 core is untouched upstream, so the substantive findings hold. To
re-check any single file before trusting its anchors:

```bash
git diff --stat 8977d31a HEAD -- <path>
```

After changing the baseline, update this block and re-verify anchors in the affected files.

---

## 1. The 30-second version

IQ-TREE is a ~500k-line C++17 codebase, organized as a dozen static libraries linked into one
executable. For substitution-model work, only four of those libraries matter:

```
utils/      Params (global CLI state), Optimization (BFGS), EigenDecomposition, Checkpoint
alignment/  Alignment: patterns, state counts, empirical frequencies
model/      ModelSubst hierarchy (Q matrices) + RateHeterogeneity hierarchy (site rates)
tree/       PhyloTree / IQTree: topology, branch lengths, likelihood kernels
main/       phyloanalysis.cpp (driver + reports), phylotesting.cpp (ModelFinder)
```

Everything optimizable — a substitution model, a rate model, a tree — inherits the same
`Optimization` interface and is driven by the same BFGS routine. **Learn that one interface and
you can extend any of them.** It is described in §5.

---

## 2. Build topology

Root `CMakeLists.txt` adds each source directory as a static library (lines 840–893), then links
them into the `iqtree3` executable (line 1019):

```
iqtree3  ←  main  model  tree  alignment  utils  nclextra  ncl  pll  pda
            lbfgsb  whtest  sprng  vectorclass  gsl  simulator  yaml-cpp
            phyloYAML  kernelsse [kernelavx kernelfma kernelavx512]
```

Facts worth knowing before you add a file:

- `model/CMakeLists.txt` is a **hand-maintained explicit source list**. There is no glob. A new
  `model/*.cpp` that is not added there simply never compiles, and you get link errors that look
  like missing symbols rather than a missing file.
- `target_link_libraries(model utils)` — the `model` library may depend on `utils` and on headers
  from `tree/` and `alignment/` (include paths are global), but not on `main/`.
- The likelihood kernels are compiled **multiple times** into separate libraries
  (`kernelsse`, `kernelavx`, `kernelfma`, `kernelavx512`) with different instruction-set flags.
  Anything included by `tree/phylokernel*.h` gets compiled once per ISA. Keep model headers out
  of the kernel include path if you can.
- Hard external dependencies: **Eigen3** and **Boost** (`find_package`, fatal if absent). Eigen is
  used directly inside `model/modelmarkov.cpp` — use it for new linear algebra.
- `iqtree3-aa` is a CMAPLE-AA build variant, unrelated to amino-acid model work despite the name.

### Dead files — do not be misled

- `model/modelnonrev.cpp` and `model/modelnonrev.h` are **0 bytes** and are not in
  `model/CMakeLists.txt`. Non-reversible logic lives in `modelmarkov.cpp` behind the
  `is_reversible` flag, plus `modelunrest.*` and `modelliemarkov.*`.
- `model/modelgtr.cpp` (887 lines) is **not** in `model/CMakeLists.txt`, and it `#include`s
  `modelgtr.h`, **a header that does not exist in the repository**. It therefore cannot compile
  even if added to the build. It is an unbuilt historical copy of the reversible code that now
  lives in `modelmarkov.cpp`, and reading it will teach you things that are no longer true.
- `cmaple/` contains a complete **second** model hierarchy with its own `NONREV`, `GTR20`, and
  amino-acid handling. It is a different engine. Greps constantly land there. Ignore it.

### There is no non-reversible module, only a flag

Non-reversible inference works and is fully supported, but it has no home of its own. There is
no `ModelNonrev` class and no non-reversible source file. The capability is **98
`is_reversible` / `isReversible()` / `useRevKernel()` branch sites** spread across 25 live
files:

```
30  model/modelmarkov.cpp        7  main/phyloanalysis.cpp     4  tree/phylotreesse.cpp
10  model/modelpomo.cpp          6  model/modelmixture.cpp     4  model/modelliemarkov.cpp
 7  model/partitionmodel.cpp     3  model/modelprotein.cpp     ... 16 more files
```

`NONREV` is likewise **a name string, not a class.** The class is `ModelProtein`, which handles
every amino-acid model through an `if/else` chain over the name, and `NONREV` is one branch of
it at `modelprotein.cpp:1207`. `GTR20` is the identical construction with `is_reversible` left
true. By contrast `UNREST` **is** a class (`ModelUnrest`), registered through
`ModelMarkov::getModelByName`. Both idioms are available when adding a model.

The three shipped entry points for non-reversible amino-acid inference are `-m NONREV` (fit 379
free rates to one alignment), `-m NQ.pfam` and siblings (use a pre-estimated matrix), and
`--model-joint NONREV` (estimate one matrix across many alignments, the nQMaker workflow).

---

## 3. Runtime control flow

```
main()                                   main/main.cpp:2219
 └─ parseArg(argc, argv, Params)         utils/tools.cpp  — fills the Params singleton
 └─ runPhyloAnalysis(params, checkpoint) main/phyloanalysis.cpp
     ├─ read alignment                   → Alignment (patterns, num_states, empirical π)
     ├─ readModelsDefinition(params)      model/modelfactory.cpp:87  → ModelsBlock
     │     builtin_mixmodels_definition + builtin_prot_models + optional --mdef file
     ├─ runModelFinder (optional)        main/phylotesting.cpp  — picks the -m string
     ├─ iqtree->initializeModel(...)     tree/iqtree.cpp:1063
     │     └─ new ModelFactory / PartitionModel / PartitionModelPlen
     │           └─ createModel(model_str, models_block, freq_type, freq_params, tree)
     │                 model/modelmixture.cpp:3122   ← THE model dispatcher
     │                 └─ new ModelProtein(...)  →  ModelProtein::init()
     │                       modelprotein.cpp:1106  →  ModelMarkov::init(freq)
     │                             →  init_state_freq()  →  decomposeRateMatrix()
     ├─ tree search / optimization
     │     IQTree::optimizeModelParameters()     tree/iqtree.cpp:2307
     │       └─ ModelFactory::optimizeParameters()   modelfactory.cpp:1553
     │            ├─ model->optimizeParameters(eps)      ← ModelMarkov::optimizeParameters
     │            ├─ site_rate->optimizeParameters(eps)
     │            └─ tree->optimizeAllBranches()
     └─ reportPhyloAnalysis(...)         main/phyloanalysis.cpp:1407  → writes .iqtree
```

For partitioned (`-p` / `-q` / `-Q`) analyses the same flow runs with `PhyloSuperTree` and
`PartitionModel`, and an extra step — `PartitionModel::optimizeLinkedModels()` — optimizes a
single Q jointly across all partitions. **That is the QMaker / nQMaker path.**

---

## 4. The object graph

```
Params (singleton)  ── read by everything, everywhere, via Params::getInstance()
    │
Alignment ────────── patterns, num_states, seq_type, empirical frequencies
    │
PhyloTree (IQTree, PhyloSuperTree)
    ├── ModelFactory *model_factory      owns and drives the pair below
    │     ├── ModelSubst *model          the Q matrix
    │     └── RateHeterogeneity *site_rate   +I, +G, +R, +H
    ├── ModelSubst *model                (same pointer, cached)
    └── RateHeterogeneity *site_rate     (same pointer, cached)
```

`PhyloTree::getModel()` and `getModelFactory()` are the accessors (`tree/phylotree.h:527`, `531`).
A model holds a back-pointer `PhyloTree *phylo_tree` — the cycle is intentional, because
`targetFunk()` needs to ask the tree for a likelihood.

### Class hierarchies

```
Optimization                       utils/optimization.h
  ├── ModelSubst                   model/modelsubst.h      (also CheckpointFactory)
  │     └── ModelMarkov            model/modelmarkov.h     (also EigenDecomposition)
  │           ├── ModelBIN
  │           ├── ModelDNA ── ModelDNAError
  │           ├── ModelProtein                 ← amino acids: LG, WAG, Q.*, GTR20, NONREV, NQ.*
  │           ├── ModelCodon ── {Empirical, Parametric, SemiEmpirical}
  │           ├── ModelMorphology
  │           ├── ModelUnrest                  ← general non-reversible (UNREST)
  │           ├── ModelLieMarkov               ← constrained-subspace non-reversible
  │           ├── ModelPoMo ── ModelPoMoMixture
  │           ├── ModelMixture (virtual) ── ModelCodonMixture
  │           └── ModelSet                     ← site-specific models
  ├── RateHeterogeneity            model/rateheterogeneity.h
  │     └── RateGamma, RateInvar, RateFree, RateHeterotachy, RateKategory, ...
  ├── ModelFactory                 model/modelfactory.h    (also CheckpointFactory)
  │     ├── PartitionModel ── PartitionModelPlen
  │     └── ModelFactoryMixlen
  └── PhyloTree                    tree/phylotree.h        (also MTree, CheckpointFactory)
        └── IQTree ── PhyloSuperTree ── PhyloSuperTreePlen / PhyloSuperTreeUnlinked
```

**Note two things:** `ModelMixture` and `ModelSet` inherit from `ModelMarkov` *and* from
`vector<ModelMarkov*>` — a container that is also a model. And `PartitionModel` is a
`ModelFactory`, not a `ModelSubst`; it implements the same `Optimization` interface at the
partition level.

---

## 5. The Optimization contract — the pattern to follow

This is the single most important pattern in the codebase. Every optimizable object implements
it, and BFGS is the only solver.

```cpp
class Optimization {                                     // utils/optimization.h
    virtual int    getNDim();                            // number of free parameters
    virtual double targetFunk(double x[]);               // function to MINIMIZE
    virtual void   setVariables(double *variables);      // model state  -> x
    virtual bool   getVariables(double *variables);      // x -> model state; returns "changed"
    virtual void   setBounds(double *lo, double *hi, bool *bound_check);
    virtual bool   restartParameters(...);               // escape boundary optima (optional)

    double minimizeMultiDimen(guess, ndim, lower, upper, bound_check, gtol, hessian=nullptr);
    double L_BFGS_B(nvar, vars, lower, upper, pgtol, maxit);
};
```

### Hard conventions

1. **All parameter vectors are 1-indexed.** `variables[1] .. variables[ndim]`. `variables[0]` is
   allocated but unused. Every `memcpy` in this codebase is `variables+1`. Off-by-one here is the
   single most common bug when extending a model.
2. **`targetFunk` returns the NEGATIVE log-likelihood.** Optimizers minimize. `optimizeParameters`
   then returns `-score` to its caller, which expects a log-likelihood.
3. **`getVariables` must return whether anything changed**, and callers use it to skip expensive
   recomputation. Returning `true` unconditionally is correct but slow; returning `false`
   incorrectly produces a silently stale likelihood.
4. **After any change to model parameters you must, in order:** `decomposeRateMatrix()`, then
   `phylo_tree->clearAllPartialLH()`, then recompute the likelihood. `targetFunk` does exactly
   this (`modelmarkov.cpp:1091–1118`). Skipping `clearAllPartialLH` gives wrong numbers with no
   warning.
5. **Constraints are box constraints only.** `setBounds` fills `lower[]`, `upper[]`, and a
   `bound_check[]` flag that asks the optimizer to retry from a different start if the optimum
   lands on a boundary. There is **no equality-constraint solver in this codebase.** A hard
   constraint has to be expressed by reparameterization (the Lie-Markov approach), by projection
   inside `setRates()`, or by a penalty returned from `targetFunk`.
6. **Invalid regions are signalled by returning `1.0e+30` from `targetFunk`**, not by throwing.
   See the π-underflow guard at `modelmarkov.cpp:1104–1110`.

### The canonical `optimizeParameters` skeleton

`ModelMarkov::optimizeParameters` (`modelmarkov.cpp:1166–1244`) is the reference implementation.
Copy its shape:

```cpp
double MyModel::optimizeParameters(double gradient_epsilon) {
    if (fixed_parameters) return 0.0;
    int ndim = getNDim();
    if (ndim == 0) return 0.0;

    double *variables   = new double[ndim+1];   // 1-indexed
    double *lower_bound = new double[ndim+1];
    double *upper_bound = new double[ndim+1];
    bool   *bound_check = new bool[ndim+1];

    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);

    double score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound,
                                       bound_check, max(gradient_epsilon, TOL_RATE));

    bool changed = getVariables(variables);
    if (changed || score == -1.0e+30) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
    }
    delete [] bound_check; delete [] lower_bound; delete [] upper_bound; delete [] variables;
    return score;
}
```

### Two parameterization styles

**Direct** (`ModelMarkov`, `ModelUnrest`, `ModelProtein`): the BFGS variables *are* entries of
`rates[]`. `setVariables` / `getVariables` are memcpys. `setRates()` is trivial or unused.

**Indirect** (`ModelLieMarkov`, `ModelDNA` with a `param_spec`): the BFGS variables are an
abstract `model_parameters[]` of lower dimension, and **`setRates()` maps them onto the full
`rates[]` array**. `getVariables` writes `model_parameters` and then calls `setRates()`.

> **A constrained nQ belongs to the indirect family.** Fixing π removes degrees of freedom from
> the 380-entry rate space, so the free parameters should live in the reduced space and
> `setRates()` should be the map back into the full rate matrix. `ModelLieMarkov` is the working
> example of exactly that structure:
> `setBounds` (861→886), `setVariables` (899), `getVariables` (919), `setBasis` (1087),
> `setRates` (1194).

---

## 6. Model lifecycle: from a `-m` string to a likelihood

```
"NONREV+F+G4"
  │
  ├─ ModelFactory ctor         strips "+G4" → RateGamma; strips "+F" → StateFreqType
  │                            model/modelfactory.cpp ~150–700
  ├─ createModel("NONREV", models_block, freq_type, freq_params, tree)
  │     model/modelmixture.cpp:3122
  │     ├─ parses NAME{params}, +FQ, +F{...}, +P (PoMo), +E (seq. error)
  │     ├─ if ModelMarkov::validModelName(str)  → UNREST or Lie-Markov
  │     └─ else by seq_type                     → ModelProtein for SEQ_PROTEIN
  │
  ├─ ModelProtein::init("NONREV", ...)          model/modelprotein.cpp:1207
  │     ├─ look up name in ModelsBlock (built-in NEXUS matrices)
  │     ├─ NONREV branch: freq default FREQ_ESTIMATE; seed Q from LG (or --init-model)
  │     ├─ setReversible(false)   → reallocates rates[] to n(n-1)=380, ROOTS the tree
  │     └─ num_params = getNumRateEntries() - 1 = 379
  │
  └─ ModelMarkov::init(freq)                    model/modelmarkov.cpp:354
        ├─ init_state_freq(freq)                sets state_freq per +F / +FO / +FQ / +FU
        └─ decomposeRateMatrix()
              └─ decomposeRateMatrixNonrev()    model/modelmarkov.cpp:1246
                    ├─ unpack rates[] into full rate_matrix, diagonal = -rowsum
                    ├─ computeStateFreqFromQMatrix(rate_matrix, state_freq, n)   ← π SOLVED HERE
                    ├─ normalize so Σ πᵢ(−Qᵢᵢ) = total_num_subst
                    └─ eigendecompose (complex spectrum)
```

From then on, `computeTransMatrix(t, P)` exponentiates via the cached eigensystem, and the
likelihood kernel consumes `P` plus `state_freq` at the root.

---

## 7. Data layouts — get these exactly right

### `double *rates` — the packed rate parameters

| Mode | Length | Packing |
|---|---|---|
| Reversible | `n(n-1)/2` (AA: 190) | Upper triangle, row-major: `(0,1) (0,2) … (0,n-1) (1,2) …`. For DNA: A-C, A-G, A-T, C-G, C-T, G-T. |
| Non-reversible | `n(n-1)` (AA: **380**) | **All off-diagonal entries of Q, row-major, skipping the diagonal:** `(0,1) (0,2) … (0,n-1) (1,0) (1,2) … (n-1,n-2)`. |

The non-reversible packing is defined by the unpacking loop in `decomposeRateMatrixNonrev`
(`modelmarkov.cpp:1256–1263`) and mirrored in `setFullRateMatrix` (`858–873`):

```cpp
for (i = 0, k = 0; i < num_states; i++) {
    double *rate_row = rate_matrix + i*num_states;
    double row_sum = 0.0;
    for (j = 0; j < num_states; j++)
        if (j != i) row_sum += (rate_row[j] = rates[k++]);
    rate_row[i] = -row_sum;          // diagonal is derived, never stored in rates[]
}
```

**File-format trap:** built-in *protein* matrices in NEXUS are written as a **lower** triangle,
while `rates[]` is an **upper** triangle. `ModelProtein::readRates` (`modelprotein.cpp:1265+`)
performs the transposition with
`id = col*(2*num_states-col-1)/2 + (row-col-1)`. Non-reversible protein matrices
(`NQ.*`) are instead written as a **full 20×20 matrix with negative diagonal**, and
`ModelMarkov::readRates` (`1665–1698`) validates that each row sums to zero.

**Reversibility is auto-detected from the sign of the first number in the file**
(`readParameters`, `modelmarkov.cpp:1798`): negative ⇒ full Q ⇒ non-reversible.

### The other arrays

| Field | Size | Meaning |
|---|---|---|
| `rate_matrix` | `n*n` | Full Q, row-major, normalized. Non-reversible only. |
| `state_freq` (in `ModelSubst`) | `n` | π. **For non-reversible models this is an OUTPUT derived from Q**, not an input. |
| `eigenvalues`, `eigenvalues_imag` | `n` | Real and imaginary parts (non-reversible spectra are complex). |
| `eigenvectors`, `inv_eigenvectors`, `inv_eigenvectors_transposed` | `n*n` | Cached eigensystem. |
| `ceval`, `cevec`, `cinv_evec` | `n`, `n*n`, `n*n` | `std::complex<double>` *aliases over the same memory* as the real arrays. Non-reversible only. |

All of these are allocated with `aligned_alloc<double>` / `ensure_aligned_allocated` for SIMD.
Use the same helpers; a plain `new double[]` will crash the AVX kernels.

### Counting parameters

- `num_params` = number of free rate parameters. For `NONREV` on amino acids: `380 - 1 = 379`.
  The last entry `rates[379]` is left at its seed value and acts as the scale anchor; global
  scale is fixed afterwards by the `total_num_subst` normalization.
- `getNDim()` = number of BFGS dimensions. **Non-reversible: `num_params`, with no frequency
  dimensions** (`modelmarkov.cpp:964–976`), because π is determined by Q.
- `getNDimFreq()` = degrees of freedom *not* in `getNDim()`, added by
  `ModelFactory::getNParameters` for AIC/BIC. Returns 0 for `FREQ_ESTIMATE`, `n-1` for
  `FREQ_EMPIRICAL`.

> If a constrained nQ removes `n-1` degrees of freedom from Q and gains none in π, **both
> `getNDim()` and `getNDimFreq()` must be updated**, or every information criterion IQ-TREE
> reports will be wrong.

---

## 8. Reversible vs non-reversible — the complete difference

| Aspect | Reversible (Q) | Non-reversible (nQ) |
|---|---|---|
| Parameterization | R (exchangeabilities) ⊗ π | Q directly |
| `rates[]` length | `n(n-1)/2` | `n(n-1)` |
| π | A free parameter (or empirical/equal) | **Solved from Q** by `computeStateFreqFromQMatrix` |
| `getNDim()` | `num_params + (n-1)` if `FREQ_ESTIMATE` | `num_params` |
| Tree | Unrooted (`convertToUnrooted()`) | **Rooted** (`convertToRooted()`), forced in `setReversible` |
| Eigen spectrum | Real, symmetric after `diag(√π)` similarity | **Complex**; may be non-diagonalizable |
| Decomposition | `SelfAdjointEigenSolver`, or `eigensystem_sym` | `EigenSolver` / `eigensystem_nonrev`; fallback to scaled squaring when `nondiagonalizable` |
| Likelihood kernel | `computePartialLikelihood*` | `computeNonrevPartialLikelihood*` (`tree/phylokernelnonrev.h`); selected by `ModelSubst::useRevKernel()` |
| `ignore_state_freq` | false | **true** — suppresses the `Q *= diag(π)` step in the decomposition |
| Cost | baseline | ~2× memory and time; complex arithmetic |

`setReversible(bool, bool adapt_tree)` (`modelmarkov.cpp:75–152`) performs the whole switch,
including converting an existing reversible half-matrix into a full non-reversible one via
`Qᵢⱼ = Rᵢⱼ πⱼ`. **It also re-roots the tree as a side effect.** That is how `-m NONREV` ends up
with a rooted tree without anyone asking for one.

---

## 9. Where π comes from — the exact code path

This is the crux of the project, so here it is precisely.

**For reversible models**, π is an independent parameter. `init_state_freq()` sets it from `+F`
(empirical), `+FQ` (equal), `+FU` (matrix-supplied), or `+FO` (`FREQ_ESTIMATE`, optimized as
`n-1` extra BFGS dimensions). Detailed balance guarantees *any* positive π is stationary for a
symmetric R — which is exactly why profile mixture models (C10–C60, `+Fmix`) can pair one R with
many π vectors.

**For non-reversible models**, π is derived. Inside `decomposeRateMatrixNonrev`
(`modelmarkov.cpp:1246`):

```cpp
// 1. initialize π to uniform, unless the user pinned it
if (freq_type != FREQ_USER_DEFINED || Params::getInstance().optimize_from_given_params)
    for (i = 0; i < num_states; i++) state_freq[i] = 1.0/num_states;

// 2. unpack rates[] into the full Q with derived diagonal
...

// 3. SOLVE for the stationary distribution
if (freq_type != FREQ_USER_DEFINED || Params::getInstance().optimize_from_given_params)
    computeStateFreqFromQMatrix(rate_matrix, state_freq, num_states);

// 4. normalize Q to the target substitution rate
for (i = 0, sum = 0.0; i < num_states; i++)
    sum -= rate_matrix[i*num_states+i] * state_freq[i];
double delta = total_num_subst / sum;
// ... Q *= delta

// 5. eigendecompose
```

and `computeStateFreqFromQMatrix` (`modelmarkov.cpp:2119–2131`) is:

```cpp
int computeStateFreqFromQMatrix(double Q[], double pi[], int n) {
    MatrixXd A(n+1, n);
    A.topRows(1).setOnes();                          // Σ πᵢ = 1
    A.bottomRows(n) = Map<MatrixXd>(Q, n, n);        // Qᵀ π = 0  (column-major Map transposes)
    VectorXd b(n+1); b.setZero(); b(0) = 1.0;
    Map<VectorXd> freq(pi, n);
    freq = A.colPivHouseholderQr().solve(b);         // least-squares solve of the overdetermined system
    ASSERT(fabs(freq.sum()-1.0) < 1e-4);
    return 0;
}
```

**π is then used as the root distribution in the likelihood.** `tree/phylokernelnew.h:963` fills
the root tip's partial-likelihood vector from `model->getStateFrequency(...)`. So "π must be the
resulting stationary frequency" is equivalent to "the root distribution the likelihood uses must
equal the supplied vector, and must genuinely be stationary for the fitted Q."

**Three existing partial mechanisms, none of which is the requested feature:**

1. `freq_type == FREQ_USER_DEFINED` with `optimize_from_given_params == false` **skips** the
   solve, leaving `state_freq` at the user's values. This *asserts* a π without constraining Q,
   so the model becomes internally inconsistent: the root distribution is not stationary for Q.
2. `ModelMarkov::adaptStateFrequency` (`modelmarkov.cpp:914–931`) rescales the non-reversible
   `rates[]` by `rates[k] *= freq[j]/state_freq[j]`. A one-shot re-weighting, not a constraint
   maintained through optimization.
3. `readParameters` / `readParametersString` (`1798`, `1846`) re-derive π from a file-supplied Q
   and **print a warning** when it disagrees with the file's π by more than `1e-3`. That warning
   is the weak, after-the-fact form of the invariant this project wants to make structural.

### How many parameters a π constraint costs

Stationarity, πᵀQ = 0, is a set of **linear** equations in Q's entries. Because the rows of Q
already sum to zero, those `n` equations sum to zero identically for any π
(Σⱼ Σᵢ πᵢ Qᵢⱼ = Σᵢ πᵢ · 0 = 0), so their rank is `n-1`, not `n`. For 20 states:

| Quantity | Count |
|---|---|
| Off-diagonal entries of Q (the diagonal is derived) | 380 |
| Minus global scale, which `total_num_subst` normalization removes | 379 (this is the current `NONREV` free-parameter count) |
| Minus the rank of πᵀQ = 0 | **360** |

Fixing π therefore costs exactly `n-1 = 19` degrees of freedom, the dimension of a frequency
vector, which is the arithmetic you would expect. The feasible set is a convex polyhedral cone:
a 361-dimensional linear subspace intersected with non-negativity on the off-diagonals. The two
models are nested (360 ⊂ 379), so a likelihood-ratio test between them is available.

### A working precedent: π-constrained non-reversible estimation already exists for DNA

`ModelLieMarkov::setBasis()` (`modelliemarkov.cpp:1087`) already implements "given π, estimate
a non-reversible Q having π as its stationary distribution", for 4-state DNA. When a Lie-Markov
model receives a user-supplied or empirical π:

```cpp
int bdf = BDF[model_num];
num_params = MODEL_PARAMS[model_num] - bdf;   // free parameters REDUCED by the frequency dof
init_state_freq(getFreqType());
double tau[3];
piToTau(state_freq, tau, symmetry);           // π mapped into frequency coordinates
...
unpermuted_rates[rate] += tau[tauIndex] * transformationMatrix[rate];   // basis shifted by τ
```

The comment above the transform tables (`modelliemarkov.cpp:262`) states the intent directly:
*"Each shows how to modify a basis matrix to enforce a fixed base frequency vector."* There is
even a feasibility check that warns when the requested π is unreachable for that model family:
*"Model %s cannot achieve requested equilibrium base frequencies ... Instead it will use ..."*.

What transfers to the 20-state problem:

- The class shape. A reduced `model_parameters[]` array, with `setRates()` mapping it into the
  full `rates[]`, and `setVariables` / `getVariables` / `setBounds` / `getNDim` overridden
  around it. Everything downstream (`decomposeRateMatrixNonrev`, the eigendecomposition, the
  non-reversible kernel, checkpointing) needs no change.
- Setting `freq_type = FREQ_USER_DEFINED` then makes `decomposeRateMatrixNonrev` skip the π
  re-solve automatically, and unlike the current behaviour that skip is **sound**, because the
  parameterization guarantees the answer.
- The non-negativity trick in `setRates()` (`modelliemarkov.cpp:1194`): take the one basis
  matrix whose off-diagonals are all non-negative as an anchor, compute the worst-case ratio
  across entries, and rescale the deviation so nothing crosses zero. This generalizes given any
  strictly positive anchor matrix having the target π as its stationary distribution.

What does not transfer: the `TRANSFORM_*` tables are hand-derived algebra for 4 states and 37
named model families. At 20 states the constraint subspace has to be computed numerically,
which Eigen already supports and which `modelmarkov.cpp` already uses elsewhere. Lie-Markov's
other constraint, closure under the Lie bracket, is a separate restriction that a general nQ
does not want.

---

## 10. Normalization and scaling conventions

- `total_num_subst` (declared in `EigenDecomposition`, `utils/eigendecomposition.h:81`) is the
  target expected substitution rate. Default 1.0; mixture components set it per class
  (`modelmixture.cpp:3446`, `3519`).
- Q is scaled so that `Σᵢ πᵢ (−Qᵢᵢ) = total_num_subst` — one expected substitution per unit
  branch length. **This makes the overall scale of `rates[]` unidentifiable**, which is why
  `num_params = getNumRateEntries() - 1` and the last rate is pinned.
- Branch lengths enter as `evol_time = time / total_num_subst` in `computeTransMatrix`.
- `ZERO_FREQ = 1e-10` (`eigendecomposition.h:24`): states with π below this are **removed** from
  the matrix before decomposition and reinserted afterwards as identity rows
  (`modelmarkov.cpp:1293–1370`). Any constrained-π implementation must handle a target π with
  near-zero entries, or those states silently drop out.
- `Params::min_state_freq` (default set in `tools.cpp`, overridable with `--min-freq`) is the
  optimizer's lower bound on any frequency; `targetFunk` returns `1e+30` below it.

---

## 11. Checkpointing contract

IQ-TREE checkpoints continuously so runs can resume. Every model implements three methods and
they must be kept in sync with any new state, or a resumed run silently reverts.

```cpp
void MyModel::startCheckpoint()  { checkpoint->startStruct("MyModel"); }

void MyModel::saveCheckpoint() {
    startCheckpoint();
    if (!fixed_parameters) CKP_ARRAY_SAVE(getNumRateEntries(), rates);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();          // chain to the base class LAST on save
}

void MyModel::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();       // chain to the base class FIRST on restore
    startCheckpoint();
    if (!fixed_parameters) CKP_ARRAY_RESTORE(getNumRateEntries(), rates);
    endCheckpoint();
    decomposeRateMatrix();                  // rebuild derived state
    if (phylo_tree) phylo_tree->clearAllPartialLH();
}
```

`ModelUnrest` (`model/modelunrest.cpp:113–138`) is the reference. Macros live in
`utils/checkpoint.h:25–40`. **Note the ordering asymmetry — base-class call last on save, first
on restore — and the mandatory `decomposeRateMatrix()` + `clearAllPartialLH()` after restore.**

---

## 12. The partition / linked-model layer (QMaker and nQMaker)

`PartitionModel : public ModelFactory` (`model/partitionmodel.h:32`) estimates **one shared Q
across many partitions**, each with its own tree and branch lengths. This is how the published
`Q.*` matrices (QMaker) and `NQ.*` matrices (nQMaker) were estimated, and almost certainly how a
constrained nQ would be estimated too.

```
PartitionModel::optimizeLinkedModels()                partitionmodel.cpp:781
  for each distinct model name in linked_models:
      unfix that model's parameters in every partition
      optimizeLinkedModel()                           partitionmodel.cpp:692
        ├─ setVariables(v)          delegate to the representative model
        ├─ model->setBounds(...)
        ├─ minimizeMultiDimen(...)  → targetFunk
        │     PartitionModel::targetFunk              partitionmodel.cpp:299
        │       parallel-for over partitions sharing the model name:
        │         part_model->targetFunk(x)           each partition scores the SAME x
        │       return the sum
        ├─ getVariables(v)          push x into EVERY partition's model
        └─ decomposeRateMatrix() on each, clearAllPartialLH, recompute
      re-fix parameters; saveCheckpoint(); dump()
```

Key detail: **the same raw variable vector `x` is handed to every partition's model object.**
Any new parameterization must therefore be a pure function of `x` plus per-model constants —
no hidden per-partition state in `getVariables`.

Driven by:

- `--link-model` (`tools.cpp:4998`) and `--model-joint MODEL` / `--link-partition`
  (`tools.cpp:5003`), which set `Params::link_model` and `Params::model_joint`.
- `--model-joint NONREV` is the nQMaker invocation. `phyloanalysis.cpp:162–177` chooses the
  QMaker vs nQMaker citation by testing whether `model_joint` contains `"NONREV"`.
- `--init-model NAME|FILE|DIVMAT` (`tools.cpp:2896`) seeds the matrix;
  `DIVMAT` seeds from the empirical divergence matrix (`partitionmodel.cpp:61–120`).
- Results are written as NEXUS via `reportNexusFile` (`phyloanalysis.cpp:420–459`).

---

## 13. Params and the CLI contract

`Params` is a mutable global singleton, reached from anywhere as `Params::getInstance()`. There
is no dependency injection; model code reads `Params` directly (for example
`Params::getInstance().min_state_freq` inside `targetFunk`). Follow that convention rather than
threading new arguments through constructors.

**Adding an option requires four edits, all in `utils/`:**

1. `utils/tools.h` — declare the field in `struct Params`.
2. `utils/tools.cpp`, the `parseArg` initializer region (~7199 and ~7442) — set its default.
3. `utils/tools.cpp`, the `parseArg` argument loop — parse it, following the neighbouring
   `strcmp(argv[cnt], "--flag") == 0` pattern, including the missing-argument error.
4. `utils/tools.cpp`, `usage_iqtree()` (~5940–5990) — document it in the help text.

Note that `--model-joint` is **absent from `usage_iqtree()`**, so it does not appear in `-h`
despite being documented on the project website. If you add a target-π option, document both.

Model-name and frequency-type strings are a second, parallel CLI surface: `-m NAME{params}+F...`
is parsed by `ModelFactory`'s constructor and then by `createModel`
(`modelmixture.cpp:3122–3220`). `StateFreqType` (`tools.h:469–482`) is the enum that carries
frequency-handling intent through the whole model layer.

---

## 14. Mechanical checklist: adding a new substitution model

This is the plumbing, not the science.

1. **Create** `model/mymodel.h` and `model/mymodel.cpp`, subclassing `ModelMarkov`. Copy the
   shape of `model/modelunrest.{h,cpp}` — it is the shortest complete example.
2. **Register the source files** in `model/CMakeLists.txt`.
3. **Register the name.** Either add a `validModelName` clause reachable from
   `ModelMarkov::getModelByName` (`modelmarkov.cpp:1954`), or — for amino acids — add a branch to
   `ModelProtein::init` (`modelprotein.cpp:1106–1245`) next to the `GTR20` and `NONREV` branches.
   Confirm it is reachable from `createModel` (`modelmixture.cpp:3122`).
4. **Implement the Optimization interface**: `getNDim`, `setVariables`, `getVariables`,
   `setBounds`, and `setRates` if the parameterization is indirect. Honour the 1-indexing and the
   `changed` flag.
5. **Implement the three checkpoint methods** with the save/restore ordering from §11.
6. **Implement `writeInfo(ostream&)`** for the `.iqtree` report, and `getNameParams()` so the
   model round-trips through its own name string.
7. **Update the parameter counts** — `getNDim()` and `getNDimFreq()` — so AIC/BIC are right.
8. **Optionally register with ModelFinder** by adding the name to the tables in
   `main/phylotesting.cpp:160–230`.
9. **Add reporting** in `main/phyloanalysis.cpp` if a new output artifact is needed
   (`reportNexusFile` at 420 is the template), plus a citation block at 162–177.
10. **Smoke-test** against `example/aa_example.phy` before anything larger.

---

## 15. Conventions, idioms, and landmines

**Idioms to match:**

- `outError(msg)`, `outWarning(msg)`, `ASSERT(cond)` — never `throw`, `assert`, or `std::cerr`.
- `convertIntToString`, `convertDoubleToString`, `convert_double_with_distribution` from
  `utils/tools.h` for all string conversion.
- `verbose_mode >= VB_MED / VB_MAX / VB_DEBUG` gating for diagnostic output.
- `aligned_alloc<double>(n)` / `ensure_aligned_allocated(ptr, n)` / `aligned_free(ptr)` for any
  array the likelihood kernels will touch.
- Doxygen `/** @param ... @return ... */` comments on public methods.
- Dated change comments, for example `// BQM 2015-09-07: ...`. Follow the local style where you
  modify old code.
- 4-space indentation in newer code; older files use tabs. Match the file you are editing.

**Landmines:**

- **Off-by-one.** Parameter vectors are 1-indexed. Rate arrays are 0-indexed.
- **Forgetting `clearAllPartialLH()`** after changing Q or π — stale likelihoods, no warning.
- **Assuming π is an input for non-reversible models.** It is an output of
  `decomposeRateMatrixNonrev` on every single evaluation.
- **Assuming a non-reversible model has an unrooted tree.** `setReversible(false)` roots it.
- **Editing `cmaple/`** because a grep for `NONREV` landed there. Different engine.
- **Reading `model/modelgtr.cpp` or `model/modelnonrev.*`.** Not built.
- **`ModelMarkov::setRates()` is `ASSERT(0)`.** If your model uses an indirect parameterization
  you must override it, and you must actually call it from `getVariables`.
- **Complex eigen arrays alias the real ones.** `ceval`/`cevec`/`cinv_evec` point into the same
  allocation as `eigenvalues`/`eigenvectors`/`inv_eigenvectors`. Do not free or reallocate one
  without the other.
- **`nondiagonalizable`** must be set when the eigenvector matrix is singular; the transition-
  probability code then falls back to scaled squaring.
- **OpenMP.** `PartitionModel::targetFunk` runs partitions in parallel. Anything it touches must
  be thread-safe or per-partition.
- **`fixed_parameters`** short-circuits `getNDim()` to 0 and is toggled around linked-model
  optimization (`partitionmodel.cpp:781–807`). Respect it.

---

## 16. Insertion points for a constrained-π nQ

A factual list of the seams, for orientation only. No design is implied.

| Concern | Location |
|---|---|
| Where π is currently solved from Q | `modelmarkov.cpp:1267` inside `decomposeRateMatrixNonrev` |
| The Q→π solver itself | `computeStateFreqFromQMatrix`, `modelmarkov.cpp:2119` |
| Existing "pin π and skip the solve" switch | the `freq_type != FREQ_USER_DEFINED \|\| optimize_from_given_params` guards at `modelmarkov.cpp:1252` and `1266` |
| Existing one-shot "reweight Q by a π" | `adaptStateFrequency`, `modelmarkov.cpp:914` |
| Dimension count to adjust | `getNDim` (`964`) and `getNDimFreq` (`978`) |
| Free-parameter ↔ rate-matrix map | `setVariables` (`1018`), `getVariables` (`1038`), `setRates` (`1949`, base asserts) |
| **Working precedent for the whole feature (4-state DNA)** | `ModelLieMarkov::setBasis`, `modelliemarkov.cpp:1087`, with `piToTau` at `1032` and `BDF[]` at `247` |
| Reduced-parameterization template | `ModelLieMarkov::setRates`, `modelliemarkov.cpp:1194` |
| Live but unsound path to be left alone or warned about | `-m NONREV+F{...}`: `ModelProtein::init` sets `FREQ_USER_DEFINED`, which makes `decomposeRateMatrixNonrev` skip the stationarity solve while leaving Q unconstrained |
| Box constraints | `setBounds`, `modelmarkov.cpp:1139`; `ModelUnrest::setBounds`, `modelunrest.cpp:84` |
| Penalty / invalid-region signalling | `targetFunk`, `modelmarkov.cpp:1104–1110` |
| Where a target π could be supplied | `Params` + `parseArg` (`tools.cpp` ~2896 for the `--init-model` pattern); or a `+F{...}` string via `createModel` (`modelmixture.cpp:3205–3215`); or a NEXUS frequency vector via `ModelsBlock` |
| Amino-acid model entry point | `ModelProtein::init`, `modelprotein.cpp:1207` (`NONREV` branch) |
| Multi-partition estimation | `PartitionModel::optimizeLinkedModel`, `partitionmodel.cpp:692` |
| Emitting the estimated matrix | `reportNexusFile`, `phyloanalysis.cpp:420` |
| Consistency check that exists today | the π-mismatch warning in `readParameters`, `modelmarkov.cpp:1833–1842` |

---

## 17. Build and test loop on this machine (Windows / MSVC)

From a prior session on this machine, configuring from a clean clone needs explicit paths and a
CMake policy override:

```bash
cd build
"C:/Program Files/CMake/bin/cmake.exe" \
    -DEIGEN3_INCLUDE_DIR="C:/ProgramData/chocolatey/lib/eigen/include/eigen3" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
"C:/Program Files/CMake/bin/cmake.exe" --build . --config Release
```

Notes: CMake is not on `PATH`; Eigen3 and Boost are hard `find_package` requirements and are not
vendored; `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is required because the vendored
`zlib-1.2.7/CMakeLists.txt` declares a `cmake_minimum_required` that modern CMake refuses.

Fast functional checks:

```bash
iqtree3 -s example/aa_example.phy -m LG+G4          # reversible baseline
iqtree3 -s example/aa_example.phy -m NONREV         # non-reversible path, rooted tree
iqtree3 -s example/aa_example.phy -m NQ.pfam        # built-in nQ matrix
iqtree3 -s example/example.phy -m UNREST            # smallest non-reversible model, DNA
iqtree3 -s example/example.phy -p example/example.nex --model-joint NONREV   # nQMaker path
```

Inspect the resulting `.iqtree` report for the Q matrix and the state frequencies, and check the
π block against the stationary vector of the reported Q. Broader regression coverage lives in
`test_scripts/` (`test_iqtree.ps1`, `verify_results.ps1`, `test_configs.txt`).
