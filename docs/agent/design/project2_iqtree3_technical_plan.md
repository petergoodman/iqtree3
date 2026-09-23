# Project 2: fixed-equilibrium nonreversible estimation in IQ-TREE3

**Technical specification and independent synthesis — 22 September 2026**  
Prepared for Peter Goodman and the IQ-TREE collaboration. Mathematical notation uses row-vector distributions. This document defines a proposed implementation; it does not report a completed IQ-TREE extension.

## 0. Decision and evidence standard

Implement the complete fixed-equilibrium model as a new nonreversible model class inside IQ-TREE's existing direct-likelihood and linked-model optimization framework. Preserve the transition-matrix, pruning, branch-length, site-rate, and tree-search machinery wherever their existing semantics match the declared analysis. Change how the optimizer's 360 variables become a generator, protect the fixed target throughout the model lifecycle, and make the numerical derivative interface suitable for those variables.

**First implementation:** the jump-chain parameterization, with 20 rows of 18 free log ratios, a small stationary-vector solve, and derived holding rates. **Competing implementation:** gauge-fixed T3 balancing, with 189 pair-strength and 171 cycle coordinates. **Independent reference:** direct balanced-flux constrained optimization. Select the production default by matched numerical and runtime tests, not by the earlier reports' agreement or by an unmeasured claim that a small inner solve is negligible. Keep substitution-history EM as a subsequent option if direct likelihood optimization proves inadequate.

This prioritization differs from both original reports. At one fixed target, T3's favorable cross-composition interpretation is not required, and a reversible model can be seeded and checked perfectly well in jump-chain coordinates. Jump-chain construction eliminates an iterative inner optimizer while using linear algebra already familiar to IQ-TREE. This is an engineering reason for implementation order, **not a theorem that jump-chain fitting is faster or better conditioned**.

Three separate claims must never be conflated:

1. **Coverage:** the unrestricted coordinates represent every strictly positive normalized generator with the specified equilibrium. Proved below for both charts.
2. **Numerical optimization:** an actual finite-precision run reaches a well-verified local solution over its declared numerical domain. Requires implementation tests, gradients, convergence diagnostics, and restarts.
3. **Global statistical optimum:** the best generator, trees, roots, branch lengths, and rate parameters have been found. Neither the parameterization nor IQ-TREE's heuristic search guarantees this.

### Evidence labels

- **Derived:** proved here from the definitions, with hypotheses specified.
- **Source:** read directly in IQ-TREE3 at commit `63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de`, retrieved from `master` on 22 September 2026 through GitHub and a local checkout. Links below pin that commit.
- **Numerical:** newly written independent NumPy/SciPy checks in Appendix B, actually run for this synthesis.
- **Proposed:** implementation or scientific protocol not yet executed in IQ-TREE.
- **Unverified:** requires compilation, runtime tracing, real training data, or further scientific evidence.

The two plans, two audits, problem contract, and `AA_MODEL_INFERENCE(2).md` were leads for investigation, not authorities. Their reported historical numerical results are not adopted as results of this work. The attached nQMaker PDF was read directly, including visual inspection of printed page 1113.

## 1. Exact statistical target

### 1.1 Generator family

Let $n=20$, with amino-acid order `A R N D C Q E G H I L K M F P S T W Y V`. Specify one target vector $\pi^\star$ with strictly positive entries summing to one. For off-diagonal entries require $q_{ij}>0$, set $q_{ii}=-\sum_{j\ne i}q_{ij}$, and define

$$
\mathcal G_{\pi^\star}=\left\{Q:q_{ij}>0\ (i\ne j),\ Q\mathbf1=0,\ \pi^\star Q=0,\ -\sum_i\pi_i^\star q_{ii}=1\right\}.
$$

The full nonreversible family includes reversible points. Do not require that the estimated matrix exhibit nonzero circulation. The proposed family has 360 free matrix parameters; fixing a matrix for downstream evaluation leaves zero matrix parameters to optimize in that evaluation.

Strict positivity is the initial mathematical domain, as in the supplied contract. General valid generators can have zero rates. Section 5 handles that closure explicitly; no finite log-coordinate chart represents its boundary exactly. The class should reject a target containing zeros rather than silently change the target or pretend the same 20-state dimension theorem still applies. Tiny positive entries are mathematically allowed but require numerical support beyond a hardcoded frequency floor.

### 1.2 Likelihood and nuisance parameters

For training alignments $D_g$, pattern counts $c_{gp}$, rooted trees $T_g$, branch lengths $t_{gb}$, rate categories $r_{gk}$, and category weights $w_{gk}$, use the ordinary phylogenetic likelihood

$$
\ell(Q,\eta,T)=\sum_g\sum_p c_{gp}\log\left[\sum_k w_{gk}\sum_a\pi_a^\star U_{gpk,\mathrm{root}}(a)\right].
$$

Here $\eta$ collects continuous nuisance parameters and

$$
U_{gpk,v}(a)=\prod_{u\in\mathrm{children}(v)}\sum_b
\left[e^{r_{gk}t_{g,vu}Q}\right]_{ab}U_{gpk,u}(b).
$$

Tip vectors encode observed/ambiguous states. Additional partition multipliers multiply branch time. An invariant category has rate zero and uses the same fixed root distribution. Existing category normalizations and partition-rate gauges must be retained consistently.

The requested conditional matrix step is

$$
\sup_{Q\in\mathcal G_{\pi^\star}}\ell(Q,\eta^{\rm current},T^{\rm current}).
$$

The full training goal also updates the permitted nuisances and trees. Write a set-valued argmax only when attainment is established; a supremum is safer on an open domain. The optimizer returns an attained approximate solution with diagnostics, not a proof that this supremum has been reached.

The sum is over sites through pattern multiplicities. Equal weighting of genes instead would change the likelihood objective. Separate gene trees and a common edge-proportional tree are distinct statistical analyses; choose according to the actual experiment. Do not infer that `--model-joint` forces one topology: it links substitution parameters and can be used with different partition-tree arrangements.

Under stationary initialization the expected number of substitutions on a branch of duration $t$ is $\int_0^t\sum_i\pi_i^\star(-q_{ii})\,ds=t$. With a scalar multiplier $r$ it is $rt$. Thus unit equilibrium mean-rate normalization preserves the familiar branch-length units. If one rescales $Q$ by $c$, branch lengths must be divided by $c$ to preserve transition probabilities; this scale equivalence is subject to any branch bounds or calibration restrictions.

### 1.3 Equilibrium versus root distribution

For this stationary model, set the root distribution to $\pi^\star$ as well as imposing $\pi^\star Q=0$. These are two logically different conditions. A nonstationary root distribution with a target-stationary generator is a valid alternative model, but it is outside this baseline and requires its own parameters and interpretation.

An irreducible shared $Q$ cannot have distinct stationary targets for different genes. Sharing coordinates but rebuilding at gene-specific targets gives different matrices and reopens Project 1's cross-composition modeling choice.

### 1.4 What is scientifically being estimated

The estimand is the best cleaned-data fit within the fixed-equilibrium working model. It is not generally the old cleaned-data estimate with its frequencies replaced. In particular, fitting reversible exchangeabilities at the target is not generally equivalent to fitting them at another composition and swapping frequencies afterward.

The target's provenance must be recorded: original training data or an external source; residue, taxon and gene weighting; ambiguity and gap handling; any pseudocount; and the exact normalized numeric vector. If estimating the target from the same training genes, the optimization is conditional on that plug-in vector. It is not automatically a 379-parameter joint MLE, and adding 19 to a BIC count does not establish valid BIC for that two-stage procedure.

## 2. Mathematical reference: balanced stationary flux

### 2.1 Linear constraints and dimension — derived

For each directed edge $e=(i,j)$, $i\ne j$, define $f_{ij}=\pi_i^\star q_{ij}$. Let $B$ have column $e_j-e_i$ for that edge; let $D$ be $B$ with its last row omitted, and set

$$
C=\begin{bmatrix}D\\\mathbf1^T\end{bmatrix},\qquad b=(0,\ldots,0,1)^T.
$$

Then

$$
(Bf)_j=\sum_{i\ne j}f_{ij}-\sum_{k\ne j}f_{jk}=(\pi^\star Q)_j,
\qquad \mathbf1^Tf=-\sum_i\pi_i^\star q_{ii}.
$$

Thus the model is exactly $f>0$, $Cf=b$, followed by $q_{ij}=f_{ij}/\pi_i^\star$ and row-sum diagonals.

**Rank proof.** If $y^TB=0$, then $y_j-y_i=0$ on every edge, so $y$ is constant. Hence $\operatorname{rank}B=n-1$. The normalization row is independent: expressing its all-one entries as $y_j-y_i$ would require both $y_j-y_i=1$ and $y_i-y_j=1$. Consequently $\operatorname{rank}C=n$. Uniform directed flux $1/[n(n-1)]$ is strictly feasible, proving the affine constraints have a nonempty positive relative interior. Therefore

$$
\dim\mathcal G_{\pi^\star}=n(n-1)-n=n(n-2)=360.
$$

The drift vector lies in $\mathbf1^\perp$, not generally in $(\pi^\star)^\perp$. Its entries sum to zero because $Q\mathbf1=0$.

The flux and generator maps are inverse linear maps at fixed positive $\pi^\star$. This proves completeness without relying on either agent's preferred parameterization. $C$ is independent of the target, so its null-space basis can be reused across targets.

### 2.2 Reversible and circulation dimensions — derived

Write $F=S+J$ off diagonal, with $S=(F+F^T)/2$ symmetric and $J=(F-F^T)/2$ antisymmetric. Symmetric flux is automatically balanced. Balance imposes zero divergence on $J$; its incidence matrix on the complete undirected graph has rank $n-1$ by the same constant-vector argument.

There are $\binom n2=190$ symmetric strengths and one normalization $2\sum_{i<j}S_{ij}=1$, leaving 189. The antisymmetric divergence-free dimension is $190-19=171$. Positivity couples them through $S_{ij}>|J_{ij}|$; these are not arbitrary independent unrestricted flux blocks. Setting $J=0$ gives the reversible subfamily. Thus $360=189+171$.

For the unrestricted positive generator, 380 off-diagonal rates determine the diagonal; a single positive scale normalization removes one dimension. Its unique equilibrium varies with the matrix, yielding 379 total dimensions. The fixed-target restriction removes 19. Parameter dimension is not a proof of statistical identifiability from a given tree and alignment.

### 2.3 Simple reversible seeds — derived

For each symmetric positive empirical exchangeability matrix $R$, form

$$
\widetilde q_{ij}=R_{ij}\pi_j^\star,\quad
\widetilde q_{ii}=-\sum_{j\ne i}\widetilde q_{ij},\quad
Q=\widetilde Q/\left(-\pi^\star\operatorname{diag}\widetilde Q\right).
$$

Detailed balance follows from $\pi_i^\star R_{ij}\pi_j^\star=\pi_j^\star R_{ji}\pi_i^\star$. It implies stationarity, and scalar normalization preserves it. This makes target versions of LG, JTT, and WAG straightforward. If a custom empirical matrix contains zeros, use a disclosed positive exchangeability mixture for an interior seed, not independent clipping of a nonreversible matrix.

For an arbitrary asymmetric factor $A$, using $q_{ij}=A_{ij}\pi_j$ instead gives $(\pi Q)_j=\pi_j[(A^T-A)\pi]_j$. This vanishes for all positive normalized $\pi$ only if $A^T=A$: linearity and an open simplex subset force the linear map to vanish everywhere. An asymmetric frequency substitution still makes a valid generator when $A\ge0$; the generally lost property is the requested stationarity, not row-sum generator validity.

### 2.4 Nonreversible seeds — derived

Given any strictly positive fitted generator $Q_0$, solve its equilibrium $\mu$, normalize its mean rate, and calculate $F_0=\operatorname{diag}(\mu)Q_0$ off diagonal. Then $q_{ij}^{\rm seed}=(F_0)_{ij}/\pi_i^\star$ is positive, balanced at the target, and normalized. This is a flux-preserving seed transfer. An alternative is T3 transport, decomposing at $\mu$ and rebuilding at $\pi^\star$.

Neither rule need define the final estimator: optimize all 360 coordinates afterward. The reversible fitted target matrix, transferred nonreversible matrices, and small feasible perturbations give complementary starting points. A feasible flux mixture $(1-\epsilon)F_0+\epsilon F_1$ remains normalized and balanced; clipping entries individually generally does not.

## 3. Interior implementation A: jump-chain coordinates

Use $K$ for the jump matrix to distinguish it from the finite-time transition matrix $P(t)=e^{tQ}$.

### 3.1 Build — derived

For each state $i$, choose a fixed reference destination $r(i)\ne i$. There are 19 possible destinations. Store 18 logits $z_{ij}$ against the reference, set the reference logit to zero, and apply a numerically stable row softmax:

$$
K_{ii}=0,\qquad K_{ij}=\frac{e^{z_{ij}}}{\sum_{k\ne i}e^{z_{ik}}},\qquad z_{i,r(i)}=0.
$$

Compute the unique stationary row vector $\nu K=\nu$, $\nu\mathbf1=1$. Then set

$$
f_{ij}=\nu_iK_{ij},\qquad q_{ij}=\frac{\nu_i}{\pi_i^\star}K_{ij},\qquad q_{ii}=-\frac{\nu_i}{\pi_i^\star}.
$$

Interpretation: $K$ chooses the next amino acid conditional on a substitution occurring. The exit rate $\lambda_i=\nu_i/\pi_i^\star$ determines the waiting time. $\nu$ is the distribution of states at jump events; it is **not** the equilibrium distribution over elapsed time. These holding rates make the latter equal to $\pi^\star$.

### 3.2 Existence and exact validity — derived

All off-diagonal entries of $K$ are positive, so its directed graph is irreducible. A stationary vector exists: Cesàro averages of any initial probability vector propagated by $K$ remain in the compact simplex, and every limit point is invariant because the difference after one application of $K$ is a bounded telescoping term divided by the averaging length. An invariant nonnegative vector is strictly positive by irreducibility. The right nullspace of $I-K$ consists of constants: a maximum component of a harmonic vector must equal all neighbors with positive weights, and irreducibility propagates equality. Thus rank is $n-1$, the left nullspace is one dimensional, and normalization gives unique $\nu$. Aperiodicity is not required for this stationary solve.

Incoming flux at $j$ is $\sum_i\nu_iK_{ij}=\nu_j$; outgoing flux is $\nu_j$. Total directed flux is $\sum_i\nu_i=1$. Therefore every finite logit vector gives a positive target-stationary unit-mean-rate generator.

### 3.3 Completeness, inverse, and smoothness — derived

Given any $Q\in\mathcal G_{\pi^\star}$, set $\lambda_i=-q_{ii}>0$ and $K_{ij}=q_{ij}/\lambda_i$ for $i\ne j$, $K_{ii}=0$. Then $\nu_i=\pi_i^\star\lambda_i$ is positive, sums to one by mean-rate normalization, and satisfies $\nu K=\nu$ by stationarity. The build therefore returns $Q$. The inverse logits are $z_{ij}=\log(K_{ij}/K_{i,r(i)})$. Conversely the jump matrix of the built $Q$ is exactly $K$. The map is bijective and uses $20(19-1)=360$ coordinates.

The normalized stationary linear system is nonsingular on its rank-completed space. Its solution depends smoothly on the entries while rank is maintained. Combined with softmax and the explicit inverse, this gives a smooth global chart of the positive family. It removes representation redundancy, not possible statistical nonidentifiability.

### 3.4 Derivatives — derived

For one row,

$$
\delta K_{ij}=K_{ij}\left(\delta z_{ij}-\sum_{k\ne i}K_{ik}\delta z_{ik}\right).
$$

Differentiating stationarity and normalization gives $\delta\nu(I-K)=\nu\delta K$, $\delta\nu\mathbf1=0$. Hence

$$
\delta\nu=\nu\delta K(I-K+\mathbf1\nu)^{-1}.
$$

The inverse exists: a left null vector of this augmented matrix has zero sum after multiplying by $\mathbf1$, is then a stationary left vector of $K$, and must be zero. Implement linear solves rather than explicitly forming an inverse. Finally

$$
\delta f_{ij}=\delta\nu_iK_{ij}+\nu_i\delta K_{ij},\quad
\delta q_{ij}=\delta f_{ij}/\pi_i^\star,
\quad \delta q_{ii}=-\sum_{j\ne i}\delta q_{ij}.
$$

A full analytic likelihood gradient still requires differentiating the tree likelihood. A cheap derivative of the builder alone does not supply that gradient.

### 3.5 Numerical contract and tradeoffs

Use stable softmax with reference destinations and state order recorded in the model definition/checkpoint. Softmax underflow to an exact zero is a numerical failure for the interior model, not an accepted sparse proposal. Check stationary-solve residuals, positivity, normalization, and sensitivity near nearly reducible chains. The existing QR stationary routine can supply a prototype solve on $K-I$, but its assertion on the sum alone is insufficient validation.

The reversible subset is defined by $\nu_iK_{ij}=\nu_jK_{ji}$, not by a coordinate block set to zero. This makes a reversible-only fit better handled by existing reversible machinery. It does not prevent exact import of its $Q$, exact likelihood comparison, or complete coverage of the reversible points within the general model.

Holding these jump coordinates fixed while changing target $\pi$ preserves flux, not generally reversible exchangeabilities. That is acceptable for this fixed-target optimization; it would be a substantive choice if reused in Project 1.

## 4. Interior implementation B: T3 balancing

### 4.1 Coordinates and build — derived

Choose symmetric log strengths $a_{ij}=a_{ji}$ and antisymmetric tilts $h_{ij}=-h_{ji}$. Pin $a_{n-1,n}=0$ to remove one strength scale. Pin $h_{i,n}=0$ for all $i<n$ to remove the 19 gradient directions. These leave 189 strength coordinates and 171 independent tilts.

Define $N_{ij}=\pi_i^\star\pi_j^\star e^{a_{ij}+h_{ij}}$ off diagonal. With $g_n=0$, minimize

$$
\Phi(g)=\sum_{i\ne j}N_{ij}e^{g_j-g_i}.
$$

At its minimizer set $x_{ij}=N_{ij}e^{g_j-g_i}$, $f=x/(\mathbf1^Tx)$, and reconstruct $Q$ from flux. The 19 potentials are derived solver outputs, not 19 additional model parameters.

### 4.2 Inner solution exists, is unique, and enforces balance — derived

The reduced gradient is $Dx$ and Hessian is

$$
H=D\operatorname{diag}(x)D^T.
$$

For a reduced vector $v$, extend it by $v_n=0$. Then $v^THv=\sum_{i\ne j}x_{ij}(v_j-v_i)^2>0$ unless $v=0$. Thus $\Phi$ is strictly convex on the slice. It is also coercive: unbounded $g$ on this slice has unbounded range, and the positive edge coefficient from a minimum to a maximum multiplies an exponential of that range. With finitely many positive coefficients their minimum is positive, giving a uniform lower bound that diverges. A unique finite minimizer exists. Its zero reduced gradient gives all balance equations because the omitted row is dependent. Normalization preserves balance.

The invertible reduced Hessian and implicit differentiation give smooth dependence on the coordinates. These facts justify the mathematical map, not every floating-point Newton implementation.

### 4.3 Inverse and uniqueness — derived

For a feasible $Q$ define

$$
a^0_{ij}=\tfrac12\log\frac{q_{ij}q_{ji}}{\pi_i^\star\pi_j^\star},\qquad
h^0_{ij}=\tfrac12\log\frac{\pi_i^\star q_{ij}}{\pi_j^\star q_{ji}}.
$$

Before gauge changes, $N=F$ is balanced and $g=0$ is the solution. Subtract $a^0_{n-1,n}$ from all strengths. Put $u_i=h^0_{i,n}$ and replace $h^0_{ij}$ by $h^0_{ij}+u_j-u_i$. This sets the star tilts to zero. The strength shift rescales $x$ uniformly and cancels in normalization; the tilt gradient is absorbed by an opposite change in $g$. This proves surjectivity.

If two coordinate sets give the same normalized flux, taking pairwise log products determines their symmetric strengths up to a common additive constant. Taking pairwise log ratios determines their tilts up to a vertex gradient. The two specified gauges remove these freedoms uniquely. Thus the chart is injective, and its explicit logarithmic inverse is smooth.

At $h=0$, $N$ is symmetric, $g=0$ balances it, and $q_{ij}\propto e^{a_{ij}}\pi_j^\star$. Conversely a reversible output has a gradient tilt, which is zero in the star gauge. Hence $h=0$ is exactly the reversible slice.

### 4.4 Builder derivative — derived

Let $v_{ij}=\delta a_{ij}+\delta h_{ij}$. Differentiating $Dx=0$ yields

$$
H\delta g=-D\operatorname{diag}(x)v,\qquad
\delta x=\operatorname{diag}(x)(v+D^T\delta g).
$$

For $s=\mathbf1^Tx$,

$$
\delta f=\delta x/s-f(\mathbf1^T\delta x)/s.
$$

Divide each off-diagonal row by its fixed target entry and complete the diagonal. A common scaling of $x$ cancels from the linear equation, so normalized flux may replace $x$ on both sides. Ignoring the derivative of the balancing potentials is incorrect.

### 4.5 Numerical design and the invalid bound of 12

Use damped Newton on the reduced system, stable exponential scaling, a domain-aware line search, a residual test, and a documented iteration limit. Near rounding-level decreases, use controlled residual/decrease acceptance rather than blindly disabling the line search. Warm starts are optional and per object; the converged result must be insensitive to evaluation order at the accuracy required by the outer derivative.

A source report proposed $|h|\le12$ as covering IQ-TREE's raw-rate box. This is false after gauge fixing. In a star gauge at state 1,

$$
h'_{23}=\tfrac12\log\frac{q_{12}q_{23}q_{31}}{q_{21}q_{32}q_{13}}.
$$

Set those three forward rates to 100, reverse rates to $10^{-4}$, and every other off-diagonal in a 20-state generator to 1. Complete the diagonal. Incoming and outgoing rates balance, so equilibrium is uniform. All raw rates lie in the claimed box, but $h'_{23}=\tfrac32\log(10^6)=20.7232658\ldots$. Scalar normalization leaves it unchanged. This independently derived and executed counterexample invalidates the containment claim.

Bounds on the un-gauged pairwise tilt are not bounds on triangle/cycle coordinates. Strengths after the scale gauge are ratios; changing to logarithms without changing the physical bounds does not enlarge their feasible range. No arbitrary fixed box should be described as the entire positive family.

## 5. Boundaries, complete coverage, and alternative direct coordinates

### 5.1 Closure and existence — derived

The closed flux polytope $\overline{\mathcal F}=\{f\ge0:Cf=b\}$ is compact because each flux lies in $[0,1]$. Consequently each off-diagonal rate lies in $[0,1/\pi_i^\star]$. For fixed finite nuisances, ordinary pattern likelihoods are continuous in $Q$ by continuity of the matrix exponential and finite sums/products. The full likelihood attains a maximum on the closure; equivalently the extended log-likelihood is upper semicontinuous. This maximum can be on the boundary. The positive family need not attain it, and boundary chains need not have a unique equilibrium, although the specified target still is stationary.

The same existence reasoning extends to a compact continuous nuisance domain and a finite set of trees. It does not establish existence for noncompact nuisance domains. Nor does it prove uniqueness or identifiability.

Every closure point is approximable by positive points: mix its flux with uniform positive flux using weight $\epsilon>0$. Hence the interior likelihood supremum equals the closure maximum for fixed nuisances when finite likelihood is available. This establishes approximation of the model, not convergence of a local optimizer to its global supremum.

### 5.2 Operational policy

The first release should describe itself as optimizing the full positive family with monitored numerical limits. Start with finite log-coordinate ranges chosen to contain the actual seeds with margin, expand if bounds are active or the solution changes materially, and report the final ranges and sensitivity. Do not invent universal safe limits. An increasing sequence of boxes exhausting all real coordinates recovers the full model's supremum **if their global suprema were computed**; heuristic local fitting does not inherit that guarantee.

When solutions approach zero flux, compare against a direct flux solver, optionally with descending declared floors. A uniform flux floor $\epsilon<1/380$ preserves a nonempty full-dimensional interior but restricts the model; it means $q_{ij}\ge\epsilon/\pi_i^\star$. A raw-rate floor means a different set. In particular, retain the distinction between a 360-dimensional restricted box and coverage of the complete 360-dimensional family.

If exact zero-rate maximum-likelihood solutions become a production requirement, add a boundary-capable flux optimization/polishing backend. The log-coordinate implementations cannot satisfy that requirement exactly. A boundary fit still needs convergence and active-set diagnostics.

### 5.3 Null-space flux solver — derived

Let $Z$ have orthonormal columns spanning $\ker C$ and choose positive feasible $f_0$. Write $f=f_0+Zy$. Every feasible flux has the unique inverse $y=Z^T(f-f_0)$, so this representation is complete on its coupled-inequality domain. For a tangent direction $d=Zv$ and lower bounds $l_e$,

$$
\alpha_{\max}=\min_{e:d_e<0}\frac{f_e-l_e}{-d_e}.
$$

A fraction-to-boundary feasible line search respects the inequalities. Coordinatewise box clamping does not. SQP, an appropriate interior-point solver, or a carefully implemented reduced-space method can optimize the same nonconvex likelihood. The domain's not being all of $\mathbb R^{360}$ does **not** make this an incomplete representation.

A Hamiltonian-cycle elimination is another linear coordinate choice: designate cycle fluxes $y_j=f_{j,j+1}$, let $b_j$ be free-edge outgoing minus incoming flux, solve $y_{j-1}-y_j=b_j$, then use total flux one to fix the additive cycle constant. Positivity of all reconstructed edges remains a coupled constraint. This can reduce coding for a reference solver but provides no reason to pass infeasible finite-difference probes to the likelihood.

### 5.4 Same model, different trajectories — derived

If $B_1,B_2$ are onto the same family, the objective values $\ell(B_1(\theta_1))$ and $\ell(B_2(\theta_2))$ range over the same set; their suprema and matrix-valued maximizer sets agree. Homeomorphic charts also preserve local maxima. They need not produce the same numerical optimization trajectory, basin, step scale, or stopping behavior.

No Jacobian determinant belongs in a frequentist likelihood merely because parameters were renamed. Bounds must describe the same physical domain to preserve comparisons; a penalty defined on $Q$ can be chart invariant, whereas unrelated coordinate penalties generally are not. Two estimates being interior to different boxes does not prove those boxes are nested.

## 6. Direct likelihood optimization: numerical requirements

### 6.1 The likelihood is not generally concave — derived counterexample

Use uniform equilibrium on three states and

$$
Q(x)=\begin{pmatrix}-1&(1+x)/2&(1-x)/2\\(1-x)/2&-1&(1+x)/2\\(1+x)/2&(1-x)/2&-1\end{pmatrix},\quad |x|<1.
$$

These form an affine positive normalized stationary family. On a two-tip rooted tree with lengths $t_1,t_2$, pattern $(0,0)$ has probability

$$
L_{00}(x)=\frac19\left[1+2e^{-3(t_1+t_2)/2}\cos\left(\frac{\sqrt3(t_2-t_1)x}{2}\right)\right].
$$

To obtain this formula, diagonalize the circulant: its eigenvalues are $0,-3/2\pm i\sqrt3 x/2$. The joint matrix is $e^{t_1Q^T}\operatorname{diag}(\pi)e^{t_2Q}$, giving the cosine term. Set $(t_1,t_2)=(0.1,5.1)$ and $x=2\pi/(5\sqrt3)$. With $A=2e^{-7.8}$ and $k=5\sqrt3/2$, the log-likelihood second derivative is $Ak^2/(1-A)>0$. Thus convex feasibility does not imply a concave observed-data likelihood.

For a 20-state counterexample, combine this block with a symmetric 17-state unit-exit-rate block, choose uniform equilibrium on 20 states, and connect blocks by small positive symmetric rates. Mean-rate normalization is independent of $x$. At zero coupling the pattern probability is scaled by $3/20$ and has the same log curvature. Continuity preserves positive curvature for sufficiently small positive coupling, giving a dense 20-state example.

### 6.2 Reuse BFGS, fix the derivative policy

The inspected single-model and linked-model paths call `minimizeMultiDimen`, whose `dfpmin` performs BFGS updates. The current numerical derivative uses $h_i=10^{-4}|x_i|$, replacing the step by $10^{-4}$ only when it is exactly zero. It can therefore use arbitrarily tiny probes near zero. T3 tilts and jump logits naturally cross zero. Shifting them by a constant moves the problem to another part of the domain and changes step sizes; positive ratios are another representation but do not guarantee gradient accuracy.

**Proposed minimal targeted change:** keep the BFGS framework but give the new model a parameter-scale-aware difference policy, available to both the single-model optimizer and the `PartitionModel` optimizer object. Use $h_i=\eta\max(s_i,|x_i|)$ with positive characteristic scale, bounds-aware one-sided or central differences, and a step-size study. Do not declare one universal $\eta$ adequate. Central differences are valuable for validation; their doubled likelihood cost should be measured before adopting them for every production gradient.

A model-subclass override alone is insufficient for joint training because the optimizer's `this` is a `PartitionModel`. Add a shared derivative-policy helper or explicit delegation in that path. Evaluate the summed linked objective at each probe, never just the representative gene's objective. After probing, restore the authoritative vector, generator, decompositions and caches; the current derivative routine restores the array entries but its calls mutate model state.

Do not change legacy-model numerical behavior globally unless maintainers deliberately approve that broader change. An opt-in helper for the new class reduces regression scope. Do not infer classical BFGS convergence guarantees from the method name when using numerical gradients and clamped trials.

### 6.3 Analytic derivatives of the likelihood — derived

For a perturbation $E=\delta Q$,

$$
\delta e^{tQ}=\int_0^t e^{(t-s)Q}Ee^{sQ}\,ds.
$$

This follows by differentiating the matrix evolution equation and variation of constants. It is not generally $tEe^{tQ}$. Differentiate pruning with this Fréchet derivative, or use an adjoint implementation. Existing branch derivatives do not suffice: branch time commutes with $Q$, but arbitrary matrix perturbations do not.

If $G$ is the full-matrix gradient of $\ell$, diagonal completion gives $\partial\ell/\partial f_{ij}=(G_{ij}-G_{ii})/\pi_i^\star$. Chain this through the chosen builder. Verify full likelihood directional derivatives, not just builder derivatives, over several probe sizes and both exponential paths.

Analytic gradients are a planned performance upgrade if finite differences are too expensive or fail accuracy gates; they are not a prerequisite to deriving the model. An accurate independent derivative oracle is required before treating small reported coordinate gradients as convincing convergence evidence.

### 6.4 Convergence diagnostics

For $\varphi(f)=-\ell(Q(f))$, an interior first-order condition is $Z^T\nabla_f\varphi=0$. With floors $f\ge l$, necessary KKT equations are

$$
Cf=b,\quad \mu\ge0,\quad f\ge l,\quad
\nabla_f\varphi+C^T\lambda-\mu=0,\quad \mu_e(f_e-l_e)=0.
$$

They follow from Lagrange multipliers for affine equalities and active inequalities; suitable local qualification is understood. They are not global certificates for this nonconvex objective. Diagnose feasibility, projected/physical gradient or KKT residual, likelihood change, step size, line-search status, bound activity, and restart spread. A chart can become ill conditioned and make its coordinate gradient small even when physical flux directions remain useful.

Keep the best feasible incumbent and accept updates only after evaluating their actual likelihood under the unchanged nuisance specification. An interrupted or failed inner solver must not overwrite a better model. Stop outer iterations only after inner accuracy and nuisance-update checks have passed. Pearson correlation between rate matrices is a descriptive continuity statistic, not an optimality condition.

## 7. What the paper says and what current IQ-TREE actually exposes

### 7.1 Paper-level workflow

The attached Dang et al. nQMaker paper, printed pp. 1112–1113, initializes empirical candidates, chooses a candidate/rate model and tree for each alignment (or a clade-linked tree), then alternates fitting a shared matrix at fixed nuisance values with fitting branch lengths at fixed matrix. It adds learned matrices to the candidate pool and repeats the outer procedure.

Your description of the matrix substep 3a is correct at this level. But the text's step 4 literally assigns `QBEST := QNEW` before testing their correlation. Taken literally this compares a matrix to itself. This was confirmed on the PDF image, not assumed to be extraction noise. The implementable interpretation is to save the previous matrix and compare it with the new one before overwriting. This is a logical correction to the printed order, not a claim about the authors' hidden workflow scripts or an official erratum.

### 7.2 Source-level facts

The following are source findings at the pinned commit. They describe the inspected paths, not an exhaustive audit of every IQ-TREE model.

| Source / symbols | Verified behavior | Consequence |
|---|---|---|
| `model/modelprotein.cpp`, `ModelProtein::init` | `NONREV` initializes from LG or a supplied model, converts reversible input when needed, and uses `getNumRateEntries()-1` parameters | 379 optimized directed-rate entries under the existing convention; new class needs 360 coordinates |
| `model/modelmarkov.cpp`, `setReversible` | Reversible-to-nonreversible conversion multiplies exchangeabilities by the **current** frequencies | Apply target frequencies before constructing a reversible seed; merely changing `state_freq` after conversion is insufficient |
| `decomposeRateMatrixNonrev` | Rebuilds the full matrix from `rates[]`; solves equilibrium unless user-defined frequencies and no `optimize_from_given_params`; then normalizes | Writing only `rate_matrix` can be overwritten. Fixing `state_freq` alone leaves rates unconstrained |
| `computeStateFreqFromQMatrix` | Column-pivoted QR solves the augmented stationary system; row-major input is mapped with transpose effect | Reusable for $K-I$ after explicit layout/residual tests; do not overwrite target with the jump stationary vector |
| `getNDim`, `getNDimFreq` | Fixed-parameter gate returns zero; fixed supplied frequencies add no frequency dimension | Preserve fixed versus fitted semantics and shared parameter counts |
| `setVariables`, `getVariables`, `targetFunk` | Packs/unpacks rates, decomposes on change, clears partial likelihoods, returns negative tree log-likelihood | Override packing/unpacking for chart coordinates; preserve cache invalidation |
| `targetFunk`, `setBounds` | Small positive frequencies below `min_state_freq` yield barrier value; base rate bounds are $[10^{-4},100]$ | Validate target support explicitly; override coordinate bounds |
| `model/modelfactory.cpp`, constructor | `model_joint` replaces the substitution model string; its frequency suffix is extracted; per-partition rate strings are retained | Fixed target can be routed through the joint construction path; runtime test still required |
| `model/partitionmodel.cpp`, constructor | Links objects by model name; pools frequency counts only for estimated/empirical frequency modes | User-defined target avoids that pooling path, but equal names alone do not prove equal target/gauge metadata |
| `PartitionModel::targetFunk` | Calls each matching partition model's objective in parallel, sums stored values in a fixed serial order | Preserve summed objective and deterministic reduction; builder can run once per partition per probe |
| `optimizeLinkedModel`, `optimizeLinkedModels` | Separate BFGS entry point; unfixes linked models for update, refixes afterward | Gradient changes must reach this path; no accidental per-gene matrix refitting |
| `PartitionModel::getNParameters` | Adds linked model dimensions once after per-partition fixed-model counts | New shared matrix contributes 360 once |
| `model/partitionmodelplen.cpp`, `optimizeParameters` | Updates partition model/rate parameters, linked matrix, gene multipliers, then branch lengths | Current linked code is not the paper's strict fixed-$V$ inner loop |
| `model/partitionmodel.cpp`, `optimizeParameters` | Alternates per-partition nuisance optimization and linked-model updates | The separate-tree/unlinked branch path also needs integration tests |
| `main/phylotesting.cpp`, `testModel` vicinity | Nonempty `model_joint` suppresses the empty-model trigger for ModelFinder | One joint-fit command does not automatically reproduce outer candidate reselection |
| `getModelSubst`, `mixRevNonrev` | Recognized nonreversible candidates require opt-in; recognized mixed reversible/nonreversible sets are rejected | The paper's mixed expanding candidate pool cannot simply be assumed to work in one ModelFinder call |
| `CandidateModelSet::generate` | Builds matrix × frequency × rate candidates | Prevent unwanted `+F`/`+FO` variants from changing the target; freeze learned candidate matrices |
| `utils/optimization.cpp` | Relative forward differences, BFGS path, clamping of trial bounds | Use the targeted derivative policy in Section 6 |
| `model/modelmarkov.cpp`, `computeTransMatrixNonrev` | Eigen-based and scaling-and-squaring paths with noninvertible-eigenvector/row-sum fallback logic | Reuse initially; row sums alone are not a comprehensive exponential accuracy test |
| `tree/phylotree.cpp`, `getNBranchParameters` | Removes dummy-root branch; removes an extra root split degree for reversible trees | Root edge versus root split matters for matched testing |
| `model/modelfactory.cpp`, `optimizeParameters` | Root search is conditional on root flags and rootedness after its continuous loop | Topology-fixed input alone must not be asserted to fix all root behavior |

Pinned source links: [ModelProtein][S1], [ModelMarkov][S2], [ModelMarkov interfaces][S3], [ModelFactory][S4], [PartitionModel][S5], [PartitionModelPlen][S6], [optimizer][S7], [ModelFinder][S8], [tree/root logic][S9], [ModelUnrest][S10], [nonreversible kernels][S11], [argument parsing/defaults][S12].

Literal `+F{...}` parsing is present in ModelFactory and ModelProtein. The source trace also shows the joint string is handed to per-partition factories and user-defined frequencies bypass pooled-frequency initialization. This is stronger than the earlier documents' untraced speculation, but **it is not an executed end-to-end validation**. The new class still needs a two-partition test confirming identical targets and matrices before, during, and after optimization and restart.

### 7.3 What can remain, what must change

| Component | Proposed treatment |
|---|---|
| Data/pattern encoding and rooted pruning | Retain |
| $P(t)=e^{tQ}$ and existing branch-length derivatives | Retain, with numerical stress checks |
| Branch, rate, topology and root algorithms | Retain under explicit nuisance policy; their numerical results will change when $Q$ changes |
| Matrix update | Replace raw directed-rate variables with complete fixed-target coordinates |
| Initial empirical candidates | Rebuild at target if choosing a fully target-consistent candidate protocol |
| Candidate frequency variants | Restrict to the fixed target; no accidental fitted/empirical composition |
| Learned candidates in later rounds | Export/import as fixed valid matrices, preserve target, rooted semantics and zero matrix fit dimension |
| Mixed candidate orchestration | Explicit adapter or separate supported runs; do not merely remove the hard-error guard |
| Derivative and bounds policy | Add targeted support for the new coordinates in both optimizer entry points |
| Reporting, checkpoints, linking | Extend to store and enforce target/chart metadata |
| Stopping | Keep familiar loop structure, add feasibility/optimization checks and best-incumbent protection |

The initial candidate list is not mathematically required to satisfy the target if it is used **only** to obtain provisional guide trees. Unconstrained guide trees do not invalidate a subsequent constrained likelihood fit. Nevertheless, using target-rebuilt reversible candidates gives a coherent all-stages protocol and is simple here. If the target is imposed only after guide-tree creation, label that policy and assess initialization sensitivity rather than calling it mathematically invalid.

### 7.4 Recommended workflow levels

**Level 1 — fixed nuisance validation:** fixed rooted trees, branch lengths and rate parameters. Expose the conditional linked matrix objective and verify only $Q$ changes. This is the clean test of paper step 3a.

**Level 2 — native IQ-TREE continuous alternation:** retain the existing linked-model scheduling, including permitted continuous site-rate updates and branch/gene-rate updates. This minimizes deviation from current architecture. Report that it is broader than the paper's strict fixed-$V$ subloop. Rate-family selection and re-estimation of continuous parameters within a selected family are distinct operations.

**Level 3 — paper-style outer reselection:** explicitly orchestrate candidate-based topology/rate selection and repeat the shared matrix fit. Keep the previous shared-model solution as an incumbent; per-gene candidate scores optimize a surrogate set of different matrices and do not prove improvement in the single-shared-$Q$ objective. Evaluate and accept based on the actual common-model objective under a comparable nuisance specification.

For minimal production scope, deliver Levels 1 and 2 first plus fixed-matrix export/reuse. Add Level 3 if reproducing the paper's full outer candidate loop is a requirement of the collaborative analysis. A thin driver using supported IQ-TREE calls is a sensible first orchestration implementation; native integration can follow maintainers' preferences. Do not market `--model-joint` alone as a verified implementation of every box in the paper's flowchart.

## 8. Concrete IQ-TREE integration specification

### 8.1 Model identity and scope

Use a dedicated class, tentatively `ModelNonRevFixedEq`, derived from the nonreversible `ModelMarkov` infrastructure. This is a proposed name, **not existing syntax**. `ModelUnrest` is useful as a lifecycle example, not code to copy without review: its last-rate pin and no-op frequency setter do not express this model.

Initially support 20-state homogeneous protein generators with scalar site-rate mixtures, one target per linked matrix, and the existing single-tree/partition paths under explicit tests. Do not implicitly claim support for profile mixtures, site-specific frequencies, branch-specific targets, PoMo, or MUTSEL. Reject unsupported combinations clearly.

Rooted/nonreversible capability must be a property of the model, including at a reversible seed. Do not change to the reversible kernel just because current circulation happens to vanish; that could alter tree representation during optimization. A separate reversible-only model supplies the control.

### 8.2 Authoritative state

Maintain:

1. Immutable normalized target vector and amino-acid order.
2. Chart kind/version; jump reference destinations or T3 gauges.
3. Authoritative 360 coordinate values.
4. Derived normalized off-diagonal rates, diagonal, and numerical decomposition state.
5. Fixed/fitted flag, numerical-domain policy, tolerances and diagnostics.

Coordinates and target are authoritative for a fitted model; the normalized full $Q$ is authoritative for a fixed imported empirical model. Never silently reconstruct a different target from alignment counts or from an unconstrained proposal.

`getNDim()` returns 360 when unfixed and zero when fixed. Fixed external target adds zero to `getNDimFreq()`. The base routines' packing direction is easy to misread: `setVariables` copies model state into the optimizer array; `getVariables` reads optimizer values back. Override both. Do not copy 360 logits into the first 360 of 380 raw rates.

### 8.3 Generator installation and decomposition

Factor the existing nonreversible decomposition into generator construction/normalization and numerical decomposition, or add a narrow protected hook that allows an already validated generator to reach the existing decomposition. Prefer a maintainable shared numerical helper over copying the whole eigen/fallback implementation.

The new class constructs $Q$ from coordinates, installs consistent `rates[]` and `rate_matrix`, keeps `state_freq=pi_target`, validates invariants, and calls the shared numerical backend. Ordinary single-matrix `total_num_subst` is 1 in the inspected initialization; keep units consistent and reject unsupported mixture-scale overrides. A harmless normalization to correct roundoff is distinct from rebuilding a model at another target.

The fixed-target invariant must survive `optimize_from_given_params`; that flag means optimize the new class's coordinates, not release its equilibrium. Do not depend solely on the old `FREQ_USER_DEFINED` guard. Do not solve stationarity of the final $Q$ and replace the immutable target with a numerical approximation; use an independent equilibrium solve only as a diagnostic.

Check the constructor's virtual-call order and allocations before invoking derived packing/building methods. Ensure all 380 directed rates exist. On failed stationary/balancing solves, reject the trial without publishing partial state, or terminate with an explicit numerical diagnostic. In release builds use runtime checks, not assertions that may compile away.

### 8.4 Joint fitting and concurrency

Add explicit compatibility checks for state count/order, target values, chart version, references/gauges and numerical policy among linked objects. Existing name matching alone is insufficient. A mismatch must fail before optimization.

Initially let each partition reconstruct from the same immutable inputs, preserving the existing parallel objective. No shared mutable T3 warm-start state. Profile actual repeated builder/decomposition cost; an immutable shared prepared-matrix cache keyed by coordinates, target and numerical options is a later optimization, not necessary for correctness. Preserve the fixed-order objective reduction.

Honor `fixParameters(true)` outside the linked update so per-partition nuisance optimizers cannot drift to different generators. Test representative selection order. At a final accepted point synchronize every linked object and invalidate relevant partial/transition caches.

### 8.5 Input, candidates, and ModelFinder

Choose a maintainer-approved explicit option for target input; distinguish it from root-frequency-only semantics. Validate finite values, exact count/order, positive sum and strictly positive entries. Normalize once if input is expressed as weights, record that normalization, and never silently floor an exact requested target. If legacy `min_state_freq` prevents the target, reject with a precise limitation or implement a tested override; do not return an inexplicable constant objective.

For initial LG/JTT/WAG candidates, preserve exchangeabilities and rebuild at target. Disable competing empirical/optimized frequency variants in this protocol. Source candidate generation combines frequency suffixes independently, so target enforcement must be tested there too.

Later learned matrices are **fixed candidates**, not 360-parameter refits for each gene. Their evaluation-gene matrix dimension is zero, though training consumed parameters and may share data with reselection. Treat in-training candidate BIC as a heuristic guide-selection criterion, not independent generalization evidence.

To handle reversible and nonreversible candidates safely, first use separate supported candidate evaluations and combine comparable scores under the same declared root and nuisance policy. Alternatively build a tested common rooted candidate runner. Reversible root split redundancy must be accounted for in scoring. User-defined matrix names may evade current name-array guards; that is not evidence that mixed handling is correct. Do not bypass a guard by renaming models.

### 8.6 Persistence and output

Checkpoint target, coordinates, chart metadata, fixed/fitted status, bounds/policy, accepted likelihood, and version. Restore into a fresh object, reconstruct $Q$, invalidate caches, and verify invariant/likelihood agreement. A checkpoint with another target is incompatible, not an invitation to overwrite the current target.

Export the full normalized matrix with high precision, target vector and amino-acid order. Preserve a lossless machine-readable fitted-model checkpoint separately from a rounded human report. Test fixed-matrix import and likelihood equality. Export/import of sparse boundary matrices, if supported later, needs its own policy and target validation.

Report conditional matrix dimension 360, exact target provenance, equilibrium and row-sum residuals, normalization, minimum flux/rate, numerical restrictions, derivative/solver status, seed/run identity, attained likelihood and nuisance policy. Separate invariant cycle affinities from raw pairwise tilts and antisymmetric flux. Update model-name-dependent citation/reporting logic: current QMaker/nQMaker attribution checks the substring `NONREV` rather than a model capability.

### 8.7 Pseudocode: conditional shared matrix update

```text
input: identical target/chart definition across linked models,
       accepted Q and coordinates, fixed nuisance snapshot
validate linked compatibility and Q feasibility
save accepted coordinates, Q, likelihood, and caches/rebuild recipe

objective(theta):
    build positive normalized target-stationary Q(theta)
    reject numerical failure; do not clip entries independently
    install consistently in each linked model
    invalidate dependent caches
    return negative sum of current gene log-likelihoods

run existing BFGS through the new scaled derivative policy
restore/evaluate returned theta after the last derivative probe
check physical feasibility and first-order diagnostics
accept only a feasible nondecreasing result; otherwise retain incumbent
record convergence, numerical bounds and any failure
refreeze linked matrices before nuisance optimization
```

Retaining the incumbent prevents decreases but does not turn a failed optimization into a converged fit. Return a failure/partial status when warranted.

## 9. EM alternative: retain mathematically, defer integration

### 9.1 Hidden-history criterion — derived

If complete substitution histories were observed, at fixed branch/rate nuisance values the $Q$-dependent log-likelihood would be

$$
\sum_{i\ne j}N_{ij}\log q_{ij}-\sum_i T_i\sum_{j\ne i}q_{ij}+\mathrm{constant}.
$$

Each jump contributes its log intensity; intervals in state $i$ contribute the negative exit-rate exposure. The root term is constant because its distribution is fixed. In the E step replace counts and exposures by their posterior expectations. Exposure incorporates site/partition rate multipliers: for category multiplier $r$, use $r$ times physical dwell time.

At fixed target, let $d_{ij}=\overline T_i/\pi_i^\star$. The normalized M step maximizes

$$
M(f)=\sum_e\overline N_e\log f_e-d^Tf,\qquad Cf=b,\ f\ge0.
$$

Its Hessian is $-\operatorname{diag}(\overline N_e/f_e^2)$. With all expected counts positive, it is strictly concave, diverges to $-\infty$ at the boundary, and has a unique interior maximizer on the compact closure. With zero counts, strict concavity and interior attainment require separate analysis; do not reuse the positive-count formulas blindly.

### 9.2 Twenty-dimensional normalized dual — derived

For the maximization Lagrangian $M(f)-\lambda^T(Cf-b)$, stationarity gives

$$
f_e=\frac{\overline N_e}{d_e+(C^T\lambda)_e}.
$$

Substitution gives the dual, up to constants,

$$
\psi(\lambda)=b^T\lambda-\sum_e\overline N_e\log[d_e+(C^T\lambda)_e].
$$

It is minimized over positive denominators, with

$$
\nabla\psi=b-Cf,\qquad
\nabla^2\psi=C\operatorname{diag}\left(\frac{\overline N_e}{[d_e+(C^T\lambda)_e]^2}\right)C^T.
$$

The Hessian is positive definite because $C$ has full row rank and the diagonal is positive. The primal optimum is interior; its gradient is orthogonal to $\ker C$, so a multiplier exists by elementary linear algebra, yielding the displayed denominator equation and a dual optimum. This proves the primal/dual correspondence here without assuming that an arbitrary Newton iteration converges. The 20 multipliers represent 19 balance equations and one normalization equation, not biological frequencies.

Use domain-preserving damped Newton with verified residuals. Active floors require the appropriate modified convex solve. This M step keeps branch lengths fixed, matching paper substep 3a most directly.

### 9.3 Nineteen-dimensional cone alternative

Omit normalization and optimize balanced positive flux. With positive counts and positive state exposures, the logarithmic boundary divergence and linear penalty at infinity give a unique primal optimum. There are only 19 independent balance multipliers. If the resulting raw generator has mean rate $m$, set $Q_{\rm new}=Q_{\rm raw}/m$ and $t_{b,\rm new}=m\,t_{b,\rm old}$. Then $t_{b,\rm new}Q_{\rm new}=t_{b,\rm old}Q_{\rm raw}$ and transitions are unchanged.

This is a legitimate different conditional update when the compensating branch rescaling is permitted and stays within bounds. It is not the same fixed-branch M step. Normalizing $Q$ alone changes transitions and forfeits that argument. For the desired architecture, prefer the normalized 20-multiplier version if EM is implemented.

### 9.4 E step and likelihood derivatives

For a branch with endpoints $a,b$, effective time $t$ and $P(t)=e^{tQ}$,

$$
E[T_i\mid a,b]=\frac{\int_0^tP_{ai}(s)P_{ib}(t-s)\,ds}{P_{ab}(t)},\quad
E[N_{ij}\mid a,b]=\frac{q_{ij}\int_0^tP_{ai}(s)P_{jb}(t-s)\,ds}{P_{ab}(t)}.
$$

These follow by conditioning on occupation or a jump at time $s$, integrating and dividing by the endpoint probability. An outside-message pass plus pruning supplies posterior endpoint and category weights. Average these bridge expectations using those weights and pattern multiplicities.

The upper-right block of $\exp[t\begin{pmatrix}Q&E\\0&Q\end{pmatrix}]$ is $\int_0^t e^{(t-s)Q}Ee^{sQ}ds$, from the block evolution equation. Real block exponentials provide a reference implementation; complex arithmetic is not mathematically necessary. Production should exploit shared structure, not separately exponentiate a block for every edge/pattern/parameter.

Differentiating the likelihood sum/integral under the usual finite positive-model regularity gives the posterior-score identity $\partial\ell/\partial q_{ij}=\overline N_{ij}/q_{ij}-\overline T_i$, with diagonals completed. Thus E-step infrastructure can also supply direct gradients.

### 9.5 What EM guarantees

Let $r(Z)=p(Z\mid D,\theta_{\rm old})$. Jensen gives

$$
\ell(\theta_{\rm new})-\ell(\theta_{\rm old})\ge
E_r[\log p(D,Z\mid\theta_{\rm new})-\log p(D,Z\mid\theta_{\rm old})].
$$

An exact or increasing feasible M update is therefore likelihood nondecreasing. It does not guarantee the global observed-likelihood optimum or a universal convergence rate. After independent nuisance changes, recompute expectations unless a separately justified ECM scheme covers the updates. Expected-history machinery is a substantial integration change absent from the inspected direct matrix-update path; existing EM elsewhere in IQ-TREE is not proof that this E step already exists.

## 10. Scientific validity and inference safeguards

### 10.1 Stationary marginals cannot be changed by circulation — derived

Since $\pi^\star Q=0$, the exponential series gives $\pi^\star e^{tQ}=\pi^\star$. Starting the root at this vector makes every modeled tip marginal equal to it, also under scalar rate mixtures. If the true retained-data marginal $p_a$ differs, all candidate generators are misspecified. Factoring a joint distribution into this marginal and its conditional distribution gives

$$
\mathrm{KL}(P_0\Vert P_Q)=\mathrm{KL}(p_a\Vert\pi^\star)+
E_{p_a}\mathrm{KL}(P_0(\cdot\mid a)\Vert P_Q(\cdot\mid a))
\ge\mathrm{KL}(p_a\Vert\pi^\star).
$$

A finite sample's empirical frequency difference alone does not establish population misspecification. Under actual mismatch, structural parameters can improve pattern correlations but cannot repair the marginal mismatch. A nonreversible fitted advantage can therefore reflect a better misspecified approximation, not necessarily biological time irreversibility. Appendix A gives newly computed numerical evidence of this possibility, not a universal claim that it occurs.

Cleaning can select on variability, gaps or structural context as well as composition. A fixed equilibrium does not model that selection mechanism. For a simple retention rule, a selection-aware probability would involve $p(D,\mathrm{retained})/p(\mathrm{retained})$; general alignment cleaners can require more complicated conditioning. Standard ascertainment correction should not be presumed to cover an arbitrary cleaner.

### 10.2 Root edge, root split and tests — derived

Let $\Pi=\operatorname{diag}(\pi^\star)$ and $Q^\leftarrow=\Pi^{-1}Q^T\Pi$. For two children of the root,

$$
J(t_u,t_v)=e^{t_uQ^T}\Pi e^{t_vQ}=\Pi e^{t_uQ^\leftarrow}e^{t_vQ}.
$$

With total length fixed, its split derivative is

$$
\frac{dJ}{dt_u}=\Pi e^{t_uQ^\leftarrow}(Q^\leftarrow-Q)e^{t_vQ}.
$$

All outer factors are invertible. The derivative is zero exactly when detailed balance holds. Under reversibility $J=\Pi e^{(t_u+t_v)Q}$, so the split is unidentifiable. Fixing a root edge alone does not fix this nuisance redundancy. The derivative proves local sensitivity of the root-child joint law at known $Q$, not global or joint identification of all parameters from leaves.

For example, the three-state cyclic family in Section 6 gives a joint law whose split dependence is periodic with phase $\sqrt3 x(t_v-t_u)/2$. Distinct positive splits separated by $2\pi/(\sqrt3|x|)$ in $t_u$ can agree exactly. Unknown $x$ and split can also trade off at fixed product $x(t_v-t_u)$ on a two-tip tree.

Dimension differences 171 and 19 do not automatically supply chi-square LRT calibrations. Correct specification, local identification, interiority, nuisance regularity and adequate optimization must hold. A reversible truth may make root splitting singular even for the fixed-versus-free equilibrium comparison. Use a justified common root/split protocol or bootstrap the actual estimation procedure for inferential claims; distinguish a statistical fitted-model bootstrap from simulations testing artifacts caused by cleaning and deliberately wrong composition.

### 10.3 Experiments

| Fit / application | Purpose |
|---|---|
| Uncleaned-data reversible and unconstrained general fits | Reference working models |
| Cleaned-data unconstrained general fit | Baseline with composition free through $Q$ |
| Historical reversible fit followed by target swap | Preserve the old workflow as a distinct comparator |
| Reversible refit at target | Matched constrained reversible baseline |
| General refit at target | Project 2 main estimator |
| General cleaned fit transferred by flux and/or T3, without refitting | Separate transport from refitting; optional sensitivity arms for the larger project |

Match training/evaluation genes, target, rate-family/root policy and computational budget as appropriate. A raw likelihood comparison across differently cleaned datasets is not a comparison on identical observations. Evaluate frozen learned matrices on the same held-out data; let only the declared tree/rate nuisances change. A reference species tree is a benchmark, not certain gene-tree truth.

Simulate reversible and nonreversible generators; correct and deliberately perturbed targets; well-behaved and near-boundary rates; outcome-independent masks and realistic content-dependent cleaning; and relevant heterogeneity. Include true target, independently estimated target and training plug-in target arms. A reversible refit cannot itself reveal spurious circulation because it forbids circulation.

If uncertainty includes target estimation, resample at the scientifically appropriate unit and recompute both target and fit. A BFGS inverse Hessian is not automatically a covariance after nuisance fitting, misspecification or tree selection. Even smooth profiling changes the relevant curvature by a nuisance Schur complement.

### 10.4 Conditional asymptotic interpretation — derived under explicit assumptions

On a compact fixed model domain, suppose normalized log-likelihoods converge uniformly to a continuous population criterion $L$, and attained fits are within $o_p(1)$ of the normalized global supremum. If $L$ has a unique maximizer $Q_0$, the fits converge to it: outside any neighborhood of $Q_0$, compactness and uniqueness give a strictly positive gap from the maximum; uniform approximation and vanishing optimization error eventually prevent an approximate maximizer from remaining there. For a nonunique maximizing set, the same argument gives convergence to that set when it is separated from its complement outside each neighborhood.

For independent patterns from a fixed population distribution, expected log probability equals a constant minus KL divergence. The population maximizer is therefore a KL-minimizing working model, not necessarily the generating process under misspecification. These conclusions require the uniform approximation, domain, and optimization hypotheses; they do not establish them for heuristic local fits, growing numbers of gene-specific nuisances, dependent sites, target vectors estimated from the same data, or topology selection. Such conditions must be demonstrated separately before claiming consistency for the full pipeline.

## 11. Development milestones and acceptance gates

### Milestone 1: standalone mathematical core

Implement target validation, jump build/inverse, T3 build/inverse, flux conversion and common residuals. Prove dimension and coverage as above; execute random and adversarial checks for $n=3,4,20$. Reversible reduction and arbitrary balanced-cycle-generated flux round trips must pass. The prototype in Appendix B covers moderate cases only.

### Milestone 2: one fixed-tree IQ-TREE model

Add class/factory parsing, 360-variable packing, frozen target, decomposition integration, scaled finite differences and lifecycle tests. Compare a fixed reversible $Q$ under the existing reversible kernel and the new general class with matched root/branch treatment. Compare the same general $Q$ against independent exponentials and likelihood calculations. Test `optimize_from_given_params`, invalid targets and sub-threshold targets explicitly.

### Milestone 3: linked training

Test two partitions with intentionally different empirical compositions but identical target. Confirm one shared $Q$, 360 counted once, fixed matrix during nuisance updates and identical target after restart. Cover both separate trees and edge-proportional branches as intended. Check thread-count/order stability, finite-difference state restoration and interrupted fits. Do not require universal bitwise equality across different compilers; require justified numerical equivalence, with determinism where the same environment permits it.

### Milestone 4: optimization comparison

From the same feasible $Q$, fixed data and nuisance snapshot, compare jump, T3 and flux implementations. Require matched effective domains or explicit representability checks. Report best likelihood, physical first-order residual, feasibility, failed trials, expensive likelihood evaluations, total runtime and repeated-start variability. Extend to 20-state likelihood-level examples before choosing a default. Tiny three-state fits do not settle performance.

The full fixed-target model contains the reversible target model. Retain its fitted $Q$ as an admissible incumbent, making a lower reported final score a protocol/numerical failure. The full unconstrained model contains the constrained family mathematically; verify actual input/bounds representability before using attained-score ordering as a diagnostic.

### Milestone 5: workflow/export integration

Test target-consistent seed candidate construction, disabled frequency variants, fixed learned candidates, mixed-candidate adapter, root handling, cache-safe export/reimport, native nuisance alternation and optional outer reselection. Correct the paper's self-correlation ordering and use likelihood/convergence diagnostics in addition to matrix correlation. Native tree/branch algorithms remain in use; only their declared scheduling and model invariants are adapted.

### Milestone 6: numerical boundaries and science

Stress tiny targets, highly skewed jump matrices, large cycle affinities, nearly defective/non-normal generators and sparse limits. Compare matrix-exponential methods and full gradients. Expand boxes or lower floors and report sensitivity. Then run held-out and cleaning/mismatch simulations from Section 10, with no test-set method tuning.

### Suggested diagnostic definitions

For unit-rate $Q$, record $r_{\rm row}=\|Q\mathbf1\|_\infty/\max(1,\|Q\|_\infty)$, $r_{\rm stat}=\|\pi^\star Q\|_\infty$ (equivalently balance of normalized flux), and $r_{\rm rate}=|{-\pi^\star\operatorname{diag}Q}-1|$. Inspect absolute and relative forms for extreme targets rather than allowing a large matrix norm to hide flux imbalance. A provisional moderate-case double-precision feasibility target is around $10^{-10}$, tightened or qualified using conditioning and likelihood impact. Derivative tolerances must be scaled and established by a step-size plateau, not a single universal threshold.

## 12. Consolidated adjudication of the supplied reports

| Earlier issue | Independent conclusion retained here |
|---|---|
| Full fixed-target dimension 360; reversible 189; circulation 171 | Correct, with full incidence-rank proof |
| Both charts cover the target family | Correct for positive matrices with gauges/reference conventions specified |
| T3 must be first because it reduces reversibly and transports well | Not forced for Project 2; jump-first prototype better matches the minimal-new-solver priority, benchmark before production choice |
| Initial feasible candidates may be the hardest part | Unlikely for LG/JTT/WAG; simple reversible rebuilding suffices |
| Initial guide candidates must always obey target | Not logically necessary if used only as guides; a target-consistent workflow is the preferred declared policy |
| Fixing `state_freq` enforces equilibrium | False; construct constrained rates too |
| Everything except initial candidates and step 3a can be untouched | Correct about much reusable mathematics, incomplete about frequency variants, later candidates, roots, gradients, linking and persistence |
| Both audits establish that no earlier mathematical error remains | False; the cycle-bound claim is contradicted explicitly |
| Positive shifted coordinates are mandatory | False; repair derivative scaling or use verified analytic gradients |
| Flux elimination is an incomplete model | False on its correct feasible domain; it needs coupled-inequality handling |
| Convex feasible space or inner problem makes ML globally easy | False; observed likelihood is nonconcave |
| Nineteen- and twenty-variable EM updates are identical | False; normalized fixed-branch and scale-releasing updates differ |
| Wrong target inevitably creates circulation | Not established; it can occur, and needs simulation rather than a universal claim |
| T3 preserves every raw circulation measure | False; cycle affinities are invariant under its gradient adjustment, raw tilts/fluxes need not be |
| Nonzero root-split derivative proves global identifiability | False; it proves local sensitivity of a particular joint distribution at known $Q$ |
| Interior fitted points imply nesting of distinct bounded charts | False; nesting is a property of the feasible sets |
| Plug-in target adds exactly 19 BIC parameters | Not automatically justified for the two-stage estimator |
| Complex arithmetic is required for EM integrals | False; real block exponentials suffice as an oracle |
| A finite penalty can never yield exact stationarity | Too strong; it does not guarantee stationarity, though special optima can satisfy it |
| Small balancing solve always negligible | Unmeasured, especially with repeated per-partition reconstruction |

## 13. Remaining decisions, without blocking the design

The implementation plan above is concrete without choosing a new biological transport rule. Before scientific runs, record: target provenance and smoothing; tree-sharing design; stationary-root and root-search policy; selected rate families and which continuous parameters move; boundary/numerical policy; initialization budget; convergence tolerances; and evaluation split.

The maintainers should decide the public option/name, whether to expose both charts, and whether the paper-style outer driver belongs in the first release. Defaults proposed here are one fixed positive target, stationary root, no regularization, jump-chart prototype, native IQ-TREE nuisance scheduling, and a separate fixed-nuisance validation mode. These choices should be adjusted by evidence, not silently changed during comparisons.

## Appendix A. Verification performed for this document

Direct source inspection and GitHub commit resolution were completed. No IQ-TREE executable was compiled or run, no production patch was written, and no empirical training alignment was fitted. End-to-end flag behavior, full 20-state likelihood performance, all unsupported model interactions and scientific predictive gains remain unverified.

The independent script in Appendix B uses balanced directed-cycle mixtures to generate test matrices rather than generating every test through the chart being checked. It tests round trips, derivatives, reversible reduction, constraint rank, a normalized EM dual against a separately constrained primal, three chart fits to one toy population distribution, a nonconcavity identity, the invalid T3 bound, a Fréchet derivative, root splitting, and an example of composition-induced nonreversible improvement.

Numerical results and the runnable source follow. This was not an IQ-TREE test. Near-boundary coverage, error bounds and global optimizer correctness are **not** certified by these moderate examples.

```json
{
  "numpy": "2.3.5",
  "scipy": "1.17.0",
  "seed": 22092026,
  "dimension_3": 3,
  "dimension_4": 8,
  "dimension_20": 360,
  "jump_roundtrip": 1.6431300764452317e-14,
  "t3_roundtrip": 1.0480505352461478e-13,
  "stationarity": 5.4817261840867104e-15,
  "rows": 2.168404344971009e-15,
  "mean_rate": 8.881784197001252e-16,
  "jump_derivative": 6.338511197591968e-10,
  "t3_derivative": 1.0473954327054301e-10,
  "reversible_reduction": 7.105427357601002e-15,
  "nonconcavity_curvature": 0.015377663244794715,
  "nonconcavity_formula_error": 2.220446049250313e-16,
  "triangle_half_affinity": 20.72326583694641,
  "em_balance_error": 1.1102230246251565e-16,
  "em_primal_dual_flux_error": 2.18899448671539e-08,
  "toy_KL_t3_jump_flux": [
    [
      4.8051099602536e-17,
      1.887817378951202e-15,
      1.2725203334019324e-15
    ],
    [
      5.742431483474016e-17,
      7.837948137793633e-16,
      7.812972982161087e-16
    ]
  ],
  "frechet_error": 6.074959979507355e-12,
  "root_split_derivative_error": 3.1143385396600887e-12,
  "mismatch_example": {
    "reversible_KL": 0.06448436727458651,
    "general_KL": 0.06447760974027678,
    "general_cycle_coordinate": -0.06728303245062825,
    "cycle_score_at_reversible_fit": 0.00020064367063643562
  }
}
```

## Appendix B. Reproducible independent check code

Requires Python with NumPy and SciPy. This is verification code, not production C++, and does not implement the full IQ-TREE likelihood. The toy population fits use exact site-pattern probabilities on a three-tip star with fixed branch lengths. The composition example is numerical evidence at attained fits; it is not an interval-certified global optimum or an estimate of artifact frequency.

```python
import json
import numpy as np
import scipy
from scipy.linalg import expm, expm_frechet, null_space
from scipy.optimize import minimize
from scipy.special import softmax

rng = np.random.default_rng(22092026)
results = {'numpy': np.__version__, 'scipy': scipy.__version__, 'seed':22092026}

def geometry(n):
    i,j=np.where(~np.eye(n,dtype=bool))
    B=np.eye(n)[:,j]-np.eye(n)[:,i]
    C=np.vstack((B[:-1],np.ones(len(i))))
    return i,j,B,C

def generator(f,pi):
    q=f/pi[:,None]
    np.fill_diagonal(q,-q.sum(axis=1))
    return q

def stationary(q):
    n=len(q)
    return np.linalg.lstsq(np.vstack((q.T,np.ones(n))),np.r_[np.zeros(n),1.],rcond=None)[0]

def jump(z,pi,dz=None):
    n=len(pi); P=np.zeros((n,n)); dP=np.zeros_like(P)
    for i in range(n):
        js=[j for j in range(n) if i!=j]
        P[i,js]=softmax(np.r_[z[i],0.])
        if dz is not None:
            v=np.r_[dz[i],0.]
            dP[i,js]=P[i,js]*(v-P[i,js]@v)
    nu=stationary(P-np.eye(n))
    f=nu[:,None]*P
    q=generator(f,pi)
    if dz is None:return q
    dnu=np.linalg.solve((np.eye(n)-P+np.ones((n,1))*nu).T,(nu@dP))
    df=dnu[:,None]*P+nu[:,None]*dP
    return q,generator(df,pi)

def jump_inverse(q):
    n=len(q); P=q/(-np.diag(q))[:,None]
    return np.array([np.log(P[i,[j for j in range(n) if j!=i]][:-1]/P[i,[j for j in range(n) if j!=i]][-1]) for i in range(n)])

def t3(a,h,pi,da=None,dh=None):
    n=len(pi);i,j,B,C=geometry(n);D=B[:-1]
    z=np.log(pi[i])+np.log(pi[j])+a[i,j]+h[i,j]
    z-=z.max(); g=np.zeros(n-1)
    for iteration in range(100):
        x=np.exp(z+D.T@g); residual=D@x
        if np.max(abs(residual))/x.sum()<2e-14: break
        H=(D*x)@D.T; step=np.linalg.solve(H,-residual)
        alpha=1.
        while np.exp(z+D.T@(g+alpha*step)).sum()>x.sum()+1e-4*alpha*(residual@step)+2e-15*x.sum():
            alpha*=.5
            if alpha<1e-12:raise RuntimeError('balancing stalled')
        g+=alpha*step
    else: raise RuntimeError('balancing exhausted')
    f=np.zeros((n,n)); f[i,j]=x/x.sum();q=generator(f,pi)
    if da is None:return q
    x=f[i,j];v=da[i,j]+dh[i,j];H=(D*x)@D.T
    dg=np.linalg.solve(H,-D@(x*v));w=v+D.T@dg
    df=np.zeros_like(f);df[i,j]=x*(w-x@w)
    return q,generator(df,pi)

def t3_inverse(q,pi):
    n=len(pi);i,j,_,_=geometry(n)
    a=np.zeros_like(q);h=np.zeros_like(q)
    a[i,j]=.5*np.log(q[i,j]*q[j,i]/(pi[i]*pi[j]))
    h[i,j]=.5*np.log(pi[i]*q[i,j]/(pi[j]*q[j,i]))
    a-=a[n-2,n-1];np.fill_diagonal(a,0.)
    u=h[:,-1].copy();h+=u[None,:]-u[:,None]
    return a,h

maxerr={k:0. for k in ['jump_roundtrip','t3_roundtrip','stationarity','rows','mean_rate','jump_derivative','t3_derivative','reversible_reduction']}
for n in [3,4,20]:
    i,j,B,C=geometry(n)
    assert np.linalg.matrix_rank(C)==n
    results[f'dimension_{n}']=len(i)-np.linalg.matrix_rank(C)
    for rep in range(6):
        # Generate balanced flux independently, as a sum of directed cycles.
        f=np.zeros((n,n))
        for u in range(n):
            for v in range(u+1,n):f[u,v]=f[v,u]=rng.lognormal(0,1)
        for _ in range(3*n):
            cycle=rng.choice(n,size=min(n,3),replace=False);w=rng.lognormal()
            for u,v in zip(cycle,np.roll(cycle,-1)):f[u,v]+=w
        f/=f.sum();pi=rng.dirichlet(np.ones(n)*2);q=generator(f,pi)
        a,h=t3_inverse(q,pi);z=jump_inverse(q)
        for key,q2 in [('jump_roundtrip',jump(z,pi)),('t3_roundtrip',t3(a,h,pi))]:
            maxerr[key]=max(maxerr[key],float(np.max(abs(q2-q))))
            maxerr['stationarity']=max(maxerr['stationarity'],float(np.max(abs(pi@q2))))
            maxerr['rows']=max(maxerr['rows'],float(np.max(abs(q2.sum(1)))))
            maxerr['mean_rate']=max(maxerr['mean_rate'],float(abs(-pi@np.diag(q2)-1)))
        dz=rng.normal(size=z.shape);eps=1e-5
        _,dq=jump(z,pi,dz)
        fd=(jump(z+eps*dz,pi)-jump(z-eps*dz,pi))/(2*eps)
        maxerr['jump_derivative']=max(maxerr['jump_derivative'],float(np.max(abs(fd-dq))/max(1,np.max(abs(dq)))))
        da=rng.normal(size=(n,n));da=(da+da.T)/2
        dh=rng.normal(size=(n,n));dh=(dh-dh.T)/2
        _,dq=t3(a,h,pi,da,dh)
        fd=(t3(a+eps*da,h+eps*dh,pi)-t3(a-eps*da,h-eps*dh,pi))/(2*eps)
        maxerr['t3_derivative']=max(maxerr['t3_derivative'],float(np.max(abs(fd-dq))/max(1,np.max(abs(dq)))))
        rev=np.exp(a)*pi[None,:];np.fill_diagonal(rev,0);np.fill_diagonal(rev,-rev.sum(1));rev/=(-pi@np.diag(rev))
        maxerr['reversible_reduction']=max(maxerr['reversible_reduction'],float(np.max(abs(rev-t3(a,np.zeros_like(a),pi)))))
results.update(maxerr)
assert max(maxerr.values())<1e-7

def cyclic(x):
    return np.array([[-1,(1+x)/2,(1-x)/2],[(1-x)/2,-1,(1+x)/2],[(1+x)/2,(1-x)/2,-1.]])
x=2*np.pi/(5*np.sqrt(3));A=2*np.exp(-1.5*5.2);k=5*np.sqrt(3)/2
results['nonconcavity_curvature']=float(A*k*k/(1-A))
observed=(expm(.1*cyclic(x))[:,0]@expm(5.1*cyclic(x))[:,0])/3
results['nonconcavity_formula_error']=abs(observed-(1+A*np.cos(k*x))/9)
assert results['nonconcavity_curvature']>0

# Counterexample to the proposed T3 cycle bound of 12.
q=np.ones((20,20));np.fill_diagonal(q,0)
q[0,1]=q[1,2]=q[2,0]=100
q[1,0]=q[2,1]=q[0,2]=1e-4
np.fill_diagonal(q,-q.sum(1))
hh=.5*np.log(q[0,1]*q[1,2]*q[2,0]/(q[1,0]*q[2,1]*q[0,2]))
results['triangle_half_affinity']=float(hh);assert hh>12

# Normalized EM dual checked against an independently constrained primal.
n=4;i,j,B,C=geometry(n);b=np.r_[np.zeros(n-1),1.];base=np.ones(len(i))/len(i);Z=null_space(C)
N=rng.uniform(.2,5,len(i));d=rng.uniform(.1,3,len(i))
lam=np.r_[np.zeros(n-1),N.sum()]
def dual(l):
    den=d+C.T@l
    return b@l-N@np.log(den) if np.all(den>0) else np.inf
for _ in range(100):
    den=d+C.T@lam;f=N/den;g=b-C@f
    if np.max(abs(g))<1e-12:break
    H=(C*(N/den**2))@C.T;step=np.linalg.solve(H,-g);alpha=1.
    while dual(lam+alpha*step)>dual(lam)+1e-4*alpha*(g@step)+1e-13:alpha*=.5
    lam+=alpha*step
f=N/(d+C.T@lam)
def primal(z):
    ff=base+Z@z
    return -N@np.log(ff)+d@ff if np.all(ff>0) else 1e50
pr=minimize(primal,np.zeros(Z.shape[1]),method='SLSQP',constraints={'type':'ineq','fun':lambda z:base+Z@z-1e-10},options={'ftol':1e-12,'maxiter':1000})
results['em_balance_error']=float(np.max(abs(C@f-b)))
results['em_primal_dual_flux_error']=float(np.max(abs(f-base-Z@pr.x)))
assert results['em_balance_error']<1e-9 and results['em_primal_dual_flux_error']<1e-6

# Same observed-data likelihood with three complete coordinate representations.
pi=np.array([.23,.31,.46]);bl=[.17,.62,1.23]
def small_t3(v,pi=pi):
    a=np.array([[0,v[0],v[1]],[v[0],0,0],[v[1],0,0.]])
    h=np.array([[0,v[2],0],[-v[2],0,0],[0,0,0.]])
    return t3(a,h,pi)
def probs(q):
    return np.einsum('i,ia,ib,ic->abc',pi,*[expm(t*q) for t in bl]).ravel()
truth=probs(small_t3([.6,-.4,.45]))
def kl(q):
    p=probs(q)
    return float(truth@np.log(truth/p)) if np.all(p>0) else 1e50
i,j,B,C=geometry(3);Z=null_space(C);f0=np.ones(6)/6
def flux_q(z):
    f=np.zeros((3,3));f[i,j]=f0+Z@z;return generator(f,pi)
fit=[]
for _ in range(2):
    q0=small_t3(rng.normal(size=3)*.1)
    a0,h0=t3_inverse(q0,pi);v0=[a0[0,1],a0[0,2],h0[0,1]]
    z0=jump_inverse(q0).ravel();u0=Z.T@((pi[:,None]*q0)[i,j]-f0)
    aa=minimize(lambda v:kl(small_t3(v)),v0,method='BFGS',options={'gtol':1e-9})
    bb=minimize(lambda z:kl(jump(z.reshape(3,1),pi)),z0,method='BFGS',options={'gtol':1e-9})
    cc=minimize(lambda z:kl(flux_q(z)),u0,method='SLSQP',constraints={'type':'ineq','fun':lambda z:f0+Z@z-1e-9},options={'ftol':1e-13,'maxiter':500})
    fit.append([aa.fun,bb.fun,cc.fun])
results['toy_KL_t3_jump_flux']=fit
assert np.max(np.abs(fit))<1e-9

q=small_t3([.3,-.2,.4]);_,E=jump(jump_inverse(q),pi,rng.normal(size=(3,1)))
exact=expm_frechet(.7*q,.7*E,compute_expm=False)
fd=(expm(.7*(q+1e-5*E))-expm(.7*(q-1e-5*E)))/(2e-5)
results['frechet_error']=float(np.max(abs(exact-fd)))
Pi=np.diag(pi);qr=np.linalg.solve(Pi,q.T@Pi);u,v=.4,.9
joint=lambda s,t:expm(s*q).T@Pi@expm(t*q)
deriv=Pi@expm(u*qr)@(qr-q)@expm(v*q)
fd=(joint(u+1e-5,v-1e-5)-joint(u-1e-5,v+1e-5))/(2e-5)
results['root_split_derivative_error']=float(np.max(abs(deriv-fd)))
assert max(results['frechet_error'],results['root_split_derivative_error'])<1e-8
# Population mismatch example: numerical evidence, not a global certificate.
pi0=np.array([.60,.28,.12]);pit=np.array([.49,.30,.21])
def star(q,root):
    return np.einsum('i,ia,ib,ic->abc',root,*[expm(t*q) for t in [.2,.7,1.1]]).ravel()
truth0=star(small_t3([.7,-.5,0],pi0),pi0)
def mismatch(v):
    p=star(small_t3(v,pit),pit)
    return float(truth0@np.log(truth0/p))
rev=minimize(lambda v:mismatch([*v,0]),[.7,-.5],method='BFGS',options={'gtol':1e-10})
full=minimize(mismatch,[*rev.x,0],method='BFGS',options={'gtol':1e-10})
eps=1e-4
score=(mismatch([*rev.x,eps])-mismatch([*rev.x,-eps]))/(2*eps)
results['mismatch_example']={'reversible_KL':rev.fun,'general_KL':full.fun,'general_cycle_coordinate':full.x[2],'cycle_score_at_reversible_fit':score}
assert rev.fun-full.fun>1e-8
print(json.dumps(results,indent=2,default=lambda v:v.item()))
```

## Primary sources and input provenance

- Dang et al. (2022), *nQMaker: Estimating Time Nonreversible Amino Acid Substitution Models*, supplied `syac007(1).pdf`, especially pp. 1111–1113. [Publisher DOI](https://doi.org/10.1093/sysbio/syac007). Workflow claims here use the supplied PDF, not unverified descriptions of subsequent corrections.
- IQ-TREE3 [pinned source tree](https://github.com/iqtree/iqtree3/tree/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de), fetched from the requested master branch on the review date.
- IQ-TREE [official model-estimation tutorial](https://iqtree.github.io/doc/Estimating-amino-acid-substitution-models), consulted for its user-facing separate-alignment versus concatenated examples. Source controls implementation claims; tutorial wording is not treated as proof of internals.
- Input guidance, not sources of proved claims: `project2_constrained_estimation_methods(1).md`, `nq_t4_constrained_estimation(1).md`, `project2_comparative_audit(1).md`, `nq_t4_two_agent_comparison.md`, `nq_problem_contract_v1(3).md`, `AA_MODEL_INFERENCE(2).md`.

No novelty claim is made for Markov-chain coordinates, balancing, EM or constrained likelihood. The derivations above stand on their stated assumptions; publication-level historical attribution needs a dedicated primary-literature review.

[S1]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelprotein.cpp
[S2]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelmarkov.cpp
[S3]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelmarkov.h
[S4]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelfactory.cpp
[S5]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/partitionmodel.cpp
[S6]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/partitionmodelplen.cpp
[S7]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/optimization.cpp
[S8]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phylotesting.cpp
[S9]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylotree.cpp
[S10]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelunrest.cpp
[S11]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylokernelnonrev.h
[S12]: https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp
