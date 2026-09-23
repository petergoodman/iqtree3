# Project 2: unified discrepancy synthesis and implementation decision record

**Version:** 1.0 — 23 September 2026  
**Audience:** agents and developers implementing or reviewing fixed-equilibrium nonreversible amino-acid estimation in IQ-TREE3  
**Status:** current working design; not a completed implementation, maintainer approval, or claim of global maximum-likelihood estimation

## 0. How to use this companion

Read this document alongside both technical plans. It resolves their conflicts, preserves useful alternatives, specifies what evidence supports the choices, and defines the conditions under which the choices should be revisited. It does not replace their detailed derivations.

Use these unambiguous aliases throughout the handoff:

| Alias | Supplied document | Original principal recommendation |
|---|---|---|
| **P-positive** | `constrained_nq_plan_agent(1).md` | Jump-chain positive weights with a pinned entry per row; retain existing optimizer conventions |
| **P-log** | `project2_iqtree3_technical_plan.md` | Jump-chain log-ratios; coordinate-appropriate derivatives and bounds; broader integration specification |
| **D-review** | `discrepancy_decisions.md` | Later adjudication favoring log-ratios, supported by a Python optimizer comparison |
| **A-differences** | `idenfitied_agent_plan_differences_a.md` | Broader discrepancy list and mathematical corrections |
| **B-differences** | `identified_agent_plan_differences_b.md` | Shorter discrepancy list and corrected seed interpretation |
| **Source guide** | `AA_MODEL_INFERENCE(3).md` | Leads for source inspection, not an authoritative specification |

**Precedence for this handoff:** later explicit user decisions take priority. Otherwise use this synthesis for the conflicts it addresses; use the two plans for compatible detail. A new mathematical contradiction or source mismatch must be documented and resolved, not concealed by invoking document precedence. This record does not certify every ancillary claim in either plan. Do not revive a superseded recommendation merely because it appears repeatedly in an older document.

The user wants one fixed-target estimation project, consistency with IQ-TREE's established algorithms where appropriate, mathematical correctness, and retention of credible alternatives. The user specifically favors T3 as a serious comparison method and does not want positive ratios discarded. The log-ratio choice below is the current synthesis recommendation, not a claim that either coordinate representation is inherently superior.

### Evidence labels

- **MATH:** algebraic result under the stated assumptions.
- **SOURCE:** inspected in the pinned IQ-TREE3 checkout identified below.
- **RERUN:** executed during this synthesis using the supplied self-contained script or the explicitly described audit variation.
- **REPORTED:** present in an earlier document but not independently reproduced here.
- **DECISION:** current engineering or scientific protocol choice.
- **GATE:** evidence or implementation work required before the associated capability is accepted.

Source baseline: `63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de`, the upstream `master` retrieved earlier in this session and rechecked locally for this synthesis. Pin this commit for reproduction; do not interpret these findings as a guarantee about future upstream versions. No constrained IQ-TREE extension was compiled or benchmarked in this synthesis.

## 1. Executive decision register

| ID | Issue | Current unified decision | Alternative retained / reopening condition |
|---|---|---|---|
| D01 | Jump coordinates | **Log-ratios first**, with stable row softmax and shared reference destinations | Positive ratios remain a maintained benchmark alternative; switch or expose both if matched compiled tests support doing so |
| D02 | Derivatives | **Scale-aware, bounds-aware finite differences from the first fitted log implementation**, in single-model and linked paths | Existing relative rule is a diagnostic control; analytic gradients are a later option |
| D03 | Bounds | Explicit coordinate domains, seed containment, bound diagnostics, and domain-sensitivity checks | Expand/revise numerical limits or add a boundary-capable method when evidence requires it |
| D04 | T3 | Implement as a serious second IQ-TREE backend using the same likelihood framework | Production exposure remains subject to benchmarks and maintainer preference |
| D05 | Initial matrices | Target-rebuilt empirical exchangeabilities for the primary fitting seed and guide-candidate protocol | Ordinary guide candidates and transported seeds remain sensitivity/multistart alternatives |
| D06 | Integration | Change the coordinate-to-generator map and protect its invariant throughout model lifecycle; reuse native likelihood and nuisance algorithms | A small adapter may orchestrate outer candidate reselection; no assumed automatic support for mixed candidate sets |
| D07 | Target | Explicit immutable target; default estimator can reproduce IQ-TREE's pooled original-training-data procedure with declared smoothing | Other declared estimators and supported small positive targets remain valid; no silent floor imposed on supplied input |
| D08 | Scientific comparisons | Include reversible/unconstrained controls, target refits, and transport comparators in the planned paper | Transport arms are not a prerequisite to mathematical core or first compiled-model validation |
| D09 | Mathematical qualifications | Adopt corrected bounds, conditional nesting statements, and boundary/local-optimization qualifications | Never infer exact global fitting, guaranteed interior MLEs, or a universal coordinate winner |
| D10 | Alternatives beyond charts | Direct balanced flux is an independent reference; EM is deferred | Revisit for persistent direct-optimization failures, exact boundary requirements, or likelihood-gradient infrastructure |

D01 settles what to implement first. It does not require an implementer to relitigate coordinates before starting. D02 is a deliberate resolution of the remaining disagreement with D-review: do not defer the derivative rule until a later fitted-model milestone. D07 resolves the other residual disagreement: default threshold compatibility is not a mathematical mandate to floor every requested target.

## 2. Mathematical and statistical contract

Fix one strictly positive target vector $\pi^\star$, in IQ-TREE order `A R N D C Q E G H I L K M F P S T W Y V`, summing to one. The initial scope is a homogeneous 20-state generator shared across the designated training alignments, with stationary root distribution $\pi^\star$ and scalar site-rate models.

The intended positive family is

$$
\mathcal G_{\pi^\star}=\{Q:q_{ij}>0\ (i\ne j),\ Q\mathbf1=0,\ \pi^\star Q=0,\ -\textstyle\sum_i\pi_i^\star q_{ii}=1\}.
$$

There are 380 off-diagonal entries, 19 independent stationarity equations, and one independent scale equation: dimension 360. Reversible members have dimension 189; the remaining 171 dimensions allow nonreversibility. Reversible matrices belong to the fitted general family and must not be excluded. These are parameter dimensions, not proofs of statistical identifiability.

Define off-diagonal stationary flux $f_{ij}=\pi_i^\star q_{ij}$. Its incoming and outgoing sums match at each state and $\sum_{i\ne j}f_{ij}=1$. With $B$ the directed incidence matrix, $D$ obtained by dropping one row, and $C=[D;\mathbf1^T]$, the reference constraint is $Cf=(0,\ldots,0,1)^T$, $f>0$. The matrix $C$ has rank 20 and is independent of the target. This is the authoritative mathematical feasibility definition, whichever computational chart is used.

Use row-vector distributions and one-based mathematical state indices in this document; C++ array indices are zero-based. Fix a directed-edge ordering before vectorizing flux. A convenient convention is row-major off-diagonals with incidence column $e_j-e_i$ for edge $(i,j)$.

The conditional matrix update maximizes the ordinary summed alignment log-likelihood, with current trees, roots, branch lengths, and rate parameters fixed. Longer alignments contribute through their site-pattern multiplicities; equal weighting of per-gene mean log-likelihoods would change the objective. The broader native training loop updates its permitted nuisances. Specify the nuisance policy rather than presenting all these objectives as one identical conditional operation.

**The central new operation is the map from 360 optimizer variables to valid $Q$.** Given that $Q$, keep IQ-TREE's nonreversible transition and pruning machinery. Both imposing $\pi^\star Q=0$ and setting the root distribution to $\pi^\star$ are required for this baseline. A fixed nonstationary root with unconstrained $Q$ is a different valid model, not an implementation of this one.

The unrestricted interior charts cover the entire strictly positive family. Finite boxes restrict it. An open positive family need not attain its likelihood supremum; exact zero rates occur only in its closure. Neither BFGS nor outer tree search proves global optimality. Do not write “optimizes the full 360-dimensional space” without distinguishing theoretical coverage, actual numerical domain, and optimization success.

**Out of initial scope:** branch-specific targets, arbitrary nonstationary root fitting, profile mixtures with distinct equilibrium vectors, exact sparse-boundary fitting, and replacing IQ-TREE's likelihood engine. They require separate designs and must not enter implicitly through existing options.

## 3. D01: coordinate choice, rationale, and retained positive alternative

### 3.1 One jump-chain construction, two coordinate representations

Use $K$ for the jump matrix and reserve $P(t)=e^{tQ}$ for finite-time transitions. Fix one reference destination $r(i)\ne i$ per row. With $z_{i,r(i)}=0$, define

$$
K_{ii}=0,\qquad K_{ij}=\frac{e^{z_{ij}}}{\sum_{k\ne i}e^{z_{ik}}},\qquad j\ne i.
$$

Implement stable softmax over the 19 off-diagonal destinations. There are 18 free log-ratios in each row, hence 360 variables. Solve for the normalized stationary vector of $K$:

$$
\nu K=\nu,\quad \nu\mathbf1=1;\qquad
q_{ij}=\frac{\nu_i}{\pi_i^\star}K_{ij}\ (i\ne j),\quad q_{ii}=-\frac{\nu_i}{\pi_i^\star}.
$$

$K$ specifies substitution destinations; the exit rates $\nu_i/\pi_i^\star$ specify waiting times. $\nu$ is the jump-event distribution, not the elapsed-time equilibrium. Incoming and outgoing flux are both $\nu_i$, and total flux is one. Thus the construction is target-stationary and normalized.

Conversely, any $Q\in\mathcal G_{\pi^\star}$ yields $K_{ij}=q_{ij}/(-q_{ii})$, $\nu_i=\pi_i^\star(-q_{ii})$, and $z_{ij}=\log(K_{ij}/K_{i,r(i)})$. This proves complete, one-to-one coverage of the positive family. Irreducibility suffices for the stationary solve; do not require an aperiodic chain or an iterative limiting distribution calculation.

The positive version uses $w_{i,r(i)}=c>0$, normalizes each off-diagonal weight row, and stores the other 18 weights. The exact correspondence is

$$
w_{ij}=c e^{z_{ij}},\qquad z_{ij}=\log(w_{ij}/c).
$$

Retain $c=10$ for reproduction of the original positive proposal. It is a numerical scaling convention, not a rate or biological constraint. A different pin or preconditioner is a separately labeled benchmark variant.

### 3.2 Why the recommendation changed

The initial preference for positive ratios weighted compatibility with IQ-TREE's relative finite differences and positive-rate conventions heavily. Those are real advantages, but insufficient evidence about BFGS behavior. The later comparison supplied a concrete mechanism and experiments: identity-initialized BFGS, unequal coordinate scales, bound clamping, and line-search failure can cause premature stopping in positive coordinates.

For the negative log-likelihood $f$, the gradients satisfy $\partial f/\partial z_k=w_k\,\partial f/\partial w_k$. Equal identity metrics in the two coordinates therefore do not produce equivalent steps. Small positive weights can have large absolute derivatives relative to their size; a proposed additive step may hit the floor. Log coordinates make steps multiplicative in the weights and often offer more suitable scaling, but no universal conditioning theorem follows.

The source-confirmed line-search/stopping behavior is described in Section 5; the newly reproduced numerical evidence and its limitations are in Section 10. The evidence favors log-ratios as a starting choice, not as a universally superior method. Mapping a separate parameter vector into rates is already a live IQ-TREE pattern in `ModelLieMarkov`; either new chart needs explicit packing, so positive weights do not eliminate most integration work.

### 3.3 Preserve positive ratios as an actual option

- Maintain compatible build/inverse and packing support sufficient for a compiled positive-ratio benchmark. A public CLI selector is a maintainer decision, not a prerequisite.
- Use identical targets, references, initial physical matrices, data, nuisance snapshots, and exactly mapped boxes for log-versus-positive comparisons.
- Separate the legacy-relative-difference positive variant from any improved-gradient or preconditioned positive variant; do not lump different optimizer configurations together.
- Retain per-start scores and failure diagnostics. One chart need not win every example; the rerun includes an example where positive ratios obtain a higher score.
- Reconsider the default if positive ratios reproducibly offer better reliability or efficiency at comparable validated solution quality, or if log implementation problems persist after derivative/domain corrections.

Neither the old positive preference nor the new log preference licenses deleting the alternative or tuning a comparison unfairly.

### 3.4 Reference convention

Choose reference destinations once from the row maxima of a designated shared target-rebuilt LG seed; break ties deterministically in amino-acid order. This is the default model-definition convention. Record it and use it for every linked object and matched multistart run, even when an alternate start comes from another matrix. Custom-seed-only workflows may choose another declared reference definition, but comparisons must preserve or explicitly map it.

Changing references is mathematically harmless without bounds, but generally changes the physical domain of a fixed rectangular box. Do not reselect them independently per partition, per likelihood evaluation, or per start. Any later explicit reference change must reencode the accepted matrix, transform/reconsider domain constraints, and restart optimizer curvature state; it is not a silent numerical tweak.

## 4. D02–D03: derivatives, bounds, and numerical acceptance

### 4.1 Derivative policy is part of the first fitted implementation

The inspected rule is a forward difference with $h_k=10^{-4}|x_k|$, using $10^{-4}$ only when the computed step is exactly zero. For log-ratios, values arbitrarily close to zero are ordinary equal-weight neighborhoods. The step can then be too small to resolve a likelihood difference. Exactly zero triggers the fallback and is not itself the problematic point. The same issue applies to signed T3 tilts.

Implement an opt-in coordinate-aware helper from the first log-coordinate fitting milestone, with initial design $h_k=\eta\max(s_k,|x_k|)$ and positive characteristic scales. For dimensionless logits/tilts, $s_k=1$, $\eta=10^{-4}$ is a reproducible initial setting, not a universally validated constant. Use representable increments, bounds-aware probing, and a validated one-sided or central stencil. Central differences are a validation reference; production use depends on accuracy/cost tests. At active bounds, use feasible one-sided checks and appropriate constrained optimality diagnostics.

The helper must be reached from both `ModelMarkov` and `PartitionModel` optimization. The latter runs the optimizer on the partition object and evaluates the sum of linked likelihoods; overriding only the model's `derivativeFunk` does not cover it. Preserve legacy behavior by default. Restore the unperturbed coordinates, generator, decompositions, and likelihood caches after probes, including errors and final termination. Reevaluate the returned accepted state before reporting its score.

A constructor/build-only prototype need not have an optimizer. A fitted log-model prototype must not declare the old derivative rule adequate solely because a few ordinary runs did not approach zero. Future analytic gradients require the derivative through the stationary/balancing solve and the full tree likelihood; a builder derivative alone is not a likelihood gradient.

### 4.2 Numerical domains must be explicit

For the historical positive bounds $w\in[10^{-4},100]$, $c=10$, the exactly corresponding log bounds are

$$
z\in[\log(10^{-5}),\log(10)]\approx[-11.5129255,2.3025851].
$$

Use these as the **historical matched-benchmark domain**, not a claim of an optimal production range. For initial implementation runs, start there only if all declared starting matrices are representable with margin; otherwise expand the domain explicitly before fitting. Final scientific runs require a recorded domain and expansion-sensitivity protocol. No arbitrary universal final box is selected here.

A positive row pin of 10 with these bounds restricts ratios to the reference to $[10^{-5},10]$. This is not the old box on raw directed rates. A raw balanced triangle with forward rates 100, reverse rates $10^{-4}$, and other off-diagonals 1 needs a weight $10^{-5}$ when its row maximum is pinned to 10: below the proposed $10^{-4}$ floor. Reusing numerical constants is not domain equivalence.

For every run, record limits, reference metadata, minima of rates/fluxes, and bound activity. Use one mapped ratio/log measure for cross-chart proximity comparisons; a fixed absolute tolerance in unrelated coordinates does not identify the same near-boundary set. Exact clipping events and distance to a bound are distinct diagnostics. A convenient shared measure is distance in log-ratio space for both jump implementations.

Expand bounds when active or when sensitivity tests show material changes; checkpoint and retain the incumbent. Absence of active bounds is not proof that the restricted solution equals the unrestricted optimum. Never silently clip an imported feasible $Q$ to fit a chart box. Rebuilding from clipped coordinates changes the matrix even if it remains target-stationary.

### 4.3 Builder validation and failure policy

The stationary solve on $K-I$ can reuse IQ-TREE's augmented QR approach, with explicit row/column layout tests. Keep $\nu$ separate from immutable $\pi^\star$. Check finite outputs, positivity, stationarity, mean rate, and residuals; the existing assertion about the sum alone is insufficient.

Complete numerical diagonals from the installed off-diagonal row sums, and compare them with the theoretical $-\nu_i/\pi_i^\star$. This preserves the generator row-sum convention while making floating-point inconsistencies visible. Keep any roundoff normalization distinct from clipping, flooring, or changing the target.

For ordinary double-precision cases, $\|\pi^\star Q\|_\infty$ and $|{-\pi^\star\operatorname{diag}Q}-1|$ around $10^{-10}$ or smaller are provisional feasibility goals. Validate tolerances against conditioning and likelihood impact. Report absolute and scaled row residuals; a large norm of $Q$ must not hide flux imbalance. Near reducibility or tiny targets require stress checks, not blanket claims that all positive inputs are numerically safe.

Underflow to exact zero, a negative stationary component, or failed balancing must not publish a partially updated model. Reject the proposal with a diagnosable status and preserve accepted state. Repeated numerical failures mean the run is not validated; they must not masquerade as a flat optimum. Do not repair a general generator by independently clipping its physical entries, which can violate stationarity.

## 5. Optimizer reuse and stopping: retain algorithms, verify results

**SOURCE:** `dfpmin` initializes its inverse Hessian to identity in the usual path. `lnsrch` forms a trial, clamps it, and compares its value with sufficient decrease predicted from the original direction. When the step falls below its minimum, it restores the old parameter vector and sets `check=1`. `dfpmin` does not use this flag before testing parameter displacement; a zero displacement immediately satisfies its small-step test. This is a failure-reporting/convergence-certification risk, not proof that every such stop is away from an optimum.

There is an additional source/port difference: C++ leaves the objective output at the last tested trial when restoring the old vector; the supplied Python port returns the old objective. The targeted audit in Section 10 found no change to fitted matrices in its four examples, but the port is not literally identical on this path. Native callers may subsequently resynchronize/recompute; do not assume the transient discrepancy necessarily survives into final IQ-TREE output.

Keep BFGS as the initial search algorithm. Add opt-in diagnostics/status propagation for the new fitting paths so failed line searches, small displacement, gradient tolerance, and iteration exhaustion are distinguishable. Preserve legacy numerical trajectories unless maintainers approve a broader repair. If source-level status propagation is unavailable in an early prototype, its stop must be labeled unclassified until external convergence checks pass; “converged” is not justified by return alone.

At termination:

1. Restore/rebuild the actual returned coordinate vector and reevaluate its full objective.
2. Check all model invariants and retain the best valid incumbent.
3. Validate a first-order diagnostic in a common declared representation; native `gtol` values in different coordinates are not equivalent accuracy standards.
4. Record domain restrictions, failed searches/trials, and stopping reason.
5. If checks fail, return a valid incumbent with an explicit incomplete/failure status; do not claim a successful optimum.

For matched jump comparisons, transform gradients by $g_z=w\odot g_w$ and evaluate a common normalized diagnostic, supplemented by feasible directional checks. For comparisons with T3/flux, use the same reference coordinates or feasible flux directions for post-fit assessment. Near active constraints, assess constrained/KKT-style conditions for the actual declared domain rather than demanding zero unconstrained gradient. These are common diagnostic conventions, not coordinate-free universal norms.

Keep targeted restarts/multistart available. Do not automatically enable the legacy bound-triggered randomization of all coordinates at dimension 360. A failed search can motivate a controlled restart or a future improved line search; any algorithmic recovery must be separately labeled and tested. Log coordinates reduce some scaling risks but do not fix the underlying stopping logic.

## 6. D04–D05: T3, empirical candidates, and starts

### 6.1 T3 is a second fitting backend, not just a transport utility

Implement T3 after the log-jump core is validated, using the same IQ-TREE likelihood engine and derivative-policy interface. A Python reference is useful for correctness but cannot by itself provide a fair compiled runtime comparison. T3's iterative inner solve is a plausible cost; its practical importance must be measured, especially when repeated across partitions and finite-difference probes.

Adopt P-log's explicit gauge convention for the unified implementation: symmetric log strengths $a_{ij}=a_{ji}$, one pin $a_{19,20}=0$, antisymmetric tilts $h_{ij}=-h_{ji}$, and $h_{i,20}=0$ for $i<20$. This leaves 189 plus 171 free coordinates. P-positive's star-at-state-1 convention describes the same unrestricted family but must be converted explicitly; never import gauge-fixed coordinates as if the gauges agreed.

With $g_{20}=0$, form

$$
N_{ij}=\pi_i^\star\pi_j^\star e^{a_{ij}+h_{ij}},\qquad
\Phi(g)=\sum_{i\ne j}N_{ij}e^{g_j-g_i}.
$$

The unique minimizer balances the positive flux $x_{ij}=N_{ij}e^{g_j-g_i}$. Normalize $f=x/\sum x$, then divide row $i$ by $\pi_i^\star$ and complete the diagonal. The 19 balancing potentials are derived quantities, not extra model parameters. The reduced Hessian is a positive-definite weighted graph Laplacian. This convex inner problem does not make the outer observed-data likelihood concave.

At $h=0$, balancing gives $g=0$ and the ordinary reversible construction. Every positive normalized target-stationary $Q$ is represented after gauge fixing. Use damped Newton, stable exponential scaling, residual checks, and deterministic accuracy. Warm starts may improve speed but must not create evaluation-order-dependent likelihood noise large enough to corrupt derivatives. Keep mutable solver state per model object.

T3 and jump backends optimize the same family before numerical restrictions. Differences in attained scores reflect optimization, bounds, or implementation, not extra biological expressiveness. T3's cross-target interpretation is useful for a transport comparator but is not required to solve this fixed-target project.

### 6.2 Separate the two meanings of “initial matrices”

**Fitting seed:** for a symmetric empirical exchangeability matrix $R$, build $\widetilde q_{ij}=R_{ij}\pi_j^\star$, derive diagonals, and normalize to unit mean rate. Use target-rebuilt LG as the default first start. Import a fitted reversible target matrix as an important additional start and nesting incumbent. WAG/JTT-derived starts and feasible perturbations remain available within a declared budget.

This target-aware seed is a chosen initialization, not the only mathematically valid seed. Row-normalizing LG at LG's original frequencies and rebuilding at the target produces a valid flux-preserving transport, but generally changes its exchangeabilities. Do not call it LG with unchanged exchangeabilities at the target. In the inspected `NONREV` initialization, conversion precedes reading the literal target frequency suffix; the new class must explicitly build its seed in the correct order or through a separate builder. Do not assume inherited `NONREV+F{...}` achieves the desired initialization.

**Guide candidates:** prefer target-rebuilt LG/WAG/JTT for the primary guide-tree/model-selection protocol. This keeps the frequency convention consistent. Ordinary empirical candidates are valid guides and remain a sensitivity arm; they do not invalidate later constrained fitting. Neither protocol is proven to give better trees or a better final local optimum.

Explicitly suppress unintended empirical/optimized frequency variants in the target-consistent candidate set. Existing reversible target construction is available; the complete candidate-generator syntax and custom-model route still require runtime verification. D-review's “unverified; either works” is not evidence that both complete command workflows work.

Changing guide candidates can change topologies, root choices, and rate-family selection. Hold all of those fixed for a coordinate-only benchmark. Evaluate guide-policy sensitivity separately rather than attributing its effects to log versus positive coordinates.

### 6.3 Other starts and later rounds

For a positive normalized nonreversible $Q_0$ with equilibrium $\mu$, flux-preserving transfer uses $f^0_{ij}=\mu_iq^0_{ij}$ and $q^{\rm start}_{ij}=f^0_{ij}/\pi_i^\star$. T3 transport is a separate retained alternative. Neither defines the final constrained estimator if all 360 variables are subsequently optimized.

Reencode the previous round's feasible $Q_{\rm best}$ at the same target exactly. No biological transport is needed when the target is unchanged. Preserve shared references even for alternate starts. If importing a sparse matrix into an interior method, use a declared feasible flux mixture with a positive reference (or a symmetric exchangeability mixture for a reversible start); do not independently clip rates.

For nesting checks, retain the entire matching nuisance snapshot along with the matrix. Preserving a matrix alone does not preserve its likelihood after trees, roots, or rate parameters have changed.

## 7. D06: integration contract and native workflow

### 7.1 Authoritative state and model identity

Use a dedicated fixed-equilibrium nonreversible model identity. Names such as `NQC` in the plans are proposed syntax, not currently available commands. The maintainers choose public names/options. Rooted/nonreversible capability belongs to the model class even when its initial matrix happens to be reversible.

The public target interface must distinguish equilibrium enforcement from a root-frequency-only request. Test its parsing through the joint-model path, including numeric precision and scientific notation; do not assume generic suffix tokenization safely preserves every valid numeric string. A target file or canonical serialized vector are reasonable interface options, subject to maintainer choice.

For fitted models, store immutable normalized target and state order; chart type/version; fixed reference/gauge definitions; authoritative 360 coordinates; domain/derivative/solver policy; and fixed/fitted status. Physical `rates[]`, normalized `rate_matrix`, and decompositions are derived. For a fixed exported model, the validated full normalized $Q$ and its target are authoritative.

Use separate coordinate storage and a pure builder rather than ambiguously using `rates[]` sometimes for auxiliary jump weights and sometimes for physical intensities. The actual pinned quantity is a jump weight, not a physical rate: all 380 physical rates may change when 360 coordinates move. This supersedes a literal implementation of P-positive's auxiliary-weights-in-`rates[]` design and resolves ambiguities in D-review's illustrative pseudocode.

`setVariables` packs model state into the optimizer array; `getVariables` reads it back. Override both deliberately. `getNDim()` is 360 when fitted and zero when fixed. An externally fixed target adds zero frequency dimensions; linked matrix parameters are counted once across partitions. An imported fixed candidate contributes zero matrix fitting dimensions in its downstream evaluation.

### 7.2 Generator installation and decomposition

Prefer a narrow shared numerical-decomposition helper/hook that accepts an already built generator over duplication of the nonreversible eigen/fallback machinery. Keep derived rate storage coherent. The existing base routine rebuilds from `rates[]`, so installing only `rate_matrix` is unsafe. A harmless roundoff normalization is acceptable only if it preserves the fixed-target invariant and declared branch-length units.

`FREQ_USER_DEFINED` bypasses the linked constructor's ordinary frequency pooling, but does not alone enforce the model. Ensure `optimize_from_given_params`, setters, seed import, restart, and inherited routines cannot replace the target or release stationarity. A diagnostic equilibrium solve on final $Q$ must not overwrite the target. General root-frequency-only semantics remain a separate legacy path.

Audit constructor order and allocations before invoking derived methods. Reject unsupported model combinations clearly. Retain the existing nonreversible kernel throughout fitting rather than changing tree representation at reversible parameter values.

### 7.3 Linked fitting, concurrency, and state restoration

Before optimization, check equality/compatibility of state order, target values, chart/gauge, references, numerical domain, and derivative policy for all linked objects. Name matching alone is insufficient. The existing partition path sends the same vector to matching objects and sums their objectives; different coordinate definitions would silently fit different matrices.

Keep matrices fixed during per-partition nuisance updates and unfix only for the shared update. Preserve existing deterministic objective reduction where applicable. Initially allow each partition to build its matrix independently from identical inputs; shared prepared-matrix caching is optional later optimization and must be immutable/thread-safe.

Check final state consistency after every probe sequence, failure, accepted update, and restart. No shared mutable T3 warm-start object across threads. Test partition ordering and different thread counts with appropriate numerical tolerances.

### 7.4 Native alternation versus outer candidate orchestration

Deliver in layers:

1. **Conditional matrix fitting:** freeze all nuisance quantities to isolate Step 3A and compare numerical methods.
2. **Native continuous alternation:** reuse the existing partition scheduling for permitted branch, site-rate, and partition-rate updates. Both separate-tree and edge-proportional arrangements need tests when supported. Do not force them into one tree-sharing interpretation.
3. **Outer candidate reselection:** add explicit orchestration if the full nQMaker-style candidate/tree loop is part of the analysis. A thin driver using supported calls is acceptable; native integration can follow.

The current continuous scheduling is broader than a strict reading in which the selected site-rate model's numerical parameters are all frozen. Rate-family selection and continuous parameter updates are different operations. Document which move. Keep native algorithms unless a justified protocol requires otherwise.

The inspected ModelFinder code rejects recognized mixed reversible/nonreversible candidate sets. Do not assume adding a learned matrix to the old list works in one ordinary invocation, remove the guard without a design, or evade it by renaming a model. Use separate compatible candidate evaluations with comparable criteria or a tested common runner. Later learned matrices are fixed candidates, not per-gene 360-parameter refits.

Preserve the previous accepted common-model solution including nuisances. Improved per-gene candidate scores do not establish improvement in the shared-$Q$ objective. Compare and accept on the actual declared shared objective; if the nuisance family changes, state how comparison/incumbent eligibility is handled. Matrix correlation can be retained as an outer diagnostic, but is not sufficient evidence of likelihood convergence. Compare the previous and new matrices, not an overwritten variable against itself.

### 7.5 Roots, reporting, persistence, and regression

Specify stationary root frequencies, input rooting/outgroup conventions, root-edge search, and root split/branch optimization separately. Fixed topology alone is not a complete root policy. For coordinate benchmarks, use identical fully specified rooted trees and branch lengths. For ordinary fitting, preserve native rooting algorithms under the chosen flags.

Checkpoint target, coordinate values, chart version, references/gauges, fixed/fitted status, bounds and tolerances, accepted score, and sufficient nuisance state. Reject incompatible restarts. Export high-precision normalized $Q$, target, and state order; keep fitted-model metadata separately from a rounded human report. Validate export/reimport likelihood equality with fixed nuisances. Test actual simulation/export paths required by the study.

Regression controls include ordinary `NONREV`, its existing literal-frequency behavior, `GTR20`, `UNREST`, and Lie-Markov models. New derivative/status hooks must default to existing behavior. Require numerical equivalence under the same environment; do not promise universal bit identity across compilers or linear-algebra libraries.

## 8. D07: target provenance, flooring, and supported inputs

The fitter consumes one explicit target and must not reestimate it during matrix optimization. For the intended cleaning study, a reasonable default is pooled original-training-data composition, computed once and held constant across cleaning treatments. Reproduce IQ-TREE's actual count/ambiguity procedure rather than assuming it equals an arithmetic mean of per-gene frequency vectors. Record taxon/site/gene weighting, missing-data conventions, the eight-iteration ambiguity procedure if used, smoothing, and the exact normalized numbers.

A kept preprocessing implementation may reproduce IQ-TREE's estimator, or use a validated native extraction route. Compare against the corresponding native calculation under the same partition/missing-data settings. Do not claim exact equivalence without that cross-check.

Three distinct policies must remain separate:

| Policy | Unified treatment |
|---|---|
| Estimator floor on a newly computed target | Allowed as an explicit scientific/procedural choice, including the native default $10^{-4}$ adjustment |
| Acceptance threshold in the fitter | Must be compatible with the exact supplied target and numerical backend; a limitation may be reported rather than silently modifying input |
| Mathematical positivity requirement | Every target entry must be strictly positive for this initial 20-state interior design; no fixed $10^{-4}$ mathematical lower bound exists |

`ModelMarkov::targetFunk` rejects positive frequencies below configured `min_state_freq`, default $10^{-4}$. **SOURCE:** `--min-freq` makes that threshold configurable. Therefore “the floor is not optional” is too strong. Merely lowering it does not establish numerical support: the decomposition also uses `ZERO_FREQ = 10^{-10}` to select states in the inspected path. A tiny positive target must not silently become a reduced-state computation. Require backend-wide checks or reject unsupported values clearly.

If input is defined as weights, normalize once, explicitly, and record the resulting target. If it is defined as an already normalized exact target, use documented sum tolerance and report any normalization; do not hide semantic changes behind a warning. Reject nonfinite entries, wrong length/order, or zeros for this model. No automatic flooring or composition substitution inside the builder.

When the target is estimated from the training data, inference is conditional on a plug-in estimate. Counting 360 matrix coordinates remains correct for the conditional optimization, but a BIC calculation across fixed/free-target procedures is not justified merely by adding 19 to a parameter count. Resampling that addresses target uncertainty must recompute the target as well as the fitted matrix.

## 9. D08–D10: scientific controls, corrections, and alternatives

### 9.1 Keep numerical and scientific comparisons distinct

Plan these scientific arms, as applicable to the training/evaluation design:

| Arm | Purpose |
|---|---|
| Original-data reversible and unrestricted general fits | Reference working models before cleaning |
| Cleaned-data reversible and unrestricted general fits | Baselines with composition fitted/derived rather than constrained to the original target |
| Cleaned reversible fit followed by exchangeability-preserving frequency replacement | Historical replacement procedure |
| Reversible refit at the target | Matched reversible null/control and feasible incumbent |
| General refit at the target | Project 2 estimate |
| Cleaned unrestricted general fit transported by fixed flux and by T3, without refitting | Separate transport choices from likelihood reestimation |

Include the transport arms in the paper plan, but do not make them a prerequisite to the first model implementation. They are not interchangeable with comparing T3 fitting to jump fitting. They need no new matrix refit, but their common downstream evaluation still has computational cost. Prespecify confirmatory claims and label exploratory work honestly; a document should not retroactively declare an analysis preregistered.

Evaluate frozen learned matrices on the same held-out observations with a declared nuisance policy. Raw training scores on differently cleaned alignments are not directly comparable likelihoods on identical data. Use training/validation splits for method tuning and reserve held-out test genes for final evaluation. Record likelihood, relevant tree and branch-length outcomes, failure rates, and uncertainty at the appropriate resampling unit.

Since $\pi^\star e^{tQ}=\pi^\star$, stationary root initialization fixes all tip marginals. Nonreversibility cannot repair a genuine marginal composition mismatch, although it may improve joint-pattern fit under misspecification. A wrong target can produce apparent circulation; it does not inevitably do so. A finite-sample frequency difference alone is not proof of population mismatch. Simulate reversible and nonreversible processes with correct and perturbed targets, and include outcome-independent masks plus relevant content-dependent cleaning. A reversible-only refit cannot diagnose spurious circulation because it disallows it.

Name circulation measures precisely: antisymmetric flux, raw pairwise tilt, and cycle affinity differ. Fixed-flux transport preserves normalized flux; T3 transport preserves cycle affinities, not every raw tilt or flux asymmetry. A transport comparator is not by itself a calibrated test for artifacts.

### 9.2 Mathematical corrections that override conflicting prose

1. **T3 triangle bound exists.** In a star-at-1 gauge,
   $$h'_{ij}=\tfrac12\log\frac{q_{1i}q_{ij}q_{j1}}{q_{i1}q_{ji}q_{1j}}.$$
   If all raw off-diagonals are in $[l,u]$, then $|h'_{ij}|\le\tfrac32\log(u/l)$. For $[10^{-4},100]$, this is $20.7232658\ldots$. The triangle example attains it before/after common normalization in the ratio expression. A bound of 12 does not contain that raw box. “No bound can be derived” is false, including its recurrence in D-review Section 3. A triangle bound alone does not establish containment of a complete T3 strength/tilt box, nor does it bound the unrestricted positive family.
2. **Sharper physical-rate bound.** Normalized stationary outgoing flux from a state equals its incoming flux; the two disjoint edge sets together have total at most one. Therefore $f_{ij}\le1/2$ and $q_{ij}\le1/(2\pi_i^\star)$. Equality belongs to suitable boundary cases; the bound is strict for the fully positive 20-state family. The older $1/\pi_i^\star$ bound is valid but loose, not a false inequality.
3. **Nesting is mathematical; independently attained scores need not be ordered.** Global suprema over nested physical/nuisance domains are ordered. Local fits can reverse the ordering. Import the smaller-model solution and its nuisances into the larger model, verify representability under actual bounds, and retain it as an incumbent for an operational nondecrease check. Merely finding two interior estimates does not prove domains are nested. Active bounds do not by themselves prove nonnesting either; inspect the domains.
4. **Interior truth does not imply interior sample MLE.** Neither many sites nor generating rates away from the floor establishes interiority. Do not label E6 floor hits spurious “by construction,” or E5's optimum boundary from its sample size alone. Conversely a floor hit does not prove a true boundary optimum.
5. **Root singularities qualify LRTs.** Under reversibility, the split of an unrooted edge into two root-adjacent lengths is unidentifiable; fixing the root edge alone is insufficient. Parameter differences 171 or 19 do not automatically yield chi-square reference distributions. Match root/split policies and verify regularity, or bootstrap the actual procedure. A derivative demonstrating local root sensitivity is not a global identifiability theorem.
6. **Exact constraint versus generator validity.** Replacing frequencies in asymmetric factors can still yield positive off-diagonals and zero row sums; the generally lost property is target stationarity. Fixing root frequencies alone is likewise not the required equilibrium constraint.
7. **No universal optimum or uncertainty claim.** Chart completeness, convex flux feasibility, and convex T3 balancing do not make the observed likelihood globally concave. Agreement across charts is evidence, not proof. The internal BFGS inverse Hessian is not automatically a covariance after nuisance fitting, target estimation, and tree selection.

### 9.3 Retained but deferred methods

Direct balanced-flux fitting with $f=f_0+Zy$ for a basis $Z$ of $\ker C$ is complete on its coupled-positivity domain. Keep it as an independent small-problem oracle and potential boundary-capable backend. It cannot be made safe by independent box clamping of $y$; use an appropriate feasible constrained solver. Exact zero-rate fitting is a separate capability, not something finite logits provide.

Constrained substitution-history EM remains mathematically useful but is deferred. If implemented for fixed branch lengths and unit-rate normalization, prefer the normalized M step with 19 balance multipliers plus one normalization multiplier. The 19-multiplier unnormalized-cone variant followed by normalization and compensating branch rescaling is a different conditional update, requiring allowed branch rescaling within bounds. Normalizing $Q$ alone does not preserve transitions. Expected-history machinery may later support analytic gradients, but it is a substantial addition to the inspected direct-fitting path.

Penalty or projection methods do not replace exact constrained likelihood estimation. They may initialize or define separately labeled transport/comparison procedures. Keep novel cross-target biological structure choices outside the fixed-target implementation unless separately requested.

## 10. Audit of the supplied numerical evidence

### 10.1 What is available and what was executed

The supplied `chart_fd_bfgs_compare.py` is self-contained with NumPy/SciPy and contains E1 gradient checks and E2 one-call BFGS comparisons. It was run unchanged for this synthesis with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python chart_fd_bfgs_compare.py \
  --n 8 --ntips 8 --nsites 1000 --seed 20260922 --nrep 4
```

Environment: NumPy 2.3.5, SciPy 1.17.0. The source report used different versions. The main program first consumes the RNG in its 20-state E1 construction, then makes E2's four 8-state problems. Reproduce that sequence; calling `experiment_bfgs` immediately with a fresh identically seeded RNG makes different problems. BLAS threading was restricted to make this audit efficient; no runtime comparison against the report's different environment is claimed.

`e6_n20_big.py` imports `IQOptLog` from `e4_nonrev`, which was not supplied. Thus E6 is not reproducible **as delivered**. No replacement for that missing module was silently invented. E3/E4/E5 experiment drivers and their original output logs were also not supplied. Their numerical results remain REPORTED, not RERUN. Obtaining them is useful follow-up for benchmarking, not a blocker for this decision record or the first implementation.

### 10.2 E1 rerun: derivative problem reproduced, universal threshold not established

The E1 data had 20 states, 8 tips, 1000 sites, and 968 unique patterns. Median relative discrepancies against the script's selected central difference were approximately $8.02\times10^{-5}$ for positive coordinates and $1.44\times10^{-4}$ for log coordinates; maxima were approximately 0.00688 and 0.01969. The controlled one-logit sweep gave:

| Logit value | Legacy log-coordinate derivative relative discrepancy | Positive-coordinate discrepancy at the same physical matrix |
|---|---:|---:|
| $10^{-2}$ | $3.24\times10^{-6}$ | $7.51\times10^{-5}$ |
| $10^{-4}$ | $1.19\times10^{-4}$ | $7.72\times10^{-5}$ |
| $10^{-6}$ | $9.55\times10^{-2}$ | $7.72\times10^{-5}$ |
| $10^{-8}$ | $6.39$ | $7.71\times10^{-5}$ |
| $0$ exactly | $1.27\times10^{-4}$ | $7.71\times10^{-5}$ |

The deterioration near zero supports D02. The exact-zero recovery reflects the special fallback. This is not a proof that only $|z|<10^{-5}$ is dangerous or that the selected central reference is exact.

**Reference-quality limitation:** `central_ref` returns all requested values, but `experiment_fd` uses the second of two steps without testing agreement; its controlled sweep uses only one central step. Despite comments mentioning a plateau, the script does not enforce a step-size plateau criterion. Treat these as measured discrepancies against central-difference estimates. Production validation must compare multiple steps and assess absolute as well as relative errors, particularly for near-zero true derivatives.

### 10.3 E2 rerun: favorable but mixed evidence for log-ratios

All comparisons use the same physical seed and exactly transformed boxes for each instance. Larger log-likelihood is better.

| Instance | Positive, legacy difference: logL / stop / evaluations | Log, legacy difference: logL / stop / evaluations | Log, scaled difference: logL / stop / evaluations |
|---|---|---|---|
| 0 | -11343.902245 / TOLX / 278 | -11303.734515 / gtol / 975 | -11303.734377 / gtol / 975 |
| 1 | -12419.902818 / TOLX / 700 | -12403.624673 / gtol / 822 | -12403.624636 / gtol / 822 |
| 2 | -14544.987015 / gtol / 449 | -14543.502241 / gtol / 666 | -14543.503265 / gtol / 666 |
| 3 | -12552.097580 / TOLX / 1426 | -12552.599324 / gtol / 770 | -12552.599299 / gtol / 770 |

The rerun broadly reproduces the reported endpoint pattern, not every count: for example D-review lists 641 positive evaluations for instance 1, while this rerun has 700. Small numerical/library differences can alter a line-search trajectory; bitwise reproduction is not claimed.

Log coordinates give higher scores in instances 0–2; positive coordinates give a higher score in instance 3. One-call evaluation counts do not uniformly favor logs, because some positive runs terminate early. The native stop labels are not equivalent common-accuracy certifications. These examples are eight-state conditional fits, not compiled 20-state nQMaker training runs.

### 10.4 Targeted line-search audit and sensitivity check

A temporary instrumented subclass counted the supplied port's returned failure flags without changing its search. A second variant changed exactly the failure return from `return xold.copy(), fold, 1` to `return xold.copy(), f, 1`, matching the inspected C++ objective-output behavior. Both variants used the same E2 problems. To reconstruct their RNG sequence without redoing E1 derivatives, call `make_problem(20,8,1000,rng)` once with seed 20260922, then create the four 8-state problems in order.

For positive/legacy and log/scaled runs respectively, failed-search counts were `[1,1,0,1]` and `[0,0,0,0]`. Thus the three positive TOLX stops in this rerun coincided with failure returns, while the positive gtol stop did not. The failure-return change left all returned matrices, endpoint scores recomputed from those matrices, and iteration/evaluation counts unchanged in these examples. The largest returned-versus-recomputed objective discrepancy in the source-return variant was about $1.7\times10^{-8}$.

As a limited likelihood-backend check, eigen-based transitions at the fitted endpoints were compared with SciPy `expm` for all branch lengths in these examples. Maximum absolute entry disagreement was below $1.1\times10^{-14}$. This checks these endpoints only, not every optimization trial or difficult nonnormal/defective matrix.

These results strengthen the plausibility of the reported mechanism without establishing differential equivalence to compiled IQ-TREE. No compiled optimizer or full phylogenetic likelihood differential test was run.

### 10.5 Limitations and corrections to the evidence narrative

- The port uses NumPy least squares rather than IQ-TREE's Eigen column-pivoted QR. They solve the same mathematical system in well-conditioned cases, but are not numerically identical implementations.
- The toy likelihood uses NumPy eigenvectors and inversion without IQ-TREE's full fallback, scaling, cache, ambiguity, and rate-category machinery. It is suitable for these small ordinary examples, not a production likelihood substitute.
- The simulated reversible seed is constructed from the symmetric part of the generating flux. It is not empirical LG at an arbitrary target. This is a useful matched seed but limits ecological realism and possible generalization.
- Branch lengths, roots, and topology remain fixed; there is no full native nuisance alternation in E2. Repeated BFGS calls in a separate script would still not reproduce that alternation by themselves.
- The base supplied optimizer does not itself report failure counts; the extra instrumentation is described above. E6's original instrumentation remains unavailable.
- The script's near-bound flags use coordinate-dependent absolute/relative tolerances, so their counts are not identical physical proximity tests across charts. Distinguish exact clipping from loosely defined near-bound counts.
- Bounds can be probed outside by the supplied forward-difference routine; a production bounds-aware implementation is a deliberate improvement, not a literal reproduction of these runs.
- The statements that E6's optimum is interior and its floor hits are spurious by construction are unsupported. The analogous E5 boundary-optimum claim is also unsupported without a solution/optimality analysis.
- Fewer line-search failures or a `gtol` stop in these examples is evidence about these runs, not a guarantee that log coordinates avoid the source vulnerability generally.

**Evidence-weighted conclusion:** starting with log-ratios is justified by mathematical equivalence, a credible source mechanism, and reproducible conditional-fit evidence. Comparative superiority in compiled 20-state training remains a GATE. The unavailable E6 result is not essential to the current decision and must not be promoted to verified evidence in a paper or proposal.

## 11. Development sequence and acceptance gates

Deliver incrementally; do not expand every milestone into a new optimizer project. Preserve decision IDs and record deviations with evidence.

| Gate | Deliverable | Acceptance conditions |
|---|---|---|
| G0: protocol/configuration | Explicit target and model definition; numerical/scientific run manifest | State order, references/gauges, chart, target provenance, domain, derivative policy, root/tree/rate policy, seed set, and source commit recorded |
| G1: mathematical core | Jump log and positive build/inverse; T3 reference; flux conversions | Tests at 3, 4, 20 states; positive random and balanced-cycle examples; round trips; reversible reduction; target stationarity; normalization; adversarial conditioning diagnostics |
| G2: first compiled conditional model | Log coordinates with initial derivative helper and status diagnostics | Same-$Q$ transitions/likelihood versus independent reference; gradient step-size checks; no target overwrite; rejected invalid targets; accepted-state restoration; source-versus-port toy case comparison |
| G3: linked and lifecycle support | Shared model fitting, persistence, fixed export/import | Two distinct-composition partitions retain one target and one physical $Q$; 360 counted once; identical metadata; fixed outside shared update; restart/export likelihood agreement; thread/order and legacy regression tests |
| G4: competing implementations | Positive-ratio compiled benchmark; compiled T3 backend; small direct-flux oracle | Identical physical starts/objectives; matched jump boxes; explicit T3-domain handling; common post-fit diagnostics; failed-search/status reports; evaluation/runtime/multistart comparison |
| G5: native and outer workflow | Native continuous alternation; target guide candidates; fixed learned candidates; optional outer driver | Intended separate-tree/edge-proportional paths tested; candidate frequency controls; safe mixed-family orchestration; incumbent protection; root policy; simulation/export support |
| G6: scientific validation/default selection | 20-state simulations and real training/evaluation runs | Domain sensitivity; multiple starts; mismatch/cleaning controls; held-out evaluation; qualified inference; evidence-based production default decision |

G2 must include scale-aware derivative support for both optimizer entry points even if its first executed test is a single alignment. It can precede full G3 linked validation. T3's standalone mathematical reference may exist at G1, but a serious compiled comparison is G4; do not treat a Python-only oracle as satisfying the user's T3 comparison goal.

At least two complementary starts should be checked in final matrix-fitting validation, including a fitted reversible target incumbent where applicable. Larger multistart budgets and numerical thresholds should be justified by observed convergence, then frozen for the main comparison. Controlled domains should contain all compared starts; no clipping to conceal an incompatibility.

For T3 versus jump benchmarks, independently chosen rectangular boxes are not automatically equivalent. Prefer a common declared physical domain when feasible. Otherwise explicitly report restrictions, cross-encode endpoints/starts, expand boxes, and test sensitivity before attributing differences to coordinate geometry. Endpoint representability alone does not prove domain equality or equal global optima.

**Required run record:** source/build and dependency versions; data identifier/split; seed; target and provenance; chart/reference/gauge; numerical limits; nuisance snapshot/policy; initial/final and best accepted likelihood; reevaluated final likelihood; feasibility residuals; derivative/optimality diagnostics; stopping reason; failed trials/searches; bound/clipping diagnostics; evaluation counts including diagnostic evaluations; runtime/hardware/thread settings; and restart/export status where relevant.

## 12. Crosswalk: what to inherit, override, or retain as an alternative

| Original location | Treatment under this synthesis |
|---|---|
| P-positive §0, §3.1.3, §13b: positive first / unchanged optimizer | Override default with D01 log first and D02 derivatives now; preserve the positive chart and rationale as the required comparison alternative |
| P-positive §3.1.4 and P-log §5: chart bounds | Retain mathematical caveats; replace any blanket unchanged-box implication with D03 and mapped-domain checks |
| P-positive §3.2.5, §4.6 and P-log §6.2: derivative path | Retain both source findings; adopt P-log's targeted interface, including the linked path |
| P-positive §3.2.6: T3 external only | Override with D04 compiled comparison; retain Python oracle |
| P-log §3 and §4: chart definitions | Retain, with shared jump references and the explicitly adopted T3 gauge; import other gauges only by conversion |
| P-positive §3.1.5, P-log §2.3–2.4: seed distinction | Retain; target-aware reversible seed is primary, transported starts remain alternatives |
| P-positive §5 versus P-log §7.3/§8.5: guide candidates | Adopt target-rebuilt primary protocol; retain ordinary guides as sensitivity analysis |
| P-positive §5/§6: “only 3a,” rate storage, checkpoint unchanged | Preserve core algorithm framing but override incomplete lifecycle/storage implications with Section 7 |
| P-log §7.4/§8: workflow, invariant, persistence | Adopt with the scoped gates and maintainers' public-interface choices |
| P-positive §7 and P-log §1.4/§8.5: target estimator/floor | Combine reproducible pooled-estimator option with immutable declared-target semantics and configurable support limits |
| P-positive §9 and P-log §10.3: transport arms | Include in planned scientific comparison; do not block core implementation on them |
| Both plans' flux and EM alternatives | Retain as reference/deferred methods under D10; no new primary backend requirement for EM |
| D-review decision 1 | Adopt log-first direction; qualify evidence through Section 10 |
| D-review decision 2 | Override delayed derivative fix: first fitted log implementation includes it |
| D-review decision 6 / §5 | Override mandatory $10^{-4}$ flooring language; retain explicit estimator/provenance requirement |
| D-review §1.2–1.3 / E5–E6 | Replace claims of guaranteed interior/boundary optima with qualified experimental observations |
| D-review §3: no T3 bound derivable | Correct using Section 9.2; a finite triangle bound does follow from a specified raw-rate box |

## 13. Remaining choices and change-control rules

The coordinate starting choice, immediate derivative support, T3 comparison role, and target-aware primary initialization are resolved enough to proceed. Do not ask the user to decide them again before implementing the mathematical core.

The following remain configuration or collaboration decisions, not hidden assumptions: public model/option naming; actual target source/vector for a dataset; separate versus shared trees; selected rate families; root search and inferential root-split policy; final bounds/tolerances; computational budget and benchmark datasets; placement of outer orchestration; and eventual public exposure of alternate charts. Use this record's defaults where it supplies them, and record dataset-specific choices before scientific fitting. Do not invent unavailable data or claim maintainer approval.

Reopen a decision only with a concrete reason: failed correctness gate, changed upstream architecture, replicated numerical evidence, incompatible scientific target, or explicit user/maintainer direction. A decision revision must identify the affected ID, supporting evidence, effect on estimand versus implementation, regression implications, and whether published/ongoing comparisons need rerunning.

For a routine handoff, an implementing agent should report: the pinned source baseline; which G-gates are complete with evidence; which failures/limitations remain; any changed D-decisions; and the next concrete deliverable. An attractive score or successful program exit is not sufficient acceptance evidence.

## 14. Source anchors and reproducibility identity

The source guide was used to locate code. The following pinned source files/symbols support the implementation contracts; links are stable for the audited commit.

| File / symbols | Relevant finding |
|---|---|
| [utils/optimization.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/optimization.cpp): `derivativeFunk`, `lnsrch`, `dfpmin`, `minimizeMultiDimen` | Relative probes, clamping, ignored failure flag before small-step return, coordinate-dependent stopping, restart wrapper |
| [model/modelmarkov.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelmarkov.cpp): `setReversible`, `getVariables`, `targetFunk`, `decomposeRateMatrixNonrev`, `computeStateFreqFromQMatrix` | Seed conversion, coordinate/rate mutation, frequency threshold, reconstruction and equilibrium guards, QR stationary solve |
| [model/modelmarkov.h](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelmarkov.h) | Rate bounds and virtual model hooks |
| [model/modelprotein.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelprotein.cpp): `ModelProtein::init`, `rescaleRates` | NONREV seed order and positive-rate scale convention |
| [model/partitionmodel.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/partitionmodel.cpp): constructor, `targetFunk`, `optimizeLinkedModel`, `getNParameters`, `optimizeParameters` | Name-based linking, target pooling guards, joint optimizer object, counting and scheduling |
| [model/partitionmodelplen.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/partitionmodelplen.cpp): `optimizeParameters` | Edge-proportional nuisance scheduling |
| [model/modelfactory.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelfactory.cpp): constructor and optimization | Joint frequency parsing, rate/nuisance and conditional root updates |
| [model/modelliemarkov.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelliemarkov.cpp): packing and `setRates` call | Live precedent for separate coordinate storage |
| [model/modelunrest.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/model/modelunrest.cpp): checkpoints | Lifecycle precedent, not a complete constrained-model implementation |
| [main/phylotesting.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/main/phylotesting.cpp): `mixRevNonrev`, candidate generation | Mixed-family guard and candidate handling |
| [alignment/alignment.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/alignment/alignment.cpp): `convertCountToFreq`, `convfreq` | Iterative ambiguity treatment and default frequency adjustment |
| [utils/tools.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/tools.cpp): `--min-freq`, defaults; [utils/eigendecomposition.h](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/utils/eigendecomposition.h): `ZERO_FREQ` | Configurable acceptance threshold and separate small-frequency threshold |
| [tree/phylotree.cpp](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylotree.cpp), [tree/phylokernelnonrev.h](https://github.com/iqtree/iqtree3/blob/63c330d90dd02241dbbbaf1e9f9e9cc6dadbd1de/tree/phylokernelnonrev.h) | Root representation/counting and root-frequency likelihood semantics |

SHA-256 fingerprints identify the exact supplied snapshots audited here. If a later copy differs, reread the relevant changes instead of assuming this adjudication automatically covers them.

| Input | SHA-256 |
|---|---|
| `constrained_nq_plan_agent(1).md` | `5125f4626f86a73cdbc0965375610d8bf03cc679d16899224116c123e712de43` |
| `project2_iqtree3_technical_plan.md` | `953438a82db410ca8c474d433ddf89f5e4208d3d743d5318a824c40c4349ad5b` |
| `discrepancy_decisions.md` | `030555bb1730190d3672b89f5d17c539b18266eafcb61c4873f52aec158c75b1` |
| `idenfitied_agent_plan_differences_a.md` | `7743311af7f70794c452607a80fdb62cd9d43072b288d96f2c2c91f5f0da2ef6` |
| `identified_agent_plan_differences_b.md` | `565a84e5609a612eb316f5ec91e32ee8f12eab5cee83416b39b2a9c05726abec` |
| `AA_MODEL_INFERENCE(3).md` | `ada51b6e492c74eed69584213f7d451ccd5b2285d27fba2648089672e6699a76` |
| `chart_fd_bfgs_compare.py` | `8f2a1d095b24efcef62cff881cea213a924d4af7f14279517ee4d90aa9d0591b` |
| `e6_n20_big.py` | `4eb138f59f1545e65b5d887e868f11ee370d67d7338d5b367fff19ccd8ee2268` |

### Review record

This synthesis was reviewed in separate passes for (1) mathematical/model distinctions and boundary qualifications, (2) source/script fidelity and evidence status, and (3) implementation handoff consistency and preservation of alternatives. Revisions made during those passes clarified the index/edge convention, numerical diagonal construction, input parsing contract, and separation of validation tuning from held-out testing. The reviews specifically checked that log-first does not imply positive-ratio removal; that derivative support is not deferred; that target flooring is not silently mandated; that all linked objects share one coordinate meaning; and that unavailable or toy evidence is not presented as compiled production validation. Structural checks confirmed the decision IDs, internal section references, input fingerprints, Markdown delimiters, and E2 table values against the actual rerun output.
