# Decisions

The active rules for this project's code, tests and method. Each entry states its rule as it
stands now. A superseded rule is rewritten in place and its "Since" line updated; the earlier
wording is in git and in `docs/agent/archive/DECISIONS_2026-09-23_to_2026-10-22.md`, which also
keeps the full reasoning, evidence and source anchors of every entry up to that date. Numbers are
permanent and never reused; the number map at the end gives the home of each retired number.

An entry is added only when Peter chooses between real alternatives and the choice constrains
later code, tests or method (CLAUDE.md, "Records"). Where the code does not yet meet a rule, the
gap is listed in `CODE_PLAN.md` section 4 under the slice that closes it, not here.

The D entries carry the IDs of the design synthesis
(`design/project2_unified_discrepancy_synthesis.md`, section 1). A change to one names, in its
"Since" line, the items of the synthesis's section 13: the affected ID, the evidence, whether the
estimand or only the implementation changes, the regression implications, and whether finished
comparisons must be rerun.

Entry format: **Rule**, **Why**, **Not chosen** or **Retained**, **Since**; at most 12 lines,
with symbol names rather than line numbers (`CODE_PLAN.md` section 5 holds the anchors). The
family throughout: 20-state generators with positive off-diagonal rates, zero row sums, π*Q = 0
and unit mean rate at π*, 360 free parameters, containing the 189-parameter reversible family.

## Principles

### 001. Develop on the fork's `nq-constrained-pi` branch, for an eventual upstream pull request

**Rule.** All work happens on `nq-constrained-pi` of `petergoodman/iqtree3`; nothing is committed
to `master` or pushed to `iqtree/iqtree3`. The work is to be offered upstream after testing, so
every change to a shared IQ-TREE file is opt-in or behaviour-preserving.
**Why.** Peter's intent. A `master` that mirrors upstream keeps syncing a fast-forward, and opt-in
changes keep the pull request reviewable.
**Not chosen.** A permanent private fork, which would allow invasive edits but forfeit
contribution.
**Since** 2026-09-23. Git procedure: `CLAUDE.md`.

### D06. Change only the map from coordinates to generator, and no other IQ-TREE behaviour

**Rule.** IQ-TREE's likelihood engine, transition matrices, pruning, and branch, rate and tree
algorithms are reused. The substantive change is the map from 360 optimizer variables to a valid
Q in nQMaker's shared-matrix update (step 3a), with π* as the root distribution. Target, chart
metadata and references survive construction, linking, checkpoint and export, and linked objects
agree on them before optimization. No IQ-TREE behaviour changes beyond what the constraint
requires, even where a change would improve IQ-TREE. Outer reselection is orchestrated explicitly.
**Why.** Synthesis section 7; the maintainers' advice of caution with `ModelMarkov` and to leave
the optimizer unchanged; a small pull request. ModelFinder rejects sets that mix reversible and
non-reversible models (`mixRevNonrev`).
**Retained.** A thin external driver for the outer loop first, native integration later.
**Since** 2026-09-23; revised 2026-10-22 to state the rule Peter set on 2026-10-05 and applied in
027 and 029 (wording only; estimand unchanged, no regression or rerun).

## Method

### D01. Parameterize the family by jump-chain log-ratios, the only compiled chart

**Rule.** Per row, 18 log-ratios against a fixed reference destination; a stable row softmax
gives the jump matrix K; one stationary solve νK = ν; exit rates ν_i/π*_i. References are the
row maxima of LG rebuilt at π*, ties to the lower amino-acid index, shared by every linked object
and start (011). No positive-ratio or T3 chart, benchmark or chart selector is compiled; the
oracle's positive-ratio and T3 code stays as a Python test reference.
**Why.** Both charts cover the family one-to-one (synthesis section 3.1); log coordinates scored
higher in three of four eight-state fits (sections 10.3, 10.4), evidence for a starting choice,
not a theorem. The decks the maintainers accepted present this chart only.
**Not chosen.** Compiled positive weights as a benchmark; T3 as a second backend (former D04).
**Since** 2026-09-23; revised 2026-10-07 by 028 (D01 and D04; Peter's direction; implementation
only, estimand unchanged; no regression or rerun).

### D02. Use a scale-aware finite-difference step

**Rule.** NQC's derivative step is h_k = η max(s_k, |x_k|) with s_k = 1 and η = 1e-4 (a
reproducible starting setting, not a validated constant), reached from both the single-model and
the linked path (005). Existing models keep the legacy step.
**Why.** IQ-TREE's step 1e-4 |x_k| collapses near zero for signed coordinates: one logit's
relative derivative error was about 9.6e-2 at z = 1e-6 and 6.4 at z = 1e-8 (synthesis section
10.2). Linked training optimizes the `PartitionModel` object, which a class override misses.
**Retained.** The legacy rule as a diagnostic control; analytic gradients later.
**Since** 2026-09-23.

### D03. Declare numerical domains explicitly and report their effect

**Rule.** Coordinate bounds are explicit and recorded for every run (default and option: 007);
the default is a benchmark domain, not a chosen production range. Bound activity is reported,
final scientific runs include a domain-expansion sensitivity check, and an imported matrix is
never clipped to a box. A start is not checked against the box.
**Why.** A finite box restricts the family, and boxes in different charts are different sets of
matrices (synthesis section 4.2).
**Retained.** Wider limits when bounds are active or sensitivity checks move the solution; a
boundary-capable method if exact zero rates are needed (D10).
**Since** 2026-09-23; start clause removed 2026-10-05 by 027 (D03; Peter's direction;
implementation only, box and estimand unchanged; no regression or rerun).

### D05. Build the first seed and the guide candidates at the target

**Rule.** The default start is LG rebuilt at the target, q_ij = R_ij π*_j at unit mean rate. The
fitted `GTR20+F{π*}` is an additional start and the nesting incumbent (025). Guide candidates are
LG, WAG and JTT rebuilt at π*, with empirical and optimized frequency variants suppressed (S6).
**Why.** One composition convention across every stage; guide candidates only select models and
trees per alignment (synthesis section 6.2). LG normalized at its own frequencies is a valid
start whose exchangeabilities are not LG's.
**Retained.** Empirical guide candidates and transported starts, as sensitivity arms.
**Since** 2026-09-23.

### D07. Take one explicit, immutable target and never floor it silently

**Rule.** The fitter consumes one explicit, strictly positive, normalized π* in IQ-TREE's state
order, and never re-estimates, floors, renormalizes after construction or replaces it. Its
provenance is the pooled composition of the original training data through IQ-TREE's estimator,
computed once and held fixed across cleaning treatments, with estimator, ambiguity and gap
treatment, smoothing and exact values recorded. An estimator floor is a declared procedural
choice; values the backend cannot support (below `min_state_freq` or `ZERO_FREQ`) are refused.
**Why.** Synthesis section 8.
**Retained.** Other declared estimators, and supported small positive targets.
**Since** 2026-09-23.

### D08. Plan reversible, unrestricted and transport comparison arms

**Rule.** The scientific comparison includes original-data and cleaned-data reversible and
unrestricted fits, the cleaned reversible fit with frequency replacement, the reversible refit at
π*, the general refit at π* (this project's estimate), and the cleaned unrestricted fit
transported to π* by fixed flux and by T3 without refitting. Frozen matrices are evaluated on the
same held-out data. The transport arms need no compiled T3 and do not block the implementation.
**Why.** Since π* e^{tQ} = π*, every tip marginal is fixed and a wrong target can produce
apparent circulation; simulations at correct and perturbed targets come before circulation is
read as biology (synthesis section 9.1).
**Since** 2026-09-23.

### D09. Adopt the corrected mathematical qualifications

**Rule.** These override conflicting prose in the older design documents. (1) With raw rates in
[l, u], the gauge-fixed T3 triangle tilt satisfies |h'| ≤ (3/2) log(u/l), 20.72 for [1e-4, 100],
so a bound of 12 does not contain that box. (2) q_ij ≤ 1/(2 π*_i). (3) Suprema over nested
domains are ordered but independently attained local fits need not be, so the smaller model's fit
and nuisances are imported as an incumbent. (4) An interior true parameter does not imply an
interior sample maximum. (5) The root-split singularity under reversibility means the 171 and
19 parameter differences need not give chi-square references. (6) Frequency replacement in
asymmetric factors can give a valid generator but loses target stationarity; fixing root
frequencies alone is not the constraint. (7) No global optimum is claimed; the BFGS inverse
Hessian is not a covariance.
**Why.** Synthesis section 9.2.
**Since** 2026-09-23.

### D10. Keep direct balanced flux as an independent reference and defer EM

**Rule.** Direct balanced-flux fitting, f = f0 + Zy with a feasible constrained solver, is a
small-problem reference (S5) and a possible boundary-capable backend. Constrained
substitution-history EM is deferred; if built, the normalized M step is preferred. Penalty and
projection methods do not substitute for the constrained maximum-likelihood estimate.
**Why.** Synthesis section 9.3.
**Retained.** EM, for persistent optimization failures, exact boundary fits, or gradient
infrastructure.
**Since** 2026-09-23.

## Model class

### 002. NQC is a `ModelProtein` subclass built through a protected constructor

**Rule.** The model is `ModelNonrevFixedEq` (`model/modelnonrevfixedeq.{h,cpp}`), public name
"NQC", provisional and held in one constant. `createModel` dispatches the name at the top of its
chain, so non-protein data reach NQC's own refusal. The base is built through a protected
`ModelProtein` constructor that reads no model; NQC then reads its seed as NONREV's branch of
`ModelProtein::init` does (the `--init-model` name or file, else LG, tree adaptation off), without
`rescaleRates`. `NONREV` and `ModelProtein::init` are unchanged.
**Why.** Tip likelihoods for B, Z and J must match NONREV and GTR20 for nesting; `createModel`
holds the models block. The public constructor reads a named model with tree adaptation on, which
unroots a rooted `-te` tree (against 016). A common seed scale changes nothing beyond rounding.
**Not chosen.** A name branch in `ModelProtein::init`; a `ModelMarkov` subclass; a base built
under the name NONREV.
**Since** 2026-09-23; revised 2026-10-20 (031, 034 item 6) and 2026-10-21 (035 item 1).

### 003. The 360 coordinates are the model's state, checkpointed losslessly

**Rule.** NQC's 360 log-ratio coordinates are its authoritative state; `rates[]` and
`rate_matrix` are derived. The class saves under its own keys, at 17 digits: chart and version,
target (from `state_freq` from S3a), references, coordinates, domain, derivative policy, fixed
status, and the accepted fit's scores, bounds and residuals with fit and evaluation counts (S3a).
Restore verifies target, references, chart and domain, rejects a mismatch, recomputes the derived
bounds and residuals and keeps the larger count; with none of the keys it leaves the model as is.
**Why.** The pinned quantity is a jump weight, not a rate (synthesis section 7.1). IQ-TREE writes
checkpoint doubles at 10 digits, and the default +I+G path saves and restores the model during
ordinary fits, so without these keys the restored Q would not match the restored nuisances.
**Not chosen.** Coordinates inside `rates[]`; a global rise of `CKP_PRECISION`; waiting for S3.
**Since** 2026-09-23; brought into S2 by 033 (2026-10-20); fit record 2026-10-22 (S3 plan).

### 005. Reach the scaled derivative step through both optimizer paths

**Rule.** NQC overrides `derivativeFunk` for single-model fits, and `PartitionModel` gains a
`derivativeFunk` override that delegates to `LinkedModelHooks` when the linked model implements
it, falling back to `Optimization::derivativeFunk` for every other model (S3b). The step is D02's,
taken backward when the forward probe would leave the domain; a central stencil exists for
validation only. From `optimizeLinkedModel` the interface also takes a start call, which checks
linked compatibility, and an end call, which runs the post-fit check (009).
**Why.** Linked optimization runs `minimizeMultiDimen` on the `PartitionModel` object, which no
model override reaches; delegating leaves `utils/optimization.cpp` untouched. Only
`optimizeLinkedModel` knows when a linked fit ends and what the optimizer returned.
**Not chosen.** A step hook inside `Optimization::derivativeFunk` (`-jointopt` is refused, 008);
reporting linked fits only at the class's report and export, which can see neither.
**Since** 2026-09-23; end call kept 2026-10-22 (S3 plan).

### 006. Take the target through the existing `+F` routes and validate it strictly in the model

**Rule.** π* arrives as `NQC+F{p1,...,p20}` or as `NQC+F<NAME>` defined in an `--mdef` file. NQC
parses it itself: exactly 20 finite entries, each above 0 and at least D07's minimum; a sum
within 1e-6 of 1 is normalized once and both vectors reported; anything else is refused, and
nothing is floored or substituted. Every `+F` form without a vector, and a missing `+F`, is
refused with a pointer to `NONREV`. Reports and exports print π* at 17 significant digits. The
target never reaches `ModelMarkov::readStateFreq`.
**Why.** Both routes reach the model with no shared change. IQ-TREE's reader accepts NaN, draws
random values for non-numeric tokens and renormalizes with a warning. The 1e-6 line accepts
IQ-TREE's 8-decimal print; Peter accepted the departure because it lives only in NQC's code.
**Not chosen.** A `--target-freq FILE` option (a maintainer question, PLAN.md); `NQC{...}` braces.
**Since** 2026-09-23; sum line confirmed 2026-10-05 (026).

### 007. Domain: the benchmark box by default, a parse-only option, no start check

**Rule.** The default box is z in [log 1e-5, log 10], the per-row form of IQ-TREE's rate bounds
[1e-4, 100] with a protein start scaled to a largest rate of 10 and the reference weight at 10.
`--nqc-domain LO,HI` (provisional name, ratio units) sets it per run; a step option is added only
if a step study needs it. Both are parse-only and read through `Params`. The start is not
compared with the box: it is evaluated as given, and the line search clamps later trial points.
**Why.** D03 needs explicit, variable domains; the box reuses IQ-TREE's constants, and the chart
is total on R^360, so a start outside the box still gives a valid Q (D06).
**Not chosen.** Compiled constants; a domain per start; a start margin; [log 1e-5, log 1e5].
**Since** 2026-09-23; start check removed by 027 (2026-10-05). Whether this default suits the
log-ratio chart is a required full decision, open (PLAN.md risk 22).

### 008. Refuse every combination the first release does not support

**Rule.** NQC stops with a message for every input outside the first release; the full list is
`CODE_PLAN.md` section 2.2, "Refusals". It covers non-protein data, untargeted `+F`,
`-optfromgiven`, `-jointopt`, mixtures, tree mixtures and HMM, GHOST `+H` (also in a partition's
own model string, S3b), site-specific frequencies, `--link-exchange-rates`, `--eigen`, ModelFinder
in any form (detected at run time from S3a), AliSim, unlinked partitioned use (S3a), and until S3b
`--model-joint` and `--link-model`. Mixture components are refused where `initMixture` makes them.
**Why.** Each path bypasses the invariants or is undesigned: `-optfromgiven` releases the
equilibrium in the base decomposition (015), `-jointopt` clamps coordinates with rate bounds, and
an `--mdef` mixture never shows "MIX" in the model string. Unlinked partitions lie outside the
first release, so Peter chose refusing them over testing them (2026-10-22).
**Not chosen.** Allowing them silently; supporting them now; detecting mixtures from the string.
**Since** 2026-09-23; revised 2026-10-20 (032; 034 items 2, 4, 5) and 2026-10-22 (S3 plan).

### 009. Check every fit afterwards, report a fall below the start, restore nothing

**Rule.** After each fit NQC re-evaluates the returned point, checks the invariants, records bound
activity and reports any fall below the start, keeping the optimizer's point as every IQ-TREE
model does. At export a central-difference first-order diagnostic in log-ratio coordinates labels
the stop "unclassified" unless it passes. `utils/optimization.{h,cpp}` is never edited. A linked
round that lowers the log-likelihood by more than 0.1 stops at IQ-TREE's `ASSERT`, as for NONREV.
**Why.** A return from `dfpmin` does not certify convergence (synthesis section 5).
`ModelMarkov::optimizeParameters` keeps the optimizer's point without comparing it with the
start; the decks' "Final fit" row asks for a check, not a restore.
**Not chosen.** Restoring the start; status fields in the shared optimizer; repairing the line
search (a maintainer question).
**Since** 2026-09-23; restore and optimizer fields removed by 029 (2026-10-07).

### 010. Export the fitted matrix at 17 digits in the `-m FILE` format

**Rule.** After the final model optimization, one call at the top of `reportPhyloAnalysis` in
`main/phyloanalysis.cpp`, active only when an NQC model exists, prints NQC's end-of-run summary
(S3a) and writes `<prefix>.NQC.qmat` (normalized Q and the π* row at 17 significant digits) and
`<prefix>.NQC.info` (target, state order, chart and version, references, domain, derivative
policy, residuals, bound activity, stop classification, seed, source commit) (S4). Re-import uses
the existing `-m FILE` path, its equality tested at Level 1; the `.iqtree` block is unchanged.
**Why.** The report prints Q at 6 digits, and re-import from it reproduced the log-likelihood only
to 4e-4; the existing reader loads a full non-reversible Q with its π as root. No model prints its
state on the `+I+G` path at default verbosity. At Level 1 `accepted_score` is the state's score.
**Not chosen.** Raising the report's precision for every model; a NEXUS block; a line per fit at
default verbosity; a final-score key (Peter, 2026-10-22).
**Since** 2026-09-23; summary and test level 2026-10-22 (S3 plan).

### 011. Every NQC object recomputes its reference destinations

**Rule.** Each NQC object computes its references from LG rebuilt at π*, read from the built-in
models block, as row maxima with ties to the lower index, whatever its start. References go into
the checkpoint, export and report, and are verified on restore and across linked objects before
optimization. User-supplied references are deferred.
**Why.** D01 needs one definition for every object and start; a function of π* gives it without
shared mutable state inside the OpenMP loop.
**Not chosen.** Copying references from a first object.
**Since** 2026-09-23.

### 015. Build Q in the class's decomposition override and call the unchanged base

**Rule.** NQC's `decomposeRateMatrixNonrev` override builds Q from the trial coordinates in
per-object workspace, validates it, writes the 380 off-diagonal rates into `rates[]` and calls
the unchanged `ModelMarkov::decomposeRateMatrixNonrev`, then stops with an error unless
`state_freq` is bitwise π*. A failed build writes nothing and makes `targetFunk` return 1e30.
`getStateFrequency` returns π* without the base's renormalization. `ModelMarkov` is not edited.
**Why.** With `FREQ_USER_DEFINED` and `-optfromgiven` refused, the base skips the π reset and
solve. Every decomposition passes through the override, so a 10-digit restore of `rates[]` never
becomes Q. The maintainers advised caution with changes to `ModelMarkov`.
**Not chosen.** A shared helper split out of the base (former 004); building Q in `getVariables`.
**Since** 2026-10-03, replacing 004; root guard added 2026-10-20 (034 item 1).

### 016. Root the tree by IQ-TREE's conversion and do not search the root in S0 to S4

**Rule.** Root frequencies: π*. An unrooted input tree is rooted by IQ-TREE's conversion
(midpoint of the longest path, or the `-o` outgroup's branch); a rooted one keeps its root. No
root-edge search: fits run under `-te` without tree search or `--root-find`. The root split is
fixed in Level 1 (`-blfix`) and left to IQ-TREE's branch-length optimization in Level 2. The
policy for scientific runs is set at gate G5.
**Why.** A non-reversible likelihood depends on the root, so comparisons need a stated policy
(synthesis section 7.5); this inherits nQMaker's practice under `-te`. S0 probe (f) found the
root unchanged.
**Not chosen.** Searching the root; requiring an outgroup; fixing the split in Level 2.
**Since** 2026-10-03.

### 024. Run Level 1 with existing options

**Rule.** A Level 1 fit uses no new option. One alignment and `-S`: `-te` with the input trees,
`-blfix` and braced rate parameters. Edge-proportional: `-q` in place of `-p`, with `-te`,
`-blfix` and a partition file giving each partition its braced rate model and its rate as `{x}`,
taken from the preceding `-p` fit's checkpoint (`PartitionModelPlen!part_rates`). Level 2 keeps
`-p` and `-S`.
**Why.** S0 probe (c): under `-p` the partition rates are refitted even with `-blfix`; under `-q`
with given rates every nuisance stayed fixed, on the same edge-proportional tree.
**Not chosen.** A new option that freezes partition rates under `-p`.
**Since** 2026-10-05.

### 025. Carry the GTR20 incumbent from a fixed-tree fit's checkpoint

**Rule.** The `GTR20+F{π*}` incumbent is fitted on fixed trees (`-te`), never after a tree
search, and carried from its checkpoint (190 rates and gamma shape at 10 digits) as a PAML-format
`-m` file with `+F{π*}` exact. The test evaluates it with `--show-lh` on the same trees and rates
and starts NQC from it with `--nqc-domain` widened to contain the start (`CODE_PLAN.md` section
3.5): NQC's first log-likelihood equals the evaluation within 1e-8 × |lnL|, and its final one
falls below it by no more than that. The domain is a test setting only (007 is unchanged); the
same fit in the default box is reported.
**Why.** A tree search leaves the checkpoint's model stale (S0 probe (g)). The 10-digit carry
moves the log-likelihood by about 3e-12 relative, kept out by re-evaluation. In the default box
the first line search clamps the start's 146 of 360 outside coordinates (`aa_example`).
**Not chosen.** Carrying from the report; the printed log-likelihood; the default box alone.
**Since** 2026-10-05; domain 2026-10-22 (Peter). The default box itself is open (PLAN.md risk 22).

## Tests

### 021. Pass lines are set before their tests and never loosened

**Rule.** The lines below are fixed before the tests that use them are written, and none is
loosened after a failure. The class's own state (final coordinates, π*, references, accepted
score) is read from its 17-digit checkpoint keys, not from the log; IQ-TREE's printed final score
is reported, not compared. The class prints the seed Q (`-v`) and its first gradient (`-vv`)
only because no file holds them.

| Check | Line | Set by |
|---|---|---|
| Oracle Q on ordinary cases: ‖π*Q‖∞, \|mean rate − 1\|, scaled row sums | ≤ 1e-10 | 021 |
| Oracle round trips Q → z → Q; log chart against T3 | ≤ 1e-10 relative per entry | 021 |
| Oracle builder derivative against a central difference | ≤ 1e-7 relative to max(1, max\|dQ\|) | 021 |
| Oracle hard cases (near-reducible, skewed π*, wide rates) | recorded, not asserted | 021 |
| Documented numbers | 360, 189, 171 exact; 20.72 and 0.015378 to printed digits (019) | 021 |
| Pruning against enumeration; sum over all patterns against 1 | ≤ 1e-12 relative | 021 |
| Reversible Q under a moved root | ≤ 1e-10 relative | 021 |
| Ported gamma rates against SciPy | ≤ 1e-5 relative | 021 |
| Initial lnL (`--show-lh`, precision 17) against the oracle, legacy models or NQC | ≤ 1e-8 × \|lnL\| | 021, 030 |
| Estimator port against IQ-TREE's "Mean state frequencies" | equal at the 8 printed digits | 021 |
| Compiled residuals | ≤ 1e-10; diagnostics only, never a rejection in a fit | 026 |
| Compiled Q against fixtures and the oracle; round trips | relative max(1e-12, 1e-14 κ) per entry, κ the 2-norm condition number of [1ᵀ; (K − I)ᵀ] | 026, 030 |
| Gradient against the oracle's, same step rule | \|g_k − g_k^oracle\| ≤ 2 × 1e-12 × \|lnL\| / \|h_k\| | 030 |
| Gradient against the central reference (steps 1e-2, 1e-3, 1e-4 × max(1, \|x_k\|); reference and curvature at 1e-3) | \|g_k − g_ref\| ≤ 2 [(\|h_k\|/2)\|f″_kk\| + 2 × 1e-12 × \|lnL\| / \|h_k\|] | 030 |
| NQC at a fixed point, more threads against `-T 1` (`CODE_PLAN.md` 3.5): log-likelihood; gradient | ≤ 1e-8 × \|lnL\|; \|g_k − g_k(`-T 1`)\| ≤ 2 × 1e-12 × \|lnL\| / \|h_k\| | 038 |
| Toy problem (`CODE_PLAN.md` 3.3): compiled final score against the oracle's port (017) | max(10 × spread, 1e-8 × \|lnL\|); spread is the largest final-score difference among oracle fits differing only by rounding (every start coordinate moved one unit in the last place; an LU solve for ν in place of the pivoted QR), measured by a kept script before the comparison runs | 030 |
| Level 1 accepted score against the oracle's re-evaluation of the checkpointed coordinates | ≤ 1e-8 × \|lnL\| | 030 |
| π* in the checkpoint after every fit | bitwise equal to the parsed, normalized vector | 030 |

h_k is the step's displacement (`CODE_PLAN.md` section 2.1), and the oracle reports the spread of
its three central differences. The step-size line allows more for rounding than the rounding it
bounds at every h, so it cannot show that a small step fails; D02's point is measured instead by
|g − g_ref| / max(1, |g_ref|), reported at a near-zero start.
**Why.** A line chosen after seeing results is not a test. Each line rests on a measurement or
derivation kept in the archive: 1e-10 is the synthesis's feasibility goal (section 4.3), the
κ-scaled Q line tracks the measured gap between two correct floating-point routes, and the
gradient lines turn a likelihood difference into a gradient difference.
**Since** 2026-10-04 (021), 2026-10-05 (026), 2026-10-20 (030, 034 item 7); a claim in 030's
reasons withdrawn 2026-10-22 (036); thread row 2026-10-22 (038).

### 012. Test the pure module with googletest in a standalone CMake project

**Rule.** C++ unit tests live in `unittest/`, a standalone CMake project that compiles
`model/fixedeqchart.cpp` directly and fetches googletest at cmaple's pinned commit. The class is
tested through the binary. Fixtures are plain text at 17 significant digits, written by the
oracle with a provenance header, and never regenerated to make a test pass. Integration into the
main build is revisited only if class-level C++ tests become necessary.
**Why.** The pure module needs only Eigen, and a standalone project changes no shared build file;
cmaple's test targets do not exist on Windows.
**Not chosen.** cmaple's googletest targets; an opt-in test directory in the main build.
**Since** 2026-09-23.

### 017. Port IQ-TREE's optimizer into the oracle from the C++ source

**Rule.** The oracle's optimizer ports `minimizeMultiDimen`, `dfpmin`, `lnsrch`, `fixBound` and
the legacy `derivativeFunk` from `utils/optimization.cpp` at `63c330d9`, with IQ-TREE's constants
and control flow, the failed-line-search path included, called as
`ModelMarkov::optimizeParameters` calls it. The scaled step is a separate function, and the
stationary solve is a column-pivoted QR as in the C++. Floating-point identity is not expected.
**Why.** The design script's port returns the old objective on a failed line search and solves by
least squares, so a disagreement with it could come from the port, not the new code.
**Not chosen.** Reusing the design port; editing the design script.
**Since** 2026-10-03.

### 018. Judge baseline run 11 by its fixed parts, and add run 12

**Rule.** Run 11 (`--model-joint NONREV`, `-p`, `-T 4`) is compared on the exit code, the
REFERENCES and SEQUENCE ALIGNMENT sections exactly, and the SUBSTITUTION PROCESS section and tree
file with every number removed; its tree drawing is not compared. Run 12, run 10 at `-T 4`, is
recorded five times from the frozen binary and compared by `CODE_PLAN.md` section 3.4's rules.
**Why.** Five repeats of run 11 ended at different log-likelihoods, none equal to `-T 1`, while
wording and topology never changed; run 12 fixes tree and matrix and exercises the threaded
likelihood.
**Not chosen.** Holding run 11's numbers to their spread; dropping run 11.
**Since** 2026-10-03.

### 038. Compare thread counts at a fixed point, not after a fit

**Rule.** NQC's linked log-likelihood at the seed and its first gradient are computed at `-T 1`
and with more threads, under `-p` and `-S`, with at least as many partitions as threads so that
partitions run concurrently (on turtle, `-p` at `-T 4`, which IQ-TREE caps there at one thread
per partition, and `-S` at `-T 3`). Each multi-threaded case runs 10 times, held to 021's lines,
not bitwise, with `-T 1` as reference. Fits at more than one thread are reported, not asserted.
**Why.** On the unmodified binary (2026-10-22) evaluations at a fixed point were exact at `-T 1`
and differed by at most 1.3e-16 relative at `-T 4`, while none of 18 joint NONREV fits at `-T 4`
matched `-T 1` (spreads 0.27 to 0.66, Level 1 included): the optimizer amplifies the rounding. A
thread-safety fault in NQC is expected to miss 021's lines by orders of magnitude.
**Not chosen.** Fit equality; a tolerance from the fits' spread; no test; a ThreadSanitizer build,
the fallback if this test fails.
**Since** 2026-10-22, retiring PLAN.md risk 19 (measurements: `explore/thread_counts.py`).

### 023. Use the turtle pooled composition as the S1 to S4 test target

**Rule.** The named test target of S1 to S4 is the pooled composition of
`test_scripts/test_data/turtle_aa.fasta` over the charsets of `turtle_aa.nex`, computed by
`oracle/target.py` under `-p` conventions and recorded at 17 digits in the G0 manifest. It serves
the `aa_example` and turtle tests; random targets in the oracle's tests and fixtures are separate
and recorded by seed. S3 adds one reported, not asserted, fit at this target cubed and
renormalized, whose LG start has 1 of 360 coordinates below the default box.
**Why.** The port reproduces IQ-TREE's printed pooled vector for this data, and a pooled
composition is the kind of target the study uses (D07).
**Not chosen.** `aa_example`'s own composition; an asserted skewed-target test.
**Since** 2026-10-04; skewed case 2026-10-22 (`explore/skewed_targets.py`).

## Tooling

| No. | Rule | Since |
|---|---|---|
| 013 | The oracle, drivers and fixture generator live in `test_scripts/fixedeq/` and run in the conda environment `iqtree3-fixedeq` in WSL, built from `environment.yml` (conda-forge, pinned) with an exported lock. The system Python is never used. | 2026-09-23 |
| 014 | Drivers are Python files; one-off `wsl.exe` command scripts use LF endings. The new directories carry their own `.gitattributes`; there is no root `.gitattributes` and `core.autocrlf` is unchanged. | 2026-09-23 |
| 019, 020 | The design scripts are rerun from copies by `test_scripts/fixedeq/design_rerun.py`, each documented value matched in one of four classes fixed before any run: exact (counts, ranks, labels: equal), digits (rounded to the printed digits: equal), zero (both below 1e-12) and bound (the printed inequality holds). The 2026-10-04 rerun is the reproduced record: 109 of 133 values reproduced; the other 24 are platform-sensitive, stay "reported" and are not test targets. | 2026-10-04 |
| 022 | Sanitizer builds use IQ-TREE's `Mem` build type with clang 14 and AddressSanitizer plus UndefinedBehaviorSanitizer, in a new directory for each slice; fallbacks are a build without cmaple, then UndefinedBehaviorSanitizer alone. Findings in upstream code are recorded, not fixed. Flags and run options: `CODE_PLAN.md` section 3.6. | 2026-10-04; per-slice directory 2026-10-21 (035 item 2) |

## Number map

Retired and merged numbers, so that every citation in code, tests and records still resolves.

| Number | Now |
|---|---|
| 004 | replaced by 015 |
| 026 | pass lines in 021; target sum in 006 |
| 027 | 007, and D03's start clause |
| 028 | D01 |
| 029 | 009 |
| 030 | pass lines in 021; S2's gradient points and toy problem in `CODE_PLAN.md` section 3.3 |
| 031 | 002 |
| 032 | 008 |
| 033 | 003 |
| 034 | item 1 in 015; items 2, 4, 5 in 008; item 3 (matrix setters stop with an error) in `CODE_PLAN.md` section 2.2; item 6 in 002; item 7 in 021 |
| 035 | item 1 in 002; item 2 in 022 |
| 036 | 021 (the withdrawn claim is gone; its measure is kept) |
| 037 | not a rule: S2's shortfalls, S3-1, listed in full in `CODE_PLAN.md` section 4 at `4749f222` and closed for one alignment in S3a (section 4 now; the linked half is S3b's) |
| D04 | dropped by 028; see D01. The T3 gauge the oracle's test reference uses ("the D04 gauge") is a_{19,20} = 0 and h_{i,20} = 0 for i < 20 (189 strength and 171 tilt coordinates, 19 balancing potentials by damped Newton), synthesis section 6.1 |
