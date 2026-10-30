# Plan: π-constrained non-reversible amino-acid models

> **Status, 2026-10-30: S4 complete, so gates G0 to G3 are met and the branch is the first release
> candidate: NQC writes its fitted matrix at 17 significant digits (`<prefix>.NQC.qmat`, which
> re-imports through `-m FILE` within 7.8e-16 of the fit) and an export record with a first-order
> stop label; the release audit's blocker (threaded runs could abort in the export) is fixed and
> its other findings closed or with the maintainers; S5 (independent reference) next; approach
> accepted by the IQ-TREE maintainers (reported by Peter, 2026-10-01).** Peter owns this document.

This is the project's living core document: what we are building and why, the constraints the
work must satisfy, where it stands, and the next step. It is read at the start of every task and
rewritten as the project moves; what belongs here and what goes elsewhere is set in `CLAUDE.md`,
"Records". A task that implements or tests a slice also reads `docs/agent/CODE_PLAN.md` after
this file. Decisions are cited by number from `docs/agent/DECISIONS.md`.

Labels: (verified) means read in source at `63c330d9`, or run on 2026-09-23; (reported) means
taken from a document; (proposed) means a plan choice not yet executed.

## Goal

Add to IQ-TREE 3 a non-reversible amino-acid substitution model whose stationary distribution is
fixed to a supplied target vector π*, and estimate its rate matrix by maximum likelihood within
the existing nQMaker joint-estimation workflow (Dang et al. 2022, doi:10.1093/sysbio/syac007).
Given one strictly positive π* in IQ-TREE's state order, the model is the set of 20-state
generators Q with positive off-diagonal rates, zero row sums, π*Q = 0, and unit mean rate at π*.
It has 360 free parameters, contains every positive reversible generator with equilibrium π*
(189 parameters), and, as an unbounded model, is nested in the 379-parameter `NONREV` family.
The root distribution is π*, and π* is never re-solved from Q. The substantive change is
confined to nQMaker's shared-matrix update (step 3a): its optimizer variables become 360
jump-chain log-ratio coordinates (D01), and every trial matrix is built so that π*Q = 0 holds by
construction at every likelihood evaluation, not only at convergence. IQ-TREE's likelihood
engine, model selection, tree inference, and branch-length and rate-parameter optimization are
reused, with peripheral changes to target input, seeding, parameter counting, and output. The
estimate is the non-reversible counterpart of the reversible workflow of Wheeler et al. (2025;
doi:10.64898/2025.12.01.691663 as given in the Part 2 slide deck, not verified), which takes
exchangeabilities from filtered training alignments and composition from the target data, a
frequency replacement that is exact only for reversible models. The work is done on branch
`nq-constrained-pi` of the fork `petergoodman/iqtree3` and is intended for eventual submission
upstream as a pull request (decision 001).

## Scientific motivation and success criteria

> Drafted from the design documents on 2026-09-23; Peter confirmed the software criteria on
> 2026-10-29. The scientific criterion below is a placeholder until he replaces it.

The reversible workflow's frequency replacement is exact only because q_ij = R_ij π_j is
stationary at every positive π (reported). In an asymmetric factor the same replacement still
gives a valid generator but loses target stationarity (D09 item 6), and IQ-TREE's
`ModelMarkov::adaptStateFrequency` performs exactly that swap on a non-reversible matrix
(verified). The constrained estimate instead refits Q inside the family with π*Q = 0 and π* as
the root distribution, so the fitted matrix is stationary at π* by construction rather than
transported there.

The intended study is the Project 2 cleaning study (reported, synthesis section 8 and D08): π*
is the pooled composition of the original training data, held fixed across cleaning treatments,
and the comparison arms are those of D08. Because π* e^{tQ} = π*, every tip marginal equals π*,
so a wrong target can produce apparent circulation; simulations with reversible and
non-reversible generators at correct and perturbed targets come before any circulation is read
as biology. Manuscript title and venue: to be supplied by Peter.

Success criteria for the software (confirmed by Peter, 2026-10-29). The feature works when gates
G0 to G3 of the synthesis are met with this evidence:

- every accepted Q has ‖π*Q‖∞ and |mean rate − 1| at or below 1e-10, the synthesis's
  provisional goal, with the actual maxima reported;
- IQ-TREE's likelihood at a fixed Q equals the independent oracle's (test layer 3);
- a constrained fit seeded at the `GTR20+F{π*}` fit, with that fit's trees and rate parameters,
  never scores below it in decision 025's widened domain;
- 360 matrix parameters are counted once across linked partitions;
- π* is bitwise unchanged through construction, linking, optimization, checkpoint, restart and
  export;
- export followed by re-import reproduces the fitted likelihood, tested at Level 1 (decision 010);
- every legacy model in the regression baseline reproduces exactly;
- a linked NQC fit takes time comparable to NONREV's on the same data and settings (thresholds,
  set 2026-10-22: 1.25 per evaluation and 2 per fit to investigate, 3 per fit to block; met in
  S3c, `CHANGELOG.md`).

All eight were met by S4's close; the audit, with each criterion's evidence, is in the S4 entry of
`CHANGELOG.md`.

Scientific success criterion (proposed placeholder). Frozen matrices are evaluated on the same
held-out genes under one declared nuisance policy across the D08 arms, with uncertainty from
resampling genes. Candidate criterion: the general refit at π* exceeds both the reversible refit
at π* and the frequency-replacement arm in held-out log-likelihood by more than the resampling
uncertainty, and non-reversibility is claimed only where it exceeds what reversible simulations
at a perturbed target produce. Whatever criterion is adopted is written as a confirmatory plan in
the repository before any training or evaluation data are read.

## Scope of the first release

The first release is gates G0 to G3 (slices S0 to S4 below), which confirms the synthesis's
default: Levels 1 and 2 of P-log section 7.4 through both `-m` and `--model-joint`, protein data
only, with fixed-matrix export and re-import. Level 1 fits the shared matrix with every nuisance
held fixed (nQMaker step 3a in isolation); Level 2 is IQ-TREE's native alternation of that fit
with branch-length and rate-parameter updates; Level 3 repeats the outer candidate reselection of
nQMaker steps 1 and 2.

In the first release:

- the NQC model (provisional public name; decision 002) with the log-ratio jump chart (D01),
  scaled finite-difference derivatives in both optimizer paths (D02, decision 005), and an
  explicit, recorded coordinate domain (D03, decision 007);
- LG rebuilt at π* as the default start, and `--init-model` starts (D05);
- invocation by `-m NQC+F{...}` or `-m NQC+F<name>` for one alignment, and by `--model-joint`
  under `-S` (separate trees) and `-p` (edge-proportional), the two arrangements nQMaker uses;
- checkpoint and restart; export of the fitted matrix at 17 significant digits and re-import
  through the existing `-m FILE` path (decision 010); post-fit diagnostics (decision 009);
- rejection, with a message, of every combination listed in decision 008.

One refinement to the synthesis's default (proposed): Level 1 fits use existing options, with
`-q` and fixed partition rates in place of `-p` for the edge-proportional arrangement
(decision 024).

Tree search with NQC runs through ordinary IQ-TREE code but is not a first-release validation
target.

Deferred, still in the project: a small direct-flux reference in Python (D10; S5, gate G4 as
reduced by D01); target-rebuilt guide candidates, fixed learned candidates in ModelFinder, the
reversible/non-reversible mixing guard, an optional Level 3 outer driver, and simulation (S6,
gate G5); domain sensitivity of the scientific fits (S3a and S3b report it on test data for
risk 22), multiple starts, mismatch simulations and held-out evaluation (S7, gate G6).

Out of scope: EM (D10), exact boundary fitting, branch-specific targets, non-stationary roots,
profile mixtures with distinct equilibria, analytic likelihood gradients, compiled positive-ratio
and T3 charts (D01), and any edit to IQ-TREE's optimizer (decision 009).

## Constraints and requirements

1. Stationarity holds by construction at every likelihood evaluation, including line-search
   trials and finite-difference probes. The log-ratio chart is total on R^360, so an unclamped
   probe still gives a valid Q.
2. Off-diagonal rates stay strictly positive. Softmax underflow, a non-finite value or a
   non-positive ν rejects the proposal: `targetFunk` returns 1e30 and no state is published.
3. π* is immutable. It is never re-solved from Q, floored, renormalized after construction or
   replaced by any setter, and the root distribution equals π* bitwise.
4. Existing behaviour is unchanged by default. Every change to a shared IQ-TREE file is opt-in,
   by model type or by option, or behaviour-preserving, and the regression baseline must
   reproduce exactly at `-T 1`.
5. Parameter counts are correct: 360 when fitted, counted once across linked partitions; 0 when
   fixed; 0 frequency degrees of freedom for the external target.
6. The model's checkpoint is lossless and an incompatible restart is rejected. This matters
   beyond restarts: the default +I+G optimization path saves and restores the model through its
   checkpoint, and IQ-TREE writes checkpoint doubles at 10 significant digits (verified).
7. Work inside `PartitionModel::targetFunk`'s OpenMP loop (verified) is thread-safe: per-object
   workspace, no static mutable state, no shared warm starts, and the existing fixed-order sum
   preserved.
8. The coordinate domain is explicit and recorded, and an imported matrix is never clipped
   (D03). A start is not checked against the domain: as for every IQ-TREE model, the line search
   clamps trial points into it (decision 007).
9. `-m NONREV+F{...}` keeps its current non-stationary-root meaning.
10. Code is ready for an eventual upstream pull request (decision 001): small diffs to shared
    files in the local style, IQ-TREE's conventions (`outError`, 1-indexed optimizer vectors,
    `aligned_alloc` for kernel arrays), and every new file listed in `model/CMakeLists.txt`.

The hard rules on contacting the maintainers and on the branch are in `CLAUDE.md`.

## Where π* comes from

The target is the pooled composition of the original training data, computed once through
IQ-TREE's estimator and held fixed across cleaning treatments (D07). The estimator itself, with
its ambiguity fixed point and its optional floor, is described in `AA_MODEL_INFERENCE.md`
section 4.1. Missing sequences count as unknown states under `-p` but not under `-S` (verified),
so the partition type is part of the record, together with the gap and ambiguity treatment, any
smoothing, and the vector at 17 significant digits. For the turtle test data the difference lies
below the 8 decimals IQ-TREE prints (baseline runs 8 and 9). The real π* is computed only after
the confirmatory plan is in the repository, because computing it reads the training data. S0 to
S4 use the test target of decision 023, recorded at 17 significant digits in the G0 manifest
(`test_scripts/fixedeq/manifest/g0_manifest.md`).

The oracle's port of the estimator, `test_scripts/fixedeq/oracle/target.py`, reproduces IQ-TREE's
"Mean state frequencies" line, printed at 8 decimals, for the turtle data under both `-p` and
`-S` (verified 2026-10-05).

The model receives the target as `NQC+F{p1,...,p20}` or `NQC+F<NAME>` (decision 006). With
`--model-joint` every partition receives the same string (verified), and a check before linked
optimization requires every linked object to hold the same target bit for bit.

π* adds no degrees of freedom. Inference is conditional on a plug-in target, so adding 19 to a
BIC count when comparing fixed-target and free-target procedures is not justified (synthesis
section 8).

## Implementation slices

File-level changes, the NQC class specification, the six test layers and each slice's task list
are in `docs/agent/CODE_PLAN.md`. Each slice maps to a gate of synthesis section 11: G0 protocol
and configuration, G1 mathematical core, G2 first compiled conditional model, G3 linked and
lifecycle support, G4 competing implementations, G5 native and outer workflow, G6 scientific
validation.

| Slice | Gate | Goal | Done when | Status |
|---|---|---|---|---|
| S0 | G0 | Evidence before any source edit | environment built from its file; design scripts rerun from copies and compared with their reported values; oracle likelihood matches the unmodified binary; baseline recorded twice with the binary's SHA-256; S0 runtime checks answered; unmodified sanitizer run recorded; G0 run manifest written | done 2026-10-05 (G0 manifest) |
| S1 | G1 | Pure mathematical core | unit tests at n = 3, 4 and 20 against oracle fixtures pass in Release and sanitizer builds; residual maxima tabled | done 2026-10-19 (`6c8d3af4`) |
| S2 | G2 | First compiled model, one alignment | root, seed, gradient, target and refusal checks and the toy comparison pass under decision 021's lines; Level 1 and Level 2 fits on `example/aa_example.phy`, `+I+G4` included; regression reproduces; sanitizers clean for new code | done 2026-10-22 (`35d4a1bf`); shortfalls in S3-1 |
| S3a | G2 | One alignment complete | S3-1 done for one alignment (fit record once under `-v`, end-of-run summary at default verbosity, full checkpoint record and restore, ModelFinder and unlinked partitions refused); GTR20-seeded nesting in the widened domain (decision 025); the skewed-target fit (023) and aa_example's risk 22 runs reported; `+R4` seed against the oracle and a Level 2 fit; incompatible restarts refused; regression reproduces; sanitizers clean for new code | done 2026-10-22 (`4749f222`) |
| S3b | G3 | Linked training and lifecycle | two partitions with different compositions under `-S` and `-p`, one target and one Q in every object; each linked fit checked and reported once (S3-1's linked half); 360 counted once; the linked gradient uses the scaled step; thread counts compared at a fixed point (decision 038); permuted order; mid-run and `--undo` restarts, incompatible ones refused; turtle's risk 22 runs and the warm starts from the box fits reported; regression reproduces; sanitizers clean for new code | done 2026-10-23 (`054914de`, reruns at `fb5bf76d`) |
| S3c | G3 | Stopping, risk 22's evidence, cost against NONREV | first part: the stopping study's asserted checks and the converged risk 22 evidence (`CODE_PLAN.md` section 4); then profiling against NONREV within its blocking threshold, or, above 3 per fit, fixed in NQC's code or accepted by Peter; any ratio above 1.25 per evaluation or 2 per fit explained (risks 5, 18, 24) | done 2026-10-26 (first part 2026-10-25, domain `0e9fb159`; profiling gate at `0b0048af`) |
| S4 | G3 | Export, re-import and report | export then `-m FILE` reproduces the log-likelihood of a Level 1 fit; metadata complete; first release candidate | done 2026-10-27 (`a0ebf5cf`); first release candidate |
| S5 | G4, as reduced by D01 | Independent reference | D10's direct balanced-flux fit in Python agrees with NQC fits on small problems; comparison report | not started |
| S6 | G5 | Workflow | candidate-set and mixing-guard tests; G5 conditions | not started |
| S7 | G6 | Science | the confirmatory plan's checks, with decision 039's convergence protocol; G6 conditions | not started |

## Risks and open questions

### Open risks

| # | Item | Resolved by |
|---|---|---|
| 9 | Real training data: a small real set for S5, and the full set with a written confirmatory plan for S7 | Peter |
| 11 | Scientific success criterion and manuscript details (placeholder above) | Peter |
| 23 | Reported by the S3 plan's review, not verified: resuming an `-S` tree search with `+I+G` may restore stale top-level model keys, for every model, because per-partition refits during the search save to the top level (the model re-optimization in `IQTree::doNNISearch`) and the resumed process keeps `opt_gammai` on, so the `+I+G` starts restore from a copy of the checkpoint (`ModelFactory::optimizeParametersGammaInvar`; anchors in `CODE_PLAN.md` section 5). Outside S3's validation paths, which use `-te` | Checked before it becomes a question for the maintainers |
| 24 | Fits stop short of their peaks under IQ-TREE's stopping rules, for every model. `dfpmin` tests its gradient relative to the log-likelihood with gtol fixed at 1e-4 and starts each call with no curvature memory; the round loops stop at a gain below `-me` (0.01; ten times that under `-S`) or at a cap (10 rounds under `-S`, 99 under `-p` and `-q`). On turtle's `-q` warm start each round was one step gaining about 0.01, so where a fit stops is close to chance; NONREV's turtle `-S` Level 1 fit hit its cap still gaining 3.5 per round (maintainer question 7). Comparing two fits' scores (box against wide, NQC against NONREV, D08's arms) therefore needs a convergence protocol. IQ-TREE's optimizer stays unchanged (decision 009). S3c's first part: NQC's endpoints lie 0.77 to 4.2 below the converged peaks at default settings and 0.67 to 2.5 after 300 relaxed rounds, more than the box costs; its later rounds average 724 evaluations, one step each. S3c's profiling: NONREV's `-S` fits all ended at the 10-round cap, NQC's 10-set `-q` and `-p` fits at the 99-round cap, so time per fit mostly follows where fits stop (NQC to NONREV 0.26 to 1.09) while time per evaluation is equal | Testing: S3c's first part (done 2026-10-25). Science: decision 039 (Peter, 2026-10-25): NQC keeps NONREV's optimization process, issues in it go to the maintainers; S7's confirmatory plan declares one stopping setting for every compared fit, reports each fit's ending rule and checks its conclusions against relaxed settings or a reference climb on a subset; the climbing routine is revisited if S7's pilot shows differences between arms the size of the stopping shortfall |

Retired risk numbers: 2, 3, 4, 6, 10 and 13 to 17 were resolved by 2026-10-07 (archived
decisions and `CHANGELOG.md`); on 2026-10-22, 7 (a linked-round `ASSERT`) and 12 (the switch to
scaling and squaring) became S3 monitoring checks in `CODE_PLAN.md` section 4, 8 joined 22,
1, 20 and 21 became maintainer questions 1, 4 and 5 below, and 19 was settled by decision 038;
on 2026-10-25, 22 (NQC's production domain) was settled by decision 007's revision, default
[1e-5, 1e5] from S3c's first part's converged peaks; on 2026-10-26, 5 (cost) and 18 (no fallback)
by S3c's profiling: an NQC evaluation costs what a NONREV one does, and by labelled extrapolation
a linked NQC `+G4` fit on about a million columns at 20 taxa takes 0.8 to 5.6 hours of evaluation
time at the round caps on 32 cores, assuming the efficiency measured at `-T 4` (`CHANGELOG.md`, S3c
entry).

### Questions for the maintainers

Peter relays these if and when he chooses; no agent contacts the maintainers (`CLAUDE.md`).

1. Final public names for the model and its options.
2. Whether a dedicated target option, such as `--target-freq FILE`, would be preferred to
   `+F{...}` (decision 006).
3. Whether the failed-line-search path of `lnsrch` and `dfpmin` should be repaired; NQC reports a
   fall below the start and changes nothing (decision 009).
4. Upstream finding: `ModelMarkov`'s constructor calls `setReversible`, which reads
   `is_reversible` before anything sets it (UndefinedBehaviorSanitizer: "load of value 190,
   which is not a valid value for type 'bool'"). It appears once in every IQ-TREE sanitizer run,
   the unmodified binary's included; its effect on results was not investigated.
5. Upstream observation, S0 probe (g): after a tree search the checkpoint keeps the model from
   before the final model optimization, and a restart from it prints "Final model parameters
   restored" (the restart itself was not tested). NQC avoids it by decision 025.
6. Upstream finding, S3b's sanitizer run: at `-v` and above, `runPhyloAnalysis` ends every
   partitioned run whose partition type is not edge-unlinked (`-Q`) by casting the tree to
   `PhyloSuperTreePlen` and printing its NNI counts. An `-S` tree is not one, so
   `printNNIcasesNUM` reads past the object: AddressSanitizer reports a heap-buffer-overflow
   and stops IQ-TREE, and a Release build prints meaningless counts. Reproduced on the
   unmodified sanitizer binary with turtle, `-S`, `LG+G4{0.5}`, `-te -blfix -v`.
7. Upstream observation, S0 probe (c2), unmodified binary: NONREV's turtle `-S` Level 1 fit
   (`--model-joint NONREV -te -blfix -m LG+G4{0.5}`) ends at the default 10-round cap
   (`--loop-model`) still gaining 3.51 per round, at −4815.29, while NQC, whose unbounded family
   lies inside NONREV's, reaches −4806.09 on the same setup (risk 24). With the round rules relaxed
   (`-me 1e-6`, 300 rounds; S3c's first part, S3b's binary, NONREV's path unchanged), NONREV's turtle
   `-S` and `-q` fits stopped at rounds 41 and 47 on a round that gained exactly nothing, after
   rounds gaining 2.9e-3 and 5.1e-3; on `aa_example` it ran to the cap still gaining. In S3c's
   profiling (binary `81729f4b`, NONREV's path unchanged) its `-S` fits ended at the 10-round cap
   in all six turtle and 10-set configurations. Cause not investigated.
8. NQC's report lines (the `.iqtree` "Model of substitution", a partitioned run's
   `.best_model.nex`, the `.iqtree` AliSim command) carry neither the target nor the fitted matrix,
   so none re-runs the fit; the fork points to `<prefix>.NQC.qmat` instead (Peter, 2026-10-26).
   Would the maintainers prefer the model string to carry the target, or the report a
   full-precision matrix?
9. Not verified: with `-mstore`, `ModelFactory` serves transition matrices cached by branch length
   from the end of each model optimization, so a linked fit that follows a partition's own update
   may evaluate stale matrices, for NONREV as for any linked model. NQC refuses `-mstore`
   (decision 008).
10. Upstream finding, S4's sanitizer run: `Alignment::createBootstrapAlignment`, which
    `runStandardBootstrap` calls for each `-b` replicate, copies a `pomo_sampling_method` the
    source alignment never set (UndefinedBehaviorSanitizer: "load of value 3200171710, which is not
    a valid value for type 'SamplingType'", reported once per run). Reproduced on the unmodified
    sanitizer binary with `-m LG -b 2`; its effect on results was not investigated.
11. Upstream finding, 2026-10-28, unmodified binary: `--runs` with `-S` crashes. On the
    two-partition data with `NONREV -te -blfix --runs 2` both runs finish and the crash comes as
    the best run is restored; a backtrace of the same crash with NQC shows
    `PhyloSuperTreeUnlinked::printResultTree`, called by `runMultipleTreeReconstruction`. On turtle
    with `LG+G4 --runs 2` and a tree search it comes in run 1's parsimony start tree
    (`PhyloTree::computeParsimonyTree`). By reading, not verified: `runMultipleTreeReconstruction`
    builds each run's tree as a `PhyloSuperTreePlen` for every partition type but `-p`, `-S`
    included.
12. Upstream observation, 2026-10-29 audit: at more than one thread the non-reversible tree
    log-likelihood is not bitwise reproducible from one call to the next, for NONREV as for NQC.
    `computeNonrevLikelihoodBranchSIMD` adds its packets' sums under `schedule(dynamic,1)` with an
    OpenMP `reduction`, while the reversible kernel stores each packet's sum and adds them in packet
    order. The differences are rounding (1.3e-16 relative at `-T 4`, decision 038), but a bitwise
    recheck in NQC's export aborted 3 of 8 runs at `-T 4`; NQC now compares within decision 038's
    line. Would the maintainers want the non-reversible kernel to add its packets in order?
13. Upstream observation, 2026-10-29, unmodified binary: a rate model written inside
    `--model-joint` is not applied. On turtle with `-te -blfix`, `--model-joint "NONREV+G4{0.5}"`
    fits `NONREV+FO` with no rate model under `-S`, `-Q` and `-p`, exit 0, while a bare
    `"NONREV+G4"` stops `-p` with "File not found NONREV+G4"; the rate model takes effect only
    through `-m`. NQC now refuses one (decision 008). Should `--model-joint` reject a rate model, or
    apply it?
14. Upstream finding, 2026-10-29, unmodified binary: `--lie-markov` with a non-reversible model
    (NONREV on aa_example, UNREST on `example.phy`) stops at `ASSERT(0 && "this line should not be
    reached")` in `ModelMarkov::computeTransMatrixNonrev`; LG runs. NQC now refuses it (decision
    008). Should it be refused for non-reversible models, or fall back to scaling and squaring?
15. Upstream finding, 2026-10-29, unmodified binary: a one-alignment `-te` run that stops after its
    model fit and before its end-of-run dump (killed at `reportPhyloAnalysis`) resumes with "Model
    parameters restored", reports the fitted model on the input tree's branch lengths and exits 0:
    aa_example on baseline run 2's tree, `LG+G4` −7330.3269 against −7302.3888 for the run that
    finished, `NONREV+G4` −6820.0884 against −6811.2908; NQC the same. By reading, `computeInitialTree`
    saves the `-te` tree into the checkpoint before `runTreeReconstruction` restores the tree from
    it, and the fitted tree reaches the candidate set only after the fit's dump; resuming with the
    fitted tree as `-te` gives the finished run's score. Linked runs not tested. Script and outputs:
    `~/iqtree3-audit/followup-20261030T054825Z/`. Should the restore keep the fitted tree?
16. Upstream finding, 2026-10-30, unmodified binary: `--root-test` re-optimizes every branch length
    after the model fit, `-blfix` included, and the report then mixes two states. On aa_example with
    baseline run 2's tree, `NONREV+G4{0.5} -te -blfix --root-test` ends with tree length 3.4968
    against 3.0282 without the option, its root test reaching −6819.462 at the original root, while
    the `.iqtree` prints the pre-test −6825.6805 beside the re-optimized tree. With a reversible
    model on an unrooted tree (`LG+G4{0.5}`) `testRootPosition` writes no `.roottest.trees` and the
    run stops with "File not found" (exit 2). NQC refuses `--root-test` (decision 008), since its
    export would describe the re-optimized state. Script and outputs:
    `~/iqtree3-audit/rootfix-20261030T174542Z/step0a/`. Should the root test honour `-blfix` and
    leave the tree as it found it?

## Current state

S0 to S4 are complete, so gates G0 to G3 and the first release's scope are met, and the branch is
the first release candidate. S3 ran as three slices (the approved plan is frozen in
`docs/agent/archive/S3_plan_2026-10-22.md`). S4's gate ran from the clean tree at `a0ebf5cf` on
2026-10-27 (`~/iqtree3-runs/s4-gate/20261027T043226Z-a0ebf5cf/`, binary frozen in
`~/iqtree3-binaries/e31974dd/`); its record and the release audit are the S4 entry of
`CHANGELOG.md`. A second release audit (2026-10-29) found one blocker: NQC runs at more than one
thread could abort in the export diagnostic. It is fixed in `6b69c125` (decision 038) and gated
from that commit (binary `~/iqtree3-binaries/db4e7aff/`; record in the 2026-10-29 entry of
`CHANGELOG.md`). Its lower-priority findings were addressed on 2026-10-30 and gated from
`743b6eab` (binary `~/iqtree3-binaries/3db13f51/`; record in the 2026-10-30 entry of
`CHANGELOG.md`); the upstream behaviours they exposed are maintainer questions 13 to 15.

- **In the repository:** the NQC class `model/modelnonrevfixedeq.{h,cpp}`, with its 17-digit
  checkpoint record, its fit records, its linked hooks (`model/linkedmodelhooks.h`, two calls in
  `PartitionModel::optimizeLinkedModel` and its `derivativeFunk`), its export (`reportRun`:
  `<prefix>.NQC.qmat`, the first-order diagnostic and stop label, the end-of-run summary and
  `<prefix>.NQC.info`; S4), and its hooks in other shared files; the pure module
  `model/fixedeqchart.*` and its unit tests (S1); the oracle with its linked objective, reference
  climber and export functions, the drivers, the two-partition data, regression baseline, probes,
  frozen binaries and G0 manifest (S0).
- **Export (S4):** the re-import reproduces a Level 1 fit to rounding, under one alignment, `-S`
  and `-q`; the fitted Level 1 exports are labelled "first-order stationary at the optimizer's
  tolerance", which certifies no peak (decision 039); all other outputs are as S3c's binary wrote
  them. The diagnostic's evaluations (`CODE_PLAN.md` section 2.2) add to every NQC run.
- **Domain and comparisons:** default domain [1e-5, 1e5] (decision 007, `0e9fb159`); compared fits
  follow decision 039's protocol (risk 24). S3c's first part left S-3 recorded as failed in case
  B (relaxed `-S` run stopped by its round threshold 0.479 below the box's peak, not shown
  NQC-specific).
- **Cost (S3c's profiling):** time per evaluation equal to NONREV's at every size (0.97 to 1.02),
  so the code is not slower; time per fit 0.26 to 1.09 of NONREV's, mostly following where the
  round rules stop each model. All thresholds met; risks 5 and 18 retired. NQC's first linked fit
  from LG takes 5,000 to 13,000 evaluations and each later round 700 to 4,400, so a training
  run's cost is set by its first fit and its round count.
- **Continuous integration:** Actions is back on since 2026-10-22. Upstream's Build workflow is
  switched off on the fork's Actions page (`ci.yaml` unchanged), so the fork no longer checks the
  macOS, Windows or ARM builds. fixedeq's first run (39817703469, `76dd4301`) passed all jobs, as
  did the run on S4 and its audit (40550150265, `548dad20`).

## Next step

Plan and implement S5 in a new session, started in plan mode with this prompt:

> Plan and implement S5, the independent reference of gate G4 as reduced by D01. Read the
> start-of-session files, then `CODE_PLAN.md` sections 3.1 and 4 (S5) and `DECISIONS.md`, D10,
> 009 and 021 in particular. S5's list in section 4 is an outline: plan it in full first,
> setting each test's pass line before the test is written (decision 021), and ask me where
> the records leave a choice, including the small real data set of risk 9. Use the keep-awake
> helper for any run longer than a few minutes (CLAUDE.md).

The working rules are in `CODE_PLAN.md` section 4, "Implementation".
