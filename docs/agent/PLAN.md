# Plan: π-constrained non-reversible amino-acid models

> **Status, 2026-11-19: S0 to S5 complete, so gates G0 to G3 and G4 as reduced by D01 are met (one
> S5 comparison failed, counted and accepted by Peter; its cause, retired risk 27, settled by
> decision 039's revision) and the branch is the first release
> candidate; SIM-1, the first simulation slice (decision 040), closed 2026-11-18; SIM-2 planned and
> approved 2026-11-19, its implementation under way; approach accepted by the IQ-TREE maintainers (reported
> by Peter, 2026-10-01).**
> Peter owns this document.

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
doi:10.64898/2025.12.01.691663, the bioRxiv preprint, resolved 2026-11-07), which takes
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
as biology (decision 040). Manuscript title and venue: to be supplied by Peter.

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
uncertainty, and non-reversibility is judged as decision 041 sets. Whatever criterion is adopted
is written as a confirmatory plan in the repository before any training or evaluation data are
read.

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

Deferred, still in the project: the simulation slices SIM-1 and SIM-2 (decision 040); target-rebuilt guide
candidates, fixed learned candidates in ModelFinder, the reversible/non-reversible mixing guard,
incumbent protection and an optional Level 3
outer driver (S6, gate G5's remaining conditions); domain sensitivity of the scientific fits (S3a
and S3b report it on test data for risk 22), multiple starts, held-out evaluation, the parametric
bootstraps of decisions 040 and 041, and gene resampling of π* as a sensitivity analysis, since
D07 holds π* fixed (S7, gate G6's remaining conditions).

Out of scope: EM (D10), exact boundary fitting, branch-specific targets, non-stationary roots,
profile mixtures with distinct equilibria, analytic likelihood gradients in the compiled model (the
oracle's S5 reference uses them, decision 042), compiled positive-ratio and T3 charts (D01), any
edit to IQ-TREE's optimizer (decision 009), simulation of root placement, of non-stationary lineages
(a limitation the manuscript states) and of alignment error with realignment and cleaning (decision
040), and a 5% likelihood-ratio test of reversibility (decision 041).

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
(`test_scripts/fixedeq/manifest/g0_manifest.md`). SIM-1 and SIM-2 simulate at NQ.bird's
stationary distribution (decision 040), never at the real π*. The oracle's port
of the estimator is checked in `CODE_PLAN.md` section 3.1.

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
validation. SIM-1 delivers G5's simulation support and begins G6's 20-state simulations; SIM-2
and S7's two bootstraps complete them. The slices run in the order decision 040 sets.

| Slice | Gate | Goal | Done when | Status |
|---|---|---|---|---|
| S0 | G0 | Evidence before any source edit | environment built from its file; design scripts rerun from copies and compared with their reported values; oracle likelihood matches the unmodified binary; baseline recorded twice with the binary's SHA-256; S0 runtime checks answered; unmodified sanitizer run recorded; G0 run manifest written | done 2026-10-05 (G0 manifest) |
| S1 | G1 | Pure mathematical core | unit tests at n = 3, 4 and 20 against oracle fixtures pass in Release and sanitizer builds; residual maxima tabled | done 2026-10-19 (`6c8d3af4`) |
| S2 | G2 | First compiled model, one alignment | root, seed, gradient, target and refusal checks and the toy comparison pass under decision 021's lines; Level 1 and Level 2 fits on `example/aa_example.phy`, `+I+G4` included; regression reproduces; sanitizers clean for new code | done 2026-10-22 (`35d4a1bf`); shortfalls in S3-1 |
| S3a | G2 | One alignment complete | S3-1 done for one alignment (fit record once under `-v`, end-of-run summary at default verbosity, full checkpoint record and restore, ModelFinder and unlinked partitions refused); GTR20-seeded nesting in the widened domain (decision 025); the skewed-target fit (023) and aa_example's risk 22 runs reported; `+R4` seed against the oracle and a Level 2 fit; incompatible restarts refused; regression reproduces; sanitizers clean for new code | done 2026-10-22 (`4749f222`) |
| S3b | G3 | Linked training and lifecycle | two partitions with different compositions under `-S` and `-p`, one target and one Q in every object; each linked fit checked and reported once (S3-1's linked half); 360 counted once; the linked gradient uses the scaled step; thread counts compared at a fixed point (decision 038); permuted order; mid-run and `--undo` restarts, incompatible ones refused; turtle's risk 22 runs and the warm starts from the box fits reported; regression reproduces; sanitizers clean for new code | done 2026-10-23 (`054914de`, reruns at `fb5bf76d`) |
| S3c | G3 | Stopping, risk 22's evidence, cost against NONREV | first part: the stopping study's asserted checks and the converged risk 22 evidence (`CODE_PLAN.md` section 4); then profiling against NONREV within its blocking threshold, or, above 3 per fit, fixed in NQC's code or accepted by Peter; any ratio above 1.25 per evaluation or 2 per fit explained (risks 5, 18, 24) | done 2026-10-26 (first part 2026-10-25, domain `0e9fb159`; profiling gate at `0b0048af`) |
| S4 | G3 | Export, re-import and report | export then `-m FILE` reproduces the log-likelihood of a Level 1 fit; metadata complete; first release candidate | done 2026-10-27 (`a0ebf5cf`); first release candidate |
| S5 | G4, as reduced by D01 | Independent reference | D10's direct balanced-flux fit in Python (decision 042) passes its known-answer validation; on four small Level 1 problems it reaches, from NQC's start and from NQC's relaxed endpoint, the peak NQC's own chart reaches, within the calibrated line, with every compared fit converged and NQC's endpoints not above that peak; failures counted; comparison report (`CODE_PLAN.md` sections 3.1, 3.3, 4) | done 2026-11-09 (gate at `13799127`); aa_example's A6 failed, labelled "W stopped short", counted and accepted by Peter (risk 27) |
| SIM-1 | G5 (simulation support), G6 (part) | Simulator check, exact recovery and screen, pipeline and pilot | AliSim's two-taxon exact check (three joint-law configurations, the root composition, frozen and gate binaries byte-identical) passes, or its fallback does; the oracle's NQC objective fitted to exact pattern frequencies returns the true Q in R1 to R6 within decision 021's SIM-1 rows; the screen's 14 fits tabled with flags and follow-ups; bird scaffold committed; pipeline with its reproducibility and wiring checks; pilot run (SIM-2's plan, first listed here, moved to its own session by Peter on 2026-11-18) | done 2026-11-18 (gate record in `CHANGELOG.md`, SIM-1 close; exact recovery passed after R4's counted failure and a solver fix) |
| SIM-2 | G6 (part) | Simulation study on Puma | the approved plan's cells (`CODE_PLAN.md` section 3.8; arms in decision 045) run, every NQC fit at stage P and stage B (decision 039), failed fits counted, every result with Monte Carlo standard errors; outputs for S7's confirmatory plan: reference distributions of the circulation statistic (decision 041), the stopping setting (decision 039), how much of the available data the fit needs; recovery at P and B for the manuscript | planned and approved 2026-11-19 (`docs/agent/archive/SIM2_plan_2026-11-19.md`); implementation next |
| S6 | G5 | Workflow | candidate-set and mixing-guard tests; G5's remaining conditions | not started |
| S7 | G6 | Science | the confirmatory plan's checks, with decision 039's convergence protocol; the parametric bootstraps of decisions 040 and 041 and the gene resampling of π*; G6's remaining conditions | not started |

## Risks and open questions

### Open risks

| # | Item | Resolved by |
|---|---|---|
| 9 | Real training data: the full set, with a written confirmatory plan, for S7 (S5 uses no real set beyond aa_example, decision 042). Expected: OrthoMaM v12 as Wheeler et al. used it (190 species, 6,769 loci, about 1,000 for training, most of 100 to 2,000 columns; decision 043) | Peter |
| 11 | Scientific success criterion and manuscript details (placeholder above) | Peter, with SIM-2's reference distributions (decision 041) |
| 23 | Reported by the S3 plan's review, not verified: resuming an `-S` tree search with `+I+G` may restore stale top-level model keys, for every model (`UPSTREAM_FINDINGS.md`, Risk 23). Outside S3's validation paths, which use `-te` | Checked before it becomes a question for the maintainers |
| 24 | Fits stop short of their peaks under IQ-TREE's stopping rules, for every model, so comparing two fits' scores (box against wide, NQC against NONREV, D08's arms) needs a convergence protocol. IQ-TREE's optimizer stays unchanged (decision 009). At scale the shortfall grows: SIM-1's 100-gene fit, capped at 20 rounds, sat at least 90 below the oracle's climb from its endpoint. The stopping rules and the measurements: `UPSTREAM_FINDINGS.md`, Risk 24 | Testing: S3c's first part (done 2026-10-25). Science: decision 039 (SIM-2's stages P and B) |
| 26 | Compute: SIM-2 and S7 run on UA HPC's Puma under UA's agentic AI policy and Peter's norms (`CLAUDE.md`, Puma). SIM-2's cost model (`CODE_PLAN.md` 3.8; estimates from SIM-1's pilot and S3c, with Puma's per-core speed, thread efficiency, evaluation and climb iteration counts assumed): 308 NQC fits with stage B and 176 GTR20 fits, about 19,700 charged core-hours in one 48 h round, rising to at least 26,300 to 39,900 over two or three rounds if the climbs run long or the estimates are 1.5 times low. S7, unverified extrapolation from S3c (`docs/agent/archive/SIM_plan_2026-11-07.md`, Appendix C): one NQC `-S` fit at 190 taxa and 1000 genes about 80 core-hours; NQC's first linked fit from LG took 5,000 to 13,000 evaluations and each later round 700 to 4,400, so a training run's cost is set by its first fit and its round count. SIM-1's pilot on the laptop (`-T 6`, round cap 20): 26 to 48 min and 0.6 GB a fit at 30 genes, 2.4 h and 1.5 GB at 100; at 1000 genes IQ-TREE asked for 10.1 GB | Peter's ceiling, 40,000 charged core-hours (2026-11-19): the agent projects the total from measured times after each round and stops for Peter only above it |
| 28 | SIM-1's 100-gene NQ.bird fit (pilot, 20-round cap) pushed a few low-flux rates toward the domain floor (5 of the chart's z coordinates within 1 of the lower bound, none at it), always underestimated, the same pairs across replicates (K→G, K→P, E→P, V→P, N→A, N→G; across the pilot's fits the largest errors sit at flux ranks 119 to 309 of 380); its unweighted log-rate RMS rose to 1.11 from 0.72 to 0.97 at 30 genes while the two highest-flux quintiles improved (RMS 0.04, 0.09) (`explore/pilot_rate_errors.py`, commit `f181cf38`). Two parts, by a planning measurement on 2026-11-19 (exploratory, SIM-2 plan archive): at 100 genes mostly IQ-TREE's process stopping short, since the oracle's climb from the capped endpoint cut the third flux quintile's pair-sum log error from 0.55 to 0.22; at 30 genes the higher point pushed the lowest-flux rates further down (0.80 to 1.35) | SIM-2 (decision 045; aim 7, stage P against stage B) |

Retired risk numbers:

| Risk | Subject | Now |
|---|---|---|
| 1, 20, 21 | | maintainer questions 1, 4 and 5 below (2026-10-22) |
| 2, 3, 4, 6, 10, 13 to 17 | | resolved by 2026-10-07 (archived decisions and `CHANGELOG.md`) |
| 5, 18 | cost; no fallback | retired 2026-10-26 by S3c's profiling: an NQC evaluation costs what a NONREV one does (`CHANGELOG.md`, S3c entry, with the labelled extrapolation to real training) |
| 7, 12 | a linked-round `ASSERT`; the switch to scaling and squaring | monitoring checks in `CODE_PLAN.md` section 4 (2026-10-22) |
| 8 | | joined 22 (2026-10-22) |
| 19 | | settled by decision 038 (2026-10-22) |
| 22 | NQC's production domain | settled 2026-10-25 by decision 007's revision, default [1e-5, 1e5] from S3c's first part's converged peaks |
| 25 | inputs for the simulation slices | settled by Peter: 2026-11-16 for SIM-1 (decisions 040, 041, 043), 2026-11-18 for SIM-2 (decisions 039, 041, 045); the evidence and the SIM-1 audit's caveats in `CHANGELOG.md`, 2026-11-18 |
| 27 | the z climb (`lbfgsb_fit`) as reference climb | retired 2026-11-16 by decision 039's revision: on aa_example it passed its z-gradient test 0.30 below the peak because a near-zero flux hid the slope, and it did not converge on S5's exact problems (`explore/risk27_probe.py`); peaks are certified by `flux.certificate`, SIM-1's by the flux reference. A z climb returns as SIM-2's stage B, labelled "highest found" unless the certificate passes (decision 039, 2026-11-19) |

### Questions for the maintainers

Peter relays these if and when he chooses; no agent contacts the maintainers (`CLAUDE.md`).
Questions 4 to 7 and 9 to 17 are stated here in brief; what was run or read, with commands,
outputs and anchors, is in `docs/agent/UPSTREAM_FINDINGS.md` under the same number.

1. Final public names for the model and its options.
2. Whether a dedicated target option, such as `--target-freq FILE`, would be preferred to
   `+F{...}` (decision 006).
3. Whether the failed-line-search path of `lnsrch` and `dfpmin` should be repaired; NQC reports a
   fall below the start and changes nothing (decision 009).
4. Upstream finding: `ModelMarkov`'s constructor reads `is_reversible` before anything sets it,
   reported by UndefinedBehaviorSanitizer once in every IQ-TREE sanitizer run, the unmodified
   binary's included; its effect on results was not investigated.
5. Upstream observation, S0 probe (g): after a tree search the checkpoint keeps the model from
   before the final model optimization. NQC avoids it by decision 025.
6. Upstream finding, S3b's sanitizer run: at `-v` and above an `-S` run ends by printing NNI
   counts from a tree it casts to the wrong type: AddressSanitizer stops IQ-TREE, and a Release
   build prints meaningless counts. Reproduced on the unmodified sanitizer binary.
7. Upstream observation, unmodified binary: NONREV's `-S` fits ended at the default 10-round cap
   in every configuration tested, on turtle still gaining and below the score NQC reaches on the
   same setup, although NQC's unbounded family lies inside NONREV's (risk 24). Cause not
   investigated.
8. NQC's report lines (the `.iqtree` "Model of substitution", a partitioned run's
   `.best_model.nex`, the `.iqtree` AliSim command) carry neither the target nor the fitted matrix,
   so none re-runs the fit; the fork points to `<prefix>.NQC.qmat` instead (Peter, 2026-10-26).
   Would the maintainers prefer the model string to carry the target, or the report a
   full-precision matrix?
9. Not verified: with `-mstore`, a linked fit that follows a partition's own update may evaluate
   stale cached transition matrices, for NONREV as for any linked model. NQC refuses `-mstore`
   (decision 008).
10. Upstream finding, S4's sanitizer run: each `-b` replicate's alignment copies a
    `pomo_sampling_method` the source alignment never set. Reproduced on the unmodified sanitizer
    binary; its effect on results was not investigated.
11. Upstream finding, 2026-10-28, unmodified binary: `--runs` with `-S` crashes, as the best run
    is restored or, with a tree search, in run 1's parsimony start tree.
12. Upstream observation, 2026-10-29 audit: at more than one thread the non-reversible tree
    log-likelihood is not bitwise reproducible from one call to the next, for NONREV as for NQC,
    because the non-reversible kernel adds its packets' sums in the order threads finish them;
    NQC compares within decision 038's line. Would the maintainers want the non-reversible kernel
    to add its packets in order?
13. Upstream observation, 2026-10-29, unmodified binary: a rate model written inside
    `--model-joint` is not applied; it takes effect only through `-m`. NQC now refuses one
    (decision 008). Should `--model-joint` reject a rate model, or apply it?
14. Upstream finding, 2026-10-29, unmodified binary: `--lie-markov` with a non-reversible model
    stops at an `ASSERT` in `ModelMarkov::computeTransMatrixNonrev`. NQC now refuses it (decision
    008). Should it be refused for non-reversible models, or fall back to scaling and squaring?
15. Upstream finding, 2026-10-29, unmodified binary: a one-alignment `-te` run that stops after its
    model fit and before its end-of-run dump resumes, reports the fitted model on the input tree's
    branch lengths and exits 0, for every model. Should the restore keep the fitted tree?
16. Upstream finding, 2026-10-30, unmodified binary: `--root-test` re-optimizes every branch length
    after the model fit, `-blfix` included, and the report then mixes two states. NQC refuses
    `--root-test` (decision 008), since its export would describe the re-optimized state. Should
    the root test honour `-blfix` and leave the tree as it found it?
17. Upstream finding, 2026-11-17, unmodified binary: `GTR20` on a rooted `-te` tree segfaults under
    `--model-joint` (`-q` and `-p`) and, as a single model, stops at an assertion in
    `PhyloTree::optimizeAllBranches`; the same tree unrooted runs. SIM-1's pipeline unroots it for
    `GTR20`. Should a reversible model unroot a rooted `-te` tree, or refuse it?

## Current state

S0 to S5 are complete, so gates G0 to G3, G4 as reduced by D01 (with aa_example's A6 failure
counted and accepted by Peter, risk 27), and the first release's scope are met, and the branch is
the first release candidate. S3 ran as three slices (the approved plan is frozen in
`docs/agent/archive/S3_plan_2026-10-22.md`). The latest IQ-TREE build gated is from `36a43711`
(binary `~/iqtree3-binaries/25d78d6e/`), the reference build for the next comparison; S5 changed
only Python test code and ran its gate at `13799127` against that binary. S4's gate record,
the three release audits and their gated fixes are the `CHANGELOG.md` entries from 2026-10-26 to
2026-10-30; the upstream behaviours they exposed are maintainer questions 10 to 16, and the
audits' low findings are in `CODE_PLAN.md` section 4, "Known and left alone".

- **Code and tests:** the inventory is `CODE_PLAN.md` sections 1.1 and 1.2.
- **Export (S4):** the re-import reproduces a Level 1 fit to rounding, under one alignment, `-S`
  and `-q`; the fitted Level 1 exports are labelled "first-order stationary at the optimizer's
  tolerance", which certifies no peak (decision 039). The diagnostic's evaluations (`CODE_PLAN.md`
  section 2.2) add to every NQC run.
- **Domain and comparisons:** default domain [1e-5, 1e5] (decision 007, `0e9fb159`); compared fits
  follow decision 039's protocol (risk 24). S3c's first part left S-3 recorded as failed in case
  B (relaxed `-S` run stopped by its round threshold 0.479 below the box's peak, not shown
  NQC-specific).
- **Cost (S3c's profiling):** an NQC evaluation costs what a NONREV one does, and time per fit
  mostly follows where the round rules stop each model; all thresholds met (`CHANGELOG.md`, S3c
  entry; training cost in risk 26).
- **Continuous integration:** Actions is back on since 2026-10-22. Upstream's Build workflow is
  switched off on the fork's Actions page (`ci.yaml` unchanged), so the fork no longer checks the
  macOS, Windows or ARM builds. fixedeq's first run (39817703469, `76dd4301`) passed all jobs, as
  did the run on S4 and its audit (40550150265, `548dad20`).
- **Independent reference (S5):** the balanced-flux reference (decision 042) passed its known-answer
  validation and, on toy, toy300 and twopart, reached the peak that the z climbs from NQC's relaxed
  endpoint reached; on aa_example the z climb met its first-order test 0.299 below the reference's
  peak (A6 failed, counted; risk 27). NQC's own endpoints lay 0.92 to 10.6 below the peak, never
  above it (gate record: `CHANGELOG.md`, S5 close). That z climb's stop was an artifact of the
  log-ratio chart, so an oracle endpoint now counts as a peak only when `flux.certificate` passes
  (decision 039; retired risk 27).
- **Plans:** simulation scope recorded (decisions 040, 041; reasons and literature in
  `docs/agent/archive/SIM_plan_2026-11-07.md`); S5's approved plan, survey and review are frozen in
  `docs/agent/archive/S5_plan_2026-11-07.md`, SIM-1's in `docs/agent/archive/SIM1_plan_2026-11-16.md`,
  SIM-2's, with its planning measurements and two reviews, in `docs/agent/archive/SIM2_plan_2026-11-19.md`.
- **Simulation (SIM-1, closed 2026-11-18, no IQ-TREE source change):** AliSim reproduces the exact
  law of an exported NQC matrix; the oracle's NQC objective recovers the true Q from exact data (R1 to
  R6); the 3-taxon screen flagged site profiles (S3) and column-context masking (S4b) at the
  certificate's activity threshold Peter set (021, 039); the bird scaffold (rooted on HUMAN) and the
  pipeline are committed with their checks passed; the pilot measured cost, the spread of Ĉ and the
  20-round cap's first effect. Ĉ is not the reversible truth's Ĉ plus the true C (0.078 for NQ.bird,
  not 0.047 + 0.058), so power is measured directly (041). Gate record and counted failures:
  `CHANGELOG.md`, SIM-1 close; the audit's corrections: its fifth-session entry and commit `f181cf38`.
  SIM-2's arms: decision 045; specification `CODE_PLAN.md` section 3.8; transport 044; pass lines 021.
- **SIM-2 (planned and approved 2026-11-19; commits 2 and 3a done):** the design, cells, methods, Puma
  execution and Peter's runbook are `CODE_PLAN.md` section 3.8's "SIM-2's plan"; the task list is section
  4. Every NQC fit runs IQ-TREE capped at 20 rounds (stage P), then the oracle's climb (stage B,
  `optimize.stage_b_climb`, "highest found" unless certified; decision 039). The drivers (`sim2.py`,
  `cloak.py`, the sbatch script) are written, reviewed and smoke-tested. In the calibration no climb
  stalled within its cap, so B's caps are still provisional (section 3.8). BAliBASE 3.0 is in place
  (CLAUDE.md, CLOAK). S6's list is written at SIM-2's close.

## Next step

SIM-2's commit 3c (`CODE_PLAN.md` section 4), CLOAK's curve. `sim/cloak.py run` started 2026-11-20 at 05:17
(laptop time) from `36fc6bc8`'s clean tree, into `~/iqtree3-runs/sim2/20261120T121706Z-36fc6bc8-cloak-run`;
its estimate is 3 to 30 h.
- If it is still running, leave it.
- If it was stopped, rerun the same command with `--out` on that directory; finished cases are
  skipped.
- When every case is done, run `cloak.py curve --run <that directory>` with the same two paths, from a
  clean tree (`CLAUDE.md`, SIM-2 commands), and commit the curve as 3c.

Then a push when Peter asks, and Peter's runbook (`CODE_PLAN.md` 3.8): preflight; test job; shards 1 to 4,
which can go before the curve; shard 5 after 3c. Before submission Peter decides two things.

1. **Stage B's stall rule and caps.** The rule's low-flux condition kept every calibration climb from
   stalling, so every B runs to its cap (2, 4, 6 and 6 h, provisional).
   - At 100 genes the lnL and Ĉ conditions alone would have stopped the climb at iteration 642, after
     which it gained another 7.0 lnL.
   - At 30 genes they would have stopped it at iteration 322, 0.06 below the endpoint.
   - The choices: keep the rule and caps; loosen the low-flux line from 0.01; or drop the condition.
   - The agent recommends keeping both. Every B records its final window movements, so round 1 shows
     whether the caps cut climbs short, and the budget guard acts after each round.
2. **The budget at 1.5 times the estimates.** With every B at its cap the projection is 26,879 charged
   core-hours at the estimates and 40,394 at 1.5 times them, above the 40,000 ceiling (risk 26). The agent
   recommends submitting round 1 as designed and deciding round 2 on the measured times, as the runbook's
   budget check already provides.

The working rules are in `CODE_PLAN.md` section 4, "Implementation".
