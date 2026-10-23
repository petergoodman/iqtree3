# Plan: π-constrained non-reversible amino-acid models

> **Status, 2026-10-22: S3a complete (gate G2's remainder): NQC fits one alignment with its full
> record and report; S3b and S3c not started; approach accepted by the IQ-TREE maintainers
> (reported by Peter, 2026-10-01).** Peter owns this document.

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

> Drafted from the design documents on 2026-09-23 for Peter's confirmation. The scientific
> criterion below is a placeholder until he replaces it.

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

Success criteria for the software (proposed). The feature works when gates G0 to G3 of the
synthesis are met with this evidence:

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
- a linked NQC fit takes time comparable to NONREV's on the same data and settings (thresholds:
  `CODE_PLAN.md` section 4, S3c's profiling).

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
| S3b | G3 | Linked training and lifecycle | two partitions with different compositions under `-S` and `-p`, one target and one Q in every object; each linked fit checked and reported once (S3-1's linked half); 360 counted once; the linked gradient uses the scaled step; thread counts compared at a fixed point (decision 038); permuted order; mid-run and `--undo` restarts, incompatible ones refused; turtle's risk 22 runs reported; regression reproduces; sanitizers clean for new code | not started |
| S3c | G3 | Cost against NONREV | profiling against NONREV within its blocking threshold, or, above 3 per fit, fixed in NQC's code or accepted by Peter; any ratio above 1.25 per evaluation or 2 per fit explained (risks 5, 18); risk 22's evidence collected | not started |
| S4 | G3 | Export, re-import and report | export then `-m FILE` reproduces the log-likelihood of a Level 1 fit; metadata complete; first release candidate | not started |
| S5 | G4, as reduced by D01 | Independent reference | D10's direct balanced-flux fit in Python agrees with NQC fits on small problems; comparison report | not started |
| S6 | G5 | Workflow | candidate-set and mixing-guard tests; G5 conditions | not started |
| S7 | G6 | Science | the confirmatory plan's checks; G6 conditions | not started |

## Risks and open questions

### Open risks

| # | Item | Resolved by |
|---|---|---|
| 5 | Cost: 361 likelihood evaluations per one-sided gradient; realistic training runs (1000 or more alignments, up to about a million columns, on HPC with dozens of CPUs) may take hours | S3c's profiling against NONREV (`CODE_PLAN.md` section 4); analytic gradients stay deferred |
| 9 | Real training data: a small real set for S5, and the full set with a written confirmatory plan for S7 | Peter |
| 11 | Scientific success criterion and manuscript details (placeholder above) | Peter |
| 18 | No fallback is chosen if S3c's profiling shows training runs impractical while analytic gradients are out of scope | Peter, if S3c's profiling crosses its blocking threshold; the candidate fallbacks are listed with it |
| 22 | NQC's production coordinate domain, a required full decision (risk 8 merged here). Decision 007's box is the per-row log-ratio counterpart of IQ-TREE's rate bounds [1e-4, 100]. Upstream fits sit on that floor: NONREV leaves 189 of 380 rates at 1e-4 on `aa_example` and 230 on turtle (checkpoints from before the final model optimization), and the GTR20 incumbent 87 of 190 (`explore/thresholds.py`). Whether the chart needs a different default is open. S2's `aa_example` fits ended with 1 or 2 coordinates at the upper bound, and the jump chains of upstream's NONREV fits (baseline runs 2 and 8), re-expressed in the chart, put 6 and 7 coordinates above it; decision 025's GTR20 start has 146 of 360 outside (-15.94 to 11.65). LG at decision 023's target lies inside (smallest ratio 3.3e-4); raising that target to powers, the start first leaves the box between frequency ratios (largest to smallest) of 306 and 963, while 6.4% of very uneven random targets (median ratio 91) put a coordinate outside (`explore/skewed_targets.py`). The nesting test and the risk 22 evidence runs use a widened domain; every other test uses the default box | Peter, before S4 is planned, from S3's bound-activity and domain-sensitivity reports; any change stays in NQC's code, and no existing inference method changes (D06) |
| 23 | Reported by the S3 plan's review, not verified: resuming an `-S` tree search with `+I+G` may restore stale top-level model keys, for every model, because per-partition refits during the search save to the top level (the model re-optimization in `IQTree::doNNISearch`) and the resumed process keeps `opt_gammai` on, so the `+I+G` starts restore from a copy of the checkpoint (`ModelFactory::optimizeParametersGammaInvar`; anchors in `CODE_PLAN.md` section 5). Outside S3's validation paths, which use `-te` | Checked before it becomes a question for the maintainers |

Retired risk numbers: 2, 3, 4, 6, 10 and 13 to 17 were resolved by 2026-10-07 (archived
decisions and `CHANGELOG.md`); on 2026-10-22, 7 (a linked-round `ASSERT`) and 12 (the switch to
scaling and squaring) became S3 monitoring checks in `CODE_PLAN.md` section 4, 8 joined 22,
1, 20 and 21 became maintainer questions 1, 4 and 5 below, and 19 was settled by decision 038.

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

## Current state

S0 to S2 and S3a are complete. S3a's gate ran from the clean tree at `4749f222` on 2026-10-22
(`~/iqtree3-runs/s3a-gate/20261022-4749f222/`); its record is the S3a entry of `CHANGELOG.md`.

- **In the repository:** the NQC class `model/modelnonrevfixedeq.{h,cpp}`, with its 17-digit
  checkpoint record, its fit records and end-of-run summary (called once from
  `reportPhyloAnalysis`), and its hooks in shared files; the pure module `model/fixedeqchart.*`
  and its unit tests (S1); the oracle, drivers, regression baseline, probes, frozen binaries and
  G0 manifest (S0).
- **One alignment is complete:** S3-1 is done for a single alignment. Its linked half, each
  linked fit checked and reported once, is S3b's (`CODE_PLAN.md` section 4).
- **S3 plan, approved 2026-10-22:** three slices, each with its own gate: S3a finished the single
  alignment, S3b adds linked training, S3c profiles against NONREV. The lists are `CODE_PLAN.md`
  section 4; the approved plan, with its survey and review findings, is frozen in
  `docs/agent/archive/S3_plan_2026-10-22.md`.
- **Risk 22 evidence so far:** S3a's `aa_example` fits in the default box and at 1e-8,1e7
  (`CHANGELOG.md`, S3a entry); S3b adds turtle's.
- **Continuous integration:** Actions is back on since 2026-10-22. Upstream's Build workflow is
  switched off on the fork's Actions page (`ci.yaml` unchanged), so the fork no longer checks the
  macOS, Windows or ARM builds. fixedeq's first run (39817703469, `76dd4301`) passed all jobs.

## Next step

1. Implement S3b in a new session, started in plan mode with this prompt (S3c likewise, after
   S3b's close):

   > Implement S3b as recorded and approved on 2026-10-22. Read the start-of-session files, then
   > `CODE_PLAN.md` sections 2.2, 3.5, 4 (its S3 introduction and S3b) and 5, and `DECISIONS.md`.
   > Where a record is unclear, read `docs/agent/archive/S3_plan_2026-10-22.md`, and ask me before
   > departing from it. Stop before any run longer than a few minutes to agree keep-awake.

   The session first checks the fork workflow's run on the push of S3a's close. Its own plan
   lists commits and gate steps only; the working rules (commits, hand-off, compute estimate and
   keep-awake before long runs, subagents) are in `CODE_PLAN.md` section 4, "Implementation".
   S3c's session plans first: its per-evaluation measurement needs Peter's approval, then the
   pilot runs.
2. Peter's full decision on the production domain (risk 22), before S4 is planned.
