# Plan: π-constrained non-reversible amino-acid models

> **Status, 2026-10-22: S2 complete (gate G2); the NQC class is in IQ-TREE's build; S3 not
> started; approach accepted by the IQ-TREE maintainers (reported by Peter, 2026-10-01).** Peter
> owns this document. Its
> companion, `docs/agent/CODE_PLAN.md`, holds the file-level change map, the NQC class
> specification and the test specification, and is read after this file. Decisions are
> recorded only in `docs/agent/DECISIONS.md` and are cited here by number.

## What this document is for

This is the living core document for the project, rewritten as the project moves rather than
appended to. It answers four questions for any agent picking up work:

1. what we are building and why,
2. what constraints and requirements the work must satisfy,
3. where the project currently stands,
4. what the next step is.

It is read at the start of every task, and most tasks need nothing else from the planning
documents. A task that implements or tests a slice also reads `CODE_PLAN.md`, after this file.
Keep this file current and short enough to read in full, and put history in `CHANGELOG.md`.

Labels: (verified) means read in source at `63c330d9`, or run on 2026-09-23; (reported) means
taken from a document; (proposed) means a plan choice not yet executed. Since S2 (`cd4c6823`) the
branch differs from `63c330d9` in `model/CMakeLists.txt`, `model/modelmixture.cpp`,
`model/modelprotein.{h,cpp}` and `utils/tools.{h,cpp}`. Anchors in those files refer to
`63c330d9`; at HEAD, lines after S2's insertions are offset by up to 15.

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
upstream as a pull request (`DECISIONS.md`, entry 001).

## Scientific motivation and success criteria

> Drafted from the design documents on 2026-09-23 for Peter's confirmation. The scientific
> criterion below is a placeholder until he replaces it.

The reversible workflow's frequency replacement is exact only because q_ij = R_ij π_j is
stationary at every positive π (reported). In an asymmetric factor the same replacement still
gives a valid generator but loses target stationarity (D09 item 6), and IQ-TREE's
`adaptStateFrequency` performs exactly that swap on a non-reversible matrix (verified,
`model/modelmarkov.cpp:914-930`). The constrained estimate instead refits Q inside the family
with π*Q = 0 and π* as the root distribution, so the fitted matrix is stationary at π* by
construction rather than transported there.

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
  never scores below it;
- 360 matrix parameters are counted once across linked partitions;
- π* is bitwise unchanged through construction, linking, optimization, checkpoint, restart and
  export;
- export followed by re-import reproduces the fitted likelihood;
- every legacy model in the regression baseline reproduces exactly.

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

One refinement to the synthesis's default (proposed): Level 1 uses existing flags (decision 024).
For one alignment and under `-S`, it uses `-te` with `-blfix` and braced rate parameters. For the
edge-proportional arrangement it uses `-q` with the partition rates fixed in the partition file,
because `-p` refits them every round even under `-blfix` (S0 probe (c)).

Tree search with NQC runs through ordinary IQ-TREE code but is not a first-release validation
target.

Deferred, still in the project: a small direct-flux reference in Python (D10; S5, gate G4 as
reduced by decision 028); target-rebuilt guide candidates, fixed learned candidates in
ModelFinder, the reversible/non-reversible mixing guard, an optional Level 3 outer driver, and
simulation (S6, gate G5); domain sensitivity, multiple starts, mismatch simulations and held-out
evaluation (S7, gate G6).

Out of scope: EM (D10), exact boundary fitting, branch-specific targets, non-stationary roots,
profile mixtures with distinct equilibria, analytic likelihood gradients, compiled positive-ratio
and T3 charts (decision 028), and any edit to IQ-TREE's optimizer (decision 029).

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
   checkpoint (verified, `model/modelfactory.cpp:1395-1505`; `opt_gammai` defaults to true,
   `utils/tools.cpp:7082`), and IQ-TREE writes checkpoint doubles at 10 significant digits
   (verified, `utils/checkpoint.h:22`).
7. Work inside `PartitionModel::targetFunk`'s OpenMP loop (verified,
   `model/partitionmodel.cpp:310-333`) is thread-safe: per-object workspace, no static mutable
   state, no shared warm starts, and the existing fixed-order sum preserved.
8. The coordinate domain is explicit and recorded, and an imported matrix is never clipped
   (D03). A start is not checked against the domain: as for every IQ-TREE model, the line search
   clamps trial points into it (decision 027).
9. `-m NONREV+F{...}` keeps its current non-stationary-root meaning.
10. Code is ready for an eventual upstream pull request (entry 001): small diffs to shared files
    in the local style, IQ-TREE's conventions (`outError`, 1-indexed optimizer vectors,
    `aligned_alloc` for kernel arrays), and every new file listed in `model/CMakeLists.txt`.
11. **Hard rule:** no agent contacts the IQ-TREE maintainers or posts to `iqtree/iqtree3`
    (issues, pull requests, discussions, comments, email). Peter handles all communication;
    questions for the maintainers go to the open questions below.
12. **Hard rule:** all development happens on branch `nq-constrained-pi` of the fork, and none
    of it waits on upstream.

## Where π* comes from

The target is the pooled composition of the original training data, computed once through
IQ-TREE's estimator and held fixed across cleaning treatments (D07). The estimator (verified)
counts states, runs an 8-round fixed point that distributes ambiguity codes
(`convertCountToFreq`, `alignment/alignment.cpp:4886-4933`), then calls `convfreq`
(`alignment/alignment.cpp:5882-5909`). By default `convfreq` returns at once, because
`keep_zero_freq` defaults to true (verified, `utils/tools.cpp:7262`); only under `--inc-zero-freq`
does it raise every entry below `min_state_freq` to that value and take the deficit from the
largest entry. Missing sequences count as unknown states under `-p` but
not under `-S` (verified, `model/partitionmodel.cpp:132-133`), so the partition type is part of
the record, together with the gap and ambiguity treatment, any smoothing, and the vector at 17
significant digits. Because the fixed point runs exactly 8 rounds, unknown cells shift the
result only slightly: for the turtle test data the padding changes the pooled vector by less than
the 8 decimals IQ-TREE prints, so baseline runs 8 (`-p`) and 9 (`-S`) printed the same vector. The
real π* is computed only after the confirmatory plan is in the repository, because computing it
reads the training data. S0 to S4 use the test target of decision 023, the turtle pooled
composition, which the G0 manifest records at 17 significant digits
(`test_scripts/fixedeq/manifest/g0_manifest.md`).

The oracle's port of the estimator, `test_scripts/fixedeq/oracle/target.py`, reproduces IQ-TREE's
"Mean state frequencies" line for the turtle data under both `-p` and `-S` (verified 2026-10-05).
IQ-TREE prints that line with 8 decimal places (verified, `model/partitionmodel.cpp:155-160`).

The model receives the target as `NQC+F{p1,...,p20}` in `-m` or `--model-joint`, or as
`NQC+F<NAME>` with the vector defined in an `--mdef` file (decision 006). It parses and validates
the vector itself, states any normalization, and rejects `+FO`, `+F`, `+FQ` or a missing `+F`
with a message that points to `NONREV`. With `--model-joint`, every partition receives the same
string (verified, `model/modelfactory.cpp:284-305`), and a check before linked optimization
requires every linked object to hold the same target bit for bit.

π* adds no degrees of freedom. Inference is conditional on a plug-in target, so adding 19 to a
BIC count when comparing fixed-target and free-target procedures is not justified (synthesis
section 8).

## Code plan

The file-level change map, the NQC class specification and the test specification are in
`docs/agent/CODE_PLAN.md`. This section summarizes them.

### Change map

New files in place:
- `model/fixedeqchart.{h,cpp}`: the pure module, Eigen and the standard library only (S1; in
  IQ-TREE's build since S2);
- `model/modelnonrevfixedeq.{h,cpp}`: the NQC class (S2);
- `unittest/`: the standalone googletest project (S1);
- `test_scripts/fixedeq/`: oracle, drivers, environment file, fixture generator, and from S2 the
  NQC check driver;
- `.github/workflows/fixedeq.yaml`: the fork-only workflow, with an NQC job from S2.

Shared IQ-TREE files touched, each opt-in or behaviour-preserving. Done in S2 (`cd4c6823`):
`model/CMakeLists.txt`, one `createModel` branch and a mixture-component check in
`model/modelmixture.cpp` (decisions 032, 034), a protected constructor in
`model/modelprotein.{h,cpp}` (decision 031), and the parse-only option `--nqc-domain` in
`utils/tools.{h,cpp}`. Proposed: in S3, the step option in `utils/tools.{h,cpp}` if needed, and
derivative delegation, a linked-compatibility check and a post-fit hook in
`model/partitionmodel.{h,cpp}`; in S4, the export hook and citation block in
`main/phyloanalysis.cpp`. Deferred: `main/phylotesting.cpp` (S6).
`utils/optimization.{h,cpp}` is not touched in any slice (decision 029).

### Decisions the code plan rests on

Every programming decision and the mathematical decisions D01 to D10 are in
`docs/agent/DECISIONS.md`, listed in the index at its top.

### Test strategy

Six layers, specified in `CODE_PLAN.md` section 3: the Python oracle, C++ unit tests, the
differential likelihood test, the regression baseline, invariant and nesting checks, and
sanitizer builds. Upstream's GitHub workflow gives portability checks only and never runs NQC.
The fork workflow added in S1 runs the unit tests and the oracle's tests, and from S2 builds
IQ-TREE and runs the NQC seed case; it has not yet run on GitHub, where Actions is disabled on
the fork.

### Implementation slices

Each slice maps to a gate of synthesis section 11: G0 protocol and configuration, G1
mathematical core, G2 first compiled conditional model, G3 linked and lifecycle support, G4
competing implementations, G5 native and outer workflow, G6 scientific validation.
`CODE_PLAN.md` section 4 lists each slice's files and tests in full.

| Slice | Gate | Goal | Main files | Tests written first | Done when |
|---|---|---|---|---|---|
| S0 | G0 | Evidence before any source edit | `test_scripts/fixedeq/` only; no IQ-TREE source | oracle tests; regression driver | environment built from its file; design scripts rerun from copies and compared with their reported values; oracle likelihood matches the unmodified binary; baseline recorded twice with the binary's SHA-256; S0 runtime checks answered; unmodified sanitizer run recorded; G0 run manifest written |
| S1 | G1 | Pure mathematical core | `model/fixedeqchart.*`, `unittest/`, fixtures, fork workflow | unit tests at n = 3, 4 and 20 against oracle fixtures | all pass in Release and sanitizer builds; residual maxima tabled |
| S2 | G2 | First compiled model, one alignment | class files with the class checkpoint (decision 033), `model/modelmixture.cpp`, `model/modelprotein.*`, `model/CMakeLists.txt`, domain option | root kept as for NONREV; differential at the seed; gradient step checks; target-overwrite and invalid-target rejections; re-evaluation of the accepted state; compiled optimizer against the Python port on a toy case; all under decision 030's lines | all pass; Level 1 and Level 2 fits on `example/aa_example.phy`, `+I+G4` included; regression reproduces; sanitizers clean for new code |
| S3 | G2, G3 | Linked training and lifecycle | `model/partitionmodel.*`, class checkpoint | two partitions with different compositions; `-S` and `-p`; `-T 1` against `-T 4`; restart; GTR20-seeded nesting | all pass; the linked gradient uses the scaled step |
| S4 | G3 | Export, re-import and report | `main/phyloanalysis.cpp`, class export | export then `-m FILE` likelihood equality; metadata completeness | G3 met; first release candidate |
| S5 | G4, as reduced by decision 028 | Independent reference | `test_scripts/fixedeq/` only: D10's direct balanced-flux fit in Python | agreement with NQC fits on small problems | comparison report |
| S6 | G5 | Workflow | `main/phylotesting.cpp`, outer driver | candidate-set and guard tests | G5 conditions |
| S7 | G6 | Science | none planned | the confirmatory plan's checks | G6 conditions |

## Risks and open questions

| # | Item | Resolved by |
|---|---|---|
| 1 | Final public model and option names | Peter, with the maintainers |
| 5 | Cost: 361 likelihood evaluations per one-sided gradient; realistic training runs may take hours | S3 profiling; analytic gradients stay deferred |
| 7 | Active `ASSERT`s abort a linked round that lowers the log-likelihood by more than 0.1 (verified, `model/partitionmodel.cpp:938`, `model/partitionmodelplen.cpp:135`) | Accepted 2026-10-07 by decision 029: NQC meets these `ASSERT`s as NONREV does, with no restore; S2 and S3 tests record whether any fit reaches them. S2's fits, all on one alignment, cannot reach these linked-round `ASSERT`s, and printed none |
| 8 | The default coordinate domain may not contain LG rebuilt at a skewed π* | No construction check (decision 027): coordinates outside the domain are clamped at the first line search, as for every IQ-TREE model, and bound activity is reported after the fit (decision 009); the domain option (decision 007) widens the domain. At decision 023's target the LG seed lies inside, 3.51 above the lower bound and 2.32 below the upper (G0 manifest). S2's `aa_example` fits at that target ended with 1 or 2 coordinates at the upper bound (0 to 5 during the fits) and none at the lower, read from each run's checkpoint by `fitrecords.py` |
| 9 | Real training data: a small real set for S5, and the full set with a written confirmatory plan for S7 | Peter |
| 11 | Scientific success criterion and manuscript details (placeholder above) | Peter |
| 12 | A continuity source the chart does not remove: `computeTransMatrixNonrev` switches from the eigen path to scaling-and-squaring when P's row sums deviate by more than 1e-4 (verified, `model/modelmarkov.cpp:489-500`), and the decomposition sets `nondiagonalizable` on a singular eigenvector matrix (verified, `model/modelmarkov.cpp:1330-1336`); shared with `NONREV`, frequency unmeasured | count "INFO: Switch to scaling-squaring" lines (printed under `-v`) in S2 and S3 fits. S2: none in the 183 fits of its three `aa_example` runs (`fitrecords.py`) |
| 18 | No fallback is recorded if S3 profiling shows training runs impractical while analytic gradients are out of scope | Peter, after S3 profiling |
| 19 | Baseline run 11 (`--model-joint NONREV` under `-p` at `-T 4`) does not reproduce: five repeats of the unmodified binary ended at log-likelihoods from -4975.0215 to -4974.5963, none equal to the `-T 1` run's -4974.5432, with the same topology and 442 of the model section's 447 numbers varying (largest spread 0.046). Its tree drawing varies too, so the spread rule of `CODE_PLAN.md` section 3.4 cannot judge it, and S3's `-T 1` against `-T 4` test meets the same behaviour in legacy code | Partly resolved 2026-10-03: run 11's rule and the added run 12 are decision 018, implemented in the driver on 2026-10-04. Run 12 (run 10 at `-T 4`) gave identical reports in five repeats, identical to run 10's at `-T 1` in every compared item, so at the report's printed precision the variation arises in the joint fit, not in a fixed evaluation. The S3 thread test's definition: Peter, before S3 |
| 20 | Upstream finding, sanitizer run of 2026-10-05: `ModelMarkov`'s constructor calls `setReversible` (`model/modelmarkov.cpp:72`), which reads `is_reversible` before anything sets it (line 76; UndefinedBehaviorSanitizer "load of value 190, which is not a valid value for type 'bool'"). It appears once in every run, so the IQ-TREE sanitizer runs of S2 to S4 will show it too (S1's sanitizer run covered only the unit tests; S2's showed it once in each of 55 runs, and nothing else); its effect on results was not investigated | Known before any edit; a question for the maintainers if Peter wishes to relay it |
| 21 | Upstream observation, S0 probe (g): after a tree search, the checkpoint keeps the model from before the final model optimization (`main/phyloanalysis.cpp:3895-3899` saves the search state, not the model), and a restart from it prints "Final model parameters restored". The restart itself was not tested. S2 to S4 fits run under `-te`, where no final optimization runs, and decision 025 carries the incumbent from a fixed-tree fit | Avoided by decision 025; a question for the maintainers if Peter wishes to relay it |
| 22 | Decision 025's nesting start, the fitted `GTR20+F{π*}` of probe (g) on `aa_example`, has 146 of its 360 coordinates outside the default domain (-15.94 to 11.65; `test_scripts/fixedeq/explore/thresholds.py`, 2026-10-05), because 87 of its 190 exchangeabilities sit at GTR20's own floor. Under decision 027 the start is evaluated exactly, so NQC's first log-likelihood still equals GTR20's, but the first line search clamps those coordinates, the fit then searches inside the domain only, and since decision 029 restores nothing, the final score can fall below the start. In the default domain, S2's fits from the LG seed ended with 1 or 2 coordinates at the upper bound (risk 8) | Peter, before S3: the nesting test in the default domain, in a domain widened by the domain option, or on data large enough for GTR20 to be well determined |

Resolved risks 2, 3, 4, 6, 10 and 13 to 17 were removed on 2026-10-07. Their resolutions are in
`DECISIONS.md` and `CHANGELOG.md`.

## Current state

S2 is complete (gate G2, 2026-10-22). The NQC class is in IQ-TREE's build.

Added in S2 (decisions 030 to 035):
- `model/modelnonrevfixedeq.{h,cpp}`: the NQC class with its 17-digit checkpoint, its hooks in
  `model/modelmixture.cpp`, `model/modelprotein.{h,cpp}` and `model/CMakeLists.txt`, and the
  `--nqc-domain` option (commit `cd4c6823`).
- The oracle's NQC objective, the toy problem and its calibration (`3edf321b`).
- `test_scripts/fixedeq/nqc/nqcchecks.py`, the sanitizer driver's NQC step and the CI NQC job
  (`35d4a1bf`).

The gate ran from the clean tree at `35d4a1bf` (`~/iqtree3-runs/s2-gate/20261021-35d4a1bf/`):
- Every asserted NQC check passed: root; seed, with IQ-TREE's first log-likelihood within
  1.06e-15 relative of the oracle's in five setups; gradient at the seed and near zero; target;
  21 refusals; Level 1, Level 2 `+G4` and `+I+G4` fits; and the toy fit against the oracle's port
  (difference 1.9e-6 against a line of 3.9e-4).
- Regression runs 1 to 12 passed against the baseline, the unit tests passed in both builds, and
  all 182 Python tests passed.
- Sanitizer (`sanitizer-rerun/`, rerun on 2026-10-21 after the first run was frozen by standby):
  regression runs 2 and 4, the 16 differential cases, probe (a), and the NQC root, seed, target,
  refusal and short cases all passed, with risk 20's upstream report the only finding.

Corrections to the check driver's `results.md` (decision 037), recomputed from the gate's saved
output by `test_scripts/fixedeq/nqc/fitrecords.py` (`fitrecords/` in the gate directory):
- It counts each fit twice under `-v`; the three fit runs made 6, 6 and 171 fits.
- Its `+I+G4` residual and bound rows describe the last start, not the restored model. Over all
  171 fit records the largest residuals are 2.8e-16 (‖π*Q‖∞), 6.7e-16 (row sums) and 4.4e-16
  (mean rate). The restored model, read from the checkpoint, has 1 coordinate at the upper bound
  and none at the lower, and residuals of at most 2.5e-16 as the oracle computes them from its
  coordinates.
- Its π* check reads a checkpoint key written from the class's own copy of the target. The guard
  on the root distribution is the bitwise check of `state_freq` after every decomposition, and
  every run exited normally.
- The Level 1 score IQ-TREE prints is 0.006 below the accepted score: the round loop kept fit 6,
  which gained 0.006, and printed the previous round's score, as decision 030 describes.

Earlier records stay in place: S1's module, unit tests and fixtures; S0's conda environment,
oracle, regression baseline, probes, frozen binaries and G0 manifest. Details are in
`CHANGELOG.md`.

## Next step

1. Peter's decisions before S3: how the `-T 1` against `-T 4` test is defined (risk 19), and which
   domain the nesting test runs in (risk 22).
2. Plan S3 (`CODE_PLAN.md` section 4). Its first change, S3-1, removes the shortfalls of the S2
   build that decision 037 records, in `model/modelnonrevfixedeq.cpp` and
   `test_scripts/fixedeq/nqc/nqcchecks.py`:
   - each fit's report printed once under `-v` on every path, partitions included, and a fall
     below the start shown at every verbosity (decision 029);
   - the target vectors, bound activity and residuals visible in a default run of every supported
     kind, `+I+G` included (decisions 006, 009), consistent with S4's report;
   - checkpoint keys for the accepted state's residuals, bound counts and start and optimizer
     log-likelihoods, the target key written from `state_freq`, and the driver reading these
     keys instead of the log (decision 034);
   - ModelFinder run without `-m`, which S2's refusal probably misses, refused and tested
     (decision 008), and refusal tests for `-fs` and the HMM option;
   - unlinked partitioned NQC (`-m NQC` with `-p`, `-q` or `-S`), accepted by S2 but untested,
     either tested or refused;
   - a checkpointed final score. The `accepted_score` key is the re-scored log-likelihood after
     the model's last accepted fit, under the rate parameters and branch lengths of that moment,
     so in Level 2 it is not the run's final score: IQ-TREE's later updates in the round raise
     the final score above it, by 66.9 in the short `+I+G4` run and by 0.005 in the `aa_example`
     one (`fitrecords.py`). S3's plan chooses between a key for the final score and running item
     3's comparison at Level 1 only.
   S3's tests also add: 360 parameters counted once across linked partitions, π* through linking
   and restart, linked Level 1 (decision 024), and permuted partition order.
3. S4 adds the export-time stop classification (decision 009), and compares the re-imported
   log-likelihood with a checkpointed score, not the printed one. `accepted_score` serves only a
   Level 1 fit (item 2, last S3-1 entry).
4. Peter asks GitHub Support to re-enable Actions, which GitHub disabled on the fork between
   2026-10-08 and 2026-10-19 (`CHANGELOG.md`, 2026-10-20). The agent then checks the fork workflow's first run: the
   micromamba route to the locked environment and the NQC job are unverified.

## Checklists

The mechanical wiring of a substitution model class is described in `docs/agent/ARCHITECTURE.md`
section 14. The per-slice file and test lists are in `docs/agent/CODE_PLAN.md` section 4.
