# Changelog

Session history for the π-constrained non-reversible amino-acid model project, newest first.
Each entry records what was done, what failed and why, and what is next. The next session starts
from the top entry.

This is a fork of `iqtree/iqtree3`. Entries here describe work on the fork, not upstream
development.

## 2026-10-05 (fourth session, continued): S0 probes recorded, sanitizer build and timing; stop point

### Done

**Probe driver.** Wrote `test_scripts/fixedeq/probes/probes.py`, with tests in `tests/test_probes.py`
(commit `d54de36f`). Each probe's outcome, as predicted from the source, is written in the driver
before it runs.

**Recorded run.** Made from the clean tree at `fc3d3280` with the frozen binary (SHA-256
`9a72950b…`, unchanged after the run). Outputs are in
`~/iqtree3-runs/probes/record-20261005T174330Z/`: 41 IQ-TREE runs, 2724 s of compute. Probes (a)
to (f) came out as predicted; (g) differs on one comparison.

- **(a) Named target.** `+F<name>` from an `--mdef` file reaches the model exactly as `+F{...}`
  does. a1 matches baseline run 3 in every compared item, the literal and named routes are
  identical under `--model-joint` with `-p` and with `-S`, and every checkpoint `state_freq`
  equals the vector.
- **(b) `e+` inside `+F{...}`.** It exits with code 2 and "Close bracket not found in
  +F{0.08e", under `-m` and under `--model-joint`. `8e-02`, and the `e+` vector given through
  `--mdef`, both match run 3.
- **(c) Level 1.** `-te -blfix` with braced rate parameters freezes every nuisance:
  - for one alignment (`+I{0.2}+G4{0.5}`), under `-S`, and under `-q`;
  - under `-q` with each partition's rate given as `{x}` in a charpartition (run 8's rates
    0.5930, 1.5373, 0.8807, kept).

  Under `-p` (c3), the partition rates were refitted (to 0.5330, 1.6143, 0.8643) and every branch
  was multiplied by 1.148146. `-m LG+G4{0.5}` with `--model-joint` was accepted, so the
  charpartition fallback was not needed.
- **(f) Root.** The root never moved to another branch in 13 comparisons. Unrooted inputs were
  rooted at the midpoint of the longest path, per partition under `-S` and on the shared tree
  under `-p` and `-q`. Level 1 kept the root-adjacent lengths, except c3, where `-p` rescaled
  them. In every Level 2 fit the root slid along its branch until one root-adjacent length was
  about 2e-6, so the root sat at one end of its edge.
- **(e) `-S` against `-p`.** All four checked differences were observed:
  - partition rates only under `-p`;
  - the `-p` and `-S` pooled vectors differ, by at most 1.3e-9;
  - one rooted tree per partition under `-S`;
  - free parameters: 379 + 30 edges + 3 under `-p`, and 379 + 86 edges under `-S`.

  The +I+G restart, identical sequences and the "logl got worse" abort were not exercised.
- **(g) Carrying a fitted `GTR20+F{π*}` on aa_example** (π* from decision 023):
  - **Fixed-tree fit (g0).** IQ-TREE ran no final model optimization, and the checkpoint matches
    the report: every rate within the block's rounding, gamma shape 0.8875725344 against
    0.8876. Log-likelihoods relative to the checkpoint route gB (-7004.61047818093):
    - GTR20{...} string: +2.4e-11;
    - report block with exact π*: +1.1e-8;
    - the same with the report's 4-decimal gamma shape: +1.7e-7;
    - report block with its own 6-decimal frequency row: +2.0e-4.

    Ten draws perturbing gB's rates within their 10th digit spread over 1.8e-8.
  - **Search fit (gS).** The checkpoint is stale. It holds the model from before "Performs
    final model parameters optimization" (98 of 190 rates differ from the report block). Its
    model on the final tree gives -7004.0655, 0.0015 from that step's printed start of -7004.064
    and 0.283 below the final -7003.7824. Source: the save after the final optimization
    (`main/phyloanalysis.cpp:3895-3899`) writes the search state and tree but not the model
    (`tree/iqtree.cpp:134-150`, `tree/phylotree.cpp:184-194`). The model is saved only after
    the first fit and after re-fits during the search (`tree/iqtree.cpp:2388-2389, 3248-3253`).
    A restart from such a checkpoint prints "Final model parameters restored"; a restart was not
    tested.

**Sanitizer build.** Built in `~/iqtree3-build-asan` with decision 022's first settings
(`sanitize.py build`). Both steps exited 0, and no fallback was needed.

**Sanitizer timing test.** Baseline runs 2 and 4 with that binary
(`~/iqtree3-runs/sanitizer/timing-20261005T174330Z/`):
- run 2: 1165 s against 58 s in Release;
- run 4: 360 s against 17.5 s.

That is about 20 times slower. Both reproduced the baseline log-likelihoods. One finding, in
both runs: UndefinedBehaviorSanitizer, `model/modelmarkov.cpp:76`, "load of value 190, which is
not a valid value for type 'bool'". `ModelMarkov`'s constructor calls `setReversible` (line 72),
which reads `is_reversible` before it is set. No AddressSanitizer report. Recorded, not fixed.

### Failed

- **Unattended jobs stalled.** The first probe dry run and the sanitizer build were started
  unattended on 2026-10-04 at 23:45. The laptop entered Modern Standby at 23:49 and resumed at
  09:44, which paused WSL. The jobs computed for about 12 minutes over 10 hours. Diagnosed from
  the Windows Kernel-Power log and from IQ-TREE's CPU time against wall time (run c1: 183 s of
  CPU over 4 h 17 min).
- **Dry-run bugs in probe (g).**
  - Its perturbed rates were written as `np.float64(...)` text, which IQ-TREE rejects.
  - Its 10-digit route used the search fit's stale checkpoint.

  Probe (g) was revised in `fc3d3280` to fit on a fixed tree and to record the search fit as a
  separate finding.
- **One wrong prediction.** The revised (g) predicted that the report's 6-decimal block with
  exact π* would differ from gB by more than the 10th-digit noise. It differs by 1.1e-8 against a
  noise spread of 1.8e-8, so probe (g) reports "differs". The prediction was left as written.
- **Sanitizer scope.** At about 20 times slower, the planned sanitizer scope (all 12 baseline
  runs, the 16 differential cases and all probes) would take about 20 hours. Only the timing test
  was run.

### Next

1. Peter's choices, with recommendations in the stop-point report:
   - decision 024, Level 1 under `-p`;
   - decision 025, the GTR20 incumbent route and the nesting tolerance;
   - the sanitizer scope.
2. Then the sanitizer run at the chosen scope, the G0 manifest, the document updates and the
   push.

## 2026-10-04 (fourth session): plan for the rest of S0, decisions 022 and 023

### Done

Peter approved a plan for S0's remaining items: probes (a) to (c) and (e) to (g), the sanitizer
run of the unmodified code, and the G0 manifest. The plan has a stop point after the probes and
the sanitizer run, for two decisions the probes feed: Level 1 under `-p`, and the transfer of a
fitted `GTR20+F{π*}`. Peter made two choices, recorded as decisions:

- 022: the sanitizer build uses the `Mem` build type (`-g -O1`), not Debug (`-O0`), with the
  flags and fallbacks the entry states.
- 023: the S1 to S4 test target is the turtle pooled composition under `-p`.

`CODE_PLAN.md` section 3.6 now cites decision 022.

Corrected a precision statement. The `.iqtree` matrix block prints 6 decimal places, not 6
significant digits, because `precision(6)` at `main/phyloanalysis.cpp:633` acts on a stream left
in fixed mode at 1528-1529. Corrected in `CODE_PLAN.md` (section 1.2 and the probe (g) text),
PLAN.md risk 16 and `AA_MODEL_INFERENCE.md` section 14. Decision 010's "Why" says "6 significant
digits"; because entries are never edited, the correction is recorded here. Its conclusion
stands: a 6-decimal block loses even more precision than its "Why" assumed. The `.GTRPMIX.nex`
writer's "6 significant digits" in `ARCHITECTURE.md` and `FILE_INDEX.md` is correct, because it
writes to a fresh stream.

### Failed

Nothing failed.

### Next

1. The probe driver, its tests, and the recorded probe run.

## 2026-10-04 (third session, closing): oracle piece E, the estimator port; the oracle is complete

### Done

`test_scripts/fixedeq/oracle/target.py`, written from `alignment/alignment.cpp`,
`alignment/superalignment.cpp` and `model/partitionmodel.cpp`, read this session. It ports:

- the state codes and appearance sets;
- `removeGappySeq`, which keeps at least three sequences;
- `countStates`;
- `convertCountToFreq`, with its 8 fixed rounds;
- `convfreq`, which applies only under `--inc-zero-freq`;
- the linked pooling, which pads missing taxa as unknown cells except under `-S`.

Also decision 006's target parser. There are 14 new tests in `test_oracle_target.py`; all 147
fixedeq tests pass.

The port reproduces, under both `-p` and `-S`, the "Mean state frequencies" vector that baseline
runs 8 and 9 printed, for all 20 entries. This settles the question raised while planning. The
code does pad missing taxa under `-p`: partitions 2 and 3 lose podarcis and alligator when they
are built (`superalignment.cpp:495-499`), before pooling. Because the fixed point runs exactly 8
rounds, unknown cells barely move the result, so for this data the padding changes nothing at the
printed precision. PLAN.md's "Where π* comes from" now says so. It also gives the print precision
correctly: IQ-TREE prints 8 decimal places, not 8 significant digits. Decision 021 words the
estimator check as "equal at the 8 significant digits printed". The comparison applied decision
019's digits rule to the printed strings, which is equality at the precision printed (8 decimal
places).

Also updated PLAN.md: risk 17 is partly resolved by decision 021, the current state records that
the oracle is complete, and the next step is to plan S0's remaining items.

### Failed

One test failed on its first run because of an error in the test: its input frequency vector
summed to about 0.989, not 1. IQ-TREE's `convfreq` adds 1 − sum to the largest entry, so that
entry rose where the test expected it to fall. The input was corrected to a normalized vector
with one zero entry. The assertion and the code were not changed. Nothing else failed in
pieces A to E.

### Next

1. Plan, for Peter's approval, S0's remaining items: probes (a) to (c) and (e) to (g), the
   sanitizer build of the unmodified code, and the G0 run manifest.
2. Peter, before S1: the compiled code's thresholds (PLAN.md risk 17).

## 2026-10-04 (third session, continued): oracle piece D, the oracle matches the unmodified binary

### Done

Wrote `test_scripts/fixedeq/differential/differential.py` (commit `6ce9adee`). It runs 16 cases:

- matrices: NQ.pfam, and a random non-reversible Q (seed 20261004), each written at 17
  significant digits;
- rate models: none, `+G4{0.5}`, `+I{0.2}` and `+I{0.2}+G4{0.5}`;
- data: `example/aa_example.phy` and partition 1 of `turtle_aa` (sites 1 to 172, all 16 taxa);
- trees: the rooted tree files of baseline runs 4 and 8, taken from the committed
  `baseline.json`.

Each case runs `iqtree3 -s ALN -m QFILE[+...] -te TREE --show-lh -seed 1 -T 1` in its own
directory. The driver compares the printed "Initial log-likelihood" with the oracle's under
decision 021.

A dry run with `--allow-dirty` (2026-10-05 00:41 UTC,
`~/iqtree3-runs/differential/dryrun-20261005T004120Z/`) checked the driver before it was
committed. The recorded run, from the clean tree at `6ce9adee`, is in
`~/iqtree3-runs/differential/record-20261005T004141Z/`. Both used the frozen binary, and their
results are identical.

All 16 cases pass. Relative differences range from 0 to 9.8e-16, so the agreement is to rounding
error. IQ-TREE printed 17 decimals in every case, which answers S0 probe (d) and resolves PLAN.md
risk 4. The NQ.pfam case without rate variation on `aa_example` gives -7577.8547, the value
baseline run 4 reported. NQ.pfam's shipped frequencies differ from its solved stationary
distribution by 6.9e-7, and the oracle roots at the solved one, so the match is consistent with
the reading that `-m FILE` discards the file's frequency row. `CODE_PLAN.md` section 3.3 now
names the output, the tolerance and the driver.

### Failed

Nothing failed.

### Next

1. Piece E: the estimator port.

## 2026-10-04 (third session, continued): oracle piece C, the rooted likelihood

### Done

Three new modules in `test_scripts/fixedeq/oracle/`:

- `readers.py`: PHYLIP, FASTA, Newick and NEXUS charset readers.
- `gamma.py`: ports of PAML's `cmpLnGamma`, `cmpIncompleteGamma`, `cmpPointNormal` and
  `cmpPointChi2`, and of IQ-TREE's `computeRatesMean` and `computeRates`, from
  `model/rategamma.cpp` and `model/rategammainvar.cpp`, read this session. This includes how the
  constructors set the rates to 1, or to 1/(1−p) under +I+G, before the rescaling that keeps
  that sum.
- `likelihood.py`, with these conventions copied from IQ-TREE (read this session):
  - B = {N,D}, Z = {Q,E}, J = {I,L}, and `X ? - . ~ ! * U O` unknown;
  - the `-m FILE` convention: π solved from Q by column-pivoted QR, and Q scaled to mean rate 1;
  - the +I term: p times the frequency of the states every character in a column allows, p for
    an all-unknown column, and the constant-state rule of `alignment/alignment.cpp:1306-1338`;
  - category weights 1/K, (1−p)/K and 1−p;
  - branch lengths of 0 or less raised to 1e-6.

  Computed independently: transition matrices by SciPy's `expm`, and pruning with per-node
  scaling in log space.

There are 14 new tests in `test_oracle_likelihood.py`:

- pruning against brute-force enumeration;
- the probabilities of all 27 three-state patterns, with +I+G4, sum to 1;
- a reversible Q's likelihood is unchanged when the root moves;
- an ambiguous B site equals the N site plus the D site;
- the constant-state rule;
- the readers on the three test files;
- the gamma port against SciPy at shapes 0.05 to 20, with and without +I, within decision 021's
  1e-5.

All 133 fixedeq tests pass.

### Failed

Nothing failed.

### Next

1. Piece D: the oracle against the unmodified binary.
2. Piece E: the estimator port.

## 2026-10-04 (third session, continued): oracle piece B, the optimizer port

### Done

`test_scripts/fixedeq/oracle/optimize.py`, written from `utils/optimization.cpp` lines 23 and
145-939 at `63c330d9` (decision 017), read in full this session. It ports `fixBound`, `lnsrch`,
`dfpmin`, `restartParameters`, `minimizeMultiDimen` and the legacy forward-difference
`derivativeFunk`, with IQ-TREE's constants. It also ports the call
`ModelMarkov::optimizeParameters` makes: gtol = max(ε, 1e-4), no bound checks, and the score
recomputed when the returned point differs from the last evaluated one. The port keeps the C++
quirks:

- a failed line search restores x but returns f from the rejected trial, so the zero
  displacement then meets the TOLX stop;
- the inverse-Hessian update runs whenever fac² > EPS·sumdg·sumxi, so also on negative
  curvature;
- `dfpmin` stops silently after 200 iterations;
- restarts redraw every coordinate, and since IQ-TREE draws with `rand()`, the caller must
  supply the draws.

The port also records a stop reason and a failed-search count, which do not change the path.
Decision 005's scaled step is a separate function: h = 1e-4·max(1, |x|), exactly representable,
and backward at the upper bound. A central-difference reference is included for validation.
There are 13 new tests in `test_oracle_optimize.py`, one or more per branch; all 119 fixedeq
tests pass.

### Failed

Nothing failed.

### Next

1. Pieces C to E.

## 2026-10-04 (third session, continued): oracle plan approved, piece A (chart core) written

### Done

Peter approved the oracle plan, in five pieces, each committed separately:

- A, the chart core: the log-ratio and positive-ratio jump charts, references, seeds,
  diagnostics, T3 in the D04 gauge with star-gauge conversion, and the built-in matrices;
- B, the derivative rules and the port of IQ-TREE's optimizer (decision 017);
- C, the rooted likelihood, with IQ-TREE's character, gamma and +I conventions;
- D, the oracle against the unmodified binary through `--show-lh` (S0 probe (d));
- E, the frequency estimator port, checked against the "Mean state frequencies" that baseline
  runs 8 and 9 printed.

The plan's test pass lines are decision 021, recorded before any test was written. Also approved
was a correction to PLAN.md's "Where π* comes from": `convfreq` floors frequencies only under
`--inc-zero-freq`, because `keep_zero_freq` defaults to true (verified,
`alignment/alignment.cpp:5884-5886`, `utils/tools.cpp:7262`). While planning, two facts came up
for pieces C to E. With `-m FILE`, IQ-TREE discards the file's frequency row and solves π from Q;
this was read by a search agent and is to be confirmed by piece D. Baseline runs 8 (`-p`) and 9
(`-S`) printed identical pooled frequencies, although `model/partitionmodel.cpp:132-133` pads
missing taxa as unknown under `-p` only; piece E settles this.

Piece A: `test_scripts/fixedeq/oracle/` with `cases.py` (balanced fluxes from directed cycles),
`chart.py`, `t3.py` and `builtin.py`, and the tests `test_oracle_chart.py`, `test_oracle_t3.py`
and `test_oracle_documented.py`. The jump chart is written from P-log section 3 and synthesis
section 3, and T3 from P-log section 4 and its Appendix B code. ν is solved by column-pivoted QR.
All 106 fixedeq tests pass in `iqtree3-fixedeq`, 64 of them new. Under decision 021, the three
hard cases are printed rather than asserted. A skewed target with minimum 1e-4, two state blocks
joined by logits of -11.5, and coordinates spread over [-11.5, 2.3] all gave residuals between
6.6e-17 and 6.2e-16.

### Failed

Nothing failed.

### Next

1. Pieces B to E.

## 2026-10-04 (third session, continued): rerun accepted as the reproduced record

### Done

Explained to Peter in chat what the design scripts are, where the 133 values come from, why they
were rerun, and what the 24 misses are; his question showed the earlier report had assumed that
background. Peter then accepted the rerun as the reproduced record, recorded as decision 020: the
109 reproduced values count as verified on this platform, and the 24 misses stay "reported",
platform-sensitive, and never pass/fail test targets. PLAN.md risk 10 is resolved, and
`CODE_PLAN.md` section 3.1 cites decision 020.

### Failed

Nothing failed.

### Next

1. Plan the oracle (`CODE_PLAN.md` section 3.1) in plan mode for Peter's approval.

## 2026-10-04 (third session, continued): design scripts rerun, 24 of 133 values not reproduced

### Done

Recorded decision 019 and wrote `test_scripts/fixedeq/design_rerun.py` with
`tests/test_design_rerun.py` (commit `398449a6`; all 42 fixedeq tests pass). The driver holds 133
values transcribed from P-positive section 12 (81 values, lines 299 to 316) and synthesis
sections 10.2 and 10.3 (52 values, lines 362 to 385). Each value has a class under decision 019:
36 exact, 71 digits, 24 zero and 2 bound. The test file, not in the approved plan, was added
because the digit counting and class boundaries are numerical edge cases.

Ran the driver at 23:35 UTC from the clean tree at `398449a6` into
`~/iqtree3-runs/design-rerun/20261004T233525Z/`, with Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0.
These are the versions synthesis section 10.1 states for its rerun; P-positive does not state its
versions. Every copy's SHA-256 matched the design README. `model/modelprotein.cpp` was taken from
`63c330d9` (SHA-256 `ad2b3426…`). All three scripts exited 0 (2.1 s, 0.1 s, 18.9 s), and
`verify_constrained_nq.py` reported no hard failures. The full table is in `comparison.md` in the
run directory.

Result: 109 reproduced and 24 not. Every zero-class value and every stop label reproduced, and
every exact value except one evaluation count. Not reproduced (rerun values rounded here to one
digit more than printed; full values in `comparison.json`):

| Script | Document:line | Quantity | Printed | Rerun |
|---|---|---|---|---|
| verify | P-positive:304 | `2.c2_dnu_formula_error` | 1.1e-10 | 4.826e-11 |
| verify | P-positive:308 | signed coordinate, FD error at h = 1e-4 | 4.7e-4 | 4.649e-4 |
| verify | P-positive:308 | signed coordinate, FD error at h = 1e-6 | 0.070 | 0.07524 |
| verify | P-positive:308 | positive ratio, relative FD error at u = 1e-2 | 2.1e-6 | 1.157e-6 |
| verify | P-positive:312 | `10.frechet_vs_central_difference` | 1.2e-10 | 1.256e-10 |
| verify | P-positive:313 | cycle coordinate, correct target | -1.2e-6 | -1.092e-7 |
| check | P-positive:316 | NQ.PFAM row sums, maximum (bound) | ≤ 2e-6 | 2.0000000000245516e-6 |
| chart | synthesis:362 | E1 positive chart, median relative error | 8.02e-5 | 8.120e-5 |
| chart | synthesis:362 | E1 positive chart, maximum | 0.00688 | 0.007035 |
| chart | synthesis:362 | E1 log chart, maximum | 0.01969 | 0.019703 |
| chart | synthesis:366 | sweep, logit 1e-2, log chart | 3.24e-6 | 3.949e-6 |
| chart | synthesis:367 | sweep, logit 1e-4, log chart | 1.19e-4 | 6.112e-4 |
| chart | synthesis:367 | sweep, logit 1e-4, positive chart | 7.72e-5 | 7.711e-5 |
| chart | synthesis:368 | sweep, logit 1e-6, log chart | 9.55e-2 | 2.757e-2 |
| chart | synthesis:369 | sweep, logit 1e-8, log chart | 6.39 | 0.2309 |
| chart | synthesis:370 | sweep, logit 0, positive chart | 7.71e-5 | 7.720e-5 |
| chart | synthesis:382 | E2 instance 0, positive, legacy: logL | -11343.902245 | -11343.9025434 |
| chart | synthesis:382 | E2 instance 0, log, legacy: logL | -11303.734515 | -11303.7345135 |
| chart | synthesis:383 | E2 instance 1, positive, legacy: logL | -12419.902818 | -12419.9003954 |
| chart | synthesis:383 | E2 instance 1, positive, legacy: evaluations | 700 | 641 |
| chart | synthesis:383 | E2 instance 1, log, legacy: logL | -12403.624673 | -12403.6246880 |
| chart | synthesis:384 | E2 instance 2, positive, legacy: logL | -14544.987015 | -14544.9869813 |
| chart | synthesis:384 | E2 instance 2, log, scaled: logL | -14543.503265 | -14543.5032655 |
| chart | synthesis:385 | E2 instance 3, positive, legacy: logL | -12552.097580 | -12552.0936059 |

By kind, 14 are finite-difference errors, 9 are optimizer end points or path counts, and 1 is the
bound. The rerun's 641 evaluations for E2 instance 1 equal the "D-review" count that synthesis
line 387 quotes against its own rerun's 700. As Peter instructed, the mismatches were reported to
him and not investigated, and the oracle plan was not proposed.

### Failed

The rerun did not reproduce 24 of the 133 documented values (above). The procedure itself did
not fail: every script ran, and no rule or table entry was changed after the results were seen.

### Next

1. Peter: decide how to treat the 24 values not reproduced (PLAN.md risk 10).
2. Then propose the oracle plan (`CODE_PLAN.md` section 3.1) for Peter's approval.

## 2026-10-04 (third session): decision 018 in the regression driver, run 12 recorded

### Done

Peter approved a plan for this session with three choices: run 11 is judged literally on the five
parts decision 018 names; the design-script rerun is a kept driver in the repository; and rerun
values are matched in four classes (to be recorded as decision 019).

Implemented decision 018 in `test_scripts/fixedeq/regression/regress.py` (commit `f9bb0e96`).
`compare` judges run 11 only on its exit code, its REFERENCES and SEQUENCE ALIGNMENT exactly, and
its SUBSTITUTION PROCESS and tree file with the numbers removed. A section added to or missing
from a later run 11 report is not checked, because run 8, the same command at `-T 1`, is compared
in full. Run 12 (run 10 at `-T 4`) joined the run list. `baseline --runs` writes a baseline of
selected runs, and `compare` accepts several `--baseline` files and refuses a run found in two.
There are nine new tests in `test_regress.py`, and all 21 pass in `iqtree3-fixedeq`. On the five
existing extracts, `compare` failed run 11 in each before the change (its tree drawing is
unstable) and passed every run after it.

Recorded run 12 from the frozen binary (SHA-256 `9a72950b…c04fb501`, equal to
`~/iqtree3-build/iqtree3`): five repeats, 23:19 to 23:31 UTC, into
`~/iqtree3-runs/baseline/run12/rep1` to `rep5`, each running run 8 first to supply run 12's Q
matrix and tree. All exited 0. Run 12's baseline is in
`test_scripts/fixedeq/regression/baseline/run12/`, and the existing `baseline.json` and
`baseline.md` are unchanged. Results, from the driver's output:

- Run 12 was identical in every compared item across all five repeats: log-likelihood
  -4974.5436, 33 free parameters, tree length 0.7132.
- Run 8, rerun five more times, passed `compare` against the original baseline each time
  (-4974.5432).
- Run 12's report at `-T 4` is identical to run 10's at `-T 1` in every compared item. This was
  checked by an ad hoc read-only script in the session scratchpad, run in the project
  environment. At the report's printed precision (log-likelihood to 4 decimals), a fixed
  evaluation does not depend on the thread count, and run 11's variation arises in the joint fit.
  PLAN.md risk 19 is updated.

`CODE_PLAN.md` section 3.4 and the `CLAUDE.md` regression-driver bullet now describe the second
baseline file and the two-file `compare`.

### Failed

Nothing failed.

### Next

1. Record decision 019, then rerun the three design scripts with the kept driver and compare
   their outputs with the documents' values (PLAN.md risk 10).
2. Propose the oracle plan (`CODE_PLAN.md` section 3.1) for Peter's approval.

## 2026-10-03 (second session, closing): decision 018 and handoff to a new session

### Done

Peter accepted both recommendations from the baseline report, recorded as decision 018: run 11
is judged only on its exit code, its REFERENCES and SEQUENCE ALIGNMENT sections, and its
SUBSTITUTION PROCESS section and tree file with every number removed (its tree section is not
compared, because the ASCII drawing changes with branch lengths); and a run 12, run 10 at
`-T 4`, joins the baseline list, to be recorded five times from the frozen binary. PLAN.md risk
19 is now partly resolved; the definition of S3's thread test stays open. `CODE_PLAN.md`
sections 3.4 and 3.5 updated.

Audit of the record, at Peter's request: all work through `6c90ff25` was committed and pushed,
every piece was documented where the document roles put it, and the memory files were current.
The one gap for a new agent was that `CLAUDE.md` said nothing about the test setup; it now
states how to run Python in the `iqtree3-fixedeq` environment, how the regression driver and the
frozen binary are used, why the driver needs `--git git.exe`, and that `test_scripts/fixedeq/`
is this project's directory.

### Failed

Nothing failed.

### Next

1. Implement decision 018 in `regress.py`, with tests, and record run 12 from the frozen binary
   without overwriting the existing baseline.
2. Rerun the three design scripts, then plan and write the oracle (PLAN.md, next step).

## 2026-10-03 (second session, continued): S0 regression driver and baseline recorded

### Done

Wrote `test_scripts/fixedeq/regression/regress.py` (commands `freeze`, `run`, `baseline`,
`compare`) and `test_scripts/fixedeq/tests/test_regress.py` (12 tests, all passing in the
`iqtree3-fixedeq` environment), committed as `086fef2f`. The driver runs each IQ-TREE command in
a new directory, compares every report section except the header, "ALISIM COMMAND" and "TIME
STAMP", plus the tree file and exit code, and records provenance: driver commit (from Windows
`git.exe`, since Linux git sees 946 false changes from line endings), binary SHA-256, compiler
and flags, package versions, system, environment variables and input checksums.

A dry run into WSL `/tmp/regress-dryrun` (runs 6 and 7 twice, runs 8 and 10 once) checked the
driver end to end before recording: runs 6 and 7 matched between repeats, `compare` passed them,
run 8 reproduced the smoke test's -4974.5432 and run 10, reading `run08.Q.txt` by its relative
name, reproduced -4974.5436.

Pre-flight: upstream `master` still `63c330d9`; working tree clean at `086fef2f`; a rebuild of
`~/iqtree3-build` compiled and linked nothing and left the binary's SHA-256 at `9a72950b…`.
Froze it as `~/iqtree3-baseline/iqtree3` with `iqtree3.sha256`. Recorded rep1 and rep2 of the 11
runs and three more repeats of run 11 under `~/iqtree3-runs/baseline/` (log in `driver.log`),
06:12 to 06:31 UTC, and wrote `test_scripts/fixedeq/regression/baseline/baseline.json` and
`baseline.md`. All runs exited 0. Results, from `baseline.md`:

| Run | Model | Log-likelihood | Free parameters | Identical across repeats |
|---|---|---|---|---|
| 1 | LG+G4 | -7301.8799 | 34 | yes |
| 2 | NONREV | -6999.7817 | 413 | yes |
| 3 | NONREV+F{...} | -7018.7227 | 413 | yes |
| 4 | NQ.pfam | -7577.8547 | 34 | yes |
| 5 | GTR20 | -7088.9631 | 241 | yes |
| 6 | UNREST (DNA) | -22669.6869 | 43 | yes |
| 7 | 12.12 (DNA) | -22669.6808 | 43 | yes |
| 8 | joint NONREV, `-p` | -4974.5432 | 412 | yes |
| 9 | joint NONREV, `-S` | -4900.5282 | 465 | yes |
| 10 | run 8's Q re-imported | -4974.5436 | 33 | yes |
| 11 | run 8 at `-T 4`, 5 repeats | -4975.0215 to -4974.5963 | 412 | no |

Run 11 does not reproduce, and none of its repeats equals run 8. The topology is the same in all
five repeats, but 442 of the 447 numbers in its model section vary (largest spread 0.046), and
its ASCII tree drawing changes with the branch lengths, so the tree section is marked
"unstable". Recorded as PLAN.md risk 19. `CODE_PLAN.md` section 3.4 now states what is compared,
the relative names in run 10 and the five repeats of run 11.

### Failed

Nothing failed. One ad hoc inspection of `baseline.json` (the run 11 spreads and the tree-section
difference) was run with the Windows system Python, read-only, instead of the project environment.

### Next

1. Peter: a comparison rule for run 11 before S2, and a definition of the S3 thread test (risk 19).
2. Rerun the three design scripts and write the oracle (PLAN.md, next step).

## 2026-10-03 (second session): root and optimizer-port decisions, S0 probes, documents committed

### Done

Explained the plan, its open items, and S0 in plain terms to Peter in chat (not recorded in the
documents). Peter then made three calls, recorded as follows:

- Root policy for test runs: IQ-TREE roots an unrooted tree itself, and the root edge is not
  searched. Recorded as decision 016, which states the four root rules of synthesis section 7.5.
  In Level 2 fits the two branch lengths beside the root follow IQ-TREE's ordinary branch-length
  optimization; that reading of "not moved" is the agent's, for Peter to revise.
- Two S0 probes added to `CODE_PLAN.md` section 4: (f) whether the root moves under `-te`; (g)
  how precisely a fitted `GTR20+F{π*}` matrix carries into a new run.
- The oracle's optimizer copy matches IQ-TREE: decision 017, a port written from
  `utils/optimization.cpp`, which states that floating-point identity across languages is not
  expected and that the toy-comparison tolerance is fixed before the test is written.

Also added to `PLAN.md`: risks 15 to 18 (root policy resolved for S0 to S4; GTR20 incumbent
transfer; provisional thresholds to be made final before S1; no fallback for an impractical
runtime), and the rule that the real π* is computed only after the confirmatory plan exists.
`CODE_PLAN.md` gained a `.gitignore` in the file map, the optimizer port in layer 1, a note on the
nesting test, and one anchor row.

Read-only checks: upstream `master` is `63c330d9` (`git ls-remote`), so no sync is needed; no
`iqtree3-fixedeq` environment exists in WSL; conda warns that adding `defaults` implicitly is
deprecated; the root `.gitignore` has no Python cache rule; the root anchors in decision 016 were
read.

Commit `049d02ee` holds the three earlier sessions' documentation work, unchanged; `87eb5115`
holds the decisions and plan edits above. Both were pushed to the fork.

Created the conda environment `iqtree3-fixedeq` in WSL (decision 013) from the new
`test_scripts/fixedeq/environment.yml`: conda-forge only, with Python 3.12, NumPy 2.3.5 and SciPy
1.17.0 (the versions the design synthesis used for its reruns), and pytest 8.4. Conda 25.7.0
resolved 45 packages, all from conda-forge, and the environment reports Python 3.12.14, NumPy
2.3.5, SciPy 1.17.0 and pytest 8.4.2. The exact package list with checksums is in
`test_scripts/fixedeq/environment.lock.txt`. Added the directory's own `.gitattributes` (decision
014) and a `.gitignore` for Python caches. The environment lives at
`~/anaconda3/envs/iqtree3-fixedeq`, outside the repository. Conda still prints a notice about the
`defaults` channel and a deprecation warning about adding it implicitly; both come from the base
installation's configuration, which was not changed.

### Failed

Nothing failed.

### Next

1. Confirm upstream `master` is unchanged, then write the regression driver and record the
   baseline (PLAN.md, next step).

## 2026-10-03: decision 015 replaces 004; `ModelMarkov` is no longer edited

### Done

Peter chose to build Q in the NQC class's `decomposeRateMatrixNonrev` override, write the 380
rates into `rates[]`, and call the unchanged `ModelMarkov::decomposeRateMatrixNonrev`, which was
decision 004's recorded fallback. Recorded as decision 015, which supersedes 004, after comparing
three placements of the z-to-Q build in source:

- in `getVariables` with the unchanged base decomposition (the maintainers' description and the
  `ModelLieMarkov` pattern): rejected, because the base decomposition then trusts `rates[]`, and
  `ModelProtein::restoreCheckpoint` writes `rates[]` at 10 significant digits and decomposes
  (`model/modelprotein.cpp:1265-1276`) on the default +I+G path (`model/modelfactory.cpp:1504,
  1550`);
- in the decomposition override with 004's helper split: rejected, because it edits the routine
  every non-reversible model uses, and its only reason (`-optfromgiven`) is covered by
  decision 008;
- in the decomposition override with the unchanged base: adopted. Every decomposition passes
  through the override, and a bitwise check of `state_freq` after the base call catches any later
  change to the base routine's π guards.

Exact form of the change 004 would have made, recorded for reference: five lines inserted between
`model/modelmarkov.cpp:1282` and `1284` (call the helper, close the function, open the helper,
declare `int i, j;`), so that lines 1284-1386 became the helper's body unchanged, plus one
protected declaration in `model/modelmarkov.h`.

Files changed: `docs/agent/DECISIONS.md` (entry 015; entry 004's status line now "superseded by
015"; preamble), `docs/agent/PLAN.md` (status, change map, decision list, S2 row, risk 13
resolved, current state), `docs/agent/CODE_PLAN.md` (`ModelMarkov` row removed and listed as not
touched, decomposition and checkpoint paragraphs, S2 file list, three anchors), `AI_DISCLOSURE.md`.
No source code was touched; IQ-TREE was neither built nor run.

Also given in chat, not recorded in the documents: an explanation of C++ inheritance and virtual
overriding as used by IQ-TREE's model classes, and a step-by-step walkthrough of the nQ inference
call path.

### Failed

An earlier chat reply in this exchange treated 004's recorded fallback as the same thing as the
`getVariables` placement, and recommended that placement. Tracing the checkpoint restore path
showed they differ; the recommendation was corrected before Peter decided.

### Next

1. S0: write the regression driver and record the baseline from the unmodified binary (PLAN.md,
   next step).
2. Peter, optionally: tell the maintainers that Q is built in the decomposition override rather
   than in `getVariables`, and why (decision 015).

## 2026-10-01: maintainer meeting notes mapped to the code

### Done

Peter reported a meeting with the IQ-TREE maintainers: they accepted the plan's method (log-ratio
jump-chain coordinates converted to a valid Q that has π* as its equilibrium by construction) and
stressed that IQ-TREE's BFGS needs a continuous objective. They named functions and gave advice,
which Peter wrote down; this session mapped each to the code, read in source at HEAD `4c5f061f`:

- Transition probability matrix: `ModelMarkov::computeTransMatrix`, which calls
  `computeTransMatrixNonrev` for non-reversible models (`model/modelmarkov.cpp:463-512`), computes
  P(t) from the cached eigensystem or by scaling and squaring.
- Optimize parameters, at three levels: `ModelFactory::optimizeParameters`
  (`model/modelfactory.cpp:1570`) alternates branch lengths with `optimizeParametersOnly`, which
  calls the model's and the rate model's `optimizeParameters` (1275-1315);
  `ModelMarkov::optimizeParameters` (`model/modelmarkov.cpp:1166-1244`) fits Q for one
  alignment; `PartitionModel::optimizeLinkedModel` (`model/partitionmodel.cpp:753-840`) fits a
  linked Q and never calls the model's `optimizeParameters`.
- `getNDim` (`model/modelmarkov.cpp:964-976`) returns `num_params` for a non-reversible model and 0
  when fixed.
- `setVariables` copies model state into the 1-indexed optimizer vector and `getVariables` copies
  it back, returning whether anything changed (`model/modelmarkov.cpp:1018-1089`); the "+1" in
  `memcpy(variables+1, ...)` is the 1-indexing.
- `minimizeMultiDimen` (`utils/optimization.cpp:750-778`) wraps `dfpmin`, which takes gradients
  from the virtual `derivativeFunk` and steps with `lnsrch`, clamping trial points by `fixBound`
  (149-156, 682-686).
- Every function above is virtual except `minimizeMultiDimen`; `dfpmin` and `lnsrch` are private
  and non-virtual (`utils/optimization.h:129-176, 229-232`; `model/modelmarkov.h:189-469`).

Comparison with the plan, reported to Peter:

- Matches: a new class (decision 002), 360 dimensions with no frequency degrees of freedom,
  1-indexed vectors, `minimizeMultiDimen` unchanged (the scaled step enters through the virtual
  `derivativeFunk`, decision 005), transition matrices reused (D06), minimal shared-file changes
  (decision 001).
- Differs in placement only: the maintainers described the Q-to-z conversion inside
  `setVariables`; the plan stores z as authoritative state (decision 003), encodes the start once,
  and builds Q from z in the decomposition override. `ModelLieMarkov` is the in-tree precedent for
  the plan's layout (`model/modelliemarkov.cpp:899-933`).
- Decision 004 rests on `-optfromgiven` overwriting π* in the unchanged base decomposition, but
  decision 008 rejects `-optfromgiven` for NQC and the flag is set only from the command line
  (`utils/tools.cpp:3316`, default false at 7284). With `FREQ_USER_DEFINED` the base skips the π
  reset and solve (`model/modelmarkov.cpp:1252, 1266`), so 004's recorded fallback needs no
  `ModelMarkov` edit; recorded as PLAN.md risk 13 for Peter.
- With `num_params = 360` and `FREQ_USER_DEFINED`, the base `getNDim` and `getNDimFreq` already
  return 360 (0 when fixed) and 0, so the overrides listed in CODE_PLAN.md section 2.2 may be
  unnecessary; to settle in S2.
- S5's recording-only edit to `dfpmin` and `lnsrch` (decision 009) is in tension with the advice
  to leave the optimizer unchanged; PLAN.md risk 14.

Continuity under BFGS: the log-ratio chart is smooth on all of R^360 and every point gives
π*Q = 0 and unit mean rate exactly (P-log section 3.3; P-positive Theorem 10, read this session),
and the fixed reference destinations (D01, decision 011) keep coordinates from jumping. NQC never
reaches the `min_state_freq` guard that returns 1e30 inside `ModelMarkov::targetFunk`
(`model/modelmarkov.cpp:1103-1108`), which `NONREV` can hit because it re-solves π at every
evaluation. Box clamping in `lnsrch` is shared with every IQ-TREE model and handled by D03. One
source is not in the plan: `computeTransMatrixNonrev` switches from the eigen path to scaling and
squaring when P's row sums deviate by more than 1e-4 (`model/modelmarkov.cpp:489-500`), and the
decomposition sets `nondiagonalizable` on a singular eigenvector matrix (1330-1336). It is shared
with `NONREV`, its frequency is unmeasured, and it prints "INFO: Switch to scaling-squaring" under
`-v` (`utils/tools.cpp:1136-1137`); PLAN.md risk 12.

Edited `docs/agent/PLAN.md` (status block, current state, risks 12 to 14) and `AI_DISCLOSURE.md`.
No decision was changed and no source code was touched. IQ-TREE was neither built nor run.

### Failed

Nothing was attempted that failed. Peter's notes end mid-sentence at "make as little"; it was read
as "as little change to existing code as possible", which is unverified.

### Next

1. Peter: decide whether decision 004 is superseded by its fallback (PLAN.md risk 13), and, if
   useful, confirm with the maintainers that storing z rather than deriving it in `setVariables`
   is acceptable.
2. S0: write the regression driver and record the baseline from the unmodified binary (PLAN.md,
   next step).

## 2026-09-23 (second session): code plan written and approved

### Done

Read the planning, design and reference documents in the order the task prompt set, then read in
source every file the plan touches and verified each anchor the plan cites at HEAD `4c5f061f`,
which matches `63c330d9` outside documentation and configuration. Wrote the code plan, which
Peter approved:

- `docs/agent/PLAN.md`: every section filled, covering the scope of the first release (gates G0
  to G3, slices S0 to S4), constraints, where π* comes from, a summary of the code plan, the test
  strategy, slices S0 to S7, risks and open questions, current state and next step. The
  scientific motivation and success criteria are drafted for Peter's confirmation, and the
  scientific criterion is a placeholder.
- `docs/agent/CODE_PLAN.md` (new, at Peter's request): the companion to PLAN.md, holding the
  file-level change map, the NQC class specification, the test specification, the per-slice
  lists and the source-anchor table. It is subordinate to PLAN.md, is never read without it, and
  records no decisions.
- `docs/agent/DECISIONS.md`: programming decisions 002 to 014, and notes resolving three of the
  pending candidates.
- `CLAUDE.md`: the hard rule that agents never contact the IQ-TREE maintainers, the companion
  document in the reading order, and the opening line.

Peter's answers this session: the provisional model name is NQC; agents never contact the IQ-TREE
maintainers, and development on the fork never waits on upstream; CODE_PLAN.md is added, PLAN.md
stays the full-scope plan that points to it, and neither overlaps DECISIONS.md.

Findings from the source that the design documents lack, each read in source:

- A third optimizer entry point: `-jointopt` runs `ModelFactory::optimizeAllParameters` with the
  model's bounds hard-coded to [MIN_RATE, MAX_RATE] (`model/modelfactory.cpp:1331-1377`), so NQC
  must reject it.
- The default +I+G optimization saves and restores the model through its checkpoint
  (`model/modelfactory.cpp:1395-1505`; `opt_gammai` defaults to true, `utils/tools.cpp:7082`),
  and checkpoints store doubles at 10 significant digits (`utils/checkpoint.h:22`), so NQC stores
  its own state at 17 digits.
- `ModelProtein::computeTipLikelihood` treats B, Z and J as two-state ambiguities
  (`model/modelprotein.cpp:1350-1366`), the `ModelSubst` default treats them as unknown
  (`model/modelsubst.cpp:199-209`), and tip partials come from this function
  (`tree/phylotreesse.cpp:366-373`), so a class outside `ModelProtein` would change the
  likelihood of such data.
- `ModelMarkov::getModelByName` receives no `ModelsBlock` (`model/modelmarkov.cpp:1954-1968`);
  `createModel` has one.
- A `+` or `*` inside a number splits a `+F{...}` model string (`model/modelfactory.cpp:259,
  296`); not yet run.
- `ModelMarkov::readStateFreq(string)` accepts NaN and hands non-numeric tokens to AliSim's
  random-distribution parser (`model/modelmarkov.cpp:1761-1796`, `utils/tools.cpp:408-428`).
- `freqTypeString` returns an empty string for protein `FREQ_USER_DEFINED`
  (`model/modelmarkov.cpp:204-206`), which explains the dropped `+F{...}` in printed model names
  observed in the first session of the day.
- `searchGAMMAInvarByRestarting` (`main/phyloanalysis.cpp:4384`) has no callers.
- `--show-lh` fixes branch lengths, keeps identical sequences, ignores checkpoints and turns on
  debug output (`utils/tools.cpp:3428-3437`); `-blfix` also switches off the +I+G restart path
  (`utils/tools.cpp:3389-3395`).
- Under `-p`, partition rates are re-optimized in every round
  (`model/partitionmodelplen.cpp:150-157`), which matters for a fixed-nuisance fit.
- Active `ASSERT`s abort a linked round whose log-likelihood falls by more than 0.1
  (`model/partitionmodel.cpp:938`, `model/partitionmodelplen.cpp:135`).

Checks run, all read-only:

- `git rev-parse HEAD 63c330d9`, `git diff --stat 63c330d9 HEAD` and
  `git config --get core.autocrlf`: HEAD differs from the baseline only in documentation,
  `.claude/settings.json` and `.gitignore`; `core.autocrlf=true`; there is no `.gitattributes`.
- `MSYS_NO_PATHCONV=1 wsl.exe -d Debian -- bash -lc "ls -d ~/miniforge3 ~/miniconda3 ~/anaconda3
  ~/mambaforge; command -v conda mamba micromamba python3 git; python3 --version; ls
  ~/iqtree3-build/iqtree3 ~/iqtree3-smoke; nproc; free -g"`: `~/anaconda3` exists; `conda`,
  `mamba` and `micromamba` are not on the PATH of such a shell; the system `python3` is 3.11.2;
  the built binary and the smoke-test outputs exist; 12 processors; 6 GB of memory.
- `MSYS_NO_PATHCONV=1 wsl.exe -d Debian -- bash -lc "~/anaconda3/bin/conda --version;
  ~/anaconda3/bin/conda env list; ~/anaconda3/bin/conda config --show solver channels"`:
  environments `base`, `geodiff` and `rdkit_env`; solver libmamba; channel `defaults`. The version
  command printed a truncated path instead of a version number; not investigated.

IQ-TREE was neither built nor run in this session.

### Corrections to reference documents

Found while planning and reported to Peter; the first five were applied later in the session, as
the last paragraph of this section records:

- `AA_MODEL_INFERENCE.md` section 4.6 names an option `--opt-model-rate-joint`, which does not
  exist; the flag that sets `optimize_model_rate_joint` is `-jointopt` (`utils/tools.cpp:3381`).
- `ARCHITECTURE.md` section 17 says the fast functional checks have not been run on the project
  machine, and `AA_MODEL_INFERENCE.md` section 16 says the build is not working; both are stale
  since the previous entry.
- `ARCHITECTURE.md` section 17 says googletest targets exist in every build that integrates CMAPLE
  "(the default)"; CMAPLE is not an option on Windows (`CMakeLists.txt:255-258`).
- `FILE_INDEX.md` line 50 places `MIN_RATE`, `TOL_RATE` and `MAX_RATE` at `modelmarkov.h:29-31`;
  they are at lines 30-32.
- `FILE_INDEX.md`'s "Minimum viable context sets" gives `tools.cpp` parse sites 4998 and 5003
  (they are 5018 and 5023) and ranges in `partitionmodel.cpp` and `modelfactory.cpp` that disagree
  with its own Tier 3 anchors.
- None of the three reference documents records the current-code facts listed under the findings
  above (checkpoint precision, the +I+G checkpoint round trip, tip ambiguity handling, the missing
  `ModelsBlock`, `readStateFreq`'s behaviour, `freqTypeString`, the dead restart routine).

A read-only re-check later in the session, at Peter's request (`sed` and `grep` on the three
documents and on each cited source line), located every item above exactly and found more anchors
of the same kind. In `FILE_INDEX.md`'s "Minimum viable context sets" the corrected anchors are:
`tools.cpp` parse sites 5018 and 5023, usage 5978-6030 (the non-reversible model list is at
5998-5999); `phyloanalysis.cpp` 164-179 and 422-461; `partitionmodel.cpp` 61-169, 295-338 and
733-868 (672-709 is the end of `computeMarginalLhForPartitions`, and 807 falls inside a
commented-out block); `modelfactory.cpp` 209-228, 284-305 and 1253-1377 (1259-1273 is a
commented-out function); `modelmarkov.cpp` 2119-2131 (2093-2117 is a commented-out older version
of the same function); `modelprotein.cpp` 1107-1249. Elsewhere, the citation block given as
165-177 at `ARCHITECTURE.md:653` and `FILE_INDEX.md:319`, and as 162-177 at `ARCHITECTURE.md:719`,
is at 164-179; `reportNexusFile`, given as 422-459 at `ARCHITECTURE.md:664`, ends at 461;
`ModelProtein::init`, given as 1106-1245 at `ARCHITECTURE.md:707` and `FILE_INDEX.md:100`, is at
1107-1249; `PartitionModel::targetFunk`, given as 299-336 at `FILE_INDEX.md:284`, ends at 338.

Applied later in the session at Peter's request, as minimal corrections: the first five items
above, at six locations (`AA_MODEL_INFERENCE.md` sections 4.6 and 16; `ARCHITECTURE.md` section 17,
in two places; `FILE_INDEX.md` line 50 and its "Minimum viable context sets", which now carries
the corrected anchors of the previous paragraph). Not applied: the anchors that paragraph lists as
elsewhere, and the missing current-code facts of the sixth item.

### Failed

A read-only search subagent, launched to sweep for other code paths that change model state,
stalled after 600 seconds without reporting. The same questions were answered by direct searches
of optimizer call sites, external writers of rates and frequencies, and `createModel` callers,
which found the `-jointopt` and +I+G checkpoint paths. One question put to Peter referred to
slices S1 and S2 before the plan had defined them, so he could not answer it as asked.

### Next

1. S0: write the regression driver and record the baseline from the unmodified binary (PLAN.md,
   next step), then complete the rest of S0.
2. Peter: confirm or replace the drafted scientific motivation and success criteria in PLAN.md,
   and supply the generating tool and model for the synthesis and for P-log, for
   `AI_DISCLOSURE.md` (carried over).

## 2026-09-23: design imported, reference documents corrected, goal written

### Done

Imported the design stage's three documents into `docs/agent/design/` as byte-identical copies
(the unified discrepancy synthesis and the two independent plans), with a `README.md` giving
provenance, SHA-256 hashes, and precedence. The two plans' hashes match the fingerprints in the
synthesis's section 14.

Recorded decisions in `docs/agent/DECISIONS.md`: entry 001 (all work on the fork's
`nq-constrained-pi` branch, intended for an eventual upstream pull request, per Peter), and the
synthesis's D01 to D10 as a separate section of mathematical and methodological decisions,
adopted by Peter and explicitly revisable. Added dated notes to the pending candidates that the
design contradicts; none was ratified.

Wrote the goal section of `docs/agent/PLAN.md` and added the design folder and
`AA_MODEL_INFERENCE.md` to its document table. Every other section of the plan is still open.

Added a statement of purpose to `ARCHITECTURE.md`, `FILE_INDEX.md` and `AA_MODEL_INFERENCE.md`:
the three describe current code and prescribe nothing. Removed or neutralized the design advice
that had accumulated in them (in `ARCHITECTURE.md`: the "constrained nQ belongs to the indirect
family" callout, the `getNDimFreq` callout, the "what transfers to the 20-state problem" passage
including the claims that nothing downstream needs changing and that the `FREQ_USER_DEFINED`
skip would be sound, the π-constraint parameter derivation, and the prescriptive rows of the
code-location table; about twenty-five similar lines in `FILE_INDEX.md`). Updated `CLAUDE.md`'s
reading order and build section.

Corrections of fact, each read in source at `63c330d9` this session:

- The MSVC build route in `CLAUDE.md` and `ARCHITECTURE.md` section 17 cannot work; replaced
  with a proposed WSL2 and clang procedure, marked as not yet run.
- `-m NONREV+F{...}` is a nonstationary-root model, not an internally inconsistent one.
- In the `NONREV` branch the LG seed is converted before any `+F{...}` vector is read
  (`modelprotein.cpp:1207-1242`), and the conversion uses whatever `state_freq` holds
  (`modelmarkov.cpp:112-119`).
- Linked optimization computes gradients with `Optimization::derivativeFunk` on the
  `PartitionModel` object; no model class can override that path. The only `derivativeFunk`
  override outside `Optimization` is `PhyloTreeMixlen`.
- The `PartitionModel` constructor calls `adaptStateFrequency` for linked `FREQ_ESTIMATE` and
  `FREQ_EMPIRICAL` models and skips it otherwise (`partitionmodel.cpp:116-165`).
- `reportNexusFile` writes 6 significant digits, labels every model `GTRPMIX`, and writes a
  uniform frequency line for a non-reversible Q.
- `-mset` and `-madd` are split on commas (`utils/tools.cpp:588`), so a `+F{...}` model string
  cannot be passed through them.
- `decomposeRateMatrixNonrev` is virtual (`modelmarkov.h:345`).
- The CLAUDE.md partition smoke test used DNA data (`example/example.phy`), which cannot
  exercise the protein `NONREV` path; replaced with `test_scripts/test_data/turtle_aa`.
- The upstream regression harness runs turtle DNA and protein analyses at a tolerance of one
  log-likelihood unit; none of its commands names `NONREV`, `NQ.*`, `GTR20`, `UNREST`, a
  Lie-Markov model, or `--model-joint`.
- **`NDEBUG` correction.** Earlier entries and documents, including this changelog's 2026-09-15
  second-session entry, state that `ASSERT` vanishes in Release builds. For gcc and clang the root
  `CMakeLists.txt` replaces `CMAKE_CXX_FLAGS_RELEASE` with flags that omit `-DNDEBUG` (lines 400
  and 420), and nothing else in the build defines it, so `ASSERT` appears to stay active in
  Release. Read from the build files; to be confirmed from the first build's configure output.
  `AA_MODEL_INFERENCE.md` section 15 and `ARCHITECTURE.md` were corrected accordingly.

State of the machine, checked this session: WSL2 Debian 12 has cmake 3.25.1, gcc 12.2, make and
git; clang, libomp, Eigen and Boost are not installed; `sudo` needs a password; there is no
IQ-TREE clone in the WSL filesystem; the VM has 12 cores and 6 GB of memory. The Windows checkout
has `core.autocrlf=true`, so its shell scripts have CRLF line endings.

Copied the design stage's check scripts `verify_constrained_nq.py`, `check_builtin_matrices.py`
and `chart_fd_bfgs_compare.py` unmodified into `docs/agent/design/scripts/`, separate from
development code, with their SHA-256 recorded in the design README. The fingerprint of
`chart_fd_bfgs_compare.py` matches the synthesis. Left out `e6_n20_big.py`, which cannot run
because its `e4_nonrev` module was never supplied. None of the scripts was rerun.

Later the same day Peter installed clang, lld, libomp, Eigen and Boost in WSL2 Debian and
**built the unmodified branch successfully** (`[100%] Built target iqtree3`, clang 14.0.6,
Boost 1.74, system zlib 1.2.13), reading the source from this Windows working copy and writing
the build to `~/iqtree3-build`, with no second clone. The procedure is now the verified one in
`CLAUDE.md`. From the configure output: the Release CXX flags contain no `-DNDEBUG`, so `ASSERT`
is active; googletest is fetched by `cmaple` at configure time; `CMAKE_POLICY_VERSION_MINIMUM`
is unused with CMake 3.25. An incremental rebuild after touching `model/modelunrest.cpp` took
15.5 s, a no-op build 2.2 s. Pushes now work through the narrowed permission rule (commit
`1767d1ee`), and the fork's GitHub Actions run upstream's CI workflow on every push to the branch.

Smoke tests of the unmodified binary, run from `~/iqtree3-smoke` with `-T 1 -seed 1`, all exit 0:

| Run | Command (after `iqtree3`) | Log-likelihood | Wall time |
|---|---|---|---|
| LG+G4 | `-s example/aa_example.phy -m LG+G4` | -7301.8799 | 16.1 s |
| NQ.pfam | `-s example/aa_example.phy -m NQ.pfam` | -7577.8547 | 20.2 s |
| UNREST | `-s example/example.phy -m UNREST` | -22669.6869 | 3.5 s |
| NONREV | `-s example/aa_example.phy -m NONREV` | -6999.7817 (413 free parameters: 379 rates plus 34 branches) | 67.4 s |
| joint | `-s test_scripts/test_data/turtle_aa.fasta -p test_scripts/test_data/turtle_aa.nex --model-joint NONREV` | -4974.5432 | 155.6 s |
| re-import | same data, `-m joint.Q.txt -te joint.treefile` | -4974.5436 | 0.1 s |
| +F through joint | same data, `-te joint.treefile --model-joint "NONREV+F{0.08,0.06,0.04,0.05,0.02,0.04,0.07,0.07,0.02,0.05,0.10,0.06,0.02,0.04,0.05,0.07,0.05,0.01,0.03,0.07}"` | -4986.6358 | 48.5 s |

Findings from these runs:

- `joint.best_model.nex` records only `NONREV+FO` per partition, without the estimated rates, so
  it cannot carry a learned matrix. The matrix is printed in the `.iqtree` report under "Full Q
  matrix and state frequencies (can be used as input for IQ-TREE)" at **6 decimal places**;
  `joint.Q.txt` is those 21 lines (20 Q rows and the frequency row). Re-importing it reproduced
  the fit's log-likelihood to 4e-4 with no frequency-mismatch warnings.
- Many rates in both joint fits sit at the optimizer floor (printed as 0.000070 and 0.000014
  after normalization) on this small data set (774 sites, 379 parameters).
- A literal `+F{...}` vector in `--model-joint` reached the linked model: the report's frequency
  row equals the supplied vector exactly, Q was fitted with it as a fixed root distribution
  (the nonstationary-root path), and no "Mean state frequencies" pooling occurred. The report
  shows only the linked model's frequencies, so the per-partition claim rests on the source
  trace. The linked model's printed name drops the `+F{...}` suffix.
- `reportNexusFile` is used only for `--link-exchange-rates`, not for `--model-joint`;
  `ARCHITECTURE.md` and `FILE_INDEX.md` were corrected.

### Failed

The first configure failed with `lld not found on PATH`: the install command proposed in
`CLAUDE.md` omitted `lld`, which the build requires for clang on Linux. Fixed by installing it.
Commands run from Claude Code through `wsl.exe` first lost their `$VAR` expansions and had their
`/mnt/c` paths rewritten by Git Bash; the working pattern (a script file, `MSYS_NO_PATHCONV=1`)
is recorded in `CLAUDE.md`. `/usr/bin/time` and `bc` are not installed, so the upstream
`test_iqtree.sh` harness cannot run here as is.

### Next

1. Write the high-level code plan into `PLAN.md` (a fresh agent, from the prompt prepared this
   session).
2. Supply the generating tool and model for the synthesis and for P-log, for `AI_DISCLOSURE.md`.

## 2026-09-15 (second session): complete description of amino-acid model inference

### Done

Wrote `docs/agent/AA_MODEL_INFERENCE.md`, a 1722-line reference describing, as mathematics with
the implementing code, every method IQ-TREE 3 uses to infer amino-acid substitution matrices,
reversible and non-reversible. It covers the seven estimation procedures (fixed empirical matrix,
`GTR20`, `NONREV`, QMaker/nQMaker joint estimation across partitions, profile mixtures and
GTRpmix, PMSF, MUTSEL), the shared numerical machinery, model selection, parameter counting,
output formats, and a list of defects. Section 16 separates what was read directly from what was
established by delegated reading, and what is a derivation rather than a code fact.

Findings that were not in `ARCHITECTURE.md` or `FILE_INDEX.md` and that bear on this project:

- **Gradients for every amino-acid rate and frequency parameter are numerical one-sided forward
  differences**, step `1e-4 * |x|`, `Optimization::derivativeFunk`, `utils/optimization.cpp:916`.
  No substitution model in the tree overrides `derivativeFunk`. One gradient costs `ndim+1` full
  tree traversals, so 380 per gradient for `NONREV`. The gradient convergence tolerance is also
  1e-4, which is the truncation-error scale of the difference itself, so the stopping test sits
  at the noise floor of the quantity it tests. This must be stated in any methods section
  reporting a `NONREV` or `GTR20` fit, and it is the strongest argument for keeping any new
  parameterisation low-dimensional.
- **Box constraints are enforced by clamping the trial point inside the line search, not by
  projecting the direction.** The Armijo test then compares the objective at the clamped point
  against a predicted decrease computed from the unclamped direction. A parameter on a bound with
  an outward direction cannot move, which is why `restartParameters` exists. `ModelMarkov` sets
  `bound_check = false` everywhere, so `GTR20` and `NONREV` never restart; the GTRpmix path
  (`ModelMixture::setBounds`, protein branch) sets it true and therefore does restart randomly.
- **L-BFGS-B is vendored and present but commented out** of `ModelMarkov::optimizeParameters`
  since 2019-09-05 over NaN issues. A genuinely bound-constrained optimiser exists in the tree and
  is not used for matrices.
- **No standard errors are available.** `dfpmin` frees its approximate inverse Hessian on every
  return path with no write-back, so the asymptotic covariance of the rate parameters is
  discarded. Copying `hessin` out before `FREEALL` is the hook if that is ever wanted.
- **`total_num_subst` is divided out at exponentiation time for reversible models and not for
  non-reversible ones.** With the default value of 1.0 the two agree exactly; the asymmetry would
  only bite for a non-reversible mixture component, which is currently unreachable.
- **Root position is a genuinely estimated quantity for non-reversible models**, by local search
  over branches within `root_move_dist` (default 2) during NNI search, plus `--root-test` and
  `--rootstrap`.
- **The nQMaker objective is an unweighted sum of per-partition log-likelihoods**, summed
  serially in a fixed order after the parallel loop specifically so the finite-difference gradient
  is bit-reproducible independent of thread count. The shared matrix is counted once in the
  degrees of freedom, not once per partition.
- **`--init-model DIVMAT` is disabled** behind `ASSERT(0 && "init_by_div_mat not working")`,
  which vanishes under `NDEBUG`, so a release build silently runs a partial version that passes
  the raw row-normalised divergence matrix instead of its logarithm and applies it only to the
  representative partition. Do not use it.
- **ModelFinder is not exhaustive by default.** `filterRates` is on (`ratehet_set = AUTO`) with
  `--score-diff 10.0`, so rate-heterogeneity options are chosen from the first matrix and then
  applied to the other 55 matrix-and-frequency combinations. The default protein candidate set is
  28 x 2 x 22 = 1232 models, and BIC uses the number of sites, not patterns, as sample size.
- **`mutsel_rust/` has now been read.** It is 9 files and 10,839 lines, not the 13 files and ~13k
  lines stated in `ARCHITECTURE.md`. The model is Halpern-Bruno mutation-selection with a shared
  reversible mutation matrix and one free log-frequency vector per site; fitnesses are derived as
  `log pi^s - log pi^mut`, which makes `pi^s` stationary for `Q^s` exactly. It is **reversible**,
  it is fitted by **AdamW with automatic differentiation** rather than BFGS, and it is a **MAP
  estimate** under three quadratic log-space priors, not an ML estimate. It is gated off by
  default because `USE_MUTSEL` is tested but never declared with `option()`.

### Corrections owed to existing documents, not yet applied

Left for Peter to approve rather than edited in place, since the task was to write a new
document:

- `docs/agent/FILE_INDEX.md:202` names a `--matrix-exp` flag that does not exist anywhere in the
  source. The real flags are `--eigenlib`, `--eigen`, `--scaling-squaring`, `--lie-markov`
  (`utils/tools.cpp:5093-5108`), and the default is `MET_EIGEN3LIB_DECOMPOSITION`.
- `docs/agent/ARCHITECTURE.md:36-41` describes `mutsel_rust/` as unread and as 13 files / ~13k
  lines. Both statements are now superseded; see Section 11 of the new document.

### Failed

Nothing was attempted that failed. The build was not exercised, so no claim in the new document
rests on running IQ-TREE.

### Next

Two things follow naturally. First, apply the two documentation corrections above. Second, the
new document's Section 7 states the parameter arithmetic for a pi-constrained nQ (360 free
parameters, nested in 379) and Section 4.5 states why a low-dimensional parameterisation matters
so much given numerical gradients; both belong in `docs/agent/PLAN.md` once the goal is settled.

## 2026-09-15: codebase mapping, documentation scaffold, fork setup

### Done

Mapped the parts of IQ-TREE relevant to non-reversible amino-acid model inference and produced
two reference documents for agents, `docs/agent/ARCHITECTURE.md` (how the system is built, its
conventions, and where this project's seams are) and `docs/agent/FILE_INDEX.md` (a tiered index
of which file to open and when, including an explicit do-not-read list).

Substantive findings, all verified against source at baseline `63c330d9`:

- Non-reversible inference has no module of its own. It is 98 `is_reversible` branch sites
  across 25 live files. `NONREV` is a name string handled inside `ModelProtein::init`
  (`modelprotein.cpp:1207`), not a class. `UNREST` by contrast is a class.
- π is re-solved from Q on every likelihood evaluation, at `modelmarkov.cpp:1267`, by
  `computeStateFreqFromQMatrix` (`modelmarkov.cpp:2119`), an Eigen QR solve.
- `ModelLieMarkov::setBasis` (`modelliemarkov.cpp:1087`) **already implements this project's
  feature for 4-state DNA**: given a target π it shifts the basis matrices so every Q in the
  span has that π as its stationary distribution, reduces the free parameter count by the
  frequency degrees of freedom, and warns when a requested π is unreachable. The 20-state
  generalization is the project. What transfers and what does not is written up in
  `ARCHITECTURE.md` section 9.
- Parameter arithmetic: 380 off-diagonal entries, 379 free under current `NONREV` after removing
  global scale, 360 under a π constraint, because πᵀQ = 0 has rank n-1. The models are nested,
  so a likelihood-ratio test between them is available.
- `-m NONREV+F{...}` is a live but unsound path: it pins `state_freq` and skips the stationarity
  solve while leaving Q unconstrained, so the root distribution is not stationary for the fitted
  Q. This is an argument for a new model name rather than overloading `NONREV`.
- Dead code that misleads searches: `model/modelnonrev.{cpp,h}` are 0 bytes and unbuilt;
  `model/modelgtr.cpp` is unbuilt and includes a header that does not exist; `cmaple/` is a
  separate engine with its own duplicate `NONREV` and `GTR20` handling.

Set up the documentation structure described in `docs/agent/PLAN.md`: `PLAN.md` (living,
overwritten), `ARCHITECTURE.md` and `FILE_INDEX.md` (corrected in place), `DECISIONS.md` (append
only), and this changelog. Added `CLAUDE.md` at the repository root and `model/CLAUDE.md`, which
Claude Code loads automatically.

Fixed the git topology. `origin` now points at the fork `petergoodman/iqtree3` and `upstream` at
`iqtree/iqtree3` with its push URL disabled. Found and removed a real hazard: `branch.master.remote`
still pointed at `upstream` after the remote rename, so a bare `git push` on `master` would have
targeted the official repository. Synced `master` from `8977d31a` to `63c330d9` (tag `v3.1.4`,
64 commits) as a fast-forward, created and pushed branch `nq-constrained-pi`, and merged the
upstream changes into it without conflict.

Re-anchored every line reference in both reference documents against the new baseline.

### Failed

The Release build did not complete. Two causes, both now understood:

1. PowerShell splits the unquoted argument `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` into
   `-DCMAKE_POLICY_VERSION_MINIMUM=3` and `.5`. Confirmed with
   `cmake -E echo -DCMAKE_POLICY_VERSION_MINIMUM=3.5`. CMake then rejects the value for every
   vendored subproject that needs the policy override, which presents as a broken CMake install
   rather than a quoting bug. The fix is to quote the whole token.
2. The bad value was cached in `build/CMakeCache.txt` as
   `CMAKE_POLICY_VERSION_MINIMUM:UNINITIALIZED=3`, so correcting the quoting alone does not
   help. The build directory has to be deleted first, and it was also stale against a
   64-commit jump that changed three `CMakeLists.txt` files.

The corrected sequence is recorded in `CLAUDE.md` and in `ARCHITECTURE.md` section 17 but has
**not yet been run to completion**. Whether upstream `v3.1.4` compiles on this machine is
therefore still unknown.

### Next

1. Run the corrected configure and build. Until that succeeds, nothing else is safe to build on.
2. Read `mutsel_rust/` (13 files, roughly 13k lines of Rust) and `utils/mutsel_wrapper.{cpp,h}`,
   which arrived upstream in v3.1.4 and have not been reviewed. Mutation-selection models bear
   directly on the relationship between a mutational process and stationary amino-acid
   frequencies, so this may overlap with, conflict with, or supply machinery for this project.
   Update the Tier 6 entry in `FILE_INDEX.md` with the verdict.
3. Write `docs/agent/PLAN.md`, which is still a placeholder. In particular settle where the
   target π comes from (literal vector, empirical composition, or named profile), since that
   choice changes the degrees of freedom charged in `getNDimFreq()`.
4. Ratify or reject the four candidate decisions recorded in `docs/agent/DECISIONS.md` and
   convert them into numbered entries.
5. Replace the deny rules in `.claude/settings.json`, which still reference a `data/raw/`
   directory that does not exist in this repository and do not protect the vendored
   subdirectories.
