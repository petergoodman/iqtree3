# SIM-2's seeded GTR20 arm: the approved plan (archive, frozen 2026-11-24)

> **Not current.** The plan Peter approved on 2026-11-24, kept unedited below, with the independent
> read-only review of its first draft appended. The live specification, runbook and pass lines are
> `docs/agent/CODE_PLAN.md` section 3.8 ("SIM-2's seeded GTR20 arm") and section 4 (SIM-2, item 3d),
> the rule is `docs/agent/DECISIONS.md` decision 046 (pass lines in 021), and the risk is
> `docs/agent/PLAN.md` risk 30; those govern. Corrections made while transcribing are listed at the end.

---

# Plan: SIM-2 round 2 records, and the seeded GTR20 arm (decision 046)

## Context

SIM-2's five shards are finished: all 1,640 phases done, none failed, every check row passed, about
19,600 charged core-hours (cores × wall; Puma's accounting not read).

Reading round 2 found that cell D's reversible arms are invalid as run:
- IQ-TREE's GTR20 keeps its last exchangeability, Tyr–Val (`rates[189]`), fixed at the seed's value after
  `rescaleRates` scales the seed's largest entry to 10. The default seed is POISSON (`Params::gtr20_model`),
  so Tyr–Val is fixed at 10.
- `ModelMarkov::setBounds` holds every free rate in [1e-4, 100] (`MIN_RATE`, `MAX_RATE`).
- So no exchangeability can exceed 10 × Tyr–Val. R0's truth needs 498, LG 43, WAG 25, Q.bird 225.

Evidence (kept scripts and outputs in `~/iqtree3-runs/sim2/20261124-round2-read/`):
- In all 80 R0 fits, Tyr–Val is exactly 10 and exactly one rate sits at 100.
- Fitted over true: Tyr–Val 44 to 47; the other 189 rates mostly 0.5 to 1.4.
- Training lnL is 339 below the truth's, which lies in GTR20+F{π*}'s family.
- The relaxed refits are the same, so this is not early stopping.

By reading, without a run, NONREV has the same mechanism: `rates[379]` is fixed from the LG seed, which
allows at most about 427×.

**Peter's decision (2026-11-24):** refit with a seed that removes the limit. Three requirements come with it:
- Record the original failure, because the history matters.
- Make the new Puma work separate from the five shards, identifiable and runnable on its own.
- Give it a distinct code history.

An independent read-only review of the first draft confirmed the `--init-model` mechanics and found five
problems:
- shard assignment;
- the optimizer's scale dependence;
- the core count;
- the floor;
- the wording "rerun".

This plan includes all five.

- **Part A**, done in this session after approval: the records update. It writes Part B into the records
  in full.
- **Part B**, implemented in a new session started in plan mode (project rule: one session per unit, the
  plan recorded first).

## Part A. Records update (this session, one commit, not pushed)

1. **`CHANGELOG.md`**, a new top entry, about 10 lines: "2026-11-24 (twelfth session, continued): round 2
   read; GTR20's range limit; the seeded arm decided (046)".
   - Done: round 2 read (tarball SHA-256 `c10cd93b`); 1,640 phases, none failed; round 2 took 22.93 h at 48
     cores; about 19,600 charged.
   - Failed: round 1's GTR20 arms, kept, with the evidence above.
   - Decided: 046.
   - Next: the Part B session, then `analyze`.
2. **`docs/agent/PLAN.md`:**
   - the status header and the slice table's SIM-2 status;
   - risk 26: the measured 19,600 plus the new arm's projection;
   - **new risk 30**, the range limit, covering:
     - SIM-2's D;
     - D05's GTR20 start and incumbent and 025's carry (S3a's nesting checks hold, but their incumbent was
       probably bound-limited);
     - S7's reversible arm;
     - NONREV by reading.

     Resolved by 046 for SIM-2 and later GTR20 fits, and by question 19 upstream.
   - **Question 19, proposed**, neutral: GTR20, and NONREV by reading, fix one rate at the rescaled seed
     value while bounding the rest to [1e-4, 100], so the attainable range depends on the seed. Would the
     maintainers prefer bounds relative to the fixed rate, another default seed, or a warning at the bound?
   - Current state's SIM-2 bullet replaced.
   - The next step rewritten, with the fixed prompt for Part B's session.
3. **`docs/agent/DECISIONS.md`:**
   - **New 046** (at most 12 lines), "SIM-2 adds a seeded GTR20 arm; every later GTR20 fit uses its seed
     rule". The rule is B1 to B3. Why: the evidence. Not chosen: B2's alternatives. Since 2026-11-24 (Peter).
   - **021**: B6's pass lines as new rows, set by Peter at this approval.
   - **039**: its GTR20 clause gains "seeded per 046".
   - **D05 and 025**: their GTR20 fit is "seeded per 046".
4. **`docs/agent/CODE_PLAN.md`:**
   - Section 3.8:
     - a new block, "SIM-2's seeded GTR20 arm (decision 046)", holding B1 to B8 and Peter's runbook;
     - the Methods GTR20 bullets point to it;
     - the cell F note: its "gtr20" start was a bound-limited fit, and F stands because its three starts'
       B endpoints coincide;
     - the known limits gain risk 30.
   - Section 4, SIM-2: a round 2 read line; new item **3d, the seeded GTR20 arm**, with its task list and
     known shortfalls; item 4 (`analyze`) reads both run directories.
   - Section 5: anchors for `ModelProtein::init` (GTR20 and NONREV branches), `rescaleRates`, `setBounds`,
     `isUnstableParameters`, and the `--init-exchange` and `--init-model` parsing.
5. **`docs/agent/UPSTREAM_FINDINGS.md`**, a new section "Q19 and Risk 30":
   - the source reading, with anchors;
   - the review's trace that `--init-model` seeds every partition's GTR20 under `--model-joint`
     (`partitionmodel.cpp`, a per-partition `ModelFactory`);
   - round 1's numbers, the built-in ranges table and the floor counts;
   - the scripts' paths;
   - the reproduction on the unmodified binary, planned in B5.
6. **`docs/agent/AA_MODEL_INFERENCE.md`**, corrected in place: the "Pinning" paragraph ("harmless … only
   chooses a point on the scale ray, not the fitted model"), "Feasible set" and "Normalisation". Under the
   box, a fixed value p limits every free rate to [1e-4/p, 100/p] × the fixed rate, for GTR20 and NONREV.
   The commit message names the fix.
7. **`docs/agent/archive/SIM2_gtr20_seeded_plan_2026-11-24.md`**: this approved plan, verbatim, with the
   review's findings appended.
8. **Check:** a fresh read-only subagent compares the records with this plan and lists gaps and
   contradictions; I fix them and it re-checks. Budgets are checked: CHANGELOG 10 lines, decisions 12.
9. **Commit**, records only, on `nq-constrained-pi`, with the attribution line. Not pushed.

## Part B. The seeded GTR20 arm (next session: SIM-2 commit 3d)

### B1. What runs, and how it is named
- **A new arm, "GTR20 seeded", not a rerun.** 3.8's rule against rerunning a failed fit with a new seed
  still holds.
- **Fits:**
  - cell D's 160 GTR20 fits: GTR20+F{π*} and GTR20+F, on R0 and NQ, 100 genes, replicates 1 to 10,
    unmasked and s = 0.5, 1, 2;
  - Dr's 16 relaxed refits.
- **Unchanged from round 1:** the same datasets, held-out replicates, data seeds, binaries and settings
  (`-q`, unrooted `-te`, `-blfix`, `-keep-ident`, the 20-round cap; relaxed `-me 1e-6`, 300 rounds,
  2 threads).
- **The only change** is the exchangeability start.
- **Declared now:**
  - "GTR20 seeded" is the primary reversible arm of D's comparison.
  - Round 1's "GTR20 default seed" arm stays in round 1's run directory and is reported beside it as failed
    by risk 30.
  - A seeded fit with a rate at the upper bound is flagged and counted, never rerun.
- **Not run again:**
  - NQC fits;
  - cell F, because its three starts' B endpoints coincided (lnL within 0.003, Ĉ to 5 digits).

### B2. Seed and window
- **Mechanism:** `--init-model seed.paml`, an existing option. For GTR20, `ModelProtein::init` calls
  `readParameters`, then `rescaleRates`. `+F{π*}`, or counting and pooling, then replaces the file's
  frequencies; `readParameters` leaves the frequency type alone (review trace).
- **The seed:** PAML lower triangle, every exchangeability written as 10 and Tyr–Val as p, then π*.
  - The largest entry is already 10, so `rescaleRates` scales by exactly 1 and Tyr–Val is fixed at exactly
    p.
  - The free rates start at 10, as POISSON does today.
  - Window: [1e-4/p, 100/p] × Tyr–Val.
- **p = 0.02** (window [0.005, 5000]), confirmed or replaced by B5's pre-set rule.
  - **Headroom:** the largest unconstrained rate in round 1 reached 343 to 622 × the true Tyr–Val, highest
    at s = 2, and R0's truth needs 498, so 5,000 leaves at least 8×.
  - **Floor:** moving the lower limit from round 1's effective 0.00045 to 0.005 adds at most 0.2 floor-held
    rates per fit (measured on round 1's 160 fits: 1.4 to 6.2 per fit already near 0, depending on
    masking).
- **Caveat:** IQ-TREE's BFGS is not scale-free (`dfpmin` and `lnsrch` test against max(|x|, 1)). With
  p = 0.02 most free rates sit below 1, so rounds, convergence and time can differ from round 1. B5 measures
  this, and Dr's relaxed refits check convergence.
- **Not chosen:**
  - `--init-exchange LG`: at most 427×, below 498. `--init-exchange` takes a model name (built-in or
    `-mdef`), so it would also need an extra NEXUS file.
  - Seeding at the truth: an unfair start.
  - Editing IQ-TREE's bounds: no change to existing inference.
  - p = 0.1: 1,000× leaves under 2× headroom at s = 2.

### B3. Code (one commit, "SIM-2 commit 3d", with tests)
- **New driver `sim/sim2_gtr20.py`:** its own `SCRIPT`, runs root (`~/iqtree3-runs/sim2-gtr20`), design and
  assignment. This is the distinct history Peter asked for.
  - `design(smoke)` takes `sim2.design(smoke)`'s D and Dr GTR20 tasks unchanged except:
    - ids get a `-seeded` suffix;
    - a field `gtr20_start = {"p": P}` (not "seed", which names the random seeds);
    - a wall cap of 40 h (round 1's maximum was 28.4 h against 30; a resource limit, not the 20-round
      setting);
    - estimates from round 1's medians (gtr20pi 13.9 h, gtr20f 16.0 h, relaxed 13.7 and 15.1 h) × B5's
      measured seeded-to-default ratio.
  - It adds `study`, `p` and the seed text's SHA-256 to `design.json`.
  - **Three shards, one wave each:**
    - shard 1, R0's 80 D fits, about 84 cores;
    - shard 2, NQ's 80 D fits, about 84 cores;
    - shard 3, the 16 Dr fits, about 34 cores.
  - It refuses frozen or gate binaries whose SHA-256 differs from round 1's (`9a72950b…`, `25d78d6e…`), so
    the seed is the only change.
  - It uses `run`, `phase`, `design` and `collect` subcommands. `collect` names the tarball
    `collect-gtr20-<stamp>.tar.gz`.
  - Its `design` checks the cumulative budget: rounds 1 and 2's 19,614 plus this projection, against 40,000.
- **`sim/sim2.py`, minimal, behaviour-preserving changes:**
  - `cmd_run` and `cmd_phase` take the design function and script path as parameters, defaulting to the
    main study, so the `phase` subprocess calls the right driver. Without them, a phase raises before its
    `try` and is counted as failed.
  - `do_phase`:
    - passes `gtr20_start` through `fit_into`'s existing `init_model` on its own code path, never through
      `start_matrix` or `init_check`;
    - after a seeded fit, adds rows from the checkpoint (`gtr20_state`):
      - the fixed rate's 10-digit text equals p's, asserted;
      - no free rate above 99 (IQ-TREE's own `0.99 × MAX_RATE` in `isUnstableParameters`), asserted, with a
        hit flagged;
      - the floor count (≤ 1.0001e-4), reported;
      - the number of "Reading model parameters from file seed.paml" lines against the number of genes,
        reported, not asserted (Peter's printed-output rule).
  - `COLLECT` gains `seed.paml`.
- **`sim/pipeline.py`:** `gtr20_seed(p, target_text)`, built with `paml_rows`, and `rate_bounds(state)`.
- **`sim/slurm/sim2_gtr20.sbatch`:** a copy of `sim2.sbatch` with job name `sim2-gtr20`, arguments
  `SHARD RUNDIR FROZEN GATE CLONE`, running `sim2_gtr20.py run`.
- **Tests (`tests/test_sim2.py` and a new `tests/test_sim2_gtr20.py`, no binary):**
  - the main design's SHA-256 equals round 1's `design.json` (`f1556f0d…`), which guards every `sim2.py`
    change;
  - the design: 176 tasks (160 D, 16 Dr), GTR20 only, the same datasets as round 1's tasks, unique `-seeded`
    ids, the shards as above, and every fit starting in the first wave under `project()`;
  - `gtr20_seed`: Tyr–Val is the last PAML entry, equal to p, and the largest entry is 10;
  - `rate_bounds` on a constructed state;
  - `init_model` set only for seeded fits;
  - the binary check refuses a wrong hash.
- **Records in the same commit:**
  - `CLAUDE.md`'s Puma SIM-2 commands gain the new driver and sbatch;
  - CODE_PLAN 3.8's runbook block is finalized.

### B4. Smoke run (laptop, clean committed tree)
- `sim2_gtr20.py run --smoke --skip-checks …` at 3 genes and 2 rounds, from the main smoke design's GTR20
  tasks, end to end: fit, rows, metrics and collect.
- It shows:
  - the fixed rate is 0.02;
  - the seed is read once per gene;
  - gtr20pi's frequencies equal π*;
  - gtr20f's frequencies equal round 1's smoke run (`~/iqtree3-runs/sim2/20261120T114049Z-smoke`) for the
    same data.

### B5. Demonstration and calibration (laptop; new kept script `explore/gtr20_window.py`)
- **Exploratory, with provenance.** About 3 to 6 h under the keep-awake helper, an estimate not yet
  measured.
- **Data:** R0, 30 genes, replicate 1. The committed cross-machine fixture shows the laptop simulates it
  byte-identically to Puma.
- **Four GTR20+F{π*} fits, 1 thread each:**
  - default seed, frozen unmodified binary;
  - default seed, gate binary;
  - seeded at p = 0.02, gate binary;
  - seeded at p = 0.05, gate binary.
- **Recorded for each:**
  - training lnL against the truth's (the oracle's `wiring_objective` at R0);
  - the fixed rate, the largest rate and headroom (100 / largest), and the floor count;
  - rounds, ending rule and wall time.
- **Rule for p, set now:** p = 0.02, unless the p = 0.05 fit's training lnL is higher by more than 1 and its
  headroom is at least 4. Then 0.05.
- The seeded-to-default wall ratio scales B3's estimates.

### B6. Pass lines (021 rows, set by Peter at this approval; a failure stops the work for Peter)
- **B5:**
  - the frozen and gate default-seed fits are identical (checkpoint rates text, final lnL);
  - each seeded fit: the fixed rate's text equals p; no rate above 99; training lnL ≥ the truth's.
  - The default fits' limit is reported. If it does not bind at 30 genes, that is reported, and the
    100-gene evidence stands.
- **Every Puma seeded fit:**
  - the fixed rate equals p;
  - no rate above 99. A hit is flagged and counted, never rerun.
- **`analyze`:**
  - Data identity: simulated datasets' gene SHA-256 equal round 1's; for masked datasets, `data.json`'s
    `files` hashes and `curve_sha256` equal round 1's.
  - gtr20f frequencies equal round 1's at the checkpoint's 10 digits.
  - gtr20pi frequencies equal π*.
  - R0 unmasked gtr20pi seeded fits' training lnL ≥ the truth's. Reported; a shortfall is read against
    the same data's Dr relaxed refit.

### B7. Puma (Peter's runbook in CODE_PLAN 3.8)
1. **Push** of 3d and B4/B5's record, with Peter's approval.
2. **Login node:** in the clone, `git pull`, then check out that commit. The tree must be clean.
3. **Test job:** 40 min at 4 cores, from a scratch directory: `sim2_gtr20.sbatch 3 <scratch>/test …`.
4. **From a new run directory** (for example `…/sim2-gtr20-run`), outside the clone, for K = 1, 2, 3:
   `sbatch --account=<account> --partition=standard --cpus-per-task=<84|84|34> <clone>/test_scripts/fixedeq/sim/slurm/sim2_gtr20.sbatch K <run dir> <frozen> <gate> <clone>`
5. **After the jobs end:** `sim2_gtr20.py collect --out <run dir>` on a compute node; WinSCP; copy into
   `~/iqtree3-runs/sim2-gtr20/`; the agent reads it.
6. **Expected:** one wave of about 25 to 40 h and roughly 5,000 to 8,000 charged core-hours, giving a total
   of about 25,000 to 28,000 of 40,000. These are rough; the design step projects them from B5's ratio
   before submission.

### B8. `analyze` (commit 4, written while Puma runs)
- It takes both run directories.
- Cell D's reversible comparison uses the seeded arm, with the default-seed arm reported beside it as
  failed (risk 30).
- Rate metrics are given with and without floor-held pairs, for both arms.
- It runs B6's `analyze` checks.

## Verification
- **Part A:** the fresh subagent's records-against-plan check passes after fixes; the line budgets hold.
- **Part B:**
  - pytest passes in `iqtree3-fixedeq`, including the main design's hash;
  - the smoke run has all phases done;
  - B6's laptop lines pass;
  - the cumulative projection is under 40,000;
  - every asserted row in Puma's three reports passes;
  - `analyze`'s checks pass.

## Critical files
- **New:** `test_scripts/fixedeq/sim/sim2_gtr20.py`, `sim/slurm/sim2_gtr20.sbatch`,
  `tests/test_sim2_gtr20.py`, `explore/gtr20_window.py`.
- **Changed:** `sim/sim2.py` (parameters, `do_phase` rows, `COLLECT`), `sim/pipeline.py`,
  `tests/test_sim2.py`.
- **Records:** `CHANGELOG.md`, `CLAUDE.md`, `docs/agent/{PLAN,CODE_PLAN,DECISIONS,UPSTREAM_FINDINGS,AA_MODEL_INFERENCE}.md`,
  `docs/agent/archive/`.
- **Reused:**
  - `sim2.checks`, `sim2.startable` and `sim2.project`, `sim2.scan`, `sim2.metrics`;
  - `pipeline.fit_into` (`init_model`), `pipeline.gtr20_state`, `pipeline.paml_rows`,
    `pipeline.matrix_metrics`, `pipeline.wiring_objective`, `pipeline.simulate_into`.
- **Read-only IQ-TREE references:** `model/modelprotein.cpp` (`ModelProtein::init`, `rescaleRates`),
  `model/modelmarkov.{h,cpp}` (`MIN_RATE`, `MAX_RATE`, `setBounds`, `isUnstableParameters`,
  `readParameters`), `utils/tools.cpp` (`--init-exchange`, `--init-model`), `utils/optimization.cpp`
  (`dfpmin`, `lnsrch`).

---

## Appendix A. The independent review of the first draft (2026-11-24)

A fresh read-only reviewer that had not seen the session was given the first draft and the code and
asked for defects. Its findings, as returned (the plan above already answers them):

1. **`assign()` cannot split this study over 2 shards** (`sim2.py:assign`). D and Dr are in `PINNED`, so
   every task goes to the last shard; non-pinned tasks get `free = range(shards-1)`, which is `[0]` with
   2 shards. The study needs its own assignment, and the main study's `assign()` must stay untouched.
2. **p changes how the optimizer behaves, not just the window.** IQ-TREE's BFGS (`utils/optimization.cpp`)
   is not scale-invariant: `dfpmin` tests step size and gradient against `max(|x|,1)`, starts from an
   identity Hessian and has an absolute `stpmax`; `lnsrch` uses the same `max(|xold|,1)` test. At
   p = 0.02 almost every fitted free rate is below 1, so rounds, convergence within 20 rounds and wall
   time can differ, and round 1's durations are not valid estimates. The 30 h cap against a measured
   27 to 28 h maximum leaves 7% margin; `do_phase` counts a capped fit as failed, which blocks its metrics.
3. **The bound affects more than SIM-2.** D05 and 025 use the fitted GTR20+F{π*} as the nesting
   incumbent; 039's "S7 uses this setting" applies the same GTR20 fits to S7. NONREV pins `rates[379]`
   at LG's rescaled value (`modelprotein.cpp`, NONREV branch), a window of about 427× at most. Decision
   046 and risk 30 should cover S7 and D05/025, and the `AA_MODEL_INFERENCE.md` fix should cover NONREV
   (the pinning table and "harmless" text, and "Feasible set" and "Normalisation").
4. **Cores: 192 > 188.** With an even split each shard needs 96 cores, so two fits per shard start late.
5. **p = 0.02 rests on weak support.** Round 1's maximum of 782 includes the rate held at 100, so it is
   censored. Tyr–Val is among R0's smallest exchangeabilities, so any error in it shifts every ratio. The
   floor goes up, from about 0.00045 × the true Tyr–Val to 0.005. The NQ fits were not looked at.
6. **`--init-model` mechanics check out (evidence, not a defect).** Each partition builds its own
   `ModelFactory` (`partitionmodel.cpp:74`), and `ModelProtein::init`'s GTR20 branch reads the file for
   every partition; the linked model is partition 1's object; the bounds come from
   `ModelMarkov::setBounds`, `bound_check` false; `rescaleRates` runs only at init; `readParameters` never
   changes `freq_type`; `+F{…}` frequencies are re-read and `+F` frequencies are counted then pooled.
   NQC's handling of `--init-model` lives only in `ModelNonrevFixedEq::readSeed`. The stated reason for
   rejecting `--init-exchange` was wrong: `findModel` also finds models defined with `-mdef`.
7. **Seed text and checks.** Write the free rates as 10 and Tyr–Val as p, so `rescaleRates` scales by
   exactly 1; compare the checkpoint's 10-digit text, not floats; use IQ-TREE's own threshold,
   `> 0.99 × MAX_RATE`; check "Reading model parameters from file seed.paml" once per gene; add
   `seed.paml` to `COLLECT`, and p and the seed's SHA-256 to `design.json`.
8. **Wiring in `sim2.py`.** `start_matrix`'s else branch writes NQ.bird transposed for any `init` value,
   and `init_check` runs whenever `init` is set: keep the GTR20 seed on its own field and path. A field
   named `seed` collides with the random seeds; the metrics branch of `do_phase` returns no rows.
9. **The `--study` plumbing is fragile.** `cmd_phase` calls `design(a.smoke)` again and `cmd_run` starts
   phases as `Path(__file__) phase` without the study; a missed pass-through raises before the `try`.
10. **Binaries are not pinned.** CR may produce a new gate binary; enforce round 1's hashes.
11. **Checks planned for `analyze`.** "Training lnL ≥ truth" applies only to GTR20+F{π*}; a 20-round fit
    may stop short, so a shortfall is a convergence finding; masked datasets' records have no per-gene
    hashes, so compare `data.json`'s file hashes and `curve_sha256`.
12. **Laptop validation.** "Seeded above the truth" is not guaranteed at 20 rounds; about 3 h is
    optimistic.
13. **Budget.** The ceiling must be checked cumulatively: round 1 about 18,500, round 2 about 1,100, plus
    this study.
14. **Records.** Call the seeded runs a new arm, not reruns (3.8: failures are never rerun with a new
    seed); declare the primary arm in advance; a seeded fit at the upper bound is flagged, not rerun;
    note that cell F's "gtr20" start was a bound-limited fit.

Its recommendations, adopted: a separate driver with its own design, assignment and estimates, `sim2.py`
parameterized and guarded by round 1's design hash; three jobs (R0's D, NQ's D, Dr); a 40 h cap; a
laptop comparison of the default seed and p = 0.02 and 0.05, with the seeded-to-default time ratio
scaling the estimates and a floor count from round 1; rate metrics with and without floor-held pairs;
round 1's binary hashes enforced.

## Appendix B. Measurements and corrections after approval (2026-11-24)

- The floor count in the plan above ("adds at most 0.2 floor-held rates per fit") was misread from an
  unkept inline script. Its kept rerun (`~/iqtree3-runs/sim2/20261124-round2-read/floors.py` and
  `floors.txt`) gives 1.3 to 6.0 rates a fit at IQ-TREE's floor and 1.4 to 6.2 below 0.005, so at most
  0.3 more. The records state 0.3.
- "Largest unconstrained rate 343 to 622" in the plan is stated in the records as medians of 343 to 560
  by masking strength, at most 622.
- A second fresh read-only reviewer checked the records against this plan. Its corrections, each
  measured in `~/iqtree3-runs/sim2/20261124-round2-read/window_checks.py` (`.txt`) or read in source:
  - NONREV's window: "about 427 ×" was GTR20's figure under an LG seed. NONREV converts its rescaled LG
    seed with LG's frequencies (`setReversible(false)`), so its fixed Val→Tyr rate is about 0.0080 and
    its window 0.0125 to 12,500 × that rate, which holds NQ.bird's range (0.33 to 2,840) but not LG's
    lowest rates (0.0053). The records state this.
  - IQ-TREE already prints a boundary warning after the model fit when `isUnstableParameters` holds
    (`IQTree::optimizeModelParameters`), so question 19 no longer offers "a warning" as an option, and
    the seeded fits report the warning's count. Whether round 1's fits printed it was not checked.
  - The longest default-seed GTR20 fit over rounds 1 and 2 took 27.5 h, not 28.4 h (that figure was a
    percentile above the maximum); GTR20+F{π*}'s median over all 80 fits is 13.3 h, not 13.9.
  - All 160 of cell D's default-seed fits, masked included, have Tyr–Val at 10 and exactly one rate above
    99; "44 to 47 times too high" holds for unmasked R0, and over all 80 R0 fits the range is 43.7 to 78.2.
  - Cell F's starts agree on Ĉ within 1e-5, not to 5 digits.
  - The records also state how a bound hit acts (the row fails in the shard's report, the phase stays done,
    the fit is counted), the design's first-wave check (a calibration ratio above about 1.9 stops the work
    for Peter), the smoke run's pass condition (every phase done), and the `analyze` rows' handling.
- B6's "a failure stops the work for Peter" is applied to the Puma rows when the agent reads the tarball,
  not mid-job: a failed fixed-rate or bound row fails in the shard's report, the phase stays done, the jobs
  go on, and the analysis of D waits for Peter. Stopping running jobs would waste their other fits.
