# CR reuse review (archive, frozen 2026-11-21)

> **Not current.** CR's step 1, the read-only reuse review, as Peter accepted it on 2026-11-21, kept
> unedited for its findings and evidence. CR's live task list is `docs/agent/CODE_PLAN.md` section 4
> ("CR") and its open choices are `docs/agent/PLAN.md` risk 29; those govern. Peter's answers are at
> the end. The three reviewers' full reports are in the companion file
> `CR_reuse_review_2026-11-21_appendix.md`.

---

Exploratory, read-only review at `91a0a45d` (`nq-constrained-pi`), CR step 1 (`CODE_PLAN.md` section 4,
CR). Nothing was built or run.

## Answer

**Principle 1 (code reuses ModelMarkov).** NQC already reuses ModelMarkov's whole non-reversible
decomposition (unpacking, normalization, eigensystem), the transition matrices and likelihood kernels,
the unchanged optimizer, the `+F` routes, the checkpoint machinery, parameter counting and
PartitionModel's linked loop. What it copies falls into four groups:

- about 50 lines that could reuse IQ-TREE now with outputs bitwise unchanged (section 1.2);
- about 75 lines that are copies only because decisions 009 and 015 keep `utils/optimization.*`
  and `model/modelmarkov.*` unedited (section 1.3);
- one copy whose reuse changes the seed's last bits (section 1.4);
- the rest is new mathematics the constraints require (section 1.5).

**Principle 2 (usage matches NONREV).** A nQMaker user changes one token: `--model-joint NONREV+FO`
becomes `--model-joint NQC+F{π*}`. Of 70 command forms, 36 are the same, 10 differ, 18 are refused and
6 are new. The largest departures that the fixed π* does not require are the mandatory end-of-run
export (with the three refusals it causes) and the refusal of unlinked partitions.

**The maintainer's `-S` folder question.** By reading, a folder reaches NQC exactly as it reaches
NONREV. No NQC run has ever used a folder, or the published joint step's form (a `.best_model.nex`
whose charsets each name a file). The CR plan adds those runs, as CR step 1 requires.

## Method and labels

Three fresh reviewers worked in parallel, read-only, each given the code, IQ-TREE's base classes and
the decisions, and told not to read the CR brief, risk 29, maintainer question 18 or the CR block:

1. the model class, the hooks header and each shared-file edit;
2. the pure chart;
3. the user-facing commands, including IQ-TREE's tutorial (read from the iqtree.github.io mirror,
   since www.iqtree.org's certificate had expired).

Every claim below carries a label:

- **V**: checked against the source in this session by the main session;
- **R**: the reviewer's reading, not rechecked;
- **D**: derived from the source, not executed.

Every "bitwise" claim is by reading. The CR plan's equivalence runs against `25d78d6e` are what would
confirm one. Line counts are the reviewers' estimates.

File abbreviations: C `model/modelnonrevfixedeq.cpp`, CH `model/fixedeqchart.cpp`, MM
`model/modelmarkov.cpp`, MP `model/modelprotein.cpp`, PM `model/partitionmodel.cpp`, OPT
`utils/optimization.cpp`, TL `utils/tools.cpp`, SA `alignment/superalignment.cpp`.

## 1. Code reuse

### 1.1 What NQC reuses now

- `decomposeRateMatrixNonrev` is called unchanged: unpacking `rates[]`, the diagonal, the rescale to
  unit mean rate, the eigensystem (C:510-537 calls MM:1246-1290; V).
- Transition matrices, the scaling-and-squaring fallback and the likelihood kernels (R).
- ModelProtein's ambiguity-code tip likelihoods and its model readers (R).
- `minimizeMultiDimen`, `dfpmin`, `lnsrch` and `fixBound`, unchanged (R).
- The `+F{...}` and `+F<name>` routes in createModel and ModelFactory (R).
- The checkpoint machinery, `getNDim` and `getNDimFreq` (R).
- PartitionModel's linked loop and its fixed-order sum (V, PM:299-338).

### 1.2 Reuse possible now with outputs bitwise unchanged (tier A candidates)

| # | NQC site | IQ-TREE code to use instead | Lines | Evidence |
|---|---|---|---|---|
| 1 | `readSeed`'s hand-copied rooting block for a non-reversible seed (C:268-276) | an unconditional `setReversible(false)` (MM:75-152) | about 6 | V. The recorded reason not to call it, that it zeroes `rates[]` (`CODE_PLAN.md` 2.2 item 4), does not hold. `rates[]` has already been copied into Q0 (C:253-260), and the constructor's decomposition rewrites it from the coordinates (C:152-156, 530-533) before anything reads it. The call does the same `convertToRooted` with the same `VB_MED` message; `num_params` is overwritten at C:152; `ignore_state_freq` affects only reversible paths (MM:1503, `eigendecomposition.cpp:488`) |
| 2 | The linked-partition list built in `PartitionModel::optimizeLinkedModel` (15 added shared-file lines), which `reportRun` builds again (C:682-689) | NQC builds it in `startLinkedFit` from the PartitionModel and keeps it for `endLinkedFit` | about 12 shared-file lines | V for both loops (same name filter, same `part_order`); the refactor is not designed. It changes the hooks' signatures, so it goes with the hooks choice |
| 3 | Chart-internal repeats: K = q/(−q_ii) three times (CH:256-260, 319, 343), the exit-rate check twice (251-255, 313-317); `freeDestinations` (CH:142-148) has no caller | one private helper that keeps the division; delete `freeDestinations` or call it from `jumpMatrix` and `inverse` | about 19 | V |
| 4 | `isName` (C:106-110) | `iEquals` (TL:210-219) | 3 | V; equal for ASCII names |
| 5 | The constructor's `freq_type = FREQ_USER_DEFINED` and `decomposeRateMatrix()` (C:151, 156) | `ModelMarkov::init(FREQ_USER_DEFINED)` (MM:354-359) | 1 | V. For protein user-defined frequencies `init_state_freq` only sums them (MM:334-341); `writeInfo` prints nothing for a 20-state non-reversible model (MM:361-380) |
| 6 | `FIRST_ORDER_GTOL` at the class's call (C:787) | `TOL_RATE` (`modelmarkov.h:31`), also 1e-4 | 0 | V; the chart keeps its constant for the standalone tests |
| 7 | `number17` (C:37-42) duplicates the chart's `str` (CH:110-116) inside the fork | one function, with the classic locale | about 5 | V for the duplicate; same bytes unless a non-classic global locale is set (R: no `setlocale` found) |
| 8 | `STATES` (C:33) | `symbols_protein` (`alignment.cpp:34`) | 1 | R |
| 9 | The `models` member (`modelnonrevfixedeq.h:104`) duplicates ModelProtein's private `models_block` | pass it to `readSeed` | 2 | R |

Bitwise but not recommended:

- **`setFullRateMatrix` for the packing loop (C:530-533).** NQC overrides that function with a guard,
  so it needs a qualified call, and Eigen's column-major Q needs a 400-entry row-major copy at every
  evaluation.
- **Reusing the likelihood that `optimizeParameters` computes (C:383-387) in place of `rescore`.** It
  saves one evaluation per fit, but removes decision 009's independent re-score.

### 1.3 Copies that exist only because upstream files are off limits

Each reuse here needs an edit to a shared upstream file. Numbers stay bitwise for other models and
for NQC (D) if the edit is made as described.

| # | Copy | NQC lines | Upstream edit that would remove it | Barred by |
|---|---|---|---|---|
| 1 | `optimizeParameters` (C:357-395) repeats MM:1166-1244 statement for statement (V), minus the unused `variables2` and the reversible-only `scaleStateFreq` block, plus a kept `optimizer_score` and `checkFit` | about 30 | ModelMarkov keeps the optimizer's score in a member; NQC calls the base, then `checkFit` | decision 015 (`modelmarkov.*`) |
| 2 | `scaledDerivative` (C:327-355) repeats `Optimization::derivativeFunk` (OPT:916-939) with the step swapped (V); `STEP_ETA` equals `ERROR_X` | about 14 | a step hook in `derivativeFunk` whose default reproduces OPT:928-931 | decision 009; decision 005's "Not chosen" |
| 3 | `solveStationary` (CH:190-204): the same augmented system and Eigen `colPivHouseholderQr` as `computeStateFreqFromQMatrix` (MM:2119-2131; V) | about 10 | move the body into a header both use; `modelmarkov.cpp` keeps its `ASSERT`. A direct call is not possible: the `ASSERT` aborts in Release builds (R), turning constraint 2's rejection into a crash | decisions 015 ("a shared helper split out of the base" not chosen) and 012 |
| 4 | `completeDiagonal` (CH:57-71): the same loop and summation order as MM:1256-1264 (V) | about 8 | a diagonal helper in `modelmarkov.cpp` | decisions 015 and 012 |
| 5 | `firstOrderTest`'s scaling (CH:435) repeats `dfpmin`'s test (OPT:849-857) | about 3 | a helper in `optimization.cpp` | decision 009 |
| 6 | The 17-digit checkpoint writers `numbers17` and `number17` (C:37-49) | about 12 | a defaulted precision argument on `Checkpoint::put`/`putArray` and `convertDoubleToString` (`utils/checkpoint.h`, `utils/tools.*`) | no decision (R) |
| 7 | `readSeed`'s LG and `--init-model` read (C:222-245) repeats `ModelProtein::init`'s NONREV branch, which itself repeats the GTR20 branch (MP:1191-1198, 1213-1239) | about 9 (plus upstream's own duplicate) | a protected ModelProtein helper used by all three, a pure move for GTR20 and NONREV | no decision; `modelprotein.*` is already edited by the fork (R) |
| 8 | `linkedmodelhooks.h` and the two `dynamic_cast`s in PartitionModel | 49 plus 2 | ModelSubst virtuals for the linked derivative, start and end, with defaults that keep other models unchanged (`model/modelsubst.h`) | no decision. This is a third form for risk 29's hooks item, besides moving the interface or a direct cast (R) |

Recommendation:

- Keep decisions 009 and 015 for CR.
- Offer items 1 to 5 to the maintainers as one concrete proposal (an edit of a few lines to
  ModelMarkov and Optimization in exchange for about 65 lines of copied code removed). The copies
  exist because of their earlier advice of caution (decision 015's Why, D06), so the trade is theirs
  to make.
- Treat items 6 to 8 as Peter's choices: each touches an upstream file that no decision protects.

### 1.4 Reuse that would change numbers (tier B)

- **`reversibleSeed` (CH:280-302, about 23 lines).** It is identical to `getQMatrix`'s reversible path
  through `computeRateMatrix` (`eigendecomposition.cpp:478-519`), except that it divides by the mean
  rate where IQ-TREE multiplies by its reciprocal (V). The seed's last bits would change, and with
  them the start, every fit, the `.qmat` and the checkpoint bytes. `readSeed` also has
  `setReversible(false)` compute R_ij π*_j into `rates[]` and then discards it (C:267-270; V).
- **Checkpoint keys through `Checkpoint::getVector`/`putVector`.** The separator changes from "," to
  ", " (R).
- **The `+H` test through `posRateHeterotachy`.** This changes only edge-case strings (R).

### 1.5 Departures the constraints require

- **D01:** `jumpMatrix`, `build`, `inverse`.
- **Constraint 1:** `generatorFromJump` and `transferSeed`. IQ-TREE moves a Q to new frequencies by
  column scaling (`adaptStateFrequency`, MM:914-930), which does not keep π stationary for a
  non-reversible Q (D09 item 6).
- **Constraint 3, D07 and decision 006:**
  - `getStateFrequency` without the base's renormalization;
  - the frequency-setter guards;
  - the immutable `target` copy;
  - `parseTarget` in place of `readStateFreq`.
- **Constraint 2 and decision 003:** `targetFunk`'s 1e30 return, and NQC's own
  `setVariables`/`getVariables` with trial and accepted arrays.
- **Constraint 8 and decision 007:** `setBounds`. The base box `MIN_RATE`/`MAX_RATE` excludes every
  negative coordinate.
- **Constraint 6 and decision 003:** NQC's own 17-digit checkpoint keys. `CKP_PRECISION` is 10.
- **Decision 021:** `diagnostics`.
- **Decisions 002 and 016:** the protected constructor.
- **No existing mechanism:** the createModel branch, the `reportRun` call site and the linked hooks.
  `writeInfo` is silent for a 20-state non-reversible model and runs per fit; `report()` is reached
  only for PoMo (R).

### 1.6 Reasons in the records that do not hold, or are imprecise

1. **`CODE_PLAN.md` 2.2 item 4** gives the zeroing of `rates[]` as the reason not to call
   `setReversible(false)` on a non-reversible seed. The zeroing is harmless at that point (1.2 item
   1; V).
2. **The comment at C:313**, "pi* was checked against it". Raw entries are checked, then divided by a
   sum within 1e-6 of 1, so a normalized entry can sit just below `min_state_freq`. The base guard
   (MM:1103-1108) would then reject every evaluation, so omitting it is right, but the comment's
   reason is imprecise (V).
3. **`fixedeqchart.h:106`**, "That displacement is exactly h0 when |x| >= h0". This contradicts
   `CODE_PLAN.md` 2.1 and the unit test (`test_fixedeqchart.cpp:422-442`), which show the subtraction
   is exact but the displacement can differ from h0 (V).
4. **Decision 006's "random draws".** IQ-TREE's reader draws only for "uniform" or a defined
   distribution name; any other token is an error (V, TL:688-699). Its other reasons hold (R): NaN
   passes `strtod`, sums are renormalized with only a warning, and zeros are accepted.
5. **`fixedeqchart.h:110-111`**, "the gtol NQC's fits receive". This holds only when
   `gradient_epsilon` ≤ 1e-4 (R).

## 2. User-facing usage

### 2.1 Overview

- **Same as NONREV:** `-s`, `-S`, `-p`, `-q`, `-Q`, `--model-joint`, `--link-model`, `-te`, `-n 0`,
  `--init-model` (syntax), `--link-alpha`, `-blfix`, the stopping options, `-T`, checkpoint and
  resume, `--undo`, `-b`/`-B`, `--show-lh`, `-m FILE` re-import, the citation (Peter, 2026-10-26),
  and parameter counting (R; 11 of the 36 "same" forms have no NQC test).

### 2.2 The `-S` folder question

- **Published workflow.** The tutorial and the nQMaker paper give the folder only to the ModelFinder
  steps on ordinary models. The joint step is `-S DIR.best_model.nex -te DIR.treefile --model-joint
  NONREV+FO --init-model LG`, a NEXUS file whose charsets each name a file in the folder. For the
  paper this is V from the deck's note (`nq_frequency_replacement_part2.tex:208`); for the tutorial
  it is R, from the mirror.
- **Code.** `-S DIR` makes `readFromParams` call `readPartitionFiles`, which builds one partition per
  file with `-m`'s model string (V, SA:91-96, 689-722). From there both forms take one path.
  - NQC's refusal block tests only `partition_file`, `model_joint` and `link_model` (V, C:203-206).
  - No NQC code reads the input layout (R).
  - As for NONREV, a folder with no `-m` gives every partition no rate model, and per-gene rate models
    come only through the `.best_model.nex` form (R).
- **Tests.** No NQC test uses a folder, or charsets that name separate files. Every NQC `-S` run is
  `-s <concatenation> -S <charsets with position ranges>` (V for `nqcchecks.py`). The search did not
  cover unscripted runs kept only under `~/iqtree3-*` or on Puma.
- **For the maintainers.** It works by reading, as for NONREV, and has not yet been run. The CR plan
  adds a folder run and a tutorial-form joint run.

### 2.3 Departures the fixed π* requires

- **A target is mandatory.** Every untargeted `+F` form is refused, with a message pointing to NONREV.
- **A new name.** `NONREV+F{...}` keeps its live meaning (constraint 9).
- **Starts are moved into the π* family.** A reversible start is rebuilt at π*; a non-reversible
  start's jump chain is transferred (D05). Same syntax.
- **Strict target parsing** (D07, 006).
- **Linked partitions must hold one target.** NONREV links by name unchecked.
- **An incompatible restart is refused**, where IQ-TREE only warns (constraint 6).
- **Refused combinations:** `-optfromgiven`, `-jointopt` (its bounds would clamp the log-ratio
  coordinates), `-fs`/`-ft`, `+FMIX`, and non-protein data.

### 2.4 Departures the fixed π* does not require

1. **The mandatory end-of-run export.** Every run writes `<prefix>.NQC.qmat` and `<prefix>.NQC.info`,
   prints a summary block, and runs a diagnostic of 722 likelihood evaluations, 6 to 14% of a first
   linked fit (`CODE_PLAN.md` 2.2; decisions 009, 010).
   - It causes the `--runs` and `--root-test` refusals and part of the `-mstore` refusal.
   - The `.qmat` repeats the `.iqtree` full-Q block, which NQC also prints, at 17 digits instead of 6
     (decision 010: a 6-decimal re-import missed by 4e-4).
   - Decision 009's reason, that a return from `dfpmin` certifies nothing, applies to NONREV too
     (maintainer question 3).
2. **Unlinked partitions are refused** (scope). NONREV fits one model per partition.
3. **ModelFinder is refused** (scope, S6). nQMaker's ModelFinder steps run on ordinary models and are
   unaffected.
4. **AliSim from the model string is refused, and the AliSim command printed in `.iqtree` fails.** It
   stops at "needs its target", because the frequency check precedes the AliSim check (V, C:180-201).
   `--alisim -m <prefix>.NQC.qmat` works through the existing reader (R).
5. **`--nqc-domain LO,HI`** is new (decision 007, D03). NONREV's box is fixed at compile time. Model-
   string braces (`NQC{lo,hi}`, refused now) would avoid the `utils/tools.*` edit, at the cost of a
   syntax change (R).
6. **New prints:** the `-v`/`-vv` records and the fall-below-start warning (decisions 009, 021).

### 2.5 Refusals that protect against upstream defects NONREV users also meet

- **`--lie-markov`** (maintainer question 14).
- **`--eigen`** (an `ASSERT`; R).
- **A rate model written inside `--model-joint`** (13).
- **`-mstore`** (9).
- **`--root-test`** (16).
- **`--runs` under `-S`** (11).

### 2.6 Message and help fixes with no loss

All of these change user-visible text, so none is bitwise:

- **Unlinked refusal.** It names `-p, -q or -S`, but `-Q` is also refused, since it sets
  `partition_file` (V, TL:2076-2081).
- **ModelFinder refusal.** "give it with -m" misleads a `-m MFP --model-joint NQC...` user (R).
- **AliSim refusal.** It could name `-m <prefix>.NQC.qmat` (R).
- **Help text.** "<name> from --mdef" understates: built-in frequency vectors also resolve (R, by
  reading).
- **Overfitting warning.** NONREV's warning (MP:1214-1217) has no NQC counterpart; it would matter
  only for single-alignment runs (V). Cosmetic.

## 3. Side findings outside reuse (not acted on)

1. **`ModelMarkov::computeTrans(time, s1, s2, d1, d2)` (MM:590-604)** reads the real eigen arrays
   whatever `is_reversible` says (V).
   - For a non-reversible model those arrays hold the complex system's storage (`model/CLAUDE.md`), so
     its result would be wrong for NONREV and NQC alike.
   - Callers are in `alignmentpairwise.cpp` (R).
   - Whether any run reaches it is unverified. A candidate maintainer question after a check.
2. **`isUnstableParameters`** flags derived rates above 99.
   - With unit mean rate, q_ij ≤ 1/(2π*_i), since a state's out-flux is at most 1/2 (checked), so it
     can fire when some π*_i < 1/198.
   - The check is meaningless in NQC's coordinates; where it prints is R.
3. **NQC's `restoreCheckpoint`** parses with `convert_double_vec` and `convert_double`, which throw a
   string. A malformed key may end in `std::terminate` rather than `outError` (R; no catch found in
   `main/main.cpp`, not traced further).
4. **The MrBayes block** writes "gtr" for NONREV and NQC (R).
5. **Checked and dismissed:** the concern that a still-fixed model would make the export diagnostic's
   probes evaluate a single point. `PartitionModel::targetFunk` unfixes each partition model around
   its own `targetFunk` call (PM:318-321), so every linked partition sees each probe (V).

## 4. What this means for risk 29's open choices

- **Chart layout.** The chart duplicates IQ-TREE in about 45 lines: items 1.2.3, 1.3.3 to 1.3.5 and
  1.4's seed. Keeping it as a second pair, option (a), therefore costs little duplication. Merging
  the files removes no duplicated code.
- **Hooks.** A third form now exists beside "move the interface" and "direct cast": ModelSubst
  virtuals with defaults (1.3 item 8). It removes the header and the casts at the price of a
  `modelsubst.h` edit. Item 1.2.2 changes the hooks' signatures, so it belongs to the same choice.
- **Tier B.** Section 1.4, and the message fixes of 2.6 if they are not allowed into CR.
- **The `-S` folder runs** go into the CR plan's equivalence list (2.2).

## Answered by Peter, 2026-11-21

- **Recorded.** This review and its appendix are frozen in `docs/agent/archive/`, with `PLAN.md`,
  `CODE_PLAN.md` and `CHANGELOG.md` updated.
- **Message and help fixes (section 2.6): yes.** They come into CR as one labelled commit, and the
  refusal comparisons' expected differences are written down before the equivalence run.
- **The copies of section 1.3, items 1 to 5: open.** Peter leans towards no: he would rather keep
  `ModelMarkov` and `Optimization` unchanged (decisions 009, 015) if possible.
