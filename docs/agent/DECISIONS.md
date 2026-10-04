# Decisions

> **Status, 2026-10-04.** Programming decisions 001 to 019 are recorded; 002 to 014 were
> approved by Peter with the code plan on 2026-09-23, 015 (which supersedes 004) to 018 on
> 2026-10-03, and 019 on 2026-10-04. Ten mathematical and methodological design
> decisions, D01 to D10, are recorded from the design synthesis and are revisable (see that
> section's preamble). The items under "Pending programming candidates" were discussed on
> 2026-09-15; the first three are now resolved by the entries named in their notes. An agent must
> not treat a pending item as a decision.

## What this document is for

This file is primarily the record of **programming decisions** made during development:
load-bearing choices about code, with the reasoning that produced them, so that later work does
not re-litigate them or silently contradict them. It also holds, in a separate section, the
mathematical and methodological decisions from the design stage that the code is built against.
The two kinds are kept apart because they change for different reasons: a programming decision
changes when the code or the codebase demands it, a mathematical decision changes when the
evidence about the model or the optimizer does.

It is **append only**. Entries are never edited or deleted. A decision that turns out to be wrong
is superseded by a new entry that references the old one by its number or ID.

It is deliberately distinct from `CHANGELOG.md`. The changelog is chronological and read
backwards from the end to answer "where did we leave off and what has been tried". This file is
a flat list read at the start of a task to answer "what am I not allowed to quietly change".

If maintaining both proves to be too much friction, drop this file and write decisions into
`CHANGELOG.md` with a `DECISION:` prefix so they stay greppable. Do not keep two homes for the
same content.

## Entry format

```markdown
## NNN. Short imperative title

- **Date:**
- **Status:** accepted | superseded by NNN | reversed
- **Decision:** one or two sentences, stated as a rule.
- **Why:** the reasoning, including the evidence that supports it.
- **Alternatives rejected:** what else was considered and the concrete reason it lost.
- **Affects:** files, commands, or behaviour that depend on this.
```

Number programming entries sequentially from 001 and never reuse a number.

---

## Programming decisions

## 001. Develop on the fork's `nq-constrained-pi` branch, for an eventual upstream pull request

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** All work happens on branch `nq-constrained-pi` of the fork
  `petergoodman/iqtree3`. Nothing is committed to `master` and nothing is pushed to
  `iqtree/iqtree3`. The work is intended to be offered upstream as a pull request from the fork
  after extensive testing and revision, so changes to shared IQ-TREE code are kept opt-in, with
  existing behaviour unchanged by default.
- **Why:** Peter's statement of intent, 2026-09-23. Keeping `master` an exact mirror of upstream
  keeps syncing a fast-forward, and keeping shared-code changes opt-in keeps an eventual pull
  request reviewable.
- **Alternatives rejected:** a permanent private fork with no merge intent, which would permit
  invasive edits to shared code but forfeit contribution upstream.
- **Affects:** every code change; the git procedure in `CLAUDE.md`; how far edits to shared files
  (`utils/optimization.cpp`, `model/partitionmodel.cpp`, `main/phyloanalysis.cpp`,
  `main/phylotesting.cpp`) may go.

## 002. Implement NQC as a `ModelProtein` subclass dispatched from `createModel`

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** The π-constrained model is a new class `ModelNonrevFixedEq`, derived from
  `ModelProtein`, in `model/modelnonrevfixedeq.{h,cpp}`. `createModel` dispatches it by name in
  its protein branch, before `new ModelProtein`, and passes the `ModelsBlock`. Its public name is
  "NQC", provisional until Peter settles the public name with the maintainers, and is held in one
  constant. `NONREV` is not modified.
- **Why:** Tip partial likelihoods come from the model's `computeTipLikelihood` (verified,
  `tree/phylotreesse.cpp:366-373`). `ModelProtein` treats B, Z and J as two-state ambiguities
  (verified, `model/modelprotein.cpp:1350-1366`), while the `ModelSubst` default treats every
  non-single state as fully unknown (verified, `model/modelsubst.cpp:199-209`), so a class outside
  the `ModelProtein` hierarchy would give the same data a different likelihood than `NONREV` or
  `GTR20` and break nesting comparisons. `createModel` holds the `ModelsBlock` needed to read the
  LG, WAG and JTT seeds and `--init-model` names (verified, `model/modelmixture.cpp:3122-3267`),
  whereas `ModelMarkov::getModelByName` receives none (verified, `model/modelmarkov.cpp:1954-1968`).
  A separate name keeps the live `NONREV+F{...}` meaning intact and permits the nested comparison
  with `NONREV`. `ModelProtein`'s constructor calls `init(name)`, which cannot parse a new name
  (verified, `model/modelprotein.cpp:1079-1088`, `1107-1249`), so the subclass constructs its base
  with the seed model's name and converts in its own constructor.
- **Alternatives rejected:** a name branch inside `ModelProtein::init`, which would put constraint
  code into a file that is mostly data tables and need type switches in every `ModelProtein`
  method; a `ModelMarkov` subclass registered through `validModelName` and `getModelByName`, as
  `UNREST` is, rejected for the tip-likelihood and `ModelsBlock` reasons above.
- **Affects:** `model/modelnonrevfixedeq.{h,cpp}`, `model/modelmixture.cpp` (`createModel`),
  `model/CMakeLists.txt`; resolves the pending candidates on a new name and on one class for both
  invocation paths.

## 003. Store the 360 coordinates as authoritative state and checkpoint them losslessly

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** The model keeps its 360 log-ratio coordinates in a separate array that is its
  authoritative state; `rates[]` holds the 380 derived physical rates and `rate_matrix` the
  derived normalized Q. The class checkpoints its own state under its own keys: chart identifier
  and version, target, references, coordinates and target as 17-significant-digit strings,
  domain, derivative policy, fixed or fitted status, and accepted score. Restore verifies target,
  references, chart and domain against the constructed object and rejects any mismatch.
- **Why:** Synthesis section 7.1: the pinned quantity is a jump weight, not a physical rate, so
  every physical rate can change when the coordinates move, and using `rates[]` for both is
  ambiguous; `ModelLieMarkov` already maps a separate parameter vector into `rates[]` (verified,
  `model/modelliemarkov.cpp:899-933`). IQ-TREE writes checkpoint doubles at 10 significant digits
  (verified, `utils/checkpoint.h:22, 352`), and the default +I+G path saves and restores the model
  through its checkpoint (verified, `model/modelfactory.cpp:1395-1505`; `opt_gammai` defaults to
  true, `utils/tools.cpp:7082`), so a 10-digit checkpoint would perturb the model during ordinary
  fits, not only at restart.
- **Alternatives rejected:** coordinates or auxiliary weights inside `rates[]` (P-positive section
  6), superseded by synthesis section 7.1; `CKP_ARRAY_SAVE` at 10 digits, which is lossy; raising
  `CKP_PRECISION` globally, which changes every model's checkpoint and is not opt-in.
- **Affects:** `model/modelnonrevfixedeq.{h,cpp}`; checkpoint and restart tests; +I+G fits.

## 004. Move the numerical half of `decomposeRateMatrixNonrev` into a protected helper

- **Date:** 2026-09-23
- **Status:** superseded by 015
- **Decision:** The eigen-decomposition part of `ModelMarkov::decomposeRateMatrixNonrev`
  (verified, `model/modelmarkov.cpp:1284-1386`) moves unchanged into a protected non-virtual
  helper that the base routine then calls. The NQC class overrides `decomposeRateMatrixNonrev`:
  it builds Q from its coordinates, validates it, installs `rates[]` and `rate_matrix`, and calls
  the helper, never touching `state_freq`.
- **Why:** The base routine resets π to uniform and re-solves it from Q unless the frequency type
  is user-defined and `-optfromgiven` is off (verified, `model/modelmarkov.cpp:1252-1254`,
  `1266-1267`). `optimize_from_given_params` is a process-wide global that cannot be switched per
  object inside the OpenMP loop, so delegating to the unchanged base would let that flag overwrite
  π* with a numerical approximation. Synthesis section 7.2 prefers a narrow shared numerical hook
  over duplicating the eigen and fallback code. Moving the statements unchanged keeps legacy
  behaviour identical, which the regression baseline checks.
- **Alternatives rejected:** writing physical rates into `rates[]` and calling the unchanged base,
  which needs no shared change but is exposed to `-optfromgiven` (kept as the fallback if the
  refactor shows any regression); copying about 100 lines of decomposition code into the class.
- **Affects:** `model/modelmarkov.{h,cpp}`, `model/modelnonrevfixedeq.cpp`; the regression
  baseline.

## 005. Reach the scaled derivative step through a `PartitionModel::derivativeFunk` override

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** The NQC class overrides `derivativeFunk` for single-model optimization, and
  `PartitionModel` gains a `derivativeFunk` override that delegates to a small hook interface
  implemented by the linked model, falling back to `Optimization::derivativeFunk` for every other
  model. The step is h = η max(s, |x|) with η = 1e-4 and s = 1 (D02), exactly representable, and
  taken backward when the forward probe would leave the declared domain; a central stencil exists
  for validation only. The same interface carries the linked-compatibility check and the post-fit
  hook of entry 009.
- **Why:** Linked optimization runs `minimizeMultiDimen` on the `PartitionModel` object (verified,
  `model/partitionmodel.cpp:786`), and `PartitionModel` does not override `derivativeFunk`, so a
  model-class override alone is not reached (verified; the only override outside `Optimization`
  is `tree/phylotreemixlen.h:165`). Delegating from `PartitionModel` leaves
  `utils/optimization.cpp` untouched, and legacy models call the same base routine as before.
- **Alternatives rejected:** a virtual step hook inside `Optimization::derivativeFunk` with a
  legacy default, which is general but edits the routine every optimizer uses; positive
  coordinates to avoid any change, contrary to D01 and D02. The third optimizer entry point,
  `ModelFactory::optimizeAllParameters` under `-jointopt`, hard-codes [MIN_RATE, MAX_RATE] bounds
  on the model's coordinates (verified, `model/modelfactory.cpp:1331-1377`), so it is rejected by
  entry 008 rather than hooked.
- **Affects:** `model/partitionmodel.{h,cpp}`, `model/modelnonrevfixedeq.{h,cpp}`; gradient step
  tests in both paths.

## 006. Take the target through the existing `+F` routes and validate it strictly in the model

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** In the first release π* arrives only as `NQC+F{p1,...,p20}` or as `NQC+F<NAME>`
  with the vector defined in an `--mdef` file. The model parses `freq_params` itself: exactly 20
  finite entries, each greater than 0 and at least `min_state_freq`; a sum within 1e-6 of 1 is
  normalized once and both vectors are reported; anything else is rejected with a message, and
  nothing is floored or substituted. `+FO`, `+F`, `+FQ` and a missing `+F` are rejected with a
  pointer to `NONREV`. Every report and export prints π* at 17 significant digits.
- **Why:** Both routes already reach the model through `ModelFactory` with no shared-code change
  (verified, `model/modelfactory.cpp:489-509` and `599-605`, with the `frequency NAME = ...;`
  syntax at `model/modelmixture.cpp:724`). The existing reader is unsuitable for an immutable
  target: it accepts NaN, turns non-numeric tokens into random draws from AliSim distributions,
  and renormalizes with only a warning (verified, `model/modelmarkov.cpp:1761-1796`,
  `utils/tools.cpp:408-428`). Known limits of the literal route: a `+` or `*` inside a number
  splits the model string (verified, `model/modelfactory.cpp:259, 296`); commas prevent its use in
  `-mset` and `-madd` (verified, `utils/tools.cpp:588`, `main/phylotesting.cpp:1095`); and protein
  model names drop the vector when printed (verified, `model/modelmarkov.cpp:204-206`,
  `model/modelprotein.cpp:1333-1348`). The named route avoids the first two. The 1e-6 sum
  tolerance is a proposed setting that covers vectors printed at 8 significant digits.
- **Alternatives rejected:** a new option such as `--target-freq FILE`, left for Peter to raise
  with the maintainers because it defines public interface; `NQC{...}` braces, which conflate the
  target with model parameters.
- **Affects:** `model/fixedeqchart.cpp` (parser), `model/modelnonrevfixedeq.cpp`; the run
  manifest; the S0 parsing checks.

## 007. Hold domain and step settings as class defaults, with parse-only options

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** The coordinate domain defaults to the D03 historical benchmark domain, z in
  [log 1e-5, log 10], and the step to D02's η = 1e-4 and s = 1. Two parse-only options, declared
  in `utils/tools.h`, defaulted and parsed in `parseArg`, documented in `usage_iqtree()` and read
  at the point of use, let a run set the domain (added in S2) and the step (added when the step
  study needs it). A start not representable with margin is rejected with a message that prints
  the domain required. Option names are provisional until Peter settles public names.
- **Why:** D03 requires explicit and recorded domains, containment of every start, and
  expansion-sensitivity checks, and LG rebuilt at a skewed target can fall outside the default
  domain, so run-time settings are needed from the first compiled slice. Options follow IQ-TREE's
  convention of global settings read through `Params`.
- **Alternatives rejected:** `NQC{...}` braces, which mean initial parameter values elsewhere in
  IQ-TREE and carry commas; compiled-in constants, which cannot be varied for sensitivity checks;
  a domain chosen automatically for each start, which breaks matched comparisons across starts.
- **Affects:** `utils/tools.{h,cpp}`, `model/modelnonrevfixedeq.cpp`; D03 sensitivity runs.

## 008. Reject every model combination the first release does not support

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** NQC stops with an error message for: non-protein data; `+FO`, `+F`, `+FQ` or no
  `+F`; `-optfromgiven` (verified, `utils/tools.cpp:3315-3318`); `-jointopt` (verified,
  `utils/tools.cpp:3381-3382`); mixtures `MIX{...}` and `+FMIX` (verified,
  `model/modelmixture.cpp:3441-3518, 4499`); tree mixtures and HMM models (verified,
  `tree/iqtreemix.cpp:1382, 1612`, `tree/iqtreemixhmm.cpp:246`); site-specific frequency models
  (verified, `model/modelfactory.cpp:643-687`); `--link-exchange-rates`; `--eigen`; ModelFinder
  runs that include NQC; and AliSim simulation from the class, for which the exported fixed matrix
  is used instead. How a mixture context is detected is settled in S2.
- **Why:** Each path bypasses the class's invariants or has not been designed. `-jointopt`
  optimizes the model's coordinates with hard-coded rate bounds and the legacy derivative
  (verified, `model/modelfactory.cpp:1331-1377`), which would clamp log-ratios at 1e-4 or above;
  `-optfromgiven` releases the equilibrium in the base decomposition (entry 004); the `--eigen`
  path is already broken for non-reversible models (reported, `AA_MODEL_INFERENCE.md` section 15;
  the base branch verified at `model/modelmarkov.cpp:1284-1289`); the others need their own
  designs (synthesis section 2, out-of-scope list).
- **Alternatives rejected:** allowing them silently, unsafe for the reasons above; supporting them
  all in the first release, beyond its scope.
- **Affects:** `model/modelnonrevfixedeq.cpp`; rejection tests in S2.

## 009. Surface failed line searches by post-fit checks first, and record status later

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** From S2, the NQC paths re-evaluate the returned point after each fit, check the
  invariants, record bound activity and restore the incumbent if the score fell, and at export
  compute a central-difference first-order diagnostic in log-ratio coordinates; the stop is
  labelled "unclassified" unless that diagnostic passes. In S5, `dfpmin` and `lnsrch` gain
  recording-only fields for line-search failures and the stop reason, with no change to any
  trajectory. The line search itself is not repaired unless Peter relays maintainer approval of
  that broader change.
- **Why:** When the step falls below its minimum, `lnsrch` restores the old point, sets
  `check = 1` and leaves `*f` at the last trial (verified, `utils/optimization.cpp:687-690`);
  `dfpmin` never reads `check`, and the zero displacement then passes its small-step test
  (verified, `utils/optimization.cpp:832-846`). A return therefore does not certify convergence
  (synthesis section 5). Post-fit checks need no change to the shared optimizer, and exact stop
  reasons matter first for the G4 comparisons.
- **Alternatives rejected:** status propagation from the start, which edits
  `utils/optimization.cpp` before it is needed; repairing the line search, which changes legacy
  trajectories.
- **Affects:** `model/modelnonrevfixedeq.cpp`, the `PartitionModel` post-fit hook, and in S5
  `utils/optimization.{h,cpp}`.

## 010. Export the fitted matrix at 17 digits in the `-m FILE` format, with a metadata file

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** After the final model optimization, a hook in `main/phyloanalysis.cpp`, active
  only when an NQC model exists, writes `<prefix>.NQC.qmat`, holding the 20 rows of the normalized
  Q and the π* row at 17 significant digits, and `<prefix>.NQC.info`, holding the target, state
  order, chart and version, references, domain, derivative policy, residuals, bound activity,
  stop classification, seed identity and source commit. Re-import uses the existing `-m FILE`
  path. The `.iqtree` report block is left as it is.
- **Why:** The report prints the matrix at 6 significant digits (verified,
  `main/phyloanalysis.cpp:633`), `.best_model.nex` omits the rates, and re-import from the 6-digit
  block reproduced the log-likelihood only to 4e-4 (reported, `CHANGELOG.md`, 2026-09-23). The
  existing reader loads a full-Q file whose first entry is negative as a non-reversible model and
  keeps the file's π as the root distribution (verified, `model/modelmarkov.cpp:1809-1812`,
  `model/modelprotein.cpp:1234-1248`), so no reader change is needed.
- **Alternatives rejected:** raising the report block's precision, which edits shared report code
  for every model; a NEXUS block through `reportNexusFile`, which writes 6 digits, labels every
  model `GTRPMIX` and writes a uniform frequency line (reported, `ARCHITECTURE.md` section 12).
- **Affects:** `main/phyloanalysis.cpp`, `model/modelnonrevfixedeq.cpp`; export and re-import
  tests.

## 011. Recompute reference destinations deterministically in every model object

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** Every NQC object computes its reference destinations from LG rebuilt at π*, read
  from the built-in models block, as the row maxima with ties going to the lower amino-acid index,
  whatever start it is given. The references are stored in the checkpoint, the export metadata
  and the report, and verified on restore and across linked objects before optimization.
  User-supplied references are deferred.
- **Why:** D01 and synthesis section 3.4 fix the rule and require one shared definition for every
  linked object and every start. A deterministic function of π* gives every object the same
  references without shared mutable state, which would be unsafe inside the OpenMP loop.
- **Alternatives rejected:** copying references from a first object, which introduces shared
  state and ordering dependence; user-supplied references in the first release.
- **Affects:** `model/fixedeqchart.cpp`, `model/modelnonrevfixedeq.cpp`; linked-compatibility and
  restart tests.

## 012. Test the pure module with googletest in a standalone CMake project

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** C++ unit tests live in `unittest/`, a standalone CMake project that compiles
  `model/fixedeqchart.cpp` directly and fetches googletest at the commit cmaple pins (verified,
  `cmaple/CMakeLists.txt:281-293`). The class's behaviour is tested through the binary. Fixtures
  are plain text at 17 significant digits, written by the oracle with a provenance header, and
  never regenerated to make a test pass. Integration into the main build is revisited only if
  class-level C++ tests become necessary.
- **Why:** The pure module depends only on Eigen, so it needs no IQ-TREE libraries, and a
  standalone project changes no shared build file. cmaple's googletest targets do not exist on
  Windows, where `USE_CMAPLE` is not an option (verified, `CMakeLists.txt:255-258`), and depending
  on them ties the tests to a vendored subproject.
- **Alternatives rejected:** reusing cmaple's targets inside the main build; an opt-in test
  subdirectory in the main build, which needs a guarded edit to the root `CMakeLists.txt` that
  nothing yet requires.
- **Affects:** `unittest/`, `.github/workflows/fixedeq.yaml`.

## 013. Keep the Python oracle in `test_scripts/fixedeq/`, in a conda environment in WSL

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** The oracle, drivers and fixture generator live in `test_scripts/fixedeq/`. They
  run in a conda environment named `iqtree3-fixedeq`, created in WSL from
  `test_scripts/fixedeq/environment.yml` (conda-forge only, with pinned python, numpy, scipy and
  pytest), with a lock file exported after creation. Because `conda` is not on the PATH of shells
  started by `wsl.exe`, commands run as `~/anaconda3/bin/conda run -n iqtree3-fixedeq ...`. The
  system Python is never used or modified.
- **Why:** `docs/agent/design/scripts/` must stay unedited reference copies (design README).
  `test_scripts/` holds IQ-TREE's test harness, and a new subdirectory there adds files without
  conflicting with upstream. The oracle must run beside the IQ-TREE binary in WSL for the
  differential tests. Checked read-only on 2026-09-23: `~/anaconda3` exists with the libmamba
  solver and the `defaults` channel, `conda` is not on the PATH in `wsl.exe` shells, and the
  system `python3` is 3.11.2. Conda-forge avoids depending on the `defaults` channel.
- **Alternatives rejected:** placing derived code beside the design scripts; a top-level `oracle/`
  directory; a virtual environment built on the system Python.
- **Affects:** `test_scripts/fixedeq/`; every oracle, fixture, regression and differential
  command.

## 014. Write drivers in Python and scope line-ending rules to the new directories

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** Test and regression drivers are Python files run as `python file.py`, and one-off
  command scripts for `wsl.exe` are written with LF line endings. The new directories carry their
  own `.gitattributes` (`*.sh text eol=lf` and `*.py text eol=lf`); no root `.gitattributes` is
  added and `core.autocrlf` is not changed.
- **Why:** This checkout has `core.autocrlf=true` and no `.gitattributes` (verified), and shell
  scripts with CRLF endings fail under WSL bash, while Python reads either ending. Scoped
  attribute files fix the new directories without changing how upstream's scripts check out in
  this working copy.
- **Alternatives rejected:** a root `.gitattributes`, which also affects upstream's files;
  changing `core.autocrlf`, a machine-level setting outside the repository.
- **Affects:** `test_scripts/fixedeq/`, `unittest/`.

## 015. Build Q in the class's decomposition override and call the unchanged base decomposition

- **Date:** 2026-10-03
- **Status:** accepted
- **Decision:** Supersedes 004. The NQC class's `decomposeRateMatrixNonrev` override builds Q from
  the stored coordinates (entry 003) into per-object workspace, validates it, writes its 380
  off-diagonal rates into `rates[]`, and calls the unchanged `ModelMarkov::decomposeRateMatrixNonrev`.
  After that call it checks that `state_freq` is bitwise equal to π* and stops with an error if it
  is not. A failed build writes nothing and sets the flag that makes the class's `targetFunk`
  return 1e30. `setVariables` and `getVariables` only copy the coordinates out of and into the
  model; they do not build Q. `ModelMarkov` is not edited.
- **Why:** When the frequency type is `FREQ_USER_DEFINED` and `-optfromgiven` is off, the base
  routine skips the π reset and the stationarity solve (verified, `model/modelmarkov.cpp:1252,
  1266`). Entry 008 rejects `-optfromgiven` for NQC, and the flag is set only from the command line
  (verified, `utils/tools.cpp:3316`, default false at 7284), so the reason given for 004 no longer
  holds. The base routine then unpacks `rates[]` into `rate_matrix`, rescales by
  total_num_subst / Σ π*_i (−q_ii), a factor equal to 1 up to rounding because the built Q has unit
  mean rate at π*, and eigendecomposes. Building Q in the decomposition override rather than in
  `getVariables` makes every eigensystem come from the coordinates, because every decomposition
  passes through that override. In particular `ModelProtein::restoreCheckpoint` writes `rates[]`
  from the checkpoint at 10 significant digits and then decomposes (verified,
  `model/modelprotein.cpp:1265-1276`; `utils/checkpoint.h:22`), and the default +I+G optimization
  restores the model on each restart and at the end (verified, `model/modelfactory.cpp:1504,
  1550`); that decomposition is dispatched to the override, so the 10-digit rates never become Q.
  The IQ-TREE maintainers advised caution with changes to `ModelMarkov` (reported by Peter,
  2026-10-01). The check after the base call turns any later change to the base routine's π
  guards into an error rather than a silent change of the target.
- **Alternatives rejected:** 004's helper split, which leaves Q bitwise as built but edits the
  routine every non-reversible model uses, against the maintainers' advice and with no remaining
  reason. Building Q in `getVariables` with the unchanged base decomposition, the placement the
  maintainers described and the one `ModelLieMarkov` uses through `setRates` (verified,
  `model/modelliemarkov.cpp:919-933, 615-625`): the base decomposition would then trust `rates[]`,
  so the checkpoint restore above, or any later writer of `rates[]`, would set Q until the next
  changed coordinate vector, and the `changed` test in `targetFunk` (verified,
  `model/modelmarkov.cpp:1092-1100`) would leave that Q as the base point of the next
  finite-difference gradient.
- **Affects:** `model/modelnonrevfixedeq.{h,cpp}`; `model/modelmarkov.{h,cpp}` are no longer
  changed; entry 008's rejection of `-optfromgiven` becomes load-bearing; S2 tests (the built Q
  compared with the oracle after the base rescaling, target immutability).

## 016. Root the tree by IQ-TREE's conversion and do not search the root in S0 to S4

- **Date:** 2026-10-03
- **Status:** accepted
- **Decision:** Every NQC run and every comparison run in S0 to S4 follows four root rules, stated
  separately as synthesis section 7.5 requires. Root frequencies: π*. Input rooting: an unrooted
  input tree is rooted by IQ-TREE's own conversion, at the midpoint of the longest path, or on
  the pendant branch of the `-o` outgroup when one is given; an already rooted input tree keeps
  its root. Root-edge search: none; fits run on a fixed topology under `-te`, with no tree search
  and without `--root-find`. Root split: the two branch lengths beside the root are fixed in
  Level 1 fits (`-blfix`) and left to IQ-TREE's ordinary branch-length optimization in Level 2
  fits. S0 probe (f) checks at run time that the root edge does not move and records what happens
  to the split. The root policy for scientific runs is set at gate G5.
- **Why:** Peter's choice on 2026-10-03, following the design's default
  (`design/constrained_nq_plan_agent.md`, root-position row: inherit nQMaker, whose root position
  is not moved under `-te` during step 3). A non-reversible likelihood depends on the root, so the
  S2 and S3 fits and the layer 3 and 5 comparisons are reproducible only under a stated policy,
  and the G0 manifest records it. In source: `convertToRooted` places the root as described
  (verified, `tree/phylotree.cpp:5906` onward); the root edge is changed only by
  `optimizeRootPosition`, which is called inside NNI search when the tree is rooted and
  `root_move_dist` is positive (verified, `tree/iqtree.cpp:3257-3259`; default 2,
  `utils/tools.cpp:7117`) and after model optimization when `--root-find` is on (verified,
  `model/modelfactory.cpp:1733-1734`; `root_find` defaults to false, `utils/tools.cpp:7118`);
  `-te` sets `min_iterations = 0` and a fixed-iteration stop (verified,
  `utils/tools.cpp:4699-4706`). That a `-te` run never reaches the NNI call has been read, not
  run; probe (f) settles it.
- **Alternatives rejected:** searching the root during fits, by tree search or `--root-find`,
  which adds a discrete nuisance that comparisons on matched rooted trees cannot absorb (deferred
  to G5); requiring an outgroup, which the test data do not supply; fixing the root split in
  Level 2 fits as well, which departs from the native branch-length optimization that synthesis
  section 7.5 asks ordinary fits to preserve.
- **Affects:** S0 probe (f) and the G0 manifest; the Level 1 and Level 2 fits of S2 and S3; test
  layers 3 and 5; the trees given to `-te` in every NQC test.

## 017. Port IQ-TREE's optimizer into the oracle from the C++ source

- **Date:** 2026-10-03
- **Status:** accepted
- **Decision:** The oracle's optimizer, against which the S2 toy comparison runs the compiled
  one, is a port written from `utils/optimization.cpp` at `63c330d9`, not copied from
  `design/scripts/chart_fd_bfgs_compare.py`. It reproduces `minimizeMultiDimen`, `dfpmin`,
  `lnsrch`, `fixBound` and the legacy forward-difference `derivativeFunk` with IQ-TREE's
  constants and control flow, including the failed-line-search path: `lnsrch` restores the old
  point, sets `check` and leaves the objective at the last trial value, and `dfpmin` never reads
  `check`. It is called with the arguments `ModelMarkov::optimizeParameters` passes. The scaled
  step of entry 005 is a separate function beside the legacy rule, and the stationary solve used
  in the toy comparison is a column-pivoted QR, as in the C++ (`CODE_PLAN.md` section 2.1). The
  match is in algorithm, constants and control flow; floating-point identity is not expected,
  because NumPy and Eigen order sums differently and use different exponential and logarithm
  routines. Agreement on the toy problem is checked against a tolerance fixed before that test is
  written. The design script stays unedited.
- **Why:** Peter's choice on 2026-10-03: the oracle's copy should match IQ-TREE. The design port
  returns the old objective on a failed line search where the C++ leaves the last trial value
  (synthesis section 5), and it solves the stationary system by NumPy least squares rather than
  Eigen's column-pivoted QR (synthesis section 10.5), so a disagreement between it and the
  compiled optimizer could come from the port rather than from the new code. Anchors (verified
  2026-09-23, `CODE_PLAN.md` section 5): `fixBound`, `utils/optimization.cpp:149-156`; `lnsrch`,
  645-718, failure path 687-690; `minimizeMultiDimen`, 750-778; `dfpmin`, 793-900, small-step
  return 843-846; legacy `derivativeFunk`, 916-939; the single-model call,
  `model/modelmarkov.cpp:1199`.
- **Alternatives rejected:** reusing the design port unchanged, for the two differences above;
  editing the design script, which must stay an unedited reference (design `README.md`;
  entry 013).
- **Affects:** `test_scripts/fixedeq/oracle/`; `CODE_PLAN.md` section 3.1; the S2 test that
  compares the compiled optimizer with the port.

## 018. Judge baseline run 11 by its fixed parts, and add run 12, run 10 at four threads

- **Date:** 2026-10-03
- **Status:** accepted
- **Decision:** In regression comparisons, run 11 (`--model-joint NONREV` under `-p` at `-T 4`)
  is judged only on the parts that did not vary across its five baseline repeats: the exit code;
  the REFERENCES and SEQUENCE ALIGNMENT sections, exactly; and the SUBSTITUTION PROCESS section
  and the tree file with every number removed, so that the wording and the topology must match
  while the values may differ. Its "MAXIMUM LIKELIHOOD TREE" section is not compared, because its
  ASCII drawing changes with the branch lengths and the topology is checked through the tree
  file. A run 12, run 10 at `-T 4`, joins the baseline list; it is recorded five times from the
  frozen binary `~/iqtree3-baseline/iqtree3` (SHA-256 `9a72950b…c04fb501`) and compared under the
  general rules of `CODE_PLAN.md` section 3.4. Both are implemented in the regression driver
  before S2.
- **Why:** Peter's choice on 2026-10-03 (PLAN.md risk 19). The five repeats of run 11 ended at
  log-likelihoods from -4975.0215 to -4974.5963, none equal to the `-T 1` run's -4974.5432, so
  its numbers carry no regression signal; its wording and topology were identical in every
  repeat and still catch a crash, a changed model description or a changed tree. Run 12 fixes
  the tree and the matrix, so it exercises the multi-threaded likelihood code with far less room
  for the optimizer's path to differ; that it reproduces exactly is expected, not yet observed.
  The recommendation as first worded, "the same report text apart from the numbers", would have
  failed on the tree drawing; this entry states the precise rule.
- **Alternatives rejected:** holding run 11's numbers to their observed spread, which a genuine
  repeat would often fail; dropping run 11, which would leave the threaded joint fit unchecked.
- **Affects:** `test_scripts/fixedeq/regression/regress.py` and its tests; the baseline in
  `test_scripts/fixedeq/regression/baseline/`; `CODE_PLAN.md` section 3.4.

## 019. Rerun the design scripts with a kept driver and match documented values in four classes

- **Date:** 2026-10-04
- **Status:** accepted
- **Decision:** The three scripts in `docs/agent/design/scripts/` are rerun by the kept driver
  `test_scripts/fixedeq/design_rerun.py`. It copies them into a new directory outside the
  repository, checks each copy against the SHA-256 in the design `README.md`, runs them with the
  arguments and environment the documents state, and compares their outputs with a table of the
  values the documents print, each cited by document and line. Each documented value has one
  class. **Exact** (counts, ranks, stop labels, exit status, failure lists): the rerun value must
  be equal. **Digits** (any other number): the rerun value, rounded to the number of significant
  digits printed, must equal the printed value. **Zero** (a quantity that is mathematically zero,
  printed below 1e-12): reproduced when the rerun's magnitude is also below 1e-12, with both
  values listed. **Bound** (a claim printed as an inequality): the rerun value must satisfy it.
  An entry whose mapping from the document's prose to an output key needed a reading is marked
  as interpreted. The classes and the 1e-12 cut were fixed before any script ran and are not
  changed after seeing results.
- **Why:** Peter's choice on 2026-10-04 (PLAN.md risk 10). Of the 133 documented values, 24 are
  rounding noise around zero, whose digits depend on the processor and the linear-algebra
  library, so requiring their printed digits would test the platform rather than the
  mathematics; every other value is held to the precision the document printed. A driver under
  version control makes the rerun reproducible, and the oracle's layer 1 test (`CODE_PLAN.md`
  section 3.1) can reuse it.
- **Alternatives rejected:** printed digits for every value, which would count platform noise as
  a mismatch; one-off scripts kept beside the outputs outside the repository, which would leave
  the procedure outside version control.
- **Affects:** `test_scripts/fixedeq/design_rerun.py` and its tests; `CODE_PLAN.md` sections 1.1
  and 3.1; PLAN.md risk 10; the G0 run manifest.

---

## Mathematical and methodological decisions (D01 to D10)

These entries carry the IDs D01 to D10 of the design synthesis,
`docs/agent/design/project2_unified_discrepancy_synthesis.md` (its section 1), so that they can
be cross-referenced with it. They were adopted by Peter on 2026-09-23 as the design the code is
built against. **They are mathematical and methodological decisions, not programming decisions,
and they are subject to change.** A change follows the synthesis's change-control rule (its
section 13): name the affected ID, the evidence, whether the change affects the estimand or only
the implementation, the regression implications, and whether finished comparisons must be rerun.
It is recorded as a new entry that references the D-number, not by editing the D entry. How each
decision is realised in code is a programming decision, recorded as a numbered entry above.

The family throughout is the set of 20-state generators with positive off-diagonal rates, zero
row sums, π*Q = 0, and unit mean rate at π*: 360 free parameters, containing the 189-dimensional
reversible family at π*.

## D01. Parameterize the family by jump-chain log-ratios first

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** The first implemented chart is the jump chain in log-ratio coordinates: per row,
  18 log-ratios against a fixed reference destination, a stable row softmax giving the jump
  matrix K, one stationary solve νK = ν, and exit rates ν_i/π*_i (360 coordinates in all).
  Reference destinations are chosen once, from the row maxima of the target-rebuilt LG seed with
  ties broken in amino-acid order, and are shared by every linked object and every start.
- **Why:** Both log-ratios and positive weights cover the whole positive family one-to-one
  (synthesis section 3.1). Log coordinates make steps multiplicative in the weights. In the
  synthesis's rerun of four eight-state conditional fits, log coordinates scored higher in three,
  and the three positive-weight runs that stopped on `TOLX` coincided with failed line searches
  (synthesis sections 10.3 and 10.4). This is evidence for a starting choice, not a theorem.
- **Alternatives retained:** positive weights with one pinned weight per row (c = 10), kept as a
  maintained benchmark. Reopen if matched compiled tests show them more reliable or efficient at
  comparable solution quality.
- **Source:** synthesis section 3; P-log section 3; P-positive section 3.1.

## D02. Use scale-aware finite differences from the first fitted log-coordinate implementation

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** The first fitted log-coordinate implementation includes an opt-in derivative step
  h_k = η max(s_k, |x_k|), with initial s_k = 1 and η = 1e-4 (a reproducible starting setting,
  not a validated constant), reachable from both the single-model optimizer path and the linked
  `PartitionModel` path. Existing models keep the legacy step by default.
- **Why:** IQ-TREE's step is 1e-4 |x_k|, which collapses near zero for signed coordinates. In the
  synthesis's rerun sweep, the relative derivative discrepancy of one logit was about 9.6e-2 at
  z = 1e-6 and 6.4 at z = 1e-8 (synthesis section 10.2). In linked training the optimizer runs on
  the `PartitionModel` object, so an override in the model class alone is not reached
  (P-positive section 4.6).
- **Alternatives retained:** the legacy relative rule as a diagnostic control; analytic
  likelihood gradients later.
- **Source:** synthesis section 4.1; P-log section 6.2.

## D03. Declare numerical domains explicitly and report their effect

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** Coordinate bounds are explicit and recorded for every run. The historical
  matched-benchmark domain is z in [log 1e-5, log 10], about [-11.513, 2.303], which corresponds
  to positive weights in [1e-4, 100] with c = 10; it is a benchmark domain, not a chosen
  production range. Every declared start must be representable with margin, bound activity is
  reported, final scientific runs include a domain-expansion sensitivity check, and an imported
  matrix is never clipped to fit a box.
- **Why:** A finite box restricts the positive family, and boxes in different charts are
  different sets of matrices (synthesis section 4.2; P-positive Theorem 6, remark 2).
- **Alternatives retained:** expanding or revising the limits when bounds are active or
  sensitivity checks move the solution; a boundary-capable method if exact zero rates are ever
  required (D10).
- **Source:** synthesis section 4.2.

## D04. Implement T3 balancing as a second compiled backend

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** After the log-jump core is validated, the T3 balancing chart is implemented in
  IQ-TREE on the same likelihood engine and derivative interface, as a serious comparison
  method. Its gauge is a_{19,20} = 0 and h_{i,20} = 0 for i < 20 (189 strength and 171 tilt
  coordinates), with the 19 balancing potentials found by damped Newton. A Python version is a
  correctness reference and does not by itself meet the comparison goal.
- **Why:** T3 and the jump chart parameterize the same family, so attained differences reflect
  optimization, bounds, or implementation, which is what the comparison is meant to measure
  (synthesis section 6.1). Its production exposure depends on benchmarks and maintainer
  preference.
- **Alternatives retained:** P-positive's star-at-state-1 gauge, which describes the same family
  and is imported only by explicit conversion.
- **Source:** synthesis section 6.1; P-log section 4.

## D05. Build the first seed and the guide candidates at the target

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** The default first fitting start is LG rebuilt at the target, q_ij = R_ij π*_j,
  normalized to unit mean rate. The fitted reversible model at the target (`GTR20+F{π*}`) is an
  additional start and the nesting incumbent. The primary guide-candidate protocol uses LG, WAG
  and JTT rebuilt at π*, with empirical and optimized frequency variants suppressed.
- **Why:** One composition convention holds across every stage. Guide candidates only select
  per-alignment models and trees, so ordinary candidates would not invalidate the constrained
  fit (synthesis section 6.2; P-log section 7.3). Row-normalizing LG at LG's own frequencies
  gives a valid but different start, whose exchangeabilities are not LG's.
- **Alternatives retained:** ordinary empirical guide candidates, and flux-transported or
  T3-transported starts, as sensitivity and multistart arms.
- **Source:** synthesis sections 6.2 and 6.3; P-positive section 3.1.5.

## D06. Change only the map from coordinates to generator, and protect its invariant

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** IQ-TREE's likelihood engine, transition matrices, pruning, and branch, rate and
  tree algorithms are reused. The substantive change is the map from 360 optimizer variables to
  a valid Q in the shared-matrix update (nQMaker step 3a), with the root distribution equal to
  π*. The fixed target, chart metadata and reference destinations must survive every part of the
  model's lifecycle (construction, linking, checkpoint, export), and linked objects must agree on
  them before optimization. Outer candidate reselection is orchestrated explicitly rather than
  assumed to work through one ModelFinder call.
- **Why:** Synthesis section 7. The inspected ModelFinder code rejects candidate sets that mix
  recognized reversible and nonreversible models (`mixRevNonrev`, `main/phylotesting.cpp`).
- **Alternatives retained:** a thin external driver for the outer loop first, native integration
  later.
- **Source:** synthesis section 7; P-positive section 5; P-log sections 7.3 and 8.

## D07. Take one explicit, immutable target and never floor it silently

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** The fitter consumes one explicit, strictly positive, normalized π* in IQ-TREE's
  state order and never re-estimates it. The default provenance is the pooled composition of the
  original training data through IQ-TREE's estimator, computed once and held fixed across
  cleaning treatments, with the estimator, ambiguity and gap treatment, any smoothing, and the
  exact numbers recorded. An estimator floor is a declared procedural choice; the fitter does not
  floor supplied input, and it rejects values its numerical backend cannot support.
- **Why:** Synthesis section 8. `ModelMarkov::targetFunk` rejects positive frequencies below
  `min_state_freq` (default 1e-4, configurable with `--min-freq`), and the decomposition removes
  states whose frequency is below `ZERO_FREQ` = 1e-10, so a small target entry must be checked
  against both rather than silently changed.
- **Alternatives retained:** other declared estimators and supported small positive targets.
- **Source:** synthesis section 8; P-positive section 7; P-log section 1.4.

## D08. Plan reversible, unrestricted and transport comparison arms

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** The scientific comparison includes original-data reversible and unrestricted
  fits, cleaned-data reversible and unrestricted fits, the cleaned reversible fit followed by
  frequency replacement, the reversible refit at π*, the general refit at π* (the Project 2
  estimate), and the cleaned unrestricted fit transported to π* by fixed flux and by T3 without
  refitting. Frozen matrices are evaluated on the same held-out data. The transport arms do not
  block the first implementation.
- **Why:** Since π* e^{tQ} = π*, every tip marginal is fixed, and a wrong target can produce
  apparent circulation; simulations with reversible and nonreversible generators at correct and
  perturbed targets are needed before circulation is read as biology (synthesis section 9.1).
- **Alternatives retained:** none dropped; the arms are the comparison plan.
- **Source:** synthesis section 9.1; P-positive sections 8 and 9; P-log section 10.3.

## D09. Adopt the corrected mathematical qualifications

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** These statements override conflicting prose in the older documents: (1) with raw
  off-diagonal rates in [l, u], the gauge-fixed T3 triangle tilt satisfies
  |h'| <= (3/2) log(u/l), which is 20.72 for [1e-4, 100], so a bound of 12 does not contain that
  box; (2) q_ij <= 1/(2 π*_i); (3) global suprema over nested domains are ordered, but
  independently attained local fits need not be, so the smaller model's fit and nuisances are
  imported as an incumbent; (4) an interior true parameter does not imply an interior sample
  maximum; (5) the root-split singularity under reversibility means the 171 and 19 parameter
  differences do not automatically give chi-square reference distributions; (6) replacing
  frequencies in asymmetric factors can still give a valid generator, and what is lost is target
  stationarity; fixing root frequencies alone is not the constraint; (7) no claim of a global
  optimum is made, and the BFGS inverse Hessian is not a covariance estimate.
- **Why:** Synthesis section 9.2, with derivations in P-positive section 2 and P-log sections 4.5,
  5 and 10.2.
- **Alternatives retained:** not applicable.
- **Source:** synthesis section 9.2.

## D10. Keep direct balanced flux as an independent reference and defer EM

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** Direct balanced-flux fitting, f = f0 + Zy with a feasible constrained solver, is
  kept as a small-problem reference and a possible boundary-capable backend. Constrained
  substitution-history EM is deferred; if implemented, the normalized M step (19 balance and one
  normalization multiplier) is preferred. Penalty and projection methods do not substitute for
  the constrained maximum-likelihood estimate.
- **Why:** Synthesis section 9.3.
- **Alternatives retained:** EM, reopened for persistent direct-optimization failures, a need for
  exact boundary fits, or likelihood-gradient infrastructure.
- **Source:** synthesis section 9.3; P-positive section 3.3; P-log sections 5.3 and 9.

---

## Pending programming candidates

These need Peter's confirmation before they become numbered entries.

### Candidate: give the constrained model a new name rather than modifying `NONREV`

Supporting evidence gathered so far: `-m NONREV+F{...}` is already a live code path that pins
`state_freq` to the supplied vector and skips the stationarity solve in
`decomposeRateMatrixNonrev`, leaving Q unconstrained, so overloading that string would change
the meaning of commands that have already been run. The two models are nested at 360 and 379
free parameters, so separate names permit a likelihood-ratio test between them. Parameter
counts, bounds, and checkpoint payloads differ, so one name would mean branching inside every
one of those. `GTR20` and `NONREV` are already separate name branches in the same class, which
is direct precedent.

Open sub-question: whether to implement as a new `ModelMarkov` subclass registered through
`ModelMarkov::getModelByName`, following `ModelUnrest`, or as another name branch inside
`ModelProtein::init`. A subclass keeps the constraint machinery generic over state count and out
of a file that is mostly data tables.

> Note, 2026-09-23: both design plans also propose a new name (`NQC` provisionally), and D09(5)
> qualifies the likelihood-ratio statement above.

> Resolved, 2026-09-23: entry 002 (a new class, `ModelNonrevFixedEq`, provisional name NQC;
> `NONREV` unchanged).

### Candidate: support single-alignment and multi-partition paths through one subclass

Supporting evidence: `PartitionModel::optimizeLinkedModel` is parameterization-agnostic. It
calls `setVariables`, `setBounds`, and `targetFunk` on the model and never inspects what the
parameters mean, so a correct `ModelMarkov` subclass works under both `-m <name>` and
`--model-joint <name>` with no additional code. Two conditions attach: the parameterization must
be a pure function of the variable vector plus fixed per-model constants, with no hidden
per-partition state; and one shared Q across partitions implies one shared π, so per-partition π
with a shared constrained Q would be a different feature.

> Note, 2026-09-23: "no additional code" is contradicted by D02 and D06. The linked path computes
> gradients through `Optimization::derivativeFunk` on the `PartitionModel` object, which no model
> class can override; the `PartitionModel` constructor calls `adaptStateFrequency` for linked
> models with `FREQ_ESTIMATE` or `FREQ_EMPIRICAL` (`model/partitionmodel.cpp:116-165`); and
> linked objects need an explicit compatibility check.

> Resolved, 2026-09-23: entries 002, 005 and 009 (one class serves both paths, with the
> linked-path hooks in `PartitionModel`); the supported invocations are in `PLAN.md`, "Scope of
> the first release".

### Candidate: follow the Lie-Markov reduced-parameterization pattern

Supporting evidence: `ModelLieMarkov::setBasis` already implements this exact feature for
4-state DNA, including reducing the free parameter count by the frequency degrees of freedom and
warning when a requested π is unreachable. See `docs/agent/ARCHITECTURE.md` section 9 for what
transfers to 20 states and what does not.

> Note, 2026-09-23: superseded in substance by D01, which fixes the chart. The "what transfers"
> passage of `ARCHITECTURE.md` section 9 has been removed as design advice. `ModelLieMarkov`
> remains relevant only as a precedent for storing coordinates separately from `rates[]` and
> mapping them in `getVariables` (synthesis section 7.1).

> Resolved, 2026-09-23: D01 and entry 003.

### Candidate: baseline and branching strategy

Which upstream commit the work is based on, and how the fork is structured. Not yet decided.
Record the chosen baseline commit here once fixed, because the anchors in
`docs/agent/ARCHITECTURE.md` and `docs/agent/FILE_INDEX.md` depend on it.

> Note, 2026-09-23: the branching half is now entry 001. The baseline is currently `63c330d9`
> (tag `v3.1.4`), the commit the branch and the documentation anchors are based on; a later sync
> would move it by the procedure in `CLAUDE.md`.
