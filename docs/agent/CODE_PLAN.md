# Code plan: π-constrained non-reversible amino-acid models

> **Read `docs/agent/PLAN.md` first.** This document is its subordinate companion and is never
> read on its own. It holds only the file-level change map, the specifications of the new code,
> the test specification, and the per-slice file and test lists for the slices PLAN.md defines.
> It records no decisions: every choice it relies on is a numbered entry in
> `docs/agent/DECISIONS.md`, cited here by number, and its options and reasoning are stated only
> there. Where this document disagrees with PLAN.md or DECISIONS.md, those win and this document
> is corrected. It is overwritten as the plan moves; history goes in `CHANGELOG.md`.

Labels are those of PLAN.md. Source anchors refer to `63c330d9` and were verified on 2026-09-23
(section 5).

## 1. Change map

### 1.1 New files

| File | Slice | Content |
|---|---|---|
| `model/fixedeqchart.{h,cpp}` | S1 | Pure functions over Eigen types, with no IQ-TREE headers and no process exit (section 2.1) |
| `model/modelnonrevfixedeq.{h,cpp}` | S2 to S4 | Class `ModelNonrevFixedEq : public ModelProtein`, public name "NQC" held in one constant (decision 002; section 2.2) |
| `unittest/` | S1 | Standalone CMake project that compiles `model/fixedeqchart.cpp` directly and fetches googletest at cmaple's pinned commit (verified, `cmaple/CMakeLists.txt:281-293`; decision 012); `fixtures/` in plain text at 17 significant digits with a provenance header; its own `.gitattributes` (decision 014) |
| `test_scripts/fixedeq/` | S0 onward | `environment.yml` and its exported lock (decision 013); the `oracle/` package; `tests/`; `regression/` and `differential/` drivers; `design_rerun.py`, which reruns the design scripts (decision 019); `probes/`, the S0 runtime probes; `sanitizer/`, the decision 022 build and run driver; `manifest/`, the G0 manifest generator and its output; `explore/`, the exploratory measurements cited by decisions 026 and 027; `make_fixtures.py`; its own `.gitattributes`, and a `.gitignore` for Python caches |
| `.github/workflows/fixedeq.yaml` | S1 | Fork-only workflow that builds and runs the unit tests and the oracle's tests on Linux, and from S2 a short NQC run |

Run outputs never go into the working copy; they go under `~/iqtree3-runs/` in WSL. Fixtures and
the summarized regression baseline are derived data kept in the repository with provenance
headers, each produced by a kept script, and are never regenerated to make a test pass.

### 1.2 Existing files touched

| File and anchor (verified) | Slice | Change (proposed) | Why the file must change | How it stays opt-in or behaviour-preserving |
|---|---|---|---|---|
| `model/CMakeLists.txt` | S2 | add `modelnonrevfixedeq` and `fixedeqchart` | the source list is explicit, with no glob | additive |
| `model/modelmixture.cpp`, `createModel`, 3247-3257 | S2 | in the protein branch, dispatch the NQC name to `new ModelNonrevFixedEq(..., models_block)` before `new ModelProtein` (decision 002) | the dispatcher that holds a `ModelsBlock` | reached only by the NQC name |
| `utils/tools.{h,cpp}`, `parseArg` and `usage_iqtree` | S2 (domain), S2 or S3 (step) | two parse-only options, read at the point of use through `Params` (decision 007) | domain expansion and step studies need run-time settings | unused unless NQC is used |
| `model/partitionmodel.{h,cpp}`, `optimizeLinkedModel`, 753-840 | S3 | a `derivativeFunk` override that delegates to the linked model's hook interface when present and otherwise calls `Optimization::derivativeFunk`; a compatibility check before optimizing; a post-fit hook that reports after (decisions 005, 009, 029) | linked optimization runs `minimizeMultiDimen` on the `PartitionModel` object (786), which no model override reaches | legacy models call the unchanged base routine; `partitionmodel.h:160-165` already carries an orphan comment for such an override |
| `main/phyloanalysis.cpp`, after final optimization near 3903-3907; citation block 164-179 | S4 | write the export files when an NQC model exists (decision 010); print the nQMaker citation for NQC | the report block prints 6 decimal places (633, in the fixed mode set at 1528-1529); the citation test matches only the substring "NONREV" | guarded by model type and name |
| `main/phylotesting.cpp`, `mixRevNonrev` 448-489 | S6 | mixed-family guard and fixed learned candidates | Level 3 | deferred |

Not touched in any slice: `model/modelmarkov.{h,cpp}` (decision 015),
`utils/optimization.{h,cpp}` (decision 029), `ModelProtein::init`, the likelihood kernels, the `PartitionModel`
constructor's frequency-pooling block (NQC's
`FREQ_USER_DEFINED` already skips it, verified `model/partitionmodel.cpp:118-119`), and every
vendored directory.

## 2. Specifications

### 2.1 Pure module `fixedeqchart`

All functions take and return plain values or Eigen objects plus a status, and none keeps static
mutable state.

- Target: strict parse of a frequency string with commas, spaces or slashes as separators;
  exactly 20 entries; every entry finite, greater than 0 and at least a caller-supplied minimum
  (`min_state_freq`); a sum within 1e-6 of 1, normalized once, with the pre-normalization sum
  returned (decision 006). It rejects NaN, infinities and the non-numeric tokens that IQ-TREE's
  own reader turns into random draws (verified, `model/modelmarkov.cpp:1761-1796`,
  `utils/tools.cpp:408-428`).
- Log-ratio jump build (D01): per row, a stable softmax over the 19 off-diagonal destinations
  with the reference logit fixed at 0; ν from the augmented system [1ᵀ; (K − I)ᵀ] ν = e₁ by
  column-pivoted Householder QR, the approach of `computeStateFreqFromQMatrix` (verified,
  `model/modelmarkov.cpp:2119-2131`); q_ij = ν_i K_ij / π*_i; the diagonal completed from the
  off-diagonal row sums and compared with −ν_i / π*_i.
- Inverse: K_ij = q_ij / (−q_ii) and z_ij = log(K_ij / K_{i,r(i)}).
- References: row maxima of a seed matrix, ties to the lower amino-acid index (D01, decision
  011).
- Seeds: a reversible seed becomes q_ij = R_ij π*_j normalized to unit mean rate; a
  non-reversible seed keeps its jump chain, which is the flux-preserving transfer of synthesis
  section 6.3.
- Diagnostics: ‖π*Q‖∞, scaled row residuals, |mean rate − 1|, minimum rate and flux, and the
  residual of the ν solve. They are reported and never used to reject a proposal (decision 026).
- Finite-difference step: signed h = η max(s, |x|), exactly representable, negative when x + h
  would leave the upper bound (decision 005).
- No other chart: no positive-ratio or T3 build is compiled (decision 028).

### 2.2 NQC model class

Construction (decisions 002, 006, 011):

1. Construct the `ModelProtein` base with the seed's name, LG by default or the `--init-model`
   name or file, because the base constructor calls `init(name)`, which cannot parse "NQC"
   (verified, `model/modelprotein.cpp:1079-1088`, `1107-1249`). Pass it an empty frequency
   string, so that `init` does not hand the target to `ModelMarkov::readStateFreq` (verified,
   `model/modelprotein.cpp:1241-1244`; decision 026).
2. Reject unsupported combinations (decision 008).
3. Parse and validate π* from `freq_params` with the pure module, under decision 006's rules.
4. For a reversible seed, set `state_freq = π*`, then call `setReversible(false)`, which
   converts R into R_ij π*_j and roots the tree (verified, `model/modelmarkov.cpp:112-119`,
   `145-149`). For a non-reversible seed, keep its jump chain.
5. Build LG at π* from the built-in models block and take the references from it, whatever the
   start (decision 011).
6. Encode the start through the inverse chart. The start is not checked against the domain; as
   for every IQ-TREE model, the line search clamps trial points into it (decision 027).
7. Set `freq_type = FREQ_USER_DEFINED`, `num_params = 360` and the name, then decompose.

Optimizer interface: `getNDim` returns 0 when fixed and 360 otherwise; `getNDimFreq` returns 0;
`setVariables` copies the coordinates out; `getVariables` compares exactly, copies them in, and
returns whether anything changed; `setBounds` sets the declared domain with `bound_check` false,
so there are no random restarts at dimension 360 (synthesis section 5).

Decomposition (decision 015): the `decomposeRateMatrixNonrev` override builds Q from the
coordinates into per-object workspace, validates it, writes the 380 off-diagonal rates into
`rates[]`, then calls the unchanged `ModelMarkov::decomposeRateMatrixNonrev`. With
`FREQ_USER_DEFINED` and `-optfromgiven` rejected, the base unpacks `rates[]` into `rate_matrix`,
skips the π reset and solve, rescales by a factor equal to 1 up to rounding, and eigendecomposes.
After the call the override checks that `state_freq` is bitwise π* and stops with an error
otherwise. The class never writes `state_freq`. A failed build writes nothing and sets a flag
that makes the class's `targetFunk` return 1e30. Every decomposition, from any caller, goes
through this override, so every eigensystem is built from the coordinates.

Derivatives (decision 005): a `derivativeFunk` override runs the scaled, bounds-aware forward
difference on the model itself; the same loop runs on the `PartitionModel` object through the
hook interface; a central stencil exists only for validation.

Post-fit checks (decision 009): after each fit, in an `optimizeParameters` override and in the
linked post-fit hook, re-evaluate the returned point, check the invariants, record bound
activity, and report whether the score fell below the start, restoring nothing (decision 029).
At export, compute a central-difference
first-order diagnostic in log-ratio coordinates; the stop is labelled "unclassified" unless that
diagnostic passes.

Guards: `setStateFrequency` and `adaptStateFrequency` accept only a vector bitwise equal to π*;
`setRateMatrix` and `setFullRateMatrix` re-encode through the inverse chart after checking
stationarity at π*, or reject.

Checkpoint (decision 003): chart identifier and version, target, references, coordinates as
17-digit strings, domain, derivative policy, fixed or fitted status, and accepted score, under
the class's own keys. Restore verifies the target, references, chart and domain against the
constructed object, rejects any mismatch, then rebuilds, decomposes and clears partial
likelihoods. If the class's restore chains to `ModelProtein::restoreCheckpoint`, that routine
restores `rates[]` at 10 significant digits and decomposes (verified,
`model/modelprotein.cpp:1265-1276`); the decomposition goes through the override and builds from
the coordinates held at that moment, so the class decomposes again after restoring its own.

Linked-compatibility signature: target, references, chart, domain and derivative policy,
compared across every object with the linked name before optimization.

Report and export (decision 010): `writeInfo` prints the target at 17 digits, the references,
domain, residuals, bound activity and stop classification, and the `.iqtree` Q block stays as it
is. The export writes `<prefix>.NQC.qmat` (20 rows of Q and the π* row at 17 significant digits,
the format `-m FILE` reads) and `<prefix>.NQC.info` (the metadata above, plus source commit and
seed identity).

## 3. Test specification

### 3.1 Layer 1: Python oracle

The package is `test_scripts/fixedeq/oracle/`, and its tests use the pass lines of decision 021.

- Charts (log-ratio, positive-ratio, T3 in the D04 gauge, and conversion from P-positive's star
  gauge), seeds and residuals. The positive-ratio and T3 charts are test references only
  (decision 028).
- Reproductions of the documented numbers: dimensions 360, 189 and 171; round trips; reversible
  reduction; the T3 triangle bound 20.72; the non-concavity curvature 0.01538 (P-log Appendix B;
  `verify_constrained_nq.py`).
- The three design scripts rerun from copies outside the repository, with their outputs compared
  against the values the design documents report, by `design_rerun.py` under decision 019. The
  24 values the 2026-10-04 rerun did not reproduce are not test targets (decision 020).
- Rooted likelihood: Newick, PHYLIP and FASTA readers; ambiguity as IQ-TREE treats it (B as N or
  D, Z as Q or E, J as I or L, verified `model/modelprotein.cpp:1350-1366`; other unknown
  characters as all states); discrete Gamma with IQ-TREE's mean categories (`gamma_median`
  defaults to false, verified `utils/tools.cpp:7276`); +I; sums in log space with scaling.
- Derivative rules: IQ-TREE's legacy step and the scaled step, with a central reference checked
  over several step sizes.
- IQ-TREE's optimizer, ported from `utils/optimization.cpp` with its failed-line-search path, for
  the S2 toy comparison (decision 017).
- A port of the target estimator (PLAN.md, "Where π* comes from").
- A direct balanced-flux reference optimizer for small problems (S5).

### 3.2 Layer 2: C++ unit tests

Fixtures at n = 3, 4 and 20: random balanced fluxes built from directed cycles rather than
through the chart under test, reversible matrices, skewed and near-reducible jump chains, and
adversarial conditioning. Each fixture's provenance header records the condition number κ of its
augmented ν system, computed by the oracle. Checks, with decision 026's pass lines: Q entries
and round trips within a relative max(1e-12, 1e-14 κ) per entry;
residuals at or below 1e-10; reference ties; target validation (zero, NaN, wrong count, below the
minimum, a sum more than 1e-6 from 1); and the step function at bounds.

### 3.3 Layer 3: differential likelihood

Against the unmodified binary first: `-m <Q file> -te <rooted tree> --show-lh`, or `-blfix
-keep-ident` if `--show-lh` gives too few digits, for NQ.pfam and for a random non-reversible Q
written at 17 digits; with and without a fixed +G4 and +I; on `example/aa_example.phy` and on
one partition of `turtle_aa`. `-blfix` also switches off the +I+G restart path (verified,
`utils/tools.cpp:3389-3395`). From S2: the Q the class builds at its seed equals the oracle's,
and IQ-TREE's first log-likelihood at the seed equals the oracle's. The output compared is the
initial log-likelihood that `--show-lh` prints at precision 17 (S0 probe (d)), and the tolerance
is 1e-8 relative (decision 021). The comparison against the unmodified binary is
`test_scripts/fixedeq/differential/differential.py`.

### 3.4 Layer 4: regression baseline

The binary is `~/iqtree3-baseline/iqtree3`, copied from the unmodified build with its SHA-256
recorded. Every run starts in its own new directory under `~/iqtree3-runs/baseline/`, uses
`-seed 1` and `-T 1` except runs 11 and 12, and is repeated twice; runs 11 and 12 are repeated
five times to measure their spread:

1. `-s example/aa_example.phy -m LG+G4`
2. `-s example/aa_example.phy -m NONREV`
3. `-s example/aa_example.phy -m "NONREV+F{0.08,0.06,0.04,0.05,0.02,0.04,0.07,0.07,0.02,0.05,0.10,0.06,0.02,0.04,0.05,0.07,0.05,0.01,0.03,0.07}"`
4. `-s example/aa_example.phy -m NQ.pfam`
5. `-s example/aa_example.phy -m GTR20`
6. `-s example/example.phy -m UNREST`
7. `-s example/example.phy -m 12.12`, a valid Lie-Markov name (verified,
   `model/modelliemarkov.cpp:217-225`)
8. `-s test_scripts/test_data/turtle_aa.fasta -p test_scripts/test_data/turtle_aa.nex --model-joint NONREV`
9. run 8 with `-S` in place of `-p` (`-S` takes a partition file, verified
   `utils/tools.cpp:2107-2116`)
10. `-s test_scripts/test_data/turtle_aa.fasta -p test_scripts/test_data/turtle_aa.nex -m run08.Q.txt -te run08.treefile`,
    with the same repeat's run 8 Q block and tree file copied under these relative names,
    because the report prints the Q file's path
11. run 8 with `-T 4`
12. run 10 with `-T 4` (decision 018)

Compared per run: every section of the `.iqtree` report except the opening header, "ALISIM
COMMAND" and "TIME STAMP", which hold paths, the build date and clock times; the tree file; and
the exit code. These contain the log-likelihood, the number of free parameters, the tree length
and the Q block. Wall time is recorded, not compared. An item identical in every baseline repeat
must reproduce exactly; one whose repeats differ only in their numbers is held to the observed
spread; one whose text differs otherwise cannot be compared. Run 11 is judged by decision 018
instead; the driver applies that rule when comparing, so the stored baseline keeps run 11's items
as recorded. The driver is `test_scripts/fixedeq/regression/regress.py`. The baseline of runs 1 to
11 and its provenance are in `test_scripts/fixedeq/regression/baseline/`; run 12's, recorded
later from the same frozen binary with run 8 rerun in each repeat to supply its inputs, is in
`baseline/run12/`. A comparison reads both files (`compare --baseline` given once for each).

### 3.5 Layer 5: invariants and nesting

- Nesting: fit `GTR20+F{π*}`, then seed NQC at that fit with its trees and rate parameters.
  IQ-TREE roots the GTR20 tree through its own conversion, and the reversible likelihood does not
  depend on the root, so NQC's first log-likelihood must equal the GTR20 fit's and its final one
  must not be lower. The fit, its transfer into the NQC run, the reference log-likelihood and
  the tolerance are those of decision 025: a fixed-tree fit carried from its checkpoint, compared
  with a `--show-lh` evaluation of the carried values within 1e-8 × |lnL|.
- Target immutability: π* compared bitwise after construction, linking, each linked update,
  checkpoint, restart and export.
- Counting: 360 plus rate and branch parameters, with the matrix counted once for linked
  partitions.
- Determinism: `-T 1` against `-T 4`, and permuted partition order. How the thread comparison
  is defined is open (PLAN.md risk 19), because legacy joint fits already differ across thread
  counts.
- Export and re-import: log-likelihood equality with fixed nuisances.
- An incompatible restart, with a different target, is rejected.

### 3.6 Layer 6: sanitizers

Build in `~/iqtree3-build-asan` with the `Mem` build type and the C, C++ and linker flags of
decision 022, `-j 2`, and run with that decision's `ASAN_OPTIONS` and `UBSAN_OPTIONS`, through
`test_scripts/fixedeq/sanitizer/sanitize.py`. The fallbacks are those of decision 022: without
cmaple, then UndefinedBehaviorSanitizer alone.

The build runs about 20 times slower than Release (measured 2026-10-05). In S0 it ran on the
unmodified tree over:
- baseline runs 2 and 4;
- the 16 differential cases;
- probe (a)'s five runs.

That covers a single NONREV fit, fixed-matrix evaluation with +G and +I on fixed trees, linked
fits under `-p` and `-S` on fixed trees, and the `--mdef` target route. At the end of S1 to S4 it
runs the unit tests and short NQC runs at `-T 1`.

### 3.7 Continuous integration

Upstream's workflow gives portability checks only (PLAN.md, test strategy). The fork workflow
runs layers 1 and 2 from S1, and a short NQC run from S2. Bitwise regression comparison stays
local.

## 4. Per-slice lists

S0, with no IQ-TREE source:

- Files: `test_scripts/fixedeq/` with `environment.yml`, `oracle/`, `tests/`, `regression/`,
  `differential/`, `probes/`, `sanitizer/`, `manifest/`, `design_rerun.py` and `.gitattributes`.
- Tests first: oracle tests that reproduce the documented numbers, and the regression driver.
- Runtime checks: (a) whether `NONREV+F<name>` with an `--mdef` file defining
  `frequency <name> = ...;` (syntax verified, `model/modelmixture.cpp:724`) reaches the model
  under `-m` and `--model-joint`; (b) whether a number such as `1e+00` inside `+F{...}` breaks
  parsing; (c) whether `-S -te <trees> -blfix` with fixed rate parameters freezes every
  nuisance; (d) which output carries 10 or more significant digits of the log-likelihood; (e)
  whether `-S` with the turtle NEXUS partition file behaves as `-p` does apart from tree linkage;
  (f) whether the root moves in a Level 1 and a Level 2 fit under `-te`, run once from an unrooted
  tree and once from a rooted one, comparing the root edge and the two branch lengths beside the
  root in the input and output trees (decision 016); (g) how precisely a fitted `GTR20+F{π*}`
  matrix on `example/aa_example.phy` carries into a new run through each route that needs no
  source change (the `.iqtree` block at 6 decimal places, and the 10-digit checkpoint values written as a
  `-m FILE`), with the tree and rate parameters fixed, recording the log-likelihood difference.
- Also: the unmodified sanitizer build, and the G0 run manifest (state order, reference rule,
  chart, target provenance, domain, derivative policy, root policy (decision 016), tree and rate
  policy, starts, source commit).

S1:

- Files: `model/fixedeqchart.{h,cpp}`, `unittest/`, `test_scripts/fixedeq/make_fixtures.py`,
  fixtures, `.github/workflows/fixedeq.yaml`.
- Tests first: layer 2.

S2:

- Files: `model/modelnonrevfixedeq.{h,cpp}`, `model/CMakeLists.txt`, `model/modelmixture.cpp`,
  `utils/tools.{h,cpp}` (domain option).
- Tests first: layer 3 at the seed; IQ-TREE's scaled gradient against the oracle's central
  reference at several steps; attempts to overwrite the target; invalid targets; the reported
  score equal to the re-evaluated score at the returned coordinates; the compiled optimizer
  against the oracle's port (decision 017) on one small identical problem; rejection of unsupported combinations;
  the full layer 4 comparison.

S3:

- Files: `model/partitionmodel.{h,cpp}`, the class checkpoint and linked hooks, and
  `utils/tools.{h,cpp}` (step option, if needed).
- Tests first: two partitions with intentionally different compositions under `-S` and `-p`
  (one Q and one target in every object, fixed outside the shared update); `-T 1` against
  `-T 4`; restart; nesting; the linked gradient's step; layer 4 again.

S4:

- Files: `main/phyloanalysis.cpp` and the class export.
- Tests first: export then `-m FILE` log-likelihood equality with fixed nuisances; no π-mismatch
  warning on re-import; metadata completeness; the citation text; layer 4 again.

S5 to S7: lists written when S4 is done.

## 5. Source anchors

Verified at `63c330d9` on 2026-09-23; HEAD `4c5f061f` matches it outside documentation. Anchors
marked as added 2026-10-03 were read at HEAD `4c5f061f` that day.

| Anchor | Fact |
|---|---|
| `utils/optimization.cpp:23` | `ERROR_X = 1.0e-4` |
| `utils/optimization.cpp:149-156` | `fixBound` clamps trial points |
| `utils/optimization.cpp:645-718` | `lnsrch`; on a step below `alamin` it restores `xold`, sets `check = 1` and leaves `*f` at the last trial (687-690) |
| `utils/optimization.cpp:750-778` | `minimizeMultiDimen`; restarts only for `bound_check` |
| `utils/optimization.cpp:793-900` | `dfpmin`; `check` is never read; small-step return at 843-846 |
| `utils/optimization.cpp:916-939`; `utils/optimization.h:145` | legacy forward difference; declared virtual |
| `tree/phylotreemixlen.h:165` | the only other `derivativeFunk` override |
| `model/modelmarkov.h:30-32, 345` | `MIN_RATE`, `TOL_RATE`, `MAX_RATE`; virtual `decomposeRateMatrixNonrev` |
| `model/modelmarkov.cpp:112-119, 145-149` | `setReversible` converts with the current `state_freq`, and roots the tree |
| `model/modelmarkov.cpp:204-206` | `freqTypeString` returns "" for protein `FREQ_USER_DEFINED` |
| `model/modelmarkov.cpp:914-930` | `adaptStateFrequency` frequency swap |
| `model/modelmarkov.cpp:964-1000` | `getNDim`, `getNDimFreq` |
| `model/modelmarkov.cpp:1018-1118` | `setVariables`, `getVariables`, `targetFunk` with the `min_state_freq` guard at 1103-1108 |
| `model/modelmarkov.cpp:1139-1164, 1199` | `setBounds` with `bound_check` false; the single-model optimizer call |
| `model/modelmarkov.cpp:1246-1387` | `decomposeRateMatrixNonrev`: reset and solve guards at 1252 and 1266, normalization 1269-1282, `--eigen` branch 1284-1289, Eigen path 1292-1386 |
| `model/modelmarkov.cpp:1761-1796, 1798-1879` | `readStateFreq(string)`; `readParameters`, with sign detection at 1809-1812 and the π check at 1834-1843 |
| `model/modelmarkov.cpp:1954-1968` | `getModelByName` and `validModelName`, with no `ModelsBlock` |
| `model/modelmarkov.cpp:2119-2131` | `computeStateFreqFromQMatrix` |
| `model/modelprotein.cpp:1079-1088, 1107-1249` | the constructor calls `init`; `init`, with the NONREV branch at 1207-1233 |
| `model/modelprotein.cpp:1333-1348, 1350-1366` | `getNameParams`; two-state tip likelihoods for B, Z and J |
| `model/modelprotein.cpp:1256-1276`; `model/modelfactory.cpp:1504, 1550` | `ModelProtein` checkpoints `rates[]` and decomposes on restore; the +I+G optimization restores the model on each restart and at the end (added 2026-10-03) |
| `tree/phylotree.cpp:5906`; `tree/iqtree.cpp:3257-3259`; `model/modelfactory.cpp:1733-1734`; `utils/tools.cpp:4699-4706, 7117-7118` | `convertToRooted` (midpoint or `-o` outgroup); root search inside NNI; `--root-find` after model optimization; `-te` sets no search iterations; `root_move_dist` 2 and `root_find` false by default (added 2026-10-03) |
| `model/modelsubst.cpp:199-209`; `tree/phylotreesse.cpp:366-373` | base tip likelihoods; tip partials come from `computeTipLikelihood` |
| `model/modelmixture.cpp:724, 3122-3267` | `frequency NAME = ...;` syntax; `createModel` |
| `model/modelfactory.cpp:259, 284-305, 489-509, 599-605` | `+F` tokenization; `model_joint`; `+F{...}`; `+F<name>` lookup |
| `model/modelfactory.cpp:1253-1257, 1331-1377, 1395-1505, 1931-1938` | `getNParameters`; the `-jointopt` path with fixed bounds; the +I+G checkpoint save and restore; the factory's `targetFunk` |
| `model/partitionmodel.cpp:35-169` | constructor: linking by name 78-82, DIVMAT `ASSERT` 96, pooling 116-168 with the skip at 118-119 |
| `model/partitionmodel.cpp:182-256, 299-338, 733-868, 883-958` | checkpoints and counting; `targetFunk`; linked optimization; outer loop with the `ASSERT` at 938 |
| `model/partitionmodelplen.cpp:75-194` | edge-proportional loop; partition rates 150-157; `ASSERT` at 135 |
| `main/phyloanalysis.cpp:164-179, 581-663, 1388-1402, 3903-3907` | citations; `reportModel`, precision set at 633; linked report; `.best_model.nex` |
| `main/phylotesting.cpp:165, 448-489, 1095, 1158-1159, 1380-1382` | non-reversible names; `mixRevNonrev`; comma split; mixing guard; `model_joint` skips ModelFinder |
| `utils/tools.cpp:588, 2107-2116, 2896, 2918-2933, 3156-3164` | comma split; `-S`; `--init-model`; `-mset` and `-madd`; `--min-freq` |
| `utils/tools.cpp:3315-3318, 3381-3382, 3389-3395, 3428-3437, 5023-5030` | `-optfromgiven`; `-jointopt`; `-blfix`; `--show-lh`; `--model-joint` |
| `utils/tools.cpp:7082, 7263, 7276, 7284` | `opt_gammai` default true; `min_state_freq` default; `gamma_median` default false; `optimize_from_given_params` default false (7284 added 2026-10-03) |
| `utils/checkpoint.h:22, 328, 352` | checkpoint precision of 10 significant digits |
| `alignment/alignment.cpp:4886-4933, 5882-5909` | ambiguity fixed point; `convfreq` |
| `utils/eigendecomposition.h:24, 81` | `ZERO_FREQ = 1e-10`; `total_num_subst` |
| `cmaple/CMakeLists.txt:281-293`; `CMakeLists.txt:255-258` | googletest fetch; CMAPLE not an option on Windows |
| `.github/workflows/ci.yaml` | upstream workflow, run on every branch push |
