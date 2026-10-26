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
| `model/modelnonrevfixedeq.{h,cpp}` | S2 to S4 | Class `ModelNonrevFixedEq : public ModelProtein`, public name "NQC" held in one constant (decision 002; section 2.2); from S3b also `public LinkedModelHooks`, listed second |
| `model/linkedmodelhooks.h` | S3b | Header only: the abstract `LinkedModelHooks` through which `PartitionModel::optimizeLinkedModel` lets a linked model check, differentiate and record its fits (decision 005; section 2.2, "Linked hooks"). A separate header because NQC's header is included by `model/modelmixture.cpp` and `partitionmodel.h` pulls in the super tree, so putting the interface in either would create an include cycle |
| `unittest/` | S1 | Standalone CMake project that compiles `model/fixedeqchart.cpp` directly, with IQ-TREE's clang Release flags and the `model` target's SSE flags (verified, `CMakeLists.txt:420`, `711-718`, `984`), and fetches googletest 1.12.1 at cmaple's pinned commit, checked by the archive's SHA-256 (verified, `cmaple/CMakeLists.txt:283-293`; decision 012); `fixtures/`, 27 chart cases and `targets.txt` in plain text at 17 significant digits with a provenance header (section 3.2); its own `.gitattributes` (decision 014) |
| `test_scripts/fixedeq/` | S0 onward | `environment.yml` and lock (decision 013); `oracle/`; `tests/`; drivers `regression/`, `differential/`, `probes/`, `sanitizer/`, `unittests/`, `nqc/` and `design_rerun.py`; `manifest/` (G0 manifest); `explore/` (exploratory measurements behind decisions 007, 021, 023, 025 and 038); `make_fixtures.py`; `.gitattributes` and `.gitignore`. Added in S3: `nqc/twopart.py`, the two-partition data (S3b); `explore/convergence.py`, the stopping study and risk 22's converged evidence, with `oracle/optimize.py`'s reference climber `lbfgsb_fit` (S3c's first part); `profile/profiling.py`, the profiling driver (S3c) |
| `.github/workflows/fixedeq.yaml` | S1 | Fork-only workflow that builds and runs the unit tests and the oracle's tests on Linux, and from S2 a short NQC run |

Run outputs never go into the working copy; they go under `~/iqtree3-runs/` in WSL. Fixtures and
the summarized regression baseline are derived data kept in the repository with provenance
headers, each produced by a kept script, and are never regenerated to make a test pass.

### 1.2 Existing files touched

| File and anchor (verified) | Slice | Change (proposed) | Why the file must change | How it stays opt-in or behaviour-preserving |
|---|---|---|---|---|
| `model/CMakeLists.txt` | S2 and S3b, done | add `modelnonrevfixedeq` and `fixedeqchart`; S3b adds `linkedmodelhooks.h` | the source list is explicit, with no glob | additive |
| `model/modelmixture.cpp`, `createModel`, 3230-3257 | S2, done | at the top of the dispatch chain, before the PoMo branch, dispatch the NQC name to `new ModelNonrevFixedEq(..., models_block)`, stopping with a message for non-protein data (decision 002) | the dispatcher that holds a `ModelsBlock`; non-protein data never reach the protein branch | reached only by the NQC name |
| `model/modelmixture.cpp`, `initMixture`, 3441-3445 and 3518 | S2, done | stop with an error when a component just created is an NQC model (decision 008) | `createModel` has no mixture flag, and `initMixture` writes component fields directly | fires only for NQC |
| `model/modelprotein.{h,cpp}`, beside the constructor at 1079-1088 | S2, done | a protected constructor that builds `ModelMarkov(tree, true, false)`, stores the models block and reads no model (decision 002) | the public constructor's `init` unroots a rooted tree when it reads a named model | called only by NQC; `ModelProtein::init` is unchanged |
| `utils/tools.{h,cpp}`, `parseArg` and `usage_iqtree` | S2 (domain, done); S4 (help line) | a parse-only option, read at the point of use through `Params` (decision 007); provisional name `--nqc-domain LO,HI`, in ratio units, default `1e-5,1e5` (`1e-5,10` until S3c). No step option: decision 007 adds one only for a step study, and S3 has none. S4: one line in `usage_iqtree`'s model list for `NQC+F{...}` and `+F<name>` | domain expansion needs a run-time setting; the model needs a help entry for the release | unused unless NQC is used; help text only |
| `model/partitionmodel.{h,cpp}`, `derivativeFunk` (declared at the orphan comment 160-165) and `optimizeLinkedModel` 753-840 | S3b, done | `derivativeFunk` delegates to the linked model's `LinkedModelHooks` when `dynamic_cast` finds them, else calls `Optimization::derivativeFunk`. In `optimizeLinkedModel`, when the hooks are present: build the list of partition models with the linked name in `part_order`; call `startLinkedFit` after the `ndim == 0` return and before `minimizeMultiDimen` (786); keep the optimizer's score before 824 overwrites it; call `endLinkedFit` after the `if (changed)` block (818-825), so before the save at 862 and the round's `ASSERT` (938; `partitionmodelplen.cpp:135`). Done in 27 added lines, with comments, and one in the header (decisions 005, 009) | linked optimization runs `minimizeMultiDimen` on the `PartitionModel` object, which no model override reaches, and only this function knows when a linked fit starts and ends and what the optimizer returned | a legacy model has no hooks, so the cast fails and every call takes the unchanged path; `minimizeMultiDimen` runs on `PartitionModel` only in this function (verified), and the cast only ever sees a `linked_models` entry, and RTTI is on (no `-fno-rtti`; the class already uses `dynamic_cast`) |
| `main/phyloanalysis.cpp`, `reportPhyloAnalysis` 1409-1416; citation block 164-179 | S3a (call); S4 (export, citation) | S3a: one call to `ModelNonrevFixedEq::reportRun` at the top of `reportPhyloAnalysis`, after its master-only return and before its `OUT_IQTREE` early return (Q1 of the S3 plan). S4: the export files and diagnostic at the same call (decisions 009, 010), only its comment changed; in `reportReferences` one flag so that `--model-joint` with an NQC model takes the nQMaker block rather than QMaker's (Peter, 2026-10-26: `--model-joint` only, NONREV's rule) | IQ-TREE prints no model details on the `+I+G` path at default verbosity, so S3-1's summary needs a print site outside the fit; `reportPhyloAnalysis` runs once per analysis on the master, after `--runs` has restored the best run, whereas `runTreeReconstruction` (3544) runs once per bootstrap replicate; the report prints Q at 6 decimals (633, in the fixed mode set at 1528-1529); the citation test matches only the substring "NONREV" | the call does nothing when the run has no NQC model |
| `main/phylotesting.cpp`, `mixRevNonrev` 448-489 | S6 | mixed-family guard and fixed learned candidates | Level 3 | deferred |

Not touched in any slice: `model/modelmarkov.{h,cpp}` (decision 015),
`utils/optimization.{h,cpp}` (decision 009), `ModelProtein::init`, the likelihood kernels, the `PartitionModel`
constructor's frequency-pooling block (NQC's
`FREQ_USER_DEFINED` already skips it, verified `model/partitionmodel.cpp:118-119`), and every
vendored directory.

## 2. Specifications

### 2.1 Pure module `fixedeqchart`

Free functions in namespace `fixedeq` over Eigen types of dynamic size, with no IQ-TREE headers
and no process exit; none keeps static mutable state. A function that can fail returns a
`Status` (`ok`, `invalidInput`, `nonFiniteCoordinate`, `softmaxUnderflow`, `nonPositiveNu`,
`nonFiniteRate`, `nonPositiveRate`, `invalidTarget`) and writes its outputs only on `ok`. Every
function works for any n ≥ 3, so the unit tests run at n = 3 and 4 as well as 20.

- Coordinates: row-major; within row i, the destinations other than i and r(i) in increasing
  index (`freeDestinations`; the G0 manifest's `chart.order`). This is the checkpoint layout,
  from S2.
- Target (`parseTarget`): strict parse of a frequency string with commas, ASCII whitespace or
  slashes as separators, a run of them counting as one; exactly the requested number of entries,
  each an ASCII decimal (`[+-]?(d+(.d*)?|.d+)([eE][+-]?d+)?`) read in the classic locale, finite,
  greater than 0 and at least a caller-supplied minimum; a sum, taken left to right, within 1e-6
  of 1, then each entry divided by that sum once, with the sum returned (decision 006). It
  rejects NaN, infinities, hexadecimal and the non-numeric tokens that IQ-TREE's own reader turns
  into random draws (verified, `model/modelmarkov.cpp:1761-1796`, `utils/tools.cpp:408-428`). It
  is stricter than the oracle's `parse_target`, which through Python also accepts non-ASCII
  digits, underscores and Unicode separators.
- Log-ratio jump build (D01; `jumpMatrix`, `solveStationary`, `generatorFromJump`, `build`): per
  row, a stable softmax over the n − 1 off-diagonal destinations with the reference logit fixed
  at 0, rejecting a jump probability that is exactly 0 after the exponential or after the
  division by the row sum (PLAN.md constraint 2; the oracle checks only the first); ν from the
  augmented system [1ᵀ; (K − I)ᵀ] ν = e₁ by column-pivoted Householder QR, the approach of
  `computeStateFreqFromQMatrix` (verified, `model/modelmarkov.cpp:2119-2131`) with the transpose
  formed explicitly, and ν required finite and positive; q_ij = (ν_i / π*_i) K_ij; the diagonal
  completed from the off-diagonal row sums; every off-diagonal rate required finite and strictly
  positive.
- Inverse: K_ij = q_ij / (−q_ii) and z_ij = log(K_ij / K_{i,r(i)}), with no box applied.
- References: row maxima of a seed matrix, ties to the lower amino-acid index (D01, decision
  011).
- Seeds: a reversible seed becomes q_ij = R_ij π*_j normalized to unit mean rate
  (`reversibleSeed`); a non-reversible seed keeps its jump chain, which is the flux-preserving
  transfer of synthesis section 6.3 (`transferSeed`).
- Diagnostics: ‖π*Q‖∞, scaled row residuals, |mean rate − 1|, minimum rate and flux, the
  residual of the ν solve, and the mismatch of the exit rates with ν_i / π*_i. They are reported
  and never used to reject a proposal (decision 021).
- Finite-difference step (`scaledStep`): signed h₀ = η max(s, |x|), negative when x + h₀ would
  leave the upper bound (decision 005), returned as the probe's displacement (x + h₀) − x,
  IQ-TREE's own idiom (verified, `utils/optimization.cpp:930-931`). When |x| ≥ h₀ the
  subtraction is exact (Sterbenz), so the step is exactly the distance the probe moved; it can
  still differ from h₀ in the last bits, because x + h₀ rounds.
- First-order test (`firstOrderTest`, S4): decision 009's export label from coordinates x, a
  gradient g of −lnL, the box, the objective f and gtol, as section 2.2 item 4 states it. It
  returns the statistic, its coordinate (the lowest index on a tie), the counts at each bound and
  outside, and whether it passed; non-finite input, unequal or empty sizes, or gtol ≤ 0 are
  `invalidInput`. The oracle's `chart.first_order_test` repeats it operation for operation.
- No other chart: no positive-ratio or T3 build is compiled (D01).

### 2.2 NQC model class

Construction (decisions 002, 006, 011):

1. Construct the `ModelProtein` base through its protected constructor, which reads no model
   (the public constructor's `init` unroots a rooted tree when it reads a named model; verified,
   `model/modelprotein.cpp:1121`, `model/modelmarkov.cpp:1856`, `97-101`). Then read the seed as
   the NONREV branch of `init` does (verified, `model/modelprotein.cpp:1214-1231`): the
   `--init-model` name from the models block, or the `--init-model` file, otherwise LG, each with
   tree adaptation off, and without `rescaleRates` (decision 002). The target never reaches
   `ModelMarkov::readStateFreq` (decision 006).
2. Reject unsupported combinations (decision 008).
3. Parse and validate π* from `freq_params` with the pure module, under decision 006's rules,
   with the minimum max(`min_state_freq`, `ZERO_FREQ`), so that D07's check against both holds.
   Report both the raw and the normalized vector.
4. For a reversible seed, set `state_freq = π*`, then call `setReversible(false)`, which
   converts R into R_ij π*_j, roots an unrooted tree and keeps a rooted one (verified,
   `model/modelmarkov.cpp:112-119`, `145-149`). For a non-reversible seed, keep its jump chain
   (`transferSeed`) and do not call `setReversible(false)`, which zeroes `rates[]` on a model that
   is already non-reversible (verified, `model/modelmarkov.cpp:133-137`); apply its rooting rule
   instead (lines 145-149), since nothing else roots the tree for a model read with adaptation off.
5. Build LG at π* from the built-in models block and take the references from it, whatever the
   start (decision 011).
6. Encode the start through the inverse chart. The start is not checked against the domain; as
   for every IQ-TREE model, the line search clamps trial points into it (decision 007).
7. Set `freq_type = FREQ_USER_DEFINED`, `num_params = 360` and the name, then decompose.

Optimizer interface: `getNDim` returns 0 when fixed and 360 otherwise; `getNDimFreq` returns 0;
`setVariables` copies the accepted coordinates out; `getVariables` compares exactly with the
trial coordinates, copies them in, and returns whether anything changed; `setBounds` sets the
declared domain with `bound_check` false, so there are no random restarts at dimension 360
(synthesis section 5).

Two coordinate arrays. The *trial* array holds the coordinates `getVariables` last received; the
*accepted* array holds the last coordinates that built successfully, and is what `setVariables`,
the checkpoint and the report use, so a proposal that failed to build is never saved or reported
(PLAN.md constraint 2). A failed build leaves the accepted array, `rates[]` and the eigensystem as
they were and sets a flag tied to the trial array. The next `getVariables` call with different
coordinates rebuilds, and the flag clears once a build succeeds; the array is never reverted on
failure, so a return to the last good point always compares as a change.

Decomposition (decision 015, which states the override's steps): the base then unpacks `rates[]`
into `rate_matrix`, skips the π reset and solve, rescales by a factor equal to 1 up to rounding,
and eigendecomposes. The class never writes `state_freq` after construction.

`targetFunk` repeats the base routine's statements (verified, `model/modelmarkov.cpp:1091-1118`)
and returns 1e30 when the flag is set. `optimizeParameters` repeats the base routine's statements
in order (verified, `model/modelmarkov.cpp:1166-1244`), keeping the optimizer's own returned
score, which the base routine discards, for the post-fit checks.

Derivatives (decision 005): one loop, `scaledDerivative(Optimization &objective, x, dfx)`,
runs the scaled, bounds-aware forward difference on any objective (S3b; in S2 it is
`derivativeFunk`'s body). The single-model `derivativeFunk` passes the model itself, and
`linkedDerivative` passes the `PartitionModel` object. The loop keeps the start-score capture and
the `-vv` first-gradient print. A central stencil exists only for validation.

Post-fit checks (decision 009, which states what they do) run in `checkFit` after every
single-model fit, in `endLinkedFit` after every linked fit (S3b), and at export (S4). From S3a:
- `state_freq` is compared bitwise with the parsed target after every fit, stopping with
  `outError` on a mismatch. This is the decomposition guard's rule, applied where no
  decomposition may follow: between NQC fits IQ-TREE updates rates and branches without
  decomposing Q.
- Under `-v` the record is printed once per fit, as three lines in S2's wording (scores; bounds;
  residuals), so `nqc/fitrecords.py` still parses it. It replaces S2's `writeInfo(cout)` call,
  which IQ-TREE's `-v` round loop then repeated.
- A re-scored lnL below the start is printed after those lines with `outWarning`: in the log at
  every verbosity, on screen except under `-quiet`.
- `checkFit` runs only outside parallel regions, because unlinked partitions are refused; a
  comment records that.

Guards: `getStateFrequency` returns π* unchanged, without the base's renormalization (decision
015); `setStateFrequency` and `adaptStateFrequency` accept only a vector bitwise equal to π*;
`setRateMatrix` and `setFullRateMatrix` stop with an error rather than re-encode, until a
reachable caller exists: their only callers (rate-linking in mixtures, site-specific
frequencies) are refused.

Refusals (decision 008; this is the full list):
- non-protein data, in `createModel`;
- `+F`, `+FQ`, a bare `+FU`, no `+F`, and `NQC{...}`;
- `-optfromgiven`; `-jointopt`; `--link-exchange-rates`;
- `+T` tree mixtures and the HMM options (`-hmmster` requires `+T`, `utils/tools.cpp:5544-5549`);
- GHOST `+H`. From S3b it is read from the partition's own model string (`aln->model_name`)
  as well as `params.model_name`, because a charset written `NQC+F{...}+H4` never reaches the
  latter;
- site-specific frequencies (`-fs`, `-ft`), checked second, after ModelFinder, because the
  site-specific path builds each model without the `+F` string (`model/modelfactory.cpp:655`);
  ModelFinder is checked first, since a candidate's string may also lack it;
- `--eigen` (other matrix-exponential methods behave as for NONREV);
- ModelFinder in any form. From S3a it is detected when `Params::dating_mf` is set. S2 matched
  TEST or MF at the start of `-m`, which missed an omitted `-m` with `-mset NQC+F<name>` and a
  partition file with a partition lacking a model. `runModelFinder` sets the flag before it builds
  any candidate (`main/phylotesting.cpp:1389`, ahead of the fast tree at 1433 and of the candidate
  tests; MixtureFinder also reaches it), and no NQC model is built before ModelFinder runs. So it
  is the run-time signal, without editing ModelFinder, although named for dating; a comment says
  so. Not covered: IQ-TREE's library build, which copies `Params` (`main/main.cpp:3627`); NQC is
  not offered there;
- AliSim (simulate from the exported matrix);
- unlinked partitioned use: a partition file (`-p`, `-q`, `-S`) without `--model-joint` or
  `--link-model` (from S3a). The check comes after the AliSim check, which also sets
  `partition_file` (`main/main.cpp:3233`), so AliSim keeps its own message;
- `-mstore` (from S4): from the end of each model optimization `ModelFactory` serves transition
  matrices from a cache keyed by branch length, which the export's probes would reuse, and which
  linked fits may reuse after each partition's own update (not verified; maintainer question 9);
- in `initMixture`, any mixture component.

IQ-TREE's parser rewrites `+FO{...}` and `+FU{...}` to `+F{...}` for every model (verified,
`model/modelfactory.cpp:484-485`), so those spellings carry a target as `+F{...}` does.

Checkpoint (decision 003; implemented in S2): chart identifier and version, target,
references, coordinates as 17-digit strings, domain, derivative policy, fixed or fitted status,
and accepted score, under the class's own keys. The accepted score is the re-scored
log-likelihood after the model's last accepted fit, under the rate parameters and branch lengths
of that moment, so in Level 2 it is not the run's final score: IQ-TREE's later updates in the
round raise the final score above it. S4's re-import equality therefore runs at Level 1, where
nothing moves after the last fit, and no final-score key is added (Peter, 2026-10-22). For a
linked model the scores are those of the linked objective, the sum over its NQC partitions.

S2's keys also include `fits`, the number of fits made. From S3a the keys add `start_score` and
`optimizer_score` (17 digits, beside `accepted_score`), `accepted_fit` (the fit that produced the
accepted state), `bounds` (lower, upper, outside), `residuals` (the seven diagnostics of section
2.1 at 17 digits) and `evaluations` (calls of `targetFunk`, for S3c's profiling), and `target` is
written from `state_freq` rather than the class's copy. Restore:
- coordinates, the three scores and `accepted_fit` come back; a score missing from an older
  checkpoint becomes NaN, printed as unknown;
- bounds and residuals are recomputed from the rebuilt Q, because they are derived. In S2 the
  `+I+G` restore of the best start brought back its coordinates but left the last start's bounds,
  residuals and start and optimizer scores (found in the S3 survey);
- `fits` and `evaluations` are counters, so restore keeps the larger of the current and saved
  values. The `+I+G` restore (`model/modelfactory.cpp:1503-1505`) then no longer winds them back
  past the starts made after the best one, and a resumed run, starting from 0, continues from the
  saved counts.

Restore verifies the target, references, chart and domain against the constructed object,
rejects any mismatch, then rebuilds, decomposes and clears partial likelihoods. A restore that finds none of the class's
keys leaves the model unchanged, as IQ-TREE's own restores do: the +I+G path copies the checkpoint before any model
state is saved (verified, `tree/iqtree.cpp:2365`, `2389`). If the class's restore chains to
`ModelProtein::restoreCheckpoint`, that routine restores `rates[]` at 10 significant digits and
decomposes (verified, `model/modelprotein.cpp:1265-1276`); the decomposition goes through the
override and builds from the coordinates held at that moment, so the class decomposes again
after restoring its own.

Verbose output (decision 021): under `-v` the class prints the seed Q after the base routine's
rescaling at 17 significant digits, and each fit's record (post-fit checks, above); under `-vv`
it prints the first gradient of its first fit with the steps taken, which in a linked run is the
representative's first linked fit. None appears in a normal run. In a linked fit the
`VB_DEBUG` "proposal rejected" print runs inside `PartitionModel::targetFunk`'s parallel loop,
so it is under `#pragma omp critical`, the idiom at `model/partitionmodel.cpp:911-915`.

Linked-compatibility signature (S3b). Every object with the linked name must be an NQC object
and must equal the representative in:
- target and `state_freq`, bitwise;
- references, domain (bitwise), chart and version, and derivative policy;
- accepted coordinates, bitwise (state consistency, synthesis section 7.3).

Otherwise `startLinkedFit` stops with `outError` naming what differs. It checks before every
linked fit, so the first check comes after round 1's per-partition nuisance updates but before
the shared matrix moves; "before optimization" in D06 and decision 011 means before the shared
matrix is optimized. Refusing before those updates would need a second hook call in the
`PartitionModel` constructor. They compute nothing wrong, since each partition uses its own
target, so that extra shared edit was not made.

Linked hooks (S3b, decision 005): `LinkedModelHooks` in `model/linkedmodelhooks.h`.
- `startLinkedFit(const vector<ModelSubst*> &linked)` checks the signature above and sets the
  start flag.
- `linkedDerivative(Optimization &objective, double x[], double dfx[])` runs `scaledDerivative`
  on `objective`.
- `endLinkedFit(const vector<ModelSubst*> &linked, double optimizer_lnl)`:
  - re-scores by summing each linked object's partition lnL directly, in the list's `part_order`
    order, which is bitwise the sum `PartitionModel::targetFunk` forms, without counting an
    evaluation. Re-scoring through `targetFunk` would record K × 1e30 at a point whose build
    failed. A failed build is warned about, and the last valid state is scored, as `checkFit`
    does;
  - records the fit in every linked object, so every partition's checkpoint struct agrees, each
    object computing its own bounds and residuals and checking its own `state_freq`;
  - has the representative print the `-v` record and any fall, once, before the round's `ASSERT`.
- The end call was kept rather than having the class report only at its report and export,
  which was the question left open for S3's plan. Without it the class cannot learn that a
  linked fit ended, or the optimizer's score. Decision 009's post-fit check could then not run
  for linked fits, and S3-1's linked half could not be met.

Thread safety (constraint 7) of the linked path:
- per object, `targetFunk` touches only its own trial and accepted arrays, Q, ν, `rates[]` and
  eigen arrays;
- `fixedeq::build` is pure;
- the base decomposition already runs in this loop for NONREV;
- the hooks run outside the parallel region, and the class has no static mutable state.

Report and export (decision 010). From S3a the `writeInfo` override is removed, so NQC prints
nothing from IQ-TREE's own call sites, as NONREV does (`ModelMarkov::writeInfo` prints nothing for
a 20-state non-reversible model, `model/modelmarkov.cpp:368-376`). Those sites are the `-v` round
loop (S2's duplicate), the default-verbosity print at the end of a fit
(`model/modelfactory.cpp:1747-1749`), the linked representative's default print
(`model/partitionmodel.cpp:953-956`, `model/partitionmodelplen.cpp:186-189`) and the kernels'
numerical-failure dumps.

In their place, the static `ModelNonrevFixedEq::reportRun(IQTree&, ostream&)` prints one block
per NQC model at default verbosity:
- target as given with its sum, target normalized from `state_freq` at 17 digits, references,
  domain;
- fits, the accepted fit's three scores, bounds, residuals, likelihood evaluations;
- for a linked model (S3b), one block per linked model, from its representative and with the
  number of partitions, not one per partition object.

It is called once per analysis from the top of `reportPhyloAnalysis` (section 1.2). With `-b 1`
IQ-TREE reports the single replicate and then the original (callers 4678 and 4738, section 5),
so a summary follows each of the two reports. A run that stops at the linked `ASSERT`
never reaches it; for that case the fall line printed before the `ASSERT` and the checkpoint keys
are the record (accepted).

S4 adds the export and the stop classification to `reportRun`, in this order (decisions 009,
010):
1. **`<prefix>.NQC.qmat`**, written first so that an abort inside a probe cannot lose the fit:
   the 20 rows of `rate_matrix`, the normalized Q the likelihood used, then `state_freq`, bitwise
   π*; each value `%.17g`, separated by single spaces, LF line ends, and no header, since the
   reader chooses a full matrix from the first token's sign. IQ-TREE's reader re-solves π from Q
   and applies `rescaleRates`, so a re-import matches the fit to rounding, not bitwise.
2. **One line at default verbosity** announcing the diagnostic and its evaluations. It also
   bounds the parsers that count switch and rejection lines: those after it come from probes.
3. **The diagnostic**, on the objective the fits used: the model for one alignment; for a linked
   model the `PartitionModel`, whose public `model` member is set to the representative for the
   diagnostic and restored after, the swap `optimizeLinkedModels` makes (outside it the member is
   a dummy "JC", and `targetFunk` matches partitions by its name and unfixes each per call, keeping
   partitions parallel as in fits).
   - Before: copy x0, the accepted coordinates (sized by `accepted.size()`, since `getNDim` is 0
     for a fixed model), and every linked object's `evaluations`; check that every linked object
     holds x0 bitwise; rebuild each at x0, then `initializeAllPartialLh` and `clearAllPartialLH` on
     the tree, because `-b N`'s report follows `copyTree` and `readTreeString` with no
     reinitialization (`optimizeConTree` reinitializes after the same re-read).
   - f0 = targetFunk(x0). For each k, probes at x0 ± h with h = 1e-3 max(1, |x_k|), 021's
     central-reference step; g_k = (f⁺ − f⁻) / (h⁺ + h⁻), h± the probes' displacements. A probe
     that returns 1e30 is a failed probe.
   - After: targetFunk(x0), then the counters are put back, so `evaluations` counts fits only.
     An `ASSERT` checks that this value equals f0 and that `rate_matrix` equals its earlier copy,
     both bitwise: the `.iqtree` prints that matrix and `-z` reuses the model.
   - 722 evaluations: 6 to 14 percent of a first linked fit's.
4. **The label**, from `fixedeq::firstOrderTest` with gtol 1e-4, the gtol NQC's fits receive
   (`gradient_epsilon` defaults to 1e-4, and `TOL_RATE` is 1e-4). Per coordinate the usable part
   of g, the derivative of −lnL, is |g_k| inside the box, max(0, −g_k) at the lower bound and
   max(0, g_k) at the upper; s_k = that × max(|x_k|, 1) / max(|f0|, 1), dfpmin's scaling. "first-order
   stationary at the optimizer's tolerance" when max s_k < 1e-4, no coordinate lies outside and no
   probe failed; otherwise "unclassified", with the reason. It certifies no peak (039).
5. **The summary** gains the label and both file paths, and says that the fit is re-used with
   `-m <prefix>.NQC.qmat`.
6. **`<prefix>.NQC.info`**, one `key: value` per line.
   - Written with the checkpoint's own strings: `target`, `chart`, `chart_version`,
     `references`, `domain`, `derivative`, `coordinates`, `fixed`, `fits`, `evaluations`,
     `accepted_fit`, `start_score`, `optimizer_score`, `accepted_score` ("unknown" where the
     checkpoint has no key), `bounds`, `residuals`.
   - Also: `state_order`; `target_as_given` and `target_sum`; `reference_letters`; `start` (LG
     rebuilt at the target, or the `--init-model` name or file and whether it was rebuilt or
     transferred, recorded at construction); `random_seed`; `command` (the checkpoint's `iqtree`
     struct); `linked_partitions` (0 for one alignment); `export_lnl` (f0, the log-likelihood of
     the exported state under the run's final nuisances, not a checkpoint key); `stop`,
     `stop_reason`, `stop_statistic`, `stop_coordinate` (1-based, as the `-vv` gradient print),
     `stop_line`, `stop_step`, `stop_failed_probes`, `stop_gradient` (360 values);
     `diagnostic_evaluations`; `qmat`; `iqtree_version`.
   - No source commit: the binary does not know it, and the drivers record commit, clean-tree
     status and the binary's SHA-256 beside every run (Peter, 2026-10-26).

The `.iqtree` Q block stays as it is. So do the `.iqtree` "Model of substitution" line, a
partitioned run's `.best_model.nex` and the `.iqtree` AliSim command, none of which carries the
target or the fitted matrix, so none re-runs the fit; the summary points to the export instead
(Peter, 2026-10-26; maintainer question 8).

## 3. Test specification

### 3.1 Layer 1: Python oracle

The package is `test_scripts/fixedeq/oracle/`, and its tests use the pass lines of decision 021.

- Charts (log-ratio, positive-ratio, T3 in the D04 gauge, and conversion from P-positive's star
  gauge), seeds, residuals, and κ of the augmented ν system (`chart.nu_condition`). The positive-ratio and T3 charts are test references only
  (D01).
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
- A linked objective (S3b, `nqc.LinkedObjective`): the partitions' `Objective`s at the same
  coordinates, summed from 0 in partition order, each with its sites, the run's input tree text
  (never the run's own `.treefile` or a checkpoint newick) and its rate model, and a partition
  rate (`-p`, `-q`) folded into the rate categories; `normalize_part_rates` repeats
  `PhyloSuperTreePlen::normalizePartRate`'s operations in order, for `-q`'s `{x}` rates.
- A port of the target estimator (PLAN.md, "Where π* comes from").
- A direct balanced-flux reference optimizer for small problems (S5).

### 3.2 Layer 2: C++ unit tests

`unittest/test_fixedeqchart.cpp`, run by `test_scripts/fixedeq/unittests/unittests.py`
(section 3.6).

The fixtures are written by `test_scripts/fixedeq/make_fixtures.py` from the oracle and committed
from a clean tree (in S1, commit `c349e41a`, seed 20261019).
- **27 chart cases** at n = 3, 4 and 20.
  - 17 have ground truth built without the chart's ν solve: random balanced fluxes from directed
    cycles (9), reversible seeds of random exchangeabilities (6), two exact 3-state cases (one
    with a non-uniform ν, which a transposed solve gets wrong, and a reference tie), and LG
    rebuilt at the decision 023 target.
  - 10 are built by the oracle's chart: NQ.pfam's jump chain transferred to that target, skewed
    targets with 19 entries at 1e-4, near-reducible chains at n = 4 and 20, coordinates across the
    first default box and beyond it, and one logit near the exponential's underflow.
- **Headers.** Each records the condition number κ of the case's augmented ν system
  (`chart.nu_condition`) and the oracle's own error against the pass line.
- **Self-check.** The generator writes nothing if any case fails decision 021's lines in the
  oracle itself, and it refuses a non-empty output directory.
- **`targets.txt`.** 25 target strings with the oracle's verdicts, and each accepted vector
  normalized as the compiled parser does it.

Checks use decision 021's pass lines; everything else is compared exactly.
- **Q.** Q built from each fixture's z, and Q built from the inverse of the fixture's Q, agree
  with the fixture within a relative max(1e-12, 1e-14 κ) per entry, the diagonal included. For a
  fixed reference set the build is injective, so the two checks pin the inverse without a line on
  z.
- **Residuals.** ‖π*Q‖∞, scaled row sums and |mean rate − 1| are at most 1e-10 on every fixture,
  hard cases included.
- **Seeds and references.**
  - Reversible and transfer seeds match their fixtures under the Q line.
  - References equal the fixtures': the row maxima, or the row maxima of LG rebuilt at the case's
    target.
  - Also checked: the oracle's tie case, and the G0 manifest's LG references.
- **Rejections**, with the outputs left untouched: a NaN coordinate; softmax underflow after the
  exponential and after the division; a rate that underflows to 0; non-positive and non-finite
  rates in the inverse; invalid targets and sizes.
- **Target parser.** Every oracle verdict is reproduced, accepted vectors and sums match bitwise,
  and non-ASCII input that Python accepts is refused.
- **Step.**
  - The displacement identity holds bitwise.
  - The displacement is exact, checked by TwoSum, whenever |x| ≥ h₀.
  - The step is negative at and beyond the upper bound, and forward below the lower bound
    (decision 007).
- **Diagnostics.** On a dyadic Q, before and after one perturbation, compared exactly.
- **First-order test (S4).** Hand-built cases compared exactly: a pass and a failure on either
  side of gtol inside the box; at each bound, a gradient pointing out of the box ignored and one
  pointing in counted; a coordinate outside failing with a zero gradient; a tie; the rejections.

**Residual table.** Each fixture test records its errors and residuals with googletest's
`RecordProperty`, as strings, because googletest 1.12.1's only numeric overload takes an `int`.
The driver turns them into `residuals.tsv`.

### 3.3 Layer 3: differential likelihood

Against the unmodified binary first: `-m <Q file> -te <rooted tree> --show-lh`, or `-blfix
-keep-ident` if `--show-lh` gives too few digits, for NQ.pfam and for a random non-reversible Q
written at 17 digits; with and without a fixed +G4 and +I; on `example/aa_example.phy` and on
one partition of `turtle_aa`. `-blfix` also switches off the +I+G restart path (verified,
`utils/tools.cpp:3389-3395`). From S2: the Q the class builds at its seed equals the oracle's,
and IQ-TREE's first log-likelihood at the seed equals the oracle's. The output compared is the
initial log-likelihood that `--show-lh` prints at precision 17 (S0 probe (d)), and the tolerance
is 1e-8 relative (decision 021). The comparison against the unmodified binary is
`test_scripts/fixedeq/differential/differential.py`. The S2 seed runs add `-nparam 2`, which skips
the fitting loop (verified, `utils/tools.cpp:2860-2867`, `model/modelfactory.cpp:1636`), since
`--show-lh` alone does not stop a model with free parameters from being fitted; the class's Q is
read from its `-v` output (decision 021). The S2 checks, with decision 021's further lines, are
`test_scripts/fixedeq/nqc/nqcchecks.py`, at these points (formerly decision 030):
- **Gradient.** `example/aa_example.phy`, decision 023's target, `+G4{0.5}`, the rooted tree of
  baseline run 2 with branch lengths fixed; at the LG seed, and at a start built by the oracle
  from the seed's coordinates with one free coordinate per row set to 1e-9, written at 17
  significant digits and passed by `--init-model`. The legacy step there is reported, not
  asserted.
- **Toy problem (decision 017).** 8 taxa, 3000 sites, no rate variation, simulated from NQ.pfam's
  jump chain transferred to the target on a rooted tree with branch lengths uniform on
  [0.05, 0.3], by the frozen binary's AliSim (or a seeded simulator in the oracle if AliSim
  refuses the Q file), every seed recorded. IQ-TREE makes exactly one
  Level 1 fit (`-te -blfix -nparam 3`); the oracle's port makes one `optimize_model` call with
  the scaled step, the same bounds and gradient tolerance; both start at LG rebuilt at the
  target. Asserted: the start scores, the first gradients and the final scores (decision 021).

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
  with a `--show-lh` evaluation of the carried values within 1e-8 × |lnL|. In S3a the GTR20 fit
  is `GTR20+F{π*}+G4` on aa_example with `-te` the unrooted topology of baseline run 2, and the
  evaluation uses that fit's tree file. The NQC run is a Level 1 fit from `--init-model` the PAML
  file, on the same tree file and α, with `-blfix -vv`. Its first log-likelihood is −lnL at the
  first-gradient print, the first fit's start; the `start_score` key holds the last fit's start,
  and a Level 1 run makes several fits. Its final one is `accepted_score`. The NQC fit runs with
  `--nqc-domain LO,HI` computed by the driver from the carried start before the run: one decade
  beyond the start's smallest and largest jump-weight ratios, rounded outward to powers of ten,
  and never inside the default box (LO at most 1e-5, HI at least 1e5); for `aa_example`'s start,
  ratios of 1.19e-7 and 1.13e5 (S3a's gate) give 1e-8 and 1e7. This is a test setting; decision 007's
  default, with no start margin, is unchanged. The same fit in the default box is reported, not
  asserted.
- Target immutability: π* compared bitwise after construction, linking, each linked update,
  checkpoint, restart and export.
  - A single-alignment run dumps no checkpoint between fits, so "after every fit" is the
    class's post-fit guard, observed as exit 0, plus the final key from `state_freq`.
  - A linked fit run uses `-cptime 0 --all-checkpoint` and reads every numbered checkpoint.
    Otherwise a linked checkpoint is dumped at most once a minute (`utils/checkpoint.cpp:150`),
    so per-fit checks would see nothing.
- Counting: 360 plus IQ-TREE's own count of rate, branch and partition-rate parameters, with the
  matrix counted once for linked partitions. IQ-TREE already counts a linked matrix once, by
  unfixing it while it counts (`model/partitionmodel.cpp:250-254`), so no code changes. The
  formulas, read in `ModelFactory`, `PartitionModel` and `PartitionModelPlen::getNParameters`
  and set here before the tests:
  - single alignment and `-S`: per partition, rate dimensions plus branches (none under
    `-blfix`, else the rooted tree's branch count less one);
  - `-p` and `-q`: the super tree's `branchNum`, which counts the root-leaf edge and ignores
    `-blfix`; plus P − 1 partition rates under `-p` and none under `-q`.
  - Rate dimensions: `+G4` 1, `+I+G4` 2, free `+R4` 6 (2 × 4 − 2), braced 0.
  - Expected:

    | Run | Count |
    |---|---|
    | aa_example, Level 1 / Level 2 `+G4` / `+I+G4` / `+R4` | 360 / 395 / 396 / 400 |
    | two-partition data (12 taxa), `-S` Level 1 | 360 |
    | two-partition data, `-S` Level 2 `+G4` / `+I+G4` / `+R4` | 406 / 408 / 416 |
    | two-partition data, `-p` with `-blfix` | 384 |
    | two-partition data, `-p` Level 2 `+G4` / `+I+G4` / `+R4` | 386 / 388 / 396 |
    | two-partition data, `-q` Level 1 | 383 |

  - The count is read from the report, as in S2. S3c also reports NONREV's count minus NQC's, 19
    in every configuration.
- Linked fits (S3b), on `nqc/twopart.py`'s data: one seeded rooted 12-taxon tree; partition A
  simulated by the frozen binary's AliSim from NQ.pfam's jump chain at decision 023's target, B at
  that target cubed and renormalized; 600 sites each, with a two-charset NEXUS file; no missing
  taxon and no identical sequences in a partition (`-p` removes identical sequences), so `-p`'s
  partition trees equal the shared tree. It reuses `nqc/toy.py`'s tree and AliSim helpers, whose
  defaults keep the toy data byte-identical, since the toy case checks its SHA-256 against the
  calibration.
  - Seed and gradient: `-S` and `-p`, `+G4{0.5}` and the given `+R4`, `-te -blfix -vv`, one round
    (`--loop-model 1`; `-nparam 2`); nothing moves before the first linked gradient. Asserted:
    −lnL there against the oracle's linked objective (1e-8); printed steps equal
    fl(fl(x + d) − x) bitwise, d by decision 005's rule with `upper` from the `domain` key; the
    same-rule and central-reference lines (021); at least one coordinate where the legacy step
    differs (0 < |x| < 1; at the LG seed every coordinate is at most 0), so that the bitwise check
    tells the rules apart (had the legacy path run, there would be no print); free parameters.
    Also a `-S` start with the first free coordinate of each row at 1e-9 (`--init-model`),
    asserted to hold 20 coordinates within 1e-8 of zero.
  - Level 1: `-S -te -blfix` with braced rates; `-q` by decision 024, its tree, gamma shapes and
    `{x}` rates from the `-p` Level 2 `+G4` fit; `--link-model` under `-S`, a charpartition giving
    both charsets `NQC+F{π*}+G4{0.5}`. Level 2: `-S` and `-p`
    with `+G4`, `+I+G4` and free `+R4`. All with `-cptime 0 --all-checkpoint`; the Level 2 `+G4`
    runs use `-vv`, since the restart test reads them. Asserted at every numbered checkpoint: π*
    bitwise in the linked struct and every partition struct, coordinates bitwise equal across
    them, references, residual keys; at the end, free parameters and the bound keys; at Level 1,
    the accepted score against the oracle (1e-8). Under `-S`, `+I+G4` takes IQ-TREE's multi-start
    in every round, so its restore is asserted block by block. Reported: `--link-model` against
    `--model-joint`; falls, `ASSERT`, the "BUG" abort, scaling-squaring, bounds. No oracle score at
    Level 2, since no file holds the branch lengths in force at a fit.
  - Refused: `--link-model` with charsets at two targets, or at the same target written to 8
    decimals, at the first linked fit; a charset without a model and no `-m` (ModelFinder runs;
    with `--model-joint` it would not); a charset `NQC+F{...}+H4`.
- Determinism (decision 038), on turtle: the linked log-likelihood at the seed and the first
  gradient at `-T 1` against `-p` at `-T 4` and `-S` at `-T 3`, 10 runs of each, within decision
  021's lines, never bitwise. Under `-p`, `SuperAlignment(Params&)` caps the threads at the sum
  over partitions of max(1, patterns × states / 4000), which on turtle is 3, one per partition,
  and prints the change; other data can allow more threads per partition. Under `-S` with fewer
  partitions than threads, `PhyloSuperTree::setNumThreads` gives each partition all the threads
  and runs the partitions one at a time, so NQC's code never runs concurrently; `-S` at `-T 3`
  gives one thread per partition, as read in that function (confirmed in S3b by the CPU-to-wall
  ratio). Fits at more than one thread are reported beside NONREV's spread (decision 038): in
  S3b, 5 `-p` Level 2 fits on turtle at `-T 4` without rate variation, the setup of NONREV's
  measured spread (`explore/thread_counts.py`'s pL2), beside one at `-T 1`.
  - The runs use the turtle trees of baseline runs 8 (`-p`) and 9 (`-S`), `-te -blfix`,
    `+G4{0.5}`, `-vv`, and one round (`-nparam 2` under `-p`, `--loop-model 1` under `-S`). The
    seed lnL and gradient are those of the first-gradient print.
  - The driver's runner must pass the thread count, since it appends `-T 1` last and the last
    `-T` wins.
  - Nothing logs the allocation under `-S`, so it is reported as the CPU-to-wall ratio.
- Permuted partition order, the charsets reversed (trees reversed for `-S`):
  - on the two-partition data at Level 1, seed lnL, first gradient and final coordinates are
    bitwise equal, because a sum of two partitions is commutative;
  - on turtle at the seed, within the lines of decision 038's thread row, which bound the same
    rounding.
- Export and re-import (S4; pass lines in decision 021). Asserted from files; printed output
  (the citation, IQ-TREE's "does not match" line on re-import, the summary's lines) is reported
  (Peter, 2026-10-26).
  - Re-import at Level 1 against `accepted_score`: aa_example's `fit-level1` export with
    `-m probe.NQC.qmat+G4{0.5} -te input.tree --show-lh`, by this build and by the frozen v3.1.4
    binary; the two-partition `-S` export with `-S tp.nex -m probe.NQC.qmat+G4{0.5} -te
    partition.trees --show-lh`, its per-partition 17-decimal lines summed; the `-q` export through
    the oracle's `file_model` at the `-q` fit's nuisances, since `-q` and `-p` print the initial
    log-likelihood to 3 decimals and no file holds more.
  - Content: the `.qmat` holds 20 rows of 20 and a π row, a negative first entry, positive
    off-diagonals and the target key's values bitwise; its Q against the oracle's build from the
    checkpointed coordinates; the residuals of its Q, computed by the oracle; every `.info` field
    present with its value count, and the fields shared with the checkpoint identical to its keys.
  - Diagnostic at the one-alignment, `-S` and `-q` Level 1 exports: the gradient against the
    oracle's central differences; the statistic, coordinate and label recomputed by the oracle's
    `first_order_test` from the C++ gradient; the label against the oracle's own. The no-fit export
    at the LG start (the `root` case's run) is "unclassified", the oracle's label checked first.
    If no export is "first-order stationary", that is reported, and S5 revisits a positive
    control (Peter, 2026-10-26).
  - Restore: an `--undo` rerun of `fit-level1`'s copied checkpoint writes a byte-identical
    `.qmat` and the same `.info` apart from `command`; the S4 binary and the S3c binary give
    identical outputs on `fit-level1`, `linked-fit-s-level1` and a one-alignment Level 1 run with
    `-z`, the only likelihood computed after `reportRun`; a `-b 2` Level 1 run on the short data
    exits 0 with its export.
  - Reported: the citation; the "does not match" line, beside the oracle's max |π̂ − π*|; a
    `+F{π*}` re-import, whose root is π*/Σπ* through `getStateFrequency`; IQ-TREE's 3-decimal `-q`
    re-import; every export's label; export_lnl against the printed final score; the diagnostic's
    time.
- Restart (S3b):
  - a `-p` and a `-S` Level 2 run with `-cptime 0 --all-checkpoint`;
  - its checkpoint #1, written after the first linked fit (under `-p` before any branch update,
    under `-S` after round 1's per-partition updates), is copied with the inputs into a new
    directory as `<prefix>.ckp.gz`, and the same command is rerun; that #1 holds `fits` 1 is
    asserted first, since the test relies on it, and a failure there stops for Peter rather than
    picking another checkpoint;
  - asserted: the first linked gradient is at #1's coordinates bitwise; `fits` continues from
    #1's value (the first resumed checkpoint whose coordinates differ from #1's holds #1's count
    plus 1, where a reset counter would hold 1); π* is bitwise unchanged; exit 0.
  - Checkpoint #1 only, because the restarted run re-reads the `-te` tree before restoring.
    No likelihood comparison, because α is refitted before the first linked fit.
  - `--undo` on a finished `-te` run is checked by exit 0 only: IQ-TREE re-saves no model
    after that restore, so its keys cannot show it.
- An incompatible restart, with a different target or domain, is rejected (`--undo`, since a
  changed command line only warns, `main/main.cpp:2476-2494`).
- Single alignment (S3a), from a finished Level 1 run's copied checkpoint and inputs: `--undo`
  with another target and with another domain each refused by NQC's restore; `--undo` with the
  same command exits 0 (reported: its summary against the original's); and S2's gate run
  `~/iqtree3-runs/s2-gate/20261021-35d4a1bf/nqc/fit-level1`, whose checkpoint lacks S3a's keys,
  rerun with `--undo` (reported: exit 0, scores unknown, bounds and residuals recomputed). The
  driver's `restart` case needs that directory as `--older-run`.

### 3.6 Layer 6: sanitizers

`test_scripts/fixedeq/sanitizer/sanitize.py` builds and runs it (decision 022), in a new
directory for each slice (S0's `~/iqtree3-build-asan`, S2's `~/iqtree3-build-asan-s2`, then
`-s3a` and `-s3b`; S3c built none, since its only C++ change was a default value and a help string,
Peter, 2026-10-25).
- **Build.** `-DCMAKE_BUILD_TYPE=Mem` (`-g -O1`) with clang 14;
  `-fsanitize=address,undefined -fno-omit-frame-pointer` in both `CMAKE_C_FLAGS` and
  `CMAKE_CXX_FLAGS`, and `-fsanitize=address,undefined` in `CMAKE_EXE_LINKER_FLAGS`; `-j 2`.
  Fallbacks, in order: the same build with `-DUSE_CMAPLE=OFF` if the build or link runs out of
  memory, then UndefinedBehaviorSanitizer alone.
- **Run.** `ASAN_OPTIONS=detect_leaks=0:allow_user_segv_handler=0` and
  `UBSAN_OPTIONS=print_stacktrace=1`.
- **Coverage.** At the end of S2 to S4, the slice's own short NQC runs at `-T 1` (in S2: the seed
  case, the refusals, and Level 1 and `+I+G4` Level 2 fits on a cut-down toy alignment), since
  the build runs about 20 times slower than Release and a full `aa_example` fit would take hours.
  S0 ran, on the unmodified tree, baseline runs 2 and 4, the 16 differential cases and probe
  (a)'s five runs (a single NONREV fit, fixed-matrix evaluation with +G and +I, linked fits under
  `-p` and `-S`, the `--mdef` target route); S2 repeated that list on its binary. S3a and S3b
  repeat it and add their own short cases: S3a the new refusals and the short runs reading keys;
  S3b a cut-down two-partition Level 1 `-S` and Level 2 `-p` run (150 sites per partition; the
  `-S` run at default verbosity, since at `-v` IQ-TREE's closing NNI print overflows on an `-S`
  tree, PLAN.md maintainer question 6).
- **Unit tests.** They are a standalone project this build does not compile, so
  `unittests.py run --variant asan-ubsan` builds them with the same flags at `-g -O1` and runs
  them with the same options.

### 3.7 Continuous integration

Upstream's workflow gives portability checks only and never runs NQC. The fork workflow runs
layers 1 and 2 from S1, and a short NQC run from S2. Bitwise regression comparison stays local.

The fork workflow is `.github/workflows/fixedeq.yaml`:
- **Trigger.** It runs on pushes to `nq-constrained-pi`, guarded by
  `github.repository == 'petergoodman/iqtree3'`.
- **Unit tests.** Both variants on `ubuntu-22.04`, with the flags of `unittests.py`. A
  sanitizer report fails the job.
- **Oracle tests.** In the environment created by micromamba from `environment.lock.txt`.
- **NQC run, from S2.** Builds IQ-TREE on Linux and runs the S2 seed case.
- **Upstream.** The file is left out of the upstream pull request.

## 4. Per-slice lists

A closed slice keeps one line pointing to its gate record. S0 also keeps its runtime checks,
which define the probes that decisions and records cite.

S0 (done 2026-10-05, no IQ-TREE source; G0 manifest, `CHANGELOG.md` 2026-10-05):

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

S1 (done 2026-10-19, gate at `6c8d3af4`; `CHANGELOG.md` 2026-10-19): the pure module, its unit
tests and fixtures, and the fork workflow (sections 2.1, 3.2, 3.7).

S2 (done 2026-10-22, gate at `35d4a1bf`; gate record in `CHANGELOG.md`, S2 entry): the NQC class
with its checkpoint, its hooks in shared files and the domain option; the oracle's NQC objective
and toy problem; the check driver `test_scripts/fixedeq/nqc/nqcchecks.py`; the sanitizer step and
the CI job. Its known shortfalls are S3-1.

S3, approved 2026-10-22 as three slices, each with its own gate. The approved plan, with the
survey and review findings behind it, is frozen in `docs/agent/archive/S3_plan_2026-10-22.md`;
read it only where these lists are unclear.
- **Why three slices.** S3-1 changes the record and report that every linked test reads, so it
  comes first and the linked tests are written once.
- S3b alone edits `PartitionModel`, so its regression and sanitizer evidence isolate that change.
- S3c is days of laptop compute and may force a performance fix, which should neither hold back
  the correctness gates nor mix into them.
- The cost is three sanitizer builds rather than one, and two extra closes.
- **Implementation, one new session per slice,** started in plan mode with the prompt in
  PLAN.md's next step.
  - The lists below are complete, so the session's own plan is short (its commits and gate
    steps), and Peter can accept it and switch to auto.
  - Commits within a slice, as in S2: oracle additions with their tests; C++; driver cases; the
    close.
  - Each slice ends with its close commit and CHANGELOG entry, which is the hand-off to the next
    session. A slice that runs long stops after any commit, and the next session resumes from
    the CHANGELOG entry. The gate's long runs leave their outputs on disk and can be checked in
    a later session.
  - Before any run expected to last more than a few minutes (each gate's sanitizer step, S3c's
    profiling), the agent gives the compute estimate, starts the keep-awake helper (CLAUDE.md)
    and monitors the run.
  - Read-only subagents help twice: a code review before the gate (S2's plan review found the
    rooting defect), and the records-against-gate check at the close.
  - Not used: parallel agents or worktrees, because CLAUDE.md allows one working copy and the
    WSL build and run directories are shared; a multi-agent workflow, because each slice
    depends on the previous gate.
- **Checks that read output text,** as in S2.
  - Print behaviour, meaning whether and how often something is printed (the record once per
    fit, the summary, the fall line), is read and reported, not asserted. That follows Peter's
    rule against asserting on printed output.
  - Asserted values come from the class's checkpoint keys, except these, all carried over from
    S2:
    - the printed values decision 021 names, because no file holds them: the seed Q (`-v`) and
      the first gradient with its −lnL and steps (`-vv`);
    - `--show-lh`'s 17-digit first lnL;
    - a refusal, asserted as a non-zero exit with an "NQC" error line; a restart refusal also
      needs the restore's own error ("the checkpoint was written with a different"), the only
      sign that the restore refused it (S3a); from S3b the linked refusals name their reason too
      ("linked partitions must share one model", "ModelFinder", "heterotachy");
    - the free-parameter count, read from the report;
    - the `+I+G` restore check: the restored model scores as the best start, from the printed
      start scores at IQ-TREE's 3 decimals (S2's `restore_check`, missing from this list until
      S3a's close; Peter kept it asserted on 2026-10-22, as the only test of the restore in an
      ordinary fit). Under `-S` each partition takes the multi-start in every round, so the
      linked run is checked block by block (Peter, 2026-10-23).

S3a (done 2026-10-22, gate at `4749f222`; gate record in `CHANGELOG.md`, S3a entry): S3-1, the
S2 build's shortfalls (report once, default visibility, checkpoint record, ModelFinder and other
refusals, unlinked partitions, final score), closed for one alignment; the oracle's FreeRate
categories; the driver's `+R4`, nesting, skewed, restart and risk22 cases (sections 2.2, 3.5,
3.6). S3-1's linked half was closed in S3b.

S3b (done 2026-10-23, gate at `054914de`, the cases its fixes changed rerun at `fb5bf76d`; gate
record in `CHANGELOG.md`, S3b entry): linked training through `LinkedModelHooks` and two calls in
`PartitionModel`, with the refusal of `--model-joint` and `--link-model` lifted and each linked
fit checked and reported once; the oracle's linked objective, `nqc/twopart.py`, and the driver's
linked, thread, permuted, restart, risk 22 and warm-start cases (sections 1.2, 2.2, 3.1, 3.5,
3.6).

Known and left alone after S3b:
- a linked NQC run prints QMaker's citation, because the test matches "NONREV" (S4);
- `--link-alpha` is neither refused nor tested; nor is `-Q` with `--model-joint`, which takes
  `-S`'s `PartitionModel` path;
- after a tree search IQ-TREE checkpoints no final model (`main/phyloanalysis.cpp:3896-3899`), so
  `--undo` restores pre-final state, for every model; NQC's restart tests use `-te`;
- PLAN.md risk 23, reported by the S3 review and not verified, outside S3's paths;
- at `-v` and above, a Release `-S` run ends by printing meaningless NNI counts (PLAN.md
  maintainer question 6).

Monitoring, reported for every S3a and S3b fit and not asserted:
- whether a linked round reaches IQ-TREE's `ASSERT` on a fall of more than 0.1 (formerly risk 7;
  decision 009), and IQ-TREE's "BUG: Tree logl gets worse after model optimization!" abort
  (`tree/iqtree.cpp:2346-2351`), a second hard stop for single-model, `-p` and `-q` runs;
- the number of "INFO: Switch to scaling-squaring" lines under `-v` (formerly risk 12): the
  non-reversible transition matrix leaves the eigen path when P's row sums deviate by more than
  1e-4, and the decomposition flags a singular eigenvector matrix, both shared with NONREV
  (anchors in section 5). S2: none in the 183 fits of its three `aa_example` runs;
- falls below the start, and the optimizer's score against the re-scored one.

S3c (done 2026-10-26; first part 2026-10-24 to 2026-10-25, profiling gate at `0b0048af`; gate
record in `CHANGELOG.md`, S3c entry, and the entries of 2026-10-24 and 2026-10-25): the stopping
study and risk 22's converged evidence (`explore/convergence.py`, the oracle's `lbfgsb_fit`); the
default domain [1e-5, 1e5] (decision 007); profiling against NONREV by `profile/profiling.py`, with
gdb-timed windows of 1000 evaluations and counted full fits, under the scaled design Peter
approved on 2026-10-25 (its departures and checks are in `0b0048af`'s commit message and this
section at that commit). S-3 stays recorded as failed (case B).

Evidence for PLAN.md risk 22 (S3a's and S3b's single fits and warm starts, S3c's first part's
converged peaks; `CHANGELOG.md`): settled by decision 007's revision on 2026-10-25. The
`nqcchecks.py` cases `risk22`, `risk22-turtle` and `warm-starts` still report the default against
`1e-8,1e7`.

S4, planned and approved 2026-10-26, with ten choices by Peter recorded where they apply
(decisions 008, 009, 010, 021; sections 2.2, 3.5). The plan was checked by three surveys and an
adversarial review against the source.
- **Found while planning.** Decision 010's reason said the `-m FILE` reader keeps the file's π as
  root; it re-solves π from Q (corrected in place). Since `0e9fb159` the `restart` case's rerun of
  S2's checkpoint was refused by the domain check, so it now passes S2's `--nqc-domain 1e-5,10`.
  `upstream/master` had not moved from v3.1.4 on 2026-10-26.
- **Files.**
  - `model/fixedeqchart.{h,cpp}`: `firstOrderTest` (section 2.1), with unit tests.
  - `model/modelnonrevfixedeq.{h,cpp}`: the export, diagnostic and label in `reportRun` (section
    2.2), the start's description recorded at construction, and the `-mstore` refusal.
  - `main/phyloanalysis.cpp`: the citation flag; the call's comment.
  - `utils/tools.cpp`: the help line.
  - The oracle: `chart.first_order_test`, `nqc.export_gradient`, and `readers.read_qmat` and
    `read_info`, each with pytest.
  - `nqc/nqcchecks.py`: cases `export` and `linked-export`, the export checks in `short` and
    `linked-short` (sanitizer), the `-mstore` refusal, the `-b 2` run, and `--previous-binary`
    (`~/iqtree3-binaries/81729f4b/iqtree3`, the S3c binary).
- **Tests** (section 3.5): re-import; content; diagnostic; restore; refusal. Printed output is
  reported, which changes two items this list carried before planning (the π-mismatch warning and
  the citation text). Layers: regression 12 of 12 with no `.NQC.*` file in any regression run;
  unit tests in both variants; pytest; the sanitizer build `~/iqtree3-build-asan-s4`, run with the
  S0 list and the short cases, `-b 2` included.
- **Release candidate** (Peter): the close's gate record carries an audit of each PLAN.md software
  success criterion against its evidence, residual maxima included; NQC appears in `-h`.
- **Commits:** records and oracle; C++; driver cases; review fixes; the close.

S5 to S7: lists written when S4 is done; S7's starts from decision 039's convergence protocol.

## 5. Source anchors

Verified at `63c330d9` on 2026-09-23; HEAD `4c5f061f` matched it outside documentation. Anchors
marked as added 2026-10-03 were read at HEAD `4c5f061f` that day. Since S2 (`cd4c6823`) the
branch differs from `63c330d9` in `model/CMakeLists.txt`, `model/modelmixture.cpp`,
`model/modelprotein.{h,cpp}` and `utils/tools.{h,cpp}`; anchors in those files still refer to
`63c330d9`, and at HEAD the lines after S2's insertions are offset by up to 15. Since S3a
(`bb89e2d9`) `main/phyloanalysis.cpp` differs too, its anchors still at `63c330d9`: at HEAD, lines
after 60 are offset by 1 and lines after 1412 by 3. Since S3b (`7388dd52`) so do
`model/partitionmodel.{h,cpp}`, anchors likewise at `63c330d9`: at HEAD, in the `.cpp` lines after
25 are offset by 1, after 752 by 9, after 784 by 24, after 788 by 25 and after 825 by 27; in the
`.h` lines after 165 by 1.

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
| `CMakeLists.txt:420, 711-718, 984`; `utils/optimization.cpp:930-931` | clang Release flags `-O3 -ffunction-sections -fdata-sections`, no fast-math; SSE flags `-D__SSE3 -msse3`, set on the `model` target; the derivative step taken as the probe's displacement (added 2026-10-19) |
| `.github/workflows/ci.yaml` | upstream workflow, run on every branch push |
| `model/modelprotein.cpp:1108, 1121, 1208-1211, 1214-1231`; `model/modelmarkov.h:179`; `model/modelmarkov.cpp:97-101, 133-137, 1856` | `init` asserts 20 states; a named model is read with tree adaptation on; NONREV's warnings; NONREV reads its seed with adaptation off; `setReversible(true, true)` unroots a rooted tree; `setReversible(false)` on a non-reversible model zeroes `rates[]` (added 2026-10-20) |
| `model/modelmarkov.cpp:886-891`; `tree/phylokernelnonrev.h:647, 1123` | `getStateFrequency` renormalizes its copy; the non-reversible kernel takes root frequencies from it (added 2026-10-20) |
| `model/modelmarkov.cpp:465-479, 501-502` | scaling and squaring for non-reversible models; `ASSERT(0)` for other techniques (added 2026-10-20) |
| `model/modelfactory.cpp:484-485, 1441-1550, 1519, 1636-1704`; `utils/tools.cpp:2860-2867, 7409` | `+FO{...}`, `+FU{...}` rewritten to `+F{...}`; the +I+G checkpoint loop and its final `ASSERT`; the round loop and its stale printed score; `-nparam`; `modelEps` 0.01 (added 2026-10-20) |
| `tree/iqtree.cpp:2333, 2365, 2389` | +I+G path chosen by `opt_gammai`; the model restore before optimization and the first save after it (added 2026-10-20) |
| `model/modelmixture.cpp:3230-3257, 3441-3448, 3518-3525, 4374-4391`; `main/phyloanalysis.cpp:190, 3568, 3586` | `createModel`'s dispatch chain; `initMixture`'s component creation and direct writes; the linked-GTR `setRateMatrix`; the `+H` test; the tree read before the model (added 2026-10-20) |
| `utils/optimization.cpp:680-690, 832-855` | `lnsrch`'s failure path; `dfpmin` returns on zero displacement before its next gradient, and its scaled-gradient stop (added 2026-10-20) |
| `utils/tools.cpp:7262`; `model/partitionmodel.cpp:132-133, 155-160` | `keep_zero_freq` true by default, so `convfreq` returns at once; missing sequences counted as unknown states except under `-S`; "Mean state frequencies" printed at precision 8 (moved from PLAN.md 2026-10-22) |
| `model/modelmarkov.cpp:72, 76, 489-500, 1330-1336`; `main/phyloanalysis.cpp:3880, 3895-3899` | the constructor's `setReversible` call and the read of unset `is_reversible`; the switch to scaling and squaring on row sums off by more than 1e-4; `nondiagonalizable` on a singular eigenvector matrix; no final optimization on a fixed tree; the post-search save that writes the search state, not the model (moved from PLAN.md 2026-10-22) |
| `main/phylotesting.cpp:1340-1389, 1433`; `main/main.cpp:3233, 3627` | `runModelFinder`'s early returns, then `dating_mf` set at 1389, before the fast tree (1433) and any candidate; AliSim sets `partition_file`; the library build copies `Params` (added 2026-10-22, the S3 plan; these files and those in the rows below are unchanged since `63c330d9`, and the `utils/tools.cpp` lines are given at `63c330d9`) |
| `main/phyloanalysis.cpp:1409-1416, 3544, 3867-3870`; callers 4259, 4678, 4738, 4751, 5377, 5720 | `reportPhyloAnalysis` with its master-only and `OUT_IQTREE` returns; `runTreeReconstruction` and its MPI-worker return; `reportPhyloAnalysis` is called once per analysis (added 2026-10-22) |
| `model/partitionmodel.cpp:786, 818-825, 862-863, 888, 911-915, 953-956`; `model/partitionmodelplen.cpp:97, 100, 135, 186-189` | the linked `minimizeMultiDimen`; the re-decomposition after it; the save and dump after each linked fit; `-S`'s round count `model_opt_steps` (`--loop-model`); the critical-section print idiom; the default-verbosity linked `writeInfo`; under `-p`, the low-precision "Initial log-likelihood", the round count `num_param_iterations` (`-nparam`), the linked `ASSERT` and the linked `writeInfo` (added 2026-10-22) |
| `model/modelfactory.cpp:1379-1538, 1487-1489, 1503-1505, 1585-1597, 1678-1680, 1747-1749` | `optimizeParametersGammaInvar`, its save of the best start and its final restore; the initial likelihood timed as one evaluation with Q already decomposed; the `-v` per-round and default end-of-fit `writeInfo` (added 2026-10-22) |
| `model/modelmarkov.cpp:361-376`; `utils/tools.cpp:110-112` | `ModelMarkov::writeInfo` prints nothing for a 20-state non-reversible model; `outWarning` writes to `cout` (added 2026-10-22) |
| `model/ratefree.cpp:31-77, 286-296, 629-657` | FreeRate given values: rates divided by Σw·r, proportions renormalized when their sum is off by more than 1e-5, `fix_params = 2`; `getNDim`; the EM loop, which reaches the substitution model only through `setTree` and the base transition matrices (added 2026-10-22) |
| `tree/phylosupertreeplen.cpp:45-77, 101-119`; `tree/phylosupertree.cpp:229-233`; `alignment/superalignmentunlinked.cpp:25-34` | `-q`'s `{x}` rates and `normalizePartRate`; `setNumThreads`; no thread cap under `-S` (added 2026-10-22) |
| `utils/checkpoint.cpp:146-216`; `utils/tools.cpp:5050, 5060, 5068, 5544-5549`; `main/main.cpp:2268-2295, 2476-2494`; `tree/iqtree.cpp:2346-2351, 2363-2394` | the dump interval and `--all-checkpoint`'s numbered files; `--undo`, `-cptime`, `--all-checkpoint`, and HMM set only with `+T`; a finished run refused without `--undo`; a changed command line only warns; the "BUG: Tree logl gets worse" abort; `ensureModelParametersAreSet` restores the model, then skips or runs the fit (added 2026-10-22) |
| `tree/phylosupertreeunlinked.cpp:20-26`; `tree/mtree.cpp:687-701` | `-S` reads one `-te` tree per partition; `-p` and `-q` read only the first (added 2026-10-22) |
| `model/modelfactory.cpp:643-670` | the site-specific frequency path (`-fs`, `-ft`) builds each model with an empty frequency string (added 2026-10-22, S3a) |
| `utils/tools.cpp:2083-2100, 3389-3395`; `model/partitionmodel.cpp:901-904, 961-966` | `-p`, `-q` and `-blfix` set `opt_gammai` false, `-S` does not; `PartitionModel::optimizeParametersGammaInvar` sets `opt_gamma_invar`, so under `-S` each partition's `ModelFactory::optimizeParametersGammaInvar` runs its `+I+G` starts in every round; `ModelFactory`'s fallback for a plain `RateHeterogeneity` (`model/modelfactory.cpp:1380-1381`) applies only to the super tree's own call, which `PartitionModel` overrides (corrected 2026-10-23, S3b review; the S3a-close row said partitioned data never take the path) |
| `tree/iqtree.cpp:3248-3254`; `model/modelfactory.cpp:1395-1398` | in `IQTree::doNNISearch`, the save of each re-optimized model to the top level of the checkpoint; the copy of the checkpoint the `+I+G` starts restore from (PLAN.md risk 23, unverified; added 2026-10-22) |
| `model/modelprotein.cpp:1236-1237, 1247`; `model/modelmarkov.cpp:1252-1267, 1834-1843` | a user model file is read, then `rescaleRates`; an unknown frequency type becomes user-defined; the decomposition inside `readParameters` runs while the type is still unknown, so π is re-solved from Q, and the file's row only feeds the 1e-3 print (added 2026-10-26, S4 plan) |
| `utils/tools.cpp:3428-3438`; `model/modelfactory.cpp:1596-1599`; `model/partitionmodelplen.cpp:97`; `main/main.cpp:2469-2470` | `--show-lh` (fixed branch lengths, `VB_DEBUG`, no checkpoint, keeps identical sequences, does not stop the run); the "1. Initial log-likelihood" line, once per partition under `-S`; the 3-decimal total under `-p` and `-q`; `cout` in fixed mode, so 17 decimals (added 2026-10-26) |
| `main/phyloanalysis.cpp:86-267, 164, 172`; `main/phyloanalysis.cpp:3532, 4677-4678, 4730-4738`; `utils/tools.cpp:3507-3510` | `reportReferences`, written to the `.iqtree` only, and its QMaker and nQMaker conditions; `optimizeConTree` reinitializing after a re-read; the `-b 1` report of the single replicate; the `-b N` report after `copyTree` and `readTreeString` with no reinitialization; `-b 1` turns off the original-alignment tree (added 2026-10-26) |
| `model/partitionmodel.cpp:842-868`; `model/modelfactory.h:224`; `main/main.cpp:2473, 2494-2500` | `optimizeLinkedModels` sets the public `model` member to the linked model, unfixes, fits, refixes and restores it; the `iqtree` checkpoint struct and its `command` record (added 2026-10-26) |
| `utils/tools.cpp:3064, 7254`; `model/modelfactory.cpp:1582, 1756, 1800-1817` | `-mstore`, off by default; transition-matrix storing stopped at the start of each model optimization and started at its end; the cache keyed by branch length (added 2026-10-26) |
