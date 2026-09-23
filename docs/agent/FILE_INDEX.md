# IQ-TREE 3 — Context File Index for nQ-with-Target-π Work

**Audience:** AI coding agents working in this repository on *constrained optimization of
non-reversible amino-acid substitution matrices (nQ) under a user-supplied target stationary
frequency vector π*.

**What this file is.** A quick guide to which file to open for which question, with line
anchors into the code as it exists at the baseline commit, so that agents can go straight to the
right place instead of grepping the repository. It is one of three reference documents that
describe the current code and prescribe nothing: `docs/agent/ARCHITECTURE.md` covers the
programming architecture, `docs/agent/AA_MODEL_INFERENCE.md` the mathematics, optimization and
algorithms of amino-acid model inference, and this file the navigation. It is not a plan or a
decision record. "Load when" lines say when a file is worth opening, not what to build; the
project's design is in `docs/agent/design/`, its decisions in `docs/agent/DECISIONS.md`, and its
plan in `docs/agent/PLAN.md`. A passage here that reads as a recommendation for the project is a
defect; report it rather than follow it.

**How to use this file.** Do not `grep` the whole tree. Load files by tier. Tier 0 is mandatory
reading before writing any code. Tiers 1–4 are loaded on demand, keyed by the "Load when" line
of each entry. Tier 6 is an explicit *do-not-read* list — those directories look relevant by
name and are not.

**Baseline: `63c330d9` (upstream tag `v3.1.4`).** Every line number below was re-verified
against that commit on 2026-09-15. They drift; treat them as anchors to grep near, not as ground
truth. `git diff --stat 63c330d9 HEAD -- <path>` tells you whether a given file has moved since.

Companion document: **`docs/agent/ARCHITECTURE.md`** — read that first for how the system fits
together, then use this file to find things.

---

## Tier 0 — Mandatory. Read these before writing any code.

### `model/modelmarkov.h` (558 lines)

The base class `ModelMarkov : public ModelSubst, public EigenDecomposition` for every
reversible *and* non-reversible Markov substitution model.

Key declarations:

- `is_reversible` flag (protected) — the single switch that changes almost all behaviour.
- `double *rates` — the free rate parameters. Packing differs by reversibility (see the
  architecture doc).
- `double *rate_matrix` — full `num_states × num_states` Q, row-major, non-reversible only.
- `int num_params` — number of free rate parameters.
- The optimization contract: `getNDim()`, `setVariables()`, `getVariables()`, `setBounds()`,
  `targetFunk()`, `optimizeParameters()`, `setRates()`.
- `decomposeRateMatrix()`, `decomposeRateMatrixNonrev()`, `decomposeRateMatrixRev()`.
- `getStateFrequency() / setStateFrequency() / adaptStateFrequency()`.
- `MIN_RATE = 1e-4`, `TOL_RATE = 1e-4`, `MAX_RATE = 100` (lines 29–31) — the default box
  constraints applied to every rate parameter.

**Load when:** always.

---

### `model/modelmarkov.cpp` (2191 lines)

The single most important file for this project. Everything about how Q is parameterized,
optimized, normalized, and eigendecomposed lives here.

| Lines | Function | Why it matters |
|---|---|---|
| 43–72 | `ModelMarkov::ModelMarkov` | Allocation; default name `NonRev` / `Rev`. |
| 75–152 | `setReversible(bool, bool adapt_tree)` | Reallocates `rates`, converts a reversible half-matrix into a full non-reversible one, and **converts the tree unrooted↔rooted**. Non-reversible models require a rooted tree. |
| 154–159 | `getNumRateEntries()` | `n(n-1)/2` reversible, `n(n-1)` non-reversible. For amino acids: 190 vs **380**. |
| 299–351 | `init_state_freq(StateFreqType)` | How `+F / +FO / +FQ / +FU` populate `state_freq`. |
| 354–359 | `init(StateFreqType)` | `init_state_freq` then `decomposeRateMatrix`. |
| 847–875 | `getRateMatrix / setRateMatrix / setFullRateMatrix` | Conversion between a full Q and the packed `rates[]` array. `setFullRateMatrix` shows the exact off-diagonal packing order for the non-reversible case. |
| 875–893 | `getStateFrequency` | Returns `state_freq` renormalized to sum 1. This is the π the likelihood kernel uses at the root. |
| 894–931 | `setStateFrequency / adaptStateFrequency` | `adaptStateFrequency` rescales non-reversible `rates[]` by a new π, `q_ij *= π_new_j / π_old_j`. The `PartitionModel` constructor applies it to linked models with `FREQ_ESTIMATE` or `FREQ_EMPIRICAL` (`partitionmodel.cpp:116–165`). |
| 932–962 | `getQMatrix` | Non-reversible: memcpy of `rate_matrix`. Reversible: builds Q from `rates` and π. |
| **964–976** | **`getNDim()`** | **Non-reversible returns `num_params` only — π contributes no free dimensions.** Returns 0 when `fixed_parameters` is set. |
| 978–999 | `getNDimFreq()` | Degrees of freedom *not* counted in `getNDim()`, used for AIC/BIC. Returns `n-1` for `FREQ_EMPIRICAL`, 3 or 9 for codon types, 0 otherwise (including `FREQ_ESTIMATE` and `FREQ_USER_DEFINED`). |
| 1018–1036 | `setVariables()` | Packs model state into the BFGS vector (**1-indexed**). |
| 1038–1089 | `getVariables()` | Unpacks the BFGS vector into model state; returns a `changed` flag. |
| 1091–1118 | `targetFunk()` | `getVariables` → `decomposeRateMatrix` → `clearAllPartialLH` → `-computeLikelihood()`. Also the penalty guard returning `1.0e+30` when any π entry falls below `Params::min_state_freq`. |
| 1139–1164 | `setBounds()` | Box constraints handed to BFGS. |
| 1166–1244 | `optimizeParameters()` | The driver: allocate, `setVariables`, `setBounds`, `minimizeMultiDimen`, `getVariables`, re-decompose, recompute likelihood. The reference single-model optimization loop. |
| **1246–1387** | **`decomposeRateMatrixNonrev()`** | Virtual (`modelmarkov.h:345`). Builds full Q from `rates[]`, then at line 1267 calls `computeStateFreqFromQMatrix(rate_matrix, state_freq, num_states)` — i.e. *solves for π from Q*. Guarded by `if (freq_type != FREQ_USER_DEFINED || optimize_from_given_params)`. Then normalizes Q so that Σᵢ πᵢ(−Qᵢᵢ) = `total_num_subst`, then eigendecomposes (Eigen3 `EigenSolver` over a complex spectrum, or `eigensystem_nonrev`). |
| 1389–1603 | `decomposeRateMatrix()` | Dispatcher; reversible path symmetrizes and uses `SelfAdjointEigenSolver` with fallbacks to `decomposeRateMatrixRev()`. |
| 1637–1700 | `readRates(istream&)` | Reversible reads a triangle; **non-reversible reads the full matrix including the diagonal and throws unless every row sums to 0**. |
| 1798–1884 | `readParameters / readParametersString` | Detect reversibility from the sign of the first entry (negative ⇒ full Q ⇒ non-reversible), read rates then π, then re-derive π from Q and print a warning on mismatch greater than 1e-3. This is the only place the code compares a supplied π with the π solved from Q. |
| 1949–1951 | `setRates()` | Base implementation is `ASSERT(0)`, and the base `getVariables` never calls it. `ModelLieMarkov` overrides it and calls it from its own `getVariables`. |
| 1954–1963 | `getModelByName()` | Static factory, reachable only for `UNREST` and Lie-Markov names. |
| **2119–2131** | **`computeStateFreqFromQMatrix(Q, pi, n)`** | Free function. Solves `[1ᵀ; Q]ᵀ x = e₁` via `colPivHouseholderQr`, and asserts only that the result sums to 1 within 1e-4. The Q→π solve performed on every non-reversible decomposition unless π is pinned. |

**Load when:** always.

---

### `model/modelprotein.cpp` (1417 lines — lines 31–1105 are data tables)

Amino-acid models. Two distinct halves:

1. **Lines 31–1105:** `builtin_prot_models`, a raw-string NEXUS block holding every built-in AA
   matrix. The **non-reversible** ones are `NQ.PFAM` (line 915), `NQ.BIRD` (938), `NQ.INSECT`
   (961), `NQ.MAMMAL` (984), `NQ.PLANT` (1007), `NQ.YEAST` (1030). Each is a full 20×20 Q
   (negative diagonal, rows summing to 0) followed by a 20-entry π row.
2. **Lines 1106–1245:** `ModelProtein::init()` — the dispatcher turning a model-name string into
   parameters. The `GTR20` branch (1180–1206) and the **`NONREV` branch (1207–1231)** are the
   entry points. `NONREV` defaults to `FREQ_ESTIMATE`, seeds Q from LG, calls
   `setReversible(false)`, and sets `num_params = getNumRateEntries() - 1 = 379`.
   `--init-model` (`Params::model_name_init`) overrides the seed. A `+F{...}` vector is read
   after the seed conversion (lines 1239–1242), so the conversion uses LG's frequencies.

Also present: `rescaleRates()` (~1090), `readRates()` (1265+, which handles the protein
lower-triangle file convention), `getNameParams()`, and checkpointing.

**Load when:** adding a new AA model name, changing how `NONREV` is initialized, or emitting an
estimated matrix. Read only the line ranges you need — the data tables are huge.

---

### `model/modelunrest.cpp` + `model/modelunrest.h` (138 + 55 lines)

`ModelUnrest` — the general non-reversible model (`-m UNREST`), written for DNA but implemented
generically over `num_states`. **This is the smallest complete worked example of a
non-reversible model class in the codebase**: constructor, `validModelName`, `setBounds`,
`setRates` (called only from the constructor), `setStateFrequency` (deliberately a no-op), and
the three checkpoint methods. About 190 lines for both files.

**Load when:** creating any new `ModelMarkov` subclass. Always worth the tokens.

---

## Tier 1 — The constrained-parameterization precedents

These are the three places IQ-TREE already does "optimize a matrix inside a constrained
subspace".

### `model/modelliemarkov.h` (161) + `model/modelliemarkov.cpp` (2336)

Lie-Markov models, for 4-state DNA. A non-reversible Q
is written as a fixed linear combination of basis matrices, so the free parameters
(`double *model_parameters`) live in a *lower-dimensional subspace* of rate space rather than
in `rates[]` directly.

Read these and nothing else in the file:

- `setBounds()` — line 886. Bounds on the *parameters*, not on the rates.
- `setVariables()` — line 899; `getVariables()` — line 919. Note the non-reversible branch: it
  copies to/from `model_parameters` and calls `setRates()` whenever anything changed.
- `setBasis()` — line 1087. Where the constraint subspace is defined.
- **`setRates()` — line 1194.** The map from free parameters to a full `rates[]`, including a
  normalization step that keeps all off-diagonals non-negative.
- `restartParameters()` — overrides `Optimization::restartParameters` to escape boundary optima,
  a known failure mode for constrained non-reversible models.

**Load when:** studying how an existing class stores parameters separately from `rates[]` and
maps them in `getVariables`, or how a basis is shifted to reach a fixed π.

### `utils/tools.cpp` lines 7936–8090 and `utils/tools.h` lines 469–482, 3821–3855

The DNA constrained-frequency machinery: `freqsFromParams()`, `paramsFromFreqs()`,
`forceFreqsConform()`, `nFreqParams()`, `setBoundsForFreqType()`, and the `StateFreqType` enum
encoding constraints such as `+FRY` (π_A + π_G = ½ = π_C + π_T) and `+F1231` (π_C = π_T).
This is IQ-TREE's existing idiom for *"π is restricted to a linear subspace, so optimize in
reduced coordinates"*, which restricts π rather than Q.

**Load when:** working with frequency-type constraints and their command-line syntax.

### `model/modeldna.cpp` lines 422–600

`ModelDNA::getNDim / getVariables / setVariables` — shows how a subclass overrides variable
packing to honour a `param_spec` string (which rate entries are tied together or held fixed)
*and* a constrained `freq_type` at the same time. The canonical example of mapping `ndim` free
variables onto a larger `rates[]` array with shared and fixed entries.

**Load when:** working with tied or fixed rate entries.

---

## Tier 2 — Optimization and numerics

### `utils/optimization.h` (≈240) + `utils/optimization.cpp`

The `Optimization` base class that `ModelSubst` inherits. Everything optimizable in IQ-TREE
implements this interface.

- `getNDim()`, `targetFunk(double x[])`, `derivativeFunk()`, `restartParameters()`.
- `derivativeFunk()` — `optimization.cpp:916`, virtual at `optimization.h:145`. One-sided
  forward differences with step `1e-4 * |x|`. `PartitionModel` does not override it, so linked
  optimization uses this implementation on the summed objective; the only override in the
  inspected code is `PhyloTreeMixlen` (`tree/phylotreemixlen.h:165`).
- `minimizeMultiDimen(guess, ndim, lower, upper, bound_check, gtol, hessian)` — line 176. BFGS
  with numerical gradients (`dfpmin` / `lnsrch`). **The workhorse; called by every
  `optimizeParameters`.**
- `L_BFGS_B(nvar, vars, lower, upper, pgtol, maxit)` — line 195. Box-constrained L-BFGS-B.
  Currently commented out of `ModelMarkov::optimizeParameters` because of historical NaN issues
  (see the block comment at `modelmarkov.cpp:1197–1220`).
- `minimizeOneDimen`, `minimizeNewton`, `brent`, `dbrent` — one-dimensional routines.

**Critical convention:** every parameter vector is **1-indexed** (`variables[1..ndim]`);
`variables[0]` is unused.

**Note:** there is **no equality-constrained solver** in this codebase today — only box bounds.
Existing models that impose structure on Q do it through their parameterization.

**Load when:** working on optimization, gradients, or bounds.

### `lbfgsb/` (`lbfgsb.c`, `lbfgsb_new.h`)

Vendored L-BFGS-B translation used by `Optimization::L_BFGS_B`, supporting per-variable bound
types via `nbd[]`.

**Load when:** you need genuine box-constrained optimization or want to re-enable L-BFGS-B.

### `utils/eigendecomposition.h` (≈190) + `utils/eigendecomposition.cpp`

`EigenDecomposition`, the second base class of `ModelMarkov`.

- `eigensystem_nonrev(rate_matrix, state_freq, eval, eval_imag, evec, inv_evec, n)` — line 71.
  Real non-symmetric eigendecomposition (`elmhes` / `eltran` / `hqr2` / `luinverse`), used when
  `--eigen` selects `MET_EIGEN_DECOMPOSITION`. There is no `--matrix-exp` option; the flags are
  `--eigenlib`, `--eigen`, `--scaling-squaring`, `--lie-markov` (`utils/tools.cpp:5093-5108`),
  and the default is `MET_EIGEN3LIB_DECOMPOSITION` (`utils/tools.cpp:7533`).
- `eigensystem_sym` — the reversible path.
- `ZERO_FREQ = 1e-10` (line 24) — threshold below which a state is dropped from the matrix.
- `total_num_subst` (line 81) — the normalization target for Q.
- `ignore_state_freq` (line 84) — set true for non-reversible models; suppresses the
  `Q *= diag(π)` step.

**Load when:** touching `decomposeRateMatrixNonrev`, debugging complex eigenvalues, or handling
a non-diagonalizable Q (the `nondiagonalizable` flag triggers a scaled-squaring fallback).

### Eigen3 (external, header-only)

Located by `FindEigen3.cmake` / `-DEIGEN3_INCLUDE_DIR=...`. `MatrixXd`, `VectorXd`,
`EigenSolver`, `SelfAdjointEigenSolver`, `FullPivLU`, and `colPivHouseholderQr` are used
directly inside `modelmarkov.cpp`. Eigen is already a hard dependency of the build.

---

## Tier 3 — Model construction, CLI, orchestration, reporting

### `model/modelmixture.cpp` lines 3122–3266 — `createModel()`

The **central dispatcher** mapping a model-name string to a concrete `ModelSubst*`. Parses `+P`
(PoMo), `+E` (sequencing error), `NAME{params}` braces, `+FQ`, `+F{...}`, then branches on
`ModelMarkov::validModelName` and sequence type into `ModelBIN / ModelDNA / ModelProtein /
ModelCodon / ModelMorphology`. Declared at `model/modelmixture.h:29`.

**Load when:** working on how model names are dispatched. Every model name is dispatched from
here.

The rest of `modelmixture.cpp` (4780 lines) is profile-mixture machinery (C10–C60, `+Fmix`) plus
`builtin_mixmodels_definition`. Relevant background: profile mixtures are exactly the "many π
sharing one exchangeability matrix R" construction that has no naive nQ analogue — the
motivation for this project. Read `ModelMixture::initMixture` (3200–3600) only if mixtures enter
scope.

### `model/modelfactory.h` (326) + `model/modelfactory.cpp` (1857)

`ModelFactory` owns the `ModelSubst*` + `RateHeterogeneity*` pair for one tree and drives their
joint optimization.

- `readModelsDefinition(Params&)` — line 88. Loads `builtin_mixmodels_definition`,
  `builtin_prot_models`, and any `--mdef` file into a `ModelsBlock`.
- Constructor (~150–700): parses the full `-m` string, splits off `+I`, `+G`, `+R`, `+F...`,
  handles `Params::model_joint` (lines 209–227 and 285–288), then calls `createModel`.
- `optimizeParametersOnly()` — 1275. Alternating model / site-rate optimization.
- `optimizeAllParameters()` — 1331. Joint BFGS over model plus site-rate dimensions.
- `optimizeParameters()` — 1570. The outer loop alternating branch lengths and parameters.
- `getNParameters()` — 1253:
  `model->getNDim() + model->getNDimFreq() + site_rate->getNDim() + branch parameters`.
- `getNDim / targetFunk / setVariables / getVariables` — 1927–1950: how the model's variable
  block is concatenated with the rate-heterogeneity block.

**Load when:** working on parameter counts, or on where optimization is actually invoked.

### `model/partitionmodel.h` (209) + `model/partitionmodel.cpp` (922)

The **QMaker / nQMaker** layer: estimating one shared Q across many partitions (`--model-joint`,
`--link-model`). This is how the published NQ.* matrices were produced.

- Constructor lines 61–120: builds `linked_models`; handles `--init-model DIVMAT`, which would
  seed Q from the empirical divergence matrix via `setFullRateMatrix` but opens with
  `ASSERT(0 && "init_by_div_mat not working")` at line 96.
- Constructor lines 116–165: for linked models with `FREQ_ESTIMATE` or `FREQ_EMPIRICAL`, pools
  state counts, prints "Mean state frequencies", and calls `adaptStateFrequency` on every
  partition's model; other frequency types skip this.
- `getNDim()` — 295, delegates to the linked model.
- `targetFunk()` — 299–336: applies one parameter vector to *every* partition sharing the model
  name and sums log-likelihoods (OpenMP over partitions).
- `setVariables / getVariables` — 733–751.
- `optimizeLinkedModel()` — 753–840: the partition-wide analogue of
  `ModelMarkov::optimizeParameters`.
- `optimizeLinkedModels()` — 842–868: loops over each distinct linked model, fixes and unfixes
  parameters, checkpoints.

**Load when:** working on the `--model-joint` path, the nQMaker workflow.

### `utils/tools.h` (≈3900) and `utils/tools.cpp` (≈8000)

The global `Params` singleton and all CLI parsing. Huge; never read whole. Grep for the field
you need. Anchors:

- `enum StateFreqType` — `tools.h:469–482`.
- `Params::freq_type` (1698), `min_state_freq` (1705), `model_name_init` (1624),
  `gtr20_model` (1836), `optimize_linked_gtr` (1834), `optimize_from_given_params` (1850),
  `link_model` (2419), `model_joint` (2422), `matrix_exp_technique` (2510).
- Option parsing: `--init-model` (`tools.cpp:2896`), `--min-freq` (3156), `--link-model` (5018),
  `--model-joint` / `--link-partition` (5023).
- Default initialization of `Params` fields: the `tools.cpp:7238` and `7480` regions.
- `usage_iqtree()` help text, including the non-reversible model list: `tools.cpp:5978–6030`.
- Frequency helpers: `freqsFromParams` (7936) and, nearby, `paramsFromFreqs`,
  `forceFreqsConform`, `nFreqParams`, `setBoundsForFreqType`.

**Load when:** working on a CLI flag or a `Params` field. `-mset` and `-madd` lists are split on
commas by `convert_string_vec` (`tools.cpp:588`). **Adding an option takes four edits:** declare the field in `tools.h`, default it in the
initializer region of `tools.cpp`, parse it in the `parseArg` argument loop, and document it in
`usage_iqtree()`.

### `main/phyloanalysis.cpp` (≈5000)

Top-level analysis driver and all `.iqtree` report generation.

- QMaker / nQMaker citation block keyed on `params.model_joint` — lines 165–177.
- `reportModelSelection` — 307.
- **`reportNexusFile(ostream&, ModelSubst*, string part_name)` — 422–461.** Writes an estimated
  matrix back out as a NEXUS `model NAME = ...;` entry. The non-reversible branch prints the
  full Q via `getQMatrix` followed by **equal** frequencies rather than the model's π. Values are
  printed at 6 significant digits, and every model is labelled `GTRPMIX`.
- `reportLinkSubstMatrix` — 463; `reportModel` — 581 and 733; `reportRate` — 796.
- `reportPhyloAnalysis` — 1409; the `.GTRPMIX.nex` emission block — 2023–2060.
- `readModelsDefinition` call sites — 3585, 4288, 4857; `initializeModel` — 3586.

**Load when:** you need output, reporting, or a citation. Grep to the function; never read the
file whole.

### `main/phylotesting.cpp` (≈7300)

ModelFinder. Holds the model-name tables: `aa_model_names` (line 161),
**`aa_model_names_nonrev[] = {"NQ.bird", ...}` (line 165)**, `aa_mixture_model_names` (168),
`aa_freq_names*` (205–208), `aa_usual_nonrev_model = "NQ.pfam"` (224). `--model-joint` handling
at 1381. `readModelsDefinition` at 1429, 7082, 7340.

ModelFinder rejects candidate sets that mix recognized reversible and non-reversible names
(`mixRevNonrev`).

**Load when:** working on which model names ModelFinder or `-mset` accept.

### `nclextra/modelsblock.h` + `.cpp`

`ModelsBlock` and `NxsModel` (with a `description` string and a `flag` of `NM_ATOMIC` or
`NM_FREQ`). The NEXUS `begin models;` reader backing both the built-in matrices and user
`--mdef` files.

**Load when:** working with NEXUS model definitions (`--mdef` or the built-in matrices).

---

## Tier 4 — Read only if you touch the likelihood or the tree

### `model/modelsubst.h` (486) / `model/modelsubst.cpp` (233)

The abstract root of the model hierarchy. Holds `num_states`, `name`, `full_name`,
`fixed_parameters`, `state_freq`, `freq_type`. `isReversible()` (70), `useRevKernel()` (73),
`fixParameters()` (82), `getNumRateEntries()` (167). Small, cheap to read, and the definitive
list of what is virtual.

### `tree/phylotree.h` (≈2400) / `tree/phylotree.cpp` (≈7000)

`PhyloTree` owns `model`, `model_factory`, `site_rate` (fields at 2376 / 2382 / 2387).

- `getModel()` 523, `getModelFactory()` 527, `clearAllPartialLH()` 788, `computeLikelihood()`
  1051, `optimizeAllBranches()` 1452, `convertToRooted()` 2256, `convertToUnrooted()` 2261.
- Kernel dispatch by reversibility around `phylotree.cpp:2622`.

**Load when:** you change `decomposeRateMatrix` semantics or need to know when partial
likelihoods must be invalidated. **Any change to Q or π must be followed by
`decomposeRateMatrix()` and then `clearAllPartialLH()`, or the likelihood is silently stale.**

### `tree/phylokernelnonrev.h`, `tree/phylokernelnew.h`

The SIMD likelihood kernels. `phylokernelnew.h:963` (and `phylokernelnonrev.h:676–682`) is where
the *root* tip likelihood vector is filled from `model->getStateFrequency(...)` — that is,
**`state_freq` is literally the root distribution used in the non-reversible likelihood.**

**Load when:** confirming what π means numerically. Do not modify these casually: they are
macro-heavy, templated over state counts, and compiled once per instruction set.

### `tree/iqtree.cpp` (≈3300)

`IQTree::initializeModel()` — 1061: chooses `PartitionModel` / `PartitionModelPlen` /
`ModelFactory`. `IQTree::optimizeModelParameters()` — 2305: called during tree search.

### `alignment/alignment.h` / `.cpp`

`computeStateFreq()` (h:803) — empirical π; `convfreq()` (829) — clamps small frequencies;
`computeDivergenceMatrix()` (847) — the empirical divergence matrix used by
`--init-model DIVMAT`.

**Load when:** working with empirical frequencies or the divergence matrix.

### `utils/checkpoint.h`

`CKP_SAVE`, `CKP_RESTORE`, `CKP_ARRAY_SAVE`, `CKP_ARRAY_RESTORE` (lines 25–40); class
`Checkpoint` (69); `CheckpointFactory` (484). Every model must round-trip its parameters through
`startCheckpoint` / `saveCheckpoint` / `restoreCheckpoint`, or a resumed run silently loses the
new parameters.

**Load when:** adding any new persistent model state. Not optional for correctness.

---

## Tier 5 — Build, test, and data

| Path | Contents | Load when |
|---|---|---|
| `model/CMakeLists.txt` | The `add_library(model ...)` source list. **Every new `.cpp`/`.h` under `model/` must be added here or it will not be compiled in.** Note that `modelnonrev.cpp` and `modelnonrev.h` exist but are **0 bytes** and are *not* listed — dead placeholders; do not use them. | Adding a file. |
| `CMakeLists.txt` (root, ≈1150) | `add_subdirectory` list at 840–893; targets `iqtree3` (916–926) and `iqtree3-aa` (929, CMAPLE-AA only); `target_link_libraries` at 1019 and 1040; per-ISA kernel libraries at 900–908. Requires Eigen3 and Boost via `find_package`. | Build changes. |
| `test_scripts/test_configs.txt` | Matrix of alignments × options read by `gen_test_standard.py` (see `test_scripts/README`). | Adding regression coverage. |
| `test_scripts/test_iqtree.sh`, `verify_results.sh`, `test_data/expect_ans.txt` | The regression harness CI runs on Linux (`.ps1` variants also exist). Turtle DNA and protein analyses only; no command names `NONREV`, `NQ.*`, `GTR20`, `UNREST`, a Lie-Markov model, or `--model-joint`. | Running the upstream regression suite. |
| `example/example.phy`, `example/aa_example.phy` | Small DNA and AA alignments — fast smoke tests. | Every manual verification. |
| `example/models.nex` | Example user model-definition file in NEXUS format. | Loading a custom Q or π from a file. |
| `example/example.nex` | Example partition file for the DNA `example.phy`. `NONREV` is protein-only, so this pair cannot exercise the protein nQMaker path. | DNA partition testing. |
| `test_scripts/test_data/turtle_aa.fasta`, `turtle_aa.nex` | 16 taxa, 3 protein partitions. | Protein multi-partition (`--model-joint`) testing. |

---

## Tier 6 — Explicitly NOT relevant. Do not read these.

Listed so you do not spend tokens discovering they are irrelevant.

- **`cmaple/`** — a *separate, vendored* likelihood engine with its **own duplicate model
  hierarchy** (`cmaple/model/model_aa.cpp`, `modelbase.cpp`, including its own `NONREV` and
  `GTR20` handling). Greps for `NONREV` hit it constantly. **It is not the model layer that
  IQ-TREE's ML search uses.** Ignore unless the task explicitly says CMAPLE.
- `terrace/`, `terraphast/`, `terracetphast/` — phylogenetic terrace analysis.
- `pll/` — vendored PLL kernels, used only for PLL-based NNI.
- `ncl/` — vendored NEXUS Class Library (raw parser; `nclextra/` is the part that matters).
- `booster/` — transfer bootstrap expectation.
- `lsd2/` — least-squares dating.
- `whtest/` — Weiss–von Haeseler test of homogeneity.
- `pda/` — phylogenetic diversity analysis.
- `gsl/`, `sprng/`, `zlib-1.2.7/`, `yaml-cpp/`, `vectorclass/` — vendored third-party code.
- `nn/`, `nn_models/` — ONNX neural-network model selection.
- `simulator/` (AliSim) — relevant *only* if you want to simulate data under a known constrained
  nQ to validate parameter recovery. `simulator/alisimulator.cpp:369` shows how the AA data type
  is inferred from the `NONREV` / `GTR20` model names. Otherwise skip.
- `obsolete/` — dead code still linked for historical reasons.
- `mutsel_rust/` (9 files, 10,839 lines) and `utils/mutsel_wrapper.{cpp,h}` — a Rust
  mutation-selection subproject that arrived upstream in v3.1.4, gated off by default
  (`USE_MUTSEL`). **Read on 2026-09-15; see `docs/agent/AA_MODEL_INFERENCE.md` section 11.** It
  stays in Tier 6 for routine work: the model is reversible, per-site, fixed on the C++ side
  (`ModelSet` with `fixParameters(true)`), and fitted externally by AdamW under priors, so none
  of its machinery is reusable for a constrained nQ. It contains one closed-form construction of
  a reversible Q with a prescribed stationary distribution: the Halpern-Bruno form
  `Q_ij = M_ij g(f_j - f_i)` with `f_i = log π_i - log π^mut_i` and `g(x) = x/(1-e^-x)`
  (`mutsel_rust/src/model.rs:141-171`).
- `doc/html/`, `doc/latex/` — generated Doxygen output. Enormous and derived; read the headers
  instead.
- `lib/`, `libmac*/`, `liblinux_arm/` — prebuilt binary blobs.
- `test_scripts/iqtree2/` — prebuilt reference binaries for regression diffing, not source.
- `phylo-yaml/`, `main/timetree.cpp`, `main/terraceanalysis.cpp` — unrelated subsystems.

---

## Minimum viable context sets

**"Understand how nQ is currently inferred"** (≈3.5k lines):
`docs/agent/ARCHITECTURE.md` → `model/modelmarkov.h` → `modelmarkov.cpp` lines 964–1400 and
2094–2131 → `model/modelprotein.cpp` lines 1106–1245.

**"Add a new constrained non-reversible model class"** — add:
`model/modelunrest.{h,cpp}` (whole) → `model/modelliemarkov.cpp` lines 886–940 and 1087–1230 →
`model/modelmixture.cpp` lines 3122–3266 → `model/CMakeLists.txt`.

**"Wire it to the command line and the reports"** — add:
`utils/tools.h` 469–482 plus the specific `Params` fields → `utils/tools.cpp` parse sites 2896,
4998, 5003 and usage 5940–5990 → `main/phyloanalysis.cpp` 162–177 and 420–459 →
`main/phylotesting.cpp` 160–230.

**"Make it work across partitions (nQMaker workflow)"** — add:
`model/partitionmodel.cpp` 61–120, 294–336, 672–807 → `model/modelfactory.cpp` 208–226,
1258–1360.
