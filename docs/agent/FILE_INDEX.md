# IQ-TREE 3 — Context File Index for nQ-with-Target-π Work

**Audience:** AI coding agents working in this repository on *constrained optimization of
non-reversible amino-acid substitution matrices (nQ) under a user-supplied target stationary
frequency vector π*.

**How to use this file.** Do not `grep` the whole tree. Load files by tier. Tier 0 is mandatory
reading before writing any code. Tiers 1–4 are loaded on demand, keyed by the "Load when" line
of each entry. Tier 6 is an explicit *do-not-read* list — those directories look relevant by
name and are not.

Line numbers were verified against the working tree at commit `8977d31a`. They drift; treat
them as anchors to grep near, not as ground truth.

Companion document: **`AGENT_IQTREE_ARCHITECTURE.md`** — read that first for how the system
fits together, then use this file to find things.

---

## Tier 0 — Mandatory. Read these before writing any code.

### `model/modelmarkov.h` (558 lines)

The base class `ModelMarkov : public ModelSubst, public EigenDecomposition` for every
reversible *and* non-reversible Markov substitution model. This is the class the new
functionality will subclass or extend.

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
| 894–931 | `setStateFrequency / adaptStateFrequency` | `adaptStateFrequency` rescales non-reversible `rates[]` by a new π — an existing, closely related "make Q consistent with a given π" operation. Read it. |
| 932–962 | `getQMatrix` | Non-reversible: memcpy of `rate_matrix`. Reversible: builds Q from `rates` and π. |
| **964–976** | **`getNDim()`** | **Non-reversible returns `num_params` only — π contributes no free dimensions.** This is precisely the assumption the new feature changes. |
| 978–999 | `getNDimFreq()` | Degrees of freedom *not* counted in `getNDim()`, used for AIC/BIC. Returns 0 for `FREQ_ESTIMATE`. |
| 1018–1036 | `setVariables()` | Packs model state into the BFGS vector (**1-indexed**). |
| 1038–1089 | `getVariables()` | Unpacks the BFGS vector into model state; returns a `changed` flag. |
| 1091–1118 | `targetFunk()` | `getVariables` → `decomposeRateMatrix` → `clearAllPartialLH` → `-computeLikelihood()`. Also the penalty guard returning `1.0e+30` when any π entry falls below `Params::min_state_freq`. |
| 1139–1164 | `setBounds()` | Box constraints handed to BFGS. |
| 1166–1244 | `optimizeParameters()` | The driver: allocate, `setVariables`, `setBounds`, `minimizeMultiDimen`, `getVariables`, re-decompose, recompute likelihood. **The template any constrained optimizer must mirror.** |
| **1246–1387** | **`decomposeRateMatrixNonrev()`** | **The heart of the problem.** Builds full Q from `rates[]`, then at line 1267 calls `computeStateFreqFromQMatrix(rate_matrix, state_freq, num_states)` — i.e. *solves for π from Q*. Guarded by `if (freq_type != FREQ_USER_DEFINED || optimize_from_given_params)`. Then normalizes Q so that Σᵢ πᵢ(−Qᵢᵢ) = `total_num_subst`, then eigendecomposes (Eigen3 `EigenSolver` over a complex spectrum, or `eigensystem_nonrev`). |
| 1389–1603 | `decomposeRateMatrix()` | Dispatcher; reversible path symmetrizes and uses `SelfAdjointEigenSolver` with fallbacks to `decomposeRateMatrixRev()`. |
| 1637–1700 | `readRates(istream&)` | Reversible reads a triangle; **non-reversible reads the full matrix including the diagonal and throws unless every row sums to 0**. |
| 1798–1884 | `readParameters / readParametersString` | Detect reversibility from the sign of the first entry (negative ⇒ full Q ⇒ non-reversible), read rates then π, then re-derive π from Q and print a warning on mismatch. That warning is the current, weak version of the consistency condition the new feature must enforce *by construction*. |
| 1949–1951 | `setRates()` | Base implementation is `ASSERT(0)`. Subclasses using an indirect parameterization override it. |
| 1954–1963 | `getModelByName()` | Static factory, reachable only for `UNREST` and Lie-Markov names. |
| **2119–2131** | **`computeStateFreqFromQMatrix(Q, pi, n)`** | Free function. Solves `[1ᵀ; Q]ᵀ x = e₁` via `colPivHouseholderQr`. **The exact Q→π map the constraint must invert.** |

**Load when:** always.

---

### `model/modelprotein.cpp` (1417 lines — lines 31–1105 are data tables)

Amino-acid models. Two distinct halves:

1. **Lines 31–1105:** `builtin_prot_models`, a raw-string NEXUS block holding every built-in AA
   matrix. The **non-reversible** ones are `NQ.PFAM` (line 915), `NQ.BIRD` (938), `NQ.INSECT`
   (961), `NQ.MAMMAL` (984), `NQ.PLANT` (1007), `NQ.YEAST` (1030). Each is a full 20×20 Q
   (negative diagonal, rows summing to 0) followed by a 20-entry π row — exactly the format a
   newly estimated constrained nQ should be emitted in.
2. **Lines 1106–1245:** `ModelProtein::init()` — the dispatcher turning a model-name string into
   parameters. The `GTR20` branch (1180–1206) and the **`NONREV` branch (1207–1231)** are the
   entry points. `NONREV` defaults to `FREQ_ESTIMATE`, seeds Q from LG, calls
   `setReversible(false)`, and sets `num_params = getNumRateEntries() - 1 = 379`.
   `--init-model` (`Params::model_name_init`) overrides the seed.

Also present: `rescaleRates()` (~1090), `readRates()` (1265+, which handles the protein
lower-triangle file convention), `getNameParams()`, and checkpointing.

**Load when:** adding a new AA model name, changing how `NONREV` is initialized, or emitting an
estimated matrix. Read only the line ranges you need — the data tables are huge.

---

### `model/modelunrest.cpp` + `model/modelunrest.h` (138 + 55 lines)

`ModelUnrest` — the general non-reversible model (`-m UNREST`), written for DNA but implemented
generically over `num_states`. **This is the smallest complete worked example of a
non-reversible model class in the codebase** and is the right file to imitate for a new
subclass: constructor, `validModelName`, `setBounds`, `setRates`, `setStateFrequency`
(deliberately a no-op), and the three checkpoint methods. About 190 lines for both files.

**Load when:** creating any new `ModelMarkov` subclass. Always worth the tokens.

---

## Tier 1 — The constrained-parameterization precedents

These are the three places IQ-TREE already does "optimize a matrix inside a constrained
subspace". A new constrained optimizer should follow one of these rather than inventing a
fourth pattern.

### `model/modelliemarkov.h` (161) + `model/modelliemarkov.cpp` (2336)

Lie-Markov models. **The closest existing analogue to the target feature.** A non-reversible Q
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

**Load when:** designing the parameterization of the constrained nQ. High value.

### `utils/tools.cpp` lines 7898–8050 and `utils/tools.h` lines 469–482, 3821–3855

The DNA constrained-frequency machinery: `freqsFromParams()`, `paramsFromFreqs()`,
`forceFreqsConform()`, `nFreqParams()`, `setBoundsForFreqType()`, and the `StateFreqType` enum
encoding constraints such as `+FRY` (π_A + π_G = ½ = π_C + π_T) and `+F1231` (π_C = π_T).
This is IQ-TREE's existing idiom for *"π is restricted to a linear subspace, so optimize in
reduced coordinates"* — the mirror image of the target feature, which restricts Q given π.

**Load when:** deciding how to expose and validate a target-π constraint on the command line,
and how to reduce coordinates.

### `model/modeldna.cpp` lines 422–600

`ModelDNA::getNDim / getVariables / setVariables` — shows how a subclass overrides variable
packing to honour a `param_spec` string (which rate entries are tied together or held fixed)
*and* a constrained `freq_type` at the same time. The canonical example of mapping `ndim` free
variables onto a larger `rates[]` array with shared and fixed entries.

**Load when:** the new parameterization needs tied or fixed rate entries.

---

## Tier 2 — Optimization and numerics

### `utils/optimization.h` (≈240) + `utils/optimization.cpp`

The `Optimization` base class that `ModelSubst` inherits. Everything optimizable in IQ-TREE
implements this interface.

- `getNDim()`, `targetFunk(double x[])`, `derivativeFunk()`, `restartParameters()`.
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
Any hard constraint must be expressed by reparameterization, projection, or penalty.

**Load when:** choosing an optimizer or adding constraints.

### `lbfgsb/` (`lbfgsb.c`, `lbfgsb_new.h`)

Vendored L-BFGS-B translation used by `Optimization::L_BFGS_B`, supporting per-variable bound
types via `nbd[]`.

**Load when:** you need genuine box-constrained optimization or want to re-enable L-BFGS-B.

### `utils/eigendecomposition.h` (≈190) + `utils/eigendecomposition.cpp`

`EigenDecomposition`, the second base class of `ModelMarkov`.

- `eigensystem_nonrev(rate_matrix, state_freq, eval, eval_imag, evec, inv_evec, n)` — line 71.
  Real non-symmetric eigendecomposition (`elmhes` / `eltran` / `hqr2` / `luinverse`), used when
  `--matrix-exp` selects `MET_EIGEN_DECOMPOSITION`.
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
directly inside `modelmarkov.cpp`. **Eigen is already a hard dependency — use it for any new
linear algebra (null-space projection, QR, SVD) rather than adding a library.**

---

## Tier 3 — Model construction, CLI, orchestration, reporting

### `model/modelmixture.cpp` lines 3122–3266 — `createModel()`

The **central dispatcher** mapping a model-name string to a concrete `ModelSubst*`. Parses `+P`
(PoMo), `+E` (sequencing error), `NAME{params}` braces, `+FQ`, `+F{...}`, then branches on
`ModelMarkov::validModelName` and sequence type into `ModelBIN / ModelDNA / ModelProtein /
ModelCodon / ModelMorphology`. Declared at `model/modelmixture.h:29`.

**Load when:** registering a new model name. Any new model must be reachable from here.

The rest of `modelmixture.cpp` (4780 lines) is profile-mixture machinery (C10–C60, `+Fmix`) plus
`builtin_mixmodels_definition`. Relevant background: profile mixtures are exactly the "many π
sharing one exchangeability matrix R" construction that has no naive nQ analogue — the
motivation for this project. Read `ModelMixture::initMixture` (3200–3600) only if mixtures enter
scope.

### `model/modelfactory.h` (326) + `model/modelfactory.cpp` (1857)

`ModelFactory` owns the `ModelSubst*` + `RateHeterogeneity*` pair for one tree and drives their
joint optimization.

- `readModelsDefinition(Params&)` — line 87. Loads `builtin_mixmodels_definition`,
  `builtin_prot_models`, and any `--mdef` file into a `ModelsBlock`.
- Constructor (~150–700): parses the full `-m` string, splits off `+I`, `+G`, `+R`, `+F...`,
  handles `Params::model_joint` (lines 208–226 and 281–284), then calls `createModel`.
- `optimizeParametersOnly()` — 1258. Alternating model / site-rate optimization.
- `optimizeAllParameters()` — 1314. Joint BFGS over model plus site-rate dimensions.
- `optimizeParameters()` — 1553. The outer loop alternating branch lengths and parameters.
- `getNParameters()` — 1236:
  `model->getNDim() + model->getNDimFreq() + site_rate->getNDim() + branch parameters`.
- `getNDim / targetFunk / setVariables / getVariables` — 1834–1857: how the model's variable
  block is concatenated with the rate-heterogeneity block.

**Load when:** the new parameterization changes `getNDim()`, or you need to know where
optimization is actually invoked.

### `model/partitionmodel.h` (209) + `model/partitionmodel.cpp` (922)

The **QMaker / nQMaker** layer: estimating one shared Q across many partitions (`--model-joint`,
`--link-model`). This is how the published NQ.* matrices were produced and is the most likely
place an "estimate a constrained nQ from data" workflow plugs in.

- Constructor lines 61–120: builds `linked_models`; handles `--init-model DIVMAT`, which seeds Q
  from the empirical divergence matrix via `setFullRateMatrix`.
- `getNDim()` — 294, delegates to the linked model.
- `targetFunk()` — 299–336: applies one parameter vector to *every* partition sharing the model
  name and sums log-likelihoods (OpenMP over partitions).
- `setVariables / getVariables` — 672–690.
- `optimizeLinkedModel()` — 692–779: the partition-wide analogue of
  `ModelMarkov::optimizeParameters`.
- `optimizeLinkedModels()` — 781–807: loops over each distinct linked model, fixes and unfixes
  parameters, checkpoints.

**Load when:** the feature must work with `--model-joint NONREV` — very likely, since that is
the nQMaker workflow.

### `utils/tools.h` (≈3900) and `utils/tools.cpp` (≈8000)

The global `Params` singleton and all CLI parsing. Huge; never read whole. Grep for the field
you need. Anchors:

- `enum StateFreqType` — `tools.h:469–482`.
- `Params::freq_type` (1700), `min_state_freq` (1707), `model_name_init` (1624),
  `gtr20_model` (1838), `optimize_linked_gtr` (1836), `optimize_from_given_params` (1852),
  `link_model` (2421), `model_joint` (2424), `matrix_exp_technique` (2512).
- Option parsing: `--init-model` (`tools.cpp:2896`), `--min-freq` (3156), `--link-model` (4998),
  `--model-joint` / `--link-partition` (5003).
- Default initialization of `Params` fields: the `tools.cpp:7199` and `7442` regions.
- `usage_iqtree()` help text, including the non-reversible model list: `tools.cpp:5940–5990`.
- Frequency helpers: `freqsFromParams` (7898) and, nearby, `paramsFromFreqs`,
  `forceFreqsConform`, `nFreqParams`, `setBoundsForFreqType`.

**Load when:** adding a CLI flag (for example a `--target-freq`-style option) or a `Params`
field. **Adding an option takes four edits:** declare the field in `tools.h`, default it in the
initializer region of `tools.cpp`, parse it in the `parseArg` argument loop, and document it in
`usage_iqtree()`.

### `main/phyloanalysis.cpp` (≈5000)

Top-level analysis driver and all `.iqtree` report generation.

- QMaker / nQMaker citation block keyed on `params.model_joint` — lines 162–177. A new method
  should add its citation here.
- `reportModelSelection` — 305.
- **`reportNexusFile(ostream&, ModelSubst*, string part_name)` — 420–459.** Writes an estimated
  matrix back out as a NEXUS `model NAME = ...;` entry. The non-reversible branch prints the
  full Q via `getQMatrix` followed by **equal** frequencies — a wart worth knowing when emitting
  a constrained nQ whose π is the whole point.
- `reportLinkSubstMatrix` — 461; `reportModel` — 579 and 731; `reportRate` — 794.
- `reportPhyloAnalysis` — 1407; the `.GTRPMIX.nex` emission block — 2000–2037.
- `readModelsDefinition` call sites — 3526, 4268, 4845; `initializeModel` — 3573.

**Load when:** you need output, reporting, or a citation. Grep to the function; never read the
file whole.

### `main/phylotesting.cpp` (≈7300)

ModelFinder. Holds the model-name tables: `aa_model_names` (line 161),
**`aa_model_names_nonrev[] = {"NQ.bird", ...}` (line 165)**, `aa_mixture_model_names` (168),
`aa_freq_names*` (205–208), `aa_usual_nonrev_model = "NQ.pfam"` (224). `--model-joint` handling
at 1381. `readModelsDefinition` at 1430, 6998, 7256.

**Load when:** a new model name should be selectable by ModelFinder or `--mset`.

### `nclextra/modelsblock.h` + `.cpp`

`ModelsBlock` and `NxsModel` (with a `description` string and a `flag` of `NM_ATOMIC` or
`NM_FREQ`). The NEXUS `begin models;` reader backing both the built-in matrices and user
`--mdef` files.

**Load when:** you want to load a target π or a seed Q from a NEXUS file.

---

## Tier 4 — Read only if you touch the likelihood or the tree

### `model/modelsubst.h` (486) / `model/modelsubst.cpp` (233)

The abstract root of the model hierarchy. Holds `num_states`, `name`, `full_name`,
`fixed_parameters`, `state_freq`, `freq_type`. `isReversible()` (70), `useRevKernel()` (73),
`fixParameters()` (82), `getNumRateEntries()` (167). Small, cheap to read, and the definitive
list of what is virtual.

### `tree/phylotree.h` (≈2400) / `tree/phylotree.cpp` (≈7000)

`PhyloTree` owns `model`, `model_factory`, `site_rate` (fields at 2379 / 2385 / 2390).

- `getModel()` 527, `getModelFactory()` 531, `clearAllPartialLH()` 788, `computeLikelihood()`
  1051, `optimizeAllBranches()` 1455, `convertToRooted()` 2259, `convertToUnrooted()` 2264.
- Kernel dispatch by reversibility around `phylotree.cpp:2600–2620`.

**Load when:** you change `decomposeRateMatrix` semantics or need to know when partial
likelihoods must be invalidated. **Any change to Q or π must be followed by
`decomposeRateMatrix()` and then `clearAllPartialLH()`, or the likelihood is silently stale.**

### `tree/phylokernelnonrev.h`, `tree/phylokernelnew.h`

The SIMD likelihood kernels. `phylokernelnew.h:963` is where the *root* tip likelihood vector is
filled from `model->getStateFrequency(...)` — that is, **`state_freq` is literally the root
distribution used in the non-reversible likelihood.** That is the operational meaning of the
constraint this project wants to impose.

**Load when:** confirming what π means numerically. Do not modify these casually: they are
macro-heavy, templated over state counts, and compiled once per instruction set.

### `tree/iqtree.cpp` (≈3300)

`IQTree::initializeModel()` — 1063: chooses `PartitionModel` / `PartitionModelPlen` /
`ModelFactory`. `IQTree::optimizeModelParameters()` — 2307: called during tree search.

### `alignment/alignment.h` / `.cpp`

`computeStateFreq()` (h:803) — empirical π; `convfreq()` (829) — clamps small frequencies;
`computeDivergenceMatrix()` (847) — the empirical divergence matrix used by
`--init-model DIVMAT`.

**Load when:** you need an empirical π as a default target or as a seed.

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
| `test_scripts/test_configs.txt` | Matrix of alignments × options driven by `run_tests.sh` / `test_iqtree.ps1`. | Adding regression coverage. |
| `test_scripts/test_iqtree.ps1`, `verify_results.ps1` | Windows test drivers (this is a Windows dev machine). | Running tests locally. |
| `example/example.phy`, `example/aa_example.phy` | Small DNA and AA alignments — fast smoke tests. | Every manual verification. |
| `example/models.nex` | Example user model-definition file, in the exact NEXUS format a custom nQ must use. | Loading a custom Q or π from a file. |
| `example/example.nex` | Example partition file, for `--model-joint` testing. | nQMaker-style testing. |

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
- `doc/html/`, `doc/latex/` — generated Doxygen output. Enormous and derived; read the headers
  instead.
- `lib/`, `libmac*/`, `liblinux_arm/` — prebuilt binary blobs.
- `test_scripts/iqtree2/` — prebuilt reference binaries for regression diffing, not source.
- `phylo-yaml/`, `main/timetree.cpp`, `main/terraceanalysis.cpp` — unrelated subsystems.

---

## Minimum viable context sets

**"Understand how nQ is currently inferred"** (≈3.5k lines):
`AGENT_IQTREE_ARCHITECTURE.md` → `model/modelmarkov.h` → `modelmarkov.cpp` lines 964–1400 and
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
