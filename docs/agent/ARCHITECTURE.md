# IQ-TREE 3 — Architecture Map for Agents

**Audience:** AI coding agents working in IQ-TREE 3's substitution-model layer, in particular
on the inference of non-reversible amino-acid rate matrices (nQ).

**What this document is.** A description of IQ-TREE 3's programming architecture as it exists at
the baseline commit: the classes, their relationships, data layouts, control flow, and coding
conventions around substitution models and their inference, especially non-reversible ones. Its
purpose is to spare agents from re-deriving that architecture by searching the tree.

**What it is not.** It is not a plan, a design, or a decision record, and nothing in it is an
instruction about what this project should build or how. Where it says "must", it states a
constraint the current code imposes (for example, that a parameter vector is 1-indexed), not a
project choice. The project's design is in `docs/agent/design/`, its decisions in
`docs/agent/DECISIONS.md`, and its plan in `docs/agent/PLAN.md`. A passage here that reads as a
recommendation for the project is a defect; report it rather than follow it.

**Companions.** Three reference documents describe the current code, each from a different
angle, and none of them prescribes anything:

- `docs/agent/ARCHITECTURE.md` (this file): the programming architecture.
- `docs/agent/AA_MODEL_INFERENCE.md`: the mathematics, optimization, and algorithms by which each
  class of amino-acid model is inferred.
- `docs/agent/FILE_INDEX.md`: which file to open for which question, with line anchors, so that
  agents do not grep the whole repository.

Read this document first, then use the file index to locate specific code.

## Baseline and staleness

**Baseline: `63c330d9` (upstream tag `v3.1.4`).** All line anchors in this document and in
`FILE_INDEX.md` were re-verified against that commit on 2026-09-15; anchors added or changed on
2026-09-23 were read at the same commit. Names and structure drift far more slowly than line
numbers; if an anchor misses, grep for the symbol name nearby.

Repository layout as of that date: `origin` is the fork `petergoodman/iqtree3`, `upstream` is
`iqtree/iqtree3` with its push URL deliberately disabled, work happens on branch
`nq-constrained-pi`, and `master` is kept as an exact mirror of upstream.

To check whether a file has moved since the baseline:

```bash
git diff --stat 63c330d9 HEAD -- <path>
```

If it has, re-verify that file's anchors and update both documents in the same change.

### Arrived upstream in v3.1.4

Two additions landed in the 64 commits that preceded this baseline. The first has now been read;
the second has not, so treat any statement about it as unverified.

- `mutsel_rust/` (**9 files, 10,839 lines of Rust**, of which `data.rs` alone is 7,854) plus
  `utils/mutsel_wrapper.{cpp,h}`. **Read on 2026-09-15 and written up in
  `docs/agent/AA_MODEL_INFERENCE.md` section 11.** In summary: a Halpern-Bruno
  mutation-selection model with one shared reversible mutation matrix and one free
  log-frequency vector per site, where fitnesses are derived as `log π^s − log π^mut`, which
  makes `π^s` the exact stationary distribution of the site's Q. It is reversible, it is fitted
  by AdamW with automatic differentiation rather than BFGS, and it is a MAP estimate under
  three quadratic log-space priors, not an ML estimate. The subproject is gated off by default
  (the configure step prints `Compile Rust Dependencies: OFF`), because `USE_MUTSEL` is tested
  with `if (USE_MUTSEL STREQUAL "ON")` but never declared with `option()`, so no Rust toolchain
  is required for a standard build.
- `main/outstreambuf.h`, and a reworked `test_scripts/` measurement harness
  (`remeasure.{ps1,sh}`, `expect_memory.txt`, `expect_runtime.txt`).

`model/modelfactorymixlen.{cpp,h}` were **deleted** in the same range. `ModelFactoryMixlen` no
longer exists.

---

## 1. The 30-second version

IQ-TREE is a ~500k-line C++17 codebase, organized as a dozen static libraries linked into one
executable. For substitution-model work, only four of those libraries matter:

```
utils/      Params (global CLI state), Optimization (BFGS), EigenDecomposition, Checkpoint
alignment/  Alignment: patterns, state counts, empirical frequencies
model/      ModelSubst hierarchy (Q matrices) + RateHeterogeneity hierarchy (site rates)
tree/       PhyloTree / IQTree: topology, branch lengths, likelihood kernels
main/       phyloanalysis.cpp (driver + reports), phylotesting.cpp (ModelFinder)
```

Everything optimizable — a substitution model, a rate model, a tree — inherits the same
`Optimization` interface and is driven by the same BFGS routine. **Learn that one interface and
you can extend any of them.** It is described in §5.

---

## 2. Build topology

Root `CMakeLists.txt` adds each source directory as a static library (lines 840–893), then links
them into the `iqtree3` executable (line 1019):

```
iqtree3  ←  main  model  tree  alignment  utils  nclextra  ncl  pll  pda
            lbfgsb  whtest  sprng  vectorclass  gsl  simulator  yaml-cpp
            phyloYAML  kernelsse [kernelavx kernelfma kernelavx512]
```

Facts worth knowing before you add a file:

- `model/CMakeLists.txt` is a **hand-maintained explicit source list**. There is no glob. A new
  `model/*.cpp` that is not added there simply never compiles, and you get link errors that look
  like missing symbols rather than a missing file.
- `target_link_libraries(model utils)` — the `model` library may depend on `utils` and on headers
  from `tree/` and `alignment/` (include paths are global), but not on `main/`.
- The likelihood kernels are compiled **multiple times** into separate libraries
  (`kernelsse`, `kernelavx`, `kernelfma`, `kernelavx512`) with different instruction-set flags.
  Anything included by `tree/phylokernel*.h` gets compiled once per ISA.
- Hard external dependencies: **Eigen3** and **Boost** (`find_package`, fatal if absent). Eigen is
  used directly inside `model/modelmarkov.cpp`, for example by `computeStateFreqFromQMatrix`.
- `iqtree3-aa` is a CMAPLE-AA build variant, unrelated to amino-acid model work despite the name.

### Dead files — do not be misled

- `model/modelnonrev.cpp` and `model/modelnonrev.h` are **0 bytes** and are not in
  `model/CMakeLists.txt`. Non-reversible logic lives in `modelmarkov.cpp` behind the
  `is_reversible` flag, plus `modelunrest.*` and `modelliemarkov.*`.
- `model/modelgtr.cpp` (887 lines) is **not** in `model/CMakeLists.txt`, and it `#include`s
  `modelgtr.h`, **a header that does not exist in the repository**. It therefore cannot compile
  even if added to the build. It is an unbuilt historical copy of the reversible code that now
  lives in `modelmarkov.cpp`, and reading it will teach you things that are no longer true.
- `cmaple/` contains a complete **second** model hierarchy with its own `NONREV`, `GTR20`, and
  amino-acid handling. It is a different engine. Greps constantly land there. Ignore it.

### There is no non-reversible module, only a flag

Non-reversible inference works and is fully supported, but it has no home of its own. There is
no `ModelNonrev` class and no non-reversible source file. The capability is **98
`is_reversible` / `isReversible()` / `useRevKernel()` branch sites** spread across 25 live
files:

```
30  model/modelmarkov.cpp        7  main/phyloanalysis.cpp     4  tree/phylotreesse.cpp
10  model/modelpomo.cpp          6  model/modelmixture.cpp     4  model/modelliemarkov.cpp
 7  model/partitionmodel.cpp     3  model/modelprotein.cpp     ... 16 more files
```

`NONREV` is likewise **a name string, not a class.** The class is `ModelProtein`, which handles
every amino-acid model through an `if/else` chain over the name, and `NONREV` is one branch of
it at `modelprotein.cpp:1207`. `GTR20` is the identical construction with `is_reversible` left
true. By contrast `UNREST` **is** a class (`ModelUnrest`), registered through
`ModelMarkov::getModelByName`. Both idioms are available when adding a model.

The three shipped entry points for non-reversible amino-acid inference are `-m NONREV` (fit 379
free rates to one alignment), `-m NQ.pfam` and siblings (use a pre-estimated matrix), and
`--model-joint NONREV` (estimate one matrix across many alignments, the nQMaker workflow).

---

## 3. Runtime control flow

```
main()                                   main/main.cpp:2219
 └─ parseArg(argc, argv, Params)         utils/tools.cpp  — fills the Params singleton
 └─ runPhyloAnalysis(params, checkpoint) main/phyloanalysis.cpp
     ├─ read alignment                   → Alignment (patterns, num_states, empirical π)
     ├─ readModelsDefinition(params)      model/modelfactory.cpp:88  → ModelsBlock
     │     builtin_mixmodels_definition + builtin_prot_models + optional --mdef file
     ├─ runModelFinder (optional)        main/phylotesting.cpp  — picks the -m string
     ├─ iqtree->initializeModel(...)     tree/iqtree.cpp:1061
     │     └─ new ModelFactory / PartitionModel / PartitionModelPlen
     │           └─ createModel(model_str, models_block, freq_type, freq_params, tree)
     │                 model/modelmixture.cpp:3122   ← THE model dispatcher
     │                 └─ new ModelProtein(...)  →  ModelProtein::init()
     │                       modelprotein.cpp:1106  →  ModelMarkov::init(freq)
     │                             →  init_state_freq()  →  decomposeRateMatrix()
     ├─ tree search / optimization
     │     IQTree::optimizeModelParameters()     tree/iqtree.cpp:2305
     │       └─ ModelFactory::optimizeParameters()   modelfactory.cpp:1570
     │            ├─ model->optimizeParameters(eps)      ← ModelMarkov::optimizeParameters
     │            ├─ site_rate->optimizeParameters(eps)
     │            └─ tree->optimizeAllBranches()
     └─ reportPhyloAnalysis(...)         main/phyloanalysis.cpp:1409  → writes .iqtree
```

For partitioned (`-p` / `-q` / `-Q`) analyses the same flow runs with `PhyloSuperTree` and
`PartitionModel`, and an extra step — `PartitionModel::optimizeLinkedModels()` — optimizes a
single Q jointly across all partitions. **That is the QMaker / nQMaker path.**

---

## 4. The object graph

```
Params (singleton)  ── read by everything, everywhere, via Params::getInstance()
    │
Alignment ────────── patterns, num_states, seq_type, empirical frequencies
    │
PhyloTree (IQTree, PhyloSuperTree)
    ├── ModelFactory *model_factory      owns and drives the pair below
    │     ├── ModelSubst *model          the Q matrix
    │     └── RateHeterogeneity *site_rate   +I, +G, +R, +H
    ├── ModelSubst *model                (same pointer, cached)
    └── RateHeterogeneity *site_rate     (same pointer, cached)
```

`PhyloTree::getModel()` and `getModelFactory()` are the accessors (`tree/phylotree.h:523`, `527`).
A model holds a back-pointer `PhyloTree *phylo_tree` — the cycle is intentional, because
`targetFunk()` needs to ask the tree for a likelihood.

### Class hierarchies

```
Optimization                       utils/optimization.h
  ├── ModelSubst                   model/modelsubst.h      (also CheckpointFactory)
  │     └── ModelMarkov            model/modelmarkov.h     (also EigenDecomposition)
  │           ├── ModelBIN
  │           ├── ModelDNA ── ModelDNAError
  │           ├── ModelProtein                 ← amino acids: LG, WAG, Q.*, GTR20, NONREV, NQ.*
  │           ├── ModelCodon ── {Empirical, Parametric, SemiEmpirical}
  │           ├── ModelMorphology
  │           ├── ModelUnrest                  ← general non-reversible (UNREST)
  │           ├── ModelLieMarkov               ← constrained-subspace non-reversible
  │           ├── ModelPoMo ── ModelPoMoMixture
  │           ├── ModelMixture (virtual) ── ModelCodonMixture
  │           └── ModelSet                     ← site-specific models
  ├── RateHeterogeneity            model/rateheterogeneity.h
  │     └── RateGamma, RateInvar, RateFree, RateHeterotachy, RateKategory, ...
  ├── ModelFactory                 model/modelfactory.h    (also CheckpointFactory)
  │     └── PartitionModel ── PartitionModelPlen
  └── PhyloTree                    tree/phylotree.h        (also MTree, CheckpointFactory)
        └── IQTree ── PhyloSuperTree ── PhyloSuperTreePlen / PhyloSuperTreeUnlinked
```

**Note two things:** `ModelMixture` and `ModelSet` inherit from `ModelMarkov` *and* from
`vector<ModelMarkov*>` — a container that is also a model. And `PartitionModel` is a
`ModelFactory`, not a `ModelSubst`; it implements the same `Optimization` interface at the
partition level.

---

## 5. The Optimization contract

This is the single most important pattern in the codebase. Every optimizable object implements
it, and BFGS is the only solver.

```cpp
class Optimization {                                     // utils/optimization.h
    virtual int    getNDim();                            // number of free parameters
    virtual double targetFunk(double x[]);               // function to MINIMIZE
    virtual void   setVariables(double *variables);      // model state  -> x
    virtual bool   getVariables(double *variables);      // x -> model state; returns "changed"
    virtual void   setBounds(double *lo, double *hi, bool *bound_check);
    virtual bool   restartParameters(...);               // escape boundary optima (optional)

    double minimizeMultiDimen(guess, ndim, lower, upper, bound_check, gtol, hessian=nullptr);
    double L_BFGS_B(nvar, vars, lower, upper, pgtol, maxit);
};
```

### Hard conventions

1. **All parameter vectors are 1-indexed.** `variables[1] .. variables[ndim]`. `variables[0]` is
   allocated but unused. Every `memcpy` in this codebase is `variables+1`. Off-by-one here is the
   single most common bug when extending a model.
2. **`targetFunk` returns the NEGATIVE log-likelihood.** Optimizers minimize. `optimizeParameters`
   then returns `-score` to its caller, which expects a log-likelihood.
3. **`getVariables` must return whether anything changed**, and callers use it to skip expensive
   recomputation. Returning `true` unconditionally is correct but slow; returning `false`
   incorrectly produces a silently stale likelihood.
4. **After any change to model parameters you must, in order:** `decomposeRateMatrix()`, then
   `phylo_tree->clearAllPartialLH()`, then recompute the likelihood. `targetFunk` does exactly
   this (`modelmarkov.cpp:1091–1118`). Skipping `clearAllPartialLH` gives wrong numbers with no
   warning.
5. **Constraints are box constraints only.** `setBounds` fills `lower[]`, `upper[]`, and a
   `bound_check[]` flag that asks the optimizer to retry from a different start if the optimum
   lands on a boundary. There is **no equality-constraint solver in this codebase.** The existing
   models that impose structure on Q do it through their parameterization (`ModelLieMarkov`, and
   `ModelDNA` with a `param_spec`).
6. **Invalid regions are signalled by returning `1.0e+30` from `targetFunk`**, not by throwing.
   See the π-underflow guard at `modelmarkov.cpp:1104–1110`.

### The canonical `optimizeParameters` skeleton

`ModelMarkov::optimizeParameters` (`modelmarkov.cpp:1166–1244`) is the reference implementation
of the single-model optimization loop. Its shape:

```cpp
double MyModel::optimizeParameters(double gradient_epsilon) {
    if (fixed_parameters) return 0.0;
    int ndim = getNDim();
    if (ndim == 0) return 0.0;

    double *variables   = new double[ndim+1];   // 1-indexed
    double *lower_bound = new double[ndim+1];
    double *upper_bound = new double[ndim+1];
    bool   *bound_check = new bool[ndim+1];

    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);

    double score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound,
                                       bound_check, max(gradient_epsilon, TOL_RATE));

    bool changed = getVariables(variables);
    if (changed || score == -1.0e+30) {
        decomposeRateMatrix();
        phylo_tree->clearAllPartialLH();
        score = phylo_tree->computeLikelihood();
    }
    delete [] bound_check; delete [] lower_bound; delete [] upper_bound; delete [] variables;
    return score;
}
```

### Two parameterization styles

**Direct** (`ModelMarkov`, `ModelUnrest`, `ModelProtein`): the BFGS variables *are* entries of
`rates[]`. `setVariables` / `getVariables` are memcpys. `setRates()` is trivial or unused.

**Indirect** (`ModelLieMarkov`, `ModelDNA` with a `param_spec`): the BFGS variables are an
abstract `model_parameters[]` of lower dimension, and **`setRates()` maps them onto the full
`rates[]` array**. `getVariables` writes `model_parameters` and then calls `setRates()`. In
`ModelLieMarkov` the relevant methods are `setBounds` (861→886), `setVariables` (899),
`getVariables` (919, which calls `setRates()` only when a value changed), `setBasis` (1087), and
`setRates` (1194).

The base `ModelMarkov::getVariables` (`modelmarkov.cpp:1038–1089`) copies the variables straight
into `rates[]` and never calls `setRates()`, whose base implementation is `ASSERT(0)`.
`ModelUnrest` overrides `setRates()` to pin its last rate to 1 (`modelunrest.cpp:106–111`), but
it uses the base `getVariables`, so its override runs only at construction.

---

## 6. Model lifecycle: from a `-m` string to a likelihood

```
"NONREV+F+G4"
  │
  ├─ ModelFactory ctor         strips "+G4" → RateGamma; strips "+F" → StateFreqType
  │                            model/modelfactory.cpp ~150–700
  ├─ createModel("NONREV", models_block, freq_type, freq_params, tree)
  │     model/modelmixture.cpp:3122
  │     ├─ parses NAME{params}, +FQ, +F{...}, +P (PoMo), +E (seq. error)
  │     ├─ if ModelMarkov::validModelName(str)  → UNREST or Lie-Markov
  │     └─ else by seq_type                     → ModelProtein for SEQ_PROTEIN
  │
  ├─ ModelProtein::init("NONREV", ...)          model/modelprotein.cpp:1207
  │     ├─ look up name in ModelsBlock (built-in NEXUS matrices)
  │     ├─ NONREV branch: freq default FREQ_ESTIMATE; seed Q from LG (or --init-model)
  │     ├─ setReversible(false)   → reallocates rates[] to n(n-1)=380, ROOTS the tree
  │     └─ num_params = getNumRateEntries() - 1 = 379
  │
  └─ ModelMarkov::init(freq)                    model/modelmarkov.cpp:354
        ├─ init_state_freq(freq)                sets state_freq per +F / +FO / +FQ / +FU
        └─ decomposeRateMatrix()
              └─ decomposeRateMatrixNonrev()    model/modelmarkov.cpp:1246
                    ├─ unpack rates[] into full rate_matrix, diagonal = -rowsum
                    ├─ computeStateFreqFromQMatrix(rate_matrix, state_freq, n)   ← π SOLVED HERE
                    ├─ normalize so Σ πᵢ(−Qᵢᵢ) = total_num_subst
                    └─ eigendecompose (complex spectrum)
```

From then on, `computeTransMatrix(t, P)` exponentiates via the cached eigensystem, and the
likelihood kernel consumes `P` plus `state_freq` at the root.

---

## 7. Data layouts — get these exactly right

### `double *rates` — the packed rate parameters

| Mode | Length | Packing |
|---|---|---|
| Reversible | `n(n-1)/2` (AA: 190) | Upper triangle, row-major: `(0,1) (0,2) … (0,n-1) (1,2) …`. For DNA: A-C, A-G, A-T, C-G, C-T, G-T. |
| Non-reversible | `n(n-1)` (AA: **380**) | **All off-diagonal entries of Q, row-major, skipping the diagonal:** `(0,1) (0,2) … (0,n-1) (1,0) (1,2) … (n-1,n-2)`. |

The non-reversible packing is defined by the unpacking loop in `decomposeRateMatrixNonrev`
(`modelmarkov.cpp:1256–1263`) and mirrored in `setFullRateMatrix` (`858–873`):

```cpp
for (i = 0, k = 0; i < num_states; i++) {
    double *rate_row = rate_matrix + i*num_states;
    double row_sum = 0.0;
    for (j = 0; j < num_states; j++)
        if (j != i) row_sum += (rate_row[j] = rates[k++]);
    rate_row[i] = -row_sum;          // diagonal is derived, never stored in rates[]
}
```

**File-format trap:** built-in *protein* matrices in NEXUS are written as a **lower** triangle,
while `rates[]` is an **upper** triangle. `ModelProtein::readRates` (`modelprotein.cpp:1265+`)
performs the transposition with
`id = col*(2*num_states-col-1)/2 + (row-col-1)`. Non-reversible protein matrices
(`NQ.*`) are instead written as a **full 20×20 matrix with negative diagonal**, and
`ModelMarkov::readRates` (`1665–1698`) validates that each row sums to zero.

**Reversibility is auto-detected from the sign of the first number in the file**
(`readParameters`, `modelmarkov.cpp:1798`): negative ⇒ full Q ⇒ non-reversible.

### The other arrays

| Field | Size | Meaning |
|---|---|---|
| `rate_matrix` | `n*n` | Full Q, row-major, normalized. Non-reversible only. |
| `state_freq` (in `ModelSubst`) | `n` | π. **For non-reversible models this is an OUTPUT derived from Q**, not an input, except under `FREQ_USER_DEFINED` without `optimize_from_given_params` (§9). |
| `eigenvalues`, `eigenvalues_imag` | `n` | Real and imaginary parts (non-reversible spectra are complex). |
| `eigenvectors`, `inv_eigenvectors`, `inv_eigenvectors_transposed` | `n*n` | Cached eigensystem. |
| `ceval`, `cevec`, `cinv_evec` | `n`, `n*n`, `n*n` | `std::complex<double>` *aliases over the same memory* as the real arrays. Non-reversible only. |

All of these are allocated with `aligned_alloc<double>` / `ensure_aligned_allocated` for SIMD.
Use the same helpers; a plain `new double[]` will crash the AVX kernels.

### Counting parameters

- `num_params` = number of free rate parameters. For `NONREV` on amino acids: `380 - 1 = 379`.
  The last entry `rates[379]` is left at its seed value and acts as the scale anchor; global
  scale is fixed afterwards by the `total_num_subst` normalization.
- `getNDim()` = number of BFGS dimensions. **Non-reversible: `num_params`, with no frequency
  dimensions** (`modelmarkov.cpp:964–976`), because π is determined by Q.
- `getNDimFreq()` = degrees of freedom *not* in `getNDim()`, added by
  `ModelFactory::getNParameters` for AIC/BIC (`modelmarkov.cpp:978–1000`). Returns `n-1` for
  `FREQ_EMPIRICAL`, 3 or 9 for the codon frequency types, and 0 otherwise, including
  `FREQ_ESTIMATE` and `FREQ_USER_DEFINED`, and 0 whenever `fixed_parameters` is set.
- `PartitionModel::getNParameters` counts a linked model's parameters once across partitions
  (`partitionmodel.cpp:242–256`).

---

## 8. Reversible vs non-reversible — the complete difference

| Aspect | Reversible (Q) | Non-reversible (nQ) |
|---|---|---|
| Parameterization | R (exchangeabilities) ⊗ π | Q directly |
| `rates[]` length | `n(n-1)/2` | `n(n-1)` |
| π | A free parameter (or empirical/equal) | **Solved from Q** by `computeStateFreqFromQMatrix` |
| `getNDim()` | `num_params + (n-1)` if `FREQ_ESTIMATE` | `num_params` |
| Tree | Unrooted (`convertToUnrooted()`) | **Rooted** (`convertToRooted()`), forced in `setReversible` |
| Eigen spectrum | Real, symmetric after `diag(√π)` similarity | **Complex**; may be non-diagonalizable |
| Decomposition | `SelfAdjointEigenSolver`, or `eigensystem_sym` | `EigenSolver` / `eigensystem_nonrev`; fallback to scaled squaring when `nondiagonalizable` |
| Likelihood kernel | `computePartialLikelihood*` | `computeNonrevPartialLikelihood*` (`tree/phylokernelnonrev.h`); selected by `ModelSubst::useRevKernel()` |
| `ignore_state_freq` | false | **true** — suppresses the `Q *= diag(π)` step in the decomposition |
| Cost | baseline | ~2× memory and time; complex arithmetic |

`setReversible(bool, bool adapt_tree)` (`modelmarkov.cpp:75–152`) performs the whole switch,
including converting an existing reversible half-matrix into a full non-reversible one via
`Qᵢⱼ = Rᵢⱼ πⱼ`, using whatever `state_freq` holds at the moment of the call
(`modelmarkov.cpp:112–119`). **It also re-roots the tree as a side effect.** That is how
`-m NONREV` ends up with a rooted tree without anyone asking for one.

In the `NONREV` branch of `ModelProtein::init` (`modelprotein.cpp:1207–1231`) the LG seed is
read and converted by `setReversible(false)` first, so the conversion uses LG's frequencies. A
`+F{...}` vector in the model string is read afterwards (`readStateFreq`, then
`freq = FREQ_USER_DEFINED`, lines 1239–1242) and does not affect the converted seed.

---

## 9. Where π comes from, by code path

**For reversible models**, π is an independent parameter. `init_state_freq()` sets it from `+F`
(empirical), `+FQ` (equal), `+FU` (matrix-supplied), or `+FO` (`FREQ_ESTIMATE`, optimized as
`n-1` extra BFGS dimensions). Detailed balance guarantees *any* positive π is stationary for a
symmetric R — which is exactly why profile mixture models (C10–C60, `+Fmix`) can pair one R with
many π vectors.

**For non-reversible models**, π is derived. Inside `decomposeRateMatrixNonrev`
(`modelmarkov.cpp:1246`):

```cpp
// 1. initialize π to uniform, unless the user pinned it
if (freq_type != FREQ_USER_DEFINED || Params::getInstance().optimize_from_given_params)
    for (i = 0; i < num_states; i++) state_freq[i] = 1.0/num_states;

// 2. unpack rates[] into the full Q with derived diagonal
...

// 3. SOLVE for the stationary distribution
if (freq_type != FREQ_USER_DEFINED || Params::getInstance().optimize_from_given_params)
    computeStateFreqFromQMatrix(rate_matrix, state_freq, num_states);

// 4. normalize Q to the target substitution rate
for (i = 0, sum = 0.0; i < num_states; i++)
    sum -= rate_matrix[i*num_states+i] * state_freq[i];
double delta = total_num_subst / sum;
// ... Q *= delta

// 5. eigendecompose
```

and `computeStateFreqFromQMatrix` (`modelmarkov.cpp:2119–2131`) is:

```cpp
int computeStateFreqFromQMatrix(double Q[], double pi[], int n) {
    MatrixXd A(n+1, n);
    A.topRows(1).setOnes();                          // Σ πᵢ = 1
    A.bottomRows(n) = Map<MatrixXd>(Q, n, n);        // Qᵀ π = 0  (column-major Map transposes)
    VectorXd b(n+1); b.setZero(); b(0) = 1.0;
    Map<VectorXd> freq(pi, n);
    freq = A.colPivHouseholderQr().solve(b);         // least-squares solve of the overdetermined system
    ASSERT(fabs(freq.sum()-1.0) < 1e-4);
    return 0;
}
```

**π is then used as the root distribution in the likelihood.** The root leaf's partial-likelihood
vector is filled from `model->getStateFrequency(...)` (`tree/phylokernelnew.h:963`, and the
non-reversible kernel at `tree/phylokernelnonrev.h:676–682`), with no transition matrix applied
to it. The kernel uses whatever `state_freq` holds; nothing on the likelihood path checks that it
is stationary for Q.

**Other code paths that set or check π for a non-reversible model:**

1. `freq_type == FREQ_USER_DEFINED` with `optimize_from_given_params == false` **skips** both the
   uniform reset and the solve, leaving `state_freq` at the user's values while Q is fitted with
   nothing tying the two together. This is what `-m NONREV+F{...}` does. The result is a
   nonstationary-root model: the supplied vector is the root distribution and is, in general,
   not stationary for the fitted Q.
2. `ModelMarkov::adaptStateFrequency` (`modelmarkov.cpp:914–931`) rescales the non-reversible
   `rates[]` by `rates[k] *= freq[j]/state_freq[j]`, a one-shot re-weighting. The
   `PartitionModel` constructor calls it on every partition's copy of a linked model whose
   frequency type is `FREQ_ESTIMATE` or `FREQ_EMPIRICAL`, after pooling state counts across
   those partitions, and skips it for every other frequency type
   (`partitionmodel.cpp:116–165`). Under `NONREV+FO` the solve in the next decomposition then
   overwrites `state_freq`. The rescaled matrix is stationary at the new vector only when the
   input matrix satisfies detailed balance.
3. `readParameters` / `readParametersString` (`1798`, `1846`) re-derive π from a file-supplied Q
   and **print a warning** when it disagrees with the file's π by more than `1e-3`
   (`modelmarkov.cpp:1833–1842`). This is the only place the code compares a supplied π with the
   π solved from Q.

Parameter counts for a family of generators constrained to a supplied π are derived in
`docs/agent/AA_MODEL_INFERENCE.md` section 7 and in `docs/agent/design/`; they describe a model
the code does not yet contain, so they are not restated here.

### Lie-Markov models with a fixed frequency vector

`ModelLieMarkov::setBasis()` (`modelliemarkov.cpp:1087`), for 4-state DNA, shifts the model's
basis matrices so that every Q in their span has a supplied or empirical π as its stationary
distribution, and reduces the free-parameter count by the frequency degrees of freedom. When a
Lie-Markov model receives a user-supplied or empirical π:

```cpp
int bdf = BDF[model_num];
num_params = MODEL_PARAMS[model_num] - bdf;   // free parameters REDUCED by the frequency dof
init_state_freq(getFreqType());
double tau[3];
piToTau(state_freq, tau, symmetry);           // π mapped into frequency coordinates
...
unpermuted_rates[rate] += tau[tauIndex] * transformationMatrix[rate];   // basis shifted by τ
```

The comment above the transform tables (`modelliemarkov.cpp:262`) states their role:
*"Each shows how to modify a basis matrix to enforce a fixed base frequency vector."* A
feasibility check warns when the requested π is unreachable for that model family:
*"Model %s cannot achieve requested equilibrium base frequencies ... Instead it will use ..."*.
The `TRANSFORM_*` tables are hand-derived for 4 states and 37 named model families.

`ModelLieMarkov::setRates()` (`modelliemarkov.cpp:1194`) builds `rates[]` from
`model_parameters[]` as a combination of basis matrices, then keeps every rate non-negative by
taking `basis[0]`, the one basis matrix whose off-diagonals are all non-negative, as an anchor
and rescaling the parameter combination by the worst-case ratio across entries. Lie-Markov
models are additionally closed under the Lie bracket, a restriction separate from the frequency
constraint.

---

## 10. Normalization and scaling conventions

- `total_num_subst` (declared in `EigenDecomposition`, `utils/eigendecomposition.h:81`) is the
  target expected substitution rate. Default 1.0; mixture components set it per class
  (`modelmixture.cpp:3446`, `3519`).
- Q is scaled so that `Σᵢ πᵢ (−Qᵢᵢ) = total_num_subst` — one expected substitution per unit
  branch length. **This makes the overall scale of `rates[]` unidentifiable**, which is why
  `num_params = getNumRateEntries() - 1` and the last rate is pinned.
- Branch lengths enter as `evol_time = time / total_num_subst` in `computeTransMatrix`.
- `ZERO_FREQ = 1e-10` (`eigendecomposition.h:24`): states with π below this are **removed** from
  the matrix before decomposition and reinserted afterwards as identity rows
  (`modelmarkov.cpp:1293–1370`). Such states therefore drop out of the decomposition without a
  message.
- `Params::min_state_freq` (default set in `tools.cpp`, overridable with `--min-freq`) is the
  optimizer's lower bound on any frequency; `targetFunk` returns `1e+30` below it.

---

## 11. Checkpointing contract

IQ-TREE checkpoints continuously so runs can resume. Every model implements three methods and
they must be kept in sync with any new state, or a resumed run silently reverts.

```cpp
void MyModel::startCheckpoint()  { checkpoint->startStruct("MyModel"); }

void MyModel::saveCheckpoint() {
    startCheckpoint();
    if (!fixed_parameters) CKP_ARRAY_SAVE(getNumRateEntries(), rates);
    endCheckpoint();
    ModelMarkov::saveCheckpoint();          // chain to the base class LAST on save
}

void MyModel::restoreCheckpoint() {
    ModelMarkov::restoreCheckpoint();       // chain to the base class FIRST on restore
    startCheckpoint();
    if (!fixed_parameters) CKP_ARRAY_RESTORE(getNumRateEntries(), rates);
    endCheckpoint();
    decomposeRateMatrix();                  // rebuild derived state
    if (phylo_tree) phylo_tree->clearAllPartialLH();
}
```

`ModelUnrest` (`model/modelunrest.cpp:113–138`) is the reference. Macros live in
`utils/checkpoint.h:25–40`. **Note the ordering asymmetry — base-class call last on save, first
on restore — and the mandatory `decomposeRateMatrix()` + `clearAllPartialLH()` after restore.**

---

## 12. The partition / linked-model layer (QMaker and nQMaker)

`PartitionModel : public ModelFactory` (`model/partitionmodel.h:32`) estimates **one shared Q
across many partitions**, each with its own tree and branch lengths. This is how the published
`Q.*` matrices (QMaker) and `NQ.*` matrices (nQMaker) were estimated.

```
PartitionModel::optimizeLinkedModels()                partitionmodel.cpp:842
  for each distinct model name in linked_models:
      unfix that model's parameters in every partition
      optimizeLinkedModel()                           partitionmodel.cpp:753
        ├─ setVariables(v)          delegate to the representative model
        ├─ model->setBounds(...)
        ├─ minimizeMultiDimen(...)  → targetFunk
        │     PartitionModel::targetFunk              partitionmodel.cpp:299
        │       parallel-for over partitions sharing the model name:
        │         part_model->targetFunk(x)           each partition scores the SAME x
        │       return the sum
        ├─ getVariables(v)          push x into EVERY partition's model
        └─ decomposeRateMatrix() on each, clearAllPartialLH, recompute
      re-fix parameters; saveCheckpoint(); dump()
```

Key detail: **the same raw variable vector `x` is handed to every partition's model object.**
The linked optimization relies on identical `x` producing an identical model in every
partition; it does not check this.

**Where the gradient comes from.** `optimizeLinkedModel` calls `minimizeMultiDimen` on the
`PartitionModel` object itself, and `PartitionModel` does not override `derivativeFunk`, so the
finite-difference gradient is `Optimization::derivativeFunk` (`utils/optimization.cpp:916`)
evaluating `PartitionModel::targetFunk`, the summed objective. A `derivativeFunk` override in a
model class is therefore reached by single-model optimization (`ModelMarkov::optimizeParameters`)
but not by linked optimization. The only override of `derivativeFunk` outside `Optimization` in
the inspected code is `PhyloTreeMixlen` (`tree/phylotreemixlen.h:165`).

**Frequencies before the first linked update.** The `PartitionModel` constructor pools state
counts over the partitions of each linked model with `FREQ_ESTIMATE` or `FREQ_EMPIRICAL`, prints
them as "Mean state frequencies", and calls `adaptStateFrequency` on every partition's model
(§9, `partitionmodel.cpp:116–165`). Linked models with any other frequency type skip this block.

Driven by:

- `--link-model` (`tools.cpp:5018`) and `--model-joint MODEL` / `--link-partition`
  (`tools.cpp:5023`), which set `Params::link_model` and `Params::model_joint`.
- `--model-joint NONREV` is the nQMaker invocation. `phyloanalysis.cpp:165–177` chooses the
  QMaker vs nQMaker citation by testing whether `model_joint` contains `"NONREV"`.
- `--init-model NAME|FILE|DIVMAT` (`tools.cpp:2896`) seeds the matrix. The `DIVMAT` branch
  (`partitionmodel.cpp:61–120`) begins with `ASSERT(0 && "init_by_div_mat not working")` at line
  96. `ASSERT` compiles to nothing only when `NDEBUG` is defined (`utils/tools.h:60–68`), and
  the build does not define it for gcc or clang (§17), so in those builds the branch stops at
  the assertion (not run). Under `NDEBUG` it would run a partial version (see
  `AA_MODEL_INFERENCE.md` section 15).
- Output of a linked fit. The estimated matrix is printed in the `.iqtree` report
  (`reportModel`), and `PhyloSuperTree::printBestPartitionParams` (`tree/phylosupertree.cpp:1526`)
  writes `.best_model.nex` with each partition's `getModelNameParams(true)` string.
  `reportNexusFile` (`phyloanalysis.cpp:422–459`) is not used on this path: it is called only when
  `params.optimize_linked_gtr` is set (`--link-exchange-rates`, `phyloanalysis.cpp:2023`), and
  writes `.GTRPMIX.nex`. It prints matrix entries at 6 significant digits, labels every model
  `GTRPMIX`, and for a non-reversible model writes the full Q followed by a uniform `1/n`
  frequency line rather than the model's π.

---

## 13. Params and the CLI contract

`Params` is a mutable global singleton, reached from anywhere as `Params::getInstance()`. There
is no dependency injection; model code reads `Params` directly (for example
`Params::getInstance().min_state_freq` inside `targetFunk`). The codebase uses that convention
rather than threading options through constructors.

**Adding an option requires four edits, all in `utils/`:**

1. `utils/tools.h` — declare the field in `struct Params`.
2. `utils/tools.cpp`, the `parseArg` initializer region (~7238 and ~7480) — set its default.
3. `utils/tools.cpp`, the `parseArg` argument loop — parse it, following the neighbouring
   `strcmp(argv[cnt], "--flag") == 0` pattern, including the missing-argument error.
4. `utils/tools.cpp`, `usage_iqtree()` (~5978–6030) — document it in the help text.

Note that `--model-joint` is **absent from `usage_iqtree()`**, so it does not appear in `-h`
despite being documented on the project website.

Model-name and frequency-type strings are a second, parallel CLI surface: `-m NAME{params}+F...`
is parsed by `ModelFactory`'s constructor and then by `createModel`
(`modelmixture.cpp:3122–3220`). `StateFreqType` (`tools.h:469–482`) is the enum that carries
frequency-handling intent through the whole model layer.

---

## 14. Where a substitution model class is wired into the code

This section lists the places in the current code that a substitution model class is registered
in or read from. It describes the codebase's plumbing; it is not a plan for this project.

1. **Class and files.** Model classes subclass `ModelMarkov` and live in `model/`.
   `model/modelunrest.{h,cpp}` is the shortest complete example of a non-reversible class.
2. **Build list.** Every source file must appear in `model/CMakeLists.txt`, which has no glob.
3. **Name dispatch.** A name is reachable either through a `validModelName` clause reached from
   `ModelMarkov::getModelByName` (`modelmarkov.cpp:1954`), as `UNREST` and the Lie-Markov models
   are, or through the `if/else` chain in `ModelProtein::init` (`modelprotein.cpp:1106–1245`), as
   `GTR20` and `NONREV` are. Both routes start from `createModel` (`modelmixture.cpp:3122`).
4. **Optimization interface.** `getNDim`, `setVariables`, `getVariables`, `setBounds`, and, for
   an indirect parameterization, `setRates` (§5), with the 1-indexing and `changed`-flag
   conventions.
5. **Checkpointing.** `startCheckpoint`, `saveCheckpoint`, `restoreCheckpoint`, with the
   save/restore ordering of §11.
6. **Reporting.** `writeInfo(ostream&)` feeds the `.iqtree` report, and `getNameParams()` lets
   the model round-trip through its own name string.
7. **Parameter counts.** `getNDim()` and `getNDimFreq()` feed AIC, AICc and BIC (§7).
8. **ModelFinder.** Candidate names come from the tables in `main/phylotesting.cpp:160–230`.
9. **Output files.** `main/phyloanalysis.cpp` writes the report and the NEXUS model file
   (`reportNexusFile` at 422), and chooses a citation block at 162–177.

---

## 15. Conventions, idioms, and landmines

**Idioms to match:**

- `outError(msg)`, `outWarning(msg)`, `ASSERT(cond)` — never `throw`, `assert`, or `std::cerr`.
- `convertIntToString`, `convertDoubleToString`, `convert_double_with_distribution` from
  `utils/tools.h` for all string conversion.
- `verbose_mode >= VB_MED / VB_MAX / VB_DEBUG` gating for diagnostic output.
- `aligned_alloc<double>(n)` / `ensure_aligned_allocated(ptr, n)` / `aligned_free(ptr)` for any
  array the likelihood kernels will touch.
- Doxygen `/** @param ... @return ... */` comments on public methods.
- Dated change comments, for example `// BQM 2015-09-07: ...`. Follow the local style where you
  modify old code.
- 4-space indentation in newer code; older files use tabs. Match the file you are editing.

**Landmines:**

- **Off-by-one.** Parameter vectors are 1-indexed. Rate arrays are 0-indexed.
- **Forgetting `clearAllPartialLH()`** after changing Q or π — stale likelihoods, no warning.
- **Assuming π is an input for non-reversible models.** It is an output of
  `decomposeRateMatrixNonrev` on every evaluation, except under `FREQ_USER_DEFINED` without
  `optimize_from_given_params` (§9).
- **Assuming a non-reversible model has an unrooted tree.** `setReversible(false)` roots it.
- **Editing `cmaple/`** because a grep for `NONREV` landed there. Different engine.
- **Reading `model/modelgtr.cpp` or `model/modelnonrev.*`.** Not built.
- **`ModelMarkov::setRates()` is `ASSERT(0)`, and the base `getVariables` never calls it.** A
  class that maps its own parameters into `rates[]` through `setRates()` has to call it itself,
  as `ModelLieMarkov::getVariables` does (§5).
- **Whether `ASSERT` is active depends on `NDEBUG`,** which the gcc and clang Release flags do
  not define (§17). An `ASSERT` that a build compiles out does not check anything.
- **Complex eigen arrays alias the real ones.** `ceval`/`cevec`/`cinv_evec` point into the same
  allocation as `eigenvalues`/`eigenvectors`/`inv_eigenvectors`. Do not free or reallocate one
  without the other.
- **`nondiagonalizable`** must be set when the eigenvector matrix is singular; the transition-
  probability code then falls back to scaled squaring.
- **OpenMP.** `PartitionModel::targetFunk` runs partitions in parallel. Anything it touches must
  be thread-safe or per-partition.
- **`fixed_parameters`** short-circuits `getNDim()` to 0 and is toggled around linked-model
  optimization (`partitionmodel.cpp:842–868`).

---

## 16. Code locations that handle π, rates, and their optimization

Where the current code handles each concern that bears on π-constrained estimation. This is a
map of existing code for orientation; it implies no design.

| Concern | Location |
|---|---|
| Where π is solved from Q | `modelmarkov.cpp:1267` inside `decomposeRateMatrixNonrev` (virtual, `modelmarkov.h:345`) |
| The Q→π solver itself | `computeStateFreqFromQMatrix`, `modelmarkov.cpp:2119` |
| The switch that pins π and skips the solve | the `freq_type != FREQ_USER_DEFINED \|\| optimize_from_given_params` guards at `modelmarkov.cpp:1252` and `1266` |
| One-shot reweighting of Q by a new π | `adaptStateFrequency`, `modelmarkov.cpp:914`; called from the `PartitionModel` constructor, `partitionmodel.cpp:116–165` |
| Reversible-to-non-reversible seed conversion | `setReversible`, `modelmarkov.cpp:75–152` (uses the current `state_freq`, lines 112–119) |
| Dimension counts | `getNDim` (`964`) and `getNDimFreq` (`978`); linked count in `PartitionModel::getNParameters`, `partitionmodel.cpp:242` |
| Free-parameter ↔ rate-matrix map | `setVariables` (`1018`), `getVariables` (`1038`), `setRates` (`1949`, base asserts) |
| Basis shift to a fixed π (4-state DNA) | `ModelLieMarkov::setBasis`, `modelliemarkov.cpp:1087`, with `piToTau` at `1032` and `BDF[]` at `247` |
| Separate parameter vector mapped into `rates[]` | `ModelLieMarkov::getVariables` (`919`) and `setRates` (`1194`) |
| Nonstationary-root path `-m NONREV+F{...}` | `ModelProtein::init` sets `FREQ_USER_DEFINED` (`modelprotein.cpp:1239–1242`), so `decomposeRateMatrixNonrev` skips the stationarity solve while Q is unconstrained |
| Box constraints | `setBounds`, `modelmarkov.cpp:1139`; `ModelUnrest::setBounds`, `modelunrest.cpp:84`; bounds `MIN_RATE`, `MAX_RATE` in `modelmarkov.h` |
| Finite-difference gradient | `Optimization::derivativeFunk`, `utils/optimization.cpp:916` (virtual, `utils/optimization.h:145`) |
| Invalid-region signalling | `targetFunk`, `modelmarkov.cpp:1104–1110` |
| Where frequency vectors enter today | `+F{...}` in a model string via `createModel` (`modelmixture.cpp:3205–3215`) and `ModelFactory`; a NEXUS frequency vector via `ModelsBlock`; `--model-joint` strings are handed to each partition's `ModelFactory` |
| Candidate lists | `-mset` / `-madd` are split on commas by `convert_string_vec` (`utils/tools.cpp:588`), so a model string containing a comma cannot be passed there |
| Amino-acid model entry point | `ModelProtein::init`, `modelprotein.cpp:1207` (`NONREV` branch) |
| Multi-partition estimation | `PartitionModel::optimizeLinkedModel`, `partitionmodel.cpp:753` |
| Emitting an estimated matrix | linked fits: `.iqtree` report and `.best_model.nex` via `PhyloSuperTree::printBestPartitionParams` (`tree/phylosupertree.cpp:1526`); `reportNexusFile` (`phyloanalysis.cpp:422`) only for `--link-exchange-rates` (§12) |
| The existing π consistency check | the π-mismatch warning in `readParameters`, `modelmarkov.cpp:1833–1842` |

---

## 17. Building IQ-TREE 3

What is known about building this tree. The commands used on the project machine live in
`CLAUDE.md`; whether a build has succeeded there is recorded in `CHANGELOG.md`.

- **MSVC cannot compile this code.** `tree/phylokernelnew.h` declares variable-length arrays,
  for example `size_t mix_addr_nstates_malign[ncat_mix], ...` at line 2290 with `ncat_mix` a
  runtime value. These are a GCC and Clang extension that MSVC does not support. The same header
  parallelizes loops over `size_t` indices, which needs OpenMP 3.0, while MSVC implements
  OpenMP 2.0. MSVC can still complete the configure step, so the failure shows up only at
  compile time.
- **Upstream CI** (`.github/workflows/ci.yaml`) builds on Linux (Ubuntu 22.04, clang 22) and on
  Windows by cross-compiling with clang against a MinGW sysroot. It has no MSVC job. Its Linux
  job passes `-DUSE_MUTSEL=ON`, which needs a Rust toolchain, and `-DIQTREE_FLAGS=static`. The
  Rust subproject is off unless `USE_MUTSEL` is set to `ON` (`CMakeLists.txt:291`).
- **Dependencies.** Eigen3 and Boost are hard `find_package` requirements and are not vendored,
  unlike most third-party code here. OpenMP provides multithreading. When the compiler is clang
  on Linux, the `lld` linker is also required: `CMakeLists.txt:451–461` stops the configure with
  a fatal error if `ld.lld` is not found, because `cmaple` enables link-time optimization.
- **googletest is fetched at configure time.** `cmaple/CMakeLists.txt:281–293` downloads
  googletest at a pinned commit with `FetchContent` and builds `gtest` and `gtest_main`, which
  `cmaple/unittest/` uses. A configure therefore needs network access, and gtest targets exist in
  every build that integrates CMAPLE (the default).
- **`-DCMAKE_POLICY_VERSION_MINIMUM=3.5`** is needed only with CMake 4.x, which rejects projects
  declaring `cmake_minimum_required` below 3.5, as vendored `zlib-1.2.7` does (2.4.4). With
  CMake 3.25 and the system zlib (the configure prints "Using system zlib" when one is found),
  CMake reports the variable as unused.
- **Build type and `NDEBUG`.** The default build type is Release (`CMakeLists.txt:134–135`). For
  gcc and clang the root `CMakeLists.txt` replaces `CMAKE_CXX_FLAGS_RELEASE` with flags that do
  not include `-DNDEBUG` (lines 400 and 420), and no other part of the build defines `NDEBUG`
  (`ncl/ncl.h:72` defines it only for the Metrowerks compiler). Confirmed on 2026-09-23 from the
  configure output of a clang 14 Release build, whose CXX flags were
  `-std=c++17 -fopenmp -pthread -O3 -ffunction-sections -fdata-sections`: `ASSERT` is active in
  gcc and clang Release builds.
- **Line endings.** A Windows checkout with `core.autocrlf=true` has CRLF line endings, and the
  `test_scripts/*.sh` harness then fails under bash. A clone made inside Linux does not have
  this problem.
- **PowerShell quoting.** PowerShell splits the unquoted argument
  `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` into `-DCMAKE_POLICY_VERSION_MINIMUM=3` and `.5`, which
  is confirmable with `cmake -E echo -DCMAKE_POLICY_VERSION_MINIMUM=3.5`. The failure presents as
  `Invalid CMAKE_POLICY_VERSION_MINIMUM value "3"` reported once per vendored subproject, and the
  bad value is then cached in `build/CMakeCache.txt`, so a retry needs a fresh build directory.

Fast functional checks. None of these has yet been run on the project machine:

```bash
iqtree3 -s example/aa_example.phy -m LG+G4          # reversible, protein
iqtree3 -s example/aa_example.phy -m NONREV         # non-reversible, protein, rooted tree
iqtree3 -s example/aa_example.phy -m NQ.pfam        # built-in nQ matrix
iqtree3 -s example/example.phy -m UNREST            # smallest non-reversible model, DNA
iqtree3 -s test_scripts/test_data/turtle_aa.fasta -p test_scripts/test_data/turtle_aa.nex \
        --model-joint NONREV                        # nQMaker path, 3 protein partitions
```

`example/example.phy` is DNA, and `NONREV` is recognized only for protein data
(`modelprotein.cpp:1207`), so it cannot exercise the protein nQMaker path;
`test_scripts/test_data/turtle_aa.{fasta,nex}` (16 taxa, 3 protein partitions) can. The
`.iqtree` report prints the Q matrix and state frequencies. The upstream regression harness in
`test_scripts/` (`test_iqtree.sh`, `verify_results.sh`, expected values in
`test_data/expect_ans.txt`) runs DNA and protein analyses of the turtle data (tree search,
ModelFinder with default candidates, partitions, mixtures, concordance factors), compared at a
tolerance of 1 log-likelihood unit. None of its commands names `NONREV`, `NQ.*`, `GTR20`,
`UNREST`, a Lie-Markov model, or `--model-joint`.
