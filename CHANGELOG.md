# Changelog

Session history for the π-constrained non-reversible amino-acid model project, newest first.
Each entry records what was done, what failed and why, and what is next. The next session starts
from the top entry.

This is a fork of `iqtree/iqtree3`. Entries here describe work on the fork, not upstream
development.

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

Later the same day Peter installed clang, libomp, Eigen and Boost in WSL2 Debian.

### Failed

Nothing attempted failed. No build was attempted and IQ-TREE was not run.

### Next

1. Build the unmodified branch in WSL2, run the smoke tests, and confirm from the configure
   output whether `NDEBUG` is defined. How the WSL build reaches the working tree (a build from
   the current Windows checkout, or a separate Linux clone) is awaiting Peter's decision; the
   procedure in `CLAUDE.md` is a proposal.
2. Supply the generating tool and model for the synthesis and for P-log, for `AI_DISCLOSURE.md`.
3. Then write the high-level code plan into `PLAN.md`.

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
