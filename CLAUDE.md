# IQ-TREE 3 fork: π-constrained non-reversible amino-acid models

Fork of `iqtree/iqtree3` (C++17, CMake). The project adds constrained maximum-likelihood
estimation of non-reversible amino-acid rate matrices under a supplied target stationary
frequency vector. Goal statement is not yet final: see `docs/agent/PLAN.md` and ask rather than
infer.

Project facts only. Personal working preferences live in `~/.claude/CLAUDE.md`; personal
project-specific notes go in `CLAUDE.local.md` (add it to `.gitignore` first).

## Read these before starting

1. `docs/agent/PLAN.md`, goal, scope, current state, next step. Read every task.
2. `docs/agent/DECISIONS.md`, settled choices not to be re-litigated. Read before changing an
   approach.
3. `docs/agent/ARCHITECTURE.md`, how the codebase is built and where the project's seams are.
   Read before writing code. Start at its "Baseline and staleness" section.
4. `docs/agent/FILE_INDEX.md`, tiered index of which file to open and when. Use instead of
   grepping the tree.
5. `CHANGELOG.md`, history including approaches that failed. <!-- TODO: file does not exist yet -->

`model/CLAUDE.md` loads automatically for work under `model/`, where most of this project lives.

## Commands

<!-- TODO: unverified in the current session. These come from a prior session's notes and
     should be re-run and corrected, and a test command added, before being trusted. -->

- Configure (Windows, MSVC, from `build/`):
  `"C:/Program Files/CMake/bin/cmake.exe" -DEIGEN3_INCLUDE_DIR="C:/ProgramData/chocolatey/lib/eigen/include/eigen3" -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..`
- Build: `"C:/Program Files/CMake/bin/cmake.exe" --build . --config Release`
- CMake is not on PATH and must be called by full path.
- Eigen3 and Boost are hard `find_package` dependencies and are not vendored.
  `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is required because vendored `zlib-1.2.7` declares a
  `cmake_minimum_required` that modern CMake rejects.
- Smoke tests: `iqtree3 -s example/aa_example.phy -m LG+G4`, then `-m NONREV`, then
  `-m NQ.pfam`. Partition path: `-s example/example.phy -p example/example.nex --model-joint NONREV`.
- Regression scripts live in `test_scripts/` (`test_iqtree.ps1`, `verify_results.ps1`).

## Layout

- `model/` substitution models. Almost all project work happens here.
- `tree/` topology, branch lengths, likelihood kernels. `alignment/` sequence data.
- `utils/` `Params` singleton, CLI parsing, BFGS, eigendecomposition, checkpointing.
- `main/` analysis driver and reports (`phyloanalysis.cpp`), ModelFinder (`phylotesting.cpp`).
- `docs/agent/` documents for agents working on this project. Everything this fork adds to the
  documentation tree goes here, so upstream merges stay clean.
- `doc/` is generated Doxygen output from upstream, not ours, and is not worth reading.

## Conventions that differ from defaults

- Optimizer parameter vectors are 1-indexed. Rate arrays are 0-indexed.
- `targetFunk` returns negative log-likelihood. Invalid regions return `1.0e+30` rather than
  throwing.
- After changing model parameters: `decomposeRateMatrix()`, then `clearAllPartialLH()`, then
  recompute. Skipping the second silently produces wrong likelihoods.
- Use `outError` / `outWarning` / `ASSERT`, not `throw` / `assert` / `std::cerr`.
- Arrays read by likelihood kernels use `aligned_alloc<double>`, never `new double[]`.
- Global options are read through `Params::getInstance()` at the point of use.
- Adding a CLI option takes four edits, all in `utils/`: declare in `tools.h`, default it in
  `parseArg`, parse it in `parseArg`, document it in `usage_iqtree()`.
- A new file under `model/` must be added to `model/CMakeLists.txt` by hand. There is no glob.

## Gotchas

- `cmaple/` is a separate vendored engine with its own duplicate model hierarchy, including its
  own `NONREV` and `GTR20`. Keyword searches land there constantly. It is not the model layer
  IQ-TREE's ML search uses.
- `model/modelnonrev.{cpp,h}` are 0-byte placeholders. `model/modelgtr.cpp` includes a header
  that does not exist. Neither is in the build.
- `NONREV` is a name string handled inside `ModelProtein::init`, not a class.
- `-m NONREV+F{...}` currently pins π and skips the stationarity solve, leaving Q unconstrained.
  It is a live path and its behaviour should not be changed without a decision entry.
- `--model-joint` is absent from `usage_iqtree()` despite being documented on the website.
- Do not edit vendored directories (`cmaple/`, `pll/`, `ncl/`, `zlib-1.2.7/`, `yaml-cpp/`,
  `gsl/`, `sprng/`, `vectorclass/`, `terraphast/`, `booster/`, `lsd2/`, `whtest/`, `pda/`,
  `obsolete/`). Changes there create merge conflicts against upstream for no benefit.

## Git

<!-- TODO: remote and branch strategy not yet decided. Fill in the baseline commit, the remote
     layout, and the working branch name once settled, and record the decision in
     docs/agent/DECISIONS.md. Until then, do not push and do not change remotes. -->

At the time of writing, `origin` points at upstream `iqtree/iqtree3`, not at a personal fork,
and the working tree is 64 commits behind it.

## Status

`docs/agent/PLAN.md` holds the current state and the next step. `CHANGELOG.md` holds history.
