# IQ-TREE 3 fork: π-constrained non-reversible amino-acid models

Fork of `iqtree/iqtree3` (C++17, CMake). The project adds constrained maximum-likelihood
estimation of non-reversible amino-acid rate matrices under a supplied target stationary
frequency vector. The goal is stated in `docs/agent/PLAN.md` (written 2026-09-23); the plan's
other sections are still open, so ask rather than infer.

Project facts only. Personal working preferences live in `~/.claude/CLAUDE.md`; personal
project-specific notes go in `CLAUDE.local.md` (add it to `.gitignore` first).

## Read these before starting

1. `docs/agent/PLAN.md`, goal, scope, current state, next step. Read every task.
2. `docs/agent/DECISIONS.md`, settled choices not to be re-litigated: programming decisions as
   numbered entries, and the design's mathematical decisions D01 to D10, which are revisable
   under their own change rule. Read before changing an approach.
3. `docs/agent/design/`, the mathematical and methodological design (start with its
   `README.md`; the unified synthesis governs conflicts). Read before planning or changing the
   method.
4. `docs/agent/ARCHITECTURE.md`, IQ-TREE's current programming architecture around substitution
   models. Read before writing code. Start at its "Baseline and staleness" section.
5. `docs/agent/AA_MODEL_INFERENCE.md`, the mathematics and algorithms of how IQ-TREE currently
   infers each class of amino-acid model. Read when the behaviour of an existing model matters.
6. `docs/agent/FILE_INDEX.md`, tiered index of which file to open and when. Use instead of
   grepping the tree.
7. `CHANGELOG.md`, history including approaches that failed. The top entry is the current state.
8. `AI_DISCLOSURE.md`, appended to whenever an AI contribution here is substantive.

Items 4 to 6 describe the code as it is and prescribe nothing. Instructions about what to build
come only from `PLAN.md`, `DECISIONS.md`, and the design documents.

`model/CLAUDE.md` loads automatically for work under `model/`, where most of this project lives.

## Commands

**The build does not use MSVC.** MSVC cannot compile this code (variable-length arrays and
OpenMP 3.0 loops in `tree/phylokernelnew.h`; see `docs/agent/ARCHITECTURE.md` section 17). The
build runs under WSL2 Debian with clang, from a clone inside the Linux filesystem.

**Proposed, not yet run to completion on this machine.** Replace this note with the verified
commands once the first build succeeds.

```bash
# one-time, needs sudo: compiler, OpenMP, Eigen, Boost (cmake, git and make are already present)
sudo apt-get update && sudo apt-get install -y clang libomp-dev libeigen3-dev libboost-dev
# one-time clone; --recursive because lsd2 and cmaple are git submodules the build needs
git clone --recursive https://github.com/petergoodman/iqtree3.git ~/BioInformatics/iqtree3
cd ~/BioInformatics/iqtree3 && git switch nq-constrained-pi
git remote add upstream https://github.com/iqtree/iqtree3.git
git remote set-url --push upstream DISABLED_NO_PUSH_TO_UPSTREAM
# from the repository root of the WSL clone
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
cmake --build build -j 4
```

- Delete `build/` before re-configuring after a failed configure: a bad cached value persists in
  `build/CMakeCache.txt` and defeats retries.
- `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is required because vendored `zlib-1.2.7`, `yaml-cpp` and
  `terraphast` declare a `cmake_minimum_required` that modern CMake rejects. If CMake is ever
  called from PowerShell, quote every `-D` token whose value contains a period.
- The `mutsel_rust` subproject is off unless `-DUSE_MUTSEL=ON`, so no Rust toolchain is needed.
  Upstream CI passes that flag; a development build should not.
- Eigen3 and Boost are hard `find_package` dependencies and are not vendored.
- `-j 4` rather than `-j`: the WSL VM has 6 GB of memory, and the per-instruction-set kernel
  files are large to compile.
- Smoke tests: `build/iqtree3 -s example/aa_example.phy -m LG+G4`, then `-m NONREV`, then
  `-m NQ.pfam`. Protein partition path:
  `-s test_scripts/test_data/turtle_aa.fasta -p test_scripts/test_data/turtle_aa.nex --model-joint NONREV`.
  (`example/example.phy` is DNA and cannot exercise `NONREV`.)
- Regression scripts live in `test_scripts/`; CI runs `test_iqtree.sh` then `verify_results.sh`.
  They do not exercise `NONREV`, `NQ.*`, `GTR20`, or `--model-joint`.

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

- `origin` is the fork `petergoodman/iqtree3`. `upstream` is `iqtree/iqtree3`, and its **push
  URL is deliberately disabled** (`DISABLED_NO_PUSH_TO_UPSTREAM`) so that nothing can reach the
  official repository by accident. Do not re-enable it; contributing back is done by opening a
  pull request from the fork on GitHub.
- Work happens on `nq-constrained-pi`. **Never commit to `master`**, which is kept as an exact
  mirror of upstream so that syncing stays a fast-forward.
- Baseline is `63c330d9` (upstream tag `v3.1.4`). The anchors in `docs/agent/` are verified
  against it.
- To sync: `git switch master; git merge --ff-only upstream/master; git switch nq-constrained-pi;
  git merge master`. Afterwards, re-verify the `docs/agent/` anchors for any file the merge
  touched, in the same change.

## Status

`docs/agent/PLAN.md` holds the current state and the next step. `CHANGELOG.md` holds history.
