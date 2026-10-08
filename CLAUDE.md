# IQ-TREE 3 fork: π-constrained non-reversible amino-acid models

Fork of `iqtree/iqtree3` (C++17, CMake). The project adds constrained maximum-likelihood
estimation of non-reversible amino-acid rate matrices under a supplied target stationary
frequency vector. The goal, scope, code plan and next step are in `docs/agent/PLAN.md` (code
plan approved 2026-09-23); where it is silent, ask rather than infer.

Project facts only. Personal working preferences live in `~/.claude/CLAUDE.md`; personal
project-specific notes go in `CLAUDE.local.md` (add it to `.gitignore` first).

## Read these before starting

Every session:

1. `docs/agent/PLAN.md`: goal, scope, constraints, open risks, current state, next step.
2. The top entry of `CHANGELOG.md`: where the last session stopped.

Only when the task needs them:

3. `docs/agent/CODE_PLAN.md`, after `PLAN.md`, when implementing or testing a slice: that
   slice's sections.
4. `docs/agent/DECISIONS.md`, through the index at its top: the entries a task cites, and any
   entry before changing an approach.
5. The slide decks in `docs/agent/design/presentation/`, the maintainer-approved statement of
   the method, before planning or changing the method. The longer design documents in
   `docs/agent/design/` (start with its `README.md`; the synthesis governs conflicts among them)
   only when the decks do not settle the question.
6. `docs/agent/ARCHITECTURE.md`, `docs/agent/AA_MODEL_INFERENCE.md` and
   `docs/agent/FILE_INDEX.md`: the relevant section, before touching IQ-TREE code not yet read in
   the session.

Not read at the start: `AI_DISCLOSURE.md` and older `CHANGELOG.md` entries.

The documents in item 6 describe the code as it is and prescribe nothing. Instructions about
what to build come only from `PLAN.md` with `CODE_PLAN.md`, `DECISIONS.md`, and the design
documents.

`model/CLAUDE.md` loads automatically for work under `model/`, where most of this project lives.

## Commands

**The build does not use MSVC.** MSVC cannot compile this code (variable-length arrays and
OpenMP 3.0 loops in `tree/phylokernelnew.h`; see `docs/agent/ARCHITECTURE.md` section 17). The
build runs under WSL2 Debian 12 with clang 14, reading the source from this Windows working copy
and writing build output to the Linux filesystem. There is one working copy; do not create a
second clone. Verified 2026-09-23 (unmodified code, `[100%] Built target iqtree3`).

```bash
# one-time, needs sudo (already done on this machine):
sudo apt-get update && sudo apt-get install -y clang lld libomp-dev libeigen3-dev libboost-dev
# in an interactive Debian terminal, turn conda off first: its base env ships its own Boost,
# OpenMP and libstdc++, which CMake can pick up instead of the system copies
conda deactivate
cmake -S /mnt/c/Users/peter/BioInformatics/iqtree3 -B ~/iqtree3-build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++
cmake --build ~/iqtree3-build -j 4
# binary: ~/iqtree3-build/iqtree3 (Linux path /home/petergoodman/iqtree3-build/iqtree3)
```

- `lld` is required: `CMakeLists.txt:451-461` stops the configure without it when the compiler is
  clang. The configure downloads googletest (via `cmaple`), so it needs network access.
- `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is unused with this machine's CMake 3.25 and system zlib;
  it is needed only with CMake 4.x. If CMake is ever called from PowerShell, quote every `-D`
  token whose value contains a period.
- The `mutsel_rust` subproject is off unless `-DUSE_MUTSEL=ON`, so no Rust toolchain is needed.
  Upstream CI passes that flag; a development build should not.
- `-j 4` rather than `-j`: the WSL VM has 6 GB of memory.
- Delete `~/iqtree3-build` before re-configuring after a failed configure, so that no stale value
  in its `CMakeCache.txt` survives.
- **Running WSL from Claude Code's Bash tool.** `wsl.exe` joins its arguments into one command
  line, so `$VAR` inside `wsl.exe -d Debian -- bash -c '...'` is expanded by the outer Linux shell
  and usually comes out empty. Put multi-step commands in a script file and run it with
  `MSYS_NO_PATHCONV=1 wsl.exe -d Debian -- bash /mnt/c/.../script.sh`; without
  `MSYS_NO_PATHCONV=1`, Git Bash rewrites `/mnt/c/...` into a Windows path. Shells started this way
  do not activate conda.
- Write run outputs outside the repository (for example `~/iqtree3-smoke`), never into the
  working copy.
- Smoke tests: `-s example/aa_example.phy` with `-m LG+G4`, `-m NQ.pfam`, `-m NONREV`;
  `-s example/example.phy -m UNREST` (DNA); protein partition path
  `-s test_scripts/test_data/turtle_aa.fasta -p test_scripts/test_data/turtle_aa.nex --model-joint NONREV`.
  `example/example.phy` is DNA and cannot exercise `NONREV`.
- Regression scripts live in `test_scripts/`; CI runs `test_iqtree.sh` then `verify_results.sh`.
  They do not exercise `NONREV`, `NQ.*`, `GTR20`, or `--model-joint`, and `test_iqtree.sh` needs
  `/usr/bin/time`, which this Debian does not have installed (nor `bc`).
- **Python runs only in the project's conda environment**, `iqtree3-fixedeq` in WSL (decision 013,
  built from `test_scripts/fixedeq/environment.yml`), never in the Windows or Debian system
  Python: `~/anaconda3/bin/conda run --no-capture-output -n iqtree3-fixedeq python ...`. Tests,
  from the repository root: `... python -m pytest test_scripts/fixedeq/tests -q -p no:cacheprovider`.
  Conda's FutureWarning about the `defaults` channel is harmless.
- **Regression driver**: `test_scripts/fixedeq/regression/regress.py` (`freeze`, `run`,
  `baseline`, `compare`; `CODE_PLAN.md` section 3.4). The baseline in
  `test_scripts/fixedeq/regression/baseline/` (runs 1 to 11) and run 12's in `baseline/run12/`
  were recorded from the frozen unmodified binary `~/iqtree3-baseline/iqtree3` (SHA-256 in
  `iqtree3.sha256` beside it), which is never rebuilt or overwritten; their raw outputs are in
  `~/iqtree3-runs/baseline/` and `~/iqtree3-runs/baseline/run12/`. After a source change, `run`
  the new build into a new directory and `compare` its `extract.json` with both files:
  `compare --baseline baseline/baseline.json --baseline baseline/run12/baseline.json`. Run 11 is
  judged by the decision 018 rule built into `compare`. Pass `--git git.exe`: Linux git inside WSL
  sees hundreds of false changes in this Windows working copy (line endings), while Windows
  `git.exe` reports the true state. `run`, `probes.py`, `differential.py` and `sanitize.py run`
  refuse an unclean working tree, untracked files included, unless given `--allow-dirty`. Use
  `--allow-dirty` for checks during development. Commit first only for a gate run whose result is
  recorded (the end-of-slice regression comparison), and do not create or edit repository files
  while it runs.
- **Sanitizer build** (decision 022):
  - `test_scripts/fixedeq/sanitizer/sanitize.py build --build-dir ~/iqtree3-build-asan --git
    git.exe` builds it, in about 30 minutes of compute.
  - `sanitize.py run --binary ~/iqtree3-build-asan/iqtree3 --build-dir ~/iqtree3-build-asan
    --out DIR --git git.exe` runs it. `--regression-runs`, `--probes` and `--no-differential`
    choose what runs, and `sanitize.py scan DIR` groups the reports.
  - The build runs about 20 times slower than Release. Every run so far shows one known upstream
    finding at `model/modelmarkov.cpp:76` (PLAN.md risk 20).
  - The unmodified sanitizer binary is frozen at `~/iqtree3-baseline-asan/iqtree3` (SHA-256 in
    `iqtree3.sha256`, with its `sanitizer-build.json` and `CMakeCache.txt`), never rebuilt or
    overwritten, so a later finding can be checked against unmodified code.
- **Probe driver**: `test_scripts/fixedeq/probes/probes.py --binary BIN --out DIR --git git.exe
  [--only a,b,...]` reruns the S0 runtime probes, in about 45 minutes of compute.
  `test_scripts/fixedeq/manifest/manifest.py --out test_scripts/fixedeq/manifest --git git.exe`
  regenerates the G0 manifest. The manifest is the frozen record of S0: later decisions
  supersede it, and it is not regenerated when they change.
- **Long jobs need the laptop awake.** Windows Modern Standby pauses the WSL VM and every job in
  it; idle sleep is set to 3 minutes on battery and 15 minutes plugged in. A job that looks hung
  is usually paused. Compare IQ-TREE's CPU time with its wall-clock time, the VM's `uptime`, and
  Windows Kernel-Power events 506 and 507 before assuming a fault.

## Layout

- `model/` substitution models. Almost all project work happens here.
- `tree/` topology, branch lengths, likelihood kernels. `alignment/` sequence data.
- `utils/` `Params` singleton, CLI parsing, BFGS, eigendecomposition, checkpointing.
- `main/` analysis driver and reports (`phyloanalysis.cpp`), ModelFinder (`phylotesting.cpp`).
- `docs/agent/` documents for agents working on this project. Everything this fork adds to the
  documentation tree goes here, so upstream merges stay clean.
- `test_scripts/fixedeq/` this project's Python environment, oracle, drivers, tests and
  regression baseline; the rest of `test_scripts/` is upstream's.
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
  official repository by accident. Do not re-enable it; contributing back is done by Peter
  opening a pull request from the fork on GitHub.
- **Never contact the IQ-TREE maintainers.** No agent opens issues or pull requests, comments,
  posts to discussions, or emails anyone at `iqtree/iqtree3`. Peter handles all communication
  with them; questions meant for the maintainers go to the open questions in
  `docs/agent/PLAN.md`. Development on the fork never waits on upstream.
- Work happens on `nq-constrained-pi`. **Never commit to `master`**, which is kept as an exact
  mirror of upstream so that syncing stays a fast-forward.
- Baseline is `63c330d9` (upstream tag `v3.1.4`). The anchors in `docs/agent/` are verified
  against it.
- To sync: `git switch master; git merge --ff-only upstream/master; git switch nq-constrained-pi;
  git merge master`. Afterwards, re-verify the `docs/agent/` anchors for any file the merge
  touched, in the same change.

## Status

`docs/agent/PLAN.md` holds the current state and the next step. `CHANGELOG.md` holds history.
