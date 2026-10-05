# AI disclosure

Record of substantive contributions by AI tools to this repository, kept so that the Methods or
Acknowledgements section of any resulting manuscript, and its cover letter, can state them
accurately. Newest first.

Entries record the tool and model version, what it did, and which files or sections it touched.
No AI-generated scientific claim, numerical result, or citation in this project has been accepted
without human verification against source.

## 2026-10-04 (third session)

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5.5, 1M context, model id `claude-opus-5-5[1m]`.
**Operator:** Peter Goodman.

**What it did.** On a plan the operator approved, implemented decision 018 in the S0 regression
driver (the comparison rule for run 11 and the added run 12), wrote tests for it, recorded run 12
five times from the frozen unmodified binary, and wrote run 12's baseline to its own file without
changing the existing baseline. Then, on the operator's choice, recorded decision 019, wrote a
driver that reruns the three design scripts from fingerprint-checked copies outside the
repository and compares 133 values transcribed from the design documents, with tests, and ran
it. Twenty-four values were not reproduced, which were reported to the operator without
investigation. On the operator's acceptance of the rerun as the reproduced record, recorded
decision 020. Then, on an oracle plan the operator approved, recorded decision 021 (the oracle's
test thresholds) and corrected PLAN.md's description of IQ-TREE's frequency floor. It then wrote
the oracle in pieces: piece A is the chart core (`test_scripts/fixedeq/oracle/cases.py`,
`chart.py`, `t3.py`, `builtin.py`, with tests `test_oracle_chart.py`, `test_oracle_t3.py` and
`test_oracle_documented.py`). Piece B is the port of IQ-TREE's optimizer and the derivative rules
(`oracle/optimize.py`, `tests/test_oracle_optimize.py`).

**Verification.** The 42 tests in `test_scripts/fixedeq/tests/` passed in the project
environment. The new regression rule was checked on the five recorded run 11 repeats (failing
before the change, passing after). Every log-likelihood and count reported was produced by the
regression driver from IQ-TREE's output files, and every rerun value by the design scripts
themselves, compared by the rerun driver. The transcription of the 133 documented values was
done by the tool and has not been checked by the operator. IQ-TREE was run only from the
unmodified binary; no source code was changed.

**Files created.** `test_scripts/fixedeq/regression/baseline/run12/baseline.json` and
`baseline.md` (written by the driver), `test_scripts/fixedeq/design_rerun.py`,
`test_scripts/fixedeq/tests/test_design_rerun.py`.

**Files modified.** `test_scripts/fixedeq/regression/regress.py`,
`test_scripts/fixedeq/tests/test_regress.py`, `docs/agent/DECISIONS.md`, `docs/agent/CODE_PLAN.md`,
`docs/agent/PLAN.md`, `CLAUDE.md`, `CHANGELOG.md`, `AI_DISCLOSURE.md`.

**Source code modified.** None.

## 2026-10-03 (second session)

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5.5, 1M context, model id `claude-opus-5-5[1m]`.
**Operator:** Peter Goodman.

**What it did.** Explained the plan, its open items and slice S0 in chat. On the operator's three
choices, recorded decisions 016 (root policy for S0 to S4) and 017 (the oracle's port of IQ-TREE's
optimizer), added two S0 probes and four risks to the planning documents, and committed the
documentation work of the three preceding sessions unchanged. Wrote the conda environment file
for the oracle and created the environment in WSL, with its lock file. Then, on a plan the
operator approved, wrote the S0 regression driver and its tests, froze the unmodified IQ-TREE
binary, ran the 11 baseline runs (twice, and run 11 five times) and wrote the baseline summary.
Finally, on the operator's choice, recorded decision 018 (the comparison rule for run 11 and an
added run 12), audited the record for a handoff, and added the test setup's facts to
`CLAUDE.md`.

**Verification.** Every source anchor cited in the new text was read in this session at HEAD
`4c5f061f`. The package versions reported were printed by the created environment. Every
log-likelihood and count reported was produced by the driver from IQ-TREE's output files; the
driver's tests passed, and a dry run reproduced two values from the 2026-09-23 smoke tests. IQ-TREE
was run only from the unmodified binary; no source code was changed.

**Files created.** `test_scripts/fixedeq/environment.yml`, `test_scripts/fixedeq/environment.lock.txt`
(written by conda), `test_scripts/fixedeq/.gitattributes`, `test_scripts/fixedeq/.gitignore`,
`test_scripts/fixedeq/regression/regress.py`, `test_scripts/fixedeq/tests/test_regress.py`,
`test_scripts/fixedeq/regression/baseline/baseline.json` and `baseline.md` (written by the driver).

**Files modified.** `docs/agent/DECISIONS.md`, `docs/agent/CODE_PLAN.md`, `docs/agent/PLAN.md`,
`CLAUDE.md`, `CHANGELOG.md`, `AI_DISCLOSURE.md`.

**Source code modified.** None.

## 2026-10-03

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5.5, 1M context, model id `claude-opus-5-5[1m]`.
**Operator:** Peter Goodman.

**What it did.** Compared three placements of the coordinate-to-matrix build against the source,
recommended one, and on the operator's decision recorded it as programming decision 015, which
supersedes decision 004, with the matching updates to the plan and code plan. Explained C++
overriding and the nQ inference call path in chat.

**Verification.** Every source anchor cited was read in this session at HEAD `4c5f061f`. Only
read-only commands were run. IQ-TREE was not built or run, and no numerical result was produced.

**Files modified.** `docs/agent/DECISIONS.md`, `docs/agent/PLAN.md`, `docs/agent/CODE_PLAN.md`,
`CHANGELOG.md`, `AI_DISCLOSURE.md`.

**Source code modified.** None.

## 2026-10-01

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5.5, 1M context, model id `claude-opus-5-5[1m]`.
**Operator:** Peter Goodman.

**What it did.** Mapped the functions and advice from the operator's notes of a meeting with the
IQ-TREE maintainers to the source, compared them with the approved plan, and recorded the
comparison and three new risks (a transition-matrix method switch as a continuity source,
whether decision 004 is still needed, and the deferred optimizer edit of S5).

**Verification.** Every source anchor cited was read in this session at HEAD `4c5f061f`. Only
read-only commands were run. IQ-TREE was not built or run, and no numerical result was produced.

**Files modified.** `CHANGELOG.md`, `docs/agent/PLAN.md`, `AI_DISCLOSURE.md`.

**Source code modified.** None.

## 2026-09-23 (second session)

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5.5, 1M context, model id `claude-opus-5-5[1m]`, with one read-only search
subagent that stalled and contributed nothing.
**Operator:** Peter Goodman.

**What it did.** Wrote the project's code plan from the design documents and a reading of the
IQ-TREE source: the completed `docs/agent/PLAN.md`; a new companion, `docs/agent/CODE_PLAN.md`
(change map, class specification, test specification, per-slice lists, source anchors); and
thirteen programming decisions, entries 002 to 014 in `docs/agent/DECISIONS.md`, each with its
options, trade-offs and recommendation, all approved by the operator. Drafted the plan's
scientific motivation and success criteria for the operator's confirmation. Added to `CLAUDE.md`,
on the operator's instruction, the rule that agents never contact the IQ-TREE maintainers, and
added the companion document to its reading order. Later in the session, on the operator's
instruction, applied five minimal corrections of fact to the reference documents: a nonexistent
option name, two stale statements that the build had not been run, a missing Windows
qualification, and line anchors that pointed at the wrong code.

**Verification.** Every source anchor cited in the new text was read in this session in the
source at HEAD `4c5f061f`, which matches `63c330d9` outside documentation. Only read-only checks
were run (git state; conda and Python in WSL). IQ-TREE was not built or run, and no numerical
result was produced. Statements taken from documents are labelled as reported, and plan choices
as proposed.

**Files created.** `docs/agent/CODE_PLAN.md`.

**Files modified.** `docs/agent/PLAN.md`, `docs/agent/DECISIONS.md`, `CLAUDE.md`, `CHANGELOG.md`,
`AI_DISCLOSURE.md`, and, for the corrections of fact, `docs/agent/ARCHITECTURE.md`,
`docs/agent/AA_MODEL_INFERENCE.md` and `docs/agent/FILE_INDEX.md`.

**Source code modified.** None.

## 2026-09-23

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5.5, 1M context, model id `claude-opus-5-5[1m]`.
**Operator:** Peter Goodman.

**What it did.** Imported three AI-generated design documents into `docs/agent/design/` without
modification, and wrote a README recording their provenance. Their authorship, as stated in the
files: `constrained_nq_plan_agent.md` was prepared by Claude, model `claude-fable-5-1`, on
2026-09-22; `project2_iqtree3_technical_plan.md` (2026-09-22) and
`project2_unified_discrepancy_synthesis.md` (2026-09-23) do not name their generating tool or
model, which is to be supplied by the operator. Recorded one programming decision and the
synthesis's ten mathematical and methodological decisions in `docs/agent/DECISIONS.md`, on the
operator's instruction. Wrote the goal section of `docs/agent/PLAN.md` from those documents.
Corrected `docs/agent/ARCHITECTURE.md`, `docs/agent/FILE_INDEX.md` and
`docs/agent/AA_MODEL_INFERENCE.md` so that they describe current code without prescribing design,
and corrected statements of fact in them against source at commit `63c330d9`. Updated the reading
order and build section of `CLAUDE.md`.

**Verification.** Every source statement added or changed was read in the source at `63c330d9`
in this session. The operator built the unmodified branch in WSL2 with the procedure now in
`CLAUDE.md`; the tool then ran smoke tests of that binary (recorded in `CHANGELOG.md`, 2026-09-23)
and confirmed from the configure output that the Release build does not define `NDEBUG`. The
tool also changed the project's `.claude/settings.json` permission rule, on the operator's
instruction, so that it can run `git push origin nq-constrained-pi`.

**Files created.** `docs/agent/design/README.md`, the three copied design documents, and three
copied check scripts in `docs/agent/design/scripts/` (`verify_constrained_nq.py`, whose
docstring names the same Claude session as `constrained_nq_plan_agent.md`;
`check_builtin_matrices.py`, cited by that document but carrying no provenance line; and
`chart_fd_bfgs_compare.py`, whose generating tool is not recorded in the file). The scripts were
copied unmodified and not rerun.

**Files modified.** `CLAUDE.md`, `docs/agent/PLAN.md`, `docs/agent/DECISIONS.md`,
`docs/agent/ARCHITECTURE.md`, `docs/agent/FILE_INDEX.md`, `docs/agent/AA_MODEL_INFERENCE.md`,
`CHANGELOG.md`, `AI_DISCLOSURE.md`, `.claude/settings.json`.

**Source code modified.** None.

## 2026-09-15 (second session)

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5, 1M context, model id `claude-opus-5[1m]`, with eight subagent readers
of the same model family.
**Operator:** Peter Goodman.

**What it did.** Read the IQ-TREE 3 source tree and wrote `docs/agent/AA_MODEL_INFERENCE.md`, a
complete description of every method IQ-TREE 3 uses to infer amino-acid substitution matrices,
reversible and non-reversible, stated as mathematics alongside the implementing code. The
document covers the parameterisations and their identifiability, the likelihood as the kernels
compute it for unrooted reversible and rooted non-reversible models, rate heterogeneity, the
seven estimation procedures (fixed empirical matrix, `GTR20`, `NONREV`, QMaker and nQMaker joint
estimation across partitions, profile mixtures and GTRpmix, PMSF, MUTSEL), the optimiser and
eigendecomposition machinery, ModelFinder, parameter counting, output formats, and a list of
defects and dead code.

**How the reading was organised.** Eight subagents read separate subsystems in parallel and
returned reports with file and line anchors and verbatim code. The main agent read the core model
layer directly (`model/modelmarkov.{cpp,h}`, `model/modelprotein.cpp`, `model/modelunrest.cpp`,
`model/modelfactory.cpp`, `model/modelmixture.cpp` dispatcher and bounds,
`utils/eigendecomposition.cpp` rate-matrix construction, `utils/optimization.cpp`
`minimizeMultiDimen` / `derivativeFunk` / `restartParameters`, `alignment/alignment.cpp`
frequency computation, `tree/phylotree.cpp` parameter counting and root optimisation,
`tree/iqtree.cpp` orchestration, `main/phyloanalysis.cpp` reporting, and `utils/tools.cpp`
defaults and option parsing). Section 16 of the document records explicitly which claims were
read directly by the main agent, which rest on delegated reading, which are derivations, and
which were not verified at all.

**Independent checks performed on delegated claims.** Four claims that contradicted existing
project documentation or that were load-bearing were re-read directly and confirmed: the absence
of any `--matrix-exp` option; the size of `mutsel_rust/` (9 files, 10,839 lines, against the 13
files and ~13k lines recorded in `ARCHITECTURE.md`); the forward-difference gradient in
`Optimization::derivativeFunk` with `ERROR_X = 1e-4`; and the protein-specific
`bound_check = true` override in `ModelMixture::setBounds`.

**Claims that are derivations, not code facts, and are labelled as such in the document.** That
`pi^T Q = 0` has rank `n - 1` given zero row sums, and hence that a pi-constrained non-reversible
amino-acid model would have 360 free parameters nested inside the 379 of `NONREV`. That the
Halpern-Bruno fitness definition used in `mutsel_rust` makes the per-site frequency vector the
exact stationary distribution of the per-site Q; this one is corroborated by a symmetry assertion
in the Rust unit tests.

**Not verified.** No claim in the document rests on running IQ-TREE. The build was not exercised
in this session, so no numerical output was checked.

**Files created.** `docs/agent/AA_MODEL_INFERENCE.md`.

**Files modified.** `CHANGELOG.md` (new session entry), `AI_DISCLOSURE.md` (this entry).

**Source code modified.** None. No file under `model/`, `tree/`, `alignment/`, `utils/`, `main/`,
`mutsel_rust/`, or any vendored directory was altered.

## 2026-09-15

**Tool:** Claude Code (VS Code extension).
**Model:** Claude Opus 5, 1M context, model id `claude-opus-5[1m]`.
**Operator:** Peter Goodman.

**What it did.** Read the IQ-TREE 3 source tree and wrote reference documentation describing the
substitution-model layer, with emphasis on non-reversible model inference. Located and verified
the code paths by which a rate matrix is parameterized, optimized, normalized, and
eigendecomposed, and by which a stationary frequency vector is derived from a rate matrix.
Identified `ModelLieMarkov::setBasis` as an existing implementation of frequency-constrained
non-reversible estimation for 4-state DNA. Derived the parameter counts for constrained and
unconstrained 20-state non-reversible models. Advised on repository documentation structure and
on git remote and branch topology, and diagnosed a build failure as a PowerShell argument-quoting
problem.

All source-level claims were verified by reading the named files at commit `63c330d9`; line
references were re-checked against that commit. The mathematical statement that πᵀQ = 0 has rank
n-1 given zero row sums is a derivation, presented as such, and has not been independently
checked. The `mutsel_rust` subproject was identified but not read, and is labelled unverified in
the documentation.

**Files created.**

- `docs/agent/ARCHITECTURE.md`
- `docs/agent/FILE_INDEX.md`
- `docs/agent/PLAN.md` (placeholder, to be written by the operator)
- `docs/agent/DECISIONS.md` (placeholder, candidate decisions unratified)
- `CLAUDE.md`
- `model/CLAUDE.md`
- `CHANGELOG.md`
- `AI_DISCLOSURE.md`

**Files modified.** `.gitignore` (ignore rules for personal configuration files only).

**Source code modified.** None. No file under `model/`, `tree/`, `alignment/`, `utils/`,
`main/`, or any vendored directory was altered, apart from the addition of `model/CLAUDE.md`,
which is documentation and is not compiled.

**Commands run by the operator on the tool's recommendation.** Git remote reconfiguration, a
fast-forward sync of `master` to upstream tag `v3.1.4`, creation of branch `nq-constrained-pi`,
and a merge of upstream into that branch. No commits were pushed to the upstream IQ-TREE
repository, and the upstream push URL was deliberately disabled.
