# AI disclosure

Record of substantive contributions by AI tools to this repository, kept so that the Methods or
Acknowledgements section of any resulting manuscript, and its cover letter, can state them
accurately. Newest first.

Entries record the tool and model version, what it did, and which files or sections it touched.
No AI-generated scientific claim, numerical result, or citation in this project has been accepted
without human verification against source.

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
in this session. The `NDEBUG` finding (that gcc and clang Release builds appear not to define it)
rests on reading `CMakeLists.txt`, not on a build. The WSL build procedure in `CLAUDE.md` is
proposed and has not been run.

**Files created.** `docs/agent/design/README.md`, the three copied design documents, and three
copied check scripts in `docs/agent/design/scripts/` (`verify_constrained_nq.py`, whose
docstring names the same Claude session as `constrained_nq_plan_agent.md`;
`check_builtin_matrices.py`, cited by that document but carrying no provenance line; and
`chart_fd_bfgs_compare.py`, whose generating tool is not recorded in the file). The scripts were
copied unmodified and not rerun.

**Files modified.** `CLAUDE.md`, `docs/agent/PLAN.md`, `docs/agent/DECISIONS.md`,
`docs/agent/ARCHITECTURE.md`, `docs/agent/FILE_INDEX.md`, `docs/agent/AA_MODEL_INFERENCE.md`,
`CHANGELOG.md`, `AI_DISCLOSURE.md`.

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
