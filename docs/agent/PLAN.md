# Plan: π-constrained non-reversible amino-acid models

> **Status, 2026-10-03: code plan approved by Peter on 2026-09-23 and amended by decisions 015 to
> 017 on 2026-10-03; approach accepted by the IQ-TREE maintainers (reported by Peter, 2026-10-01); no
> source code changed yet.** Peter owns
> this document. Its companion, `docs/agent/CODE_PLAN.md`, holds the file-level change map, the
> NQC class specification and the test specification; it is subordinate to this document and is
> never read without reading this one first. Decisions are recorded only in
> `docs/agent/DECISIONS.md` and are cited here by number.

## What this document is for

This is the living core document for the project, rewritten as the project moves rather than
appended to. It answers four questions for any agent picking up work:

1. what we are building and why,
2. what constraints and requirements the work must satisfy,
3. where the project currently stands,
4. what the next step is.

It is read at the start of every task, and most tasks need nothing else from the planning
documents. A task that implements or tests a slice also reads `CODE_PLAN.md`, after this file.
Keep this file current and short enough to read in full, and put history in `CHANGELOG.md`.

## How it relates to the other documents

| Document | Lifecycle | Read when |
|---|---|---|
| `docs/agent/PLAN.md` (this file) | overwritten | start of every task |
| `docs/agent/CODE_PLAN.md` | overwritten together with this file | only after this file, when implementing or testing a slice |
| `docs/agent/design/` | fixed copies; a revision is a new file | before planning, or before changing the method |
| `docs/agent/DECISIONS.md` | append only, never edited | before changing a settled choice |
| `docs/agent/ARCHITECTURE.md` | corrected in place as findings land | before writing code |
| `docs/agent/AA_MODEL_INFERENCE.md` | corrected in place | when the method of an existing model matters |
| `docs/agent/FILE_INDEX.md` | corrected in place | when locating code |
| `CHANGELOG.md` | append only | to resume, or to check whether something was already tried |

Plans live in this file and in `CODE_PLAN.md`; decisions live only in `DECISIONS.md`
(programming decisions as numbered entries, the design's mathematical decisions as D01 to D10).
Neither planning document restates a decision's options or reasoning. `ARCHITECTURE.md`,
`AA_MODEL_INFERENCE.md`, and `FILE_INDEX.md` describe IQ-TREE's code as it currently is and
prescribe nothing.

Labels: (verified) means read in source at `63c330d9`, which HEAD `4c5f061f` matches outside
documentation, or run on 2026-09-23; (reported) means taken from a document; (proposed) means a
plan choice not yet executed.

## Goal

Add to IQ-TREE 3 a non-reversible amino-acid substitution model whose stationary distribution is
fixed to a supplied target vector π*, and estimate its rate matrix by maximum likelihood within
the existing nQMaker joint-estimation workflow (Dang et al. 2022, doi:10.1093/sysbio/syac007).
Given one strictly positive π* in IQ-TREE's state order, the model is the set of 20-state
generators Q with positive off-diagonal rates, zero row sums, π*Q = 0, and unit mean rate at π*.
It has 360 free parameters, contains every positive reversible generator with equilibrium π*
(189 parameters), and, as an unbounded model, is nested in the 379-parameter `NONREV` family.
The root distribution is π*, and π* is never re-solved from Q. The substantive change is
confined to nQMaker's shared-matrix update (step 3a): its optimizer variables become 360
jump-chain log-ratio coordinates (D01), and every trial matrix is built so that π*Q = 0 holds by
construction at every likelihood evaluation, not only at convergence. IQ-TREE's likelihood
engine, model selection, tree inference, and branch-length and rate-parameter optimization are
reused, with peripheral changes to target input, seeding, parameter counting, and output. The
estimate is the non-reversible counterpart of the reversible workflow of Wheeler et al. (2025;
doi:10.64898/2025.12.01.691663 as given in the Part 2 slide deck, not verified), which takes
exchangeabilities from filtered training alignments and composition from the target data, a
frequency replacement that is exact only for reversible models. The work is done on branch
`nq-constrained-pi` of the fork `petergoodman/iqtree3` and is intended for eventual submission
upstream as a pull request (`DECISIONS.md`, entry 001).

## Scientific motivation and success criteria

> Drafted from the design documents on 2026-09-23 for Peter's confirmation. The scientific
> criterion below is a placeholder until he replaces it.

The reversible workflow's frequency replacement is exact only because q_ij = R_ij π_j is
stationary at every positive π (reported). In an asymmetric factor the same replacement still
gives a valid generator but loses target stationarity (D09 item 6), and IQ-TREE's
`adaptStateFrequency` performs exactly that swap on a non-reversible matrix (verified,
`model/modelmarkov.cpp:914-930`). The constrained estimate instead refits Q inside the family
with π*Q = 0 and π* as the root distribution, so the fitted matrix is stationary at π* by
construction rather than transported there.

The intended study is the Project 2 cleaning study (reported, synthesis section 8 and D08): π*
is the pooled composition of the original training data, held fixed across cleaning treatments,
and the comparison arms are those of D08. Because π* e^{tQ} = π*, every tip marginal equals π*,
so a wrong target can produce apparent circulation; simulations with reversible and
non-reversible generators at correct and perturbed targets come before any circulation is read
as biology. Manuscript title and venue: to be supplied by Peter.

Success criteria for the software (proposed). The feature works when gates G0 to G3 of the
synthesis are met with this evidence:

- every accepted Q has ‖π*Q‖∞ and |mean rate − 1| at or below 1e-10, the synthesis's
  provisional goal, with the actual maxima reported;
- IQ-TREE's likelihood at a fixed Q equals the independent oracle's (test layer 3);
- a constrained fit seeded at the `GTR20+F{π*}` fit, with that fit's trees and rate parameters,
  never scores below it;
- 360 matrix parameters are counted once across linked partitions;
- π* is bitwise unchanged through construction, linking, optimization, checkpoint, restart and
  export;
- export followed by re-import reproduces the fitted likelihood;
- every legacy model in the regression baseline reproduces exactly.

Scientific success criterion (proposed placeholder). Frozen matrices are evaluated on the same
held-out genes under one declared nuisance policy across the D08 arms, with uncertainty from
resampling genes. Candidate criterion: the general refit at π* exceeds both the reversible refit
at π* and the frequency-replacement arm in held-out log-likelihood by more than the resampling
uncertainty, and non-reversibility is claimed only where it exceeds what reversible simulations
at a perturbed target produce. Whatever criterion is adopted is written as a confirmatory plan in
the repository before any training or evaluation data are read.

## Scope of the first release

The first release is gates G0 to G3 (slices S0 to S4 below), which confirms the synthesis's
default: Levels 1 and 2 of P-log section 7.4 through both `-m` and `--model-joint`, protein data
only, with fixed-matrix export and re-import. Level 1 fits the shared matrix with every nuisance
held fixed (nQMaker step 3a in isolation); Level 2 is IQ-TREE's native alternation of that fit
with branch-length and rate-parameter updates; Level 3 repeats the outer candidate reselection of
nQMaker steps 1 and 2.

In the first release:

- the NQC model (provisional public name; decision 002) with the log-ratio jump chart (D01),
  scaled finite-difference derivatives in both optimizer paths (D02, decision 005), and an
  explicit, recorded coordinate domain (D03, decision 007);
- LG rebuilt at π* as the default start, and `--init-model` starts (D05);
- invocation by `-m NQC+F{...}` or `-m NQC+F<name>` for one alignment, and by `--model-joint`
  under `-S` (separate trees) and `-p` (edge-proportional), the two arrangements nQMaker uses;
- checkpoint and restart; export of the fitted matrix at 17 significant digits and re-import
  through the existing `-m FILE` path (decision 010); post-fit diagnostics (decision 009);
- rejection, with a message, of every combination listed in decision 008.

Three refinements to the synthesis's default (proposed):

1. Positive-ratio build and inverse live in the pure module from S1, because D01 retains the
   positive chart and they cost little; the compiled positive benchmark comes in S5.
2. The compiled T3 backend (D04) comes after the first release but before any scientific
   comparison, matching D04's "after the log-jump core is validated".
3. Level 1 uses existing flags if S0 shows they freeze every nuisance: under `-S`, `-te` with
   `-blfix` and fixed rate parameters. Under `-p`, partition rates are re-optimized in every
   round (verified, `model/partitionmodelplen.cpp:150-157`), so Level 1 there needs `-q` or a
   small opt-in option; the S0 result decides.

Tree search with NQC runs through ordinary IQ-TREE code but is not a first-release validation
target.

Deferred, still in the project: the positive-ratio compiled benchmark, the compiled T3 backend,
a small direct-flux oracle (D10) and recording-only line-search status (S5, gate G4);
target-rebuilt guide candidates, fixed learned candidates in ModelFinder, the
reversible/non-reversible mixing guard, an optional Level 3 outer driver, and simulation (S6,
gate G5); domain sensitivity, multiple starts, mismatch simulations and held-out evaluation (S7,
gate G6).

Out of scope: EM (D10), exact boundary fitting, branch-specific targets, non-stationary roots,
profile mixtures with distinct equilibria, and analytic likelihood gradients.

## Constraints and requirements

1. Stationarity holds by construction at every likelihood evaluation, including line-search
   trials and finite-difference probes. The log-ratio chart is total on R^360, so an unclamped
   probe still gives a valid Q.
2. Off-diagonal rates stay strictly positive. Softmax underflow, a non-finite value or a
   non-positive ν rejects the proposal: `targetFunk` returns 1e30 and no state is published.
3. π* is immutable. It is never re-solved from Q, floored, renormalized after construction or
   replaced by any setter, and the root distribution equals π* bitwise.
4. Existing behaviour is unchanged by default. Every change to a shared IQ-TREE file is opt-in,
   by model type or by option, or behaviour-preserving, and the regression baseline must
   reproduce exactly at `-T 1`.
5. Parameter counts are correct: 360 when fitted, counted once across linked partitions; 0 when
   fixed; 0 frequency degrees of freedom for the external target.
6. The model's checkpoint is lossless and an incompatible restart is rejected. This matters
   beyond restarts: the default +I+G optimization path saves and restores the model through its
   checkpoint (verified, `model/modelfactory.cpp:1395-1505`; `opt_gammai` defaults to true,
   `utils/tools.cpp:7082`), and IQ-TREE writes checkpoint doubles at 10 significant digits
   (verified, `utils/checkpoint.h:22`).
7. Work inside `PartitionModel::targetFunk`'s OpenMP loop (verified,
   `model/partitionmodel.cpp:310-333`) is thread-safe: per-object workspace, no static mutable
   state, no shared warm starts, and the existing fixed-order sum preserved.
8. The coordinate domain is explicit and recorded, every start is representable with margin,
   and an imported matrix is never clipped (D03).
9. `-m NONREV+F{...}` keeps its current non-stationary-root meaning.
10. Code is ready for an eventual upstream pull request (entry 001): small diffs to shared files
    in the local style, IQ-TREE's conventions (`outError`, 1-indexed optimizer vectors,
    `aligned_alloc` for kernel arrays), and every new file listed in `model/CMakeLists.txt`.
11. **Hard rule:** no agent contacts the IQ-TREE maintainers or posts to `iqtree/iqtree3`
    (issues, pull requests, discussions, comments, email). Peter handles all communication;
    questions for the maintainers go to the open questions below.
12. **Hard rule:** all development happens on branch `nq-constrained-pi` of the fork, and none
    of it waits on upstream.

## Where π* comes from

The target is the pooled composition of the original training data, computed once through
IQ-TREE's estimator and held fixed across cleaning treatments (D07). The estimator (verified)
counts states, runs an 8-round fixed point that distributes ambiguity codes
(`convertCountToFreq`, `alignment/alignment.cpp:4886-4933`), then raises every entry below
`min_state_freq` to that value and takes the deficit from the largest entry (`convfreq`,
`alignment/alignment.cpp:5882-5909`). Missing sequences count as unknown states under `-p` but
not under `-S` (verified, `model/partitionmodel.cpp:132-133`), so the partition type is part of
the record, together with the gap and ambiguity treatment, any smoothing, and the vector at 17
significant digits. The real π* is computed only after the confirmatory plan is in the
repository, because computing it reads the training data; S0 to S4 use test targets.

A kept script in the oracle reproduces the estimator. It is cross-checked against IQ-TREE's own
"Mean state frequencies" line for a linked `+FO` run under the same partition type, which
IQ-TREE prints at 8 significant digits (verified, `model/partitionmodel.cpp:155-160`).

The model receives the target as `NQC+F{p1,...,p20}` in `-m` or `--model-joint`, or as
`NQC+F<NAME>` with the vector defined in an `--mdef` file (decision 006). It parses and validates
the vector itself, states any normalization, and rejects `+FO`, `+F`, `+FQ` or a missing `+F`
with a message that points to `NONREV`. With `--model-joint`, every partition receives the same
string (verified, `model/modelfactory.cpp:284-305`), and a check before linked optimization
requires every linked object to hold the same target bit for bit.

π* adds no degrees of freedom. Inference is conditional on a plug-in target, so adding 19 to a
BIC count when comparing fixed-target and free-target procedures is not justified (synthesis
section 8).

## Code plan

The file-level change map, the NQC class specification and the test specification are in
`docs/agent/CODE_PLAN.md`. This section summarizes them.

### Change map

New files (proposed): `model/fixedeqchart.{h,cpp}` (pure module, Eigen and the standard
library only), `model/modelnonrevfixedeq.{h,cpp}` (the NQC class), `unittest/` (standalone
googletest project), `test_scripts/fixedeq/` (oracle, drivers, environment file, fixture
generator) and `.github/workflows/fixedeq.yaml` (fork-only workflow).

Shared IQ-TREE files touched, each opt-in or behaviour-preserving (proposed): in S2,
`model/CMakeLists.txt` and one `createModel` branch in `model/modelmixture.cpp`; in S2 and S3,
parse-only options in `utils/tools.{h,cpp}`; in S3, derivative delegation, a linked-compatibility
check and a post-fit hook in `model/partitionmodel.{h,cpp}`; in S4, the export hook and citation
block in `main/phyloanalysis.cpp`. Deferred: `utils/optimization.{h,cpp}` (S5) and
`main/phylotesting.cpp` (S6).

### Decisions the code plan rests on

Programming decisions in `DECISIONS.md`: 002 model class, name and registration; 003
authoritative coordinates and lossless checkpoint; 005 derivative-step delegation; 006 target
input; 007 domain and step settings; 008 unsupported combinations; 009 surfacing failed line
searches; 010 export format; 011 reference destinations; 012 C++ test framework; 013 oracle
location and environment; 014 line endings; 015 Q built in the class's decomposition, which
calls the unchanged base (superseding 004); 016 root policy for S0 to S4; 017 the oracle's port
of IQ-TREE's optimizer. Mathematical decisions: D01 to D10.

### Test strategy

Six layers, specified in `CODE_PLAN.md` section 3: (1) a Python oracle in
`test_scripts/fixedeq/`, derived from the unedited `docs/agent/design/scripts/` and P-log
Appendix B, whose own tests reproduce the design documents' numbers; (2) C++ unit tests of the
pure module against oracle fixtures; (3) a differential test of IQ-TREE's likelihood at a fixed Q
on a fixed rooted tree against the oracle, validated first against the unmodified binary; (4) a
regression baseline of NONREV, `NONREV+F{...}`, NQ.pfam, GTR20, UNREST, a Lie-Markov model and
`--model-joint NONREV` under `-p` and `-S` (list in `CODE_PLAN.md` section 3.4), recorded twice
from the unmodified build before the first source edit, which must reproduce exactly at `-T 1`,
with any other run held to the spread seen between the two repeats; (5) invariant and nesting
checks (GTR20-seeded non-decrease, target immutability, parameter counts, thread count and
partition order, export and re-import); (6) Debug builds with AddressSanitizer and
UndefinedBehaviorSanitizer, run first on the unmodified tree so that existing findings are known.

GitHub Actions runs upstream's workflow on every push to any branch (verified,
`.github/workflows/ci.yaml`). It supplies cross-platform compile checks (clang 22 and gcc on
Linux, Apple clang, clang cross-compiling for MinGW) and upstream's turtle regression, but never
runs NQC, and its compilers differ from the local build, so no bitwise comparison happens there.
The fork-only workflow added in S1 runs the unit tests and the oracle's tests.

### Implementation slices

Each slice maps to a gate of synthesis section 11: G0 protocol and configuration, G1
mathematical core, G2 first compiled conditional model, G3 linked and lifecycle support, G4
competing implementations, G5 native and outer workflow, G6 scientific validation.
`CODE_PLAN.md` section 4 lists each slice's files and tests in full.

| Slice | Gate | Goal | Main files | Tests written first | Done when |
|---|---|---|---|---|---|
| S0 | G0 | Evidence before any source edit | `test_scripts/fixedeq/` only; no IQ-TREE source | oracle tests; regression driver | environment built from its file; design scripts rerun from copies and compared with their reported values; oracle likelihood matches the unmodified binary; baseline recorded twice with the binary's SHA-256; S0 runtime checks answered; unmodified sanitizer run recorded; G0 run manifest written |
| S1 | G1 | Pure mathematical core | `model/fixedeqchart.*`, `unittest/`, fixtures, fork workflow | unit tests at n = 3, 4 and 20 against oracle fixtures | all pass in Release and sanitizer builds; residual maxima tabled |
| S2 | G2 | First compiled model, one alignment | class files, `model/modelmixture.cpp`, `model/CMakeLists.txt`, domain option | differential at the seed; gradient step checks; target-overwrite and invalid-target rejections; re-evaluation of the accepted state; compiled optimizer against the Python port on a toy case | all pass; Level 1 and Level 2 fits on `example/aa_example.phy`; regression reproduces; sanitizers clean for new code |
| S3 | G2, G3 | Linked training and lifecycle | `model/partitionmodel.*`, class checkpoint | two partitions with different compositions; `-S` and `-p`; `-T 1` against `-T 4`; restart; GTR20-seeded nesting | all pass; the linked gradient uses the scaled step |
| S4 | G3 | Export, re-import and report | `main/phyloanalysis.cpp`, class export | export then `-m FILE` likelihood equality; metadata completeness | G3 met; first release candidate |
| S5 | G4 | Competing implementations | pure module (T3, positive ratios), chart selector, optional status recording | matched starts and domains; common diagnostics | comparison report |
| S6 | G5 | Workflow | `main/phylotesting.cpp`, outer driver | candidate-set and guard tests | G5 conditions |
| S7 | G6 | Science | none planned | the confirmatory plan's checks | G6 conditions |

## Risks and open questions

| # | Item | Resolved by |
|---|---|---|
| 1 | Final public model and option names | Peter, with the maintainers |
| 2 | A `+` or `*` inside a number, for example `1e+00`, splits a `+F{...}` model string (verified by reading, `model/modelfactory.cpp:259, 296`) | S0 run; documented; `+F<name>` avoids it |
| 3 | Whether existing flags give a fully fixed-nuisance Level 1 fit | S0 run |
| 4 | An output carrying 10 or more significant digits of the log-likelihood for the differential test; first candidate `--show-lh` (verified, `utils/tools.cpp:3428-3437`), then `-wsl`, then a debug print from the class | S0 |
| 5 | Cost: 361 likelihood evaluations per one-sided gradient; realistic training runs may take hours | S3 profiling; analytic gradients stay deferred |
| 6 | AddressSanitizer with link-time optimization and libomp inside the 6 GB WSL VM | S0 attempt; UndefinedBehaviorSanitizer-only fallback |
| 7 | Active `ASSERT`s abort a linked round that lowers the log-likelihood by more than 0.1 (verified, `model/partitionmodel.cpp:938`, `model/partitionmodelplen.cpp:135`) | incumbent protection and tests in S2 and S3 |
| 8 | The default coordinate domain may not contain LG rebuilt at a skewed π* | construction check and the domain option (decision 007) |
| 9 | Real training data: a small real set for S5, and the full set with a written confirmatory plan for S7 | Peter |
| 10 | The design scripts have not been rerun in this repository, and `e6_n20_big.py` cannot run | S0 reruns; mismatches reported |
| 11 | Scientific success criterion and manuscript details (placeholder above) | Peter |
| 12 | A continuity source the chart does not remove: `computeTransMatrixNonrev` switches from the eigen path to scaling-and-squaring when P's row sums deviate by more than 1e-4 (verified, `model/modelmarkov.cpp:489-500`), and the decomposition sets `nondiagonalizable` on a singular eigenvector matrix (verified, `model/modelmarkov.cpp:1330-1336`); shared with `NONREV`, frequency unmeasured | count "INFO: Switch to scaling-squaring" lines (printed under `-v`) in S2 and S3 fits |
| 13 | Resolved 2026-10-03: whether decision 004's `ModelMarkov` edit was still needed | decision 015, which supersedes 004; `ModelMarkov` is no longer edited |
| 14 | S5's recording-only edit to `dfpmin` and `lnsrch` (decision 009) touches the optimizer the maintainers advised leaving unchanged, and cannot be made by override (private, non-virtual, verified `utils/optimization.h:229-232`) | Peter, with the maintainers, before S5 |
| 15 | Resolved 2026-10-03 for S0 to S4: the root policy, on which a non-reversible likelihood depends | decision 016; the policy for scientific runs at G5 |
| 16 | How precisely a fitted `GTR20+F{π*}` incumbent reaches an NQC run: the report prints Q at 6 digits, decision 010 exports only NQC, and checkpoints hold 10 digits | S0 probe (g), then a decision entry before S3 |
| 17 | Thresholds marked provisional (1e-10 residuals, 1e-12 relative Q entries, the 1e-6 target sum) must be final, or declared reported rather than asserted, before their tests are written, because a threshold is not relaxed after a failure | Peter, before S1 |
| 18 | No fallback is recorded if S3 profiling shows training runs impractical while analytic gradients are out of scope | Peter, after S3 profiling |
| 19 | Baseline run 11 (`--model-joint NONREV` under `-p` at `-T 4`) does not reproduce: five repeats of the unmodified binary ended at log-likelihoods from -4975.0215 to -4974.5963, none equal to the `-T 1` run's -4974.5432, with the same topology and 442 of the model section's 447 numbers varying (largest spread 0.046). Its tree drawing varies too, so the spread rule of `CODE_PLAN.md` section 3.4 cannot judge it, and S3's `-T 1` against `-T 4` test meets the same behaviour in legacy code | Peter: a comparison rule for run 11 before S2, and the S3 thread test's definition before S3 |

## Current state

The design is settled (D01 to D10) and programming decisions 001 to 017 are recorded, 015
superseding 004; 016 sets the root policy for S0 to S4 and 017 the oracle's optimizer port.
The planning documents are committed on the branch. The IQ-TREE maintainers accepted the log-ratio jump-chain approach and stressed
that BFGS needs a continuous objective (reported by Peter, 2026-10-01; the mapping of their notes
to the code is in `CHANGELOG.md`). The unmodified branch builds in WSL2 and passes the smoke tests recorded in
`CHANGELOG.md`. No source code has changed. The conda environment `iqtree3-fixedeq` exists in WSL,
built from `test_scripts/fixedeq/environment.yml` with its lock in `environment.lock.txt`. The
regression baseline is recorded (2026-10-03) from the frozen unmodified binary
`~/iqtree3-baseline/iqtree3` by `test_scripts/fixedeq/regression/regress.py`, with its summary and
provenance in `test_scripts/fixedeq/regression/baseline/`: runs 1 to 10 reproduce exactly between
repeats, run 11 does not (risk 19). The oracle and fixtures do not exist yet.

## Next step

Rerun the three supplied design scripts from copies outside the repository, in the
`iqtree3-fixedeq` environment, and compare their outputs with the values the design documents
report (risk 10); then write the oracle and its tests (`CODE_PLAN.md` section 3.1). Before S2,
Peter sets the comparison rule for run 11 (risk 19).

## Checklists

The mechanical wiring of a substitution model class is described in `docs/agent/ARCHITECTURE.md`
section 14. The per-slice file and test lists are in `docs/agent/CODE_PLAN.md` section 4.
