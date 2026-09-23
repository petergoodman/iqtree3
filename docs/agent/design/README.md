# Design documents for π-constrained nonreversible estimation

This folder holds the mathematical and methodological design for estimating a nonreversible
amino-acid rate matrix at a fixed target stationary distribution π* inside IQ-TREE 3. The
documents were written before any code was changed and before IQ-TREE had been built on the
project machine, so every statement about IQ-TREE in them comes from reading source at commit
`63c330d9` (tag `v3.1.4`), not from running it.

These files are copies, byte-identical to the operator's working files at the time of import
(2026-09-23). Do not edit them. A revised design is added as a new file and recorded in
`CHANGELOG.md`.

## Files

| File | Alias used in the synthesis | Date | Author as stated in the file | SHA-256 |
|---|---|---|---|---|
| `project2_unified_discrepancy_synthesis.md` | (the synthesis) | 2026-09-23, version 1.0 | not stated | `7b37d0255bde46f73ea223de28cecdb64e45466ed417b5f19a29b57e8a326d41` |
| `constrained_nq_plan_agent.md` | P-positive | 2026-09-22 | Claude, model `claude-fable-5-1` | `5125f4626f86a73cdbc0965375610d8bf03cc679d16899224116c123e712de43` |
| `project2_iqtree3_technical_plan.md` | P-log | 2026-09-22 | not stated | `953438a82db410ca8c474d433ddf89f5e4208d3d743d5318a824c40c4349ad5b` |

The synthesis cites the two plans as `constrained_nq_plan_agent(1).md` and
`project2_iqtree3_technical_plan.md` and gives their SHA-256 fingerprints in its section 14. Both
copies here match those fingerprints.

## Precedence

The synthesis sets the order (its section 0): explicit later decisions by Peter come first, then
the synthesis for the conflicts it resolves, then the two plans for detail compatible with it.
The synthesis's ten decisions, D01 to D10, are recorded in `docs/agent/DECISIONS.md` as
mathematical and methodological decisions. They are revisable, and a revision is recorded there
as a new entry.

## What these documents are not

They are not a code plan. The code plan, once written, lives in `docs/agent/PLAN.md`.
Programming decisions made during development are numbered entries in
`docs/agent/DECISIONS.md`. The descriptions of IQ-TREE's current code live in
`docs/agent/ARCHITECTURE.md`, `docs/agent/AA_MODEL_INFERENCE.md`, and
`docs/agent/FILE_INDEX.md`, and where those differ from a statement here about current code, the
source at the baseline commit decides.

## Scripts cited by these documents

`scripts/` holds the design stage's numerical check scripts, copied unmodified on 2026-09-23.
They are reference material for agents: records of how the design's numbers were produced, and a
starting point for the independent test reference that development code will be checked
against. They are not development code, are not built or run by IQ-TREE, and must not be edited
in place; working code derived from them goes elsewhere, as the code plan specifies.

| Script | Backs | SHA-256 |
|---|---|---|
| `scripts/verify_constrained_nq.py` | P-positive section 12, checks 1 to 12 | `fef4315dd728712197a708893816cd5ef9c25b558df0736336095b0166d56b4b` |
| `scripts/check_builtin_matrices.py` | P-positive section 12, LG and `NQ.PFAM` figures | `71c120b94d1e68019faa927f3c32c3243ed33871de1073e8fa9f3573c5f4233c` |
| `scripts/chart_fd_bfgs_compare.py` | synthesis section 10 (E1, E2) and Part 2 slide deck Backup B | `8f2a1d095b24efcef62cff881cea213a924d4af7f14279517ee4d90aa9d0591b` |

The fingerprint of `chart_fd_bfgs_compare.py` matches the one in the synthesis's section 14. The
other two have no published fingerprint to compare against.

Conventions differ between scripts. `verify_constrained_nq.py` implements the positive-ratio
jump chart (C2) and the T3 chart with its potential gauge at state 1, which D04 does not adopt;
`chart_fd_bfgs_compare.py` implements both the positive-ratio and log-ratio jump charts and a
Python port of IQ-TREE's `dfpmin`, `lnsrch` and `derivativeFunk`. The check code printed in P-log
Appendix B implements the log-ratio jump chart and the T3 gauge that D04 adopts; it exists only
inside `project2_iqtree3_technical_plan.md`.

Running `verify_constrained_nq.py` or `check_builtin_matrices.py` writes a `*_results.json` file
next to the script, so run a copy placed outside the repository rather than the file here.
`check_builtin_matrices.py` takes the path to `model/modelprotein.cpp` as its first argument.
None of the scripts has been rerun in this repository, so their numbers remain as reported in
the documents until they are.

Not included: `e6_n20_big.py` (synthesis section 10.1), which imports a module `e4_nonrev` that
was never supplied and so cannot run, and the experiment drivers for E3 to E5, which were not
supplied.
