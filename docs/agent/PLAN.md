# Plan: π-constrained non-reversible amino-acid models

> **Status, 2026-09-23: goal written, everything else placeholder.** Peter owns this document.
> The goal section below was written from the design documents in `docs/agent/design/` and the
> decisions D01 to D10 in `docs/agent/DECISIONS.md`. Every other section is still unwritten; an
> agent should treat those questions as open and ask rather than infer.

## What this document is for

This is the living core document for the project, rewritten as the project moves rather than
appended to. It answers four questions for any agent picking up work:

1. what we are building and why,
2. what constraints and requirements the work must satisfy,
3. where the project currently stands,
4. what the next step is.

It is read at the start of every task. Keep it current, keep it short enough to read in full,
and put history in `CHANGELOG.md` rather than here.

## How it relates to the other documents

| Document | Lifecycle | Read when |
|---|---|---|
| `docs/agent/PLAN.md` (this file) | overwritten | start of every task |
| `docs/agent/design/` | fixed copies; a revision is a new file | before planning, or before changing the method |
| `docs/agent/DECISIONS.md` | append only, never edited | before changing a settled choice |
| `docs/agent/ARCHITECTURE.md` | corrected in place as findings land | before writing code |
| `docs/agent/AA_MODEL_INFERENCE.md` | corrected in place | when the method of an existing model matters |
| `docs/agent/FILE_INDEX.md` | corrected in place | when locating code |
| `CHANGELOG.md` | append only | to resume, or to check whether something was already tried |

Plans live in this file and decisions in `DECISIONS.md` (programming decisions as numbered
entries, the design's mathematical decisions as D01 to D10). `ARCHITECTURE.md`,
`AA_MODEL_INFERENCE.md`, and `FILE_INDEX.md` describe IQ-TREE's code as it currently is and
prescribe nothing.

## Sections to fill in

### Goal

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

### Scientific motivation and success criteria

Why the constrained estimate is the right object, what result would count as the feature
working, and what comparison or validation demonstrates it. Include the intended manuscript or
question if there is one.

### Scope

What is in and out. Two scope questions have provisional answers recorded in `DECISIONS.md` and
should be restated here once settled: which invocation paths are supported (single alignment,
multi-partition, or both), and whether the feature gets a new model name or modifies `NONREV`.

### Constraints and requirements

Hard requirements the implementation must satisfy. Candidates to consider, not yet decided:

- the constraint holds exactly at every optimizer step rather than approximately or at
  convergence only,
- off-diagonal rates remain non-negative throughout,
- existing `-m NONREV` and `-m NQ.*` behaviour is unchanged,
- parameter counts reported for AIC and BIC are correct under the constraint,
- the feature checkpoints and resumes correctly,
- the feature works under the same options as unconstrained nQ estimation.

### Where π comes from

Decide and record whether the target π is supplied as a literal vector, taken as the empirical
composition of the alignment, or read from a named profile, and how that choice affects the
degrees of freedom charged in `getNDimFreq()`. For the multi-partition case, state which
alignment's composition "empirical" refers to.

### Current state

Overwritten each session. What exists, what builds, what is verified, what is broken.

### Next step

One concrete next action, not a list.

### Checklists

Implementation and verification checklists, added once scope is settled. A mechanical checklist
for adding a substitution model already exists at `docs/agent/ARCHITECTURE.md` section 14 and
should be referenced rather than copied.

### Open questions

Questions blocking progress, each with who or what resolves it.
