# Plan: π-constrained non-reversible amino-acid models

> **Status: placeholder.** Peter is writing this document. Nothing below the next heading has
> been filled in yet. An agent reading this should treat the project goal as not yet formally
> specified, and should ask rather than infer.

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
| `docs/agent/ARCHITECTURE.md` | corrected in place as findings land | before writing code |
| `docs/agent/FILE_INDEX.md` | corrected in place | when locating code |
| `docs/agent/DECISIONS.md` | append only, never edited | before changing a settled choice |
| `CHANGELOG.md` | append only | to resume, or to check whether something was already tried |

## Sections to fill in

### Goal

One paragraph stating precisely what the feature does, in terms a reviewer would accept. The
working statement so far, to be replaced: given a target stationary frequency vector π, estimate
by maximum likelihood the non-reversible amino-acid rate matrix Q that has π as its stationary
distribution, under otherwise the same conditions as existing nQ estimation.

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
