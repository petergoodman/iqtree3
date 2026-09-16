# Decisions

> **Status: placeholder, format defined, no decisions ratified.** The entries listed under
> "Pending" below were discussed in session on 2026-09-15 but have not been confirmed by Peter.
> They are recorded here so the reasoning is not lost, not because they are settled. An agent
> must not treat a pending item as a decision.

## What this document is for

Load-bearing choices, with the reasoning that produced them, so that later work does not
re-litigate them or silently contradict them. It is **append only**. Entries are never edited
or deleted. A decision that turns out to be wrong is superseded by a new entry that references
the old one by number.

It is deliberately distinct from `CHANGELOG.md`. The changelog is chronological and read
backwards from the end to answer "where did we leave off and what has been tried". This file is
a flat list read at the start of a task to answer "what am I not allowed to quietly change".

If maintaining both proves to be too much friction, drop this file and write decisions into
`CHANGELOG.md` with a `DECISION:` prefix so they stay greppable. Do not keep two homes for the
same content.

## Entry format

```markdown
## NNN. Short imperative title

- **Date:**
- **Status:** accepted | superseded by NNN | reversed
- **Decision:** one or two sentences, stated as a rule.
- **Why:** the reasoning, including the evidence that supports it.
- **Alternatives rejected:** what else was considered and the concrete reason it lost.
- **Affects:** files, commands, or behaviour that depend on this.
```

Number entries sequentially and never reuse a number.

---

## Pending, not yet decided

These need Peter's confirmation before they become entries 001 onward.

### Candidate: give the constrained model a new name rather than modifying `NONREV`

Supporting evidence gathered so far: `-m NONREV+F{...}` is already a live code path that pins
`state_freq` to the supplied vector and skips the stationarity solve in
`decomposeRateMatrixNonrev`, leaving Q unconstrained, so overloading that string would change
the meaning of commands that have already been run. The two models are nested at 360 and 379
free parameters, so separate names permit a likelihood-ratio test between them. Parameter
counts, bounds, and checkpoint payloads differ, so one name would mean branching inside every
one of those. `GTR20` and `NONREV` are already separate name branches in the same class, which
is direct precedent.

Open sub-question: whether to implement as a new `ModelMarkov` subclass registered through
`ModelMarkov::getModelByName`, following `ModelUnrest`, or as another name branch inside
`ModelProtein::init`. A subclass keeps the constraint machinery generic over state count and out
of a file that is mostly data tables.

### Candidate: support single-alignment and multi-partition paths through one subclass

Supporting evidence: `PartitionModel::optimizeLinkedModel` is parameterization-agnostic. It
calls `setVariables`, `setBounds`, and `targetFunk` on the model and never inspects what the
parameters mean, so a correct `ModelMarkov` subclass works under both `-m <name>` and
`--model-joint <name>` with no additional code. Two conditions attach: the parameterization must
be a pure function of the variable vector plus fixed per-model constants, with no hidden
per-partition state; and one shared Q across partitions implies one shared π, so per-partition π
with a shared constrained Q would be a different feature.

### Candidate: follow the Lie-Markov reduced-parameterization pattern

Supporting evidence: `ModelLieMarkov::setBasis` already implements this exact feature for
4-state DNA, including reducing the free parameter count by the frequency degrees of freedom and
warning when a requested π is unreachable. See `docs/agent/ARCHITECTURE.md` section 9 for what
transfers to 20 states and what does not.

### Candidate: baseline and branching strategy

Which upstream commit the work is based on, and how the fork is structured. Not yet decided.
Record the chosen baseline commit here once fixed, because the anchors in
`docs/agent/ARCHITECTURE.md` and `docs/agent/FILE_INDEX.md` depend on it.
