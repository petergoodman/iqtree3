# Decisions

> **Status, 2026-09-23.** One programming decision is recorded (001). Ten mathematical and
> methodological design decisions, D01 to D10, are recorded from the design synthesis and are
> revisable (see that section's preamble). The items under "Pending programming candidates" were
> discussed on 2026-09-15 and have not been confirmed by Peter. An agent must not treat a pending
> item as a decision.

## What this document is for

This file is primarily the record of **programming decisions** made during development:
load-bearing choices about code, with the reasoning that produced them, so that later work does
not re-litigate them or silently contradict them. It also holds, in a separate section, the
mathematical and methodological decisions from the design stage that the code is built against.
The two kinds are kept apart because they change for different reasons: a programming decision
changes when the code or the codebase demands it, a mathematical decision changes when the
evidence about the model or the optimizer does.

It is **append only**. Entries are never edited or deleted. A decision that turns out to be wrong
is superseded by a new entry that references the old one by its number or ID.

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

Number programming entries sequentially from 001 and never reuse a number.

---

## Programming decisions

## 001. Develop on the fork's `nq-constrained-pi` branch, for an eventual upstream pull request

- **Date:** 2026-09-23
- **Status:** accepted
- **Decision:** All work happens on branch `nq-constrained-pi` of the fork
  `petergoodman/iqtree3`. Nothing is committed to `master` and nothing is pushed to
  `iqtree/iqtree3`. The work is intended to be offered upstream as a pull request from the fork
  after extensive testing and revision, so changes to shared IQ-TREE code are kept opt-in, with
  existing behaviour unchanged by default.
- **Why:** Peter's statement of intent, 2026-09-23. Keeping `master` an exact mirror of upstream
  keeps syncing a fast-forward, and keeping shared-code changes opt-in keeps an eventual pull
  request reviewable.
- **Alternatives rejected:** a permanent private fork with no merge intent, which would permit
  invasive edits to shared code but forfeit contribution upstream.
- **Affects:** every code change; the git procedure in `CLAUDE.md`; how far edits to shared files
  (`utils/optimization.cpp`, `model/partitionmodel.cpp`, `main/phyloanalysis.cpp`,
  `main/phylotesting.cpp`) may go.

---

## Mathematical and methodological decisions (D01 to D10)

These entries carry the IDs D01 to D10 of the design synthesis,
`docs/agent/design/project2_unified_discrepancy_synthesis.md` (its section 1), so that they can
be cross-referenced with it. They were adopted by Peter on 2026-09-23 as the design the code is
built against. **They are mathematical and methodological decisions, not programming decisions,
and they are subject to change.** A change follows the synthesis's change-control rule (its
section 13): name the affected ID, the evidence, whether the change affects the estimand or only
the implementation, the regression implications, and whether finished comparisons must be rerun.
It is recorded as a new entry that references the D-number, not by editing the D entry. How each
decision is realised in code is a programming decision, recorded as a numbered entry above.

The family throughout is the set of 20-state generators with positive off-diagonal rates, zero
row sums, π*Q = 0, and unit mean rate at π*: 360 free parameters, containing the 189-dimensional
reversible family at π*.

## D01. Parameterize the family by jump-chain log-ratios first

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** The first implemented chart is the jump chain in log-ratio coordinates: per row,
  18 log-ratios against a fixed reference destination, a stable row softmax giving the jump
  matrix K, one stationary solve νK = ν, and exit rates ν_i/π*_i (360 coordinates in all).
  Reference destinations are chosen once, from the row maxima of the target-rebuilt LG seed with
  ties broken in amino-acid order, and are shared by every linked object and every start.
- **Why:** Both log-ratios and positive weights cover the whole positive family one-to-one
  (synthesis section 3.1). Log coordinates make steps multiplicative in the weights. In the
  synthesis's rerun of four eight-state conditional fits, log coordinates scored higher in three,
  and the three positive-weight runs that stopped on `TOLX` coincided with failed line searches
  (synthesis sections 10.3 and 10.4). This is evidence for a starting choice, not a theorem.
- **Alternatives retained:** positive weights with one pinned weight per row (c = 10), kept as a
  maintained benchmark. Reopen if matched compiled tests show them more reliable or efficient at
  comparable solution quality.
- **Source:** synthesis section 3; P-log section 3; P-positive section 3.1.

## D02. Use scale-aware finite differences from the first fitted log-coordinate implementation

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** The first fitted log-coordinate implementation includes an opt-in derivative step
  h_k = η max(s_k, |x_k|), with initial s_k = 1 and η = 1e-4 (a reproducible starting setting,
  not a validated constant), reachable from both the single-model optimizer path and the linked
  `PartitionModel` path. Existing models keep the legacy step by default.
- **Why:** IQ-TREE's step is 1e-4 |x_k|, which collapses near zero for signed coordinates. In the
  synthesis's rerun sweep, the relative derivative discrepancy of one logit was about 9.6e-2 at
  z = 1e-6 and 6.4 at z = 1e-8 (synthesis section 10.2). In linked training the optimizer runs on
  the `PartitionModel` object, so an override in the model class alone is not reached
  (P-positive section 4.6).
- **Alternatives retained:** the legacy relative rule as a diagnostic control; analytic
  likelihood gradients later.
- **Source:** synthesis section 4.1; P-log section 6.2.

## D03. Declare numerical domains explicitly and report their effect

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** Coordinate bounds are explicit and recorded for every run. The historical
  matched-benchmark domain is z in [log 1e-5, log 10], about [-11.513, 2.303], which corresponds
  to positive weights in [1e-4, 100] with c = 10; it is a benchmark domain, not a chosen
  production range. Every declared start must be representable with margin, bound activity is
  reported, final scientific runs include a domain-expansion sensitivity check, and an imported
  matrix is never clipped to fit a box.
- **Why:** A finite box restricts the positive family, and boxes in different charts are
  different sets of matrices (synthesis section 4.2; P-positive Theorem 6, remark 2).
- **Alternatives retained:** expanding or revising the limits when bounds are active or
  sensitivity checks move the solution; a boundary-capable method if exact zero rates are ever
  required (D10).
- **Source:** synthesis section 4.2.

## D04. Implement T3 balancing as a second compiled backend

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** After the log-jump core is validated, the T3 balancing chart is implemented in
  IQ-TREE on the same likelihood engine and derivative interface, as a serious comparison
  method. Its gauge is a_{19,20} = 0 and h_{i,20} = 0 for i < 20 (189 strength and 171 tilt
  coordinates), with the 19 balancing potentials found by damped Newton. A Python version is a
  correctness reference and does not by itself meet the comparison goal.
- **Why:** T3 and the jump chart parameterize the same family, so attained differences reflect
  optimization, bounds, or implementation, which is what the comparison is meant to measure
  (synthesis section 6.1). Its production exposure depends on benchmarks and maintainer
  preference.
- **Alternatives retained:** P-positive's star-at-state-1 gauge, which describes the same family
  and is imported only by explicit conversion.
- **Source:** synthesis section 6.1; P-log section 4.

## D05. Build the first seed and the guide candidates at the target

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** The default first fitting start is LG rebuilt at the target, q_ij = R_ij π*_j,
  normalized to unit mean rate. The fitted reversible model at the target (`GTR20+F{π*}`) is an
  additional start and the nesting incumbent. The primary guide-candidate protocol uses LG, WAG
  and JTT rebuilt at π*, with empirical and optimized frequency variants suppressed.
- **Why:** One composition convention holds across every stage. Guide candidates only select
  per-alignment models and trees, so ordinary candidates would not invalidate the constrained
  fit (synthesis section 6.2; P-log section 7.3). Row-normalizing LG at LG's own frequencies
  gives a valid but different start, whose exchangeabilities are not LG's.
- **Alternatives retained:** ordinary empirical guide candidates, and flux-transported or
  T3-transported starts, as sensitivity and multistart arms.
- **Source:** synthesis sections 6.2 and 6.3; P-positive section 3.1.5.

## D06. Change only the map from coordinates to generator, and protect its invariant

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** IQ-TREE's likelihood engine, transition matrices, pruning, and branch, rate and
  tree algorithms are reused. The substantive change is the map from 360 optimizer variables to
  a valid Q in the shared-matrix update (nQMaker step 3a), with the root distribution equal to
  π*. The fixed target, chart metadata and reference destinations must survive every part of the
  model's lifecycle (construction, linking, checkpoint, export), and linked objects must agree on
  them before optimization. Outer candidate reselection is orchestrated explicitly rather than
  assumed to work through one ModelFinder call.
- **Why:** Synthesis section 7. The inspected ModelFinder code rejects candidate sets that mix
  recognized reversible and nonreversible models (`mixRevNonrev`, `main/phylotesting.cpp`).
- **Alternatives retained:** a thin external driver for the outer loop first, native integration
  later.
- **Source:** synthesis section 7; P-positive section 5; P-log sections 7.3 and 8.

## D07. Take one explicit, immutable target and never floor it silently

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** The fitter consumes one explicit, strictly positive, normalized π* in IQ-TREE's
  state order and never re-estimates it. The default provenance is the pooled composition of the
  original training data through IQ-TREE's estimator, computed once and held fixed across
  cleaning treatments, with the estimator, ambiguity and gap treatment, any smoothing, and the
  exact numbers recorded. An estimator floor is a declared procedural choice; the fitter does not
  floor supplied input, and it rejects values its numerical backend cannot support.
- **Why:** Synthesis section 8. `ModelMarkov::targetFunk` rejects positive frequencies below
  `min_state_freq` (default 1e-4, configurable with `--min-freq`), and the decomposition removes
  states whose frequency is below `ZERO_FREQ` = 1e-10, so a small target entry must be checked
  against both rather than silently changed.
- **Alternatives retained:** other declared estimators and supported small positive targets.
- **Source:** synthesis section 8; P-positive section 7; P-log section 1.4.

## D08. Plan reversible, unrestricted and transport comparison arms

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** The scientific comparison includes original-data reversible and unrestricted
  fits, cleaned-data reversible and unrestricted fits, the cleaned reversible fit followed by
  frequency replacement, the reversible refit at π*, the general refit at π* (the Project 2
  estimate), and the cleaned unrestricted fit transported to π* by fixed flux and by T3 without
  refitting. Frozen matrices are evaluated on the same held-out data. The transport arms do not
  block the first implementation.
- **Why:** Since π* e^{tQ} = π*, every tip marginal is fixed, and a wrong target can produce
  apparent circulation; simulations with reversible and nonreversible generators at correct and
  perturbed targets are needed before circulation is read as biology (synthesis section 9.1).
- **Alternatives retained:** none dropped; the arms are the comparison plan.
- **Source:** synthesis section 9.1; P-positive sections 8 and 9; P-log section 10.3.

## D09. Adopt the corrected mathematical qualifications

- **Date:** 2026-09-23
- **Status:** adopted; mathematical, revisable
- **Decision:** These statements override conflicting prose in the older documents: (1) with raw
  off-diagonal rates in [l, u], the gauge-fixed T3 triangle tilt satisfies
  |h'| <= (3/2) log(u/l), which is 20.72 for [1e-4, 100], so a bound of 12 does not contain that
  box; (2) q_ij <= 1/(2 π*_i); (3) global suprema over nested domains are ordered, but
  independently attained local fits need not be, so the smaller model's fit and nuisances are
  imported as an incumbent; (4) an interior true parameter does not imply an interior sample
  maximum; (5) the root-split singularity under reversibility means the 171 and 19 parameter
  differences do not automatically give chi-square reference distributions; (6) replacing
  frequencies in asymmetric factors can still give a valid generator, and what is lost is target
  stationarity; fixing root frequencies alone is not the constraint; (7) no claim of a global
  optimum is made, and the BFGS inverse Hessian is not a covariance estimate.
- **Why:** Synthesis section 9.2, with derivations in P-positive section 2 and P-log sections 4.5,
  5 and 10.2.
- **Alternatives retained:** not applicable.
- **Source:** synthesis section 9.2.

## D10. Keep direct balanced flux as an independent reference and defer EM

- **Date:** 2026-09-23
- **Status:** adopted; methodological, revisable
- **Decision:** Direct balanced-flux fitting, f = f0 + Zy with a feasible constrained solver, is
  kept as a small-problem reference and a possible boundary-capable backend. Constrained
  substitution-history EM is deferred; if implemented, the normalized M step (19 balance and one
  normalization multiplier) is preferred. Penalty and projection methods do not substitute for
  the constrained maximum-likelihood estimate.
- **Why:** Synthesis section 9.3.
- **Alternatives retained:** EM, reopened for persistent direct-optimization failures, a need for
  exact boundary fits, or likelihood-gradient infrastructure.
- **Source:** synthesis section 9.3; P-positive section 3.3; P-log sections 5.3 and 9.

---

## Pending programming candidates

These need Peter's confirmation before they become numbered entries.

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

> Note, 2026-09-23: both design plans also propose a new name (`NQC` provisionally), and D09(5)
> qualifies the likelihood-ratio statement above.

### Candidate: support single-alignment and multi-partition paths through one subclass

Supporting evidence: `PartitionModel::optimizeLinkedModel` is parameterization-agnostic. It
calls `setVariables`, `setBounds`, and `targetFunk` on the model and never inspects what the
parameters mean, so a correct `ModelMarkov` subclass works under both `-m <name>` and
`--model-joint <name>` with no additional code. Two conditions attach: the parameterization must
be a pure function of the variable vector plus fixed per-model constants, with no hidden
per-partition state; and one shared Q across partitions implies one shared π, so per-partition π
with a shared constrained Q would be a different feature.

> Note, 2026-09-23: "no additional code" is contradicted by D02 and D06. The linked path computes
> gradients through `Optimization::derivativeFunk` on the `PartitionModel` object, which no model
> class can override; the `PartitionModel` constructor calls `adaptStateFrequency` for linked
> models with `FREQ_ESTIMATE` or `FREQ_EMPIRICAL` (`model/partitionmodel.cpp:116-165`); and
> linked objects need an explicit compatibility check.

### Candidate: follow the Lie-Markov reduced-parameterization pattern

Supporting evidence: `ModelLieMarkov::setBasis` already implements this exact feature for
4-state DNA, including reducing the free parameter count by the frequency degrees of freedom and
warning when a requested π is unreachable. See `docs/agent/ARCHITECTURE.md` section 9 for what
transfers to 20 states and what does not.

> Note, 2026-09-23: superseded in substance by D01, which fixes the chart. The "what transfers"
> passage of `ARCHITECTURE.md` section 9 has been removed as design advice. `ModelLieMarkov`
> remains relevant only as a precedent for storing coordinates separately from `rates[]` and
> mapping them in `getVariables` (synthesis section 7.1).

### Candidate: baseline and branching strategy

Which upstream commit the work is based on, and how the fork is structured. Not yet decided.
Record the chosen baseline commit here once fixed, because the anchors in
`docs/agent/ARCHITECTURE.md` and `docs/agent/FILE_INDEX.md` depend on it.

> Note, 2026-09-23: the branching half is now entry 001. The baseline is currently `63c330d9`
> (tag `v3.1.4`), the commit the branch and the documentation anchors are based on; a later sync
> would move it by the procedure in `CLAUDE.md`.
