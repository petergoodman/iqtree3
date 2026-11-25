# Upstream findings

Evidence behind the questions for the maintainers in `docs/agent/PLAN.md` and behind the
upstream-behaviour risks 23, 24 and 30: what was run or read, on unmodified binaries and, where a
question compares them, on NQC builds. It describes and never prescribes. Entries are keyed by
PLAN.md's numbers, which stay the citation; questions 1 to 3 and 8 carry no evidence and live
only in PLAN.md. Not read at session start; read when relaying a question or rechecking an
upstream behaviour. Source anchors are in `CODE_PLAN.md` section 5, named here by symbol.

## Q4

Upstream finding: `ModelMarkov`'s constructor calls `setReversible`, which reads
`is_reversible` before anything sets it (UndefinedBehaviorSanitizer: "load of value 190, which is
not a valid value for type 'bool'"). It appears once in every IQ-TREE sanitizer run, the
unmodified binary's included; its effect on results was not investigated.

Anchors: the `ModelMarkov` constructor's `setReversible` call and the read of unset
`is_reversible`.

## Q5

Upstream observation, S0 probe (g): after a tree search the checkpoint keeps the model from
before the final model optimization, and a restart from it prints "Final model parameters
restored" (the restart itself was not tested). NQC avoids it by decision 025.

Anchors: the post-search save that writes the search state, not the model (`runTreeReconstruction`).

## Q6

Upstream finding, S3b's sanitizer run: at `-v` and above, `runPhyloAnalysis` ends every
partitioned run whose partition type is not edge-unlinked (`-Q`) by casting the tree to
`PhyloSuperTreePlen` and printing its NNI counts. An `-S` tree is not one, so `printNNIcasesNUM`
reads past the object: AddressSanitizer reports a heap-buffer-overflow and stops IQ-TREE, and a
Release build prints meaningless counts. Reproduced on the unmodified sanitizer binary with
turtle, `-S`, `LG+G4{0.5}`, `-te -blfix -v`.

## Q7

Upstream observation, S0 probe (c2), unmodified binary: NONREV's turtle `-S` Level 1 fit
(`--model-joint NONREV -te -blfix -m LG+G4{0.5}`) ends at the default 10-round cap
(`--loop-model`) still gaining 3.51 per round, at −4815.29, while NQC, whose unbounded family
lies inside NONREV's, reaches −4806.09 on the same setup (risk 24). With the round rules relaxed
(`-me 1e-6`, 300 rounds; S3c's first part, S3b's binary, NONREV's path unchanged), NONREV's
turtle `-S` and `-q` fits stopped at rounds 41 and 47 on a round that gained exactly nothing,
after rounds gaining 2.9e-3 and 5.1e-3; on `aa_example` it ran to the cap still gaining. In
S3c's profiling (binary `81729f4b`, NONREV's path unchanged) its `-S` fits ended at the 10-round
cap in all six turtle and 10-set configurations. Cause not investigated.

Anchors: `dfpmin`, `lnsrch` and their failure and stopping paths; the round counts
`model_opt_steps` (`--loop-model`) and `num_param_iterations` (`-nparam`).

## Q9

Not verified: with `-mstore`, `ModelFactory` serves transition matrices cached by branch length
from the end of each model optimization, so a linked fit that follows a partition's own update
may evaluate stale matrices, for NONREV as for any linked model. NQC refuses `-mstore`
(decision 008).

Anchors: `-mstore`; transition-matrix storing stopped and started around each model
optimization; the cache keyed by branch length.

## Q10

Upstream finding, S4's sanitizer run: `Alignment::createBootstrapAlignment`, which
`runStandardBootstrap` calls for each `-b` replicate, copies a `pomo_sampling_method` the source
alignment never set (UndefinedBehaviorSanitizer: "load of value 3200171710, which is not a valid
value for type 'SamplingType'", reported once per run). Reproduced on the unmodified sanitizer
binary with `-m LG -b 2`; its effect on results was not investigated.

Anchors: `createBootstrapAlignment`; `runStandardBootstrap`.

## Q11

Upstream finding, 2026-10-28, unmodified binary: `--runs` with `-S` crashes. On the
two-partition data with `NONREV -te -blfix --runs 2` both runs finish and the crash comes as the
best run is restored; a backtrace of the same crash with NQC shows
`PhyloSuperTreeUnlinked::printResultTree`, called by `runMultipleTreeReconstruction`. On turtle
with `LG+G4 --runs 2` and a tree search it comes in run 1's parsimony start tree
(`PhyloTree::computeParsimonyTree`). By reading, not verified: `runMultipleTreeReconstruction`
builds each run's tree as a `PhyloSuperTreePlen` for every partition type but `-p`, `-S`
included.

## Q12

Upstream observation, 2026-10-29 audit: at more than one thread the non-reversible tree
log-likelihood is not bitwise reproducible from one call to the next, for NONREV as for NQC.
`computeNonrevLikelihoodBranchSIMD` adds its packets' sums under `schedule(dynamic,1)` with an
OpenMP `reduction`, while the reversible kernel stores each packet's sum and adds them in packet
order. The differences are rounding (1.3e-16 relative at `-T 4`, decision 038), but a bitwise
recheck in NQC's export aborted 3 of 8 runs at `-T 4`; NQC now compares within decision 038's
line. Would the maintainers want the non-reversible kernel to add its packets in order?

Anchors: the non-reversible and reversible packet sums; packets per thread; the `-p` and `-q`
thread cap.

## Q13

Upstream observation, 2026-10-29, unmodified binary: a rate model written inside
`--model-joint` is not applied. On turtle with `-te -blfix`, `--model-joint "NONREV+G4{0.5}"`
fits `NONREV+FO` with no rate model under `-S`, `-Q` and `-p`, exit 0, while a bare
`"NONREV+G4"` stops `-p` with "File not found NONREV+G4"; the rate model takes effect only
through `-m`. NQC now refuses one (decision 008). Should `--model-joint` reject a rate model, or
apply it?

Anchors: `model_joint` in `ModelFactory` and `--model-joint`'s parsing.

## Q14

Upstream finding, 2026-10-29, unmodified binary: `--lie-markov` with a non-reversible model
(NONREV on aa_example, UNREST on `example.phy`) stops at `ASSERT(0 && "this line should not be
reached")` in `ModelMarkov::computeTransMatrixNonrev`; LG runs. NQC now refuses it (decision
008). Should it be refused for non-reversible models, or fall back to scaling and squaring?

Anchors: scaling and squaring for non-reversible models; the `ASSERT(0)` for other techniques.

## Q15

Upstream finding, 2026-10-29, unmodified binary: a one-alignment `-te` run that stops after its
model fit and before its end-of-run dump (killed at `reportPhyloAnalysis`) resumes with "Model
parameters restored", reports the fitted model on the input tree's branch lengths and exits 0:
aa_example on baseline run 2's tree, `LG+G4` −7330.3269 against −7302.3888 for the run that
finished, `NONREV+G4` −6820.0884 against −6811.2908; NQC the same. By reading,
`computeInitialTree` saves the `-te` tree into the checkpoint before `runTreeReconstruction`
restores the tree from it, and the fitted tree reaches the candidate set only after the fit's
dump; resuming with the fitted tree as `-te` gives the finished run's score. Linked runs not
tested. Script and outputs: `~/iqtree3-audit/followup-20261030T054825Z/`. Should the restore keep
the fitted tree?

Anchors: `computeInitialTree` and `initializeParams`; the candidate set after the fit's dump;
the end-of-run forced dump.

## Q16

Upstream finding, 2026-10-30, unmodified binary: `--root-test` re-optimizes every branch length
after the model fit, `-blfix` included, and the report then mixes two states. On aa_example with
baseline run 2's tree, `NONREV+G4{0.5} -te -blfix --root-test` ends with tree length 3.4968
against 3.0282 without the option, its root test reaching −6819.462 at the original root, while
the `.iqtree` prints the pre-test −6825.6805 beside the re-optimized tree. With a reversible model
on an unrooted tree (`LG+G4{0.5}`) `testRootPosition` writes no `.roottest.trees` and the run
stops with "File not found" (exit 2). NQC refuses `--root-test` (decision 008), since its export
would describe the re-optimized state. Script and outputs:
`~/iqtree3-audit/rootfix-20261030T174542Z/step0a/`. Should the root test honour `-blfix` and
leave the tree as it found it?

Anchors: `root_test`; `testRootPosition`, its only caller and its branch-length re-optimization.

## Q17

Upstream finding, 2026-11-17, unmodified binary (`~/iqtree3-baseline/iqtree3`, v3.1.4): `GTR20` on a
rooted `-te` tree fails. On turtle (`turtle_aa.fasta`, the `-m LG -n 0` tree rooted at the midpoint of
Anolis's pendant edge), `-p turtle_aa.nex --model-joint GTR20 -te rooted.nwk` segfaults at "Estimate model
parameters" (exit 139, printed by the script, not saved; the log shows the segmentation fault), and the same
tree unrooted runs to completion; `NONREV` on both trees passes that point and was still optimizing when the
script's 20 s timeout stopped it. A first turtle attempt (`gtr20-turtle-20261117T194343/`) used a tree whose
taxa did not match the alignment and ran nothing. A single model, `-m GTR20
-te rooted.nwk`, stops at `ASSERT(tree_lh-new_tree_lh < max_delta_lh)` in `PhyloTree::optimizeAllBranches`
(new_tree_lh −12638.85 against tree_lh −9639.93), Anolis's root edge grown from 0.027 to 2.36. On SIM-1's
52-taxon bird data (30 genes) the segfault reproduces under `-q` and `-p`, with and without `+F{π*}` and
`-blfix`, at `-T 1` and `-T 6`, on the NQC gate binary as on the unmodified one; on the rooted tree
`LG+F{π*}` (no free parameters) runs to completion, and `NONREV+F{π*}`, like `GTR20` on the tree unrooted,
passes the crash point and was still optimizing at the script's 30 s timeout (`GTR20` on the bird tree
unrooted later ran to completion in SIM-1's pipeline). The scripts call their helpers by a path in another
session's scratchpad; copies sit beside them (`unroot_tree.py`, `turtle_rooted.py`). Scripts and outputs:
`~/iqtree3-scratch/gtr20-turtle-20261117T194402/`, `gtr20-narrow-20261117T194153/`,
`gtr20-crash-20261117T194123/`. Cause not investigated. SIM-1's pipeline gives `GTR20` the tree unrooted,
on which a reversible likelihood is the same. Should a reversible model unroot a rooted `-te` tree, or
refuse it?

Anchors: `PhyloTree::optimizeAllBranches`'s assertion; the segfault's site not located.

## Q19 and Risk 30

Read in source at `63c330d9`: `ModelProtein::init`'s GTR20 branch seeds `rates[]` from `--init-model` (a
model name, or a file through `readParameters`) or from `Params::gtr20_model` (`--init-exchange`, default
POISSON), calls `rescaleRates`, which scales the largest rate to 10 (`AA_SCALE`), and frees 189 rates, so
`ModelMarkov::setVariables` leaves `rates[189]`, Tyr–Val, at the seed's value. `ModelMarkov::setBounds`
holds every free rate in [`MIN_RATE`, `MAX_RATE`] = [1e-4, 100], with `bound_check` false, so without
restarts. From POISSON every rate starts at 10, so Tyr–Val is fixed at 10 and every other exchangeability
lies in [1e-5, 10] × Tyr–Val. `ModelMarkov::isUnstableParameters` flags a rate above 0.99 × `MAX_RATE`;
after the model fit `IQTree::optimizeModelParameters` prints "Estimated model parameters are at boundary
that can cause numerical instability!" when it is true for any partition (`PartitionModel::
isUnstableParameters`) and leaves the fit as it is. Whether SIM-2's fits printed it was not checked:
their logs were not collected. NONREV's branch, by reading only, converts its rescaled LG seed with LG's
frequencies (`setReversible(false)`, q_ij = R_ij π_j) and keeps `rates[379]` (Val→Tyr), about 0.0080,
so its window is 0.0125 to 12,500 times that rate: NQ.bird's rates span 0.33 to 2,840 times its own
Val→Tyr, inside it, and LG's 0.0053 to 86.5, the lowest below the floor (computed from the built-in
matrices, not from a NONREV run). Under `--model-joint` each partition builds its
own `ModelFactory` (`PartitionModel`'s constructor), so a seed file reaches every partition (traced by a
read-only review on 2026-11-24, which also found that `readParameters` leaves the frequency type alone).
IQ-TREE's BFGS tests its step and gradient against max(|x|, 1) (`dfpmin`, `lnsrch`), so where the fixed
rate sits also changes the free rates' scale for the optimizer.

Seen in SIM-2's round 1 (gate binary `25d78d6e`, whose GTR20 code is upstream's; `-q --model-joint`,
unrooted `-te`, `-blfix`, `-keep-ident`, 20-round cap; cell D, 100 genes, 52 taxa):
- Every one of cell D's 160 fits (R0 and NQ, every masking strength): checkpoint Tyr–Val 10 exactly and
  exactly one rate above 99 (at 100 in the 40 unmasked fits checked for it); 0 to 11 rates a fit at
  the floor on R0, 1 to 12 on NQ.
- The 20 unmasked R0 fits, fitted over true exchangeability, the truth scaled by the fit's median ratio:
  Tyr–Val 43.7 to 46.7; the other 189 with 5th percentiles 0.48 to 0.66 and 95th 1.33 to 1.46 a fit.
  Over all 80 R0 fits, masked included, Tyr–Val 43.7 to 78.2.
- Unmasked R0, training lnL minus the truth's (the oracle's; the truth lies in GTR20+F{π*}'s family),
  n = 10: GTR20+F{π*} mean −339.33 (−397.11 to −270.20), GTR20+F −339.17; relaxed refits (`-me 1e-6`,
  300 rounds, n = 2) −338.70 and −334.05; NQC on the same data +18.22 at P and +171.52 at B.
- Held-out lnL minus the truth's, unmasked R0: GTR20+F{π*} −55.57; NQC at B −18.48.
- Largest exchangeability over Tyr–Val in the built-in matrices (the oracle's
  `builtin.exchangeabilities`): LG 42.7, WAG 24.9, JTT 60.1, Q.pfam 34.9, Q.bird 224.7, Q.mammal 103.6,
  Q.insect 57.2, Q.plant 82.0, Q.yeast 47.0, mtREV 643.7, HIVb 498.1; R0's truth 497.7.
- Round 1's fits on the true Tyr–Val's scale: the largest rate off the bound had medians 343 to 560 by
  masking strength, at most 622; 1.3 to 6.0 rates a fit sat at IQ-TREE's floor, rising with masking.

Scripts and outputs (each `.txt` headed by its command, commit and date):
`~/iqtree3-runs/sim2/20261124-round2-read/` (`extra.py`, `bounds.py`, `ranges.py`, `floors.py`,
`window_checks.py`, `firstlook.py`, `status2.py`), on round 2's collect (`collect-20261124T164042Z.tar.gz`, SHA-256
`c10cd93b…`).

Reproduced on the unmodified binary (decision 046's laptop demonstration, 2026-11-24, `explore/gtr20_window.py
fits` at `182a0dd7`; R0, 30 genes, replicate 1, the same settings, 1 thread): the default-seed GTR20+F{π*} fits on
the frozen binary `9a72950b` and the gate binary gave identical checkpoint rates and final lnL (−234155.85); Tyr–Val
10, exactly one rate (I–V) at 100, 8 at the floor; training lnL 32.09 below the truth's (oracle). Seeded fits on the
gate binary (Tyr–Val fixed at p): p = 0.02 and 0.05 gave 98.76 and 100.10 above the truth's, each with I–V at 100
too; with Tyr–Val alone scaled by 0.01 to 10 the oracle's lnL moved by at most 0.63, so that rate is weakly
identified at 30 genes and the free rates drift to the ceiling at little cost. Outputs:
`~/iqtree3-runs/sim2-gtr20/20261124T192220Z-182a0dd7-window/`; diagnosis
`~/iqtree3-runs/sim2-gtr20/20261124-smoke-diagnosis/` (`window_bounds.py`, `bounds.py` for the 3-gene smoke fits).

Anchors: `rescaleRates`; `ModelProtein::init`'s GTR20 and NONREV branches; `ModelMarkov::setBounds`,
`isUnstableParameters`, `readParameters`, `setReversible`'s conversion; `IQTree::optimizeModelParameters`'s
boundary warning; `PartitionModel::isUnstableParameters`; `--init-exchange` and `--init-model` parsing;
`dfpmin`, `lnsrch`; `PartitionModel`'s per-partition `ModelFactory`.

## Risk 23

Reported by the S3 plan's review, not verified: resuming an `-S` tree search with `+I+G` may
restore stale top-level model keys, for every model, because per-partition refits during the
search save to the top level (the model re-optimization in `IQTree::doNNISearch`) and the resumed
process keeps `opt_gammai` on, so the `+I+G` starts restore from a copy of the checkpoint
(`ModelFactory::optimizeParametersGammaInvar`; anchors in `CODE_PLAN.md` section 5). Outside S3's
validation paths, which use `-te`.

## Risk 24

Fits stop short of their peaks under IQ-TREE's stopping rules, for every model. `dfpmin` tests
its gradient relative to the log-likelihood with gtol fixed at 1e-4 and starts each call with no
curvature memory; the round loops stop at a gain below `-me` (0.01; ten times that under `-S`)
or at a cap (10 rounds under `-S`, 99 under `-p` and `-q`). On turtle's `-q` warm start each
round was one step gaining about 0.01, so where a fit stops is close to chance; NONREV's turtle
`-S` Level 1 fit hit its cap still gaining 3.5 per round (maintainer question 7). Comparing two
fits' scores (box against wide, NQC against NONREV, D08's arms) therefore needs a convergence
protocol. IQ-TREE's optimizer stays unchanged (decision 009). S3c's first part: NQC's endpoints
lie 0.77 to 4.2 below the converged peaks at default settings and 0.67 to 2.5 after 300 relaxed
rounds, more than the box costs; its later rounds average 724 evaluations, one step each. S3c's
profiling: NONREV's `-S` fits all ended at the 10-round cap, NQC's 10-set `-q` and `-p` fits at
the 99-round cap, so time per fit mostly follows where fits stop (NQC to NONREV 0.26 to 1.09)
while time per evaluation is equal.

SIM-1's pilot, measured while planning SIM-2 on 2026-11-19 (exploratory; commands in
`docs/agent/archive/SIM2_plan_2026-11-19.md`). An NQC Level 1 `-q` fit of 100 simulated genes (52 taxa,
NQ.bird truth), capped at 20 linked rounds (`-nparam 21`), still gained 1.06 in its last round. The
oracle evaluates the same objective and re-scored the accepted endpoint exactly (−852797.834). From that
endpoint, its L-BFGS-B with analytic gradients rose 90.2 in 90 iterations (20 minutes) and was still
rising. The capped fit's training lnL exceeded the true matrix's by 42.6 against 132.8 after the climb,
and Ĉ moved from 0.0610 to 0.0592 (truth 0.0577). At 30 genes the climb rose 15.6 in 124 iterations.
This is NONREV's process at scale, which bears on maintainer question 7. SIM-2 measures it for every fit
(decision 039, stages P and B) before any question is put.

Anchors: `dfpmin`'s scaled-gradient stop; the round loops, `modelEps` and the round counts.
