# CR reuse review, the reviewers' reports (archive, frozen 2026-11-21)

> **Not current.** Companion to `CR_reuse_review_2026-11-21.md`, which governs over these reports
> where they differ; CR's live records are `docs/agent/CODE_PLAN.md` section 4 ("CR") and
> `docs/agent/PLAN.md` risk 29.

These are the three read-only reviewers' final reports as delivered, at `91a0a45d`. They are the R
claims of `CR_reuse_review_2026-11-21.md`. The main file's V labels mark the claims rechecked against
the source, and its sections 1.6 and 3 correct or dismiss three of the reviewers' points:

- the chart's step comment (confirmed wrong);
- the C:313 comment (imprecise, the omission right);
- the fixed-model concern (dismissed).

---

## Appendix A: the model class, the hooks and the shared-file edits

**Answer.** NQC reuses most of IQ-TREE:

- ModelMarkov's whole non-reversible decomposition: unpacking, normalisation and the eigensystem.
- The transition matrices and likelihood kernels.
- ModelProtein's B/Z/J tip likelihoods and its model readers.
- The optimizer, unchanged (`minimizeMultiDimen`, `dfpmin`, `lnsrch`).
- The `+F` parsing routes, the checkpoint machinery, parameter counting (inherited
  `getNDim`/`getNDimFreq`) and PartitionModel's linked-fit loop.

It re-implements roughly 100 lines of base logic:

- the body of `optimizeParameters`, the derivative loop and the body of `targetFunk`;
- the `setVariables`/`getVariables`/`setBounds` pattern;
- NONREV's seed-reading control flow and the rooting block of `setReversible`;
- small formatters and constants.

Most of these copies are forced, either by decisions 009 and 015 (no edits to `optimization.*` or
`modelmarkov.*`) or by constraints 2, 3, 6 and 8. About 25 lines could reuse IQ-TREE now with
bitwise-identical output.

Three of the stated reasons for not reusing the base do not fully hold:

- **CODE_PLAN 2.2 item 4 (does not hold as a reason).** It says not to call `setReversible(false)`
  on a non-reversible seed because that zeroes `rates[]`. The zeroing is real but harmless at that
  point.
- **The comment at `C:313` (partly holds).** It justifies omitting the base `min_state_freq` guard,
  but π* is checked against that minimum before normalisation, not after.
- **Decision 006's "random draws" (partly holds).** IQ-TREE's reader draws random values only for
  tokens naming a defined distribution or "uniform"; any other token is an error.

Abbreviations: H `model/modelnonrevfixedeq.h`, C `model/modelnonrevfixedeq.cpp`, LH
`model/linkedmodelhooks.h`, MM/MMh `model/modelmarkov.cpp`/`.h`, MP `model/modelprotein.cpp`, MS/MSh
`model/modelsubst.cpp`/`.h`, OPT `utils/optimization.cpp`, PM `model/partitionmodel.cpp`, MF
`model/modelfactory.cpp`, CK `utils/checkpoint.h`, PA `main/phyloanalysis.cpp`, TL `utils/tools.cpp`,
CH `model/fixedeqchart.cpp`.

Effect classes:

- (i) bitwise-preserving;
- (ii) changes numbers or behaviour;
- (iii) needs a shared upstream edit;
- (iv) infeasible under a named constraint;
- (v) uncertain;
- n/a: nothing to reuse.

### Table A: members of ModelNonrevFixedEq and LinkedModelHooks

| # | NQC symbol | What it does | Nearest counterpart | Relation | Concrete differences | Reason given, and whether it holds | Could it reuse now? | Effect |
|---|---|---|---|---|---|---|---|---|
| 1 | `~LinkedModelHooks` LH:25 | virtual destructor | none | no counterpart | – | – | – | n/a |
| 2 | `LinkedModelHooks::startLinkedFit` LH:32 | call before each linked fit | none; `setBounds` is called once per linked fit (PM:787) but cannot see other partitions | no counterpart | – | Decision 005: only `optimizeLinkedModel` knows a fit starts. Holds | no | n/a |
| 3 | `LinkedModelHooks::linkedDerivative` LH:39 | derivative of the PartitionModel objective | `Optimization::derivativeFunk` (optimization.h:145, OPT:916-939) | no counterpart (interface) | takes the objective as a parameter | Decision 005: `minimizeMultiDimen` runs on PartitionModel (PM:810), so no model override is reached. Holds | only as a new `ModelSubst` virtual | (iii) `modelsubst.h`; a default keeps other models unchanged |
| 4 | `LinkedModelHooks::endLinkedFit` LH:46 | call after each linked fit, with the optimizer's score | none | no counterpart | – | Decisions 005, 009: the score is a local that PM:849 overwrites. Holds | no | n/a |
| 5 | `NAME` H:32, C:28 | "NQC" | `ModelSubst::name` (MSh:436), set from it at C:153 | no counterpart | – | – | – | n/a |
| 6 | `isName` H:35, C:106-110 | case-insensitive name match | `iEquals` (tools.h:3148, TL:210-219) | reimplements | `toupper` then `==` against `tolower` per character; same result for ASCII | none given | yes | (i) |
| 7 | `refuseAsMixtureComponent` H:39, C:112-115 | `outError` if a mixture component is NQC | none (`initMixture` has no component veto) | no counterpart | – | Decision 008: an `--mdef` mixture never shows "MIX". Unverified | no | n/a |
| 8 | `reportRun` H:43, C:663-701 | once-per-run export and summary per NQC model | `ModelSubst::report` (MSh:361; called only for POMO at PA:671-673, and from `PartitionModel::reportLinkedModel` PM:897-904, which has no callers); `writeInfo` (MM:361-380); linked-model loop PA:1392-1405 | no counterpart | needs the objective and tree; runs before the `OUT_IQTREE` return; writes files and cout | CODE_PLAN 1.2: no print site runs once per analysis; `writeInfo` runs per fit or round (MF:1748, PM:980-983). Holds | only if `reportModel` called `report()` for NQC and could reach the objective | (iii) PA; output would move into `.iqtree` and be lost when `.iqtree` is suppressed |
| 9 | constructor H:53, C:117-164 | protected base, refusals, π* parse, seed, encode, decompose | public ModelProtein constructor plus `init` (MP:1079-1088, 1113-1255) | wraps base (protected constructor MP:1090-1094) | C:151 and C:156 inline what `ModelMarkov::init(FREQ_USER_DEFINED)` (MM:354-359) does; `init` itself ends with that call (MP:1254) | Decision 002 holds: the public `init` reads a named model with tree adaptation on and unroots a rooted tree (MP:1127 → MM:1856 → MM:97-101). The other branches have side effects too: warnings (MP:1185-1188, 1214-1217); GTR20 refuses a non-reversible `--init-model` (MP:1200-1201); NONREV sends the target through `readStateFreq` (MP:1246-1249) | partly: call `init(FREQ_USER_DEFINED)` after `constructed = true` | (i): for user-defined frequencies `init_state_freq` only sums `state_freq` (MM:334-341); `writeInfo` at `VB_MAX` prints nothing for a 20-state non-reversible model (MM:368-376) |
| 10 | `targetFunk` H:57, C:304-317 | −lnL at x; 1e30 on a failed build | `ModelMarkov::targetFunk` (MM:1091-1118) | reimplements | adds a counter; skips `clearAllPartialLH` on a failed build; returns 1e30; omits the `min_state_freq` guard (MM:1103-1108) | C:313's reason partly holds: `parseTarget` checks raw entries (CH:384), then divides by the sum (CH:399). A normalised entry can fall below `min_state_freq` by a factor down to 1/(1+1e-6); the base guard would then return 1e30 at every evaluation. Constraint 2 needs a 1e30 return without computing, which the base cannot give. Holds | no pre-likelihood validity hook exists | (iv) constraint 2 |
| 11 | `derivativeFunk` H:60, C:319-321 | single-model gradient | OPT:916-939 | reimplements (via #30) | step 1e-4·max(1,\|x\|), backward at the upper bound; base uses 1e-4·\|x\|, or 1e-4 at 0 | D02, decision 005: the base step collapses near 0 (OPT:928-929). Mechanism holds; quoted error figures unverified | only with a step hook in Optimization | (iii) `utils/optimization.{h,cpp}`; barred by decision 009 |
| 12 | `setBounds` H:62, C:295-302 | box in log-ratio units, `bound_check` false | MM:1139-1164; same 5-line pattern in ModelUnrest (`modelunrest.cpp:85-92`) and ModelLieMarkov (`modelliemarkov.cpp:886-897`) | reimplements (IQ-TREE's per-class pattern) | constants | Decisions 007, D03. Holds: the base hard-codes `MIN_RATE`/`MAX_RATE` (MM:1146-1147), which exclude every negative coordinate | no | (iv) constraint 8 |
| 13 | `optimizeParameters` H:65, C:357-395 | the base routine statement for statement, plus `fit_starting`, the optimizer score and `checkFit` | MM:1166-1244 | reimplements | drops the unused `variables2` (MM:1182, 1196, 1240); drops the reversible estimated-frequency block (MM:1227-1230), unreachable for a non-reversible model | CODE_PLAN 2.2: the base discards the optimizer's score. Holds (assigned at MM:1199, overwritten at MM:1234) | only if ModelMarkov kept the score in a member | (iii) `model/modelmarkov.{h,cpp}`; barred by decision 015. Without the score, (ii): the `optimizer_score` key is lost |
| 14 | `decomposeRateMatrixNonrev` H:68, C:510-537 | build Q from trial coordinates, pack into `rates[]`, call base, check π | MM:1246-1387 | wraps base | packing loop C:530-533 equals `setFullRateMatrix`'s non-reversible branch (MM:867-871) | Decision 015 holds: with user-defined frequencies and `-optfromgiven` refused, the base skips the π reset and solve (MM:1252, 1266) | packing could be `ModelMarkov::setFullRateMatrix(rowMajorCopy.data(), state_freq)` | (i), at a 400-entry copy; marginal |
| 15 | `getStateFrequency` H:71, C:539-545 | π* copied exactly | MM:875-892 | reimplements (base without renormalisation) | the base multiplies by 1/sum (MM:887-891) | Decision 015 and constraint 3. Holds: the kernel takes root frequencies from this function (`phylokernelnonrev.h:682, 1153`) | no | (iv) constraint 3 |
| 16 | `setStateFrequency` H:74, C:547-551 | accept only π* | MM:894-912 | wraps base | guard; precedent ModelUnrest's no-op (`modelunrest.cpp:113`) | constraint 3. Holds | already calls base | n/a |
| 17 | `adaptStateFrequency` H:75, C:553-560 | base before construction; after, accept only π* and do nothing | MM:914-930 | wraps base | after construction never calls base | constraint 3. Holds: the base rescales `rates[]` by f/f (MM:917-928), not bitwise identity. Callers PM:165 (estimated or empirical only, PM:119), `modelmixture.cpp:3981` (refused) | no | (iv) constraint 3 |
| 18 | `setRateMatrix` H:78, C:562-566 | stop after construction | MM:852-856 | wraps base | guard | CODE_PLAN 2.2: every caller is refused. Holds: MF:661 and 682 (site-specific, refused); `modelmixture.cpp:4399` (refused); PA:4429, 4456 (no callers); `iqtree.cpp:2007` (4 states only); MF:1265 (comment) | calling the base would be silently discarded at the next decomposition | (ii) an error becomes a silent no-op |
| 19 | `setFullRateMatrix` H:79, C:568-572 | same | MM:858-873 | wraps base | guard | only caller PM:111 sits behind `ASSERT(0)` (PM:97). Holds | as #18 | (ii) |
| 20 | `startCheckpoint` H:81, C:574-576 | own struct name | MP:1257-1259 | reimplements (one-line pattern) | name | IQ-TREE convention | dropping it moves every key under "ModelProtein" | (ii) key paths change |
| 21 | `saveCheckpoint` H:82, C:578-606 | 17-digit string keys, then `ModelProtein::saveCheckpoint` | MP:1262-1269; ModelLieMarkov (`modelliemarkov.cpp:606-613`) | wraps base | the base call also writes `rates` (380 values at 10 digits) and `state_freq` into the NQC struct; both redundant | Decision 003 holds: `CKP_PRECISION` is 10 (CK:22, 328, 352), and the +I+G path saves and restores the model (MF:1488, 1504) | only with a precision argument on `put`/`putArray` | (iii) `utils/checkpoint.h`; a defaulted argument keeps other models unchanged |
| 22 | `restoreCheckpoint` H:83, C:608-656 | base restore, verify, rebuild from 17-digit coordinates | MP:1271-1282; ModelLieMarkov (615-625) | wraps base | the base decomposes twice from the pre-restore coordinates (MS:58, MP:1278), then discarded; parses with `convert_double_vec`/`convert_double` (TL:430-446, 383-393) | Decision 015 holds: every decomposition goes through the override, so 10-digit `rates[]` never become Q | parsing could use `Checkpoint::getVector` (CK:230-257) if keys used ", " | (ii) checkpoint and `.info` strings change |
| 23 | `startLinkedFit` C:414-442 | linked-compatibility signature check | none | no counterpart | – | D06, decisions 005, 011. Holds | see Table C, C11 | n/a |
| 24 | `linkedDerivative` H:90, C:323-325 | scaled differences on PartitionModel | OPT:916-939 | reimplements (via #30) | as #11 | as #11 | as #11 | (iii) |
| 25 | `endLinkedFit` H:93, C:444-458 | re-score the linked set, record in every object, report once | partly the base re-score PM:843-850 | no counterpart | the base sums the whole super tree, only when the point changed; NQC sums only linked NQC partitions in `part_order` | C:445-449 matches `PartitionModel::targetFunk`'s summation order (PM:314-333). Holds | no: a different quantity | (ii) |
| 26 | `setVariables` H:96, C:279-283 | copy accepted coordinates out | MM:1018-1036; ModelLieMarkov (899-913) | reimplements (pattern, own array) | which array | Decision 003 holds: the base is wired to `rates[]` (MM:1030), which the override overwrites with Q | no | (iv) constraint 2 |
| 27 | `getVariables` H:97, C:285-293 | exact compare, copy into trial | MM:1038-1089; ModelLieMarkov (919-930) | reimplements | same comparison semantics (MM:1056-1058) | a failed proposal must not reach `rates[]` | no | (iv) constraint 2 |
| 28 | `refuseUnsupported` C:166-220 | 17 refusals | none overall; `posRateHeterotachy` (`modelfactory.h:41`, MF:76-78) for the +H test | no counterpart | raw `find("+H")` (C:190-191); IQ-TREE's helper also catches "*H" and skips "+H" followed by a letter | Decision 008 as policy; nothing given for skipping the helper | +H test: yes | (ii) edge strings only |
| 29 | `readSeed` C:222-277 | LG read for references; `--init-model` read; seed at π*; set π*; root the tree | NONREV branch of `ModelProtein::init` (MP:1213-1239; same init-model read in the GTR20 branch, MP:1191-1198); `setReversible` rooting (MM:145-149) | wraps base readers; reimplements the branch's control flow and the rooting block | no `rescaleRates`; C:271-276 copies the rooting rule by hand | Decision 002 (no `rescaleRates`) holds. CODE_PLAN 2.2 item 4's fact (MM:133-137 zeroes `rates[]`) is true but does not hold as a reason: `rates[]` has been consumed into Q0 (C:253-260), and the constructor's decomposition rewrites it from z (C:156 → C:530-533) before any read | yes: an unconditional `setReversible(false)`; the init-model read could become a protected ModelProtein helper shared with the GTR20 and NONREV branches | rooting: (i). Helper: (iii) `modelprotein.cpp`, a pure move for GTR20 and NONREV |
| 30 | `scaledDerivative` C:327-355 | forward differences, start-score capture, `-vv` print | OPT:916-939 | reimplements | step rule; objective as a parameter | as #11 | as #11 | (iii) |
| 31 | `checkFit` C:404-412 | re-score, warn, record, report | none | no counterpart | – | decision 009 | – | n/a |
| 32 | `rescore` C:460-468 | π guard, decompose, clear partials, compute likelihood | base re-score MM:1231-1235, also copied at C:383-387 | reimplements | adds the guard | decision 009 | when `changed` or the score is −1e30, C:386 has already computed the same likelihood | (i) at `-T 1`; saves one evaluation per fit |
| 33 | `recordFit` C:470-477 | update the fit record | none | no counterpart | – | – | – | n/a |
| 34 | `reportFit` C:479-488 | `-v` record and fall warning | role of `writeInfo` (MM:361-380) | no counterpart | – | CODE_PLAN 2.2: replaced S2's `writeInfo` because the `-v` round loop repeated it. Holds (MF:1748, PM:982) | as a `writeInfo` override | (ii) printed at every `writeInfo` site |
| 35 | `measureState` C:490-500 | bound counts and residuals | none; the counts duplicate `firstOrderTest`'s (CH:423-432) | no counterpart | – | – | – | n/a |
| 36 | `writeState` C:502-508 | print the record | none | no counterpart | – | – | – | n/a |
| 37 | `printMatrix` C:658-661 | 17-digit rows of Q | `reportModel`'s full-Q rows (PA:654-658), 6 decimals | reimplements | precision, destination, separator | decision 021 | no | (ii) |
| 38 | `ExportDiagnostic` H:141-148 | result struct | none | no counterpart | – | – | – | n/a |
| 39 | `exportRun` C:703-715 | runs the export steps in order | none | no counterpart | – | – | – | n/a |
| 40 | `writeQmat` C:717-731 | `.qmat` file | `reportModel`'s "Full Q matrix ... (can be used as input for IQ-TREE)" block (PA:650-665); the ofstream-exception idiom (PA:2600-2640) | reimplements the content; reuses the idiom | 17 digits, own file, no header | Decision 010 holds: the report prints Q at 6 decimals (PA:637, fixed mode at PA:1535) | only by raising `reportModel`'s precision | (iii) `main/phyloanalysis.cpp`; changes every model's `.iqtree` |
| 41 | `exportDiagnostic` C:733-803 | central differences, then rebuild | none: IQ-TREE has only forward differences | no counterpart | the call at C:787 passes `FIRST_ORDER_GTOL`; `TOL_RATE` (MMh:31) has the same value | – | constant only | n/a (the `TOL_RATE` swap is (i)) |
| 42 | `writeSummary` C:805-827 | end-of-run summary | none | no counterpart | letters loop, see #45 | – | – | n/a |
| 43 | `writeExportInfo` C:829-898 | `.info` file | none; reuses the file idiom | no counterpart | – | – | – | n/a |
| 44 | `CHART`, `CHART_VERSION`, `DERIVATIVE` C:30-32 | metadata | none | no counterpart | – | – | – | n/a |
| 45 | `STATES` C:33, used at 807-809 and 836-838 | amino-acid letters | `symbols_protein` (`alignment.cpp:34`; also copied in `modelcodon.cpp:579`) via `Alignment::convertStateBack` (`alignment.cpp:1916-1917`) | reimplements (duplicate constant) | none | none given | yes | (i) |
| 46 | `RECHECK_LNL` C:35 | decision 038's line | none | no counterpart | – | – | – | n/a |
| 47 | `number17` C:37-42 | 17-digit double | `convertDoubleToString` (TL:204-208, precision 6); the chart's private `str` (CH:110-116, 17 digits, classic locale) | reimplements | precision; no locale set | decisions 003, 010 need 17 digits | upstream only with a precision argument; inside the fork it duplicates CH:110-116 | (iii) `utils/tools.{h,cpp}` |
| 48 | `numbers17` C:44-49 | join doubles | `Checkpoint::putArray` (CK:345-358) | reimplements | 17 digits, "," separator | decision 003. Holds | only with a precision argument | (iii) `utils/checkpoint.h` |
| 49 | `integers` C:51-56 | join ints | `Checkpoint::putVector` (CK:366-379), lossless for ints | reimplements | "," against ", " | none given | yes, for `references` and `bounds` | (ii) format change |
| 50 | `score17` C:59-62 | NaN printed as "unknown" | none | no counterpart | – | – | – | n/a |
| 51 | `oneLine` C:66-78 | collapse whitespace | `trimString` (TL:7081-7084) only trims the ends | no counterpart | – | – | – | n/a |
| 52 | `jointRateModel` C:82-104 | non-`+F` component in `--model-joint` | the joint string's `+F` extraction (MF:284-304), same `find_first_of("+*")` idiom, no callable function | reimplements the idiom | skips braces | decision 008 policy | no function exists | n/a |
| 53 | `BUILD_FAILED_AFTER_FIT` C:397-398 | message | none | no counterpart | – | – | – | n/a |
| 54 | `sameBits` C:400-402 | `memcmp` wrapper | – | no counterpart | – | – | – | n/a |

Data members that duplicate base state:

| # | Member | Base state it duplicates | Assessment | Effect |
|---|---|---|---|---|
| 55 | `models` H:104 | `ModelProtein::models_block` (`modelprotein.h:101`, private; set at MP:1093) | Used only in `readSeed` (C:225, 239), which only the constructor calls, where `models_block` is a parameter. Pass it to `readSeed` and drop the member, or make the base member protected ((iii), `modelprotein.h`, no behaviour change) | (i) |
| 56 | `target` H:107 | `state_freq` (MSh:450) | the immutable reference for the guards at C:463, 535, 548, 558. Holds | (iv) constraint 3 |
| 57 | `Q` H:116 | `rates[]` (MMh:494) and `rate_matrix` (MMh:540) | a third copy of Q, the chart's output buffer, used only at C:515 and 533; `rate_matrix` cannot stand in, because the base rewrites and rescales it. Could be a local | (i), at one heap allocation per evaluation |
| 58 | `trial`, `accepted` H:114 | `rates[]` as the optimizer's state in the base (MM:1058) | Decision 003 and constraint 2. Holds | (iv) constraint 2 |
| 59 | `lower`, `upper` H:110 | `Params::nqc_domain_lo/hi` (`tools.h:1709`) | log of the option values; needed for the checkpoint and the linked check | n/a |
| 60 | `constructed` H:102 | none | Holds: during `readSeed`, `readParametersString` reaches the virtual decomposition (MM:1873) | n/a |

### Table B: base-class virtuals a fitted NQC object inherits unchanged

| Inherited behaviour | Behaviour for NQC | Correct? | Blocked or bypassed? | Should NQC rely on it? |
|---|---|---|---|---|
| `ModelMarkov::getNDim` MM:964-976 | non-reversible: `num_params` (360), 0 when fixed | yes (constraint 5) | – | yes |
| `ModelMarkov::getNDimFreq` MM:978-1000 | 0 for user-defined frequencies | yes | – | yes |
| `ModelMarkov::getName` MM:241-245 | "NQC" (`freqTypeString` returns "" for protein user-defined, MM:204-206) | yes for linking; every NQC partition shares the name whatever its target; `startLinkedFit` catches a mismatch | – | yes |
| `ModelProtein::getNameParams` MP:1339-1354 | "NQC", without the target | lossy: the `.iqtree` model line cannot re-run the fit (maintainer question 8) | – | by design |
| `isReversible`, `getNumRateEntries` MMh:63, MM:154-159 | false; 380 | yes | – | yes |
| `ModelMarkov::decomposeRateMatrix` MM:1389-1395 | dispatches to NQC's override | yes; the hook NQC relies on | – | yes |
| `computeTransMatrix`, `computeTransMatrixNonrev`, `computeTransDerv`, `computeTrans(t,i,j)` MM:507-511, 463-505, 732-768, 570-588 | built from `rate_matrix` and the complex eigensystem, with the scaling-and-squaring fallback | yes (D06) | – | yes |
| `computeTrans(t,i,j,d1,d2)` MM:590-604 | uses the real eigen arrays regardless of `is_reversible`; stale for non-reversible models | wrong, identically for NONREV; callers `alignmentpairwise.cpp:320-382`; whether NQC runs reach it unverified | not blocked | upstream issue |
| `getRateMatrix` MM:847-850 | copies `rates[]` | adequate; `reportModel` switches to `getQMatrix` above 20 parameters (PA:632-652) | – | yes |
| `getQMatrix` MM:932-937 | copies normalised `rate_matrix` | yes | – | yes |
| `getFreqType` MMh:272 | user-defined | yes | – | yes |
| `ModelProtein::computeTipLikelihood` MP:1356-1372 | B/Z/J ambiguity codes | yes; the reason NQC subclasses ModelProtein (decision 002) | – | yes |
| `isUnstableParameters` MM:1120-1137 | flags any derived rate above 99; with unit mean rate, q_ij ≤ 1/(2π*_i), so it can fire when π*_i < 1/198; warning at `iqtree.cpp:2342-2346` | meaningless in chart coordinates | not blocked | override, or leave it; either way (ii) |
| `writeInfo` MM:361-380 | prints nothing for 20 states | yes; NQC relies on the silence | – | yes |
| `ModelSubst::report` MSh:361 | empty; reached only for POMO | n/a | – | – |
| `ModelProtein::printMrBayesModelText` MP:1374-1422 | writes "gtr" with free frequencies | wrong model in the MrBayes block, as for NONREV; not refused | no | user-facing note |
| `setReversible` MM:75-152 | called only in `readSeed` | yes | – | yes (#29) |
| `ModelProtein::init` MP:1113-1255 | never called | – | bypassed through the protected constructor | correct to bypass |
| `ModelProtein::readRates(istream&)`, `readStateFreq(istream&)` MP:1284-1336, MM:1739-1759 | used to read LG and `--init-model`; the target never passes through them | yes | – | yes |
| `readRates(string)`, `readStateFreq(string)` MM:1702-1737, 1761-1796 | unreachable | – | refused upstream of them | – |
| `ModelMarkov::setRates` MM:1949-1952 | `ASSERT(0)`, never called | – | unused | (i), saves no code |
| `ModelSubst::fixParameters` MSh:82-86 | linking uses it | yes | – | yes |
| `Optimization::restartParameters` OPT:726-748 | never restarts (`bound_check` false) | yes | – | yes |
| `optimFunc`, `optimGradient` OPT:1261-1268 | L-BFGS-B path, commented out at its callers | fine | – | – |
| `writeParameters`, `update_eigen_pointers`, `get_num_states_total`, mixture virtuals, `getMutationModel`, `isFused` | defaults | yes | mixtures refused | – |
| non-virtual `scaleStateFreq` MM:1002-1016 | every caller gated on reversible plus estimated frequencies (MM:1227, PM:839, PA:551), or mixtures, or codons | unreachable | the guards at C:463 and C:535 would stop the run | – |
| non-virtual `init(StateFreqType)` MM:354-359 | not called | – | inlined by NQC | could be called, (i) (#9) |

### Table C: each shared-file edit

| Edit | What it adds | Existing extension point, or new? | Could an existing IQ-TREE mechanism do it? |
|---|---|---|---|
| C1 PA:61 include | header for C2 and C3 | – | – |
| C2 PA:165-176 citation flag | NQC under `--model-joint` gets nQMaker's citation | modifies an existing condition | No helper exists (`aa_model_names_nonrev`, `phylotesting.cpp:165`, lists only NQ.*; `validModelName`, MM:1965-1968, covers UNREST and Lie-Markov). Behaviour-preserving without NQC. Its prefix parse differs from `jointRateModel`'s, so the fork has two parsers |
| C3 PA:1417-1418 `reportRun` call | once-per-run export | new call site | No. `report()` is reached only for POMO and writes `.iqtree` after the `OUT_IQTREE` return; `writeInfo` runs per fit or round |
| C4 `model/CMakeLists.txt` +3 | sources | required: no glob | – |
| C5 `modelmixture.cpp:18` include | – | – | – |
| C6 `createModel` 3231-3235 | NQC name dispatch, protein check | new branch; reuses `createModel`'s own brace and `+F` parsing (3186-3219) | No. `validModelName`/`getModelByName` (MM:1954-1968) has no models block, sits after the PoMo branch and is in `modelmarkov.cpp` (decision 015). ModelProtein's name handling can only construct a `ModelProtein` (3264), so NQC's overrides would become flag branches inside the base: (iii)/(iv) |
| C7 `initMixture` 3453, 3527 | refuse NQC as a component | new | No component veto exists; string detection was rejected by decision 008 (premise unverified) |
| C8 `modelprotein.h:91-98`, MP:1090-1094 | protected constructor | new | No: the public constructor always runs `init`, and every branch adapts the tree, warns, or refuses a non-reversible `--init-model`. The same edit could make `models_block` protected (#55) |
| C9 PM:26 include | – | – | – |
| C10 `partitionmodel.h:166`, PM:754-760 | `derivativeFunk` override delegating to the hooks | overrides an existing virtual on PartitionModel; new interface on the model side | Only by a new `ModelSubst` virtual, a larger shared edit. For other models it calls `Optimization::derivativeFunk` exactly as before, at one `dynamic_cast` per gradient |
| C11 PM:794-807 start hook and list | pre-fit check | new | `setBounds` (PM:787) could carry the start flag but cannot see other partitions. The list loop copies `targetFunk`'s name filter (PM:314-320) and already exists in `reportRun` (C:682-689). NQC could build the list itself: (i), same list and order, about 12 fewer shared lines |
| C12 PM:813 `optimizer_score` | keeps the optimizer's score | local variable | needed because PM:849 overwrites `score` |
| C13 PM:851-852 end hook | post-fit check | new | no post-fit hook carries the score; `writeInfo` at PM:980-983 runs once per call and only at `VB_MIN` |
| C14 `tools.h:1706-1709`, TL:3167-3178, 6013-6016, 7279-7280 | `--nqc-domain` | existing option mechanism (copying `--min-freq`, TL:3156-3165; reuses `convert_double_vec`) | No model-parameter bound option exists. Model-string braces (`modelmixture.cpp:3186-3205`, refused at C:178-179) would avoid editing `tools.*`, at the cost of a user-visible syntax change: (ii) |

### Section D: duplication outside the tables

- **Parsing.**
  - `parseTarget` against `readStateFreq(string)` (MM:1761-1796). Decision 006's reason holds for
    NaN (`strtod` accepts it) and for renormalisation with only a warning (MM:1789-1795); it partly
    holds for "random values" (TL:684-699). The base also accepts zeros and has no minimum.
    Reusing it is (iv): D07, constraint 3.
  - `restoreCheckpoint` calls `convert_double_vec` and `convert_double` without a try. They throw
    `std::string`, so a malformed key ends in `std::terminate` rather than `outError` (C:639, 647).
  - The +H test should use `posRateHeterotachy`. Raw `"+T"` matches IQ-TREE's own practice
    (PA:4910).
- **Printing.**
  - `number17` duplicates the chart's `str`.
  - `STATES` and the letters loop appear twice (C:807-809, 836-838).
  - The `diagnostics[7]` array is built twice (C:582-584, 833-835), `domain[2]` five times and
    `score_keys` twice.
  - `printMatrix` and `writeQmat` print the same rows with different separators.
  - The version string copies `main.cpp:132`.
- **File writing.** Reuses IQ-TREE's idiom; the `.qmat` content duplicates `reportModel`'s full-Q
  block at higher precision.
- **Parameter counting.** Inherited, no duplication.
- **Constants.** `FIRST_ORDER_GTOL` equals `TOL_RATE`; `STEP_ETA` equals `ERROR_X` (file-static in
  OPT).
- **Layout conversions.**
  - Q to `rates[]` (C:530-533) equals MM:867-871.
  - Packed upper triangle to a symmetric R is done twice (C:228-231, 248-250), copying base unpack
    loops (MM:1612-1619, 946-952). `getQMatrix` (reversible branch, `state_freq` = π*) returns
    R·π* at unit mean rate, but divides where `reversibleSeed` multiplies by the reciprocal: (ii).
  - Off-diagonal `rates[]` to Q0 (C:253-260) copies MM:1256-1264. After the non-reversible
    `readParameters`, `rate_matrix` already holds Q0 rescaled at its own π (MM:1269-1282,
    1838/1873); reading it is (ii).
- **Partition loops.** The list loop is duplicated (PM:799-805, C:682-689). `reportRun`'s swap of
  the `model` member (C:694-699) copies `optimizeLinkedModels` (PM:873-891).
- **Redundant base work.**
  - Each save writes the derived `rates` and `state_freq` keys at 10 digits.
  - Each restore runs two decompositions that are then discarded.
  - Each single-model fit evaluates the likelihood twice (#32).
- **Internal duplication.** `measureState`'s bound counts duplicate `firstOrderTest`'s.

### Section E: summary

- **Counts by relation** (54 rows): calls base unchanged 0; wraps base 9; reimplements 19; no
  counterpart 26.
- **Counts by effect** (54 rows): (i) 6; (ii) 9; (iii) 10; (iv) 6; (v) 0; n/a 23.
- **Data members:** (i) 2, (iv) 2, n/a 2.

Largest opportunities:

1. **The `optimizeParameters` copy**, about 30 lines: (iii), barred by decision 015.
2. **The derivative loop**, about 14 lines: (iii), barred by decision 009.
3. **`readSeed`:** the rooting block (i); the LG and `--init-model` reads as a ModelProtein helper
   (iii).
4. **The linked list built inside NQC:** (i), about 12 shared lines.
5. **Small (i) items**, about 15 lines in total:
   - `iEquals`;
   - `convertStateBack`;
   - pass `models_block`;
   - `init(FREQ_USER_DEFINED)`;
   - a qualified `setFullRateMatrix`;
   - `TOL_RATE`;
   - the re-score.

Departures judged necessary:

- **Constraint 3 and D07:** `getStateFrequency`, the immutable target and the guards, and
  `parseTarget`.
- **Constraint 2 and decision 003:** the 1e30 path and the trial/accepted arrays.
- **Constraint 2:** no `computeStateFreqFromQMatrix` inside the build (its `ASSERT` aborts in
  Release).
- **Constraint 8:** `setBounds`.
- **Constraint 6:** the 17-digit keys.
- **Decisions 009 and 015 only:** the derivative step and the `optimizeParameters` copy.
- **Decision 016:** the protected constructor.
- **No existing mechanism:** the `createModel` branch, the hooks and the `reportRun` call site.

Unverified:

- decision 008's premise that an `--mdef` mixture never shows "MIX";
- the include-cycle reason for a separate hooks header (LH:9-11);
- D02's error figures;
- whether NQC reaches the non-reversible `computeTrans` with derivatives;
- whether NONREV with a non-reversible `--init-model` leaves an unrooted tree unrooted;
- whether alternative seed rounding flips a reference tie;
- whether the `isUnstableParameters` warning has fired;
- that the global locale is never changed;
- that ModelFinder builds no NQC model before `dating_mf` is set;
- whether a single-alignment NQC model can be fixed at export. The main session dismissed the
  linked case: `PartitionModel::targetFunk` unfixes each partition model, PM:318-321.

---

## Appendix B: the chart (`model/fixedeqchart.{h,cpp}`)

**Answer.** About 45 of the module's roughly 580 lines repeat work IQ-TREE already does: one
stationary solve, one diagonal fill, one reversible seed, two constants, and the formulas for the
derivative step and the gradient test. The target parser and its helpers (about 70 lines) differ on
purpose (decision 006). Everything else has no counterpart or would break a requirement.

NQC already reuses IQ-TREE's numerical back end:

- the packed `rates[]` layout;
- the diagonal fill and unit-mean-rate normalization, through the unchanged
  `decomposeRateMatrixNonrev`;
- the eigensystem, P(t) and the likelihood;
- the optimizer and `fixBound`;
- the displacement idiom;
- the QR formulation of the stationary solve.

Only `reversibleSeed` can be reused by changing chart code alone, and that changes the seed's last
bits. Moving `computeStateFreqFromQMatrix` or the diagonal loop into a shared helper would keep every
number bitwise, but needs a `modelmarkov.cpp` edit, which decision 015 rules out.

The counterpart files (`modelmarkov.*`, `eigendecomposition.*`, `optimization.*`, `alisimulator.cpp`,
`ratefree.cpp`, `lbfgsb_new.cpp`) are byte-identical to `63c330d9`. Eigen 3.4.0 was read in WSL.

### Table A: every function in fixedeqchart

| Chart function | What it computes | Nearest counterpart | Relation | Differences | Reason, and whether it holds | Reuse now? | Effect |
|---|---|---|---|---|---|---|---|
| `freeDestinations` cpp:142-148 | row i's destinations other than i and ref[i], ascending | none (IQ-TREE's packing, skip diagonal, MM:1256-1262, 868-871, with the reference also skipped) | no counterpart | nothing calls it; `jumpMatrix` and `inverse` inline the rule | D01 order. Holds | dead code | n/a |
| `jumpMatrix` cpp:150-188 | per-row softmax with reference logit 0, max-shifted; rejects an exact 0 | related: `RateFree::getVariables` (ratefree.cpp:440-451); FREQ_ESTIMATE (MM:1061-1067, 1227-1230); max-shifted exp (PM:702-703) | related, different | RateFree uses no log or exp, ratio bounds, member state, no underflow check | D01; constraint 2. Holds | no | (iv) |
| `solveStationary` cpp:190-204 | x M = 0, Σx = 1, from [1ᵀ; Mᵀ]x = e₁ by `colPivHouseholderQr`; rejects non-finite | `computeStateFreqFromQMatrix` MM:2119-2131 (free function, declared only at MM:38) | reimplements | same system, right-hand side and solver; upstream maps row-major as column-major (implicit transpose), the chart writes `M.transpose()`; upstream `ASSERT`s abs(Σπ−1) < 1e-4, which aborts; equal M gives bitwise-equal output (derived from Eigen 3.4.0 `ColPivHouseholderQR.h:587-607`) | decision 012 holds as a project choice; no process exit holds through the `ASSERT`; thread safety and static state do not hold as reasons | the unit-test project cannot link it | (iii): a shared header, keeping the `ASSERT` upstream; pure move, edits `modelmarkov.cpp` (decision 015). A direct call is (iv) |
| `generatorFromJump` cpp:206-231 | ν of K, ν > 0; q_ij = (ν_i/π_i)K_ij; diagonal; rate checks | none for the row scaling; IQ-TREE column-scales (MM:114-119, 914-930) | related, different | column scaling does not keep π stationary for a non-reversible Q (D09 item 6) | constraint 1, D01. Holds | no | (iv) |
| `build` cpp:233-239 | `jumpMatrix` then `generatorFromJump` | NONREV's map: raw rates, π re-solved (MM:1256-1267) | related, different | different chart; NONREV re-solves π | constraints 1, 3. Holds | no | (iv) |
| `inverse` cpp:241-268 | K_ij = q_ij/(−q_ii); z = log(K_ij/K_ir) | none for log-ratios; the K step also in `AliSimulator::extractRatesJMatrix` (alisimulator.cpp:2728-2739) | related | AliSim multiplies by 1/(−q_ii), a different rounding, and is a member | D01; decision 012. Holds | no | (iv) |
| `references` cpp:270-278 | row argmax of off-diagonals, ties to the lower index | the inline argmax for `highest_freq_state` (MM:1188-1192; copied at C:373-375) | related (idiom) | one-dimensional | D01, decision 011 | no callable unit | (iii), not worth it |
| `reversibleSeed` cpp:280-302 | q_ij = R_ij π_j, diagonal, divided by Σπ_i(−q_ii) | `EigenDecomposition::computeRateMatrix` (eigendecomposition.cpp:478-519) via `getQMatrix`'s reversible branch (MM:932-962); `setReversible(false)` (114-119); Eigen3 path (1503-1514) | reimplements | R·π against π·R bitwise equal; same row-sum and mean-rate order; the chart divides, upstream multiplies by delta = 1/sum; `ignore_state_freq` skips the π multiplication once `setReversible(false)` has run; `getQMatrix` zeroes rates where π ≤ `ZERO_FREQ` | decision 012 holds for the chart; from the class no constraint forbids `getQMatrix` | yes from the class: set `state_freq` to π*, call `getQMatrix` | (ii): seed's last bits change, then start coordinates, every fit, `.qmat` and checkpoint bytes |
| `transferSeed` cpp:304-323 | K = q0_ij/(−q0_ii), then `generatorFromJump` at π | the K step (AliSim); `adaptStateFrequency` (914-930) | related | `adaptStateFrequency` does not preserve stationarity; for reversible Q0 = Rπ it returns `reversibleSeed(R, π)` exactly in exact arithmetic (derived) | constraint 1; D05. Holds | no | (iv) |
| `diagnostics` cpp:325-350 | ‖πQ‖∞, scaled row sums, mean-rate error, minimum rate and flux, ν residual, exit-rate mismatch | none; related checks MM:2128-2129, 1696-1697, 1834-1843, 474-478, 489-500 | related, different | IQ-TREE checks other quantities and aborts, warns or switches method | decision 021. Holds | no | (iv); its K step (343) is a third copy in the chart |
| `parseTarget` cpp:352-403 | 20 ASCII decimals; > 0 and ≥ minimum; sum within 1e-6; divided once | `readStateFreq(string)` (1761-1796); `convert_double_with_distribution` (tools.cpp:408-428); `normalize_frequencies` (772-795) | reimplements (different rules) | `strtod` accepts NaN and hex; non-numeric tokens become draws through the global RNG (tools.cpp:684-726); allows 0; renormalizes by the reciprocal when abs(sum−1) ≥ 1e-7 with a warning; `outError` | decision 006, D07, no exit, decision 012. Hold | no | (iv); would also change π*'s bits |
| `scaledStep` cpp:405-411 | h0 = η·max(s, abs(x)), negated past the upper bound; returns (x + h0) − x | `Optimization::derivativeFunk`'s inline step (OPT:927-931) | related (deliberate variant) | floor at 1; backward at the upper bound; same idiom; for abs(x) ≥ 1 below the bound equal to the legacy probe bitwise (derived) | D02, decisions 005, 009. Holds | no | (iii), forbidden by decision 009 |
| `firstOrderTest` cpp:413-444 | decision 009's label; projected gradient; s_k scaling; argmax and counts | `dfpmin`'s stop test (OPT:849-857); `projgr` (lbfgsb_new.cpp:3908-3951) | related (dfpmin's formula done again) | same expression inside the box; adds bound projection, outside count, argmax; `projgr` is a different statistic | decision 009. Holds | no | (iii), forbidden by 009; `projgr` (ii) |
| `statusText` cpp:120-140 | text per Status | none | no counterpart | n/a | constraint 2. Holds | n/a | n/a |
| `validReferences` cpp:23-31 | size and range checks | none | no counterpart | n/a | n/a | n/a | n/a |
| `validTarget` cpp:33-38 | π finite and > 0 | the `targetFunk` guard (MM:1103-1107) allows 0; `init_state_freq` sum check (334-341) | related | accepts zeros | D07. Holds | no | (iv) |
| `square` cpp:40 | shape check | none | no counterpart | n/a | n/a | n/a | n/a |
| `checkRates` cpp:43-54 | off-diagonals finite and > 0 | `readRates` throws or `outError`s (1686-1687, 1723-1724); NONREV bounds at `MIN_RATE` (1146) | related | those allow 0 or exit | constraint 2; no exit. Holds | no | (iv) |
| `completeDiagonal` cpp:57-71 | `checkRates`, then q_ii = −Σ q_ij | unpack loop of `decomposeRateMatrixNonrev` (1256-1264); `computeRateMatrix` (495-499, 511) | reimplements | same order, bitwise equal; in the per-evaluation path the class copies only off-diagonals (C:530-533), so the chart's diagonal goes unused there | decisions 012, 015, as project choices | no callable unit | (iii); rejected by decision 015 |
| `separator` cpp:73-75 | `,` `/` and six ASCII whitespace characters | `readStateFreq` picks `,` or `/` once and accepts a space (1766-1768, 1779) | related | different set | decision 006 | no | (iv) |
| `digit` cpp:77 | ASCII digit | none | no counterpart | n/a | n/a | n/a | n/a |
| `decimal` cpp:80-108 | ASCII decimal grammar | `is_number` (tools.cpp:813-818), `convert_double` (383-393), both `strtod` | related | `strtod` accepts hex, nan, leading whitespace | decisions 006, 012 | no | (iv) |
| `str` cpp:110-116 | 17 digits, classic locale | `convertDoubleToString` (tools.cpp:204-208, precision 6); the class's `number17` (C:37-42, no `imbue`) | related | precision, locale | decision 012 for tools | inside the fork, one function | (iv) for tools; merging with `number17` same bytes unless a non-classic global locale is installed (no `setlocale` found) |
| `STEP_ETA`, `STEP_SCALE` h:47-48 | 1e-4, 1 | `ERROR_X` (OPT:23), internal linkage | reimplements (constant) | none | decision 009 | no | (iii) |
| `FIRST_ORDER_ETA` h:112 | 1e-3 | none | no counterpart | n/a | n/a | n/a | n/a |
| `FIRST_ORDER_GTOL` h:113 | 1e-4 | `TOL_RATE` (modelmarkov.h:31); fits receive max(gradient_epsilon, TOL_RATE), default gradient_epsilon 1e-4 (modelfactory.h:205, 214; partitionmodel.h:91) | reimplements (value) | equal only when gradient_epsilon ≤ 1e-4 | decision 012 for the chart; the class could pass `TOL_RATE` | yes, from the class | (v) |

### Section B: conventions

- **Rate layout.** Dense n×n Eigen matrices with (i,j) as the rate from i to j, IQ-TREE's meaning. The
  coordinate order is IQ-TREE's packing with the reference also skipped. The class converts, and
  repeats loops IQ-TREE has:
  - packing Q into `rates[]` (C:530-533) equals `setFullRateMatrix`'s non-reversible branch
    (868-871): (i);
  - unpacking R (229-231, 248-250) repeats 114-118, 946-952, 1613-1618;
  - unpacking Q0 (253-260) repeats 1256-1264, with `Q0.row(i).sum()`'s order unverified. After the
    non-reversible read, `rate_matrix` already holds the normalized Q0 (1838, 1873): (ii);
  - reading `rate_matrix` uses a zero-copy `Map` (498).
- **State order.** The class's `STATES` (C:33) duplicates `symbols_protein` (alignment.cpp:34)
  without X.
- **Normalization.** Same definition as ModelMarkov's (Σπ*_i(−q_ii) = 1), reached implicitly through
  Σν = 1; the base then multiplies by delta ≈ 1 (MM:1269-1282). The only explicit copy is
  `reversibleSeed`'s division, which does nothing in production, because `inverse` and `references`
  are scale-invariant in exact arithmetic; the seed fixtures test it.
- **Stationary solve.** Same formulation, solver and right-hand side as IQ-TREE, on K − I.
  IQ-TREE's own π solve is skipped for NQC by `FREQ_USER_DEFINED` (1252-1254, 1266-1267), as
  constraint 3 requires.
- **Bounds.** IQ-TREE's mechanism is reused (`setBounds`, `fixBound` in `lnsrch` OPT:684,
  `bound_check` false). The box differs: [log 1e-5, log 1e5] on log-ratios, against
  [`MIN_RATE`, `MAX_RATE`] = [1e-4, 100] on raw rates for NONREV and GTR20 (MM:1139-1164).
- **Finite-difference step.** Its own rule (D02), on IQ-TREE's displacement idiom.

### Section C: the reverse direction

1. **`computeStateFreqFromQMatrix` as a header helper.** It edits `modelmarkov.cpp`; upstream numbers
   are unchanged if moved verbatim (derived); saves about 8 lines; needs decision 015 revisited.
2. **A diagonal fill helper from `decomposeRateMatrixNonrev`.** Bitwise; rejected by decision 015.
3. **`EigenDecomposition::computeRateMatrix` as a free function.** It edits
   `eigendecomposition.{h,cpp}` (overridden by PoMo, `modelpomo.h:321`); upstream bitwise if the
   member forwards; NQC changes, (ii), unless the chart keeps its division.
4. **`derivativeFunk`'s step and `dfpmin`'s test as helpers.** Forbidden by decision 009 and decision
   005's "Not chosen". The class's `scaledDerivative` (C:327-355) is a statement-level copy of
   `derivativeFunk` (OPT:916-939) with the step swapped: the largest copy of upstream code in NQC's
   numerical path.
5. **Without upstream edits.** `setFullRateMatrix` for packing, (i). In `readSeed`, feed the
   R_ij π*_j that `setReversible(false)` puts in `rates[]` to `transferSeed`, which equals
   `reversibleSeed` in exact arithmetic, (ii).
6. **Nothing to factor** in the matrix exponential, eigensystem, P(t) checks or likelihood.

### Section D: summary

- **By relation:** 6 reimplement, 14 related, 6 no counterpart (26 rows).
- **By effect:** (i) 0, (ii) 1, (iii) 6, (iv) 12, (v) 1, n/a 6.

Largest opportunities:

1. **`reversibleSeed`:** about 23 lines, (ii).
2. **Chart-internal repeats** (K three times, the exit check twice, dead `freeDestinations`): about
   19 lines, (i).
3. **`solveStationary`:** about 10 lines, (iii), bitwise.
4. **`completeDiagonal`:** about 8 lines, (iii), bitwise.
5. **The `firstOrderTest` and `scaledStep` formulas:** about 9 lines, (iii), forbidden by decision
   009.

Doc errors found:

- `fixedeqchart.h:106` "exactly h0 when abs(x) ≥ h0" contradicts CODE_PLAN 2.1 (lines 90-94) and
  the test at `test_fixedeqchart.cpp:424-426, 441-442` (a 2.2e-17 difference at x = 2);
- `fixedeqchart.h:110-113` "the gtol NQC's fits receive" holds only when gradient_epsilon ≤ 1e-4.

Not verified:

- bitwise equality of the `computeStateFreqFromQMatrix` route and of the pure moves;
- that its `ASSERT` never fires for a finite, positive K;
- that division and reciprocal multiplication differ on actual LG-at-π* entries;
- Eigen's `row().sum()` order;
- that every caller passes gradient_epsilon ≤ 1e-4.

Searched `model/`, `utils/`, `main/`, `tree/` and `alignment/`, plus targeted searches in
`simulator/` and `lbfgsb/`. Not searched: `nn/`, `terrace/`, `pll/` and vendored directories.

---

## Appendix C: user-facing usage

**Answer.** A nQMaker user's commands carry over almost unchanged: `--model-joint NONREV+FO` becomes
`--model-joint NQC+F{π*}` (or `NQC+F<name>`). The departures fall into three groups:

1. **Required by π*:** a target is mandatory, a new name, the start moved into the π* family, one
   target across linked partitions, and `-optfromgiven`, `-fs` and `-ft` refused.
2. **Not required by π*, caused by NQC's own additions:**
   - a mandatory end-of-run export (two files, a summary block and a 722-evaluation diagnostic);
   - the `--runs`, `--root-test` and `-mstore` refusals that follow from it;
   - `--nqc-domain`.
3. **Scope or protection:**
   - scope refusals: unlinked partitions, ModelFinder, AliSim from the model string, mixtures, `+H`,
     `+T`;
   - protective refusals where NONREV meets an upstream defect: `--lie-markov`, `--eigen`, and a
     rate model inside `--model-joint`.

On `-S` with a folder: by the code, it should work exactly as for NONREV, but no NQC run has used it.

Abbreviations: NQC.cpp `model/modelnonrevfixedeq.cpp`; MF `model/modelfactory.cpp`; MM
`model/modelmarkov.cpp`; MP `model/modelprotein.cpp`; MIX `model/modelmixture.cpp`; PM
`model/partitionmodel.cpp`; SA `alignment/superalignment.cpp`; PA `main/phyloanalysis.cpp`; PT
`main/phylotesting.cpp`; T `utils/tools.cpp`; chk `test_scripts/fixedeq/nqc/nqcchecks.py`.

Test status: the last full `nqcchecks.py` gate was at `36a43711` (binary `25d78d6e`), 442 asserted
checks passed (CHANGELOG.md:337-339). The C++ has not changed since. CI runs only `--cases seed`.

### Table A: workflow steps and command forms

| # | Step or form | GTR20 / NONREV | NQC | Class | What differs | Reason, needed for π*? | Existing form? | Tests |
|---|---|---|---|---|---|---|---|---|
| 1 | general non-reversible model, one alignment | `-m NONREV[+G4]`, LG seed, FREQ_ESTIMATE, 379 parameters (MP:1213-1239); π re-solved each decomposition (MM:1252-1268) | `-m NQC+F{π*}[+G4]` (MIX:3231-3235; NQC.cpp:117-164) | differs | new name; target mandatory; 360 parameters | decision 002 (name provisional, maintainer question 1); constraint 9. Yes | no (constraints 4, 9) | chk `seed-*`, `fit-*`, `root-*` |
| 2 | `-m NONREV+F{π}` | π pinned (MP:1246-1250), solve skipped (MM:1252, 1266); Q unconstrained | `NQC+F{π*}`: π*Q = 0 by construction | differs | same syntax, different model | constraint 9 | – | `probes.py` (a), (b), unmodified binary |
| 3 | `-m GTR20+F{π*}` | 189 parameters | `NQC+F{π*}` | same | same `+F{...}` route; 0 frequency degrees of freedom | 006 | – | chk `nesting-*` |
| 4 | joint step `--model-joint NONREV+FO` / `GTR20+FO` | MF:284-305 | `--model-joint NQC+F{π*}` | differs | `+FO` → `+F{π*}` or `+F<name>` | 006; required | no | chk `linked-*` |
| 5 | `+F`, `+FO`, `+FQ`, bare `+FU`, none | accepted; π solved from Q; `+F` books 19 frequency degrees of freedom (by reading); bare `+FU` pins LG's π (by reading) | refused (NQC.cpp:180-182) | refused | message points to NONREV; exit 2 | 006, 008; required | no | chk `refuse-freq-*` |
| 6 | literal `+F{...}` | IQ-TREE's reader (MM:1761-1797): distribution tokens, zeros, renormalised with a warning | NQC's parser: 20 entries, finite, > 0 and ≥ max(min_state_freq, ZERO_FREQ), sum within 1e-6, normalised once | differs | stricter acceptance | 006, D07; the 1e-6 line is a choice | syntax reused | chk `target-sum`, `target-minimum` |
| 7 | `+FO{...}`, `+FU{...}` | rewritten to `+F{...}` (MF:483-485) | same | same | – | – | – | none |
| 8 | `+F<name>` | models-block lookup (MF:599-604) | same; built-in vectors such as `C10pi1` also resolve (by reading) | same | – | – | – | chk `seed-aa-G4-named` |
| 9 | `NAME{params}` | `NONREV{379 values}` read as a fixed matrix (MP:1140-1183) | refused (NQC.cpp:178-179) | refused | – | 008, scope | `-m FILE` does this job | chk `refuse-braces` |
| 10 | fixed matrix by name or file | MP:1122-1138, 1240-1244; full matrix read as non-reversible, π re-solved (MM:1798-1840) | the fit re-used as `-m <prefix>.NQC.qmat` | same | – | 010 | – | chk `export-reimport*`, `linked-export-reimport-*`, `short-reimport` |
| 11 | `-s` | – | same | same | – | – | – | chk `fit-*`, `root-*` |
| 12 | `-S DIR` | T:2107-2116; SA:91-96, 689-722 | same path | same, by reading | rate model only from `-m`, as for NONREV | – | – | none found |
| 13 | `-S FILE` | SA:97-106, 586-588, 637-645 | same | same | – | – | – | chk `linked-*` (all `-s concat -S charsets` with position ranges) |
| 14 | `-p FILE` or `-p DIR` | T:2083-2090 | same | same | – | 043 | – | chk `-p` Level 2, `options-link-alpha`; DIR: none |
| 15 | `-q` | T:2092-2100 | same; Level 1 uses `-q` with `{x}` rates | same | – | 024 | – | chk `linked-fit-q-level1`; `sim/pipeline.py:400-405` |
| 16 | `-Q` | T:2076-2081 | same | same | – | – | – | chk `options-edge-unlinked` |
| 17 | `--link-model` | T:5030-5033; links by `getName` (PM:79-83) | same, plus the target check | same | – | – | – | chk `linked-fit-link-model` |
| 18 | partition file without linking | one NONREV per partition (by reading) | refused (NQC.cpp:203-204) | refused | the message names `-p, -q or -S` but `-Q` is also refused | 008, scope; not required by π* | not without code | chk `refuse-unlinked-*` |
| 19 | `--init-model NAME` | NONREV: read, rescale, `setReversible(false)` at the start's π (MP:1220-1230); GTR20 rejects a non-reversible start | reversible start rebuilt at π*; non-reversible start's jump chain transferred; no rescale (NQC.cpp:237-265) | differs | same syntax, different conversion | D05, 002; the start must lie in the π* family | syntax reused | none found (by name) |
| 20 | `--init-model FILE` | as 19 | as 19 | differs | as 19 | as 19 | – | chk `nesting-*`, `gradient-nearzero`, `linked-gradient-nearzero`, `warm-*`; `sim2.py:523-533` |
| 21 | default start | LG rescaled and converted at LG's π (MP:1232-1237) | LG rebuilt at π* (NQC.cpp:225-235) | differs | internal | D05; required | – | chk `seed-*` |
| 22 | `-te`, rooted or unrooted | roots an unrooted tree, keeps a rooted one (MM:145-149) | same rule (NQC.cpp:267-276) | same | – | 016 | – | chk `root-*` |
| 23 | `-n 0` | T:2848-2858 | not refused | same, by reading | – | – | – | none |
| 24 | tree search | ordinary | ordinary; not a validation target | same | – | – | – | chk `options-ufboot` |
| 25 | rate heterogeneity in `-m` | – | same | same | – | – | – | chk `seed-*`, `fit-level2-*`, linked Level 2 |
| 26 | rate heterogeneity inside `--model-joint` | NONREV: `+G4{0.5}` silently dropped; bare `+G4` stops with "File not found NONREV+G4" (UPSTREAM_FINDINGS Q13) | refused (NQC.cpp:205-206) | refused | – | 008; protective | rates via `-m` or charsets accepted | chk `linked-refuse-joint-rate-*` |
| 27 | rate heterogeneity from charset strings (`.best_model.nex`) | PM:67-70; under `--model-joint` the charset's `+F` is discarded, its rate kept (MF:258-305) | same | same | – | – | – | chk `linked-fit-q-level1`, `linked-fit-link-model` |
| 28 | ModelFinder involving NQC | NONREV: no refusal; NQ.* candidates need `--nonrev-model` (PT:479-485) | refused once `dating_mf` is set (PT:1389; NQC.cpp:174-175) | refused | "give it with -m" misleads a `-m MFP --model-joint` user (by reading) | 008; deferred to S6 | ModelFinder on ordinary models unchanged | chk `refuse-modelfinder*`, `linked-refuse-modelfinder`; `-madd` and the `--model-joint` case: none |
| 29 | `--link-alpha` | T:5025-5028 | same | same | – | – | – | chk `options-link-alpha` |
| 30 | `-T N`, `-T AUTO` | – | same; non-bitwise threads as NONREV (Q12) | same | – | 038 | – | chk `options-threads-auto`, `threads`, `export-threads` |
| 31 | checkpoint and resume | 10-digit rates (MP:1262-1282); a changed command line only warns (`main/main.cpp:2476-2494`) | adds 17-digit struct; restore refuses a different target, references, chart or domain (NQC.cpp:578-628) | differs | stricter; exit 2 | 003, constraint 6; holds | – | chk `restart-*` |
| 32 | `-redo` | T:5050-5055 | same | same | – | – | – | none |
| 33 | `--undo` | T:5062-5065 | same | same | – | – | – | chk `restart-undo`, `export-undo`, `restart-linked-undo` |
| 34 | `--runs N > 1` | one alignment: not established; `-S` crashes (Q11) | refused (NQC.cpp:216-217) | refused | suggests one run per `-seed` | 008: the restored tree cannot compute the likelihood the export diagnostic needs; not required by π* | only by readying the tree or skipping the diagnostic | chk `refuse-runs*`, `linked-refuse-runs-*`; `runs-one` accepted |
| 35 | `-b` | – | same; `-S` refuses `-b` for every model (T:5596-5597) | same | – | – | – | chk `short-bootstrap` |
| 36 | `-B` | – | same | same | – | – | – | chk `options-ufboot` |
| 37 | `-blfix` | T:3401-3408 | same | same | – | 016, 024 | – | chk `fit-level1` |
| 38 | `--root-test` | accepted; re-fits branches even with `-blfix`; report mixes states (Q16) | refused (NQC.cpp:212-213) | refused | – | 008; not required by π* | – | chk `refuse-root-test` |
| 39 | `-mstore` | accepted; possible stale cache (Q9, unverified) | refused (NQC.cpp:208-209) | refused | – | 008; partly protective | – | chk `refuse-mstore` |
| 40 | `--lie-markov`, `--eigen` | stop at the `ASSERT` in `computeTransMatrixNonrev` (MM:501-503); Lie-Markov recorded (Q14), eigen by reading | refused (NQC.cpp:195-199) | refused | – | 008; protective | – | chk `refuse-lie-markov`, `refuse-eigen` |
| 41 | `--scaling-squaring` | MM:464-466 | as NONREV | same, by reading | – | – | – | none |
| 42 | `-jointopt` | accepted, `MIN_RATE`/`MAX_RATE` bounds (MF:1347-1352) | refused (NQC.cpp:185-186) | refused | – | 008: the bounds would clamp log-ratio coordinates; required by the chart | – | chk `refuse-jointopt` |
| 43 | `-optfromgiven` | forces the π solve (MM:1252, 1266) | refused (NQC.cpp:183-184) | refused | – | 008, 015; required | – | chk `refuse-optfromgiven` |
| 44 | `--show-lh` | T:3440-3450 | same | same | – | 025 | – | chk `seed-*`, `export-reimport*` |
| 45 | `-o` outgroup | – | same rooting rule | same, by reading | – | 016 | – | none |
| 46 | `--loop-model`, `-nparam`, `-me` | – | same | same | – | 039 | – | used across chk; `sim/pipeline.py:405` |
| 47 | AliSim with the model string | accepted; simulates from the start, not a fit (by reading) | refused (NQC.cpp:200-201) | refused | – | 008, scope | row 48 does the job | chk `refuse-alisim` |
| 48 | `--alisim -m FILE` | same reader | same | same | – | – | – | format exercised (`sim/alisim_check.py:159`, `sim/pipeline.py:197`, `nqc/toy.py:53`); an actual `.NQC.qmat` into AliSim: none |
| 49 | `.iqtree` model line | `NONREV+FO+G4` (MM:241-246, 201-239) | `NQC+G4` | same | target not shown, as for `NONREV+F{}` and `GTR20+F{}` | maintainer question 8 | – | not asserted |
| 50 | `.iqtree` full Q block | PA:632-666, 6 decimals | printed for NQC; last row π* (NQC.cpp:539-545) | same | – | 010 adds `.qmat` | – | `regress.py:62` for NONREV only |
| 51 | `.iqtree` state frequencies | NONREV+FO: "(estimated with maximum likelihood)" plus pi lines (PA:684-713) | "(model)", no pi lines | same mechanism | – | – | – | none |
| 52 | `.best_model.nex` | `tree/phylosupertree.cpp:1526-1565`; NONREV's string carries π but not Q | `NQC+G4{a}`, no target | same mechanism | re-use without `--model-joint` refused ("needs its target") | maintainer question 8 | – | none |
| 53 | `.treefile` | – | same (rooted) | same | – | – | – | chk `root-*` |
| 54 | `.iqtree` ALISIM COMMAND | PA:2094-2187; NONREV's printed command simulates from a fresh LG-based NONREV (by reading) | prints `-m "NQC+G4{a}"` or `-p best_model.nex`, which is refused with "needs its target" (frequency check before the AliSim check, NQC.cpp:180 vs 200; by reading) | differs | the printed command fails | maintainer question 8 | pointing to `-m <prefix>.NQC.qmat` would need a shared change | none |
| 55 | `<prefix>.NQC.qmat` | – | NQC.cpp:717-731 | new | existing `-m FILE` format | 010; not required by π* | overlaps row 50 | chk export rows |
| 56 | `<prefix>.NQC.info` | – | NQC.cpp:829-898 | new | – | 009, 010; not required by π* | partly duplicates the checkpoint | chk `.info` rows |
| 57 | end-of-run summary and export diagnostic | NONREV prints nothing (MM:368-376) | NQC.cpp:663-827; 722 extra evaluations | new | – | 009, 010; not required by π* | – | chk export and diagnostic rows |
| 58 | `--nqc-domain LO,HI` | none (compile-time box) | T:3167-3178, 6016; NQC.cpp:140-141 | new | – | 007, D03; not required by π* | not without a rebuild per domain | chk `fits`, `nesting-*`, `risk22-*`, restart |
| 59 | `-v`/`-vv` prints and the fall warning | none | NQC.cpp:160-163, 344-351, 479-488 | new | – | 009, 021 | – | `fitrecords.py` |
| 60 | "might be overfitting" warning | NONREV, when not linked (MP:1214-1217) | absent | differs | cosmetic | none stated | yes | – |
| 61 | "Linking ..." and "Mean state frequencies" lines | printed only for FREQ_ESTIMATE/EMPIRICAL (PM:119-160) | not printed | same mechanism | – | – | – | – |
| 62 | help text | NONREV listed; `--model-joint`, `--link-model`, `--init-model` absent (upstream) | two added lines (T:6013-6016) | new | – | – | – | none |
| 63 | citation | nQMaker block when `model_joint` contains "NONREV" (PA:165-183) | same block when NQC is the joint's first component | same | – | Peter, 2026-10-26 | – | chk L1866 (reported) |
| 64 | parameter count | 379, once when linked (PM:251-255) | 360, 0 frequency degrees of freedom | same mechanism | – | constraint 5 | – | chk count rows |
| 65 | linked partitions with different targets | both link as "NONREV", no check (by reading) | refused at the first linked fit (NQC.cpp:414-442) | refused | – | required | – | chk `linked-compat-*` |
| 66 | `MIX{...}`, `+FMIX`, `--mdef` mixture | not established | refused (MIX:3453, 3527; NQC.cpp:112-115) | refused | – | `+FMIX` required; `MIX` scope | – | chk `refuse-mix/-fmix/-mdef-mix` |
| 67 | `+H`, `+T`, `-hmmster` | not established | refused (NQC.cpp:187-192) | refused | – | scope | – | chk `refuse-ghost/-tree-mixture/-hmm`, `linked-refuse-ghost` |
| 68 | `-fs`, `-ft` | not established | refused (NQC.cpp:176-177) | refused | – | required | – | chk `refuse-site-freq*` |
| 69 | `--link-exchange-rates` | IQ-TREE refuses it without GTR (T:5600-5601) | refused (NQC.cpp:193-194) | refused | – | concept does not apply | – | chk `refuse-link-exchange` |
| 70 | non-protein data | not established | refused (MIX:3232-3233) | refused | – | inherent | – | chk `refuse-dna` |

### Section B: the `-S` folder question

- **Published workflow** (tutorial on the iqtree.github.io mirror, updated May 26, 2026):
  - step 1 runs `iqtree3 -S train_plant -mset LG,WAG,JTT -cmax 4 -pre train_plant`;
  - the joint fit runs `iqtree3 -S train_plant.best_model.nex -te train_plant.treefile --model-joint
    NONREV+FO --init-model LG -pre train_plant.NONREV`;
  - the deck's note (`nq_frequency_replacement_part2.tex:208`) gives the same pattern for the paper.
- **Code path.**
  - `parseArg` sets `partition_file` and `TOPO_UNLINKED` (T:2107-2116).
  - `readFromParams` sees a directory and calls `readPartitionFiles` (SA:91-96), which builds one
    alignment per file with `params.model_name` (SA:689-722).
  - A partition file goes through `readPartitionNexus` (SA:586-588, 637-645).
  - From there both take one path:
    - the `PartitionModel` constructor (PM:66-90);
    - LG as the default for an empty model string under `--model-joint` (MF:167-178);
    - `--model-joint` replacing the substitution and frequency strings and keeping the rate string
      (MF:284-305);
    - no ModelFinder for empty models under `--model-joint` (PT:1380-1386);
    - the createModel dispatch (MIX:3231-3235);
    - `link_model` (PM:79-83).
- **NQC's path tells them apart nowhere.**
  - `refuseUnsupported` tests `partition_file`, `model_joint` and `link_model` (NQC.cpp:203-206).
  - `startLinkedFit` uses the partition name only in its message.
  - `reportRun` and the hook select partitions by model name and `part_order`.
- **Consequences, the same as for NONREV.**
  - A folder with no `-m` fits with no rate heterogeneity.
  - With `-m LG+G4` every partition gets a free `+G4`.
  - Per-gene rate models come only through the `.best_model.nex` form.
  - `-S DIR -m NQC+F{π*}` without linking is refused.
- **Tests.** No test, driver case or recorded run uses a folder with NQC, or a partition file whose
  charsets name separate files. Every NQC `-S` run is `-s <concatenated> -S <charsets>`:
  - chk L493, L1424-1425, L1481-1486, L1499-1501, L1516-1534, L1627, L1797-1807, L1862, L2043;
  - `profile/profiling.py:350-356`;
  - `reference/reference.py:395`.

  Upstream's baseline run09 and probes a4, a5, c2, f8 and f9 use the same layout with NONREV.
- **Missed by the search:** unscripted runs under `~/iqtree3-runs`, `~/iqtree3-audit` or
  `~/iqtree3-scratch`, Puma runs, and commit messages without the search terms.

### Section C: new user-facing surface

1. **Model name `NQC`** (NQC.cpp:28, 153-154). Needed: constraint 9; final name is maintainer
   question 1.
2. **`--nqc-domain LO,HI`** (T:3167-3178). Not needed for π*; it serves D03/007 and decision 025.
3. **`<prefix>.NQC.qmat`**, 17 digits in the `-m FILE` format. Existing: the `.iqtree` full Q block
   at 6 decimals; `.GTRPMIX.nex` is a precedent. Not needed for π*.
4. **`<prefix>.NQC.info`**. Existing: the checkpoint, `.iqtree`, `.log`. Not needed.
5. **End-of-run block and the 722-evaluation diagnostic**, 6 to 14% of a first linked fit. Not needed
   for π*; decision 009's reason applies to NONREV (maintainer question 3).
6. **Verbose prints.** Not needed for π*.
7. **Default-verbosity warnings:** the fall below the start (not needed for π*), and the build-failure
   warning (specific to NQC's chart).
8. **Linked-compatibility refusal.** Needed.
9. **Incompatible-restart refusal.** Needed.
10. **Strict target parsing.** Needed in part.
11. **The `--init-model` conversion.** Needed.
12. **Two help lines and the citation routing.** Both reuse existing mechanisms.
13. **Twenty refusal and validation messages.**

### Section D: refusals

All exit through `outError` (status 2, "ERROR: " prefix), at construction, in 0.1 s or less with no
export (CHANGELOG.md:372-376). Messages are as in `refuseUnsupported` (NQC.cpp:166-220), plus
"NQC: invalid target +F{...}: <reason>" (136-138), the mixture refusal (112-115), the
linked-compatibility refusal (438-439) and the restore refusal (628). What NONREV does with each
combination, and whether each refusal is necessary or protective, is given in Table A, rows 5, 6, 9,
18, 26, 28, 34, 38 to 40, 42, 43, 47 and 65 to 70.

### Section E: summary

**Counts over the 70 rows:** 36 same (11 without an NQC test: rows 7, 12, 23, 32, 41, 45, 49, 50, 51,
52, 61), 10 differ, 18 refused, 6 new.

Departures ranked by how far they move NQC from nQMaker:

1. **A target is mandatory.** Needed.
2. **Unlinked partitions refused.** Scope.
3. **Mandatory export**, and the three refusals it drives. Not needed.
4. **ModelFinder refused.** Not needed; deferred to S6.
5. **AliSim from the model string refused, and the printed AliSim command fails.** Not needed.
6. **`--runs`, `--root-test` and `-mstore` refused.** From the export; partly protective.
7. **`--nqc-domain`.** Needed for D03.
8. **The `--init-model` and default-start conversions.** Needed.
9. **Strict parsing, the restart refusal and the linked-target refusal.** Needed.
10. **Protective refusals.**
11. **Cosmetic differences.**

No-loss adoptions:

- NONREV's construction warning;
- the unlinked refusal naming `-Q`;
- the ModelFinder message for `--model-joint` users;
- the AliSim refusal naming the `.qmat`;
- the help text's "<name> from --mdef";
- optionally, gating the end-of-run summary at `-v`, which would reverse Peter's S3-1 choice.

Not verified at runtime:

- `-S DIR`, `-p DIR` and the per-file charset form;
- `--init-model NAME`, `-n 0`, `-redo`, `-o`, `--scaling-squaring` and `-madd`;
- built-in frequency names as targets;
- the order of failure of the printed AliSim command;
- NONREV's behaviour where marked "by reading" or "not established".

Test outcomes came from CHANGELOG gate summaries. The nQMaker paper was not read; its commands come
from the deck's note.
