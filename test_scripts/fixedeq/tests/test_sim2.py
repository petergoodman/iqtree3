"""SIM-2's drivers (sim/sim2.py, sim/cloak.py and the pipeline additions; CODE_PLAN.md 3.8 and 4) on constructed
inputs and committed files only, so CI runs them without a binary or downloaded data."""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

FIXEDEQ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIXEDEQ))
sys.path.insert(0, str(FIXEDEQ / "sim"))
from oracle import builtin, chart, simulate  # noqa: E402
import cloak  # noqa: E402
import common  # noqa: E402
import pipeline as pl  # noqa: E402
import sim2  # noqa: E402

TREE = "((a:0.1,b:0.2):0.05,(c:0.3,(d:0.1,e:0.15):0.1):0.2);"


def small_dataset(path, seed=3, lengths=(40, 25)):
    """A replicate's data as simulate_into writes it, for five taxa and two genes of random residues."""
    rng = np.random.default_rng(seed)
    path.mkdir()
    order = pl.leaf_order(TREE)
    alns = [{t: "".join(rng.choice(list(pl.STATES), n)) for t in order} for n in lengths]
    genes = [{"length": n, "shape": 0.5, "rate": 1.0} for n in lengths]
    for k, aln in enumerate(alns, start=1):
        (path / f"gene{k:04d}").mkdir()
        pl.write_phylip(path / f"gene{k:04d}" / f"g{k:04d}.phy", aln, order)
    text, sets = pl.concatenate(alns, order)
    (path / "alignment.phy").write_text(text)
    (path / "sets.nex").write_text(sets)
    (path / "parts.nex").write_text(pl.q_file(sets, genes, [1.0, 1.0]))
    (path / "tree.nwk").write_text(TREE + "\n")
    (path / pl.TREE_RATE1).write_text(TREE + "\n")
    pi, truths = simulate.truths()
    pl.write_matrix(path / "truth.qmat", truths["NQ.bird"], pi)
    (path / "data.json").write_text(json.dumps({"genes": genes, "m": 1.0, "seed": 7, "replicate": 2,
                                                "target": common.target_text(), "taxa": 5, "files": {}}))
    return alns


# ---------------------------------------------------------------- pipeline additions

def test_profile_counts_sum_to_n_and_differ_from_the_genes_own_stream():
    w = builtin.frequency_profiles("C20")[1]
    counts = pl.component_counts(sim2.BASE, 1, 7, 523, w)
    assert sum(counts) == 523 and len(counts) == 20 and counts == pl.component_counts(sim2.BASE, 1, 7, 523, w)
    own = np.random.default_rng(np.random.SeedSequence([sim2.BASE, 1, 7])).multinomial(523, w).tolist()
    assert counts != own


def test_draw_is_prefix_stable():
    loci = pl.read_scaffold(sim2.SCAFFOLD)[1]
    assert pl.draw(loci, 30, sim2.BASE, 4) == pl.draw(loci, 100, sim2.BASE, 4)[:30]


def test_shift_truths_are_frequency_replacement_of_r0():
    pi = common.target_values()
    r0 = simulate.dial(simulate.truths()[1]["NQ.bird"], pi, 0.0)
    R = r0 / pi[None, :]
    np.fill_diagonal(R, 0.0)
    for name, tv in (("shift:LG:0.06", 0.06), ("shift:Q.MAMMAL:own", None)):
        Q, pi_d = pl.truth(name, pi)
        assert np.max(np.abs(Q - chart.reversible_seed(R, pi_d)) / np.abs(chart.reversible_seed(R, pi_d))) <= 1e-10
        if tv is not None:
            assert abs(simulate.total_variation(pi_d, pi) - tv) <= 1e-12
    assert not pl.at_target("shift:LG:0.06") and not pl.at_target(sim2.PROFILES) and pl.at_target("dial:0.1")


def test_profile_reference_is_reversible_and_balanced_at_pi():
    pi = common.target_values()
    comps, ref = pl.profile_components(pi)
    assert len(comps) == 20 and simulate.circulation(ref, pi) <= 1e-12
    assert np.max(np.abs(pi @ ref)) <= 1e-12


def test_risk28_split_is_zero_at_the_truth():
    pi, truths = simulate.truths()
    m = pl.rate_metrics(truths["NQ.bird"], truths["NQ.bird"], pi)
    assert m["log_error_rms"] == 0 and m["tenfold"] == 0 and m["log_error_rms_flux_weighted"] == 0
    assert m["pair_sum_log_error_by_quintile"] == [0.0] * 5 and m["direction_share_error_by_quintile"] == [0.0] * 5


def test_pruning_keeps_the_charsets(tmp_path):
    alns = small_dataset(tmp_path / "src")
    (tmp_path / "dst").mkdir()
    man, same = pl.prune_dataset(tmp_path / "src", tmp_path / "dst", {"a", "c", "e"})
    assert same and man["taxa"] == 3 and man["pruned"] == ["a", "c", "e"]
    assert sorted(pl.leaf_order((tmp_path / "dst" / pl.TREE_RATE1).read_text().strip())) == ["a", "c", "e"]
    gene = pl.gene_alignments_of(tmp_path / "dst", 2)[1]
    assert gene == {t: alns[1][t] for t in ("a", "c", "e")}


def test_masked_datasets_nest_across_strengths(tmp_path):
    small_dataset(tmp_path / "src", lengths=(300, 200))
    edges, p = np.array([0.0, 0.5, 1.0]), np.array([0.6, 0.2])
    masks = []
    for s in sim2.STRENGTHS:
        (tmp_path / f"s{s}").mkdir()
        man, same = pl.mask_dataset(tmp_path / "src", tmp_path / f"s{s}", s, edges, p)
        masks.append(np.array([[c == "-" for c in g[t]] for g in pl.gene_alignments_of(tmp_path / f"s{s}", 2)
                               for t in "abcde"], dtype=object))
        assert same and 0 < man["masking"]["masked_fraction"] < 1
    flat = [np.concatenate([np.array(row, dtype=bool) for row in m]) for m in masks]
    assert all(np.all(a <= b) for a, b in zip(flat, flat[1:])) and flat[0].sum() < flat[2].sum()


def test_gtr20_state_and_paml_carry():
    rates = [f"{0.5 + k / 100:.10g}" for k in range(190)]
    freqs = ", ".join(["0.05"] * 20)
    ck = {"PartitionModelPlen!GTR20!ModelProtein!rates": ", ".join(rates),
          "PartitionModelPlen!GTR20!ModelProtein!state_freq": freqs,
          "PartitionModelPlen!g0001!ModelProtein!state_freq": freqs}
    st = pl.gtr20_state(ck)
    assert np.array_equal(st["R"], st["R"].T) and st["R"][0, 1] == 0.5 and st["R"][18, 19] == float(rates[-1])
    assert list(st["part_freqs"]) == ["g0001"] and np.allclose(st["freq"], 0.05)
    rows = pl.gtr20_paml(st, common.target_text()).splitlines()
    assert [len(r.split()) for r in rows[:19]] == list(range(1, 20)) and len(rows[19].split()) == 20
    with pytest.raises(ValueError):
        pl.gtr20_state({k: v for k, v in ck.items() if "rates" not in k})
    counted = {k.replace("!GTR20!", "!GTR20+F!"): v for k, v in ck.items()}   # GTR20+F's linked name
    assert np.array_equal(pl.gtr20_state(counted)["R"], st["R"]) and list(pl.gtr20_state(counted)["part_freqs"]) == \
        ["g0001"]


def test_level2_state_reads_the_checkpoint_tree_by_alignment_order(tmp_path):
    ckp = ("---\nPhyloTree:\n newick: ((0:0.1000000001,1:0.2):0.05,(2:0.3,(3:0.1,4:0.15):0.1):0.2);\n"
           "PartitionModelPlen:\n part_rates: 0.5, 1.5\n g0001!RateGamma!gamma_shape: 0.7\n"
           " g0002!RateGamma!gamma_shape: 1.2\n")
    with gzip.open(tmp_path / f"{pl.PREFIX}.ckp.gz", "wt") as fh:
        fh.write(ckp)
    (tmp_path / f"{pl.PREFIX}.treefile").write_text("((a:0.1000000001,b:0.2):0.05,(c:0.3,(d:0.1,e:0.15):0.1):0.2);\n")
    st = pl.level2_state(tmp_path, list("abcde"))
    assert st["rates"] == [0.5, 1.5] and st["shapes"] == [0.7, 1.2]
    assert pl.leaf_depths(pl.readers.parse_newick(st["tree"]))["a"] == 0.05 + 0.1000000001
    with pytest.raises(ValueError):
        pl.level2_state(tmp_path, list("bacde"))


def test_two_cpus_takes_the_lowest_allowed():
    assert pl.two_cpus({9, 5, 3}) == "3,5"
    with pytest.raises(ValueError):
        pl.two_cpus({4})


def test_read_curve(tmp_path):
    path = tmp_path / "curve.tsv"
    path.write_text("# script: x\nbin_lo\tbin_hi\tresidues\tremoved\tp\n0\t0.5\t10\t4\t0.4\n0.5\t1\t10\t1\t0.1\n")
    edges, p = pl.read_curve(path)
    assert edges.tolist() == [0.0, 0.5, 1.0] and p.tolist() == [0.4, 0.1]


# ---------------------------------------------------------------- the design and the scheduler

def test_design_is_deterministic_and_has_the_planned_cells():
    des = sim2.design()
    assert sim2.design_body(des) == sim2.design_body(sim2.design())
    ids = [t["id"] for t in des["tasks"]]
    assert len(ids) == len(set(ids))
    nqc = [t for t in des["tasks"] if t["method"] == "nqc"]
    assert len(nqc) == 308 and len(des["tasks"]) - len(nqc) == 176
    held = {d["replicate"] for d in des["datasets"].values() if d["genes"] == sim2.HELDOUT_GENES}
    fitted = {des["datasets"][t["dataset"]]["replicate"] for t in des["tasks"]}
    assert held == {sim2.HELDOUT_REPLICATE} and sim2.HELDOUT_REPLICATE not in fitted


def test_main_design_is_round_1s():
    """Every change to sim2.py keeps the five shards' design: its body hashes as round 1's design.json (SIM-2
    commit 3d, CODE_PLAN 3.8)."""
    body = sim2.design_body(sim2.design())
    assert hashlib.sha256(body.encode()).hexdigest() == "f1556f0d38256ca1e01c88742403c1238270c9e8c96e05160e8c8a4b1c39fa41"


def test_smoke_design_has_unique_tasks():
    des = sim2.design(smoke=True)
    ids = [t["id"] for t in des["tasks"]]
    assert len(ids) == len(set(ids)) and {t["cell"] for t in des["tasks"]} == {"A", "B", "C", "D", "Dr", "E", "F",
                                                                               "G", "H"}


def test_rate_metrics_of_a_reversible_fit_at_its_own_frequencies():
    """A reversible matrix scored at its own stationary distribution has no net flux and equal direction
    shares, so its net-flux error is the truth's net flux and its share errors are the truth's departures."""
    pi, truths = simulate.truths()
    rev = chart.reversible_seed(builtin.exchangeabilities("LG")[0], pi * 0 + 0.05)
    m = pl.rate_metrics(rev, truths["Q.bird"], pi, pi_fit=pi * 0 + 0.05)
    assert m["net_flux_error_max"] <= 1e-15 and max(m["direction_share_error_by_quintile"]) <= 1e-15
    at_pi = pl.rate_metrics(rev, truths["Q.bird"], pi)     # its fluxes taken at pi* instead: not balanced there
    assert at_pi["net_flux_error_max"] > 1e-4


def test_read_curve_refuses_an_empty_bin(tmp_path):
    path = tmp_path / "curve.tsv"
    path.write_text("bin_lo\tbin_hi\tp\n0\t0.5\tnan\n0.5\t1\t0.1\n")
    with pytest.raises(ValueError):
        pl.read_curve(path)


def test_shards_balanced_with_d_and_f_together():
    des = sim2.design()
    tasks = {t["id"]: t for t in des["tasks"]}
    for t in des["tasks"]:   # the last shard holds the pinned cells and nothing else
        assert (des["assignment"][t["id"]] == sim2.SHARDS) == (t["cell"] in sim2.PINNED)
        if t["init_from"]:
            assert des["assignment"][t["init_from"]] == des["assignment"][t["id"]]
    load = [0.0] * sim2.SHARDS
    for tid, k in des["assignment"].items():
        load[k - 1] += sim2.task_core_hours(tasks[tid])
    free = [t for t in des["tasks"] if t["cell"] not in sim2.PINNED]
    assert max(load[:-1]) - min(load[:-1]) <= max(sim2.task_core_hours(t) for t in free)


def fake(phases):
    out = {pid: {"id": pid, "cores": c, "est_h": e, "deps": d, "tail": e} for pid, (c, e, d) in phases.items()}
    for p in out.values():   # tails for these short chains
        p["tail"] = p["est_h"] + sum(q["est_h"] for q in out.values() if p["id"] in q["deps"])
    return out


def test_scheduler_core_budget_deadline_and_blocking():
    ph = fake({"a": (4, 2.0, []), "b": (4, 1.0, []), "c": (4, 3.0, []), "d": (16, 1.0, []), "e": (1, 1.0, ["x"]),
               "x": (1, 1.0, [])})
    picked = [p["id"] for p in sim2.startable(ph, {}, {}, 8, 0.0, 100.0, 12)]
    assert picked == ["c", "a"]                      # longest tail first, 8 cores
    alone = fake({"d": (16, 1.0, [])})
    assert [p["id"] for p in sim2.startable(alone, {}, {}, 12, 0.0, 100.0, 12)] == ["d"]   # 16 cores take the 12
    assert sim2.startable(alone, {}, {}, 11, 0.0, 100.0, 12) == []
    late = sim2.startable(ph, {}, {}, 100, 10.0, 12.0, 100)
    assert {p["id"] for p in late} == {"b", "d", "x"}  # 10 + 1.5 x 2 and 10 + 1.5 x 3 pass the deadline 12
    state = {"x": "failed"}
    sim2.startable(ph, state, {}, 100, 0.0, 100.0, 100)
    assert state["e"] == "blocked"


def test_backfilling_reserves_for_the_first_wide_phase():
    """r holds 6 cores until 5, and w (8 cores, first in line) reserves time 5. With 8 cores in all nothing is
    spare then: short (ends at 4) starts, long (ends at 6) would delay w and waits. With 10 cores, 2 are spare at
    5, so long starts on them and short still fits before 5."""
    ph = fake({"w": (8, 9.0, []), "long": (2, 6.0, []), "short": (2, 4.0, []), "r": (6, 5.0, [])})
    ph["w"]["tail"], ph["long"]["tail"], ph["short"]["tail"] = 10.0, 6.0, 4.0
    assert [p["id"] for p in sim2.startable(ph, {}, {"r": 5.0}, 2, 0.0, 100.0, 8)] == ["short"]
    assert [p["id"] for p in sim2.startable(ph, {}, {"r": 5.0}, 4, 0.0, 100.0, 10)] == ["long", "short"]


def test_metrics_waits_for_b_but_runs_after_its_failure():
    ph = fake({"P": (1, 1.0, []), "B": (1, 1.0, ["P"]), "metrics": (1, 0.1, ["P"])})
    ph["metrics"]["after"] = ["B"]
    assert [p["id"] for p in sim2.startable(ph, {"P": "done"}, {}, 10, 0.0, 100.0, 10)] == ["B"]
    state = {"P": "done", "B": "failed"}
    assert [p["id"] for p in sim2.startable(ph, state, {}, 10, 0.0, 100.0, 10)] == ["metrics"]


def test_projection_respects_cores_and_the_deadline():
    """late holds a core from 0 to 30, so c (all 8 cores, after a) waits for it: makespan 31. With deadline 20,
    1.5 x 30 does not fit, late never starts, and c runs at 2 to 3."""
    ph = fake({"a": (4, 2.0, []), "b": (4, 2.0, []), "c": (8, 1.0, ["a"]), "late": (1, 30.0, [])})
    r = sim2.project(ph, 8, 50.0)
    assert r["makespan_h"] == 31.0 and r["phases_left"] == 0 and r["charged_core_h"] == 8 * 31.0
    r = sim2.project(ph, 8, 20.0)
    assert r["phases_left"] == 1 and r["makespan_h"] == 3.0


def test_scan_renames_leftover_partials_and_seeds_the_ratios(tmp_path):
    p = {"id": "T/B", "task": "T", "kind": "B", "size": 100, "est_h": 2.0}
    final = sim2.phase_dir(tmp_path, p)
    final.with_name("B.partial").mkdir(parents=True)
    state, interrupted, _ = sim2.scan(tmp_path, {"T/B": p})
    assert interrupted == ["T/B"] and state == {} and final.with_name("B.interrupted-1").is_dir()
    final.with_name("B.partial").mkdir()
    sim2.scan(tmp_path, {"T/B": p})
    final.mkdir()
    (final / "record.json").write_text(json.dumps({"status": "done", "record": {}, "seconds": 3 * 3600}))
    state, _, ratios = sim2.scan(tmp_path, {"T/B": p})
    assert final.with_name("B.interrupted-2").is_dir() and state == {"T/B": "done"} and ratios == {("B", 100): [1.5]}


def test_measured_counts_capped_fits_and_skips_crashes():
    p = {"est_h": 2.0}
    assert sim2.measured({"status": "failed", "record": {"at_cap": True}, "seconds": 7200}, p) == 1.0
    assert sim2.measured({"status": "failed", "record": {"at_cap": False}, "seconds": 7200}, p) is None
    assert sim2.measured({"status": "failed", "error": "trace", "seconds": 5}, p) is None
    assert sim2.measured({"status": "failed", "exit_code": -9, "phase": "x"}, p) is None


# ---------------------------------------------------------------- CLOAK's curve

def test_cloak_r_and_join(tmp_path):
    """Columns of the reference: (A, A, gap), (C, D, C), (E, E, E). r: the As 1 (one other residue present), the
    Cs 0.5, D 0, the Es 1."""
    (tmp_path / "c.tfa").write_text(">s1\nACE\n>s2\nADE\n>s3\nCE\n")
    (tmp_path / "c.msf").write_text(" Name: s1 oo Len: 3\n Name: s2 oo Len: 3\n Name: s3 oo Len: 3\n//\n"
                                    "s1 ACE\ns2 ADE\ns3 .CE\n")
    (tmp_path / "c.xml").write_text("<colsco-name>coreblock</colsco-name><colsco-data>1 0 1</colsco-data>")
    (tmp_path / "r.tsv").write_text("label\tpos\tresidue\tmsa0_col\tkept\tout_col\n"
                                    "s1\t1\tA\t1\t1\t1\ns1\t2\tC\t2\t0\t\ns1\t3\tE\t3\t1\t2\n"
                                    "s2\t1\tA\t1\t1\t1\ns2\t2\tD\t2\t0\t\ns2\t3\tE\t3\t1\t2\n"
                                    "s3\t1\tC\t2\t1\t3\ns3\t2\tE\t3\t1\t2\n")
    files = {"tfa": tmp_path / "c.tfa", "msf": tmp_path / "c.msf", "xml": tmp_path / "c.xml"}
    res, why, alone = cloak.case_residues(files, tmp_path / "r.tsv")
    assert why is None and alone == 0 and len(res) == 8
    assert sorted((round(r, 6), kept, core) for r, kept, core in res) == sorted(
        [(1, True, True), (0.5, False, False), (1, True, True), (1, True, True), (0, False, False), (1, True, True),
         (0.5, True, False), (1, True, True)])
    counts = cloak.bin_counts(res)
    assert counts[0].sum() == 8 and counts[1].sum() == 2 and counts[1][0] == 1
    total, p, lo, hi, core_p = cloak.curve_table([counts, counts], draws=50)
    assert total[0].sum() == 16 and p[0] == 1.0
    (tmp_path / "r2.tsv").write_text((tmp_path / "r.tsv").read_text().replace("s2\t2\tD", "s2\t2\tK"))
    assert cloak.case_residues(files, tmp_path / "r2.tsv")[0] is None
