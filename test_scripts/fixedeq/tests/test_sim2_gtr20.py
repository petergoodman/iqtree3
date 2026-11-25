"""SIM-2's seeded GTR20 arm (sim/sim2_gtr20.py, the seeded path of sim2.do_phase, pipeline.gtr20_seed and
rate_bounds, explore/gtr20_window.choose_p; decision 046, CODE_PLAN.md 3.8) on constructed inputs, without a
binary."""
import gzip
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

FIXEDEQ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FIXEDEQ))
sys.path.insert(0, str(FIXEDEQ / "sim"))
sys.path.insert(0, str(FIXEDEQ / "explore"))
import common  # noqa: E402
import pipeline as pl  # noqa: E402
import sim2  # noqa: E402
import sim2_gtr20 as g  # noqa: E402
import gtr20_window  # noqa: E402


def main_gtr20(smoke=False):
    return {t["id"]: t for t in sim2.design(smoke)["tasks"] if t["cell"] in ("D", "Dr") and t["method"] != "nqc"}


def test_design_is_cell_d_gtr20_seeded_in_three_one_wave_shards():
    des, main = g.design(), main_gtr20()
    assert sim2.design_body(des) == sim2.design_body(g.design())
    tasks = des["tasks"]
    assert len(tasks) == 176 and sum(t["cell"] == "D" for t in tasks) == 160 and sum(t["cell"] == "Dr" for t in tasks) == 16
    assert {t["method"] for t in tasks} == {"gtr20pi", "gtr20f"}
    assert sorted(t["id"].removesuffix("-seeded") for t in tasks) == sorted(main) and all(
        t["id"].endswith("-seeded") for t in tasks)
    full = sim2.design()["datasets"]
    for t in tasks:
        old = main[t["id"].removesuffix("-seeded")]
        assert t["gtr20_start"] == {"p": g.P} and (t["dataset"], t["heldout"]) == (old["dataset"], old["heldout"])
        assert {k: v for k, v in t.items() if k not in ("id", "gtr20_start", "phases")} == \
            {k: v for k, v in old.items() if k not in ("id", "phases")}
        fit = t["phases"][0]
        assert fit["cap_h"] == g.FIT_CAP_H and fit["est_h"] == g.EST_H[(t["method"], t["relaxed"])] * g.RATIO
        assert t["phases"][1] == old["phases"][1]
    assert all(des["datasets"][d] == full[d] for d in des["datasets"])
    by_shard = {k: [t for t in tasks if des["assignment"][t["id"]] == k] for k in (1, 2, 3)}
    assert len(by_shard[1]) == len(by_shard[2]) == 80 and len(by_shard[3]) == 16
    assert {full[t["dataset"]]["truth"] for t in by_shard[1]} == {sim2.TRUTHS["R0"]}
    assert {full[t["dataset"]]["truth"] for t in by_shard[2]} == {sim2.TRUTHS["NQ"]}
    assert {t["cell"] for t in by_shard[3]} == {"Dr"}
    for k, cores in g.SHARD_CORES.items():
        assert g.first_wave(des, k, cores)[0]
    assert not g.first_wave(des, 1, 60)[0]     # 80 fits at one core each cannot all start on 60


def test_first_wave_fails_above_the_ratio_the_deadline_allows(monkeypatch):
    """GTR20+F's 16.01 h estimate times 1.5 passes the 46 h deadline above a ratio of about 1.9."""
    monkeypatch.setattr(g, "RATIO", 1.85)
    assert g.first_wave(g.design(), 2, 84)[0]
    monkeypatch.setattr(g, "RATIO", 1.95)
    assert not g.first_wave(g.design(), 2, 84)[0]


def test_smoke_design_keeps_the_main_smoke_phases():
    des, main = g.design(smoke=True), main_gtr20(smoke=True)
    assert len(des["tasks"]) == 5 and des["shards"] == 1 and set(des["assignment"].values()) == {1}
    assert all(t["phases"] == main[t["id"].removesuffix("-seeded")]["phases"] for t in des["tasks"])


def test_gtr20_seed_fixes_tyr_val_at_p_and_scales_by_one():
    rows = [r.split() for r in pl.gtr20_seed(0.02, common.target_text()).splitlines()]
    assert [len(r) for r in rows[:19]] == list(range(1, 20)) and rows[19] == common.target_text().split(",")
    entries = [v for r in rows[:19] for v in r]
    assert (pl.STATES[19], pl.STATES[18]) == ("V", "Y") and rows[18][18] == "0.02"   # PAML's last entry: Val-Tyr
    assert entries.count("10") == 189 and max(float(v) for v in entries) == 10.0


def test_rate_bounds_counts_the_upper_and_the_floor():
    st = {"rates_text": ["1"] * 186 + ["0.0001", "0.00010001", "99.5", "0.05"]}
    b = pl.rate_bounds(st)
    assert b["fixed_text"] == "0.05" and b["largest_free"] == 99.5 and b["above_upper"] == 1 and b["at_floor"] == 2
    assert b["headroom"] == pytest.approx(100 / 99.5)


def fake_fit(rates, calls):
    """pipeline.fit_into's stand-in: records its init_model and writes a linked GTR20 checkpoint and stdout."""
    def fit_into(out, data, binary, threads, timeout, **kw):
        calls.append(kw.get("init_model"))
        d = Path(out) / "fit"
        d.mkdir(parents=True)
        freqs = ", ".join(["0.05"] * 20)
        with gzip.open(d / f"{pl.PREFIX}.ckp.gz", "wt") as fh:
            fh.write(f"---\nPartitionModelPlen:\n GTR20!ModelProtein!rates: {', '.join(rates)}\n"
                     f" GTR20!ModelProtein!state_freq: {freqs}\n")
        (d / "stdout.txt").write_text(f"{sim2.SEED_READ} seed.paml\n" * 3)
        return {"exit_code": 0, "at_cap": False, "assertion_or_abort": False, "data_match": True}
    return fit_into


def run_fit_phase(des, tid, tmp_path, monkeypatch, rates):
    calls = []
    monkeypatch.setattr(sim2.pl, "fit_into", fake_fit(rates, calls))
    t = next(t for t in des["tasks"] if t["id"] == tid)
    data = tmp_path / "data" / t["dataset"]
    data.mkdir(parents=True)
    (data / "data.json").write_text(json.dumps({"genes": [{}, {}, {}]}))
    phases = sim2.shard_phases(des, 1)
    out = tmp_path / "out"
    res, rows = sim2.do_phase(des, phases, phases[f"{tid}/fit"], tmp_path, out, SimpleNamespace(gate=str(tmp_path),
                                                                                               cores=12))
    return calls[0], res, rows


def test_only_seeded_fits_get_the_seed(tmp_path, monkeypatch):
    rates = ["10"] * 189 + [f"{g.P:.10g}"]
    init, res, rows = run_fit_phase(sim2.design(smoke=True), "D-gtr20pi-R0-g3-r1", tmp_path / "main", monkeypatch, rates)
    assert init is None and rows == [] and res["status"] == "done"
    init, res, rows = run_fit_phase(g.design(smoke=True), "D-gtr20pi-R0-g3-r1-seeded", tmp_path / "seeded", monkeypatch,
                                    rates)
    assert init == ("seed.paml", pl.gtr20_seed(g.P, common.target_text())) and res["status"] == "done"
    assert [r["passed"] for r in rows] == [True, None, None, None] and not rows[1]["flagged"]
    assert rows[3]["seed_reads"] == rows[3]["genes"] == 3


def test_a_bound_hit_is_flagged_and_a_wrong_fixed_rate_fails_with_the_phase_done(tmp_path, monkeypatch):
    """Peter, 2026-11-24: a rate above 99 is reported, flagged and counted; the fixed Tyr-Val stays asserted."""
    rates = ["10"] * 188 + ["99.5", "0.03"]
    _, res, rows = run_fit_phase(g.design(smoke=True), "D-gtr20f-R0-g3-r1-seeded", tmp_path, monkeypatch, rates)
    assert res["status"] == "done" and [r["passed"] for r in rows[:2]] == [False, None]   # fixed 0.03, not P
    assert rows[1]["flagged"] and rows[1]["above_upper"] == 1


def test_binary_check_refuses_other_binaries(tmp_path):
    (tmp_path / "a").write_bytes(b"not iqtree")
    (tmp_path / "b").write_bytes(b"nor this")
    assert [k for k, _ in g.binary_check(tmp_path / "a", tmp_path / "b")] == ["frozen", "gate"]


def test_choose_p_rule():
    assert gtr20_window.choose_p({0.02: -100.0, 0.05: -98.5}, {0.02: 10.0, 0.05: 4.0}) == 0.05
    assert gtr20_window.choose_p({0.02: -100.0, 0.05: -99.5}, {0.02: 10.0, 0.05: 9.0}) == 0.02    # gain not above 1
    assert gtr20_window.choose_p({0.02: -100.0, 0.05: -90.0}, {0.02: 10.0, 0.05: 3.9}) == 0.02    # headroom below 4
