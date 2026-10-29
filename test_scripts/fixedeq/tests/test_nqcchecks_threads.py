"""The driver's thread checks (decision 038): the export rows of fixed_point_against, the NaN-safe
seed_against and linked residual line, and the -p and -q thread cap (CODE_PLAN.md section 3.5)."""
import itertools
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nqc"))
import nqcchecks as nc  # noqa: E402

STATIONARY = nc.STATIONARY
# f = 100 and |x| <= 1: the export_lnl line is 1e-6, the gradient line about 1e-7, the label band about 1e-9
INFO = {"command": "iqtree3 -s aln.phy", "coordinates": "0.5,-1", "export_lnl": "-100",
        "stop_gradient": "0.001,0.002", "stop_statistic": "1.0000000000000001e-05", "stop": STATIONARY,
        "stop_reason": "the projected gradient passes the optimizer's test", "stop_coordinate": "2", "fits": "3"}


class FakeRun:
    """What the checks read from a Runner's run: its exit code, stdout and export files."""

    def __init__(self, root, name, code=0, info=INFO, qmat="-1 1\n1 -1\n0.5 0.5\n", stdout=""):
        self.name, self.code, self.dir = name, code, root / name
        self.dir.mkdir()
        (self.dir / "stdout.txt").write_text(stdout, encoding="utf-8")
        if qmat is not None:
            self.file(".NQC.qmat").write_text(qmat, encoding="utf-8")
        if info is not None:
            text = "# NQC export record\n" + "".join(f"{k}: {v}\n" for k, v in info.items())
            self.file(".NQC.info").write_text(text, encoding="utf-8")

    def file(self, suffix):
        return self.dir / f"{nc.PREFIX}{suffix}"

    def stdout(self):
        return (self.dir / "stdout.txt").read_text(encoding="utf-8")

    def report(self):
        return ""


def rows(tmp_path, *changes):
    """fixed_point_against for a -T 1 reference and one run per change, each a dict of keyword
    arguments for FakeRun; returns {row's label after the colon: passed}."""
    ref = FakeRun(tmp_path, "ref")
    runs = [FakeRun(tmp_path, f"run{k}", **change) for k, change in enumerate(changes)]
    return {row["check"].split(": ", 1)[1]: row["passed"] for row in nc.fixed_point_against("case", ref, runs)}


def info(**values):
    return {**INFO, **values}


def asserted(result):
    return {k: v for k, v in result.items() if v is not None}


def test_equal_runs_pass_every_asserted_row(tmp_path):
    result = rows(tmp_path, {}, {"info": info(command="iqtree3 -T 4")}, {})
    assert all(asserted(result).values()) and len(asserted(result)) == 6


def test_an_aborted_run_fails(tmp_path):
    result = rows(tmp_path, {}, {"code": -6, "info": None, "stdout": "Assertion `restored == d.f0' failed."})
    assert not result["every run finishes with both export files (decision 038)"]
    assert not any(asserted(result).values())


def test_a_missing_info_fails(tmp_path):
    assert not rows(tmp_path, {"info": None})["every run finishes with both export files (decision 038)"]


def test_a_nan_export_lnl_fails_even_after_a_good_run(tmp_path):
    result = rows(tmp_path, {}, {"info": info(export_lnl="nan")})
    assert not result["export_lnl within 1e-8 |lnL| of the reference's (decision 038)"]


def test_a_nan_gradient_fails_even_after_a_good_run(tmp_path):
    result = rows(tmp_path, {}, {"info": info(stop_gradient="0.001,nan")})
    assert not result["stop_gradient within the same-rule line of the reference's (decision 021)"]


def test_a_different_qmat_fails(tmp_path):
    assert not rows(tmp_path, {"qmat": "-1 1\n1 -1\n0.5 0.50000000000000011\n"})[
        ".qmat byte-identical to the -T 1 reference's (decision 021)"]


def test_another_info_field_fails(tmp_path):
    assert not rows(tmp_path, {"info": info(fits="4")})[
        ".info as the reference's apart from the command and the diagnostic (decision 021)"]


def test_export_lnl_against_the_line(tmp_path):
    key = "export_lnl within 1e-8 |lnL| of the reference's (decision 038)"
    inside, outside = tmp_path / "inside", tmp_path / "outside"
    inside.mkdir()
    outside.mkdir()
    assert rows(inside, {"info": info(export_lnl="-100.0000009")})[key]
    assert not rows(outside, {"info": info(export_lnl="-100.0000011")})[key]


def test_the_label_is_asserted_away_from_the_band(tmp_path):
    other = info(stop="unclassified", stop_reason="the projected gradient fails the optimizer's test")
    result = rows(tmp_path, {"info": other})
    assert result["stop label equal to the reference's"] is False


def test_the_label_is_reported_within_the_band(tmp_path):
    near = {**INFO, "stop_statistic": "0.0001"}
    ref = FakeRun(tmp_path, "ref", info=near)
    other = FakeRun(tmp_path, "run", info={**near, "stop": "unclassified", "stop_reason": "fails"})
    result = {row["check"].split(": ", 1)[1]: row["passed"] for row in nc.fixed_point_against("case", ref, [other])}
    assert result["stop label equal to the reference's (report: within the rounding band)"] is None
    assert result[".info as the reference's apart from the command and the diagnostic (decision 021)"]


def test_a_reference_without_its_export_stops_for_peter(tmp_path):
    ref = FakeRun(tmp_path, "ref", code=-6, info=None)
    out = nc.fixed_point_against("case", ref, [FakeRun(tmp_path, "run")])
    assert len(out) == 1 and out[0]["passed"] is False and "stop for Peter" in out[0]["check"]


def gradient_text(lnl, g):
    return (f"NQC first gradient at -lnL {lnl} (coordinate index, coordinate, step, derivative of -lnL):\n" +
            "".join(f"NQC gradient {k} 0.5 0.0001 {v}\n" for k, v in enumerate(g, 1)))


def seed_rows(tmp_path, texts):
    ref = FakeRun(tmp_path, "ref", stdout=gradient_text(100.0, [0.5, 0.25]))
    runs = [FakeRun(tmp_path, f"run{k}", stdout=t) for k, t in enumerate(texts)]
    return [row["passed"] for row in nc.seed_against("seed", ref, runs)][:2]


def test_seed_against_passes_equal_runs(tmp_path):
    assert seed_rows(tmp_path, [gradient_text(100.0, [0.5, 0.25])] * 2) == [True, True]


def test_seed_against_fails_a_nan_lnl_in_a_later_run(tmp_path):
    # the built-in max kept the first run's value and passed this
    ok = gradient_text(100.0, [0.5, 0.25])
    assert seed_rows(tmp_path, [ok, gradient_text("nan", [0.5, 0.25])])[0] is False


def test_seed_against_fails_a_nan_gradient_in_a_later_run(tmp_path):
    ok = gradient_text(100.0, [0.5, 0.25])
    assert seed_rows(tmp_path, [ok, gradient_text(100.0, [0.5, "nan"])])[1] is False


def test_beyond_counts_nan_as_outside_the_line():
    assert nc.beyond(math.nan, 1.0) and nc.beyond(1.0 + 1e-12, 1.0)
    assert not nc.beyond(1.0, 1.0) and not nc.beyond(0.5, 1.0)


class FakeTwoPart:
    def __init__(self, columns):
        """columns: one list of two-taxon columns per partition."""
        self.aln = {"a": "", "b": ""}
        self.charsets, start = {}, 0
        for k, cols in enumerate(columns):
            self.aln["a"] += "".join(c[0] for c in cols)
            self.aln["b"] += "".join(c[1] for c in cols)
            self.charsets[f"p{k}"] = np.arange(start, start + len(cols))
            start += len(cols)


def test_thread_cap_follows_superalignment():
    pairs = ["".join(p) for p in itertools.product("ARNDCQEGHILKMFPSTWYV", repeat=2)]   # 400 distinct columns
    cap, patterns = nc.thread_cap(FakeTwoPart([pairs, pairs[:399] + pairs[:5]]))
    # 400 x 20 // 4000 = 2; 399 distinct x 20 // 4000 = 1
    assert patterns == [400, 399] and cap == 3
