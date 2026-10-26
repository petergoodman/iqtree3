import re
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "profile"))
import profiling as pf  # noqa: E402

COUNT_OUTPUT = """Breakpoint 1 at 0x66b7f0
IQ-TREE output
[Inferior 1 (process 134) exited normally]
PROFILE-EXITCODE
$1 = 0
PROFILE-EXITSIGNAL
$2 = void
PROFILE-BREAKPOINTS
Num     Type           Disp Enb Address            What
1       breakpoint     keep y   0x0000555555bbf7f0 <PartitionModel::targetFunk(double*)>
\tbreakpoint already hit 8343 times
\tignore next 1999991657 hits
2       breakpoint     keep y   0x0000555555aa1000 <ModelMarkov::targetFunk(double*)+4>
\tbreakpoint already hit 1 time
\tignore next 1999999999 hits
"""


def test_gdb_counts_exit_code_and_void_signal_are_read():
    g = pf.parse_gdb(COUNT_OUTPUT)
    assert g["exit_code"] == 0 and g["exit_signal"] is None
    assert g["hits"] == {pf.COUNTED: 8343, pf.PER_PARTITION: 1}


def test_an_abort_reads_as_a_signal_and_an_unhit_breakpoint_as_zero():
    text = ("PROFILE-EXITCODE\n$1 = void\nPROFILE-EXITSIGNAL\n$2 = 6\nPROFILE-BREAKPOINTS\n"
            "1       breakpoint     keep y   0x0000555555bbf7f0 <PartitionModel::targetFunk(double*)>\n")
    g = pf.parse_gdb(text)
    assert g["exit_code"] is None and g["exit_signal"] == 6 and g["hits"] == {pf.COUNTED: 0}


def test_a_paused_run_fails_the_cpu_guard_at_one_thread_only():
    text = "21.50 1000 20\n25.70 1400 30\n"   # 4.20 s of wall, 410 ticks of CPU
    assert pf.window_cpu(text) == pytest.approx(4.10)
    assert pf.window_cpu("21.50 246.78\n25.70 292.97\n") is None
    assert pf.cpu_row("w", 4.10, 4.20, 1)["passed"] is True
    assert pf.cpu_row("f", 42.9, 97.8, 1)["passed"] is False   # the paused development run
    assert pf.cpu_row("f", None, 97.8, 1)["passed"] is False
    assert pf.cpu_row("w", 1.0, 4.20, 4)["passed"] is None


def test_window_needs_both_readings_and_ignores_the_crossings_between():
    assert pf.window_seconds("21.50 246.78\n25.70 292.97\n") == pytest.approx(4.20)
    assert pf.window_seconds("21.50 246.78\n") is None
    script = pf.gdb_script("window")
    assert f"ignore 1 {pf.WINDOW - 1}" in script
    # inside the breakpoint's commands nothing may follow continue
    lines = [line.strip() for line in script.splitlines()]
    assert lines[lines.index("continue") + 1] == "end"


@pytest.mark.parametrize("mode, extra", [("window", ()), ("count", ()), ("count", (pf.PER_PARTITION,))])
def test_every_gdb_echo_ends_its_line_and_the_fit_breakpoint_is_numbered_last(mode, extra):
    lines = [line.strip() for line in pf.gdb_script(mode, extra).splitlines()]
    echoes = [line for line in lines if line.startswith("echo")]
    assert echoes and all(line.endswith("\\n") for line in echoes)
    if mode == "count":
        assert f"commands {len(extra) + 2}" in lines


def test_scaled_tree_keeps_topology_and_scales_every_length():
    text = "((T1:0.1000,T2:0.2000):0.0500,T3:0.3000);\n"
    scaled = pf.scale_tree(text, 2.5)
    assert scaled.replace(":", " ").split()[0] == "((T1"
    got = [float(v) for v in re.findall(r":([0-9.eE+-]+)", scaled)]
    assert got == pytest.approx([0.25, 0.5, 0.125, 0.75], rel=1e-15)


def test_partition_rates_are_seeded_with_mean_one():
    r = pf.partition_rates()
    assert len(r) == pf.COLLECTION and np.mean(r) == pytest.approx(1.0, rel=1e-14)
    assert np.array_equal(r, pf.partition_rates())


@pytest.mark.parametrize("n", pf.SETS)
def test_set_rates_survive_iqtree_renormalization_and_reproduce_simulated_lengths(n):
    r = pf.partition_rates()
    m, q = pf.set_rates(r, n)
    # equal-length partitions: normalizePartRate divides by the plain mean, which must already be 1
    assert np.mean(q) == pytest.approx(1.0, rel=1e-14)
    assert np.array(q) * m == pytest.approx(r[:n], rel=1e-14)


def test_charsets_tile_the_concatenation_and_sets_nest():
    text = pf.nexus_sets(3, sites=4)
    assert "charset p001 = 1-4;" in text and "charset p003 = 9-12;" in text
    alns = [{"T1": c * 4, "T2": c.lower() * 4, "T10": "-" * 4} for c in "ACDE"]
    whole = pf.concatenate(alns).splitlines()
    first = pf.concatenate(alns[:2]).splitlines()
    assert whole[0] == "3 16" and first[0] == "3 8"
    assert [line.split()[0] for line in whole[1:]] == ["T1", "T2", "T10"]
    assert all(w.split()[1].startswith(f.split()[1]) for w, f in zip(whole[1:], first[1:]))


@pytest.mark.parametrize("x, passed, verdict", [(3.5, False, "blocks"), (2.5, True, "investigate"),
                                                (1.0, True, "within")])
def test_per_fit_verdicts(x, passed, verdict):
    row = pf.per_fit_row("k", {"per_fit": x})
    assert row["passed"] is passed and row["verdict"].startswith(verdict)


def test_per_evaluation_verdict_only_at_one_thread():
    pe = {("10", "-S", 1, "nqc"): [0.13], ("10", "-S", 1, "nonrev"): [0.1],
          ("10", "-S", 4, "nqc"): [0.13], ("10", "-S", 4, "nonrev"): [0.1]}
    rows = pf.per_eval_rows(pe)
    by_threads = {r["check"].split(":")[0].endswith("-T 1"): r for r in rows}
    assert by_threads[True]["verdict"].startswith("investigate") and by_threads[False]["verdict"] is None
    assert by_threads[True]["ratio"] == pytest.approx(1.3)


def test_a_crossing_hidden_by_the_gdb_stop_cost_is_still_flagged():
    # raw 6 / 5 = 1.2 is within 1.25, but with 1.5 ms per stop removed from both it is 4.5 / 3.5
    pe = {("turtle", "-S", 1, "nqc"): [0.0060], ("turtle", "-S", 1, "nonrev"): [0.0050]}
    row = pf.per_eval_rows(pe, stop=0.0015)[0]
    assert row["ratio"] < pf.PER_EVAL_LINE < row["ratio_without_stop_cost"]
    assert row["verdict"].startswith("investigate")
    nqc, nonrev = {"complete": True, "wall_s": 29.0, "hits": 20000}, {"complete": True, "wall_s": 10.0, "hits": 1000}
    parts = pf.split(nqc, nonrev, 0.001, 0.001, stop=0.001)
    assert parts["per_fit"] == pytest.approx(2.9) and parts["per_fit_without_stop_cost"] == pytest.approx(1.0)
    # raw 10 / 4 = 2.5; with 2 s of stops removed from each, 8 / 2 = 4
    parts = pf.split({"complete": True, "wall_s": 10.0, "hits": 1000}, {"complete": True, "wall_s": 4.0, "hits": 1000},
                     0.001, 0.001, stop=0.002)
    row = pf.per_fit_row("k", parts)
    assert parts["per_fit"] < pf.BLOCK_LINE < parts["per_fit_without_stop_cost"] and row["passed"] is False


def test_per_fit_ratio_splits_into_per_evaluation_and_count_ratios():
    nqc, nonrev = {"complete": True, "wall_s": 30.0, "hits": 1000}, {"complete": True, "wall_s": 50.0, "hits": 2500}
    parts = pf.split(nqc, nonrev, 0.03, 0.02)
    assert parts["per_fit"] == pytest.approx(parts["per_evaluation"] * parts["evaluation_count"])
    assert parts["outside_evaluations_s"]["nqc"] == pytest.approx(0.0, abs=1e-12)
    assert parts["outside_evaluations_s"]["nonrev"] == pytest.approx(0.0, abs=1e-12)


def test_extrapolation_arithmetic():
    a, b = 0.002, 0.004   # seconds, seconds per alignment
    pe = {(str(n), "-S", 1, m): [a + b * n] for n in pf.SETS for m in ("nqc", "nonrev")}
    pe.update({("100", "-S", 4, m): [(a + 100 * b) / 2] for m in ("nqc", "nonrev")})   # efficiency 0.5
    # a long first linked fit, then short ones: the cap repeats the short ones, not the average
    fits = {"10-L1-S-G4": {"nqc": {"complete": True, "per_fit_evaluations": [9000, 700, 900]}},
            "10-L2-S-G4": {"nqc": {"complete": False, "per_fit_evaluations": [9000, 800]}}}
    got = pf.extrapolate(pe, fits, {10: 15, 30: 35, 100: 105})
    row = got["nqc -S"]
    e = a + 2000 * b
    n = 9000 + (pf.LOOP_DEFAULT - 1) * 800
    assert row["per_evaluation_s_million_columns"] == pytest.approx(e)
    assert row["efficiency_T4"] == pytest.approx(0.5)
    assert row["L1_at_cap"]["evaluations"] == pytest.approx(n)
    assert row["L1_at_cap"]["hours_32_cores"] == pytest.approx(n * e / 3600 / 16)
    assert "L2_at_cap" not in row   # an incomplete fit is not extrapolated
    assert got["memory_mb_million_columns"] == pytest.approx(5 + 2000)
    assert pf.linear_fit(pf.SETS, [a + b * n for n in pf.SETS])[2] == pytest.approx(0, abs=1e-12)
    assert pf.at_cap([9000], 10) is None


def test_evaluations_per_linked_fit_from_the_counts_at_each_start():
    text = ("x\nPROFILE-FIT\nNum     Type           Disp Enb Address            What\n"
            "1       breakpoint     keep y   0x1 <PartitionModel::targetFunk(double*)>\n"
            "IQ-TREE line\nhalf a linePROFILE-FIT\nNum     Type\n"
            "1       breakpoint     keep y   0x1 <PartitionModel::targetFunk(double*)>\n"
            "\tbreakpoint already hit 9000 times\n\tignore next 5 hits\nmore\n")
    g = pf.parse_gdb(text)
    assert g["fit_starts"] == [0, 9000]
    assert pf.evaluations_per_fit(g["fit_starts"], 9724) == [9000, 724]
    assert pf.evaluations_per_fit([], 9724) == [] and pf.evaluations_per_fit([0], None) == []


@pytest.mark.parametrize("mode, g, timed_out, want", [
    ("window", {"window_s": 4.2, "hits": {pf.COUNTED: 1001}}, False, 0),
    ("window", {"window_s": None, "hits": {pf.COUNTED: 1001}}, False, "window incomplete"),
    ("window", {"window_s": 4.2, "hits": {pf.COUNTED: 600}}, False, "window incomplete"),
    ("count", {"exit_code": None, "exit_signal": 6, "hits": {}}, False, "signal 6"),
    ("count", {"exit_code": None, "exit_signal": None, "hits": {}}, False, "no exit code"),
    ("count", {"exit_code": 2, "exit_signal": None, "hits": {}}, False, 2),
    ("count", {"exit_code": 0, "exit_signal": None, "hits": {}}, True, "timeout")])
def test_run_code(mode, g, timed_out, want):
    assert pf.run_code(mode, g, 1000, timed_out) == want


def test_nonrev_bounds_under_p_and_s_by_row():
    rates = ["0.0001"] + ["0.5"] * 18 + ["100"] + ["0.5"] * 360
    for struct in ("PartitionModel", "PartitionModelPlen"):
        got = pf.nonrev_bounds({f"{struct}!NONREV+FO!ModelProtein!rates": ", ".join(rates)})
        assert got == {"at_lower": {"A": 1}, "at_upper": {"R": 1}, "lower": 1, "upper": 1}
    assert pf.nonrev_bounds({"ModelProtein!rates": "1"}) is None


def test_missing_results_never_pass():
    rows = pf.per_eval_rows({("10", "-S", 1, "nqc"): [None], ("10", "-S", 1, "nonrev"): [0.02]})
    assert rows[0]["verdict"] == "missing"
    parts = pf.split({"complete": True, "wall_s": 1.0, "hits": 10}, {"complete": False, "exit_code": "timeout"},
                     0.1, 0.1)
    assert parts["per_fit"] is None and pf.per_fit_row("k", parts)["passed"] is False


def test_window_keys_survive_the_round_trip_and_windows_runs_are_checked():
    key = ("10", "-q", 4, "nqc")
    assert pf.decode_key(pf.encode_key(key)) == key
    w = {"provenance": {"binary": {"sha256": "b"}},
         "results": {"data": [{"passed": True, "data_json_sha256": "d"}], "windows": [{"check": "x", "passed": True}]}}
    assert pf.windows_problems(w, "b", "d") == []
    w["results"]["windows"].append({"check": "w-10 window complete", "passed": False})
    assert pf.windows_problems(w, "c", "e") == ["another binary", "other data", "failed rows: ['w-10 window complete']"]
