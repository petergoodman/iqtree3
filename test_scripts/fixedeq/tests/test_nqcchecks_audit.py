"""The driver checks added after the 2026-10-29 release audit: refusals by their reason, the +I+G
restore line, and linked runs given --all-checkpoint that wrote no numbered checkpoint; and after
the 2026-10-30 audit: verdicts counted whatever their type, NaN residuals, the NQC error line and
the summation order IQ-TREE prints."""
import ast
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "nqc"))
import nqcchecks as nc  # noqa: E402
from oracle import chart, nqc, target  # noqa: E402


class FakeRun:
    """What the checks read from a Runner's run: exit code, stdout, arguments and final checkpoint."""

    def __init__(self, root, code=0, stdout="", args=(), checkpoint=None):
        self.code, self.args, self._checkpoint = code, list(args), checkpoint or {}
        self.dir = root / "run"
        self.dir.mkdir()
        (self.dir / "stdout.txt").write_text(stdout, encoding="utf-8")

    def stdout(self):
        return (self.dir / "stdout.txt").read_text(encoding="utf-8")

    def checkpoint(self):
        return self._checkpoint


def test_a_refusal_must_name_its_reason(tmp_path):
    run = FakeRun(tmp_path, code=2, stdout="ERROR: NQC is not supported with --eigen\n")
    assert nc.refused("x", run, "--eigen")["passed"]
    assert not nc.refused("x", run, "--lie-markov")["passed"]


def test_a_run_that_finishes_is_not_a_refusal(tmp_path):
    run = FakeRun(tmp_path, code=0, stdout="ERROR: NQC is not supported with --eigen\n")
    assert not nc.refused("x", run, "--eigen")["passed"]


def starts_text(starts, restored):
    lines = [f"Init pinv, alpha: 0.1, 1.0 / Estimate: 0.1, 1.0 / LogL: {s}" for s in starts]
    return "\n".join(lines + [f"Optimal pinv,alpha: 0.1, 1.0 / LogL: {restored}"]) + "\n"


def test_restore_line_is_decision_038s_capped_at_the_earlier_line():
    assert nc.restore_line(-7000.0) == 1e-8 * 7000.0
    assert nc.restore_line(-1e5) == 5e-4


def test_restore_check_holds_the_restored_score_to_the_line():
    starts = [-7000.5, -7000.0, -7000.2]
    assert nc.restore_check("x", starts_text(starts, -7000.0 + 1e-5))["passed"]
    assert not nc.restore_check("x", starts_text(starts, -7000.0 + 1e-4))["passed"]
    assert not nc.restore_check("x", starts_text(starts, "nan"))["passed"]
    assert not nc.restore_check("x", starts_text([-1e5 - 0.5, -1e5], -1e5 + 6e-4))["passed"]
    assert nc.restore_check("x", starts_text([-1e5 - 0.5, -1e5], -1e5 + 4e-4))["passed"]


def test_restore_check_reports_whether_the_best_start_was_the_last():
    assert nc.restore_check("x", starts_text([-7000.5, -7000.0, -7000.2], -7000.0))["informative"]
    assert not nc.restore_check("x", starts_text([-7000.5, -7000.0], -7000.0))["informative"]


def test_block_restore_check_judges_each_block():
    text = starts_text([-7000.5, -7000.0, -7000.2], -7000.0) + starts_text([-7000.1, -7000.3], -7000.1 + 1e-3)
    row = nc.block_restore_check("x", text)
    assert not row["passed"] and row["failing"] == [1] and row["blocks"] == 2 and row["informative"] == 2


def linked_checkpoint():
    keys = {"target": "0.5,0.5", "coordinates": "0.1,0.2", "references": "1,0", "residuals": "0,0,0,0,0,0,0"}
    return {f"PartitionModel!{name}!ModelNonrevFixedEq!{k}": v for name in ("NQC", "partA", "partB")
            for k, v in keys.items()}


def first_row(tmp_path, args):
    c = SimpleNamespace(pi=np.array([0.5, 0.5]), ref=[1, 0])
    run = FakeRun(tmp_path, args=args, checkpoint=linked_checkpoint())
    return nc.linked_invariants("x", run, c, 2)[0]


def test_a_run_without_all_checkpoint_is_judged_at_its_final_checkpoint(tmp_path):
    row = first_row(tmp_path, ["-S", "parts.nex"])
    assert row["passed"] and "the final checkpoint" in row["check"]


def test_a_run_given_all_checkpoint_must_have_written_numbered_checkpoints(tmp_path):
    row = first_row(tmp_path, ["-S", "parts.nex", "-cptime", "0", "--all-checkpoint"])
    assert not row["passed"] and row["numbered_expected"]


# ---------------------------------------------------------------- 2026-10-30 audit fixes

def test_a_numpy_verdict_is_stored_as_a_bool_so_that_main_counts_it():
    assert nc.check("x", np.False_)["passed"] is False
    assert nc.check("x", np.True_)["passed"] is True
    assert nc.check("x", None)["passed"] is None


DRIVERS = [Path(__file__).resolve().parents[1] / p
           for p in ("nqc/nqcchecks.py", "profile/profiling.py", "explore/convergence.py")]


def unwrapped_chains(path):
    """Lines of check() calls whose verdict is an and/or chain outside bool(), at the top level or
    as a branch of a conditional verdict: a chain returns its first falsy operand, so a None would
    turn a failing row into a reported one and a numpy False would escape main's count."""
    lines = []
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "check" \
                and len(node.args) > 1:
            verdict = node.args[1]
            parts = [verdict.body, verdict.orelse] if isinstance(verdict, ast.IfExp) else [verdict]
            if any(isinstance(p, ast.BoolOp) for p in parts):
                lines.append(node.lineno)
    return sorted(lines)


@pytest.mark.parametrize("path", DRIVERS, ids=lambda p: p.name)
def test_and_or_verdicts_are_wrapped_in_bool(path):
    assert unwrapped_chains(path) == []


class ExportRun(FakeRun):
    def file(self, suffix):
        return self.dir / f"{nc.PREFIX}{suffix}"


def export_rows(tmp_path, positive_first_entry=False):
    """export_checks on a .qmat the oracle builds from the LG seed's coordinates at the turtle target;
    rows 0 and 3 are the .qmat's form and its residuals."""
    pi = target.parse_target(nc.turtle_target())[0]
    seed, ref = nqc.lg_seed(pi)
    z = chart.inverse(seed, ref)
    Q = chart.build(z, pi, ref)[0]
    if positive_first_entry:
        Q[0, 0] = -Q[0, 0]
    run = ExportRun(tmp_path)
    run.file(".NQC.qmat").write_text("".join(" ".join(f"{v:.17g}" for v in row) + "\n"
                                             for row in [*Q, pi]), encoding="utf-8")
    run.file(".NQC.info").write_text("model: NQC\n", encoding="utf-8")
    ck = {"coordinates": ",".join(f"{v:.17g}" for v in z), "target": ",".join(f"{v:.17g}" for v in pi)}
    return nc.export_checks("x", run, ck, SimpleNamespace(pi=pi, ref=ref), 0)


def test_the_qmat_form_row_passes_on_a_built_matrix(tmp_path):
    rows = export_rows(tmp_path)
    assert rows[0]["passed"] is True and rows[3]["passed"] is True


def test_a_qmat_with_a_positive_first_entry_is_counted_as_failed(tmp_path):
    assert export_rows(tmp_path, positive_first_entry=True)[0]["passed"] is False


def test_a_nan_residual_after_the_first_fails_the_residual_row(tmp_path, monkeypatch):
    residuals = {"stationarity": 1e-13, "row_sums": float("nan"), "mean_rate_error": 1e-13}
    monkeypatch.setattr(nc.chart, "diagnostics", lambda Q, pi: residuals)
    assert export_rows(tmp_path)[3]["passed"] is False


def test_nqc_error_is_a_line_nqc_wrote(tmp_path):
    echo = FakeRun(tmp_path, code=2, stdout="ERROR: File not found NQC+G4\n")
    assert nc.nqc_error(echo) is None
    for i, line in enumerate(("ERROR: NQC: invalid target +F{...}: entry 18",
                              "ERROR: NQC is not supported with --eigen")):
        (tmp_path / str(i)).mkdir()
        run = FakeRun(tmp_path / str(i), code=2, stdout="Reading alignment\n" + line + "\n")
        assert nc.nqc_error(run) == line


# as printed by the 743b6eab gate's permuted-turtle-p run at -vv
ORDER_BLOCK = """0.0054754364255107848,-1.3174503022034523
Partitions ordered by computation costs:
#nexus
begin sets;
  charset ENSGALG00000010671.macse_AA_gb = 389-774;
  charset ENSGALG00000016996.macse_AA_gb = 173-388;
  charset ENSGALG00000016122.macse_AA_gb = 1-172;
end;
NQC first gradient at -lnL 4806.1 (coordinate index, coordinate, step, derivative of -lnL):
"""


def test_summation_order_reads_the_printed_block():
    assert nc.summation_order(ORDER_BLOCK) == ["ENSGALG00000010671.macse_AA_gb", "ENSGALG00000016996.macse_AA_gb",
                                               "ENSGALG00000016122.macse_AA_gb"]
    assert nc.summation_order("NQC first gradient at -lnL 4806.1\n") is None
