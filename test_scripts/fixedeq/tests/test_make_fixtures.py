"""The S1 fixture generator: reproducible, read back bitwise, and all-or-nothing (CODE_PLAN.md 3.2)."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import make_fixtures as mf  # noqa: E402

PROV = {"script": mf.SCRIPT, "git_commit": "test", "git_clean": False, "generated_utc": "test", "seed": mf.SEED}
TURTLE = json.loads(mf.MANIFEST.read_text())["manifest"]["target"]["values_17_digits"]


@pytest.fixture(scope="module")
def fixtures():
    return mf.build_cases()


def test_cases_are_reproducible_from_the_seed(fixtures):
    again = mf.build_cases()
    assert [c["name"] for c in again] == [c["name"] for c in fixtures]
    for a, b in zip(fixtures, again):
        assert a.keys() == b.keys()
        assert all(np.array_equal(a[k], b[k]) for k in a if isinstance(a[k], np.ndarray)), a["name"]


def test_every_case_and_target_passes_its_checks(fixtures):
    assert len(fixtures) == 27 and len({c["name"] for c in fixtures}) == 27
    for c in fixtures:
        assert mf.check(c)[1] == [], c["name"]
    assert mf.render_targets(TURTLE, PROV)[1] == []


def test_written_fixtures_read_back_bitwise(fixtures):
    for c in fixtures:
        header, arrays = mf.read_fixture(mf.render(c, PROV)[0])
        assert header["name"] == c["name"] and int(header["n"]) == len(c["pi"])
        for label in ("pi", "ref", "z", "nu", "Q", "R", "Q0", "R_ref"):
            if label in c:
                assert np.array_equal(arrays[label], np.ravel(c[label])), (c["name"], label)


def test_a_failing_case_writes_nothing(fixtures, tmp_path):
    bad = dict(fixtures[0], name="bad", Q=fixtures[0]["Q"].copy())
    bad["Q"][0, 1] *= 1 + 1e-9
    out = tmp_path / "fixtures"
    summary, failures = mf.write_all(out, [fixtures[1], bad], TURTLE, PROV)
    assert summary == [] and failures and all(f.startswith("bad:") for f in failures)
    assert not out.exists()


def test_a_non_empty_output_directory_is_refused(tmp_path, monkeypatch):
    (tmp_path / "old.txt").write_text("x")
    monkeypatch.setattr(sys, "argv", ["make_fixtures.py", "--out", str(tmp_path), "--allow-dirty", "--git", "no-such-git"])
    with pytest.raises(SystemExit, match="never regenerated"):
        mf.main()
    assert [p.name for p in tmp_path.iterdir()] == ["old.txt"]
