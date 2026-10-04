import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "regression"))
import regress  # noqa: E402

REPORT = """IQ-TREE 3.1.4 built Sep 23 2026

Input file name: /some/path/aln.phy
Random seed number: 1

REFERENCES
----------

Cite this.

SUBSTITUTION PROCESS
--------------------

Model of substitution: LG+G4
Gamma shape alpha: 0.6588

USER TREE
---------

Log-likelihood of the tree: -7301.8799 (s.e. 292.5146)

ALISIM COMMAND
--------------
--alisim simulated_MSA -t out.treefile

TIME STAMP
----------

Date and time: Wed Sep 23 11:43:08 2026
"""

ROW = " " + " ".join(["0.1"] * 20)


def q_report(rows):
    return regress.Q_HEADER + " \n\n" + "\n".join(rows) + "\n\nState frequencies: (estimated)\n"


def baseline_of(*item_sets):
    extracts = [(f"rep{i}", {"runs": {"run01": {"items": items}}}) for i, items in enumerate(item_sets, 1)]
    return regress.make_baseline(extracts)["run01"]


def test_compared_sections_exclude_paths_and_times():
    sections = regress.split_sections(REPORT)
    assert [k for k in sections if k not in regress.SKIPPED] == ["REFERENCES", "SUBSTITUTION PROCESS", "USER TREE"]
    assert "/some/path" in sections["(header)"]
    assert "Gamma shape alpha: 0.6588" in sections["SUBSTITUTION PROCESS"]


def test_repeated_title_is_kept_separately():
    sections = regress.split_sections(REPORT + "\nUSER TREE\n---------\nsecond\n")
    assert sections["USER TREE #2"].endswith("second")


def test_q_block_returns_the_21_rows():
    assert regress.q_block(q_report([ROW] * 21)).splitlines() == [ROW] * 21


@pytest.mark.parametrize("text", [
    q_report([ROW] * 20),
    q_report([" " + " ".join(["0.1"] * 19)] * 21),
    q_report([ROW] * 20 + [ROW.replace("0.1", "abc", 1)]),
    "no matrix here\n",
])
def test_q_block_rejects_a_malformed_block(text):
    with pytest.raises(ValueError):
        regress.q_block(text)


def test_split_numbers():
    skeleton, numbers = regress.split_numbers("lnL -7301.8799 (s.e. 292.5) alpha 1e-05 model 12.12")
    assert numbers == ["-7301.8799", "292.5", "1e-05", "12.12"]
    assert skeleton == "lnL # (s.e. #) alpha # model #"


def test_identical_repeats_are_compared_exactly():
    base = baseline_of({"T": "lnL -10.5"}, {"T": "lnL -10.5"})
    assert base["items"]["T"]["mode"] == "exact"
    assert regress.compare_run(base, {"items": {"T": "lnL -10.5"}}) == []
    assert regress.compare_run(base, {"items": {"T": "lnL -10.6"}})


def test_numbers_that_vary_get_their_spread():
    base = baseline_of({"T": "lnL -10.50"}, {"T": "lnL -10.52"})
    assert base["items"]["T"]["mode"] == "spread"
    assert regress.compare_run(base, {"items": {"T": "lnL -10.51"}}) == []
    assert regress.compare_run(base, {"items": {"T": "lnL -10.60"}})
    assert regress.compare_run(base, {"items": {"T": "LL -10.50"}})


def test_text_that_varies_is_unstable_and_never_passes():
    base = baseline_of({"T": "topology A"}, {"T": "topology B"})
    assert base["items"]["T"]["mode"] == "unstable"
    assert regress.compare_run(base, {"items": {"T": "topology A"}})


def test_missing_and_extra_items_fail():
    base = baseline_of({"T": "a", "U": "b"}, {"T": "a", "U": "b"})
    assert regress.compare_run(base, {"items": {"T": "a"}})
    assert regress.compare_run(base, {"items": {"T": "a", "U": "b", "V": "c"}})


RULE = regress.RULES["run11"]


def run11_items(x="0.1", drawing="drawing 1"):
    return {"REFERENCES": "Cite this.",
            "SEQUENCE ALIGNMENT": "Input data: 16 sequences with 2 sites",
            "SUBSTITUTION PROCESS": f"Linked model of substitution: NONREV+FO\nrate {x}",
            "MAXIMUM LIKELIHOOD TREE": drawing,
            "(treefile)": f"((A:{x},B:{x}):{x},C:{x});",
            "(exit code)": "0"}


def run11_baseline():
    return baseline_of(run11_items("0.1", "drawing 1"), run11_items("0.2", "drawing 2"))


def test_run11_rule_ignores_numbers_and_the_tree_drawing():
    new = {"items": run11_items("0.9", "drawing 3")}
    assert regress.compare_run(run11_baseline(), new, RULE) == []
    assert regress.compare_run(run11_baseline(), new)


@pytest.mark.parametrize("changes", [
    {"REFERENCES": "Cite that."},
    {"SEQUENCE ALIGNMENT": "Input data: 16 sequences with 3 sites"},
    {"SUBSTITUTION PROCESS": "Linked model of substitution: GTR20+FO\nrate 0.1"},
    {"(treefile)": "((A:0.1,C:0.1):0.1,B:0.1);"},
    {"(exit code)": "1"},
    {"SUBSTITUTION PROCESS": None},
])
def test_run11_rule_still_catches_wording_topology_and_exit_code(changes):
    items = {k: v for k, v in {**run11_items(), **changes}.items() if v is not None}
    assert regress.compare_run(run11_baseline(), {"items": items}, RULE)


def test_run11_rule_fails_on_an_item_unstable_in_the_baseline():
    base = baseline_of(run11_items(), {**run11_items("0.2"), "(treefile)": "((A:0.2,C:0.2):0.2,B:0.2);"})
    assert base["items"]["(treefile)"]["mode"] == "unstable"
    assert regress.compare_run(base, {"items": run11_items()}, RULE)


def test_baselines_merge_but_a_run_in_two_is_refused(tmp_path):
    for name, runs in [("a.json", {"run11": {}}), ("b.json", {"run12": {}}), ("c.json", {"run12": {}})]:
        (tmp_path / name).write_text(json.dumps({"runs": runs}))
    assert regress.load_baselines([tmp_path / "a.json", tmp_path / "b.json"]) == {"run11": {}, "run12": {}}
    with pytest.raises(SystemExit):
        regress.load_baselines([tmp_path / "b.json", tmp_path / "c.json"])
