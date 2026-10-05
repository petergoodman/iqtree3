"""The frequency estimator port and the target parser (decisions 006 and 021)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import design_rerun  # noqa: E402
from oracle import readers, target  # noqa: E402

REPO = Path(__file__).resolve().parents[3]

# "Mean state frequencies:" as printed, identically, by baseline runs 8 (-p) and 9 (-S) of the
# frozen binary (~/iqtree3-runs/baseline/rep1/run08/stdout.txt and run09/stdout.txt, recorded
# 2026-10-04 at driver commit 086fef2f); fixed format with 8 decimals (partitionmodel.cpp:155-160).
PRINTED = ("0.05952707 0.04208893 0.03825436 0.05788368 0.02903314 0.04053684 0.08655163 0.05925317 "
           "0.02665936 0.06135305 0.08554734 0.07340455 0.02793755 0.05587510 0.04208893 0.04446270 "
           "0.04455400 0.00876472 0.04281932 0.07340455").split()


@pytest.mark.parametrize("partition_type", ["p", "S"])
def test_pooled_turtle_frequencies_match_iqtree(partition_type):
    aln = readers.read_fasta(REPO / "test_scripts" / "test_data" / "turtle_aa.fasta")
    sets = readers.read_charsets(REPO / "test_scripts" / "test_data" / "turtle_aa.nex")
    freq = target.pooled_frequencies(aln, sets, partition_type)
    assert all(design_rerun.reproduced("digits", p, f) for p, f in zip(PRINTED, freq)), list(freq)


def test_empty_sequences_are_removed_but_three_are_kept():
    aln = {"a": "--", "b": "AC", "c": "??", "d": "X-"}
    assert list(target.remove_gappy(aln)) == ["a", "b", "c"]
    assert list(target.remove_gappy({"a": "A", "b": "-", "c": "C", "d": "D"})) == ["a", "c", "d"]


def test_codes_and_appearance_follow_iqtree():
    assert [target.state_code(c) for c in "AVBZJX-?"] == [0, 19, 20, 21, 22, 23, 23, 23]
    assert list(np.nonzero(target.appearance(20))[0]) == [2, 3] and target.appearance(23).sum() == 20


def test_definite_states_give_their_proportions_and_zero_stays_zero_by_default():
    counts = np.zeros(24, dtype=int)
    counts[[0, 1, 2]] = [2, 1, 1]
    freq = target.convert_count_to_freq(counts)
    assert freq[3] == 0.0 and list(freq[:3]) == [0.5, 0.25, 0.25]


def test_floor_applies_only_when_zero_frequencies_are_not_kept():
    freq = np.array([0.5, 0.3] + [0.2 / 17] * 18)
    freq[5] = 0.0  # a normalized vector with one zero entry
    floored = target.convfreq(freq.copy(), 1e-4)
    assert floored[5] == 1e-4 and floored[0] < 0.5
    assert np.array_equal(np.delete(floored, [0, 5]), np.delete(freq, [0, 5]))
    assert abs(floored.sum() - 1.0) < 1e-15


def test_target_parser_accepts_and_normalizes_once():
    text = "/".join(["0.05"] * 10) + " " + ",".join(["0.05"] * 9) + ",0.0500005"
    pi, raw = target.parse_target(text)
    assert raw == pytest.approx(1.0000005, abs=1e-15) and pi.sum() == pytest.approx(1.0, abs=1e-15)


@pytest.mark.parametrize("text", [
    ",".join(["0.05"] * 19),
    ",".join(["0.05"] * 19 + ["nan"]),
    ",".join(["0.05"] * 19 + ["inf"]),
    ",".join(["0.05"] * 19 + ["abc"]),
    ",".join(["0.1"] + ["0.05"] * 18 + ["0.0"]),
    ",".join(["0.05"] * 19 + ["0.00005"]),
    ",".join(["0.05"] * 19 + ["0.051"]),
])
def test_target_parser_rejects(text):
    with pytest.raises(ValueError):
        target.parse_target(text)
