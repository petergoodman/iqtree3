import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import design_rerun as dr  # noqa: E402


@pytest.mark.parametrize("printed, digits", [
    ("0.070", 2), ("-11343.902245", 11), ("1.6e4", 2), ("6146", 4), ("8.02e-5", 3), ("1.0000002", 8), ("-0.0595", 3),
])
def test_significant_digits_count_trailing_zeros_and_ignore_sign(printed, digits):
    assert dr.significant_digits(printed) == digits


@pytest.mark.parametrize("cls, printed, value, ok", [
    ("digits", "0.070", 0.0698, True),
    ("digits", "0.070", 0.0706, False),
    ("digits", "-11343.902245", -11343.9022454, True),
    ("digits", "-11343.902245", -11343.902246, False),
    ("digits", "1.6e4", 16420.0, True),
    ("zero", "2.0e-16", -3e-13, True),
    ("zero", "2.0e-16", 2e-12, False),
    ("bound", "2e-6", 1.9e-6, True),
    ("bound", "2e-6", 2.1e-6, False),
    ("bound", "0", 0.0, True),
    ("bound", "0", 0.01, False),
    ("exact", "TOLX", "gtol", False),
    ("exact", 278, 278, True),
])
def test_the_four_classes(cls, printed, value, ok):
    assert dr.reproduced(cls, printed, value) is ok


def test_the_transcribed_table_follows_decision_019():
    for script, path, printed, cls, doc, line, note in dr.REPORTED:
        assert script in dr.RUNS and cls in ("exact", "digits", "zero", "bound")
        if cls == "zero":
            assert abs(float(printed)) < dr.ZERO, path
        elif cls == "digits":
            assert abs(float(printed)) >= dr.ZERO, path
