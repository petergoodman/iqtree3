"""The unit-test driver's reading of googletest's XML (test_scripts/fixedeq/unittests/unittests.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "unittests"))
import unittests  # noqa: E402

XML = """<?xml version="1.0" encoding="UTF-8"?>
<testsuites tests="3" failures="{failures}" disabled="0" errors="0" name="AllTests">
  <testsuite name="Oracle/ChartFixture" tests="2" failures="{failures}" errors="0">
    <testcase name="MatchesOracle/a" value_param="&quot;a&quot;" status="run" result="completed" time="0" classname="Oracle/ChartFixture">
      <properties>
        <property name="fixture" value="a"/>
        <property name="build_error" value="2e-16"/>
        <property name="q_limit" value="1e-12"/>
      </properties>
    </testcase>
    <testcase name="MatchesOracle/b" status="run" result="completed" time="0" classname="Oracle/ChartFixture" fixture="b" build_error="5e-13" q_limit="1e-12">
      {failure}
    </testcase>
  </testsuite>
  <testsuite name="Step" tests="1" failures="0" errors="0">
    <testcase name="StepsBack" status="run" result="completed" time="0" classname="Step"/>
  </testsuite>
</testsuites>
"""


def test_rows_come_from_property_children_and_attributes(tmp_path):
    path = tmp_path / "gtest.xml"
    path.write_text(XML.format(failures=0, failure=""))
    problems, rows = unittests.results(path, ["a", "b"])
    assert problems == []
    assert [(r["fixture"], r["status"], r["build_error"]) for r in rows] == [("a", "passed", "2e-16"), ("b", "passed", "5e-13")]
    assert unittests.maxima(rows)["build_error_over_limit"] == 0.5


def test_failures_and_missing_rows_are_problems(tmp_path):
    path = tmp_path / "gtest.xml"
    path.write_text(XML.format(failures=1, failure='<failure message="x" type=""/>'))
    problems, rows = unittests.results(path, ["a", "b", "c"])
    assert [r["status"] for r in rows] == ["passed", "failed"]
    assert any("1 of 3 tests failed" in p for p in problems) and any("c" in p for p in problems if "no row" in p)


def test_a_missing_or_truncated_xml_is_a_problem(tmp_path):
    assert unittests.results(tmp_path / "absent.xml", ["a"])[0]
    cut = tmp_path / "cut.xml"
    cut.write_text(XML.format(failures=0, failure="")[:300])
    assert unittests.results(cut, ["a"])[0]
