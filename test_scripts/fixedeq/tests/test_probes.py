import gzip
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "probes"))
import probes  # noqa: E402
from oracle import readers, target  # noqa: E402

ROOTED = "((A:1.0000000000,B:2.0000000000):0.5000000000,(C:1.0000000000,D:1.0000000000):1.5000000000):0.0000000000;"


def test_checkpoint_reader_joins_struct_and_key(tmp_path):
    p = tmp_path / "x.ckp.gz"
    with gzip.open(p, "wt", encoding="utf-8") as fh:
        fh.write("--- # IQ-TREE Checkpoint ver >= 1.6\nModelProtein:\n rates: 1, 2, 3\nPartitionModelPlen:\n"
                 " part1!ModelProtein!state_freq: 0.25, 0.75\n part_rates: 0.5, 1.5\nfinished: true\n")
    ck = probes.read_checkpoint(p)
    assert ck == {"ModelProtein!rates": "1, 2, 3", "PartitionModelPlen!part1!ModelProtein!state_freq": "0.25, 0.75",
                  "PartitionModelPlen!part_rates": "0.5, 1.5", "finished": "true"}
    assert probes.ckp_values(ck["PartitionModelPlen!part_rates"]) == [0.5, 1.5]


def test_root_info_and_unroot_keep_the_unrooted_tree():
    t = readers.parse_newick(ROOTED)
    split, lengths = probes.root_info(t)
    assert split == frozenset("CD") and lengths == (1.5, 0.5)
    u = probes.unroot(t)
    assert probes.root_info(u) is None and len(u.children) == 3
    assert probes.edge_map(u) == probes.edge_map(t)
    assert probes.to_newick(u) == "(A:1.0000000000,B:2.0000000000,(C:1.0000000000,D:1.0000000000):2.0000000000);"


@pytest.mark.parametrize("start", ["A", "B", "D"])
def test_midpoint_follows_convert_to_rooted(start):
    # longest path A-C = 1 + 1 + 4 = 6; its midpoint lies on the edge to C, 3 from C and 1 from the centre
    t = readers.parse_newick("((A:1,B:0.5):1,C:4,D:1);")
    split, lengths = probes.midpoint(t, start)
    assert split == frozenset("C")
    assert lengths == pytest.approx((3.0, 1.0))


def test_tree_comparison_against_the_midpoint():
    unrooted = "((A:1,B:0.5):1,C:4,D:1);"
    rooted_out = "(((A:1,B:0.5):1,D:1):1,C:3):0;"
    c = probes.tree_comparison(unrooted, rooted_out, start="D")
    assert c["same_topology"] and c["max_abs_change"] == 0.0
    assert c["root_split_kept"] and c["root_lengths_kept"]
    moved = probes.tree_comparison(unrooted, "(((A:1,B:0.5):1,C:4):0.5,D:0.5):0;", start="D")
    assert moved["same_topology"] and not moved["root_split_kept"]


def test_paml_rows_invert_the_internal_pair_order():
    # internal order for n = 4: (0,1) (0,2) (0,3) (1,2) (1,3) (2,3)
    assert probes.paml_rows(list(range(6)), n=4) == [[0], [1, 3], [2, 4, 5]]


def test_charpartition_and_mdef_text():
    nex = "#nexus\nbegin sets;\ncharset p1 = 1-10;\ncharset p2 = 11-20;\nend;\n"
    text = probes.with_charpartition(nex, [("LG+G4{0.5}", "p1", "0.8"), ("LG+G4{0.5}", "p2", None)])
    assert "charpartition probe = LG+G4{0.5}: p1{0.8}, LG+G4{0.5}: p2;\nend;" in text
    assert probes.mdef("X", "0.5,0.5") == "#nexus\nbegin models;\nfrequency X = 0.5,0.5;\nend;\n"


def test_partition_table():
    report = ("SUBSTITUTION PROCESS\n--------------------\n\nTopology-unlinked partition model\n\n"
              "  ID  Model         TreeLen  Parameters\n   1  NONREV+FO      0.4336  NONREV+FO+G4{0.5}\n"
              "   2  NONREV+FO      1.0772  NONREV+FO\n\nLinked model of substitution: NONREV+FO\n")
    title, rows = probes.partition_table(report)
    assert title == "TreeLen"
    assert rows == [("NONREV+FO", "0.4336", "NONREV+FO+G4{0.5}"), ("NONREV+FO", "1.0772", "NONREV+FO")]


def test_turtle_target_is_a_valid_plain_decimal_vector():
    text = probes.turtle_target()
    assert "e+" not in text and "E+" not in text
    pi, raw = target.parse_target(text)
    assert len(pi) == 20 and abs(raw - 1.0) < 1e-14
