import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "sanitizer"))
import sanitize  # noqa: E402

ASAN_REPORT = """some output
==1234==ERROR: AddressSanitizer: heap-buffer-overflow on address 0x602000000014 at pc 0x55 bp 0x7f sp 0x7f
READ of size 8 at 0x602000000014 thread T0
    #0 0x55d5c in __interceptor_memcpy /usr/src/compiler-rt/lib/asan/asan_interceptors.cpp:22:3
    #1 0x55d6d in std::vector<double, std::allocator<double> >::at(unsigned long) const /mnt/c/x/tree/phylotree.cpp:123:45
    #2 0x55d7e in main /mnt/c/x/main/main.cpp:10:2
"""


def test_asan_report_is_grouped_at_its_first_project_frame():
    found = sanitize.scan_text(ASAN_REPORT)
    assert found == [("AddressSanitizer", "/mnt/c/x/tree/phylotree.cpp:123 "
                      "(std::vector<double, std::allocator<double> >::at(unsigned long) const)", "heap-buffer-overflow")]


def test_ubsan_values_are_masked_so_repeats_group_together():
    a = sanitize.scan_text("/mnt/c/x/model/a.cpp:7:9: runtime error: signed integer overflow: 2147483647 + 1 cannot be represented in type 'int'")
    b = sanitize.scan_text("/mnt/c/x/model/a.cpp:7:9: runtime error: signed integer overflow: 2147483600 + 99 cannot be represented in type 'int'")
    assert a == b and a[0][:2] == ("UndefinedBehaviorSanitizer", "/mnt/c/x/model/a.cpp:7")


def test_iqtree_crash_message_and_clean_output():
    assert sanitize.scan_text("*** IQ-TREE CRASHES WITH SIGNAL SEGMENTATION FAULT") == [
        ("IQ-TREE crash", None, "signal SEGMENTATION")]
    assert sanitize.scan_text("Total wall-clock time used: 1.2 seconds\n") == []


def test_build_environment_drops_conda(monkeypatch):
    monkeypatch.setenv("PATH", "/home/u/anaconda3/envs/x/bin:/usr/bin:/home/u/anaconda3/condabin")
    monkeypatch.setenv("CONDA_PREFIX", "/home/u/anaconda3/envs/x")
    monkeypatch.setenv("CMAKE_PREFIX_PATH", "/home/u/anaconda3/envs/x")
    env = sanitize.build_env()
    assert env["PATH"] == "/usr/bin" and "CONDA_PREFIX" not in env and "CMAKE_PREFIX_PATH" not in env
