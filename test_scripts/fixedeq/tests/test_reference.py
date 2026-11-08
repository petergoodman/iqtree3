"""S5's comparison driver (reference/reference.py): the same-peak line and the outcome table's labels
on constructed inputs (CODE_PLAN.md 3.3)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reference"))
import reference as ref  # noqa: E402

L = 1e-3
H, W = -100.0, -100.01


def test_the_same_peak_line_has_a_floor():
    assert ref.same_peak_line(-5.0, -5.0) == 1e-3
    assert ref.same_peak_line(-5.0, -5.00005) == 1e-3
    assert abs(ref.same_peak_line(-5.0, -5.01) - 0.1) < 1e-12


def test_a_lower_endpoint_that_stopped_short():
    assert ref.label_pair(H, W, (H + 2e-4, True), (H - 1e-4, True), W, L) == "W stopped short"


def test_a_higher_endpoint_not_at_a_peak():
    assert ref.label_pair(H, W, (W, True), (H + 0.5, True), W, L) == "H not at a peak"


def test_distinct_stationary_points_on_one_ridge_or_two_peaks():
    assert ref.label_pair(H, W, (W + 1e-4, True), (H, True), W - 5e-4, L) == "distinct stationary points, one ridge"
    assert ref.label_pair(H, W, (W + 1e-4, True), (H, True), W - 0.2, L) == \
        "distinct stationary points, distinct peaks"


def test_unconverged_or_falling_cross_starts():
    assert ref.label_pair(H, W, (H, False), (H, True), W, L) == "cross-start unconverged"
    assert ref.label_pair(H, W, (W - 0.1, True), (H, True), W, L) == "cross-start unconverged"
    assert ref.label_pair(H, W, (H, True), (H - 0.1, True), W, L) == "cross-start unconverged"


def test_cross_starts_take_fluxes_from_both_fitters():
    import numpy as np
    from oracle import cases, chart, flux, nqc, readers

    rng = np.random.default_rng(4)
    pi = cases.random_target(4, rng)
    Qs = chart.reversible_seed(cases.random_exchangeabilities(4, rng), pi)
    ref_ = chart.references(Qs)
    E = np.eye(4)
    cols = [{name: E[rng.integers(4)] for name in "abcd"} for _ in range(60)]
    obj = nqc.Objective.from_patterns(readers.parse_newick("((a:0.1,b:0.2):0.05,(c:0.3,d:0.1):0.2);"), cols,
                                      rng.integers(1, 6, 60), pi, ref_)
    dom, lf = flux.matched(ref_), flux.LikelihoodFlux(obj, pi)
    pool = ref.Pool(obj, 1)
    r = flux.fit(lf, flux.inward(Qs, pi, dom)[0], dom)
    start = flux.inward(Qs, pi, dom)[0]
    rows = ref.cross_starts("A6", (r.f, r.value, "flux"), (start, lf(start)[0], "z"), obj, pi, ref_, dom, lf, pool, L)
    assert sum(row["check"].startswith("A5 ") for row in rows) == 5   # 3 for the flux fit, 2 for the z climb
    assert rows[-1]["label"] in {"W stopped short", "H not at a peak", "cross-start unconverged", "undetermined",
                                 "distinct stationary points, one ridge", "distinct stationary points, distinct peaks"}


def test_anything_else_is_undetermined():
    # H's fitter from W climbs halfway to H: neither W stopped short nor distinct points
    assert ref.label_pair(H, W, (H - 0.005, True), (H, True), W, L) == "undetermined"
