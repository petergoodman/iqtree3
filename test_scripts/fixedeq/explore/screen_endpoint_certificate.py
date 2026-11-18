"""Exploratory: flux.certificate (activity threshold 1e-6, Peter 2026-11-17, decisions 021 and 039) at the endpoints a
screen run saved (exact.py screen records endpoint_flux), with no refit: for a screen run made before the threshold
changed, such as the 4-taxon S4b repeat at 72077ade.

  python test_scripts/fixedeq/explore/screen_endpoint_certificate.py RESULTS_JSON
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "sim"))
import exact as ex  # noqa: E402
import common  # noqa: E402
from oracle import flux, nqc, readers, simulate  # noqa: E402

res = json.load(open(sys.argv[1]))
prov = res["provenance"]
pi = common.target_values()
_, truths = simulate.truths()
ref = nqc.lg_seed(pi)[1]
dom = flux.matched(ref)
tree_text = prov["tree"]
tree = readers.parse_newick(tree_text)
names = tuple(sorted(leaf.name for leaf in tree.leaves()))
loci = ex.read_loci(prov["loci"]["path"]) if prov.get("loci") else None
for rec in res["fits"]:
    if "endpoint_flux" not in rec:
        continue
    f = np.array(rec["endpoint_flux"])
    w, gaps, _, _ = ex.screen_data(rec["scenario"], truths[rec["base"]], pi, tree_text, loci)
    obj = nqc.Objective.from_patterns(tree, simulate.columns(names, gaps=gaps), w, pi, ref)
    value, gf = flux.LikelihoodFlux(obj, pi)(f)
    cert = flux.certificate(f, gf, dom)
    line = ex.CERTIFICATE * max(1.0, abs(value))
    print(f"{rec['scenario']} {rec['base']}: lnL {value:.6f} (record {rec.get('lnl')}), certificate {cert['residual']:.3g} "
          f"at 1e-6 (line {line:.3g}), active {cert['active']}, passes {cert['residual'] <= line}, C {rec.get('C_fit')}")
