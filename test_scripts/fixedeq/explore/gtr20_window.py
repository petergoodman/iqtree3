"""Decision 046's laptop demonstration and calibration, and the read of the seeded arm's smoke run (CODE_PLAN.md 3.8,
"SIM-2's seeded GTR20 arm", laptop runs); exploratory.

  python test_scripts/fixedeq/explore/gtr20_window.py fits --frozen BIN --gate BIN --git git.exe [--cap-hours 12]
      [--out DIR]
  python test_scripts/fixedeq/explore/gtr20_window.py smoke --run RUNDIR --reference RUNDIR --git git.exe [--out DIR]

fits: R0 at 30 genes, replicate 1, simulated by the frozen binary (seed base sim2.BASE; its genes checked against
the committed cross-machine fixture), then four GTR20+F{pi*} fits at SIM-2's settings (-q, the tree unrooted,
-blfix, -keep-ident, 20 rounds, -seed 1), 1 thread each, run at the same time so that they share one load: the
default seed on the frozen and on the gate binary, and seeded (pipeline.gtr20_seed) at p = 0.02 and 0.05 on the
gate binary. For each: the oracle's training lnL (sim2.objective_for, the trees AliSim used) against the truth's,
the fixed rate, the largest free rate, headroom and floor count (pipeline.rate_bounds), rounds, ending rule and
wall time. Asserted, decision 021's laptop rows: the two default fits' checkpoint rates text and final lnL
identical; each seeded fit's fixed rate at p and no free rate above 99; each seeded fit's training lnL not below
the truth's. Reported: whether the limit binds for the default fits at 30 genes; p by the rule set before the run
(choose_p); the seeded-to-default wall ratio at that p, which scales sim2_gtr20.py's estimates.
smoke: a seeded smoke run (sim2_gtr20.py run --smoke) against round 1's smoke run on the same data: every phase
done; the shard's seeded rows passed; every dataset's file hashes equal; each gtr20f fit's linked frequencies equal
round 1's at the checkpoint's 10 digits; each gtr20pi fit's equal pi* at 10 digits; the seed's printed reads
against the gene count, reported.
"""
import argparse
import concurrent.futures
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "sim"))
import sim2  # noqa: E402  (BLAS threads at 1; puts the drivers on sys.path)
import sim2_gtr20  # noqa: E402
import common  # noqa: E402
import pipeline as pl  # noqa: E402
from oracle import chart, readers  # noqa: E402
from nqcchecks import check  # noqa: E402
from probes import PREFIX, read_checkpoint  # noqa: E402
from regress import REPO, read  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
GENES, REPLICATE = 30, 1
PS = (0.02, 0.05)
LNL_MARGIN, HEADROOM_MIN = 1.0, 4.0     # the rule for p (CODE_PLAN 3.8), set before the run


def choose_p(lnl, headroom):
    """0.02, unless the 0.05 fit's training lnL is higher by more than 1 and its headroom is at least 4; then 0.05.
    lnl and headroom map each p to its seeded fit's values."""
    return 0.05 if lnl[0.05] - lnl[0.02] > LNL_MARGIN and headroom[0.05] >= HEADROOM_MIN else 0.02


def linked_freq_text(ck):
    """The linked GTR20 model's frequencies as the checkpoint writes them (10 significant digits)."""
    hits = [v for k, v in ck.items() if k.endswith("!ModelProtein!state_freq") and k.split("!")[1] in ("GTR20", "GTR20+F")]
    if len(hits) != 1:
        raise ValueError(f"{len(hits)} linked GTR20 frequency entries")
    return [v.strip() for v in hits[0].split(",")]


def cmd_fits(a):
    prov = common.provenance(a, SCRIPT, "window", "exploratory", genes=GENES, replicate=REPLICATE, seed_base=sim2.BASE,
                             iqtree_seed=1, threads=1, rounds=sim2.ROUNDS, cap_hours=a.cap_hours, ps=PS,
                             rule={"lnl_margin": LNL_MARGIN, "headroom_min": HEADROOM_MIN})
    out = common.new_out(a.out, prov, runs=sim2_gtr20.RUNS)
    bad = sim2_gtr20.binary_check(a.frozen, a.gate)
    if bad:
        sys.exit(f"not round 1's binaries (decision 046): {bad}")
    frozen, gate, prov["binaries"] = sim2.binaries(a)
    title = "SIM-2 seeded GTR20: laptop demonstration and calibration"
    tree, loci = pl.read_scaffold(sim2.SCAFFOLD)
    data = out / "data"
    data.mkdir()
    man, rows = pl.simulate_into(data, frozen, tree, loci, sim2.TRUTHS["R0"], GENES, sim2.BASE, REPLICATE)
    genes = [g["sha256"] for g in man["genes"]]
    rows.append(check("data: R0's genes byte-identical to the committed cross-machine fixture (021)",
                      genes == json.loads(read(sim2.FIXTURE))["genes"]["R0"]))
    if any(r["passed"] is False for r in rows):
        common.finish(out, prov, f"{title}: the data failed", rows)

    target = common.target_text()
    specs = {"default-frozen": (frozen, None), "default-gate": (gate, None),
             **{f"seeded-{p:g}": (gate, ("seed.paml", pl.gtr20_seed(p, target))) for p in PS}}

    def run(name):
        (out / name).mkdir()
        binary, init = specs[name]
        return name, pl.fit_into(out / name, data, binary, 1, a.cap_hours * 3600, gtr20="pi", rounds=sim2.ROUNDS,
                                 init_model=init)

    with concurrent.futures.ThreadPoolExecutor(len(specs)) as pool:
        recs = dict(pool.map(run, specs))

    pi = common.target_values()
    obj = sim2.objective_for(data, read(data / pl.TREE_RATE1).strip(), pi)[0]
    truth = obj.log_likelihood(readers.read_qmat(read(data / "truth.qmat"))[0])
    texts, fits = {}, {}
    for name, rec in recs.items():
        stdout = out / name / "fit" / "stdout.txt"
        worse = stdout.exists() and "Tree logl gets worse" in read(stdout)
        ok = rec["exit_code"] == 0 and not rec["at_cap"] and not rec.get("assertion_or_abort") and not worse
        rows.append(check(f"{name}: the fit exits 0 within its cap, with no abort and no fall (sim2.do_phase's failures)",
                          ok, exit_code=rec["exit_code"], at_cap=rec["at_cap"],
                          assertion_or_abort=rec.get("assertion_or_abort"), logl_worse=worse))
        if not ok:
            continue
        st = pl.gtr20_state(read_checkpoint(out / name / "fit" / f"{PREFIX}.ckp.gz"))
        lnl = obj.log_likelihood(chart.reversible_seed(st["R"], pi))
        texts[name] = st["rates_text"]
        fits[name] = {**pl.rate_bounds(st), "lnl_oracle": lnl, "lnl_truth_oracle": truth, "excess_over_truth": lnl - truth,
                      "final_lnl_printed": rec["final_lnl"], "rounds": rec["rounds"], "ended_by": rec["ended_by"],
                      "wall_s": rec["wall_s"], "cpu_s": rec["cpu_s"], "peak_rss_mb": rec["peak_rss_mb"]}
        rows.append(check(f"{name}: record", None, **fits[name]))
    if "default-frozen" in fits and "default-gate" in fits:
        same = texts["default-frozen"] == texts["default-gate"]
        finals = [fits[k]["final_lnl_printed"] for k in ("default-frozen", "default-gate")]
        rows.append(check("default seed: frozen and gate binaries give identical checkpoint rates text and final lnL "
                          "(021)", same and finals[0] is not None and finals[0] == finals[1], rates_text_identical=same,
                          final_lnl=finals))
        g = fits["default-gate"]
        rows.append(check("default seed: the range limit (reported; if it does not bind at 30 genes the 100-gene "
                          "evidence stands)", None, fixed_text=g["fixed_text"], above_upper=g["above_upper"],
                          binds=g["above_upper"] > 0, excess_over_truth=g["excess_over_truth"]))
    for p in PS:
        f = fits.get(f"seeded-{p:g}")
        if f is None:
            continue
        rows.append(check(f"seeded p = {p:g}: the fixed Tyr-Val's checkpoint text equals p's, and no free rate above 99 "
                          "(021)", f["fixed_text"] == f"{p:.10g}" and f["above_upper"] == 0, fixed_text=f["fixed_text"],
                          above_upper=f["above_upper"], largest_free=f["largest_free"], headroom=f["headroom"]))
        rows.append(check(f"seeded p = {p:g}: training lnL not below the truth's, by the oracle (021)",
                          f["lnl_oracle"] >= truth, excess_over_truth=f["excess_over_truth"]))
    extra = {"fits": fits}
    if all(f"seeded-{p:g}" in fits for p in PS) and "default-gate" in fits:
        chosen = choose_p({p: fits[f"seeded-{p:g}"]["lnl_oracle"] for p in PS},
                          {p: fits[f"seeded-{p:g}"]["headroom"] for p in PS})
        ratio = fits[f"seeded-{chosen:g}"]["wall_s"] / fits["default-gate"]["wall_s"]
        extra.update({"p": chosen, "ratio": ratio})
        rows.append(check("p by the rule set before the run: 0.02 unless 0.05's lnL is higher by more than 1 and its "
                          "headroom at least 4", None, p=chosen,
                          lnl_gain_005=fits["seeded-0.05"]["lnl_oracle"] - fits["seeded-0.02"]["lnl_oracle"],
                          headroom_005=fits["seeded-0.05"]["headroom"]))
        rows.append(check("seeded-to-default wall ratio at that p (scales sim2_gtr20.py's estimates)", None, ratio=ratio,
                          seeded_wall_s=fits[f"seeded-{chosen:g}"]["wall_s"],
                          default_wall_s=fits["default-gate"]["wall_s"]))
    common.finish(out, prov, title, rows, extra)


def cmd_smoke(a):
    run, ref = (Path(p).expanduser().resolve() for p in (a.run, a.reference))
    prov = common.provenance(a, SCRIPT, "smoke-read", "development", run=str(run), reference=str(ref))
    out = common.new_out(a.out, prov, runs=sim2_gtr20.RUNS)
    des = json.loads(read(run / "design.json"))
    shard = run / "shard-1"
    last = max(shard.glob("run-*"), key=lambda p: int(p.name.split("-")[1]))
    res = json.loads(read(last / "results.json"))
    phases = sim2.shard_phases(des, 1)
    rows = [check("every phase done", all(res["state"].get(pid) == "done" for pid in phases),
                  not_done={pid: res["state"].get(pid) for pid in phases if res["state"].get(pid) != "done"})]
    seeded = [r for r in res["results"] if "seeded GTR20" in r["check"]]
    rows.append(check("the shard's seeded rows: four for each task and none failed",
                      len(seeded) == 4 * len(des["tasks"]) and all(r["passed"] is not False for r in seeded),
                      rows=len(seeded), tasks=len(des["tasks"]),
        failed=[r["check"] for r in seeded if r["passed"] is False],
        fixed_text={r["check"].split(":")[0]: r.get("fixed_text") for r in seeded if "fixed Tyr-Val" in r["check"]}))
    rows.append(check("the seed's printed reads against the genes, and boundary warnings (printed output)", None,
                      **{r["check"].split(":")[0]: {k: r[k] for k in ("seed_reads", "genes", "boundary_warnings")}
                         for r in seeded if "printed seed reads" in r["check"]}))
    for d in des["datasets"]:
        mine, theirs = (json.loads(read(x / "shard-1" / "data" / d / "data.json"))["files"] for x in (run, ref))
        rows.append(check(f"data/{d}: file hashes equal round 1's smoke run", mine == theirs))
    pi_text = [f"{v:.10g}" for v in common.target_values()]
    for t in des["tasks"]:
        old = t["id"].removesuffix("-seeded")
        try:
            mine = linked_freq_text(read_checkpoint(shard / "tasks" / t["id"] / "fit" / "fit" / f"{PREFIX}.ckp.gz"))
            if t["method"] == "gtr20pi":
                rows.append(check(f"{t['id']}: linked frequencies equal pi* at 10 digits", mine == pi_text))
            else:
                theirs = linked_freq_text(read_checkpoint(ref / "shard-1" / "tasks" / old / "fit" / "fit" /
                                                          f"{PREFIX}.ckp.gz"))
                rows.append(check(f"{t['id']}: linked frequencies equal round 1's {old} at the checkpoint's 10 digits",
                                  mine == theirs))
        except Exception as e:  # noqa: BLE001  a missing or unreadable checkpoint fails its row; the report is kept
            rows.append(check(f"{t['id']}: linked frequencies readable", False, error=repr(e)[:500]))
    common.finish(out, prov, "SIM-2 seeded GTR20: smoke run against round 1's", rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    f = sub.add_parser("fits")
    f.add_argument("--frozen", required=True)
    f.add_argument("--gate", required=True)
    f.add_argument("--cap-hours", type=float, default=12.0)
    s = sub.add_parser("smoke")
    s.add_argument("--run", required=True, help="the seeded smoke run's directory")
    s.add_argument("--reference", required=True, help="round 1's smoke run directory")
    for p in (f, s):
        p.add_argument("--out", help=f"new directory outside the repository (default under {sim2_gtr20.RUNS}/)")
        p.add_argument("--git", default="git")
        p.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    {"fits": cmd_fits, "smoke": cmd_smoke}[a.step](a)


if __name__ == "__main__":
    main()
