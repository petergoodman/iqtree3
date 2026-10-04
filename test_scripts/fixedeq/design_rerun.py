#!/usr/bin/env python3
"""Rerun the three design scripts from copies outside the repository and compare their outputs
with the values the design documents report (decision 019; docs/agent/CODE_PLAN.md, section 3.1).

  python design_rerun.py --out DIR [--git git.exe]

DIR must be new and outside the repository. The scripts in docs/agent/design/scripts/ are copied
there and run from the copies, never in place. Exit status 1 if a script exits nonzero or any
documented value is not reproduced.
"""
import argparse
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "regression"))
from regress import REPO, capture, sha256, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
DESIGN = REPO / "docs" / "agent" / "design"
BASE_COMMIT = "63c330d9"
ZERO = 1e-12  # decision 019: a documented value below this is rounding noise around zero

# from docs/agent/design/README.md
FINGERPRINTS = {
    "verify_constrained_nq.py": "fef4315dd728712197a708893816cd5ef9c25b558df0736336095b0166d56b4b",
    "check_builtin_matrices.py": "71c120b94d1e68019faa927f3c32c3243ed33871de1073e8fa9f3573c5f4233c",
    "chart_fd_bfgs_compare.py": "8f2a1d095b24efcef62cff881cea213a924d4af7f14279517ee4d90aa9d0591b",
}
V, C, E = FINGERPRINTS
# script: (arguments, added environment, file its results are read from), as the documents state
RUNS = {
    V: (["--n", "20", "--seed", "20260922"], {}, "verify_constrained_nq_results.json"),
    C: (["modelprotein.cpp"], {}, "check_builtin_matrices_results.json"),
    E: (["--n", "8", "--ntips", "8", "--nsites", "1000", "--seed", "20260922", "--nrep", "4"],
        {"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}, "chart_fd_bfgs_compare.stdout"),
}

P = "constrained_nq_plan_agent.md"
S = "project2_unified_discrepancy_synthesis.md"
READ = "interpreted: "
# (script, output path, printed value, class, document, line, note). An "exact" entry holds the
# expected value itself; the others hold the number as printed, whose digits set the precision.
REPORTED = [
    (V, ("(exit code)",), 0, "exact", P, 299, "'no hard failures'"),
    (V, ("_meta", "failures"), [], "exact", P, 299, "'no hard failures'"),
    (V, ("1.rank_stationarity_map",), 19, "exact", P, 303, ""),
    (V, ("1.rank_with_mean_rate_row",), 20, "exact", P, 303, ""),
    (V, ("1.dim_fiber",), 360, "exact", P, 303, ""),
    (V, ("1.dim_reversible",), 189, "exact", P, 303, ""),
    (V, ("1.dim_surplus",), 171, "exact", P, 303, ""),
    (V, ("2.c2_stationarity_residual",), "2.0e-16", "zero", P, 304, ""),
    (V, ("2.c2_mean_rate_minus_1",), "5.6e-16", "zero", P, 304, ""),
    (V, ("2.c2_roundtrip_P",), "1.1e-16", "zero", P, 304, ""),
    (V, ("2.c2_nu_equals_pi_times_exit",), "1.4e-17", "zero", P, 304, ""),
    (V, ("2.c2_row_scale_is_gauge",), "7.1e-15", "zero", P, 304, ""),
    (V, ("2.c2_jacobian_rank",), 360, "exact", P, 304, ""),
    (V, ("2.c2_dnu_formula_error",), "1.1e-10", "digits", P, 304, ""),
    (V, ("3.c1_stationarity_residual",), "4.5e-17", "zero", P, 305, ""),
    (V, ("3.c1_mean_rate_minus_1",), "2.2e-16", "zero", P, 305, ""),
    (V, ("3.c1_newton_iterations_cold",), 5, "exact", P, 305, ""),
    (V, ("3.c1_gauge_fixed_roundtrip",), "4.4e-16", "zero", P, 305, ""),
    (V, ("3.c1_gauge_orbit_invariance",), "6.7e-16", "zero", P, 305, ""),
    (V, ("3.c1_h0_equals_R_diag_pi",), "0.0", "zero", P, 305, READ + "one printed 'h = 0 slice' value, two script keys"),
    (V, ("3.c1_h0_gstar_is_zero",), "0.0", "zero", P, 305, READ + "one printed 'h = 0 slice' value, two script keys"),
    (V, ("3.c1_jacobian_rank",), 360, "exact", P, 305, ""),
    (V, ("3.c1_c2_agree_on_fiber_point",), "4.5e-14", "zero", P, 305, ""),
    (V, ("4.counterexample_uniform_stationary",), "5e-16", "zero", P, 306, ""),
    (V, ("4.raw_h_max",), "6.91", "digits", P, 306, ""),
    (V, ("4.gauge_fixed_h_max",), "20.72", "digits", P, 306, ""),
    (V, ("4.max_flux_entry_over_samples",), "0.40", "digits", P, 306, ""),
    (V, ("5.circulant_formula_error",), "2.2e-16", "zero", P, 307, ""),
    (V, ("5.second_derivative_numeric",), "0.015378", "digits", P, 307, ""),
    (V, ("5.second_derivative_analytic",), "0.015378", "digits", P, 307, ""),
    (V, ("6.fd_rule_on_signed_coordinate", "h=0.01", "iqtree_fd_error"), "4.1e-5", "digits", P, 308, ""),
    (V, ("6.fd_rule_on_signed_coordinate", "h=0.0001", "iqtree_fd_error"), "4.7e-4", "digits", P, 308, ""),
    (V, ("6.fd_rule_on_signed_coordinate", "h=1e-06", "iqtree_fd_error"), "0.070", "digits", P, 308, ""),
    (V, ("6.fd_rule_on_signed_coordinate", "h=1e-08", "iqtree_fd_error"), "6.5", "digits", P, 308, ""),
    (V, ("6.fd_rule_on_signed_coordinate", "h=0.01", "reference_derivative"), "-12.87", "digits", P, 308,
     READ + "'reference -12.87, then -13.75': first h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=0.0001", "reference_derivative"), "-13.75", "digits", P, 308,
     READ + "'then -13.75': the later h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=1e-06", "reference_derivative"), "-13.75", "digits", P, 308,
     READ + "'then -13.75': the later h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=1e-08", "reference_derivative"), "-13.75", "digits", P, 308,
     READ + "'then -13.75': the later h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=0.01", "shifted_by_12_error"), "0.053", "digits", P, 308,
     READ + "'0.053 throughout': every h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=0.0001", "shifted_by_12_error"), "0.053", "digits", P, 308,
     READ + "'0.053 throughout': every h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=1e-06", "shifted_by_12_error"), "0.053", "digits", P, 308,
     READ + "'0.053 throughout': every h"),
    (V, ("6.fd_rule_on_signed_coordinate", "h=1e-08", "shifted_by_12_error"), "0.053", "digits", P, 308,
     READ + "'0.053 throughout': every h"),
    (V, ("6.fd_rule_on_positive_ratio_coordinate", "u=1", "iqtree_fd_relative_error"), "1.6e-4", "digits", P, 308, ""),
    (V, ("6.fd_rule_on_positive_ratio_coordinate", "u=0.01", "iqtree_fd_relative_error"), "2.1e-6", "digits", P, 308, ""),
    (V, ("6.fd_rule_on_positive_ratio_coordinate", "u=0.0001", "iqtree_fd_relative_error"), "2.7e-5", "digits", P, 308, ""),
    (V, ("7.swap_reversible_seed_stationary_at_new_pi",), "3.5e-17", "zero", P, 309, ""),
    (V, ("7.swap_nonreversible_seed_residual_at_new_pi",), "3.09", "digits", P, 309, ""),
    (V, ("7.c2_seed_from_foreign_Q_stationary",), "1.6e-16", "zero", P, 309, ""),
    (V, ("8.tip_marginal_equals_pi",), "1.4e-16", "zero", P, 310, ""),
    (V, ("9.mstep_normalized_residual",), "2.8e-17", "zero", P, 311, ""),
    (V, ("9.mstep_cone_balance_residual",), "1.6e-14", "zero", P, 311, ""),
    (V, ("9.mstep_cone_mean_rate",), "2.36", "digits", P, 311, ""),
    (V, ("9.mstep_normalized_beats_rescaled_cone_at_fixed_T",), "0.181", "digits", P, 311, ""),
    (V, ("9.mstep_kkt_projected_gradient",), "1.7e-14", "zero", P, 311, ""),
    (V, ("9.compensated_rescale_transition_error",), "5.6e-17", "zero", P, 311, ""),
    (V, ("10.frechet_vs_central_difference",), "1.2e-10", "digits", P, 312, ""),
    (V, ("10.naive_tE_expm_error",), "0.052", "digits", P, 312, ""),
    (V, ("11.wrong_target_population_example", "correct", "gain"), "-3e-15", "zero", P, 313, ""),
    (V, ("11.wrong_target_population_example", "wrong", "gain"), "4.6e-6", "digits", P, 313, ""),
    (V, ("11.wrong_target_population_example", "correct", "cycle_coordinate"), "-1.2e-6", "digits", P, 313, ""),
    (V, ("11.wrong_target_population_example", "wrong", "cycle_coordinate"), "-0.0595", "digits", P, 313, ""),
    (V, ("12.flux_transport_stationary_at_pi_b",), "2.0e-16", "zero", P, 314, ""),
    (V, ("12.flux_transport_mean_rate",), "1.0", "digits", P, 314, ""),
    (C, ("LG", "stationarity_residual_at_shipped_pi"), "2e-17", "zero", P, 316, ""),
    (C, ("LG", "exchangeability_min_max", 0), "0.0035", "digits", P, 316, ""),
    (C, ("LG", "exchangeability_min_max", 1), "10.65", "digits", P, 316, ""),
    (C, ("LG", "within_row_dynamic_range_max"), "6146", "digits", P, 316, ""),
    (C, ("LG", "within_row_dynamic_range_median"), "295", "digits", P, 316, ""),
    (C, ("LG", "global_dynamic_range"), "1.6e4", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "row_sums_max_abs"), "2e-6", "bound", P, 316, "printed as '<= 2e-6'"),
    (C, ("NQ.PFAM", "mean_rate_at_solved_pi"), "1.0000002", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "shipped_pi_vs_solved_pi_max_abs"), "6.9e-7", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "within_row_dynamic_range_max"), "4587", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "global_dynamic_range"), "1.3e4", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "c2_min_ratio_to_row_max"), "2.2e-4", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "c2_fraction_ratios_below_1e-4"), "0", "bound", P, 316, "printed as 'no ratio below 1e-4'"),
    (C, ("NQ.PFAM", "c1_raw_h_max_abs"), "0.94", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "c1_gauge_fixed_h_max_abs"), "0.87", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "fraction_pairs_flux_asymmetry_gt_0.2"), "0.26", "digits", P, 316, "printed as '26%'"),
    (C, ("NQ.PFAM", "c2_retarget_to_LG_pi_residual"), "1.2e-8", "digits", P, 316, ""),
    (C, ("NQ.PFAM", "c2_retarget_max_rel_change_vs_original"), "0.20", "digits", P, 316,
     READ + "'up to 20%' read as 0.20 to two digits"),
    (E, ("E1_fd_accuracy_n20", "n"), 20, "exact", S, 362, ""),
    (E, ("E1_fd_accuracy_n20", "npatterns"), 968, "exact", S, 362, ""),
    (E, ("E1_fd_accuracy_n20", "ratio_chart_rel_err_median"), "8.02e-5", "digits", S, 362, ""),
    (E, ("E1_fd_accuracy_n20", "log_chart_rel_err_median"), "1.44e-4", "digits", S, 362, ""),
    (E, ("E1_fd_accuracy_n20", "ratio_chart_rel_err_max"), "0.00688", "digits", S, 362, ""),
    (E, ("E1_fd_accuracy_n20", "log_chart_rel_err_max"), "0.01969", "digits", S, 362, ""),
]
SWEEP = [("0.01", "3.24e-6", "7.51e-5", 366), ("0.0001", "1.19e-4", "7.72e-5", 367),
         ("1e-06", "9.55e-2", "7.72e-5", 368), ("1e-08", "6.39", "7.71e-5", 369),
         ("0.0", "1.27e-4", "7.71e-5", 370)]
for z, log_err, ratio_err, line in SWEEP:
    REPORTED += [(E, ("E1_fd_accuracy_n20", "log_coordinate_near_zero_sweep", z, "rel_err"), log_err, "digits", S, line, ""),
                 (E, ("E1_fd_accuracy_n20", "log_coordinate_near_zero_sweep", z, "ratio_chart_same_point_rel_err"),
                  ratio_err, "digits", S, line, "")]
E2 = [  # synthesis section 10.3: (logL, stop, evaluations) for each chart and rule
    (382, ("-11343.902245", "TOLX", 278), ("-11303.734515", "gtol", 975), ("-11303.734377", "gtol", 975)),
    (383, ("-12419.902818", "TOLX", 700), ("-12403.624673", "gtol", 822), ("-12403.624636", "gtol", 822)),
    (384, ("-14544.987015", "gtol", 449), ("-14543.502241", "gtol", 666), ("-14543.503265", "gtol", 666)),
    (385, ("-12552.097580", "TOLX", 1426), ("-12552.599324", "gtol", 770), ("-12552.599299", "gtol", 770)),
]
for rep, (line, *variants) in enumerate(E2):
    for variant, (logl, stop, nfev) in zip(("ratio_iqtree_rule", "log_iqtree_rule", "log_scale_aware_rule"), variants):
        REPORTED += [(E, ("E2_bfgs", rep, variant, "logL_final"), logl, "digits", S, line, ""),
                     (E, ("E2_bfgs", rep, variant, "stop"), stop, "exact", S, line, ""),
                     (E, ("E2_bfgs", rep, variant, "nfev"), nfev, "exact", S, line, "")]

NOT_COMPARED = [
    "wall times (P-positive section 12 reports 1.7 s for verify_constrained_nq.py)",
    "the Part 2 slide deck's Backup B, which is not in the repository",
    "synthesis section 10.4, which used modified copies of chart_fd_bfgs_compare.py",
    "the 'D-review' evaluation counts quoted in synthesis section 10.3, which come from another run",
    "P-log Appendix A, whose numbers come from its own Appendix B code (oracle work, CODE_PLAN.md 3.1)",
    "script outputs the documents do not quote",
]


def significant_digits(printed):
    return len(printed.lower().lstrip("+-").split("e")[0].replace(".", "").lstrip("0"))


def reproduced(cls, printed, value):
    """Decision 019's four classes."""
    if cls == "exact":
        return value == printed
    x = float(value)
    if cls == "zero":
        return abs(x) < ZERO
    if cls == "bound":
        return x <= float(printed)
    d = significant_digits(printed) - 1
    return f"{x:.{d}e}" == f"{float(printed):.{d}e}"


def lookup(data, path):
    for key in path:
        data = data[key]
    return data


def compare(results):
    rows = []
    for script, path, printed, cls, doc, line, note in REPORTED:
        try:
            value = lookup(results[script], path)
            ok = reproduced(cls, printed, value)
        except (KeyError, IndexError, TypeError, ValueError) as e:
            value, ok = f"missing ({e!r})", False
        rows.append({"script": script, "path": list(path), "printed": printed, "class": cls,
                     "document": f"{doc}:{line}", "rerun": value, "reproduced": ok, "note": note})
    return rows


def provenance(a, out):
    import numpy
    import scipy
    return {
        "script": SCRIPT,
        "command": sys.argv,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
        "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
        "out": str(out),
        "python": sys.version,
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "numpy_build": numpy.show_config(mode="dicts").get("Build Dependencies"),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
        "environment": {k: v for k, v in os.environ.items()
                        if k == "LANG" or k.startswith(("LC_", "OMP_", "OPENBLAS_", "MKL_"))},
        "fingerprints": FINGERPRINTS,
    }


def markdown(prov, runs, rows):
    failed = [r for r in rows if not r["reproduced"]]
    lines = [
        "# Design-script rerun compared with the documents",
        "",
        f"Generated by `{SCRIPT}` on {prov['started_utc']} at commit `{prov['git_commit']}` "
        f"(clean: {prov['git_clean']}), in `{prov['out']}`. Rule: decision 019 (zero below {ZERO:g}).",
        f"Python {prov['python'].split()[0]}, NumPy {prov['numpy']}, SciPy {prov['scipy']}; {prov['platform']}.",
        "",
        "| Script | Arguments | Added environment | Exit code | Wall time (s) |",
        "|---|---|---|---|---|",
    ]
    for name, r in runs.items():
        env = " ".join(f"{k}={v}" for k, v in r["environment"].items()) or "none"
        lines.append(f"| {name} | `{' '.join(r['arguments'])}` | {env} | {r['exit_code']} | {r['wall_s']} |")
    lines += ["", f"{len(rows)} documented values compared; {len(rows) - len(failed)} reproduced, "
              f"{len(failed)} not reproduced.", "",
              "| Script | Output key | Document:line | Printed | Rerun | Class | Reproduced | Note |",
              "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['script']} | `{'/'.join(map(str, r['path']))}` | {r['document']} | {r['printed']} "
                     f"| {json.dumps(r['rerun'])} | {r['class']} | {'yes' if r['reproduced'] else 'NO'} | {r['note']} |")
    lines += ["", "Not compared:", ""] + [f"- {item}" for item in NOT_COMPARED]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="new directory outside the repository")
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    out = Path(a.out).expanduser().resolve()
    if out.exists():
        sys.exit(f"{out} exists; every rerun goes into a new directory")
    if out.is_relative_to(REPO):
        sys.exit("the rerun directory must be outside the repository")
    prov = provenance(a, out)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    out.mkdir(parents=True)
    for name, digest in FINGERPRINTS.items():
        shutil.copyfile(DESIGN / "scripts" / name, out / name)
        if sha256(out / name) != digest:
            sys.exit(f"{name}: the copy's SHA-256 differs from the design README's {digest}")
    blob = subprocess.run([a.git, "cat-file", "blob", f"{BASE_COMMIT}:model/modelprotein.cpp"],
                          cwd=REPO, capture_output=True, check=True).stdout
    (out / "modelprotein.cpp").write_bytes(blob)
    prov["inputs"] = {f"model/modelprotein.cpp at {BASE_COMMIT}": sha256(out / "modelprotein.cpp")}
    runs, results = {}, {}
    for name, (args, env, output) in RUNS.items():
        stem = Path(name).stem
        print(f"{name}: started {datetime.datetime.now():%H:%M:%S}", flush=True)
        t0 = time.perf_counter()
        with open(out / f"{stem}.stdout", "w", encoding="utf-8") as so, \
                open(out / f"{stem}.stderr", "w", encoding="utf-8") as se:
            code = subprocess.run([sys.executable, name, *args], cwd=out, env={**os.environ, **env},
                                  stdout=so, stderr=se).returncode
        runs[name] = {"arguments": args, "environment": env, "exit_code": code,
                      "wall_s": round(time.perf_counter() - t0, 1)}
        try:
            data = json.loads((out / output).read_text())
        except (OSError, ValueError) as e:
            data, runs[name]["output_error"] = {}, repr(e)
        results[name] = {"(exit code)": code, **data}
        print(f"{name}: exit {code}, {runs[name]['wall_s']} s", flush=True)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    rows = compare(results)
    write(out / "provenance.json", json.dumps(dict(prov, runs=runs), indent=1))
    write(out / "comparison.json", json.dumps(rows, indent=1))
    write(out / "comparison.md", markdown(prov, runs, rows))
    failed = [r for r in rows if not r["reproduced"]]
    print(f"{len(rows)} documented values; {len(failed)} not reproduced; wrote {out / 'comparison.md'}")
    for r in failed:
        print(f"    {r['script']} {'/'.join(map(str, r['path']))}: printed {r['printed']} "
              f"({r['document']}), rerun {r['rerun']}")
    sys.exit(1 if failed or any(r["exit_code"] for r in runs.values()) else 0)


if __name__ == "__main__":
    main()
