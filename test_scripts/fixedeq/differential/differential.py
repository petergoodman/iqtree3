#!/usr/bin/env python3
"""The oracle's likelihood against an IQ-TREE binary (docs/agent/CODE_PLAN.md 3.3; S0 probe (d)).

  python differential.py --binary ~/iqtree3-baseline/iqtree3 --out DIR [--git git.exe]

Sixteen cases: two Q matrices written at 17 significant digits (NQ.pfam and a random
non-reversible Q), four fixed rate models, and two data sets with rooted trees taken from the
regression baseline. Each case runs
`iqtree3 -s ALN -m QFILE[+...] -te TREE --show-lh -seed 1 -T 1` in its own directory, and the
"Initial log-likelihood" it prints at precision 17 is compared with the oracle's under decision
021: |dlnL| <= 1e-8 |lnL|. Exit status 1 if a case fails or prints fewer than 10 significant
digits.
"""
import argparse
import datetime
import json
import platform
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "regression"))
from oracle import builtin, cases, likelihood as lk, readers  # noqa: E402
from regress import REPO, capture, read, sha256, write  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
TOL = 1e-8
SEED = 20261004
ALPHA, PINV = 0.5, 0.2
RATE_MODELS = {"none": "", "G4": f"+G4{{{ALPHA}}}", "I": f"+I{{{PINV}}}", "IG4": f"+I{{{PINV}}}+G4{{{ALPHA}}}"}
BASELINE = REPO / "test_scripts" / "fixedeq" / "regression" / "baseline" / "baseline.json"
INITIAL = re.compile(r"^1\. Initial log-likelihood: (\S+)", re.M)


def significant_digits(printed):
    return len(printed.lstrip("+-").replace(".", "").lstrip("0"))


def matrices():
    Q_nq, pi_nq = builtin.full_matrix("NQ.PFAM")
    rng = np.random.default_rng(SEED)
    pi = cases.random_target(20, rng)
    return {"nqpfam": (Q_nq, pi_nq), "random": (cases.generator_from_flux(cases.cycle_flux(20, rng), pi), pi)}


def write_matrix(path, Q, pi):
    rows = [" ".join(f"{v:.17g}" for v in row) for row in Q] + [" ".join(f"{v:.17g}" for v in pi)]
    write(path, "\n".join(rows) + "\n")


def read_matrix(path):
    return np.array(read(path).split(), dtype=float)[:400].reshape(20, 20)


def datasets(out):
    """Alignment path and rooted tree path for each data set, written into out where needed."""
    trees = json.loads(read(BASELINE))["runs"]
    turtle = readers.read_fasta(REPO / "test_scripts" / "test_data" / "turtle_aa.fasta")
    first = next(iter(readers.read_charsets(REPO / "test_scripts" / "test_data" / "turtle_aa.nex").values()))
    write(out / "turtle_part1.fasta", "".join(f">{n}\n{''.join(s[i] for i in first)}\n" for n, s in turtle.items()))
    sets = {"aa_example": (REPO / "example" / "aa_example.phy", "run04"),
            "turtle_part1": (out / "turtle_part1.fasta", "run08")}
    result = {}
    for name, (aln, run) in sets.items():
        tree = out / f"{run}.treefile"
        write(tree, trees[run]["items"]["(treefile)"]["text"])
        if len(readers.read_newick(tree).children) != 2:
            sys.exit(f"{tree}: the baseline tree is not rooted")
        result[name] = (aln, tree)
    return result


def oracle_lnl(aln_path, tree_path, qmat, rate_model):
    aln = readers.read_phylip(aln_path) if aln_path.suffix == ".phy" else readers.read_fasta(aln_path)
    Q, root = lk.file_model(read_matrix(qmat))
    return lk.alignment_log_likelihood(aln, readers.read_newick(tree_path), Q, root,
                                       shape=ALPHA if "G4" in rate_model else None,
                                       p_invar=PINV if rate_model.startswith("I") else 0.0)


def provenance(a, binary, out):
    import numpy
    import scipy
    return {
        "script": SCRIPT, "command": sys.argv,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "git_commit": capture(a.git, "rev-parse", "HEAD", cwd=REPO),
        "git_clean": capture(a.git, "status", "--porcelain", cwd=REPO) == "",
        "binary": {"path": str(binary), "sha256": sha256(binary)},
        "out": str(out), "python": sys.version, "numpy": numpy.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(),
        "settings": {"tolerance": TOL, "seed": SEED, "alpha": ALPHA, "p_invar": PINV},
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--binary", required=True)
    ap.add_argument("--out", required=True, help="new directory outside the repository")
    ap.add_argument("--git", default="git", help="git executable (git.exe under WSL)")
    ap.add_argument("--allow-dirty", action="store_true")
    ap.add_argument("--timeout", type=float, default=600)
    a = ap.parse_args()
    out = Path(a.out).expanduser().resolve()
    if out.exists() or out.is_relative_to(REPO):
        sys.exit(f"{out}: needs a new directory outside the repository")
    binary = Path(a.binary).expanduser().resolve()
    prov = provenance(a, binary, out)
    if not prov["git_clean"] and not a.allow_dirty:
        sys.exit(f"the repository is not clean ({prov['git_commit']}); commit first or pass --allow-dirty")
    out.mkdir(parents=True)
    qmats = {}
    for name, (Q, pi) in matrices().items():
        qmats[name] = out / f"{name}.qmat"
        write_matrix(qmats[name], Q, pi)
    data = datasets(out)
    prov["inputs"] = {str(p.relative_to(REPO) if p.is_relative_to(REPO) else p.name): sha256(p)
                      for p in [*qmats.values(), *(x for pair in data.values() for x in pair)]}
    rows = []
    for dname, (aln, tree) in data.items():
        for mname, qmat in qmats.items():
            for rname, suffix in RATE_MODELS.items():
                case = f"{dname}-{mname}-{rname}"
                cdir = out / case
                cdir.mkdir()
                shutil.copyfile(qmat, cdir / qmat.name)
                argv = [str(binary), "-s", str(aln), "-m", qmat.name + suffix, "-te", str(tree),
                        "--show-lh", "-seed", "1", "-T", "1", "--prefix", case]
                t0 = time.perf_counter()
                with open(cdir / "stdout.txt", "w", encoding="utf-8") as fh:
                    code = subprocess.run(argv, cwd=cdir, stdout=fh, stderr=subprocess.STDOUT, timeout=a.timeout).returncode
                m = INITIAL.search(read(cdir / "stdout.txt"))
                row = {"case": case, "model": qmat.name + suffix, "exit_code": code,
                       "iqtree": m.group(1) if m else None, "wall_s": round(time.perf_counter() - t0, 1)}
                row["oracle"] = oracle_lnl(aln, tree, qmat, rname)
                if m:
                    iq = float(m.group(1))
                    row.update(digits=significant_digits(m.group(1)), abs_diff=abs(row["oracle"] - iq),
                               rel_diff=abs(row["oracle"] - iq) / abs(iq))
                    row["passed"] = code == 0 and row["digits"] >= 10 and row["abs_diff"] <= TOL * abs(iq)
                else:
                    row["passed"] = False
                rows.append(row)
                print(f"{case}: IQ-TREE {row['iqtree']}, oracle {row['oracle']!r}, "
                      f"relative {row.get('rel_diff', float('nan')):.2e}, {'pass' if row['passed'] else 'FAIL'}", flush=True)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    write(out / "results.json", json.dumps({"provenance": prov, "cases": rows}, indent=1))
    lines = ["# Oracle against IQ-TREE", "",
             f"Generated by `{SCRIPT}` on {prov['started_utc']} at commit `{prov['git_commit']}` "
             f"(clean: {prov['git_clean']}); binary SHA-256 `{prov['binary']['sha256']}`.",
             f"Rule: decision 021, |dlnL| <= {TOL:g} |lnL|; alpha {ALPHA}, p_invar {PINV}; "
             f"NumPy {prov['numpy']}, SciPy {prov['scipy']}.", "",
             "| Case | Exit | IQ-TREE initial lnL | Oracle lnL | Relative difference | Digits | Passed |",
             "|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['case']} | {r['exit_code']} | {r['iqtree']} | {r['oracle']!r} | "
                     f"{r.get('rel_diff', float('nan')):.2e} | {r.get('digits')} | {'yes' if r['passed'] else 'NO'} |")
    write(out / "results.md", "\n".join(lines) + "\n")
    failed = [r["case"] for r in rows if not r["passed"]]
    print(f"{len(rows) - len(failed)} of {len(rows)} cases passed; wrote {out / 'results.md'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
