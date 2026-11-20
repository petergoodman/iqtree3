"""CLOAK's removal curve on BAliBASE 3.0 (CODE_PLAN.md 3.8, "Laptop work before submission", step 2; decision 045):
the frequency with which CLOAK v0.3.0 removes a residue, binned by r, the share of the other residues in its
reference column identical to it, which sets cell D's masking probability (simulate.rarity_mask).

  cloak.py run --balibase TARBALL --cloak-dir DIR --git git.exe [--out RUNDIR] [--jobs 3] [--threads 4] [--derived DIR]
  cloak.py curve --run RUNDIR --cloak-dir DIR --balibase TARBALL --git git.exe [--into FILE]

run: extracts the 218 full-length BB cases' .tfa, .msf and .xml (tarfile, data filter) into --derived once, the
tarball read-only; then per case, largest first, CLOAK's stratified ensemble of the unaligned .tfa (-align
-stratified, 16 alignments) and cloak_residue_map.py, which says for each residue (label, ungapped position)
whether CLOAK kept it. Resumable: a rerun on the same RUNDIR skips cases with a residue map. curve: joins each
residue to its reference .msf column by label and ungapped position (a case whose letters disagree is excluded and
counted); r = (identical residues in the column - 1) / (non-gap residues in the column - 1), a residue with no other
residue in its column excluded and counted (Peter, 2026-11-19); removal by bin of r (edges 0, 1/51, 0.05, 0.1, 0.2,
..., 0.9, 1; the first bin is r = 0 in simulated columns), with 95% case-bootstrap intervals, the overall removed
fraction and the same for core columns (.xml coreblock 1) beside them; written with provenance lines to FILE
(default sim/masking/cloak_curve.tsv), never over an existing file.
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
import tarfile
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
for _p in ("", "regression", "nqc"):
    sys.path.insert(0, str(HERE.parent / _p))
sys.path.insert(0, str(HERE))
import common  # noqa: E402
from regress import REPO, read, sha256, write  # noqa: E402
from nqcchecks import check  # noqa: E402

SCRIPT = Path(__file__).resolve().relative_to(REPO).as_posix()
BINARY, MAPPER = "muscle-cloak-v0.3.0", "cloak_residue_map.py"
CASE = re.compile(r"bb3_release/(RV\d\d)/(BB\d{5})\.(tfa|msf|xml)")
EDGES = np.array([0.0, 1 / 51, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
BOOTSTRAP, SEED = 2000, 20261119
CURVE = HERE / "masking" / "cloak_curve.tsv"


# ---------------------------------------------------------------- inputs

def extract(tarball, derived):
    """The BB cases' three files, extracted once (a marker file records a finished extraction; a directory
    without it is refused, not overwritten); returns {case: {suffix: path}}."""
    derived = Path(derived)
    marker = derived / ".extracted"
    if derived.exists() and not marker.exists():
        sys.exit(f"{derived} exists without {marker.name}: an interrupted extraction; remove it to extract again")
    if not derived.exists():
        with tarfile.open(tarball) as tar:
            members = [m for m in tar.getmembers() if CASE.fullmatch(m.name)]
            tar.extractall(derived, members=members, filter="data")
        write(marker, f"{Path(tarball).name} {sha256(tarball)}\n")
    cases = {}
    for path in derived.rglob("BB*"):
        m = CASE.fullmatch(path.relative_to(derived).as_posix())
        if m:
            cases.setdefault(m.group(2), {})[m.group(3)] = path
    return {c: f for c, f in sorted(cases.items()) if set(f) == {"tfa", "msf", "xml"}}


def fasta(path):
    out, name = {}, None
    for line in read(path).splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = []
        elif name is not None and line.strip():
            out[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in out.items()}


def msf(path):
    """A GCG MSF alignment: {label: row}, '.' and '-' as gaps, labels in Name: order."""
    text = read(path)
    names = re.findall(r"^\s*Name:\s+(\S+)", text, re.M)
    rows = {n: [] for n in names}
    body = text.split("//", 1)[1]
    for line in body.splitlines():
        parts = line.split()
        if parts and parts[0] in rows:
            rows[parts[0]].append("".join(parts[1:]))
    return {n: "".join(r).upper().replace("-", ".") for n, r in rows.items()}


def coreblocks(path):
    m = re.search(r"<colsco-name>\s*coreblock\s*</colsco-name>.*?<colsco-data>(.*?)</colsco-data>", read(path), re.S)
    return np.array([int(v) for v in m.group(1).split()]) if m else None


# ---------------------------------------------------------------- run

def run_case(case, files, out, binary, mapper, threads):
    d = out / "cases" / case
    if (d / "map" / "residues.tsv").exists() and (d / "case.json").exists():
        return json.loads(read(d / "case.json"))
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    t0 = time.monotonic()
    a = subprocess.run([str(binary), "-align", str(files["tfa"]), "-stratified", "-output", "ens.efa", "-threads",
                        str(threads)], cwd=d, capture_output=True, text=True)
    rec = {"case": case, "align_exit": a.returncode, "align_s": round(time.monotonic() - t0, 1)}
    if a.returncode == 0:
        t1 = time.monotonic()
        m = subprocess.run([sys.executable, str(mapper), str(binary), "ens.efa", "map"], cwd=d, capture_output=True,
                           text=True)
        rec.update({"map_exit": m.returncode, "map_s": round(time.monotonic() - t1, 1), "summary": m.stdout.strip(),
                    "map_error": m.stderr[-1500:] if m.returncode else ""})
        for p in (d / "map").glob("plane*"):
            p.unlink()
    rec["error"] = (a.stderr if a.returncode else "")[-1500:]
    if not (d / "map" / "residues.tsv").exists():
        rec["failed"] = True
    write(d / "case.json", json.dumps(rec, indent=1))
    return rec


def cmd_run(a):
    cloak_dir = Path(a.cloak_dir).expanduser().resolve()
    tarball = Path(a.balibase).expanduser().resolve()
    binary, mapper = cloak_dir / BINARY, cloak_dir / MAPPER
    prov = common.provenance(a, SCRIPT, "cloak-run", "measurement of a setting", jobs=a.jobs, threads=a.threads,
                             binary={"path": str(binary), "sha256": sha256(binary)},
                             mapper={"path": str(mapper), "sha256": sha256(mapper)},
                             balibase={"path": str(tarball), "sha256": sha256(tarball)})
    if a.out and Path(a.out).expanduser().exists():
        out = Path(a.out).expanduser().resolve()
        prov["out"], prov["resumed"] = str(out), True
    else:
        out = common.new_out(a.out, prov, runs="~/iqtree3-runs/sim2")
    cases = extract(tarball, Path(a.derived).expanduser())
    order = sorted(cases, key=lambda c: (-sum(len(s) for s in fasta(cases[c]["tfa"]).values()), c))
    rows, records = [check("BB cases extracted (BAliBASE 3.0's 218)", len(cases) == 218, cases=len(cases))], {}
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for rec in pool.map(lambda c: run_case(c, cases[c], out, binary, mapper, a.threads), order):
            records[rec["case"]] = rec
            print(f"{rec['case']}: align {rec['align_exit']} ({rec['align_s']} s), map {rec.get('map_exit')}", flush=True)
            common.save(out, prov, "CLOAK on BAliBASE 3.0, in progress", rows, {"cases": records})
    failed = sorted(c for c, r in records.items() if r.get("failed"))
    rows.append(check("cases with a residue map", None, done=len(records) - len(failed), failed=failed))
    common.finish(out, prov, "CLOAK on BAliBASE 3.0", rows, {"cases": records})


# ---------------------------------------------------------------- curve

def case_residues(files, residues_tsv):
    """Each residue of a case as (r, kept, core), or the reason the case is excluded; and the count of residues
    alone in their column."""
    ref = msf(files["msf"])
    seqs = fasta(files["tfa"])
    if sorted(ref) != sorted(seqs) or any(ref[n].replace(".", "") != seqs[n] for n in ref):
        return None, "the .msf rows differ from the .tfa sequences", 0
    core = coreblocks(files["xml"])
    names = list(ref)
    cols = np.array([list(ref[n]) for n in names])
    length = cols.shape[1]
    if core is None or len(core) != length:
        return None, "no coreblock scores of the alignment's length", 0
    where = {}
    for i, n in enumerate(names):
        nongap = np.flatnonzero(cols[i] != ".")
        where.update({(n, k + 1): (i, j) for k, j in enumerate(nongap)})
    lines = [line.split("\t") for line in read(residues_tsv).splitlines()]
    head, body = lines[0], lines[1:]
    li, pi_, ri, ki = (head.index(c) for c in ("label", "pos", "residue", "kept"))
    if len(body) != len(where):
        return None, f"{len(body)} mapped residues against {len(where)} in the reference", 0
    out, alone = [], 0
    for row in body:
        key = (row[li].split()[0], int(row[pi_]))
        if key not in where:
            return None, f"residue {key} not in the reference", 0
        i, j = where[key]
        if cols[i, j] != row[ri].upper():
            return None, f"residue {key} is {row[ri]} here and {cols[i, j]} in the reference", 0
        column = cols[:, j]
        present = int(np.sum(column != "."))
        if present < 2:
            alone += 1
            continue
        out.append(((np.sum(column == cols[i, j]) - 1) / (present - 1), row[ki] == "1", core[j] == 1))
    return out, None, alone


def bin_counts(residues):
    """Per bin of r: residues, removed, core residues, core removed."""
    r = np.array([x[0] for x in residues])
    removed = np.array([not x[1] for x in residues])
    core = np.array([x[2] for x in residues])
    k = np.minimum(np.searchsorted(EDGES, r, side="right") - 1, len(EDGES) - 2)
    m = len(EDGES) - 1
    return np.stack([np.bincount(k, minlength=m), np.bincount(k, removed, minlength=m),
                     np.bincount(k, core, minlength=m), np.bincount(k, removed & core, minlength=m)])


def curve_table(per_case, draws=BOOTSTRAP, seed=SEED):
    """Pooled removal per bin with 95% case-bootstrap intervals (cases resampled with replacement, seeded)."""
    counts = np.array(per_case, dtype=float)              # cases x 4 x bins
    total = counts.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)   # a bin no case reaches has no interval
        p = total[1] / total[0]
        rng = np.random.default_rng(seed)
        boot = np.array([(lambda s: s[1] / s[0])(counts[rng.integers(len(counts), size=len(counts))].sum(axis=0))
                         for _ in range(draws)])
        lo, hi = np.nanpercentile(boot, [2.5, 97.5], axis=0)
        core_p = total[3] / total[2]
    return total, p, lo, hi, core_p


def cmd_curve(a):
    run = Path(a.run).expanduser().resolve()
    into = Path(a.into).expanduser().resolve()
    if into.exists():
        sys.exit(f"{into} exists; derived data are never regenerated in place")
    cloak_dir = Path(a.cloak_dir).expanduser().resolve()
    tarball = Path(a.balibase).expanduser().resolve()
    prov = common.provenance(a, SCRIPT, "cloak-curve", "derived data", run=str(run))
    cases = extract(tarball, Path(a.derived).expanduser())
    per_case, excluded, alone, used = [], {}, 0, []
    for case, files in cases.items():
        tsv = run / "cases" / case / "map" / "residues.tsv"
        if not tsv.exists():
            excluded[case] = "no residue map"
            continue
        res, why, n_alone = case_residues(files, tsv)
        if res is None:
            excluded[case] = why
            continue
        per_case.append(bin_counts(res))
        alone += n_alone
        used.append(case)
    total, p, lo, hi, core_p = curve_table(per_case)
    header = {"content": "CLOAK v0.3.0's removal frequency on BAliBASE 3.0's full-length BB cases by r, the share "
                         "of the other residues in a residue's reference column identical to it (decision 045)",
              "script": SCRIPT, "git_commit": prov["git_commit"], "git_clean": prov["git_clean"],
              "generated_utc": prov["started_utc"], "run": str(run),
              "cloak_sha256": sha256(cloak_dir / BINARY), "mapper_sha256": sha256(cloak_dir / MAPPER),
              "balibase_sha256": sha256(tarball), "cases_used": len(used), "cases_excluded": json.dumps(excluded),
              "residues_alone_excluded": alone, "removed_overall": float(total[1].sum() / total[0].sum()),
              "removed_core_overall": float(total[3].sum() / total[2].sum()),
              "bootstrap": f"{BOOTSTRAP} case resamples, seed {SEED}, 95% percentile intervals"}
    lines = [f"# {k}: {v}" for k, v in header.items()]
    lines.append("\t".join(("bin_lo", "bin_hi", "residues", "removed", "p", "ci_lo", "ci_hi", "core_residues",
                            "core_removed", "core_p")))
    for b in range(len(EDGES) - 1):
        lines.append("\t".join([f"{EDGES[b]:.17g}", f"{EDGES[b + 1]:.17g}", str(int(total[0, b])), str(int(total[1, b])),
                                f"{p[b]:.6g}", f"{lo[b]:.6g}", f"{hi[b]:.6g}", str(int(total[2, b])),
                                str(int(total[3, b])), f"{core_p[b]:.6g}"]))
    into.parent.mkdir(parents=True, exist_ok=True)
    write(into, "\n".join(lines) + "\n")
    print(f"wrote {into}: {len(used)} cases, {len(excluded)} excluded")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="step", required=True)
    r = sub.add_parser("run")
    r.add_argument("--out", help="the run directory: new, or an earlier run's to resume (default under "
                                 "~/iqtree3-runs/sim2/)")
    r.add_argument("--jobs", type=int, default=3)
    r.add_argument("--threads", type=int, default=4)
    c = sub.add_parser("curve")
    c.add_argument("--run", required=True)
    c.add_argument("--into", default=str(CURVE))
    for p in (r, c):
        p.add_argument("--balibase", required=True, help="BAliBASE_R1-5.tar.gz, read-only")
        p.add_argument("--cloak-dir", required=True, help="the folder holding muscle-cloak-v0.3.0 and its mapper")
        p.add_argument("--derived", default="~/iqtree3-data/derived/balibase")
        p.add_argument("--git", default="git")
        p.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args()
    {"run": cmd_run, "curve": cmd_curve}[a.step](a)


if __name__ == "__main__":
    main()
