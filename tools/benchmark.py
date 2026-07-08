#!/usr/bin/env python3
"""Speed benchmark: Fortran CPQMC (code/src, src-small build) vs the Python port
(pyqmc) on the same small two-orbital problem. Reports wall time and throughput
(walker-updates / second), normalizing for walker count and step count.

    source tools/env.sh
    python tools/benchmark.py            # default 2x2, 4+4, uxx=2

Note: the Fortran build is MPI + compiled with NWLKRS=1000 (compile-time); the
Python port is interpreted numpy with brute-force determinant ratios. The
throughput ratio is the apples-to-apples comparison.
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NWLKRS_FORT = 1000          # src-small compile-time walker count


def run_fortran(lx, ly, nup, ndn, uxx, uxy, v, nblk, nblkstps, ranks):
    """Run src-small; returns (wall_sec, walker_steps). Assumes it is built."""
    indat = f"""0,0,{nblk},{nblkstps}
10,10,10,10,10
0.010000,-50.000000
in.dat
-1.000000,-1.00000,-1.000000,-1.00000,0.1
{uxx}.0,{uxy}.0,{v}.0
0.000000,0.000000
1.000000
1.000000
"""
    tmp = REPO / "results" / "_bench_in.dat"
    tmp.parent.mkdir(exist_ok=True); tmp.write_text(indat)
    r = subprocess.run(["bash", str(REPO / "tools/run.sh"), "src-small", str(tmp), str(ranks)],
                       capture_output=True, text=True, timeout=900)
    # newest run dir
    runs = sorted((REPO / "results/src-small").glob("*/"), key=os.path.getmtime)
    meta = json.loads((runs[-1] / "run_meta.json").read_text())
    return float(meta["wall_sec"]), NWLKRS_FORT * nblk * nblkstps


def run_python(lx, ly, nup, ndn, uxx, uxy, v, nw, nsteps):
    sys.path.insert(0, str(REPO / "pyqmc")); sys.path.insert(0, str(REPO / "ed"))
    from cpqmc import CPMC
    from altermagnet_ed import build_hopping
    K = build_hopping(lx, ly, -1, -1, -1, -1)
    qmc = CPMC(lx, ly, nup, ndn, U=uxx, dt=0.01, nwalkers=nw, K=K, uxy=uxy, v=v)
    t0 = time.time()
    qmc.run(nequil=0, nmeas=nsteps, ortho=10, pc=10, meas_every=nsteps)
    wall = time.time() - t0
    return wall, nw * nsteps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=4); ap.add_argument("--ndn", type=int, default=4)
    ap.add_argument("--uxx", type=float, default=2); ap.add_argument("--uxy", type=float, default=0)
    ap.add_argument("--v", type=float, default=0)
    ap.add_argument("--fort-nblk", type=int, default=2); ap.add_argument("--fort-nblkstps", type=int, default=20)
    ap.add_argument("--py-nw", type=int, default=100); ap.add_argument("--py-nsteps", type=int, default=20)
    ap.add_argument("--ranks", type=int, default=2)
    ap.add_argument("-o", "--out")
    a = ap.parse_args()

    fw, fsteps = run_fortran(a.lx, a.ly, a.nup, a.ndn, int(a.uxx), int(a.uxy), int(a.v),
                             a.fort_nblk, a.fort_nblkstps, a.ranks)
    pw, psteps = run_python(a.lx, a.ly, a.nup, a.ndn, a.uxx, a.uxy, a.v, a.py_nw, a.py_nsteps)
    fort_tp = fsteps / fw if fw else 0
    py_tp = psteps / pw if pw else 0
    res = {
        "config": {"lx": a.lx, "ly": a.ly, "nup": a.nup, "ndn": a.ndn,
                   "uxx": a.uxx, "uxy": a.uxy, "v": a.v},
        "fortran": {"wall_sec": round(fw, 3), "walker_steps": fsteps,
                    "throughput": round(fort_tp, 1), "walkers": NWLKRS_FORT, "ranks": a.ranks},
        "python": {"wall_sec": round(pw, 3), "walker_steps": psteps,
                   "throughput": round(py_tp, 1), "walkers": a.py_nw},
        "speedup_fortran_over_python": round(fort_tp / py_tp, 1) if py_tp else None,
    }
    out = json.dumps(res, indent=2)
    if a.out:
        Path(a.out).write_text(out)
    print(out)
    print(f"\nThroughput (walker-updates/s): Fortran {fort_tp:,.0f}  |  Python {py_tp:,.0f}"
          f"  ->  Fortran is {res['speedup_fortran_over_python']}x faster")


if __name__ == "__main__":
    main()
