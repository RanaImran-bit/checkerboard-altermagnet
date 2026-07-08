#!/usr/bin/env python3
"""Phase 6 WS3 follow-up: how does the bias vs ED scale with the number of walkers,
for the FIXED vs the ADAPTIVE trial?

Hypothesis: the constrained-path bias is systematic (NOT statistical), so for the
FIXED trial dE(ED) is independent of the walker count (only the error bar shrinks
~1/sqrt(Nw)). For the ADAPTIVE trial the trial is rebuilt from the ensemble 1-RDM,
whose accuracy improves with more walkers -> a better trial -> the bias itself
should coherently DECREASE with the walker count.

    source tools/env.sh
    python pyqmc/scan_walkers.py --lx 4 --ly 4 --nup 2 --ndn 2 --U 4 \
        --nws 100 250 500 1000
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def ed_energy(lx, ly, nup, ndn, t, U):
    from hubbard_ed import build
    _, H, _, _ = build(lx, ly, nup, ndn, t, U)
    return float(np.linalg.eigvalsh(H.toarray())[0])


def one(mode, nw, a, e0):
    from cpqmc import CPMC
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=nw,
             seed=a.seed, trial=mode, trial_every=a.trial_every, trial_mix=a.trial_mix)
    e, err = q.run(nequil=a.nequil, nmeas=a.nmeas)
    return e, err, e - e0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nmeas", type=int, default=400)
    ap.add_argument("--trial-every", type=int, default=10); ap.add_argument("--trial-mix", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nws", type=int, nargs="+", default=[100, 250, 500, 1000])
    ap.add_argument("-o", "--out", help="write JSON of the scan")
    a = ap.parse_args()

    e0 = ed_energy(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}; ED E0 = {e0:.6f}")
    print(f"{'Nw':>6} | {'fixed E':>10} {'+/-':>7} {'dE(ED)':>8} | "
          f"{'adapt E':>10} {'+/-':>7} {'dE(ED)':>8}")
    scan = {"e0": e0, "nws": a.nws, "fixed": [], "adaptive": []}
    for nw in a.nws:
        ef, errf, df = one("fixed", nw, a, e0)
        ea, erra, da = one("adaptive", nw, a, e0)
        scan["fixed"].append({"nw": nw, "E": ef, "err": errf, "dE": df})
        scan["adaptive"].append({"nw": nw, "E": ea, "err": erra, "dE": da})
        print(f"{nw:>6} | {ef:>10.4f} {errf:>7.4f} {df:>+8.4f} | "
              f"{ea:>10.4f} {erra:>7.4f} {da:>+8.4f}")
    if a.out:
        import json
        with open(a.out, "w") as f:
            json.dump(scan, f, indent=2)
        print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
