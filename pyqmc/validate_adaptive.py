#!/usr/bin/env python3
"""Phase 6 WS3 test: does the adaptive (self-consistent natural-orbital) trial
reduce the constrained-path bias vs ED, relative to the fixed free-electron trial,
in a hard / sign-problematic regime?

Compares fixed vs adaptive trials (mixed + back-propagated energy) against exact
ED on the single-band Hubbard model.

    source tools/env.sh
    python pyqmc/validate_adaptive.py --lx 4 --ly 4 --nup 2 --ndn 2 --U 4
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def ed_energy(lx, ly, nup, ndn, t, U):
    from hubbard_ed import build
    basis, H, K, V = build(lx, ly, nup, ndn, t, U)
    w = np.linalg.eigvalsh(H.toarray())
    return float(w[0])


def run(label, mode, a, e0, **kw):
    from cpqmc import CPMC
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw,
             seed=a.seed, trial=mode, trial_every=a.trial_every, trial_mix=a.trial_mix)
    em, errm = q.run(nequil=a.nequil, nmeas=a.nmeas)
    q2 = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw,
              seed=a.seed, trial=mode, trial_every=a.trial_every, trial_mix=a.trial_mix)
    eb, errb = q2.run_bp(nequil=a.nequil, nblocks=max(a.nmeas // a.bp, 12), bp=a.bp)
    for est, e, err in (("mixed", em, errm), (f"bp{a.bp}", eb, errb)):
        z = abs(e - e0) / err if err > 0 else float("nan")
        print(f"  {label:9s} {est:7s} E = {e:9.4f} +/- {err:6.4f}   "
              f"dE(ED) = {e - e0:+7.4f}   z = {z:5.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=300)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nmeas", type=int, default=300)
    ap.add_argument("--bp", type=int, default=14); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--trial-every", type=int, default=10)
    ap.add_argument("--trial-mix", type=float, default=0.5)
    ap.add_argument("--modes", nargs="+", default=["fixed", "adaptive", "sample1"])
    a = ap.parse_args()

    e0 = ed_energy(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U}, t={a.t}")
    print(f"ED ground state E0 = {e0:.6f}\n")
    for mode in a.modes:
        run(mode, mode, a, e0)


if __name__ == "__main__":
    main()
