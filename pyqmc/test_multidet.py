#!/usr/bin/env python3
"""Phase 6 WS3 variant: multi-determinant trial. Compares fixed / adaptive /
sample1 / multidet(k) against ED on the single-band Hubbard model (mixed-estimator
energy), and scans the number of determinants k for the multidet trial.

Hypothesis: a weighted superposition of k walkers (sampled by weight) is the
principled middle ground between the single-determinant natural-orbital average
(low variance, but capped bias) and the single-walker trial (huge variance). As k
grows the multidet trial should LOWER both the bias (it converges toward the
correlated state) and the variance (relative to one walker).

    source tools/env.sh
    python pyqmc/test_multidet.py --lx 4 --ly 4 --nup 2 --ndn 2 --U 4 --ks 4 8 16
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


def run(label, mode, a, e0, k=8):
    from cpqmc import CPMC
    q = CPMC(a.lx, a.ly, a.nup, a.ndn, t=a.t, U=a.U, dt=a.dt, nwalkers=a.nw,
             seed=a.seed, trial=mode, trial_every=a.trial_every,
             trial_mix=a.trial_mix, trial_k=k)
    e, err = q.run(nequil=a.nequil, nmeas=a.nmeas)
    z = abs(e - e0) / err if err > 0 else float("nan")
    print(f"  {label:14s} E = {e:9.4f} +/- {err:6.4f}   dE(ED) = {e - e0:+7.4f}   z = {z:5.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--t", type=float, default=1.0)
    ap.add_argument("--dt", type=float, default=0.01); ap.add_argument("--nw", type=int, default=400)
    ap.add_argument("--nequil", type=int, default=200); ap.add_argument("--nmeas", type=int, default=300)
    ap.add_argument("--trial-every", type=int, default=10); ap.add_argument("--trial-mix", type=float, default=0.5)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--ks", type=int, nargs="+", default=[4, 8, 16])
    a = ap.parse_args()

    e0 = ed_energy(a.lx, a.ly, a.nup, a.ndn, a.t, a.U)
    print(f"# single-band {a.lx}x{a.ly}, nup={a.nup} ndn={a.ndn}, U={a.U} (mixed estimator)")
    print(f"ED ground state E0 = {e0:.6f}  (nw={a.nw})\n")
    run("fixed", "fixed", a, e0)
    run("adaptive", "adaptive", a, e0)
    run("sample1", "sample1", a, e0)
    for k in a.ks:
        run(f"multidet k={k}", "multidet", a, e0, k=k)


if __name__ == "__main__":
    main()
