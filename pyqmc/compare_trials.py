#!/usr/bin/env python3
"""Head-to-head: fixed vs adaptive vs CASCI trial, ACCURACY and TIME, same config.

For one cluster/coupling, run all three trials at matched statistics, report the
energy, the bias vs ED, the bias reduction relative to fixed, and the wall time of
the (steady-state) measurement run. The CASCI setup (seed fixed run -> 1-RDM ->
natural orbitals -> CAS diag) is one-time and reported separately; the per-run time
is the fair cost comparison.

    python pyqmc/compare_trials.py --lx 4 --ly 4 --nup 2 --ndn 2 --U 4 --ed -11.5303 --nact 4
"""
from __future__ import annotations
import argparse, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))
from cpqmc import CPMC, square_hopping
from casci import ensemble_rdm, natural_orbitals, casci_trial


def timed_run(q, nequil, nmeas):
    t0 = time.time(); e, err = q.run(nequil=nequil, nmeas=nmeas); wall = time.time() - t0
    return e, err, wall, wall / (nequil + nmeas) * 1e3      # ms/step


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--nup", type=int, default=2); ap.add_argument("--ndn", type=int, default=2)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--ed", type=float, required=True)
    ap.add_argument("--ncore", type=int, default=0); ap.add_argument("--nact", type=int, default=4)
    ap.add_argument("--nw", type=int, default=240); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nequil", type=int, default=140); ap.add_argument("--nmeas", type=int, default=220)
    ap.add_argument("--trial-every", type=int, default=10); ap.add_argument("--trial-mix", type=float, default=0.5)
    a = ap.parse_args()
    n = a.lx * a.ly; K = square_hopping(a.lx, a.ly, 1.0)
    base = dict(lx=a.lx, ly=a.ly, nup=a.nup, ndn=a.ndn, t=1.0, U=a.U, dt=0.01, nwalkers=a.nw, seed=a.seed)
    print(f"# {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} U={a.U}  ED={a.ed:.4f}  "
          f"CAS(core {a.ncore}, active {a.nact})  nw={a.nw} steps={a.nequil}+{a.nmeas}")

    # fixed
    ef, errf, wf, msf = timed_run(CPMC(**base, trial="fixed"), a.nequil, a.nmeas)
    # adaptive
    ea, erra, wa, msa = timed_run(CPMC(**base, trial="adaptive", trial_every=a.trial_every,
                                       trial_mix=a.trial_mix), a.nequil, a.nmeas)
    # CASCI: one-time setup (seed RDM -> natural orbitals -> CAS), then timed run
    t0 = time.time()
    qs = CPMC(**base, trial="fixed"); qs.run(nequil=a.nequil, nmeas=0)
    ru, rd = ensemble_rdm(qs); occ, W = natural_orbitals(ru + rd)
    _, ups, dns, coef, ncfg = casci_trial(K, a.U, W, a.nup, a.ndn, a.ncore, a.nact)
    qc = CPMC(**base, trial="multidet"); qc.trial.set_multidet(ups, dns, coef); qc.trial._refresh_olp(qc.walkers)
    if qc.walkers.olp[0] < 0: qc.trial.coef = -qc.trial.coef; qc.trial._refresh_olp(qc.walkers)
    setup = time.time() - t0
    ec, errc, wc, msc = timed_run(qc, a.nequil, a.nmeas)

    def row(name, e, err, ms, w, extra=""):
        dE = e - a.ed
        cut = (1 - abs(dE) / abs(ef - a.ed)) * 100 if abs(ef - a.ed) > 0 else 0
        return (f"{name:>10}  E={e:8.4f}+/-{err:.4f}  dE={dE:+.4f}  "
                f"bias_cut={cut:5.1f}%   {ms:7.1f} ms/step  {w/wf:5.2f}x time {extra}")
    print(f"\n{'trial':>10}  {'energy':>16}  {'dE(ED)':>9}  {'accuracy':>13}   {'speed':>20}")
    print(row("fixed", ef, errf, msf, wf))
    print(row("adaptive", ea, erra, msa, wa))
    print(row("CASCI", ec, errc, msc, wc, f"(k={ncfg} dets; +{setup:.1f}s one-time setup)"))


if __name__ == "__main__":
    main()
