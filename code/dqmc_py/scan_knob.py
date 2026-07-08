#!/usr/bin/env python3
"""Scan an altermagnet knob (tam or t1) at fixed mu/beta and emit the k-resolved d-wave
pair susceptibility reductions (maxk/k0/r0) + <sign> + density. CSV to stdout. Used to
build the finite-T-DQMC vs T=0-CPQMC cross-comparison (Phase 3) on a small lattice."""
import os, sys, argparse
sys.path.insert(0, os.path.dirname(__file__))
from dqmc import DQMC


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=4); ap.add_argument("--ly", type=int, default=4)
    ap.add_argument("--U", type=float, default=4.0); ap.add_argument("--mu", type=float, default=1.5)
    ap.add_argument("--beta", type=float, default=4.0); ap.add_argument("--dt", type=float, default=0.05)
    ap.add_argument("--knob", choices=["tam", "t1", "tp"], default="tam")
    ap.add_argument("--vals", type=str, default="0.0,0.1,0.2,0.3")
    ap.add_argument("--nwarm", type=int, default=200); ap.add_argument("--nmeas", type=int, default=1500)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    vals = [float(v) for v in a.vals.split(",")]
    print(f"# DQMC scan {a.lx}x{a.ly} U={a.U} mu={a.mu} beta={a.beta} dt={a.dt} knob={a.knob}")
    print("knob,dens,sign,sc_maxk,sc_k0,sc_r0,cr_maxk,cr_k0,cr_r0")
    for v in vals:
        kw = dict(tam=0.0, t1=0.0, tp=0.0); kw[a.knob] = v
        q = DQMC(a.lx, a.ly, a.U, a.mu, a.beta, a.dt, seed=a.seed, **kw)
        r = q.run(a.nwarm, a.nmeas, chi=True, kres=True)
        sc = r["susc_d"]; cr = r["corr_d"]
        print(f"{v:.3f},{r['dens']:.4f},{r['sign']:.4f},"
              f"{sc['maxk']:.4f},{sc['k0']:.4f},{sc['r0']:.4f},"
              f"{cr['maxk']:.4f},{cr['k0']:.4f},{cr['r0']:.4f}", flush=True)


if __name__ == "__main__":
    main()
