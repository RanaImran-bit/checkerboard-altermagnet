#!/usr/bin/env python3
"""2D (t1, tam) map of the d-wave VERTEX and FULL susceptibility/correlation (all 4
reductions maxk/k0/r0/rgt) from the validated finite-T DQMC, at the appeal filling.
Output CSV for the appeal heatmaps."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from dqmc import DQMC
lx = ly = int(sys.argv[1]) if len(sys.argv) > 1 else 6
mu = float(sys.argv[2]) if len(sys.argv) > 2 else 1.8
beta = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
vals = [0.0, 0.1, 0.2, 0.3, 0.4]
cols = ["maxk", "k0", "r0", "rgt"]
hdr = ["tam", "t1", "dens", "sign"]
for f in ("suscV", "suscF", "corrV", "corrF"):
    hdr += [f"{f}_{c}" for c in cols]
print(",".join(hdr))
for tam in vals:
    for t1 in vals:
        q = DQMC(lx, ly, 4.0, mu, beta, 0.0625, tam=tam, t1=t1, seed=11)
        r = q.run(120, 600, chi=True, kres=True)
        row = [f"{tam:.2f}", f"{t1:.2f}", f"{r['dens']:.4f}", f"{r['sign']:.4f}"]
        for key in ("suscV_d", "susc_d", "corrV_d", "corr_d"):
            row += [f"{r[key][c]:.4f}" for c in cols]
        print(",".join(row), flush=True)
