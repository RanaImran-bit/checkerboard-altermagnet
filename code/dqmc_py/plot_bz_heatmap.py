#!/usr/bin/env python3
"""Brillouin-zone heatmaps of the d-wave pair SUSCEPTIBILITY P(q) (pair center-of-mass
momentum) from the finite-T DQMC, 6x6 U=4 doped (beta=2), for three knob settings:
baseline, anisotropic t' on, altermagnet tam on. Shows WHERE in q-space the d-wave weight
sits (the 'maxk' peak) and how t'/tam move it. P(q)=FFT2 of S(R) (reduce_mat)."""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from dqmc import DQMC, _shift_index, reduce_mat

L = 6; beta = 2.0; mu = 1.0; U = 4.0; dt = 0.0625
settings = [(0.0, 0.0, r"baseline ($t'{=}t_{am}{=}0$)"),
            (0.0, 0.4, r"$t'{=}0.4$ (anisotropic NNN)"),
            (0.4, 0.0, r"$t_{am}{=}0.4$ (altermagnet)")]
shift = _shift_index(L, L)
maps = []
for tam, t1, lab in settings:
    q = DQMC(L, L, U, mu, beta, dt, tam=tam, t1=t1, seed=7)
    r = q.run(120, 500, chi=True, kres=True)
    red = r["susc_d"]
    Pq = np.fft.fftshift(red["Pq"])         # center q=0
    maps.append((Pq, r["dens"], r["sign"], lab))
    print(f"{lab}: n={r['dens']:.3f} sign={r['sign']:.3f} maxk={red['maxk']:.2f} k0={red['k0']:.2f}", flush=True)

vmax = max(m[0].max() for m in maps)
fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
ext = [-1, 1, -1, 1]    # q in units of pi
for a, (Pq, dens, sg, lab) in zip(ax, maps):
    im = a.imshow(Pq, origin="lower", extent=ext, cmap="inferno", vmin=0, vmax=vmax,
                  interpolation="bilinear", aspect="equal")
    a.set_title(lab + f"\n$n{{=}}{dens:.2f}$, $\\langle$sign$\\rangle{{=}}{sg:.2f}$", fontsize=9.5)
    a.set_xlabel(r"$q_x/\pi$"); a.set_ylabel(r"$q_y/\pi$")
    a.plot(0, 0, "c+", ms=8, mew=1.5)
    fig.colorbar(im, ax=a, fraction=0.046, pad=0.04)
fig.suptitle(r"Finite-$T$ DQMC: $d$-wave pair susceptibility $P(q)$ over the Brillouin zone "
             r"(6$\times$6, $U{=}4$, $\beta{=}2$)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = os.path.join(os.path.dirname(__file__), "..", "..", "docs", "dqmc_bz_susc_dwave.png")
fig.savefig(os.path.normpath(out), dpi=130); print("wrote", os.path.normpath(out))
