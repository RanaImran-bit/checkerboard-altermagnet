#!/usr/bin/env python3
"""PILOT: momentum-resolved magnetic susceptibility chi_zz(q) of the COMBINED altermagnet
model (tam=tA NN anisotropy + t1=t' NNN spin-dependent hopping), finite-T DQMC, 4x4 U=4.
Four corners of the (tam, t1) plane at half filling -> where does the dominant chi_zz(q)
peak sit, and does the combination differ from either knob alone? Writes one CSV per
point (full q-grid + summary) so interrupted scans resume; --plot assembles the figure."""
from __future__ import annotations
import os, sys, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from spin_susc import SpinDQMC

OUT = os.environ.get("PILOT_OUT", os.path.join(os.path.dirname(__file__), "..", "..", "results", "chi_spin_pilot"))
POINTS = [(0.0, 0.0), (0.3, 0.0), (0.0, 0.3), (0.3, 0.3)]   # (tam, t1)


def run_point(tam, t1, L=4, U=4.0, mu=2.0, beta=3.0, dt=0.0625, nwarm=200, nmeas=600, seed=11):
    q = SpinDQMC(L, L, U, mu, beta, dt, tam=tam, t1=t1, seed=seed)
    r = q.run_spin(nwarm, nmeas)
    os.makedirs(OUT, exist_ok=True)
    f = os.path.join(OUT, f"chi_tam{tam:g}_t1{t1:g}.csv")
    with open(f, "w") as fh:
        fh.write(f"# 4x4 U={U} mu={mu} beta={beta} dt={dt} tam={tam} t1={t1} "
                 f"nmeas={nmeas} dens={r['dens']:.5f} sign={r['sign']:.4f}\n")
        fh.write(f"# chi_q0={r['chi_q0']:.6f} chi_max={r['chi_max']:.6f} S_max={r['S_max']:.6f}\n")
        for tag, grid in (("chi", r["chi_q"]), ("S", r["S_q"])):
            for kx in range(L):
                fh.write(tag + "," + str(kx) + "," + ",".join(f"{v:.6f}" for v in grid[kx]) + "\n")
    print(f"done tam={tam} t1={t1}: dens={r['dens']:.4f} sign={r['sign']:.3f} "
          f"chi_max={r['chi_max']:.4f} chi_q0={r['chi_q0']:.4f}", flush=True)


def plot(L=4):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(16, 4))
    for a, (tam, t1) in zip(ax, POINTS):
        f = os.path.join(OUT, f"chi_tam{tam:g}_t1{t1:g}.csv")
        lines = open(f).read().splitlines()
        hdr = lines[0] + " " + lines[1]
        grid = np.array([[float(v) for v in l.split(",")[2:]]
                         for l in lines if l.startswith("chi,")])
        P = np.fft.fftshift(grid)                     # center q=0
        im = a.imshow(P, origin="lower", cmap="magma",
                      extent=[-np.pi, np.pi, -np.pi, np.pi])
        a.set_title(rf"$t_A={tam}$, $t'={t1}$", fontsize=11)
        a.set_xlabel(r"$q_x$")
        plt.colorbar(im, ax=a, shrink=0.8)
        kmax = np.unravel_index(grid.argmax(), grid.shape)
        qx = 2 * np.pi * kmax[0] / L; qy = 2 * np.pi * kmax[1] / L
        a.text(0.03, 0.94, rf"$Q^*=({qx/np.pi:.1f}\pi,{qy/np.pi:.1f}\pi)$",
               transform=a.transAxes, color="w", fontsize=10)
    ax[0].set_ylabel(r"$q_y$")
    fig.suptitle(r"DQMC $\chi_{zz}(\mathbf{q})$, $4{\times}4$, $U=4$, $\beta=3$, half filling "
                 r"-- combined $t_A$ (NN) + $t'$ (NNN) spin-dependent hopping", fontsize=12)
    fig.tight_layout()
    out = os.path.join(OUT, "chi_spin_pilot.png")
    fig.savefig(out, dpi=160)
    print("wrote", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--point", type=int, default=-1, help="index into POINTS; -1 = all")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--beta", type=float, default=3.0)
    ap.add_argument("--nmeas", type=int, default=600)
    a = ap.parse_args()
    if a.plot:
        plot()
    elif a.point >= 0:
        run_point(*POINTS[a.point], beta=a.beta, nmeas=a.nmeas)
    else:
        for p in POINTS:
            run_point(*p, beta=a.beta, nmeas=a.nmeas)
