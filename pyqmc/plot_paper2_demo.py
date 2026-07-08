#!/usr/bin/env python3
"""Paper-2 demo figure: chi_zz(q) and chi_pm(q) maps for the combined altermagnet
model (tA + t') at the six phase-representative (n, tA, t') points, plus a
per-point method comparison (DQMC / CP-DQMC / CPQMC[x2 T=0]).

Reads results/paper2_demo/{method}_{tag}.json (from paper2_chi_point.py runs).

    python pyqmc/plot_paper2_demo.py       # -> docs/paper2_chi_demo.png
"""
from __future__ import annotations
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.join(os.path.dirname(__file__), "..")
DEMO = os.path.join(ROOT, "results", "paper2_demo")
OUT = os.path.join(ROOT, "docs", "paper2_chi_demo.png")

POINTS = [("1.0", "0.3", "0.0"), ("1.0", "0.5", "0.5"), ("0.9", "0.1", "0.7"),
          ("0.9", "0.5", "0.1"), ("0.75", "0.3", "0.3"), ("1.0", "0.9", "0.1")]
METHODS = ["dqmc", "cpdqmc", "cpqmc"]
MLABEL = {"dqmc": "DQMC", "cpdqmc": "CP-DQMC", "cpqmc": r"CPQMC ($2\chi^{T=0}$)"}


def load(m, n, ta, tt):
    p = os.path.join(DEMO, f"{m}_n{n}_ta{ta}_tt{tt}.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def qmap(ax, grid, title, cmap="magma"):
    g = np.fft.fftshift(np.array(grid))
    im = ax.imshow(g.T, origin="lower", cmap=cmap,
                   extent=[-1, 1, -1, 1], aspect="equal")
    ax.set_title(title, fontsize=8)
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
    ax.tick_params(labelsize=6)
    return im


def main():
    ncol = len(POINTS)
    fig, axes = plt.subplots(3, ncol, figsize=(2.3 * ncol, 7.6))
    for c, (n, ta, tt) in enumerate(POINTS):
        recs = {m: load(m, n, ta, tt) for m in METHODS}
        base = recs.get("cpdqmc") or recs.get("dqmc") or recs.get("cpqmc")
        hdr = f"n={n}  $t_A$={ta}  $t'$={tt}"
        if base is None:
            axes[0][c].set_title(hdr + "\n(missing)", fontsize=8)
            continue
        # rows 1-2: CP-DQMC maps (workhorse); chi_zz and chi_pm on the same footing
        im0 = qmap(axes[0][c], base["chi_q"], hdr + "\n" + r"$\chi_{zz}(q)$")
        im1 = qmap(axes[1][c], base["chi_pm_q"], r"$\chi_{+-}(q)$", cmap="viridis")
        fig.colorbar(im0, ax=axes[0][c], fraction=0.046, pad=0.02).ax.tick_params(labelsize=5)
        fig.colorbar(im1, ax=axes[1][c], fraction=0.046, pad=0.02).ax.tick_params(labelsize=5)
        # row 3: method comparison of the peak values (zz solid, pm hatched)
        ax = axes[2][c]
        w = 0.35
        for k, m in enumerate(METHODS):
            r = recs.get(m)
            if r is None:
                continue
            pz, pp = r["chi_q_peak"], r["chi_pm_q_peak"]
            ax.bar(k - w / 2, pz["value"], w, color="C0")
            ax.bar(k + w / 2, pp["value"], w, color="C2", hatch="//")
            ax.text(k, 0.02, f"({pz['qx_over_pi']:+.2f},{pz['qy_over_pi']:+.2f})$\\pi$",
                    rotation=90, ha="center", va="bottom", fontsize=5.5)
            sgn = r.get("sign", 1.0)
            lab = MLABEL[m] + (f"\ns={sgn:.2f}" if sgn < 0.995 else "")
            ax.text(k, -0.02, lab, ha="center", va="top", fontsize=5.5,
                    transform=ax.get_xaxis_transform())
        ax.set_xticks([])
        if c == 0:
            ax.set_ylabel("peak value\n(zz blue / pm green)", fontsize=7)
        ax.tick_params(labelsize=6)
        for r0 in (0, 1):
            if c == 0:
                axes[r0][c].set_ylabel(r"$q_y/\pi$", fontsize=7)
            axes[r0][c].set_xlabel(r"$q_x/\pi$", fontsize=6)
    fig.suptitle("Combined altermagnet model ($t_A$ + $t'$, U=4, L=8, $\\beta$=5): "
                 r"$\chi_{zz}(q)$ vs $\chi_{+-}(q)$ [CP-DQMC maps; bars: all methods]",
                 fontsize=10, y=0.995)
    fig.tight_layout()
    fig.savefig(OUT, dpi=170)
    print(f"wrote {OUT}")
    # peak table to stdout
    print(f"\n{'point':<22} {'method':<8} {'dens':>6} {'sign':>6} "
          f"{'zz peak (q/pi)':>18} {'pm peak (q/pi)':>18} {'pm/2zz':>7}")
    for (n, ta, tt) in POINTS:
        for m in METHODS:
            r = load(m, n, ta, tt)
            if r is None:
                continue
            pz, pp = r["chi_q_peak"], r["chi_pm_q_peak"]
            print(f"n={n} tA={ta} t'={tt:<6} {m:<8} {r.get('dens', float('nan')):>6.3f} "
                  f"{r.get('sign', 1.0):>6.3f} "
                  f"{pz['value']:>7.3f} ({pz['qx_over_pi']:+.2f},{pz['qy_over_pi']:+.2f}) "
                  f"{pp['value']:>7.3f} ({pp['qx_over_pi']:+.2f},{pp['qy_over_pi']:+.2f}) "
                  f"{pp['value'] / (2 * pz['value']):>7.2f}")


if __name__ == "__main__":
    main()
