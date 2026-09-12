#!/usr/bin/env python -u
"""Figure 6: the polarization against symmetry-breaking field, and its
extrapolation to vanishing field against system size.

Sections V of the v1 and v3 manuscripts argue the absence of long-range order
entirely from this extrapolation, and no figure existed for it.

Panel (a): Delta_tot/N against h at each L with the h -> 0 fit.

The fitted quantity is Delta_tot/N, following code/analyze_hscan.py, NOT m_stag.
The two differ by about a factor of two and only Delta_tot/N reproduces the
intercepts quoted in the text.
Panel (b): the intercepts against 1/L, with the two candidate fits made on the
L = 8 to 14 points only and extended to L = 16, so the reader can see that the
measured L = 16 value falls far below both.

Data: data/hscan/polarization_L{8,10,12,14,16}.csv, half filling, U = 4,
delta = 0.3, periodic boundaries.
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D    = os.path.join(REPO, "data", "hscan")
OUT  = os.path.join(REPO, "manuscript", "figures")

plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif",
    "font.size":9,"axes.labelsize":10,"legend.fontsize":8,
    "xtick.labelsize":9,"ytick.labelsize":9,"axes.linewidth":0.9})
# Okabe-Ito, colour-vision safe
COL = {8:"#0072B2", 10:"#009E73", 12:"#D55E00", 14:"#CC79A7", 16:"#000000"}

SIZES = [8, 10, 12, 14, 16]
fig, ax = plt.subplots(1, 2, figsize=(6.75, 2.8))

inter = {}
for L in SIZES:
    f = os.path.join(D, f"polarization_L{L}.csv")
    if not os.path.exists(f): continue
    d = pd.read_csv(f)
    d["dtot_N"] = d.delta_tot / (d.L ** 2)
    g = d.groupby("h")["dtot_N"].agg(["mean", "std", "count"]).reset_index()
    g["err"] = g["std"] / np.sqrt(g["count"].clip(lower=1))
    ax[0].errorbar(g.h, g["mean"], yerr=g.err.fillna(0), marker="o", ms=3.2, lw=1.2,
                   capsize=2, color=COL[L], label=f"$L={L}$")
    # linear fit in h over the smallest fields, which is what the intercept means
    g = g.sort_values("h")
    sel = g.iloc[:3]                       # the three smallest fields, as in analyze_hscan.py
    c = np.polyfit(sel.h, sel["mean"], 1)
    hh = np.linspace(0, sel.h.max(), 20)
    ax[0].plot(hh, np.polyval(c, hh), ls=":", lw=0.9, color=COL[L])
    # intercept uncertainty by reseeding the fit inside the error bars
    rng = np.random.default_rng(0)
    e = sel.err.fillna(0).values
    boot = [np.polyfit(sel.h.values, rng.normal(sel["mean"].values, e), 1)[-1]
            for _ in range(2000)]
    inter[L] = (float(np.polyval(c, 0.0)), float(np.std(boot)))

ax[0].set_xlabel(r"$h$"); ax[0].set_ylabel(r"$\Delta_{\rm tot}/N$")
ax[0].set_xlim(left=0); ax[0].legend(frameon=False, handlelength=1.3, labelspacing=.25)
ax[0].text(.04, .93, "(a)", transform=ax[0].transAxes, fontweight="bold")

Ls  = np.array(sorted(inter))
y   = np.array([inter[L][0] for L in Ls])
ye  = np.array([inter[L][1] for L in Ls])
fitL = Ls[Ls <= 14]; fity = y[Ls <= 14]
c1 = np.polyfit(1/fitL, fity, 1)            # linear in 1/L
c2 = np.polyfit(1/fitL**2, fity, 1)         # linear in 1/L^2
xs = np.linspace(0, 1/7.5, 100)
ax[1].plot(xs, np.polyval(c1, xs),    lw=1.2, color="#0072B2", label=r"fit in $1/L$, $L\leq14$")
ax[1].plot(xs, np.polyval(c2, xs**2), lw=1.2, ls="--", color="#D55E00", label=r"fit in $1/L^2$, $L\leq14$")
ax[1].errorbar(1/fitL, fity, yerr=ye[Ls <= 14], ls="none", marker="o", ms=4.5,
               capsize=2.5, color="k", label=r"$L=8$ to $14$")
if 16 in inter:
    ax[1].errorbar([1/16], [inter[16][0]], yerr=[inter[16][1]], ls="none", marker="s",
                   ms=6, capsize=2.5, mfc="w", mec="k", mew=1.4, color="k", label=r"$L=16$")
ax[1].axhline(0, lw=.7, color="0.6")
ax[1].set_xlabel(r"$1/L$"); ax[1].set_ylabel(r"$\Delta_{\rm tot}/N\,(h\to0)$")
ax[1].set_xlim(0, 1/7.5); ax[1].legend(frameon=False, handlelength=1.6, labelspacing=.25)
ax[1].text(.04, .93, "(b)", transform=ax[1].transAxes, fontweight="bold")

fig.tight_layout(pad=0.4)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"fig06_fss_hscan.{ext}"), dpi=300, bbox_inches="tight")

print(f"{'L':>4}{'dtot_N(h->0)':>14}{'err':>10}")
for L in Ls: print(f"{L:>4}{inter[L][0]:>14.5f}{inter[L][1]:>10.5f}")
if 16 in inter:
    print(f"\n  fits from L<=14 predict at L=16: {np.polyval(c1,1/16):.5f} (1/L), "
          f"{np.polyval(c2,(1/16)**2):.5f} (1/L^2)")
    print(f"  measured at L=16: {inter[16][0]:.5f}")
    print(f"  ratio to L=14: {inter[16][0]/inter[14][0]:.3f}   (1/L would give {14/16:.3f})")
print("\nwrote manuscript/figures/fig06_fss_hscan.pdf and .png")
