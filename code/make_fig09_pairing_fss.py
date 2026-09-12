"""Figure 9: pairing vertex at half filling, and its finite-size behaviour.

The vertex is the interaction-induced part of the pair-field susceptibility, the
full correlator minus its uncorrelated Wick piece. It is the meaningful measure of
whether U builds pairing, since the full correlator is large and positive already
for free electrons. All three sizes use byte-identical drivers and settings
(beta = 32, N_w = 500, tau window 0.8, six seeds).

(a) All four channels at L = 12. The on-site s vertex is negative throughout, as a
    repulsive U requires, and d_xy leads for delta >= 0.3.
(b) d_xy at L = 8, 10, 12. It GROWS with system size at every delta where it leads,
    which is the opposite of the magnetic order parameter's behaviour and is what
    separates a pairing tendency from a finite-size artifact.
(c) The same growth at the d_xy peak, delta = 0.4, against 1/N.

  python make_fig09_pairing_fss.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "font.size": 8, "axes.labelsize": 9, "legend.fontsize": 7,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "lines.linewidth": 1.3, "lines.markersize": 4.2,
    "figure.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data")
SRC = {8: os.path.join(D, "chi_fss_L8_U4.csv"),
       10: os.path.join(D, "chi_fss_L10_U4.csv"),
       12: os.path.join(D, "chi_L12", "eqtime_L12_U4.csv")}
dat = {L: pd.read_csv(f) for L, f in SRC.items() if os.path.exists(f)}
Ls = sorted(dat)

CH = [("son", r"on-site $s$", "#E69F00", "D"),
      ("sext", r"extended $s$", "#009E73", "^"),
      ("d", r"$d_{x^2-y^2}$", "#0072B2", "o"),
      ("dxy", r"$d_{xy}$", "#D55E00", "s")]
CL = {8: "#0072B2", 10: "#D55E00", 12: "#009E73"}
MK = {8: "o", 10: "s", 12: "^"}

def curve(df, ch):
    g = df.groupby("delta")[f"chi_{ch}_vertex"]
    x = np.array(sorted(df.delta.unique()))
    return x, g.mean().reindex(x).to_numpy(), (g.std(ddof=1)/np.sqrt(g.count())).reindex(x).to_numpy()

fig, axs = plt.subplots(1, 3, figsize=(6.75, 2.35),
                        gridspec_kw={"wspace": 0.36})

a = axs[0]
for ch, lab, c, m in CH:
    x, y, e = curve(dat[12], ch)
    a.errorbar(x, y, yerr=e, fmt=m + "-", color=c, mfc="white", mew=1.1,
               capsize=1.6, elinewidth=0.7, label=lab)
a.axhline(0, color="0.6", lw=0.7, ls=":")
a.set_xlabel(r"$\delta$"); a.set_ylabel(r"$\chi^{\mathrm{vertex}}_{\alpha}$")
a.set_title("(a)", loc="left", fontsize=9)
a.legend(frameon=False, handlelength=1.3, labelspacing=0.18, borderaxespad=0.3,
         loc="lower right", ncol=2, columnspacing=0.9)
a.tick_params(direction="in", top=True, right=True)
a.text(0.96, 0.94, r"$L=12$", transform=a.transAxes, fontsize=7, color="0.3",
       va="top", ha="right")

b = axs[1]
for L in Ls:
    x, y, e = curve(dat[L], "dxy")
    b.errorbar(x, y, yerr=e, fmt=MK[L] + "-", color=CL[L], mfc="white", mew=1.1,
               capsize=1.6, elinewidth=0.7, label=rf"$L={L}$")
b.axhline(0, color="0.6", lw=0.7, ls=":")
b.set_xlabel(r"$\delta$"); b.set_ylabel(r"$\chi^{\mathrm{vertex}}_{d_{xy}}$")
b.set_title("(b)", loc="left", fontsize=9)
b.legend(frameon=False, handlelength=1.4, labelspacing=0.18, borderaxespad=0.3,
         loc="lower right")
b.tick_params(direction="in", top=True, right=True)

c = axs[2]
for dl, col, mk in [(0.4, "#D55E00", "s"), (0.5, "#009E73", "^"), (0.7, "#E69F00", "D")]:
    xs, ys, es = [], [], []
    for L in Ls:
        s = dat[L][np.isclose(dat[L].delta, dl)]["chi_dxy_vertex"]
        if len(s):
            xs.append(1.0/(L*L)); ys.append(s.mean()); es.append(s.std(ddof=1)/np.sqrt(len(s)))
    c.errorbar(xs, ys, yerr=es, fmt=mk + "-", color=col, mfc="white", mew=1.1,
               capsize=1.6, elinewidth=0.7, label=rf"$\delta={dl}$")
c.set_xlabel(r"$1/N$"); c.set_ylabel(r"$\chi^{\mathrm{vertex}}_{d_{xy}}$")
c.set_title("(c)", loc="left", fontsize=9)
c.legend(frameon=False, handlelength=1.4, labelspacing=0.18, borderaxespad=0.3)
c.tick_params(direction="in", top=True, right=True)
c.set_xlim(0, 0.0185)
c.margins(y=0.16)
top = dat  # label each size once, above the delta = 0.4 series
for L in Ls:
    yv = dat[L][np.isclose(dat[L].delta, 0.4)]["chi_dxy_vertex"].mean()
    c.annotate(rf"$L={L}$", xy=(1.0/(L*L), yv), xytext=(0, 7),
               textcoords="offset points", ha="center", fontsize=6.4, color="0.4")

out = os.path.join(HERE, "..", "manuscript", "figures", "fig09_pairing_fss")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png")
print(f"wrote {out}.pdf / .png   sizes={Ls}")
for dl in (0.3, 0.4, 0.5, 0.7):
    v = [dat[L][np.isclose(dat[L].delta, dl)]["chi_dxy_vertex"].mean() for L in Ls]
    print(f"  delta={dl}: " + "  ".join(f"L{L}={x:.3f}" for L, x in zip(Ls, v))
          + f"   growth {v[-1]/v[0]:.2f}x")
