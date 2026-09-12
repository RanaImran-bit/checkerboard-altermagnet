"""The altermagnetic order parameter: the d_xy anisotropy of the spin correlations.

The AFM structure factor is NOT the altermagnetic quantity. An altermagnet is
spin-compensated, so it has no net polarisation, and S_AFM only measures how
strongly spins are correlated, not the symmetry of that correlation.

What distinguishes an altermagnet is that the spin response carries a d-wave
form factor fixed by the lattice. On the checkerboard the C4-odd harmonic is
sin(qx) sin(qy), so the order parameter is

    Psi_dxy = (1/N) sum_q  sin(qx) sin(qy)  S_zz(q)

already computed by checkerboard_mag_scan.py as psi_harmonic(Sq, "dxy").
Psi_dx2y2, built from cos(qx) - cos(qy), is the wrong symmetry here and acts as
the control: it should stay zero.

Symmetry requires Psi_dxy = 0 at delta = 0, where the two sublattices are still
related by a translation and the state is an ordinary antiferromagnet.
"""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
    "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "figure.dpi": 130})

D = "/Users/liujiaxin/Desktop/checkerboard-altermagnet/data"
d = pd.read_csv(f"{D}/magnetic_master.csv"); d = d[(d.L == 12) & (d.U > 0)]
# All black. With colour gone the four U values are told apart by marker shape,
# with filled/open alternating and a distinct line style, which is how Huang's
# own figures do it and stays readable in greyscale print.
MK = {2.: "o", 4.: "s", 6.: "^", 8.: "D"}
FC = {2.: "k", 4.: "white", 6.: "k", 8.: "white"}
LS = {2.: "-", 4.: "--", 6.: "-.", 8.: ":"}
US = [2., 4., 6., 8.]

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))

# (a) Psi_dxy against delta, half filling
h = d[np.isclose(d.n, 1.0)]
g = h.groupby(["U", "delta"])[["Psi_dxy", "Psi_dx2y2"]].agg(["mean", "sem"])
a = ax[0]
for U in US:
    s = g.loc[U]
    a.errorbar(s.index, np.abs(s[("Psi_dxy", "mean")]) * 1e3,
               yerr=s[("Psi_dxy", "sem")] * 1e3, color="k", marker=MK[U], mfc=FC[U], mec="k",
               ls=LS[U], ms=5, lw=1.3, capsize=2.5,
               label=f"$U={U:g}$")
a.axhline(0, color="k", lw=1.0, ls="--")
a.set_xlabel(r"anisotropy  $\delta$")
a.set_ylabel(r"$|\Psi_{d_{xy}}|\times 10^{3}$")
a.set_title(r"(a) altermagnetic order parameter", fontsize=15)
a.legend(fontsize=10, loc="best")

# (b) the wrong-symmetry control
a = ax[1]
for U in US:
    s = g.loc[U]
    a.errorbar(s.index, np.abs(s[("Psi_dx2y2", "mean")]) * 1e3,
               yerr=s[("Psi_dx2y2", "sem")] * 1e3, color="k", marker=MK[U], mfc=FC[U], mec="k",
               ls=LS[U], ms=5, lw=1.3, capsize=2.5,
               label=f"$U={U:g}$")
a.axhline(0, color="k", lw=1.0, ls="--")
a.set_xlabel(r"anisotropy  $\delta$")
a.set_ylabel(r"$|\Psi_{d_{x^2-y^2}}|\times 10^{3}$")
a.set_title(r"(b) control: wrong symmetry", fontsize=15)
# Shared scale so the control reads against panel (a), but wide enough that no
# error bar is clipped. Clipped bars ran off the top and looked like stray
# vertical rules rather than uncertainties.
lo = min(ax[0].get_ylim()[0], ax[1].get_ylim()[0])
hi = max(ax[0].get_ylim()[1], ax[1].get_ylim()[1])
for _a in (ax[0], ax[1]):
    _a.set_ylim(lo, hi)
a.legend(fontsize=10, loc="best")

# (c) against filling at fixed delta
a = ax[2]
for U in US:
    rows = []
    for nv in sorted(d.n.unique()):
        s = d[(d.U == U) & np.isclose(d.n, nv) & (d.delta == 0.7)]
        rows.append((nv, np.abs(s.Psi_dxy.mean()) * 1e3, s.Psi_dxy.sem() * 1e3))
    r = np.array(rows)
    a.errorbar(r[:, 0], r[:, 1], yerr=r[:, 2], color="k", marker=MK[U], mfc=FC[U], mec="k",
               ls=LS[U], ms=5, lw=1.3, capsize=2.5,
               label=f"$U={U:g}$")
a.axhline(0, color="k", lw=1.0, ls="--")
a.set_xlabel(r"filling  $n$")
a.set_ylabel(r"$|\Psi_{d_{xy}}|\times 10^{3}$")
a.set_title(r"(c) filling dependence at $\delta=0.7$", fontsize=15)
a.legend(fontsize=10, loc="best")
fig.tight_layout()
fig.savefig(f"{D}/../figs/fig_psi_am.png", dpi=600, bbox_inches="tight",
            facecolor="white")

print("Psi_dxy at half filling: significance in sigma from zero")
print(f"{'delta':>6}" + "".join(f"{'U='+str(int(u)):>10}" for u in US))
for dl in sorted(h.delta.unique()):
    row = f"{dl:6.1f}"
    for U in US:
        m = g.loc[(U, dl), ("Psi_dxy", "mean")]
        e = g.loc[(U, dl), ("Psi_dxy", "sem")]
        row += f"{abs(m)/e if e > 0 else 0:10.1f}"
    print(row)
print("\ncontrol Psi_dx2y2, same units:")
print(f"{'delta':>6}" + "".join(f"{'U='+str(int(u)):>10}" for u in US))
for dl in sorted(h.delta.unique()):
    row = f"{dl:6.1f}"
    for U in US:
        m = g.loc[(U, dl), ("Psi_dx2y2", "mean")]
        e = g.loc[(U, dl), ("Psi_dx2y2", "sem")]
        row += f"{abs(m)/e if e > 0 else 0:10.1f}"
    print(row)
print("\ngrowth with U at delta=0.7 (half filling), |Psi_dxy| x10^3:")
for U in US:
    print(f"  U={U:g}: {abs(g.loc[(U,0.7),('Psi_dxy','mean')])*1e3:.2f}")
