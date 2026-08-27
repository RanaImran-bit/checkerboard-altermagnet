# --- the small-U peak in Delta_tot, EVERY lattice size we have. ---
# Paste-and-run on 251.
#
# L = 8, 10, 12   complete 6x7 (U, delta) grids, no error bars recorded
# L = 14          delta = 0.1, 0.2, 0.3 only, and U=2 exists ONLY at delta=0.2
# L = 16          two cells total, U=4 at delta = 0.2 and 0.3, drawn as bare markers
# L = 6           DOES NOT EXIST for Delta_tot at half filling. There are L=6 pairing
#                 and chi_zz files, but no polarisation run. It would need new cells.
#
# The gap that matters: at delta = 0.1 and 0.3 the U=2 point was never run at L=14,
# which is precisely where the peak sits. So the largest-lattice check of the peak
# exists at delta = 0.2 and nowhere else. Two cells would close it.
#
# Delta_tot is constrained-path and therefore biased by the trial wave function.
# See the printout: the peak lives where Delta_tot stops tracking m_trial.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt

ORIENT = "U-y"          # "U-y" matches the heat maps; "U-x" is the response curve

def find(name):
    for c in [name, f"../data/{name}", f"data/{name}",
              os.path.expanduser(f"~/results/figdata/{name}"),
              os.path.expanduser(f"~/results/{name}"),
              os.path.expanduser(f"~/checkerboard-altermagnet/data/{name}")]:
        if os.path.exists(c): return c
    h = glob.glob(os.path.expanduser(f"~/results/**/{name}"), recursive=True)
    if h: return h[0]
    raise SystemExit(f"cannot find {name}. Try:  find ~ -name {name}")

plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif",
    "font.size":15,"axes.labelsize":19,"legend.fontsize":13,
    "xtick.labelsize":15,"ytick.labelsize":15,"axes.linewidth":1.4})

f  = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv")); f = f[np.isclose(f.n, 1.0)]
dn = pd.read_csv(find("fortran_dnk.csv"));              dn = dn[np.isclose(dn.n, 1.0)]

DELTAS = [0.1, 0.2, 0.3]
COL = {6:"#999999", 8:"#0072B2", 10:"#009E73", 12:"#D55E00", 14:"#000000", 16:"#CC79A7"}
SWAP = (ORIENT == "U-y")

fig, ax = plt.subplots(1, len(DELTAS),
                       figsize=((4.8 if SWAP else 5.3)*len(DELTAS), 5.6 if SWAP else 4.9),
                       sharey=True, sharex=SWAP, gridspec_kw={"wspace":.08})
for a, dl in zip(ax, DELTAS):
    for L in (8, 10, 12):
        g = f[(f.L == L) & np.isclose(f.delta, dl) & (f.U > 0)].sort_values("U")
        xy = (g.dtot_N, g.U) if SWAP else (g.U, g.dtot_N)
        a.plot(*xy, marker="o", ms=6, lw=2.0, color=COL[L], label=f"$L={L}$")
    for L, mk in ((14, "s"), (16, "D")):
        g = dn[(dn.L == L) & np.isclose(dn.delta, dl) & (dn.U >= 2)].sort_values("U")
        if not len(g): continue
        e  = g.delta_tot_err.values / (L*L)          # per site, matching dtot_N
        xy = (g.dtot_N, g.U) if SWAP else (g.U, g.dtot_N)
        err = {"xerr": e} if SWAP else {"yerr": e}
        a.errorbar(*xy, marker=mk, ms=7, capsize=3, color=COL[L],
                   lw=2.0 if len(g) > 1 else 0, label=f"$L={L}$", **err)
    a.set_title(rf"$\delta = {dl:g}$", fontsize=19, pad=10)
    if SWAP:
        a.set_xlabel(r"$\Delta_{\rm tot}/N$", labelpad=7)
        a.set_ylim(1.8, 5.2); a.set_yticks([2, 3, 4, 5])
        a.set_xlim(0, .098); a.set_xticks([0, .02, .04, .06, .08])
    else:
        a.set_xlabel(r"$U/t$", labelpad=7)
        a.set_xlim(1.8, 5.2); a.set_xticks([2, 3, 4, 5])
    a.legend(loc="best", frameon=True, framealpha=.92, edgecolor="0.75",
             handlelength=1.6, handletextpad=.6, labelspacing=.32, borderpad=.5)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$" if SWAP else r"$\Delta_{\rm tot}/N$", labelpad=8)
fig.suptitle(r"half filling, every available lattice size."
             r"  $L=6$ has no polarisation data.", fontsize=16, y=1.01)
plt.show()

print("=" * 62)
print("WHAT EXISTS, per (L, delta), at U > 0")
print("=" * 62)
print(f"{'delta':>7}" + "".join(f"{'L='+str(L):>22}" for L in (8, 10, 12, 14, 16)))
for dl in DELTAS:
    cells = []
    for L in (8, 10, 12, 14, 16):
        src = f if L <= 12 else dn
        g = src[(src.L == L) & np.isclose(src.delta, dl) & (src.U > 0)]
        us = sorted(float(u) for u in g.U.unique())
        cells.append(("-" if not us else
                      ("U=" + ",".join(f"{u:g}" for u in us))[:21]))
    print(f"{dl:>7.1f}" + "".join(f"{c:>22}" for c in cells))

print("\n  L=6: no Delta_tot at half filling anywhere, local or on 250-256.")
print("       (L=6 pairing and chi_zz files exist, but not the polarisation.)")
print("  L=14 is MISSING U=2 at delta=0.1 and 0.3 -- exactly where the peak is.")
print("       Two cells would let the largest lattice test the peak at all three delta.")

print("\n" + "=" * 62)
print("PEAK OF Delta_tot/N OVER U <= 3")
print("=" * 62)
print(f"{'delta':>7}" + "".join(f"{'L='+str(L):>10}" for L in (8, 10, 12, 14)))
for dl in DELTAS:
    row = []
    for L in (8, 10, 12, 14):
        src = f if L <= 12 else dn
        g = src[(src.L == L) & np.isclose(src.delta, dl) & (src.U > 0) & (src.U <= 3)]
        row.append(g.dtot_N.max() if len(g) else np.nan)
    print(f"{dl:>7.1f}" + "".join(f"{v:>10.4f}" for v in row))
print("  delta=0.2 is the only column where L=14 tests the U=2 peak: 0.0738")
print("  against 0.0804 at L=12. Slightly lower, not gone.")
