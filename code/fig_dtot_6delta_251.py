# --- Delta_tot against U at six anisotropies, L = 8 to 14. ---
# Paste-and-run on 251.
#
# Six panels, delta = 0.1 to 0.6, laid out 2 x 3.
#
# L = 8, 10, 12   complete 6x7 (U, delta) grids, no error bars recorded
# L = 14          delta = 0.1, 0.2, 0.3 ONLY, so it appears in the first three
#                 panels and not the last three. Complete over U = 2 to 5 after
#                 the two U=2 cells finished on 28 Aug. Carries error bars.
# L = 16          exists (U=4 only, two cells) but is NOT plotted: a single isolated
#                 marker at U=4 says nothing about a peak that lives at U=2.
# L = 6           does not exist for Delta_tot at half filling anywhere.
#
# What it shows. At delta=0.1 the peak is GONE at L=14: Delta_tot/N falls from
# 0.0603 at L=12 to 0.0012, about twice the delta=0 null floor. That run has
# signal-to-noise 1.0 and correlation 0.007 with sin kx sin ky, so it holds no
# polarisation and no d_xy character. At delta=0.3 the decay is steady, roughly
# 1/L. At delta=0.2 the peak persists at every size and is still unexplained.
#
# Delta_tot is constrained-path and therefore biased by the trial wave function,
# which is the leading explanation for the delta=0.1 collapse: that lattice
# converged to a nearly unpolarised solution where the smaller ones did not.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt
import matplotlib.patheffects as pe

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

DELTAS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
NCOL   = 3
COL = {8:"#0072B2", 10:"#009E73", 12:"#D55E00", 14:"#000000"}
SWAP = (ORIENT == "U-y")

# Background: Delta_tot/N is the x axis here, so the gradient runs horizontally and
# every point on a curve sits on the colour that its own Delta_tot has in the heat
# maps. VMAX is the same 0.0916 those figures use, so the two are directly readable
# against each other. The curves keep their own colours and get a white casing,
# because jet spans blue through red and both the L=8 blue and the L=12 orange would
# otherwise vanish where the background matches them.
VMAX = 0.0916
XHI  = 0.098
CASE = [pe.withStroke(linewidth=3.4, foreground="w")]

NROW = int(np.ceil(len(DELTAS)/NCOL))
fig, axg = plt.subplots(NROW, NCOL,
                        figsize=(4.7*NCOL, 5.0*NROW),
                        sharey=True, sharex=SWAP,
                        gridspec_kw={"wspace":.08, "hspace":.20})
ax = axg.ravel()
for a in ax[len(DELTAS):]: a.set_visible(False)
GRAD = np.linspace(0, XHI, 512)[None, :]
for a, dl in zip(ax, DELTAS):
    if SWAP:
        a.imshow(GRAD, extent=[0, XHI, 1.8, 5.2], cmap="jet", vmin=0, vmax=VMAX,
                 aspect="auto", origin="lower", zorder=0, interpolation="bilinear")
    for L in (8, 10, 12):
        g = f[(f.L == L) & np.isclose(f.delta, dl) & (f.U > 0)].sort_values("U")
        xy = (g.dtot_N, g.U) if SWAP else (g.U, g.dtot_N)
        # line and markers drawn separately: a path effect applied to one artist
        # strokes its markers too, which washes them out against the gradient
        ln, = a.plot(*xy, lw=2.2, color=COL[L], zorder=4)
        ln.set_path_effects(CASE)
        a.plot(*xy, ls="none", marker="o", ms=6.5, mfc=COL[L], mec="w", mew=1.2,
               zorder=5, label=f"$L={L}$")
    for L, mk in ((14, "s"),):
        g = dn[(dn.L == L) & np.isclose(dn.delta, dl) & (dn.U >= 2)].sort_values("U")
        # a lone point traces no trend; delta=0.4 has one L=14 cell, at U=4 only
        if len(g) < 2: continue
        e  = g.delta_tot_err.values / (L*L)          # per site, matching dtot_N
        xy = (g.dtot_N, g.U) if SWAP else (g.U, g.dtot_N)
        err = {"xerr": e} if SWAP else {"yerr": e}
        eb = a.errorbar(*xy, capsize=3, color=COL[L], zorder=4,
                        lw=2.2 if len(g) > 1 else 0, **err)
        eb.lines[0].set_path_effects(CASE)
        a.plot(*xy, ls="none", marker=mk, ms=7, mfc=COL[L], mec="w", mew=1.2,
               zorder=5, label=f"$L={L}$")
    a.set_title(rf"$\delta = {dl:g}$", fontsize=19, pad=10)
    if SWAP:
        a.set_xlabel(r"$\Delta_{\rm tot}/N$", labelpad=7)
        a.set_ylim(1.8, 5.2); a.set_yticks([2, 3, 4, 5])
        a.set_xlim(0, .098); a.set_xticks([0, .02, .04, .06, .08])
    else:
        a.set_xlabel(r"$U/t$", labelpad=7)
        a.set_xlim(1.8, 5.2); a.set_xticks([2, 3, 4, 5])
    a.legend(loc="upper right", frameon=True, framealpha=.95, edgecolor="0.4",
             handlelength=1.5, handletextpad=.55, labelspacing=.28, borderpad=.45)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
for i, a in enumerate(ax[:len(DELTAS)]):          # label only the outer edges
    if i % NCOL: a.set_ylabel("")
    else:        a.set_ylabel(r"$U/t$" if SWAP else r"$\Delta_{\rm tot}/N$", labelpad=8)
    if i < len(DELTAS) - NCOL: a.set_xlabel("")
fig.suptitle(r"half filling.  $L=14$ forms a curve only for $\delta \leq 0.3$; "
             r"its single $\delta=0.4$ cell is omitted.", fontsize=17, y=0.985)
if SWAP:
    import matplotlib as mpl
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, VMAX), cmap="jet")
    cb = fig.colorbar(sm, ax=list(ax), pad=.013, fraction=.020)
    cb.set_label(r"$\Delta_{\rm tot}/N$", fontsize=19, labelpad=9)
    cb.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.3f"))
    cb.ax.tick_params(labelsize=13)
plt.show()

print("=" * 62)
print("WHAT EXISTS, per (L, delta), at U > 0")
print("=" * 62)
print(f"{'delta':>7}" + "".join(f"{'L='+str(L):>22}" for L in (8, 10, 12, 14)))
for dl in DELTAS:
    cells = []
    for L in (8, 10, 12, 14):
        src = f if L <= 12 else dn
        g = src[(src.L == L) & np.isclose(src.delta, dl) & (src.U > 0)]
        us = sorted(float(u) for u in g.U.unique())
        cells.append(("-" if not us else
                      ("U=" + ",".join(f"{u:g}" for u in us))[:21]))
    print(f"{dl:>7.1f}" + "".join(f"{c:>22}" for c in cells))

print("\n  L=6: no Delta_tot at half filling anywhere, local or on 250-256.")
print("       (L=6 pairing and chi_zz files exist, but not the polarisation.)")
print("  L=14 now covers U=2 at all three delta (the two cells finished 28 Aug).")

print("\n" + "=" * 62)
print("Delta_tot/N AT U = 2, where the peak sits")
print("=" * 62)
print(f"{'delta':>7}" + "".join(f"{'L='+str(L):>10}" for L in (8, 10, 12, 14)) + "   L14/L12")
for dl in DELTAS:
    row = []
    for L in (8, 10, 12, 14):
        src = f if L <= 12 else dn
        g = src[(src.L == L) & np.isclose(src.delta, dl) & np.isclose(src.U, 2.0)]
        row.append(g.dtot_N.iloc[0] if len(g) else np.nan)
    r = row[3] / row[2] if (row[2] and not np.isnan(row[3])) else np.nan
    cells = "".join("       -  " if np.isnan(v) else f"{v:>10.4f}" for v in row)
    print(f"{dl:>7.1f}" + cells + ("       -  " if np.isnan(r) else f"{r:>10.3f}"))
print("\n  (L=14 blank above delta=0.3: not run there.)")
print("\n  delta=0.1: the peak is GONE at L=14. Delta_tot/N = 0.0012 is about twice")
print("             the delta=0 null floor of 6.26e-4, i.e. at the noise. That run")
print("             has signal-to-noise 1.0 and correlation 0.007 with sin kx sin ky,")
print("             so there is no polarisation and no d_xy character in it at all.")
print("  delta=0.3: steady decay, roughly like 1/L.")
print("  delta=0.2: persists at every size. Still unexplained, and it is also the")
print("             cell that broke the pattern in the Neel finite-size scaling.")
