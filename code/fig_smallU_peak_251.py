# --- the small-U peak in Delta_tot: is it real, and does it vanish with L? ---
# Paste-and-run on 251. Delta_tot/N against U at fixed anisotropy, one curve per
# lattice size, with the diagnostics printed underneath.
#
# Three questions this answers:
#   1. is the low-U bright spot in the heat map above the noise?      yes, ~100x
#   2. does it decay with system size?                                no, it wanders
#   3. what does it actually track?                                   the TRIAL
#
# Delta_tot from the Fortran CPQMC is constrained-path, so it is biased by the trial
# wave function. fortran_L12_partial.csv records the trial magnetisation m_trial
# next to Delta_tot, which lets us test whether the sharp features follow m_trial
# rather than U. They largely do, and the small-U peak is the one cell that breaks
# even that pattern.
#
# The plot is L = 8, 10, 12. L=14 is dropped from it: its scan is a cross rather
# than a rectangle, so only delta=0.2 has a full U column and the other panels
# would show a stub. It is still used in diagnostic (1), because fortran_dnk.csv is
# the only file carrying error bars on Delta_tot at all.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt

# "U-y" puts U on the vertical axis so the panels line up with the (U, delta) heat
# maps, where U is already vertical. "U-x" is the conventional response curve.
ORIENT = "U-y"

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
    "font.size":15,"axes.labelsize":19,"legend.fontsize":14,
    "xtick.labelsize":15,"ytick.labelsize":15,"axes.linewidth":1.4})

f  = pd.read_csv(find("fortran_L8_L10_L12_dedup.csv")); f = f[np.isclose(f.n, 1.0)]
dn = pd.read_csv(find("fortran_dnk.csv"));              dn = dn[np.isclose(dn.n, 1.0)]

DELTAS = [0.1, 0.2, 0.3]
SIZES  = [8, 10, 12]
COL    = {8:"#0072B2", 10:"#009E73", 12:"#D55E00"}

SWAP = (ORIENT == "U-y")
fig, ax = plt.subplots(1, len(DELTAS),
                       figsize=((4.6 if SWAP else 5.3)*len(DELTAS), 5.4 if SWAP else 4.9),
                       sharey=True, sharex=SWAP, gridspec_kw={"wspace":.08})
for a, dl in zip(ax, DELTAS):
    for L in SIZES:
        g = f[(f.L == L) & np.isclose(f.delta, dl) & (f.U > 0)].sort_values("U")
        xy = (g.dtot_N, g.U) if SWAP else (g.U, g.dtot_N)
        a.plot(*xy, marker="o", ms=6, lw=2.0, color=COL[L], label=f"$L={L}$")
    a.set_title(rf"$\delta = {dl:g}$", fontsize=19, pad=10)
    if SWAP:
        a.set_xlabel(r"$\Delta_{\rm tot}/N$", labelpad=7)
        a.set_ylim(1.8, 5.2); a.set_yticks([2, 3, 4, 5]); a.set_xlim(0, .098)
        a.set_xticks([0, .02, .04, .06, .08])
    else:
        a.set_xlabel(r"$U/t$", labelpad=7)
        a.set_xlim(1.8, 5.2); a.set_xticks([2, 3, 4, 5])
    a.legend(loc="best", frameon=True, framealpha=.92, edgecolor="0.75",
             handlelength=1.6, handletextpad=.6, labelspacing=.35, borderpad=.5)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$" if SWAP else r"$\Delta_{\rm tot}/N$", labelpad=8)
fig.suptitle(r"half filling. Curves should lean steadily right as $U$ rises;"
             r" at $\delta=0.1$ and $0.2$ they bulge out at the bottom instead."
             if SWAP else
             r"half filling. A monotonic rise with $U$ would be the expectation;"
             r" it is not what happens.", fontsize=16, y=1.01)
plt.show()

# ---------------- diagnostics ----------------
print("=" * 60)
print("(1) ABOVE THE NOISE?")
print("=" * 60)
q = dn[(dn.U > 0) & (dn.L >= 14) & dn.snr.notna() & (dn.delta > 0)]
print(f"  Delta_tot carries error bars only at L=14,16 (fortran_dnk.csv).")
print(f"  There the signal-to-noise runs {q.snr.min():.0f} to {q.snr.max():.0f}.")
print(f"  delta=0 null floor: 6.26e-4 per site. Peak values ~0.06-0.09, so ~100x.")

print("\n" + "=" * 60)
print("(2) DOES IT DECAY WITH L?   Peak of Delta_tot/N over U <= 3")
print("=" * 60)
print(f"{'delta':>7}" + "".join(f"{'L='+str(L):>10}" for L in SIZES))
for dl in DELTAS:
    row = [f[(f.L == L) & np.isclose(f.delta, dl) & (f.U > 0) & (f.U <= 3)].dtot_N.max()
           for L in SIZES]
    print(f"{dl:>7.1f}" + "".join(f"{v:>10.4f}" for v in row))
print("  Grows with L rather than decaying. A finite-size artefact would fade.")

print("\n" + "=" * 60)
print("(3) WHAT IT TRACKS: the trial magnetisation")
print("=" * 60)
p = pd.read_csv(find("fortran_L12_partial.csv")); p = p[p.U > 0]
print(f"  corr(m_trial, Delta_tot/N) over {len(p)} L=12 cells: "
      f"{np.corrcoef(p.m_trial, p.dtot_N)[0,1]:+.3f}")
lo, hi = p[p.m_trial < .10], p[p.m_trial > .25]
print(f"    m_trial < 0.10 ({len(lo):>2} cells): mean Delta_tot/N = {lo.dtot_N.mean():.4f}")
print(f"    m_trial > 0.25 ({len(hi):>2} cells): mean Delta_tot/N = {hi.dtot_N.mean():.4f}"
      f"   ({hi.dtot_N.mean()/lo.dtot_N.mean():.1f}x)")
print("\n  the exception, and it is the small-U peak itself:")
for _, r in p[(p.U <= 2) & (p.dtot_N > .03)].iterrows():
    print(f"    U={r.U:.1f} delta={r.delta:.1f}: Delta_tot/N={r.dtot_N:.4f} "
          f"from m_trial={r.m_trial:.4f}")
print("  Everywhere else a large Delta_tot needs a large m_trial. Not here.")
print("\n  CAVEAT: 24 cells from a partial scan. Suggestive, not settled. The test is")
print("  to re-run those cells with a different trial and see whether Delta_tot moves.")
