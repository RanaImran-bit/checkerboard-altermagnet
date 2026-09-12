# --- which d channel wins: time-integrated vertex, L = 8, 10, 12. Paste-and-run on 251. ---
# Self-contained: reads the RAW per-U run files and builds the summary itself, so it
# needs no prepared CSV. Sources searched, in order:
#     ~/collected/eqt_b32/eqtime_L*_U*.csv     (L=8 and L=10, full U grid)
#     ~/results/chi_L12/eqtime_L12_U*.csv      (L=12)
# All three sizes share beta=32, tau window 0.8, N_w=500, six seeds, and a complete
# 7x7 (U, delta) grid. The script asserts that before merging them.
#
# Red: d_xy is the more enhanced channel. Blue: d_x2-y2. Black line: they are equal.
# One shared symmetric colour scale across panels, so L=8 reads paler because its
# signal really is ~2.3x weaker, not because of per-panel rescaling.
import os, glob, numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib
from scipy.interpolate import RegularGridInterpolator

CH   = ["son", "sext", "d", "dxy"]
DIRS = ["~/collected/eqt_b32", "~/results/chi_L12", "~/eqt_b32", "~/chi_L12",
        "../data/chi_fss_b32", "../data/chi_L12", "chi_fss_b32", "chi_L12"]

files = []
for d in DIRS:
    files += glob.glob(os.path.join(os.path.expanduser(d), "eqtime_L*_U*.csv"))
if not files:
    raise SystemExit("no eqtime_L*_U*.csv found. Try:  find ~ -name 'eqtime_L*_U*.csv'")
raw = pd.concat([pd.read_csv(f) for f in sorted(set(files))], ignore_index=True)
raw = raw.drop_duplicates(subset=["L", "U", "delta", "seed"])
print(f"{len(set(files))} files, {len(raw)} rows, L = {sorted(raw.L.unique())}")

for c, want in (("beta_proj", 32.0), ("tau_max", 0.8), ("nw", 500)):
    got = sorted(raw[c].unique())
    assert got == [want], f"{c} differs across the files: {got}. Do not merge these."
print(f"settings uniform: beta=32, tau_max=0.8, nw=500, "
      f"seeds={raw.seed.nunique()}, U={sorted(raw.U.unique())}")

rows = []
for (L, U, dl), s in raw.groupby(["L", "U", "delta"]):
    r = dict(L=int(L), U=U, delta=dl, nseed=len(s))
    for c in CH:
        v = s[f"chi_{c}_vertex"].to_numpy(float)
        r[c] = v.mean()
        r[f"{c}_err"] = v.std(ddof=1)/np.sqrt(len(v)) if len(v) > 1 else np.nan
    rows.append(r)
S = pd.DataFrame(rows).sort_values(["L","U","delta"]).reset_index(drop=True)
S = S[S.U >= 2]
LS, UMAX = sorted(S.L.unique()), S.U.max()
for L in LS:
    p = S[S.L == L].pivot_table(index="U", columns="delta", values="d")
    assert not p.isna().any().any(), f"L={L}: incomplete (U, delta) grid"
    print(f"  L={L}: {p.shape[0]}x{p.shape[1]} grid, complete")

plt.rcParams.update({"font.family":"serif","mathtext.fontset":"dejavuserif",
    "font.size":15,"axes.labelsize":20,"xtick.labelsize":16,"ytick.labelsize":16,
    "axes.linewidth":1.4})
gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,UMAX,400))

def grid(s, col):
    """Bilinear on the measured rectangle. NOT scipy.griddata: that triangulates,
    and on a rectangular lattice the triangulation is neither unique nor stable
    under adding rows, which moves contours far from the new data."""
    p = s.pivot_table(index="U", columns="delta", values=col)
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method="linear", bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

def nice_step(span, target=6):
    raw_ = span / target; mag = 10.0 ** np.floor(np.log10(raw_))
    for m in (1, 2, 2.5, 5, 10):
        if raw_ <= m * mag: return m * mag
    return 10 * mag

Z = {L: grid(S[S.L == L].assign(_v=lambda t: t.dxy - t.d), "_v") for L in LS}
M = max(np.abs(z).max() for z in Z.values())
step = nice_step(2 * M)
ticks = np.arange(np.ceil(-M/step)*step, M + step/2, step)
ticks = ticks[(ticks >= -M) & (ticks <= M)]

fig, ax = plt.subplots(1, len(LS), figsize=(5.6*len(LS), 5.2), sharex=True,
                       sharey=True, squeeze=False, gridspec_kw={"wspace":.09})
ax = ax[0]
for a, L in zip(ax, LS):
    im = a.contourf(gx, gy, Z[L], levels=np.linspace(-M, M, 81), cmap="RdBu_r",
                    vmin=-M, vmax=M)
    a.contour(gx, gy, Z[L], levels=[0.], colors="k", linewidths=3.2)
    a.set_title(rf"$L={int(L)}$", fontsize=21, pad=10)
    a.set_xlabel(r"anisotropy  $\delta$", labelpad=8)
    a.set_xticks([.2,.3,.4,.5,.6]); a.set_xlim(.1,.7); a.set_ylim(2., UMAX)
    a.tick_params(direction="out", top=False, right=False, length=6, width=1.4)
ax[0].set_ylabel(r"$U/t$", labelpad=8)
cb = fig.colorbar(im, ax=ax.tolist(), pad=.014, fraction=.019, ticks=ticks)
cb.set_label(r"$d_{xy}-d_{x^2-y^2}$", fontsize=20, labelpad=10)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.0f"))
cb.ax.tick_params(labelsize=15)
fig.suptitle(r"time-integrated vertex, half filling."
             r"   black line: $d_{xy}=d_{x^2-y^2}$", fontsize=18, y=1.00)
plt.show()

def zc(x, y):
    s = np.sign(y); k = np.where(s[:-1]*s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]; return x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])
print("\ncrossover delta, from the measured grid (no interpolation):")
print(f"{'U':>5}" + "".join(f"{'L='+str(L):>9}" for L in LS))
A = []
for U in sorted(S.U.unique()):
    r = [zc(g.delta.values, (g.dxy - g.d).values)
         for g in (S[(S.L == L) & (S.U == U)].sort_values("delta") for L in LS)]
    A.append(r); print(f"{U:>5.1f}" + "".join(f"{x:>9.3f}" for x in r))
A = np.array(A)
print("  spread over U" + "".join(f"{np.nanmax(A[:,j])-np.nanmin(A[:,j]):>9.3f}"
                                  for j in range(A.shape[1])))
