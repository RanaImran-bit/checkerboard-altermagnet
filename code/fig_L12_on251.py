"""L=12 half-filling grid, straight from the Fortran run folders.

Run it on a node and it rescans dir-kVals every time, so rerunning after more
jobs land picks them up with no bookkeeping.

    python fig_L12_on251.py

A run is used only when nohup.out ends with the timing line AND both n_up.dat
and n_dn.dat exist. Partial runs write valid-looking files with unconverged
statistics, so they are skipped rather than averaged in.

NOTE ON COVERAGE: each node holds only its own share of the grid, so this sees
~6 of the 28 finished cells when run on 251 alone. To plot all of them, gather
the nodes into one CSV first and point CSV_IN at it:

    for h in 250 251 252 253 254; do ssh $h 'python ~/fig_L12_on251.py --dump'; done > ~/L12_all.csv

Delta_tot = sum_k |n_up(k) - n_dn(k)| / L^2, PRL Eq. (2). The absolute value is
inside the sum, so it stays finite even though the net magnetisation vanishes.
Pairing numbers are the vertex contributions at k = (0,0).
"""
import os, re, sys, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib as mpl, matplotlib.pyplot as plt

BASE   = os.path.expanduser("~/Checkerboard_Model")
OUTDIR = os.path.expanduser("~/analysis")
CSV_IN = os.path.expanduser("~/L12_all.csv")   # used if it exists, else scan BASE
CHAN   = {"sowave": "son", "swave": "sext", "dwave": "d", "dd12wave": "dxy"}


def scan(base):
    rows = []
    for d in sorted(glob.glob(f"{base}/L12n*/")):
        nm = os.path.basename(d.rstrip("/"))
        m = re.match(r"L(\d+)n([0-9.]+)u([0-9.]+)tA-([0-9.]+)tt", nm)
        if not m:
            continue
        try:
            if "time_sec" not in "".join(open(d + "nohup.out").readlines()[-3:]):
                continue
        except OSError:
            continue
        k = d + "dir-kVals/"
        try:
            up = pd.read_csv(k + "n_up.dat", sep=r"\s+", skiprows=1, header=None,
                             names=["kx", "ky", "v", "e"])
            dn = pd.read_csv(k + "n_dn.dat", sep=r"\s+", skiprows=1, header=None,
                             names=["kx", "ky", "v", "e"])
        except OSError:
            continue
        L = int(m.group(1))
        r = dict(L=L, n=float(m.group(2)), U=float(m.group(3)),
                 delta=float(m.group(4)),
                 dtot_N=float(np.abs(up.v.values - dn.v.values).sum()) / L**2)
        for f, key in CHAN.items():
            try:
                with open(f"{k}Vertex_{f}.dat") as fh:
                    fh.readline()
                    r[key] = float(fh.readline().split()[2])
            except (OSError, IndexError, ValueError):
                r[key] = np.nan
        rows.append(r)
    return pd.DataFrame(rows)


if "--dump" in sys.argv:                       # for gathering nodes into one CSV
    print(scan(BASE).to_csv(index=False, header=False), end="")
    sys.exit()

if os.path.exists(CSV_IN):
    d = pd.read_csv(CSV_IN, header=None,
                    names=["L", "n", "U", "delta", "dtot_N", "son", "sext", "d", "dxy"])
    print(f"read {CSV_IN}")
else:
    d = scan(BASE)
    print(f"scanned {BASE}")
d = d[np.isclose(d.n, 1.0)].sort_values(["U", "delta"])
print(f"{len(d)} complete cells")

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["DejaVu Serif"],
    "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "legend.frameon": False,
    "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12})

US, DS = sorted(d.U.unique()), sorted(d.delta.unique())
CU = plt.cm.viridis(np.linspace(0.05, 0.88, len(US)))
CD = plt.cm.plasma(np.linspace(0.05, 0.82, len(DS)))
MK = ["o", "s", "^", "D", "v", "P", "X"]

# ===================== FIGURE 1: Delta_tot on the x-axis =====================
fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.9))

for i, U in enumerate(US):                             # (a) delta against Dtot
    s = d[d.U == U].sort_values("delta")
    ax[0].plot(s.dtot_N, s.delta, marker=MK[i % len(MK)], ms=6.5, lw=1.6,
               color=CU[i], label=f"$U={U:g}$")
ax[0].set_xlabel(r"$\Delta_{\mathrm{tot}}$")
ax[0].set_ylabel(r"anisotropy  $\delta$")
ax[0].legend(fontsize=10, ncol=2, loc="best")
ax[0].set_title("(a)", loc="left", fontsize=15)

for i, dl in enumerate(DS):                                # (b) U against Dtot
    s = d[d.delta == dl].sort_values("U")
    if len(s) < 2:
        continue
    ax[1].plot(s.dtot_N, s.U, marker=MK[i % len(MK)], ms=6.5, lw=1.6,
               color=CD[i], label=rf"$\delta={dl:g}$")
ax[1].set_xlabel(r"$\Delta_{\mathrm{tot}}$"); ax[1].set_ylabel(r"$U/t$")
ax[1].legend(fontsize=10, ncol=2, loc="best")
ax[1].set_title("(b)", loc="left", fontsize=15)

fig.tight_layout()
os.makedirs(OUTDIR, exist_ok=True)
fig.savefig(f"{OUTDIR}/fig_L12_dtot.png", dpi=600, bbox_inches="tight",
            facecolor="white")

# ================= FIGURE 2: pairing, blue dx2-y2 vs red dxy =================
fig2, ax2 = plt.subplots(1, 2, figsize=(11.6, 4.9))

for i, U in enumerate([2.0, 3.0, 3.5, 4.0]):
    s = d[d.U == U].sort_values("delta")
    if len(s) < 2:
        continue
    al = 0.35 + 0.65 * i / 3
    ax2[0].plot(s.delta, s.d, "-o", ms=5.5, lw=1.6, color="tab:blue", alpha=al,
                label=rf"$d_{{x^2-y^2}}$, $U={U:g}$")
    ax2[0].plot(s.delta, s.dxy, "-s", ms=5.5, lw=1.6, color="tab:red", alpha=al,
                label=rf"$d_{{xy}}$, $U={U:g}$")
ax2[0].axhline(0, color="k", lw=0.8, ls=":")
ax2[0].set_xlabel(r"anisotropy  $\delta$")
ax2[0].set_ylabel(r"pairing vertex at $\mathbf{k}=0$")
ax2[0].legend(fontsize=8.5, ncol=2, loc="best")
ax2[0].set_title("(a)", loc="left", fontsize=15)

w = d[d.U > 0].copy()
w["dif"] = w.dxy - w.d
v = np.abs(w.dif).max()
sc = ax2[1].scatter(w.delta, w.U, c=w.dif, cmap="coolwarm", vmin=-v, vmax=v,
                    s=280, marker="s", edgecolor="k", linewidth=0.7)
cb = fig2.colorbar(sc, ax=ax2[1], pad=0.02)
cb.set_label(r"$d_{xy} - d_{x^2-y^2}$", fontsize=13)
ax2[1].set_xlabel(r"anisotropy  $\delta$"); ax2[1].set_ylabel(r"$U/t$")
ax2[1].set_xlim(0.03, 0.67); ax2[1].set_ylim(1.5, 5.5)
ax2[1].set_title(r"(b)  red = $d_{xy}$ wins", loc="left", fontsize=15)

fig2.tight_layout()
fig2.savefig(f"{OUTDIR}/fig_L12_pairing.png", dpi=600, bbox_inches="tight",
             facecolor="white")

print(f"wrote {OUTDIR}/fig_L12_dtot.png and {OUTDIR}/fig_L12_pairing.png\n")
print("Delta_tot (rows U, cols delta):")
print(d.pivot_table(index="U", columns="delta", values="dtot_N").round(4).to_string())
