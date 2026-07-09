#!/usr/bin/env python3
"""Path-A analysis: longitudinal-vs-transverse spin ANISOTROPY across the
combined-altermagnet phase diagram, straight from the EXISTING production
CPQMC output (no new simulations).

Your Fortran (mc2duph.f90) already writes, per run folder, the equal-time
momentum structure factors in dir-kVals/*.dat  (cols: k_x  k_y  value  error):

    sdwz.dat  ->  S^z(q)     longitudinal (the channel your phase diagram uses)
    sdwx.dat  ->  S^x(q)     transverse   (same 0.25 prefactor / Wick structure)
    pmdf.dat  ->  S^{+-}(q)  spin-flip (transverse) correlator

Because sdwx and sdwz are built symmetrically, their ratio is pinned to 1 in the
SU(2)-isotropic limit (tA = t' = 0) and departs when the spin-dependent hopping
breaks spin isotropy.  We therefore define, at the dominant magnetic ordering
wavevector Q* (= argmax of S^z(q), consistent with the phase-diagram classifier):

    R_aniso(point) = S^x(Q*) / S^z(Q*)          (also reports S^{+-}(Q*)/2 S^z(Q*))

    R < 1  -> longitudinal / Ising-Neel  (t_A-driven)
    R ~ 1  -> isotropic (SU(2) intact)
    R > 1  -> transverse / easy-plane / spiral  (t'-driven)

This is the equal-time analogue of the chi_+-/2 chi_zz "meter" from the demo,
computed for the FULL L=8/10/14 grid you already have on the cluster.

Output: one 2x3 phase-diagram figure (panels = t'), points colored by R_aniso and
shaped by ordering-wavevector class, plus a CSV of every (n, tA, t', Q*, R).

    python3 plot_paper2_anisotropy.py --base /home/phd25imran/CPQMC/Imran/L=14 \
        --L 14 --U 4.0 --out-dir /home/phd25imran/CPQMC/Imran/L=14
"""
from __future__ import annotations
import os, re, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D

PI = np.pi
FOLDER_RE = re.compile(
    r'L(?P<L>\d+)n(?P<n>[\d.]+)u(?P<u>[\d.]+)tA(?P<tA>[\d.]+)tt(?P<tt>[\d.]+)N(?P<N>\d+)')


# ----------------------------------------------------------------------------- IO
def read_kfile(path):
    """Read a dir-kVals/*.dat file -> DataFrame[kx, ky, val, err]. None on failure."""
    try:
        d = pd.read_csv(path, skiprows=1, header=None, sep=r'\s+')
        if d.empty or d.shape[1] < 3:
            return None
        d = d.iloc[:, :4] if d.shape[1] >= 4 else d.iloc[:, :3]
        d.columns = ["kx", "ky", "val", "err"][:d.shape[1]]
        return d
    except Exception:
        return None


def find_folders(base, L, U):
    out = []
    for folder in os.listdir(base):
        m = FOLDER_RE.match(folder)
        if not m:
            continue
        if int(m.group('L')) != L or abs(float(m.group('u')) - U) > 1e-6:
            continue
        kdir = os.path.join(base, folder, "dir-kVals")
        if not os.path.isfile(os.path.join(kdir, "sdwz.dat")):
            continue
        out.append(dict(folder=folder, n=float(m.group('n')),
                        tA=float(m.group('tA')), tt=float(m.group('tt')), kdir=kdir))
    out.sort(key=lambda x: (x['tt'], x['tA'], x['n']))
    return out


# ------------------------------------------------------------------- ordering-Q class
def classify_Q(kx, ky, tol=0.25, delta=0.1):
    """Label the ordering wavevector: M=(pi,pi) Neel, X=(pi,0)/(0,pi) stripe,
    N=(pi,q)/(q,pi) spiral, Q=(q,q) diagonal.  Returns (label, marker, q_frac)."""
    dM = min(np.hypot(kx - sx, ky - sy) for sx in (PI, -PI) for sy in (PI, -PI))
    if dM < tol:
        return "Neel", "o", None
    dX = min(np.hypot(kx - sx, ky - sy) for sx, sy in ((PI, 0), (-PI, 0), (0, PI), (0, -PI)))
    if dX < tol:
        return "stripe", "s", None
    if abs(abs(ky) - PI) < tol and delta < abs(kx) < PI - delta:
        return "spiral", "^", abs(kx) / PI
    if abs(abs(kx) - PI) < tol and delta < abs(ky) < PI - delta:
        return "spiral", "^", abs(ky) / PI
    if abs(kx - ky) < tol and abs(kx) < PI - delta:
        return "diagonal", "*", abs(kx) / PI
    return "other", "D", None


# ------------------------------------------------------------------------- analysis
def analyze_point(kdir):
    """Peak of S^z(q) -> Q*; evaluate S^x and S^{+-} at the SAME Q*."""
    sz = read_kfile(os.path.join(kdir, "sdwz.dat"))
    if sz is None:
        return None
    ip = sz["val"].idxmax()
    kx, ky = float(sz.at[ip, "kx"]), float(sz.at[ip, "ky"])
    sz_peak = float(sz.at[ip, "val"])

    def at_Q(fname):
        d = read_kfile(os.path.join(kdir, fname))
        if d is None:
            return np.nan
        j = ((d["kx"] - kx) ** 2 + (d["ky"] - ky) ** 2).idxmin()
        return float(d.at[j, "val"])

    sx_peak = at_Q("sdwx.dat")
    pm_peak = at_Q("pmdf.dat")
    label, marker, qfrac = classify_Q(kx, ky)
    r_x = sx_peak / sz_peak if sz_peak else np.nan          # S^x / S^z  (->1 isotropic)
    r_pm = pm_peak / (2 * sz_peak) if sz_peak else np.nan   # S^{+-}/2S^z (->1 isotropic)
    return dict(qx=kx, qy=ky, qx_pi=kx / PI, qy_pi=ky / PI,
                sz=sz_peak, sx=sx_peak, pm=pm_peak,
                R=r_x, R_pm=r_pm, order=label, marker=marker, qfrac=qfrac)


# ---------------------------------------------------------------------------- driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/home/phd25imran/CPQMC/Imran/L=14",
                    help="folder holding the L..n..u..tA..tt..N.. run directories")
    ap.add_argument("--L", type=int, default=14)
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--tt", type=float, nargs="*",
                    default=[0.0, 0.1, 0.3, 0.5, 0.7, 0.9], help="t' panels to show")
    ap.add_argument("--ratio", choices=["x", "pm"], default="x",
                    help="x: S^x/S^z (default);  pm: S^{+-}/2S^z")
    ap.add_argument("--vmax", type=float, default=2.0, help="upper end of the ratio color scale")
    ap.add_argument("--out-dir", default=".")
    a = ap.parse_args()

    folders = find_folders(a.base, a.L, a.U)
    if not folders:
        raise SystemExit(f"no matching folders in {a.base} for L={a.L}, U={a.U}")
    print(f"found {len(folders)} folders")

    rows = []
    for f in folders:
        if f["tt"] not in a.tt:
            continue
        res = analyze_point(f["kdir"])
        if res is None:
            continue
        rows.append({**{k: f[k] for k in ("n", "tA", "tt", "folder")}, **res})
    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("no analyzable points (missing sdwz/sdwx?)")
    rcol = "R" if a.ratio == "x" else "R_pm"
    rlabel = (r"$S^x(Q^*)/S^z(Q^*)$" if a.ratio == "x"
              else r"$S^{+-}(Q^*)/2\,S^z(Q^*)$")

    os.makedirs(a.out_dir, exist_ok=True)
    csv = os.path.join(a.out_dir, f"anisotropy_L{a.L}_U{a.U}.csv")
    df.to_csv(csv, index=False)
    print(f"wrote {csv}")

    # ---- phase-diagram figure: color = anisotropy ratio, shape = ordering class ----
    tt_vals = [t for t in sorted(df["tt"].unique()) if t in a.tt]
    ncol = 3
    nrow = int(np.ceil(len(tt_vals) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 4.3 * nrow),
                             sharex=True, sharey=True, facecolor="white", squeeze=False)
    norm = mcolors.TwoSlopeNorm(vmin=0.0, vcenter=1.0, vmax=a.vmax)
    cmap = plt.cm.RdYlGn                       # red<1 (Ising/Neel), green>1 (transverse)
    panel = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
    tA_all = sorted(df["tA"].unique()); n_all = sorted(df["n"].unique())

    for i, tt in enumerate(tt_vals):
        ax = axes[i // ncol][i % ncol]
        sub = df[df["tt"] == tt]
        for _, r in sub.iterrows():
            ax.scatter(r["n"], r["tA"], marker=r["marker"],
                       c=[cmap(norm(np.clip(r[rcol], 0, a.vmax)))],
                       s=150, edgecolors="black", linewidth=0.7, zorder=3)
        ax.set_title(f"{panel[i]}  $t' = {tt:.1f}$", fontsize=18)
        ax.set_xlim(min(n_all) - 0.03, max(n_all) + 0.03)
        ax.set_ylim(min(tA_all) - 0.05, max(tA_all) + 0.05)
        ax.tick_params(labelsize=12)
        if i % ncol == 0:
            ax.set_ylabel(r"Anisotropy $t_A$", fontsize=16)
        if i // ncol == nrow - 1:
            ax.set_xlabel(r"Filling $n$", fontsize=16)
    for j in range(len(tt_vals), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")

    fig.subplots_adjust(left=0.07, right=0.86, top=0.90, bottom=0.09, wspace=0.12, hspace=0.28)
    cax = fig.add_axes([0.88, 0.30, 0.017, 0.45])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label(rlabel + "\n(<1 longitudinal / Neel,  >1 transverse)", fontsize=13)
    cb.ax.axhline(1.0, color="k", lw=1.2)
    cb.ax.tick_params(labelsize=11)

    handles = [Line2D([0], [0], marker=m, color="w", markerfacecolor="0.6",
                      markeredgecolor="k", markersize=12, label=lab)
               for lab, m in [(r"$(\pi,\pi)$ Neel", "o"), (r"$(\pi,0)$ stripe", "s"),
                              (r"$(\pi,q)$ spiral", "^"), (r"$(q,q)$ diag", "*")]]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.865, 0.92),
               title="ordering $Q^*$", frameon=True, fontsize=11, title_fontsize=12)
    fig.suptitle(f"Spin anisotropy across the combined-altermagnet phase diagram "
                 f"(L={a.L}, U={a.U}) -- equal-time structure factors", fontsize=15, y=0.965)

    out = os.path.join(a.out_dir, f"anisotropy_phase_L{a.L}_U{a.U}.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"wrote {out}")

    # ---- console summary ----
    print(f"\n{'t':>5} {'tA':>5} {'n':>6} {'order':>9} {'Q*/pi':>14} "
          f"{'S^z':>8} {'S^x':>8} {rcol:>7}")
    for _, r in df.sort_values(["tt", "tA", "n"]).iterrows():
        print(f"{r['tt']:>5.1f} {r['tA']:>5.2f} {r['n']:>6.3f} {r['order']:>9} "
              f"({r['qx_pi']:+.2f},{r['qy_pi']:+.2f}) {r['sz']:>8.4f} {r['sx']:>8.4f} "
              f"{r[rcol]:>7.2f}")


if __name__ == "__main__":
    main()
