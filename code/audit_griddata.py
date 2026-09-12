"""Which figure scripts are affected by the griddata triangulation problem?

For every (data file, x, y) parameter map used by the figure scripts, check that
the points form a complete rectangle, then measure how much the interpolated field
moves when a row of the scan is dropped. A local interpolant should not move at
all away from the dropped row. griddata does.
"""
import os, glob, numpy as np, pandas as pd
from scipy.interpolate import griddata
from gridinterp import is_regular, reggrid

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

CASES = [
    # (label, csv, filter, xcol, ycol, value columns)
    ("L=12 half filling (U,delta)", "fortran_L8_L10_L12_dedup.csv",
     lambda d: d[np.isclose(d.n,1.0) & (d.L==12) & (d.U>=2)], "delta", "U",
     ["d","dxy","dtot_N","son","sext"]),
    ("L=10 half filling (U,delta)", "fortran_L8_L10_L12_dedup.csv",
     lambda d: d[np.isclose(d.n,1.0) & (d.L==10) & (d.U>=2)], "delta", "U",
     ["d","dxy","dtot_N"]),
    ("L=8 half filling (U,delta)", "fortran_L8_L10_L12_dedup.csv",
     lambda d: d[np.isclose(d.n,1.0) & (d.L==8) & (d.U>=2)], "delta", "U",
     ["d","dxy","dtot_N"]),
    ("chi vertex L=12 (U,delta)", "chi_L12_vertex_summary.csv",
     lambda d: d[d.U>=2], "delta", "U", ["d","dxy","son","sext"]),
]

print(f"{'case':<32}{'col':>8}{'regular':>9}{'griddata drift':>16}{'reggrid drift':>15}")
print("-"*80)
for lab, csv, filt, xc, yc, cols in CASES:
    path = os.path.join(D, csv)
    if not os.path.exists(path):
        print(f"{lab:<32}  MISSING {csv}"); continue
    d = filt(pd.read_csv(path))
    ys = np.sort(d[yc].unique())
    drop = ys[-1]                       # drop the top row, as the U=5 case did
    lo, hi = float(d[xc].min()), float(d[xc].max())
    ylo, yhi = float(ys[0]), float(ys[-2])
    gx, gy = np.meshgrid(np.linspace(lo,hi,220), np.linspace(ylo,yhi,220))
    for c in cols:
        if c not in d.columns: continue
        reg = is_regular(d, xc, yc, c)
        sub = d[d[yc] != drop]
        P, Ps = (d[xc], d[yc]), (sub[xc], sub[yc])
        a = griddata(P, d[c], (gx,gy), "linear")
        b = griddata(Ps, sub[c], (gx,gy), "linear")
        m = ~(np.isnan(a) | np.isnan(b)); rng = np.nanmax(a)-np.nanmin(a)
        g1 = 100*np.max(np.abs(a-b)[m])/rng if rng else 0.0
        if reg:
            a2 = reggrid(d, xc, yc, c, gx, gy); b2 = reggrid(sub, xc, yc, c, gx, gy)
            g2 = 100*np.max(np.abs(a2-b2))/(a2.max()-a2.min())
        else:
            g2 = float("nan")
        print(f"{lab:<32}{c:>8}{str(reg):>9}{g1:>15.1f}%{g2:>14.2f}%")
print("\ndrift = biggest change in the interpolated field, away from the dropped row,")
print("as a percentage of the panel's range. A local interpolant must give 0.")
