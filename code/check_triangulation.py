"""Why does adding a U=5 row change the interpolated field at U=3?

Linear interpolation on a triangulation is local: a point inside a triangle
depends only on that triangle's three corners. So adding data far away should
change nothing. It does. This finds out why.
"""
import os, numpy as np, pandas as pd
from scipy.interpolate import griddata
from scipy.spatial import Delaunay

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.L == 12) & (f.U >= 2)]
A = f[["delta", "U"]].values
B = f[f.U <= 4.5][["delta", "U"]].values
print(f"points: with U=5 {len(A)}, without {len(B)}")
print(f"U grid: {sorted(f.U.unique())}   (spacing 1.0, 0.5, 0.5, 0.5, 0.5 -- NOT uniform)")
print(f"delta grid: {sorted(f.delta.unique())}  (uniform 0.1)\n")

tA, tB = Delaunay(A), Delaunay(B)
def tris(t, P):
    return {tuple(sorted(map(tuple, np.round(P[s], 6)))) for s in t.simplices}
sA, sB = tris(tA, A), tris(tB, B)
low = lambda S: {t for t in S if max(p[1] for p in t) <= 4.5}
print(f"triangles entirely below U=4.5:  with U=5 row {len(low(sA))},  without {len(low(sB))}")
print(f"  identical between the two: {len(low(sA) & low(sB))}")
print(f"  present only WITH the U=5 row: {len(low(sA) - low(sB))}")
print(f"  present only WITHOUT it:      {len(low(sB) - low(sA))}")

print("\nThe grid is a rectangular lattice, so every cell can be split into two")
print("triangles along either diagonal and both are valid Delaunay solutions.")
print("Qhull's tie-break is not stable when a point is added, so diagonals flip")
print("far from the change. That is the whole effect.\n")

print("Check: does rescaling the axes before triangulating stabilise it?")
gx, gy = np.meshgrid(np.linspace(.1,.7,120), np.linspace(2.,4.5,120))
for rs in (False, True):
    out = []
    for col in ("d", "dxy", "dtot_N"):
        a = griddata((f.delta, f.U), f[col], (gx,gy), "linear", rescale=rs)
        g = f[f.U <= 4.5]
        b = griddata((g.delta, g.U), g[col], (gx,gy), "linear", rescale=rs)
        m = ~(np.isnan(a) | np.isnan(b))
        rng = np.nanmax(a) - np.nanmin(a)
        out.append(f"{col} {100*np.max(np.abs(a-b)[m])/rng:5.1f}%")
    print(f"  rescale={str(rs):>5}:  max |difference| as % of range   " + "   ".join(out))
