"""The fix: the data is on a regular grid, so stop using scattered-data interpolation.

griddata() triangulates. On a rectangular lattice every cell can be split along
either diagonal, both are valid Delaunay triangulations, and Qhull's choice is not
stable when points are added. That is why adding a U=5 row moved contours at U=3.

RegularGridInterpolator works on the rectangle directly. Bilinear inside a cell
depends only on that cell's four corners, so it is unique, stable, and local.
"""
import os, numpy as np, pandas as pd, scipy
from scipy.interpolate import RegularGridInterpolator

print("scipy", scipy.__version__)
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.L == 12) & (f.U >= 2)]

def rg(df, col, how, gx, gy):
    p = df.pivot_table(index="U", columns="delta", values=col)
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method=how, bounds_error=False, fill_value=None)
    return fn(np.stack([gy.ravel(), gx.ravel()], -1)).reshape(gx.shape)

gx, gy = np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,4.5,400))
for how in ("linear", "cubic"):
    print(f"\n--- RegularGridInterpolator, {how} ---")
    ok = True
    for col in ("d", "dxy", "dtot_N"):
        A = rg(f,             col, how, gx, gy)     # U=5 row included
        B = rg(f[f.U <= 4.5], col, how, gx, gy)     # U=5 row excluded
        dm = np.max(np.abs(A - B)); rng = A.max() - A.min()
        print(f"  {col:>7}: max |difference| = {dm:.3e}  ({100*dm/rng:.3f}% of range)")
        ok &= dm < 1e-12
    print(f"  -> shared window {'IDENTICAL' if ok else 'still moves'} with and without the U=5 row")

print("\n--- and does it stay inside the data, unlike the cubic overshoot? ---")
for how in ("linear", "cubic"):
    A = rg(f, "d", how, *np.meshgrid(np.linspace(.1,.7,400), np.linspace(2.,5.,400)))
    print(f"  {how:>7}: {A.min():+.4f} .. {A.max():+.4f}   "
          f"(measured {f.d.min():+.4f} .. {f.d.max():+.4f})")
