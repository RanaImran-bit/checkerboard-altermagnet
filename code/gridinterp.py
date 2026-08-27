"""Grid interpolation for the (U, delta) and (n, delta) parameter maps.

Why this module exists
----------------------
These figures were all built with scipy.griddata, which treats the points as
scattered and triangulates them. But the data is a rectangular lattice, and on a
rectangle every cell can be split along either diagonal with both splits equally
valid as a Delaunay triangulation. Qhull breaks the tie arbitrarily and the choice
is NOT stable when points are added: adding a U=5 row to the L=12 half-filling grid
flipped 22 of the 48 triangles lying below U=4.5, moving drawn contours by up to
20% of the panel range at U=3, nowhere near the new data.

A figure whose contours depend on which rows you happened to include is not a
function of the data alone.

RegularGridInterpolator works on the rectangle directly. With method="linear"
(bilinear) each cell depends only on its own four corners, so the result is unique,
local, reproducible, and bounded by the measurements. Verified: identical to machine
precision with and without the U=5 row, and the interpolated range matches the
measured range to 1e-4.

Do not switch to method="cubic" for a paper figure. It is non-local by
construction, so it still moves U=3 contours by 7.9% when U=5 data is added, and it
undershoots to -0.05 where the smallest measurement is +0.02.

See check_umax_effect.py, check_triangulation.py and check_reggrid_fix.py.
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator

__all__ = ["is_regular", "reggrid"]


def is_regular(df, xcol, ycol, vcol):
    """True when (xcol, ycol) form a complete rectangle with one value per cell."""
    p = df.pivot_table(index=ycol, columns=xcol, values=vcol)
    return (not p.isna().any().any()) and len(p) > 1 and p.shape[1] > 1


def reggrid(df, xcol, ycol, vcol, gx, gy, method="linear", fill_value=None):
    """Interpolate vcol onto (gx, gy). gx varies along xcol, gy along ycol.

    fill_value=None extrapolates outside the measured rectangle, matching what the
    figures relied on before. Pass np.nan to blank the outside instead, which some
    scripts use to mask the region beyond the scan.

    Raises if the parameter grid has holes, rather than silently falling back to a
    triangulation. A hole means the scan is incomplete and the figure would be
    inventing the missing cell either way, which is worth knowing about.
    """
    p = df.pivot_table(index=ycol, columns=xcol, values=vcol)
    if p.isna().any().any():
        miss = [(float(p.index[i]), float(p.columns[j]))
                for i, j in zip(*np.where(p.isna().values))]
        raise ValueError(f"{vcol}: grid has {len(miss)} missing cell(s), "
                         f"first few (y, x) = {miss[:5]}")
    fn = RegularGridInterpolator((p.index.values, p.columns.values), p.values,
                                 method=method, bounds_error=False,
                                 fill_value=fill_value)
    return fn(np.stack([np.asarray(gy).ravel(), np.asarray(gx).ravel()], -1)
              ).reshape(np.shape(gx))
