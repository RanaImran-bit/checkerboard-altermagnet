"""Finite-size scaling of the Neel structure factor across the half-filling grid.

S(Q)/N is the squared staggered moment. Long-range order means it approaches a
nonzero constant as N grows; short-ranged correlations mean it decays as 1/N.

R_p, the peak-sharpness ratio, cannot make this distinction: it is computed at one
lattice size and normalised by the peak height, so a correlation length merely
comparable to L already gives R_p near 1.

PROVENANCE. Every run directory was verified to be the checkerboard model before
inclusion, not by folder name but by reading its own copy of the source and inputs:
  in.dat line 5 must have t0 < 0, t1 = +0.3, tam = 0 (delta enters through t2)
  mc2duph.f90 must contain the parity branch mod( ixv(i)+iyv(i), 2 ) and tp = t1 + t2
The spin-dependent-hopping model of the sibling paper uses the SAME tA/tt folder
naming but has t1 = -0.3, t2 = 0, tam = -delta and no parity branch. 15 such runs
were found mixed in and are quarantined in data/sdwz_half_rejected/. Folder names
alone cannot separate the two models.

Data: dir-kVals/sdwz.dat, columns k_x k_y value error, one header line. The
convention is value = (1/N) sum_ij e^{iQr} <SiSj>, confirmed by the Brillouin-zone
mean being size independent at 0.184 across L = 12 to 18. So value is extensive and
value/N is the order parameter.
"""
import os, glob, re
import numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "sdwz_half")
PAT = re.compile(r"U([\d.]+)_d(-?[\d.]+)_L(\d+)\.dat")

rows = []
for fp in sorted(glob.glob(os.path.join(D, "*.dat"))):
    m = PAT.search(os.path.basename(fp))
    if not m: continue
    U, d, L = float(m.group(1)), float(m.group(2)), int(m.group(3))
    try: a = np.loadtxt(fp, skiprows=1)
    except Exception: continue
    if a.ndim != 2 or a.shape[1] < 3 or not len(a): continue
    i = int(np.argmax(a[:, 2]))
    rows.append(dict(U=U, delta=d, L=L, N=L*L, kx=a[i,0], ky=a[i,1],
                     S=a[i,2], err=a[i,3] if a.shape[1] > 3 else np.nan,
                     m2=a[i,2]/(L*L)))
t = pd.DataFrame(rows)
print(f"loaded {len(t)} cells, L = {sorted(t.L.unique())}\n")

def fit(g):
    """Least squares of m^2 against 1/L. Returns intercept, its error, and the
    ratio of the largest-L value to the smallest, which 1/N decay would push toward
    (L_small/L_large)^2."""
    g = g.sort_values("L")
    x, y = 1.0/g.L.values, g.m2.values
    A = np.column_stack([np.ones(len(x)), x])
    c, s = np.linalg.lstsq(A, y, rcond=None)[0]
    r = y - (c + s*x); dof = max(len(x)-2, 1)
    e = np.sqrt((r**2).sum()/dof) * np.sqrt(np.linalg.inv(A.T@A)[0,0]) if len(x) > 2 else np.nan
    return c, e, y[-1]/y[0], (g.L.iloc[0]/g.L.iloc[-1])**2

print("Cells with three or more lattice sizes. 'expected' is the ratio pure 1/N")
print("decay would give between the smallest and largest L present.\n")
print(f"{'U':>5}{'delta':>7}{'sizes':>18}{'m2 (largest L)':>16}{'ratio':>8}{'expected':>10}{'intercept':>12}   verdict")
out = []
for (U, d), g in t.groupby(["U", "delta"]):
    g = g[g.L >= 8]
    if len(g) < 3: continue
    c, e, ratio, exp = fit(g)
    g = g.sort_values("L")
    if ratio > 0.80:      v = "ORDER"
    elif ratio < 1.45*exp: v = "1/N decay"
    else:                  v = "intermediate"
    print(f"{U:>5.1f}{d:>7.2f}{','.join(str(x) for x in g.L):>18}"
          f"{g.m2.iloc[-1]:>16.5f}{ratio:>8.3f}{exp:>10.3f}{c:>12.5f}   {v}")
    out.append(dict(U=U, delta=d, nL=len(g), m2_max=g.m2.iloc[-1],
                    ratio=ratio, expected=exp, intercept=c, verdict=v))
o = pd.DataFrame(out)
o.to_csv(os.path.join(D, "..", "neel_fss_summary.csv"), index=False)
print(f"\nwrote data/neel_fss_summary.csv  ({len(o)} cells)")
print(f"  ORDER: {(o.verdict=='ORDER').sum()}   intermediate: {(o.verdict=='intermediate').sum()}"
      f"   1/N decay: {(o.verdict=='1/N decay').sum()}")
