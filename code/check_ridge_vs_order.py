"""Why does the Delta_tot ridge run diagonally in the (U, delta) plane?

Delta_tot needs two things at once and they pull against each other:

  MAGNETIC ORDER. With no ordered moment there is no spin splitting to resolve.
  Order is favoured by large U. Anisotropy frustrates the Neel state, so the
  interaction needed to sustain order grows with delta.

  ANISOTROPY. At delta = 0 the two sublattices are related by a translation, the
  splitting is forbidden by symmetry, and Delta_tot vanishes identically. The
  splitting can only grow with delta.

So Delta_tot is largest where delta is as large as the order can tolerate. If that
is right, the ridge should sit at or just inside the boundary of the ordered region
measured independently from the structure-factor scaling.

This compares the two, both from data: the ridge from Delta_tot, the order boundary
from S(Q)/N finite-size scaling in neel_fss_summary.csv.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
f = pd.read_csv(f"{D}/fortran_L8_L10_L12_dedup.csv")
f = f[np.isclose(f.n, 1.0) & (f.U >= 2)]
fss = pd.read_csv(f"{D}/neel_fss_summary.csv")

def peak(x, y):
    i = int(np.argmax(y))
    if i in (0, len(y)-1): return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + .5*(y[i-1]-y[i+1])/den*(x[i+1]-x[i])

def boundary(U, thresh=0.70):
    """Largest delta at that U whose S(Q)/N ratio still exceeds `thresh`,
    i.e. the outer edge of the ordered region. Interpolated to the first delta
    that falls below it."""
    g = fss[np.isclose(fss.U, U)].sort_values("delta")
    g = g[g.delta > 0]
    if len(g) < 2: return np.nan
    d, r = g.delta.values, g.ratio.values
    ok = r >= thresh
    if not ok.any(): return 0.0
    last = np.where(ok)[0].max()
    if last == len(d) - 1: return d[last]
    d0, d1, r0, r1 = d[last], d[last+1], r[last], r[last+1]
    return d0 + (r0 - thresh) * (d1 - d0) / (r0 - r1) if r0 != r1 else d0

print(f"{'U':>5}{'Dtot ridge (L=8)':>19}{'(L=10)':>10}{'(L=12)':>10}"
      f"{'order edge':>13}{'ridge-edge':>12}")
rows = []
for U in sorted(f.U.unique()):
    rid = []
    for L in (8, 10, 12):
        g = f[(f.L == L) & (f.U == U)].sort_values("delta")
        rid.append(peak(g.delta.values, g.dtot_N.values))
    b = boundary(U)
    rows.append((U, *rid, b))
    print(f"{U:>5.1f}{rid[0]:>19.3f}{rid[1]:>10.3f}{rid[2]:>10.3f}"
          f"{b:>13.3f}{rid[2]-b:>12.3f}")

a = np.array(rows, float)
m = ~np.isnan(a[:, 4])
print(f"\n  correlation of the L=12 ridge with the order edge: "
      f"{np.corrcoef(a[m,3], a[m,4])[0,1]:+.3f}   (n = {m.sum()})")
print(f"  mean offset ridge - edge: {np.mean(a[m,3]-a[m,4]):+.3f}"
      f"   rms {np.sqrt(np.mean((a[m,3]-a[m,4])**2)):.3f}")
print(f"\n  both quantities move outward together as U rises:")
print(f"    ridge  {a[0,3]:.3f} -> {a[-1,3]:.3f}")
print(f"    edge   {a[m,4][0]:.3f} -> {a[m,4][-1]:.3f}")
