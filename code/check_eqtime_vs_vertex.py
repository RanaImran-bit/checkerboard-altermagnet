"""Equal-time pair correlations vs the time-integrated vertex, L=12.

fortran_L8_L10_L12_dedup.csv holds EQUAL-TIME pair correlations. The Fortran
driver computes nothing else (see checkerboard_eqtime.py). chi_L12_vertex_summary.csv
holds the TIME-INTEGRATED susceptibility vertex, the full correlator minus its
uncorrelated Wick piece, at beta = 32, N_w = 500, tau window 0.8, six seeds.

Both live on the same (U, delta) grid at L=12, so the two can be compared cell by
cell and the d_xy / d_x2-y2 crossover can be located in each.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.L == 12)]
v = pd.read_csv(os.path.join(D, "chi_L12_vertex_summary.csv"))

print("=== coverage ===")
print(f"equal-time L=12 : U {sorted(f.U.unique())}")
print(f"vertex     L=12 : U {sorted(v.U.unique())}")
print(f"vertex other L  : only U=4 (chi_fss_L8_U4.csv, chi_fss_L10_U4.csv)"
      f" -> no (U,delta) map at L=8 or 10")

print("\n=== the U=0 null, vertex only ===")
z = v[v.U == 0]
print(f"  d   max |value| = {z.d.abs().max():.2e}")
print(f"  dxy max |value| = {z.dxy.abs().max():.2e}   (vertex vanishes for free electrons)")

def zero_cross(x, y):
    s = np.sign(y); k = np.where(s[:-1]*s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]; return x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])

def peak(x, y):
    i = int(np.argmax(y))
    if i in (0, len(y)-1): return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + .5*(y[i-1]-y[i+1])/den*(x[i+1]-x[i])

print("\n=== where d_xy overtakes d_x2-y2, and the Delta_tot ridge ===")
print(f"{'U':>5} {'eq-time':>10} {'vertex':>10} {'ridge':>10}   {'eq-vtx':>8}")
a, b = [], []
for U in sorted(set(f.U.unique()) & set(v.U.unique()) - {0.0}):
    gf = f[f.U == U].sort_values("delta"); gv = v[v.U == U].sort_values("delta")
    xe = zero_cross(gf.delta.values, (gf.dxy - gf.d).values)
    xv = zero_cross(gv.delta.values, (gv.dxy - gv.d).values)
    xr = peak(gf.delta.values, gf.dtot_N.values)
    print(f"{U:>5.1f} {xe:>10.3f} {xv:>10.3f} {xr:>10.3f}   {xe-xv:>8.3f}")
    if np.isfinite(xe) and np.isfinite(xv): a.append(xe); b.append(xv)
a, b = np.array(a), np.array(b)
print(f"\n  equal-time vs vertex crossover: r = {np.corrcoef(a,b)[0,1]:+.3f}, "
      f"mean offset {np.mean(a-b):+.3f}, rms {np.sqrt(np.mean((a-b)**2)):.3f}, n = {len(a)}")

print("\n=== sign of each channel (U >= 2) ===")
for nm, t in (("equal-time", f[f.U >= 2]), ("vertex", v[v.U >= 2])):
    print(f"  {nm:>11}: d_x2-y2 negative in {(t.d<0).sum():>2}/{len(t)}, "
          f"d_xy negative in {(t.dxy<0).sum():>2}/{len(t)}"
          + (f"  (delta <= {t[t.dxy<0].delta.max():.1f})" if (t.dxy<0).any() else ""))

print("\n=== does the vertex sign survive its own error bar? ===")
w = v[v.U >= 2].copy()
w["d_sig"]   = w.d / w.d_err
w["dxy_sig"] = w.dxy / w.dxy_err
neg = w[w.dxy < 0]
print(f"  d_xy negative cells: {len(neg)}, of which beyond 3 sigma: {(neg.dxy_sig < -3).sum()}")
print(f"  d_x2-y2 cells below 3 sigma from zero: {(w.d/w.d_err < 3).sum()}/{len(w)}")
