"""Does the Python equal-time vertex agree with the Fortran equal-time?

The raw eqtime_L12_U*.csv carry BOTH eq_*_vertex (the tau=0 slice) and
chi_*_vertex (its trapezoidal integral over the tau window). Comparing the Python
tau=0 slice against the Fortran columns tells us whether the equal-time to
integrated difference seen in the two-panel figure is physical, or an artefact of
two codes normalising differently.
"""
import os, glob, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
r = pd.concat([pd.read_csv(x) for x in sorted(glob.glob(os.path.join(D, "chi_L12", "eqtime_L12_U*.csv")))],
              ignore_index=True)
print("raw chi_L12 files: U =", sorted(r.U.unique()), " rows =", len(r))
print("seeds per U:", r.groupby("U").size().to_dict())

g = r.groupby(["U","delta"]).agg(eq_d=("eq_d_vertex","mean"), eq_dxy=("eq_dxy_vertex","mean"),
                                 chi_d=("chi_d_vertex","mean"), chi_dxy=("chi_dxy_vertex","mean"),
                                 nseed=("seed","size")).reset_index()

f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n,1.0) & (f.L==12)][["U","delta","d","dxy"]]
m = g.merge(f, on=["U","delta"], suffixes=("","_F"))
m = m[m.U >= 2]
print(f"\noverlapping cells (U>=2): {len(m)}")
for ch, py, fo in (("d_x2-y2","eq_d","d"), ("d_xy","eq_dxy","dxy")):
    r_ = np.corrcoef(m[py], m[fo])[0,1]
    rat = m[py]/m[fo]
    print(f"  {ch:>8}: Python eq vs Fortran   r = {r_:+.4f}   "
          f"ratio median {np.median(rat):.3f}  spread {np.percentile(rat,5):.3f}..{np.percentile(rat,95):.3f}")

print("\n=== crossover delta, all three quantities ===")
def zc(x, y):
    s = np.sign(y); k = np.where(s[:-1]*s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]; return x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])

print(f"{'U':>5} {'Fortran eq':>12} {'Python eq':>12} {'Python int':>12}")
for U in sorted(m.U.unique()):
    s = m[m.U==U].sort_values("delta"); x = s.delta.values
    print(f"{U:>5.1f} {zc(x,(s.dxy-s.d).values):>12.3f} "
          f"{zc(x,(s.eq_dxy-s.eq_d).values):>12.3f} {zc(x,(s.chi_dxy-s.chi_d).values):>12.3f}")
