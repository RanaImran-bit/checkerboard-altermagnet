"""Why do the two L=12 equal-time panels look different?

The three-size script sets UMAX = 5.0. The two-panel script sets
UMAX = min(fortran, vertex) = 4.5, because the vertex data stops there. Same CSV,
same channels, same contour levels. This checks two things:

  1. how much of the difference is simply the extra U = 4.5 to 5 band, and
  2. whether the SHARED window 2 to 4.5 changes too, which it can: cubic
     interpolation matches slopes across triangle edges, so adding a row of data
     at U = 5 alters the surface below it as well.
"""
import os, numpy as np, pandas as pd
from scipy.interpolate import griddata

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.L == 12) & (f.U >= 2)]

# evaluate both on the SAME display mesh, the shared window only
gx, gy = np.meshgrid(np.linspace(.1, .7, 400), np.linspace(2., 4.5, 400))

def field(df, col, how):
    P = (df.delta.values, df.U.values)
    a = griddata(P, df[col].values, (gx, gy), how)
    b = griddata(P, df[col].values, (gx, gy), "linear")
    a[np.isnan(a)] = b[np.isnan(a)]
    return a

print("U rows present:")
print(f"  UMAX=5.0 script -> {sorted(f.U.unique())}")
print(f"  UMAX=4.5 script -> {sorted(f[f.U <= 4.5].U.unique())}")

for how in ("cubic", "linear"):
    print(f"\n--- {how} interpolation, compared on the shared window U in [2, 4.5] ---")
    for col in ("d", "dxy", "dtot_N"):
        A = field(f,             col, how)      # U=5 row included
        B = field(f[f.U <= 4.5], col, how)      # U=5 row excluded
        dmax = np.nanmax(np.abs(A - B))
        rng  = np.nanmax(A) - np.nanmin(A)
        print(f"  {col:>7}: max |difference| = {dmax:.4f}   "
              f"({100*dmax/rng:.1f}% of the panel's range)")

print("\n--- and the crossover line itself, on the shared window ---")
def zc(x, y):
    s = np.sign(y); k = np.where(s[:-1]*s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]; return x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])
print(f"{'U':>5}{'delta_c':>10}   (from the measured grid, no interpolation)")
for U in sorted(f[f.U <= 4.5].U.unique()):
    g = f[f.U == U].sort_values("delta")
    print(f"{U:>5.1f}{zc(g.delta.values, (g.dxy-g.d).values):>10.3f}")
print("\nThe measured crossover does not depend on UMAX at all: it is computed")
print("column by column. Only the DRAWN contours shift, and only via interpolation.")


# --- where does the difference actually live? ---
print("\n" + "="*64)
print("WHERE the difference sits (shared window, U in [2, 4.5])")
print("="*64)
for how in ("cubic", "linear"):
    print(f"\n--- {how} ---")
    for col in ("d", "dxy", "dtot_N"):
        A = field(f, col, how); B = field(f[f.U <= 4.5], col, how)
        Dm = np.abs(A - B); rng = np.nanmax(A) - np.nanmin(A)
        big = Dm > 0.02*rng                       # 2% of range
        if not big.any():
            print(f"  {col:>7}: nothing above 2% of range"); continue
        Us = gy[big]
        print(f"  {col:>7}: {100*big.mean():>5.1f}% of the panel above 2% of range,"
              f"  those points span U = {Us.min():.2f} to {Us.max():.2f},"
              f"  median U = {np.median(Us):.2f}")
        # how much is confined to the top band?
        frac_top = (Us >= 4.0).mean()
        print(f"           {100*frac_top:.0f}% of them sit at U >= 4.0"
              f" (the row adjacent to the added U=5 data)")
