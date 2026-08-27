"""Finite-size behaviour of the vertex, now that all three sizes exist.

Until today the vertex was thought to be an L=12-only result. The L=8 and L=10 runs
on the same grid and settings were sitting uncollected on 251. This asks the two
questions that matter: does the d_xy / d_x2-y2 crossover sit at the same delta at
every size, and does the d_xy vertex still grow with system size.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
S = pd.read_csv(os.path.join(D, "chi_vertex_summary_3L.csv"))
S = S[S.U >= 2]

def zc(x, y):
    s = np.sign(y); k = np.where(s[:-1]*s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]; return x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])

print("=== crossover delta where d_xy overtakes d_x2-y2 ===")
print(f"{'U':>5}" + "".join(f"{'L='+str(L):>10}" for L in sorted(S.L.unique())))
tab = {}
for U in sorted(S.U.unique()):
    row = []
    for L in sorted(S.L.unique()):
        g = S[(S.L == L) & (S.U == U)].sort_values("delta")
        row.append(zc(g.delta.values, (g.dxy - g.d).values))
    tab[U] = row
    print(f"{U:>5.1f}" + "".join(f"{x:>10.3f}" for x in row))
A = np.array(list(tab.values()))
print("\n  mean over U   " + "".join(f"{np.nanmean(A[:,j]):>10.3f}" for j in range(A.shape[1])))
print("  spread over U " + "".join(f"{np.nanmax(A[:,j])-np.nanmin(A[:,j]):>10.3f}" for j in range(A.shape[1])))
print(f"\n  spread ACROSS sizes at fixed U: {np.nanmax(np.nanmax(A,1)-np.nanmin(A,1)):.3f} at most")

print("\n=== does the d_xy vertex grow with size? peak value at each U ===")
print(f"{'U':>5}" + "".join(f"{'L='+str(L):>12}" for L in sorted(S.L.unique())) + f"{'L12/L8':>10}")
for U in sorted(S.U.unique()):
    vals = []
    for L in sorted(S.L.unique()):
        g = S[(S.L == L) & (S.U == U)]
        vals.append(g.dxy.max())
    print(f"{U:>5.1f}" + "".join(f"{v:>12.3f}" for v in vals) + f"{vals[-1]/vals[0]:>10.2f}")

print("\n=== where does d_xy peak in delta, per size ===")
for L in sorted(S.L.unique()):
    g = S[S.L == L]
    pk = {float(U): float(t.loc[t.dxy.idxmax(), "delta"]) for U, t in g.groupby("U")}
    print(f"  L={int(L):>3}: " + ", ".join(f"U={U:.1f}->{d:.1f}" for U, d in pk.items()))

print("\n=== signs, all cells U>=2 ===")
for L in sorted(S.L.unique()):
    g = S[S.L == L]
    print(f"  L={int(L):>3}: d_x2-y2 negative in {(g.d<0).sum():>2}/{len(g)}, "
          f"d_xy negative in {(g.dxy<0).sum():>2}/{len(g)}"
          + (f" (delta <= {g[g.dxy<0].delta.max():.1f})" if (g.dxy<0).any() else ""))
