"""Is there one optimal path in (U, delta), and do polarisation and pairing share it?

Three candidate "paths" through the (U, delta) plane, each a delta as a function of U:
  1. the Delta_tot ridge      -- where the magnetic polarisation is largest
  2. the d_xy optimum         -- where the d_xy pairing vertex is largest
  3. the channel crossover    -- where d_xy overtakes d_x2-y2
If magnetism and pairing shared a single optimal path, all three would coincide.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
f = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
f = f[np.isclose(f.n, 1.0) & (f.U >= 2)]
v = pd.read_csv(os.path.join(D, "chi_vertex_summary_3L.csv")); v = v[v.U >= 2]

def peak(x, y):
    i = int(np.argmax(y))
    if i in (0, len(y)-1): return x[i]
    den = y[i-1] - 2*y[i] + y[i+1]
    return x[i] if den == 0 else x[i] + .5*(y[i-1]-y[i+1])/den*(x[i+1]-x[i])

def zc(x, y):
    s = np.sign(y); k = np.where(s[:-1]*s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]; return x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])

for L in (8, 10, 12):
    fs, vs = f[f.L == L], v[v.L == L]
    print(f"\n=== L = {L} ===")
    print(f"{'U':>5}{'Dtot ridge':>13}{'d_xy optimum':>15}{'crossover':>12}")
    R, X, C = [], [], []
    for U in sorted(vs.U.unique()):
        g = fs[fs.U == U].sort_values("delta")
        h = vs[vs.U == U].sort_values("delta")
        r = peak(g.delta.values, g.dtot_N.values)
        x = peak(h.delta.values, h.dxy.values)
        c = zc(h.delta.values, (h.dxy - h.d).values)
        R.append(r); X.append(x); C.append(c)
        print(f"{U:>5.1f}{r:>13.3f}{x:>15.3f}{c:>12.3f}")
    R, X, C = np.array(R), np.array(X), np.array(C)
    print(f"  range over U:  ridge {R.min():.3f}-{R.max():.3f}   "
          f"d_xy opt {X.min():.3f}-{X.max():.3f}   crossover {C.min():.3f}-{C.max():.3f}")
    print(f"  ridge vs d_xy optimum : gap {R[0]-X[0]:+.3f} at U=2  ->  {R[-1]-X[-1]:+.3f} at U=5")
    print(f"  ridge vs crossover    : gap {R[0]-C[0]:+.3f} at U=2  ->  {R[-1]-C[-1]:+.3f} at U=5")

print("\n" + "="*64)
print("DOES POLARISATION PREDICT PAIRING AT FIXED GEOMETRY?")
print("="*64)
def resid(y, X_):
    A = np.column_stack([np.ones(len(y))] + list(X_))
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    return y - A @ b
for L in (8, 10, 12):
    m = v[v.L == L].merge(f[f.L == L][["U","delta","dtot_N"]], on=["U","delta"])
    dl, U = m.delta.values, m.U.values
    rx = np.corrcoef(resid(m.dxy.values, [dl, U]), resid(m.dtot_N.values, [dl, U]))[0,1]
    rd = np.corrcoef(resid(m.d.values,   [dl, U]), resid(m.dtot_N.values, [dl, U]))[0,1]
    print(f"  L={L}: partial corr with Delta_tot, delta and U removed -- "
          f"d_xy {rx:+.3f},  d_x2-y2 {rd:+.3f}   ({len(m)} cells)")
