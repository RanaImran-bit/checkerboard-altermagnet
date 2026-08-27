"""Does the d_xy / d_x2-y2 crossover follow the ridge of Delta_tot?

For each L and each measured U, find (a) the delta where dxy - d changes sign,
by linear interpolation between the bracketing grid points, and (b) the delta
where Delta_tot/N is largest, refined by a parabola through the peak point and
its two neighbours. Then compare the two curves.
"""
import os, numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
d = pd.read_csv(os.path.join(D, "fortran_L8_L10_L12_dedup.csv"))
d = d[np.isclose(d.n, 1.0) & (d.U >= 2)]

def zero_cross(x, y):
    s = np.sign(y)
    k = np.where(s[:-1] * s[1:] < 0)[0]
    if len(k) != 1: return np.nan
    i = k[0]
    return x[i] - y[i] * (x[i+1] - x[i]) / (y[i+1] - y[i])

def peak(x, y):
    i = int(np.argmax(y))
    if i in (0, len(y)-1): return x[i]
    d0, d1, d2 = y[i-1], y[i], y[i+1]
    den = d0 - 2*d1 + d2
    if den == 0: return x[i]
    return x[i] + 0.5 * (d0 - d2) / den * (x[i+1] - x[i])

if __name__ == "__main__":
    rows = []
    for L in sorted(d.L.unique()):
        s = d[d.L == L]
        print(f"\n=== L = {L} ===")
        print(f"{'U':>5} {'delta_cross':>12} {'delta_peak':>11} {'diff':>8}")
        a, b = [], []
        for U in sorted(s.U.unique()):
            g = s[s.U == U].sort_values("delta")
            x = g.delta.values
            xc = zero_cross(x, (g.dxy - g.d).values)
            xp = peak(x, g.dtot_N.values)
            print(f"{U:>5.1f} {xc:>12.3f} {xp:>11.3f} {xc-xp:>8.3f}"
                  if np.isfinite(xc) else
                  f"{U:>5.1f} {'none':>12} {xp:>11.3f} {'-':>8}")
            if np.isfinite(xc): a.append(xc); b.append(xp); rows.append((L,U,xc,xp))
        a, b = np.array(a), np.array(b)
        if len(a) > 2:
            r = np.corrcoef(a, b)[0,1]
            print(f"  crossover vs ridge: r = {r:+.3f}, "
                  f"mean offset {np.mean(a-b):+.3f}, rms {np.sqrt(np.mean((a-b)**2)):.3f}, "
                  f"n = {len(a)}")

    t = pd.DataFrame(rows, columns=["L","U","delta_cross","delta_peak"])
    t.to_csv(os.path.join(D, "crossover_vs_ridge.csv"), index=False)
    r = np.corrcoef(t.delta_cross, t.delta_peak)[0,1]
    print(f"\npooled over all L: r = {r:+.3f}, n = {len(t)}")
    print(f"wrote {os.path.join(D,'crossover_vs_ridge.csv')}")

    print("\n--- sign of each channel, all cells ---")
    for L in sorted(d.L.unique()):
        s = d[d.L == L]
        print(f"L={L}: d_x2-y2 negative in {(s.d<0).sum():>2}/{len(s)} cells, "
              f"d_xy negative in {(s.dxy<0).sum():>2}/{len(s)} cells "
              f"(delta <= {s[s.dxy<0].delta.max() if (s.dxy<0).any() else float('nan'):.1f})")
