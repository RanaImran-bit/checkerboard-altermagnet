"""Does R_p indicate long-range order? No. Here is the test that does.

R_p = (S(Q) - mean of the neighbouring k points) / S(Q) is a PEAK SHARPNESS ratio,
computed at one lattice size. It says the structure-factor peak stands well above its
immediate neighbours on that lattice's k grid. It cannot distinguish long-range order
from a correlation length comparable to L, because both give a sharp peak at fixed L.
It is also normalised by S(Q) itself, so the peak's magnitude cancels out entirely.

Long-range order is a statement about how the peak WEIGHT survives the thermodynamic
limit: S(Q)/N must approach a nonzero constant as N grows. That needs several lattice
sizes at the same parameters, which these series provide.

  n=1, U=4, delta=0.2 : L = 4, 12, 14, 16, 18
  n=1, U=4, delta=0.3 : L = 12, 14, 16, 18
    (the only L=8 run at delta=0.3 is beta16_..., projected to beta=16 rather than 32,
     so it is not part of the series and is excluded)
"""
import os, glob, re
import numpy as np, pandas as pd

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "sdwz_fss")

def nearest_idx(d, kx, ky, mi, num):
    r = np.sqrt((d[:,0]-kx)**2 + (d[:,1]-ky)**2)
    r[mi] = np.inf
    return np.argsort(r)[:num]

def analyse(fp, L):
    d = np.loadtxt(fp, skiprows=1)
    mi = int(np.argmax(d[:,2]))
    kx, ky, S, err = d[mi,0], d[mi,1], d[mi,2], d[mi,3]
    ex = np.isclose(abs(kx), np.abs(d[:,0]).max())
    ey = np.isclose(abs(ky), np.abs(d[:,1]).max())
    num = 1 if (ex and ey) else (2 if (ex or ey) else 4)
    avg = d[nearest_idx(d, kx, ky, mi, num), 2].mean()
    rp = (S - avg) / S if S else 0.0
    return dict(L=L, N=L*L, kx=kx, ky=ky, S=S, err=err, Rp=rp,
                S_over_N=S/(L*L), npts=len(d))

for tag in ("d0.2", "d0.3"):
    fs = sorted(glob.glob(os.path.join(D, f"{tag}_L*.dat")),
                key=lambda p: int(re.search(r"_L(\d+)", p).group(1)))
    if not fs: continue
    rows = [analyse(p, int(re.search(r"_L(\d+)", p).group(1))) for p in fs]
    t = pd.DataFrame(rows).sort_values("L").reset_index(drop=True)
    print("=" * 74)
    print(f"n = 1, U = 4, delta = {tag[1:]}   ({len(t)} lattice sizes)")
    print("=" * 74)
    print(f"{'L':>4}{'N':>6}{'Q/pi':>14}{'S(Q)':>11}{'error':>10}{'R_p':>8}{'S(Q)/N':>11}")
    for _, r in t.iterrows():
        print(f"{int(r.L):>4}{int(r.N):>6}"
              f"{'(' + format(r.kx/np.pi,'.2f') + ',' + format(r.ky/np.pi,'.2f') + ')':>14}"
              f"{r.S:>11.4f}{r.err:>10.5f}{r.Rp:>8.3f}{r.S_over_N:>11.6f}")

    big = t[t.L >= 12]
    print(f"\n  R_p      over L>=12: {big.Rp.min():.3f} to {big.Rp.max():.3f}"
          f"   -- flat, and it stays high at every size")
    print(f"  S(Q)     over L>=12: {big.S.min():.4f} to {big.S.max():.4f}"
          f"   -- grows with L")
    print(f"  S(Q)/N   over L>=12: {big.S_over_N.min():.6f} to {big.S_over_N.max():.6f}"
          f"   -- falls by {big.S_over_N.iloc[0]/big.S_over_N.iloc[-1]:.2f}x")

    x = 1.0 / big.L.values
    y = big.S_over_N.values
    A = np.column_stack([np.ones(len(x)), x])
    c, m = np.linalg.lstsq(A, y, rcond=None)[0]
    res = y - (c + m*x)
    dof = max(len(x) - 2, 1)
    sig = np.sqrt((res**2).sum()/dof) * np.sqrt(np.linalg.inv(A.T@A)[0,0])
    print(f"\n  linear fit of S(Q)/N against 1/L, using L >= 12:")
    print(f"    intercept = {c:+.6f} +/- {sig:.6f}   ({abs(c)/sig if sig else 0:.1f} sigma from zero)")
    print(f"    -> {'consistent with ZERO: no long-range order' if abs(c) < 2*sig else 'nonzero, would indicate order'}")
    print()


# ---------------------------------------------------------------------------
# Independent check at quarter filling, where L=10 and L=14 both exist.
#
# If S(Q)/N falls as 1/N the ratio between the two sizes is N_10/N_14 = 0.510.
# If the order parameter is finite the ratio sits near 1.
# ---------------------------------------------------------------------------
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "sdwz_pairs")
if os.path.isdir(P):
    print("=" * 74)
    print("n = 0.5, U = 4   L=10 vs L=14   (pure 1/N would give ratio 0.510)")
    print("=" * 74)
    print(f"{'delta':>7}{'S/N (L=10)':>13}{'S/N (L=14)':>13}{'ratio':>9}   verdict")
    for dl in ("0.00", "0.10", "0.20", "0.30", "0.40"):
        try:
            a = analyse(os.path.join(P, f"n0.500_U4.0_d{dl}_L10.dat"), 10)
            b = analyse(os.path.join(P, f"n0.500_U4.0_d{dl}_L14.dat"), 14)
        except OSError:
            continue
        r = b["S_over_N"] / a["S_over_N"]
        v = "1/N decay" if r < 0.62 else ("plateau" if r > 0.85 else "intermediate")
        print(f"{float(dl):>7.2f}{a['S_over_N']:>13.6f}{b['S_over_N']:>13.6f}{r:>9.3f}   {v}")
