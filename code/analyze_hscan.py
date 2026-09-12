"""h -> 0 extrapolation of Delta_tot: is the spin splitting spontaneous?

Delta_tot is measured with a staggered pinning field h, which makes the trial
state symmetry broken. Because the constrained path is defined by the trial, that
biases the walkers toward the broken state, so a finite Delta_tot at finite h
proves nothing on its own. The physical question is the h -> 0 intercept:

    intercept > 0, holding or growing with L   ->  spontaneous order
    intercept consistent with 0                ->  the splitting is field-induced

This is the test flagged as mandatory in checkerboard_polarization.py and never
run until now. It should be read against the finite-size scaling of S(pi,pi),
which falls as 1/N at every delta, filling and U up to 8 (no long-range order).

Both a linear and a quadratic fit in h are reported. The linear fit uses the three
smallest fields, where the response should be closest to linear; disagreement
between the two fits is itself information about how far the extrapolation reaches.

  python analyze_hscan.py [dir]
"""
import sys, glob, os
import numpy as np
import pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "../data/hscan"
files = sorted(glob.glob(os.path.join(src, "polarization_L*.csv")))
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
d["dtot_N"] = d.delta_tot / (d.L ** 2)
print(f"U={sorted(d.U.unique())} delta={sorted(d.delta.unique())} "
      f"h={sorted(d.h.unique())} L={sorted(d.L.unique())} "
      f"seeds={len(d.seed.unique())}\n")

for L in sorted(d.L.unique()):
    s = d[d.L == L]
    g = s.groupby("h")["dtot_N"]
    h = np.array(sorted(s.h.unique()))
    y = g.mean().reindex(h).to_numpy()
    e = (g.std(ddof=1) / np.sqrt(g.count())).reindex(h).to_numpy()
    print(f"L = {L}  ({L*L} sites)")
    for hh, yy, ee in zip(h, y, e):
        print(f"    h = {hh:5.2f}   Delta_tot/N = {yy:.6f} +/- {ee:.6f}")
    # linear fit on the three smallest fields, weighted
    k = 3
    w = 1.0 / np.maximum(e[:k], 1e-12) ** 2
    c1 = np.polyfit(h[:k], y[:k], 1, w=np.sqrt(w))
    # quadratic on all five
    c2 = np.polyfit(h, y, 2, w=np.sqrt(1.0 / np.maximum(e, 1e-12) ** 2))
    # intercept uncertainty by reseeding the fit within the error bars
    rng = np.random.default_rng(0)
    b = [np.polyfit(h[:k], rng.normal(y[:k], e[:k]), 1)[-1] for _ in range(2000)]
    print(f"    linear  (h <= {h[k-1]:.2f}) intercept = {c1[-1]:+.6f} +/- {np.std(b):.6f}")
    print(f"    quadratic (all h)   intercept = {c2[-1]:+.6f}")
    print(f"    -> {'NONZERO' if c1[-1] > 3*np.std(b) else 'consistent with ZERO'}"
          f"  (floor from the delta=0 control: 0.000626 +/- 0.000214)\n")

print("Summary: intercept vs L (growing or steady => order survives the limit)")
for L in sorted(d.L.unique()):
    s = d[d.L == L]; g = s.groupby("h")["dtot_N"]
    h = np.array(sorted(s.h.unique())); y = g.mean().reindex(h).to_numpy()
    e = (g.std(ddof=1)/np.sqrt(g.count())).reindex(h).to_numpy()
    rng = np.random.default_rng(0)
    b = [np.polyfit(h[:3], rng.normal(y[:3], e[:3]), 1)[-1] for _ in range(2000)]
    print(f"  L={L:3d}: {np.polyfit(h[:3], y[:3], 1)[-1]:+.6f} +/- {np.std(b):.6f}")
