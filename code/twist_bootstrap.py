"""Bootstrap the delta at which the d_xy vertex changes sign.

Twist average uses the four boundary-condition corners with weights 1:2:1 --
(+1,+1) periodic, (-1,+1) and (+1,-1) which are degenerate by the lattice x<->y
symmetry, and (-1,-1). Seeds are resampled with replacement INSIDE each leg,
so the error bar reflects Monte Carlo noise, not scatter between legs.
"""
import numpy as np, pandas as pd

NB = 4000
rng = np.random.default_rng(12345)
DEL = np.arange(0, 0.75, 0.1)

t = pd.read_csv('twist_merged.csv')
p = pd.read_csv('/Users/liujiaxin/Desktop/checkerboard-altermagnet/data/pairing_master.csv')
p = p[(p.L == 12) & np.isclose(p.n, 1.0)]

# leg -> (weight, {(U,delta): array of per-seed chi_dxy})
legs = {}
for (ax, ay), g in t.groupby(['apx', 'apy']):
    w = 2.0 if ax != ay else 1.0
    legs[(ax, ay)] = (w, {k: v.chi_dxy_vertex.values for k, v in g.groupby(['U', 'delta'])})
legs[(1, 1)] = (1.0, {k: v.chi_dxy.values for k, v in p.groupby(['U', 'delta'])})

US = [2., 4., 6., 8.]


def crossing(x, y):
    """First sign change of y(x), linearly interpolated. nan if none."""
    i = np.where(np.diff(np.sign(y)))[0]
    if not len(i):
        return np.nan
    k = i[0]
    return x[k] + (-y[k] / (y[k + 1] - y[k])) * (x[k + 1] - x[k])


def curve(U, resample):
    num = np.zeros(len(DEL)); den = 0.0
    for w, tab in legs.values():
        vals = []
        for d in DEL:
            a = tab[(U, round(d, 1))]
            vals.append(rng.choice(a, len(a), replace=True).mean() if resample else a.mean())
        num += w * np.array(vals); den += w
    return num / den


print('d_xy sign change in delta -- L=12, half filling, 4000 bootstrap resamples')
print(f"{'U':>3} {'twist-avg dc':>16} {'68% CI':>16} {'no-crossing':>12}   {'periodic dc':>12}")
for U in US:
    dc = crossing(DEL, curve(U, False))
    bs = np.array([crossing(DEL, curve(U, True)) for _ in range(NB)])
    ok = bs[~np.isnan(bs)]
    lo, hi = np.percentile(ok, [16, 84])
    pv = p[p.U == U].groupby('delta').chi_dxy.mean()
    dp = crossing(pv.index.values, pv.values)
    print(f'{U:3g} {dc:16.3f} {f"[{lo:.3f}, {hi:.3f}]":>16} {1-len(ok)/NB:11.1%}   {dp:12.3f}')

print()
print('Is dc the same at every U?  spread of the twist-averaged values:')
v = [crossing(DEL, curve(U, False)) for U in US]
print(f'  min {np.nanmin(v):.3f}  max {np.nanmax(v):.3f}  range {np.nanmax(v)-np.nanmin(v):.3f}')
