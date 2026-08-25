"""Finite-size behaviour of the pairing vertex at U = 4.

The paper's pairing claim rested on L = 12 alone. These runs add L = 8 and L = 10
at settings byte-matched to the committed L = 12 run (beta = 32, N_w = 500,
tau window 0.8, six seeds), so the three sizes form one series.

Two questions. Does d_xy remain the leading channel at every size, and does its
vertex grow with L? Growth with size is what distinguishes a genuine pairing
tendency from a finite-size artifact, and it is the opposite of what the magnetic
order parameter does.

Note the multi-size data already on disk from earlier work is FULL chi, which is
dominated by the uncorrelated bubble and cannot answer either question.

  python analyze_chi_fss.py
"""
import glob, os
import numpy as np
import pandas as pd

CH = {"son": "on-site s", "sext": "ext s", "d": "dx2-y2", "dxy": "dxy"}
frames = {8: "../data/chi_fss_L8_U4.csv", 10: "../data/chi_fss_L10_U4.csv",
          12: "../data/chi_L12/eqtime_L12_U4.csv"}
d = {}
for L, f in frames.items():
    if os.path.exists(f):
        d[L] = pd.read_csv(f)
Ls = sorted(d)
print(f"U = 4, sizes {Ls}, "
      f"beta={sorted(set(np.concatenate([d[L].beta_proj.unique() for L in Ls])))}, "
      f"nw={sorted(set(np.concatenate([d[L].nw.unique() for L in Ls])))}\n")

print("chi_dxy VERTEX (mean +/- standard error over seeds)")
print(f"{'delta':>7} | " + "".join(f"{'L=%d'%L:>18}" for L in Ls) + "   L8 -> L12")
for dl in sorted(d[Ls[0]].delta.unique()):
    row, vals = f"{dl:7.1f} | ", []
    for L in Ls:
        s = d[L][np.isclose(d[L].delta, dl)]["chi_dxy_vertex"]
        if len(s) == 0:
            row += f"{'--':>18}"; vals.append(np.nan); continue
        m, e = s.mean(), s.std(ddof=1)/np.sqrt(len(s))
        row += f"{m:11.4f}+/-{e:5.4f}"; vals.append(m)
    ok = [v for v in vals if np.isfinite(v)]
    if len(ok) >= 2:
        row += f"   {ok[-1]/ok[0]:5.2f}x " + ("grows" if ok[-1] > ok[0] else "falls")
    print(row)

print("\nleading channel at each (L, delta)")
for L in Ls:
    out = []
    for dl in sorted(d[L].delta.unique()):
        s = d[L][np.isclose(d[L].delta, dl)]
        v = {c: s[f"chi_{c}_vertex"].mean() for c in CH}
        out.append(f"{dl:.1f}:{CH[max(v, key=v.get)]}")
    print(f"  L={L:3d}: " + "  ".join(out))

print("\nall four channels at the d_xy peak (delta = 0.4)")
print(f"{'L':>4} | " + "".join(f"{v:>12}" for v in CH.values()))
for L in Ls:
    s = d[L][np.isclose(d[L].delta, 0.4)]
    if len(s) == 0: continue
    print(f"{L:4d} | " + "".join(f"{s[f'chi_{c}_vertex'].mean():12.4f}" for c in CH))
