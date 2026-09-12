"""Our pairing channels against those of Pan et al. (arXiv:2409.16523).

Evaluated at delta = 0.3, half filling, the single point where the two Hamiltonians
coincide (their t'/t = 0.6; identical spectra to 4e-15; the sign of t' is fixed by
particle-hole symmetry, which is exact at half filling).

Their basis has no d_xy. The comparison therefore reports the leading channel three
ways: in the full basis, in ours with d_xy removed, and in theirs alone.

CAVEAT that must travel with these numbers: they work at n = 0.9 and T/t = 1/6,
we at half filling in the ground state. Under particle-hole their n = 0.9 maps to
n = 1.1, so their doped results are not ours even in principle.
"""
import glob, os, sys
import numpy as np, pandas as pd
src = sys.argv[1] if len(sys.argv) > 1 else "../data/pan"
d = pd.concat([pd.read_csv(x) for x in sorted(glob.glob(os.path.join(src, "pan_*.csv")))],
              ignore_index=True)
CH = {"son":"on-site s","sext":"ext s","d":"dx2-y2","dxy":"dxy","pan_s":"Pan s","pan_d":"Pan d"}
print(f"{'U':>5} | " + "".join(f"{v:>12}" for v in CH.values()))
for U in sorted(d.U.unique()):
    s = d[d.U == U]
    print(f"{U:5.1f} | " + "".join(f"{s[f'chi_{c}_vertex'].mean():12.4f}" for c in CH))
print()
for U in sorted(d.U.unique()):
    s = d[d.U == U]
    o = {c: s[f"chi_{c}_vertex"].mean() for c in ("son","sext","d","dxy")}
    p = {c: s[f"chi_{c}_vertex"].mean() for c in ("pan_s","pan_d")}
    n = {c: v for c, v in o.items() if c != "dxy"}
    print(f"U={U:4.1f}: full -> {CH[max(o,key=o.get)]:9s} | minus d_xy -> {CH[max(n,key=n.get)]:9s}"
          f" | their basis -> {CH[max(p,key=p.get)]}")
