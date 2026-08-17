"""Pairing channels for EVERY Fortran run, on the same (U, delta) / (n, delta)
axes as Delta_tot -- the two live in the same run folders.

Channel map from mc2duph.f90 (lines 1353-1355, 1386):
    sowave -> on-site s   swave -> extended s
    dwave  -> d_x2-y2     dd12wave -> d_xy
Vertex_*.dat is the connected part, at k = (0,0).
"""
import numpy as np, pandas as pd, glob, os, re

MAP = {"sowave": "son", "swave": "sext", "dwave": "d", "dd12wave": "dxy"}
rows = []
for p in sorted(glob.glob("/home/phd25imran/Checkerboard_Model/L*/dir-kVals")):
    nm = os.path.basename(os.path.dirname(p))
    m = re.match(r"L(\d+)n([0-9.]+)u([0-9.]+)tA-?([0-9.]+)tt", nm)
    if not m:
        continue
    r = dict(L=int(m.group(1)), n=float(m.group(2)),
             U=float(m.group(3)), delta=float(m.group(4)))
    ok = True
    for f, key in MAP.items():
        try:                              # line 2 of each file is k = (0,0)
            with open(f"{p}/Vertex_{f}.dat") as fh:
                fh.readline()
                r[key] = float(fh.readline().split()[2])
        except Exception:
            ok = False; break
    if ok:
        rows.append(r)

d = pd.DataFrame(rows).sort_values(["L", "U", "delta", "n"])
d.to_csv("/home/phd25imran/analysis/data/fortran_pairing_all.csv", index=False)
print(f"wrote {len(d)} runs")
for L in sorted(d.L.unique()):
    s = d[d.L == L]
    print(f"  L={L}: {len(s)} runs, U={sorted(s.U.unique())}")
h = d[(d.L == 14) & np.isclose(d.n, 1.0)]
print("\nleading channel, L=14 half filling (rows U, cols delta):")
K = ["son", "sext", "d", "dxy"]
w = pd.Series([K[i] for i in h[K].values.argmax(1)],
              index=pd.MultiIndex.from_arrays([h.U, h.delta]))
print(w.unstack().iloc[::-1].to_string())
