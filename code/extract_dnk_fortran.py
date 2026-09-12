"""Delta_tot and the Delta n(k) maps, straight from the Fortran output.

Every run directory holds dir-kVals/n_up.dat and n_dn.dat, the momentum-resolved
occupations. The altermagnetic order parameter of the PRL is

    Delta_tot = sum_k |n_up(k) - n_dn(k)|

The Fortran uses a symmetry-broken trial (wfup.txt / wfdn.txt are separate), so
the spin degeneracy is already lifted and the difference is physical rather than
zero by symmetry as it is in the spin-restricted Python runs.
"""
import numpy as np, pandas as pd, glob, os, re

rows, maps, meta = [], [], []
for p in sorted(glob.glob("/home/phd25imran/Checkerboard_Model/L*/dir-kVals")):
    nm = os.path.basename(os.path.dirname(p))
    m = re.match(r"L(\d+)n([0-9.]+)u([0-9.]+)tA-?([0-9.]+)tt", nm)
    if not m:
        continue
    L, n, U, dl = int(m.group(1)), float(m.group(2)), float(m.group(3)), float(m.group(4))
    try:
        u = pd.read_csv(p + "/n_up.dat", sep=r"\s+", skiprows=1, header=None,
                        names=["kx", "ky", "v", "e"])
        w = pd.read_csv(p + "/n_dn.dat", sep=r"\s+", skiprows=1, header=None,
                        names=["kx", "ky", "v", "e"])
    except Exception:
        continue
    dn = u.v.values - w.v.values
    err = np.sqrt(u.e.values ** 2 + w.e.values ** 2)
    kx, ky = u.kx.values, u.ky.values
    f = np.sin(kx) * np.sin(ky)
    rows.append(dict(L=L, n=n, U=U, delta=dl,
                     delta_tot=float(np.abs(dn).sum()),
                     delta_tot_err=float(np.sqrt((err ** 2).sum())),
                     net=float(dn.sum()),                 # must be ~0: compensated
                     snr=float(np.abs(dn / err).mean()),
                     corr_dxy=float(np.corrcoef(dn, f)[0, 1])))
    maps.append(np.column_stack([kx, ky, dn, err]))
    meta.append([L, n, U, dl])

d = pd.DataFrame(rows).sort_values(["L", "U", "delta", "n"])
d.to_csv("/home/phd25imran/analysis/data/fortran_dnk.csv", index=False)
np.savez_compressed("/home/phd25imran/analysis/data/fortran_dnk_maps.npz",
                    maps=np.array(maps, dtype=object), meta=np.array(meta),
                    cols=np.array(["L", "n", "U", "delta"]))
print(f"wrote {len(d)} runs")
print("\nnet magnetisation sum_k dn(k) -- should be ~0 if compensated:")
print(f"  |net| max = {d.net.abs().max():.3e}   median = {d.net.abs().median():.3e}")
print("\nL=14 half filling, the (U, delta) grid:")
h = d[(d.L == 14) & (np.isclose(d.n, 1.0))]
print(h.pivot_table(index="U", columns="delta", values="delta_tot").round(3).to_string())
