"""Re-extract the L=12 half-filling grid KEEPING the k-resolved arrays.

extract49.py computed n_up(k) and n_dn(k) and then threw them away, saving only
the scalar Delta_tot. This keeps both: the scalar AND the full Delta n(k) map
per run, so the Brillouin-zone panels can be drawn at delta = 0.1, 0.4, 0.7.

Writes two files to ~/analysis:
    fortran_L12_all49.csv       scalars (same as before, plus pairing channels)
    fortran_L12_maps.npz        maps[i] = (nk, 3) array of kx, ky, dn ; meta rows L,n,U,delta
"""
import numpy as np, pandas as pd, glob, os, re
MAP = {"sowave": "son", "swave": "sext", "dwave": "d", "dd12wave": "dxy"}
rows, maps, meta = [], [], []
for d in sorted(glob.glob(os.path.expanduser("~/Checkerboard_Model/L12n1.000*/"))):
    nm = os.path.basename(d.rstrip("/"))
    m = re.match(r"L(\d+)n([0-9.]+)u([0-9.]+)tA-([0-9.]+)tt", nm)
    if not m: continue
    try:
        if "time_sec" not in "".join(open(d + "nohup.out").readlines()[-3:]): continue
    except OSError: continue
    k = d + "dir-kVals/"
    try:
        up = pd.read_csv(k+"n_up.dat", sep=r"\s+", skiprows=1, header=None,
                         names=["kx","ky","v","e"])
        dn = pd.read_csv(k+"n_dn.dat", sep=r"\s+", skiprows=1, header=None,
                         names=["kx","ky","v","e"])
    except OSError: continue
    L=int(m.group(1)); n_=float(m.group(2)); U=float(m.group(3)); dl=float(m.group(4))
    diff = up.v.values - dn.v.values
    row=[L,n_,U,dl,float(np.abs(diff).sum())/L**2,float(diff.sum())]
    for f in MAP:
        try:
            with open(f"{k}Vertex_{f}.dat") as fh:
                fh.readline(); row.append(float(fh.readline().split()[2]))
        except (OSError, IndexError, ValueError): row.append(np.nan)
    rows.append(row)
    # KEEP the momentum-resolved difference, not just its sum
    maps.append(np.column_stack([up.kx.values, up.ky.values, diff]))
    meta.append([L, n_, U, dl])

cols="L,n,U,delta,dtot_N,net,son,sext,d,dxy"
out=os.path.expanduser("~/analysis")
with open(f"{out}/fortran_L12_all49.csv","w") as fh:
    fh.write(cols+"\n")
    for r in rows: fh.write(",".join(str(x) for x in r)+"\n")
np.savez_compressed(f"{out}/fortran_L12_maps.npz",
                    maps=np.array(maps), meta=np.array(meta),
                    cols=np.array(["L","n","U","delta"]))
print(f"{len(rows)} runs; maps array {np.array(maps).shape}")
