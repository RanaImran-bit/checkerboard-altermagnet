import numpy as np, pandas as pd, glob, os, re
MAP={"sowave":"son","swave":"sext","dwave":"d","dd12wave":"dxy"}
for d in sorted(glob.glob(os.path.expanduser("~/Checkerboard_Model/L12n0.847*/"))):
    nm=os.path.basename(d.rstrip("/"))
    m=re.match(r"L(\d+)n([0-9.]+)u([0-9.]+)tA-([0-9.]+)tt",nm)
    if not m: continue
    try:
        if "time_sec" not in "".join(open(d+"nohup.out").readlines()[-3:]): continue
    except OSError: continue
    k=d+"dir-kVals/"
    try:
        up=pd.read_csv(k+"n_up.dat",sep=r"\s+",skiprows=1,header=None,names=["kx","ky","v","e"])
        dn=pd.read_csv(k+"n_dn.dat",sep=r"\s+",skiprows=1,header=None,names=["kx","ky","v","e"])
    except OSError: continue
    L=int(m.group(1)); diff=up.v.values-dn.v.values
    row=[L,float(m.group(2)),float(m.group(3)),float(m.group(4)),
         float(np.abs(diff).sum())/L**2,float(diff.sum())]
    for f in MAP:
        try:
            with open(f"{k}Vertex_{f}.dat") as fh:
                fh.readline(); row.append(float(fh.readline().split()[2]))
        except (OSError,IndexError,ValueError): row.append(np.nan)
    print(",".join(str(x) for x in row))
