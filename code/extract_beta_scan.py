import numpy as np, pandas as pd, glob, os, re
snap=lambda a: np.where(np.abs(np.asarray(a))>np.pi-1e-4,-np.pi,np.asarray(a))
for d in sorted(glob.glob(os.path.expanduser("~/Checkerboard_Model/beta*/"))):
    nm=os.path.basename(d.rstrip("/"))
    m=re.match(r"beta(\d+)_L(\d+)n[0-9.]+u([0-9.]+)tA-([0-9.]+)tt",nm)
    if not m: continue
    try:
        if "time_sec" not in "".join(open(d+"nohup.out").readlines()[-3:]): continue
    except OSError: continue
    k=d+"dir-kVals/"
    try:
        up=pd.read_csv(k+"n_up.dat",sep=r"\s+",skiprows=1,header=None,names=["kx","ky","v","e"])
        dn=pd.read_csv(k+"n_dn.dat",sep=r"\s+",skiprows=1,header=None,names=["kx","ky","v","e"])
    except OSError: continue
    L=int(m.group(2)); diff=up.v.values-dn.v.values
    agg={}
    for a,b,v in zip(snap(up.kx.values),snap(up.ky.values),diff):
        agg.setdefault((round(a,4),round(b,4)),[]).append(v)
    vals=np.array([np.mean(v) for v in agg.values()])
    print(f"{m.group(1)},{L},{m.group(3)},{m.group(4)},{np.abs(vals).sum()/L**2:.6f},{len(vals)}")
