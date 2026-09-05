"""Staggered moment of the UHF TRIAL itself, so the QMC Delta_tot can be read
against what the trial handed it. A near-zero QMC value under a strongly
magnetic trial means the constrained-path projection suppressed the splitting
(physics); under a weak trial it would mean nothing (artifact)."""
import numpy as np, glob, os, re
for d in sorted(glob.glob(os.path.expanduser("~/Checkerboard_Model/L12n1.000*/"))):
    nm = os.path.basename(d.rstrip("/"))
    try:
        if "time_sec" not in "".join(open(d+"nohup.out").readlines()[-3:]): continue
    except Exception: continue
    try:
        P = {s: np.loadtxt(d+f"wf{s}.txt") for s in ("up","dn")}
    except Exception: continue
    n = {}
    for s, M in P.items():
        M = M.reshape(144, -1)
        G = M @ np.linalg.solve(M.T @ M, M.T)
        n[s] = np.diag(G)
    i = np.arange(144); sg = (-1.0)**((i//12)+(i%12))
    m = float(np.abs((n["up"]-n["dn"])*sg).sum()/144)
    U = re.search(r"u([0-9.]+)tA", nm).group(1)
    dl = re.search(r"tA-([0-9.]+)tt", nm).group(1)
    print(f"{U},{dl},{m:.4f}")
