# Full Lieb altermagnet summary across the finished 48-site runs.
#
# Reports the sublattice-resolved moments and the altermagnetic order parameter
# m_B - m_C. B and C are related by the C4 rotation about an A site and by no
# translation or inversion, so antialigned B and C is altermagnetic rather than
# antiferromagnetic. A carries no interaction (U lives only on B and C) and
# should stay unpolarised.
import sys, os, glob, re, numpy as np

root = sys.argv[1]
rows = []
for d in sorted(glob.glob(os.path.join(root, "L8*/"))):
    name = os.path.basename(d.rstrip("/"))
    m = re.search(r"u([\d.]+)tA([\d.]+)", name)
    U, t2 = float(m.group(1)), float(m.group(2))
    L = int(re.search(r"L(\d+)", name).group(1)); cells = L // 2; NS = 3 * cells * cells
    sites = [(ix, iy) for ix in range(1, L+1) for iy in range(1, L+1)
             if not (ix % 2 == 0 and iy % 2 == 0)]
    def diag(f):
        p = os.path.join(d, "dir-rVals", f)
        if not os.path.exists(p): return None, None
        a = np.loadtxt(p, skiprows=1)
        return np.diag(a[:, 2].reshape(NS, NS)), np.diag(a[:, 3].reshape(NS, NS))
    nu, eu = diag("gx_up.dat"); nd, ed = diag("gx_dn.dat")
    if nu is None: continue
    mm = nu - nd; err = np.hypot(eu, ed)
    lab = ["A" if (ix % 2 and iy % 2) else ("B" if iy % 2 else "C") for ix, iy in sites]
    g = {s: np.array([mm[i] for i, t in enumerate(lab) if t == s]) for s in "ABC"}
    ge = {s: np.array([err[i] for i, t in enumerate(lab) if t == s]) for s in "ABC"}
    e_bc = float(np.sqrt((ge["B"]**2).sum() + (ge["C"]**2).sum()) / len(g["B"]))
    rows.append((U, t2, g["A"].mean(), g["B"].mean(), g["C"].mean(),
                 g["B"].mean() - g["C"].mean(), e_bc, mm.mean()))

print(f"{'U':>6}{'t2':>6}{'m_A':>10}{'m_B':>10}{'m_C':>10}{'m_B-m_C':>11}{'err':>9}{'net/site':>10}")
for r in sorted(rows, key=lambda r: (r[1], r[0])):
    print(f"{r[0]:>6.1f}{r[1]:>6.1f}{r[2]:>10.4f}{r[3]:>10.4f}{r[4]:>10.4f}"
          f"{r[5]:>11.4f}{r[6]:>9.4f}{r[7]:>10.4f}")
print(f"\n  runs: {len(rows)}   sites: {NS}   n = 4 electrons per unit cell (compensated, NUP = NDN)")
mA = max(abs(r[2]) for r in rows)
print(f"  largest |m_A| anywhere: {mA:.5f}  (A should be non-magnetic: U acts only on B and C)")
