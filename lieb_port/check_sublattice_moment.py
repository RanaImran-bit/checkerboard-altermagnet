# Sublattice-resolved moments on the Lieb lattice.
#
# The diagonal of gx_up / gx_dn in dir-rVals is the per-site density, so
# m_i = <n_i^up> - <n_i^dn> resolves A, B and C separately. This is the
# measurement that decides whether the Lieb lattice can carry altermagnetism at
# all: B and C are related by the C4 rotation, not by a translation, so opposite
# moments on B and C would be altermagnetic. Lieb's theorem puts B and C in the
# same bipartite sublattice at half filling and so forbids it there, which makes
# this a null test at n=1 and the real question away from half filling.
#
# Site order matches mc2duph.f90: ix outermost, iy innermost, (even,even) skipped
# in 1-based fine-grid coordinates.
#
#   python check_sublattice_moment.py <run directory>
import sys, os, numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else "."
name = os.path.basename(os.path.abspath(d))
L = int(name.split("L")[1].split("n")[0])
cells = L // 2
NS = 3 * cells * cells

sites = [(ix, iy) for ix in range(1, L + 1) for iy in range(1, L + 1)
         if not (ix % 2 == 0 and iy % 2 == 0)]
assert len(sites) == NS

def diag(fname):
    a = np.loadtxt(os.path.join(d, "dir-rVals", fname), skiprows=1)
    assert len(a) == NS * NS, f"{fname}: {len(a)} rows, expected {NS*NS}"
    v = a[:, 2].reshape(NS, NS)
    e = a[:, 3].reshape(NS, NS)
    return np.diag(v), np.diag(e)

nu, eu = diag("gx_up.dat")
nd, ed = diag("gx_dn.dat")
m, em = nu - nd, np.hypot(eu, ed)

lab = ["A" if (ix % 2 and iy % 2) else ("B" if iy % 2 else "C")
       for ix, iy in sites]
print(f"{name}:  {NS} sites")
print(f"{'sub':>5}{'count':>7}{'<n_up>':>12}{'<n_dn>':>12}{'m':>12}{'err':>12}")
for s in "ABC":
    k = [i for i, t in enumerate(lab) if t == s]
    print(f"{s:>5}{len(k):>7}{nu[k].mean():>12.6f}{nd[k].mean():>12.6f}"
          f"{m[k].mean():>12.6f}{np.sqrt((em[k]**2).sum())/len(k):>12.6f}")
mA = m[[i for i, t in enumerate(lab) if t == "A"]].mean()
mB = m[[i for i, t in enumerate(lab) if t == "B"]].mean()
mC = m[[i for i, t in enumerate(lab) if t == "C"]].mean()
print(f"\n  net moment per site = {m.mean():+.6f}")
print(f"  m_B - m_C           = {mB - mC:+.6f}   "
      f"(nonzero would break the C4 that relates B and C)")
print(f"  m_A vs m_B + m_C    = {mA:+.6f} vs {mB + mC:+.6f}")
