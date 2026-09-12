# Unbiased control trials for the Lieb altermagnet, to test whether the order
# seen at t2 = 0 is physics or constrained-path trial bias.
#
# The production runs use an ALTERMAGNETIC trial (A neutral, B up, C down). CPQMC
# is constrained path, so that trial can hold the walkers in the altermagnetic
# state. The 48-site scan gives its LARGEST order parameter at t2 = 0
# (m_B - m_C = 1.96), while the published mean-field plus small-cluster ED work
# says t2 is required to stabilise the state. Before claiming a disagreement we
# have to start from a trial that carries no magnetic order at all.
#
#   MODE=free   phiT = lowest NUP eigenvectors of the one-body K. No interaction,
#               no spin dependence, identical for both spins. At n = 4 per cell
#               and t2 = 0 this is closed shell (lower dispersive band plus the
#               full flat band), so there is no degeneracy ambiguity.
#   MODE=pm     the self-consistent solution reached from the PARAMAGNETIC seed,
#               taken regardless of its energy, i.e. the SCF is not allowed to
#               fall into the altermagnetic basin.
#
# Writes wfup.txt / wfdn.txt exactly like wf_lieb.py.
import os, re, sys, numpy as np, importlib.util

root = os.getcwd()
mode = os.environ.get("MODE", "free")

src = os.path.join(root, "wf_lieb.py")
body = open(src).read().split("# ================================================================\n# Main")[0]
open("/tmp/_wfl.py", "w").write("EPS_A = 0.0\n" + body)
spec = importlib.util.spec_from_file_location("_wfl", "/tmp/_wfl.py")
w = importlib.util.module_from_spec(spec); spec.loader.exec_module(w)

folder = os.path.basename(root)
L  = int(re.search(r"L(\d+)", folder).group(1))
u  = float(re.search(r"u([\d.]+)", folder).group(1))
tA = float(re.search(r"tA(-?[\d.]+)", folder).group(1))
tt = float(re.search(r"tt(-?[\d.]+)", folder).group(1))
p  = open(os.path.join(root, "parameter.f90")).read()
m  = re.search(r"parameter *\( *NUP *= *(\d+) *, *NDN *= *(\d+)", p)
NUP, NDN = int(m.group(1)), int(m.group(2))
EPS = float(open(os.path.join(root, "in.dat")).readlines()[4].split(",")[2])
w.EPS_A = EPS

print(f"CONTROL TRIAL  mode={mode}  L={L} U={u} t2={tA} NUP={NUP} NDN={NDN} eps_A={EPS}")

K = w.GetK(L, L, 1, tt, tA)
if EPS != 0.0:
    for i, (ix, iy) in enumerate(w.lieb_sites(L, L)):
        if ix % 2 == 1 and iy % 2 == 1:
            K[i, i] += EPS

if mode == "free":
    ev, vec = np.linalg.eigh(K)
    gap_u = ev[NUP] - ev[NUP - 1] if NUP < len(ev) else float("nan")
    print(f"  one-body spectrum: HOMO={ev[NUP-1]:.6f} LUMO={ev[NUP]:.6f} gap={gap_u:.6f}")
    if abs(gap_u) < 1e-8:
        print("  WARNING: open shell at this filling, the free trial is ambiguous")
    phiT_up = vec[:, :NUP].T.copy()
    phiT_dn = vec[:, :NDN].T.copy()
elif mode == "pm":
    nup, ndn, E, vu, vd = w.Iteration(L, L, u, tA, tt, NUP=NUP, NDN=NDN, seed="pm")
    print(f"  paramagnetic SCF energy (variational): {E:.6f}")
    phiT_up = vu[:, :NUP].T.copy()
    phiT_dn = vd[:, :NDN].T.copy()
else:
    sys.exit(f"unknown MODE={mode}")

# report the magnetic content of the trial itself, so the control is auditable
sites = w.lieb_sites(L, L)
nu = (phiT_up.T ** 2).sum(1); nd = (phiT_dn.T ** 2).sum(1)
mm = nu - nd
lab = ["A" if (ix % 2 and iy % 2) else ("B" if iy % 2 else "C") for ix, iy in sites]
g = {s: float(np.mean([mm[i] for i, t in enumerate(lab) if t == s])) for s in "ABC"}
print(f"  trial moments: m_A={g['A']:+.5f} m_B={g['B']:+.5f} m_C={g['C']:+.5f} "
      f"m_B-m_C={g['B']-g['C']:+.5f}   <-- must be ~0 for an unbiased control")

phiT_up[np.abs(phiT_up) < 1e-10] = 0.0
phiT_dn[np.abs(phiT_dn) < 1e-10] = 0.0
np.savetxt(os.path.join(root, "wfup.txt"), phiT_up, fmt="%.12e", delimiter=" ")
np.savetxt(os.path.join(root, "wfdn.txt"), phiT_dn, fmt="%.12e", delimiter=" ")
print("Saved: wfup.txt  wfdn.txt")
