# k-space sum rules and the magnetic peak for a finished Lieb run.
#
# The writer duplicates zone-boundary rows so that a plot closes on itself, and
# it writes each duplicate immediately after its parent, so the mesh is NOT the
# first NKPTS rows. Every duplicate carries kx = -pi or ky = -pi, which the
# folded mesh 2*pi*i/L - 2*pi*[i > L/2] never produces, so dropping those rows
# leaves exactly the mesh. Sum rules: FourierTransform divides by
# the number of SITES while summing over the fine grid's NKPTS points, and
# sum_k cos(k.(r_m - r_n)) = NKPTS * delta_mn, so
#     sum_k n_sigma(k) = (NKPTS / NSTATES) * N_sigma
#
#   python check_lieb_kspace.py <run directory>
import sys, os, numpy as np

d = sys.argv[1] if len(sys.argv) > 1 else "."
cells = None
for tok in os.path.basename(os.path.abspath(d)).split("n"):
    if tok.startswith("L") or tok[:1].isdigit():
        pass
L = int(os.path.basename(os.path.abspath(d)).split("L")[1].split("n")[0])
cells = L // 2
NK, NS = L * L, 3 * cells * cells
NUP, NDN = 2 * cells * cells, cells * cells

def load(name):
    a = np.loadtxt(os.path.join(d, "dir-kVals", name), skiprows=1)
    keep = (a[:, 0] > -np.pi + 1e-6) & (a[:, 1] > -np.pi + 1e-6)
    a = a[keep]
    assert len(a) == NK, f"{name}: kept {len(a)} rows, expected {NK}"
    return a

u, dn, s = load("n_up.dat"), load("n_dn.dat"), load("sdwz.dat")
print(f"{os.path.basename(os.path.abspath(d))}:  sites={NS}  k points={NK}  "
      f"NUP={NUP} NDN={NDN}")
print(f"  sum n_up(k) = {u[:,2].sum():.6f}   expected {NK/NS*NUP:.6f}")
print(f"  sum n_dn(k) = {dn[:,2].sum():.6f}   expected {NK/NS*NDN:.6f}")
dt = np.abs(u[:, 2] - dn[:, 2]).sum()
print(f"  Delta_tot   = {dt:.6f}   Delta_tot/N = {dt/NS:.6f}")
i = int(np.argmax(s[:, 2]))
print(f"  sdwz peaks at k = ({s[i,0]:+.4f}, {s[i,1]:+.4f}) with {s[i,2]:.6f}")
g = s[(np.abs(s[:, 0]) < 1e-9) & (np.abs(s[:, 1]) < 1e-9), 2]
print(f"  sdwz at Gamma = {g[0]:.6f}" if len(g) else "  no Gamma row")
