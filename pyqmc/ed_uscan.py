#!/usr/bin/env python3
"""EXACT (ED) two-orbital altermagnet: scan the on-site Hubbard U at FIXED anisotropic
bands (the leeb2024spontaneous / manuscript Fig.S9 setup: t1,t2,t3,t4 fixed, anisotropy
built into the orbitals, NOT scanned). Interaction is what is varied. We track whether
altermagnetic order and d-wave pairing turn on together as U grows.

Observables (ground state):
  S_d, S_s : d_{x2-y2} / ext-s intra-orbital NN pair structure factors  <O_a^dag O_a>
  S_AM     : altermagnetic order structure factor = < (sum_i g_i S^z_i)^2 >,
             g_i = +1 (orbital 1) / -1 (orbital 2)  -- orbital-staggered magnetization
  S_Neel   : <(sum_i (-1)^{x+y} S^z_i)^2>  (columnar Neel, reference)
  m2       : uniform <(sum_i S^z_i)^2> (net-moment reference; AM has zero net moment)

  source tools/env.sh
  python pyqmc/ed_uscan.py --lx 2 --ly 2 --nup 4 --ndn 4 \
      --t1 1 --t2 1.75 --t3 0.85 --t4 0.65 --Ulist 0 1 2 4 6 8 12
"""
from __future__ import annotations
import os, sys, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")
import numpy as np
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ed"))


def pair_terms(lx, ly, fkind):
    """intra-orbital NN pair operator terms; fkind='d' -> (+x,-x:+1;+y,-y:-1), 's' -> all +1."""
    lxy = lx * ly
    def idx(x, y, orb): return orb * lxy + (x % lx) * ly + (y % ly)
    out = []
    for orb in (0, 1):
        for x in range(lx):
            for y in range(ly):
                m = idx(x, y, orb)
                for (dx, dy, fd) in ((1,0,+1.0),(-1,0,+1.0),(0,1,-1.0),(0,-1,-1.0)):
                    f = fd if fkind == "d" else 1.0
                    out.append((m, idx(x+dx, y+dy, orb), f))
    return out


def scan_point(lx, ly, nup, ndn, t1, t2, t3, t4, U, uxy=0.0, v=0.0):
    from quspin.operators import hamiltonian
    from altermagnet_ed import _build_static
    static, basis, nsites = _build_static(lx, ly, nup, ndn, t1, t2, t3, t4, U, uxy, v)
    nc = dict(check_pcon=False, check_symm=False, check_herm=False)
    H = hamiltonian(static, [], basis=basis, dtype=np.float64, **nc)
    w, V = np.linalg.eigh(H.toarray()); psi0 = V[:, 0]; E0 = float(w[0]); p2 = psi0 ** 2
    lxy = lx * ly
    # ---- pairing structure factors S_a = <O_a^dag O_a> ----
    def pairSF(fkind):
        bonds = pair_terms(lx, ly, fkind)
        terms = [[fa*fb, m, mm, j, k] for (m,j,fa) in bonds for (mm,k,fb) in bonds]
        O = hamiltonian([["+-|+-", terms]], [], basis=basis, dtype=np.float64, **nc)
        return float(psi0 @ O.dot(psi0))
    Sd = pairSF("d"); Ss = pairSF("s")
    # ---- diagonal S^z vectors for spin structure factors ----
    du = [np.real(hamiltonian([["n|", [[1.0, i]]]], [], basis=basis, dtype=np.float64, **nc).diagonal()) for i in range(nsites)]
    dd = [np.real(hamiltonian([["|n", [[1.0, i]]]], [], basis=basis, dtype=np.float64, **nc).diagonal()) for i in range(nsites)]
    sz = [0.5*(du[i]-dd[i]) for i in range(nsites)]
    def xyorb(i):
        orb = i // lxy; loc = i % lxy; return loc // ly, loc % ly, orb
    gAM = np.array([+1.0 if (i//lxy)==0 else -1.0 for i in range(nsites)])      # orbital sign
    gNe = np.array([(-1.0)**(xyorb(i)[0]+xyorb(i)[1]) for i in range(nsites)])  # (pi,pi)
    OAM = sum(gAM[i]*sz[i] for i in range(nsites))
    ONe = sum(gNe[i]*sz[i] for i in range(nsites))
    Ouni = sum(sz[i] for i in range(nsites))
    S_AM = float((p2*OAM*OAM).sum()); S_Ne = float((p2*ONe*ONe).sum()); m2 = float((p2*Ouni*Ouni).sum())
    return dict(E0=E0, Sd=Sd, Ss=Ss, S_AM=S_AM, S_Neel=S_Ne, m2=m2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lx", type=int, default=2); ap.add_argument("--ly", type=int, default=2)
    ap.add_argument("--nup", type=int, default=4); ap.add_argument("--ndn", type=int, default=4)
    ap.add_argument("--t1", type=float, default=1.0); ap.add_argument("--t2", type=float, default=1.75)
    ap.add_argument("--t3", type=float, default=0.85); ap.add_argument("--t4", type=float, default=0.65)
    ap.add_argument("--uxy", type=float, default=0.0); ap.add_argument("--v", type=float, default=0.0)
    ap.add_argument("--Ulist", type=float, nargs="+", default=[0,1,2,3,4,6,8,12])
    a = ap.parse_args()
    print(f"# ED two-orbital U-scan {a.lx}x{a.ly} nup={a.nup} ndn={a.ndn} "
          f"t1={a.t1} t2={a.t2} t3={a.t3} t4={a.t4} (fixed bands; scanning U)")
    print(f"# {'U':>5} {'E0':>10} {'S_d':>9} {'S_s':>9} {'S_AM':>9} {'S_Neel':>9} {'m2':>8}")
    for U in a.Ulist:
        r = scan_point(a.lx, a.ly, a.nup, a.ndn, a.t1, a.t2, a.t3, a.t4, U, a.uxy, a.v)
        print(f"  {U:>5.1f} {r['E0']:>10.4f} {r['Sd']:>9.4f} {r['Ss']:>9.4f} "
              f"{r['S_AM']:>9.4f} {r['S_Neel']:>9.4f} {r['m2']:>8.4f}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
