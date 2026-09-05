# Bring the Lieb port to the published altermagnetic model.
#
# Kaushal and Franz, PRL (arXiv:2412.16421), Eq. (1):
#
#   H = t1 sum_<ij> c+c  +  t2 sum_<<ij>> c+c
#       + U sum_{i in B,C} n_up n_dn  +  eps_A sum_{i in A} n_i  -  mu sum_i n_i
#
# Three things separate this from a plain Lieb Hubbard, and the port had none of
# them:
#
#   1. U acts ONLY on B and C. A is the oxygen and carries no interaction. This
#      is the whole mechanism: the moments live on B and C, which are related by
#      C4 and not by translation, so their antialignment is altermagnetic.
#   2. t2, the B-C next-nearest hopping, is REQUIRED to stabilise the state. It
#      is already wired as the in.dat t1 slot.
#   3. eps_A, an on-site energy on A. Taken from the previously unused in.dat t2
#      slot. DEFAULT 0, and it should stay 0 until the walker weights are checked
#      at eps_A /= 0: a large tk diagonal is exactly what killed W_Down in the
#      first Lieb port. eps_A = 0 already shows a large altermagnetic region at
#      n = 4 (their Fig. 2), so nothing is lost by starting there.
#
# The altermagnet is COMPENSATED, so it lives at NUP = NDN. A run with
# NUP - NDN fixed at cells^2 imposes a net moment and cannot find it.
import re, sys

mc = open("mc2duph.f90").read()

def sub1(s, old, new, tag):
    if s.count(old) != 1:
        sys.exit(f"[{tag}] expected 1 occurrence, found {s.count(old)}")
    return s.replace(old, new)

# -- U on B and C only, eps_A on A -------------------------------------------
mc = sub1(mc, """do i=1,nsites_cu
hub_u(i)=ud            ! every indexed site is a real Lieb site now
epsil(i)=zero
end do""",
"""do i=1,nsites_cu
if ( mod(ixv(i),2)==0 .and. mod(iyv(i),2)==0 ) then
   hub_u(i)=zero       ! sublattice A: the oxygen, no Hubbard U
else
   hub_u(i)=ud         ! sublattices B and C: the magnetic sites
end if
epsil(i)=zero
end do""", "hub_u")

# eps_A rides the unused in.dat t2 slot. epsil() is declared but never read by
# this code, so it goes on the tk diagonal, which tred2 and the energy sum both
# handle. Guarded because a large diagonal broke the walker weights before.
mc = sub1(mc, "do 120 i=1,nsites_cu\n   lbpx = mod( ixv(i), 2 )",
"""if ( abs(t2) > 1.0e-8 ) then
   do i=1,nsites_cu
      if ( mod(ixv(i),2)==0 .and. mod(iyv(i),2)==0 ) then
         tk(i,i,1)=t2   ! eps_A on sublattice A. CHECK W_Down before trusting
         tk(i,i,2)=t2   ! any run with this nonzero: see patch_kaushal_franz.py
      end if
   end do
end if

do 120 i=1,nsites_cu
   lbpx = mod( ixv(i), 2 )""", "eps_A")
open("mc2duph.f90","w").write(mc)

# -- trial wave function ------------------------------------------------------
wf = open("wf_lieb.py").read()

wf = sub1(wf, """def GetHup(Lx, Ly, nup, ndn, K, U):
    return K + np.diag(ndn) * U

def GetHdn(Lx, Ly, nup, ndn, K, U):
    return K + np.diag(nup) * U""",
"""def hubbard_u(Lx, Ly, U):
    \"\"\"Per-site U. Zero on A, U on B and C, matching mc2duph.f90 and Eq. (1) of
    the published model. A uniform U is a DIFFERENT Hamiltonian and does not
    show the altermagnet.\"\"\"
    return np.array([0.0 if (ix % 2 == 1 and iy % 2 == 1) else U
                     for ix, iy in lieb_sites(Lx, Ly)])

def GetHup(Lx, Ly, nup, ndn, K, U):
    return K + np.diag(ndn * U)

def GetHdn(Lx, Ly, nup, ndn, K, U):
    return K + np.diag(nup * U)""", "GetH")

wf = sub1(wf, "    nup *= NUP / nup.sum(); ndn *= NDN / ndn.sum()\n    return nup, ndn",
"""    nup *= NUP / nup.sum(); ndn *= NDN / ndn.sum()
    return nup, ndn


def am_seed(Lx, Ly, NUP, NDN):
    \"\"\"Altermagnetic seed: A non-magnetic, B up, C down.

    This is the state the published work finds at n = 2 and n = 4 electrons per
    unit cell. B and C are related by the C4 rotation about an A site and NOT by
    any translation or by inversion, which is exactly the altermagnet criterion,
    so antialigning them is altermagnetic rather than antiferromagnetic. It is
    compensated, so it only exists at NUP = NDN, and the ferrimagnetic seed
    cannot reach it.
    \"\"\"
    sites = lieb_sites(Lx, Ly)
    N = len(sites)
    nup = np.zeros(N); ndn = np.zeros(N)
    for i, (ix, iy) in enumerate(sites):
        if ix % 2 == 1 and iy % 2 == 1:
            nup[i], ndn[i] = 0.5, 0.5       # A, non-magnetic
        elif iy % 2 == 1:
            nup[i], ndn[i] = 0.9, 0.1       # B, up
        else:
            nup[i], ndn[i] = 0.1, 0.9       # C, down
    nup *= NUP / nup.sum(); ndn *= NDN / ndn.sum()
    return nup, ndn""", "am_seed")

wf = sub1(wf, """    if seed == 'neel':
        nup, ndn = neel_seed(Lx, Ly, NUP, NDN)""",
"""    if seed == 'neel':
        nup, ndn = neel_seed(Lx, Ly, NUP, NDN)
    elif seed == 'am':
        nup, ndn = am_seed(Lx, Ly, NUP, NDN)""", "seed branch")

wf = sub1(wf, "    K_up   = GetK(Lx, Ly, t, tt, tA)\n    K_down = K_up                        # spin independent",
"""    K_up   = GetK(Lx, Ly, t, tt, tA)
    if EPS_A != 0.0:                     # on-site energy on sublattice A
        for i, (ix, iy) in enumerate(lieb_sites(Lx, Ly)):
            if ix % 2 == 1 and iy % 2 == 1:
                K_up[i, i] += EPS_A
    K_down = K_up                        # spin independent
    U = hubbard_u(Lx, Ly, U)             # per site: 0 on A, U on B and C""", "K + eps")

wf = sub1(wf, "    seeds = ['neel', 'pm'] + ['random'] * 5",
          "    seeds = ['am', 'neel', 'pm'] + ['random'] * 5", "seed list")

# eps_A comes from in.dat, the same slot the Fortran reads it from
wf = sub1(wf, 'print(f"Folder : {current_folder}")',
"""# eps_A is in.dat line 5 field 3, the slot the Fortran reads as t2.
EPS_A = float(open(os.path.join(root, "in.dat")).readlines()[4].split(",")[2])

print(f"Folder : {current_folder}")""", "eps read")
wf = sub1(wf, 'print(f"Params : L={L}, n={n}, U={u}, tA={tA}, tt={tt}, N={int(N)}")',
          'print(f"Params : L={L}, n={n}, U={u}, tA={tA}, tt={tt}, N={int(N)}, eps_A={EPS_A}")',
          "params print")
open("wf_lieb.py","w").write(wf)
print("patched mc2duph.f90 and wf_lieb.py")
