import numpy as np
import os
import re
try:
    from numba import njit          # optional speedup
except ImportError:                 # server without numba: plain python is fine
    def njit(f):
        return f

# ================================================================
# Geometry helpers — UNCHANGED from wf_tt.py
# ================================================================
def cart_coord(i, Lx):
    ix = i % Lx
    iy = i // Lx + 1
    if i % Lx == 0:
        iy -= 1
        ix = Lx
    return ix, iy

def linear_index(ix, iy, Lx):
    return int((iy - 1) * Lx + ix)

def Neighbor(index, Lx, Ly):
    ix, iy = cart_coord(index, Lx)

    up_y = iy - 1
    if up_y < 1:
        up_y = Ly
    up = linear_index(ix, up_y, Lx)

    down_y = iy + 1
    if down_y > Ly:
        down_y = 1
    down = linear_index(ix, down_y, Lx)

    left_x = ix - 1
    if left_x < 1:
        left_x = Lx
    left = linear_index(left_x, iy, Lx)

    right_x = ix + 1
    if right_x > Lx:
        right_x = 1
    right = linear_index(right_x, iy, Lx)

    # NNN diagonal neighbors
    up_right_x, up_right_y = ix + 1, iy - 1
    if up_right_x > Lx: up_right_x = 1
    if up_right_y < 1:  up_right_y = Ly
    up_right = linear_index(up_right_x, up_right_y, Lx)

    down_right_x, down_right_y = ix + 1, iy + 1
    if down_right_x > Lx: down_right_x = 1
    if down_right_y > Ly: down_right_y = 1
    down_right = linear_index(down_right_x, down_right_y, Lx)

    up_left_x, up_left_y = ix - 1, iy - 1
    if up_left_x < 1:  up_left_x = Lx
    if up_left_y < 1:  up_left_y = Ly
    up_left = linear_index(up_left_x, up_left_y, Lx)

    down_left_x, down_left_y = ix - 1, iy + 1
    if down_left_x < 1:  down_left_x = Lx
    if down_left_y > Ly: down_left_y = 1
    down_left = linear_index(down_left_x, down_left_y, Lx)

    return up, down, left, right, up_right, down_right, up_left, down_left


# ================================================================
# GetK — UNCHANGED from wf_tt.py
#
# One function handles BOTH spins via sign convention:
#   spin-up:   GetK(Lx, Ly, t, +tA, +tt)
#   spin-down: GetK(Lx, Ly, t, -tA, -tt)
#
# Resulting dispersion:
#   ε↑(k) = -2(1+tA)cos(kx) - 2(1-tA)cos(ky) + 4tt·sin(kx)sin(ky)
#   ε↓(k) = -2(1-tA)cos(kx) - 2(1+tA)cos(ky) - 4tt·sin(kx)sin(ky)
#
# Limiting cases:
#   tt=0       → CPL (NN only) model
#   tA=0       → PRB (NNN only) model
#   tA≠0,tt≠0  → Unified model  ← NEW
# ================================================================
# ================================================================
# Lieb coordinate helpers -- FORTRAN ordering, deliberately NOT cart_coord
#
# cart_coord/linear_index use  i = (iy-1)*Lx + ix   (ix is the FAST index).
# mc2duph.f90 uses             i = (ix-1)*Ly + iy   (iy is the FAST index).
# The two are transposes of each other. That never mattered for the
# checkerboard, whose sublattice test (ix+iy)%2 is symmetric under swapping
# the axes, so a transposed trial was still a valid checkerboard trial. On
# Lieb it matters: A, B and C distinguish x from y, and a transposed trial
# would swap the B and C sublattices. So the Lieb code uses the Fortran
# ordering everywhere and never calls cart_coord.
# ================================================================
def lieb_sites(Lx, Ly):
    """The real sites, in the Fortran loop order.

    mc2duph.f90 runs ixx over 1..lx outermost and iyy over 1..ly innermost and
    increments the site counter only when the position is not (odd, odd), so the
    fine-grid quarter that carries no site is simply never indexed. Anything that
    builds a trial wave function has to walk the lattice in the same order, or
    every row of phiT lands on the wrong site.

    Returns a list of (ix, iy), 1-based, whose position in the list is the
    0-based Fortran site index.
    """
    return [(ix, iy) for ix in range(1, Lx + 1) for iy in range(1, Ly + 1)
            if not (ix % 2 == 0 and iy % 2 == 0)]


def lieb_maps(Lx, Ly):
    """(sites, pos->index). pos is wrapped, so lookups may cross the boundary.
    A position with no site is absent from the map, which is the Python echo of
    the Fortran iposit==0 sentinel."""
    sites = lieb_sites(Lx, Ly)
    pos = {p: i for i, p in enumerate(sites)}
    return sites, pos


def lieb_nsites(Lx, Ly):
    return 3 * (Lx // 2) * (Ly // 2)


def GetK(Lx, Ly, t, tprime, delta):
    """Lieb-lattice kinetic matrix. MUST match mc2duph.f90 element for element.

    The Lieb lattice is a square lattice with a quarter of the sites removed.
    On a fine grid of spacing 1/2 (Lx, Ly are TWICE the cell count and EVEN):

        (odd ,odd ) = A     (even,odd ) = B
        (odd ,even) = C     (even,even) = no site

    with 1-based ix, iy, so A-B and A-C are the fine-grid AXIAL bonds and the
    next-nearest B-C bonds are the fine-grid DIAGONALS.

        A : +-x is B, +-y is C          diagonals empty -> none
        B : +-x is A, +-y empty         all four diagonals are C
        C : +-y is A, +-x empty         all four diagonals are B

    Argument names are kept from the checkerboard version so the caller and the
    folder-name parsing do not change:
        t       -> unused (the amplitude comes from tprime, matching in.dat t0)
        tprime  -> NEAREST neighbour, in.dat slot t0
        delta   -> NEXT nearest (B-C only), in.dat slot t1

    There is no vacancy term of any kind. An earlier version kept the missing
    quarter as orbitals pushed above the Fermi level by a large diagonal, and
    that diagonal broke the walker weights, so the sites are gone instead.
    """
    sites, pos = lieb_maps(Lx, Ly)
    N = len(sites)
    K = np.zeros((N, N))

    def put(i, ix, iy, amp):
        j = pos.get((((ix - 1) % Lx) + 1, ((iy - 1) % Ly) + 1))
        if j is not None:
            K[i, j] = amp

    for i, (ix, iy) in enumerate(sites):
        px, py = ix % 2, iy % 2         # 1 on the odd (A-carrying) rows/columns

        if py == 1:                     # A or B: the +-x neighbour is a site
            put(i, ix + 1, iy, tprime)
            put(i, ix - 1, iy, tprime)

        if px == 1:                     # A or C: the +-y neighbour is a site
            put(i, ix, iy + 1, tprime)
            put(i, ix, iy - 1, tprime)

        if px + py == 1:                # B or C: all four diagonals are B-C
            for dx, dy in ((1, 1), (-1, -1), (1, -1), (-1, 1)):
                put(i, ix + dx, iy + dy, delta)

    return K


# ================================================================
# Hartree-Fock helpers — UNCHANGED from wf_tt.py
# ================================================================
def hubbard_u(Lx, Ly, U):
    """Per-site U. Zero on A, U on B and C, matching mc2duph.f90 and Eq. (1) of
    the published model. A uniform U is a DIFFERENT Hamiltonian and does not
    show the altermagnet."""
    return np.array([0.0 if (ix % 2 == 1 and iy % 2 == 1) else U
                     for ix, iy in lieb_sites(Lx, Ly)])

def GetHup(Lx, Ly, nup, ndn, K, U):
    return K + np.diag(ndn * U)

def GetHdn(Lx, Ly, nup, ndn, K, U):
    return K + np.diag(nup * U)

@njit
def GetEigen(H):
    eig_values, eig_vecs = np.linalg.eigh(H)
    idx_sorted = np.argsort(eig_values)
    return eig_values[idx_sorted], eig_vecs[:, idx_sorted]

def Getnup(upvecs, NUP):
    return np.sum(upvecs[:, :NUP] ** 2, axis=1)

def Getndn(dnvecs, NDN):
    return np.sum(dnvecs[:, :NDN] ** 2, axis=1)


# ================================================================
# Iteration — UNCHANGED from wf_tt.py
# ================================================================
def neel_seed(Lx, Ly, NUP, NDN):
    """Ferrimagnetic seed for the Lieb lattice: A up, B and C down.

    Lieb's theorem fixes this. The Lieb lattice is bipartite with sublattices
    {A} and {B, C}, of sizes Lc^2 and 2*Lc^2, so the half-filled Hubbard ground
    state carries S = ||A| - |B u C||/2 = Lc^2/2. The ordered state is therefore
    FERRIMAGNETIC, not the compensated Neel of the checkerboard, and the seed
    has to reflect that or the SCF starts in the wrong sector.
    """
    sites = lieb_sites(Lx, Ly)
    N = len(sites)
    nup = np.zeros(N); ndn = np.zeros(N)
    for i, (ix, iy) in enumerate(sites):
        if ix % 2 == 1 and iy % 2 == 1:
            nup[i], ndn[i] = 0.9, 0.1       # A sublattice, up
        else:
            nup[i], ndn[i] = 0.1, 0.9       # B and C sublattice, down
    nup *= NUP / nup.sum(); ndn *= NDN / ndn.sum()
    return nup, ndn


def am_seed(Lx, Ly, NUP, NDN):
    """Altermagnetic seed: A non-magnetic, B up, C down.

    This is the state the published work finds at n = 2 and n = 4 electrons per
    unit cell. B and C are related by the C4 rotation about an A site and NOT by
    any translation or by inversion, which is exactly the altermagnet criterion,
    so antialigning them is altermagnetic rather than antiferromagnetic. It is
    compensated, so it only exists at NUP = NDN, and the ferrimagnetic seed
    cannot reach it.
    """
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
    return nup, ndn


def Iteration(Lx, Ly, U, tA, tt, NUP=None, NDN=None, seed='random'):
    """One HF SCF attempt. seed: 'neel' | 'pm' | 'random'.
    CHECKERBOARD frustration makes random seeds wander at U>=2, so the SCF is
    capped: at most MAX_RESETS resets of 1000 iterations, then the best state
    reached is returned (never hangs)."""
    N     = lieb_nsites(Lx, Ly)
    t     = 1
    alpha = 0.5
    times = 0
    MAX_RESETS = 1
    resets = 0

    # LIEB: hopping is SPIN INDEPENDENT -> identical K for both spins.
    # Folder fields REPURPOSED: tt -> tprime (mean diagonal t'), tA -> delta.
    # Collinear Neel / d-wave AM is EMERGENT from the HF U*n term, not the band.
    # folder fields: tt -> NN amplitude (in.dat t0), tA -> B-C amplitude (t1)
    K_up   = GetK(Lx, Ly, t, tt, tA)
    if EPS_A != 0.0:                     # on-site energy on sublattice A
        for i, (ix, iy) in enumerate(lieb_sites(Lx, Ly)):
            if ix % 2 == 1 and iy % 2 == 1:
                K_up[i, i] += EPS_A
    K_down = K_up                        # spin independent
    Uscalar = U                          # keep the scalar for messages
    U = hubbard_u(Lx, Ly, U)             # per site: 0 on A, U on B and C

    if seed == 'neel':
        nup, ndn = neel_seed(Lx, Ly, NUP, NDN)
    elif seed == 'am':
        nup, ndn = am_seed(Lx, Ly, NUP, NDN)
    elif seed == 'pm':
        nup = np.full(N, NUP / N); ndn = np.full(N, NDN / N)
    else:
        nup = np.random.random(N); ndn = np.random.random(N)
    E1 = E2 = 0

    while True:
        if times == 1000:
            resets += 1
            if resets > MAX_RESETS:
                print(f"Warning: U={Uscalar} seed={seed}: {MAX_RESETS} resets exhausted, "
                      f"returning best-so-far E={E2:.6f}", flush=True)
                break
            nup = np.random.random(N)
            ndn = np.random.random(N)
            E1 = E2 = 0
            print(f"Warning: U={Uscalar} seed={seed} did not converge in 1000 steps. "
                  f"Resetting ({resets}/{MAX_RESETS}).", flush=True)
            times = 0

        Hup = GetHup(Lx, Ly, nup, ndn, K_up,   U)
        Hdn = GetHdn(Lx, Ly, nup, ndn, K_down,  U)

        eig_eup, eig_vup = GetEigen(Hup)
        eig_edn, eig_vdn = GetEigen(Hdn)

        # Convergence is still tracked on the eigenvalue sum, which is what the
        # SCF actually stationarises. E_hf below is the variational energy and is
        # what Run() compares across seeds.
        E  = sum(eig_eup[:NUP]) + sum(eig_edn[:NDN])
        E1 = E2
        E2 = E

        nup_new = (1 - alpha) * nup + alpha * Getnup(eig_vup, NUP)
        ndn_new = (1 - alpha) * ndn + alpha * Getndn(eig_vdn, NDN)

        nup = nup_new
        ndn = ndn_new

        if abs(E2 - E1) < 0.0001:
            break
        times += 1

    # Variational energy of the determinant the orbitals actually describe.
    # Deliberately NOT E - sum(U n_up n_dn): that identity assumes exact self
    # consistency, which the mixed-density SCF does not reach for states that
    # only converge loosely, and it then flatters exactly those states.
    r_u = eig_vup[:, :NUP] @ eig_vup[:, :NUP].T
    r_d = eig_vdn[:, :NDN] @ eig_vdn[:, :NDN].T
    n_u = Getnup(eig_vup, NUP)
    n_d = Getndn(eig_vdn, NDN)
    E_hf = float(np.sum(K_up * r_u) + np.sum(K_up * r_d) + np.sum(U * n_u * n_d))
    return nup, ndn, E_hf, eig_vup, eig_vdn


# ================================================================
# Run — UNCHANGED from wf_tt.py
# ================================================================
def Run(Lx, Ly, U=None, NUP=None, NDN=None, tA=None, tt=None):
    N    = lieb_nsites(Lx, Ly)
    # Attempt schedule: the physical Neel seed first (converges fast on the
    # checkerboard and is the expected ground state), a paramagnetic seed as
    # the unbiased control, then random restarts. Best energy wins, as before.
    seeds = ['am', 'neel', 'pm'] + ['random'] * 2
    nitr = len(seeds)
    phiT_ups, phiT_dns = [], []
    Es = np.zeros(nitr)

    for i in range(nitr):
        nup, ndn, E, eig_vup, eig_vdn = Iteration(
            Lx, Ly, U, tA, tt, NUP=NUP, NDN=NDN, seed=seeds[i]
        )
        print(f"attempt {i+1}/{nitr} (seed={seeds[i]}): E_HF = {E:.6f}", flush=True)
        Es[i] = E
        phiT_ups.append(eig_vup[:, :NUP].T)
        phiT_dns.append(eig_vdn[:, :NDN].T)

    idx    = np.argmin(Es)
    E_gs   = Es[idx]
    phiT_up = phiT_ups[idx]
    phiT_dn = phiT_dns[idx]
    return phiT_up, phiT_dn, E_gs


# ================================================================
# extract_variables — CHANGED: now reads BOTH tA and tt
#
# Folder naming convention for unified model:
#   L{L}n{n}u{u}tA{tA}tt{tt}N{N}
#   e.g.  L18n0.800u4.0tA0.30tt0.20N260
#
# Backward compatible:
#   CPL folders (no tt):  tt defaults to 0.0
#   PRB folders (no tA):  tA defaults to 0.0
# ================================================================
def extract_variables(folder_name):
    L = int(re.search(r'-?\d+(\.\d+)?',   folder_name.split("L")[1]).group())
    n = float(re.search(r'-?\d+(\.\d+)?', folder_name.split("n")[1]).group())
    u = float(re.search(r'-?\d+(\.\d+)?', folder_name.split("u")[1]).group())
    N = float(re.search(r'-?\d+(\.\d+)?', folder_name.split("N")[1]).group())

    # CHANGE 1: read tA (default 0 for PRB-only folders)
    tA = float(re.search(r'-?\d+(\.\d+)?',
               folder_name.split("tA")[1]).group()) if "tA" in folder_name else 0.0

    # CHANGE 2: read tt (default 0 for CPL-only folders)
    tt = float(re.search(r'-?\d+(\.\d+)?',
               folder_name.split("tt")[1]).group()) if "tt" in folder_name else 0.0

    return L, n, u, tA, tt, N


# ================================================================
# Main
# ================================================================
root           = os.getcwd()
current_folder = os.path.basename(root)
variables      = extract_variables(current_folder)

if variables:
    L, n, u, tA, tt, N = variables   # CHANGE 3: unpack tA from variables
else:
    raise RuntimeError(f"Could not parse folder name: {current_folder}")

# eps_A is in.dat line 5 field 3, the slot the Fortran reads as t2.
EPS_A = float(open(os.path.join(root, "in.dat")).readlines()[4].split(",")[2])

print(f"Folder : {current_folder}")
print(f"Params : L={L}, n={n}, U={u}, tA={tA}, tt={tt}, N={int(N)}, eps_A={EPS_A}")

# CHANGE 4: pass actual tA (not hardcoded 0) to Run
# NUP and NDN come from parameter.f90, NOT from the folder name. The name
# carries a single N, but the Lieb ground state is ferrimagnetic (Lieb's
# theorem: S = ||A|-|BuC||/2 = cells^2/2), so NUP != NDN. Parsing one N from the
# name silently built a trial with NUP=NDN while the Fortran expected NUP!=NDN,
# and the mismatched down determinant made the down-spin weight collapse.
import re as _re
_p = open(os.path.join(root, "parameter.f90")).read()
_m = _re.search(r"parameter *\( *NUP *= *(\d+) *, *NDN *= *(\d+)", _p)
if not _m:
    raise RuntimeError("could not read NUP, NDN from parameter.f90")
NUP_F, NDN_F = int(_m.group(1)), int(_m.group(2))
print(f"Filling : NUP={NUP_F}, NDN={NDN_F} (read from parameter.f90; folder N={int(N)})")
if NUP_F != NDN_F:
    print(f"          ferrimagnetic sector, S = {(NUP_F-NDN_F)/2}")

phiT_up, phiT_dn, E_gs = Run(L, L, U=u, NUP=NUP_F, NDN=NDN_F, tA=tA, tt=tt)
print(f"HF ground state energy: {E_gs:.6f}   (variational, not the eigenvalue sum)")

phiT_up[np.abs(phiT_up) < 1e-10] = 0.0
phiT_dn[np.abs(phiT_dn) < 1e-10] = 0.0

np.savetxt(os.path.join(root, 'wfup.txt'), phiT_up,
           fmt='%.12e', delimiter=' ', header='', comments='')
np.savetxt(os.path.join(root, 'wfdn.txt'), phiT_dn,
           fmt='%.12e', delimiter=' ', header='', comments='')

print("Saved: wfup.txt  wfdn.txt")
