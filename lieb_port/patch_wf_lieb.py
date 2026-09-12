# Convert wf_lieb.py from the embedded 4*cells^2 grid to the reduced Lieb basis.
#
# The Fortran now indexes only the 3*cells^2 real sites, in the order
#     do ixx=1,lx; do iyy=1,ly; skip (odd,odd); i=i+1
# The trial wave function has to use exactly that ordering, or every row of
# phiT lands on the wrong site. EVAC is gone with the vacancies.
#
# Run from the directory holding wf_lieb.py.
import sys

w = open("wf_lieb.py").read()

old_helpers_start = w.index("def lieb_coord(i, Ly):")
old_helpers_end   = w.index("def GetK(Lx, Ly, t, tprime, delta):")
w = w[:old_helpers_start] + '''def lieb_sites(Lx, Ly):
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


''' + w[old_helpers_end:]

# ---- GetK on the reduced basis ----
gk_start = w.index("def GetK(Lx, Ly, t, tprime, delta):")
gk_end   = w.index("# ================================================================\n"
                   "# Hartree-Fock helpers")
w = w[:gk_start] + '''def GetK(Lx, Ly, t, tprime, delta):
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


''' + w[gk_end:]

# ---- ferrimagnetic seed on the reduced basis ----
ns_start = w.index("def neel_seed(Lx, Ly, NUP, NDN):")
ns_end   = w.index("def Iteration(Lx, Ly, U, tA, tt,")
w = w[:ns_start] + '''def neel_seed(Lx, Ly, NUP, NDN):
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


''' + w[ns_end:]

# ---- Iteration and Run: N is the reduced count, no vacancy masking ----
w = w.replace("""    N     = Lx * Ly
    t     = 1""", """    N     = lieb_nsites(Lx, Ly)
    t     = 1""")
w = w.replace("""        nup = np.random.random(N); ndn = np.random.random(N)
        # keep random seeds off the vacant sublattice
        for i in range(N):
            ix, iy = lieb_coord(i + 1, Ly)
            if (ix - 1) % 2 == 1 and (iy - 1) % 2 == 1:
                nup[i] = ndn[i] = 0.0""",
"""        nup = np.random.random(N); ndn = np.random.random(N)""")
w = w.replace("""def Run(Lx, Ly, U=None, NUP=None, NDN=None, tA=None, tt=None):
    N    = Lx * Ly""",
"""def Run(Lx, Ly, U=None, NUP=None, NDN=None, tA=None, tt=None):
    N    = lieb_nsites(Lx, Ly)""")

for bad in ("EVAC", "lieb_coord(", "lieb_index("):
    if bad in w:
        sys.exit(f"still referenced after patch: {bad}")

open("wf_lieb.py", "w").write(w)
print("wf_lieb.py converted to the reduced Lieb basis")
