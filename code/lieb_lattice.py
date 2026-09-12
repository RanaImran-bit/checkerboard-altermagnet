"""Lieb lattice: build the hopping matrix and verify against known band structure.

Lieb = square lattice of A sites with B sites on the horizontal bond midpoints and
C sites on the vertical ones. 3 sites per cell, NN hopping connects A-B and A-C
only. Bloch Hamiltonian:

    H(k) = [[0,            2t cos(kx/2), 2t cos(ky/2)],
            [2t cos(kx/2), 0,            0           ],
            [2t cos(ky/2), 0,            0           ]]

Eigenvalues 0 and +/- 2t sqrt(cos^2(kx/2)+cos^2(ky/2)) -- so a FLAT BAND at E=0
with one state per unit cell. That flat band is the signature to check against;
if the real-space matrix reproduces it, the construction is right.

Why this lattice for us: B and C are related by a 90 degree ROTATION (x <-> y),
which is exactly the altermagnet criterion, and the hopping is spin-independent,
so there is NO spin splitting at U=0. Any splitting must come from interaction.
"""
import numpy as np

def lieb_hopping(Lx, Ly, t=-1.0, apx=1, apy=1):
    """3*Lx*Ly hopping matrix. Site order per cell: A=0, B=1 (x-bond), C=2 (y-bond)."""
    n = 3*Lx*Ly
    K = np.zeros((n, n))
    def idx(x, y, s): return 3*((x % Lx)*Ly + (y % Ly)) + s
    def sgn(x, d, L, ap): return ap if (x + d < 0 or x + d >= L) else 1
    for x in range(Lx):
        for y in range(Ly):
            A = idx(x, y, 0); B = idx(x, y, 1); C = idx(x, y, 2)
            K[A, B] = K[B, A] = t                     # A - B within the cell
            K[A, C] = K[C, A] = t                     # A - C within the cell
            s = t*sgn(x, 1, Lx, apx)
            K[B, idx(x+1, y, 0)] = K[idx(x+1, y, 0), B] = s   # B - A next cell in x
            s = t*sgn(y, 1, Ly, apy)
            K[C, idx(x, y+1, 0)] = K[idx(x, y+1, 0), C] = s   # C - A next cell in y
    return K

for L in [6, 8]:
    K = lieb_hopping(L, L)
    E = np.linalg.eigvalsh(K)
    nflat = int(np.sum(np.abs(E) < 1e-10))
    # analytic spectrum on the same k-mesh
    ks = 2*np.pi*np.arange(L)/L
    an = [0.0]*(L*L)
    for kx in ks:
        for ky in ks:
            r = 2*abs(np.sqrt(np.cos(kx/2)**2 + np.cos(ky/2)**2))
            an += [-r, +r]
    an = np.sort(np.array(an))
    print(f"L={L}: sites={3*L*L}  spectrum {E.min():+.4f}..{E.max():+.4f}")
    # L^2 flat-band states, PLUS 2 from the dispersive bands touching zero at
    # (pi,pi) where cos(kx/2)=cos(ky/2)=0. So L^2+2, not L^2.
    print(f"   states at E=0: {nflat}  (expect {L*L}+2 = {L*L+2})  -> "
          f"{'OK' if nflat==L*L+2 else 'MISMATCH'}")
    print(f"   matches analytic 3-band form: {np.allclose(np.sort(E), an, atol=1e-8)}")
