# Sanity check + U=0 reference energies for the reduced Lieb build.
#
# wf_lieb.py runs its driver at import, so the functions are pulled out by
# exec'ing the file up to the driver rather than importing it.
#
#   python check_lieb_K.py /path/to/wf_lieb.py
import sys, numpy as np

src = open(sys.argv[1] if len(sys.argv) > 1 else "wf_lieb.py").read()
cut = src.index("# ================================================================\n"
                "# extract_variables")
ns = {}
exec(compile(src[:cut], "wf_lieb.py", "exec"), ns)
GetK, lieb_sites = ns["GetK"], ns["lieb_sites"]

print(f"{'cells':>6}{'t1':>6}{'N':>5}{'bonds/site':>16}{'zero modes':>12}"
      f"{'E0(U=0)':>14}{'E0/N':>12}")
for cells in (2, 3, 4):
    for t1 in (0.0, 0.3):
        L = 2 * cells
        K = GetK(L, L, 1, -1.0, t1)
        N = K.shape[0]
        assert N == 3 * cells * cells, (N, cells)
        assert np.allclose(K, K.T), "K is not symmetric"
        sites = lieb_sites(L, L)
        # A has 4 neighbours (2 in x, 2 in y); B and C have 2 axial plus, when
        # t1 is on, 4 diagonals
        deg = (np.abs(K) > 1e-12).sum(1)
        nA = sum(1 for ix, iy in sites if ix % 2 == 1 and iy % 2 == 1)
        assert nA == cells * cells, (nA, cells)
        assert set(deg[[i for i, (ix, iy) in enumerate(sites)
                        if ix % 2 == 1 and iy % 2 == 1]]) == {4}, "A site degree"
        NUP, NDN = 2 * cells * cells, cells * cells
        e = np.linalg.eigvalsh(K)
        E0 = e[:NUP].sum() + e[:NDN].sum()
        nz = int(np.sum(np.abs(e) < 1e-10))
        print(f"{cells:>6}{t1:>6}{N:>5}{str(sorted(set(deg))):>16}{nz:>12}"
              f"{E0:>14.8f}{E0/N:>12.8f}")

print("\n  At t1 = 0 the flat band gives cells^2 exact zero modes, which is the")
print("  Lieb count. E0 is the U=0 ground state for NUP = 2*cells^2 up and")
print("  NDN = cells^2 down, the ferrimagnetic sector Lieb's theorem picks out,")
print("  and CPQMC at U=0 must reproduce it exactly.")
