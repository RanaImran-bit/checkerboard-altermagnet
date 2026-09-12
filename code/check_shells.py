"""Which fillings does each code use, and are they open or closed shell?

Both codes pick the electron number the same way, nup = round(n * L^2 / 2)
(Python: checkerboard_pr.py N_TARGETS; Fortran: setup_run.sh line 27). The
difference is WHICH n values were chosen.

A filling is CLOSED shell if the free-electron spectrum has a gap at the Fermi
level, E[nup] - E[nup-1] > 0. Then the trial determinant is unique. Open shell
means a degenerate manifold sits at the Fermi level and the trial is one
arbitrary choice within it, which biases the constrained path.
"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import checkerboard as cb

DELTAS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
TOL = 1e-8


def gap(L, nup, delta):
    K = cb.checkerboard_hopping(L, L, -1.0, 0.3, -delta)
    e = np.linalg.eigvalsh(K)
    return e[nup] - e[nup - 1]


def report(tag, L, nups, deltas):
    print(f"\n=== {tag}:  L = {L} ===")
    print(f"{'nup':>5}{'n':>8}  " + "".join(f"{'d='+str(d):>9}" for d in deltas))
    for nup in nups:
        n = 2 * nup / (L * L)
        row = f"{nup:5d}{n:8.3f}  "
        for d in deltas:
            g = gap(L, nup, d)
            row += f"{('open' if g < TOL else f'{g:.3f}'):>9}"
        print(row)


# Python: n = 2k/36 for k in [9,10,12,14,16,18] -> round fractions, chosen for
# convenience, NOT for shell structure.
report("PYTHON fillings (round fractions 2k/36)", 12,
       [36, 40, 48, 56, 64, 72], DELTAS)

# Fortran: the n values that actually appear as run directories.
FOR10 = [25, 29, 33, 37, 39, 41, 43, 45, 47, 49]
FOR14 = [49, 57, 61, 69, 73, 75, 77, 83, 87, 91, 95, 98]
report("FORTRAN fillings", 10, FOR10, [0.0, 0.1, 0.2, 0.3, 0.4])
report("FORTRAN fillings", 14, FOR14, [0.0, 0.1, 0.2, 0.3, 0.4])

print("\n" + "=" * 72)
for tag, L, nups, ds in [("PYTHON L=12", 12, [36, 40, 48, 56, 64, 72], DELTAS),
                         ("FORTRAN L=10", 10, FOR10, [0.0, .1, .2, .3, .4]),
                         ("FORTRAN L=14", 14, FOR14, [0.0, .1, .2, .3, .4])]:
    tot = len(nups) * len(ds)
    closed = sum(1 for nup in nups for d in ds if gap(L, nup, d) > TOL)
    print(f"{tag:14s}: closed shell in {closed}/{tot} (filling, delta) combinations"
          f"  = {100*closed/tot:.0f}%")
