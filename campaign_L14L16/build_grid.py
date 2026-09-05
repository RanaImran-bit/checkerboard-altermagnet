# The L=14 / L=16 campaign grid, rebuilt with the corrected trial.
#
# Every cell is re-run because the audit showed the trial selection was wrong in
# 36% of runs and, at half filling, always in the direction of OVERSTATING the
# staggered moment. Comparisons against the old numbers are only meaningful if
# the whole grid is regenerated the same way.
#
#   sizes    L = 14 (196 sites), L = 16 (256 sites)
#   U        0, 2, 3, 3.5, 4, 4.5, 5
#   delta    0.1 .. 0.7          (in.dat tam = -delta, tt = t1 = 0.3)
#   filling  n = 1 plus 8 doped targets spanning 0.5 to 0.95
#
# NUP = NDN = round(n * nsites / 2), so the achieved n is snapped to the lattice.
import csv, sys

L_SET     = [14, 16]
U_SET     = [0.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0]
DELTA_SET = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
N_TARGET  = [1.00, 0.95, 0.90, 0.82, 0.75, 0.70, 0.65, 0.60, 0.50]

rows = []
for L in L_SET:
    ns = L * L
    for nt in N_TARGET:
        N = round(nt * ns / 2)
        n_act = 2 * N / ns
        for U in U_SET:
            for de in DELTA_SET:
                rows.append(dict(L=L, nsites=ns, n_target=nt, n=round(n_act, 4),
                                 N=N, U=U, delta=de))
with open("grid.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print(f"total runs: {len(rows)}")
print(f"  per (L, filling) block: {len(U_SET)*len(DELTA_SET)}")
print(f"\n{'L':>4}{'n target':>10}{'N':>6}{'n achieved':>12}{'runs':>7}")
for L in L_SET:
    for nt in N_TARGET:
        N = round(nt*L*L/2); print(f"{L:>4}{nt:>10.2f}{N:>6}{2*N/(L*L):>12.4f}"
                                   f"{len(U_SET)*len(DELTA_SET):>7}")
