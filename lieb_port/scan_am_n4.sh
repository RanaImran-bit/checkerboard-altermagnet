#!/bin/bash
# Altermagnetic scan on the Lieb lattice at n = 4 electrons per unit cell.
#
# cells = 4, so 48 sites. Kaushal and Franz (PRL, arXiv:2412.16421) had
# unrestricted Hartree-Fock plus ED on 12- and 18-site clusters and no QMC at
# all, so this is 2.7 times their largest cluster and unbiased where HF is not.
#
# Two cuts through their Fig. 2(a), the t2 vs U phase diagram at eps_A = 0:
#   t2 scan at fixed U = 10, including t2 = 0 where they say the altermagnet
#     should NOT survive (the next-nearest hopping is required to stabilise it)
#   U scan at fixed t2 = 0.3, down to weak coupling
#
# Read out with check_sublattice_moment.py: m_A should stay at zero, and
# m_B - m_C is their staggered magnetisation, the altermagnetic order parameter.
set -u
cd "$(dirname "$0")/fortran" || exit 1
NP=${NP:-6}
run () { echo "=== cells=4 U=$1 t2=$2 ==="; bash run_lieb.sh 4 "$1" "$2" 4 "$NP"; }

for t2 in 0.0 0.1 0.2 0.3 0.4 0.5; do run 10.0 "$t2"; done   # t2 cut at U=10
for u  in 2.0 4.0 6.0 8.0;         do run "$u"  0.3;     done   # U cut at t2=0.3
echo; echo "launched. watch with:  tail -f */nohup.out"
