#!/bin/bash
# Launch the two missing L=14 half-filling cells at U=2, delta=0.1 and 0.3.
#
# Why: the small-U peak in Delta_tot is strongest at delta=0.1 and 0.2, but L=14
# was only ever run at U=2 for delta=0.2. So the largest-lattice test of the peak
# exists in one column out of three. These two cells close that.
#
# Uses the existing ~/Checkerboard_Model/run_one.sh, which creates the folder,
# edits parameter.f90 (lx, ly, NUP, NDN) and in.dat (t2 = -delta, ud = U),
# generates the trial wave function, builds with the Intel toolchain, and starts
# mpirun in the background.
#
#   folder names produced:
#     L14n1.000u2.0tA-0.1000tt0.3N98
#     L14n1.000u2.0tA-0.3000tt0.3N98
#
# Comparable finished cells took 20.6 h (U=3, delta=0.1) and 21.4 h (U=3,
# delta=0.3), and run time falls slightly as U rises, so expect 21-23 h.
#
# 256 is deliberately not used: it has no `make`.
set -u
NP=${NP:-8}
run() {   # run <host> <delta>
  local h=$1 d=$2
  echo ">>> host $h   delta=$d   NP=$NP"
  ssh -n -o BatchMode=yes "$h" "
    cd ~/Checkerboard_Model || exit 1
    setsid nohup bash run_one.sh 14 2.0 $d 1.0 $NP \
      > launch_L14U2d${d}.log 2>&1 < /dev/null &
    echo '    detached, log: ~/Checkerboard_Model/launch_L14U2d${d}.log'"
}
run 250 0.1
run 251 0.3
