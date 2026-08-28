#!/bin/bash
# Pull the multi-lattice-size sdwz.dat series for the finite-size test of magnetic order.
# R_p (peak sharpness at one L) cannot distinguish long-range order from a correlation
# length comparable to L. S(Q)/N against 1/L can. These series make that test possible:
#   n=1, U=4, delta=0.2 : L = 4, 12, 14, 16, 18
#   n=1, U=4, delta=0.3 : L = 8, 12, 14, 16, 18
set -u
OUT=${1:-../data/sdwz_fss}
mkdir -p "$OUT"
get() {  # get <host> <remote path> <label>
  scp -q -o BatchMode=yes "$1:$2" "$OUT/$3.dat" 2>/dev/null \
    && echo "  ok  $3" || echo "  MISSING  $3  ($1:$2)"
}
B=Checkerboard_Model
get 250 "~/$B/L4n1.000u4.0tA-0.2000tt0.3N8/dir-kVals/sdwz.dat"    d0.2_L4
get 251 "~/$B/L12n1.000u4tA-0.2000tt0.3N72/dir-kVals/sdwz.dat"    d0.2_L12
get 251 "~/$B/L14n1.000u4.0tA-0.2tt0.3N98/dir-kVals/sdwz.dat"     d0.2_L14
get 251 "~/$B/L16n1.000u4.0tA-0.2tt0.3N128/dir-kVals/sdwz.dat"    d0.2_L16
get 252 "~/$B/L18n1.000u4.0tA-0.2000tt0.3N162/dir-kVals/sdwz.dat" d0.2_L18

# delta = 0.3 series. NOTE: the only L=8 run is beta16_L8..., projected to beta=16
# rather than 32, so it is NOT part of this series and is deliberately omitted.
get 251 "~/$B/L12n1.000u4tA-0.3000tt0.3N72/dir-kVals/sdwz.dat"    d0.3_L12
get 251 "~/L14_archive/L14n1.000u4.0tA-0.3tt0.3N98/dir-kVals/sdwz.dat" d0.3_L14
get 250 "~/$B/L16n1.000u4.0tA-0.3000tt0.3N128/dir-kVals/sdwz.dat" d0.3_L16
get 250 "~/$B/L18n1.000u4.0tA-0.3000tt0.3N162/dir-kVals/sdwz.dat" d0.3_L18
