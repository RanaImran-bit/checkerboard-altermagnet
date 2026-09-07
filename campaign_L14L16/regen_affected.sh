#!/bin/bash
# Re-run the 26 half-filling cells the trial audit found were constrained to a
# variationally worse trial, now with the corrected selection. These feed the
# finite-size extrapolation behind the "no long-range order" claim, and the old
# bias was uneven across L (7 cells at L=8, 7 at L=10, 8 at L=12, 3 at L=14, 1 at
# L=18), which is exactly what can corrupt an extrapolation.
# 4 seeds each, one core per run.
set -u
D=$HOME/campaign_L14L16
cd $D || exit 1
: > /tmp/regen_work.txt
tail -n +2 /tmp/regen_cells.csv | while IFS=, read -r L U DE; do
  L=${L%.0}; N=$((L*L/2))
  for S in 1 2 3 4; do echo "$L $U $DE $N $S" >> /tmp/regen_work.txt; done
done
echo "runs queued: $(wc -l < /tmp/regen_work.txt)"
# build any binary we do not have yet
for L in 8 10 12 14 18; do
  N=$((L*L/2))
  [ -f bin/CPMC_L${L}_N${N}.exe ] || bash campaign.sh build $L $N
done
ls bin/*.exe | wc -l | xargs echo "binaries:"
