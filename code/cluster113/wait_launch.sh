#!/bin/bash
# Wait for the L=8 test job to finish, verify it produced real output, then
# submit the L=8 and L=10 grids. Aborts rather than submitting 98 jobs if the
# test failed, since every job would fail the same way.
B=/opt/data/hpc/ZONES/all/PROJECTS/a1/PRIVATE/jianyu/Imran/checkerboard
JOB=52566
LOG=$B/wait_launch.log
echo "[$(date +%H:%M)] waiting for test job $JOB" > $LOG
while squeue -j $JOB -h -o %T 2>/dev/null | grep -qE "RUNNING|PENDING"; do sleep 120; done
D=$(ls -d $B/runs/L8n1.000* 2>/dev/null | head -1)
echo "[$(date +%H:%M)] test job left the queue; dir=$D" >> $LOG
ok=1
tail -3 "$D/nohup.out" 2>/dev/null | grep -q time_sec || { echo "  FAIL: no time_sec" >> $LOG; ok=0; }
[ -d "$D/dir-kVals" ] && [ -d "$D/dir-rVals" ] || { echo "  FAIL: missing dir-kVals/dir-rVals" >> $LOG; ok=0; }
if [ "$ok" = "1" ]; then
  k=$(ls "$D/dir-kVals" | wc -l); r=$(ls "$D/dir-rVals" | wc -l)
  t=$(tail -3 "$D/nohup.out" | grep -o "time_sec[^0-9]*[0-9.]*" | grep -o "[0-9.]*$")
  echo "  OK: dir-kVals=$k dir-rVals=$r time_sec=$t" >> $LOG
  cd $B
  bash launch_grid.sh 8  >> $LOG 2>&1
  bash launch_grid.sh 10 >> $LOG 2>&1
  echo "[$(date +%H:%M)] queue now: $(squeue -u $USER -h | wc -l) jobs" >> $LOG
else
  echo "[$(date +%H:%M)] ABORTED, grids NOT submitted" >> $LOG
fi
