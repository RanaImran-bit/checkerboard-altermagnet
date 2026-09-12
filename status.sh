#!/bin/bash
# status.sh -- one-shot status of every checkerboard / Lieb job on the lab machines.
#
#   bash ~/Desktop/checkerboard-altermagnet/status.sh
#
# Reads only: ps, squeue, and our own log/CSV files. Changes nothing, touches no
# other user's data. Safe to run as often as you like.
#
# Progress markers it knows about:
#   Fortran cell (delta=0 control) : blocks written out of 10; a finished cell ends
#                                    its nohup.out with a "time_sec" line
#   Delta_tot python driver        : rows in polarization_L12*.csv out of 36
#   Lieb m_am scan                 : "[n/42]" lines in mam_*.log
#   113new chi jobs                : squeue elapsed time; they emit NOTHING until done

NODES="250 251 252 253 254 255 256"
SSH="ssh -n -o ConnectTimeout=15 -o BatchMode=yes -o ControlMaster=no"
echo "=================== $(date '+%Y-%m-%d %H:%M') ==================="

echo
echo "### 1. Fortran  delta=0 control column  (L=12, 10 blocks per cell, ~11-13 h each)"
for h in $NODES; do
  out=$($SSH $h '
    cd ~/Checkerboard_Model 2>/dev/null || exit 0
    n=$(ps -u $USER -o comm= | grep -c "^CPMC.exe$")
    any=0
    for d in L12n1.000u*tA-0.0000*; do
      [ -d "$d" ] || continue
      any=1
      U=$(basename $d | sed "s/.*u//;s/tA.*//")
      blk=$(grep -cE "^ +[0-9]+ +-[0-9]" $d/nohup.out 2>/dev/null)
      if tail -3 $d/nohup.out 2>/dev/null | grep -q time_sec; then
        secs=$(tail -3 $d/nohup.out | grep time_sec | awk "{print \$2}")
        printf "    U=%-4s DONE (%.1f h)\n" "$U" "$(awk -v s=$secs "BEGIN{print s/3600}")"
      else
        age=$(( ($(date +%s) - $(stat -c %Y $d/nohup.out 2>/dev/null || echo 0)) / 60 ))
        printf "    U=%-4s running  blocks=%s/10  last write %s min ago\n" "$U" "$blk" "$age"
      fi
    done
    [ $any -eq 1 ] && echo "    (CPMC ranks alive on this node: $n)"
  ' 2>/dev/null)
  [ -n "$out" ] && { echo "  [$h]"; echo "$out"; }
done

echo
echo "### 2. Delta_tot python drivers  (36 jobs per run; CSV written only at the END)"
for h in $NODES; do
  out=$($SSH $h '
    cd ~/apbc_test 2>/dev/null || exit 0
    n=$(ps -u $USER -o args= | grep -c "[c]heckerboard_polarization.py")
    for c in polarization_L12*.csv; do
      [ -f "$c" ] && printf "    %-34s %s/36 rows\n" "$c" "$(( $(wc -l < $c) - 1 ))"
    done
    [ "$n" -gt 0 ] && echo "    still running: $n worker processes"
  ' 2>/dev/null)
  [ -n "$out" ] && { echo "  [$h]"; echo "$out"; }
done

echo
echo "### 3. Lieb m_am scans  (42 jobs each; these DO report progress live)"
for h in $NODES; do
  out=$($SSH $h '
    cd ~/lieb-altermagnet/code 2>/dev/null || exit 0
    for f in mam_*.log; do
      [ -f "$f" ] || continue
      d=$(grep -c "^  \[" $f)
      if grep -q "^wrote " $f; then
        printf "    %-20s DONE  (%s)\n" "$f" "$(grep "^wrote " $f | tail -1)"
      else
        printf "    %-20s %s/42  last: %s\n" "$f" "$d" "$(grep "^  \[" $f | tail -1 | sed "s/^  //")"
      fi
    done
  ' 2>/dev/null)
  [ -n "$out" ] && { echo "  [$h]"; echo "$out"; }
done

echo
echo "### 3b. Lieb S(q) structure factor  (unbroken; writes lieb_sq_partial.csv as it goes)"
for h in $NODES; do
  out=$($SSH $h '
    cd ~/lieb-altermagnet/code 2>/dev/null || exit 0
    for f in sq_*.log; do
      [ -f "$f" ] || continue
      tot=$(grep -o "= [0-9]* jobs" $f | head -1 | tr -dc 0-9)
      d=$(grep -c "^  \[" $f)
      if grep -q "^wrote " $f; then
        printf "    %-22s DONE  %s\n" "$f" "$(grep "^wrote " $f | tail -1)"
      else
        printf "    %-22s %s/%s  last: %s\n" "$f" "$d" "${tot:-?}" "$(grep "^  \[" $f | tail -1 | sed "s/^  //")"
      fi
    done
  ' 2>/dev/null)
  [ -n "$out" ] && { echo "  [$h]"; echo "$out"; }
done

echo
echo "### 3c. h-scan  (Delta_tot vs pinning field; CSV written only at the END)"
for h in $NODES; do
  out=$($SSH $h '
    cd ~/hscan 2>/dev/null || exit 0
    n=$(ps -u $USER -o args= | grep -c "[c]heckerboard_polarization.py")
    for f in hscan_L*.log; do
      [ -f "$f" ] || continue
      tot=$(grep -o "= [0-9]* jobs" "$f" | head -1 | tr -dc 0-9)
      if grep -q "^wrote " "$f"; then
        printf "    %-16s DONE  %s\n" "$f" "$(grep "^wrote " "$f" | tail -1)"
      else
        printf "    %-16s running (%s cells total, no progress output until done)\n" "$f" "${tot:-?}"
      fi
    done
    [ "$n" -gt 0 ] && echo "    workers: $n"
  ' 2>/dev/null)
  [ -n "$out" ] && { echo "  [$h]"; echo "$out"; }
done

echo
echo "### 3d. chi vertex finite-size  (L=8/10 to match the committed L=12 run)"
for h in $NODES; do
  out=$($SSH $h '
    cd ~/chi_fss 2>/dev/null || exit 0
    n=$(ps -u $USER -o args= | grep -c "[c]heckerboard_eqtime.py")
    for f in chifss_L*.log; do
      [ -f "$f" ] || continue
      if grep -q "^wrote " "$f"; then
        printf "    %-16s DONE  %s\n" "$f" "$(grep "^wrote " "$f" | tail -1)"
      else
        printf "    %-16s running\n" "$f"
      fi
    done
    [ "$n" -gt 0 ] && echo "    workers: $n"
  ' 2>/dev/null)
  [ -n "$out" ] && { echo "  [$h]"; echo "$out"; }
done

echo
echo "### 4. 113new SLURM  (our jobs are named py12u* / calib113)"
ssh -n -o ConnectTimeout=15 -o BatchMode=yes 113new '
  squeue -u $USER -o "%.9i %.14j %.2t %.11M %.5C %R" 2>/dev/null | grep -E "JOBID|py12|calib" || echo "  none of ours queued"
  B=/opt/data/hpc/ZONES/all/PROJECTS/a1/PRIVATE/jianyu/Imran/checkerboard
  echo "  --- finished chi CSVs ---"
  ls -la $B/pychi/eqtime_L12_*.csv 2>/dev/null | awk "{print \"    \"\$9\"  \"\$5\" bytes  \"\$6\" \"\$7\" \"\$8}" || echo "    (none yet)"
' 2>/dev/null

echo
echo "### 5. Node load  (free = cores - load; negative means oversubscribed)"
for h in $NODES; do
  $SSH $h 'printf "  %-6s cores=%s load=%-6s free=%s\n" "$(hostname -I | awk "{print \$1}" | cut -d. -f4)" "$(nproc)" "$(cut -d" " -f1 /proc/loadavg)" "$(awk -v c=$(nproc) -v l=$(cut -d" " -f1 /proc/loadavg) "BEGIN{printf \"%d\", c-l}")"' 2>/dev/null | sed "s/^  [0-9]*/  $h/"
done
echo "==============================================================="
