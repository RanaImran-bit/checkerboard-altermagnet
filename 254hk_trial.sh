for d in $HOME/Checkerboard_Model/L12n1.000*/; do
  nm=$(basename "$d")
  tail -3 "$d/nohup.out" 2>/dev/null | grep -q time_sec || continue
  [ -f "$d/wfup.txt" ] && [ -f "$d/wfdn.txt" ] || { echo "$nm NOWF"; continue; }
  if cmp -s "$d/wfup.txt" "$d/wfdn.txt"; then echo "$nm TRIAL_NONMAGNETIC"
  else echo "$nm TRIAL_BROKEN"; fi
done
