#!/bin/bash
# Archive every FINISHED L=14 run onto 251, resumably.
#
# Only runs with corr.dat are copied, so a job still writing is never half
# archived. Runs already on 251 are skipped, so this can be re-run after an
# interruption (the first attempt died with the session at 105 of 592).
#
# The receiving ssh must NOT carry "< /dev/null": that replaces the piped tar
# stream with empty stdin and gzip dies with "unexpected end of file". Use
# ssh -n on the SENDING side instead.
set -u
DEST='$HOME/collected/L14_campaign/fortran'
have=$(ssh -n -o BatchMode=yes 251 "for d in $DEST/*/; do [ -f \"\$d/corr.dat\" ] && basename \"\$d\"; done 2>/dev/null" | sort)
echo "complete on 251: $(echo "$have" | grep -c . )"
ssh -n -o BatchMode=yes 251 "for d in $DEST/*/; do [ -f \"\$d/corr.dat\" ] || rm -rf \"\$d\"; done; echo pruned-incomplete"

copy_box () {                       # $1 = box, $2 = remote runs dir
  local b=$1 R=$2
  local done_list
  done_list=$(ssh -n -o BatchMode=yes "$b" "cd $R 2>/dev/null && for d in */; do [ -f \"\$d/corr.dat\" ] && basename \"\$d\"; done" | sort)
  local todo
  todo=$(comm -23 <(echo "$done_list") <(echo "$have"))
  local n=$(echo "$todo" | grep -c .)
  echo "=== $b: finished=$(echo "$done_list" | grep -c .)  to copy=$n  $(date +%H:%M:%S) ==="
  [ "$n" -eq 0 ] && return 0
  echo "$todo" | ssh -o BatchMode=yes "$b" "cat > /tmp/tocopy.txt; cd $R && tar czf - -T /tmp/tocopy.txt" \
    | ssh -o BatchMode=yes 251 "tar xzf - -C $DEST"
  echo "=== $b done $(date +%H:%M:%S)  on 251 now: $(ssh -n -o BatchMode=yes 251 "ls -d $DEST/*/ 2>/dev/null | wc -l") ==="
  have=$(printf '%s\n%s\n' "$have" "$todo" | sort)
}

for b in 250 252 253 254; do copy_box "$b" '$HOME/campaign_L14L16/runs'; done
copy_box 113new '$HOME/private/a1/campaign_L14L16/runs'
# 251's own runs are local
ssh -n -o BatchMode=yes 251 "cd \$HOME/campaign_L14L16/runs && for d in */; do
  [ -f \"\$d/corr.dat\" ] || continue
  [ -d \"$DEST/\$d\" ] || cp -r \"\$d\" $DEST/
done; echo '251 local merged'"
ssh -n -o BatchMode=yes 251 "echo TOTAL: \$(ls -d $DEST/*/ | wc -l) dirs, \$(ls $DEST/*/corr.dat | wc -l) with corr.dat, \$(du -sh $DEST | cut -f1)"
