#!/bin/bash
# Move the finished L=14 runs onto 251 (the archive box).
# NOTE: the receiving ssh must NOT have "< /dev/null" -- that replaces the piped
# tar stream with an empty stdin and gzip dies with "unexpected end of file".
# Use ssh -n on the SENDING side instead when stdin must be detached.
set -u
for b in 250 252 253 254; do
  echo "=== $b start $(date +%H:%M:%S) ==="
  ssh -n -o BatchMode=yes $b 'cd $HOME/campaign_L14L16/runs && tar czf - */' \
    | ssh -o BatchMode=yes 251 'tar xzf - -C $HOME/collected/L14_campaign/fortran'
  echo "=== $b done $(date +%H:%M:%S)  total on 251: $(ssh -n -o BatchMode=yes 251 'ls -d $HOME/collected/L14_campaign/fortran/*/ 2>/dev/null | wc -l') ==="
done
ssh -n -o BatchMode=yes 251 'cd $HOME/campaign_L14L16/runs && for d in */; do [ -d "$HOME/collected/L14_campaign/fortran/$d" ] || cp -r "$d" "$HOME/collected/L14_campaign/fortran/"; done; echo "251 local merged"'
ssh -n -o BatchMode=yes 251 'echo "TOTAL run dirs: $(ls -d $HOME/collected/L14_campaign/fortran/*/ | wc -l)   with corr.dat: $(ls $HOME/collected/L14_campaign/fortran/*/corr.dat | wc -l)   size: $(du -sh $HOME/collected/L14_campaign/fortran | cut -f1)"'
