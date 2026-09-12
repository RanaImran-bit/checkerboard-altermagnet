#!/bin/bash
# The nine points that were stopped on 251, resubmitted to 113new Slurm.
# The tenth, U=10 t2=0.3, is still running on 251 and is deliberately not here.
set -eu
cd "$(dirname "$0")"
for t in 0.0 0.1 0.2 0.4 0.5; do bash submit_lieb.sh 4 10.0 "$t" 4; done   # t2 cut at U=10
for u in 2.0 4.0 6.0 8.0;     do bash submit_lieb.sh 4 "$u"  0.3 4; done   # U cut at t2=0.3
