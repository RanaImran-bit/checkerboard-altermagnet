set -u
D=$HOME/campaign_L14L16; cd $D
# trials first (one per cell, shared by that cell's seeds)
awk '{print $1, $2, $3, $4}' /tmp/regen_work.txt | sort -u > /tmp/regen_cellargs.txt
echo "cells needing trials: $(wc -l < /tmp/regen_cellargs.txt)"
cat /tmp/regen_cellargs.txt | xargs -P 20 -n 4 sh -c 'cd $HOME/campaign_L14L16 && bash campaign.sh trial $0 $1 $2 $3' > /tmp/regen_trials.log 2>&1
echo "trials done: $(grep -c "trial ok\|trial exists" /tmp/regen_trials.log) failed: $(grep -c FAILED /tmp/regen_trials.log)"
# then the runs, one core each
awk '{L=$1;U=$2;DE=$3;N=$4;S=$5; t2=sprintf("%.4f",-DE); n=sprintf("%.3f",2*N/(L*L));
      printf "L%dn%su%stA%stt0.3N%d %d\n", L, n, U, t2, N, S}' /tmp/regen_work.txt > /tmp/regen_runs.txt
echo "runs: $(wc -l < /tmp/regen_runs.txt)"; head -2 /tmp/regen_runs.txt
nohup sh -c 'xargs -P 50 -n 2 sh -c "cd $HOME/campaign_L14L16 && bash campaign.sh runfg \$0 \$1" < /tmp/regen_runs.txt > /tmp/regen_run.log 2>&1; echo DONE > /tmp/regen.flag' >/dev/null 2>&1 &
sleep 25; echo "active: $(pgrep -cxu $USER CPMC.exe)"
