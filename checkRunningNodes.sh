# check the nodes where you running jobs is locating in tscc,ucsd.

qstat -u shg047 | awk '{print $1,$10}' | grep R| grep local | awk '{print $1} | xargs -n 1 checkjob | grep tscc
