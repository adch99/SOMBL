#!/bin/bash
tmpfile="logs/running_jobnos.log"
qstat -u mursalin | tail -n+6 | grep " R " | cut -f 1 -d" " > $tmpfile
for runno in $(cat $tmpfile); do
    echo $(qstat -f $runno | grep "Job_Name" | cut -d"=" -f 2);
done
