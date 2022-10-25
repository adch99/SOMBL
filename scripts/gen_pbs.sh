#!/bin/bash

# Author: Aditya Chincholi
# All rights reserved.

# Script to generate a PBS file for input
# in the cluster

size=40
runs=100
hopup=1.0
hopdn=1.0

dis_high=18
dis_low=8.0
dis_step=1.0

coup_high=2.0
coup_low=0.0
coup_step=0.1

couplings=$(seq $coup_low $coup_step $coup_high)
disorders=$(seq $dis_high -$dis_step $dis_low)

for coupling in $couplings; do
    for disorder in $disorders; do
        jobname="mbl_${size}x${size}_W${disorder}_C${coupling}_TU${hopup}_TD${hopdn}_N${runs}"
        pbs_filename="jobs/${jobname}.pbs"
        log_filename="logs/${jobname}.log"
        cat <<END_OF_PROGRAM > ${pbs_filename}
#!/bin/sh
# The following line specifies the maximum cpu utilization
# in terms of cpu time (3 days)
#PBS -l cput=72:00:00
# The following line specifies the maximum cpu utilization
# in terms of wall time (3 days)
#PBS -l walltime=72:00:00
# The following line merges the stdout and stderr together
#PBS -j oe
#PBS -l mem=8gb
# The following line specifies the name of the job
#PBS -N ${jobname}
    job=${jobname}
cd \${PBS_O_WORKDIR}
#########################################################################

echo \`hostname\` \`date\` \`pwd\` >> ./trial\$job
(time ./build/exact_diag_simulation -n \$runs -c \$coupling -w \$disorder -s \$size -u \$hopup -d \$hopdn) >> ./trial\$job
echo \`date\` >> ./trial\$job

cd \${PBS_O_WORKDIR}
date >> ./trial\$job
echo Job ended at \`date\` >> ./trial\$job
exit 0

END_OF_PROGRAM

    done
done
