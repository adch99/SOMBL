#!/bin/bash

# Author: Aditya Chincholi
# All rights reserved.

# Script to generate a PBS file for input
# in the cluster.
# Average the green's functions from the
# different batches.

size=60
runs=100
batchsize=10
numbatches=10
bins=50
hopup=1.0
hopdn=1.0

dis_high=12
dis_low=10
dis_step=1.0

coup_high=1.7
coup_low=1.5
coup_step=0.1

couplings=$(seq $coup_low $coup_step $coup_high)
disorders=$(seq $dis_high -$dis_step $dis_low)
batchnums=$(seq 1 1 $numbatches)

for coupling in $couplings; do
    for disorder in $disorders; do
        # jobname="mbl_average_${size}x${size}_W${disorder}_C${coupling}_TU${hopup}_TD${hopdn}_N${runs}_BS${batchsize}"
        jobname="mbl_density_${size}x${size}_W${disorder}_C${coupling}_TU${hopup}_TD${hopdn}_N${runs}"
        pbs_filename="jobs/${jobname}.pbs"
        log_filename="logs/${jobname}.log"
        # execname="build/batch_average"
        # execname="build/keldysh_energy_batch_average"
        execname="build/keldysh_densities"
        echo "Out: $pbs_filename"
        cat <<END_OF_PROGRAM > ${pbs_filename}
#!/bin/bash
# The following line specifies the maximum cpu utilization
# in terms of cpu time (3 days)
#PBS -l cput=72:00:00
# The following line specifies the maximum cpu utilization
# in terms of wall time (3 days)
#PBS -l walltime=72:00:00
# The following line merges the stdout and stderr together
#PBS -j oe
#PBS -o ${log_filename}
#PBS -l mem=4gb
# The following line specifies the name of the job
#PBS -N ${jobname}
#########################################################################

job=${jobname}
cd \${PBS_O_WORKDIR}
hostname
date
echo "Running job ${jobname}"
pwd
time ${execname} -n $runs --batchsize $batchsize --bins $bins -c $coupling -w $disorder -s $size -u $hopup -d $hopdn
date

cd \${PBS_O_WORKDIR}
date
echo Job ended at \`date\`
# Use this only if there's a separate executable
# for each job.
# rm ${execname}
exit 0

END_OF_PROGRAM
    done
done
