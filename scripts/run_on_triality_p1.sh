#!/bin/bash

# Author: Aditya Chincholi
# All rights reserved.

# Script to generate a PBS file for input
# in the cluster.
# Average the green's functions from the
# different batches.

size=100
runs=100
batchsize=10
numbatches=10
bins=50
hopup=1.0
hopdn=1.0

dis_high=8
dis_low=8
dis_step=1.0

coup_high=1.0
coup_low=0.0
coup_step=0.1

couplings=$(seq $coup_low $coup_step $coup_high)
disorders=$(seq $dis_high -$dis_step $dis_low)
batchnums=$(seq 1 1 $numbatches)

for coupling in $couplings; do
    for disorder in $disorders; do
        jobname="mbl_density_error_${size}x${size}_W${disorder}_C${coupling}_TU${hopup}_TD${hopdn}_N${runs}"
        log_filename="logs/${jobname}.log"
        execname="build/keldysh_densities_error"
        /usr/bin/time -v ${execname} -n $runs --batchsize $batchsize --bins $bins -c $coupling -w $disorder -s $size -u $hopup -d $hopdn &> $log_filename
    done
done
