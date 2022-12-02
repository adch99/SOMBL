#!/bin/bash

size=40
runs=100
batchsize=15
numbatches=6
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
batchnums=$(seq 1 1 $numbatches)

for coupling in $couplings; do
    for disorder in $disorders; do
        for batchnum in $batchnums; do
            jobfilename="jobs/mbl_${size}x${size}_W${disorder}_C${coupling}_TU${hopup}_TD${hopdn}_N${runs}_BS${batchsize}_B${batchnum}.pbs"
            qsub $jobfilename
        done
    done
done