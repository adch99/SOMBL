#!/bin/bash

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

chunksize=20
numchunks=$((runs / chunksize))
if [ "$((runs % chunksize))" != 0 ]; then
    echo "We need the number of runs to be a multiple of the chunksize!"
    echo "$runs is not a multiple of $chunksize"
    exit 1
fi



# nproc=4
# procnum=1
for coupling in $couplings; do
    for disorder in $disorders; do
        jobfilename="jobs/mbl_${size}x${size}_W${disorder}_C${coupling}_TU${hopup}_TD${hopdn}_N${runs}.pbs"
        qsub $jobfilename
        # ./build/exact_diag_simulation -n $runs -c $coupling -w $disorder -s $size -u $hopup -d $hopdn &
        # rem=$((procnum % nproc))
        # if [ "$rem" = "0" ]; then
        #     wait
        # fi
        # procnum=$((procnum+1))
    done
done
wait


function run_chunked
# {
#     # $1 = number of runs
#     # $2 = chunk size
#     # $3 = coupling strength
#     # $4 = disorder
#     numchunks=$(($1 / $2))
#     if [ "$(($1 % $2))" != 0 ]; then
#         echo "We need the chunksize to be a multiple of number of runs!"
#         echo "$1 is not a multiple of $2"
#         exit 1
#     fi
#     for chunkidx in $(seq 1 $numchunks); do
#         simple "./build/exact_diag_chunk -k $chunkidx -n $2 -c $3 -w $4 -s $size -u $hopup -d $hopdn" 5 &
#     done
# }