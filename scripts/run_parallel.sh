#!/bin/bash

logfile="/dev/null"

size=10
coupling=0.0
hopup=1.0
hopdn=1.0
hopping=1.0
numruns=100
nospin=false
disorder_low=5.0
disorder_high=20.0
disorder_step=1.5

run_sim()
{
    echo "Running -c $1 -w $2"
    build/exact_diag_simulation -c $1 -w $2 \
        -s $size -u $hopup -d $hopdn -n $numruns
}

export -f run_sim

parallel -j$(nproc) run_sim ::: 0 ::: $(seq $disorder_low $disorder_step $disorder_high)

# for coupling in 0 0.1 0.5 1.0 2.0 5.0 10.0; do
# for coupling in 0; do
#     for disorder in $(seq $disorder_low $disorder_step $disorder_high); do
#         # echo "Running -c $coupling -w $disorder"
#     done
# done

