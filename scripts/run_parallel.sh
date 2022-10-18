#!/bin/bash

# logfile="/dev/null"

size=40
# coupling=0.0
runs=100
# hopup=1.0
# hopdn=1.0
# hopping=1.0
# nospin=false
# disorder_low=6.5
# disorder_high=18.5
# disorder_step=-1.5

basename="run_diag_params"
srcbase="src/${basename}.c"
for coupling in 0 0.1 0.5 1; do
    srcname="src/${basename}_c${coupling}.c"
    execname="build/${basename}_c${coupling}"
    sed -e "s/#SIZE#/${size}/g" $srcbase > src/tmp1.c
    sed -e "s/#COUP#/${coupling}/g" src/tmp1.c > src/tmp2.c
    sed -e "s/#RUNS#/${runs}/g" src/tmp2.c > $srcname
    make $execname
    rm src/tmp{1,2}.c
    rm $srcname
    echo "c=${coupling} s=${size} n=${runs}"
    # echo "$execname"
    nohup $execname &
done

# run_basic()
# {
#     echo "Running -c $1 -w $2"
#     build/exact_diag_simulation -c $1 -w $2 \
#         -s $size -u $hopup -d $hopdn -n $numruns
# }

# export -f run_basic

# parallel -j$(nproc) run_basic ::: 0 ::: $(seq $disorder_low $disorder_step $disorder_high)

# for coupling in 0 0.1 0.5 1.0 2.0 5.0 10.0; do
# for coupling in 0; do
#     for disorder in $(seq $disorder_low $disorder_step $disorder_high); do
#         # echo "Running -c $coupling -w $disorder"
#     done
# done
