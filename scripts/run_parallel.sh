#!/bin/bash

# logfile="/dev/null"

size=40
runs=100
hopup=1.0
hopdn=1.0
hopping=1.0
# nospin=false

coup_low=2.0
coup_high=20.0
coup_step=2.0
couplings=($(seq $coup_low $coup_step $coup_high))

disorder=18.5
# dis_low=6.5
# dis_high=18.5
# dis_step=-1.5
# disorders=($(seq $dis_low $dis_step $dis_high))

# If we have a 2D space of params
# declare -a coup_array
# declare -a dis_array
# i=0
# for c in ${couplings[@]}; do
#     for w in ${disorders[@]}; do
#         coup_array[i]=$c
#         dis_array[i]=$w
#         i++
#     done
# done

# declare -p coup_array
# declare -p dis_array

# echo ${coup_array[*]}
# echo ${dis_array[*]}

coup_array=(${couplings[*]})
# dis_array=(${disorders[*]})

total_samples=10

run_basic()
{
    echo "Running -c $1 -w $2"
    build/exact_diag_simulation -c $1 -w $2 \
        -s $size -u $hopup -d $hopdn -n $numruns
}

export -f run_basic


nproc=4
batches=$((total_samples / nproc))
leftover=$((total_samples % nproc))

echo $batches

for (( batch=1; batch<=batches; batch++ ))
do
    echo "Batch $batch starting"
    echo $(date --iso-8601=seconds)
    for (( proc=1; proc<=$nproc; proc++ ))
    do
        num=$(($nproc*($batch - 1) + ($proc - 1)))
        coupling=${coup_array[num]}
        # disorder=${dis_array[num]}
        # run_basic $coupling $disorder
        simple "-c $coupling -w $disorder" 5 &
    done
    wait
done

# Leftover processes
batch=$batches
echo "Batch $batch starting"
echo $(date --iso-8601=seconds)
for (( proc=1; proc<=$leftover; proc++ ))
do
    num=$(($nproc*($batch - 1) + ($proc - 1)))
    coupling=${coup_array[num]}
    # disorder=${dis_array[num]}
    # run_basic $coupling $disorder
    simple "-c $coupling -w $disorder" 5 &
done
wait



# basename="run_diag_params"
# srcbase="src/${basename}.c"
# for coupling in 0 0.1 0.5 1; do
#     srcname="src/${basename}_c${coupling}.c"
#     execname="build/${basename}_c${coupling}"
#     sed -e "s/#SIZE#/${size}/g" $srcbase > src/tmp1.c
#     sed -e "s/#COUP#/${coupling}/g" src/tmp1.c > src/tmp2.c
#     sed -e "s/#RUNS#/${runs}/g" src/tmp2.c > $srcname
#     make $execname
#     rm src/tmp{1,2}.c
#     rm $srcname
#     echo "c=${coupling} s=${size} n=${runs}"
#     # echo "$execname"
#     nohup $execname &
# done

# parallel -j$(nproc) run_basic ::: 0 ::: $(seq $disorder_low $disorder_step $disorder_high)

# for coupling in 0 0.1 0.5 1.0 2.0 5.0 10.0; do
# for coupling in 0; do
#     for disorder in $(seq $disorder_low $disorder_step $disorder_high); do
#         # echo "Running -c $coupling -w $disorder"
#     done
# done
