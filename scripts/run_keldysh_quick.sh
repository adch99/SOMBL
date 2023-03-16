#!/bin/bash

W=15
C=0.2
size=10
bins=10
batchsize=5
numbatches=5
N=25
logfile=
for i in $(seq 1 $numbatches); do
    echo "Batch $i"
    #./build/keldysh_window_batch -s $size -w $W -c $C --batch $i --batchsize $batchsize --bins $bins >> $logfile;
    seed=$(od -vAn -N4 -tu4 < /dev/urandom)
    echo $seed
    ./build/keldysh_window_batch_nobins_error -s $size -n $N -w $W -c $C --batch $i --batchsize $batchsize --bins $bins --seed $seed;
    # sleep 1;
done
#./build/keldysh_energy_batch_average -s $size -w $W -c $C --batchsize $batchsize --bins $bins -a ffffffffff >> $logfile
./build/keldysh_energy_batch_average_error -s $size -n $N -w $W -c $C --batchsize $batchsize -a ffffffffff
# ./build/keldysh_densities -s $size -w $W -c $C --bins $bins >> $logfile
./build/keldysh_densities_error_random_sum -s $size -n $N -w $W -c $C
