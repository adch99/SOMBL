#!/bin/bash

W=15
C=0.2
size=10
bins=10
batchsize=10
numbatches=10
N=100
logfile=logs/output.log
for i in $(seq 1 $numbatches); do
    echo "Batch $i"
    ./build/keldysh_window_batch -s $size -w $W -c $C --batch $i --batchsize $batchsize --bins $bins >> $logfile;
done
./build/keldysh_energy_batch_average -s $size -w $W -c $C --batchsize $batchsize --bins $bins >> $logfile
./build/keldysh_densities -s $size -w $W -c $C --bins $bins >> $logfile
