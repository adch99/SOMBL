#!/bin/bash

# couplings=(2.0 5.0 7.5 10.0 12.5 15.0 17.5 20.0)
couplings=($(seq 2.0 2.0 20.0))
echo ${couplings[3]}

for i in {0..1}; do
    num=$((4*i))
    simple "$num -c ${couplings[num]}" 5 &
    num=$((4*i + 1))
    simple "$num -c ${couplings[num]}" 5 &
    num=$((4*i + 2))
    simple "$num -c ${couplings[num]}" 5 &
    num=$((4*i + 3))
    simple "$num -c ${couplings[num]}" 5 &
    wait   
done