#!/bin/bash

initialcond="pnjunction"
disorders="$(seq 8 1 18)"
couplings="0 $(seq 0.1 0.1 0.9) 1 $(seq 1.1 0.1 1.9) 2 $(seq 2.1 0.1 2.9) 3"
# echo $disorders
# echo $couplings

echo -n "   "
for c in $couplings; do
       printf "%-3.1f " $c
done
echo
for w in $disorders; do
        printf "%-2d " $w
        for c in $couplings; do
                filename="data/mbl_density_100x100_W${w}_C${c}_TU1_TD1_N100_a0_b0_full_${initialcond}.dat"
                if [ -f $filename ]; then
                        echo -n "Y   "
                        # echo -n "(${c},${w}), "
                else
                        echo -n "N   "
                        # echo $filename
                fi
        done
        echo
done

echo "Variance"
echo "--------"
echo -n "   "
for c in $couplings; do
       printf "%-3.1f " $c
done
echo
for w in $disorders; do
        # printf "%-2d " $w
        for c in $couplings; do
                filename="data/mbl_density_100x100_W${w}_C${c}_TU1_TD1_N100_a0_b0_full_${initialcond}.dat.variance"
                if [ -f $filename ]; then
                        # echo -n "Y   "
                        echo -n "(${c},${w}), "
                # else
                        # echo -n "N   "
                        # echo $filename
                fi
        done
        echo
done


