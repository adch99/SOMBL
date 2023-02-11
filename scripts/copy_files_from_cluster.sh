#!bash

w=15.0
c=0.2
batch=1
batchsize=10
n=100

path=""
dir=""
fname="data/mbl_100x100_W$w_C$c_TU1_TD1_N$n_BS$batchsize_B$batch_grgrstar_a0_b0_full.dat"
fullname="$path/$dir/data/$name"
rsync mursalin@flock:$fullname $fname

fname="data/mbl_100x100_W$w_C$c_TU1_TD1_N$n_BS$batchsize_B$batch_grgrstar_a1_b1_full.dat"
fullname="$path/$dir/data/$name"
rsync mursalin@flock:$fullname $fname

fname="data/mbl_100x100_W$w_C$c_TU1_TD1_N$n_BS$batchsize_B$batch_grgrstar_a0_b1_full.dat"
fullname="$path/$dir/data/$name"
rsync mursalin@flock:$fullname $fname
