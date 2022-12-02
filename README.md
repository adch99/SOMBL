# SOMBL

We have a 2D Anderson model with fermions on a $LxL$ lattice.
The fermion spins are coupled to the motion via spin-orbit
coupling.

## Compiling and Running
```
$ mkdir build
$ mkdir build/{utils,diag,ham_gen,io,params,gfunc,tests,extern}; mkdir build/extern/unity
$ make
```

To check if everything runs correctly do

```
$ make tests
```
