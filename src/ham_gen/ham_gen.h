#ifndef MBL_HAM_GEN_H
#define MBL_HAM_GEN_H

#include "../constants.h"

#define NEIGHS 4

int hamiltonian(CDTYPE * ham, int len, int width,
                DTYPE coupling_const, DTYPE disorder_strength,
                DTYPE hop_strength, int (* neighbours)[NEIGHS]);
int hamiltonian_nospin(CDTYPE * ham, int len, int width,
                DTYPE disorder_strength, DTYPE hop_strength,
                int (*neighbours)[NEIGHS]);
int get_neighbour_lists(int (*neighbours)[NEIGHS], int len, int width);
int get_neighbours(int index, int len, int width, int * nlist);
int check_neighbour(int index, int * nlist);

#endif /* MBL_HAM_GEN_H*/