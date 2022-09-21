#ifndef MBL_HAM_GEN_H
#define MBL_HAM_GEN_H

#include "../constants.h"

int hamiltonian(CDTYPE * ham, int len, int width,
                DTYPE coupling_const, DTYPE disorder_strength,
                DTYPE hop_strength_upup, DTYPE hop_strength_dndn,
                int (*neighbours)[NEIGHS]);
int hamiltonian_nospin(CDTYPE * ham, int len, int width,
                DTYPE disorder_strength, DTYPE hop_strength,
                int (*neighbours)[NEIGHS]);
int get_neighbour_lists(int (*neighbours)[NEIGHS], int len, int width);
int get_neighbours_edges(int (*neighbours)[NEIGHS], int len, int width);
int get_neighbours_corners(int (*neighbours)[NEIGHS], int len, int width);
int get_neighbours_bulk(int index, int len, int * nlist);
int check_neighbour(int index, int * nlist);

#define LATIDX(x,y,len) ((x) + (len)*(y))

#endif /* MBL_HAM_GEN_H*/