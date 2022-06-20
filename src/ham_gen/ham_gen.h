#ifndef MBL_HAM_GEN_H
#define MBL_HAM_GEN_H

#include <complex.h>

#define CDTYPE complex double

int hamiltonian(CDTYPE * ham, int len, int width, double coupling_const, double disorder_strength);

#endif /* MBL_HAM_GEN_H*/