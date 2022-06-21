#ifndef MBL_UTILS_H
#define MBL_UTILS_H

#include "../constants.h"

double utils_loc_len(double energy, double * eigenvals, double hop_strength, int len, int eigenfunc_num);
int utils_preprocess_lapack(CDTYPE * matrix, int size, CDTYPE * preprocd);
int utils_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals);

#endif //MBL_UTILS_H