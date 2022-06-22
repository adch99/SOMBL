#ifndef MBL_UTILS_H
#define MBL_UTILS_H

#include "../constants.h"

#define RTC(i1, i2, size) (i1 + i2*size) // row major to col major

DTYPE utils_loc_len(DTYPE energy, DTYPE * eigenvals, DTYPE hop_strength, int len, int eigenfunc_num);
int utils_preprocess_lapack(CDTYPE * matrix, int size, CDTYPE * preprocd);
int utils_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals);
int utils_row_to_col(int index1, int index2, int size);
int utils_uniform_dist(double low, double high, int num_samples, double * samples);

#endif //MBL_UTILS_H