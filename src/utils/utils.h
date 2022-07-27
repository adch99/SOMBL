#ifndef MBL_UTILS_H
#define MBL_UTILS_H

#include "../constants.h"

#define RTC(i1, i2, size) ((i1) + (i2)*(size)) // row major to col major

DTYPE utils_loc_len(DTYPE energy, const DTYPE * eigenvals, DTYPE hop_strength,
                    int len, int eigenfunc_num);
int utils_preprocess_lapack(CDTYPE * matrix, int size, CDTYPE * preprocd);
int utils_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals);
int utils_row_to_col(int index1, int index2, int size);
int utils_uniform_dist(double low, double high, int num_samples,
                        double * samples, int seed_with_time);
int utils_print_matrix_complex_F(CDTYPE * matrix, int m, int n);
int utils_print_matrix_real_C(DTYPE * matrix, int m, int n);
int utils_print_matrix(void * matrix, int m, int n,
                    char type, char ordering);
int utils_get_eigh(CDTYPE * matrix, int size, DTYPE * eigvals);
int utils_fit_exponential(DTYPE * x, DTYPE * y, int length, DTYPE * exponent,
                        DTYPE * mantissa, DTYPE * residuals);
int utils_get_green_func_lim(CDTYPE * eigenvectors, int size,
                        DTYPE * green_func);
int utils_construct_data_vs_dist(DTYPE * matrix, int size,
                                DTYPE * dists, DTYPE * func);
int utils_compare_datapoints(const void * a, const void * b);



#endif //MBL_UTILS_H