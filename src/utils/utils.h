#ifndef MBL_UTILS_H
#define MBL_UTILS_H

#include "../constants.h"

#define RTC(i1, i2, size) ((i1) + (i2)*(size)) // row major to col major

DTYPE utils_loc_len(DTYPE energy, const DTYPE * eigenvals, DTYPE hop_strength,
                    int len, int eigenfunc_num);
int utils_uniform_dist(double low, double high, int num_samples,
                        double * samples, int seed_with_time);
int utils_print_matrix(void * matrix, int m, int n,
                    char type, char ordering);
int utils_save_matrix(void * matrix, int m, int n,
                    char type, char ordering, FILE * ofile);
int utils_get_green_func_lim(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy);
DTYPE utils_compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy);
int utils_get_lattice_index(int index, int length, int nospin,
                        int * x, int * y, unsigned int * spin);
int utils_get_matrix_index(int x, int y, unsigned int spin,
                        int length, int nospin);
int utils_construct_data_vs_dist(DTYPE * matrix, int size, int length,
                            int bins, DTYPE ** dists, DTYPE ** func,
                            DTYPE ** funcerr);
int utils_bin_data(DTYPE index, DTYPE value, int bins, int * counts,
                DTYPE lowest, DTYPE bin_width, DTYPE * value_hist,
                DTYPE * errors);
DTYPE utils_get_charge_imbalance(DTYPE * gfuncsq, int * occupied_set_up,
                        int set_length_up, int * occupied_set_dn,
                        int set_length_dn, int num_states);
DTYPE utils_pbc_chord_length_sq(int index1, int length1,
                                int index2, int length2);
int utils_reflect_upper_to_lower(DTYPE * matrix, int n);

// Bitwise manipulations
// From SO community wiki post
/* a=target variable, b=bit number to act upon 0-n */
#define BIT_SET(a,b) ((a) |= (1ULL<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1ULL<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1ULL<<(b)))
#define BIT_CHECK(a,b) (!!((a) & (1ULL<<(b))))        // '!!' to make sure this returns 0 or 1

#define BITMASK_SET(x, mask) ((x) |= (mask))
#define BITMASK_CLEAR(x, mask) ((x) &= (~(mask)))
#define BITMASK_FLIP(x, mask) ((x) ^= (mask))
#define BITMASK_CHECK_ALL(x, mask) (!(~(x) & (mask)))
#define BITMASK_CHECK_ANY(x, mask) ((x) & (mask))

// Can't believe this isn't a std library function
// Do NOT use these for floats/doubles
#define INTMIN(a,b) (((a) > (b))?(b):(a))
#define INTMAX(a,b) (((a) < (b))?(b):(a))

// Flags for convenience
#define DEGEN_EIGVALS 1
#define NONDEGEN_EIGVALS 0


#endif //MBL_UTILS_H