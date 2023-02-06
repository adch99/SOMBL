#ifndef MBL_UTILS_H
#define MBL_UTILS_H

#include "../constants.h"

#define RTC(i1, i2, size) ((i1) + (i2)*(size)) // row major to col major

DTYPE utils_loc_len(DTYPE energy, const DTYPE * eigenvals, DTYPE hop_strength,
                    int len, int eigenfunc_num);
int utils_uniform_dist(DTYPE low, DTYPE high, int num_samples,
                        DTYPE * samples, int seed_with_time);
int utils_print_matrix(void * matrix, int m, int n,
                    char type, char ordering);
int utils_save_matrix(void * matrix, int m, int n,
                    char type, char ordering, FILE * ofile);
int utils_get_green_func_lim(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy);
DTYPE utils_compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy);
DTYPE utils_gfuncsq_sigma(int i, uint alpha, int j, uint beta,
                        CDTYPE * sigma, CDTYPE * eigvecs, int num_states);
int utils_gfuncsq_sigma_matrix(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length);
int utils_gfuncsq_sigma_matrix_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length);
int utils_gfuncsq_sigma_matrix_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length);
int utils_gfuncsq_sigma_matrix_nondeg_restr(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length, int nmin, int nmax);
int utils_gfuncsq_sigma_matrix_deg_restr(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length, int nmin, int nmax);
int utils_gfuncsq_sigma_coeff_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length);
int utils_gfuncsq_sigma_coeff_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length);
int utils_gfuncsq_sigma_coeff(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length);
int utils_multiply_restricted(CDTYPE * mat1, int m1, int n1min, int n1max,
                            CDTYPE * mat2, int m2, int n2min, int n2max,
                            CDTYPE * prod);
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
int utils_bin_energy_range(DTYPE * energies, int len, int num_bins,
                        DTYPE core_min, DTYPE core_max, int * bin_edges);
int utils_add_to_matrix_real(DTYPE * matrix1, DTYPE * matrix2,
                            int m, int n);
int utils_add_to_matrix_complex(CDTYPE * matrix1, CDTYPE * matrix2,
                            int m, int n);


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