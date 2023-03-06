#include <complex.h>
#include "../constants.h"

int gfuncsq_std(CDTYPE * eigenvectors, int size,
                    DTYPE * green_func, int degeneracy);
DTYPE gfuncsq_std_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy);
// int gfuncsq_restr_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
//                                 CDTYPE * eigvecs, int length,
//                                 int nmin, int nmax);
// int gfuncsq_restr_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
//                                     CDTYPE * eigvecs, int length,
//                                     int nmin, int nmax);
// int gfuncsq_restr_sigma_op_deg(DTYPE * gfuncsq, CDTYPE * sigma,
//                                 CDTYPE * eigvecs, int length,
//                                 int nmin, int nmax);
int gfuncsq_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length);
int gfuncsq_sigma_op_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length);
int gfuncsq_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length);
int gfuncsq_sym_GR_GRstar_nondeg(CDTYPE * eigvecs, DTYPE * gfuncsq,
                                int length, int nmin, int nmax,
                                uint alpha);
int gfuncsq_asym_GR_GRstar_nondeg(CDTYPE * eigvecs, CDTYPE * gfuncsq,
                                int length, int nmin, int nmax,
                                uint alpha, uint beta);
int gfuncsq_sym_GR_GRstar_deg(CDTYPE * eigvecs, DTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha);
int gfuncsq_asym_GR_GRstar_deg(CDTYPE * eigvecs, CDTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha, uint beta);
int gfunc_direct_full(CDTYPE * eigvecs, CDTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha, uint beta);


int gfuncsq_sym_GR_GRstar_nondeg_error(CDTYPE * eigvecs, DTYPE * gfuncsq,
                                        DTYPE * sqsum, int length,
                                        int nmin, int nmax, uint alpha);
int gfuncsq_sym_GR_GRstar_deg_error(CDTYPE * eigvecs, DTYPE * gfuncsq,
                                    DTYPE * sqsum, int length, int nmin,
                                    int nmax, uint alpha);
int gfuncsq_asym_GR_GRstar_nondeg_error(CDTYPE * eigvecs, CDTYPE * gfuncsq,
                                CDTYPE * sqsum, int length, int nmin,
                                int nmax, uint alpha, uint beta);
int gfuncsq_asym_GR_GRstar_deg_error(CDTYPE * eigvecs, CDTYPE * gfuncsq,
                            CDTYPE * sqsum, int length, int nmin,
                            int nmax, uint alpha, uint beta);
