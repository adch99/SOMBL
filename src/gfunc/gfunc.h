#include <complex.h>
#include "../constants.h"

int gfuncsq_std(CDTYPE * eigenvectors, int size,
                    DTYPE * green_func, int degeneracy);
DTYPE gfuncsq_std_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy);
int gfuncsq_restr_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length,
                                int nmin, int nmax);
int gfuncsq_restr_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                                    CDTYPE * eigvecs, int length,
                                    int nmin, int nmax);
int gfuncsq_restr_sigma_op_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length,
                                int nmin, int nmax);
int gfuncsq_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length);
int gfuncsq_sigma_op_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length);
int gfuncsq_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length);

