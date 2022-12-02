#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
// #include <lapacke.h>
#include "mkl.h"
#include "diag.h"
#include "../constants.h"


/*
    Returns the eigenvalues of the hermitian matrix given in
    the array of eigvals given. We assume the array has
    already been preprocessed for the used library
    (currently LAPACK). NOTE: The matrix given to
    diagonalize will be destroyed by this function.
*/
int diag_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals)
{    
    // Diagonalize with LAPACK
    // printf("Entered utils_get_eigvalsh\n");
    int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', size,
                            matrix, size, eigvals);
    if (info != 0)
    {
        printf("LAPACKE_zheev error! Code: %d", info);
        return info; // Some error has occured.
    }
    
    return 0;
}

/*
    Returns the eigenvalues of the hermitian matrix given in
    the array of eigvals given. We assume the array has
    already been preprocessed for the used library
    (currently LAPACK). NOTE: The matrix given to
    diagonalize will be destroyed by this function.
    Orthonormalized eigenvectors will be written into the
    matrix.
*/
int diag_get_eigh(CDTYPE * matrix, int size, DTYPE * eigvals)
{
    // printf("Running zheev...");
    // fflush(stdout);
    lapack_int info, lda, n;

    lda = size;
    n = size;

    info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', n,
                        matrix, lda, eigvals);
    if (info != 0)
    {
        printf("LAPACKE_zheev error! Code: %lld", info);
        return(info); // Some error has occured.
    }
    return(0);
}
