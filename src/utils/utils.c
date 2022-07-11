#include "utils.h"
#include <math.h>
#include <complex.h>
#include <lapacke.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "../constants.h"
 

DTYPE utils_loc_len(DTYPE energy, const DTYPE * eigenvals, DTYPE hop_strength, int len, int eigenfunc_num)
{
    /*
    Calculates the localization length from the so-called "Lyapunov exponent"
    which is computed from the given spectrum. We assume the eigenvalues are
    sorted in ascending order and are real.
    We also assume the hopping strength is constant, however this can be easily
    extended. This can also be done for a complex hopping strength even though
    the phase doesn't actually matter, only the absolute values do.
    If eigenfunc_num is in [0,len-1], then we assume it means we are calculating it for
    that eigenfunction, otherwise, we assume the energy given is not an eigenvalue.
    */

    DTYPE lambda = 0;
    DTYPE eig;
    int i;
    int skip_one;
    skip_one = (eigenfunc_num >=0 && eigenfunc_num < len);
    if (skip_one)
        energy = *(eigenvals + eigenfunc_num);

    for(i=0;i<len;i++)
    {
        if(skip_one && i == eigenfunc_num)
            continue;
        eig = *(eigenvals + i);
        lambda += log(fabs(energy - eig)); /* + log(cabsd(*(hop_strengths+i))) if needed*/
    }

    if (skip_one)
        lambda /= (DTYPE) (len - 1);
    else
        lambda /= (DTYPE) len;
    lambda += -log(fabs(hop_strength));

    return(1.0 / lambda);    
}

int utils_preprocess_lapack(CDTYPE * matrix, int size, CDTYPE * preprocd)
{
    /*
        This function takes a sizexsize matrix and returns the column major
        version of it with whatever other preprocessing is needed to get it to
        work with LAPACK.
        We assume a SQUARE MATRIX and that preprocd has already been allocated
        a size^2 array of CDTYPE.
    */
    int i, j;

    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            *(preprocd + i + size*j) = *(matrix + i*size + j);
        }
    }

    return(0);
}

int utils_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals)
{
    /*
    Returns the eigenvalues of the hermitian matrix given in the
    array of eigvals given.
    We assume the array has already been preprocessed for the used
    library (currently LAPACK).
    NOTE: The matrix given to diagonalize will be destroyed by this
    function.
    */

    // Convert first to column major
    // double * colMajorMatrix = malloc(sizeof(complex double)*(size*size));
    

    // Diagonalize with LAPACK
    /*
        If this isn't done by LAPACKE, then we need to call zheev twice.
        First time to get the optimal number of eigenvalues that are to
        be calculated (WORK, LWORK).
    */
    // printf("Entered utils_get_eigvalsh\n");
    int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', size, matrix, size, eigvals);
    if (info != 0)
    {
        printf("LAPACKE_zheev error! Code: %d", info);
        return info; // Some error has occured.
    }
    
    // Extract and return eigenvalues


    return 0;
}

int utils_row_to_col(int index1, int index2, int size)
{
    // index1 should be in [0, size-1]
    return(index1 + index2*size);
}

int utils_uniform_dist(double low, double high, int num_samples, double * samples,
                        int seed_with_time)
{
    /*
        Generates 'num_samples' numbers in the range [low, high]
        using a uniform distribution. Please ensure that you
    */
    int randint, i;
    double div, scale = (high - low);
    if (seed_with_time)
        srandom((unsigned) time(NULL));
    
    for(i = 0; i < num_samples; i++)
    {
        randint = random();
        div = ((double) randint) / ((double) RAND_MAX);
        *(samples + i) = low + scale * div;
    }
    return 0;
}

int utils_print_matrix(CDTYPE * matrix, int m, int n)
{
    // Prints a matrix stored in Column Major Form
    int i, j;
    CDTYPE element;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            element = *(matrix + RTC(i,j,m));
            printf("%08.4lf+%08.4lfj\t", creal(element), cimag(element));
        }
        printf("\n");
    }
    return 0;
}