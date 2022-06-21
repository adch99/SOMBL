#include <math.h>
#include <complex.h>
#include <lapacke.h>
#include "../constants.h"

double utils_loc_len(DTYPE energy, DTYPE * eigenvals, DTYPE hop_strength, int len, int eigenfunc_num)
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
        if(i == eigenfunc_num)
            continue;
        eig = *(eigenvals + i);
        lambda += log(fabs(energy - eig)); /* + log(cabsd(*(hop_strengths+i))) if needed*/
    }

    if (skip_one)
        lambda /= len - 1;
    else
        lambda /= len;
    lambda += log(fabs(hop_strength));

    return(1/lambda);    
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
    int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', size, matrix, size, eigvals);
    if (info != 0)
        return info; // Some error has occured.

    
    // Extract and return eigenvalues


    return 0;
}