#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <lapacke.h>
#include <gsl/gsl_fit.h>
#include <omp.h>
#include "utils.h"
#include "../constants.h"
 

/*
    Calculates the localization length from the so-called
    "Lyapunov exponent" which is computed from the given
    spectrum. We assume the eigenvalues are sorted in
    ascending order and are real. We also assume the hopping
    strength is constant, however this can be easily
    extended. This can also be done for a complex hopping
    strength even though the phase doesn't actually matter,
    only the absolute values do. If eigenfunc_num is in
    [0,len-1], then we assume it means we are calculating it
    for that eigenfunction, otherwise, we assume the energy
    given is not an eigenvalue.
*/
DTYPE utils_loc_len(DTYPE energy, const DTYPE * eigenvals, DTYPE hop_strength,
                    int len, int eigenfunc_num)
{
    DTYPE lambda = 0;
    DTYPE eig;
    DTYPE loc_len;
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
        lambda += log(fabs(energy - eig));
        /* + log(cabsd(*(hop_strengths+i))) if needed*/
    }

    if (skip_one)
        lambda /= (DTYPE) (len - 1);
    else
        lambda /= (DTYPE) len;
    lambda += -log(fabs(hop_strength));

    loc_len = 1.0 / lambda;
    return(loc_len);    
}

/*
    This function takes a sizexsize matrix and returns the
    column major version of it with whatever other
    preprocessing is needed to get it to work with LAPACK.
    We assume a SQUARE MATRIX and that preprocd has already
    been allocated a size^2 array of CDTYPE.
*/
int utils_preprocess_lapack(CDTYPE * matrix, int size, CDTYPE * preprocd)
{
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

/*
    Returns the eigenvalues of the hermitian matrix given in
    the array of eigvals given. We assume the array has
    already been preprocessed for the used library
    (currently LAPACK). NOTE: The matrix given to
    diagonalize will be destroyed by this function.
*/
int utils_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals)
{

    // Convert first to column major
    // double * colMajorMatrix = malloc(sizeof(complex double)*(size*size));
    

    // Diagonalize with LAPACK
    /*
        If this isn't done by LAPACKE, then we need to call
        zheev twice. First time to get the optimal number of
        eigenvalues that are to be calculated (WORK, LWORK).
    */
    // printf("Entered utils_get_eigvalsh\n");
    int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'N', 'U', size,
                            matrix, size, eigvals);
    if (info != 0)
    {
        printf("LAPACKE_zheev error! Code: %d", info);
        return info; // Some error has occured.
    }
    
    // Extract and return eigenvalues


    return 0;
}

// index1 should be in [0, size-1]
int utils_row_to_col(int index1, int index2, int size)
{
    return(index1 + index2*size);
}

/*
    Generates 'num_samples' numbers in the range [low, high]
    using a uniform distribution. Please ensure that you
*/
int utils_uniform_dist(double low, double high, int num_samples,
                     double * samples, int seed_with_time)
{
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

// Prints a matrix stored in Column Major Form
int utils_print_matrix_complex_F(CDTYPE * matrix, int m, int n)
{
    int i, j;
    CDTYPE element;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            element = *(matrix + RTC(i,j,m));
            printf("%08.4lf+%08.4lfj\t", creal(element),
                    cimag(element));
        }
        printf("\n");
    }
    return 0;
}

int utils_print_matrix_real_C(DTYPE * matrix, int m, int n)
{
    int i, j;
    DTYPE element;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            element = *(matrix + n*i + j);
            printf("%08.4lf\t", element);
        }
        printf("\n");
    }
    return 0;
}

int utils_print_matrix(void * matrix, int m, int n,
                    char type, char ordering)
{
    int i, j, index;
    DTYPE elemr;
    CDTYPE elemc;


    if(type != 'C' && type != 'R')
    {
        printf("Invalid type passed to utils_print_matrix.\n");
        printf("type must be either 'C' or 'R'.\n");
    }


    if(ordering != 'C' && ordering != 'F')
    {
        printf("Invalid ordering passed to utils_print_matrix.\n");
        printf("ordering must be either 'C' or 'F'.\n");
    }


    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = (ordering=='C')?(i*n + j):RTC(i, j, m);
            if(type == 'C')
            {
                elemc = *((CDTYPE*)matrix + index);
                printf("%08.4lf+%08.4lfj\t", creal(elemc),
                        cimag(elemc));
            }
            else
            {
                elemr = *((DTYPE*)matrix + index);
                printf("%08.4lf\t", elemr);
            }
        }
        printf("\n");
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
int utils_get_eigh(CDTYPE * matrix, int size, DTYPE * eigvals)
{
    int info = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'U', size,
                            matrix, size, eigvals);
    if (info != 0)
    {
        printf("LAPACKE_zheev error! Code: %d", info);
        return info; // Some error has occured.
    }
    return 0;
}

/*
    This function takes the eigenvectors, calculates the
    long time limit of the Green's function squared and ADDS
    it to given matrix. If you want only the Green's
    function for the given eigenvectors, then please ensure
    that the green_func array is initialized to zeroes. This
    behaviour helps us calculate the disorder averaged
    Green's function in-place which saves memory. 
*/
int utils_get_green_func_lim(CDTYPE * eigenvectors, int size, DTYPE * green_func)
{

    // Construct green function limit squared
    int i, j;
    // int num_threads;
    // #pragma omp critical
    // {
    //     num_threads = omp_get_num_threads();
    // }
    
    #pragma omp parallel for collapse(2)
    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            DTYPE sum = utils_compute_gfsq_elem(i, j, eigenvectors, size);
            #pragma omp critical
            {
                *(green_func + i*size + j) += sum;
                if(i != j)
                    *(green_func + j*size + i) += sum;
            }

        }
    }    
    // printf("Threads: %d\n", num_threads);
    return 0;
}

/*
    Computes the (i,j)th element of the Green's function
    squared in the long time limit.
*/
DTYPE utils_compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors, int size)
{
    DTYPE sum = 0.0;
    int k;
    #pragma omp parallel for reduction (+:sum) schedule(auto)
    for(k = 0; k < size; k++)
    {
        int index1, index2;
        CDTYPE psi1, psi2;
        DTYPE mod1, mod2;
        index1 = RTC(k, i, size);
        index2 = RTC(k, j, size);
        psi1 = *(eigenvectors + index1);
        psi2 = *(eigenvectors + index2);
        mod1 = creal(psi1)*creal(psi1) + cimag(psi1)*cimag(psi1);
        mod2 = creal(psi2)*creal(psi2) + cimag(psi2)*cimag(psi2);
        sum += mod1 * mod2;
    }

    return(sum);
}


int utils_fit_exponential(DTYPE * x, DTYPE * y, int length, DTYPE * exponent,
                        DTYPE * mantissa, DTYPE * residuals)
{
    DTYPE * logdata = malloc(length*sizeof(DTYPE));
    DTYPE diff;
    int i;
    
    // Take log of the data
    for(i = 0; i < length; i++)
    {
        *(logdata + i) = log(*(y + i));
        if(isnan(*(logdata + i)))
            printf("NaN detected in log output.\n");
    }

    // Construct a weight function for the fitting
    // We use 1 for the middle 80% of the data
    // We use a linearly decaying rate for the ends.

    int start_of_mid = ceil(0.25 * length);
    int end_of_mid = length - ceil(0.25 * length);
    DTYPE * weights = malloc(length * sizeof(DTYPE));
    DTYPE slope = 1.0 / (DTYPE) start_of_mid;
    for(i = start_of_mid; i < end_of_mid; i++)
        *(weights + i) = 1;
    
    for(i = 0; i < start_of_mid; i++)
        *(weights + i) = slope * i;

    for(i = end_of_mid; i < length; i++)
        *(weights + i) = 1 - slope * (i - end_of_mid);

    // Fit line to data using GSL.
    DTYPE c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_wlinear(x, 1, weights, 1, logdata, 1, length,
                &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    if(isnan(c0) || isnan(c1))
        printf("gsl returned NaN.\n");  
    
    *mantissa = exp(c0);
    *exponent = c1;

    if(isnan(*mantissa) || isnan(*exponent))
        printf("mantissa or exponent are NaN.\n");  


    // Calculate residuals
    *residuals = 0;
    for(i = 0; i < length; i++)
    {   
        diff = (*mantissa) * exp((*exponent) * i);
        *residuals += *(weights + i)  * diff*diff/length;
    }

    free(logdata);
    return 0;
}


/*
    Takes the matrix M and calculates the function
    f(r) = average of {M_ij : |i-j|=r}.
    
    This function will initialize `dists` and `func` to the
    appropriate lengths. Do not initialize these beforehand.
    And remember to free them after use.
*/
int utils_construct_data_vs_dist(DTYPE * matrix, int size, int length,
                            int bins, DTYPE ** dists, DTYPE ** func)
{
    int i, j;
    int nospin = (size==length*length)?1:0;

    unsigned int spin1down, spin2down;
    spin1down = spin2down = 0;
    BIT_SET(spin1down, 1);
    BIT_SET(spin2down, 0);

    DTYPE lowest = 0;
    DTYPE highest = length / sqrt(2);
    DTYPE bin_width = (highest - lowest) / (DTYPE) bins;
    int * counts = calloc(bins, sizeof(int));
    *dists = calloc(bins, sizeof(DTYPE));
    *func = calloc(bins, sizeof(DTYPE));

    // printf("bin_width = %e\n", bin_width);

    // Initialize the x values
    // to the midpoints of the bins
    for(i = 0; i < bins; i++)
    {
        *(*dists + i) = ((DTYPE) i + 0.5) * bin_width; 
        // printf("i=%d dist=%e\n", i, *(*dists + i));
    }

    // printf("size = %d\tlength = %d\n", size, length);

    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            int site1x, site2x, site1y, site2y;
            int dx, dy;
            unsigned int spin_index1, spin_index2;
            unsigned int spins;
            DTYPE len, value;

            utils_get_lattice_index(i, length, nospin,
                            &site1x, &site1y, &spin_index1);
            utils_get_lattice_index(j, length, nospin,
                            &site2x, &site2y, &spin_index2);
 
            dx = abs(site1x - site2x);
            dy = abs(site1y - site2y);
            // printf("i = %d\tj = %d\n", i, j);
            // printf("x1 = %d\ty1 = %d\tx2 = %d\ty2 = %d\n", site1x, site1y, site2x, site2y);
            // printf("dx = %d\tdy = %d\n", dx, dy);
            dx = INTMIN(abs(length/2-dx), dx);
            dy = INTMIN(abs(length/2-dy), dy);
            len = sqrt((DTYPE) (dx*dx + dy*dy));
            // printf("dx = %d\tdy = %d\tlen = %e\n", dx, dy, len);
            spins = spin_index1*spin1down + spin_index2*spin2down;

            if(i == j)
            {
                value = *(matrix + RTC(i, j, size));
                utils_bin_data(len, value, bins, counts,
                                lowest, bin_width, *func);
            }
            else
            {
                value = *(matrix + RTC(i, j, size));
                utils_bin_data(len, value, bins, counts,
                                lowest, bin_width, *func);
                value = *(matrix + RTC(j, i, size));
                utils_bin_data(len, value, bins, counts,
                                lowest, bin_width, *func);
            }
        }
    }

    // Calculate the averages
    for(i = 0; i < bins; i++)
    {
        if(isnan(*(*func + i)))
        {
            printf("NaN detected at site i = %d\n", i);
        }

        if(*(counts + i) == 0)
        {
            *(*func + i) = DBL_MIN;
            // printf("bin i=%d has counts=%d\n", i, *(counts + i));
        }
        else
            *(*func + i) /= (DTYPE) *(counts + i);
    }

    return(bins);
}

/*
    Given the matrix index `index`, returns the (x,y)
    position index on the lattice as well the appropriate
    spin index. For nospin=1, spin is set to 0. 
*/
int utils_get_lattice_index(int index, int length, int nospin,
                        int * x, int * y, unsigned int * spin)
{
    int site_index;
    if(nospin == 1)
    {
        site_index = index;
        *spin = 0;
    }
    else
    {
        site_index = index / 2;
        *spin = (unsigned int) (index % 2);
    }
    
    *x = site_index % length;
    *y = site_index / length;

    return(0);
}

int utils_bin_data(DTYPE index, DTYPE value, int bins, int * counts,
                DTYPE lowest, DTYPE bin_width, DTYPE * value_hist)
{
    /*
        Each data point is an index and a value.
        The index determines which bin the data
        point goes into. The value is added to
        the sum corresponding to that bin.
        At the end, this will be averaged out. 
    */

    if(isnan(value) || isnan(index))
    {
        printf("NaN value passed\nvalue = %e\nindex = %e\n", value, index);
        return(-1);
    }

    int bin_num = (int) floor((index - lowest) / bin_width);

    if (bin_num >= bins || bin_num < 0)
    {
        printf("Value %e doesn't lie in range: 0 - %e\n", index, bin_width*bins);
        return(-1);
    }

    *(counts + bin_num) += 1;
    *(value_hist + bin_num) += value;

    return(0);
}


DTYPE utils_pbc_chord_length(int index1, int length1, int index2, int length2)
{
    DTYPE radius1 = (DTYPE) length1 / (2.0 * M_PI);
    DTYPE radius2 = (DTYPE) length2 / (2.0 * M_PI);

    DTYPE angle1 = 2.0 * M_PI * ((DTYPE) index1 / (DTYPE) length1);
    DTYPE angle2 = 2.0 * M_PI * ((DTYPE) index2 / (DTYPE) length2);



    return 0;   
}