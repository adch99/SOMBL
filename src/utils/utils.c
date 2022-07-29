#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <lapacke.h>
#include <gsl/gsl_fit.h>
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
    int i, j, k;
    int index1, index2;
    DTYPE value;
    CDTYPE psi1, psi2;
    DTYPE mod1, mod2;
    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            for(k = 0; k < size; k++)
            {
                index1 = RTC(k, i, size);
                index2 = RTC(k, j, size);
                // value = cabs(*(eigenvectors + index1)
                //         * conj(*(eigenvectors + index2)));
                psi1 = *(eigenvectors + index1);
                psi2 = *(eigenvectors + index2);
                mod1 = creal(psi1)*creal(psi1) + cimag(psi1)*cimag(psi1);
                mod2 = creal(psi2)*creal(psi2) + cimag(psi2)*cimag(psi2);
                value = mod1 * mod2;
                *(green_func + i*size + j) += value;
                if(i != j)
                    *(green_func + j*size + i) += value;
            }
        }
    }    

    return 0;
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
    }
    // Fit line to data using GSL or with our own formulas.
    DTYPE c0, c1, cov00, cov01, cov11, sumsq;
    gsl_fit_linear(x, 1, logdata, 1, length, &c0, &c1, &cov00,
                    &cov01, &cov11, &sumsq);
    *mantissa = exp(c0);
    *exponent = c1;

    // Calculate residuals
    *residuals = 0;
    for(i = 0; i < length; i++)
    {   
        diff = (*mantissa) * exp((*exponent) * i);
        *residuals += diff*diff/length;
    }

    free(logdata);
    return 0;
}



struct FuncDataPoint
{
    DTYPE func_val;
    int dist; // Keeping ints allows us to compare easily
    int counter;
    // spins: up,up = 0b00, up,down=0b01, down,up=0b10=2
    // down,down=0b11=3
    unsigned int spins; 
};

/*
    Takes the matrix M and calculates the function
    f(r) = average of {M_ij : |i-j|=r}.
    
    This function will initialize `dists` and `func` to the
    appropriate lengths. Do not initialize these beforehand.
    And remember to free them after use.
*/
int utils_construct_data_vs_dist(DTYPE * matrix, int size, int length,
                                DTYPE ** dists, DTYPE ** func)
{
    // Construct a list of all possible lengths
    // Max possible lengths is (L/2)^2 choose 2, so we'll
    // use that as our baseline for array length.

    int i, j, k, total_lens, match;
    int len, dx, dy;
    int nospin = (size==length*length)?0:1;

    // int max_poss_lens = (size*size)*((size*size) - 4) / 32;

    // Let D(n) = number of distinct distances on a lattice
    // (w/o PBCs) of side n. We computed D(n) numerically
    // for n=0-100 and found D(n) ~ 0.33n^2 + 3.38n - 15.22
    // D(n) ~ n^1.83 * 10^(-0.106) For lattice with PBC
    // D_pbc(n) ≈ D(n/2) and also by the way we define our
    // lattice, length = n + 1. So we use a safe upper bound
    // by doing the following:    
    int max_poss_lens;
    if(nospin)
        max_poss_lens = 0.4*length*length + 5*length - 14;
    else
        max_poss_lens = 0.4*length*length + 5*length - 14;

    if(max_poss_lens < 32)
        max_poss_lens = 32; 
    match = 0;
    total_lens = 0;
    struct FuncDataPoint * data;
    struct FuncDataPoint * datapoint;
    data = calloc(max_poss_lens, sizeof(struct FuncDataPoint));

    // While the property is named dist, within this
    // function it will actually store dist^2 values
    // because for the sake of comparison, it is far better
    // to avoid the extra cost of calculating square roots.
    // At the end, however, we will calculate the square
    // roots and the passed value is of the distances.


    int site_index1, site_index2;
    unsigned int spin_index1, spin_index2;
    unsigned int spins;
    unsigned int spin1down, spin2down;
    spin1down = spin2down = 0;
    BIT_SET(spin1down,1);
    BIT_SET(spin2down, 0);
    int site1x, site2x, site1y, site2y;

    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            // I know a linear search is not the best, but
            // the alternative is to sort at each stage and
            // I really don't want to write an insertion
            // sort.

            if(nospin)
            {
                site_index1 = i;
                site_index2 = j;
                spin_index1 = spin_index2 = 0;
            }
            else
            {
                site_index1 = i / 2;
                site_index2 = j / 2;
                spin_index1 = i % 2;
                spin_index2 = j % 2;
            }

            site1x = site_index1 % length;
            site1y = site_index1 / length;
            site2x = site_index2 % length;
            site2y = site_index2 / length;
 

            dx = abs(site1x - site2x);
            dy = abs(site1y - site2y);
            dx = INTMIN(abs(length/2-dx), dx);
            dy = INTMIN(abs(length/2-dy), dy);
            len = dx*dx + dy*dy;
            spins = spin_index1*spin1down + spin_index2*spin2down;
            match = 0;

            for(k = 0; k < total_lens; k++)
            {
                datapoint = data + k;
                if((datapoint->dist == len) && (datapoint->spins == spins))
                {   
                    match = 1;
                    break;
                }
            }
            if(match == 0)
            {
                datapoint = data + total_lens; // This is redundant but
                // present for clarity
                datapoint->dist = len;
                datapoint->spins = spins;
                datapoint->counter = 0;
                datapoint->func_val = 0;
                total_lens += 1;
            }

            if(i == j)
            {
                datapoint->counter += 1;
                datapoint->func_val += *(matrix + RTC(i,j,size));
            }
            else
            {
                datapoint->counter += 2;
                datapoint->func_val += *(matrix + RTC(i,j,size));
                datapoint->func_val += *(matrix + RTC(j,i,size));
            }
        }
    }

    // Calculate the averages
    for(i = 0; i < total_lens; i++)
    {
        datapoint = data + i;
        datapoint->func_val /= (DTYPE) datapoint->counter;
    }

    // Sort the array
    qsort(data, total_lens, sizeof(struct FuncDataPoint), utils_compare_datapoints);

    *dists = malloc(total_lens * sizeof(DTYPE));
    *func = malloc(total_lens * sizeof(DTYPE));
    
    for(i = 0; i < total_lens; i++)
    {
        datapoint = data + i;
        *(*dists + i) = sqrt((DTYPE) datapoint->dist);
        *(*func + i) = datapoint->func_val;
    }
    
    free(data);

    return(total_lens);
}

int utils_compare_datapoints(const void * a, const void * b)
{
  const struct FuncDataPoint * da = (const struct FuncDataPoint *) a;
  const struct FuncDataPoint * db = (const struct FuncDataPoint *) b;
  return (da->dist > db->dist) - (da->dist < db->dist);
}
