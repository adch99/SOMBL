#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
// #include "mkl.h"
// #include <cblas.h>
// #include <omp.h>
#include "utils.h"
#include "../constants.h"

#define TOL 1e-6
#define THRESHOLD TOL

enum Spin {UP=0, DOWN};


/*
    Generates 'num_samples' numbers in the range [low, high]
    using a uniform distribution. Please ensure that you
    give a reasonable number for low and high i.e.
    low < high.
*/
int utils_uniform_dist(DTYPE low, DTYPE high, int num_samples,
                    DTYPE * samples, int seed_with_time)
{
    int randint, i;
    DTYPE div, scale = (high - low);
    if (seed_with_time)
        srandom((unsigned) time(NULL));
    
    for(i = 0; i < num_samples; i++)
    {
        randint = random();
        div = ((DTYPE) randint) / ((DTYPE) RAND_MAX);
        *(samples + i) = low + scale * div;
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


    if(ordering != 'R' && ordering != 'C')
    {
        printf("Invalid ordering passed to utils_print_matrix.\n");
        printf("ordering must be either 'C' or 'R'.\n");
    }


    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = (ordering=='R')?(i*n + j):RTC(i, j, m);
            if(type == 'C')
            {
                elemc = *((CDTYPE*)matrix + index);
                printf("%10.4lg%+11.4lgj ", crealf(elemc),
                        cimagf(elemc));
            }
            else
            {
                elemr = *((DTYPE*)matrix + index);
                printf("%10.4lg\t", elemr);
            }
        }
        printf("\n");
    }
    return 0;
}






/*
    Divides up the energy values into `num_bins` bins.
    `numbins-2` bins are in [core_min, core_max)
    The other two bins will be on either side of the core.
    `bin_edges` contains the starting index of the energy
    values in that bin.
*/
int utils_bin_energy_range(DTYPE * energies, int len, int num_bins,
                        DTYPE core_min, DTYPE core_max, int * bin_edges)
{
    if(num_bins <= 2)
    {
        printf("More than 2 bins needed to function correctly!\n");
        return(-1);
    }

    DTYPE bin_width = (core_max - core_min) / (num_bins - 2.0);

    *(bin_edges + 0) = 0;
    int i;
    int bin_ctr = 1;
    DTYPE cur_bin_max_energy;
    for(i = 0; i < len; i++)
    {
        cur_bin_max_energy = core_min + (bin_ctr - 1)*bin_width;
        if(*(energies + i) >= cur_bin_max_energy)
        {
            *(bin_edges + bin_ctr) = i;
            bin_ctr++;
            if(bin_ctr == num_bins)
                break;
        }     
    }

    return(0);
}




/*
    Takes the matrix M and calculates the function
    f(r) = average of {M_ij : |i-j|=r}.
    
    This function will initialize `dists` and `func` to the
    appropriate lengths. Do not initialize these beforehand.
    And remember to free them after use.
*/
int utils_construct_data_vs_dist(DTYPE * matrix, int size, int length,
                            int bins, DTYPE ** dists, DTYPE ** func,
                            DTYPE ** funcerr)
{
    int i, j;
    int nospin = (size==length*length)?1:0;

    unsigned int spin1down, spin2down;
    spin1down = spin2down = 0;
    BIT_SET(spin1down, 1);
    BIT_SET(spin2down, 0);

    DTYPE lowest = 0;

    // For our purposes, we will ignore the
    // range outside length/2. Change the
    // line below to include the entire range
    DTYPE highest = length;
    // DTYPE highest = (DTYPE) length / 2.0;
    // DTYPE highest = (DTYPE) length * sqrt(2.0) / M_PI;
    int * counts;
    DTYPE bin_width = (highest - lowest) / (DTYPE) bins;
    if(nospin == 1)
    {
        counts = calloc(bins, sizeof(int));
        *dists = calloc(bins, sizeof(DTYPE));
        *func = calloc(bins, sizeof(DTYPE));
        *funcerr = calloc(bins, sizeof(DTYPE));
    }
    else
    {
        *dists = calloc(bins, sizeof(DTYPE));
        counts = calloc(4*bins, sizeof(int));
        *func = calloc(4*bins, sizeof(DTYPE));
        *funcerr = calloc(4*bins, sizeof(DTYPE));
    }

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
            // dx = INTMIN(abs(length/2-dx), dx);
            // dy = INTMIN(abs(length/2-dy), dy);
            len = sqrt((DTYPE) (dx*dx + dy*dy));
            // len = sqrt(utils_pbc_chord_length_sq(dx, length, dy, length));
            // printf("dx = %d\tdy = %d\tlen = %e\n", dx, dy, len);
            spins = spin_index1*spin1down + spin_index2*spin2down;

            if(i == j)
            {
                value = *(matrix + RTC(i, j, size));
                utils_bin_data(len, value, bins, counts + bins*spins,
                                lowest, bin_width, *func + bins*spins,
                                *funcerr + bins*spins);
            }
            else
            {
                value = *(matrix + RTC(i, j, size));
                utils_bin_data(len, value, bins, counts + bins*spins,
                                lowest, bin_width, *func + bins*spins,
                                *funcerr + bins*spins);
                value = *(matrix + RTC(j, i, size));
                utils_bin_data(len, value, bins, counts + bins*spins,
                                lowest, bin_width, *func + bins*spins,
                                *funcerr + bins*spins);
            }
        }
    }

    // Calculate the averages
    unsigned int spins;
    unsigned int total_spins;
    int index;
    if(nospin == 1)
        total_spins = 1;
    else
        total_spins = 4;

    for(spins = 0; spins < total_spins; spins++)
    {
        for(i = 0; i < bins; i++)
        {
            index = spins*bins + i;
            if(isnan(*(*func + index)))
            {
                printf("NaN detected at site i = %d\n", i);
            }

            if(*(counts + index) == 0)
            {
                *(*func + index) = DBL_MIN;
                // printf("bin i=%d has counts=%d\n", i, *(counts + i));
            }
            else
            {
                // Divide by the counts
                *(*func + index) /= (DTYPE) *(counts + index);
                // Near zero values can be negative
                // Make them positive for plotting
                // and convenience.
                DTYPE value = *(*func + index);
                if(value < 0 && fabs(value) < 1e-15)
                    *(*func + index) *= -1;

                // Divide by counts to get <x^2>
                // Then calculate sqrt(<x^2> - <x>^2)
                *(*funcerr + index) /= (DTYPE) *(counts + index);
                *(*funcerr + index) -= value*value;
                *(*funcerr + index) = sqrt(*(*funcerr + index));
            }
        }
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

/*
    Given the lattice position (x, y) and spin,
    returns the index in the gfunc/ham matrix.
*/
int utils_get_matrix_index(int x, int y, unsigned int spin,
                        int length, int nospin)
{
    int pbc_x, pbc_y;
    if(x < 0)
        pbc_x = length - (abs(x) % length);
    else
        pbc_x = x % length;
    
    if(y < 0)
        pbc_y = length - (abs(y) % length);
    else
        pbc_y = y % length;

    int lattice_index = pbc_x + length*pbc_y;
    int index;
    if(nospin == 1)
        index = lattice_index;
    else
        index = 2*lattice_index + spin;
    return(index);
}

/*
    Each data point is an index and a value.
    The index determines which bin the data
    point goes into. The value is added to
    the sum corresponding to that bin.
    If the index value is outside the range
    [lowest, highest] then it is ignored.
    At the end, this will be averaged out. 
*/
int utils_bin_data(DTYPE index, DTYPE value, int bins, int * counts,
                DTYPE lowest, DTYPE bin_width, DTYPE * value_hist,
                DTYPE * errors)
{
    if(isnan(value) || isnan(index))
    {
        printf("NaN value passed\nvalue = %e\nindex = %e\n", value, index);
        return(-1);
    }

    int bin_num = (int) floor((index - lowest) / bin_width);

    (void) bins;
    // if (bin_num >= bins || bin_num < 0)
    // {
    //     printf("Value %e doesn't lie in range: 0 - %e\n", index, bin_width*bins);
    //     return(-1);
    // }

    *(counts + bin_num) += 1;
    *(value_hist + bin_num) += value;
    *(errors + bin_num) += value*value; // We store sum of squares here
    // we will subtract mean^2 later to get stddev

    return(0);
}

/*
    Given the green function matrix for long times, calculates
    the long time imbalance between the initially occupied sites
    and the initially empty sites. The `occupied_set` array must
    be sorted according to numerical value and contains the lattice
    index of the occupied sites.
*/
DTYPE utils_get_charge_imbalance(DTYPE * gfuncsq, int * occupied_set_up,
                        int set_length_up, int * occupied_set_dn,
                        int set_length_dn, int num_states)
{
    int i, j, ctr_up, ctr_dn;
    int site_final, site_initial;
    int spin_final, spin_initial;
    DTYPE imbalance;
    int sign, index1, index2;
    DTYPE elem;

    imbalance = 0;
    ctr_up = 0;
    ctr_dn = 0;

    for(i = 0; i < num_states; i++)
    {
        site_final = i / 2;
        spin_final = i % 2;

        // Check if site in A (occupied set)
        int next_occ_up = *(occupied_set_up + ctr_up);
        int next_occ_dn = *(occupied_set_dn + ctr_dn);

        if(site_final == next_occ_up)
        {
            sign = 1;
            ctr_up += 1;
        }
        else if (site_final == next_occ_dn)
        {
            sign = 1;
            ctr_dn += 1;
        }
        else
            sign = -1;
        
        for(j = 0; j < set_length_up; j++)
        {
            site_initial = *(occupied_set_up + j);
            spin_initial = 0;
            index1 = 2*site_final + spin_final;
            index2 = 2*site_initial + spin_initial;
            elem = *(gfuncsq + RTC(index1, index2, num_states));
            printf("(%d,%d):%d:%e\n", index1, index2, sign, elem);
            imbalance += sign * elem; 
        }

        for(j = 0; j < set_length_dn; j++)
        {
            site_initial = *(occupied_set_dn + j);
            spin_initial = 1;
            index1 = 2*site_final + spin_final;
            index2 = 2*site_initial + spin_initial;
            elem = *(gfuncsq + RTC(index1, index2, num_states));
            printf("(%d,%d):%d:%e\n", index1, index2, sign, elem);
            imbalance += sign * elem; 
        }
    }
    return imbalance / (DTYPE) (set_length_up + set_length_dn);
}




/*
    Creates symmetric matrix from the upper half.
    A lot of symmetric blas routines may only
    write to one half of the matrix. We need to
    reflect the values to the lower half also
    for other parts of the program to work. 
*/
int utils_reflect_upper_to_lower(DTYPE * matrix, int n)
{
    int i, j;
    for(i = 0; i < n; i++)
    {
        for(j = 0; j < i; j++)
        {
            *(matrix + RTC(i,j,n)) = *(matrix + RTC(j,i,n));
        }
    }
    return(0);
}



/*
    Adds matrix2 to matrix1 and stores the result back
    in matrix1. Matrices should be real (DTYPE).
    Assumes column major format.
*/
int utils_add_to_matrix_real(DTYPE * matrix1, DTYPE * matrix2, int m, int n)
{
    int i, j, index;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            *(matrix1 + index) = *(matrix1 + index) + *(matrix2 + index);
        }
    }
    return(0);
}

/*
    Adds matrix2 to matrix1 and stores the result back
    in matrix1. Matrices should be real (DTYPE).
    Assumes column major format.
*/
int utils_add_to_matrix_real_error(DTYPE * matrix1, DTYPE * matrix2, int m, int n, DTYPE * sqsum)
{
    int i, j, index;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            *(matrix1 + index) = *(matrix1 + index) + *(matrix2 + index);
            *(sqsum + index) = *(sqsum + index) + *(matrix2 + index) * (*(matrix2 + index));
        }
    }
    return(0);
}

/*
    Adds matrix2*w2 to matrix1*w1 and stores the result back
    in matrix1. Matrices should be real (DTYPE).
    Assumes column major format.
*/
int utils_add_to_matrix_real_weighted(DTYPE * matrix1, DTYPE * matrix2,
                                    int m, int n, int w1, int w2)
{
    int i, j, index;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            *(matrix1 + index) = *(matrix1 + index) * w1 + *(matrix2 + index) * w2;
        }
    }
    return(0);
}


/*
    Adds matrix2 to matrix1 and stores the result back
    in matrix1. Matrices should be complex (CDTYPE)
    Assumes column major format.
*/
int utils_add_to_matrix_complex(CDTYPE * matrix1, CDTYPE * matrix2, int m, int n)
{
    int i, j, index;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            *(matrix1 + index) = *(matrix1 + index) + *(matrix2 + index);
        }
    }
    return(0);
}


/*
    Adds matrix2 to matrix1 and stores the result back
    in matrix1. Matrices should be complex (CDTYPE)
    Assumes column major format.
*/
int utils_add_to_matrix_complex_error(CDTYPE * matrix1, CDTYPE * matrix2, int m, int n,
                                    CDTYPE * sqsum)
{
    int i, j, index;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            CDTYPE elem = *(matrix2 + index);
            *(matrix1 + index) = *(matrix1 + index) + elem;
            *(sqsum + index) = *(sqsum + index) + (creal(elem)*creal(elem) + I*cimag(elem)*cimag(elem));
        }
    }
    return(0);
}



/*
    Change v to 1-2*v[k,γ] for the initial vector v.
*/
int utils_create_keldysh_vector(void * initial_cond, char type, int num_states)
{
    int i;
    CDTYPE * ptr_c;
    DTYPE *ptr_r;
    for(i = 0; i < num_states; i++)
    {
        if(type == 'C')
        {
            ptr_c = (CDTYPE *) initial_cond + i;
            *ptr_c = 1 - 2* (*ptr_c);
        }
        else if(type == 'R')
        {
            ptr_r = (DTYPE *) initial_cond + i;
            *ptr_r = 1 - 2* (*ptr_r);
        }
        else
        {
            fprintf(stderr, "type should be 'R' or 'C'!\n");
            exit(EXIT_FAILURE);
        }
    }
    return(0);
}

int utils_diag_energysubspace_to_updown(CDTYPE * eigvecs, int Lsq)
{
    // We have v_1 and v_2 as the eigenvectors for an energy
    // subspace
    // Construct v_down = (v_1 + i λ v_2) / (1 + λ^2)^0.5
    // and       v_up   = (v_1 - i λ v_2) / (1 + λ^2)^0.5
    
    int twoLsq = 2 * Lsq;
    CDTYPE lambda = 0;
    DTYPE norm; 
    int n, i, p;
    CDTYPE elem1, elem2;
    int samples = 0;
    
    
    for(n = 0; n < Lsq; n++)
    {
        for(p = 0; p < Lsq; p++)
        {
            elem1 = *(eigvecs + RTC(2*p+DOWN,2*n,twoLsq));
            elem2 = *(eigvecs + RTC(2*p+DOWN,2*n+1,twoLsq));
            if (cabs(elem2) > TOL)
            {
                lambda += I * elem1 / elem2;
                samples++;
            }
        }
        lambda /= samples;
        norm = sqrt(1 + lambda*lambda);

        for(i = 0; i < twoLsq; i++)
        {
            elem1 = *(eigvecs + RTC(i,2*n,twoLsq));
            elem2 = *(eigvecs + RTC(i,2*n+1,twoLsq));
            *(eigvecs + RTC(i,2*n,twoLsq)) = (elem1 + I*lambda*elem2) / norm;
        }


    }
    return(0);
}