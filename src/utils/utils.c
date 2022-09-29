#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <lapacke.h>
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
    given is not an eigenvalue. This works only for 1D systems.
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
    Returns the eigenvalues of the hermitian matrix given in
    the array of eigvals given. We assume the array has
    already been preprocessed for the used library
    (currently LAPACK). NOTE: The matrix given to
    diagonalize will be destroyed by this function.
*/
int utils_get_eigvalsh(CDTYPE * matrix, int size, DTYPE * eigvals)
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
int utils_get_eigh(CDTYPE * matrix, int size, DTYPE * eigvals)
{
    // printf("Running zheev...");
    // fflush(stdout);
    lapack_int info, lda, n, lwork;
    double * rwork;
    lapack_complex_double * work, query;

    lda = size;
    n = size;
    rwork = calloc(3*n-2, sizeof(double));

    lwork = -1;
    info = LAPACKE_zheev_work(LAPACK_COL_MAJOR, 'V', 'U', n,
                        matrix, lda, eigvals, &query, lwork, rwork);
    
    lwork = query;
    work = calloc(query, sizeof(lapack_complex_double));
    info = LAPACKE_zheev_work(LAPACK_COL_MAJOR, 'V', 'U', n,
                        matrix, lda, eigvals, work, lwork, rwork);
    
    if (info != 0)
    {
        printf("LAPACKE_zheev error! Code: %d", info);
        return(info); // Some error has occured.
    }

    free(work);
    free(rwork);

    return(0);
}


/*
    Generates 'num_samples' numbers in the range [low, high]
    using a uniform distribution. Please ensure that you
    give a reasonable number for low and high i.e.
    low < high.
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

int utils_save_matrix(void * matrix, int m, int n,
                    char type, char ordering, FILE * ofile)
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
                fprintf(ofile, "%e+%ej ", creal(elemc), cimag(elemc));
            }
            else
            {
                elemr = *((DTYPE*)matrix + index);
                fprintf(ofile, "%e ", elemr);
            }
        }
        fprintf(ofile, "\n");
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
int utils_get_green_func_lim(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy)
{

    // Construct green function limit squared
    int i, j;
    // int num_threads;
    // #pragma omp critical
    // {
    //     num_threads = omp_get_num_threads();
    // }
    
    // #pragma omp parallel for collapse(2)
    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            DTYPE sum = utils_compute_gfsq_elem(i, j, eigenvectors, size,
                                                degeneracy);
            // #pragma omp critical
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
DTYPE utils_compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy)
{
    DTYPE sum = 0.0;
    int k, l;
    // #pragma omp parallel for reduction (+:sum) schedule(auto)
    for(k = 0; k < size; k++)
    {
        int index1, index2;
        CDTYPE psi1, psi2;
        DTYPE mod1, mod2;
        index1 = RTC(i, k, size);
        index2 = RTC(j, k, size);
        psi1 = *(eigenvectors + index1);
        psi2 = *(eigenvectors + index2);
        mod1 = creal(psi1)*creal(psi1) + cimag(psi1)*cimag(psi1);
        mod2 = creal(psi2)*creal(psi2) + cimag(psi2)*cimag(psi2);
        sum += mod1 * mod2;
    }

    if(degeneracy == DEGEN_EIGVALS)
    {
        // We assume that the degeneracy occurs in
        // pairs and that we have all the eigenvalues
        // being degenerate (Kramer degeneracy).
        for(k = 0; k < size; k += 2)
        {
            l = k + 1;
            CDTYPE psi_k_i, psi_l_i;
            CDTYPE psi_k_j, psi_l_j;
            CDTYPE term;
            psi_k_i = *(eigenvectors + RTC(i, k, size));
            psi_l_i = *(eigenvectors + RTC(i, l, size));
            psi_k_j = *(eigenvectors + RTC(j, k, size));
            psi_l_j = *(eigenvectors + RTC(j, l, size));

            term = psi_k_i * conj(psi_l_i) * conj(psi_k_j) * psi_l_j;
            sum += term + conj(term);
        }
    }
    return(sum);
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

    if (bin_num >= bins || bin_num < 0)
    {
        printf("Value %e doesn't lie in range: 0 - %e\n", index, bin_width*bins);
        return(-1);
    }

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


DTYPE utils_pbc_chord_length_sq(int index1, int length1, int index2, int length2)
{
    DTYPE radius1 = (DTYPE) length1 / (2.0 * M_PI);
    DTYPE radius2 = (DTYPE) length2 / (2.0 * M_PI);

    DTYPE angle1 = 2.0 * M_PI * ((DTYPE) index1 / (DTYPE) length1);
    DTYPE angle2 = 2.0 * M_PI * ((DTYPE) index2 / (DTYPE) length2);

    DTYPE chord1 = 2 * radius1 * sin(angle1 / 2.0);
    DTYPE chord2 = 2 * radius2 * sin(angle2 / 2.0);

    DTYPE chordsq = chord1*chord1 + chord2*chord2;
    return chordsq;   
}

