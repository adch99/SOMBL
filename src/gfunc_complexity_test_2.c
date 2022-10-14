#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <cblas.h>
#include "utils/utils.h"


// #define RTC(i,j,size) ((i) + (j)*(size))

int init_matrix(double complex * A, int size);
int check_matrix_colmajor(double * A, double * B, int m,
                        int n, double tol);
int the_right_way(double * G, double complex * A, int size);
int the_easier_way(double * G, double complex * A, int size);

int the_reference_way(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy);
DTYPE reference_way_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy);


int main(int argc, char ** argv)
{
    (void) argc; (void) argv;
    int size = 500;
    double complex * A = calloc(size*size, sizeof(double complex));
    // double A[4] = {1, 3, 2, 4};
    init_matrix(A, size);
    time_t start_time;


    double * G1 = calloc(size*size, sizeof(double));
    start_time = time(NULL);
    the_easier_way(G1, A, size);
    printf("G1 took %lds.\n", (time(NULL) - start_time));

    double * G2 = calloc(size*size, sizeof(double));
    start_time = time(NULL);
    the_right_way(G2, A, size);
    printf("G2 took %lds.\n", (time(NULL) - start_time));

    double * G3 = calloc(size*size, sizeof(double));
    start_time = time(NULL);
    the_reference_way(A, size, G3, DEGEN_EIGVALS);
    printf("G3 took %lds.\n", (time(NULL) - start_time));

    double * G4 = calloc(size*size, sizeof(double));
    // the_current_way
    start_time = time(NULL);
    utils_get_green_func_lim(A, size, G4, DEGEN_EIGVALS);
    printf("G4 took %lds.\n", (time(NULL) - start_time));

    double diff = 0;
    int i;
    for(i = 0; i < size*size; i++)
    {
        diff += fabs(*(G1 + i) - *(G2 + i));
    }
    printf("Diff: %e\n", diff);

    int success = check_matrix_colmajor(G1, G2, size, size, 1e-3);
    success = success && check_matrix_colmajor(G1, G3, size, size, 1e-3);
    success = success && check_matrix_colmajor(G1, G4, size, size, 1e-3);
    if(success == 1)
        printf("All matrices are equal.\n");
    else
    {
        printf("G1:\n");
        utils_print_matrix(G1, size, size, 'R', 'F');

        printf("\nG2:\n");
        utils_print_matrix(G2, size, size, 'R', 'F');

        printf("\nG3:\n");
        utils_print_matrix(G3, size, size, 'R', 'F');

        printf("\nG4:\n");
        utils_print_matrix(G4, size, size, 'R', 'F');
    }

    free(G1);
    free(G2);
    // free(A);
    return(0);
}

int init_matrix(double complex * A, int size)
{
    int i, j;
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            *(A + RTC(i,j,size)) = 0.2*(i + 0.34*j) + 0.3*I*(-0.254*i + 0.52*j);
        }
    }
    return(0);
}

int check_matrix_colmajor(double * A, double * B, int m,
                        int n, double tol)
{
    int i, j;
    int success = 1;
    double diff; 
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            diff = *(A + i + j*n) - *(B + i + j*n);
            if(fabs(diff) > tol)
            {
                printf("Issue at (%d,%d)! ", i, j);
                printf("Diff = %12.5e\n", diff);
                success = 0;
            }
        }
    }
    return(success);
}

int the_right_way(double * G, double complex * A, int size)
{
    int i, j, k;
    double mod1, mod2, sum = 0;
    double complex term;
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            sum = 0;
            for(k = 0; k < size; k++)
            {
                mod1 = cabs(*(A + RTC(i,k,size)));
                mod2 = cabs(*(A + RTC(j,k,size)));
                sum += mod1 * mod1 * mod2 * mod2;   
            }

            // Degenerate subspace part
            for(k = 0; k < size; k+=2)
            {
                term = conj(*(A + RTC(i,k,size))) * (*(A + RTC(j,k,size))) * conj(*(A + RTC(j,k+1,size))) * (*(A + RTC(i,k+1,size)));
                sum += 2 * creal(term);
            }

            *(G + RTC(i,j,size)) += sum;
        }
    }
    return(0);
}

int the_easier_way(double * G, double complex * A, int size)
{
    double * Asq = calloc(size*size, sizeof(double));
    int i, j, n;
    double mod;
    for(i = 0; i < size*size; i++)
    {
        mod = cabs(*(A + i));
        *(Asq + i) = mod*mod;
    }

    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, size, size,
            1.0, Asq, size, 1.0, G, size);
    // Only the upper half matrix is updated.
    // Copy that to lower half
    // for(i = 0; i < size; i++)
    // {
    //     for(j = 0; j < i; j++)
    //     {
    //         *(G + RTC(i, j, size)) = *(G + RTC(j, i, size));
    //     }
    // }

    free(Asq);

    printf("Hello\n");

    double complex * B = calloc(size*size/2, sizeof(double complex));
    double complex * BBh = calloc(size*size, sizeof(double complex));

    for(i = 0; i < size; i++)
    {
        for(n = 0; n < size/2; n++)
        {
            *(B + RTC(i,n,size)) = conj(*(A + RTC(i,2*n,size))) * (*(A + RTC(i,2*n+1,size)));
        }
    }
    cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,
                size, size/2, 1.0, B, size, 0.0, BBh, size);
    printf("Bello\n");
    for(i = 0; i < size; i++)
    {
        for(j = i; j < size; j++)
        {
            *(G + RTC(i,j,size)) += 2*creal(*(BBh + RTC(i, j, size)));
        }
    }

    // Only the upper half matrix is updated.
    // Copy that to lower half
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < i; j++)
        {
            *(G + RTC(i, j, size)) = *(G + RTC(j, i, size));
        }
    }

    free(B);
    free(BBh);
    return(0);
}

int the_reference_way(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy)
{

    // Construct green function limit squared
    int i, j;
    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            DTYPE sum = reference_way_elem(i, j, eigenvectors, size,
                                                degeneracy);
            *(green_func + i*size + j) += sum;
            if(i != j)
                *(green_func + j*size + i) += sum;
        }
    }
    return 0;
}

/*
    Computes the (i,j)th element of the Green's function
    squared in the long time limit.
*/
DTYPE reference_way_elem(int i, int j, CDTYPE * eigenvectors,
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
