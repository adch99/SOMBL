#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cblas.h"
#include "constants.h"
#include "utils/utils.h"

#define SIZE 40

int get_green_func_lim(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy);
DTYPE compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy);

int main(int argc, char ** argv)
{
    (void) argc; (void) argv;

    int num_states = SIZE*SIZE*2;
    CDTYPE * eigvecs = calloc(num_states*num_states, sizeof(CDTYPE));
    DTYPE * gfuncsq = calloc(num_states*num_states, sizeof(DTYPE));

    int i, j;
    printf("Initialization started\n");
    for(i = 0; i < num_states; i++)
    {
        for(j = 0; j < num_states; j++)
        {
            *(eigvecs + RTC(i, j, num_states)) = 0.01*(i + I*j);
        }
    }
    printf("Done.\n");

    printf("Gfunc construction started.\n");
    time_t start_time = time(NULL);   
    get_green_func_lim(eigvecs, num_states, gfuncsq, DEGEN_EIGVALS);
    time_t end_time = time(NULL);

    printf("Done. Time taken is %lds.\n", (end_time - start_time));

    return(0);
}

int get_green_func_lim(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy)
{

    int i, j;

    // First let's square the eigenvectors matrix
    // element by element
    DTYPE * eigvec_sq = calloc(size*size, sizeof(DTYPE));
    CDTYPE elem;

    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            elem = *(eigenvectors + RTC(i, j, size));
            *(eigvec_sq + RTC(i, j, size)) = creal(elem)*creal(elem) + cimag(elem)*cimag(elem); 
        }
    }

    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, size, size,
                1.0, eigvec_sq, size, 1.0, green_func, size);

    
    // #pragma omp parallel for collapse(2)
    int index_ij = 0;
    int index_ji;
    for(i = 0; i < size; i++)
    {
        index_ji = i;
        for(j = 0; j <= i; j++)
        {
            DTYPE sum = compute_gfsq_elem(i, j, eigenvectors, size,
                                                degeneracy);
            // #pragma omp critical
            *(green_func + index_ij) += sum;
            if(i != j)
                *(green_func + index_ji) += sum;
            index_ji += size;
            index_ij++;
        }
    }    
    // free(eigvec_sq);
    return 0;
}

/*
    Computes the (i,j)th element of the Green's function
    squared in the long time limit.
*/
DTYPE compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy)
{
    DTYPE sum = 0.0;
    int k, l;

    // int index1, index2;
    // index1 = RTC(i, 0, size);
    // index2 = RTC(j, 0, size);
    // CDTYPE psi1, psi2;
    // DTYPE mod1, mod2;

    // // #pragma omp parallel for reduction (+:sum) schedule(auto)
    // for(k = 0; k < size; k++)
    // {
    //     psi1 = *(eigenvectors + index1);
    //     psi2 = *(eigenvectors + index2);
    //     mod1 = creal(psi1)*creal(psi1) + cimag(psi1)*cimag(psi1);
    //     mod2 = creal(psi2)*creal(psi2) + cimag(psi2)*cimag(psi2);
    //     sum += mod1 * mod2;
    //     index1 += size;
    //     index2 += size;
    // }

    if(degeneracy == DEGEN_EIGVALS)
    {
        CDTYPE * ptr_ik = eigenvectors + RTC(i, 0, size);
        CDTYPE * ptr_il = eigenvectors + RTC(i, 1, size);
        CDTYPE * ptr_jk = eigenvectors + RTC(j, 0, size);
        CDTYPE * ptr_jl = eigenvectors + RTC(j, 1, size);
        CDTYPE term;

        // #pragma omp parallel for reduction (+:sum) schedule(auto)
        for(k = 0; k < size; k += 2)
        {
            // l = k + 1;

            term = (*ptr_ik) * conj(*ptr_il) * conj(*ptr_jk) * (*ptr_jl);
            sum += term + conj(term);

            ptr_ik += 2*size;            
            ptr_il += 2*size;            
            ptr_jk += 2*size;            
            ptr_jl += 2*size;            
        }
    }
    return(sum);
}
