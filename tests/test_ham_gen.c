#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/utils/utils.h"
#include "../src/ham_gen/ham_gen.h"
#include "../src/constants.h"
#define TOL 1e-6

// Helper Function
int tester(int (*test_func)(), char * name);

int test_hamiltonian(); 
int test_hamiltonian_nospin();
int test_get_neighbour_list();

int main(int argc, char ** argv)
{
    srandom(17);
    // tester(test_hamiltonian, "test_hamiltonian");
    tester(test_hamiltonian_nospin, "test_hamiltonian_nospin");
    tester(test_get_neighbour_list, "test_get_neighbour_list");
    return 0;
}

int tester(int (*test_func)(), char * name)
{
    int result = test_func();
    if (result)
        printf("Test Passed:\t%s\n", name);
    else
        printf("Test Failed:\t%s\n", name);

    return(0);
}

int test_hermitian(CDTYPE * matrix, int size)
{
    int success = 1;
    int i, j, index, index_flip;

    for(i = 0; i < size; i++)
    {
        for(j = 0; j <= i; j++)
        {
            index = RTC(i, j, size);
            index_flip = RTC(j, i, size);
            if(*(matrix + index) != conj(*  (matrix + index_flip)))
                success = 0;
        }
    }
    return(success);
}

int test_hamiltonian()
{
    int length, width;
    length = width = 3;
    int num_sites = length*width;
    int num_states = 2*num_sites;
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));
    int (*neighbour)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbour, length, width);
    hamiltonian(ham, length, width, 0, 15.0, 1.0, neighbour);

    printf("H (with spin):\n");
    utils_print_matrix(ham, num_states, num_states);

    printf("\nneighbours:\n");
    int i, j;
    for(i = 0; i < num_sites; i++)
    {
        printf("%d: ", i);
        for(j = 0; j < NEIGHS; j++)
        {
            printf("%d ", *(*(neighbour + i) + j));
        }
        printf("\n");
    }

    DTYPE * eigvals = malloc(num_states*sizeof(DTYPE));
    utils_get_eigvalsh(ham, num_states, eigvals);
    printf("\nEigvals:\n");
    for(i = 0; i < num_states; i++)
    {
        printf("%-.2g ", *(eigvals + i));
    }
    return 0;
}

int test_hamiltonian_nospin()
{
    int success = 1;
    int length, width;
    length = width = 10;
    int num_sites = length*width;
    int num_states = num_sites;
    int i, j, index, site1, site2;
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));
    CDTYPE value;
    int (*neighbour)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbour, length, width);
    hamiltonian_nospin(ham, length, width, 15.0, 1.0, neighbour);

    if (test_hermitian(ham, num_states) == 0)
    {
        printf("Matrix not Hermitian!\n");
        success = 0;
    }

    for (i=0;i<num_sites;i++)
    {
        for(j=0;j<NEIGHS;j++)
        {   
            site1 = i;
            site2 = *(*(neighbour + i) + j);
            index = RTC(site1, site2, num_states);
            if (*(ham + index) != -1.0)
            {
                printf("Problem at site (%d,%d) index=%d", site1, site2, index);
                success = 0;
            }
        }
    }

    // printf("H (spinless):\n");
    // utils_print_matrix(ham, num_states, num_states);

    // printf("\nneighbours:\n");
    // for(i = 0; i < num_sites; i++)
    // {
    //     printf("%d: ", i);
    //     for(j = 0; j < NEIGHS; j++)
    //     {
    //         printf("%d ", *(*(neighbour + i) + j));
    //     }
    //     printf("\n");
    // }

    for(i=0;i<num_states;i++)
    {
        for(j=0;j<num_states;j++)
        {
            index = RTC(i, j, num_states);
            value = *(ham + index);

            if(i == j)
            {
                if(value == 0+0*I)
                {
                    printf("Problem at (%d,%d) index=%d!\n", i, j, index);
                    printf("H(i,j) = %e+%ej\n", creal(value), cimag(value));
                    success = 0;
                }
            }
            else if(check_neighbour(j, *(neighbour + i)) >= 0)
            {
                if(value != -1+0*I)
                {
                    printf("Problem at (%d,%d) index=%d!\n", i, j, index);
                    printf("H(i,j) = %e+%ej =/= -1\n", creal(value), cimag(value));
                    success = 0;
                }
            }
            else
            {
                if(value != 0+0*I)
                {
                    printf("Problem at (%d,%d) index=%d\n!", i, j, index);
                    printf("H(i,j) = %e+%ej =/= 0\n", creal(value), cimag(value));
                    success = 0;
                }
            }
        }
    }

    DTYPE * eigvals = malloc(num_states*sizeof(DTYPE));
    utils_get_eigvalsh(ham, num_states, eigvals);
    // printf("\nEigvals:\n");
    // for(i = 0; i < num_states; i++)
    // {
    //     printf("%-.4g ", *(eigvals + i));
    // }
    // printf("\n");
    return success;
}

int test_get_neighbour_list()
{
    int success = 1;
    int length, width;
    length = width = 3;
    int num_sites = length * width;
    int (*neighbours)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbours, length, width);

    int expected[][NEIGHS] = {
        {6, 3, 2, 1},
        {7, 4, 0, 2},
        {8, 5, 1, 0},
        {0, 6, 5, 4},
        {1, 7, 3, 5},
        {2, 8, 4, 3},
        {3, 0, 8, 7},
        {4, 1, 6, 8},
        {5, 2, 7, 6}
    };

    int i, j;
    for(i=0;i<num_sites;i++)
    {
        for(j=0;j<NEIGHS;j++)
        {
            if(expected[i][j] != *(*(neighbours + i) + j))
                success = 0;
        }
    }

    return(success);
}