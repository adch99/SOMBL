#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "../src/utils/utils.h"
#include "../src/ham_gen/ham_gen.h"
#include "../src/constants.h"
#define TOL 1e-6

// Helper Function
int tester(int (*test_func)(), char * name);

// Smaller functions
int test_ham_lattice_site(int index, CDTYPE * ham, int length, DTYPE hopup,
                        DTYPE hopdn, DTYPE disorder, DTYPE coupling);
int test_ceq(CDTYPE a, CDTYPE b, DTYPE tol);

// Functions to test
int test_hamiltonian(); 
int test_hamiltonian_nospin();
int test_get_neighbour_list_pbc();
int test_get_neighbour_list_obc();

int main(int argc, char ** argv)
{
    (void) argc;
    (void) argv;
    // srandom(17);
    // tester(test_hamiltonian, "test_hamiltonian");
    // tester(test_hamiltonian_nospin, "test_hamiltonian_nospin");
    tester(test_get_neighbour_list_obc, "test_get_neighbour_list_obc");
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

int test_ceq(CDTYPE a, CDTYPE b, DTYPE tol)
{
    DTYPE diff_real = fabs(creal(a) - creal(b));
    DTYPE diff_imag = fabs(cimag(a) - cimag(b));
    return((diff_real < tol) && (diff_imag < tol));
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
            if(*(matrix + index) != conj(*(matrix + index_flip)))
                success = 0;
        }
    }
    return(success);
}

int test_ham_lattice_site(int index, CDTYPE * ham, int length, DTYPE hopup,
                        DTYPE hopdn, DTYPE disorder, DTYPE coupling)
{
    int success = 1;
    int num_states = 2 * length * length;
    int x, y, cond, nb_index;
    unsigned int spin;
    CDTYPE value, expected, hopping;
    utils_get_lattice_index(index, length, 0, &x, &y, &spin);

    CDTYPE spin_orbit_exp[4][2] = {{-1, 1}, {1, -1}, {I, I}, {-I, -I}};

    if(spin == 0)
        hopping = hopup;
    else if(spin == 1)
        hopping = hopdn;
    else
    {
        printf("Error: spin is not 0 or 1! spin = %d\n", spin);
        success = 0;
    }

    // Check if disorder strength is in range
    value = *(ham + RTC(index, index, num_states));
    cond = ((cabs(value) <= disorder/2) && (cabs(value) >= -disorder/2));
    if(cond == 0)
        printf("Error at (%d, %d) diagonal\n", index, index);
    success = success && cond; 

    // (x-1, y) same spin
    expected = -hopping;
    nb_index = utils_get_matrix_index(x-1, y, spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x-1, y, spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x-1, y) opp spin
    expected = coupling * spin_orbit_exp[0][1-spin];
    nb_index = utils_get_matrix_index(x-1, y, 1-spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x-1, y, 1-spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x+1, y) same spin
    expected = -hopping;
    nb_index = utils_get_matrix_index(x+1, y, spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x+1, y, spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x+1, y) opp spin
    expected = coupling * spin_orbit_exp[1][1-spin];
    nb_index = utils_get_matrix_index(x+1, y, 1-spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x+1, y, 1-spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x, y-1) same spin
    expected = -hopping;
    nb_index = utils_get_matrix_index(x, y-1, spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x, y-1, spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x, y-1) opp spin
    expected = coupling * spin_orbit_exp[2][1-spin];
    nb_index = utils_get_matrix_index(x, y-1, 1-spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x, y-1, 1-spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x, y+1) same spin
    expected = -hopping;
    nb_index = utils_get_matrix_index(x, y+1, spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x, y+1, spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    //(x, y+1) opp spin
    expected = coupling * spin_orbit_exp[3][1-spin];
    nb_index = utils_get_matrix_index(x, y+1, 1-spin, length, 0);
    value = *(ham + RTC(nb_index, index, num_states));
    cond = (cabs(value - expected) < TOL);
    if(cond == 0)
    {
        printf("Error in ham:\n");
        printf("i = %d = (%d, %d, %d)\n", nb_index, x, y+1, 1-spin);
        printf("j = %d = (%d, %d, %d)\n", index, x, y, spin);
        printf("value = %.3e + i%.3e\nexpected = %.3e + i%.3e\n\n", creal(value), cimag(value),
                creal(expected), cimag(expected));
    }
    success = success && cond;

    return(success);
}

int test_hamiltonian()
{
    int success = 1;
    int length, width;
    length = width = 3;
    int num_sites = length*width;
    int num_states = 2*num_sites;
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));
    int (*neighbour)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbour, length, width);
    hamiltonian(ham, length, width, 0.1, 15.0, 1.0, 1.5, neighbour);

    // printf("H (with spin):\n");
    // utils_print_matrix(ham, num_states, num_states, 'C', 'F');

    // int i, result;
    // for(i = 0; i < num_states; i++)
    // {
    //     result = test_ham_lattice_site(i, ham, length, 1.0, 1.5, 15.0, 0.1);
    //     success = (success && result);
    // }

    // printf("\nneighbours:\n");
    // int i, j;
    // for(i = 0; i < num_sites; i++)
    // {
    //     printf("%d: ", i);
    //     for(j = 0; j < NEIGHS; j++)
    //     {
    //         printf("%d ", *(*(neighbour + i) + j));
    //     }
    //     printf("\n");
    // }

    // DTYPE * eigvals = malloc(num_states*sizeof(DTYPE));
    // utils_get_eigvalsh(ham, num_states, eigvals);
    // printf("\nEigvals:\n");
    // for(i = 0; i < num_states; i++)
    // {
    //     printf("%-.2g ", *(eigvals + i));
    // }

    free(ham);
    free(neighbour);
    return(success);
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
            if (test_ceq(*(ham + index), -1+0*I, TOL) == 0)
            {
                printf("Problem at site (%d,%d) index=%d\n", site1, site2, index);
                success = 0;
            }
        }
    }

    for(i=0;i<num_states;i++)
    {
        for(j=0;j<num_states;j++)
        {
            index = RTC(i, j, num_states);
            value = *(ham + index);

            if(i == j)
            {
                if(test_ceq(value, 0+0*I, TOL) == 1)
                {
                    printf("Problem at (%d,%d) index=%d!\n", i, j, index);
                    printf("Result: %d\n", test_ceq(value, 0+0*I, TOL));
                    printf("H(i,j) = %e+%ej\n", creal(value), cimag(value));
                    success = 0;
                }
            }
            else if(check_neighbour(j, *(neighbour + i)) >= 0)
            {
                if(test_ceq(value, -1+0*I, TOL) == 0)
                {
                    printf("Problem at (%d,%d) index=%d!\n", i, j, index);
                    printf("H(i,j) = %e+%ej =/= -1\n", creal(value), cimag(value));
                    success = 0;
                }
            }
            else
            {
                if(test_ceq(value, 0+0*I, TOL) == 0)
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
    free(ham);
    free(neighbour);
    free(eigvals);
    return success;
}

int test_get_neighbour_list_pbc()
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

    free(neighbours);
    return(success);
}

int test_get_neighbour_list_obc()
{
    int success = 1;
    int length, width;
    length = width = 3;
    int num_sites = length * width;
    int (*neighbours)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbours, length, width);

    int expected[][NEIGHS] = {
        {-1,  1, -1,  3},
        { 0,  2, -1,  4},
        { 1, -1, -1,  5},
        {-1,  4,  0,  6},
        { 3,  5,  1,  7},
        { 4, -1,  2,  8},
        {-1,  7,  3, -1},
        { 6,  8,  4, -1},
        { 7, -1,  5, -1}
    };

    int i, j;
    for(i=0;i<num_sites;i++)
    {
        // printf("%d: ", i);
        for(j=0;j<NEIGHS;j++)
        {
            // printf("%3d ", *(*(neighbours + i) + j));
            if(expected[i][j] != *(*(neighbours + i) + j))
                success = 0;
        }
        // printf("\n");
    }
    // printf("\n");

    free(neighbours);
    return(success);
}


int test_check_neighbour()
{
    int success = 1;
    return(success);
}