#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "constants.h"
#include "utils/utils.h"
#include "ham_gen/ham_gen.h"
#include "diag/diag.h"
// #include "io/io.h"

#define SIZE 40
#define SIZESTR "40"
#define UP 0
#define DOWN 1

int main(int argc, char ** argv)
{
    (void) argc;
    (void) argv;
    // srandom(17);
    srandom(time(NULL));
    int length, width;
    length = width = SIZE;
    int num_sites = length*width;
    int num_states = 2*num_sites;
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));
    int (*neighbour)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbour, length, width);
    hamiltonian(ham, length, width, 2.0, 18.0, 1.0, 1.0, neighbour);
    FILE * ofile = fopen("data/matrix_"SIZESTR"x"SIZESTR"_sample.dat", "w");
    utils_save_matrix(ham, num_states, num_states, 'C', 'F', ofile);
    fclose(ofile);
    DTYPE * eigvals = malloc(num_states * sizeof(DTYPE));
    // diag_get_eigh(ham, num_states, eigvals);
    diag_get_eigvalsh(ham, num_states, eigvals);

    int i;
    printf("\nEigvals: ");
    for(i = 0; i < num_states; i++)
        printf("%e ", *(eigvals + i));
    printf("\n");

    // ofile = fopen("data/matrix_"SIZESTR"x"SIZESTR"_eigvecs.dat", "w");
    // utils_save_matrix(ham, num_states, num_states, 'C', 'F', ofile);
    // fclose(ofile);

    // DTYPE * gfuncsq = calloc(num_states*num_states, sizeof(DTYPE));
    // int site1 = 15;
    // unsigned int spin1 = DOWN;
    // int site2 = 11;
    // unsigned int spin2 = DOWN;
    // int index1 = 2*site1+spin1;
    // int index2 = 2*site2+spin2;
    // DTYPE value = utils_compute_gfsq_elem(index1, index2, ham,
    //                                     num_states, DEGEN_EIGVALS);
    // printf("G(%d : %d, %d : %d) = %e\n", site1, spin1, site2, spin2, value);

    // utils_get_green_func_lim(ham, num_states, gfuncsq, DEGEN_EIGVALS);
    // ofile = fopen("data/matrix_"SIZESTR"x"SIZESTR"_gfuncsq.dat", "w");
    // utils_save_matrix(gfuncsq, num_states, num_states, 'R', 'F', ofile);
    // fclose(ofile);

    // int j;
    // int murs_scheme[4] = {1, 3, 0, 2};
    // printf("Neighbours:\n");
    // printf("     x+1\t  y+1\t  x-1\t  y-1\t\n");
    // for(i = 0; i < num_sites; i++)
    // {
    //     printf("%d: ", i);
    //     for(j = 0; j < 4; j++)
    //     {
    //         printf("%5d\t", *(*(neighbour + i) + murs_scheme[j]));
    //     }
    //     printf("\n");
    // }
    // CDTYPE x = I;
    // printf("%e + %ej\n", crealf(x), cimagf(x));

    // free(gfuncsq);
    free(eigvals);
    free(ham);
    return(0);
}