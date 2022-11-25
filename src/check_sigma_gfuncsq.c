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
    srandom(17);
    int length, width;
    length = width = SIZE;
    int num_sites = length*width;
    int num_states = 2*num_sites;
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));
    int (*neighbour)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbour, length, width);
    hamiltonian(ham, length, width, 2.0, 18.0, 1.0, 1.0, neighbour);

    DTYPE * eigvals = malloc(num_states * sizeof(DTYPE));
    diag_get_eigh(ham, num_states, eigvals);

    // DTYPE * gfuncsq = calloc(num_states*num_states, sizeof(DTYPE));
    CDTYPE * sigma = calloc(2*2, sizeof(CDTYPE));
    *(sigma + RTC(0, 0, 2)) = 1.1;
    *(sigma + RTC(1, 1, 2)) = 0.6;
    *(sigma + RTC(0, 1, 2)) = -0.1;
    *(sigma + RTC(1, 0, 2)) = 0.2;


    int site1 = 15;
    int site2 = 13;
    uint spin1 = UP;
    uint spin2 = DOWN;

    printf("Starting timer...\n");
    time_t start_time = time(NULL);
    DTYPE value = utils_gfuncsq_sigma(site1, spin1, site2, spin2,
                        sigma, ham, num_states);
    time_t duration = time(NULL) - start_time;
    printf("G(%d : %d, %d : %d) = %e\n", site1, spin1, site2, spin2, value);
    printf("Time taken t = %d\n", duration);

    // utils_get_green_func_lim(ham, num_states, gfuncsq, DEGEN_EIGVALS);
    // ofile = fopen("data/matrix_sigma_"SIZESTR"x"SIZESTR"_gfuncsq.dat", "w");
    // utils_save_matrix(gfuncsq, num_states, num_states, 'R', 'F', ofile);
    // fclose(ofile);


    // free(gfuncsq);
    free(eigvals);
    free(ham);
    return(0);
}