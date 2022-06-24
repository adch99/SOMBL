#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "ham_gen/ham_gen.h"
#include "utils/utils.h"

double run();

int main(int argc, char ** argv)
{
    int ctr;
    int numRuns = 50;
    char outfilename[] = "data/loclengths.dat";
    FILE * outfile = fopen(outfilename, "w");
    DTYPE avg_loc_len;

    printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = 1; ctr <= numRuns; ctr++)
    {
        printf("Run %d started...", ctr);
        /* Call run */
        avg_loc_len = run();
        fprintf(outfile, "%lf\n", avg_loc_len);
        printf("\tDone\n");
    }
    fclose(outfile);
    return(0);
}

double run()
{
    /* Set parameters */
    /* Create hamiltonian */
    int len, width;
    DTYPE coupling_const, disorder_strength, hop_strength, energy;

    len = width = 100;
    coupling_const = 0;
    disorder_strength = 15;
    hop_strength = 1;
    energy = 1;

    int num_sites = len * width;
    int num_states = 2*num_sites;

    int (*neighbours)[NEIGHS] = malloc(sizeof(DTYPE[num_sites][NEIGHS])); 
    CDTYPE * ham = malloc(sizeof(CDTYPE)*(num_states*num_states));
    
    get_neighbour_lists(neighbours, len, width);
    hamiltonian(ham, len, width, coupling_const, disorder_strength,
                hop_strength, neighbours);

    // utils_print_matrix(ham, num_states, num_states);
    /* Calculate eigenvalues */
    // DTYPE * eigvals;
    DTYPE * eigvals = calloc(num_states, sizeof(DTYPE));
    utils_get_eigvalsh(ham, num_states, eigvals);

    /* Calculate localization lengths */
    int i;
    DTYPE avg_loc_len = 0;
    for(i = 0; i < num_states; i++)
    {
        avg_loc_len += utils_loc_len(-1, eigvals, hop_strength, len, i); 
    }
    avg_loc_len /= num_states;
    
    /* Return the localization lengths */

    free(ham);
    free(neighbours); // Remove this later. We should reuse neighbour list.
    return(avg_loc_len);
}

