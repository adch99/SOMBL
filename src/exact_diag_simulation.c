#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "ham_gen/ham_gen.h"
#include "utils/utils.h"

struct SystemParams {
    int len;
    int width;
    DTYPE coupling_const;
    DTYPE disorder_strength;
    DTYPE hop_strength;
    int (*neighbours)[NEIGHS];
};


double run(struct SystemParams * params, int create_neighbours);

int main(int argc, char ** argv)
{
    struct SystemParams params;
    params.len = params.width = 10;
    params.coupling_const = 0;
    params.disorder_strength = 15;
    params.hop_strength = 1;

    int ctr;
    int numRuns = 50;
    char outfilename[32];
    sprintf(outfilename, "data/loclengths%dx%d.dat", params.len, params.width);
    FILE * outfile = fopen(outfilename, "w");
    DTYPE avg_loc_len;
    int create_neighbours = 1;

    printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = 1; ctr <= numRuns; ctr++)
    {
        printf("Run %d started...", ctr);
        fflush(stdout);
        /* Call run */
        avg_loc_len = run(&params, create_neighbours);
        create_neighbours = 0;
        fprintf(outfile, "%lf\n", avg_loc_len);
        printf("\tDone\n");
    }
    fclose(outfile);
    free(params.neighbours);
    return(0);
}

double run(struct SystemParams * params, int create_neighbours)
{
    /* Set parameters */
    /* Create hamiltonian */

    DTYPE energy = 1;

    int num_sites = params->len * params->width;
    int num_states = 2*num_sites;

    if (create_neighbours)
    {
        params->neighbours = malloc((num_sites*NEIGHS)*sizeof(int)); 
        get_neighbour_lists(params->neighbours, params->len, params->width);
    }

    CDTYPE * ham = malloc(sizeof(CDTYPE)*(num_states*num_states));
    
    hamiltonian(ham, params->len, params->width, params->coupling_const,
                params->disorder_strength, params->hop_strength,
                params->neighbours);

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
        avg_loc_len += utils_loc_len(-1, eigvals, params->hop_strength,
                                    params->len, i); 
    }
    avg_loc_len /= num_states;
    
    /* Return the localization lengths */

    free(eigvals);
    free(ham);
    // free(neighbours); // Remove this later-> We should reuse neighbour list->
    return(avg_loc_len);
}

