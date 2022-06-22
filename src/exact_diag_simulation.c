#include <stdio.h>
#include <stdlib.h>
#include "ham_gen/ham_gen.h"

double run();

int main(int argc, char ** argv)
{
    int ctr;
    int numRuns = 10;

    printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = 1; ctr <= numRuns; ctr++)
    {
        printf("Run %d started...", ctr);
        /* Call run */
        run();
        printf("\tDone\n");
    }
    return(0);
}

double run()
{
    /* Set parameters */
    /* Create hamiltonian */
    int len, width;
    DTYPE coupling_const, disorder_strength, hop_strength;

    len = width = 3;
    coupling_const = 0;
    disorder_strength = 15;
    hop_strength = 1;

    int num_sites = len * width;

    int ** neighbours = malloc(sizeof(DTYPE)*(num_sites*NEIGHS)); 
    CDTYPE * ham = malloc(sizeof(DTYPE)*(num_sites*num_sites));
    
    hamiltonian(ham, len, width, coupling_const, disorder_strength,
                hop_strength, neighbours);
    /* Calculate eigenvalues */
    /* Calculate localization lengths */
    /* Return the localization lengths */

    free(ham);
    free(neighbours); // Remove this later. We should reuse neighbour list.
    return(0);
}

