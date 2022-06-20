#include <stdio.h>
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
    /* Calculate eigenvalues */
    /* Calculate localization lengths */
    /* Return the localization lengths */
    return(0);
}