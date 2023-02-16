#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#define L (4)
#define BUFFSIZE (L*L*L*L*2)
int main(int argc, char ** argv)
{
    // char filename[] = "data/mbl_10x10_W15_C0.2_TU1_TD1_N25_BS5_B1_grgrstar_a0_b0_full.dat";
    DTYPE buffer[BUFFSIZE];
    FILE * fp = fopen(*argv, "rb");
    fread(buffer, sizeof(DTYPE), BUFFSIZE, fp);
    fclose(fp);
    int i, j;
    // for(i = 0; i < BUFFSIZE; i++)
    //     printf("%e\n", buffer[i]);
    for(i = 0; i < L*L; i++)
    {
        for(j = 0; j < 2*L*L; j++)
        {
            printf("%le ", buffer[i + j*L*L]);
        }
        printf("\n");
    }

    return(0);
}