#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "io/io.h"
#include "utils/utils.h"

int main(int argc, char ** argv)
{
    if(argc != 6)
    {
        fprintf(stderr, "ERROR: 5 cmd line args are needed: fname1 n1 fname2 n2 ofname\n");
        fprintf(stderr, "Only %d args recieved.\n", argc);        
        exit(1);
    }
    int L = 100;
    int Lsq = L*L;
    int twoLsq = 2*Lsq;
    twoLsq = 10;
    printf("Running with wrong L value: twoLsq = %d\n", twoLsq);
    
    char fname1[128];
    char fname2[128];
    char ofname[128];
    strcpy(fname1, *(argv + 1));
    int n1 = atoi(*(argv + 2));
    strcpy(fname2, *(argv + 3));
    int n2 = atoi(*(argv + 4));
    strcpy(ofname, *(argv + 5));


    // char fname1[] = "data/mbl_100x100_W12_C1.3_TU1_TD1_N90_BS15_B3_greenfuncsq.dat";
    // int n1 = 5;
    // char fname2[] = "data/mbl_100x100_W12_C1.3_TU1_TD1_N90_BS15_B3_greenfuncsq.dat";
    // int n2 = 5;
    // char ofname[] = "data/mbl_100x100_W12_C1.3_TU1_TD1_N90_BS15_B3_greenfuncsq.dat.new";

    DTYPE * array = calloc(twoLsq*twoLsq, sizeof(DTYPE));
    DTYPE * buffer = calloc(twoLsq*twoLsq, sizeof(DTYPE));

    // Read the array and add it to the sum
    printf("Reading file %s...\n", fname1);
    io_read_array('R', 'C', buffer, twoLsq, twoLsq, fname1);
    utils_add_to_matrix_real_weighted(array, buffer, twoLsq, twoLsq, 0, n1);

    // Read the array and add it to the sum
    printf("Reading file %s...\n", fname2);
    io_read_array('R', 'C', buffer, twoLsq, twoLsq, fname2);
    utils_add_to_matrix_real_weighted(array, buffer, twoLsq, twoLsq, 1, n2);

    int i, j;
    for(i = 0; i < twoLsq; i++)
    {
        for(j = 0; j < twoLsq; j++)
        {
            *(array + RTC(i,j,twoLsq)) /= (DTYPE) (n1 + n2);
            printf("%8.2f ", *(array + RTC(i,j,twoLsq)));
        }
        printf("\n");
    }



    printf("Writing file %s...\n", ofname);
    io_write_array('R', 'C', array, twoLsq, twoLsq, ofname);

    free(array);
    free(buffer);
}