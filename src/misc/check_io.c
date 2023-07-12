#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "io/io.h"
#include "utils/utils.h"

#define L 4

int main(void)
{
    // char fname[] = "data/mbl_10x10_W15.1_C0.04_TU1_TD1_N100_BS10_B3_grgrstar_a0_b1_bin0.dat";
    // char fname[] = "data/mbl_10x10_W15_C0.02_TU1_TD1_N100_BS10_B3_grgrstar_a0_b1_bin0.dat";
    char fname[] = "data/mbl_4x4_W15_C0.02_TU1_TD1_N100_BS10_B3_grgrstar_a0_b1_bin0.dat";
    // char fname[] = "data/test_matrix.dat";
    // int m = 3;
    // int n = 4;
    int m = L*L;
    int n = 2*L*L;
    printf("m = %d n = %d\n", m, n);
    CDTYPE * matrix = calloc(m * n, sizeof(CDTYPE));
    CDTYPE * matrix2 = calloc(m * n, sizeof(CDTYPE));
    io_read_array('C', 'C', matrix2, m, n, fname);
    // io_read_array('C', 'C', matrix, m, n, fname);
    FILE * ifile = io_safely_open('r', fname);
    int i, j, info;
    DTYPE zreal, zimag;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            info = fscanf(ifile, " (%le%lej) ", &zreal, &zimag);
            if(info != 2)
            {
                fclose(ifile);
                fprintf(stderr, "Wrong input format (%d,%d) info = %d\n", i, j, info);
                exit(EXIT_FAILURE);
            }
            *(matrix + i + j*m) = zreal + I*zimag;
            // printf("%d ", info);
            printf("%6.1e%+7.1ej ", zreal, zimag);
            if(cabs(*(matrix + i + j*m) - *(matrix2 + i + j*m)) > 1e-6)
                printf("Difference at %d,%d\n", i, j);
        }
        // io_read_until_char('\n', ifile);
        // fscanf(ifile, "\n");
        printf("\n");
    }
    fclose(ifile);

    // FILE * ofile = io_safely_open('w', "data/test_matrix_copy.dat");
    io_write_array('C', 'C', matrix, m, n, "data/test_matrix_copy.dat");
    // utils_print_matrix(matrix, m, n, 'C', 'C');
    // fclose(ofile);
    free(matrix);
    return(0);
}