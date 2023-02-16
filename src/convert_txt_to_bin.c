#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "io/io.h"
#include "constants.h"

int main(int argc, char ** argv)
{
    char filename[64] = "data/sample_matrix_10x10.dat";
    char newfilename[64] = "data/sample_matrix_bin_10x10.dat";
    CDTYPE * bufferc;
    int m, n;
    m = 10;
    n = 10;
    bufferc = calloc(m*n, sizeof(CDTYPE));
    io_read_array('C', 'C', bufferc, m, n, filename);
    io_write_array_bin('C', bufferc, m, n, newfilename);
    free(bufferc);

    // Verification
    // CDTYPE * arrayc = calloc(m*n, sizeof(CDTYPE));
    // CDTYPE elem;
    // io_read_array_bin('C', arrayc, m, n, newfilename);
    // int i, j;
    // for(i = 0; i < m; i++)
    // {
    //     for(j = 0; j < n; j++)
    //     {
    //         elem = *(arrayc + i + m*j);
    //         printf("(%e+%ej) ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // free(arrayc);

    return 0;
}