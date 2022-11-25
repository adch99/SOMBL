#ifndef MBL_IO_H
#define MBL_IO_H

#include "../constants.h"

int io_get_gfuncsq_from_file(DTYPE * matrix, struct OutStream outfiles,
                    struct SystemParams params);
int io_get_initial_condition(int ** occupied_set_up, int * set_length_up,
                            int ** occupied_set_dn, int * set_length_dn,
                            char * filename);
int io_output_function_data(DTYPE * dists, DTYPE * gfuncsq,
                        DTYPE * errors, char * filename,
                        int data_len);
int io_zread_2d(CDTYPE * array, int m, int n, FILE * ifile);
int io_dread_2d(DTYPE * array, int m, int n, FILE * ifile);

#endif