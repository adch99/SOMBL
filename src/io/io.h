#ifndef MBL_IO_H
#define MBL_IO_H

#include "../constants.h"
#include "../params/params.h"


FILE * io_safely_open(char purpose, char * filename);
FILE * io_safely_open_binary(char purpose, char * filename);
int io_read_array(char type, char ordering, void * array,
                int m, int n, char * filename);
int io_read_array_bin(char type, void * array,
                int m, int n, char * filename);
int io_read_array_bin_d2f(char type, void * array,
                int m, int n, char * filename);
int io_read_array_real(char ordering, DTYPE * array,
                        int m, int n, FILE * ifile);
int io_read_array_complex(char ordering, CDTYPE * array,
                        int m, int n, FILE * ifile);
int io_read_array_int(char ordering, int * array,
                        int m, int n, FILE * ifile);
int io_write_array(char type, char ordering, void * array,
                    int m, int n, char * filename);
int io_write_array_bin(char type, void * array,
                    int m, int n, char * filename);

int io_read_until(char * stopstr, FILE * ifile);
int io_read_until_char(char stopchar, FILE * ifile);
int io_get_array_from_file(int * array, int length, FILE * ifile);


int io_get_gfuncsq_from_file(DTYPE * matrix, struct OutStream outfiles,
                    struct SystemParams params, int sigma);
int io_get_initial_condition(int ** occupied_set_up, int * set_length_up,
                            int ** occupied_set_dn, int * set_length_dn,
                            char * filename);
int io_output_function_data(DTYPE * dists, DTYPE * gfuncsq,
                        DTYPE * errors, char * filename,
                        int data_len);
int io_zread_2d(CDTYPE * array, int m, int n, FILE * ifile);
int io_dread_2d(DTYPE * array, int m, int n, FILE * ifile);
int io_append_array(char type, void * array, int n, char * filename);
#endif