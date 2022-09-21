#include <stdlib.h>
#include <string.h>
#include "../constants.h"
#include "../params/params.h"
#include "../utils/utils.h"

#define MAXCHAR 1000

int io_read_until(char * stopstr, FILE * ifile);
int io_get_array_from_file(int * array, int length, FILE * ifile);


int io_get_gfuncsq_from_file(DTYPE * matrix, struct OutStream outfiles,
                    struct SystemParams params)
{
    FILE * datafile = fopen(outfiles.gfuncsq, "r");

    if(datafile == NULL)
    {
        printf("Cannot open file %s!\n", outfiles.gfuncsq);
        exit(-1);
    }

    int i, j, index;
    DTYPE elem;
    for(i = 0; i < params.num_states; i++)
    {
        for(j = 0; j < params.num_states; j++)
        {
            index = RTC(i, j, params.num_states);
            fscanf(datafile, "%lf", &elem);
            *(matrix + index) = elem;
        }
    }
    fclose(datafile);
    return 0;
}

int io_get_initial_condition(int ** occupied_set_up, int * set_length_up,
                            int ** occupied_set_dn, int * set_length_dn,
                            char * filename)
{
    FILE * ifile = fopen(filename, "r");
    if(ifile == NULL)
    {
        printf("Cannot open file %s!\n", filename);
        return(-1);
    }

    char row[MAXCHAR];

    if(io_read_until("#NUMUP\n", ifile) == -1)
    {
        printf("Format not valid! Cannot detect #NUMUP\n");
        return(-1);
    }

    fgets(row, MAXCHAR, ifile);
    *set_length_up = atoi(row);
    *occupied_set_up = malloc(*set_length_up * sizeof(int));

    if(io_read_until("#NUMDN\n", ifile) == -1)
    {
        printf("Format not valid! Cannot detect #NUMDN\n");
        return(-1);
    }

    fgets(row, MAXCHAR, ifile);
    *set_length_dn = atoi(row);
    *occupied_set_dn = malloc(*set_length_dn * sizeof(int));

    if(io_read_until("#SPINSUP\n", ifile) == -1)
    {
        printf("Format not valid! Cannot detect #SPINSUP\n");
        return(-1);
    }

    io_get_array_from_file(*occupied_set_up, *set_length_up, ifile);

    if(io_read_until("#SPINSDN\n", ifile) == -1)
    {
        printf("Format not valid! Cannot detect #SPINSDN\n");
        return(-1);
    }

    io_get_array_from_file(*occupied_set_dn, *set_length_dn, ifile);

    fclose(ifile);
    return(0);
}


int io_read_until(char * stopstr, FILE * ifile)
{
    char row[MAXCHAR];
    while (!feof(ifile))
    {
        fgets(row, MAXCHAR, ifile);
        // printf("Row: %s", row);
        if(strcmp(row, stopstr) == 0)
            break;
    }
    if(feof(ifile))
        return(-1);
    else
        return(0);
}

int io_get_array_from_file(int * array, int length, FILE * ifile)
{
    char row[MAXCHAR];
    char * token;
    int counter = 0;
    while (!feof(ifile) && counter < length)
    {
        fgets(row, MAXCHAR, ifile);
        // printf("Row: %s", row);

        if(row[0] == '#')
        {
            fseek(ifile, -strlen(row), SEEK_CUR);
            break;
        }

        token = strtok(row, ",");

        while(token != NULL)
        {
            *(array + counter) = atoi(token);
            counter += 1; 
            // printf("Token: %s\n", token);
            token = strtok(NULL, ",");
        }
    }

    return(0);
}

int io_output_function_data(DTYPE * dists, DTYPE * gfuncsq,
                        char * filename, int data_len)
{
    // Write the values to a file
    FILE * ofile = fopen(filename, "w");
    if (ofile == NULL)
    {
        printf("Cannot open file %s!", filename);
        return(-1);
    }
    int i;
    for(i = 0; i < data_len; i++)
    {
        // if(isnan(*(dists + i)) || isnan(*(gfuncsq + i)))
        //     printf("NaN detected in postprocess.");

        fprintf(ofile, "%e %e\n", *(dists + i), *(gfuncsq + i));
        printf("%e %e\n", *(dists + i), *(gfuncsq + i));

    }
    fclose(ofile);
    return 0;
}
