#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../constants.h"
#include "../params/params.h"
#include "../utils/utils.h"
#include "io.h"

#define MAXCHAR 1000


/*
    Opens the file for the given purpose (must be either 'r'
    or 'w') and checks for errors while opening.
*/
FILE * io_safely_open(char purpose, char * filename)
{
    FILE * openfile;
    if(purpose == 'w')
    {
        openfile = fopen(filename, "w");
        if (openfile == NULL)
        {
            printf("Error in opening file %s: %s\n",
                    filename, strerror(errno));
            fclose(openfile);
            exit(EXIT_FAILURE);
        }
        return(openfile);
    }

    if(purpose == 'r')
    {
        openfile = fopen(filename, "r");
        if (openfile == NULL)
        {
            printf("Error in opening file %s: %s\n",
                    filename, strerror(errno));
            fclose(openfile);
            exit(EXIT_FAILURE);
        }
        return(openfile);
    }
    else
    {
        printf("Argument 'purpose' must be either 'r' or 'w'\n");
        exit(EXIT_FAILURE);
    }
}

/*
    Reads the array of the given type from the file and
    stores it in the array in the given ordering. type must
    be 'R' for real, 'C' for complex and 'I' for int.
    ordering must be 'R' for row major and 'C' for column
    major.
*/
int io_read_array(char type, char ordering, void * array,
                int m, int n, char * filename)
{
    FILE * ifile = io_safely_open('r', filename);
    int info;

    if(ordering != 'C' && ordering != 'R')
    {
        printf("Error: ordering must be eiher 'C' for"
                "column major or 'R' for row major.");
        exit(EXIT_FAILURE);
    }

    if(type == 'C')
    {
        info = io_read_array_complex(ordering, (CDTYPE *) array,
                                    m, n, ifile);
    }
    else if(type == 'R')
    {        
        info = io_read_array_real(ordering, (DTYPE *) array,
                                    m, n, ifile);
    }
    else if(type == 'I')
    {
        info = io_read_array_int(ordering, (int *) array,
                                    m, n, ifile);
    }
    else
    {
        printf("ERROR: type given to io_read_array must be one of "
                "the following: 'C', 'R' or 'I'!\n");
        exit(EXIT_FAILURE);
    }
    
    fclose(ifile);
    return(info);
}

/*
    Reads the complex array from the file and stores it in
    the array in the given ordering. ordering must be
    'R' for row major and 'C' for column major.
*/
int io_read_array_complex(char ordering, CDTYPE * array,
                        int m, int n, FILE * ifile)
{
    int i, j, found;
    DTYPE zreal, zimag;
    // printf("Beginning\n");
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            // printf("Pointer at %d\n", ftell(ifile));
            found = io_read_until_char('(', ifile);
            if(found != 0)
            {
                printf("Something is wrong! Cannot parse file!\n");
                return(-1);
            }
            // printf("Pointer at %d\n", ftell(ifile));
            fseek(ifile, 1, SEEK_CUR);
            // printf("Pointer at %d\n", ftell(ifile));
            fscanf(ifile, "%le", &zreal);
            // printf("Pointer at %d\n", ftell(ifile));
            fscanf(ifile, "%le", &zimag);
            // printf("Pointer at %d\n\n", ftell(ifile));
            // printf("z = %le+%lej\n", zreal, zimag);

            if (ordering == 'C')
                *(array + RTC(i, j, m)) = zreal + I*zimag;
            else if (ordering == 'R')
                *(array + i*n + j) = zreal + I*zimag;
            else
                return(-1);
        }
        io_read_until_char('\n', ifile);
        // printf("Pointer at %d\n", ftell(ifile));
        fseek(ifile, 1, SEEK_CUR);
        // printf("Pointer at %d\n-------\n", ftell(ifile));
    }
    return(0);
}

/*
    Reads the real array from the file and stores it in
    the array in the given ordering. ordering must be
    'R' for row major and 'C' for column major.
*/
int io_read_array_real(char ordering, DTYPE * array,
                        int m, int n, FILE * ifile)
{
    int i, j, match;
    DTYPE elem;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            match = fscanf(ifile, "%le", &elem);
            if(match == 0)
            {
                printf("An error occured while parsing the file!\n");
                return(-1);
            }
            if (ordering == 'C')
                *(array + RTC(i, j, m)) = elem;
            else if (ordering == 'R')
                *(array + i*n + j) = elem;
            else
                return(-1);
        }
        io_read_until_char('\n', ifile);
        fseek(ifile, 1, SEEK_CUR);
        // printf("\n");
    }
    // printf("\n");
    return(0);
}


/*
    Reads the int array from the file and stores it in
    the array in the given ordering. ordering must be
    'R' for row major and 'C' for column major.
*/
int io_read_array_int(char ordering, int * array,
                        int m, int n, FILE * ifile)
{
    return(0);
}

int io_write_array(char type, char ordering, void * array,
                    int m, int n, char * filename)
{
    FILE * ofile = io_safely_open('w', filename);
    int i, j, index;
    DTYPE elemr;
    CDTYPE elemc;


    if(type != 'C' && type != 'R')
    {
        printf("Invalid type passed to io_write_array.\n");
        printf("type must be either 'C' or 'R'.\n");
    }


    if(ordering != 'C' && ordering != 'F')
    {
        printf("Invalid ordering passed to io_write_array.\n");
        printf("ordering must be either 'C' or 'F'.\n");
    }


    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = (ordering=='C')?(i*n + j):RTC(i, j, m);
            if(type == 'C')
            {
                elemc = *((CDTYPE*)array + index);
                fprintf(ofile, "(%le%+lej) ", creal(elemc), cimag(elemc));
            }
            else
            {
                elemr = *((DTYPE*)array + index);
                fprintf(ofile, "%le ", elemr);
            }
        }
        fprintf(ofile, "\n");
    }
    return 0;
}


int io_get_gfuncsq_from_file(DTYPE * matrix, struct OutStream outfiles,
                    struct SystemParams params, int sigma)
{
    FILE * datafile = fopen(outfiles.gfuncsq, "r");

    if(datafile == NULL)
    {
        printf("Cannot open file %s!\n", outfiles.gfuncsq);
        exit(-1);
    }

    int i, j, index;
    int size;
    if (sigma == 0)
        size = params.num_states;
    else
        size = params.num_sites;
    DTYPE elem;
    
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            index = RTC(i, j, size);
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

int io_read_until_char(char stopchar, FILE * ifile)
{
    char row[MAXCHAR];
    int offset, read_len;
    while (!feof(ifile))
    {
        if(fgets(row, MAXCHAR, ifile) == NULL)
            break;
        read_len = strlen(row);
        // printf("Row: %s", row);
        char * loc = strchr(row, stopchar);
        if(loc == NULL)
            continue;
        else
        {
            offset = -read_len + (int)(loc - row);
            fseek(ifile, offset, SEEK_CUR);
            break;
        }
    }
    if(feof(ifile))
        return(-1);
    else
        return(0);
}

/*
    Reads 2d real array of size m x n from given file
    into the array.
*/
int io_dread_2d(DTYPE * array, int m, int n, FILE * ifile)
{
    int i, j, match;
    DTYPE elem;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            match = fscanf(ifile, "%le", &elem);
            if(match == 0)
            {
                printf("An error occured while parsing the file!\n");
                return(-1);
            }
            *(array + RTC(i, j, m)) = elem;
        }
        io_read_until_char('\n', ifile);
        fseek(ifile, 1, SEEK_CUR);
        // printf("\n");
    }
    // printf("\n");
    return(0);
}

/*
    Reads 2d complex array of size m x n from given file
    into the array.
*/
int io_zread_2d(CDTYPE * array, int m, int n, FILE * ifile)
{
    int i, j, found;
    DTYPE zreal, zimag;
    // printf("Beginning\n");
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            // printf("Pointer at %d\n", ftell(ifile));
            found = io_read_until_char('(', ifile);
            if(found != 0)
            {
                printf("Something is wrong! Cannot parse file!");
                return(-1);
            }
            // printf("Pointer at %d\n", ftell(ifile));
            fseek(ifile, 1, SEEK_CUR);
            // printf("Pointer at %d\n", ftell(ifile));
            fscanf(ifile, "%le", &zreal);
            // printf("Pointer at %d\n", ftell(ifile));
            fscanf(ifile, "%le", &zimag);
            // printf("Pointer at %d\n\n", ftell(ifile));
            // printf("z = %le+%lej\n", zreal, zimag);

            *(array + RTC(i, j, m)) = zreal + I*zimag;
        }
        io_read_until_char('\n', ifile);
        // printf("Pointer at %d\n", ftell(ifile));
        fseek(ifile, 1, SEEK_CUR);
        // printf("Pointer at %d\n-------\n", ftell(ifile));
    }
    return(0);
}


int io_output_function_data(DTYPE * dists, DTYPE * gfuncsq,
                        DTYPE * errors, char * filename,
                        int data_len)
{
    // Write the values to a file
    FILE * ofile = io_safely_open('w', filename);
    int i;
    for(i = 0; i < data_len; i++)
    {
        // if(isnan(*(dists + i)) || isnan(*(gfuncsq + i)))
        //     printf("NaN detected in postprocess.");

        fprintf(ofile, "%e %e %e\n",
                *(dists + i), *(gfuncsq + i), *(errors + i));
        printf("%e %e %e\n",
                *(dists + i), *(gfuncsq + i), *(errors + i));

    }
    fclose(ofile);
    return 0;
}
