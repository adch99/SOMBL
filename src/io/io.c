#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../constants.h"
#include "../params/params.h"
#include "../utils/utils.h"
#include "io.h"

#define MAXCHAR 1024
enum SpinPair {UPUP=0, UPDN, DNUP, DNDN};
enum Spin {UP=0, DOWN};

int check_terminated(const char * str, int max_len)
{
    int i;
    for(i = 0; i < max_len; i++)
    {
        if(*(str + i) == '\0')
            return(1);
    }
    return(0);
}

/*
    Opens the file for the given purpose (must be either 'r'
    or 'w') and checks for errors while opening.
*/
FILE * io_safely_open(char purpose, char * filename)
{
    if(purpose != 'w' && purpose != 'r' && purpose != 'a')
    {
        printf("\nArgument 'purpose' must be either 'r' or 'w' or 'a'\n");
        printf("purpose = %c\n", purpose);
        exit(EXIT_FAILURE);
    }

    FILE * openfile;
    char mode[3];
    sprintf(mode, "%c", purpose);
    // printf("Opening %s in mode %s\n", filename, mode);
    openfile = fopen(filename, mode);
    if (openfile == NULL)
    {
        printf("Error in opening file '%s': %s\n",
                filename, strerror(errno));
        fclose(openfile);
        exit(EXIT_FAILURE);
    }
    return(openfile);
}

/*
    Opens the binary file for the given purpose (must be either 'r'
    or 'w') and checks for errors while opening.
*/
FILE * io_safely_open_binary(char purpose, char * filename)
{
    if(purpose != 'w' && purpose != 'r' && purpose != 'a')
    {
        printf("\nArgument 'purpose' must be either 'r' or 'w' or 'a'\n");
        printf("purpose = %c\n", purpose);
        exit(EXIT_FAILURE);
    }

    FILE * openfile;
    char mode[3];
    sprintf(mode, "%cb", purpose);
    // printf("Opening %s in mode %s\n", filename, mode);
    openfile = fopen(filename, mode);
    if (openfile == NULL)
    {
        printf("Error in opening file '%s': %s\n",
                filename, strerror(errno));
        fclose(openfile);
        exit(EXIT_FAILURE);
    }
    return(openfile);
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
        // printf("io_read_array: m = %d n = %d\n", m, n);
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
    if(info != 0)
    {
        printf("Erring Filename: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    return(info);
}


/*
    Reads the array of the given type from the file and
    stores it in the array in the given ordering. type must
    be 'R' for real, 'C' for complex and 'I' for int.
    ordering must be 'R' for row major and 'C' for column
    major.
*/
int io_read_array_bin(char type, void * array,
                int m, int n, char * filename)
{
    FILE * ifile = io_safely_open_binary('r', filename);
    size_t count;

    if(type == 'C')
    {
        // printf("io_read_array: m = %d n = %d\n", m, n);
        count = fread((CDTYPE *) array, sizeof(CDTYPE), m*n, ifile);
    }
    else if(type == 'R')
    {    
        count = fread((DTYPE *) array, sizeof(DTYPE), m*n, ifile);
    }
    else if(type == 'I')
    {
        count = fread((int *) array, sizeof(int), m*n, ifile);
    }
    else
    {
        printf("ERROR: type given to io_read_array_bin must be one of "
                "the following: 'C', 'R' or 'I'!\n");
        exit(EXIT_FAILURE);
    }
    
    fclose(ifile);
    if(count != m*n)
    {
        printf("Read Error: Only %d out of %d read\n", count, m*n);
        printf("Erring Filename: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    return(count);
}


/*
    Reads the complex array from the file and stores it in
    the array in the given ordering. ordering must be
    'R' for row major and 'C' for column major.
*/
int io_read_array_complex(char ordering, CDTYPE * array,
                        int m, int n, FILE * ifile)
{
    int i, j, info;
    DTYPE zreal, zimag;

    if (ordering == 'C')
    {
        for(i = 0; i < m; i++)
        {
            for(j = 0; j < n; j++)
            {
                info = fscanf(ifile, " (%e%ej) ", &zreal, &zimag);
                *(array + RTC(i, j, m)) = zreal + I*zimag;
            }
        }

    }
    else if (ordering == 'R')
    {
        for(i = 0; i < m; i++)
        {
            for(j = 0; j < n; j++)
            {
                info = fscanf(ifile, " (%e%ej) ", &zreal, &zimag);
                *(array + i*n + j) = zreal + I*zimag;
            }
        }
    }
    else
        return(-1);

    if(info != 2)
    {
        fprintf(stderr, "Input file in wrong format at (%d,%d) info = %d!\nstrerror: %s\n",
                i, j, info, strerror(errno));
        char line[32];
        fgets(line, 32, ifile);
        fprintf(stderr, "Error line: '%s'\n", line);
        fprintf(stderr, "feof: %d ferror: %d\n", feof(ifile), ferror(ifile));
        fprintf(stderr, "m: %d n: %d\n", m, n);
        return(-1);
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
    if(ordering == 'C')
    {
        for(i = 0; i < m; i++)
        {
            for(j = 0; j < n; j++)
            {
                match = fscanf(ifile, "%e", &elem);
                *(array + RTC(i, j, m)) = elem;
            }
        }

    }
    else if(ordering == 'R')
    {
        for(i = 0; i < m; i++)
        {
            for(j = 0; j < n; j++)
            {
                match = fscanf(ifile, "%e", &elem);
                *(array + i*n + j) = elem;
            }
        }

    }
    else
    {
        return(-1);
    }

    if(match == 0)
    {
        printf("An error occured while parsing the file!\n");
        return(-1);
    }

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
    int i, j, match;
    int elem;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            match = fscanf(ifile, "%d", &elem);
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
    Writes the array given to the given file. type must be
    'R' for real, 'C' for complex. ordering must be 'R' for
    row major and 'C' for column major.
*/
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


    if(ordering != 'C' && ordering != 'R')
    {
        printf("Invalid ordering passed to io_write_array.\n");
        printf("ordering must be either 'C' or 'R'.\n");
    }

    // printf("Writing %dx%d\n", m, n);
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = (ordering=='R') ? (i*n + j) : RTC(i, j, m);
            if(type == 'C')
            {
                elemc = *((CDTYPE*)array + index);
                fprintf(ofile, " (%le%+lej) ", crealf(elemc), cimagf(elemc));
            }
            else
            {
                elemr = *((DTYPE*)array + index);
                fprintf(ofile, "%le ", elemr);
            }
        }
        fprintf(ofile, "\n");
    }
    fclose(ofile);
    return 0;
}


/*
    Writes the array given to the given binary file. type must be
    'R' for real, 'C' for complex. ordering must be 'R' for
    row major and 'C' for column major.
*/
int io_write_array_bin(char type, void * array,
                    int m, int n, char * filename)
{
    FILE * ofile = io_safely_open_binary('w', filename);
    size_t count;

    if(type != 'C' && type != 'R')
    {
        printf("Invalid type passed to io_write_array.\n");
        printf("type must be either 'C' or 'R'.\n");
    }

    if(type == 'R')
    {
        count = fwrite((DTYPE *) array, sizeof(DTYPE), m*n, ofile);
    }
    else if(type == 'C')
    {
        count = fwrite((CDTYPE *) array, sizeof(CDTYPE), m*n, ofile);
    }

    if(count != m*n)
    {
        printf("Write Error: Only %d out of %d written\n", count, m*n);
        printf("Erring Filename: %s\n", filename);
        exit(EXIT_FAILURE);        
    }

    fclose(ofile);
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
            fscanf(datafile, "%f", &elem);
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
    fscanf(ifile, "%d", array);
    int ctr = 1;
    for(ctr = 1; ctr < length; ctr++)
    {
        fscanf(ifile, ",%d", array + ctr);
    }



    // char row[MAXCHAR];
    // char * token;
    // int counter = 0;
    // while (!feof(ifile) && counter < length)
    // {
    //     fgets(row, MAXCHAR, ifile);
    //     // printf("Row: %s", row);

    //     if(row[0] == '#')
    //     {
    //         fseek(ifile, -strlen(row), SEEK_CUR);
    //         break;
    //     }

    //     token = strtok(row, ",");

    //     while(token != NULL)
    //     {
    //         *(array + counter) = atoi(token);
    //         counter += 1; 
    //         // printf("Token: %s\n", token);
    //         token = strtok(NULL, ",");
    //     }
    // }

    return(0);
}

int io_read_until_char(char stopchar, FILE * ifile)
{
    char value;
    int found = 0;
    do {
        value = fgetc(ifile);
        if(feof(ifile))
        {
            printf("Reached EOF!\n");
            return(-1);
        }
        else if(ferror(ifile))
        {
            fprintf(stderr, strerror(errno));
            printf("File Reading Error!\n");
            return(-1);
        }
        else if(value == stopchar)
        {
            found = 1;
            break;
        }
    } while(found == 0);

    
    fseek(ifile, -1, SEEK_CUR);
    // clearerr(ifile);
    
    if(found != 1)
    {
        printf("Not found! Last value = %c\n", value);
        return(-1);
    }
    
    return(0);
}


int io_old_read_until_char(char stopchar, FILE * ifile)
{
    // char row[MAXCHAR];
    // int offset, read_len;
    // char * loc = NULL;
    // int passes = 0;

    // while (fgets(row, MAXCHAR, ifile) != NULL)
    // {
    //     passes += 1;
    //     if(check_terminated(row, MAXCHAR) == 0)
    //     {
    //         printf("Something is wrong!\nstring not terminated!\n");
    //         return(-1);
    //     }
    //     read_len = strlen(row);
    //     loc = strchr(row, stopchar);
    //     if(loc == NULL)
    //         continue;
    //     else
    //     {
    //         offset = -read_len + (int)(loc - row);
    //         fseek(ifile, offset, SEEK_CUR);
    //         break;
    //     }
    // }
    // if(ferror(ifile))
    // {
    //     perror("Ferror in reading file!\n");
    //     return(-1);
    // }
    // else if(feof(ifile) && loc == NULL)
    // {
    //     if(stopchar != '\n')
    //         printf("Finding %c in Row: %s FAILED len = %d passes = %d\n", stopchar, row, read_len, passes);
    //     fflush(stdout);
    //     return(-1);
    // }
    // else
    //     return(0);
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
            match = fscanf(ifile, "%e", &elem);
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
            fscanf(ifile, "%e", &zreal);
            // printf("Pointer at %d\n", ftell(ifile));
            fscanf(ifile, "%e", &zimag);
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
        // printf("%e %e %e\n",
        //         *(dists + i), *(gfuncsq + i), *(errors + i));

    }
    fclose(ofile);
    return 0;
}

/*
    Appends the given array to the end of the given file
    in a single line, ending with a newline.
    `type` must be either 'R' for real or 'C' for complex.
    `array` must be a 1D array.
*/
int io_append_array(char type, void * array, int n, char * filename)
{
    if(type != 'C' && type != 'R')
    {
        printf("Invalid type passed to io_write_array.\n");
        printf("type must be either 'C' or 'R'.\n");
    }

    FILE * ofile = io_safely_open('a', filename);
    int i;
    CDTYPE elemc;
    DTYPE elemr;
    for(i = 0; i < n; i++)
    {
        if(type == 'C')
        {
            elemc = *((CDTYPE*)array + i);
            fprintf(ofile, "(%le%+lej) ", crealf(elemc), cimagf(elemc));
        }
        else
        {
            elemr = *((DTYPE*)array + i);
            fprintf(ofile, "%le ", elemr);
        }
    }
    fprintf(ofile, "\n");
    fclose(ofile);
    return(0);
}

/*
    Reads the array of the given type from the file and
    stores it in the array. type must be 'R' for real and
    'C' for complex. Assumes the file to be read is in
    binary double or complex double whereas the array
    is float or complex float.
*/
int io_read_array_bin_d2f(char type, void * array,
                int m, int n, char * filename)
{
    FILE * ifile = io_safely_open_binary('r', filename);
    size_t count;
    double rbuffer;
    complex double cbuffer;
    int i;

    if(type == 'C')
    {
        // printf("io_read_array: m = %d n = %d\n", m, n);
        for(i = 0; i < m*n; i++)
        {
            count = fread(&cbuffer, sizeof(complex double), 1, ifile);
            if(count != 1)
            {
                printf("Read Error: Only %d out of %d read\n", i, m*n);
                printf("Erring Filename: %s\n", filename);
                exit(EXIT_FAILURE);
            }
            *((CDTYPE *) array + i) = (CDTYPE) cbuffer;
        }
    }
    else if(type == 'R')
    {    
        for(i = 0; i < m*n; i++)
        {
            count = fread(&rbuffer, sizeof(double), 1, ifile);
            if(count != 1)
            {
                printf("Read Error: Only %d out of %d read\n", i, m*n);
                printf("Erring Filename: %s\n", filename);
                exit(EXIT_FAILURE);
            }
            *((DTYPE *) array + i) = (DTYPE) rbuffer;
        }
    }
    else
    {
        printf("ERROR: type given to io_read_array_bin_d2f must be one of "
                "the following: 'C' or 'R'!\n");
        exit(EXIT_FAILURE);
    }
    
    fclose(ifile);

    return(count);
}

int io_get_initial_cond_vector(void * initial_cond, char type, char * filename)
{
    int * occupied_set_up;
    int * occupied_set_down;
    int set_length_up, set_length_down;
    int i, index;
    char err_msg[] = "get_initial_cond requires type "
    "as either 'C' for complex or 'R' for real.\n";
    io_get_initial_condition(&occupied_set_up, &set_length_up,
                            &occupied_set_down, &set_length_down,
                            filename);
    printf("set_length_up: %d set_length_down: %d\n", set_length_up, set_length_down);

    if(set_length_up > 0)
    {
        for(i = 0; i < set_length_up; i++)
        {
            index = *(occupied_set_up + i);
            if(type == 'R')
                *((DTYPE *)initial_cond + 2*index + UP) = 1;
            else if(type == 'C')
                *((CDTYPE *)initial_cond + 2*index + UP) = 1;
            else
            {
                fprintf(stderr, err_msg);
                exit(EXIT_FAILURE);
            }
        }
        free(occupied_set_up);
    }
    if(set_length_down > 0)
    {
        for(i = 0; i < set_length_down; i++)
        {
            index = *(occupied_set_down + i);
            if(type == 'R')
                *((DTYPE *)initial_cond + 2*index + DOWN) = 1;
            else if(type == 'C')
                *((CDTYPE *)initial_cond + 2*index + DOWN) = 1;
            else
            {
                fprintf(stderr, err_msg);
                exit(EXIT_FAILURE);
            }
        }
        free(occupied_set_down);
    }
    return(0);
}
