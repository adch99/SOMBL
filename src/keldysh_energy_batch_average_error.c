#include <stdio.h>
#include <stdlib.h>
#include "utils/utils.h"
#include "params/params.h"
#include "io/io.h"

/* Constant Declarations */
#define MAXLEN 64
const char *argp_program_version =
    "keldysh_energy_batch_average 1.0";
const char *argp_program_bug_address =
    "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
    "keldysh_energy_batch_average -- a simulation of spin-orbit coupled" 
    "2d many-body localized systems.";
// A description of the arguments we accept.
static char args_doc[] = "-s <size> -c <coupling_const>"
                        "-w <disorder_strength>"
                        " -t <hop_strength> -n <num_runs>"
                        " [-p]";
// The options we understand.
static struct argp_option options[] = {
    {"size",      's', "SIZE",       0, "Length and width of the lattice",        0},
    {"coupling",  'c', "COUPLING",   0, "Spin-orbit coupling constant",           0},
    {"disorder",  'w', "DISORDER",   0, "Strength of the disorder",               0},
    {"hopping",   't', "HOPPING",    0, "Strength of the hopping",                0},
    {"hopup",     'u', "HOPUP",      0, "Strength of the hopping for up spins",   0},
    {"hopdn",     'd', "HOPDN",      0, "Strength of the hopping for down spins", 0},
    {"runs",      'n', "NUMRUNS",    0, "Number of runs in the disorder average", 0},
    {"nospin",    'p', 0,            0, "Use a spinless model hamiltonian.",      0},
    {"batch",     'b', "BATCHNUM",   0, "Batch number of this set of runs.",      0},
    {"batchsize", 'j', "BATCHSIZE",  0, "Number of runs in each batch.",          0},
    {"seed",      'x', "SEED",       0, "Seed for random number generation.",     0},
    {"bins",      'e', "ENERGYBINS", 0, "Number of bins for energy binning.",     0},
    {"dtypes",    'a', "DTYPELIST",  0, "Word having dtypes for each batch",      0},
    { 0 }
};
// Our argp parser.
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};

//------------------------------------------------------------------------
int average_real_matrix(char * ofilename, char ** ifilenames,
                        int num_batches, int m, int n,
                        char * dtypelist);
int average_complex_matrix(char * ofilename, char ** ifilenames,
                        int num_batches, int m, int n,
                        char * dtypelist);
int get_real_matrix(char * filename, char dtype,
                    DTYPE * array, int m, int n);
int get_complex_matrix(char * filename, char dtype,
                    CDTYPE * array, int m, int n);



int main(int argc, char ** argv)
{
    struct SystemParams params;
    struct OutStream outfiles;
    uint alpha, beta;
    int bin, num_batches, L, Lsq;
    int batchnum;
    char ofilename[MAXLEN];
    char ** ifilenames;

    params_setup(argc, argv, &params, &outfiles, &argp, 0);
    params.energybins = 0;
    L = params.len;
    Lsq = L*L;
    num_batches = params.numRuns / params.batch_size;
    // allocate the filenames
    ifilenames = malloc(num_batches*sizeof(char *));
    for(batchnum = 1; batchnum <= num_batches; batchnum++)
        *(ifilenames + batchnum - 1) = malloc(MAXLEN*sizeof(char));

    printf("dtypes: %s\n", params.dtypelist);

    for(alpha = 0; alpha < 2; alpha++)
    {        
        printf("α = %d β = %d\n", alpha, alpha);
        printf("-------------\n");
        // We start from bin = -1 as that accounts for full
        // window case as well.
        printf("energybins: %d\n", params.energybins);
        for(bin = -1; bin < params.energybins; bin++)
        {
            printf("bin: %d\n", bin);
            params.batch_num = -1;
            params_gr_grstar_filename(ofilename, params,
                                    bin, alpha, alpha);
            printf("Output: %s\n", ofilename);

            printf("Input:\n");
            for(batchnum = 1; batchnum <= num_batches; batchnum++)
            {
                params.batch_num = batchnum;
                params_gr_grstar_filename(*(ifilenames + batchnum - 1),
                                        params, bin, alpha, alpha);
                printf("%s\n", *(ifilenames + batchnum - 1));
            }
            printf("\n");
            average_real_matrix(ofilename, ifilenames, num_batches,
                                Lsq, 2*Lsq, params.dtypelist);
        }
    }

    // alpha != beta case
    alpha = 0;
    beta = 1;
    printf("α = %d β = %d\n", alpha, beta);
    printf("-------------\n");

    // We start from bin = -1 as that accounts for full
    // window case as well.
    for(bin = -1; bin < params.energybins; bin++)
    {
        params.batch_num = -1;
        params_gr_grstar_filename(ofilename, params,
                                bin, alpha, beta);
        printf("Output: %s\n", ofilename);

        printf("Input:\n");
        for(batchnum = 1; batchnum <= num_batches; batchnum++)
        {
            params.batch_num = batchnum;
            params_gr_grstar_filename(*(ifilenames + batchnum - 1),
                                    params, bin, alpha, beta);
            printf("%s\n", *(ifilenames + batchnum - 1));
        }
        printf("\n");
        average_complex_matrix(ofilename, ifilenames, num_batches,
                            Lsq, 2*Lsq, params.dtypelist);
    }


    for(batchnum = 1; batchnum <= num_batches; batchnum++)
        free(*(ifilenames + batchnum - 1));
    free(ifilenames);

    return(0);
}

int average_real_matrix(char * ofilename, char ** ifilenames,
                        int num_batches, int m, int n,
                        char * dtypelist)
{
    DTYPE * avg_gfuncsq = calloc(m*n, sizeof(DTYPE));
    DTYPE * avg_sq_gfuncsq = calloc(m*n, sizeof(DTYPE));
    DTYPE * batch_gfuncsq = calloc(m*n, sizeof(DTYPE));
    char * ifilename;
    int batchnum, i, j;

    // For each batch
    for(batchnum = 1; batchnum <= num_batches; batchnum++)
    {
        // Generate the correct filename
        ifilename = *(ifilenames + batchnum - 1);
        // Read the array and add it to the sum
        printf("Reading file %s...\n", ifilename);
        char dtype = *(dtypelist + batchnum - 1);
        get_real_matrix(ifilename, dtype, batch_gfuncsq, m, n);
        // io_read_array('R', 'C', batch_gfuncsq, m, n, ifilename);
        utils_add_to_matrix_real_error(avg_gfuncsq, batch_gfuncsq,
                                    m, n, avg_sq_gfuncsq);
    }

    // Divide to get the average
    int index;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            *(avg_gfuncsq + index) /= num_batches;
            *(avg_sq_gfuncsq + index) = *(avg_sq_gfuncsq + index)/num_batches
                                    - (*(avg_gfuncsq + index))*(*(avg_gfuncsq + index));
        }
    }


    // Save this array in a file
    printf("Writing to file %s...\n", ofilename);
    io_write_array_bin('R', avg_gfuncsq, m, n, ofilename);

    char varfilename[128];
    sprintf(varfilename, "%s.variance", ofilename);
    io_write_array_bin('R', avg_sq_gfuncsq, m, n, varfilename);

    free(avg_gfuncsq);
    free(avg_sq_gfuncsq);
    free(batch_gfuncsq);
    return(0);
}

int average_complex_matrix(char * ofilename, char ** ifilenames,
                            int num_batches, int m, int n,
                            char * dtypelist)
{
    CDTYPE * avg_gfuncsq = calloc(m*n, sizeof(CDTYPE));
    CDTYPE * avg_sq_gfuncsq = calloc(m*n, sizeof(CDTYPE));
    CDTYPE * batch_gfuncsq = calloc(m*n, sizeof(CDTYPE));
    char * ifilename;
    int batchnum, i, j;

    // For each batch
    for(batchnum = 1; batchnum <= num_batches; batchnum++)
    {
        // Generate the correct filename
        ifilename = *(ifilenames + batchnum - 1);
        // Read the array and add it to the sum
        printf("Reading file %s...\n", ifilename);
        // io_read_array('C', 'C', batch_gfuncsq, m, n, ifilename);
        char dtype = *(dtypelist + batchnum - 1);
        get_complex_matrix(ifilename, dtype, batch_gfuncsq, m, n);

        utils_add_to_matrix_complex_error(avg_gfuncsq, batch_gfuncsq, m, n, avg_sq_gfuncsq);
    }

    // Divide to get the average
    int index;
    CDTYPE elem;
    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = RTC(i, j, m);
            *(avg_gfuncsq + index) /= num_batches;
            elem = *(avg_gfuncsq + index);
            *(avg_sq_gfuncsq + index) = *(avg_sq_gfuncsq + index)/num_batches - (creal(elem)*creal(elem) + I*cimag(elem)*cimag(elem));

        }
    }

    // Save this array in a file
    printf("Writing to file %s...\n", ofilename);
    io_write_array_bin('C', avg_gfuncsq, m, n, ofilename);

    char varfilename[128];
    sprintf(varfilename, "%s.variance", ofilename);
    io_write_array_bin('C', avg_sq_gfuncsq, m, n, varfilename);


    free(avg_gfuncsq);
    free(avg_sq_gfuncsq);
    free(batch_gfuncsq);
    return(0);
}

int get_real_matrix(char * filename, char dtype, DTYPE * array,
                    int m, int n)
{
    if (dtype == 't')
        io_read_array('R', 'C', array, m, n, filename);
    else if (dtype == 'd')
        io_read_array_bin_d2f('R', array, m, n, filename);
    else if(dtype == 'f')
        io_read_array_bin('R', array, m, n, filename);
    else
    {
        fprintf(stderr, "dtype passed to get_real_matrix should be either 't', 'd' or 'f'!\n");
        fprintf(stderr, "passed value is %c\n", dtype);
        return(-1);
    }
    return(0);
}

int get_complex_matrix(char * filename, char dtype, CDTYPE * array,
                    int m, int n)
{
    if (dtype == 't')
        io_read_array('C', 'C', array, m, n, filename);
    else if (dtype == 'd')
        io_read_array_bin_d2f('C', array, m, n, filename);
    else if(dtype == 'f')
        io_read_array_bin('C', array, m, n, filename);
    else
    {
        fprintf(stderr, "dtype passed to get_real_matrix should be either 't', 'd' or 'f'!\n");
        fprintf(stderr, "passed value is %c\n", dtype);
        return(-1);
    }
    return(0);
}