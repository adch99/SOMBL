#include <stdlib.h>
#include <stdio.h>
#include <argp.h>
#include <math.h>
#include "constants.h"
#include "params/params.h"
#include "utils/utils.h"

/* Constant Declarations */
const char *argp_program_version =
  "calculate_dist_vs_gfuncsq 1.0";
const char *argp_program_bug_address =
  "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
  "calculate_dist_vs_gfuncsq -- converts the green function " 
  "squared matrix from data file to a function of the "
  "distances between lattice sites";
// A description of the arguments we accept.
static char args_doc[] = "-s <size> -c <coupling_const>"
                        "-w <disorder_strength>"
                        " -t <hop_strength> -n <num_runs>"
                        " [-p]";
// The options we understand.
static struct argp_option options[] = {
  {"size",     's', "SIZE",     0, "Length and width of the lattice",        0},
  {"coupling", 'c', "COUPLING", 0, "Spin-orbit coupling constant",           0},
  {"disorder", 'w', "DISORDER", 0, "Strength of the disorder",               0},
  {"hopping",  't', "HOPPING",  0, "Strength of the hopping",                0},
  {"runs",     'n', "NUMRUNS",  0, "Number of runs in the disorder average", 0},
  {"nospin",   'p', 0,          0, "Use a spinless model hamiltonian.",      0},
  { 0 }
};
// Our argp parser.
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};


int setup(int argc, char ** argv, struct SystemParams * params,
        struct OutStream * outfiles);
int get_gfuncsq_data(DTYPE * matrix, struct OutStream outfiles,
                    struct SystemParams params);
int output_function_data(DTYPE * dists, DTYPE * gfuncsq,
                        struct OutStream outfiles, int data_len);


int main(int argc, char ** argv)
{
    // Get SystemParams from defaults or cmd line args
    struct SystemParams params;
    struct OutStream outfiles;

    setup(argc, argv, &params, &outfiles);

    printf("Calculate G^2(r)\n");
    printf("----------------\n");
    printf("len: %d\t nospin: %d\t coupling: %.2e\n"
        "disorder: %.2e\t hopping: %.2e\n",
        params.len, params.nospin, params.coupling_const,
        params.disorder_strength, params.hop_strength);


    // Get the gfuncsq matrix from the data file
    DTYPE * matrix = malloc(params.num_states*params.num_states*sizeof(DTYPE));
    get_gfuncsq_data(matrix, outfiles, params);
    // utils_print_matrix(matrix, params.num_states, params.num_states, 'R', 'F');
    // Bin the gfunsq matrix and create a function of
    // distance G^2(|r_i - r_j|)
    DTYPE * dists;
    DTYPE * gfuncsq;
    int data_len;
    // Let's bin the data into bins of width 1.
    int bins = ceil((DTYPE) params.len / 2.0);
    data_len = utils_construct_data_vs_dist(matrix, params.num_states, params.len,
                                        bins, &dists, &gfuncsq);
    output_function_data(dists, gfuncsq, outfiles, data_len);

    free(matrix);
    free(dists);
    free(gfuncsq);
    return 0;
}

int setup(int argc, char ** argv, struct SystemParams * params,
        struct OutStream * outfiles)
{
    /* Default Values of Parameters */
    params->len = params->width = 20;
    params->coupling_const = 0;
    params->disorder_strength = 10;
    params->hop_strength = 1;
    params->numRuns = 100;
    params->nospin = 0;

    /* Parse our arguments; every option seen by parse_opt will
        be reflected in params. */
    argp_parse (&argp, argc, argv, 0, 0, params);
    if(params->nospin)
        params->num_states = params->len*params->width;
    else
        params->num_states = 2*params->len*params->width;

    // Set up the data files for the system params.
    *outfiles = params_set_up_datastream(*params);

    return 0;
}

int get_gfuncsq_data(DTYPE * matrix, struct OutStream outfiles,
                    struct SystemParams params)
{
    FILE * datafile = fopen(outfiles.gfuncsq, "r");

    if(datafile == NULL)
    {
        printf("Cannot open file %s!\n", outfiles.gfuncsq);
        return(-1);
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

int output_function_data(DTYPE * dists, DTYPE * gfuncsq,
                        struct OutStream outfiles, int data_len)
{
    // Write the values to a file
    FILE * ofile = fopen(outfiles.dist_vs_gfuncsq, "w");
    if (ofile == NULL)
    {
        printf("Cannot open file %s!", outfiles.dist_vs_gfuncsq);
        return(-1);
    }
    int i;
    for(i = 0; i < data_len; i++)
    {
        if(isnan(*(dists + i)) || isnan(*(gfuncsq + i)))
            printf("NaN detected in postprocess.");

        fprintf(ofile, "%e %e\n", *(dists + i), *(gfuncsq + i));
        printf("%e %e\n", *(dists + i), *(gfuncsq + i));

    }
    fclose(ofile);
    return 0;
}