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
  {"spinless", 'p', 0,          0, "Use a spinless model hamiltonian.",      0},
  { 0 }
};
// Our argp parser.
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};


int setup(int argc, char ** argv, struct SystemParams * params,
        struct OutStream * outfiles);




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
    // Bin the gfunsq matrix and create a function of
    // distance G^2(|r_i - r_j|)

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

int get_gfuncsq_data()
{
    return 0;
}

DTYPE post_process(struct SystemParams params, struct OutStream outfiles,
                DTYPE * gfunc)
{
    int i;
    int num_states = params.num_states;


    // Construct gfunc vs distance datapoints
    DTYPE * dists;
    DTYPE * gfuncsq_vals;
    int data_len;
    int bins = 10;
    DTYPE exponent, mantissa, residuals;
    // exponent = mantissa = residuals = 0;
    data_len = utils_construct_data_vs_dist(gfunc, num_states, params.len,
                                        bins, &dists, &gfuncsq_vals);


    // Write the values to a file
    FILE * ofile = fopen(outfiles.dist_vs_gfuncsq, "w");
    for(i = 0; i < data_len; i++)
    {
        if(isnan(*(dists + i)) || isnan(*(gfuncsq_vals + i)))
            printf("NaN detected in postprocess.");

        fprintf(ofile, "%e %e\n", *(dists + i), *(gfuncsq_vals + i));
        printf("%e %e\n", *(dists + i), *(gfuncsq_vals + i));

    }
    fclose(ofile);

    // Fit exponential to the data
    utils_fit_exponential(dists, gfuncsq_vals, data_len,
                        &exponent, &mantissa, &residuals);
    
    printf("Residuals: %e\n", residuals);
    printf("Exponent: %e\tMantissa: %e\n", exponent, mantissa);
    // Get localization length from exponent
    DTYPE loclen = -2 / exponent;

    free(dists);
    free(gfuncsq_vals);
    return(loclen);
}
