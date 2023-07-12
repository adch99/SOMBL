#include <stdlib.h>
#include <stdio.h>
#include <argp.h>
#include <math.h>
#include "constants.h"
#include "params/params.h"
#include "utils/utils.h"
#include "io/io.h"

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
    {"hopup",    'u', "HOPUP",    0, "Strength of the hopping for up spins",   0},
    {"hopdn",    'd', "HOPDN",    0, "Strength of the hopping for down spins", 0},
    {"runs",     'n', "NUMRUNS",  0, "Number of runs in the disorder average", 0},
    {"nospin",   'p', 0,          0, "Use a spinless model hamiltonian.",      0},
    { 0 }
};
// Our argp parser.
static struct argp argp = {options, params_parse_opt, args_doc, doc, 0, 0, 0};

int main(int argc, char ** argv)
{
    // Get SystemParams from defaults or cmd line args
    struct SystemParams params;
    struct OutStream outfiles;

    params_setup(argc, argv, &params, &outfiles, &argp, 0);

    printf("Calculate G^2(r)\n");
    printf("----------------\n");
    printf("len: %d\t nospin: %d\t coupling: %.2e\n"
        "disorder: %.2e\t hop_up: %.2e\t hop_dn: %.2e\n",
        params.len, params.nospin, params.coupling_const,
        params.disorder_strength, params.hop_strength_upup,
        params.hop_strength_dndn);


    // Get the gfuncsq matrix from the data file
    DTYPE * matrix = malloc(params.num_states*params.num_states*sizeof(DTYPE));
    io_get_gfuncsq_from_file(matrix, outfiles, params, 0);
    // utils_print_matrix(matrix, params.num_states, params.num_states, 'R', 'F');
    // Bin the gfunsq matrix and create a function of
    // distance G^2(|r_i - r_j|)
    DTYPE * dists;
    DTYPE * gfuncsq;
    DTYPE * gfuncsq_err;
    int data_len;
    // Let's bin the data into bins of width 1.
    // int bins = floor(params.len * sqrt(2) / M_PI);
    // int bins = floor((DTYPE) params.len * sqrt(2));
    int bins = params.len;
    printf("Bins: %d\n", bins);
    // int bins = 30;
    data_len = utils_construct_data_vs_dist(matrix, params.num_states, params.len,
                                        bins, &dists, &gfuncsq, &gfuncsq_err);

    if(params.nospin == 1)
    {
        io_output_function_data(dists, gfuncsq, gfuncsq_err,
                                outfiles.dist_vs_gfuncsq, data_len);
    }        
    else
    {
        int i;
        for(i = 0; i < 4; i++)
        {
            char * fname = outfiles.dist_vs_gfuncsq_spin[i];
            printf("Outputting to %s\n", fname);
            io_output_function_data(dists, gfuncsq + i*bins,
                                    gfuncsq_err + i*bins, fname, data_len);
        }
    }
    params_cleanup(&params, &outfiles);
    free(matrix);
    free(dists);
    free(gfuncsq);
    return 0;
}

