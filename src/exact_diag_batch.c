#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <argp.h>
#include "constants.h"
#include "params/params.h"
#include "ham_gen/ham_gen.h"
#include "utils/utils.h"
#include "diag/diag.h"

// For NaN detection
// #define _GNU_SOURCE
// #include <fenv.h>

//------------------------------------------------------------------------

/* Function Declarations */ 
int run(struct SystemParams * params, DTYPE * gfunc);
int output_gfuncsq_matrix(int runs_done, DTYPE * gfuncsq,
                struct SystemParams params, struct OutStream outfiles);
DTYPE post_process(struct SystemParams params, struct OutStream outfiles,
                DTYPE * gfunc);
//------------------------------------------------------------------------

/* Constant Declarations */
const char *argp_program_version =
    "exact_diag_batch 1.0";
const char *argp_program_bug_address =
    "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
    "exact_diag_batch -- a simulation of spin-orbit coupled" 
    "2d many-body localized systems.";
// A description of the arguments we accept.
static char args_doc[] = "-s <size> -c <coupling_const>"
                        "-w <disorder_strength>"
                        " -t <hop_strength> -n <num_runs>"
                        " [-p]";
// The options we understand.
static struct argp_option options[] = {
    {"size",      's', "SIZE",      0, "Length and width of the lattice",        0},
    {"coupling",  'c', "COUPLING",  0, "Spin-orbit coupling constant",           0},
    {"disorder",  'w', "DISORDER",  0, "Strength of the disorder",               0},
    {"hopping",   't', "HOPPING",   0, "Strength of the hopping",                0},
    {"hopup",     'u', "HOPUP",     0, "Strength of the hopping for up spins",   0},
    {"hopdn",     'd', "HOPDN",     0, "Strength of the hopping for down spins", 0},
    {"runs",      'n', "NUMRUNS",   0, "Number of runs in the disorder average", 0},
    {"nospin",    'p', 0,           0, "Use a spinless model hamiltonian.",      0},
    {"batch",     'b', "BATCHNUM",  0, "Batch number of this set of runs.",      0},
    {"batchsize", 'j', "BATCHSIZE", 0, "Number of runs in each batch.",          0},
    {"seed",      'x', "SEED",      0, "Seed for random number generation.",     0},
    { 0 }
};
// Our argp parser.
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};

//------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    // For NaN detection
    // feenableexcept(FE_INVALID | FE_OVERFLOW);

    struct SystemParams params;
    struct OutStream outfiles;

    params_setup(argc, argv, &params, &outfiles, &argp, 0);

    int ctr, start, stop;

    // Initialize Random Seed
    // ----------------------
    // If seed not given, use current time to generate
    // the seed. Otherwise, use the given seed.
    if (params.seed == -1)
        srandom((unsigned) time(NULL));
    else
        srandom((unsigned int) params.seed);


    // Allocate Green's Function Matrix
    // --------------------------------
    // Represents long time limit of the green's function squared.
    DTYPE * gfunc = calloc(params.num_states*params.num_states, sizeof(DTYPE));
    if(gfunc == NULL)
    {
        printf("Could not allocate memory for gfunc! Exiting...\n");
        return(1);
    }

    // Create neighbour list
    // ---------------------
    params.neighbours = malloc((params.num_sites*NEIGHS)*sizeof(int));
    get_neighbour_lists(params.neighbours, params.len, params.width);



    // Set the iteration conditions for the runs
    if (params.batch_num == -1 || params.batch_size == -1)
    {
        start = 1;
        stop = params.numRuns;
    }
    else
    {
        start = 1 + params.batch_size * (params.batch_num - 1);
        stop = params.batch_size * params.batch_num;
    }

    printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = start; ctr <= stop; ctr++)
    {
        printf("Run %d started...", ctr);
        fflush(stdout);
        run(&params, gfunc);
        output_gfuncsq_matrix(ctr - start + 1, gfunc, params, outfiles);
        printf("Done\n");
    }

    // DTYPE avg_loc_len = post_process(params, outfiles, gfunc);
    // printf("Disorder Averaged Loc Len: %e\n", avg_loc_len);

    free(params.neighbours);
    params_cleanup(&params, &outfiles);
    free(gfunc);
    return(0);
}

/*
    Creates a hamiltonian based on params, calculates the
    long time green's function squared and ADDS it to gfunc.
    Please ensure gfunc is initialized properly for your
    purpose.
*/
int run(struct SystemParams * params, DTYPE * gfunc)
{
    int num_sites = params->num_sites;
    int num_states = params->num_states;

    // Create hamiltonian
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));

    if(ham == NULL)
    {
        printf("Could not allocate memory for ham! Exiting...\n");
        exit(1);
    }
    
    // printf("Creating Ham...");
    // fflush(stdout);
    if(params->nospin == 1)
        hamiltonian_nospin(ham, params->len, params->width,
                params->disorder_strength, params->hop_strength_upup,
                params->neighbours);
    else
        hamiltonian(ham, params->len, params->width, params->coupling_const,
                params->disorder_strength, params->hop_strength_upup,
                params->hop_strength_dndn, params->neighbours);


    // Calculate eigenvectors
    // printf("Eigh...");
    // fflush(stdout);
    DTYPE * eigvals = calloc(num_states, sizeof(DTYPE));
    diag_get_eigh(ham, num_states, eigvals);

    // Calculate and add green's function long time limit squared
    // printf("Gfunc...");
    // fflush(stdout);
    int degeneracy;
    if(params->nospin == 0)
        degeneracy = DEGEN_EIGVALS;
    else
        degeneracy = NONDEGEN_EIGVALS;
    utils_get_green_func_lim(ham, num_states, gfunc, degeneracy);
    // printf("Wrapping up...");
    // fflush(stdout);

    free(eigvals);
    free(ham);
    return(0);
}


/*
    Outputs the gfuncsq matrix to the corresponding file in outfiles.
*/
int output_gfuncsq_matrix(int runs_done, DTYPE * gfuncsq,
                struct SystemParams params, struct OutStream outfiles)
{
    int i, j;
    DTYPE elem;
    FILE * ofile = fopen(outfiles.gfuncsq, "w");
    if (ofile == NULL)
    {
        printf("Cannot open file %s!\n", outfiles.gfuncsq);
        return(-1);
    }
    // Divide the gfunc matrix to get disorder avged Green's function
    for(i = 0; i < params.num_states; i++)
    {
        for(j = 0; j < params.num_states; j++)
        {
            elem = *(gfuncsq + RTC(i, j, params.num_states)) / (DTYPE) runs_done;
            fprintf(ofile, "%le ", elem);
            // if(isnan(elem))
            // {
            //     printf("We have a problem!!\n");
            //     printf("NaN detected at i=%d in gfunc", i);
            //     return(-1);
            // }
        }
        fprintf(ofile, "\n");
    }

    fclose(ofile);
    return 0;
}