#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <argp.h>
#include "constants.h"
#include "params/params.h"
#include "ham_gen/ham_gen.h"
#include "utils/utils.h"

// For NaN detection
// #define _GNU_SOURCE
// #include <fenv.h>

//------------------------------------------------------------------------

/* Function Declarations */ 
int run(struct SystemParams * params, int create_neighbours,
            DTYPE * gfunc);
int output_gfuncsq_matrix(int runs_done, DTYPE * gfuncsq,
                struct SystemParams params, struct OutStream outfiles);
DTYPE post_process(struct SystemParams params, struct OutStream outfiles,
                DTYPE * gfunc);
//------------------------------------------------------------------------

/* Constant Declarations */
const char *argp_program_version =
  "exact_diag_simulation 1.0";
const char *argp_program_bug_address =
  "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
  "exact_diag_simulation -- a simulation of spin-orbit coupled" 
  "2d many-body localized systems.";
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
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};

//------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    // For NaN detection
    // feenableexcept(FE_INVALID | FE_OVERFLOW);

    struct SystemParams params;
    struct OutStream outfiles;

    params_setup(argc, argv, &params, &outfiles, &argp);

    // printf("Exact Diagonalization\n");
    // printf("---------------------\n");
    // printf("len: %d\t nospin: %d\t coupling: %.2e\n"
    //     "disorder: %.2e\t hop_up: %.2e\t hop_dn: %.2e\n",
    //     params.len, params.nospin, params.coupling_const,
    //     params.disorder_strength, params.hop_strength_upup,
    //     params.hop_strength_dndn);

    // if(params.nospin == 1)
    //     params.num_states = params.len*params.width;
    // else
    //     params.num_states = 2*params.len*params.width;

    int ctr;
    int create_neighbours = 1;

    // Seed all random numbers generated by the time 
    srandom((unsigned) time(NULL));

    // Allocate Green's Function Matrix
    // Represents long time limit of the green's function squared.
    DTYPE * gfunc = calloc(params.num_states*params.num_states, sizeof(DTYPE));

    if(gfunc == NULL)
    {
        printf("Could not allocate memory for gfunc! Exiting...\n");
        return(1);
    }

    // printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = 1; ctr <= params.numRuns; ctr++)
    {
        // printf("Run %d started...", ctr);
        // fflush(stdout);
        /* Call run */
        run(&params, create_neighbours, gfunc);
        create_neighbours = 0;
        output_gfuncsq_matrix(ctr, gfunc, params, outfiles);
        // printf("Done\n");
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
int run(struct SystemParams * params, int create_neighbours,
            DTYPE * gfunc)
{
    int num_sites = params->len * params->width;
    int num_states = params->num_states;

    // Create neighbour list if not present
    if (create_neighbours)
    {
        params->neighbours = malloc((num_sites*NEIGHS)*sizeof(int));
        get_neighbour_lists(params->neighbours, params->len, params->width);
    }

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
    utils_get_eigh(ham, num_states, eigvals);

    // Calculate and add green's function long time limit squared
    // printf("Gfunc...");
    // fflush(stdout);
    int degeneracy;
    if(params->nospin == 0)
        degeneracy = DEGEN_EIGVALS;
    else
        degeneracy = NONDEGEN_EIGVALS;
    utils_get_green_func_lim(ham, num_states, gfunc, degeneracy);

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
            fprintf(ofile, "%e ", elem);
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