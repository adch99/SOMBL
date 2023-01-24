#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <argp.h>
#include "constants.h"
#include "params/params.h"
#include "ham_gen/ham_gen.h"
#include "utils/utils.h"
#include "diag/diag.h"
#include "gfunc/gfunc.h"
#include "io/io.h"

// For NaN detection
// #define _GNU_SOURCE
// #include <fenv.h>

//------------------------------------------------------------------------

struct GFunc {
    DTYPE * gfunc_sym[2];
    CDTYPE * gfunc_asym;
};

/* Function Declarations */ 
int run(struct SystemParams * params, struct GFunc gfuncs,
        int run_num, bool isLastRun);
int output_gfuncsq_matrix(int runs_done, DTYPE * gfuncsq,
                struct SystemParams params, struct OutStream outfiles);
DTYPE post_process(struct SystemParams params, struct OutStream outfiles,
                DTYPE * gfunc);
int process_gr_grstar(struct SystemParams * params, struct GFunc gfuncs, 
                    CDTYPE * eigvecs, uint alpha, uint beta, int run_num,
                    bool isLastRun);

//------------------------------------------------------------------------

/* Constant Declarations */
const char *argp_program_version =
    "keldysh_window_batch 1.0";
const char *argp_program_bug_address =
    "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
    "keldysh_window_batch_nobins -- a simulation of spin-orbit coupled" 
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
    {"start",     'f', "STARTFROM",  0, "Run number to start from.",              0},
    { 0 }
};
// Our argp parser.
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};

//  --------------------
// | Relevant Constants |
//  --------------------
// #define ENERGYCOREMIN (-1.0)
// #define ENERGYCOREMAX (1.0)
// #define NUMBINS (50)

//------------------------------------------------------------------------


int main(int argc, char ** argv)
{
    // For NaN detection
    // feenableexcept(FE_INVALID | FE_OVERFLOW);

    struct SystemParams params;
    struct OutStream outfiles;
    printf("Hello!\n");
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


    // Create neighbour list
    // ---------------------
    params.neighbours = malloc((params.num_sites*NEIGHS)*sizeof(int));
    get_neighbour_lists(params.neighbours, params.len, params.width);



    // Set the iteration conditions for the runs
    if (params.batch_num == -1 || params.batch_size == -1)
    {
        start = params.startfrom;
        stop = params.numRuns;
    }
    else
    {
        start = params.startfrom + params.batch_size * (params.batch_num - 1);
        stop = params.batch_size * params.batch_num;
    }

    // Initialize the green functions
    int Lsq = params.len * params.len;
    struct GFunc gfuncs;
    gfuncs.gfunc_sym[0] = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
    gfuncs.gfunc_sym[1] = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
    gfuncs.gfunc_asym = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));

    printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = start; ctr <= stop; ctr++)
    {
        printf("Run %d started...\n", ctr);
        fflush(stdout);
        run(&params, gfuncs, ctr - start + 1, (ctr == stop));
        printf("Run %d done\n", ctr);
    }

    // DTYPE avg_loc_len = post_process(params, outfiles, gfunc);
    // printf("Disorder Averaged Loc Len: %e\n", avg_loc_len);

    free(params.neighbours);
    params_cleanup(&params, &outfiles);
    return(0);
}

/*
    Creates a hamiltonian based on params, calculates the
    long time green's function squared and ADDS it to gfunc.
    Please ensure gfunc is initialized properly for your
    purpose.
*/
int run(struct SystemParams * params, struct GFunc gfuncs,
        int run_num, bool isLastRun)
{
    // int num_sites = params->num_sites;
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

    // Append eigenvalues to file
    char filename[100];
    params_basefilename(*params, "data/mbl", "_eigvals.dat",
                    filename);
    io_append_array('R', eigvals, num_states, filename);

    // For alpha = beta cases
    uint alpha, beta;
    for(alpha = 0; alpha < 2; alpha++)
    {
        process_gr_grstar(params, gfuncs, ham, alpha, alpha,
                        run_num, isLastRun);
    }

    // Now for the alpha != beta case
    // We only need one case
    alpha = 0;
    beta = 1;
    process_gr_grstar(params, gfuncs, ham, alpha, beta,
                    run_num, isLastRun);

    free(eigvals);
    free(ham);
    return(0);
}

int process_gr_grstar(struct SystemParams * params, struct GFunc gfuncs, 
                    CDTYPE * eigvecs, uint alpha, uint beta, int run_num,
                    bool isLastRun)
{
    int L = params->len;
    int Lsq = L*L;
    char filename[100];
    int nmin, nmax, i;
    bool isCheckpoint, isFirst, startMiddle;
    int runs_done = run_num - params->startfrom;

    // Now for the full set
    nmin = 0;
    nmax = 2*Lsq;
    printf("    alpha:%d beta:%d nmin:%d nmax:%d\n", alpha,
            beta, nmin, nmax);
    // Get the filename
    params_gr_grstar_filename(filename, *params, -1, alpha, beta);

    if(alpha == beta)
    {
        DTYPE * gfunc_sym = gfuncs.gfunc_sym[alpha];
        isFirst = (run_num == params->startfrom);
        startMiddle = (params->startfrom > 1);
        if(isFirst && startMiddle)
        {
            // Get the matrix from the file
            io_read_array('R', 'C', gfunc_sym, Lsq, 2*Lsq, filename);

            // Multiply matrix by the number of runs already done
            for(i = 0; i < Lsq * 2*Lsq; i++)
                *(gfunc_sym + i) *= run_num - 1;
        }
        gfuncsq_sym_GR_GRstar_nondeg(eigvecs, gfunc_sym, L,
                                    nmin, nmax, alpha);
        gfuncsq_sym_GR_GRstar_deg(eigvecs, gfunc_sym, L,
                                nmin, nmax, alpha);

        isCheckpoint = ((runs_done % 5 == 0) && (runs_done > 0));
        if(isCheckpoint || isLastRun)
        {
            // Divide matrix by the number of runs already done
            for(i = 0; i < Lsq * 2*Lsq; i++)
                *(gfunc_sym + i) /= run_num;
            
            io_write_array_bin('R', gfunc_sym, Lsq, 2*Lsq, filename);

            // Multiply matrix by the number of runs already done
            for(i = 0; i < Lsq * 2*Lsq; i++)
                *(gfunc_sym + i) *= (CDTYPE) run_num;

        }


        // free(gfunc_sym);
    }
    else // alpha != beta
    {
        CDTYPE * gfunc_asym = gfuncs.gfunc_asym;
        isFirst = (run_num == params->startfrom);
        startMiddle = (params->startfrom > 1);
        if(isFirst && startMiddle)
        {
            // Get the matrix from the file
            io_read_array('C', 'C', gfunc_asym, Lsq, 2*Lsq, filename);

            // Multiply matrix by the number of runs already done
            for(i = 0; i < Lsq * 2*Lsq; i++)
                *(gfunc_asym + i) *= run_num - 1;
        }
        gfuncsq_asym_GR_GRstar_nondeg(eigvecs, gfunc_asym, L,
                                    nmin, nmax, alpha, beta);
        gfuncsq_asym_GR_GRstar_deg(eigvecs, gfunc_asym, L,
                                    nmin, nmax, alpha, beta);
        
        isCheckpoint = ((runs_done % 5 == 0) && (runs_done > 0));
        if(isCheckpoint || isLastRun)
        {
            // Divide matrix by the number of runs already done
            for(i = 0; i < Lsq * 2*Lsq; i++)
                *(gfunc_asym + i) /= (CDTYPE) run_num;

            // Write the new matrix back to the file.
            io_write_array_bin('C', gfunc_asym, Lsq, 2*Lsq, filename);

            // Multiply matrix by the number of runs already done
            for(i = 0; i < Lsq * 2*Lsq; i++)
                *(gfunc_asym + i) *= (CDTYPE) run_num;

        }
        // free(gfunc_asym);
    }
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
        }
        fprintf(ofile, "\n");
    }

    fclose(ofile);
    return 0;
}
