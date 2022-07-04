#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <argp.h>
#include "constants.h"
#include "ham_gen/ham_gen.h"
#include "utils/utils.h"

//------------------------------------------------------------------------

/* Struct Declarations */
struct SystemParams {
    int len;
    int width;
    DTYPE coupling_const;
    DTYPE disorder_strength;
    DTYPE hop_strength;
    int (*neighbours)[NEIGHS];
};

struct OutStream {
    FILE * eigvals;
    FILE * loc_lens;
};

//------------------------------------------------------------------------

/* Function Declarations */ 
double run(struct SystemParams * params, int create_neighbours,
            struct OutStream outfiles);
struct OutStream set_up_datastream(struct SystemParams params);
int analysis(DTYPE * eigvals, int num_eigvals, struct SystemParams *params,
            struct OutStream outfiles);
static error_t parse_opt (int key, char *arg, struct argp_state *state);

//------------------------------------------------------------------------

/* Constant Declarations */
const char *argp_program_version =
  "exact_diag_simulation 1.0";
const char *argp_program_bug_address =
  "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
  "exact_diag_simulation -- a simulation of spin-orbit coupled 2d many-body localized systems.";
// A description of the arguments we accept.
static char args_doc[] = "-s <size> -c <coupling_const> -w <disorder_strength> -t <hop_strength>";
// The options we understand.
static struct argp_option options[] = {
  {"size",     's', "SIZE",     0, "Length and width of the lattice" },
  {"coupling", 'c', "COUPLING", 0, "Spin-orbit coupling constant" },
  {"disorder", 'w', "DISORDER", 0, "Strength of the disorder"},
  {"hopping",  't', "HOPPING",  0, "Strength of the hopping"},
  { 0 }
};
// Our argp parser.
static struct argp argp = { options, parse_opt, args_doc, doc };

//------------------------------------------------------------------------

int main(int argc, char ** argv)
{
    struct SystemParams params;
    /* Default Values of Parameters */
    params.len = params.width = 20;
    params.coupling_const = 0;
    params.disorder_strength = 10;
    params.hop_strength = 1;


    /* Parse our arguments; every option seen by parse_opt will
        be reflected in params. */
    argp_parse (&argp, argc, argv, 0, 0, &params);

    printf("len: %d coupling: %.2e\ndisorder: %.2e hopping: %.2e\n",
        params.len, params.coupling_const, params.disorder_strength,
        params.hop_strength);


    int ctr;
    int numRuns = 100;
    struct OutStream outfiles = set_up_datastream(params);

    DTYPE avg_loc_len = 0;
    int create_neighbours = 1;

    printf("Starting Simulation for Exact Diagonalization...\n");
    for(ctr = 1; ctr <= numRuns; ctr++)
    {
        printf("Run %d started...", ctr);
        fflush(stdout);
        /* Call run */
        avg_loc_len += run(&params, create_neighbours, outfiles);
        create_neighbours = 0;
        
        printf("\tDone\n");
    }
    avg_loc_len /= (double) numRuns;
    printf("Disorder Averaged Loc Len: %lf\n", avg_loc_len);

    fclose(outfiles.eigvals);
    fclose(outfiles.loc_lens);
    free(params.neighbours);
    return(0);
}

struct OutStream set_up_datastream(struct SystemParams params)
{
    struct OutStream outfiles;
    char basename[64];
    char eigvalsname[64];
    char loclensname[64];
    sprintf(basename, "data/mbl%dx%d_W%.2g_C%.2g_T%.2g_", params.len, params.width,
            params.disorder_strength, params.coupling_const, params.hop_strength);
    sprintf(eigvalsname, "%seigvals.dat", basename);
    sprintf(loclensname, "%sloclens.dat", basename);

    outfiles.eigvals = fopen(eigvalsname, "w");
    outfiles.loc_lens = fopen(loclensname, "w");
    
    if (outfiles.eigvals == NULL)
    {
        printf("Failed to open %s", eigvalsname);
        exit(1);
    }

    if (outfiles.loc_lens == NULL)
    {
        printf("Failed to open %s", loclensname);
        exit(1);
    }

    return(outfiles);
}

double run(struct SystemParams * params, int create_neighbours, struct OutStream outfiles)
{
    /* Set parameters */
    /* Create hamiltonian */

    DTYPE energy = 1;

    int num_sites = params->len * params->width;
    int num_states = 2*num_sites;

    if (create_neighbours)
    {
        params->neighbours = malloc((num_sites*NEIGHS)*sizeof(int)); 
        get_neighbour_lists(params->neighbours, params->len, params->width);
    }

    CDTYPE * ham = malloc(sizeof(CDTYPE)*(num_states*num_states));
    
    hamiltonian(ham, params->len, params->width, params->coupling_const,
                params->disorder_strength, params->hop_strength,
                params->neighbours);

    /* Calculate eigenvalues */
    // DTYPE * eigvals;
    DTYPE * eigvals = calloc(num_states, sizeof(DTYPE));
    utils_get_eigvalsh(ham, num_states, eigvals);

    // printf("\nmin: %lf\tmax: %lf\n", *eigvals, *(eigvals + num_states-1));

    /* Analysis */
    DTYPE avg_loc_len = analysis(eigvals, num_states, params, outfiles);

    free(eigvals);
    free(ham);

    /* Return the localization lengths */
    return(avg_loc_len);
}

int analysis(DTYPE * eigvals, int num_eigvals, struct SystemParams *params,
            struct OutStream outfiles)
{
    /* Calculate localization lengths */
    int i;
    DTYPE avg_loc_len = 0;
    DTYPE loc_len;
    for(i = 0; i < num_eigvals; i++)
    {
        loc_len = utils_loc_len(-1, eigvals, params->hop_strength,
                                    params->len, i); 
        avg_loc_len += loc_len;
        fprintf(outfiles.eigvals, "%e", *(eigvals + i));
        fprintf(outfiles.loc_lens, "%e", loc_len);

        if (i != num_eigvals-1)
        {
            fprintf(outfiles.eigvals, " ");
            fprintf(outfiles.loc_lens, " ");
        }
    }
    fprintf(outfiles.eigvals, "\n");
    fprintf(outfiles.loc_lens, "\n");
    avg_loc_len /= num_eigvals;
    return(avg_loc_len);
}

/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct SystemParams *params = state->input;

  switch (key)
    {
    case 's':
      params->len = params->width = atoi(arg);
      break;
    case 'c':
      params->coupling_const = atof(arg);
      break;
    case 'w':
      params->disorder_strength = atof(arg);
      break;
    case 't':
      params->hop_strength = atof(arg);
      break;

    // case ARGP_KEY_ARG:
    //   if (state->arg_num >= 2)
    //     /* Too many arguments. */
    //     argp_usage (state);
    //   break;

    // case ARGP_KEY_END:
    //   if (state->arg_num < 2)
    //     /* Not enough arguments. */
    //     argp_usage (state);
    //   break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}



