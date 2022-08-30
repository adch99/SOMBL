#include <stdlib.h>
#include <stdio.h>
#include <argp.h>
#include "constants.h"
#include "params/params.h"
#include "utils/utils.h"
#include "io/io.h"

/* Constant Declarations */
const char *argp_program_version =
  "calculate_imbalance 1.0";
const char *argp_program_bug_address =
  "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
  "Calculate imbalance -- calculates the imbalance "
  "based on initial conditions and green function "
  "squared data from input files";
// A description of the arguments we accept.
static char args_doc[] = "-s <size> -c <coupling_const>"
                        " -w <disorder_strength>"
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


int main(int argc, char ** argv)
{
    // Get SystemParams from defaults or cmd line args
    struct SystemParams params;
    struct OutStream outfiles;

    params_setup(argc, argv, &params, &outfiles, &argp);

    printf("Calculate Imbalance\n");
    printf("-------------------\n");
    printf("len: %d\t nospin: %d\t coupling: %.2e\n"
        "disorder: %.2e\t hop_up: %.2e\t hop_dn: %.2e\n",
        params.len, params.nospin, params.coupling_const,
        params.disorder_strength, params.hop_strength_upup,
        params.hop_strength_dndn);


    int * occupied_set_up, * occupied_set_dn;
    int set_length_up, set_length_dn;
    DTYPE * gfuncsq = malloc(params.num_states*params.num_states * sizeof(DTYPE));
    DTYPE charge_imb;

    io_get_initial_condition(&occupied_set_up, &set_length_up,
                        &occupied_set_dn, &set_length_dn,
                        "data/alt_up_empty_L3.dat");

    int i;
    printf("Up spins: %d\n", set_length_up);
    for(i = 0; i < set_length_up; i++)
        printf("%d ", *(occupied_set_up + i));
    printf("\n");

    printf("Dn spins: %d\n", set_length_dn);
    for(i = 0; i < set_length_dn; i++)
        printf("%d ", *(occupied_set_dn + i));
    printf("\n");

    io_get_gfuncsq_from_file(gfuncsq, outfiles, params);

    charge_imb = utils_get_charge_imbalance(gfuncsq, occupied_set_up,
                        set_length_up, occupied_set_dn, set_length_dn,
                        params.num_states);

    printf("Charge Imbalance: %e\n", charge_imb);
    free(occupied_set_up);
    free(occupied_set_dn);
    return 0;
}