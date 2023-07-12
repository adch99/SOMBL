#include <stdio.h>
#include <stdlib.h>
#include "utils/utils.h"
#include "params/params.h"
#include "io/io.h"

/* Constant Declarations */
const char *argp_program_version =
    "batch_average 1.0";
const char *argp_program_bug_address =
    "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
    "batch_average -- a simulation of spin-orbit coupled" 
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
    struct SystemParams params;
    struct OutStream outfiles;
    int size, num_batches, batchnum;
    char * ofilename, *ifilename;
    int i, j;

    params_setup(argc, argv, &params, &outfiles, &argp, 0);
    size = params.num_states;
    num_batches = params.numRuns / params.batch_size;
    DTYPE * avg_gfuncsq = calloc(size*size, sizeof(DTYPE));
    DTYPE * batch_gfuncsq = calloc(size*size, sizeof(DTYPE));
    ofilename = outfiles.gfuncsq;

    // For each batch
    for(batchnum = 1; batchnum <= num_batches; batchnum++)
    {
        // Generate the correct filename
        params.batch_num = batchnum;
        params_set_up_datastream(params, &outfiles, 0);
        ifilename = outfiles.gfuncsq;
        // Read the array and add it to the sum
        printf("Reading file %s...\n", ifilename);
        io_read_array('R', 'C', batch_gfuncsq, size, size, ifilename);
        utils_add_to_matrix_real(avg_gfuncsq, batch_gfuncsq,
                                size, size);
    }

    // Divide to get the average
    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            *(avg_gfuncsq + RTC(i,j,size)) /= num_batches;
        }
    }

    // Save this array in a file
    printf("Writing to file %s...\n", ofilename);
    io_write_array('R', 'C', avg_gfuncsq, size, size, ofilename);

    return(0);
}


