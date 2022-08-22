#include <argp.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "params.h"
#include "../constants.h"


/* Parse a single option. */
error_t params_parse_opt(int key, char *arg, struct argp_state *state)
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
        case 'n':
            params->numRuns = atoi(arg);
            break;
        case 'p':
            params->nospin = 1;
            break;
        default:
            return ARGP_ERR_UNKNOWN;
        }
    return 0;
}

struct OutStream params_set_up_datastream(struct SystemParams params)
{
    struct OutStream outfiles;
    char base[16];
    char basename[64];
    // int gfuncsq_check;
    // int dist_vs_gfuncsq_check;

    if (params.nospin)
      strcpy(base, "data/mbl_nospin");
    else
      strcpy(base, "data/mbl");

    sprintf(basename, "%s_%dx%d_W%.4g_C%.4g_T%.4g_N%d_", base, params.len,
            params.width, params.disorder_strength, params.coupling_const,
            params.hop_strength, params.numRuns);
    sprintf(outfiles.gfuncsq, "%sgreenfuncsq.dat", basename);
    sprintf(outfiles.dist_vs_gfuncsq, "%sdistvsgfsq.dat", basename);

    // gfuncsq_check = access(outfiles.gfuncsq, W_OK);
    // dist_vs_gfuncsq_check = access(outfiles.dist_vs_gfuncsq, W_OK);
    
    // if (gfuncsq_check != 0)
    // {
    //     printf("Cannot open %s!\n", outfiles.gfuncsq);
    //     exit(1);
    // }
    // if (dist_vs_gfuncsq_check != 0)
    // {
    //     printf("Cannot open %s!\n", outfiles.dist_vs_gfuncsq);
    //     exit(1);
    // }

    // fclose(dist_vs_gfuncsq_file);
    // fclose(gfuncsq_file);

    return(outfiles);
}
