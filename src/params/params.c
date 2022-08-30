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
            params->hop_strength_upup = atof(arg);
            params->hop_strength_dndn = atof(arg);
            break;
        case 'u':
            params->hop_strength_upup = atof(arg);
            break;
        case 'd':
            params->hop_strength_dndn = atof(arg);
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

int params_setup(int argc, char ** argv, struct SystemParams * params,
        struct OutStream * outfiles, struct argp * argp)
{
    /* Default Values of Parameters */
    params->len = params->width = 20;
    params->coupling_const = 0;
    params->disorder_strength = 10;
    params->hop_strength_upup = 1;
    params->hop_strength_dndn = 1;
    params->numRuns = 100;
    params->nospin = 0;

    /* Parse our arguments; every option seen by parse_opt will
        be reflected in params. */
    argp_parse(argp, argc, argv, 0, 0, params);
    if(params->nospin)
        params->num_states = params->len*params->width;
    else
        params->num_states = 2*params->len*params->width;

    // Set up the data files for the system params.
    params_set_up_datastream(*params, outfiles);

    printf("Filenames:\n");
    printf("gfuncsq: %s\n", outfiles->gfuncsq);
    printf("dist_vs_gfuncsq: %s\n", outfiles->dist_vs_gfuncsq);
    int i;
    for(i = 0; i < 4; i++)
        printf("dist_vs_gfuncsq_spin[%d]: %s\n", i, outfiles->dist_vs_gfuncsq_spin[i]);

    return 0;
}

int params_set_up_datastream(struct SystemParams params, struct OutStream * outfiles)
{
    char base[16];
    char basename[64];
    // int gfuncsq_check;
    // int dist_vs_gfuncsq_check;

    if (params.nospin)
      strcpy(base, "data/mbl_nospin");
    else
      strcpy(base, "data/mbl");

    sprintf(basename, "%s_%dx%d_W%.4g_C%.4g_TU%.4g_TD%.4g_N%d", base, params.len,
            params.width, params.disorder_strength, params.coupling_const,
            params.hop_strength_upup, params.hop_strength_dndn, params.numRuns);
    sprintf(outfiles->gfuncsq, "%s_greenfuncsq.dat", basename);
    sprintf(outfiles->dist_vs_gfuncsq_spin[0], "%s_upup_distvsgfsq.dat", basename);
    sprintf(outfiles->dist_vs_gfuncsq_spin[1], "%s_updn_distvsgfsq.dat", basename);
    sprintf(outfiles->dist_vs_gfuncsq_spin[2], "%s_dnup_distvsgfsq.dat", basename);
    sprintf(outfiles->dist_vs_gfuncsq_spin[3], "%s_dndn_distvsgfsq.dat", basename);


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
