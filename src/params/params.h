#ifndef MBL_PARAMS_H
#define MBL_PARAMS_H

#include <argp.h>
#include "../constants.h"

/* Struct Declarations */
struct SystemParams {
    int len;
    int width;
    DTYPE coupling_const;
    DTYPE disorder_strength;
    DTYPE hop_strength_upup;
    DTYPE hop_strength_dndn;
    int (*neighbours)[NEIGHS];
    int numRuns; // Technically not a system parameter but it's convenient
    // to have it here.
    int nospin;
    int num_states;
    int num_sites;
    int batch_num;
    int batch_size;
    int seed;
};


struct OutStream {
    char * gfuncsq;
    char * dist_vs_gfuncsq;
    char * dist_vs_gfuncsq_spin[4];
};

error_t params_parse_opt(int key, char *arg, struct argp_state *state);
int params_setup(int argc, char ** argv, struct SystemParams * params,
        struct OutStream * outfiles, struct argp * argp, int sigma);
int params_set_up_datastream(struct SystemParams params,
                            struct OutStream * outfiles, int sigma);
int params_cleanup(struct SystemParams * params, struct OutStream * outfiles);
#endif