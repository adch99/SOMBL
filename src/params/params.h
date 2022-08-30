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
};


struct OutStream {
    char gfuncsq[90];
    char dist_vs_gfuncsq[90];
    char dist_vs_gfuncsq_spin[4][100];
};

error_t params_parse_opt(int key, char *arg, struct argp_state *state);
int params_setup(int argc, char ** argv, struct SystemParams * params,
        struct OutStream * outfiles, struct argp * argp);
int params_set_up_datastream(struct SystemParams params,
                            struct OutStream * outfiles);
#endif