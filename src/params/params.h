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
    DTYPE hop_strength;
    int (*neighbours)[NEIGHS];
    int numRuns; // Technically not a system parameter but it's convenient
    // to have it here.
    int nospin;
    int num_states;
};


struct OutStream {
    char gfuncsq[80];
    char dist_vs_gfuncsq[80];
};

error_t params_parse_opt(int key, char *arg, struct argp_state *state);
struct OutStream params_set_up_datastream(struct SystemParams params);

#endif