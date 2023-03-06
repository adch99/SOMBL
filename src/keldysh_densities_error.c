#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include <string.h>
#include "constants.h"
#include "params/params.h"
#include "utils/utils.h"
#include "io/io.h"
#include "mkl.h"
//#include "cblas.h"
//#include "lapacke.h"

/* Constant Declarations */
const char *argp_program_version =
    "keldysh_window_batch 1.0";
const char *argp_program_bug_address =
    "<aditya.chincholi@students.iiserpune.ac.in>";
// Program documentation.
static char doc[] =
    "exact_diag_batch -- a simulation of spin-orbit coupled" 
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
    { 0 }
};
// Our argp parser.
static struct argp argp = { options, params_parse_opt, args_doc, doc, 0, 0, 0};

// --------------------------
// Program specific constants
// --------------------------
enum SpinPair {UPUP=0, UPDN, DNUP, DNDN};
enum Spin {UP=0, DOWN};

// ----------------------
// Function Declarations
// ----------------------
int get_initial_cond(void * initial_cond, char type, char * filename);
int get_density(void * density, void * initial_cond,
                struct SystemParams params, uint alpha, uint beta, int bin,
                void * variance);
int save_density(void * density, struct SystemParams params,
                uint alpha, uint beta, int bin, char * basename,
                void * variance);
int create_initial_vector(void * initial_cond, char type, int num_states);

int main(int argc, char ** argv)
{
    struct SystemParams params;
    struct OutStream outfiles;
    uint alpha, beta;
    int bin;
    char filename[256];
    char basename[] = "alt_up_down";
    params_setup(argc, argv, &params, &outfiles, &argp, 0);
    params.energybins = 0;
    snprintf(filename, 256, "data/%s_L%d.dat", basename, params.len);
    printf("Initial Cond File: %s\n", filename);

    // bin = -1 gives the full matrix.
    for(bin = -1; bin < params.energybins; bin++)
    {
        for(alpha = 0; alpha < 2; alpha++)
        {
            for(beta = 0; beta < 2; beta++)
            {
                if (alpha == beta)
                {
                    DTYPE * initial_cond = calloc(params.num_states, sizeof(DTYPE));
                    DTYPE * density = calloc(params.num_sites, sizeof(DTYPE));
                    DTYPE * variance = calloc(params.num_sites, sizeof(DTYPE));
                    io_get_initial_cond_vector(initial_cond, 'R', filename);
                    utils_create_keldysh_vector(initial_cond, 'R', params.num_states);
                    get_density(density, initial_cond, params, alpha, beta, bin, variance);
                    save_density(density, params, alpha, beta, bin, basename, variance);
                    free(variance);
                    free(initial_cond);
                    free(density);        
                }
                else if (alpha != beta)
                {
                    CDTYPE * initial_cond = calloc(params.num_states, sizeof(CDTYPE));
                    CDTYPE * density = calloc(params.num_sites, sizeof(CDTYPE));
                    CDTYPE * variance = calloc(params.num_sites, sizeof(CDTYPE));
                    io_get_initial_cond_vector(initial_cond, 'C', filename);
                    utils_create_keldysh_vector(initial_cond, 'C', params.num_states);
                    get_density(density, initial_cond, params, alpha, beta, bin, variance);
                    save_density(density, params, alpha, beta, bin, basename, variance);
                    free(variance);
                    free(initial_cond);
                    free(density);        
                }
            }
        }
    }

    return(0);
}



int get_density(void * density, void * initial_cond,
                struct SystemParams params, uint alpha,
                uint beta, int bin, void * variance)
{
    if(alpha > 2 || beta > 2)
    {
        fprintf(stderr, "α and β should be either 0 (UP) or 1 (DOWN).\n");
        exit(EXIT_FAILURE);
    }

    int L = params.len;
    int Lsq = L*L;

    if(alpha == beta)
    {
        DTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
        DTYPE * gfuncvar = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
        char filename[64];
        char varfilename[128];
        params_gr_grstar_filename(filename, params, bin, alpha, beta);
        sprintf(varfilename, "%s.variance", filename);
        io_read_array_bin('R', gfuncsq, Lsq, 2*Lsq, filename);
        io_read_array_bin('R', gfuncvar, Lsq, 2*Lsq, varfilename);
        cblas_sgemv(CblasColMajor, CblasNoTrans, Lsq, 2*Lsq,
                    1.0, gfuncsq, Lsq, (DTYPE *) initial_cond, 1,
                    0.0, (DTYPE *) density, 1);
        cblas_sgemv(CblasColMajor, CblasNoTrans, Lsq, 2*Lsq,
                    1.0/(double)(2*Lsq), gfuncvar, Lsq, (DTYPE *) initial_cond, 1,
                    0.0, (DTYPE *) variance, 1);
        
        free(gfuncsq);
    }
    else if(alpha == 0 && beta == 1)
    {
        CDTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
        CDTYPE * gfuncvar = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
        char filename[64];
        char varfilename[128];
        params_gr_grstar_filename(filename, params, bin, alpha, beta);
        sprintf(varfilename, "%s.variance", filename);
        io_read_array_bin('C', gfuncsq, Lsq, 2*Lsq, filename);
        io_read_array_bin('C', gfuncvar, Lsq, 2*Lsq, varfilename);
        lapack_complex_float blas_alpha, blas_beta;
        blas_alpha = 1;
        blas_beta = 0;
        // blas_alpha.real = 1;
        // blas_alpha.imag = 0;
        // blas_beta.real = 0;
        // blas_beta.imag = 0;
        cblas_cgemv(CblasColMajor, CblasNoTrans, Lsq, 2*Lsq,
                    &blas_alpha, gfuncsq, Lsq, (CDTYPE *) initial_cond, 1,
                    &blas_beta, (CDTYPE *) density, 1);
        blas_alpha = 1/(double)(2*Lsq);
        blas_beta = 0;
        cblas_cgemv(CblasColMajor, CblasNoTrans, Lsq, 2*Lsq,
                    &blas_alpha, gfuncvar, Lsq, (CDTYPE *) initial_cond, 1,
                    &blas_beta, (CDTYPE *) variance, 1);
        free(gfuncsq);
    }
    else if(alpha == 1 && beta == 0)
    {
        // We need the conjugate value here.
        // G_R(i, α; k, γ)G_R^*(i, β; k, γ) = [G_R(i, β; k, γ)G_R^*(i, α; k, γ)]^*
        CDTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
        CDTYPE * gfuncvar = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
        char filename[128];
        char varfilename[256];
        params_gr_grstar_filename(filename, params, bin, beta, alpha);
        sprintf(varfilename, "%s.variance", filename);
        io_read_array_bin('C', gfuncsq, Lsq, 2*Lsq, filename);
        io_read_array_bin('C', gfuncvar, Lsq, 2*Lsq, varfilename);
        // Take the conjugate of the initial_cond
        int i;
        for(i = 0; i < 2*Lsq; i++)
            *((CDTYPE *)initial_cond + i) = conj(*((CDTYPE *)initial_cond + i));
        lapack_complex_float blas_alpha, blas_beta;
        blas_alpha = 1;
        blas_beta = 0;
        // blas_alpha.real = 1;
        // blas_alpha.imag = 0;
        // blas_beta.real = 0;
        // blas_beta.imag = 0;
        cblas_cgemv(CblasColMajor, CblasNoTrans, Lsq, 2*Lsq,
                    &blas_alpha, gfuncsq, Lsq, (CDTYPE *) initial_cond, 1,
                    &blas_beta, (CDTYPE *) density, 1);
        blas_alpha = 1/(double)(2*Lsq);
        blas_beta = 0;
        cblas_cgemv(CblasColMajor, CblasNoTrans, Lsq, 2*Lsq,
                    &blas_alpha, gfuncvar, Lsq, (CDTYPE *) initial_cond, 1,
                    &blas_beta, (CDTYPE *) variance, 1);
        // Take the conjugate of the output
        for(i = 0; i < Lsq; i++)
            *((CDTYPE *)density + i) = conj(*((CDTYPE *)density + i));
        free(gfuncsq);
    }

    return(0);
}

int save_density(void * density, struct SystemParams params,
                uint alpha, uint beta, int bin, char * basename,
                void * variance)
{
    char filename[128];
    char varfilename[256];
    char suffix[64];
    int i;
    if(bin >= 0)
        snprintf(suffix, 64, "_a%d_b%d_bin%d_%s.dat", alpha, beta, bin, basename);
    else
        snprintf(suffix, 64, "_a%d_b%d_full_%s.dat", alpha, beta, basename);
    params_basefilename(params, "data/mbl_density", suffix, filename);

    sprintf(varfilename, "%s.variance", filename);

    printf("Writing density alpha=%d beta=%d bin=%d to %s\n",
            alpha, beta, bin, filename);
    printf("Writing variance alpha=%d beta=%d bin=%d to %s\n",
            alpha, beta, bin, varfilename);


    FILE * ifile = io_safely_open('w', filename);
    for(i = 0; i < params.num_sites; i++)
    {
        if(alpha == beta)
            fprintf(ifile, "%e ", *((DTYPE *)density + i));
        else
        {
            CDTYPE elemc = *((CDTYPE *) density + i);
            fprintf(ifile, " (%e+%ej) ", crealf(elemc), cimagf(elemc));
        }
    }
    fclose(ifile);

    ifile = io_safely_open('w', varfilename);
    for(i = 0; i < params.num_sites; i++)
    {
        if(alpha == beta)
            fprintf(ifile, "%e ", *((DTYPE *)variance + i));
        else
        {
            CDTYPE elemc = *((CDTYPE *) variance + i);
            fprintf(ifile, " (%e+%ej) ", crealf(elemc), cimagf(elemc));
        }
    }
    fclose(ifile);


    return(0);
}