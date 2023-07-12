#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "constants.h"
#include "gfunc/gfunc.h"
#include "utils/utils.h"

int main()
{
    srandom((unsigned) time(NULL));
    int L = 20;
    int Lsq = L*L;
    CDTYPE * eigvecs = calloc(2*Lsq * 2*Lsq, sizeof(CDTYPE));
    int i, j;
    DTYPE norm;
    CDTYPE elem;


    for(i = 0; i < 2*Lsq; i++)
    {
        for(j = 0; j < 2*Lsq; j++)
        {
            *(eigvecs + RTC(i, j, 2*Lsq)) = ((double) random()) / ((double) RAND_MAX) - 0.5;
            *(eigvecs + RTC(i, j, 2*Lsq)) += I*(((double) random()) / ((double) RAND_MAX) - 0.5);
        }
    }
    for(i = 0; i < 2*Lsq; i++)
    {
        norm = 0;
        for(j = 0; j < 2*Lsq; j++)
        {
            elem = *(eigvecs + RTC(i, j, 2*Lsq));
            norm += crealf(elem)*crealf(elem) + cimagf(elem)*cimagf(elem);
        }
        for(j = 0; j < 2*Lsq; j++)
            *(eigvecs + RTC(i, j, 2*Lsq)) /= sqrt(norm);
    }

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            elem = *(eigvecs + RTC(i, j, 2*Lsq));
            printf("%10.3e%+10.3ej ", crealf(elem), cimagf(elem));
        }
        printf("\n");
    }

    int nmin = 0;
    int nmax = 2*Lsq;
    uint alpha = 0;
    uint beta = 0;

    DTYPE * gfunc_sym = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
    time_t std_time = time(NULL);
    gfuncsq_sym_GR_GRstar_nondeg(eigvecs, gfunc_sym, L,
                                nmin, nmax, alpha);
    gfuncsq_sym_GR_GRstar_deg(eigvecs, gfunc_sym, L,
                            nmin, nmax, alpha);
    std_time = time(NULL) - std_time;
    printf("Time taken by standard: %ld\n", std_time);

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            printf("%-10.3e ", *(gfunc_sym + RTC(i, j, Lsq)));
        }
        printf("\n");
    }

    free(gfunc_sym);

    CDTYPE * gfunc_direct = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
    size_t direct_time = time(NULL);
    gfunc_direct_full(eigvecs, gfunc_direct, L, nmin, nmax, alpha, beta);
    direct_time = time(NULL) - direct_time;
    printf("Time taken by direct: %ld\n", direct_time);

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            elem = *(gfunc_direct + RTC(i, j, 2*Lsq));
            printf("%-10.3e+i%-10.3e ", crealf(elem), cimagf(elem));
        }
        printf("\n");
    }


    free(gfunc_direct);

    nmin = 0;
    nmax = 2*Lsq;
    alpha = 0;
    beta = 1;

    CDTYPE * gfunc_asym = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
    std_time = time(NULL);
    gfuncsq_asym_GR_GRstar_nondeg(eigvecs, gfunc_asym, L,
                                nmin, nmax, alpha, beta);
    gfuncsq_asym_GR_GRstar_deg(eigvecs, gfunc_asym, L,
                            nmin, nmax, alpha, beta);
    std_time = time(NULL) - std_time;
    printf("Time taken by standard: %ld\n", std_time);

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            elem =  *(gfunc_asym + RTC(i, j, Lsq));
            printf("%10.3e%+10.3ej ", crealf(elem), cimagf(elem));
        }
        printf("\n");
    }

    free(gfunc_asym);

    gfunc_direct =  calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
    direct_time = time(NULL);
    gfunc_direct_full(eigvecs, gfunc_direct, L, nmin, nmax, alpha, beta);
    direct_time = time(NULL) - direct_time;
    printf("Time taken by direct: %ld\n", direct_time);

    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            elem = *(gfunc_direct + RTC(i, j, 2*Lsq));
            printf("%10.3e%+10.3ej ", crealf(elem), cimagf(elem));
        }
        printf("\n");
    }


    free(gfunc_direct);



}