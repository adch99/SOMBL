#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "unity.h"
#include "../src/utils/utils.h"
#include "../src/gfunc/gfunc.h"
#include "../src/io/io.h"
#include "../src/constants.h"

#define TOL 1e-6

// For NaN detection
// #define _GNU_SOURCE
// #include <fenv.h>


void setUp(void) {}

void tearDown(void) {}

void assert_equal_complex_double_array(const double complex * expected,
                                        const double complex * actual,
                                        int length)
{
    int i;
    double complex elem1, elem2;
    char template[] = "Element %d %s Part";
    char message[40];
    for(i = 0; i < length; i++)
    {
        elem1 = *(actual + i);
        elem2 = *(expected + i);
        snprintf(message, 40, template, i, "Real");
        TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(creal(elem2), creal(elem1), message);
        snprintf(message, 40, template, i, "Imag");
        TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(cimag(elem2), cimag(elem1), message);

    }
}

void assert_equal_complex_float_array(const float complex * expected,
                                        const float complex * actual,
                                        int length)
{
    int i;
    double complex elem1, elem2;
    char template[] = "Element %d %s Part";
    char message[40];
    for(i = 0; i < length; i++)
    {
        elem1 = *(actual + i);
        elem2 = *(expected + i);
        snprintf(message, 40, template, i, "Real");
        TEST_ASSERT_FLOAT_WITHIN_MESSAGE(TOL, crealf(elem2), crealf(elem1), message);
        snprintf(message, 40, template, i, "Imag");
        TEST_ASSERT_FLOAT_WITHIN_MESSAGE(TOL, cimagf(elem2), cimagf(elem1), message);

    }
}


void test_utils_gfuncsq_sigma_matrix_nondeg()
{
    int L = 4;
    int Lsq = L * L;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *sigma = calloc(2 * 2, sizeof(CDTYPE));
    CDTYPE *eigvecs = calloc(2 * Lsq * 2 * Lsq, sizeof(CDTYPE));
    DTYPE *exp_gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));
    DTYPE *gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));

    *(sigma + RTC(0, 0, 2)) = 1.1 + 0*I;
    *(sigma + RTC(0, 1, 2)) = 0.2 + 0*I;
    *(sigma + RTC(1, 0, 2)) = -0.1 + 0*I;
    *(sigma + RTC(1, 1, 2)) = 0.6 + 0*I;

    // printf("Sigma:\n");
    // for (i = 0; i < 2; i++)
    // {
    //     for (j = 0; j < 2; j++)
    //     {
    //         elem = *(sigma + RTC(i, j, 2));
    //         printf("%.2g+%.2gj  ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *eigfile = fopen("data/fake_eigenvectors_4x4.dat", "r");
    io_zread_2d(eigvecs, 2 * Lsq, 2 * Lsq, eigfile);
    fclose(eigfile);

    // printf("First few values of eigvecs:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         elem = *(eigvecs + RTC(i, j, 2 * Lsq));
    //         printf("(%7.2e+%7.2ej) ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *gfuncfile = fopen("data/nondeg_sigma_gfuncsq_python_4x4.dat", "r");
    io_dread_2d(exp_gfuncsq, Lsq, Lsq, gfuncfile);
    fclose(gfuncfile);

    // printf("First few values of EXPECTED gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(exp_gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    utils_gfuncsq_sigma_matrix_nondeg(gfuncsq, sigma, eigvecs, L);

    // printf("First few values of ACTUAL gfuncsq:\n");
    // for (i = 0; i < Lsq; i++)
    // {
    //     for (j = 0; j < Lsq; j++)
    //     {
    //         printf("%.3g ", *(gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // DTYPE diff = 0;
    // for (i = 0; i < Lsq * Lsq; i++)
    //     diff += fabs(*(gfuncsq + i) - *(exp_gfuncsq + i));

    // printf("Diff = %.2g\n", diff);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(exp_gfuncsq, gfuncsq, Lsq * Lsq);

    // char message[64];
    // CDTYPE actual, expected;
    // for(i = 0; i < Lsq; i++)
    // {
    //     for(j = 0; j < Lsq; j++)
    //     {
    //         actual = *(gfuncsq + RTC(i, j, Lsq));
    //         expected = *(exp_gfuncsq + RTC(i, j, Lsq));
    //         sprintf(message, "Error at (%d,%d) real part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(crealf(actual), crealf(expected), message);
    //         sprintf(message, "Error at (%d,%d) imag part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(cimagf(actual), cimagf(expected), message);
    //     }
    // }

    free(sigma);
    free(eigvecs);
    free(exp_gfuncsq);
    free(gfuncsq);
}

void test_utils_gfuncsq_sigma_matrix_deg()
{
    int L = 4;
    int Lsq = L * L;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *sigma = calloc(2 * 2, sizeof(CDTYPE));
    CDTYPE *eigvecs = calloc(2 * Lsq * 2 * Lsq, sizeof(CDTYPE));
    DTYPE *exp_gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));
    DTYPE *gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));

    *(sigma + RTC(0, 0, 2)) = 1.1 + 0*I;
    *(sigma + RTC(0, 1, 2)) = 0.2 + 0*I;
    *(sigma + RTC(1, 0, 2)) = -0.1 + 0*I;
    *(sigma + RTC(1, 1, 2)) = 0.6 + 0*I;

    // printf("Sigma:\n");
    // for (i = 0; i < 2; i++)
    // {
    //     for (j = 0; j < 2; j++)
    //     {
    //         elem = *(sigma + RTC(i, j, 2));
    //         printf("%.2g+%.2gj  ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *eigfile = fopen("data/fake_eigenvectors_4x4.dat", "r");
    io_zread_2d(eigvecs, 2 * Lsq, 2 * Lsq, eigfile);
    fclose(eigfile);

    // printf("First few values of eigvecs:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         elem = *(eigvecs + RTC(i, j, 2 * Lsq));
    //         printf("(%7.2e+%7.2ej) ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *gfuncfile = fopen("data/deg_sigma_gfuncsq_python_4x4.dat", "r");
    io_dread_2d(exp_gfuncsq, Lsq, Lsq, gfuncfile);
    fclose(gfuncfile);

    // printf("First few values of EXPECTED gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(exp_gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    utils_gfuncsq_sigma_matrix_deg(gfuncsq, sigma, eigvecs, L);

    // printf("First few values of ACTUAL gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // DTYPE diff = 0;
    // for (i = 0; i < Lsq * Lsq; i++)
    //     diff += fabs(*(gfuncsq + i) - *(exp_gfuncsq + i));

    // printf("Diff = %.2g\n", diff);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(exp_gfuncsq, gfuncsq, Lsq * Lsq);

    // char message[64];
    // CDTYPE actual, expected;
    // for(i = 0; i < Lsq; i++)
    // {
    //     for(j = 0; j < Lsq; j++)
    //     {
    //         actual = *(gfuncsq + RTC(i, j, Lsq));
    //         expected = *(exp_gfuncsq + RTC(i, j, Lsq));
    //         sprintf(message, "Error at (%d,%d) real part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(crealf(actual), crealf(expected), message);
    //         sprintf(message, "Error at (%d,%d) imag part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(cimagf(actual), cimagf(expected), message);
    //     }
    // }

    free(sigma);
    free(eigvecs);
    free(exp_gfuncsq);
    free(gfuncsq);
}

void test_utils_gfuncsq_sigma_coeff_nondeg()
{
    int L = 4;
    int Lsq = L * L;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *sigma = calloc(2 * 2, sizeof(CDTYPE));
    CDTYPE *eigvecs = calloc(2 * Lsq * 2 * Lsq, sizeof(CDTYPE));
    DTYPE *exp_gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));
    DTYPE *gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));

    *(sigma + RTC(0, 0, 2)) = 1.1 + 0*I;
    *(sigma + RTC(0, 1, 2)) = 0.2 + 0*I;
    *(sigma + RTC(1, 0, 2)) = -0.1 + 0*I;
    *(sigma + RTC(1, 1, 2)) = 0.6 + 0*I;

    // printf("Sigma:\n");
    // for (i = 0; i < 2; i++)
    // {
    //     for (j = 0; j < 2; j++)
    //     {
    //         elem = *(sigma + RTC(i, j, 2));
    //         printf("%.2g+%.2gj  ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *eigfile = fopen("data/fake_eigenvectors_4x4.dat", "r");
    io_zread_2d(eigvecs, 2 * Lsq, 2 * Lsq, eigfile);
    fclose(eigfile);

    // printf("First few values of eigvecs:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         elem = *(eigvecs + RTC(i, j, 2 * Lsq));
    //         printf("(%7.2e+%7.2ej) ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *gfuncfile = fopen("data/nondeg_sigma_coeff_gfuncsq_python_4x4.dat", "r");
    io_dread_2d(exp_gfuncsq, Lsq, Lsq, gfuncfile);
    fclose(gfuncfile);

    // printf("First few values of EXPECTED gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(exp_gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    utils_gfuncsq_sigma_coeff_nondeg(gfuncsq, sigma, eigvecs, L);

    // printf("First few values of ACTUAL gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // DTYPE diff = 0;
    // for (i = 0; i < Lsq * Lsq; i++)
    //     diff += fabs(*(gfuncsq + i) - *(exp_gfuncsq + i));

    // printf("Diff = %.2g\n", diff);

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(exp_gfuncsq, gfuncsq, Lsq * Lsq);

    free(sigma);
    free(eigvecs);
    free(exp_gfuncsq);
    free(gfuncsq);
}

void test_utils_gfuncsq_sigma_coeff_deg()
{
    int L = 4;
    int Lsq = L * L;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *sigma = calloc(2 * 2, sizeof(CDTYPE));
    CDTYPE *eigvecs = calloc(2 * Lsq * 2 * Lsq, sizeof(CDTYPE));
    DTYPE *exp_gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));
    DTYPE *gfuncsq = calloc(Lsq * Lsq, sizeof(DTYPE));

    *(sigma + RTC(0, 0, 2)) = 1.1 + 0*I;
    *(sigma + RTC(0, 1, 2)) = 0.2 + 0*I;
    *(sigma + RTC(1, 0, 2)) = -0.1 + 0*I;
    *(sigma + RTC(1, 1, 2)) = 0.6 + 0*I;

    // printf("Sigma:\n");
    // for (i = 0; i < 2; i++)
    // {
    //     for (j = 0; j < 2; j++)
    //     {
    //         elem = *(sigma + RTC(i, j, 2));
    //         printf("%.2g+%.2gj  ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *eigfile = fopen("data/fake_eigenvectors_4x4.dat", "r");
    io_zread_2d(eigvecs, 2 * Lsq, 2 * Lsq, eigfile);
    fclose(eigfile);

    // printf("First few values of eigvecs:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         elem = *(eigvecs + RTC(i, j, 2 * Lsq));
    //         printf("(%7.2e+%7.2ej) ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    FILE *gfuncfile = fopen("data/deg_sigma_coeff_gfuncsq_python_4x4.dat", "r");
    io_dread_2d(exp_gfuncsq, Lsq, Lsq, gfuncfile);
    fclose(gfuncfile);

    // printf("First few values of EXPECTED gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(exp_gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    utils_gfuncsq_sigma_coeff_deg(gfuncsq, sigma, eigvecs, L);

    // printf("First few values of ACTUAL gfuncsq:\n");
    // for (i = 0; i < 4; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3g ", *(gfuncsq + RTC(i, j, Lsq)));
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    // DTYPE diff = 0;
    // for (i = 0; i < Lsq * Lsq; i++)
    //     diff += fabs(*(gfuncsq + i) - *(exp_gfuncsq + i));

    // printf("Diff = %.2g\n", diff);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY(exp_gfuncsq, gfuncsq, Lsq * Lsq);

    free(sigma);
    free(eigvecs);
    free(exp_gfuncsq);
    free(gfuncsq);
}


void test_utils_multiply_restricted(void)
{
    CDTYPE A[25] = {(-3+5*I), (-5-4*I), (-7+1*I), (7-9*I), (9-10*I),
                    (-7-2*I), (3-7*I), (7+5*I), (7-5*I), (5+7*I),
                    (7+3*I), (-3-7*I), (-5-8*I), (-6-1*I), (5+5*I),
                    (4-8*I), (-3-4*I), (-10+9*I), (-2-2*I), (-6-8*I),
                    (7-7*I), (-9+6*I), (5+2*I), (-5-5*I), (-3+9*I)};
    CDTYPE B[25] = {(9+2*I), (-4+1*I), (-4+0*I), (-7+2*I), (-10+5*I),
                    (4+3*I), (-6+6*I), (4+9*I), (2+8*I), (2-7*I),
                    (-2-9*I), (-6-2*I), (9-5*I), (1+5*I), (-10-8*I),
                    (-5+6*I), (-8+6*I), (-6+4*I), (-8-1*I), (-3-5*I),
                    (4+1*I), (-9-3*I), (-4-10*I), (-9-4*I), (5-2*I)};
    CDTYPE exp_C[25] = {(-9-93*I), (-43+40*I), (-54+116*I), (14+35*I), (-27-81*I),
                        (-50+18*I), (21+101*I), (55-17*I), (66+56*I), (68-16*I),
                        (34+26*I), (42+40*I), (-52-68*I), (60+18*I), (94-78*I),
                        (-63+71*I), (74+77*I), (102-62*I), (75+12*I), (41+13*I),
                        -6*I, (37-20*I), (-13+24*I), (17-41*I), (-54-42*I)};

    CDTYPE * actual_C = calloc(5*5, sizeof(CDTYPE));

    utils_multiply_restricted(A, 5, 2, 4, B, 5, 3, 5, actual_C);

    int i;
    CDTYPE actual, expected;
    for(i = 0; i < 5*5; i++)
    {
        actual = *(actual_C + i);
        expected = *(exp_C + i);
        TEST_ASSERT_EQUAL_FLOAT(crealf(expected), crealf(actual));
        TEST_ASSERT_EQUAL_FLOAT(cimagf(expected), cimagf(actual));
    }

    free(actual_C);
}

void test_utils_bin_energy_range(void)
{
    DTYPE energies[100];
    int num_bins = 12;
    int actual_bin_edges[12];
    int exp_bin_edges[12] = {0, 20, 26, 32, 38, 44, 50, 56, 62, 68, 74, 80};
    int i;
    
    for(i = 0; i < 100; i++)
        energies[i] = 0.01 * i;
        // energies[i] = 0.01 * (i / 2) + 0.001;

    utils_bin_energy_range(energies, 100, num_bins, 0.20, 0.80, actual_bin_edges);
    // utils_bin_energy_range(energies, 100, num_bins, 0.10, 0.40, actual_bin_edges);
    // printf("\n");
    // for(i = 0; i < 12; i++)
    //     printf("%d ", actual_bin_edges[i]);
    // printf("\n");
    TEST_ASSERT_INT_ARRAY_WITHIN(1, exp_bin_edges, actual_bin_edges, 12);
}

void test_gfuncsq_GR_GRstar_nondeg(void)
{
    int L = 4;
    int Lsq = L * L;
    uint alpha, beta;
    int i;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *eigvecs = calloc(2*Lsq * 2*Lsq, sizeof(CDTYPE));

    // Get eigenvectors
    io_read_array('C', 'C', eigvecs, 2*Lsq, 2*Lsq,
                "data/fake_eigenvectors_4x4.dat");
    char message[32];
    for(alpha = 0; alpha < 2; alpha++)
    {
        for(beta = 0; beta < 2; beta++)
        {
            // printf("α = %d β = %d\n", alpha, beta);
            // Get expected value from file
            char filename[50];
            snprintf(filename, 50,
                    "data/nondeg_gfunc_gfuncstar_a%d_b%d_python_4x4.dat",
                    alpha, beta);

            if(alpha == beta)
            {
                DTYPE *exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
                DTYPE *gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));

                io_read_array('R', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                // Zero out the gfuncsq for the run.
                for(i = 0; i < Lsq*2*Lsq; i++)
                    *(gfuncsq + i) = 0;
                gfuncsq_sym_GR_GRstar_nondeg(eigvecs, gfuncsq, L,
                                        0, 2*Lsq, alpha);
                snprintf(message, 32, "alpha = %d beta = %d", alpha, beta);
                TEST_ASSERT_FLOAT_ARRAY_WITHIN_MESSAGE(TOL, exp_gfuncsq, gfuncsq,
                                                        Lsq*2*Lsq, message);
                // TEST_ASSERT_FLOAT_WITHIN(TOL, 1.2323432e-3, 1.2323442e-3);

                free(exp_gfuncsq);
                free(gfuncsq);

            }
            else
            {
                CDTYPE * exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                CDTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                io_read_array('C', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                gfuncsq_asym_GR_GRstar_nondeg(eigvecs, gfuncsq, L,
                                        0, 2*Lsq, alpha, beta);
                // snprintf(message, 40, "alpha = %d beta = %d gamma = %d",
                //         alpha, beta, gamma);
                assert_equal_complex_float_array(exp_gfuncsq, gfuncsq, Lsq*2*Lsq);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
        }
    }
    free(eigvecs);
}

void test_gfuncsq_GR_GRstar_deg(void)
{
    int L = 4;
    int Lsq = L * L;
    uint alpha, beta;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *eigvecs = calloc(2*Lsq * 2*Lsq, sizeof(CDTYPE));

    // Get eigenvectors
    io_read_array('C', 'C', eigvecs, 2*Lsq, 2*Lsq,
                "data/fake_eigenvectors_4x4.dat");

    char message[40];
    for(alpha = 0; alpha < 2; alpha++)
    {
        for(beta = 0; beta < 2; beta++)
        {
            // Get expected value from file
            char filename[60];
            snprintf(filename, 60,
                    "data/deg_gfunc_gfuncstar_a%d_b%d_python_4x4.dat",
                    alpha, beta);
            if(alpha == beta)
            {
                DTYPE * exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
                DTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
                io_read_array('R', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                gfuncsq_sym_GR_GRstar_deg(eigvecs, gfuncsq, L,
                                        0, 2*Lsq, alpha);
                // printf("alpha = %d beta = %d\n", alpha, beta);
                // printf("Expected:\n");
                // for(i = 0; i < 4; i++)
                // {
                //     for(j = 0; j < 4; j++)
                //     {
                //         printf("%12.3e ", *(exp_gfuncsq + RTC(i, j, Lsq)));
                //     }
                //     printf("\n");
                // }
                // printf("\n");
                // printf("Actual:\n");
                // for(i = 0; i < 4; i++)
                // {
                //     for(j = 0; j < 4; j++)
                //     {
                //         printf("%12.3e ", *(gfuncsq + RTC(i, j, Lsq)));
                //     }
                //     printf("\n");
                // }
                // printf("\n");

                snprintf(message, 40, "alpha = %d beta = %d",
                        alpha, beta);
                TEST_ASSERT_FLOAT_ARRAY_WITHIN_MESSAGE(TOL, exp_gfuncsq, gfuncsq,
                                                        Lsq*2*Lsq, message);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
            else
            {
                CDTYPE * exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                CDTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                io_read_array('C', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                gfuncsq_asym_GR_GRstar_deg(eigvecs, gfuncsq, L,
                                        0, 2*Lsq, alpha, beta);
                // snprintf(message, 40, "alpha = %d beta = %d gamma = %d",
                //         alpha, beta, gamma);
                // printf("alpha = %d beta = %d\n", alpha, beta);
                // printf("Expected:\n");
                // for(i = 0; i < 4; i++)
                // {
                //     for(j = 0; j < 4; j++)
                //     {
                //         elem = *(exp_gfuncsq + RTC(i, j, Lsq));
                //         printf("(%10.3e%+10.3ej) ",crealf(elem), cimagf(elem));
                //     }
                //     printf("\n");
                // }
                // printf("\n");
                // printf("Actual:\n");
                // for(i = 0; i < 4; i++)
                // {
                //     for(j = 0; j < 4; j++)
                //     {
                //         elem = *(gfuncsq + RTC(i, j, Lsq));
                //         printf("(%10.3e%+10.3ej) ",crealf(elem), cimagf(elem));
                //     }
                //     printf("\n");
                // }

                assert_equal_complex_float_array(exp_gfuncsq, gfuncsq, Lsq*2*Lsq);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
        }
    }

    free(eigvecs);
}

void test_gfuncsq_GR_GRstar_nondeg_restr(void)
{
    int L = 4;
    int nmin = 8;
    int nmax = 16;
    int Lsq = L * L;
    uint alpha, beta;
    int i;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *eigvecs = calloc(2*Lsq * 2*Lsq, sizeof(CDTYPE));
    
    // Get eigenvectors
    io_read_array('C', 'C', eigvecs, 2*Lsq, 2*Lsq,
                "data/fake_eigenvectors_4x4.dat");

    char message[32];
    for(alpha = 0; alpha < 2; alpha++)
    {
        for(beta = 0; beta < 2; beta++)
        {
            // Get expected value from file
            char filename[70];
            snprintf(filename, 70,
                    "data/nondeg_gfunc_gfuncstar_restr_a%d_b%d_min%d_max%d_python_4x4.dat",
                    alpha, beta, nmin, nmax);
            if(alpha == beta)
            {
                DTYPE *exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
                DTYPE *gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));

                io_read_array('R', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                // Zero out the gfuncsq for the run.
                for(i = 0; i < Lsq*2*Lsq; i++)
                    *(gfuncsq + i) = 0;
                gfuncsq_sym_GR_GRstar_nondeg(eigvecs, gfuncsq, L,
                                        nmin, nmax, alpha);
                snprintf(message, 32, "alpha = %d beta = %d", alpha, beta);
                TEST_ASSERT_FLOAT_ARRAY_WITHIN_MESSAGE(TOL, exp_gfuncsq, gfuncsq,
                                                        Lsq*2*Lsq, message);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
            else
            {
                CDTYPE * exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                CDTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                io_read_array('C', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                gfuncsq_asym_GR_GRstar_nondeg(eigvecs, gfuncsq, L,
                                        nmin, nmax, alpha, beta);
                // snprintf(message, 40, "alpha = %d beta = %d gamma = %d",
                //         alpha, beta, gamma);
                assert_equal_complex_float_array(exp_gfuncsq, gfuncsq, Lsq*2*Lsq);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
        }
    }
    free(eigvecs);
}


void test_gfuncsq_GR_GRstar_deg_restr(void)
{
    int L = 4;
    int nmin = 8;
    int nmax = 16;
    int Lsq = L * L;
    uint alpha, beta;
    // int i, j;
    // CDTYPE elem;

    CDTYPE *eigvecs = calloc(2*Lsq * 2*Lsq, sizeof(CDTYPE));

    // Get eigenvectors
    io_read_array('C', 'C', eigvecs, 2*Lsq, 2*Lsq,
                "data/fake_eigenvectors_4x4.dat");

    char message[40];
    for(alpha = 0; alpha < 2; alpha++)
    {
        for(beta = 0; beta < 2; beta++)
        {
            // Get expected value from file
            char filename[70];
            snprintf(filename, 70,
                    "data/deg_gfunc_gfuncstar_restr_a%d_b%d_min%d_max%d_python_4x4.dat",
                    alpha, beta, nmin, nmax);
            if(alpha == beta)
            {
                DTYPE * exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
                DTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(DTYPE));
                io_read_array('R', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                gfuncsq_sym_GR_GRstar_deg(eigvecs, gfuncsq, L,
                                        nmin, nmax, alpha);
                snprintf(message, 40, "alpha = %d beta = %d",
                        alpha, beta);
                TEST_ASSERT_FLOAT_ARRAY_WITHIN_MESSAGE(TOL, exp_gfuncsq, gfuncsq,
                                                        Lsq*2*Lsq, message);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
            else
            {
                CDTYPE * exp_gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                CDTYPE * gfuncsq = calloc(Lsq * 2*Lsq, sizeof(CDTYPE));
                io_read_array('C', 'C', exp_gfuncsq, Lsq, 2*Lsq, filename);
                gfuncsq_asym_GR_GRstar_deg(eigvecs, gfuncsq, L,
                                        nmin, nmax, alpha, beta);
                assert_equal_complex_float_array(exp_gfuncsq, gfuncsq, Lsq*2*Lsq);
                free(exp_gfuncsq);
                free(gfuncsq);
            }
        }
    }
    free(eigvecs);
}


int main(void)
{
    // For NaN detection
    // feenableexcept(FE_INVALID | FE_OVERFLOW);
    UNITY_BEGIN();
    // RUN_TEST(test_utils_gfuncsq_sigma_matrix_nondeg);
    // RUN_TEST(test_utils_gfuncsq_sigma_matrix_deg);
    // RUN_TEST(test_utils_gfuncsq_sigma_coeff_nondeg);
    // RUN_TEST(test_utils_gfuncsq_sigma_coeff_deg);
    RUN_TEST(test_utils_multiply_restricted);
    RUN_TEST(test_utils_bin_energy_range);
    RUN_TEST(test_gfuncsq_GR_GRstar_nondeg);
    RUN_TEST(test_gfuncsq_GR_GRstar_deg);
    RUN_TEST(test_gfuncsq_GR_GRstar_nondeg_restr);
    RUN_TEST(test_gfuncsq_GR_GRstar_deg_restr);
    return (UNITY_END());
}