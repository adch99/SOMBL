#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../extern/unity/unity.h"
#include "../src/utils/utils.h"
#include "../src/io/io.h"

void setUp(void) {}

void tearDown(void) {}

void test_utils_gfuncsq_sigma_matrix_nondeg()
{
    int L = 4;
    int Lsq = L * L;
    int i, j;
    CDTYPE elem;

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
    //         printf("%.2g+%.2gj  ", creal(elem), cimag(elem));
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
    //         printf("(%7.2e+%7.2ej) ", creal(elem), cimag(elem));
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

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(exp_gfuncsq, gfuncsq, Lsq * Lsq);

    // char message[64];
    // CDTYPE actual, expected;
    // for(i = 0; i < Lsq; i++)
    // {
    //     for(j = 0; j < Lsq; j++)
    //     {
    //         actual = *(gfuncsq + RTC(i, j, Lsq));
    //         expected = *(exp_gfuncsq + RTC(i, j, Lsq));
    //         sprintf(message, "Error at (%d,%d) real part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(creal(actual), creal(expected), message);
    //         sprintf(message, "Error at (%d,%d) imag part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(cimag(actual), cimag(expected), message);
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
    int i, j;
    CDTYPE elem;

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
    //         printf("%.2g+%.2gj  ", creal(elem), cimag(elem));
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
    //         printf("(%7.2e+%7.2ej) ", creal(elem), cimag(elem));
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

    TEST_ASSERT_EQUAL_DOUBLE_ARRAY(exp_gfuncsq, gfuncsq, Lsq * Lsq);

    // char message[64];
    // CDTYPE actual, expected;
    // for(i = 0; i < Lsq; i++)
    // {
    //     for(j = 0; j < Lsq; j++)
    //     {
    //         actual = *(gfuncsq + RTC(i, j, Lsq));
    //         expected = *(exp_gfuncsq + RTC(i, j, Lsq));
    //         sprintf(message, "Error at (%d,%d) real part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(creal(actual), creal(expected), message);
    //         sprintf(message, "Error at (%d,%d) imag part", i, j);
    //         TEST_ASSERT_EQUAL_DOUBLE_MESSAGE(cimag(actual), cimag(expected), message);
    //     }
    // }

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
        TEST_ASSERT_EQUAL_DOUBLE(creal(expected), creal(actual));
        TEST_ASSERT_EQUAL_DOUBLE(cimag(expected), cimag(actual));
    }

    free(actual_C);
}

int main(void)
{
    UNITY_BEGIN();
    RUN_TEST(test_utils_gfuncsq_sigma_matrix_nondeg);
    RUN_TEST(test_utils_gfuncsq_sigma_matrix_deg);
    RUN_TEST(test_utils_multiply_restricted);
    return (UNITY_END());
}