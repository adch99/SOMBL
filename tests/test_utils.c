#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/utils/utils.h"
#include "../src/constants.h"
#define TOL 1e-6

// Helper Functions
int tester(int (*test_func)(), char * name);

// Test Functions
int test_loc_len();
int test_get_eigvalsh();
int test_uniform_dist();
int test_get_eigh();
int test_fit_exponential();
int test_get_green_func_lim();

int main(int argc, char ** argv)
{
    (void) argc;
    (void) argv;
    tester(test_loc_len, "test_loc_len");
    tester(test_get_eigvalsh, "test_get_eigvalsh");
    tester(test_uniform_dist, "test_uniform_dist");
    tester(test_get_eigh, "test_get_eigh");
    tester(test_fit_exponential, "test_fit_exponential");
    tester(test_get_green_func_lim, "test_get_green_func_lim");
    return(0);
}

int tester(int (*test_func)(), char * name)
{
    int result = test_func();
    if (result)
        printf("Test Passed:\t%s\n", name);
    else
        printf("Test Failed:\t%s\n", name);

    return(0);
}

int test_loc_len()
{
    int len = 30;
    DTYPE * eigvals;
    DTYPE hop_strength = 1.45;
    DTYPE result;
    int i, success;

    /* Independent energy */
    DTYPE energy = 14.68;

    eigvals = malloc(sizeof(DTYPE)*len);

    for(i=0;i<len;i++)
    {
        *(eigvals + i) = i;
    }
    
    result = utils_loc_len(energy, eigvals, hop_strength, len, -1);

    // if(fabs(1/result - 2.096965023993401) > TOL)
    if(fabs(1/result - 1.3538379111284349) > TOL)
    {
        success = 0;
        printf("Failed test for utils_loc_len independent energy.\n");
        printf("Diff: %g\n", 1/result - 1.3538379111284349);
    }
    else
        success = 1;

    /* For specific eigenfunction */
    energy = *(eigvals + 20);
    result = utils_loc_len(energy, eigvals, hop_strength, len, 20);

    // if(fabs(1/result - 2.2728547268061026) > TOL)
    if(fabs(1/result - 1.5297276139411364) > TOL)
    {
        success = 0;
        printf("Failed test for utils_loc_len particular eigenfunc.\n");
        printf("Diff: %g\n", 1/result - 1.5297276139411364);
    }
    else
        success *= 1;


    free(eigvals);

    return(success);
}

int test_get_eigvalsh()
{
    CDTYPE matrix[4][4] = {{0, 5, 0, 0}, {5, 0, 0, 0}, {0, 0, 0, -3*I}, {0, 0, 3*I, 0}};
    CDTYPE colMajorMat[16];
    DTYPE eigvals[4];    
    utils_preprocess_lapack(matrix[0], 4, colMajorMat);
    int info = utils_get_eigvalsh(colMajorMat, 4, eigvals);
    if (info != 0)
    {
        printf("Error occured: Code %d", info);
        return(0);
    }
    else
    {
        // printf("eigvals: ");
        // int i;
        // for(i = 0; i < 4; i++)
        // {
        //     printf("%lf ", *(eigvals + i));
        // }
        // printf("\n");
        DTYPE answer[4] = {-5, -3, 3, 5};
        int i, success = 1;
        for(i = 0; i < 4; i++)
        {
            if(answer[i] != *(eigvals + i))
                success = 0;
        }
        return(success);
    }

}

int test_uniform_dist()
{
    int success = 1;
    DTYPE * array = calloc(100, sizeof(DTYPE));
    utils_uniform_dist(40, 65, 100, array, 1);
    int i;
    // printf("Uniform Dist: ");
    for(i = 0; i < 100; i++)
        if(*(array + i) <= 40 || *(array + i) >= 65)
            success = 0;
        // printf("%lf ", *(array + i));
    
    // printf("\n");

    return(success);
}

int test_get_eigh()
{
    int i, j, index, success = 1;
    const CDTYPE expected_eigvecs[4][4] = {
        {-0.70710678, 0,             0,            0.70710678},
        { 0.70710678, 0,             0,            0.70710678},
        { 0,          0.70710678*I, -0.70710678*I,   0},
        { 0,          0.70710678,    0.70710678,   0}};
    const DTYPE expected_eigvals[4] = {-5, -3, 3, 5};


    CDTYPE matrix[4][4] = { {0, 5, 0,   0  },
                            {5, 0, 0,   0  },
                            {0, 0, 0,  -3*I},
                            {0, 0, 3*I, 0  }};
    CDTYPE colMajorMat[16];
    DTYPE eigvals[4];    
    utils_preprocess_lapack(matrix[0], 4, colMajorMat);
    int info = utils_get_eigh(colMajorMat, 4, eigvals);
    if (info != 0)
    {
        printf("Error occured: Code %d", info);
        return(0);
    }
    else
    { 
        for(i = 0; i < 4; i++)
        {
            if(expected_eigvals[i] != *(eigvals + i))
                success = 0;
        }

        // Check Eigenvectors
        DTYPE diff;
        for(i = 0; i < 4; i++)
        {
            for(j = 0; j < 4; j++)
            {
                // check the jth value of the ith eigenvector
                index = RTC(j, i, 4);
                diff = cabs(*(colMajorMat + index) - expected_eigvecs[j][i]);
                if(diff > TOL)
                    success = 0;
            }
        }

        if (success == 0)
        {
            // print the values
            printf("eigvals: ");
            int i;
            for(i = 0; i < 4; i++)
            {
                printf("%lf ", *(eigvals + i));
            }
            printf("\n");

            CDTYPE elem;
            printf("Expected:\n");
            for(i = 0; i < 4; i++)
            {
                for(j = 0; j < 4; j++)
                {
                    elem = expected_eigvecs[i][j];
                    printf("%08.4lf+%08.4lfj\t", creal(elem), cimag(elem));
                }
                printf("\n");
            }
            printf("\nOutput:\n");
            utils_print_matrix_complex_F(colMajorMat, 4, 4);
        }

        return(success);
    }
}

int test_get_green_func_lim()
{
    int success = 1;
    // CDTYPE matrix[4][4] = {{0.5896361835771357, 0.44745597581725566, -0.5189882797525668, 0.7517877523715167},
    //                         {0.8623924559345224, -1.180358328933656, -1.6474429460441498, -0.3026617290848911},
    //                         {-0.2920270049953856, -0.36570382127221435, 0.9146900109670734, -0.39089851247025076},
    //                         {0.20841774473573638, -0.9977012020932581, -1.718062870761319, -1.7469738133322454}};
    
    // eigenvecs is already in columnwise form.
    // ith eigenvector is in the ith column.
    CDTYPE eigenvecs[4][4] = {{0.12841189352385896, 0.4175634422815194, -0.8477995373902318, 0.30064447720278203},
                            {-0.500110400748712, -0.7515171157585533, -0.4269000372902363, 0.05355343076315152},
                            {-0.3586900050619972, 0.2979874471254812, -0.21211115881446221, -0.8588095353049181},
                            {-0.7776512411956948, 0.41480872429059784, 0.23238149424790583, 0.4113284702647177}};
    DTYPE exp_gfunc[4][4] = {{0.5554646436568518, 0.23384788281967755, 0.11660739068618453, 0.09408008283728597},
                            {0.23384788281967755, 0.4147503502608533, 0.09264378557053818, 0.25875798134893274},
                            {0.11660739068618453, 0.09264378557053818, 0.5704476737777641, 0.22030114996551375},
                            {0.09408008283728597, 0.25875798134893274, 0.22030114996551375, 0.4268607858482684}};
    DTYPE * gfunc = calloc(16, sizeof(DTYPE));
    utils_get_green_func_lim(eigenvecs[0], 4, gfunc);

    int i, j, index;
    DTYPE diff;
    for(i = 0; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            index = RTC(i, j, 4);
            diff = fabs(exp_gfunc[i][j] - *(gfunc + index));
            if(diff > TOL && isnormal(diff))
                success = 0;
        }
    }

    if(success == 0)
    {
        printf("Output:\n");
        utils_print_matrix(gfunc, 4, 4, 'R', 'F');
        printf("Expected:\n");
        utils_print_matrix(exp_gfunc, 4, 4, 'R', 'C');
    }

    return(success);
}

int test_fit_exponential()
{
    int success = 1;
    DTYPE exponent = -1.23;
    DTYPE mantissa = 113.56;
    DTYPE x[10], y[10];
    DTYPE noise[10];
    DTYPE est_exp, est_mant, residuals;
    int i;
    utils_uniform_dist(-1e-4, 1e-4, 10, noise, 1);
    for(i = 0; i < 10; i++)
    {
        x[i] = (DTYPE) i;
        y[i] = mantissa * exp(x[i] * exponent) + noise[i];
        // printf("(%.3e, %.3e)\n", x[i], y[i]);
    }

    utils_fit_exponential(y, 10, &est_exp, &est_mant, &residuals);
    // printf("Estimates: %e * exp(%e x)\n", est_mant, est_exp);;
    // printf("Residuals per datapoint: %e\n", residuals/10);

    DTYPE errp_exp = fabs((est_exp - exponent) / exponent);
    DTYPE errp_mant = fabs((est_mant - mantissa) / mantissa);

    if(errp_exp > 1e-2 || errp_mant > 1e-2)
        success = 0;


    return(success);
}