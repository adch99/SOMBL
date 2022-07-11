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

int main(int argc, char ** argv)
{
    tester(test_loc_len, "test_loc_len");
    tester(test_get_eigvalsh, "test_get_eigvalsh");
    tester(test_uniform_dist, "test_uniform_dist");
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