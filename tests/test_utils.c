#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/utils/utils.h"
#define TOL 1e-6

int test_loc_len();

int main(int argc, char ** argv)
{
    int result = test_loc_len();
    if (result)
        printf("Test Passed:\ttest_loc_len\n");
    else
        printf("Test Failed:\ttest_loc_len\n");
    return(0);
}

int test_loc_len()
{
    int len = 30;
    double * eigvals;
    double hop_strength = 1.45;
    double result;
    int i, success;

    /* Independent energy */
    double energy = 14.68;

    eigvals = malloc(sizeof(double)*len);

    for(i=0;i<len;i++)
    {
        *(eigvals + i) = i;
    }
    
    result = utils_loc_len(energy, eigvals, hop_strength, len, -1);

    if(abs(1/result - 2.096965023993401) > TOL)
    {
        success = 0;
        printf("Failed test for utils_loc_len independent energy.");
    }
    else
        success = 1;

    /* For specific eigenfunction */
    energy = *(eigvals + 20);
    result = utils_loc_len(energy, eigvals, hop_strength, len, 20);

    if(abs(1/result - 2.2728547268061026) > TOL)
    {
        success = 0;
        printf("Failed test for utils_loc_len particular eigenfunc.");
    }
    else
        success *= 1;


    free(eigvals);

    return(success);
}