#include <stdio.h>
#include <stdlib.h>
#include "utils/utils.h"
#include "constants.h"

int main(int argc, char ** argv)
{
    int length = 10;
    DTYPE hop_strength = 1;
    DTYPE disorder_strength = 10;
    CDTYPE * ham = calloc(length*length, sizeof(CDTYPE));
    ham_gen_1d(ham, length, hop_strength, disorder_strength);
    
    CDTYPE * eigvals
    

    return(0);
}

int ham_gen_1d(CDTYPE * ham, int length, DTYPE hop_strength, DTYPE disorder_strength)
{
    DTYPE * disorder;
    utils_uniform_dist(0, disorder_strength, length, disorder, 0);

    int i;
    for(i = 0; i < length; i++)
    {
        *(ham + i*length + i) = *(disorder + i);
    }

    for(i = 0; i < length-1; i++)
    {
        *(ham + i*length + (i+1)) = -hop_strength;
        *(ham + (i+1)*length + i) = -hop_strength;
    }

    return 0;
}