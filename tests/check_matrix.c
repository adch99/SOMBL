#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../src/utils/utils.h"
#include "../src/ham_gen/ham_gen.h"

int main(int argc, char ** argv)
{
    (void) argc;
    (void) argv;
    int length, width;
    length = width = 50;
    int num_sites = length*width;
    int num_states = num_sites;
    int i, j, index;
    CDTYPE value;
    CDTYPE * ham = calloc(num_states*num_states, sizeof(CDTYPE));
    int (*neighbour)[NEIGHS] = malloc((num_sites * NEIGHS) * sizeof(int));
    get_neighbour_lists(neighbour, length, width);
    hamiltonian_nospin(ham, length, width, 15.0, 1.0, neighbour);

    for(i = 0; i < num_states; i++)
    {
        for(j = 0; j < num_states; j++)
        {
            index = RTC(i, j, num_states);
            value = *(ham + index);
            printf("%e+%ej ", creal(value), cimag(value));
        }
        printf("\n");
    }
    printf("\n");

    return 0;
}