#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "utils/utils.h"
#include "io/io.h"

enum Spin {UP=0, DOWN};

void spiral(DTYPE qx,  DTYPE qy, int L, CDTYPE * density);
void print_init_cond(CDTYPE * density, int L);

int main()
{
    int L = 4;
    char pattern[] = "data/spiral_qx0.25_qy0.25";
    DTYPE qx = 0.25;
    DTYPE qy = 0.25;

    CDTYPE * density;
    char filename[64];
    FILE * ofile;
    size_t writesize;

    density = calloc(2*L*L, sizeof(CDTYPE));
    // Calculate density here
    spiral(qx, qy, L, density);
    sprintf(filename, "%s_L%d.dat", pattern, L);
    ofile = io_safely_open_binary('w', filename);
    writesize = fwrite(density, sizeof(CDTYPE), 2*L*L, ofile);
    if(writesize != (size_t) 2*L*L)
    {
        fprintf(stderr, "ERROR: Complete array could not be written out to file.\n");
        exit(EXIT_FAILURE);
    }
    print_init_cond(density, L);
    fclose(ofile);
    free(density);    
}

void print_init_cond(CDTYPE * density, int L)
{
    int x, y, pos;
    CDTYPE elemc;
    for(y = 0; y < L; y++)
    {
        for(x = 0; x < L; x++)
        {
            pos = x + y*L;
            elemc = *(density + 2*pos + UP);
            // printf("%10.4le%+11.4lej ", crealf(elemc),
            //         cimagf(elemc));
            printf("%-10.4lg ", crealf(elemc));
        }
        printf("\n");
        for(x = 0; x < L; x++)
        {
            pos = x + y*L;
            elemc = *(density + 2*pos + DOWN);
            // printf("%10.4le%+11.4lej ", crealf(elemc),
            //         cimagf(elemc));
            printf("%-10.4lg ", crealf(elemc));
        }
        printf("\n");
        printf("\n");
    }
}

void spiral(DTYPE qx,  DTYPE qy, int L, CDTYPE * density)
{
    int x, y, pos;
    DTYPE phase;
    DTYPE shift = 0;
    CDTYPE alpha, beta;
    for(x = 0; x < L; x++)
    {
        for(y = 0; y < L; y++)
        {
            phase = qx * x + qy * y + shift; 
            pos = x + y*L;
            // |Ψ> = α |up> + β |down>
            alpha = cos(phase);
            beta = sin(phase);
            *(density + 2*pos + UP) = alpha * conj(alpha);
            *(density + 2*pos + DOWN) = beta * conj(beta);
        }
    }
}
