#include <math.h>

double utils_loc_len(double energy, double * eigenvals, double hop_strength, int len, int eigenfunc_num)
{
    /*
    Calculates the localization length from the so-called "Lyapunov exponent"
    which is computed from the given spectrum. We assume the eigenvalues are
    sorted in ascending order and are real.
    We also assume the hopping strength is constant, however this can be easily
    extended. This can also be done for a complex hopping strength even though
    the phase doesn't actually matter, only the absolute values do.
    If eigenfunc_num is in [0,len-1], then we assume it means we are calculating it for
    that eigenfunction, otherwise, we assume the energy given is not an eigenvalue.
    */

    double lambda = 0;
    double eig;
    int i;
    
    if(eigenfunc_num >=0 && eigenfunc_num < len)
        energy = *(eigenvals + eigenfunc_num);

    for(i=0;i<len;i++)
    {
        if(i == eigenfunc_num)
            continue;
        eig = *(eigenvals + i);
        lambda += log(fabs(energy - eig)); /* + log(cabsd(*(hop_strengths+i))) if needed*/
    }

    lambda /= len;
    lambda += log(fabs(hop_strength));

    return(lambda);    
}