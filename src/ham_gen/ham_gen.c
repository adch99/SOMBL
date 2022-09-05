/*
Functions to create Spin-Orbit coupled Hamiltonians
for a 2d spin lattice
*/
#include <stdio.h>
#include <lapacke.h>
#include <omp.h>
#include "ham_gen.h"
#include "../utils/utils.h"

int hamiltonian(CDTYPE * ham, int len, int width,
                DTYPE coupling_const, DTYPE disorder_strength,
                DTYPE hop_strength_upup, DTYPE hop_strength_dndn,
                int (*neighbours)[NEIGHS])
{
    int num_sites = len*width;
    // Produce disorder
    DTYPE * disorder;
    disorder = malloc(sizeof(DTYPE)*num_sites);
    DTYPE low = -disorder_strength / 2.0;
    DTYPE high = disorder_strength / 2.0;
    utils_uniform_dist(low, high, num_sites, disorder, 0);
    // Produce matrix
    int site1, site2;
    int index_up_up, index_dn_dn, index_up_dn, index_dn_up;
    int i, j, loc;

    for(i = 0; i < num_sites*2; i++)
    {
        for(j = 0; j < num_sites*2; j++)
            *(ham + RTC(i, j, 2*num_sites)) = 0;
    }

    for(site1 = 0; site1 < num_sites; site1++)
    {
        for(loc = 0; loc < NEIGHS; loc++)
        {
            site2 = *(*(neighbours + site1) + loc);
            // RTC is row major to column major
            index_up_up = RTC(2*site1, 2*site2, 2*num_sites); 
            index_dn_dn = RTC(2*site1+1, 2*site2+1, 2*num_sites); 
            index_up_dn = RTC(2*site1, 2*site2+1, 2*num_sites); 
            index_dn_up = RTC(2*site1+1, 2*site2, 2*num_sites); 

            switch(loc)
            {
                case 0:
                    *(ham + index_up_up) = -hop_strength_upup;
                    *(ham + index_dn_dn) = -hop_strength_dndn;
                    *(ham + index_up_dn) = -coupling_const;
                    *(ham + index_dn_up) = coupling_const;
                    break; 

                case 1:
                    *(ham + index_up_up) = -hop_strength_upup;
                    *(ham + index_dn_dn) = -hop_strength_dndn;
                    *(ham + index_up_dn) = coupling_const;
                    *(ham + index_dn_up) = -coupling_const;
                    break; 

                case 2:
                    *(ham + index_up_up) = -hop_strength_upup;
                    *(ham + index_dn_dn) = -hop_strength_dndn;
                    *(ham + index_up_dn) = I*coupling_const;
                    *(ham + index_dn_up) = I*coupling_const;
                    break; 

                case 3:
                    *(ham + index_up_up) = -hop_strength_upup;
                    *(ham + index_dn_dn) = -hop_strength_dndn;
                    *(ham + index_up_dn) = -I*coupling_const;
                    *(ham + index_dn_up) = I*coupling_const;
                    break; 
            }
        }
    }

    // // #pragma omp parallel for
    // for(site1 = 0; site1 < num_sites; site1++)
    // {
    //     // #pragma omp parallel for
    //     for(site2 = 0; site2 <= site1; site2++)
    //     {
    //         // RTC is row major to column major
    //         index_up_up = RTC(2*site1, 2*site2, 2*num_sites); 
    //         index_dn_dn = RTC(2*site1+1, 2*site2+1, 2*num_sites); 
    //         index_up_dn = RTC(2*site1, 2*site2+1, 2*num_sites); 
    //         index_dn_up = RTC(2*site1+1, 2*site2, 2*num_sites); 

    //         if (site1 == site2)
    //         {
    //             // #pragma omp critical
    //             {
    //                 *(ham + index_up_up) = *(disorder + site1);
    //                 *(ham + index_dn_dn) = *(disorder + site1);
    //                 *(ham + index_up_dn) = 0;
    //                 *(ham + index_dn_up) = 0;
    //             }
    //         }
    //         else
    //         {
    //             int loc = check_neighbour(site1, *(neighbours + site2));
                
    //             switch(loc)
    //             {
    //                 case -1:
    //                     *(ham + index_up_up) = 0;
    //                     *(ham + index_dn_dn) = 0;
    //                     *(ham + index_up_dn) = 0;
    //                     *(ham + index_dn_up) = 0;
    //                     break;

    //                 case 0:
    //                     *(ham + index_up_up) = -hop_strength_upup;
    //                     *(ham + index_dn_dn) = -hop_strength_dndn;
    //                     *(ham + index_up_dn) = coupling_const;
    //                     *(ham + index_dn_up) = -coupling_const;
    //                     break; 

    //                 case 1:
    //                     *(ham + index_up_up) = -hop_strength_upup;
    //                     *(ham + index_dn_dn) = -hop_strength_dndn;
    //                     *(ham + index_up_dn) = -coupling_const;
    //                     *(ham + index_dn_up) = coupling_const;
    //                     break; 

    //                 case 2:
    //                     *(ham + index_up_up) = -hop_strength_upup;
    //                     *(ham + index_dn_dn) = -hop_strength_dndn;
    //                     *(ham + index_up_dn) = -I*coupling_const;
    //                     *(ham + index_dn_up) = -I*coupling_const;
    //                     break; 

    //                 case 3:
    //                     *(ham + index_up_up) = -hop_strength_upup;
    //                     *(ham + index_dn_dn) = -hop_strength_dndn;
    //                     *(ham + index_up_dn) = I*coupling_const;
    //                     *(ham + index_dn_up) = -I*coupling_const;
    //                     break; 
    //             }
    //         }
    //     }
    // }

    free(disorder);
    return(0);
}

/*
    Creates a neighbour list for each index
    of a 2D lattice of len x width.
    Order of neighbours: x-1, x+1, y-1, y+1
*/
int get_neighbour_lists(int (*neighbours)[NEIGHS], int len, int width)
{


    int index;

    // #pragma omp parallel for
    for(index = 0; index < len*width; index++)
    {
        get_neighbours(index, len, width, *(neighbours + index));
    }

    return(0);
}

/*
    Given a lattice index (= i + len*j) and an array,
    writes the lattice neighbours of the site in the
    array. Order of neighbours: [x-1, x+1, y-1, y+1]
*/
int get_neighbours(int index, int len, int width, int * nlist)
{
    /*
        Assumes 2D square lattice.
        Assumes Periodic Boundary Conditions.
        Assumes index = i + len*j

            |---length-----|
        -   .    .    .    .
        |   .    .    .    . 
        |   .    .    .    .
        |   .    .(i,j)    .
      width .    .    .    .
        |   .    .    .    .
        |   .    .    .    .
        |   .    .    .    .
        -   .    .    .    .
        (i,j) = (1, 3)
        index = 1 + 3*4 = 13
        i = index % len
        j = index // len
    */

   int i, j, inew, jnew;

    i = index % len;
    j = index / len;

    // x-1 => (i, j-1)
    inew = i;
    jnew = ((j-1) >= 0)? j-1 : j-1+width;
    *(nlist) = inew + len*jnew;

    // x+1 => (i, j+1)
    inew = i;
    jnew = ((j+1) < width)? j+1 : j+1-width;
    *(nlist + 1) = inew + len*jnew;

    // y-1 => (i-1, j)
    inew = ((i-1) >= 0)? i-1 : i-1+len;
    jnew = j;
    *(nlist + 2) = inew + len*jnew;

    // y+1 => (i+1, j)
    inew = ((i+1) < len)? i+1 : i+1-len;
    jnew = j;
    *(nlist + 3) = inew + len*jnew;

    return(0);
}

// Check if index is in the nlist.
// If yes, then return the index.
// 0->x-1, 1->x+1, 2->y-1, 3->y+1
// -1->Not neighbour
int check_neighbour(int index, int * nlist)
{

    int i, loc = -1;
    for(i = 0; i < NEIGHS; i++)
    {
        if(*(nlist + i) == index)
            loc = i;
    } 
    return(loc);
}


int hamiltonian_nospin(CDTYPE * ham, int len, int width,
                DTYPE disorder_strength, DTYPE hop_strength,
                int (*neighbours)[NEIGHS])
{
    int num_sites = len*width;
    // Produce disorder
    DTYPE * disorder;
    disorder = malloc(sizeof(DTYPE)*num_sites);
    DTYPE low = -disorder_strength / 2.0;
    DTYPE high = disorder_strength / 2.0;
    utils_uniform_dist(low, high, num_sites, disorder, 0);

    // Produce matrix
    int site1, site2, index;

    for(site1 = 0; site1 < num_sites; site1++)
    {
        for(site2 = 0; site2 < num_sites; site2++)
        {
            // RTC is row to column major
            index = RTC(site1, site2, num_sites);

            if (site1 == site2)
                *(ham + index) = *(disorder + site1);

            else
            {
                int loc = check_neighbour(site1, *(neighbours + site2));

                if (loc == -1)
                    *(ham + index) = 0;
                else
                    *(ham + index) = -hop_strength;
            }
        }
    }

    free(disorder);
    return(0);
}
