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
                DTYPE hop_strength_upup, DTYPE  hop_strength_dndn,
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

    // printf("\nDisorder: [");
    // for(i = 0; i < num_sites; i++)
    // {
    //     printf("%e,", *(disorder + i));
    // }
    // printf("]\n");


    //                                      updn dnup
    DTYPE spinorbit_coeffs[NEIGHS][2] = {   {-1, 1},  // x-1
                                            {1, -1},  // x+1
                                            {I,  I},  // y-1 
                                            {-I, I}}; // y+1

    for(i = 0; i < num_sites*2; i++)
    {
        for(j = 0; j < num_sites*2; j++)
            *(ham + RTC(i, j, 2*num_sites)) = 0;
    }

    for(site1 = 0; site1 < num_sites; site1++)
    {
        *(ham + RTC(2*site1, 2*site1, 2*num_sites)) = *(disorder + site1);
        *(ham + RTC(2*site1+1, 2*site1+1, 2*num_sites)) = *(disorder + site1);

        for(loc = 0; loc < NEIGHS; loc++)
        {
            site2 = *(*(neighbours + site1) + loc);

            // Check if the neighbour actually exists
            // If not, then index will be -1,
            // so check and skip updating
            if(site2 == -1)
                continue;

            // RTC is row major to column major
            index_up_up = RTC(2*site1, 2*site2, 2*num_sites); 
            index_dn_dn = RTC(2*site1+1, 2*site2+1, 2*num_sites); 
            index_up_dn = RTC(2*site1, 2*site2+1, 2*num_sites); 
            index_dn_up = RTC(2*site1+1, 2*site2, 2*num_sites); 

            *(ham + index_up_up) = -hop_strength_upup;
            *(ham + index_dn_dn) = -hop_strength_dndn;
            *(ham + index_up_dn) = coupling_const * spinorbit_coeffs[loc][0];
            *(ham + index_dn_up) = coupling_const * spinorbit_coeffs[loc][1];
        }
    }

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
    int x, y, index;
    get_neighbours_corners(neighbours, len, width);
    get_neighbours_edges(neighbours, len, width);
    for(x = 1; x < len-1; x++)
    {
        for(y = 1; y < width-1; y++)
        {
            index = LATIDX(x, y, len);
            get_neighbours_bulk(index, len, *(neighbours + index));
        }
    }
    return(0);
}

/*
    Given the neighbour lists, updates the values
    for the neighbours of the 4 corners.
    Order: x-1, x+1, y-1, y+1
    A value of -1 indicates no neighbour exists.
    Uses open boundary condition.
*/
int get_neighbours_corners(int (*neighbours)[NEIGHS], int len, int width)
{
    int x, y, index;
    int * nlist;

    // We have 4 corners: (0,0), (0,width-1),
    // (len-1,0), (len-1,width-1)
        // (0,W-1)-----(L-1,W-1)
        //        |   |
        //        |   |
        //   (0,0)-----(L-1,0)

    // (0, 0) has only x+1 and y+1 neighbours
    x = 0;
    y = 0;
    index = LATIDX(x, y, len);
    nlist = *(neighbours + index);
    *(nlist + 0) = -1;
    *(nlist + 1) = LATIDX(x+1, y, len);
    *(nlist + 2) = -1;
    *(nlist + 3) = LATIDX(x, y+1, len);

    // (0,L-1) has only x+1 and y-1
    x = 0;
    y = width-1;
    index = LATIDX(x, y, len);
    nlist = *(neighbours + index);
    *(nlist + 0) = -1;
    *(nlist + 1) = LATIDX(x+1, y, len);
    *(nlist + 2) = LATIDX(x, y-1, len);
    *(nlist + 3) = -1;

    // (L-1,0) has only x-1 and y+1
    x = len-1;
    y = 0;
    index = LATIDX(x, y, len);
    nlist = *(neighbours + index);
    *(nlist + 0) = LATIDX(x-1, y, len);
    *(nlist + 1) = -1;
    *(nlist + 2) = -1;
    *(nlist + 3) = LATIDX(x, y+1, len);

    // (L-1,L-1) has only x-1 and y-1
    x = len-1;
    y = width-1;
    index = LATIDX(x, y, len);
    nlist = *(neighbours + index);
    *(nlist + 0) = LATIDX(x-1, y, len);
    *(nlist + 1) = -1;
    *(nlist + 2) = LATIDX(x, y-1, len);
    *(nlist + 3) = -1;

    return(0);
}
/*
    Given the neighbour lists, updates the values
    for the neighbours of the 4 edges. Does NOT
    include the 4 corners.
    Order: x-1, x+1, y-1, y+1
    A value of -1 indicates no neighbour exists.
    Uses open boundary condition.
*/
int get_neighbours_edges(int (*neighbours)[NEIGHS], int len, int width)
{
    int x, y, index;
    int * nlist;
    // Bottom edge (x, 0) has no y-1 neighbour
    y = 0;
    for(x = 1; x < len-1; x++)
    {
        index = LATIDX(x, y, len);
        nlist = *(neighbours + index);
        *(nlist + 0) = LATIDX(x-1, y, len);
        *(nlist + 1) = LATIDX(x+1, y, len);
        *(nlist + 2) = -1;
        *(nlist + 3) = LATIDX(x, y+1, len);
    }

    // Top edge (x, W-1) has no y+1 neighbour
    y = width-1;
    for(x = 1; x < len-1; x++)
    {
        index = LATIDX(x, y, len);
        nlist = *(neighbours + index);
        *(nlist + 0) = LATIDX(x-1, y, len);
        *(nlist + 1) = LATIDX(x+1, y, len);
        *(nlist + 2) = LATIDX(x, y-1, len);
        *(nlist + 3) = -1;
    }

    // Left edge (0, y) has no x-1 neighbour
    x = 0;
    for(y = 1; y < width-1; y++)
    {
        index = LATIDX(x, y, len);
        nlist = *(neighbours + index);
        *(nlist + 0) = -1;
        *(nlist + 1) = LATIDX(x+1, y, len);
        *(nlist + 2) = LATIDX(x, y-1, len);
        *(nlist + 3) = LATIDX(x, y+1, len);
    }

    // Right edge (L-1, y) has no x+1 neighbour
    x = len-1;
    for(y = 1; y < width-1; y++)
    {
        index = LATIDX(x, y, len);
        nlist = *(neighbours + index);
        *(nlist + 0) = LATIDX(x-1, y, len);
        *(nlist + 1) = -1;
        *(nlist + 2) = LATIDX(x, y-1, len);
        *(nlist + 3) = LATIDX(x, y+1, len);
    }
    return(0);
}

/*
    Given a lattice index (= i + len*j) and an array,
    writes the lattice neighbours of the site in the
    array. Order of neighbours: [x-1, x+1, y-1, y+1]
*/
int get_neighbours_bulk(int index, int len, int * nlist)
{
    /*
        Assumes 2D square lattice.
        Assumes index = x + len*y

            |---length-----|
        -   .    .    .    .
        |   .    .    .    . 
        |   .    .    .    .
        |   .    .(x,y)    .
      width .    .    .    .
        |   .    .    .    .
        |   .    .    .    .
        |   .    .    .    .
        -   .    .    .    .
        (x, y) = (1, 3)
        index = 1 + 3*4 = 13
        x = index % len
        y = index // len
    */

    int i, j;
    // int inew, jnew;

    i = index % len;
    j = index / len;

    // x-1 => (i, j-1)
    // inew = i;
    // jnew = j-1;
    *(nlist) = LATIDX(i-1, j, len);

    // x+1 => (i, j+1)
    // inew = i;
    // jnew = j+1;
    *(nlist + 1) = LATIDX(i+1, j, len);

    // y-1 => (i-1, j)
    // inew = i-1;
    // jnew = j;
    *(nlist + 2) = LATIDX(i, j-1, len);

    // y+1 => (i+1, j)
    // inew = i+1;
    // jnew = j;
    *(nlist + 3) = LATIDX(i, j+1, len);

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
    int site1, site2, index, loc;

    for(site1 = 0; site1 < num_sites; site1++)
    {
        // RTC is row to column major
        // First we set the diagonal element
        index = RTC(site1, site1, num_sites);
        *(ham + index) = *(disorder + site1);

        // And now the hopping terms
        for(loc = 0; loc < NEIGHS; loc++)
        {
            site2 = *(*(neighbours + site1) + loc);
            if(site2 == -1) // This neighbour doesn't exist.
                continue; // So skip it.
            index = RTC(site1, site2, num_sites);
            *(ham + index) = -hop_strength;
        }
    }

    free(disorder);
    return(0);
}
