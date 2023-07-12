#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include "mkl.h"
// #include <cblas.h>
// #include <omp.h>
#include "utils.h"
#include "../constants.h"

#define TOL 1e-6
#define THRESHOLD TOL

enum Spin {UP=0, DOWN};

/*
    Calculates the localization length from the so-called
    "Lyapunov exponent" which is computed from the given
    spectrum. We assume the eigenvalues are sorted in
    ascending order and are real. We also assume the hopping
    strength is constant, however this can be easily
    extended. This can also be done for a complex hopping
    strength even though the phase doesn't actually matter,
    only the absolute values do. If eigenfunc_num is in
    [0,len-1], then we assume it means we are calculating it
    for that eigenfunction, otherwise, we assume the energy
    given is not an eigenvalue. This works only for 1D systems.
*/
DTYPE utils_loc_len(DTYPE energy, const DTYPE * eigenvals, DTYPE hop_strength,
                    int len, int eigenfunc_num)
{
    DTYPE lambda = 0;
    DTYPE eig;
    DTYPE loc_len;
    int i;
    int skip_one;
    skip_one = (eigenfunc_num >=0 && eigenfunc_num < len);
    if (skip_one)
        energy = *(eigenvals + eigenfunc_num);

    for(i=0;i<len;i++)
    {
        if(skip_one && i == eigenfunc_num)
            continue;
        eig = *(eigenvals + i);
        lambda += log(fabs(energy - eig));
        /* + log(cabsd(*(hop_strengths+i))) if needed*/
    }

    if (skip_one)
        lambda /= (DTYPE) (len - 1);
    else
        lambda /= (DTYPE) len;
    lambda += -log(fabs(hop_strength));

    loc_len = 1.0 / lambda;
    return(loc_len);    
}

int utils_save_matrix(void * matrix, int m, int n,
                    char type, char ordering, FILE * ofile)
{
    int i, j, index;
    DTYPE elemr;
    CDTYPE elemc;


    if(type != 'C' && type != 'R')
    {
        printf("Invalid type passed to utils_print_matrix.\n");
        printf("type must be either 'C' or 'R'.\n");
    }


    if(ordering != 'C' && ordering != 'F')
    {
        printf("Invalid ordering passed to utils_print_matrix.\n");
        printf("ordering must be either 'C' or 'F'.\n");
    }


    for(i = 0; i < m; i++)
    {
        for(j = 0; j < n; j++)
        {
            index = (ordering=='C')?(i*n + j):RTC(i, j, m);
            if(type == 'C')
            {
                elemc = *((CDTYPE*)matrix + index);
                fprintf(ofile, "%e+%ej ", crealf(elemc), cimagf(elemc));
            }
            else
            {
                elemr = *((DTYPE*)matrix + index);
                fprintf(ofile, "%e ", elemr);
            }
        }
        fprintf(ofile, "\n");
    }
    return 0;
}

/*
    This function takes the eigenvectors, calculates the
    long time limit of the Green's function squared and ADDS
    it to given matrix. If you want only the Green's
    function for the given eigenvectors, then please ensure
    that the green_func array is initialized to zeroes. This
    behaviour helps us calculate the disorder averaged
    Green's function in-place which saves memory. 
*/
int utils_get_green_func_lim(CDTYPE * eigenvectors, int size,
                            DTYPE * green_func, int degeneracy)
{

    // Construct green function limit squared
    int i, j, n;

    // First let's square the eigenvectors matrix
    // element by element
    DTYPE * eigvec_sq = calloc(size*size, sizeof(DTYPE));
    CDTYPE elem;


    for(i = 0; i < size; i++)
    {
        for(j = 0; j < size; j++)
        {
            elem = *(eigenvectors + RTC(i, j, size));
            *(eigvec_sq + RTC(i, j, size)) = crealf(elem)*crealf(elem) + cimagf(elem)*cimagf(elem); 
        }
    }

    cblas_ssyrk(CblasColMajor, CblasUpper, CblasNoTrans, size, size,
                1.0, eigvec_sq, size, 1.0, green_func, size);
    free(eigvec_sq);
    
    // for(i = 0; i < size; i++)
    // {
    //     for(j = 0; j <= i; j++)
    //     {
    //         DTYPE sum = utils_compute_gfsq_elem(i, j, eigenvectors, size,
    //                                             degeneracy);
    //         *(green_func + RTC(i,j,size)) += sum;
    //         if(i != j)
    //             *(green_func + RTC(j,i,size)) += sum;
    //     }
    // }

    if (degeneracy == DEGEN_EIGVALS)
    {
        CDTYPE * B = calloc(size*size/2, sizeof(CDTYPE));
        CDTYPE * BBh = calloc(size*size, sizeof(CDTYPE));

        // Create B = (N x N/2) matrix
        // B_in = psi_{2n}(i)^*  psi_{2n+1}(i)
        for(i = 0; i < size; i++)
        {
            for(n = 0; n < size/2; n++)
            {
                *(B + RTC(i,n,size)) = conj(*(eigenvectors + RTC(i,2*n,size))) * (*(eigenvectors + RTC(i,2*n+1,size)));
            }
        }
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,
                    size, size/2, 1.0, B, size, 0.0, BBh, size);
        // printf("Bello\n");
        int index_ij;
        for(i = 0; i < size; i++)
        {
            for(j = i; j < size; j++)
            {
                index_ij = RTC(i, j, size);
                *(green_func + index_ij) += 2*crealf(*(BBh + index_ij));
            }
        }

        free(B);
        free(BBh);

    }
    // Only upper half matrix is updated by BLAS.
    // Copying to lower half.
    utils_reflect_upper_to_lower(green_func, size);
    return 0;
}

/*
    Computes the (i,j)th element of the Green's function
    squared in the long time limit.
*/
DTYPE utils_compute_gfsq_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy)
{
    DTYPE sum = 0.0;
    int k, l;

    int index1, index2;
    index1 = RTC(i, 0, size);
    index2 = RTC(j, 0, size);
    // // We can eliminate some operations
    // // by adding 'size' to index1
    // // rather than calling RTC each time.
    // // This removes one multiply computation
    // // from the L^6 block.

    // // If you want to parallelize this block,
    // // move this back down into the block as
    // // index1 and index2 will need to be local
    // // for parallelization.
    // // index1 = RTC(i, k, size);
    // // index2 = RTC(j, k, size);

    CDTYPE psi1, psi2;
    DTYPE mod1, mod2;

    // // #pragma omp parallel for reduction (+:sum) schedule(auto)
    for(k = 0; k < size; k++)
    {
        psi1 = *(eigenvectors + index1);
        psi2 = *(eigenvectors + index2);
        mod1 = crealf(psi1)*crealf(psi1) + cimagf(psi1)*cimagf(psi1);
        mod2 = crealf(psi2)*crealf(psi2) + cimagf(psi2)*cimagf(psi2);
        sum += mod1 * mod2;
        index1 += size;
        index2 += size;
    }

    if(degeneracy == DEGEN_EIGVALS)
    {
        // We assume that the degeneracy occurs in
        // pairs and that we have all the eigenvalues
        // being degenerate (Kramer degeneracy).
        
        int index_ik, index_il, index_jk, index_jl;

        CDTYPE psi_k_i, psi_l_i;
        CDTYPE psi_k_j, psi_l_j;
        CDTYPE term;

        // #pragma omp parallel for reduction (+:sum) schedule(auto)
        for(k = 0; k < size; k += 2)
        {
            l = k + 1;
            index_ik = RTC(i, k, size);
            index_il = RTC(i, l, size);
            index_jk = RTC(j, k, size);
            index_jl = RTC(j, l, size);

            psi_k_i = *(eigenvectors + index_ik);
            psi_l_i = *(eigenvectors + index_il);
            psi_k_j = *(eigenvectors + index_jk);
            psi_l_j = *(eigenvectors + index_jl);

            term = psi_k_i * conj(psi_l_i) * conj(psi_k_j) * psi_l_j;
            sum += term + conj(term);
        }
    }
    return(sum);
}

/*
    Computes the pauli-generalized green's function square
    matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} <i,alpha| sigma exp(-iHt)|j,beta>|^2
*/
int utils_gfuncsq_sigma_matrix(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length)
{
    utils_gfuncsq_sigma_matrix_nondeg(gfuncsq, sigma, eigvecs, length);
    utils_gfuncsq_sigma_matrix_deg(gfuncsq, sigma, eigvecs, length);
    return(0);
}

/*
    Computes the non-degenerate terms of pauli-generalized green's
    function square matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} <i,alpha| sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
*/
int utils_gfuncsq_sigma_matrix_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                                    CDTYPE * eigvecs, int length)
{
    uint alpha, alpha_p, gamma, gamma_p;
    CDTYPE coeff;
    int Lsq = length*length; 
    int i, n;

    CDTYPE * matrix1 = calloc(Lsq*(2*Lsq), sizeof(CDTYPE));
    CDTYPE * matrix2 = calloc(Lsq*(2*Lsq), sizeof(CDTYPE));
    CDTYPE * product = calloc(Lsq*Lsq, sizeof(CDTYPE));
    CDTYPE elem;

    for(gamma = 0; gamma < 2; gamma++)
    {
        for(gamma_p = 0; gamma_p < 2; gamma_p++)
        {
            // First we create the coefficients
            coeff = 0;
            for(alpha = 0; alpha < 2; alpha++)
            {
                for(alpha_p = 0; alpha_p < 2; alpha_p++)
                {
                    coeff += conj(*(sigma + RTC(alpha_p,gamma,2))) * *(sigma + RTC(alpha,gamma_p,2));
                }
            }
            // printf("g=%d,g'=%d  coeff=%e+%e\n", gamma, gamma_p, crealf(coeff), cimagf(coeff));
            
            // Check if coeff is zero
            if(cabs(coeff) < TOL)
                continue; // Just skip this matrix

            // Adds the matrices with appropriate coeffs
            for(i = 0; i < Lsq; i++)
            {
                for(n = 0; n < 2*Lsq; n++)
                {
                    elem = conj(*(eigvecs + RTC(2*i+gamma,n,2*Lsq))) * *(eigvecs + RTC(2*i+gamma_p,n,2*Lsq));
                    *(matrix1 + RTC(i,n,Lsq)) += coeff * elem;
                    *(matrix2 + RTC(i,n,Lsq)) += elem;
                }
            }
        }
    }

    // int j;
    // printf("First few values of matrix1:\n");
    // for (i = 0; i < Lsq; i++)
    // {
    //     for (j = 0; j < 2*Lsq; j++)
    //     {
    //         elem = *(matrix1 + RTC(i,j,Lsq));
    //         printf("%.3g+%.3gj ", crealf(elem), cimagf(elem));
    //     }
    //     printf("\n");
    // }
    // printf("\n");




    CDTYPE cblas_alpha = 1.0;
    CDTYPE cblas_beta = 0.0;

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
                Lsq, 2*Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
                Lsq, &cblas_beta, product, Lsq);
    
    // printf("G(0,0) = %le+%lej\n", crealf(*(product + 1)), cimagf(*(product + 1)));

    for(i = 0; i < Lsq*Lsq; i++)
    {
        *(gfuncsq + i) += crealf(*(product + i));
    }

    free(matrix1);
    free(matrix2);
    free(product);
    return(0);
}

/*
    Computes the degenerate terms of pauli-generalized green's
    function square matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} <i,alpha| sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
*/
int utils_gfuncsq_sigma_matrix_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length)
{
    uint alpha, alpha_p, gamma, gamma_p;
    CDTYPE coeff;
    int Lsq = length*length; 
    int i, j, p;

    CDTYPE * matrix1 = calloc(Lsq*Lsq, sizeof(CDTYPE));
    CDTYPE * matrix2 = calloc(Lsq*Lsq, sizeof(CDTYPE));
    CDTYPE * product = calloc(Lsq*Lsq, sizeof(CDTYPE));
    CDTYPE elem;


    for(gamma = 0; gamma < 2; gamma++)
    {
        for(gamma_p = 0; gamma_p < 2; gamma_p++)
        {
            // First we create the coefficients
            coeff = 0;
            for(alpha = 0; alpha < 2; alpha++)
            {
                for(alpha_p = 0; alpha_p < 2; alpha_p++)
                {
                    coeff += conj(*(sigma + RTC(alpha_p,gamma,2))) * *(sigma + RTC(alpha,gamma_p,2));
                }
            }
            // printf("g=%d,g'=%d  coeff=%e+%e\n", gamma, gamma_p, crealf(coeff), cimagf(coeff));
            
            // Check if coeff is zero
            if(cabs(coeff) < TOL)
                continue; // Just skip this matrix

            // Adds the matrices with appropriate coeffs
            for(i = 0; i < Lsq; i++)
            {
                for(p = 0; p < Lsq; p++)
                {
                    elem = conj(*(eigvecs + RTC(2*i+gamma,2*p,2*Lsq))) * *(eigvecs + RTC(2*i+gamma_p,2*p+1,2*Lsq));
                    *(matrix1 + RTC(i,p,Lsq)) += coeff * elem;
                    *(matrix2 + RTC(i,p,Lsq)) += elem;
                }
            }
        }
    }

    CDTYPE cblas_alpha = 1.0;
    CDTYPE cblas_beta = 0.0;

    
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
                Lsq, Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
                Lsq, &cblas_beta, product, Lsq);

    for(i = 0; i < Lsq; i++)
    {
        for(j = 0; j < Lsq; j++)
        {
            *(gfuncsq + RTC(i, j, Lsq)) += 2 * crealf(*(product + RTC(i, j, Lsq)));
        }
    }

    free(matrix1);
    free(matrix2);
    free(product);
    return(0);
}

/*
    Computes the pauli-generalized green's function square
    matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} <i,alpha| sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
*/
// int utils_gfuncsq_sigma_matrix_restr(DTYPE * gfuncsq, CDTYPE * sigma,
//                                 CDTYPE * eigvecs, int length, int nmin, int nmax)
// {
//     utils_gfuncsq_sigma_matrix_nondeg_restr(gfuncsq, sigma, eigvecs, length, nmin, nmax);
//     utils_gfuncsq_sigma_matrix_deg_restr(gfuncsq, sigma, eigvecs, length, nmin, nmax);
//     return(0);
// }

/*
    Computes the non-degenerate terms of pauli-generalized green's
    function square matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} <i,alpha| sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
*/
// int utils_gfuncsq_sigma_matrix_nondeg_restr(DTYPE * gfuncsq, CDTYPE * sigma,
//                                     CDTYPE * eigvecs, int length, int nmin, int nmax)
// {
//     uint alpha, alpha_p, gamma, gamma_p;
//     CDTYPE coeff;
//     int Lsq = length*length; 
//     int i, n;

//     CDTYPE * matrix1 = calloc(Lsq*(2*Lsq), sizeof(CDTYPE));
//     CDTYPE * matrix2 = calloc(Lsq*(2*Lsq), sizeof(CDTYPE));
//     CDTYPE * product = calloc(Lsq*Lsq, sizeof(CDTYPE));
//     CDTYPE elem;

//     for(gamma = 0; gamma < 2; gamma++)
//     {
//         for(gamma_p = 0; gamma_p < 2; gamma_p++)
//         {
//             // First we create the coefficients
//             coeff = 0;
//             for(alpha = 0; alpha < 2; alpha++)
//             {
//                 for(alpha_p = 0; alpha_p < 2; alpha_p++)
//                 {
//                     coeff += conj(*(sigma + RTC(alpha_p,gamma,2))) * *(sigma + RTC(alpha,gamma_p,2));
//                 }
//             }
//             // printf("g=%d,g'=%d  coeff=%e+%e\n", gamma, gamma_p, crealf(coeff), cimagf(coeff));
            
//             // Check if coeff is zero
//             if(cabs(coeff) < TOL)
//                 continue; // Just skip this matrix

//             // Adds the matrices with appropriate coeffs
//             for(i = 0; i < Lsq; i++)
//             {
//                 for(n = 0; n < 2*Lsq; n++)
//                 {
//                     elem = conj(*(eigvecs + RTC(2*i+gamma,n,2*Lsq))) * *(eigvecs + RTC(2*i+gamma_p,n,2*Lsq));
//                     *(matrix1 + RTC(i,n,Lsq)) += coeff * elem;
//                     *(matrix2 + RTC(i,n,Lsq)) += elem;
//                 }
//             }
//         }
//     }

//     CDTYPE cblas_alpha = 1.0;
//     CDTYPE cblas_beta = 0.0;

//     cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
//                 Lsq, 2*Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
//                 Lsq, &cblas_beta, product, Lsq);

//     for(i = 0; i < Lsq*Lsq; i++)
//     {
//         *(gfuncsq + i) += crealf(*(product + i));
//     }

//     free(matrix1);
//     free(matrix2);
//     free(product);
//     return(0);
// }

// /*
//     Computes the degenerate terms of pauli-generalized green's
//     function square matrix and adds it to gfuncsq.
//     |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} <i,alpha| sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
// */
// int utils_gfuncsq_sigma_matrix_deg_restr(DTYPE * gfuncsq, CDTYPE * sigma,
//                                 CDTYPE * eigvecs, int length, int nmin, int nmax)
// {
//     uint alpha, alpha_p, gamma, gamma_p;
//     CDTYPE coeff;
//     int Lsq = length*length; 
//     int i, j, p;

//     CDTYPE * matrix1 = calloc(Lsq*Lsq, sizeof(CDTYPE));
//     CDTYPE * matrix2 = calloc(Lsq*Lsq, sizeof(CDTYPE));
//     CDTYPE * product = calloc(Lsq*Lsq, sizeof(CDTYPE));
//     CDTYPE elem;


//     for(gamma = 0; gamma < 2; gamma++)
//     {
//         for(gamma_p = 0; gamma_p < 2; gamma_p++)
//         {
//             // First we create the coefficients
//             coeff = 0;
//             for(alpha = 0; alpha < 2; alpha++)
//             {
//                 for(alpha_p = 0; alpha_p < 2; alpha_p++)
//                 {
//                     coeff += conj(*(sigma + RTC(alpha_p,gamma,2))) * *(sigma + RTC(alpha,gamma_p,2));
//                 }
//             }
//             // printf("g=%d,g'=%d  coeff=%e+%e\n", gamma, gamma_p, crealf(coeff), cimagf(coeff));
            
//             // Check if coeff is zero
//             if(cabs(coeff) < TOL)
//                 continue; // Just skip this matrix

//             // Adds the matrices with appropriate coeffs
//             for(i = 0; i < Lsq; i++)
//             {
//                 for(p = 0; p < Lsq; p++)
//                 {
//                     elem = conj(*(eigvecs + RTC(2*i+gamma,2*p,2*Lsq))) * *(eigvecs + RTC(2*i+gamma_p,2*p+1,2*Lsq));
//                     *(matrix1 + RTC(i,p,Lsq)) += coeff * elem;
//                     *(matrix2 + RTC(i,p,Lsq)) += elem;
//                 }
//             }
//         }
//     }

//     CDTYPE cblas_alpha = 1.0;
//     CDTYPE cblas_beta = 0.0;

    
//     cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
//                 Lsq, Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
//                 Lsq, &cblas_beta, product, Lsq);

//     for(i = 0; i < Lsq; i++)
//     {
//         for(j = 0; j < Lsq; j++)
//         {
//             *(gfuncsq + RTC(i, j, Lsq)) += 2 * crealf(*(product + RTC(i, j, Lsq)));
//         }
//     }

//     free(matrix1);
//     free(matrix2);
//     free(product);
//     return(0);
// }



DTYPE utils_gfuncsq_sigma(int i, uint alpha, int j, uint beta,
                        CDTYPE * sigma, CDTYPE * eigvecs, int num_states)
{
    DTYPE elem = 0;
    int n;
    uint gamma, gamma_p;
    
    CDTYPE sigma_a_g, sigma_a_gp, psi_ig_n, psi_igp_n, psi_jb_n;
    CDTYPE term1, term2, term3;

    for(n = 0; n < num_states; n++)
    {
        for(gamma = 0; gamma < 2; gamma++)
        {
            for(gamma_p = 0; gamma_p < 2; gamma_p++)
            {
                sigma_a_g = *(sigma + RTC(alpha, gamma, 2));
                sigma_a_gp = *(sigma + RTC(alpha, gamma_p, 2));
                psi_ig_n = *(eigvecs + RTC(2*i+gamma, n, num_states));
                psi_igp_n = *(eigvecs + RTC(2*i+gamma_p, n, num_states));
                psi_jb_n = *(eigvecs + RTC(2*j+beta, n, num_states));

                term1 = conj(sigma_a_g) * sigma_a_gp;
                term2 = conj(psi_ig_n) * psi_igp_n;
                term3 = conj(psi_jb_n) * psi_jb_n;
                elem += term1 * term2 * term3;
            }
        }
    }
    return(elem);
}

DTYPE utils_pbc_chord_length_sq(int index1, int length1, int index2, int length2)
{
    DTYPE radius1 = (DTYPE) length1 / (2.0 * M_PI);
    DTYPE radius2 = (DTYPE) length2 / (2.0 * M_PI);

    DTYPE angle1 = 2.0 * M_PI * ((DTYPE) index1 / (DTYPE) length1);
    DTYPE angle2 = 2.0 * M_PI * ((DTYPE) index2 / (DTYPE) length2);

    DTYPE chord1 = 2 * radius1 * sin(angle1 / 2.0);
    DTYPE chord2 = 2 * radius2 * sin(angle2 / 2.0);

    DTYPE chordsq = chord1*chord1 + chord2*chord2;
    return chordsq;   
}

/*
    Given matrix mat1 of size m1 x n1
    and matrix mat2 of size m2 x n2
    Computes A @ B^H where
    A = mat1[:, n1min:n1max]
    B = mat2[:, n2min:n2max]
*/
int utils_multiply_restricted(CDTYPE * mat1, int m1, int n1min, int n1max,
                            CDTYPE * mat2, int m2, int n2min, int n2max,
                            CDTYPE * prod)
{
    if((n1max-n1min) != (n2max-n2min))
    {
        printf("Error! Mismatch in sizes of restricted ranges! %d != %d\n", (n1max-n1min), 
                (n2max-n2min));
        return(-1);
    }

    int nlen = n1max - n1min;
    CDTYPE * A = mat1 + n1min*m1;
    CDTYPE * B = mat2 + n2min*m2;
    CDTYPE alpha = 1.0;
    CDTYPE beta = 1.0;

    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, 
                m1, m2, nlen, &alpha,
                A, m1, B, m2, &beta, prod, m1);

    return(0);
}

int utils_gfuncsq_sigma_coeff_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length)
{
    int Lsq = length * length;
    CDTYPE * matrix1 = calloc(Lsq*2*Lsq, sizeof(CDTYPE));
    CDTYPE * matrix2 = calloc(Lsq*2*Lsq, sizeof(CDTYPE));
    CDTYPE * product = calloc(Lsq*Lsq, sizeof(CDTYPE));
    // CDTYPE * inter = calloc(Lsq*Lsq, sizeof(CDTYPE));

    int index1, index2;
    uint alpha, alpha_p, beta, beta_p;
    int i, j, n;
    CDTYPE coeff = 0;
    CDTYPE elem, b;

    for(index1 = 0; index1 < 4; index1++)
    {
        alpha = index1 / 2;
        alpha_p = index1 % 2;

        // Build matrix1
        for(i = 0; i < Lsq; i++)
        {
            for(n = 0; n < 2*Lsq; n++)
            {
                *(matrix1 + RTC(i,n,Lsq)) = conj(*(eigvecs + RTC(2*i+alpha,n,2*Lsq))) 
                                            * (*(eigvecs + RTC(2*i+alpha_p,n,2*Lsq)));
            }
        }

        // Calculate self product: M_aa' @ M_aa'^H
        beta = alpha;
        beta_p = alpha_p;
        coeff = conj(*(sigma + RTC(alpha,beta,2)))
                * (*(sigma + RTC(alpha_p,beta_p,2)));
        if(cabs(coeff) > THRESHOLD)
        {
            // printf("%dx%d a=%d a'=%d b=%d b'=%d\n",
            //     index1, index1, alpha, alpha_p, beta, beta_p);
            // printf("index1 = %d\n", index1);
            cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,
                        Lsq, 2*Lsq, coeff, matrix1, Lsq,
                        0.0, matrix2, Lsq);
            // We are using matrix2 as just a placeholder for the
            // result because the output is going to be only filled
            // in the upper triangular part.
            // Even though matrix2 is intended to be a L^2 x 2L^2
            // matrix, it is ultimately just a 1d array the way it
            // is allocated.

            for(i = 0; i < Lsq; i++)
            {
                *(product + RTC(i,i,Lsq)) += *(matrix2 + RTC(i,i,Lsq));
                for(j = i+1; j < Lsq; j++)
                {
                    elem = *(matrix2 + RTC(i,j,Lsq));
                    *(product + RTC(i,j,Lsq)) += elem;
                    *(product + RTC(j,i,Lsq)) += elem;
                }
            }
        }
        for(index2 = index1+1; index2 < 4; index2++)
        {
            beta = index2 / 2;
            beta_p = index2 % 2;

            // printf("%dx%d a=%d a'=%d b=%d b'=%d\n",
            //         index1, index2, alpha, alpha_p, beta, beta_p);
            // printf("%dx%d a=%d a'=%d b=%d b'=%d\n",
            //         index2, index1, beta, beta_p, alpha, alpha_p);
            // printf("index1 = %d index2 = %d\n", index1, index2);

            // Build matrix2
            for(i = 0; i < Lsq; i++)
            {
                for(n = 0; n < 2*Lsq; n++)
                {
                    *(matrix2 + RTC(i,n,Lsq)) = conj(*(eigvecs + RTC(2*i+beta,n,2*Lsq)))
                                                * (*(eigvecs + RTC(2*i+beta_p,n,2*Lsq)));
                }
            }

            // Compute matrix1 @ matrix2
            // M_aa' @ M_bb'^H
            coeff = conj(*(sigma + RTC(alpha, beta, 2)))
                    * (*(sigma + RTC(alpha_p, beta_p, 2)));
            b = 1.0;
            if(cabs(coeff) > THRESHOLD)
            {
                // Do the computation
                cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                            Lsq, Lsq, 2*Lsq, &coeff, matrix1, Lsq,
                            matrix2, Lsq, &b, product, Lsq);
            }

            // Compute matrix2 @ matrix1
            // M_bb'^H @ M_aa'
            coeff = conj(*(sigma + RTC(beta, alpha, 2)))
                    * (*(sigma + RTC(beta_p, alpha_p, 2)));
            b = 1.0;
            if(cabs(coeff) > THRESHOLD)
            {
                // Do the computation
                cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                            Lsq, Lsq, 2*Lsq, &coeff, matrix2, Lsq,
                            matrix1, Lsq, &b, product, Lsq);
            }
        }
    }

    for(i = 0; i < Lsq; i++)
    {
        for(j = 0; j < Lsq; j++)
        {
            elem = *(product + RTC(i,j,Lsq));
            *(gfuncsq + RTC(i,j,Lsq)) += crealf(elem);
            // printf("(%d,%d) = %.2g + %.2gj\n", i, j, crealf(elem), cimagf(elem));
            // if(fabs(cimagf(elem)) > TOL)
            //     printf("Problem at (%d,%d), imag part is %.3g", i, j, cimagf(elem));
        }
    }


    return(0);
}

int utils_gfuncsq_sigma_coeff_deg(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length)
{
    int Lsq = length * length;
    CDTYPE * matrix1 = calloc(Lsq*Lsq, sizeof(CDTYPE));
    CDTYPE * matrix2 = calloc(Lsq*Lsq, sizeof(CDTYPE));
    CDTYPE * product = calloc(Lsq*Lsq, sizeof(CDTYPE));
    // CDTYPE * inter = calloc(Lsq*Lsq, sizeof(CDTYPE));

    int index1, index2;
    uint alpha, alpha_p, beta, beta_p;
    int i, j, p;
    CDTYPE coeff = 0;
    CDTYPE elem, b;

    for(index1 = 0; index1 < 4; index1++)
    {
        alpha = index1 / 2;
        alpha_p = index1 % 2;

        // Build matrix1
        for(i = 0; i < Lsq; i++)
        {
            for(p = 0; p < Lsq; p++)
            {
                *(matrix1 + RTC(i,p,Lsq)) = conj(*(eigvecs + RTC(2*i+alpha,2*p,2*Lsq))) 
                                            * (*(eigvecs + RTC(2*i+alpha_p,2*p+1,2*Lsq)));
            }
        }

        // Calculate self product: M_aa' @ M_aa'^H
        beta = alpha;
        beta_p = alpha_p;
        coeff = conj(*(sigma + RTC(alpha,beta,2)))
                * (*(sigma + RTC(alpha_p,beta_p,2)));
        if(cabs(coeff) > THRESHOLD)
        {
            // printf("%dx%d a=%d a'=%d b=%d b'=%d\n",
            //     index1, index1, alpha, alpha_p, beta, beta_p);
            // printf("index1 = %d\n", index1);
            cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,
                        Lsq, Lsq, coeff, matrix1, Lsq,
                        0.0, matrix2, Lsq);
            // We are using matrix2 as just a placeholder for the
            // result because the output is going to be only filled
            // in the upper triangular part.
            // Even though matrix2 is intended to be a L^2 x 2L^2
            // matrix, it is ultimately just a 1d array the way it
            // is allocated.

            for(i = 0; i < Lsq; i++)
            {
                *(product + RTC(i,i,Lsq)) += *(matrix2 + RTC(i,i,Lsq));
                for(j = i+1; j < Lsq; j++)
                {
                    elem = *(matrix2 + RTC(i,j,Lsq));
                    *(product + RTC(i,j,Lsq)) += elem;
                    *(product + RTC(j,i,Lsq)) += elem;
                }
            }
        }
        for(index2 = index1+1; index2 < 4; index2++)
        {
            beta = index2 / 2;
            beta_p = index2 % 2;

            // printf("%dx%d a=%d a'=%d b=%d b'=%d\n",
            //         index1, index2, alpha, alpha_p, beta, beta_p);
            // printf("%dx%d a=%d a'=%d b=%d b'=%d\n",
            //         index2, index1, beta, beta_p, alpha, alpha_p);
            // printf("index1 = %d index2 = %d\n", index1, index2);

            // Build matrix2
            for(i = 0; i < Lsq; i++)
            {
                for(p = 0; p < Lsq; p++)
                {
                    *(matrix2 + RTC(i,p,Lsq)) = conj(*(eigvecs + RTC(2*i+beta,2*p,2*Lsq)))
                                                * (*(eigvecs + RTC(2*i+beta_p,2*p+1,2*Lsq)));
                }
            }

            // Compute matrix1 @ matrix2
            // M_aa' @ M_bb'^H
            coeff = conj(*(sigma + RTC(alpha, beta, 2)))
                    * (*(sigma + RTC(alpha_p, beta_p, 2)));
            b = 1.0;
            if(cabs(coeff) > THRESHOLD)
            {
                // Do the computation
                cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                            Lsq, Lsq, Lsq, &coeff, matrix1, Lsq,
                            matrix2, Lsq, &b, product, Lsq);
            }

            // Compute matrix2 @ matrix1
            // M_bb'^H @ M_aa'
            coeff = conj(*(sigma + RTC(beta, alpha, 2)))
                    * (*(sigma + RTC(beta_p, alpha_p, 2)));
            b = 1.0;
            if(cabs(coeff) > THRESHOLD)
            {
                // Do the computation
                cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                            Lsq, Lsq, Lsq, &coeff, matrix2, Lsq,
                            matrix1, Lsq, &b, product, Lsq);
            }
        }
    }

    for(i = 0; i < Lsq; i++)
    {
        for(j = 0; j < Lsq; j++)
        {
            elem = *(product + RTC(i,j,Lsq));
            *(gfuncsq + RTC(i,j,Lsq)) += 2*crealf(elem);
            // printf("(%d,%d) = %.2g + %.2gj\n", i, j, crealf(elem), cimagf(elem));
            // if(fabs(cimagf(elem)) > TOL)
            //     printf("Problem at (%d,%d), imag part is %.3g", i, j, cimagf(elem));
        }
    }


    return(0);
}

/*
    Computes the pauli-generalized green's function square
    matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta} sigma_{alpha,beta} <i,alpha| exp(-iHt)|j,beta>|^2
*/
int utils_gfuncsq_sigma_coeff(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length)
{
    utils_gfuncsq_sigma_coeff_nondeg(gfuncsq, sigma, eigvecs, length);
    utils_gfuncsq_sigma_coeff_deg(gfuncsq, sigma, eigvecs, length);
    return(0);
}
