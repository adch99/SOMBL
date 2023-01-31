#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <complex.h>
// #include <cblas.h>
#include "mkl.h"
#include "../constants.h"
#include "../utils/utils.h"
#include "gfunc.h"

#define TOL 1e-6

/*
    This function takes the eigenvectors, calculates the
    long time limit of the Green's function squared and ADDS
    it to given matrix. If you want only the Green's
    function for the given eigenvectors, then please ensure
    that the green_func array is initialized to zeroes. This
    behaviour helps us calculate the disorder averaged
    Green's function in-place which saves memory. 
*/
int gfuncsq_std(CDTYPE * eigenvectors, int size,
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
            *(eigvec_sq + RTC(i, j, size)) = creal(elem)*creal(elem)
                                            + cimag(elem)*cimag(elem); 
        }
    }

    cblas_dsyrk(CblasColMajor, CblasUpper, CblasNoTrans, size, size,
                1.0, eigvec_sq, size, 1.0, green_func, size);
    free(eigvec_sq);
    
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
                *(B + RTC(i,n,size)) = conj(*(eigenvectors + RTC(i,2*n,size)))
                                    * (*(eigenvectors + RTC(i,2*n+1,size)));
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
                *(green_func + index_ij) += 2*creal(*(BBh + index_ij));
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
DTYPE gfuncsq_std_elem(int i, int j, CDTYPE * eigenvectors,
                            int size, int degeneracy)
{
    DTYPE sum = 0.0;
    int k;
    int l;

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
        mod1 = creal(psi1)*creal(psi1) + cimag(psi1)*cimag(psi1);
        mod2 = creal(psi2)*creal(psi2) + cimag(psi2)*cimag(psi2);
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

// /*
//     Computes the pauli-generalized green's function square
//     matrix and adds it to gfuncsq.
//     |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
//     <i,alpha| sigma exp(-iHt)|j,beta>|^2
// */
// int gfuncsq_restr_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
//                                 CDTYPE * eigvecs, int length,
//                                 int nmin, int nmax)
// {
//     gfuncsq_restr_sigma_op_nondeg(gfuncsq, sigma, eigvecs,
//                                 length, nmin, nmax);
//     gfuncsq_restr_sigma_op_deg(gfuncsq, sigma, eigvecs,
//                                 length, nmin, nmax);
//     return(0);
// }

// /*
//     Computes the non-degenerate terms of pauli-generalized green's
//     function square matrix and adds it to gfuncsq.
//     |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
//     <i,alpha| sigma exp(-iHt)|j,beta>|^2
// */
// int gfuncsq_restr_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
//                                     CDTYPE * eigvecs, int length,
//                                     int nmin, int nmax)
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
//                     coeff += conj(*(sigma + RTC(alpha_p,gamma,2)))
//                             * *(sigma + RTC(alpha,gamma_p,2));
//                 }
//             }
            
//             // Check if coeff is zero
//             if(cabs(coeff) < TOL)
//                 continue; // Just skip this matrix

//             // Adds the matrices with appropriate coeffs
//             for(i = 0; i < Lsq; i++)
//             {
//                 for(n = 0; n < 2*Lsq; n++)
//                 {
//                     elem = conj(*(eigvecs + RTC(2*i+gamma,n,2*Lsq)))
//                             * *(eigvecs + RTC(2*i+gamma_p,n,2*Lsq));
//                     *(matrix1 + RTC(i,n,Lsq)) += coeff * elem;
//                     *(matrix2 + RTC(i,n,Lsq)) += elem;
//                 }
//             }
//         }
//     }

//     CDTYPE cblas_alpha = 1.0;
//     CDTYPE cblas_beta = 0.0;

//     cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
//                 Lsq, 2*Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
//                 Lsq, &cblas_beta, product, Lsq);

//     for(i = 0; i < Lsq*Lsq; i++)
//     {
//         *(gfuncsq + i) += creal(*(product + i));
//     }

//     free(matrix1);
//     free(matrix2);
//     free(product);
//     return(0);
// }

// /*
//     Computes the degenerate terms of pauli-generalized green's
//     function square matrix and adds it to gfuncsq.
//     |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
//     <i,alpha| sigma exp(-iHt)|j,beta>|^2
// */
// int gfuncsq_restr_sigma_op_deg(DTYPE * gfuncsq, CDTYPE * sigma,
//                                 CDTYPE * eigvecs, int length,
//                                 int nmin, int nmax)
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
//                     coeff += conj(*(sigma + RTC(alpha_p,gamma,2)))
//                             * *(sigma + RTC(alpha,gamma_p,2));
//                 }
//             }

//             // Check if coeff is zero
//             if(cabs(coeff) < TOL)
//                 continue; // Just skip this matrix

//             // Adds the matrices with appropriate coeffs
//             for(i = 0; i < Lsq; i++)
//             {
//                 for(p = 0; p < Lsq; p++)
//                 {
//                     elem = conj(*(eigvecs + RTC(2*i+gamma,2*p,2*Lsq)))
//                         * *(eigvecs + RTC(2*i+gamma_p,2*p+1,2*Lsq));
//                     *(matrix1 + RTC(i,p,Lsq)) += coeff * elem;
//                     *(matrix2 + RTC(i,p,Lsq)) += elem;
//                 }
//             }
//         }
//     }

//     CDTYPE cblas_alpha = 1.0;
//     CDTYPE cblas_beta = 0.0;

    
//     cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
//                 Lsq, Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
//                 Lsq, &cblas_beta, product, Lsq);

//     for(i = 0; i < Lsq; i++)
//     {
//         for(j = 0; j < Lsq; j++)
//         {
//             *(gfuncsq + RTC(i, j, Lsq)) += 2 * creal(*(product + RTC(i, j, Lsq)));
//         }
//     }

//     free(matrix1);
//     free(matrix2);
//     free(product);
//     return(0);
// }


/*
    Computes the pauli-generalized green's function square
    matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty
    |Sum_{alpha,beta} <i,alpha| sigma exp(-iHt)|j,beta>|^2
*/
int gfuncsq_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
                    CDTYPE * eigvecs, int length)
{
    gfuncsq_sigma_op_nondeg(gfuncsq, sigma, eigvecs, length);
    gfuncsq_sigma_op_deg(gfuncsq, sigma, eigvecs, length);
    return(0);
}

/*
    Computes the non-degenerate terms of pauli-generalized
    green's function square matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
    <i,alpha|sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
*/
int gfuncsq_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
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
                    coeff += conj(*(sigma + RTC(alpha_p,gamma,2))) 
                            * *(sigma + RTC(alpha,gamma_p,2));
                }
            }
            
            // Check if coeff is zero
            if(cabs(coeff) < TOL)
                continue; // Just skip this matrix

            // Adds the matrices with appropriate coeffs
            for(i = 0; i < Lsq; i++)
            {
                for(n = 0; n < 2*Lsq; n++)
                {
                    elem = conj(*(eigvecs + RTC(2*i+gamma,n,2*Lsq)))
                            * *(eigvecs + RTC(2*i+gamma_p,n,2*Lsq));
                    *(matrix1 + RTC(i,n,Lsq)) += coeff * elem;
                    *(matrix2 + RTC(i,n,Lsq)) += elem;
                }
            }
        }
    }

    CDTYPE cblas_alpha = 1.0;
    CDTYPE cblas_beta = 0.0;

    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
                Lsq, 2*Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
                Lsq, &cblas_beta, product, Lsq);
    
    for(i = 0; i < Lsq*Lsq; i++)
    {
        *(gfuncsq + i) += creal(*(product + i));
    }

    free(matrix1);
    free(matrix2);
    free(product);
    return(0);
}

/*
    Computes the degenerate terms of pauli-generalized green's
    function square matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
    <i,alpha| sigma_{alpha,beta} exp(-iHt)|j,beta>|^2
*/
int gfuncsq_sigma_op_deg(DTYPE * gfuncsq, CDTYPE * sigma,
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
                    coeff += conj(*(sigma + RTC(alpha_p,gamma,2)))
                            * *(sigma + RTC(alpha,gamma_p,2));
                }
            }
            
            // Check if coeff is zero
            if(cabs(coeff) < TOL)
                continue; // Just skip this matrix

            // Adds the matrices with appropriate coeffs
            for(i = 0; i < Lsq; i++)
            {
                for(p = 0; p < Lsq; p++)
                {
                    elem = conj(*(eigvecs + RTC(2*i+gamma,2*p,2*Lsq)))
                            * *(eigvecs + RTC(2*i+gamma_p,2*p+1,2*Lsq));
                    *(matrix1 + RTC(i,p,Lsq)) += coeff * elem;
                    *(matrix2 + RTC(i,p,Lsq)) += elem;
                }
            }
        }
    }

    CDTYPE cblas_alpha = 1.0;
    CDTYPE cblas_beta = 0.0;

    
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
                Lsq, Lsq, &cblas_alpha, matrix1, Lsq, matrix2,
                Lsq, &cblas_beta, product, Lsq);

    for(i = 0; i < Lsq; i++)
    {
        for(j = 0; j < Lsq; j++)
        {
            *(gfuncsq + RTC(i, j, Lsq)) += 2 * creal(*(product + RTC(i, j, Lsq)));
        }
    }

    free(matrix1);
    free(matrix2);
    free(product);
    return(0);
}

int gfuncsq_sym_GR_GRstar_nondeg(CDTYPE * eigvecs, DTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha)
{
    int num_sites = length * length;
    int num_states = 2*num_sites;
    int nwindow = nmax - nmin;
    // Case where alpha = beta
    int i; // Positional index
    int n, m; // Pos+Spin index or eigenvector index
    int q; // index that is summed over
    CDTYPE elem;

    DTYPE * Z = calloc(num_states * nwindow, sizeof(DTYPE));
    DTYPE * A = calloc(num_sites * nwindow, sizeof(DTYPE));
    // printf("A = %d x %d\n", num_sites, nwindow);

    // Create Z_{2i+alpha,n} = |psi_n(i,alpha)|^2
    for(m = 0; m < num_states; m++)
    {
        for(n = nmin; n < nmax; n++)
        {

            // printf("Index m = %d\tn = %d\n", m, n);
            q = n - nmin;
            elem = *(eigvecs + RTC(m, n, num_states));
            *(Z + RTC(m, q, num_states)) = creal(elem)*creal(elem)
                                        + cimag(elem)*cimag(elem);            
        }
    }
    // printf("Z Initialized!\n\n");

    // Create A_{alpha,alpha} =     
    for(i = 0; i < num_sites; i++)
    {
        for(n = nmin; n < nmax; n++)
        {
            // printf("(%d, %d) ", i, n);
            q = n - nmin;
            *(A + RTC(i, q, num_sites)) = *(eigvecs + RTC(2*i+alpha, n, num_states))
                                        * conj(*(eigvecs + RTC(2*i+alpha, n, num_states)));
        }
        // printf("\n");
    }

    lapack_int Lsq = num_sites;
    lapack_int twoLsq = num_states;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                Lsq, twoLsq, nwindow, 1.0, A, Lsq,
                Z, twoLsq, 1.0, gfuncsq, Lsq);
    free(Z);
    free(A);

    return(0);
}

int gfuncsq_asym_GR_GRstar_nondeg(CDTYPE * eigvecs, CDTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha, uint beta)
{
    int num_sites = length * length;
    int num_states = 2*num_sites;
    int nwindow = nmax - nmin;
    // Case where alpha = beta
    int i; // Positional index
    int n, m; // Pos+Spin index or eigenvector index
    int q; // index that is summed over
    CDTYPE elem;

    CDTYPE * Z = calloc(num_states * nwindow, sizeof(CDTYPE));
    CDTYPE * A = calloc(num_sites * nwindow, sizeof(CDTYPE));
    // printf("A = %d x %d\n", num_sites, nwindow);

    // Create Z_{2i+alpha,n} = |psi_n(i,alpha)|^2
    for(m = 0; m < num_states; m++)
    {
        for(n = nmin; n < nmax; n++)
        {

            // printf("Index m = %d\tn = %d\n", m, n);
            q = n - nmin;
            elem = *(eigvecs + RTC(m, n, num_states));
            *(Z + RTC(m, q, num_states)) = creal(elem)*creal(elem)
                                        + cimag(elem)*cimag(elem);            
        }
    }
    // printf("Z Initialized!\n\n");

    // Create A_{alpha,alpha} =     
    for(i = 0; i < num_sites; i++)
    {
        for(n = nmin; n < nmax; n++)
        {
            // printf("(%d, %d) ", i, n);
            q = n - nmin;
            *(A + RTC(i, q, num_sites)) = *(eigvecs + RTC(2*i+alpha, n, num_states))
                                        * conj(*(eigvecs + RTC(2*i+beta, n, num_states)));
        }
        // printf("\n");
    }

    lapack_int Lsq = num_sites;
    lapack_int twoLsq = num_states;
    lapack_complex_double blas_alpha, blas_beta;
    blas_alpha.real = 1.0;
    blas_alpha.imag = 0.0;
    blas_beta.real = 1.0;
    blas_beta.imag = 0.0;
    cblas_cgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                Lsq, twoLsq, nwindow, &blas_alpha, A, Lsq,
                Z, twoLsq, &blas_beta, gfuncsq, Lsq);
    free(Z);
    free(A);

    return(0);
}


int gfuncsq_asym_GR_GRstar_deg(CDTYPE * eigvecs, CDTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha, uint beta)
{
    int Lsq = length * length;
    CDTYPE * C_ab, * C_ba_star, * C_gg, * output;
    int i, k, p, q;
    int pmin = nmin / 2;
    int pmax = nmax / 2;
    uint gamma;

    lapack_int M = Lsq;
    lapack_int N = Lsq;
    lapack_int K = (nmax - nmin) / 2;
    lapack_complex_double blas_alpha, blas_beta;

    output = calloc(Lsq*Lsq, sizeof(CDTYPE));

    for(gamma = 0; gamma < 2; gamma++)
    {
        // Create C_gg
        C_gg = calloc(Lsq*Lsq, sizeof(CDTYPE));
        for(i = 0; i < Lsq; i++)
        {
            for(p = pmin; p < pmax; p++)
            {
                q = p - pmin;
                *(C_gg + RTC(i, q, Lsq)) = *(eigvecs + RTC(2*i+gamma,2*p,2*Lsq))
                                            * conj(*(eigvecs + RTC(2*i+gamma, 2*p+1, 2*Lsq)));
            }
        }
        // Create C_ab
        C_ab = calloc(Lsq*Lsq, sizeof(CDTYPE));
        for(i = 0; i < Lsq; i++)
        {
            for(p = pmin; p < pmax; p++)
            {
                q = p - pmin;
                *(C_ab + RTC(i, q, Lsq)) = *(eigvecs + RTC(2*i+alpha,2*p,2*Lsq))
                                            * conj(*(eigvecs + RTC(2*i+beta, 2*p+1, 2*Lsq)));
            }
        }

        // compute C_ab C_gg^T and add to output
        // blas_beta = 0 is important to clear out
        // any pre-existing values.
        blas_alpha.real = 1.0;
        blas_alpha.imag = 0.0;
        blas_beta.real = 0.0;
        blas_beta.imag = 0.0;
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    M, N, K, &blas_alpha, C_ab, M, C_gg, N,
                    &blas_beta, output, M);
        free(C_ab);

        // Create C_ba_star
        C_ba_star = calloc(Lsq*Lsq, sizeof(CDTYPE));
        for(i = 0; i < Lsq; i++)
        {
            for(p = pmin; p < pmax; p++)
            {
                q = p - pmin;
                *(C_ba_star + RTC(i, q, Lsq)) = conj(*(eigvecs + RTC(2*i+beta, 2*p, 2*Lsq)))
                                        * (*(eigvecs + RTC(2*i+alpha, 2*p+1, 2*Lsq)));
            }
        }

        // compute C_ba_star C_gg^H and add to output
        blas_alpha.real = 1.0;
        blas_alpha.imag = 0.0;
        blas_beta.real = 1.0;
        blas_beta.imag = 0.0;
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
                    M, N, K, &blas_alpha, C_ba_star, M, C_gg, N,
                    &blas_beta, output, M);
        free(C_ba_star);
        free(C_gg);

        // Add the output to the gfuncsq
        for(i = 0; i < Lsq; i++)
        {
            for(k = 0; k < Lsq; k++)
            {
                *(gfuncsq + RTC(i, 2*k+gamma, Lsq)) += *(output + RTC(i, k, Lsq));
            }
        }
    }
    free(output);
    return(0);
}

int gfuncsq_sym_GR_GRstar_deg(CDTYPE * eigvecs, DTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha)
{
    int Lsq = length * length;
    CDTYPE * C_ab, * C_gg, * output;
    int i, k, p, q;
    uint gamma;
    uint beta = alpha;
    int pmin = nmin / 2;
    int pmax = nmax / 2;

    lapack_int M = Lsq;
    lapack_int N = Lsq;
    lapack_int K = (nmax - nmin) / 2;
    // lapack_complex_double blas_alpha, blas_beta;
    CDTYPE blas_alpha, blas_beta;

    output = calloc(Lsq*Lsq, sizeof(CDTYPE));

    for(gamma = 0; gamma < 2; gamma++)
    {
        // for(i = 0; i < Lsq*Lsq; i++)
        // {
        //     *(output + i) = 0;
        // }

        // Create C_gg
        C_gg = calloc(Lsq*Lsq, sizeof(CDTYPE));
        for(i = 0; i < Lsq; i++)
        {
            for(p = pmin; p < pmax; p++)
            {
                q = p - pmin;
                *(C_gg + RTC(i, q, Lsq)) = *(eigvecs + RTC(2*i+gamma,2*p,2*Lsq))
                                            * conj(*(eigvecs + RTC(2*i+gamma, 2*p+1, 2*Lsq)));
            }
        }
        // Create C_ab
        C_ab = calloc(Lsq*Lsq, sizeof(CDTYPE));
        for(i = 0; i < Lsq; i++)
        {
            for(p = pmin; p < pmax; p++)
            {
                q = p - pmin;
                *(C_ab + RTC(i, q, Lsq)) = *(eigvecs + RTC(2*i+alpha,2*p,2*Lsq))
                                            * conj(*(eigvecs + RTC(2*i+beta, 2*p+1, 2*Lsq)));
            }
        }

        // compute C_ab C_gg^T and add to output
        // blas_beta = 0 is important to clear out
        // any pre-existing values.
        // blas_alpha.real = 1.0;
        // blas_alpha.imag = 0.0;
        // blas_beta.real = 0.0;
        // blas_beta.imag = 0.0;
        blas_alpha = 1.0;
        blas_beta = 0.0;
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    M, N, K, &blas_alpha, C_ab, M, C_gg, N,
                    &blas_beta, output, M);
        free(C_ab);
        free(C_gg);

        // if(gamma == 1)
        // {
        //     CDTYPE elem;
        //     printf("gamma = 1 output:\n");
        //     for(i = 0; i < 4; i++)
        //     {
        //         for(k = 0; k < 4; k++)
        //         {
        //             elem = *(output + RTC(i, k, Lsq));
        //             printf("(%10.3e%+10.3ej) ",creal(elem), cimag(elem));
        //         }
        //         printf("\n");
        //     }
        //     printf("\n");
        // }

        // Add the output to the gfuncsq
        for(i = 0; i < Lsq; i++)
        {
            for(k = 0; k < Lsq; k++)
            {
                *(gfuncsq + RTC(i, 2*k+gamma, Lsq)) += 2*creal(*(output + RTC(i, k, Lsq)));
            }
        }
    }
    free(output);
    return(0);
}



int gfunc_calc_direct(CDTYPE * eigvecs, int length, int i, int k, uint alpha,
                    uint beta, uint gamma, int nmin, int nmax)
{
    int n;
    int Lsq = length*length;
    CDTYPE sum = 0;
    CDTYPE * ptr_n_ialpha = eigvecs + RTC(2*i+alpha, 0, 2*Lsq);
    CDTYPE * ptr_n_ibeta = eigvecs + RTC(2*i+beta, 0, 2*Lsq);
    CDTYPE * ptr_n_kgamma = eigvecs + RTC(2*k+gamma, 0, 2*Lsq);
    size_t next_col = RTC(2*i+alpha, 1, 2*Lsq) - RTC(2*i+alpha, 0, 2*Lsq);

    for(n = nmin; n < nmax; n++)
    {
        sum += (*ptr_n_ialpha) * conj(*ptr_n_ibeta) * (creal(*ptr_n_kgamma)*creal(*ptr_n_kgamma)
                + cimag(*ptr_n_kgamma)*cimag(*ptr_n_kgamma));
        ptr_n_ialpha += next_col;
        ptr_n_ibeta += next_col;
        ptr_n_kgamma += next_col;
    }
    return(sum);
}

int gfunc_direct_full(CDTYPE * eigvecs, CDTYPE * gfuncsq, int length,
                        int nmin, int nmax, uint alpha, uint beta)
{
    int i, k;
    int Lsq = length * length;
    uint gamma;
    CDTYPE elem;
    for(i = 0; i < Lsq; i++)
    {
        for(k = 0; k < Lsq; k++)
        {
            for(gamma = 0; gamma < 2; gamma++)
            {
                elem = gfunc_calc_direct(eigvecs, length, i, k, alpha, beta,
                                        gamma, nmin, nmax);
                *(gfuncsq + RTC(i, 2*k+gamma, Lsq)) += elem;
            }
        }
    }
    return 0;
}   