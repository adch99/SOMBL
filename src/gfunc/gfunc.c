#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <cblas.h>
#include "../constants.h"
#include "../utils/utils.h"

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
            // l = k + 1;
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
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
    <i,alpha| sigma exp(-iHt)|j,beta>|^2
*/
int gfuncsq_restr_sigma_op(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length,
                                int nmin, int nmax)
{
    gfuncsq_restr_sigma_op_nondeg(gfuncsq, sigma, eigvecs,
                                length, nmin, nmax);
    gfuncsq_restr_sigma_op_deg(gfuncsq, sigma, eigvecs,
                                length, nmin, nmax);
    return(0);
}

/*
    Computes the non-degenerate terms of pauli-generalized green's
    function square matrix and adds it to gfuncsq.
    |G(i,j)|^2 = lim t to infty |Sum_{alpha,beta}
    <i,alpha| sigma exp(-iHt)|j,beta>|^2
*/
int gfuncsq_restr_sigma_op_nondeg(DTYPE * gfuncsq, CDTYPE * sigma,
                                    CDTYPE * eigvecs, int length,
                                    int nmin, int nmax)
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

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
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
    <i,alpha| sigma exp(-iHt)|j,beta>|^2
*/
int gfuncsq_sigma_op_deg_restr(DTYPE * gfuncsq, CDTYPE * sigma,
                                CDTYPE * eigvecs, int length,
                                int nmin, int nmax)
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

    
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
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

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
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

    
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans, Lsq,
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
