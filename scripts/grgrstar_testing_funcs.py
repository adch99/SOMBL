import numpy as np

def gfunc_gfuncstar_nondeg_restr(i, alpha, beta, k, gamma, eigvecs, L, nmin, nmax):
    elem = 0
    Lsq = L*L
    for n in range(nmin, nmax):
        elem += eigvecs[2*i+alpha, n] * np.conj(eigvecs[2*i+beta, n]) * np.abs(eigvecs[2*k+gamma,n])**2    
    return elem

def gfunc_gfuncstar_deg(i, alpha, beta, k, gamma, eigvecs, L):
    elem = 0
    Lsq = L*L   
    for p in range(Lsq):
        term1 = eigvecs[2*i+alpha,2*p] * np.conj(eigvecs[2*i+beta,2*p+1])
        term2 = np.conj(eigvecs[2*k+gamma, 2*p+1]) * eigvecs[2*k+gamma,2*p]
        term3 = eigvecs[2*i+alpha,2*p+1] * np.conj(eigvecs[2*i+beta,2*p])
        term4 = np.conj(eigvecs[2*k+gamma, 2*p]) * eigvecs[2*k+gamma,2*p+1]
        elem += term1*term2 + term3*term4
    return elem

def gfunc_gfuncstar_deg_restr(i, alpha, beta, k, gamma, eigvecs, L, nmin, nmax):
    elem = 0
    Lsq = L*L
    pmin = nmin//2
    pmax = nmax//2
    for p in range(pmin, pmax):
        term1 = eigvecs[2*i+alpha,2*p] * np.conj(eigvecs[2*i+beta,2*p+1])
        term2 = np.conj(eigvecs[2*k+gamma, 2*p+1]) * eigvecs[2*k+gamma,2*p]
        term3 = eigvecs[2*i+alpha,2*p+1] * np.conj(eigvecs[2*i+beta,2*p])
        term4 = np.conj(eigvecs[2*k+gamma, 2*p]) * eigvecs[2*k+gamma,2*p+1]
        elem += term1*term2 + term3*term4
    return elem

def gfunc_gfuncstar_mat_nondeg(alpha, beta, eigvecs, L):
    Lsq = L*L
    # Create Z
    Z = np.abs(eigvecs)**2
    # Create A
    A = np.empty((Lsq,2*Lsq), dtype=complex)
    for i in range(Lsq):
        for n in range(2*Lsq):
            A[i,n] = eigvecs[2*i+alpha,n] * np.conj(eigvecs[2*i+beta, n])
    
    return A @ Z.T

def gfunc_gfuncstar_mat_deg(alpha, beta, gamma, eigvecs, L):
    Lsq = L*L
    C_ab = np.empty((Lsq,Lsq), dtype=complex)
    C_ba = np.empty((Lsq,Lsq), dtype=complex)
    C_gg = np.empty((Lsq,Lsq), dtype=complex)
    for i in range(Lsq):
        for p in range(Lsq):
            C_ab[i,p] = eigvecs[2*i+alpha,2*p] * np.conj(eigvecs[2*i+beta,2*p+1])
    for i in range(Lsq):
        for p in range(Lsq):
            C_ba[i,p] = eigvecs[2*i+beta,2*p] * np.conj(eigvecs[2*i+alpha,2*p+1])
    for i in range(Lsq):
        for p in range(Lsq):
            C_gg[i,p] = eigvecs[2*i+gamma,2*p] * np.conj(eigvecs[2*i+gamma,2*p+1])
    output = (C_ab @ C_gg.T) + np.conj(C_ba @ C_gg.T)
    gfuncsq = np.zeros((Lsq, 2*Lsq), dtype=complex)
    for i in range(Lsq):
        for k in range(Lsq):
            gfuncsq[i,2*k+gamma] = output[i,k]
    return gfuncsq



def check_deg_gfuncsq(L, eigvecs, alpha, beta):
    direct_deg_gfuncsq = np.zeros((L*L,2*L*L), dtype=complex)
    mat_deg_gfuncsq = np.zeros((L*L,2*L*L), dtype=complex)
    for gamma in range(2):
        for i in range(L*L):
            for j in range(L*L):
                direct_deg_gfuncsq[i, 2*j+gamma] = gfunc_gfuncstar_deg(i, alpha, beta, j,
                                                                       gamma, eigvecs, L)

        addend = gfunc_gfuncstar_mat_deg(alpha, beta, gamma, eigvecs, L)
        mat_deg_gfuncsq += addend

    diff = direct_deg_gfuncsq - mat_deg_gfuncsq
    filename = f"../data/deg_gfunc_gfuncstar_a{alpha}_b{beta}_python_{L}x{L}.dat"
    if alpha == beta:
        np.savetxt(filename, mat_deg_gfuncsq.real)
    else:
        np.savetxt(filename, mat_deg_gfuncsq)
    print(f"α = {alpha} β = {beta} γ = {gamma}")
    print(f"Diff Sum = {np.sum(np.abs(diff)):.3e}")
    print(f"Diff Max = {np.max(np.abs(diff)):.3e}")
    return direct_deg_gfuncsq, mat_deg_gfuncsq

def check_nondeg_gfuncsq(L, eigvecs, alpha, beta):
    direct_nondeg_gfuncsq = np.empty((L*L,2*L*L), dtype=complex)

    for alpha in range(2):
        for beta in range(2):
            for i in range(L*L):
                for j in range(L*L):
                    for gamma in range(2):
                        direct_nondeg_gfuncsq[i, 2*j+gamma] = gfunc_gfuncstar_nondeg(i, alpha, beta, j,
                                                                                    gamma, eigvecs, L)

            mat_nondeg_gfuncsq = gfunc_gfuncstar_mat_nondeg(alpha, beta, eigvecs, L)

            diff = direct_nondeg_gfuncsq - mat_nondeg_gfuncsq
            if alpha == beta:
                np.savetxt(f"../data/nondeg_gfunc_gfuncstar_a{alpha}_b{beta}_python_{L}x{L}.dat",
                        mat_nondeg_gfuncsq.real)
            else:
                np.savetxt(f"../data/nondeg_gfunc_gfuncstar_a{alpha}_b{beta}_python_{L}x{L}.dat",
                        mat_nondeg_gfuncsq)
                
            print(f"alpha: {alpha} beta: {beta}")
            print("Shape:", mat_nondeg_gfuncsq.shape)
            print(np.sum(np.abs(diff)))
            print(np.max(np.abs(diff)))
            print()

def write_gfuncsq_deg_restr(L, eigvecs, alpha, beta, nmin, nmax):
    direct_deg_gfuncsq = np.zeros((L*L,2*L*L), dtype=complex)
    for gamma in range(2):
        for i in range(L*L):
            for j in range(L*L):
                direct_deg_gfuncsq[i, 2*j+gamma] = gfunc_gfuncstar_deg_restr(i, alpha, beta, j,
                                                                            gamma, eigvecs, L,
                                                                             nmin, nmax)
    filename = f"../data/deg_gfunc_gfuncstar_restr_a{alpha}_b{beta}_min{nmin}_max{nmax}_python_{L}x{L}.dat"
    if alpha == beta:
        np.savetxt(filename, direct_deg_gfuncsq.real)
    else:
        np.savetxt(filename, direct_deg_gfuncsq)

def write_gfuncsq_nondeg_restr(L, eigvecs, alpha, beta, nmin, nmax):
    direct_nondeg_gfuncsq = np.zeros((L*L,2*L*L), dtype=complex)
    for gamma in range(2):
        for i in range(L*L):
            for j in range(L*L):
                direct_nondeg_gfuncsq[i, 2*j+gamma] = gfunc_gfuncstar_nondeg_restr(i, alpha, beta, j, gamma,
                                                                          eigvecs, L, nmin, nmax)
    filename = f"../data/nondeg_gfunc_gfuncstar_restr_a{alpha}_b{beta}_min{nmin}_max{nmax}_python_{L}x{L}.dat"
    if alpha == beta:
        np.savetxt(filename, direct_nondeg_gfuncsq.real)
    else:
        np.savetxt(filename, direct_nondeg_gfuncsq)

