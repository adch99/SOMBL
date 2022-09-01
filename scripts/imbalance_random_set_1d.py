import numpy as np
from tqdm import tqdm

def hamiltonian(length, disorderStr, hopStr):
    diagonal = np.random.uniform(low=-disorderStr/2,
                                high=disorderStr/2,
                                size=length)
    offDiag = -hopStr * np.ones(length-1)
    ham = np.diag(diagonal) + np.diag(offDiag, k=-1) + np.diag(offDiag, k=1)
    return ham


def gfuncsq_long_time(ham, gfuncsq):
    eigvals, eigvecs = np.linalg.eigh(ham)
    for i in range(ham.shape[0]):
        for j in range(ham.shape[1]):
            terms = eigvecs[i, :] * np.conj(eigvecs[j, :])
            gfuncsq[i, j] += np.sum(np.abs(terms)**2)
    return gfuncsq


def imbalance(gfuncsq, initialState):
    length = gfuncsq.shape[0]
    mask = np.zeros(length, dtype=bool)
    mask[initialState] = True
    sigma = 2*mask - 1
    imbalance = np.sum(sigma[:, np.newaxis] * gfuncsq[:, mask])
    return imbalance


def get_initial_state(length, numParticles):
    rng = np.random.default_rng()
    initialState = rng.choice(a=length, size=numParticles,
                            replace=False)
    return initialState


def main():
    length = 1000
    numParticles = length // 2
    disorderStr = 15.0
    hopStr = 1.0
    disorderSamples = 100
    initialStateSamples = 100
    imbVals = np.empty(shape=(initialStateSamples,))

    for i in range(initialStateSamples):
        initialState = get_initial_state(length, numParticles)
        gfuncsq = np.zeros(shape=(length, length), dtype=complex)
        print(f"Run {i+1}")
        for j in tqdm(range(disorderSamples)):
            ham = hamiltonian(length, disorderStr, hopStr)
            gfuncsq_long_time(ham, gfuncsq)
        gfuncsq /= disorderSamples
        imbVals[i] = imbalance(gfuncsq, initialState)        
    
    imbMean = np.mean(imbVals)
    imbStd = np.std(imbVals)

    print(f"Mean: {imbMean}")
    print(f"Stddev: {imbStd}")
    print("Dist:", imbVals)

    filename = f"data/random_initial_1d_L{length}" \
            + f"_W{disorderStr}_T{hopStr}_DS{disorderSamples}" \
            + f"_IS{initialStateSamples}.dat"

    np.savetxt(filename, imbVals)

if __name__ == "__main__":
    main()
