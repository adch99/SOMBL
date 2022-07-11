import numpy as np
import matplotlib.pyplot as plt

"""
    This is a test for determining if our formula
    and procedure is indeed correct or not.
    For this purpose, we will be using a spinless
    1d Anderson model with uniform disorder. 
"""


def ham_gen(length, hop_strength, disorder_strength):
    ham = np.zeros((length, length))
    disorder = np.random.uniform(0, disorder_strength, size=length)

    for i in range(length-1):
        ham[i, i+1] = -hop_strength
        ham[i+1, i] = -hop_strength

    return ham

def loc_len(spectrum, hop_strength, energy=None, index=None):
    length = spectrum.shape[0]
    mask = np.full(length, True)

    if index is None and energy is None:
        raise ValueError("At least one of index and energy have to be given.")
    elif index is not None and energy is not None:
        raise ValueError("Both index and energy cannot be given.")
    elif index is not None:
        energy = spectrum[index]
        mask[index] = False
        denom = length - 1
    else:
        denom = length

    return np.sum(np.log(np.abs(spectrum[mask] - energy))) / denom - np.log(np.abs(hop_strength))

def run():
    hop_strength = 1
    disorder_strength = 15
    length = 1000
    ham = ham_gen(length, hop_strength, disorder_strength)
    spectrum = np.linalg.eigvalsh(ham)

    locLens = np.empty(length)
    for i in range(length):
        locLens[i] = loc_len(spectrum, hop_strength, index=i)
    
    avgLocLen = np.mean(locLens)

    print(f"Avg Loc Len: {avgLocLen:e}")

    plt.plot(np.arange(length), locLens, marker=".")
    plt.xlabel("Eigenvalue Index")
    plt.ylabel(r"$\xi_{loc}$")
    
    fig = plt.figure()
    plt.plot(spectrum, locLens, marker=".")
    plt.xlabel("Energy")
    plt.ylabel(r"$\xi_{loc}$")

if __name__ == "__main__":
    run()
    plt.show()