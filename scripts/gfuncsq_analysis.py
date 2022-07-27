import numpy as np


def getData():
    fname = "data/mbl_10x10_W16_C0_T1_greenfuncsq.dat"
    gfuncsq = np.loadtxt(f"data/mbl_{}")
    return data

def analyse(data):
    dists = []
    numDists = 0
    size = data.shape[0]
    for i in range(size):
        for j in range(size):
            dist = i*i + j*j
            try:
                index = dists.index(dist)
                dists.append((dist, 1))
            except ValueError:
                index = numDists
                dists[index] = (dists[index])
