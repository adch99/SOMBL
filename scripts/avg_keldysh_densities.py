import numpy as np

def load(L, fname, *args, **kwargs):
    return np.loadtxt(fname, *args, **kwargs).reshape((L,L)).T

def save(L, fname, X, *args, **kwargs):
    return np.savetxt(fname, X.T.reshape(L*L), *args, **kwargs)


dirname = "densities_100x100/"
# dirname = "data"
L = 100
maxruns = 100
N = 100
disorders = np.arange(8, 18+1)
couplings = np.arange(0, 2.01, 0.1)
pattern = "altn_random_updown"

for w in disorders:
    for c in couplings:
        basename = f"mbl_density_{L}x{L}_W{w}_C{c:.1g}_TU1_TD1_N{25}_"
        for alpha in range(2):
            for beta in range(2):
                if alpha == beta:
                    avg = np.zeros((L,L))
                    var = np.zeros((L,L))
                    bname = dirname + "/" + basename + f"a{alpha}_b{beta}_full_" + pattern
                    for n in range(1, maxruns+1):
                        fname = bname + f"_n{n}.dat"
                        avg += load(L, fname)
                        var += load(L, fname + ".variance")
                    avg /= maxruns
                    var /= maxruns
                    save(L, bname + ".dat" + ".average", avg.T.reshape(100))
                    save(L, bname + ".dat" + ".variance" + ".average", var)
                if alpha != beta:
                    avg = np.zeros((L, L), dtype=complex)
                    var = np.zeros((L, L), dtype=complex)
                    bname = dirname + "/" + basename + f"a{alpha}_b{beta}_full_" + pattern
                    for n in range(1, maxruns+1):
                        fname = bname + f"_n{n}.dat"
                        avg += load(L, fname, dtype=complex)
                        var += load(L, fname + ".variance", dtype=complex)
                    avg /= maxruns
                    var /= maxruns
                    save(L, bname + ".dat" + ".average", avg)
                    save(L, bname + ".dat" + ".variance" + ".average", var)
