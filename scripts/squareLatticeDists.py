import bisect as bi
import numpy as np
import matplotlib.pyplot as plt
import scripts.plot_utils as putils
from tqdm import tqdm

class DistSq(object):
    def __init__(self, value, count):
        self.value = value
        self.count = count

def get_series_sum(phi, latAvals, latAcounts, latBvals, latBcounts):
    seriesL = np.sum(np.exp(-np.sqrt(latAvals) / phi) * latAcounts)
    seriesM = np.sum(np.exp(-np.sqrt(latBvals) / phi) * latBcounts)

    return seriesL, seriesM

def get_series(Lx, Ly):
    latA = []
    latB = []

    keyfunc = lambda x: x.value
    for x in range(-Lx, Lx+1):
        for y in range(-Ly, Ly+1):
            value = x*x + y*y
            if (x + y) % 2 == 0:
                # check if the value is already in the list
                index = bi.bisect(latA, value, key=keyfunc)
                if len(latA) > 0 and latA[index-1].value == value:
                    latA[index-1].count += 1
                else:
                    latA.insert(index, DistSq(value, 1))
            else:
                # check if the value is already in the list
                index = bi.bisect(latB, value, key=keyfunc)
                if len(latB) > 0 and latB[index-1].value == value:
                    latB[index-1].count += 1
                else:
                    latB.insert(index, DistSq(value, 1))
    return latA, latB
# print("A:", [(x.count, (x.value)) for x in latA])
# print("B:", [(x.count, (x.value)) for x in latB])

def imbalance_nup(L, SL, SM, init_cond):
    upList, downList = init_cond
    f = SL - SM
    imb = 0
    denom = L*L
    for x in range(L):
        for y in range(L):
            for gamma in range(2):
                index = x + y * L
                if (x + y) % 2 == 0:
                    b = -1
                else:
                    b = 1

                if gamma == 0 and index in upList:
                    p = -1
                elif gamma == 1 and index in downList:
                    p = -1
                else:
                    p = 1
                
                # print(b, p, f)
                imb += b * p * f
    
    imb /= denom

    return imb


                
def main():
    Lx = Ly = 750
    L = 100
    latA, latB = get_series(Lx, Ly)
    latAvals = np.array([np.sqrt(x.value) for x in latA])
    latAcounts = np.array([x.count for x in latA])
    latBvals = np.array([np.sqrt(x.value) for x in latB])
    latBcounts = np.array([x.count for x in latB])

    phi = np.linspace(0.1, 10, 100)
    # phi = np.array([1.0])
    SL = np.empty(phi.shape[0])
    SM = np.empty(phi.shape[0])
    I_nup = np.empty(phi.shape[0])

    for i, p in tqdm(enumerate(phi)):
        sL, sM = get_series_sum(p, latAvals, latAcounts, latBvals, latBcounts)
        SL[i] = sL
        SM[i] = sM
        f = sL - sM
        psi = sL + sM
        init_cond = putils.get_initial_condition("altn_altupdown_updown", L)
        I_nup[i] = imbalance_nup(L, SL[i], SM[i], init_cond)

    print("L for series:", Lx, Ly)
    print("L for imbalance:", L)
    plt.plot(phi, I_nup, label=r"$I_{n_{\uparrow}}$")
    print(I_nup)
    # print("SL:", SL)
    # print("SM:", SM)

    plt.figure()
    plt.plot(phi, SL + SM, label="SL+SM")
    # plt.plot(phi, SL, label="SL")
    # plt.plot(phi, SM, label="SM")
    # plt.xlabel(r"$\phi$")
    # plt.ylabel(r"$S$")
    # plt.legend()
    # plt.savefig("plots/PNGs/square_lattice_series_sum.png")


    # plt.figure()
    # plt.plot(phi, SL-SM, label="Diff")
    # # plt.plot(phi, (SL-SM)/(SL+SM))
    # plt.xlabel(r"$\phi$")
    # plt.ylabel(r"$S_L - S_M$")
    # plt.legend()
    # plt.savefig("plots/PNGs/square_lattice_series_diff.png")
    

if __name__ == "__main__":
    main()
    plt.show()