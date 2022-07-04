import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic


def getData(length, width):
    basename = f"data/mbl{length}x{width}_"
    eigvals = np.loadtxt(basename+"eigvals.dat", np.float64, delimiter=" ")
    loc_lens = np.loadtxt(basename+"loclens.dat", np.float64, delimiter=" ")

    print(f"eigvals.shape: {eigvals.shape}")
    print(f"loc_lens.shape: {loc_lens.shape}")
    return eigvals, loc_lens

def processData(data, length, width):
    eigvals, loc_lens = data
    print(f"min: {eigvals.min()}")
    eigvalsf = eigvals.ravel() # flattened version
    # mask = np.logical_and(eigvalsf > -1, eigvalsf < 1)
    # mask = (eigvalsf < 1)
    # eigvalsfm = eigvalsf[mask]
    # loc_lensfm = loc_lens.ravel()[mask] # flattened version
    eigvalsfm = eigvalsf
    loc_lensfm = loc_lens.ravel() # flattened version
    
    bins = 50
    hist, bin_edges = np.histogram(eigvalsfm, bins=bins)
    mean_result = binned_statistic(eigvalsfm, loc_lensfm, statistic="mean", bins=bins)
    std_result = binned_statistic(eigvalsfm, loc_lensfm, statistic="std", bins=bins)
    mean_loc_lens = mean_result[0]
    std_loc_lens = std_result[0]
    return (bins, bin_edges, hist, mean_loc_lens, std_loc_lens)

def plotData(proc_data, length, width):
    bins, bin_edges, hist, mean_loc_lens, std_loc_lens = proc_data

    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
    
    fig.suptitle(f"Spin Orbit Coupled MBL {length}x{width}")

    axes[0].hist(bin_edges[:-1], bin_edges, weights=hist)
    axes[0].set_xlabel("Energy")
    axes[0].set_ylabel("Disorder Avgd DOS")
    # axes[0].set_xscale("symlog")

    xvals = (bin_edges[:-1] + bin_edges[1:]) / 2
    axes[1].errorbar(xvals, mean_loc_lens, yerr=std_loc_lens, marker=".", capsize=2)
    axes[1].set_xlabel("Energy")
    axes[1].set_ylabel("Disorder Avgd Loc Length")
    # axes[1].set_xscale("symlog")
    # axes[1].set_yscale("log")

    # axes[1].set_xlim(-1, 1)

    fig.tight_layout()
    fig.savefig(f"plots/mbl{length}x{width}.pdf")
    fig.savefig(f"plots/mbl{length}x{width}.png")

def main(length, width):
    data = getData(length, width)
    proc_data = processData(data, length, width)
    plotData(proc_data, length, width)

if __name__ == "__main__":
    length = width = 20
    main(length, width)
    plt.show()
