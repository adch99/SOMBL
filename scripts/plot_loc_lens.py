import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic


def getData(params):
    basename = f"data/mbl{params['length']}x{params['width']}" \
                + f"_W{params['disorder']:.2g}" \
                + f"_C{params['coupling']:.2g}" \
                + f"_T{params['hopping']:.2g}_"
    eigvals = np.loadtxt(basename+"eigvals.dat", np.float64, delimiter=" ")
    loc_lens = np.loadtxt(basename+"loclens.dat", np.float64, delimiter=" ")

    print(f"eigvals.shape: {eigvals.shape}")
    print(f"loc_lens.shape: {loc_lens.shape}")
    return eigvals, loc_lens

def processData(data, params):
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

def plotData(proc_data, params):
    bins, bin_edges, hist, mean_loc_lens, std_loc_lens = proc_data

    fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True)
    
    fig.suptitle(f"Spin Orbit Coupled MBL {params['length']}x{params['width']}")

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
    fig.savefig(f"plots/mbl{params['length']}x{params['width']}.pdf")
    fig.savefig(f"plots/mbl{params['length']}x{params['width']}.png")

def main(params):
    data = getData(params)
    proc_data = processData(data, params)
    plotData(proc_data, params)

if __name__ == "__main__":
    params = {
        "length": 30,
        "width": 30,
        "disorder": 15,
        "coupling": 0.0,
        "hopping": 1.0
    }
    main(params)
    plt.show()
