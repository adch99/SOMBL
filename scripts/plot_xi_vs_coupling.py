import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# import pandas as pd
import scripts.plot_utils as putils


def main():
    plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(figsize=(16,9))

    ax.set_xlabel(r"$\alpha$", fontsize="x-large")
    ax.set_ylabel(r"$\xi$", fontsize="x-large")
    ax.set_ylim(-0.5, 20.5)

    spin = "upup"
    for disorder in range(11, 17, 1):
        couplings = []
        loclens = []
        loclenerrs = []
        for coupling in np.arange(0, 2.1, 0.1):
            kwargs = {
                "size": 40,
                "coupling": coupling,
                "disorder": disorder,
                "hopup": 1.0,
                "hopdn": 1.0,
                "runs": 100,
                "nospin": False
            }
            params = putils.SystemParams(**kwargs)
            data = putils.getData(params, spin=spin)
            cleanedData = putils.cleanData(data)
            exp, mant, residpp, cutoff, experrsq = putils.fitData(cleanedData, params, cov=True)
            fitData = exp, mant, residpp, cutoff
            xi = -2/exp
            xi_error = 2 * np.sqrt(experrsq) / exp**2
            couplings.append(coupling)
            loclens.append(xi)
            loclenerrs.append(xi_error)
        
        kwargs = {
            "yerr": loclenerrs,
            "label": f"W = {disorder}",
            "marker": "o",
            "capsize": 5
        }


        ax.errorbar(couplings, loclens, **kwargs)

    # ax.plot(np.linspace(0,2,100), 40*np.ones(100), linestyle="--", color="black")
    # ax.plot(np.linspace(0,2,100), 20*np.ones(100), linestyle="--", color="blue",  linewidth=4)

    # ax.legend()
    fig.tight_layout()
    filename = "xi_vs_couplings_more_restricted"
    fig.savefig("plots/PDFs/" + filename + ".pdf")
    fig.savefig("plots/PNGs/" + filename + ".png")
    fig.savefig("plots/PDFs/" + filename + ".svg")



if __name__ == "__main__":
    main()
    plt.show()