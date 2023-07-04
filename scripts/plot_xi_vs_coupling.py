import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# import pandas as pd
import scripts.plot_utils as putils


def main():
    # plt.rcParams.update({'font.size': 22})
    fig, ax = plt.subplots(figsize=(5,4))

    ax.set_xlabel(r"$\alpha$", fontsize="xx-large")
    # ax.set_ylabel(r"$\xi^{-1}$", fontsize="xx-large", rotation=0, labelpad=20)
    ax.set_ylabel(r"$\xi}$", fontsize="xx-large",
                rotation=0, labelpad=5)
    ax.tick_params(axis="both", labelsize="x-large")
    # ax.set_ylim(-0.5, 20.5)
    ax.set_ylim(-0.5, 40.5)

    spins = ["upup", "updn"]

    for spin in spins:
        # for i, disorder in enumerate(range(11, 17, 1)):
        for disorder in [13]:
            i = 0
            couplings = []
            loclens = []
            loclenerrs = []
            # locleninvs = []
            # locleninverrs = []
            for coupling in np.arange(0.1, 2.1, 0.1):
            # for coupling in [0.3, 1.0, 1.5]:
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
                # xi_inv = -exp / 2
                # xi_inv_error = 2 * np.sqrt(experrsq)
                # couplings.append(coupling)
                # locleninvs.append(xi_inv)
                # locleninverrs.append(xi_inv_error)

            if spin == "upup":
                s = 0
            if spin == "updn":
                s = 1

            kwargs = {
                "yerr": loclenerrs,
                # "yerr": locleninverrs,
                "label": f"W = {disorder} s = {spin}",
                # "marker": "o",
                "marker": ".",
                "capsize": 5,
                "linestyle": "-",
                "color": plt.cm.tab10(2*i+s)
            }


            ax.errorbar(couplings, loclens, **kwargs)
            # ax.errorbar(couplings, locleninvs, **kwargs)

    # ax.plot(np.linspace(0,2,100), 40*np.ones(100), linestyle="--", color="black")
    ax.plot(np.linspace(0,2,100), 20*np.ones(100), linestyle="--",
            color="green", linewidth=4)

    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    filename = "xi_vs_couplings_upup_vs_updn_thesis"
    fig.savefig("plots/PDFs/" + filename + ".pdf")
    fig.savefig("plots/PNGs/" + filename + ".png")
    fig.savefig("plots/PDFs/" + filename + ".svg")



if __name__ == "__main__":
    main()
    plt.show()