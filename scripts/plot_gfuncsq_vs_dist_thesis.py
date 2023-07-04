import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import ticker
# import pandas as pd
import scripts.plot_utils as putils


def main():
    sns.set()
    sns.set_style("ticks")
    sns.set_context("paper")
    # plt.rcParams.update({'font.size': 22})
    title=""
    outfilename = "delocalized_and_intermediate_phase_example_thesis"
    spin = "upup"
    numPhasePoints = 2
    disorders = np.array([13, 13])
    couplings = np.array([1.0, 1.5])
    emphasis = [False, ]*numPhasePoints
    # emphasis[2] = True
    # phasePointNum = 0
    fig, ax, palette = createCanvas(numPhasePoints, title=title)
    for phasePointNum in range(numPhasePoints):
        kwargs = {
            "size": 40,
            "coupling": couplings[phasePointNum],
            "disorder": disorders[phasePointNum],
            "hopup": 1.0,
            "hopdn": 1.0,
            "runs": 100,
            "nospin": False
        }
        params = putils.SystemParams(**kwargs)
        data = putils.getData(params, spin=spin)
        cleanedData = putils.cleanData(data)
        fitData = putils.fitData(cleanedData, params)
        dataColor = palette[2*phasePointNum+1]  # Use with "Paired"
        fitColor = palette[2*phasePointNum]  # Use with "Paired"
        # dataColor = palette[phasePointNum]
        # fitColor = palette[phasePointNum]
        plotPhasePoint(params, spin, cleanedData, fitData, ax,
                    color=dataColor, em=emphasis[phasePointNum])
        plotPhasePointFit(params, spin, cleanedData, fitData, ax,
                        color=fitColor)

        # if phasePointNum == 0:
        #     ax.text(5, 1e-6, r"Localized Phase $\alpha = 0.3$",
        #         fontsize="x-large", color=dataColor)
        # if phasePointNum == 1:
        #     ax.text(20, 1e-3, r"Delocalized Phase $\alpha = 1.5$",
        #         fontsize="x-large", color=dataColor)
        # if phasePointNum == 2:
        #     ax.text(25,4e-6, r"Intermediate Phase $\alpha = 1.0$",
        #         fontsize="x-large", color=dataColor)


    postPlot(fig, ax, filename=outfilename)


def createCanvas(numPhasePoints, title):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.set_xlabel(r"$r$", fontsize="xx-large")
    ax.set_ylabel(r"$\log_{10}\left\langle|G_R|^2(r)\right\rangle$",
                fontsize="xx-large", labelpad=15)
    ax.set_title(title)
    ax.set_yscale("log")
    ax.tick_params(axis="both", labelsize="x-large")
    ax.yaxis.set_major_formatter(ticker.LogFormatterExponent())
    palette = sns.color_palette(palette="Paired", n_colors=2*numPhasePoints)
    # palette = sns.color_palette(palette="tab10", n_colors=numPhasePoints)
    return fig, ax, palette


def postPlot(fig, ax, filename=None):
    # ax.legend()
    # ax.tick_params(rotation=0, size=25)

    fig.tight_layout()
    if filename is None:
        filename = "fitting_diffs_W10_C0.2to0.5_em0.5"
    fig.savefig("plots/PDFs/" + filename + ".pdf", bbox_inches="tight")
    fig.savefig("plots/PNGs/" + filename + ".png", bbox_inches="tight")


def plotPhasePoint(params, spin, cleanedData, fitData, ax,
                    color=None, em=True):
    dists, gfuncsq, errors = cleanedData
    exp, mant, resid, cutoff = fitData

    # Plot the points actually fit to the data
    # y = gfuncsq[cutoff:]
    # x = dists[cutoff:]
    y = gfuncsq
    x = dists
    # label = f"W = {params.disorder:.2g} α = {params.coupling:.2g} s={spin}"
    label = f"α = {params.coupling:.2g}"
    kwargs = {
        "label": label,
        "marker": ".",
        "linestyle": "",
        "linewidth": 2,
    }
    if color is not None:
        kwargs["markeredgecolor"] = color
        # kwargs["markerfacecolor"] = color
        kwargs["color"] = color
    # if em is False:
    #     kwargs["alpha"] = 0.5
    # if em is True:
    #     kwargs = {**kwargs,
    #         "marker": "o",
    #         "markersize": 13,
    #         # "markersize": 5,
    #         # "markerfacecolor": "white",
    #         # "markeredgewidth": 3
        # }
        # kwargs["zorder"] = 3
        # kwargs["markersize"] += 1
        # kwargs["linewidth"] += 1

    ax.plot(x, y, **kwargs)

    # # Plot the points excluded from the fit
    # y = gfuncsq[:cutoff]
    # x = dists[:cutoff]
    # # label = f"Excluded W = {params.disorder:.2g}"
    # # label +=" α = {params.coupling:.2g} s={spin}"
    # label = f"α = {params.coupling:.2g} Excluded"

    # kwargs = {
    #     "label": label,
    #     "marker": "s",
    #     "markersize": 13,
    #     "linestyle": "-",
    #     "linewidth": 2,
    #     "markerfacecolor": "white",
    #     "markeredgewidth": 3
    # }
    # if color is not None:
    #     kwargs["markeredgecolor"] = color
    #     kwargs["color"] = color
    # if em is False:
    #     kwargs["alpha"] = 0.5
    # if em is True:
    #     kwargs["zorder"] = 3
    #     kwargs["linewidth"] += 1

    # ax.plot(x, y, **kwargs)


def plotPhasePointFit(params, spin, cleanedData, fitData, ax, color=None):
    exp, mant, resid, cutoff = fitData
    dists, gfuncsq, errors = cleanedData

    x = np.linspace(dists.min(), dists.max(), 100)
    y = mant*np.exp(exp*x)

    # label = f"W = {params.disorder:.2g} α = {params.coupling:.2g} s={spin}"
    kwargs = {
        # "label": label
        "linewidth": 2,
        "alpha": 1,
        "zorder": 1,
    }
    if color is not None:
        kwargs["color"] = color

    ax.plot(x, y, **kwargs)


if __name__ == "__main__":
    main()
    plt.show()
