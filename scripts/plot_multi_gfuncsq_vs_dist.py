import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# import pandas as pd
import scripts.plot_utils as putils


def main():
    sns.set()
    sns.set_style("whitegrid")
    sns.set_context("talk")
    plt.rcParams.update({'font.size': 22})
    # filename = "data/loc_len_vals_60percent.csv"
    # df = putils.createDataFrame(filename, minFraction=0.6)
    # df = pd.read_csv(filename)
    # title = r"Fitting Issues $s = \uparrow\uparrow$"
    # outfilename = "fitting_diffs_W10_C0.2to0.5_upup_em0.5.pdf"
    # title = r"Qualitative Differences $W = 13$"
    title=""
    outfilename = "localized_phase_example"
    spin = "upup"
    numPhasePoints = 1
    disorders = np.array([13])
    couplings = np.array([0.3])
    emphasis = [False, ]*numPhasePoints
    emphasis[0] = True
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

        if phasePointNum == 0:
            ax.text(5, 1e-6, r"Localized Phase $\alpha = 0.3$",
                fontsize="x-large", color=dataColor)
        if phasePointNum == 1:
            ax.text(20, 1e-3, r"Delocalized Phase $\alpha = 1.5$",
                fontsize="x-large", color=dataColor)
        if phasePointNum == 2:
            ax.text(25,4e-6, r"Intermediate Phase $\alpha = 1.0$",
                fontsize="x-large", color=dataColor)


    postPlot(fig, ax, filename=outfilename)


def createCanvas(numPhasePoints, title):
    fig, ax = plt.subplots(figsize=(16, 9))
    ax.set_xlabel(r"$r$", fontsize="x-large")
    ax.set_ylabel(r"$|G_R|^2(r)$", rotation=0, fontsize="x-large")
    ax.set_title(title)
    ax.set_yscale("log")
    palette = sns.color_palette(palette="Paired", n_colors=2*numPhasePoints)
    return fig, ax, palette


def postPlot(fig, ax, filename=None):
    # ax.legend()
    ax.tick_params(rotation=0, size=25)

    fig.tight_layout()
    if filename is None:
        filename = "fitting_diffs_W10_C0.2to0.5_em0.5"
    fig.savefig("plots/PDFs/" + filename + ".pdf")
    fig.savefig("plots/PNGs/" + filename + ".png")


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
        "marker": "o",
        "markersize": 8,
        # "markerfacecolor": "white",
        "markeredgewidth": 3,
        "linestyle": "-",
        "linewidth": 2,
    }
    if color is not None:
        kwargs["markeredgecolor"] = color
        kwargs["color"] = color
    if em is False:
        kwargs["alpha"] = 0.5
    if em is True:
        kwargs = {**kwargs,
            "marker": "o",
            "markersize": 13,
            "markerfacecolor": "white",
            "markeredgewidth": 3
        }
        kwargs["zorder"] = 3
        kwargs["markersize"] += 1
        kwargs["linewidth"] += 1

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
        "linewidth": 5,
        "alpha": 0.9,
        "zorder": 1,
    }
    if color is not None:
        kwargs["color"] = color

    ax.plot(x, y, **kwargs)


if __name__ == "__main__":
    main()
    plt.show()
