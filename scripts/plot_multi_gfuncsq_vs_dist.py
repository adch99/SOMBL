import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# import pandas as pd
import scripts.plot_utils as putils


def main():
    sns.set()
    sns.set_context("talk")
    # filename = "data/loc_len_vals_60percent.csv"
    # df = putils.createDataFrame(filename, minFraction=0.6)
    # df = pd.read_csv(filename)
    title = r"Qualitative Differences $W = 13$"
    outfilename = "qualitative_diffs_upup_W13_C0.3n1.0n1.5_em1.5"
    spin = "upup"
    numPhasePoints = 3
    disorders = np.array([13, 13, 13])
    couplings = np.array([0.3, 1.0, 1.5])
    emphasis = [False, False, False, False]
    emphasis[2] = True
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
        plotPhasePoint(params, spin, cleanedData, fitData, ax,
                        color=palette[2*phasePointNum+1], em=emphasis[phasePointNum])
        plotPhasePointFit(params, spin, cleanedData, fitData, ax, color=palette[2*phasePointNum])


    postPlot(fig, ax, filename=outfilename)


def createCanvas(numPhasePoints, title):
    fig, ax = plt.subplots(figsize=(16,9))
    ax.set_xlabel(r"$r$")
    ax.set_ylabel(r"$|G_R|^2(r)$")
    ax.set_title(title)
    ax.set_yscale("log")
    palette = sns.color_palette(palette="Paired", n_colors=2*numPhasePoints)
    return fig, ax, palette


def postPlot(fig, ax, filename=None):
    ax.legend()
    fig.tight_layout()
    if filename is None:
        filename = "fitting_diffs_W10_C0.2to0.5_em0.5"
    fig.savefig("plots/PDFs/" + filename + ".pdf")   
    fig.savefig("plots/PNGs/" + filename + ".png")   


def plotPhasePoint(params, spin, cleanedData, fitData, ax, color=None, em=True):
    dists, gfuncsq, errors = cleanedData
    exp, mant, resid, cutoff = fitData

    # Plot the points actually fit to the data
    y = gfuncsq[cutoff:]
    x = dists[cutoff:]
    # label = f"W = {params.disorder:.2g} α = {params.coupling:.2g} s={spin}"
    label = f"α = {params.coupling:.2g}"
    kwargs = {
        "label": label,
        "marker": "o",
        "markersize": 8,
        "linestyle": "-",
        "linewidth": 2,
        "markerfacecolor": "white",
        "markeredgewidth": 3
    }
    if color is not None:
        kwargs["markeredgecolor"] = color
        kwargs["color"] = color
    if em is False:
        kwargs["alpha"] = 0.5
    if em is True:
        kwargs["zorder"] = 3
        kwargs["markersize"] += 1
        kwargs["linewidth"] += 1

    ax.plot(x, y, **kwargs)


    # Plot the points excluded from the fit
    y = gfuncsq[:cutoff]
    x = dists[:cutoff]
    # label = f"Excluded W = {params.disorder:.2g} α = {params.coupling:.2g} s={spin}"
    label = f"α = {params.coupling:.2g} Excluded"

    kwargs = {
        "label": label,
        "marker": "o",
        "markersize": 8,
        "linestyle": "-",
        "linewidth": 2,
        "markerfacecolor": "black",
        "markeredgewidth": 3
    }
    if color is not None:
        kwargs["markeredgecolor"] = color
        kwargs["color"] = color
    if em is False:
        kwargs["alpha"] = 0.5
    if em is True:
        kwargs["zorder"] = 3
        kwargs["linewidth"] += 1

    ax.plot(x, y, **kwargs)

def plotPhasePointFit(params, spin, cleanedData, fitData, ax, color=None):
    exp, mant, resid, cutoff = fitData
    dists, gfuncsq, errors = cleanedData

    x = np.linspace(dists.min(), dists.max(), 100)
    y = mant*np.exp(exp*x)

    label = f"W = {params.disorder:.2g} α = {params.coupling:.2g} s={spin}"
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
