import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns
# import scripts.plot_utils as putils


def main():
    sns.set()
    sns.set_style("white")
    sns.set_context("paper")
    # plt.rcParams.update({'font.size': 16})
    # sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    filename = "data/loc_len_vals_70percent.csv"
    # df = putils.createDataFrame(filename, minFraction=0.6)
    df = pd.read_csv(filename)
    plotHeatMap(df, 40)
    spins = ["upup", "dndn", "updn", "dnup"]
    for spin in spins:
        plotHeatMap(df.loc[df["spin"] == spin], 40, spins=spin)

def organize_ticks(numticks, plot, ax, data, axis="x"):
    if axis.lower() == "x":
        plotted_ticks = plot.get_xticklabels()
        pos = 0
    elif axis.lower() == "y":
        plotted_ticks = plot.get_yticklabels()
        pos = 1
    else:
        raise ValueError("axis must be either 'x' or 'y'")


    totalticks = len(data)  
    skip = (totalticks) // (numticks - 1)

    toShow = [0, totalticks-1]
    for n in range(skip, totalticks-1):
        if n % skip == 0:
            toShow.append(n)
    toShow.sort()

    ticks = data[toShow]
    tickpos = [plotted_ticks[i].get_position()[pos] for i in toShow]
    ticklabels = [f"{coup:.1f}" for coup in ticks]

    if axis == "x":
        ax.set_xticks(tickpos, labels=ticklabels)
    if axis == "y":
        ax.set_yticks(tickpos, labels=ticklabels)


def plotHeatMap(df, size, spins=None):
    fig, ax = plt.subplots(figsize=(5, 4))
    pivoted = df.pivot_table(index="disorder", columns="coupling", values="xi")

    # We need to sort the disorder values so that they increase as we go
    # from bottom to top (as it is on the y-axis)
    pivoted = pivoted.sort_values(by="disorder", ascending=False)

    cbar_kws = {"shrink": .85, "fraction": 0.05}
    plot = sns.heatmap(data=pivoted, annot=False, fmt=".1f",
                        vmin=0, vmax=40, cmap="icefire", ax=ax, cbar=True,
                        cbar_kws=cbar_kws)
    title = f"Loc Lengths for {size}x{size}"
    if spins is not None:
        title += f" s = {spins}"
    # plot.set(title=title)
    plot.set_xlabel(r"$\alpha$", fontsize="xx-large")
    plot.set_ylabel(r"$W$", fontsize="xx-large", rotation=0,
                    labelpad=15)

    plt.annotate("Localized\nPhase", (1.3, 2.5), color="white",
                fontsize="xx-large")

    plt.annotate("Delocalized\nPhase", (11.5, 10), color="black",
                fontsize="xx-large")


    organize_ticks(5, plot, ax, pivoted.columns, axis="x")
    organize_ticks(5, plot, ax, pivoted.index, axis="y")

    ax.xaxis.set_tick_params(rotation=0, labelsize="x-large")
    ax.yaxis.set_tick_params(rotation=0, labelsize="x-large")
    cbar_ax = fig.axes[-1]
    cbar_ax.tick_params(labelsize="x-large")

    
    fig.tight_layout()
    filename = "loc_lens_disorder_vs_coupling_heatmap_70p_new"
    if spins is not None:
        filename += f"_{spins}"
    fig.savefig("plots/PDFs/" + filename + ".pdf")
    fig.savefig("plots/PNGs/" + filename + ".png")


if __name__ == "__main__":
    main()
    plt.show()
