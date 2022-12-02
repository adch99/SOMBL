# import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import seaborn as sns
# import scripts.plot_utils as putils


def main():
    sns.set()
    sns.set_style("white")
    sns.set_context("talk")
    plt.rcParams.update({'font.size': 16})
    # sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    filename = "data/loc_len_vals_60percent.csv"
    # df = putils.createDataFrame(filename, minFraction=0.6)
    df = pd.read_csv(filename)
    plotHeatMap(df, 40)
    spins = ["upup", "dndn", "updn", "dnup"]
    for spin in spins:
        plotHeatMap(df.loc[df["spin"] == spin], 40, spins=spin)


def plotHeatMap(df, size, spins=None):
    fig, ax = plt.subplots(figsize=(16, 9))
    pivoted = df.pivot_table(index="disorder", columns="coupling", values="xi")

    # We need to sort the disorder values so that they increase as we go
    # from bottom to top (as it is on the y-axis)
    pivoted = pivoted.sort_values(by="disorder", ascending=False)

    plot = sns.heatmap(data=pivoted, annot=False, fmt=".1f",
                        vmin=0, vmax=40, cmap="icefire", ax=ax, cbar=True,
                        cbar_kws={"shrink": .85})
    title = f"Loc Lengths for {size}x{size}"
    if spins is not None:
        title += f" s = {spins}"
    # plot.set(title=title)
    plot.set_xlabel(r"$\alpha$", fontsize="xx-large")
    plot.set_ylabel(r"$W$", fontsize="xx-large", rotation=0)
    xticklabels = [f"{coup:.1f}" for coup in pivoted.columns]
    plot.set_xticklabels(xticklabels, rotation=0, fontsize="large")
    plot.set_yticklabels(plot.get_yticklabels(), rotation=0, fontsize="x-large")
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    fig.tight_layout()
    filename = "loc_lens_disorder_vs_coupling_heatmap_60p"
    if spins is not None:
        filename += f"_{spins}"
    fig.savefig("plots/PDFs/" + filename + ".pdf")
    fig.savefig("plots/PNGs/" + filename + ".png")


if __name__ == "__main__":
    main()
    plt.show()
