# import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# import scripts.plot_utils as putils


def main():
    sns.set()
    sns.set_context("paper")
    # sns.set_context("paper", rc={"font.size":8,"axes.titlesize":8,"axes.labelsize":5})
    filename = "data/loc_len_vals_60percent.csv"
    # df = putils.createDataFrame(filename, minFraction=0.6)
    df = pd.read_csv(filename)
    plotHeatMap(df, 40)
    spins = ["upup", "dndn", "updn", "dnup"]
    for spin in spins:
        plotHeatMap(df.loc[df["spin"] == spin], 40, spins=spin)


def plotHeatMap(df, size, spins=None):
    fig, ax = plt.subplots(figsize=(12, 6))
    pivoted = df.pivot_table(index="disorder", columns="coupling", values="xi")
    plot = sns.heatmap(data=pivoted, annot=True, fmt=".1f",
                        vmin=0, vmax=40, cmap="icefire", ax=ax, cbar=True,
                        cbar_kws={"shrink": .85})
    title = f"Loc Lengths for {size}x{size}"
    if spins is not None:
        title += f" s = {spins}"
    plot.set(title=title)
    xticklabels = [f"{coup:.1f}" for coup in pivoted.columns]
    plot.set_xticklabels(xticklabels, rotation=0)
    plot.set_yticklabels(plot.get_yticklabels(), rotation=0)
    fig.tight_layout()
    filename = "loc_lens_disorder_vs_coupling_heatmap_60p"
    if spins is not None:
        filename += f"_{spins}"
    fig.savefig("plots/PDFs/" + filename + ".pdf")
    fig.savefig("plots/PNGs/" + filename + ".png")


if __name__ == "__main__":
    main()
    plt.show()
