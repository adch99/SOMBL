import subprocess
import numpy as np
import matplotlib.pyplot as plt


def getData():
    coupling_const = 0.0
    hop_strength = 1.0
    disorder_vals = np.linspace(10, 20, 11)
    size = 70
    loc_lens = np.zeros(disorder_vals.shape)

    for i, disorder_strength in enumerate(disorder_vals):
        args = ["./build/exact_diag_simulation",
                "-s", f"{size}",
                "-c", f"{coupling_const}",
                "-t", f"{hop_strength}",
                "-w", f"{disorder_strength}"]

        result = subprocess.run(args, capture_output=True, text=True)
        result.check_returncode()
        loc_lens[i] = float((result.stdout.split("\n")[-2]).split(":")[1])
        print(f"For W = {disorder_strength}, Xi = {loc_lens[i]}")

    return disorder_vals, loc_lens

def plotData(disorder_vals, loc_lens, name=None):
    fig, ax = plt.subplots()
    ax.plot(disorder_vals, loc_lens, marker=".")
    ax.set_xlabel("W")
    ax.set_ylabel(r"$\langle\Xi_{loc}\rangle_{disorder}$")

    if name is None:
        name = "plots/loc_lens_vs_W"
    fig.savefig(name + ".pdf")
    fig.savefig(name + ".png")

def main():
    data = getData()
    plotData(*data, name="plots/loc_lens_vs_W_30x30_C0_T1")

if __name__ == "__main__":
    main()
    plt.show()