import numpy as np
import matplotlib.pyplot as plt

params = {
    "coupling_const":   0.0,
    "hop_strength":     1.0,
    "disorder":         10,
    "size":             40,
    "num_runs":         100,
    "nospin":           True
}

def getData(params):
    filename = "data/" + getFilename(params) + ".dat"
    data = np.loadtxt(filename)
    return data.T

def getFilename(params):
    basename = "mbl_nospin" if params["nospin"] else "mbl"
    basename += f"_{params['size']}x{params['size']}"
    basename += f"_W{params['disorder']}"
    basename += f"_T{params['hop_strength']}"
    basename += f"_C{params['coupling_const']}"
    basename += f"_N{params['num_runs']}"
    return basename

def fitData(data, params):
    return None

def plotData(data, params, fit):
    fig, ax = plt.subplots()
    # Scatter plot the data
    dists = data[0, :]
    loggfuncsq = np.log(data[1, :])
    ax.plot(dists, gfuncsq, marker="o", label="sim")
    ax.set_xlabel(r"$|r_i - r_j|$")
    ax.set_ylabel(r"$ln\left[G^2(|r_i - r_j|)\right]$")
    ax.legend()

    # Plot the fit

    # Print the residual

    filename = "plots/" + getFilename(params)
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")

def main():
    data = getData(params)
    fit = fitData(data, params)
    plotData(data, params, fit)

if __name__ == "__main__":
    main()
    plt.show()
