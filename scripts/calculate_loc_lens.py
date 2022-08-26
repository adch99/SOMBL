#!/usr/bin/env python3
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt

def main():
    params = getParams()
    data = getData(params)
    fit_vals = fitData(data, params)
    plotData(fit_vals, data, params)
    exp, mant, resid, cutoff = fit_vals
    # print(f"xi:\t{-2/exp:e}")
    # print(f"residpp:\t{resid:e}")
    # print(f"cutoff:\t{cutoff}")
    output = {
        "xi": -2 / exp,
        "residpp": resid,
        "cutoff": cutoff
    }
    print(json.dumps(output))
    if not params.silent:
        plt.show()

def getParams():
    desc = """
    Reads the corresponding file from data/ and fits an exponential
    function to the G^2(r) data. Prints the localization length
    extracted and plots the fit and the datapoints.
    """
    parser = argparse.ArgumentParser(prog="calculate_loc_lens", description=desc)
    parser.add_argument("-s", "--size", help="Length and width of the lattice",
                        type=int, default=20)
    parser.add_argument("-c", "--coupling", help="Spin-orbit coupling constant",
                        type=float, default=0.0)
    parser.add_argument("-w", "--disorder", help="Strength of the disorder",
                        type=float, default=10.0)
    parser.add_argument("-t", "--hopping", help="Strength of the hopping",
                        type=float, default=1.0)
    parser.add_argument("-n", "--runs", help="Number of runs in the disorder average",
                        type=int, default=100)
    parser.add_argument("-p", "--nospin", help="Use a spinless model hamiltonian.",
                        action="store_true")
    parser.add_argument("--silent", help="Do not show the plot interactively.",
                        action="store_true")
    params = parser.parse_args()

    return params

def getData(params):
    filename = "data/" + getFilename(params) + ".dat"
    dists, gfuncsq = np.loadtxt(filename).T
    return dists, gfuncsq

def getFilename(params):
    basename = "mbl_nospin" if params.nospin else "mbl"
    basename += f"_{params.size}x{params.size}"
    basename += f"_W{params.disorder:.4g}"
    basename += f"_C{params.coupling:.4g}"
    basename += f"_T{params.hopping:.4g}"
    basename += f"_N{params.runs}"
    basename += "_distvsgfsq"
    return basename

def fitData(data, params):
    poly = np.polynomial.polynomial.Polynomial
    dists, gfuncsq = data
    num_points = dists.shape[0]
    diff = -1000
    residuals = 10000
    cutoff = -1
    while diff < 0 and (num_points - cutoff) >= 4: 
        # We need at least 3 points, hence the second condition
        cutoff += 1
        last_residuals = residuals
        x = dists[cutoff:]
        y = np.log(gfuncsq[cutoff:])
        series, extras = poly.fit(x, y, deg=1, full=True)
        residuals = extras[0][0] / (num_points - cutoff)
        diff = residuals - last_residuals
   
    c0, c1 = series.convert().coef
    exponent = c1
    mantissa = np.exp(c0)
    return exponent, mantissa, residuals, cutoff

def plotData(fit_vals, data, params):
    exp, mant, resid, cutoff = fit_vals
    dists, gfuncsq = data
    x = np.linspace(dists[cutoff], dists[-1], 100)
    y = mant * np.exp(exp*x)
    fig, axes = plt.subplots(figsize=(12, 8))
    if cutoff > 0:
        axes.scatter(dists[:cutoff], gfuncsq[:cutoff], label="Excluded")
    axes.scatter(dists[cutoff:], gfuncsq[cutoff:], label="Main Data")
    axes.plot(x, y)
    axes.set_yscale("log")
    axes.set_xlabel(r"$r$")
    axes.set_ylabel(r"$G_R^2(r)$")
    axes.set_title(f"Exponential Fit of y = {mant:.3e}*exp({exp:.3e}x)")
    axes.legend()

    filename = "plots/" + getFilename(params)
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")


if __name__ == "__main__":
    main()