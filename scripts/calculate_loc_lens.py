#!/usr/bin/env python3
import argparse
import json
import numpy as np
import matplotlib.pyplot as plt


def main():
    params = getParams()
    if params.nospin:
        spins = [None, ]
    else:
        spins = ["upup", "updn", "dnup", "dndn"]

    for spin in spins:
        data = getData(params, spin=spin)
        data = cleanData(data)
        num_points = data[0].shape[0]
        if num_points == 0:
            print(f"All points for spin = {spin} are too close to zero!")
            print("Cannot evaluate anything.\n")
            continue
        fit_vals = fitData(data, params)
        plotData(fit_vals, data, params, spin=spin)
        exp, mant, resid, cutoff = fit_vals
        output = {
            "xi": -2 / exp,
            "residpp": resid,
            "cutoff": cutoff,
            "spin": spin
        }
        print(json.dumps(output))

    if not params.silent:
        plt.show()


def getData(params, spin=None):
    filename = "data/" + getFilename(params, spin=spin) + ".dat"
    dists, gfuncsq, errors = np.loadtxt(filename).T
    return dists, gfuncsq, errors


def getFilename(params, spin=None):
    allowedSpins = [None, "upup", "updn", "dnup", "dndn"]
    if spin not in allowedSpins:
        raise TypeError(f"spin must be one of {allowedSpins}")
    basename = "mbl_nospin" if params.nospin else "mbl"
    basename += f"_{params.size}x{params.size}"
    basename += f"_W{params.disorder:.4g}"
    basename += f"_C{params.coupling:.4g}"
    # basename += f"_T{params.hopping:.4g}"
    basename += f"_TU{params.hopup:.4g}"
    basename += f"_TD{params.hopdn:.4g}"
    basename += f"_N{params.runs}"
    if spin is not None:
        basename += f"_{spin}"
    basename += "_distvsgfsq"
    return basename


def cleanData(data, tol=1e-10):
    dists, gfuncsq, errors = data
    cond = (gfuncsq > tol)
    return dists[cond], gfuncsq[cond], errors[cond]


def fitData(data, params):
    poly = np.polynomial.polynomial.Polynomial
    dists, gfuncsq, errors = data
    num_points = dists.shape[0]

    diff = -1000
    cost = 10000
    cutoff = -1
    minPoints = np.ceil(0.6 * num_points)

    while diff < 0 and (num_points - cutoff) >= minPoints:
        # We need at least 3 points, hence the second condition
        cutoff += 1
        last_cost = cost
        x = dists[cutoff:]
        y = np.log(gfuncsq[cutoff:])
        series, extras = poly.fit(x, y, deg=1, full=True)#, w=1/errors)
        
        used_points = num_points - cutoff
        cost = extras[0][0] / used_points #+ alpha*used_points
        diff = cost - last_cost

    c0, c1 = series.convert().coef
    exponent = c1
    mantissa = np.exp(c0)
    return exponent, mantissa, cost, cutoff


def checkData(data):
    dists, gfuncsq, errors = data
    tol = 1e-6
    probsum = 0
    dr = dists[1] - dists[0]
    for i in range(len(dists)):
        r = dists[i]
        probsum += r**2 * dr * gfuncsq[i]

    diff = abs(probsum - 1)

    if diff < tol:
        return True
    else:
        print("probsum:", probsum)
        return False


def plotData(fit_vals, data, params, spin=None):
    exp, mant, resid, cutoff = fit_vals
    dists, gfuncsq, errors = data
    x = np.linspace(dists[cutoff], dists[-1], 100)
    y = mant * np.exp(exp*x)
    fig, axes = plt.subplots(figsize=(12, 8))

    style_args = {
        "marker": "o",
        "capsize": 5,
        "linestyle": "none"
    }

    if cutoff > 0:
        data_x = dists[:cutoff]
        data_y = gfuncsq[:cutoff]
        data_yerr = errors[:cutoff]
        label = "Excluded"
        axes.errorbar(x=data_x, y=data_y, yerr=data_yerr,
                    label=label, **style_args)

    data_x = dists[cutoff:]
    data_y = gfuncsq[cutoff:]
    data_yerr = errors[cutoff:]
    label = "Main Data"
    # axes.errorbar(x=data_x, y=data_y, yerr=data_yerr,
    #             label=label, **style_args)
    axes.plot(data_x, data_y, label=label)
    # axes.plot(x, y)
    axes.set_yscale("log")
    axes.set_xscale("log")
    axes.set_xlabel(r"$r$")
    axes.set_ylabel(r"$G_R^2(r)$")

    title = f"s = {spin} "
    title += r"$\xi$" + f" = {-2/exp:.3e} "
    title += f"c={params.coupling} w={params.disorder} "
    title += f"l={params.size} n={params.runs}"

    axes.set_title(title)
    axes.legend()

    filename = "plots/" + getFilename(params, spin=spin)
    fig.savefig(filename + ".pdf")
    fig.savefig(filename + ".png")


def getParams():
    desc = """
    Reads the corresponding file from data/ and fits an exponential
    function to the G^2(r) data. Prints the localization length
    extracted and plots the fit and the datapoints.
    """
    parser = argparse.ArgumentParser(
                                    prog="calculate_loc_lens",
                                    description=desc)
    parser.add_argument(
                        "-s", "--size",
                        help="Length and width of the lattice",
                        type=int, default=20)
    parser.add_argument(
                        "-c", "--coupling",
                        help="Spin-orbit coupling constant",
                        type=float, default=0.0)
    parser.add_argument(
                        "-w", "--disorder",
                        help="Strength of the disorder",
                        type=float, default=10.0)
    parser.add_argument(
                        "-t", "--hopping",
                        help="Strength of the hopping",
                        type=float, default=1.0)
    parser.add_argument(
                        "-u", "--hopup",
                        help="Strength of the hopping of up spins",
                        type=float, default=1.0)
    parser.add_argument(
                        "-d", "--hopdn",
                        help="Strength of the hopping of down spins",
                        type=float, default=1.0)
    parser.add_argument(
                        "-n", "--runs",
                        help="Number of runs in the disorder average",
                        type=int, default=100)
    parser.add_argument(
                        "-p", "--nospin",
                        help="Use a spinless model hamiltonian.",
                        action="store_true")
    parser.add_argument(
                        "--silent",
                        help="Do not show the plot interactively.",
                        action="store_true")
    params = parser.parse_args()
    return params


if __name__ == "__main__":
    main()
