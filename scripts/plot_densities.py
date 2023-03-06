#!python

# Author: Aditya Chincholi
# License: All Rights Reserved
# Purpose: To plot the particle spin densities

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm

import scripts.plot_utils as putils


def convert_to_2d(array, length):
    matrix = np.empty((length, length), dtype=array.dtype)
    for x in range(length):
        for y in range(length):
            matrix[x, y] = array[x + y*length]
    return matrix


def get_final_densities(GR_GRstar, binnum, errorExists=False, variance=None):
    n_up = 0.5 * (1 - GR_GRstar[0, 0, binnum])
    n_down = 0.5 * (1 - GR_GRstar[1, 1, binnum])
    S_plus = -0.5 * GR_GRstar[1, 0, binnum]
    S_minus = -0.5 * GR_GRstar[0, 1, binnum]
    densities = (n_up, n_down, S_plus, S_minus)
    
    if errorExists and variance is None:
        raise TypeError("if errorExists is True, variance must not be None")
    
    if errorExists and variance is not None:
        n_up_var = 0.25 * variance[0, 0, binnum]
        n_down_var = 0.25 * variance[1, 1, binnum]
        S_plus_var = 0.25 * variance[1, 0, binnum]
        S_minus_var = 0.25 * variance[0, 1, binnum]
        density_vars = (n_up_var, n_down_var, S_plus_var, S_minus_var)
        return densities, density_vars

    return densities

def plot_set(densities, title="", colorbar_kwargs={}):
    n_up, n_down, S_plus, S_minus = densities
    figures = []

    # Plot n_up and n_down
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4.5))
    figures.append(fig)
    fig.suptitle(title)

    vmin = min(n_up.min(), n_down.min())
    vmax = max(n_up.max(), n_down.max())
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    ax = axes[0]
    ax.set_title(r"$\langle n_{\uparrow} \rangle$")
    mappable = ax.imshow(n_up, origin="lower", norm=norm)
    fig.colorbar(mappable, ax=ax, **colorbar_kwargs)

    ax = axes[1]
    ax.set_title(r"$\langle n_{\downarrow} \rangle$")
    mappable = ax.imshow(n_down, origin="lower", norm=norm)
    fig.colorbar(mappable, ax=ax, **colorbar_kwargs)

    fig.tight_layout()

    # Plot S_plus and S_minus
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 4.5))
    figures.append(fig)
    fig.suptitle(title)

    vmin = min(np.abs(S_plus).min(), np.abs(S_minus).min())
    vmax = max(n_up.max(), n_down.max())
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    ax = axes[0, 0]
    ax.set_title(r"$|\langle S^+ \rangle|$")
    mappable = ax.imshow(np.abs(S_plus), origin="lower",
                        cmap="RdBu", norm=norm)
    fig.colorbar(mappable, ax=ax, **colorbar_kwargs)

    ax = axes[0, 1]
    ax.set_title(r"$|\langle S^- \rangle|$")
    mappable = ax.imshow(np.abs(S_minus), origin="lower",
                        cmap="RdBu", norm=norm)
    fig.colorbar(mappable, ax=ax, **colorbar_kwargs)

    vmin = min(np.angle(S_plus).min(), np.angle(S_minus).min())
    vmax = max(n_up.max(), n_down.max())
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    ax = axes[1, 0]
    ax.set_title(r"$\theta [\langle S^+ \rangle]$")
    mappable = ax.imshow(np.angle(S_plus), origin="lower",
                        cmap="RdBu", norm=norm)
    fig.colorbar(mappable, ax=ax, **colorbar_kwargs)

    ax = axes[1, 1]
    ax.set_title(r"$\theta [\langle S^- \rangle]$")
    mappable = ax.imshow(np.angle(S_minus), origin="lower",
                        cmap="RdBu", norm=norm)
    fig.colorbar(mappable, ax=ax, **colorbar_kwargs)

    fig.tight_layout()
    return figures


def calculate(kwargs, totalbins, prefix, pattern, errorExists=False):
    GR_GRstar = np.empty((2,2,totalbins+1), dtype=object)
    if errorExists:
        GR_GRstar_var = np.empty((2,2,totalbins+1), dtype=object)
    for alpha in range(2):
        for beta in range(2):
            kwargs["alpha"] = alpha
            kwargs["beta"] = beta
            kwargs["binnum"] = -1
            params = putils.SystemParams(**kwargs)
            filename = putils.getFilename(params, endname=f"_{pattern}.dat",
                                        prefix=prefix)
            varfilename = filename + ".variance"
            # print(filename)
            if(alpha == beta):
                matrix = convert_to_2d(np.loadtxt(filename), params.size)
                GR_GRstar[alpha, beta, totalbins] = matrix
                if errorExists:
                    var = convert_to_2d(np.loadtxt(varfilename), params.size)
                    GR_GRstar_var[alpha, beta, totalbins] = var
                
            else: #(alpha != beta)
                matrix = convert_to_2d(np.loadtxt(filename, dtype=complex),
                                    params.size)
                GR_GRstar[alpha, beta, totalbins] = matrix
                if errorExists:
                    var = convert_to_2d(np.loadtxt(varfilename, dtype=complex), params.size)
                    GR_GRstar_var[alpha, beta, totalbins] = var


            for binnum in range(totalbins):
                kwargs["alpha"] = alpha
                kwargs["beta"] = beta
                kwargs["binnum"] = binnum
                params = putils.SystemParams(**kwargs)
                filename = putils.getFilename(params, endname=f"_{pattern}.dat",
                                            prefix=prefix)
                varfilename = filename + ".variance"
                # print("varfilename:", varfilename)
                # print(filename)
                if(alpha == beta):
                    matrix = convert_to_2d(np.loadtxt(filename), params.size)
                    GR_GRstar[alpha, beta, binnum] = matrix
                    if errorExists:
                        var = convert_to_2d(np.loadtxt(varfilename), params.size)
                        GR_GRstar_var[alpha, beta, binnum] = var
                else: # (alpha != beta)
                    matrix = convert_to_2d(np.loadtxt(filename, dtype=complex),
                                        params.size)
                    GR_GRstar[alpha, beta, binnum] = matrix
                    if errorExists:
                        var = convert_to_2d(np.loadtxt(varfilename, dtype=complex),
                                            params.size)
                        GR_GRstar_var[alpha, beta, binnum] = var

    if errorExists:
        return GR_GRstar, GR_GRstar_var
    return GR_GRstar


def plot_all(params, totalbins, pattern, GR_GRstar, colorbar_kwargs):
    for binnum in tqdm(range(totalbins)):
        densities = get_final_densities(GR_GRstar, binnum)
        fig1, fig2 = plot_set(densities, title=f"Bin {binnum}", colorbar_kwargs=colorbar_kwargs)
        params.alpha = params.beta = None
        params.binnum = binnum
        fname1 = putils.getFilename(params, endname=f"_{pattern}_n_updown", prefix="mbl_densities")
        fig1.savefig("plots/PDFs/"+fname1+".pdf")
        fig1.savefig("plots/PNGs/"+fname1+".png")
        fname2 = putils.getFilename(params, endname=f"_{pattern}_S_plusminus", prefix="mbl_densities")
        fig2.savefig("plots/PDFs/"+fname2+".pdf")
        fig2.savefig("plots/PNGs/"+fname2+".png")
        plt.close(fig1)
        plt.close(fig2)
        
    densities = get_final_densities(GR_GRstar, totalbins)
    fig1, fig2 = plot_set(densities, title="Full Energy Window", colorbar_kwargs=colorbar_kwargs)
    params.alpha = params.beta = None
    params.binnum = -1
    fname1 = putils.getFilename(params, endname=f"_{pattern}_n_updown", prefix="mbl_densities")
    fig1.savefig("plots/PDFs/"+fname1+".pdf")
    fig1.savefig("plots/PNGs/"+fname1+".png")
    fname2 = putils.getFilename(params, endname=f"_{pattern}_S_plusminus", prefix="mbl_densities")
    fig2.savefig("plots/PDFs/"+fname2+".pdf")
    fig2.savefig("plots/PNGs/"+fname2+".png")
    plt.close(fig1)
    plt.close(fig2)


def main():
    colorbar_kwargs = {
        "location": "right",
        "fraction": 0.046,
        "pad": 0.04,
        "format": "%.2e"
    }

    prefix = "data/mbl_density"
    pattern = "alt_up_down"

    for c in [0.3, 0.4, 0.5]:
        for w in [13, 14, 15]:
            print(f"C = {c} W = {w}")
            totalbins = 50
            kwargs = {
                "size": 60,
                "coupling": c,
                "disorder": w,
                "hopup": 1.0,
                "hopdn": 1.0,
                "runs": 100,
                "nospin": False,
                "binnum": 2,
                "alpha": 0,
                "beta": 0
            }
            params = putils.SystemParams(**kwargs)
            GR_GRstar = calculate(kwargs, totalbins, prefix, pattern)
            plot_all(params, totalbins, pattern, GR_GRstar, colorbar_kwargs)

if __name__ == "__main__":
    main()
