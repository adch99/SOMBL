import numpy as np
import matplotlib.pyplot as plt

params = {
    "coupling_const":   0.0,
    "hop_strength":     1.0,
    "disorder":         10,
    "size":             10,
    "num_runs":         100,
    "nospin":           True
}

def getData(params):
    filename = "data/" + getFilename(params) + ".dat"
    

def getFilename(params):
    basename = "mbl_nospin" if params["nospin"] else "mbl"
    basename += f"_{params['size']}x{params['size']}"
    basename += f"_W{params['disorder']}"
    basename += f"_T{params['hop_strength']}"
    basename += f"_C{params['coupling_const']}"
    basename += f"_N{params['num_runs']}"
    return basename


def plotData(data, params, fitData):
    fig, ax = plt.subplots()
    # Scatter plot the data

    # Plot the fit

    # Print the residual