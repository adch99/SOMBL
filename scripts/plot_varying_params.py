import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt

params = {
    "coupling_const":   0.0,
    "hop_strength":     1.0,
    "disorder_vals":    np.linspace(5, 20, 11),
    "size":             10,
    "num_runs":         100,
    "nospin":           True
}

def getFilename(params):
    W_min = params['disorder_vals'][0]
    W_max = params['disorder_vals'][-1]

    basename = "mbl_nospin" if params["nospin"] else "mbl"
    basename += f"_{params['size']}x{params['size']}"
    basename += f"_W{W_min}to{W_max}"
    basename += f"_t{params['hop_strength']}"
    basename += f"_c{params['coupling_const']}"
    basename += f"_n{params['num_runs']}"
    basename += "_Xi_vs_W"
    return basename


def getData(params):

    loc_lens = np.zeros(params['disorder_vals'].shape)
    print(f"{'W':5} {'Xi':11} {'Time':5}")

    for i, disorder_strength in enumerate(params['disorder_vals']):
        args = ["./build/exact_diag_simulation",
                "-s", f"{params['size']}",
                "-c", f"{params['coupling_const']}",
                "-t", f"{params['hop_strength']}",
                "-w", f"{disorder_strength}",
                "-n", f"{params['num_runs']}",
                "-p"]
        start_time = time.time()
        result = subprocess.run(args, capture_output=True, text=True)
        result.check_returncode()
        time_taken = time.time() - start_time
        loc_lens[i] = float((result.stdout.split("\n")[-2]).split(":")[1])
        print(f"{disorder_strength:-5} {loc_lens[i]:-.5e} {time_taken:-.3f}")

    return params['disorder_vals'], loc_lens

def outputData(data, outfilename):
    np.savetxt("data/" + outfilename + ".dat", np.array(data).T)
    
def plotData(disorder_vals, loc_lens, name=None):
    fig, ax = plt.subplots()
    ax.plot(disorder_vals, loc_lens, marker=".")
    ax.set_ylim(0, 100)
    ax.set_xlim(0, 21)
    ax.set_xlabel("W")
    ax.set_ylabel(r"$\langle\xi_{loc}\rangle_{disorder}$")

    if name is None:
        name = "loc_lens_vs_W"
    fig.savefig("plots/" + name + ".pdf")
    fig.savefig("plots/" + name + ".png")

def main(params):
    data = getData(params)
    plotData(*data, name=getFilename(params))

if __name__ == "__main__":
    main(params)
    plt.show()