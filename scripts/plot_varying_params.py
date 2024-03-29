#!/usr/bin/env python3
import subprocess
import time
import json
import numpy as np
import matplotlib.pyplot as plt


def getParams(paramsFilename):
    params = {
        "coupling_const":   0.0,
        "hopup":            1.0,
        "hopdn":            1.0,
        "disorder_low":     5,
        "disorder_high":    20,
        "disorder_samples": 11,
        "size":             10,
        "num_runs":         100,
        "nospin":           False
    }
    return params


def getFilename(params):
    W_min = params['disorder_low']
    W_max = params['disorder_high']

    basename = "mbl_nospin" if params["nospin"] else "mbl"
    basename += f"_{params['size']}x{params['size']}"
    basename += f"_W{W_min}to{W_max}"
    # basename += f"_T{params['hop_strength']}"
    basename += f"_TU{params['hopup']}"
    basename += f"_TD{params['hopdn']}"
    basename += f"_C{params['coupling_const']}"
    basename += f"_N{params['num_runs']}"
    basename += "_Xi_vs_W"
    return basename


def runExactDiag(params):
    print("Starting exact diagonalizations")
    iterator = np.linspace(
                        params["disorder_low"],
                        params["disorder_high"],
                        params["disorder_samples"]
                        )
    for i, disorder_strength in enumerate(iterator):
        args = ["./build/exact_diag_simulation",
                "-s", f"{params['size']}",
                "-c", f"{params['coupling_const']}",
                # "-t", f"{params['hop_strength']}",
                "-u", str(params["hopup"]),
                "-d", str(params["hopdn"]),
                "-w", f"{disorder_strength}",
                "-n", f"{params['num_runs']}"]
        if params["nospin"]:
            args.append("-p")

        command = " ".join(args)

        print(
            f"DIAG: W = {disorder_strength:.2f} Started...",
            end="", flush=True)
        start_time = time.time()
        result = subprocess.run(
                                [command, ],
                                shell=True,
                                capture_output=True,
                                text=True
                                )
        try:
            result.check_returncode()
        except subprocess.CalledProcessError:
            print("\nCall to build/exact_diag_simulation failed!")
            print("PARAMS:", params)
            print("CALL: ", " ".join(args))
            print("STDERR:", result.stderr)
            print("STDOUT:", result.stdout)
            print("RETURNCODE:", result.returncode)
            exit(1)
        time_taken = time.time() - start_time
        print(f"Done in {time_taken:.1f}s")


def runFuncCalc(params):
    print("Starting calculation of G(r)")
    iterator = np.linspace(
                        params["disorder_low"],
                        params["disorder_high"],
                        params["disorder_samples"]
                        )
    for i, disorder_strength in enumerate(iterator):
        args = ["./build/calculate_dist_vs_gfuncsq",
                "-s", f"{params['size']}",
                "-c", f"{params['coupling_const']}",
                # "-t", f"{params['hop_strength']}",
                "-u", str(params["hopup"]),
                "-d", str(params["hopdn"]),
                "-w", f"{disorder_strength}",
                "-n", f"{params['num_runs']}"]
        if params["nospin"]:
            args.append("-p")

        print(
            f"FUNC: W = {disorder_strength:.2e} Started...",
            end="", flush=True)
        start_time = time.time()
        result = subprocess.run(args, capture_output=True, text=True)
        try:
            result.check_returncode()
        except subprocess.CalledProcessError:
            print("\nCALL TO build/calculate_dist_vs_gfuncsq FAILED!")
            print("PARAMS:", params)
            print("CALL: ", " ".join(args))
            print("STDERR:", result.stderr)
            print("STDOUT:", result.stdout)
            print("RETURNCODE:", result.returncode)
            exit(1)
        time_taken = time.time() - start_time
        print(f"Done in {time_taken:.1f}s")


def runLocLens(params):
    loc_lens = np.zeros(params["disorder_samples"])
    print(f"{'W':5} {'Xi':11} {'Time':5} {'Residpp':11} {'Cutoff':5}")

    iterator = np.linspace(
                        params["disorder_low"],
                        params["disorder_high"],
                        params["disorder_samples"]
                        )
    for i, disorder_strength in enumerate(iterator):
        args = ["./scripts/calculate_loc_lens.py",
                "-s", f"{params['size']}",
                "-c", f"{params['coupling_const']}",
                # "-t", f"{params['hop_strength']}",
                "-u", str(params["hopup"]),
                "-d", str(params["hopdn"]),
                "-w", f"{disorder_strength}",
                "-n", f"{params['num_runs']}",
                "--silent"]

        if params["nospin"]:
            args.append("-p")

        start_time = time.time()
        result = subprocess.run(args, capture_output=True, text=True)
        try:
            result.check_returncode()
        except subprocess.CalledProcessError:
            print("Call to scripts/calculate_loc_lens.py failed!")
            print("params:", params)
            print("stderr:", result.stderr)
            exit(1)
        time_taken = time.time() - start_time
        # data = json.loads(result.stdout)
        # loc_lens[i] = data["xi"]
        # output_string = f"{disorder_strength:-5} {loc_lens[i]:-.5e}"
        # output_string += f" {time_taken:-.3f} {data['residpp']:-5e}"
        # output_string += f" {data['cutoff']}"
        # print(output_string)
    return params["disorder_vals"], loc_lens


def outputData(data, outfilename):
    np.savetxt("data/" + outfilename + ".dat", np.array(data).T)


def plotData(disorder_vals, loc_lens, name=None):
    fig, ax = plt.subplots()
    ax.plot(disorder_vals, loc_lens, marker=".")
    # ax.set_ylim(0, 100)
    # ax.set_xlim(0, 21)
    ax.set_xlabel("W")
    ax.set_ylabel(r"$\langle\xi_{loc}\rangle_{disorder}$")

    if name is None:
        name = "loc_lens_vs_W"
    fig.savefig("plots/" + name + ".pdf")
    fig.savefig("plots/" + name + ".png")


def main(params):
    # runExactDiag(params)
    # print("")
    runFuncCalc(params)
    print("")
    data = runLocLens(params)
    plotData(*data, name=getFilename(params))


if __name__ == "__main__":
    params = getParams("data/params.json")
    print(params)
    main(params)
    plt.show()
