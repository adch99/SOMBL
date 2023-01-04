import numpy as np
import pandas as pd


class SystemParams(object):
    def __init__(self, size=20, coupling=0.0, disorder=10.0,
                hopup=1.0, hopdn=1.0, runs=100, nospin=False,
                batch=None, batchsize=None, binnum=None,
                alpha=None, beta=None):
        self.size = size
        self.coupling = coupling
        self.disorder = disorder
        self.hopup = hopup
        self.hopdn = hopdn
        self.runs = runs
        self.nospin = nospin
        self.batchsize = batchsize
        self.batch = batch
        self.binnum = binnum
        self.alpha = alpha
        self.beta = beta


def createDataFrame(filename, minFraction=0.6):
    dataList = []
    numCoup = 21
    couplings = np.linspace(0, 2, numCoup)
    numDis = 11
    disorder = np.linspace(8, 18, numDis)
    CC, WW = np.meshgrid(couplings, disorder, indexing="ij")
    for x in range(numCoup):
        for y in range(numDis):
            coup = CC[x, y]
            dis = WW[x, y]
            kwargs = {
                "size": 40,
                "coupling": coup,
                "disorder": dis,
                "hopup": 1.0,
                "hopdn": 1.0,
                "runs": 100,
                "nospin": False
            }
            params = SystemParams(**kwargs)
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
                fit_vals = fitData(data, params, minFraction=minFraction)
                exp, mant, resid, cutoff = fit_vals
                output = {
                    "xi": -2 / exp,
                    "residpp": resid,
                    "cutoff": cutoff,
                    "spin": spin
                }
                dataLine = {**kwargs, **output}
                print(dataLine)
                dataList.append(dataLine)

    df = pd.DataFrame(dataList)
    print(df)
    df.to_csv(filename)
    return df


def getFilename(params, spin=None, endname="_distvsgfsq", prefix=None):
    allowedSpins = [None, "upup", "updn", "dnup", "dndn"]
    if spin not in allowedSpins:
        raise TypeError(f"spin must be one of {allowedSpins}")
    if prefix is not None:
        basename = prefix
    else:
        basename = "mbl_nospin" if params.nospin else "mbl"
    basename += f"_{params.size}x{params.size}"
    basename += f"_W{params.disorder:.4g}"
    basename += f"_C{params.coupling:.4g}"
    # basename += f"_T{params.hopping:.4g}"
    basename += f"_TU{params.hopup:.4g}"
    basename += f"_TD{params.hopdn:.4g}"
    basename += f"_N{params.runs}"
    if params.batch is not None and params.batchsize is not None:
        basename += f"_BS{params.batchsize}"
        basename += f"_B{params.batch}"
    if params.alpha is not None:
        basename += f"_a{params.alpha}"
    if params.beta is not None:
        basename += f"_b{params.beta}"
    if params.binnum is not None:
        if params.binnum == -1:
            basename += "_full"
        else:
            basename += f"_bin{params.binnum}"
    if spin is not None:
        basename += f"_{spin}"
    basename += f"{endname}"
    return basename


def getData(params, spin=None):
    filename = "data/" + getFilename(params, spin=spin) + ".dat"
    dists, gfuncsq, errors = np.loadtxt(filename).T
    return dists, gfuncsq, errors


def cleanData(data, tol=1e-10):
    dists, gfuncsq, errors = data
    cond = (gfuncsq > tol)
    return dists[cond], gfuncsq[cond], errors[cond]


def fitData(data, params, minFraction=0.6, cov=False):
    poly = np.polynomial.polynomial.Polynomial
    dists, gfuncsq, errors = data
    num_points = dists.shape[0]

    diff = -1000
    cost = 10000
    cutoff = -1
    minPoints = np.ceil(minFraction * num_points)

    while diff < 0 and (num_points - cutoff) >= minPoints:
        # We need at least 3 points, hence the second condition
        cutoff += 1
        last_cost = cost
        x = dists[cutoff:]
        y = np.log(gfuncsq[cutoff:])
        series, extras = poly.fit(x, y, deg=1, full=True)  # , w=1/errors)

        used_points = num_points - cutoff
        cost = extras[0][0] / used_points  # + alpha*used_points
        diff = cost - last_cost

    c0, c1 = series.convert().coef
    exponent = c1
    mantissa = np.exp(c0)
    if cov:
        p, V = np.polyfit(x, y, deg=1, cov=True)
        return exponent, mantissa, cost, cutoff, V[0, 0] 
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
