import numpy as np
import pandas as pd
import scripts.plot_utils as putils
import scripts.plot_densities as pdens
from tqdm import tqdm

def get_densities(length, coupling, disorder, binnum, pattern):
    prefix = "data/mbl_density"
    totalbins = 0
    kwargs = {
        "size": length,
        "coupling": coupling,
        "disorder": disorder,
        "hopup": 1.0,
        "hopdn": 1.0,
        "runs": 100,
        "nospin": False,
        "binnum": 2,
        "alpha": 0,
        "beta": 0
    }
    # params = putils.SystemParams(**kwargs)
    GR_GRstar, GR_GRstar_var = pdens.calculate(kwargs, totalbins, prefix, pattern, errorExists=True)
    # print("GRGR*var", GR_GRstar_var)
    densities, density_vars = pdens.get_final_densities(GR_GRstar, binnum, errorExists=True, variance=GR_GRstar_var)
    return densities, density_vars

def calculate_imbalance(density, variance, length, setA, norm=None):
    imb = 0
    imb_var = 0
    sorted_setA = sorted(setA)
    for x in range(length):
        for y in range(length):
            index = x + length*y
            if index not in setA:
                imb -= density[x, y]
                imb_var += variance[x, y]

    for index in setA:
        x = index % length
        y = index // length
        imb += density[x, y]
        imb_var += variance[x, y]

    if norm is None:
        imb /= np.sum(density)
        imb_var /= np.sum(density)**2
    else:
        imb /= norm
        imb_var /= norm**2

    return imb, imb_var

def calculate_spin_imbalances_df(pattern, length, pspace):
    upList, downList = putils.get_initial_condition(pattern, length)
    # For charge imbalance
    if upList is None:
        upList = []
    if downList is None:
        downList = []
    setAup = list(upList)
    setAdown = list(downList)
    
    mainSubLattice = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 2 == 0):
                mainSubLattice.append(x + y*length)


    data = []
    for coupling, disorder in tqdm(pspace):
        # print(f"α = {coupling} W = {disorder}")
        binnum = -1
        densities, density_vars = get_densities(length, coupling, disorder, binnum, pattern)
        # print(density_vars)
        n_up, n_down, S_plus, S_minus = densities
        n_up_var, n_down_var, S_plus_var, S_minus_var = density_vars
        
        S_x = (S_plus + S_minus)/2
        S_y = -1j*(S_plus - S_minus)/2
        S_z = n_up - n_down
        charge = n_up + n_down

        S_x_var = (S_plus_var + S_minus_var)/4
        S_y_var = (S_plus_var + S_minus_var)/4
        S_z_var = n_up_var + n_down_var
        charge_var = n_up_var + n_down_var
        
        norm = len(mainSubLattice) / 2

        imb_n_up = calculate_imbalance(n_up, n_up_var, length, mainSubLattice)
        imb_n_down = calculate_imbalance(n_down, n_down_var, length, mainSubLattice)
        imb_S_plus = calculate_imbalance(S_plus, S_plus_var, length, mainSubLattice, norm=norm)
        imb_S_minus = calculate_imbalance(S_minus, S_minus_var, length, mainSubLattice, norm=norm)
        imb_charge = calculate_imbalance(charge, charge_var, length, mainSubLattice)
        imb_S_x = calculate_imbalance(S_x, S_x_var, length, setAup, norm=norm)
        imb_S_y = calculate_imbalance(S_y, S_y_var, length, setAup, norm=norm)
        imb_S_z = calculate_imbalance(S_z, S_z_var, length, setAup, norm=norm)
        
        datapoint = {
            "coupling": coupling,
            "disorder": disorder,
            "binnum": binnum,
            "spin_up_imb_n_up": imb_n_up[0],
            "spin_up_imb_n_down": imb_n_down[0],
            "spin_up_imb_S_plus": imb_S_plus[0],
            "spin_up_imb_S_minus": imb_S_minus[0],
            "spin_up_imb_charge": imb_charge[0],
            "spin_up_imb_S_x": imb_S_x[0],
            "spin_up_imb_S_y": imb_S_y[0],
            "spin_up_imb_S_z": imb_S_z[0],

            "spin_up_imb_n_up_var": imb_n_up[1],
            "spin_up_imb_n_down_var": imb_n_down[1],
            "spin_up_imb_S_plus_var": imb_S_plus[1],
            "spin_up_imb_S_minus_var": imb_S_minus[1],
            "spin_up_imb_charge_var": imb_charge[1],
            "spin_up_imb_S_x_var": imb_S_x[1],
            "spin_up_imb_S_y_var": imb_S_y[1],
            "spin_up_imb_S_z_var": imb_S_z[1]

        }
        data.append(datapoint)

    df = pd.DataFrame(data)
    return df

def main():
    pattern = "altn_altupdown_updown"
    length = 60
    # pspace = [(c, w) for c in np.arange(0, 2.01, 0.1) for w in range(8, 18+1)]
    pspace = [(c, w) for c in np.arange(0, 3.01, 0.1) for w in range(8, 18+1)]
    # pspace += [(c, w) for c in np.arange(2.0, 3.01, 0.1) for w in range(12, 18+1)]
    spindf = calculate_spin_imbalances_df(pattern, length, pspace)
    spindf.to_csv(f"data/spin_imbalances_error_L{length}_{pattern}.dat")

if __name__ == "__main__":
    main()