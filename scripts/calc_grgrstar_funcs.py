import numpy as np
import pandas as pd
from scipy.fft import fft2, fftfreq, fftshift
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

def calculate_imbalance_twoset(density, variance, length, setA, setB, norm=None):
    imb = 0
    imb_var = 0
    # sorted_setA = sorted(setA)
    for index in setB:
        x = index % length
        y = index // length
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


def calculate_staggered_magnetization(S, variance, length, setA, setB):
    mag = 0 + 0j
    mag_var = 0 + 0j

    for index in setB:
        x = index % length
        y = index // length
        mag -= S[x, y]
        mag_var += variance[x, y]

    for index in setA:
        x = index % length
        y = index // length
        mag += S[x, y]
        mag_var += variance[x, y]

    norm = length*length / 4
    mag /= norm
    mag_var /= norm**2

    return mag, mag_var

def get_oneset_mask(length, setA):
    mask = np.full((length, length), -1)
    setAX = [(index % length) for index in setA]
    setAY = [(index // length) for index in setA]
    mask[setAX, setAY] = 1
    var_mask = np.abs(mask)
    return mask, var_mask


def get_twoset_mask(length, setPlus, setMinus):
    mask = np.full((length, length), 0)
    setPlusX = [(index % length) for index in setPlus]
    setPlusY = [(index // length) for index in setPlus]
    setMinusX = [(index % length) for index in setMinus]
    setMinusY = [(index // length) for index in setMinus]

    mask[setPlusX, setPlusY] = 1
    mask[setMinusX, setMinusY] = -1
    var_mask = np.abs(mask)
    return mask, var_mask

def calc_imb_mask(density, variance, mask, var_mask, norm):
    imb = np.sum(density * mask) / norm
    imb_var = np.sum(variance * var_mask) / norm**2
    return imb, imb_var

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

def calculate_spin_imbalances_twoset_df(pattern, length, pspace):
    upList, downList = putils.get_initial_condition(pattern, length)
    # For charge imbalance
    if upList is None:
        upList = []
    if downList is None:
        downList = []
    
    setPlus = []
    setMinus = []
    for x in range(length):
        for y in range(length):
            if ((x + y) % 4 == 1):
                setPlus.append(x + y*length)
            elif ((x + y) % 4 == 3):
                setMinus.append(x + y*length)
            else:
                continue

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
        S_z = (n_up - n_down) / 2
        charge = n_up + n_down

        S_x_var = (S_plus_var + S_minus_var) / 4
        S_y_var = (S_plus_var + S_minus_var) / 4
        S_z_var = (n_up_var + n_down_var) / 4
        charge_var = n_up_var + n_down_var
        
        norm = (len(setPlus) + len(setMinus)) / 2

        mask, var_mask = get_twoset_mask(length, setPlus, setMinus)

        # imb_n_up = calculate_imbalance_twoset(n_up, n_up_var, length, setPlus, setMinus)
        # imb_n_down = calculate_imbalance_twoset(n_down, n_down_var, length, setPlus, setMinus)
        # imb_S_plus = calculate_imbalance_twoset(S_plus, S_plus_var, length, setPlus, setMinus, norm=norm)
        # imb_S_minus = calculate_imbalance_twoset(S_minus, S_minus_var, length, setPlus, setMinus, norm=norm)
        # imb_charge = calculate_imbalance_twoset(charge, charge_var, length, setPlus, setMinus)
        # imb_S_x = calculate_imbalance_twoset(S_x, S_x_var, length, setPlus, setMinus, norm=norm)
        # imb_S_y = calculate_imbalance_twoset(S_y, S_y_var, length, setPlus, setMinus, norm=norm)
        # imb_S_z = calculate_imbalance_twoset(S_z, S_z_var, length, setPlus, setMinus, norm=norm)
        # st_mag_z = calculate_staggered_magnetization(S_z, S_z_var, length, setPlus, setMinus) 

        imb_n_up = calc_imb_mask(n_up, n_up_var, mask, var_mask, np.sum(n_up))
        imb_n_down = calc_imb_mask(n_down, n_down, mask, var_mask, np.sum(n_down))
        imb_S_plus = calc_imb_mask(S_plus, S_plus_var, mask, var_mask, norm)
        imb_S_minus = calc_imb_mask(S_minus, S_minus_var, mask, var_mask, norm)
        imb_charge = calc_imb_mask(charge, charge_var, mask, var_mask, norm)
        imb_S_x = calc_imb_mask(S_x, S_x_var, mask, var_mask, norm)
        imb_S_y = calc_imb_mask(S_y, S_y_var, mask, var_mask, norm)
        imb_S_z = calc_imb_mask(S_z, S_z_var, mask, var_mask, norm)
        st_mag_z = calc_imb_mask(S_z, S_z_var, mask, var_mask, norm)

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
            "staggered_mag_z": st_mag_z[0],

            "spin_up_imb_n_up_var": imb_n_up[1],
            "spin_up_imb_n_down_var": imb_n_down[1],
            "spin_up_imb_S_plus_var": imb_S_plus[1],
            "spin_up_imb_S_minus_var": imb_S_minus[1],
            "spin_up_imb_charge_var": imb_charge[1],
            "spin_up_imb_S_x_var": imb_S_x[1],
            "spin_up_imb_S_y_var": imb_S_y[1],
            "spin_up_imb_S_z_var": imb_S_z[1],
            "staggered_mag_z_var": st_mag_z[1]
        }
        data.append(datapoint)

    df = pd.DataFrame(data)
    return df


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
        S_z = (n_up - n_down) / 2
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
        imb_S_x = calculate_imbalance(S_x, S_x_var, length, mainSubLattice, norm=norm)
        imb_S_y = calculate_imbalance(S_y, S_y_var, length, mainSubLattice, norm=norm)
        imb_S_z = calculate_imbalance(S_z, S_z_var, length, mainSubLattice, norm=norm)
        
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

def fourier_transform(density, length):
    density_k = fftshift(fft2(density, (length, length), norm="ortho"))
    kx = fftshift(fftfreq(length)) * 2 * np.pi
    ky = fftshift(fftfreq(length)) * 2 * np.pi
    return kx, ky, density_k

def fourier_component(kx, ky, density, length, phase=0):
    x = y = np.arange(0, length)
    X, Y = np.meshgrid(x, y)
    weight = np.exp(-1j*(kx*X + ky*Y + phase))
    val = np.sum(density * weight) / length
    return val

def main():
    pattern = "altn_altupdown_updown"
    length = 100
    pspace = [(c, w) for c in np.arange(0, 2.01, 0.1) for w in range(8, 18+1)]
    pspace += [(c, w) for c in np.arange(2.0, 3.01, 0.1) for w in range(12, 18+1)]
    pspace = [x for x in pspace if x != (0.8, 11)]
    spindf = calculate_spin_imbalances_twoset_df(pattern, length, pspace)
    spindf.to_csv(f"data/spin_imbalances_twoset_error_sample_L{length}_{pattern}.dat")

if __name__ == "__main__":
    main()