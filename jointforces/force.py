import numpy as np
import pandas as pd
from tqdm import tqdm
from glob import glob
from natsort import natsorted
from .simulation import load_lookup_functions
from .utils import load

def reconstruct(folder, lookupfile, muperpixel, outfile=None, r_min=2):
    # get filepaths for PIV results
    dis_files = natsorted(glob(folder+'/dis*.npy'))
    seg_files = natsorted(glob(folder+'/seg*.npy'))

    # load lookup table
    get_displacement, get_pressure = load_lookup_functions(lookupfile)

    # initial spheroid radius and surface (used for force estimation)
    r0 = load(seg_files[0])['radius']
    A0 = 4 * np.pi * (r0 * (10 ** -6)) ** 2

    # initialize result dictionary
    results = {'pressure_mean': [], 'pressure_median': [], 'pressure_std': [],
               'contractility_mean': [], 'contractility_median': [], 'contractility_std': []}

    u_sum = None
    v_sum = None

    # loop over series of PIV results
    for (dis_file, seg_file) in tqdm(zip(dis_files, seg_files)):
        dis = load(dis_file)
        seg = load(seg_file)

        x_rav = np.ravel(dis['x'])
        y_rav = np.ravel(dis['y'])

        try:
            u_sum += np.ravel(dis['u'])
            v_sum += np.ravel(dis['v'])
        except:
            u_sum = np.ravel(dis['u'])
            v_sum = np.ravel(dis['v'])

        cx, cy = seg['centroid']

        distance, displacement, pressure = infer_pressure(x_rav, y_rav, u_sum, v_sum, cx, cy, r0, get_pressure)

        mask = distance > r_min

        pressure_mean = np.nanmean(pressure[mask])
        pressure_median = np.nanmedian(pressure[mask])
        pressure_std = np.nanstd(pressure[mask], ddof=1)

        contractility_mean = pressure_mean*A0*(muperpixel**2.)*(10**6)  # unit: µN
        contractility_median = pressure_median*A0*(muperpixel**2.)*(10**6)  # unit: µN
        contractility_std = pressure_std*A0*(muperpixel**2.)*(10**6)  # unit: µN

        results['pressure_mean'].append(pressure_mean)
        results['pressure_median'].append(pressure_median)
        results['pressure_std'].append(pressure_std)

        results['contractility_mean'].append(contractility_mean)
        results['contractility_median'].append(contractility_median)
        results['contractility_std'].append(contractility_std)

    df = pd.DataFrame.from_dict(results)
    df.columns = ['Mean Pressure (Pa)',
                  'Median Pressure (Pa)',
                  'St.dev. Pressure (Pa)',
                  'Mean Contractility (µN)',
                  'Median Contractility (µN)',
                  'St.dev. Contractility (µN)']

    if outfile is not None:
        df.to_excel(outfile)

    return df


def infer_pressure(x_rav, y_rav, u_rav, v_rav, x_sph, y_sph, r_sph, get_pressure):
    dx = (x_rav[~np.isnan(u_rav)] - x_sph) / r_sph
    dy = (y_rav[~np.isnan(u_rav)] - y_sph) / r_sph
    distance = np.sqrt(dx ** 2 + dy ** 2)

    u_rav2 = u_rav[~np.isnan(u_rav)] / r_sph
    v_rav2 = v_rav[~np.isnan(u_rav)] / r_sph

    m = (np.array([dx, dy]).T) / np.expand_dims(distance, axis=1)
    d = np.array([u_rav2, v_rav2]).T

    displacement = np.array([-np.dot(di, mi) for di, mi in zip(d, m)])

    abs = np.sqrt(u_rav2**2. + v_rav2**2.)
    mask = displacement/abs > 0.94 # cos(20deg)
    distance = distance[mask]
    displacement = displacement[mask]

    pressure = get_pressure(distance, displacement)

    return [distance, displacement, pressure]
