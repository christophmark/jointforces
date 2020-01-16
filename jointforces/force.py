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
               'contractility_mean': [], 'contractility_median': [], 'contractility_std': [],
               'pressure_min': [], 'pressure_max': [], 'contractility_min': [], 'contractility_max': [],
               'angle_min': [], 'angle_max': []}

    # create dict with all angles
    angles_dict = {str(a): [] for a in range(-175, 175, 5)}
    
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

        distance, displacement, angle, pressure = infer_pressure(x_rav, y_rav, u_sum, v_sum, cx, cy, r0, get_pressure)
        mask = distance > r_min

        pr_angle = []
        pr_median = []
        for alpha in range(-175, 175, 5):
            mask2 = (angle >= (alpha-5)*np.pi/180.) & (angle < (alpha+5)*np.pi/180.)
            pr_angle.append(alpha)
            pr_median.append(np.nanmedian(pressure[mask & mask2]))

        i_min = np.nanargmin(pr_median)
        alpha_min = pr_angle[i_min]
        pressure_min = pr_median[i_min]

        i_max = np.nanargmax(pr_median)
        alpha_max = pr_angle[i_max]
        pressure_max = pr_median[i_max]

        pressure_mean = np.nanmean(pressure[mask])
        pressure_median = np.nanmedian(pressure[mask])
        pressure_std = np.nanstd(pressure[mask], ddof=1)

        contractility_mean = pressure_mean*A0*(muperpixel**2.)*(10**6)  # unit: µN
        contractility_median = pressure_median*A0*(muperpixel**2.)*(10**6)  # unit: µN
        contractility_std = pressure_std*A0*(muperpixel**2.)*(10**6)  # unit: µN

        contractility_min = pressure_min*A0*(muperpixel**2.)*(10**6)  # unit: µN
        contractility_max = pressure_max*A0*(muperpixel**2.)*(10**6)  # unit: µN

        results['pressure_mean'].append(pressure_mean)
        results['pressure_median'].append(pressure_median)
        results['pressure_std'].append(pressure_std)

        results['contractility_mean'].append(contractility_mean)
        results['contractility_median'].append(contractility_median)
        results['contractility_std'].append(contractility_std)

        results['pressure_min'].append(pressure_min)
        results['pressure_max'].append(pressure_max)

        results['contractility_min'].append(contractility_min)
        results['contractility_max'].append(contractility_max)

        results['angle_min'].append(alpha_min)
        results['angle_max'].append(alpha_max)
        
        # append pressures for all angle data
        for i,a in enumerate(angles_dict):
            angles_dict[a].append(pr_median[i])

    df = pd.DataFrame.from_dict(results)
    df.columns = ['Mean Pressure (Pa)',
                  'Median Pressure (Pa)',
                  'St.dev. Pressure (Pa)',
                  'Mean Contractility (µN)',
                  'Median Contractility (µN)',
                  'St.dev. Contractility (µN)',
                  'Minimal Median Pressure (Pa)',
                  'Maximal Median Pressure (Pa)',
                  'Minimal Median Contractility (µN)',
                  'Maximal Median Contractility (µN)',
                  'Angle of minimal Pr./Contr. (deg)',
                  'Angle of maximal Pr./Contr. (deg)']

    if outfile is not None:
        df.to_excel(outfile)
    else:
        df.to_excel(folder+'//result.xlsx')
     
    # save pressures for all angles
    an = pd.DataFrame.from_dict(angles_dict)
    an.columns = [a for a in range(-175, 175, 5)]
    if outfile is not None:
        an.to_excel(outfile[:-5]+'_angles.xlsx')
    else:
        an.to_excel(folder+'//result_angles.xlsx')    

    return df    
 



def infer_pressure(x_rav, y_rav, u_rav, v_rav, x_sph, y_sph, r_sph, get_pressure):
    dx = (x_rav[~np.isnan(u_rav)] - x_sph) / r_sph
    dy = (y_rav[~np.isnan(u_rav)] - y_sph) / r_sph
    distance = np.sqrt(dx ** 2 + dy ** 2)
    angle = np.arctan2(dy, dx)

    u_rav2 = u_rav[~np.isnan(u_rav)] / r_sph
    v_rav2 = v_rav[~np.isnan(u_rav)] / r_sph

    m = (np.array([dx, dy]).T) / np.expand_dims(distance, axis=1)
    d = np.array([u_rav2, v_rav2]).T

    displacement = np.array([-np.dot(di, mi) for di, mi in zip(d, m)])

    abs = np.sqrt(u_rav2**2. + v_rav2**2.)
    mask = displacement/abs > 0.94  # cos(20deg)

    distance = distance[mask]
    angle = angle[mask]
    displacement = displacement[mask]

    pressure = get_pressure(distance, displacement)

    return [distance, displacement, angle, pressure]
