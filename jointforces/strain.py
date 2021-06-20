import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
from natsort import natsorted
from scipy.ndimage import generic_filter
from skimage.morphology import disk
from scipy.interpolate import LinearNDInterpolator
from saenopy import macro
from saenopy.materials import SemiAffineFiberMaterial, LinearMaterial
from scipy import interpolate
import matplotlib 
import os


def get_displ(x_rav, y_rav, u_rav, v_rav, x_sph, y_sph):
    dx = (x_rav - x_sph)
    dy = (y_rav - y_sph)

    distance = np.sqrt(dx ** 2 + dy ** 2)
    angle = np.arctan2(dy, dx)
    m = (np.array([dx, dy]).T) / np.expand_dims(distance, axis=1)
    d = np.array([u_rav, v_rav]).T
    displacement = np.array([-np.dot(di, mi) for di, mi in zip(d, m)])

    return [distance, angle, displacement]


def create_strain_maps(folder, delta, outfolder='strain-maps', radius=2, i_max=None, vmin=1, vmax=1.5, ticks_strains=[ 1, 1.1 , 1.2 , 1.3, 1.4, 1.5]):
    """
    Creates 2D Strain maps for a time series. 
    (and single plot with mean and maximal strain for eacht time step)
    
    folder: containing the defromation series (output folder of compute_displacement_series)
    outfolder: output of Strain maps 
    delta: Windowsize (px) as it was used in the evaluation of deformations for correct strain
    radius: used to smooth the deformation field    
    i_max: can be used to evaluate only first i timesteps
    """
    # creates folder if it doesn't exist
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    # creates folder if it doesn't exist
    if not os.path.exists(os.path.join(outfolder,"plots")):
        os.makedirs(os.path.join(outfolder,"plots"))   
        
    
    # load segmentations
    seg_files = natsorted(glob(folder+'/seg*.npy'))
    

    # load deformations
    # look for accumulated deformation files (new standard)
    d_accumulated_files = natsorted(glob(folder+'/def*.npy'))
    # look also for not-accummulated deformations (old standard)
    d_notaccumulated_files = natsorted(glob(folder+'/dis*.npy'))   
    # if not-accumulated deformations are found chose different mode
    if len(d_notaccumulated_files) > len(d_accumulated_files):
        accumulated = False
        dis_files = d_notaccumulated_files
        print("Found not-accumulated deformation files (old standard) and will conduct calculations accordingly.")
    # else do the calcualtion with accumulated d eformations already
    else:
        accumulated = True
        dis_files = d_accumulated_files
    

    dis = np.load(dis_files[0], allow_pickle="True").item()

    x_rav, y_rav = np.ravel(dis['x']), np.ravel(dis['y'])
    u_rav, v_rav = np.zeros_like(np.ravel(dis['u'])), np.zeros_like(np.ravel(dis['v']))

    if not i_max:
        i_max= len(dis_files)
    median_strain = []
    max_strain = []
    max_defo = []
    median_defo = []
    for i in tqdm(range(i_max)):
        # get displacement field
        dis = np.load(dis_files[i], allow_pickle="True").item()
        
        # get deformations
        # sum up if we have not-accummulated deformations (old standard)
        if accumulated == False:
            try:
                u_rav += np.ravel(dis['u'])
                v_rav += np.ravel(dis['v'])
            except:
                u_rav = np.ravel(dis['u'])
                v_rav = np.ravel(dis['v'])
        # else read in accummulated deformations directly (new standard)
        else:
            u_rav = np.ravel(dis['u'])
            v_rav = np.ravel(dis['v'])
        

        seg = np.load(seg_files[i], allow_pickle="True").item()
        x_sph, y_sph = seg['centroid']

        dist, angle, displ = get_displ(x_rav, y_rav, u_rav, v_rav, x_sph, y_sph)

        # make sure spheroid mask is visible at the end
        dist[np.isnan(u_rav)] = np.nan

        # smooth displacement field
        displ = displ.reshape(dis['x'].shape).copy()
        displ = generic_filter(displ, np.nanmean, footprint=disk(radius), mode='nearest')

        # linear interpolation of discrete displacement field
        displ = np.ravel(displ)
        x, y = np.ravel(dis['x']), np.ravel(dis['y'])

        mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(displ))

        x2 = x[mask]
        y2 = y[mask]
        displ2 = displ[mask]

        f = LinearNDInterpolator(np.array([x2, y2]).T, displ2)

        # compute strain
        x_old = (dist * np.cos(angle) + x_sph).reshape(dis['x'].shape)
        y_old = (dist * np.sin(angle) + y_sph).reshape(dis['y'].shape)
        dist_old = f(x_old, y_old)

        x_new = ((dist + delta) * np.cos(angle) + x_sph).reshape(dis['x'].shape)
        y_new = ((dist + delta) * np.sin(angle) + y_sph).reshape(dis['y'].shape)
        dist_new = f(x_new, y_new)

        strain = 1 + (dist_old - dist_new) / delta

        #print (strain)
        print ("")
        print ("Max. Strain: "+str(np.nanmax(strain)))
        max_defo.append(np.nanmax(np.sqrt(u_rav**2 + v_rav**2 )))
        median_defo.append(np.nanmedian(np.sqrt(u_rav**2 + v_rav**2 )))
        max_strain.append(np.nanmax(strain))
        median_strain.append(np.nanmedian(strain))
        plt.figure(figsize=(8, 6))
        height = strain.shape[0]
        width = strain.shape[1]
        # plot strain map            
        plt.imshow(strain, vmin=vmin, vmax=vmax, cmap='magma', extent=[0, width, height, 0], origin='lower')
        cb = plt.colorbar(fraction=0.025, label='Strain', ticks=ticks_strains)
            
        #cb.ax.set_yticklabels(['$\mathregular{\leq}$ 300 Pa', '400 Pa', '500 Pa', '$\mathregular{\geq}$ 600 Pa'])
        plt.xlim((0, width))
        plt.ylim((0, height))
        plt.axis('off')
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.tight_layout()
        plt.savefig(outfolder+'/strain'+str(i+1).zfill(6)+'.png', dpi=250)
        plt.close()
   
    # save last strain map + deformations interpolated
    np.save(os.path.join(outfolder,"strain_end.npy"),strain)
    np.save(os.path.join(outfolder,"x_end.npy"),x2)
    np.save(os.path.join(outfolder,"y_end.npy"),y2)
    np.save(os.path.join(outfolder,"displ_end.npy"),displ2)
    # raw accumulated data
    np.save(os.path.join(outfolder,"x_rav_end.npy"),x_rav)
    np.save(os.path.join(outfolder,"y_rav_end.npy"),y_rav)
    np.save(os.path.join(outfolder,"u_rav_end.npy"),u_rav)
    np.save(os.path.join(outfolder,"v_rav_end.npy"),v_rav)  
    
    
    # plot and save overview
    np.savetxt(os.path.join(outfolder,"plots")+'/median_strain.txt', median_strain)
    np.savetxt(os.path.join(outfolder,"plots")+'/max_strain.txt', max_strain)  
    np.savetxt(os.path.join(outfolder,"plots")+'/max_defo.txt', max_defo)  
    np.savetxt(os.path.join(outfolder,"plots")+'/median_defo.txt', median_defo)  
    
    # deformations  
    plt.figure(figsize=(6,2))
    plt.grid()
    plt.xlabel("timesteps")
    plt.ylabel("deformations (px)")
    plt.plot(max_defo)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,"plots")+'/max_defo.png', dpi=300)
    plt.close()
    # median strain    
    plt.figure(figsize=(6,2))
    plt.grid()
    plt.xlabel("timesteps", fontsize=7)
    plt.ylabel("median deformations (px)", fontsize=7)
    plt.plot(median_defo)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,"plots")+'/median_deformation.png', dpi=300)
    plt.close()
    # median strain    
    plt.figure(figsize=(6,2))
    plt.grid()
    plt.xlabel("timesteps")
    plt.ylabel("median strain")
    plt.plot(median_strain)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,"plots")+'/median_strain.png', dpi=300)
    plt.close()
    # max strain
    # max strain
    plt.figure(figsize=(6,2))
    plt.grid()
    plt.xlabel("timesteps")
    plt.ylabel("max strain")
    plt.plot(max_strain)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,"plots")+'/max_strain.png', dpi=300)
    plt.close()
    return 


def create_stiffness_maps(folder, delta, outfolder='stiffness-maps', k0=1645, eps_s=0.0075, d_s=0.033, d_0=0.0008,
                           radius=2, i_max=None):
    
    """
    Creates 2D stiffness maps for a time series. 
    (and single plot with mean and maximal stiffness for eacht time step)
    
    folder: containing the defromation series (output folder of compute_displacement_series)
    outfolder: output of Strain maps 
    delta: Windowsize (px) as it was used in the evaluation of deformations for correct strain
    radius: used to smooth the deformation field    
    i_max: can be used to evaluate only first i timesteps
    k0, eps_s,d_s,_do: material parameters to compute the stiffness

    """
    # creates folder if it doesn't exist
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    # creates folder if it doesn't exist
    if not os.path.exists(os.path.join(outfolder,"plots")):
        os.makedirs(os.path.join(outfolder,"plots"))   
        
    # fiber material
    material = SemiAffineFiberMaterial(k0, d_0, eps_s, d_s)
    # loead extensional stress data and get stiffness - might need to adjust range
    epsilon = np.arange(-1, 1.9, 0.0001)
    x, y = macro.getExtensionalRheometerStress(epsilon, material)
    get_stiffness = interpolate.interp1d(x[1:], (y[1:]-y[:-1])/(x[1:]-x[:-1]))
        
    # load in displacement data and ssheroid mask
    # load segmentations
    seg_files = natsorted(glob(folder+'/seg*.npy'))

    # load deformations
    # look for accumulated deformation files (new standard)
    d_accumulated_files = natsorted(glob(folder+'/def*.npy'))
    # look also for not-accummulated deformations (old standard)
    d_notaccumulated_files = natsorted(glob(folder+'/dis*.npy'))   
    # if not-accumulated deformations are found chose different mode
    if len(d_notaccumulated_files) > len(d_accumulated_files):
        accumulated = False
        dis_files = d_notaccumulated_files
        print("Found not-accumulated deformation files (old standard) and will conduct calculations accordingly.")
    # else do the calcualtion with accumulated d eformations already
    else:
        accumulated = True
        dis_files = d_accumulated_files
    
    

    dis = np.load(dis_files[0], allow_pickle="True").item()

    x_rav, y_rav = np.ravel(dis['x']), np.ravel(dis['y'])
    u_rav, v_rav = np.zeros_like(np.ravel(dis['u'])), np.zeros_like(np.ravel(dis['v']))
    
    # take all timesteps if not specified
    if not i_max:
        i_max= len(dis_files)
    
    # loop over all times
    median_stiffness = []
    max_stiffness = []    
    for i in tqdm(range(i_max)):
        # get displacement field
        dis = np.load(dis_files[i], allow_pickle="True").item()
      
        # get deformations
        # sum up if we have not-accummulated deformations (old standard)
        if accumulated == False:
            try:
                u_rav += np.ravel(dis['u'])
                v_rav += np.ravel(dis['v'])
            except:
                u_rav = np.ravel(dis['u'])
                v_rav = np.ravel(dis['v'])
        # else read in accummulated deformations directly (new standard)
        else:
            u_rav = np.ravel(dis['u'])
            v_rav = np.ravel(dis['v'])
   

        seg = np.load(seg_files[i], allow_pickle="True").item()
        x_sph, y_sph = seg['centroid']

        dist, angle, displ = get_displ(x_rav, y_rav, u_rav, v_rav, x_sph, y_sph)

        # make sure spheroid mask is visible at the end
        dist[np.isnan(u_rav)] = np.nan

        # smooth displacement field
        displ = displ.reshape(dis['x'].shape).copy()
        displ = generic_filter(displ, np.nanmean, footprint=disk(radius), mode='nearest')

        # linear interpolation of discrete displacement field
        displ = np.ravel(displ)
        x, y = np.ravel(dis['x']), np.ravel(dis['y'])

        mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(displ))

        x2 = x[mask]
        y2 = y[mask]
        displ2 = displ[mask]

        f = LinearNDInterpolator(np.array([x2, y2]).T, displ2)

        # compute strain
        x_old = (dist * np.cos(angle) + x_sph).reshape(dis['x'].shape)
        y_old = (dist * np.sin(angle) + y_sph).reshape(dis['y'].shape)
        dist_old = f(x_old, y_old)

        x_new = ((dist + delta) * np.cos(angle) + x_sph).reshape(dis['x'].shape)
        y_new = ((dist + delta) * np.sin(angle) + y_sph).reshape(dis['y'].shape)
        dist_new = f(x_new, y_new)

        strain = 1 + (dist_old - dist_new) / delta

        # compute stiffening
        stiffness = get_stiffness(strain)
        log_stiffness = np.log10(stiffness)
        #stiffness = np.zeros_like(strain) + np.nan  # initalize with linear stiffness
        #stiffness[strain <= 0] = k0 * np.exp(strain[strain <= 0] / d_0)  # buckling
        #stiffness[(strain > 0) & (strain < eps_s)] = k0  # linear regime
        #stiffness[strain >= eps_s] = k0 * np.exp((strain[strain >= eps_s] - eps_s) / d_s)  # strain stiffening
        # plot result
        plt.figure(figsize=(8, 6))
        height = log_stiffness.shape[0]
        width = log_stiffness.shape[1]
        plt.imshow(log_stiffness, vmin=2, vmax=4, cmap='inferno', extent=[0, width, height, 0], origin='lower')
        cb = plt.colorbar(fraction=0.025, label='Matrix stiffness', ticks=np.log10([100, 250, 500, 1000, 2500, 5000, 10000]))
        cb.ax.set_yticklabels(['$\mathregular{\leq}$ 100 Pa', '250 Pa', '500 Pa', '1.0 kPa', '2.5 kPa', '5.0 kPa', '$\mathregular{\geq}$ 10.0 kPa'])
        plt.xlim((0, width))
        plt.ylim((0, height))
        plt.axis('off')
        plt.gca().get_xaxis().set_visible(False)
        plt.gca().get_yaxis().set_visible(False)
        plt.tight_layout()
        plt.savefig(outfolder+'/stiffening'+str(i+1).zfill(6)+'.png', dpi=250)
        print('')
        print("Max. Stiffness: "+str(np.nanmax(stiffness)))
        print('')
        median_stiffness.append(np.nanmedian(stiffness))
        max_stiffness.append(np.nanmax(stiffness))
        plt.close()
        
     # plot and save overview    
    np.savetxt(os.path.join(outfolder,"plots")+'/median_stiffness.txt', median_stiffness)
    np.savetxt(os.path.join(outfolder,"plots")+'/max_stiffness.txt', max_stiffness)
    # median stiffness
    plt.figure(figsize=(6,2))
    plt.grid()
    plt.xlabel("timesteps")
    plt.ylabel("median stiffness")
    plt.plot(median_stiffness)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,"plots")+'/median_stiffness.png', dpi=300)
    plt.close()
    # max stiffness
    plt.figure(figsize=(6,2))
    plt.grid()
    plt.xlabel("timesteps")
    plt.ylabel("max. stiffness")
    plt.plot(max_stiffness)
    plt.tight_layout()
    plt.savefig(os.path.join(outfolder,"plots")+'/max_stiffness.png', dpi=300)
    plt.close()
    return 
