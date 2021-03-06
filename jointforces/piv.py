import warnings
import os
import gc
import numpy as np
import scipy.ndimage.morphology as scipy_morph
import scipy.ndimage.measurements as scipy_meas
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import remove_small_objects
from skimage.exposure import adjust_gamma   
from skimage import io
from glob import glob
from tqdm import tqdm
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import openpiv.pyprocess
#import openpiv.process
import openpiv.filters
from roipoly import RoiPoly



def enhance_contrast(img, gauss=False, gamma=None):
    img = img.astype(np.float)
    # apply gaussian filter 
    if gauss==True:
        img = gaussian(img, 2)
    # contrast enhancement
    img -= np.percentile(img, 0.5)
    img /= np.percentile(img, 99.5)
    img[img < 0.] = 0.
    img[img > 1.] = 1.
    img /= 256
    # gamma correction 
    if gamma is not None:
        img = adjust_gamma(img, gamma)  
    return img
        


def segment_spheroid(img, enhance=True, thres = 0.9):
    """
    Image segmentation function to create spheroid mask, radius, and position of a spheroid in a grayscale image.

    Args:
        img(array): Grayscale image as a Numpy array
        thres(float): To adjust the segmentation

    Returns:
        dict: Dictionary with keys: mask, radius, centroid (x/y)
    """
    height = img.shape[0]
    width = img.shape[1]
    
    # contrast enhancement
    if enhance:
        img = enhance_contrast(img, gauss=True)
        
    
    # flip y (to match x/y coordinates) and binarize
    mask = img[::-1] < threshold_otsu(img) * thres    
    
    # remove other objects
    mask = scipy_morph.binary_dilation(mask, iterations=3)
    mask = scipy_morph.binary_fill_holes(mask)
    mask = remove_small_objects(mask, min_size=1000)
    mask = scipy_morph.binary_closing(mask, iterations=3)
    mask = scipy_morph.binary_fill_holes(mask)

    # identify spheroid as the most centered object
    labeled_mask, max_lbl = scipy_meas.label(mask)
    center_of_mass = np.array(scipy_meas.center_of_mass(mask, labeled_mask, range(1, max_lbl + 1)))
    distance_to_center = np.sqrt(np.sum((center_of_mass - np.array([height / 2, width / 2])) ** 2, axis=1))

    mask = (labeled_mask == distance_to_center.argmin() + 1)

    # determine radius of spheroid
    radius = np.sqrt(np.sum(mask) / np.pi)

    # determine center position of spheroid
    cy, cx = center_of_mass[distance_to_center.argmin()]

    # return dictionary containing spheroid information
    return {'mask': mask, 'radius': radius, 'centroid': (cx, cy)}



def custom_mask(img):
    """
    Image segmentation function to create a custom polygon mask, and evalute radius and position of the masked object.
    Need to use %matplotlib qt in jupyter notebook
    Args:
        img(array): Grayscale image as a Numpy array
    Returns:
        dict: Dictionary with keys: mask, radius, centroid (x/y)
    """
    height = img.shape[0]
    width  = img.shape[1]
    # click polygon mask interactive
    plt.ion()
    plt.imshow(img, extent=[0, width, height, 0])
    plt.text(0.5, 1.05,'Click Polygon Mask with left click, finish with right click',  fontsize=12,
         horizontalalignment='center',
         verticalalignment='center',c='darkred', transform= plt.gca().transAxes)#     transform = ax.transAxes)
    my_roi = RoiPoly(color='r')
    # Extract mask and segementation details
    mask = np.flipud(my_roi.get_mask(img))   # flip mask due to imshow 
    # determine radius of spheroid
    radius = np.sqrt(np.sum(mask) / np.pi)
    # determine center of mass
    cy, cx = scipy_meas.center_of_mass(mask)
    # hide pop-up windows   
    plt.ioff()
    # return dictionary containing spheroid information
    return {'mask': mask, 'radius': radius, 'centroid': (cx, cy)} 


def compute_displacements(window_size, img0, img1, mask1=None, cutoff=None, drift_correction=True):
    # get image size
    height = img1.shape[0]
    width = img1.shape[1]

    # ignore ubiquitous warning messages from OpenPIV
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # compute displacements
        ut, vt, sig2noise = openpiv.pyprocess.extended_search_area_piv(img0,   
                                                  img1,
                                                  window_size=window_size,
                                                  search_area_size = window_size,
                                                  overlap=window_size // 2,
                                                  dt=1,
                                                  sig2noise_method='peak2peak')

        # ToDo  Signal to noise filter option by threshold
        # u, v, mask = validation.sig2noise_val( ut, vt, sig2noise, threshold = 2.5 )
        # u, v = filters.replace_outliers( u, v, method='localmean', max_iter=10, kernel_size=2)

        # replace outliers
        ut, vt = openpiv.filters.replace_outliers(ut, vt, method='localmean', max_iter=3, kernel_size=1)

        # get coordinates corresponding to displacement
        x, y = openpiv.pyprocess.get_coordinates(image_size=img1.shape,
                                               search_area_size=window_size,
                                               overlap=window_size // 2)
    # flip y-values of PIV-arrows for correct orientation
    y = y[::-1]
    # if mask is specified, replace displacements within mask with NaN
    # if cutoff is specified, replace displacements further out from the center than cutoff value (px) with NaN
    for i, xi in enumerate(x[0, :]):
        for j, yj in enumerate(y[:, 0]):
            if mask1 is not None:
                if np.any(mask1[int(yj - 0.5 * window_size):int(yj + 0.5 * window_size),
                          int(xi - 0.5 * window_size):int(xi + 0.5 * window_size)]):
                    ut[j, i] = np.nan
                    vt[j, i] = np.nan
            if cutoff is not None:   
                if np.sqrt(((xi - 0.5 * width) ** 2. + (yj - 0.5 * height) ** 2.)) > cutoff:
                    ut[j, i] = np.nan
                    vt[j, i] = np.nan

    # subtracting mean displacements acts like a drift correction
    if drift_correction:
        ut -= np.nanmean(ut)
        vt -= np.nanmean(vt)
    return {'x': x, 'y': y , 'u': ut, 'v': vt}


def displacement_plot(img, segmentation, displacements, quiver_scale=1, color_norm=75., cmap=cm.jet, s=50, **kwargs):
    # get image size
    height = img.shape[0]
    width = img.shape[1]
    
    x, y, u, v = displacements['x'], displacements['y'], displacements['u'], displacements['v']
    mask = segmentation['mask']
    cx, cy = segmentation['centroid']

    plt.imshow(enhance_contrast(img), cmap='Greys_r', extent=[0, width, height, 0], origin='lower')
   
    if cmap is None:
        p = plt.quiver(x, y, u, v,
                       alpha=0.8,
                       scale=quiver_scale,
                       units='xy',
                       pivot='mid',
                       **kwargs)
    else:
        d = (u ** 2 + v ** 2) ** 0.5
        p = plt.quiver(x, y, u, v, d,
                       clim=[0, color_norm],
                       cmap=cmap,
                       alpha=0.8,
                       scale=quiver_scale,
                       units='xy',
                       pivot='mid',
                       **kwargs)

    overlay = mask.astype(int) - scipy_morph.binary_erosion(mask, iterations=4).astype(int)
    overlay = np.array([overlay.T,
                        np.zeros_like(overlay).T,
                        np.zeros_like(overlay).T,
                        overlay.T]).T
 
    plt.imshow(overlay.astype(np.float), extent=[0, width, height, 0], zorder=1000)
    plt.scatter([cx], [cy], lw=1, edgecolors='k', c='r', s=s)

    plt.xlim((0, width))
    plt.ylim((0, height))

    plt.axis('off')
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    
    return p


def save_displacement_plot(filename, img, segmentation, displacements, quiver_scale=1, color_norm=75., cmap=cm.jet,
                           s=50, **kwargs):
    height = img.shape[0]
    width = img.shape[1]
    fig = plt.figure(figsize=(10 * (width / height), 10))

    displacement_plot(img, segmentation, displacements, quiver_scale, color_norm, cmap, **kwargs)

    plt.savefig(filename, bbox_inches='tight', pad_inches=0, dpi=150)
    plt.clf()
    plt.close(fig)





def compute_displacement_series(folder, filter, outfolder, n_max=None, n_min=None,
                                enhance=True, window_size=70, cutoff=None, drift_correction=True,
                                plot=True, continous_segmentation = False, quiver_scale=1, color_norm=75., 
                                draw_mask = False, gamma=None, gauss=False, load_mask=None, thres_segmentation = 0.9,
                                cut_img = False, cut_img_val = (None,None,None,None)):
  
    img_files = natsorted(glob(folder+'/'+filter))
    # mainimal and maximal timesteps for analysis
    if n_max is not None:
        img_files = img_files[:n_max+1]
        
    if n_min is not None:
        img_files = img_files[n_min:]   
    # create output folder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    # read in image - (mask certain area in image if specified)    
    print(plt.imread(img_files[0]).dtype)

    if cut_img:
        img0 = io.imread(img_files[0], as_gray='True')[cut_img_val[0]:cut_img_val[1],  cut_img_val[2]:cut_img_val[3] ] 
    else:
        img0 = io.imread(img_files[0], as_gray='True')
        
    
      
    # segment , draw or load in mask
    if draw_mask == False:  # normal segmentation
        seg0 = segment_spheroid(img0, enhance=enhance, thres=thres_segmentation)
    else:                   # draw manual
        seg0 = custom_mask(img0)
    if load_mask is not None:   # load in mask
        seg0 = np.load(load_mask, allow_pickle=True).item()
        
    # save initial segmentation
    np.save(outfolder+'/seg000000.npy', seg0)

    u_sum = None
    v_sum = None
    
    # loop through all images
    for i in tqdm(range(1, len(img_files))):
        # read in new image (and cut if specified)
        if cut_img:
            img1= io.imread(img_files[i], as_gray='True')[cut_img_val[0]:cut_img_val[1],  cut_img_val[2]:cut_img_val[3] ] 
        else:
            img1 = io.imread(img_files[i], as_gray='True')
        
        # take the mask of the current timestep if continous_segmentation is active -
        # if not always take the mask from t0 
        # both options have advantages and disadvantages: continous_segmentation might
        # determine the center over time more accurately and offers information over growth,
        # however wrongly detected masks here may lead to force fluctuations 
        # and can reduce the FoV too much. Since the standard approach for the 
        # later force reconstruction is unsing the initial timestep anyway, 
        # the approach of always using the t0 mask here as well is quite robust 
        if (draw_mask == False) & (continous_segmentation == True):
            seg1 = segment_spheroid(img1, enhance=enhance, thres=thres_segmentation)
        else:
            seg1 = seg0.copy()
        
        # compute and save the matrx deformations and mask    
        dis = compute_displacements(window_size, img0, img1, mask1=seg1['mask'],
                                cutoff=cutoff, drift_correction=drift_correction)
        np.save(outfolder + '/seg'+str(i).zfill(6)+'.npy', seg1)
        np.save(outfolder + '/dis'+str(i).zfill(6)+'.npy', dis)
                     
        if plot:
            if u_sum is None:
                u_sum = dis['u']
                v_sum = dis['v']
            else:
                u_sum += dis['u']
                v_sum += dis['v']

            dis_sum = {'x': dis['x'], 'y': dis['y'], 'u': u_sum, 'v': v_sum}

            # plot same mask for all if continous segmentation is False or a custom mask is drawn
            if (draw_mask == False) & (continous_segmentation == True):
                save_displacement_plot(outfolder+'/plot'+str(i).zfill(6)+'.png', img1, seg1, dis_sum,
                                       quiver_scale=quiver_scale, color_norm=color_norm)
            # plot individual mask for each timestep ( however for forcereconstruction later we do use
            # the first timestep only to avoid errornous force fluctuations)    
            else:
                save_displacement_plot(outfolder+'/plot'+str(i).zfill(6)+'.png', img1, seg0, dis_sum,
                                       quiver_scale=quiver_scale, color_norm=color_norm)
        # next round we update the images        
        img0 = img1.copy()



def compute_noise_level(folder, filter, n_max=10, enhance=True,
                        window_size=70, cutoff=None, drift_correction=True):

    img_files = natsorted(glob(folder + '/' + filter))

    if n_max is not None:
        img_files = img_files[:n_max+1]

    img0 = plt.imread(img_files[0])
    mask = np.zeros_like(img0).astype(np.bool)

    displ = []

    for i in tqdm(range(1, len(img_files))):
        img1 = plt.imread(img_files[i])
        seg1 = segment_spheroid(img1, enhance=enhance)

        dis = compute_displacements(window_size, img0, img1, mask1=mask,
                                    cutoff=cutoff, drift_correction=drift_correction)

        displ.extend(list(np.ravel(np.sqrt(dis['u']**2. + dis['v']**2.))))

        img0 = img1.copy()

    return np.percentile(displ, 95)
