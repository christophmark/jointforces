import warnings
import os
import gc
import numpy as np
import scipy.ndimage.morphology as scipy_morph
import scipy.ndimage.measurements as scipy_meas
from skimage.filters import gaussian, threshold_otsu
from skimage.morphology import remove_small_objects
from glob import glob
from tqdm import tqdm
from natsort import natsorted
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import openpiv.pyprocess
import openpiv.process
import openpiv.filters


def enhance_contrast(img):
    img = img.astype(np.float)
    img = gaussian(img, 2)
    img -= np.percentile(img, 0.5)
    img /= np.percentile(img, 99.5)
    img[img < 0.] = 0.
    img[img > 1.] = 1.
    return img


def segment_spheroid(img, enhance=True):
    """
    Image segmentation function to create spheroid mask, radius, and position of a spheroid in a grayscale image.

    Args:
        img(array): Grayscale image as a Numpy array

    Returns:
        dict: Dictionary with keys: mask, radius, centroid (x/y)
    """
    height = img.shape[0]
    width = img.shape[1]

    # contrast enhancement
    if enhance:
        img = enhance_contrast(img)

    # flip y (to match x/y coordinates) and binarize
    mask = img[::-1] < threshold_otsu(img) * 0.9

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


def compute_displacements(window_size, img0, img1, mask1=None, cutoff=None, drift_correction=True):
    # get image size
    height = img1.shape[0]
    width = img1.shape[1]

    # ignore ubiquitous warning messages from OpenPIV
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # compute displacements
        ut, vt, sig2noise = openpiv.pyprocess.piv(img0,
                                                  img1,
                                                  window_size=window_size,
                                                  overlap=window_size // 2,
                                                  dt=1,
                                                  sig2noise_method='peak2peak')

        # replace outliers
        ut, vt = openpiv.filters.replace_outliers(ut, vt, method='localmean', max_iter=3, kernel_size=1)

        # get coordinates corresponding to displacement
        x, y = openpiv.process.get_coordinates(image_size=img1.shape,
                                               window_size=window_size,
                                               overlap=window_size // 2)

    # if mask is specified, replace displacements within mask with NaN
    # if cutoff is specified, replace displacements further out from the center than cutoff value with NaN
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

    return {'x': x, 'y': y, 'u': ut, 'v': vt}


def save_displacement_plot(filename, img, segmentation, displacements, quiver_scale=1, color_norm=75.):
    # get image size
    height = img.shape[0]
    width = img.shape[1]

    x, y, u, v = displacements['x'], displacements['y'], displacements['u'], displacements['v']
    mask = segmentation['mask']
    cx, cy = segmentation['centroid']

    fig = plt.figure(figsize=(10 * (width / height), 10))
    plt.imshow(img, cmap='Greys_r', extent=[0, width, height, 0], origin='lower')

    d = (u ** 2 + v ** 2) ** 0.5

    p = plt.quiver(x, y, u, v, d,
                   clim=[0, color_norm],
                   cmap=cm.jet,
                   alpha=0.8,
                   scale=quiver_scale,
                   units='xy',
                   pivot='mid')

    overlay = mask.astype(int) - scipy_morph.binary_erosion(mask, iterations=4).astype(int)
    overlay = np.array([overlay.T,
                        np.zeros_like(overlay).T,
                        np.zeros_like(overlay).T,
                        overlay.T]).T

    plt.imshow(overlay.astype(np.float), extent=[0, width, height, 0], zorder=1000)
    plt.scatter([cx], [cy], s=50, lw=1, edgecolors='k', c='r')

    plt.xlim((0, width))
    plt.ylim((0, height))

    plt.axis('off')
    plt.gca().get_xaxis().set_visible(False)
    plt.gca().get_yaxis().set_visible(False)

    plt.savefig(filename, bbox_inches='tight', pad_inches=0, dpi=150)
    plt.clf()
    plt.close(fig)


def compute_displacement_series(folder, filter, outfolder, n_max=None, enhance=True,
                                window_size=70, cutoff=None, drift_correction=True,
                                plot=True, quiver_scale=1, color_norm=75.):
    img_files = natsorted(glob(folder+'/'+filter))

    if n_max is not None:
        img_files = img_files[:n_max+1]

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    img0 = plt.imread(img_files[0])
    seg0 = segment_spheroid(img0, enhance=enhance)

    np.save(outfolder+'/seg000000.npy', seg0)

    u_sum = None
    v_sum = None

    for i in tqdm(range(1, len(img_files))):
        img1 = plt.imread(img_files[i])
        seg1 = segment_spheroid(img1, enhance=enhance)

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

            save_displacement_plot(outfolder+'/plot'+str(i).zfill(6)+'.png', img1, seg1, dis_sum,
                                   quiver_scale=quiver_scale, color_norm=color_norm)

        img0 = img1.copy()
