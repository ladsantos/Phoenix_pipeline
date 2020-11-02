import numpy as np
from astropy import visualization as aviz
from astropy.nddata.utils import block_reduce, Cutout2D
from matplotlib import pyplot as plt
import matplotlib.colors as clr
import os

def save_image(image, name = 'Image', path = os.getcwd()):
	"""
	image = image data
	name = name of the image, default is Image
	path = path of the resultant image,
		default is the the present working directory
	"""
	fig, ax = plt.subplots(1, 1, figsize = (5,10))
	im = ax.imshow(image, origin = 'lower', norm = clr.LogNorm(vmin=3, vmax=20000),cmap='PuBu_r')
	fig.colorbar(im, ax = ax, fraction = 0.046, pad = 0.04)
	fig.savefig(path, name + '.png')



def find_nearest_dark_exposure(image, dark_exposure_times, tolerance=0.5):
    """
    Find the nearest exposure time of a dark frame to the exposure time of the image,
    raising an error if the difference in exposure time is more than tolerance.

    This function is taken from the CCD Data reduction guide:
    https://mwcraig.github.io/ccd-as-book/00-00-Preface.html
    
    Parameters
    ----------
    
    image : astropy.nddata.CCDData
        Image for which a matching dark is needed.
    
    dark_exposure_times : list
        Exposure times for which there are darks.
    
    tolerance : float or ``None``, optional
        Maximum difference, in seconds, between the image and the closest dark. Set
        to ``None`` to skip the tolerance test.
    
    Returns
    -------
    
    float
        Closest dark exposure time to the image.
    """

    dark_exposures = np.array(list(dark_exposure_times))
    idx = np.argmin(np.abs(dark_exposures - image.header['exptime']))
    closest_dark_exposure = dark_exposures[idx]

    if (tolerance is not None and 
        np.abs(image.header['exptime'] - closest_dark_exposure) > tolerance):
        return None
        #raise RuntimeError('Closest dark exposure time is {} for flat of exposure '
                           #'time {}.'.format(closest_dark_exposure, image.header['exptime']))
        
    
    return closest_dark_exposure


def inverse_median(a):
	return 1/np.median(a)
