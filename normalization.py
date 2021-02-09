import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy.stats import mad_std
import utils as utl
import astropy.units as u
import shutil

def normal(flux_file, degree=1):
    """
    Parameters:
    -----------
    flux_file : .dat file
        flux file generated from flux_extraction.py file
    degree : int
        degree of the fitted polynomial
        default is 1
    -----------
    returns
    -----------
    normal_flux : .dat file
        data file containing pixel, normalized flux and error in normalized flux.
    """
    pix, fl, fle = np.loadtxt(flux_file, usecols=(0,1,2), unpack=True)
    