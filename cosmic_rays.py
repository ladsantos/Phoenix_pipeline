import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
import os
from pathlib import Path
from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy.stats import mad_std
import utils as utl
import astropy.units as u
import shutil

def cosmic_removal(loc, files, loc_out):
    """
    Parameters:
    -----------
    loc : str
       location of normalized flux files
    files : tuple
       tuple containing the names of normalized 
       spectrum files
    loc_out : str
        location of the output files
    -----------
    returns
    -----------
    files : .dat files
        data files containing the pixel number,
        normalized flux and error in it, with
        removal of cosmic rays.
    """
    # Creating Master flux
    master_pix, all_fl, all_fle = np.loadtxt(loc + files[0], usecols=(0,1,2), unpack=True)
    for i in range(len(files)-1):
        fl1, fle1 = np.loadtxt(loc + files[i+1], usecols=(1,2), unpack=True)
        all_fl = np.vstack((all_fl, fl1))
        all_fle = np.vstack((all_fle, fle1))
    
    master_fl = np.mean(all_fl, axis=0)
    master_fle = np.array([])
    for i in range(len(all_fle[0])):
        fel11 = np.sum(all_fle[:,i]**2)/len(all_fle[:,i])
        master_fle = np.hstack((master_fle, fel11))


    # master pixels - master_pix
    # master flux - master_fl
    # error in master flux - master_fle
    # Replacing 5-sigma outliers with average values of four pixels around it
    for i in range(len(files)):
        ff = open(loc_out + 'cosmic_' + files[i], 'w')
        pixx, fll, flle = np.loadtxt(loc + files[i], usecols=(0, 1, 2), unpack=True)
        for j in range(len(fll)):
            if j != len(fll)-2 or 0 or 1 or len(fll)-1:
                if np.abs(fll[j] - master_fl[j]) > 5*flle[j]:
                    fll[j] = (fll[j-1] + fll[j-2] + fll[j+1] + fll[j+2])/4
                    flle[j] = (flle[j-1]**2 + flle[j-2]**2 +
                            flle[j+1]**2 + flle[j+2]**2)/4
                else:
                    continue
            elif j == 0 or j == 1:
                if np.abs(fll[j] - master_fl[j]) > 5*flle[j]:
                    fll[j] = fll[2]
                    flle[j] = flle[2]
                else:
                    continue
            elif j == len(fll)-2 or j == len(fll)-1:
                if np.abs(fll[j] - master_fl[j]) > 5*flle[j]:
                    fll[j] = fll[j-2]
                    flle[j] = flle[j-2]
                else:
                    continue
        for k in range(len(pixx)):
            ff.write(str(pixx[k]) + '\t' + str(fll[k]) +
                    '\t' + str(flle[k]) + '\n')
        ff.close()