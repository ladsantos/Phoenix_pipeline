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

def normal_spectrum(flux_file, out_path, degree=3):
    """
    Parameters:
    -----------
    flux_file : .dat file
        flux file generated from flux_extraction.py file
    out_path : str
        Path of the output data file
    degree : int
        degree of the fitted polynomial
        default is 3
    -----------
    returns
    -----------
    normal_flux : .dat file
        data file containing pixel, normalized flux and error in normalized flux.
    """
    pix, fl, fle = np.loadtxt(flux_file, usecols=(0,1,2), unpack=True)
    sigmas = np.copy(fle)
    # ---------------------------------------
    # To select the region of spectral lines
    # ---------------------------------------
    selection = []
    fig = plt.figure(figsize=(16, 12))
    ax = fig.add_subplot(211)

    ax.errorbar(pix, fl, yerr=fle, '-')
    #ax.set_ylim(-2, 2)
    ax.set_title('Press left mouse button and drag to select the region of the spectrum you want to mask.')

    ax2 = fig.add_subplot(212)
    line2, = ax2.plot(pix, fl, '-')

    def onselect(xmin, xmax):
        indmin, indmax = np.searchsorted(x, (xmin, xmax))
        indmax = min(len(x) - 1, indmax)

        thisx = x[indmin:indmax]
        thisy = y[indmin:indmax]
        line2.set_data(thisx, thisy)
        ax2.set_xlim(thisx[0], thisx[-1])
        ax2.set_ylim(thisy.min(), thisy.max())
        fig.canvas.draw_idle()
        ab = [indmin, indmax]
        selection.append(ab)

    # set useblit True on gtkagg for enhanced performance
    span = SpanSelector(ax, onselect, 'horizontal', useblit=True, 
                         rectprops=dict(alpha=0.5, facecolor='red'))

    plt.show()
    print(selection)

    for i in range(len(selection)):
        sigmas[int(selection[0]):int(selection[1])] = 0
    
    bb = np.polyfit(pix, fl, sigmas, deg=degree)
    # Remeber that polyfit gives you an array of the
    # coefficient of the polynomial, starting with the highest degree.
    # while the function in the utilities uses the array
    # starting with a coefficients of lowest degree.

    coefs = np.flip(bb)
    polynom = utl.arbi_poly(pix, coefs)
    normal_flux = fl/polynom
    normal_flux_err = fle/polynom
    flux_file1 = 'normal_' + flux_file
    f11 = open(out_path, flux_file1, 'w')
    for i in range(len(pix)):
        f11.write(str(pix[i]) + '\t' + str(normal_flux[i]) + '\t' + str(normal_flux_err[i]) + '\n')
    f11.close()

pt1 = os.getcwd()
ff1 = 'sky_sub_2009oct30_0009.fits_flux.dat'
normal_spectrum(ff1, pt1, degree=3)