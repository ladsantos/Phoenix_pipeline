import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize as mz
from scipy.optimize import curve_fit as cft
from matplotlib.widgets import TextBox
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Button
import os
import utils as utl

def wave_soln(path, fname):
    """
    Parameters:
    -----------
    path : str
        path of the input file
    fname : str
        name of the input telluric file.
        Ideally this should be the file from which cosmic rays 
        are removed. 3 columns, pixel, flux, error in flux
    -----------
    returns
    -----------
    popt : array like
        optimized values of parameters of linear fit between pixel 
        space and wavelength space.
    pcov : 2x2 matrix
        covariance matrix of optimized values
    -----------
    """
    fig, ax = plt.subplots(figsize=(16,12))
    fig.subplots_adjust(bottom=0.2)

    pix, fl, fle = np.loadtxt(path + fname, usecols=(0,1,2), unpack=True)
    plts = ax.plot(pix, fl)
    ax.set_title('Press the left mouse button to select the region where the line is.')
    
    new_data = []    # Wavelength data provided by user
    indices = []

    # Dummy variables to save data
    new_data1 = []
    indices1 = []

    def onselect(pmin, pmax):
        global indices1
        indmin, indmax = np.searchsorted(pix, (pmin, pmax))
        indmax = min(len(pix) - 1, indmax)

        ab = [indmin, indmax]
        indices1 = []
        indices1.append(ab)

    def submit(expression):
        global new_data1
        ydata = expression
        new_data1 = []
        new_data1.append(ydata)

    def enter(self):
        global new_data1, indices1
        new_data.append(new_data1[0])
        indices.append(indices1[0])
        text_box1.set_val("")

    # To select a region
    span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                        rectprops=dict(alpha=0.5, facecolor='red'))

    # Adding a box to enter values
    axbox1 = fig.add_axes([0.2, 0.05, 0.4, 0.075])
    text_box1 = TextBox(axbox1, "Enter the corresponding\n wavelength here (in Angstrom)")
    text_box1.on_submit(submit)
    #text_box1.stop_typing()

    # Enter button
    axenter = plt.axes([0.7, 0.05, 0.1, 0.075])
    bnext = Button(axenter, 'Enter')
    bnext.on_clicked(enter)

    plt.show()

    # new_data has values of mid-wavelength of lines in Angstrom
    # indices has indices of the starting and ending point of lines
    # 
    # We want to fit a Gaussian to pixel-flux (with flux-error) data
    # to find out mid-wavelength. And then assign this wavelength (in
    # pixel space) to user provided wavelength stored in new_data.
    # 
    # mid_pix stores the value of pixel at mid-wavelength in wavelength space

    mid_pix = []

    for i in range(len(indices)):
        aa = int(indices[i][0])
        bb = int(indices[i][1])
        pix1 = pix[aa:bb]
        fl1 = fl[aa:bb]
        fle1 = fle[aa:bb]
        xinit = np.array([(pix1[0] + pix1[-1])/2, 1, 1, 1])
        def min_log_likelihood(x):
            model = utl.neg_gaus(pix1, x[0], x[1], x[2], x[3])
            chi2 = (fl1 - model)/fle1
            chi22 = np.sum(chi2**2)
            yy = np.sum(np.log(fle1)) + 0.5*chi22
            return yy
        soln = mz(min_log_likelihood, xinit, method='L-BFGS-B')
        mid_pix.append(soln.x[0])

    mid_wave_pix = np.asarray(mid_pix)
    mid_wave_lam = np.asarray(new_data)

    popt, pcov = cft(utl.line, mid_wave_pix, mid_wave_lam)

    return popt, pcov