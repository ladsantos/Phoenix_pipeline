import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Button
import os

def wave_soln(path, fname, out_path):
    """
    Parameters:
    -----------
    path : str
        path of the input file
    fname : str
        name of the input telluric file.
        Ideally this should be the file from which cosmic rays 
        are removed. 3 columns, pixel, flux, error in flux
    out_path : str
        path of the output file
    -----------
    returns
    -----------
    wavelength solution : data file
        .dat file containing pixel number and corresponding
        wavelength in angstrom.
    -----------
    """
    fig, ax = plt.subplots(figsize=(16,12))
    fig.subplots_adjust(bottom=0.2)

    pix, fl, fle = np.loadtxt(path + fname, usecols=(0,1,2), unpack=True)
    plts = ax.plot(pix, fl)
    ax.set_title('Press the left mouse button to select the region where the line is.')

    new_data = [] # wavelength data provided by user
    indices = [] # corresponding indices of starting and ending of line

    # Dummy variables
    new_data1 = []
    indices1 = []

    def onselect(tmin, tmax):
        global indices1
        indmin, indmax = np.searchsorted(t, (tmin, tmax))
        indmax = min(len(t) - 1, indmax)

        ab = [indmin, indmax]
        indices1 = []
        indices1.append(ab)

    def submit(expression):
        global new_data1
        ydata = expression
        new_data1 = []
        new_data1.append(ydata)

    def enter(self):
        new_data.append(new_data1[0])
        indices.append(indices1[0])
        text_box.set_val("")

    # To select a region
    span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                        rectprops=dict(alpha=0.5, facecolor='red'))

    # Adding a box to enter values
    axbox = fig.add_axes([0.2, 0.05, 0.4, 0.075])
    text_box = TextBox(axbox, "Enter the corresponding\n wavelength here (in Angstrom)")
    text_box.on_submit(submit)
    text_box.stop_typing()

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
