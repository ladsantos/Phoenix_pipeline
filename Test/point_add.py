import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Button
import os
from scipy.optimize import minimize as mz
import utils as utl

fig, ax = plt.subplots(figsize=(16, 12))
fig.subplots_adjust(bottom=0.2)

pt = os.getcwd()

t, y, ye = np.loadtxt(pt + '/Test/cosmic_normal_sky_sub_2009oct30_0009.fits_flux.dat', usecols=(0,1, 2), unpack=True)

l = ax.plot(t, y, lw=2)
ax.set_title('Press left mouse button and drag to test')

new_data = []
indices = []

# Dummy variables to save data
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

print(new_data)
print(indices)

# new_data has values of mid-wavelength of lines in Angstrom
# indices has indices of the starting and ending point of lines
# 
# We want to fit a Gaussian to pixel-flux (with flux-error) data
# to find out mid-wavelength. And then assign this wavelength (in
# pixel space) to user provided wavelength stored in new_data.

mid_pix = []


for i in range(len(indices)):
    pix1 = t[int(indices[i][0]):int(indices[i][1])]
    fl1 = y[int(indices[i][0]):int(indices[i][1])]
    fle1 = ye[int(indices[i][0]):int(indices[i][1])]
    xinit = np.array([(pix1[0] + pix1[-1])/2, 1, 1, 1])
    def min_log_likelihood(x):
        global pix1, fl1, fle1
        model = utl.neg_gaus(pix1, x[0], x[1], x[2], x[3])
        chi2 = (fl1 - model)/fle1
        chi22 = np.sum(chi2**2)
        yy = np.sum(np.log(fle1)) + 0.5*chi22
        return yy
    soln = mz(min_log_likelihood, xinit, method='L-BFGS-B')
    plt.plot(pix1, fl1)
    plt.plot(pix1, utl.neg_gaus(pix1, *soln.x))
    plt.show()
    mid_pix.append(soln.x[0])

print(mid_pix)