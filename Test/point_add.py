import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from matplotlib.widgets import SpanSelector
from matplotlib.widgets import Button
import os

fig, ax = plt.subplots(figsize=(16, 12))
fig.subplots_adjust(bottom=0.2)

pt = os.getcwd()

t, y = np.loadtxt(pt + '/Test/cosmic_normal_sky_sub_2009oct30_0009.fits_flux.dat', usecols=(0,1), unpack=True)

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