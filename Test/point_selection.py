import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.widgets import Cursor

pt = os.getcwd()
print(pt)

pix, fl , fle = np.loadtxt(pt + '/Test/cosmic_normal_sky_sub_2009oct30_0009.fits_flux.dat', usecols=(0,1,2), unpack=True)


plt.errorbar(pix, fl, yerr=fle)


data = plt.ginput(n=-1)
print(data)