import numpy as np
import astropy.nddata import CCDData
import ccdproc as ccdp
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from scipy.optimize import curve_fit as cft

file_name = input('Enter the name of fits file: ')

pt = Path('.')
f1 = ccdp.ImageFileCollection(pt)
ccd = CCDData.read(file_name + '.fits')

# Trimming the Image
trimmed = ccdp.trim_image(ccd, fits_section = '[1:256, 100:1000]')
trimmed.meta['TRIM'] = True
trimmed.header = ccd.header
trimmed.write(file_name + '_trim.fits')

# Reading the data from Trimmed image
data1 = trimmed.data
# Cutting the intial and last 50 columns from the data
data2 = np.transpose(data1)
data3 = data2[50:190]
data = np.transpose(data3)

# Detecting a line where spectrum could reside
ys = np.array([100, 200, 700, 800])
xs = np.array([])
for i in range(len(ys)):
	xd = data[ys[i]]
	ma = np.max(xd)
	ab = np.where(xd == ma)
	xs = np.hstack((xs, ab[0][0]))

def line(x, m, c):
	return m*x + c

popt, pcov = cft(line, xs, ys)

# Finding maximum flux at each point according to detected spectrum
def inv_line(x, m, c):
	bc = (x - c)/m
	return bc

# Finding total flux
def total_flux(lam, xlim=25):
	ydata = data[lam]
	xmid = inv_line(lam, popt[0], popt[1])
	xlow = xmid - xlim
	xup = xmid + xlim
	total_flux1 = 0
	xdata = np.arange(int(xlow), int(xup+1), 1)
	for i in range(len(xdata)):
		total_flux1 = total_flux1 + ydata[xdata[i]]
	return total_flux

# Flux as a function of pixel
flux = np.array([])
y11 = np.arange(0, 901, 1)
for i in range(len(y11)):
	f11 = total_flux(y11[i])
	flux = np.hstack((flux, f11))
