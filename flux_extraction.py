import numpy as np
from astropy.nddata import CCDData
import ccdproc as ccdp
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from scipy.optimize import curve_fit as cft


def flux_extraction(file_name, path, out_path, images=True):
	"""
	Parameters
	----------
	----------
	file_name : str
			Name of the image/telluric file
			from which flux has to be extracted
	path : str
			Path of the desired image file
	out_path : str
			Path of the output data and/or image file
	images : bool
			True if one wants to save visualization of flux data
			False if not.
			Default is True
	----------
	returns
	----------
	flux : data file
			.dat file containing the flux at
			various pixel values
			Path of this file would be similar to
			that of image file.
	----------
	"""
	pt = Path(path)
	f1 = ccdp.ImageFileCollection(pt)
	ccd = CCDData.read(path + file_name)# + '.fits')
	
	# Trimming the Image
	trimmed = ccdp.trim_image(ccd, fits_section = '[1:256, 100:1000]')
	trimmed.meta['TRIM'] = True
	trimmed.header = ccd.header
	#trimmed.write(file_name + '_trim.fits')
	
	# Reading the data from Trimmed image
	data = trimmed.data
	
	# Creating a function to detect the edges of slit
	# For lower edge
	def xlow(raw_data):
		"""
		Parameters
		----------
		----------
		raw_data : numpy.ndarray
				Array containing flux at some particular wavelength
		----------
		returns
		----------
		number : float
				A pixel number showing the lower edge of slit
		----------
		"""
		j = 0
		for i in range(int(len(raw_data)/5)):
			st = np.std(raw_data[j:j+5])
			xlw = 0
			if st < 2:
				xlw = j
			if xlw != 0:
				break
			j = j + 5
		return xlw
	# For upper edge
	def xup(raw_data):
		"""
		Parameters
		----------
		----------
		raw_data : numpy.ndarray
				Array containing flux at some particular wavelength
		----------
		returns
		----------
		number : float
				A pixel number showing the upper edge of slit
		----------
		"""
		j = 255
		for i in range(int(len(raw_data)/5)):
			st = np.std(raw_data[j-5:j])
			xup = 0
			if st < 2:
				xup = j
			if xup != 0:
				break
			j = j - 5
		return xup
	
	# Defining line and inverse line
	def line(x, m, c):
		return m*x + c
	def inv_line(x, m, c):
		bc = (x-c)/m
		return bc
	
	# Detecting the edges of the spectrum
	ys = np.array([150, 300, 450, 600, 750])
	xs_left = np.array([])
	xs_right = np.array([])
	xs_mid = np.array([])
	for i in range(len(ys)):
		dd1 = data[ys[i]]
		xll = xlow(dd1)
		xs_left = np.hstack((xs_left, xll))
		xuu = xup(dd1)
		xs_right = np.hstack((xs_right, xuu))
	
	popt_l, pcov_l = cft(line, xs_left, ys)
	popt_r, pcov_r = cft(line, xs_right, ys)

	# Detecting a line where spectrum could reside
	for i in range(len(ys)):
		ran_l = inv_line(ys[i], popt_l[0], popt_l[1])
		ran_r = inv_line(ys[i], popt_r[0], popt_r[1])
		xd1 = data[ys[i]]
		xd = xd1[int(ran_l):int(ran_r)]
		ma = np.max(xd)
		ab = np.where(xd == ma)
		xs_mid = np.hstack((xs_mid, ab[0][0] + ran_l))
	
	popt_m, pcov_m = cft(line, xs_mid, ys)
	
	# Finding total flux
	def total_flux(lam, xlim=20):
		ydata = data[lam]
		xmid = inv_line(lam, popt_m[0], popt_m[1])
		xlow = xmid - xlim
		xup = xmid + xlim
		total_flux1 = 0
		xdata = np.arange(int(xlow), int(xup+1), 1)
		for i in range(len(xdata)):
			total_flux1 = total_flux1 + ydata[xdata[i]]
		return total_flux1
	
	# Flux as a function of pixel
	flux = np.array([])
	y11 = np.arange(0, 900, 1)
	for i in range(len(y11)):
		f11 = total_flux(y11[i])
		flux = np.hstack((flux, f11))

	# Saving the image file for flux
	if images == True:
		fig1 = plt.figure(figsize = (20,10))
		plt.plot(flux)
		plt.xlabel('Pixel Number')
		plt.ylabel('Total Flux')
		plt.title('Total flux for ' + file_name + ' observation')
		plt.grid()
		plt.savefig(out_path + '/' + file_name + '_flux.png')
		plt.close(fig1)
	
	# Saving Data file of the flux
	f1 = open(out_path + '/' + file_name + '_flux.dat', 'w')
	f1.write('#Pixel\t\tFlux\n')
	for i in range(len(y11)):
		f1.write(str(y11[i]) + '\t\t' + str(flux[i]) + '\n')
	f1.close()
