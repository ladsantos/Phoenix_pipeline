import numpy as np
from astropy.nddata import CCDData
import ccdproc as ccdp
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from scipy.optimize import curve_fit as cft


def flux_extraction(file_name, file_err_name, path, path_err, out_path, images=True):
	"""
	Parameters
	----------
	----------
	file_name : str
			Name of the image/telluric file
			from which flux has to be extracted
	file_err_name: str
			Name of the image/telluric variance file
	path : str
			Path of the desired image file
	path_err : str
			Path of the variance of desired image file
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
	# Reading Data File
	ccd = CCDData.read(path + file_name)# + '.fits')
	
	# Trimming the Image
	trimmed = ccdp.trim_image(ccd, fits_section = '[1:256, 100:1000]')
	trimmed.meta['TRIM'] = True
	trimmed.header = ccd.header
	#trimmed.write(file_name + '_trim.fits')
	
	# Reading the data from Trimmed image
	data = trimmed.data
	
	# Reading Variance File
	ccd_err = CCDData.read(path_err + file_err_name)# + '.fits')
	
	# Trimming the Image
	trimmed_err = ccdp.trim_image(ccd_err, fits_section = '[1:256, 100:1000]')
	trimmed_err.meta['TRIM'] = True
	trimmed_err.header = ccd.header
	#trimmed.write(file_name + '_trim.fits')
	
	# Reading the data from Trimmed image
	data_err = trimmed_err.data
	data_err[data_err == 0] = 1
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
	
	# Creating xdata and ydata in range of ccd
	xall = np.arange(0,256,1)
	yall = np.arange(0, 901, 1)
	
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

	# Defining a Gaussian to create Spatial Image
	def gaus(x, mu, sigma):
		a1 = np.sqrt(2*np.pi*sigma*sigma)**-1
		a2 = np.exp(-0.5*((x-mu)/sigma)**2)
		return a1*a2
	
	# Finding spatial profile
	def spatial(data1, data1_err, lam, xlim=25):
		ydata = data1[lam]
		xmid = inv_line(lam, popt_m[0], popt_m[1])
		xlow = xmid - xlim
		xup = xmid + xlim
		p2 = ydata[int(xlow):int(xup)]
		p1 = p2/np.sum(p2)
		xdata = np.arange(1, len(p1)+1, 1)
		poptg, pcovg = cft(gaus, xdata, p1)
		fwhm = np.sqrt(poptg[1]*poptg[1]*np.log(256))
		mu1 = poptg[0] + inv_line(lam, *popt_m) - xlim
		return mu1, poptg[1], fwhm
	
	# Finding total flux
	def total_flux(data1, data1_err, lam):
		ydata = data1[lam]
		ydata_err = data1_err[lam]
		mu1, sig1, fm1 = spatial(data1, data1_err, lam)
		p_x = gaus(xall, mu1, sig1)
		a1 = 0
		a2 = 0
		for i in range(len(xall)):
			a11 = p_x[i]*ydata[i]/ydata_err[i]
			a1 = a1 + a11
			a22 = p_x[i]*p_x[i]/ydata_err[i]
			a2 = a2 + a22
		fopt = a1/a2
		var = 1/a2
		return fopt, var

	# Cosmic Rays Removal
	def cosmic_ray(data1, data1_err, lam, threshold=16):
		fopt, var = total_flux(data1, data1_err, lam)
		d_s = data1[lam]
		d_s_err = data1_err[lam]
		mu2, sig2, fm2 = spatial(data1, data1_err, lam)
		p_x = gaus(xall, mu2, sig2)
		data_wo_cr = np.array([])
		data_wo_cr_err = np.array([])
		for i in range(len(xall)):
			xx = (d_s[i] - fopt*p_x[i])**2
			yy = xx/d_s_err[i]
			if yy>threshold:
				if i == len(xall)-1:
					xxx = data1[lam][i-1]
				else:
					xxx = (data1[lam][i-1] + data1[lam][i+1])/2
				data_wo_cr = np.hstack((data_wo_cr, xxx))
				data_wo_cr_err = np.hstack((data_wo_cr_err, 1))
			else:
				data_wo_cr = np.hstack((data_wo_cr, data1[lam][i]))
				data_wo_cr_err = np.hstack((data_wo_cr_err, data1_err[lam][i]))
		return data_wo_cr, data_wo_cr_err

	# Data Without Cosmic Rays
	final_data = np.array([])
	final_data_err = np.array([])
	final_data, final_data_err = cosmic_ray(data, data_err, yall[0], threshold=10)
	
	for i in range(len(yall)-1):
		fda, fdae = cosmic_ray(data, data_err, yall[i+1])
		final_data = np.vstack((final_data, fda))
		final_data_err = np.vstack((final_data_err, fdae))
	
	# Flux as a function of pixel
	flux = np.array([])
	flux_err = np.array([])
	for i in range(len(yall)):
		f11, v11 = total_flux(final_data, final_data_err, yall[i])
		flux = np.hstack((flux, f11))
		flux_err = np.hstack((flux_err, v11))

	# Saving the image file for flux
	if images == True:
		fig1 = plt.figure(figsize = (20,10))
		plt.errorbar(yall, flux, yerr=flux_err)
		plt.xlabel('Pixel Number')
		plt.ylabel('Total Flux')
		plt.title('Total flux for ' + file_name + ' observation')
		plt.grid()
		plt.savefig(out_path + '/' + file_name + '_flux.png')
		plt.close(fig1)
	
	# Saving Data file of the flux
	f1 = open(out_path + '/' + file_name + '_flux.dat', 'w')
	f1.write('#Pixel\t\tFlux\n')
	for i in range(len(yall)):
		f1.write(str(yall[i]) + '\t\t' + str(flux[i]) + '\t' + str(flux_err[i]) + '\n')
	f1.close()
