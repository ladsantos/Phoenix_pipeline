import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy.stats import mad_std
import utils as utl
import astropy.units as u
import shutil

#-----------------------------------
#
#-------Loading data files----------
#
#-----------------------------------
#"""
x_b = ''
x_d = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/30_Oct_2009/dark/'
x_f = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/30_Oct_2009/flat/'
#x_s = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/30_Oct_2009/telluric_standard/'
x_s = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/30_Oct_2009/WASP-7/'
it_s = 'object'
"""
x_b = input('Enter the path of Bias images, if there are none then press Enter: ')
x_d = input('Enter the path of Dark images: ')
x_f = input('Enter the path of Flat-field images: ')
x_s = input('Enter the path of Science images: ')
it_s = input('Enter the Science Image file type: ')
"""

path_b = Path(x_b)
path_d = Path(x_d)
path_f = Path(x_f)
path_s = Path(x_s)

#-----Bias
if x_b == '' and path_b == Path(''):
	print('Be aware: You did not provide the Bias files; the process will still continue though.')
	files_b = None
elif not path_b.is_dir():
	raise RuntimeError('The path you provided for the Bias files does not exist.')
else:
	files_b = ccdp.ImageFileCollection(path_b)
#-----Dark
if x_d == '' or not path_d.is_dir():
	raise RuntimeError('You must provide Dark files for processing.\n Or the path you provided does not exist.')
else:
	files_d = ccdp.ImageFileCollection(path_d)
#-----Flat
if x_f == '' or not path_f.is_dir():
	raise RuntimeError('You must provide Flatfield files for processing.\n Or the path you provided does not exist.')
else:
	files_f = ccdp.ImageFileCollection(path_f)
#-----Science
if x_s == '' or not path_s.is_dir():
	raise RuntimeError('You must provide Science images for processing.\n Or the path you provided does not exist.')
else:
	files_s = ccdp.ImageFileCollection(path_s)

#-----------------------------------
#
#--------Calibrating Images---------
#
#-----------------------------------

if files_b is not None:
	#-------------------------------
	#------Creating Master-bias-----
	#-------------------------------
	cali_bias_path = Path(path_b / 'cali_bias')
	cali_bias_path.mkdir(exist_ok = True)
	files_b_cali = files_b.files_filtered(imagetyp = 'bias', include_path = True)
	combined_bias = ccdp.combine(files_b_cali, method='average', sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, mem_limit=350e6)
	combined_bias.meta['combined'] = True
	combined_bias.write(cali_bias_path / 'master_bias.fits')
	# Reading master bias
	master_bias = CCDData.read(cali_bias_path / 'master_bias.fits')
else:
	master_bias = None


#-------------------------------
#-------Calibrating Darks-------
#-------------------------------
cali_dark_path = Path(path_d / 'cali_dark')
cali_dark_path.mkdir(exist_ok = True)
files_d_cali = files_d.files_filtered(imagetyp = 'DARK', include_path = True)
for ccd, file_name in files_d.ccds(imagetyp = 'DARK', return_fname = True, ccd_kwargs = {'unit':'adu'}):
	if master_bias is not None:
		# Subtract bias
		ccd = ccdp.subtract_bias(ccd, master_bias)
	else:
		ccd = ccd
	# Save the result 
	ccd.write(cali_dark_path / file_name)


#--------------------------------
#------Creating Master-Dark------
#--------------------------------
red_dark = ccdp.ImageFileCollection(cali_dark_path)
# Calculating exposure times of DARK images
dark_times = set(red_dark.summary['exptime'][red_dark.summary['imagetyp'] == 'DARK'])
for exposure in sorted(dark_times):
	cali_darks = red_dark.files_filtered(imagetyp = 'dark', exptime = exposure, include_path = True)
	combined_dark = ccdp.combine(cali_darks, method='average', sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, mem_limit=350e6)
	combined_dark.meta['combined'] = True
	com_dark_name = 'combined_dark_{:6.3f}.fits'.format(exposure)
	combined_dark.write(cali_dark_path / com_dark_name)
# Reading master dark of various exposure times
red_dark = ccdp.ImageFileCollection(cali_dark_path)
combined_darks = {ccd.header['exptime']: ccd for ccd in red_dark.ccds(imagetyp = 'DARK', combined = True)}

#--------------------------------
#-------Calibrating Flats--------
#--------------------------------
cali_flat_path = Path(path_f / 'cali_flat')
cali_flat_path.mkdir(exist_ok = True)
files_f_cali = files_f.files_filtered(imagetyp = 'FLAT', include_path = True)
for ccd, file_name in files_f.ccds(imagetyp = 'FLAT', ccd_kwargs = {'unit' : 'adu'}, return_fname = True):
	# Subtract bias
	if master_bias is not None:
		ccd = ccdp.subtract_bias(ccd, master_bias)
	else:
		ccd = ccd
	closest_dark = utl.find_nearest_dark_exposure(ccd, dark_times)
	if closest_dark is None:
		closest_dark1 = utl.find_nearest_dark_exposure(ccd, dark_times, tolerance = 100)
		# Subtract scaled Dark
		ccd = ccdp.subtract_dark(ccd, combined_darks[closest_dark1], exposure_time = 'exptime', exposure_unit = u.second, scale = True)
		ccd.write(cali_flat_path / ('flat-' + file_name))
	else:
		closest_dark2 = utl.find_nearest_dark_exposure(ccd, dark_times)
		# Subtracting Darks
		ccd = ccdp.subtract_dark(ccd, combined_darks[closest_dark2], exposure_time = 'exptime', exposure_unit = u.second)
		ccd.write(cali_flat_path / ('flat-' + file_name))


#--------------------------------
#-----Creating Master-Flat-------
#--------------------------------
red_flats = ccdp.ImageFileCollection(cali_flat_path)
cali_flats = red_flats.files_filtered(imagetyp = 'FLAT', include_path = True)
combined_flat = ccdp.combine(cali_flats, method='average', scale = utl.inverse_median, sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, mem_limit=350e6)
combined_flat.meta['combined'] = True
combined_flat.write(cali_flat_path / 'master_flat.fits')
# Reading master flat
red_flats = ccdp.ImageFileCollection(cali_flat_path)
combined_flat = CCDData.read(cali_flat_path / 'master_flat.fits')#{ccd.header['imagetyp'] : ccd for ccd in red_flats.ccds(imagetyp = 'FLAT', combined = True)}


#--------------------------------
#---Calibrating Science Images---
#--------------------------------

# Creating a list of spectrum images
files_spec = files_s.summary['file', 'view_pos']
files_spec_list = np.array([])
for i in range(len(files_spec)):
	xxx = files_spec['view_pos'][i]
	if xxx[0:4] == 'open':
		files_spec_list = np.hstack((files_spec_list, files_spec['file'][i]))


# Sky subtracting images
cali_science_path1 = Path(path_s / 'Sky_subtracted_science')
cali_science_path1.mkdir(exist_ok = True)
j = 0
for i in range(int(len(files_spec_list)/2)):
	ccd1 = CCDData.read(x_s + files_spec_list[j], unit='adu')
	ccd2 = CCDData.read(x_s + files_spec_list[j+1], unit = 'adu')
	sky_sub1 = ccd1.data - ccd2.data
	ss1 = CCDData(sky_sub1, unit='adu')
	ss1.header = ccd1.header
	ss1.meta['sky_sub'] = True
	name1 = 'sky_sub_' + files_spec_list[j] + '.fits'
	ss1.write(cali_science_path1 / name1)
	sky_sub2 = ccd2.data - ccd1.data
	ss2 = CCDData(sky_sub2, unit='adu')
	ss2.header = ccd2.header
	ss2.meta['sky_sub'] = True
	name2 = 'sky_sub_' + files_spec_list[j+1] + '.fits'
	ss2.write(cali_science_path1 / name2)
	j = j+2

files_s1 = ccdp.ImageFileCollection(cali_science_path1)
final_calibrated = Path(path_s / 'Final_calibrated_science')
final_calibrated.mkdir(exist_ok = True)

# Correcting for flat
for ccd, file_name in files_s1.ccds(imagetyp = it_s, ccd_kwargs = {'unit' : 'adu'}, return_fname = True):
	# Subtract scaled Dark
	#ccd = ccdp.subtract_dark(ccd, combined_darks[closest_dark1], exposure_time = 'exptime', exposure_unit = u.second, scale = True)
	ccd = ccdp.flat_correct(ccd, combined_flat)#['FLAT'])
	ccd.write(final_calibrated / file_name)
"""
	else:
		closest_dark2 = utl.find_nearest_dark_exposure(ccd, dark_times)
		# Subtracting Darks
		ccd = ccdp.subtract_dark(ccd, combined_darks[closest_dark2], exposure_time = 'exptime', exposure_unit = u.second)
		ccd = ccdp.flat_correct(ccd, combined_flat)#['FLAT'])
		ccd.write(cali_science_path / file_name)
		ccd.write(final_calibrated / file_name)
"""
