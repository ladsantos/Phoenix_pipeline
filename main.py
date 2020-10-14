import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from astropy.nddata import CCDData
import ccdproc as ccdp


#-----------------------------------
#
#-------Loading data files----------
#
#-----------------------------------
x_b = input('Enter the path of Bias images, if there are none then press Enter: ')
x_d = input('Enter the path of Dark images: ')
x_f = input('Enter the path of Flat-field images: ')
x_s = input('Enter the path of Science images: ')

path_b = Path(x_b)
path_d = Path(x_d)
path_f = Path(x_f)
path_s = Path(x_s)

#-----Bias
if x_b == '' and path_b == Path(''):
	print('Be aware: You did not provide the Bias files; the process will still continue though.')
	files_b = None
else if not path_b.is_dir():
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
#------Checking for overscan--------
#
#-----------------------------------

oscn = np.array([])
for f1 in files_s.hdus():
	xx = 'biassec' in a.header
	if xx:
		oscn = np.hstack((oscn, a.header['biassec']))
	else:
		oscn = np.hstack((oscn, xx))

if oscn[0] == 0:
	print('No Overscan detected!')
	osn1 = input('If there is overscan, but not detected, then please type the region of overscan: ')
	osn = True
	if osn = '':
		osn = False
else:
	osn1 = input('Please input the useful region of overscan: ')
	osn = True


#-----------------------------------
#
#--------Calibrating Images---------
#
#-----------------------------------
if file_b is not None:
	if osn:
		
