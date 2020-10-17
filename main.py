import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy.stats import mad_std


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
#--------Calibrating Images---------
#
#-----------------------------------
#---If bias file exists---
if file_b is not None:
	#-------------------------------
	#------Creating Master-bias-----
	#-------------------------------
	cali_bias_path = Path(path_b / 'master_bias')
	cali_bias_path.mkdir(exist_ok=True)
	combined_bias = ccdp.combine(files_b, method='average', sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std, mem_limit=350e6)
	combined_bias.meta['combined'] = True
	combined_bias.write(cali_bias_path / 'master_bias.fit')
	
