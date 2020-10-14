import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
from astropy.nddata import CCDData
import ccdproc as ccdp

x_b = input('Enter the path of Bias images, if there are not one then press Enter: ')
x_d = input('Enter the path of Dark images: ')
x_f = input('Enter the path of Flat-field images: ')
x_s = input('Enter the path of Science images: ')

path_b = Path(x_b)
path_d = Path(x_d)
path_f = Path(x_f)
path_s = Path(x_s)

if x_b == '' and path_b == Path(''):
	print('Be aware: You did not provide the Bias files; the process will still continue though.')
	continue
if x_d == '':
	raise RuntimeError('You must provide Dark files for processing.')
if x_f == '':
	raise RuntimeError('You must provide Flatfield files for processing.')
if x_s == '':
	raise RuntimeError('You must provide Science images for processing.')
