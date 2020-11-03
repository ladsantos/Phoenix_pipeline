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
x_d = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/20_Oct_2009/dark/'
x_f = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/20_Oct_2009/flat/'
#x_s = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/30_Oct_2009/telluric_standard/'
x_s = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/20_Oct_2009/WASP-7/'
it_s = 'object'
"""
x_b = input('Enter the path of Bias images, if there are none then press Enter: ')
x_d = input('Enter the path of Dark images: ')
x_f = input('Enter the path of Flat-field images: ')
x_s = input('Enter the path of Science images: ')
it_s = input('Enter the Science Image file type: ')
"""
