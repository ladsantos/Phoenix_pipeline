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
import calibration as cbr
import flux_extraction as fx


obs = np.array(['30_Oct_2009'])

path1 = os.getcwd()

for i in range(len(obs)):
    x_d = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/' + obs[i] + '/dark/'
    x_f = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/' + obs[i] + '/flat/'
    x_s = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/' + obs[i] + '/telluric_standard/'
    #x_s = '/home/jayshil/Documents/UNIGE/APL/APL1/WASP-7b_Phoenix/' + obs[i] + '/WASP-7/'
    it_s = 'object'
    cbr.calibrate_images(x_d = x_d, x_f = x_f, x_s = x_s, it_s = it_s)
    p1 = x_s + 'Final_calibrated_science/'
    p2 = x_s + 'Error_final_calibrated_science/'
    os.mkdir(x_s + 'Flux/')
    file_cali = os.listdir(p1)
    file_cali.sort(key = utl.natural_keys)
    file_cali_err = os.listdir(p2)
    file_cali_err.sort(key = utl.natural_keys)
    p3 = x_s + 'Flux/'
    for j in range(len(file_cali)):
        fx.flux_extraction(file_name = file_cali[j], file_err_name = file_cali_err[j], path = p1, path_err=p2, out_path = p3)
