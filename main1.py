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
import normalization as nml
import cosmic_rays as csr
from tkinter import filedialog
from tkinter import *

root = Tk()
root.withdraw()

print('***************************************************')
print('                                                   ')
print('          Phoenix Instrument Pipeline              ')
print('          ---------------------------              ')
print('                                                   ')
print('                  Created by                       ')
print('      Jayshil A Patel and Leonardo Dos Santos      ')
print('                                                   ')
print('       Email: Jayshil.Patel at etu.unige.ch        ')
print('              Leonardo.DosSantos at unige.ch       ')
print('                                                   ')
print('***************************************************')
print('                                                   ')

#-------------------------------------------
#-------- For Bias files -------------------
#-------------------------------------------

print('For Bias files:')
print('---------------')
print('               ')
xx = input('Do you want to provide Bias files (y/n)?: ')
print('                                          ')
if xx == 'Y' or xx == 'y':
    input('Press Enter to choose the directory of Bias files...')
    pt_b1 = filedialog.askdirectory(title='Choose Bias files directory')
    pt_b = x_b1 + '/'
else:
    pt_b = ''

#-------------------------------------------
#-------- For Dark files -------------------
#-------------------------------------------

print('                                                   ')
print('For Dark files:')
print('---------------')
print('               ')
input('Press Enter to choose the directory of Dark files...')
pt_d1 = filedialog.askdirectory(title='Choose Dark files directory')
pt_d = pt_d1 + '/'

#-------------------------------------------
#-------- For Flat files -------------------
#-------------------------------------------

print('                                                   ')
print('For Flat files:')
print('---------------')
print('               ')
input('Press Enter to choose the directory of Flat files...')
pt_f1 = filedialog.askdirectory(title='Choose Flat files directory')
pt_f = pt_f1 + '/'

#-------------------------------------------
#-------- For Science files ----------------
#-------------------------------------------

print('                                                   ')
print('For Science files:')
print('------------------')
print('                  ')
input('Press Enter to choose the directory of Science files...')
pt_s1 = filedialog.askdirectory(title='Choose Science files directory')
pt_s = pt_s1 + '/'
it_s = input('Please enter the type of the Science files (Object/Science): ')

#-------------------------------------------
#--------- For Output files ----------------
#-------------------------------------------

print('                                                   ')
print('For Output files:')
print('-----------------')
print('                  ')
print('To save output files, the program will create')
print('a folder named Output within your choosen directory.')
input('Press Enter to choose the directory to dump output files...')
pt_out1 = filedialog.askdirectory()
os.mkdir(pt_out1 + '/Output')
pt_out = pt_out1 + '/Output'

#-------------------------------------------
#------------- Calibration -----------------
#-------------------------------------------
print('                                                   ')
print('Calibration Starting...')
print('-----------------------')
print('                  ')
cbr.calibrate_images(x_d=pt_d, x_f=pt_f, x_s=pt_s, it_s=it_s, x_b=pt_b)
print('Calibration done...')
print('-------------------')
print('                   ')

# Path of calibrated Science files
p1 = pt_s + 'Final_calibrated_science/'
p2 = pt_s + 'Error_final_calibrated_science/'

#-------------------------------------------
#----------- Flux Extraction ---------------
#-------------------------------------------
print('Flux Extraction Starting...')
print('---------------------------')
print('                           ')
# Flux directory
os.mkdir(pt_s + 'Flux/')
file_cali = os.listdir(p1)
file_cali.sort(key = utl.natural_keys)
file_cali_err = os.listdir(p2)
file_cali_err.sort(key = utl.natural_keys)
# Simple Flux Extraction
p3 = pt_s + 'Flux/'
for j in range(len(file_cali)):
    fx.flux_extraction(file_name = file_cali[j], file_err_name = file_cali_err[j], path = p1, path_err=p2, out_path = p3, images=False)

print('Flux Extraction done...')
print('-----------------------')
print('                       ')

#-------------------------------------------
# Flux Normalization and Cosmic Ray removal 
#-------------------------------------------
print('Flux Normalization and Cosmic Rays removal starting...')
print('------------------------------------------------------')
print('                                                      ')
flux_files = os.listdir(p3)
flux_files.sort(key=utl.natural_keys)
for i in range(len(flux_files)):
    nml.normal_spectrum(p3, flux_files[i], p3)

nml_files1 = os.listdir(p3)
nml_files = []
for i in range(len(nml_files1)):
    if nml_files1[i][0:6] == 'normal':
        nml_files.append(nml_files1[i])

nml_files.sort(key=utl.natural_keys)
csr.cosmic_removal(p3, nml_files, p3)

print('Flux Normalization and Cosmic Rays removal done...')
print('--------------------------------------------------')
print('                                                  ')

#----------------------------------------------
# Dumping all of the files in output directory 
#----------------------------------------------
if pt_b != '':
    os.system('mv ' + pt_b + 'cali_bias ' + pt_out)


os.system('mv ' + pt_d + 'cali_dark ' + pt_out)
os.system('mv ' + pt_f + 'cali_flat ' + pt_out)
os.system('mv ' + pt_s + 'cali_science ' + pt_out)
os.system('mv ' + p1 + ' ' + pt_out)
os.system('mv ' + p2 + ' ' + pt_out)
os.system('mv ' + p3 + ' ' + pt_out)

print('This pipleline is currently under developement...')
print('The latest stage is the extraction of normalized')
print('spectrum with cosmic rays removed...')
print('These files you can access in the Output folder, Flux directory...')