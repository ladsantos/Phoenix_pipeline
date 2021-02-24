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
import wavelength_soln as wsoln
from tkinter import filedialog
from tkinter import *
import glob
from astropy.io import fits

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
list3 = glob.glob(pt_s + '*.fits')
imgtyp = []
for i in range(len(list3)):
    hdul = fits.open(list3[i])
    h1 = hdul[0].header
    if h1['view_pos'][0:4] == 'open':
        h22 = h1['imagetyp']
        imgtyp.append(h22)

it_s = imgtyp[0]

#-------------------------------------------
#-------- For Telluric files ---------------
#-------------------------------------------

print('                                                   ')
print('For Telluric files:')
print('-------------------')
print('                   ')
input('Press Enter to choose the directory of Telluric files...')
pt_t1 = filedialog.askdirectory(title='Choose Telluric files directory')
pt_t = pt_t1 + '/'

'''
list5 = os.listdir(pt_t)
print(list5)
for i in range(len(list5)):
    os.system('mv ' + pt_t + list5[i] + ' ' + pt_t + 'telluric_' + list5[i])
'''
list4 = glob.glob(pt_t + '*.fits')
imgtyp1 = []
for i in range(len(list4)):
    hdul = fits.open(list4[i])
    h1 = hdul[0].header
    if h1['view_pos'][0:4] == 'open':
        h22 = h1['imagetyp']
        imgtyp1.append(h22)

it_t = imgtyp1[0]

#-------------------------------------------
#--------- For Output files ----------------
#-------------------------------------------

print('                                                   ')
print('For Output files:')
print('-----------------')
print('                  ')
print('To save output files, the program will create')
print('a folder named Output within your chosen directory.')
input('Press Enter to choose the directory to dump output files...')
pt_out1 = filedialog.askdirectory()
try:
    os.mkdir(pt_out1 + '/Output')
except:
    pt_out = pt_out1 + '/Output'
pt_out = pt_out1 + '/Output'

#-------------------------------------------
#------------- Calibration -----------------
#-------------------------------------------
print('                                                   ')
print('Calibration Starting...')
print('For Tellurics...')
print('-----------------------')
print('                  ')
cbr.calibrate_images(x_d=pt_d, x_f=pt_f, x_s=pt_t, it_s=it_t, x_b=pt_b)
os.system('rm -r ' + pt_d + 'cali_dark')
os.system('rm -r ' + pt_f + 'cali_flat')

if pt_b != '':
    os.system('rm -r ' + pt_b + 'cali_bias')

print('For Science...')
print('-----------------------')
print('                  ')
cbr.calibrate_images(x_d=pt_d, x_f=pt_f, x_s=pt_s, it_s=it_s, x_b=pt_b)
print('                   ')
print('Calibration done...')
print('-------------------')
print('                   ')

# Path of calibrated Science files
p1 = pt_s + 'Final_calibrated_science/'
p2 = pt_s + 'Error_final_calibrated_science/'

# Path of calibrated Telluric files
p1_t = pt_t + 'Final_calibrated_science/'
p2_t = pt_t + 'Error_final_calibrated_science/'

#-------------------------------------------
#----------- Flux Extraction ---------------
#-------------------------------------------
print('Flux Extraction Starting...')
print('---------------------------')
print('                           ')
# For Tellurics
print('For Tellurics...')
# Flux directory
os.mkdir(pt_t + 'Flux_t/')
file_cali_t = os.listdir(p1_t)
file_cali_t.sort(key = utl.natural_keys)
file_cali_err_t = os.listdir(p2_t)
file_cali_err_t.sort(key = utl.natural_keys)
# Simple Flux Extraction
p3_t = pt_t + 'Flux_t/'
for j in range(len(file_cali_t)):
    fx.flux_extraction(file_name = file_cali_t[j], file_err_name = file_cali_err_t[j], path = p1_t, path_err=p2_t, out_path = p3_t, images=False)

# For Science
print('For Science...')
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
print('For Tellurics...')
flux_files_t = os.listdir(p3_t)
flux_files_t.sort(key=utl.natural_keys)
for i in range(len(flux_files_t)):
    nml.normal_spectrum(p3_t, flux_files_t[i], p3_t)

nml_files1_t = os.listdir(p3_t)
nml_files_t = []
for i in range(len(nml_files1_t)):
    if nml_files1_t[i][0:6] == 'normal':
        nml_files_t.append(nml_files1_t[i])

nml_files_t.sort(key=utl.natural_keys)
csr.cosmic_removal(p3_t, nml_files_t, p3_t)

print('                                                      ')
print('For Science...')
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
os.system('mv ' + pt_t + 'cali_science ' + pt_out + '/cali_telluric')
os.system('mv ' + p1 + ' ' + pt_out)
os.system('mv ' + p2 + ' ' + pt_out)
os.system('mv ' + p3 + ' ' + pt_out)
os.system('mv ' + p1_t + ' ' + pt_out + '/Final_calibrated_telluric')
os.system('mv ' + p2_t + ' ' + pt_out + '/Error_final_calibrated_telluric')
os.system('mv ' + p3_t + ' ' + pt_out)


#-------------------------------------------
#--- Wavelength Calibration in Tellurics ---
#-------------------------------------------

pt_tel = pt_out +'/Flux_t/'
list5 = os.listdir(pt_tel)
csr_tel = []
for i in range(len(list5)):
    if list5[i][0:6] == 'cosmic':
        csr_tel.append(list5[i])

soln_m = np.array([])
soln_c = np.array([])

for i in range(len(csr_tel)):
    popt, pcov = wsoln.wave_soln(pt_tel, csr_tel[i])
    soln_m = np.hstack((soln_m, popt[0]))
    soln_c = np.hstack((soln_c, popt[1]))

# Linear Solution of the Wavelength Calibration

avg_soln_m = np.mean(soln_m)
avg_soln_c = np.mean(soln_c)

# Applying this to Science files
pt_sci = pt_out + '/Flux/'
list6 = os.listdir(pt_sci)
csr_sci = []
for i in range(len(list6)):
    if list6[i][0:6] == 'cosmic':
        csr_sci.append(list6[i])

# Creating a directory to save Spectrum files
os.mkdir(pt_out + '/Spectrum')
pt_spectrum = pt_out + '/Spectrum/'

for i in range(len(csr_sci)):
    pix, fl, fle = np.loadtxt(pt_sci + csr_sci[i], usecols=(0,1,2), unpack=True)
    wave = utl.line(pix, avg_soln_m, avg_soln_c)
    f1 = open(pt_spectrum + 'spectra_' + csr_sci[i][22:], 'w')
    f1.write('#Wavelength (in A)\t Flux \t Flux Error\n')
    for j in range(len(wave)):
        f1.write(str(wave[j]) + '\t' + str(fl[j]) + '\t' + str(fle[j]) + '\n')
    f1.close()

print('This pipleline is currently under developement...')
print('The latest stage is the extraction of normalized')
print('spectrum with cosmic rays removed...')
print('These files you can access in the Output folder, Flux directory...')