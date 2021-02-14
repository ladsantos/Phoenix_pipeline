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
import cosmic_rays as cs
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
print('         Email: Jayshil.Patel@etu.unige.ch         ')
print('                Leonardo.DosSantos@unige.ch        ')
print('                                                   ')
print('***************************************************')
print('                                                   ')
print('For Bias files:')
print('---------------')
print('               ')
xx = input('Do you want to provide Bias files (y/n)?:')
print('                                          ')
if xx == 'Y' or 'y':
    yy = input('Press Enter to choose the directory of Bias file')
    x_b = filedialog.askdirectory()
else:
    x_b = ''