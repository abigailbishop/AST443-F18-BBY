"""

AST 443
Lab2
Data Analysis
7.2a2 Dark Current Correction

"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Load Lab constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Get Master Flat
flatDir = info['fitsFiles']
masterflatHDU = fits.open(flatDir + 'masterFlat.fits')
masterFlat = masterflatHDU[0].data
masterflatHDU.close()

# Get dcc Science Images Names
sciDirNames = ['dccNebulaData', 'dccRefstarData']
SciNames = []
for subdir in sciDirNames:
    sciDir = info['dataDir'] + info[subdir]
    files = open(sciDir + 'names.txt', 'r')
    for line in files:
        SciNames.append(sciDir + line.strip('\n'))

# Calibrate Images
for fname in SciNames:
    hdu = fits.open(fname)
    head = hdu[0].header
    sciExpTime = head['EXPTIME']
    data = hdu[0].data
    hdu.close()
    cordata = data/ masterFlat
    calfname = fname.strip('_dcc.FIT') + '_cal.fits'
    fits.writeto(calfname,cordata,header=head)
