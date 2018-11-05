"""

AST 443
Lab2
Data Analysis
7.2a Master Flat

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

# Open Files
flat_dir = info['dataDir'] + info['flatSubdir']
files = open(flat_dir+'names.txt', 'r')
flats = []
for line in files:
    flats.append(fits.open(flat_dir+line.strip('\n')))
flats_data = []
for flat in flats:
    flats_data.append(flat[0].data)

# Average the flats and remove the bias
flat_avg = np.median( flats_data, axis=0 )
flat_avg_mode = stats.mode(flat_avg.flatten())[0][0]
flat_avg_shape = flat_avg.shape
flat_avgnorm = np.zeros(flat_avg_shape)
for column in range(flat_avg_shape[1]):
    for row in range(flat_avg_shape[0]):
        flat_avgnorm[row][column] = flat_avg[row][column] / flat_avg_mode

# Save the Master flat
mfFname = info['fitsFiles'] + 'masterFlat.fits'
fits.writeto(mfFname,flat_avgnorm)

