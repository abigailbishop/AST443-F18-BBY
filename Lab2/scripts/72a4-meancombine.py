# Data Analysis for part 4.1
# Mean combines the data files
# November 4, 2018

# Imports
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
data_dir = info['dataDir'] + info['flatCalibrated']
files = ['nebula_00000000_dcc_flat.fits', 
         'nebula_00000001_dcc_flat.fits',
         'nebula_00000002_dcc_flat.fits']
datas = []
for fileName in files:
    datas.append(fits.open(data_dir+fileName))
datas_data = []
for data in datas:
    datas_data.append(data[0].data)

# Create the Dark Master Frame
data_mean = np.median(datas_data, axis=0)
master_write = fits.PrimaryHDU(data_mean)
master_write.writeto(info['fitsFiles'] + 'nebula_avg.fits')
