# Data Analysis for part 4.1
# Calibrates Science Image (SI - Master_Dark) / Master_Flat
# October 9, 2018
# NOTES: Need to check work on imports and save calibrated images

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

# Open Master Dark and Flat
mdark_dir = info['fitsSubdir'] + info['masterDark']
mdark = open(mdarkdir)
mdark_data = mdark[0].data

mflat_dir = info['fitsSubdir'] + info['masterFlat']
mflat = open(mflat_dir)
mflat_data = mflat[0].data

# Open Image Files
sci_dir = info['dataDir'] + info[]
files = open(sci_dir+'names.txt', 'r')
images = []
for line in files:
    images.append(fits.open(sci_dir_line.strip('\n')))

# Calibrate and Save Images
for image in images:
    images_data = image[0].data
    cal_image = (image_data - mdark_data) / mflat_data


