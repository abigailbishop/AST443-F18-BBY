# Data Analysis for part 4.1
# Calibrates Science Image (SI - Master_Dark) / Master_Flat
# October 9, 2018
# NOTES: Need calImage in constants, change to Emily's directory

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

# Open Master Dark, Master Flat, Bad Pixel Map
mdark_dir = info['fitsSubdir'] + info['masterDark']
mdark = fits.open(mdark_dir)
mdark_data = mdark[0].data

mflat_dir = info['fitsSubdir'] + info['masterFlat']
mflat = fits.open(mflat_dir)
mflat_data = mflat[0].data

bpm_dir = infor['fitsSubdir'] + info['badPixelMap']
bpm = fits.open(bpm_dir)
bpm_data = bpm[0].data

# Open Image Files
sci_dir = info['dataDir']
files = open(sci_dir+'names.txt', 'r')
images = []
for line in files:
    images.append(fits.open(sci_dir_line.strip('\n')))
files.close()

# Calibrate and Save Images
i=0
for image in images:
    fname = files[i]
    images_data = image[0].data
    cal_image = (image_data - mdark_data) * bpm_data / mflat_data
    cal_image_write = fits.PrimaryHDU(cal_image)
    cal_image_write.writeto(info['calImage'] + fname
    i=i+1


