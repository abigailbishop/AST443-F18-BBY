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

bpm_dir = info['fitsSubdir'] + info['badPixelMap']
bpm = fits.open(bpm_dir)
bpm_data = bpm[0].data

# Open Image Files
sci_dir = info['dataDir'] + info['rawDataSubdir']
files = open(sci_dir+'names.txt', 'r')
image_data = []
fnames = []
for line in files:
    image = fits.open(sci_dir + line.strip('\n'))
    image_data.append(image[0].data)
    fnames.append(line.strip('\n'))
    image.close()

# Calibrate and Save Images
cal_image_dir = info['dataDir'] + info['calImageSubdir']
i=0
for image in image_data:
    fname = fnames[i]
    cal_image = (image - mdark_data) * bpm_data / mflat_data
    cal_image_write = fits.PrimaryHDU(cal_image)
    cal_image_write.writeto(cal_image_dir + fname)
    i=i+1
    
files.close()


