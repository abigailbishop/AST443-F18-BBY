# Lab 0 Data Analysis for part 4.1
# September 7, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Open files
masterflat_f = fits.open('3.3_SpectrumFlat/masterflat.fits')
masterflat = masterflat_f[0].data
masterflat_crop = masterflat[175:223, 20:529]

# Save the cropped flat
masterflat_write = fits.PrimaryHDU(masterflat_crop)
#masterflat_write.writeto('3.3_SpectrumFlat/masterflat_50E-6m.fits')
