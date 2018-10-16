# Lab 0 Data Analysis for part 4.4
# Makes a bad pixel Map
# September 7, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Load constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Open files
dark_f = fits.open(info['fitsSubdir']+info['masterDark'])
dark = dark_f[0].data
flat_f = fits.open(info['fitsSubdir']+info['masterFlat'])
flat = flat_f[0].data

# Identify hot and dead pixels
print('Creating bad pixel map')
dead_c = []
dead_r = []
hot_c = []
hot_r = []
for column in range( len( flat[0] ) ):
    for row in range( len( flat[1] ) ):
        if flat[column][row] < 0.8:
            dead_c.append(column)
            dead_r.append(row)
        elif dark[column][row] > 8000:
            hot_c.append(column)
            hot_r.append(row)
bad_pixel_map = np.ones( (len(flat[0]), len(flat[1]) ) )
for i in range( len( dead_c ) ):
    bad_pixel_map[dead_c[i]][dead_r[i]] = 0
for i in range( len( hot_c ) ): 
    bad_pixel_map[hot_c[i]][hot_r[i]] = 0
print( 'Number of bad pixels: ', (len(dead_c) + len(hot_c)) )
print( 'Percent of pixels that are bad: %.10f' % ( 100 * 
    (len(dead_c) + len(hot_c)) / ( len(flat[0])**2 ) ) )

# Saving the Bad Pixel map
print('Saving bad pixel map')
map_write = fits.PrimaryHDU(bad_pixel_map)
map_write.writeto(info['fitsSubdir']+info['badPixelMap'])
