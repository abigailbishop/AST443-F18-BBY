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

# Obtain information for dead pixels
print('Analyzing flat field for dead pixel info')
flat_flat = flat.flatten()
flat_flat_sigclip = stats.sigmaclip(flat_flat, 5.0, 5.0)
flat_flat_cut = flat_flat_sigclip[0]
flat_cut_min = flat_flat_sigclip[1]
flat_cut_max = flat_flat_sigclip[2]
flat_cut_mean = np.mean(flat_flat_cut)
flat_cut_mode = stats.mode(flat_flat_cut)[0][0]
flat_cut_stddev = np.std(flat_flat_cut)
flat_cut_5sigma = 5 * flat_cut_stddev

# Obtain information for hot pixels
print('Analyzing dark field for hot pixel info')
dark_adjusted = dark - flat
dark_adjusted_flat = dark_adjusted.flatten()
dark_flat_sigclip = stats.sigmaclip(dark_adjusted_flat, 5.0, 5.0)
dark_flat_cut = dark_flat_sigclip[0]
dark_cut_min = dark_flat_sigclip[1]
dark_cut_max = dark_flat_sigclip[2]
dark_cut_mean = np.mean(dark_flat_cut)
dark_cut_mode = stats.mode(dark_flat_cut)[0][0]
dark_cut_stddev = np.std(dark_flat_cut)
dark_cut_5sigma = 5 * dark_cut_stddev

# Identify hot and dead pixels
print('Creating bad pixel map')
dead_c = []
dead_r = []
hot_c = []
hot_r = []
for column in range( len( flat[0] ) ):
    for row in range( len( flat[1] ) ):
        if flat[column][row] < flat_cut_5sigma:
            dead_c.append(column)
            dead_r.append(row)
        elif dark_adjusted[column][row] > dark_cut_5sigma:
            hot_c.append(column)
            hot_r.append(row)
bad_pixel_map = np.ones( (len(flat[0]), len(flat[1]) ) )
for i in range( len( dead_c ) ):
    bad_pixel_map[dead_c[i]][dead_r[i]] = 0
for i in range( len( hot_c ) ): 
    bad_pixel_map[hot_c[i]][hot_r[i]] = 0
print( 'Percent of pixels that are bad: %.2f' % ( 100 * 
    (len(dead_c) + len(hot_c)) / ( len(flat[0])**2 ) ) )

# Save the master DARK minus the BIAS
print('Saving bad pixel map')
map_write = fits.PrimaryHDU(bad_pixel_map)
map_write.writeto(info['fitsSubdir']+info['badPixelMap'])
