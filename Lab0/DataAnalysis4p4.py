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

# Open files
bias_f = fits.open('neg10BIAS_master.fits')
bias = bias_f[0].data
dark_f = fits.open('neg10DARK_adjusted.fits')
dark = dark_f[0].data
flat_f = fits.open('flat_master.fits')
flat = flat_f[0].data

# Obtain information for dead pixels
bias_flat = bias.flatten()
bias_flat_sigclip = stats.sigmaclip(bias_flat, 5.0, 5.0)
bias_flat_cut = bias_flat_sigclip[0]
bias_cut_min = bias_flat_sigclip[1]
bias_cut_max = bias_flat_sigclip[2]
bias_cut_mean = np.mean(bias_flat_cut)
bias_cut_mode = stats.mode(bias_flat_cut)[0][0]
bias_cut_stddev = np.std(bias_flat_cut)
bias_cut_5sigma = 5 * bias_cut_stddev

# Obtain information for hot pixels
dark_adjusted = dark - bias
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
dead_c = []
dead_r = []
hot_c = []
hot_r = []
for column in range( len( bias[0] ) ):
    for row in range( len( bias[1] ) ):
        if bias[column][row] < bias_cut_5sigma:
            dead_c.append(column)
            dead_r.append(row)
        elif dark_adjusted[column][row] > dark_cut_5sigma:
            hot_c.append(column)
            hot_r.append(row)
bad_pixel_map = np.zeros( (len(bias[0]), len(bias[1]) ) )
for i in range( len( dead_c ) ):
    bad_pixel_map[dead_c[i]][dead_r[i]] = 1
for i in range( len( hot_c ) ): 
    bad_pixel_map[hot_c[i]][hot_r[i]] = 1
print( 'Percent of pixels that are bad: %.2f' % ( 100 * 
    (len(dead_c) + len(hot_c)) / ( len(bias[0])**2 ) ) )

# Save the master DARK minus the BIAS
map_write = fits.PrimaryHDU(bad_pixel_map)
#map_write.writeto('bad_pixel_map.fits')
