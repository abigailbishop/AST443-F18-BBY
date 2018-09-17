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
flat0_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000067.FIT')
flat0 = flat0_f[0].data
flat1_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000068.FIT')
flat1 = flat1_f[0].data
flat2_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000069.FIT')
flat2 = flat2_f[0].data
flat3_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000070.FIT')
flat3 = flat3_f[0].data
flat4_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000071.FIT')
flat4 = flat4_f[0].data
flat5_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000072.FIT')
flat5 = flat5_f[0].data
flat6_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000073.FIT')
flat6 = flat6_f[0].data
flat7_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000074.FIT')
flat7 = flat7_f[0].data
flat8_f = fits.open('3.3_SpectrumFlat/3.3_spectrumFlats.00000075.FIT')
flat8 = flat8_f[0].data

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
