# AST443 Lab 0 Data Analysis for part 4.3
# September 7, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Open Files
flat_f_0 = fits.open('3.2_CCDFlats/3.2_flat_real.00000000.FIT')
flat_0 = flat_f_0[0].data
flat_f_1 = fits.open('3.2_CCDFlats/3.2_flat_real.00000001.FIT')
flat_1 = flat_f_1[0].data
flat_f_2 = fits.open('3.2_CCDFlats/3.2_flat_real.00000002.FIT')
flat_2 = flat_f_2[0].data
flat_f_3 = fits.open('3.2_CCDFlats/3.2_flat_real.00000003.FIT')
flat_3 = flat_f_3[0].data
flat_f_4 = fits.open('3.2_CCDFlats/3.2_flat_real.00000004.FIT')
flat_4 = flat_f_4[0].data
flat_f_5 = fits.open('3.2_CCDFlats/3.2_flat_real.00000005.FIT')
flat_5 = flat_f_5[0].data
flat_f_6 = fits.open('3.2_CCDFlats/3.2_flat_real.00000006.FIT')
flat_6 = flat_f_6[0].data
flat_f_7 = fits.open('3.2_CCDFlats/3.2_flat_real.00000007.FIT')
flat_7 = flat_f_7[0].data
flat_f_8 = fits.open('3.2_CCDFlats/3.2_flat_real.00000008.FIT')
flat_8 = flat_f_8[0].data
flat_f_9 = fits.open('3.2_CCDFlats/3.2_flat_real.00000009.FIT')
flat_9 = flat_f_9[0].data

# Average the flats and remove the bias
flat_all = [ flat_0, flat_1, flat_2, flat_3, flat_4, flat_5,
    flat_5, flat_6, flat_7, flat_8, flat_9]
flat_avg = np.mean( flat_all, axis=0 )
flat_f_bias = fits.open('neg10BIAS_master.fits')
flat_bias = flat_f_bias[0].data
flat_avgnorm = flat_avg - flat_bias

# Save the Master flat
master_write = fits.PrimaryHDU(flat_avgnorm)
#master_write.writeto('flat_master.fits')

# NOTE: The difference between the maximum and minimum points are 
#     1.59E4 - 1.18E4 = 0.41E4 = 4100

# Histogram of counts in the master flat field
flat_avgnorm_flat = flat_avgnorm.flatten()
flat_avgnorm_mean = np.mean(flat_avgnorm_flat)
flat_avgnorm_median = np.median(flat_avgnorm_flat)
flat_avgnorm_mode = stats.mode(flat_avgnorm_flat)[0][0]
flat_avgnorm_stddev = np.std(flat_avgnorm_flat)
fig, ax = plt.subplots()
flat_avgnorm_bins = 100
flat_avgnorm_min = min( flat_avgnorm_flat )
flat_avgnorm_max = max( flat_avgnorm_flat )
flat_avgnorm_norm = ( (flat_avgnorm_max-flat_avgnorm_min) 
    / flat_avgnorm_bins * len(flat_avgnorm_flat) )
xgauss = np.linspace( flat_avgnorm_min, flat_avgnorm_max, 
    10 * flat_avgnorm_bins )
ygauss = flat_avgnorm_norm * norm.pdf(xgauss,loc=flat_avgnorm_mean, scale=flat_avgnorm_stddev)
xmode = [flat_avgnorm_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f$' % (flat_avgnorm_mean, ),
    'Median={}'.format(flat_avgnorm_median),
    'Mode={}'.format(flat_avgnorm_mode) ,
    r'$\sigma=%.2f$' % (flat_avgnorm_stddev, )))
ax.hist(flat_avgnorm_flat, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Raw Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e6])
gauss = norm.pdf(xgauss,loc=flat_avgnorm_mean, scale=flat_avgnorm_stddev)
plt.plot(xgauss, ygauss, color="red", linewidth=1.0)
plt.plot(xmode, ymode, color="yellow", linewidth=1.0)
plt.savefig('flat_master_raw.pdf', ppi=300)
plt.clf()
