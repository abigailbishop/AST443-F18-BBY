# Lab 0 Data Analysis for part 4.1
# Creates and analyzes BIAS frames set at +10C
# September 7, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

pos10BIAS_f = fits.open('3.1_p10_BIAS/3.1_p10.00000000.BIAS.FIT')
pos10BIAS_header = pos10BIAS_f[0].header
pos10BIAS_gain = pos10BIAS_header['EGAIN']
print( 'Gain = ', pos10BIAS_gain )
pos10BIAS_d = pos10BIAS_f[0].data
pos10BIAS_0 = pos10BIAS_d.flatten()

# Unclipped Data
fig, ax = plt.subplots()
data_mean = np.mean(pos10BIAS_0)
data_median = np.median(pos10BIAS_0)
data_mode = float( stats.mode(pos10BIAS_0)[0][0] )
data_stddev = np.std(pos10BIAS_0)
data_readnoise = data_stddev * pos10BIAS_gain
print( 'Unclipped Read Noise: ', data_readnoise)
data_bins = 100
data_min = min( pos10BIAS_0 )
data_max = max( pos10BIAS_0 )
data_norm = (data_max-data_min) / data_bins * len(pos10BIAS_0)
xgauss = np.linspace( data_min, data_max, 10 * data_bins )
ygauss = data_norm * norm.pdf(xgauss,loc=data_mean, scale=data_stddev)
xmode = [data_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f' % (data_mean, ),
    'Median={}'.format(data_median),
    'Mode={}'.format(data_mode) ,
    r'$\sigma=%.2f$' % (data_stddev, )))
ax.hist(pos10BIAS_0, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Raw Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e6])
gauss = norm.pdf(xgauss,loc=data_mean, scale=data_stddev)
plt.plot(xgauss, ygauss, color="red", linewidth=1.0)
plt.plot(xmode, ymode, color="yellow", linewidth=1.0)
plt.savefig('pos10BIAS_raw.pdf', ppi=300)
plt.clf()

# Clipped Data
fig, ax = plt.subplots()
pos10BIAS_c = pos10BIAS_0[pos10BIAS_0>960]
pos10BIAS_c = pos10BIAS_c[pos10BIAS_c<1040]
clip_mean = np.mean(pos10BIAS_c)
clip_median = np.median(pos10BIAS_c)
clip_mode = float( stats.mode(pos10BIAS_c)[0][0] )
clip_stddev = np.std(pos10BIAS_c)
clip_readnoise = clip_stddev * pos10BIAS_gain
print( 'Clipped Read Noise: ', clip_readnoise)
clip_bins = 100
clip_min = min( pos10BIAS_c )
clip_max = max( pos10BIAS_c )
clip_norm = (clip_max-clip_min) / clip_bins * len(pos10BIAS_c)
xgauss = np.linspace( clip_min, clip_max, 10 * clip_bins )
ygauss = clip_norm * norm.pdf(xgauss,loc=clip_mean, scale=clip_stddev)
xmode = [clip_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f' % (clip_mean, ),
    'Median={}'.format(clip_median),
    'Mode={}'.format(clip_mode) ,
    r'$\sigma=%.2f$' % (clip_stddev, )))
ax.hist(pos10BIAS_c, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Raw Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e6])
gauss = norm.pdf(xgauss,loc=data_mean, scale=data_stddev)
plt.plot(xgauss, ygauss, color="red", linewidth=3.0)
plt.plot(xmode, ymode, color="yellow", linewidth=3.0)
plt.savefig('pos10BIAS_clipped.pdf', ppi=300)
plt.clf()

print( 'Percent Rejected = {}'.format( 100 *
    ( len(pos10BIAS_0) - len(pos10BIAS_c) ) / len (pos10BIAS_0) ) )

# Creating a master BIAS Frame
pos10BIAS_1= fits.open('3.1_p10_BIAS/3.1_p10.00000001.BIAS.FIT')[0].data
pos10BIAS_2= fits.open('3.1_p10_BIAS/3.1_p10.00000002.BIAS.FIT')[0].data
pos10BIAS_3= fits.open('3.1_p10_BIAS/3.1_p10.00000003.BIAS.FIT')[0].data
pos10BIAS_4= fits.open('3.1_p10_BIAS/3.1_p10.00000004.BIAS.FIT')[0].data
pos10BIAS_5= fits.open('3.1_p10_BIAS/3.1_p10.00000005.BIAS.FIT')[0].data
pos10BIAS_6= fits.open('3.1_p10_BIAS/3.1_p10.00000006.BIAS.FIT')[0].data
pos10BIAS_7= fits.open('3.1_p10_BIAS/3.1_p10.00000007.BIAS.FIT')[0].data
pos10BIAS_8= fits.open('3.1_p10_BIAS/3.1_p10.00000008.BIAS.FIT')[0].data
pos10BIAS_9= fits.open('3.1_p10_BIAS/3.1_p10.00000009.BIAS.FIT')[0].data
pos10BIAS_array = [ pos10BIAS_d, pos10BIAS_1, pos10BIAS_2, pos10BIAS_3, 
    pos10BIAS_4, pos10BIAS_5, pos10BIAS_6,
    pos10BIAS_7, pos10BIAS_8, pos10BIAS_9 ]
pos10BIAS_master = np.mean(pos10BIAS_array, axis=0)
pos10BIAS_flatmaster = pos10BIAS_master.flatten()
#pos10BIAS_sum = pos10BIAS_sum[ pos10BIAS_sum >clip_min]
#pos10BIAS_sum = pos10BIAS_sum[pos10BIAS_sum<clip_max]
master_mean = np.mean(pos10BIAS_flatmaster)
print( 'Master mean = ', master_mean)
master_stddev = np.std(pos10BIAS_flatmaster)
print( 'Master Standard Deviation = ', master_stddev)
print( 'STD Decreased by: ', ( 1 - master_stddev / data_stddev ) )
print( '1/sqrt{N_images} = ', 1 / np.sqrt( len(pos10BIAS_array) ) )

# Save the master BIAS
master_write = fits.PrimaryHDU(pos10BIAS_master)
#master_write.writeto('pos10BIAS_master.fits')
