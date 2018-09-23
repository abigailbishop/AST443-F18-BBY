# Lab 0 Data Analysis for part 4.1
# Creates and analyzes BIAS frames set at -10C
# September 7, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

neg10BIAS_f = fits.open('3.1_m10_BIAS/3.1_m10.00000000.BIAS.FIT')
neg10BIAS_header = neg10BIAS_f[0].header
neg10BIAS_gain = neg10BIAS_header['EGAIN']
print( 'Gain = ', neg10BIAS_gain )
neg10BIAS_d = neg10BIAS_f[0].data
neg10BIAS_0 = neg10BIAS_d.flatten()

# Unclipped Data
fig, ax = plt.subplots()
data_mean = np.mean(neg10BIAS_0)
data_median = np.median(neg10BIAS_0)
data_mode = float( stats.mode(neg10BIAS_0)[0][0] )
data_stddev = np.std(neg10BIAS_0)
data_readnoise = data_stddev * neg10BIAS_gain
print( 'Unclipped Read Noise: ', data_readnoise)
data_bins = 100
data_min = min( neg10BIAS_0 )
data_max = max( neg10BIAS_0 )
data_norm = (data_max-data_min) / data_bins * len(neg10BIAS_0)
xgauss = np.linspace( data_min, data_max, 10 * data_bins )
ygauss = data_norm * norm.pdf(xgauss,loc=data_mean, scale=data_stddev)
xmode = [data_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f' % (data_mean, ),
    'Median={}'.format(data_median),
    'Mode={}'.format(data_mode) ,
    r'$\sigma=%.2f$' % (data_stddev, )))
ax.hist(neg10BIAS_0, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Raw Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,3e6])
gauss = norm.pdf(xgauss,loc=data_mean, scale=data_stddev)
plt.plot(xmode, ymode, color="yellow", linewidth=1.0)
plt.plot(xgauss, ygauss, color="red", linewidth=1.0)
plt.savefig('neg10BIAS_raw.pdf', ppi=300)
plt.clf()

# Clipped Data
fig, ax = plt.subplots()
neg10BIAS_c = neg10BIAS_0[neg10BIAS_0>960]
neg10BIAS_c = neg10BIAS_c[neg10BIAS_c<1040]
clip_mean = np.mean(neg10BIAS_c)
clip_median = np.median(neg10BIAS_c)
clip_mode = float( stats.mode(neg10BIAS_c)[0][0] )
clip_stddev = np.std(neg10BIAS_c)
clip_readnoise = clip_stddev * neg10BIAS_gain
print( 'Clipped Read Noise: ', clip_readnoise)
clip_bins = 100
clip_min = min( neg10BIAS_c )
clip_max = max( neg10BIAS_c )
clip_norm = (clip_max-clip_min) / clip_bins * len(neg10BIAS_c)
xgauss = np.linspace( clip_min, clip_max, 10 * clip_bins )
ygauss = clip_norm * norm.pdf(xgauss,loc=clip_mean, scale=clip_stddev)
xmode = [clip_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f' % (clip_mean, ),
    'Median={}'.format(clip_median),
    'Mode={}'.format(clip_mode) ,
    r'$\sigma=%.2f$' % (clip_stddev, )))
ax.hist(neg10BIAS_c, bins=100, color='black')
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
plt.savefig('neg10BIAS_clipped.pdf', ppi=300)
plt.clf()

print( 'Percent Rejected = {}'.format( 100 * 
    ( len(neg10BIAS_0) - len(neg10BIAS_c) ) / len (neg10BIAS_0) ) )

# Creating a master BIAS Frame
neg10BIAS_1= fits.open('3.1_m10_BIAS/3.1_m10.00000001.BIAS.FIT')[0].data
neg10BIAS_2= fits.open('3.1_m10_BIAS/3.1_m10.00000002.BIAS.FIT')[0].data
neg10BIAS_3= fits.open('3.1_m10_BIAS/3.1_m10.00000003.BIAS.FIT')[0].data
neg10BIAS_4= fits.open('3.1_m10_BIAS/3.1_m10.00000004.BIAS.FIT')[0].data
neg10BIAS_5= fits.open('3.1_m10_BIAS/3.1_m10.00000005.BIAS.FIT')[0].data
neg10BIAS_6= fits.open('3.1_m10_BIAS/3.1_m10.00000006.BIAS.FIT')[0].data
neg10BIAS_7= fits.open('3.1_m10_BIAS/3.1_m10.00000007.BIAS.FIT')[0].data
neg10BIAS_8= fits.open('3.1_m10_BIAS/3.1_m10.00000008.BIAS.FIT')[0].data
neg10BIAS_9= fits.open('3.1_m10_BIAS/3.1_m10.00000009.BIAS.FIT')[0].data
neg10BIAS_array = [ neg10BIAS_d, neg10BIAS_1, neg10BIAS_2, neg10BIAS_3, 
    neg10BIAS_4, neg10BIAS_5, neg10BIAS_6,
    neg10BIAS_7, neg10BIAS_8, neg10BIAS_9 ]
neg10BIAS_master = np.mean(neg10BIAS_array, axis=0)
neg10BIAS_flatmaster = neg10BIAS_master.flatten()
#neg10BIAS_sum = neg10BIAS_sum[ neg10BIAS_sum >clip_min]
#neg10BIAS_sum = neg10BIAS_sum[neg10BIAS_sum<clip_max]
master_mean = np.mean(neg10BIAS_flatmaster)
print( 'Master mean = ', master_mean)
master_stddev = np.std(neg10BIAS_flatmaster)
print( 'Master Standard Deviation = ', master_stddev)
print( 'STD Decreased by: ', ( 1 - master_stddev / data_stddev ) )
print( '1/sqrt{N_images} = ', 1 / np.sqrt( len(neg10BIAS_array) ) )

# Save the master BIAS
master_write = fits.PrimaryHDU(neg10BIAS_master)
#master_write.writeto('neg10BIAS_master.fits')
