# Lab 0 Data Analysis for part 4.2
# Creates and analyzes DARK frames set at -10C
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
dark_f_0 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000017.DARK.FIT')
dark_0 = dark_f_0[0].data
dark_f_1 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000018.DARK.FIT')
dark_1 = dark_f_1[0].data
dark_f_2 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000019.DARK.FIT')
dark_2 = dark_f_2[0].data
dark_f_3 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000020.DARK.FIT')
dark_3 = dark_f_3[0].data
dark_f_4 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000021.DARK.FIT')
dark_4 = dark_f_4[0].data
dark_f_5 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000022.DARK.FIT')
dark_5 = dark_f_5[0].data
dark_f_6 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000023.DARK.FIT')
dark_6 = dark_f_6[0].data
dark_f_7 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000024.DARK.FIT')
dark_7 = dark_f_7[0].data
dark_f_8 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000025.DARK.FIT')
dark_8 = dark_f_8[0].data
dark_f_9 = fits.open('3.1_m10_DARK_1minexpo/3.1_m10.00000026.DARK.FIT')
dark_9 = dark_f_9[0].data

# Create the Dark Master Frame
dark_master_array = [dark_0, dark_1, dark_2, dark_3, dark_4, dark_5,
    dark_6, dark_7, dark_8, dark_9 ]
dark_master = np.median(dark_master_array, axis=0)

# Unclipped Data Statistical Properties
dark_master_flat = dark_master.flatten()
dark_master_flat = dark_master_flat[dark_master_flat<15000]
dark_mean = np.mean(dark_master_flat)
dark_median = np.median(dark_master_flat)
dark_mode = stats.mode(dark_master_flat)[0][0]
dark_stddev = np.std(dark_master_flat)
fig, ax = plt.subplots()
dark_bins = 100
dark_min = min( dark_master_flat )
dark_max = max( dark_master_flat )
dark_norm = (dark_max-dark_min) / dark_bins * len(dark_master_flat)
xgauss = np.linspace( dark_min, dark_max, 10 * dark_bins )
ygauss = dark_norm * norm.pdf(xgauss,loc=dark_mean, scale=dark_stddev)
xmode = [dark_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f' % (dark_mean, ),
    'Median={}'.format(dark_median),
    'Mode={}'.format(dark_mode) ,
    r'$\sigma=%.2f$' % (dark_stddev, )))
ax.hist(dark_master_flat, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Raw Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e6])
gauss = norm.pdf(xgauss,loc=dark_mean, scale=dark_stddev)
plt.plot(xgauss, ygauss, color="red", linewidth=1.0)
plt.plot(xmode, ymode, color="yellow", linewidth=1.0)
plt.savefig('neg10DARK_noOutliers.pdf', ppi=300)
plt.clf()

# Clipped Data Statistical Properties
fig, ax = plt.subplots()
dark_master_sigclip = stats.sigmaclip(dark_master_flat, 5.0, 5.0)
dark_master_cut = dark_master_sigclip[0]
dark_min = dark_master_sigclip[1]
dark_max = dark_master_sigclip[2]
dark_mean = np.mean(dark_master_cut)
dark_median = np.median(dark_master_cut)
dark_mode = stats.mode(dark_master_cut)[0][0]
dark_stddev = np.std(dark_master_cut)
fig, ax = plt.subplots()
dark_bins = 100
dark_norm = (dark_max-dark_min) / dark_bins * len(dark_master_cut)
xgauss = np.linspace( dark_min, dark_max, 10 * dark_bins )
ygauss = dark_norm * norm.pdf(xgauss,loc=dark_mean, scale=dark_stddev)
xmode = [dark_mode] * 100
ymode = np.linspace( 0, max(ygauss), len(xmode) )
textstr = '\n'.join((
    'Mean=%.2f' % (dark_mean, ),
    'Median={}'.format(dark_median),
    'Mode={}'.format(dark_mode) ,
    r'$\sigma=%.2f$' % (dark_stddev, )))
ax.hist(dark_master_cut, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Sigma Clipped Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e6])
gauss = norm.pdf(xgauss,loc=dark_mean, scale=dark_stddev)
plt.plot(xgauss, ygauss, color="red", linewidth=1.0)
plt.plot(xmode, ymode, color="yellow", linewidth=1.0)
plt.savefig('neg10DARK_cut.pdf', ppi=300)
plt.clf()
print( 'Percent Rejected = {}'.format( 100 * 
    ( len(dark_master_flat.flatten()) - len(dark_master_cut) ) / 
      len(dark_master_flat.flatten()) ) )

# Identifying Hot Pixels
# Hot pixels have high counts, above the BIAS. 
# Counts beyond 5\sigma of the standard distribution are hot
# I want to subtract the master BIAS from the master DARK and then analyze
bias_master = fits.open('neg10BIAS_master.fits')
bias_master = bias_master[0].data
dark_adjusted = dark_master - bias_master
dark_adjusted_flat = dark_adjusted.flatten()
dark_adjusted_flat_sigclip = stats.sigmaclip( dark_adjusted_flat, 5.0, 5.0)
dark_adjusted_clip = dark_adjusted_flat_sigclip[0]
dark_adjusted_min = dark_adjusted_flat_sigclip[1]
dark_adjusted_max = dark_adjusted_flat_sigclip[2]
dark_adjusted_mean = np.mean(dark_adjusted_flat)
dark_adjusted_median = np.median(dark_adjusted_flat)
dark_adjusted_mode = stats.mode(dark_adjusted_flat)[0][0]
dark_adjusted_stddev = np.std(dark_adjusted_flat)
fig, ax = plt.subplots()
ax.set_title("Adjusted Data")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e7])
textstr = '\n'.join((
    'Mean=%.2f' % (dark_adjusted_mean, ),
    'Median=%.2f' % (dark_adjusted_median, ),
    'Mode=%.2f' % (dark_adjusted_mode, ) ,
    r'$\sigma=%.2f$' % (dark_adjusted_stddev, )))
ax.hist(dark_adjusted_flat, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.65, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
plt.savefig('neg10DARK_adjusted.pdf', ppi=300)
hot_column = []
hot_row = []
dark_5sigma = (dark_max - dark_min) / 2.0
for column in range( dark_adjusted.shape[0] ) : 
    for row in range( dark_adjusted.shape[1] ) :
        if ( dark_adjusted[column][row] > dark_adjusted_max ):
            hot_column.append(column)
            hot_row.append(row)
print( 'Percent Hot = {}'.format( 100 * 
    ( len(hot_column)  / ( len(dark_adjusted[0])*len(dark_adjusted[1]) ) ) ) )

# Save the master DARK minus the BIAS
master_write = fits.PrimaryHDU(dark_master)
#master_write.writeto('neg10dark_adjusted.fits')

# Calculating the Dark Current
dark_adjusted_array = []
for dark in dark_master_array:
    dark_adjusted_array.append(dark - bias_master)
dark_adjusted_all_flat = np.asarray(dark_adjusted_array)
dark_adjusted_all_flat = dark_adjusted_all_flat.flatten()
fig, ax = plt.subplots()
dark_adjusted_all_sigclip = stats.sigmaclip(dark_adjusted_all_flat, 5.0, 5.0)
dark_adjusted_all_cut = dark_adjusted_all_sigclip[0]
dark_adjusted_all_min = dark_adjusted_all_sigclip[1]
dark_adjusted_all_max = dark_adjusted_all_sigclip[2]
dark_adjusted_all_mean = np.mean(dark_adjusted_all_cut)
dark_adjusted_all_median = np.median(dark_adjusted_all_cut)
dark_adjusted_all_mode = stats.mode(dark_adjusted_all_cut)[0][0]
dark_adjusted_all_stddev = np.std(dark_adjusted_all_cut)
dark_adjusted_all_uncert = np.sqrt( dark_adjusted_all_mean )
fig, ax = plt.subplots()
textstr = '\n'.join((
    'Mean=%.2f' % (dark_adjusted_all_mean, ),
    'Median=%.2f' % (dark_adjusted_all_median),
    'Mode=%.2f' % (dark_adjusted_all_mode) ,
    r'$\sigma=%.2f$' % (dark_adjusted_all_stddev, ),
    'Uncertainty=%.2f' % (dark_adjusted_all_uncert, ) ))
ax.hist(dark_adjusted_all_flat, range=[0, 5000], bins=500, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.6, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
ax.set_title("Sigma Adjusted Dark Current")
ax.set_yscale("log", nonposy='clip')
ax.set_xlabel('Number of Counts in a Bin')
ax.set_ylabel('Number of Bins')
ax.set_ylim([0.1,1e7])
gauss = norm.pdf(xgauss,loc=dark_mean, scale=dark_stddev)
plt.savefig('neg10DARK_darkcurrent.pdf', ppi=300)
plt.clf()

# Typical Counts as a function of Exposure time
dark_f_t0 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000010.DARK.FIT')
dark_t0 = (dark_f_t0[0].data - bias_master).flatten()
dark_t0_exp = dark_f_t0[0].header['EXPTIME']
dark_f_t1 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000011.DARK.FIT')
dark_t1 = (dark_f_t1[0].data - bias_master).flatten()
dark_t1_exp = dark_f_t1[0].header['EXPTIME']
dark_f_t2 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000012.DARK.FIT')
dark_t2 = (dark_f_t2[0].data - bias_master).flatten()
dark_t2_exp = dark_f_t2[0].header['EXPTIME']
dark_f_t3 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000013.DARK.FIT')
dark_t3 = (dark_f_t3[0].data - bias_master).flatten()
dark_t3_exp = dark_f_t3[0].header['EXPTIME']
dark_f_t4 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000014.DARK.FIT')
dark_t4 = (dark_f_t4[0].data - bias_master).flatten()
dark_t4_exp = dark_f_t4[0].header['EXPTIME']
dark_f_t5 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000015.DARK.FIT')
dark_t5 = (dark_f_t5[0].data - bias_master).flatten()
dark_t5_exp = dark_f_t5[0].header['EXPTIME']
dark_f_t6 = fits.open('3.1_m10_DARK_varyexpo/3.1_m10.00000016.DARK.FIT')
dark_t6 = (dark_f_t6[0].data - bias_master).flatten()
dark_t6_exp = dark_f_t6[0].header['EXPTIME']
dark_modes = [ 
    stats.mode(dark_t0)[0][0], 
    stats.mode(dark_t1)[0][0], 
    stats.mode(dark_t2)[0][0], 
    stats.mode(dark_t3)[0][0], 
    stats.mode(dark_t4)[0][0], 
    stats.mode(dark_t5)[0][0], 
    stats.mode(dark_t6)[0][0] 
    ]
dark_means = [
    np.mean(dark_t0),
    np.mean(dark_t1),
    np.mean(dark_t2),
    np.mean(dark_t3),
    np.mean(dark_t4),
    np.mean(dark_t5),
    np.mean(dark_t6)
    ]
dark_modes_uncerts = []
for mean in dark_means: 
    dark_modes_uncerts.append( np.sqrt( abs(mean) ) ) 
dark_exps = [
    dark_t0_exp,
    dark_t1_exp,
    dark_t2_exp,
    dark_t3_exp,
    dark_t4_exp,
    dark_t5_exp,
    dark_t6_exp
    ]
fig, ax = plt.subplots()
ax.set_title('"Typical" Counts as a function of Exposure Time')
ax.set_xlabel('Exposure Time (seconds)')
ax.set_ylabel('Number of "Typical" Counts (e-/pixel)')
ax.errorbar( dark_exps, dark_modes, yerr=dark_modes_uncerts, fmt='o')
darks_linearfit = np.polynomial.polynomial.polyfit(
    dark_exps, dark_modes, deg=1 )
def linear_function( x, b, m ):
    return m * x + b
darks_linregression = linear_function(np.asarray(dark_exps), *darks_linearfit)
ax.plot(dark_exps, darks_linregression, linestyle = '--')
dark_darkcurrent = ( ( darks_linregression[1] - darks_linregression[0] ) /
    ( dark_exps[1] - dark_exps[0] ) )
textstr = 'Dark Current=%.2f e-/px/sec' % (dark_darkcurrent, )
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.06, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
plt.savefig('neg10DARK_typical-exposure.pdf', ppi=300)
plt.clf()
