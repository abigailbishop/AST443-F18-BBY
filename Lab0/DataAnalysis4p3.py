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
flat_avg_mode = stats.mode(flat_avg.flatten())[0][0]
flat_avg_shape = flat_avg.shape
flat_avgnorm = np.zeros(flat_avg_shape)
for column in range(flat_avg_shape[1]):
    for row in range(flat_avg_shape[0]):
        flat_avgnorm[row][column] = flat_avg[row][column] / flat_avg_mode

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
    'Mean=%.2f' % (flat_avgnorm_mean, ),
    'Median=%.2f' % (flat_avgnorm_median, ),
    'Mode=%.2f' % (flat_avgnorm_mode, ) ,
    r'$\sigma=%.2f$' % (flat_avgnorm_stddev, )))
ax.hist(flat_avgnorm_flat, bins=100, color='black')
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
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

# Plot the relative sensitivity
sensitivity_0deg = np.loadtxt('sensitivity_data_0deg.txt', skiprows = 2)
absx_0deg = sensitivity_0deg[:,0]
absy_0deg = sensitivity_0deg[:,1]
absbright_0deg = sensitivity_0deg[:,2]
dist_0deg = []
relbright_0deg = []
sensitivity_90deg = np.loadtxt('sensitivity_data_90deg.txt', skiprows = 2)
absx_90deg = sensitivity_90deg[:,0]
absy_90deg = sensitivity_90deg[:,1]
absbright_90deg = sensitivity_90deg[:,2]
dist_90deg = []
relbright_90deg = []
for item in range( len( absx_0deg ) ):
    dist_0deg.append( np.sqrt( 
        ( absx_0deg[item] - absx_0deg[0] )**2 + 
        ( absy_0deg[item] - absy_0deg[0] )**2 
        ) )
    relbright_0deg.append( absbright_0deg[item] / absbright_0deg[0] )
    dist_90deg.append( np.sqrt( 
        ( absx_90deg[item] - absx_90deg[0] )**2 + 
        ( absy_90deg[item] - absy_90deg[0] )**2 
        ) )
    relbright_90deg.append( absbright_90deg[item] / absbright_90deg[0] )
fig, ax = plt.subplots()
ax.set_title("Relative Brightness")
ax.set_xlabel('Relative Distance from Center (pixels)')
ax.set_ylabel('Relative Brightness (sum/pix**2)')
plt.plot(dist_0deg, relbright_0deg, label='Rotation= 0 degrees')
plt.plot(dist_90deg, relbright_90deg, label='Rotation= 90 degrees')
plt.legend(loc='best')
plt.savefig('brightness-distance.pdf', ppi=300)
plt.clf()
# The relative brightness of a star will appear to decrease if it is
# closer to the edges of the CCD's field of view
# If you forgot to take flatfields on the night of observations, you 
# could not retake them later because the flat field depends on the orientation
# of the CCD in the telescope.

# Save the Master flat
master_write = fits.PrimaryHDU(flat_avgnorm)
#master_write.writeto('flat_master.fits')
