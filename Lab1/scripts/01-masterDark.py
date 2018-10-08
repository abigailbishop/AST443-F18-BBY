# Data Analysis for part 4.1
# Creates the Master Dark Frame
# October 8, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Load Lab constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Open Files
dark_dir = info['dataDir'] + info['darkSubdir']
files = open(dark_dir+'names.txt', 'r')
darks = []
for line in files:
    darks.append(fits.open(dark_dir+line.strip('\n')))
darks_data = []
for dark in darks:
    darks_data.append(dark[0].data)

# Create the Dark Master Frame
dark_master = np.median(darks_data, axis=0)
master_write = fits.PrimaryHDU(dark_master)
master_write.writeto(info['fitsFiles']+info['masterDark'])

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
plt.savefig(info['images']+'masterDark-dist.pdf', ppi=300)
plt.clf()
