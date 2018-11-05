# Lab 2 normalizing the flat fields. Part 7.2.a
# Applies the normalized flat field to the arc star spectrum
# September 10, 2018

# Imports
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm
from scipy.optimize import curve_fit

# Load constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Open files
masterflat_f = fits.open(info['fitsFiles'] + 'masterFlat.fits')
masterflat = masterflat_f[0].data
masterflat_np = np.asarray(masterflat)
masterflat_shape =  masterflat.shape
starflat_f = fits.open(info['dataDir'] + info['dccRefstarData'] +
     'data_.00000008_dcc.fits')
starflat = starflat_f[0].data
starflat_np = np.asarray(starflat)[46:94, :]
starflat_shape =  starflat.shape

# Obtain initial data arrays
pixel = []
master_avg = []
star_avg = []
for column in range(masterflat_shape[1]):
    # Looping over every column
    pixel.append(column)
    master_colavg = np.mean( masterflat[:, column].flatten() )
    master_avg.append( master_colavg )
    star_colavg = np.mean( starflat[:, column].flatten())
    star_avg.append( star_colavg )
x = np.asarray(pixel)
y = np.asarray(master_avg)
ystar = np.asarray(star_avg)

# Define the funtion to model our data after
def unknown_function(x,a,b,c):
    return a*x**2 + b*x + c

# Fit the Data
popt, pcov = curve_fit(unknown_function, x, y)
yfit = unknown_function(x,*popt)

# Plot the Avg. value as a function of pixel position
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Average Value')
ax.plot(x, y, color="red", linewidth=2.0, label='Data')
plt.plot(x, yfit, linestyle = "--", 
    color = "blue", label = "Quadratic Fit")
plt.legend()
plt.savefig(info['images'] + 'spectrograph_crop_fit.pdf', ppi=300)
plt.clf()

# Normalize the Data
ymax = max(y)
y_norm_raw = []
ystar_norm_raw = []
for i in range(len(x)):
    y_norm_raw.append( y[i] / yfit[i] )
    ystar_norm_raw.append( ystar[i] / yfit[i] )
y_norm = np.asarray(y_norm_raw)
ystar_norm = np.asarray(ystar_norm_raw)

# Plot the Normalized Data
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Normalized Average Value per Y Value')
ax.plot(x, y_norm, color="red", linewidth=2.0, label='Normalized Data')
plt.legend()
plt.savefig(info['images'] + 'spectrograph_cropnorm_fit.pdf', ppi=300)
plt.clf()

# Plot the adjusted arc spectrum
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Normalized Average Value per column')
ax.plot(x, ystar_norm, color="red", linewidth=2.0, 
    label='Adjusted Arc Lamp Data')
plt.legend()
plt.savefig(info['images'] + 'spectrograph_staradj_fit.pdf', ppi=300)
plt.clf()

