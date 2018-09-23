# Lab 0 Data Analysis for part 4.5
# Applies the normalized flat field to the arc lamp spectrum
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

# Open files
masterflat_f = fits.open('3.3_SpectrumFlats_all/master_DomeFlat_50E-6m.fits')
masterflat = masterflat_f[0].data
masterflat_np = np.asarray(masterflat)
masterflat_shape =  masterflat.shape
lampflat_f = fits.open('3.3_SpectrumFlat/masterflat_lamp.fits')
lampflat = lampflat_f[0].data
lampflat_np = np.asarray(lampflat)[175:223, :]
lampflat_shape =  lampflat.shape

# Obtain initial data arrays
pixel = []
master_avg = []
lamp_avg = []
for column in range(masterflat_shape[1]):
    # Looping over every column
    pixel.append(column)
    master_colavg = np.mean( masterflat[:, column].flatten() )
    master_avg.append( master_colavg )
    lamp_colavg = np.mean( lampflat[:, column].flatten())
    lamp_avg.append( lamp_colavg )
x = np.asarray(pixel)
y = np.asarray(master_avg)
ylamp = np.asarray(lamp_avg)

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
plt.savefig('spectrograph_crop_fit.pdf', ppi=300)
plt.clf()

# Normalize the Data
ymax = max(y)
y_norm_raw = []
ylamp_norm_raw = []
for i in range(len(x)):
    y_norm_raw.append( y[i] / yfit[i] )
    ylamp_norm_raw.append( ylamp[i] / yfit[i] )
y_norm = np.asarray(y_norm_raw)
ylamp_norm = np.asarray(ylamp_norm_raw)

# Plot the Normalized Data
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Normalized Average Value per Y Value')
ax.plot(x, y_norm, color="red", linewidth=2.0, label='Normalized Data')
plt.legend()
plt.savefig('spectrograph_cropnorm_fit.pdf', ppi=300)
plt.clf()

# Plot the adjusted arc spectrum
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Normalized Average Value per column')
ax.plot(x, ylamp_norm, color="red", linewidth=2.0, 
    label='Adjusted Arc Lamp Data')
plt.legend()
plt.savefig('spectrograph_lampadj_fit.pdf', ppi=300)
plt.clf()

# Plot the adjusted arc spectrum without normalization
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Normalized Average Value per column')
ax.plot(x, ylamp-y, color="red", linewidth=2.0, 
    label='Adjusted Arc Lamp Data')
plt.legend()
plt.savefig('spectrograph_lampadjnonorm_fit.pdf', ppi=300)
plt.clf()
