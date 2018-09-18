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
from scipy.optimize import curve_fit

# Open files
masterflat_f = fits.open('3.3_SpectrumFlats_all/master_DomeFlat_50E-6m.fits')
masterflat = masterflat_f[0].data
masterflat_np = np.asarray(masterflat)
masterflat_shape =  masterflat.shape

# Obtain initial data arrays
pixel = []
avg_value = []
for column in range(masterflat_shape[1]):
    # Looping over every column
    pixel.append(column)
    avg_value.append( np.mean( masterflat[:, column].flatten() ))
x = np.asarray(pixel)
y = np.asarray(avg_value)

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
    color = "blue", label = "curve fit")
plt.legend()
plt.savefig('spectrograph_crop_fit.pdf', ppi=300)
plt.clf()

# Normalize the Data
ymax = max(y)
y_norm_raw = []
yfit_norm_raw = []
for i in range(len(x)):
    y_norm_raw.append( y[i] / ymax )
    yfit_norm_raw.append( yfit[i] / ymax )
y_norm = np.asarray(y_norm_raw)
yfit_norm = np.asarray(yfit_norm_raw)

# Plot the Normalized Data
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Normalized Average Value per Y Value')
ax.plot(x, y_norm, color="red", linewidth=2.0, label='Data')
plt.plot(x, yfit_norm, linestyle = "--",
    color = "blue", label = "Curve Fit")
plt.legend()
plt.savefig('spectrograph_cropnorm_fit.pdf', ppi=300)
plt.clf()
