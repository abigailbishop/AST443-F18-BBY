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
from scipy import interpolate

# Open files
masterflat_f = fits.open('3.3_SpectrumFlat/masterflat_50E-6m_lamp.fits')
masterflat = masterflat_f[0].data
masterflat_np = np.asarray(masterflat)
masterflat_shape =  masterflat.shape

# Obtain arrays for plotting
pixel = []
avg_value = []
for column in range(masterflat_shape[1]):
    # Looping over every column
    pixel.append(column)
    avg_value.append( np.mean( masterflat[:, column].flatten() ))

# Plot the Avg. value as a function of pixel position
fig, ax = plt.subplots()
ax.set_title("Pixel Value as a function of Pixel Position")
ax.set_xlabel('Horizontal Pixel Index')
ax.set_ylabel('Average Value of Verticle Cells')
ax.plot(pixel, avg_value, color="red", linewidth=1.0, label='Data')
plt.legend()
plt.savefig('spectrograph_lampcrop_fit.pdf', ppi=300)
plt.clf()

