# Analyze the single dish observations for Lab 3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from astropy.stats import gaussian_sigma_to_fwhm
from operator import itemgetter

# Load Constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Import the Sun data
slew_dir = info['dataDir'] + info['raw2dish'] + 'Sun/'
files = open(slew_dir + 'names.txt')
slews = []
# for slew[i], inex 0 is the file name, 1 is time, 2 is current
for line in files:
    slews.append( [line.strip('\n'), np.loadtxt(slew_dir+line.strip('\n'), 
                    skiprows=1, delimiter=',') ] )
delta_az = 20 * np.pi / 180.     # change in azimuth in radians
wavelength = 2.7     # cm
cm2in = 1. / 2.54      # inches / cm
wavelength = wavelength * cm2in     # wavelength of light in inches

# Important data output file
output = open('../baselines.txt', 'w')

# Defines a process to flip the signal arrays to get maxima instead of minima
def centerflip(array):
    center = np.mean(array[-5:])
    for i in range(len(array)):
        array[i] = -1 * (array[i] - center)
    return(array)

# Defines a process that finds the index and value of an array closest to a 
#     provided number
def nearest(array, value):
    array = np.asarray(array)
    i = (np.abs(array - value)).argmin()
    return [array[i], i]

# Loop over every slew across the object and save its plot of signal vs azimuth
for i in range(len(slews)):
    # So this is the analysis for a singular slew across an object
    print( 'Analyzing ', slews[i][0])
    # Load Data
    times = slews[i][1][:,0]
    currents = slews[i][1][:,1]
    currents = centerflip(currents)
    alt_deg = float(slews[i][0][31:33])
    alt_rad = alt_deg * np.pi / 180.
    # Convert times to azimuthal angles adjusted for altitude in sky
    slewrate = delta_az / times[-1]
    times = [j*slewrate for j in times]
    #times = [j*np.cos(alt_rad) for j in times]
    # Plot 
    plt.plot(times, currents)
    plt.xlabel(r'$\Delta$ Azimuth (radians)')
    plt.ylabel('Current (A)')
    plt.minorticks_on()
    plt.title('Interferometer - Sun - %.1f degrees Alt' % alt_deg)
    plt.savefig(info['sun2dishplots'] + slews[i][0][:-4] + '.pdf' , ppi=300)
    plt.clf()
