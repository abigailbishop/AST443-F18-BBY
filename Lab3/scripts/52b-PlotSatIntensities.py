# Analyze the satelite iterferometer observations for Lab 3

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
slew_dir = info['dataDir'] + info['raw2dish'] 
files = open(slew_dir + 'names.txt')
slews = []
fileNames = []
# for slew[i], inex 0 is the file name, 1 is time, 2 is current
for line in files:
    slews.append( [line.strip('\n'), np.loadtxt(slew_dir+line.strip('\n'), 
                    skiprows=1, delimiter=',') ] )
    fileNames.append(line.strip('\n'))
delta_az = 20 * np.pi / 180.     # change in azimuth in radians
wavelength = 2.7     # cm
cm2in = 1. / 2.54      # inches / cm
wavelength = wavelength * cm2in     # wavelength of light in inches

# Important data output file
output = open('../baselines-sat.txt', 'w')

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

# Initialize Figures
plt.figure(0)
plt.figure(1)
plt.figure(2)
plt.figure(3)
plt.figure(4)
plt.figure(5)
plt.figure(6)
plt.figure(7)
plt.figure(8)
plt.figure(9)
plt.figure(10)
plt.figure(11)
plt.figure(12)
plt.figure(13)
plt.figure(14)
plt.figure(15)
plt.figure(16)
plt.figure(17)
plt.figure(18)
plt.figure(19)
plt.figure(20)

# Loop over every slew across the object and save its plot of signal vs azimuth
for i in range(len(slews)):
    plt.figure(i)
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
    times = [j*np.cos(alt_rad) for j in times]
    # Center plot around maximum
    center = times[np.argmax(currents)]
    baseline_approx = wavelength / 2. / float(fileNames[i][24:26])
    # Plot 
    plt.plot(times, currents)
    plt.xlabel(r'$\Delta$ Azimuth (radians)')
    plt.ylabel('Current (A)')
    plt.xlim(center - 3*baseline_approx, center + 3*baseline_approx)
    plt.minorticks_on()
    plt.title('Interferometer - Sat - %.1f degrees Alt' % alt_deg)
    plt.savefig(info['sat2dishplots'] + slews[i][0][:-4] + '.pdf' , ppi=300)
    plt.clf()
