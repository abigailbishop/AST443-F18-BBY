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

# Important data output file
output = open('../baselines.txt', 'w')

# Flip and center the currents ALSO plot the curves
def centerflip(array):
    center = np.mean(array[-5:])
    for i in range(len(array)):
        array[i] = -1 * (array[i] - center)
    return(array)
for i in range(len(slews)):
    print( 'Analyzing ', slews[i][0])
    times = slews[i][1][:,0]
    currents = slews[i][1][:,1]
    currents = centerflip(currents)
    alt_deg = float(slews[i][0][31:33])
    alt_rad = alt_deg * np.pi / 180.
    times = [j*np.cos(alt_rad) for j in times]
    slewrate = delta_az / times[-1]
    times = [j*slewrate for j in times]
    plt.plot(times, currents)
    plt.xlabel(r'$\Delta$ Azimuth (radians)')
    plt.ylabel('Current (A)')
    plt.title('Interferometer - Sun - %.1f degrees Alt' % alt_deg)
    plt.savefig(info['sun2dishplots'] + slews[i][0][:-4] , ppi=300)
    plt.clf()
    # Get the big points
    bigMaxIdx = np.argmax(currents) 
    rightMinIdx = 0
    for j in range(bigMaxIdx+1, len(currents)):
        if currents[j] <= currents[j-1]:
            rightMinIdx = j
        else:
            break
    rightMaxIdx = np.argmax(currents[rightMinIdx:]) + rightMinIdx
    rightMinIdx = np.argmin(currents[bigMaxIdx:rightMaxIdx]) + bigMaxIdx
    deltaAngle = (times[rightMaxIdx] - times[bigMaxIdx])
    baseline = wavelength / deltaAngle      # cm 
    output.write('%s Baseline = %.6f\n' %  (slews[i][0], baseline) )
