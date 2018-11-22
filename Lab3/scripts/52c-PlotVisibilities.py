# Analyze the sun and satellite data for Lab 3

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

# Load in the sun and satellite data from plots
# Key: file,BigMaxX,BigMaxY,NextMinX,NextMinY,NextMaxX,NextMaxY
sun = np.loadtxt('../maxmin-sun.txt', skiprows=3, delimiter=',')
sat = np.loadtxt('../maxmin-sat.txt', skiprows=3, delimiter=',')
types = ['Sun', 'Satellite']
files = [sun[:,0], sat[:,0]]
bigMaxX = [sun[:1], sat[:1]]
bigMaxY = [sun[:2], sat[:2]]
nextMinX = [sun[:3], sat[:3]]
nextMinY = [sun[:4], sat[:4]]
nextMaxX = [sun[:5], sat[:5]]
nextMaxY = [sun[:6], sat[:6]]

# Important constants
delta_az = 20 * np.pi / 180.     # change in azimuth in radians
wavelength = 2.7     # cm
cm2in = 1. / 2.54      # inches / cm
wavelength = wavelength * cm2in     # wavelength of light in inches

# Defines a process to flip the signal arrays to get maxima instead of minima
def centerflip(array):
    center = np.mean(array[-5:])
    for i in range(len(array)):
        array[i] = -1 * (array[i] - center)
    return(array)

# Calculate Visibilities. Index 0 = sun. Index 1 = satellite
visibilites = [[], []]
baselines = [[], []]
for slew in range(len(files)):
    for i in range(len(files[i])):
        #Calculate visibility
        visibilites[slew].append(
            (bigMaxY[slew][i] - nextMinY[slew][i]) /
            (bigMaxY[slew][i] + nextMinY[slew][i]) 
        )
        baselines[slew].append(1. / (bigMaxX[slew][i] + nextMinX[slew][i]) )
    plt.plot(baselines[slew], visibilities[slew])
    plt.xlabel(r'$B_{\lambda}$')
    plt.ylabel(r'Visibility, $V_0(B_{\lambda}$')
    plt.minorticks_on()
    plt.title('%s Interferometer Visibility' % types[slew])
    plt.savefig(info['images'] + 'visibilities-%s.pdf' % types[slew] , ppi=300)
    plt.clf()
