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
sun = [np.loadtxt('../maxmin-sun.txt', skiprows=3, 
                  delimiter=',', usecols=[0], dtype=str),
       np.loadtxt('../maxmin-sun.txt', skiprows=3, 
                  delimiter=',', usecols=range(1,9)) ]
sat = [np.loadtxt('../maxmin-sat.txt', skiprows=3, 
                  delimiter=',', usecols=[0], dtype=str),
       np.loadtxt('../maxmin-sat.txt', skiprows=3, 
                  delimiter=',', usecols=range(1,9)) ]
types = ['Sun', 'Satellite']
files = [sun[0], sat[0]]
bigMaxX = [sun[1][:,0], sat[1][:,0]]
bigMaxY = [sun[1][:,1], sat[1][:,1]]
nextMinX = [sun[1][:,2], sat[1][:,2]]
nextMinY = [sun[1][:,3], sat[1][:,3]]
nextMaxX = [sun[1][:,4], sat[1][:,4]]
nextMaxY = [sun[1][:,5], sat[1][:,5]]
errorX = [sun[1][:,6], sat[1][:,6]]
errorY = [sun[1][:,7], sat[1][:,7]]

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
visibilities = [[], []]
baselines = [[], []]
baselines_exp = [[], []]
for slew in range(len(files)):
    print('Analyzing %s data' % types[slew])
    errorV = []
    errorB = []
    for i in range(len(files[slew])):
        maxX = bigMaxX[slew][i]
        maxY = bigMaxY[slew][i]
        minX = nextMinX[slew][i]
        minY = nextMinY[slew][i]
        errX = errorX[slew][i]
        errY = errorY[slew][i]
        visibilities[slew].append( (maxY - minY) / (maxY + minY) )
        errorV.append( np.sqrt( (2 * errY**2) * ( (1. / (maxY + minY)**2) + 
                     ( (maxY - minY)**2 / (maxY + minY)**4)) ) )
        baselines[slew].append(.5 / (maxX + minX) )
        baselines_exp[slew].append(float(files[slew][i][24:26]))
        errorB.append( np.sqrt((2 * errX**2) / (minX - maxX)**2 ) ) 
    plt.errorbar(baselines[slew], visibilities[slew], 
    #plt.errorbar(baselines_exp[slew], visibilities[slew], 
    #       xerr = errorB, yerr = errorV,
    #       xerr = [0.5]*len(errorV), yerr = errorV,
    fmt = '.')
    plt.xlabel(r'$B_{\lambda}$')
    plt.ylabel(r'Visibility, $V_0(B_{\lambda}$)')
    plt.minorticks_on()
    plt.title('%s Interferometer Visibility' % types[slew])
    plt.savefig(info['images'] + 'visibilities-%s.pdf' % types[slew] , ppi=300)
    plt.clf()
