# Analyze the sun and satellite data for Lab 3

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
from scipy.stats import norm
from scipy import asarray as ar,exp
from scipy.optimize import curve_fit
from astropy.stats import gaussian_sigma_to_fwhm
from operator import itemgetter

AU2km = 1.496e+8 # NASA , km
alpha_litVal = 1919.0*(1/3600.0)*(np.pi/180) # NASA, 1919 seconds of arc
d_litVal = alpha_litVal * AU2km

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
                  delimiter=',', usecols=range(1,5)) ]
sat = [np.loadtxt('../maxmin-sat.txt', skiprows=3, 
                  delimiter=',', usecols=[0], dtype=str),
       np.loadtxt('../maxmin-sat.txt', skiprows=3, 
                  delimiter=',', usecols=range(1,5)) ]
types = ['Sun', 'Satellite']
files = [sun[0], sat[0]]
bigMaxX = [sun[1][:,0], sat[1][:,0]]
bigMaxY = [sun[1][:,1], sat[1][:,1]]
nextMinY = [sun[1][:,2], sat[1][:,2]]
nextMaxX = [sun[1][:,3], sat[1][:,3]]

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

# Define a sinc function
def sincFunc(BLambda, alpha):
    return abs( np.sinc(BLambda * alpha) )

# Calculate Visibilities. Index 0 = sun. Index 1 = satellite
baselines_exp = [[], []]
for slew in range(len(files)):
    print('Analyzing %s data' % types[slew])
    plotV = []
    errorV = []
    plotB = []
    errorB = []
    tmpV = []
    tmpB = []
    j = 0   # Counts how many measurements from the same plot we've considered
    for i in range(len(files[slew])):
        j += 1
        maxX = bigMaxX[slew][i]
        maxY = bigMaxY[slew][i]
        minY = nextMinY[slew][i]
        maxXNext = nextMaxX[slew][i]
        tmpV.append((maxY - minY) / (maxY + minY))
        tmpB.append(abs(maxXNext - maxX))
        if j == 3:    # Number of measurements per plot (must all be the same!)
            plotV.append(np.mean(tmpV))
            errorV.append(np.std(tmpV)/np.sqrt(3.0))
            plotB.append(np.mean(tmpB))
            errorB.append(np.std(tmpB)/np.sqrt(3.0))
            baselines_exp[slew].append(
                                2.0 * float(files[slew][i][24:26])/wavelength)
            tmpV= []
            tmpB= []
            j = 0
    # Fit this data
    plotB = baselines_exp[slew]
    if slew == 0:
        popt, pcov = curve_fit(sincFunc, plotB, plotV, p0=[1.0e+6/AU2km], 
                               sigma=errorV, absolute_sigma=True)
    else:
        popt, pcov = curve_fit(sincFunc, plotB, plotV)
    fittedB = np.arange(min(plotB), max(plotB),
                        step=( (max(plotB)-min(plotB))/100. ) )
    fittedV = []
    trueV = []
    for b in fittedB:
        fittedV.append(sincFunc(b, *popt))
        trueV.append(sincFunc(b, alpha_litVal))

   # Calculate Diameter of sun
    if slew == 0:
        alpha = popt[0]
        sigAlpha = np.sqrt(pcov[0])
        d = AU2km * alpha # small angle approx to find diameter
        sigd = AU2km * sigAlpha # prop. uncertainty
        
        # Literature Agreement
        litAgree = abs(d - d_litVal)/sigd[0]
        
        print('Xi^2 = ', pcov[0])
        print('alpha = {:e}'.format(alpha))
        print('d =  {:e} pm {:e} km'.format(d, sigd[0]))
        print('Agreement = {:e} sigma'.format(litAgree))

    # Plot this data
    plt.figure(slew)
    plt.errorbar(plotB, plotV, 
    #plt.errorbar(plotB, plotV, 
           xerr = errorB, yerr = errorV,
    #       xerr = [0.5]*len(errorV), yerr = errorV,
    fmt = '.', label='Analyzed Data'
    )
    if slew == 0:
        plt.plot(fittedB, fittedV, label='Fitted Sinc Function')
        plt.plot(fittedB, trueV,'--', label='Expected Sinc Function')
    plt.xlabel(r'$B_{\lambda}$')
    #plt.xlabel(r'$B$ (inches)')
    plt.ylabel(r'Visibility, $V_0(B_{\lambda})$')
    plt.minorticks_on()
    plt.legend()
    plt.title('%s Interferometer Visibility' % types[slew])
    plt.savefig(info['images'] + 'visibilities-%s.pdf' % types[slew] , ppi=300)
    #plt.savefig(info['images'] + 'visibilities-%s-noAltAdjustments.pdf' % types[slew] , ppi=300)
    #plt.savefig(info['images'] + 'visibilities-%s-BRaw.pdf'%types[slew],ppi=300)
    #plt.clf()

    # Save data
    fname = '../52c-' + types[slew] + '_prelimPlotdata.txt'
    np.savetxt(fname, np.c_[plotB, errorB, plotV, errorV])

