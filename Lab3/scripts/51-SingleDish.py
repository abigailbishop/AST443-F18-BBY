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

# Analyze Sun and Sat Data
sun = np.loadtxt(
    info['dataDir'] + info['raw1dish'] + 'single_01_140_170_sun.csv', 
    skiprows=1, delimiter=',')
sun_times = sun[:,0]
sun_currents = sun[:,1]
sat = np.loadtxt(
    info['dataDir'] + info['raw1dish'] + 'single_05_120_140_satellite_new.csv', 
    skiprows=1, delimiter=',')
sat_times=sat[:,0]
sat_currents = sat[:,1]

# Flip and center the currents
def centerflip(array):
    center = np.mean(array[-5:])
    for i in range(len(array)):
        array[i] = -1 * (array[i] - center)
    return(array)
sun_currents = centerflip(sun_currents)
sat_currents = centerflip(sat_currents)

# Resize Graphs
max_current = max([max(sun_currents), max(sat_currents)])
sun_currents = [i*max_current/max(sun_currents) for i in sun_currents]
sat_currents = [i*max_current/max(sat_currents) for i in sat_currents]
sat_times = [i+11 for i in sat_times]
sat_alt = 26. * np.pi / 180.
sun_alt = 30. * np.pi / 180.
sun_times = [i*(np.cos(sat_alt)/np.cos(sun_alt)) for i in sun_times]

# Signal to Noise Ratios
sun_max = 0
sat_max = 0
sun_min = 20
sat_min = 20
for t in range(len(sun_times)):
    if sun_times[t] < 18 or sun_times[t] > 32:
        sun_max = sun_currents[t] if sun_currents[t] > sun_max else sun_max
        sun_min = sun_currents[t] if sun_currents[t] < sun_min else sun_min
for t in range(len(sat_times)):
    if sat_times[t] < 18 or sat_times[t] > 32:
        sat_max = sat_currents[t] if sat_currents[t] > sat_max else sat_max
        sat_min = sat_currents[t] if sat_currents[t] < sat_min else sat_min
sat_SNR = min(sat_currents) / (sat_max - sat_min)
sun_SNR = min(sun_currents) / (sun_max - sun_min)
print('Satellite SNR: ', sat_SNR)
print('Sun SNR:       ', sun_SNR)

# Full width and half the max distance from original current
def fwhm(x, y):
    diff = max(y) - min(y)
    hm = 0.5 * diff
    maxIndex = 0
    hmIndex = 0
    for i in range(len(y)):
        if y[i] == max(y):  maxIndex = i
        if abs(y[i] - hm) < abs(y[hmIndex] - hm): hmIndex = i
    fwhm = 2 * abs( x[maxIndex] - x[hmIndex] )
    return( fwhm, maxIndex, hmIndex)
fwhmSun, SunMaxIndex, SunHMIndex = fwhm(sun_times, sun_currents)
fwhmSat, SatMaxIndex, SatHMIndex = fwhm(sat_times, sat_currents)
print( 'Sun fwhm ', fwhmSun)
print( 'Sat fwhm ', fwhmSat)

# Plot data
plt.plot(sat_times, sat_currents, label='Satellite')
plt.plot(sun_times, sun_currents, label='Sun')
plt.xlabel('Time (s)')
plt.ylabel('Normalized Current (A)')
plt.title('Single Dish Observations')
plt.legend()
plt.savefig(info['images'] + 'profile_single_sat05.pdf', ppi=300)

# Getting Sigma
sunSigma = fwhmSun / sun_SNR
satSigma = fwhmSat / sat_SNR
print('SunSigma = ', sunSigma, 'SatSigma = ', satSigma)
sunSigma = fwhmSun / sun_SNR
