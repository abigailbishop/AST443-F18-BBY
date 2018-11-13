# Analyze the single dish observations for Lab 3

import numpy as np
import matplotlib.pyplot as plt
import math

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

# Defines gaussian function with standard deviation (full-width as half maximum) sigma
def gaussian(x,sigma):
    return 1/(2*math.sqrt(2*math.pi)) * math.exp(-(x/(2*sigma))**2)
n=1000 # number of points
t=0.001 # sampling interval
hm = min(sun_currents) * 0.5
time_at_hm = 26
ws = [0, 0]
for c in range(len(sun_currents)):
    if abs(sun_currents[c]-hm) < abs(ws[0] - hm) and sun_times[c] < time_at_hm:
        ws[0] = sun_currents[c]
    elif abs(sun_currents[c]-hm)<abs(ws[1] - hm) and sun_times[c] > time_at_hm:
        ws[1] = sun_currents[c]
fwhm = ws[1] - ws[0]
resolved_source_width = 1.0 # radians (make it an even # of t's)
                            # exaggerated: actual solar value is
                            # 0.01 rad
# define n dimensional vector of the gaussian function evaluated over the sampling interval 
gaus = [gaussian( (i- n/2.0 + 1)*t, fwhm ) for i in range(n)]
# define top hat over sampling interval
tophat = [0]*n
for i in range(n):
    if(i**2<(resolved_source_width/t/2)**2):
        tophat[i] = 1.0
# Convolution of guassian and top hat
gauss_hat = numpy.convolve(gaus, tophat)
# plot results
y = [(i- n/2.0 + 1)*t/2 for i in range(2*n-1)]
print gauss_hat

#print x# Plot data
plt.plot(sat_times, sat_currents, label='Satellite')
plt.plot(sun_times, sun_currents, label='Sun')
plt.plot(y,gauss_hat)
plt.xlabel('Time (s)')
plt.ylabel('Normalized Current (A)')
plt.title('Single Dish Observations')
plt.legend()
plt.savefig(info['images'] + 'profile_single_sat05.pdf', ppi=300)
