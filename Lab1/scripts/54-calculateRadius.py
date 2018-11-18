# Data Analysis for part 5 or 6 - Calculating the star's radius

import matplotlib.pyplot as plt
import numpy as np

# Import Constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Load Data
data = np.loadtxt(info['normFluxSubdir'] + 'dataFinal-normToSource.txt', 
                  delimiter=',', skiprows=1)
imageNums = data[:,0]
times = data[:,1]
ratios = data[:,2]

# Calculate radius by averaging values in key regions
radiusTimes = [120, 140, 215, 215]
    # index: 0=preTransitEnd, 1=beginTransit, 2=endTransit, 3=postTransitBegin
radiusNumValues = [0]*( len(radiusTimes) - 1 )
    # index: 0=preTransit, 1=InTransit, 2=afterTransit
radiusFluxs = [0]*( len(radiusTimes) - 1 )
radiusFluxStds = [0]*( len(radiusTimes) - 1 )
# Sum all ratios, count exposures in each set of data
for image in range(len(imageNums)):
    if times[image] > 250:
        continue
    if times[image] < radiusTimes[0]:
        radiusNumValues[0] += 1
        radiusFluxs[0] += ratios[image]
    if times[image] > radiusTimes[1] and times[image] < radiusTimes[2]:
        radiusNumValues[1] += 1
        radiusFluxs[1] += ratios[image]
    if times[image] > radiusTimes[3]:
        radiusNumValues[2] += 1
        radiusFluxs[2] += ratios[image]
# Calculate mean ratio in each zone
for i in range(len(radiusFluxs)):
    radiusFluxs[i] = radiusFluxs[i] / radiusNumValues[i]
# Calculate (science - mean)^2 and sum them for each zone
for image in range(len(imageNums)):
    if times[image] > 250:
        continue
    if times[image] < radiusTimes[0]:
        radiusFluxStds[0] += (ratios[image] - radiusFluxs[0])**2
    if times[image] > radiusTimes[1] and times[image] < radiusTimes[2]:
        radiusFluxStds[1] += (ratios[image] - radiusFluxs[1])**2
    if times[image] > radiusTimes[3]:
        radiusFluxStds[2] += (ratios[image] - radiusFluxs[2])**2
# Calculate STD then uncertainty on the mean for each zone
for i in range(len(radiusFluxs)):
    # Standard Deviation
    radiusFluxStds[i] = (radiusFluxStds[i] / radiusNumValues[i])**0.5
    print(radiusFluxStds[i])
    # Uncertainty on the mean
    radiusFluxStds[i] = radiusFluxStds[i] / radiusNumValues[i]**0.5
    
# NOTE the below calculations are based on only the flux before the transit
#   as the data after the transit looks unreliable
depth = radiusFluxs[0] - radiusFluxs[1]
depthErr = ( radiusFluxStds[0]**2 + radiusFluxStds[1]**2 )**0.5
radiusFrac = depth**0.5
radiusFracErr = depthErr / depth**0.5 * 0.5 
save = open(info['images'] + 'finalRadius.txt', 'w')
save.write('Depth: %f +- %f\n' % (depth, depthErr) )
save.write('R_planet / R_star: %f +- %f\n' % (radiusFrac, radiusFracErr) )
save.write('preTransit: %f +- %f\n' % (radiusFluxs[0], radiusFluxStds[0]) )
save.write('inTransit: %f +- %f\n' % (radiusFluxs[1], radiusFluxStds[1]) )
save.write('afterTransit: %f +- %f\n' % (radiusFluxs[2], radiusFluxStds[2]) )
save.close()
