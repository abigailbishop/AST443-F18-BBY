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
ratioErrs = data[:,3]

# Calculate radius by averaging values in key regions
radiusTimes = [100, 145, 175, 200]
    # index: 0=preTransit, 1=beginTransit, 2=endTransit, 3=postTransit
radiusNumValues = [0]*( len(radiusTimes) - 1 )
    # index: 0=preTransit, 1=InTransit, 2=afterTransit
radiusFluxs = [0]*( len(radiusTimes) - 1 )
radiusFluxErrs = [0]*( len(radiusTimes) - 1 )
for image in range(len(imageNums)):
    if times[image] < radiusTimes[0]:
        radiusNumValues[0] += 1
        radiusFluxs[0] += ratios[image]
        radiusFluxErrs[0] += ratioErrs[image]**2
    if times[image] > radiusTimes[1] and times[image] < radiusTimes[2]:
        radiusNumValues[1] += 1
        radiusFluxs[1] += ratios[image]
        radiusFluxErrs[1] += ratioErrs[image]**2
    if times[image] > radiusTimes[3]:
        radiusNumValues[2] += 1
        radiusFluxs[2] += ratios[image]
        radiusFluxErrs[2] += ratioErrs[image]**2
print(radiusFluxs)
for i in range(len(radiusNumValues)):
    radiusFluxs[i] = radiusFluxs[i] / radiusNumValues[i]
    radiusFluxErrs[i] = radiusFluxErrs[i]**0.5 / radiusNumValues[i]
# NOTE the below calculations are based on only the flux before the transit
#   as the data after the transit looks unreliable
radiusFrac = ( 1 - radiusFluxs[1]/radiusFluxs[0] )**0.5
radiusFracErr = ( 0.5 * radiusFrac**0.5 * 
    ( (radiusFluxErrs[0]/radiusFluxs[0])**2 + 
      (radiusFluxErrs[1]/radiusFluxs[1])**2 )**0.5 ) 
print( 'Radius Fraction: ', radiusFrac )
print( 'Radius Fraction Error: ', radiusFracErr )
