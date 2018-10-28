# Data Analysis for 4.5 - Plot literal light curve

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

# Bin data
dt = 10       # Binning time interval in Minutes
maxTime = int( max(times) + 2*dt - max(times)%dt )
timesBinned = range(0,maxTime, dt)
ratiosBinned = [0]*len(timesBinned)
ratioErrsBinned = [0]*len(timesBinned)
numImagesBinned = [0]*len(timesBinned)
for image in range(len(imageNums)):
    index = int(times[image]/dt)
    ratiosBinned[index] = ratiosBinned[index] + ratios[image] 
    ratioErrsBinned[index] = ratioErrsBinned[index] + ratioErrs[image]**2
    numImagesBinned[index] = numImagesBinned[index] + 1
for ratio in range(len(ratiosBinned)):
    if numImagesBinned[ratio] > 1:
        ratiosBinned[ratio] = ratiosBinned[ratio] / numImagesBinned[ratio]
        ratioErrsBinned[ratio] = ( ratioErrsBinned[ratio]**0.5
            / numImagesBinned[ratio] )
plt.errorbar(timesBinned, ratiosBinned, yerr=ratioErrsBinned, fmt='x')


# Plot the plot
#plt.errorbar(times, ratios, fmt='x')
#plt.errorbar(imageNums, ratios, yerr=ratioErrs, fmt='x')
plt.ylim(0.93,1.07)
plt.title('Light Curve of Target')
plt.xlabel('Time from start (minutes)')
#plt.xlabel('File Number')
plt.ylabel('Relative Brightness')
plt.savefig('../lightCurve-Final.pdf', ppi=300)
