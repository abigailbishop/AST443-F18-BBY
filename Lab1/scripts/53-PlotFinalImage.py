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
dt = 5       # Binning time interval in Minutes
maxTime = int( max(times) + 2*dt - max(times)%dt )
timesBinned = np.arange(0.,maxTime, dt)
ratiosBinned = [[] for _ in range(len(timesBinned))]
ratioErrsBinned = [[] for _ in range(len(timesBinned))]
numImagesBinned = [0]*len(timesBinned)
for image in range(len(imageNums)):
    index = int(times[image]/dt)
    ratiosBinned[index].append(ratios[image])
    ratioErrsBinned[index].append(1/ratioErrs[image]**2)
    numImagesBinned[index] = numImagesBinned[index] + 1
timesBinnedPlot = []
ratiosBinnedPlot = []
ratioErrsBinnedPlot = []
for ratio in range(len(ratiosBinned)):
    if numImagesBinned[ratio] > 1:
        if timesBinned[ratio] > 250:
            continue
        timesBinnedPlot.append(timesBinned[ratio])
        ratiosBinnedPlot.append(np.mean(ratiosBinned[ratio]))
        sumErrors = 0
        for r in ratiosBinned[ratio]:
            sumErrors += (ratiosBinnedPlot[ratio] - r)**2
        ratioErrsBinnedPlot.append((1./(len(ratiosBinned[ratio])-1)
             *sumErrors)**0.5)
plt.errorbar(timesBinnedPlot, ratiosBinnedPlot, 
             yerr=ratioErrsBinnedPlot, fmt=',')

xs = np.arange(0,300)
ys = [1]*len(xs)

# Plot the plot
#plt.errorbar(times, ratios, fmt='x')
#plt.errorbar(imageNums, ratios, yerr=ratioErrs, fmt=',')
#plt.errorbar(times, ratios, yerr=ratioErrs, fmt=',')
plt.errorbar(xs, ys)
#plt.ylim(0.85,1.05)
plt.title('Light Curve of Target')
plt.xlabel('Time from start (minutes)')
#plt.xlabel('File Number')
plt.ylabel('Relative Brightness')
plt.savefig(info['images']+'lightCurve-Final.pdf', ppi=300)
