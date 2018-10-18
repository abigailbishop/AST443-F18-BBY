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

# Plot the plot
plt.errorbar(times, ratios, yerr=ratioErrs, fmt='x')
plt.title('Light Curve of Target')
plt.xlabel('Time from start (minutes)')
plt.ylabel('Relative Brightness')
plt.savefig('../lightCurve-Final.pdf', ppi=300)
