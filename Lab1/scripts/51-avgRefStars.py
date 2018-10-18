# Data Analysis for part 4.4
# Rescales flux and error of ref. stars, plots light curves
# October 16, 2018

# Imports
import numpy as np
import matplotlib.pyplot as plt
import math

# Load lab constants
print('Load Constants')
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Open light curve tables

dataSubDir = info['fluxSubdir']

# get fnames
files = open(dataSubDir + 'names.txt', 'r')
fnames = []
for line in files:
    fnames.append(line.strip('\n'))

print('Pulling fluxes and times from other files')
imageNums = [[] for i in range(len(fnames))]
times = [[] for i in range(len(fnames))]
fluxs = [[] for i in range(len(fnames))]
fluxErrs = [[] for i in range(len(fnames))]
for fname1 in range(len(fnames)):
    data = np.loadtxt(dataSubDir + fnames[fname1], delimiter=',', skiprows=1)
    imageNums[fname1].append(data[:,0])
    times[fname1].append(data[:,1])
    fluxs[fname1].append(data[:,2])
    fluxErrs[fname1].append(data[:,3])

# Weighted Mean and std dev
print('Calculating for each exposure')
saveFile = open(info['normFluxSubdir'] + 'weightedAvgs.txt', 'w')
saveFile.write('#Time,flux,errFlux,weightMean,errWeightMean,ratio,ratioErr\n')
print(len(times[0]))
for image in range(len(times[0])):
    time = times[0][image]
    if image % 50 == 0:
        print(time)
    fluxSource = fluxs[0][image]
    fluxErrSource = fluxErrs[0][image]
    meanWeight = 0
    meanStdDev = 0
    meanWeightNum = 0
    meanWeightDen = 0
    for i in range(len(fluxs)):
        meanWeightNum = meanWeightNum + fluxs[i][image]/(fluxErrs[i][image]**2)
        meanWeightDen = meanWeightDen + 1/(fluxErrs[i][image]**2)
    print(meanWeightNum)
    print(meanWeightDen)
    meanWeight = meanWeightNum / meanWeightDen
    meanStdDev = math.sqrt(1 / meanWeightDen)
    ratio = fluxSource / meanWeight
    ratioErr = (1/meanWeight) * math.sqrt( fluxErrSource**2 + 
                 (ratio * meanStdDev)**2 )

    saveFile.write("%i,%.6f,%i,%i,%.6f,%.6f,%.6f,%.6f\n" % (
               imageNums[image], time, fluxSource, fluxErrSource, meanWeight, 
               meanStdDev, ratio, ratioErr ) )
saveFile.close
