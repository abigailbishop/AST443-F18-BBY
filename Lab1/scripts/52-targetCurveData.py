# Data Analysis for part 4.5
# Rescales flux and error of ref. stars, plots light curves
# October 17, 2018

# Imports
import numpy as np
import matplotlib.pyplot as plt
import math

# Load lab constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Open scaled info
weighted = np.loadtxt(info['normFluxSubdir'] + 'weightedAvgs.txt', 
                      delimiter=',', skiprows=1)
fileNums = weighted[:,0]
times = weighted[:,1]
fluxs = weighted[:,2]
fluxErrs = weighted[:,3]
means = weighted[:,4]
stdDevs = weighted[:,5]
ratios = weighted[:,6]
ratioErrs = weighted[:,7]

# Weight the target flux
numFiles = 0
avgRatios = 0
avgRatioErrs = 0
for image in range(len(times)):
    if (fileNums[image] < 105 or fileNums[image] > 310):
        avgRatios = avgRatios + ratios[image]
        avgRatioErrs = avgRatioErrs + ratios[image]
        numFiles = numFiles + 1
avgFluxSource = avgFluxSource/numFiles
meanWeight = meanWeightNum / meanWeightDen
meanStdDev = math.sqrt(1 / meanWeightDen)
for image in range(len(times)):
    fluxs[image] = fluxs[image] / avgFluxSource
    fluxErrs[image] = fluxErrs[image] / avgFluxSource
print(fluxs[3])
print(fluxErrs[3])
print(avgFluxSource)

# Normalize ratios to the baseline flux
saveFile = open(info['normFluxSubdir'] + 'dataFinal-normToSource.txt', 'w')
saveFile.write('#ImageNum,Time,NormRi,NormSigRi\n')
for image in range(len(times)):
    ratiosNew = (ratios[image] / fluxs[image])
    ratioErrsNew = ( ratiosNew * math.sqrt( 
        (ratioErrs[image]/ratios[image])**2 + 
        (fluxErrs[image]/fluxs[image])**2 ) )
    if ratiosNew < 0.0001:
        saveFile.write('%d,%.6f,%.6f,%.6f\n' % (
                   fileNums[image], times[image], ratiosNew, ratioErrsNew ) )
saveFile.close()
