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

# Average the ratios
numFiles = 0
avgRatios = 0
#badImages = ( [2] + range(35, 45) + [66, 67] + range(77, 95) + range(109,134) +
#     range(161, 179) + range(226, 229) + [243] + range(251, 356) + 
#     range(368, 379) ) 
badImages = ( [2, 66, 67] + range(77,95) + range(226,229) + [243] )
for image in range(len(fileNums)):
    if ((fileNums[image] < 105 or fileNums[image] > 310) and 
        (fileNums[image] not in badImages) ):
        avgRatios = avgRatios + ratios[image]
        numFiles = numFiles + 1
avgRatios = avgRatios / numFiles
for image in range(len(times)):
    ratios[image] = ratios[image] / avgRatios
    ratioErrs[image] = ratioErrs[image] / avgRatios
print(ratios[3])
print(ratioErrs[3])

# Normalize ratios to the baseline ratio
saveFile = open(info['normFluxSubdir'] + 'dataFinal-normToSource.txt', 'w')
saveFile.write('#ImageNum,Time,NormRi,NormSigRi\n')
for image in range(len(fileNums)):
    if (ratios[image] < 2 and (fileNums[image] not in badImages)):
        saveFile.write('%d,%.6f,%.6f,%.6f\n' % (
            fileNums[image], times[image], ratios[image], ratioErrs[image]))
saveFile.close()
