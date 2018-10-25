# Data Analysis for part 4.4
# Rescales flux and error of ref. stars, plots light curves
# October 16, 2018

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Load lab constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Open light curve tables

dataSubDir = info['fluxSubdir']
outSubDir = info['normFluxSubdir'] + 'NORM_'

# get fnames
files = open(dataSubDir + 'names.txt', 'r')
fnames = []
for line in files:
    fnames.append(line.strip('\n'))

i=0
for fname1 in fnames:
    data = np.loadtxt(dataSubDir + fname1, delimiter=',', skiprows=1)
    imageNumber = data[:,0]
    time = data[:,1]
    flux = data[:,2]
    sigFlux = data[:,3]

    avgFlux = np.mean(flux)
    flux = flux/avgFlux
    
    sigFlux = sigFlux/avgFlux
    
    #fname2 = outSubDir + fname1
    #np.savetxt(fname2, np.c_[imageNumber,time,flux,sigFlux])
    
    fname3 = outSubDir + fname1.strip('.txt') + '.png'
    title = fname1.strip('.txt')
    plt.figure(i+1)
    plt.title(title)
    plt.errorbar(time, flux, yerr=sigFlux, fmt='x', label=r'data')
    plt.xlabel(r'Time from Start [m]')
    plt.ylabel(r'Flux [count]')
    plt.savefig(fname3,format='png',dpi=1000,bbox_inches='tight')
    i=i+1

