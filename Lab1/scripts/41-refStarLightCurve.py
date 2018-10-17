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
N = 6 # number of ref stars
dataSubDir = 
for i in range(N):
    fname1 = FNAME1 # PUT FILENAME HERE
    data = np.loadtxt(fname1)
    time = data[0]
    flux = data[1]
    sigFlux = data[2]

    avgFlux = np.mean(flux)
    flux = flux/avgFlux
    
    avgSig = np.mean(sigFlux)
    sigFlux = sigFlux/avgSig
    
    fname2 = FNAME2 # PUT NEW FILE NAME HERE
    numpy.savetxt(fname2, np.c_[time,flux,sigFlux])
    
    fname3 = FNAME3 # PUT IMAGE FNAME HERE
    plt.figure(i+1)
    plt.title(r'Lightcurve of star', + str(i))
    plt.errorbar(time, flux, yerr=sigFlux, fmt='x', label=r'data')
    plt.xlabel(r'Time from Start [m]')
    plt.ylabel(r'Flux [count]')
    plt.savefig(fname3,format='png',dpi=1000,bbox_inches='tight')

