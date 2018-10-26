# This code will save the flux and error in flux of all notable stars and
#   track the time that each flux is measured at relative to the first time

import numpy as np
import os
from astropy.io import fits
import math

# Import path and file name constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Import star names, RA, and Dec
refs = []
for line in open('33-ourStars.txt'):
    li=line.strip()
    if not li.startswith("#"):
        stars = [x.strip() for x in line.split(',')]
        stars.append(stars[0])
        stars.append(stars[1])
        stars.append(stars[2])
    refs.append(stars)

#Loop over every stat in 33-ourStars.txt
numExposures = 2700
fileNames = open(info['fluxSubdir'] + "names.txt","w")
time0 = 0 
saveData = [[] for i in range(numExposures)]
openFiles = []
for star in range(len(refs)):
    # Get file name and save it. Open save-to file
    starFileName = "star_%s.txt" % (refs[star][0])
    fileNames.write(starFileName+'\n')
    print('Saving info from: ' + refs[star][0])
    openFiles.append(open(info['fluxSubdir'] + starFileName,"w"))

    # Loop over every exposure for flux and time
    numFiles = 0
    j = 26
    for i in range(len(saveData)):
        j=j+1
        # Determines path to catalog and fits file
        catalog = None
        pathToFits = None
        pathToCats = None
        #if j < 10 :
        #    pathToFits = (info['dataDir']+info['rawDataSubdir']+
        #              'hd189733.0000000{}.FIT'.format(i)) 
        #    pathToCats = (info['dataDir']+info['catFileSubdir']+
        #              'hd189733_0000000{}.cat'.format(i) )
        if j < 100 :
            pathToFits = (info['dataDir']+info['rawDataSubdir']+
                      'hd189733.000000{}.FIT'.format(i)) 
            pathToCats = (info['dataDir']+info['catFileSubdir']+
                      'hd189733_000000{}.cat'.format(i) )
        elif j < 1000 :
            pathToFits = (info['dataDir']+info['rawDataSubdir']+
                      'hd189733.00000{}.FIT'.format(i)) 
            pathToCats = (info['dataDir']+info['catFileSubdir']+
                      'hd189733_00000{}.cat'.format(i) )
        else: 
             pathToFits = (info['dataDir']+info['rawDataSubdir']+
                      'hd189733.0000{}.FIT'.format(i)) 
             pathToCats = (info['dataDir']+info['catFileSubdir']+
                      'hd189733_0000{}.cat'.format(i) )
        # If both the fits file and the catalog exist, analyze it
        if os.path.isfile(pathToFits) and os.path.isfile(pathToCats):
            # Save Time of exposure
            date = fits.open(pathToFits)[0].header['TIME-OBS']
            time = float(date[:2])*60. + float(date[3:5]) + float(date[6:])/60.
            # Save flux info for every thing close to where we expect
            #   the star to be
            catalog = np.loadtxt(pathToCats)
            fluxes = []
            fluxErres = []
            for thing in catalog:
                if math.sqrt( ((thing[3]-float(refs[star][1]))*
                                  math.cos(thing[4]/180*math.pi))**2 + 
                     (thing[4]-float(refs[star][2]))**2) < 0.01 :
                    fluxes.append(thing[5])
                    fluxErres.append(thing[6])
            # If we found more than one thing close to where we expected the 
            #   star to be, then we save information from this exposure
            if len(fluxes) > 0 :
                if (numFiles == 0 and star==0):
                    time0 = time
                maxFlux = 0
                maxFluxErr = 0
                # If more than one objects were found, pull the one with
                #   higher flux
                for n in range(len(fluxes)):
                    if fluxes[n] > maxFlux:
                        maxFlux = fluxes[n]
                        maxFluxErr = fluxErres[n]
                saveData[i].append( [i, time-time0, maxFlux, maxFluxErr] )
                numFiles = numFiles + 1

for fileI in openFiles:
    fileI.write("#file number,time from start,flux,flux error\n")

for i in range(len(saveData)):
     if len(saveData[i]) == len(refs):
         for star in range(len(refs)):
             openFiles[star].write("%s,%.6f,%i,%i\n" % (saveData[i][star][0], 
               saveData[i][star][1],saveData[i][star][2],saveData[i][star][3]))
for fileI in openFiles:
    fileI.close()
fileNames.close()
