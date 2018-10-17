import numpy as np
import os
from astropy.io import fits

info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

refs = []
for line in open('33-ourStars.txt'):
    li=line.strip()
    if not li.startswith("#"):
        stars = [x.strip() for x in line.split(',')]
        stars.append(stars[0])
        stars.append(stars[1])
        stars.append(stars[2])
    refs.append(stars)

fileNames = open(info['fluxSubdir'] + "names.txt","w")
for star in range(len(refs)):
    starFileName = "star_%s.txt" % (refs[star][0])
    fileNames.write(starFileName+'\n')
    print('Saving info from: ' + refs[star][0])
    file = open(info['fluxSubdir'] + starFileName,"w")
    numFiles = 0
    time0 = 0 
    for i in range(379):
        catalog = None
        pathToFits = None
        pathToCats = None
        if i < 10 :
            pathToFits = (info['dataDir']+info['rawDataSubdir']+
                      'data_.0000000{}.FIT'.format(i)) 
            pathToCats = (info['dataDir']+info['catFileSubdir']+
                      'data_0000000{}.cat'.format(i) )
        elif i < 100 :
            pathToFits = (info['dataDir']+info['rawDataSubdir']+
                      'data_.000000{}.FIT'.format(i)) 
            pathToCats = (info['dataDir']+info['catFileSubdir']+
                      'data_000000{}.cat'.format(i) )
        else: 
            pathToFits = (info['dataDir']+info['rawDataSubdir']+
                      'data_.00000{}.FIT'.format(i)) 
            pathToCats = (info['dataDir']+info['catFileSubdir']+
                      'data_00000{}.cat'.format(i) )
        if os.path.isfile(pathToFits) and os.path.isfile(pathToCats):
            # Save Time Info
            date = fits.open(pathToFits)[0].header['TIME-OBS']
            time = float(date[:2])*60. + float(date[3:5]) + float(date[6:])/60.
            # Save Flux Info
            catalog = np.loadtxt(pathToCats)
            fluxes = []
            fluxErres = []
            for thing in catalog:
                if ( (abs(thing[3]-float(refs[star][1])) < 0.05) and
                     (abs(thing[4]-float(refs[star][2])) < 0.05) ):
                    fluxes.append(thing[5])
                    fluxErres.append(thing[6])
            if len(fluxes) > 0 :
                if numFiles == 0:
                    time0 = time
                    file.write("#file number,time from start,flux,flux error\n")
                maxFlux = 0
                maxFluxErr = 0
                for n in range(len(fluxes)):
                    if fluxes[n] > maxFlux:
                        maxFlux = fluxes[n]
                        maxFluxErr = fluxErres[n]
                file.write("%s,%.6f,%i,%i\n" % (
                       i, time-time0, maxFlux, maxFluxErr) )
                numFiles = numFiles + 1
    file.close()
fileNames.close()
