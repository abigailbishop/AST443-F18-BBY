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
for line in open('ourStars.txt'):
    li=line.strip()
    if not li.startswith("#"):
        stars = [x.strip() for x in line.split(',')]
        stars.append(stars[0])
        stars.append(stars[1])
        stars.append(stars[2])
    refs.append(stars)

for star in range(len(refs)):
    file = open(info['fluxSubdir'] + "star_%s.txt" % (refs[star][0]),"w")
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
            date = fits.open(pathToFits)[0].header['TIME-OBS']
            time = float(date[:2])*60. + float(date[3:5]) + float(date[6:])/60.
            if numFiles == 0:
                time0 = time
            catalog = np.loadtxt(pathToCats)
            print(date)
            np.loadtxt(pathToCats)
            flux=0
            fluxErr=0
            file.write("%s,%.6f,%i,%i\n" % (
                       refs[star][0], time-time0, flux, fluxErr) )
            numFiles = numFiles + 1
    file.close()
