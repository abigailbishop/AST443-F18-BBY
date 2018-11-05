"""

AST 443
Lab2
Data Analysis
7.1 Dark Current Correction

"""

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import stats
from scipy.stats import norm

# Load Lab constants
info = {}
for line in open('inputs.txt'):
    li=line.strip()
    if not li.startswith("#"):
        data = [x.strip() for x in line.split(',')]
        info[data[0]] = data[1]

# Get Raw Science Images Names
sciDirNames = ['rawNebulaData', 'rawRefstarData','arcSpecSubdir','flatSubdir']
SciNames = []
for subdir in sciDirNames:
  sciDir = info['dataDir'] + info[subdir]
  files = open(sciDir + 'names.txt', 'r')
  for line in files:
    SciNames.append(sciDir + line.strip('\n'))

# Make master dark for each exposure time

subdirNames = ['darkSubdir10s','darkSubdir30s','darkSubdir2min',
               'darkSubdir3min','darkSubdir5min','darkSubdir10min']

for subdir in subdirNames:
  dark_dir = info['dataDir'] + info[subdir]
  # get fnames
  files = open(dark_dir + 'names.txt', 'r')
  fnames = []
  for line in files:
      fnames.append(line.strip('\n'))
  
  # make master dark
  allimages = []
  i=0
  for file in fnames:
    fname = dark_dir + file
    hdu = fits.open(fname)
    head = hdu[0].header
    if i==0:
      darkExpTime = head['EXPTIME']
      i=i+1
    allimages.append(hdu[0].data)
    hdu.close()
  masterdark = np.median(allimages,axis=0)
  path = info['fitsFiles']
  mdname = path + 'masterdark_' + str(darkExpTime) + '.fits'
  #fits.writeto(mdname,masterdark,header=head)
  

  for fname in SciNames:
    hdu = fits.open(fname)
    head = hdu[0].header
    sciExpTime = head['EXPTIME']
    data = hdu[0].data
    hdu.close()
    if sciExpTime == darkExpTime:
      cordata = data - masterdark
      calfname = fname.strip('.FIT') + '_dcc.fits'
      fits.writeto(calfname,cordata,header=head)

  print(subdir, 'done')
