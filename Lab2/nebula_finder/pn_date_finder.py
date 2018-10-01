# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 18:16:50 2018

@author: emily
"""
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation

#from astropy.io.votable import parse_single_table
#table = parse_single_table("5246234-Planetaries.xls")
data = np.loadtxt("pn_data.csv",skiprows=1,usecols = (2,3,4,5),
                  dtype='str', delimiter=',')

NAME = data[:,0]
RA_tmp = data[:,1]  # In hours and minutes
RA = []
DEC_tmp = data[:,2] # In degrees but instead of a . it's a space
DEC = []
MAG_tmp = data[:,3] # As a string
MAG = []

for ra in RA_tmp:
    ra_hour = ra[0:2]
    ra_min = ra[3:]
    RA.append( ( float(ra_hour) + float(ra_min) / 60 ) * 360 / 24 )

for dec in DEC_tmp:
    dec_int = dec[1:3]
    dec_decimal = dec[4:]
    if dec[0] == '+':
        DEC.append( float(dec_int) + float(dec_decimal) / 100 )
    else:
        DEC.append( 0 - float(dec_int) + float(dec_decimal) / 100 )

for mag in MAG_tmp:
    MAG.append( float(mag) )

JD = 2458408.625
            
output='transit_aau.csv'
# location of Stony Brook
observing_location = EarthLocation(lat=40.914224*u.deg, lon=-73.11623*u.deg, height=0)
# read in the input file
#transits=np.loadtxt(transit_rdj.txt)
# define the Ra, Dec pairs
coords=SkyCoord(RA[:]*u.deg,DEC[:]*u.deg, frame='icrs')
# define the observing times
times=Time(JD, format='jd')
# do the transformation to AltAz
aa = AltAz(location=observing_location, obstime=times)
aacoords=coords.transform_to(aa)
# specify output timeformat to be in UTC
times.format='fits'
# write output file
#np.savetxt(output,np.transpose([aacoords.az,aacoords.alt,times]), fmt="%.6f %.6f %.30s")

NAME_obs = []
RA_obs = []
DEC_obs = []
Alt_obs = []
Az_obs = []
V_obs = []
for i in range(len(aacoords.alt.deg)):
    if aacoords.alt.deg[i]>40 and MAG[i]<12:
        NAME_obs.append(NAME[i])
        RA_obs.append(RA[i])
        DEC_obs.append(DEC[i])
        Alt_obs.append(aacoords.alt.deg[i])
        Az_obs.append(aacoords.az.deg[i])
        V_obs.append(MAG[i])

#np.savetxt('obs_transit_aau.txt',np.transpose([ID_obs,Az_obs,Alt_obs,Time_obs]), fmt="%.30s %.6f %.6f %.30s")
with open("nebula_options.txt", "w") as param_file:
    header = "DATE:" + str(JD) + "\n NAME RA DEC AZ ALT V \n"
    param_file.write(header)
    for i in range(len(NAME_obs)):
        outstring = str(NAME_obs[i]) + " " 
        outstring += str(RA_obs[i]) + " " 
        outstring += str(DEC_obs[i]) + " " 
        outstring += str(Az_obs[i]) + " " 
        outstring += str(Alt_obs[i]) + " " 
        outstring += str(V_obs[i]) + " " + "\n"
        param_file.write(outstring)
        
        
        
        
    
