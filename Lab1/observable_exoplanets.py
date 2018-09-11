# September 3 2018
# Saving stars with transiting exoplanets observable from SBU to file

import numpy as np
import sys
from astropy.io.votable import parse
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation

# Load exoplanet catalog and initialize array
votable = parse("exoplanet_catalog.vot.xml")
table = votable.get_first_table()
data = table.array
data_sorted = []
names = []

# Select candidates from magnitude, transit duration, and depth
count = 0
for i in range( len( data ) ):
    transit_depth = (data['radius'][i] * 0.1028 / data['star_radius'][i])**2
    magnitude = data['mag_v'][i]
    transit_duration = ( (data['star_radius'][i] * data['orbital_period'][i]) / 
                         (np.pi * data['semi_major_axis'][i] * 215) ) 
    JD_planet = data['tzero_tr'][i]
    orbital_period = data['orbital_period'][i]
    if( ( transit_depth >= 0.01 )
        and   
        ( magnitude <= 12 )
        and
        ( transit_duration <= 0.125 )
        and 
        ( JD_planet > 0 )
        ): 
        count = count + 1
#        print( 'duration: {}'.format(transit_duration * 24))
#        print( 'depth: {}'.format(transit_depth))
#        print( 'magnitude: {}'.format(magnitude))
#        print( 'JD: {}'.format(JD_planet))
        # Updates planet crossings to after Sept. 4
        while (JD_planet < 2458367.500000 ) :
            JD_planet = JD_planet + orbital_period
        # All planet transits before Sept. 22
        while (JD_planet < 2458383.500000 ) :
            data_sorted.append( [
                i,
                data['ra'][i],
                data['dec'][i],
                JD_planet,
                transit_depth,
                magnitude,
                transit_duration
                ] )
            names.append( [
                data['name'][i]
                ] )
            JD_planet = JD_planet + orbital_period
        # Updates planet crossings to after Oct. 1
        while (JD_planet < 2458392.500000 ) :
            JD_planet = JD_planet + orbital_period
        # All planet transits before Oct. 6
        while (JD_planet < 2458397.500000 ) :
            data_sorted.append( [
                i,
                data['ra'][i],
                data['dec'][i],
                JD_planet,
                transit_depth,
                magnitude,
                transit_duration
                ] )
            names.append( [
                data['name'][i]
                ] )
            JD_planet = JD_planet + orbital_period

# See if the transit is observable from SBU
output = sys.argv[1]
transits = np.asarray(data_sorted)
observing_location = EarthLocation(
    lat=40.914224*u.deg, lon=-73.11623*u.deg, height=0)
coords=SkyCoord(transits[:,1]*u.deg, transits[:,2]*u.deg, frame='icrs')
times=Time(transits[:,3], format='jd')
aa = AltAz(location=observing_location, obstime=times)
aacoords=coords.transform_to(aa)
times.format='fits'
obs_names = []
obs_az = []
obs_alt = []
obs_RA = []
obs_dec = []
obs_times = []
obs_depth = []
obs_mag = []
obs_duration = []
for i in range( len (transits) ) :
    if ( 
        int( str(times[i])[11:13] ) < 6 
        and
        int( str(times[i])[11:13] ) > 1 
        and
        aacoords.alt[i].deg > 40 
        ) :
        obs_names.append( names[i] )
        obs_az.append( aacoords.az[i].deg )
        obs_alt.append( aacoords.alt[i].deg )
        obs_RA.append( transits[i,1] )
        obs_dec.append( transits[i,2] )
        obs_times.append( times[i] )
        obs_depth.append( transits[i,4] )
        obs_mag.append( transits[i,5] )
        obs_duration.append( transits[i,6]*24. )
 
print('Number of Candidates: {}'.format( count ) )
print('Number of Opportunities: {}'.format( len(obs_names) ) )

np.savetxt(output,
    np.transpose( [ obs_names, obs_az, obs_alt, obs_RA, obs_dec, obs_times, 
        obs_depth, obs_mag, obs_duration ] ), 
    header = 'Names\tAz(deg)\tAlt(deg)\tRA(deg)\tDec(deg)\tTime\tDepth\tMag\tDuration(hrs)',
    fmt="%s %.6f %.6f %.6f %.6f %.30s %.6f %.6f %.6f")
