import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from itertools import repeat
import pandas as pd
import pyoorb as oo


import time
def dtime(time_prev):
   return (time.time() - time_prev, time.time())


def pack_oorbArray(orbits):
    """Translate orbital element dictionary (easy for humans) into pyoorb-suitable input orbit array."""
    # Translate orbital elements into array that pyoorb will like.
    # PyOrb wants ::
    # 0: orbitId
    # 1 - 6: orbital elements, using radians for angles
    # 7: element type code, where 2 = cometary - means timescale is TT, too
    # 8: epoch
    # 9: timescale for the epoch; 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
    # 10: magHv
    # 11: G
    elem_type = np.zeros(len(orbits)) + 2
    epoch_type = np.zeros(len(orbits)) + 3
    gval = np.zeros(len(orbits)) + 0.15
    # Also, the orbitID has to be a float, rather than a string, so substitute if needed.
    if ((isinstance(orbits['objid'][0], float) == True) |
        (isinstance(orbits['objid'][0], int) == True)):
        orbids = orbits['objid']
    else:
        orbids = np.arange(0, len(orbits['objid']), 1)
    # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
    oorbArray = np.column_stack((orbids, orbits['q'], orbits['e'], np.radians(orbits['i']),
                                    np.radians(orbits['node']), np.radians(orbits['argperi']),
                                    orbits['t_p'], elem_type, orbits['t_0'], epoch_type, orbits['H'], gval))
    return oorbArray


if __name__ == "__main__":

    t = time.time()
    print "starting..."

    # check against OpenOrb command line
    timespan = 1000
    timestep = 10
    command = "/bin/rm test_out ; oorb --code=807 --task=ephemeris --orb-in=test_orbits.des --timespan=%d --step=%d > test_out" %(timespan, timestep)
    print command
    subprocess.call(command, shell=True)
    dt, t = dtime(t)
    print "Calculating ephemerides by command line took %f s "%(dt)

    # Read the command line version back, to look for differences.
    data = pd.read_table('test_out', sep="\s*", engine='python')
    dt, t = dtime(t)
    print "Reading data back from file %f s" %(dt)
    print "Read %d ephemerides" %(len(data['RA']))
    ctimes = data['MJD_UTC/UT1'][data['#Designation'] == 1]
    print "Created %d unique times; %d times total" %(len(np.unique(data['MJD_UTC/UT1'])), len(ctimes))

    # Read the orbits from disk.
        # Read the orbits from disk.
    dt, t = dtime(t)
    orbits = pd.read_table('test_orbits.des', sep='\s*', engine='python')
    newcols = orbits.columns.values
    newcols[0] = 'objid'
    orbits.columns = newcols
    dt, t = dtime(t)
    print "Reading %d orbits required %f s" %(len(orbits['q']), dt)
    # convert orbit array to 'packed' form needed in oorb.
    oorbArray = pack_oorbArray(orbits)

    # set up oorb.
    ephfile = os.path.join(os.getenv('OORB_DATA'), 'de430.dat')
    oo.pyoorb.oorb_init(ephemeris_fname=ephfile)

    # set observatory code
    obscode = 807
    # Set up dates to predict ephemerides.
    timestart = orbits['t_0'][0]
    times = np.arange(timestart, timestart + timespan + timestep/2.0, timestep)
    times = np.array(ctimes)
    # For pyoorb, we need to tag times with timescales;
    # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
    ephTimes = np.array(zip(times, repeat(1, len(times))), dtype='double')
    print times
    dt, t = dtime(t)
    print "Ready for ephemeris generation .. %f s" %(dt)

    # Generate ephemerides.
    oorbephs, err = oo.pyoorb.oorb_ephemeris(in_orbits = oorbArray, in_obscode=obscode, in_date_ephems=ephTimes)
    dt, t = dtime(t)
    print "Calculating ephemerides by python required %f s" %(dt)

    # Returned ephems contain a 3-D Fortran array of ephemerides, the axes are:
    #   [objid][time][ephemeris information element]
    # the ephemeris information elements are (in order):
    # distance, ra, dec, mag, ephem mjd, ephem mjd timescale, dradt(sky), ddecdt(sky)
    # per object, per date, 8 elements (array shape is OBJ(s)/DATE(s)/VALUES)
    # Note that ra/dec, dradt, etc. are all in DEGREES.
    # First: (to arrange ephems for easier later use)
    # Swap the order of the axes: DATE / Objs / values

    # Unpack ephemerides.
    times = np.ravel(oorbephs.swapaxes(0, 1).swapaxes(0, 2)[4])
    ra = np.ravel(oorbephs.swapaxes(0, 1).swapaxes(0, 2)[1])
    dec = np.ravel(oorbephs.swapaxes(0, 1).swapaxes(0, 2)[2])

    radiff = data['RA'] - ra
    decdiff = data['Dec'] - dec

    radiff *= 3600
    decdiff *= 3600

    print "min/max ra offsets", radiff.min(), radiff.max()
    print "min/max dec offsets", decdiff.min(), decdiff.max()

    plt.figure()
    plt.plot(radiff, decdiff, 'k.')
    plt.xlabel('Difference in RA (arcsec)')
    plt.ylabel('Difference in Dec (arcsec)')

    #plt.show()
