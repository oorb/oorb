from __future__ import print_function
from __future__ import absolute_import
import os
import subprocess
import warnings
import numbers
import numpy as np
from itertools import repeat
import pandas as pd
import pyoorb as oo


import time


def dtime(time_prev):
    return (time.time() - time_prev, time.time())


def pack_oorbArray(orbits):
    """Translate orbital element dataframe (easy for humans) into pyoorb-suitable input orbit array."""
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
    if isinstance(orbits['objid'][0], numbers.Real):
        orbids = orbits['objid']
    else:
        orbids = np.arange(0, len(orbits['objid']), 1)
    # Convert to format for pyoorb, INCLUDING converting inclination, node, argperi to RADIANS
    oorbArray = np.column_stack((orbids, orbits['q'], orbits['e'], np.radians(orbits['i']),
                                 np.radians(orbits['node']), np.radians(orbits['argperi']),
                                 orbits['t_p'], elem_type, orbits['t_0'], epoch_type, orbits['H'], gval))
    return oorbArray


def unpack_oorbArray(orbits):
    """Translate pyoorb-style orbital element array back into dataframe."""
    orbits = pd.DataFrame(oorbArray, columns=['objId', 'q', 'e', 'i', 'node', 'argperi',
                                              't_p', 'elem_type', 't_0', 'epoch_type',
                                              'H', 'gval'])
    return orbits


if __name__ == "__main__":

    t = time.time()
    print("starting...")

    # check against OpenOrb command line
    newepoch = 49385
    try:
        os.unlink('test_orbits_new.des')
    except OSError:
        pass
    command = "oorb --task=propagation --orb-in=test_orbits.des" \
        "--orb-out=test_orbits_new.des --epoch-mjd-tt=%f" % newepoch
    print(command)
    subprocess.call(command, shell=True)
    dt, t = dtime(t)
    print("Propagating orbits by command line took %f s " % (dt))

    # Read the command line version back, to look for differences.
    data = pd.read_table('test_orbits_new.des', sep="\s*", engine='python')
    dt, t = dtime(t)
    print("Reading data back from file %f s" % (dt))
    print("Read %d orbits" % (len(data['q'])))

    # Read the orbits from disk.
    dt, t = dtime(t)
    orbits = pd.read_table('test_orbits.des', sep='\s*', engine='python')
    newcols = orbits.columns.values
    newcols[0] = 'objid'
    orbits.columns = newcols
    dt, t = dtime(t)
    print("Reading %d orbits required %f s" % (len(orbits['q']), dt))
    # convert orbit array to 'packed' form needed in oorb.
    oorbArray = pack_oorbArray(orbits)

    # set up oorb.
    ephfile = os.path.join(os.getenv('OORB_DATA'), 'de405.dat')
    oo.pyoorb.oorb_init(ephemeris_fname=ephfile)

    # For pyoorb, we need to tag times with timescales;
    # 1= MJD_UTC, 2=UT1, 3=TT, 4=TAI
    epochTime = np.array(zip([newepoch], repeat(3, len([newepoch]))), dtype='double')
    dt, t = dtime(t)
    print("Ready for orbit propagation .. %f s" % (dt))

    newOorbElems, err = oo.pyoorb.oorb_propagation_nb(in_orbits=oorbArray, in_epoch=epochTime)
    if err != 0:
        warnings.warn('Orbit propagation returned error %d' % err)
    print("Calculating new orbits by python required %f s" % (dt))
    orbits = unpack_oorbArray(newOorbElems)
    print(orbits)
