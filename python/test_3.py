# ================================================================== #
#                                                                    #
# Copyright 2002,2003,2004,2005,2006,2007,2008,2009                  #
# Mikael Granvik, Jenni Virtanen, Karri Muinonen, Teemu Laakso,      #
# Dagmara Oszkiewicz                                                 #
#                                                                    #
# This file is part of OpenOrb.                                      #
#                                                                    #
# OpenOrb is free software: you can redistribute it and/or modify it #
# under the terms of the GNU General Public License as published by  #
# the Free Software Foundation, either version 3 of the License, or  #
# (at your option) any later version.                                #
#                                                                    #
# OpenOrb is distributed in the hope that it will be useful, but     #
# WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                           #
#                                                                    #
# You should have received a copy of the GNU General Public License  #
# along with OpenOrb. If not, see <http://www.gnu.org/licenses/>.    #
#                                                                    #
# ================================================================== #
#
#

from __future__ import print_function
from __future__ import absolute_import

import os
import numpy as np
import pylab
from itertools import repeat
import pyoorb as oo

import time


def dtime(time_prev):
    return (time.time() - time_prev, time.time())


def read_orbits(orbitfilename):
    """ Read orbits that will be using for orbit prediction """
    """ set up for reading SSM des files"""
    try:
        f = open(orbitfilename, 'r')
    except IOError:
        print("Couldn't open file %s" % (orbitfilename))
        exit()
    q = []
    e = []
    i = []
    node = []
    argperi = []
    timeperi = []
    H = []
    epoch = []
    # read values from file
    for line in f:
        if line.startswith("!"):
            continue
        if line.startswith("#"):
            continue
        values = line.split()
        q.append(values[2])
        e.append(values[3])
        i.append(values[4])
        node.append(values[5])
        argperi.append(values[6])
        timeperi.append(values[7])
        H.append(values[8])
        epoch.append(values[9])
    # convert to numpy arrays
    q = np.array(q, dtype=np.double)
    e = np.array(e, dtype=np.double)
    i = np.array(i, dtype=np.double)
    node = np.array(node, dtype=np.double)
    argperi = np.array(argperi, dtype=np.double)
    timeperi = np.array(timeperi, dtype=np.double)
    epoch = np.array(epoch, dtype=np.double)
    H = np.array(H, dtype=np.double)
    # add some extra columns required for ephemeris calculation
    objid = np.arange(0.0, len(q), 1.0, dtype=np.double)
    G = np.ones((len(q),), dtype=np.double)
    G = G * 0.15
    # convert values to radians
    i = np.radians(i)
    node = np.radians(node)
    argperi = np.radians(argperi)
    # and add the bookkeeping info for cometary type (2) and timescale (1)
    timescale = np.ones((len(q),), dtype=np.double)
    element_type = np.ones((len(q),), dtype=np.double)
    element_type = element_type * 2.0

    # convert to the proper sort of array for pyoorb
    orbits = np.column_stack((objid, q, e, i, node, argperi, timeperi, element_type, epoch, timescale, H, G))
    orbital_epoch = epoch[0]
    return orbits, orbital_epoch


if __name__ == "__main__":

    singledate = False  # calculate ephemerides at a single time at some separation from epoch.
    # otherwise - calculate ephemerides at multiple times between epoch and future date.

    t = time.time()
    print("Starting oorb ephemeris generation timing test")
    dt, t = dtime(t)
    ephfile = os.path.join(os.getenv('OORB_DATA'), 'de405.dat')
    oo.pyoorb.oorb_init(ephemeris_fname=ephfile, error_verbosity=5, info_verbosity=1)
    dt, t = dtime(t)
    print("calling oorb_init() took %f s" % (dt))
    # read in orbit DES file
    orbits, orbital_epoch = read_orbits('test_orbits.des')
    dt, t = dtime(t)
    print("Reading %d orbits required %f s" % (len(orbits), dt))

    # set observatory code
    obscode = 807
    # Generate ephemerides for ten years, at timesteps of 30 days.
    dates = np.arange(orbital_epoch, orbital_epoch + 2 * 365.10 + 30.0 / 2., 30.0)
    dates = np.arange(orbital_epoch, orbital_epoch + 30. + 1. / 24., 1. / 24. * .75)
    # timescale: 1 = UTC, 2 = UT1, 3=TT, 4=TAI
    # Create an array to store the timing information.
    tReq = np.zeros(len(dates), dtype='float')
    scaling = 1000
    for i in range(1, len(dates)):
        if singledate:
            # Calculate ephemerides at only the 'date'
            ephem_dates = np.zeros([1, 2], dtype=np.double)
            ephem_dates[0][:] = [dates[i], 1]
        else:
            # Calculate ephemerides at all times up to date
            d = dates[0:i]
            ephem_dates = np.array(zip(d, repeat(1, len(d))), dtype='double')
        dt, t = dtime(t)
        eph, err = oo.pyoorb.oorb_ephemeris(in_orbits=orbits,
                                            in_obscode=obscode,
                                            in_date_ephems=ephem_dates)
        dt, t = dtime(t)
        if singledate:
            print("Calculating ephemerides @ date %f required %f s" % (dates[i], dt))
        else:
            print("Calculating ephemerides up to date %f required %f s" % (dates[i], dt))
        tReq[i] = dt / float(len(orbits)) * scaling

    deltatime = dates - dates[0]
    pylab.plot(deltatime, tReq, 'bo')
    p = np.polyfit(deltatime, tReq, 1)
    dt = np.arange(deltatime[0], deltatime[-1], 0.1)
    y = np.polyval(p, dt)
    if singledate:
        pylab.plot(dt, y, 'b-',
                   label='Single ephemeride for %d objects @ date: y=%.1f + %.2f dt'
                   % (scaling, p[1], p[0]))
    else:
        pylab.plot(dt, y, 'b-',
                   label='Multiple ephemerides for %d objects up to date: y=%.1f + %.2f dt'
                   % (scaling, p[1], p[0]))
    pylab.legend(loc='lower right', fontsize='smaller', fancybox=True)
    pylab.xlabel("Days from epoch of orbit")
    pylab.ylabel("Approx time required to calculate ephemerides for 1000 objects")
    pylab.title("OpenOrb ephemeride calculation requirements")
    pylab.savefig('timing_test.png', format='png')
    # pylab.show()
