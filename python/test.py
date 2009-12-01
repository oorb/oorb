#====================================================================#
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
#====================================================================#
#
# Script for testing pyoorb library using Python. The results should
# be identical to the results by test.f90 if everything is okay.
#
# @author  MG, JM
# @version 2009-12-01
#

#!/usr/bin/env python

import pyoorb
import numpy

if __name__ == "__main__":
    print "starting..."
    print "calling oorb_init():"
    pyoorb.pyoorb.oorb_init(error_verbosity=5, info_verbosity=5)
    #orb is id, 6 elements, epoch_mjd, H, G, element type index
    #keplerian appears to be element type index 3
    #orbits = numpy.array([0.,1.,2.,3.,4.,5.,6.,5373.,1.,1.,3.])
    #using the first item in PS-SSM, 1/100th density, s1.01 file.
    orbits = numpy.zeros([2,12], dtype=numpy.double, order='F')
    print orbits
    orbits[0][0] = 1.0
    orbits[0][1] = 1.853233422926951E+00
    orbits[0][2] = 3.871022198610788E-01
    orbits[0][3] = 6.474397461946572E+00
    orbits[0][4] = 2.122138159851688E+02
    orbits[0][5] = 1.602711715213041E+02
    orbits[0][6] = 3.372535412461616E+01
    for i in range(3,7):
        orbits[0][i] = numpy.radians(orbits[0][i])
    orbits[0][7] = 3.0
    orbits[0][8] = 55200.0
    orbits[0][9] = 3.0
    orbits[0][10]= 24.541
    orbits[0][11]= .15

    orbits[1][0] = 2.0
    for i in range(1,12):
        orbits[1][i] = orbits[0][i]
    
    print orbits[0]
    print orbits[1]
    covariances = numpy.zeros([2,6,6], order='F') #zero covariance is okay for now...

    covariances[0][0][:] = [ 1.59357530E-06, 4.54646999E-07, 5.59222797E-06, \
                              -3.87304158E-06, -3.54135866E-07, -4.31574921E-05 ]

    covariances[0][1][:] = [ 4.54646999E-07, 1.29710940E-07, 1.59544037E-06, \
                              -1.10495903E-06, -1.00707885E-07, -1.23129696E-05 ]

    covariances[0][2][:] = [ 5.59222797E-06, 1.59544037E-06, 1.96320645E-05, \
                                 -1.36008250E-05, -1.34577989E-06, -1.51404682E-04 ]

    covariances[0][3][:] = [ -3.87304158E-06, -1.10495903E-06, -1.36008250E-05, \
                                9.42766927E-06, 9.60529697E-07, 1.04844563E-04 ]

    covariances[0][4][:] = [ -3.54135866E-07, -1.00707885E-07, -1.34577989E-06, \
                                9.60529697E-07, 2.72099434E-06, 8.49016456E-06 ]

    covariances[0][5][:] = [ -4.31574921E-05, -1.23129696E-05, -1.51404682E-04, \
                                1.04844563E-04, 8.49016456E-06, 1.16925907E-03 ]

    for i in range(0,6):
        for j in range(0,6):
            covariances[1][i][j] = covariances[0][i][j]
    
    
    obscode = "500"
    ephem_dates=numpy.zeros([2,2], dtype=numpy.double, order='F')
    ephem_dates[0][:] = [ 55148.0, 1.0 ]
    ephem_dates[1][:] = [ 55150.0, 1.0 ]

    print "calling oorb_ephemeris with covariances"
    eph, err = pyoorb.pyoorb.oorb_ephemeris(in_orbits=orbits, \
                                                in_covariances=covariances, \
                                                in_obscode=obscode, \
                                                in_date_ephems=ephem_dates)
    print "error code:", err
    print eph[0][0][:]
    print eph[0][1][:]
    print eph[1][0][:]
    print eph[1][1][:]

    print "calling oorb_ephemeris without covariances"
    eph, err = pyoorb.pyoorb.oorb_ephemeris(in_orbits=orbits, \
                                                in_obscode=obscode, \
                                                in_date_ephems=ephem_dates)
    print "error code:", err
    print eph[0][0][:]
    print eph[0][1][:]
    print eph[1][0][:]
    print eph[1][1][:]
