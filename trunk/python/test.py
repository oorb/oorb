#!/usr/bin/env python

import pyoorb
import numpy

if __name__ == "__main__":
    print "starting..."
    print "calling oorb_init():"
    pyoorb.pyoorb.oorb_init(error_verbosity=5, info_verbosity=5)
    #print "calling popagation_2b():"
    #pyoorb.pyoorb.propagation_2b(numpy.array([1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,0.]), 
    #                             5300, 5600)
    print "calling oorb_ephemeris with singleton orbit"
    #orb is id, 6 elements, epoch_mjd, H, G, element type index
    #keplerian appears to be element type index 3
    #orbits = numpy.array([0.,1.,2.,3.,4.,5.,6.,5373.,1.,1.,3.])
    #using the first item in PS-SSM, 1/100th density, s1.01 file.
    orbits = numpy.zeros([2,12], dtype=numpy.double, order='F')
    print orbits
    #numbers taken from test.f90
    orbits[0][0] = 1.0
    orbits[0][1] = 1.853233422926951E+00
    orbits[0][2] = 3.871022198610788E-01
    orbits[0][3] = 6.474397461946572E+00
    orbits[0][4] = 2.122138159851688E+02
    orbits[0][5] = 1.602711715213041E+02
    orbits[0][6] = 3.372535412461616E+01
    for i in range(1,7):
        orbits[0][i] = numpy.radians(orbits[0][i])
    orbits[0][7] = 3.0
    orbits[0][8] = 55200.0
    orbits[0][9] = 3.0
    orbits[0][10]= 24.541
    orbits[0][11]= .15

    orbits[1][0] = 2.0
    orbits[1][1] = 1.153233422926951E+00
    orbits[1][2] = 3.171022198610788E-01
    orbits[1][3] = 6.174397461946572E+00
    orbits[1][4] = 2.122138159851688E+02
    orbits[1][5] = 1.102711715213041E+02
    orbits[1][6] = 3.172535412461616E+01
    for i in range(1,7):
        orbits[1][i] = numpy.radians(orbits[1][i])
    orbits[1][7] = 3.0
    orbits[1][8] = 55200.0
    orbits[1][9] = 3.0
    orbits[1][10]= 24.541
    orbits[1][11]= .15
    
    print orbits[0]
    print orbits[1]
    covariances = numpy.zeros([2,6,6], order='F') #zero covariance is okay for now...

    covariances[0][0][0:6] = [ 1.59357530E-06, 4.54646999E-07, 5.59222797E-06, \
                              -3.87304158E-06, -3.54135866E-07, -4.31574921E-05 ]

    covariances[0][1][0:6] = [ 4.54646999E-07, 1.29710940E-07, 1.59544037E-06, \
                              -1.10495903E-06, -1.00707885E-07, -1.23129696E-05 ]

    covariances[0][2][0:6] = [ 5.59222797E-06, 1.59544037E-06, 1.96320645E-05, \
                                 -1.36008250E-05, -1.34577989E-06, -1.51404682E-04 ]

    covariances[0][3][0:6] = [ -3.87304158E-06, -1.10495903E-06, -1.36008250E-05, \
                                9.42766927E-06, 9.60529697E-07, 1.04844563E-04 ]

    covariances[0][4][0:6] = [ -3.54135866E-07, -1.00707885E-07, -1.34577989E-06, \
                                9.60529697E-07, 2.72099434E-06, 8.49016456E-06 ]

    covariances[0][5][0:6] = [ -4.31574921E-05, -1.23129696E-05, -1.51404682E-04, \
                                1.04844563E-04, 8.49016456E-06, 1.16925907E-03 ]

    covariances[1][0][0:6] = [ 1.59357530E-06, 4.54646999E-07, 5.59222797E-06, \
                              -3.87304158E-06, -3.54135866E-07, -4.31574921E-05 ]

    covariances[1][1][0:6] = [ 4.54646999E-07, 1.29710940E-07, 1.59544037E-06, \
                              -1.10495903E-06, -1.00707885E-07, -1.23129696E-05 ]

    covariances[1][2][0:6] = [ 5.59222797E-06, 1.59544037E-06, 1.96320645E-05, \
                                 -1.36008250E-05, -1.34577989E-06, -1.51404682E-04 ]

    covariances[1][3][0:6] = [ -3.87304158E-06, -1.10495903E-06, -1.36008250E-05, \
                                9.42766927E-06, 9.60529697E-07, 1.04844563E-04 ]

    covariances[1][4][0:6] = [ -3.54135866E-07, -1.00707885E-07, -1.34577989E-06, \
                                9.60529697E-07, 2.72099434E-06, 8.49016456E-06 ]

    covariances[1][5][0:6] = [ -4.31574921E-05, -1.23129696E-05, -1.51404682E-04, \
                                1.04844563E-04, 8.49016456E-06, 1.16925907E-03 ]
    
    obscode = "500"
    ephem_dates=numpy.zeros([2,2], dtype=numpy.double, order='F')
    ephem_dates[0][0:2] = [ 55148.0, 1.0 ]
    ephem_dates[1][0:2] = [ 55148.0, 1.0 ]
    print "calling ephemeris_static"
    eph, err = pyoorb.pyoorb.oorb_ephemeris(orbits, \
                                                obscode, \
                                                ephem_dates)
    print "succesfully returned!"
    print eph, err
#     print "calling oorb_ephemeris (plural version)"
#     orbits = numpy.array( [ orbits] )
#     covariances = numpy.array( [ covariances ] )
#     ephem_dates = numpy.array( [ ephem_dates ] )
#     eph, err = pyoorb.pyoorb.oorb_ephemeris(orbits, covariances, obscode, ephem_dates)
#     print eph,err
