#
# LSST Data Management System
# Copyright 2008, 2009 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
# F. Pierfederici <fpierfed@gmail.com>
import os
import ctypes
import ctypes.util
import numpy

from constants import *



# Find the oorb library
libPath = ctypes.util.find_library('oorb')
if(not libPath):
    raise(Exception('Could not find liboorb!'))
ffi = ctypes.CDLL(libPath)





def init(ephemFile, verbosity=0):
    """
    Init the oorb module.
    
    @param ephemFile: full path of the JPL ephemeris de405.dat file
    @param verbosity: integer indicating verbosity level. 0=quiet.
    
    Return
    None
    """
    # Init the error code.
    err = ctypes.c_int(0)
    
    # Do a simple sanity check in the input.
    if(not ephemFile):
        raise(Exception('Bad input.'))
    if(verbosity < 0):
        verbosity = 0
    elif(verbosity > 5):
        verbosity = 5
    
    # Invoke the fortran init function.
    ffi.oorb_init_(ctypes.c_char_p(ephemFile), 
                   ctypes.byref(ctypes.c_int(verbosity)), 
                   ctypes.byref(err),
                   ctypes.c_int(len(ephemFile)))
    if(err.value):
        raise(Exception('Error: oorb_init failed (error=%d)' \
                        %(err.value)))
    return



def calendardate_to_mjd(y, m, d, timescale='TT'):
    """
    Convert a calendar date expressed as year, month, day.fraction to MJD.
    timescale is either UTC, UT1 or TT. Output MJD is TAI
    
    @param y: year number
    @param m: month number
    @param d: day number
    @param timescale: either "UTC", "UT1", "TT" or "TAI"
    
    Return
    MJD TAI
    """
    # Init.
    mjd = ctypes.c_double(0.)
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(y < MIN_YEAR or m < 0 or m > 12 or d < 0 or d >= 32 or not timescale):
        raise(Exception('Bad input.'))
    
    # Enter Fortran land.
    ffi.oorb_calendardate_to_mjd_(ctypes.byref(ctypes.c_int(y)), 
                                  ctypes.byref(ctypes.c_int(m)), 
                                  ctypes.byref(ctypes.c_double(d)), 
                                  ctypes.byref(mjd), 
                                  ctypes.c_char_p(timescale), 
                                  ctypes.byref(err), 
                                  ctypes.c_int(len(timescale)))
    if(err.value):
        raise(Exception('Error: oorb_calendardate_to_mjd failed (error=%d)' \
                        %(err.value)))
    return(mjd.value)



def mjdutc_to_mjdtai(mjdIn):
    """
    Convert a calendar date expressed as MJD UTC to MJD TAI.
    
    @param mjdIn: MJD UTC
    
    Return
    MJD TAI
    """
    return(_convertTimescale(mjdIn, fctn=ffi.oorb_mjdutc_to_mjdtai_))



def mjdtai_to_mjdutc(mjdIn):
    """
    Convert a calendar date expressed as MJD TAI to MJD UTC.
    
    @param mjdIn: MJD TAI
    
    Return
    MJD UTC
    """
    return(_convertTimescale(mjdIn, fctn=ffi.oorb_mjdtai_to_mjdutc_))



def mjdtt_to_mjdtai(mjdIn):
    """
    Convert a calendar date expressed as MJD TT to MJD TAI. Remember that, in 
    first approximation, TT = TAI + 32.184 seconds.
    
    @param mjdIn: MJD TT
    
    Return
    MJD TAI
    """
    return(_convertTimescale(mjdIn, fctn=ffi.oorb_mjdtt_to_mjdtai_))



def mjdtai_to_mjdtt(mjdIn):
    """
    Convert a calendar date expressed as MJD TAI to MJD TT. Remember that, in 
    first approximation, TT = TAI + 32.184 seconds.
    
    @param mjdIn: MJD TAI
    
    Return
    MJD TT
    """
    return(_convertTimescale(mjdIn, fctn=ffi.oorb_mjdtai_to_mjdtt_))



def ranging_fast(trackId, coords, mjds, mags, obscodes, filters, elementType, 
                 numOrbits):
    """
    Statistical ranging. Given a track (meaning an ID, coordinates, times, mags,
    filters, obscodes), produce an orbit cloud of the given size.
    
    Use a 2-body dynamical model.
    
    
    @param trackId: track ID. Just a number, really.
    @param coords: RA, Dec array. RA and Dec in decimal degrees.
    @param mjds: MJD TAI array
    @param mags: mag array
    @param obscodes: obscode array. Obscodes are strings of length 4
    @param filters: filter array. Filters are strings of length 1
    @param elementType: desired type of orbital elements: keplerian or cartesian
    @param numOrbits: number of ranging orbits to generate.
    
    Return
    List of numOrbits flattened orbits. Each orbits is a 15-element array of the
    form
        [trackId, *elements, epochMJD, Un-normalized p.d.f., Reduced chi2, 
         Regularized apr, *Jacobian, elementTypeIdx]
    where
        elements is the 6-element orbital element array.
        Jacobian is the 3-element Jacobian array (diagonal of the 3x3 matrix).
    """
    return(_ranging(trackId, 
                    coords, 
                    mjds, 
                    mags, 
                    obscodes, 
                    filters, 
                    elementType, 
                    numOrbits,
                    fctn=ffi.oorb_ranging_fast_))



def ranging(trackId, coords, mjds, mags, obscodes, filters, elementType, 
            numOrbits):
    """
    Statistical ranging. Given a track (meaning an ID, coordinates, times, mags,
    filters, obscodes), produce an orbit cloud of the given size.
    
    Use the full n-body dynamical model.
    
    
    @param trackId: track ID. Just a number, really.
    @param coords: RA, Dec array. RA and Dec in decimal degrees.
    @param mjds: MJD TAI array
    @param mags: mag array
    @param obscodes: obscode array. Obscodes are strings of length 4
    @param filters: filter array. Filters are strings of length 1
    @param elementType: desired type of orbital elements: keplerian or cartesian
    @param numOrbits: number of ranging orbits to generate.
    
    Return
    List of numOrbits flattened orbits. Each orbits is a 15-element array of the
    form
        [trackId, *elements, epochMJD, Un-normalized p.d.f., Reduced chi2, 
         Regularized apr, *Jacobian, elementTypeIdx]
    where
        elements is the 6-element orbital element array.
        Jacobian is the 3-element Jacobian array (diagonal of the 3x3 matrix).
    """
    return(_ranging(trackId, 
                    coords, 
                    mjds, 
                    mags, 
                    obscodes, 
                    filters, 
                    elementType, 
                    numOrbits,
                    fctn=ffi.oorb_ranging_))



def lsl_fast(trackId, coords, mjds, mags, obscodes, filters, rangingOrbits):
    """
    Orbital inversion using least squares with linearized covariances, that is, 
    fixing the resulting shape of the orbital-element pdf to a Gaussian. 
    
    Given a track (meaning an ID, coordinates, times, mags, filters, obscodes) 
    and a cloud of ranging orbits, produce a single orbit. The orbit cloud is a
    list of lists of the same form as the output of the ranging routines.
    
    Use a 2-body dynamical model.
    
    
    @param trackId: track ID. Just a number, really.
    @param coords: RA, Dec array. RA and Dec in decimal degrees.
    @param mjds: MJD TAI array
    @param mags: mag array
    @param obscodes: obscode array. Obscodes are strings of length 4
    @param filters: filter array. Filters are strings of length 1
    @param rangingOrbits: orbit cloud generated by statistical ranging.
    
    Return
    Tuple of the form
        (orb, cov)
    where
        orb is a single flatted orbit of the form
            [track_id, *elements, epoch, H, G, el_type_index]
        cov is the 6x6 covariance matrix
        elements is the 6-element orbital element array of the same type as the
            ranging orbit cloud elements.
        el_type_index is an internal index. Leave it alone, unless you know what
            you are doing.
    """
    return(_lsl(trackId, 
                coords, 
                mjds, 
                mags, 
                obscodes, 
                filters, 
                rangingOrbits,
                fctn=ffi.oorb_lsl_fast_))



def lsl(trackId, coords, mjds, mags, obscodes, filters, rangingOrbits):
    """
    Orbital inversion using least squares with linearized covariances, that is, 
    fixing the resulting shape of the orbital-element pdf to a Gaussian. 
    
    Given a track (meaning an ID, coordinates, times, mags, filters, obscodes) 
    and a cloud of ranging orbits, produce a single orbit. The orbit cloud is a
    list of lists of the same form as the output of the ranging routines.
    
    Use the full n-body dynamical model.
    
    
    @param trackId: track ID. Just a number, really.
    @param coords: RA, Dec array. RA and Dec in decimal degrees.
    @param mjds: MJD TAI array
    @param mags: mag array
    @param obscodes: obscode array. Obscodes are strings of length 4
    @param filters: filter array. Filters are strings of length 1
    @param rangingOrbits: orbit cloud generated by statistical ranging.
    
    Return
    Tuple of the form
        (orb, cov)
    where
        orb is a single flatted orbit of the form
            [track_id, *elements, epoch, H, G, el_type_index]
        cov is the 6x6 covariance matrix
        elements is the 6-element orbital element array of the same type as the
            ranging orbit cloud elements.
        el_type_index is an internal index. Leave it alone, unless you know what
            you are doing.
    """
    return(_lsl(trackId, 
                coords, 
                mjds, 
                mags, 
                obscodes, 
                filters, 
                rangingOrbits,
                fctn=ffi.oorb_lsl_))



def propagate_orbit_fast(orb, cov, mjd):
    """
    Propagate a single orbit to a given MJD TAI.
    
    Use a 2-body dynamical model.
    
    
    @param orb: flattened orbit in the same format as the output of lsl
    @param cov: covariance 6x6 matrix (same format as lsl)
    @param mjd: MJD TAI to propagate the input orbit to
    
    Return
    Tuple of the form
        (newOrb, newCov)
    where 
        newOrb is the new propagated orbit
        newCov is the new propagated 6x6 covariance matrix
    """
    return(_propagate_orbit(orb, 
                            cov, 
                            mjd,
                            fctn=ffi.oorb_propagate_orbit_fast_))



def propagate_orbit(orb, cov, mjd):
    """
    Propagate a single orbit to a given MJD TAI.
    
    Use the full n-body dynamical model.
    
    
    @param orb: flattened orbit in the same format as the output of lsl
    @param cov: covariance 6x6 matrix (same format as lsl)
    @param mjd: MJD TAI to propagate the input orbit to
    
    Return
    Tuple of the form
        (newOrb, newCov)
    where 
        newOrb is the new propagated orbit
        newCov is the new propagated 6x6 covariance matrix
    """
    return(_propagate_orbit(orb, 
                            cov, 
                            mjd,
                            fctn=ffi.oorb_propagate_orbit_))



def ephemeris_fast(orb, cov, obscode, n, step):
    """
    Compute n ephemeris step (fractional) days apart for the given orbit and 
    covariance matrix at the orbit epoch. To compute ephemeris from an arbitrary
    MJD, propagate the orbit to that MJD first.
    
    Use a 2-body dynamical model.
    
    
    @param orb: flattened orbit in the same format as the output of lsl
    @param cov: covariance 6x6 matrix (same format as lsl)
    @param obscode: 4-letter MPC observatory code
    @param n: number of ephemeris to compute
    @param step: time separation in fractional days for successive ephemeris
    
    Return
    A list of n ephemeris. Each ephemeris is a list of the form
        [dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa]
    where
        dist is the distance from Earth in AU
        ra is RA in decimal degrees
        dec is Dec in decimal degrees
        mag is the predicted magnitude (NOT COMPUTED YET)
        mjd is the MJD TAI of the ephemeris
        raErr is the error on RA (NOT COMPUTED YET)
        decErr is the error on Dec (NOT COMPUTED YET)
        smaa is the positional error ellipse semi major axis (NOT COMPUTED YET)
        smia is the positional error ellipse semi minor axis (NOT COMPUTED YET)
        pa is the positional error ellipse position angle (NOT COMPUTED YET)
    """
    return(_ephemeris(orb, 
                      cov, 
                      obscode, 
                      n, 
                      step, 
                      fctn=ffi.oorb_ephemeris_fast_))



def ephemeris(orb, cov, obscode, n, step):
    """
    Compute n ephemeris step (fractional) days apart for the given orbit and 
    covariance matrix at the orbit epoch. To compute ephemeris from an arbitrary
    MJD, propagate the orbit to that MJD first.
    
    Use the full n-body dynamical model.
    
    
    @param orb: flattened orbit in the same format as the output of lsl
    @param cov: covariance 6x6 matrix (same format as lsl)
    @param obscode: 4-letter MPC observatory code
    @param n: number of ephemeris to compute
    @param step: time separation in fractional days for successive ephemeris
    
    Return
    A list of n ephemeris. Each ephemeris is a list of the form
        [dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa]
    where
        dist is the distance from Earth in AU
        ra is RA in decimal degrees
        dec is Dec in decimal degrees
        mag is the predicted magnitude (NOT COMPUTED YET)
        mjd is the MJD TAI of the ephemeris
        raErr is the error on RA (NOT COMPUTED YET)
        decErr is the error on Dec (NOT COMPUTED YET)
        smaa is the positional error ellipse semi major axis (NOT COMPUTED YET)
        smia is the positional error ellipse semi minor axis (NOT COMPUTED YET)
        pa is the positional error ellipse position angle (NOT COMPUTED YET)
    """
    return(_ephemeris(orb, 
                      cov, 
                      obscode, 
                      n, 
                      step, 
                      fctn=ffi.oorb_ephemeris_))



def moid(orb, cov):
    """
    Compute the MOID given an orbit and its covariance matrix.
    
    
    @param orb: flattened orbit in the same format as the output of lsl
    @param cov: covariance 6x6 matrix (same format as lsl)
    
    Return
    moid
    """
    # Init.
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(orb == None or cov == None):
        raise(Exception('Bad input.'))
    
    # Allocate space for the output moid: a simple float.
    moid = ctypes.c_double(-1.)
    
    # Enter Fortran land.
    ffi.oorb_moid_(orb.ctypes.data_as(ctypes.c_void_p),
                   cov.ctypes.data_as(ctypes.c_void_p),
                   ctypes.byref(moid),
                   ctypes.byref(err))
    if(err.value):
        raise(Exception('Error: oorb_moid failed (error=%d)' \
                        %(err.value)))
    return(moid.value)



def classification(rangingOrbits):
    """
    Compute the likelihood that a given ranging orbit cloud could represent a
    given asteroid class.
    
    
    @param rangingOrbits: orbit cloud generated by statistical ranging.
    
    Return
    Tuple of the form
        (asteroid_classes, likelihoods)
    where
        asteroid_classes is a list of asteroid class names
        likelihoods is a list of probability values
    both lists are ordered so that likelihoods[j] corresponds to 
    asteroid_classes[j]
    """
    # Init.
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(rangingOrbits == None):
        raise(Exception('Bad input.'))
    
    # Allocate space for the output weights: a 1-d array of 17 doubles.
    weights = numpy.zeros(shape=(17,), dtype='d', order='Fortran')
    
    # Enter Fortran land.
    ffi.priv_classification_(rangingOrbits.ctypes.data_as(ctypes.c_void_p),
                             ctypes.byref(ctypes.c_int(len(rangingOrbits))),
                             weights.ctypes.data_as(ctypes.c_void_p),
                             ctypes.byref(err))
    return(ASTEROID_CLASSES, weights)



# Internal routines.
def _ranging(trackId, coords, mjds, mags, obscodes, filters, elementType, 
            numOrbits, fctn):
    # Init.
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(coords == None or obscodes == None or filters == None or 
       elementType == None):
        raise(Exception('Bad input.'))
    
    # Make sure that each obscode element is indeed 4-char long and that filters
    # are 1-char long.
    obscodes = ['%-4s' %(o[:4]) for o in obscodes]
    filters = ['%s' %(f[:1]) for f in filters]
    # And make sure that elementType is 11 char long.
    elementType = '%-11s' %(elementType[:11])
    
    # Allocate space for the output orbits: an array of the form
    # [[id, elements(0:5), epoch, Un-normalized p.d.f., Reduced chi2, 
    #   Regularized apr, Jacobian det(0:2), elementTypeIdx], ] with numOrbits elements.
    rawOrbits = numpy.zeros(shape=(numOrbits, 15), dtype='d', order='Fortran')
    # orbits = numpy.zeros(shape=(15, numOrbits), dtype='d')
    mjds = numpy.array(mjds, dtype='d', order='Fortran')
    mags = numpy.array(mags, dtype='d', order='Fortran')
    coords = numpy.array(coords, dtype='d', order='Fortran')
    numCoords = len(mjds)
    
    # Enter Fortran land.
    fctn(ctypes.c_char_p(''.join(obscodes)),
         ctypes.c_char_p(''.join(filters)),
         ctypes.byref(ctypes.c_int(trackId)),
         ctypes.byref(ctypes.c_int(numCoords)),
         coords.ctypes.data_as(ctypes.c_void_p),
         mjds.ctypes.data_as(ctypes.c_void_p),
         mags.ctypes.data_as(ctypes.c_void_p),
         ctypes.c_char_p(elementType),
         ctypes.byref(ctypes.c_int(numOrbits)),
         rawOrbits.ctypes.data_as(ctypes.c_void_p),
         ctypes.byref(err),
         ctypes.c_int(len(elementType)))
    if(err.value):
        raise(Exception('Error: %s failed (error=%d)' \
              %(fctn.__name__[5:-1], err.value)))
    
    # Now transpose the array and create a matrix.
#     for j in range(numOrbits):
#         for k in range(15):
#             orbits[j][k]  = rawOrbits[numOrbits*k + j]
    return(rawOrbits)


    
def _lsl(trackId, coords, mjds, mags, obscodes, filters, rangingOrbits, fctn):
    # Init.
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(coords == None or obscodes == None or filters == None or 
       rangingOrbits == None):
        raise(Exception('Bad input.'))
    
    # Make sure that each obscode element is indeed 4-char long and that filters
    # are 1-char long.
    obscodes = ['%-4s' %(o[:4]) for o in obscodes]
    filters = ['%-1s' %(f[:1]) for f in filters]
    
    # Allocate space for the output orbit: a 1-d array of the form
    # [track_id, elements(1:6), epoch, H, G, el_type_index]
    orbit = numpy.zeros(shape=(11,), dtype='d', order='Fortran')
    # covariance is a 6x6 matrix of doubles:
    cov = numpy.zeros(shape=(6, 6), dtype='d', order='Fortran')
    sigmas = numpy.zeros(shape=(6,), dtype='d', order='Fortran')
    corr = numpy.zeros(shape=(6, 6), dtype='d', order='Fortran')
    mjds = numpy.array(mjds, dtype='d', order='Fortran')
    mags = numpy.array(mags, dtype='d', order='Fortran')
    coords = numpy.array(coords, dtype='d', order='Fortran')
    numCoords = len(mjds)
    numOrbits = len(rangingOrbits)
    
    # Enter Fortran land.
    fctn(ctypes.c_char_p(''.join(obscodes)),
         ctypes.c_char_p(''.join(filters)),
         ctypes.byref(ctypes.c_int(trackId)),
         ctypes.byref(ctypes.c_int(numCoords)),
         coords.ctypes.data_as(ctypes.c_void_p),
         mjds.ctypes.data_as(ctypes.c_void_p),
         mags.ctypes.data_as(ctypes.c_void_p),
         rangingOrbits.ctypes.data_as(ctypes.c_void_p),
         ctypes.byref(ctypes.c_int(numOrbits)),
         orbit.ctypes.data_as(ctypes.c_void_p),
         cov.ctypes.data_as(ctypes.c_void_p),
         sigmas.ctypes.data_as(ctypes.c_void_p),
         corr.ctypes.data_as(ctypes.c_void_p),
         ctypes.byref(err))
    if(err.value):
        raise(Exception('Error: %s failed (error=%d)' \
              %(fctn.__name__[5:-1], err.value)))
    return(orbit, cov)



def _propagate_orbit(orb, cov, mjd, fctn):
    # Init.
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(orb == None or cov == None or not mjd):
        raise(Exception('Bad input.'))
    
    # Allocate space for the output orbit: a 1-d array of the form
    # [track_id, elements(1:6), epoch, H, G, el_type_index]
    newOrbit = numpy.zeros(shape=(11,), dtype='d', order='Fortran')
    # covariance is a 6x6 matrix of doubles:
    newCov = numpy.zeros(shape=(6, 6), dtype='d', order='Fortran')
    newSigmas = numpy.zeros(shape=(6,), dtype='d', order='Fortran')
    newCorr = numpy.zeros(shape=(6, 6), dtype='d', order='Fortran')
    
    # Enter Fortran land.
    fctn(orb.ctypes.data_as(ctypes.c_void_p),
         cov.ctypes.data_as(ctypes.c_void_p),
         ctypes.byref(ctypes.c_double(mjd)),
         newOrbit.ctypes.data_as(ctypes.c_void_p),
         newCov.ctypes.data_as(ctypes.c_void_p),
         newSigmas.ctypes.data_as(ctypes.c_void_p),
         newCorr.ctypes.data_as(ctypes.c_void_p),
         ctypes.byref(err))
    if(err.value):
        raise(Exception('Error: %s failed (error=%d)' \
              %(fctn.__name__[5:-1], err.value)))
    return(newOrbit, newCov)



def _ephemeris(orb, cov, obscode, n, step, fctn):
    # Init.
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(orb == None or cov == None or not n or not step):
        raise(Exception('Bad input.'))
    
    # Allocate space for the output ephems: a n-d array of the form
    # [[dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa], ]
    ephems = numpy.zeros(shape=(n, 10), dtype='d', order='Fortran')
    
    # Make sure obscode is 4-char long.
    obscode = '%-4s' %(obscode[:4])
    
    # Enter Fortran land.
    fctn(orb.ctypes.data_as(ctypes.c_void_p),
         cov.ctypes.data_as(ctypes.c_void_p),
         ctypes.c_char_p(obscode),
         ctypes.byref(ctypes.c_int(n)),
         ctypes.byref(ctypes.c_double(step)),
         ephems.ctypes.data_as(ctypes.c_void_p),
         ctypes.byref(err),
         ctypes.c_int(len(obscode)))
    if(err.value):
        raise(Exception('Error: %s failed (error=%d)' \
              %(fctn.__name__[5:-1], err.value)))
    return(ephems)



def _convertTimescale(mjdIn, fctn):
    # Init.
    mjdOut = ctypes.c_double(-1.)
    err = ctypes.c_int(0)
    
    # Simple sanity check.
    if(not mjdIn):
        raise(Exception('Bad input.'))
    
    # Enter Fortran land.
    fctn(ctypes.byref(ctypes.c_double(mjdIn)), 
         ctypes.byref(mjdOut), 
         ctypes.byref(err))
    if(err.value):
        raise(Exception('Error: oorb_calendardate_to_mjd failed (error=%d)' \
                        %(err.value)))
    return(mjdOut.value)








if(__name__ == '__main__'):
    import sys
    # Do not use sys.stdout.write|flush because it confuses the fortran code :-)
    # print is OK.
    
    # All times in/out need to be TAI. The MJD constants used below are UTC and
    # therefore need to be converted...
    
    sys.stdout.write('oorb.init:                        ')
    ef = os.path.join(os.environ['OORB_DATA'], 'JPL_ephemeris', 'de405.dat')
    init(ef, verbosity=2)
    sys.stdout.write('.')
    sys.stdout.write(' (all tests ok)\n')
    
    sys.stdout.write('oorb.calendardate_to_mjd:         ')
    mjd = calendardate_to_mjd(2008, 5, 31.35234, 'UTC')
    assert(mjd == mjdtt_to_mjdtai(54617.353094444443))
    sys.stdout.write('.')
    mjd = calendardate_to_mjd(2008, 5, 31.39302, 'UTC')
    assert(mjd == mjdtt_to_mjdtai(54617.393774444448))
    sys.stdout.write('.')
    mjd = calendardate_to_mjd(2008, 5, 31.43442, 'UTC')
    assert(mjd == mjdtt_to_mjdtai(54617.435174444443))
    sys.stdout.write('.')
    sys.stdout.write(' (all tests ok)\n')
    
    coords = [[253.64316666666659, 8.33333333333333306E-05, 19.381388888888889, 8.33333333333333306E-05],
              [253.64175000000000, 8.33333333333333306E-05, 19.381833333333333, 8.33333333333333306E-05],
              [253.64033333333330, 8.33333333333333306E-05, 19.382222222222218, 8.33333333333333306E-05]]
    mjds = [mjdtt_to_mjdtai(54617.353094444443), 
            mjdtt_to_mjdtai(54617.393774444448), 
            mjdtt_to_mjdtai(54617.435174444443)]
    mags = [23.699999999999999, 
            23.699999999999999, 
            23.800000000000001]
    obscodes = ['568', ] * len(mags)
    filters = ['r', ] * len(mags)
        
    
    
    try:
        rangingOrbits = ranging(0, coords, mjds, mags, obscodes, filters,
                                'keplerian', 5000)
        for i in range(10):
            print(rangingOrbits[i])
    except:
        print('ranging failed (which is normal for this data set).')
        


    rangingOrbits = ranging_fast(0, coords, mjds, mags, obscodes, filters,
                                 'keplerian', 5000)
    for i in range(10):
        print(rangingOrbits[i])
    
    
    
    
    
    coords = [[253.64316666666659, 8.3333333333333331e-05, 19.381388888888889, 8.3333333333333331e-05], 
              [253.64175000000000, 8.3333333333333331e-05, 19.381833333333333, 8.3333333333333331e-05], 
              [253.64033333333330, 8.3333333333333331e-05, 19.382222222222222, 8.3333333333333331e-05], 
              [253.37720833333330, 8.3333333333333331e-05, 19.449750000000002, 8.3333333333333331e-05], 
              [253.37575000000001, 8.3333333333333331e-05, 19.450166666666664, 8.3333333333333331e-05], 
              [253.34449999999998, 8.3333333333333331e-05, 19.456777777777777, 8.3333333333333331e-05], 
              [252.87362500000003, 8.3333333333333331e-05, 19.518361111111112, 8.3333333333333331e-05], 
              [252.87212500000004, 8.3333333333333331e-05, 19.518388888888889, 8.3333333333333331e-05], 
              [252.84200000000004, 8.3333333333333331e-05, 19.519861111111108, 8.3333333333333331e-05], 
              [252.84029166666667, 8.3333333333333331e-05, 19.519944444444445, 8.3333333333333331e-05], 
              [252.42179166666665, 8.3333333333333331e-05, 19.507138888888889, 8.3333333333333331e-05], 
              [252.42158333333333, 8.3333333333333331e-05, 19.507111111111112, 8.3333333333333331e-05], 
              [252.42124999999999, 8.3333333333333331e-05, 19.507055555555556, 8.3333333333333331e-05], 
              [252.42108333333331, 8.3333333333333331e-05, 19.507055555555556, 8.3333333333333331e-05], 
              [252.42091666666667, 8.3333333333333331e-05, 19.507027777777779, 8.3333333333333331e-05]]
    mjds = [mjdtt_to_mjdtai(54617.353094444443), 
            mjdtt_to_mjdtai(54617.393774444448), 
            mjdtt_to_mjdtai(54617.435174444443), 
            mjdtt_to_mjdtai(54625.212504444447), 
            mjdtt_to_mjdtai(54625.256664444445), 
            mjdtt_to_mjdtai(54626.181794444448), 
            mjdtt_to_mjdtai(54640.367084444442), 
            mjdtt_to_mjdtai(54640.414734444443), 
            mjdtt_to_mjdtai(54641.353904444448), 
            mjdtt_to_mjdtai(54641.406674444443), 
            mjdtt_to_mjdtai(54655.124204444444), 
            mjdtt_to_mjdtai(54655.131814444445), 
            mjdtt_to_mjdtai(54655.145554444447), 
            mjdtt_to_mjdtai(54655.150364444446), 
            mjdtt_to_mjdtai(54655.155144444449)]
    mags = [23.699999999999999, 
            23.699999999999999, 
            23.800000000000001, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006, 
            99.900000000000006]
    obscodes = ['568',
                '568',
                '568',
                '807',
                '807',
                '807',
                '696',
                '696',
                '696',
                '696',
                '807',
                '807',
                '807',
                '807',
                '807']
    filters = ['r',
               'r',
               'r',
               '',
               '',
               '',
               '',
               '',
               '',
               '',
               '',
               '',
               '',
               '',
               '']
    
    
    
    sys.stdout.write('oorb.lsl:                         ')
    sys.stdout.write('\n')
    orb, cov = lsl(0, coords, mjds, mags, obscodes, filters, rangingOrbits)
    
    
    
    
    
    sys.stdout.write('oorb.lsl_fast:                    ')
    sys.stdout.write('\n')
    orb, cov = lsl_fast(0, coords, mjds, mags, obscodes, filters, rangingOrbits)
    # print(orb)
    # print(cov)
    
    
    
    
    
    mjd = mjdtt_to_mjdtai(54994.)
    sys.stdout.write('oorb.propagate_orbit:             ')
    sys.stdout.write('\n')
    newOrb, newCov = propagate_orbit(orb, cov, mjd)
    # print(orb)
    # print(cov)
    
    
    
    
    
    mjd = mjdtt_to_mjdtai(54994.)
    sys.stdout.write('oorb.propagate_orbit_fast:        ')
    sys.stdout.write('\n')
    newOrb, newCov = propagate_orbit_fast(orb, cov, mjd)
    # print(orb)
    # print(cov)

    
    
    
    
    
    obscode = '500'
    n = 10
    step = 1.
    sys.stdout.write('oorb.ephemeris                    ')
    sys.stdout.write('\n')
    ephems = ephemeris(newOrb, newCov, obscode, n, step)
    print(ephems)
    
    
    
    
    
    obscode = '500'
    n = 10
    step = 1.
    sys.stdout.write('oorb.ephemeris_fast               ')
    sys.stdout.write('\n')
    ephems = ephemeris_fast(newOrb, newCov, obscode, n, step)
    print(ephems)
    
    
    
    
    
    moid = moid(newOrb, newCov)
    print('moid: %.6f' %(moid))
    
    
    
    
    
    
    classes, weights = classification(rangingOrbits)
    # Write classification out.
    print("Class             Probability")
    for j in range(len(weights)):
        print(" %s  %f" %(classes[j], weights[j]))























