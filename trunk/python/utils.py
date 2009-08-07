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
from constants import *
from api import *



def sexagesimalHoursToDecimalDegrees(h, m, s):
    """
    Convert an angle given in sexagesimal hour form to decimal degrees.
    """
    return(15. * sexagesimalDegreesToDecimalDegrees(h, m, s))


def sexagesimalDegreesToDecimalDegrees(d, m, s):
    """
    Convert an angle given in sexagesimal degree form to decimal degrees.
    """
    if(d >= 0.):
        return(float(d) + float(m) / 60. + float(s) / 3600.)
    return(float(d) - float(m) / 60. - float(s) / 3600.)


def _updateConfig(config, var, val):
    if(config.has_key(var) and not isinstance(config[var], list)):
        config[var] = [config[var], val]
    elif(config.has_key(var)):
        config[var].append(val)
    else:
        config[var] = val
    return(config)


def _guessType(val):
    # remove any comments.
    val = val.split('#', 1)[0]
    val = val.strip()
    if(not val):
        return(None)
    
    # Is it an array or a simple type (including strings)? Arrays are space 
    # separated and are homogeneus collections.
    tokens = val.split()
    if(len(tokens) == 1):
        # Nope: single value.
        return(_guessSingleType(val))
    # Either a string with spaces of an array of int/float/Bool.
    # Guess the type of the first token and assume that if it is not a string, 
    # we have an array (not that safe, actually).
    # TODO: Improve array parsing.
    typedToken = _guessSingleType(tokens[0])
    if(isinstance(typedToken, str)):
        return(val)
    # We have an array!
    typedVal = [_guessSingleType(t) for t in tokens]
    # Simple sanity check.
    for v in typedVal:
        if(not isinstance(v, typedToken.__class__)):
            # Ops! We made a mistake in assuming that val was an array.
            return(val)
    return(typedVal)


def _guessSingleType(val):
    if(not val):
        return(None)
    if(val.upper() == 'T'):
        return(True)
    elif(val.upper() == 'F'):
        return(False)
    
    try:
        return(float(val))
    except:
        pass
    try:
        return(int(val))
    except:
        pass
    return(val)


def readConfigurationFile(fileName):
    config = {}
    
    f = open(fileName)
    for line in f:
        line = line.split('#', 1)[0]
        line = line.strip()
        if(not line):
            continue
        
        try:
            (var, val) = line.split(':', 1)
        except:
            raise(SyntaxError('Unable to parse line "%s"' %(line)))
        val = _guessType(val)
        config = _updateConfig(config, var, val)
    f.close()
    return(config)


def _parseMPCLine(line, defaultStdDev):
    # Remove comments and leading/trailing blank spaces.
    line = line.split('#', 1)[0]
    line = line.strip()
    if(not line or line.startswith("#")):
        return
    
    tokens = line.split()
    # FIXME: Improve MPC format parsing.
    if(len(tokens) not in (13, 11)):
        raise(NotImplementedError("FIXME: Improve MPC format parsing."))
    
    objName = tokens[0]
    if(objName.endswith('*')):
        objName = objName[:-1]
    
    obsTypeYear = tokens[1]
    if(len(obsTypeYear) >= 5):
        obsType = obsTypeYear[0]
        year = obsTypeYear[1:]
    else:
        year = obsTypeYear
    year = int(year)
    
    month = int(tokens[2])
    day = float(tokens[3])
    
    # The timescale used by MPC prior to year 1972 is UT1 and since then UTC has
    # been used.
    if(year < 1972):
        timescale = "UT1"
    else:
        timescale = "UTC"
    mjd = calendardate_to_mjd(year, month, day, timescale)
    
    h = int(tokens[4])
    m = int(tokens[5])
    s = float(tokens[6])
    ra = sexagesimalHoursToDecimalDegrees(h, m, s)
    raErr = defaultStdDev[1]
    
    d = int(tokens[7])
    m = int(tokens[8])
    s = float(tokens[9])
    dec = sexagesimalDegreesToDecimalDegrees(d, m, s)
    decErr = defaultStdDev[2]
    
    # Filters need to be 1-character long.
    if(len(tokens) == 13):
        # We have mag and filter info.
        mag = float(tokens[10])
        filter = tokens[11]
    else:
        # No mag/filter info.
        mag = 99.9
        filter = ' '

    # The obscode are the last 3 characters from the last token. We need to pad
    # them to 4 character lomg.
    obscode = tokens[-1][-3:] + ' '
    return((objName, mjd, ra, raErr, dec, decErr, mag, filter, obscode))
    
    

def readObsFile(fileName, defaultStdDev):
    """
    Parse a detection file in the (old) MPC format. If no coordinate uncertainty
    is profided use the input defaultStdDev values.
    
    @param fileName: full path of the file to parse.
    @param defaultStdDev: 6 element array with uncertainties on distance, ra, 
           dec and their time derivatives. Since we typically only measure ra 
           and dec, defaultStdDev has the form:
                [None, raErr, decErr, None, None, None]
    Return
        [(objectName, MJD, RA, RAErr, Dec, DecErr, mag, filterName, obscode), ]
    
    Angles in decimal degrees. Times in TT.
    """
    data = []
    
    f = open(fileName)
    for line in f:
        values = _parseMPCLine(line, defaultStdDev)
        if(not values):
            continue
        data.append(values)
    f.close()
    return(data)


def splitObservations(obs):
    """
    Group observations per object name.
    
    @param obs: list of the form [(objectName, MJD, RA, Dec, mag, filterName), ]
    
    Return
        [obs1, obs2, ...]
    Where
        obsi = [(objectNamei, MJD, RA, Dec, mag, filterName), ]
    """
    if(not obs):
        return(obs)
    
    refName = obs[0][0]
    currentObs = []
    separatedObs = []
    for o in obs:
        name = o[0]
        if(name == refName):
            # Another observation from the current group.
            currentObs.append(o)
            continue
        else:
            # Another group.
            separatedObs.append(currentObs)
            refName = name
            currentObs = [o, ]
    # Make sure that we got all groups.
    if(currentObs not in separatedObs):
        separatedObs.append(currentObs)
    return(separatedObs)