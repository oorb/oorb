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
"""
OpenOrb Python Module
"""
import os

from api import *
import constants
import defaults
import utils



# Constants
__verbosity__ = 0
__version__ = '1.0a1'
__author__ = 'F. Pierfederici <fpierfed@gmail.com>'
__all__ = ('calendardate_to_mjd',
           'ranging_fast',
           'ranging',
           'lsl_fast',
           'lsl',
           'propagate_orbit_fast',
           'propagate_orbit',
           'ephemeris_fast',
           'ephemeris',
           'moid',
           'classification')


# Get the path to the data dir and infer the JPL ephemeris location.
dataDir = os.environ.get('OORB_DATA', None)
if(not dataDir):
    raise(Exception('Please point $OORB_DATA to the OpenOrb data directory.'))
ephemFilePath = os.path.join(dataDir, 'JPL_ephemeris', 'de405.dat')

# Init the module.
err = init(ephemFilePath, __verbosity__)
if(err):
    raise(Exception('OpenOrb initialization failed with error code %d' %(err)))




