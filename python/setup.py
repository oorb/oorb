#!/usr/bin/env python
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
from distutils.core import setup



if __name__ == "__main__": 
    setup(name = 'oorb', 
          description = "OpenOrb module", 
          author = "Francesco Pierfederici", 
          author_email = "fpierfed@lsst.org",
          version='1.0a1',
          package_dir={'oorb': '.'},
          packages=['oorb', ],
) 
