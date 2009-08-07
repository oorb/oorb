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



initModel = "2-body"
fullModel = "n-body"
method = "continued fraction"
method_sw = "continued fraction"
step = 5.0
initStep = 5.0
a_max = 1000.0
a_min = 0.00465424
periapsis_max = -1.0
periapsis_min = -1.0
apoapsis_max = -1.0
apoapsis_min = -1.0
rho_max = -1.0
rho_min = -1.0
sor_type = 2
norb = 5000
norb_sw = 500
ntrial = 10000000
ntrial_sw = 200000
niter = 3
genwin_mult = 6.0
accwin_mult = 4.0
i = 0
j = 0
elementMask = True
outlierRejection = False
uniform = False
regularized = False
random_obs = False
gaussian_rho = False
pdf_ml_init = -1.0
ephemStep = 0.0
obscode = "500"
timespan = 10.0
mjd = 54994.0
moid = -1.0
sor_rho_init = [1.79769313486231571E+308, ] * 4
sor_genwin_offset = [0.0, ] * 4
stdev_arr = [0.3 * arcsec_to_deg, ] * 6
outlierMultiplier = 3.0
elementType = "keplerian"
sor_norb = 5000
sor_ntrial = 10000000
sor_norb_sw = 500
sor_ntrial_sw = 200000
correctionFactor = 0.2
niterMajorMax = 20
niterMajorMin = 2
niterMinor = 10
