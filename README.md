# OpenOrb (or OOrb), an open-source orbit-computation package.

[![CI](https://github.com/oorb/oorb/workflows/CI/badge.svg)]()
[![Test Coverage](https://codecov.io/gh/mjuric/oorb/branch/master/graph/badge.svg)](https://codecov.io/gh/mjuric/oorb)


# Introduction #

OOrb contains, e.g., the statistical orbital ranging method (hereafter
referred to as Ranging). Ranging is used to solve the orbital inverse
problem of computing non-Gaussian orbital-element probability density
functions (p.d.f.s) based on input astrometry.

Ranging is optimized for cases where the amount of astrometry is
scarce or spans a relatively short time span. Ranging-based methods
have successfully been applied to a variety of different topics such
as rigorous ephemeris prediction, orbital-element-distribution studies
for trans-neptunian objects, the computation of invariant collision
probabilities between NEOs and the Earth, detecting linkages between
astrometric asteroid observations within an apparition as well as
between apparitions, and in the rigorous analysis of the impact of
orbital arc-length and/or astrometric uncertainty on the uncertainty
of the resulting orbits.

In OOrb, tools for making ephemeris predictions and classification of
objects (i.e., NEO-MBO-TNO) are also available.


Documentation on usage and installation is available on the [oorb wiki](https://github.com/oorb/oorb/wiki).


# Citation information #

When using this software please cite

	Granvik, M., Virtanen, J., Oszkiewicz, D., Muinonen, K. (2009).
	OpenOrb: Open-source asteroid orbit computation software including statistical ranging.
	Meteoritics & Planetary Science 44(12), 1853-1861.


# Licensing information #

OpenOrb is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenOrb is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should receive a copy of the GNU General Public License along
with OpenOrb. If not, see <http://www.gnu.org/licenses/>.
