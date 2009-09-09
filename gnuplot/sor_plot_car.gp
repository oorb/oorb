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
# Gnuplot script for plotting the orbital-element probability-density
# function in Cartesian phase space.
#
# Author:  MG
# Version: 2009-08-11
#
reset
set terminal postscript eps enhanced linewidth 1.0 8.0
set out 'sor_results.eps'
unset key
set pointsize 0.2
set size 1.0,1.0
set origin 0.0,0.0
set multiplot
set lmargin 12
set rmargin 2
set tmargin 1
set bmargin 3
set size 0.5,0.33
set origin 0.0,0.66
set xlabel 'x [AU]'
set ylabel 'y [AU]'
plot 'sor_orbits.out' using 1:2 7
set size 0.5,0.33
set origin 0.5,0.66
set xlabel 'x [AU]'
set ylabel 'z [AU]'
plot 'sor_orbits.out' using 1:3 7
set size 0.5,0.33
set origin 0.0,0.33
set xlabel 'x [AU]'
set ylabel 'dx/dt [AU/d]'
plot 'sor_orbits.out' using 1:4 7
set size 0.5,0.33
set origin 0.5,0.33
set xlabel 'x [AU]'
set ylabel 'dy/dt [AU/d]'
plot 'sor_orbits.out' using 1:5 7
set size 0.5,0.33
set origin 0.0,0.0
set xlabel 'x [AU]'
set ylabel 'dz/dt [AU/d]'
plot 'sor_orbits.out' using 1:6 7
set size 0.5,0.33
set origin 0.5,0.0
set xlabel 'x [AU]'
set ylabel 'Unnormalized discrete p.d.f.'
set logscale y
set format y '10^{%L}'
plot 'sor_orbits.out' using 1:7 7
unset multiplot
reset
