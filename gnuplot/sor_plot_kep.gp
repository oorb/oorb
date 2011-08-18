#====================================================================#
#                                                                    #
# Copyright 2002,2003,2004,2005,2006,2007,2008,2009,2010,2011        #
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
# function in Keplerian phase space.
#
# Author:  MG
# Version: 2011-08-18
#
reset
set terminal postscript eps enhanced linewidth 1.0 12.0
set out 'sor_results.eps'
set logscale x
unset key
set pointsize 0.3
set size 1.0,1.0
set origin 0.0,0.0
set multiplot
set lmargin 11
set rmargin 2
set tmargin 1
set bmargin 3
set size 0.5,0.33
set origin 0.0,0.66
set xlabel 'a [AU]'
set ylabel 'e'
plot 'sor_orbits.out' using 1:2 pt 7
set size 0.5,0.33
set origin 0.5,0.66
set xlabel 'a [AU]'
set ylabel '{/Symbol i} [deg]'
plot 'sor_orbits.out' using 1:3 pt 7
set size 0.5,0.33
set origin 0.0,0.33
set xlabel 'a [AU]'
set ylabel '{/Symbol O} [deg]'
plot 'sor_orbits.out' using 1:4 pt 7
set size 0.5,0.33
set origin 0.5,0.33
set xlabel 'a [AU]'
set ylabel '{/Symbol o} [deg]'
plot 'sor_orbits.out' using 1:5 pt 7
set size 0.5,0.33
set origin 0.0,0.0
set xlabel 'a [AU]'
set ylabel 'M [deg]'
plot 'sor_orbits.out' using 1:6 pt 7
set logscale y
set format y '10^{%L}'
set size 0.5,0.33
set origin 0.5,0.0
set xlabel 'a [AU]'
set ylabel 'Unnormalized discrete p.d.f.'
plot 'sor_orbits.out' using 1:7 pt 7
unset multiplot
reset
