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
# Makefile for building the OOrb library.
#
# Author: MG
# Date:   2009-11-13

include ../make.config
include ../Makefile.include


# Make everything:
all: liboorb

liboorb: io
	$(FC) $(FCOPTIONS) $(FC_SHARED) -o liboorb.$(LIBEXT) \
	`ls ../$(MODULEPATH)/*.o ../$(CLASSPATH)/*.o ../$(MAINPATH)/*.o`

io:
	cd ../$(MAINPATH) ; make io.o

clean:
	rm -f *~ *.$(LIBEXT)
