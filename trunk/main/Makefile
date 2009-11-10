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
# Makefile for main programs.
#
# Author: MG
# Date:   2009-11-09

include ../make.config
include ../Makefile.include

IOMOD = io

all: $(ORBITBIN)

# Compile and link classes to the orbit program:
$(ORBITBIN): classes $(IOMOD).o
	$(FC) $(FCOPTIONS) $(ORBITSRC) -o $(ORBITBIN) \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	$(IOMOD).o ../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o

tico: classes
	$(FC) $(FCOPTIONS) tico.f90 -o tico \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	../$(CLASSPATH)/*.o ../$(MODULEPATH)/*.o

oorb4mops: classes $(IOMOD).o
	$(FC) $(FCOPTIONS) oorb4mops.f90 -o oorb4mops \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	$(IOMOD).o ../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o

oorb4neossat: classes $(IOMOD).o
	$(FC) $(FCOPTIONS) oorb4neossat.f90 -o oorb4neossat \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	$(IOMOD).o ../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o

oorb4research: classes $(IOMOD).o
	$(FC) $(FCOPTIONS) oorb4research.f90 -o oorb4research \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	$(IOMOD).o ../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o

$(IOMOD).o: classes $(IOMOD).f90
	$(FC) -c $(FCOPTIONS) $(FC_INC)../$(CLASSPATH) \
	$(FC_INC)../$(MODULEPATH) $(IOMOD).f90

classes:
	cd ../$(CLASSPATH) ; $(MAKE) all


modules:
	cd ../$(MODULEPATH) ; $(MAKE) all


# Remove program binaries and trash:
clean:
	rm -f *~ *.o *.mod $(ORBITBIN) \
	tico oorb4mops oorb4research oorb4neossat
