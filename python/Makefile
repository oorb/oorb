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
include ../make.config
include ../Makefile.include

all: liboorb.$(LIBEXT)

liboorb.o: liboorb.f90
	$(FC) -c $(FCOPTIONS) $(FC_INC)../$(CLASSPATH) \
	$(FC_INC)../$(MODULEPATH) $(FC_INC)../main liboorb.f90

wrapper.o: wrapper.f90 liboorb.o
	$(FC) -c $(FCOPTIONS) $(FC_INC)../$(CLASSPATH) \
	$(FC_INC)../$(MODULEPATH) \
	$(FC_INC)../main \
	wrapper.f90

wrapper4mops.o: wrapper4mops.f90 liboorb.o
	$(FC) -c $(FCOPTIONS) $(FC_INC)../$(CLASSPATH) \
	$(FC_INC)../$(MODULEPATH) \
	$(FC_INC)../main \
	wrapper4mops.f90


liboorb.$(LIBEXT): modules classes main wrapper.o liboorb.o
	$(FC) $(FCOPTIONS) -shared -o liboorb.$(LIBEXT) \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	liboorb.o wrapper.o \
	../main/io.o \
	../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o

liboorb4mops.$(LIBEXT): modules classes main wrapper4mops.o liboorb.o
	$(FC) $(FCOPTIONS) -shared -o liboorb.$(LIBEXT) \
	$(FC_INC)../$(CLASSPATH) $(FC_INC)../$(MODULEPATH) \
	liboorb.o wrapper4mops.o \
	../main/io.o \
	../$(CLASSPATH)/*.o \
	../$(MODULEPATH)/*.o

mops: liboorb4mops.$(LIBEXT)

main:
	cd ../main/ ; make io.o

classes:
	cd ../$(CLASSPATH) ; $(MAKE) all

modules:
	cd ../$(MODULEPATH) ; $(MAKE) all


# Remove library and modules:
clean:
	rm -f *~ *.mod *.o *.$(LIBEXT) *.pyc *.pyc.$(OS)

