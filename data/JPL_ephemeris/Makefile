#====================================================================#
#                                                                    #
# Copyright 2002-2014,2015                                           #
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
# This is a makefile for Jet Propulsion Laboratory's planetary ephemeris.
#
# Author:  MG
# Version: 2015-11-04

include ../../make.config
include ../../Makefile.include

# TYPE OF EPHEMERIS [405|406|430|431]
EPH_TYPE ?= 430
# FTP client
#FTP_CLIENT = wget -nv -c -O
FTP_CLIENT = curl -L -s -C - -f --output
# Output name of the program:
EPH_BIN = de$(EPH_TYPE).dat
# Name of ascii-to-binary converter program source file:
CONV_SRC = asc2eph.f90
# Name of ascii-to-binary converter program:
CONVERTER = asc2eph
# Name of test program source file:
TESTER_SRC = ephtester.f90
# Name of test program:
TESTER = ephtester
# ASCII versions of de405 ephemeris files:
EPH_ASCII_DE405 = header.405 ascp1600.405 ascp1620.405 ascp1640.405 ascp1660.405 ascp1680.405 \
                       ascp1700.405 ascp1720.405 ascp1740.405 ascp1760.405 ascp1780.405 \
                       ascp1800.405 ascp1820.405 ascp1840.405 ascp1860.405 ascp1880.405 \
                       ascp1900.405 ascp1920.405 ascp1940.405 ascp1960.405 ascp1980.405 \
                       ascp2000.405 ascp2020.405 ascp2040.405 ascp2060.405 ascp2080.405 \
                       ascp2100.405 ascp2120.405 ascp2140.405 ascp2160.405 ascp2180.405 \
                       ascp2200.405
# ASCII versions of de406 ephemeris files:
EPH_ASCII_DE406 = header.406 ascm3000.406 ascm2900.406 ascm2800.406 ascm2700.406 ascm2600.406 \
                             ascm2500.406 ascm2400.406 ascm2300.406 ascm2200.406 ascm2100.406 \
                             ascm2000.406 ascm1900.406 ascm1800.406 ascm1700.406 ascm1600.406 \
                             ascm1500.406 ascm1400.406 ascm1300.406 ascm1200.406 ascm1100.406 \
                             ascm1000.406 ascm0900.406 ascm0800.406 ascm0700.406 ascm0600.406 \
                             ascm0500.406 ascm0400.406 ascm0300.406 ascm0200.406 ascm0100.406 \
                             ascp0000.406 ascp0100.406 ascp0200.406 ascp0300.406 ascp0400.406 \
                             ascp0500.406 ascp0600.406 ascp0700.406 ascp0800.406 ascp0900.406 \
                             ascp1000.406 ascp1100.406 ascp1200.406 ascp1300.406 ascp1400.406 \
                             ascp1500.406 ascp1600.406 ascp1700.406 ascp1800.406 ascp1900.406 \
                             ascp2000.406 ascp2100.406 ascp2200.406 ascp2300.406 ascp2400.406 \
                             ascp2500.406 ascp2600.406 ascp2700.406 ascp2800.406 ascp2900.406
# ASCII versions of de430 ephemeris files:
EPH_ASCII_DE430 = header.430_229 ascp1550.430 ascp1650.430 ascp1750.430 ascp1850.430 \
                                 ascp1950.430 ascp2050.430 ascp2150.430 ascp2250.430 \
                                 ascp2350.430 ascp2450.430 ascp2550.430
# ASCII versions of de431 ephemeris files:
EPH_ASCII_DE431 = header.431_229 ascm13000.431 ascm12000.431 ascm11000.431 ascm10000.431 \
                                 ascm09000.431 ascm08000.431 ascm07000.431 ascm06000.431 \
                                 ascm05000.431 ascm04000.431 ascm03000.431 ascm02000.431 \
                                 ascm01000.431 ascp00000.431 ascp01000.431 ascp02000.431 \
                                 ascp03000.431 ascp04000.431 ascp05000.431 ascp06000.431 \
                                 ascp07000.431 ascp08000.431 ascp09000.431 ascp10000.431 \
                                 ascp11000.431 ascp12000.431 ascp13000.431 ascp14000.431 \
                                 ascp15000.431 ascp16000.431
# Input file for ascii-to-binary converter:
EPH_INPUT  = ascii.$(EPH_TYPE)
# Project name:
PROJNAME = JPL_Ephemeris

# Compile converter, build input file and convert ascii ephemeris files to a single binary file:
.PHONY: eph
eph : $(EPH_BIN)
	cp $(EPH_BIN) ..

$(EPH_BIN) : $(CONVERTER) $(EPH_INPUT)
	rm -f $(EPH_BIN)
	DYLD_LIBRARY_PATH="../lib:$$CONDA_PREFIX/lib" ./$(CONVERTER) --eph-type=$(EPH_TYPE) < $(EPH_INPUT)

test : $(TESTER) $(EPH_BIN) testpo.$(EPH_TYPE)
	./$(TESTER) --eph-type=$(EPH_TYPE)

# If the library hasn't been built yet, here's the rule to build it
../../lib/liboorb.a:
	$(MAKE) -C ../.. liboorb

# Compile ascii-to-binary converter:
$(CONVERTER): $(CONV_SRC) ../../lib/liboorb.a
	$(FC) -o $(CONVERTER) $(FC_INC)../../build $(CONV_SRC) ../../lib/liboorb.a $(ADDLIBS) 

# Compile tester:
$(TESTER): $(TEST_SRC) ../../lib/liboorb.a
	$(FC) -o $(TESTER) $(FC_INC)../../build $(TESTER_SRC) ../../lib/liboorb.a $(ADDLIBS) 

# Download individual ephemeris files (via a temp file, so that partial downloads aren't
# seen as successful by make
$(EPH_ASCII_DE$(EPH_TYPE)) testpo.$(EPH_TYPE):
	$(FTP_CLIENT) $@.tmp ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de$(EPH_TYPE)/$@ && mv $@.tmp $@

# Build input file for ascii-to-binary converter:
$(EPH_INPUT): $(EPH_ASCII_DE$(EPH_TYPE)) testpo.$(EPH_TYPE)
	cat $(EPH_ASCII_DE$(EPH_TYPE)) > $(EPH_INPUT)

# Write back-up:
backup:
	$(SHELL) -c "if test -d ../backup_$(PROJNAME); then true; else mkdir ../backup_$(PROJNAME); fi"
	cp -a * ../backup_$(PROJNAME)

# Make compressed tar-ball:
tar: clean
	rm -f $(PROJNAME)-*.tgz
	tar cvzf $(PROJNAME)-$(EPH_TYPE).tgz *

# Clean up directory:
clean:
	rm -f *.o *~ core $(CONVERTER) $(TESTER) *.$(EPH_TYPE) \
              $(EPH_ASCII_DE$(EPH_TYPE)) header.$(EPH_TYPE)* \
              $(EPH_BIN)

