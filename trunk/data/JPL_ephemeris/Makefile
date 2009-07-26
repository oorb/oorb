#====================================================================#
#                                                                    #
# Copyright 2009 Mikael Granvik, Jenni Virtanen, Karri Muinonen,     #
#                Teemu Laakso, Dagmara Oszkiewicz                    #
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
# Version: 2008-08-12

include ../../make.config
include ../../Makefile.include

# FTP client
FTP_CLIENT = ncftpget
# Output name of the program: 
EPH_BIN = de405.dat
# Name of ascii-to-binary converter program source file:
CONV_SRC = asc2eph.f90
# Name of ascii-to-binary converter program:
CONVERTER = asc2eph
# Name of test program source file:
TESTER_SRC = de405_test.f90
# Name of test program:
TESTER = de405_test
# ASCII versions of ephemeris files:
EPH_ASCII = header.405 ascp1600.405 ascp1620.405 ascp1640.405 ascp1660.405 ascp1680.405 \
                       ascp1700.405 ascp1720.405 ascp1740.405 ascp1760.405 ascp1780.405 \
                       ascp1800.405 ascp1820.405 ascp1840.405 ascp1860.405 ascp1880.405 \
                       ascp1900.405 ascp1920.405 ascp1940.405 ascp1960.405 ascp1980.405 \
                       ascp2000.405 ascp2020.405 ascp2040.405 ascp2060.405 ascp2080.405 \
                       ascp2100.405 ascp2120.405 ascp2140.405 ascp2160.405 ascp2180.405 \
                       ascp2200.405
# Input file for ascii-to-binary converter:
EPH_INPUT  = ascii.405 
# Project name:
PROJNAME = JPL_Ephemeris
# Version:
VERSION = 405

# Compile converter, build input file and convert ascii ephemeris files to a single binary file:
de405 : $(CONVERTER) $(EPH_INPUT)
	rm -f $(EPH_BIN)
	nice ./$(CONVERTER) < $(EPH_INPUT)
	rm -f $(EPH_INPUT)
	cp $(EPH_BIN) ../

test : $(TESTER)
	$(TESTER)

# Compile ascii-to-binary converter:
$(CONVERTER): $(CONV_SRC)
	$(FC) -o $(CONVERTER) $(CONV_SRC)

# Compile tester:
$(TESTER): modules $(TEST_SRC)
	$(FC) -o $(TESTER) $(FC_INC)../../modules ../../$(MODULEPATH)/*.o $(TESTER_SRC)

# Build input file for ascii-to-binary converter:
$(EPH_INPUT):
	$(FTP_CLIENT) ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/*
	cat $(EPH_ASCII) > $(EPH_INPUT)
	rm -f $(EPH_ASCII)

# Write back-up:
backup:
	$(SHELL) -c "if test -d ../backup_$(PROJNAME); then true; else mkdir ../backup_$(PROJNAME); fi"
	cp -a * ../backup_$(PROJNAME)

# Make compressed tar-ball:
tar: clean
	rm -f $(PROJNAME)-*.tgz
	tar cvzf $(PROJNAME)-$(VERSION).tgz *

modules:
	cd ../../$(MODULEPATH) ; $(MAKE) all

# Clean up directory:
clean:
	rm -f *.o *~ core $(CONVERTER) $(TESTER) *.405
