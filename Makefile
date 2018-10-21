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
# This is a makefile for the classes, modules, executables,
# and documentation in the OpenOrb-project.
#
# Author:  MG
# Version: 2009-11-09

include make.config
include Makefile.include

PREFIX ?= /opt/oorb

# Write back-up:
backup:
	$(SHELL) -c "if test -d ../backup_$(PROJNAME); then true; else mkdir ../backup_$(PROJNAME); fi"
	cp -a * ../backup_$(PROJNAME)

# Build binary and Python module
all:
	cd $(MAINPATH)   && $(MAKE) oorb
	cd $(PYTHONPATH) && $(MAKE) pyoorb

# Make tar-ball:
tar: all_clean
	cd .. ; tar cvf $(PROJNAME)_v$(VERSION).tar $(PROJNAME) ; gzip $(PROJNAME)_v$(VERSION).tar

all_clean: clean
	cd $(DOCPATH)    ; $(MAKE) clean
	cd $(CLASSPATH)  ; $(MAKE) clean
	cd $(MODULEPATH) ; $(MAKE) clean
	cd $(MAINPATH)   ; $(MAKE) clean
	cd $(PYTHONPATH) ; $(MAKE) clean
	cd $(LIBPATH)    ; $(MAKE) clean

install:
	@echo "Installing into $(PREFIX)"
	mkdir -p $(PREFIX)/bin $(PREFIX)/etc $(PREFIX)/lib $(PREFIX)/data $(PREFIX)/python
	cp -a main/oorb $(PREFIX)/bin/
	cp -a main/oorb.conf $(PREFIX)/etc/
	cp -a lib/liboorb* $(PREFIX)/lib/
	cp -a data/* $(PREFIX)/data/ && rm -rf "$(PREFIX)/data/JPL_ephemeris"
	cp -a python/pyoorb*.so $(PREFIX)/python/

# Remove library and modules:
clean:
	rm -f *~ *.mod *.o

