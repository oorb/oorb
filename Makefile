############################################################################
#
# OpenOrb Build System
#
# This file and its companion in build/Makefile define build targets for all
# command line executables (see PROGRAMS in ../make.config), liboorb
# libraries, and the pyoorb Python module.  Default rule (named `all`)
# builds all of them.
#
# Example usage (with -j4 to take advantage of multi-threaded builds):
#
#	$ make -j4
#
# or (for example) to just build oorb, run:
#
#	$ make -j4 oorb
#
# To install the code, run:
#
#	$ make install
#
# This will install OpenOrb to whatever you specified as --prefix when
# running configure, and pyoorb to the standard modules location (i.e., the
# "site-path") of python used to build pyoorb.  The location where pyoorb is
# installed can be given explicitly via the PYTHON_SITE_PATH environmental
# variable.
#
# FOR DEVELOPERS:
#
# * The build is executed in the build/ subdirectory.  All intermediate
#   files end up there (.o, .mod, etc.)
#
# * To add new programs, add the .f90 file into main/ and list it in
#   PROGRAMS variable in make.config. Then rerun `make depends`.
#
# * To add new files, add them to modules/ or classes/ and list the file
#   name in the apropriate variable in make.config.  Then rerun `make
#   depends`
#
# * If there's a change in depdendency of any FORTRAN file (e.g., you've
#   USEd another module), run:
#
#       $ make depends
#
#   to rebuild the dependencies file (build/make.depends).  You must be
#   using gfortran for this rebuild to work.
#
# Author: mjuric@astro.washington.edu (http://github.com/mjuric)
#
#############################################################################


include make.config
include Makefile.include

# Verify $PREFIX has been set -- this may happen if the developer
# just `git pull`-ed the version of the code that introduced $(PREFIX).
ifeq ("$(PREFIX)", "")
    $(error "Install prefix is not set; please rerun ./configure and try again.")
endif

.PHONY: all
all:
	@ $(MAKE) -C build $@

ephem:
	@ $(MAKE) -C data/JPL_ephemeris

# Forward everything we don't recognize to the makefile in build/
%:
	$(MAKE) -C build $@

.PHONY: clean
clean:
	$(MAKE) -C build clean
	$(MAKE) -C doc clean

# Make tar-ball:
.PHONY: tar
tar: clean
	cd .. && tar czvf $(PROJNAME)_v$(VERSION).tar.gz --exclude $(PROJNAME)/.git $(PROJNAME)

.PHONY: install
install:
	install -d $(PREFIX)/{bin,etc,lib,share/oorb}

	for P in $(PROGRAMS); do \
		install bin/$$P $(PREFIX)/bin; \
	done

	install -m644 lib/liboorb.a         $(PREFIX)/lib
	install       lib/liboorb.$(LIBEXT) $(PREFIX)/lib
	install -m644 main/oorb.conf        $(PREFIX)/etc
	cp -a $(shell find data -maxdepth 1 -type f) $(PREFIX)/share/oorb

ifeq ("$(PYOORB)","1")
	install       python/$(shell cat python/pyoorb.name) $(shell echo $${PYTHON_SITE_PATH:-$$(cat python/pyoorb.sp_dir)})
endif

.PHONY: test
test: all
	@hash $(PYTEST) 2>/dev/null || { echo "You need to have pytest installed to run the tests." && exit -1; }
	PYTHONPATH="lib:$$PYTHONPATH" PYOORB=$(PYOORB) $(PYTEST) tests
ifeq ("$(PYOORB)","1")
# integration tests, will run only if JPL ephemeris data has been downloaded
ifneq ("$(wildcard data/de430.dat)","")
	PYTHONPATH="lib:$$PYTHONPATH" OORB_DATA=data $(PYTHON) python/test.py
else
	@ echo WARNING: Not running pyoorb integration tests as data/de430.dat is not present. Run "'make ephem'" to build it first.
endif
endif
