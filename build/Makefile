############################################################################
#
# OpenOrb Makefile
#
# See the main `Makefile` for the documentation.
#
# Author: mjuric@astro.washington.edu (http://github.com/mjuric)
#
#############################################################################

include ../make.config
include ../Makefile.include

#
# Default rule builds all command-line executables and the pyoorb Python
# module (if configured with --with-pyoorb).
#
# All build targets
ALLTARGETS=$(PROGRAMS) liboorb
ifeq ("$(PYOORB)","1")
ALLTARGETS+=pyoorb
endif

.PHONY: all
all: $(ALLTARGETS)

#
# Clean all temp files
#
.PHONY: clean
clean:
	rm -f *.o *.mod *.tmp
	rm -f version.h
	rm -f $(addprefix ../main/, $(PROGRAMS)) 
	rm -f ../lib/liboorb.{$(LIBEXT),a}
	rm -rf *.pyf ../python/pyoorb.{so,dylib} ../python/pyoorb.*.so ../python/pyoorb.{name,sp_dir} _pyoorb_build
	rm -rf ../bin ../lib

#
# Convenience so we can run 'make oorb' (instead of 'make ../bin/oorb')
#
.PHONY: $(PROGRAMS)
$(PROGRAMS): %: ../bin/%

#
# Build programs
#
# These rules leave binaries both in bin/ and main/, the latter for
# backwards compatibility.
#
# See https://stackoverflow.com/a/34983297 on why static pattern
# rules are needed here (otherwise, deleting main/oorb and running
# make oorb would not rebuild it, as it's considered an intermediate
# file by make).

$(addprefix ../bin/, $(PROGRAMS)): ../bin/%: ../main/%
	mkdir -p ../bin
	cp -a $< $@

$(addprefix ../main/, $(PROGRAMS)): ../main/%: %.o ../lib/liboorb.a
	$(FC) $(FCOPTIONS) $^ -o $@.tmp
	fix-rpath.sh $@.tmp
	mv $@.tmp $@

#
# Build version include file. This is included by oorb.f90 to
# compile-in the current version
#
version.h: ../VERSION
	generate-version-h.sh $$(cat $<) > $@

#
# Build liboorb, static and dynamic (in lib/ directory)
#
.PHONY: liboorb
liboorb: ../lib/liboorb.$(LIBEXT) ../lib/liboorb.a

LIBOBJECTS=$(notdir ${LIBSOURCES:.f90=.o})
../lib/liboorb.$(LIBEXT): $(LIBOBJECTS)
	mkdir -p ../lib
	$(FC) $(FCOPTIONS) $(FC_SHARED) $(LIBOBJECTS) -o $@

../lib/liboorb.a: $(LIBOBJECTS)
	mkdir -p ../lib
ifeq ($(shell uname),Darwin)
	libtool -static -no_warning_for_no_symbols $(LIBOBJECTS) -o $@
else
	ar rcs $@ $(LIBOBJECTS)
endif

ifeq ("$(PYOORB)","1")
#
# Build pyoorb
#
# These rules drop the compiled module both in python/ and lib/, the latter
# for backwards compatibility.
#

# Full name of the pyoorb library (e.g., pyoorb.cpython-36m-darwin.so)
PYOORB_SO            = pyoorb$(shell PYTHON=$(PYTHON) $(ROOT_DIR)/build-tools/get-python-extension-suffix.sh)

.PHONY: pyoorb
pyoorb: ../lib/$(PYOORB_SO)

../lib/$(PYOORB_SO): ../python/$(PYOORB_SO)
	mkdir -p ../lib && cp -a $< $@

../python/$(PYOORB_SO): ../lib/liboorb.a pyoorb.pyf pyoorb.o pyoorb.mod
	@ # Check for prerequisites and output a helpful message if they're not installed
	@ hash $(F2PY) 2>/dev/null && hash $(PYTHON) 2>/dev/null || \
		{ echo "****************************************************************************"; \
		  echo "  You need python and f2py (part of numpy) to build pyoorb."; \
		  echo "  If you don't care about building pyoorb, run 'make oorb' instead."; \
		  echo "  Otherwise, see python/README.rst for how to set up the build environment."; \
		  echo "****************************************************************************"; \
		  exit -1; \
		}
	@ [[ ! -z "$(LDFLAGS)" ]] && echo "$$LDFLAGS" | grep -E -v -q -- "(-bundle|-shared)" && { \
	  [[ `uname` == "Darwin" ]] && LDFLAGS="$$LDFLAGS -bundle" || LDFLAGS="$$LDFLAGS -shared"; \
	}; \
	: \ "# We override whatever is in F77= as that confuses f2py when used with Intel Compilers (conda compiler packages set it)"
	F77= $(F2PY) --quiet -m pyoorb pyoorb.o pyoorb.pyf ../lib/liboorb.a --build-dir ./_pyoorb_build -c --noarch \
	        $(F2PY_FCOMPILER) --f90exec=$(FC) --f90flags="$(FC_INC)../build" && \
	mv $(PYOORB_SO) ../python
	@ echo
	@ # FIXME: this exists so `sudo make install` knows what to install. It's hacky this way
	@ #        and should be made to be generated/updated like the VERSION file is.
	echo $(PYOORB_SO) > ../python/pyoorb.name
	$(PYTHON) -c "import site; print(site.getsitepackages()[0])" > ../python/pyoorb.sp_dir
	@ echo "Built the $(PYOORB_SO) module."

pyoorb.pyf: ../VERSION ../python/pyoorb.f90
	$(F2PY) --quiet -m pyoorb ../python/pyoorb.f90 --overwrite-signature -h pyoorb.pyf.tmp
	@ echo
	@ # generate __version__ string code
	generate-version-pyf.sh $$(cat ../VERSION) > version.pyf
	@ # splice the generated version information so pyoorb.__version__ works
	sed '/    interface  ! in :pyoorb/ r version.pyf' pyoorb.pyf.tmp > pyoorb.pyf.tmp2
	@ # make sure it worked
	@ grep -q '__version__' pyoorb.pyf.tmp2 || { echo "Failed to splice pyoorb.__version__ information into pyoorb.pyf. Verify the sed stanza in the Makefile is not out of date."; exit -1; }
	mv pyoorb.pyf.tmp2 pyoorb.pyf
	rm -f pyoorb.pyf.tmp version.pyf
else
.PHONY: pyoorb
pyoorb:
	@echo "Not configured to build pyoorb. Rerun ./configure with --with-pyoorb switch and try again."
endif

#
# Rules to generate dependency information using gfortran's preprocessor.
#
# These are not intented to be run by the users, just the developers
# whenever the dependency information changes (i.e., if a source file is
# added, or module USE information is changed).
#
# Note: they're only known to work with gfortran; the generated file
# will work with other compilers, assuming they follow gfortran
# conventions on naming of .mod files.
#
.PHONY: depends
depends: clean version.h
	make-fortran-deps.sh $(ALLSOURCES) > make.depends

#
# Recipes to build individual files, which are called by the rules generated
# in make.depends
#
# We can't use pattern rules here because something like '%.mod: %.f90' will
# fail when the module name differs from .f90 file name.
#
# The complex form of MODBUILDCMD is to correctly handle parallel builds
# when more than one .mod file is generated from the same .f90 file, e.g.:
#
#    foo.mod bar.mod: baz.f90
#            $(MODBUILDCMD)
#
# make will execute $(MODBUILDCMD) twice and in parallel (once for foo.mod,
# once for bar.mod), but each of the invocations will try to write _both_
# files (and thus we get issues with overwriting).  The stanza below avoids
# this by running each build in a temporary directory, and then just copying
# back the target .mod file.
#
OBJBUILDCMD=$(FC) -c $(FCOPTIONS) $(FC_INC). $<
MODBUILDCMD=@test -f $@ || ( TMPDIR=$$(mktemp -d) && cd "$$TMPDIR" && $(FC) -c $(FCOPTIONS) $(FC_INC)$(PWD) $(abspath $<) && mv $@ $(PWD) && rm -r "$$TMPDIR"  )
-include make.depends
