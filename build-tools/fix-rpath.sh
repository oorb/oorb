#!/bin/bash
#
# If conda's gfortran compiler has been used, convert the relative linker
# paths that it inserts into absolute paths.  This makes the as-built binary
# non-relocatable, but it makes it possible to run it from the build
# directory. This is typically what a developer expects, as well as the
# odd user who will build oorb from source using conda's gfortran.
#
# When creating conda packages, conda-build will re-establish
# realocatability when it scans the binary for libs residing in
# $CONDA_PREFIX.
#

set -e

#
# This is relevant only on macOS
#
[[ $(uname) == Darwin ]] || exit 0

#
# See if there is the tell-tale signs of conda's compiler
#
if otool -L "$1" | grep -q '@rpath/libgfortran'; then
	#
	# compute the location of lib directory
	#
	F=$(which gfortran)
	D=$(dirname "$F")
	D="$D/../lib"
	
	for L in $(otool -L "$1" | grep '@rpath' | awk '{print $1}'); do
		## verify the library is there
		ABSPATH="$D/${L#@rpath/}"
		if test -f "$ABSPATH"; then
			install_name_tool -change "$L" "$ABSPATH" "$1"
		fi
	done
fi
