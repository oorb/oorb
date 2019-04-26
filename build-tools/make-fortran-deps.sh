#!/bin/bash
#
#  Inspired by the solution by FX Coudert
#       https://gcc.gnu.org/ml/fortran/2007-11/msg00013.html
#  and further discussed at
#       http://lagrange.mechse.illinois.edu/f90_mod_deps/
#

# Verify we're using gfortran -- ifort is known to be buggy with dependency
# generation

($FC --version 2>/dev/null | grep -q "GNU Fortran") || {
	echo "*** Running 'make depends' requires the compiler to be gfortran." >&2
	exit -1
}

$FC -cpp -I. -M "$@" | while read line; do
	IFS=":" read -r TARGETS DEPS <<< "$line"

	# FIXME: I assume gfortran will output .o first
	IFS=" " read -r OBJ MODULES <<< "$TARGETS"

	# FIXME: I assume gfortran will output the corresponding .f90 first
	IFS=" " read -r SRC <<< "$DEPS"

	# Separate dependency for the object file
	echo "$OBJ: $DEPS"
	echo $'\t''$(OBJBUILDCMD)'

	if [[ ! -z $MODULES ]]; then
		# Separate dependency for the module file(s), if any
		echo "$MODULES: $SRC $OBJ"
		echo $'\t''$(MODBUILDCMD)'
	fi
done
