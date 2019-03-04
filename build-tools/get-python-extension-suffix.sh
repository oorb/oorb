#!/bin/bash

PYTHON=${PYTHON:-python}

# figure out the major python version
PYVER=$($PYTHON -c 'import sys; print(sys.version_info[0])')

if [[ $PYVER == 2 ]]; then
	echo ".so"
elif [[ $PYVER == 3 ]]; then
	python3-config --extension-suffix
fi
