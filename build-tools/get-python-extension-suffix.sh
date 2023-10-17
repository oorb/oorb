#!/bin/bash

# Prefer '$PYTHON', or 'python3', or 'python', in that order
PYTHON=${PYTHON:-python3}
if ! command -v $PYTHON >/dev/null 2>&1; then
	PYTHON=python
fi

# figure out the major python version
PYVER=$($PYTHON -c 'import sys; print(sys.version_info[0])')

if [[ $PYVER == 2 ]]; then
	echo ".so"
elif [[ $PYVER == 3 ]]; then
	$PYTHON -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))'
fi
