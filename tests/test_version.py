#
# Check that the version is correctly compiled into
# the binaries.
#

from __future__ import print_function
from subprocess import check_output
import os
from pkg_resources import parse_version

import pytest

def shell(cmd):
	return check_output(cmd, shell=True).decode('utf-8')

@pytest.mark.skipif(os.environ.get("PYOORB", None) == "0", reason="pyoorb not built")
def test_version_pyoorb():
	expectver = check_output("./build-tools/compute-version.sh", shell=True).decode('utf-8').rstrip()

	import pyoorb
	assert pyoorb.__version__ == expectver

def test_version_oorb():
	expectver = 'v' + check_output("./build-tools/compute-version.sh", shell=True).decode('utf-8').rstrip()
	ver = check_output('./main/oorb  --version | grep "^OpenOrb v" | cut -d\\  -f 2', shell=True).decode('utf-8').rstrip()

	assert ver == expectver
