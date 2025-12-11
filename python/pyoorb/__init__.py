# Licensed under GNU General Public License Version 3 - see COPYING

from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass

# import pyoorb namespace from FORTRAN library
from . import pyoorb as _pyoorb
from .pyoorb import pyoorb
