#!/usr/bin/env python
# F. Pierfederici (fpierfed@lsst.org)
from distutils.core import setup



if __name__ == "__main__": 
    setup(name = 'oorb', 
          description = "OpenOrb module", 
          author = "Francesco Pierfederici", 
          author_email = "fpierfed@lsst.org",
          version='1.0a1',
          package_dir={'oorb': '.'},
          packages=['oorb', ],
) 
