pyoorb - A Python Wrapper for OpenOrb
=====================================

This submodule builds a Python wrapper for OpenOrb. Most of the code
has been developed by the LSST Data Management team (see below). It
has been slightly modified to be of use to a broader community.

Current State
-------------

In its current state, the wrapper only provides orbit propagation
functionality - ranging and other orbit fitting capabilities will be
added in the near future.

Documentation will be added in the near future, too. For now, we refer
to the ``test.py`` file for an example usecase, as well as `Lynne
Jones' Jupyter notebook
<https://github.com/rhiannonlynne/notebooks/blob/master/PyOorb%20Demo.ipynb>`_.


Installation
------------

0. This is not a requirement but highly recommended: get the latest
   version of `Anaconda <https://www.anaconda.com/download>`_ Python
   3.6 and make it your default Python before you install OpenOrb.

1. Follow the OpenOrb installation guidelines including the generation
   and updating of the data files (no need to update ``ET minus UT``
   and ``TAI minus UTC``)

2. From the ``oorb/python``
   directory, run the following commands to build pyoorb:

       >>> make clean
       >>> make
       >>> make pyoorb
       
3. A few environment variables have to be set before pyoorb can be
   used from any location on your machine. In order to do so, open
   your ``.bashrc``, ``.cshrc``, or ``.profile`` file in your home
   directory with a text editor and add the ``oorb/python`` directory
   on your machine to the variables ``LD_LIBRARY_PATH`` and
   ``PYTHONPATH``. Using bash, the corresponding lines would look like
   this:
   
       >>> export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:<oorb/python>/"
       >>> export PYTHONPATH="$PYTHONPATH:<oorb/python>/"

   where you replace ``<oorb/python>`` with the path on your system.

4. Now, pyoorb is installed on your system and you can test it by running:

       >>> python test.py

   in your ``oorb/python`` directory. This will test different
   functions and should generate a lot of numerical output. Each
   function should complete with error code 0.


Acknowledgements and License Information
----------------------------------------

LSST Data Management System
Copyright 2008, 2009 LSST Corporation.

This product includes software developed by the
LSST Project (http://www.lsst.org/).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the LSST License Statement and
the GNU General Public License along with this program.  If not,
see <http://www.lsstcorp.org/LegalNotices/>.

original wrapper developer: F. Pierfederici <fpierfed@gmail.com>

This code has been slightly modified to be of use to a broader
community in the framework of the `sbpy project <http://sbpy.org>`_.
