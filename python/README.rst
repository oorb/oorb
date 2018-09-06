pyoorb - A Python Wrapper for OpenOrb
=====================================

This submodule builds a Python wrapper for OpenOrb. Most of the code
has been developed by the LSST Data Management team (see below). It
has been slightly modified to be of use to a broader community.

Current State
-------------

In its current state, the wrapper provides functionality for
transformations between orbital element definitions, the computation
of ephemerides, and orbit propagation - ranging and other orbit
fitting capabilities will be added in the near future.

Some documentation and installation guides can be found below - more
documentation will be added in the near future. Please refer to the
``test.py`` file for some example usecases, as well as `Lynne Jones'
Jupyter notebook
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
   ``PYTHONPATH``; additional variables are listed below that have to
   be set as well. Using bash, the corresponding lines would look like
   this:
   
       >>> export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:<oorb>/python/"
       >>> export PYTHONPATH="$PYTHONPATH:<oorb>/python/"
       >>> export OORB_DATA="<oorb>/data"
       >>> export OORB_CONF="<oorb>/main/oorb.conf"
   
   where you replace ``<oorb>`` with the corresponding path to the
   oorb directory on your system.

4. Now, pyoorb is installed on your system and you can test it by running:

       >>> python test.py

   in your ``oorb/python`` directory. This will test different
   functions and should generate a lot of numerical output. Each
   function should complete with error code 0.


Documentation/API
-----------------

Defining Orbits 
^^^^^^^^^^^^^^^^

`pyoorb` uses the same input orbit format in all functions. Orbits
have to be in the form of a `numpy.array` (one-dimensional for one
orbit, two-dimensional for multiple orbits). very orbit array has to
contain the following 12 properties (in this order), some of which
depend on the type of the elements provided (``COM`` for cometary
orbits, ``KEP`` for Keplerian orbits, and ``CART`` for cartesian
coordinates):

1. `orbit id`: an integer number to identify the orbit; usually ranges
   from 0 to `n`-1, where `n` is the number of orbits
2. `perihelion distance` (au) for ``COM``, semimajor axis `a` (au) for
   ``KEP``, `x` (au) for ``CART``
3. `eccentricity` for ``COM`` or ``KEP``, `y` (au) for ``CART``
4. `inclination` (deg) for ``COM`` or ``KEP``, `z` (au) for ``CART``
5. `longitude of the ascending node` (deg) for ``COM`` and ``KEP``,
   `dx/dt` (au/day) for ``CART``
6. `argument of the periapsis` (deg) for ``COM`` and ``KEP``, `dy/dt`
   (au/day) for ``CARqT``
7. `epoch of perihelion` (modified Julian date) for ``COM``, mean
   anomaly (deg) for ``KEP``, `dz/dt` for ``CART``
8. `orbital element type`; integer value: ``CART``: 1, ``COM``: 2,
   ``KEP``: 3
9. `epoch` of the osculating elements (modified Julian date)
10. `timescale type` of the epochs provided; integer value: ``UTC``:
    1, ``UT1``: 2, ``TT``: 3, ``TAI``: 4
11. `absolute magnitude` (``KEP`` or ``CART``) or `M1` parameter
    (``COM``)
12. `photometric slope parameter` (``KEP`` or ``CART``) or `K1`
    parameter (``COM``)

In order to provide compatibility with the underlying FORTRAN library,
only use double values and set the corresponding order for the array
(``order='F'``). The following lines provide an example orbit array:

    >>> import numpy as np
    >>> orbits = np.array([[0,1,2],  # id
    ...                    [2.45123, 1.2343, 7.1235],  # a 
    ...                    [0.1234, 0.05453, 0.0000234],  # e
    ...                    [3.6, 14.5322, 0.002323],  # i
    ...                    [231.34534543, 23.345345, 45.42342],  # longnode
    ...                    [345.4324324, 125.243324, 45.34242],  # argper
    ...                    [13.234234, 56.234234, 184.234324],  # mean anom
    ...                    [3, 3, 3],  # orbit type
    ...                    [53124.0, 51624.0, 52623.0],  # epoch
    ...                    [1, 1, 1],  # timescale type
    ...                    [12.5, 6.3, 20.5],  # absolute magnitude
    ...                    [0.15, 0.2, 0.15]],  # slope parameter
    ...                   dtype=np.double, order='F')
    

Orbital Element Transformation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Orbit Propagation
^^^^^^^^^^^^^^^^^


Ephemerides Computation
^^^^^^^^^^^^^^^^^^^^^^^
   

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

This code has been modified by Michael Mommert
(<mommermiscience@gmail.com>) to be of use to a broader community in the
framework of the `sbpy project <http://sbpy.org>`_.
