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
documentation below, the included ``test.py`` file, as well as `Lynne
Jones' Jupyter notebook
<https://github.com/rhiannonlynne/notebooks/blob/master/PyOorb%20Demo.ipynb>`_
for some example usecases.

Convenience functions for the use of this module are provided in the
framework of the `sbpy <http://sbpy.org>`_ project. Please refer to
the `sbpy.data` module for `more information
<https://sbpy.readthedocs.io/en/latest/sbpy/data.html>`_.

Installation
------------

``pyoorb`` is automatically installed alongside ``oorb`` if the `conda
installer <https://github.com/conda-forge/openorb-feedstock>`_ has
been used. This is the **recommended installation method** and by far
the easiest method, as well.

Manual installation procedures are provided here in case you do not
want to use the conda installer:

0. This is not a requirement but highly recommended: get the latest
   version of `Anaconda <https://www.anaconda.com/download>`_ Python
   3.6 and make it your default Python before you install OpenOrb. Ensure
   you have ``numpy`` and ``pytest`` installed by running:

       >>> conda install numpy pytest

1. Follow the OpenOrb installation guidelines including the generation
   and updating of the data files (no need to update ``ET minus UT``
   and ``TAI minus UTC``)

2. From the root directory, run:

       >>> make pyoorb

   The shared library will be built in ``oorb/python``.
       
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

   For a more comprehensive list of tests, run:

       >>> make test

   from the root ``oorb`` directory.

Documentation/API
-----------------

Defining Orbits 
^^^^^^^^^^^^^^^^

pyoorb uses the same input orbit format in all functions. Orbits
have to be in the form of a numpy.array (one-dimensional for one
orbit, two-dimensional for multiple orbits). Every orbit array has to
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
11. `absolute magnitude` of the target 
12. `photometric slope parameter` of the target 

In order to provide compatibility with the underlying FORTRAN library,
only use double values (``dtype=np.double``) and set the corresponding
order for the array (``order='F'``). Angles have to be provided in
radians, not in degrees; epochs have to be provided as modified Julian
dates.

The following example defines a single orbit:

    >>> import numpy as np
    >>> orbits = np.array([[0, 1.46905, 0.33435, np.deg2rad(14.3024),
    ...                     np.deg2rad(224.513), np.deg2rad(27.5419),
    ...                     np.deg2rad(324.697), 3, 51544.5, 1, 12.5, 0.15]],
    ...                   dtype=np.double, order='F')

Note that ``orbits`` is an array of arrays, where the latter ones
define the orbits. If you want to define more than one orbit, you can
simply add more orbit arrays. A more convenient way is shown in the
next example that uses an array transposition and defines a total of
three orbits:

    >>> import numpy as np
    >>> orbits = np.array(
    ...    np.array([[0, 1, 2],  # id
    ...              [1.46905, 2.49241, 2.43064],  # a
    ...              [0.33435, 0.78954, 0.70177],  # e
    ...              np.deg2rad([14.3024, 7.99749, 0.52734]),  # incl
    ...              np.deg2rad([224.513, 124.587, 334.990]),  # longnode
    ...              np.deg2rad([27.5419, 267.104, 205.376]),  # argper
    ...              np.deg2rad([324.697, 297.655, 5.56391]),  # mean anom
    ...              [3, 3, 3],  # orbit type
    ...              [51544.5, 51544.5, 51544.5],  # epoch
    ...              [1, 1, 1],  # timescale type
    ...              [12.5, 6.3, 20.5],  # absolute magnitude
    ...              [0.15, 0.15, 0.15]]).transpose(),  # slope parameter
    ...    dtype=np.double, order='F')



Initializing pyoorb
^^^^^^^^^^^^^^^^^^^

Before any pyoorb functionality can be used, the module has to be
initialized using the following two lines:

    >>> import pyoorb as oo
    >>> oo.pyoorb.oorb_init()

In case you installed ``pyoorb`` manually (i.e., you did not use the
conda installer), you have to manually define which ephemerides to
use:

    >>> import os
    >>> ephfile = os.path.join(os.getenv('OORB_DATA'), 'de430.dat')
    >>> oo.pyoorb.oorb_init(ephfile)

This initialization requires the ``'OORB_DATA'`` environment variable
to be properly defined (see installation guide above). Note that in
this example the ``DE430`` planetary and lunar ephemerides are used;
other definition files can be used, but those have to be present in
the ``'OORB_DATA'`` directory.


Orbital Element Transformation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Function ``pyoorb.oorb_element_transformation`` provides
transformations between different orbital element
definitions. Required parameters are ``in_orbits`` (an orbit array as
defined above) and ``in_element_type`` (the orbital element type
integer code: ``CART``: 1, ``COM``: 2, ``KEP``: 3). The function
outputs the orbit array as defined above using the element scheme
defined through ``in_element_type``, as well as the corresponding
error code.

The following example transforms the orbit array defined above from a
Keplerian to a cartesian definition:

    >>> new_orbits, err = oo.pyoorb.oorb_element_transformation(
    ...     in_orbits=orbits,
    ...     in_element_type=1)
    >>> print(err)
    0
    >>> print(new_orbits)
    [[ 0.00000000e+00 -1.13248995e+00 -1.21090780e-01 -1.80398368e-01
       5.81046365e-03 -1.61918824e-02  3.98214412e-03  1.00000000e+00
       5.15445000e+04  1.00000000e+00  1.25000000e+01  1.50000000e-01]
                                     ...
     [ 2.00000000e+00 -6.02921087e-01 -5.37976414e-01 -6.83341536e-03
       1.03814468e-02 -2.24243823e-02 -1.46645896e-04  1.00000000e+00
       5.15445000e+04  1.00000000e+00  2.05000000e+01  1.50000000e-01]]    

The definitions of the individual columns are provided above.

Orbit Propagation
^^^^^^^^^^^^^^^^^

The function ``pyoorb.oorb_propagation`` propagates one or more orbits
(``in_orbits``) to a desired epoch (``in_epoch``). The epoch has to be
provided as an array with ``dtype=np.double, order='F'`` (see example
below) and has to contain two elements: the epoch as modified Julian
date and the corresponding timescale type (``UTC``: 1, ``UT1``: 2,
``TT``: 3, ``TAI``: 4). The user can decide between an N-body
integration (``in_dynmodel='N'``) or a faster but less accurate
two-body integration (``in_dynmodel='2'``).

The following example creates a target epoch and propagates all three
orbits defined above using an N-body integration:

    >>> epoch = np.array([51232.23, 3], dtype=np.double, order='F')
    >>> orb, err = oo.pyoorb.oorb_propagation(in_orbits=orbits,
    ...                                       in_epoch=epoch,
    ...                                       in_dynmodel='N')
    >>> print(err)
    0
    >>> print(orb)
    [[0.00000000e+00 1.46902259e+00 3.34239111e-01 2.49612224e-01
      3.91849109e+00 4.81057800e-01 2.64937688e+00 3.00000000e+00
      5.12322300e+04 3.00000000e+00 0.00000000e+00 0.00000000e+00]
     [1.00000000e+00 2.49303082e+00 7.89455978e-01 1.39563343e-01
      2.17499921e+00 4.66102654e+00 3.83095531e+00 3.00000000e+00
      5.12322300e+04 3.00000000e+00 0.00000000e+00 0.00000000e+00]
     [2.00000000e+00 2.43234703e+00 7.01646868e-01 9.21820573e-03
      5.85002081e+00 3.58026888e+00 4.96413849e+00 3.00000000e+00
      5.12322300e+04 3.00000000e+00 0.00000000e+00 0.00000000e+00]]

Ephemeris Computation
^^^^^^^^^^^^^^^^^^^^^

The function ``pyoorb.oorb_ephemeris_full`` computes ephemeris for
orbits ``in_orbits`` relative to observer location ``in_obscode`` (the
official `Minor Planet Center observatory code
<https://minorplanetcenter.net/iau/lists/ObsCodesF.html>`_) and for
epochs ``in_date_ephems``. Epochs are defined as arrays with
``dtype=np.double, order='F'`` (see example below) containing
two-element arrays with the epoch as modified Julian date and the
corresponding timescale type (``UTC``: 1, ``UT1``: 2, ``TT``: 3,
``TAI``: 4). The user can decide between an N-body
integration (``in_dynmodel='N'``) or a faster but less accurate
two-body integration (``in_dynmodel='2'``).

The following example computes ephemeris for the orbits defined above,
as seen from Maunakea, and for a range of epochs, using an N-body
integration:

    >>> mjds = np.arange(51232, 51233, 1/24)
    >>> epochs = np.array(list(zip(mjds, [1]*len(mjds))), dtype=np.double, order='F')
    >>> eph, err = oo.pyoorb.oorb_ephemeris_full(in_orbits=orbits,
    ...                                          in_obscode='568',
    ...                                          in_date_ephems=epochs,
    ...                                          in_dynmodel='N')
    >>> print(err)
    0
    >>> print(eph)
    [[[ 5.12320000e+04  2.97420305e+01  9.39898382e+00 ... -8.88643901e-01
        4.34680947e-01  1.86372479e-05]
      [ 5.12320417e+04  2.97612057e+01  9.40381070e+00 ... -8.88970865e-01
        4.34043644e-01  1.44636163e-05]
      [ 5.12320833e+04  2.97803611e+01  9.40862679e+00 ... -8.89300047e-01
        4.33405588e-01  1.05174851e-05]
      ...
      [ 5.12328750e+04  1.71222795e+01  7.66784905e+00 ... -8.95406065e-01
        4.21023347e-01  2.93152579e-05]
      [ 5.12329167e+04  1.71324739e+01  7.67202457e+00 ... -8.95715537e-01
        4.20379156e-01  2.60843567e-05]
      [ 5.12329583e+04  1.71426400e+01  7.67619777e+00 ... -8.96026864e-01
        4.19736086e-01  2.22713554e-05]]]

``eph`` is a nested array with one element per input orbit, one
element per epoch, and 33 properties that are calculated by pyoorb. In
the case of ``pyoorb.oorb_ephemeris_full``, these properties are:

0. modified julian date
1. right ascension (deg)
2. declination (deg)
3. dra/dt sky-motion (deg/day, including cos(dec) factor)
4. ddec/dt sky-motion (deg/day)
5. solar phase angle (deg)
6. solar elongation angle (deg)
7. heliocentric distance (au)
8. geocentric distance (au)
9. predicted apparent V-band magnitude
10. position angle for direction of motion (deg)
11. topocentric ecliptic longitude (deg)
12. topocentric ecliptic latitude (deg)
13. opposition-centered topocentric ecliptic longitude (deg)
14. opposition-centered topocentric ecliptic latitude (deg)
15. heliocentric ecliptic longitude (deg)
16. heliocentric ecliptic latitude (deg)
17. opposition-centered heliocentric ecliptic longitude (deg)
18. opposition-centered heliocentric ecliptic latitude (deg)
19. topocentric object altitude (deg)
20. topocentric solar altitude (deg)
21. topocentric lunar altitude (deg)
22. lunar phase [0...1]
23. lunar elongation (deg, distance between the target and the Moon)
24. heliocentric ecliptic cartesian x coordinate for the object (au)
25. heliocentric ecliptic cartesian y coordinate for the object (au)
26. heliocentric ecliptic cartesian z coordinate for the objects (au)
27. heliocentric ecliptic cartesian x rate for the object (au/day)
28. heliocentric ecliptic cartesian y rate for the object (au/day)
29. heliocentric ecliptic cartesian z rate for the objects (au/day)
30. heliocentric ecliptic cartesian x coordinate for the observatory (au)
31. heliocentric ecliptic cartesian y coordinate for the observatory (au)
32. heliocentric ecliptic cartesian z coordinate for the observatory (au)
33. true anomaly (deg)

``pyoorb.oorb_ephemeris_basic`` only provides a subset of these
properties, enabling fast computations and requiring less memory:

0. modified julian date
1. right ascension (deg)
2. declination (deg)
3. dra/dt sky-motion (deg/day, including cos(dec) factor)
4. ddec/dt sky-motion (deg/day)
5. solar phase angle (deg)
6. solar elongation angle (deg)
7. heliocentric distance (au)
8. geocentric distance (au)
9. predicted apparent V-band magnitude
10. true anomaly (deg)


     
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

Original wrapper developer: F. Pierfederici <fpierfed@gmail.com>

This code has been modified by Michael Mommert to be of use to a
broader community in the framework of the `sbpy project
<http://sbpy.org>`_.
