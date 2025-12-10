cmake_minimum_required(VERSION 3.15...3.30)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)

find_package(LAPACK)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# Grab the variables from a local Python installation F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_custom_command(
  OUTPUT prefix.h
  COMMAND echo 'CHARACTER(len=0) :: PREFIX = \"\"' > ${CMAKE_CURRENT_SOURCE_DIR}/prefix.h)

# Generate the wrapper files
add_custom_command(
  OUTPUT pyoorbmodule.c pyoorb-f2pywrappers2.f90
  COMMAND Python::Interpreter -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/python/pyoorb.f90" -m pyoorb
  VERBATIM
  DEPENDS python/pyoorb.f90)

Python_add_library(
  pyoorb
  MODULE
  WITH_SOABI
  modules/linal.f90
  modules/cl_options.f90
  modules/data_structures.f90
  modules/integrators.f90
  modules/estimators.f90
  modules/random.f90
  modules/functions.f90
  modules/sort.f90
  modules/planetary_data.f90
  modules/utilities.f90
  modules/parameters.f90
  modules/statistics.f90
  classes/Orbit_class.f90
  classes/Observation_class.f90
  classes/Time_class.f90
  classes/CartesianCoordinates_class.f90
  classes/File_class.f90
  classes/SphericalCoordinates_class.f90
  classes/Observatory_class.f90
  classes/Unit_class.f90
  classes/Observations_class.f90
  classes/PhysicalParameters_class.f90
  classes/Base_class.f90
  classes/Observatories_class.f90
  classes/StochasticOrbit_class.f90
  prefix.h
  "${CMAKE_CURRENT_BINARY_DIR}/pyoorbmodule.c"
  "${CMAKE_CURRENT_BINARY_DIR}/pyoorb-f2pywrappers2.f90"
  "${F2PY_INCLUDE_DIR}/fortranobject.c"
  python/pyoorb.f90
)

target_include_directories(pyoorb PUBLIC ${F2PY_INCLUDE_DIR})
target_link_libraries(
  pyoorb
  PUBLIC
  Python::NumPy
  LAPACK::LAPACK
)

install(TARGETS pyoorb DESTINATION .)