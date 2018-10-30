!====================================================================!
!                                                                    !
! Copyright 2002-2011,2012                                           !
! Mikael Granvik, Jenni Virtanen, Karri Muinonen, Teemu Laakso,      !
! Dagmara Oszkiewicz, Lynne Jones                                    !
!                                                                    !
! This file is part of OpenOrb.                                      !
!                                                                    !
! OpenOrb is free software: you can redistribute it and/or modify it !
! under the terms of the GNU General Public License as published by  !
! the Free Software Foundation, either version 3 of the License, or  !
! (at your option) any later version.                                !
!                                                                    !
! OpenOrb is distributed in the hope that it will be useful, but     !
! WITHOUT ANY WARRANTY; without even the implied warranty of         !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  !
! General Public License for more details.                           !
!                                                                    !
! You should have received a copy of the GNU General Public License  !
! along with OpenOrb. If not, see <http://www.gnu.org/licenses/>.    !
!                                                                    !
!====================================================================!
!
!! Module for which f2py builds Python wrappers.
!!
!! @author  MG, FP, LJ
!! @version 2012-09-21
!!
!! updated 2018-09-05 by Michael Mommert (sbpy.org)

MODULE pyoorb

  USE planetary_data
  USE Base_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  USE Observatories_cl
  USE Orbit_cl
  USE StochasticOrbit_cl
  USE PhysicalParameters_cl
  IMPLICIT NONE
  TYPE (Observatories), SAVE, PRIVATE :: obsies
  CHARACTER(len=11), DIMENSION(6), PUBLIC :: element_types = (/ &
       "cartesian  ",&
       "cometary   ",&
       "keplerian  ",&
       "delaunay   ",&
       "poincare   ",&
       "equinoctial" &
       /)
  CHARACTER(len=3), DIMENSION(4), PUBLIC :: timescales = (/     &
       "UTC",&
       "UT1",&
       "TT ",&
       "TAI" &
       /)

CONTAINS

  SUBROUTINE oorb_init(ephemeris_fname, error_verbosity, &
       info_verbosity, error_code)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: ephemeris_fname
    INTEGER, INTENT(in), OPTIONAL :: error_verbosity
    INTEGER, INTENT(in), OPTIONAL :: info_verbosity
    INTEGER, INTENT(out) :: error_code

    TYPE (Time) :: t

    ! Set path to data files by reading env parameter OORB_DATA_DIR:
    CALL setAccessToDataFiles()

    ! Read OBSCODE.dat
    CALL NEW(obsies)
    IF (error) THEN
       error_code = 1
       error = .FALSE.
       RETURN       
    END IF

    ! Read TAI-UTC.dat and ET-UT.dat
    CALL NEW(t)
    IF (error) THEN
       error_code = 2
       error = .FALSE.
       RETURN       
    END IF

    ! Read de430.dat
    IF (PRESENT(ephemeris_fname) .AND. LEN_TRIM(ephemeris_fname) /= 0) THEN
       CALL JPL_ephemeris_init(error, ephemeris_fname)
       IF (error) THEN
          error_code = 3
          error = .FALSE.
          RETURN
       END IF
    ELSE
       CALL JPL_ephemeris_init(error, &
            filename=TRIM(OORB_DATA_DIR) // "/" // TRIM(EPH_FNAME)) 
       IF (error) THEN
          error_code = 4
          error = .FALSE.
          RETURN
       END IF
    END IF

    ! Set verbosity of error messages (0=nothing,5=everything,1=default)
    IF (PRESENT(error_verbosity)) THEN
       err_verb = error_verbosity
    ELSE
       err_verb = 1
    END IF

    ! Set verbosity of info messages (0=nothing,5=everything,1=default)
    IF (PRESENT(info_verbosity)) THEN
       info_verb = info_verbosity
    ELSE
       info_verb = 1
    END IF

  END SUBROUTINE oorb_init


  SUBROUTINE oorb_memfree()

    IMPLICIT NONE

    CALL nullifyTime()
    CALL NULLIFY(obsies)
    CALL JPL_ephemeris_nullify()

  END SUBROUTINE oorb_memfree


  SUBROUTINE oorb_element_transformation(in_norb, &
       in_orbits,       &
       in_element_type, &
       out_orbits,      &
       error_code)

    ! Input/Output variables.
    ! in_norb: number of input orbits
    INTEGER, INTENT(in)                                 :: in_norb
    ! in_orbits: input flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)           :: in_orbits ! (1:norb,1:12)
    ! in_element_type: type of output elements ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    INTEGER, INTENT(in)                                 :: in_element_type
    ! out_orbits: output flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(out)          :: out_orbits ! (1:norb,1:12)
    ! error_code: output error code
    INTEGER, INTENT(out)                                :: error_code

    TYPE (Time) :: t
    TYPE (Orbit) :: orb
    INTEGER :: i

    ! Loop over orbits:
    DO i=1,SIZE(in_orbits,dim=1)

       ! Get the element type from the input flattened orbit.
       IF(NINT(in_orbits(i,8)) < 0 .OR.                                       &
            NINT(in_orbits(i,8)) > SIZE(element_types)) THEN
          ! Error: unsupported orbital elements.
          error_code = 58
          RETURN
       END IF

       ! Get each flattened orbit and create an Orbit instance. First
       ! create a Time instance:
       CALL NEW(t, in_orbits(i,9), timescales(NINT(in_orbits(i,10))))
       IF (error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an Orbit instance:
       CALL NEW(orb, &
            in_orbits(i,2:7), &
            element_types(NINT(in_orbits(i,8))), &
            "ecliptic", &
            copy(t))
       CALL NULLIFY(t)

       out_orbits(i,1) = in_orbits(i,1)
       out_orbits(i,2:7) = getElements(orb, element_types(in_element_type), &
            "ecliptic")
       IF (error) THEN
          ! Error in transformation
          error_code = 58
          RETURN
       END IF
       out_orbits(i,8) = REAL(in_element_type,8)
       out_orbits(i,9:) = in_orbits(i,9:)
       CALL NULLIFY(orb)

    END DO

  END SUBROUTINE oorb_element_transformation


  SUBROUTINE oorb_propagation(in_norb, &
       in_orbits,       &
       in_epoch,        &
       in_dynmodel,     &
       out_orbits,      &
       error_code)

    ! Input/Output variables.
    ! in_norb: number of input orbits
    INTEGER, INTENT(in)                                 :: in_norb
    ! in_orbits: input flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)           :: in_orbits ! (1:norb,1:12)
    ! in_epoch: epoch and timescale of output orbits:
    ! (1) modified Julian date
    ! (2) timescale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    REAL(8), DIMENSION(2), INTENT(in)                   :: in_epoch
    ! in_model: "2"=2-body dynamical model, "N"=n-body dynamical model
    CHARACTER(len=1), INTENT(in)                        :: in_dynmodel
    ! out_orbits: output flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(out)          :: out_orbits ! (1:norb,1:12)
    ! error_code: output error code
    INTEGER, INTENT(out)                                :: error_code

    TYPE (Orbit) :: orb
    TYPE (Time) :: t0, t1
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    CHARACTER(len=6) :: dyn_model
    REAL(8) :: integration_step
    INTEGER :: i
    LOGICAL, DIMENSION(10) :: perturbers

    ! Init
    errstr = ""
    error_code = 0
    IF (in_dynmodel .EQ. "2") THEN
       dyn_model = "2-body"
    END IF
    IF (in_dynmodel .EQ. "N") THEN
       dyn_model = "n-body"
    END IF
    integrator = "bulirsch-stoer"
    integration_step = 5.0_8
    perturbers = .TRUE.

    ! Loop over orbits:
    DO i=1,SIZE(in_orbits,dim=1)

       ! Get the element type from the input flattened orbit.
       IF(NINT(in_orbits(i,8)) < 0 .OR. &
            NINT(in_orbits(i,8)) > SIZE(element_types)) THEN
          ! Error: unsupported orbital elements.
          error_code = 58
          RETURN
       END IF

       ! Get each flattened orbit and create an Orbit instance.
       ! First create a Time instance corresponding to the orbit epoch.
       CALL NEW(t0, in_orbits(i,9), timescales(NINT(in_orbits(i,10))))
       IF(error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an Orbit instance.
       CALL NEW(orb, &
            in_orbits(i,2:7), &
            element_types(NINT(in_orbits(i,8))), &
            "ecliptic", &
            copy(t0))
       CALL NULLIFY(t0)
       CALL setParameters(orb, &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
          error_code = 37
          RETURN
       END IF
       ! Create a Time instance based on the requested output epoch
       CALL NEW(t1, in_epoch(1), timescales(NINT(in_epoch(2))))
       IF (error) THEN
          ! Error in transformation
          error_code = 58
          RETURN
       END IF
       ! Compute topocentric ephemerides
       CALL propagate(orb, &
            t1)
       IF (error) THEN
          ! Error in getEphemerides()
          error_code = 40
          RETURN
       END IF
       out_orbits(i,1) = in_orbits(i,1)
       out_orbits(i,2:7) = getElements(orb, element_types(NINT(in_orbits(i,8))), &
            "ecliptic")
       IF (error) THEN
          ! Error in transformation
          error_code = 58
          RETURN
       END IF

       ! convert angles to degrees
       IF (element_types(NINT(in_orbits(i,8))) .EQ. "keplerian" .OR.  &
           element_types(NINT(in_orbits(i,8))) .EQ. "cometary") THEN
          out_orbits(i,4) = out_orbits(i,4)/rad_deg
          out_orbits(i,5) = out_orbits(i,5)/rad_deg
          out_orbits(i,6) = out_orbits(i,6)/rad_deg
          out_orbits(i,7) = out_orbits(i,7)/rad_deg
       END IF
          
       ! Set element type (same as input element type)
       out_orbits(i,8) = in_orbits(i,8)
       out_orbits(i,9) = in_epoch(1)
       out_orbits(i,10) = in_epoch(2)
       out_orbits(i,11) = in_orbits(i,11)
       out_orbits(i,12) = in_orbits(i,12)       
       CALL NULLIFY(orb)
       CALL NULLIFY(t1)
    END DO

  END SUBROUTINE oorb_propagation

  
  SUBROUTINE oorb_ephemeris_full(in_norb, &
       in_orbits,                    &
       in_obscode,                   &
       in_ndate,                     &
       in_date_ephems,               &
       in_dynmodel,                  &
       out_ephems,                   &
       error_code)
 
    ! in_norb: number of input orbits
    INTEGER, INTENT(in)                                 :: in_norb
    ! in_orbits: input flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)           :: in_orbits ! (1:norb,1:12)
    ! in_obscode: observatory code as defined by the Minor Planet Center
    CHARACTER(len=*), INTENT(in)                        :: in_obscode
    ! in_ndate: number of input epochs
    INTEGER, INTENT(in)                                 :: in_ndate
    ! in_date_ephems: epochs for ephemeris:
    ! (1) modified Julian date
    ! (2) timescale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    REAL(8), DIMENSION(in_ndate,2), INTENT(in)          :: in_date_ephems ! (1:ndate,1:2)
    ! in_model: "2"=2-body dynamical model, "N"=n-body dynamical model
    CHARACTER(len=1), INTENT(in)                        :: in_dynmodel
    ! out_ephemn: output ephemeris, 33 columns per target and epoch:
    ! (1) modified julian date                                                                
    ! (2) right ascension (deg)                                                               
    ! (3) declination (deg)                                                                   
    ! (4) dra/dt sky-motion (deg/day, including cos(dec) factor)                              
    ! (5) ddec/dt sky-motion (deg/day)                                                        
    ! (6) solar phase angle (deg)                                                             
    ! (7) solar elongation angle (deg)                                                        
    ! (8) heliocentric distance (au)                                                          
    ! (9) geocentric distance (au)                                                            
    ! (10) predicted apparent V-band magnitude                                                 
    ! (11) position angle for direction of motion (deg)                                        
    ! (12) topocentric ecliptic longitude (deg)                                                
    ! (13) topocentric ecliptic latitude (deg)                                                 
    ! (14) opposition-centered topocentric ecliptic longitude (deg)                            
    ! (15) opposition-centered topocentric ecliptic latitude (deg)                             
    ! (16) heliocentric ecliptic longitude (deg)                                               
    ! (17) heliocentric ecliptic latitude (deg)                                                
    ! (18) opposition-centered heliocentric ecliptic longitude (deg)                           
    ! (19) opposition-centered heliocentric ecliptic latitude (deg)                            
    ! (20) topocentric object altitude (deg)                                                   
    ! (21) topocentric solar altitude (deg)                                                    
    ! (22) topocentric lunar altitude (deg)                                                    
    ! (23) lunar phase [0...1]                                                                 
    ! (24) lunar elongation (deg, distance between the target and the Moon)                    
    ! (25) heliocentric ecliptic cartesian x coordinate for the object (au)                    
    ! (26) heliocentric ecliptic cartesian y coordinate for the object (au)                    
    ! (27) heliocentric ecliptic cartesian z coordinate for the objects (au)                   
    ! (28) heliocentric ecliptic cartesian x rate for the object (au/day))                     
    ! (29) heliocentric ecliptic cartesian y rate for the object (au/day)                      
    ! (30) heliocentric ecliptic cartesian z rate for the objects (au/day)                     
    ! (31) heliocentric ecliptic cartesian coordinates for the observatory (au)                
    ! (32) heliocentric ecliptic cartesian coordinates for the observatory (au)                
    ! (33) heliocentric ecliptic cartesian coordinates for the observatory (au)
    ! (34) true anomaly
    REAL(8), DIMENSION(in_norb,in_ndate,34), INTENT(out) :: out_ephems ! (1:norb,1:ndate,1:34)
    ! error_code: output error code
    INTEGER, INTENT(out)                                :: error_code

    ! Internal variables.  
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr
    TYPE (Orbit), DIMENSION(1) :: orb_arr
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: observers
    TYPE (CartesianCoordinates) :: ccoord
    TYPE (CartesianCoordinates) :: obsy_ccoord
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides
    TYPE (SphericalCoordinates) :: scoord
    TYPE (Time) :: t
    REAL(8), DIMENSION(:,:), POINTER :: planeph
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    CHARACTER(len=6) :: dyn_model
    REAL(8), DIMENSION(6) :: coordinates, &
         elements, &
         comp_coord, &
         h_ecl_car_coord_obsy, &
         h_ecl_car_coord_obj         
    REAL(8), DIMENSION(3) :: obsy_obj, &
         obsy_pos, &
         geoc_obsy, &
         obsy_sun, &
         vec3, &
         pos, &
         sun_moon, &
         obsy_moon, &
         pos_opp
    REAL(8) :: cos_phase, &
         ephemeris_r2, &
         heliocentric_r2, &
         integration_step, &
         mjd, &
         mjd_tt, &
         observer_r2, &
         phase, &
         solar_elongation, &
         vmag, &
         hlon, &
         hlat, &
         ta_s, &
         ta_c, &
         fak, &
         ecc_anom, &
         true_anom, &
         obsy_moon_r2, &
         cos_obj_phase, &
         lunar_phase, &
         lunar_elongation, &
         lunar_alt, &
         solar_alt, &
         obj_alt, &
         pa, &
         tlon, &
         tlat, &
         opplon, &
         opplat, &
         toclon, &
         toclat, &
         hoclon, &
         hoclat
    INTEGER :: i, &
         j
    LOGICAL, DIMENSION(10) :: perturbers

    ! Init
    errstr = ""
    error_code = 0
    IF (in_dynmodel .EQ. "2") THEN
       dyn_model = "2-body"
    END IF
    IF (in_dynmodel .EQ. "N") THEN
       dyn_model = "n-body"
    END IF
    integrator = "bulirsch-stoer"
    integration_step = 5.0_8
    perturbers = .TRUE.
    ALLOCATE(observers(SIZE(in_date_ephems,dim=1)))
    DO i=1,SIZE(in_date_ephems,dim=1)
       CALL NEW(t, in_date_ephems(i,1), timescales(NINT(in_date_ephems(i,2))))
       IF (error) THEN
          ! Error in creating a new Time object.
          error_code = 35
          RETURN
       END IF
       ! Compute heliocentric observatory coordinates
       observers(i) = getObservatoryCCoord(obsies, in_obscode, t)
       IF (error) THEN
          ! Error in getObservatoryCCoord()
          error_code = 36
          RETURN
       END IF
       CALL rotateToEquatorial(observers(i))
       CALL NULLIFY(t)
    END DO

    ! Loop over orbits:
    DO i=1,SIZE(in_orbits,dim=1)
    
       ! Get the element type from the input flattened orbit.
       IF(NINT(in_orbits(i,8)) < 0 .OR.                                       &
            NINT(in_orbits(i,8)) > SIZE(element_types)) THEN
          ! Error: unsupported orbital elements.
          error_code = 58
          RETURN
       END IF

       ! Get each flattened orbit and create an Orbit instance.
       elements(1:6) = in_orbits(i,2:7)
       ! Create a Time instance.
       CALL NEW(t, in_orbits(i,9), timescales(NINT(in_orbits(i,10))))
       IF(error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an instant of StochasticOrbit via an Orbit instance.
       CALL NEW(orb_arr(1), &
            elements(1:6), &
            element_types(NINT(in_orbits(i,8))), &
            "ecliptic", &
            copy(t))
       CALL NULLIFY(t)
       CALL setParameters(orb_arr(1), &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
          error_code = 37
          RETURN
       END IF
       ! Use cartesian equatorial coordinates:
       CALL toCartesian(orb_arr(1), "equatorial")
       !CALL toCometary(orb_arr(1))
       IF (error) THEN
          ! Error in setParameters()
          error_code = 39
          RETURN
       END IF
       ! Compute topocentric ephemerides
       CALL getEphemerides(orb_arr, &
            observers, &
            ephemerides, &
            this_lt_corr_arr=orb_lt_corr_arr)
       IF (error) THEN
          ! Error in getEphemerides()
          error_code = 40
          RETURN
       END IF
       !CALL NULLIFY(orb_arr(1))
   
       ! Now export the ephem_arr to a flat array.
       DO j=1,SIZE(observers)

          ! Make sure that the ephemeris is equatorial:
          CALL rotateToEquatorial(ephemerides(1,j))
          coordinates = getCoordinates(ephemerides(1,j))
          ! coordinates(1) = Delta
          ! coordinates(2) = ra
          ! coordinates(3) = dec
          ! coordinates(4) = dDelta 
          ! coordinates(5) = dra (to be multiplied with cos(coordinates(3))
          ! coordinates(6) = ddec

          ! Compute position angle for direction of motion
          pa = ATAN2(coordinates(5)*COS(coordinates(3)), coordinates(6))
              IF (pa < 0.0_bp) THEN
                 pa = two_pi + pa
              END IF

          ! Extract topocentric ecliptic lon and lat
          CALL rotateToEcliptic(ephemerides(1,j))        
          comp_coord = getCoordinates(ephemerides(1, j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (40)',1)
             STOP
          END IF
          tlon = comp_coord(2)
          tlat = comp_coord(3)

          ! Compute topocentric opposition coordinates (actually,
          ! get heliocentric observatory coordinates)
          scoord = getSCoord(observers(j))
          CALL rotateToEcliptic(scoord)
          pos_opp = getPosition(scoord)
          opplon = pos_opp(2)
          opplat = pos_opp(3)
          CALL NULLIFY(scoord)

          ! Compute opposition-centered topocentric ecliptic coordinates
          toclon = tlon - opplon
          toclat = tlat - opplat
          IF (toclon > pi) THEN
             toclon = toclon - two_pi
          ELSE IF (toclon < -pi) THEN
             toclon = toclon + two_pi
          END IF

          ! Extract heliocentric ecliptic lon and lat
          ccoord = getCCoord(orb_lt_corr_arr(1,j), "ecliptic")
          scoord = getSCoord(ccoord)
          comp_coord = getCoordinates(scoord)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (45)',1)
             STOP
          END IF
          hlon = comp_coord(2)
          hlat = comp_coord(3)
          CALL NULLIFY(ccoord)
          CALL NULLIFY(scoord)

          ! Compute opposition-centered heliocentric ecliptic coordinates
          hoclon = hlon - opplon
          hoclat = hlat - opplat
          IF (hoclon > pi) THEN
             hoclon = hoclon - two_pi
          ELSE IF (hoclon < -pi) THEN
             hoclon = hoclon + two_pi
          END IF

          ! r, Delta, phase angle
          CALL NEW(ccoord, ephemerides(1,j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (50)',1)
             STOP
          END IF
          CALL rotateToEquatorial(ccoord)
          obsy_obj = getPosition(ccoord)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (55)',1)
             STOP
          END IF
          CALL NULLIFY(ccoord)
          ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
          CALL toCartesian(orb_lt_corr_arr(1,j), frame='equatorial')
          pos = getPosition(orb_lt_corr_arr(1,j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (60)',1)
             STOP
          END IF
          heliocentric_r2 = DOT_PRODUCT(pos,pos)
          obsy_pos = getPosition(observers(j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (65)',1)
             STOP
          END IF
          ! geoc_obsy = obsy_pos - pos  not in oorb.f90
          observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)
          cos_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
               observer_r2) / (SQRT(heliocentric_r2 * ephemeris_r2))
          phase = ACOS(cos_phase)

          ! apparent brightness
          IF (in_orbits(i,8) .EQ. 2) THEN
             ! if the target is a comet (elements type == 'COM')
             vmag = in_orbits(i,11)+5*LOG10(coordinates(1))+2.5*in_orbits(i,12)*LOG10(SQRT(heliocentric_r2))
          ELSE
             ! if the target is not a comet (any other elements type)
             vmag = getApparentHGMagnitude(H=in_orbits(i,11), &
               G=in_orbits(i,12), r=SQRT(heliocentric_r2), &
               Delta=coordinates(1), phase_angle=phase)
          END IF

          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (70)',1)
             STOP
          END IF
         
          ! Parameters relevant for Earth-based observers

          ! ephem date
          t = getTime(observers(j))
          mjd = getMJD(t, timescales(NINT(in_date_ephems(j,2))))
          mjd_tt = getMJD(t, "TT")

          obsy_ccoord = getGeocentricObservatoryCCoord(obsies, in_obscode, t)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (75)',1)
             STOP
          END IF

          CALL rotateToEquatorial(obsy_ccoord)        
          geoc_obsy = getPosition(obsy_ccoord)

          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (80)',1)
             STOP
          END IF

          CALL NULLIFY(t)
          
          vec3 = cross_product(geoc_obsy,obsy_obj)
          obj_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_obj))
          
          ! Compute (approximate) altitude of the Sun
          ! Position of the geocenter as seen from the Sun:
          planeph => JPL_ephemeris(mjd_tt, 3, 11, error)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (85)',1)
             STOP
          END IF
          ! Position of the Sun as seen from the observatory:
          obsy_sun = -(planeph(1,1:3) + geoc_obsy)
          DEALLOCATE(planeph)
          vec3 = cross_product(geoc_obsy,obsy_sun)
          solar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_sun))

          ! Compute the solar elongation:
          vec3 = cross_product(obsy_obj,obsy_sun)
          solar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_sun))
          
          ! Compute phase of the Moon:
          ! Position of the Moon as seen from the Sun:
          planeph => JPL_ephemeris(mjd_tt, 10, 11, error)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (95)',1)
             STOP
          END IF
          sun_moon = planeph(1,1:3)
          DEALLOCATE(planeph)
          ! Angle between Sun and Moon as seen from the observatory:
          obsy_moon = sun_moon - obsy_pos
          obsy_moon_r2 = DOT_PRODUCT(obsy_moon,obsy_moon)
          cos_obj_phase = DOT_PRODUCT(obsy_moon,-obsy_pos) / &
               (SQRT(observer_r2) * SQRT(obsy_moon_r2))
          lunar_phase = (1.0_bp-cos_obj_phase)/2.0_bp

          ! Compute (approximate) distance between the target and the Moon:
          vec3 = cross_product(obsy_obj,obsy_moon)
          lunar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_moon))

          ! Compute (approximate) altitude of the Moon:
          vec3 = cross_product(geoc_obsy,obsy_moon)
          lunar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_moon))

          ! Extract heliocentric ecliptic cartesian coordinates for the object
          h_ecl_car_coord_obj = getElements(orb_lt_corr_arr(1,j), "cartesian", "ecliptic")
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (90)',1)
             STOP
          END IF

          ! Extract heliocentric ecliptic cartesian coordinates for the observer
          CALL rotateToEcliptic(observers(j))
          h_ecl_car_coord_obsy = getCoordinates(observers(j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (95)',1)
             STOP
          END IF

          ! calculate true anomaly
          elements = getElements(orb_arr(1), "cometary")
          IF (elements(2) < 1.0_bp) THEN
             t = getTime(observers(j))
             ! Compute eccentric and true anomalies if orbit is elliptic:
             CALL solveKeplerEquation(orb_arr(1), t, ecc_anom)
             CALL NULLIFY(t)
             ta_s = SIN(ecc_anom)
             ta_c = COS(ecc_anom)
             fak = SQRT(1 - elements(2) * elements(2))
             true_anom = MODULO(ATAN2(fak * ta_s, ta_c - elements(2)),two_pi)
          ELSE
             ecc_anom = -99.0_bp
             true_anom = -99.0_bp
          END IF
          
          ! Write the output ephem array.
          out_ephems(i,j,1) = mjd                                        ! modified julian date
          out_ephems(i,j,2) = coordinates(2)/rad_deg                     ! right ascension (deg)
          out_ephems(i,j,3) = coordinates(3)/rad_deg                     ! declination (deg)
          out_ephems(i,j,4) = coordinates(5)*COS(coordinates(3))/rad_deg ! dra/dt sky-motion (deg/day, including cos(dec) factor)
          out_ephems(i,j,5) = coordinates(6)/rad_deg                     ! ddec/dt sky-motion (deg/day)
          out_ephems(i,j,6) = phase/rad_deg                              ! solar phase angle (deg)
          out_ephems(i,j,7) = solar_elongation/rad_deg                   ! solar elongation angle (deg)
          out_ephems(i,j,8) = SQRT(heliocentric_r2)                      ! heliocentric distance (au)
          out_ephems(i,j,9) = coordinates(1)                             ! geocentric distance (au)
          out_ephems(i,j,10) = vmag                                      ! predicted apparent V-band magnitude
          out_ephems(i,j,11) = pa/rad_deg                                ! position angle for direction of motion (deg)
          out_ephems(i,j,12) = tlon/rad_deg                              ! topocentric ecliptic longitude (deg)
          out_ephems(i,j,13) = tlat/rad_deg                              ! topocentric ecliptic latitude (deg)
          out_ephems(i,j,14) = toclon/rad_deg                            ! opposition-centered topocentric ecliptic longitude (deg)
          out_ephems(i,j,15) = toclat/rad_deg                            ! opposition-centered topocentric ecliptic latitude (deg)
          out_ephems(i,j,16) = hlon/rad_deg                              ! heliocentric ecliptic longitude (deg)
          out_ephems(i,j,17) = hlat/rad_deg                              ! heliocentric ecliptic latitude (deg)
          out_ephems(i,j,18) = hoclon/rad_deg                            ! opposition-centered heliocentric ecliptic longitude (deg)
          out_ephems(i,j,19) = hoclat/rad_deg                            ! opposition-centered heliocentric ecliptic latitude (deg)
          out_ephems(i,j,20) = obj_alt/rad_deg                           ! topocentric object altitude (deg)
          out_ephems(i,j,21) = solar_alt/rad_deg                         ! topocentric solar altitude (deg)
          out_ephems(i,j,22) = lunar_alt/rad_deg                         ! topocentric lunar altitude (deg)
          out_ephems(i,j,23) = lunar_phase                               ! lunar phase [0...1]
          out_ephems(i,j,24) = lunar_elongation/rad_deg                  ! lunar elongation (deg, distance between the target and the Moon)
          out_ephems(i,j,25) = h_ecl_car_coord_obj(1)                    ! heliocentric ecliptic cartesian x coordinate for the object (au)
          out_ephems(i,j,26) = h_ecl_car_coord_obj(2)                    ! heliocentric ecliptic cartesian y coordinate for the object (au)
          out_ephems(i,j,27) = h_ecl_car_coord_obj(3)                    ! heliocentric ecliptic cartesian z coordinate for the objects (au)
          out_ephems(i,j,28) = h_ecl_car_coord_obj(4)                    ! heliocentric ecliptic cartesian x rate for the object (au/day))
          out_ephems(i,j,29) = h_ecl_car_coord_obj(5)                    ! heliocentric ecliptic cartesian y rate for the object (au/day)
          out_ephems(i,j,30) = h_ecl_car_coord_obj(6)                    ! heliocentric ecliptic cartesian z rate for the objects (au/day)
          out_ephems(i,j,31) = h_ecl_car_coord_obsy(1)                   ! heliocentric ecliptic cartesian coordinates for the observatory (au)
          out_ephems(i,j,32) = h_ecl_car_coord_obsy(2)                   ! heliocentric ecliptic cartesian coordinates for the observatory (au)
          out_ephems(i,j,33) = h_ecl_car_coord_obsy(3)                   ! heliocentric ecliptic cartesian coordinates for the observatory (au)
          out_ephems(i,j,34) = true_anom/rad_deg                         ! true anomaly (deg)                                                           

          CALL NULLIFY(ephemerides(1,j))
          CALL NULLIFY(orb_lt_corr_arr(1,j))
          
       END DO
       CALL NULLIFY(orb_arr(1))
       DEALLOCATE(ephemerides, orb_lt_corr_arr)

    END DO
    DO i=1,SIZE(observers)
      CALL NULLIFY(observers(i))
    END DO
    DEALLOCATE(observers)

  END SUBROUTINE oorb_ephemeris_full


  SUBROUTINE oorb_ephemeris_basic(in_norb, &
       in_orbits,                    &
       in_obscode,                   &
       in_ndate,                     &
       in_date_ephems,               &
       in_dynmodel,                  &
       out_ephems,                   &
       error_code)
 
    ! in_norb: number of input orbits
    INTEGER, INTENT(in)                                 :: in_norb
    ! in_orbits: input flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)           :: in_orbits ! (1:norb,1:12)
    ! in_obscode: observatory code as defined by the Minor Planet Center
    CHARACTER(len=*), INTENT(in)                        :: in_obscode
    ! in_ndate: number of input epochs
    INTEGER, INTENT(in)                                 :: in_ndate
    ! in_date_ephems: epochs for ephemeris:
    ! (1) modified Julian date
    ! (2) timescale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    REAL(8), DIMENSION(in_ndate,2), INTENT(in)          :: in_date_ephems ! (1:ndate,1:2)
    ! in_model: "2"=2-body dynamical model, "N"=n-body dynamical model
    CHARACTER(len=1), INTENT(in)                        :: in_dynmodel
    ! out_ephemn: output ephemeris, 33 columns per target and epoch:
    ! (1) modified julian date                                                                
    ! (2) right ascension (deg)                                                               
    ! (3) declination (deg)                                                                   
    ! (4) dra/dt sky-motion (deg/day, including cos(dec) factor)                              
    ! (5) ddec/dt sky-motion (deg/day)                                                        
    ! (6) solar phase angle (deg)                                                             
    ! (7) solar elongation angle (deg)                                                        
    ! (8) heliocentric distance (au)                                                          
    ! (9) geocentric distance (au)                                                            
    ! (10) predicted apparent V-band magnitude
    ! (11) true anomaly

    
    REAL(8), DIMENSION(in_norb,in_ndate,11), INTENT(out) :: out_ephems ! (1:norb,1:ndate,1:11)
    ! error_code: output error code
    INTEGER, INTENT(out)                                :: error_code

    ! Internal variables.  
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr
    TYPE (Orbit), DIMENSION(1) :: orb_arr
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: observers
    TYPE (CartesianCoordinates) :: ccoord
    TYPE (CartesianCoordinates) :: obsy_ccoord
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides
    TYPE (SphericalCoordinates) :: scoord
    TYPE (Time) :: t
    REAL(8), DIMENSION(:,:), POINTER :: planeph
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    CHARACTER(len=6) :: dyn_model
    REAL(8), DIMENSION(6) :: coordinates, &
         elements, &
         comp_coord, &
         h_ecl_car_coord_obsy, &
         h_ecl_car_coord_obj         
    REAL(8), DIMENSION(3) :: obsy_obj, &
         obsy_pos, &
         geoc_obsy, &
         obsy_sun, &
         vec3, &
         pos
    REAL(8) :: cos_phase, &
         ephemeris_r2, &
         heliocentric_r2, &
         integration_step, &
         mjd, &
         mjd_tt, &
         observer_r2, &
         phase, &
         solar_elongation, &
         vmag, &
         ta_s, &
         ta_c, &
         fak, &
         ecc_anom, &
         true_anom
    INTEGER :: i, &
         j
    LOGICAL, DIMENSION(10) :: perturbers

    ! Init
    errstr = ""
    error_code = 0
    IF (in_dynmodel .EQ. "2") THEN
       dyn_model = "2-body"
    END IF
    IF (in_dynmodel .EQ. "N") THEN
       dyn_model = "n-body"
    END IF
    integrator = "bulirsch-stoer"
    integration_step = 5.0_8
    perturbers = .TRUE.
    ALLOCATE(observers(SIZE(in_date_ephems,dim=1)))
    DO i=1,SIZE(in_date_ephems,dim=1)
       CALL NEW(t, in_date_ephems(i,1), timescales(NINT(in_date_ephems(i,2))))
       IF (error) THEN
          ! Error in creating a new Time object.
          error_code = 35
          RETURN
       END IF
       ! Compute heliocentric observatory coordinates
       observers(i) = getObservatoryCCoord(obsies, in_obscode, t)
       IF (error) THEN
          ! Error in getObservatoryCCoord()
          error_code = 36
          RETURN
       END IF
       CALL rotateToEquatorial(observers(i))
       CALL NULLIFY(t)
    END DO

    ! Loop over orbits:
    DO i=1,SIZE(in_orbits,dim=1)
    
       ! Get the element type from the input flattened orbit.
       IF(NINT(in_orbits(i,8)) < 0 .OR.                                       &
            NINT(in_orbits(i,8)) > SIZE(element_types)) THEN
          ! Error: unsupported orbital elements.
          error_code = 58
          RETURN
       END IF

       ! Get each flattened orbit and create an Orbit instance.
       elements(1:6) = in_orbits(i,2:7)
       ! Create a Time instance.
       CALL NEW(t, in_orbits(i,9), timescales(NINT(in_orbits(i,10))))
       IF(error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an instant of StochasticOrbit via an Orbit instance.
       CALL NEW(orb_arr(1), &
            elements(1:6), &
            element_types(NINT(in_orbits(i,8))), &
            "ecliptic", &
            copy(t))
       CALL NULLIFY(t)
       CALL setParameters(orb_arr(1), &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
          error_code = 37
          RETURN
       END IF
       ! Use cartesian equatorial coordinates:
       CALL toCartesian(orb_arr(1), "equatorial")
       !CALL toCometary(orb_arr(1))
       IF (error) THEN
          ! Error in setParameters()
          error_code = 39
          RETURN
       END IF
       ! Compute topocentric ephemerides
       CALL getEphemerides(orb_arr, &
            observers, &
            ephemerides, &
            this_lt_corr_arr=orb_lt_corr_arr)
       IF (error) THEN
          ! Error in getEphemerides()
          error_code = 40
          RETURN
       END IF

   
       ! Now export the ephem_arr to a flat array.
       DO j=1,SIZE(observers)

          ! Make sure that the ephemeris is equatorial:
          CALL rotateToEquatorial(ephemerides(1,j))
          coordinates = getCoordinates(ephemerides(1,j))
          ! coordinates(1) = Delta
          ! coordinates(2) = ra
          ! coordinates(3) = dec
          ! coordinates(4) = dDelta 
          ! coordinates(5) = dra (to be multiplied with cos(coordinates(3))
          ! coordinates(6) = ddec

          ! ephem date
          t = getTime(observers(j))
          mjd = getMJD(t, timescales(NINT(in_date_ephems(j,2))))
          mjd_tt = getMJD(t, "TT")
          CALL NULLIFY(t)

          ! r, Delta, phase angle
          CALL NEW(ccoord, ephemerides(1,j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (50)',1)
             STOP
          END IF
          CALL rotateToEquatorial(ccoord)
          obsy_obj = getPosition(ccoord)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (55)',1)
             STOP
          END IF
          CALL NULLIFY(ccoord)
          ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
          CALL toCartesian(orb_lt_corr_arr(1,j), frame='equatorial')
          pos = getPosition(orb_lt_corr_arr(1,j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (60)',1)
             STOP
          END IF
          heliocentric_r2 = DOT_PRODUCT(pos,pos)
          obsy_pos = getPosition(observers(j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (65)',1)
             STOP
          END IF
          ! geoc_obsy = obsy_pos - pos  not in oorb.f90
          observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)
          cos_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
               observer_r2) / (SQRT(heliocentric_r2 * ephemeris_r2))
          phase = ACOS(cos_phase)

          ! apparent brightness
          IF (in_orbits(i,8) .EQ. 2) THEN
             ! if the target is a comet (elements type == 'COM')
             vmag = in_orbits(i,11)+5*LOG10(coordinates(1))+2.5*in_orbits(i,12)*LOG10(SQRT(heliocentric_r2))
          ELSE
             ! if the target is not a comet (any other elements type)
             vmag = getApparentHGMagnitude(H=in_orbits(i,11), &
               G=in_orbits(i,12), r=SQRT(heliocentric_r2), &
               Delta=coordinates(1), phase_angle=phase)
          END IF

          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (70)',1)
             STOP
          END IF
         
          ! Compute (approximate) altitude of the Sun
          ! Position of the geocenter as seen from the Sun:
          planeph => JPL_ephemeris(mjd_tt, 3, 11, error)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (85)',1)
             STOP
          END IF
          ! Position of the Sun as seen from the observatory:
          obsy_sun = -(planeph(1,1:3) + geoc_obsy)
          DEALLOCATE(planeph)

          ! Compute the solar elongation:
          vec3 = cross_product(obsy_obj,obsy_sun)
          solar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_sun))

          ! calculate true anomaly
          elements = getElements(orb_arr(1), "cometary")
          IF (elements(2) < 1.0_bp) THEN
             t = getTime(observers(j))
             ! Compute eccentric and true anomalies if orbit is elliptic:
             CALL solveKeplerEquation(orb_arr(1), t, ecc_anom)
             CALL NULLIFY(t)
             ta_s = SIN(ecc_anom)
             ta_c = COS(ecc_anom)
             fak = SQRT(1 - elements(2) * elements(2))
             true_anom = MODULO(ATAN2(fak * ta_s, ta_c - elements(2)),two_pi)
          ELSE
             ecc_anom = -99.0_bp
             true_anom = -99.0_bp
          END IF
          
          ! Write the output ephem array.
          out_ephems(i,j,1) = mjd                                        ! modified julian date
          out_ephems(i,j,2) = coordinates(2)/rad_deg                     ! right ascension (deg)
          out_ephems(i,j,3) = coordinates(3)/rad_deg                     ! declination (deg)
          out_ephems(i,j,4) = coordinates(5)*COS(coordinates(3))/rad_deg ! dra/dt sky-motion (deg/day, including cos(dec) factor)
          out_ephems(i,j,5) = coordinates(6)/rad_deg                     ! ddec/dt sky-motion (deg/day)
          out_ephems(i,j,6) = phase/rad_deg                              ! solar phase angle (deg)
          out_ephems(i,j,7) = solar_elongation/rad_deg                   ! solar elongation angle (deg)
          out_ephems(i,j,8) = SQRT(heliocentric_r2)                      ! heliocentric distance (au)
          out_ephems(i,j,9) = coordinates(1)                             ! geocentric distance (au)
          out_ephems(i,j,10) = vmag                                      ! predicted apparent V-band magnitude
          out_ephems(i,j,11) = true_anom/rad_deg                         ! true anomaly (deg)                                                           
          CALL NULLIFY(ephemerides(1,j))
          CALL NULLIFY(orb_lt_corr_arr(1,j))

          WRITE (*,*) true_anom/rad_deg
          
       END DO
       CALL NULLIFY(orb_arr(1))
       DEALLOCATE(ephemerides, orb_lt_corr_arr)

    END DO
    DO i=1,SIZE(observers)
      CALL NULLIFY(observers(i))
    END DO
    DEALLOCATE(observers)

  END SUBROUTINE oorb_ephemeris_basic
  

  SUBROUTINE oorb_ephemeris_covariance(in_norb, &
       in_orbits,                    &
       in_covariances,               &
       in_obscode,                   &
       in_ndate,                     &
       in_date_ephems,               &
       out_ephems,                   &
       error_code)

    ! this function is working but has not yet been updated

    ! Lynne's comment: "I would also like to request an update to
    ! oorb_ephemeris_covariance to match oorb_ephemeris_basic
    ! (i.e. return the same set of values, although now they'd include
    ! the uncertainties). Since it looks like only the Ra/Dec values
    ! have errors, maybe it's reasonable to just update the returned
    ! values so that they look the same as the values returned by
    ! oorb_ephemeris_basic, but then still add the uncertainty ellipse
    ! at the end."
    
    ! in_norb: number of input orbits
    INTEGER, INTENT(in)                                  :: in_norb
    ! in_orbits: input flattened orbits, 12 columns per target:
    ! (1) object id (integer value)
    ! (2-7) orbital elements:
    !   * (q, e, i, longnode, argper, perihelion epoch in mjd) for comet format
    !   * (a, e, i, longnode, argper, mean anomaly) for keplerian format
    !   * (x, y, z, dx, dy, dz) for cartesian format
    ! (8) orbital elements type ('CART': 1, 'COM': 2, 'KEP': 3, 'DEL': 4, 'EQX': 5)
    ! (9) epoch in mjd
    ! (10) time scale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    ! (11) target absolute magnitude H/M1 parameter for comets ('COM' elements type)
    ! (12) target slope parameter G/K1 parameter for comets ('COM' elements type)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)             :: in_orbits ! (1:norb,1:12)
    ! in_covariances: input covariance matrices for orbital elements
    REAL(8), DIMENSION(in_norb,6,6), INTENT(in)           :: in_covariances ! (1:norb,1:6,1:6)
    ! in_obscode: observatory code as defined by the Minor Planet Center
    CHARACTER(len=*), INTENT(in)                          :: in_obscode
    ! in_ndate: number of input epochs
    INTEGER, INTENT(in)                                   :: in_ndate
    ! in_date_ephems: epochs for ephemeris:
    ! (1) modified Julian date
    ! (2) timescale type ('UTC': 1, 'UT1': 2, 'TT': 3, 'TAI': 4)
    REAL(8), DIMENSION(in_ndate,2), INTENT(in)            :: in_date_ephems ! (1:ndate,1:2)
    ! out_ephem: output ephemeris, 14 columns per target and epoch:
    ! (1) distance from the observer (au)
    ! (2) right ascension (deg)                                                               
    ! (3) declination (deg)                                                                   
    ! (4) predicted apparent V-band magnitude                                                 
    ! (5) epoch (modified Julian date) 
    ! (6) timescale type
    ! (7) dra/dt sky-motion (deg/day, including cos(dec) factor)                              
    ! (8) ddec/dt sky-motion (deg/day)                                                        
    ! (9) solar phase angle (deg)                                                             
    ! (10) ra uncertainty (arcsec)
    ! (11) dec uncertainty (arcsec)
    ! (12) semi-major axis of the uncertainty ellipse (arcmin)
    ! (13) semi-minor axis of the uncertainty ellipse (arcmin)
    ! (14) position angle of the uncertainty ellipse (deg)
    REAL(8), DIMENSION(in_norb,in_ndate,14), INTENT(out)  :: out_ephems ! (1:norb,1:ndate,1:14)
    ! error_code: output error code
    INTEGER, INTENT(out)                                  :: error_code

    ! Internal variables.  
    TYPE (StochasticOrbit) :: storb
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr
    TYPE (Orbit), DIMENSION(1) :: orb_arr
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: observers
    TYPE (CartesianCoordinates) :: ccoord
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides
    REAL(bp), DIMENSION(:,:), POINTER :: planeph
    TYPE (Time) :: t
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    CHARACTER(len=6) :: dyn_model
    REAL(8), DIMENSION(:,:,:), POINTER :: cov_arr
    REAL(8), DIMENSION(6,6) :: eigenvec
    REAL(8), DIMENSION(6) :: coordinates, &
         eigenval, &
         elements
    REAL(8), DIMENSION(3) :: obsy_obj, &
         obsy_pos, &
         pos
    REAL(8) :: cos_phase, &
         ephemeris_r2, &
         heliocentric_r2, &
         integration_step, &
         mjd, &
         observer_r2, &
         pa, &
         phase, &
         sigma_dec, &
         sigma_ra, &
         smaa, &
         smia, &
         vmag
    INTEGER :: ismaa, &
         ismia, &
         i, &
         j, &
         nrot
    LOGICAL, DIMENSION(10) :: perturbers

    ! Init
    errstr = ""
    error_code = 0
    eigenval = 0.0_8
    eigenvec = 0.0_8
    dyn_model = "n-body"
    integrator = "bulirsch-stoer"
    integration_step = 5.0_8
    perturbers = .TRUE.
    ALLOCATE(observers(SIZE(in_date_ephems,dim=1)))
    DO i=1,SIZE(in_date_ephems,dim=1)
       CALL NEW(t, in_date_ephems(i,1), timescales(NINT(in_date_ephems(i,2))))
       IF (error) THEN
          ! Error in creating a new Time object.
          error_code = 35
          RETURN
       END IF
       ! Compute heliocentric observatory coordinates
       observers(i) = getObservatoryCCoord(obsies, in_obscode, t)
       IF (error) THEN
          ! Error in getObservatoryCCoord()
          error_code = 36
          RETURN
       END IF
       CALL rotateToEquatorial(observers(i))
       CALL NULLIFY(t)
    END DO

    ! Loop over orbits:
    DO i=1,SIZE(in_orbits,dim=1)

       ! Get the element type from the input flattened orbit.
       IF(NINT(in_orbits(i,8)) < 0 .OR.                                       &
            NINT(in_orbits(i,8)) > SIZE(element_types)) THEN
          ! Error: unsupported orbital elements.
          error_code = 58
          RETURN
       END IF

       ! Get each flattened orbit and create an Orbit instance.
       elements(1:6) = in_orbits(i,2:7)
       ! Create a Time instance.
       CALL NEW(t, in_orbits(i,9), timescales(NINT(in_orbits(i,10))))
       IF(error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an instant of StochasticOrbit via an Orbit instance.
       CALL NEW(orb_arr(1), &
            elements(1:6), &
            element_types(NINT(in_orbits(i,8))), &
            "ecliptic", &
            copy(t))
       CALL NULLIFY(t)
       CALL setParameters(orb_arr(1), &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
          error_code = 37
          RETURN
       END IF
       CALL NEW(storb, orb_arr(1), in_covariances(i,:,:), &
            cov_type=element_types(NINT(in_orbits(i,8))), &
            element_type=element_types(NINT(in_orbits(i,8))))
       IF (error) THEN
          error_code = 34
          RETURN
       END IF
       ! Use cartesian equatorial coordinates:
       CALL toCartesian(storb, "equatorial")
       !CALL toCometary(storb)
       IF (error) THEN
          ! Error in setParameters()
          error_code = 36
          RETURN
       END IF

       ! Set integration parameters
       CALL setParameters(storb, &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
          ! Error in setParameters()
          error_code = 37
          RETURN
       END IF

       ! Compute topocentric ephemerides
       CALL getEphemerides(storb, &
            observers, &
            ephemerides, &
            cov_arr=cov_arr, &
            this_lt_corr_arr=orb_lt_corr_arr)
       IF (error) THEN
          ! Error in getEphemerides()
          error_code = 38
          RETURN
       END IF
       CALL NULLIFY(storb)
       CALL NULLIFY(orb_arr(1))

       ! Now export the ephem_arr to a flat array.
       DO j=1,SIZE(observers)

          ! Make sure that the ephemeris is equatorial:
          CALL rotateToEquatorial(ephemerides(1,j))
          ! Extract RA & Dec
          coordinates = getCoordinates(ephemerides(1,j))

          ! Calculate apparent brightness
          CALL NEW(ccoord, ephemerides(1,j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (50)',1)
             STOP
          END IF
          CALL rotateToEquatorial(ccoord)
          obsy_obj = getPosition(ccoord)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (55)',1)
             STOP
          END IF
          CALL NULLIFY(ccoord)
          ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
          CALL toCartesian(orb_lt_corr_arr(1,j), frame='equatorial')
          pos = getPosition(orb_lt_corr_arr(1,j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (60)',1)
             STOP
          END IF
          heliocentric_r2 = DOT_PRODUCT(pos,pos)
          obsy_pos = getPosition(observers(j))
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (65)',1)
             STOP
          END IF
          observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)
          cos_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
               observer_r2) / (SQRT(heliocentric_r2) * &
               SQRT(ephemeris_r2))
          phase = ACOS(cos_phase)
          vmag = getApparentHGMagnitude(H=in_orbits(i,11), &
               G=in_orbits(i,12), r=SQRT(heliocentric_r2), &
               Delta=coordinates(1), phase_angle=phase)
          IF (error) THEN
             CALL errorMessage('oorb / ephemeris', &
                  'TRACE BACK (70)',1)
             STOP
          END IF

          ! ephem date
          t = getTime(observers(j))
          mjd = getMJD(t, timescales(NINT(in_date_ephems(j,2))))
          CALL NULLIFY(t)


          cov_arr(:,:,j) = SIGN(SQRT(ABS(cov_arr(:,:,j))),cov_arr(:,:,j))
          sigma_ra = cov_arr(2,2,j)/rad_asec
          sigma_dec = cov_arr(3,3,j)/rad_asec
          CALL eigen_decomposition_jacobi(cov_arr(2:3,2:3,j), &
               eigenval(2:3), eigenvec(2:3,2:3), nrot, errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             error_code = 7
             WRITE(*,*) TRIM(errstr)
             RETURN
          END IF
          ismaa = 1 + MAXLOC(ABS(eigenval(2:3)),dim=1)
          ismia = 1 + MINLOC(ABS(eigenval(2:3)),dim=1)
          smaa = ABS(eigenval(ismaa))/rad_amin
          smia = ABS(eigenval(ismia))/rad_amin
          pa = ATAN2(eigenvec(3,ismaa),eigenvec(2,ismaa))/rad_deg
          ! Write the output ephem array.
          out_ephems(i,j,1) = coordinates(1)                             ! distance
          out_ephems(i,j,2) = coordinates(2)/rad_deg                     ! ra
          out_ephems(i,j,3) = coordinates(3)/rad_deg                     ! dec
          out_ephems(i,j,4) = vmag                                       ! mag
          out_ephems(i,j,5) = mjd                                        ! ephem mjd
          out_ephems(i,j,6) = NINT(in_date_ephems(j,2))                  ! ephem mjd timescale
          out_ephems(i,j,7) = coordinates(5)*COS(coordinates(3))/rad_deg ! dra/dt  sky-motion
          out_ephems(i,j,8) = coordinates(6)/rad_deg                     ! ddec/dt sky-motion
          out_ephems(i,j,9) = phase/rad_deg                              ! phase angle
          out_ephems(i,j,10) = sigma_ra                                  ! raErr
          out_ephems(i,j,11) = sigma_dec                                 ! decErr
          out_ephems(i,j,12) = smaa                                      ! semi-major axis
          out_ephems(i,j,13) = smia                                      ! semi-minor axis
          out_ephems(i,j,14) = pa                                        ! position angle

          CALL NULLIFY(ephemerides(1,j))
          CALL NULLIFY(orb_lt_corr_arr(1,j))

       END DO

       DEALLOCATE(ephemerides, orb_lt_corr_arr, cov_arr)

    END DO
    DO i=1,SIZE(observers)
       CALL NULLIFY(observers(i))
    END DO
    DEALLOCATE(observers)

  END SUBROUTINE oorb_ephemeris_covariance



END MODULE pyoorb
