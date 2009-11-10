!====================================================================!
!                                                                    !
! Copyright 2002,2003,2004,2005,2006,2007,2008,2009                  !
! Mikael Granvik, Jenni Virtanen, Karri Muinonen, Teemu Laakso,      !
! Dagmara Oszkiewicz                                                 !
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
!! @author  MG, FP
!! @version 2009-11-09
!!
MODULE pyoorb

  USE planetary_data
  USE Base_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  USE Observatories_cl
  USE Orbit_cl
  USE StochasticOrbit_cl
  IMPLICIT NONE
  CHARACTER(len=11), DIMENSION(6), PUBLIC :: element_types = (/             &
       "cartesian  ",&
       "cometary   ",&
       "keplerian  ",&
       "delaunay   ",&
       "poincare   ",&
       "equinoctial" &
       /)
  CHARACTER(len=8), PARAMETER :: internal_timescale = "TAI"


CONTAINS

  SUBROUTINE oorb_init(ephemeris_fname, error_verbosity, &
       info_verbosity, error_code)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: ephemeris_fname
    INTEGER, INTENT(in), OPTIONAL :: error_verbosity
    INTEGER, INTENT(in), OPTIONAL :: info_verbosity
    INTEGER, INTENT(out) :: error_code

    TYPE (Observatories) :: obsies
    TYPE (Time) :: t


    ! Set path to data files:
    CALL setAccessToDataFiles()
    CALL JPL_ephemeris_init(error, &
         filename=TRIM(OORB_DATA_DIR) // "/" // TRIM(EPH_FNAME)) 
    IF (error) THEN
       CALL errorMessage("oorb", &
            "Could not initialize planetary ephemerides using the " // &
            TRIM(OORB_DATA_DIR) // "/" // TRIM(EPH_FNAME) // " file.", 1)
       STOP
    END IF
    ! Read OBSCODE.dat
    CALL NEW(obsies)
    ! Read TAI-UTC.dat and ET-UT.dat
    CALL NEW(t)
    ! Read de405.dat
    IF (PRESENT(ephemeris_fname)) THEN
       CALL JPL_ephemeris_init(error, ephemeris_fname)
    ELSE
       CALL JPL_ephemeris_init(error, ephemeris_fname)
    END IF
    IF (error) THEN
       error_code = 1
       RETURN
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



  SUBROUTINE propagation_2b(elements, mjd_tt0, mjd_tt1)

    IMPLICIT NONE
    REAL(8), DIMENSION(:,:), INTENT(inout) :: elements
    REAL(8), INTENT(in) :: mjd_tt0, mjd_tt1

    TYPE (Orbit) :: orb
    TYPE (Time) :: t0, t1

    CALL NEW(t0, mjd_tt0, "TT")
    CALL NEW(t1, mjd_tt1, "TT")
    CALL NEW(orb, elements(1,:), "cartesian", "ecliptic", t0)
    CALL propagate(orb, t1)
    elements(1,:) = getElements(orb, "cartesian", "ecliptic")
    CALL NULLIFY(t0)
    CALL NULLIFY(t1)
    CALL NULLIFY(orb)

  END SUBROUTINE propagation_2b



  SUBROUTINE oorb_ephemeris(in_orbits, &
       in_covariances,                 &
       in_obscode,                     &
       in_date_ephems,                 &
       out_ephems,                     &
       error_code)

    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    REAL(8),DIMENSION(:,:), INTENT(in)                  :: in_orbits ! (1:norb,1:11)
    ! Uncertainty matrices:
    REAL(8), DIMENSION(:,:,:), INTENT(in)               :: in_covariances ! (1:norb,1:6,1:6)
    ! Observatory code as defined by the Minor Planet Center
    CHARACTER(len=4), INTENT(in)                        :: in_obscode
    ! Ephemeris dates.
    REAL(8), DIMENSION(:), INTENT(in)                   :: in_date_ephems ! (1:ndate)
    ! Output ephemeris
    ! out_ephems = ((dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa), )
    !REAL(8), DIMENSION(SIZE(in_orbits,dim=1), &
    !     SIZE(in_date_ephems,dim=1),10), INTENT(out)    :: out_ephems ! (1:norb,1:ndate,1:10)
    REAL(8), DIMENSION(:,:,:), INTENT(out)              :: out_ephems ! (1:norb,1:ndate,1:10)
    ! Output error code
    INTEGER, INTENT(out)                                :: error_code

    ! Internal variables.  
    TYPE (StochasticOrbit) :: storb
    TYPE (Orbit) :: orb
    TYPE (Observatories) :: obsies
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: observers
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides
    TYPE (Time) :: t
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    CHARACTER(len=6) :: dyn_model
    REAL(8), DIMENSION(:,:,:), POINTER :: cov_arr
    REAL(8), DIMENSION(6,6) :: corr
    REAL(8), DIMENSION(6) :: coordinates, &
         elements, &
         stdev
    REAL(8) :: integration_step, &
         mjd
    INTEGER :: i, &
         j, &
         k, &
         l
    LOGICAL, DIMENSION(10) :: perturbers

    ! Init
    error_code = 0
    dyn_model = "n-body"
    integration_step = 5.0_8
    perturbers = .TRUE.
    CALL NEW(obsies)
    ALLOCATE(observers(SIZE(in_date_ephems)))
    DO i=1,SIZE(in_date_ephems)
       CALL NEW(t, in_date_ephems(i), internal_timescale)
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
       IF(NINT(in_orbits(i,11)) < 0 .OR.                                       &
            NINT(in_orbits(i,11)) > SIZE(element_types)) THEN
          ! Error: unsupported orbital elements.
          error_code = 58
          RETURN
       END IF

       ! Get each flattened orbit and create an Orbit instance.
       ! Just to beat on a dead horse:
       ! in_orbits(1):      id
       ! in_orbits(2:7):    elements(1:6)
       ! in_orbits(8):      epoch_mjd
       ! in_orbits(9):      H
       ! in_orbits(10):     G
       ! in_orbits(j,11):   element_type_index
       ! Convert angles to radians, if needed.
       elements(1:6) = in_orbits(i,2:7)
       IF (element_types(NINT(in_orbits(i,11))) == "keplerian") THEN
          elements(3:6) = elements(3:6) * rad_deg
       ELSE IF (element_types(NINT(in_orbits(i,11))) == "delaunay") THEN
          elements(1:3) = elements(1:3) * rad_deg
       ELSE IF (element_types(NINT(in_orbits(i,11))) == "poincare") THEN
          elements(4) = elements(4) * rad_deg
       ELSE IF (element_types(NINT(in_orbits(i,11))) == "equinoctial") THEN
          elements(6) = elements(6) * rad_deg
       END IF

       ! Create a Time instance.
       CALL NEW(t, in_orbits(i,8), internal_timescale)
       IF(error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an Orbit instance.
       CALL NEW(orb, &
            elements(1:6), &
            element_types(NINT(in_orbits(i,11))), &
            "ecliptic", &
            copy(t))
       CALL NULLIFY(t)
       CALL setParameters(orb, &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
          ! Error in setParameters()
          error_code = 37
          RETURN
       END IF

       ! Initialize stochasticorbit:
       CALL NEW(storb, orb, in_covariances(i,:,:), &
            cov_type=element_types(NINT(in_orbits(i,11))), &
            element_type=element_types(NINT(in_orbits(i,11))))
       IF (error) THEN
          ! Error in StochasticOrbit instantiation
          error_code = 34
          RETURN
       END IF
       CALL NULLIFY(orb)

       ! EPHEMS:
       ! Use equatorial coordinates:
       CALL toCartesian(storb, "equatorial")

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
            cov_arr=cov_arr)
       IF (error) THEN
          ! Error in getEphemerides()
          error_code = 38
          RETURN
       END IF

       ! Now Export the ephem_arr to a flar array.
       DO j=1,SIZE(observers)
          t = getTime(observers(j))
          mjd = getMJD(t, internal_timescale)
          CALL NULLIFY(t)

          ! Input orbits correspond to one or more single-point estimates of the pdf.
          ! Make sure that the ephemeris is equatorial:
          CALL rotateToEquatorial(ephemerides(1,j))
          coordinates = getCoordinates(ephemerides(1,j))
          DO k=1,6
             stdev(k) = SQRT(cov_arr(k,k,j)) 
          END DO
          DO k=1,6
             DO l=1,6
                corr(k,l) = cov_arr(k,l,j) / &
                     (stdev(k)*stdev(l))
             END DO
          END DO

          ! Write the output ephem array.
          out_ephems(i,j,1) = coordinates(1)                ! distance
          out_ephems(i,j,2:3) = coordinates(2:3)/rad_deg    ! ra/dec
          ! FIXME: compute predicted magnitude!
          out_ephems(i,j,4) = 99.0_8                       ! mag
          out_ephems(i,j,5) = mjd                           ! ephem mjd
          ! FIXME: compute positional uncertainties.
          out_ephems(i,j,6) = 99.0_8                       ! raErr
          out_ephems(i,j,7) = 99.0_8                       ! decErr
          out_ephems(i,j,8) = 99.0_8                       ! semi-major axis
          out_ephems(i,j,9) = 99.0_8                       ! semi-minor axis
          out_ephems(i,j,10) = 99.0_8                      ! position angle

          CALL NULLIFY(ephemerides(1,j))
       END DO

       CALL NULLIFY(storb)
       DEALLOCATE(ephemerides, cov_arr)

    END DO
    DO i=1,SIZE(observers)
       CALL NULLIFY(observers(i))
    END DO
    DEALLOCATE(observers)

  END SUBROUTINE oorb_ephemeris




!!$  SUBROUTINE ephemeris(elements, mjd_utc, in_obscode, ephemeris)
!!$
!!$    IMPLICIT NONE
!!$    REAL(8), DIMENSION(:,:), INTENT(in) :: elements
!!$    REAL(8), INTENT(in) :: mjd_utc
!!$    REAL(8), INTENT(in) :: in_obscode
!!$    REAL(8), dimension(:,:) INTENT(out) :: ephemeris
!!$
!!$    TYPE (Orbit) :: orb
!!$    TYPE (Time) :: t0, t
!!$
!!$    CALL NEW(t0, elements(1,7), "TT")
!!$    CALL NEW(t, mjd_utc, "UTC")
!!$    CALL NEW(orb, elements(1,1:6), "cartesian", "ecliptic", t0)
!!$    CALL NULLIFY(t0)
!!$    CALL NULLIFY(t)
!!$    CALL NULLIFY(orb)
!!$
!!$  END SUBROUTINE ephemeris




END MODULE pyoorb
