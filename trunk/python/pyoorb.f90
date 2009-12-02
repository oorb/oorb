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
!! @version 2009-12-01
!!
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

    ! Read de405.dat
    IF (PRESENT(ephemeris_fname)) THEN
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





  SUBROUTINE oorb_ephemeris(in_norb, &
       in_orbits,                    &
       in_obscode,                   &
       in_ndate,                     &
       in_date_ephems,               &
       out_ephems,                   &
       error_code)

    ! Input/Output variables.
    INTEGER, INTENT(in)                                 :: in_norb
    ! Input flattened orbit:
    !  (track_id, elements(1:6), element_type_index, epoch, timescale, H, G)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)           :: in_orbits ! (1:norb,1:12)
    ! Observatory code as defined by the Minor Planet Center
    CHARACTER(len=*), INTENT(in)                        :: in_obscode
    INTEGER, INTENT(in)                                 :: in_ndate
    ! Ephemeris dates.
    ! (mjd, timescale)
    REAL(8), DIMENSION(in_ndate,2), INTENT(in)          :: in_date_ephems ! (1:ndate,1:2)
    ! Output ephemeris
    ! out_ephems = ((dist, ra, dec, mag, mjd, timescale), )
    REAL(8), DIMENSION(in_norb,in_ndate,6), INTENT(out) :: out_ephems ! (1:norb,1:ndate,1:6)
    ! Output error code
    INTEGER, INTENT(out)                                :: error_code

    ! Internal variables.  
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr
    TYPE (Orbit), DIMENSION(1) :: orb_arr
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: observers
    TYPE (CartesianCoordinates) :: ccoord
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides
    TYPE (Time) :: t
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    CHARACTER(len=6) :: dyn_model
    REAL(8), DIMENSION(6) :: coordinates, &
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
         phase, &
         vmag
    INTEGER :: i, &
         j
    LOGICAL, DIMENSION(10) :: perturbers

    ! Init
    errstr = ""
    error_code = 0
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
          vmag = getApparentMagnitude(H=in_orbits(i,11), &
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

          ! Write the output ephem array.
          out_ephems(i,j,1) = coordinates(1)            ! distance
          out_ephems(i,j,2) = coordinates(2)/rad_deg    ! ra
          out_ephems(i,j,3) = coordinates(3)/rad_deg    ! dec
          out_ephems(i,j,4) = vmag                      ! mag
          out_ephems(i,j,5) = mjd                       ! ephem mjd
          out_ephems(i,j,6) = NINT(in_date_ephems(j,2)) ! ephem mjd timescale

          CALL NULLIFY(ephemerides(1,j))
          CALL NULLIFY(orb_lt_corr_arr(1,j))

       END DO

       DEALLOCATE(ephemerides, orb_lt_corr_arr)

    END DO
    DO i=1,SIZE(observers)
       CALL NULLIFY(observers(i))
    END DO
    DEALLOCATE(observers)

  END SUBROUTINE oorb_ephemeris




  SUBROUTINE oorb_ephemeris_covariance(in_norb, &
       in_orbits,                    &
       in_covariances,               &
       in_obscode,                   &
       in_ndate,                     &
       in_date_ephems,               &
       out_ephems,                   &
       error_code)

    ! Input/Output variables.
    INTEGER, INTENT(in)                                   :: in_norb
    ! Input flattened orbit:
    !  (track_id, elements(1:6), element_type_index, epoch, timescale, H, G)
    REAL(8),DIMENSION(in_norb,12), INTENT(in)             :: in_orbits ! (1:norb,1:12)
    ! Covariance matrices:
    REAL(8), DIMENSION(in_norb,6,6), INTENT(in)           :: in_covariances ! (1:norb,1:6,1:6)
    ! Observatory code as defined by the Minor Planet Center
    CHARACTER(len=*), INTENT(in)                          :: in_obscode
    INTEGER, INTENT(in)                                   :: in_ndate
    ! Ephemeris dates.
    ! (mjd, timescale)
    REAL(8), DIMENSION(in_ndate,2), INTENT(in)            :: in_date_ephems ! (1:ndate,1:2)
    ! Output ephemeris
    ! out_ephems = ((dist, ra, dec, mag, mjd, timescale, raErr, decErr, smaa, smia, pa), )
    REAL(8), DIMENSION(in_norb,in_ndate,11), INTENT(out)  :: out_ephems ! (1:norb,1:ndate,1:11)
    ! Output error code
    INTEGER, INTENT(out)                                  :: error_code

    ! Internal variables.  
    TYPE (StochasticOrbit) :: storb
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr
    TYPE (Orbit), DIMENSION(1) :: orb_arr
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: observers
    TYPE (CartesianCoordinates) :: ccoord
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides
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
          vmag = getApparentMagnitude(H=in_orbits(i,11), &
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
          out_ephems(i,j,1) = coordinates(1)            ! distance
          out_ephems(i,j,2) = coordinates(2)/rad_deg    ! ra
          out_ephems(i,j,3) = coordinates(3)/rad_deg    ! dec
          out_ephems(i,j,4) = vmag                      ! mag
          out_ephems(i,j,5) = mjd                       ! ephem mjd
          out_ephems(i,j,6) = NINT(in_date_ephems(j,2)) ! ephem mjd timescale
          out_ephems(i,j,7) = sigma_ra                  ! raErr
          out_ephems(i,j,8) = sigma_dec                 ! decErr
          out_ephems(i,j,9) = smaa                      ! semi-major axis
          out_ephems(i,j,10) = smia                     ! semi-minor axis
          out_ephems(i,j,11) = pa                       ! position angle

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
