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
  CHARACTER(len=11), DIMENSION(6), PUBLIC :: ORBITAL_ELEMENTS = (/             &
       "cartesian  ",&
       "cometary   ",&
       "keplerian  ",&
       "delaunay   ",&
       "poincare   ",&
       "equinoctial" &
       /)
  CHARACTER(len=8), PARAMETER :: internal_timescale = "TAI"


CONTAINS

  SUBROUTINE init(ephemeris_fname, error_verbosity, info_verbosity)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in), OPTIONAL :: ephemeris_fname
    INTEGER, INTENT(in), OPTIONAL :: error_verbosity
    INTEGER, INTENT(in), OPTIONAL :: info_verbosity

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

  END SUBROUTINE init



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



  SUBROUTINE oorb_ephemeris(in_orbit,                                         &
       in_covariance,                                    &
       in_obscode,                                       &
       in_num_ephems,                                    &
       in_step,                                          &
       out_ephems,                                       &
       error_code)

    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    REAL(8),DIMENSION(11), INTENT(in)                   :: in_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), INTENT(in)                 :: in_covariance
    CHARACTER(len=4), INTENT(in)                        :: in_obscode
    ! Compute epehemeris from the orbit epoch to epoch+in_step*in_num_ephems
    ! Number of ephemeris to compute.
    INTEGER, INTENT(in)                                 :: in_num_ephems
    ! Ephemeris step in fractional days.
    REAL(8), INTENT(in)                                 :: in_step
    ! Output ephemeris
    ! out_ephems = ((dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa), )
    REAL(8), DIMENSION(in_num_ephems,10), INTENT(out)   :: out_ephems
    ! Output error code
    INTEGER, INTENT(out)                                :: error_code

    ! Internal variables.
    TYPE (StochasticOrbit)                              :: storb
    TYPE (Orbit), DIMENSION(:), POINTER         :: orb_arr_in
    TYPE (Observatories) :: obsies
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: observers
    TYPE (SphericalCoordinates), DIMENSION(:,:),POINTER :: ephemerides
    TYPE (Time)                                         :: t
    CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: id_arr_in
    CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: element_type_pdf_arr_in
    CHARACTER(len=DESIGNATION_LEN)              :: id
    CHARACTER(len=INTEGRATOR_LEN)              :: integrator
    CHARACTER(len=6)                                    :: dyn_model
    CHARACTER(len=11)                           :: element_type
    REAL(8), DIMENSION(:,:,:), POINTER         :: cov_arr_in
    REAL(8), DIMENSION(:,:,:), POINTER                  :: cov_arr
    REAL(8), DIMENSION(:,:), POINTER                    :: pdfs_arr
    REAL(8), DIMENSION(:,:), POINTER           :: HG_arr_in
    REAL(8), DIMENSION(:,:), ALLOCATABLE       :: jac_arr_in
    REAL(8), DIMENSION(:), ALLOCATABLE         :: rchi2_arr_in, &
         pdf_arr_in, &
         reg_apr_arr_in
    REAL(8), DIMENSION(1,6,6)                            :: cov
    REAL(8), DIMENSION(6,6)                            :: corr
    REAL(8), DIMENSION(6)                      :: elements, coordinates, stdev
    REAL(8)                                            :: integration_step, mjd, step, timespan 
    INTEGER                                             :: element_type_index, i, j, k, l, norb, nstep
    LOGICAL, DIMENSION(:), POINTER                      :: perturbers


    CALL NEW(obsies)


    ! Init
    error_code = 0
    dyn_model = "n-body"
    integration_step = 5.0_8
    ALLOCATE(perturbers(10), stat=error_code)
    IF(error_code /= 0) THEN
       ! Error in allocating memory!
       error_code = 57
       RETURN
    END IF
    perturbers = .TRUE.


    ! ORBITS:
    ! Get the element type from the input flattened orbit.
    element_type_index = in_orbit(11)
    IF(element_type_index .LE. 0 .OR.                                       &
         element_type_index .GT. SIZE(ORBITAL_ELEMENTS)) THEN
       ! Error: unsupported orbital elements.
       error_code = 58
       RETURN
    END IF
    element_type = ORBITAL_ELEMENTS(element_type_index)

    ! Compute timespan given step and num_ephems.
    timespan = in_step * in_num_ephems
    step = in_step          ! Might be modified by the code below.

    ! Init vars.
    error = .FALSE.
    error_code = 0
    ! We try to keep these two subroutines similar so that we can merge in 
    ! the future, hence the norb = 1 thing here.
    norb = 1

    ! Get the element type from the input flattened in_orbit.
    element_type_index = in_orbit(11)
    IF(element_type_index .LE. 0 .OR.                               &
         element_type_index .GT. SIZE(ORBITAL_ELEMENTS)) THEN
       ! Error: unsupported orbital elements.
       error_code = 58
       RETURN
    END IF
    element_type = ORBITAL_ELEMENTS(element_type_index)

    ! Allocate memory for the orbit arrays.
    ALLOCATE(id_arr_in(norb), orb_arr_in(norb), &
         element_type_pdf_arr_in(norb), cov_arr_in(norb,6,6), &
         HG_arr_in(norb,2), pdf_arr_in(norb), &
         rchi2_arr_in(norb), jac_arr_in(norb,3), &
         reg_apr_arr_in(norb), stat=error_code)
    IF (error_code /= 0) THEN
       ! Error in memory allocation
       error_code = 33
       RETURN
    END IF
    jac_arr_in = -1.0_8
    pdf_arr_in = -1.0_8
    cov_arr_in = 0.0_8
    id_arr_in = " "
    ! There are no H/G in the input in_orbit.
    FORALL (j=1:norb)
       HG_arr_in(j,1:2) = (/ 99.0_8, 9.9_8 /)
    END FORALL
    element_type_pdf_arr_in = element_type

    ! Get each flattened orbit and create an Orbit instance.
    DO j=1,norb
       ! Just to beat on a dead horse:
       ! in_orbit(1):      id
       ! in_orbit(2:7):    elements(1:6)
       ! in_orbit(8):      epoch_mjd
       ! in_orbit(9):      H
       ! in_orbit(10):     G
       ! in_orbit(j,11):   element_type_index
       ! Convert angles to radians, if needed.
       WRITE(id_arr_in(j), fmt="(A,I10.10)") "TRK", IDINT(in_orbit(1))

       elements(1:6) = in_orbit(2:7)
       IF (element_type == "keplerian") THEN
          elements(3:6) = elements(3:6) * rad_deg
       ELSE IF (element_type == "delaunay") THEN
          elements(1:3) = elements(1:3) * rad_deg
       ELSE IF (element_type == "poincare") THEN
          elements(4) = elements(4) * rad_deg
       ELSE IF (element_type == "equinoctial") THEN
          elements(6) = elements(6) * rad_deg
       END IF

       ! Create a Time instance.
       CALL NEW(t, in_orbit(8), internal_timescale)
       IF(error) THEN
          ! Error in creating a Time instance.
          error_code = 57
          RETURN
       END IF

       ! Now create an Orbit instance.
       CALL NEW(orb_arr_in(j),                                 &
            elements(1:6),                                 &
            element_type,                                  &
            "ecliptic",                                    &
            copy(t))
       CALL NULLIFY(t)

       ! Covariance. We do not need to re-compute it since we already have 
       ! it.
       cov_arr_in(j,:,:) = in_covariance(:,:)
    END DO

    ! Initialize stochasticorbits if uncertainty information available:
    IF (norb > 1 .AND. &
         ALL(pdf_arr_in > 0.0_8) .AND. &
         ALL(jac_arr_in > 0.0_8)) THEN
       CALL NULLIFY(storb)
       CALL NEW(storb, orb_arr_in, pdf_arr_in, &
            element_type_pdf_arr_in(1), jac_arr=jac_arr_in, &
            reg_apr_arr=reg_apr_arr_in, &
            rchi2_arr=rchi2_arr_in)
       id = id_arr_in(1)
       DEALLOCATE(id_arr_in)
       ALLOCATE(id_arr_in(1))
       id_arr_in(1) = id
    ELSE IF (norb == 1 .AND. ALL(cov_arr_in(:,1,1) > 0.0_8)) THEN
       CALL NEW(storb, orb_arr_in(1), cov_arr_in(1,:,:), &
            cov_type=element_type, element_type=element_type)
    ELSE
       ! Error in StochasticOrbit instantiation
       error_code = 34
       RETURN
    END IF

    ! Cleanup.
    DO j=1,norb
       CALL NULLIFY(orb_arr_in(j))
    END DO
    DEALLOCATE(orb_arr_in, pdf_arr_in, rchi2_arr_in, jac_arr_in, &
         reg_apr_arr_in, element_type_pdf_arr_in)

    ! EPHEMS:
    ! Init vars
    error = .FALSE.
    error_code = 0
    IF (step <= 0.0_8) THEN
       nstep = 1        
    ELSE
       step = SIGN(ABS(step),timespan)
       IF (ABS(timespan) > 10.0_8*EPSILON(timespan) .AND. &
            ABS(timespan) < ABS(step)) THEN
          step = timespan
       END IF
       nstep = NINT(timespan/step) + 1
    END IF
    integration_step = MIN(ABS(step),integration_step)

    ! Use equatorial coordinates:
    CALL toCartesian(storb, "equatorial")

    t = getTime(storb)
    mjd = getMJD(t, internal_timescale)
    CALL NULLIFY(t)
    ALLOCATE(observers(nstep))
    DO j=1,nstep
       CALL NEW(t, mjd+(j-1)*step, internal_timescale)
       IF (error) THEN
          ! Error in creating a new Time object.
          error_code = 35
          RETURN
       END IF
       ! Compute heliocentric observatory coordinates
       observers(j) = getObservatoryCCoord(obsies, in_obscode, t)
       IF (error) THEN
          ! Error in getObservatoryCCoord()
          error_code = 36
          RETURN
       END IF
       CALL rotateToEquatorial(observers(j))
       CALL NULLIFY(t)
    END DO

    ! Set integration parameters
    IF (dyn_model == "n-body") THEN
       CALL setParameters(storb, &
            dyn_model=dyn_model, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
    ELSE
       CALL setParameters(storb, &
            dyn_model=dyn_model)
    END IF
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
         pdfs_arr=pdfs_arr)
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
       IF (containsSampledPDF(storb)) THEN
          ! We do not support exporting ephems from ranging orbits yet.
          ! FIXME: support ephems from ranging orbits?
          error_code = 59
          RETURN
       END IF

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
       out_ephems(j, 1) = coordinates(1)                ! distance
       out_ephems(j, 2:3) = coordinates(2:3)/rad_deg    ! ra/dec
       ! FIXME: compute predicted magnitude!
       out_ephems(j, 4) = 99.0_8                       ! mag
       out_ephems(j, 5) = mjd                           ! ephem mjd
       ! FIXME: compute positional uncertainties.
       out_ephems(j, 6) = 99.0_8                       ! raErr
       out_ephems(j, 7) = 99.0_8                       ! decErr
       out_ephems(j, 8) = 99.0_8                       ! semi-major axis
       out_ephems(j, 9) = 99.0_8                       ! semi-minor axis
       out_ephems(j, 10) = 99.0_8                      ! position angle
    END DO
    IF(error_code /= 0) THEN
       RETURN
    END IF

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
