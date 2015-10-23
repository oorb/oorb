!====================================================================!
!                                                                    !
! Copyright 2002-2014,2015                                           !
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
!! *Class*description*:
!!
!! Type and routines for observatories.
!! 
!! @author  MG, JV
!! @version 2015-10-23
!!
MODULE Observatories_cl

  USE Base_cl
  USE File_cl
  USE Observatory_cl
  USE Time_cl
  USE CartesianCoordinates_cl

  USE utilities
  USE planetary_data
  USE linal

  IMPLICIT NONE
  PRIVATE :: new_Obsies
  PRIVATE :: new_Obsies_file
  PRIVATE :: nullify_Obsies
  PRIVATE :: copy_Obsies
  PRIVATE :: exist_Obsies
  PRIVATE :: getPosition_Obsies
  PRIVATE :: getName_Obsies
  PRIVATE :: getObservatoryCCoord_code
  PRIVATE :: getObservatoryCCoord_obsy

  TYPE Observatories
     PRIVATE
     TYPE (Observatory), DIMENSION(:), POINTER :: observatory_arr => NULL()
     CHARACTER(len=FNAME_LEN)                  :: code_fname          =  ""
     INTEGER                                   :: no_of_observatories =  0
     LOGICAL                                   :: is_initialized      = .FALSE.
  END TYPE Observatories

  INTERFACE NEW
     MODULE PROCEDURE new_Obsies
     MODULE PROCEDURE new_Obsies_file
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_Obsies
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_Obsies
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_Obsies
  END INTERFACE exist

  INTERFACE getPosition
     MODULE PROCEDURE getPosition_Obsies
  END INTERFACE getPosition

  INTERFACE getName
     MODULE PROCEDURE getName_Obsies
  END INTERFACE getName

  INTERFACE getObservatory
     MODULE PROCEDURE getObservatory_Obsies
  END INTERFACE getObservatory

  INTERFACE getObservatoryCCoord
     MODULE PROCEDURE getObservatoryCCoord_code
     MODULE PROCEDURE getObservatoryCCoord_obsy
  END INTERFACE getObservatoryCCoord


CONTAINS





  !! *Description*:
  !!
  !! Default initialization routine. Reads the observatory data from
  !! the file to which the parameter CODE_FNAME points.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Obsies(this)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(inout) :: this

    CALL NEW(this, TRIM(OORB_DATA_DIR) // "/" &
         // TRIM(CODE_FNAME))

  END SUBROUTINE new_Obsies





  !! *Description*:
  !!
  !! Reads the observatory data from the given file.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Obsies_file(this, filename)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(inout) :: this
    CHARACTER(len=*), INTENT(in)        :: filename
    TYPE (File)                         :: code_file
    CHARACTER(len=96)                   :: name, form
    CHARACTER(len=OBSY_CODE_LEN)        :: code
    REAL(bp), DIMENSION(3)              :: coordinates, position
    INTEGER                             :: i, err, length

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    length = LEN_TRIM(filename)
    IF (length /= 0 .AND. length <= FNAME_LEN) THEN
       this%code_fname = TRIM(filename)
    ELSE IF (length > FNAME_LEN) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / new", &
            "Filename too long; adjust parameters.", 1)
       RETURN
    ELSE
       this%code_fname = TRIM(CODE_FNAME)
    END IF

    CALL NEW(code_file, TRIM(this%code_fname))
    IF (error) THEN
       CALL errorMessage("Observatories / new", &
            "TRACE BACK 1", 1)
       RETURN
    END IF

    CALL setActionRead(code_file)
    CALL setStatusOld(code_file)
    CALL OPEN(code_file)
    IF (error) THEN
       CALL errorMessage("Observatories / new", &
            "TRACE BACK 4", 1)
       RETURN
    END IF

    this%no_of_observatories = getNrOfLines(code_file)

    ALLOCATE(this%observatory_arr(this%no_of_observatories), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / new", &
            "Could not allocate memory for observatories.", 1)
       RETURN
    END IF

    form = "(A3,1X,F9.5,F8.6,F9.6,A96)"
    DO i=1, this%no_of_observatories

       READ(getUnit(code_file), TRIM(form), iostat=err) code, coordinates, name
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observatories / new", &
               "Error while reading from observatory datafile.", 1)
          RETURN
       END IF

       ! Make the transformation to Cartesian geocentric equatorial
       ! coordinates in AUs (body-fixed coordinates):
       coordinates(1) = coordinates(1)*rad_deg
       position(1)    = r_earth*coordinates(2)*COS(coordinates(1))
       position(2)    = r_earth*coordinates(2)*SIN(coordinates(1))
       position(3)    = r_earth*coordinates(3)

       CALL NEW(this%observatory_arr(i), code, name, position)
       IF (error) THEN
          CALL errorMessage("Observatories / new", &
               "TRACE BACK 6", 1)
          RETURN
       END IF

    END DO

    CALL NULLIFY(code_file)
    IF (error) THEN
       CALL errorMessage("Observatories / new", &
            "TRACE BACK 7", 1)
       RETURN
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE new_Obsies_file





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_Obsies(this)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(inout) :: this
    INTEGER :: i, err

    IF (ASSOCIATED(this%observatory_arr)) THEN
       DO i=1,SIZE(this%observatory_arr,dim=1)
          CALL NULLIFY(this%observatory_arr(i))
       END DO
       DEALLOCATE(this%observatory_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observatories / nullify", &
               "Could not deallocate memory.", 1)
          RETURN
       END IF
    END IF
    this%no_of_observatories =  0
    this%code_fname          =  ""
    this%is_initialized      =  .FALSE.

  END SUBROUTINE nullify_Obsies





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  !! Returns error.
  !!
  FUNCTION copy_Obsies(this)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this
    TYPE (Observatories)             :: copy_Obsies
    INTEGER :: err, i, nobsies

    nobsies = 0
    IF (ASSOCIATED(this%observatory_arr)) THEN
       nobsies = SIZE(this%observatory_arr,dim=1)
       ALLOCATE(copy_Obsies%observatory_arr(nobsies), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observatories / copy", &
               "Could not allocate memory.", 1)
          RETURN
       END IF
       DO i=1,nobsies
          copy_Obsies%observatory_arr(i) = &
               copy(this%observatory_arr(i))
       END DO
    END IF
    copy_Obsies%no_of_observatories = nobsies
    copy_Obsies%observatory_arr     = this%observatory_arr
    copy_Obsies%code_fname          = this%code_fname
    copy_Obsies%is_initialized      = this%is_initialized

  END FUNCTION copy_Obsies





  !! *Description*:
  !!
  !! Returns the status of this object, i.e. whether
  !! it exists or not.
  !!
  LOGICAL FUNCTION exist_Obsies(this)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this

    exist_Obsies = this%is_initialized

  END FUNCTION exist_Obsies





  !! *Description*:
  !!
  !! Returns the equation of the equinoxes (difference between
  !! apparent sidereal time and mean sidereal time) in radians.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION equationOfEquinoxes(t)

    IMPLICIT NONE
    TYPE(Time), INTENT(inout) :: t
    REAL(bp)                  :: mjd_tt, oblm, dpsi, deps

    IF (.NOT. exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / equationOfEquinoxes", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd_tt = getMJD(t, "TT")
    oblm = meanObliquity(mjd_tt)
    CALL getNutationAngles(t, dpsi, deps)
    IF (error) THEN
       CALL errorMessage("Observatories / equationOfEquinoxes", &
            "TRACE BACK", 1)
       RETURN
    END IF

    equationOfEquinoxes = rad_asec*dpsi*COS(oblm)

  END FUNCTION equationOfEquinoxes






  !! *Description*:
  !!
  !! Returns name of observatory based on the given observatory code.
  !!
  !! Returns error.
  !!
  CHARACTER(len=128) FUNCTION getName_Obsies(this, code)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this
    CHARACTER(len=*), INTENT(in)     :: code
    CHARACTER(len=OBSY_CODE_LEN)     :: trial_code
    INTEGER                          :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getName", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    DO i=1, this%no_of_observatories
       trial_code = getCode(this%observatory_arr(i))
       IF (TRIM(code) == TRIM(trial_code)) THEN
          getName_Obsies = getName(this%observatory_arr(i))
          RETURN
       END IF
    END DO

    ! The matching code could not be found, which indicates it is not
    ! included in the observatory code file:
    error = .TRUE.
    CALL errorMessage("Observatories / getName", &
         "Code " // TRIM(code) // " does not refer to any observatory listed in " &
         // this%code_fname, 1)

  END FUNCTION getName_Obsies





  !! *Description*:
  !!
  !! Returns observatory object based on the given observatory code.
  !!
  !! Returns error.
  !!
  FUNCTION getObservatory_Obsies(this, code)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this
    CHARACTER(len=*), INTENT(in)     :: code
    TYPE (Observatory)               :: getObservatory_Obsies
    CHARACTER(len=OBSY_CODE_LEN)     :: trial_code
    REAL(bp), DIMENSION(6)           :: position
    INTEGER                          :: i, indx

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getObservatory", &
            "Object has not been initialized yet.", 1)
       RETURN
    END IF

    indx = INDEX(code,"-")
    IF (indx == 0) THEN
       DO i=1, this%no_of_observatories
          trial_code = getCode(this%observatory_arr(i))
          IF (error) THEN
             CALL errorMessage("Observatories / getObservatory", &
                  "TRACE BACK (5)", 1)
             RETURN
          END IF
          ! If the matching code is found, then 
          ! collect the data and return:
          IF (TRIM(code) == TRIM(trial_code)) THEN
             getObservatory_Obsies = copy(this%observatory_arr(i))
             IF (error) THEN
                CALL errorMessage("Observatories / getObservatory", &
                     "TRACE BACK (10)", 1)
                CALL NULLIFY(getObservatory_Obsies)
                RETURN
             END IF
             RETURN
          END IF
       END DO
    ELSE
       CALL toInt(code(indx+1:), i, error)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatory", &
               "TRACE BACK (15)", 1)
          CALL NULLIFY(getObservatory_Obsies)
          RETURN
       END IF
       position = 0.0_bp
       CALL NEW(getObservatory_Obsies, code, planetary_locations(i), position)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatory", &
               "TRACE BACK (20)", 1)
          CALL NULLIFY(getObservatory_Obsies)
          RETURN
       END IF
       RETURN
    END IF

    ! The matching code could not be found, which indicates it is not
    ! included in the observatory code file:
    error = .TRUE.
    CALL errorMessage("Observatories / getObservatory", &
         "Code " // TRIM(code) // &
         " does not refer to any observatory listed in " &
         // this%code_fname, 1)

  END FUNCTION getObservatory_Obsies





  !! *Description*:
  !!
  !! Returns heliocentric equatorial Cartesian coordinates for the
  !! given observatory code and epoch.
  !!
  !! Returns error
  !!
  FUNCTION getObservatoryCCoord_code(this, code, t)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in)  :: this
    CHARACTER(len=*), INTENT(in)      :: code
    TYPE (Time), INTENT(in)           :: t
    TYPE (CartesianCoordinates)       :: getObservatoryCCoord_code

    TYPE (Time)                       :: t_
    TYPE (CartesianCoordinates)       :: geocenter_ccoord, geocentric_obs_ccoord
    REAL(bp), DIMENSION(:,:), POINTER :: coordinates => NULL()
    REAL(bp)                          :: mjd_tt
    INTEGER                           :: err, indx, i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "t has not yet been initialized.", 1)
       RETURN
    END IF

    t_ = copy(t)
    indx = INDEX(code,"-")
    IF (indx == 0) THEN
       ! We are dealing with an Earth based observer or a satellite
       ! with its coordinates given relative to the Earth.

       ! Equatorial Cartesian coordinates for the Earth (geocenter):
       mjd_tt = getMJD(t_, "TT")
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       coordinates => JPL_ephemeris(mjd_tt, 3, 11, error)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "Could not get planetary ephemeris (5).", 1)
          RETURN
       END IF
       CALL NEW(geocenter_ccoord, coordinates(1,:), "equatorial", copy(t_))
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       DEALLOCATE(coordinates, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
       ! Equatorial geocentric Cartesian state for the observatory
       geocentric_obs_ccoord = getGeocentricObservatoryCCoord(this, code, t_)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       ! Equatorial heliocentric Cartesian state for the observatory
       getObservatoryCCoord_code = copy(geocenter_ccoord + geocentric_obs_ccoord)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
       CALL rotateToEquatorial(getObservatoryCCoord_code)
       CALL NULLIFY(t_)
       CALL NULLIFY(geocenter_ccoord)
       CALL NULLIFY(geocentric_obs_ccoord)

    ELSE
       ! We are dealing with a planetocentric observer.

       ! Equatorial Cartesian planetocentric coordinates:
       mjd_tt = getMJD(t_, "TT")
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
       CALL toInt(TRIM(code(indx+1:)), i, error)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "Could not convert string to integer.", 1)
          RETURN
       END IF
       IF (i < 1 .OR. i>11) THEN
          error = .TRUE.
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "Observatory code invalid: " // TRIM(code(indx+1:)), 1)
          RETURN
       END IF
       coordinates => JPL_ephemeris(mjd_tt, i, 11, error)
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "Could not get planetary ephemeris (10).", 1)
          RETURN
       END IF
       CALL NEW(getObservatoryCCoord_code, coordinates(1,:), "equatorial", copy(t_))
       IF (error) THEN
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "TRACE BACK (30)", 1)
          RETURN
       END IF
       DEALLOCATE(coordinates, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observatories / getObservatoryCCoord", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF
       CALL rotateToEquatorial(getObservatoryCCoord_code)
       CALL NULLIFY(t_)

    END IF

  END FUNCTION getObservatoryCCoord_code





  !! *Description*:
  !!
  !! Returns heliocentric equatorial Cartesian coordinates for the
  !! given observatory and epoch.
  !!
  !! Returns error
  !!
  FUNCTION getObservatoryCCoord_obsy(this, obsy, t)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this
    TYPE (Observatory), INTENT(in)   :: obsy
    TYPE (Time), INTENT(in)          :: t
    TYPE (CartesianCoordinates)      :: getObservatoryCCoord_obsy
    CHARACTER(len=OBSY_CODE_LEN)     :: code

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obsy)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "obsy has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "t has not yet been initialized.", 1)
       RETURN
    END IF

    code = getCode(obsy)
    IF (error) THEN
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    getObservatoryCCoord_obsy = getObservatoryCCoord(this, code, t)
    IF (error) THEN
       CALL errorMessage("Observatories / getObservatoryCCoord", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

  END FUNCTION getObservatoryCCoord_obsy





  !! *Description*:
  !!
  !! Returns geocentric equatorial Cartesian coordinates for the given
  !! observatory code and epoch.
  !!
  !! Returns error.
  !!
  FUNCTION getGeocentricObservatoryCCoord(this, code, t)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this
    CHARACTER(len=*), INTENT(in)     :: code
    TYPE (Time), INTENT(inout)       :: t

    TYPE (CartesianCoordinates)      :: getGeocentricObservatoryCCoord
    TYPE (Time)                      :: t_obs
    ! Diurnal rotation matrix (transformation from body-fixed to
    ! true-of-date frames):
    REAL(bp), DIMENSION(3,3)         :: diurnal_matrix
    ! Rotation matrix from true-of-date frame to the reference
    ! system in which observations are given (i.e., for
    ! nutation and precession):
    REAL(bp), DIMENSION(3,3)         :: pn_matrix
    REAL(bp), DIMENSION(3)           :: omega, bf_position, &
         true_position, true_velocity, position, velocity
    REAL(bp)                         :: epoch_obs, gast
    CHARACTER(len=20)                :: frame_obs, refsys_obs, &
         frame_true, refsys_true

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "t has not yet been initialized.", 1)
       RETURN
    END IF

    ! Reference system in which the observations are given
    frame_obs = "equatorial"
    refsys_obs = "mean"
    epoch_obs = 51544.5_bp ! J2000.0
    CALL NEW(t_obs, epoch_obs, "TT")
    IF (error) THEN
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    ! True-of-date reference system
    frame_true = "equatorial"
    refsys_true = "true-of-date"

    ! Earth angular velocity (rad/d)
    omega(1:2) = 0.0_bp
    omega(3) = two_pi*1.00273790934_bp

    ! Greenwich Apparent Sidereal Time = Greenwich Mean Sidereal Time +
    ! Equation of the Equinoxes
    gast = getGMST(t) + equationOfEquinoxes(t)
    IF (error) THEN
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

    ! Diurnal rotation matrix (transformation from body-fixed to
    ! true-of-date frames), neglecting polar motion
    diurnal_matrix = rotationMatrix(-gast,3)

    ! Coordinates of the observatory in the body-fixed frame
    bf_position = getPosition(this, code)
    IF (error) THEN
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF

    ! Coordinates of the observatory in the true-of-date frame
    true_position = MATMUL(diurnal_matrix, bf_position)

    ! Observatory velocity in the true-of-date frame: V = OMEGA X R
    true_velocity = cross_product(omega, true_position)

    ! Rotation matrix from the true-of-date reference system
    ! to the reference system of the observations. Notice that
    ! the EPOCH of the true-of-date frame is the time at which
    ! the observation is obtained: therefore this statement
    ! must be within the loop on observations
    pn_matrix = precessionAndNutationMatrix(frame_true, refsys_true, &
         t, frame_obs, refsys_obs, t_obs)
    IF (error) THEN
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF

    ! Coordinates of the observatory in the frame of observations
    position = MATMUL(pn_matrix, true_position)
    velocity = MATMUL(pn_matrix, true_velocity)

    CALL NEW(getGeocentricObservatoryCCoord, position, velocity, &
         frame_true, t)
    IF (error) THEN
       CALL errorMessage("Observatories / getGeocentricObservatoryCCoord", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF

    CALL NULLIFY(t_obs)

  END FUNCTION getGeocentricObservatoryCCoord





  !! *Description*:
  !!
  !! Returns the Cartesian geocentric equatorial coordinates in AUs
  !! (body-fixed coordinates).
  !!
  !! Returns error.
  !!
  FUNCTION getPosition_Obsies(this, code)

    IMPLICIT NONE
    TYPE (Observatories), INTENT(in) :: this
    CHARACTER(len=*), INTENT(in)     :: code
    REAL(bp), DIMENSION(3)           :: getPosition_Obsies
    CHARACTER(len=OBSY_CODE_LEN)     :: trial_code
    INTEGER                          :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getPosition", &
            "Object has not been initialized yet.", 1)
       RETURN
    END IF

    DO i=1, this%no_of_observatories
       trial_code = getCode(this%observatory_arr(i))
       ! If the matching code is found, then collect the data and return:
       IF (TRIM(code) == TRIM(trial_code)) THEN
          getPosition_Obsies = getPosition(this%observatory_arr(i))
          RETURN
       END IF
    END DO

    ! The matching code could not be found, which indicates it is not included in the
    ! observatory code file:
    error = .TRUE.
    CALL errorMessage("Observatories / getPosition", &
         "Code " // TRIM(code) // " does not refer to any observatory listed in " &
         // this%code_fname, 1)

  END FUNCTION getPosition_Obsies





  !! *Description*:
  !!
  !! Computation of the nutation matrix from nutation angles 
  !!
  !! Uses all 106 terms of the series for Nutation in Longitude and Obliquity
  !! as listed in the Expl. Suppl. to AA, pp. 112-113.
  !!
  !! Note 1. No corrections (Herring 1987) as listed in p. 116 (presumable due to the
  !! difference of the real Earth from the Wahr model). Effect?
  !! Note 2. No planetary terms (Vondrak 1983), p. 118-119 with 85 terms. Should be used 
  !! if mas accuracy is required.
  !!
  !! Based on f77 routine by Mario Carpino (nutn80.f).
  !!
  !! Returns error.
  !!
  SUBROUTINE getNutationAngles(t0, dpsi, deps)

    IMPLICIT NONE
    TYPE(Time), INTENT(inout) :: t0
    REAL(bp),INTENT(out)      :: dpsi, deps
    REAL(bp) :: mjd_tt, dl, dp, df, dd, dn, t1
    REAL(lp) :: l, n, cp, sp, cx, sx, cd, sd, cn, sn, cl, &
         sl, cp2, sp2, cd2, sd2, cn2, sn2, cl2, sl2, ca, &
         sa, cb, sb, cc, sc, cv, sv, ce, se, cf, sf, cg, &
         sg, ch, sh, cj, sj, ck, sk, cm, sm, cq, sq, cr, &
         sr, cs, ss, ct, st, cu, su, cw, sw, t,  &
         t2, t3, p, x, d

    IF (.NOT. exist(t0)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getNutationAngles", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd_tt = getMJD(t0, "TT")
    IF (error) THEN
       CALL errorMessage("Observatories / getNutationAngles", &
            "TRACE BACK", 1)
       RETURN
    END IF
    t1 = (mjd_tt-51544.5_bp)/36525.0_bp
    t = REAL(t1,kind=lp)
    t2 = t*t
    t3 = t2*t

    ! Fundamental arguments (IAU 1980, Expl.Suppl. to AA p. 114):
    dl = ( 485866.733_bp +1717915922.633_bp*t1 +31.310*t2 +0.064*t3)*rad_asec
    dp = (1287099.804_bp + 129596581.224_bp*t1 - 0.577*t2 -0.012*t3)*rad_asec
    df = ( 335778.877_bp +1739527263.137_bp*t1 -13.257*t2 +0.011*t3)*rad_asec
    dd = (1072261.307_bp +1602961601.328_bp*t1 - 6.891*t2 +0.019*t3)*rad_asec
    dn = ( 450160.280_bp -   6962890.539_bp*t1 + 7.455*t2 +0.008*t3)*rad_asec
    l = REAL(MOD(dl,two_pi),kind=lp)
    p = REAL(MOD(dp,two_pi),kind=lp)
    x = REAL(MOD(df,two_pi)*2.0_bp,kind=lp)
    d = REAL(MOD(dd,two_pi),kind=lp)
    n = REAL(MOD(dn,two_pi),kind=lp)
    cl = COS(l)
    sl = SIN(l)
    cp = COS(p)
    sp = SIN(p)
    cx = COS(x)
    sx = SIN(x)
    cd = COS(d)
    sd = SIN(d)
    cn = COS(n)
    sn = SIN(n)
    cp2 = 2.0_lp*cp*cp - 1.0_lp
    sp2 = 2.0_lp*sp*cp
    cd2 = 2.0_lp*cd*cd - 1.0_lp
    sd2 = 2.0_lp*sd*cd
    cn2 = 2.0_lp*cn*cn - 1.0_lp
    sn2 = 2.0_lp*sn*cn
    cl2 = 2.0_lp*cl*cl - 1.0_lp
    sl2 = 2.0_lp*sl*cl
    ca = cx*cd2 +sx*sd2
    sa = sx*cd2 -cx*sd2
    cb = ca*cn -sa*sn
    sb = sa*cn +ca*sn
    cc = cb*cn -sb*sn
    sc = sb*cn +cb*sn
    cv = cx*cd2 -sx*sd2
    sv = sx*cd2 +cx*sd2
    ce = cv*cn -sv*sn
    se = sv*cn +cv*sn
    cf = ce*cn -se*sn
    sf = se*cn +ce*sn
    cg = cl*cd2 +sl*sd2
    sg = sl*cd2 -cl*sd2
    ch = cx*cn2 -sx*sn2
    sh = sx*cn2 +cx*sn2
    cj = ch*cl -sh*sl
    sj = sh*cl +ch*sl
    ck = cj*cl -sj*sl
    sk = sj*cl +cj*sl
    cm = cx*cl2 +sx*sl2
    sm = sx*cl2 -cx*sl2
    cq = cl*cd +sl*sd
    sq = sl*cd -cl*sd
    cr = 2.0_lp*cq*cq - 1.0_lp
    sr = 2.0_lp*sq*cq
    cs = cx*cn -sx*sn
    ss = sx*cn +cx*sn
    ct = cs*cl -ss*sl
    st = ss*cl +cs*sl
    cu = cf*cl +sf*sl
    su = sf*cl -cf*sl
    cw = cp*cg -sp*sg
    sw = sp*cg +cp*sg

    ! Series for dpsi
    ! Terms 1-54:
    dpsi = -(171996.0_lp+174.2_lp*t)*sn      +(2062.0_lp+0.2_lp*t)*sn2 &
         +46.0_lp*(sm*cn+cm*sn)              -11.0_lp*sm &
         -3.0_lp*(sm*cn2+cm*sn2)             -3.0_lp*(sq*cp-cq*sp) &
         -2.0_lp*(sb*cp2-cb*sp2)             +(sn*cm-cn*sm) &
         -(13187.0_lp+1.6_lp*t)*sc           +(1426.0_lp-3.4_lp*t)*sp &
         -(517.0_lp-1.2_lp*t)*(sc*cp+cc*sp)  +(217.0_lp-0.5_lp*t)*(sc*cp-cc*sp) &
         +(129.0_lp+0.1_lp*t)*sb             +48.0_lp*sr &
         -22.0_lp*sa                         +(17.0_lp-0.1_lp*t)*sp2 &
         -15.0_lp*(sp*cn+cp*sn)              -(16.0_lp-0.1_lp*t)*(sc*cp2+cc*sp2) &
         -12.0_lp*(sn*cp-cn*sp)              -6.0_lp*(sn*cr-cn*sr) &
         -5.0_lp*(sb*cp-cb*sp)               +4.0_lp*(sr*cn+cr*sn) &
         +4.0_lp*(sb*cp+cb*sp)               -4.0_lp*sq &
         +(sr*cp+cr*sp)                      +(sn*ca-cn*sa) &
         -(sp*ca-cp*sa)                      +(sp*cn2+cp*sn2) &
         +(sn*cq-cn*sq)                      -(sp*ca+cp*sa) &
         -(2274.0_lp+0.2*t)*sh               +(712.0_lp+0.1_lp*t)*sl &
         -(386.0_lp+0.4*t)*ss                -301.0_lp*sj &
         -158.0_lp*sg                        +123.0_lp*(sh*cl-ch*sl) &
         +63.0_lp*sd2                        +(63.0_lp+0.1_lp*t)*(sl*cn+cl*sn) &
         -(58.0_lp+0.1_lp*t)*(sn*cl-cn*sl)   -59.0_lp*su &
         -51.0_lp*st                         -38.0_lp*sf &
         +29.0_lp*sl2                        +29.0_lp*(sc*cl+cc*sl) &
         -31.0_lp*sk                         +26.0_lp*sx &
         +21.0_lp*(ss*cl-cs*sl)              +16.0_lp*(sn*cg-cn*sg) &
         -13.0_lp*(sn*cg+cn*sg)              -10.0_lp*(se*cl-ce*sl) &
         -7.0_lp*(sg*cp+cg*sp)               +7.0_lp*(sh*cp+ch*sp) &
         -7.0_lp*(sh*cp-ch*sp)               -8.0_lp*(sf*cl+cf*sl)
    ! Terms 55-106
    dpsi = dpsi +6.0_lp*(sl*cd2+cl*sd2)      +6.0_lp*(sc*cl2+cc*sl2) &
         -6.0_lp*(sn*cd2+cn*sd2)             -7.0_lp*se &
         +6.0_lp*(sb*cl+cb*sl)               -5.0_lp*(sn*cd2-cn*sd2) &
         +5.0_lp*(sl*cp-cl*sp)               -5.0_lp*(ss*cl2+cs*sl2) &
         -4.0_lp*(sp*cd2-cp*sd2)             +4.0_lp*(sl*cx-cl*sx) &
         -4.0_lp*sd                          -3.0_lp*(sl*cp+cl*sp) &
         +3.0_lp*(sl*cx+cl*sx)               -3.0_lp*(sj*cp-cj*sp) &
         -3.0_lp*(su*cp-cu*sp)               -2.0_lp*(sn*cl2-cn*sl2) &
         -3.0_lp*(sk*cl+ck*sl)               -3.0_lp*(sf*cp-cf*sp) &
         +2.0_lp*(sj*cp+cj*sp)               -2.0_lp*(sb*cl-cb*sl) &
         +2.0_lp*(sn*cl2+cn*sl2)             -2.0_lp*(sl*cn2+cl*sn2) &
         +2.0_lp*(sl*cl2+cl*sl2)             +2.0_lp*(sh*cd+ch*sd) &
         +(sn2*cl-cn2*sl)                    -(sg*cd2-cg*sd2) &
         +(sf*cl2-cf*sl2)                    -2.0_lp*(su*cd2+cu*sd2) &
         -(sr*cd2-cr*sd2)                    +(sw*ch+cw*sh) &
         -(sl*ce+cl*se)                      -(sf*cr-cf*sr) &
         +(su*ca+cu*sa)                      +(sg*cp-cg*sp) &
         +(sb*cl2+cb*sl2)                    -(sf*cl2+cf*sl2) &
         -(st*ca-ct*sa)                      +(sc*cx+cc*sx) &
         +(sj*cr+cj*sr)                      -(sg*cx+cg*sx) &
         +(sp*cs+cp*ss)                      +(sn*cw-cn*sw) &
         -(sn*cx-cn*sx)                      -(sh*cd-ch*sd) &
         -(sp*cd2+cp*sd2)                    -(sl*cv-cl*sv) &
         -(ss*cp-cs*sp)                      -(sw*cn+cw*sn) &
         -(sl*ca-cl*sa)                      +(sl2*cd2+cl2*sd2) &
         -(sf*cd2+cf*sd2)                    +(sp*cd+cp*sd)

    ! Series for deps:
    deps = (92025.0_lp+8.9*t)*cn             -(895.0_lp-0.5*t)*cn2 &
         -24.0_lp*(cm*cn-sm*sn)              +(cm*cn2-sm*sn2) &
         +(cb*cp2+sb*sp2)                    +(5736.0_lp-3.1_lp*t)*cc &
         +(54.0_lp-0.1_lp*t)*cp              +(224.0_lp-0.6*t)*(cc*cp-sc*sp) &
         -(95.0_lp-0.3*t)*(cc*cp+sc*sp)      -70.0_lp*cb &
         +cr                                 +9.0_lp*(cp*cn-sp*sn) &
         +7.0_lp*(cc*cp2-sc*sp2)             +6.0_lp*(cn*cp+sn*sp) &
         +3.0_lp*(cn*cr+sn*sr)               +3.0_lp*(cb*cp+sb*sp) &
         -2.0_lp*(cr*cn-sr*sn)               -2.0_lp*(cb*cp-sb*sp) &
         +(977.0_lp-0.5*t)*ch                -7.0_lp*cl &
         +200.0_lp*cs                        +(129.0_lp-0.1_lp*t)*cj &
         -cg                                 -53.0_lp*(ch*cl+sh*sl) &
         -2.0_lp*cd2                         -33.0_lp*(cl*cn-sl*sn) &
         +32.0_lp*(cn*cl+sn*sl)              +26.0_lp*cu &
         +27.0_lp*ct                         +16.0_lp*cf &
         -cl2                                -12.0_lp*(cc*cl-sc*sl)
    deps = deps &
         +13.0_lp*ck                         -cx &
         -10.0_lp*(cs*cl+ss*sl)              -8.0_lp*(cn*cg+sn*sg) &
         +7.0_lp*(cn*cg-sn*sg)               +5.0_lp*(ce*cl+se*sl) &
         -3.0_lp*(ch*cp-sh*sp)               +3.0_lp*(ch*cp+sh*sp) &
         +3.0_lp*(cf*cl-sf*sl)               -3.0_lp*(cc*cl2-sc*sl2) &
         +3.0_lp*(cn*cd2-sn*sd2)             +3.0_lp*ce &
         -3.0_lp*(cb*cl-sb*sl)               +3.0_lp*(cn*cd2+sn*sd2) &
         +3.0_lp*(cs*cl2-ss*sl2)             +(cj*cp+sj*sp) &
         +(cu*cp+su*sp)                      +(cn*cl2+sn*sl2) &
         +(ck*cl-sk*sl)                      +(cf*cp+sf*sp) &
         -(cj*cp-sj*sp)                      +(cb*cl+sb*sl) &
         -(cn*cl2-sn*sl2)                    +(cl*cn2-sl*sn2) &
         -(ch*cd-sh*sd)                      -(cn2*cl+sn2*sl) &
         -(cf*cl2+sf*sl2)                    +(cu*cd2-su*sd2) &
         -(cw*ch-sw*sh)                      +(cl*ce-sl*se) &
         +(cf*cr+sf*sr)                      -(cb*cl2-sb*sl2) 

    dpsi = dpsi*0.0001_bp
    deps = deps*0.0001_bp

  END SUBROUTINE getNutationAngles






  !! *Description*:
  !!
  !! Computes nutation matrix according to Wahr (IAU-1980) theory.
  !! The nutation matrix (RNUT) transforms mean coordinates into true
  !! coordinates:
  !!
  !!     Xtrue = RNUT Xmean
  !!
  !! Use transpose to get from Xtrue to Xmean:
  !!
  !!     Xmean = transpose(RNUT) Xtrue
  !!
  !! Based on f77 routine by Mario Carpino (rnut80.f).
  !!
  !! Returns error.
  !!
  FUNCTION getNutationMatrix(t)

    IMPLICIT NONE

    TYPE (Time), INTENT(inout) :: t
    REAL(bp), DIMENSION(3,3) :: getNutationMatrix
    REAL(bp), DIMENSION(3,3) :: r1, r2, r3, r23
    REAL(bp) :: epsm, epst, dpsi, deps, mjd_tt

    IF (.NOT. exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getNutationMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd_tt = getMJD(t, "TT")
    IF (error) THEN
       CALL errorMessage("Observatories / getNutationMatrix", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    epsm = meanObliquity(mjd_tt)
    CALL getNutationAngles(t, dpsi, deps)
    IF (error) THEN
       CALL errorMessage("Observatories / getNutationMatrix", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    dpsi = dpsi * rad_asec
    epst = epsm + deps * rad_asec

    r3 = rotationMatrix(  epsm, 1)
    r2 = rotationMatrix( -dpsi, 3)
    r1 = rotationMatrix( -epst, 1)

    r23 = MATMUL(r2,r3)
    getNutationMatrix = MATMUL(r1,r23)

  END FUNCTION getNutationMatrix






  !! *Description*:
  !!
  !! Returns precession matrix which transforms equatorial coordinates
  !! from J2000.0 into mean equatorial coordinates at chosen epoch t:
  !!
  !!     Xepoch = RPREC Xj2000.0
  !!
  !! Use transpose to get from Xepoch to Xj2000.0:
  !!
  !!     Xj2000.0 = transpose(RPREC) Xepoch
  !!
  !! Based on f77 routine by Mario Carpino (prec.f)
  !!
  !! Returns error.
  !!
  FUNCTION getPrecessionMatrix(t)

    IMPLICIT NONE
    TYPE(Time), INTENT(inout) :: t
    REAL(bp), DIMENSION(3,3)  :: getPrecessionMatrix
    REAL(bp), DIMENSION(3,3) :: r1, r2, r3, r12
    REAL(bp) :: mjd_tt, cent, zeta, theta, z
    REAL(bp), PARAMETER :: &
                                ! Calcolo costanti usate (vedi Astronomical Almanac 1987, B18)
                                ! Linear terms:
         zed = 0.6406161_bp * rad_deg, &
         zd  = 0.6406161_bp * rad_deg, &
         thd = 0.5567530_bp * rad_deg, &
                                ! Quadratical terms:
         zedd= 0.0000839_bp * rad_deg, &
         zdd = 0.0003041_bp * rad_deg, &
         thdd = -0.0001185_bp * rad_deg, &
                                ! Cubical terms:
         zeddd = 0.0000050_bp * rad_deg, &
         zddd = 0.0000051_bp * rad_deg, &
         thddd = - 0.0000116_bp * rad_deg

    IF (.NOT. exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / getPrecessionMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd_tt = getMJD(t, "TT")
    IF (error) THEN
       CALL errorMessage("Observatories / getPrecessionMatrix", &
            "TRACE BACK", 1)
       RETURN
    END IF

    ! Fundamental arguments:
    cent  = ( mjd_tt - 51544.5_bp ) / 36525.0_bp ! centuries since J2000.0
    zeta  = ( ( zeddd * cent + zedd ) * cent + zed ) * cent
    z     = ( (  zddd * cent +  zdd ) * cent +  zd ) * cent
    theta = ( ( thddd * cent + thdd ) * cent + thd ) * cent

    r3 = rotationMatrix(- zeta, 3)
    r2 = rotationMatrix( theta, 2)
    r1 = rotationMatrix(-    z, 3)

    ! Rotation matrix:
    r12 = MATMUL(r1,r2)
    getPrecessionMatrix = MATMUL(r12,r3)

  END FUNCTION getPrecessionMatrix






  !! *Description*:
  !! 
  !! Construction of rotation matrix for transformations between
  !! different reference systems.
  !!
  !! INPUT:    frame1    -  starting frame ("ecliptical" or "equatorial")
  !!           rsys1     -  starting reference system ("mean" or "true-of-date")
  !!           t1        -  starting reference time (mjd, tt)
  !!           frame2    -  final frame ("ecliptical" or "equatorial")
  !!           rsys2     -  final reference system ("mean" or "true-of-date")
  !!           t2        -  final reference time (mjd, tt)
  !!
  !! OUTPUT:   ROT(3,3)  -  rotation matrix giving the transformation from
  !!                        starting to final reference systems:
  !!                        x2 = ROT x1
  !!
  !! Returns error.
  !!
  !! Based on f77 routine by Mario Carpino (rotpn.f).
  !!
  FUNCTION precessionAndNutationMatrix(frame1, rsys1, t1, frame2, rsys2, t2)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout) :: t1, t2
    CHARACTER(len=*), INTENT(in) :: frame1, rsys1, frame2, rsys2
    REAL(bp), DIMENSION(3,3) :: precessionAndNutationMatrix
    TYPE (Time) :: t_
    CHARACTER(len=20) :: frame, rsys
    REAL(bp), DIMENSION(3,3) :: rot, r
    REAL(bp), PARAMETER :: eps = 1.0e-6_bp
    REAL(bp) :: tt_, tt1, tt2, obl

    IF (frame1 /= "ecliptical" .AND. frame1 /= "equatorial") THEN
       error = .TRUE.
       CALL errorMessage("Observatories / precessionAndNutationMatrix", &
            "Unkown starting frame.", 1)
       RETURN
    END IF

    IF (rsys1 /= "mean" .AND. rsys1 /= "true-of-date") THEN
       error = .TRUE.
       CALL errorMessage("Observatories / precessionAndNutationMatrix", &
            "Unkown starting reference system.", 1)
       RETURN
    END IF

    IF (frame2 /= "ecliptical" .AND. frame2 /= "equatorial") THEN
       error = .TRUE.
       CALL errorMessage("Observatories / precessionAndNutationMatrix", &
            "Unkown final frame.", 1)
       RETURN
    END IF

    IF (rsys2 /= "mean" .AND. rsys2 /= "true-of-date") THEN
       error = .TRUE.
       CALL errorMessage("Observatories / precessionAndNutationMatrix", &
            "Unkown final reference system.", 1)
       RETURN
    END IF

    IF (.NOT. exist(t1)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / precessionAndNutationMatrix", &
            "Unkown starting epoch.", 1)
       RETURN
    END IF

    IF (.NOT. exist(t2)) THEN
       error = .TRUE.
       CALL errorMessage("Observatories / precessionAndNutationMatrix", &
            "Unkown final epoch.", 1)
       RETURN
    END IF

    frame = frame1
    rsys = rsys1
    t_ = copy(t1)

    tt1 = getMJD(t1, "TT")
    IF (error) THEN
       CALL errorMessage("Observatories / " // &
            "precessionAndNutationMatrix", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    tt2 = getMJD(t2, "TT")
    IF (error) THEN
       CALL errorMessage("Observatories / " // &
            "precessionAndNutationMatrix", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    tt_ = tt1

    ! Initialization of the rotation matrix (equal to the unit matrix)
    rot = identity_matrix(3)

    IF (info_verb >= 5) THEN
       CALL matrix_print(rot,stdout,errstr)
       WRITE(stdout,*)
    END IF

    ! Building of the rotation matrix
    DO

       ! Different epochs
       IF (ABS(tt_ - tt2) > eps) THEN

          IF (frame == "ecliptical") THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,*) "Transformation of FRAME (ecliptical --> equatorial)"
             END IF
             obl = meanObliquity(tt_)
             r = rotationMatrix(-obl,1)
             rot = MATMUL(r,rot)
             frame = "equatorial"
             IF (error) THEN
                CALL errorMessage("Observatories / " // &
                     "precessionAndNutationMatrix", &
                     "TRACE BACK (15)", 1)
                RETURN
             END IF
          ELSE IF (frame == "equatorial") THEN
             IF (TRIM(rsys) == "true-of-date") THEN
                IF (info_verb >= 5) THEN
                   WRITE(stdout,*) "Transformation of RSYS (True-of-date --> mean)"
                END IF
                r = getNutationMatrix(t_)
                IF (error) THEN
                   CALL errorMessage("Observatories / " // &
                        "precessionAndNutationMatrix", &
                        "TRACE BACK (20)", 1)
                   RETURN
                END IF
                rot = MATMUL(TRANSPOSE(r),rot)
                rsys = "mean"
             ELSE IF (rsys == "mean") THEN
                IF (info_verb >= 5) THEN
                   WRITE(stdout,*) "Transformation of T (precession)"
                END IF
                r = getPrecessionMatrix(t_)
                rot = MATMUL(TRANSPOSE(r),rot)
                r = getPrecessionMatrix(t2)
                rot = MATMUL(r,rot)
                tt_ = tt2
                t_ = copy(t2)
                IF (error) THEN
                   CALL errorMessage("Observatories / " // &
                        "precessionAndNutationMatrix", &
                        "TRACE BACK (25)", 1)
                   RETURN
                END IF
             END IF
          END IF
          !
          ! Already at the same epoch
          !
       ELSE IF (TRIM(rsys) /= TRIM(rsys2)) THEN

          IF (TRIM(rsys) == "true-of-date") THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,*) "Transformation of RSYS (True-of-date --> mean)"
             END IF
             r = getNutationMatrix(t_)
             IF (error) THEN
                CALL errorMessage("Observatories / " // &
                     "precessionAndNutationMatrix", &
                     "TRACE BACK (30)", 1)
                RETURN
             END IF
             rot = MATMUL(TRANSPOSE(r),rot)
             rsys = "mean"
          ELSE IF (rsys == "mean") THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,*) "Transformation of RSYS (Mean --> True-of-date)"
             END IF
             IF (TRIM(rsys2) == "true-of-date") THEN
                r = getNutationMatrix(t_)
                IF (error) THEN
                   CALL errorMessage("Observatories / " // &
                        "precessionAndNutationMatrix", &
                        "TRACE BACK (35)", 1)
                   RETURN
                END IF
                rot = MATMUL(r,rot)
                rsys = rsys2
             ELSE
                error = .TRUE.
                CALL errorMessage("Observatories / " // &
                     "precessionAndNutationMatrix", &
                     "Internal error (5).", 1)
                RETURN
             END IF
          END IF

       ELSE IF (frame /= frame2) THEN

          IF (frame == "ecliptical") THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,*) "Transformation of FRAME (ecliptical --> equatorial)"
             END IF
             obl = meanObliquity(tt_)
             r = rotationMatrix(-obl,1)
             rot = MATMUL(r,rot)
             frame = "equatorial"
          ELSE IF (frame == "equatorial") THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,*) "Transformation of FRAME (equatorial --> ecliptical)"
             END IF
             IF (frame2 == "ecliptical") THEN
                obl = meanObliquity(tt_)
                r = rotationMatrix(obl,1)
                rot = MATMUL(r,rot)
                frame = "ecliptical"
             ELSE
                error = .TRUE.
                CALL errorMessage("Observatories / " // &
                     "precessionAndNutationMatrix", &
                     "Internal error (10).", 1)
                RETURN
             END IF
          END IF

       ELSE

          EXIT

       END IF

       IF (info_verb >= 5) THEN
          CALL matrix_print(rot,stdout,errstr)
          WRITE(stdout,*)
       END IF

    END DO

    precessionAndNutationMatrix = rot

  END FUNCTION precessionAndNutationMatrix




END MODULE Observatories_cl











