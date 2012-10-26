!====================================================================!
!                                                                    !
! Copyright 2002-2011,2012                                           !
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
!! Generic type and routines for an astrometric observation, which can
!! be presented with spherical coordinates and their time derivatives.
!!  
!! @see Observations_class 
!!  
!! @author  MG, JV 
!! @version 2012-10-26
!!  
MODULE Observation_cl

  USE Base_cl
  USE Time_cl
  USE Observatory_cl
  USE SphericalCoordinates_cl
  USE CartesianCoordinates_cl

  USE utilities

  IMPLICIT NONE
  PRIVATE :: NEW_Obs
  PRIVATE :: NULLIFY_Obs
  PRIVATE :: copy_Obs
  PRIVATE :: exist_Obs
  PRIVATE :: reallocate_a1_Obs
  PRIVATE :: reallocate_a2_Obs
  PRIVATE :: addMultinormalDeviate_Obs
  PRIVATE :: addUniformDeviate_Obs
  PRIVATE :: setDesignation_Obs
  PRIVATE :: equal_Obs
  PRIVATE :: getRA_Obs
  PRIVATE :: getDec_Obs
  PRIVATE :: getNumber_Obs
  PRIVATE :: getTime_Obs
  PRIVATE :: getCode_Obs
  PRIVATE :: getCovarianceMatrix_Obs
  PRIVATE :: getObservationMask_Obs
  PRIVATE :: getObservatory_Obs
  PRIVATE :: getObservatoryCCoord_Obs
  PRIVATE :: getObservatoryCode_Obs
  PRIVATE :: getStandardDeviations_Obs
  PRIVATE :: setNumber_Obs

  TYPE Observation
     PRIVATE
     INTEGER                        :: number
     CHARACTER(len=DESIGNATION_LEN) :: designation
     LOGICAL                        :: discovery
     CHARACTER(len=1)               :: note1, note2, mode
     CHARACTER(len=2)               :: filter
     TYPE (SphericalCoordinates)    :: obs_scoord
     REAL(bp), DIMENSION(6,6)       :: covariance
     LOGICAL, DIMENSION(6)          :: obs_mask
     REAL(bp)                       :: mag
     REAL(bp)                       :: mag_unc
     REAL(bp)                       :: s2n
     CHARACTER(len=DESIGNATION_LEN) :: secret_name
     TYPE (Observatory)             :: obsy
     TYPE (CartesianCoordinates)    :: obsy_ccoord
     TYPE (CartesianCoordinates)    :: satellite_ccoord
     INTEGER                        :: coord_unit
     LOGICAL                        :: is_initialized = .FALSE.
  END TYPE Observation

  INTERFACE NEW
     MODULE PROCEDURE NEW_Obs
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE NULLIFY_Obs
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_Obs
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_Obs
  END INTERFACE exist

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_a1_Obs
     MODULE PROCEDURE reallocate_a2_Obs
  END INTERFACE reallocate

  INTERFACE addMultinormalDeviate
     MODULE PROCEDURE addMultinormalDeviate_Obs
  END INTERFACE addMultinormalDeviate

  INTERFACE addUniformDeviate
     MODULE PROCEDURE addUniformDeviate_Obs
  END INTERFACE addUniformDeviate

  INTERFACE getDesignation
     MODULE PROCEDURE getDesignation_Obs
  END INTERFACE getDesignation

  INTERFACE setDesignation
     MODULE PROCEDURE setDesignation_Obs
  END INTERFACE setDesignation

  INTERFACE equal
     MODULE PROCEDURE equal_Obs
  END INTERFACE equal

  INTERFACE getDec
     MODULE PROCEDURE getDec_Obs
  END INTERFACE getDec

  INTERFACE getID
     MODULE PROCEDURE getID_Obs
  END INTERFACE getID

  INTERFACE getNumber
     MODULE PROCEDURE getNumber_Obs
  END INTERFACE getNumber

  INTERFACE getRA
     MODULE PROCEDURE getRA_Obs
  END INTERFACE getRA

  INTERFACE getTime
     MODULE PROCEDURE getTime_Obs
  END INTERFACE getTime

  INTERFACE getCode
     MODULE PROCEDURE getCode_Obs
  END INTERFACE getCode

  INTERFACE getCovarianceMatrix
     MODULE PROCEDURE getCovarianceMatrix_Obs
  END INTERFACE getCovarianceMatrix

  INTERFACE getObservationMask
     MODULE PROCEDURE getObservationMask_Obs
  END INTERFACE getObservationMask

  INTERFACE getObservatory
     MODULE PROCEDURE getObservatory_Obs
  END INTERFACE getObservatory

  INTERFACE getObservatoryCCoord
     MODULE PROCEDURE getObservatoryCCoord_Obs
  END INTERFACE getObservatoryCCoord

  INTERFACE getObservatoryCode
     MODULE PROCEDURE getObservatoryCode_Obs
  END INTERFACE getObservatoryCode

  INTERFACE getStandardDeviations
     MODULE PROCEDURE getStandardDeviations_Obs
  END INTERFACE getStandardDeviations

  INTERFACE setNumber
     MODULE PROCEDURE setNumber_Obs
  END INTERFACE setNumber

CONTAINS




  !! *Description*:
  !!
  !! Initializes this object with the given values.
  !!
  !! Returns error.
  !!
  SUBROUTINE NEW_Obs(this, number, designation, discovery, note1, &
       note2, obs_scoord, covariance, obs_mask, mag, mag_unc, filter, &
       s2n, obsy, obsy_ccoord, satellite_ccoord, coord_unit, secret_name)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout)                 :: this
    INTEGER, INTENT(in)                               :: number
    CHARACTER(len=*), INTENT(in)                      :: designation
    LOGICAL, INTENT(in)                               :: discovery
    CHARACTER(len=*), INTENT(in)                      :: note1, note2
    CHARACTER(len=*), INTENT(in)                      :: filter
    TYPE (SphericalCoordinates), INTENT(in)           :: obs_scoord
    REAL(bp), DIMENSION(6,6), INTENT(in)              :: covariance
    LOGICAL, DIMENSION(6), INTENT(in)                 :: obs_mask
    REAL(bp), INTENT(in)                              :: mag
    REAL(bp), INTENT(in), OPTIONAL                    :: mag_unc
    REAL(bp), INTENT(in), OPTIONAL                    :: s2n
    TYPE (Observatory), INTENT(in)                    :: obsy
    TYPE (CartesianCoordinates), INTENT(in)           :: obsy_ccoord
    TYPE (CartesianCoordinates), INTENT(in), OPTIONAL :: satellite_ccoord
    INTEGER, INTENT(in), OPTIONAL                     :: coord_unit
    CHARACTER(len=*), INTENT(in), OPTIONAL            :: secret_name

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%number             = number
    IF (LEN_TRIM(designation) <= designation_len) THEN
       this%designation = TRIM(designation)
    ELSE
       error = .TRUE.
       CALL errorMessage("Observation / new", &
            "Too many characters in designation:" // TRIM(designation), 1)
       RETURN
    END IF
    this%discovery         = discovery
    this%note1             = note1
    this%note2             = note2
    this%obs_scoord        = copy(obs_scoord)
    this%covariance        = covariance
    this%obs_mask          = obs_mask
    this%mag               = mag
    IF (PRESENT(mag_unc)) THEN
       this%mag_unc        = mag_unc
    END IF
    this%filter              = filter
    IF (PRESENT(s2n)) THEN
       this%s2n            = s2n
    END IF
    this%obsy              = copy(obsy)
    this%obsy_ccoord       = copy(obsy_ccoord)
    IF ((note2 == "S" .NEQV. PRESENT(satellite_ccoord)) .OR. &
         (note2 == "S" .NEQV. PRESENT(coord_unit))) THEN
       error = .TRUE.
       CALL errorMessage("Observation / new", &
            "Input of satellite (?) observation is inconsistent.", 1)
       RETURN
    END IF
    IF (PRESENT(satellite_ccoord)) THEN
       this%satellite_ccoord = copy(satellite_ccoord)
    END IF
    IF (PRESENT(coord_unit)) THEN
       this%coord_unit = coord_unit
    END IF
    IF (PRESENT(secret_name)) THEN
       this%secret_name = secret_name
    END IF
    this%is_initialized = .TRUE.

  END SUBROUTINE NEW_Obs





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE NULLIFY_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout) :: this

    this%number            = 0
    this%designation       = ""
    this%discovery         = .FALSE.
    this%note1             = ""
    this%note2             = ""
    CALL NULLIFY(this%obs_scoord)
    this%covariance = 0.0_bp
    this%obs_mask          = .FALSE.
    this%mag               = 0.0_bp
    this%mag_unc           = -1.0_bp
    this%filter              = "  "
    this%s2n               = -1.0_bp
    CALL NULLIFY(this%obsy)
    CALL NULLIFY(this%obsy_ccoord)
    CALL NULLIFY(this%satellite_ccoord)
    this%coord_unit        = 0
    this%secret_name       = " "
    this%is_initialized    = .FALSE.

  END SUBROUTINE NULLIFY_Obs





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    TYPE (Observation)             :: copy_Obs

    copy_Obs%number            = this%number
    copy_Obs%designation       = this%designation
    copy_Obs%discovery         = this%discovery
    copy_Obs%note1             = this%note1
    copy_Obs%note2             = this%note2
    copy_Obs%obs_scoord        = copy(this%obs_scoord)
    IF (error) THEN
       CALL errorMessage("Observation / copy", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    copy_Obs%covariance = this%covariance
    copy_Obs%obs_mask          = this%obs_mask
    copy_Obs%mag               = this%mag
    copy_Obs%mag_unc           = this%mag_unc
    copy_Obs%filter              = this%filter
    copy_Obs%s2n               = this%s2n
    copy_Obs%obsy              = copy(this%obsy)
    IF (error) THEN
       CALL errorMessage("Observation / copy", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    copy_Obs%obsy_ccoord       = copy(this%obsy_ccoord)
    IF (error) THEN
       CALL errorMessage("Observation / copy", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    copy_Obs%satellite_ccoord  = copy(this%satellite_ccoord)
    IF (error) THEN
       CALL errorMessage("Observation / copy", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF
    copy_Obs%coord_unit        = this%coord_unit
    copy_Obs%secret_name       = this%secret_name
    copy_Obs%is_initialized    = this%is_initialized

  END FUNCTION copy_Obs





  !! *Description*:
  !!
  !! Returns the status of this object, i.e. whether it exists or not.
  !!
  LOGICAL FUNCTION exist_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    exist_Obs = this%is_initialized

  END FUNCTION exist_Obs





  !! *Description*:
  !!
  !! Reallocates a pointer array of Observation-objects and
  !! copies the existing (non-nullified) data from the old array to
  !! the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myobservations => reallocate(myobservations,4)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_a1_Obs(array,n)

    IMPLICIT NONE
    TYPE (Observation), DIMENSION(:), POINTER :: reallocate_a1_Obs, array
    INTEGER, INTENT(in), OPTIONAL             :: n
    INTEGER                                   :: i, nn, nold, err, nonexistent

    nonexistent = 0
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array,dim=1)
    ELSE
       nold = 0
    END IF
    IF (PRESENT(n)) THEN
       nn = n
    ELSE
       nn = 0
       DO i=1,nold
          IF (array(i)%is_initialized) THEN
             nn = nn + 1
          END IF
       END DO
    END IF
    IF (nn /= 0) THEN
       ALLOCATE(reallocate_a1_Obs(nn), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / reallocate", &
               "Could not allocate memory.", 1)
          reallocate_a1_Obs => NULL()
          RETURN
       END IF
    ELSE
       reallocate_a1_Obs => NULL()
       DO i=1,SIZE(array)
          CALL NULLIFY(array(i))
       END DO
       DEALLOCATE(array, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / reallocate", &
               "Could not deallocate memory.", 1)
          RETURN
       END IF
    END IF

    ! Return if there isn't anything to copy:
    IF (.NOT.ASSOCIATED(array)) THEN
       RETURN
    END IF

    ! Scan the whole initial array if n isn't present:
    IF (.NOT.PRESENT(n)) THEN
       nn = nold
    END IF
    DO i=1,MIN(nn,nold)
       IF (.NOT.array(i)%is_initialized) THEN
          nonexistent = nonexistent + 1
          CYCLE
       END IF
       reallocate_a1_Obs(i-nonexistent) = copy(array(i))
    END DO
    DO i=1,SIZE(array)
       CALL NULLIFY(array(i))
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observation / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_a1_Obs





  !! *Description*:
  !!
  !! Reallocates a pointer array of Observation-objects and
  !! copies existing data from the old array to
  !! the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myobservations => reallocate(myobservations,4,5)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_a2_Obs(array,n,m)

    IMPLICIT NONE
    TYPE (Observation), DIMENSION(:,:), POINTER :: reallocate_a2_Obs, array
    INTEGER, INTENT(in)                         :: n, m
    INTEGER                                     :: i, j, nold, mold, err

    nold = SIZE(array,dim=1)
    mold = SIZE(array,dim=2)
    ALLOCATE(reallocate_a2_Obs(n,m), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observation / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_a2_Obs => NULL()
       RETURN
    END IF

    ! Return if there isn't anything to copy:
    IF (.NOT.ASSOCIATED(array)) THEN
       RETURN
    END IF
    DO i=1,MIN(n,nold)
       DO j=1,MIN(m,mold)
          reallocate_a2_Obs(i,j) = copy(array(i,j))
       END DO
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observation / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_a2_Obs





  !! *Description*:
  !!
  !! Adds multinormal deviates to the observation. The mean relative
  !! to the original coordinates and the (optional) covariance matrix
  !! should be given in radians for the angular space coordinates, AUs
  !! for distance, radians per day for angular velocities and AUs per
  !! day for line-of-sight velocity. Note that the covariance matrix
  !! for the observational uncertainties is not correctly updated if
  !! the original observation has a non-zero covariance matrix.
  !!
  !! Returns error.
  !!
  SUBROUTINE addMultinormalDeviate_Obs(this, mean, covariance)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout)    :: this
    REAL(bp), DIMENSION(6), INTENT(in)   :: mean
    REAL(bp), DIMENSION(6,6), INTENT(in) :: covariance

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / addMultinormalDeviate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    CALL addMultinormalDeviate(this%obs_scoord, mean, covariance)
    IF (error) THEN
       CALL errorMessage("Observation / addMultinormalDeviate", &
            "TRACE BACK", 1)
       RETURN
    END IF
    ! New noise estimate = old noise estimate + added noise
    this%covariance = this%covariance + covariance

  END SUBROUTINE addMultinormalDeviate_Obs





  !! *Description*:
  !!
  !! Adds uniform deviates to the observation. The center relative to
  !! the original coordinates and the absolute values of the boundary
  !! values ((i,1)=center, (i,2)=abs(boundary), position i=1:3 and
  !! velocity i=4:6) should be given in radians for the angular space
  !! coordinates, AUs for distance, radians per day for angular
  !! velocities and AUs per day for line-of-sight velocity. Note that
  !! the covariance matrix for the observational uncertainties is not
  !! updated.
  !!
  !! Returns error.
  !!
  SUBROUTINE addUniformDeviate_Obs(this, center_and_absbound)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout)    :: this
    REAL(bp), DIMENSION(6,2), INTENT(in) :: center_and_absbound

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / addUniformDeviate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    CALL addUniformDeviate(this%obs_scoord, center_and_absbound)
    IF (error) THEN
       CALL errorMessage("Observation / addUniformDeviate", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END SUBROUTINE addUniformDeviate_Obs





  LOGICAL FUNCTION equal_Obs(this, that)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    TYPE (Observation), INTENT(in) :: that

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / equal", &
            "1st object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. that%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / equal", &
            "2nd object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Designation:
    IF (.NOT.(this%designation == that%designation)) THEN
       equal_Obs = .FALSE.
       RETURN
    END IF

    ! Epoch and coordinates:
    IF (.NOT.equal(this%obs_scoord,that%obs_scoord)) THEN
       equal_Obs = .FALSE.
       RETURN
    END IF

    ! Observatory
    IF (.NOT.equal(this%obsy,that%obsy)) THEN
       equal_Obs = .FALSE.
       RETURN
    END IF

    ! Assuming that all of the above comparisons are true,
    ! it can be concluded that the two objects are the same:
    equal_Obs = .TRUE.

  END FUNCTION equal_Obs





  !! *Description*:
  !!
  !! Returns the code of the observatory where the observation was
  !! made.
  !!
  !! Returns error.
  !!
  CHARACTER(len=OBSY_CODE_LEN) FUNCTION getCode_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObsCode", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getCode_Obs = getCode(this%obsy)

  END FUNCTION getCode_Obs





  !! *Description*:
  !!
  FUNCTION getCovarianceMatrix_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    REAL(bp), DIMENSION(6,6)       :: getCovarianceMatrix_Obs

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getCovarianceMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getCovarianceMatrix_Obs = this%covariance

  END FUNCTION getCovarianceMatrix_Obs





  !! *Description*:
  !!
  !! Returns Declination [rad].
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getDec_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getDec", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getDec_Obs = getLatitude(this%obs_scoord)
    IF (error) THEN
       CALL errorMessage("Observation / getDec", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END FUNCTION getDec_Obs





  !! *Description*:
  !!
  !! Returns designation of this object.
  !!
  !! Returns error.
  !!
  CHARACTER(len=DESIGNATION_LEN) FUNCTION getDesignation_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getDesignation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getDesignation_Obs = TRIM(this%designation)

  END FUNCTION getDesignation_Obs





  !! *Description*:
  !!
  !! Returns filter.
  !!
  !! Returns error.
  !!
  CHARACTER(len=2) FUNCTION getFilter(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getFilter", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getFilter = this%filter

  END FUNCTION getFilter





  !! *Description*:
  !!
  !! Returns the ID (number or designation) for an single object object.
  !!
  !! Returns error.
  !!
  FUNCTION getID_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    CHARACTER(len=DESIGNATION_LEN) :: getID_Obs

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getID", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%number == 0) THEN
       getID_Obs = this%designation
    ELSE
       CALL toString(this%number, getID_Obs, error)
       IF (error) THEN
          CALL errorMessage("Observation / getID", &
               "Could not convert integer to string.", 1)
          RETURN
       END IF
    END IF

  END FUNCTION getID_Obs





  !! *Description*:
  !!
  !! Returns magnitude [mag].
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getMagnitude(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getMagnitude", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getMagnitude = this%mag

  END FUNCTION getMagnitude





  !! *Description*:
  !!
  !! Returns magnitude uncertainty [mag].
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getMagnitudeUncertainty(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getMagnitudeUncertainty", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getMagnitudeUncertainty = this%mag_unc

  END FUNCTION getMagnitudeUncertainty





  !! *Description*:
  !!
  !! Returns number of the asteroid corresponding to this object. If
  !! the asteroid is unnumbered, the function returns the default
  !! value (usually 0).
  !!
  !! Returns error.
  !!
  INTEGER FUNCTION getNumber_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getNumber", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getNumber_Obs = this%number

  END FUNCTION getNumber_Obs





  !! *Description*:
  !!
  !! Returns the observation records as character strings of a given
  !! format (old and new MPC, Bowell) containing the data of this
  !! object.
  !!
  !! Returns error.
  !!
  FUNCTION getObservationRecords(this, frmt, number) RESULT(records)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in)                       :: this
    CHARACTER(len=*), INTENT(in)                         :: frmt
    INTEGER, INTENT(in), OPTIONAL                        :: number
    CHARACTER(len=OBS_RECORD_LEN), DIMENSION(:), POINTER :: records

    TYPE (Time)                       :: t
    TYPE (CartesianCoordinates)       :: ccoord
    CHARACTER(len=DESIGNATION_LEN)    :: designation, secret_name
    CHARACTER(len=12)                 :: day_str
    CHARACTER(len=9)                  :: ra_unc, dec_unc
    CHARACTER(len=7)                  :: number_str
    CHARACTER(len=6)                  :: s_str, mag_str
    CHARACTER(len=5)                  :: as_str
    CHARACTER(len=4)                  :: obsy_code
    CHARACTER(len=2)                  :: month_str, h_str, m_str, deg_str, am_str, filter_
    CHARACTER(len=1)                  :: discovery, sign_dec, sign_x, sign_y, sign_z
    REAL(bp), DIMENSION(6)            :: coord
    REAL(bp)                          :: day, ra, s, dec, as, correlation
    INTEGER(ihp)                      :: day_integer
    INTEGER                           :: number_, err, year, month, h, m, &
         deg, am, n, i, indx

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObservationRecords", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    t = getTime(this)
    IF (error) THEN
       CALL errorMessage("Observation / getObservationRecords", &
            "TRACE BACK (1)", 1)
       RETURN       
    END IF
    ra = getRA(this)
    IF (error) THEN
       CALL errorMessage("Observation / getObservationRecords", &
            "TRACE BACK (2)", 1)
       RETURN       
    END IF
    dec = getDec(this)
    IF (error) THEN
       CALL errorMessage("Observation / getObservationRecords", &
            "TRACE BACK (3)", 1)
       RETURN       
    END IF
    obsy_code = getCode(this%obsy)
    IF (error) THEN
       CALL errorMessage("Observation / getObservationRecords", &
            "TRACE BACK (4)", 1)
       RETURN       
    END IF

    IF (TRIM(frmt) /= "des") THEN

       IF (PRESENT(number)) THEN
          number_ = number
       ELSE
          number_ = this%number
       END IF
       IF (this%discovery) THEN
          discovery = "*"
       ELSE
          discovery = " "
       END IF
       ! The timescale used by MPC prior to year 1972 is UT1 and since
       ! then UTC has been used. The maximum difference between UT1 and 
       ! UTC is 0.9 seconds.
       CALL getCalendarDate(t, "TT", year, month, day)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "TRACE BACK (5)", 1)
          RETURN       
       END IF
       IF (year >= 1972) THEN
          CALL getCalendarDate(t, "UTC", year, month, day)
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
       ELSE
          CALL getCalendarDate(t, "UT1", year, month, day)
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "TRACE BACK (15)", 1)
             RETURN       
          END IF
       END IF
       CALL toString(month, month_str, error)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not convert month to string.", 1)
          WRITE(stderr,*) year, month, day, getMJD(t,"utc")
          RETURN
       END IF
       IF (month < 10) THEN
          month_str = "0" // TRIM(month_str)
       END IF
       CALL radiansToHMS(ra, h, m, s)
       IF (ABS(60.0_bp-s) < 0.00001) THEN
          m = m + 1
          s = 0.0_bp
       END IF
       IF (m == 60) THEN
          m = 0
          h = h + 1
       END IF
       IF (h == 24) THEN
          h = 0
       END IF
       CALL radiansToDAMAS(ABS(dec), deg, am, as)
       IF (ABS(60.0_bp-as) < 0.00001) THEN
          am = am + 1
          as = 0.0_bp
       END IF
       IF (am == 60) THEN
          am = 0
          deg = deg + 1
       END IF
       IF (dec >= 0.0_bp) THEN
          sign_dec = "+"
       ELSE
          sign_dec = "-"
       END IF
       CALL toString(h, h_str, error)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (5).", 1)
          RETURN       
       END IF
       IF (h < 10) h_str = "0" // TRIM(h_str)
       CALL toString(m, m_str, error)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (10).", 1)
          RETURN       
       END IF
       IF (m < 10) m_str = "0" // TRIM(m_str)
       CALL toString(s, s_str, error, frmt="(F6.3)")
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (15).", 1)
          RETURN       
       END IF
       IF (s < 10.0_bp) s_str = "0" // TRIM(s_str)
       CALL toString(deg, deg_str, error)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (20).", 1)
          RETURN       
       END IF
       IF (deg < 10) deg_str = "0" // TRIM(deg_str)
       CALL toString(am, am_str, error)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (25).", 1)
          RETURN       
       END IF
       IF (am < 10) am_str = "0" // TRIM(am_str)
       CALL toString(as, as_str, error, frmt="(F5.2)")
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (30).", 1)
          RETURN       
       END IF
       IF (as < 10.0_bp) as_str = "0" // TRIM(as_str)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Conversion error (35).", 1)
          RETURN       
       END IF
       IF (this%mag < 99.0_bp) THEN
          CALL toString(this%mag, mag_str, error, frmt="(F6.3)")
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Conversion error (40).", 1)
             RETURN
          END IF
          IF (this%mag < 10.0_bp) mag_str = " " // TRIM(mag_str)
       ELSE
          mag_str = " "
       END IF

    END IF

    ALLOCATE(records(5))
    records = " "
    SELECT CASE (TRIM(frmt))

    CASE ("mpc","sor")

       IF (number_ == 0) THEN
          number_str(1:5) = "     "
       ELSE
          CALL toString(number_, number_str, error)
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert number to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          DO WHILE (LEN_TRIM(number_str) < 5)
             number_str = "0" // TRIM(number_str)
          END DO
       END IF
       n = LEN_TRIM(this%designation)
       IF (n <= 7) THEN
          designation = this%designation
       ELSE
          designation = this%designation
          CALL MPC3DesToMPCDes(designation)
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "TRACE BACK (20)", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
       END IF
       day_integer = NINT(day*10_ihp**6,kind=ihp)
       day = day_integer/10.0_bp**6
       CALL toString(day, day_str, error, frmt="(F12.9)")
       IF (day < 10.0_bp) day_str = "0" // TRIM(day_str)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not convert day to string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       err = 0
       WRITE(records(1),"(A5,A7,3A1,I4,1X,A2,1X,A8,1X,A2,1X,A2,1X,A6," // &
            "A1,A2,1X,A2,1X,A5,9X,A5,A2,5X,A3)", iostat=err) number_str(1:5), &
            designation(1:7), discovery, this%note1, this%note2, year, &
            month_str, day_str(1:8), h_str, m_str, s_str, sign_dec, &
            deg_str, am_str, as_str, mag_str(1:5), this%filter, obsy_code
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not write output string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       SELECT CASE (this%note2)
       CASE ("S") ! Satellite-based observation
          ccoord = copy(this%satellite_ccoord)
          CALL rotateToEquatorial(ccoord)
          coord = getCoordinates(ccoord)
          IF (coord(1) >= 0.0_bp) THEN
             sign_x = "+"
          ELSE
             sign_x = "-"
          END IF
          IF (coord(2) >= 0.0_bp) THEN
             sign_y = "+"
          ELSE
             sign_y = "-"
          END IF
          IF (coord(3) >= 0.0_bp) THEN
             sign_z = "+"
          ELSE
             sign_z = "-"
          END IF
          SELECT CASE (this%coord_unit)
          CASE (1)
             coord(1:3) = coord(1:3)*km_au
             coord = ABS(coord)
             WRITE(records(2),"(A5,A7,3A1,I4,1X,A2,1X,A8,1X,I1,3(1X,A1,F10.4)," // &
                  "8X,A3)", iostat=err) number_str(1:5), designation(1:7), " ", &
                  this%note1, "s", year, month_str, day_str(1:8), &
                  this%coord_unit, sign_x, coord(1), sign_y, coord(2), &
                  sign_z, coord(3), obsy_code
          CASE (2)
             coord = ABS(coord)
             WRITE(records(2), &
                  "(A5,A7,3A1,I4,1X,A2,1X,A8,1X,I1,3(1X,A1,F10.8),8X,A3)",iostat=err) &
                  number_str(1:5), designation(1:7), " ", &
                  this%note1, "s", year, month_str, day_str(1:8), &
                  this%coord_unit, sign_x, coord(1), sign_y, coord(2), sign_z, coord(3), obsy_code
          CASE default
             error = .TRUE.
             CALL errorMessage("Observation / getObservationRecords", &
                  "Unit type not available.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END SELECT
          records => reallocate(records,2)
       CASE default
          records => reallocate(records,1)
       END SELECT
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not write output string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF

    CASE ("mpc2")

       DO WHILE (LEN_TRIM(obsy_code) < 4)
          obsy_code = "0" // TRIM(obsy_code)
       END DO

       IF (number_ == 0) THEN
          number_str(1:7) = "       "
       ELSE
          CALL toString(number_, number_str, error)
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert number to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          DO WHILE (LEN_TRIM(number_str) < 7)
             number_str = "0" // TRIM(number_str)
          END DO
       END IF
       n = LEN_TRIM(this%designation)
       IF (n <= 9) THEN
          designation = this%designation
       ELSE
          designation = this%designation(n-8:n)
       END IF
       day_integer = NINT(day*10_ihp**9_ihp,kind=ihp)
       day = day_integer/10.0_bp**9.0_bp
       CALL toString(day, day_str, error, frmt="(F12.9)")
       IF (day < 10.0_bp) day_str = "0" // TRIM(day_str)
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not convert day to string (10).", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       IF (SQRT(this%covariance(2,2))*SQRT(this%covariance(3,3)) > EPSILON(correlation)) THEN
          correlation = this%covariance(2,3) / &
               (SQRT(this%covariance(2,2))*SQRT(this%covariance(3,3)))
       ELSE
          correlation = 0.0_bp
       END IF
       ! Optical:
       IF (SQRT(this%covariance(2,2))/rad_asec > 0.01_bp .AND. &
            SQRT(this%covariance(3,3))/rad_asec > 0.01_bp) THEN
          WRITE(records(1),"(A7,A9,3A1,I4,1X,A2,1X,A12,1X,A2,1X,A2,1X,A6," // &
               "3X,A1,A2,1X,A2,1X,A5,4X,A6,A2,A1,1X,A3,A1,A12,1X,A12," // & 
               "1X,A12,A1,2A4)", iostat=err) number_str(1:7), designation(1:9), "1", &
               discovery, this%note2, year, month_str, day_str, h_str, m_str, &
               s_str, sign_dec, deg_str, am_str, as_str, mag_str, this%filter, &
               " ", "   ", " ", "            ", "            ", "            ", &
               " ", "    ", obsy_code

          ! Optical (assumes zero-correlation and error estimates):
          WRITE(records(2),"(A7,A9,3A1,I4,1X,A2,1X,A12,1X,F6.3,1X,F6.3,3X," // &
               "F8.5,A1,A7,F12.10,F5.3,A4,A4,F5.1,I1,3X,A1,2A8,A1,2A4)", iostat=err) &
               number_str(1:7), designation(1:9), "2", " ", this%note2, year, month_str, day_str, &
               SQRT(this%covariance(2,2))/rad_asec, &
               SQRT(this%covariance(3,3))/rad_asec, &
               correlation, "X", "       ", 0.0_bp, 0.0_bp, "    ", "    ", 0.0_bp, 0, &
               " ", "        ", "        ", " ", "    ", obsy_code
       ELSE
          WRITE(records(1),"(A7,A9,3A1,I4,1X,A2,1X,A12,1X,F14.12,1X," // &
               "A1,F14.12,1X,A6,A2,A1,1X,A3,A1,A12,1X,A12," // & 
               "1X,A12,A1,2A4)", iostat=err) number_str(1:7), designation(1:9), "1", &
               discovery, this%note2, year, month_str, day_str, ra,  &
               sign_dec, ABS(dec), mag_str, this%filter, &
               " ", "   ", " ", "            ", "            ", "            ", &
               " ", "    ", obsy_code
          CALL toString(SQRT(this%covariance(2,2))*10000000, ra_unc, error, frmt="(F9.7)")
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert RA uncertainty to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          ra_unc(2:8) = ra_unc(3:9)
          IF (ra_unc == "00000000") THEN
             ra_unc = " "
             ra_unc(1:1) = "0"
          ELSE
             DO i=1,8
                IF (ra_unc(i:i) /= "0") THEN
                   indx = i
                   EXIT
                END IF
             END DO
             ra_unc(indx+2:8) = " "
          END IF
          CALL toString(SQRT(this%covariance(3,3))*10000000, dec_unc, error, frmt="(F9.7)")
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert Dec uncertainty to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          dec_unc(2:8) = dec_unc(3:9)
          IF (dec_unc == "00000000") THEN
             dec_unc = " "
             dec_unc(1:1) = "0"
          ELSE
             DO i=1,8
                IF (dec_unc(i:i) /= "0") THEN
                   indx = i
                   EXIT
                END IF
             END DO
             dec_unc(indx+2:8) = " "
          END IF
          ! Optical (assumes zero-correlation and error estimates):
          WRITE(records(2),"(A7,A9,3A1,I4,1X,A2,1X,A12,1X,A8,A8," // &
               "F8.5,A1,A7,F12.10,F5.3,A4,A4,F5.1,I1,3X,A1,2A8,A1,2A4)", iostat=err) &
               number_str(1:7), designation(1:9), "2", " ", this%note2, &
               year, month_str, day_str, ra_unc(1:8), dec_unc(1:8), correlation, "X", &
               "       ", 0.0_bp, 0.0_bp, "    ", "    ", 0.0_bp, 0, &
               " ", "        ", "        ", " ", "    ", obsy_code
       END IF

       SELECT CASE (this%note2)
       CASE ("S") ! Satellite-based observation
          records(3) = records(1)
          records(3)(17:17) = "3"
          ccoord = copy(this%satellite_ccoord)
          CALL rotateToEquatorial(ccoord)
          coord = getCoordinates(ccoord)
          SELECT CASE (this%coord_unit)
          CASE (1)
             ! in km:
             coord(1:3) = coord(1:3)*km_au
             WRITE(records(3)(41:123),"(A,3(F18.12))", iostat=err) "0", coord(1:3)
          CASE (2)
             ! in AU:
             WRITE(records(3)(41:123),"(A,3(F18.12))", iostat=err) "1", coord(1:3)
          CASE default
             error = .TRUE.
             CALL errorMessage("Observation / getObservationRecords", &
                  "Unit type not available.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END SELECT
          records(3)(18:18) = "+"
          records => reallocate(records,3)
       CASE default
          records(2)(18:18) = "+"
          records => reallocate(records,2)
       END SELECT
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not write output string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF

    CASE ("mpc3")

       DO WHILE (LEN_TRIM(obsy_code) < 4)
          obsy_code = "0" // TRIM(obsy_code)
       END DO

       IF (number_ == 0) THEN
          number_str(1:7) = "       "
       ELSE
          CALL toString(number_, number_str, error)
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert number to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          DO WHILE (LEN_TRIM(number_str) < 7)
             number_str = "0" // TRIM(number_str)
          END DO
       END IF
       designation = " "
       n = LEN_TRIM(this%designation)
       IF (simulated_observations) THEN
          designation = this%designation
       ELSE
          IF (n == 9) THEN
             designation = this%designation
          ELSE IF (n < 9 .AND. n > 0) THEN
             designation = this%designation
             CALL MPCDesToMPC3Des(designation)
          ELSE IF (n > 9) THEN
             designation = this%designation(n-8:n)
             number_str(1:6) = number_str(2:7)
             number_str(7:7) = this%designation(n-9:n-9)
          END IF
       END IF
       day_integer = NINT(day*10_ihp**9_ihp,kind=ihp)
       day = day_integer/10.0_bp**9.0_bp
       CALL toString(day, day_str, error, frmt="(F12.9)")
       IF (error) THEN
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not convert day to string (10).", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       IF (day < 10.0_bp) THEN
          day_str = "0" // TRIM(day_str)
       END IF
       IF (SQRT(this%covariance(2,2))*SQRT(this%covariance(3,3)) > EPSILON(correlation)) THEN
          correlation = this%covariance(2,3) / &
               (SQRT(this%covariance(2,2))*SQRT(this%covariance(3,3)))
       ELSE
          correlation = 0.0_bp
       END IF
       ! Optical:
       IF (SQRT(this%covariance(2,2))/rad_asec > 0.01_bp .AND. &
            SQRT(this%covariance(3,3))/rad_asec > 0.01_bp) THEN
          WRITE(records(1),"(A7,A9,2X,3A1,I4,A2,A12,1X,A2,1X,A2,1X,A6," // &
               "3X,A1,A2,1X,A2,1X,A5,4X,A6,A2,A1,1X,A3,A1,A12,1X,A12," // & 
               "1X,A12,A1,2A4)", iostat=err) number_str(1:7), designation(1:9), "1", &
               discovery, this%note2, year, month_str, day_str, h_str, m_str, &
               s_str, sign_dec, deg_str, am_str, as_str, mag_str, this%filter, &
               " ", "   ", " ", "            ", "            ", "            ", &
               " ", "    ", obsy_code

          ! Optical (assumes zero-correlation and error estimates):
          WRITE(records(2),"(A7,A9,2X,3A1,I4,A2,A12,1X,F6.3,1X,F6.3,3X," // &
               "F8.5,A1,A7,F12.10,F5.3,A4,A4,F5.1,I1,3X,A1,2A8,A1,2A4)", iostat=err) &
               number_str(1:7), designation(1:9), "2", " ", this%note2, &
               year, month_str, day_str, &
               SQRT(this%covariance(2,2))/rad_asec, &
               SQRT(this%covariance(3,3))/rad_asec, &
               correlation, "X", "       ", 0.0_bp, 0.0_bp, "    ", "    ", 0.0_bp, 0, &
               " ", "        ", "        ", " ", "    ", obsy_code
       ELSE
          WRITE(records(1),"(A7,A9,2X,3A1,I4,A2,A12,1X,F14.12,1X," // &
               "A1,F14.12,1X,A6,A2,A1,1X,A3,A1,A12,1X,A12," // & 
               "1X,A12,A1,2A4)", iostat=err) number_str(1:7), &
               designation(1:9), "1", discovery, this%note2, year, month_str, &
               day_str, ra, sign_dec, ABS(dec), mag_str, this%filter, &
               " ", "   ", " ", "            ", "            ", &
               "            ", " ", "    ", obsy_code
          CALL toString(SQRT(this%covariance(2,2))*10000000, ra_unc, error, frmt="(F9.7)")
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert RA uncertainty to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          ra_unc(2:8) = ra_unc(3:9)
          IF (ra_unc == "00000000") THEN
             ra_unc = " "
             ra_unc(1:1) = "0"
          ELSE
             DO i=1,8
                IF (ra_unc(i:i) /= "0") THEN
                   indx = i
                   EXIT
                END IF
             END DO
             ra_unc(indx+2:8) = " "
          END IF
          CALL toString(SQRT(this%covariance(3,3))*10000000, dec_unc, error, frmt="(F9.7)")
          IF (error) THEN
             CALL errorMessage("Observation / getObservationRecords", &
                  "Could not convert Dec uncertainty to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
          dec_unc(2:8) = dec_unc(3:9)
          IF (dec_unc == "00000000") THEN
             dec_unc = " "
             dec_unc(1:1) = "0"
          ELSE
             DO i=1,8
                IF (dec_unc(i:i) /= "0") THEN
                   indx = i
                   EXIT
                END IF
             END DO
             dec_unc(indx+2:8) = " "
          END IF
          ! Optical (assumes zero-correlation and error estimates):
          WRITE(records(2),"(A7,A9,2X,3A1,I4,A2,A12,1X,A8,A8," // &
               "F8.5,A1,A7,F12.10,F5.3,A4,A4,F5.1,I1,3X,A1,2A8,A1,2A4)", &
               iostat=err) number_str(1:7), designation(1:9), "2", &
               " ", this%note2, year, month_str, day_str, ra_unc(1:8), &
               dec_unc(1:8), correlation, "X", "       ", 0.0_bp, 0.0_bp, &
               "    ", "    ", 0.0_bp, 0, " ", "        ", "        ", &
               " ", "    ", obsy_code
       END IF

       SELECT CASE (this%note2)
       CASE ("S") ! Satellite-based observation
          records(3) = records(1)
          records(3)(19:19) = "3"
          ccoord = copy(this%satellite_ccoord)
          CALL rotateToEquatorial(ccoord)
          coord = getCoordinates(ccoord)
          SELECT CASE (this%coord_unit)
          CASE (1)
             ! in km:
             coord(1:3) = coord(1:3)*km_au
             WRITE(records(3)(41:123),"(A,3(F18.12))", iostat=err) "0", coord(1:3)
          CASE (2)
             ! in AU:
             WRITE(records(3)(41:123),"(A,3(F18.12))", iostat=err) "1", coord(1:3)
          CASE default
             error = .TRUE.
             CALL errorMessage("Observation / getObservationRecords", &
                  "Unit type not available.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END SELECT
          records(3)(20:20) = "+"
          records => reallocate(records,3)
       CASE default
          records(2)(20:20) = "+"
          records => reallocate(records,2)
       END SELECT
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / getObservationRecords", &
               "Could not write output string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF

    CASE ("elgb")

       IF (number_ == 0) THEN
          n = LEN_TRIM(this%designation)
          IF (n <= 8) THEN
             designation = this%designation
          ELSE
             designation = this%designation(n-7:n)
          END IF
       ELSE
          CALL toString(number_, designation, error)
          IF (error) THEN
             CALL errorMessage("Observation / getELGBString", &
                  "Could not convert number to string.", 1)
             DEALLOCATE(records, stat=err)
             RETURN
          END IF
       END IF

       err = 0
       WRITE(records(1),"(I4,1X,A2,2X,A8,2X,A2,1X,A2,1X,A6,2X,A1,A2," // &
            "1X,A2,1X,A5,2X,A4,A1,1X,A8,A1,1X,A3,2X,A1,A6)",iostat=err) &
            year, month_str, day_str(1:8), h_str, m_str, s_str, sign_dec, &
            deg_str, am_str, as_str, mag_str(1:4), TRIM(this%filter), &
            designation(1:8), discovery, obsy_code, this%note2, "      "
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / getELGBString", &
               "Could not write output string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       records => reallocate(records,1)

    CASE ("des")

       filter_ = this%filter
       IF (LEN_TRIM(filter_) == 0) THEN
          filter_ = "X"
       END IF
       IF (LEN_TRIM(this%secret_name) == 0) THEN
          secret_name = "X"
       ELSE
          secret_name = this%secret_name
       END IF
       IF (this%number /= 0) THEN
          WRITE(records(1),"(I0,1X,F16.10,1X,A,3(1X,F14.10)," // &
               "2(1X,A),2(1X,F14.10),1X,F14.10,1X,E14.7,1X,A)", iostat=err) & 
               this%number, getMJD(t, "UTC"), "O", &
               ra/rad_deg, dec/rad_deg, this%mag, filter_, &
               TRIM(obsy_code), SQRT(this%covariance(2,2))/rad_asec, &
               SQRT(this%covariance(3,3))/rad_asec, this%mag_unc, &
               this%s2n, TRIM(secret_name) 
       ELSE IF (LEN_TRIM(this%designation) /= 0) THEN
          WRITE(records(1),"(A,1X,F16.10,1X,A,3(1X,F14.10)," // &
               "2(1X,A),2(1X,F14.10),1X,F14.10,1X,E14.7,1X,A)", iostat=err) & 
               TRIM(this%designation), getMJD(t, "UTC"), "O", &
               ra/rad_deg, dec/rad_deg, this%mag, filter_, &
               TRIM(obsy_code), SQRT(this%covariance(2,2))/rad_asec, &
               SQRT(this%covariance(3,3))/rad_asec, this%mag_unc, &
               this%s2n, TRIM(secret_name) 
       ELSE
          error = .TRUE.
          CALL errorMessage("Observation / getObservationRecords / des", &
               "Number and designation missing for this record.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observation / getObservationRecords / des", &
               "Could not write output string.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF
       records => reallocate(records,1)

    END SELECT

  END FUNCTION getObservationRecords





  !! *Description*:
  !!
  !! Returns the observation mask to be used with the spherical
  !! coordinates. Elements set to true corresponds to coordinates that
  !! have been initialized.
  !!
  !! Returns error.
  !!
  FUNCTION getObservationMask_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    LOGICAL, DIMENSION(6)          :: getObservationMask_Obs

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObservationMask", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getObservationMask_Obs = this%obs_mask

  END FUNCTION getObservationMask_Obs





  !! *Description*:
  !!
  !! Returns the spherical observation coordinates (can be any, or
  !! all, of rho, RA, Dec, drho/dt, dRA/dt, dDec/dt) and epoch of
  !! this observation.
  !!
  !! Returns error.
  !!
  FUNCTION getObservationSCoord(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    TYPE (SphericalCoordinates)    :: getObservationSCoord

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObservationSCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getObservationSCoord = copy(this%obs_scoord)

  END FUNCTION getObservationSCoord





  !! *Description*:
  !!
  !! Returns the Observatory object.
  !!
  !! Returns error.
  !!
  FUNCTION getObservatory_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    TYPE (Observatory)             :: getObservatory_Obs

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObservatory", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getObservatory_Obs = copy(this%obsy)

  END FUNCTION getObservatory_Obs





  !! *Description*:
  !!
  !! Returns the Cartesian heliocentric coordinates of the observatory 
  !! at the epoch when the observation was made.
  !!
  !! Returns error.
  !!
  FUNCTION getObservatoryCCoord_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    TYPE (CartesianCoordinates)    :: getObservatoryCCoord_Obs

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObservatoryCCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getObservatoryCCoord_Obs = copy(this%obsy_ccoord)

  END FUNCTION getObservatoryCCoord_Obs





  !! *Description*:
  !!
  !! Returns the IAU designated observatory code for the observatory
  !! where this particular observation was made.
  !!
  !! Returns error.
  !!
  CHARACTER(len=OBSY_CODE_LEN) FUNCTION getObservatoryCode_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getObservatoryCode", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getObservatoryCode_Obs = getCode(this%obsy)

  END FUNCTION getObservatoryCode_Obs





  !! *Description*:
  !!
  !! Returns Right Ascension [rad].
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getRA_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getRA", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getRA_Obs = getLongitude(this%obs_scoord)
    IF (error) THEN
       CALL errorMessage("Observation / getRA", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END FUNCTION getRA_Obs





  !! *Description*:
  !!
  !! Returns solar elongation [rad].
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getSolarElongation(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this

    TYPE (SphericalCoordinates) :: scoord
    REAL(bp), DIMENSION(6) :: coord1, coord2

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getSolarElongation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    coord1 = getCoordinates(this%obs_scoord)
    scoord = getSCoord(opposite(this%obsy_ccoord))
    IF (getFrame(this%obs_scoord) == "equatorial") THEN
       CALL rotateToEquatorial(scoord)
    ELSE
       CALL rotateToEcliptic(scoord)
    END IF
    coord2 = getCoordinates(scoord)
    CALL NULLIFY(scoord)
    getSolarElongation = angularDistance(coord1(2), coord1(3), &
         coord2(2), coord2(3))

  END FUNCTION getSolarElongation





  !! *Description*:
  !!
  FUNCTION getStandardDeviations_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    REAL(bp), DIMENSION(6)         :: getStandardDeviations_Obs
    INTEGER                        :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getStandardDeviations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1,6
       getStandardDeviations_Obs(i) = SQRT(this%covariance(i,i))
    END DO
    WHERE (.NOT.this%obs_mask)
       getStandardDeviations_Obs = 0.0_bp
    END WHERE

  END FUNCTION getStandardDeviations_Obs




  !! *Description*:
  !!
  !! Returns epoch as a Time object.
  !!
  !! Returns error.
  !!
  FUNCTION getTime_Obs(this)

    IMPLICIT NONE
    TYPE (Observation), INTENT(in) :: this
    TYPE (Time)                    :: getTime_Obs

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / getTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getTime_Obs = getTime(this%obs_scoord)
    IF (error) THEN
       CALL errorMessage("Observation / getTime", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END FUNCTION getTime_Obs





  !! *Description*:
  !!
  SUBROUTINE setCovarianceMatrix(this, covariance)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout)    :: this
    REAL(bp), DIMENSION(6,6), INTENT(in) :: covariance

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / setCovarianceMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    this%covariance = covariance

  END SUBROUTINE setCovarianceMatrix





  !! *Description*:
  !!
  !! Set designation of object. The length of the designation
  !! string must be smaller or equal to the DESIGNATION_LEN parameter.
  !!
  !! Returns error.
  !!
  SUBROUTINE setDesignation_Obs(this, designation)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout) :: this
    CHARACTER(len=*), INTENT(in)      :: designation

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / setDesignation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (LEN_TRIM(designation) > DESIGNATION_LEN) THEN
       error = .TRUE.
       CALL errorMessage("Observations / setDesignation", &
            "New designation too long.", 1)
       RETURN
    END IF

    this%designation = TRIM(designation)

  END SUBROUTINE setDesignation_Obs





  !! *Description*:
  !!
  SUBROUTINE setNumber_Obs(this, number)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout) :: this
    INTEGER, INTENT(in)               :: number

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / setNumber", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    this%number = number

  END SUBROUTINE setNumber_Obs





  !! *Description*:
  !!
  !! Set the spherical observation coordinates (can be any, or
  !! all, of rho, RA, Dec, drho/dt, dRA/dt, dDec/dt) and epoch of
  !! this observation.
  !!
  !! Returns error.
  !!
  SUBROUTINE setObservationSCoord(this, obs_scoord)

    IMPLICIT NONE
    TYPE (Observation), INTENT(inout)       :: this
    TYPE (SphericalCoordinates), INTENT(in) :: obs_scoord

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observation / setObservationSCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    CALL NULLIFY(this%obs_scoord)
    this%obs_scoord = copy(obs_scoord)

  END SUBROUTINE setObservationSCoord





END MODULE Observation_cl





