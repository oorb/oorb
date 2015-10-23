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
!! This class provides a datatype and routines for handling
!! observations. Only optical astrometric observations given in the
!! old and new <a
!! href="http://cfa-www.harvard.edu/iau/info/OpticalObs.html">MPC</a>
!! formats, the Data Exchange Standard (DES) format, the Lowell
!! format, and the preliminary GAIA formats are supported for the time
!! being. Besides the observations file, indirectly this class also
!! needs a file containing codes and coordinates for observatories
!! (usually OBSCODE.dat).
!!  
!! *Example*: 
!!  
!! <pre> 
!! use Base_cl 
!! use Observations_cl 
!! implicit none 
!! type (Observations) :: obss 
!! type (File) :: myfile 
!! .  
!! .  
!! .  
!! call new(myfile,"myfilen.ame") 
!! call open(myfile) 
!! call new(obss,myfile) 
!! if (error) stop 
!! call nullify(obss) 
!! </pre> 
!!  
!! @see StochasticOrbit_class 
!!  
!! @author  MG, JV 
!! @version 2015-10-23
!!  
MODULE Observations_cl

  USE Base_cl
  USE File_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  USE CartesianCoordinates_cl
  USE Observation_cl
  USE Observatories_cl
  USE Orbit_cl

  USE utilities
  USE sort
  USE linal
  !$ use omp_lib

  IMPLICIT NONE
  PRIVATE :: NEW_Obss
  PRIVATE :: NEW_Obss_file
  PRIVATE :: NEW_Obss_obs
  PRIVATE :: NEW_Obss_obs_arr
  PRIVATE :: NULLIFY_Obss
  PRIVATE :: copy_Obss
  PRIVATE :: exist_Obss
  PRIVATE :: addMultinormalDeviates_Obss
  PRIVATE :: addUniformDeviates_Obss
  PRIVATE :: addition_Obss
  PRIVATE :: setDesignation_Obss
  PRIVATE :: clean_Obss
  PRIVATE :: getNumber_Obss
  PRIVATE :: getDeclinations_Obss
  PRIVATE :: setNumber_Obss
  PRIVATE :: sortObservations
  PRIVATE :: getDates_Obss
  PRIVATE :: getMinAndMaxValues_Obss
  PRIVATE :: getObservationMasks_Obss
  PRIVATE :: getObservatoryCodes_Obss
  PRIVATE :: getStandardDeviations_Obss
  !  PRIVATE :: readGaiaFile
  PRIVATE :: reallocate_a1_Obss

  TYPE Observations
     PRIVATE
     !! Observations.
     TYPE (Observation), DIMENSION(:), POINTER :: obs_arr => NULL()
     !! Different objects in this Observations object.
     CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER  :: objects => NULL()
     !! Sorting criteria (epochs as tdt and mjd).
     REAL(bp), DIMENSION(:), POINTER           :: criteria => NULL()  
     !! Index vector which puts the observations in ascending order.
     INTEGER, DIMENSION(:), POINTER            :: ind => NULL()  
     !! Number of observations.
     INTEGER                                   :: nobs = 0
     !! Number of different objects in this Observations object.
     INTEGER                                   :: nobjects = 0
     !! If _true_, then this object has been initialized.
     LOGICAL                                   :: is_initialized = .FALSE.
     !! Place for additional information related to observations,
     !! E.g., notes from observation file: additional mask + obs.pair flags
     CHARACTER(len=4), DIMENSION(:), POINTER   :: obs_note_arr => NULL()
  END TYPE Observations

  !! Initializes an Observations-object.
  INTERFACE NEW
     MODULE PROCEDURE NEW_Obss
     MODULE PROCEDURE NEW_Obss_file
     MODULE PROCEDURE NEW_Obss_obs
     MODULE PROCEDURE NEW_Obss_obs_arr
  END INTERFACE NEW

  !! Nullifies an Observations-object.
  INTERFACE NULLIFY
     MODULE PROCEDURE NULLIFY_Obss
  END INTERFACE NULLIFY

  !! Makes a copy of an Observations-object:
  INTERFACE copy
     MODULE PROCEDURE copy_Obss
  END INTERFACE copy

  !! Returns state of an Observations-object:
  INTERFACE exist
     MODULE PROCEDURE exist_Obss
  END INTERFACE exist

  !! Adds Gaussian noise to the observations. 
  INTERFACE addMultinormalDeviates
     MODULE PROCEDURE addMultinormalDeviates_Obss
  END INTERFACE addMultinormalDeviates

  !! Adds uniform noise to the observations. 
  INTERFACE addUniformDeviates
     MODULE PROCEDURE addUniformDeviates_Obss
  END INTERFACE addUniformDeviates

  !! Returns MJDs in UTC for all observations in this object.
  INTERFACE getDates
     MODULE PROCEDURE getDates_Obss
  END INTERFACE getDates

  !! Returns the declinations of all observations in this object.
  INTERFACE getDeclinations
     MODULE PROCEDURE getDeclinations_Obss
  END INTERFACE getDeclinations

  !! Changes the designation for all observations in this object.
  INTERFACE getDesignation
     MODULE PROCEDURE getDesignation_Obss
  END INTERFACE getDesignation

  !! Changes the designation for all observations in this object.
  INTERFACE setDesignation
     MODULE PROCEDURE setDesignation_Obss
  END INTERFACE setDesignation

  !! Removes multiple observations from this object:
  INTERFACE clean
     MODULE PROCEDURE clean_Obss
  END INTERFACE clean

  INTERFACE getID
     MODULE PROCEDURE getID_Obss
  END INTERFACE getID

  !! Gets the end (read documentation) values of this object. 
  INTERFACE getMinAndMaxValues
     MODULE PROCEDURE getMinAndMaxValues_Obss
  END INTERFACE getMinAndMaxValues

  !! Returns the number of this object assuming that the set contains only one object.. 
  INTERFACE getNumber
     MODULE PROCEDURE getNumber_Obss
  END INTERFACE getNumber

  INTERFACE getObservationMasks
     MODULE PROCEDURE getObservationMasks_Obss
  END INTERFACE getObservationMasks

  INTERFACE getObservatoryCodes
     MODULE PROCEDURE getObservatoryCodes_Obss
  END INTERFACE getObservatoryCodes

  !! Returns an Observations-object containing the combined 
  !! observations from two given Observations-objects.
  INTERFACE OPERATOR (+) 
     MODULE PROCEDURE addition_Obss
  END INTERFACE OPERATOR (+)

  INTERFACE getStandardDeviations
     MODULE PROCEDURE getStandardDeviations_Obss
  END INTERFACE getStandardDeviations

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_a1_Obss
  END INTERFACE reallocate

  INTERFACE setNumber
     MODULE PROCEDURE setNumber_Obss
  END INTERFACE setNumber

CONTAINS




  !! *Description*:
  !!
  !! Initializes a new Observations instance with default values.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object already has been initialized or 
  !!               something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! .
  !! .
  !! .
  !! call new(obss)
  !! </pre>
  SUBROUTINE NEW_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%obs_arr => NULL()
    this%objects => NULL()
    this%criteria => NULL()
    this%ind => NULL()
    this%nobs = 0
    this%nobjects = 0
    this%obs_note_arr => NULL()
    this%is_initialized = .TRUE.

  END SUBROUTINE NEW_Obss





  !! *Description*:
  !!
  !! Initializes a new Observations instance. Reads the given 
  !! observations file, transforms the dates to a Julian dates, 
  !! target coordinates to radians and observatory codes to 
  !! geocentric equatorial coordinates in AUs of the observatories.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object already has been initialized or 
  !!               something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! type (File)         :: myfile
  !! .
  !! .
  !! .
  !! call new(myfile,"myfilen.ame")
  !! call open(myfile)
  !! call new(obss, myfile)
  !! </pre>
  !!
  SUBROUTINE NEW_Obss_file(this, f1, stdev, orb_sim)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)           :: this
    TYPE (File), INTENT(in)                      :: f1
    REAL(bp), DIMENSION(6), OPTIONAL, INTENT(in) :: stdev
    TYPE (Orbit), OPTIONAL, INTENT(inout)        :: orb_sim

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(stdev)) THEN
       IF (PRESENT(orb_sim)) THEN
          CALL readObservationFile(this, f1, stdev, orb_sim=orb_sim)
       ELSE
          CALL readObservationFile(this, f1, stdev)
       END IF
    ELSE
       IF (PRESENT(orb_sim)) THEN
          CALL readObservationFile(this, f1, orb_sim=orb_sim)
       ELSE
          CALL readObservationFile(this, f1)
       END IF
    END IF
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    ! Update number and names of different objects:
    CALL sortObservations(this)
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE NEW_Obss_file





  !! *Description*:
  !!
  !! Initializes a new Observations instance using the given
  !! Observation object.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object already has been initialized or
  !!               something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! type (Observation)  :: myobs
  !! .
  !! .
  !! .
  !! call new(obss, myobs)
  !! </pre>
  !!
  SUBROUTINE NEW_Obss_obs(this, obs)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    TYPE (Observation), INTENT(in)     :: obs
    TYPE (Time)                        :: t
    CHARACTER(len=DESIGNATION_LEN)     :: number
    INTEGER                            :: err


    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%nobs = 1 
    ALLOCATE(this%obs_arr(this%nobs), this%objects(this%nobs), &
         this%criteria(this%nobs), this%ind(this%nobs), &
         this%obs_note_arr(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    this%obs_arr(1) = copy(obs)
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "TRACE BACK (5)", 1)
       CALL NULLIFY(this)
       RETURN
    END IF
    this%nobjects = 1
    number = getNumber(obs)
    IF (LEN_TRIM(number) /= 0) THEN
!!$       DO WHILE (LEN_TRIM(number) < 7)
!!$          number = "0" // TRIM(number)
!!$       END DO
       this%objects(this%nobjects) = TRIM(number)
    ELSE
       this%objects(this%nobjects) = getDesignation(obs)
       IF (error) THEN
          CALL errorMessage("Observations / new", &
               "TRACE BACK (20)", 1)
          CALL NULLIFY(this)
          RETURN
       END IF
    END IF
    t = getTime(obs)
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "TRACE BACK (25)", 1)
       CALL NULLIFY(this)
       RETURN
    END IF
    this%criteria(1) = getMJD(t, "tdt")
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "TRACE BACK (30)", 1)
       CALL NULLIFY(this)
       RETURN
    END IF
    this%ind(1) = 1
    this%obs_note_arr = " "
    this%is_initialized = .TRUE.

  END SUBROUTINE NEW_Obss_obs





  !! *Description*:
  !!
  !! Initializes a new Observations instance using the given
  !! array of Observation objects.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object already has been initialized or
  !!               something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! type (Observation), dimension(7) :: myobs_arr
  !! .
  !! .
  !! .
  !! call new(obss, myobs_arr)
  !! </pre>
  !!
  SUBROUTINE NEW_Obss_obs_arr(this, obs_arr)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    TYPE (Observation), DIMENSION(:), INTENT(in) :: obs_arr
    INTEGER                            :: err, i

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%nobs = SIZE(obs_arr,dim=1)
    ALLOCATE(this%obs_arr(this%nobs), this%objects(this%nobs), &
         this%criteria(this%nobs), this%ind(this%nobs), &
         this%obs_note_arr(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Could not allocate memory.", 1)
       CALL NULLIFY(this)
       RETURN
    END IF
    DO i=1,this%nobs
       this%obs_arr(i) = copy(obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / new", &
               "TRACE BACK (5)", 1)
          CALL NULLIFY(this)
          RETURN
       END IF
    END DO

    ! Update number and names of different objects:
    CALL sortObservations(this, force_full=.TRUE.)
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "TRACE BACK (10)", 1)
       CALL NULLIFY(this)
       RETURN
    END IF
    this%obs_note_arr = " "

    this%is_initialized = .TRUE.

  END SUBROUTINE NEW_Obss_obs_arr





  !! *Description*:
  !!
  !! Nullifies this object, i.e. deallocates arrays and 
  !! sets default values to parameters.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! .
  !! .
  !! .
  !! call nullify(obss)
  !! </pre>
  SUBROUTINE NULLIFY_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    INTEGER :: i, err

    IF (ASSOCIATED(this%obs_arr)) THEN
       DO i=1,SIZE(this%obs_arr,dim=1)
          CALL NULLIFY(this%obs_arr(i))
       END DO
       DEALLOCATE(this%obs_arr, stat=err)
       IF (err /= 0) THEN
          NULLIFY(this%obs_arr)
       END IF
    END IF
    IF (ASSOCIATED(this%ind)) THEN 
       DEALLOCATE(this%ind, stat=err)
       IF (err /= 0) THEN
          NULLIFY(this%ind)
       END IF
    END IF
    IF (ASSOCIATED(this%criteria)) THEN
       DEALLOCATE(this%criteria, stat=err)
       IF (err /= 0) THEN
          NULLIFY(this%criteria)
       END IF
    END IF
    IF (ASSOCIATED(this%objects)) THEN
       DEALLOCATE(this%objects, stat=err)
       IF (err /= 0) THEN
          NULLIFY(this%objects)
       END IF
    END IF
    IF (ASSOCIATED(this%obs_note_arr)) THEN
       DEALLOCATE(this%obs_note_arr, stat=err)
       IF (err /= 0) THEN
          NULLIFY(this%obs_note_arr)
       END IF
    END IF
    this%nobjects =  0
    this%nobs =  0
    this%is_initialized = .FALSE.

  END SUBROUTINE NULLIFY_Obss





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  !! Returns error if allocation fails.
  !!
  FUNCTION copy_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    TYPE (Observations)             :: copy_Obss
    INTEGER                         :: i, err

    copy_Obss%nobs = this%nobs
    copy_Obss%nobjects = this%nobjects
    IF (ASSOCIATED(this%obs_arr)) THEN
       ALLOCATE(copy_Obss%obs_arr(this%nobs), &
            copy_Obss%objects(this%nobjects), &
            copy_Obss%ind(this%nobs), &
            copy_Obss%criteria(this%nobs), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / copy", &
               "Could not allocate memory (5).", 1)
          CALL NULLIFY(copy_Obss)
          RETURN
       END IF
       DO i=1,SIZE(this%obs_arr)
          copy_Obss%obs_arr(i) = copy(this%obs_arr(i))
          IF (error) THEN
             CALL errorMessage("Observations / Observations", &
                  "TRACE BACK", 1)
             CALL NULLIFY(copy_Obss)
             RETURN
          END IF
       END DO
       copy_Obss%objects = this%objects
       copy_Obss%ind = this%ind
       copy_Obss%criteria = this%criteria
    END IF
    IF (ASSOCIATED(this%obs_note_arr)) THEN
       ALLOCATE(copy_Obss%obs_note_arr(this%nobs), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / copy", &
               "Could not allocate memory (10).", 1)
          CALL NULLIFY(copy_Obss)
          RETURN
       END IF
       copy_Obss%obs_note_arr = this%obs_note_arr
    END IF
    copy_Obss%is_initialized = this%is_initialized

  END FUNCTION copy_Obss





  !! *Description*:
  !!
  !! Returns the status of the object, i.e. whether it exists or not.
  !!
  LOGICAL FUNCTION exist_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this

    exist_Obss = this%is_initialized

  END FUNCTION exist_Obss





  !! *Description*:
  !!
  !! Adds Gaussian deviates to the observations, optionally using a
  !! mask. The new center relative to the original coordinates and the
  !! covariance matrices should be given in radians for the angular
  !! space coordinates, AUs for distance, radians per day for angular
  !! velocities and AUs per day for line-of-sight velocity. If mask is
  !! present, apply only for observations for which mask is true.
  !!
  !! Returns error.
  !!
  SUBROUTINE addMultinormalDeviates_Obss(this, mean_arr, covariance_arr, mask_arr)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)          :: this
    REAL(bp), DIMENSION(:,:), INTENT(in)        :: mean_arr
    REAL(bp), DIMENSION(:,:,:), INTENT(in)      :: covariance_arr
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask_arr
    LOGICAL, DIMENSION(:), ALLOCATABLE          :: mask_arr_
    INTEGER                                     :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addMultinormalDeviates", &
            "Object has not been initialized.", 1)
       RETURN       
    END IF

    IF (SIZE(this%obs_arr) /= SIZE(mean_arr,dim=1)) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addMultinormalDeviates", &
            "Number of observations and length of mean_arr are not compatible.", 1)
       RETURN
    END IF
    IF (SIZE(this%obs_arr) /= SIZE(covariance_arr,dim=1)) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addMultinormalDeviates", &
            "Number of observations and length of covariance_arr are not compatible.", 1)
       RETURN
    END IF

    ALLOCATE(mask_arr_(SIZE(this%obs_arr)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addMultinormalDeviate", &
            "Could not allocate memory.", 1)
       DEALLOCATE(mask_arr_, stat=err)
       RETURN
    END IF

    IF (PRESENT(mask_arr)) THEN
       IF (SIZE(this%obs_arr) /= SIZE(mask_arr,dim=1)) THEN
          error = .TRUE.
          CALL errorMessage("Observations / addMultinormalDeviates", &
               "Number of observations and length of mask_arr are not compatible.", 1)
          RETURN
       END IF
       mask_arr_ = mask_arr 
    ELSE 
       mask_arr_(:) = .TRUE.
    END IF

    DO i=1,SIZE(this%obs_arr)
       IF (mask_arr_(i)) THEN
          CALL addMultinormalDeviate(this%obs_arr(this%ind(i)), mean_arr(i,:), covariance_arr(i,:,:))
          IF (error) THEN
             CALL errorMessage("Observations / addMultinormalDeviate", &
                  "TRACE BACK", 1)
             DEALLOCATE(mask_arr_, stat=err)
             RETURN
          END IF
       END IF
    END DO

    DEALLOCATE(mask_arr_, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addMultinormalDeviate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END SUBROUTINE addMultinormalDeviates_Obss





  !! *Description*:
  !!
  !! Adds uniform deviates to the observations (i), optionally using a
  !! mask. The center relative to the original state and the absolute
  !! values of the boundary values ((i,j,1)=center,
  !! (i,j,2)=abs(boundary), position j=1:3 and velocity j=4:6) should
  !! be given in radians for the angular space coordinates, AUs for
  !! distance, radians per day for angular velocities and AUs per day
  !! for line-of-sight velocity.
  !!
  !! Returns error.
  !!
  SUBROUTINE addUniformDeviates_Obss(this, center_and_absbound_arr, mask_arr)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)          :: this
    REAL(bp), DIMENSION(:,:,:), INTENT(in)      :: center_and_absbound_arr
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask_arr
    LOGICAL, DIMENSION(:), ALLOCATABLE          :: mask_arr_
    INTEGER                                     :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addUniformDeviates", &
            "Object has not been initialized.", 1)
       RETURN       
    END IF
    IF (SIZE(this%obs_arr) /= SIZE(center_and_absbound_arr,dim=1)) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addUniformDeviates", &
            "Number of observations and length of center_and_absbound_arr are not compatible.", 1)
       RETURN
    END IF

    ALLOCATE(mask_arr_(SIZE(this%obs_arr)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addUniformDeviates", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    IF (PRESENT(mask_arr)) THEN
       IF (SIZE(this%obs_arr) /= SIZE(mask_arr,dim=1)) THEN
          error = .TRUE.
          CALL errorMessage("Observations / addMultinormalDeviates", &
               "Number of observations and length of mask_arr are not compatible.", 1)
          RETURN
       END IF
       mask_arr_ = mask_arr 
    ELSE 
       mask_arr_(:) = .TRUE.
    END IF

    DO i=1,SIZE(this%obs_arr)
       IF (mask_arr_(i)) THEN
          CALL addUniformDeviate(this%obs_arr(this%ind(i)), center_and_absbound_arr(i,:,:))
          IF (error) THEN
             CALL errorMessage("Observations / addUniformDeviates", &
                  "TRACE BACK", 1)
             DEALLOCATE(mask_arr_, stat=err)
             RETURN       
          END IF
       END IF
    END DO

    DEALLOCATE(mask_arr_, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addUniformDeviates", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END SUBROUTINE addUniformDeviates_Obss





  !! *Description*:
  !!
  !! Returns the sum of two objects, i.e. one object that contains all
  !! observations from both original objects.
  !!
  !! Returns error.
  !!
  FUNCTION addition_Obss(this, that)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    TYPE (Observations), INTENT(in) :: that
    TYPE (Observations)             :: addition_Obss
    INTEGER                         :: err, i

    IF (.NOT. this%is_initialized) THEN
       CALL errorMessage("Observations / + (addition)", &
            "Left object has not been initialized.", 1)
       CALL NULLIFY(addition_Obss)
       RETURN       
    END IF

    IF (.NOT. that%is_initialized) THEN
       CALL errorMessage("Observations / + (addition)", &
            "Right object has not been initialized.", 1)
       CALL NULLIFY(addition_Obss)
       RETURN       
    END IF

    ! Add observations:
    ALLOCATE(addition_Obss%obs_arr(this%nobs+that%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addition", &
            "Could not allocate astrometric array.", 1)
       CALL NULLIFY(addition_Obss)
       RETURN
    END IF
    DO i=1, this%nobs
       addition_Obss%obs_arr(i) = copy(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / addition", &
               "TRACE BACK (5)", 1)
          CALL NULLIFY(addition_Obss)
          RETURN
       END IF
    END DO
    DO i=1, that%nobs
       addition_Obss%obs_arr(this%nobs+i) = copy(that%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / addition", &
               "TRACE BACK (10)", 1)
          CALL NULLIFY(addition_Obss)
          RETURN
       END IF
    END DO

    ! Use number of observations, and number and names 
    ! of different objects corresponding to "this" as such:
    addition_Obss%nobs = this%nobs
    addition_Obss%nobjects = this%nobjects
    ALLOCATE(addition_Obss%objects(this%nobjects), &
         addition_Obss%ind(this%nobs), &
         addition_Obss%criteria(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addition", &
            "Could not allocate memory.", 1)
       CALL NULLIFY(addition_Obss)
       RETURN
    END IF
    addition_Obss%objects  = this%objects
    addition_Obss%ind      = this%ind
    addition_Obss%criteria = this%criteria

    ! Add obs_note_arr if present in either of the objects:
    IF (ASSOCIATED(this%obs_note_arr) .OR. &
         ASSOCIATED(that%obs_note_arr)) THEN
       ALLOCATE(addition_Obss%obs_note_arr(this%nobs+that%nobs), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / addition", &
               "Could not allocate obs_note array.", 1)
          CALL NULLIFY(addition_Obss)
          RETURN
       END IF
       IF (ASSOCIATED(this%obs_note_arr)) THEN
          addition_Obss%obs_note_arr(1:this%nobs) = &
               this%obs_note_arr(1:this%nobs)
       ELSE
          addition_Obss%obs_note_arr(1:this%nobs) = " "
       END IF
       IF (ASSOCIATED(that%obs_note_arr)) THEN
          addition_Obss%obs_note_arr(this%nobs+1:) = &
               that%obs_note_arr(1:that%nobs)
       ELSE
          addition_Obss%obs_note_arr(this%nobs+1:) = " "
       END IF
    END IF

    ! Add number of observations, and number and names 
    ! of different objects corresponding to "that":
    CALL sortObservations(addition_Obss)
    IF (error) THEN
       CALL errorMessage("Observations / addition", &
            "TRACE BACK", 1)
       CALL NULLIFY(addition_Obss)
       RETURN
    END IF

    addition_Obss%is_initialized = .TRUE.

  END FUNCTION addition_Obss





  !! *Description*:
  !!
  !! Adds a single observation to this object.
  !!
  !! Returns error.
  !!
  SUBROUTINE addObservation(this, obs, sort, note)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    TYPE (Observation), INTENT(in)     :: obs
    LOGICAL, INTENT(in), OPTIONAL      :: sort
    CHARACTER(4), INTENT(in), OPTIONAL :: note
    INTEGER :: nobs
    LOGICAL                            :: sort_


    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addObservation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obs)) THEN
       error = .TRUE.
       CALL errorMessage("Observations / addObservation", &
            "Observation to be added has not yet been initialized.", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%obs_arr)) THEN
       nobs = SIZE(this%obs_arr)
    ELSE
       nobs = 0
    END IF
    this%obs_arr => reallocate(this%obs_arr, nobs+1)
    this%obs_arr(nobs+1) = copy(obs)
    this%obs_note_arr => reallocate(this%obs_note_arr, nobs+1)
    IF (PRESENT(note)) THEN
       this%obs_note_arr(nobs+1) = note
    ELSE
       this%obs_note_arr(nobs+1) = " "
    END IF
    ! Update number of observations, and number and names 
    ! of different objects as default behaviour:
    IF (PRESENT(sort)) THEN
       sort_ = sort
    ELSE
       sort_ = .TRUE.
    END IF
    IF (sort_) THEN
       CALL sortObservations(this)
       IF (error) THEN
          CALL errorMessage("Observations / addObservation", &
               "TRACE BACK", 1)
          RETURN
       END IF
    ELSE
       this%nobs = nobs + 1
       this%nobjects = 0
    END IF

  END SUBROUTINE addObservation





  !! *Description*:
  !!
  !! Set designation for all observations within this object.
  !!
  !! Returns error.
  !!
  SUBROUTINE setDesignation_Obss(this, designation)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    CHARACTER(len=*), INTENT(in)       :: designation
    INTEGER                            :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / setDesignation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1, SIZE(this%obs_arr)
       CALL setDesignation(this%obs_arr(i), TRIM(designation))
       IF (error) THEN
          CALL errorMessage("Observations / setDesignation", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    END DO

    ! Update number and names of different objects:
    this%nobjects = 0
    CALL sortObservations(this)
    IF (error) THEN
       CALL errorMessage("Observations / setDesignation", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

  END SUBROUTINE setDesignation_Obss





  !! *Description*:
  !!
  !! Removes multiple instances of the identical observations,
  !! reallocates the observation arrays, and sorts the data according
  !! to (1) designations and (2) numbers..
  !!
  !! Returns error.
  !!
  SUBROUTINE clean_Obss(this)

    TYPE (Observations), INTENT(inout) :: this
    INTEGER :: i, j

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / clean", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    clean1:DO i=1,SIZE(this%obs_arr,dim=1)
       IF (.NOT.exist(this%obs_arr(this%ind(i)))) CYCLE clean1
       clean2:DO j=i+1,SIZE(this%obs_arr,dim=1)
          IF (.NOT.exist(this%obs_arr(this%ind(j)))) CYCLE clean2
          IF (equal(this%obs_arr(this%ind(i)),this%obs_arr(this%ind(j)))) THEN
             CALL NULLIFY(this%obs_arr(this%ind(j)))
          END IF
          IF (error) THEN
             CALL errorMessage("Observations / clean", &
                  "TRACE BACK (5)", 1)
             RETURN
          END IF
       END DO clean2
    END DO clean1
    this%obs_arr => reallocate(this%obs_arr)

    ! Update number and names of different objects:
    this%nobjects = 0
    CALL sortObservations(this,primary_sort="designation")
    IF (error) THEN
       CALL errorMessage("Observations / clean", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

  END SUBROUTINE clean_Obss





  !! *Description*:
  !!
  !! Returns true, if there are different coordinates for simultaneous
  !! observations from the same observatory code.
  !!
  LOGICAL FUNCTION containsInconsistencies(this)

    TYPE (Observations), INTENT(in) :: this
    TYPE (SphericalCoordinates), DIMENSION(:), ALLOCATABLE :: scoord_arr
    TYPE (Time), DIMENSION(:), ALLOCATABLE :: t_arr
    TYPE (Observatory), DIMENSION(:), ALLOCATABLE :: obsy_arr
    INTEGER :: i, j, nobs, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / containsInconsistencies", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    nobs = SIZE(this%obs_arr,dim=1)
    ALLOCATE(scoord_arr(nobs), t_arr(nobs), obsy_arr(nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / containsInconsistencies", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    DO i=1,nobs
       scoord_arr(i) = getObservationSCoord(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / containsInconsistencies", &
               "TTRACE BACK (5)", 1)
          DEALLOCATE(scoord_arr, stat=err)
          DEALLOCATE(t_arr, stat=err)
          DEALLOCATE(obsy_arr, stat=err)
          RETURN
       END IF
       t_arr(i) = getTime(scoord_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / containsInconsistencies", &
               "TTRACE BACK (10)", 1)
          DEALLOCATE(scoord_arr, stat=err)
          DEALLOCATE(t_arr, stat=err)
          DEALLOCATE(obsy_arr, stat=err)
          RETURN
       END IF
       obsy_arr(i) = getObservatory(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / containsInconsistencies", &
               "TTRACE BACK (15)", 1)
          DEALLOCATE(scoord_arr, stat=err)
          DEALLOCATE(t_arr, stat=err)
          DEALLOCATE(obsy_arr, stat=err)
          RETURN
       END IF
    END DO

    containsInconsistencies = .FALSE.
    ci:DO i=1,nobs
       DO j=i+1,nobs
          IF (.NOT.equal(scoord_arr(i), scoord_arr(j)) .AND. &
               equal(t_arr(i), t_arr(j)) .AND. &
               equal(obsy_arr(i), obsy_arr(j))) THEN
             containsInconsistencies = .TRUE.
             EXIT ci
          END IF
          IF (error) THEN
             CALL errorMessage("Observations / containsInconsistencies", &
                  "TTRACE BACK (20)", 1)
             DEALLOCATE(scoord_arr, stat=err)
             DEALLOCATE(t_arr, stat=err)
             DEALLOCATE(obsy_arr, stat=err)
             RETURN
          END IF
       END DO
    END DO ci
    DEALLOCATE(scoord_arr, t_arr, obsy_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(scoord_arr, stat=err)
       DEALLOCATE(t_arr, stat=err)
       DEALLOCATE(obsy_arr, stat=err)
       CALL errorMessage("Observations / containsInconsistencies", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

  END FUNCTION containsInconsistencies





  !! *Description*:
  !!
  !! Returns the blocks of the block diagonal information matrix
  !! containing all observations in this object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! block_diagonal_information_matrix =>
  !! getBlockDiagonalInformationMatrix(myobservations)
  !!
  FUNCTION getBlockDiagInformationMatrix(this)

    TYPE (Observations), INTENT(in)     :: this
    REAL(bp), DIMENSION(:,:,:), POINTER :: getBlockDiagInformationMatrix

    REAL(bp), DIMENSION(:,:,:), POINTER :: covariance_matrices => NULL()
    INTEGER :: i, j, err
    LOGICAL, DIMENSION(6) :: obs_mask

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getBlockDiagInformationMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getBlockDiagInformationMatrix", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getBlockDiagInformationMatrix(this%nobs,6,6), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getBlockDiagInformationMatrix", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    covariance_matrices => getCovarianceMatrices(this)
    IF (error) THEN
       CALL errorMessage("Observations / getBlockDiagInformationMatrix", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(covariance_matrices, stat=err)
       DEALLOCATE(getBlockDiagInformationMatrix, stat=err)
       RETURN
    END IF

    getBlockDiagInformationMatrix = 0.0_bp
    DO i=1,this%nobs
       ! The covariance matrix can be inverted separately for
       ! each observation by using the properties of block diagonal
       ! matrices:
       obs_mask = getObservationMask(this%obs_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getBlockDiagInformationMatrix", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(covariance_matrices, stat=err)
          DEALLOCATE(getBlockDiagInformationMatrix, stat=err)
          RETURN
       END IF
       DO j=1,6
          IF (.NOT.obs_mask(j)) THEN
             covariance_matrices(i,j,j) = 1.0_bp
          END IF
       END DO
       getBlockDiagInformationMatrix(i,1:6,1:6) = &
            matinv(covariance_matrices(i,:,:), errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / " // &
               "getBlockDiagInformationMatrix", &
               "Could not invert covariance matrix for observations " // &
               TRIM(errstr), 1)
          DEALLOCATE(covariance_matrices, stat=err)
          DEALLOCATE(getBlockDiagInformationMatrix, stat=err)
          RETURN
       END IF
       DO j=1,6
          IF (.NOT.obs_mask(j)) THEN
             getBlockDiagInformationMatrix(i,j,j) = 0.0_bp
          END IF
       END DO
    END DO
    DEALLOCATE(covariance_matrices, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getBlockDiagInformationMatrix", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF


  END FUNCTION getBlockDiagInformationMatrix





  !! *Description*:
  !!
  !! Returns an array containing covariance matrices for all
  !! observations in this object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! covariance_matrices => getCovarianceMatrices(myobservations)
  !!
  FUNCTION getCovarianceMatrices(this)

    TYPE (Observations), INTENT(in)     :: this
    REAL(bp), DIMENSION(:,:,:), POINTER :: getCovarianceMatrices

    INTEGER                             :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getCovarianceMatrices", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getCovarianceMatrices", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getCovarianceMatrices(this%nobs,6,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getCovarianceMatrices", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getCovarianceMatrices(i,:,:) = getCovarianceMatrix(this%obs_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getCovarianceMatrices", &
               "TRACE BACK", 1)
          DEALLOCATE(getCovarianceMatrices, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getCovarianceMatrices





  !! *Description*:
  !!
  !! Returns MJDs in UTC for all observations within this object.
  !!
  !! Returns error.
  !!
  FUNCTION getDates_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER :: getDates_Obss

    TYPE (Time)                     :: t
    INTEGER                         :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDates", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDates", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getDates_Obss(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDates", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1,this%nobs
       t = getTime(this%obs_arr(i))
       getDates_Obss(i) = getMJD(t, "UTC")
       IF (error) THEN
          CALL errorMessage("Observations / getDates", &
               "TRACE BACK", 1)
          DEALLOCATE(getDates_Obss, stat=err)
          RETURN
       END IF
       CALL NULLIFY(t)
    END DO

  END FUNCTION getDates_Obss





  !! *Description*:
  !!
  !! Returns declinations for all observations within this object.
  !!
  !! Returns error.
  !!
  FUNCTION getDeclinations_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER :: getDeclinations_Obss

    INTEGER                         :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDeclinations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDeclinations", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getDeclinations_Obss(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDeclinations", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1,this%nobs
       getDeclinations_Obss(i) = getDec(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / getDeclinations", &
               "TRACE BACK", 1)
          DEALLOCATE(getDeclinations_Obss, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getDeclinations_Obss





  !! *Description*:
  !!
  !! Returns designations for all observations within this object.
  !!
  !! Returns error.
  !!
  FUNCTION getDesignations(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)                       :: this
    CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: getDesignations

    INTEGER                                               :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDesignations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDesignations", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getDesignations(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDesignations", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1,this%nobs
       getDesignations(i) = getDesignation(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / getDesignations", &
               "TRACE BACK", 1)
          DEALLOCATE(getDesignations, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getDesignations





  !! *Description*:
  !!
  !! Returns the _only_ designation for all observations within this
  !! object. If more than one designations exist, an error occurs.
  !!
  !! Returns error.
  !!
  FUNCTION getDesignation_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    CHARACTER(len=DESIGNATION_LEN)  :: getDesignation_Obss

    CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: &
         designations => NULL()
    INTEGER                                               :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDesignation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getDesignation", &
            "Observations missing.", 1)
       RETURN
    END IF

    designations => getDesignations(this)
    IF (error) THEN
       CALL errorMessage("Observations / getDesignation", &
            "TRACE BACK", 1)
       RETURN
    END IF

    DO i=2,this%nobs
       IF (designations(1) /= designations(i)) THEN
          error = .TRUE.
          CALL errorMessage("Observations / getDesignation", &
               "More than one designation.", 1)
          DEALLOCATE(designations, stat=err)
          RETURN
       END IF
    END DO
    getDesignation_Obss = designations(1)
    DEALLOCATE(designations, stat=err)

  END FUNCTION getDesignation_Obss





  !! *Description*:
  !!
  !! Returns filters for observations.
  !!
  !! *Usage*:
  !!
  !! filters => getFilters(myobservations)
  !!
  !! Returns error.
  !!
  FUNCTION getFilters(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)         :: this
    CHARACTER(len=2), DIMENSION(:), POINTER :: getFilters

    INTEGER                                 :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getFilters", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getFilters", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getFilters(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getFilters", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getFilters(i) = getFilter(this%obs_arr(this%ind(i))) 
       IF (error) THEN
          CALL errorMessage("Observations / getFilters", &
               "TRACE BACK", 1)
          DEALLOCATE(getFilters, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getFilters





  !! *Description*:
  !!
  !! Returns the ID (number or designation) for an single object object.
  !!
  !! Returns error.
  !!
  FUNCTION getID_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    CHARACTER(len=DESIGNATION_LEN)  :: getID_Obss

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getID", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1 .OR. .NOT.ASSOCIATED(this%objects)) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getID", &
            "Observations missing.", 1)
       RETURN
    END IF

    IF (SIZE(this%objects, dim=1) == 1) THEN
       getID_Obss = this%objects(1)
    ELSE
       error = .TRUE.
       CALL errorMessage("Observations / getID", &
            "Contains more than one objects.", 1)
       RETURN
    END IF

  END FUNCTION getID_Obss





  !! *Description*:
  !!
  !! Returns the complete information matrix containing all
  !! observations in this object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! information_matrix => getInformationMatrix(myobservations)
  !!
  FUNCTION getInformationMatrix(this)

    TYPE (Observations), INTENT(in)     :: this
    REAL(bp), DIMENSION(:,:), POINTER   :: getInformationMatrix

    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         covariance_matrices => NULL()
    INTEGER                             :: i, j, k, err
    LOGICAL, DIMENSION(6)               :: obs_mask

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getInformationMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getInformationMatrix", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getInformationMatrix(this%nobs*6,this%nobs*6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getInformationMatrix", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    covariance_matrices => getCovarianceMatrices(this)
    IF (error) THEN
       CALL errorMessage("Observations / getInformationMatrix", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(getInformationMatrix, stat=err)
       DEALLOCATE(covariance_matrices, stat=err)
       RETURN
    END IF
    getInformationMatrix = 0.0_bp
    DO i=1,this%nobs
       ! The covariance matrix can be inverted separately for
       ! each observation by using the properties of block diagonal
       ! matrices:
       j = (i-1)*6
       obs_mask = getObservationMask(this%obs_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getInformationMatrix", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(getInformationMatrix, stat=err)
          DEALLOCATE(covariance_matrices, stat=err)
          RETURN
       END IF
       DO k=1,6
          IF (.NOT.obs_mask(k)) THEN
             covariance_matrices(i,k,k) = 1.0_bp
          END IF
       END DO
       getInformationMatrix(j+1:j+6,j+1:j+6) = &
            matinv(covariance_matrices(i,:,:), errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / " // &
               "getInformationMatrix", &
               "Could not invert covariance matrix for observations " // &
               TRIM(errstr), 1)
          DEALLOCATE(getInformationMatrix, stat=err)
          DEALLOCATE(covariance_matrices, stat=err)
          RETURN
       END IF
       DO k=1,6
          IF (.NOT.obs_mask(k)) THEN
             getInformationMatrix(j+k,j+k) = 0.0_bp
          END IF
       END DO
    END DO
    DEALLOCATE(covariance_matrices, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getInformationMatrix", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF


  END FUNCTION getInformationMatrix





  !! *Description*:
  !!
  !! Returns magnitudes for observations [mag].
  !!
  !! *Usage*:
  !!
  !! mags => getMagnitudes(myobservations)
  !!
  !! Returns error.
  !!
  FUNCTION getMagnitudes(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER :: getMagnitudes

    INTEGER                         :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMagnitudes", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMagnitudes", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getMagnitudes(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMagnitudes", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getMagnitudes(i) = getMagnitude(this%obs_arr(this%ind(i))) 
       IF (error) THEN
          CALL errorMessage("Observations / getMagnitudes", &
               "TRACE BACK", 1)
          DEALLOCATE(getMagnitudes, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getMagnitudes





  !! *Description*:
  !!
  !! Returns magnitude uncertainties for observations [mag].
  !!
  !! *Usage*:
  !!
  !! mag_uncs => getMagnitudeUncertainties(myobservations)
  !!
  !! Returns error.
  !!
  FUNCTION getMagnitudeUncertainties(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER :: getMagnitudeUncertainties

    INTEGER                         :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMagnitudeUncertainties", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMagnitudeUncertainties", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getMagnitudeUncertainties(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMagnitudeUncertainties", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getMagnitudeUncertainties(i) = getMagnitudeUncertainty(this%obs_arr(this%ind(i))) 
       IF (error) THEN
          CALL errorMessage("Observations / getMagnitudeUncertainties", &
               "TRACE BACK", 1)
          DEALLOCATE(getMagnitudeUncertainties, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getMagnitudeUncertainties





  !! *Description*:
  !!
  !! Returns two onedimensional arrays containing the minimum and maximum values
  !! for R.A., Dec. and Julian date (UTC). The values does not necessarily correspond
  !! to the same observation. They are just the minima and maxima of this particular
  !! Observations object.
  !!
  !! *OBS*! _For_ R.A. (index=2) the extremes are given as
  !! a minimum for the rightmost coordinate, and maximum for the leftmost. In other
  !! words: if the object contains R.A. coordinates from 23h to 1h, then 23h will 
  !! be the minimum and 1h the maximum. If there are only 2 coordinates, the minimum
  !! and maximum value will be defined by requiring ABS(&Delta(R.A.)) to be as small
  !! as possible.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object has not been initialized or something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! real(bp), dimension(:,:), pointer :: minmax
  !! .
  !! .
  !! .
  !! minmax => getMinAndMaxValues(obss)
  !! </pre>
  FUNCTION getMinAndMaxValues_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)       :: this
    REAL(bp), DIMENSION(:,:), POINTER     :: getMinAndMaxValues_Obss

    TYPE (Time)                           :: t
    TYPE (CartesianCoordinates)           :: ccoord
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: obs_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE   :: tmp
    REAL(bp)                              :: ra_min1, ra_max1, &
         ra_min2, ra_max2
    INTEGER                               :: i, no, err
    LOGICAL, DIMENSION(:), ALLOCATABLE    :: mask_array

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMinAndMaxValues", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMinAndMaxValues", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getMinAndMaxValues_Obss(3,2), &
         obs_arr(this%nobs,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMinAndMaxValues", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs

       ! Observation date MJD (UTC):
       t = getTime(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / getMinAndMaxValues", &
               "TRACE BACK 1", 1)
          DEALLOCATE(getMinAndMaxValues_Obss, stat=err)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(t)
          RETURN
       END IF
       obs_arr(i,1) = getMJD(t, "utc")
       IF (error) THEN
          CALL errorMessage("Observations / getMinAndMaxValues", &
               "TRACE BACK 2", 1)
          DEALLOCATE(getMinAndMaxValues_Obss, stat=err)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(t)
          RETURN
       END IF

       ! Right Ascension (rad):
       obs_arr(i,2) = getRA(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / getMinAndMaxValues", &
               "TRACE BACK 3", 1)
          DEALLOCATE(getMinAndMaxValues_Obss, stat=err)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(t)
          RETURN
       END IF

       ! Declination (rad):
       obs_arr(i,3) = getDec(this%obs_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / getMinAndMaxValues", &
               "TRACE BACK 4", 1)
          DEALLOCATE(getMinAndMaxValues_Obss, stat=err)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(t)
          RETURN
       END IF

       ! Observatory coordinates (heliocentric equatorial Cartesian coordinates
       ! in AU's) 
       ccoord = getObservatoryCCoord(this%obs_arr(i))
       obs_arr(i,4:6) = getPosition(ccoord)
       CALL NULLIFY(ccoord)
       CALL NULLIFY(t)
       IF (error) THEN
          CALL errorMessage("Observations / getMinAndMaxValues", &
               "TRACE BACK 5", 1)
          DEALLOCATE(getMinAndMaxValues_Obss, stat=err)
          DEALLOCATE(obs_arr, stat=err)
          RETURN
       END IF

    END DO

    ! Find the extreme values:
    ! Julian date (UTC):
    getMinAndMaxValues_Obss(1,1) = MINVAL(obs_arr(:,1))
    getMinAndMaxValues_Obss(1,2) = MAXVAL(obs_arr(:,1))

    ! Right ascension:
    ! If there are more than one observation and R.A. goes from, for 
    ! instance, 23h -> 1h, we have a problem since it seems like 1h
    ! would be the minimum instead of 23h. So we ca not make an easy
    ! substitution, but must make a work-around: 
    ALLOCATE(mask_array(SIZE(obs_arr,1)), &
         tmp(SIZE(obs_arr,1)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getMinAndMaxValues", &
            "Could not allocate memory for temporary arrays.", 1)
       DEALLOCATE(getMinAndMaxValues_Obss, stat=err)
       DEALLOCATE(obs_arr, stat=err)
       DEALLOCATE(tmp, stat=err)
       DEALLOCATE(mask_array, stat=err)
       RETURN
    END IF
    mask_array = .FALSE.
    ! Choose observations with R.A. in 
    ! the interval [0,pi[ ([0h,12h[):
    tmp = obs_arr(:,2) + pi
    WHERE (tmp < two_pi) mask_array = .TRUE.
    ! Count them:
    no = 0
    DO i=1, SIZE(mask_array)
       IF (mask_array(i)) no = no + 1
    END DO
    ! If there are observations, find the extremes in 
    ! the interval [0,pi[ ([0h,12h[):
    IF (no > 0) THEN
       ra_min1 = -pi + MINVAL(tmp, mask=mask_array)
       ra_max1 = -pi + MAXVAL(tmp, mask=mask_array)
    ELSE ! Else put a mark:
       ra_min1 = -1.0_bp
    END IF
    mask_array = .FALSE.
    ! Choose observations with R.A. in 
    ! the interval [pi,2pi[ ([12h,24h[):
    tmp = obs_arr(:,2) - pi
    WHERE (tmp >= 0.0_bp) mask_array = .TRUE.
    ! Count them:
    no = 0
    DO i=1, SIZE(mask_array)
       IF (mask_array(i)) no = no + 1
    END DO
    ! If there are observations, find the extremes in 
    ! the interval [pi,2pi[ ([12h,24h[):
    IF (no > 0) THEN
       ra_min2 = pi + MINVAL(tmp, mask=mask_array)
       ra_max2 = pi + MAXVAL(tmp, mask=mask_array)
    ELSE ! Else put a mark:
       ra_min2 = -1.0_bp
    END IF

    ! Choose the extremes:
    IF (ra_min1 < 0.0_bp) THEN
       ! There are observations only in 
       ! the interval [pi,2pi[: 
       getMinAndMaxValues_Obss(2,1) = ra_min2
       getMinAndMaxValues_Obss(2,2) = ra_max2
    ELSE IF (ra_min2 < 0.0_bp) THEN
       ! There are observations only in 
       ! the interval [0,pi[: 
       getMinAndMaxValues_Obss(2,1) = ra_min1
       getMinAndMaxValues_Obss(2,2) = ra_max1
    ELSE
       ! There are observations in both intervals:
       IF ((ra_min2 - ra_max1) <= &
            ((two_pi + ra_min1) - ra_max2)) THEN
          ! The distance ra_max1 -> ra_min2 is
          ! shorter than the distance from
          ! ra_max2 -> ra_min1:
          getMinAndMaxValues_Obss(2,1) = ra_min1
          getMinAndMaxValues_Obss(2,2) = ra_max2
       ELSE
          ! The distance ra_max2 -> ra_min1 is
          ! shorter than the distance from
          ! ra_max1 -> ra_min2:
          getMinAndMaxValues_Obss(2,1) = ra_min2
          getMinAndMaxValues_Obss(2,2) = ra_max1
       END IF
    END IF

    ! Declination:
    getMinAndMaxValues_Obss(3,1) = MINVAL(obs_arr(:,3))
    getMinAndMaxValues_Obss(3,2) = MAXVAL(obs_arr(:,3))

    DEALLOCATE(tmp, mask_array, obs_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(obs_arr, stat=err)
       DEALLOCATE(tmp, stat=err)
       DEALLOCATE(mask_array, stat=err)
       CALL errorMessage("Observations / getMinAndMaxValues", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION getMinAndMaxValues_Obss





  !! *Description*:
  !!
  !! Returns the number of different objects in an Observations
  !! object.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object has not been initialized or 
  !!               something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! integer :: different_objects
  !! .
  !! .
  !! .
  !! different_objects = getNrOfObjects(obss)
  !! </pre>
  INTEGER FUNCTION getNrOfObjects(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getNrOfObjects", &
            "Astrometric observations (MPC) have already been initialized.", 1)
       RETURN
    END IF

    getNrOfObjects = this%nobjects

  END FUNCTION getNrOfObjects





  !! *Description*:
  !!
  !! Returns the total number of observations.
  !!
  !! Returns error.
  !!
  INTEGER FUNCTION getNrOfObservations(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getNrOfObservations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getNrOfObservations = this%nobs

  END FUNCTION getNrOfObservations





  !! *Description*:
  !!
  !! Returns the number of this object.
  !!
  !! Returns error if there are more than one numbers.
  !!
  CHARACTER(len=DESIGNATION_LEN) FUNCTION getNumber_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)  :: this
    INTEGER :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getNumber", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getNumber", &
            "Observations missing.", 1)
       RETURN
    END IF

    getNumber_Obss = getNumber(this%obs_arr(1))
    DO i=2,this%nobs
       IF (getNumber_Obss /= getNumber(this%obs_arr(i))) THEN
          error = .TRUE.
          CALL errorMessage("Observations / getNumber", &
               "Different numbers in this set of observations.", 1)
          RETURN          
       END IF
    END DO

  END FUNCTION getNumber_Obss





  !! *Description*:
  !!
  !! Returns the number of observations.
  !!
  !! Returns error.
  !!
  INTEGER FUNCTION getNumberOfObservations(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)  :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getNumberOfObservations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getNumberOfObservations = this%nobs

  END FUNCTION getNumberOfObservations





  !! *Description*:
  !!
  !! Returns the id's of different objects in an Observations
  !! object.
  !!
  !! Returns error.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! character(len=16), dimension(:), pointer :: objects
  !! .
  !! .
  !! .
  !! objects => getObjects(obss)
  !! </pre>
  FUNCTION getObjects(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)                    :: this
    CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: getObjects 

    INTEGER                                               :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObjects", &
            "Astrometric observations (MPC) have already been initialized.", 1)
       RETURN
    END IF

    IF (this%nobjects < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObjects", &
            "Objects missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObjects(this%nobjects), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObjects", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    getObjects = this%objects

  END FUNCTION getObjects





  !! *Description*:
  !!
  !! Returns the observation indicated by index i in this
  !! object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! obs = getObservation(myobservations, 2)
  !!
  FUNCTION getObservation(this, i)

    TYPE (Observations), INTENT(in) :: this
    INTEGER, INTENT(in)             :: i
    TYPE (Observation)              :: getObservation

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < i) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservations", &
            "Number of observation is too large.", 1)
       RETURN
    END IF

    getObservation = copy(this%obs_arr(this%ind(i)))
    IF (error) THEN
       CALL errorMessage("Observations / getObservations", &
            "TRACE BACK", 1)
       CALL NULLIFY(getObservation)
       RETURN
    END IF

  END FUNCTION getObservation





  !! *Description*:
  !!
  !! Returns the the angular arc (radians) between 
  !! the first and the last observation.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getObservationalAngularArc(this)

    TYPE (Observations), INTENT(in) :: this
    REAL(bp)                        :: ra1, dec1, ra2, dec2, cos_alpha

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationalAngularArc", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ra1  = getRA(this%obs_arr(this%ind(1)))
    dec1 = getDec(this%obs_arr(this%ind(1)))
    ra2  = getRA(this%obs_arr(this%ind(this%nobs)))
    dec2 = getDec(this%obs_arr(this%ind(this%nobs)))
    cos_alpha = COS(ra2-ra1)*COS(dec1)*COS(dec2) + SIN(dec1)*SIN(dec2)
    IF (ABS(cos_alpha) > 1.0_bp) THEN
       cos_alpha = SIGN(1.0_bp, cos_alpha)
    END IF
    getObservationalAngularArc = ACOS(cos_alpha)

  END FUNCTION getObservationalAngularArc





  !! *Description*:
  !!
  !! Returns the observational timespan
  !! (&DeltaT) from the first to the last observation
  !! in [Julian] days.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getObservationalTimespan(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationalTimespan", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationalTimespan", &
            "Observations missing.", 1)
       RETURN
    END IF

    getObservationalTimespan = this%criteria(this%ind(this%nobs)) - this%criteria(this%ind(1))
    IF (getObservationalTimespan < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationalTimespan", &
            "Negative time arc: check code.", 1)
       RETURN
    END IF

  END FUNCTION getObservationalTimespan





  !! *Description*:
  !!
  !! Returns an array containing masks for all observations in this
  !! object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! obs_masks => getObservationMasks(myobservations)
  !!
  !FUNCTION getObservationMasks_Obss(this)
  FUNCTION getObservationMasks_Obss(this, use_notes)

    TYPE (Observations), INTENT(in)  :: this
    LOGICAL, INTENT(in), OPTIONAL    :: use_notes
    LOGICAL, DIMENSION(:,:), POINTER :: getObservationMasks_Obss
    LOGICAL                          :: use_notes_

    INTEGER                          :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationMasks", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationMasks", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObservationMasks_Obss(this%nobs,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationMasks", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    IF (PRESENT(use_notes)) THEN
       use_notes_ = use_notes
    ELSE
       use_notes_ = .FALSE.
    END IF

    DO i=1, this%nobs
       !IF (this%obs_note_arr(this%ind(i)) == "#") THEN
       IF (use_notes_ .AND. this%obs_note_arr(this%ind(i)) == "#") THEN
          getObservationMasks_Obss(i,1:6) = .FALSE.
       ELSE
          getObservationMasks_Obss(i,1:6) = getObservationMask(this%obs_arr(this%ind(i)))
          IF (error) THEN
             CALL errorMessage("Observations / getObservationMasks", &
                  "TRACE BACK", 1)
             DEALLOCATE(getObservationMasks_Obss, stat=err)
             RETURN
          END IF
       END IF
    END DO

  END FUNCTION getObservationMasks_Obss






  !! *Description*:
  !!
  !! Returns an array containing notes for all observations in this
  !! object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! obs_notes => getObservationNotes(myobservations)
  !!
  FUNCTION getObservationNotes(this)

    TYPE (Observations), INTENT(in)  :: this
    CHARACTER(4), DIMENSION(:), POINTER :: getObservationNotes

    INTEGER                          :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationNotes", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationNotes", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObservationNotes(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationNotes", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getObservationNotes(i) = this%obs_note_arr(this%ind(i))
    END DO

  END FUNCTION getObservationNotes





  !! *Description*:
  !!
  !! Returns an array containing all observations in this
  !! object.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! obs => getObservations(myobservations)
  !!
  FUNCTION getObservations(this)

    TYPE (Observations), INTENT(in)           :: this
    TYPE (Observation), DIMENSION(:), POINTER :: getObservations

    INTEGER                                   :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / get", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObservations(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservations", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getObservations(i) = copy(this%obs_arr(this%ind(i)))
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / getObservations", &
               "TRACE BACK", 1)
          DEALLOCATE(getObservations, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getObservations





  !! *Description*:
  !!

  !! Returns observation weights for all 6 observables (distance, ra,
  !! dec, and their time derivatives) based on the relative standard
  !! deviations (stdev): the observation with the smallest stdev gets
  !! weight=1, while the largest stdevs corresponds to smaller
  !! weights. A negative value indicates that a measurement does not
  !! exist (cfg. obs_mask = .false.).
  !!
  !! Example:
  !!
  !! obs1_stdev_ra = 0.10''  obs1_stdev_dec = 0.05''
  !! obs2_stdev_ra = 0.50''  obs2_stdev_dec = 0.10''
  !!
  !! ==> 
  !!
  !! obs1_weight_ra = 1.0   obs1_weight_dec = 1.0
  !! obs2_weight_ra = 0.2   obs2_weight_dec = 0.5
  !!
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! weights => getObservationWeights(myobservations)
  !!
  FUNCTION getObservationWeights(this)

    TYPE (Observations), INTENT(in)   :: this
    REAL(bp), DIMENSION(:,:), POINTER :: getObservationWeights

    REAL(bp)                          :: stdev_min
    INTEGER                           :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationWeights", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationWeights", &
            "Observations missing.", 1)
       RETURN
    END IF

    getObservationWeights => getStandardDeviations(this)
    IF (error) THEN
       CALL errorMessage("Observations / getObservationWeights", &
            "TRACE BACK", 1)
       DEALLOCATE(getObservationWeights, stat=err)
       RETURN
    END IF

    ! For each observable..
    DO i=1,6
       ! Find the smallest stdev: 
       stdev_min = MINVAL(ABS(getObservationWeights(:,i)))
       ! The weight is equal to the inverse of the ratio of the stdev
       ! and the smallest stdev:
       getObservationWeights(:,i) = stdev_min / getObservationWeights(:,i)
    END DO

  END FUNCTION getObservationWeights





  !! *Description*:
  !!
  !! Returns an array of spherical coordinates corresponding to the
  !! observations.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! scoords => getObservationSCoords(myobservations)
  !!
  FUNCTION getObservationSCoords(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)                    :: this
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: getObservationSCoords

    INTEGER                                            :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationSCoords", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationSCoords", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObservationSCoords(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservationSCoords", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getObservationSCoords(i) = getObservationSCoord(this%obs_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getObservationSCoords", &
               "TRACE BACK", 1)
          DEALLOCATE(getObservationSCoords, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getObservationSCoords





  !! *Description*:
  !!
  !! Returns an array of Cartesian coordinates corresponding to the
  !! heliocentric locations of the observers at times corresponding to
  !! the observations.
  !!
  !! Returns error.
  !!
  !!
  !! *Usage*:
  !!
  !! ccoords => getObservatoryCCoords(myobservations)
  !!
  FUNCTION getObservatoryCCoords(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)                    :: this
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: getObservatoryCCoords

    INTEGER                                            :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservatoryCCoords", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservatoryCCoords", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObservatoryCCoords(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservatoryCCoords", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getObservatoryCCoords(i) = getObservatoryCCoord(this%obs_arr(this%ind(i))) 
       IF (error) THEN
          CALL errorMessage("Observations / getObservatoryCCoords", &
               "TRACE BACK", 1)
          DEALLOCATE(getObservatoryCCoords, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getObservatoryCCoords




  !! *Description*:
  !!
  !! Returns the observatory codes for the observatories where these observations were made.
  !!
  !! Returns error.
  !!
  !!
  !! *Usage*:
  !!
  !! codes => getObservatoryCodes(myobservations)
  !!
  FUNCTION getObservatoryCodes_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)                    :: this
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: getObservatoryCodes_Obss

    INTEGER                                            :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservatoryCodes", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservatoryCodes", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getObservatoryCodes_Obss(this%nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getObservatoryCodes", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, this%nobs
       getObservatoryCodes_Obss(i) = getObservatoryCode(this%obs_arr(this%ind(i))) 
       IF (error) THEN
          CALL errorMessage("Observations / getObservatoryCodes", &
               "TRACE BACK", 1)
          DEALLOCATE(getObservatoryCodes_Obss, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getObservatoryCodes_Obss




  !! *Description*:
  !!
  !! Returns an Observations array containing the
  !! input Observations object separated into Observations objects
  !! with different (1) numbers or (2) designations.
  !!
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if this object has not been initialized or something goes wrong.
  !! 
  !! *Usage*:
  !!
  !! <pre>
  !! type (Observations) :: obss
  !! type (Observations), dimension(:), pointer :: separated_sets
  !! .
  !! .
  !! .
  ! separated_sets => getSeparatedSets(obss)
  !! </pre>
  !!
  FUNCTION getSeparatedSets(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)         :: this
    TYPE (Observations), DIMENSION(:), POINTER :: getSeparatedSets

    CHARACTER(len=DESIGNATION_LEN)             :: number, id
    INTEGER                                    :: i, j, err

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getSeparatedSets", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getSeparatedSets", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getSeparatedSets(this%nobjects), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getSeparatedSets", &
            "Could not allocate memory for separated objects.", 1)
       RETURN
    END IF

    DO i=1, this%nobjects
       CALL NULLIFY(getSeparatedSets(i))
       CALL NEW(getSeparatedSets(i))
       IF (error) THEN
          CALL errorMessage("Observations / getSeparatedSets", &
               "TRACE BACK 2", 1)
          DEALLOCATE(getSeparatedSets, stat=err)
          RETURN
       END IF
    END DO

    DO i=1, this%nobs

       ! (1) number or (2) designation: 
       number = getNumber(this%obs_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getSeparatedSets", &
               "TRACE BACK 3", 1)
          DEALLOCATE(getSeparatedSets, stat=err)
          RETURN
       END IF
       IF (LEN_TRIM(number) /= 0) THEN
          id = number
       ELSE
          id = getDesignation(this%obs_arr(this%ind(i)))
          IF (error) THEN
             CALL errorMessage("Observations / getSeparatedSets", &
                  "TRACE BACK 5", 1)
             DEALLOCATE(getSeparatedSets, stat=err)
             RETURN
          END IF
       END IF
       DO j=1, this%nobjects
          IF (id == this%objects(j)) EXIT
       END DO
       CALL addObservation(getSeparatedSets(j), &
            this%obs_arr(this%ind(i)), &
            note=this%obs_note_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getSeparatedSets", &
               "TRACE BACK 6", 1)
          DEALLOCATE(getSeparatedSets, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getSeparatedSets





  !! *Description*:
  !!
  !! 
  !!
  !! Returns error.
  !!
  FUNCTION getStandardDeviations_Obss(this)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)   :: this
    REAL(bp), DIMENSION(:,:), POINTER :: getStandardDeviations_Obss

    INTEGER                           :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getStandardDeviations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getStandardDeviations", &
            "Observations missing.", 1)
       RETURN
    END IF

    ALLOCATE(getStandardDeviations_Obss(this%nobs,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / getStandardDeviations", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1,this%nobs
       getStandardDeviations_Obss(i,1:6) = getStandardDeviations(this%obs_arr(this%ind(i)))
       IF (error) THEN
          CALL errorMessage("Observations / getStandardDeviations", &
               "TRACE BACK", 1)
          DEALLOCATE(getStandardDeviations_Obss, stat=err)
          RETURN
       END IF
    END DO

  END FUNCTION getStandardDeviations_Obss





  !! *Description*:
  !!
  !! Groups a set of observations into subgroups based on observation
  !! date and number of observations.
  !!
  !! Returns error.
  !!
  SUBROUTINE groupObservations(this, time_interval, nobs_max, obss_arr, tdt_arr)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in)   :: this
    REAL(bp), INTENT(in) :: time_interval
    INTEGER, INTENT(in) :: nobs_max
    TYPE (Observations), DIMENSION(:), POINTER :: obss_arr
    REAL(bp), DIMENSION(:,:), POINTER :: tdt_arr

    REAL(bp) :: tdt_, tdt_sum 
    INTEGER  :: i, j, nobs, nobs_max_, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / groupObservations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs < 1) THEN
       error = .TRUE.
       CALL errorMessage("Observations / groupObservations", &
            "Observations missing.", 1)
       RETURN
    END IF

    IF (nobs_max < 0) THEN
       nobs_max_ = HUGE(nobs_max_)
    ELSE
       nobs_max_ = nobs_max
    END IF

    tdt_ = this%criteria(this%ind(1))
    i = 0
    j = 0
    DO
       IF (j+1 > this%nobs) THEN
          EXIT
       END IF
       i = i + 1
       tdt_arr => reallocate(tdt_arr,i,2)
       tdt_arr(i,2) = tdt_
       obss_arr => reallocate(obss_arr,i)
       CALL NULLIFY(obss_arr(i))
       CALL NEW(obss_arr(i))
       IF (error) THEN
          CALL errorMessage("Observations / groupObservations", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(tdt_arr, stat=err)
          DEALLOCATE(obss_arr, stat=err)
          RETURN
       END IF
       nobs = 0
       tdt_sum = 0.0_bp
       DO
          IF (j+1 > this%nobs) THEN
             EXIT
          ELSE IF (this%criteria(this%ind(j+1)) - tdt_ > time_interval &
               .OR. nobs+1 > nobs_max_) THEN
             tdt_ = this%criteria(this%ind(j+1))
             EXIT
          ELSE
             nobs = nobs + 1
             j = j + 1
             tdt_sum = tdt_sum + this%criteria(this%ind(j))
             CALL addObservation(obss_arr(i), this%obs_arr(this%ind(j)))
             IF (error) THEN
                CALL errorMessage("Observations / groupObservations", &
                     "TRACE BACK (5)", 1)
                DEALLOCATE(tdt_arr, stat=err)
                DEALLOCATE(obss_arr, stat=err)
                RETURN
             END IF
          END IF
       END DO
       tdt_arr(i,1) = tdt_sum/REAL(nobs,bp)
       tdt_arr(i,2) = this%criteria(this%ind(j)) - tdt_arr(i,2)
    END DO
    tdt_arr(:,1) = tdt_arr(:,1) - tdt_arr(1,1)

  END SUBROUTINE groupObservations





  !! *Description*:
  !!
  !! Intrinsic routine that reads a file of the given format.
  !!
  !! Returns error.
  !!
  SUBROUTINE readObservationFile(this, obsf, stdev, orb_sim)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)           :: this
    TYPE (File), INTENT(in)                      :: obsf
    REAL(bp), DIMENSION(6), INTENT(in), OPTIONAL :: stdev
    TYPE (Orbit), INTENT(inout), OPTIONAL        :: orb_sim

    TYPE (File) :: orbfile
    TYPE (Orbit) :: orb
    TYPE (Observatories) :: obsies
    TYPE (Observatory) :: obsy
    TYPE (Time) :: t, t0
    TYPE (SphericalCoordinates) :: obs_scoord, ephemeris
    TYPE (CartesianCoordinates) :: obsy_ccoord, geocenter_ccoord, &
         satellite_ccoord
    CHARACTER(len=512) :: line
    CHARACTER(len=FNAME_LEN) :: fname, suffix
    CHARACTER(len=DESIGNATION_LEN) :: number, designation
    CHARACTER(len=132), DIMENSION(5) :: records
    CHARACTER(len=132) :: record
    CHARACTER(len=124) :: line1, line2, line_, str, filter, obstype
    CHARACTER(len=7) :: nr_str
    CHARACTER(len=4) :: timescale, obsy_code, strk_class
    REAL(bp), DIMENSION(:,:), POINTER :: planeph => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: element_arr
    REAL(bp), DIMENSION(6,6) :: covariance
    REAL(bp), DIMENSION(6) :: coordinates, stdev_, mean
    REAL(bp), DIMENSION(3) :: position, velocity, pos1, pos2
    REAL(bp) :: day, sec, arcsec, mag, ra, dec, jd, mjd_utc, dt, &
         ecl_lon, ecl_lat, angscan, pos_unc_along, pos_unc_across, &
         vel_unc_along, vel_unc_across, rot_angle, correlation, &
         rmin, rarcmin, mag_unc, s2n, mjd_tt, mjd_tcb, strk_len, &
         strk_len_unc, strk_direction, strk_direction_unc
    INTEGER :: i, j, err, year, month, hour, min, deg, arcmin, &
         nlines, coord_unit, indx, iobs, irecord, norb, ccd, &
         border1, border2
    LOGICAL, DIMENSION(6) :: obs_mask
    LOGICAL :: discovery, converttonewformat

    converttonewformat = .FALSE.

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readObservationFile", &
            "Astrometric observations (MPC) have already been initialized.", 1)
       RETURN
    END IF

    nlines = getNrOfLines(obsf)
    IF (error) THEN
       CALL errorMessage("Observations / readObservationFile", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    ALLOCATE(this%obs_arr(nlines), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readObservationFile", &
            "Could not allocate memory for the observations.", 1)
       RETURN
    END IF

    CALL NEW(obsies)
    IF (error) THEN
       CALL errorMessage("Observations / readObservationFile", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

    ALLOCATE(this%obs_note_arr(nlines), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readObservationFile", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    this%obs_note_arr = " "

    ! Decide which routine should read the input 
    ! file by checking the suffix:
    fname = getFileName(obsf)
    indx = INDEX(fname, ".", back=.TRUE.)
    IF (indx == 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / new", &
            "Suffix required in file name.",1)
    END IF
    suffix = fname(indx+1:LEN_TRIM(fname))
    CALL locase(suffix, error)
    IF (error) THEN
       CALL errorMessage("Observations / new", &
            "The suffix string contains forbidden characters.", 1)
       RETURN
    END IF

    IF (info_verb >= 1 .AND. PRESENT(stdev) .AND. &
         suffix == "mpc2") THEN
       WRITE(stdout,"(1X,A)") "WARNING! Overriding " // &
            "observation uncertainties given in the " // &
            "observation file with the ones given in " // &
            "the option file."
    END IF

    SELECT CASE (TRIM(suffix))

    CASE ("des")

       ! This is the Data Exchange Standard adopted by MOPS.

       ! Read first line to check if the header exists:
       READ(getUnit(obsf), "(A)", iostat=err) line1
       IF (.NOT.(line1(1:2) == "!!" .OR. line1(1:1) == "#")) THEN
          REWIND(getUnit(obsf))
       END IF
       i = 0
       DO
          designation = " "
          covariance = 0.0
          READ(getUnit(obsf), *, iostat=err) designation, mjd_utc, &
               obstype, ra, dec, mag, filter, obsy_code, &
               covariance(2,2), covariance(3,3), mag_unc, s2n, str
          IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading DES observations from file.", 1)
             RETURN
          ELSE IF (err < 0) THEN ! end-of-file
             EXIT
          ELSE IF (designation(1:1) == "!" .OR. designation(1:1) == "#") THEN
             CYCLE
          END IF
          CALL NEW(t, mjd_utc, "UTC")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (15)", 1)
             RETURN
          END IF
          CALL NEW(obs_scoord, ra*rad_deg, dec*rad_deg, t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
          IF (PRESENT(stdev)) THEN
             DO j=1,6
                covariance(j,j) = stdev(j)**2.0_bp
             END DO
          ELSE
             covariance(2,2) = (covariance(2,2)*rad_asec)**2.0_bp
             covariance(3,3) = (covariance(3,3)*rad_asec)**2.0_bp
          END IF
          IF (obstype == "O" .OR. obstype == "o") THEN
             obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)
          END IF
          obsy = getObservatory(obsies, obsy_code)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (25)", 1)
             RETURN
          END IF
          obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (30)", 1)
             RETURN
          END IF
          IF (filter == "X") THEN
             filter = " "
          END IF
          i = i + 1
          CALL NEW(this%obs_arr(i), number="", designation=designation, &
               discovery=.FALSE., note1="", note2="", &
               obs_scoord=obs_scoord, covariance=covariance, &
               obs_mask=obs_mask, mag=mag, mag_unc=mag_unc, &
               filter=filter, s2n=s2n, obsy=obsy, &
               obsy_ccoord=obsy_ccoord, secret_name=TRIM(str))
          CALL NULLIFY(t)
          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)
       END DO

    CASE ("mpc") ! Minor Planet Center old format

       ! This is the file format (http://cfa-www.harvard.edu/iau/info/OpticalObs.html):
       ! to be read from the observations file:
       !
       ! Nr.,Des.,New       Time                R.A.                 
       !   ---15---- --------17--------- --------12---------      
       ! "(I5,A7,3A1,I4,1X,I2,1X,F8.5,1X,I2,1X,I2,1X,F5.2,1X," // &
       !         Dec.           Mag.,Filter   Obs. 
       !  --------12---------    ---6---    -3-
       ! "I3,1X,I2,1X,F4.1,1X,9X,F5.3,A1,6X,A3)"

       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)
       i = 0
       DO

          READ(getUnit(obsf), "(A124)", iostat=err) line_
          IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file.", 1)
             RETURN
          ELSE IF (err < 0) THEN ! end-of-file
             EXIT
          END IF

          IF (info_verb >= 4) THEN
             WRITE(stdout,"(A80)") line_
          END IF

          ! If the line is empty or marked with "#", read the next line:
          IF (LEN_TRIM(line_) == 0 .OR. line_(1:1) == "#") THEN
             CYCLE
          END IF
          ! Skip observation if line contains radar observations
          IF (line_(15:15) == "r" .OR. line_(15:15) == "R") THEN
             CYCLE
          END IF
          ! Skip observation if line contains "Roving Observer" observations
          IF (line_(15:15) == "v" .OR. line_(15:15) == "V") THEN
             CYCLE
          END IF

          i = i + 1

          ! Read additional note from file
          IF (line_(1:1) == "!") THEN
             str = line_(2:2)
             CALL removeLeadingBlanks(str)
             this%obs_note_arr(i) = TRIM(str)
             line1 = line_(3:)
          ELSE
             this%obs_note_arr(i) = " "
             line1 = line_
          END IF

          number = line1(1:5)
          designation = line1(6:12)
          IF (line1(13:13) == "*") THEN
             discovery = .TRUE.
          ELSE
             discovery = .FALSE.
          END IF
          CALL toInt(line1(16:19), year, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (5).", 1)
             RETURN
          END IF
          CALL toInt(line1(21:22), month, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (10).", 1)
             RETURN
          END IF
          CALL toReal(line1(24:32), day, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (15).", 1)
             RETURN
          END IF
          CALL toInt(line1(33:34), hour, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (20).", 1)
             RETURN
          END IF
          CALL toInt(line1(36:37), min, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (25).", 1)
             RETURN
          END IF
          CALL toReal(line1(39:44), sec, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (30).", 1)
             RETURN
          END IF
          CALL toInt(line1(46:47), deg, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (35).", 1)
             RETURN
          END IF
          CALL toInt(line1(49:50), arcmin, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (40).", 1)
             RETURN
          END IF
          CALL toReal(line1(52:56), arcsec, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (45).", 1)
             RETURN
          END IF
          IF (LEN_TRIM(line1(66:70)) /= 0) THEN
             CALL toReal(line1(66:70), mag, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (50).", 1)
                RETURN
             END IF
          ELSE
             mag = 99.9_bp
          END IF

          ! R.A. and Dec. + epoch:
          !
          ! The timescale used by MPC prior to year 1972 is UT1 and since
          ! then UTC has been used:  
          IF (year < 1972) THEN
             timescale = "UT1"
          ELSE
             timescale = "UTC"
          END IF
          ! Create epoch object:
          CALL NEW(t, year, month, day, timescale)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (35)", 1)
             RETURN
          END IF
          IF (line1(45:45) == "-" .AND. deg /= 0) THEN
             deg = -1 * deg
          ELSE IF (line1(45:45) == "-" .AND. deg == 0 .AND. arcmin /= 0) THEN
             arcmin = -1 * arcmin
          ELSE IF (line1(45:45) == "-" .AND. deg == 0 .AND. arcmin == 0) THEN
             arcsec = -1.0_bp * arcsec
          END IF
          ! Create observation object containing R.A. and Dec. + epoch:
          CALL NEW(obs_scoord, hour, min, sec, deg, arcmin, arcsec, t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (40)", 1)
             RETURN
          END IF
          IF (PRESENT(stdev)) THEN
             stdev_ = stdev
          ELSE
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Standard deviation for MPC format must be explicitly given.", 1)             
             RETURN
          END IF
          covariance = 0.0_bp
          covariance(2,2) = stdev_(2)**2.0_bp
          covariance(3,3) = stdev_(3)**2.0_bp

          ! Compute the heliocentric position of the observer at epoch t:
          obsy_code = " "
          obsy_code(1:3) = line1(78:80)
          SELECT CASE (obsy_code(1:3))
          CASE (" -3")
             obsy_code = "500"
          END SELECT
          obsy = getObservatory(obsies, TRIM(obsy_code))
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (45)", 1)
             RETURN
          END IF

          ! Change old type descriptions to new ones
          SELECT CASE (line1(15:15))
          CASE ("c", " ")
             line1(15:15) = "C"
          END SELECT

          ! Create observation object:
          CALL NULLIFY(this%obs_arr(i))
          SELECT CASE (line1(15:15))
          CASE ("S")
             READ(getUnit(obsf), "(A80)", iostat=err) line2
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Error while reading satellite coordinates from file.", 1)
                RETURN
             END IF
             CALL toReal(line2(36:46), position(1), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (55).", 1)
                RETURN
             END IF
             IF (line2(35:35) == "-") THEN
                position(1) = -1.0_bp*position(1)
             END IF
             CALL toReal(line2(48:58), position(2), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (60).", 1)
                RETURN
             END IF
             IF (line2(47:47) == "-") THEN
                position(2) = -1.0_bp*position(2)
             END IF
             CALL toReal(line2(60:70), position(3), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (50).", 1)
                RETURN
             END IF
             IF (line2(59:59) == "-") THEN
                position(3) = -1.0_bp*position(3)
             END IF
             SELECT CASE (line2(33:33))
             CASE ("1", "2")
                ! 
                IF (line2(33:33) == "1") THEN
                   coord_unit = 1
                   position(1:3) = position(1:3)/km_au
                ELSE
                   coord_unit = 2
                END IF
                geocenter_ccoord = getObservatoryCCoord(obsies, "500", t)
                CALL rotateToEquatorial(geocenter_ccoord)
                coordinates = getCoordinates(geocenter_ccoord)
                coordinates(1:3) = coordinates(1:3) + position(1:3)
                CALL NEW(obsy_ccoord, coordinates, "equatorial", t)
                coordinates = 0.0_bp
                coordinates(1:3) = position(1:3)
                CALL NEW(satellite_ccoord, coordinates, "equatorial", t)
             CASE ("3")
                coordinates = 0.0_bp
                coordinates(1:3) = position(1:3)
                CALL NEW(obsy_ccoord, coordinates, "ecliptic", t)
                CALL rotateToEquatorial(obsy_ccoord)
                geocenter_ccoord = getObservatoryCCoord(obsies, "500", t)
                CALL rotateToEcliptic(geocenter_ccoord)
                coordinates = getCoordinates(geocenter_ccoord)
                coordinates(1:3) = position(1:3) - coordinates(1:3)
                CALL NEW(satellite_ccoord, coordinates, "ecliptic", t)
                CALL rotateToEquatorial(satellite_ccoord)
                coord_unit = 2
             END SELECT
             CALL NEW(this%obs_arr(i), number, designation, discovery, &
                  line1(14:14), line1(15:15), obs_scoord, covariance, &
                  obs_mask, mag, -1.0_bp, line1(71:71), -1.0_bp, &
                  obsy, obsy_ccoord, &
                  satellite_ccoord=satellite_ccoord, coord_unit=coord_unit)
          CASE default
             obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)
             IF (PRESENT(orb_sim)) THEN
                CALL NULLIFY(obs_scoord)
                CALL getEphemeris(orb_sim, obsy_ccoord, obs_scoord)
                mean = 0.0_bp
                CALL addMultinormalDeviate(obs_scoord, mean, covariance)
                this%obs_note_arr(i) = "S"
             END IF
             CALL NEW(this%obs_arr(i), number=number, &
                  designation=designation, discovery=discovery, &
                  note1=line1(14:14), note2=line1(15:15), &
                  obs_scoord=obs_scoord, covariance=covariance, &
                  obs_mask=obs_mask, mag=mag, filter=line1(71:71), &
                  obsy=obsy, obsy_ccoord=obsy_ccoord)
          END SELECT
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (55)", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(geocenter_ccoord)
          CALL NULLIFY(obsy_ccoord)
          CALL NULLIFY(satellite_ccoord)

       END DO

    CASE ("mpc2") ! Minor Planet Center new format (20060601->)

       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)
       iobs = 0
       obsloop: DO

          records = " "
          irecord = 0
          recordloop: DO
             record = " "
             READ(getUnit(obsf), "(A132)", iostat=err) record
             IF (err > 0) THEN ! error
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Error while reading observations from file: " // TRIM(fname), 1)
                RETURN
             ELSE IF (err < 0 .AND. irecord > 0) THEN ! problematic end-of-file
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Observation file " // TRIM(fname) // " ended abruptly.", 1)
                RETURN
             ELSE IF (err < 0 .AND. irecord == 0) THEN ! end-of-file
                EXIT obsloop
             END IF
             ! If the line is empty or marked with "#", read the next line:
             IF (LEN_TRIM(record) == 0 .OR. record(1:1) == "#") THEN
                CYCLE recordloop
             END IF
             irecord = irecord + 1
             IF (irecord /= 1 .AND. (record(17:17) == "1" .OR. record(17:17) == "A")) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Records missing from observation file.", 1)
                RETURN
             END IF
             IF (irecord > 5) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "End-of-record sign missing from observation record.", 1)
                RETURN
             END IF
             IF (record(17:17) == "1") iobs = iobs + 1
             records(irecord) = record 
             IF (records(irecord)(18:18) == "+") EXIT ! all available records for this observation read  
          END DO recordloop

          IF (info_verb >= 4) THEN
             DO irecord=1,5
                IF (LEN_TRIM(records(irecord)) /= 0) WRITE(stdout,"(A132)") records(irecord)
             END DO
          END IF

          number = records(1)(1:7)
          designation = records(1)(8:16)
          IF (records(1)(18:18) == "*") THEN
             discovery = .TRUE.
          ELSE
             discovery = .FALSE.
          END IF
          CALL toInt(records(1)(20:23), year, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (5).", 1)
             RETURN
          END IF
          CALL toInt(records(1)(25:26), month, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (10).", 1)
             RETURN
          END IF
          CALL toReal(records(1)(28:40), day, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (15).", 1)
             RETURN
          END IF
          ! Create epoch object:
          ! The timescale used by MPC prior to year 1972 is UT1 and since
          ! then UTC has been used:  
          IF (year < 1972) THEN
             timescale = "UT1"
          ELSE
             timescale = "UTC"
          END IF
          CALL NEW(t, year, month, day, timescale)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (60)", 1)
             RETURN
          END IF
          IF (records(1)(43:43) == " ") THEN
             CALL toInt(records(1)(41:42), hour, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (20).", 1)
                RETURN
             END IF
             CALL toInt(records(1)(44:45), min, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (25).", 1)
                RETURN
             END IF
             CALL toReal(records(1)(47:52), sec, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (30).", 1)
                RETURN
             END IF
             CALL toInt(records(1)(57:58), deg, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (35).", 1)
                RETURN
             END IF
             CALL toInt(records(1)(60:61), arcmin, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (40).", 1)
                RETURN
             END IF
             CALL toReal(records(1)(63:67), arcsec, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (45).", 1)
                RETURN
             END IF
             IF (records(1)(56:56) == "-" .AND. deg /= 0) THEN
                deg = -1 * deg
             ELSE IF (records(1)(56:56) == "-" .AND. deg == 0 .AND. arcmin /= 0) THEN
                arcmin = -1 * arcmin
             ELSE IF (records(1)(56:56) == "-" .AND. deg == 0 .AND. arcmin == 0) THEN
                arcsec = -1.0_bp * arcsec
             END IF
             ! Create observation object containing R.A. and Dec. + epoch:
             CALL NEW(obs_scoord, hour, min, sec, deg, arcmin, arcsec, t)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "TRACE BACK (65)", 1)
                RETURN
             END IF
          ELSE
             CALL toReal(records(1)(41:55), ra, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (50).", 1)
                RETURN
             END IF
             CALL toReal(records(1)(56:71), dec, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (55).", 1)
                RETURN
             END IF
             ! Create observation object containing R.A. and Dec. + epoch:
             CALL NEW(obs_scoord, longitude=ra, latitude=dec, t=t)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "TRACE BACK (70)", 1)
                RETURN
             END IF
          END IF
          IF (LEN_TRIM(records(1)(72:77)) /= 0) THEN
             CALL toReal(records(1)(72:77), mag, error)
          ELSE
             mag = 99.9_bp
          END IF

          ! Uncertainties:
          stdev_ = 0.0_bp
          IF (PRESENT(stdev)) THEN
             stdev_ = stdev
             correlation = 0.0_bp
          ELSE
             IF (records(2)(43:43) == ".") THEN
                CALL toReal(records(2)(41:46), stdev_(2), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (60).", 1)
                   RETURN
                END IF
                CALL toReal(records(2)(48:53), stdev_(3), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (65).", 1)
                   RETURN
                END IF
                stdev_ = stdev_*rad_asec
             ELSE
                CALL toReal(records(2)(41:41) // "." // &
                     records(2)(42:48), stdev_(2), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (70).", 1)
                   RETURN
                END IF
                CALL toReal(records(2)(49:49) // "." // &
                     records(2)(50:56), stdev_(3), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (75).", 1)
                   RETURN
                END IF
                ! Transform units to radians from 10**(-7) radians:
                stdev_ = stdev_*0.0000001_bp
             END IF
             CALL toReal(records(2)(57:64), correlation, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (77).", 1)
                RETURN
             END IF
          END IF
          covariance = 0.0_bp
          covariance(2,2) = stdev_(2)**2.0_bp
          covariance(3,3) = stdev_(3)**2.0_bp
          covariance(2,3) = correlation*stdev_(2)*stdev_(3)
          covariance(3,2) = covariance(2,3)

          ! Compute the heliocentric position of the observer at epoch t:
          obsy_code = records(1)(129:132)
          IF (obsy_code(1:1) == "0") THEN
             obsy_code = obsy_code(2:4)
             obsy_code(4:4) = " "
          END IF
          obsy = getObservatory(obsies, obsy_code)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (80)", 1)
             RETURN
          END IF

          ! Create observation object:
          CALL NULLIFY(this%obs_arr(iobs))
          SELECT CASE (records(1)(19:19))
          CASE (" ", "C", "A", "T", "M", "P")
             obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)
             CALL NEW(this%obs_arr(iobs), number=number, designation=designation, &
                  discovery=discovery, note1=records(1)(124:124), &
                  note2=records(1)(19:19), obs_scoord=obs_scoord, &
                  covariance=covariance, obs_mask=obs_mask, &
                  mag=mag, filter=records(1)(78:79), obsy=obsy, &
                  obsy_ccoord=obsy_ccoord)
          CASE ("S")
             IF (LEN_TRIM(records(3)) == 0) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "No satellite coordinates available in observation file.", 1)
                RETURN
             END IF
             CALL toReal(records(3)(42:59), position(1), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (80).", 1)
                RETURN
             END IF
             CALL toReal(records(3)(60:77), position(2), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (85).", 1)
                RETURN
             END IF
             CALL toReal(records(3)(78:95), position(3), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (90).", 1)
                RETURN
             END IF
             SELECT CASE (records(3)(41:41))
             CASE ("0")
                coord_unit = 1
                position(1:3) = position(1:3)/km_au
             CASE ("1")
                coord_unit = 2
             CASE default
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Unit type "" // records(3)(41:41) // "" not available.", 1)
                RETURN
             END SELECT
             obsy_ccoord = getObservatoryCCoord(obsies, "500", t)
             CALL rotateToEquatorial(obsy_ccoord)
             coordinates = getCoordinates(obsy_ccoord)
             CALL NULLIFY(obsy_ccoord)
             coordinates(1:3) = coordinates(1:3) + position(1:3)
             CALL NEW(obsy_ccoord, coordinates, "equatorial", t)
             coordinates = 0.0_bp
             coordinates(1:3) = position(1:3)
             CALL NEW(satellite_ccoord, coordinates, "equatorial", t)
             CALL NEW(this%obs_arr(iobs), number=number, designation=designation, &
                  discovery=discovery, note1=records(1)(124:124), &
                  note2=records(1)(19:19), obs_scoord=obs_scoord, &
                  covariance=covariance, obs_mask=obs_mask, &
                  mag=mag, filter=records(1)(78:79), obsy=obsy, &
                  obsy_ccoord=obsy_ccoord, satellite_ccoord=satellite_ccoord, &
                  coord_unit=coord_unit)
          CASE default
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "No such (" // records(1)(19:19) // ") option available.", 1)
             IF (err_verb >= 1) THEN
                WRITE(stderr,"(A)") records(1)
             END IF
             RETURN
          END SELECT
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (85)", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(geocenter_ccoord)
          CALL NULLIFY(obsy_ccoord)
          CALL NULLIFY(satellite_ccoord)

       END DO obsloop

       i = iobs

    CASE ("mpc3","sor","vov","kep") ! Minor Planet Center new format (20060601->)

       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)
       iobs = 0
       obsloop_mpc3: DO

          records = " "
          irecord = 0
          recordloop_mpc3: DO
             READ(getUnit(obsf), "(A132)", iostat=err) record
             IF (err > 0) THEN ! error
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Error while reading observations from file: " // TRIM(fname), 1)
                RETURN
             ELSE IF (err < 0 .AND. irecord > 0) THEN ! problematic end-of-file
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Observation file " // TRIM(fname) // " ended abruptly.", 1)
                RETURN
             ELSE IF (err < 0 .AND. irecord == 0) THEN ! end-of-file
                EXIT obsloop_mpc3
             END IF
             ! If the line is empty or marked with "#", read the next line:
             IF (LEN_TRIM(record) == 0 .OR. record(1:1) == "#" .OR. &
                  record(7:7) == "X") THEN
                records = " "
                irecord = 0
                CYCLE recordloop_mpc3
             END IF
             irecord = irecord + 1
             IF (irecord /= 1 .AND. (record(19:19) == "1" .OR. record(19:19) == "A")) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Records missing from observation file.", 1)
                RETURN
             END IF
             IF (irecord > 5) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "End-of-record sign missing from observation record.", 1)
                RETURN
             END IF
             IF (record(19:19) == "1") THEN
                iobs = iobs + 1
             END IF
             records(irecord) = record 
             IF (records(irecord)(20:20) == "+") THEN
                ! All available records for this observation read  
                EXIT
             END IF
          END DO recordloop_mpc3

          IF (info_verb >= 4) THEN
             DO irecord=1,5
                IF (LEN_TRIM(records(irecord)) /= 0) THEN
                   WRITE(stdout,"(A132)") records(irecord)
                END IF
             END DO
          END IF

          SELECT CASE (records(1)(7:7))
          CASE ("C", "P", "D", "X", "A")
             number = records(1)(1:6)
             designation = records(1)(7:16)
          CASE default
             number = records(1)(1:7)
             designation = records(1)(8:16)
          END SELECT
          IF (records(1)(20:20) == "*") THEN
             discovery = .TRUE.
          ELSE
             discovery = .FALSE.
          END IF
          CALL toInt(records(1)(22:25), year, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string (" // records(1)(22:25) &
                  // ") to number.", 1)
             RETURN
          END IF
          CALL toInt(records(1)(26:27), month, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (10).", 1)
             RETURN
          END IF
          CALL toReal(records(1)(28:39), day, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (15).", 1)
             RETURN
          END IF
          ! Create epoch object:
          ! The timescale used by MPC prior to year 1972 is UT1 and since
          ! then UTC has been used:  
          IF (year < 1972) THEN
             timescale = "UT1"
          ELSE
             timescale = "UTC"
          END IF
          CALL NEW(t, year, month, day, timescale)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (90)", 1)
             RETURN
          END IF
          IF (records(1)(43:43) == " ") THEN
             CALL toInt(records(1)(41:42), hour, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (20).", 1)
                RETURN
             END IF
             IF (LEN_TRIM(records(1)(47:52)) == 0 .OR. records(1)(46:46) == ".") THEN
                CALL toReal(records(1)(44:52), rmin, error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (22).", 1)
                   RETURN
                END IF
                min = FLOOR(rmin)
                sec = 60.0_bp*(rmin-min)
             ELSE
                CALL toInt(records(1)(44:45), min, error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (25).", 1)
                   RETURN
                END IF
                CALL toReal(records(1)(47:52), sec, error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (30).", 1)
                   RETURN
                END IF
             END IF
             CALL toInt(records(1)(57:58), deg, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (35).", 1)
                RETURN
             END IF
             IF (LEN_TRIM(records(1)(63:67)) == 0 .OR. records(1)(62:62) == ".") THEN
                CALL toReal(records(1)(60:67), rarcmin, error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (32).", 1)
                   RETURN
                END IF
                arcmin = FLOOR(rarcmin)
                arcsec = 60.0_bp*(rarcmin-arcmin)
             ELSE
                CALL toInt(records(1)(60:61), arcmin, error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (40).", 1)
                   RETURN
                END IF
                CALL toReal(records(1)(63:67), arcsec, error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (45).", 1)
                   RETURN
                END IF
             END IF
             IF (records(1)(56:56) == "-" .AND. deg /= 0) THEN
                deg = -1 * deg
             ELSE IF (records(1)(56:56) == "-" .AND. deg == 0 .AND. arcmin /= 0) THEN
                arcmin = -1 * arcmin
             ELSE IF (records(1)(56:56) == "-" .AND. deg == 0 .AND. arcmin == 0) THEN
                arcsec = -1.0_bp * arcsec
             END IF
             ! Create observation object containing R.A. and Dec. + epoch:
             CALL NEW(obs_scoord, hour, min, sec, deg, arcmin, arcsec, t)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "TRACE BACK (95)", 1)
                RETURN
             END IF
          ELSE
             CALL toReal(records(1)(41:55), ra, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (50).", 1)
                RETURN
             END IF
             CALL toReal(records(1)(56:71), dec, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (55).", 1)
                RETURN
             END IF
             ! Create observation object containing R.A. and Dec. + epoch:
             CALL NEW(obs_scoord, longitude=ra, latitude=dec, t=t)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "TRACE BACK (100)", 1)
                RETURN
             END IF
          END IF
          IF (LEN_TRIM(records(1)(72:77)) /= 0) THEN
             CALL toReal(records(1)(72:77), mag, error)
          ELSE
             mag = 99.9_bp
          END IF
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "Could not convert string to number (57).", 1)
             RETURN
          END IF

          ! Uncertainties:
          stdev_ = 0.0_bp
          IF (PRESENT(stdev)) THEN
             stdev_ = stdev
             correlation = 0.0_bp
          ELSE
             IF (records(2)(43:43) == ".") THEN
                CALL toReal(records(2)(41:46), stdev_(2), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (60).", 1)
                   RETURN
                END IF
                CALL toReal(records(2)(48:53), stdev_(3), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (65).", 1)
                   RETURN
                END IF
                stdev_ = stdev_*rad_asec
             ELSE
                CALL toReal(records(2)(41:41) // "." // &
                     records(2)(42:48), stdev_(2), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (70).", 1)
                   RETURN
                END IF
                CALL toReal(records(2)(49:49) // "." // &
                     records(2)(50:56), stdev_(3), error)
                IF (error) THEN
                   CALL errorMessage("Observations / readObservationFile", &
                        "Could not convert string to number (75).", 1)
                   RETURN
                END IF
                ! Transform units to radians from 10**(-7) radians:
                stdev_ = stdev_*0.0000001_bp
             END IF
             CALL toReal(records(2)(57:64), correlation, error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (77).", 1)
                RETURN
             END IF
          END IF
          covariance = 0.0_bp
          covariance(2,2) = stdev_(2)**2.0_bp
          covariance(3,3) = stdev_(3)**2.0_bp
          covariance(2,3) = correlation*stdev_(2)*stdev_(3)
          covariance(3,2) = covariance(2,3)

          ! Compute the heliocentric position of the observer at epoch t:
          obsy_code = records(1)(129:132)
          IF (obsy_code(1:1) == "0") THEN
             obsy_code = obsy_code(2:4)
             obsy_code(4:4) = " "
          END IF
          obsy = getObservatory(obsies, obsy_code)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (105)", 1)
             RETURN
          END IF

          ! Create observation object:
          CALL NULLIFY(this%obs_arr(iobs))
          SELECT CASE (records(1)(21:21))
          CASE (" ", "C", "A", "T", "M", "P", "c")
             obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)
             CALL NEW(this%obs_arr(iobs), number=number, designation=designation, &
                  discovery=discovery, note1=records(1)(124:124), &
                  note2=records(1)(21:21), obs_scoord=obs_scoord, &
                  covariance=covariance, obs_mask=obs_mask, &
                  mag=mag, filter=records(1)(78:79), obsy=obsy, &
                  obsy_ccoord=obsy_ccoord)
          CASE ("S")
             IF (LEN_TRIM(records(3)) == 0) THEN
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "No satellite coordinates available in observation file.", 1)
                RETURN
             END IF
             CALL toReal(records(3)(42:59), position(1), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (80).", 1)
                RETURN
             END IF
             CALL toReal(records(3)(60:77), position(2), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (85).", 1)
                RETURN
             END IF
             CALL toReal(records(3)(78:95), position(3), error)
             IF (error) THEN
                CALL errorMessage("Observations / readObservationFile", &
                     "Could not convert string to number (90).", 1)
                RETURN
             END IF
             SELECT CASE (records(3)(41:41))
             CASE ("0")
                coord_unit = 1
                position(1:3) = position(1:3)/km_au
             CASE ("1")
                coord_unit = 2
             CASE default
                error = .TRUE.
                CALL errorMessage("Observations / readObservationFile", &
                     "Unit type "" // records(3)(41:41) // "" not available.", 1)
                RETURN
             END SELECT
             obsy_ccoord = getObservatoryCCoord(obsies, "500", t)
             CALL rotateToEquatorial(obsy_ccoord)
             coordinates = getCoordinates(obsy_ccoord)
             CALL NULLIFY(obsy_ccoord)
             coordinates(1:3) = coordinates(1:3) + position(1:3)
             CALL NEW(obsy_ccoord, coordinates, "equatorial", t)
             coordinates = 0.0_bp
             coordinates(1:3) = position(1:3)
             CALL NEW(satellite_ccoord, coordinates, "equatorial", t)
             CALL NEW(this%obs_arr(iobs), number=number, designation=designation, &
                  discovery=discovery, note1=records(1)(124:124), &
                  note2=records(1)(21:21), obs_scoord=obs_scoord, &
                  covariance=covariance, obs_mask=obs_mask, &
                  mag=mag, filter=records(1)(78:79), obsy=obsy, &
                  obsy_ccoord=obsy_ccoord, satellite_ccoord=satellite_ccoord, &
                  coord_unit=coord_unit)
          CASE ("R")
             iobs = iobs - 1
          CASE default
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "No such (" // records(1)(21:21) // ") option available.", 1)
             IF (err_verb >= 1) THEN
                WRITE(stderr,"(A)") records(1)
             END IF
             RETURN
          END SELECT
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (110)", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(geocenter_ccoord)
          CALL NULLIFY(obsy_ccoord)
          CALL NULLIFY(satellite_ccoord)

       END DO obsloop_mpc3

       i = iobs


    CASE ("gaia3")

       covariance = 0.0_bp
       position = 0.0_bp
       velocity = 0.0_bp
       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)

       ! Origin of observation dates: 1.0 Jan 2010 ( = JD 2455197.5 = MJD 55197.0)
       mjd_tcb = 55197.0_bp
       i = 0 
       DO

          READ(getUnit(obsf), "(A)", iostat=err) line
          IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file (1).", 1)
             RETURN
          ELSE IF (err < 0) THEN ! end-of-file
             EXIT
          ELSE IF (line(1:1) == "#") THEN
             CYCLE
          END IF
          READ(line, *, iostat=err) number, ccd, dt, position(2), &
               position(3), covariance(2,2), covariance(2,3), covariance(3,3), &
               coordinates(1:6)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file (2).", 1)
             RETURN
          END IF
          covariance(3,2) = covariance(2,3)

          i = i + 1
          IF (i == 1) THEN
             discovery = .TRUE.
          ELSE
             discovery = .FALSE.
          END IF

          CALL NEW(t, mjd_tcb + dt, "TCB")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (115)", 1)
             RETURN
          END IF
          obsy = getObservatory(obsies, "247")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (120)", 1)
             RETURN
          END IF
          ! Transform barycentric spacecraft coordinates to
          ! heliocentric spacecraft coordinates
          mjd_tt = getMJD(t, "TT")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (121)", 1)
             RETURN
          END IF
          planeph => JPL_ephemeris(mjd_tt, 12, 11, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (122)", 1)
             RETURN
          END IF
          CALL NEW(obsy_ccoord, coordinates + planeph(1,:), "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (125)", 1)
             RETURN
          END IF
          DEALLOCATE(planeph, stat=err)
          CALL NEW(obs_scoord, position, velocity, "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (130)", 1)
             RETURN
          END IF
          satellite_ccoord = getObservatoryCCoord(obsies, "500", t)
          CALL rotateToEquatorial(satellite_ccoord)
          CALL rotateToEquatorial(obsy_ccoord)
          coordinates = getCoordinates(obsy_ccoord) - getCoordinates(satellite_ccoord)
          CALL NULLIFY(satellite_ccoord)
          CALL NEW(satellite_ccoord, coordinates, "equatorial", t)
          CALL NULLIFY(this%obs_arr(i))
          CALL NEW(this%obs_arr(i), number=number, designation=" ", &
               discovery=discovery, note1=" ", note2="S", &
               obs_scoord=obs_scoord, covariance=covariance, &
               obs_mask=obs_mask, mag=99.9_bp, filter=" ", &
               obsy=obsy, obsy_ccoord=obsy_ccoord, &
               satellite_ccoord=satellite_ccoord, &
               coord_unit=2)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (135)", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)
          CALL NULLIFY(satellite_ccoord)

       END DO

    CASE ("geo")

       covariance = 0.0_bp
       position = 0.0_bp
       velocity = 0.0_bp
       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)

       i = 0
       DO

          READ(getUnit(obsf), "(A)", iostat=err) line
          IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file (1).", 1)
             RETURN
          ELSE IF (err < 0) THEN ! end-of-file
             EXIT
          ELSE IF (line(1:1) == "#") THEN
             CYCLE
          END IF
          READ(line, *, iostat=err) number, obsy_code, mjd_utc, position(2), &
               position(3), covariance(2,2), covariance(3,3), coordinates(1:6)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file (2).", 1)
             RETURN
          END IF

          i = i + 1
          IF (i == 1) THEN
             discovery = .TRUE.
          ELSE
             discovery = .FALSE.
          END IF

          CALL NEW(t, mjd_utc, "UTC")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (115)", 1)
             RETURN
          END IF
          obsy = getObservatory(obsies, "247")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (120)", 1)
             RETURN
          END IF
          ! Transform geocentric spacecraft coordinates to
          ! heliocentric spacecraft coordinates
          mjd_tt = getMJD(t, "TT")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (121)", 1)
             RETURN
          END IF
          planeph => JPL_ephemeris(mjd_tt, 3, 11, error)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (122)", 1)
             RETURN
          END IF
          CALL NEW(obsy_ccoord, coordinates + planeph(1,:), "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (125)", 1)
             RETURN
          END IF
          DEALLOCATE(planeph, stat=err)
          CALL NEW(obs_scoord, position, velocity, "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (130)", 1)
             RETURN
          END IF
          satellite_ccoord = getObservatoryCCoord(obsies, "500", t)
          CALL rotateToEquatorial(satellite_ccoord)
          CALL rotateToEquatorial(obsy_ccoord)
          coordinates = getCoordinates(obsy_ccoord) - getCoordinates(satellite_ccoord)
          CALL NULLIFY(satellite_ccoord)
          CALL NEW(satellite_ccoord, coordinates, "equatorial", t)
          CALL NULLIFY(this%obs_arr(i))
          CALL NEW(this%obs_arr(i), number=number, designation=" ", &
               discovery=discovery, note1=" ", note2="V", &
               obs_scoord=obs_scoord, covariance=covariance, &
               obs_mask=obs_mask, mag=99.9_bp, filter=" ", &
               obsy=obsy, obsy_ccoord=obsy_ccoord, &
               satellite_ccoord=satellite_ccoord, &
               coord_unit=2)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (135)", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)
          CALL NULLIFY(satellite_ccoord)

       END DO

    CASE ("strk")


       DEALLOCATE(this%obs_arr, this%obs_note_arr, stat=err)
       ALLOCATE(this%obs_arr(2*nlines), this%obs_note_arr(2*nlines), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readObservationFile", &
               "Could not allocate memory.", 1)
          RETURN
       END IF
       this%obs_note_arr = " "

       number = ""
       covariance = 0.0_bp
       position = 0.0_bp
       velocity = 0.0_bp
       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)

       i = 0
       DO

          READ(getUnit(obsf), "(A)", iostat=err) line
          IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file (1).", 1)
             RETURN
          ELSE IF (err < 0) THEN ! end-of-file
             EXIT
          ELSE IF (line(1:1) == "#") THEN
             CYCLE
          END IF
          READ(line, *, iostat=err) designation, mjd_utc, dt, &
               obsy_code, position(2), covariance(2,2), position(3), &
               covariance(3,3), covariance(2,3), strk_len, &
               strk_len_unc, strk_direction, strk_direction_unc, &
               mag, mag_unc, pos1(2), pos1(3), border1, pos2(2), &
               pos2(3), border2, strk_class
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file (2).", 1)
             RETURN
          END IF

          ! Change correlation to covariance
          covariance(2,3) = covariance(2,3) * covariance(2,2) * covariance(3,3)
          covariance(3,2) = covariance(2,3)
          ! Change standard deviations to variances
          covariance(2,2) = covariance(2,2)**2
          covariance(3,3) = covariance(3,3)**2
          ! Change degrees to radians
          position = position*rad_deg
          covariance = covariance*(rad_deg**2)
          pos1 = pos1*rad_deg
          pos2 = pos2*rad_deg
          strk_len = strk_len*rad_deg
          strk_len_unc = strk_len_unc*rad_deg
          strk_direction = strk_direction*rad_deg
          strk_direction_unc = strk_direction_unc*rad_deg

          ! Streak start point
          i = i + 1
          discovery = .FALSE.
          CALL NEW(t, mjd_utc-dt/2.0_bp, "UTC")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (115)", 1)
             RETURN
          END IF
          CALL NEW(obs_scoord, pos1, velocity, "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (130)", 1)
             RETURN
          END IF
          obsy = getObservatory(obsies, obsy_code)
          obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)
          CALL NULLIFY(this%obs_arr(i))
          CALL NEW(this%obs_arr(i), number=number, designation=designation, &
               discovery=discovery, note1=" ", note2="C", &
               obs_scoord=obs_scoord, covariance=covariance, &
               obs_mask=obs_mask, mag=mag, filter=" ", &
               obsy=obsy, obsy_ccoord=obsy_ccoord)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (135)", 1)
             RETURN
          END IF
          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)

          ! Streak end point
          i = i + 1
          discovery = .FALSE.
          CALL NEW(t, mjd_utc+dt/2.0_bp, "UTC")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (115)", 1)
             RETURN
          END IF
          CALL NEW(obs_scoord, pos2, velocity, "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (130)", 1)
             RETURN
          END IF
          obsy = getObservatory(obsies, obsy_code)
          obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)
          CALL NULLIFY(this%obs_arr(i))
          CALL NEW(this%obs_arr(i), number=number, designation=designation, &
               discovery=discovery, note1=" ", note2="C", &
               obs_scoord=obs_scoord, covariance=covariance, &
               obs_mask=obs_mask, mag=mag, filter=" ", &
               obsy=obsy, obsy_ccoord=obsy_ccoord)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (135)", 1)
             RETURN
          END IF
          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)

       END DO




    CASE ("sentinel")

       covariance = 0.0_bp
       position = 0.0_bp
       velocity = 0.0_bp
       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)

       ! Origin of observation dates: 1.0 Jan 2010 ( = JD 2455197.5 = MJD 55197.0)
       mjd_utc = 55197.0_bp
       i = 0
       DO

          READ(getUnit(obsf), "(A)", iostat=err) line
          IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file.", 1)
             RETURN
          ELSE IF (err < 0) THEN ! end-of-file
             EXIT
          ELSE IF (line(1:1) == "#") THEN
             CYCLE
          END IF
          READ(line, *, iostat=err) number, ccd, dt, position(2), &
               position(3), covariance(2,2), covariance(2,3), covariance(3,3), &
               coordinates(1:6)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Observations / readObservationFile", &
                  "Error while reading observations from file.", 1)
             RETURN
          END IF
          covariance(3,2) = covariance(2,3)

          i = i + 1
          IF (i == 1) THEN
             discovery = .TRUE.
          ELSE
             discovery = .FALSE.
          END IF

          CALL NEW(t, mjd_utc + dt, "utc")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (115)", 1)
             RETURN
          END IF
          obsy = getObservatory(obsies, "247")
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (120)", 1)
             RETURN
          END IF
          CALL NEW(obsy_ccoord, coordinates, "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (125)", 1)
             RETURN
          END IF
          CALL NEW(obs_scoord, position, velocity, "equatorial", t)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (130)", 1)
             RETURN
          END IF
          satellite_ccoord = getObservatoryCCoord(obsies, "500", t)
          CALL rotateToEquatorial(satellite_ccoord)
          CALL rotateToEquatorial(obsy_ccoord)
          coordinates = getCoordinates(obsy_ccoord) - getCoordinates(satellite_ccoord)
          CALL NULLIFY(satellite_ccoord)
          CALL NEW(satellite_ccoord, coordinates, "equatorial", t)
          CALL NULLIFY(this%obs_arr(i))
          CALL NEW(this%obs_arr(i), number=number, designation=" ", &
               discovery=discovery, note1=" ", note2="S", &
               obs_scoord=obs_scoord, covariance=covariance, &
               obs_mask=obs_mask, mag=99.9_bp, filter=" ", &
               obsy=obsy, obsy_ccoord=obsy_ccoord, &
               satellite_ccoord=satellite_ccoord, &
               coord_unit=2)
          IF (error) THEN
             CALL errorMessage("Observations / readObservationFile", &
                  "TRACE BACK (135)", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)
          CALL NULLIFY(satellite_ccoord)

       END DO

    CASE ("elgb")

       CALL readELGBFile(this, obsf, stdev)
       i = getNumberOfObservations(this)

    CASE default

       error = .TRUE.
       CALL errorMessage("Observations / readObservationFile", &
            "Cannot recognize the format of the following file: " // TRIM(fname),1)

    END SELECT

    CALL NULLIFY(obsies)
    IF (error) THEN
       CALL errorMessage("Observations / readObservationFile", &
            "TRACE BACK (140)", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%obs_arr) .AND. i > 0) THEN
       this%obs_arr => reallocate(this%obs_arr,i)
    END IF
    IF (ASSOCIATED(this%obs_note_arr) .AND. i > 0) THEN
       this%obs_note_arr => reallocate(this%obs_note_arr,i)
    END IF

  END SUBROUTINE readObservationFile





  !! *Description*:
  !!
  !! Intrisic routine that reads a file of elgb-format (used by Ted B.).
  !!
  !! Returns error.
  !!
  SUBROUTINE readELGBFile(this, elgb_file, stdev)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    TYPE (File), INTENT(in)            :: elgb_file
    REAL(bp), DIMENSION(6), INTENT(in) :: stdev
    TYPE (Observatories)               :: obsies
    TYPE (Observatory)                 :: obsy
    TYPE (Time)                        :: t
    TYPE (SphericalCoordinates)        :: obs_scoord
    TYPE (CartesianCoordinates)        :: obsy_ccoord, &
         geocenter_ccoord, satellite_ccoord
    CHARACTER(len=DESIGNATION_LEN)     :: number
    CHARACTER(len=96)                  :: line1
    CHARACTER(len=3)                   :: timescale, obsy_code
    REAL(bp), DIMENSION(6,6)           :: covariance
    REAL(bp)                           :: day, sec, arcsec, mag
    INTEGER                            :: i, err, year, &
         month, hour, min, deg, arcmin, nlines
    LOGICAL, DIMENSION(6)              :: obs_mask
    LOGICAL                            :: discovery

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readELGBFile", &
            "Astrometric observations (MPC) have already been initialized.", 1)
       RETURN
    END IF

    nlines = getNrOfLines(elgb_file)
    IF (error) THEN
       CALL errorMessage("Observations / readELGBFile", &
            "TRACE BACK 4", 1)
       RETURN
    END IF

    ALLOCATE(this%obs_arr(nlines), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readELGBFile", &
            "Could not allocate memory for the observations.", 1)
       RETURN
    END IF

    CALL NEW(obsies)
    IF (error) THEN
       CALL errorMessage("Observations / readELGBFile", &
            "TRACE BACK 5", 1)
       RETURN
    END IF

    obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)

    ! This is the file format (http://asteroid.lowell.edu/asteroid/loneos/public_obs.html):
    ! to be read from the observations file (motions are ignored):
    !
    !         Date              R.A.                  Dec.          Mag.&Color
    !  --------17------- --------13--------- ----------14---------- -----8----      
    ! "(I4,1X,I2,2X,F8.5,2X,I2,1X,I2,1X,F6.3,2X,A1,I2,1X,I2,1X,F5.2,2X,F4.1,A1,
    !
    ! LONEOS Discovery  Obs.           LONEOS RA motion Dec motion
    !   ID   Asterisk   Code  Detector Region  s/day      am/day
    ! --9---  ----1---- --4-- ---3---- --6--- ---14---- ----6-----
    ! 1X,A8,    A1,     1X,A3, 2X,A1,    A6,   7X,F7.1,   F6.2)"

    i = 0
    DO

       READ(getUnit(elgb_file), "(A75)", iostat=err) line1
       IF (err > 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readELGBFile", &
               "Error while reading observations from file.", 1)
          RETURN
       ELSE IF (err < 0) THEN ! end-of-file
          EXIT
       END IF

       ! If the line is empty or marked with "#", read the next line:
       IF (LEN_TRIM(line1) == 0 .OR. line1(1:1) == "#") THEN
          CYCLE
       END IF

       i = i + 1
       !       CALL toInt(line1(1:5), nr, error)
       number = ""
       IF (line1(62:62) == "*") THEN
          discovery = .TRUE.
       ELSE
          discovery = .FALSE.
       END IF
       CALL toInt(line1(1:4), year, error)
       CALL toInt(line1(6:7), month, error)
       CALL toReal(line1(10:17), day, error)
       CALL toInt(line1(20:21), hour, error)
       CALL toInt(line1(23:24), min, error)
       CALL toReal(line1(26:31), sec, error)
       CALL toInt(line1(35:36), deg, error)
       CALL toInt(line1(38:39), arcmin, error)
       CALL toReal(line1(41:45), arcsec, error)
       IF (LEN_TRIM(line1(48:51)) /= 0) THEN
          CALL toReal(line1(48:51), mag, error)
       ELSE
          mag = 99.9_bp
       END IF
       IF (error) THEN
          CALL errorMessage("Observations / readELGBFile", &
               "Could not convert string to number.", 1)
          RETURN
       END IF

       ! R.A. and Dec. + epoch:
       !
       ! The timescale used by MPC prior to year 1972 is UT1 and since
       ! then UTC has been used:  
       IF (year < 1972) THEN
          timescale = "UT1"
       ELSE
          timescale = "UTC"
       END IF
       ! Create epoch object:
       CALL NEW(t, year, month, day, timescale)
       IF (error) THEN
          CALL errorMessage("Observations / readELGBFile", &
               "TRACE BACK 17", 1)
          RETURN
       END IF
       IF (line1(34:34) == "-" .AND. deg /= 0) THEN
          deg = -1 * deg
       ELSE IF (line1(34:34) == "-" .AND. deg == 0 .AND. arcmin /= 0) THEN
          arcmin = -1 * arcmin
       ELSE IF (line1(34:34) == "-" .AND. deg == 0 .AND. arcmin == 0) THEN
          arcsec = -1.0_bp * arcsec
       END IF
       ! Create observation object containing R.A. and Dec. + epoch:
       CALL NEW(obs_scoord, hour, min, sec, deg, arcmin, arcsec, t)
       IF (error) THEN
          CALL errorMessage("Observations / readELGBFile", &
               "TRACE BACK 18", 1)
          RETURN
       END IF
       covariance = 0.0_bp
       covariance(2,2) = stdev(2)**2.0_bp
       covariance(3,3) = stdev(3)**2.0_bp

       ! Compute the heliocentric position of the observer at epoch t:
       obsy_code = line1(64:66)
       SELECT CASE (obsy_code)
       CASE (" -3")
          obsy_code = "500"
       END SELECT
       obsy = getObservatory(obsies, obsy_code)
       IF (error) THEN
          CALL errorMessage("Observations / readELGBFile", &
               "Could not convert observatory code string to integer.", 1)
          RETURN
       END IF

       ! Create observation object:
       CALL NULLIFY(this%obs_arr(i))
       obsy_ccoord = getObservatoryCCoord(obsies, obsy_code, t)

       CALL NEW(this%obs_arr(i), number=number, &
            designation=line1(54:61), discovery=discovery, & 
            note1=" ", note2=line1(69:69), obs_scoord=obs_scoord, &
            covariance=covariance, obs_mask=obs_mask, mag=mag, &
            filter=line1(52:52), obsy=obsy, obsy_ccoord=obsy_ccoord)

       IF (error) THEN
          CALL errorMessage("Observations / readELGBFile", &
               "TRACE BACK 20", 1)
          RETURN
       END IF

       CALL NULLIFY(obs_scoord)
       CALL NULLIFY(t)
       CALL NULLIFY(obsy)
       CALL NULLIFY(geocenter_ccoord)
       CALL NULLIFY(obsy_ccoord)
       CALL NULLIFY(satellite_ccoord)

    END DO

    CALL NULLIFY(obsies)
    IF (error) THEN
       CALL errorMessage("Observations / readELGBFile", &
            "TRACE BACK 21", 1)
       RETURN
    END IF

    this%obs_arr => reallocate(this%obs_arr,i)

  END SUBROUTINE readELGBFile





  !! *Description*:
  !!
  !!
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE readGAIAFile(this, gaia_file)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)    :: this
    TYPE (File), INTENT(in)               :: gaia_file
    TYPE (Observatories)                  :: obsies
    TYPE (Observatory)                    :: obsy
    TYPE (Time)                           :: t
    TYPE (SphericalCoordinates)           :: obs_scoord
    TYPE (CartesianCoordinates)           :: obsy_ccoord
    CHARACTER(len=DESIGNATION_LEN)        :: number
    CHARACTER(len=93)                     :: line, frst_line
    CHARACTER(len=14)                     :: empty_form, full_form
    CHARACTER(len=13)                     :: designation
    CHARACTER(len=5)                      :: form_, line_length, nobs_char
    CHARACTER(len=3), PARAMETER           :: rec_length="14"
    INTEGER, PARAMETER                    :: nrecords=17
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: stdevs
    REAL(bp), DIMENSION(6,6)              :: covariance
    REAL(bp), DIMENSION(:), ALLOCATABLE   :: dates, magnitudes, &
         distances, longitudes, latitudes, rot_angles
    REAL(bp), DIMENSION(3)                :: position, velocity
    REAL(bp)                              :: mjd_utc
    INTEGER                               :: i, err, nlines, nobs, nobs_, len
    LOGICAL, DIMENSION(6)                 :: obs_mask
    LOGICAL                               :: discovery

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readGAIAFile", &
            "This object has already been initialized.", 1)
       RETURN
    END IF

    nlines = getNrOfLines(gaia_file)
    IF (error) THEN
       CALL errorMessage("Observations / readGAIAFile", &
            "TRACE BACK 2", 1)
       RETURN
    END IF

    IF (MOD(nlines,nrecords) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / readGAIAFile", &
            "Number of records not what's expected for GAIA data.", 1)
       RETURN
    END IF

    CALL NEW(obsies)
    IF (error) THEN
       CALL errorMessage("Observations / readGAIAFile", &
            "TRACE BACK 5", 1)
       RETURN
    END IF

    nobs = 0
    DO

       READ(getUnit(gaia_file), "(A36)", iostat=err) frst_line
       IF (err > 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 4", 1)
          RETURN
       ELSE IF (err < 0) THEN
          EXIT
       END IF

       ! Maximum number of observations is assumed to be 9999
       CALL toInt(frst_line(33:36), nobs_, error)
       IF (error) THEN
          CALL errorMessage("Observations / readGAIAFile", &
               "Could not convert character string to integer (1).", 1)
          RETURN       
       END IF
       designation = frst_line(17:29)
       number = frst_line(11:15)

       nobs = nobs + nobs_
       this%obs_arr => reallocate(this%obs_arr, nobs)
       !IF (err /= 0) THEN
       !   error = .TRUE.
       !   CALL errorMessage("Observations / readGAIAFile", &
       !        "Could not allocate memory for the observations.", 1)
       !   RETURN
       !END IF

       ALLOCATE(dates(nobs_),magnitudes(nobs_),distances(nobs_),&
            longitudes(nobs_),latitudes(nobs_),rot_angles(nobs_), &
            stdevs(nobs_,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "Could not allocate memory.", 1)
          RETURN
       END IF

       len = 10*nobs_+1
       CALL toString(len, form_, error)
       IF (error) THEN
          CALL errorMessage("Observations / readGAIAFile", &
               "Could not convert integer to character string.", 1)
          RETURN       
       END IF

       line_length = ADJUSTL(form_)
       empty_form = "(A"//TRIM(line_length)//")"

       nobs_char = ADJUSTL(frst_line(33:36))
       full_form = "("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".4)"

       ! Which field of view
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 10", 1)
          RETURN
       END IF

       ! Observation date
       READ(getUnit(gaia_file), full_form, iostat=err) dates
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 15", 1)
          RETURN
       END IF

       ! Transverse ordinate in degrees
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 20", 1)
          RETURN
       END IF

       ! Apparent magnitudes
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".1)"
       READ(getUnit(gaia_file), full_form, iostat=err) magnitudes
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 25", 1)
          RETURN
       END IF
       WHERE (ABS(magnitudes) < 10.0_bp*EPSILON(magnitudes)) magnitudes = 99.9_bp

       ! Angular diameter in mas
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 30", 1)
          RETURN
       END IF

       ! Distance to earth in au
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".3)"
       READ(getUnit(gaia_file), full_form, iostat=err) distances
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 35", 1)
          RETURN
       END IF

       ! Radial velocity to earth in km/s
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 40", 1)
          RETURN
       END IF
       ! Distance to sun in au   
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 45", 1)
          RETURN
       END IF
       ! Phase angle in degrees
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 50", 1)
          RETURN
       END IF
       ! Inertial speed along-scan in mas/s
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 55", 1)
          RETURN
       END IF
       ! Inertial speed across-scan in mas/s
       READ(getUnit(gaia_file), empty_form, iostat=err) line
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 60", 1)
          RETURN
       END IF

       ! Orientation of the scan (deg)
       !    READ(getUnit(gaia_file), empty_form, iostat=err) line
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".2)"
       READ(getUnit(gaia_file), full_form, iostat=err) rot_angles
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 65", 1)
          RETURN
       END IF
       rot_angles = rot_angles*rad_deg

       ! Ecliptic longitude (deg)
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".9)"
       READ(getUnit(gaia_file), full_form, iostat=err) longitudes
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 70", 1)
          RETURN
       END IF
       longitudes = longitudes*rad_deg

       ! Ecliptic latitude  (deg)
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".9)"
       READ(getUnit(gaia_file), full_form, iostat=err) latitudes
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 75", 1)
          RETURN
       END IF
       latitudes = latitudes*rad_deg

       ! Along scan stdev (milli-as)
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".2)"
       READ(getUnit(gaia_file), full_form, iostat=err) stdevs(:,1)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 80", 1)
          RETURN
       END IF
       ! Across scan stdev (milli-as)
       full_form="("//TRIM(nobs_char)//"F"//TRIM(rec_length)//".2)"
       READ(getUnit(gaia_file), full_form, iostat=err) stdevs(:,2)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / readGAIAFile", &
               "TRACE BACK 85", 1)
          RETURN
       END IF
       stdevs = stdevs*1.0e-3_bp*rad_asec


       ! True rotation angle is: 
       ! direction of the scan w.r.t. equator = (direction of the scan
       ! w.r.t. ecliptic pole) - (obliquity of the ecliptic)
       rot_angles = rot_angles - eps

       ! Do not include the distances
       distances = 0.0_bp

       obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)

       DO i=1,nobs_
          discovery = .FALSE.
          ! Origin of dates: 3.5 Jan 2010, JD 2455200.d0
          mjd_utc = 55199.5_bp
          CALL NEW(t, mjd_utc+dates(i), "utc")
          IF (error) THEN
             CALL errorMessage("Observations / readGAIAFile", &
                  "TRACE BACK 90", 1)
             RETURN
          END IF
          obsy = getObservatory(obsies, "500")
          IF (error) THEN
             CALL errorMessage("Observations / readGaiaFile", &
                  "TRACE BACK 95", 1)
             RETURN
          END IF
          obsy_ccoord = getObservatoryCCoord(obsies, "500", t)
          IF (error) THEN
             CALL errorMessage("Observations / readGaiaFile", &
                  "TRACE BACK 100", 1)
             RETURN
          END IF
          ! Note: ECLIPTIC coordinates
          position = (/ distances(i), longitudes(i), latitudes(i) /)
          velocity = 0.0_bp
          CALL NEW(obs_scoord, position, velocity, "ecliptic", t)
          CALL rotateToEquatorial(obs_scoord)
          IF (error) THEN
             CALL errorMessage("Observations / readGAIAFile", &
                  "TRACE BACK 105", 1)
             RETURN
          END IF

          covariance = 0.0_bp
          covariance(2,2) = (stdevs(i,1)*COS(rot_angles(i)))**2.0_bp + &
               (stdevs(i,2)*SIN(rot_angles(i)))**2.0_bp
          covariance(2,3) = (stdevs(i,1)**2.0_bp - stdevs(i,2)**2.0_bp)* &
               COS(rot_angles(i))*SIN(rot_angles(i))
          covariance(3,2) = (stdevs(i,1)**2.0_bp - stdevs(i,2)**2.0_bp)* &
               COS(rot_angles(i))*SIN(rot_angles(i))
          covariance(3,3) = (stdevs(i,1)*SIN(rot_angles(i)))**2.0_bp + &
               (stdevs(i,2)*COS(rot_angles(i)))**2.0_bp

          CALL NULLIFY(this%obs_arr(nobs-nobs_+i))
          CALL NEW(this%obs_arr(nobs-nobs_+i), number=number, &
               designation=designation, discovery=discovery, & 
               note1=" ", note2="C", obs_scoord=obs_scoord, &
               covariance=covariance, obs_mask=obs_mask, &
               mag=magnitudes(i), filter=" ", obsy=obsy, &
               obsy_ccoord=obsy_ccoord)
          IF (error) THEN
             CALL errorMessage("Observations / readGAIAFile", &
                  "TRACE BACK 110", 1)
             RETURN
          END IF

          CALL NULLIFY(obs_scoord)
          CALL NULLIFY(t)
          CALL NULLIFY(obsy)
          CALL NULLIFY(obsy_ccoord)

       END DO

       DEALLOCATE(dates, magnitudes, distances, longitudes, latitudes, &
            rot_angles, stdevs, stat=err)

    END DO

    this%obs_arr => reallocate(this%obs_arr,nobs)
    CALL NULLIFY(obsies)
    IF (error) THEN
       CALL errorMessage("Observations / readGaiaFile", &
            "TRACE BACK 115", 1)
       RETURN
    END IF


  END SUBROUTINE readGAIAFile




  !! *Description*:
  !!
  !! Reallocates a pointer array of Observations-objects and
  !! copies the existing (not-nullified) data from the old array to
  !! the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myobservations => reallocate(myobservations,4)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_a1_Obss(array,n)

    IMPLICIT NONE
    TYPE (Observations), DIMENSION(:), POINTER :: reallocate_a1_Obss, array
    INTEGER, INTENT(in)                        :: n
    INTEGER                                    :: i, nold, err, nexist

    nexist = 0
    nold = SIZE(array,dim=1)
    ALLOCATE(reallocate_a1_Obss(n), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_a1_Obss => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    DO i=1, MIN(n,nold)
       IF (.NOT.array(i)%is_initialized) THEN
          CYCLE
       END IF
       nexist = nexist + 1
       reallocate_a1_Obss(nexist) = copy(array(i))
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_a1_Obss





  !! *Description*:
  !!
  !! 
  !!
  !! Returns error.
  !!
  SUBROUTINE setCovarianceMatrices(this, covariance)

    TYPE (Observations), INTENT(inout)   :: this
    REAL(bp), DIMENSION(6,6), INTENT(in) :: covariance
    INTEGER                              :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / setCovarianceMatrices", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1,this%nobs
       CALL setCovarianceMatrix(this%obs_arr(i), covariance)
       IF (error) THEN
          CALL errorMessage("Observations / setCovarianceMatrices", &
               "TRACE BACK", 1)
          RETURN
       END IF
    END DO

  END SUBROUTINE setCovarianceMatrices





  !! *Description*:
  !!
  !! 
  !!
  !! Returns error.
  !!
  SUBROUTINE setNumber_Obss(this, number)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout) :: this
    INTEGER, INTENT(in)                :: number
    INTEGER                            :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / setNumber", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1,SIZE(this%obs_arr,dim=1)
       CALL setNumber(this%obs_arr(i), number)
       IF (error) THEN
          CALL errorMessage("Observations / setNumber", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    END DO

    this%nobjects = 0
    CALL sortObservations(this)
    IF (error) THEN
       CALL errorMessage("Observations / setNumber", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

  END SUBROUTINE setNumber_Obss





  !! *Description*:
  !!
  !! Internal routine. Searches for different objects (different (1)
  !! number or (2) designation) and updates the total number of 
  !! observations. Updates the index vector, which is used to
  !! sort the observations in order of ascending epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE sortObservations(this, primary_sort, force_full)

    IMPLICIT NONE
    TYPE (Observations), INTENT(inout)     :: this
    CHARACTER(len=*), INTENT(in), OPTIONAL :: primary_sort
    LOGICAL, INTENT(in), OPTIONAL          :: force_full
    TYPE (Time)                            :: t
    CHARACTER(len=DESIGNATION_LEN)         :: number, id
    CHARACTER(len=16)                      :: primarysort
    REAL(bp), DIMENSION(:), ALLOCATABLE    :: tmp
    REAL(bp)                               :: mjd_tdt
    INTEGER                                :: i, iobs, iobj

    id = " "
    IF (this%nobjects == 0) THEN
       ! Since the number of different objects is 0, 
       ! all observations must be scanned:
       this%nobs = 0
    END IF

    ! Force full sort if requested:
    IF (PRESENT(force_full)) THEN
       IF (force_full) THEN
          this%nobjects = 0
          this%nobs = 0
       END IF
    END IF

    IF (PRESENT(primary_sort)) THEN
       primarysort = primary_sort
       CALL locase(primarysort, error)
       IF (error) THEN
          CALL errorMessage("Observations / sortObservations", &
               "The primary sort string contains forbidden characters.", 1)
          RETURN
       END IF
       IF (primarysort /= "designation" .AND. &
            primarysort /= "number") THEN
          error = .TRUE.
          CALL errorMessage("Observations / sortObservations", &
               "Cannot sort by " // TRIM(primarysort) // ".", 1)
          RETURN
       END IF
    ELSE
       ! Default:
       primarysort = "number"
    END IF

    ! Find total number of different objects:
    DO iobs=this%nobs+1, SIZE(this%obs_arr)

       number = getNumber(this%obs_arr(iobs))
       IF (error) THEN
          CALL errorMessage("Observations / sortObservations", &
               "TRACE BACK 1", 1)
          RETURN
       END IF

       ! Use (1) number or (2) designation:
       IF (LEN_TRIM(number) /= 0 .AND. primarysort == "number") THEN
          id = number
!!$          DO WHILE (LEN_TRIM(id) < 7)
!!$             id = "0" // TRIM(id)
!!$          END DO
       ELSE
          id = getDesignation(this%obs_arr(iobs))
          IF (error) THEN
             CALL errorMessage("Observations / sortObservations", &
                  "TRACE BACK 2", 1)
             RETURN
          END IF
       END IF

       ! ...look for it among already registered objects and
       ! exit when you find it or when you have searched through all... 
       iobj = this%nobjects
       DO WHILE (iobj>0)
          IF (this%objects(iobj) == TRIM(id)) EXIT
          iobj = iobj - 1
       END DO

       ! ...and if turns out that you did not find it, put it on the list:
       IF (iobj == 0) THEN
          this%nobjects = this%nobjects + 1
          this%objects => reallocate(this%objects,this%nobjects)
          this%objects(this%nobjects) = TRIM(id)
       END IF

    END DO

    IF (ASSOCIATED(this%ind)) THEN
       DEALLOCATE(this%ind)
    END IF
    ALLOCATE(this%ind(SIZE(this%obs_arr,dim=1)))
    IF (ASSOCIATED(this%criteria)) THEN
       ALLOCATE(tmp(SIZE(this%criteria,dim=1)))
       tmp = this%criteria
       DEALLOCATE(this%criteria)
    END IF
    ALLOCATE(this%criteria(SIZE(this%obs_arr,dim=1)))
    IF (ALLOCATED(tmp)) THEN
       this%criteria(1:SIZE(tmp)) = tmp
       DEALLOCATE(tmp)
    END IF

    ! Update ascending epoch index vector:
    DO i=this%nobs+1,SIZE(this%obs_arr)
       t = getTime(this%obs_arr(i))
       mjd_tdt = getMJD(t, "tdt")
       IF (error) THEN
          CALL errorMessage("Observations / sortObservations", &
               "TRACE BACK 4", 1)
          RETURN
       END IF
       this%criteria(i) = mjd_tdt
       CALL NULLIFY(t)
    END DO
    CALL quicksort(this%criteria, this%ind, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / sortObservations", &
            "Could not quicksort observations " // &
            TRIM(errstr), 1)
       RETURN
    END IF

    ! Update number of observations:
    this%nobs = SIZE(this%obs_arr,dim=1)

  END SUBROUTINE sortObservations





  !! *Description*:
  !!
  !! Writes the observations in specified format to the given logical unit.
  !!
  !! Returns error.
  !!
  SUBROUTINE writeObservationFile(this, lu, frmt, number)

    IMPLICIT NONE
    TYPE (Observations), INTENT(in) :: this
    CHARACTER(len=*)                :: frmt
    INTEGER, INTENT(in)             :: lu         
    CHARACTER(len=*), INTENT(in), OPTIONAL :: number

    CHARACTER(len=OBS_RECORD_LEN), DIMENSION(:), POINTER :: &
         records => NULL()
    CHARACTER(len=2)                :: note
    INTEGER                         :: i, j, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observations / writeObservationFile", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%nobs == 0) THEN
       error = .TRUE.
       CALL errorMessage("Observations / writeObservationFile", &
            "Object does not contain any observations.", 1)
       RETURN
    END IF

    DO i=1,this%nobs

       ! Get formatted observation records:
       IF (PRESENT(number)) THEN
          records => getObservationRecords(this%obs_arr(this%ind(i)), TRIM(frmt), number=number)
       ELSE
          records => getObservationRecords(this%obs_arr(this%ind(i)), TRIM(frmt))
       END IF
       IF (error) THEN
          CALL errorMessage("Observations / writeObservationFile", &
               "TRACE BACK 5", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF

       ! Write a different portion of the records depending on the
       ! observation format:
       SELECT CASE (TRIM(frmt))
       CASE ("mpc")
          DO j=1,SIZE(records,dim=1)
             IF (LEN_TRIM(this%obs_note_arr(this%ind(i))) == 0) THEN
                WRITE(lu, "(A80)", iostat=err) records(j)(1:80)
             ELSE IF (LEN_TRIM(this%obs_note_arr(this%ind(i))) /= 0 &
                  .AND. j == 1) THEN
                note = "!" // TRIM(this%obs_note_arr(this%ind(i)))
                WRITE(lu, "(A82)", iostat=err) note // records(j)(1:80)
             END IF
          END DO
       CASE ("mpc2")
          DO j=1,SIZE(records,dim=1)
             WRITE(lu, "(A132)", iostat=err) records(j)(1:132)
          END DO
       CASE ("mpc3")
          DO j=1,SIZE(records,dim=1)
             WRITE(lu, "(A132)", iostat=err) records(j)(1:132)
          END DO
       CASE ("elgb")
          DO j=1,SIZE(records,dim=1)
             WRITE(lu, "(A75)", iostat=err) records(j)(1:75)
          END DO
       CASE ("des")
          DO j=1,SIZE(records,dim=1)
             WRITE(lu, "(A)", iostat=err) TRIM(records(j)(:))
          END DO
       CASE default
          error = .TRUE.
          CALL errorMessage("Observations / writeObservationFile", &
               "Unknown observation format: " // TRIM(frmt), 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END SELECT
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / writeObservationFile", &
               "Error while writing records.", 1)
          DEALLOCATE(records, stat=err)
          RETURN
       END IF

       DEALLOCATE(records, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Observations / writeObservationFile", &
               "Could not deallocate memory.", 1)
          RETURN
       END IF

    END DO

  END SUBROUTINE writeObservationFile





END MODULE Observations_cl
