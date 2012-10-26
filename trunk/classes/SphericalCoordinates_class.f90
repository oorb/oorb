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
!! Type and routines for spherical coordinates.
!!
!! @see Observations_class
!! 
!! @author  MG
!! @version 2012-10-26
!!
MODULE SphericalCoordinates_cl

  USE Base_cl
  USE Time_cl

  USE random
  USE utilities
  USE linal

  IMPLICIT NONE

  PRIVATE :: new_SC
  PRIVATE :: new_SC_coordinates
  PRIVATE :: new_SC_hour
  PRIVATE :: new_SC_hourAndDistance
  PRIVATE :: new_SC_rad
  PRIVATE :: new_SC_radAndDistance
  PRIVATE :: new_SC_posAndVel
  PRIVATE :: nullify_SC
  PRIVATE :: copy_SC
  PRIVATE :: exist_SC
  PRIVATE :: addMultinormalDeviate_SC
  PRIVATE :: addUniformDeviate_SC
  PRIVATE :: checkAngles
  PRIVATE :: equal_SC
  PRIVATE :: getPosition_SC
  PRIVATE :: getVelocity_SC
  PRIVATE :: getFrame_SC
  PRIVATE :: getTime_SC
  PRIVATE :: rotateToEquatorial_SC
  PRIVATE :: rotateToEcliptic_SC

  TYPE SphericalCoordinates
     PRIVATE
     REAL(bp), DIMENSION(3)   :: position ! distance   [AU]
     !                                      longitude  [rad]
     !                                      latitude   [rad]
     REAL(bp), DIMENSION(3)   :: velocity ! ddistance  [AU/d]
     !                                      dlongitude [rad/d]
     !                                      dlatitude  [rad/d]
     CHARACTER(len=FRAME_LEN) :: frame = "equatorial"
     TYPE (Time)              :: t
     LOGICAL                  :: is_initialized = .FALSE.
  END TYPE SphericalCoordinates

  INTERFACE NEW
     MODULE PROCEDURE new_SC
     MODULE PROCEDURE new_SC_coordinates
     MODULE PROCEDURE new_SC_hour
     MODULE PROCEDURE new_SC_hourAndDistance
     MODULE PROCEDURE new_SC_rad
     MODULE PROCEDURE new_SC_radAndDistance
     MODULE PROCEDURE new_SC_posAndVel
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_SC
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_SC
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_SC
  END INTERFACE exist

  INTERFACE addMultinormalDeviate
     MODULE PROCEDURE addMultinormalDeviate_SC
  END INTERFACE addMultinormalDeviate

  INTERFACE addUniformDeviate
     MODULE PROCEDURE addUniformDeviate_SC
  END INTERFACE addUniformDeviate

  INTERFACE equal
     MODULE PROCEDURE equal_SC
  END INTERFACE equal

  INTERFACE getCoordinates
     MODULE PROCEDURE getCoordinates_SC
  END INTERFACE getCoordinates

  INTERFACE getPosition
     MODULE PROCEDURE getPosition_SC
  END INTERFACE getPosition

  INTERFACE getVelocity
     MODULE PROCEDURE getVelocity_SC
  END INTERFACE getVelocity

  INTERFACE getFrame
     MODULE PROCEDURE getFrame_SC
  END INTERFACE getFrame

  INTERFACE getTime
     MODULE PROCEDURE getTime_SC
  END INTERFACE getTime

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_SC_1
     MODULE PROCEDURE reallocate_SC_2
  END INTERFACE reallocate

  INTERFACE rotateToEquatorial
     MODULE PROCEDURE rotateToEquatorial_SC
  END INTERFACE rotateToEquatorial

  INTERFACE rotateToEcliptic
     MODULE PROCEDURE rotateToEcliptic_SC
  END INTERFACE rotateToEcliptic


CONTAINS




  !! *Description*:
  !!
  !! Initializes object using default values:
  !!   - position and velocity vectors are zero
  !!   - coordinate frame is equatorial
  !!   - epoch is the default Time object 
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%position       = 0.0_bp
    this%velocity       = 0.0_bp
    this%frame          = "equatorial"
    CALL NEW(this%t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_SC





  !! *Description*:
  !!
  !! Initializes object using values given as spherical coordinates
  !! and their time derivatives, the coordinate frame, and the epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC_coordinates(this, coordinates, frame, t)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6), INTENT(in)         :: coordinates
    CHARACTER(len=*), INTENT(in)               :: frame
    TYPE (Time), INTENT(in)                    :: t
    CHARACTER(len=FRAME_LEN)                   :: frame_

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%position       = coordinates(1:3)
    this%velocity       = coordinates(4:6)
    frame_ = frame
    CALL locase(frame_, error)
    IF (error) THEN
       CALL errorMessage("SphericalCoordinates / new", &
            "The frame string contains forbidden characters.", 1)
       RETURN
    END IF
    this%frame          = TRIM(frame_)
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_SC_coordinates





  !! *Description*:
  !!
  !! Initializes object using values given as the traditional
  !! astrometric coordinates, and a Time object. Coordinate 
  !! frame is equatorial, and the distance and the velocity 
  !! vector are set to zero.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC_hour(this, h, min, sec, deg, arcmin, arcsec, t)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), INTENT(in)                 :: sec, arcsec
    INTEGER, INTENT(in)                  :: h, min, deg, arcmin
    TYPE (Time), INTENT(in)              :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    ! Distance:
    this%position(1)    = 0.0_bp

    ! Transform longitude:
    ! Transform to hours:
    this%position(2)    = h + min/60.0_bp + sec/3600.0_bp
    ! Transform to radians:
    this%position(2)    = this%position(2) * rad_hour

    ! Transform latitude:
    ! Transform to degrees:
    IF (deg /= 0) THEN
       this%position(3) = SIGN((ABS(deg) + arcmin/60.0_bp + arcsec/3600.0_bp), REAL(deg, bp))
    ELSE IF (deg == 0 .AND. arcmin /= 0) THEN
       this%position(3) = SIGN((ABS(arcmin)/60.0_bp + arcsec/3600.0_bp), REAL(arcmin, bp))
    ELSE IF (deg == 0 .AND. arcmin == 0) THEN
       this%position(3) = arcsec/3600.0_bp
    END IF
    ! Transform to radians:
    this%position(3)    = this%position(3) * rad_deg

    CALL checkAngles(this)

    this%velocity       = 0.0_bp
    this%frame          = "equatorial"
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_SC_hour





  !! *Description*:
  !!
  !! Initializes object using values given as the traditional 
  !! astrometric coordinates, the distance, and a Time object. 
  !! Coordinate frame is equatorial, and the velocity vector 
  !! is set to zero.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC_hourAndDistance(this, &
       distance, h, min, sec, deg, arcmin, arcsec, t)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), INTENT(in)                 :: distance, sec, arcsec
    INTEGER, INTENT(in)                  :: h, min, deg, arcmin    
    TYPE (Time), INTENT(in)              :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    CALL NEW(this, h, min, sec, deg, arcmin, arcsec, t)
    IF (error) THEN
       CALL errorMessage("SphericalCoordinates / new", &
            "TRACE BACK", 1)
       RETURN
    END IF

    ! Distance:
    IF (distance < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Negative distance is not allowed.", 1)
       this%is_initialized = .FALSE.
       RETURN
    END IF
    this%position(1)    = distance

  END SUBROUTINE new_SC_hourAndDistance





  !! *Description*:
  !!
  !! Initializes object using values given as the longitude,
  !! latitude, and a Time object. Angles are given in radians.
  !! Coordinate frame is equatorial, and the distance and the 
  !! velocity vector are set to zero.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC_rad(this, longitude, latitude, t)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), INTENT(in)                 :: longitude, latitude
    TYPE (Time), INTENT(in)              :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%position(1)    = 0.0_bp
    this%position(2)    = longitude
    this%position(3)    = latitude
    this%velocity       = 0.0_bp
    CALL checkAngles(this)
    this%frame          = "equatorial"
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_SC_rad





  !! *Description*:
  !!
  !! Initializes object using values given as the distance [AU], 
  !! longitude, latitude, and a Time object. Angles are given in 
  !! radians. Coordinate frame is equatorial, and the 
  !! velocity vector is set to zero.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC_radAndDistance(this, distance, longitude, &
       latitude, t)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), INTENT(in)                 :: distance, longitude, latitude
    TYPE (Time), INTENT(in)              :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (distance < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Negative distance is not allowed.", 1)
       RETURN
    END IF
    this%position(1)    = distance
    this%position(2)    = longitude
    this%position(3)    = latitude
    this%velocity       = 0.0_bp
    CALL checkAngles(this)
    this%frame          = "equatorial"
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_SC_radAndDistance





  !! *Description*:
  !!
  !! Initializes object using values given as the position and 
  !! velocity vectors, the coordinate frame, and a Time object
  !! for the epoch. 
  !!
  !! The vectors contain:
  !!    - position(1): distance   [AU]
  !!    - position(2): longitude  [rad]
  !!    - position(3): latitude   [rad]
  !!    - velocity(1): ddistance  [AU/d]
  !!    - velocity(2): dlongitude [rad/d]
  !!    - velocity(3): dlatitude  [rad/d]
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SC_posAndVel(this, position, velocity, frame, t)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(3), INTENT(in)   :: position, velocity
    CHARACTER(len=*), INTENT(in)         :: frame
    TYPE (Time), INTENT(in)              :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (position(1) < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / new", &
            "Negative distance is not allowed.", 1)
       RETURN
    END IF
    this%position       = position
    this%velocity       = velocity
    CALL checkAngles(this)
    this%frame          = TRIM(frame)
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_SC_posAndVel





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this

    this%position       = 0.0_bp
    this%velocity       = 0.0_bp
    this%frame          = ""
    CALL NULLIFY(this%t)
    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_SC





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this
    TYPE (SphericalCoordinates)             :: copy_SC

    copy_SC%position       = this%position
    copy_SC%velocity       = this%velocity
    copy_SC%frame          = this%frame
    copy_SC%t              = copy(this%t)
    IF (error) THEN
       CALL errorMessage("SphericalCoordinates / copy", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    copy_SC%is_initialized = this%is_initialized

  END FUNCTION copy_SC





  !! *Description*:
  !!
  !! Returns the status of this object, i.e. whether
  !! it exists or not.
  !!
  LOGICAL FUNCTION exist_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this

    exist_SC = this%is_initialized

  END FUNCTION exist_SC





  !! *Description*:
  !!
  !! Adds Gaussian deviates to the coordinates. The new center
  !! relative to the original coordinates and the standard deviations
  !! should be given in radians for the angular space coordinates, AUs
  !! for distance, radians per day for angular velocities and AU per
  !! day for line-of-sight velocity. Outputs the final values for
  !! center (which should not change) and deviates.
  !!
  !! Returns error.
  !!
  SUBROUTINE addMultinormalDeviate_SC(this, mean, covariance)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6), INTENT(in)         :: mean
    REAL(bp), DIMENSION(6,6), INTENT(in)       :: covariance

    REAL(bp), DIMENSION(6,6)                   :: A, covariance_
    REAL(bp), DIMENSION(6)                     :: norm, p, deviates
    REAL(bp)                                   :: cosdelta
    INTEGER                                    :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / addMultinormalDeviate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (ALL(mean == 0.0_bp) .AND. ALL(covariance == 0)) THEN
       RETURN
    END IF

    ! Definitions:
    !
    ! x      :: multivariate deviate vector
    ! v      :: normally distributed (0,1) random variables
    ! SIGMA  :: covariance matrix

    ! sigma RA is larger for a higher declination:
    covariance_ = covariance
    cosdelta = COS(this%position(3))
    IF (ABS(cosdelta) < 10.0_bp*EPSILON(cosdelta)) THEN
       cosdelta = 10.0_bp*EPSILON(cosdelta)
    END IF
    covariance_(2,:) = covariance_(2,:)/cosdelta
    covariance_(:,2) = covariance_(:,2)/cosdelta
    DO i=1,6
       IF (covariance_(i,i) == 0.0_bp) THEN
          covariance_(i,i) = 1.0_bp
       END IF
    END DO

    A = 0.0_bp
    DO i=1,6
       A(i,i:6) = covariance_(i,i:6)
    END DO
    CALL cholesky_decomposition(A, p, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / addMultinormalDeviate:", &
            "Cholesky decomposition unsuccessful:", 1)
       WRITE(stderr,"(A)") TRIM(errstr)
       RETURN
    END IF
    DO i=1,6
       A(i,i) = p(i)
    END DO

    ! dx = A v + mean:
    CALL randomGaussian(norm)
    deviates = MATMUL(A,norm) + mean

    ! x_new = x_old + dx:
    this%position = this%position + deviates(1:3)
    this%velocity = this%velocity + deviates(4:6)

    CALL checkAngles(this)

  END SUBROUTINE addMultinormalDeviate_SC





  !! *Description*:
  !!
  !! Adds uniform deviates to the coordinates. The center relative to
  !! the original coordinates and the absolute values of the boundary
  !! values ((i,1)=center, (i,2)=abs(boundary), position i=1:3 and
  !! velocity i=4:6) should be given in radians for the angular space
  !! coordinates, AUs for distance, radians per day for angular 
  !! velocities and AUs per day for line-of-sight velocity.
  !!
  !! Returns error.
  !!
  SUBROUTINE addUniformDeviate_SC(this, center_and_absbound)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6,2), INTENT(in)       :: center_and_absbound
    REAL(bp), DIMENSION(6,2) :: center_and_absbound_
    REAL(bp), DIMENSION(6) :: ran
    REAL(bp) :: tmp

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / addUniformDeviate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    center_and_absbound_ = center_and_absbound

    ! Error value in RA is larger for a higher declination:
    tmp = COS(this%position(3))
    IF (ABS(tmp) < 10.0_bp*EPSILON(tmp)) THEN
       tmp = 10.0_bp*EPSILON(tmp)
    END IF
    center_and_absbound_(2,2) = center_and_absbound_(2,2)/tmp

    ! Multiply the random numbers with the given limits to get the
    ! final deviates:
    CALL randomNumber(ran)
    center_and_absbound_(1:6,2) = &
         (2.0_bp*ran - 1.0_bp)*center_and_absbound_(1:6,2)

    ! New coordinates = old coordinates + center relative to old
    ! coordinates + deviates:
    this%position = this%position + center_and_absbound_(1:3,1) + &
         center_and_absbound_(1:3,2)
    this%velocity = this%velocity + center_and_absbound_(4:6,1) + &
         center_and_absbound_(4:6,2)

    CALL checkAngles(this)

  END SUBROUTINE addUniformDeviate_SC





  ! Description:
  !
  ! Checks that longitude and latitude are within [0,2pi] and 
  ! [-pi/2,pi/2], respectively.
  !
  SUBROUTINE checkAngles(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this

    ! 0 <= angle < 2pi :
    this%position(2)    = MODULO(this%position(2), two_pi)

    ! -pi/2 <= angle <= pi/2 :
    this%position(3)    = ASIN(SIN(this%position(3)))

  END SUBROUTINE checkAngles





  LOGICAL FUNCTION equal_SC(this, that)

    TYPE (SphericalCoordinates), INTENT(in) :: this
    TYPE (SphericalCoordinates), INTENT(in) :: that

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / equal", &
            "1st object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. that%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / equal", &
            "2nd object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Designation:
    IF (.NOT.(this%frame == that%frame)) THEN
       equal_SC = .FALSE.
       RETURN
    END IF

    ! Position:
    IF (ANY(ABS(this%position - that%position) > EPSILON(this%position(1)))) THEN
       equal_SC = .FALSE.
       RETURN
    END IF

    ! Velocity:
    IF (ANY(ABS(this%velocity - that%velocity) > EPSILON(this%velocity(1)))) THEN
       equal_SC = .FALSE.
       RETURN
    END IF

    ! Epoch
    IF (.NOT.equal(this%t,that%t)) THEN
       equal_SC = .FALSE.
       RETURN
    END IF

    ! Assuming that all of the above comparisons are true,
    ! it can be concluded that the two objects are the same:
    equal_SC = .TRUE.    

  END FUNCTION equal_SC





  !! *Description*:
  !!
  !! Returns the position vector as [AU,rad,rad,AU/day,rad/day,rad/day].
  !!
  !! Returns error if not initialized.
  !!
  FUNCTION getCoordinates_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(6)                  :: getCoordinates_SC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getCoordinates", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getCoordinates_SC(1:3) = this%position
    getCoordinates_SC(4:6) = this%velocity

  END FUNCTION getCoordinates_SC





  !! *Description*:
  !!
  !! Returns the position vector as [AU,rad,rad].
  !!
  !! Returns error if not initialized.
  !!
  FUNCTION getPosition_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(3)            :: getPosition_SC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getPosition", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getPosition_SC = this%position

  END FUNCTION getPosition_SC





  !! *Description*:
  !!
  !! Returns the velocity vector as [AU/day,rad/day,rad/day].
  !!
  !! Returns error.
  !!
  FUNCTION getVelocity_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(3)            :: getVelocity_SC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getVelocity", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getVelocity_SC = this%velocity

  END FUNCTION getVelocity_SC





  !! *Description*:
  !!
  !! Returns coordinate frame.
  !!
  !! Returns error.
  !!
  CHARACTER(len=FRAME_LEN) FUNCTION getFrame_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getFrame", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getFrame_SC = TRIM(this%frame)

  END FUNCTION getFrame_SC





  !! *Description*:
  !!
  !! Returns the epoch as a Time object.
  !!
  !! Returns error.
  !!
  FUNCTION getTime_SC(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this
    TYPE (Time)                       :: getTime_SC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getTime_SC = copy(this%t)

  END FUNCTION getTime_SC





  !! *Description*:
  !!
  !! Returns distance.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getDistance(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getDistance = this%position(1)

  END FUNCTION getDistance





  !! *Description*:
  !!
  !! Returns longitude (or Right Ascension).
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getLongitude(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getLongitude", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getLongitude = this%position(2)

  END FUNCTION getLongitude





  !! *Description*:
  !!
  !! Returns latitude (or Declination).
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getLatitude(this)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / getLatitude", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getLatitude = this%position(3)

  END FUNCTION getLatitude





  !! *Description*:
  !!
  !! Reallocates a pointer array of SphericalCoordinates-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myscoords => reallocate(myscoords,4)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_SC_1(array,n)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: reallocate_SC_1, array
    INTEGER, INTENT(in)                   :: n
    INTEGER                               :: i, nold, err

    ALLOCATE(reallocate_SC_1(n), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_SC_1 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,dim=1)
    DO i=1, MIN(n,nold)
       reallocate_SC_1(i) = copy(array(i))
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_SC_1





  !! *Description*:
  !!
  !! Reallocates a pointer array of SphericalCoordinates-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myscoords => reallocate(myscoords,4,2)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_SC_2(array,n,m)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: reallocate_SC_2, array
    INTEGER, INTENT(in)                   :: n, m
    INTEGER                               :: i, j, nold, mold, err

    ALLOCATE(reallocate_SC_2(n,m), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_SC_2 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,dim=1)
    mold = SIZE(array,dim=2)
    DO i=1, MIN(n,nold)
       DO j=1, MIN(m,mold)
          reallocate_SC_2(i,j) = copy(array(i,j))
       END DO
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_SC_2





  !! *Description*:
  !!
  !!
  !! Rotates position coordinates to equatorial coordinates.
  !! Optionally, returns partial derivatives between final and initial
  !! coordinates. An identity matrix is returned if the coordinates
  !! are not changed, otherwise:
  !! 
  !! dr_eq/dr_ec    dr_eq/dlon_ec    dr_eq/dlat_ec    dr_eq/ddr_ec    dr_eq/ddlon_ec    dr_eq/ddlat_ec  
  !!
  !! dlon_eq/dr_ec  dlon_eq/dlon_ec  dlon_eq/dlat_ec  dlon_eq/ddr_ec  dlon_eq/ddlon_ec  dlon_eq/ddlat_ec  
  !!
  !! dlat_eq/dr_ec  dlat_eq/dlon_ec  dlat_eq/dlat_ec  dlat_eq/ddr_ec  dlat_eq/ddlon_ec  dlat_eq/ddlat_ec  
  !!
  !! ddr_eq/dr_ec   ddr_eq/dlon_ec   ddr_eq/dlat_ec   ddr_eq/ddr_ec   ddr_eq/ddlon_ec   ddr_eq/ddlat_ec  
  !!
  !! ddlon_eq/dr_ec ddlon_eq/dlon_ec ddlon_eq/dlat_ec ddlon_eq/ddr_ec ddlon_eq/ddlon_ec ddlon_eq/ddlat_ec  
  !!
  !! ddlat_eq/dr_ec ddlat_eq/dlon_ec ddlat_eq/dlat_ec ddlat_eq/ddr_ec ddlat_eq/ddlon_ec ddlat_eq/ddlat_ec  
  !!
  SUBROUTINE rotateToEquatorial_SC(this, partials)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: partials
    REAL(bp) :: lon_ec, lat_ec, sin_lon_ec, cos_lon_ec, sin_lat_ec, &
         cos_lat_ec, sin_lat_eq, cos_lat_eq, sin_lon_eq, cos_lon_eq, &
         dlon_ec, dlat_ec, sin_eps, cos_eps

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / rotateToEquatorial", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%frame /= "equatorial") THEN

       ! Position
       lon_ec = this%position(2)
       lat_ec  = this%position(3)
       sin_lon_ec = SIN(lon_ec)
       cos_lon_ec = COS(lon_ec)
       sin_lat_ec = SIN(lat_ec)
       cos_lat_ec = COS(lat_ec)
       sin_eps = SIN(eps)
       cos_eps = COS(eps)
       sin_lat_eq = sin_lat_ec*cos_eps + &
            cos_lat_ec*sin_eps*sin_lon_ec
       IF (ABS(sin_lat_eq) > 1.0_bp) THEN
          sin_lat_eq = SIGN(1.0_bp, sin_lat_eq)
       END IF
       this%position(3) = ASIN(sin_lat_eq)
       cos_lat_eq = COS(this%position(3))
       sin_lon_eq = (-sin_lat_ec*sin_eps + &
            cos_lat_ec*cos_eps*sin_lon_ec) / &
            cos_lat_eq
       cos_lon_eq = cos_lat_ec*cos_lon_ec / &
            cos_lat_eq
       this%position(2) = angle(cos_lon_eq, sin_lon_eq)

       ! Velocity
       dlon_ec = this%velocity(2)
       dlat_ec  = this%velocity(3)
       this%velocity(3) = (dlat_ec*cos_lat_ec*cos_eps - &
            dlat_ec*sin_lat_ec*sin_eps*sin_lon_ec + &
            dlon_ec*cos_lat_ec*sin_eps*cos_lon_ec) / &
            cos_lat_eq
       this%velocity(2) = (dlat_ec*sin_lat_ec*cos_lon_ec + &
            dlon_ec*cos_lat_ec*sin_lon_ec - &
            this%velocity(3)*cos_lon_eq*sin_lat_eq) / &
            (sin_lon_eq*cos_lat_eq)

       CALL checkAngles(this)
       this%frame = "equatorial"

       IF (PRESENT(partials)) THEN
          partials(1,:) = (/ 1.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp /)
          partials(2,:) = (/ 0.0_bp, &
               cos_lat_ec*sin_lon_ec/(SQRT(1.0_bp-cos_lon_eq**2)*cos_lat_eq), &
               sin_lat_ec*cos_lon_ec/(SQRT(1.0_bp-cos_lon_eq**2)*cos_lat_eq), &
               0.0_bp, 0.0_bp, 0.0_bp /)
          partials(3,:) = (/ 0.0_bp, &
               cos_lat_ec*sin_eps*cos_lon_ec/SQRT(1.0_bp-sin_lat_eq**2), &
               (cos_lat_ec*cos_eps-sin_lat_ec*sin_eps*sin_lon_ec)/SQRT(1.0_bp-sin_lat_eq**2), &
               0.0_bp, 0.0_bp, 0.0_bp /)
          partials(4,:) = (/ 0.0_bp, 0.0_bp, 0.0_bp, 1.0_bp, 0.0_bp, 0.0_bp /)
          partials(5,:) = (/ 0.0_bp, &
               (-dlat_ec*sin_lat_ec*sin_lon_ec+dlon_ec*cos_lat_ec*cos_lon_ec)/(sin_lon_eq*cos_lat_eq), &
               (dlat_ec*cos_lat_ec*cos_lon_ec-dlon_ec*sin_lat_ec*sin_lon_ec)/(sin_lon_eq*cos_lat_eq), &
               0.0_bp, &
               (cos_lat_ec*sin_lon_ec)/(sin_lon_eq*cos_lat_eq), &
               (sin_lat_ec*cos_lon_ec)/(sin_lon_eq*cos_lat_eq) /)
          partials(6,:) = (/ 0.0_bp, &
               (-dlat_ec*sin_lat_ec*sin_eps*cos_lon_ec-dlon_ec*cos_lat_ec*sin_eps*sin_lon_ec)/cos_lat_eq, &
               (-dlat_ec*sin_lat_ec*cos_eps - &
               dlat_ec*cos_lat_ec*sin_eps*sin_lon_ec - &
               dlon_ec*sin_lat_ec*sin_eps*cos_lon_ec)/cos_lat_eq, &
               0.0_bp, &
               (cos_lat_ec*sin_eps*cos_lon_ec)/cos_lat_eq, &
               (cos_lat_ec*cos_eps-sin_lat_ec*sin_eps*sin_lon_ec)/cos_lat_eq /)
       END IF

    ELSE

       IF (PRESENT(partials)) THEN
          partials = identity_matrix(6)
       END IF

    END IF

  END SUBROUTINE rotateToEquatorial_SC





  !! *Description*:
  !!
  !! Rotates position coordinates to ecliptic coordinates.
  !! Optionally, returns partial derivatives between final and initial
  !! coordinates. An identity matrix is returned if the coordinates
  !! are not changed, otherwise:
  !! 
  !! dr_ec/dr_eq    dr_ec/dlon_eq    dr_ec/dlat_eq    dr_ec/ddr_eq    dr_ec/ddlon_eq    dr_ec/ddlat_eq  
  !!
  !! dlon_ec/dr_eq  dlon_ec/dlon_eq  dlon_ec/dlat_eq  dlon_ec/ddr_eq  dlon_ec/ddlon_eq  dlon_ec/ddlat_eq  
  !!
  !! dlat_ec/dr_eq  dlat_ec/dlon_eq  dlat_ec/dlat_eq  dlat_ec/ddr_eq  dlat_ec/ddlon_eq  dlat_ec/ddlat_eq  
  !!
  !! ddr_ec/dr_eq   ddr_ec/dlon_eq   ddr_ec/dlat_eq   ddr_ec/ddr_eq   ddr_ec/ddlon_eq   ddr_ec/ddlat_eq  
  !!
  !! ddlon_ec/dr_eq ddlon_ec/dlon_eq ddlon_ec/dlat_eq ddlon_ec/ddr_eq ddlon_ec/ddlon_eq ddlon_ec/ddlat_eq  
  !!
  !! ddlat_ec/dr_eq ddlat_ec/dlon_eq ddlat_ec/dlat_eq ddlat_ec/ddr_eq ddlat_ec/ddlon_eq ddlat_ec/ddlat_eq  
  !!
  SUBROUTINE rotateToEcliptic_SC(this, partials)

    IMPLICIT NONE
    TYPE (SphericalCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: partials
    REAL(bp) :: lon_eq, lat_eq, sin_lon_eq, cos_lon_eq, sin_lat_eq, &
         cos_lat_eq, sin_lat_ec, cos_lat_ec, sin_lon_ec, cos_lon_ec, &
         dlon_eq, dlat_eq, sin_eps, cos_eps

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("SphericalCoordinates / rotateToEcliptic", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%frame /= "ecliptic") THEN

       ! Position
       lon_eq = this%position(2)
       lat_eq  = this%position(3)
       sin_lon_eq = SIN(lon_eq)
       cos_lon_eq = COS(lon_eq)
       sin_lat_eq = SIN(lat_eq)
       cos_lat_eq = COS(lat_eq)
       sin_eps = SIN(eps)
       cos_eps = COS(eps)
       sin_lat_ec = sin_lat_eq*cos_eps - &
            cos_lat_eq*sin_eps*sin_lon_eq
       IF (ABS(sin_lat_ec) > 1.0_bp) THEN
          sin_lat_ec = SIGN(1.0_bp, sin_lat_ec)
       END IF
       this%position(3) = ASIN(sin_lat_ec)
       cos_lat_ec = COS(this%position(3))
       sin_lon_ec = (sin_lat_eq*sin_eps + &
            cos_lat_eq*cos_eps*sin_lon_eq) / &
            cos_lat_ec
       cos_lon_ec = cos_lat_eq*cos_lon_eq / &
            cos_lat_ec
       this%position(2) = angle(cos_lon_ec, sin_lon_ec)

       ! Velocity
       dlon_eq = this%velocity(2)
       dlat_eq  = this%velocity(3)
       this%velocity(3) = (dlat_eq*cos_lat_eq*cos_eps + &
            dlat_eq*sin_lat_eq*sin_eps*sin_lon_eq - &
            dlon_eq*cos_lat_eq*sin_eps*cos_lon_eq) / &
            cos_lat_ec
       this%velocity(2) = (dlat_eq*sin_lat_eq*cos_lon_eq + &
            dlon_eq*cos_lat_eq*sin_lon_eq - &
            this%velocity(3)*cos_lon_ec*sin_lat_ec) / &
            (sin_lon_ec*cos_lat_ec)

       CALL checkAngles(this)
       this%frame = "ecliptic"

       IF (PRESENT(partials)) THEN
          partials(1,:) = (/ 1.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp /)
          partials(2,:) = (/ 0.0_bp, &
               cos_lat_eq*sin_lon_eq/(SQRT(1.0_bp-cos_lon_ec**2)*cos_lat_ec), &
               sin_lat_eq*cos_lon_eq/(SQRT(1.0_bp-cos_lon_ec**2)*cos_lat_ec), &
               0.0_bp, 0.0_bp, 0.0_bp /)
          partials(3,:) = (/ 0.0_bp, &
               -cos_lat_eq*sin_eps*cos_lon_eq/SQRT(1.0_bp-sin_lat_ec**2), &
               (cos_lat_eq*cos_eps+sin_lat_eq*sin_eps*sin_lon_eq)/SQRT(1.0_bp-sin_lat_ec**2), &
               0.0_bp, 0.0_bp, 0.0_bp /)
          partials(4,:) = (/ 0.0_bp, 0.0_bp, 0.0_bp, 1.0_bp, 0.0_bp, 0.0_bp /)
          partials(5,:) = (/ 0.0_bp, &
               (-dlat_eq*sin_lat_eq*sin_lon_eq+dlon_eq*cos_lat_eq*cos_lon_eq)/(sin_lon_ec*cos_lat_ec), &
               (dlat_eq*cos_lat_eq*cos_lon_eq-dlon_eq*sin_lat_eq*sin_lon_eq)/(sin_lon_ec*cos_lat_ec), &
               0.0_bp, &
               (cos_lat_eq*sin_lon_eq)/(sin_lon_ec*cos_lat_ec), &
               (sin_lat_eq*cos_lon_eq)/(sin_lon_ec*cos_lat_ec) /)
          partials(6,:) = (/ 0.0_bp, &
               (dlat_eq*sin_lat_eq*sin_eps*cos_lon_eq+dlon_eq*cos_lat_eq*sin_eps*sin_lon_eq)/cos_lat_ec, &
               (-dlat_eq*sin_lat_eq*cos_eps + &
               dlat_eq*cos_lat_eq*sin_eps*sin_lon_eq + &
               dlon_eq*sin_lat_eq*sin_eps*cos_lon_eq)/cos_lat_ec, &
               0.0_bp, &
               (-cos_lat_eq*sin_eps*cos_lon_eq)/cos_lat_ec, &
               (cos_lat_eq*cos_eps+sin_lat_eq*sin_eps*sin_lon_eq)/cos_lat_ec /)
       END IF

    ELSE

       IF (PRESENT(partials)) THEN
          partials = identity_matrix(6)
       END IF

    END IF

  END SUBROUTINE rotateToEcliptic_SC




END MODULE SphericalCoordinates_cl
