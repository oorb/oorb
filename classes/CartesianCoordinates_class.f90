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
!! *Class*description*: 
!!  
!! Type and routines for static Cartesian coordinates.
!!
!! @see Orbit_class 
!!  
!! @author  MG, JV, KM, TL
!! @version 2009-12-31
!!
MODULE CartesianCoordinates_cl

  USE Base_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  IMPLICIT NONE

  PRIVATE :: new_CC
  PRIVATE :: new_CC_coordinates
  PRIVATE :: new_CC_spherical
  PRIVATE :: new_CC_values
  PRIVATE :: nullify_CC
  PRIVATE :: copy_CC
  PRIVATE :: exist_CC
  PRIVATE :: getCoordinates_CC
  PRIVATE :: getSCoord_CC
  PRIVATE :: getPosition_CC
  PRIVATE :: getVelocity_CC
  PRIVATE :: getFrame_CC
  PRIVATE :: getTime_CC
  PRIVATE :: opposite_CC
  PRIVATE :: rotateToEcliptic_CC
  PRIVATE :: rotateToEcliptic_simple
  PRIVATE :: rotateToEquatorial_CC
  PRIVATE :: rotateToEquatorial_simple
  PRIVATE :: addition_CC
  PRIVATE :: subtraction_CC
  PRIVATE :: partialsScoordWrtCcoord_CC

  TYPE CartesianCoordinates
     PRIVATE
     REAL(bp), DIMENSION(3)   :: position
     REAL(bp), DIMENSION(3)   :: velocity
     ! Frame can be ecliptic or equatorial.
     CHARACTER(len=FRAME_LEN) :: frame = "ecliptic"
     TYPE (Time)              :: t
     LOGICAL                  :: is_initialized = .FALSE.
  END TYPE CartesianCoordinates

  INTERFACE NEW
     MODULE PROCEDURE new_CC
     MODULE PROCEDURE new_CC_coordinates
     MODULE PROCEDURE new_CC_spherical
     MODULE PROCEDURE new_CC_values
  END INTERFACE

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_CC
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_CC
  END INTERFACE

  INTERFACE exist
     MODULE PROCEDURE exist_CC
  END INTERFACE

  INTERFACE getCoordinates
     MODULE PROCEDURE getCoordinates_CC
  END INTERFACE

  INTERFACE getPosition
     MODULE PROCEDURE getPosition_CC
  END INTERFACE

  INTERFACE getVelocity
     MODULE PROCEDURE getVelocity_CC
  END INTERFACE

  INTERFACE getFrame
     MODULE PROCEDURE getFrame_CC
  END INTERFACE

  INTERFACE getSCoord
     MODULE PROCEDURE getSCoord_CC
  END INTERFACE

  INTERFACE getTime
     MODULE PROCEDURE getTime_CC
  END INTERFACE

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_CC
     MODULE PROCEDURE reallocate_CC_1
     MODULE PROCEDURE reallocate_CC_2
  END INTERFACE

  INTERFACE rotateToEquatorial
     MODULE PROCEDURE rotateToEquatorial_CC
     MODULE PROCEDURE rotateToEquatorial_simple
  END INTERFACE

  INTERFACE rotateToEcliptic
     MODULE PROCEDURE rotateToEcliptic_CC
     MODULE PROCEDURE rotateToEcliptic_simple
  END INTERFACE

  !! Returns a Cartesian state, which is the sum of the two given
  !! Cartesian states (ecliptic or equatorial).
  INTERFACE OPERATOR (+) 
     MODULE PROCEDURE addition_CC
  END INTERFACE

  !! Returns a Cartesian state, which is the difference of the two
  !! given Cartesian states (ecliptic or equatorial).
  INTERFACE OPERATOR (-) 
     MODULE PROCEDURE subtraction_CC
  END INTERFACE

  INTERFACE opposite
     MODULE PROCEDURE opposite_CC
  END INTERFACE

  INTERFACE partialsScoordWrtCcoord
     MODULE PROCEDURE partialsScoordWrtCcoord_CC
  END INTERFACE

CONTAINS





  !! *Description*:
  !!
  !! Initializes a CartesianCoordinates-object. Position and velocity
  !! are zero-vectors. Coordinate frame is ecliptic.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%position       = (/ 0.0_bp, 0.0_bp, 0.0_bp /)
    this%velocity       = (/ 0.0_bp, 0.0_bp, 0.0_bp /)
    this%frame          = "ecliptic"
    CALL NULLIFY(this%t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_CC





  !! *Description*:
  !!
  !! Initializes a CartesianCoordinates-object using a given
  !! coordinates, coordinate frame, and epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_CC_coordinates(this, coordinates, frame, t)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6), INTENT(in)         :: coordinates
    CHARACTER(len=*), INTENT(in)               :: frame
    TYPE (Time), INTENT(in)                    :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%position = coordinates(1:3)
    this%velocity = coordinates(4:6)

    IF (TRIM(frame) /= "equatorial" .AND. &
         TRIM(frame) /= "ecliptic") THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / new", &
            "Wrong frame of reference. " // &
            "Choose either 'equatorial' or 'ecliptic'.", 1)
       RETURN
    ELSE
       this%frame       = TRIM(frame)
    END IF
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_CC_coordinates





  !! *Description*:
  !!
  !! Initializes a CartesianCoordinates-object using a given
  !! SphericalCoordinates-object.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_CC_spherical(this, scoord)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this
    TYPE (SphericalCoordinates), INTENT(in)    :: scoord
    REAL(bp), DIMENSION(3)                     :: position, velocity

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    position = getPosition(scoord)
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / new",&
            "TRACE BACK 1", 1)
       RETURN
    END IF
    this%position(1) = position(1) * COS(position(2)) * COS(position(3))
    this%position(2) = position(1) * SIN(position(2)) * COS(position(3))
    this%position(3) = position(1) * SIN(position(3))

    velocity = getVelocity(scoord)
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / new",&
            "TRACE BACK 2", 1)
       RETURN
    END IF
    this%velocity(1) = COS(position(2))*COS(position(3))*velocity(1) - &
         position(1)*SIN(position(2))*COS(position(3))*velocity(2)   - &
         position(1)*COS(position(2))*SIN(position(3))*velocity(3)
    this%velocity(2) = SIN(position(2))*COS(position(3))*velocity(1) + &
         position(1)*COS(position(2))*COS(position(3))*velocity(2)   - &
         position(1)*SIN(position(2))*SIN(position(3))*velocity(3)
    this%velocity(3) = SIN(position(3))*velocity(1) + &
         position(1)*COS(position(3))*velocity(3)

    this%frame = getFrame(scoord)
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / new",&
            "TRACE BACK 3", 1)
       RETURN
    END IF

    this%t = getTime(scoord)
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / new",&
            "TRACE BACK 4", 1)
       RETURN
    END IF
    this%is_initialized = .TRUE.

  END SUBROUTINE new_CC_spherical





  !! *Description*:
  !!
  !! Initializes a CartesianCoordinates-object using a given position,
  !! velocity, coordinate frame, and epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_CC_values(this, position, velocity, frame, t)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(3), INTENT(in)         :: position
    REAL(bp), DIMENSION(3), INTENT(in)         :: velocity
    CHARACTER(len=*), INTENT(in)               :: frame
    TYPE (Time), INTENT(in)                    :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%position       = position
    this%velocity       = velocity

    IF (TRIM(frame) /= "equatorial" .AND. &
         TRIM(frame) /= "ecliptic") THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / new", &
            "Wrong frame of reference. " // &
            "Choose either 'equatorial' or 'ecliptic'.", 1)
       RETURN
    ELSE
       this%frame       = TRIM(frame)
    END IF
    this%t              = copy(t)
    this%is_initialized = .TRUE.

  END SUBROUTINE new_CC_values





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this

    this%position       = 0.0_bp
    this%velocity       = 0.0_bp
    this%frame          = " "
    CALL NULLIFY(this%t)
    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_CC





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    TYPE (CartesianCoordinates)             :: copy_CC

    copy_CC%position       = this%position
    copy_CC%velocity       = this%velocity
    copy_CC%frame          = TRIM(this%frame)
    copy_CC%t              = copy(this%t)
    copy_CC%is_initialized = this%is_initialized

  END FUNCTION copy_CC





  !! *Description*:
  !!
  !! Returns the status of the object, i.e. whether it exists or not.
  !!
  LOGICAL FUNCTION exist_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this

    exist_CC = this%is_initialized

  END FUNCTION exist_CC





  !! *Description*:
  !!
  !! Adds "that" to "this", and returns the sum as a new object. 
  !! Makes sure they are given in the same frame of reference
  !! and that they are in the same epoch before adding. Returns
  !! a nullified object if an error occurs.
  !!
  !! *NOTE*:
  !! *Do*not*use*this*function*by*name,*but*by*using*the*addition*operator*!
  !!
  FUNCTION addition_CC(this, that)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    TYPE (CartesianCoordinates), INTENT(in) :: that
    TYPE (CartesianCoordinates)             :: addition_CC, this_, that_

    IF (.NOT. this%is_initialized) THEN
       CALL errorMessage("CartesianCoordinates / + (addition)", &
            "Left object has not been initialized.", 1)
       CALL NULLIFY(addition_CC)
       RETURN       
    END IF

    IF (.NOT. that%is_initialized) THEN
       CALL errorMessage("CartesianCoordinates / + (addition)", &
            "Right object has not been initialized.", 1)
       CALL NULLIFY(addition_CC)
       RETURN       
    END IF

    IF (this%frame /= that%frame) THEN
       CALL errorMessage("CartesianCoordinates / + (addition)", &
            "Coordinates are not given in the same frame of reference.", 1)
       CALL NULLIFY(addition_CC)
       RETURN       
    END IF

    this_ = copy(this)
    that_ = copy(that)

    IF (equal(this_%t,that_%t)) THEN
       addition_CC = copy(this_)
       addition_CC%position = this_%position + that_%position
       addition_CC%velocity = this_%velocity + that_%velocity
    ELSE
       CALL errorMessage("CartesianCoordinates / + (addition)", &
            "Coordinatess at different epochs cannot be used in addition.", 1)
       CALL NULLIFY(addition_CC)
       RETURN
    END IF

  END FUNCTION addition_CC





  !! *Description*:
  !!
  !! Computes light-time-corrected epoch for object's
  !! Cartesian coordinates given the distance to the observer.
  !!
  !! Returns error.
  !!
  SUBROUTINE estimateLightTime(this, distance)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this
    REAL(bp), INTENT(in)                       :: distance
    REAL(bp)                                   :: tdt

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / estimateLightTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    tdt = getMJD(this%t, "tdt")
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / estimateLightTime", &
            "TRACE BACK 1", 1)
       RETURN
    END IF
    CALL NULLIFY(this%t)

    ! Subtract the time it takes for light to travel from the target
    ! to the observer:
    CALL NEW(this%t, tdt - distance/sol, "tdt")
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / estimateLightTime", &
            "TRACE BACK 2", 1)
       RETURN
    END IF

  END SUBROUTINE estimateLightTime





  !! *Description*:
  !!
  !! Rotates coordinates from the ecliptical to the equatorial
  !! coordinate frame.
  !!
  SUBROUTINE rotateToEquatorial_simple(w)

    IMPLICIT NONE
    REAL(bp), DIMENSION(6), INTENT(inout) :: w
    REAL(bp), DIMENSION(3,3)              :: R

    R(1,:) = (/ 1.0_bp,   0.0_bp,    0.0_bp /)
    R(2,:) = (/ 0.0_bp, COS(eps), -SIN(eps) /) 
    R(3,:) = (/ 0.0_bp, SIN(eps),  COS(eps) /)
    w(1:3) = MATMUL(R,w(1:3))
    w(4:6) = MATMUL(R,w(4:6))

  END SUBROUTINE rotateToEquatorial_simple





  !! *Description*:
  !!
  !! Rotates coordinates from the equatorial to the ecliptical
  !! coordinate frame.
  !!
  SUBROUTINE rotateToEcliptic_simple(w)

    IMPLICIT NONE
    REAL(bp), DIMENSION(6), INTENT(inout) :: w
    REAL(bp), DIMENSION(3,3)              :: R

    R(1,:) = (/ 1.0_bp,    0.0_bp,   0.0_bp /)
    R(2,:) = (/ 0.0_bp,  COS(eps), SIN(eps) /) 
    R(3,:) = (/ 0.0_bp, -SIN(eps), COS(eps) /)
    w(1:3) = MATMUL(R,w(1:3))
    w(4:6) = MATMUL(R,w(4:6))

  END SUBROUTINE rotateToEcliptic_simple





  !! *Description*:
  !!
  !! Returns the partial derivatives of spherical coordinates 
  !! (ra,dec) wrt Cartesian coordinates.
  !!
  !! @author TL
  !!
  SUBROUTINE partialsSCoordWrtCCoord_CC(this, partials)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(2,6), INTENT(out)   :: partials
    REAL(bp), DIMENSION(3) :: pos
    REAL(bp) :: tmp1, tmp2, tmp3

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / " // &
            "partialsSCoordWrtCCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Common terms.
    pos  = getPosition(this)
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / " // &
            "partialsSCoordWrtCCoord", &
            "TRACE BACK", 1)
       RETURN
    END IF
    tmp1 = 1.0_bp/SUM(pos(1:2)**2)
    tmp2 = SQRT(tmp1)
    tmp3 = 1.0_bp/SUM(pos**2)

    ! First row: ra
    partials(1,1)   = -pos(2)*tmp1
    partials(1,2)   = pos(1)*tmp1
    partials(1,3:6) = 0.0_bp
    ! Second row: dec
    partials(2,1)   = -pos(1)*pos(3)*tmp2*tmp3
    partials(2,2)   = -pos(2)*pos(3)*tmp2*tmp3
    partials(2,3)   = tmp2*(1-pos(3)**2*tmp3)
    partials(2,4:6) = 0.0_bp

  END SUBROUTINE partialsSCoordWrtCCoord_CC





  !! *Description*:
  !!
  !! Returns the coordinate frame.
  !!
  !! Returns error.
  !!
  CHARACTER(len=FRAME_LEN) FUNCTION getFrame_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / getFrame", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getFrame_CC = TRIM(this%frame)

  END FUNCTION getFrame_CC





  !! *Description*:
  !!
  !! Returns the position vector as [AU,AU,AU,AU/day,AU/day,AU/day].
  !!
  !! Returns error if not initialized.
  !!
  FUNCTION getCoordinates_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(6)                  :: getCoordinates_CC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / getCoordinates", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getCoordinates_CC(1:3) = this%position
    getCoordinates_CC(4:6) = this%velocity

  END FUNCTION getCoordinates_CC





  !! *Description*:
  !!
  !! Returns the Cartesian position vector as AUs.
  !!
  !! Returns error.
  !!
  FUNCTION getPosition_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(3)                  :: getPosition_CC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / getPosition", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getPosition_CC = this%position

  END FUNCTION getPosition_CC





  !! *Description*:
  !!
  !! Transforms Cartesian coordinates to spherical coordinates and
  !! returns them as a SphericalCoordinates object with the same frame
  !! as the CartesianCoordinates object.
  !!
  !! Returns error.
  !!
  FUNCTION getSCoord_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    TYPE (SphericalCoordinates)             :: getSCoord_CC
    REAL(bp), DIMENSION(3)                  :: position, velocity
    REAL(bp)                                :: sin_theta, cos_theta, &
         sin_phi, r, rv

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / getSCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ! rho
    position(1) = SQRT(DOT_PRODUCT(this%position,this%position))
    ! theta
    sin_theta = this%position(2) / &
         SQRT(DOT_PRODUCT(this%position(1:2),this%position(1:2)))
    cos_theta = this%position(1) / &
         SQRT(DOT_PRODUCT(this%position(1:2),this%position(1:2)))
    IF (ABS(cos_theta) > 1.0_bp) cos_theta = SIGN(1.0_bp, cos_theta)
    position(2) = ACOS(cos_theta)
    IF (sin_theta < 0.0_bp) position(2) = two_pi - position(2)
    ! astronomical phi (= pi/2 - mathematical phi)
    sin_phi = this%position(3)/position(1)
    IF (ABS(sin_phi) > 1.0_bp) sin_phi = SIGN(1.0_bp, sin_phi)
    position(3) = ASIN(sin_phi)

    r = position(1)
    rv = DOT_PRODUCT(this%position,this%velocity)
    velocity(1) = rv/r
    velocity(2) = (this%position(1)*this%velocity(2) - &
         this%velocity(1)*this%position(2)) / &
         (this%position(1)**2.0_bp + &
         this%position(2)**2.0_bp)
    velocity(3) = (r**2.0_bp*this%velocity(3) - rv*this%position(3)) / &
         (r**2.0_bp * SQRT(this%position(1)**2.0_bp + &
         this%position(2)**2.0_bp))

    CALL NEW(getSCoord_CC, position, velocity, this%frame, copy(this%t))
    IF (error) THEN
       CALL errorMessage("CartesianCoordinates / getSCoord", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

  END FUNCTION getSCoord_CC





  !! *Description*:
  !!
  !! Returns the Time-object corresponding to this state.
  !!
  !! Returns error.
  !!
  FUNCTION getTime_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    TYPE (Time)                             :: getTime_CC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / getTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getTime_CC = copy(this%t)

  END FUNCTION getTime_CC





  !! *Description*:
  !!
  !! Returns the Cartesian velocity vector as AUs per day.
  !!
  !! Returns error.
  !!
  FUNCTION getVelocity_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    REAL(bp), DIMENSION(3)                  :: getVelocity_CC

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / getVelocity", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getVelocity_CC = this%velocity

  END FUNCTION getVelocity_CC






  !! *Description*:
  !!
  !! Computes the opposite Cartesian coordinates, e.g., position and
  !! velocity multiplied by -1.
  !!
  !! Returns error. 
  !!
  FUNCTION opposite_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    TYPE (CartesianCoordinates)             :: opposite_CC

    opposite_CC = copy(this)
    opposite_CC%position = -1.0_bp * this%position
    opposite_CC%velocity = -1.0_bp * this%velocity

  END FUNCTION opposite_CC





  !! *Description*:
  !!
  !! Reallocates a pointer array of CartesianCoordinates-objects by
  !! copying only the initialized objects to the new array.
  !!
  !! *Usage*:
  !!
  !! myccoords => reallocate(myccoords)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_CC(array)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: reallocate_CC, array
    INTEGER :: i, j, n, err

    n = 0
    DO i=1,SIZE(array)
       IF (exist(array(i))) THEN
          n = n + 1
       END IF
    END DO
    IF (n == 0) THEN
       reallocate_CC => NULL()
    ELSE
       ALLOCATE(reallocate_CC(n), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("CartesianCoordinates / reallocate", &
               "Could not allocate memory.", 1)
          DEALLOCATE(array)
          reallocate_CC => NULL()
          RETURN
       END IF
       j = 0
       DO i=1,SIZE(array)
          IF (exist(array(i))) THEN
             j = j + 1
             reallocate_CC(j) = copy(array(i))
             IF (error) THEN
                DEALLOCATE(array)
                CALL errorMessage("CartesianCoordinates / reallocate", &
                     "TRACE BACK", 1)
                RETURN          
             END IF
          END IF
       END DO
    END IF
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_CC





  !! *Description*:
  !!
  !! Reallocates a pointer array of CartesianCoordinates-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myccoords => reallocate(myccoords,4)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_CC_1(array, n)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: reallocate_CC_1, array
    INTEGER, INTENT(in)                   :: n
    INTEGER                               :: i, nold, err

    ALLOCATE(reallocate_CC_1(n), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_CC_1 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,dim=1)
    DO i=1, MIN(n,nold)
       reallocate_CC_1(i) = copy(array(i))
       IF (error) THEN
          CALL errorMessage("CartesianCoordinates / reallocate", &
               "TRACE BACK", 1)
          RETURN          
       END IF
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_CC_1





  !! *Description*:
  !!
  !! Reallocates a pointer array of CartesianCoordinates-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myccoords => reallocate(myccoords,4,2)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_CC_2(array, n, m)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), DIMENSION(:,:), POINTER :: reallocate_CC_2, array
    INTEGER, INTENT(in)                   :: n, m
    INTEGER                               :: i, j, nold, mold, err

    ALLOCATE(reallocate_CC_2(n,m), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_CC_2 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,dim=1)
    mold = SIZE(array,dim=2)
    DO i=1, MIN(n,nold)
       DO j=1, MIN(m,mold)
          reallocate_CC_2(i,j) = copy(array(i,j))
       END DO
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_CC_2





  !! *Description*:
  !!
  !! Rotates coordinates to the ecliptic coordinate frame, if not
  !! already the case.
  !!
  SUBROUTINE rotateToEcliptic_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6)                     :: coord

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / rotateToEcliptic", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%frame /= "ecliptic") THEN
       coord(1:3) = this%position
       coord(4:6) = this%velocity
       CALL rotateToEcliptic(coord)
       this%position = coord(1:3)
       this%velocity = coord(4:6)
       this%frame = "ecliptic"
    END IF

  END SUBROUTINE rotateToEcliptic_CC





  !! *Description*:
  !!
  !! Rotates coordinates to the equatorial coordinate frame, if not
  !! already the case.
  !!
  SUBROUTINE rotateToEquatorial_CC(this)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(inout) :: this
    REAL(bp), DIMENSION(6)                     :: coord

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("CartesianCoordinates / rotateToEquatorial", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%frame /= "equatorial") THEN
       coord(1:3) = this%position
       coord(4:6) = this%velocity
       CALL rotateToEquatorial(coord)
       this%position = coord(1:3)
       this%velocity = coord(4:6)
       this%frame = "equatorial"
    END IF

  END SUBROUTINE rotateToEquatorial_CC





  !! *Description*:
  !!
  !! Subtracts "that" from "this", and returns the difference as a new
  !! object. Makes sure they are given in the same frame of reference
  !! and that they are in the same epoch before subtracting. Returns
  !! a nullified object if an error occurs.
  !!
  !! *NOTE*:
  !! *Do*not*use*this*function*by*name,*but*by*using*the*subtraction*operator*!
  !!
  FUNCTION subtraction_CC(this, that)

    IMPLICIT NONE
    TYPE (CartesianCoordinates), INTENT(in) :: this
    TYPE (CartesianCoordinates), INTENT(in) :: that
    TYPE (CartesianCoordinates)             :: subtraction_CC, this_, that_

    IF (.NOT. this%is_initialized) THEN
       CALL errorMessage("CartesianCoordinates / - (subtraction)", &
            "Left object has not been initialized.", 1)
       CALL NULLIFY(subtraction_CC)
       RETURN       
    END IF

    IF (.NOT. that%is_initialized) THEN
       CALL errorMessage("CartesianCoordinates / - (subtraction)", &
            "Right object has not been initialized.", 1)
       CALL NULLIFY(subtraction_CC)
       RETURN       
    END IF

    IF (this%frame /= that%frame) THEN
       CALL errorMessage("CartesianCoordinates / - (subtraction)", &
            "Coordinates are not given in the same frame of reference.", 1)
       CALL NULLIFY(subtraction_CC)
       RETURN       
    END IF

    this_ = copy(this)
    that_ = copy(that)

    IF (equal(this_%t,that_%t)) THEN
       subtraction_CC = copy(this_)
       subtraction_CC%position = this_%position - that_%position
       subtraction_CC%velocity = this_%velocity - that_%velocity
    ELSE
       CALL errorMessage("CartesianCoordinates / - (subtraction)", &
            "Coordinatess at different times can not be used in subtraction.", 1)
       CALL NULLIFY(subtraction_CC)
       RETURN          
    END IF

  END FUNCTION subtraction_CC





END MODULE CartesianCoordinates_cl
