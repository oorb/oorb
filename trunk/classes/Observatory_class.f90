!====================================================================!
!                                                                    !
! Copyright 2009 Mikael Granvik, Jenni Virtanen, Karri Muinonen,     !
!                Teemu Laakso, Dagmara Oszkiewicz                    !
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
!! Type and routines for an observatory.
!!
!! @author  MG
!! @version 2008-08-12
!!
MODULE Observatory_cl

  USE Base_cl

  IMPLICIT NONE

  PRIVATE :: new_Obsy
  PRIVATE :: nullify_Obsy
  PRIVATE :: copy_Obsy
  PRIVATE :: exist_Obsy
  PRIVATE :: equal_Obsy
  PRIVATE :: getCode_Obsy
  PRIVATE :: getName_Obsy
  PRIVATE :: getPosition_Obsy

  TYPE Observatory
     PRIVATE
     CHARACTER(len=96)            :: name           =  ""
     CHARACTER(len=OBSY_CODE_LEN) :: code           =  ""
     REAL(bp), DIMENSION(3)       :: position       =  0.0_bp
     LOGICAL                      :: is_initialized = .FALSE.
  END TYPE Observatory

  INTERFACE NEW
     MODULE PROCEDURE new_Obsy
  END INTERFACE

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_Obsy
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_Obsy
  END INTERFACE

  INTERFACE exist
     MODULE PROCEDURE exist_Obsy
  END INTERFACE

  INTERFACE equal
     MODULE PROCEDURE equal_Obsy
  END INTERFACE

  INTERFACE getCode
     MODULE PROCEDURE getCode_Obsy
  END INTERFACE

  INTERFACE getName
     MODULE PROCEDURE getName_Obsy
  END INTERFACE

  INTERFACE getPosition
     MODULE PROCEDURE getPosition_Obsy
  END INTERFACE


CONTAINS




  !! *Description*:
  !!
  !! Initializes a new object based on given information. The position
  !! contains the body-fixed Cartesian geocentric equatorial
  !! coordinates of the observatory given in AUs.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Obsy(this, code, name, position)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(out)    :: this
    CHARACTER(len=*), INTENT(in)       :: code
    CHARACTER(len=*), INTENT(in)       :: name
    REAL(bp), DIMENSION(3), INTENT(in) :: position

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (LEN_TRIM(code) > OBSY_CODE_LEN) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / new", &
            "Observatory code ("// TRIM(code) // &
            ")too long.", 1)
       RETURN
    END IF

    this%code           =  code
    this%name           =  name
    this%position       =  position
    this%is_initialized = .TRUE.

  END SUBROUTINE new_Obsy





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_Obsy(this)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(inout) :: this

    this%code           =  ""
    this%name           =  ""
    this%position       =  0.0_bp
    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_Obsy





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_Obsy(this)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(in) :: this
    TYPE (Observatory)             :: copy_Obsy

    copy_Obsy%code           = this%code
    copy_Obsy%name           = this%name
    copy_Obsy%position       = this%position
    copy_Obsy%is_initialized = this%is_initialized

  END FUNCTION copy_Obsy





  !! *Description*:
  !!
  !! Returns the status of this object, i.e. whether
  !! it exists or not.
  !!
  LOGICAL FUNCTION exist_Obsy(this)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(in) :: this

    exist_Obsy = this%is_initialized

  END FUNCTION exist_Obsy





  LOGICAL FUNCTION equal_Obsy(this, that)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(in) :: this
    TYPE (Observatory), INTENT(in) :: that

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / equal", &
            "1st object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. that%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / equal", &
            "2nd object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Observatory code:
    IF (.NOT.(this%code == that%code)) THEN
       equal_Obsy = .FALSE.
       RETURN
    END IF

    ! Position wrt. the geocenter:
    IF (ANY(ABS(this%position - that%position) > EPSILON(this%position(1)))) THEN
       equal_Obsy = .FALSE.
       RETURN
    END IF

    ! Assuming that all of the above comparisons are true,
    ! it can be concluded that the two objects are the same:
    equal_Obsy = .TRUE.

  END FUNCTION equal_Obsy





  !! *Description*:
  !!
  !! Returns observatory code as a character string.
  !!
  !! Returns error.
  !!
  CHARACTER(len=OBSY_CODE_LEN) FUNCTION getCode_Obsy(this)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / getCode", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getCode_Obsy = this%code

  END FUNCTION getCode_Obsy





  !! *Description*:
  !!
  !! Returns name of observatory as a character string.
  !!
  !! Returns error.
  !!
  CHARACTER(len=128) FUNCTION getName_Obsy(this)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / getName", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getName_Obsy = this%name

  END FUNCTION getName_Obsy





  !! *Description*:
  !!
  !! Returns geocentric coordinates of the observatory.
  !!
  !! Returns error.
  !!
  FUNCTION getPosition_Obsy(this)

    IMPLICIT NONE
    TYPE (Observatory), INTENT(in) :: this
    REAL(bp), DIMENSION(3)         :: getPosition_Obsy

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Observatory / getPosition", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getPosition_Obsy = this%position

  END FUNCTION getPosition_Obsy





END MODULE Observatory_cl

