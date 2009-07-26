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
!! Type and routines for logical unit handling.
!!
!! @see File_class
!! 
!! @author  MG, JV
!! @version 2008-08-12
!!
MODULE Unit_cl

  USE Base_cl
  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: min_lu = 10 !! Minimum logical unit.
  INTEGER, PARAMETER, PRIVATE :: max_lu = 99 !! Maximum logical unit.
  INTEGER, SAVE, PRIVATE :: next_lu = min_lu ! Next supposedly free unit.

  PRIVATE :: new_U
  PRIVATE :: nullify_U
  PRIVATE :: copy_U
  PRIVATE :: exist_U
  PRIVATE :: safeLogicalUnit
  PRIVATE :: isOpen

  TYPE Unit
     PRIVATE
     INTEGER :: lu                       ! Logical unit.
     LOGICAL :: is_initialized = .FALSE. ! Has this unit been initialized? 
     !                                     (true=yes, false=no)
  END TYPE Unit

  !! Initializes a new instance of Unit.
  INTERFACE NEW
     MODULE PROCEDURE new_U
  END INTERFACE

  !! Nullifies a Unit-object, e.g. closes the unit.
  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_U
  END INTERFACE

  !! Makes a copy of a Unit-object.
  INTERFACE copy
     MODULE PROCEDURE copy_U
  END INTERFACE

  !!
  INTERFACE exist
     MODULE PROCEDURE exist_U
  END INTERFACE

  !! Returns the logical unit number of a Unit instance.
  INTERFACE getUnit
     MODULE PROCEDURE getUnit_U
  END INTERFACE





CONTAINS





  !! *Description*:
  !!
  !! Initializes a new object and finds a free 
  !! logical unit. Returns _error_ =
  !!     - _false_, if a free unit is found.
  !!     - _true_, if this unit has already been opened, 
  !!       or a free logical unit ca not be found.
  !!
  !! *Usage*:
  !!
  !! type (Unit) :: myunit
  !!
  !! ...
  !!
  !! call new(myunit) 
  !!
  SUBROUTINE new_U(this)

    IMPLICIT NONE
    TYPE (Unit), INTENT(inout) :: this
    INTEGER                    :: lu

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Unit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    lu = next_lu
    CALL safeLogicalUnit(lu)
    IF (error) THEN
       CALL errorMessage("Unit / new", &
            "TRACE BACK", 1)
       RETURN
    END IF

    this%lu             = lu
    this%is_initialized = .TRUE.


  END SUBROUTINE new_U





  !! *Description*:
  !!
  !! Closes logical unit. Returns _error_ =
  !!     - _false_, if the unit is properly closed, or if the unit 
  !!                has not been opened at all. 
  !!     - _true_, if an error occurs during the closing procedure. 
  !!
  !! *Usage*:
  !!
  !! type (Unit) :: myunit
  !!
  !! ...
  !!
  !! call nullify(myunit)
  !!
  SUBROUTINE nullify_U(this)

    IMPLICIT NONE
    TYPE (Unit), INTENT(inout) :: this
    INTEGER                    :: err

    IF (isOpen(this)) THEN
       CLOSE(this%lu, iostat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Unit / nullify", &
               "Could not close unit.", 1)
          RETURN
       END IF
    END IF

    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_U





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_U(this) RESULT(a_copy)

    IMPLICIT NONE
    TYPE (Unit), INTENT(in)  :: this
    TYPE (Unit)              :: a_copy

    a_copy%lu             = this%lu
    a_copy%is_initialized = this%is_initialized

  END FUNCTION copy_U





  !! *Description*:
  !!
  !! Returns the status of the object, i.e. whether
  !! it exists or not.
  !!
  LOGICAL FUNCTION exist_U(this)

    IMPLICIT NONE
    TYPE (Unit), INTENT(in)  :: this

    exist_U = this%is_initialized

  END FUNCTION exist_U





  !! *Description*:
  !!
  !! Returns
  !!     - logical unit number, if the instance has 
  !!                            been initialized.
  !!     - -1, if the unit has not been initialized.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! type (Unit) :: myunit
  !!
  !! integer :: logical_unit
  !!
  !! ...
  !!
  !! logical_unit = getUnit(myunit) 
  !!
  INTEGER FUNCTION getUnit_U(this)

    IMPLICIT NONE
    TYPE (Unit), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Unit / getUnit", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getUnit_U = this%lu

  END FUNCTION getUnit_U





  !! *Description*:
  !!
  !! Examines a Unit instance. Returns
  !!     - _true_, if the unit is initialized and in use.
  !!     - _false_, if the unit is not used or not initialized. 
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! type (Unit) :: myunit
  !!
  !! logical :: open
  !!
  !! ...
  !!
  !! open = isOpen(myunit)
  !!
  LOGICAL FUNCTION isOpen(this)

    IMPLICIT NONE
    TYPE (Unit), INTENT(in) :: this
    LOGICAL                 :: used
    INTEGER                 :: err

    INQUIRE(this%lu, opened=used, iostat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Unit / isOpen", &
            "Inquiry returned error.", 1)
       RETURN
    ELSE
       IF (this%is_initialized) THEN
          isOpen = used
       ELSE
          isOpen = .FALSE.
       END IF
    END IF

  END FUNCTION isOpen





  !! *Description*:
  !!
  !! Examines a logical unit number. If it is unused,
  !! the value is returned unchanged. If it is in use, 
  !! the routine finds an unused logical unit number 
  !! between parameters "min_lu" and "max_lu" (incl.) 
  !! and returns that value. Returns _error_ =
  !!     - _false_, if the input unit is ok, or a free 
  !!                unit was found.
  !!     - _true_,  if the input unit is not ok, and a 
  !!                free unit ca not be found.
  !!
  !! *Usage*:
  !!
  !! integer :: unit_nr = 10
  !!
  !! ...
  !!
  !! call safeLogicalUnit(unit_nr) 
  !!
  SUBROUTINE safeLogicalUnit(lu)

    IMPLICIT NONE
    INTEGER, INTENT(inout) :: lu
    INTEGER                :: count, err
    LOGICAL                :: done, used

    done = .FALSE.
    count = min_lu
    DO WHILE (.NOT. done)
       ! Figure out whether this unit is taken or not:
       INQUIRE(unit=lu, opened=used, iostat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Unit / safeLogicalUnit", &
               "Inquiry failed.", 1)
          RETURN
       END IF
       IF (used) THEN
          count = count + 1
          ! If more than max_lu units have been tried,
          ! every available unit has been tried at least once.
          ! A free unit could not be found: 
          IF (count > max_lu) THEN
             error = .TRUE.
             CALL errorMessage("Unit / safeLogicalUnit", &
                  "Could not find a free logical unit.", 1)
             RETURN
          END IF
          ! Try next supposedly free unit:
          lu = next_lu
          ! Push next supposedly free unit one unit further:
          next_lu = next_lu + 1
          ! Back to beginning if top is reached:
          IF (next_lu > max_lu) next_lu = min_lu
       ELSE
          done = .TRUE.
       END IF
    END DO

  END SUBROUTINE safeLogicalUnit





  !! *Description*:
  !!
  !! Sets the logical unit of a Unit instance to 
  !! the debug logical unit. Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if the unit has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (Unit) :: myunit
  !!
  !! ...
  !!
  !! call setDebugUnit(myunit)
  !! 
  SUBROUTINE setDebugUnit(this)

    IMPLICIT NONE
    TYPE (Unit), INTENT(inout) :: this

    IF (isOpen(this)) THEN
       error = .TRUE.
       CALL errorMessage("Unit / setDebugUnit", &
            "Unit has already been opened.", 1)
       RETURN
    END IF

    this%lu = debug

  END SUBROUTINE setDebugUnit





END MODULE Unit_cl
