!====================================================================!
!                                                                    !
! Copyright 2002-2016,2017                                           !
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
!! *Module*description*:
!!
!! Contains generic functions that return command line options.
!!
!! @author  MG
!! @version 2017-10-26
!!
MODULE cl_options

  USE parameters
  IMPLICIT NONE
  INTEGER, PARAMETER :: CL_LEN = 4096

  PRIVATE :: get_cl_option_char
  PRIVATE :: get_cl_option_real8
  PRIVATE :: get_cl_option_real
  PRIVATE :: get_cl_option_integer
  PRIVATE :: get_cl_option_logical

  INTERFACE get_cl_option
     MODULE PROCEDURE get_cl_option_char
     MODULE PROCEDURE get_cl_option_real8
     MODULE PROCEDURE get_cl_option_real
     MODULE PROCEDURE get_cl_option_integer
     MODULE PROCEDURE get_cl_option_logical
  END INTERFACE


CONTAINS





  CHARACTER(len=256) FUNCTION get_cl_option_char(option, default) RESULT(value)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: option, default
    CHARACTER(len=CL_LEN) :: cl = ""
    INTEGER :: i, di, j, err

    CALL get_command_line(cl)
    i = INDEX(TRIM(cl), TRIM(option))
    di = LEN_TRIM(option)
    IF (i /= 0) THEN
       READ(cl(i+di:),'(a)',iostat=err) value
       IF (err /= 0) value = TRIM(default)
    ELSE
       value = TRIM(default)
    END IF
    ! Remove blanks in the beginning:
    j = 0
    DO WHILE (value(1:1) == " ")
       j = j + 1
       IF (j>256) RETURN
       value = TRIM(value(2:))
    END DO
    ! Return empty value, if no value is given:
    IF (value(2:2) == "-") THEN
       value = ""
       RETURN
    END IF
    ! Remove a blank found in the middle and everything
    ! coming after it:
    i = INDEX(TRIM(value)," ")
    IF (i /= 0) value = TRIM(value(1:i-1))

  END FUNCTION get_cl_option_char





  REAL(rprec8) FUNCTION get_cl_option_real8(option, default) RESULT(value)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: option
    REAL(rprec8), INTENT(in) :: default
    CHARACTER(len=CL_LEN) :: cl = ""
    INTEGER :: i, di, err

    CALL get_command_line(cl)
    i = INDEX(TRIM(cl), TRIM(option))
    di = LEN_TRIM(option) - 1
    IF (i /= 0) THEN
       READ(cl(i+di+1:),*,iostat=err) value
       IF (err /= 0) value = default
    ELSE
       value = default
    END IF

  END FUNCTION get_cl_option_real8





  REAL FUNCTION get_cl_option_real(option, default) RESULT(value)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: option
    REAL, INTENT(in) :: default
    CHARACTER(len=CL_LEN) :: cl = ""
    INTEGER :: i, di, err

    CALL get_command_line(cl)
    i = INDEX(TRIM(cl), TRIM(option))
    di = LEN_TRIM(option) - 1
    IF (i /= 0) THEN
       READ(cl(i+di+1:),*,iostat=err) value
       IF (err /= 0) value = default
    ELSE
       value = default
    END IF

  END FUNCTION get_cl_option_real





  INTEGER FUNCTION get_cl_option_integer(option, default) RESULT(value)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: option
    INTEGER, INTENT(in) :: default
    CHARACTER(len=CL_LEN) :: cl = ""
    INTEGER :: i, di, err

    CALL get_command_line(cl)
    i = INDEX(TRIM(cl), TRIM(option))
    di = LEN_TRIM(option) - 1
    IF (i /= 0) THEN
       READ(cl(i+di+1:),*,iostat=err) value
       IF (err /= 0) value = default
    ELSE
       value = default
    END IF

  END FUNCTION get_cl_option_integer





  ! Returns true if the option is found, and the given default value otherwise.
  !
  LOGICAL FUNCTION get_cl_option_logical(option, default) RESULT(value)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: option
    LOGICAL, INTENT(in) :: default
    CHARACTER(len=CL_LEN) :: cl = " "
    INTEGER :: i

    CALL get_command_line(cl)
    i = INDEX(TRIM(cl), TRIM(option))
    IF (i /= 0) THEN
       value = .TRUE.
    ELSE
       value = default
    END IF

  END FUNCTION get_cl_option_logical





  SUBROUTINE get_command_line(cl)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(out) :: cl
    EXTERNAL :: getarg
    INTEGER, EXTERNAL :: iargc
    CHARACTER(len=256) :: arg
    INTEGER :: iarg, narg

    cl = " "
    ! If an error occurs during compilation, comment the section
    ! between the ='s in/out and the section between the +'s out/in:
    !=======================================
!!$    narg = iargc()
!!$    DO iarg=1, narg
!!$       CALL getarg(iarg, arg)
!!$       cl = TRIM(cl) // " " // TRIM(arg)
!!$    END DO
    !=======================================
    !+++++++++++++++++++
    CALL get_command(cl)
    !+++++++++++++++++++

  END SUBROUTINE get_command_line





END MODULE cl_options

