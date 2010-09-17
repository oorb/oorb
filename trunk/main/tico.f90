!====================================================================!
!                                                                    !
! Copyright 2002,2003,2004,2005,2006,2007,2008,2009,2010             !
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
!!
!! DESCRIPTION
!!
!! tico performs conversions between timescales (TAI, TT, UTC).  
!!
!!
!! USAGE
!!
!!  For example, to convert MJD 53653.36413 TAI to UTC do
!!
!!     tico --tai_in=53653.36413 --utc_out
!!
!!
!! REQUIREMENTS
!!
!! Environment variable $OORB_DATA has to point to a directory which
!! contains files ET-UT.dat and TAI-UTC.dat. Make sure these files are
!! updated regularly.  
!!
!!
!! AUTHOR
!!
!! Mikael Granvik 2010-09-17
!!
PROGRAM tico

  USE cl_options
  USE Base_cl
  USE Time_cl

  IMPLICIT NONE
  TYPE (Time) :: t
  CHARACTER(len=64) :: str
  REAL(bp) :: mjd_in, mjd_out

  ! Set path to data files:
  CALL setAccessToDataFiles()

  ! Input
  IF (get_cl_option("--tai_in=", .FALSE.)) THEN
     str = get_cl_option("--tai_in=", "")
     CALL readDate(str, "TAI", t)
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
  ELSE IF (get_cl_option("--utc_in=", .FALSE.)) THEN
     str = get_cl_option("--utc_in=", "")
     CALL readDate(str, "UTC", t)
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (10)", 1)
        STOP
     END IF
  ELSE IF (get_cl_option("--tt_in=", .FALSE.)) THEN
     str = get_cl_option("--tt_in=", "")
     CALL readDate(str, "TT", t)
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (15)", 1)
        STOP
     END IF
  ELSE
     CALL errorMessage("tico", &
          "Input time/timescale not specified. Use the '--[tai|utc|tt]_in=' option.", 1)     
     STOP
  END IF

  ! Output
  IF (get_cl_option("--utc_out", .FALSE.)) THEN
     mjd_out = getMJD(t, "UTC")
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (20)", 1)
        STOP
     END IF
     WRITE(stdout,'(F15.8)') mjd_out
  ELSE IF (get_cl_option("--tai_out", .FALSE.)) THEN
     mjd_out = getMJD(t, "TAI")
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (25)", 1)
        STOP
     END IF
     WRITE(stdout,'(F15.8)') mjd_out
  ELSE IF (get_cl_option("--tt_out", .FALSE.)) THEN
     mjd_out = getMJD(t, "TT")
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (30)", 1)
        STOP
     END IF
     WRITE(stdout,'(F15.8)') mjd_out
  ELSE
     CALL errorMessage("tico", &
          "Output timescale not specified. Use the '--[tai|utc|tt]_out' option.", 1)
  END IF
  CALL NULLIFY(t)

CONTAINS

  SUBROUTINE readDate(str, timescale, t)

    IMPLICIT NONE
    TYPE (Time), INTENT(out) :: t
    CHARACTER(len=*), INTENT(in) :: str, timescale
    REAL(bp) :: mjd_in, day
    INTEGER, DIMENSION(0:2) :: indx_arr
    INTEGER :: i, j, year, month

    ! Count hyphens to figure out if this is a calendar date or MJD:
    indx_arr = 0
    i = 0
    j = 1
    DO WHILE (j > 0)
       j = INDEX(str(indx_arr(i)+1:),"-")
       IF (j > 0) THEN
          indx_arr(i+1) = indx_arr(i) + j
          i = i + 1
       END IF
    END DO
    IF (i > 1) THEN
       CALL toInt(str(1:indx_arr(1)-1), year, error)
       CALL toInt(str(indx_arr(1)+1:indx_arr(2)-1), month, error)
       CALL toReal(TRIM(str(indx_arr(2)+1:)), day, error)
       CALL NEW(t, year, month, day, TRIM(timescale))
    ELSE
       CALL toReal(str, mjd_in, error)
       CALL NEW(t, mjd_in, TRIM(timescale))
    END IF
    IF (error) THEN
       CALL errorMessage("tico / readDate", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END SUBROUTINE readDate


END PROGRAM tico
