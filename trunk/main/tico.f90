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
!! Mikael Granvik 2009-06-18
!!
PROGRAM tico

  USE cl_options
  USE Base_cl
  USE Time_cl

  IMPLICIT NONE
  TYPE (Time) :: t
  REAL(bp) :: mjd_in, mjd_out

  ! Set path to data files:
  CALL setAccessToDataFiles()

  ! Input
  IF (get_cl_option("--tai_in=", .FALSE.)) THEN
     mjd_in = get_cl_option("--tai_in=", 0.0_bp)
     CALL NEW(t, mjd_in, "TAI")
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
  ELSE IF (get_cl_option("--utc_in=", .FALSE.)) THEN
     mjd_in = get_cl_option("--utc_in=", mjd_in)
     CALL NEW(t, mjd_in, "UTC")
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (10)", 1)
        STOP
     END IF
  ELSE IF (get_cl_option("--tt_in=", .FALSE.)) THEN
     mjd_in = get_cl_option("--tt_in=", mjd_in)
     CALL NEW(t, mjd_in, "TT")
     IF (error) THEN
        CALL errorMessage("tico", &
             "TRACE BACK (15)", 1)
        STOP
     END IF
  ELSE
     CALL errorMessage("tico", &
          "Input time/timescale not specified.", 1)     
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
          "Output timescale not specified.", 1)
  END IF
  CALL NULLIFY(t)

END PROGRAM tico
