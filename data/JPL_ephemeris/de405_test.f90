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
!!*Program*description*:
!!
!! This program tests the 'de405.dat' JPL planetary ephemeris file and
!! the routine 'JPL_ephemeris' that extracts interpolated positions.
!! Positions and velocities of the Sun, the 8 planets, the Moon, and
!! Pluto are included.
!!
!! @author  MG
!! @version 2009-08-19
!!
PROGRAM de405_tester

  USE parameters
  USE planetary_data

  CHARACTER(len=256) :: str
  REAL(rprec8), DIMENSION(:,:), POINTER :: ephemeris
  REAL(rprec8) :: mjd_tt, jd_tt, value
  INTEGER :: i, ntarget, ncenter, ncoord, err
  LOGICAL :: error

  error = .FALSE.
  OPEN(12,file="testpo.405",status="old")
  str = ""
  DO WHILE (str(1:3) /= "EOT")
     READ(12,"(A)") str
  END DO
  CALL JPL_ephemeris_init(error, "de405.dat")
  IF (error) THEN
     WRITE(0,*) "***** INITIALIZATION ERROR OCCURRED *****"     
     STOP
  END IF
  DO
     READ(12,*,iostat=err) i, str, jd_tt, ntarget, ncenter, ncoord, value
     IF (err < 0) THEN
        EXIT
     ELSE IF (err > 0) THEN
        IF (error) THEN
           WRITE(0,*) "***** READ ERROR OCCURRED *****"     
           STOP
        END IF
     END IF
     IF (ncenter == 0 .OR. ncenter > 11 .OR. ntarget > 11) THEN
        CYCLE
     END IF
     mjd_tt = jd_tt - 2400000.5_rprec8
     ephemeris => JPL_ephemeris(mjd_tt, ntarget, ncenter, error)
     IF (ABS(ephemeris(1,ncoord)-value) > 10.0E-13_rprec8) THEN
        WRITE(0,*) "***** COMPARISON ERROR OCCURRED *****"
        STOP
     END IF
     DEALLOCATE(ephemeris)
  END DO
  CLOSE(12)
  CALL JPL_ephemeris_nullify()
  WRITE(*,*) "*********** EVERYTHING OK ***************"

END PROGRAM de405_tester
