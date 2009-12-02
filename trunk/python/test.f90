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
!! Test program for pyoorb module.
!!
!! @author  MG, JM
!! @version 2009-12-01
!!
PROGRAM test

  USE pyoorb
  IMPLICIT NONE
  CHARACTER(len=4) :: in_obscode
  INTEGER, PARAMETER :: in_norb = 3
  INTEGER, PARAMETER :: in_ndate = 2
  REAL(8), DIMENSION(in_norb,6,6) :: in_covariances
  REAL(8), DIMENSION(in_norb,12) :: in_orbits
  REAL(8), DIMENSION(in_norb,in_ndate,11) :: out_ephems
  REAL(8), DIMENSION(in_ndate,2) :: in_date_ephems
  INTEGER :: error_code, i, j

  ! Initialize 
  error_code = 0
  CALL oorb_init(error_code=error_code, info_verbosity=1)
  IF (error_code /= 0) THEN
     WRITE(stderr,*) "Error in oorb_init. Code: ", error_code
     STOP
  END IF

  ! 2009 TL4 on 2009-11-12 (74 obs in 2009) from (Horizons & AstDyS)
  in_orbits(1,1) = 1.0_8
  in_orbits(1,2:7) = (/ &
       1.853233422926951E+00_8, &
       3.871022198610788E-01_8, &
       6.474397461946572E+00_8, &
       2.122138159851688E+02_8, &
       1.602711715213041E+02_8, &
       3.372535412461616E+01_8 /)
  in_orbits(1,4:7) = in_orbits(1,4:7)*rad_deg
  in_orbits(1,8) = 3.0_8
  in_orbits(1,9) = 55200.0_8
  in_orbits(1,10) = 3.0_8
  in_orbits(1,11:12) = (/ 21.541_8, 0.15_8 /)
  ! Copy the same orbit for testing purposes
  in_orbits(2,:) = in_orbits(1,:)
  in_orbits(2,1) = 2.0_8
  in_orbits(3,:) = in_orbits(1,:)
  in_orbits(3,1) = 3.0_8

  ! From AstDyS
  in_covariances(1,1,1:6) = (/ 1.59357530E-06_8, 4.54646999E-07_8, 5.59222797E-06_8, &
       -3.87304158E-06_8, -3.54135866E-07_8, -4.31574921E-05_8 /)
  in_covariances(1,2,1:6) = (/ 4.54646999E-07_8, 1.29710940E-07_8, 1.59544037E-06_8, &
       -1.10495903E-06_8, -1.00707885E-07_8, -1.23129696E-05_8 /)
  in_covariances(1,3,1:6) = (/ 5.59222797E-06_8, 1.59544037E-06_8, 1.96320645E-05_8, &
       -1.36008250E-05_8, -1.34577989E-06_8, -1.51404682E-04_8 /)
  in_covariances(1,4,1:6) = (/ -3.87304158E-06_8, -1.10495903E-06_8, -1.36008250E-05_8, &
       9.42766927E-06_8, 9.60529697E-07_8, 1.04844563E-04_8 /)
  in_covariances(1,5,1:6) = (/ -3.54135866E-07_8, -1.00707885E-07_8, -1.34577989E-06_8, &
       9.60529697E-07_8, 2.72099434E-06_8, 8.49016456E-06_8 /)
  in_covariances(1,6,1:6) = (/ -4.31574921E-05_8, -1.23129696E-05_8, -1.51404682E-04_8, &
       1.04844563E-04_8, 8.49016456E-06_8, 1.16925907E-03_8 /)
  DO i=1,6
     IF (i >= 3) THEN
        in_covariances(1,i,:) = in_covariances(1,i,:)*rad_deg
        in_covariances(1,:,i) = in_covariances(1,:,i)*rad_deg
     END IF
  END DO
  ! Copy the covariance matrix for testing purposes
  DO i=2,3
     in_covariances(i,:,:) = in_covariances(1,:,:)
  END DO

  ! Observatory code
  in_obscode = "500"

  ! Dates for output ephemerides
  in_date_ephems(1,1:2) = (/ 55148.0_8, 1.0_8 /)
  in_date_ephems(2,1:2) = (/ 55150.0_8, 1.0_8 /)

  ! Compute ephemerides with covariance
  CALL oorb_ephemeris_covariance(in_norb, &
       in_orbits,                         &
       in_covariances,                    &
       in_obscode,                        &
       in_ndate,                          & 
       in_date_ephems,                    &
       out_ephems,                        &
       error_code)
  IF (error_code /= 0) THEN
     WRITE(stderr,*) "Error in oorb_ephemeris. Code: ", error_code
     STOP
  END IF
  WRITE(stdout,"(12(2X,A12))") "ID", "Delta", "RA", "Dec", &
       "Mag", "MJD", "Timescale", "sigma_RA", "sigma_Dec", &
       "Uncert smaa", "Uncert smia", "Uncert PA"
  DO i=1,SIZE(out_ephems,dim=1)
     DO j=1,SIZE(out_ephems,dim=2)
        WRITE(stdout,"(12(2X,E12.5))") in_orbits(i,1), out_ephems(i,j,:)
     END DO
  END DO

  ! Compute ephemerides without covariance
  out_ephems = -1.0_8
  CALL oorb_ephemeris(in_norb, &
       in_orbits,              &
       in_obscode,             &
       in_ndate,               & 
       in_date_ephems,         &
       out_ephems(:,:,1:6),    &
       error_code)
  IF (error_code /= 0) THEN
     WRITE(stderr,*) "Error in oorb_ephemeris. Code: ", error_code
     STOP
  END IF
  WRITE(stdout,"(7(2X,A12))") "ID", "Delta", "RA", "Dec", &
       "Mag", "MJD", "Timescale"
  DO i=1,SIZE(out_ephems,dim=1)
     DO j=1,SIZE(out_ephems,dim=2)
        WRITE(stdout,"(7(2X,E12.5))") in_orbits(i,1), out_ephems(i,j,1:6)
     END DO
  END DO

!!$ Ephemeris for 2009 TL4 at epoch 55148 MJD (from AstDyS)
!!$ Right Ascension:      17.801   deg
!!$ Declination:          -0.08723 deg
!!$ V Magnitude:          19.98    mag
!!$ Solar Elongation:   -145.14    deg
!!$ Phase:                28.48    deg
!!$ Galactic Latitude:   -62.54    deg
!!$ Distance from Earth:   0.23    AU
!!$ Distance from Sun:     1.186   AU
!!$ Apparent motion 
!!$ Rate:                  0.848   deg/day
!!$ Direction:            95.121   deg
!!$ Uncertainty Ellipse 
!!$ Semimajor axis:        0.021   arcmin
!!$ Semiminor axis:        0.017   arcmin
!!$ Orientation:        -161.75    deg
!!$ From JPL:
!!$ RA  3sigma unc:        3.097   arcsec
!!$ Dec 3sigma unc:        3.672   arcsec 

  CALL oorb_memfree()

END PROGRAM test
