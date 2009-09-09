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
!! *Module*description*: 
!!  
!! Random number generators. 
!!
!! @see SphericalCoordinates_class 
!!  
!! @author  MG
!! @version 2008-08-12
!!
MODULE random

  USE parameters
  IMPLICIT NONE
  INTEGER, PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed 
  INTEGER :: idum, idum_prm = -1
  LOGICAL :: first_ran = .TRUE.

  INTERFACE ran_pmm
     MODULE PROCEDURE ran_pmm_r4
     MODULE PROCEDURE ran_pmm_r8
     MODULE PROCEDURE ran_pmm_r16
  END INTERFACE

  INTERFACE randomNumber
     MODULE PROCEDURE randomNumber_single_r4
     MODULE PROCEDURE randomNumber_single_r8
     MODULE PROCEDURE randomNumber_single_r16
     MODULE PROCEDURE randomNumber_array_r4
     MODULE PROCEDURE randomNumber_2array_r4
     MODULE PROCEDURE randomNumber_array_r8
     MODULE PROCEDURE randomNumber_2array_r8
     MODULE PROCEDURE randomNumber_array_r16
  END INTERFACE

  INTERFACE randomGaussian
     MODULE PROCEDURE randomGaussian_single_r4
     MODULE PROCEDURE randomGaussian_single_r8
     MODULE PROCEDURE randomGaussian_single_r16
     MODULE PROCEDURE randomGaussian_array_r4
     MODULE PROCEDURE randomGaussian_array_r8
     MODULE PROCEDURE randomGaussian_array_r16
  END INTERFACE


CONTAINS



  SUBROUTINE initializeRandomNumberGenerator(idum)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: idum

    idum_prm = idum

  END SUBROUTINE initializeRandomNumberGenerator





  !! Minimal random number generator of Park and Miller combined with a
  !! Marsaglia shift sequence. Returns a uniform random deviate between
  !! 0.0 and 1.0 (exclusive of the endpoint values). This fully portable,
  !! scalar generator has the traditional (not Fortran 90)
  !! calling sequence with a random deviate as the returned function value:
  !! call with idum a negative integer to initialize; thereafter, do not
  !! alter idum except to reinitialize. The period of this generator is
  !! about 3.1 × 10^18.
  !!
  SUBROUTINE ran_pmm_r4(idum, ran)

    IMPLICIT NONE 
    INTEGER, INTENT(inout) :: idum 
    REAL(rprec4), INTENT(out) :: ran
    REAL(rprec4), SAVE :: am
    INTEGER, SAVE :: ix=-1, iy=-1, k

    IF (idum <= 0 .OR. iy < 0) THEN 
       ! Initialize.
       am=NEAREST(1.0_rprec4,-1.0_rprec4)/IM
       iy=IOR(IEOR(888889999,ABS(idum)),1)
       ix=IEOR(777755555,ABS(idum))
       ! Set idum positive.
       idum=ABS(idum)+1
    END IF

    ! Marsaglia shift sequence with period 2^32 - 1.
    ix=IEOR(ix,ISHFT(ix,13))
    ix=IEOR(ix,ISHFT(ix,-17))
    ix=IEOR(ix,ISHFT(ix,5))

    ! Park-Miller sequence by Schrage's method, period 2^31 - 2
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    IF (iy < 0) THEN
       iy=iy+IM
    END IF

    !Combine the two generators with masking to ensure nonzero value
    ran = am*IOR(IAND(IM,IEOR(ix,iy)),1)

  END SUBROUTINE ran_pmm_r4





  !! Minimal random number generator of Park and Miller combined with a
  !! Marsaglia shift sequence. Returns a uniform random deviate between
  !! 0.0 and 1.0 (exclusive of the endpoint values). This fully portable,
  !! scalar generator has the traditional (not Fortran 90)
  !! calling sequence with a random deviate as the returned function value:
  !! call with idum a negative integer to initialize; thereafter, do not
  !! alter idum except to reinitialize. The period of this generator is
  !! about 3.1 × 10^18.
  !!
  SUBROUTINE ran_pmm_r8(idum, ran)

    IMPLICIT NONE 
    INTEGER, INTENT(inout) :: idum 
    REAL(rprec8), INTENT(out) :: ran
    REAL(rprec8), SAVE :: am
    INTEGER, SAVE :: ix=-1, iy=-1, k

    IF (idum <= 0 .OR. iy < 0) THEN 
       ! Initialize.
       am=NEAREST(1.0_rprec8,-1.0_rprec8)/IM
       iy=IOR(IEOR(888889999,ABS(idum)),1)
       ix=IEOR(777755555,ABS(idum))
       ! Set idum positive.
       idum=ABS(idum)+1
    END IF

    ! Marsaglia shift sequence with period 2^32 - 1.
    ix=IEOR(ix,ISHFT(ix,13))
    ix=IEOR(ix,ISHFT(ix,-17))
    ix=IEOR(ix,ISHFT(ix,5))

    ! Park-Miller sequence by Schrage's method, period 2^31 - 2
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    IF (iy < 0) THEN
       iy=iy+IM
    END IF

    !Combine the two generators with masking to ensure nonzero value
    ran = am*IOR(IAND(IM,IEOR(ix,iy)),1)

  END SUBROUTINE ran_pmm_r8




  !! Minimal random number generator of Park and Miller combined with a
  !! Marsaglia shift sequence. Returns a uniform random deviate between
  !! 0.0 and 1.0 (exclusive of the endpoint values). This fully portable,
  !! scalar generator has the traditional (not Fortran 90)
  !! calling sequence with a random deviate as the returned function value:
  !! call with idum a negative integer to initialize; thereafter, do not
  !! alter idum except to reinitialize. The period of this generator is
  !! about 3.1 × 10^18.
  !!
  SUBROUTINE ran_pmm_r16(idum, ran)

    IMPLICIT NONE 
    INTEGER, INTENT(inout) :: idum 
    REAL(rprec16), INTENT(out) :: ran
    REAL(rprec16), SAVE :: am
    INTEGER, SAVE :: ix=-1, iy=-1, k

    IF (idum <= 0 .OR. iy < 0) THEN 
       ! Initialize.
       am=NEAREST(1.0_rprec16,-1.0_rprec16)/IM
       iy=IOR(IEOR(888889999,ABS(idum)),1)
       ix=IEOR(777755555,ABS(idum))
       ! Set idum positive.
       idum=ABS(idum)+1
    END IF

    ! Marsaglia shift sequence with period 2^32 - 1.
    ix=IEOR(ix,ISHFT(ix,13))
    ix=IEOR(ix,ISHFT(ix,-17))
    ix=IEOR(ix,ISHFT(ix,5))

    ! Park-Miller sequence by Schrage's method, period 2^31 - 2
    k=iy/IQ
    iy=IA*(iy-k*IQ)-IR*k
    IF (iy < 0) THEN
       iy=iy+IM
    END IF

    !Combine the two generators with masking to ensure nonzero value
    ran = am*IOR(IAND(IM,IEOR(ix,iy)),1)

  END SUBROUTINE ran_pmm_r16





  SUBROUTINE randomNumber_single_r4(ran)

    IMPLICIT NONE
    REAL(rprec4), INTENT(out) :: ran

    IF (first_ran) THEN
       idum = idum_prm
       first_ran = .FALSE.
    END IF
    CALL ran_pmm(idum, ran)

  END SUBROUTINE randomNumber_single_r4





  SUBROUTINE randomNumber_single_r8(ran)

    IMPLICIT NONE
    REAL(rprec8), INTENT(out) :: ran

    IF (first_ran) THEN
       idum = idum_prm
       first_ran = .FALSE.
    END IF
    CALL ran_pmm(idum, ran)

  END SUBROUTINE randomNumber_single_r8





  SUBROUTINE randomNumber_single_r16(ran)

    IMPLICIT NONE
    REAL(rprec16), INTENT(out) :: ran

    IF (first_ran) THEN
       idum = idum_prm
       first_ran = .FALSE.
    END IF
    CALL ran_pmm(idum, ran)

  END SUBROUTINE randomNumber_single_r16





  SUBROUTINE randomNumber_array_r4(ran)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:), INTENT(out) :: ran
    INTEGER :: i

    DO i=1,SIZE(ran)
       CALL randomNumber(ran(i))
    END DO

  END SUBROUTINE randomNumber_array_r4





  SUBROUTINE randomNumber_2array_r4(ran)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:,:), INTENT(out) :: ran
    INTEGER :: i, j

    DO i=1,SIZE(ran,dim=1)
       DO j=1,SIZE(ran,dim=2)
          CALL randomNumber(ran(i,j))
       END DO
    END DO

  END SUBROUTINE randomNumber_2array_r4





  SUBROUTINE randomNumber_array_r8(ran)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(out) :: ran
    INTEGER :: i

    DO i=1,SIZE(ran)
       CALL randomNumber(ran(i))
    END DO

  END SUBROUTINE randomNumber_array_r8





  SUBROUTINE randomNumber_2array_r8(ran)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(out) :: ran
    INTEGER :: i, j

    DO i=1,SIZE(ran,dim=1)
       DO j=1,SIZE(ran,dim=2)
          CALL randomNumber(ran(i,j))
       END DO
    END DO

  END SUBROUTINE randomNumber_2array_r8





  SUBROUTINE randomNumber_array_r16(ran)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(out) :: ran
    INTEGER :: i

    DO i=1,SIZE(ran)
       CALL randomNumber(ran(i))
    END DO

  END SUBROUTINE randomNumber_array_r16





  SUBROUTINE randomGaussian_single_r4(rangauss)

    IMPLICIT NONE
    INTEGER, PARAMETER :: nran = 12
    REAL(rprec4), INTENT(out) :: rangauss
    REAL(rprec4), DIMENSION(nran) :: ran
    REAL(rprec4) :: x1, sum 
    REAL(rprec4), SAVE :: x2, sln
    LOGICAL, SAVE :: second = .FALSE.

    IF (second) THEN
       ! Use the second random number generated on last call:
       rangauss = x2*sln
       second = .FALSE.
    ELSE
       ! Generate a pair of random normals:
       second = .TRUE.
       sum = HUGE(sum)
       DO WHILE (sum >= 1.0_rprec4)
          CALL randomNumber(ran)
          x1 = 2.0_rprec4 * ran(1) - 1.0_rprec4
          ! Because congruential algorithms do not strictly speaking
          ! produce independent random numbers, more randomness has
          ! been created by not taking the two first numbers, but the
          ! first and one of the numbers coming after the two first
          ! numbers.
          x2 = 2.0_rprec4 * ran(2+CEILING((nran-2)*ran(2))) - 1.0_rprec4
          ! Tiny added to prevent LOG(zero)/zero:
          sum = x1**2.0_rprec4 + x2**2.0_rprec4 + TINY(sum)
       END DO
       sln = SQRT((-2.0_rprec4*LOG(sum))/sum)
       rangauss = x1 * sln
    END IF

  END SUBROUTINE randomGaussian_single_r4





  SUBROUTINE randomGaussian_single_r8(rangauss)

    IMPLICIT NONE
    INTEGER, PARAMETER :: nran = 12
    REAL(rprec8), INTENT(out) :: rangauss
    REAL(rprec8), DIMENSION(nran) :: ran
    REAL(rprec8) :: x1, sum 
    REAL(rprec8), SAVE :: x2, sln
    LOGICAL, SAVE :: second = .FALSE.

    IF (second) THEN
       ! Use the second random number generated on last call:
       rangauss = x2*sln
       second = .FALSE.
    ELSE
       ! Generate a pair of random normals:
       second = .TRUE.
       sum = HUGE(sum)
       DO WHILE (sum >= 1.0_rprec8)
          CALL randomNumber(ran)
          x1 = 2.0_rprec8 * ran(1) - 1.0_rprec8
          ! Because congruential algorithms do not strictly speaking
          ! produce independent random numbers, more randomness has
          ! been created by not taking the two first numbers, but the
          ! first and one of the numbers coming after the two first
          ! numbers.
          x2 = 2.0_rprec8 * ran(2+CEILING((nran-2)*ran(2))) - 1.0_rprec8
          ! Tiny added to prevent LOG(zero)/zero:
          sum = x1**2.0_rprec8 + x2**2.0_rprec8 + TINY(sum)
       END DO
       sln = SQRT((-2.0_rprec8*LOG(sum))/sum)
       rangauss = x1 * sln
    END IF

  END SUBROUTINE randomGaussian_single_r8





  SUBROUTINE randomGaussian_single_r16(rangauss)

    IMPLICIT NONE
    INTEGER, PARAMETER :: nran = 12
    REAL(rprec16), INTENT(out) :: rangauss
    REAL(rprec16), DIMENSION(nran) :: ran
    REAL(rprec16) :: x1, sum 
    REAL(rprec16), SAVE :: x2, sln
    LOGICAL, SAVE :: second = .FALSE.

    IF (second) THEN
       ! Use the second random number generated on last call:
       rangauss = x2*sln
       second = .FALSE.
    ELSE
       ! Generate a pair of random normals:
       second = .TRUE.
       sum = HUGE(sum)
       DO WHILE (sum >= 1.0_rprec16)
          CALL randomNumber(ran)
          x1 = 2.0_rprec16 * ran(1) - 1.0_rprec16
          ! Because congruential algorithms do not strictly speaking
          ! produce independent random numbers, more randomness has
          ! been created by not taking the two first numbers, but the
          ! first and one of the numbers coming after the two first
          ! numbers.
          x2 = 2.0_rprec16 * ran(2+CEILING((nran-2)*ran(2))) - 1.0_rprec16
          ! Tiny added to prevent LOG(zero)/zero:
          sum = x1**2.0_rprec16 + x2**2.0_rprec16 + TINY(sum)
       END DO
       sln = SQRT((-2.0_rprec16*LOG(sum))/sum)
       rangauss = x1 * sln
    END IF

  END SUBROUTINE randomGaussian_single_r16





  SUBROUTINE randomGaussian_array_r4(rangauss)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:), INTENT(out) :: rangauss 
    INTEGER :: i

    DO i=1,SIZE(rangauss)
       CALL randomGaussian(rangauss(i))
    END DO

  END SUBROUTINE randomGaussian_array_r4





  SUBROUTINE randomGaussian_array_r8(rangauss)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(out) :: rangauss 
    INTEGER :: i

    DO i=1,SIZE(rangauss)
       CALL randomGaussian(rangauss(i))
    END DO

  END SUBROUTINE randomGaussian_array_r8





  SUBROUTINE randomGaussian_array_r16(rangauss)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(out) :: rangauss 
    INTEGER :: i

    DO i=1,SIZE(rangauss)
       CALL randomGaussian(rangauss(i))
    END DO

  END SUBROUTINE randomGaussian_array_r16





  SUBROUTINE setRandomSeed(seed)

    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(in) :: seed
    INTEGER, DIMENSION(8) :: dt 
    INTEGER, DIMENSION(1) :: seed_array
    INTEGER :: k

    CALL RANDOM_SEED(size=k)
    IF (PRESENT(seed)) THEN
       seed_array = seed
       CALL RANDOM_SEED(put=seed_array(1:k))
    ELSE
       CALL DATE_AND_TIME(values=dt)
       seed_array = dt(8)
       CALL RANDOM_SEED(put=seed_array(1:k))
    END IF

  END SUBROUTINE setRandomSeed





END MODULE random
