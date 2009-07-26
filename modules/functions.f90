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
!! *Module*description*:
!!
!! Contains different mathematical functions.
!!
!! @author  MG
!! @version 2008-08-12
!!
MODULE functions

  USE parameters
  IMPLICIT NONE

  INTERFACE stumpff
     MODULE PROCEDURE stumpff_r8_series
     MODULE PROCEDURE stumpff_r16_series
     MODULE PROCEDURE stumpff_r8_recursive
     MODULE PROCEDURE stumpff_r16_recursive
  END INTERFACE



CONTAINS




  !! *Description*:
  !!
  !! Determines Stumpff-function c_k(x) by using series. The number of
  !! terms is dictated by the precision: the series is truncated when
  !! the accuracy given by additional terms can not be presented with
  !! the used precision. 
  !!
  !! Input values: x=R and k=0,1,2,...
  !!
  !! If k<0, the function returns 0.
  !!
  REAL(rprec8) FUNCTION stumpff_r8_series(x, k)

    IMPLICIT NONE
    REAL(rprec8), INTENT(in) :: x
    INTEGER, INTENT(in) :: k
    REAL(rprec8) :: term, stumpff_old
    INTEGER :: i

    ! k = 0,1,2,...
    IF (k < 0) THEN
       stumpff_r8_series = 0.0_rprec8
       RETURN
    END IF

    ! Compute the first term of the Stumpff-"function":
    term = 1.0_rprec8/factorial(k)
    stumpff_r8_series = term

    ! Compute terms of the Stumpff-"function" until the 
    ! limit of the computer accuracy is reached and 
    ! further improvement of the accuracy is pointless:  
    stumpff_old = HUGE(stumpff_old)
    i = 0
    DO
       i = i + 1
       term = term * x/REAL((k+i*2)*(k+i*2-1),rprec8)
       stumpff_r8_series = stumpff_r8_series + &
            REAL((-1)**i,rprec8) * term
       IF (ABS(stumpff_r8_series - stumpff_old) < &
            10.0_rprec8*EPSILON(stumpff_r8_series)) THEN
          EXIT
       END IF
       stumpff_old = stumpff_r8_series
    END DO

  END FUNCTION stumpff_r8_series





  !! *Description*:
  !!
  !! Determines Stumpff-function c_k(x) by using series. The number of
  !! terms is dictated by the precision: the series is truncated when
  !! the accuracy given by additional terms can not be presented with
  !! the used precision. 
  !!
  !! Input values: x=R and k=0,1,2,...
  !!
  !! If k<0, the function returns 0.
  !!
  REAL(rprec16) FUNCTION stumpff_r16_series(x, k)

    IMPLICIT NONE
    REAL(rprec16), INTENT(in) :: x
    INTEGER, INTENT(in) :: k
    REAL(rprec16) :: term, stumpff_old
    INTEGER :: i

    ! k = 0,1,2,...
    IF (k < 0) THEN
       stumpff_r16_series = 0.0_rprec16
       RETURN
    END IF

    ! Compute the first term of the Stumpff-"function":
    term = 1.0_rprec16/factorial(k)
    stumpff_r16_series = term

    ! Compute terms of the Stumpff-"function" until the 
    ! limit of the computer accuracy is reached and 
    ! further improvement of the accuracy is pointless:  
    stumpff_old = HUGE(stumpff_old)
    i = 0
    DO
       i = i + 1
       term = term * x/REAL((k+i*2)*(k+i*2-1),rprec16)
       stumpff_r16_series = stumpff_r16_series + &
            REAL((-1)**i,rprec16) * term
       IF (ABS(stumpff_r16_series - stumpff_old) < &
            10.0_rprec16*EPSILON(stumpff_r16_series)) THEN
          EXIT
       END IF
       stumpff_old = stumpff_r16_series
    END DO

  END FUNCTION stumpff_r16_series





  !! *Description*:
  !!
  !! Determines Stumpff-function c_k(x) by using the recursion
  !! formula.
  !!
  !! If 'to_higher' is abscent or .false., the function returns c_k(x)
  !! assuming that 'stumpff' is c_(k+2)(x). If 'to_higher' is .true.,
  !! the function returns c_(k+2)(x) assuming that 'stumpff' is
  !! c_(k)(x).
  !!
  !! Input values: x=R, k=0,1,2,..., stumpff>0
  !!
  !! If k<0 or stumpff<=0, the function returns 0.
  !!
  REAL(rprec8) FUNCTION stumpff_r8_recursive(x, k, stumpff, to_higher)

    IMPLICIT NONE
    REAL(rprec8), INTENT(in) :: x, stumpff
    INTEGER, INTENT(in) :: k
    LOGICAL, OPTIONAL :: to_higher

    ! k = 0,1,2,... and stumpff > 0
    IF (k < 0 .OR. stumpff <= 0) THEN
       stumpff_r8_recursive = 0.0_rprec8
       RETURN
    END IF

    IF (PRESENT(to_higher)) THEN
       IF (to_higher) THEN
          ! c_k -> c_(k+2)
          stumpff_r8_recursive = (1.0_rprec8/factorial(k-2) - stumpff)/x
          RETURN
       END IF
    END IF

    ! c_(k+2) -> c_k
    stumpff_r8_recursive = 1.0_rprec8/factorial(k) - x*stumpff

  END FUNCTION stumpff_r8_recursive





  !! *Description*:
  !!
  !! Determines Stumpff-function c_k(x) by using the recursion
  !! formula.
  !!
  !! If 'to_higher' is abscent or .false., the function returns c_k(x)
  !! assuming that 'stumpff' is c_(k+2)(x). If 'to_higher' is .true.,
  !! the function returns c_(k+2)(x) assuming that 'stumpff' is
  !! c_(k)(x).
  !!
  !! Input values: x=R, k=0,1,2,..., stumpff>0
  !!
  !! If k<0 or stumpff<=0, the function returns 0.
  !!
  REAL(rprec16) FUNCTION stumpff_r16_recursive(x, k, stumpff, to_higher)

    IMPLICIT NONE
    REAL(rprec16), INTENT(in) :: x, stumpff
    INTEGER, INTENT(in) :: k
    LOGICAL, OPTIONAL :: to_higher

    ! k = 0,1,2,... and stumpff > 0
    IF (k < 0 .OR. stumpff <= 0) THEN
       stumpff_r16_recursive = 0.0_rprec16
       RETURN
    END IF

    IF (PRESENT(to_higher)) THEN
       IF (to_higher) THEN
          ! c_(k+2) <- c_k
          stumpff_r16_recursive = (1.0_rprec16/factorial(k-2) - stumpff)/x
          RETURN
       END IF
    END IF

    ! c_k <- c_(k+2)
    stumpff_r16_recursive = 1.0_rprec16/factorial(k) - x*stumpff

  END FUNCTION stumpff_r16_recursive





  INTEGER FUNCTION factorial(k)

    INTEGER, INTENT(in) :: k
    INTEGER             :: l

    ! k = 0,1,2,...
    IF (k < 0) THEN
       factorial = 0
       RETURN
    END IF

    l = k
    factorial = 1
    DO WHILE (l > 0)
       factorial = factorial * l
       l = l - 1
    END DO

  END FUNCTION factorial





END MODULE functions
