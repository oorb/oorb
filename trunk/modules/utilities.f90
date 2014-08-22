!====================================================================!
!                                                                    !
! Copyright 2002-2013,2014                                           !
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
!! Independent utilities.
!!
!! @author  MG
!! @version 2014-08-22
!!
MODULE utilities

  USE parameters

  IMPLICIT NONE

  PRIVATE :: integerToArray_ri8
  PRIVATE :: arrayToInteger_r8i8
  PRIVATE :: arrayToReal_r8r16
  PRIVATE :: cumulativeDistribution_2d_r8
  PRIVATE :: realToArray_r16r8
  PRIVATE :: reallocate_ch_1
  PRIVATE :: reallocate_ch_2
  PRIVATE :: reallocate_r4_1
  PRIVATE :: reallocate_r4_2
  PRIVATE :: reallocate_r8_1
  PRIVATE :: reallocate_r8_2
  PRIVATE :: reallocate_r8_3
  PRIVATE :: reallocate_r16_1
  PRIVATE :: reallocate_r16_2
  PRIVATE :: reallocate_i1_1
  PRIVATE :: reallocate_i4_1
  PRIVATE :: reallocate_i4_2
  PRIVATE :: reallocate_i8_1
  PRIVATE :: reallocate_i8_2
  PRIVATE :: reallocate_l_1
  PRIVATE :: reallocate_l_2
  PRIVATE :: secToHMS_r4
  PRIVATE :: secToHMS_r8
  PRIVATE :: toInt_4
  PRIVATE :: toReal_8
  PRIVATE :: toReal_16
  PRIVATE :: toString_i4
  PRIVATE :: toString_i8
  PRIVATE :: toString_r4
  PRIVATE :: toString_r8

  INTERFACE arrayToInteger
     MODULE PROCEDURE arrayToInteger_r8i8
  END INTERFACE arrayToInteger

  INTERFACE arrayToReal
     MODULE PROCEDURE arrayToReal_r8r16
  END INTERFACE arrayToReal

  INTERFACE cumulativeDistribution
     MODULE PROCEDURE cumulativeDistribution_2d_r8
  END INTERFACE cumulativeDistribution

  INTERFACE integerToArray
     MODULE PROCEDURE integerToArray_ri8
  END INTERFACE integerToArray

  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_i4
     MODULE PROCEDURE imaxloc_i8
     MODULE PROCEDURE imaxloc_r4
     MODULE PROCEDURE imaxloc_r8
     MODULE PROCEDURE imaxloc_r16
  END INTERFACE imaxloc

  INTERFACE iminloc
     MODULE PROCEDURE iminloc_i4
     MODULE PROCEDURE iminloc_i8
     MODULE PROCEDURE iminloc_r4
     MODULE PROCEDURE iminloc_r8
     MODULE PROCEDURE iminloc_r16
  END INTERFACE iminloc

  INTERFACE realToArray
     MODULE PROCEDURE realToArray_r16r8
  END INTERFACE realToArray

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_ch_1
     MODULE PROCEDURE reallocate_ch_2
     MODULE PROCEDURE reallocate_r4_1
     MODULE PROCEDURE reallocate_r4_2
     MODULE PROCEDURE reallocate_r8_1
     MODULE PROCEDURE reallocate_r8_2
     MODULE PROCEDURE reallocate_r8_3
     MODULE PROCEDURE reallocate_r16_1
     MODULE PROCEDURE reallocate_r16_2
     MODULE PROCEDURE reallocate_i1_1
     MODULE PROCEDURE reallocate_i4_1
     MODULE PROCEDURE reallocate_i4_2
     MODULE PROCEDURE reallocate_i8_1
     MODULE PROCEDURE reallocate_i8_2
     MODULE PROCEDURE reallocate_l_1
     MODULE PROCEDURE reallocate_l_2
  END INTERFACE reallocate

  INTERFACE swap
     MODULE PROCEDURE swap_i4
     MODULE PROCEDURE swap_i8
     MODULE PROCEDURE swap_scalar_r8
     MODULE PROCEDURE swap_vector_r8
     MODULE PROCEDURE swap_r16
     MODULE PROCEDURE swap_ch
  END INTERFACE swap

  INTERFACE toInt
     MODULE PROCEDURE toInt_4
  END INTERFACE toInt

  INTERFACE secToHMS
     MODULE PROCEDURE secToHMS_r4
     MODULE PROCEDURE secToHMS_r8
  END INTERFACE secToHMS

  INTERFACE toReal
     MODULE PROCEDURE toReal_8
     MODULE PROCEDURE toReal_16
  END INTERFACE toReal

  INTERFACE toString
     MODULE PROCEDURE toString_i4
     MODULE PROCEDURE toString_i8
     MODULE PROCEDURE toString_r4
     MODULE PROCEDURE toString_r8
  END INTERFACE toString



CONTAINS





  INTEGER(iprec8) FUNCTION arrayToInteger_r8i8(array, elements, bin_size, nbins, bounds, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)     :: array    ! input array
    LOGICAL, DIMENSION(:), INTENT(in)         :: elements ! used elements
    REAL(rprec8), DIMENSION(:), INTENT(in)     :: bin_size ! bin sizes
    INTEGER, DIMENSION(:), INTENT(in)         :: nbins     ! number of bins
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: bounds   ! variable bounds
    LOGICAL, INTENT(out)                      :: error
    INTEGER(iprec8), DIMENSION(:), ALLOCATABLE :: kp
    INTEGER, DIMENSION(:), ALLOCATABLE        :: help, nbins_
    INTEGER(iprec8)                            :: bins_coeff, intgr
    INTEGER                                   :: ncolumn, i, j, err

    ! Make an index (help) of the elements that are used:
    ncolumn = COUNT(elements)
    ALLOCATE(kp(ncolumn), help(ncolumn), nbins_(ncolumn), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(help, stat=err)
       DEALLOCATE(nbins_, stat=err)
       RETURN
    END IF
    j = 0
    DO i=1,SIZE(elements)
       IF (elements(i)) THEN
          j = j + 1
          help(j) = i
          nbins_(j) = nbins(i)
       END IF
    END DO

    ! Compute coordinates of the bin:
    DO i=1,SIZE(kp)
       ! Compute a value describing a single elements position:
       kp(i) = INT((array(help(i)) - bounds(help(i),1)) / &
            bin_size(help(i)),4) + 1
    END DO

    ! If one or more of the elements are out of bounds, skip the
    ! array:
    IF (ANY(kp>nbins_) .OR. ANY(kp<0)) THEN
       arrayToInteger_r8i8 = 0_iprec8
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(help, stat=err)
       DEALLOCATE(nbins_, stat=err)
       RETURN
    END IF

    ! Transform bin coordinates to a single integer value
    ! (size(kp) == size(help)):
    intgr = 0_iprec8
    DO i=1,SIZE(kp)
       bins_coeff = 1_iprec8
       DO j=2,i
          bins_coeff = bins_coeff*INT(nbins_(j-1),iprec8)
       END DO
       intgr = intgr + INT((kp(i) - 1),iprec8)*bins_coeff
    END DO
    intgr = intgr + 1_iprec8

    DEALLOCATE(kp, help, nbins_, stat=err)
    IF (err /= 0) THEN
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(help, stat=err)
       DEALLOCATE(nbins_, stat=err)
    END IF
    arrayToInteger_r8i8 = intgr

  END FUNCTION arrayToInteger_r8i8





  REAL(rprec16) FUNCTION arrayToReal_r8r16(array, elements, bin_size, nbins, bounds, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)   :: array    ! input array
    LOGICAL, DIMENSION(:), INTENT(in)       :: elements ! used elements
    REAL(rprec8), DIMENSION(:), INTENT(in)   :: bin_size ! bin sizes
    INTEGER, DIMENSION(:), INTENT(in)       :: nbins     ! number of bins
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: bounds   ! variable bounds
    LOGICAL, INTENT(out)                    :: error
    INTEGER(iprec8), DIMENSION(:), ALLOCATABLE   :: kp
    INTEGER, DIMENSION(:), ALLOCATABLE      :: help, nbins_
    REAL(rprec16)                            :: bins_coeff, rreal
    REAL(rprec8)                             :: tmp
    INTEGER                                 :: ncolumn, i, j, err

    ! Make an index (help) of the elements that are used:
    ncolumn = COUNT(elements)
    ALLOCATE(kp(ncolumn), help(ncolumn), nbins_(ncolumn), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(help, stat=err)
       DEALLOCATE(nbins_, stat=err)
       RETURN
    END IF
    j = 0
    DO i=1,SIZE(elements,dim=1)
       IF (elements(i)) THEN
          j = j + 1
          help(j) = i
          nbins_(j) = nbins(i)
       END IF
    END DO

    ! Compute coordinates of the bin:
    DO i=1,SIZE(kp,dim=1)
       ! Compute a value describing a single elements position:
       tmp = array(help(i)) - bounds(help(i),1)
       tmp = tmp / bin_size(help(i))
       kp(i) = INT(tmp,4) + 1
    END DO

    ! If one or more of the elements are out of bounds, skip the
    ! array:
    IF (ANY(kp > nbins_) .OR. ANY(kp < 0)) THEN
       arrayToReal_r8r16 = 0.0_rprec16
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(help, stat=err)
       DEALLOCATE(nbins_, stat=err)
       RETURN
    END IF

    ! Transform bin coordinates to a single integer value
    ! (size(kp) == size(help)):
    rreal = 0.0_rprec16
    DO i=1,SIZE(kp,dim=1)
       bins_coeff = 1.0_rprec16
       DO j=2,i
          bins_coeff = bins_coeff*REAL(nbins_(j-1),rprec16)
       END DO
       rreal = rreal + REAL((kp(i) - 1),rprec16)*bins_coeff
    END DO
    rreal = rreal + 1.0_rprec16

    DEALLOCATE(kp, help, nbins_, stat=err)
    IF (err /= 0) THEN
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(help, stat=err)
       DEALLOCATE(nbins_, stat=err)
    END IF
    arrayToReal_r8r16 = rreal

  END FUNCTION arrayToReal_r8r16





  INTEGER FUNCTION imaxloc_i4(array)

    IMPLICIT NONE
    INTEGER(iprec4), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imax

    imax = MAXLOC(array)
    imaxloc_i4 = imax(1)

  END FUNCTION imaxloc_i4





  INTEGER FUNCTION imaxloc_i8(array)

    IMPLICIT NONE
    INTEGER(iprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imax

    imax = MAXLOC(array)
    imaxloc_i8 = imax(1)

  END FUNCTION imaxloc_i8





  INTEGER FUNCTION imaxloc_r4(array)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imax

    imax = MAXLOC(array)
    imaxloc_r4 = imax(1)

  END FUNCTION imaxloc_r4





  INTEGER FUNCTION imaxloc_r8(array)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imax

    imax = MAXLOC(array)
    imaxloc_r8 = imax(1)

  END FUNCTION imaxloc_r8





  INTEGER FUNCTION imaxloc_r16(array)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imax

    imax = MAXLOC(array)
    imaxloc_r16 = imax(1)

  END FUNCTION imaxloc_r16





  INTEGER FUNCTION iminloc_i4(array)

    IMPLICIT NONE
    INTEGER(iprec4), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imin

    imin = MINLOC(array)
    iminloc_i4 = imin(1)

  END FUNCTION iminloc_i4





  INTEGER FUNCTION iminloc_i8(array)

    IMPLICIT NONE
    INTEGER(iprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imin

    imin = MINLOC(array)
    iminloc_i8 = imin(1)

  END FUNCTION iminloc_i8





  INTEGER FUNCTION iminloc_r4(array)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imin

    imin = MINLOC(array)
    iminloc_r4 = imin(1)

  END FUNCTION iminloc_r4





  INTEGER FUNCTION iminloc_r8(array)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imin

    imin = MINLOC(array)
    iminloc_r8 = imin(1)

  END FUNCTION iminloc_r8





  INTEGER FUNCTION iminloc_r16(array)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(1) :: imin

    imin = MINLOC(array)
    iminloc_r16 = imin(1)

  END FUNCTION iminloc_r16





  FUNCTION integerToArray_ri8(intgr, elements, bin_size, nbins, bounds, error)

    IMPLICIT NONE
    INTEGER(iprec8), INTENT(in)              :: intgr    ! input integer
    LOGICAL, DIMENSION(:), INTENT(in)        :: elements ! used elements
    REAL(rprec8), DIMENSION(:), INTENT(in)   :: bin_size ! bin sizes
    INTEGER, DIMENSION(:), INTENT(in)        :: nbins     ! number of bins
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: bounds   ! variable bounds
    LOGICAL, INTENT(out)                     :: error
    REAL(rprec8), DIMENSION(SIZE(elements))  :: integerToArray_ri8

    INTEGER(iprec8)                    :: intgr_coeff, bins_coeff
    INTEGER, DIMENSION(:), ALLOCATABLE :: kp, nbins_
    INTEGER                            :: ncolumn, i, j, k, err

    ! Make an index (help) of the elements that are used:
    ncolumn = COUNT(elements)
    ALLOCATE(kp(ncolumn), nbins_(ncolumn), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(nbins_, stat=err)
       RETURN
    END IF
    j = 0
    DO i=1,SIZE(elements)
       IF (elements(i)) THEN
          j = j + 1
          nbins_(j) = nbins(i)
       END IF
    END DO

    ! Transform integer to bin coordinates
    ! (size(kp) == size(help)):
    bins_coeff = 1_rprec8
    DO i=ncolumn, 1, -1
       intgr_coeff = intgr
       DO j=i+1, ncolumn
          bins_coeff = 1_rprec8
          DO k=1,j-1
             bins_coeff = bins_coeff*INT(nbins_(k),rprec8)
          END DO
          intgr_coeff = intgr_coeff - INT((kp(j)-1)*bins_coeff,rprec8)
       END DO
       bins_coeff = 1_rprec8
       DO j=1,i-1
          bins_coeff = bins_coeff*INT(nbins_(j),rprec8)
       END DO
       kp(i) = 1 + INT(intgr_coeff/REAL(bins_coeff),iprec4)
    END DO
    kp(1) = kp(1) - 1

    ! Transform bin coordinates to real values: 
    j = 0
    DO i=1,SIZE(elements)
       IF (elements(i)) THEN
          j = j + 1
          integerToArray_ri8(i) = bounds(i,1) + kp(j)*bin_size(i)
       ELSE
          integerToArray_ri8(i) = -1.0_rprec8
       END IF
    END DO

    DEALLOCATE(kp, nbins_, stat=err)
    IF (err /= 0) THEN
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(nbins_, stat=err)
    END IF

  END FUNCTION integerToArray_ri8





  LOGICAL FUNCTION is_downcase(d)

    ! Description:
    !
    ! This function returns TRUE if the character 'd' is a lowercase letter.
    ! It returns FALSE otherwise.
    !
    ! Programming interface:
    !
    ! Name         In/out  Type          Structure  Meaning
    !
    ! is_downcase  out     logical       scalar     TRUE if d is lowercase letter 
    !                                               FALSE otherwise
    ! d            in      character(1)  scalar     The character to be examined
    !
    ! Errors:
    !
    ! No errors should occur and none are trapped.
    ! This routine assumes the ASCII character set.

    IMPLICIT NONE
    CHARACTER(len=1) :: d

    is_downcase = ((d >= 'a') .AND. (d <= 'z'))

  END FUNCTION is_downcase





  !! Description:
  !!
  !! This subroutine takes a 2D distribution and computes its
  !! cumulative version as a function of rows and columns.
  !!
  !! Programming interface:
  !!
  !! Name         In/out  Type          Structure  Meaning
  !!
  !! in           in      real(8)       array      Incremental 2D distribution
  !! out          out     real(8)       array      Cumulative 2D distribution
  !!
  !! Errors:
  !!
  !! No errors should occur and none are trapped.
  !!
  SUBROUTINE cumulativeDistribution_2d_r8(in, out)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: in
    REAL(rprec8), DIMENSION(:,:), INTENT(out) :: out

    INTEGER :: i, j

    DO i=1,SIZE(in,dim=1)
       DO j=1,SIZE(in,dim=2)
          IF (i == 1) THEN
             out(i,j) = SUM(in(i,:j))
          ELSE
             out(i,j) = SUM(in(i,:j)) + out(i-1,j)
          END IF
       END DO
    END DO

  END SUBROUTINE cumulativeDistribution_2d_r8





  FUNCTION outerand(a,b)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL, DIMENSION(SIZE(a),SIZE(b)) :: outerand

    outerand = SPREAD(a,dim=2,ncopies=SIZE(b)) .AND. &
         SPREAD(b,dim=1,ncopies=SIZE(a))

  END FUNCTION outerand





  FUNCTION realToArray_r16r8(rreal, elements, bin_size, nbins, bounds, error)

    IMPLICIT NONE
    REAL(rprec16), INTENT(in)               :: rreal    ! input integer (real)
    LOGICAL, DIMENSION(:), INTENT(in)      :: elements ! used elements
    REAL(rprec8), DIMENSION(:), INTENT(in)   :: bin_size ! bin sizes
    INTEGER, DIMENSION(:), INTENT(in)      :: nbins     ! number of bins
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: bounds   ! variable bounds
    LOGICAL, INTENT(out)                   :: error
    REAL(rprec8), DIMENSION(:), POINTER      :: realToArray_r16r8
    REAL(rprec16)                          :: rreal_coeff, bins_coeff
    INTEGER, DIMENSION(:), ALLOCATABLE     :: kp, nbins_
    INTEGER                                :: ncolumn, i, j, k, err

    ! Make an index (help) of the elements that are used:
    ncolumn = COUNT(elements)
    ALLOCATE(realToArray_r16r8(SIZE(elements)), kp(ncolumn), &
         nbins_(ncolumn), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(realToArray_r16r8, stat=err)
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(nbins_, stat=err)
       RETURN
    END IF
    j = 0
    DO i=1,SIZE(elements)
       IF (elements(i)) THEN
          j = j + 1
          nbins_(j) = nbins(i)
       END IF
    END DO

    ! Transform integer to bin coordinates
    ! (size(kp) == size(help)):
    bins_coeff = 1.0_rprec16
    DO i=ncolumn, 1, -1
       rreal_coeff = rreal
       DO j=i+1, ncolumn
          bins_coeff = 1.0_rprec16
          DO k=1,j-1
             bins_coeff = bins_coeff*REAL(nbins_(k),rprec16)
          END DO
          rreal_coeff = rreal_coeff - REAL((kp(j)-1)*bins_coeff,rprec16)
       END DO
       bins_coeff = 1.0_rprec16
       DO j=1,i-1
          bins_coeff = bins_coeff*REAL(nbins_(j),rprec16)
       END DO
       kp(i) = 1 + INT(rreal_coeff/bins_coeff,iprec4)
    END DO
    kp(1) = kp(1) - 1

    ! Transform bin coordinates to real values: 
    j = 0
    DO i=1,SIZE(elements)
       IF (elements(i)) THEN
          j = j + 1
          realToArray_r16r8(i) = bounds(i,1) + kp(j)*bin_size(i)
       ELSE
          realToArray_r16r8(i) = -1.0_rprec8
       END IF
    END DO

    DEALLOCATE(kp, nbins_, stat=err)
    IF (err /= 0) THEN
       DEALLOCATE(kp, stat=err)
       DEALLOCATE(nbins_, stat=err)
    END IF

  END FUNCTION realToArray_r16r8





  FUNCTION reallocate_ch_1(array,n)

    IMPLICIT NONE
    CHARACTER(len=*), DIMENSION(:), POINTER :: array
    CHARACTER(len=LEN(array(1))), DIMENSION(:), POINTER :: reallocate_ch_1
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_ch_1(n))
    reallocate_ch_1 = " "
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array, dim=1)
       reallocate_ch_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_ch_1





  FUNCTION reallocate_ch_2(array,n,m)

    IMPLICIT NONE
    CHARACTER(len=*), DIMENSION(:,:), POINTER :: array
    CHARACTER(len=LEN(array(1,1))), DIMENSION(:,:), POINTER :: reallocate_ch_2
    INTEGER, INTENT(in) :: n, m
    INTEGER :: nold, mold, i, j

    ALLOCATE(reallocate_ch_2(n,m))
    reallocate_ch_2 = " "
    IF (.NOT. ASSOCIATED(array)) THEN
       RETURN
    END IF
    nold = SIZE(array, dim=1)
    mold = SIZE(array, dim=2)
    DO i=1,MIN(n,nold)
       DO j=1,MIN(m,mold)
          reallocate_ch_2(i,j) = array(i,j)
       END DO
    END DO
    DEALLOCATE(array)

  END FUNCTION reallocate_ch_2





  FUNCTION reallocate_ch18_1(array,n)

    IMPLICIT NONE
    CHARACTER(len=18), DIMENSION(:), POINTER :: reallocate_ch18_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_ch18_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array, dim=1)
    reallocate_ch18_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_ch18_1





  FUNCTION reallocate_ch24_1(array,n)

    IMPLICIT NONE
    CHARACTER(len=24), DIMENSION(:), POINTER :: reallocate_ch24_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_ch24_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array, dim=1)
    reallocate_ch24_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_ch24_1





  FUNCTION reallocate_i1_1(array,n)

    IMPLICIT NONE
    INTEGER(iprec1), DIMENSION(:), POINTER :: reallocate_i1_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_i1_1(n))
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array, dim=1)
       reallocate_i1_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_i1_1





  FUNCTION reallocate_i4_1(array,n)

    IMPLICIT NONE
    INTEGER(iprec4), DIMENSION(:), POINTER :: reallocate_i4_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_i4_1(n))
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array, dim=1)
       reallocate_i4_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_i4_1





  FUNCTION reallocate_i4_2(array,n,m)

    IMPLICIT NONE
    INTEGER(iprec4), DIMENSION(:,:), POINTER :: reallocate_i4_2, array
    INTEGER, INTENT(in) :: n, m
    INTEGER :: nold, mold

    ALLOCATE(reallocate_i4_2(n,m))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,1)
    mold = SIZE(array,2)
    reallocate_i4_2(1:MIN(n,nold),1:MIN(m,mold)) = &
         array(1:MIN(n,nold),1:MIN(m,mold))
    DEALLOCATE(array)

  END FUNCTION reallocate_i4_2





  FUNCTION reallocate_i8_1(array,n)

    IMPLICIT NONE
    INTEGER(rprec8), DIMENSION(:), POINTER :: reallocate_i8_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_i8_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    reallocate_i8_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_i8_1





  FUNCTION reallocate_i8_2(array,n,m)

    IMPLICIT NONE
    INTEGER(rprec8), DIMENSION(:,:), POINTER :: reallocate_i8_2, array
    INTEGER, INTENT(in) :: n, m
    INTEGER :: nold, mold

    ALLOCATE(reallocate_i8_2(n,m))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,1)
    mold = SIZE(array,2)
    reallocate_i8_2(1:MIN(n,nold),1:MIN(m,mold)) = &
         array(1:MIN(n,nold),1:MIN(m,mold))
    DEALLOCATE(array)

  END FUNCTION reallocate_i8_2





  FUNCTION reallocate_r4_1(array,n)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:), POINTER :: reallocate_r4_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_r4_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    reallocate_r4_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_r4_1





  FUNCTION reallocate_r4_2(array,n,m)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:,:), POINTER :: reallocate_r4_2, array
    INTEGER, INTENT(in) :: n, m
    INTEGER :: nold, mold

    ALLOCATE(reallocate_r4_2(n,m))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,1)
    mold = SIZE(array,2)
    reallocate_r4_2(1:MIN(n,nold),1:MIN(m,mold)) = &
         array(1:MIN(n,nold),1:MIN(m,mold))
    DEALLOCATE(array)

  END FUNCTION reallocate_r4_2





  FUNCTION reallocate_r8_1(array,n)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), POINTER :: reallocate_r8_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_r8_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    reallocate_r8_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_r8_1





  FUNCTION reallocate_r8_2(array,n,m,nmask)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), POINTER :: reallocate_r8_2, array
    INTEGER, INTENT(in) :: n, m
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: nmask
    INTEGER :: i, in, nold, mold

    ALLOCATE(reallocate_r8_2(n,m))
    IF (ASSOCIATED(array) .AND. .NOT.PRESENT(nmask)) THEN
       nold = SIZE(array,1)
       mold = SIZE(array,2)
       reallocate_r8_2(1:MIN(n,nold),1:MIN(m,mold)) = &
            array(1:MIN(n,nold),1:MIN(m,mold))
       DEALLOCATE(array)
    ELSE IF (ASSOCIATED(array) .AND. PRESENT(nmask)) THEN
       in = 0
       DO i=1,nold
          IF (in == n) THEN
             EXIT
          END IF
          IF (nmask(i)) THEN
             in = in + 1
             reallocate_r8_2(in,1:MIN(m,mold)) = array(i,1:MIN(m,mold))
          END IF
       END DO
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_r8_2





  FUNCTION reallocate_r8_3(array,n,m,o)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:,:), POINTER :: reallocate_r8_3, array
    INTEGER, INTENT(in) :: n, m, o
    INTEGER :: nold, mold, oold

    ALLOCATE(reallocate_r8_3(n,m,o))
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array,1)
       mold = SIZE(array,2)
       oold = SIZE(array,3)
       reallocate_r8_3(1:MIN(n,nold),1:MIN(m,mold),1:MIN(o,oold)) = &
            array(1:MIN(n,nold),1:MIN(m,mold),1:MIN(o,oold))
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_r8_3





  FUNCTION reallocate_r16_1(array,n)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), POINTER :: reallocate_r16_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_r16_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    reallocate_r16_1(1:MIN(n,nold)) = array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_r16_1





  FUNCTION reallocate_r16_2(array,n,m)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:,:), POINTER :: reallocate_r16_2, array
    INTEGER, INTENT(in) :: n, m
    INTEGER :: nold, mold

    ALLOCATE(reallocate_r16_2(n,m))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,1)
    mold = SIZE(array,2)
    reallocate_r16_2(1:MIN(n,nold),1:MIN(m,mold)) = &
         array(1:MIN(n,nold),1:MIN(m,mold))
    DEALLOCATE(array)

  END FUNCTION reallocate_r16_2





  FUNCTION reallocate_l_1(array,n)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:), POINTER :: reallocate_l_1, array
    INTEGER, INTENT(in) :: n
    INTEGER :: nold

    ALLOCATE(reallocate_l_1(n))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    reallocate_l_1(1:MIN(n,nold)) = &
         array(1:MIN(n,nold))
    DEALLOCATE(array)

  END FUNCTION reallocate_l_1





  FUNCTION reallocate_l_2(array,n,m)

    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:), POINTER :: reallocate_l_2, array
    INTEGER, INTENT(in) :: n, m
    INTEGER :: nold, mold

    ALLOCATE(reallocate_l_2(n,m))
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,dim=1)
    mold = SIZE(array,dim=2)
    reallocate_l_2(1:MIN(n,nold),1:MIN(m,mold)) = &
         array(1:MIN(n,nold),1:MIN(m,mold))
    DEALLOCATE(array)

  END FUNCTION reallocate_l_2





  SUBROUTINE removeLeadingBlanks(str)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: str

    DO WHILE (INDEX(TRIM(str),CHAR(0)) == 1 .OR. &
         INDEX(TRIM(str),CHAR(32)) == 1)
       str(1:LEN_TRIM(str)) = str(2:LEN_TRIM(str)) // CHAR(32)
    END DO

  END SUBROUTINE removeLeadingBlanks





!!$  FUNCTION TRIM(str_in) RESULT(str_out)
!!$
!!$    IMPLICIT NONE
!!$    CHARACTER(len=*), INTENT(in) :: str_in
!!$    CHARACTER(len=LEN_TRIM(str_in)) :: str_out
!!$
!!$    str_out = str_in(1:LEN_TRIM(str_in))
!!$
!!$  END FUNCTION trim





  CHARACTER(len=64) FUNCTION secToHMS_r4(sec, error)

    IMPLICIT NONE
    REAL(rprec4), INTENT(in)    :: sec
    LOGICAL, INTENT(inout) :: error
    CHARACTER(len=16)      :: str
    INTEGER                :: h, min, ind

    secToHMS_r4 = ' '
    ! Hours:
    h = FLOOR(sec/3600.0_4)
    IF (h /= 0) THEN
       CALL toString(h, str, error)
       secToHMS_r4 = TRIM(secToHMS_r4) // TRIM(str) // 'h '
    END IF
    ! Minutes:
    min = FLOOR((sec - h*3600.0_4)/60.0_4)
    IF (min /= 0 .OR. h /= 0) THEN
       CALL toString(min, str, error)
       secToHMS_r4 = TRIM(secToHMS_r4) // TRIM(str) // 'min '
    END IF
    ! Seconds with 3 decimals:
    CALL toString((sec-h*3600.0_4-min*60.0_4), str, error)
    ind = INDEX(str,'.') + 6
    secToHMS_r4 = TRIM(secToHMS_r4) // str(1:ind) // 'sec'

  END FUNCTION secToHMS_r4





  CHARACTER(len=64) FUNCTION secToHMS_r8(sec, error)

    IMPLICIT NONE
    REAL(rprec8), INTENT(in)    :: sec
    LOGICAL, INTENT(inout) :: error
    CHARACTER(len=32)      :: str
    INTEGER                :: h, min, ind

    secToHMS_r8 = ' '
    ! Hours:
    h = FLOOR(sec/3600.0_rprec8)
    IF (h /= 0) THEN
       CALL toString(h, str, error)
       secToHMS_r8 = TRIM(secToHMS_r8) // TRIM(str) // 'h '
    END IF
    ! Minutes:
    min = FLOOR((sec - h*3600.0_rprec8)/60.0_rprec8)
    IF (min /= 0 .OR. h /= 0) THEN
       CALL toString(min, str, error)
       secToHMS_r8 = TRIM(secToHMS_r8) // TRIM(str) // 'min '
    END IF
    ! Seconds with 3 decimals:
    CALL toString((sec-h*3600.0_rprec8-min*60.0_rprec8), str, error)
    ind = INDEX(str,'.') + 6
    secToHMS_r8 = TRIM(secToHMS_r8) // str(1:ind) // 'sec'

  END FUNCTION secToHMS_r8





  SUBROUTINE swap_i4(first, second, mask)

    IMPLICIT NONE
    INTEGER(iprec4), INTENT(inout) :: first, second
    LOGICAL, OPTIONAL, INTENT(in) :: mask
    INTEGER(iprec4) :: dum

    IF (PRESENT(mask)) THEN
       IF (mask) THEN
          dum = first
          first = second
          second = dum
       END IF
    ELSE
       dum = first
       first = second
       second = dum
    END IF

  END SUBROUTINE swap_i4





  SUBROUTINE swap_i8(first, second, mask)

    IMPLICIT NONE
    INTEGER(iprec8), INTENT(inout) :: first, second
    LOGICAL, OPTIONAL, INTENT(in) :: mask
    INTEGER(iprec8) :: dum

    IF (PRESENT(mask)) THEN
       IF (mask) THEN
          dum = first
          first = second
          second = dum
       END IF
    ELSE
       dum = first
       first = second
       second = dum
    END IF

  END SUBROUTINE swap_i8





  SUBROUTINE swap_scalar_r8(first, second, mask)

    IMPLICIT NONE
    REAL(rprec8), INTENT(inout) :: first, second
    LOGICAL, OPTIONAL, INTENT(in) :: mask
    REAL(rprec8) :: dum

    IF (PRESENT(mask)) THEN
       IF (mask) THEN
          dum = first
          first = second
          second = dum
       END IF
    ELSE
       dum = first
       first = second
       second = dum
    END IF

  END SUBROUTINE swap_scalar_r8





  SUBROUTINE swap_vector_r8(first, second, mask)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: first, second
    LOGICAL, OPTIONAL, INTENT(in) :: mask
    REAL(rprec8), DIMENSION(SIZE(first,dim=1)) :: dum

    IF (PRESENT(mask)) THEN
       IF (mask) THEN
          dum = first
          first = second
          second = dum
       END IF
    ELSE
       dum = first
       first = second
       second = dum
    END IF

  END SUBROUTINE swap_vector_r8





  SUBROUTINE swap_r16(first, second, mask)

    IMPLICIT NONE
    REAL(rprec16), INTENT(inout) :: first, second
    LOGICAL, OPTIONAL, INTENT(in) :: mask
    REAL(rprec16) :: dum

    IF (PRESENT(mask)) THEN
       IF (mask) THEN
          dum = first
          first = second
          second = dum
       END IF
    ELSE
       dum = first
       first = second
       second = dum
    END IF

  END SUBROUTINE swap_r16





  SUBROUTINE swap_ch(first, second, mask)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: first, second
    LOGICAL, OPTIONAL, INTENT(in) :: mask
    CHARACTER(len=MAX(LEN(first),LEN(second))) :: dum

    IF (PRESENT(mask)) THEN
       IF (mask) THEN
          dum = first
          first = second
          second = dum
       END IF
    ELSE
       dum = first
       first = second
       second = dum
    END IF

  END SUBROUTINE swap_ch





  !! *Description*:
  !!
  !! Converts a string to a integer(4) number. An error is reported if
  !! the string cannot be converted to number.
  !!
  SUBROUTINE toInt_4(str, k, error)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: str
    INTEGER(iprec4), INTENT(out) :: k
    LOGICAL, INTENT(inout)       :: error
    INTEGER                      :: err

    READ(str,'(I30)',iostat=err) k
    IF (err /= 0) THEN
       error = .TRUE.
    END IF

  END SUBROUTINE toInt_4





  !! *Description*:
  !!
  !! Converts a string to a real(prec8) number. An error is reported if
  !! the string cannot be converted to number.
  !!
  SUBROUTINE toReal_8(str, r, error)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)   :: str
    REAL(rprec8), INTENT(out)        :: r
    LOGICAL, INTENT(inout)         :: error
    CHARACTER(len=LEN_TRIM(str)+2) :: lstr
    REAL(rprec8)                     :: plusminus
    INTEGER                        :: err

    plusminus = 1.0_rprec8

    lstr = TRIM(str)
    IF (INDEX(lstr, '+') == 1) THEN
       lstr = lstr(2:LEN_TRIM(lstr)) 
    ELSE IF (INDEX(lstr, '-') == 1) THEN
       plusminus = -1.0_rprec8
       lstr = lstr(2:LEN_TRIM(lstr)) 
    END IF
    IF (INDEX(str, '.') == 0) lstr = TRIM(str) // '.0'
    READ(lstr, '(F35.16)', iostat=err) r
    IF (err /= 0) THEN
       error = .TRUE.
    END IF
    r = plusminus * r

  END SUBROUTINE toReal_8





  !! *Description*:
  !!
  !! Converts a string to a real(16) number. An error is reported if
  !! the string cannot be converted to number.
  !!
  SUBROUTINE toReal_16(str, r, error)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)   :: str
    REAL(rprec16), INTENT(out)        :: r
    LOGICAL, INTENT(inout)         :: error
    CHARACTER(len=LEN_TRIM(str)+2) :: lstr
    REAL(rprec16)                     :: plusminus
    INTEGER                        :: err

    plusminus = 1.0_rprec16

    lstr = TRIM(str) 
    IF (INDEX(lstr, '+') == 1) THEN
       lstr = lstr(2:LEN_TRIM(lstr)) 
    ELSE IF (INDEX(lstr, '-') == 1) THEN
       plusminus = -1.0_rprec16
       lstr = lstr(2:LEN_TRIM(lstr)) 
    END IF
    IF (INDEX(str, '.') == 0) lstr = TRIM(str) // '.0' 
    READ(lstr, '(F64.32)', iostat=err) r
    IF (err /= 0) THEN
       error = .TRUE.
    END IF
    r = plusminus * r

  END SUBROUTINE toReal_16





  !! *Description*:
  !!
  !! This subroutine returns the input string with all the lowercase
  !! letters converted to uppercase letters.
  !!
  SUBROUTINE upcase(str)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: str
    INTEGER :: i, limit

    limit = LEN_TRIM(str)
    DO i=1,limit
       IF ((str(i:i) >= 'a') .AND. (str(i:i) <= 'z')) THEN
          str(i:i) = CHAR(ICHAR(str(i:i))-32)
       END IF
    END DO

  END SUBROUTINE upcase





  !! *Description*:
  !!
  !! This subroutine returns the input string with all the uppercase
  !! letters converted to lowercase letters.
  !!
  SUBROUTINE locase(str, error)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: str
    LOGICAL, INTENT(inout) :: error
    INTEGER :: i, limit, ilo, ihi

    limit = LEN(str)
    ilo = IACHAR('A')
    ihi = IACHAR('Z')
    DO i=1,limit
       IF (IACHAR(str(i:i)) > 127 .OR. IACHAR(str(i:i)) < 0) THEN
          error = .TRUE.
          RETURN
       END IF
       IF (IACHAR(str(i:i)) >= ilo .AND. IACHAR(str(i:i)) <= ihi) THEN
          str(i:i) = CHAR(ICHAR(str(i:i))+32)
       END IF
    END DO

  END SUBROUTINE locase





  !!
  !!
  SUBROUTINE toString_i4(k, str, error)

    ! Description:
    !
    ! This routine converts an integer number to a string.
    !
    ! Programming interface:
    !
    ! Name   In/Out  Type       Structure  Meaning
    ! 
    ! str    out    character  Scalar     The string
    ! k      in     integer(4) Scalar     The integer number to be converted
    ! error  inout  logical    Scalar     Set TRUE if an error occurs
    !
    ! Errors:
    !
    ! An error is reported if the number cannot be converted to an integer.

    IMPLICIT NONE
    INTEGER(iprec4), INTENT(in)   :: k
    CHARACTER(len=*), INTENT(out) :: str
    LOGICAL, INTENT(inout)        :: error
    INTEGER                       :: err

    WRITE(str,'(I0)',iostat=err) k
    IF (err /= 0) THEN
       error = .TRUE.
       str = '***** Error in utilities / toString *****'
    END IF

  END SUBROUTINE toString_i4





  !!
  !!
  SUBROUTINE toString_i8(k, str, error)

    ! Description:
    !
    ! This routine converts an integer number to a string.
    !
    ! Programming interface:
    !
    ! Name   In/Out  Type       Structure  Meaning
    ! 
    ! str    out    character  Scalar     The string
    ! k      in     integer(4) Scalar     The integer number to be converted
    ! error  inout  logical    Scalar     Set TRUE if an error occurs
    !
    ! Errors:
    !
    ! An error is reported if the number cannot be converted to an integer.

    IMPLICIT NONE
    INTEGER(iprec8), INTENT(in)   :: k
    CHARACTER(len=*), INTENT(out) :: str
    LOGICAL, INTENT(inout)        :: error
    INTEGER                       :: err

    WRITE(str,'(I0)',iostat=err) k
    IF (err /= 0) THEN
       error = .TRUE.
       str = '***** Error in utilities / toString *****'
    END IF

  END SUBROUTINE toString_i8





  !!
  !!
  SUBROUTINE toString_r4(r, str, error, frmt)

    ! Description:
    !
    ! This routine converts a real number to a string.
    !
    ! Programming interface:
    !
    ! Name   In/Out  Type       Structure  Meaning
    ! 
    ! str    out    character  Scalar     The string
    ! r4     in     real(4)    Scalar     The real number to be converted
    ! error  inout  logical    Scalar     Set TRUE if an error occurs
    !
    ! Errors:
    !
    ! An error message is written to the output string, if the number
    ! cannot be converted to an integer.

    IMPLICIT NONE
    REAL(rprec4), INTENT(in)               :: r
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frmt
    CHARACTER(len=*), INTENT(out)          :: str
    LOGICAL, INTENT(inout)                 :: error

    INTEGER :: err

    IF (PRESENT(frmt)) THEN
       WRITE(str,TRIM(frmt),iostat=err) r
    ELSE
       WRITE(str,'(F0.8)',iostat=err) r
    END IF
    IF (err /= 0) THEN
       error = .TRUE.
       str = '***** Error in utilities / toString *****'
    END IF
    CALL removeLeadingBlanks(str)

  END SUBROUTINE toString_r4





  !!
  !!
  SUBROUTINE toString_r8(r, str, error, frmt)

    ! Description:
    !
    ! This routine converts a real number to a string.
    !
    ! Programming interface:
    !
    ! Name   In/Out  Type       Structure  Meaning
    ! 
    ! str    out    character  Scalar     The string (len > 19)
    ! r4     in     real(8)    Scalar     The real number to be converted
    ! error  inout  logical    Scalar     Set TRUE if an error occurs
    !
    ! Errors:
    !
    ! An error message is written to the output string, if the number
    ! cannot be converted to an integer.

    IMPLICIT NONE
    REAL(rprec8), INTENT(in)               :: r
    CHARACTER(len=*), INTENT(out)          :: str
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frmt
    LOGICAL, INTENT(inout)                 :: error

    INTEGER :: err

    IF (PRESENT(frmt)) THEN
       WRITE(str,TRIM(frmt),iostat=err) r
    ELSE
       WRITE(str,'(F0.16)',iostat=err) r
    END IF
    IF (err /= 0) THEN
       error = .TRUE.
       str = '***** Error in utilities / toString *****'
    END IF
    str = ADJUSTL(str)

  END SUBROUTINE toString_r8





END MODULE utilities
