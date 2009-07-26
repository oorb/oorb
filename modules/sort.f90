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
!! Contains sorting and searching routines.
!!
!! @author  MG
!! @version 2008-06-05
!!
MODULE sort

  USE parameters
  USE utilities

  IMPLICIT NONE
  INTEGER, PARAMETER :: NN = 15, NSTACK = 100

  INTERFACE insertionSort
     MODULE PROCEDURE insertionSort_i8
     MODULE PROCEDURE insertionSort_r16
     MODULE PROCEDURE insertionSort_ch
  END INTERFACE

  INTERFACE quickSort
     MODULE PROCEDURE quickSort_i4
     MODULE PROCEDURE quickSort_i4_index
     MODULE PROCEDURE quickSort_i8
     MODULE PROCEDURE quickSort_i8_index
     MODULE PROCEDURE quickSort_r8
     MODULE PROCEDURE quickSort_r8_index
     MODULE PROCEDURE quickSort_r16
     MODULE PROCEDURE quickSort_ch
     MODULE PROCEDURE quickSort_ch_index
  END INTERFACE

  INTERFACE binarySearch
     MODULE PROCEDURE binarySearch_i8
     MODULE PROCEDURE binarySearch_r8
     MODULE PROCEDURE binarySearch_r16
     MODULE PROCEDURE binarySearch_ch_index
  END INTERFACE

  INTERFACE findLocation
     MODULE PROCEDURE findLocation_r8
     MODULE PROCEDURE findLocation_r8_indx
     MODULE PROCEDURE findLocation_r16
  END INTERFACE





CONTAINS





  !! *Description*:
  !!
  !! The insertion sort technique is the same as the one used for
  !! sorting cards: Pick out the second card and put it in order with
  !! the first card. Then pick out the third card and put it into the
  !! sequence among the first two. Continue until the last card has
  !! been picked out and inserted. In a random situation the algorithm
  !! scales as N**2, but assuming an almost sorted array it can even
  !! reach N.  Generally, the insertion sort should not be used
  !! for N >= 100. This implementation is an in-place sort.
  !!
  SUBROUTINE insertionSort_i8(array)

    IMPLICIT NONE
    INTEGER(iprec8), DIMENSION(:), INTENT(inout) :: array
    INTEGER(iprec8)                              :: element
    INTEGER                                      :: i, j

    DO i=2, SIZE(array)
       element = array(i)
       j = i
       DO WHILE (j > 1)
          IF (array(j-1) <= element) THEN
             EXIT
          END IF
          array(j) = array(j-1)
          j = j - 1
       END DO
       array(j) = element
    END DO

  END SUBROUTINE insertionSort_i8





  !! *Description*:
  !!
  !! The insertion sort technique is the same as the one used for
  !! sorting cards: Pick out the second card and put it in order with
  !! the first card. Then pick out the third card and put it into the
  !! sequence among the first two. Continue until the last card has
  !! been picked out and inserted. In a random situation the algorithm
  !! scales as N**2, but assuming an almost sorted array it can even
  !! reach N.  Generally, the insertion sort should not be used
  !! for N >= 100. This implementation is an in-place sort.
  !!
  SUBROUTINE insertionSort_ch(array)

    IMPLICIT NONE
    CHARACTER(len=*), DIMENSION(:), INTENT(inout) :: array
    CHARACTER(len=LEN(array(1)))                  :: element
    INTEGER                                       :: i, j

    DO i=2,SIZE(array)
       element = array(i)
       j = i
       DO WHILE (j > 1)
          IF (array(j-1) <= element) THEN
             EXIT
          END IF
          array(j) = array(j-1)
          j = j - 1
       END DO
       array(j) = element
    END DO

  END SUBROUTINE insertionSort_ch





  !! *Description*:
  !!
  !! The insertion sort technique is the same as the one used for
  !! sorting cards: Pick out the second card and put it in order with
  !! the first card. Then pick out the third card and put it into the
  !! sequence among the first two. Continue until the last card has
  !! been picked out and inserted. In a random situation the algorithm
  !! scales as N**2, but assuming an almost sorted array it can even
  !! reach N.  Generally, the insertion sort should not be used
  !! for N >= 100. This implementation is an in-place sort.
  !!
  SUBROUTINE insertionSort_r16(array)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(inout) :: array
    REAL(rprec16)                              :: element
    INTEGER                                    :: i, j

    DO i=2, SIZE(array)
       element = array(i)
       j = i
       DO WHILE (j > 1)
          IF (array(j-1) <= element) THEN
             EXIT
          END IF
          array(j) = array(j-1)
          j = j - 1
       END DO
       array(j) = element
    END DO

  END SUBROUTINE insertionSort_r16





  SUBROUTINE quickSort_i4(array, error)

    IMPLICIT NONE
    INTEGER(iprec4), DIMENSION(:), INTENT(inout) :: array
    LOGICAL, INTENT(inout) :: error
    INTEGER(iprec4) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(i)
             DO j=i-1,left,-1
                IF (array(j) <= element) EXIT
                array(j+1) = array(j)
             END DO
             array(j+1) = element
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(array(center),array(left+1))
          CALL swap(array(left),array(right),array(left)>array(right))
          CALL swap(array(left+1),array(right),array(left+1)>array(right))
          CALL swap(array(left),array(left+1),array(left)>array(left+1))
          i = left + 1
          j = right
          element = array(left+1)
          DO
             DO
                i = i + 1
                IF (array(i) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(j) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(array(i), array(j))
          END DO
          array(left+1) = array(j)
          array(j)      = element
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_i4





  SUBROUTINE quickSort_i4_index(array, ind, error)

    IMPLICIT NONE
    INTEGER(iprec4), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(SIZE(array)), INTENT(out) :: ind
    LOGICAL, INTENT(inout) :: error
    INTEGER(iprec4) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack, tmp

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n
    DO i=1,SIZE(ind)
       ind(i) = i
    END DO

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(ind(i))
             tmp = ind(i)
             DO j=i-1,left,-1
                IF (array(ind(j)) <= element) THEN
                   EXIT
                END IF
                ind(j+1) = ind(j)
             END DO
             ind(j+1) = tmp
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(ind(center),ind(left+1))
          CALL swap(ind(left),ind(right),array(ind(left))>array(ind(right)))
          CALL swap(ind(left+1),ind(right),array(ind(left+1))>array(ind(right)))
          CALL swap(ind(left),ind(left+1),array(ind(left))>array(ind(left+1)))
          i = left + 1
          j = right
          element = array(ind(left+1))
          tmp = ind(left+1)
          DO
             DO
                i = i + 1
                IF (array(ind(i)) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(ind(j)) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(ind(i), ind(j))
          END DO
          ind(left+1) = ind(j)
          ind(j)      = tmp
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_i4_index





  SUBROUTINE quickSort_i8(array, error)

    IMPLICIT NONE
    INTEGER(rprec8), DIMENSION(:), INTENT(inout) :: array
    LOGICAL, INTENT(inout) :: error
    INTEGER(rprec8) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(i)
             DO j=i-1,left,-1
                IF (array(j) <= element) EXIT
                array(j+1) = array(j)
             END DO
             array(j+1) = element
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(array(center), array(left+1))
          CALL swap(array(left),   array(right),  array(left)   > array(right))
          CALL swap(array(left+1), array(right),  array(left+1) > array(right))
          CALL swap(array(left),   array(left+1), array(left)   > array(left+1))
          i = left + 1
          j = right
          element = array(left+1)
          DO
             DO
                i = i + 1
                IF (array(i) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(j) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(array(i), array(j))
          END DO
          array(left+1) = array(j)
          array(j)      = element
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_i8





  SUBROUTINE quickSort_i8_index(array, ind, error)

    IMPLICIT NONE
    INTEGER(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(SIZE(array)), INTENT(out) :: ind
    LOGICAL, INTENT(inout) :: error
    INTEGER(rprec8) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack, tmp

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n
    DO i=1,SIZE(ind)
       ind(i) = i
    END DO

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(ind(i))
             tmp = ind(i)
             DO j=i-1,left,-1
                IF (array(ind(j)) <= element) THEN
                   EXIT
                END IF
                ind(j+1) = ind(j)
             END DO
             ind(j+1) = tmp
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(ind(center),ind(left+1))
          CALL swap(ind(left),ind(right),array(ind(left))>array(ind(right)))
          CALL swap(ind(left+1),ind(right),array(ind(left+1))>array(ind(right)))
          CALL swap(ind(left),ind(left+1),array(ind(left))>array(ind(left+1)))
          i = left + 1
          j = right
          element = array(ind(left+1))
          tmp = ind(left+1)
          DO
             DO
                i = i + 1
                IF (array(ind(i)) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(ind(j)) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(ind(i), ind(j))
          END DO
          ind(left+1) = ind(j)
          ind(j)      = tmp
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_i8_index





  SUBROUTINE quickSort_r8(array, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: array
    LOGICAL, INTENT(inout) :: error
    REAL(rprec8) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(i)
             DO j=i-1,left,-1
                IF (array(j) <= element) EXIT
                array(j+1) = array(j)
             END DO
             array(j+1) = element
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(array(center),array(left+1))
          CALL swap(array(left),array(right),array(left)>array(right))
          CALL swap(array(left+1),array(right),array(left+1)>array(right))
          CALL swap(array(left),array(left+1),array(left)>array(left+1))
          i = left + 1
          j = right
          element = array(left+1)
          DO
             DO
                i = i + 1
                IF (array(i) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(j) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(array(i), array(j))
          END DO
          array(left+1) = array(j)
          array(j)      = element
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_r8





  SUBROUTINE quickSort_r8_index(array, ind, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(SIZE(array)), INTENT(out) :: ind
    LOGICAL, INTENT(inout) :: error
    REAL(rprec8) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack, tmp

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n
    DO i=1,SIZE(ind)
       ind(i) = i
    END DO

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(ind(i))
             tmp = ind(i)
             DO j=i-1,left,-1
                IF (array(ind(j)) <= element) THEN
                   EXIT
                END IF
                ind(j+1) = ind(j)
             END DO
             ind(j+1) = tmp
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(ind(center),ind(left+1))
          CALL swap(ind(left),ind(right),array(ind(left))>array(ind(right)))
          CALL swap(ind(left+1),ind(right),array(ind(left+1))>array(ind(right)))
          CALL swap(ind(left),ind(left+1),array(ind(left))>array(ind(left+1)))
          i = left + 1
          j = right
          element = array(ind(left+1))
          tmp = ind(left+1)
          DO
             DO
                i = i + 1
                IF (array(ind(i)) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(ind(j)) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(ind(i), ind(j))
          END DO
          ind(left+1) = ind(j)
          ind(j)      = tmp
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_r8_index





  SUBROUTINE quickSort_r16(array, error)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(inout) :: array
    LOGICAL, INTENT(inout) :: error
    REAL(rprec16) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(i)
             DO j=i-1,left,-1
                IF (array(j) <= element) EXIT
                array(j+1) = array(j)
             END DO
             array(j+1) = element
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(array(center),array(left+1))
          CALL swap(array(left),array(right),array(left)>array(right))
          CALL swap(array(left+1),array(right),array(left+1)>array(right))
          CALL swap(array(left),array(left+1),array(left)>array(left+1))
          i = left + 1
          j = right
          element = array(left+1)
          DO
             DO
                i = i + 1
                IF (array(i) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(j) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(array(i), array(j))
          END DO
          array(left+1) = array(j)
          array(j)      = element
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_r16





  SUBROUTINE quickSort_ch(array, error)

    IMPLICIT NONE
    CHARACTER(len=*), DIMENSION(:), INTENT(inout) :: array
    LOGICAL, INTENT(inout) :: error
    CHARACTER(len=LEN(array(1))) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack

    n = SIZE(array,dim=1)
    istack = 0
    left  = 1
    right = n

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(i)
             DO j=i-1,left,-1
                IF (array(j) <= element) EXIT
                array(j+1) = array(j)
             END DO
             array(j+1) = element
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(array(center),array(left+1))
          CALL swap(array(left),array(right),array(left)>array(right))
          CALL swap(array(left+1),array(right),array(left+1)>array(right))
          CALL swap(array(left),array(left+1),array(left)>array(left+1))
          i = left + 1
          j = right
          element = array(left+1)
          DO
             DO
                i = i + 1
                IF (array(i) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(j) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(array(i), array(j))
          END DO
          array(left+1) = array(j)
          array(j)      = element
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_ch





  SUBROUTINE quickSort_ch_index(array, ind, error)

    IMPLICIT NONE
    CHARACTER(len=*), DIMENSION(:), INTENT(inout) :: array
    INTEGER, DIMENSION(SIZE(array)), INTENT(out) :: ind
    LOGICAL, INTENT(inout) :: error
    CHARACTER(len=LEN(array(1))) :: element
    INTEGER, DIMENSION(NSTACK) :: stack
    INTEGER :: i, j, center, left, right, n, istack, tmp

    n = SIZE(array)
    istack = 0
    left  = 1
    right = n
    DO i=1,SIZE(ind)
       ind(i) = i
    END DO

    DO
       IF (right-left < NN) THEN
          DO i=left+1,right
             element = array(ind(i))
             tmp = ind(i)
             DO j=i-1,left,-1
                IF (array(ind(j)) <= element) THEN
                   EXIT
                END IF
                ind(j+1) = ind(j)
             END DO
             ind(j+1) = tmp
          END DO
          IF (istack == 0) RETURN
          right  = stack(istack)
          left   = stack(istack-1)
          istack = istack - 2
       ELSE
          center = (left + right)/2
          CALL swap(ind(center),ind(left+1))
          CALL swap(ind(left),ind(right),array(ind(left))>array(ind(right)))
          CALL swap(ind(left+1),ind(right),array(ind(left+1))>array(ind(right)))
          CALL swap(ind(left),ind(left+1),array(ind(left))>array(ind(left+1)))
          i = left + 1
          j = right
          element = array(ind(left+1))
          tmp = ind(left+1)
          DO
             DO
                i = i + 1
                IF (array(ind(i)) >= element) EXIT
             END DO
             DO
                j = j - 1
                IF (array(ind(j)) <= element) EXIT
             END DO
             IF (j < i) EXIT
             CALL swap(ind(i), ind(j))
          END DO
          ind(left+1) = ind(j)
          ind(j)      = tmp
          istack = istack + 2
          IF (istack > NSTACK) THEN
             error = .TRUE.
             WRITE(0,*) 'Error: sort / quicksort: NSTACK too small...'
             RETURN
          END IF
          IF (right-i+1 >= j-1) THEN
             stack(istack)   = right
             stack(istack-1) = i
             right           = j - 1
          ELSE
             stack(istack)   = j - 1
             stack(istack-1) = left
             left            = i
          END IF
       END IF
    END DO

  END SUBROUTINE quickSort_ch_index





  !! *Description*:
  !!
  !! Binary search routine. Returns the position of the search key in
  !! the one-dimensional input array. If the search key is not found,
  !! it returns -1.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the search
  !! key by repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the key, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the value is found, or the interval
  !! is empty.
  !!
  INTEGER FUNCTION binarySearch_i8(key, array)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    INTEGER(iprec8), INTENT(in)               :: key
    INTEGER(iprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER                              :: n, left, right, center, i

    n = SIZE(array, dim=1)
    ! Return immediately
    ! - if the length of the array is zero:
    IF (n == 0) THEN
       binarySearch_i8 = -1
       RETURN
    END IF
    ! - if the key (value to be searched for) is smaller or 
    !   larger than the minimum or maximum values of
    !   the array:
    IF (key < array(1) .OR. key > array(n)) THEN
       binarySearch_i8 = -1
       RETURN
    END IF

    i = 1
    left = 1
    right = n
    DO WHILE (left <= right .AND. i<=n+2)
       i = i + 1
       center = CEILING((left+right)/2.0)
       IF (key == array(center)) THEN
          binarySearch_i8 = center
          RETURN
       ELSE IF (key < array(center)) THEN
          right = center - 1 
       ELSE IF (key > array(center)) THEN
          left = center + 1
       END IF
    END DO

    IF (i == n+2) THEN
       WRITE(0,*) "ERROR! Binary search stuck in a loop with key:"
       WRITE(0,*) key
       WRITE(0,*) "and array:"
       WRITE(0,*) array
    END IF

    ! The key could not be found in the array. 
    ! Return a negative index:
    binarySearch_i8 = -1

  END FUNCTION binarySearch_i8





  !! *Description*:
  !!
  !! Binary search routine. Returns the position of the search key in
  !! the one-dimensional input array. If the search key is not found,
  !! it returns -1.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the search
  !! key by repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the key, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the value is found, or the interval
  !! is empty.
  !!
  INTEGER FUNCTION binarySearch_r8(key, array)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    REAL(rprec8), INTENT(in)               :: key
    REAL(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER                           :: n, left, right, center, i

    n = SIZE(array, dim=1)
    ! Return immediately
    ! - if the length of the array is zero:
    IF (n == 0) THEN
       binarySearch_r8 = -1
       RETURN
    END IF
    ! - if the key (value to be searched for) is smaller or 
    !   larger than the minimum or maximum values of
    !   the array:
    IF (key < array(1) .OR. key > array(n)) THEN
       binarySearch_r8 = -1
       RETURN
    END IF

    i = 1
    left = 1
    right = n
    DO WHILE (left <= right .AND. i<=n+2)
       i = i + 1
       center = CEILING((left+right)/2.0)
       IF (ABS(key - array(center)) < 10.0_rprec8*EPSILON(key)) THEN
          binarySearch_r8 = center
          RETURN
       ELSE IF (key < array(center)) THEN
          right = center - 1 
       ELSE IF (key > array(center)) THEN
          left = center + 1
       END IF
    END DO

    IF (i == n+2) THEN
       WRITE(0,*) "ERROR! Binary search stuck in a loop with key:"
       WRITE(0,*) key
       WRITE(0,*) "and array:"
       WRITE(0,*) array
    END IF

    ! The key could not be found in the array. 
    ! Return a negative index:
    binarySearch_r8 = -1

  END FUNCTION binarySearch_r8





  !! *Description*:
  !!
  !! Binary search routine. Returns the position of the search key in
  !! the one-dimensional input array. If the search key is not found,
  !! it returns -1.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the search
  !! key by repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the key, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the value is found, or the interval
  !! is empty.
  !!
  INTEGER FUNCTION binarySearch_r16(key, array)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    REAL(rprec16), INTENT(in)               :: key
    REAL(rprec16), DIMENSION(:), INTENT(in) :: array
    INTEGER                            :: n, left, right, center, i

    n = SIZE(array, dim=1)
    ! Return immediately
    ! - if the length of the array is zero:
    IF (n == 0) THEN
       binarySearch_r16 = -1
       RETURN
    END IF
    ! - if the key (value to be searched for) is smaller or 
    !   larger than the minimum or maximum values of
    !   the array:
    IF (key < array(1) .OR. key > array(n)) THEN
       binarySearch_r16 = -1
       RETURN
    END IF

    i = 1
    left = 1
    right = n
    DO WHILE (left <= right .AND. i<=n+2)
       i = i + 1
       center = CEILING((left+right)/2.0)
       IF (ABS(key - array(center)) < 10.0_rprec16*EPSILON(key)) THEN
          binarySearch_r16 = center
          RETURN
       ELSE IF (key < array(center)) THEN
          right = center - 1 
       ELSE IF (key > array(center)) THEN
          left = center + 1
       END IF
    END DO

    IF (i == n+2) THEN
       WRITE(0,*) "ERROR! Binary search stuck in a loop with key:"
       WRITE(0,*) key
       WRITE(0,*) "and array:"
       WRITE(0,*) array
    END IF

    ! The key could not be found in the array. 
    ! Return a negative index:
    binarySearch_r16 = -1

  END FUNCTION binarySearch_r16





  !! *Description*:
  !!
  !! Binary search routine. Returns the position of the search key in
  !! the one-dimensional input array. If the search key is not found,
  !! it returns -1.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the search
  !! key by repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the key, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the value is found, or the interval
  !! is empty.
  !!
  INTEGER FUNCTION binarySearch_ch_index(key, array, indx)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)               :: key
    CHARACTER(len=*), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(:), INTENT(in) :: indx
    INTEGER                            :: n, left, right, center, i

    n = SIZE(array, dim=1)
    ! Return immediately
    ! - if the length of the array is zero:
    IF (n == 0) THEN
       binarySearch_ch_index = -1
       RETURN
    END IF
    ! - if the key (value to be searched for) is smaller or 
    !   larger than the minimum or maximum values of
    !   the array:
    IF (key < array(indx(1)) .OR. key > array(indx(n))) THEN
       binarySearch_ch_index = -1
       RETURN
    END IF

    i = 1
    left = 1
    right = n
    DO WHILE (left <= right .AND. i<=n+2)
       i = i + 1
       center = CEILING((left+right)/2.0)
       IF (key == array(indx(center))) THEN
          binarySearch_ch_index = indx(center)
          RETURN
       ELSE IF (key < array(indx(center))) THEN
          right = center - 1 
       ELSE IF (key > array(indx(center))) THEN
          left = center + 1
       END IF
    END DO

    IF (i == n+2) THEN
       WRITE(0,*) "ERROR! Binary search stuck in a loop with key:"
       WRITE(0,*) key
       WRITE(0,*) "and array:"
       WRITE(0,*) array
    END IF

    ! The key could not be found in the array. 
    ! Return a negative index:
    binarySearch_ch_index = -1

  END FUNCTION binarySearch_ch_index





  !! *Description*:
  !!
  !! Binary search routine. Returns the position where the value
  !! should be found (starting from left) in the one-dimensional input
  !! array, if the array is sorted in ascending order.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the value by
  !! repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the value, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the interval is empty.
  !!
  INTEGER FUNCTION findLocation_r8(value, array)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    REAL(rprec8), INTENT(in)               :: value
    REAL(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER                              :: n, left, right, center

    n = SIZE(array)
    ! Return immediately, if the value is smaller or larger than the
    ! minimum or maximum values of the array:
    IF (value < array(1)) THEN
       findLocation_r8 = 1
       RETURN
    ELSE IF (value >= array(n)) THEN
       findLocation_r8 = n+1
       RETURN
    END IF

    left = 1
    right = n
    center = CEILING((left+right)/2.0)
    DO WHILE (right-left > 1)
       center = CEILING((left+right)/2.0)
       IF (value < array(center)) THEN
          right = center
       ELSE IF (value >= array(center)) THEN
          left = center
       END IF
    END DO
    IF (value < array(center)) THEN
       findLocation_r8 = center
    ELSE
       findLocation_r8 = center + 1
    END IF

  END FUNCTION findLocation_r8





  !! *Description*:
  !!
  !! Binary search routine. Returns the position where the value
  !! should be found (starting from left) in the one-dimensional input
  !! array, given an index array providing the ascending order of the
  !! array elements.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the value by
  !! repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the value, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the interval is empty.
  !!
  INTEGER FUNCTION findLocation_r8_indx(value, array, indx_array)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    INTEGER, PARAMETER                   :: prec = 8
    REAL(rprec8), INTENT(in)               :: value
    REAL(rprec8), DIMENSION(:), INTENT(in) :: array
    INTEGER, DIMENSION(:), INTENT(in)    :: indx_array
    INTEGER                              :: n, left, right, center

    n = SIZE(array)
    ! Return immediately, if the value is smaller or larger than the
    ! minimum or maximum values of the array:
    IF (value < array(indx_array(1))) THEN
       findLocation_r8_indx = 1
       RETURN
    ELSE IF (value >= array(indx_array(n))) THEN
       findLocation_r8_indx = n+1
       RETURN
    END IF

    left = 1
    right = n
    center = CEILING((left+right)/2.0)
    DO WHILE (right-left > 1)
       center = CEILING((left+right)/2.0)
       IF (value < array(indx_array(center))) THEN
          right = center
       ELSE IF (value >= array(indx_array(center))) THEN
          left = center
       END IF
    END DO
    IF (value < array(indx_array(center))) THEN
       findLocation_r8_indx = center
    ELSE
       findLocation_r8_indx = center + 1
    END IF

  END FUNCTION findLocation_r8_indx





  !! *Description*:
  !!
  !! Binary search routine. Returns the position where the value
  !! should be found (starting from left) in the one-dimensional input
  !! array, if the array is sorted in ascending order.
  !!
  !! As a prerequisite, the array which is to be scanned has to be sorted
  !! in ascending order. The algorithm searches for the value by
  !! repeatedly dividing the search interval in half. Initially
  !! the interval covers the whole array. If the item in the middle of
  !! the interval is larger than the value, the interval is narrowed to
  !! the lower half. Otherwise it is narrowed to the upper half. The
  !! procedure is continued until the interval is empty.
  !!
  INTEGER FUNCTION findLocation_r16(value, array)

    ! The input array needs to be sorted in ascending order.

    IMPLICIT NONE
    REAL(rprec16), INTENT(in)               :: value
    REAL(rprec16), DIMENSION(:), INTENT(in) :: array
    INTEGER                              :: n, left, right, center

    n = SIZE(array)
    ! Return immediately, if the value is smaller or larger than the
    ! minimum or maximum values of the array:
    IF (value < array(1)) THEN
       findLocation_r16 = 1
       RETURN
    ELSE IF (value >= array(n)) THEN
       findLocation_r16 = n+1
       RETURN
    END IF

    left = 1
    right = n
    center = CEILING((left+right)/2.0)
    DO WHILE (right-left > 1)
       center = CEILING((left+right)/2.0)
       IF (value < array(center)) THEN
          right = center
       ELSE IF (value >= array(center)) THEN
          left = center
       END IF
    END DO
    IF (value < array(center)) THEN
       findLocation_r16 = center
    ELSE
       findLocation_r16 = center + 1
    END IF

  END FUNCTION findLocation_r16





END MODULE sort
