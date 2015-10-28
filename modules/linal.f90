!====================================================================!
!                                                                    !
! Copyright 2002-2014,2015                                           !
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
!! This is a linear algebra module written in Fortran 90. For additional   
!! information see the README.txt-file and the linal_doc.ps-file, which are
!! distributed together with this sourcefile (linal.f90).                  
!!                                                                         
!! Following routines are included:                                        
!!                                                                         
!! eigen_decomposition_jacobi - Solves the eigenproblem.                   
!!                                                                         
!! sq_matrix_check        - Checks wheter the given matrix is a square     
!!                          matrix, if the ranges are ok etc.              
!!                                                                         
!! tri_matrix_check       - Checks wheter the given matrix is a compressed 
!!                          tridiagonal matrix, if the ranges are ok etc.  
!!                                                                         
!! vector_check           - Checks wheter the vector ranges are ok.        
!!                                                                         
!! LU_factor              - Finds the LU factorization of a given square   
!!                          matrix.                                        
!!                                                                         
!! LU_solve               - Solves a linear system by using the LU factori-
!!                          zation.                                        
!!                                                                         
!! LU_separation          - Separates an LU factorization matrix to a lower
!!                          and an upper triangular matrix.                
!!                                                                         
!! matinv                 - Finds the inverse matrix of a given square     
!!                          matrix by means of the LU factorization.       
!!                                                                         
!! Gauss_elimination      - Solves a linear system by means of the Gauss   
!!                          elimination procedure.                         
!!                                                                         
!! tridiagonal_solve      - Solves a tridiagonal linear system.            
!!                                                                         
!! Jacobi_iteration       - Solves a linear system by means of the Jacobi  
!!                          iteration procedure.                           
!!                                                                         
!! Gauss_Seidel_iteration - Solves a linear system by means of the Gauss-  
!!                          Seidel iteration procedure.                    
!!                                                                         
!! determinant            - Finds the determinant of a given matrix.       
!!                                                                         
!! cond_number            - Finds the condition number of a given matrix.  
!!                                                                         
!! matrix_norm            - Finds the matrix norm of a given matrix.       
!!                                                                         
!! matrix_print          - Writes the given matrix to the given i/o-unit  
!!                          using a nice layout.                           
!!                                                                         
!! vector_print          - Writes the given vector to the given i/o-unit  
!!                          using a nice layout.                           
!!                                                                         
!! outer_product          - Computes the outer product of two given vectors.
!!                                                                         
!! cross_product          - Computes the cross product of two given vectors.
!!                                                                         
!! identity_matrix        - Returns the identity matrix of a given         
!!                          dimension.                                     
!!                                                                         
!!                                                                         
!!                                                                         
!!                                                                         
!!   THE tri_matrix_check ROUTINE SHOULD ALWAYS BE CALLED BEFORE CALLING   
!!                                                                         
!!   tridiagonal_solve                                                     
!!                                                                         
!!   THE sq_matrix_check ROUTINE SHOULD ALWAYS BE CALLED BEFORE CALLING    
!!                                                                         
!!   LU_factor                                                             
!!   LU_solve                                                              
!!   LU_separation                                                         
!!   Gauss_elimination                                                     
!!   Jacobi_iteration                                                      
!!   Gauss_Seidel_iteration                                                
!!   matinv                                                                
!!   determinant                                                           
!!   cond_number                                                           
!!                                                                         
!!   YOU DON'T NEED TO CALL THE CHECKING ROUTINES, BUT NOT CALLING THEM IS 
!!   ADVISABLE ONLY IN CASE YOU KNOW EXACTLY WHAT YOU ARE DOING.           
!!
!!
!! @author  MG
!! @version 2015-10-28
!!
MODULE linal

  USE parameters
  USE utilities

  IMPLICIT NONE

  !! Parameters:
  !!
  !!  max_iter         - Maximum number of iteration loops.
  !!
  !!  layout           - Layout of one number.
  !!

  INTEGER, PARAMETER           :: max_iter = 1000000
  CHARACTER(len=30), PARAMETER :: frmt_default = '(E22.14,1X)'

  PRIVATE :: matrix_print_r8
  PRIVATE :: matrix_print_r16
  PRIVATE :: outer_product_r4
  PRIVATE :: outer_product_r8
  PRIVATE :: outer_product_r16

  INTERFACE diagonal_multiplication
     MODULE PROCEDURE diagonal_multiplication_vec_r8
     MODULE PROCEDURE diagonal_multiplication_sca_r8
  END INTERFACE diagonal_multiplication

  INTERFACE identity_matrix
     MODULE PROCEDURE identity_matrix_r8
  END INTERFACE identity_matrix

  INTERFACE LU_factor
     MODULE PROCEDURE LU_factor_r8
     MODULE PROCEDURE LU_factor_r16
  END INTERFACE LU_factor

  INTERFACE LU_solve
     MODULE PROCEDURE LU_solve_r8
     MODULE PROCEDURE LU_solve_r16
  END INTERFACE LU_solve

  INTERFACE matrix_print
     MODULE PROCEDURE matrix_print_r8
     MODULE PROCEDURE matrix_print_r16
  END INTERFACE matrix_print

  INTERFACE matinv
     MODULE PROCEDURE matinv_r8
     MODULE PROCEDURE matinv_r16
  END INTERFACE matinv

  INTERFACE outer_product
     MODULE PROCEDURE outer_product_r4
     MODULE PROCEDURE outer_product_r8
     MODULE PROCEDURE outer_product_r16
  END INTERFACE outer_product

  INTERFACE triple_product
     MODULE PROCEDURE triple_product_r8
  END INTERFACE triple_product

  INTERFACE cross_product
     MODULE PROCEDURE cross_product_r8
  END INTERFACE cross_product

CONTAINS



  !! *Description*:
  !!
  !! Cholesky decomposition.
  !!
  !! Given an N x N positive-definite symmetric matrix 'a', this routine
  !! constructs its Cholesky decomposition, A = L · Transpose(L) . On
  !! input, only the upper triangle of 'a' needs to be given; it is not
  !! modified. The Cholesky factor L is returned in the lower triangle
  !! of 'a', except for its diagonal elements, which are returned in 'p',
  !! a vector of length N.
  !!
  SUBROUTINE cholesky_decomposition(a, p, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(inout) :: a
    REAL(rprec8), DIMENSION(:), INTENT(out) :: p
    CHARACTER(len=*), INTENT(inout) :: error

    REAL(rprec8) :: summ
    INTEGER(iprec4) :: i, n

    n = SIZE(a,dim=1)
    IF (n /= SIZE(a,dim=2) .OR. n /= SIZE(p)) THEN
       error = " -> linal : cholesky_decomposition : Matrix and vector sizes are not compatible." // &
            TRIM(error)
       RETURN
    END IF
    DO i=1,n
       summ = a(i,i) - DOT_PRODUCT(a(i,1:i-1),a(i,1:i-1))
       IF (summ <= 0.0_rprec8) THEN 
          ! a, WITH rounding errors, is not positive definite
          error = " -> linal : cholesky_decomposition : Matrix is not positive definite." // &
               TRIM(error)
          RETURN
       END IF
       p(i) = SQRT(summ)
       a(i+1:n,i) = (a(i,i+1:n) - MATMUL(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
    END DO

  END SUBROUTINE cholesky_decomposition





  !! *Description*:
  !!
  !! Solves the set of N linear equations A · x = b, where a is a
  !! positive-definite symmetric matrix. a (N × N) and p (of length N)
  !! are input as the output of the routine choldc.  Only the lower
  !! subdiagonal portion of a is accessed. b is the input right-hand-
  !! side vector, of length N. The solution vector, also of length N,
  !! is returned in x. a and p are not modified and can be left in
  !! place for successive calls with different right-hand sides b. b
  !! is not modified unless you identify b and x in the calling
  !! sequence, which is allowed.
  !!
  SUBROUTINE cholesky_solve(a, p, b, x, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: a
    REAL(rprec8), DIMENSION(:), INTENT(in) :: p, b
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: x
    CHARACTER(len=*), INTENT(inout) :: error

    INTEGER(iprec4) :: i, n

    n = SIZE(a,dim=1)
    IF (n /= SIZE(a,dim=2) .OR. n /= SIZE(p) .OR. &
         n /= SIZE(b) .OR. n /= SIZE(x)) THEN
       error = " -> linal : cholesky_solve : Matrix and/or vector sizes are not compatible." // &
            TRIM(error)
       RETURN
    END IF
    DO i=1,n ! Solve L · y = b, storing y in x.
       x(i) = (b(i) - DOT_PRODUCT(a(i,1:i-1),x(1:i-1))) / p(i)
    END DO
    DO i=n,1,-1 ! Solve LT · x = y.
       x(i) = (x(i) - DOT_PRODUCT(a(i+1:n,i),x(i+1:n))) / p(i)
    END DO

  END SUBROUTINE cholesky_solve





  FUNCTION diagonal_multiplication_vec_r8(matrix, diagonal, error, mask)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: matrix
    REAL(rprec8), DIMENSION(:), INTENT(in) :: diagonal
    REAL(rprec8), DIMENSION(SIZE(matrix,dim=1),SIZE(matrix,dim=2)) :: diagonal_multiplication_vec_r8
    CHARACTER(len=*), INTENT(inout) :: error
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: mask 
    INTEGER :: i, n

    n = MIN(SIZE(matrix,dim=1),SIZE(matrix,dim=2))
    IF (SIZE(diagonal) < n) THEN
       error = " -> linal : diagonal_multiplication : Matrix and diagonal vector are not compatible." // &
            TRIM(error)
       RETURN
    END IF
    IF (PRESENT(mask)) THEN
       IF (SIZE(mask) < n) THEN
          error = " -> linal : diagonal_multiplication : Matrix and mask vector are not compatible." // &
               TRIM(error)
          RETURN
       END IF
    END IF
    diagonal_multiplication_vec_r8 = matrix
    DO i=1,n
       IF (PRESENT(mask)) THEN
          IF (.NOT.mask(i)) THEN
             CYCLE
          END IF
       END IF
       diagonal_multiplication_vec_r8(i,i) = matrix(i,i) * diagonal(i)
    END DO

  END FUNCTION diagonal_multiplication_vec_r8





  FUNCTION diagonal_multiplication_sca_r8(matrix, diagonal, error, mask)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: matrix
    REAL(rprec8), INTENT(in) :: diagonal
    REAL(rprec8), DIMENSION(SIZE(matrix,dim=1),SIZE(matrix,dim=2)) :: diagonal_multiplication_sca_r8
    CHARACTER(len=*), INTENT(inout) :: error
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in) :: mask 
    INTEGER :: i, n

    n = MIN(SIZE(matrix,dim=1),SIZE(matrix,dim=2))
    IF (PRESENT(mask)) THEN
       IF (SIZE(mask) < n) THEN
          error = " -> linal : diagonal_multiplication : Matrix and mask vector are not compatible." // &
               TRIM(error)
          RETURN
       END IF
    END IF
    diagonal_multiplication_sca_r8 = matrix
    DO i=1,n
       IF (PRESENT(mask)) THEN
          IF (.NOT.mask(i)) THEN
             CYCLE
          END IF
       END IF
       diagonal_multiplication_sca_r8(i,i) = matrix(i,i) * diagonal
    END DO

  END FUNCTION diagonal_multiplication_sca_r8





  !! *Description*:
  !!
  !! Computes all eigenvalues and eigenvectors of a real symmetric
  !! N x N matrix 'a'. 'd' is a vector of length N that returns the
  !! eigenvalues of 'a'. 'v' is an N x N matrix whose columns contain,
  !! on output, the normalized eigenvectors of 'a'. 'nrot' returns the
  !! number of Jacobi rotations that were required.
  !! 
  !! Reference:
  !! Press, Teukolsky, Vetterling, and Flannery: Numerical Recipes in
  !! Fortran 90 -- The Art of Parallel Scientific Computing, 2nd Ed.,
  !! 1999, Cambridge University Press, pp. 1225--1227.
  !!
  SUBROUTINE eigen_decomposition_jacobi(a, d, v, nrot, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)  :: a
    REAL(rprec8), DIMENSION(:), INTENT(out)   :: d
    REAL(rprec8), DIMENSION(:,:), INTENT(out) :: v
    INTEGER, INTENT(out)                      :: nrot
    CHARACTER(len=*), INTENT(inout)           :: error

    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: aa
    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: b, z
    REAL(rprec8) :: c, g, h, s, sm, t, tau, theta, tresh
    INTEGER :: i, j, ip, iq, n, err
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: upper_triangle

    n = SIZE(a,dim=1)
    IF (n /= SIZE(a,dim=2) .OR. n /= SIZE(d) .OR. &
         n /= SIZE(v,dim=1) .OR. n /= SIZE(v,dim=2)) THEN
       error = " -> linal : eigen_decomposition_jacobi : Matrix and/or vector sizes are not compatible." // &
            TRIM(error)
       RETURN
    END IF
    ALLOCATE(aa(n,n), b(n), z(n), upper_triangle(n,n), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : eigen_decomposition_jacobi : Could not allocate memory." // &
            TRIM(error)
       DEALLOCATE(aa, stat=err)
       DEALLOCATE(b, stat=err)
       DEALLOCATE(z, stat=err)
       DEALLOCATE(upper_triangle, stat=err)
       RETURN
    END IF
    aa = a
    upper_triangle = .FALSE.
    DO i=1,n-1
       upper_triangle(i,i+1:n) = .TRUE.
       b(i) = aa(i,i)
    END DO
    v = identity_matrix(n)
    d(:) = b(:)
    z(:) = 0.0_rprec8
    nrot = 0
    DO i=1, 100
       sm = SUM(ABS(aa), mask=upper_triangle)
       IF (sm == 0.0_rprec8) THEN
          DEALLOCATE(aa, b, z, upper_triangle, stat=err)
          IF (err /= 0) THEN
             error = " -> linal : eigen_decomposition_jacobi : Could not deallocate memory (1)." // &
                  TRIM(error)
             DEALLOCATE(aa, stat=err)
             DEALLOCATE(b, stat=err)
             DEALLOCATE(z, stat=err)
             DEALLOCATE(upper_triangle, stat=err)
          END IF
          ! Make sure that the eigenvalues are non-negative:
          d = ABS(d)
          RETURN
       END IF
       tresh = MERGE((0.2_rprec8*sm)/(n**2), 0.0_rprec8, i < 4)
       DO ip=1,n-1
          DO iq=ip+1,n
             g = 100.0_rprec8*ABS(aa(ip,iq))
             IF((i > 4) .AND. (ABS(d(ip))+g == ABS(d(ip))) &
                  .AND. (ABS(d(iq))+g == ABS(d(iq)))) THEN
                aa(ip,iq) = 0.0_rprec8
             ELSE IF (ABS(aa(ip,iq)) > tresh) THEN
                h = d(iq) - d(ip)
                IF(ABS(h)+g == ABS(h)) THEN
                   t = aa(ip,iq)/h
                ELSE
                   theta = 0.5_rprec8*h/aa(ip,iq)
                   t = 1.0_rprec8/(ABS(theta)+SQRT(1.0_rprec8+theta**2.0_rprec8))
                   IF(theta < 0.0) THEN
                      t = -t
                   END IF
                END IF
                c = 1.0_rprec8 / SQRT(1.0_rprec8+t**2.0_rprec8)
                s = t*c
                tau = s/(1.0_rprec8 + c)
                h = t*aa(ip,iq)
                z(ip) = z(ip) - h
                z(iq) = z(iq) + h
                d(ip) = d(ip) - h
                d(iq) = d(iq) + h
                aa(ip,iq) = 0.0_rprec8
                CALL jrotate(aa(1:ip-1, ip), aa(1:ip-1, iq), error)
                IF (LEN_TRIM(error) /= 0) THEN
                   error = " -> linal : eigen_decomposition_jacobi : (1)" // TRIM(error) 
                   RETURN
                END IF
                CALL jrotate(aa(ip, ip+1:iq-1), aa(ip+1:iq-1, iq), error)
                IF (LEN_TRIM(error) /= 0) THEN
                   error = " -> linal : eigen_decomposition_jacobi : (2)" // TRIM(error) 
                   RETURN
                END IF
                CALL jrotate(aa(ip, iq+1:n), aa(iq, iq+1:n), error)
                IF (LEN_TRIM(error) /= 0) THEN
                   error = " -> linal : eigen_decomposition_jacobi : (3)" // TRIM(error) 
                   RETURN
                END IF
                CALL jrotate(v(:,ip), v(:,iq), error)
                IF (LEN_TRIM(error) /= 0) THEN
                   error = " -> linal : eigen_decomposition_jacobi : (4)" // TRIM(error) 
                   RETURN
                END IF
                IF (LEN_TRIM(error) /= 0) THEN
                   DEALLOCATE(aa, stat=err)
                   DEALLOCATE(b, stat=err)
                   DEALLOCATE(z, stat=err)
                   DEALLOCATE(upper_triangle, stat=err)
                   RETURN
                END IF
                nrot = nrot + 1
             END IF
          END DO
       END DO
       b(:) = b(:) + z(:)
       d(:) = b(:)
       z(:) = 0.0_rprec8
    END DO
    error = " -> linal : eigen_decomposition_jacobi : Too many iterations." // &
         TRIM(error)
    DEALLOCATE(aa, b, z, upper_triangle, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : eigen_decomposition_jacobi : Could not deallocate memory (2)." // &
            TRIM(error)
       DEALLOCATE(aa, stat=err)
       DEALLOCATE(b, stat=err)
       DEALLOCATE(z, stat=err)
       DEALLOCATE(upper_triangle, stat=err)
       RETURN
    END IF

  CONTAINS

    SUBROUTINE jrotate(a1, a2, error)

      IMPLICIT NONE
      REAL(rprec8), DIMENSION(:), INTENT(inout) :: a1, a2
      CHARACTER(len=*), INTENT(inout)           :: error

      REAL(rprec8), DIMENSION(:), ALLOCATABLE :: wk1
      INTEGER :: err = 0

      ALLOCATE(wk1(SIZE(a1)), stat=err)
      IF (err /= 0) THEN
         error = " -> linal : jrotate : Could not allocate memory." // &
              TRIM(error)
         RETURN
      END IF

      wk1(:) = a1(:)
      a1(:) = a1(:) - s*(a2(:) + a1(:) * tau)
      a2(:) = a2(:) + s*(wk1(:) - a2(:) * tau)

      DEALLOCATE(wk1, stat=err)
      IF (err /= 0) THEN
         error = " -> linal : jrotate : Could not deallocate memory." // &
              TRIM(error)
         RETURN
      END IF

    END SUBROUTINE jrotate

  END SUBROUTINE eigen_decomposition_jacobi





  !! *Description*:
  !!
  !! Returns an identity matrix of a given dimension. 
  !!
  FUNCTION identity_matrix_r8(n)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    REAL(rprec8), DIMENSION(n,n) :: identity_matrix_r8
    INTEGER :: i

    identity_matrix_r8 = 0.0_rprec8
    DO i=1,n
       identity_matrix_r8(i,i) = 1.0_rprec8
    END DO

  END FUNCTION identity_matrix_r8





  SUBROUTINE sq_matrix_check(A, n, error)

    !! Finds out the following things about the given matrix:
    !!    - is it a square matrix (n,n)? (yes=ok)
    !!
    !! Input:  - Matrix                                      A(:,:)
    !! Output: - Number of columns/rows (if square, else 0)  n
    !!         - Error (.false., if ok)                      error

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: A
    INTEGER, INTENT(out) :: n
    CHARACTER(len=*), INTENT(inout) :: error
    INTEGER, DIMENSION(2) :: up_bound

    n = 0
    up_bound = UBOUND(A) 

    ! Is A a square matrix?
    IF (up_bound(1) /= up_bound(2)) THEN
       error = " -> linal : sq_matrix_check : A is not a square matrix." // &
            TRIM(error)
       RETURN
       ! Is it bigger than 1x1?
    ELSE IF (up_bound(1) < 2) THEN
       error = " -> linal : sq_matrix_check : A is a 1x1 matrix." // &
            TRIM(error)
       RETURN
       ! If everything is ok, give n the right value:
    ELSE
       n = up_bound(1)
    END IF

  END SUBROUTINE sq_matrix_check





  SUBROUTINE tri_matrix_check(A, n, error)

    !! Finds out if the given matrix is a compressed tridiagonal
    !! matrix (n,3)? (yes=ok)
    !!
    !! Input:  - Matrix                                      A(:,:)
    !! Output: - Number of columns/rows (if square, else 0)  n
    !!         - Error (.false., if ok)                      error

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: A
    INTEGER, INTENT(out) :: n
    CHARACTER(len=*), INTENT(inout) :: error
    INTEGER, DIMENSION(2) :: up_bound

    n = 0
    up_bound = UBOUND(A) 

    ! Is A a tridiagonal matrix?
    IF (up_bound(2) /= 3) THEN
       error = " -> linal : tri_matrix_check : A is not a compressed tridiagonal matrix." // &
            TRIM(error)
       RETURN
       ! If everything is ok, give n the right value:
    ELSE
       n = up_bound(1)
    END IF

  END SUBROUTINE tri_matrix_check





  SUBROUTINE vector_size(b, n)

    !! Finds out the number of elements of the given vector.
    !!
    !! Input:  - Vector                     b(:)
    !! Output: - Size of vector             n

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: b
    INTEGER, INTENT(out)                   :: n

    INTEGER, DIMENSION(1) :: up_bound

    up_bound = UBOUND(b) 
    n = up_bound(1)

  END SUBROUTINE vector_size





  SUBROUTINE LU_factor_r8(A2LU, indx, error)

    !! Finds the LU factorization of a given square matrix.
    !!
    !! Input:  - Square matrix              A2LU(n,n)
    !! Output: - LU factorization           A2LU(n,n)
    !!         - Index vector               indx(n)
    !!         - Error (.false., if ok)     error

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(inout)       :: A2LU
    INTEGER, DIMENSION(SIZE(A2LU,dim=1)), INTENT(out) :: indx
    CHARACTER(len=*), INTENT(inout)                            :: error

    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: vv, helpvec
    INTEGER :: n, i, imax, err

    n = SIZE(A2LU,dim=1)
    ALLOCATE(vv(SIZE(indx)), helpvec(SIZE(indx)), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : LU_factor : Could not allocate memory." // &
            TRIM(error)
       DEALLOCATE(vv, stat=err)
       DEALLOCATE(helpvec, stat=err)
       RETURN
    END IF

    ! Find scaling factor for each row so that the maximum 
    ! element of the row is equal to one:
    vv = MAXVAL(ABS(A2LU),dim=2)
    IF (MINVAL(vv) < 1.0e-32_rprec8*EPSILON(vv)) THEN
       error = " -> linal : LU_factor : Maximum element of a row too small." // &
            TRIM(error)
       DEALLOCATE(vv, stat=err)
       DEALLOCATE(helpvec, stat=err)
       RETURN
    ELSE
       vv = 1.0_rprec8 / vv
    END IF

    DO i=1,n
       imax = (i-1) + imaxloc(vv(i:n)*ABS(A2LU(i:n,i)))
       IF (i /= imax) THEN
          helpvec = A2LU(imax,1:n)
          A2LU(imax,1:n) = A2LU(i,1:n)
          A2LU(i,1:n) = helpvec
          vv(imax) = vv(i)
       END IF
       indx(i) = imax
       ! Avoid division with zero:
       IF (ABS(A2LU(i,i)) < 1.0e-32_rprec8*EPSILON(A2LU(i,i))) THEN
          error = " -> linal : LU_factor : Division by almost zero." // &
               TRIM(error)
          DEALLOCATE(vv, stat=err)
          DEALLOCATE(helpvec, stat=err)
          RETURN
       END IF
       A2LU(i+1:n,i) = A2LU(i+1:n,i) / A2LU(i,i)
       A2LU(i+1:n,i+1:n) = A2LU(i+1:n,i+1:n) - &
            outer_product(A2LU(i+1:n,i),A2LU(i,i+1:n))
    END DO

    DEALLOCATE(vv, helpvec, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : LU_factor : Could not deallocate memory." // &
            TRIM(error)
       DEALLOCATE(vv, stat=err)
       DEALLOCATE(helpvec, stat=err)
       RETURN
    END IF

  END SUBROUTINE LU_factor_r8





  SUBROUTINE LU_factor_r16(A2LU, indx, error)

    !! Finds the LU factorization of a given square matrix.
    !!
    !! Input:  - Square matrix              A2LU(n,n)
    !! Output: - LU factorization           A2LU(n,n)
    !!         - Index vector               indx(n)
    !!         - Error (.false., if ok)     error

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:,:), INTENT(inout)    :: A2LU
    INTEGER, DIMENSION(SIZE(A2LU,dim=1)), INTENT(out) :: indx
    CHARACTER(len=*), INTENT(inout)                            :: error

    REAL(rprec16), DIMENSION(:), ALLOCATABLE :: vv, helpvec
    INTEGER :: n, i, imax, err

    n = SIZE(A2LU,dim=1)
    ALLOCATE(vv(SIZE(indx)), helpvec(SIZE(indx)), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : LU_factor : Could not allocate memory." // &
            TRIM(error)
       DEALLOCATE(vv, stat=err)
       DEALLOCATE(helpvec, stat=err)
       RETURN
    END IF

    ! Find scaling factor for each row so that the maximum 
    ! element of the row is equal to one:
    vv = MAXVAL(ABS(A2LU),dim=2)
    IF (MINVAL(vv) < 1.0e-32_rprec16*EPSILON(vv)) THEN
       error = " -> linal : LU_factor : Maximum element of a row too small." // &
            TRIM(error)
       DEALLOCATE(vv, stat=err)
       DEALLOCATE(helpvec, stat=err)
       RETURN
    ELSE
       vv = 1.0_rprec16 / vv
    END IF

    DO i=1,n
       imax = (i-1) + imaxloc(vv(i:n)*ABS(A2LU(i:n,i)))
       IF (i /= imax) THEN
          helpvec = A2LU(imax,1:n)
          A2LU(imax,1:n) = A2LU(i,1:n)
          A2LU(i,1:n) = helpvec
          vv(imax) = vv(i)
       END IF
       indx(i) = imax
       ! Avoid division with zero:
       IF (ABS(A2LU(i,i)) < 1.0e-32_rprec16*EPSILON(A2LU(i,i))) THEN
          error = " -> linal : LU_factor : Division by almost zero." // &
               TRIM(error)
          DEALLOCATE(vv, stat=err)
          DEALLOCATE(helpvec, stat=err)
          RETURN
       END IF
       A2LU(i+1:n,i) = A2LU(i+1:n,i) / A2LU(i,i)
       A2LU(i+1:n,i+1:n) = A2LU(i+1:n,i+1:n) - &
            outer_product(A2LU(i+1:n,i),A2LU(i,i+1:n))
    END DO

    DEALLOCATE(vv, helpvec, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : LU_factor : Could not deallocate memory." // &
            TRIM(error)
       DEALLOCATE(vv, stat=err)
       DEALLOCATE(helpvec, stat=err)
       RETURN
    END IF

  END SUBROUTINE LU_factor_r16





  SUBROUTINE LU_solve_r8(LU, indx, b)

    !! Solves a linear system by means of the LU factorization.
    !! 
    !! Input:  - LU factorization square matrix  LU(n,n)
    !!         - Row permutation index vector    indx(n)
    !!         - LUx=b                           b(n)
    !! Output: - LUx=b => solve x => b=x         b(n)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)  :: LU
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: b
    INTEGER, DIMENSION(:), INTENT(in)         :: indx

    REAL(rprec8) :: summ
    INTEGER :: n, i, ii, ll

    n = SIZE(LU,dim=1)
    ii = 0

    ! Solve 'a' (=Ux) from La=b, by means of forward
    ! substitution:
    DO i=1,n
       ll = indx(i)
       summ = b(ll)
       b(ll) = b(i)
       IF (ii /= 0) THEN
          summ = summ - DOT_PRODUCT(LU(i,ii:i-1), b(ii:i-1))
       ELSE IF (ABS(summ) > 0.0_rprec8) THEN
          ii = i
       END IF
       b(i) = summ
    END DO

    ! Solve 'x' from Ux=a, by means of backward
    ! substitution:
    DO i=n,1,-1
       b(i) = (b(i) - DOT_PRODUCT(LU(i,i+1:n), b(i+1:n)))/LU(i,i)
    END DO

  END SUBROUTINE LU_solve_r8





  SUBROUTINE LU_solve_r16(LU, indx, b)

    !! Solves a linear system by means of the LU factorization.
    !! 
    !! Input:  - LU factorization square matrix  LU(n,n)
    !!         - Row permutation index vector    indx(n)
    !!         - LUx=b                           b(n)
    !! Output: - LUx=b => solve x => b=x         b(n)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:,:), INTENT(in)  :: LU
    REAL(rprec16), DIMENSION(:), INTENT(inout) :: b
    INTEGER, DIMENSION(:), INTENT(in)          :: indx

    REAL(rprec16) :: summ
    INTEGER :: n, i, ii, ll

    n = SIZE(LU,dim=1)
    ii = 0

    ! Solve 'a' (=Ux) from La=b, by means of forward
    ! substitution:
    DO i=1,n
       ll = indx(i)
       summ = b(ll)
       b(ll) = b(i)
       IF (ii /= 0) THEN
          summ = summ - DOT_PRODUCT(LU(i,ii:i-1), b(ii:i-1))
       ELSE IF (ABS(summ) > 0.0_rprec16) THEN
          ii = i
       END IF
       b(i) = summ
    END DO

    ! Solve 'x' from Ux=a, by means of backward
    ! substitution:
    DO i=n,1,-1
       b(i) = (b(i) - DOT_PRODUCT(LU(i,i+1:n), b(i+1:n)))/LU(i,i)
    END DO

  END SUBROUTINE LU_solve_r16





  SUBROUTINE LU_separation(LU2L, U)

    !! Separates an LU factorization square matrix to a lower 
    !! and an upper triangular square matrix.
    !!
    !! Input:  - LU factorization square matrix  LU2L(n,n)
    !!         - Number of columns/rows          n
    !! Output: - Lower triangular square matrix  LU2L(n,n)
    !!         - Upper triangular square matrix  U(n,n)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(inout) :: LU2L
    REAL(rprec8), DIMENSION(:,:), INTENT(out)   :: U

    INTEGER :: n, i, j

    n = SIZE(LU2L,dim=1)

    U = 0.0_rprec8
    ! Copy the upper triangular matrix from LU2L to U:
    FORALL(i=1:n, j=1:n, i<=j) 
       U(i,j) = LU2L(i,j)
    END FORALL
    ! Transform the LU2L matrix to a lower triangular matrix:
    FORALL(i=1:n)
       LU2L(i,i) = 1.0_rprec8
    END FORALL
    FORALL(i=1:n, j=1:n, i<j)
       LU2L(i,j) = 0.0_rprec8
    END FORALL

  END SUBROUTINE LU_separation





  FUNCTION matinv_r8(A, error, method)

    !! Returns the inverse of square matrix A using either the LU
    !! factorization algorithm or the Cholesky decomposition algorithm
    !! 
    !! Input:   - Square matrix              A(n,n)
    !!          - Inverse square matrix      inv_A(n,n)
    !! Output:  - Error (.false., if ok)     error 

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)             :: A
    REAL(rprec8), DIMENSION(SIZE(A,dim=2),SIZE(A,dim=1)) :: matinv_r8
    CHARACTER(len=*), INTENT(inout)                      :: error
    CHARACTER(len=*), INTENT(in), OPTIONAL               :: method

    CHARACTER(len=16) :: method_
    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: LU, L
    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: p
    INTEGER, DIMENSION(:), ALLOCATABLE :: indx
    INTEGER :: i, n, err

    n = SIZE(A,dim=1)
    IF (n /= SIZE(A,dim=2)) THEN
       error = " -> linal : matinv : Input not a square matrix." // &
            TRIM(error)
       RETURN
    END IF

    method_ = " "
    IF (PRESENT(method)) THEN
       method_ = method
    ELSE
       method_ = "LU"
    END IF

    IF (method_ == "LU") THEN

       ALLOCATE(LU(n,n), indx(n), stat=err)
       IF (err /= 0) THEN
          error = " -> linal : matinv : Could not allocate memory." // &
               TRIM(error)
          DEALLOCATE(LU, stat=err)
          DEALLOCATE(indx, stat=err)
          RETURN
       END IF

       ! A => LU :
       LU(:,:) = A(:,:)
       CALL LU_factor(LU(:,:), indx, error)
       IF (LEN_TRIM(error) /= 0) THEN
          error = " -> linal : matinv : LU factorization unsuccessful." // &
               TRIM(error)
          DEALLOCATE(LU, stat=err)
          DEALLOCATE(indx, stat=err)
          RETURN
       END IF

       ! Initialize the identity matrix:
       matinv_r8(:,:) = 0.0_rprec8
       FORALL(i=1:n)
          matinv_r8(i,i) = 1.0_rprec8
       END FORALL

       ! Solve b from LUb = e_i and put A(:,i) = b:
       DO i=1,n
          CALL LU_solve(LU, indx, matinv_r8(1:n,i))
       END DO

       DEALLOCATE(LU, indx, stat=err)
       IF (err /= 0) THEN
          error = " -> linal : matinv : Could not deallocate memory." // &
               TRIM(error)
          DEALLOCATE(LU, stat=err)
          DEALLOCATE(indx, stat=err)
          RETURN
       END IF

    ELSE IF (method_ == "Cholesky") THEN

       ALLOCATE(L(n,n), p(n), stat=err)
       IF (err /= 0) THEN
          error = " -> linal : matinv : Could not allocate memory." // &
               TRIM(error)
          DEALLOCATE(L, stat=err)
          DEALLOCATE(p, stat=err)
          RETURN
       END IF

       L(:,:) = A(:,:)
       CALL cholesky_decomposition(L, p, error)
       IF (LEN_TRIM(error) /= 0) THEN
          error = " -> linal : matinv : ." // &
               TRIM(error)
          DEALLOCATE(L, stat=err)
          DEALLOCATE(p, stat=err)
          RETURN
       END IF

       ! Initialize the identity matrix:
       matinv_r8(:,:) = 0.0_rprec8
       FORALL(i=1:n)
          matinv_r8(i,i) = 1.0_rprec8
       END FORALL

       ! Solve b from L * Transpose(L) * b = e_i and put A(:,i) = b:
       DO i=1,n
          CALL cholesky_solve(L, p, matinv_r8(1:n,i), matinv_r8(1:n,i), error)
          IF (LEN_TRIM(error) /= 0) THEN
             error = " -> linal : matinv : ." // &
                  TRIM(error)
             DEALLOCATE(L, stat=err)
             DEALLOCATE(p, stat=err)
             RETURN
          END IF
       END DO

       DEALLOCATE(L, p, stat=err)
       IF (err /= 0) THEN
          error = " -> linal : matinv : Could not deallocate memory." // &
               TRIM(error)
          DEALLOCATE(L, stat=err)
          DEALLOCATE(p, stat=err)
          RETURN
       END IF

    ELSE

       ! No such option...
       error = " -> linal : matinv : No such option." // &
            TRIM(error)

    END IF

  END FUNCTION matinv_r8





  FUNCTION matinv_r16(A, error)

    !! Returns the inverse matrix of square matrix. 
    !! 
    !! Input:   - Square matrix              A(n,n)
    !!          - Inverse square matrix      inv_A(n,n)
    !! Output:  - Error (.false., if ok)     error 

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:,:), INTENT(in)             :: A
    REAL(rprec16), DIMENSION(SIZE(A,dim=2),SIZE(A,dim=1)) :: matinv_r16
    CHARACTER(len=*), INTENT(inout)                                :: error

    REAL(rprec16), DIMENSION(:,:), ALLOCATABLE :: LU
    INTEGER, DIMENSION(:), ALLOCATABLE :: indx
    INTEGER :: i, n, err

    n = SIZE(A,dim=1)
    IF (n /= SIZE(A,dim=2)) THEN
       error = " -> linal : matinv : Input not a square matrix." // &
            TRIM(error)
       RETURN
    END IF

    ALLOCATE(LU(n,n), indx(n), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : matinv : Could not allocate memory." // &
            TRIM(error)
       DEALLOCATE(LU, stat=err)
       DEALLOCATE(indx, stat=err)
       RETURN
    END IF

    ! A => LU :
    LU(:,:) = A(:,:)
    CALL LU_factor(LU(:,:), indx, error)
    IF (LEN_TRIM(error) /= 0) THEN
       error = " -> linal : matinv : LU factorization unsuccessful." // &
            TRIM(error)
       DEALLOCATE(LU, stat=err)
       DEALLOCATE(indx, stat=err)
       RETURN
    END IF

    ! Initialize the identity matrix:
    matinv_r16(:,:) = 0.0_rprec16
    FORALL(i=1:n)
       matinv_r16(i,i) = 1.0_rprec16
    END FORALL

    ! Solve b from LUb = e_i and put A(:,i) = b:
    DO i=1,n
       CALL LU_solve(LU, indx, matinv_r16(1:n,i))
    END DO

    DEALLOCATE(LU, indx, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : matinv : Could not deallocate memory." // &
            TRIM(error)
       DEALLOCATE(LU, stat=err)
       DEALLOCATE(indx, stat=err)
       RETURN
    END IF

  END FUNCTION matinv_r16





  SUBROUTINE Gauss_elimination(A, b, error)

    !! Solves a linear system by means of the Gauss elimination 
    !! procedure.
    !!
    !! Input:  - Square matrix              A(n,n)
    !!         - Ax=b                       b(n)
    !!         - Number of columns/rows     n
    !! Output: - Ax=b => solve x => b=x     b(n)
    !!         - Error (.false., if ok)     error

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)             :: A
    REAL(rprec8), DIMENSION(:), INTENT(inout)            :: b
    CHARACTER(len=*), INTENT(inout)                                  :: error
    REAL(rprec8), DIMENSION(SIZE(A,dim=1),SIZE(A,dim=2)) :: AA
    REAL(rprec8), DIMENSION(SIZE(b))                     :: help_b
    INTEGER, DIMENSION(SIZE(b))                             :: ind
    REAL(rprec8)                                         :: smax, s, c 
    INTEGER                                                 :: i, ii, j, k, kk, l, m, n

    n = SIZE(A,dim=1)

    AA = A
    ! Initialize the row permutation vector ind:
    FORALL(i=1:n)
       ind(i) = i
    END FORALL

    DO i=1, n
       m = i
       ! Find maximum element among the remaining
       ! elements of this column.
       smax = AA(ind(i),i)
       DO k=i+1, n
          s = ABS(AA(ind(k),i))
          IF (s > smax) THEN
             smax = s
             m = k
          END IF
       END DO
       l = ind(i)
       ind(i) = ind(m)
       ind(m) = l
       ii = ind(i)
       IF (ABS(AA(ii,i)) < EPSILON(AA(ii,i))) THEN
          error = " -> linal : Gauss_elimination : Attempted division by zero." // &
               TRIM(error)
          RETURN
       END IF
       ! Calculate the upper triangular matrix by
       ! means of the Gauss procedure:
       DO k=i+1, n
          kk = ind(k)
          c = AA(kk,i) / AA(ii,i)
          DO j=i, n
             AA(kk,j) = AA(kk,j) - c * AA(ii,j)
          END DO
          b(kk) = b(kk) - c * b(ii)
       END DO
    END DO

    ! Solve the upper triangular matrix by means
    ! of backsubstitution:
    DO i=n, 1, -1
       ii = ind(i)
       b(ii) = b(ii) / AA(ii,i)
       DO k=i-1, 1, -1
          kk = ind(k)
          b(kk) = b(kk) - AA(kk,i) * b(ii)
       END DO
    END DO

    ! "Unpermutate" the solution in order to get
    ! a solution, which corresponds to the original,
    ! given matrix:  
    help_b = b
    FORALL(i=1:n)
       b(i) = help_b(ind(i))
    END FORALL

  END SUBROUTINE Gauss_elimination





  SUBROUTINE gauss_jordan(a, b, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(inout) :: a,b
    CHARACTER(len=*), INTENT(inout) :: error
    INTEGER, DIMENSION(SIZE(a,1)) :: ipiv,indxr,indxc
    LOGICAL, DIMENSION(SIZE(a,1)) :: lpiv
    REAL(rprec8) :: pivinv
    REAL(rprec8), DIMENSION(SIZE(a,1)) :: dumc
    INTEGER, DIMENSION(2), TARGET :: irc
    INTEGER :: i,l,n
    INTEGER, POINTER :: irow => NULL(), icol => NULL()

    n = SIZE(a,dim=1)
    irow => irc(1)
    icol => irc(2)
    ipiv = 0
    DO i=1,n
       lpiv = (ipiv == 0)
       irc = MAXLOC(ABS(a),outerand(lpiv,lpiv))
       ipiv(icol) = ipiv(icol) + 1
       IF (ipiv(icol) > 1) THEN
          error = " -> linal : gauss_jordan : Singular matrix (1)." // &
               TRIM(error)
          RETURN
       END IF
       IF (irow /= icol) THEN
          CALL swap(a(irow,:),a(icol,:))
          CALL swap(b(irow,:),b(icol,:))
       END IF
       indxr(i) = irow
       indxc(i) = icol
       IF (a(icol,icol) == 0.0) THEN
          error = " -> linal : gauss_jordan : Singular matrix (2)." // &
               TRIM(error)
          RETURN
       END IF
       pivinv = 1.0_rprec8/a(icol,icol)
       a(icol,icol) = 1.0_rprec8
       a(icol,:) = a(icol,:)*pivinv
       b(icol,:) = b(icol,:)*pivinv
       dumc = a(:,icol)
       a(:,icol) = 0.0
       a(icol,icol) = pivinv
       a(1:icol-1,:) = a(1:icol-1,:) - outer_product(dumc(1:icol-1),a(icol,:))
       b(1:icol-1,:) = b(1:icol-1,:) - outer_product(dumc(1:icol-1),b(icol,:))
       a(icol+1:,:) = a(icol+1:,:) - outer_product(dumc(icol+1:),a(icol,:))
       b(icol+1:,:) = b(icol+1:,:) - outer_product(dumc(icol+1:),b(icol,:))
    END DO
    DO l=n,1,-1
       CALL swap(a(:,indxr(l)),a(:,indxc(l)))
    END DO

  END SUBROUTINE gauss_jordan





  SUBROUTINE tridiagonal_solve(A, b, error)

    !! Solves a tridiagonal linear system. 
    !!
    !! Input:  - Tridigonal matrix (the non-null diagonals)  A(n,3)
    !!         - Ax=b                                        b(n)
    !! Output: - Ax=b => x=bA^(-1) => b=x                    b(n)
    !!         - Error (.false., if ok)                      error

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)  :: A
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: b
    CHARACTER(len=*), INTENT(inout)                       :: error
    REAL(rprec8), DIMENSION(SIZE(A,dim=1),3)  :: AA
    INTEGER                                      :: i, n

    n = SIZE(A,dim=1)
    AA = A

    DO i=2, n
       IF (ABS(AA(i-1,2)) < EPSILON(AA(i-1,2))) THEN
          error = " -> linal : tridiagonal_solve : Attempted division by almost zero." // &
               TRIM(error)
          RETURN
       END IF
       AA(i,2) = AA(i,2) - (AA(i,1) * AA(i-1,3)) / AA(i-1,2)
       b(i) = b(i) - (AA(i,1) * b(i-1)) / AA(i-1,2)
    END DO

    b(n) = b(n) / AA(n,2) 
    DO i=n-1, 1, -1
       b(i) = (b(i) - AA(i,3) * b(i+1)) / AA(i,2)
    END DO

  END SUBROUTINE tridiagonal_solve





  SUBROUTINE Jacobi_iteration(A, b, error)

    !! Solves a linear system by means of the Jacobi iteration
    !! procedure.
    !!
    !! Input:  - Square matrix                            A(n,n)
    !!         - Ax=b                                     b(n)
    !! Output: - Ax=b => solve x => b=x                   b(n)
    !!         - Error (.false., if ok)                   error
    !!
    !! Other:  - Old solution vector                      x(n)
    !!         - New solution vector                      y(n)
    !!         - Difference between solution vectors      diff
    !!         - Old difference between solution vectors  old_diff

    IMPLICIT NONE
    REAL(rprec8), PARAMETER                              :: iteration_factor = 1.0_rprec8
    REAL(rprec8), DIMENSION(:,:), INTENT(in)             :: A
    REAL(rprec8), DIMENSION(:), INTENT(inout)            :: b
    CHARACTER(len=*), INTENT(inout)                                  :: error
    REAL(rprec8), DIMENSION(SIZE(A,dim=1),SIZE(A,dim=2)) :: LU
    REAL(rprec8), DIMENSION(SIZE(b))                     :: D, x, y
    REAL(rprec8)                                         :: diff, old_diff
    INTEGER                                                 :: iter, i, n

    n = SIZE(A,dim=1)

    ! A = L + U + D = (L + U) + D :
    LU = A
    FORALL(i=1:n)
       D(i) = LU(i,i)
       LU(i,i) = 0.0_rprec8
    END FORALL
    x = 0.0_rprec8
    diff = HUGE(diff)
    iter = 1

    DO WHILE (diff > iteration_factor*EPSILON(diff))
       DO i=1, n
          IF (ABS(D(i)) < EPSILON(D(i))) THEN
             error = " -> linal : Jacobi_iteration : Attempted division by almost zero." // &
                  TRIM(error)
             RETURN
          END IF
          y(i) = (b(i) - DOT_PRODUCT(LU(i,:), x)) / D(i)
       END DO
       old_diff = diff
       diff = SQRT(DOT_PRODUCT(x-y, x-y))
       IF (diff > old_diff) THEN
          error = " -> linal : Jacobi_iteration : Iteration divergent." // &
               TRIM(error)
          RETURN
       END IF
       x = y
       iter = iter + 1
       IF (iter > max_iter) THEN
          error = " -> linal : Jacobi_iteration : Reached maximum number of iterations." // &
               TRIM(error)
          RETURN
       END IF
    END DO

    b = x

  END SUBROUTINE Jacobi_iteration




  SUBROUTINE Gauss_Seidel_iteration(A, b, error)

    !! Solves a linear system by means of the Gauss-Seidel iteration
    !! procedure.
    !!
    !! Input:  - Square matrix                            A(n,n)
    !!         - Ax=b                                     b(n)
    !! Output: - Ax=b => solve x => b=x                   b(n)
    !!         - Error (.false., if ok)                   error
    !!
    !! Other:  - Old solution vector                      x(n)
    !!         - New solution vector                      y(n)
    !!         - Difference between solution vectors      diff
    !!         - Old difference between solution vectors  old_diff

    IMPLICIT NONE
    REAL(rprec8), PARAMETER                               :: iteration_factor = 1.0_rprec8
    REAL(rprec8), DIMENSION(:,:), INTENT(in)              :: A
    REAL(rprec8), DIMENSION(SIZE(A,dim=1)), INTENT(inout) :: b
    CHARACTER(len=*), INTENT(inout)                                   :: error
    REAL(rprec8), DIMENSION(SIZE(A,dim=1),SIZE(A,dim=2))  :: LU
    REAL(rprec8), DIMENSION(SIZE(b))                      :: D, x, y
    REAL(rprec8)                                          :: diff, old_diff
    INTEGER                                                  :: n, iter, i
    INTRINSIC sqrt

    n = SIZE(A,dim=1)

    ! A = L + U + D = (L + U) + D :
    LU = A
    FORALL(i=1:n)
       D(i) = LU(i,i)
       LU(i,i) = 0.0_rprec8
    END FORALL

    x = 0.0_rprec8
    y = 0.0_rprec8
    diff = HUGE(diff)
    iter = 1

    DO WHILE (diff > iteration_factor*EPSILON(diff))
       DO i=1, n
          IF (ABS(D(i)) < EPSILON(D(i))) THEN
             error = " -> linal : Gauss_Seidel_iteration : Attempted division by almost zero." // &
                  TRIM(error)
             RETURN
          END IF
          y(i) = (b(i) - DOT_PRODUCT(LU(i,1:i-1), y(1:i-1)) - &
               DOT_PRODUCT(LU(i,i+1:n), x(i+1:n))) / D(i)
       END DO
       old_diff = diff
       diff = SQRT(DOT_PRODUCT(x - y, x - y))
       IF (diff > old_diff) THEN
          error = " -> linal : Gauss_Seidel_iteration : Iteration divergent." // &
               TRIM(error)
          RETURN
       END IF
       x = y
       iter = iter + 1
       IF (iter > max_iter) THEN
          error = " -> linal : Gauss_Seidel_iteration : Reached maximum number of iterations." // &
               TRIM(error)
          RETURN
       END IF
    END DO

    b = x

  END SUBROUTINE Gauss_Seidel_iteration





  REAL(rprec8) FUNCTION determinant(A, error)

    !! Returns the determinant of a square matrix.
    !!
    !! Input:  - Square matrix                   A(n,n)
    !! Output: - Error message (.false., if ok)  error

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)             :: A
    CHARACTER(len=*), INTENT(inout)                                  :: error
    REAL(rprec8), DIMENSION(SIZE(A,dim=1),SIZE(A,dim=2)) :: B
    INTEGER, DIMENSION(SIZE(A,dim=1))                       :: indx
    INTEGER                                                 :: i

    B = A

    ! A => LU :
    CALL LU_factor(B, indx, error)
    IF (LEN_TRIM(error) /= 0) THEN
       error = " -> linal : determinant : ." // TRIM(error)
       determinant = 0.0_rprec8
       RETURN
    END IF
    ! det = U(1,1)*U(2,2)*U(3,3)*...*U(n,n)
    determinant = 1.0_rprec8
    DO i=1, SIZE(B,dim=1)
       determinant = determinant * B(i,i)
    END DO

  END FUNCTION determinant





  REAL(rprec8) FUNCTION cond_nr(A, error)

    !! Finds the condition number of a given matrix. A failure in the
    !! matrix inversion results in a huge condition number.
    !!
    !! Input:  - Square matrix              A(n,n)
    !! Output: - Error (.false., if ok)     error
    !!
    !! Other:  - Matrix norm of A           A_norm
    !!         - Matrix norm of A^(-1)      inv_A_norm

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)             :: A
    CHARACTER(len=*), INTENT(inout)                               :: error
    REAL(rprec8), DIMENSION(SIZE(A,dim=2),SIZE(A,dim=1)) :: inv_A
    REAL(rprec8)                                         :: A_norm, inv_A_norm

    A_norm = matnorm(A)
    inv_A = matinv(A, error)
    IF (LEN_TRIM(error) /= 0) THEN
       error = " -> linal : cond_nr : ." // TRIM(error)
       cond_nr = HUGE(cond_nr)
       RETURN
    END IF
    inv_A_norm = matnorm(inv_A)
    cond_nr = A_norm * inv_A_norm

  END FUNCTION cond_nr





  REAL(rprec8) FUNCTION matnorm(A)

    !! Finds the matrix norm of a given matrix.
    !!
    !! Input:  - Matrix       A(:,:)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: A
    REAL(rprec8), DIMENSION(SIZE(A,dim=1))   :: s
    INTEGER                                     :: i

    ! Find the sum of each row:
    FORALL(i=1:SIZE(A,dim=1))
       s(i) = SUM(ABS(A(i,1:SIZE(A,dim=2))))
    END FORALL
    ! Matrix norm is equal to the greatest sum:
    matnorm = MAXVAL(s)

  END FUNCTION matnorm




  SUBROUTINE matrix_print_r8(A, lu, error, frmt, row_indx, column_indx)

    !! Writes the given matrix to the given i/o-unit using a nice
    !! layout.
    !!
    !! Input:  - Square matrix           A(n,n)
    !!         - I/O-unit                iounit

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)    :: A
    INTEGER, INTENT(in)                         :: lu
    CHARACTER(len=*), INTENT(inout)             :: error
    CHARACTER(len=*), INTENT(in), OPTIONAL      :: frmt
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: row_indx
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: column_indx
    CHARACTER(len=32)                           :: form
    INTEGER, DIMENSION(:), ALLOCATABLE          :: indx1, indx2
    INTEGER                                     :: i, j, siz1, siz2, err

    siz1 = SIZE(A, dim=1)
    siz2 = SIZE(A, dim=2)
    ALLOCATE(indx1(siz1), indx2(siz2), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : matrix_print : Could not allocate memory." // &
            TRIM(error)
       RETURN
    END IF

    IF (PRESENT(frmt)) THEN
       form = TRIM(frmt)
    ELSE
       form = TRIM(frmt_default)
    END IF
    IF (PRESENT(row_indx)) THEN
       indx1 = row_indx
    ELSE
       DO i=1,siz1
          indx1(i) = i
       END DO
    END IF
    IF (PRESENT(column_indx)) THEN
       indx2 = column_indx
    ELSE
       DO i=1,siz2
          indx2(i) = i
       END DO
    END IF

    ! Write the matrix row-wise:
    DO i=1,siz1
       DO j=1,siz2
          WRITE(lu,TRIM(form),advance='no') A(indx1(i),indx2(j))
       END DO
       WRITE(lu,*)
    END DO

    DEALLOCATE(indx1, indx2, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : matrix_print : Could not deallocate memory." // &
            TRIM(error)
       DEALLOCATE(indx1, stat=err)
       DEALLOCATE(indx2, stat=err)
       RETURN       
    END IF

  END SUBROUTINE matrix_print_r8





  SUBROUTINE matrix_print_r16(A, lu, error, frmt)

    !! Writes the given matrix to the given i/o-unit using a nice
    !! layout.
    !!
    !! Input:  - Square matrix           A(n,n)
    !!         - I/O-unit                lu

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:,:), INTENT(in)  :: A
    INTEGER, INTENT(in)                          :: lu
    CHARACTER(len=*), INTENT(inout)       :: error
    CHARACTER(len=*), INTENT(in), OPTIONAL       :: frmt
    REAL(rprec16), DIMENSION(:,:), ALLOCATABLE :: AA
    CHARACTER(len=30)                            :: form, c
    REAL                                         :: x, y
    INTEGER, DIMENSION(2)                        :: lo_bound, up_bound
    INTEGER                                      :: i, n, p, err

    INTRINSIC lbound, ubound, fraction, nint, char, trim 

    lo_bound = LBOUND(A)
    up_bound = UBOUND(A)
    ALLOCATE(AA(lo_bound(1):up_bound(1),lo_bound(2):up_bound(2)), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : matrix_print : Could not allocate memory." // &
            TRIM(error)
       RETURN       
    END IF
    AA = A

    ! Change integer n to character c (= n):
    n = up_bound(2) - lo_bound(2) + 1
    c = ''
    DO WHILE (n > 0)
       x = n / 10.0
       y = FRACTION(x)
       n = NINT(x - y)
       p = NINT(10 * y)
       c = CHAR(48 + p) // TRIM(c)
    END DO

    IF (PRESENT(frmt)) THEN
       form = '(' // TRIM(c) // TRIM(frmt) // ')'
    ELSE
       form = '(' // TRIM(c) // TRIM(frmt_default) // ')'
    END IF

    ! Write the matrix row-wise:
    DO i=lo_bound(1), up_bound(1)
       WRITE(lu,TRIM(form)) AA(i,:)
    END DO

    DEALLOCATE(AA, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : matrix_print : Could not deallocate memory." // &
            TRIM(error)
       RETURN       
    END IF

  END SUBROUTINE matrix_print_r16





  !! *Description*:
  !!
  !! Ref. Numerical Recipes for Fortran 90, Sect. 23.5 
  !!
  FUNCTION outer_product_r4(a, b)

    IMPLICIT NONE
    REAL(rprec4), DIMENSION(:), INTENT(in) :: a, b
    REAL(rprec4), DIMENSION(SIZE(a),SIZE(b)) :: outer_product_r4

    outer_product_r4 = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
         SPREAD(b,dim=1,ncopies=SIZE(a))

  END FUNCTION outer_product_r4





  !! *Description*:
  !!
  !! Ref. Numerical Recipes for Fortran 90, Sect. 23.5 
  !!
  FUNCTION outer_product_r8(a, b)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: a, b
    REAL(rprec8), DIMENSION(SIZE(a),SIZE(b)) :: outer_product_r8

    outer_product_r8 = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
         SPREAD(b,dim=1,ncopies=SIZE(a))

  END FUNCTION outer_product_r8





  !! *Description*:
  !!
  !! Ref. Numerical Recipes for Fortran 90, Sect. 23.5 
  !!
  FUNCTION outer_product_r16(a, b)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:), INTENT(in) :: a, b
    REAL(rprec16), DIMENSION(SIZE(a),SIZE(b)) :: outer_product_r16

    outer_product_r16 = SPREAD(a,dim=2,ncopies=SIZE(b)) * &
         SPREAD(b,dim=1,ncopies=SIZE(a))

  END FUNCTION outer_product_r16





  !! *Description*:
  !!
  !! Returns the triple product (x dot (y cross z)) of three vectors
  !! with three elements.
  !!
  REAL(rprec8) FUNCTION triple_product_r8(x,y,z)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(3), INTENT(in) :: x, y, z

    triple_product_r8 = DOT_PRODUCT(x,cross_product(y,z))

  END FUNCTION triple_product_r8





  !! *Description*:
  !!
  !! Returns cross product of two 3-vectors.
  !!
  FUNCTION cross_product_r8(x,y) RESULT(cross_prod)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(3), INTENT(in) :: x, y 
    REAL(rprec8), DIMENSION(3)             :: cross_prod

    cross_prod(1) = x(2)*y(3) - x(3)*y(2)
    cross_prod(2) = x(3)*y(1) - x(1)*y(3)
    cross_prod(3) = x(1)*y(2) - x(2)*y(1)

  END FUNCTION cross_product_r8





  SUBROUTINE vector_print(b, lu, column, error, frmt)

    !! Writes the given vector to the given i/o-unit using a nice
    !! layout.
    !!
    !! Input:  - Vector                          b(n)
    !!         - Column = .true., row = .false.  column
    !!         - I/O-unit                        lu

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)  :: b
    INTEGER, INTENT(in)                        :: lu
    LOGICAL, INTENT(in)                        :: column
    CHARACTER(len=*), INTENT(inout)     :: error
    CHARACTER(len=*), INTENT(in), OPTIONAL     :: frmt
    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: bb
    CHARACTER(len=30)                          :: form, c
    REAL                                       :: x, y
    INTEGER, DIMENSION(1)                      :: lo_bound, up_bound
    INTEGER                                    :: i, n, p, indx, err

    INTRINSIC lbound, ubound, fraction, nint, char, trim 

    lo_bound = LBOUND(b)
    up_bound = UBOUND(b)
    ALLOCATE(bb(lo_bound(1):up_bound(1)), stat=err)
    IF (err /= 0) THEN
       error = " -> linal : vector_print : Could not allocate memory." // &
            TRIM(error)
       RETURN
    END IF
    bb = b

    WHERE (ABS(bb) < EPSILON(bb)) bb = 0

    IF (column) THEN
       ! Writes a column vector:
       IF (PRESENT(frmt)) THEN
          ! Remove empty spaces:
          indx = INDEX(frmt,',',back=.TRUE.)
          form = '(' // frmt(1:indx-1) // '))'
       ELSE
          indx = INDEX(frmt_default,',',back=.TRUE.)
          form = '(' // frmt_default(1:indx-1) // '))'
       END IF
       DO i=lo_bound(1), up_bound(1)
          WRITE(lu,TRIM(form)) bb(i)
       END DO
    ELSE
       ! Writes a row vector:
       ! Change integer n to character c (= n):
       n = up_bound(1) - lo_bound(1) + 1
       c = ''
       DO WHILE (n > 0)
          x = n / 10.0
          y = FRACTION(x)
          n = NINT(x - y)
          p = NINT(10 * y)
          c = CHAR(48 + p) // TRIM(c)
       END DO
       IF (PRESENT(frmt)) THEN
          form = '(' // TRIM(c) // TRIM(frmt) // ')'
       ELSE
          form = '(' // TRIM(c) // TRIM(frmt_default) // ')'
       END IF
       WRITE(lu,TRIM(form)) bb
    END IF

    DEALLOCATE(bb, stat=err)
    IF (err /= 0) THEN
       error = " -> linal : vector_print : Could not deallocate memory." // &
            TRIM(error)
       RETURN
    END IF

  END SUBROUTINE vector_print



END MODULE linal
