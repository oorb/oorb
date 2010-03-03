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
!! *Module*description*:
!!
!! Contains generic routines for parameter estimation.
!!
!! @author  MG
!! @version 2010-03-02
!!
MODULE estimators

  USE parameters
  USE linal
  USE statistics
  IMPLICIT NONE

  PRIVATE :: LevenbergMarquardt_r8_blockdiag

  INTERFACE leastSquares
     MODULE PROCEDURE leastSquares_r8_matrix
     MODULE PROCEDURE leastSquares_r8_blockdiag
  END INTERFACE

  INTERFACE LevenbergMarquardt
     MODULE PROCEDURE LevenbergMarquardt_r8_blockdiag
  END INTERFACE

CONTAINS



  SUBROUTINE leastSquares_r8_matrix(indata, inform_mat_indata, &
       mask_indata, computed, design_mat, correction_factor, mask_param, &
       cov_mat_param, parameters, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: indata ! incl. cos(dec_obs)
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: inform_mat_indata ! incl. cos(dec_obs)
    LOGICAL, DIMENSION(:,:), INTENT(in)   :: mask_indata
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: computed ! incl. cos(dec_comp)
    REAL(rprec8), DIMENSION(:,:,:), INTENT(in) :: design_mat ! (type,data,param) incl. cos(dec_obs)
    REAL(rprec8), INTENT(in)                   :: correction_factor
    LOGICAL, DIMENSION(:), INTENT(in)        :: mask_param
    REAL(rprec8), DIMENSION(:,:), INTENT(out)  :: cov_mat_param
    REAL(rprec8), DIMENSION(:), INTENT(inout)  :: parameters
    CHARACTER(len=*), INTENT(inout)                   :: error
    REAL(rprec8), DIMENSION(SIZE(cov_mat_param,dim=1),SIZE(cov_mat_param,dim=2)) :: &
         inform_mat_param!, identity, test, matrix
    REAL(rprec8), DIMENSION(SIZE(indata,dim=1),SIZE(indata,dim=2)) :: data_minus_computed 
    REAL(rprec8), DIMENSION(SIZE(parameters)) :: param_corrections, tmp, d
    INTEGER :: i, j, ndata, nparam, nmultidata

    ndata = SIZE(indata,dim=1)
    nmultidata = SIZE(indata,dim=2)
    nparam = SIZE(parameters) 

    ! Check shapes of matrices:
    IF (nmultidata*ndata /= SIZE(inform_mat_indata,dim=1) .OR. &
         nmultidata*ndata /= SIZE(inform_mat_indata,dim=2) .OR. &
         nmultidata /= SIZE(design_mat,dim=1) .OR. &
         ndata /= SIZE(design_mat,dim=2) .OR. &
         nparam /= SIZE(design_mat,dim=3) .OR. &
         ndata < SIZE(mask_indata,dim=1) .OR. &
         nmultidata < SIZE(mask_indata,dim=2) .OR. &
         ndata /= SIZE(computed,dim=1) .OR. &
         nmultidata /= SIZE(computed,dim=2) .OR. &
         nparam < SIZE(mask_param) .OR. &
         nparam /= SIZE(cov_mat_param,dim=1) .OR. &
         nparam /= SIZE(cov_mat_param,dim=2)) THEN
       error = " -> estimators : leastSquares : Shapes of arrays do not conform." // &
            TRIM(error)
       RETURN
    END IF

    ! (1) Compute the parameters' covariance matrix COV_MAT_PARAM and
    !     its inverse INFORM_MAT_PARAM.

    ! Sigma_param^(-1) = A^T Sigma_data^(-1) A:
    inform_mat_param(:,:) = 0.0_rprec8
    DO i=1,nmultidata
       inform_mat_param(:,:) = inform_mat_param(:,:) + &
            MATMUL(MATMUL(TRANSPOSE(design_mat(i,:,:)), &
            inform_mat_indata(i:nmultidata*ndata:nmultidata, &
            i:nmultidata*ndata:nmultidata)), design_mat(i,:,:))
    END DO

    ! Check which parameters are not used:
    DO i=1,nparam
       IF (.NOT.mask_param(i)) THEN
          ! Nullify the corresponding row..
          inform_mat_param(i,:) = 0.0_rprec8
          ! ..and column,
          inform_mat_param(:,i) = 0.0_rprec8
          ! and set the diagonal element to 1.
          inform_mat_param(i,i) = 1.0_rprec8
       END IF
    END DO

    ! Renormalize to gain accuracy:
    DO i=1,nparam
       tmp(i) = SQRT(inform_mat_param(i,i)) 
    END DO
    DO i=1,nparam
       DO j=i,nparam
          inform_mat_param(i,j) = inform_mat_param(i,j)/(tmp(i)*tmp(j))
          ! Due to symmetry:
          inform_mat_param(j,i) = inform_mat_param(i,j)
       END DO
    END DO

    ! Sigma_param = (Sigma_param^(-1))^(-1) 
    cov_mat_param(:,:) = matinv(inform_mat_param(:,:), error, method="Cholesky")
    IF (LEN_TRIM(error) /= 0) THEN
       WRITE(0,*) "Could not find inverse of inverse covariance matrix."
       RETURN
    END IF
    !CALL matrix_print(MATMUL(cov_mat_param(:,:),inform_mat_param(:,:)),0)

    ! Enforce symmetry:
    DO i=1,nparam
       DO j=i,nparam
          cov_mat_param(i,j) = 0.5_rprec8*(cov_mat_param(i,j) + cov_mat_param(j,i)) 
          cov_mat_param(j,i) = cov_mat_param(i,j)
       END DO
    END DO

    ! Renormalize back:
    DO i=1,nparam
       DO j=i,nparam
          cov_mat_param(i,j) = cov_mat_param(i,j)/(tmp(i)*tmp(j))
          ! Due to symmetry:
          cov_mat_param(j,i) = cov_mat_param(i,j)
          inform_mat_param(i,j) = inform_mat_param(i,j)*(tmp(i)*tmp(j))
          ! Due to symmetry
          inform_mat_param(j,i) = inform_mat_param(i,j)
       END DO
    END DO

!!$    ! Is AA^(-1)=I ?
!!$    identity(:,:) = 0.0_rprec8
!!$    FORALL(i=1:nparam)
!!$       identity(i,i) = 1.0_rprec8
!!$    END FORALL
!!$    test(:,:) = MATMUL(inform_mat_param, cov_mat_param)
!!$    IF (ANY(ABS(test(:,:) - identity(:,:)) > 10*EPSILON(test(:,:)))) THEN
!!$       error = .TRUE.
!!$       WRITE(0,*) 'Error: Identity criterion not fulfilled:'
!!$       call matrix_print(test,0)
!!$       RETURN
!!$    END IF

    ! (2) Compute parameter corrections and make the scaled
    !     corrections. If COMPUTED and PARAMETERS are zero,
    !     and CORRECTION_FACTOR is one, this corresponds to the 
    !     linear least-squares fit.

    data_minus_computed(:,:) = indata(:,:) - computed(:,:)
    WHERE (.NOT.mask_indata(:,:))
       data_minus_computed(:,:) = 0.0_rprec8
    END WHERE

    ! d = A^T Sigma^(-1) y:
    d(:) = 0.0_rprec8
    DO i=1,nmultidata
       d(:) = d(:) + MATMUL(MATMUL(TRANSPOSE(design_mat(i,:,:)), &
            inform_mat_indata(i:nmultidata*ndata:nmultidata, &
            i:nmultidata*ndata:nmultidata)), data_minus_computed(:,i))
    END DO
    WHERE (.NOT.mask_param(:))
       d(:) = 0.0_rprec8
    END WHERE

    ! param = Sigma_param d:
    param_corrections(:) = MATMUL(cov_mat_param(:,:), d(:))
    WHERE (.NOT.mask_param(:))
       param_corrections(:) = 0.0_rprec8
    END WHERE
    parameters(:) = parameters(:) + correction_factor*param_corrections(:)

  END SUBROUTINE leastSquares_r8_matrix





  SUBROUTINE leastSquares_r8_blockdiag(indata, block_diag_inform_mat_indata, &
       mask_indata, computed, design_mat, correction_factor, mask_param, &
       cov_mat_param, parameters, error)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: indata ! incl. cos(dec_obs)
    REAL(rprec8), DIMENSION(:,:,:), INTENT(in) :: block_diag_inform_mat_indata ! incl. cos(dec_obs)
    LOGICAL, DIMENSION(:,:), INTENT(in)   :: mask_indata
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: computed ! incl. cos(dec_comp)
    REAL(rprec8), DIMENSION(:,:,:), INTENT(in) :: design_mat ! (type,data,param) incl. cos(dec_obs)
    REAL(rprec8), INTENT(in)                   :: correction_factor
    LOGICAL, DIMENSION(:), INTENT(in)        :: mask_param
    REAL(rprec8), DIMENSION(:,:), INTENT(out)  :: cov_mat_param
    REAL(rprec8), DIMENSION(:), INTENT(inout)  :: parameters
    CHARACTER(len=*), INTENT(inout)                   :: error
    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: &
         inform_mat_param!, identity, test, matrix
    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: data_minus_computed 
    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: param_corrections, tmp, d
    INTEGER :: i, j, ndata, nparam, nmultidata, err

    ndata = SIZE(indata,dim=1)
    nmultidata = SIZE(indata,dim=2)
    nparam = SIZE(parameters)

    ! Check shapes of matrices:
    IF (ndata /= SIZE(block_diag_inform_mat_indata,dim=1) .OR. &
         nmultidata /= SIZE(block_diag_inform_mat_indata,dim=2) .OR. &
         nmultidata /= SIZE(block_diag_inform_mat_indata,dim=3) .OR. &
         nmultidata /= SIZE(design_mat,dim=1) .OR. &
         ndata /= SIZE(design_mat,dim=2) .OR. &
         nparam /= SIZE(design_mat,dim=3) .OR. &
         ndata < SIZE(mask_indata,dim=1) .OR. &
         nmultidata < SIZE(mask_indata,dim=2) .OR. &
         ndata /= SIZE(computed,dim=1) .OR. &
         nmultidata /= SIZE(computed,dim=2) .OR. &
         nparam < SIZE(mask_param) .OR. &
         nparam /= SIZE(cov_mat_param,dim=1) .OR. &
         nparam /= SIZE(cov_mat_param,dim=2)) THEN
       error = " -> estimators : leastSquares : Shapes of arrays do not conform." // &
            TRIM(error)
       RETURN
    END IF

    ALLOCATE(inform_mat_param(nparam,nparam), &
         data_minus_computed(ndata,nmultidata), &
         param_corrections(nparam), &
         tmp(nparam), &
         d(nparam), &
         stat=err)
    IF (err /= 0) THEN
       error = " -> estimators : leastSquares : Could not allocate memory." // &
            TRIM(error)
       DEALLOCATE(inform_mat_param, stat=err)
       DEALLOCATE(data_minus_computed, stat=err)
       DEALLOCATE(param_corrections, stat=err)
       DEALLOCATE(tmp, stat=err)
       DEALLOCATE(d, stat=err)
       RETURN
    END IF

    ! (1) Compute the parameters' covariance matrix COV_MAT_PARAM and
    !     its inverse INFORM_MAT_PARAM.

    ! Sigma_param^(-1) = A^T Sigma_data^(-1) A:
    inform_mat_param(:,:) = 0.0_rprec8
    DO i=1,ndata
       inform_mat_param(:,:) = inform_mat_param(:,:) + &
            MATMUL(MATMUL(TRANSPOSE(design_mat(1:nmultidata,i,1:nparam)), &
            block_diag_inform_mat_indata(i,1:nmultidata,1:nmultidata)), &
            design_mat(1:nmultidata,i,1:nparam))
    END DO

    ! Check which parameters are not used:
    DO i=1,nparam
       IF (.NOT.mask_param(i)) THEN
          ! Nullify the corresponding row..
          inform_mat_param(i,:) = 0.0_rprec8
          ! ..and column,
          inform_mat_param(:,i) = 0.0_rprec8
          ! and set the diagonal element to 1.
          inform_mat_param(i,i) = 1.0_rprec8
       END IF
    END DO

    ! Renormalize to gain accuracy:
    DO i=1,nparam
       tmp(i) = SQRT(inform_mat_param(i,i)) 
    END DO
    DO i=1,nparam
       DO j=i,nparam
          inform_mat_param(i,j) = inform_mat_param(i,j)/(tmp(i)*tmp(j))
          ! Due to symmetry:
          inform_mat_param(j,i) = inform_mat_param(i,j)
       END DO
    END DO

    ! Sigma_param = (Sigma_param^(-1))^(-1) 
    cov_mat_param(:,:) = matinv(inform_mat_param(:,:), error, method="Cholesky")
    IF (LEN_TRIM(error) /= 0) THEN
       WRITE(0,*) "Could not find inverse of inverse covariance matrix."
       DEALLOCATE(inform_mat_param, stat=err)
       DEALLOCATE(data_minus_computed, stat=err)
       DEALLOCATE(param_corrections, stat=err)
       DEALLOCATE(tmp, stat=err)
       DEALLOCATE(d, stat=err)
       RETURN
    END IF

    ! Enforce symmetry:
    DO i=1,nparam
       DO j=i,nparam
          cov_mat_param(i,j) = 0.5_rprec8*(cov_mat_param(i,j) + cov_mat_param(j,i)) 
          cov_mat_param(j,i) = cov_mat_param(i,j)
       END DO
    END DO

    ! Renormalize back:
    DO i=1,nparam
       DO j=i,nparam
          cov_mat_param(i,j) = cov_mat_param(i,j)/(tmp(i)*tmp(j))
          ! Due to symmetry:
          cov_mat_param(j,i) = cov_mat_param(i,j)
          inform_mat_param(i,j) = inform_mat_param(i,j)*(tmp(i)*tmp(j))
          ! Due to symmetry
          inform_mat_param(j,i) = inform_mat_param(i,j)
       END DO
    END DO

!!$    ! Is AA^(-1)=I ?
!!$    identity(:,:) = 0.0_rprec8
!!$    FORALL(i=1:nparam)
!!$       identity(i,i) = 1.0_rprec8
!!$    END FORALL
!!$    test(:,:) = MATMUL(inform_mat_param, cov_mat_param)
!!$    IF (ANY(ABS(test(:,:) - identity(:,:)) > 10*EPSILON(test(:,:)))) THEN
!!$       error = .TRUE.
!!$       WRITE(0,*) 'Error: Identity criterion not fulfilled:'
!!$       call matrix_print(test,0)
!!$       RETURN
!!$    END IF

    ! (2) Compute parameter corrections and make the scaled
    !     corrections. If COMPUTED and PARAMETERS are zero,
    !     and CORRECTION_FACTOR is one, this corresponds to the 
    !     linear least-squares fit.

    data_minus_computed(:,:) = indata(:,:) - computed(:,:)
    WHERE (.NOT.mask_indata(:,:))
       data_minus_computed(:,:) = 0.0_rprec8
    END WHERE

    ! d = A^T Sigma^(-1) y:
    d(:) = 0.0_rprec8
    DO i=1,ndata
       d(:) = d(:) + MATMUL(MATMUL(TRANSPOSE(design_mat(1:nmultidata,i,1:nparam)), &
            block_diag_inform_mat_indata(i,1:nmultidata,1:nmultidata)), &
            data_minus_computed(i,1:nmultidata))
    END DO
    WHERE (.NOT.mask_param(:))
       d(:) = 0.0_rprec8
    END WHERE

    ! param = Sigma_param d:
    param_corrections(:) = MATMUL(cov_mat_param(:,:), d(:))
    WHERE (.NOT.mask_param(:))
       param_corrections(:) = 0.0_rprec8
    END WHERE
    parameters(:) = parameters(:) + correction_factor*param_corrections(:)

    DEALLOCATE(inform_mat_param, data_minus_computed, &
         param_corrections, tmp, d, stat=err)
    IF (err /= 0) THEN
       error = " -> estimators : leastSquares : Could not deallocate memory." // &
            TRIM(error)
       DEALLOCATE(inform_mat_param, stat=err)
       DEALLOCATE(data_minus_computed, stat=err)
       DEALLOCATE(param_corrections, stat=err)
       DEALLOCATE(tmp, stat=err)
       DEALLOCATE(d, stat=err)
       RETURN
    END IF

  END SUBROUTINE leastSquares_r8_blockdiag





  SUBROUTINE LevenbergMarquardt_r8_blockdiag(x, measur, &
       information_matrix_measur, mask_measur, sub, &
       niter, chi2_frac, cov_mat_param, params, &
       residuals, chi2, errstr)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: x
    REAL(rprec8), DIMENSION(:,:), INTENT(in)   :: measur ! incl. cos(dec_obs)
    REAL(rprec8), DIMENSION(:,:,:), INTENT(in) :: information_matrix_measur ! incl. cos(dec_obs)
    LOGICAL, DIMENSION(:,:), INTENT(in)        :: mask_measur
    INTERFACE
       SUBROUTINE sub(x, measur, params, residuals, jacobians, chi2, errstr)
         USE parameters
         IMPLICIT NONE
         REAL(rprec8), DIMENSION(:,:), INTENT(in)    :: x
         REAL(rprec8), DIMENSION(:,:), INTENT(in)    :: measur         
         REAL(rprec8), DIMENSION(:), INTENT(inout)   :: params
         REAL(rprec8), DIMENSION(:,:), INTENT(out)   :: residuals
         REAL(rprec8), DIMENSION(:,:,:), INTENT(out) :: jacobians 
         REAL(rprec8), INTENT(out)                   :: chi2
         CHARACTER(len=*), INTENT(inout)             :: errstr
       END SUBROUTINE sub
    END INTERFACE
    REAL(rprec8), INTENT(in)                   :: chi2_frac
    INTEGER, INTENT(inout)                     :: niter
    REAL(rprec8), DIMENSION(:,:), INTENT(out)  :: cov_mat_param
    REAL(rprec8), DIMENSION(:), INTENT(inout)  :: params
    REAL(rprec8), DIMENSION(:,:), INTENT(out)  :: residuals
    REAL(rprec8), INTENT(out)                  :: chi2   
    CHARACTER(len=*), INTENT(inout)            :: errstr

    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: alpha
    REAL(rprec8) :: chi2_previous, lambda
    INTEGER :: i, ndata, nparam, nmultidata, nfit, err


    ndata = SIZE(residuals,dim=1)
    nmultidata = SIZE(residuals,dim=2)
    nparam = SIZE(params)

    ALLOCATE(alpha(nparam,nparam), stat=err)

    chi2_previous = HUGE(chi2_previous)
    lambda = -1.0_rprec8
    DO i=1,niter
       CALL LevenbergMarquardt_private
       IF (ABS(chi2 - chi2_previous)/chi2 < chi2_frac .AND. &
            chi2 < chi2_previous) THEN
          EXIT
       ELSE
          chi2_previous = chi2
       END IF
    END DO
    lambda = 0.0_rprec8
    CALL LevenbergMarquardt_private

    DEALLOCATE(alpha, stat=err)

  CONTAINS


    SUBROUTINE LevenbergMarquardt_private

      REAL(rprec8), SAVE :: chi2_
      REAL(rprec8), DIMENSION(:), ALLOCATABLE, SAVE :: params_, beta
      REAL(rprec8), DIMENSION(:,:), ALLOCATABLE, SAVE :: param_corrections

      IF (lambda < 0.0_rprec8) THEN
         ALLOCATE(params_(nparam),beta(nparam),param_corrections(nparam,1))
         lambda = 0.001_rprec8
         params_ = params
         CALL coefficients(params_, alpha, beta)
         chi2_ = chi2
      END IF
      cov_mat_param = alpha
      cov_mat_param = diagonal_multiplication(cov_mat_param,1.0_rprec8+lambda,errstr)
      IF (LEN_TRIM(errstr) /= 0) THEN
         errstr = " -> estimators : LevenbergMarquardt : LevenbergMarquardt_private : ." // &
              TRIM(errstr)
         RETURN
      END IF
      param_corrections(:,1) = beta
      CALL gauss_jordan(cov_mat_param, param_corrections, errstr)
      IF (LEN_TRIM(errstr) /= 0) THEN
         errstr = " -> estimators : LevenbergMarquardt : LevenbergMarquardt_private : ." // &
              TRIM(errstr)
         RETURN         
      END IF
      IF (lambda == 0.0_rprec8) THEN
         DEALLOCATE(params_, beta, param_corrections)
         RETURN
      END IF
      params_ = params + param_corrections(:,1)
      CALL coefficients(params_, cov_mat_param, param_corrections(:,1))
      IF (chi2 < chi2_) THEN
         lambda = 0.1_rprec8*lambda
         chi2_ = chi2
         alpha = cov_mat_param
         beta = param_corrections(:,1)
         params = params_
      ELSE
         lambda = 10.0_rprec8*lambda
         chi2 = chi2_
      END IF

    END SUBROUTINE LevenbergMarquardt_private


    SUBROUTINE coefficients(params, alpha, beta)

      IMPLICIT NONE
      REAL(rprec8), DIMENSION(:), INTENT(inout) :: params
      REAL(rprec8), DIMENSION(:), INTENT(out) :: beta
      REAL(rprec8), DIMENSION(:,:), INTENT(out) :: alpha

      REAL(rprec8), DIMENSION(nmultidata,ndata,nparam) :: jacobians
      REAL(rprec8), DIMENSION(nparam,nmultidata) :: tmp
      INTEGER :: i

      CALL sub(x, measur, params, residuals, jacobians, chi2, errstr)
      IF (LEN_TRIM(errstr) /= 0) THEN
         errstr = " -> estimators : LevenbergMarquardt : coefficients : ." // &
              TRIM(errstr)
         RETURN
      END IF
      ! Approximate Hessian by multiplying Jacobians
      ! alpha = cov_param^(-1) = J^T Sigma_obs^(-1) J:
      ! beta = J^T Sigma_obs^(-1) y:
      alpha = 0.0_rprec8
      beta = 0.0_rprec8
      DO i=1,ndata
         tmp = MATMUL(TRANSPOSE(jacobians(1:nmultidata,1:nparam,i)), &
              information_matrix_measur(i,1:nmultidata,1:nmultidata))
         alpha = alpha + MATMUL(tmp, jacobians(1:nmultidata,1:nparam,i))
         beta = beta + MATMUL(tmp, residuals(i,1:nmultidata))
      END DO

    END SUBROUTINE coefficients

  END SUBROUTINE LevenbergMarquardt_r8_blockdiag





  !! *Description*:
  !!
  !! Minimization of the function func in N dimensions by the
  !! downhill simplex method of Nelder and Mead. The (N + 1) × N
  !! matrix p is input. Its N + 1 rows are N-dimensional vectors that
  !! are the vertices of the starting simplex. Also input is the
  !! vector y of length N + 1, whose components must be preinitialized
  !! to the values of func evaluated at the N + 1 vertices (rows) of
  !! p and ftol the fractional convergence tolerance to be achieved in
  !! the function value (n.b.!). On output, p and y will have been
  !! reset to N+1 new points all within ftol of a minimum function
  !! value, and iter gives the number of function evaluations taken.
  !!
  !! Parameters: The maximum allowed number of function evaluations,
  !! and a small number.
  !!
  SUBROUTINE simplex(func, param, thresh, iter, p_matrix, y_vector, verb, errstr)

    IMPLICIT NONE
    INTERFACE
       REAL(8) FUNCTION func(x, param, errstr)
         IMPLICIT NONE
         REAL(8), DIMENSION(:), INTENT(in) :: x
         REAL(8), DIMENSION(:), INTENT(in) :: param
         CHARACTER(len=*), INTENT(inout) :: errstr
       END FUNCTION func
    END INTERFACE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: param
    REAL(rprec8), INTENT(in) :: thresh
    INTEGER, INTENT(inout) :: iter
    REAL(rprec8), DIMENSION(:,:), INTENT(inout) :: p_matrix
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: y_vector
    INTEGER(iprec4), INTENT(in) :: verb
    CHARACTER(len=*), INTENT(inout) :: errstr

    CHARACTER(len=8) :: str
    REAL(rprec8), DIMENSION(SIZE(p_matrix,dim=2)) :: psum
    REAL(rprec8) :: ysave, ytry, ytmp
    INTEGER(iprec4) :: i, ilo, ihi, inhi, ndim, iter_max
    LOGICAL :: error

    ndim = SIZE(p_matrix,dim=2)
    error = .FALSE.
    CALL toString(ndim-1, str, error)
    IF (error) THEN
       errstr = " -> estimators : simplex : Could not convert integer to character string." // &
            TRIM(errstr)
       RETURN       
    END IF
    IF (ndim /= SIZE(p_matrix,dim=1) - 1 .OR. ndim /= SIZE(y_vector) - 1) THEN
       errstr = " -> estimators : simplex : Input vectors are not compatible." // &
            TRIM(errstr)
       RETURN
    END IF
    iter_max = iter
    iter = 0
    psum(:) = SUM(p_matrix(:,:),dim=1) 
    DO !Iteration loop.  
       ! Determine which point is the highest (worst),
       ! next-highest, and lowest (best)
       ilo = iminloc(y_vector(:))
       ihi = imaxloc(y_vector(:))
       ytmp = y_vector(ihi) 
       y_vector(ihi) = y_vector(ilo)
       inhi = imaxloc(y_vector(:))
       y_vector(ihi) = ytmp
       IF (verb >= 2) THEN
          DO i=1,SIZE(y_vector,dim=1)
             IF (i == ilo) THEN
                WRITE(*,"('ILO:',1X,"//TRIM(str)//"(F20.15,1X))") p_matrix(i,:), y_vector(i) 
             ELSE IF (i == ihi) THEN
                WRITE(*,"('IHI:',1X,"//TRIM(str)//"(F20.15,1X))") p_matrix(i,:), y_vector(i) 
             ELSE
                WRITE(*,"(5X,"//TRIM(str)//"(F20.15,1X))") p_matrix(i,:), y_vector(i) 
             END IF
          END DO
          WRITE(*,"(1X)")
       END IF
       ! Compute the fractional range from highest to lowest and
       ! return if satisfactory.
       IF (y_vector(ilo) < thresh) THEN 
          ! If returning, put best point and value in slot 1. 
          CALL swap(y_vector(1),y_vector(ilo))
          CALL swap(p_matrix(1,:),p_matrix(ilo,:)) 
          RETURN 
       END IF
       IF (iter > iter_max) THEN
          errstr = " -> estimators : simplex_try : Maximum number of iterations exceeded." // &
               TRIM(errstr)
          RETURN
       END IF
       ! Begin a new iteration. First extrapolate by a factor -1
       ! through the face of the simplex across from the high
       ! point, i.e., reflect the simplex from the high point.
       ytry = simplex_try(func, param, -1.0_rprec8, ihi, psum, p_matrix, y_vector, errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          errstr = " -> estimators : simplex : ." // &
               TRIM(errstr)
          RETURN
       END IF
       iter = iter + 1 
       IF (ytry <= y_vector(ilo)) THEN
          ! Gives a result better than the best point, so try an
          ! additional extrapolation by a factor of 2.  
          ytry = simplex_try(func, param, 2.0_rprec8, ihi, psum, p_matrix, y_vector, errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             errstr = " -> estimators : simplex : ." // &
                  TRIM(errstr)
             RETURN
          END IF
          iter = iter + 1 
       ELSE IF (ytry >= y_vector(inhi)) THEN 
          ! The reflected point is worse than the second highest,
          ! so look for an intermediate lower point, i.e., do a
          ! one-dimensional contraction.
          ysave = y_vector(ihi) 
          ytry = simplex_try(func, param, 0.5_rprec8, ihi, psum, p_matrix, y_vector, errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             errstr = " -> estimators : simplex : ." // &
                  TRIM(errstr)
             RETURN
          END IF
          iter = iter + 1
          IF (ytry >= ysave) THEN
             ! Can't seem to get rid of that high point. Better
             ! contract around the lowest (best) point.
             p_matrix(:,:) = 0.5_rprec8*(p_matrix(:,:)+SPREAD(p_matrix(ilo,:),1,SIZE(p_matrix,1)))
             DO i=1,ndim+1
                IF (i /= ilo) THEN
                   y_vector(i) = func(p_matrix(i,:),param,errstr)
                   IF (LEN_TRIM(errstr) /= 0) THEN
                      errstr = " -> estimators : simplex : ." // &
                           TRIM(errstr)
                      RETURN
                   END IF
                END IF
             END DO
             iter = iter + ndim ! Keep track of function evaluations.
             psum(:) = SUM(p_matrix(:,:),dim=1)
          END IF
       END IF
    END DO ! Go back for the test of doneness and the next iteration.

  END SUBROUTINE simplex


  !! *Description*:
  !!
  !! Extrapolates by a factor fac through the face of the
  !! simplex across from the high point, tries it, and replaces
  !! the high point if the new point is better.
  !!
  REAL(rprec8) FUNCTION simplex_try(func, param, fac, ihi, psum, p_matrix, y_vector, errstr)

    IMPLICIT NONE
    INTERFACE
       REAL(8) FUNCTION func(x, param, errstr)
         IMPLICIT NONE
         REAL(8), DIMENSION(:), INTENT(in) :: x
         REAL(8), DIMENSION(:), INTENT(in) :: param
         CHARACTER(len=*), INTENT(inout) :: errstr
       END FUNCTION func
    END INTERFACE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: param
    REAL(rprec8), INTENT(IN) :: fac
    INTEGER, INTENT(in) :: ihi
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: psum
    REAL(rprec8), DIMENSION(:,:), INTENT(inout) :: p_matrix
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: y_vector
    CHARACTER(len=*), INTENT(inout) :: errstr

    REAL(rprec8), DIMENSION(SIZE(p_matrix,2)) :: ptry 
    REAL(rprec8) :: fac1, fac2, ytry

    fac1 = (1.0_rprec8-fac)/SIZE(p_matrix,dim=2)
    fac2 = fac1 - fac 
    ptry(:) = psum(:)*fac1-p_matrix(ihi,:)*fac2 
    ! Evaluate the function at the trial point.
    ytry = func(ptry,param,errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       errstr = " -> estimators : simplex_try : ." // &
            TRIM(errstr)
       RETURN
    END IF
    IF (ytry < y_vector(ihi)) THEN 
       ! If it's better than the highest, then replace
       ! the highest.
       y_vector(ihi) = ytry
       psum(:) = psum(:) - p_matrix(ihi,:) + ptry(:)
       p_matrix(ihi,:) = ptry(:)
    END IF
    simplex_try = ytry 

  END FUNCTION simplex_try


END MODULE estimators
