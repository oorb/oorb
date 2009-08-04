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
!! Contains generic routines for parameter estimation.
!!
!! @author  MG
!! @version 2009-08-03
!!
MODULE estimators

  USE parameters
  USE linal
  IMPLICIT NONE

  INTERFACE leastSquares
     MODULE PROCEDURE leastSquares_r8_matrix
     MODULE PROCEDURE leastSquares_r8_blockdiag
     MODULE PROCEDURE leastSquares_r16_matrix
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
    LOGICAL, INTENT(inout)                   :: error
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
       error = .TRUE.
       WRITE(0,*) "fitdata: leastSquares: Shapes of arrays do not conform."
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
    IF (error) THEN
       WRITE(0,*) "Could not find inverse of inverse covariance matrix."
       RETURN
    END IF
    CALL matrix_print(MATMUL(cov_mat_param(:,:),inform_mat_param(:,:)),0)

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
    LOGICAL, INTENT(inout)                   :: error
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
       error = .TRUE.
       WRITE(0,*) "fitdata: leastSquares: Shapes of arrays do not conform."
       RETURN
    END IF

    ALLOCATE(inform_mat_param(nparam,nparam), &
         data_minus_computed(ndata,nmultidata), &
         param_corrections(nparam), &
         tmp(nparam), &
         d(nparam), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
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
    IF (error) THEN
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
       error = .TRUE.
       DEALLOCATE(inform_mat_param, stat=err)
       DEALLOCATE(data_minus_computed, stat=err)
       DEALLOCATE(param_corrections, stat=err)
       DEALLOCATE(tmp, stat=err)
       DEALLOCATE(d, stat=err)
       RETURN
    END IF

  END SUBROUTINE leastSquares_r8_blockdiag





  SUBROUTINE leastSquares_r16_matrix(indata, inform_mat_indata, &
       mask_indata, computed, design_mat, correction_factor, mask_param, &
       cov_mat_param, parameters, error)

    IMPLICIT NONE
    REAL(rprec16), DIMENSION(:,:), INTENT(in)   :: indata ! incl. cos(dec_obs)
    REAL(rprec16), DIMENSION(:,:,:), INTENT(in) :: inform_mat_indata ! incl. cos(dec_obs)
    LOGICAL, DIMENSION(:,:), INTENT(in)   :: mask_indata
    REAL(rprec16), DIMENSION(:,:), INTENT(in)   :: computed ! incl. cos(dec_comp)
    REAL(rprec16), DIMENSION(:,:,:), INTENT(in) :: design_mat ! (type,data,param) incl. cos(dec_obs)
    REAL(rprec16), INTENT(in)                   :: correction_factor
    LOGICAL, DIMENSION(:), INTENT(in)        :: mask_param
    REAL(rprec16), DIMENSION(:,:), INTENT(out)  :: cov_mat_param
    REAL(rprec16), DIMENSION(:), INTENT(inout)  :: parameters
    LOGICAL, INTENT(inout)                   :: error
    REAL(rprec16), DIMENSION(SIZE(parameters),SIZE(inform_mat_indata,dim=2)) :: tmp_mat
    REAL(rprec16), DIMENSION(SIZE(cov_mat_param,dim=1),SIZE(cov_mat_param,dim=2)) :: &
         inform_mat_param, identity, test
    REAL(rprec16), DIMENSION(SIZE(indata,dim=1),SIZE(indata,dim=2)) :: data_minus_computed 
    REAL(rprec16), DIMENSION(SIZE(parameters)) :: param_corrections, tmp, d
    INTEGER :: i, j, ndata, nparam, nmultidata

    ndata = SIZE(indata,dim=1)
    nmultidata = SIZE(indata,dim=2)
    nparam = SIZE(parameters)

    ! Check shapes of matrices:
    IF (nmultidata /= SIZE(inform_mat_indata,dim=1) .OR. &
         ndata /= SIZE(inform_mat_indata,dim=2) .OR. &
         ndata /= SIZE(inform_mat_indata,dim=3) .OR. &
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
       error = .TRUE.
       RETURN
    END IF

    ! (1) Compute the parameters' covariance matrix COV_MAT_PARAM and
    !     its inverse INFORM_MAT_PARAM.

    ! Sigma_param^(-1) = A^T Sigma_data^(-1) A:
    inform_mat_param(:,:) = 0.0_rprec16
    DO i=1, nmultidata
       tmp_mat = MATMUL(TRANSPOSE(design_mat(i,:,:)),inform_mat_indata(i,:,:))
       inform_mat_param(:,:) = inform_mat_param(:,:) + &
            MATMUL(tmp_mat(:,:),design_mat(i,:,:))
    END DO
    DO i=1,nparam
       ! Check which parameters are not used:
       IF (.NOT.mask_param(i)) THEN
          ! Nullify the corresponding row..
          inform_mat_param(i,:) = 0.0_rprec16
          ! ..and column,
          inform_mat_param(:,i) = 0.0_rprec16
          ! and set the diagonal element to 1.
          inform_mat_param(i,i) = 1.0_rprec16
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
    cov_mat_param(:,:) = matinv(inform_mat_param(:,:), error)
    IF (error) THEN
       WRITE(0,*) "Could not find inverse of inverse covariance matrix."
       RETURN
    END IF

    ! Enforce symmetry:
    DO i=1,nparam
       DO j=i,nparam
          cov_mat_param(i,j) = 0.5_rprec16*(cov_mat_param(i,j) + cov_mat_param(j,i)) 
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

    ! Is AA^(-1)=I ?
    identity(:,:) = 0.0_rprec16
    FORALL(i=1:nparam)
       identity(i,i) = 1.0_rprec16
    END FORALL
    test(:,:) = MATMUL(inform_mat_param, cov_mat_param)
    IF (ANY(ABS(test(:,:) - identity(:,:)) > 100*EPSILON(test(:,:)))) THEN
       error = .TRUE.
       WRITE(0,*) 'Error: Identity criterion not fulfilled:'
       CALL matrix_print(test,0)
       RETURN
    END IF

    ! (2) Compute parameter corrections and make the scaled
    !     corrections. If COMPUTED and PARAMETERS are zero,
    !     and CORRECTION_FACTOR is one, this corresponds to the 
    !     basic least-squares fit.

    data_minus_computed(:,:) = indata(:,:) - computed(:,:)
    WHERE (mask_indata(:,:))
       data_minus_computed(:,:) = 0.0_rprec16
    END WHERE
    ! d = A^T Sigma^(-1) y:
    d(:) = 0.0_rprec16
    DO i=1,nmultidata
       tmp_mat = MATMUL(TRANSPOSE(design_mat(i,:,:)),inform_mat_indata(i,:,:))
       d(:) = d(:) + MATMUL(tmp_mat, data_minus_computed(:,i))
    END DO
    WHERE (.NOT.mask_param(:))
       d(:) = 0.0_rprec16
    END WHERE

    ! param = Sigma_param d:
    param_corrections(:) = MATMUL(cov_mat_param(:,:), d(:))
    WHERE (.NOT.mask_param(:))
       param_corrections(:) = 0.0_rprec16
    END WHERE
    parameters(:) = parameters(:) + correction_factor*param_corrections(:)

  END SUBROUTINE leastSquares_r16_matrix





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
  SUBROUTINE amoeba(p, y, ftol, func, iter, error)

    IMPLICIT NONE
    INTERFACE
       REAL(8) FUNCTION func(x, error)
         IMPLICIT NONE
         REAL(8), DIMENSION(:), INTENT(in) :: x
         LOGICAL, INTENT(inout) :: error
       END FUNCTION func
    END INTERFACE
    REAL(rprec8), DIMENSION(:,:), INTENT(inout) :: p
    REAL(rprec8), DIMENSION(:), INTENT(inout) :: y
    REAL(rprec8), INTENT(in) :: ftol
    !real(rprec8), external :: func
    INTEGER(iprec4), INTENT(out) :: iter
    LOGICAL, INTENT(inout) :: error


    INTEGER(iprec4), PARAMETER :: ITMAX = 5000
    REAL(rprec8), PARAMETER :: TINY = 1.0e-10_rprec8

    INTEGER(iprec4) :: ihi, ndim ! Global variables.  
    REAL(rprec8), DIMENSION(SIZE(p,2)) :: psum 

    CALL amoeba_private 

  CONTAINS 

    SUBROUTINE amoeba_private 

      IMPLICIT NONE
      INTEGER(iprec4) :: i, ilo, inhi 
      REAL(rprec8) :: rtol, ysave, ytry, ytmp

      ndim = SIZE(p,dim=2)
      IF (ndim /= SIZE(p,dim=1) - 1 .OR. ndim /= SIZE(y) - 1) THEN
         error = .TRUE.
         RETURN
      END IF
      iter = 0 
      psum(:) = SUM(p(:,:),dim=1) 
      DO !Iteration loop.  
         ! Determine which point is the highest (worst),
         ! next-highest, and lowest (best)
         ilo = iminloc(y(:))
         ihi = imaxloc(y(:))
         ytmp = y(ihi) 
         y(ihi) = y(ilo)
         inhi = imaxloc(y(:))
         y(ihi) = ytmp
         rtol = 2.0_rprec8*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo))+TINY)
         ! Compute the fractional range from highest to lowest and
         ! return if satisfactory.  
         IF (rtol < ftol) THEN 
            ! If returning, put best point and value in slot 1. 
            CALL swap(y(1),y(ilo))
            CALL swap(p(1,:),p(ilo,:)) 
            RETURN 
         END IF
         IF (iter >= ITMAX) THEN
            ! TMAX exceeded in amoeba
            error = .TRUE.
            RETURN
         END IF
         ! Begin a new iteration. First extrapolate by a factor -1
         ! through the face of the simplex across from the high
         ! point, i.e., reflect the simplex from the high point.
         ytry = amotry(-1.0_rprec8) 
         iter = iter + 1 
         IF (ytry <= y(ilo)) THEN
            ! Gives a result better than the best point, so try an
            ! additional extrapolation by a factor of 2.  
            ytry=amotry(2.0_rprec8)
            iter = iter + 1 
         ELSE IF (ytry >= y(inhi)) THEN 
            ! The reflected point is worse than the second highest,
            ! so look for an intermediate lower point, i.e., do a
            ! one-dimensional contraction.
            ysave = y(ihi) 
            ytry = amotry(0.5_rprec8)
            iter = iter + 1
            IF (ytry >= ysave) THEN
               ! Can't seem to get rid of that high point. Better
               ! contract around the lowest (best) point.
               p(:,:) = 0.5_rprec8*(p(:,:)+SPREAD(p(ilo,:),1,SIZE(p,1)))
               DO i=1,ndim+1
                  IF (i /= ilo) y(i) = func(p(i,:), error) 
               END DO
               iter=iter+ndim ! Keep track of function evaluations.
               psum(:) = SUM(p(:,:),dim=1)
            END IF
         END IF
      END DO ! Go back for the test of doneness and the next iteration.

    END SUBROUTINE amoeba_private



    FUNCTION amotry(fac)

      IMPLICIT NONE
      REAL(rprec8), INTENT(IN) :: fac
      REAL(rprec8) :: amotry
      ! Extrapolates by a factor fac through the face of the
      ! simplex across from the high point, tries it, and replaces
      ! the high point if the new point is better.
      REAL(rprec8) :: fac1, fac2, ytry
      REAL(rprec8), DIMENSION(SIZE(p,2)) :: ptry 

      fac1 = (1.0_rprec8-fac)/ndim
      fac2 = fac1 - fac 
      ptry(:) = psum(:)*fac1-p(ihi,:)*fac2 
      ! Evaluate the function at the trial point.
      ytry = func(ptry, error)
      IF (ytry < y(ihi)) THEN 
         ! If it's better than the highest, then replace
         ! the highest.
         y(ihi) = ytry
         psum(:) = psum(:) - p(ihi,:) + ptry(:)
         p(ihi,:) = ptry(:)
      END IF
      amotry = ytry 

    END FUNCTION amotry

  END SUBROUTINE amoeba




END MODULE estimators
