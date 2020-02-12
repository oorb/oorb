!====================================================================!
!                                                                    !
! Copyright 2002-2018,2019                                           !
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
!! Tools for statistics.
!!
!! @author  MG
!! @version 2019-10-29
!!
MODULE statistics

  USE parameters
  USE utilities
  USE sort
  USE linal

  IMPLICIT NONE

  PRIVATE :: chi_square_blockdiag
  PRIVATE :: moments_r8
  PRIVATE :: confidence_limits_hist_r8
  PRIVATE :: confidence_limits_sample_r8

  INTERFACE chi_square
     MODULE PROCEDURE chi_square_blockdiag
  END INTERFACE chi_square

  INTERFACE moments
     MODULE PROCEDURE moments_r8
  END INTERFACE moments

  INTERFACE confidence_limits
     MODULE PROCEDURE confidence_limits_hist_r8
     MODULE PROCEDURE confidence_limits_sample_r8
  END INTERFACE confidence_limits

  INTERFACE mahalanobis_distance
    MODULE PROCEDURE mahalanobis_distance_r8
  END INTERFACE mahalanobis_distance

CONTAINS




  !! *Description*:
  !!
  !! Tested.
  !!
  !! @author  MG
  !! @version 2009-10-13
  !!
  REAL(rprec8) FUNCTION chi_square_blockdiag(residuals, information_matrix, mask, errstr)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: residuals
    REAL(rprec8), DIMENSION(:,:,:), INTENT(in) :: information_matrix
    LOGICAL, DIMENSION(:,:), INTENT(in), OPTIONAL :: mask
    CHARACTER(len=*), INTENT(inout) :: errstr

    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: residuals_
    REAL(rprec8), DIMENSION(1) :: chi_square_blockdiag_
    INTEGER :: i, nobs, nmulti, err

    chi_square_blockdiag = 0.0_rprec8
    chi_square_blockdiag_ = 0.0_rprec8
    nobs = SIZE(residuals,dim=1)
    nmulti = SIZE(residuals,dim=2)

    IF (SIZE(information_matrix,dim=1) /= nobs .OR. &
         SIZE(information_matrix,dim=2) /= nmulti .OR. &
         SIZE(information_matrix,dim=3) /= nmulti) THEN
       errstr = " -> statistics : chi_square : Shape of input matrices do not conform." // &
            TRIM(errstr)
       RETURN
    END IF

    ALLOCATE(residuals_(nobs,nmulti), stat=err)
    IF (err /= 0) THEN
       errstr = " -> statistics : chi_square : Could not allocate memory." // &
            TRIM(errstr)
       RETURN
    END IF
    residuals_ = residuals
    IF (PRESENT(mask)) THEN
       WHERE (.NOT. mask)
          residuals_ = 0.0_rprec8
       END WHERE
    END IF
    DO i=1,nobs
       chi_square_blockdiag_ = chi_square_blockdiag_ + &
            MATMUL(MATMUL(residuals_(i,1:nmulti), &
            information_matrix(i,1:nmulti,1:nmulti)), &
            TRANSPOSE(residuals_(i:i,1:nmulti)))
    END DO
    chi_square_blockdiag = chi_square_blockdiag_(1)
    DEALLOCATE(residuals_, stat=err)
    IF (err /= 0) THEN
       WRITE(0,*) "chi_square" // &
            "Could not deallocate memory."
       RETURN
    END IF

  END FUNCTION chi_square_blockdiag





  !! *Description*:
  !!
  !! Computes a histogram (histo) for given data (indata) either as
  !! number density (without pdf) or from pdf (with pdf).
  !!
  SUBROUTINE histogram(indata, histo, xmin_in, xmax_in, pdf, xmax, dx)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)                      :: indata
    REAL(rprec8), DIMENSION(:,:), INTENT(out)                   :: histo
    REAL(rprec8), INTENT(in), OPTIONAL                          :: xmin_in
    REAL(rprec8), INTENT(in), OPTIONAL                          :: xmax_in
    REAL(rprec8), DIMENSION(SIZE(indata)), INTENT(in), OPTIONAL :: pdf
    REAL(rprec8), INTENT(out), OPTIONAL                         :: xmax
    REAL(rprec8), INTENT(out), OPTIONAL                         :: dx

    REAL(rprec8), DIMENSION(SIZE(indata)) :: pdf_
    REAL(rprec8) :: dx_, x_min, x_max
    INTEGER :: ngrid, nvalue, i, j

    ngrid = SIZE(histo,dim=1)
    nvalue = SIZE(indata)
    IF (PRESENT(xmin_in)) THEN
       x_min = xmin_in
    ELSE
       x_min = MINVAL(indata)
    END IF
    IF (PRESENT(xmax_in)) THEN
       x_max = xmax_in
    ELSE
       x_max = MAXVAL(indata)
    END IF
    dx_ = (x_max-x_min) / ngrid
    DO i=1, ngrid
       histo(i,1) = x_min + 0.5_rprec8*dx_ + dx_*(i-1)
    END DO

    IF (PRESENT(pdf)) THEN
       ! Use given pdf the sum of which is normalized to unity:
       pdf_ = pdf/SUM(pdf)
    ELSE
       ! Use number density:
       pdf_ = 1.0_rprec8
    END IF

    histo(:,2) = 0.0_rprec8
    DO i=1,nvalue
       DO j=1,ngrid
          IF (ABS(indata(i)-histo(j,1)) < 0.5_rprec8*dx_) THEN
             histo(j,2) = histo(j,2) + pdf_(i)
             EXIT
          END IF
       END DO
    END DO
    IF (PRESENT(xmax)) THEN
       xmax = histo(MAXLOC(histo(:,2),dim=1),1)
    END IF
    IF (PRESENT(dx)) THEN
       dx = dx_
    END IF

  END SUBROUTINE histogram





  !! *Description*:
  !!
  !! Computes various statistical quantities for a given data set.
  !! Current output: arithmetic mean, standard deviation, skewness,
  !! and kurtosis.
  !!
  SUBROUTINE moments_r8(indata, pdf, mask, mean, std_dev, skew, kurt, errstr)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)           :: indata
    REAL(rprec8), DIMENSION(:), INTENT(in), OPTIONAL :: pdf
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in)      :: mask
    REAL(rprec8), OPTIONAL, INTENT(out)              :: mean, std_dev, skew, kurt
    CHARACTER(len=*), INTENT(inout)                  :: errstr

    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: pdf_
    REAL(rprec8) :: mean_, std_dev_
    INTEGER :: ndata, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask_

    ndata = SIZE(indata)
    ALLOCATE(mask_(ndata), pdf_(ndata), stat=err)
    IF (err /= 0) THEN
       RETURN
    END IF
    mask_ = .TRUE.
    IF (PRESENT(mask)) THEN
       mask_ = mask
    END IF
    ndata = COUNT(mask_)
    IF (PRESENT(pdf)) THEN
       ! Make sure the distribution is normalized:
       pdf_ = pdf/SUM(pdf, mask=mask_)
    END IF

    ! Lowest moments:
    ! Mean:
    IF (PRESENT(pdf)) THEN
       mean_ = SUM(indata*pdf_, mask=mask_)
    ELSE
       mean_ = SUM(indata, mask=mask_)/ndata
    END IF
    IF (PRESENT(mean)) THEN
       mean = mean_
    END IF

    ! Standard deviation:
    IF (PRESENT(std_dev) .OR. PRESENT(skew) .OR. PRESENT(kurt)) THEN
       IF (PRESENT(pdf)) THEN
          std_dev_ = SQRT(SUM(((indata - &
               indata(MAXLOC(pdf_,dim=1)))**2.0_rprec8)*pdf_, &
               mask=mask_))
       ELSE
          std_dev_ = SQRT(SUM(((indata - mean_)**2.0_rprec8), &
               mask=mask_) / (ndata-1))
       END IF
       IF (PRESENT(std_dev)) THEN
          std_dev = std_dev_
       END IF
    END IF

    ! Skewness:
    IF (PRESENT(skew)) THEN
       IF (PRESENT(pdf)) THEN
          skew = SUM(((indata - indata(MAXLOC(pdf_,dim=1)))**3.0_rprec8)*pdf_, &
               mask=mask_) / std_dev_
       ELSE
          skew = SUM(((indata - mean_)**3.0_rprec8), &
               mask=mask_) / (std_dev_*ndata)
       END IF
    END IF

    ! Kurtosis:
    IF (PRESENT(kurt)) THEN
       IF (PRESENT(pdf)) THEN
          kurt = SUM(((indata - indata(MAXLOC(pdf_,dim=1)))**4.0_rprec8)*pdf_, &
               mask=mask_) / std_dev_
       ELSE
          kurt = SUM(((indata - mean_)**4.0_rprec8), &
               mask=mask_) / (std_dev_*ndata)
       END IF
    END IF

    !r = sum((data_set-mean_)*(data_set)

    DEALLOCATE(pdf_, mask_, stat=err)
    IF (err /= 0) THEN
       DEALLOCATE(pdf_, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF

  END SUBROUTINE moments_r8





  !! *Description*:
  !!
  !! This routine calculates the end points of the confidence interval
  !! (or credible interval) by incorporating bins of a histogram until
  !! the sum of normalized weights is larger than the requested
  !! probability mass.
  !!
  SUBROUTINE confidence_limits_hist_r8(indata, pdf, nhist, &
       mask, probability_mass, peak, bounds, errstr)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)            :: indata
    REAL(rprec8), DIMENSION(:), INTENT(in)            :: pdf
    INTEGER, INTENT(in)                               :: nhist
    REAL(rprec8), INTENT(in), OPTIONAL                :: probability_mass
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in)       :: mask
    REAL(rprec8), OPTIONAL, INTENT(out)               :: peak
    REAL(rprec8), DIMENSION(2), OPTIONAL, INTENT(out) :: bounds
    CHARACTER(len=*), INTENT(inout)                            :: errstr

    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: histo, histo_
    REAL(rprec8) :: probability_mass_
    INTEGER :: ndata, err, imax, ilo, ihi
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask_

    ndata = SIZE(indata)
    ALLOCATE(mask_(ndata), histo(nhist,2), histo_(nhist,2), stat=err)
    IF (err /= 0) THEN
       errstr = " -> statistics : confidence_limits : Could not allocate memory." // &
            TRIM(errstr)
       DEALLOCATE(histo, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF
    mask_ = .TRUE.
    IF (PRESENT(mask)) THEN
       mask_ = mask
    ELSE
       mask_ = .TRUE.
    END IF

    IF (PRESENT(bounds) .AND. .NOT.PRESENT(probability_mass)) THEN
       errstr = " -> statistics : confidence_limits : Probability mass not given." // &
            TRIM(errstr)
       DEALLOCATE(histo, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF

    ! Produce histogram based on pdf:
    CALL histogram(indata, histo, pdf=pdf)
    IF (PRESENT(peak)) THEN
       ! ML solution
       imax = imaxloc(histo(:,2))
       peak = histo(imax,1)
    END IF

    IF (PRESENT(bounds)) THEN
       histo(:,2) = histo(:,2)/SUM(histo(:,2))
       histo_ = histo
       bounds(1) = HUGE(bounds(1))
       bounds(2) = -HUGE(bounds(2))
       ilo = nhist
       ihi = 1
       probability_mass_ = 0.0_rprec8
       DO WHILE (probability_mass_ < probability_mass)
          imax = imaxloc(histo_(:,2))
          IF (imax < ilo) THEN
             ilo = imax
          END IF
          IF (imax > ihi) THEN
             ihi = imax
          END IF
          probability_mass_ = SUM(histo(ilo:ihi,2))
          IF (histo(imax,1) < bounds(1)) THEN
             bounds(1) = histo(imax,1)
          END IF
          IF (histo(imax,1) > bounds(2)) THEN
             bounds(2) = histo(imax,1)
          END IF
          histo_(imax,2) = 0.0_rprec8
       END DO
    END IF

    DEALLOCATE(histo, histo_, mask_, stat=err)
    IF (err /= 0) THEN
       errstr = " -> statistics : confidence_limits : Could not deallocate memory." // &
            TRIM(errstr)
       DEALLOCATE(histo, stat=err)
       DEALLOCATE(histo_, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF

  END SUBROUTINE confidence_limits_hist_r8





  !! *Description*:
  !!
  !! This routine calculates the end points of the confidence interval
  !! (or credible interval) by descending from the ML solution until
  !! the sum of normalized weights is larger than the requested
  !! probability mass.
  !!
  SUBROUTINE confidence_limits_sample_r8(indata, pdf, mask, &
       probability_mass, peak, bounds, errstr)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in)            :: indata
    REAL(rprec8), DIMENSION(:), INTENT(in)            :: pdf
    REAL(rprec8), INTENT(in), OPTIONAL                :: probability_mass
    LOGICAL, DIMENSION(:), OPTIONAL, INTENT(in)       :: mask
    REAL(rprec8), OPTIONAL, INTENT(out)               :: peak
    REAL(rprec8), DIMENSION(2), OPTIONAL, INTENT(out) :: bounds
    CHARACTER(len=*), INTENT(inout)                   :: errstr

    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: pdf_
    REAL(rprec8) :: probability_mass_
    INTEGER, DIMENSION(:), ALLOCATABLE :: indx_arr
    INTEGER :: ndata, imax, i, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask_

    ndata = SIZE(indata)
    IF (ndata /= SIZE(pdf)) THEN
       errstr = " -> statistics : confidence_limits : Size of vectors does not conform." // &
            TRIM(errstr)
       RETURN
    END IF
    ALLOCATE(pdf_(ndata), mask_(ndata), stat=err)
    IF (err /= 0) THEN
       errstr = " -> statistics : confidence_limits : Could not allocate memory." // &
            TRIM(errstr)
       DEALLOCATE(pdf_, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF
    mask_ = .TRUE.
    IF (PRESENT(mask)) THEN
       IF (ndata /= SIZE(mask)) THEN
          errstr = " -> statistics : confidence_limits : Size of mask does not conform with data." // &
               TRIM(errstr)
          DEALLOCATE(pdf_, stat=err)
          DEALLOCATE(mask_, stat=err)
          RETURN
       END IF
       mask_ = mask
    ELSE
       mask_ = .TRUE.
    END IF
    pdf_ = pdf
    WHERE (.NOT.mask_)
       pdf_ = 0.0_rprec8
    END WHERE
    pdf_ = pdf_/SUM(pdf_)

    IF (PRESENT(bounds) .AND. .NOT.PRESENT(probability_mass)) THEN
       errstr = " -> statistics : confidence_limits : Probability mass not given." // &
            TRIM(errstr)
       DEALLOCATE(pdf_, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF

    IF (PRESENT(peak)) THEN
       ! Find peak of phase angle histogram:
       imax = imaxloc(pdf_)
       peak = indata(imax)
    END IF
    IF (PRESENT(bounds)) THEN
       ALLOCATE(indx_arr(ndata), stat=err)
       IF (err /= 0) THEN
          errstr = " -> statistics : confidence_limits : Could not allocate memory." // &
               TRIM(errstr)
          DEALLOCATE(pdf_, stat=err)
          DEALLOCATE(mask_, stat=err)
          DEALLOCATE(indx_arr, stat=err)
          RETURN
       END IF
       CALL quicksort(pdf_, indx_arr, errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          errstr = " -> statistics : confidence_limits : ." // &
               TRIM(errstr)
          DEALLOCATE(pdf_, stat=err)
          DEALLOCATE(mask_, stat=err)
          DEALLOCATE(indx_arr, stat=err)
          RETURN
       END IF
       pdf_ = pdf_/SUM(pdf_)
       bounds(1) = HUGE(bounds(1))
       bounds(2) = -HUGE(bounds(2))
       probability_mass_ = 0.0_rprec8
       DO i=ndata,1,-1
          IF (mask_(indx_arr(i))) THEN
             probability_mass_ = probability_mass_ + pdf_(indx_arr(i))
             IF (indata(indx_arr(i)) < bounds(1)) THEN
                bounds(1) = indata(indx_arr(i))
             END IF
             IF (indata(indx_arr(i)) > bounds(2)) THEN
                bounds(2) = indata(indx_arr(i))               
             END IF
             IF (probability_mass_ >= probability_mass) THEN
                EXIT
             END IF
          END IF
       END DO
       DEALLOCATE(indx_arr, stat=err)
       IF (err /= 0) THEN
          errstr = " -> statistics : confidence_limits : Could not deallocate memory (1)." // &
               TRIM(errstr)
          DEALLOCATE(pdf_, stat=err)
          DEALLOCATE(mask_, stat=err)
          DEALLOCATE(indx_arr, stat=err)
          RETURN
       END IF
    END IF

    DEALLOCATE(pdf_, mask_, stat=err)
    IF (err /= 0) THEN
       errstr = " -> statistics : confidence_limits : Could not deallocate memory (2)." // &
            TRIM(errstr)
       DEALLOCATE(pdf_, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF

  END SUBROUTINE confidence_limits_sample_r8





  SUBROUTINE credible_region(pdf_arr, probability_mass, indx_arr, errstr, repetition_arr)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:), INTENT(in) :: pdf_arr
    REAL(rprec8), INTENT(in) :: probability_mass
    INTEGER, DIMENSION(:), INTENT(out) :: indx_arr
    CHARACTER(len=*), INTENT(inout) :: errstr
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: repetition_arr

    REAL(rprec8), DIMENSION(:), ALLOCATABLE :: pdf_arr_
    REAL(rprec8) :: probability_mass_
    INTEGER :: i, j, irepet, nrepet

    ! Sort the pdf in ascending order:
    CALL quickSort(pdf_arr, indx_arr, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       errstr = " -> statistics : credible_region : TRACE BACK." // &
            TRIM(errstr)
       RETURN
    END IF
    ! Start from the largest pdf value
    i = SIZE(pdf_arr)
    IF (PRESENT(repetition_arr)) THEN
       ! MCMC 
       ! Which number of repetitions corresponds to the wanted
       ! probability mass?
       nrepet = CEILING(SUM(repetition_arr)*probability_mass)
       ! Add distinct sampling points until reaching the desired
       ! number of (typically non-unique) sample points:
       irepet = 0
       DO WHILE (irepet < nrepet)
          irepet = irepet + repetition_arr(indx_arr(i))
          i = i - 1
       END DO
    ELSE
       ! MC
       ALLOCATE(pdf_arr_(SIZE(pdf_arr)))
       ! Normalize the pdf distribution
       pdf_arr_ = pdf_arr/SUM(pdf_arr)
       ! Add sample points until the desired probability mass is
       ! reached:
       probability_mass_ = 0.0_rprec8
       DO WHILE (probability_mass_ < probability_mass)
          probability_mass_ = probability_mass_ + pdf_arr_(indx_arr(i))
          i = i - 1
       END DO
       DEALLOCATE(pdf_arr_)
    END IF
    ! Fill the resulting index array with negative values for the
    ! sample points that do not fit within the desired probability
    ! mass:
    DO j=i,1,-1
       indx_arr(j) = -1
    END DO

  END SUBROUTINE credible_region





  SUBROUTINE population_covariance(indata, covariance, mean, mask)

    IMPLICIT NONE
    REAL(rprec8), DIMENSION(:,:), INTENT(in) :: indata
    REAL(rprec8), DIMENSION(:,:), INTENT(out) :: covariance
    REAL(rprec8), DIMENSION(:), INTENT(out) :: mean
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: mask

    INTEGER :: i, j
    LOGICAL,DIMENSION(SIZE(indata,dim=2)) :: mask_

    IF (PRESENT(mask)) THEN
       mask_ = mask
       covariance = 0.0_rprec8
    ELSE
       mask_ = .TRUE.
    END IF

    DO i=1,SIZE(indata,dim=2)
       ! Compute mean for variable i
       mean(i) = SUM(indata(:,i))/SIZE(indata,dim=1)
       ! Compute elements of covariance matrix from 1,1 to i,i
       DO j=1,i
          IF (mask_(i) .AND. mask_(j)) THEN
             covariance(i,j) = SUM((indata(:,i)-mean(i)) * (indata(:,j)-mean(j))) / &
                  SIZE(indata,dim=1)
             covariance(j,i) = covariance(i,j) 
          END IF
       END DO
    END DO

  END SUBROUTINE population_covariance

  ! Given the covariance matrix and residuals of a given observation,
  ! computes the Mahalanobis distance.

  REAL(rprec8) FUNCTION mahalanobis_distance_r8(cov_matrix,residuals, errstr)

    IMPLICIT NONE

    REAL(rprec8), DIMENSION(:,:), INTENT(in)  :: cov_matrix, residuals
    CHARACTER(len=*), INTENT(inout) :: errstr

    REAL(rprec8), DIMENSION(:,:), ALLOCATABLE :: inverse_cov_matrix
    REAL(rprec8), DIMENSION(1,1)                :: mahalanobis
    INTEGER n

    n = size(cov_matrix,dim=1)
    ALLOCATE(inverse_cov_matrix(n,n))

    inverse_cov_matrix = matinv(cov_matrix,errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       errstr = " -> statistics : mahalanobis_distance : TRACE BACK." // &
            TRIM(errstr)
       RETURN
    END IF

    mahalanobis = SQRT(MATMUL(MATMUL(TRANSPOSE(residuals),inverse_cov_matrix),(residuals)))
    mahalanobis_distance_r8 = mahalanobis(1,1)
    DEALLOCATE(inverse_cov_matrix)

  END FUNCTION mahalanobis_distance_r8



END MODULE statistics
