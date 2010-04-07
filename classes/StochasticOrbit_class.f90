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
!! *Class*description*: 
!!  
!! Type and routines for stochastic orbits (that is, orbital elements
!! with their uncertainties). Contains the algorithms for the
!! [statistical orbital] ranging method and the least-squares method.
!!
!! @author MG, JV, KM 
!! @version 2010-04-06
!!  
MODULE StochasticOrbit_cl

  USE Base_cl
  USE Orbit_cl
  USE Observations_cl
  USE SphericalCoordinates_cl
  USE CartesianCoordinates_cl
  USE Time_cl

  USE planetary_data
  USE linal
  USE estimators
  USE statistics

  IMPLICIT NONE
  !! (Initial) maximum probability density value.
  !! Note: case-dependent, you should try different values
  !!       by changing SO%pdf_ml_init_prm.
  REAL(bp), PARAMETER, PRIVATE :: pdf_ml_init  = 50.0_bp
  !! Default bound for equivalent chi2-difference. Change it through
  !! setParameters!
  !! E.g., 50.0 corresponds to probability mass (1 - 4.7e-9).  
  !! E.g., 30.0 corresponds to probability mass (1 - ?).
  REAL(bp), PARAMETER, PRIVATE :: dchi2       = 30.0_bp    
  !! Maximum end value for histogram.
  REAL(bp), PARAMETER, PRIVATE :: histo_end   = 0.3_bp     
  !! Maximum number of orbits to be processed simultaneously.
  INTEGER, PARAMETER, PRIVATE  :: norb_simult_max = 5000

  PRIVATE :: new_SO
  PRIVATE :: new_SO_observations
  PRIVATE :: new_SO_orb_arr
  PRIVATE :: nullify_SO
  PRIVATE :: copy_SO
  PRIVATE :: exist_SO
  PRIVATE :: getApoapsisDistance_SO
  PRIVATE :: getEphemerides_SO
  PRIVATE :: getEphemeris_SO
  PRIVATE :: getCovarianceMatrix_SO
  PRIVATE :: getParameters_SO
  PRIVATE :: getPeriapsisDistance_SO
  PRIVATE :: getPhaseAngle_SO_pdf
  PRIVATE :: getPhaseAngle_SO_point
  PRIVATE :: getPhaseAngles_SO
  PRIVATE :: getObservationMasks_SO
  PRIVATE :: getRangeBounds_SO
  PRIVATE :: getResults_SO
  PRIVATE :: getSampleOrbit_SO
  PRIVATE :: getStandardDeviations_SO
  PRIVATE :: propagate_SO
  PRIVATE :: setObservationMask_one
  PRIVATE :: setObservationMask_all
  PRIVATE :: setObservationMask_all_notes
  PRIVATE :: setObservationPair_default
  PRIVATE :: setObservationPair_pair
  PRIVATE :: setParameters_SO
  PRIVATE :: toCartesian_SO
  PRIVATE :: toKeplerian_SO

  TYPE StochasticOrbit
     PRIVATE
     TYPE (Time)                         :: t_inv_prm
     TYPE (Orbit), DIMENSION(:), POINTER :: orb_arr_cmp            => NULL()
     TYPE (Orbit)                        :: orb_ml_cmp
     TYPE (Observations)                 :: obss
     CHARACTER(len=ELEMENT_TYPE_LEN)     :: cov_type_prm           = " "
     CHARACTER(len=ELEMENT_TYPE_LEN)     :: element_type_prm       = "cartesian"
     REAL(bp), DIMENSION(:,:,:), POINTER :: res_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: cov_ml_cmp             => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: rms_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: res_accept_prm         => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: jac_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: rchi2_arr_cmp      => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: pdf_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: reg_apr_arr_cmp        => NULL()
     REAL(bp)                            :: pdf_ml_init_prm        = -1.0_bp
     REAL(bp)                            :: dchi2_prm              = dchi2
     REAL(bp)                            :: pdf_ml_prm             = -1.0_bp
     REAL(bp)                            :: pdf_ml_cmp             = -1.0_bp
     REAL(bp)                            :: accept_multiplier_prm  = -1.0_bp
     REAL(bp)                            :: outlier_multiplier_prm = -1.0_bp
     LOGICAL, DIMENSION(:,:), POINTER    :: obs_masks_prm          => NULL()
     LOGICAL                             :: outlier_rejection_prm  = .FALSE.
     LOGICAL                             :: regularization_prm     = .TRUE.
     LOGICAL                             :: jacobians_prm          = .TRUE.
     LOGICAL                             :: multiple_obj_prm       = .FALSE.
     LOGICAL                             :: is_initialized_prm     = .FALSE.
     LOGICAL                             :: uniform_pdf_prm        = .FALSE.

     ! Bayesian informative apriori assumptions
     REAL(bp)                            :: apriori_a_max_prm                    = 500.0
     REAL(bp)                            :: apriori_a_min_prm                    = planetary_radii(11)
     REAL(bp)                            :: apriori_periapsis_max_prm           = -1.0_bp
     REAL(bp)                            :: apriori_periapsis_min_prm           = -1.0_bp
     REAL(bp)                            :: apriori_apoapsis_max_prm             = -1.0_bp
     REAL(bp)                            :: apriori_apoapsis_min_prm             = -1.0_bp
     REAL(bp)                            :: apriori_rho_max_prm                  = 200.0_bp
     REAL(bp)                            :: apriori_rho_min_prm                  = 0.0_bp
     REAL(bp)                            :: apriori_hcentric_dist_max_prm = 200.0_bp
     REAL(bp)                            :: apriori_hcentric_dist_min_prm = planetary_radii(11)
     REAL(bp)                            :: apriori_velocity_max_prm             = 0.12_bp
     LOGICAL                             :: informative_apriori_prm              = .TRUE.

     ! Parameters for propagation:
     CHARACTER(len=DYN_MODEL_LEN)        :: dyn_model_prm        = "2-body"
     CHARACTER(len=INTEGRATOR_LEN)       :: integrator_prm       = "gauss-radau"
     REAL(bp), DIMENSION(:), POINTER     :: finite_diff_prm      => NULL()
     REAL(bp)                            :: integration_step_prm = 1.0_bp
     LOGICAL, DIMENSION(10)              :: perturbers_prm       = .FALSE.

     ! Parameters for statistical ranging:
     CHARACTER(len=64)                   :: sor_2point_method_prm      = "continued fraction"
     CHARACTER(len=64)                   :: sor_2point_method_sw_prm   = "continued fraction"
     REAL(bp), DIMENSION(:,:,:), POINTER :: sor_deviates_prm           => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: sor_rho_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: sor_pair_histo_prm         => NULL()
     REAL(bp), DIMENSION(2,2)            :: sor_rho_prm                = -1.0_bp
     REAL(bp), DIMENSION(2,2)            :: sor_rho_cmp                = -1.0_bp
     REAL(bp)                            :: sor_generat_multiplier_prm = -1.0_bp
     INTEGER, DIMENSION(:,:), POINTER    :: sor_pair_arr_prm           => NULL()
     INTEGER                             :: sor_ntrial_prm             = -1
     INTEGER                             :: sor_ntrial_cmp             = -1
     INTEGER                             :: sor_norb_prm               = -1
     INTEGER                             :: sor_norb_cmp               = -1
     INTEGER                             :: sor_ntrial_sw_prm          = -1
     INTEGER                             :: sor_ntrial_sw_cmp          = -1
     INTEGER                             :: sor_norb_sw_prm            = -1
     INTEGER                             :: sor_norb_sw_cmp            = -1
     INTEGER                             :: sor_niter_cmp              = -1
     INTEGER                             :: sor_niter_prm              = -1
     INTEGER                             :: sor_rho_histo_cmp          = 1
     LOGICAL, DIMENSION(4)               :: sor_iterate_bounds_prm     = .TRUE.
     LOGICAL                             :: sor_random_obs_prm         = .FALSE.
     ! Allow generation of rho from Gaussian p.d.f.
     ! NOTE: only applicable to first statisticalRanging in autoSR!
     LOGICAL                             :: sor_gaussian_pdf_prm       = .FALSE.

     ! Parameters for VoV:
     REAL(bp), DIMENSION(:,:), POINTER   :: vov_map_cmp           => NULL()
     REAL(bp), DIMENSION(6,2)            :: vov_scaling_prm       = -1.0_bp
     REAL(bp), DIMENSION(6,2)            :: vov_scaling_cmp       = -1.0_bp
     INTEGER                             :: vov_norb_prm          = -1
     INTEGER                             :: vov_norb_iter_prm     = -1
     INTEGER                             :: vov_norb_cmp          = -1
     INTEGER                             :: vov_ntrial_prm        = -1
     INTEGER                             :: vov_ntrial_iter_prm   = -1
     INTEGER                             :: vov_ntrial_cmp        = -1
     INTEGER                             :: vov_niter_prm         = -1
     INTEGER                             :: vov_niter_cmp         = -1
     INTEGER                             :: vov_nmap_prm          = -1
     LOGICAL, DIMENSION(6)               :: vov_mapping_mask_prm  = &
          (/ .TRUE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /)
     LOGICAL, DIMENSION(6,2)             :: vov_scaling_ready_cmp = &
          .FALSE.

     ! Parameters for the least-squares fitting:
     REAL(bp)                            :: ls_corr_fac_prm         = 1.0_bp
     REAL(bp)                            :: ls_rchi2_diff_tresh_prm = 0.00001_bp
     INTEGER                             :: ls_niter_major_max_prm  = 10
     INTEGER                             :: ls_niter_major_min_prm  = 2
     INTEGER                             :: ls_niter_minor_prm      = 20
     LOGICAL, DIMENSION(6)               :: ls_elem_mask_prm        = .TRUE.

     ! Parameters for the covariance sampling:
     REAL(bp)                            :: cos_nsigma_prm   = 7.5_bp
     INTEGER                             :: cos_norb_prm     = 100
     INTEGER                             :: cos_ntrial_prm   = 100000
     LOGICAL                             :: cos_gaussian_prm = .FALSE.

     ! Parameters for simplex optimization:
     REAL(bp)                            :: smplx_tol_prm  = 1.05_bp
     INTEGER                             :: smplx_niter_prm = 1000
     INTEGER                             :: smplx_niter_cmp

  END TYPE StochasticOrbit


  INTERFACE NEW
     MODULE PROCEDURE new_SO
     MODULE PROCEDURE new_SO_observations
     MODULE PROCEDURE new_SO_orb_cov
     MODULE PROCEDURE new_SO_orb_arr
  END INTERFACE

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_SO
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_SO
  END INTERFACE

  INTERFACE exist
     MODULE PROCEDURE exist_SO
  END INTERFACE

  INTERFACE getApoapsisDistance
     MODULE PROCEDURE getApoapsisDistance_SO
  END INTERFACE

  INTERFACE getChi2
     MODULE PROCEDURE getChi2_matrix
     MODULE PROCEDURE getChi2_this_orb
  END INTERFACE

  INTERFACE getCovarianceMatrix
     MODULE PROCEDURE getCovarianceMatrix_SO
  END INTERFACE

  INTERFACE getEphemerides
     MODULE PROCEDURE getEphemerides_SO
  END INTERFACE

  INTERFACE getEphemeris
     MODULE PROCEDURE getEphemeris_SO
  END INTERFACE

  INTERFACE getParameters
     MODULE PROCEDURE getParameters_SO
  END INTERFACE

  INTERFACE getPeriapsisDistance
     MODULE PROCEDURE getPeriapsisDistance_SO
  END INTERFACE

  INTERFACE getPhaseAngle
     MODULE PROCEDURE getPhaseAngle_SO_pdf
     MODULE PROCEDURE getPhaseAngle_SO_point
  END INTERFACE

  INTERFACE getPhaseAngles
     MODULE PROCEDURE getPhaseAngles_SO
  END INTERFACE

  INTERFACE getObservationMasks
     MODULE PROCEDURE getObservationMasks_SO
  END INTERFACE

  INTERFACE getRangeBounds
     MODULE PROCEDURE getRangeBounds_SO
  END INTERFACE

  INTERFACE getResiduals
     MODULE PROCEDURE getResiduals_single
  END INTERFACE

  INTERFACE getResults
     MODULE PROCEDURE getResults_SO
  END INTERFACE

  INTERFACE getRMS
     MODULE PROCEDURE getRMS_single
  END INTERFACE

  INTERFACE getSampleOrbit
     MODULE PROCEDURE getSampleOrbit_SO
  END INTERFACE

  INTERFACE getStandardDeviations
     MODULE PROCEDURE getStandardDeviations_SO
  END INTERFACE

  INTERFACE leastSquares
     MODULE PROCEDURE leastSquares_SO
  END INTERFACE

  INTERFACE levenbergMarquardt
     MODULE PROCEDURE levenbergMarquardt_SO
  END INTERFACE

  INTERFACE getTime
     MODULE PROCEDURE getTime_SO
  END INTERFACE

  INTERFACE propagate
     MODULE PROCEDURE propagate_SO
  END INTERFACE

  INTERFACE setAcceptanceWindow
     MODULE PROCEDURE setAcceptanceWindow_sigma
  END INTERFACE

  INTERFACE setObservationMask
     MODULE PROCEDURE setObservationMask_one
     MODULE PROCEDURE setObservationMask_all
     MODULE PROCEDURE setObservationMask_all_notes
  END INTERFACE

  INTERFACE setObservationPair
     MODULE PROCEDURE setObservationPair_default
     MODULE PROCEDURE setObservationPair_pair
  END INTERFACE

  INTERFACE setParameters
     MODULE PROCEDURE setParameters_SO
  END INTERFACE

  INTERFACE setRangeBounds
     MODULE PROCEDURE setRangeBounds_3sigma
     MODULE PROCEDURE setRangeBounds_values
  END INTERFACE

  INTERFACE toCartesian
     MODULE PROCEDURE toCartesian_SO
  END INTERFACE

  INTERFACE toCometary
     MODULE PROCEDURE toCometary_SO
  END INTERFACE

  INTERFACE toKeplerian
     MODULE PROCEDURE toKeplerian_SO
  END INTERFACE

CONTAINS





  !! *Description*:
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this

    IF (this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%is_initialized_prm = .TRUE.

  END SUBROUTINE new_SO





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SO_observations(this, obss)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Observations), INTENT(in)       :: obss

    IF (this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF(.NOT. exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Observations not intialized.", 1)
       RETURN
    END IF

    this%obss = copy(obss)
    this%obs_masks_prm => getObservationMasks(this%obss)

    ! Initialize other variables needed in the computation
    this%is_initialized_prm = .TRUE.
    CALL setObservationPair(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / new", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END SUBROUTINE new_SO_observations





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SO_orb_cov(this, orb, cov, cov_type, &
       element_type, obss)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)     :: this
    TYPE (Orbit), INTENT(in)                  :: orb
    CHARACTER(len=*), INTENT(in)              :: cov_type
    REAL(bp), DIMENSION(6,6), INTENT(in)      :: cov
    CHARACTER(len=*), INTENT(in), OPTIONAL    :: element_type
    TYPE (Observations), INTENT(in), OPTIONAL :: obss
    INTEGER                                   :: err

    IF (this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(obss)) THEN
       IF(.NOT. exist(obss)) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / new", &
               "Observations not intialized.", 1)
          RETURN
       END IF
       this%obss = copy(obss)
       this%obs_masks_prm => getObservationMasks(this%obss)
    END IF
    ALLOCATE(this%cov_ml_cmp(6,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    ! Initialize other variables needed in the computation
    this%orb_ml_cmp = copy(orb)
    this%cov_ml_cmp = cov
    this%cov_type_prm = cov_type
    IF (PRESENT(element_type)) THEN
       this%element_type_prm = element_type
    ELSE
       this%element_type_prm = "cartesian"
    END IF
    this%is_initialized_prm = .TRUE.
    IF (PRESENT(obss)) THEN
       CALL setObservationPair(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit /  new", &
               "TRACE BACK", 1)
          RETURN
       END IF
    END IF

  END SUBROUTINE new_SO_orb_cov





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SO_orb_arr(this, orb_arr, pdf_arr, element_type, jac_arr, reg_apr_arr, rchi2_arr)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)          :: this
    TYPE (Orbit), DIMENSION(:), INTENT(in)         :: orb_arr
    CHARACTER(len=*), INTENT(in)                   :: element_type
    REAL(bp), DIMENSION(:), INTENT(in)             :: pdf_arr
    REAL(bp), DIMENSION(:,:), INTENT(in), OPTIONAL :: jac_arr
    REAL(bp), DIMENSION(:), INTENT(in), OPTIONAL   :: reg_apr_arr
    REAL(bp), DIMENSION(:), INTENT(in), OPTIONAL   :: rchi2_arr

    INTEGER :: i, err

    IF (this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%sor_norb_cmp = SIZE(orb_arr, dim=1)
    ALLOCATE(this%orb_arr_cmp(this%sor_norb_cmp), &
         this%pdf_arr_cmp(this%sor_norb_cmp), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / new", &
            "Could not allocate pointer (5).", 1)
       RETURN
    END IF
    DO i=1,this%sor_norb_cmp
       this%orb_arr_cmp(i) = copy(orb_arr(i))
       this%pdf_arr_cmp(i) = pdf_arr(i)
    END DO
    IF (PRESENT(jac_arr)) THEN
       ALLOCATE(this%jac_arr_cmp(this%sor_norb_cmp,3), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / new", &
               "Could not allocate pointer (10).", 1)
          RETURN
       END IF
       DO i=1,this%sor_norb_cmp
          this%jac_arr_cmp(i,1) = jac_arr(i,1)
          this%jac_arr_cmp(i,2) = jac_arr(i,2)
          this%jac_arr_cmp(i,3) = jac_arr(i,3)
       END DO
    END IF
    IF (PRESENT(rchi2_arr)) THEN
       ALLOCATE(this%rchi2_arr_cmp(this%sor_norb_cmp), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / new", &
               "Could not allocate pointer (15).", 1)
          RETURN
       END IF
       DO i=1,this%sor_norb_cmp
          this%rchi2_arr_cmp(i) = rchi2_arr(i)
       END DO
    END IF
    IF (PRESENT(reg_apr_arr)) THEN
       ALLOCATE(this%reg_apr_arr_cmp(this%sor_norb_cmp), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / new", &
               "Could not allocate pointer (20).", 1)
          RETURN
       END IF
       DO i=1,this%sor_norb_cmp
          this%reg_apr_arr_cmp(i) = reg_apr_arr(i)
       END DO
    END IF
    this%element_type_prm = element_type
    this%is_initialized_prm = .TRUE.

  END SUBROUTINE new_SO_orb_arr






  !! *Description*:
  !!
  SUBROUTINE nullify_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    INTEGER :: i, err

    CALL NULLIFY(this%t_inv_prm)
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DO i=1,SIZE(this%orb_arr_cmp)
          CALL NULLIFY(this%orb_arr_cmp(i))
       END DO
       DEALLOCATE(this%orb_arr_cmp, stat=err)
    END IF
    CALL NULLIFY(this%orb_ml_cmp)
    CALL NULLIFY(this%obss)
    this%element_type_prm = "cartesian"
    this%cov_type_prm = " "
    IF (ASSOCIATED(this%res_arr_cmp)) THEN
       DEALLOCATE(this%res_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%res_accept_prm)) THEN
       DEALLOCATE(this%res_accept_prm, stat=err)
    END IF
    IF (ASSOCIATED(this%cov_ml_cmp)) THEN
       DEALLOCATE(this%cov_ml_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%rms_arr_cmp)) THEN
       DEALLOCATE(this%rms_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%reg_apr_arr_cmp)) THEN
       DEALLOCATE(this%reg_apr_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%jac_arr_cmp)) THEN
       DEALLOCATE(this%jac_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%obs_masks_prm)) THEN
       DEALLOCATE(this%obs_masks_prm, stat=err)
    END IF
    this%pdf_ml_prm = -1.0_bp
    this%pdf_ml_init_prm = -1.0_bp
    this%pdf_ml_cmp = -1.0_bp
    this%accept_multiplier_prm = -1.0_bp
    this%outlier_rejection_prm = .FALSE.
    this%regularization_prm = .TRUE.
    this%jacobians_prm = .TRUE.
    this%multiple_obj_prm = .FALSE.

    ! Bayesian apriori parameters
    this%apriori_a_max_prm = 500.0
    this%apriori_a_min_prm = planetary_radii(11)
    this%apriori_periapsis_max_prm = -1.0_bp
    this%apriori_periapsis_min_prm = -1.0_bp
    this%apriori_apoapsis_max_prm = -1.0_bp
    this%apriori_apoapsis_min_prm = -1.0_bp
    this%apriori_rho_max_prm = 200.0_bp
    this%apriori_rho_min_prm = 0.0_bp
    this%informative_apriori_prm = .TRUE.

    ! Propagation parameters
    this%dyn_model_prm = "2-body"
    this%integrator_prm = "gauss-radau"
    IF (ASSOCIATED(this%finite_diff_prm)) THEN
       DEALLOCATE(this%finite_diff_prm, stat=err)
    END IF
    this%integration_step_prm = 1.0_bp 
    this%perturbers_prm = .FALSE.

    ! Ranging variables
    this%sor_2point_method_prm = "continued fraction"
    this%sor_2point_method_sw_prm = "continued fraction"
    IF (ASSOCIATED(this%sor_deviates_prm)) THEN
       DEALLOCATE(this%sor_deviates_prm, stat=err)
    END IF
    IF (ASSOCIATED(this%sor_rho_arr_cmp)) THEN
       DEALLOCATE(this%sor_rho_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%sor_pair_histo_prm)) THEN
       DEALLOCATE(this%sor_pair_histo_prm, stat=err)
    END IF
    this%sor_rho_prm = -1.0_bp
    this%sor_rho_cmp = -1.0_bp
    this%sor_generat_multiplier_prm = -1.0_bp
    IF (ASSOCIATED(this%sor_pair_arr_prm)) THEN
       DEALLOCATE(this%sor_pair_arr_prm, stat=err)
    END IF
    this%sor_ntrial_prm = -1
    this%sor_ntrial_cmp = -1
    this%sor_norb_prm = -1
    this%sor_norb_cmp = -1
    this%sor_niter_cmp = -1
    this%sor_niter_prm = -1
    this%sor_rho_histo_cmp = 1
    this%sor_random_obs_prm = .FALSE.
    this%sor_gaussian_pdf_prm = .FALSE.
    this%uniform_pdf_prm = .FALSE.

    ! Parameters for VoV:
    IF (ASSOCIATED(this%vov_map_cmp)) THEN
       DEALLOCATE(this%vov_map_cmp, stat=err)
    END IF
    this%vov_scaling_prm = -1.0_bp
    this%vov_scaling_cmp = -1.0_bp
    this%vov_norb_prm = -1
    this%vov_norb_cmp = -1
    this%vov_ntrial_prm = -1
    this%vov_ntrial_cmp = -1
    this%vov_niter_prm = -1
    this%vov_niter_cmp = -1
    this%vov_nmap_prm = -1
    this%vov_mapping_mask_prm = &
         (/ .TRUE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /)
    this%vov_scaling_ready_cmp = .FALSE.

    ! Parameters for the least-squares fit:
    this%ls_corr_fac_prm = 1.0_bp
    this%ls_niter_major_max_prm = 10
    this%ls_niter_major_min_prm = 2
    this%ls_niter_minor_prm = 20
    this%ls_elem_mask_prm = .TRUE.

    this%is_initialized_prm = .FALSE.

  END SUBROUTINE nullify_SO





  !! *Description*:
  !!
  !!
  !! Returns error if allocation fails.
  !!
  FUNCTION copy_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit) :: this
    TYPE (StochasticOrbit) :: copy_SO
    INTEGER :: i, err, nobs, norb

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / copy", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    copy_SO%t_inv_prm = copy(this%t_inv_prm) 
    IF (exist(this%obss)) THEN
       copy_SO%obss = copy(this%obss)
       nobs = getNrOfObservations(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / copy", &
               "TRACE BACK (5).", 1)
          RETURN
       END IF
    END IF
    norb = SIZE(this%orb_arr_cmp)
    IF (ASSOCIATED(copy_SO%orb_arr_cmp)) THEN
       DO i=1,SIZE(copy_SO%orb_arr_cmp)
          CALL NULLIFY(copy_SO%orb_arr_cmp(i))
       END DO
       DEALLOCATE(copy_SO%orb_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       ALLOCATE(copy_SO%orb_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (5).", 1)
          RETURN
       END IF
       DO i=1,norb
          copy_SO%orb_arr_cmp(i) = copy(this%orb_arr_cmp(i))       
       END DO
    END IF
    copy_SO%orb_ml_cmp = copy(this%orb_ml_cmp)
    copy_SO%cov_type_prm = this%cov_type_prm
    IF (ASSOCIATED(this%res_arr_cmp)) THEN
       ALLOCATE(copy_SO%res_arr_cmp(norb,nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (10).", 1)
          RETURN
       END IF
       copy_SO%res_arr_cmp = this%res_arr_cmp
    END IF
    IF (ASSOCIATED(this%res_accept_prm)) THEN
       ALLOCATE(copy_SO%res_accept_prm(nobs,6))
       copy_SO%res_accept_prm = this%res_accept_prm
    END IF
    IF (ASSOCIATED(this%cov_ml_cmp)) THEN
       ALLOCATE(copy_SO%cov_ml_cmp(6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (15).", 1)
          RETURN
       END IF
       copy_SO%cov_ml_cmp = this%cov_ml_cmp
    END IF
    IF (ASSOCIATED(this%rms_arr_cmp)) THEN
       ALLOCATE(copy_SO%rms_arr_cmp(norb,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (25).", 1)
          RETURN
       END IF
       copy_SO%rms_arr_cmp = this%rms_arr_cmp
    END IF
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       ALLOCATE(copy_SO%rchi2_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (30).", 1)
          RETURN
       END IF
       copy_SO%rchi2_arr_cmp = this%rchi2_arr_cmp
    END IF
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       ALLOCATE(copy_SO%pdf_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (35).", 1)
          RETURN
       END IF
       copy_SO%pdf_arr_cmp = this%pdf_arr_cmp
    END IF
    IF (ASSOCIATED(this%reg_apr_arr_cmp)) THEN
       ALLOCATE(copy_SO%reg_apr_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (40).", 1)
          RETURN
       END IF
       copy_SO%reg_apr_arr_cmp = this%reg_apr_arr_cmp
    END IF
    IF (ASSOCIATED(this%jac_arr_cmp)) THEN
       ALLOCATE(copy_SO%jac_arr_cmp(norb,3), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (45).", 1)
          RETURN
       END IF
       copy_SO%jac_arr_cmp = this%jac_arr_cmp
    END IF
    IF (ASSOCIATED(this%obs_masks_prm)) THEN
       ALLOCATE(copy_SO%obs_masks_prm(nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (20).", 1)
          RETURN
       END IF
       copy_SO%obs_masks_prm = this%obs_masks_prm
    END IF
    copy_SO%pdf_ml_init_prm = this%pdf_ml_init_prm
    copy_SO%pdf_ml_prm = this%pdf_ml_prm
    copy_SO%pdf_ml_cmp = this%pdf_ml_cmp
    copy_SO%accept_multiplier_prm = this%accept_multiplier_prm
    copy_SO%outlier_rejection_prm = this%outlier_rejection_prm
    copy_SO%regularization_prm = this%regularization_prm
    copy_SO%jacobians_prm = this%jacobians_prm
    copy_SO%multiple_obj_prm = this%multiple_obj_prm

    ! Bayesian apriori parameters
    copy_SO%apriori_a_max_prm = this%apriori_a_max_prm
    copy_SO%apriori_a_min_prm = this%apriori_a_min_prm
    copy_SO%apriori_periapsis_max_prm = this%apriori_periapsis_max_prm
    copy_SO%apriori_periapsis_min_prm = this%apriori_periapsis_min_prm
    copy_SO%apriori_apoapsis_max_prm = this%apriori_apoapsis_max_prm
    copy_SO%apriori_apoapsis_min_prm = this%apriori_apoapsis_min_prm
    copy_SO%apriori_rho_min_prm = this%apriori_rho_min_prm
    copy_SO%informative_apriori_prm = this%informative_apriori_prm

    ! Propagation parameters
    copy_SO%dyn_model_prm = this%dyn_model_prm
    copy_SO%integrator_prm = this%integrator_prm
    IF (ASSOCIATED(this%finite_diff_prm)) THEN
       ALLOCATE(copy_SO%finite_diff_prm(6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (50).", 1)
          RETURN
       END IF
       copy_SO%finite_diff_prm = this%finite_diff_prm
    END IF
    copy_SO%integration_step_prm = this%integration_step_prm
    copy_SO%perturbers_prm = this%perturbers_prm

    ! Ranging parameters
    copy_SO%element_type_prm = this%element_type_prm
    copy_SO%sor_2point_method_prm = this%sor_2point_method_prm
    copy_SO%sor_2point_method_sw_prm = this%sor_2point_method_sw_prm
    IF (ASSOCIATED(this%sor_deviates_prm)) THEN
       ALLOCATE(copy_SO%sor_deviates_prm(nobs,6,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (55).", 1)
          RETURN
       END IF
       copy_SO%sor_deviates_prm = this%sor_deviates_prm
    END IF
    IF (ASSOCIATED(this%sor_rho_arr_cmp)) THEN
       ALLOCATE(copy_SO%sor_rho_arr_cmp(norb,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (60).", 1)
          RETURN
       END IF
       copy_SO%sor_rho_arr_cmp = this%sor_rho_arr_cmp
    END IF
    IF (ASSOCIATED(this%sor_pair_histo_prm)) THEN
       ALLOCATE(copy_SO%sor_pair_histo_prm(SIZE(this%sor_pair_histo_prm)), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (65).", 1)
          RETURN
       END IF
       copy_SO%sor_pair_histo_prm = this%sor_pair_histo_prm
    END IF
    copy_SO%sor_rho_prm = this%sor_rho_prm
    copy_SO%sor_rho_cmp = this%sor_rho_cmp
    copy_SO%sor_generat_multiplier_prm = this%sor_generat_multiplier_prm
    IF (ASSOCIATED(this%sor_pair_arr_prm)) THEN
       ALLOCATE(copy_SO%sor_pair_arr_prm(SIZE(this%sor_pair_arr_prm,dim=1), &
            SIZE(this%sor_pair_arr_prm,dim=2)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (70).", 1)
          RETURN
       END IF
       copy_SO%sor_pair_arr_prm = this%sor_pair_arr_prm
    END IF
    copy_SO%sor_ntrial_prm = this%sor_ntrial_prm
    copy_SO%sor_ntrial_cmp = this%sor_ntrial_cmp
    copy_SO%sor_norb_prm = this%sor_norb_prm
    copy_SO%sor_norb_cmp = this%sor_norb_cmp
    copy_SO%sor_niter_cmp = this%sor_niter_cmp
    copy_SO%sor_niter_prm = this%sor_niter_prm
    copy_SO%sor_rho_histo_cmp = this%sor_rho_histo_cmp
    copy_SO%sor_random_obs_prm = this%sor_random_obs_prm
    copy_SO%sor_gaussian_pdf_prm = this%sor_gaussian_pdf_prm
    copy_SO%uniform_pdf_prm = this%uniform_pdf_prm

    ! VoV parameters
    IF (ASSOCIATED(this%vov_map_cmp)) THEN
       ALLOCATE(copy_SO%vov_map_cmp(SIZE(this%vov_map_cmp,dim=1), &
            SIZE(this%vov_map_cmp,dim=2)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (75).", 1)
          RETURN
       END IF
       copy_SO%vov_map_cmp = this%vov_map_cmp
    END IF
    copy_SO%vov_scaling_prm = this%vov_scaling_prm
    copy_SO%vov_scaling_cmp = this%vov_scaling_cmp
    copy_SO%vov_norb_prm = this%vov_norb_prm
    copy_SO%vov_norb_cmp = this%vov_norb_cmp
    copy_SO%vov_ntrial_prm = this%vov_ntrial_prm
    copy_SO%vov_ntrial_cmp = this%vov_ntrial_cmp
    copy_SO%vov_niter_prm = this%vov_niter_prm
    copy_SO%vov_niter_cmp = this%vov_niter_cmp
    copy_SO%vov_nmap_prm = this%vov_nmap_prm
    copy_SO%vov_mapping_mask_prm = this%vov_mapping_mask_prm
    copy_SO%vov_scaling_ready_cmp = this%vov_scaling_ready_cmp

    ! Least-squares parameters
    copy_SO%ls_corr_fac_prm = this%ls_corr_fac_prm
    copy_SO%ls_niter_major_max_prm = this%ls_niter_major_max_prm
    copy_SO%ls_niter_major_min_prm = this%ls_niter_major_min_prm
    copy_SO%ls_niter_minor_prm = this%ls_niter_minor_prm
    copy_SO%ls_elem_mask_prm = this%ls_elem_mask_prm

    copy_SO%is_initialized_prm = this%is_initialized_prm

  END FUNCTION copy_SO





  !! *Description*:
  !!
  LOGICAL FUNCTION exist_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this

    exist_SO = this%is_initialized_prm

  END FUNCTION exist_SO





  !! *Description*:
  !!
  !! Automated statistical orbital ranging. By gradually increasing
  !! the number of sample orbits (10->200->norb_init), the topocentric
  !! range intervals are iterated until they correspond to an unbiased 
  !! phase space region of possible orbits.
  !!
  !! Output: Probability density function for Cartesian position and
  !!         velocity at specified epoch in two-body approximation.
  !!
  !!
  !! Returns error.
  !!
  !! @author  JV, MG
  !! @version 2009-08-04
  !!
  SUBROUTINE autoStatisticalRanging(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this

    REAL(bp), DIMENSION(2,2) :: rho_bounds_
    REAL(bp) :: ddchi2, & ! difference to the given dchi2 limit
         pdf_ml_final, pdf_ml_
    INTEGER :: niter_, norb_prm 
    LOGICAL :: uniform_pdf_prm

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    this%sor_niter_cmp = 0
    norb_prm = this%sor_norb_prm
    uniform_pdf_prm = this%uniform_pdf_prm    

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "***************"
       WRITE(stdout,"(2X,A)") "FIRST ITERATION"
       WRITE(stdout,"(2X,A)") "***************"
    END IF

    ! Force 10 orbits and uniform pdf in the first step
    this%sor_norb_prm = 10
    this%uniform_pdf_prm = .TRUE.
    CALL statisticalRanging(this)
    IF (error .OR. this%sor_norb_cmp <= 1) THEN
       CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
            "First iteration failed.", 1)
       RETURN
    END IF
    this%sor_niter_cmp = 1
    CALL updateRanging(this, automatic=.TRUE.)

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(1X)")
       WRITE(stdout,"(2X,A)") "*****************"
       WRITE(stdout,"(2X,A)") "SECOND ITERATION "
       WRITE(stdout,"(2X,A)") "*****************"
    END IF

    ! Force 200 orbits and uniform pdf in the second step
    this%sor_norb_prm = 200
    ! Draw rho from uniform p.d.f. regardless of initial choice
    this%sor_gaussian_pdf_prm = .FALSE.
    CALL statisticalRanging(this)
    IF (error .OR. this%sor_norb_cmp <= 1) THEN
       CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
            "Second iteration failed.", 1)
       RETURN
    END IF
    this%sor_niter_cmp = 2
    CALL updateRanging(this, automatic=.TRUE.)

    ! Return to user-defined pdf mode 
    this%uniform_pdf_prm = uniform_pdf_prm

    ! Additional iteration, if requested nr of orbits very large
    IF (norb_prm/this%sor_norb_prm >= 50) THEN
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(1X)")
          WRITE(stdout,"(2X,A)") "********************"
          WRITE(stdout,"(2X,A)") "ADDITIONAL ITERATION"
          WRITE(stdout,"(2X,A)") "********************"
       END IF
       this%sor_norb_prm = 2000
       CALL statisticalRanging(this)
       IF (error .OR. this%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
               "Additional iteration failed.", 1)
          RETURN
       END IF

       this%sor_niter_cmp = 3
       CALL updateRanging(this, automatic=.TRUE.)
    END IF

    ! Final iteration with requested nr of orbits, subsequent iterations if necessary
    IF (info_verb >= 2) THEN 
       WRITE(stdout,"(1X)")
       WRITE(stdout,"(2X,A)") "****************"
       WRITE(stdout,"(2X,A)") "FINAL ITERATIONS"
       WRITE(stdout,"(2X,A)") "****************"
    END IF

    niter_ = this%sor_niter_cmp
    this%sor_norb_prm  = norb_prm
    this%sor_rho_histo_cmp = 1
    ddchi2             = 10.0_bp
    DO WHILE (this%sor_rho_histo_cmp > 0 .OR. ddchi2 > 2.0_bp)

       IF (this%sor_niter_cmp >= (this%sor_niter_prm + niter_)) THEN
          EXIT
       END IF

       CALL statisticalRanging(this)
       IF (error .OR. this%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
               "Subsequent iteration failed.", 1)
          RETURN
       END IF

       this%sor_niter_cmp = this%sor_niter_cmp + 1
       IF (this%sor_norb_cmp < this%sor_norb_prm .AND. err_verb >= 2) THEN
          WRITE(stderr,"(A,I0)") " Warning by autoStatisticalRanging:" // &
               " Number of sample orbits too small: ", this%sor_norb_cmp
          WRITE(stderr,"(1X)")
       END IF

       pdf_ml_final = MAXVAL(this%pdf_arr_cmp)
       ddchi2 = ABS(2.0_bp*LOG(pdf_ml_final/this%pdf_ml_prm))

       pdf_ml_ = this%pdf_ml_prm
       rho_bounds_ = this%sor_rho_prm
       CALL updateRanging(this, automatic = .TRUE.)
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,I0,A,F10.5)") "Histogram flag = ", &
               this%sor_rho_histo_cmp,", d(dchi2) = ", ddchi2
          WRITE(stdout,"(1X)")
       END IF

    END DO
    ! Save the last _used_ rho bounds
    this%sor_rho_prm = rho_bounds_
    this%pdf_ml_prm = pdf_ml_

  END SUBROUTINE autoStatisticalRanging





  SUBROUTINE comparePropagationParameters(this, orb, parameters_agree)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Orbit), INTENT(in) :: orb
    LOGICAL, INTENT(out) :: parameters_agree
    CHARACTER(len=DYN_MODEL_LEN) :: storb_dyn_model, orb_dyn_model
    CHARACTER(len=INTEGRATOR_LEN) :: storb_integrator, orb_integrator
    REAL(bp), DIMENSION(6) :: storb_finite_diff, orb_finite_diff
    REAL(bp) :: storb_integration_step, orb_integration_step

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / comparePropagationParameters", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    parameters_agree = .TRUE.
    CALL getParameters(this, &
         dyn_model=storb_dyn_model, &
         integration_step=storb_integration_step, &
         integrator=storb_integrator, &
         finite_diff=storb_finite_diff)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / comparePropagationParameters", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    CALL getParameters(orb, &
         dyn_model=orb_dyn_model, &
         integration_step=orb_integration_step, &
         integrator=orb_integrator, &
         finite_diff=orb_finite_diff)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / comparePropagationParameters", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    IF (storb_dyn_model == "n-body" .AND. &
         orb_dyn_model == "n-body") THEN
       IF (storb_integrator /= orb_integrator .OR. &
            storb_integration_step /= orb_integration_step) THEN
          parameters_agree = .FALSE.
          CALL infoMessage("StochasticOrbit / comparePropagationParameters", &
               "Integrators or step-sizes disagree.", stdout, 2)
       END IF
       IF (ANY(ABS(storb_finite_diff-orb_finite_diff) > &
            EPSILON(storb_finite_diff))) THEN
          parameters_agree = .FALSE.
          CALL infoMessage("StochasticOrbit / comparePropagationParameters", &
               "Finite difference values disagree.", stdout, 2)
       END IF
    ELSE IF (storb_dyn_model /= orb_dyn_model) THEN
       parameters_agree = .FALSE.
       CALL infoMessage("StochasticOrbit / comparePropagationParameters", &
            "Propagation models disagree.", stdout, 2)
    END IF

  END SUBROUTINE comparePropagationParameters





  !! *Description*:
  !!
  !! Returns .TRUE. if this StochasticOrbit instance contains an
  !! orbital-element pdf in the form of sample orbits and their
  !! weights. Otherwise returns .FALSE.
  !!
  LOGICAL FUNCTION containsSampledPDF(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this

    ! The default result is .FALSE.:
    containsSampledPDF = .FALSE.

    IF (.NOT. this%is_initialized_prm) THEN
       RETURN
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp) .AND. &
         ASSOCIATED(this%pdf_arr_cmp) .AND. &
         ASSOCIATED(this%jac_arr_cmp)) THEN
       containsSampledPDF = .TRUE.
    END IF

  END FUNCTION containsSampledPDF





  SUBROUTINE covarianceSampling(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this

    TYPE (Time) :: t0
    TYPE (Orbit), DIMENSION(:), ALLOCATABLE :: &
         orb_arr
    TYPE (Orbit) :: &
         orb_nominal, &
         orb
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
         comp_scoords, obs_scoords
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoords
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=64) :: &
         frmt = "(F20.15,1X)", &
         efrmt = "(E10.4,1X)"
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         residuals3, &
         partials_arr, &
         information_matrix_obs
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         jacobian_arr, &
         obs_coords, &
         residuals2
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         reg_apriori_arr, &
         pdf_arr, &
         rchi2_arr, &
         cosdec0
    REAL(bp), DIMENSION(6,6) :: &
         cov, &
         A, &
         information_matrix_elem, &
         jacobian_matrix
    REAL(bp), DIMENSION(6) :: &
         elements_nominal, &
         ran, &
         deviates, &
         mean, &
         elements, &
         comp_coord, &
         stdev
    REAL(bp) :: &
         apriori, &
         pdf_ml_global_ls, &
         chi2, &
         pdf_val, &
         pdf_relative_bound, &
         obs_, &
         comp_
    INTEGER, DIMENSION(:), ALLOCATABLE :: &
         failed_flag
    INTEGER, DIMENSION(6) :: &
         n0
    INTEGER :: &
         i, j, &
         itrial, &
         iorb, &
         err, &
         nobs, &
         nfailed
    LOGICAL, DIMENSION(:,:), POINTER :: &
         mask_arr

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1,6
       n0(i) = COUNT(this%obs_masks_prm(:,i))
    END DO
    orb_nominal = getNominalOrbit(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    frame = getFrame(orb_nominal)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    elements_nominal = getElements(orb_nominal, this%element_type_prm, frame)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    cov = getCovarianceMatrix(this, this%element_type_prm, frame)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF
    DO i=1,6
       stdev(i) = SQRT(cov(i,i))
    END DO
    t0 = getTime(orb_nominal)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF


    ! Observations
    information_matrix_obs => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (30)", 1)
       RETURN
    END IF
    obs_scoords => getObservationSCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (35)", 1)
       RETURN
    END IF
    nobs = SIZE(obs_scoords,dim=1)
    ALLOCATE(orb_arr(this%cos_norb_prm), cosdec0(nobs), &
         residuals2(nobs,6), residuals3(this%cos_norb_prm,nobs,6), &
         failed_flag(4), reg_apriori_arr(this%cos_norb_prm), &
         pdf_arr(this%cos_norb_prm), rchi2_arr(this%cos_norb_prm), &
         jacobian_arr(this%cos_norb_prm,3), obs_coords(nobs,6), &
         mask_arr(nobs,6), stat=err)
    DO i=1,nobs
       obs_coords(i,:) = getCoordinates(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (40)",1)
          RETURN
       END IF
       cosdec0(i) = COS(obs_coords(i,3))
    END DO
    obsy_ccoords => getObservatoryCCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (45)",1)
       RETURN
    END IF

    ! ML PDF value
    CALL getEphemerides(orb_nominal, obsy_ccoords, comp_scoords, &
         partials_arr=partials_arr)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (50)",1)
       RETURN
    END IF
    ! Multiply RA partials with cosine of observed declination:
    DO i=1,nobs
       partials_arr(2,:,i) = partials_arr(2,:,i)*cosdec0(i)
    END DO
    ! Residuals and chi2 for the global fit:
    DO i=1,nobs
       comp_coord = getCoordinates(comp_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (55)",1)
          RETURN
       END IF
       residuals2(i,1:6) = obs_coords(i,1:6) - comp_coord(1:6)
       residuals2(i,2) = residuals2(i,2) * cosdec0(i)
       IF (ABS(residuals2(i,2)) > pi) THEN
          obs_ = obs_coords(i,2)
          comp_ = comp_coord(2)
          IF (obs_ < comp_) THEN
             obs_ = obs_ + two_pi
          ELSE
             comp_ = comp_ + two_pi
          END IF
          residuals2(i,2) = (obs_ - comp_) * cosdec0(i)
       END IF
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,3"//TRIM(frmt)//")") "observed pos.", obs_coords(i,1:3)
          WRITE(stdout,"(2X,A,3"//TRIM(frmt)//")") "computed pos.", comp_coord(1:3)
       END IF
    END DO
    ! Compute chi2:
    chi2 = chi_square(residuals2, information_matrix_obs, this%obs_masks_prm, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (60)" // TRIM(errstr),1)
       RETURN
    END IF
    DEALLOCATE(residuals2)
    IF (chi2 < 0.0_bp) THEN
       WRITE(stdout,*) "Negative chi2 (", chi2, ") for ML solution."
    END IF
    ! Jeffrey's apriori:
    IF (this%regularization_prm) THEN
       ! Sigma_elements^(-1) = A^T Sigma_obs^(-1) A, where A is the
       ! partial derivatives matrix of ephemerides wrt elements:
       information_matrix_elem(:,:) = 0.0_bp
       DO i=1,nobs
          information_matrix_elem = information_matrix_elem + &
               MATMUL(MATMUL(TRANSPOSE(partials_arr(1:6,1:6,i)), &
               information_matrix_obs(i,1:6,1:6)), &
               partials_arr(1:6,1:6,i))
       END DO
       ! Jeffrey's apriori:
       apriori = SQRT(ABS(determinant(information_matrix_elem, errstr)))
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (65)", 1)
          DO i=1,nobs
             CALL matrix_print(partials_arr(1:6,1:6,i),stderr,errstr)
             WRITE(stderr,*)
             CALL matrix_print(information_matrix_obs(i,1:6,1:6),stderr,errstr)
             WRITE(stderr,*)
          END DO
          CALL matrix_print(information_matrix_elem,stderr,errstr)
          RETURN
       END IF
    ELSE
       apriori = 1.0_bp
    END IF
    pdf_ml_global_ls = apriori*EXP(-0.5_bp*(chi2 - COUNT(this%obs_masks_prm)))
    ! If initial estimate for maximum likelihood pdf remains
    ! undetermined (indicated by a negative value), set it equal to
    ! the pdf computed for the global least squares:
    IF (this%pdf_ml_prm < 0.0_bp) THEN
       this%pdf_ml_prm = pdf_ml_global_ls
    END IF

    ! From Devroye: "Non-Uniform Random Variate Generation", 1986,
    ! chapter 11
    A = 0.0_bp
    DO i=1,6
       DO j=1,i
          IF (j == 1) THEN
             A(i,1) = cov(i,1)/SQRT(cov(1,1))
          ELSE IF (i == j) THEN
             A(i,i) = SQRT(cov(i,i) - SUM(A(i,1:i-1)**2))
          ELSE IF (j < i) THEN
             A(i,j) = (cov(i,j) - SUM(A(i,1:j-1)*A(j,1:j-1)))/A(j,j)
          END IF
       END DO
    END DO

    ! Relative bound for probability density:
    pdf_relative_bound = EXP(-0.5_bp*this%dchi2_prm)
    failed_flag = 0
    iorb = 0
    itrial = 0
    mean = 0.0_bp
    IF (this%element_type_prm == "keplerian") THEN
       WRITE(stdout,"(A20,6(1X,F15.10))") "Nominal elements:", elements_nominal(1:2), elements_nominal(3:6)/rad_deg
    ELSE
       WRITE(stdout,"(A20,6(1X,F15.10))") "Nominal elements:", elements_nominal(1:6)
    END IF
    DO WHILE (iorb < this%cos_norb_prm .AND. itrial < this%cos_ntrial_prm)

       itrial = itrial + 1

       IF (itrial == 1) THEN
          orb = copy(orb_nominal)
       ELSE

          WRITE(stdout,"(1X,A,I0)") "Trial orbit #", itrial
          WRITE(stdout,"(1X,A,I0)") "Accepted orbits so far: ", iorb

          IF (this%cos_gaussian_prm) THEN
             CALL randomGaussian(ran)
          ELSE
             CALL randomNumber(ran)
             ran = 2.0_bp*ran - 1.0_bp
          END IF
          deviates = this%cos_nsigma_prm*MATMUL(A,ran) + mean
          ! New coordinates = old coordinates + deviates:
          elements = elements_nominal + deviates
          !WRITE(50,"(12(F20.15,1X))") elements, ran
          IF (info_verb >= 2) THEN
             IF (this%element_type_prm == "keplerian") THEN
                WRITE(stdout,"(A20,6(1X,F15.10))") "Element uncertainty:", stdev(1:2), stdev(3:6)/rad_deg
                WRITE(stdout,"(A20,6(1X,F15.10))") "Deviates:", deviates(1:2), deviates(3:6)/rad_deg
                IF (elements(2) < 0.0_bp .OR. elements(2) > 1.0_bp) THEN
                   ! Eccentricity out of bounds (or should non-elliptic orbits be accepted?).
                   failed_flag(1) = failed_flag(1) + 1
                   CYCLE
                END IF
                IF (elements(3) < 0.0_bp .OR. elements(3) > pi) THEN
                   ! Inclination not defined.
                   failed_flag(2) = failed_flag(2) + 1
                   CYCLE
                END IF
                ! 0 <= angle < 2pi :
                elements(4:6) = MODULO(elements(4:6), two_pi)
             ELSE
                WRITE(stdout,"(A20,6(1X,F15.10))") "Element uncertainty:", stdev(1:6)
                WRITE(stdout,"(A20,6(1X,F15.10))") "Deviates:", deviates(1:6)          
             END IF
             IF (this%element_type_prm == "keplerian") THEN
                WRITE(stdout,"(A20,6(1X,F15.10))") "Generated elements:", elements(1:2), elements(3:6)/rad_deg
             ELSE
                WRITE(stdout,"(A20,6(1X,F15.10))") "Generated elements:", elements(1:6)          
             END IF
          END IF
          IF (ANY(ABS(deviates) > this%cos_nsigma_prm*stdev)) THEN
             WRITE(stdout,*) ABS(deviates) > this%cos_nsigma_prm*stdev
             !stop
          END IF
          CALL NEW(orb, elements, this%element_type_prm, frame, t0)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "TRACE BACK (70)", 1)
             RETURN
          END IF
       END IF

       CALL setParameters(orb, dyn_model=this%dyn_model_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm, &
            perturbers=this%perturbers_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (75)",1)
          RETURN
       END IF

       !!
       !! 5) ACCEPTANCE / REJECTION OF GENERATED ORBIT
       !!
       CALL getEphemerides(orb, obsy_ccoords, comp_scoords, &
            partials_arr=partials_arr)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (80)",1)
          RETURN
       END IF
       ! Multiply RA partials with cosine of observed declination:
       DO i=1,nobs
          partials_arr(2,:,i) = partials_arr(2,:,i)*cosdec0(i)
       END DO

       ! Sky-plane residuals and chi-squares:
       DO i=1,nobs
          comp_coord = getCoordinates(comp_scoords(i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "TRACE BACK (85)",1)
             RETURN
          END IF
          residuals3(iorb+1,i,1:6) = obs_coords(i,1:6) - comp_coord(1:6)
          residuals3(iorb+1,i,2) = residuals3(iorb+1,i,2) * cosdec0(i)
          IF (ABS(residuals3(iorb+1,i,2)) > pi) THEN
             obs_ = obs_coords(i,2)
             comp_ = comp_coord(2)
             IF (obs_ < comp_) THEN
                obs_ = obs_ + two_pi
             ELSE
                comp_ = comp_ + two_pi
             END IF
             residuals3(iorb+1,i,2) = (obs_ - comp_) * cosdec0(i)
          END IF
          IF (info_verb >= 4) THEN
             WRITE(stdout,"(2X,A,3"//TRIM(frmt)//")") "observed pos.", obs_coords(i,1:3)
             WRITE(stdout,"(2X,A,3"//TRIM(frmt)//")") "computed pos.", comp_coord(1:3)
          END IF
       END DO

!!$       mask_arr = .FALSE.
!!$       WHERE (this%obs_masks_prm .AND. ABS(residuals3(iorb+1,:,:)) > this%res_accept_prm)
!!$          mask_arr = .TRUE.
!!$       END WHERE
!!$       IF (info_verb >= 4) THEN
!!$          DO i=1,nobs
!!$             WRITE(stdout,"(2X,A,2"//TRIM(frmt)//")") "O-C residuals (RA, Dec):", &
!!$                  residuals3(iorb+1,i,2:3)/rad_asec
!!$          END DO
!!$          WRITE(stdout,"(2X,A,I0,A,I0)") &
!!$               "No of omitted obs/included obs: ", &
!!$               COUNT(mask_arr),"/",n0(2)
!!$       END IF
!!$       IF (COUNT(mask_arr) > 0) THEN
!!$          ! Residuals are too large for at least one observation.
!!$          failed_flag(3) = failed_flag(3) + 1
!!$          IF (info_verb >= 5) THEN
!!$             WRITE(stdout,"(2X,A)") &
!!$                  "Failed (residuals are too large)"
!!$          END IF
!!$          CALL NULLIFY(orb)
!!$          DEALLOCATE(comp_scoords, partials_arr, stat=err)
!!$          IF (err /= 0) THEN
!!$             error = .TRUE.
!!$             DEALLOCATE(comp_scoords, stat=err)
!!$             DEALLOCATE(partials_arr, stat=err)
!!$             CALL errorMessage("StochasticOrbit / covarianceSampling", &
!!$                  "Could not deallocate memory (5)", 1)
!!$             RETURN
!!$          END IF
!!$          CYCLE
!!$       END IF

       ! Compute chi2:
       chi2 = chi_square(residuals3(iorb+1,:,:), information_matrix_obs, this%obs_masks_prm, errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (90)" // TRIM(errstr),1)
          RETURN
       END IF
       IF (chi2 < 0.0_bp) THEN
          WRITE(stdout,*) "Negative chi2 (", chi2, ") at trial ", itrial
       END IF
       ! Jeffrey's apriori:
       ! Monitor the matrix inversion?
       IF (this%regularization_prm) THEN
          ! Sigma_elements^(-1) = A^T Sigma_obs^(-1) A, where A is the
          ! partial derivatives matrix of ephemerides wrt elements:
          information_matrix_elem(:,:) = 0.0_bp
          DO i=1,nobs
             information_matrix_elem = information_matrix_elem + &
                  MATMUL(MATMUL(TRANSPOSE(partials_arr(1:6,1:6,i)), &
                  information_matrix_obs(i,1:6,1:6)), &
                  partials_arr(1:6,1:6,i))
          END DO
          ! Jeffrey's apriori:
          apriori = SQRT(ABS(determinant(information_matrix_elem, errstr)))
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "TRACE BACK (95)", 1)
             DO i=1,nobs
                CALL matrix_print(partials_arr(1:6,1:6,i),stderr,errstr)
                WRITE(stderr,*)
                CALL matrix_print(information_matrix_obs(i,1:6,1:6),stderr,errstr)
                WRITE(stderr,*)
             END DO
             CALL matrix_print(information_matrix_elem,stderr,errstr)
             RETURN
          END IF
       ELSE
          apriori = 1.0_bp
       END IF

       ! Probability density function (note that the '- nobs' term
       ! is there for practical reasons):
       pdf_val = apriori*EXP(-0.5_bp*(chi2 - COUNT(this%obs_masks_prm)))
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A)") "Sample information matrix:"
          CALL matrix_print(information_matrix_elem,stdout,errstr)
          WRITE(stdout,"(2X,A,1X,2"//TRIM(frmt)//")") "Sample chi2, rchi2:", chi2, chi2-SUM(n0(1:6))
          WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Sample apriori:", apriori
          WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Sample pdf:", pdf_val
          WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "ML pdf:", this%pdf_ml_prm
          WRITE(stdout,*) "Uniform sampling: ", this%uniform_pdf_prm 
          WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Sample PDF / ML PDF:", pdf_val / this%pdf_ml_prm
          WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Relative bound for pdf:", pdf_relative_bound        
       END IF
       IF (.NOT.this%uniform_pdf_prm .AND. &
            pdf_val / this%pdf_ml_prm < pdf_relative_bound) THEN
          ! The PDF is used and its value is not acceptable.
          failed_flag(4) = failed_flag(4) + 1
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(2X,A,1X,E10.5)") &
                  "Failed (PDF value not acceptable)", pdf_val
          END IF
          CALL NULLIFY(orb)
          DEALLOCATE(comp_scoords, partials_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             DEALLOCATE(comp_scoords, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "Could not deallocate memory (5)", 1)
             RETURN
          END IF
          CYCLE
       END IF

       iorb = iorb + 1
       WRITE(stdout,"(1X,A,I0,1X,A)") "Orbit #", iorb, "accepted."
       orb_arr(iorb) = copy(orb)
       reg_apriori_arr(iorb) = apriori
       pdf_arr(iorb) = pdf_val
!!$       if (this%pdf_ml_prm < pdf_val) then
!!$          this%pdf_ml_prm = pdf_val
!!$       end if
       rchi2_arr(iorb) = chi2 - SUM(n0(1:6))
       ! P.d.f.:
       ! - exponential part
       ! - a priori p.d.f. for invariance (optionally)
       ! ***************************** debiasing_factor_arr ****************************

       IF (this%jacobians_prm) THEN
          ! Determinant of Jacobian between topocentric spherical
          ! coordinates of the first and the last observation and
          ! orbital parameters required for output ("Topocentric Wrt
          ! Cartesian/Keplerian"):
          jacobian_matrix(1:3,:) = partials_arr(1:3,:,1) / &
               cosdec0(1)
          jacobian_matrix(4:6,:) = partials_arr(1:3,:,nobs) / &
               cosdec0(nobs)
          jacobian_arr(iorb,1) = ABS(determinant(jacobian_matrix, errstr))
          IF (LEN_TRIM(errstr) /= 0) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "Unsuccessful computation of determinant of orbital element " // &
                  "jacobian matrix " // TRIM(errstr), 1)
             errstr = ""
             IF (err_verb >= 1) THEN
                CALL matrix_print(jacobian_matrix, stderr, errstr)
             END IF
             errstr = ""
             CYCLE
          END IF

          ! Determinant of Jacobian between Cartesian and Keplerian
          ! orbital elements ("Cartesian Wrt Keplerian"):
          CALL partialsCartesianWrtKeplerian(orb, &
               jacobian_matrix, "equatorial")
          jacobian_arr(iorb,2) = ABS(determinant(jacobian_matrix, errstr)) 
          IF (LEN_TRIM(errstr) /= 0) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "Unsuccessful computation of " // &
                  "jacobian matrix " // TRIM(errstr), 1)
             errstr = ""
             IF (err_verb >= 1) THEN
                CALL matrix_print(jacobian_matrix, stderr, errstr)
             END IF
             errstr = ""
             CALL NULLIFY(orb)
             RETURN
          END IF

          ! Determinant of Jacobian between equinoctial and
          ! Keplerian orbital elements ("Equinoctial Wrt
          ! Keplerian"):
          elements = getElements(orb, "keplerian")
          jacobian_arr(iorb,3) = 0.5_bp*elements(2) * &
               SIN(0.5_bp*elements(3)) / COS(0.5_bp*elements(3))**3
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(2X,A,F15.12)") "Chi2 new: ", chi2
             WRITE(stdout,"(2X,A,3(1X,F10.5))") "PDF values:", pdf_val, &
                  pdf_val/this%pdf_ml_prm, pdf_relative_bound
          END IF
       ELSE
          jacobian_arr(iorb,:) = -1.0_bp
       END IF

       CALL NULLIFY(orb)
       DEALLOCATE(comp_scoords, partials_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "Could not deallocate memory (15).", 1)
          RETURN
       END IF

    END DO

    IF (info_verb >= 1) THEN
       WRITE(stdout,*) "Final number of orbits and the required trials:"
       WRITE(stdout,"(2(I0,2X))") iorb, itrial
       WRITE(stdout,*) "Total failure percentage (1), and failure " // &
            "due to (2) eccentricity, (3) inclination, " // &
            "(4) residuals, and (5) pdf:"
       nfailed = SUM(failed_flag)
       nfailed = MAX(1,nfailed)
       WRITE(stdout,"(5"//TRIM(frmt)//")") &
            100.0_bp*REAL(SUM(failed_flag),bp)/itrial, &
            100.0_bp*REAL(failed_flag(1),bp)/nfailed, &
            100.0_bp*REAL(failed_flag(2),bp)/nfailed, &
            100.0_bp*REAL(failed_flag(3),bp)/nfailed, &
            100.0_bp*REAL(failed_flag(4),bp)/nfailed
    END IF

    ALLOCATE(this%reg_apr_arr_cmp(iorb), this%res_arr_cmp(iorb,nobs,6))


    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DO i=1,SIZE(this%orb_arr_cmp)
          CALL NULLIFY(this%orb_arr_cmp(i))
       END DO
       DEALLOCATE(this%orb_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "Could not deallocate memory (20).", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(this%orb_arr_cmp(iorb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "Could not allocate memory (5).", 1)
       RETURN
    END IF
    DO i=1,iorb
       this%orb_arr_cmp(i) = copy(orb_arr(i))
    END DO
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "Could not deallocate memory (25).", 1)
          RETURN       
       END IF
    END IF
    ALLOCATE(this%pdf_arr_cmp(iorb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "Could not allocate memory (10).", 1)
       RETURN
    END IF
    this%pdf_arr_cmp = pdf_arr
    IF (ASSOCIATED(this%jac_arr_cmp)) THEN
       DEALLOCATE(this%jac_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "Could not deallocate memory (25).", 1)
          RETURN       
       END IF
    END IF
    ALLOCATE(this%jac_arr_cmp(iorb,3), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "Could not allocate memory (10).", 1)
       RETURN
    END IF
    this%jac_arr_cmp = jacobian_arr
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "Could not deallocate memory (25).", 1)
          RETURN       
       END IF
    END IF
    ALLOCATE(this%rchi2_arr_cmp(iorb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "Could not allocate memory (10).", 1)
       RETURN
    END IF
    this%rchi2_arr_cmp = rchi2_arr
    DO i=1,SIZE(orb_arr)
       CALL NULLIFY(orb_arr(i))
    END DO
    this%reg_apr_arr_cmp = reg_apriori_arr
    this%res_arr_cmp = residuals3

    CALL NULLIFY(orb_nominal)

    DEALLOCATE(orb_arr, failed_flag, pdf_arr, &
         information_matrix_obs, jacobian_arr, residuals3, stat=err)

    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "Could not deallocate memory (30).", 1)
       RETURN       
    END IF

  END SUBROUTINE covarianceSampling






  !! *Description*:
  !!
  !! Returns a discrete sampling of the apoapsis distance
  !! probability-density function for this orbit.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getApoapsisDistance_SO(this, Q)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), DIMENSION(:,:), POINTER     :: Q

    TYPE (Orbit) :: orb
    REAL(bp), DIMENSION(6,6) :: partials, covariance1, covariance2
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    REAL(bp) :: jac
    INTEGER :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getApoapsisDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (containsSampledPDF(this)) THEN

       pdf_arr => getPDFValues(this, "keplerian")
       IF (error) THEN
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF

       ALLOCATE(Q(SIZE(this%orb_arr_cmp,dim=1),2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "Could not allocate memory (5).", 1)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
       partials = identity_matrix(6)
       DO i=1,SIZE(this%orb_arr_cmp,dim=1)
          CALL getApoapsisDistance(this%orb_arr_cmp(i), &
               Q(i,1), partials=partials(1,1:6))
          IF (error) THEN
             CALL errorMessage("Orbit / getApoapsisDistance", &
                  "TRACE BACK (5)", 1)
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          jac = ABS(determinant(partials, errstr))
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getApoapsisDistance", &
                  "Could not compute determinant for Jacobian. " // &
                  TRIM(errstr), 1)
             errstr = ""
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          Q(i,2) = pdf_arr(i)*jac
       END DO

       DEALLOCATE(pdf_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getApoapsisDistance", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF

    ELSE

       ALLOCATE(Q(1,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "Could not allocate memory (5).", 1)
          RETURN
       END IF
       orb = copy(this%orb_ml_cmp)
       CALL toKeplerian(orb)
       partials = identity_matrix(6)
       CALL getApoapsisDistance(orb, &
            Q(1,1), &
            partials=partials(1,1:6))
       IF (error) THEN
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       CALL NULLIFY(orb)

       covariance1 = getCovarianceMatrix(this, "keplerian")
       covariance2 = MATMUL(MATMUL(partials,covariance1),TRANSPOSE(partials))
       Q(1,2) = SQRT(covariance2(1,1))

    END IF

  END SUBROUTINE getApoapsisDistance_SO





  SUBROUTINE getAPrioriWeights(this, apriori, apriori_pdf)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)   :: this
    REAL(bp), DIMENSION(:,:), INTENT(in) :: apriori
    REAL(bp), DIMENSION(:), POINTER      :: apriori_pdf
    REAL(bp), DIMENSION(6)               :: elements
    INTEGER                              :: i, j, norb, nbin, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getAPrioriWeights", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nbin = SIZE(apriori,dim=1)
    IF (.NOT. ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getAPrioriWeights", &
            "No sample orbits available.", 1)
       RETURN
    END IF
    norb = SIZE(this%orb_arr_cmp,dim=1)
    ALLOCATE(apriori_pdf(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getAPrioriWeights", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    apriori_pdf = 0.0_bp

    DO i=1,norb
       elements = getElements(this%orb_arr_cmp(i), "keplerian")
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getAPrioriWeights", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       elements(3) = elements(3)/rad_deg
       DO j=1,nbin
          IF (elements(1) >= apriori(j,1) .AND. elements(1) < apriori(j,2) .AND. &
               elements(2) >= apriori(j,3) .AND. elements(2) < apriori(j,4) .AND. &
               elements(3) >= apriori(j,5) .AND. elements(3) < apriori(j,6)) THEN
             apriori_pdf(i) = apriori(j,7)
             EXIT
          END IF
       END DO
    END DO

  END SUBROUTINE getAPrioriWeights





  !! Description:
  !!
  !! Returns the sample orbit with the highest PDF value. Optionally,
  !! the PDF of all sample orbits is multiplied with a given priori
  !! distribution before choosing the highest value.
  !!
  FUNCTION getBestFittingSampleOrbit(this, apriori)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this
    REAL(bp), DIMENSION(:,:), OPTIONAL  :: apriori
    TYPE (Orbit)                        :: getBestFittingSampleOrbit
    REAL(bp), DIMENSION(:), POINTER     :: apriori_pdf
    REAL(bp), DIMENSION(:), ALLOCATABLE :: pdf
    INTEGER                             :: norb, err, indx 

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    norb = SIZE(this%orb_arr_cmp,dim=1)
    ALLOCATE(pdf(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    pdf = this%pdf_arr_cmp
    pdf = pdf / SUM(pdf)
    IF (PRESENT(apriori)) THEN
       CALL getAPrioriWeights(this, apriori, apriori_pdf)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "TRACE BACK (5).", 1)
          DEALLOCATE(pdf, stat=err)
          IF (err /= 0) CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "Could not deallocate memory (5).", 1)
          RETURN
       ELSE IF (SIZE(apriori_pdf,dim=1) /= norb) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "Size of a priori array not consistent with number of orbits.", 1)
          DEALLOCATE(pdf, apriori_pdf, stat=err)
          IF (err /= 0) CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
       pdf(1:norb) = pdf(1:norb)*apriori_pdf(1:norb)
       IF (SUM(pdf) < EPSILON(pdf(1))) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "Use of priori results in zero probability.", 1)
          DEALLOCATE(pdf, apriori_pdf, stat=err)
          IF (err /= 0) CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "Could not deallocate memory (10).", 1)
          RETURN
       ELSE
          pdf = pdf / SUM(pdf)
       END IF
       DEALLOCATE(apriori_pdf, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "Could not deallocate memory (15).", 1)
          RETURN
       END IF
    END IF
    indx = MAXLOC(pdf,dim=1)
    DEALLOCATE(pdf, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
            "Could not deallocate memory (20).", 1)
       RETURN
    END IF
    getBestFittingSampleOrbit = copy(this%orb_arr_cmp(indx))

  END FUNCTION getBestFittingSampleOrbit





  REAL(bp) FUNCTION getChi2_this_orb(this, orb, obs_masks)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Orbit), INTENT(in) :: orb
    LOGICAL, DIMENSION(:,:), INTENT(in), OPTIONAL :: obs_masks

    REAL(bp), DIMENSION(:,:,:), POINTER :: information_matrix
    REAL(bp), DIMENSION(:,:), POINTER :: residuals
    INTEGER :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getChi2", &
            "Object has not been initialized.", 1)
       RETURN
    END IF
    IF (.NOT.exist(orb)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getChi2", &
            "Orbit has not been initialized.", 1)
       RETURN
    END IF
    residuals => getResiduals(this, orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getChi2", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(residuals, stat=err)
       RETURN
    END IF
    information_matrix => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getChi2", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(information_matrix, stat=err)
       RETURN
    END IF
    IF (PRESENT(obs_masks)) THEN
       getChi2_this_orb = chi_square(residuals, information_matrix, obs_masks, errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getChi2", &
               "TRACE BACK (15) " // TRIM(errstr), 1)
          errstr = ""
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(information_matrix, stat=err)
          RETURN
       END IF
    ELSE
       getChi2_this_orb = chi_square(residuals, information_matrix, this%obs_masks_prm, errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getChi2", &
               "TRACE BACK (25) " // TRIM(errstr), 1)
          errstr = ""
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(information_matrix, stat=err)
          RETURN
       END IF
    END IF
    DEALLOCATE(residuals, information_matrix, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getChi2", &
            "Could not deallocate memory (10).", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(information_matrix, stat=err)
       RETURN
    END IF

  END FUNCTION getChi2_this_orb





  !! *Description*:
  !!
  !! Returns chi2 based on residuals and information matrix. Has been tested.
  !!
  !! @author MG, JV
  !! @version 22.5.2006
  !!
  FUNCTION getChi2_matrix(residuals, information_matrix, mask)

    IMPLICIT NONE
    REAL(bp), DIMENSION(:,:), INTENT(in)          :: residuals
    REAL(bp), DIMENSION(:,:), INTENT(in)          :: information_matrix
    LOGICAL, DIMENSION(:,:), INTENT(in), OPTIONAL :: mask
    REAL(bp)                                      :: getChi2_matrix
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: residuals_
    REAL(bp), DIMENSION(:), ALLOCATABLE :: residual_vector, tmp
    INTEGER :: i, j, nobs, nmulti,err

    ! chi2 = residuals*information_matrix*transpose(residuals), where
    ! "residuals" is (nmulti*nobs x 1) and "information_matrix" is
    ! (nmulti*nobs x nmulti*nobs):

    nobs = SIZE(residuals,dim=1)
    nmulti = SIZE(residuals,dim=2)
    ALLOCATE(residuals_(nobs,nmulti), residual_vector(nobs*nmulti), &
         tmp(nobs*nmulti), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getChi2", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    residuals_ = residuals
    IF (PRESENT(mask)) THEN
       WHERE (.NOT. mask)
          residuals_ = 0.0_bp
       END WHERE
    END IF

    DO i=1,nobs
       j = (i-1)*nmulti
       residual_vector(j+1:j+nmulti) = residuals_(i,:)
    END DO
    DO i=1,nobs*nmulti
       tmp(i) = DOT_PRODUCT(residual_vector(:), information_matrix(:,i))
    END DO
    getChi2_matrix = DOT_PRODUCT(tmp(:), residual_vector(:))

    DEALLOCATE(residuals_, residual_vector, tmp, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getChi2", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION getChi2_matrix





  !! *Description*:
  !!
  !! Returns orbital element covariance matrix for a given orbital
  !! element type (Keplerian/Cartesian).
  !!
  !! Tested:
  !!
  !! cov_type="keplerian" has been tested using Keplerian ls and
  !! Cartesian ls + partials, where partials are partials
  !! Keplerian wrt Cartesian. Relative errors of sqrt(diagonal
  !! elements) are  ~10^-8 sigma.
  !!
  !! cov_type="cartesian" has been tested using Keplerian ls +
  !! inv(partials) and Cartesian ls, where partials are partials
  !! Keplerian wrt Cartesian. Relative errors of sqrt(diagonal
  !! elements) are ~10^-7 sigma.
  !!
  !! Returns error.
  !!
  FUNCTION getCovarianceMatrix_SO(this, cov_type, frame)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)     :: this
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cov_type
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame
    REAL(bp), DIMENSION(6,6)               :: getCovarianceMatrix_SO

    CHARACTER(len=ELEMENT_TYPE_LEN)        :: cov_type_
    CHARACTER(len=FRAME_LEN)               :: frame_
    REAL(bp), DIMENSION(6,6)               :: partials

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. ASSOCIATED(this%cov_ml_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
            "Covariances are not available.", 1)
       RETURN
    END IF

    IF (PRESENT(cov_type)) THEN
       cov_type_ = cov_type
    ELSE
       cov_type_ = this%cov_type_prm
    END IF
    CALL locase(cov_type_, error)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
            "The element type string contains forbidden characters.", 1)
       RETURN
    END IF

    IF (PRESENT(frame)) THEN
       frame_ = frame
    ELSE IF (cov_type_ == "keplerian" .OR. cov_type_ == "cometary") THEN
       frame_ = "ecliptic"
    ELSE IF (.NOT.PRESENT(frame) .AND. .NOT.PRESENT(cov_type)) THEN
       frame_ = getFrame(this%orb_ml_cmp)
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
            "Frame must be given for non-Keplerian elements.", 1)       
       RETURN
    END IF

    IF (cov_type_ == this%cov_type_prm) THEN
       getCovarianceMatrix_SO = this%cov_ml_cmp
       RETURN
    ELSE IF (this%cov_type_prm == "cartesian" .AND. &
         cov_type_ == "cometary") THEN
       CALL partialsCometaryWrtCartesian(this%orb_ml_cmp, partials)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    ELSE IF (this%cov_type_prm == "cartesian" .AND. &
         cov_type_ == "keplerian") THEN
       CALL partialsKeplerianWrtCartesian(this%orb_ml_cmp, partials)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    ELSE IF (this%cov_type_prm == "cometary" .AND. &
         cov_type_ == "cartesian") THEN
       CALL partialsCartesianWrtCometary(this%orb_ml_cmp, partials, frame_)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    ELSE IF (this%cov_type_prm == "cometary" .AND. &
         cov_type_ == "keplerian") THEN
       CALL partialsKeplerianWrtCometary(this%orb_ml_cmp, partials)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    ELSE IF (this%cov_type_prm == "keplerian" .AND. &
         cov_type_ == "cartesian") THEN
       CALL partialsCartesianWrtKeplerian(this%orb_ml_cmp, partials, frame_)
    ELSE IF (this%cov_type_prm == "keplerian" .AND. &
         cov_type_ == "cometary") THEN
       CALL partialsCometaryWrtKeplerian(this%orb_ml_cmp, partials)
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getCovarianceMatrix", &
            "No such option: " // TRIM(this%cov_type_prm) // &
            " in and " // TRIM(cov_type_) // " out.", 1)
       RETURN
    END IF
    getCovarianceMatrix_SO = MATMUL(MATMUL(partials, this%cov_ml_cmp), TRANSPOSE(partials))

  END FUNCTION getCovarianceMatrix_SO





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE getEphemerides_SO(this, observers, ephemerides_arr, &
       lt_corr, cov_arr, pdfs_arr, this_lt_corr_arr)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)                 :: this
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: observers
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER  :: ephemerides_arr
    LOGICAL, INTENT(in), OPTIONAL                         :: lt_corr
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL         :: cov_arr
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL           :: pdfs_arr
    TYPE (Orbit), DIMENSION(:,:), POINTER, OPTIONAL       :: this_lt_corr_arr

    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: ephemerides_arr_
    TYPE (Orbit), DIMENSION(:), POINTER                :: this_lt_corr_arr_
    REAL(bp), DIMENSION(:,:,:,:), POINTER              :: partials_arr4 
    REAL(bp), DIMENSION(:,:,:), POINTER                :: partials_arr3
    REAL(bp), DIMENSION(6,6)                           :: cov_elm
    REAL(bp)                                           :: det
    INTEGER                                            :: i, j, err, norb, nobs
    LOGICAL                                            :: lt_corr_

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getEphemeris", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nobs = SIZE(observers,dim=1)
    DO i=1,nobs
       IF (.NOT. exist(observers(i))) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "One or more of the observers has not been initialized.", 1)
          RETURN
       END IF
    END DO

    IF (PRESENT(lt_corr)) THEN
       lt_corr_ = lt_corr
    ELSE
       lt_corr_ = .TRUE.
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp) .AND. PRESENT(pdfs_arr)) THEN
       ! Sampled p.d.f.:
       IF (PRESENT(this_lt_corr_arr)) THEN
          CALL getEphemerides(this%orb_arr_cmp, observers, &
               ephemerides_arr, lt_corr=lt_corr_, &
               partials_arr=partials_arr4, &
               this_lt_corr_arr=this_lt_corr_arr)
       ELSE
          CALL getEphemerides(this%orb_arr_cmp, observers, &
               ephemerides_arr, lt_corr=lt_corr_, &
               partials_arr=partials_arr4)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(partials_arr4, stat=err)
          RETURN
       END IF
       norb = SIZE(this%orb_arr_cmp, dim=1)
       ALLOCATE(pdfs_arr(norb,nobs), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not allocate memory (5).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(partials_arr4, stat=err)
          DEALLOCATE(pdfs_arr, stat=err)
          RETURN
       END IF
       DO i=1,norb
          DO j=1,nobs
             det = determinant(partials_arr4(i,:,:,j), errstr)
             IF (LEN_TRIM(errstr) /= 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / getEphemerides", &
                     "Computation of determinant of partials between " // &
                     "ephemeris and orbital elements failed. " // &
                     TRIM(errstr), 1)
                errstr = ""
                IF (err_verb >= 1) THEN
                   CALL matrix_print(partials_arr4(i,:,:,j),stderr, errstr)
                END IF
                errstr = ""
                DEALLOCATE(ephemerides_arr, stat=err)
                DEALLOCATE(partials_arr4, stat=err)
                DEALLOCATE(pdfs_arr, stat=err)
                RETURN
             END IF
             ! Changed 2008-12-14
             !pdfs_arr(i,j) = this%pdf_arr_cmp(i)/ABS(det)
             pdfs_arr(i,j) = this%pdf_arr_cmp(i) * ABS(det)
          END DO
       END DO
       DEALLOCATE(partials_arr4, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not deallocate memory (5).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(pdfs_arr, stat=err)
          RETURN
       END IF
    ELSE IF (exist(this%orb_ml_cmp) .AND. PRESENT(cov_arr)) THEN
       ! Least squares solution + covariance matrix:
       IF (getElementType(this%orb_ml_cmp) /= &
            this%element_type_prm) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Element types for ML Orbit and StochasticOrbit are not compatible.", 1)
          DEALLOCATE(ephemerides_arr_, stat=err)
          DEALLOCATE(partials_arr3, stat=err)
          RETURN
       END IF
       IF (getElementType(this%orb_ml_cmp) /= &
            this%cov_type_prm) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Element types for ML orbit and covariance are not compatible.", 1)
          DEALLOCATE(ephemerides_arr_, stat=err)
          DEALLOCATE(partials_arr3, stat=err)
          RETURN
       END IF
       IF (PRESENT(this_lt_corr_arr)) THEN
          CALL getEphemerides(this%orb_ml_cmp, observers, &
               ephemerides_arr_, lt_corr=lt_corr_, &
               partials_arr=partials_arr3, &
               this_lt_corr_arr=this_lt_corr_arr_)
          ALLOCATE(this_lt_corr_arr(1,nobs))
          DO i=1,nobs
             this_lt_corr_arr(1,i) = copy(this_lt_corr_arr_(i))
          END DO
          DEALLOCATE(this_lt_corr_arr_)
       ELSE
          CALL getEphemerides(this%orb_ml_cmp, observers, ephemerides_arr_, &
               lt_corr=lt_corr_, partials_arr=partials_arr3)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "TRACE BACK (10)", 1)
          DEALLOCATE(ephemerides_arr_, stat=err)
          DEALLOCATE(partials_arr3, stat=err)
          RETURN
       END IF
       cov_elm = getCovarianceMatrix(this, &
            getElementType(this%orb_ml_cmp), &
            getFrame(this%orb_ml_cmp))
       ALLOCATE(ephemerides_arr(1,nobs), cov_arr(6,6,nobs), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(ephemerides_arr_, stat=err)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(cov_arr, stat=err)
          DEALLOCATE(partials_arr3, stat=err)
          RETURN
       END IF
       DO i=1,nobs
          ephemerides_arr(1,i) = copy(ephemerides_arr_(i))
          CALL NULLIFY(ephemerides_arr_(i))
          cov_arr(:,:,i) = MATMUL(MATMUL(partials_arr3(:,:,i), cov_elm), &
               TRANSPOSE(partials_arr3(:,:,i)))
       END DO
       DEALLOCATE(partials_arr3, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not deallocate memory (10).", 1)
          DEALLOCATE(ephemerides_arr_, stat=err)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(cov_arr, stat=err)
          RETURN
       END IF
       DEALLOCATE(ephemerides_arr_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not deallocate memory (15).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(cov_arr, stat=err)
          RETURN
       END IF
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getEphemerides", &
            "pdfs_arr or cov_arr must be specified.", 1)
       RETURN
    END IF

  END SUBROUTINE getEphemerides_SO





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE getEphemeris_SO(this, observer, ephemeris_arr, &
       lt_corr, cov, pdf_arr, this_lt_corr_arr, &
       cov_lt_corr, pdf_lt_corr_arr)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)              :: this
    TYPE (CartesianCoordinates), INTENT(in)            :: observer
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: ephemeris_arr
    LOGICAL, INTENT(in), OPTIONAL                      :: lt_corr
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL    :: cov
    REAL(bp), DIMENSION(:), POINTER, OPTIONAL          :: pdf_arr
    TYPE (Orbit), DIMENSION(:), POINTER, OPTIONAL      :: this_lt_corr_arr
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL    :: cov_lt_corr
    REAL(bp), DIMENSION(:), POINTER, OPTIONAL          :: pdf_lt_corr_arr

    REAL(bp), DIMENSION(:,:,:), POINTER :: partials_arr, &
         jacobian_lt_corr_arr, jacobian_prop_arr
    REAL(bp), DIMENSION(6,6)            :: cov_elm, partials, &
         jacobian_lt_corr, jacobian_prop, jacobian
    REAL(bp)                            :: det
    INTEGER                             :: i, err, norb
    LOGICAL                             :: lt_corr_

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getEphemeris", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getEphemeris", &
            "'observer' has not been initialized.", 1)
       RETURN
    END IF

    ! Light-time correction is applied by default:
    IF (PRESENT(lt_corr)) THEN
       lt_corr_ = lt_corr
    ELSE
       lt_corr_ = .TRUE.
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp) .AND. PRESENT(pdf_arr)) THEN
       ! Sampled p.d.f.:
       IF (PRESENT(this_lt_corr_arr) .AND. PRESENT(pdf_lt_corr_arr)) THEN
          CALL getEphemeris(this%orb_arr_cmp, observer, ephemeris_arr, &
               lt_corr=lt_corr_, partials_arr=partials_arr, &
               this_lt_corr_arr=this_lt_corr_arr, &
               jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
               jacobian_prop_arr=jacobian_prop_arr)
       ELSE IF (.NOT.(PRESENT(this_lt_corr_arr) .OR. PRESENT(pdf_lt_corr_arr))) THEN
          CALL getEphemeris(this%orb_arr_cmp, observer, ephemeris_arr, &
               lt_corr=lt_corr_, partials_arr=partials_arr)
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "Are you sure you want to do this?", 1)
          RETURN          
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "TRACE BACK", 1)
          DEALLOCATE(ephemeris_arr, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          DEALLOCATE(jacobian_prop_arr, stat=err)
          RETURN
       END IF
       norb = SIZE(this%orb_arr_cmp, dim=1)
       ALLOCATE(pdf_arr(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "Could not allocate memory.", 1)
          DEALLOCATE(ephemeris_arr, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          DEALLOCATE(jacobian_prop_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
       DO i=1,norb
          det = determinant(partials_arr(i,:,:), errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getEphemeris", &
                  "Computation of determinant of partials between " // &
                  "ephemeris and orbital elements failed. " // &
                  TRIM(errstr), 1)
             errstr = ""
             DEALLOCATE(ephemeris_arr, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             IF (PRESENT(this_lt_corr_arr)) THEN
                DEALLOCATE(this_lt_corr_arr, stat=err)
             END IF
             DEALLOCATE(jacobian_lt_corr_arr, stat=err)
             DEALLOCATE(jacobian_prop_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          ! Changed 2008-12-14
          !pdf_arr(i) = this%pdf_arr_cmp(i)/ABS(det)
          pdf_arr(i) = this%pdf_arr_cmp(i) * ABS(det)
       END DO
       DEALLOCATE(partials_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "Could not deallocate memory.", 1)
          DEALLOCATE(ephemeris_arr, stat=err)
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          DEALLOCATE(jacobian_prop_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
       IF (PRESENT(this_lt_corr_arr) .AND. PRESENT(pdf_lt_corr_arr)) THEN
          ALLOCATE(pdf_lt_corr_arr(norb), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getEphemeris", &
                  "Could not allocate memory.", 1)
             DEALLOCATE(ephemeris_arr, stat=err)
             DEALLOCATE(this_lt_corr_arr, stat=err)
             DEALLOCATE(jacobian_lt_corr_arr, stat=err)
             DEALLOCATE(jacobian_prop_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          DO i=1,norb
             det = determinant(MATMUL(jacobian_lt_corr_arr(i,:,:), &
                  jacobian_prop_arr(i,:,:)), errstr)
             IF (LEN_TRIM(errstr) /= 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / getEphemeris", &
                     "Computation of determinant of partials between " // &
                     "orbital elements at different epochs failed. " // &
                     TRIM(errstr), 1)
                errstr = ""
                DEALLOCATE(ephemeris_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr, stat=err)
                DEALLOCATE(jacobian_lt_corr_arr, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                RETURN
             END IF
             ! Changed 2008-12-14
             !pdf_lt_corr_arr(i) = this%pdf_arr_cmp(i)/ABS(det)
             pdf_lt_corr_arr(i) = this%pdf_arr_cmp(i) * ABS(det)
          END DO
          DEALLOCATE(jacobian_lt_corr_arr, jacobian_prop_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getEphemeris", &
                  "Could not deallocate memory.", 1)
             DEALLOCATE(ephemeris_arr, stat=err)
             DEALLOCATE(this_lt_corr_arr, stat=err)
             DEALLOCATE(jacobian_lt_corr_arr, stat=err)
             DEALLOCATE(jacobian_prop_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
       END IF
    ELSE IF (exist(this%orb_ml_cmp) .AND. PRESENT(cov)) THEN
       ! Least squares solution + covariance matrix:
       ALLOCATE(ephemeris_arr(1), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "Could not allocate memory.", 1)
          DEALLOCATE(ephemeris_arr, stat=err)
          RETURN
       END IF
       cov_elm = getCovarianceMatrix(this, &
            getElementType(this%orb_ml_cmp), &
            getFrame(this%orb_ml_cmp))
       IF (PRESENT(this_lt_corr_arr) .AND. PRESENT(cov_lt_corr)) THEN
          ALLOCATE(this_lt_corr_arr(1), stat=err)
          CALL getEphemeris(this%orb_ml_cmp, observer, ephemeris_arr(1), &
               lt_corr=lt_corr_, partials=partials, &
               this_lt_corr=this_lt_corr_arr(1), &
               jacobian_lt_corr=jacobian_lt_corr, &
               jacobian_prop=jacobian_prop)
          jacobian = MATMUL(jacobian_lt_corr,jacobian_prop)
          cov_lt_corr = MATMUL(MATMUL(jacobian, cov_elm), &
               TRANSPOSE(jacobian))
       ELSE IF (.NOT.(PRESENT(this_lt_corr_arr) .OR. PRESENT(pdf_lt_corr_arr))) THEN
          CALL getEphemeris(this%orb_ml_cmp, observer, ephemeris_arr(1), &
               lt_corr=lt_corr_, partials=partials)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "TRACE BACK", 1)
          DEALLOCATE(ephemeris_arr, stat=err)
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          RETURN
       END IF
       cov = MATMUL(MATMUL(partials, cov_elm), &
            TRANSPOSE(partials))
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getEphemeris", &
            "Ephemeris cannot be requested without proper uncertainty estimate.", 1)
       RETURN
    END IF

  END SUBROUTINE getEphemeris_SO






  !! *Description*:
  !!
  !! Returns error.
  !! 
  SUBROUTINE getGroupWeights(this, weights, groups, apriori_pdf)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)           :: this
    REAL(bp), DIMENSION(:), POINTER              :: weights
    CHARACTER(len=*), DIMENSION(:), POINTER      :: groups
    REAL(bp), DIMENSION(:), INTENT(in), OPTIONAL :: apriori_pdf

    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: elements
    REAL(bp), DIMENSION(:), ALLOCATABLE :: periapsis, apoapsis, pdf
    LOGICAL, DIMENSION(:), POINTER :: mask_array, mask_array_tot
    INTEGER :: norb, err, i, ngroup

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ngroup = 17
    IF (.NOT. ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "No sample orbits available.", 1)
       RETURN
    END IF
    norb = SIZE(this%orb_arr_cmp,dim=1)
    ALLOCATE(weights(ngroup), groups(ngroup), elements(norb,6), &
         apoapsis(norb), periapsis(norb), mask_array(norb), &
         mask_array_tot(norb), pdf(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    mask_array_tot = .TRUE. ! which are left?
    pdf = this%pdf_arr_cmp
    pdf = pdf / SUM(pdf)
    DO i=1,norb
       elements(i,:) = getElements(this%orb_arr_cmp(i), "keplerian")
    END DO
    apoapsis = elements(:,1) * (1.0_bp + elements(:,2))
    periapsis = elements(:,1) * (1.0_bp - elements(:,2))
    IF (PRESENT(apriori_pdf)) THEN
       IF (SIZE(apriori_pdf,dim=1) /= norb) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getGroupWeights", &
               "Size of a priori array not consistent with number of orbits.", 1)
          RETURN
       END IF
       pdf(1:norb) = pdf(1:norb)*apriori_pdf(1:norb)
       IF (SUM(pdf) < EPSILON(pdf(1))) THEN
          pdf = 0.0_bp
       ELSE
          pdf = pdf / SUM(pdf)
       END IF
    END IF

    ! Atira (AKA IEO AKA Apohele)
    groups(1) = "Atira"
    mask_array = .FALSE.
    WHERE (apoapsis <= 0.983_bp)
       mask_array = .TRUE.
    END WHERE
    weights(1) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Aten
    groups(2) = "NEO/Aten"
    mask_array = .FALSE.
    WHERE (elements(:,1) <= 1.0_bp .AND. &
         apoapsis > 0.983_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(2) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Apollo
    groups(3) = "NEO/Apollo"
    mask_array = .FALSE.
    WHERE (elements(:,1) > 1.0_bp .AND. &
         periapsis <= 1.017_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(3) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Amor
    groups(4) = "NEO/Amor"
    mask_array = .FALSE.
    WHERE (periapsis > 1.017_bp .AND. &
         periapsis <= 1.3_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(4) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Aethra (Mars crosser)
    groups(5) = "Aethra"
    mask_array = .FALSE.
    WHERE (periapsis > 1.3_bp .AND. &
         periapsis <= 5.0_bp/3.0_bp .AND. & ! 1.6666... AU
         mask_array_tot) mask_array = .TRUE.
    weights(5) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Hungaria
    groups(6) = "Hungaria"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 1.8_bp .AND. &
         elements(:,1) <= 2.0_bp .AND. &
         elements(:,2) >= 0.0_bp .AND. &
         elements(:,2) <= 0.16_bp .AND. &
         elements(:,3) >= 15.0_bp .AND. &
         elements(:,3) <= 35.0_bp .AND. &
         periapsis > 5.0_bp/3.0_bp .AND. & ! 1.6666... AU
         mask_array_tot) mask_array = .TRUE.
    weights(6) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Phocaea (part of the mainbelt)
    groups(7) = "MBO/Phocaea"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 2.26_bp .AND. &
         elements(:,1) <= 2.50_bp .AND. &
         elements(:,2) >= 0.1_bp .AND. &
         elements(:,2) <= 0.3_bp .AND. &
         elements(:,3) >= 10.0_bp .AND. &
         elements(:,3) <= 30.0_bp .AND. &
         periapsis > 5.0_bp/3.0_bp .AND. & ! 1.6666... AU
         mask_array_tot) mask_array = .TRUE.
    weights(7) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Mainbelt (without Phocaea)
    groups(8) = "MBO exc. Phocaea"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 2.1_bp .AND. &
         elements(:,1) <= 3.5_bp .AND. &
         elements(:,2) >= 0.0_bp .AND. &
         elements(:,2) <= 0.35_bp .AND. &
         elements(:,3) >= 0.0_bp .AND. &
         elements(:,3) <= 35.0_bp .AND. &
         periapsis > 5.0_bp/3.0_bp .AND. & ! 1.6666... AU
         periapsis <= 4.95_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(8) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Hilda
    groups(9) = "Hilda"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 3.75_bp .AND. &
         elements(:,1) <= 4.02_bp .AND. &
         elements(:,2) >= 0.0_bp .AND. &
         elements(:,2) <= 0.26_bp .AND. &
         elements(:,3) >= 0.0_bp .AND. &
         elements(:,3) <= 18.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(9) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Jupiter Trojan
    groups(10) = "Jupiter Trojan"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 5.08_bp .AND. &
         elements(:,1) <= 5.33_bp .AND. &
         elements(:,2) >= 0.0_bp .AND. &
         elements(:,2) <= 0.28_bp .AND. &
         elements(:,3) >= 0.0_bp .AND. &
         elements(:,3) <= 40.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(10) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Centaur
    groups(11) = "Centaur"
    mask_array = .FALSE.
    WHERE (periapsis > 5.023_bp .AND. &
         apoapsis < 30.06_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(11) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! TNO Inner belt
    groups(12) = "TNO/Inner belt"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 36.0_bp .AND. &
         elements(:,1) <= 39.0_bp .AND. &
         periapsis > 35.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(12) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Plutino
    groups(13) = "TNO/Plutino"
    mask_array = .FALSE.
    WHERE (elements(:,1) > 39.0_bp .AND. &
         periapsis < 40.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(13) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Classical TNO
    groups(14) = "TNO/Classical"
    mask_array = .FALSE.
    WHERE (elements(:,1) >= 40.0_bp .AND. &
         elements(:,1) <= 48.0_bp .AND. &
         periapsis > 35.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(14) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! TNO Outer belt
    groups(15) = "TNO/Outer belt"
    mask_array = .FALSE.
    WHERE (elements(:,1) > 48.0_bp .AND. &
         periapsis > 36.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(15) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Scattered TNO
    groups(16) = "TNO/Scattered"
    mask_array = .FALSE.
    WHERE (elements(:,1) > 49.0_bp .AND. &
         periapsis < 36.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(16) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! The rest
    groups(17) = "The rest"
    weights(17) = SUM(pdf,mask=mask_array_tot)

    DEALLOCATE(elements, pdf, mask_array, mask_array_tot, &
         periapsis, apoapsis, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END SUBROUTINE getGroupWeights





  !! *Description*:
  !!
  !! Computes an rms value using mean residuals computed using the
  !! sample orbits.
  !!
  !! Returns error.
  !!
  FUNCTION getMeanRMS(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), DIMENSION(6)                :: getMeanRMS
    INTEGER                               :: norb, i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getMeanRMS", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. ASSOCIATED(this%rms_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getMeanRMS", &
            "RMSs are missing -> make an orbit distribution.", 1)
       RETURN       
    END IF

    norb = SIZE(this%rms_arr_cmp,dim=1)
    getMeanRMS = 0.0_bp
    DO i=1,6
       getMeanRMS(i) = SUM(this%rms_arr_cmp(:,i))/norb
    END DO

  END FUNCTION getMeanRMS





  !! *Description*:
  !!
  !! Depending on the orbital solution returns either 
  !! 
  !!  - a discrete sampling of the periapsis distance
  !!    probability-density function, or
  !!
  !!  - the nominal periapsis distance and its 1-sigma uncertainty.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getPeriapsisDistance_SO(this, q)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), DIMENSION(:,:), POINTER     :: q

    TYPE (Orbit) :: orb
    REAL(bp), DIMENSION(6,6) :: partials, covariance1, covariance2
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    REAL(bp) :: jac
    INTEGER :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPeriapsisDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (containsSampledPDF(this)) THEN

       pdf_arr => getPDFValues(this, "keplerian")
       IF (error) THEN
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF

       ALLOCATE(q(SIZE(this%orb_arr_cmp,dim=1),2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "Could not allocate memory (5).", 1)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
       partials = identity_matrix(6)
       DO i=1,SIZE(this%orb_arr_cmp,dim=1)
          CALL getPeriapsisDistance(this%orb_arr_cmp(i), &
               q(i,1), partials=partials(1,1:6))
          IF (error) THEN
             CALL errorMessage("Orbit / getPeriapsisDistance", &
                  "TRACE BACK (5)", 1)
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          jac = ABS(determinant(partials, errstr))
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getPeriapsisDistance", &
                  "Could not compute determinant for Jacobian. " // &
                  TRIM(errstr), 1)
             errstr = ""
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          q(i,2) = pdf_arr(i)*jac
       END DO

       DEALLOCATE(pdf_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPeriapsisDistance", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF

    ELSE

       ALLOCATE(q(1,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "Could not allocate memory (5).", 1)
          RETURN
       END IF
       orb = copy(this%orb_ml_cmp)
       CALL toKeplerian(orb)
       partials = identity_matrix(6)
       CALL getPeriapsisDistance(orb, &
            q(1,1), &
            partials=partials(1,1:6))
       IF (error) THEN
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       CALL NULLIFY(orb)

       covariance1 = getCovarianceMatrix(this, "keplerian")
       covariance2 = MATMUL(MATMUL(partials,covariance1),TRANSPOSE(partials))
       q(1,2) = SQRT(covariance2(1,1))

    END IF

  END SUBROUTINE getPeriapsisDistance_SO





  !! Computes the probability of the target described by the sample
  !! orbits being a PHA (MOID wrt. the Earth <= 0.05 AU) at a given  
  !! moment (determined by the epoch of the sample orbits). 
  !! Computations are done using the 2-body approximation.
  !!
  REAL(bp) FUNCTION getPHAProbability(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this
    TYPE (Orbit)                        :: orbit_earth
    TYPE (Time)                         :: t
    REAL(bp), DIMENSION(:), ALLOCATABLE :: moid, pdf
    REAL(bp), DIMENSION(:,:), POINTER   :: elem
    REAL(bp)                            :: mjd_tdt
    INTEGER                             :: norb, err, i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "No sample orbits available.", 1)
       RETURN
    END IF

    t = getTime(this%orb_arr_cmp(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    mjd_tdt = getMJD(t, "tdt")
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    elem => JPL_ephemeris(mjd_tdt, 3, 11, error)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    CALL NULLIFY(orbit_earth)
    CALL NEW(orbit_earth, elem(1,:), "cartesian", "equatorial", t)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF
    DEALLOCATE(elem, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF
    CALL toKeplerian(orbit_earth)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF

    norb = SIZE(this%orb_arr_cmp,dim=1)
    ALLOCATE(moid(norb), pdf(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    pdf = this%pdf_arr_cmp
    pdf = pdf / SUM(pdf)

    getPHAProbability = 0.0_bp
    DO i=1,norb
       moid(i) = getMOID(orbit_earth, this%orb_arr_cmp(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getPHAProbability", &
               "TRACE BACK (35)", 1)
          RETURN
       END IF
       IF (moid(i) <= 0.05_bp) THEN
          getPHAProbability = getPHAProbability + pdf(i)
       END IF
    END DO

    DEALLOCATE(moid, pdf, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION getPHAProbability





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getNominalOrbit(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Orbit)                       :: getNominalOrbit

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getNominalOrbit", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(this%orb_ml_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getNominalOrbit", &
            "The nominal orbit is not available.", 1)
       RETURN
    END IF

    getNominalOrbit = copy(this%orb_ml_cmp)

  END FUNCTION getNominalOrbit





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  INTEGER FUNCTION getNrOfSampleOrbits(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getNrOfSampleOrbits", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getNrOfSampleOrbits", &
            "Sample orbits do not exist.", 1)       
       RETURN
    END IF

    getNrOfSampleOrbits = SIZE(this%orb_arr_cmp,dim=1)

  END FUNCTION getNrOfSampleOrbits






  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getObservationMasks_SO(this) RESULT(obs_masks)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    LOGICAL, DIMENSION(:,:), POINTER   :: obs_masks
    INTEGER                            :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getObservationMasks", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%obs_masks_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getObservationMasks", &
            "Observation masks have not been allocated.", 1)
       RETURN       
    END IF

    ALLOCATE(obs_masks(SIZE(this%obs_masks_prm,dim=1), &
         SIZE(this%obs_masks_prm,dim=2)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getObservationMasks", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    obs_masks = this%obs_masks_prm

  END FUNCTION getObservationMasks_SO




  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getObservationPairs(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    INTEGER, DIMENSION(:,:), POINTER   :: getObservationPairs

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getObservationPair", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    ALLOCATE(getObservationPairs(SIZE(this%sor_pair_arr_prm,dim=1), &
         SIZE(this%sor_pair_arr_prm,dim=2)))

    getObservationPairs = this%sor_pair_arr_prm

  END FUNCTION getObservationPairs







  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE getParameters_SO(this, &
       dyn_model, &
       perturbers, &
       integration_step, &
       integrator, &
       finite_diff, &
       t_inv, &
       element_type, &
       multiple_objects, &
       outlier_rejection, &
       uniform_pdf, regularized_pdf, &
       accept_multiplier, &
       outlier_multiplier, &
       res_accept, &
       gaussian_pdf, &
       pdf_ml_init_prm, &
       pdf_ml_prm, &
       prob_mass, &
       apriori_a_max, apriori_a_min, apriori_periapsis_max, apriori_periapsis_min, &
       apriori_apoapsis_max, apriori_apoapsis_min, apriori_rho_max, apriori_rho_min, &
       sor_2point_method, &
       sor_norb, sor_norb_sw, &
       sor_ntrial, sor_ntrial_sw, &
       sor_rho1_l, sor_rho1_u, &
       sor_rho2_l, sor_rho2_u, &
       sor_random_obs_selection, & 
       sor_niter, &
       sor_generat_multiplier, &
       sor_deviates, &
       vov_norb, vov_ntrial, vov_norb_iter, vov_ntrial_iter, &
       vov_nmap, vov_niter, vov_scaling, vov_mapping_mask, &
       ls_correction_factor, ls_niter_major_max, ls_niter_major_min, ls_niter_minor, &
       ls_element_mask)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Time), INTENT(out), OPTIONAL :: t_inv
    CHARACTER(len=*), INTENT(out), OPTIONAL :: &
         dyn_model, &
         integrator, &
         element_type, &
         sor_2point_method
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL :: &
         sor_deviates
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL :: &
         res_accept
    REAL(bp), DIMENSION(6,2), INTENT(out), OPTIONAL :: &
         vov_scaling
    REAL(bp), DIMENSION(6), INTENT(out), OPTIONAL :: &
         finite_diff
    REAL(bp), INTENT(out), OPTIONAL :: &
         integration_step, &
         accept_multiplier, &
         outlier_multiplier, &
         pdf_ml_init_prm, &
         pdf_ml_prm, &
         apriori_a_max, &
         apriori_a_min, &
         apriori_periapsis_max, &
         apriori_periapsis_min, &
         apriori_apoapsis_max, &
         apriori_apoapsis_min, &
         apriori_rho_max, &
         apriori_rho_min, &
         prob_mass, &
         sor_rho1_l, &
         sor_rho1_u, &
         sor_rho2_l, &
         sor_rho2_u, &
         sor_generat_multiplier, &
         ls_correction_factor
    INTEGER, INTENT(out), OPTIONAL :: &
         sor_norb, &
         sor_norb_sw, &
         sor_ntrial, &
         sor_ntrial_sw, &
         sor_niter, &
         vov_norb, &
         vov_ntrial, &
         vov_norb_iter, &
         vov_ntrial_iter, &
         vov_nmap, &
         vov_niter, &
         ls_niter_major_max, &
         ls_niter_major_min, &
         ls_niter_minor
    LOGICAL, DIMENSION(:), INTENT(out), OPTIONAL :: &
         perturbers, &
         vov_mapping_mask, &
         ls_element_mask
    LOGICAL, INTENT(out), OPTIONAL :: &
         multiple_objects, &
         regularized_pdf, &
         uniform_pdf, &
         sor_random_obs_selection, &
         gaussian_pdf, &
         outlier_rejection
    INTEGER :: err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getParameters", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(dyn_model)) THEN
       dyn_model = this%dyn_model_prm
    END IF
    IF (PRESENT(perturbers)) THEN
       perturbers = this%perturbers_prm
    END IF
    IF (PRESENT(integration_step)) THEN
       integration_step = this%integration_step_prm
    END IF
    IF (PRESENT(integrator)) THEN
       integrator = this%integrator_prm
    END IF
    IF (PRESENT(finite_diff)) THEN
       IF (ASSOCIATED(this%finite_diff_prm)) THEN
          finite_diff = this%finite_diff_prm
       ELSE
          finite_diff = -1.0_bp
       END IF
    END IF
    IF (PRESENT(t_inv)) THEN
       t_inv = copy(this%t_inv_prm)
    END IF
    IF (PRESENT(element_type)) THEN
       element_type = this%element_type_prm 
    END IF
    IF (PRESENT(multiple_objects)) THEN
       multiple_objects = this%multiple_obj_prm 
    END IF
    IF (PRESENT(outlier_rejection)) THEN
       outlier_rejection = this%outlier_rejection_prm
    END IF
    IF (PRESENT(outlier_multiplier)) THEN
       outlier_multiplier = this%outlier_multiplier_prm 
    END IF
    IF (PRESENT(uniform_pdf)) THEN
       uniform_pdf = this%uniform_pdf_prm 
    END IF
    IF (PRESENT(regularized_pdf)) THEN
       regularized_pdf = this%regularization_prm 
    END IF
    IF (PRESENT(accept_multiplier)) THEN
       accept_multiplier = this%accept_multiplier_prm 
    END IF
    IF (PRESENT(res_accept)) THEN
       IF (ASSOCIATED(this%res_accept_prm)) THEN
          ALLOCATE(res_accept(SIZE(this%res_accept_prm,dim=1), &
               SIZE(this%res_accept_prm,dim=2)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getParameters", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          res_accept = this%res_accept_prm
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getParameters", &
               "Acceptance windows not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(pdf_ml_init_prm)) THEN
       pdf_ml_init_prm = this%pdf_ml_init_prm 
    END IF
    IF (PRESENT(pdf_ml_prm)) THEN
       pdf_ml_prm = this%pdf_ml_prm 
    END IF
    IF (PRESENT(prob_mass)) THEN
       ! prob_mass = getProbabilityMass(this%dchi2_prm)
       ! temporarily
       prob_mass = this%dchi2_prm
    END IF
    IF (PRESENT(apriori_a_max)) THEN
       apriori_a_max = this%apriori_a_max_prm
    END IF
    IF (PRESENT(apriori_a_min)) THEN
       apriori_a_min = this%apriori_a_min_prm
    END IF
    IF (PRESENT(apriori_periapsis_max)) THEN
       apriori_periapsis_max = this%apriori_periapsis_max_prm
    END IF
    IF (PRESENT(apriori_periapsis_min)) THEN
       apriori_periapsis_min = this%apriori_periapsis_min_prm
    END IF
    IF (PRESENT(apriori_apoapsis_max)) THEN
       apriori_apoapsis_max = this%apriori_apoapsis_max_prm
    END IF
    IF (PRESENT(apriori_apoapsis_min)) THEN
       apriori_apoapsis_min = this%apriori_apoapsis_min_prm
    END IF
    IF (PRESENT(apriori_rho_max)) THEN
       apriori_rho_max = this%apriori_rho_max_prm
    END IF
    IF (PRESENT(apriori_rho_min)) THEN
       apriori_rho_min = this%apriori_rho_min_prm
    END IF
    IF (PRESENT(sor_2point_method)) THEN
       sor_2point_method = this%sor_2point_method_prm
    END IF
    IF (PRESENT(sor_norb)) THEN
       sor_norb = this%sor_norb_prm
    END IF
    IF (PRESENT(sor_norb_sw)) THEN
       sor_norb_sw = this%sor_norb_sw_prm 
    END IF
    IF (PRESENT(sor_ntrial)) THEN
       sor_ntrial = this%sor_ntrial_prm 
    END IF
    IF (PRESENT(sor_ntrial_sw)) THEN
       sor_ntrial_sw = this%sor_ntrial_sw_prm
    END IF
    IF (PRESENT(sor_rho1_l)) THEN
       sor_rho1_l = this%sor_rho_prm(1,1) 
    END IF
    IF (PRESENT(sor_rho1_u)) THEN
       sor_rho1_u = this%sor_rho_prm(1,2) 
    END IF
    IF (PRESENT(sor_rho2_l)) THEN
       sor_rho2_l = this%sor_rho_prm(2,1)
    END IF
    IF (PRESENT(sor_rho2_u)) THEN
       sor_rho2_u = this%sor_rho_prm(2,2)
    END IF
    IF (PRESENT(sor_random_obs_selection)) THEN
       sor_random_obs_selection = this%sor_random_obs_prm
    END IF
    IF (PRESENT(sor_niter)) THEN
       sor_niter = this%sor_niter_prm 
    END IF

    IF (PRESENT(gaussian_pdf)) THEN
       gaussian_pdf = this%sor_gaussian_pdf_prm
    END IF
    IF (PRESENT(sor_generat_multiplier)) THEN
       sor_generat_multiplier = this%sor_generat_multiplier_prm 
    END IF
    IF (PRESENT(sor_deviates)) THEN
       IF (ASSOCIATED(this%sor_deviates_prm)) THEN
          ALLOCATE(sor_deviates(SIZE(this%sor_deviates_prm,dim=1), &
               SIZE(this%sor_deviates_prm,dim=2), &
               SIZE(this%sor_deviates_prm,dim=3)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getParameters", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          sor_deviates = this%sor_deviates_prm
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getParameters", &
               "Generation windows not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(vov_norb)) THEN
       vov_norb = this%vov_norb_prm
    END IF
    IF (PRESENT(vov_ntrial)) THEN
       vov_ntrial = this%vov_ntrial_prm 
    END IF
    IF (PRESENT(vov_norb_iter)) THEN
       vov_norb_iter = this%vov_norb_iter_prm 
    END IF
    IF (PRESENT(vov_ntrial_iter)) THEN
       vov_ntrial_iter = this%vov_ntrial_iter_prm 
    END IF
    IF (PRESENT(vov_nmap)) THEN
       vov_nmap = this%vov_nmap_prm 
    END IF
    IF (PRESENT(vov_niter)) THEN
       vov_niter = this%vov_niter_prm 
    END IF
    IF (PRESENT(vov_scaling)) THEN
       vov_scaling = this%vov_scaling_prm
    END IF
    IF (PRESENT(vov_mapping_mask)) THEN
       vov_mapping_mask = this%vov_mapping_mask_prm 
    END IF
    IF (PRESENT(ls_correction_factor)) THEN
       ls_correction_factor = this%ls_corr_fac_prm 
    END IF
    IF (PRESENT(ls_niter_major_max)) THEN
       ls_niter_major_max = this%ls_niter_major_max_prm
    END IF
    IF (PRESENT(ls_niter_major_min)) THEN
       ls_niter_major_min = this%ls_niter_major_min_prm
    END IF
    IF (PRESENT(ls_niter_minor)) THEN
       ls_niter_minor = this%ls_niter_minor_prm
    END IF
    IF (PRESENT(ls_element_mask)) THEN
       ls_element_mask = this%ls_elem_mask_prm 
    END IF

  END SUBROUTINE getParameters_SO





  !! *Description*:
  !!
  !! Returns error.
  !!
  FUNCTION getPDFValues(this, element_type)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER    :: getPDFValues
    CHARACTER(len=*), INTENT(in), OPTIONAL :: element_type
    REAL(bp), DIMENSION(:), ALLOCATABLE :: jac
    INTEGER                            :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPDFValues", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPDFValues", &
            "PDF values do not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getPDFValues(SIZE(this%pdf_arr_cmp)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPDFValues", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    IF (PRESENT(element_type)) THEN
       ALLOCATE(jac(SIZE(this%pdf_arr_cmp)))
       jac = 1.0_bp
       IF (this%element_type_prm == "cartesian" .AND. &
            element_type == "keplerian") THEN
          ! Changed 2008-12-14
          !jac = this%jac_arr_cmp(:,2)
          jac = 1.0_bp/this%jac_arr_cmp(:,2)
       ELSE IF (this%element_type_prm == "keplerian" .AND. &
            element_type == "cartesian") THEN
          ! Changed 2008-12-14
          !jac = 1.0_bp/this%jac_arr_cmp(:,2)
          jac = this%jac_arr_cmp(:,2)
       ELSE IF (this%element_type_prm == "cometary" .OR. &
            element_type == "cometary") THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPDFValues", &
               "Not yet implemented for cometary elements.", 1)
          RETURN
       END IF
       getPDFValues = this%pdf_arr_cmp*jac
       DEALLOCATE(jac, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPDFValues", &
               "Could not deallocate memory.", 1)
          RETURN
       END IF
    ELSE
       getPDFValues = this%pdf_arr_cmp
    END IF

  END FUNCTION getPDFValues





  !! *Description*:
  !!
  !! Returns a discrete sampling of phase-angle probability-density
  !! function for an object on this orbit at a given epoch. The
  !! observer coordinates are assumed to be exact.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getPhaseAngle_SO_pdf(this, observer, phase_angle_pdf)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)         :: this
    TYPE (CartesianCoordinates), INTENT(in)       :: observer
    REAL(bp), DIMENSION(:,:), POINTER             :: phase_angle_pdf

    REAL(bp), DIMENSION(6,6) :: partials
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    REAL(bp) :: jac
    INTEGER :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Observer object has not been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%orb_arr_cmp) .OR. &
         .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Orbital-element pdf missing.", 1)
       RETURN
    END IF

    pdf_arr => getPDFValues(this, "cartesian")
    IF (error) THEN
       CALL errorMessage("Orbit / getPhaseAngle", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(pdf_arr, stat=err)
       RETURN
    END IF

    ALLOCATE(phase_angle_pdf(SIZE(this%orb_arr_cmp,dim=1),2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(pdf_arr, stat=err)
       RETURN
    END IF
    partials = identity_matrix(6)
    DO i=1,SIZE(this%orb_arr_cmp,dim=1)
       CALL getPhaseAngle(this%orb_arr_cmp(i), observer, &
            phase_angle_pdf(i,1), partials(1,1:6))
       IF (error) THEN
          CALL errorMessage("Orbit / getPhaseAngle", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
       jac = ABS(determinant(partials, errstr))
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPhaseAngle", &
               "Could not compute determinant for Jacobian.", 1)
          errstr = ""
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
       phase_angle_pdf(i,2) = pdf_arr(i)*jac
    END DO

    DEALLOCATE(pdf_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPhaseAngle", &
            "Could not deallocate memory (10).", 1)
       RETURN
    END IF

  END SUBROUTINE getPhaseAngle_SO_pdf





  !! *Description*:
  !!
  !! Returns a phase angles corresponding to one or more orbits as
  !! seen by one or more observers.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getPhaseAngles_SO(this, observers, phase_angles)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)                 :: this
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: observers
    REAL(bp), DIMENSION(:,:), POINTER                     :: phase_angles

    REAL(bp), DIMENSION(:), POINTER :: phase_angles_
    INTEGER :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPhaseAngles", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1,SIZE(observers)
       IF (.NOT.exist(observers(i))) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPhaseAngles", &
               "All observer objects have not been initialized.", 1)
          RETURN
       END IF
    END DO

    IF (containsSampledPDF(this)) THEN
       ALLOCATE(phase_angles(SIZE(this%orb_arr_cmp),SIZE(observers)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPhaseAngles", &
               "Could not allocate memory (5).", 1)
          RETURN
       END IF
       DO i=1,SIZE(this%orb_arr_cmp)
          CALL getPhaseAngles(this%orb_arr_cmp(i), observers, phase_angles_)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / getPhaseAngles", &
                  "TRACE BACK (5).", 1)
             RETURN
          END IF
          phase_angles(i,:) = phase_angles_
          DEALLOCATE(phase_angles_, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getPhaseAngles", &
                  "Could not deallocate memory (5).", 1)
             RETURN
          END IF
       END DO
    ELSE
       CALL getPhaseAngles(this%orb_ml_cmp, observers, phase_angles_)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getPhaseAngles", &
               "TRACE BACK (10).", 1)
          RETURN
       END IF
       ALLOCATE(phase_angles(1,SIZE(phase_angles_)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPhaseAngles", &
               "Could not allocate memory (10).", 1)
          RETURN
       END IF
       phase_angles(1,:) = phase_angles_
       DEALLOCATE(phase_angles_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getPhaseAngles", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF
    END IF

  END SUBROUTINE getPhaseAngles_SO





  !! *Description*:
  !!
  !! Returns the phase angle and its uncertainty for an object on
  !! this orbit at a given epoch. The observer coordinates are assumed
  !! to be exact.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getPhaseAngle_SO_point(this, observer, phase_angle, sigma)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)   :: this
    TYPE (CartesianCoordinates), INTENT(in) :: observer
    REAL(bp), INTENT(out) :: phase_angle
    REAL(bp), INTENT(out), OPTIONAL :: sigma

    REAL(bp), DIMENSION(6,6) :: cov
    REAL(bp), DIMENSION(1,6) :: partials
    REAL(bp), DIMENSION(1,1) :: variance

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Observer object has not been initialized.", 1)
       RETURN
    END IF

    ! Phase angle
    IF (.NOT.ASSOCIATED(this%orb_arr_cmp) .AND. &
         exist(this%orb_ml_cmp)) THEN
       CALL getPhaseAngle(this%orb_ml_cmp, observer, phase_angle, &
            partials(1,:))
       IF (error) THEN
          CALL errorMessage("Orbit / getPhaseAngle", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       IF (PRESENT(sigma)) THEN
          ! Compute phase-angle uncertainty (sigma):
          cov = getCovarianceMatrix(this, &
               getElementType(this%orb_ml_cmp), &
               getFrame(this%orb_ml_cmp))
          IF (error) THEN
             CALL errorMessage("Orbit / getPhaseAngle", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
          variance(1:1,1:1) = MATMUL(MATMUL(partials(1:1,1:3), &
               cov(1:3,1:3)), TRANSPOSE(partials(1:1,1:3)))
          sigma = SQRT(variance(1,1))
       END IF
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Either the orbital-element p.d.f. exist, " // &
            "or the maximum likelihood orbit does not exist.", 1)
       RETURN
    END IF

  END SUBROUTINE getPhaseAngle_SO_point





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getPositionDistribution(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(:,:), POINTER  :: getPositionDistribution
    INTEGER                            :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPositionDistribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPositionDistribution", &
            "Sample orbits do not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getPositionDistribution(SIZE(this%orb_arr_cmp),3), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPositionDistribution", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, SIZE(this%orb_arr_cmp)
       getPositionDistribution(i,:) = getPosition(this%orb_arr_cmp(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getPositionDistribution", &
               "TRACE BACK", 1)
          RETURN
       END IF
    END DO

  END FUNCTION getPositionDistribution







  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getRangeBounds_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(4)             :: getRangeBounds_SO

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRangeBounds", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    getRangeBounds_SO(1) = this%sor_rho_prm(1,1) 
    getRangeBounds_SO(2) = this%sor_rho_prm(1,2)
    getRangeBounds_SO(3) = this%sor_rho_prm(2,1)
    getRangeBounds_SO(4) = this%sor_rho_prm(2,2)

  END FUNCTION getRangeBounds_SO





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getReducedChi2Distribution(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER    :: getReducedChi2Distribution
    INTEGER                            :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getReducedChi2Distribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%rchi2_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getReducedChi2Distribution", &
            "Reduced chi2 array does not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getReducedChi2Distribution(SIZE(this%rchi2_arr_cmp,dim=1)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getReducedChi2Distribution", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    getReducedChi2Distribution(:) = this%rchi2_arr_cmp

  END FUNCTION getReducedChi2Distribution





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getResidualDistribution(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this
    REAL(bp), DIMENSION(:,:,:), POINTER :: getResidualDistribution
    INTEGER                             :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResidualDistribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%res_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResidualDistribution", &
            "Residuals do not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getResidualDistribution(SIZE(this%res_arr_cmp,dim=1), &
         SIZE(this%res_arr_cmp,dim=2), SIZE(this%res_arr_cmp,dim=3)), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResidualDistribution", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    getResidualDistribution(:,:,:) = this%res_arr_cmp

  END FUNCTION getResidualDistribution





  !! *Description*:
  !!
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getResiduals_single(this, orb)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)                 :: this
    TYPE (Orbit), INTENT(in)                           :: orb
    REAL(bp), DIMENSION(:,:), POINTER                  :: getResiduals_single
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: computed_scoords, &
         observed_scoords
    REAL(bp), DIMENSION(:,:), ALLOCATABLE              :: observed_coords, &
         computed_coords
    INTEGER                                            :: err, i, nobs

    nobs = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    ALLOCATE(getResiduals_single(nobs,6), observed_coords(nobs,6), &
         computed_coords(nobs,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    observed_scoords => getObservationSCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

    obsy_ccoords => getObservatoryCCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF

    CALL getEphemerides(orb, obsy_ccoords, computed_scoords)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF

    observed_coords = 0.0_bp
    computed_coords = 0.0_bp
    DO i=1,nobs
       observed_coords(i,:) = getCoordinates(observed_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getResiduals", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
       computed_coords(i,:) = getCoordinates(computed_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getResiduals", &
               "TRACE BACK (30)", 1)
          RETURN
       END IF
    END DO

    getResiduals_single(1:nobs,1:6) = observed_coords(1:nobs,1:6) - &
         computed_coords(1:nobs,1:6)        
    getResiduals_single(1:nobs,2) = getResiduals_single(1:nobs,2) * &
         COS(observed_coords(1:nobs,3))
    DO i=1,nobs
       IF (ABS(getResiduals_single(i,2)) > pi) THEN
          IF (observed_coords(i,2) < computed_coords(i,2)) THEN
             observed_coords(i,2) = observed_coords(i,2) + two_pi
          ELSE
             computed_coords(i,2) = computed_coords(i,2) + two_pi
          END IF
          getResiduals_single(i,2) = (observed_coords(i,2) - &
               computed_coords(i,2)) * COS(observed_coords(i,3))
       END IF
    END DO
    DO i=1,SIZE(observed_scoords)
       CALL NULLIFY(observed_scoords(i))
    END DO
    DO i=1,SIZE(obsy_ccoords)
       CALL NULLIFY(obsy_ccoords(i))
    END DO
    DO i=1,SIZE(computed_scoords)
       CALL NULLIFY(computed_scoords(i))
    END DO
    DEALLOCATE(observed_coords, computed_coords, observed_scoords, &
         obsy_ccoords, computed_scoords, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION getResiduals_single





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE getResults_SO(this,&
       reg_apr_arr, jac_arr, &
       sor_norb_cmp, sor_ntrial_cmp,& 
       sor_rho_cmp, sor_niter_cmp, &
       sor_rho_arr_cmp, sor_rho_histo_cmp, &
       vov_norb_cmp, vov_ntrial_cmp, &
       vov_niter_cmp, vov_scaling_cmp, &
       vov_map_cmp, vov_scaling_ready_cmp)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)              :: this
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL     :: &
         jac_arr
    REAL(bp), DIMENSION(:), POINTER, OPTIONAL       :: reg_apr_arr
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL     :: sor_rho_arr_cmp
    REAL(bp), DIMENSION(2,2), INTENT(out), OPTIONAL :: &
         sor_rho_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: sor_norb_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: sor_ntrial_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: sor_niter_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: sor_rho_histo_cmp
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL     :: vov_map_cmp    
    REAL(bp), DIMENSION(6,2), INTENT(out), OPTIONAL :: vov_scaling_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: vov_norb_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: vov_ntrial_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: vov_niter_cmp
    LOGICAL, DIMENSION(6,2), INTENT(out), OPTIONAL  :: vov_scaling_ready_cmp
    INTEGER :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResults", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(reg_apr_arr)) THEN
       IF (ASSOCIATED(this%reg_apr_arr_cmp)) THEN
          ALLOCATE(reg_apr_arr(SIZE(this%reg_apr_arr_cmp)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getResults", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          reg_apr_arr = this%reg_apr_arr_cmp
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getResults", &
               "Apriori array not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(jac_arr)) THEN
       IF (ASSOCIATED(this%jac_arr_cmp)) THEN
          ALLOCATE(jac_arr(SIZE(this%jac_arr_cmp,dim=1), &
               SIZE(this%jac_arr_cmp,dim=2)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getResults", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          jac_arr = this%jac_arr_cmp
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getResults", &
               "Jacobians not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(sor_norb_cmp)) THEN
       sor_norb_cmp = this%sor_norb_cmp
    END IF
    IF (PRESENT(sor_ntrial_cmp)) THEN
       sor_ntrial_cmp = this%sor_ntrial_cmp
    END IF

    IF (PRESENT(sor_rho_cmp)) THEN
       sor_rho_cmp = this%sor_rho_cmp
    END IF
    IF (PRESENT(sor_niter_cmp)) THEN
       sor_niter_cmp = this%sor_niter_cmp
    END IF
    IF (PRESENT(sor_rho_arr_cmp)) THEN
       IF (ASSOCIATED(this%sor_rho_arr_cmp)) THEN
          ALLOCATE(sor_rho_arr_cmp(this%sor_norb_cmp,2), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getResults", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          sor_rho_arr_cmp = this%sor_rho_arr_cmp
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getResults", &
               "Generated rho values not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(sor_rho_histo_cmp)) THEN
       sor_rho_histo_cmp = this%sor_rho_histo_cmp
    END IF
    IF (PRESENT(vov_norb_cmp)) THEN
       vov_norb_cmp = this%vov_norb_cmp
    END IF
    IF (PRESENT(vov_ntrial_cmp)) THEN
       vov_ntrial_cmp = this%vov_ntrial_cmp
    END IF
    IF (PRESENT(vov_niter_cmp)) THEN
       vov_niter_cmp = this%vov_niter_cmp
    END IF
    IF (PRESENT(vov_scaling_cmp)) THEN
       vov_scaling_cmp = this%vov_scaling_cmp
    END IF
    IF (PRESENT(vov_map_cmp)) THEN
       IF (ASSOCIATED(this%vov_map_cmp)) THEN
          ALLOCATE(vov_map_cmp(this%vov_nmap_prm,12), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getResults", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          vov_map_cmp = this%vov_map_cmp
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getResults", &
               "VoV map not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(vov_scaling_ready_cmp)) THEN
       vov_scaling_ready_cmp = this%vov_scaling_ready_cmp 
    END IF

  END SUBROUTINE getResults_SO





  !! *Description*:
  !!
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getRhoDistribution(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(:,:), POINTER  :: getRhoDistribution
    INTEGER                            :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRhoDistribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%sor_rho_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRhoDistribution", &
            "Rho distribution does not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getRhoDistribution(SIZE(this%sor_rho_arr_cmp,dim=1), &
         SIZE(this%sor_rho_arr_cmp,dim=2)), stat=err)
    getRhoDistribution = this%sor_rho_arr_cmp

  END FUNCTION getRhoDistribution





  !! *Description*:
  !!
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getRMS_single(this, orb)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Orbit), INTENT(inout)        :: orb
    REAL(bp), DIMENSION(6)             :: getRMS_single
    REAL(bp), DIMENSION(:,:), POINTER  :: residuals
    INTEGER                            :: i, nobs, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRMS", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    residuals => getResiduals(this, orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getRMS", &
            "TRACE BACK", 1)
       RETURN       
    END IF

    nobs = SIZE(residuals,dim=1)
    getRMS_single = 0.0_bp
    DO i=1,6
       getRMS_single(i) = SQRT(SUM(residuals(1:nobs,i)**2.0_bp)/nobs)
    END DO

    DEALLOCATE(residuals, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRMS", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF


  END FUNCTION getRMS_single





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getRMSDistribution(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this
    REAL(bp), DIMENSION(:,:), POINTER   :: getRMSDistribution
    INTEGER                             :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRMSDistribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%rms_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRMSDistribution", &
            "rms distribution does not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getRMSDistribution(SIZE(this%rms_arr_cmp,dim=1), &
         SIZE(this%rms_arr_cmp,dim=2)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRMSDistribution", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    getRMSDistribution(:,:) = this%rms_arr_cmp

  END FUNCTION getRMSDistribution






  !! *Description*:
  !!
  !! Returns the total O-C residual rms in the Pythagorean sense
  !! (rms_tot = sqrt(rms_RA^2 + rms_Dec^2)) for each sample orbit.
  !!
  !! Returns error.
  !!
  FUNCTION getRMSValues(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER  :: getRMSValues
    INTEGER                            :: i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getRMSValues", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ALLOCATE(getRMSValues(SIZE(this%rms_arr_cmp,dim=1)))
    DO i=1,SIZE(this%rms_arr_cmp,dim=1)
       getRMSValues(i) = SQRT(SUM(this%rms_arr_cmp(i,2:3)**2.0_bp))
    END DO

  END FUNCTION getRMSValues





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getSampleOrbit_SO(this, index)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    INTEGER, INTENT(in)                :: index
    TYPE (Orbit)                       :: getSampleOrbit_SO

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getSampleOrbit", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getSampleOrbit", &
            "Sample orbits do not exist.", 1)
       RETURN
    END IF

    getSampleOrbit_SO = copy(this%orb_arr_cmp(index))

  END FUNCTION getSampleOrbit_SO





  !! *Description*:
  !!
  !! Returns orbits corresponding to the orbital-element PDF.
  !!
  !! Returns error.
  !!
  FUNCTION getSampleOrbits(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this
    TYPE (Orbit), DIMENSION(:), POINTER :: getSampleOrbits
    INTEGER                             :: err, i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getSampleOrbits", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getSampleOrbits", &
            "Sample orbits do not exist.", 1)
       RETURN
    END IF

    ALLOCATE(getSampleOrbits(SIZE(this%orb_arr_cmp)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getSampleOrbits", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    DO i=1, SIZE(this%orb_arr_cmp)
       getSampleOrbits(i) = copy(this%orb_arr_cmp(i))
    END DO

  END FUNCTION getSampleOrbits





  !! *Description*:
  !!
  !! Returns error.
  !!
  FUNCTION getStandardDeviations_SO(this, element_type)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(6)             :: getStandardDeviations_SO, stdev
    CHARACTER(len=*), INTENT(in), OPTIONAL :: element_type
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: element_arr
    REAL(bp), DIMENSION(:), POINTER :: pdf
    CHARACTER(len=ELEMENT_TYPE_LEN)    :: element_type_
    INTEGER                            :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPDFValues", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getPDFValues", &
            "PDF values do not exist.", 1)
       RETURN
    END IF

    IF (PRESENT(element_type)) THEN
       element_type_ = element_type
    ELSE
       element_type_ = this%element_type_prm
    END IF

    ALLOCATE(element_arr(SIZE(this%orb_arr_cmp,dim=1),6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getStandardDeviations", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    IF (PRESENT(element_type)) THEN
       pdf => getPDFValues(this, element_type)
    ELSE
       pdf = this%pdf_arr_cmp
    ENDIF

    DO i=1,SIZE(this%orb_arr_cmp,dim=1)
       !       IF (element_type_ == "cartesian") THEN
       !          CALL rotateToEcliptic(this%orb_arr_cmp(l))
       !       END IF
       element_arr(i,1:6) = getElements(this%orb_arr_cmp(i), element_type_)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getStandardDeviations", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
    END DO
    DO i=1,6
       CALL moments(element_arr(:,i), pdf=pdf, std_dev=stdev(i), &
            error=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getStandardDeviations", &
               "TRACE BACK (15) " // TRIM(errstr), 1)
          RETURN
       END IF
    END DO
    getStandardDeviations_SO = stdev

    DEALLOCATE(pdf, element_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getStandardDeviations", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION getStandardDeviations_SO





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getTime_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Time)                        :: getTime_SO
    TYPE (Time)                        :: t
    INTEGER                            :: i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (exist(this%orb_ml_cmp) .AND. .NOT. ASSOCIATED(this%orb_arr_cmp)) THEN
       t = getTime(this%orb_ml_cmp)
    ELSE IF (.NOT.exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%orb_arr_cmp)) THEN
       t = getTime(this%orb_arr_cmp(1))
       DO i=2,SIZE(this%orb_arr_cmp,dim=1)
          IF (.NOT.equal(getTime(this%orb_arr_cmp(i)),t)) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getTime", &
                  "Sample orbits have different epochs.", 1)
             CALL NULLIFY(t)
             RETURN
          END IF
       END DO
    ELSE IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%orb_arr_cmp)) THEN
       t = getTime(this%orb_ml_cmp)
       DO i=1,SIZE(this%orb_arr_cmp,dim=1)
          IF (.NOT.equal(getTime(this%orb_arr_cmp(i)),t)) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getTime", &
                  "Sample orbits and the nominal orbit have different epochs.", 1)
             CALL NULLIFY(t)
             RETURN
          END IF
       END DO
    END IF

    getTime_SO = copy(t)

  END FUNCTION getTime_SO





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  INTEGER FUNCTION getTrials(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getTrials", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getTrials = this%sor_ntrial_cmp

  END FUNCTION getTrials






  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE includeObservations(this, obss)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Observations), INTENT(in)       :: obss

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / includeObservations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / includeObservations", &
            "Observations not intialized.", 1)
       RETURN
    END IF

    IF (exist(this%obss)) THEN
       CALL NULLIFY(this%obss)
       DEALLOCATE(this%obs_masks_prm)
    END IF

    this%obss = copy(obss)
    this%obs_masks_prm => getObservationMasks(this%obss)

    ! Initialize other variables needed in the computation
    CALL setObservationPair(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / includeObservations", &
            "TRACE BACK", 1)
       RETURN
    END IF

  END SUBROUTINE includeObservations






  !! *Description*:
  !!
  !! Computation of the least-squares orbital elements (Cartesian or
  !! Keplerian) and their covariances in the two-body and full
  !! many-body approaches.
  !!
  !! Returns error.
  !!
  SUBROUTINE leastSquares_SO(this, preliminary_orbit)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), INTENT(in)              :: preliminary_orbit    

    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: obs_scoords, ephemerides
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords
    TYPE (Time) :: t, t_
    TYPE (Orbit) :: orb
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: obsy_codes
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type, element_type_ls
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=DYN_MODEL_LEN) :: dyn_model, dyn_model_
    CHARACTER(len=INTEGRATOR_LEN) :: orb_integrator
    CHARACTER(len=32) :: str
    REAL(bp), DIMENSION(:,:,:), POINTER :: partials_arr, &
         cov_mat_obs, inform_mat_obs_bd
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE :: design_mat ! (data,obs,parameter)
    REAL(bp), DIMENSION(:,:), POINTER :: stdev_arr_obs, &
         elements_iter_arr, rmss_iter_arr
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: observed, & ! incl. cos(dec)
         computed, & ! incl. cos(dec)
         residuals
    REAL(bp), DIMENSION(6,6) :: cov
    REAL(bp), DIMENSION(6) :: elements, elements_pre, &
         elem_diff, stdev, stdev_max, rms_pre, orb_finite_diff
    REAL(bp) :: sigma_multiplier_rms, &
         correction_factor, orb_integration_step, obs_, comp_
    INTEGER, DIMENSION(6) :: nobs_arr
    INTEGER, DIMENSION(3) :: approach_direction
    INTEGER :: i, j, k, err, iiter, approach, &
         final_iterations, nobs, err_verb_
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: obs_masks_
    LOGICAL, DIMENSION(6) :: element_mask
    LOGICAL :: elements_converge, rmss_converge, elements_exist

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "LEAST SQUARES"
    END IF

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(preliminary_orbit)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "Preliminary orbit has not been initialized.", 1)
       RETURN
    END IF

    IF (this%accept_multiplier_prm < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "Acceptance criterion has not been defined.", 1)
       RETURN
    END IF

    IF (this%outlier_rejection_prm) THEN
       IF (this%outlier_multiplier_prm < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / leastSquares", &
               "Outlier criterion has not been defined.", 1)
          RETURN
       END IF
    END IF
    sigma_multiplier_rms = this%accept_multiplier_prm
    stdev_max = 0.0_bp
    ! Number of previous iteration rounds to consider
    iiter = 3

    IF (.NOT.ASSOCIATED(this%obs_masks_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "Observation mask is missing.", 1)
       RETURN
    END IF

    IF (this%ls_niter_minor_prm <= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "ls_niter_minor_prm undefined.", 1)
       RETURN       
    END IF

    IF (this%ls_niter_major_min_prm <= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "ls_niter_major_min_prm undefined.", 1)
       RETURN       
    END IF

    IF (this%ls_niter_major_max_prm <= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "ls_niter_major_max_prm undefined.", 1)
       RETURN       
    END IF

    ! Allocate memory:
    nobs = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    ALLOCATE(observed(nobs,6), computed(nobs,6), residuals(nobs,6), &
         design_mat(6,nobs,6), elements_iter_arr(iiter,8), &
         rmss_iter_arr(iiter,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    orb = copy(preliminary_orbit)

    ! Inversion epoch is equal to epoch of initial orbit:
    t = getTime(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF

    ! Find parameters of the peliminary orbit:
    CALL getParameters(orb, dyn_model=dyn_model)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    IF (dyn_model /= this%dyn_model_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "Inconsistent propagation schemes: " // &
            "orb=" // TRIM(dyn_model) // " and storb=" // &
            TRIM(this%dyn_model_prm) // ".", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    dyn_model_ = dyn_model

    ! Initialize elements:
    frame = getFrame(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (15)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    element_type = this%element_type_prm
    CALL locase(element_type, error)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "The element type string contains forbidden characters.", 1)
       RETURN
    END IF
    SELECT CASE (TRIM(element_type))
    CASE ("cartesian")
       CALL toCartesian(orb, frame=frame)
    CASE ("keplerian")
       CALL toKeplerian(orb)
    CASE default
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / leastSquares", &
            "Can not use elements of type: " // TRIM(element_type), 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END SELECT
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (25)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    elements_pre(1:6) = getElements(orb, element_type)
    DO i=1,iiter
       elements_iter_arr(i,1:6) = elements_pre(1:6)
    END DO
    element_type_ls = element_type
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (27)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(1X,3(1X,A,1X,A))") "Using", TRIM(element_type), &
            "elements."
    END IF

    ! Observations and observers:
    obs_scoords => getObservationSCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (30)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    ALLOCATE(obs_masks_(SIZE(this%obs_masks_prm,dim=1),SIZE(this%obs_masks_prm,dim=2)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "Could not allocate memory.", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    obs_masks_ = this%obs_masks_prm
    DO i=1,6
       nobs_arr(i) = COUNT(obs_masks_(:,i))
    END DO
    obsy_codes => getObservatoryCodes(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (45)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF
    obsy_ccoords => getObservatoryCCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (45)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF
    DO i=1,nobs
       observed(i,:) = getCoordinates(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "leastSquares", &
               "TRACE BACK (50)", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(ephemerides, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(design_mat, stat=err) 
          DEALLOCATE(stdev_arr_obs, stat=err) 
          DEALLOCATE(elements_iter_arr, stat=err)
          DEALLOCATE(rmss_iter_arr, stat=err) 
          DEALLOCATE(observed, stat=err)
          DEALLOCATE(computed, stat=err) 
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(obs_masks_, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          RETURN
       END IF
    END DO
    observed(:,2) = observed(:,2)*COS(observed(:,3))
    DEALLOCATE(obs_scoords, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "Could not deallocate memory (5).", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF

    ! Covariance matrix for observations:
    stdev_arr_obs => getStandardDeviations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (51)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF
    stdev_max(2) = MAXVAL(stdev_arr_obs(:,2))
    stdev_max(3) = MAXVAL(stdev_arr_obs(:,3))

    inform_mat_obs_bd => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (53)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(inform_mat_obs_bd, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF

    ! Compute positions and partial derivatives:
    CALL getEphemerides(orb, obsy_ccoords, ephemerides, &
         partials_arr=partials_arr)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (55)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(inform_mat_obs_bd, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF
    CALL propagate(orb, t)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "TRACE BACK (56)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(inform_mat_obs_bd, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF

    ! Write design matrix and add cosine term to it and the computed positions:
    design_mat = 0.0_bp
    DO j=1,nobs
       computed(j,:) = getCoordinates(ephemerides(j))
       computed(j,2) = computed(j,2)*COS(observed(j,3))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "leastSquares", &
               "TRACE BACK (60)", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(ephemerides, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(design_mat, stat=err) 
          DEALLOCATE(stdev_arr_obs, stat=err) 
          DEALLOCATE(elements_iter_arr, stat=err)
          DEALLOCATE(rmss_iter_arr, stat=err) 
          DEALLOCATE(inform_mat_obs_bd, stat=err) 
          DEALLOCATE(observed, stat=err)
          DEALLOCATE(computed, stat=err) 
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(obs_masks_, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          RETURN
       END IF
       design_mat(2,j,:) = partials_arr(2,:,j)*COS(observed(j,3))
       design_mat(3,j,:) = partials_arr(3,:,j)
    END DO
    DEALLOCATE(ephemerides, partials_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "Could not deallocate memory (10).", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(ephemerides, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(cov_mat_obs, stat=err)
       DEALLOCATE(design_mat, stat=err) 
       DEALLOCATE(stdev_arr_obs, stat=err) 
       DEALLOCATE(elements_iter_arr, stat=err)
       DEALLOCATE(rmss_iter_arr, stat=err) 
       DEALLOCATE(inform_mat_obs_bd, stat=err) 
       DEALLOCATE(observed, stat=err)
       DEALLOCATE(computed, stat=err) 
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(obs_masks_, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       RETURN
    END IF


    ! Output orbital elements and sky-plane residuals and rms if needed:
    IF (info_verb >= 2) THEN
       t_ = getTime(orb)
       WRITE(stdout,"(2X,A,I0,A)") "Initial elements, residuals, and rms:"
       str = getCalendarDateString(t_, "TT")
       SELECT CASE (element_type)
       CASE ("keplerian")
          elements = getElements(orb, element_type)
          WRITE(stdout,"(2X,A,6(1X,F20.13),1X,A)") "Kep: ", &
               elements(1:2), elements(3:6)/rad_deg, &
               TRIM(str)
       CASE ("cartesian")
          elements = getElements(orb, element_type, frame=frame)
          WRITE(stdout,'(2X,A,6(1X,F20.13),1X,A)') "Car: ", &
               elements(1:6),TRIM(str)
       END SELECT
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "leastSquares", &
               "TRACE BACK (65)", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(ephemerides, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(design_mat, stat=err) 
          DEALLOCATE(stdev_arr_obs, stat=err) 
          DEALLOCATE(elements_iter_arr, stat=err)
          DEALLOCATE(rmss_iter_arr, stat=err) 
          DEALLOCATE(inform_mat_obs_bd, stat=err) 
          DEALLOCATE(observed, stat=err)
          DEALLOCATE(computed, stat=err) 
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(obs_masks_, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          RETURN
       END IF
       CALL NULLIFY(t_)
       WRITE(stdout,"(2X,A)") "Residuals RA & Dec [as]:"
    END IF
    DO j=1,nobs
       residuals(j,1:6) = observed(j,1:6) - computed(j,1:6)
       IF (ABS(residuals(j,2)) > pi) THEN
          obs_ = observed(j,2)
          comp_ = computed(j,2)
          IF (obs_ < comp_) THEN
             obs_ = obs_ + two_pi
          ELSE
             comp_ = comp_ + two_pi
          END IF
          residuals(j,2) = obs_ - comp_
       END IF
       IF (info_verb >= 2) THEN
          IF (ALL(obs_masks_(j,2:3))) THEN
             WRITE(stdout,"(2X,A,2(F20.7,1X),A,1X,A)") " ", &
                  residuals(j,2:3)/rad_asec, " ", TRIM(obsy_codes(j))
          ELSE
             WRITE(stdout,"(2X,A,2(F15.7,1X),A,1X,A)") "(", &
                  residuals(j,2:3)/rad_asec, ")", TRIM(obsy_codes(j))
          END IF
       END IF
    END DO
    rmss_iter_arr = 0.0_bp
    rms_pre = 0.0_bp
    rms_pre(2) = SQRT(SUM(residuals(:,2)**2.0_bp,mask=obs_masks_(:,2))/nobs_arr(2))
    rms_pre(3) = SQRT(SUM(residuals(:,3)**2.0_bp,mask=obs_masks_(:,3))/nobs_arr(3))
    DO i=1,iiter
       rmss_iter_arr(i,:) = rms_pre
    END DO
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,2(1X,F15.7))") "RA & Dec RMS [as]: ", &
            rmss_iter_arr(iiter,2)/rad_asec, rmss_iter_arr(iiter,3)/rad_asec
       WRITE(stdout,"(2X,A)") ""
    END IF

    final_iterations = 0
    approach = 11
    approach_direction(1:3) = 1
    correction_factor = this%ls_corr_fac_prm
    element_mask(1:6) = this%ls_elem_mask_prm

    ! Main LSL loop
    ls_i: DO i=1,this%ls_niter_major_max_prm

       ! Perform linear least-squares fit:
       elements_converge = .FALSE.
!!$       rmss_converge = .FALSE.
       rmss_converge = .TRUE.
       ! "Correction factor" loop
       ls_j: DO j=1,this%ls_niter_minor_prm

          ! Save results from last two rounds for convergence checks
          elements_iter_arr(1:iiter-1,1:6) = elements_iter_arr(2:iiter,1:6)
          rmss_iter_arr(1:iiter-1,1:6) = rmss_iter_arr(2:iiter,1:6)
          rmss_iter_arr(iiter,1:6) = 0.0_bp

          ! Perform least squares fit:
          CALL leastSquares(observed(1:nobs,1:6), &
               inform_mat_obs_bd(1:nobs,1:6,1:6), &
               obs_masks_(1:nobs,1:6), &
               computed(1:nobs,1:6), &
               design_mat(1:6,1:nobs,1:6), &
               correction_factor, &
               element_mask(1:6), &
               cov(1:6,1:6), &
               elements_iter_arr(iiter,1:6), &
               errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "Could not find a linear least-squares solution. " // &
                  TRIM(errstr), 1)
             errstr = ""
             elements_converge = .FALSE.
             elements_iter_arr(iiter,1:6) = elements_iter_arr(iiter-1,1:6)
             EXIT ls_j
          END IF

          ! Update orbit:
          CALL getParameters(orb, &
               integration_step=orb_integration_step, &
               integrator=orb_integrator, &
               finite_diff=orb_finite_diff)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "TRACE BACK (69)", 1)
             RETURN
          END IF
          CALL NULLIFY(orb)
          err_verb_ = err_verb
          err_verb = err_verb! + 1
          IF (element_type_ls == "keplerian") THEN
             elements_iter_arr(iiter,4:6) = &
                  MODULO(elements_iter_arr(iiter,4:6), two_pi)
          ELSE IF (element_type_ls == "cometary") THEN
             elements_iter_arr(iiter,4:5) = &
                  MODULO(elements_iter_arr(iiter,4:5), two_pi)
          ELSE IF (element_type_ls == "cartesian" .AND. &
               SQRT(DOT_PRODUCT(elements_iter_arr(iiter,4:6), &
               elements_iter_arr(iiter,4:6))) > sol) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "Speed of object is larger than speed of light.", 1)
             elements_converge = .FALSE.
             elements_iter_arr(iiter,1:6) = elements_iter_arr(iiter-1,1:6)
             EXIT ls_j
          END IF
          CALL NEW(orb, elements_iter_arr(iiter,1:6), &
               TRIM(element_type_ls), TRIM(frame), copy(t))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "TRACE BACK (70)", 1)
             error = .FALSE.
             elements_converge = .FALSE.
             rmss_converge = .FALSE.
             elements_exist = .FALSE.
             elements_iter_arr(iiter,1:6) = elements_iter_arr(iiter-1,1:6)
             err_verb = err_verb_
             EXIT ls_j
          END IF
          elements_exist = .TRUE.
          err_verb = err_verb_
          CALL setParameters(orb, &
               dyn_model=dyn_model_, &
               perturbers=this%perturbers_prm, &
               integration_step=orb_integration_step, &
               integrator=orb_integrator, &
               finite_diff=orb_finite_diff)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "TRACE BACK (75)", 1)
             DEALLOCATE(obs_scoords, stat=err) 
             DEALLOCATE(ephemerides, stat=err) 
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             DEALLOCATE(cov_mat_obs, stat=err)
             DEALLOCATE(design_mat, stat=err) 
             DEALLOCATE(stdev_arr_obs, stat=err) 
             DEALLOCATE(elements_iter_arr, stat=err)
             DEALLOCATE(rmss_iter_arr, stat=err) 
             DEALLOCATE(inform_mat_obs_bd, stat=err) 
             DEALLOCATE(observed, stat=err)
             DEALLOCATE(computed, stat=err) 
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(obs_masks_, stat=err)
             DEALLOCATE(obsy_codes, stat=err)
             RETURN
          END IF

          ! Compute positions and partial derivatives:
          CALL getEphemerides(orb, obsy_ccoords, ephemerides, &
               partials_arr=partials_arr)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "TRACE BACK (80)", 1)
             error = .FALSE.
             elements_converge = .FALSE.
             elements_iter_arr(iiter,1:6) = elements_iter_arr(iiter-1,1:6)
             DEALLOCATE(ephemerides, partials_arr, stat=err)
             EXIT ls_j
          END IF
          CALL propagate(orb, t)

          ! Write design matrix and add cosine term to it and the computed positions:
          design_mat = 0.0_bp
          DO k=1,nobs
             computed(k,:) = getCoordinates(ephemerides(k))
             computed(k,2) = computed(k,2)*COS(observed(k,3))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / " // &
                     "leastSquares", &
                     "TRACE BACK (85)", 1)
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                RETURN
             END IF
             design_mat(2,k,:) = partials_arr(2,:,k)*COS(observed(k,3))
             design_mat(3,k,:) = partials_arr(3,:,k)
          END DO
          DEALLOCATE(ephemerides, partials_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / " // &
                  "leastSquares", &
                  "Could not deallocate memory (15).", 1)
             DEALLOCATE(obs_scoords, stat=err) 
             DEALLOCATE(ephemerides, stat=err) 
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             DEALLOCATE(cov_mat_obs, stat=err)
             DEALLOCATE(design_mat, stat=err) 
             DEALLOCATE(stdev_arr_obs, stat=err) 
             DEALLOCATE(elements_iter_arr, stat=err)
             DEALLOCATE(rmss_iter_arr, stat=err) 
             DEALLOCATE(inform_mat_obs_bd, stat=err) 
             DEALLOCATE(observed, stat=err)
             DEALLOCATE(computed, stat=err) 
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(obs_masks_, stat=err)
             DEALLOCATE(obsy_codes, stat=err)
             RETURN
          END IF

          ! Output orbital elements if needed:
          IF (info_verb >= 3) THEN
             t_ = getTime(orb)
             WRITE(stdout,"(1X,A,2(1X,I0,1X,A))") &
                  "Elements, residuals, and rms after major iteration", i, &
                  "and minor iteration", j, ":"
             SELECT CASE (element_type)
             CASE ("keplerian")
                elements = getElements(orb, element_type)
                WRITE(stdout,"(2X,A,6(1X,F20.13),1X,A)") "Kep: ", &
                     elements(1:2), elements(3:6)/rad_deg, &
                     getCalendarDateString(t_, "tdt")
             CASE ("cartesian")
                elements = getElements(orb, element_type, frame=frame)
                WRITE(stdout,"(2X,A,6(1X,F20.13),1X,A)") "Car: ", &
                     elements(1:6), &
                     getCalendarDateString(t_, "tdt")
             END SELECT
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / " // &
                     "leastSquares", &
                     "TRACE BACK (90)", 1)
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                RETURN
             END IF
             CALL NULLIFY(t_)
          END IF
          IF (info_verb >= 4) THEN 
             WRITE(stdout,"(2X,A)") "Covariance matrix:"
             CALL matrix_print(cov, stdout, errstr)
             errstr = ""
             WRITE(stdout,"(2X,A)") "Residuals RA & Dec [as]:"
          END IF

          ! Compute sky-plane residuals and rmss, and output them if needed:
          DO k=1,nobs
             residuals(k,1:6) = observed(k,1:6) - computed(k,1:6)
             IF (ABS(residuals(k,2)) > pi) THEN
                obs_ = observed(k,2)
                comp_ = computed(k,2)
                IF (obs_ < comp_) THEN
                   obs_ = obs_ + two_pi
                ELSE
                   comp_ = comp_ + two_pi
                END IF
                residuals(k,2) = obs_ - comp_
             END IF
             IF (info_verb >= 4) THEN
                WRITE(stdout,"(2X,2(F15.7,1X),A)") residuals(k,2:3)/rad_asec, TRIM(obsy_codes(k))
             END IF
          END DO
          rmss_iter_arr(iiter,2) = SQRT(SUM(residuals(:,2)**2.0_bp,mask=obs_masks_(:,2))/nobs_arr(2))
          rmss_iter_arr(iiter,3) = SQRT(SUM(residuals(:,3)**2.0_bp,mask=obs_masks_(:,3))/nobs_arr(3))
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(2X,A,2(1X,F15.7))") "RA & Dec RMS [as]: ", &
                  rmss_iter_arr(iiter,2)/rad_asec, rmss_iter_arr(iiter,3)/rad_asec
          END IF

          ! Exit iteration loop if the difference between the elements of
          ! the new and the old orbit are sufficiently small:
          DO k=1,SIZE(cov,dim=1)
             IF (cov(k,k) < 0.0_bp) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / " // &
                     "leastSquares", &
                     "Negative diagonal in covariance matrix:", 1)
                IF (err_verb >= 1) THEN
                   CALL matrix_print(cov,stderr,errstr)
                   errstr = ""
                   WRITE(stderr,*) "Orbital elements:"
                   WRITE(stderr,*) elements_iter_arr(iiter,1:6)                  
                END IF
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                RETURN
             END IF
             stdev(k) = SQRT(cov(k,k))
          END DO
          elem_diff = ABS(elements_iter_arr(iiter,1:6) - &
               elements_iter_arr(iiter-1,1:6)) / correction_factor
          IF (ALL(elem_diff < stdev) .AND. j >= 3) THEN
             elements_converge = .TRUE.
          END IF

          IF (elements_converge .AND. approach == 11 .AND. &
               ((ALL(rmss_iter_arr(iiter,2:3) < &
               sigma_multiplier_rms*stdev_max(2:3))) .OR. &
               (j > 1 .AND. ALL(ABS(rmss_iter_arr(iiter,2:3) - &
               rmss_iter_arr(iiter-1,2:3)) < &
               1000.0_bp*EPSILON(rmss_iter_arr(iiter,2:3)))))) THEN
             IF (info_verb >= 3) THEN
                WRITE(stdout,"(1X,A,I0)") "Exiting @", j
             END IF
             EXIT ls_j
          ELSE IF (elements_converge .AND. approach /= 11) THEN
             IF (info_verb >= 3) THEN
                WRITE(stdout,"(1X,A,I0)") "Exiting @", j
             END IF
             EXIT ls_j
          END IF

       END DO ls_j

!!$       IF (iiter >= 3) THEN
!!$          IF (elements_exist .AND. &
!!$               ((SQRT(SUM((rmss_iter_arr(iiter-2,:) - rmss_iter_arr(iiter-1,:))**2.0_bp)) >= &
!!$               SQRT(SUM((rmss_iter_arr(iiter-1,:) - rmss_iter_arr(iiter,:))**2.0_bp)) .OR. &
!!$               SQRT(SUM((rmss_iter_arr(iiter-1,:) - rmss_iter_arr(iiter,:))**2.0_bp)) < &
!!$               0.000001_bp*SQRT(SUM((rmss_iter_arr(iiter,:))**2.0_bp))) )) THEN
!!$             rmss_converge = .TRUE.
!!$          ELSE
!!$             rmss_converge = .FALSE.                
!!$          END IF
!!$       END IF

       IF (elements_converge .AND. rmss_converge .AND. &
            final_iterations == this%ls_niter_major_min_prm) THEN
          IF (info_verb >= 2) THEN
             WRITE(stdout,*) "LS solution found (1)."
          END IF
          IF (this%outlier_rejection_prm) THEN
             DO k=1,6
                nobs_arr(k) = COUNT(obs_masks_(:,k))
             END DO
             IF (ANY(nobs_arr(2:3) < 0.2_bp*nobs)) THEN
                IF (info_verb >= 2) THEN
                   WRITE(*,*) "WARNING from LS: More than 20% of the observations excluded!", &
                        " Check the assumptions used for outlier rejection (e.g., noise)."
                END IF
             END IF
             IF (ALL(nobs_arr < 4)) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / leastSquares", &
                     "Least squares approach not suitable when " // &
                     "possible outliers have been excluded.", 1)
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                RETURN
             END IF
             IF (ANY(rmss_iter_arr(iiter,2:3) > sigma_multiplier_rms*stdev_max(2:3))) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / leastSquares", &
                     "Least squares solution not reliable when " // &
                     "rms values are so large as compared to assumed " // &
                     "observational uncertainty. (Check noise assumption!)", 1)
                RETURN
             END IF
          END IF
          EXIT ls_i
       ELSE IF (elements_converge .AND. rmss_converge .AND. approach == 11) THEN
          final_iterations = final_iterations + 1
          IF (this%outlier_rejection_prm .AND. &
               ALL(rmss_iter_arr(iiter,2:3) <= sigma_multiplier_rms*stdev_max(2:3))) THEN
             obs_masks_ = this%obs_masks_prm
             DO k=1,nobs
                ! The following assumes that the stdevs are equal for
                ! both coordinates:
                IF (SQRT(SUM(residuals(k,2:3)**2)) > &
                     this%outlier_multiplier_prm*MAXVAL(stdev_arr_obs(k,2:3))) THEN
                   obs_masks_(k,:) = .FALSE.
                END IF
             END DO
             DO k=1,6
                nobs_arr(k) = COUNT(obs_masks_(:,k))
             END DO
             IF (ANY(nobs_arr(2:3) < 0.2*nobs)) THEN
                IF (info_verb >= 2) THEN
                   WRITE(stdout,*) "More than 20% of the observations excluded! Keeping all."
                END IF
                obs_masks_ = this%obs_masks_prm
                DO k=1,6
                   nobs_arr(k) = COUNT(obs_masks_(:,k))
                END DO
             END IF
             IF (ALL(nobs_arr < 4)) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / leastSquares", &
                     "Least squares approach not suitable when " // &
                     "possible outliers have been excluded.", 1)
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                RETURN
             END IF
          END IF
       ELSE IF (elements_converge .AND. rmss_converge .AND. approach /= 11) THEN
          final_iterations = 0
          approach_direction(1:2) = approach_direction(2:3)
          approach_direction(3) = -1
          approach = approach - 1
          IF (approach == 40) THEN
             approach = 12
          END IF
       ELSE IF (elements_converge .AND. .NOT.rmss_converge) THEN
          final_iterations = 0
          approach_direction(1:2) = approach_direction(2:3)
          approach_direction(3) = 1
          IF (approach < 21) THEN
             approach = 21
          ELSE
             approach = approach + 1
          END IF
       ELSE IF (.NOT.elements_converge .AND. rmss_converge) THEN
          final_iterations = 0
          approach_direction(1:2) = approach_direction(2:3)
          approach_direction(3) = 1
          approach = approach + 1
       ELSE IF (.NOT.elements_converge .AND. .NOT.rmss_converge) THEN
          final_iterations = 0
          CALL NULLIFY(orb)
          orb = copy(preliminary_orbit)
          DO k=1,iiter
             elements_iter_arr(k,1:6) = elements_pre(1:6)
          END DO
          DO k=1,iiter
             rmss_iter_arr(k,1:6) = rms_pre(1:6)
          END DO
          IF (approach >= 41) THEN
             approach_direction(1:2) = approach_direction(2:3)
             approach_direction(3) = 1
             approach = approach + 1
          ELSE IF (approach == 11) THEN
             approach = 12
          ELSE IF (approach < 41) THEN
             approach = 41
          END IF
       END IF
       !IF (approach > 48) THEN
       IF (approach >= 43) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / leastSquares", &
               "No converging solution found (5).", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(ephemerides, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(design_mat, stat=err) 
          DEALLOCATE(stdev_arr_obs, stat=err) 
          DEALLOCATE(elements_iter_arr, stat=err)
          DEALLOCATE(rmss_iter_arr, stat=err) 
          DEALLOCATE(inform_mat_obs_bd, stat=err) 
          DEALLOCATE(observed, stat=err)
          DEALLOCATE(computed, stat=err) 
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(obs_masks_, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          RETURN
       END IF
       IF (approach_direction(3) < 0 .AND. &
            i > 25 .AND. approach < 44) THEN
          IF (info_verb >= 2) THEN
             WRITE(stdout,*) "LS solution found (2)."
          END IF
          EXIT ls_i
       END IF

       DO
          SELECT CASE (approach)
          CASE (11)
             dyn_model_ = dyn_model
             correction_factor = this%ls_corr_fac_prm
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (12)
             dyn_model_ = dyn_model
             correction_factor = this%ls_corr_fac_prm
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (21)
             dyn_model_ = dyn_model
             correction_factor = this%ls_corr_fac_prm
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (41)
             dyn_model_ = dyn_model
             correction_factor = 0.1_bp*this%ls_corr_fac_prm
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (42)
             dyn_model_ = dyn_model
             correction_factor = 0.01_bp*this%ls_corr_fac_prm
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (43)
             dyn_model_ = dyn_model
             correction_factor = 0.001_bp*this%ls_corr_fac_prm
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (44)
             dyn_model_ = dyn_model
             correction_factor = 0.0001_bp*this%ls_corr_fac_prm
             IF (element_type_ls /= element_type) THEN
                element_type_ls = element_type
                IF (element_type_ls == "keplerian") THEN
                   CALL toKeplerian(orb)
                ELSE IF (element_type_ls == "cartesian") THEN
                   CALL toCartesian(orb, frame=frame)
                END IF
                elements_iter_arr(iiter,1:6) = getElements(orb, element_type_ls, frame=frame)
             END IF
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / " // &
                     "leastSquares", &
                     "TRACE BACK (105)", 1)
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                RETURN
             END IF
             element_mask(1:6) = this%ls_elem_mask_prm
             EXIT
          CASE (45)
             dyn_model_ = dyn_model
             correction_factor = 0.001_bp*this%ls_corr_fac_prm
             IF (element_type_ls == "keplerian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1) = .FALSE.
             ELSE IF (element_type_ls == "cartesian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1) = .FALSE.
             END IF
             EXIT
          CASE (46)
             dyn_model_ = dyn_model
             correction_factor = 0.001_bp*this%ls_corr_fac_prm
             IF (element_type_ls == "keplerian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(2) = .FALSE.
             ELSE IF (element_type_ls == "cartesian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1) = .FALSE.
                element_mask(4) = .FALSE.
             END IF
             EXIT
          CASE (47)
             dyn_model_ = dyn_model
             correction_factor = 0.001_bp*this%ls_corr_fac_prm
             IF (element_type_ls == "keplerian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1:2) = .FALSE.
             ELSE IF (element_type_ls == "cartesian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1) = .FALSE.
                element_mask(4:5) = .FALSE.
             END IF
             EXIT
          CASE (48)
             dyn_model_ = dyn_model
             correction_factor = 0.0001_bp*this%ls_corr_fac_prm
             IF (element_type_ls == "keplerian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1:2) = .FALSE.
             ELSE IF (element_type_ls == "cartesian") THEN
                element_mask(1:6) = .TRUE.
                element_mask(1) = .FALSE.
                element_mask(4:5) = .FALSE.
             END IF
             EXIT
          CASE default
             IF (approach > 48) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / leastSquares", &
                     "No converging solution found (10).", 1)
                DEALLOCATE(obs_scoords, stat=err) 
                DEALLOCATE(ephemerides, stat=err) 
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(cov_mat_obs, stat=err)
                DEALLOCATE(design_mat, stat=err) 
                DEALLOCATE(stdev_arr_obs, stat=err) 
                DEALLOCATE(elements_iter_arr, stat=err)
                DEALLOCATE(rmss_iter_arr, stat=err) 
                DEALLOCATE(inform_mat_obs_bd, stat=err) 
                DEALLOCATE(observed, stat=err)
                DEALLOCATE(computed, stat=err) 
                DEALLOCATE(residuals, stat=err) 
                DEALLOCATE(obs_masks_, stat=err)
                DEALLOCATE(obsy_codes, stat=err)
                WRITE(stderr,*) "Returning from leastSquares"
                RETURN
             ELSE
                approach = approach + approach_direction(3)
             END IF
          END SELECT
       END DO
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,L1)") "Elements converge: ", elements_converge
          WRITE(stdout,"(2X,A,L1)") "RMSs converge: ", rmss_converge
          WRITE(stdout,"(2X,A,I0)") "Approach to be used: ", approach
          WRITE(stdout,"(2X,A,3(I0,1X))") "Approach directions: ", approach_direction
          WRITE(stdout,"(2X,A,6(I0,1X))") "Observation included: ",nobs_arr 
       END IF

    END DO ls_i

    ! Orbital elements and sky-plane residuals and rms":
    IF (info_verb >= 2) THEN
       t_ = getTime(orb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "leastSquares", &
               "TRACE BACK (130)", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(ephemerides, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(design_mat, stat=err) 
          DEALLOCATE(stdev_arr_obs, stat=err) 
          DEALLOCATE(elements_iter_arr, stat=err)
          DEALLOCATE(rmss_iter_arr, stat=err) 
          DEALLOCATE(inform_mat_obs_bd, stat=err) 
          DEALLOCATE(observed, stat=err)
          DEALLOCATE(computed, stat=err) 
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(obs_masks_, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          RETURN
       END IF
       WRITE(stdout,"(2X,A)") "# Final results:"
       err_verb_ = err_verb
       err_verb = 0
       elements = getElements(orb, "keplerian")
       err_verb = err_verb_
       IF (error) THEN
          error = .FALSE.
       ELSE
          str = getCalendarDateString(t_, "TT")
          WRITE(stdout,"(2X,A,2X,A,6(1X,F20.13),1X,A)") "#", "Kep: ", &
               elements(1:2), elements(3:6)/rad_deg, &
               TRIM(str)
       END IF
       elements = getElements(orb, "cartesian", frame="ecliptic")
       IF (error) THEN
          error = .FALSE.
       ELSE
          str = getCalendarDateString(t_, "TT")
          WRITE(stdout,"(2X,A,2X,A,6(1X,F20.13),1X,A)") "#", "Car: ", &
               elements(1:6), TRIM(str)
       END IF
       CALL NULLIFY(t_)
       WRITE(stdout,"(2X,A,2X,A)") "#", "Residuals RA & Dec [as]:"
       DO j=1,nobs
          IF (ALL(obs_masks_(j,2:3))) THEN
             WRITE(stdout,"(2X,A,2X,A,2(F15.7,1X),A,1X,A)") "#", " ", &
                  residuals(j,2:3)/rad_asec, " ", TRIM(obsy_codes(j))
          ELSE
             WRITE(stdout,"(2X,A,2X,A,2(F15.7,1X),A,1X,A)") "#", "(", &
                  residuals(j,2:3)/rad_asec, ")", TRIM(obsy_codes(j))
          END IF
       END DO
       WRITE(stdout,"(2X,A,2X,A,2(1X,F15.7))") "#", "RA & Dec RMS [as]: ", &
            rmss_iter_arr(iiter,2)/rad_asec, rmss_iter_arr(iiter,3)/rad_asec
    END IF

    CALL NULLIFY(this%orb_ml_cmp)
    this%orb_ml_cmp = copy(orb)
    IF (.NOT.ASSOCIATED(this%cov_ml_cmp)) THEN
       ALLOCATE(this%cov_ml_cmp(6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / " // &
               "leastSquares", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(ephemerides, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(design_mat, stat=err) 
          DEALLOCATE(stdev_arr_obs, stat=err) 
          DEALLOCATE(elements_iter_arr, stat=err)
          DEALLOCATE(rmss_iter_arr, stat=err) 
          DEALLOCATE(inform_mat_obs_bd, stat=err) 
          DEALLOCATE(observed, stat=err)
          DEALLOCATE(computed, stat=err) 
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(obs_masks_, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          RETURN
       END IF
    END IF
    this%cov_ml_cmp = cov
    this%cov_type_prm = element_type_ls
    this%obs_masks_prm = obs_masks_

    CALL NULLIFY(orb)
    CALL NULLIFY(t)
    CALL NULLIFY(t_)
    DEALLOCATE(obs_scoords, stat=err) 
    DEALLOCATE(ephemerides, stat=err) 
    DEALLOCATE(obsy_ccoords, stat=err)
    DEALLOCATE(partials_arr, stat=err)
    DEALLOCATE(cov_mat_obs, stat=err)
    DEALLOCATE(design_mat, stat=err) 
    DEALLOCATE(stdev_arr_obs, stat=err) 
    DEALLOCATE(elements_iter_arr, stat=err)
    DEALLOCATE(rmss_iter_arr, stat=err) 
    DEALLOCATE(inform_mat_obs_bd, stat=err) 
    DEALLOCATE(observed, stat=err)
    DEALLOCATE(computed, stat=err) 
    DEALLOCATE(residuals, stat=err)
    DEALLOCATE(obs_masks_, stat=err)
    DEALLOCATE(obsy_codes, stat=err)

    IF (final_iterations /= this%ls_niter_major_min_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "leastSquares", &
            "Could not reach final iteration rounds.", 1)
    END IF

  END SUBROUTINE leastSquares_SO





  !! *Description*:
  !!
  !! Computation of the least-squares orbital elements (Cartesian or
  !! Keplerian) and their covariances in the two-body and full
  !! many-body approaches.
  !!
  !! Returns error.
  !!
  SUBROUTINE levenbergMarquardt_SO(this, preliminary_orbit)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), INTENT(in)              :: preliminary_orbit    

    INTEGER, PARAMETER :: niter_size = 3
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: obs_scoords
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords
    TYPE (Time) :: t, t_
    TYPE (Orbit) :: orb
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: obsy_codes
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type_
    CHARACTER(len=FRAME_LEN) :: frame_
    CHARACTER(len=DYN_MODEL_LEN) :: dyn_model_
    CHARACTER(len=INTEGRATOR_LEN) :: integrator_
    CHARACTER(len=32) :: str
    REAL(bp), DIMENSION(:,:,:), POINTER :: information_matrix_measur
    REAL(bp), DIMENSION(:,:), POINTER :: stdev_arr_measur
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE :: jacobians
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: measur, & ! incl. cos(dec)
         residuals, &
         alpha
    REAL(bp), DIMENSION(6,6) :: cov_mat_param
    REAL(bp), DIMENSION(niter_size,6) :: elements_iter_arr
    REAL(bp), DIMENSION(6) :: params, finite_diff_
    REAL(bp) :: rchi2, integration_step_, rchi2_previous, rchi2_old, lambda
    INTEGER :: i, j, k, err, ndata, nparam, nmultidata
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: mask_measur
    LOGICAL, DIMENSION(6) :: mask_param

    REAL(bp), PARAMETER :: rchi2_max_prm = 1.5_bp

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "LEAST SQUARES"
    END IF

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(preliminary_orbit)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "Preliminary orbit has not been initialized.", 1)
       RETURN
    END IF

    IF (this%outlier_rejection_prm) THEN
       IF (this%outlier_multiplier_prm < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
               "Outlier criterion has not been defined.", 1)
          RETURN
       END IF
    END IF

    IF (.NOT.ASSOCIATED(this%obs_masks_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "Observation mask is missing.", 1)
       RETURN
    END IF

    IF (this%ls_niter_minor_prm <= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "ls_niter_minor_prm undefined.", 1)
       RETURN       
    END IF

    IF (this%ls_niter_major_max_prm <= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "ls_niter_major_max_prm undefined.", 1)
       RETURN       
    END IF

    ! Allocate memory:
    ndata = SIZE(this%obs_masks_prm,dim=1)
    nmultidata = SIZE(this%obs_masks_prm,dim=2)
    nparam = 6
    ALLOCATE(measur(ndata,nmultidata), residuals(ndata,nmultidata), &
         alpha(nparam,nparam), jacobians(nmultidata,nparam,ndata), &
         mask_measur(ndata,nmultidata), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    orb = copy(preliminary_orbit)

    ! Inversion epoch is equal to epoch of initial orbit:
    t = getTime(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF

    ! Find parameters of the peliminary orbit:
    CALL getParameters(orb, dyn_model=dyn_model_)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    IF (dyn_model_ /= this%dyn_model_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "Inconsistent propagation schemes: " // &
            "orb=" // TRIM(dyn_model_) // " and storb=" // &
            TRIM(this%dyn_model_prm) // ".", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    CALL getParameters(orb, &
         integration_step=integration_step_, &
         integrator=integrator_, &
         finite_diff=finite_diff_)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (15)", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF

    ! Initialize elements:
    frame_ = getFrame(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (20)", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    element_type_ = this%element_type_prm
    CALL locase(element_type_, error)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "The element type string contains forbidden characters.", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    SELECT CASE (TRIM(element_type_))
    CASE ("cartesian")
       CALL toCartesian(orb, frame=frame_)
    CASE ("keplerian")
       CALL toKeplerian(orb)
    CASE ("cometary")
       CALL toCometary(orb)
    CASE default
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / levenbergMarquardt", &
            "Can not use elements of type: " // TRIM(element_type_), 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END SELECT
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (25)", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    params(1:6) = getElements(orb, element_type_)
    DO i=1,niter_size
       elements_iter_arr(i,1:6) = params(1:6)
    END DO
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (30)", 1)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(1X,3(1X,A,1X,A))") "Using", TRIM(element_type_), &
            "elements."
    END IF

    ! Observations and observers:
    obs_scoords => getObservationSCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (35)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    mask_measur = this%obs_masks_prm
    obsy_codes => getObservatoryCodes(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (40)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(obsy_codes, stat=err)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    obsy_ccoords => getObservatoryCCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (45)", 1)
       DEALLOCATE(obs_scoords, stat=err) 
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    DO i=1,ndata
       measur(i,:) = getCoordinates(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "levenbergMarquardt", &
               "TRACE BACK (50)", 1)
          DEALLOCATE(obs_scoords, stat=err) 
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(obsy_codes, stat=err)
          DEALLOCATE(measur, stat=err)
          DEALLOCATE(residuals, stat=err) 
          DEALLOCATE(mask_measur, stat=err)
          DEALLOCATE(alpha, stat=err)
          DEALLOCATE(jacobians, stat=err)
          RETURN
       END IF
    END DO
    measur(:,2) = measur(:,2)*COS(measur(:,3))
    DEALLOCATE(obs_scoords, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "Could not deallocate memory (5).", 1)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(obsy_codes, stat=err)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    information_matrix_measur => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (55)", 1)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(information_matrix_measur, stat=err) 
       DEALLOCATE(obsy_codes, stat=err)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    stdev_arr_measur => getStandardDeviations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (60)", 1)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(information_matrix_measur, stat=err) 
       DEALLOCATE(obsy_codes, stat=err)
       DEALLOCATE(stdev_arr_measur, stat=err)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF
    mask_param = this%ls_elem_mask_prm

    rchi2_old = HUGE(rchi2_old)
    rchi2_previous = HUGE(rchi2_previous)
    lambda = -1.0_bp
    DO i=1,this%ls_niter_major_max_prm

       DO j=1,this%ls_niter_minor_prm
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(2X,A)") "Call to 'LevenbergMarquardt_private' from " // &
                  "loop in 'levenbergMarquardt_SO'."
          END IF
          CALL LevenbergMarquardt_private
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "levenbergMarquardt", &
                  "Could not find a least-squares solution (1).", 1)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(information_matrix_measur, stat=err) 
             DEALLOCATE(obsy_codes, stat=err)
             DEALLOCATE(stdev_arr_measur, stat=err)
             DEALLOCATE(measur, stat=err)
             DEALLOCATE(residuals, stat=err) 
             DEALLOCATE(mask_measur, stat=err)
             DEALLOCATE(alpha, stat=err)
             DEALLOCATE(jacobians, stat=err)
             RETURN
          END IF
          IF (ABS(rchi2 - rchi2_previous) < this%ls_rchi2_diff_tresh_prm .AND. &
               rchi2 <= rchi2_previous) THEN
             !             rchi2_previous = rchi2
             EXIT
             !          ELSE
             !             rchi2_previous = rchi2
          END IF
       END DO

       ! Outlier rejection
       IF (this%outlier_rejection_prm .AND. ABS(rchi2 - rchi2_old) > this%ls_rchi2_diff_tresh_prm) THEN
          mask_measur = this%obs_masks_prm
          DO k=1,ndata
             IF (ANY(ABS(residuals(k,2:3)) > &
                  this%outlier_multiplier_prm*stdev_arr_measur(k,2:3))) THEN
                mask_measur(k,:) = .FALSE.
             END IF
          END DO
          rchi2_old = rchi2
       ELSE IF (.NOT.this%outlier_rejection_prm) THEN
          EXIT
       END IF

    END DO

    lambda = 0.0_bp
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "Final call to 'LevenbergMarquardt_private' from " // &
            "'levenbergMarquardt_SO'."
    END IF
    CALL levenbergMarquardt_private
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "Could not find a least-squares solution (2).", 1)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(information_matrix_measur, stat=err) 
       DEALLOCATE(obsy_codes, stat=err)
       DEALLOCATE(stdev_arr_measur, stat=err)
       DEALLOCATE(measur, stat=err)
       DEALLOCATE(residuals, stat=err) 
       DEALLOCATE(mask_measur, stat=err)
       DEALLOCATE(alpha, stat=err)
       DEALLOCATE(jacobians, stat=err)
       RETURN
    END IF

    IF (element_type_ == "keplerian") THEN
       params(4:6) = MODULO(params(4:6), two_pi)
    ELSE IF (element_type_ == "cometary") THEN
       params(4:5) = MODULO(params(4:5), two_pi)
    ELSE IF (element_type_ == "cartesian" .AND. &
         SQRT(DOT_PRODUCT(params(4:6),params(4:6))) > sol) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "Speed of object is larger than speed of light.", 1)
       DEALLOCATE(mask_measur, stat=err)
       RETURN
    END IF

    ! Output final orbital elements and sky-plane residuals and rms if needed:
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "Call to 'ephemeris_lsl' from " // &
            "'levenbergMarquardt_SO'."
    END IF
    CALL ephemeris_lsl(params, residuals, jacobians, rchi2)

    DEALLOCATE(obsy_ccoords, stat=err)
    DEALLOCATE(information_matrix_measur, stat=err) 
    DEALLOCATE(obsy_codes, stat=err)
    DEALLOCATE(stdev_arr_measur, stat=err)
    DEALLOCATE(measur, stat=err)
    DEALLOCATE(residuals, stat=err) 
    DEALLOCATE(alpha, stat=err)
    DEALLOCATE(jacobians, stat=err)

    CALL NULLIFY(this%orb_ml_cmp)
    CALL NEW(this%orb_ml_cmp, params, TRIM(element_type_), TRIM(frame_), copy(t))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (65)", 1)
       DEALLOCATE(mask_measur, stat=err)
       RETURN
    END IF
    CALL setParameters(this%orb_ml_cmp, &
         dyn_model=dyn_model_, &
         perturbers=this%perturbers_prm, &
         integration_step=integration_step_, &
         integrator=integrator_, &
         finite_diff=finite_diff_)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "TRACE BACK (70)", 1)
       DEALLOCATE(mask_measur, stat=err)
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(this%cov_ml_cmp)) THEN
       ALLOCATE(this%cov_ml_cmp(6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / " // &
               "leastSquares", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(mask_measur, stat=err)
          RETURN
       END IF
    END IF
    DO i=1,nparam
       IF (.NOT.mask_param(i)) THEN
          cov_mat_param(i,i) = 0.0_bp
       END IF
    END DO
    this%cov_ml_cmp = cov_mat_param
    this%cov_type_prm = element_type_
    this%obs_masks_prm = mask_measur
    DEALLOCATE(mask_measur, stat=err)

    ! Check whether acceptable solution based on rchi2
    IF (rchi2 > rchi2_max_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "levenbergMarquardt", &
            "Could not find an acceptable least-squares solution.", 1)
       RETURN
    END IF

  CONTAINS

    ! LevenbergMarquardt_private and coefficients are essentially
    ! identical to the subroutines found in estimators.f90.

    SUBROUTINE levenbergMarquardt_private

      !      REAL(bp), SAVE :: rchi2_
      REAL(bp), DIMENSION(:), ALLOCATABLE, SAVE :: params_, beta
      REAL(bp), DIMENSION(:,:), ALLOCATABLE, SAVE :: param_corrections

      IF (lambda < 0.0_bp) THEN
         ALLOCATE(params_(nparam),beta(nparam),param_corrections(nparam,1))
         !lambda = 0.001_bp
         lambda = 0.000001_bp
         params_ = params
         IF (info_verb >= 2) THEN
            WRITE(stdout,"(2X,A)") "First call to 'coefficients' from " // &
                 "'levenbergMarquardt_private' when lambda < 0."
         END IF
         CALL coefficients(params_, alpha, beta)
         IF (error) THEN
            CALL errorMessage("StochasticOrbit / " // &
                 "levenbergMarquardt / levenbergMarquardt_private", &
                 "TRACE BACK (5)", 1)
            DEALLOCATE(params_, beta, param_corrections)
            RETURN
         END IF
         !         rchi2_ = rchi2
         rchi2_previous = rchi2
      ELSE
         rchi2_previous = rchi2         
      END IF
      cov_mat_param = alpha
      cov_mat_param = diagonal_multiplication(cov_mat_param, 1.0_bp+lambda, errstr, mask_param)
      IF (LEN_TRIM(errstr) /= 0) THEN
         error = .TRUE.
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / levenbergMarquardt_private", &
              TRIM(errstr), 1)
         errstr = ""
         DEALLOCATE(params_, beta, param_corrections)
         RETURN
      END IF
      param_corrections(:,1) = beta
      CALL gauss_jordan(cov_mat_param, param_corrections, errstr)
      IF (LEN_TRIM(errstr) /= 0) THEN
         error = .TRUE.
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / levenbergMarquardt_private", &
              TRIM(errstr), 1)
         errstr = ""
         DEALLOCATE(params_, beta, param_corrections)
         RETURN
      END IF
      IF (ABS(lambda) < EPSILON(lambda)) THEN
         DEALLOCATE(params_, beta, param_corrections)
         RETURN
      END IF
      params_ = params + param_corrections(:,1)
      IF (info_verb >= 2) THEN
         WRITE(stdout,"(2X,A)") "Second call to 'coefficients' from " // &
              "'levenbergMarquardt_private' when lambda > 0."
      END IF
      CALL coefficients(params_, cov_mat_param, param_corrections(:,1))
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / levenbergMarquardt_private", &
              "TRACE BACK (10)", 1)
         DEALLOCATE(params_, beta, param_corrections)
         RETURN
      END IF
      !      IF (rchi2 < rchi2_) THEN
      IF (rchi2 < rchi2_previous) THEN
         lambda = 0.5_bp*lambda
         !         rchi2_ = rchi2
         alpha = cov_mat_param
         beta = param_corrections(:,1)
         params = params_
      ELSE
         !lambda = 10.0_bp*lambda
         lambda = 10.0_bp*lambda
         !         rchi2 = rchi2_
      END IF
      IF (lambda < EPSILON(lambda)) THEN
         lambda = EPSILON(lambda)
      END IF

    END SUBROUTINE levenbergMarquardt_private


    SUBROUTINE coefficients(params, alpha, beta)

      !IMPLICIT NONE
      REAL(bp), DIMENSION(:), INTENT(inout) :: params
      REAL(bp), DIMENSION(:), INTENT(out) :: beta
      REAL(bp), DIMENSION(:,:), INTENT(out) :: alpha

      REAL(bp), DIMENSION(nparam,nmultidata) :: tmp
      INTEGER :: i

      IF (info_verb >= 2) THEN
         WRITE(stdout,"(2X,A)") "Call to 'ephemeris_lsl' from " // &
              "'coefficients'."
      END IF
      CALL ephemeris_lsl(params, residuals, jacobians, rchi2)
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / coefficients", &
              "TRACE BACK (5)", 1)
         RETURN
      END IF
      ! Approximate Hessian by multiplying Jacobians
      ! alpha = cov_param^(-1) = J^T Sigma_obs^(-1) J:
      ! beta = J^T Sigma_obs^(-1) y:
      alpha = 0.0_bp
      beta = 0.0_bp
      DO i=1,ndata
         tmp = MATMUL(TRANSPOSE(jacobians(1:nmultidata,1:nparam,i)), &
              information_matrix_measur(i,1:nmultidata,1:nmultidata))
         alpha = alpha + MATMUL(tmp, jacobians(1:nmultidata,1:nparam,i))
         beta = beta + MATMUL(tmp, residuals(i,1:nmultidata))
      END DO
      DO i=1,nparam
         IF (.NOT.mask_param(i)) THEN
            alpha(i,:) = 0.0_bp
            alpha(:,i) = 0.0_bp
            alpha(i,i) = 1.0_bp
            beta(i) = 0.0_bp
         END IF
      END DO

    END SUBROUTINE coefficients


    SUBROUTINE ephemeris_lsl(elements, residuals, jacobians, rchi2)

      !implicit none
      REAL(bp), DIMENSION(:), INTENT(inout) :: elements ! nparam
      REAL(bp), DIMENSION(:,:), INTENT(out) :: residuals ! ndata,nmultidata
      REAL(bp), DIMENSION(:,:,:), INTENT(out) :: jacobians ! nmultidata,nparam,ndata
      REAL(bp), INTENT(out) :: rchi2

      TYPE (Orbit) :: orb
      TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: ephemerides
      REAL(bp), DIMENSION(:,:,:), POINTER :: partials_arr
      REAL(bp), DIMENSION(SIZE(residuals,dim=1),SIZE(residuals,dim=2)) :: computed
      INTEGER :: j

      IF (element_type_ == "keplerian") THEN
         elements(4:6) = MODULO(elements(4:6), two_pi)
      ELSE IF (element_type_ == "cometary") THEN
         elements(4:5) = MODULO(elements(4:5), two_pi)
      ELSE IF (element_type_ == "cartesian" .AND. &
           SQRT(DOT_PRODUCT(elements(4:6),elements(4:6))) > sol) THEN
         error = .TRUE.
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / ephemeris_lsl", &
              "Speed of object is larger than speed of light.", 1)
         RETURN
      END IF
      CALL NEW(orb, elements, TRIM(element_type_), TRIM(frame_), copy(t))
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / ephemeris_lsl", &
              "TRACE BACK (5)", 1)
         RETURN
      END IF
      CALL setParameters(orb, &
           dyn_model=dyn_model_, &
           perturbers=this%perturbers_prm, &
           integration_step=integration_step_, &
           integrator=integrator_, &
           finite_diff=finite_diff_)
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / ephemeris_lsl", &
              "TRACE BACK (10)", 1)
         RETURN
      END IF

      ! Compute positions and partial derivatives:
      CALL getEphemerides(orb, obsy_ccoords, ephemerides, &
           partials_arr=partials_arr)
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / ephemeris_lsl", &
              "TRACE BACK (15)", 1)
         RETURN
      END IF

      ! Write design matrix and add cosine term to it and the computed positions:
      jacobians = 0.0_bp
      DO j=1,ndata
         computed(j,:) = getCoordinates(ephemerides(j))
         computed(j,2) = computed(j,2)*COS(measur(j,3))
         IF (error) THEN
            CALL errorMessage("StochasticOrbit / " // &
                 "levenbergMarquardt / ephemeris_lsl", &
                 "TRACE BACK (20)", 1)
            RETURN
         END IF
         jacobians(:,:,j) = partials_arr(:,:,j)
         jacobians(2,:,j) = jacobians(2,:,j)*COS(measur(j,3))
      END DO
      DEALLOCATE(ephemerides, partials_arr, stat=err)
      IF (err /= 0) THEN
         error = .TRUE.
         CALL errorMessage("StochasticOrbit / " // &
              "levenbergMarquardt / ephemeris_lsl", &
              "Could not deallocate memory.", 1)
         RETURN
      END IF
      ! Output orbital elements and sky-plane residuals and rms if needed:
      IF (info_verb >= 2) THEN
         t_ = getTime(orb)
         WRITE(stdout,"(2X,A,I0,A)") "Elements, residuals, RMS, and reduced chi2:"
         str = getCalendarDateString(t_, "TT")
         SELECT CASE (element_type_)
         CASE ("keplerian")
            elements = getElements(orb, element_type_)
            WRITE(stdout,"(2X,A,6(1X,F20.13),1X,A)") "Kep: ", &
                 elements(1:2), elements(3:6)/rad_deg, &
                 TRIM(str)
         CASE ("cartesian")
            elements = getElements(orb, element_type_, frame=frame_)
            WRITE(stdout,'(2X,A,6(1X,F20.13),1X,A)') "Car: ", &
                 elements(1:6),TRIM(str)
         END SELECT
         IF (error) THEN
            CALL errorMessage("StochasticOrbit / " // &
                 "levenbergMarquardt / ephemeris_lsl", &
                 "TRACE BACK (25)", 1)
            RETURN
         END IF
         CALL NULLIFY(t_)
         WRITE(stdout,"(2X,A)") "Residuals RA & Dec [as]:"
      END IF
      CALL NULLIFY(orb)
      residuals = 0.0_bp
      DO j=1,ndata
         residuals(j,1:6) = measur(j,1:6) - computed(j,1:6)
         IF (ABS(residuals(j,2)) > pi) THEN
            residuals(j,2) = two_pi - residuals(j,2)
         END IF
         IF (info_verb >= 2) THEN
            IF (ALL(mask_measur(j,2:3))) THEN
               WRITE(stdout,"(2X,A,2(F20.7,1X),A,1X,A)") " ", &
                    residuals(j,2:3)/rad_asec, " ", TRIM(obsy_codes(j))
            ELSE
               WRITE(stdout,"(2X,A,2(F15.7,1X),A,1X,A)") "(", &
                    residuals(j,2:3)/rad_asec, ")", TRIM(obsy_codes(j))
            END IF
         END IF
      END DO

      ! Compute reduced chi2:
      rchi2 = chi_square(residuals, information_matrix_measur, mask_measur, errstr) / &
           REAL(COUNT(mask_measur)-COUNT(mask_param),bp)

      IF (info_verb >= 2) THEN
         WRITE(stdout,"(2X,A,2(1X,F15.7))") "RMS RA & Dec [arcsec]: ", &
              SQRT(SUM(residuals(:,2)**2,mask=mask_measur(:,2))/COUNT(mask_measur(:,2)))/rad_asec, &
              SQRT(SUM(residuals(:,3)**2,mask=mask_measur(:,3))/COUNT(mask_measur(:,3)))/rad_asec
         WRITE(stdout,"(2X,A,1X,F20.5)") "Reduced chi2: ", &
              rchi2
         WRITE(stdout,"(2X,A,1X,F20.15)") "Lambda: ", lambda
         WRITE(stdout,"(2X,A)") ""
      END IF

    END SUBROUTINE ephemeris_lsl

  END SUBROUTINE levenbergMarquardt_SO





  !! *Description*:
  !!
  !! Arrange observations in pairs. Observation mask used.
  !!
  !! Returns error.
  !!
  !! Modifications needed:
  !! To make ranging more efficient, too close pairs should be avoided.
  !! Also, it might be useful to use only observations from the far ends of
  !! the data set. Where should this filtering be done?
  !!
  SUBROUTINE makePairsOfObservations(this, idx_pair)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout) :: this
    INTEGER, DIMENSION(:,:), POINTER     :: idx_pair
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: idx_pair_
    REAL(bp)                             :: x
    INTEGER                              :: nobs, nobs_orig, n1, n2, ncomb
    INTEGER                              :: i1, i2, nomit, k, err
    INTEGER, DIMENSION(:), ALLOCATABLE   :: ip, ip2

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / makePairsOfObservations", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nomit = 0
    nobs_orig = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / makePairsOfObservations", &
            "TRACE BACK (1", 1)
       RETURN
    END IF

    nobs = COUNT(this%obs_masks_prm(:,2))
    IF (nobs < 2) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / makePairsOfObservations", &
            "Too few observations for pairing.",1)
       RETURN
    END IF
    nomit = nobs_orig - nobs

    ncomb = INT(0.5_bp*nobs*(nobs-1))
    ALLOCATE(idx_pair(ncomb,2),idx_pair_(ncomb,2))

    ALLOCATE(ip(nobs_orig),ip2(nobs))
    ip = (/ (k,k = 1,nobs_orig) /)
    ip2 = PACK(ip,this%obs_masks_prm(:,2))
    x = 1.0_bp*(ncomb-nobs)/nobs
    n1 = 1+INT(x)
    n2 = nobs*(x-INT(x))
    DO k = 1,n1 
       idx_pair(((k-1)*nobs+1):k*nobs,1) = ip2
    END DO
    i1 = 1
    i2 = nobs
    DO k = 1,n1
       idx_pair(i1:i2,2) = CSHIFT(ip2,k)
       i1 = i1+nobs
       i2 = i2+nobs
    END DO

    IF (n2 /= 0) THEN
       idx_pair(i1:i1+n2-1,1) = ip2(1:n2)
       ip = CSHIFT(ip2,n1+1)
       idx_pair(i1:i1+n2-1,2) = CSHIFT(ip(1:n2),n1+1)
    END IF

    ! Order so that first observation number in the pair is always
    ! smaller than the second.
    idx_pair_ = idx_pair
    WHERE(idx_pair(:,2) < idx_pair(:,1))
       idx_pair(:,1) = idx_pair_(:,2)
       idx_pair(:,2) = idx_pair_(:,1)
    END WHERE

    DEALLOCATE(idx_pair_, ip, ip2, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / makePairsOfObservations", &
            "Could not deallocate memory.",1)
       DEALLOCATE(idx_pair, stat=err)
       DEALLOCATE(ip, stat=err)
       DEALLOCATE(ip2, stat=err)       
       RETURN
    END IF

  END SUBROUTINE makePairsOfObservations





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE propagate_SO(this, t, encounters)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout)          :: this
    TYPE (Time), INTENT(inout)                    :: t
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL :: encounters 

    REAL(bp), DIMENSION(:,:,:), POINTER    :: jacobians
    REAL(bp), DIMENSION(6,6)               :: jacobian, cov
    REAL(bp)                               :: det
    INTEGER                                :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / propagate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       !CALL propagate(this%orb_arr_cmp, t)
       IF (PRESENT(encounters)) THEN
          CALL propagate(this%orb_arr_cmp, t, jacobian=jacobians, encounters=encounters)
       ELSE
          CALL propagate(this%orb_arr_cmp, t, jacobian=jacobians)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / propagate", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(jacobians, stat=err)
          RETURN
       END IF
       ! Propagation of pdf
       DO i=1,SIZE(this%orb_arr_cmp)
          det = determinant(jacobians(i,:,:), errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             CALL errorMessage("StochasticOrbit / propagate", &
                  "TRACE BACK (10) " // TRIM(errstr), 1)
             errstr = ""
             DEALLOCATE(jacobians, stat=err)
             RETURN
          END IF
          ! Changed 2008-12-14
          !this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i) / det
          !this%jac_arr_cmp(i,1) = this%jac_arr_cmp(i,1) / det
          this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i) * det
          this%jac_arr_cmp(i,1) = this%jac_arr_cmp(i,1) * det
       END DO
       DEALLOCATE(jacobians,stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / propagate", &
               "Could not deallocate memory.", 1)
          RETURN
       END IF
    END IF
    IF (exist(this%orb_ml_cmp)) THEN
       cov = getCovarianceMatrix(this, & 
            getElementType(this%orb_ml_cmp), & 
            getFrame(this%orb_ml_cmp))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / propagate", &
               "Missing covariances.", 1)
          RETURN
       END IF
       IF (PRESENT(encounters)) THEN
          CALL propagate(this%orb_ml_cmp, t, jacobian=jacobian, encounters=encounters)
       ELSE
          CALL propagate(this%orb_ml_cmp, t, jacobian=jacobian)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / propagate", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       this%cov_ml_cmp = MATMUL(MATMUL(jacobian, cov), &
            TRANSPOSE(jacobian))
    END IF

  END SUBROUTINE propagate_SO





  !! *Description*:
  !!
  !! The default values for the angular deviations ("generation window")
  !! and the maximum residuals ("acceptance window") are adjusted
  !! according to the given factor. These values can be overridden
  !! with the appropriate routine if needed.
  !!
  !! Returns error.
  !!
  SUBROUTINE setAcceptanceWindow_sigma(this, accept_multiplier)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), INTENT(in), OPTIONAL :: accept_multiplier
    REAL(bp), DIMENSION(:,:), POINTER :: stdevs
    REAL(bp) :: accept_multiplier_
    INTEGER :: nobs, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
            "Object has not been initialized.", 1)
       RETURN
    END IF
    IF (PRESENT(accept_multiplier)) THEN
       IF (accept_multiplier < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
               "Acceptance multiplier not accepted.", 1)
          RETURN
       END IF
       accept_multiplier_ = accept_multiplier
       this%accept_multiplier_prm = accept_multiplier
    ELSE IF (this%accept_multiplier_prm > 0) THEN
       accept_multiplier_ = this%accept_multiplier_prm
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
            "Acceptance multiplier not available.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (ASSOCIATED(this%res_accept_prm)) THEN
       IF (SIZE(this%res_accept_prm,dim=1) /= nobs) THEN
          DEALLOCATE(this%res_accept_prm, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
                  "Could not deallocate memory (5).", 1)
             RETURN
          END IF
       END IF
    END IF
    IF (.NOT.ASSOCIATED(this%res_accept_prm)) THEN    
       ALLOCATE(this%res_accept_prm(nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
               "Could not allocate memory.", 1)
          RETURN
       END IF
    END IF
    stdevs => getStandardDeviations(this%obss)
    ! Acceptance windows for residuals:
    this%res_accept_prm(1:nobs,1:6) = accept_multiplier_*stdevs(1:nobs,1:6)
    DEALLOCATE(stdevs, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
            "Could not deallocate memory (10).", 1)
       RETURN
    END IF

  END SUBROUTINE setAcceptanceWindow_sigma








  !! *Description*:
  !!
  !! The default values for the angular deviations ("generation window")
  !! and the maximum residuals ("acceptance window") are adjusted
  !! according to the given factor. These values can be overridden
  !! with the appropriate routine if needed.
  !!
  !! Returns error.
  !!
  SUBROUTINE setGenerationWindow(this, offset)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), DIMENSION(:), INTENT(in), OPTIONAL :: offset
    REAL(bp), DIMENSION(:,:), POINTER :: stdevs
    REAL(bp) :: mean
    INTEGER :: i, j, nobs, err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setGenerationWindow", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (this%sor_generat_multiplier_prm < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setGenerationWindow", &
            "Generation multiplier not acceptable.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (ASSOCIATED(this%sor_deviates_prm)) THEN
       IF (SIZE(this%sor_deviates_prm,dim=1) /= nobs) THEN
          DEALLOCATE(this%sor_deviates_prm, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / setGenerationWindow", &
                  "Could not deallocate memory (5).", 1)
             RETURN
          END IF
       END IF
    END IF
    IF (.NOT.ASSOCIATED(this%sor_deviates_prm)) THEN
       ALLOCATE(this%sor_deviates_prm(nobs,6,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setGenerationWindow", &
               "Could not allocate memory.", 1)
          RETURN
       END IF
    END IF

    this%sor_deviates_prm = 0.0_bp
    ! Angular deviates:
    IF (ASSOCIATED(this%res_arr_cmp)) THEN
       ! Use knowledge of previous results, i.e., move the center of
       ! the sampling window to the mean of the previous residuals:
       IF (SIZE(this%res_arr_cmp,dim=1) /= SIZE(this%orb_arr_cmp,dim=1) .OR. &
            SIZE(this%res_arr_cmp,dim=2) /= nobs) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setGenerationWindow", &
               "Inconsistencies in array dimensions", 1)
          RETURN
       END IF
       DO i=1,nobs
          DO j=1,6
             CALL moments(this%res_arr_cmp(:,i,j), mean=mean, error=errstr)!, std_dev=stdev)!, pdf=this%pdf_arr_cmp)
             IF (LEN_TRIM(errstr) /= 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / setGenerationWindow", &
                     "Could not compute moments. " // TRIM(errstr), 1)
                errstr = ""
                RETURN
             END IF
             this%sor_deviates_prm(i,j,1) = mean
             !this%sor_deviates_prm(i,j,2) = this%sor_generat_multiplier_prm*stdev
          END DO
       END DO
    ELSE
       stdevs => getStandardDeviations(this%obss)
       this%sor_deviates_prm(1:nobs,1:6,1) = 0.0_bp          
       IF (PRESENT(offset)) THEN
          this%sor_deviates_prm(1,2:3,1) = offset(1:2)
          this%sor_deviates_prm(nobs,2:3,1) = offset(3:4)
       END IF
       this%sor_deviates_prm(1:nobs,1:6,2) = this%sor_generat_multiplier_prm*stdevs(1:nobs,1:6)
       DEALLOCATE(stdevs, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setGenerationWindow", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF
    END IF

  END SUBROUTINE setGenerationWindow






  !! *Description*:
  !!
  !! Set parameters.
  !!
  !! Returns error.
  !!
  SUBROUTINE setParameters_SO(this, &
       dyn_model, perturbers, integration_step, integrator, &
       finite_diff, &
       t_inv, element_type, multiple_objects, outlier_rejection, &
       uniform_pdf, regularized_pdf, jacobians_pdf, &
       accept_multiplier, outlier_multiplier, &
       gaussian_pdf, pdf_ml, &
       apriori_a_max, apriori_a_min, apriori_periapsis_max, apriori_periapsis_min, &
       apriori_apoapsis_max, apriori_apoapsis_min, apriori_rho_max, apriori_rho_min, &
       sor_norb, sor_norb_sw, sor_ntrial, sor_ntrial_sw, &
       sor_rho1_l, sor_rho1_u, sor_rho2_l, sor_rho2_u, &
       sor_random_obs_selection, sor_niter, sor_generat_multiplier, &
       sor_generat_offset, sor_2point_method, sor_2point_method_sw, &
       sor_iterate_bounds, &
       vov_norb, vov_ntrial, vov_norb_iter, vov_ntrial_iter, &
       vov_nmap, vov_niter, vov_scaling, vov_mapping_mask, &
       ls_correction_factor, ls_niter_major_max, ls_niter_major_min, ls_niter_minor, &
       ls_element_mask, &
       cos_nsigma, cos_norb, cos_ntrial, cos_gaussian, &
       smplx_tol, smplx_niter)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Time), INTENT(in), OPTIONAL :: t_inv
    CHARACTER(len=*), INTENT(in), OPTIONAL :: &
         dyn_model, &
         integrator, &
         element_type, &
         sor_2point_method, &
         sor_2point_method_sw
    REAL(bp), DIMENSION(6,2), INTENT(in), OPTIONAL :: &
         vov_scaling
    REAL(bp), DIMENSION(6), INTENT(in), OPTIONAL :: &
         finite_diff
    REAL(bp), DIMENSION(4), INTENT(in), OPTIONAL :: &
         sor_generat_offset
    REAL(bp), INTENT(in), OPTIONAL :: &
         integration_step, &
         accept_multiplier, &
         outlier_multiplier, &
         pdf_ml, &
         apriori_a_max, &
         apriori_a_min, &
         apriori_periapsis_max, &
         apriori_periapsis_min, &
         apriori_apoapsis_max, &
         apriori_apoapsis_min, &
         apriori_rho_max, &
         apriori_rho_min, &
         sor_rho1_l, &
         sor_rho1_u, &
         sor_rho2_l, &
         sor_rho2_u, &
         sor_generat_multiplier, &
         ls_correction_factor, &
         cos_nsigma, &
         smplx_tol
    INTEGER, INTENT(in), OPTIONAL :: &
         sor_norb, &
         sor_norb_sw, &
         sor_ntrial, &
         sor_ntrial_sw, &
         sor_niter, &
         vov_norb, &
         vov_ntrial, &
         vov_norb_iter, &
         vov_ntrial_iter, &
         vov_nmap, &
         vov_niter, &
         ls_niter_major_max, &
         ls_niter_major_min, &
         ls_niter_minor, &
         cos_norb, &
         cos_ntrial, &
         smplx_niter
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: &
         perturbers, &
         sor_iterate_bounds, &
         vov_mapping_mask, &
         ls_element_mask
    LOGICAL, INTENT(in), OPTIONAL :: &
         multiple_objects, &
         regularized_pdf, &
         uniform_pdf, &
         jacobians_pdf, &
         sor_random_obs_selection, &
         gaussian_pdf, &
         outlier_rejection, &
         cos_gaussian
    CHARACTER(len=256) :: str
    INTEGER :: i, err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setParameters", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(dyn_model)) THEN
       IF (LEN_TRIM(dyn_model) > DYN_MODEL_LEN) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Parameter 'dyn_model' too long.", 1)
          RETURN
       END IF
       str = dyn_model
       CALL locase(str, error)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / setParameters", &
               "The dynamical model string contains forbidden characters.", 1)
          RETURN
       END IF
       IF (str /= "2-body" .AND. &
            str /= "n-body" .AND. str /= "pseudo-n-body") THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Option " // TRIM(dyn_model) // " not available.", 1)
          RETURN
       END IF
       this%dyn_model_prm = TRIM(str)
       IF (str == "2-body" .OR. str == "n-body") THEN
          IF (ASSOCIATED(this%orb_arr_cmp)) THEN
             DO i=1,SIZE(this%orb_arr_cmp)
                CALL setParameters(this%orb_arr_cmp(i), dyn_model=this%dyn_model_prm)
             END DO
          END IF
          IF (exist(this%orb_ml_cmp)) THEN
             CALL setParameters(this%orb_ml_cmp, dyn_model=this%dyn_model_prm)
          END IF
       END IF
    END IF
    IF (PRESENT(perturbers)) THEN
       IF (SIZE(perturbers) > SIZE(this%perturbers_prm)) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Perturber array too small.", 1)
          RETURN                       
       END IF
       this%perturbers_prm = perturbers
       IF (ASSOCIATED(this%orb_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL setParameters(this%orb_arr_cmp(i), perturbers=this%perturbers_prm)
          END DO
       END IF
       IF (exist(this%orb_ml_cmp)) THEN
          CALL setParameters(this%orb_ml_cmp, perturbers=this%perturbers_prm)
       END IF
    END IF
    IF (PRESENT(integration_step)) THEN
       this%integration_step_prm = integration_step
       IF (ASSOCIATED(this%orb_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL setParameters(this%orb_arr_cmp(i), integration_step=this%integration_step_prm)
          END DO
       END IF
       IF (exist(this%orb_ml_cmp)) THEN
          CALL setParameters(this%orb_ml_cmp, integration_step=this%integration_step_prm)
       END IF
    END IF
    IF (PRESENT(integrator)) THEN
       str = integrator
       CALL locase(str, error)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / setParameters", &
               "The integrator string contains forbidden characters.", 1)
          RETURN
       END IF
       this%integrator_prm = TRIM(str)
       IF (ASSOCIATED(this%orb_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL setParameters(this%orb_arr_cmp(i), integrator=this%integrator_prm)
          END DO
       END IF
       IF (exist(this%orb_ml_cmp)) THEN
          CALL setParameters(this%orb_ml_cmp, integrator=this%integrator_prm)
       END IF
    END IF
    IF (PRESENT(finite_diff)) THEN
       IF (ALL(finite_diff > 0.0_bp)) THEN
          IF (.NOT.ASSOCIATED(this%finite_diff_prm)) THEN
             ALLOCATE(this%finite_diff_prm(6), stat=err)
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / setParameters", &
                     "Could not allocate memory.", 1)
                RETURN
             END IF
          END IF
          this%finite_diff_prm = finite_diff
          IF (ASSOCIATED(this%orb_arr_cmp)) THEN
             DO i=1,SIZE(this%orb_arr_cmp)
                CALL setParameters(this%orb_arr_cmp(i), finite_diff=this%finite_diff_prm)
             END DO
          END IF
          IF (exist(this%orb_ml_cmp)) THEN
             CALL setParameters(this%orb_ml_cmp, finite_diff=this%finite_diff_prm)
          END IF
       ENDIF
    END IF
    IF (PRESENT(t_inv)) THEN
       this%t_inv_prm = copy(t_inv)
    END IF
    IF (PRESENT(element_type)) THEN
       this%element_type_prm = element_type
       DO i=1,LEN(this%element_type_prm)
          IF (IACHAR(this%element_type_prm(i:i)) > 127) THEN
             error = .TRUE.
             EXIT
          END IF
       END DO
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / setParameters", &
               "The element type string contains forbidden characters.", 1)
          RETURN
       END IF
       CALL locase(this%element_type_prm, error)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / setParameters", &
               "The element type string contains forbidden characters.", 1)
          RETURN
       END IF
       IF (this%element_type_prm /= "keplerian" .AND. &
            this%element_type_prm /= "cartesian" .AND. &
            this%element_type_prm /= "cometary") THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Orbital element type is not recognized: " // &
               TRIM(element_type) // ".", 1)
          RETURN          
       END IF
    END IF
    IF (PRESENT(multiple_objects)) THEN
       this%multiple_obj_prm = multiple_objects
    END IF
    IF (PRESENT(outlier_rejection)) THEN
       this%outlier_rejection_prm = outlier_rejection
    END IF
    IF (PRESENT(outlier_multiplier)) THEN
       IF (outlier_multiplier <= 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Multiplier for the outlier criterion is zero or negative.", 1)
          RETURN
       END IF
       this%outlier_multiplier_prm = outlier_multiplier
    END IF
    IF (PRESENT(uniform_pdf)) THEN
       this%uniform_pdf_prm = uniform_pdf
    END IF
    IF (PRESENT(regularized_pdf)) THEN
       this%regularization_prm = regularized_pdf
    END IF
    IF (PRESENT(jacobians_pdf)) THEN
       this%jacobians_prm = jacobians_pdf
    END IF
    IF (PRESENT(accept_multiplier)) THEN
       IF (accept_multiplier <= 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Scaling of acceptance window is zero or negative.", 1)
          RETURN
       END IF
       this%accept_multiplier_prm = accept_multiplier
       CALL setAcceptanceWindow(this)
    END IF
    IF (PRESENT(pdf_ml)) THEN
       this%pdf_ml_init_prm = pdf_ml
    END IF
    IF (PRESENT(apriori_a_max)) THEN
       this%apriori_a_max_prm = apriori_a_max
       this%informative_apriori_prm = .TRUE.
    END IF
    IF (PRESENT(apriori_a_min)) THEN
       this%apriori_a_min_prm = apriori_a_min
       this%informative_apriori_prm = .TRUE.
    END IF
    IF (PRESENT(apriori_periapsis_max)) THEN
       this%apriori_periapsis_max_prm = apriori_periapsis_max
       this%informative_apriori_prm = .TRUE.
    END IF
    IF (PRESENT(apriori_periapsis_min)) THEN
       this%apriori_periapsis_min_prm = apriori_periapsis_min
       this%informative_apriori_prm = .TRUE.
    END IF
    IF (PRESENT(apriori_apoapsis_max)) THEN
       this%apriori_apoapsis_max_prm = apriori_apoapsis_max
       this%informative_apriori_prm = .TRUE.
    END IF
    IF (PRESENT(apriori_apoapsis_min)) THEN
       this%apriori_apoapsis_min_prm = apriori_apoapsis_min
       this%informative_apriori_prm = .TRUE.
    END IF
    IF (PRESENT(apriori_rho_max)) THEN
       this%apriori_rho_max_prm = apriori_rho_max
    END IF
    IF (PRESENT(apriori_rho_min)) THEN
       this%apriori_rho_min_prm = apriori_rho_min
    END IF
    IF (PRESENT(sor_norb)) THEN
       IF (sor_norb < 0 .OR. sor_norb > 100000) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Number of sample orbits must be positive and " // &
               "less than 100001.", 1)
          RETURN
       ELSE
          this%sor_norb_prm = sor_norb
       END IF
    END IF
    IF (PRESENT(sor_norb_sw)) THEN
       IF (sor_norb_sw < 0 .OR. sor_norb_sw > 10000) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Number of sample orbits for the stepwise " // &
               "approach must be positive and " // &
               "less than 10001.", 1)
          RETURN
       ELSE
          this%sor_norb_sw_prm = sor_norb_sw
       END IF
    END IF
    IF (PRESENT(sor_ntrial)) THEN
       IF (sor_ntrial < 0 .OR. sor_ntrial > 100000000) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Number of trial orbits must be positive and " // &
               "less than 100,000,001.", 1)
          RETURN
       ELSE
          this%sor_ntrial_prm = sor_ntrial
       END IF
    END IF
    IF (PRESENT(sor_ntrial_sw)) THEN
       IF (sor_ntrial_sw < 0 .OR. sor_ntrial_sw > 1000000) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Number of trial orbits for the stepwise " // & 
               "must be positive and less than 1,000,001.", 1)
          RETURN
       ELSE
          this%sor_ntrial_sw_prm = sor_ntrial_sw
       END IF
    END IF
    IF (PRESENT(sor_rho1_l)) THEN
       this%sor_rho_prm(1,1) = sor_rho1_l
    END IF
    IF (PRESENT(sor_rho1_u)) THEN
       this%sor_rho_prm(1,2) = sor_rho1_u
    END IF
    IF (PRESENT(sor_rho2_l)) THEN
       this%sor_rho_prm(2,1) = sor_rho2_l
    END IF
    IF (PRESENT(sor_rho2_u)) THEN
       this%sor_rho_prm(2,2) = sor_rho2_u
    END IF
    IF (PRESENT(sor_random_obs_selection)) THEN
       IF (sor_random_obs_selection) THEN
          CALL setRandomObservationSelection(this)
       END IF
    END IF
    IF (PRESENT(sor_niter)) THEN
       this%sor_niter_prm = sor_niter
    END IF
    IF (PRESENT(gaussian_pdf)) THEN
       IF (gaussian_pdf) THEN
          this%sor_gaussian_pdf_prm = .TRUE.
       ELSE
          this%sor_gaussian_pdf_prm = .FALSE.
       END IF
       !this%sor_gaussian_pdf_prm = gaussian_pdf
    END IF
    IF (PRESENT(sor_generat_multiplier) .AND. PRESENT(sor_generat_offset)) THEN
       IF (sor_generat_multiplier <= 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Scaling of generation window is zero or negative.", 1)
          RETURN
       END IF
       this%sor_generat_multiplier_prm = sor_generat_multiplier
       CALL setGenerationWindow(this, sor_generat_offset)
    ELSE IF (PRESENT(sor_generat_multiplier)) THEN
       IF (sor_generat_multiplier <= 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Scaling of generation window is zero or negative.", 1)
          RETURN
       END IF
       this%sor_generat_multiplier_prm = sor_generat_multiplier
       CALL setGenerationWindow(this)
    END IF
    IF (PRESENT(sor_2point_method)) THEN
       this%sor_2point_method_prm = sor_2point_method
    END IF
    IF (PRESENT(sor_2point_method_sw)) THEN
       this%sor_2point_method_sw_prm = sor_2point_method_sw
    END IF
    IF (PRESENT(sor_iterate_bounds)) THEN
       this%sor_iterate_bounds_prm = sor_iterate_bounds
    END IF
    IF (PRESENT(vov_norb)) THEN
       this%vov_norb_prm = vov_norb
    END IF
    IF (PRESENT(vov_ntrial)) THEN
       this%vov_ntrial_prm = vov_ntrial
    END IF
    IF (PRESENT(vov_norb_iter)) THEN
       this%vov_norb_iter_prm = vov_norb_iter
    END IF
    IF (PRESENT(vov_ntrial_iter)) THEN
       this%vov_ntrial_iter_prm = vov_ntrial_iter
    END IF
    IF (PRESENT(vov_nmap)) THEN
       this%vov_nmap_prm = vov_nmap
    END IF
    IF (PRESENT(vov_niter)) THEN
       this%vov_niter_prm = vov_niter
    END IF
    IF (PRESENT(vov_scaling)) THEN
       this%vov_scaling_prm = vov_scaling
    END IF
    IF (PRESENT(vov_mapping_mask)) THEN
       this%vov_mapping_mask_prm = vov_mapping_mask
    END IF
    IF (PRESENT(ls_correction_factor)) THEN
       this%ls_corr_fac_prm = ls_correction_factor
    END IF
    IF (PRESENT(ls_niter_major_max)) THEN
       this%ls_niter_major_max_prm = ls_niter_major_max
    END IF
    IF (PRESENT(ls_niter_major_min)) THEN
       this%ls_niter_major_min_prm = ls_niter_major_min
    END IF
    IF (PRESENT(ls_niter_minor)) THEN
       this%ls_niter_minor_prm = ls_niter_minor
    END IF
    IF (PRESENT(ls_element_mask)) THEN
       this%ls_elem_mask_prm = ls_element_mask
    END IF
    IF (PRESENT(cos_nsigma)) THEN
       this%cos_nsigma_prm = cos_nsigma
    END IF
    IF (PRESENT(cos_norb)) THEN
       this%cos_norb_prm = cos_norb
    END IF
    IF (PRESENT(cos_ntrial)) THEN
       this%cos_ntrial_prm = cos_ntrial
    END IF
    IF (PRESENT(cos_gaussian)) THEN
       this%cos_gaussian_prm = cos_gaussian
    END IF
    IF (PRESENT(smplx_tol)) THEN
       this%smplx_tol_prm = smplx_tol
    END IF
    IF (PRESENT(smplx_niter)) THEN
       this%smplx_niter_prm = smplx_niter
    END IF

  END SUBROUTINE setParameters_SO





  !! *Description*:
  !!
  !! Sets the ranging parameters that are suitable for NEOs.
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE setNEORanging(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setNEORanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    CALL setRangeBounds(this, 0.1_bp, 6.0_bp, -1.0_bp, 1.0_bp)
    CALL setParameters(this, sor_norb=2000, sor_ntrial=1000000)

  END SUBROUTINE setNEORanging



  !! *Description*:
  !!
  !! Updates the observation mask of the observation specified, i.e.,
  !! the inclusion or omission of observations in orbit inversion is
  !! changed.
  !!
  !! Returns error.
  !!
  SUBROUTINE setObservationMask_one(this, i, obs_mask)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    LOGICAL, DIMENSION(6), INTENT(in)     :: obs_mask
    INTEGER, INTENT(in)                   :: i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationMask", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (i > getNrOfObservations(this%obss)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationMask", &
            "Unrealistic input values.", 1)
       RETURN
    END IF

    this%obs_masks_prm(i,:) = obs_mask

  END SUBROUTINE setObservationMask_one





  !! *Description*:
  !!
  !! Updates the observation mask of the observations, i.e.,
  !! the inclusion or omission of observations in orbit inversion is
  !! changed.
  !!
  !! Returns error.
  !!
  SUBROUTINE setObservationMask_all(this, obs_masks)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    LOGICAL, DIMENSION(:,:), INTENT(in)   :: obs_masks

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationMask", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    this%obs_masks_prm = obs_masks

  END SUBROUTINE setObservationMask_all





  !! *Description*:
  !!
  !! Updates the observation mask of the observation specified using
  !! user-specified information from observation notes, i.e.,
  !! the inclusion or omission of observations in orbit inversion is
  !! changed.
  !!
  !! Returns error.
  !!
  SUBROUTINE setObservationMask_all_notes(this, use_notes)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    LOGICAL, INTENT(in)                   :: use_notes
    INTEGER :: err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationMask", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF(.NOT. exist(this%obss)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationMask", &
            "Observations not intialized.", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%obs_masks_prm)) THEN
       DEALLOCATE(this%obs_masks_prm, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setObservationMask", &
               "Could not deallocate memory.", 1)
          RETURN          
       END IF
    END IF
    this%obs_masks_prm => getObservationMasks(this%obss, use_notes)

  END SUBROUTINE setObservationMask_all_notes




  !! Description:
  !!
  !! Sets observation pair to default values:
  !! 1) if none excluded, first and last observation of the set
  !! 2) if some excluded, first and last observation NOT excluded
  !!
  !! Returns error.
  !!
  SUBROUTINE setObservationPair_default(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    INTEGER                               :: nobs, i, err, iobs1, iobs2 
    INTEGER, DIMENSION(:), ALLOCATABLE    :: intarr

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (ASSOCIATED(this%sor_pair_arr_prm)) THEN
       DEALLOCATE(this%sor_pair_arr_prm, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setObservationPair", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(intarr(nobs), this%sor_pair_arr_prm(1,2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    DO i=1,nobs
       intarr(i) = i
    END DO
    iobs1 = MINLOC(intarr, dim=1,mask=(this%obs_masks_prm(:,2) .AND. &
         this%obs_masks_prm(:,3)))
    iobs2 = MAXLOC(intarr, dim=1,mask=(this%obs_masks_prm(:,2) .AND. &
         this%obs_masks_prm(:,3)))
    this%sor_pair_arr_prm(1,1) = iobs1
    this%sor_pair_arr_prm(1,2) = iobs2
    IF (this%sor_pair_arr_prm(1,1) == this%sor_pair_arr_prm(1,2)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Only one observation available.", 1)
       RETURN
    END IF
    DEALLOCATE(intarr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Could not deallocate memory (10).", 1)
       RETURN
    END IF

  END SUBROUTINE setObservationPair_default





  !! *Description*:
  !!
  !! Sets observation pair as given.
  !!
  !! Returns error.
  !!
  SUBROUTINE setObservationPair_pair(this, iobs1, iobs2)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    INTEGER, INTENT(in)                   :: iobs1, iobs2
    INTEGER                               :: nobs, err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (iobs1 <= 0 .OR. iobs2 <= 0 .OR. iobs1 > nobs .OR. iobs2 > nobs &
         .OR. iobs1 == iobs2) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Unrealistic input values.", 1)
       RETURN
    END IF
    IF (.NOT. this%obs_masks_prm(iobs1,2) .OR. .NOT. this%obs_masks_prm(iobs2,2)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "One or both observations excluded.", 1)
       RETURN
    END IF
    IF (ASSOCIATED(this%sor_pair_arr_prm)) THEN
       DEALLOCATE(this%sor_pair_arr_prm, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setObservationPair", &
               "Could not deallocate memory.", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(this%sor_pair_arr_prm(1,2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setObservationPair", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    this%sor_pair_arr_prm(1,1) = iobs1
    this%sor_pair_arr_prm(1,2) = iobs2

  END SUBROUTINE setObservationPair_pair





  !! *Description*:
  !!
  !! Sets observation pair selection to random selection.
  !! Constructs a list of possible pairs
  !! (number of pairs = number of combinations of two).
  !!
  !! Returns error.
  !!
  SUBROUTINE setRandomObservationSelection(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    INTEGER, DIMENSION(:,:), POINTER      :: idx_pair
    INTEGER                               :: err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRandomObservationSelection", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    this%sor_random_obs_prm = .TRUE.
    CALL makePairsOfObservations(this,idx_pair)
    DEALLOCATE(this%sor_pair_arr_prm, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRandomObservationSelection", &
            "Could not deallocate memory (5).", 1)
       RETURN
    END IF
    ALLOCATE(this%sor_pair_arr_prm(SIZE(idx_pair,dim=1), SIZE(idx_pair,dim=2)), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRandomObservationSelection", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    this%sor_pair_arr_prm = idx_pair
    DEALLOCATE(idx_pair, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRandomObservationSelection", &
            "Could not deallocate memory (10).", 1)
       RETURN
    END IF

  END SUBROUTINE setRandomObservationSelection





  !! *Description*:
  !!
  !! Determines the new range bounds from the 3-sigma cutoff values
  !! of the a posteriori range probability density.
  !! 
  !!
  !! Returns error.
  !!
  SUBROUTINE setRangeBounds_3sigma(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), DIMENSION(2)                :: mean, stdev
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: histo
    REAL(bp), DIMENSION(:), ALLOCATABLE   :: topo_range
    REAL(bp), DIMENSION(2,2)              :: rho_
    REAL(bp)                              :: rhomin1, rhomax1, &
         rhomin2, rhomax2
    INTEGER                               :: grid_num, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRangeBounds", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    rho_ = this%sor_rho_prm
    this%sor_rho_histo_cmp = 0
    grid_num   = 40

    IF (ANY(this%sor_iterate_bounds_prm(1:2))) THEN

       ! Compute mean and std for first topocentric range,
       ! create histogram:
       CALL moments(this%sor_rho_arr_cmp(:,1), mean=mean(1), &
            std_dev=stdev(1), error=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "Could not compute moments (5). " // TRIM(errstr), 1)
          errstr = ""
          RETURN
       END IF
       rhomin1 = MINVAL(this%sor_rho_arr_cmp(:,1))
       rhomax1 = MAXVAL(this%sor_rho_arr_cmp(:,1))
       ALLOCATE(histo(grid_num,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "Could not allocate memory.", 1)
          RETURN
       END IF
       CALL histogram(this%sor_rho_arr_cmp(:,1), histo)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "TRACE BACK", 1)
          RETURN
       END IF
       histo(:,2) = histo(:,2)/MAXVAL(histo(:,2),dim=1)
       ! Check the end points of the histogram, raise flag if the wings do not
       ! fall off as expected. Typically, histo_end = 0.3_bp
       IF (histo(1,2) > histo_end) THEN ! .OR. histo(grid_num,2) > histo_end) THEN
          this%sor_rho_histo_cmp = 1
       END IF
       IF (histo(grid_num,2) > histo_end) THEN
          this%sor_rho_histo_cmp = this%sor_rho_histo_cmp + 2
       END IF
       DEALLOCATE(histo, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
       ! 3-sigma bounds for topocentric ranges:
       IF (this%sor_iterate_bounds_prm(1)) THEN
          rho_(1,1) = mean(1) - 3*stdev(1)
       END IF
       IF (this%sor_iterate_bounds_prm(2)) THEN
          rho_(1,2) = mean(1) + 3*stdev(1)
       END IF
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,2(1X,F10.5))") &
               "rho1 mean and standard deviation: ", &
               mean(1), stdev(1)
       END IF

    END IF

    IF (ANY(this%sor_iterate_bounds_prm(3:4))) THEN

       ! Second topocentric range (with respect to the first)
       ALLOCATE(topo_range(this%sor_norb_cmp), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "Could not allocate memory for arrays.", 1)
          RETURN
       END IF
       topo_range = this%sor_rho_arr_cmp(:,2) - this%sor_rho_arr_cmp(:,1)
       rhomin2 = MINVAL(topo_range)
       rhomax2 = MAXVAL(topo_range)
       CALL moments(topo_range, mean=mean(2), std_dev=stdev(2), &
            error=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "Could not compute moments (10). " // TRIM(errstr), 1)
          errstr = ""
          RETURN
       END IF
       DEALLOCATE(topo_range, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setRangeBounds", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF
       ! 3-sigma bounds for topocentric ranges:
       IF (this%sor_iterate_bounds_prm(3)) THEN
          rho_(2,1) = mean(2) - 3*stdev(2)
       END IF
       IF (this%sor_iterate_bounds_prm(4)) THEN
          rho_(2,2) = mean(2) + 3*stdev(2)
       END IF
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,2(1X,F10.5))") &
               "rho2-rho1 mean and standard deviation: ", &
               mean(2), stdev(2)
       END IF

    END IF

    IF (rho_(1,1) < MAX(0.0_bp,this%apriori_rho_min_prm)) THEN
       rho_(1,1) = MAX(0.0_bp,this%apriori_rho_min_prm)
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,4(1X,F10.5))") &
            "Updated rho bounds: ", rho_(1,:), rho_(2,:)
    END IF
    this%sor_rho_prm = rho_

  END SUBROUTINE setRangeBounds_3sigma





  !! *Description*:
  !!
  !! Sets new range bounds given the lower and upper values for rho1 and rho2.
  !!
  !! Returns error.
  !!
  SUBROUTINE setRangeBounds_values(this, rho1_lower, rho1_upper, &
       rho2_lower, rho2_upper)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), INTENT(in)                             :: rho1_lower, rho1_upper, &
         rho2_lower, rho2_upper

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRangeBounds", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    this%sor_rho_prm(1,1) = rho1_lower
    this%sor_rho_prm(1,2) = rho1_upper
    this%sor_rho_prm(2,1) = rho2_lower
    this%sor_rho_prm(2,2) = rho2_upper

  END SUBROUTINE setRangeBounds_values



  !! *Description*:
  !!
  !! Sets the regularization of the statistical treatment:
  !! ON (.true.), OFF (.false.). By default, should always
  !! be ON.
  !!
  !!
  !! Returns error.
  !!

  SUBROUTINE setRegularization(this,regularization)
    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    LOGICAL, INTENT(in)                   :: regularization

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setRegularizationOff", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    this%regularization_prm = regularization

  END SUBROUTINE setRegularization



  !! *Description*:
  !!
  !! Sets the ranging parameters are suitable for TNOs.
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE setTNORanging(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setTNORanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    CALL setRangeBounds(this, 20.0_bp, 60.0_bp, -1.0_bp, 1.0_bp)
    CALL setParameters(this, sor_norb=2000, sor_ntrial=1000000)
    CALL setAcceptanceWindow(this, 4.0_bp)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / setTNORanging", &
            "TRACE BACK (1", 1)
       RETURN
    END IF

  END SUBROUTINE setTNORanging





  !! *Description*:
  !!
  !! Minimization of the function func in N dimensions by the
  !! downhill simplex method of Nelder and Mead. The (N + 1) × N
  !! matrix p is input. Its N + 1 rows are N-dimensional vectors that
  !! are the vertices of the starting simplex. Also input is the
  !! vector y of length N + 1, whose components must be preinitialized
  !! to the values of func evaluated at the N + 1 vertices (rows) of
  !! p and ftol the fractional convergence tolerance to be achieved in
  !! the function value (n.b.!). The input value of iter is the
  !! maximum number of function evaluations, while on output the same
  !! parameter gives the actual number of function evaluations
  !! perfomed. Also on output, p and y will have been reset to N+1 new
  !! points all within this%smplx_tol_prm of a minimum function value.
  !!
  !! Parameters: The maximum allowed number of function evaluations,
  !! and a small number.
  !!
  SUBROUTINE simplexOrbits(this, orb_arr)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), DIMENSION(:), INTENT(inout) :: orb_arr

    REAL(bp), PARAMETER :: TINY = 1.0e-10_bp

    CHARACTER(len=FRAME_LEN) :: frame
    REAL(bp), DIMENSION(7,6) :: p, p_init
    REAL(bp), DIMENSION(7) :: y
    REAL(bp), DIMENSION(6) :: psum, p_best
    REAL(bp) :: y_best, rchi2_denom
    INTEGER(ibp) :: ihi, ilo, ndim ! Global variables (within this subroutine).  
    INTEGER :: err, i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / simplexOrbits", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    CALL NULLIFY(this%orb_ml_cmp)
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DO i=1,SIZE(this%orb_arr_cmp)
          CALL NULLIFY(this%orb_arr_cmp(i))
       END DO
       DEALLOCATE(this%orb_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp, stat=err)
    END IF

    ! Denominator for reduced chi2:
    rchi2_denom = 1.0_bp/REAL(COUNT(this%obs_masks_prm)-6,bp)

    frame = getFrame(orb_arr(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / simplexOrbits", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    CALL propagate(orb_arr, this%t_inv_prm)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / simplexOrbits", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    DO i=1,7
       p(i,1:6) = getElements(orb_arr(i), this%element_type_prm, &
            frame=frame)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / simplexOrbits", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       p_init(i,1:6) = p(i,1:6)
       y(i) = getChi2(this, orb_arr(i))*rchi2_denom
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / simplexOrbits", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
    END DO
    p_best = p(MINLOC(y,dim=1),:)
    y_best = MINVAL(y)
    CALL simplex_private
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / simplexOrbits", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       DO i=1,7
          IF (i == ilo) THEN
             WRITE(stdout,"('ILO:',1X,6(F20.15,1X),1(F25.5,1X))") &
                  p(i,:), y(i)
          ELSE IF (i == ihi) THEN
             WRITE(stdout,"('IHI:',1X,6(F20.15,1X),1(F25.5,1X))") &
                  p(i,:), y(i)
          ELSE
             WRITE(stdout,"(5X,6(F20.15,1X),1(F25.5,1X))") &
                  p(i,:), y(i)
          END IF
       END DO
       WRITE(stdout,"(1X)")
    END IF
    IF (this%smplx_niter_cmp >= this%smplx_niter_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / simplexOrbits", &
            "Maximum number of iterations exceeded.", 1)
       RETURN
    END IF
    ALLOCATE(this%orb_arr_cmp(7), this%rchi2_arr_cmp(7), stat=err)
    DO i=1,7
       CALL NEW(this%orb_arr_cmp(i), p(i,:), this%element_type_prm, &
            frame, this%t_inv_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / simplexOrbits", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       CALL setParameters(this%orb_arr_cmp(i), &
            dyn_model=this%dyn_model_prm, &
            perturbers=this%perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / simplexOrbits", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       this%rchi2_arr_cmp(i) = y(i)
    END DO
    this%orb_ml_cmp = copy(this%orb_arr_cmp(1))

  CONTAINS 

    SUBROUTINE simplex_private 

      IMPLICIT NONE
      TYPE (Orbit) :: orb
      REAL(bp), DIMENSION(6) :: ptry_ref, ptry_exp, ptry_co
      REAL(bp) :: rtol, ytry_ref, ytry_exp, ytry_co, ytmp
      INTEGER(ibp) :: i, j, inhi
      LOGICAL :: contraction_successful

      ndim = SIZE(p,dim=2)
      this%smplx_niter_cmp = 0 
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
         IF (y_best > y(ilo)) THEN
            p_best = p(ilo,:)
            y_best = y(ilo)
         END IF
         ! Compute the fractional range from highest to lowest and
         ! reinitialize if very similar values
         IF (2.0_bp*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo))+TINY) &
              < 0.000001_bp) THEN
            j = 0
            DO i=1,ndim+1
               IF (i == ilo) THEN
                  CYCLE
               END IF
               j = j + 1
               CALL NULLIFY(orb_arr(i))
               p(i,:) = p(ilo,:)
               p(i,j) = p(i,j)*1.1_bp
               CALL NEW(orb_arr(i), p(i,:), this%element_type_prm, &
                    frame, this%t_inv_prm)
               CALL setParameters(orb_arr(i), &
                    dyn_model=this%dyn_model_prm, &
                    perturbers=this%perturbers_prm, &
                    integrator=this%integrator_prm, &
                    integration_step=this%integration_step_prm)
               y(i) = getChi2(this, orb_arr(i))*rchi2_denom
            END DO
            IF (error) THEN
               CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                    "TRACE BACK (20)", 1)
               RETURN
            END IF
            ilo = iminloc(y(:))
            ihi = imaxloc(y(:))
            ytmp = y(ihi) 
            y(ihi) = y(ilo)
            inhi = imaxloc(y(:))
            y(ihi) = ytmp
            IF (info_verb >= 2) THEN
               WRITE(stdout,"(5X,A)") "Reinitialization performed."
            END IF
         END IF
         IF (info_verb >= 3) THEN
            DO i=1,7
               IF (i == ilo) THEN
                  WRITE(stdout,"('ILO:',1X,6(F20.15,1X),1(F15.3,1X))") &
                       p(i,:), y(i)
               ELSE IF (i == ihi) THEN
                  WRITE(stdout,"('IHI:',1X,6(F20.15,1X),1(F15.3,1X))") &
                       p(i,:), y(i)
               ELSE
                  WRITE(stdout,"(5X,6(F20.15,1X),1(F15.3,1X))") &
                       p(i,:), y(i)
               END IF
            END DO
            WRITE(stdout,"(1X)")
         END IF
         rtol = y(ihi)
         ! return if satisfactory.  
         IF (rtol < this%smplx_tol_prm) THEN 
            ! If returning, put best point and value in slot 1. 
            CALL swap(y(1),y(ilo))
            CALL swap(p(1,:),p(ilo,:)) 
            ilo = 1
            RETURN 
         END IF
         IF (this%smplx_niter_cmp >= this%smplx_niter_prm) THEN
            RETURN
         END IF
         ! Begin a new iteration. First extrapolate by a factor -1
         ! through the face of the simplex across from the high
         ! point, i.e., reflect the simplex from the high point.
         CALL simtry(this, 1.0_bp, ptry_ref, ytry_ref) ! Reflect
         IF (error) THEN
            CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                 "TRACE BACK (30)", 1)
            RETURN
         END IF
         this%smplx_niter_cmp = this%smplx_niter_cmp + 1 
         IF (ytry_ref >= y(ilo) .AND. ytry_ref < y(inhi)) THEN
            ! If it's better than the highest, then replace
            ! the highest.
            y(ihi) = ytry_ref
            psum(:) = psum(:) - p(ihi,:) + ptry_ref(:)
            p(ihi,:) = ptry_ref(:)
         ELSE IF (ytry_ref < y(ilo)) THEN
            ! Gives a result better than the best point, so try an
            ! additional extrapolation by a factor of 2.  
            CALL simtry(this, 2.0_bp, ptry_exp, ytry_exp) ! Expand
            IF (error) THEN
               CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                    "TRACE BACK (35)", 1)
               RETURN
            END IF
            this%smplx_niter_cmp = this%smplx_niter_cmp + 1
            IF (ytry_exp < ytry_ref) THEN
               y(ihi) = ytry_exp
               psum(:) = psum(:) - p(ihi,:) + ptry_exp(:)
               p(ihi,:) = ptry_exp(:)
            ELSE
               y(ihi) = ytry_ref
               psum(:) = psum(:) - p(ihi,:) + ptry_ref(:)
               p(ihi,:) = ptry_ref(:)
            END IF
         ELSE IF (ytry_ref >= y(inhi)) THEN 
            ! The reflected point is worse than the second highest,
            ! so look for an intermediate lower point, i.e., do a
            ! one-dimensional contraction.
            contraction_successful = .FALSE.
            IF (ytry_ref >= y(inhi) .AND. ytry_ref < y(ihi)) THEN
               CALL simtry(this, 0.5_bp, ptry_co, ytry_co) ! Outside contraction
               IF (error) THEN
                  CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                       "TRACE BACK (40)", 1)
                  RETURN
               END IF
               IF (ytry_co < ytry_ref) THEN
                  y(ihi) = ytry_co
                  psum(:) = psum(:) - p(ihi,:) + ptry_co(:)
                  p(ihi,:) = ptry_co(:)
                  contraction_successful = .TRUE. 
               END IF
            ELSE ! ytry_ref >= y(ihi) 
               CALL simtry(this, -0.5_bp, ptry_co, ytry_co) ! Inside contraction
               IF (error) THEN
                  CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                       "TRACE BACK (40)", 1)
                  RETURN
               END IF
               IF (ytry_co < y(ihi)) THEN
                  y(ihi) = ytry_co
                  psum(:) = psum(:) - p(ihi,:) + ptry_co(:)
                  p(ihi,:) = ptry_co(:)
                  contraction_successful = .TRUE. 
               END IF
            END IF
            this%smplx_niter_cmp = this%smplx_niter_cmp + 1
            IF (.NOT.contraction_successful) THEN
               ! Can't seem to get rid of that high point. Better
               ! contract around the lowest (best) point.
               p(:,:) = 0.5_bp*(p(:,:) + SPREAD(p(ilo,:),1,SIZE(p,dim=1)))
               DO i=1,ndim+1
                  IF (i /= ilo) THEN
                     IF (this%element_type_prm == "keplerian") THEN
                        IF (p(i,1) < this%apriori_a_min_prm .OR. &
                             p(i,1) > this%apriori_a_max_prm .OR. &
                             p(i,2) < 0.0_bp .OR. &
                             p(i,2) > 1.0_bp .OR. &
                             p(i,3) < 0.0_bp .OR. &
                             p(i,3) > pi) THEN
                           p(i,:) = p_init(i,:)
                        END IF
                        p(i,4:6) = MODULO(p(i,4:6), two_pi)
                     ELSE IF (this%element_type_prm == "cartesian") THEN
!!!$                        IF (SQRT(SUM(p(i,1:3)**2)) > this%apriori_hcentric_dist_max_prm .or. &
!!!$                             SQRT(SUM(p(i,4:6)**2)) > this%apriori_velocity_max_prm) THEN
!!!$                           p(i,:) = p_init(i,:)
!!!$                        END IF
                     END IF
                     CALL NEW(orb, p(i,:), this%element_type_prm, &
                          frame, this%t_inv_prm)
                     IF (error) THEN
                        CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     CALL setParameters(orb, &
                          dyn_model=this%dyn_model_prm, &
                          perturbers=this%perturbers_prm, &
                          integrator=this%integrator_prm, &
                          integration_step=this%integration_step_prm)
                     IF (error) THEN
                        CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     y(i) = getChi2(this, orb)*rchi2_denom
                     IF (error) THEN
                        CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                             "TRACE BACK (50)", 1)
                        RETURN
                     END IF
                     CALL NULLIFY(orb)
                  END IF
               END DO
               this%smplx_niter_cmp = this%smplx_niter_cmp + ndim ! Keep track of function evaluations.
               psum(:) = SUM(p(:,:),dim=1)
            END IF
         END IF
      END DO ! Go back for the test of doneness and the next iteration.

    END SUBROUTINE simplex_private


    !! *Description*:
    !!
    !! Extrapolates by a factor fac through the face of the
    !! simplex across from the high point, tries it, and replaces
    !! the high point if the new point is better.
    !!
    SUBROUTINE simtry(this, fac, ptry, ytry)

      IMPLICIT NONE
      TYPE (StochasticOrbit), INTENT(in) :: this
      REAL(bp), INTENT(IN) :: fac
      REAL(bp), DIMENSION(ndim), INTENT(out) :: ptry 
      REAL(bp), INTENT(out) :: ytry

      TYPE (Orbit) :: orb
      REAL(bp) :: ran

      ptry(:) = (1.0_bp+fac)*(psum(:)-p(ihi,:))/ndim - fac*p(ihi,:)
      ! Evaluate the function at the trial point.
      IF (this%element_type_prm == "keplerian") THEN
         IF (ptry(1) < this%apriori_a_min_prm .OR. &
              ptry(1) > this%apriori_a_max_prm .OR. &
              ptry(2) < 0.0_bp .OR. &
              ptry(2) > 1.0_bp .OR. &
              ptry(3) < 0.0_bp .OR. &
              ptry(3) > pi) THEN
            CALL randomNumber(ran)
            ptry = p_init(1+NINT(6*ran),:)
         END IF
         ptry(4:6) = MODULO(ptry(4:6), two_pi)
      ELSE IF (this%element_type_prm == "cartesian") THEN
!!$         IF (SQRT(SUM(ptry(1:3)**2)) > this%apriori_hcentric_dist_max_prm .or. &
!!$              SQRT(SUM(ptry(4:6)**2)) > this%apriori_velocity_max_prm) THEN
!!$            CALL randomNumber(ran)
!!$            ptry = p_init(1+NINT(6*(ran)),:)            
!!$         END IF
      END IF
      CALL NEW(orb, ptry, this%element_type_prm, &
           frame, this%t_inv_prm)
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / simplexOrbits / simtry", &
              "TRACE BACK (55)", 1)
         RETURN
      END IF
      CALL setParameters(orb, &
           dyn_model=this%dyn_model_prm, &
           perturbers=this%perturbers_prm, &
           integrator=this%integrator_prm, &
           integration_step=this%integration_step_prm)
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / simplexOrbits / simtry", &
              "TRACE BACK (45)", 1)
         RETURN
      END IF
      ytry = getChi2(this, orb)*rchi2_denom
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / simplexOrbits / simtry", &
              "TRACE BACK (60)", 1)
         RETURN
      END IF
      CALL NULLIFY(orb)

    END SUBROUTINE simtry

  END SUBROUTINE simplexOrbits





  !! *Description*:
  !!
  !! Statistical orbital ranging using two astrometric observations.
  !! Make use of the multiple orbit propagation scheme. CONTINUE TESTING!
  !! 
  !! Output: Probability density function for Cartesian position and
  !!         velocity at specified epoch in two-body approximation.
  !!
  !! Simplified:
  !!             - No random observation selection
  !!
  !! Returns error.
  !!
  SUBROUTINE statisticalRanging(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)                :: this
    TYPE (Orbit), DIMENSION(:), POINTER                  :: orb_arr_
    TYPE (Orbit), DIMENSION(:), POINTER                  :: orb_arr
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER   :: obsies_ccoords
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: scoords
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER   :: obs_scoords
    TYPE (CartesianCoordinates)                          :: obs_ccoord_topo1, &
         obs_ccoord_helio1, obs_ccoord_topo2, obs_ccoord_helio2
    TYPE (SphericalCoordinates)                          :: obs_scoord1, &
         obs_scoord2, obs_scoord_
    TYPE (Time)                                          :: t1, t2
    CHARACTER(len=DYN_MODEL_LEN)                         :: dyn_model
    CHARACTER(len=INTEGRATOR_LEN)                        :: integrator
    REAL(bp), DIMENSION(:,:,:,:), POINTER                :: partials_arr
    REAL(bp), DIMENSION(:,:,:), POINTER                  :: cov_matrices, &
         sphdev
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE              :: residuals
    REAL(bp), DIMENSION(:,:,:), POINTER                  :: information_matrix_obs
    REAL(bp), DIMENSION(:,:), ALLOCATABLE                :: rms, &
         rho_distribution, jacobians
    REAL(bp), DIMENSION(:), POINTER                      :: rho1, rho2, mjd_lt
    REAL(bp), DIMENSION(:), ALLOCATABLE                  :: cosdec0_arr, &
         pdf_arr, reg_apriori_arr, rchi2_arr, residual_vector
    REAL(bp), DIMENSION(6,6)                             :: information_matrix_elem, &
         jacobian_matrix
    REAL(bp), DIMENSION(6,2)                             :: bounds
    REAL(bp), DIMENSION(6)                               :: observed_coord, &
         computed_coord, elements
    REAL(bp), DIMENSION(3)                               :: position1, &
         position2, acc, r_ra_dec
    REAL(bp)                                             :: pdf_val, &
         pdf_relative_bound, apriori, jac_sph_inv, rho_comp1, rho_comp2, tdt, &
         rho_mid, a, ran, tmp, chi2, integration_step, &
         trials_per_orb, jac_car_kep, jac_equ_kep, obs_arc, &
         rang, distance, q, ftol
    INTEGER, DIMENSION(:,:), POINTER                     :: obs_pair_arr
    INTEGER, DIMENSION(:), ALLOCATABLE                   :: pair_histogram
    INTEGER, DIMENSION(6)                                :: n0, n0_
    INTEGER, DIMENSION(7)                                :: failed_flag
    INTEGER                                              :: err, i, j, k, &
         nobs, iorb, itrial, nobj, ncomb, ipair, err_verb_, &
         naccepted, norb
    LOGICAL, DIMENSION(:,:), ALLOCATABLE                 :: maskarr
    LOGICAL                                              :: first
    CHARACTER(len=DESIGNATION_LEN)                       :: object_id

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nobj = getNrOfObjects(this%obss)
    IF (nobj > 1 .AND. .NOT.this%multiple_obj_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "More than one object!", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    IF (nobs < 2) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Less than two observations available.", 1)
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(this%sor_deviates_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Generation window not set.", 1)
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(this%res_accept_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Window for accepted residuals not set.", 1)
       RETURN
    END IF
    IF (ALL(this%sor_rho_prm < 0.0_bp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Topocentric range intervals not set.", 1)
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(this%sor_pair_arr_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Observation pair not set.", 1)
       RETURN
    END IF
    IF (this%sor_norb_prm < 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Required number of sample orbits not set.", 1)
       RETURN
    END IF

    IF (this%pdf_ml_prm < 0.0_bp) THEN
       first = .TRUE.
       IF (this%pdf_ml_init_prm <= 0.0_bp) THEN
          this%pdf_ml_prm = pdf_ml_init
       ELSE
          this%pdf_ml_prm = this%pdf_ml_init_prm
       END IF
    ELSE
       first = .FALSE.
    END IF
    IF (.NOT.exist(this%t_inv_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Inversion epoch not set.", 1)
       RETURN
    END IF

    CALL getParameters(this, &
         dyn_model=dyn_model, &
         integration_step=integration_step, &
         integrator=integrator)!, &
    !         finite_diff=finite_diff)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF

    IF (integrator == "gauss-radau" .AND. dyn_model /= "2-body") THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Cannot use variational-equations approach for "//TRIM(integrator), 1)
       RETURN
    END IF

    object_id = getID(this%obss)
    obs_arc   = getObservationalTimespan(this%obss)

    ! ------------
    ! Initializing
    ! ------------

    ALLOCATE(orb_arr_(this%sor_norb_prm), &
         rho_distribution(this%sor_norb_prm,2), &
         residuals(this%sor_norb_prm,nobs,6), &
         pdf_arr(this%sor_norb_prm), &
         reg_apriori_arr(this%sor_norb_prm), &
         jacobians(this%sor_norb_prm,3), &
         rchi2_arr(this%sor_norb_prm), &
         rms(this%sor_norb_prm,6), &
         pair_histogram(SIZE(this%sor_pair_arr_prm,1)), &
         cosdec0_arr(nobs), &
         residual_vector(nobs*6), &
         maskarr(nobs,6), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "Could not allocate memory (1).", 1)
       RETURN
    END IF
    pdf_arr              = 0.0_bp
    rho_distribution = 0.0_bp
    residuals        = 0.0_bp
    reg_apriori_arr      = 0.0_bp
    jacobians        = 0.0_bp
    rms              = 0.0_bp
    rchi2_arr           = 0.0_bp
    pair_histogram   = 0
    acc              = 0.0_bp

    this%sor_rho_cmp(:,1) = HUGE(this%sor_rho_cmp)
    this%sor_rho_cmp(:,2) = -HUGE(this%sor_rho_cmp)

    ! Spherical toposentric observation coordinates:
    obs_scoords => getObservationSCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "TRACE BACK (20)", 1)
       DO i=1,SIZE(orb_arr_)
          CALL NULLIFY(orb_arr_(i))
       END DO
       DEALLOCATE(orb_arr_, stat=err)
       DEALLOCATE(rho_distribution, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(jacobians, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms, stat=err)
       DEALLOCATE(pair_histogram, stat=err)
       DEALLOCATE(cosdec0_arr, stat=err)
       DEALLOCATE(residual_vector, stat=err)
       DEALLOCATE(maskarr, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       RETURN
    END IF

    information_matrix_obs => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "TRACE BACK (30)", 1)
       DO i=1,SIZE(orb_arr_)
          CALL NULLIFY(orb_arr_(i))
       END DO
       DEALLOCATE(orb_arr_, stat=err)
       DEALLOCATE(rho_distribution, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(jacobians, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms, stat=err)
       DEALLOCATE(pair_histogram, stat=err)
       DEALLOCATE(cosdec0_arr, stat=err)
       DEALLOCATE(residual_vector, stat=err)
       DEALLOCATE(maskarr, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       RETURN
    END IF

    cov_matrices => getCovarianceMatrices(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "TRACE BACK (35)", 1)
       DO i=1,SIZE(orb_arr_)
          CALL NULLIFY(orb_arr_(i))
       END DO
       DEALLOCATE(orb_arr_, stat=err)
       DEALLOCATE(rho_distribution, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(jacobians, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms, stat=err)
       DEALLOCATE(pair_histogram, stat=err)
       DEALLOCATE(cosdec0_arr, stat=err)
       DEALLOCATE(residual_vector, stat=err)
       DEALLOCATE(maskarr, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(cov_matrices, stat=err)
       RETURN
    END IF

    ! Relative bound for probability density:
    pdf_relative_bound = EXP(-0.5_bp*this%dchi2_prm)

    ! Observation number counter (observation mask must be up-to-date!), 
    ! construct cosine array:
    DO i=1,6
       n0(i) = COUNT(this%obs_masks_prm(:,i))
    END DO

    DO i=1,nobs
       CALL rotateToEquatorial(obs_scoords(i))
       r_ra_dec = getPosition(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / statisticalRanging", &
               "TRACE BACK (40)", 1) 
          DO j=1,SIZE(orb_arr_)
             CALL NULLIFY(orb_arr_(j))
          END DO
          DEALLOCATE(orb_arr_, stat=err)
          DEALLOCATE(rho_distribution, stat=err)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(jacobians, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms, stat=err)
          DEALLOCATE(pair_histogram, stat=err)
          DEALLOCATE(cosdec0_arr, stat=err)
          DEALLOCATE(residual_vector, stat=err)
          DEALLOCATE(maskarr, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(cov_matrices, stat=err)
          RETURN
       END IF
       cosdec0_arr(i) = COS(r_ra_dec(3))
    END DO

    ! Extract equatorial observatory coordinates
    obsies_ccoords => getObservatoryCCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "TRACE BACK (45)", 1)
       DO i=1,SIZE(orb_arr_)
          CALL NULLIFY(orb_arr_(i))
       END DO
       DEALLOCATE(orb_arr_, stat=err)
       DEALLOCATE(rho_distribution, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(jacobians, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms, stat=err)
       DEALLOCATE(pair_histogram, stat=err)
       DEALLOCATE(cosdec0_arr, stat=err)
       DEALLOCATE(residual_vector, stat=err)
       DEALLOCATE(maskarr, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(cov_matrices, stat=err)
       DEALLOCATE(obsies_ccoords, stat=err)
       RETURN
    END IF

    IF (info_verb >= 1) THEN
       IF (this%uniform_pdf_prm) THEN
          WRITE(stdout,"(2X,A)") "WARNING: NOT using PDF in acceptance!"
       END IF
       IF (.NOT.this%regularization_prm) THEN
          WRITE(stdout,"(2X,A)") "WARNING: REGULARIZATION is OFF!"
       END IF
       IF (.NOT.this%jacobians_prm) THEN
          WRITE(stdout,"(2X,A)") "WARNING: JACOBIANS are NOT USED!"
       END IF
    END IF

    IF (first .AND. info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "- - - - - - - - - - - - ------------"
       WRITE(stdout,"(2X,A)") "STARTING MULTI-ORBIT RANGING . . . ."
       WRITE(stdout,"(2X,A)") "- - - - - - - - - - - - ------------"
       WRITE(stdout,"(1X)")
       WRITE(stdout,"(2X,A,A)")     "ID                      : ", &
            TRIM(object_id)
       !WRITE(stdout,"(2X,A,1X,I0)") "Nr of observations     :", &
       !     nobs
       WRITE(stdout,"(2X,A,1X,I0,A,I0,A,I0)") "Nr of incl. observations :", &
            n0(2),'+',n0(3),'/',nobs
       WRITE(stdout,"(2X,A,1X,F10.4,A)") "Observational time arc   :", &
            obs_arc," days"
       WRITE(stdout,"(2X,A)") "OPTIONS:"
       IF (.NOT. this%sor_random_obs_prm) THEN
          WRITE(stdout,"(2X,A,2(1X,I0))") "Observation pair used    :", &
               this%sor_pair_arr_prm
       ELSE 
          WRITE(stdout,"(2X,A)") "Random observation selection is used."
       END IF
       WRITE(stdout,"(2X,A,1X,I0)") "Nr of sample orbits      :", &
            this%sor_norb_prm
       WRITE(stdout,"(2X,A)") "Stdevs (R.A. & Dec. [arcsec]) and correlation :"
       DO i=1,nobs
          WRITE(stdout,"(4X,3(3X,F12.7))") &
               SQRT(cov_matrices(i,2,2))/rad_asec, &
               SQRT(cov_matrices(i,3,3))/rad_asec, &
               cov_matrices(i,2,3) / &
               SQRT(cov_matrices(i,2,2))*SQRT(cov_matrices(i,3,3))

       END DO
       WRITE(stdout,"(2X,A,A,A)") "PDF evaluated in ", &
            TRIM(this%element_type_prm)," elements"
       IF (this%uniform_pdf_prm) THEN
          WRITE(stdout,"(2X,A)") "NOT using PDF in acceptance!"
       ENDIF
       IF (.NOT. this%regularization_prm) WRITE(stdout,"(2X,A)",advance="no") "NOT"
       WRITE(stdout,"(2X,A,L2)") "Using regularization by Jeffreys"
       WRITE(stdout,"(2X,A,A)") "Epoch                   : ", &
            TRIM(getCalendarDateString(this%t_inv_prm,"tdt"))
       WRITE(stdout,'(2X,A,E10.4)') "Initial PDF maximum     : ", &
            this%pdf_ml_prm
       WRITE(stdout,"(2X,A,A,A)") "Using ",TRIM(dyn_model)," dynamical model"
       WRITE(stdout,"(2X,A,A,A)") "Using ", &
            TRIM(this%sor_2point_method_prm), " method"
       IF (this%sor_gaussian_pdf_prm) THEN
          WRITE(stdout,"(2X,A)") "Using Gaussian deviate for rho!"
       END IF
       WRITE(stdout,"(2X,A,4(1X,E10.4))") "Initial rho bounds (AU) :", &
            this%sor_rho_prm(1,1:2), this%sor_rho_prm(2,1:2)

       WRITE(stdout,"(1X)")
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A)") "Observations in MPC format: "
          CALL writeObservationFile(this%obss, stdout, "mpc")
          WRITE(stdout,"(2X,A)") "Generation window for each observation in seconds of arc"
          WRITE(stdout,"(2X,A)") "(RA_mean, RA_abs, Dec_mean, Dec_abs):"
          IF (.NOT. this%sor_random_obs_prm) THEN
             DO i=1,2
                j = this%sor_pair_arr_prm(1,i)
                WRITE(stdout,"(2X,A,I0,A,4(1X,F10.5))") "", j, ". obs:", &
                     this%sor_deviates_prm(j,2,1)/rad_asec, &
                     this%sor_deviates_prm(j,2,2)/rad_asec, &
                     this%sor_deviates_prm(j,3,1)/rad_asec, &
                     this%sor_deviates_prm(j,3,2)/rad_asec
             END DO
          ELSE
             DO i=1,nobs
                WRITE(stdout,"(2X,A,I0,A,4(1X,F10.5))") "", i, ". obs:", &
                     this%sor_deviates_prm(i,2,1)/rad_asec, &
                     this%sor_deviates_prm(i,2,2)/rad_asec, &
                     this%sor_deviates_prm(i,3,1)/rad_asec, &
                     this%sor_deviates_prm(i,3,2)/rad_asec
             END DO
          ENDIF
          WRITE(stdout,"(2X,A)") "Accepted residual window for each observation in seconds of arc"
          WRITE(stdout,"(2X,A)") "(RA_min, RA_max, Dec_min, Dec_max):"
          DO i=1,nobs
             WRITE(stdout,"(2X,A,I0,A,4(1X,F10.5))") "", i, ". obs:", &
                  -this%res_accept_prm(i,2)/rad_asec, &
                  this%res_accept_prm(i,2)/rad_asec, &
                  -this%res_accept_prm(i,3)/rad_asec, &
                  this%res_accept_prm(i,3)/rad_asec
          END DO
          DO i=1,nobs
             WRITE(stdout,"(2X,A,I0,A)") "Covariance matrix for ", i, &
                  ". observation: "
             CALL matrix_print(cov_matrices(i,:,:)/rad_asec**2.0_bp, &
                  stdout, errstr)
             IF (LEN_TRIM(errstr) /= 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "Could not print covariance matrix " // &
                     TRIM(errstr), 1) 
                RETURN
             END IF
          END DO
          WRITE(stdout,"(2X,A,6(I0,1X))") "Number of included observations " // &
               "(r,ra,dec,dr,dra,ddec): ", n0(1:6)
       END IF
    END IF

    !------------------------
    ! Monte Carlo selection :
    !------------------------

    iorb        = 0
    itrial      = 0
    trials_per_orb = 1.0_bp
    failed_flag = 0
    norb = 10*this%sor_norb_prm
    IF (norb > norb_simult_max) THEN
       norb = norb_simult_max
    END IF

    sor_main: DO WHILE (iorb < this%sor_norb_prm .AND. itrial < this%sor_ntrial_prm)

       naccepted = 0
       IF (.NOT.ASSOCIATED(orb_arr)) THEN
          ALLOCATE(orb_arr(norb), obs_pair_arr(norb,2), rho1(norb), &
               rho2(norb), sphdev(norb,2,3), mjd_lt(norb), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "Could not allocate memory (10).", 1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             RETURN
          END IF
       END IF

       i = 0
       sor_orb_gen: DO WHILE (i < norb .AND. itrial < this%sor_ntrial_prm)

          itrial = itrial + 1
          CALL NULLIFY(orb_arr(i+1))
          CALL NULLIFY(obs_ccoord_topo1)
          CALL NULLIFY(obs_ccoord_helio1)
          CALL NULLIFY(obs_ccoord_topo2)
          CALL NULLIFY(obs_ccoord_helio2)
          CALL NULLIFY(obs_scoord1)

          ! Random selection of observation pair
          IF (this%sor_random_obs_prm) THEN
             ncomb = SIZE(this%sor_pair_arr_prm,1)
             CALL randomNumber(ran)
             ipair = INT(ran*ncomb)+1
             obs_pair_arr(i+1,1:2) = this%sor_pair_arr_prm(ipair,1:2)
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,2(1X,I0))") &
                     "Random observation pair used:", obs_pair_arr(i+1,1:2)
             END IF
             ! Use predefined observation pair if not using random selection
          ELSE
             obs_pair_arr(i+1,1:2) = RESHAPE(this%sor_pair_arr_prm, (/ 2/))
          END IF
          obs_scoord1 = copy(obs_scoords(obs_pair_arr(i+1,1)))
          sphdev(i+1,1,:) = getPosition(obs_scoord1)

          ! Generate topocentric ranges and add deviates 
          ! to angular coordinates:
          bounds = 0.0_bp
          bounds(1,1) = 0.5_bp*(this%sor_rho_prm(1,1) + this%sor_rho_prm(1,2))
          bounds(1,2) = 0.5_bp*ABS(this%sor_rho_prm(1,2) - this%sor_rho_prm(1,1))
          bounds(2:3,1:2) = this%sor_deviates_prm(obs_pair_arr(i+1,1),2:3,1:2)
          CALL addUniformDeviate(obs_scoord1, bounds)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (50)", 1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             CALL NULLIFY(obs_scoord1)
             RETURN
          END IF
          IF (this%sor_gaussian_pdf_prm) THEN
             obs_scoord_ = copy(obs_scoord1)
             CALL randomGaussian(rang)
             distance = this%sor_rho_prm(1,1) + rang*this%sor_rho_prm(1,2)
             CALL NULLIFY(obs_scoord1)
             CALL NEW(obs_scoord1, distance, getLongitude(obs_scoord_), &
                  getLatitude(obs_scoord_), getTime(obs_scoord_))
             CALL NULLIFY(obs_scoord_)
          END IF

          ! addUniformDeviate does not return the actual deviations,
          ! accuracies computed in the end ar wrong!
          sphdev(i+1,1,:) = getPosition(obs_scoord1) - sphdev(i+1,1,:)
          sphdev(i+1,1,2) = sphdev(i+1,1,2)*cosdec0_arr(obs_pair_arr(i+1,1))
          rho1(i+1) = getDistance(obs_scoord1)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (55)", 1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             CALL NULLIFY(obs_scoord1)
             RETURN
          END IF
          IF (rho1(i+1) <= planetary_radii(3)) THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,F10.7,A)") &
                     "Failed (rho1 smaller than the Earth radius: ", &
                     rho1(i+1), " AU)"
             END IF
             CYCLE
          END IF
          rho_mid = 0.5_bp*(this%sor_rho_prm(2,1) + this%sor_rho_prm(2,2))
          bounds(1,1) = rho1(i+1) + rho_mid 
          IF (bounds(1,1) <= planetary_radii(3)) THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,F10.7,A)") &
                     "Failed (Mean of rho2 smaller than the Earth radius: ", &
                     bounds(1,1), " AU)"
             END IF
             CYCLE sor_orb_gen
          END IF
          bounds(2:3,1:2) = this%sor_deviates_prm(obs_pair_arr(i+1,2),2:3,1:2)
          rho2(i+1) = 0.0_bp
          ! Make sure rho2 is larger than zero length units
          DO WHILE (rho2(i+1) <= EPSILON(rho2(i+1)))
             bounds(1,2) = 0.5_bp*ABS(this%sor_rho_prm(2,2) - this%sor_rho_prm(2,1))
             CALL NULLIFY(obs_scoord2)
             obs_scoord2 = copy(obs_scoords(obs_pair_arr(i+1,2)))
             sphdev(i+1,2,:) = getPosition(obs_scoord2)
             CALL addUniformDeviate(obs_scoord2, bounds)
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "TRACE BACK (60)", 1)
                DO i=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(i))
                END DO
                DO i=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(i))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                CALL NULLIFY(obs_scoord1)
                CALL NULLIFY(obs_scoord2)
                RETURN
             END IF

             IF (this%sor_gaussian_pdf_prm) THEN
                obs_scoord_ = copy(obs_scoord2)
                CALL randomGaussian(rang)
                distance = rho1(i+1) + rang*this%sor_rho_prm(2,2)
                CALL NULLIFY(obs_scoord2)
                CALL NEW(obs_scoord2, distance, getLongitude(obs_scoord_), &
                     getLatitude(obs_scoord_), getTime(obs_scoord_))
             END IF

             rho2(i+1) = getDistance(obs_scoord2)
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "TRACE BACK (65)", 1)
                DO i=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(i))
                END DO
                DO i=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(i))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                CALL NULLIFY(obs_scoord1)
                CALL NULLIFY(obs_scoord2)
                RETURN
             END IF
          END DO
          sphdev(i+1,2,:) = getPosition(obs_scoord2) - sphdev(i+1,2,:)
          sphdev(i+1,2,2) = sphdev(i+1,2,2)*cosdec0_arr(obs_pair_arr(i+1,2))
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(2X,A,1X,2(F15.10,1X))") &
                  "Generated ranges are ", rho1(i+1), rho2(i+1)
          END IF

          CALL NEW(obs_ccoord_topo1, obs_scoord1)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (70)",1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             CALL NULLIFY(obs_scoord1)
             CALL NULLIFY(obs_scoord2)
             CALL NULLIFY(obs_ccoord_topo1)
             RETURN
          END IF

          CALL rotateToEcliptic(obs_ccoord_topo1)
          IF (info_verb >= 5) THEN
             position1 = getPosition(obs_scoord1)
             position2 = getPosition(obs_scoords(obs_pair_arr(i+1,1)))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "TRACE BACK (75)",1)
                DO i=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(i))
                END DO
                DO i=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(i))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                CALL NULLIFY(obs_scoord1)
                CALL NULLIFY(obs_scoord2)
                CALL NULLIFY(obs_ccoord_topo1)
                RETURN
             END IF
             WRITE(stdout,"(2X,A,1X,2(F15.10,1X))") &
                  "Spurious coordinates - Original observation 1:", & 
                  (position1(2:3)-position2(2:3))/rad_asec * &
                  cosdec0_arr(obs_pair_arr(i+1,1))
          END IF

          CALL NEW(obs_ccoord_topo2, obs_scoord2)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (80)",1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             CALL NULLIFY(obs_scoord1)
             CALL NULLIFY(obs_scoord2)
             CALL NULLIFY(obs_ccoord_topo1)
             CALL NULLIFY(obs_ccoord_topo2)
             RETURN
          END IF
          CALL rotateToEcliptic(obs_ccoord_topo2)
          IF (info_verb >= 5) THEN
             position1 = getPosition(obs_scoord2)
             position2 = getPosition(obs_scoords(obs_pair_arr(i+1,2)))
             WRITE(stdout,"(2X,A,1X,2(F15.10,1X))") &
                  "Spurious coordinates - Original observation 2:", &
                  (position1(2:3)-position2(2:3))/rad_asec * &
                  cosdec0_arr(obs_pair_arr(i+1,1))
             position1 = getPosition(obs_ccoord_topo2)
             WRITE(stdout,"(2X,A,1X,3(F15.10,1X))") &
                  "Cartesian coordinates, topocentric ecliptic 2:", position1
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "TRACE BACK (85)",1)
                DO i=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(i))
                END DO
                DO i=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(i))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                CALL NULLIFY(obs_scoord1)
                CALL NULLIFY(obs_scoord2)
                CALL NULLIFY(obs_ccoord_topo1)
                CALL NULLIFY(obs_ccoord_topo2)
                RETURN
             END IF
          END IF

          ! Ecliptic observatory coordinates
          CALL rotateToEcliptic(obsies_ccoords(obs_pair_arr(i+1,1)))
          CALL rotateToEcliptic(obsies_ccoords(obs_pair_arr(i+1,2)))

          ! Cartesian heliocentric coordinates:
          obs_ccoord_helio1 = &
               copy(obsies_ccoords(obs_pair_arr(i+1,1)) + obs_ccoord_topo1)
          obs_ccoord_helio2 = &
               copy(obsies_ccoords(obs_pair_arr(i+1,2)) + obs_ccoord_topo2)

          IF (info_verb >= 5) THEN
             t1 = getTime(obsies_ccoords(obs_pair_arr(i+1,1)))
             t2 = getTime(obs_ccoord_topo1)
             WRITE(stdout,"(2X,A,2(1X,F15.7))") &
                  "Heliocentric times 1:", &
                  getMJD(t1,"tdt"), getMJD(t2,"tdt")
             position1 = getPosition(obs_ccoord_helio1)
             position2 = getPosition(obsies_ccoords(obs_pair_arr(i+1,1)))
             WRITE(stdout,"(2X,A,1X,3(F15.10,1X))") &
                  "Heliocentric coordinates for observatory 1", position2
             CALL NULLIFY(t1)
             CALL NULLIFY(t2)
          END IF

          ! ------------------------------------------------
          ! Cartesian ecliptic position and velocity from 
          ! two cartesian ecliptic coordinates. P-iteration.
          ! ------------------------------------------------

          CALL estimateLightTime(obs_ccoord_helio1, rho1(i+1))
          CALL estimateLightTime(obs_ccoord_helio2, rho2(i+1))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (90)",1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             CALL NULLIFY(obs_scoord1)
             CALL NULLIFY(obs_scoord2)
             CALL NULLIFY(obs_ccoord_topo1)
             CALL NULLIFY(obs_ccoord_topo2)
             CALL NULLIFY(obs_ccoord_helio1)
             CALL NULLIFY(obs_ccoord_helio2)
             RETURN
          END IF
          t1 = getTime(obs_ccoord_helio1)
          mjd_lt(i+1) = getMJD(t1,"TDT")
          CALL NULLIFY(t1)
          t1 = getTime(obs_scoords(obs_pair_arr(i+1,1)))

          IF (info_verb >= 5) THEN
             !             t2 = getTime(obs_ccoord_helio2)
             !             tdt = getMJD(t2,"TDT")
             WRITE(stdout,"(2X,A,1X,F15.7)") "lt-time 1", mjd_lt(i+1)
             !             t1 = getTime(obs_scoords(obs_pair_arr(i+1,2)))
             !             t1 = getTime(obs_scoords(obs_pair_arr(i+1,1)))
             tdt = getMJD(t1,"TDT")
             WRITE(stdout,"(2X,A,1X,F15.7)") "original time 1", tdt
          END IF

          IF (this%sor_2point_method_prm == "n-body amoeba") THEN
             ! The tolerance is equal to one tenth of the tranverse
             ! uncertainty bounds (assuming that the transverse
             ! uncertainty is very small as compared to the
             ! topocentric distance of the object):
             position1 = getPosition(obs_scoord1)
             ftol = 0.05_bp * position1(1) * &
                  SQRT(SUM(ABS(bounds(2,1:2)))**2 + SUM(ABS(bounds(3,1:2)))**2)
          END IF

          ! Do not show errors:
          err_verb_ = err_verb
          err_verb = 0
          ! Find orbit candidate at the epoch of the first observation by
          ! using the chosen method to solve the 2-point boundary value
          ! problem:
          CALL NEW(orb_arr(i+1), obs_ccoord_helio1, obs_ccoord_helio2, &
               this%sor_2point_method_prm, this%apriori_a_max_prm, &
               ftol=ftol, perturbers=this%perturbers_prm, &
               integrator=this%integrator_prm, &
               integration_step=this%integration_step_prm)
          err_verb = err_verb_
          IF (error) THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A)") "Failed solving 2-point boundary value problem." 
             END IF
             failed_flag(1) = failed_flag(1) + 1
             error = .FALSE.
             CYCLE sor_orb_gen
          END IF

          IF (this%informative_apriori_prm) THEN
             ! Semimajor axis:
             a = -1.0_bp
             IF (this%apriori_a_min_prm >= 0.0_bp .OR. &
                  this%apriori_a_max_prm >= 0.0_bp) THEN
                a = getSemimajorAxis(orb_arr(i+1))
                IF (this%apriori_a_min_prm >= 0.0_bp .AND. &
                     a < this%apriori_a_min_prm) THEN
                   ! Semimajor axis too small
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F13.7,A)") &
                           "Failed (semimajor axis too small: ", a, " AU)"
                   END IF
                   failed_flag(2) = failed_flag(2) + 1
                   CYCLE sor_orb_gen
                END IF
                IF (this%apriori_a_max_prm >= 0.0_bp .AND. &
                     a > this%apriori_a_max_prm) THEN
                   ! Semimajor axis too large
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (semimajor axis too large: ", a, " AU)"
                   END IF
                   failed_flag(3) = failed_flag(3) + 1
                   CYCLE sor_orb_gen
                END IF
             END IF
             ! Periapsis distance:
             IF (this%apriori_periapsis_min_prm >= 0.0_bp .OR. &
                  this%apriori_periapsis_max_prm >= 0.0_bp) THEN
                IF (a >= 0.0_bp) THEN
                   CALL getPeriapsisDistance(orb_arr(i+1), q, a)
                ELSE
                   CALL getPeriapsisDistance(orb_arr(i+1), q)
                END IF
                IF (error) THEN
                   !call errorMessage("StochasticOrbit / statisticalRanging", &
                   !     "Error when computing periapsis distance.", 1)
                   !return
                   ! q not computed as orbit is hyperbolic
                   error = .FALSE.
                   CYCLE sor_orb_gen
                END IF
                ! Periapsis distance too small:
                IF (this%apriori_periapsis_min_prm >= 0.0_bp .AND. &
                     q < this%apriori_periapsis_min_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F13.7,A)") &
                           "Failed (periapsis distance too small: ", a, " AU)"
                   END IF
                   failed_flag(4) = failed_flag(4) + 1
                   CYCLE sor_orb_gen
                END IF
                ! Periapsis distance too large:
                IF (this%apriori_periapsis_max_prm >= 0.0_bp .AND. &
                     q > this%apriori_periapsis_max_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (periapsis distance too large: ", a, " AU)"
                   END IF
                   failed_flag(5) = failed_flag(5) + 1
                   CYCLE sor_orb_gen
                END IF
             END IF
             ! Apoapsis distance:
             IF (this%apriori_apoapsis_min_prm >= 0.0_bp .OR. &
                  this%apriori_apoapsis_max_prm >= 0.0_bp) THEN
                IF (a >= 0.0_bp) THEN
                   CALL getApoapsisDistance(orb_arr(i+1), Q, a)
                ELSE
                   CALL getApoapsisDistance(orb_arr(i+1), Q)
                END IF
                ! Apoapsis distance too small:
                IF (this%apriori_apoapsis_min_prm >= 0.0_bp .AND. &
                     Q < this%apriori_apoapsis_min_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F13.7,A)") &
                           "Failed (apoapsis distance too small: ", a, " AU)"
                   END IF
                   failed_flag(4) = failed_flag(4) + 1
                   CYCLE sor_orb_gen
                END IF
                ! Apoapsis distance too large:
                IF (this%apriori_apoapsis_max_prm >= 0.0_bp .AND. &
                     Q > this%apriori_apoapsis_max_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (apoapsis distance too large: ", a, " AU)"
                   END IF
                   failed_flag(5) = failed_flag(5) + 1
                   CYCLE sor_orb_gen
                END IF
             END IF
          END IF

          ! Propagate orbits first individually to a close-by joint epoch.
          CALL setParameters(orb_arr(i+1), &
               dyn_model=dyn_model, &
               perturbers=this%perturbers_prm, &
               integration_step=integration_step, &
               integrator=integrator)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (101)",1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             CALL NULLIFY(obs_scoord1)
             CALL NULLIFY(obs_scoord2)
             CALL NULLIFY(obs_ccoord_topo1)
             CALL NULLIFY(obs_ccoord_topo2)
             CALL NULLIFY(obs_ccoord_helio1)
             CALL NULLIFY(obs_ccoord_helio2)
             RETURN
          END IF
          CALL propagate(orb_arr(i+1), t1)
          IF (error) THEN 
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (95)",1)
             error = .FALSE.
             CYCLE sor_orb_gen
          END IF

          ! To allow comparison (e.g., derivatives) against the old code
          IF (TRIM(this%element_type_prm) == "keplerian") THEN
             CALL toKeplerian(orb_arr(i+1))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "Failed to change to Keplerian elements.",1)
                DO i=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(i))
                END DO
                DO i=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(i))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                CALL NULLIFY(obs_scoord1)
                CALL NULLIFY(obs_scoord2)
                CALL NULLIFY(obs_ccoord_topo1)
                CALL NULLIFY(obs_ccoord_topo2)
                CALL NULLIFY(obs_ccoord_helio1)
                CALL NULLIFY(obs_ccoord_helio2)
                RETURN
             END IF
          ELSE
             CALL rotateToEquatorial(orb_arr(i+1))
          END IF
          i = i + 1
       END DO sor_orb_gen
       CALL NULLIFY(obs_scoord1)
       CALL NULLIFY(obs_scoord2)
       CALL NULLIFY(obs_ccoord_topo1)
       CALL NULLIFY(obs_ccoord_topo2)
       CALL NULLIFY(obs_ccoord_helio1)
       CALL NULLIFY(obs_ccoord_helio2)

       IF (i < norb .AND. i /= 0) THEN
          ! All trials used without finding required number of orbits: 
          norb = i
          orb_arr => reallocate(orb_arr,norb)
          obs_pair_arr => reallocate(obs_pair_arr,norb,2)
          rho1 => reallocate(rho1,norb)
          rho2 => reallocate(rho2,norb)
          sphdev => reallocate(sphdev,norb,2,3)
          mjd_lt => reallocate(mjd_lt,norb)
       ELSE IF (i == 0) THEN
          ! All trials used and no orbits found:
          EXIT sor_main
       END IF

       ! Orbital elements at the specified epoch:
       ! DO THE MAIN PROPAGATION TO THE DESIRED EPOCH HERE
       CALL propagate(orb_arr, this%t_inv_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / statisticalRanging", &
               "TRACE BACK (105)",1)
          DO i=1,SIZE(orb_arr_)
             CALL NULLIFY(orb_arr_(i))
          END DO
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr_, stat=err)
          DEALLOCATE(rho_distribution, stat=err)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(jacobians, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms, stat=err)
          DEALLOCATE(pair_histogram, stat=err)
          DEALLOCATE(cosdec0_arr, stat=err)
          DEALLOCATE(residual_vector, stat=err)
          DEALLOCATE(maskarr, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(cov_matrices, stat=err)
          DEALLOCATE(obsies_ccoords, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(obs_pair_arr, stat=err)
          DEALLOCATE(rho1, stat=err)
          DEALLOCATE(rho2,stat=err)
          DEALLOCATE(sphdev, stat=err)
          DEALLOCATE(mjd_lt, stat=err)
          RETURN
       END IF

       ! ---------------------------------------------
       ! Topocentric positions and partial derivatives
       ! in two-body approximation.
       ! ---------------------------------------------

       ! Compute residuals and PDF value (via chi2) in order to either
       ! accept the orbit or discard it:
       IF (ASSOCIATED(scoords)) THEN
          DEALLOCATE(scoords, stat=err)
       END IF
       IF (this%regularization_prm .OR. this%jacobians_prm) THEN
          IF (ASSOCIATED(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          CALL getEphemerides(orb_arr, obsies_ccoords, scoords, &
               partials_arr=partials_arr)
       ELSE
          CALL getEphemerides(orb_arr, obsies_ccoords, scoords)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / statisticalRanging", &
               "TRACE BACK (110)", 1)
          DO i=1,SIZE(orb_arr_)
             CALL NULLIFY(orb_arr_(i))
          END DO
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr_, stat=err)
          DEALLOCATE(rho_distribution, stat=err)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(jacobians, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms, stat=err)
          DEALLOCATE(pair_histogram, stat=err)
          DEALLOCATE(cosdec0_arr, stat=err)
          DEALLOCATE(residual_vector, stat=err)
          DEALLOCATE(maskarr, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(cov_matrices, stat=err)
          DEALLOCATE(obsies_ccoords, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(obs_pair_arr, stat=err)
          DEALLOCATE(rho1, stat=err)
          DEALLOCATE(rho2,stat=err)
          DEALLOCATE(sphdev, stat=err)
          DEALLOCATE(mjd_lt, stat=err)
          DEALLOCATE(scoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          RETURN
       END IF

       sor_orb_acc: DO i=1,norb
          ! Sky-plane residuals and chi-squares:
          DO j=1,nobs
             observed_coord = getCoordinates(obs_scoords(j))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "TRACE BACK (115)",1)
                DO k=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(k))
                END DO
                DO k=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(k))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                DEALLOCATE(scoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                RETURN
             END IF
             computed_coord = getCoordinates(scoords(i,j))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "TRACE BACK (120)",1)
                DO k=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(k))
                END DO
                DO k=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(k))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                DEALLOCATE(scoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                RETURN
             END IF
             residuals(iorb+1,j,1:6) = observed_coord(1:6) - computed_coord(1:6)
             residuals(iorb+1,j,2) = residuals(iorb+1,j,2) * cosdec0_arr(j)

             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,3(1X,F15.10))") "observed pos.", observed_coord(1:3)
                WRITE(stdout,"(2X,A,3(1X,F15.10))") "computed pos.", computed_coord(1:3)
             END IF
             IF (this%regularization_prm .OR. this%jacobians_prm) THEN
                ! Multiply RA partials with cosine of observed declination:
                partials_arr(i,2,:,j) = partials_arr(i,2,:,j)*cosdec0_arr(j)
             END IF
          END DO

          maskarr = .FALSE.
          WHERE (this%obs_masks_prm .AND. ABS(residuals(iorb+1,:,:)) > this%res_accept_prm) 
             maskarr = .TRUE.
          END WHERE
          IF (info_verb >= 5) THEN
             DO j=1,nobs
                WRITE(stdout,"(2X,A,2(1X,F10.5))") "O-C residuals (RA, Dec):", &
                     residuals(iorb+1,j,2:3)/rad_asec
             END DO
             WRITE(stdout,"(2X,A,I0,A,I0)") &
                  "No of omitted obs/included obs: ", &
                  COUNT(maskarr),"/",n0(2)
          END IF
          IF (COUNT(maskarr) > 0) THEN
             ! Residuals are too large for at least one observation.
             failed_flag(6) = failed_flag(6) + 1
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A)") &
                     "Failed (residuals are too large)"
             END IF
             CYCLE sor_orb_acc
          END IF

          ! Compute chi2:
          chi2 = chi_square(residuals(iorb+1,:,:), information_matrix_obs, this%obs_masks_prm, errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "TRACE BACK (125)", 1)
             DO j=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(j))
             END DO
             DO j=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(j))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             DEALLOCATE(scoords, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             RETURN          
          END IF

          ! Jeffrey's apriori:
          ! Monitor the matrix inversion?
          IF (this%regularization_prm) THEN
             ! Sigma_elements^(-1) = A^T Sigma_obs^(-1) A, where A is the
             ! partial derivatives matrix of ephemerides wrt elements:
             information_matrix_elem = 0.0_bp
             DO j=1,nobs
                information_matrix_elem = information_matrix_elem + &
                     MATMUL(MATMUL(TRANSPOSE(partials_arr(i,1:6,1:6,j)), &
                     information_matrix_obs(j,1:6,1:6)), &
                     partials_arr(i,1:6,1:6,j))
             END DO
             apriori = SQRT(ABS(determinant(information_matrix_elem, errstr)))
             IF (LEN_TRIM(errstr) /= 0) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "Unsuccessful computation of determinant of orbital element " // &
                     "information matrix " // TRIM(errstr), 1)
                errstr = ""
                IF (err_verb >= 1) THEN
                   CALL matrix_print(information_matrix_elem, stderr, errstr)
                END IF
                errstr = ""
                CYCLE
             END IF
          ELSE
             apriori = 1.0_bp
          END IF

          IF (this%jacobians_prm) THEN
             ! Determinant of Jacobian between topocentric
             ! coordinates (inverse problem coordinates)
             ! and orbital parameters required for output
             ! ("Topocentric Wrt Cartesian/Keplerian"):
             jacobian_matrix(1:3,:) = partials_arr(i,1:3,:,obs_pair_arr(i,1)) / &
                  cosdec0_arr(obs_pair_arr(i,1))
             jacobian_matrix(4:6,:) = partials_arr(i,1:3,:,obs_pair_arr(i,2)) / &
                  cosdec0_arr(obs_pair_arr(i,2))
             jac_sph_inv = ABS(determinant(jacobian_matrix, errstr))
             IF (LEN_TRIM(errstr) /= 0) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "Unsuccessful computation of determinant of orbital element " // &
                     "jacobian matrix " // TRIM(errstr), 1)
                errstr = ""
                IF (err_verb >= 1) THEN
                   CALL matrix_print(jacobian_matrix, stderr, errstr)
                END IF
                errstr = ""
                CYCLE
             END IF

             ! Determinant of Jacobian between Cartesian and Keplerian
             ! orbital elements ("Cartesian Wrt Keplerian"):
             CALL partialsCartesianWrtKeplerian(orb_arr(i), &
                  jacobian_matrix, "equatorial")
             IF (error .AND. TRIM(this%element_type_prm) /= "keplerian") THEN
                error = .FALSE.
                jac_car_kep = -1.0_bp
             ELSE IF (error) THEN
                CALL errorMessage("StochasticOrbit / statisticalRanging", &
                     "Unsuccessful computation of jacobian matrix " // &
                     "between Cartesian and Keplerian elements.", 1)
                DO j=1,SIZE(orb_arr_)
                   CALL NULLIFY(orb_arr_(j))
                END DO
                DO j=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(j))
                END DO
                DEALLOCATE(orb_arr_, stat=err)
                DEALLOCATE(rho_distribution, stat=err)
                DEALLOCATE(residuals, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(jacobians, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms, stat=err)
                DEALLOCATE(pair_histogram, stat=err)
                DEALLOCATE(cosdec0_arr, stat=err)
                DEALLOCATE(residual_vector, stat=err)
                DEALLOCATE(maskarr, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(cov_matrices, stat=err)
                DEALLOCATE(obsies_ccoords, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(obs_pair_arr, stat=err)
                DEALLOCATE(rho1, stat=err)
                DEALLOCATE(rho2,stat=err)
                DEALLOCATE(sphdev, stat=err)
                DEALLOCATE(mjd_lt, stat=err)
                DEALLOCATE(scoords, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                RETURN
             ELSE
                jac_car_kep = ABS(determinant(jacobian_matrix, errstr)) 
                IF (LEN_TRIM(errstr) /= 0) THEN
                   CALL errorMessage("StochasticOrbit / statisticalRanging", &
                        "Unsuccessful computation of determinant of" // &
                        " jacobian matrix between Cartesian and Keplerian" // &
                        " elements:" // TRIM(errstr), 1)
                   errstr = ""
                   IF (err_verb >= 1) THEN
                      CALL matrix_print(jacobian_matrix, stderr, errstr)
                   END IF
                   errstr = ""
                   DO j=1,SIZE(orb_arr_)
                      CALL NULLIFY(orb_arr_(j))
                   END DO
                   DO j=1,SIZE(orb_arr)
                      CALL NULLIFY(orb_arr(j))
                   END DO
                   DEALLOCATE(orb_arr_, stat=err)
                   DEALLOCATE(rho_distribution, stat=err)
                   DEALLOCATE(residuals, stat=err)
                   DEALLOCATE(pdf_arr, stat=err)
                   DEALLOCATE(reg_apriori_arr, stat=err)
                   DEALLOCATE(jacobians, stat=err)
                   DEALLOCATE(rchi2_arr, stat=err)
                   DEALLOCATE(rms, stat=err)
                   DEALLOCATE(pair_histogram, stat=err)
                   DEALLOCATE(cosdec0_arr, stat=err)
                   DEALLOCATE(residual_vector, stat=err)
                   DEALLOCATE(maskarr, stat=err)
                   DEALLOCATE(obs_scoords, stat=err)
                   DEALLOCATE(information_matrix_obs, stat=err)
                   DEALLOCATE(cov_matrices, stat=err)
                   DEALLOCATE(obsies_ccoords, stat=err)
                   DEALLOCATE(orb_arr, stat=err)
                   DEALLOCATE(obs_pair_arr, stat=err)
                   DEALLOCATE(rho1, stat=err)
                   DEALLOCATE(rho2,stat=err)
                   DEALLOCATE(sphdev, stat=err)
                   DEALLOCATE(mjd_lt, stat=err)
                   DEALLOCATE(scoords, stat=err)
                   DEALLOCATE(partials_arr, stat=err)
                   RETURN
                END IF
             END IF

             ! Determinant of Jacobian between equinoctial and
             ! Keplerian orbital elements ("Equinoctial Wrt
             ! Keplerian"):
             elements = getElements(orb_arr(i), "keplerian")
             jac_equ_kep = 0.5_bp*elements(2) * &
                  SIN(0.5_bp*elements(3)) / COS(0.5_bp*elements(3))**3
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,F15.12)") "Chi2 new: ", chi2
                WRITE(stdout,"(2X,A,3(1X,F10.5))") "PDF values:", pdf_val, &
                     pdf_val/this%pdf_ml_prm, pdf_relative_bound
             END IF
          ELSE
             jac_sph_inv = 1.0_bp
             jac_car_kep = 1.0_bp
             jac_equ_kep = 1.0_bp
          END IF

          ! Probability density function (note that the '- nobs' term
          ! is there for practical reasons):
          pdf_val = apriori*EXP(-0.5_bp*(chi2 - COUNT(this%obs_masks_prm)))/jac_sph_inv

          IF (.NOT.this%uniform_pdf_prm .AND. &
               pdf_val / this%pdf_ml_prm < pdf_relative_bound) THEN
             ! The PDF is used and its value is not acceptable.
             failed_flag(7) = failed_flag(7) + 1
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,1X,E10.5)") &
                     "Failed (PDF value not acceptable)", pdf_val
             END IF
             CYCLE sor_orb_acc
          END IF

          ! Update metrics for computational accuracy:
          rho_comp1 = getDistance(scoords(i,obs_pair_arr(i,1)))
          rho_comp2 = getDistance(scoords(i,obs_pair_arr(i,2)))

          ! 1) Topocentric ranges
          IF (ABS(rho_comp1 - rho1(i)) > acc(1)) THEN
             acc(1) = ABS(rho_comp1 - rho1(i))
          END IF
          IF (ABS(rho_comp2 - rho2(i)) > acc(1)) THEN
             acc(1) = ABS(rho_comp2 - rho2(i))
          END IF

          ! 2) R.A. and Dec.
          tmp = ABS(residuals(iorb+1,obs_pair_arr(i,1),2) + sphdev(i,1,2))
          IF (tmp > acc(2)) THEN
             acc(2) = tmp
          END IF
          tmp = ABS(residuals(iorb+1,obs_pair_arr(i,2),2) + sphdev(i,2,2))
          IF (tmp > acc(2)) THEN
             acc(2) = tmp
          END IF
          tmp = ABS(residuals(iorb+1,obs_pair_arr(i,1),3) + sphdev(i,1,3))
          IF (tmp > acc(3)) THEN
             acc(3) = tmp
          END IF
          tmp = ABS(residuals(iorb+1,obs_pair_arr(i,2),3) + sphdev(i,2,3))
          IF (tmp > acc(3)) THEN
             acc(3) = tmp
          END IF

          ! Update topocentric distance bounds:
          IF (rho1(i) < this%sor_rho_cmp(1,1)) THEN
             this%sor_rho_cmp(1,1) = rho1(i)
          END IF
          IF (rho1(i) > this%sor_rho_cmp(1,2)) THEN
             this%sor_rho_cmp(1,2) = rho1(i)
          END IF
          IF (rho2(i) - rho1(i) < this%sor_rho_cmp(2,1)) THEN
             this%sor_rho_cmp(2,1) = rho2(i) - rho1(i)
          END IF
          IF (rho2(i) - rho1(i) > this%sor_rho_cmp(2,2)) THEN
             this%sor_rho_cmp(2,2) = rho2(i) - rho1(i)
          END IF

          ! Orbital element storage:
          iorb = iorb + 1
          naccepted = naccepted + 1
          ! - accepted sample orbit:
          orb_arr_(iorb) = copy(orb_arr(i))
          ! - p.d.f.:
          pdf_arr(iorb) = pdf_val
          ! - regularizing apriori:
          reg_apriori_arr(iorb) = apriori
          ! - determinant of jacobian between computed positions and
          !   orbital elements:
          jacobians(iorb,1) = jac_sph_inv
          ! - determinant of jacobian between Cartesian and Keplerian
          !   elements:
          jacobians(iorb,2) = jac_car_kep
          ! - determinant of jacobian between Equinoctial and Keplerian
          !   elements:
          jacobians(iorb,3) = jac_equ_kep
          ! Nonlinear reduced chi2 (= chi2/(number of measurements))
          rchi2_arr(iorb) = chi2/REAL(COUNT(this%obs_masks_prm),bp)

          n0_ = n0
          WHERE (n0_ == 0)
             n0_ = 1
          END WHERE
          ! - rms:
          rms(iorb,:) = SQRT(SUM(residuals(iorb,:,:)**2.0_bp,dim=1,mask=this%obs_masks_prm)/n0_)
          ! - generated distance to object at first epoch 
          rho_distribution(iorb,1) = rho_comp1
          ! - generated distance to object at second epoch (NB! Relative
          !   to the distance at the first epoch.)
          rho_distribution(iorb,2) = rho_comp2
          ! Update histogram of observation pairs used as boundary values
          ! when solving the 2-point boundary value problem:
          DO j=1,SIZE(this%sor_pair_arr_prm,dim=1)
             IF (obs_pair_arr(i,1) == this%sor_pair_arr_prm(j,1) .AND. &
                  obs_pair_arr(i,2) == this%sor_pair_arr_prm(j,2)) THEN
                pair_histogram(j) = pair_histogram(j) + 1
             END IF
          END DO
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(2X,A,I0)") "Sample orbit ", iorb
          END IF
          IF (iorb == this%sor_norb_prm) THEN
             EXIT
          END IF

       END DO sor_orb_acc

       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,2(1X,I0))") &
               "Nr of accepted orbits and accepted trial orbits:", &
               naccepted, norb
       END IF
       IF (naccepted == 0) THEN
          norb = norb*100
       ELSE
          norb = NINT((this%sor_norb_prm - iorb)*(1.1_bp*norb/naccepted))
       END IF
       IF (norb > norb_simult_max) THEN
          norb = norb_simult_max
       END IF
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,1X,I0)") "Nr of orbits to be generated:", norb
       END IF
       IF (norb /= SIZE(orb_arr,dim=1) .OR. &
            iorb == this%sor_norb_prm .OR. &
            itrial == this%sor_ntrial_prm) THEN
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, obs_pair_arr, rho1, rho2, sphdev, mjd_lt, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "Could not deallocate memory (10).", 1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             DEALLOCATE(scoords, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             RETURN
          END IF
       END IF
       IF (this%regularization_prm .OR. this%jacobians_prm) THEN
          DEALLOCATE(scoords, partials_arr, stat=err)
       ELSE
          DEALLOCATE(scoords, stat=err)
       END IF
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / statisticalRanging", &
               "Could not deallocate memory (12).", 1)
          DO i=1,SIZE(orb_arr_)
             CALL NULLIFY(orb_arr_(i))
          END DO
          IF (ASSOCIATED(orb_arr)) THEN
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
          END IF
          DEALLOCATE(orb_arr_, stat=err)
          DEALLOCATE(rho_distribution, stat=err)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(jacobians, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms, stat=err)
          DEALLOCATE(pair_histogram, stat=err)
          DEALLOCATE(cosdec0_arr, stat=err)
          DEALLOCATE(residual_vector, stat=err)
          DEALLOCATE(maskarr, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(cov_matrices, stat=err)
          DEALLOCATE(obsies_ccoords, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(obs_pair_arr, stat=err)
          DEALLOCATE(rho1, stat=err)
          DEALLOCATE(rho2,stat=err)
          DEALLOCATE(sphdev, stat=err)
          DEALLOCATE(mjd_lt, stat=err)
          DEALLOCATE(scoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          RETURN
       END IF

    END DO sor_main

    ! Final number of orbits and unnormalized a posteriori 
    ! probability density values:
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,I0)") "Final number of orbits : ", iorb
       WRITE(stdout,"(2X,A,I0)") "Final number of trials : ", itrial
       WRITE(stdout,"(2X,A,F6.2)") "Failed-%  - total       : ", &
            SUM(failed_flag(1:5))/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - 2-point     : ", &
            failed_flag(1)/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - a too small : ", & 
            failed_flag(2)/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - a too large : ", &
            failed_flag(3)/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - q too small : ", & 
            failed_flag(4)/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - q too large : ", & 
            failed_flag(5)/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - residuals   : ", & 
            failed_flag(6)/REAL(itrial)*100.0_bp
       WRITE(stdout,"(2X,A,F6.2)") "          - p.d.f.      : ", & 
            failed_flag(7)/REAL(itrial)*100.0_bp
    END IF

    this%sor_norb_cmp   = iorb
    this%sor_ntrial_cmp = itrial
    IF (this%sor_norb_cmp > 0) THEN
       IF (ASSOCIATED(this%orb_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL NULLIFY(this%orb_arr_cmp(i))
          END DO
          DEALLOCATE(this%orb_arr_cmp, stat=err)
       END IF
       DEALLOCATE(this%sor_rho_arr_cmp, stat=err)
       DEALLOCATE(this%res_arr_cmp, stat=err)
       DEALLOCATE(this%pdf_arr_cmp, stat=err)
       DEALLOCATE(this%reg_apr_arr_cmp, stat=err)
       DEALLOCATE(this%jac_arr_cmp, stat=err)
       DEALLOCATE(this%rchi2_arr_cmp, stat=err)
       DEALLOCATE(this%rms_arr_cmp, stat=err)
       DEALLOCATE(this%sor_pair_histo_prm, stat=err)
       ALLOCATE(this%orb_arr_cmp(this%sor_norb_cmp), &
            this%sor_rho_arr_cmp(this%sor_norb_cmp,2), &
            this%res_arr_cmp(this%sor_norb_cmp,nobs,6), &
            this%pdf_arr_cmp(this%sor_norb_cmp), &
            this%reg_apr_arr_cmp(this%sor_norb_cmp), &
            this%jac_arr_cmp(this%sor_norb_cmp,3), &
            this%rchi2_arr_cmp(this%sor_norb_cmp), &
            this%rms_arr_cmp(this%sor_norb_cmp,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / statisticalRanging", &
               "Could not allocate memory (3).", 1)
          DO i=1,SIZE(orb_arr_)
             CALL NULLIFY(orb_arr_(i))
          END DO
          IF (ASSOCIATED(orb_arr)) THEN
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
          END IF
          DEALLOCATE(orb_arr_, stat=err)
          DEALLOCATE(rho_distribution, stat=err)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(jacobians, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms, stat=err)
          DEALLOCATE(pair_histogram, stat=err)
          DEALLOCATE(cosdec0_arr, stat=err)
          DEALLOCATE(residual_vector, stat=err)
          DEALLOCATE(maskarr, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(cov_matrices, stat=err)
          DEALLOCATE(obsies_ccoords, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(obs_pair_arr, stat=err)
          DEALLOCATE(rho1, stat=err)
          DEALLOCATE(rho2,stat=err)
          DEALLOCATE(sphdev, stat=err)
          DEALLOCATE(mjd_lt, stat=err)
          DEALLOCATE(scoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          RETURN
       END IF

       DO i=1, this%sor_norb_cmp
          this%orb_arr_cmp(i) = copy(orb_arr_(i))
       END DO
       this%pdf_arr_cmp       = pdf_arr(1:iorb)
       i = MAXLOC(this%pdf_arr_cmp,dim=1)
       CALL NULLIFY(this%orb_ml_cmp)
       this%orb_ml_cmp        = copy(this%orb_arr_cmp(i)) 
       this%reg_apr_arr_cmp   = reg_apriori_arr(1:iorb)
       this%jac_arr_cmp       = jacobians(1:iorb,1:3)
       this%rchi2_arr_cmp = rchi2_arr(1:iorb)
       this%rms_arr_cmp       = rms(1:iorb,:)
       this%sor_rho_arr_cmp   = rho_distribution(1:iorb,:)
       this%res_arr_cmp       = residuals(1:iorb,:,:)
       IF (this%sor_random_obs_prm) THEN
          ALLOCATE(this%sor_pair_histo_prm(SIZE(this%sor_pair_arr_prm,1)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / statisticalRanging", &
                  "Could not allocate memory (4).", 1)
             DO i=1,SIZE(orb_arr_)
                CALL NULLIFY(orb_arr_(i))
             END DO
             IF (ASSOCIATED(orb_arr)) THEN
                DO i=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(i))
                END DO
             END IF
             DEALLOCATE(orb_arr_, stat=err)
             DEALLOCATE(rho_distribution, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(jacobians, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms, stat=err)
             DEALLOCATE(pair_histogram, stat=err)
             DEALLOCATE(cosdec0_arr, stat=err)
             DEALLOCATE(residual_vector, stat=err)
             DEALLOCATE(maskarr, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(cov_matrices, stat=err)
             DEALLOCATE(obsies_ccoords, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(obs_pair_arr, stat=err)
             DEALLOCATE(rho1, stat=err)
             DEALLOCATE(rho2,stat=err)
             DEALLOCATE(sphdev, stat=err)
             DEALLOCATE(mjd_lt, stat=err)
             DEALLOCATE(scoords, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             RETURN
          END IF
          this%sor_pair_histo_prm = pair_histogram
       END IF
       DO i=1,SIZE(orb_arr_)
          CALL NULLIFY(orb_arr_(i))
       END DO
       IF (ASSOCIATED(orb_arr)) THEN
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
       END IF
       DEALLOCATE(orb_arr_, stat=err)
       DEALLOCATE(rho_distribution, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(jacobians, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms, stat=err)
       DEALLOCATE(pair_histogram, stat=err)
       DEALLOCATE(cosdec0_arr, stat=err)
       DEALLOCATE(residual_vector, stat=err)
       DEALLOCATE(maskarr, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(cov_matrices, stat=err)
       DEALLOCATE(obsies_ccoords, stat=err)
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(obs_pair_arr, stat=err)
       DEALLOCATE(rho1, stat=err)
       DEALLOCATE(rho2, stat=err)
       DEALLOCATE(sphdev, stat=err)
       DEALLOCATE(mjd_lt, stat=err)
       DEALLOCATE(scoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / statisticalRanging", &
            "No sample orbits found!", 1)
       DO i=1,SIZE(orb_arr_)
          CALL NULLIFY(orb_arr_(i))
       END DO
       IF (ASSOCIATED(orb_arr)) THEN
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
       END IF
       DEALLOCATE(orb_arr_, stat=err)
       DEALLOCATE(rho_distribution, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(jacobians, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms, stat=err)
       DEALLOCATE(pair_histogram, stat=err)
       DEALLOCATE(cosdec0_arr, stat=err)
       DEALLOCATE(residual_vector, stat=err)
       DEALLOCATE(maskarr, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(cov_matrices, stat=err)
       DEALLOCATE(obsies_ccoords, stat=err)
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(obs_pair_arr, stat=err)
       DEALLOCATE(rho1, stat=err)
       DEALLOCATE(rho2,stat=err)
       DEALLOCATE(sphdev, stat=err)
       DEALLOCATE(mjd_lt, stat=err)
       DEALLOCATE(scoords, stat=err)
       DEALLOCATE(partials_arr, stat=err)
       RETURN
    END IF

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,2(1X,E10.4))") "Init/max pdf values    :", &
            this%pdf_ml_prm, MAXVAL(this%pdf_arr_cmp)
       WRITE(stdout,"(2X,A,4(1X,E10.4))") "Min/max rhos, @ end    :" , &
            this%sor_rho_cmp(1,1:2), this%sor_rho_cmp(2,1:2)
       acc(2:3) = acc(2:3)/rad_asec
       WRITE(stdout,"(2X,A,3(1X,E10.4),A)") "Accur. (rho/R.A./Dec.) :", acc, " (AU/as/as)"
       WRITE(stdout,"(2X,A)") ""
    END IF

  END SUBROUTINE statisticalRanging





  !!  *Description*:
  !!
  !! Finds initial values for the full inversion by including the
  !! observations one at a time and perfoming a partial inversion.
  !! A limited amount of required sample orbits and trial orbits are
  !! used during the initiation phase.
  !! 
  !!
  SUBROUTINE stepwiseRanging(this, nobs_max)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)     :: this
    INTEGER, INTENT(in)                       :: nobs_max

    TYPE (Observations)                       :: obss
    TYPE (Observation), DIMENSION(:), POINTER :: obs_arr
    TYPE (StochasticOrbit)                    :: storb
    CHARACTER(len=4)                          :: str1, str2
    REAL(bp), DIMENSION(4)                    :: rho
    REAL(bp)                                  :: pdf_ml_, pdf_ml_final, ddchi2
    INTEGER(ibp)                              :: nobs, nobs_max_, err, &
         i, j, k

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    obs_arr => getObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    nobs = SIZE(obs_arr,dim=1)
    IF (nobs < 2) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "Number of observations is not sufficient. At least two required.", 1)
       DEALLOCATE(obs_arr, stat=err)
       RETURN
    END IF
    IF (nobs_max < 0 .OR. nobs_max > nobs) THEN
       nobs_max_ = nobs
    ELSE
       nobs_max_ = nobs_max
    END IF
    CALL NEW(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(obs_arr, stat=err)
       CALL NULLIFY(obss)
       RETURN
    END IF
    CALL addObservation(obss, obs_arr(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "TRACE BACK (15)", 1)
       DEALLOCATE(obs_arr, stat=err)
       CALL NULLIFY(obss)
       RETURN
    END IF

    i = 1
    j = 2
    k = nobs
    rho(1) = this%sor_rho_prm(1,1)
    rho(2) = this%sor_rho_prm(1,2)
    rho(3) = this%sor_rho_prm(2,1)
    rho(4) = this%sor_rho_prm(2,2)
    DO WHILE (k-j /= -1 .AND. i < nobs_max_ .AND. i < nobs)
       i = i + 1
       IF ((-1)**i < 0) THEN
          CALL addObservation(obss, obs_arr(j))
          j = j + 1
       ELSE
          CALL addObservation(obss, obs_arr(k))
          k = k - 1
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / stepwiseRanging", &
               "TRACE BACK (20)", 1)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(obss)
          RETURN
       END IF
       CALL NULLIFY(storb)
       CALL NEW(storb, obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / stepwiseRanging", &
               "TRACE BACK (22)", 1)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(obss)
          CALL NULLIFY(storb)
          RETURN
       END IF
       CALL setParameters(storb, &
            dyn_model=this%dyn_model_prm, &
            perturbers=this%perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm, &
            t_inv=this%t_inv_prm, &
            element_type=this%element_type_prm, &
            multiple_objects=this%multiple_obj_prm, &
            outlier_rejection=this%outlier_rejection_prm, &
            uniform_pdf=this%uniform_pdf_prm, &
            regularized_pdf=this%regularization_prm, &
            jacobians_pdf=this%jacobians_prm, &
            accept_multiplier=this%accept_multiplier_prm, &
            apriori_a_max=this%apriori_a_max_prm, &
            apriori_a_min=this%apriori_a_min_prm, &
            apriori_periapsis_max=this%apriori_periapsis_max_prm, &
            apriori_periapsis_min=this%apriori_periapsis_min_prm, &
            apriori_apoapsis_max=this%apriori_apoapsis_max_prm, &
            apriori_apoapsis_min=this%apriori_apoapsis_min_prm, &
            apriori_rho_min=this%apriori_rho_min_prm, &
            sor_norb=this%sor_norb_sw_prm, &
                                ! Number of trial orbits increases with number of observations:
                                !sor_ntrial=this%sor_ntrial_sw_prm*i, &
            sor_ntrial=this%sor_ntrial_sw_prm, &
            sor_rho1_l=rho(1), &
            sor_rho1_u=rho(2), &
            sor_rho2_l=rho(3), &
            sor_rho2_u=rho(4), &
            sor_niter=1, &
            sor_generat_multiplier=this%sor_generat_multiplier_prm, &
            sor_2point_method=this%sor_2point_method_sw_prm, &
            sor_iterate_bounds=this%sor_iterate_bounds_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / stepwiseRanging", &
               "TRACE BACK (27)", 1)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(obss)
          CALL NULLIFY(storb)
          RETURN
       END IF
       IF (i == 2) THEN
          CALL autoStatisticalRanging(storb)
       ELSE
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(2X,A,1X,I0,A,I0)") "Nr of observations     :", i, "/", nobs
          END IF
          storb%pdf_ml_prm = pdf_ml_
          CALL statisticalRanging(storb)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / stepwiseRanging", &
               "TRACE BACK (30)", 1)
          DEALLOCATE(obs_arr, stat=err)
          CALL NULLIFY(obss)
          CALL NULLIFY(storb)
          RETURN
       ELSE
          this%sor_norb_sw_cmp = storb%sor_norb_cmp
          this%sor_ntrial_sw_cmp = storb%sor_ntrial_cmp
          IF (this%sor_norb_sw_cmp < 2) THEN
             CALL toString(this%sor_norb_sw_cmp, str1, error)
             !             CALL toString(i, str2, error)
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  TRIM(str1) // "sample orbits found when " // &
                  TRIM(str2) // " observations were included." , 1)
             error = .TRUE.
             DEALLOCATE(obs_arr, stat=err)
             CALL NULLIFY(obss)
             CALL NULLIFY(storb)
             RETURN
          END IF
          CALL setRangeBounds(storb)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  "TRACE BACK (33", 1)
             DEALLOCATE(obs_arr, stat=err)
             CALL NULLIFY(obss)
             CALL NULLIFY(storb)
             RETURN
          END IF
          rho(1:4) = getRangeBounds(storb)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  "TRACE BACK (35)", 1)
             DEALLOCATE(obs_arr, stat=err)
             CALL NULLIFY(obss)
             CALL NULLIFY(storb)
             RETURN
          END IF
          pdf_ml_ = MAXVAL(storb%pdf_arr_cmp)
       END IF
    END DO
    DEALLOCATE(obs_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "Could not deallocate memory.", 1)
       CALL NULLIFY(obss)
       CALL NULLIFY(storb)
       RETURN
    END IF
    CALL NULLIFY(obss)

    this%sor_niter_cmp = 0
    this%sor_rho_histo_cmp = 1
    IF (nobs == nobs_max_) THEN
       ddchi2 = 10.0_bp
       DO WHILE (this%sor_niter_cmp < this%sor_niter_prm .AND. &
            (this%sor_rho_histo_cmp > 0 .OR. ddchi2 > 2.0_bp .OR. &
            this%sor_norb_cmp < this%sor_norb_prm))
          CALL setParameters(this, &
               sor_rho1_l=rho(1), &
               sor_rho1_u=rho(2), &
               sor_rho2_l=rho(3), &
               sor_rho2_u=rho(4))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  "TRACE BACK (40)", 1)
             CALL NULLIFY(storb)
             RETURN
          END IF
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(2X,A,1X,I0)") "Nr of observations     :", i
          END IF
          this%pdf_ml_prm = pdf_ml_
          CALL statisticalRanging(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  "Subsequent iteration failed.", 1)
             CALL NULLIFY(storb)
             RETURN
          END IF
          this%sor_niter_cmp = this%sor_niter_cmp + 1
          IF (this%sor_norb_cmp < this%sor_norb_prm .AND. err_verb >= 2) THEN
             WRITE(stderr,"(A,I0)") "Warning by stepwiseRanging:" // &
                  " Number of sample orbits too small: ", this%sor_norb_cmp
             WRITE(stderr,"(1X)")
          END IF
          pdf_ml_final = MAXVAL(this%pdf_arr_cmp)
          ddchi2  = ABS(2.0_bp*LOG(pdf_ml_final/this%pdf_ml_prm))
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(2X,A,1X,F12.4)") "Delta dchi2            :", &
                  ddchi2
          END IF
          CALL setRangeBounds(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  "TRACE BACK (45)", 1)
             CALL NULLIFY(storb)
             RETURN
          END IF
          rho(1:4) = getRangeBounds(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  "TRACE BACK (50)", 1)
             CALL NULLIFY(storb)
             RETURN
          END IF
          pdf_ml_ = MAX(MAXVAL(this%pdf_arr_cmp),this%pdf_ml_prm)
       END DO
       !WRITE(stderr,*) this%sor_norb_cmp
       !WRITE(stderr,*) this%sor_ntrial_cmp
    ELSE
       CALL NULLIFY(this)
       this = copy(storb)
    END IF
    CALL NULLIFY(storb)

  END SUBROUTINE stepwiseRanging





  SUBROUTINE toCartesian_SO(this, frame)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this
    CHARACTER(len=*), INTENT(in) :: frame

    CHARACTER(len=FRAME_LEN) :: frame_
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    INTEGER :: i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / toCartesian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    frame_ = TRIM(frame)
    CALL locase(frame_, error)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / toCartesian", &
            "The frame string contains forbidden characters.", 1)
       RETURN
    END IF
    IF (frame_ /= "equatorial" .AND. &
         frame_ /= "ecliptic") THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / toCartesian", &
            "Frame " // TRIM(frame_) // " not recognized.", 1)
    END IF

    IF (this%element_type_prm == "cartesian") THEN
       IF (ASSOCIATED(this%orb_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL toCartesian(this%orb_arr_cmp(i), frame_)
          END DO
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "cartesian", frame_)
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          i = MAXLOC(this%pdf_arr_cmp,dim=1)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
          CALL toCartesian(this%orb_ml_cmp, frame_)
       END IF
    ELSE IF (this%element_type_prm == "cometary" .OR. &
         this%element_type_prm == "keplerian") THEN
       IF (ASSOCIATED(this%orb_arr_cmp) .AND. &
            ASSOCIATED(this%pdf_arr_cmp) .AND. &
            ASSOCIATED(this%jac_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL toCartesian(this%orb_arr_cmp(i), frame_)
             !this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i)/this%jac_arr_cmp(i,2)
          END DO
          pdf_arr => getPDFValues(this, "cartesian")
          this%pdf_arr_cmp = pdf_arr
          DEALLOCATE(pdf_arr)
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "cartesian", frame_)
          this%cov_type_prm = "cartesian"
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          i = MAXLOC(this%pdf_arr_cmp,dim=1)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
          CALL toCartesian(this%orb_ml_cmp, frame_)
       END IF
       this%element_type_prm = "cartesian"
    ELSE
       error =.TRUE.
       CALL errorMessage("StochasticOrbit / toCartesian", &
            "Cannot convert from " // TRIM(this%element_type_prm) // &
            " to Cartesian elements.", 1)
       RETURN
    END IF

  END SUBROUTINE toCartesian_SO





  SUBROUTINE toCometary_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this

    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    INTEGER :: i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / toCometary", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%element_type_prm == "cometary") THEN
       RETURN
    END IF

    IF (this%element_type_prm == "cartesian" .OR. &
         this%element_type_prm == "keplerian") THEN
       IF (ASSOCIATED(this%orb_arr_cmp) .AND. &
            ASSOCIATED(this%pdf_arr_cmp) .AND. &
            ASSOCIATED(this%jac_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL toCometary(this%orb_arr_cmp(i))
             !this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i)*this%jac_arr_cmp(i,2)
          END DO
          pdf_arr => getPDFValues(this, "cometary")
          this%pdf_arr_cmp = pdf_arr
          DEALLOCATE(pdf_arr)
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "cometary")
          this%cov_type_prm = "cometary"
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          i = MAXLOC(this%pdf_arr_cmp,dim=1)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
          CALL toCometary(this%orb_ml_cmp)          
       END IF
       this%element_type_prm = "cometary"
    ELSE
       error =.TRUE.
       CALL errorMessage("StochasticOrbit / toCometary", &
            "Cannot convert from " // TRIM(this%element_type_prm) // &
            " to Cometary elements.", 1)
       RETURN
    END IF

  END SUBROUTINE toCometary_SO





  SUBROUTINE toKeplerian_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this

    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    INTEGER :: i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / toKeplerian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%element_type_prm == "keplerian") THEN
       RETURN
    END IF

    IF (this%element_type_prm == "cartesian" .OR. &
         this%element_type_prm == "cometary") THEN
       IF (ASSOCIATED(this%orb_arr_cmp) .AND. &
            ASSOCIATED(this%pdf_arr_cmp) .AND. &
            ASSOCIATED(this%jac_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL toKeplerian(this%orb_arr_cmp(i))
             !this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i)*this%jac_arr_cmp(i,2)
          END DO
          pdf_arr => getPDFValues(this, "keplerian")
          this%pdf_arr_cmp = pdf_arr
          DEALLOCATE(pdf_arr)
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "keplerian")
          this%cov_type_prm = "keplerian"
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          i = MAXLOC(this%pdf_arr_cmp,dim=1)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
          CALL toKeplerian(this%orb_ml_cmp)          
       END IF
       this%element_type_prm = "keplerian"
    ELSE
       error =.TRUE.
       CALL errorMessage("StochasticOrbit / toKeplerian", &
            "Cannot convert from " // TRIM(this%element_type_prm) // &
            " to Keplerian elements.", 1)
       RETURN
    END IF

  END SUBROUTINE toKeplerian_SO




  !! *Description*:
  !!
  !! Updates the StochasticOrbit between subsequent calls
  !! of statistical ranging.
  !!
  !! Output: - New range bounds, currently from the 3-sigma cutoff values 
  !!           of the range probability density.
  !!         - Maximum probability density value
  !!         - Optional outlier detection: recognizes outlier observations
  !!           and optionally flags them (obs_masks -> .false.).
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE updateRanging(this, automatic)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout)  :: this
    LOGICAL, INTENT(in), OPTIONAL         :: automatic

    CHARACTER(len=128) :: str, frmt
    REAL(bp), DIMENSION(:,:), POINTER :: stdevs
    REAL(bp), DIMENSION(:), ALLOCATABLE :: ra_mean, dec_mean
    REAL(bp) :: ra_std, dec_std
    INTEGER, DIMENSION(:), ALLOCATABLE :: nr_array
    INTEGER :: i, nobs, err, nr_of_omitted
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks
    LOGICAL, DIMENSION(:), ALLOCATABLE :: maskarr
    LOGICAL :: automatic_

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / updateRanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    ALLOCATE(ra_mean(nobs), dec_mean(nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / updateRanging", &
            "Could not allocate memory (5).", 1)
       RETURN
    END IF

    IF (PRESENT(automatic)) THEN
       automatic_ = automatic
    ELSE
       automatic_ = .FALSE.
    END IF

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "- - - - - - - - - "
       WRITE(stdout,"(2X,A)") "UPDATING . . . . ."
       WRITE(stdout,"(2X,A)") "- - - - - - - - - "
    END IF

    ! Update range intervals
    CALL setRangeBounds(this)
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,4(1X,E10.4))") "Real rho bounds:", &
            this%sor_rho_cmp(1,1:2), this%sor_rho_cmp(2,1:2)
       WRITE(stdout,"(2X,A,4(1X,E10.4))") "New rho bounds: ", &
            this%sor_rho_prm(1,1:2), this%sor_rho_prm(2,1:2)
    END IF

    DO i=1,nobs
       CALL moments(this%res_arr_cmp(:,i,2), mean=ra_mean(i), &
            std_dev=ra_std, error=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / updateRanging", &
               "Error in moment computation for RA " // &
               TRIM(errstr), 1)
          RETURN
       END IF
       CALL moments(this%res_arr_cmp(:,i,3), mean=dec_mean(i), &
            std_dev=dec_std, error=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / updateRanging", &
               "Error in moment computation for Dec " // &
               TRIM(errstr), 1)
          RETURN
       END IF
    END DO

    IF (info_verb >= 3) THEN
       WRITE(stdout,"(2X,A)") "Residual offsets, R.A. & Dec:"
       DO i=1,nobs
          WRITE(stdout,"(2X,A,I0,A,2(1X,F10.5))") "Obs no ", i, ":", &
               ra_mean(i)/rad_asec, dec_mean(i)/rad_asec
       END DO
    END IF

    ! Update maximum probability density value and generating windows
    this%pdf_ml_prm = MAXVAL(this%pdf_arr_cmp)

    ! Recognize outlier observations
    IF (this%outlier_multiplier_prm > 0.0_bp) THEN
       ALLOCATE(nr_array(nobs), maskarr(nobs), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / updateRanging", &
               "Could not allocate memory (10).", 1)
          RETURN
       END IF
       nr_array = RESHAPE((/ (i, i = 1,nobs) /), (/ nobs /))
       obs_masks => getObservationMasks(this%obss, use_notes=.TRUE.)
       stdevs => getStandardDeviations(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / updateRanging", &
               "TRACE BACK ()", 1)
          RETURN
       END IF
       maskarr = .FALSE.
       DO i=1,nobs
          IF (ABS(ra_mean(i)) > this%outlier_multiplier_prm*stdevs(i,2) .OR. &
               ABS(dec_mean(i)) > this%outlier_multiplier_prm*stdevs(i,3)) THEN
             obs_masks(i,1:6) = .FALSE.
             maskarr(i) = .TRUE.
          END IF
       END DO

       nr_of_omitted = COUNT(maskarr)
       IF (this%outlier_rejection_prm) THEN
          ! Remove outlier observations
          this%obs_masks_prm = obs_masks
          IF (info_verb >= 2) THEN
             IF (nr_of_omitted == 0) THEN
                WRITE(stdout,"(2X,A)") "No of omitted: 0"
             ELSE IF (nr_of_omitted == 1) THEN
                WRITE(stdout,"(2X,A,2(I0,A))") "No of omitted: ", &
                     nr_of_omitted,"(",PACK(nr_array,maskarr),")"
             ELSE IF (nr_of_omitted >= 2) THEN
                str = " "
                CALL toString(nr_of_omitted-1, str, error)
                IF (error) THEN
                   CALL errorMessage("StochasticOrbit / updateRanging", &
                        "Error in conversion from integer to character.", 1)
                   RETURN
                END IF
                frmt = "(A,I0,A,I0," // TRIM(str) // "(1X,I0),A)"
                WRITE(stdout,TRIM(frmt)) "No of omitted: ", &
                     nr_of_omitted," (",PACK(nr_array,maskarr),")"
             END IF
          END IF
       ELSE IF (nr_of_omitted > 0) THEN
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(1X)")
             WRITE(stdout,"(2X,A,I0)") "WARNING from updateRanging: " // &
                  "potential outlier observations detected! No: ", nr_of_omitted
          END IF
       END IF
       IF (info_verb >= 3 .AND. this%outlier_rejection_prm) THEN
          WRITE(stdout,"(2X,A)") &
               "Obs mask after update + mean RA and Dec:"
          DO i=1,SIZE(this%obs_masks_prm,dim=1)
             WRITE(stdout,"(2X,A,I0,A,6(1X,L1),2(1X,F10.5))") &
                  "Obs #", i, ":", this%obs_masks_prm(i,:), &
                  ra_mean(i)/rad_asec, dec_mean(i)/rad_asec
          END DO
       END IF
       DEALLOCATE(nr_array, obs_masks, maskarr, stdevs, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / updateRanging", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
    END IF

    DEALLOCATE(ra_mean, dec_mean, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(ra_mean, stat=err)
       DEALLOCATE(dec_mean, stat=err)
       CALL errorMessage("StochasticOrbit / updateRanging", &
            "Could not deallocate memory (10).", 1)
       RETURN
    END IF

  END SUBROUTINE updateRanging




END MODULE StochasticOrbit_cl
