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
!! *Class*description*: 
!!  
!! Type and routines for stochastic orbits (that is, orbital elements
!! with their uncertainties). Contains the algorithms for the
!! [statistical orbital] ranging method and the least-squares method.
!!
!! @author MG, JV, KM, DO 
!! @version 2019-10-29
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
  !! Default bound for equivalent chi2-difference. Change it through
  !! setParameters!
  !! E.g., 50.0 corresponds to probability mass (1 - 4.7e-9).  
  !! E.g., 30.0 corresponds to probability mass (1 - ?).
  !! E.g., 20.1 corresponds to probability mass 0.9973 which
  !!  corresponds to the 1D 3-sigma confidence limits
  REAL(bp), PARAMETER, PRIVATE :: dchi2       = 20.1_bp    
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
  PRIVATE :: getID_SO
  PRIVATE :: getCovarianceMatrix_SO
  PRIVATE :: getParameters_SO
  PRIVATE :: getPeriapsisDistance_SO
  PRIVATE :: getPhaseAngle_SO_pdf
  PRIVATE :: getPhaseAngle_SO_point
  PRIVATE :: getPhaseAngles_SO
  PRIVATE :: getObservationMasks_SO
  PRIVATE :: getObservations_SO
  PRIVATE :: getRangeBounds_SO
  PRIVATE :: getResiduals_SO_obss
  PRIVATE :: getResiduals_SO_orb
  PRIVATE :: getResiduals_SO_orb_arr
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
     TYPE (Orbit)                        :: orb_ml_prm
     TYPE (Observations)                 :: obss
     CHARACTER(len=DESIGNATION_LEN)      :: id_prm                 = " "
     CHARACTER(len=ELEMENT_TYPE_LEN)     :: cov_type_prm           = " "
     CHARACTER(len=ELEMENT_TYPE_LEN)     :: element_type_prm       = "cartesian"
     REAL(bp), DIMENSION(:,:,:), POINTER :: res_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: cov_ml_cmp             => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: rms_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: res_accept_prm         => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: jac_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: rchi2_arr_cmp          => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: pdf_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: reg_apr_arr_cmp        => NULL()
     REAL(bp)                            :: chi2_min_init_prm      = -1.0_bp
     REAL(bp)                            :: dchi2_prm              = dchi2
     REAL(bp)                            :: chi2_min_prm           = -1.0_bp
     REAL(bp)                            :: chi2_min_cmp           = -1.0_bp
     REAL(bp)                            :: accept_multiplier_prm  = -1.0_bp
     REAL(bp)                            :: outlier_multiplier_prm = -1.0_bp
     INTEGER, DIMENSION(:), POINTER      :: repetition_arr_cmp     => NULL()
     INTEGER                             :: center_prm       = 11

     LOGICAL, DIMENSION(:,:), POINTER    :: obs_masks_prm          => NULL()
     LOGICAL                             :: outlier_rejection_prm  = .FALSE.
     LOGICAL                             :: regularization_prm     = .TRUE.
     LOGICAL                             :: jacobians_prm          = .TRUE.
     LOGICAL                             :: multiple_obj_prm       = .FALSE.
     LOGICAL                             :: is_initialized_prm     = .FALSE.
     LOGICAL                             :: dchi2_rejection_prm        = .TRUE.
     LOGICAL                             :: mh_acceptance_prm        = .FALSE.
     LOGICAL                             :: generat_gaussian_deviates_prm = .TRUE.
     REAL(bp)                            :: generat_multiplier_prm = -1.0_bp

     ! Bayesian informative apriori assumptions
     REAL(bp)                            :: apriori_a_max_prm             = -1.0_bp
     REAL(bp)                            :: apriori_a_min_prm             = -1.0_bp
     REAL(bp)                            :: apriori_periapsis_max_prm     = -1.0_bp
     REAL(bp)                            :: apriori_periapsis_min_prm     = -1.0_bp
     REAL(bp)                            :: apriori_apoapsis_max_prm      = -1.0_bp
     REAL(bp)                            :: apriori_apoapsis_min_prm      = -1.0_bp
     REAL(bp)                            :: apriori_rho_max_prm           = -1.0_bp
     REAL(bp)                            :: apriori_rho_min_prm           = -1.0_bp
     REAL(bp)                            :: apriori_hcentric_dist_max_prm = -1.0_bp
     REAL(bp)                            :: apriori_hcentric_dist_min_prm = -1.0_bp
     REAL(bp)                            :: apriori_velocity_max_prm      = -1.0_bp
     LOGICAL                             :: informative_apriori_prm       = .FALSE.

     ! Parameters for propagation:
     CHARACTER(len=DYN_MODEL_LEN)        :: dyn_model_prm        = "2-body"
     CHARACTER(len=INTEGRATOR_LEN)       :: integrator_prm       = "gauss-radau"
     REAL(bp), DIMENSION(:), POINTER     :: finite_diff_prm      => NULL()
     REAL(bp)                            :: integration_step_prm = 1.0_bp
     LOGICAL, DIMENSION(10)              :: perturbers_prm       = .FALSE.
     LOGICAL                             :: ast_perturbers_prm   = .FALSE.

     ! Parameters for statistical ranging:
     CHARACTER(len=64)                   :: sor_2point_method_prm      = "continued fraction"
     CHARACTER(len=64)                   :: sor_2point_method_sw_prm   = "continued fraction"
     REAL(bp), DIMENSION(:,:,:), POINTER :: sor_deviates_prm           => NULL()
     REAL(bp), DIMENSION(:,:), POINTER   :: sor_rho_arr_cmp            => NULL()
     REAL(bp), DIMENSION(:), POINTER     :: sor_pair_histo_prm         => NULL()
     REAL(bp), DIMENSION(2,2)            :: sor_rho_prm                = -1.0_bp
     REAL(bp), DIMENSION(2,2)            :: sor_rho_cmp                = -1.0_bp
     !REAL(bp), DIMENSION(2)              :: sor_rho_mean_prm           = -1.0_bp
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
     ! Define number of burn-in orbits required in MCMC ranging
     INTEGER                             :: sor_burnin_prm             = 1
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

     ! Parameters for VOMCMC:
     REAL(bp), DIMENSION(:,:), POINTER   :: vomcmc_map_cmp           => NULL()
     REAL(bp), DIMENSION(6,2)            :: vomcmc_scaling_prm       = -1.0_bp
     REAL(bp), DIMENSION(6,2)            :: vomcmc_scaling_cmp       = -1.0_bp
     INTEGER                             :: vomcmc_norb_prm          = -1
     INTEGER                             :: vomcmc_norb_iter_prm     = -1
     INTEGER                             :: vomcmc_norb_cmp          = -1
     INTEGER                             :: vomcmc_ntrial_prm        = -1
     INTEGER                             :: vomcmc_ntrial_iter_prm   = -1
     INTEGER                             :: vomcmc_ntrial_cmp        = -1
     INTEGER                             :: vomcmc_niter_prm         = -1
     INTEGER                             :: vomcmc_niter_cmp         = -1
     INTEGER                             :: vomcmc_nmap_prm          = -1
     LOGICAL, DIMENSION(6)               :: vomcmc_mapping_mask_prm  = &
          (/ .TRUE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /)
     LOGICAL, DIMENSION(6,2)             :: vomcmc_scaling_ready_cmp = &
          .FALSE.

     ! Parameters for the least-squares fitting:
     REAL(bp)                            :: ls_corr_fac_prm         = 1.0_bp
     REAL(bp)                            :: ls_rchi2_diff_tresh_prm = 0.00001_bp
     REAL(bp)                            :: ls_rchi2_acceptable_prm = 1.5_bp
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
     REAL(bp)                            :: smplx_similarity_tol_prm = 0.0001
     INTEGER                             :: smplx_niter_prm = 1000
     INTEGER                             :: smplx_niter_cmp
     LOGICAL                             :: smplx_force_prm = .FALSE.

     ! Parameters for MCMC observation sampling:
     INTEGER                             :: os_norb_prm = 500
     INTEGER                             :: os_ntrial_prm = 5000
     INTEGER                             :: os_sampling_type_prm = 1

  END TYPE StochasticOrbit


  INTERFACE NEW
     MODULE PROCEDURE new_SO
     MODULE PROCEDURE new_SO_observations
     MODULE PROCEDURE new_SO_orb_cov
     MODULE PROCEDURE new_SO_orb_arr
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_SO
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_SO
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_SO
  END INTERFACE exist

  INTERFACE getApoapsisDistance
     MODULE PROCEDURE getApoapsisDistance_SO
  END INTERFACE getApoapsisDistance

  INTERFACE getChi2
     MODULE PROCEDURE getChi2_matrix
     MODULE PROCEDURE getChi2_this_orb
  END INTERFACE getChi2

  INTERFACE getCovarianceMatrix
     MODULE PROCEDURE getCovarianceMatrix_SO
  END INTERFACE getCovarianceMatrix

  INTERFACE getEphemerides
     MODULE PROCEDURE getEphemerides_SO
  END INTERFACE getEphemerides

  INTERFACE getEphemeris
     MODULE PROCEDURE getEphemeris_SO
  END INTERFACE getEphemeris

  INTERFACE getID
     MODULE PROCEDURE getID_SO
  END INTERFACE getID

  INTERFACE getParameters
     MODULE PROCEDURE getParameters_SO
  END INTERFACE getParameters

  INTERFACE getPeriapsisDistance
     MODULE PROCEDURE getPeriapsisDistance_SO
  END INTERFACE getPeriapsisDistance

  INTERFACE getPhaseAngle
     MODULE PROCEDURE getPhaseAngle_SO_pdf
     MODULE PROCEDURE getPhaseAngle_SO_point
  END INTERFACE getPhaseAngle

  INTERFACE getPhaseAngles
     MODULE PROCEDURE getPhaseAngles_SO
  END INTERFACE getPhaseAngles

  INTERFACE getObservationMasks
     MODULE PROCEDURE getObservationMasks_SO
  END INTERFACE getObservationMasks

  INTERFACE getObservations
     MODULE PROCEDURE getObservations_SO
  END INTERFACE getObservations

  INTERFACE getRangeBounds
     MODULE PROCEDURE getRangeBounds_SO
  END INTERFACE getRangeBounds

  INTERFACE getResiduals
     MODULE PROCEDURE getResiduals_SO_obss
     MODULE PROCEDURE getResiduals_SO_orb
     MODULE PROCEDURE getResiduals_SO_orb_arr
  END INTERFACE getResiduals

  INTERFACE getResults
     MODULE PROCEDURE getResults_SO
  END INTERFACE getResults

  INTERFACE getRMS
     MODULE PROCEDURE getRMS_single
  END INTERFACE getRMS

  INTERFACE getSampleOrbit
     MODULE PROCEDURE getSampleOrbit_SO
  END INTERFACE getSampleOrbit

  INTERFACE getStandardDeviations
     MODULE PROCEDURE getStandardDeviations_SO
  END INTERFACE getStandardDeviations

  INTERFACE leastSquares
     MODULE PROCEDURE leastSquares_SO
  END INTERFACE leastSquares

  INTERFACE levenbergMarquardt
     MODULE PROCEDURE levenbergMarquardt_SO
  END INTERFACE levenbergMarquardt

  INTERFACE getTime
     MODULE PROCEDURE getTime_SO
  END INTERFACE getTime

  INTERFACE propagate
     MODULE PROCEDURE propagate_SO
  END INTERFACE propagate

  INTERFACE setAcceptanceWindow
     MODULE PROCEDURE setAcceptanceWindow_sigma
  END INTERFACE setAcceptanceWindow

  INTERFACE setObservationMask
     MODULE PROCEDURE setObservationMask_one
     MODULE PROCEDURE setObservationMask_all
     MODULE PROCEDURE setObservationMask_all_notes
  END INTERFACE setObservationMask

  INTERFACE setObservationPair
     MODULE PROCEDURE setObservationPair_default
     MODULE PROCEDURE setObservationPair_pair
  END INTERFACE setObservationPair

  INTERFACE setParameters
     MODULE PROCEDURE setParameters_SO
  END INTERFACE setParameters

  INTERFACE setRangeBounds
     MODULE PROCEDURE setRangeBounds_3sigma
     MODULE PROCEDURE setRangeBounds_values
  END INTERFACE setRangeBounds

  INTERFACE toCartesian
     MODULE PROCEDURE toCartesian_SO
  END INTERFACE toCartesian

  INTERFACE toCometary
     MODULE PROCEDURE toCometary_SO
  END INTERFACE toCometary

  INTERFACE toKeplerian
     MODULE PROCEDURE toKeplerian_SO
  END INTERFACE toKeplerian

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

    INTEGER :: err_verb_

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
    err_verb_ = err_verb
    err_verb = 0
    CALL setObservationPair(this)
    err_verb = err_verb_
    IF (error) THEN
       CALL warningMessage("StochasticOrbit / new", &
            "Observation pair could not be selected. Possibly " // &
            "because only one observation has been provided.", 1)
       error = .FALSE.
    END IF

  END SUBROUTINE new_SO_observations





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE new_SO_orb_cov(this, orb, cov, cov_type, &
       element_type, obss, id)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)     :: this
    TYPE (Orbit), INTENT(in)                  :: orb
    CHARACTER(len=*), INTENT(in)              :: cov_type
    REAL(bp), DIMENSION(6,6), INTENT(in)      :: cov
    CHARACTER(len=*), INTENT(in), OPTIONAL    :: element_type
    TYPE (Observations), INTENT(in), OPTIONAL :: obss
    CHARACTER(len=*), INTENT(in), OPTIONAL    :: id
    INTEGER                                   :: err, err_verb_

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
       err_verb_ = err_verb
       err_verb = 0
       CALL setObservationPair(this)
       err_verb = err_verb_
       IF (error) THEN
          CALL warningMessage("StochasticOrbit / new", &
               "Observation pair could not be selected. Possibly " // &
               "because only one observation has been provided.", 1)
          error = .FALSE.
       END IF
    END IF
    IF (PRESENT(id)) THEN
       this%id_prm = id 
    END IF

  END SUBROUTINE new_SO_orb_cov





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!

  SUBROUTINE new_SO_orb_arr(this, orb_arr, pdf_arr, element_type, &
       jac_arr, reg_apr_arr, rchi2_arr, repetition_arr, obss, id)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)          :: this
    TYPE (Orbit), DIMENSION(:), INTENT(in)         :: orb_arr
    TYPE (Observations), INTENT(in), OPTIONAL      :: obss
    CHARACTER(len=*), INTENT(in)                   :: element_type
    CHARACTER(len=*), INTENT(in), OPTIONAL         :: id
    REAL(bp), DIMENSION(:), INTENT(in)             :: pdf_arr
    REAL(bp), DIMENSION(:,:), INTENT(in), OPTIONAL :: jac_arr
    REAL(bp), DIMENSION(:), INTENT(in), OPTIONAL   :: reg_apr_arr
    REAL(bp), DIMENSION(:), INTENT(in), OPTIONAL   :: rchi2_arr
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL    :: repetition_arr

    INTEGER :: i, err, err_verb_

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
       this%rchi2_arr_cmp = rchi2_arr
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
       this%reg_apr_arr_cmp = reg_apr_arr
    END IF
    IF (PRESENT(repetition_arr)) THEN
       IF (.NOT.ALL(repetition_arr == 0)) THEN
          ALLOCATE(this%repetition_arr_cmp(this%sor_norb_cmp), &
               stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / new", &
                  "Could not allocate pointer (25).", 1)
             RETURN
          END IF
          this%repetition_arr_cmp = repetition_arr
       END IF
    END IF
    this%element_type_prm = element_type
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
    this%is_initialized_prm = .TRUE.

    IF (PRESENT(obss)) THEN
       err_verb_ = err_verb
       err_verb = 0
       CALL setObservationPair(this)
       err_verb = err_verb_
       IF (error) THEN
          CALL warningMessage("StochasticOrbit / new", &
               "Observation pair could not be selected. Possibly " // &
               "because only one observation has been provided.", 1)
          error = .FALSE.
       END IF
    END IF
    IF (PRESENT(id)) THEN
       this%id_prm = id 
    END IF

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
    CALL NULLIFY(this%orb_ml_prm)
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
    IF (ASSOCIATED(this%repetition_arr_cmp)) THEN
       DEALLOCATE(this%repetition_arr_cmp, stat=err)
    END IF
    IF (ASSOCIATED(this%obs_masks_prm)) THEN
       DEALLOCATE(this%obs_masks_prm, stat=err)
    END IF
    this%chi2_min_init_prm = -1.0_bp
    this%chi2_min_prm = -1.0_bp
    this%dchi2_prm = dchi2
    this%accept_multiplier_prm = -1.0_bp
    this%outlier_rejection_prm = .FALSE.
    this%regularization_prm = .TRUE.
    this%jacobians_prm = .TRUE.
    this%multiple_obj_prm = .FALSE.

    ! Bayesian apriori parameters
    this%apriori_a_max_prm = -1.0_bp
    this%apriori_a_min_prm = -1.0_bp
    this%apriori_periapsis_max_prm = -1.0_bp
    this%apriori_periapsis_min_prm = -1.0_bp
    this%apriori_apoapsis_max_prm = -1.0_bp
    this%apriori_apoapsis_min_prm = -1.0_bp
    this%apriori_rho_max_prm = -1.0_bp
    this%apriori_rho_min_prm = -1.0_bp
    this%informative_apriori_prm = .FALSE.

    ! Propagation parameters
    this%dyn_model_prm = "2-body"
    this%integrator_prm = "gauss-radau"
    IF (ASSOCIATED(this%finite_diff_prm)) THEN
       DEALLOCATE(this%finite_diff_prm, stat=err)
    END IF
    this%integration_step_prm = 1.0_bp
    this%perturbers_prm = .FALSE.
    this%ast_perturbers_prm = .FALSE.

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
    this%generat_multiplier_prm = -1.0_bp
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
    this%sor_burnin_prm = 1
    this%sor_random_obs_prm = .FALSE.
    this%sor_gaussian_pdf_prm = .FALSE.
    this%dchi2_rejection_prm = .TRUE.

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

    ! Parameters for VOMCMC:
    IF (ASSOCIATED(this%vomcmc_map_cmp)) THEN
       DEALLOCATE(this%vomcmc_map_cmp, stat=err)
    END IF
    this%vomcmc_scaling_prm = -1.0_bp
    this%vomcmc_scaling_cmp = -1.0_bp
    this%vomcmc_norb_prm = -1
    this%vomcmc_norb_cmp = -1
    this%vomcmc_ntrial_prm = -1
    this%vomcmc_ntrial_cmp = -1
    this%vomcmc_niter_prm = -1
    this%vomcmc_niter_cmp = -1
    this%vomcmc_nmap_prm = -1
    this%vomcmc_mapping_mask_prm = &
         (/ .TRUE., .FALSE., .FALSE., .FALSE., .FALSE., .FALSE. /)
    this%vomcmc_scaling_ready_cmp = .FALSE.

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
    copy_SO%orb_ml_prm = copy(this%orb_ml_prm)
    copy_SO%cov_type_prm = this%cov_type_prm
    IF (ASSOCIATED(this%res_arr_cmp)) THEN
       ALLOCATE(copy_SO%res_arr_cmp(norb,nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (10).", 1)
          RETURN
       END IF
       copy_SO%res_arr_cmp = this%res_arr_cmp(1:norb,1:nobs,1:6)
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
       copy_SO%rms_arr_cmp = this%rms_arr_cmp(1:norb,1:6)
    END IF
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       ALLOCATE(copy_SO%rchi2_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (30).", 1)
          RETURN
       END IF
       copy_SO%rchi2_arr_cmp = this%rchi2_arr_cmp(1:norb)
    END IF
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       ALLOCATE(copy_SO%pdf_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (35).", 1)
          RETURN
       END IF
       copy_SO%pdf_arr_cmp = this%pdf_arr_cmp(1:norb)
    END IF
    IF (ASSOCIATED(this%reg_apr_arr_cmp)) THEN
       ALLOCATE(copy_SO%reg_apr_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (40).", 1)
          RETURN
       END IF
       copy_SO%reg_apr_arr_cmp = this%reg_apr_arr_cmp(1:norb)
    END IF
    IF (ASSOCIATED(this%jac_arr_cmp)) THEN
       ALLOCATE(copy_SO%jac_arr_cmp(norb,3), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (45).", 1)
          RETURN
       END IF
       copy_SO%jac_arr_cmp = this%jac_arr_cmp(1:norb,1:3)
    END IF
    IF (ASSOCIATED(this%repetition_arr_cmp)) THEN
       ALLOCATE(copy_SO%repetition_arr_cmp(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (46).", 1)
          RETURN
       END IF
       copy_SO%repetition_arr_cmp = this%repetition_arr_cmp(1:norb)
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
    copy_SO%chi2_min_init_prm = this%chi2_min_init_prm
    copy_SO%chi2_min_prm = this%chi2_min_prm
    copy_SO%center_prm = this%center_prm
    copy_SO%dchi2_prm = this%dchi2_prm
    copy_SO%accept_multiplier_prm = this%accept_multiplier_prm
    copy_SO%outlier_rejection_prm = this%outlier_rejection_prm
    copy_SO%regularization_prm = this%regularization_prm
    copy_SO%jacobians_prm = this%jacobians_prm
    copy_SO%multiple_obj_prm = this%multiple_obj_prm
    copy_SO%generat_multiplier_prm = this%generat_multiplier_prm
    copy_SO%generat_gaussian_deviates_prm = this%generat_gaussian_deviates_prm

    ! Bayesian apriori parameters
    copy_SO%apriori_a_max_prm = this%apriori_a_max_prm
    copy_SO%apriori_a_min_prm = this%apriori_a_min_prm
    copy_SO%apriori_periapsis_max_prm = this%apriori_periapsis_max_prm
    copy_SO%apriori_periapsis_min_prm = this%apriori_periapsis_min_prm
    copy_SO%apriori_apoapsis_max_prm = this%apriori_apoapsis_max_prm
    copy_SO%apriori_apoapsis_min_prm = this%apriori_apoapsis_min_prm
    copy_SO%apriori_rho_min_prm = this%apriori_rho_min_prm
    copy_SO%apriori_rho_max_prm = this%apriori_rho_max_prm
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
    copy_SO%ast_perturbers_prm = this%ast_perturbers_prm

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
       ALLOCATE(copy_SO%sor_rho_arr_cmp(SIZE(this%sor_rho_arr_cmp,dim=1),2), stat=err)
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
    copy_SO%sor_burnin_prm = this%sor_burnin_prm
    copy_SO%sor_random_obs_prm = this%sor_random_obs_prm
    copy_SO%sor_gaussian_pdf_prm = this%sor_gaussian_pdf_prm
    copy_SO%dchi2_rejection_prm = this%dchi2_rejection_prm

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

    ! VOMCMC parameters
    IF (ASSOCIATED(this%vomcmc_map_cmp)) THEN
       ALLOCATE(copy_SO%vomcmc_map_cmp(SIZE(this%vomcmc_map_cmp,dim=1), &
            SIZE(this%vomcmc_map_cmp,dim=2)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / copy", &
               "Could not allocate pointer (75).", 1)
          RETURN
       END IF
       copy_SO%vomcmc_map_cmp = this%vomcmc_map_cmp
    END IF
    copy_SO%vomcmc_scaling_prm = this%vomcmc_scaling_prm
    copy_SO%vomcmc_scaling_cmp = this%vomcmc_scaling_cmp
    copy_SO%vomcmc_norb_prm = this%vomcmc_norb_prm
    copy_SO%vomcmc_norb_cmp = this%vomcmc_norb_cmp
    copy_SO%vomcmc_ntrial_prm = this%vomcmc_ntrial_prm
    copy_SO%vomcmc_ntrial_cmp = this%vomcmc_ntrial_cmp
    copy_SO%vomcmc_niter_prm = this%vomcmc_niter_prm
    copy_SO%vomcmc_niter_cmp = this%vomcmc_niter_cmp
    copy_SO%vomcmc_nmap_prm = this%vomcmc_nmap_prm
    copy_SO%vomcmc_mapping_mask_prm = this%vomcmc_mapping_mask_prm
    copy_SO%vomcmc_scaling_ready_cmp = this%vomcmc_scaling_ready_cmp

    ! Least-squares parameters
    copy_SO%ls_corr_fac_prm = this%ls_corr_fac_prm
    copy_SO%ls_niter_major_max_prm = this%ls_niter_major_max_prm
    copy_SO%ls_niter_major_min_prm = this%ls_niter_major_min_prm
    copy_SO%ls_niter_minor_prm = this%ls_niter_minor_prm
    copy_SO%ls_elem_mask_prm = this%ls_elem_mask_prm

    ! Simplex parameters
    copy_SO%smplx_tol_prm = this%smplx_tol_prm
    copy_SO%smplx_niter_prm = this%smplx_niter_prm
    copy_SO%smplx_niter_cmp = this%smplx_niter_cmp
    copy_SO%smplx_force_prm = this%smplx_force_prm
    copy_SO%smplx_similarity_tol_prm = this%smplx_similarity_tol_prm

    ! Observation sampling parameters
    copy_SO%os_norb_prm = this%os_norb_prm
    copy_SO%os_ntrial_prm = this%os_ntrial_prm
    copy_SO%os_sampling_type_prm = this%os_sampling_type_prm

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
  !! Returns error.
  !!
  !! @author  JV, MG
  !! @version 2011-08-22
  !!
  SUBROUTINE autoStatisticalRanging(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this

    REAL(bp), DIMENSION(2,2) :: rho_bounds_
    REAL(bp) :: dchi2, & 
         ddchi2, & ! difference to the given dchi2 limit
         chi2_min_final, chi2_min_
    INTEGER :: niter_, norb_prm 
    LOGICAL :: dchi2_rejection_prm

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       CALL constrainRangeDistributions(this, this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       CALL setRangeBounds(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
    END IF

    this%sor_niter_cmp = 0
    norb_prm = this%sor_norb_prm
    dchi2_rejection_prm = this%dchi2_rejection_prm    

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "***************"
       WRITE(stdout,"(2X,A)") "FIRST ITERATION"
       WRITE(stdout,"(2X,A)") "***************"
    END IF

    ! Force 10 orbits and uniform pdf in the first step
    this%sor_norb_prm = 10
    this%dchi2_rejection_prm = .FALSE.
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

    ! TEST: Return to user-defined pdf mode earlier
    !this%dchi2_rejection_prm = dchi2_rejection_prm

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
    this%dchi2_rejection_prm = dchi2_rejection_prm

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
    DO WHILE (this%sor_rho_histo_cmp > 0 .OR. (this%dchi2_rejection_prm .AND. ddchi2 > 2.0_bp))

       IF (this%sor_niter_cmp >= (this%sor_niter_prm + niter_)) THEN
          EXIT
       END IF

       CALL statisticalRanging(this)
       IF (error .OR. this%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoStatisticalRanging", &
               "Subsequent iteration failed.", 1)
          RETURN
       END IF
       dchi2 = MAXVAL(this%rchi2_arr_cmp) - MINVAL(this%rchi2_arr_cmp)
       ddchi2 = ABS(dchi2 - this%dchi2_prm)

       this%sor_niter_cmp = this%sor_niter_cmp + 1
       IF (this%sor_norb_cmp < this%sor_norb_prm .AND. err_verb >= 2) THEN
          WRITE(stderr,"(A,I0)") " Warning by autoStatisticalRanging:" // &
               " Number of sample orbits too small: ", this%sor_norb_cmp
          WRITE(stderr,"(1X)")
       END IF

       chi2_min_ = this%chi2_min_prm
       rho_bounds_ = this%sor_rho_prm
       CALL updateRanging(this, automatic = .TRUE.)
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,I0,2(A,F10.5))") "Rho histogram flag = ", &
               this%sor_rho_histo_cmp,", dchi2 = ", dchi2, &
               ", d(dchi2) = ", ddchi2
          WRITE(stdout,"(1X)")
       END IF

    END DO
    ! Save the last _used_ rho bounds
    this%sor_rho_prm = rho_bounds_
    this%chi2_min_prm = chi2_min_

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
  LOGICAL FUNCTION containsDiscretePDF(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this

    ! The default result is .FALSE.:
    containsDiscretePDF = .FALSE.

    IF (.NOT. this%is_initialized_prm) THEN
       RETURN
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp) .AND. &
         ASSOCIATED(this%pdf_arr_cmp)) THEN
       containsDiscretePDF = .TRUE.
    END IF

  END FUNCTION containsDiscretePDF





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
         comp_scoords => NULL(), &
         obs_scoords => NULL()
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoords => NULL()
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=12) :: str1, str2
    CHARACTER(len=64) :: &
         frmt = "(F20.15,1X)", &
         efrmt = "(E11.4,1X)"
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         residuals3 => NULL(), &
         partials_arr => NULL(), &
         information_matrix_obs => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: &
         stdev_arr_measur => NULL()
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
         A, &
         cov, &
         eigenvectors, &
         information_matrix_elem, &
         jacobian_matrix, &
         sqrt_eigenvalues
    REAL(bp), DIMENSION(6) :: &
         elements_nominal, &
         ran, &
         deviates, &
         mean, &
         eigenvalues, &
         elements, &
         comp_coord, &
         p, &
         stdev
    REAL(bp) :: &
         sma, &
         apriori, &
         chi2_min_global_ls, &
         chi2, &
         dchi2, &
         pdf_val, &
         obs_, &
         q, &
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
         nfailed, &
         nrotation
    LOGICAL, DIMENSION(:,:), POINTER :: &
         mask_arr => NULL(), &
         mask_measur => NULL()

    iorb = 0
    itrial = 0
    mean = 0.0_bp

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
    stdev_arr_measur => getStandardDeviations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (31)", 1)
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
         failed_flag(11), reg_apriori_arr(this%cos_norb_prm), &
         pdf_arr(this%cos_norb_prm), rchi2_arr(this%cos_norb_prm), &
         jacobian_arr(this%cos_norb_prm,3), obs_coords(nobs,6), &
         mask_arr(nobs,6), mask_measur(nobs,6), stat=err)
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
    failed_flag = 0

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
    ! Residuals for the input orbit:
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
          WRITE(stdout,"(2X,A,1X,A,3"//TRIM(frmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               "observed pos.", obs_coords(i,1:3)
          WRITE(stdout,"(2X,A,1X,A,3"//TRIM(frmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               "computed pos.", comp_coord(1:3)
       END IF
    END DO

    ! Observation mask
    mask_measur = this%obs_masks_prm
    ! Outlier rejection
    IF (this%outlier_rejection_prm) THEN
       DO i=1,nobs
          IF (ANY(ABS(residuals2(i,2:3)) > &
               this%outlier_multiplier_prm*stdev_arr_measur(i,2:3))) THEN
             mask_measur(i,:) = .FALSE.
          END IF
       END DO
    END IF

    ! Reference chi2 from the input orbit which is assumed to be the
    ! maximum-likelihood orbit:
    chi2 = chi_square(residuals2, information_matrix_obs, mask_measur, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / covarianceSampling", &
            "TRACE BACK (60)" // TRIM(errstr),1)
       RETURN
    END IF
    DEALLOCATE(residuals2)
    IF (chi2 < 0.0_bp) THEN
       WRITE(stdout,"(2X,A,1X,A,F10.5,A)") &
            "StochasticOrbit / covarianceSampling:", &
            "Negative chi2 (", chi2, ") for ML solution."
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
    chi2_min_global_ls = chi2
    ! If initial estimate for chi2 corresponding to maximum likelihood
    ! solution remains undetermined (indicated by a negative value),
    ! set it equal to the chi2 computed for the global least squares:
    IF (this%chi2_min_prm < 0.0_bp) THEN
       this%chi2_min_prm = chi2_min_global_ls
    END IF

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,1X,A)") &
            "StochasticOrbit / covarianceSampling:", &
            "Orbital-element covariance matrix:"
       CALL matrix_print(cov,stdout,errstr)
       WRITE(stdout,"(2X,2(A,1X),E12.4)") &
            "StochasticOrbit / covarianceSampling:", &
            "Condition number for orbital-element covariance matrix:", cond_nr(cov, errstr)
    END IF

    A = 0.0_bp
    ! Code needs only the upper diagonal matrix as input. Lower
    ! diagonal matrix will be populated with the Cholesky
    ! decomposition and p will have the diagonal elements.
    DO i=1,6
       A(i,i:6) = cov(i,i:6)
    END DO
    CALL cholesky_decomposition(A, p, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / covarianceSampling:", &
            "Cholesky decomposition unsuccessful:", 1)
       WRITE(stderr,"(A)") TRIM(errstr)
       RETURN
    END IF
    DO i=1,6
       ! Get rid of the upper diagonal matrix that contains the
       ! covariance matrix
       A(i,i:6) = 0.0_bp
       ! Insert diagonal elements to the Cholesky decomposition
       A(i,i) = p(i)
    END DO

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,1X,A)") &
            "StochasticOrbit / covarianceSampling:", &
            "A matrix:"
       CALL matrix_print(A,stdout,errstr)
       IF (this%element_type_prm == "keplerian") THEN
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
               "StochasticOrbit / covarianceSampling:", &
               "Element uncertainty:", stdev(1:2), stdev(3:6)/rad_deg
       ELSE
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
               "StochasticOrbit / covarianceSampling:", &
               "Element uncertainty:", stdev(1:6)
       END IF
       DO i=1,6
          ran = 0.0_bp
          ran(i) = -1.0_bp
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") "StochasticOrbit / covarianceSampling:", &
               "1-sigma deviates", this%cos_nsigma_prm*MATMUL(A,ran) + mean
          ran(i) = 1.0_bp
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") "StochasticOrbit / covarianceSampling:", &
               "1-sigma deviates", this%cos_nsigma_prm*MATMUL(A,ran) + mean
       END DO
       IF (this%element_type_prm == "keplerian") THEN
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
               "StochasticOrbit / covarianceSampling:", &
               "Nominal elements:", elements_nominal(1:2), elements_nominal(3:6)/rad_deg
       ELSE
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
               "StochasticOrbit / covarianceSampling:", &
               "Nominal elements:", elements_nominal(1:6)
       END IF
       WRITE(stdout,"(2X,A,1X,A)") &
            "StochasticOrbit / covarianceSampling:", &
            "Starting sampling..."
       WRITE(stdout,*)
    END IF

    DO WHILE (iorb < this%cos_norb_prm .AND. itrial < this%cos_ntrial_prm)

       itrial = itrial + 1

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A,I0)") &
               "StochasticOrbit / covarianceSampling:", &
               "Trial orbit #", itrial
          WRITE(stdout,"(2X,A,1X,A,I0)") &
               "StochasticOrbit / covarianceSampling:", &
               "Accepted orbits so far: ", iorb
          WRITE(stdout,"(2X,A,1X,2A)") &
               "StochasticOrbit / covarianceSampling:", &
               "Element type: ", TRIM(this%element_type_prm)
          IF (this%element_type_prm == "keplerian") THEN
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Nominal elements:", elements_nominal(1:2), elements_nominal(3:6)/rad_deg
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Element uncertainty:", stdev(1:2), stdev(3:6)/rad_deg
          ELSE
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Nominal elements:", elements_nominal(1:6)
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Element uncertainty:", stdev(1:6)
          END IF
       END IF

       IF (itrial == 1) THEN
          ! Use the nominal orbit as first trial orbit.
          IF (info_verb >= 2) THEN
             deviates = 0.0_bp
             elements = elements_nominal
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Deviates:", deviates
          END IF
          CALL NULLIFY(orb)
          orb = copy(orb_nominal)
       ELSE
          IF (this%cos_gaussian_prm) THEN
             CALL randomGaussian(ran)
          ELSE
             CALL randomNumber(ran)
             ran = 2.0_bp*ran - 1.0_bp
          END IF
          WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
               "StochasticOrbit / covarianceSampling:", &
               "Random numbers:", ran
          deviates = this%cos_nsigma_prm*MATMUL(A,ran) + mean
          ! New coordinates = old coordinates + deviates:
          elements = elements_nominal + deviates
          IF (this%element_type_prm == "keplerian") THEN
             ! 0 <= angle < 2pi :
             elements(4:6) = MODULO(elements(4:6), two_pi)
          END IF
          IF (info_verb >= 2) THEN
             IF (this%element_type_prm == "keplerian") THEN
                WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                     "StochasticOrbit / covarianceSampling:", &
                     "Deviates:", deviates(1:2), deviates(3:6)/rad_deg
             ELSE
                WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                     "StochasticOrbit / covarianceSampling:", &
                     "Deviates:", deviates(1:6)          
             END IF
          END IF
          CALL NULLIFY(orb)
          CALL NEW(orb, elements, this%element_type_prm, frame, t0)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "TRACE BACK (70)", 1)
             RETURN
          END IF
       END IF
       IF (info_verb >= 2) THEN
          IF (this%element_type_prm == "keplerian") THEN
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Generated elements:", elements(1:2), elements(3:6)/rad_deg
          ELSE
             WRITE(stdout,"(2X,A,1X,A19,6(1X,F15.10))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Generated elements:", elements(1:6)          
          END IF
       END IF

       ! Check whether there are a priori requirements on the orbits
       IF (this%informative_apriori_prm) THEN
          elements = getElements(orb,"cometary")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "TRACE BACK (71)", 1)
             RETURN
          END IF
          ! Semimajor axis:
          IF (this%apriori_a_min_prm >= 0.0_bp .OR. &
               this%apriori_a_max_prm >= 0.0_bp) THEN
             sma = elements(1)/(1.0_bp-elements(2))
             IF (this%apriori_a_min_prm >= 0.0_bp .AND. &
                  sma < this%apriori_a_min_prm) THEN
                ! Semimajor axis too small
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (semimajor axis too small: ", sma, " au)"
                END IF
                failed_flag(1) = failed_flag(1) + 1
                CYCLE
             END IF
             IF (this%apriori_a_max_prm >= 0.0_bp .AND. &
                  sma > this%apriori_a_max_prm) THEN
                ! Semimajor axis too large
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (semimajor axis too large: ", sma, " au)"
                END IF
                failed_flag(2) = failed_flag(2) + 1
                CYCLE
             END IF
          END IF
          ! Eccentricity:
          IF (elements(2) < 0.0_bp) THEN
             ! Eccentricity too small.
             failed_flag(3) = failed_flag(3) + 1
             CYCLE
          END IF
          ! Inclination:
          IF (elements(3) < 0.0_bp) THEN
             ! Inclination not defined.
             failed_flag(5) = failed_flag(5) + 1
             CYCLE
          END IF
          IF (elements(3) > pi) THEN
             ! Inclination not defined.
             failed_flag(6) = failed_flag(6) + 1
             CYCLE
          END IF
          ! Periapsis distance:
          IF (this%apriori_periapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_periapsis_max_prm >= 0.0_bp) THEN
             ! Periapsis distance too small:
             IF (this%apriori_periapsis_min_prm >= 0.0_bp .AND. &
                  elements(1) < this%apriori_periapsis_min_prm) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (periapsis distance too small: ", q, " au)"
                END IF
                failed_flag(7) = failed_flag(7) + 1
                CYCLE
             END IF
             ! Periapsis distance too large:
             IF (this%apriori_periapsis_max_prm >= 0.0_bp .AND. &
                  elements(1) > this%apriori_periapsis_max_prm) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (periapsis distance too large: ", q, " au)"
                END IF
                failed_flag(8) = failed_flag(8) + 1
                CYCLE
             END IF
          END IF
          ! Apoapsis distance:
          IF (this%apriori_apoapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_apoapsis_max_prm >= 0.0_bp) THEN
             Q = sma*(1.0_bp+elements(2))
             ! Apoapsis distance too small:
             IF (this%apriori_apoapsis_min_prm >= 0.0_bp .AND. &
                  Q < this%apriori_apoapsis_min_prm) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (apoapsis distance too small: ", Q, " au)"
                END IF
                failed_flag(9) = failed_flag(9) + 1
                CYCLE
             END IF
             ! Apoapsis distance too large:
             IF (this%apriori_apoapsis_max_prm >= 0.0_bp .AND. &
                  Q > this%apriori_apoapsis_max_prm) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (apoapsis distance too large: ", Q, " au)"
                END IF
                failed_flag(10) = failed_flag(10) + 1
                CYCLE
             END IF
          END IF
       END IF

       CALL setParameters(orb, dyn_model=this%dyn_model_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (75)",1)
          RETURN
       END IF

       !!
       !! 5) ACCEPTANCE / REJECTION OF GENERATED ORBIT
       !!
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,1X,A)") &
               "StochasticOrbit / covarianceSampling:", &
               "Start computing ephemerides..."
       END IF
       CALL getEphemerides(orb, obsy_ccoords, comp_scoords, &
            partials_arr=partials_arr)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (80)",1)
          RETURN
       END IF
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,1X,A)") &
               "StochasticOrbit / covarianceSampling:", &
               "Ephemerides ready..."
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
             WRITE(stdout,"(2X,A,1X,A21,3"//TRIM(frmt)//")") &
                  "StochasticOrbit / covarianceSampling:", &
                  "observed pos.", obs_coords(i,1:3)
             WRITE(stdout,"(2X,A,1X,A21,3"//TRIM(frmt)//")") &
                  "StochasticOrbit / covarianceSampling:", &
                  "computed pos.", comp_coord(1:3)
             WRITE(stdout,"(2X,A,1X,A21,3"//TRIM(frmt)//")") &
                  "StochasticOrbit / covarianceSampling:", &
                  "residuals [au,arcsec]", residuals3(iorb+1,i,1), &
                  residuals3(iorb+1,i,2:3)/rad_asec
          END IF
       END DO

!!$       mask_arr = .FALSE.
!!$       WHERE (this%obs_masks_prm .AND. ABS(residuals3(iorb+1,:,:)) > this%res_accept_prm)
!!$          mask_arr = .TRUE.
!!$       END WHERE
!!$       IF (info_verb >= 4) THEN
!!$          DO i=1,nobs
!!$             WRITE(stdout,"(2X,A,1X,A,2"//TRIM(frmt)//")") "O-C residuals (RA, Dec):", &
!!$                  residuals3(iorb+1,i,2:3)/rad_asec
!!$          END DO
!!$          WRITE(stdout,"(2X,A,1X,A,I0,A,I0)") &
!!$               "No of omitted obs/included obs: ", &
!!$               COUNT(mask_arr),"/",n0(2)
!!$       END IF
!!$       IF (COUNT(mask_arr) > 0) THEN
!!$          ! Residuals are too large for at least one observation.
!!$          failed_flag(3) = failed_flag(3) + 1
!!$          IF (info_verb >= 5) THEN
!!$             WRITE(stdout,"(2X,A,1X,A)") &
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
       chi2 = chi_square(residuals3(iorb+1,:,:), information_matrix_obs, mask_measur, errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / covarianceSampling", &
               "TRACE BACK (90)" // TRIM(errstr),1)
          RETURN
       END IF
       IF (chi2 < 0.0_bp) THEN
          WRITE(stdout,"(2X,A,1X,A,F10.5,A,I0)") &
               "StochasticOrbit / covarianceSampling:", &
               "Negative chi2 (", chi2, ") at trial ", itrial
       END IF
       dchi2 = chi2 - this%chi2_min_prm
       IF (this%dchi2_rejection_prm .AND. &
            dchi2 > this%dchi2_prm) THEN
          ! chi2 filtering is used and dchi2 is too large.
          failed_flag(11) = failed_flag(11) + 1
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(2X,A,1X,A,1X,E10.5)") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Failed (dchi2 too large)", dchi2
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

       ! This trial orbit fulfills all criteria so add it to the pile
       ! of (accepted) sample orbits
       iorb = iorb + 1
       orb_arr(iorb) = copy(orb)
       reg_apriori_arr(iorb) = apriori
       pdf_arr(iorb) = apriori*EXP(-0.5_bp*(chi2/COUNT(mask_measur)))
       rchi2_arr(iorb) = chi2/(SUM(n0(1:6)) - 6)
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
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
                  "Unsuccessful computation of determinant of orbital element " // &
                  "jacobian matrix " // TRIM(errstr), 1)
             errstr = ""
             IF (err_verb >= 1) THEN
                WRITE(stderr,"(A)") "cosdec0(1):"
                WRITE(stderr,"(F15.10)") cosdec0(1)
                WRITE(stderr,"(A)") "partials_arr(1:3,:,1):"
                CALL matrix_print(partials_arr(1:3,:,1), stderr, errstr)
                WRITE(stderr,"(A)") "cosdec0(nobs):"
                WRITE(stderr,"(F15.10)") cosdec0(nobs)
                WRITE(stderr,"(A)") "partials_arr(1:3,:,nobs):"
                CALL matrix_print(partials_arr(1:3,:,nobs), stderr, errstr)
                WRITE(stderr,"(A)") "jacobian_matrix:"
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
             CALL errorMessage("StochasticOrbit / covarianceSampling", &
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
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(2X,A,1X,A,F25.15)") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Chi2 new: ", chi2
             WRITE(stdout,"(2X,A,1X,A,3(1X,F10.5))") &
                  "StochasticOrbit / covarianceSampling:", &
                  "dchi2:", chi2 - this%chi2_min_prm
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

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A,I0,1X,A)") &
               "StochasticOrbit / covarianceSampling:", &
               "Orbit #", iorb, "accepted."
          IF (this%regularization_prm) THEN
             WRITE(stdout,"(2X,A,1X,A)") &
                  "StochasticOrbit / covarianceSampling:", &
                  "Sample information matrix:"
             CALL matrix_print(information_matrix_elem,stdout,errstr)
          END IF
          WRITE(stdout,"(2X,A,1X,A,1X,2"//TRIM(efrmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               "Sample chi2, rchi2:", chi2, chi2/(SUM(n0(1:6)) - 6)
          WRITE(stdout,"(2X,A,1X,A,1X,1"//TRIM(efrmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               "Sample apriori:", apriori
          WRITE(stdout,"(2X,A,1X,A,1X,1"//TRIM(efrmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               "Sample pdf:", pdf_arr(iorb)
          WRITE(stdout,"(2X,A,1X,A,1X,1"//TRIM(efrmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               "dchi2:", dchi2
          WRITE(stdout,"(2X,A,1X,A,L1)") &
               "StochasticOrbit / covarianceSampling:", &
               "dchi2 filtering: ", this%dchi2_rejection_prm 
          WRITE(stdout,*)
       END IF

    END DO

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,1X,A)") &
            "StochasticOrbit / covarianceSampling:", &
            "Final number of orbits and the required trials:"
       WRITE(stdout,"(2X,A,1X,2(I0,2X))") &
            "StochasticOrbit / covarianceSampling:", &
            iorb, itrial
       WRITE(stdout,"(2X,A,1X,A)") &
            "StochasticOrbit / covarianceSampling:", &
            "Total failure percentage (1), and failure due to " // &
            "(2) min a, " // &
            "(3) max a, " // &
            "(4) min e, " // &
            "(5) max e, " // &
            "(6) min i, " // &
            "(7) max i, " // &
            "(8) min q, " // &
            "(9) max q, " // &
            "(10) min Q, " // &
            "(11) max Q, " // &
            "(12) pdf:"
       nfailed = SUM(failed_flag)
       nfailed = MAX(1,nfailed)
       WRITE(stdout,"(2X,A,1X,1"//TRIM(frmt)//")") &
            "StochasticOrbit / covarianceSampling:", &
            100.0_bp*REAL(SUM(failed_flag),bp)/itrial
       DO i=1,SIZE(failed_flag)
          WRITE(stdout,"(2X,A,1X,1"//TRIM(frmt)//")") &
               "StochasticOrbit / covarianceSampling:", &
               100.0_bp*REAL(failed_flag(i),bp)/nfailed
       END DO
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
         information_matrix_obs, jacobian_arr, residuals3, &
         stdev_arr_measur, mask_measur, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(jacobian_arr, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(stdev_arr_measur, stat=err)
       DEALLOCATE(mask_measur, stat=err)
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
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr => NULL()
    REAL(bp) :: jac
    INTEGER :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getApoapsisDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (containsDiscretePDF(this)) THEN

       pdf_arr => getDiscretePDF(this, "keplerian")
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

    REAL(bp), DIMENSION(:), POINTER     :: apriori_pdf => NULL(), &
         pdf => NULL()
    INTEGER                             :: norb, err, indx 

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    norb = getNrOfSampleOrbits(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
            "TRACE BACK (5).", 1)
       RETURN
    END IF

    pdf => getDiscretePDF(this)
    pdf = pdf / SUM(pdf)
    IF (PRESENT(apriori)) THEN
       CALL getAPrioriWeights(this, apriori, apriori_pdf)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getBestFittingSampleOrbit", &
               "TRACE BACK (10).", 1)
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

    REAL(bp), DIMENSION(:,:,:), POINTER :: information_matrix => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: residuals => NULL()
    INTEGER :: err, i

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
  REAL(bp) FUNCTION getChi2_matrix(residuals, information_matrix, mask)

    IMPLICIT NONE
    REAL(bp), DIMENSION(:,:), INTENT(in)          :: residuals ! (1:nobs,1:nmulti)
    REAL(bp), DIMENSION(:,:), INTENT(in)          :: information_matrix ! (1:nobs*nmulti,1:nobs*nmulti)
    LOGICAL, DIMENSION(:,:), INTENT(in), OPTIONAL :: mask ! (1:nobs,1:nmulti)

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
       lt_corr, cov_arr, pdf_arr, this_lt_corr_arr)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)                 :: this
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: observers
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER  :: ephemerides_arr
    LOGICAL, INTENT(in), OPTIONAL                         :: lt_corr
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL         :: cov_arr
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL           :: pdf_arr
    TYPE (Orbit), DIMENSION(:,:), POINTER, OPTIONAL       :: this_lt_corr_arr

    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: ephemerides_arr_ => NULL()
    TYPE (Orbit), DIMENSION(:), POINTER                :: this_lt_corr_arr_ => NULL()
    REAL(bp), DIMENSION(:,:,:,:), POINTER              :: partials_arr4 => NULL()
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE            :: coord_arr
    REAL(bp), DIMENSION(:,:,:), POINTER                :: partials_arr3 => NULL()
    REAL(bp), DIMENSION(:,:), POINTER                  :: pdf_arr_2 => NULL()
    REAL(bp), DIMENSION(:), POINTER                    :: pdf_arr_1 => NULL()
    REAL(bp), DIMENSION(6,6)                           :: cov_elm
    REAL(bp), DIMENSION(6)                             :: mean
    REAL(bp)                                           :: det, w
    INTEGER                                            :: i, j, k, err, norb, nobs
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

    IF (containsDiscretePDF(this)) THEN
       ! Discrete orbital-element pdf
       pdf_arr_1 => getDiscretePDF(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "TRACE BACK (3)", 1)
          DEALLOCATE(pdf_arr_1, stat=err)
          RETURN          
       END IF
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
          DEALLOCATE(pdf_arr_1, stat=err)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(partials_arr4, stat=err)
          RETURN
       END IF
       norb = SIZE(this%orb_arr_cmp, dim=1)
       ALLOCATE(pdf_arr_2(norb,nobs), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not allocate memory (5).", 1)
          DEALLOCATE(pdf_arr_1, stat=err)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(partials_arr4, stat=err)
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
                   CALL matrix_print(partials_arr4(i,:,:,j), stderr, errstr)
                END IF
                errstr = ""
                DEALLOCATE(pdf_arr_1, stat=err)
                DEALLOCATE(pdf_arr_2, stat=err)
                DEALLOCATE(ephemerides_arr, stat=err)
                DEALLOCATE(partials_arr4, stat=err)
                RETURN
             END IF
             ! Changed 2008-12-14
             !pdf_arr_2(i,j) = pdf_arr_1(i)/ABS(det)
             pdf_arr_2(i,j) = pdf_arr_1(i) * ABS(det)
          END DO
       END DO
       ! Normalize resulting pdf for each observatory and date
       DO i=1,nobs
          pdf_arr_2(:,i) = pdf_arr_2(:,i)/SUM(pdf_arr_2(:,i))
       END DO
       IF (PRESENT(cov_arr)) THEN
          ! estimate covariance matrices based on the discrete pdf
          norb = SIZE(ephemerides_arr,dim=1)
          ALLOCATE(coord_arr(norb,nobs,6), cov_arr(6,6,nobs), stat=err)
          DO i=1,norb
             DO j=1,nobs
                coord_arr(i,j,:) = getCoordinates(ephemerides_arr(i,j))
             END DO
          END DO
          cov_arr = 0.0_bp
          DO i=1,nobs
             w = SUM(pdf_arr_2(:,i)**2)
             DO j=1,6
                mean(j) = SUM(pdf_arr_2(:,i)*coord_arr(:,i,j))
             END DO
             DO j=1,6
                DO k=j,6
                   cov_arr(j,k,i) = (1.0_bp/(1.0_bp-w)) * &
                        SUM(pdf_arr_2(:,i) * &
                        (coord_arr(:,i,j)-mean(j)) * &
                        (coord_arr(:,i,k)-mean(k)))
                   cov_arr(k,j,i) = cov_arr(j,k,i)
                END DO
             END DO
          END DO
          DEALLOCATE(coord_arr, stat=err)
       END IF
       IF (PRESENT(pdf_arr)) THEN
          ALLOCATE(pdf_arr(norb,nobs), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getEphemerides", &
                  "Could not allocate memory (5).", 1)
             DEALLOCATE(pdf_arr_1, stat=err)
             DEALLOCATE(pdf_arr_2, stat=err)
             DEALLOCATE(ephemerides_arr, stat=err)
             DEALLOCATE(partials_arr4, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             RETURN
          END IF
          pdf_arr = pdf_arr_2
       END IF
       DEALLOCATE(partials_arr4, pdf_arr_1, pdf_arr_2, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Could not deallocate memory (5).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          RETURN
       END IF
    ELSE IF (exist(this%orb_ml_cmp) .AND. PRESENT(cov_arr)) THEN
       ! Least squares solution + covariance matrix:
       IF (getElementType(this%orb_ml_cmp) /= &
            this%element_type_prm) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Element types for ML Orbit and StochasticOrbit are not compatible.", 1)
          RETURN
       END IF
       IF (getElementType(this%orb_ml_cmp) /= &
            this%cov_type_prm) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemerides", &
               "Element types for ML orbit and covariance are not compatible.", 1)
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
            "pdf_arr or cov_arr must be specified.", 1)
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

    REAL(bp), DIMENSION(:,:,:), POINTER :: partials_arr => NULL(), &
         jacobian_lt_corr_arr => NULL(), &
         jacobian_prop_arr => NULL()
    REAL(bp), DIMENSION(:), POINTER     :: pdf_arr_ => NULL()
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
       pdf_arr_ => getDiscretePDF(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(pdf_arr_, stat=err)
          RETURN          
       END IF
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
               "TRACE BACK (10)", 1)
          DEALLOCATE(pdf_arr_, stat=err)
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
          DEALLOCATE(pdf_arr_, stat=err)
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
             DEALLOCATE(pdf_arr_, stat=err)
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
          !pdf_arr(i) = pdf_arr_(i)/ABS(det)
          pdf_arr(i) = pdf_arr_(i) * ABS(det)
       END DO
       DEALLOCATE(partials_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getEphemeris", &
               "Could not deallocate memory.", 1)
          DEALLOCATE(pdf_arr_, stat=err)
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
             DEALLOCATE(pdf_arr_, stat=err)
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
                DEALLOCATE(pdf_arr_, stat=err)
                DEALLOCATE(ephemeris_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr, stat=err)
                DEALLOCATE(jacobian_lt_corr_arr, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                RETURN
             END IF
             ! Changed 2008-12-14
             !pdf_lt_corr_arr(i) = pdf_arr_(i)/ABS(det)
             pdf_lt_corr_arr(i) = pdf_arr_(i) * ABS(det)
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
       DEALLOCATE(pdf_arr_, stat=err)
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
    REAL(bp), DIMENSION(:), POINTER :: pdf => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: periapsis, apoapsis
    LOGICAL, DIMENSION(:), POINTER :: mask_array => NULL(), &
         mask_array_tot => NULL()
    INTEGER :: norb, err, i, ngroup

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ngroup = 18
    IF (.NOT. ASSOCIATED(this%orb_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "No sample orbits available.", 1)
       RETURN
    END IF
    norb = SIZE(this%orb_arr_cmp,dim=1)
    ALLOCATE(weights(ngroup), groups(ngroup), elements(norb,6), &
         apoapsis(norb), periapsis(norb), mask_array(norb), &
         mask_array_tot(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "Could not allocate memory.", 1)
       RETURN
    END IF
    mask_array_tot = .TRUE. ! which are left?
    pdf => getDiscretePDF(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getGroupWeights", &
            "TRACE BACK", 1)          
       RETURN
    END IF
    pdf = pdf / SUM(pdf)
    DO i=1,norb
       elements(i,:) = getElements(this%orb_arr_cmp(i), "cometary")
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getGroupWeights", &
               "TRACE BACK", 1)          
          RETURN
       END IF
       !       write(stdout,"(6(E15.5,1X))") elements(i,:)
    END DO
    apoapsis = elements(:,1) * (1.0_bp + elements(:,2)) / (1.0_bp - elements(:,2))
    periapsis = elements(:,1)
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

    ! Hyperbolic orbits (e.g., interstellar comets and asteroids)
    groups(18) = "Hyperbolic"
    mask_array = .FALSE.
    WHERE (elements(:,2) > 1.0_bp)
       mask_array = .TRUE.
    END WHERE
    weights(18) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Atira (AKA IEO AKA Apohele)
    groups(1) = "Atira"
    mask_array = .FALSE.
    WHERE (apoapsis <= 0.983_bp .AND. &
         mask_array_tot)
       mask_array = .TRUE.
    END WHERE
    weights(1) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Aten
    groups(2) = "NEO/Aten"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) <= 1.0_bp .AND. &
         apoapsis > 0.983_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(2) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Apollo
    groups(3) = "NEO/Apollo"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) > 1.0_bp .AND. &
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
         periapsis <= 5.0_bp/3.0_bp .AND. & ! 1.6666... au
         mask_array_tot) mask_array = .TRUE.
    weights(5) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Hungaria
    groups(6) = "Hungaria"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 1.8_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 2.0_bp .AND. &
         elements(:,2) >= 0.0_bp .AND. &
         elements(:,2) <= 0.16_bp .AND. &
         elements(:,3) >= 15.0_bp .AND. &
         elements(:,3) <= 35.0_bp .AND. &
         periapsis > 5.0_bp/3.0_bp .AND. & ! 1.6666... au
         mask_array_tot) mask_array = .TRUE.
    weights(6) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Phocaea
    groups(7) = "Phocaea"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 2.26_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 2.50_bp .AND. &
         elements(:,2) >= 0.1_bp .AND. &
         elements(:,2) <= 0.3_bp .AND. &
         elements(:,3) >= 10.0_bp .AND. &
         elements(:,3) <= 30.0_bp .AND. &
         periapsis > 5.0_bp/3.0_bp .AND. & ! 1.6666... au
         mask_array_tot) mask_array = .TRUE.
    weights(7) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Mainbelt
    groups(8) = "MBO"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 2.1_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 3.5_bp .AND. &
         elements(:,2) >= 0.0_bp .AND. &
         elements(:,2) <= 0.35_bp .AND. &
         elements(:,3) >= 0.0_bp .AND. &
         elements(:,3) <= 35.0_bp .AND. &
         periapsis > 5.0_bp/3.0_bp .AND. & ! 1.6666... au
         periapsis <= 4.95_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(8) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Hilda
    groups(9) = "Hilda"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 3.75_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 4.02_bp .AND. &
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
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 5.08_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 5.33_bp .AND. &
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
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 36.0_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 39.0_bp .AND. &
         periapsis > 35.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(12) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Plutino
    groups(13) = "TNO/Plutino"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) > 39.0_bp .AND. &
         periapsis < 40.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(13) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Classical TNO
    groups(14) = "TNO/Classical"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) >= 40.0_bp .AND. &
         elements(:,1)/(1.0_bp-elements(:,2)) <= 48.0_bp .AND. &
         periapsis > 35.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(14) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! TNO Outer belt
    groups(15) = "TNO/Outer belt"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) > 48.0_bp .AND. &
         periapsis > 36.0_bp .AND. &
         mask_array_tot) mask_array = .TRUE.
    weights(15) = SUM(pdf,mask=mask_array)
    WHERE (mask_array) mask_array_tot = .FALSE.

    ! Scattered TNO
    groups(16) = "TNO/Scattered"
    mask_array = .FALSE.
    WHERE (elements(:,1)/(1.0_bp-elements(:,2)) > 49.0_bp .AND. &
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





  CHARACTER(len=DESIGNATION_LEN) FUNCTION getID_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this

    IF (TRIM(this%id_prm) == "") THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getID", &
         "Object does not contain an ID.", 1)
       RETURN
    END IF

    getID_SO = this%id_prm

  END FUNCTION getID_SO





  !! *Description*:
  !!
  !! Computes the mean of residual rms values using the sample orbits.
  !!
  !! Returns error.
  !!
  FUNCTION getResidualMeanRMS(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    REAL(bp), DIMENSION(6)                :: getResidualMeanRMS
    INTEGER                               :: norb, i

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResidualMeanRMS", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. ASSOCIATED(this%rms_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResidualMeanRMS", &
            "RMSs are missing -> make an orbit distribution.", 1)
       RETURN       
    END IF

    norb = SIZE(this%rms_arr_cmp,dim=1)
    getResidualMeanRMS = 0.0_bp
    DO i=1,6
       getResidualMeanRMS(i) = SUM(this%rms_arr_cmp(:,i))/norb
    END DO

  END FUNCTION getResidualMeanRMS





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
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr => NULL()
    REAL(bp) :: jac
    INTEGER :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPeriapsisDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (containsDiscretePDF(this)) THEN

       pdf_arr => getDiscretePDF(this, "keplerian")
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
  !! orbits being a PHA (MOID wrt. the Earth <= 0.05 au) at a given  
  !! moment (determined by the epoch of the sample orbits). 
  !! Computations are done using the 2-body approximation.
  !!
  REAL(bp) FUNCTION getPHAProbability(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this

    TYPE (Orbit)                        :: orbit_earth
    TYPE (Time)                         :: t
    REAL(bp), DIMENSION(:), POINTER     :: pdf => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: moid
    REAL(bp), DIMENSION(:,:), POINTER   :: elem => NULL()
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
    pdf => getDiscretePDF(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getPHAProbability", &
            "TRACE BACK (30)", 1)
       RETURN
    END IF
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
  FUNCTION getNominalOrbit(this, ml_orbit)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Orbit)                       :: getNominalOrbit
    LOGICAL, INTENT(in), OPTIONAL      :: ml_orbit

    REAL(bp), DIMENSION(:), POINTER    :: pdf => NULL()
    LOGICAL                            :: ml_orbit_
    INTEGER                            :: indx_ml, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getNominalOrbit", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(ml_orbit)) THEN
       ml_orbit_ = ml_orbit
    ELSE
       ml_orbit_ = .FALSE.
    ENDIF

    IF (.NOT. exist(this%orb_ml_cmp) .AND. .NOT. ml_orbit_) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getNominalOrbit", &
            "The nominal orbit is not available.", 1)
       RETURN
    ELSE IF (ml_orbit_ .AND. containsDiscretePDF(this)) THEN
       pdf => getDiscretePDF(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getNominalOrbit", &
               "TRACE BACK", 1)
          RETURN
       END IF
       indx_ml = MAXLOC(pdf,dim=1)
       DEALLOCATE(pdf, stat=err)
       getNominalOrbit = copy(this%orb_arr_cmp(indx_ml))
    ELSE
       getNominalOrbit = copy(this%orb_ml_cmp)
    END IF

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
  FUNCTION getObservations_SO(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (Observations) :: getObservations_SO

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getObservations", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    getObservations_SO = copy(this%obss)

  END FUNCTION getObservations_SO





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
       asteroid_perturbers, &
       integration_step, &
       integrator, &
       finite_diff, &
       t_inv, &
       element_type, &
       multiple_objects, &
       outlier_rejection, &
       dchi2_rejection, regularized_pdf, &
       accept_multiplier, &
       outlier_multiplier, &
       res_accept, &
       gaussian_pdf, &
       chi2_min_init_prm, &
       chi2_min_prm, &
       center, &
       dchi2_prm, &
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
       generat_multiplier, &
       sor_deviates, &
       vov_norb, vov_ntrial, vov_norb_iter, vov_ntrial_iter, &
       vov_nmap, vov_niter, vov_scaling, vov_mapping_mask, &
       vomcmc_norb, vomcmc_ntrial, vomcmc_norb_iter, vomcmc_ntrial_iter, &
       vomcmc_nmap, vomcmc_niter, vomcmc_scaling, vomcmc_mapping_mask, &
       ls_correction_factor, ls_niter_major_max, ls_niter_major_min, ls_niter_minor, &
       ls_element_mask, &
       smplx_tol, smplx_niter, smplx_force, smplx_similarity_tol, &
       os_norb, os_ntrial, os_sampling_type, generat_gaussian_deviates)

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
         vov_scaling, &
         vomcmc_scaling
    REAL(bp), DIMENSION(6), INTENT(out), OPTIONAL :: &
         finite_diff
    REAL(bp), INTENT(out), OPTIONAL :: &
         integration_step, &
         accept_multiplier, &
         outlier_multiplier, &
         chi2_min_init_prm, &
         chi2_min_prm, &
         dchi2_prm, &
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
         generat_multiplier, &
         ls_correction_factor, &
         smplx_tol, &
         smplx_similarity_tol
    INTEGER, INTENT(out), OPTIONAL :: &
         center, &
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
         vomcmc_norb, &
         vomcmc_ntrial, &
         vomcmc_norb_iter, &
         vomcmc_ntrial_iter, &
         vomcmc_nmap, &
         vomcmc_niter, &
         ls_niter_major_max, &
         ls_niter_major_min, &
         ls_niter_minor, &
         smplx_niter, &
         os_norb, &
         os_ntrial, &
         os_sampling_type
    LOGICAL, DIMENSION(:), INTENT(out), OPTIONAL :: &
         perturbers, &
         vov_mapping_mask, &
         vomcmc_mapping_mask, &
         ls_element_mask
    LOGICAL, INTENT(out), OPTIONAL :: &
         multiple_objects, &
         regularized_pdf, &
         dchi2_rejection, &
         sor_random_obs_selection, &
         gaussian_pdf, &
         outlier_rejection, &
         smplx_force, &
         generat_gaussian_deviates, &
         asteroid_perturbers
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
    IF (PRESENT(asteroid_perturbers)) THEN
       asteroid_perturbers = this%ast_perturbers_prm
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
    IF (PRESENT(dchi2_rejection)) THEN
       dchi2_rejection = this%dchi2_rejection_prm 
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
    IF (PRESENT(chi2_min_init_prm)) THEN
       chi2_min_init_prm = this%chi2_min_init_prm 
    END IF
    IF (PRESENT(chi2_min_prm)) THEN
       chi2_min_prm = this%chi2_min_prm 
    END IF
    IF (PRESENT(dchi2_prm)) THEN
       dchi2_prm = this%dchi2_prm 
    END IF
    IF (PRESENT(center)) THEN
       center = this%center_prm 
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
    IF (PRESENT(generat_multiplier)) THEN
       generat_multiplier = this%generat_multiplier_prm 
    END IF
    IF (PRESENT(generat_gaussian_deviates)) THEN
       generat_gaussian_deviates = this%generat_gaussian_deviates_prm 
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
    IF (PRESENT(vomcmc_norb)) THEN
       vomcmc_norb = this%vomcmc_norb_prm
    END IF
    IF (PRESENT(vomcmc_ntrial)) THEN
       vomcmc_ntrial = this%vomcmc_ntrial_prm 
    END IF
    IF (PRESENT(vomcmc_norb_iter)) THEN
       vomcmc_norb_iter = this%vomcmc_norb_iter_prm 
    END IF
    IF (PRESENT(vomcmc_ntrial_iter)) THEN
       vomcmc_ntrial_iter = this%vomcmc_ntrial_iter_prm 
    END IF
    IF (PRESENT(vomcmc_nmap)) THEN
       vomcmc_nmap = this%vomcmc_nmap_prm 
    END IF
    IF (PRESENT(vomcmc_niter)) THEN
       vomcmc_niter = this%vomcmc_niter_prm 
    END IF
    IF (PRESENT(vomcmc_scaling)) THEN
       vomcmc_scaling = this%vomcmc_scaling_prm
    END IF
    IF (PRESENT(vomcmc_mapping_mask)) THEN
       vomcmc_mapping_mask = this%vomcmc_mapping_mask_prm 
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
    IF (PRESENT(smplx_tol)) THEN
       smplx_tol = this%smplx_tol_prm 
    END IF
    IF (PRESENT(smplx_niter)) THEN
       smplx_niter = this%smplx_niter_prm 
    END IF
    IF (PRESENT(smplx_force)) THEN
       smplx_force = this%smplx_force_prm 
    END IF
    IF (PRESENT(smplx_similarity_tol)) THEN
       smplx_similarity_tol = this%smplx_similarity_tol_prm
    END IF
    IF (PRESENT(os_norb)) THEN
       os_norb = this%os_norb_prm 
    END IF
    IF (PRESENT(os_ntrial)) THEN
       os_ntrial = this%os_ntrial_prm 
    END IF
    IF (PRESENT(os_sampling_type)) THEN
       os_sampling_type = this%os_sampling_type_prm 
    END IF


  END SUBROUTINE getParameters_SO





  !! *Description*:
  !!
  !! Returns error.
  !!
  FUNCTION getDiscretePDF(this, element_type)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    REAL(bp), DIMENSION(:), POINTER    :: getDiscretePDF
    CHARACTER(len=*), INTENT(in), OPTIONAL :: element_type

    REAL(bp), DIMENSION(:), ALLOCATABLE :: jac, pdf
    INTEGER                             :: norb, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getDiscretePDF", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.containsDiscretePDF(this)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getDiscretePDF", &
            "PDF values do not exist.", 1)
       RETURN
    END IF

    norb = SIZE(this%orb_arr_cmp)
    IF (norb /= SIZE(this%pdf_arr_cmp) .AND. &
         norb /= SIZE(this%repetition_arr_cmp)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getDiscretePDF", &
            "Sizes of arrays do not agree with each other.", 1)
       RETURN       
    END IF

    ALLOCATE(pdf(norb), getDiscretePDF(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getDiscretePDF", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%repetition_arr_cmp) .AND. &
         .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
       pdf = this%repetition_arr_cmp
    ELSE
       pdf = this%pdf_arr_cmp
    END IF

    IF (PRESENT(element_type)) THEN
       ALLOCATE(jac(SIZE(pdf)), stat=err)
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
          CALL errorMessage("StochasticOrbit / getDiscretePDF", &
               "Not yet implemented for cometary elements.", 1)
          RETURN
       END IF
       getDiscretePDF = pdf*jac
       DEALLOCATE(jac, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getDiscretePDF", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
    ELSE
       getDiscretePDF = pdf
    END IF
    ! Normalize sum of pdf to unity
    getDiscretePDF = getDiscretePDF/SUM(getDiscretePDF)

    DEALLOCATE(pdf, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getDiscretePDF", &
            "Could not deallocate memory (10).", 1)
       RETURN
    END IF

  END FUNCTION getDiscretePDF





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
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr => NULL()
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

    IF (.NOT.containsDiscretePDF(this)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Orbital-element pdf missing.", 1)
       RETURN
    END IF

    pdf_arr => getDiscretePDF(this, "cartesian")
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

    REAL(bp), DIMENSION(:), POINTER :: phase_angles_ => NULL()
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

    IF (containsDiscretePDF(this)) THEN
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





  FUNCTION getResiduals_SO_obss(this, obss)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)                   :: this
    TYPE (Observations), INTENT(in)                      :: obss
    REAL(bp), DIMENSION(:,:,:), POINTER                  :: getResiduals_SO_obss

    TYPE (CartesianCoordinates), DIMENSION(:), POINTER   :: obsy_ccoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: computed_scoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER   :: observed_scoords => NULL()
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE              :: computed_coords
    REAL(bp), DIMENSION(:,:), ALLOCATABLE                :: observed_coords
    INTEGER                                              :: err, i, j, nobs, norb

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Object 'obss' has not yet been initialized.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    norb = SIZE(this%orb_arr_cmp)
    ALLOCATE(getResiduals_SO_obss(nobs,norb,6), observed_coords(nobs,6), &
         computed_coords(nobs,norb,6), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    observed_scoords => getObservationSCoords(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

    obsy_ccoords => getObservatoryCCoords(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF

    CALL getEphemerides(this%orb_arr_cmp, obsy_ccoords, computed_scoords)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF

    observed_coords = 0.0_bp
    computed_coords = 0.0_bp
    DO i=1,nobs
       CALL rotateToEquatorial(observed_scoords(i))
       observed_coords(i,:) = getCoordinates(observed_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getResiduals", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
       DO j=1,norb
          CALL rotateToEquatorial(computed_scoords(j,i))
          computed_coords(i,j,:) = getCoordinates(computed_scoords(j,i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / getResiduals", &
                  "TRACE BACK (30)", 1)
             RETURN
          END IF
       END DO
    END DO

    DO j=1,norb
       getResiduals_SO_obss(1:nobs,j,1:6) = observed_coords(1:nobs,1:6) - &
            computed_coords(1:nobs,j,1:6)        
       getResiduals_SO_obss(1:nobs,j,2) = getResiduals_SO_obss(1:nobs,j,2) * &
            COS(observed_coords(1:nobs,3))
       DO i=1,nobs
          IF (ABS(getResiduals_SO_obss(i,j,2)) > pi) THEN
             IF (observed_coords(i,2) < computed_coords(i,j,2)) THEN
                observed_coords(i,2) = observed_coords(i,2) + two_pi
             ELSE
                computed_coords(i,j,2) = computed_coords(i,j,2) + two_pi
             END IF
             getResiduals_SO_obss(i,j,2) = (observed_coords(i,2) - &
                  computed_coords(i,j,2)) * COS(observed_coords(i,3))
          END IF
       END DO
    END DO
    DO i=1,SIZE(observed_scoords)
       CALL NULLIFY(observed_scoords(i))
    END DO
    DO i=1,SIZE(obsy_ccoords)
       CALL NULLIFY(obsy_ccoords(i))
    END DO
    DO i=1,SIZE(computed_scoords,dim=1)
       DO j=1,SIZE(computed_scoords,dim=2)
          CALL NULLIFY(computed_scoords(i,j))
       END DO
    END DO
    DEALLOCATE(observed_coords, computed_coords, observed_scoords, &
         obsy_ccoords, computed_scoords, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION getResiduals_SO_obss





  !! *Description*:
  !!
  !!
  !!
  !! Returns error.
  !!
  FUNCTION getResiduals_SO_orb(this, orb)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)                 :: this
    TYPE (Orbit), INTENT(in)                           :: orb
    REAL(bp), DIMENSION(:,:), POINTER                  :: getResiduals_SO_orb

    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: computed_scoords => NULL(), &
         observed_scoords => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE              :: observed_coords, &
         computed_coords
    INTEGER                                            :: err, i, nobs

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(orb)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "Object 'orb' has not yet been initialized.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / getResiduals", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    ALLOCATE(getResiduals_SO_orb(nobs,6), observed_coords(nobs,6), &
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

    getResiduals_SO_orb(1:nobs,1:6) = observed_coords(1:nobs,1:6) - &
         computed_coords(1:nobs,1:6)        
    getResiduals_SO_orb(1:nobs,2) = getResiduals_SO_orb(1:nobs,2) * &
         COS(observed_coords(1:nobs,3))
    DO i=1,nobs
       IF (ABS(getResiduals_SO_orb(i,2)) > pi) THEN
          getResiduals_SO_orb(i,2) = two_pi - getResiduals_SO_orb(i,2)
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

  END FUNCTION getResiduals_SO_orb

  FUNCTION getResiduals_SO_orb_arr(this_arr, orb_arr) result(residuals)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in), DIMENSION(:)       :: this_arr
    TYPE (Orbit), INTENT(in), DIMENSION(:)                 :: orb_arr
    TYPE (SparseArray)                                     :: residuals
    TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: obsy_ccoords
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER     :: temp_ccoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER   :: computed_scoords, sorted_scoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER     :: observed_scoords => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE                  :: observed_coords, computed_coords
    REAL(bp), DIMENSION(:), ALLOCATABLE                    :: mjd_tt_arr, mjd_tt_arr_sorted, mjd_tt_arr_new
    INTEGER                                                :: err, i,j,k, nstorb
    INTEGER, DIMENSION(:), ALLOCATABLE                     :: nobs_arr,indx_arr

    TYPE (Time)                                            :: t

    nstorb = SIZE(this_arr)
    ALLOCATE(nobs_arr(nstorb))
    ! Initialize amount of objects in sparse array.
    ALLOCATE(residuals%vectors(nstorb))

    DO i=1,nstorb
       nobs_arr(i) = getNrOfObservations(this_arr(i)%obss)
       ! Initialize amount of observations for each object in sparse array + amount of residuals.
       ALLOCATE(residuals%vectors(i)%elements(nobs_arr(i),6))
    END DO
    ALLOCATE(obsy_ccoords(SUM(nobs_arr(:))))

    temp_ccoords => getObservatoryCCoords(this_arr(1)%obss)
    obsy_ccoords(1:nobs_arr(1)) = temp_ccoords(:)
    DEALLOCATE(temp_ccoords)
    NULLIFY(temp_ccoords)

    DO i=2,nstorb
       temp_ccoords => getObservatoryCCoords(this_arr(i)%obss)
       obsy_ccoords(1+SUM(nobs_arr(1:i-1)):SUM(nobs_arr(1:i-1))&
            + nobs_arr(i)) = temp_ccoords(:)
       DEALLOCATE(temp_ccoords)
       NULLIFY(temp_ccoords)
    END DO

    CALL getEphemerides(orb_arr, obsy_ccoords, computed_scoords)

    ! Here we compute the residuals for each object and appropriately store them in the result variable.
    DO i=1,nstorb
       ALLOCATE(observed_coords(nobs_arr(i),6), &
            computed_coords(nobs_arr(i),6), stat=err)
       observed_coords = 0.0_bp
       computed_coords = 0.0_bp
       observed_scoords => getObservationSCoords(this_arr(i)%obss)
       DO j=1,nobs_arr(i)
          observed_coords(j,:) = getCoordinates(observed_scoords(j))
          computed_coords(j,:) = getCoordinates(computed_scoords(i,j+SUM(nobs_arr(1:i-1))))
       END DO
       residuals%vectors(i)%elements(:,1:6) = &
            observed_coords(:,1:6) - computed_coords(:,1:6)
       residuals%vectors(i)%elements(1:nobs_arr(i),2) = &
            residuals%vectors(i)%elements(1:nobs_arr(i),2) * &
            COS(observed_coords(1:nobs_arr(i),3))
       DO j=1,nobs_arr(i)
          IF (ABS(residuals%vectors(i)%elements(j,2)) > pi) THEN
             residuals%vectors(i)%elements(j,2) = two_pi - &
                  residuals%vectors(i)%elements(j,2)
          END IF
       END DO
       DEALLOCATE(observed_scoords,observed_coords,computed_coords)

    END DO

    DEALLOCATE(obsy_ccoords, computed_scoords,nobs_arr)
  END FUNCTION getResiduals_SO_orb_arr




  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE getResults_SO(this,&
       reg_apr_arr, jac_arr, repetition_arr_cmp, &
       sor_norb_cmp, sor_ntrial_cmp,& 
       sor_rho_cmp, sor_niter_cmp, &
       sor_rho_arr_cmp, sor_rho_histo_cmp, &
       vov_norb_cmp, vov_ntrial_cmp, &
       vov_niter_cmp, vov_scaling_cmp, &
       vov_map_cmp, vov_scaling_ready_cmp, &
       vomcmc_norb_cmp, vomcmc_ntrial_cmp, &
       vomcmc_niter_cmp, vomcmc_scaling_cmp, &
       vomcmc_map_cmp, vomcmc_scaling_ready_cmp)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)              :: this
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL     :: &
         jac_arr
    REAL(bp), DIMENSION(:), POINTER, OPTIONAL       :: reg_apr_arr
    INTEGER, DIMENSION(:), POINTER, OPTIONAL       :: repetition_arr_cmp
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
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL     :: vomcmc_map_cmp    
    REAL(bp), DIMENSION(6,2), INTENT(out), OPTIONAL :: vomcmc_scaling_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: vomcmc_norb_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: vomcmc_ntrial_cmp
    INTEGER, INTENT(out), OPTIONAL                  :: vomcmc_niter_cmp
    LOGICAL, DIMENSION(6,2), INTENT(out), OPTIONAL  :: vomcmc_scaling_ready_cmp

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
    IF (PRESENT(repetition_arr_cmp)) THEN
       IF (ASSOCIATED(this%repetition_arr_cmp)) THEN
          ALLOCATE(repetition_arr_cmp(SIZE(this%repetition_arr_cmp)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getResults", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          repetition_arr_cmp = this%repetition_arr_cmp
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getResults", &
               "Repetition array not available.", 1)
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
          ALLOCATE(sor_rho_arr_cmp(SIZE(this%sor_rho_arr_cmp,dim=1),2), stat=err)
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
    IF (PRESENT(vomcmc_norb_cmp)) THEN
       vomcmc_norb_cmp = this%vomcmc_norb_cmp
    END IF
    IF (PRESENT(vomcmc_ntrial_cmp)) THEN
       vomcmc_ntrial_cmp = this%vomcmc_ntrial_cmp
    END IF
    IF (PRESENT(vomcmc_niter_cmp)) THEN
       vomcmc_niter_cmp = this%vomcmc_niter_cmp
    END IF
    IF (PRESENT(vomcmc_scaling_cmp)) THEN
       vomcmc_scaling_cmp = this%vomcmc_scaling_cmp
    END IF
    IF (PRESENT(vomcmc_map_cmp)) THEN
       IF (ASSOCIATED(this%vomcmc_map_cmp)) THEN
          ALLOCATE(vomcmc_map_cmp(this%vomcmc_nmap_prm,12), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / getResults", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          vomcmc_map_cmp = this%vomcmc_map_cmp
       ELSE
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getResults", &
               "VoV map not available.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(vomcmc_scaling_ready_cmp)) THEN
       vomcmc_scaling_ready_cmp = this%vomcmc_scaling_ready_cmp 
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

    REAL(bp), DIMENSION(:,:), POINTER  :: residuals => NULL()
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
    REAL(bp), DIMENSION(:), POINTER :: getRMSValues

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
  FUNCTION getSampleOrbits(this, probability_mass)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in)  :: this
    REAL(bp), INTENT(in), OPTIONAL :: probability_mass
    TYPE (Orbit), DIMENSION(:), POINTER :: getSampleOrbits

    REAL(bp), DIMENSION(:), POINTER :: pdf_arr
    INTEGER, DIMENSION(:), ALLOCATABLE :: indx_arr
    INTEGER :: err, i, j, norb

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

    IF (PRESENT(probability_mass)) THEN
       pdf_arr => getDiscretePDF(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / getSampleOrbits", &
               "TRACE BACK", 1)
          RETURN
       END IF
       ALLOCATE(indx_arr(SIZE(pdf_arr)))
       CALL credible_region(pdf_arr, probability_mass, indx_arr, errstr)
       IF (LEN_TRIM(errstr) > 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getSampleOrbits", &
               "Error in computation of credible region: " // TRIM(errstr), 1)
          RETURN
       END IF
       norb = 0
       DO i=1,SIZE(indx_arr)
          IF (indx_arr(i) > 0) THEN
             norb = norb + 1
          END IF
       END DO
       ALLOCATE(getSampleOrbits(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / getSampleOrbits", &
               "Could not allocate memory.", 1)
          RETURN
       END IF
       j = 0
       DO i=1, SIZE(indx_arr)
          IF (indx_arr(i) > 0) THEN
             j = j + 1
             getSampleOrbits(j) = copy(this%orb_arr_cmp(indx_arr(i)))
          END IF
       END DO
       DEALLOCATE(indx_arr, pdf_arr, stat=err)
    ELSE
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
    END IF

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
    REAL(bp), DIMENSION(:), POINTER :: pdf => NULL()
    CHARACTER(len=ELEMENT_TYPE_LEN)    :: element_type_
    INTEGER                            :: i, err

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getStandardDeviations", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.containsDiscretePDF(this)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / getStandardDeviations", &
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
       pdf => getDiscretePDF(this, element_type)
    ELSE
       pdf => getDiscretePDF(this)
    END IF

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
            errstr=errstr)
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
  !! Automated version of MCMC ranging. Stable version (pre-2015 version, see autoMCMCRanging2)
  !!
  !! Calls MCMCRanging (earlier: MCMCRanging3).
  !!
  !! For initialization using either
  !!
  !!  - no a priori orbit information: Ranging (i.e. uniform proposal)
  !!
  !!  - orbital distribution from input file
  !!
  SUBROUTINE autoMCMCRanging(this)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout) :: this
    TYPE(StochasticOrbit) :: storb
    REAL(bp), DIMENSION(4) :: sor_rho_init

    IF (.NOT.containsDiscretePDF(this)) THEN

       this%sor_niter_cmp = 0
       ! First iteration with Ranging and uniform sampling
       storb = copy(this)
       storb%sor_norb_prm = 20
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm
       storb%dchi2_rejection_prm = .FALSE.

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A)") "***************"
          WRITE(stdout,"(2X,A)") "FIRST ITERATION"
          WRITE(stdout,"(2X,A)") "***************"
       END IF
       CALL statisticalRanging(storb)
       IF (error .OR. storb%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "First iteration failed.", 1)
          RETURN
       END IF
       this%sor_niter_cmp = 1
       CALL updateRanging(storb, automatic=.TRUE.)

       ! Second iteration
       storb%sor_norb_prm = this%sor_norb_sw_prm
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm    
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(1X)")
          WRITE(stdout,"(2X,A)") "*****************"
          WRITE(stdout,"(2X,A)") "SECOND ITERATION "
          WRITE(stdout,"(2X,A)") "*****************"
       END IF
       CALL statisticalRanging(storb)
       IF (error .OR. storb%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "Second iteration failed.", 1)
          RETURN
       END IF
       this%sor_niter_cmp = 2
       CALL updateRanging(storb, automatic=.TRUE.)

    ELSE

       storb = copy(this)
       storb%sor_norb_prm = 50
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm
       ! Compute range distribution from orbit distribution
       ! - this version updates this%orb_ml_cmp to the best-fit orbit 
       ! (based on the chi2 values of the new observations)
       ! - needed because for mcmc, orbital p.d.f = 1
       CALL constrainRangeDistributions2(storb, storb%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "TRACE BACK (15)", 1)
          STOP
       END IF
       CALL setRangeBounds(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "TRACE BACK (20)", 1)
          STOP
       END IF
       !       storb%orb_ml_prm=getNominalOrbit(storb, ml_orbit=.true.)
       ! Optionally, use ML-orbit to start the next chain (TO BE TESTED!)
       IF (this%sor_iterate_bounds_prm(1)) THEN
          storb%orb_ml_prm = copy(storb%orb_ml_cmp)
       END IF

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A)") "***************"
          WRITE(stdout,"(2X,A)") "FIRST ITERATION"
          WRITE(stdout,"(2X,A)") "***************"
       END IF
       CALL MCMCRanging(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "TRACE BACK (10)", 1)
          STOP
       END IF
       this%sor_niter_cmp = 1

       ! Is this needed, since burn in in MCMC discards orbits based on their rms-values? 
       !CALL constrainRangeDistributions(storb, storb%obss)
       !IF (error) THEN
       !   CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
       !        "TRACE BACK (15)", 1)
       !   STOP
       !END IF
       CALL setRangeBounds(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "TRACE BACK (20)", 1)
          STOP
       END IF
       ! Optionally, use ML-orbit to start the next chain (TO BE TESTED!)
       IF (this%sor_iterate_bounds_prm(1)) THEN
          storb%orb_ml_prm = copy(storb%orb_ml_cmp)
       END IF
       !storb%orb_ml_prm = getNominalOrbit(storb, ml_orbit=.true.)

       ! Second iteration
       storb%sor_norb_prm = this%sor_norb_sw_prm
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(1X)")
          WRITE(stdout,"(2X,A)") "*****************"
          WRITE(stdout,"(2X,A)") "SECOND ITERATION "
          WRITE(stdout,"(2X,A)") "*****************"
       END IF
       CALL MCMCRanging(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "TRACE BACK (10)", 1)
          STOP
       END IF
       this%sor_niter_cmp = 2
       CALL setRangeBounds(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
               "TRACE BACK (20)", 1)
          STOP
       END IF
    END IF

    sor_rho_init = getRangeBounds(storb)
    this%sor_rho_prm(1,1:2) = sor_rho_init(1:2)
    this%sor_rho_prm(2,1:2) = sor_rho_init(3:4)
    ! Optionally, use ML-orbit to start the next chain (TO BE TESTED!)
    IF (this%sor_iterate_bounds_prm(1)) THEN
       this%orb_ml_prm = copy(storb%orb_ml_cmp)
    END IF
    CALL NULLIFY(storb)

    ! Start actual Markov chain
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "*****************"
       WRITE(stdout,"(2X,A,1X,I0,1X,A)") "FINAL SAMPLING FOR", this%sor_norb_prm, "SAMPLE ORBITS..."
       WRITE(stdout,"(2X,A)") "*****************"
    END IF
    CALL MCMCRanging(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / autoMCMCRanging", &
            "TRACE BACK (50)", 1)
       STOP
    END IF
    this%sor_niter_cmp = 3

  END SUBROUTINE autoMCMCRanging





  !! *Description*:
  !!
  !! Automated version of random walk. Calls randomWalkRanging.
  !!
  !! For initialization using either
  !!
  !! - no a priori orbit information: Ranging (i.e. uniform proposal)
  !!
  !! - orbital distribution from input file
  !!
  SUBROUTINE autoRandomWalkRanging(this)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout) :: this

    TYPE(StochasticOrbit) :: storb
    REAL(bp), DIMENSION(4) :: sor_rho_init
    LOGICAL, DIMENSION(:,:), POINTER :: &
         obs_masks => NULL()

    IF (.NOT.containsDiscretePDF(this)) THEN

       this%sor_niter_cmp = 0
       ! First iteration with Ranging and uniform sampling
       storb = copy(this)
       storb%sor_norb_prm = 20
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm
       storb%dchi2_rejection_prm = .FALSE.

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A)") "***************"
          WRITE(stdout,"(2X,A)") "FIRST ITERATION"
          WRITE(stdout,"(2X,A)") "***************"
       END IF
       CALL statisticalRanging(storb)
       IF (error .OR. storb%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "First iteration failed.", 1)
          RETURN
       END IF
       this%sor_niter_cmp = 1
       CALL updateRanging(storb, automatic=.TRUE.)

       !CALL NULLIFY(storb)

       ! Second iteration
       !storb = copy(this)
       storb%sor_norb_prm = this%sor_norb_sw_prm
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm    
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(1X)")
          WRITE(stdout,"(2X,A)") "*****************"
          WRITE(stdout,"(2X,A)") "SECOND ITERATION "
          WRITE(stdout,"(2X,A)") "*****************"
       END IF
       CALL statisticalRanging(storb)
       IF (error .OR. storb%sor_norb_cmp <= 1) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "Second iteration failed.", 1)
          RETURN
       END IF
       this%sor_niter_cmp = 2
       CALL updateRanging(storb, automatic=.TRUE.)

    ELSE
       storb = copy(this)
       storb%sor_norb_prm = 50
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm
       obs_masks => getObservationMasks(storb%obss)
       storb%chi2_min_prm = REAL(COUNT(obs_masks), bp)
       ! Compute range distribution from orbit distribution
       CALL constrainRangeDistributions2(storb, storb%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "TRACE BACK (15)", 1)
          STOP
       END IF
       CALL setRangeBounds(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "TRACE BACK (20)", 1)
          STOP
       END IF
       !       storb%orb_ml_prm=getNominalOrbit(storb, ml_orbit=.true.)
       ! Optionally, use ML-orbit to start the next chain (TO BE TESTED!)
       IF (this%sor_iterate_bounds_prm(1)) storb%orb_ml_prm=copy(storb%orb_ml_cmp)

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A)") "***************"
          WRITE(stdout,"(2X,A)") "FIRST ITERATION"
          WRITE(stdout,"(2X,A)") "***************"
       END IF
       CALL randomWalkRanging(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "TRACE BACK (10)", 1)
          STOP
       END IF
       this%sor_niter_cmp = 1

       ! Is this needed, since burn in in MCMC discards orbits based on their rms-values? 
       !CALL constrainRangeDistributions(storb, storb%obss)
       !IF (error) THEN
       !   CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
       !        "TRACE BACK (15)", 1)
       !   STOP
       !END IF
       CALL setRangeBounds(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "TRACE BACK (20)", 1)
          STOP
       END IF
       IF (storb%chi2_min_prm > storb%chi2_min_cmp) storb%chi2_min_prm = storb%chi2_min_cmp
       ! Optionally, use ML-orbit to start the next chain (TO BE TESTED!)
       IF (this%sor_iterate_bounds_prm(1)) storb%orb_ml_prm=copy(storb%orb_ml_cmp)

       !Second iterarion
       storb%sor_norb_prm = this%sor_norb_sw_prm
       storb%sor_ntrial_prm = this%sor_ntrial_sw_prm
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(1X)")
          WRITE(stdout,"(2X,A)") "*****************"
          WRITE(stdout,"(2X,A)") "SECOND ITERATION "
          WRITE(stdout,"(2X,A)") "*****************"
       END IF
       CALL randomWalkRanging(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "TRACE BACK (10)", 1)
          STOP
       END IF
       this%sor_niter_cmp = 2
       CALL setRangeBounds(storb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
               "TRACE BACK (20)", 1)
          STOP
       END IF

    ENDIF
    sor_rho_init=getRangeBounds(storb)
    this%sor_rho_prm(1,1:2) = sor_rho_init(1:2)
    this%sor_rho_prm(2,1:2) = sor_rho_init(3:4)
    ! Optionally, use ML-orbit to start the next chain (TO BE TESTED!)
    IF (this%sor_iterate_bounds_prm(1)) THEN
       this%orb_ml_prm=copy(storb%orb_ml_cmp)
    END IF
    this%chi2_min_prm = storb%chi2_min_cmp
    CALL NULLIFY(storb)

    ! Start actual Markov chain
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "*****************"
       WRITE(stdout,"(2X,A,1X,I0,1X,A)") "FINAL SAMPLING FOR", this%sor_norb_prm, "SAMPLE ORBITS..."
       WRITE(stdout,"(2X,A)") "*****************"
    END IF
    CALL randomWalkRanging(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / autoRandomWalkRanging", &
            "TRACE BACK (50)", 1)
       STOP
    END IF
    this%sor_niter_cmp = 3

  END SUBROUTINE autoRandomWalkRanging





  !! *Description:*
  !!
  !! Random-walk ranging method.
  !!
  !!   - Burn in phase: discard accepted orbits until first chi2 in the phase-space area 
  !!     defined by dchi2.
  !!
  !!   - Random-walk option (ran_walk) can be turned off, which reverts the routine to 
  !!     a version of MCMC ranging
  !!
  SUBROUTINE randomWalkRanging(this)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout) :: this

    TYPE(Orbit) :: orb
    TYPE(SphericalCoordinates), DIMENSION(:), POINTER :: &
         obs_scoords => NULL(), &
         comp_scoords => NULL()
    TYPE(SphericalCoordinates) :: obs_scoord1, obs_scoord2, obs_help
    TYPE(CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords => NULL()
    TYPE(CartesianCoordinates) :: obs_ccoord_helio1, &
         obs_ccoord_helio2, obs_ccoord_topo1, obs_ccoord_topo2
    TYPE(Time) :: tt
    REAL(bp), DIMENSION(:,:,:), POINTER :: information_matrix_obs => NULL(), &
         partials_arr => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: residuals, &
         comp_coords, &
         obs_coords, &
         jacobians
    REAL(bp), DIMENSION(:), ALLOCATABLE :: cosdec0
    REAL(bp), DIMENSION(6,6):: information_matrix_elem, &
         jacobian_matrix
    REAL(bp), DIMENSION(6) :: elements, rans, state, state_, &
         proposal_density
    REAL(bp), DIMENSION(3) :: pos, rms
    REAL(bp) :: chi2, ran, t1, t2, obs_coord, jac_sph_inv, &
         jac_sph_inv_, chi2_, chi2min, avalue, a, q, e, dchi2, chi2min_prm
    INTEGER, DIMENSION(:,:), POINTER :: obs_pair_arr => NULL()
    INTEGER, DIMENSION(6) :: n0, n0_
    INTEGER, DIMENSION(2) :: obs_pair
    INTEGER :: ndof, i, itrial, err, nobs, j, norb_acc
    INTEGER, PARAMETER :: burnin_max = 50    ! Make user-defined parameter? 
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks => NULL()
    LOGICAL :: accepted, first = .TRUE., burn_in = .TRUE., frst_found = .FALSE., burnin_case=.FALSE.

    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Starting MCMC ranging using random walk"
    END IF

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Allocate memory for solution storing - orbits
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DEALLOCATE(this%orb_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(this%orb_arr_cmp(this%sor_norb_prm), stat=err) 
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "Could not allocate memory (5).", 1)
       RETURN
    END IF
    nobs = getNrOfObservations(this%obss)
    ALLOCATE(this%res_arr_cmp(this%sor_norb_prm,nobs,6), &
         this%pdf_arr_cmp(this%sor_norb_prm), &
         residuals(nobs,6), obs_coords(nobs,6), &
         comp_coords(nobs, 6), cosdec0(nobs), &
         this%reg_apr_arr_cmp(this%sor_norb_prm), &
         this%jac_arr_cmp(this%sor_norb_prm,3), &
         this%rchi2_arr_cmp(this%sor_norb_prm), &
         this%rms_arr_cmp(this%sor_norb_prm,6), &
         this%repetition_arr_cmp(this%sor_norb_prm), &
         this%sor_rho_arr_cmp(this%sor_norb_prm,2), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "Could not allocate memory (10).", 1)
       RETURN 
    END IF
    this%repetition_arr_cmp = 0
    ! Allocate memory for solution storing - pdf 
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(this%pdf_arr_cmp(this%sor_norb_prm), stat=err) 
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "Could not allocate memory (15).", 1)
       RETURN
    END IF
    this%sor_rho_cmp(:,1) = HUGE(this%sor_rho_cmp)
    this%sor_rho_cmp(:,2) = -HUGE(this%sor_rho_cmp)

    ! Get observed sky positions (original observational data) 
    obs_scoords => getObservationSCoords(this%obss) 
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    obs_coords = 0.0_bp
    DO i=1,nobs
       CALL rotateToEquatorial(obs_scoords(i))
       obs_coords(i,:) = getCoordinates(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (10)", 1)
          CALL NULLIFY(obs_scoords(i))
          RETURN
       END IF
       ! CALL NULLIFY(obs_scoords(i))
       cosdec0(i) = COS(obs_coords(i,3))
    END DO

    obs_masks => getObservationMasks(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    ! Observation number counter (observation mask must be up-to-date!), 
    ! construct cosine array:
    DO i=1,6
       n0(i) = COUNT(obs_masks(:,i))
    END DO
    ndof = COUNT(obs_masks) - 6

    information_matrix_obs => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "TRACE BACK (20)", 1)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(comp_coords, stat=err)
       RETURN
    END IF

    obsy_ccoords => getObservatoryCCoords(this%obss) 
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF

    ! Use predefined observation pair if not using random selection
    IF (.NOT. this%sor_random_obs_prm) THEN
       obs_pair = RESHAPE(this%sor_pair_arr_prm, (/ 2 /))
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "Use of random observation pairs not accepted.", 1)
       RETURN
    END IF
    ! Copy observation pair
    obs_scoord1 = copy(obs_scoords(obs_pair(1)))
    obs_scoord2 = copy(obs_scoords(obs_pair(2)))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "TRACE BACK (30)", 1)
       RETURN
    END IF

    ! proposal density; assuming 3sigma rho bounds
    proposal_density(1) = (this%sor_rho_prm(1,2) - this%sor_rho_prm(1,1))/6.0_bp
    proposal_density(2) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(1),2,2))
    proposal_density(3) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(1),3,3))
    proposal_density(4) = (this%sor_rho_prm(2,2) - this%sor_rho_prm(2,1))/6.0_bp 
    proposal_density(5) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(2),2,2))
    proposal_density(6) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(2),3,3))
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(A)") "Proposal density for the first observation (y_l [au], RA [as], Dec [as]): "
       WRITE(stdout,"(F15.7)") proposal_density(1)
       WRITE(stdout,"(F15.7)") proposal_density(2)/rad_asec
       WRITE(stdout,"(F15.7)") proposal_density(3)/rad_asec
       WRITE(stdout,"(A)") "Proposal density for the second observation (y_r [au], RA [as], Dec [as]): "
       WRITE(stdout,"(F15.7)") proposal_density(4)
       WRITE(stdout,"(F15.7)") proposal_density(5)/rad_asec
       WRITE(stdout,"(F15.7)") proposal_density(6)/rad_asec
    END IF

    ! rho1, RA1, Dec1, drho12, RA2, Dec2 of the proposed orbit
    state_(1:3) = getPosition(obs_scoord1)
    state_(4:6) = getPosition(obs_scoord2)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "TRACE BACK (50)", 1)
       DEALLOCATE(obs_masks, stat=err)
       RETURN
    END IF

    ! Center the first trial to the ml orbit if available
    IF (exist(this%orb_ml_prm)) THEN
       CALL getEphemerides(this%orb_ml_prm, obsy_ccoords, comp_scoords)
       state_(1)=getDistance(comp_scoords(obs_pair(1)))
       state_(4)=getDistance(comp_scoords(obs_pair(2)))
    ELSE
       state_(1) = (this%sor_rho_prm(1,2) + this%sor_rho_prm(1,1))/2.0_bp
       state_(4) = state_(1) + (this%sor_rho_prm(2,2) + this%sor_rho_prm(2,1))/2.0_bp
    ENDIF

    ! set counter some parameters
    this%sor_norb_cmp = 0
    this%sor_ntrial_cmp = -1 ! first loop is setting the parameters
    chi2min = HUGE(chi2min)
    first = .TRUE.
    itrial = 0
    norb_acc = 0
    burn_in = .TRUE.

    DO WHILE (this%sor_norb_cmp < this%sor_norb_prm .AND. &
         this%sor_ntrial_cmp < this%sor_ntrial_prm)

       ! Get next state(topocentric distances and sky positions for
       ! the two observations) around the previous computed state
       IF (itrial == 0 .AND. info_verb >= 3) THEN
          WRITE(stdout,"(A)") "Generating new state..."
       END IF

       this%sor_ntrial_cmp = this%sor_ntrial_cmp + 1
       itrial = itrial + 1
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(A,I0)") "Trial state #", itrial
       END IF

       IF (first) THEN
          state = state_
       ELSE
          CALL randomGaussian(rans)
          state(1) = state_(1) + rans(1) * proposal_density(1) + rans(4) * proposal_density(4)
          state(2:3) = state_(2:3) + rans(2:3) * proposal_density(2:3)
          state(4) = state_(4) + rans(1) * proposal_density(1) - rans(4) * proposal_density(4)
          state(5:6) = state_(5:6) + rans(5:6) * proposal_density(5:6)
       END IF

       ! Check if the distances are negative
       IF (state(1) <= 0.0_bp .OR. state(4) <= 0.0_bp) THEN
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(2X,A,F10.7,A)") &
                  "Failed (One or both topocentric distances smaller than the Earth radius.)" 
          END IF
          CYCLE
       END IF

       ! Create new spherical coordinates with generated states in equatorial frame
       CALL NULLIFY(obs_scoord1)
       CALL NULLIFY(obs_scoord2)
       CALL NEW(obs_scoord1, state(1), state(2), state(3), getTime(obsy_ccoords(obs_pair(1))))
       CALL NEW(obs_scoord2, state(4), state(5), state(6), getTime(obsy_ccoords(obs_pair(2))))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (65)", 1)
          RETURN
       END IF

       ! Create new topocentric cartesian coordinates of the observations in equatorial frame
       CALL NULLIFY(obs_ccoord_topo1)
       CALL NULLIFY(obs_ccoord_topo2)
       CALL NEW(obs_ccoord_topo1, obs_scoord1) 
       CALL NEW(obs_ccoord_topo2, obs_scoord2) 
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (70)", 1)
          RETURN
       END IF

       ! Rotate to ecliptic
       CALL rotateToEcliptic(obs_ccoord_topo1)
       CALL rotateToEcliptic(obs_ccoord_topo2)
       CALL rotateToEcliptic(obsy_ccoords(obs_pair(1)))
       CALL rotateToEcliptic(obsy_ccoords(obs_pair(2)))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (75)", 1)
          RETURN
       END IF

       ! Cartesian heliocentric coordinates
       CALL NULLIFY(obs_ccoord_helio1)
       CALL NULLIFY(obs_ccoord_helio2)
       obs_ccoord_helio1 = copy(obsy_ccoords(obs_pair(1)) + obs_ccoord_topo1)
       obs_ccoord_helio2 = copy(obsy_ccoords(obs_pair(2)) + obs_ccoord_topo2)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (80)", 1)
          RETURN
       END IF

       CALL estimateLightTime(obs_ccoord_helio1, state(1))! changed from 3 and 6
       CALL estimateLightTime(obs_ccoord_helio2, state(1)+state(4))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (85)", 1)
          RETURN
       END IF

       ! Find orbit candidate at the epoch of the first observation by
       ! using the chosen method to solve the 2-point boundary value
       ! problem:
       CALL NULLIFY(orb)
       CALL NEW(orb, obs_ccoord_helio1, obs_ccoord_helio2, &
            this%sor_2point_method_prm, this%apriori_a_max_prm, &
            center=this%center_prm, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (86)", 1)
          error = .FALSE.
          CYCLE
       END IF

       IF (this%informative_apriori_prm .AND. .NOT. first) THEN
          ! Semimajor axis:
          a = -1.0_bp
          IF (this%apriori_a_min_prm >= 0.0_bp .OR. &
               this%apriori_a_max_prm >= 0.0_bp) THEN
             a = getSemimajorAxis(orb)
             IF (this%apriori_a_min_prm >= 0.0_bp .AND. &
                  a < this%apriori_a_min_prm) THEN
                ! Semimajor axis too small
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (semimajor axis too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             IF (this%apriori_a_max_prm >= 0.0_bp .AND. &
                  a > this%apriori_a_max_prm) THEN
                ! Semimajor axis too large
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (semimajor axis too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
          IF (error) THEN 
             CALL errorMessage("StochasticOrbit / randomWalkRanging", &
                  "TRACE BACK (xx)", 1)
             error = .FALSE.
             CYCLE
          END IF
          ! Periapsis distance:
          IF (this%apriori_periapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_periapsis_max_prm >= 0.0_bp) THEN
             IF (a >= 0.0_bp) THEN
                CALL getPeriapsisDistance(orb, q, a)
             ELSE
                CALL getPeriapsisDistance(orb, q)
             END IF
             IF (error) THEN
                error = .FALSE.
                CYCLE
             END IF
             ! Periapsis distance too small:
             IF (this%apriori_periapsis_min_prm >= 0.0_bp .AND. &
                  q < this%apriori_periapsis_min_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (periapsis distance too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             ! Periapsis distance too large:
             IF (this%apriori_periapsis_max_prm >= 0.0_bp .AND. &
                  q > this%apriori_periapsis_max_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (periapsis distance too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
          ! Apoapsis distance:
          IF (this%apriori_apoapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_apoapsis_max_prm >= 0.0_bp) THEN
             IF (a >= 0.0_bp) THEN
                CALL getApoapsisDistance(orb, Q, a)
             ELSE
                CALL getApoapsisDistance(orb, Q)
             END IF
             ! Apoapsis distance too small:
             IF (this%apriori_apoapsis_min_prm >= 0.0_bp .AND. &
                  Q < this%apriori_apoapsis_min_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (apoapsis distance too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             ! Apoapsis distance too large:
             IF (this%apriori_apoapsis_max_prm >= 0.0_bp .AND. &
                  Q > this%apriori_apoapsis_max_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (apoapsis distance too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
       END IF

       CALL setParameters(orb, &
            dyn_model=this%dyn_model_prm, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)

       ! checking if the epochs are the same
       IF (.NOT. equal(this%t_inv_prm,getTime(orb))) THEN
          CALL propagate(orb, this%t_inv_prm)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / randomWalkRanging", &
                  "TRACE BACK (90)", 1)
             error = .FALSE.
             CYCLE
          END IF
       END IF

       IF (info_verb >= 3) THEN
          WRITE(stdout,"(A,1X,I0,1X,A)") "New state generated using", itrial, "trials."
       END IF
       itrial = 0

       !! 
       !! COMPUTE CHI2 AND PDF FOR PROPOSED ORBIT
       !!

       ! get computed sky positions - spherical coordinates
       CALL getEphemerides(orb, obsy_ccoords, comp_scoords, partials_arr=partials_arr)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (45)", 1)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          !DEALLOCATE(obsy_ccoords, stat=err)
          !DEALLOCATE(obs_masks, stat=err)
          !DEALLOCATE(residuals, stat=err)
          !DEALLOCATE(information_matrix_obs, stat=err)
          !DEALLOCATE(obs_coords, stat=err)         
          !DEALLOCATE(comp_coords, stat=err)
          !RETURN
          error = .FALSE.
          CYCLE
       END IF

       ! Residuals
       comp_coords = 0.0_bp
       DO i=1,nobs
          CALL rotateToEquatorial(comp_scoords(i))
          comp_coords(i,:) = getCoordinates(comp_scoords(i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / randomWalkRanging", &
                  "TRACE BACK (55)", 1)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(obs_masks, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(comp_coords, stat=err)
             RETURN
          END IF
          CALL NULLIFY(comp_scoords(i))
       END DO
       DEALLOCATE(comp_scoords)
       residuals(1:nobs,1:6) = obs_coords(1:nobs,1:6) - &
            comp_coords(1:nobs,1:6)        
       residuals(1:nobs,2) = residuals(1:nobs,2) * &
            COS(obs_coords(1:nobs,3))
       DO i=1,nobs
          IF (ABS(residuals(i,2)) > pi) THEN
             obs_coord = obs_coords(i,2)
             IF (obs_coord < comp_coords(i,2)) THEN
                obs_coord = obs_coord + two_pi
             ELSE
                comp_coords(i,2) = comp_coords(i,2) + two_pi
             END IF
             residuals(i,2) = (obs_coord - &
                  comp_coords(i,2)) * COS(obs_coords(i,3))
          END IF
       END DO

       chi2 = chi_square(residuals, information_matrix_obs, errstr=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "TRACE BACK (60)", 1)
          DEALLOCATE(obs_masks, stat=err)
          RETURN
       END IF

       ! Jacobians
       ! Determinant of Jacobian between topocentric
       ! coordinates (inverse problem coordinates)
       ! and orbital parameters required for output
       ! ("Topocentric Wrt Cartesian/Keplerian"):
       jacobian_matrix(1:3,:) = partials_arr(1:3,:,obs_pair(1)) / &
            cosdec0(obs_pair(1))
       jacobian_matrix(4:6,:) = partials_arr(1:3,:,obs_pair(2)) / &
            cosdec0(obs_pair(2))
       DEALLOCATE(partials_arr, stat=err)
       jac_sph_inv = ABS(determinant(jacobian_matrix, errstr))
       IF (LEN_TRIM(errstr) /= 0) THEN
          CALL errorMessage("StochasticOrbit / randomWalkRanging", &
               "Unsuccessful computation of determinant of orbital element " // &
               "jacobian matrix: " // TRIM(errstr), 1)
          CALL matrix_print(jacobian_matrix, stderr, errstr)
          errstr = ""
          CYCLE ! main loop
       END IF

       IF (first) THEN
          chi2_ = chi2
          jac_sph_inv_ = jac_sph_inv
          chi2min_prm = chi2_
          first = .FALSE.
          CYCLE ! start MCMC loop
       END IF

       !!
       !! MAKE DECISION WHETHER PROPOSED ORBIT IS ACCEPTED
       !! 
       dchi2 = chi2 - this%chi2_min_prm
       IF (burn_in) THEN
          avalue = EXP(-0.5_bp*(chi2-chi2_))*jac_sph_inv_/jac_sph_inv
       ELSE
          avalue = this%dchi2_prm/ABS(dchi2)
       END IF

       ! Decision making accept or reject depending on pdf
       IF (avalue >= 1.0_bp) THEN
          accepted = .TRUE.
          norb_acc = norb_acc + 1
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(A,1X,I0,1X,A,1X,E12.4,2X,A,6(1X,F20.10))") &
                  "Orbit", this%sor_norb_cmp  + 1, &
                  "accepted. pdf higher:", avalue, &
                  "State:", state(1), state(2:3)/rad_deg, &
                  state(4), state(5:6)/rad_deg
          END IF
       END IF

       IF (accepted) THEN
          !!
          !! PROPOSED ORBIT ACCEPTED
          !!
          n0_ = n0
          WHERE (n0_ == 0)
             n0_ = 1
          END WHERE
          ! - rms:
          rms = SQRT(SUM(residuals(:,:)**2.0_bp,dim=1,mask=this%obs_masks_prm)/n0_)

          burnin_case = .FALSE.
          IF (burn_in) THEN
             IF (info_verb >= 2) THEN
                IF (norb_acc == 1) THEN
                   WRITE(stdout,"(A)") "Burn-in: norb_acc, chi2, chi2min, dchi2"
                END IF
                WRITE(stdout,"(10X,I2,1X,3(F5.1,1X))") norb_acc, chi2, this%chi2_min_prm, dchi2 !, chi2min_prm
             END IF
             ! RANDOM-WALK BURN IN: Require one accepted orbit in the area defined by dchi2_prm.
             burnin_case = (dchi2 <= this%dchi2_prm)
             IF (chi2 <= chi2min_prm) THEN
                chi2min_prm = chi2
             END IF
             IF (burnin_case) THEN
                burn_in = .FALSE.
                IF (info_verb >= 2) THEN
                   WRITE(stdout,"(A,2(F5.1,1X))") 'BURN-IN OVER (dchi2, dchi2_prm): ', &
                        dchi2, this%dchi2_prm
                   WRITE(stdout,"(A,I0,1X,A)") ' after checking ',norb_acc,'accepted orbits '
                END IF
             END IF
          END IF

          ! Only save orbits after burn_in is over:
          IF (.NOT.burn_in) THEN
             this%sor_norb_cmp = this%sor_norb_cmp  + 1
             this%rchi2_arr_cmp(this%sor_norb_cmp) = chi2 - REAL(ndof,bp)
             this%pdf_arr_cmp(this%sor_norb_cmp) =  &
                  EXP(-0.5_bp * (chi2 - REAL(ndof,bp))) / &
                  jac_sph_inv
             this%res_arr_cmp(this%sor_norb_cmp,1:nobs,1:6) = residuals(1:nobs,1:6)
             this%orb_arr_cmp(this%sor_norb_cmp) = copy(orb)
             this%sor_rho_arr_cmp(this%sor_norb_cmp,1) = state(1)
             ! From ranging:
             ! "- generated distance to object at second epoch (NB! Relative
             !   to the distance at the first epoch.)"
             this%sor_rho_arr_cmp(this%sor_norb_cmp,2) = state(4) ! - state(4)
          ELSE IF (norb_acc > burnin_max) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / randomWalkRanging", &
                  "Burn in failed!", 1)
             RETURN
          END IF
          chi2_ = chi2
          state_ = state
          jac_sph_inv_ = jac_sph_inv

          IF (chi2min > chi2) THEN
             chi2min = chi2
          END IF
          accepted = .FALSE.

       END IF

       IF (this%sor_norb_cmp >= 1) THEN
          this%repetition_arr_cmp(this%sor_norb_cmp) = this%repetition_arr_cmp(this%sor_norb_cmp) + 1
       END IF

       IF (info_verb >= 2 .AND. MOD(this%sor_ntrial_cmp,5000) == 0) THEN
          WRITE(stdout,"(2X,A,3(I0,1X))") "Nr of accepted orbits and trials: ", &
               this%sor_norb_cmp, this%sor_ntrial_cmp
       END IF

    END DO
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(A)") "END OF MCMC RANGING"
    END IF

    this%orb_arr_cmp => reallocate(this%orb_arr_cmp, this%sor_norb_cmp)
    this%repetition_arr_cmp => reallocate(this%repetition_arr_cmp, this%sor_norb_cmp)
    this%res_arr_cmp => reallocate(this%res_arr_cmp, this%sor_norb_cmp, nobs, 6)
    this%rchi2_arr_cmp => reallocate(this%rchi2_arr_cmp, this%sor_norb_cmp)
    this%pdf_arr_cmp => reallocate(this%pdf_arr_cmp, this%sor_norb_cmp)
    this%sor_rho_arr_cmp => reallocate(this%sor_rho_arr_cmp, this%sor_norb_cmp, 2)
    this%sor_rho_cmp(1,1) = MINVAL(this%sor_rho_arr_cmp(:,1))
    this%sor_rho_cmp(1,2) = MAXVAL(this%sor_rho_arr_cmp(:,1))
    this%sor_rho_cmp(2,1) = MINVAL(this%sor_rho_arr_cmp(:,2)-this%sor_rho_arr_cmp(:,1))
    this%sor_rho_cmp(2,2) = MAXVAL(this%sor_rho_arr_cmp(:,2)-this%sor_rho_arr_cmp(:,1))
    ! Compute pdf based on repetitions and exp(chi2) and normalize the
    ! sum of the pdf
    this%pdf_arr_cmp = this%pdf_arr_cmp * this%repetition_arr_cmp
    this%pdf_arr_cmp = this%pdf_arr_cmp/SUM(this%pdf_arr_cmp)
    IF (this%sor_norb_cmp > 0) THEN
       CALL propagate(this%orb_arr_cmp, this%t_inv_prm)
       i = MAXLOC(this%pdf_arr_cmp,dim=1)
       this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       this%chi2_min_cmp = chi2min
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(3(A,1X))") "Element types for covariance matrix and ML orbit:", &
               this%cov_type_prm, getElementType(this%orb_ml_cmp)
          WRITE(stdout,"(A,1X,I0)") "Final number of orbits:", this%sor_norb_cmp
          WRITE(stdout,"(A,1X,I0)") "Final number of trials:", this%sor_ntrial_cmp
          WRITE(stdout,"(A,1X,E7.2)") "Acceptance rate: ", &
               REAL(this%sor_norb_cmp)/REAL(this%sor_ntrial_cmp)
          WRITE(stdout,"(A,1X,E8.3)") "chi2 min: ", chi2min
       END IF
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "No sample orbits found!", 1)
    END IF

    DO i=1,nobs
       CALL NULLIFY(obs_scoords(i))
       CALL NULLIFY(obsy_ccoords(i))
    END DO
    DEALLOCATE(obs_masks, obs_coords, comp_coords, cosdec0, &
         residuals, obs_scoords, obsy_ccoords, &
         information_matrix_obs, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / randomWalkRanging", &
            "Could not deallocate memory (20).", 1)
       RETURN
    END IF

  END SUBROUTINE randomWalkRanging





  !! *Description:*
  !!
  !! Markov Chain Monte Carlo ranging method. 
  !!   - Currently stable version! (Pre-2015 version, see MCMCRanging2)
  !!   - Burn in phase: discard accepted orbits until (R.A., Dec.) residual rms are reasonable.
  !!
  !! Returns error.
  !!
  SUBROUTINE MCMCRanging(this)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout) :: this

    TYPE(Orbit) :: orb
    TYPE(SphericalCoordinates), DIMENSION(:), POINTER :: &
         obs_scoords => NULL(), comp_scoords => NULL() !observed and
    !computed sky
    !positions
    TYPE(SphericalCoordinates) :: obs_scoord1, obs_scoord2, obs_help
    TYPE(CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords => NULL()
    TYPE(CartesianCoordinates) :: obs_ccoord_helio1, &
         obs_ccoord_helio2, obs_ccoord_topo1, obs_ccoord_topo2
    TYPE(Time) :: tt
    REAL(bp), DIMENSION(:,:,:), POINTER :: information_matrix_obs => NULL(), &
         partials_arr => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: residuals, &
         comp_coords, &
         obs_coords, &
         jacobians
    REAL(bp), DIMENSION(:), ALLOCATABLE :: cosdec0
    REAL(bp), DIMENSION(6,6):: information_matrix_elem, &
         jacobian_matrix
    REAL(bp), DIMENSION(6) :: elements, rans, state, state_, &
         proposal_density, rms
    REAL(bp), DIMENSION(3) :: pos
    REAL(bp) :: chi2, ran, t1, t2, obs_coord, jac_sph_inv, &
         jac_sph_inv_, chi2_, chi2min, avalue, a, q
    INTEGER, DIMENSION(:,:), POINTER :: obs_pair_arr => NULL()
    INTEGER, DIMENSION(6) :: n0, n0_
    INTEGER, DIMENSION(2) :: obs_pair
    INTEGER :: ndof, i, itrial, err, nobs, j, count_rms, norb_acc
    INTEGER, PARAMETER :: rms_prm = 3    ! 3 or 5? Make user-defined parameter? this%rms_prm
    INTEGER, PARAMETER :: burnin_max = 50    ! Make user-defined parameter?
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks => NULL()
    LOGICAL :: accepted, first = .TRUE., burn_in = .TRUE., rms_found = .FALSE.

    IF (info_verb >= 2) THEN
       WRITE(*,*)"Starting MCMC ranging"
    END IF

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Allocate memory for solution storing - orbits
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DEALLOCATE(this%orb_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "Could not deallocate memory (5).", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(this%orb_arr_cmp(this%sor_norb_prm), stat=err) 
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "Could not allocate memory (5).", 1)
       RETURN
    END IF
    nobs = getNrOfObservations(this%obss)
    ALLOCATE(this%res_arr_cmp(this%sor_norb_prm,nobs,6), &
         this%pdf_arr_cmp(this%sor_norb_prm), &
         residuals(nobs,6), obs_coords(nobs,6), &
         comp_coords(nobs, 6), cosdec0(nobs), &
         this%reg_apr_arr_cmp(this%sor_norb_prm), &
         this%jac_arr_cmp(this%sor_norb_prm,3), &
         this%rchi2_arr_cmp(this%sor_norb_prm), &
         this%rms_arr_cmp(this%sor_norb_prm,6), &
         this%repetition_arr_cmp(this%sor_norb_prm), &
         this%sor_rho_arr_cmp(this%sor_norb_prm,2), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "Could not allocate memory (10).", 1)
       RETURN 
    END IF
    this%repetition_arr_cmp = 0
    ! Allocate memory for solution storing - pdf 
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "Could not deallocate memory (10).", 1)
          RETURN
       END IF
    END IF
    ALLOCATE(this%pdf_arr_cmp(this%sor_norb_prm), stat=err) 
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "Could not allocate memory (15).", 1)
       RETURN
    END IF
    this%sor_rho_cmp(:,1) = HUGE(this%sor_rho_cmp)
    this%sor_rho_cmp(:,2) = -HUGE(this%sor_rho_cmp)

    ! Get observed sky positions (original observational data) 
    obs_scoords => getObservationSCoords(this%obss) 
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    obs_coords = 0.0_bp
    DO i=1,nobs
       CALL rotateToEquatorial(obs_scoords(i))
       obs_coords(i,:) = getCoordinates(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (10)", 1)
          CALL NULLIFY(obs_scoords(i))
          RETURN
       END IF
       ! CALL NULLIFY(obs_scoords(i))
       cosdec0(i) = COS(obs_coords(i,3))
    END DO

    obs_masks => getObservationMasks(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    ! Observation number counter (observation mask must be up-to-date!), 
    ! construct cosine array:
    DO i=1,6
       n0(i) = COUNT(obs_masks(:,i))
    END DO
    ndof = COUNT(obs_masks) - 6

    information_matrix_obs => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "TRACE BACK (20)", 1)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(comp_coords, stat=err)
       RETURN
    END IF

    obsy_ccoords => getObservatoryCCoords(this%obss) 
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF

    ! Use predefined observation pair if not using random selection
    IF (.NOT. this%sor_random_obs_prm) THEN
       obs_pair = RESHAPE(this%sor_pair_arr_prm, (/ 2 /))
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "Use of random observation pairs not accepted.", 1)
       RETURN
    END IF
    ! Copy observation pair
    obs_scoord1 = copy(obs_scoords(obs_pair(1)))
    obs_scoord2 = copy(obs_scoords(obs_pair(2)))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "TRACE BACK (30)", 1)
       RETURN
    END IF

    ! proposal density; assuming 3sigma rho bounds
    proposal_density(1) = (this%sor_rho_prm(1,2) - this%sor_rho_prm(1,1))/6.0_bp
    proposal_density(2) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(1),2,2))
    proposal_density(3) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(1),3,3))
    proposal_density(4) = (this%sor_rho_prm(2,2) - this%sor_rho_prm(2,1))/6.0_bp 
    proposal_density(5) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(2),2,2))
    proposal_density(6) = 1.0_bp / SQRT(information_matrix_obs(obs_pair(2),3,3))
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(A)") "Proposal density for the first observation (y_l [au], RA [as], Dec [as]): "
       WRITE(stdout,"(F15.7)") proposal_density(1)
       WRITE(stdout,"(F15.7)") proposal_density(2)/rad_asec
       WRITE(stdout,"(F15.7)") proposal_density(3)/rad_asec
       WRITE(stdout,"(A)") "Proposal density for the second observation (y_r [au], RA [as], Dec [as]): "
       WRITE(stdout,"(F15.7)") proposal_density(4)
       WRITE(stdout,"(F15.7)") proposal_density(5)/rad_asec
       WRITE(stdout,"(F15.7)") proposal_density(6)/rad_asec
    END IF

    ! rho1, RA1, Dec1, drho12, RA2, Dec2 of the proposed orbit
    state_(1:3) = getPosition(obs_scoord1)
    state_(4:6) = getPosition(obs_scoord2)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "TRACE BACK (50)", 1)
       DEALLOCATE(obs_masks, stat=err)
       RETURN
    END IF
    ! Center the first trial to the ml orbit if available
    IF (exist(this%orb_ml_prm)) THEN
       CALL getEphemerides(this%orb_ml_prm, obsy_ccoords, comp_scoords)
       state_(1) = getDistance(comp_scoords(obs_pair(1)))
       state_(4) = getDistance(comp_scoords(obs_pair(2)))
    ELSE
       state_(1) = (this%sor_rho_prm(1,2) + this%sor_rho_prm(1,1))/2.0_bp
       state_(4) = state_(1) + (this%sor_rho_prm(2,2) + this%sor_rho_prm(2,1))/2.0_bp
    END IF

    ! set counter some parameters
    this%sor_norb_cmp = 0
    this%sor_ntrial_cmp = -1 ! first loop is setting the parameters
    chi2min = HUGE(chi2min)
    first = .TRUE.
    itrial = 0

    norb_acc=0
    count_rms=0
    burn_in = .TRUE.
    !    burn_in = .false.
    DO WHILE (this%sor_norb_cmp < this%sor_norb_prm .AND. &
         this%sor_ntrial_cmp < this%sor_ntrial_prm)

       ! Get next state(topocentric distances and sky positions for
       ! the two observations) around the previous computed state
       IF (itrial == 0 .AND. info_verb >= 3) THEN
          WRITE(stdout,"(A)") "Generating new state..."
       END IF

       this%sor_ntrial_cmp = this%sor_ntrial_cmp + 1
       itrial = itrial + 1
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(A,I0)") "Trial state #", itrial
       END IF

       IF (first) THEN
          state = state_
       ELSE
          CALL randomGaussian(rans)
          state(1) = state_(1) + rans(1) * proposal_density(1) + rans(4) * proposal_density(4)
          state(2:3) = state_(2:3) + rans(2:3) * proposal_density(2:3)
          state(4) = state_(4) + rans(1) * proposal_density(1) - rans(4) * proposal_density(4)
          state(5:6) = state_(5:6) + rans(5:6) * proposal_density(5:6)
       END IF

       ! Check if the distances are negative
       IF (state(1) <= 0.0_bp .OR. state(4) <= 0.0_bp) THEN
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(2X,A,F10.7,A)") &
                  "Failed (One or both topocentric distances smaller than the Earth radius.)" 
          END IF
          CYCLE
       END IF

       ! Create new spherical coordinates with generated states in equatorial frame
       CALL NULLIFY(obs_scoord1)
       CALL NULLIFY(obs_scoord2)
       CALL NEW(obs_scoord1, state(1), state(2), state(3), getTime(obsy_ccoords(obs_pair(1))))
       CALL NEW(obs_scoord2, state(4), state(5), state(6), getTime(obsy_ccoords(obs_pair(2))))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (65)", 1)
          RETURN
       END IF

       ! Create new topocentric cartesian coordinates of the observations in equatorial frame
       CALL NULLIFY(obs_ccoord_topo1)
       CALL NULLIFY(obs_ccoord_topo2)
       CALL NEW(obs_ccoord_topo1, obs_scoord1) 
       CALL NEW(obs_ccoord_topo2, obs_scoord2) 
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (70)", 1)
          RETURN
       END IF

       ! Rotate to ecliptic
       CALL rotateToEcliptic(obs_ccoord_topo1)
       CALL rotateToEcliptic(obs_ccoord_topo2)
       CALL rotateToEcliptic(obsy_ccoords(obs_pair(1)))
       CALL rotateToEcliptic(obsy_ccoords(obs_pair(2)))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (75)", 1)
          RETURN
       END IF

       ! Cartesian heliocentric coordinates
       CALL NULLIFY(obs_ccoord_helio1)
       CALL NULLIFY(obs_ccoord_helio2)
       obs_ccoord_helio1 = copy(obsy_ccoords(obs_pair(1)) + obs_ccoord_topo1)
       obs_ccoord_helio2 = copy(obsy_ccoords(obs_pair(2)) + obs_ccoord_topo2)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (80)", 1)
          RETURN
       END IF

       CALL estimateLightTime(obs_ccoord_helio1, state(1))! changed from 3 and 6
       CALL estimateLightTime(obs_ccoord_helio2, state(1)+state(4))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (85)", 1)
          RETURN
       END IF

       ! Find orbit candidate at the epoch of the first observation by
       ! using the chosen method to solve the 2-point boundary value
       ! problem:
       CALL NULLIFY(orb)
       CALL NEW(orb, obs_ccoord_helio1, obs_ccoord_helio2, &
            this%sor_2point_method_prm, this%apriori_a_max_prm, &
            center=this%center_prm, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (86)", 1)
          error = .FALSE.
          CYCLE
       END IF

       IF (this%informative_apriori_prm .AND. .NOT. first) THEN
          ! Semimajor axis:
          a = -1.0_bp
          IF (this%apriori_a_min_prm >= 0.0_bp .OR. &
               this%apriori_a_max_prm >= 0.0_bp) THEN
             a = getSemimajorAxis(orb)
             IF (this%apriori_a_min_prm >= 0.0_bp .AND. &
                  a < this%apriori_a_min_prm) THEN
                ! Semimajor axis too small
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (semimajor axis too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             IF (this%apriori_a_max_prm >= 0.0_bp .AND. &
                  a > this%apriori_a_max_prm) THEN
                ! Semimajor axis too large
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (semimajor axis too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
          ! Periapsis distance:
          IF (this%apriori_periapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_periapsis_max_prm >= 0.0_bp) THEN
             IF (a >= 0.0_bp) THEN
                CALL getPeriapsisDistance(orb, q, a)
             ELSE
                CALL getPeriapsisDistance(orb, q)
             END IF
             IF (error) THEN
                error = .FALSE.
                CYCLE
             END IF
             ! Periapsis distance too small:
             IF (this%apriori_periapsis_min_prm >= 0.0_bp .AND. &
                  q < this%apriori_periapsis_min_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (periapsis distance too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             ! Periapsis distance too large:
             IF (this%apriori_periapsis_max_prm >= 0.0_bp .AND. &
                  q > this%apriori_periapsis_max_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (periapsis distance too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
          ! Apoapsis distance:
          IF (this%apriori_apoapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_apoapsis_max_prm >= 0.0_bp) THEN
             IF (a >= 0.0_bp) THEN
                CALL getApoapsisDistance(orb, Q, a)
             ELSE
                CALL getApoapsisDistance(orb, Q)
             END IF
             ! Apoapsis distance too small:
             IF (this%apriori_apoapsis_min_prm >= 0.0_bp .AND. &
                  Q < this%apriori_apoapsis_min_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (apoapsis distance too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             ! Apoapsis distance too large:
             IF (this%apriori_apoapsis_max_prm >= 0.0_bp .AND. &
                  Q > this%apriori_apoapsis_max_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (apoapsis distance too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
       END IF

       CALL setParameters(orb, &
            dyn_model=this%dyn_model_prm, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)

       ! checking if the epochs are the same
       IF (.NOT. equal(this%t_inv_prm,getTime(orb))) THEN
          CALL propagate(orb, this%t_inv_prm)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / MCMCRanging", &
                  "TRACE BACK (90)", 1)
             error = .FALSE.
             CYCLE
          END IF
       END IF

       IF (info_verb >= 3) THEN
          WRITE(stdout,"(A,1X,I0,1X,A)") "New state generated using", itrial, "trials."
       END IF
       itrial = 0

       !! 
       !! COMPUTE CHI2 AND PDF FOR PROPOSED ORBIT
       !!

       ! get computed sky positions - spherical coordinates
       CALL getEphemerides(orb, obsy_ccoords, comp_scoords, partials_arr=partials_arr)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (45)", 1)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          !DEALLOCATE(obsy_ccoords, stat=err)
          !DEALLOCATE(obs_masks, stat=err)
          !DEALLOCATE(residuals, stat=err)
          !DEALLOCATE(information_matrix_obs, stat=err)
          !DEALLOCATE(obs_coords, stat=err)         
          !DEALLOCATE(comp_coords, stat=err)
          !RETURN
          error = .FALSE.
          CYCLE
       END IF

       ! Residuals
       comp_coords = 0.0_bp
       DO i=1,nobs
          CALL rotateToEquatorial(comp_scoords(i))
          comp_coords(i,:) = getCoordinates(comp_scoords(i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / MCMCRanging", &
                  "TRACE BACK (55)", 1)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(obs_masks, stat=err)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(comp_coords, stat=err)
             RETURN
          END IF
          CALL NULLIFY(comp_scoords(i))
       END DO
       DEALLOCATE(comp_scoords)
       residuals(1:nobs,1:6) = obs_coords(1:nobs,1:6) - &
            comp_coords(1:nobs,1:6)        
       residuals(1:nobs,2) = residuals(1:nobs,2) * &
            COS(obs_coords(1:nobs,3))
       DO i=1,nobs
          IF (ABS(residuals(i,2)) > pi) THEN
             obs_coord = obs_coords(i,2)
             IF (obs_coord < comp_coords(i,2)) THEN
                obs_coord = obs_coord + two_pi
             ELSE
                comp_coords(i,2) = comp_coords(i,2) + two_pi
             END IF
             residuals(i,2) = (obs_coord - &
                  comp_coords(i,2)) * COS(obs_coords(i,3))
          END IF
       END DO

       chi2 = chi_square(residuals, information_matrix_obs, errstr=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "TRACE BACK (60)", 1)
          DEALLOCATE(obs_masks, stat=err)
          RETURN
       END IF

       ! Jacobians
       ! Determinant of Jacobian between topocentric
       ! coordinates (inverse problem coordinates)
       ! and orbital parameters required for output
       ! ("Topocentric Wrt Cartesian/Keplerian"):
       jacobian_matrix(1:3,:) = partials_arr(1:3,:,obs_pair(1)) / &
            cosdec0(obs_pair(1))
       jacobian_matrix(4:6,:) = partials_arr(1:3,:,obs_pair(2)) / &
            cosdec0(obs_pair(2))
       DEALLOCATE(partials_arr, stat=err)
       jac_sph_inv = ABS(determinant(jacobian_matrix, errstr))
       IF (LEN_TRIM(errstr) /= 0) THEN
          CALL errorMessage("StochasticOrbit / MCMCRanging", &
               "Unsuccessful computation of determinant of orbital element " // &
               "jacobian matrix: " // TRIM(errstr), 1)
          CALL matrix_print(jacobian_matrix, stderr, errstr)
          errstr = ""
          CYCLE ! main loop
       END IF

       IF (first) THEN
          chi2_ = chi2
          jac_sph_inv_ = jac_sph_inv
          first = .FALSE.
          CYCLE ! start MCMC loop
       END IF

       !!
       !! MAKE DECISION WHETHER PROPOSED ORBIT IS ACCEPTED
       !! 
       avalue = EXP(-0.5_bp*(chi2-chi2_))*jac_sph_inv_/jac_sph_inv

       ! JV 16.2.2015: do we ever go here???
       IF (avalue > HUGE(avalue)) THEN
          IF (info_verb >= 4) THEN
             WRITE(stdout,*) "huge avalue", MIN(chi2,chi2_)
          END IF
          avalue = chi2_/chi2
          burn_in = .TRUE.
       END IF

       ! Decision making accept or reeject depending on pdf
       IF (avalue >= 1.0_bp) THEN
          accepted = .TRUE.
          norb_acc = norb_acc + 1
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(A,1X,I0,1X,A,1X,E12.4,2X,A,6(1X,F20.10))") &
                  "Orbit", this%sor_norb_cmp  + 1, &
                  "accepted. pdf higher:", avalue, &
                  "State:", state(1), state(2:3)/rad_deg, &
                  state(4), state(5:6)/rad_deg
          END IF
       ELSE IF (.NOT. burn_in) THEN
          CALL randomNumber(ran)
          IF (avalue > ran) THEN
             accepted = .TRUE.
             norb_acc = norb_acc + 1
             IF (info_verb >= 3) THEN
                WRITE(stdout,"(A,1X,I0,1X,A,1X,E12.4,2X,A,6(1X,F20.10))") &
                     "Orbit", this%sor_norb_cmp  + 1, &
                     "accepted. pdf random:", avalue, &
                     "State:", state(1), state(2:3)/rad_deg, &
                     state(4), state(5:6)/rad_deg
             END IF
          END IF
       END IF

       IF (accepted) THEN

          !!
          !! PROPOSED ORBIT ACCEPTED
          !!

          n0_ = n0
          WHERE (n0_ == 0)
             n0_ = 1
          END WHERE
          ! - rms:
          rms = SQRT(SUM(residuals(:,:)**2.0_bp,dim=1,mask=this%obs_masks_prm)/n0_)

          ! CHECK residuals: if larger than acceptance windows, burn in is not over yet.
          ! Require RMS_PRM (typically, 3-5) subsequently accepted orbits.
          IF (burn_in .AND. rms(2) < this%res_accept_prm(1,2) .AND. rms(3) < this%res_accept_prm(1,3) &
               .AND. count_rms < rms_prm) THEN
             count_rms = count_rms + 1
             SELECT CASE (count_rms)
             CASE(1) 
                rms_found = .TRUE.
             CASE (rms_prm)
                burn_in = .FALSE.
                IF (info_verb >= 2) THEN
                   WRITE(stdout,"(A,2(L1,1X),I0,1X,2(F6.3))") 'BURN-IN OVER: ', burn_in,rms_found, &
                        count_rms,rms(2:3)/rad_asec
                   WRITE(stdout,"(A,I0,1X,A)") ' after checking ',norb_acc,'accepted orbits'
                   WRITE(stdout,"(A,2(F5.2))") ' against residual acceptance windows ',&
                        this%res_accept_prm(1,2:3)/rad_asec
                END IF
                !               CASE default
                !                  IF (rms_found) count_rms = count_rms + 1
             END SELECT
          ELSE
             rms_found = .FALSE.
             count_rms = 0
          END IF
          IF (burn_in) THEN
             WRITE(stdout,"(A,3(I0,1X),2(F6.3,1X))") "Burn in: ", norb_acc, &
                  this%sor_ntrial_cmp, count_rms, rms(2:3)/rad_asec
          END IF

          ! Only save orbits after burn_in is over:
          IF (.NOT.burn_in) THEN
             this%sor_norb_cmp = this%sor_norb_cmp  + 1
             this%rchi2_arr_cmp(this%sor_norb_cmp) = chi2 - REAL(ndof,bp)
             this%pdf_arr_cmp(this%sor_norb_cmp) =  &
                  EXP(-0.5_bp * (chi2 - REAL(ndof,bp))) / &
                  jac_sph_inv
             this%res_arr_cmp(this%sor_norb_cmp,1:nobs,1:6) = residuals(1:nobs,1:6)
             this%orb_arr_cmp(this%sor_norb_cmp) = copy(orb)
             this%sor_rho_arr_cmp(this%sor_norb_cmp,1) = state(1)
             ! From ranging:
             ! "- generated distance to object at second epoch (NB! Relative
             !   to the distance at the first epoch.)"
             this%sor_rho_arr_cmp(this%sor_norb_cmp,2) = state(4) ! - state(1)
          ELSE IF (norb_acc > burnin_max) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / MCMCRanging", &
                  "Burn in failed!", 1)
             RETURN
          END IF

          chi2_ = chi2
          state_ = state
          jac_sph_inv_ = jac_sph_inv

          IF (chi2min > chi2) THEN
             chi2min = chi2
          END IF
          accepted = .FALSE.

       END IF

       IF (this%sor_norb_cmp >= 1) THEN
          this%repetition_arr_cmp(this%sor_norb_cmp) = this%repetition_arr_cmp(this%sor_norb_cmp) + 1 
       END IF

       IF (info_verb >= 2 .AND. MOD(this%sor_ntrial_cmp,5000) == 0) THEN
          WRITE(stdout,"(2X,A,3(I0,1X))") "Nr of accepted orbits and trials: ", &
               this%sor_norb_cmp, this%sor_ntrial_cmp
       END IF

    END DO
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(A)") "END OF MCMC RANGING"
    END IF

    this%orb_arr_cmp => reallocate(this%orb_arr_cmp, this%sor_norb_cmp)
    this%repetition_arr_cmp => reallocate(this%repetition_arr_cmp, this%sor_norb_cmp)
    this%res_arr_cmp => reallocate(this%res_arr_cmp, this%sor_norb_cmp, nobs, 6)
    this%rchi2_arr_cmp => reallocate(this%rchi2_arr_cmp, this%sor_norb_cmp)
    this%pdf_arr_cmp => reallocate(this%pdf_arr_cmp, this%sor_norb_cmp)
    this%sor_rho_arr_cmp => reallocate(this%sor_rho_arr_cmp, this%sor_norb_cmp, 2)
    this%sor_rho_cmp(1,1) = MINVAL(this%sor_rho_arr_cmp(:,1))
    this%sor_rho_cmp(1,2) = MAXVAL(this%sor_rho_arr_cmp(:,1))
    this%sor_rho_cmp(2,1) = MINVAL(this%sor_rho_arr_cmp(:,2)-this%sor_rho_arr_cmp(:,1))
    this%sor_rho_cmp(2,2) = MAXVAL(this%sor_rho_arr_cmp(:,2)-this%sor_rho_arr_cmp(:,1))
    ! Compute pdf based on repetitions
    this%pdf_arr_cmp = REAL(this%repetition_arr_cmp,bp)/SUM(this%repetition_arr_cmp)
    IF (this%sor_norb_cmp > 0) THEN
       CALL propagate(this%orb_arr_cmp, this%t_inv_prm)
       !    i = MAXLOC(this%pdf_arr_cmp,dim=1)
       ! USE best fitting i.e. min chi2 orbit?
       i = MINLOC(this%rchi2_arr_cmp,dim=1)
       this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(3(A,1X))") "Element types for covariance matrix and ML orbit:", &
               this%cov_type_prm, getElementType(this%orb_ml_cmp)
          WRITE(stdout,"(A,1X,I0)") "Final number of orbits:", this%sor_norb_cmp
          WRITE(stdout,"(A,1X,I0)") "Final number of trials:", this%sor_ntrial_cmp
          WRITE(stdout,"(A,1X,E7.2)") "Acceptance rate: ", &
               REAL(this%sor_norb_cmp)/REAL(this%sor_ntrial_cmp)
          WRITE(stdout,"(A,1X,E8.3)") "chi2 min: ", chi2min
       END IF
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "No sample orbits found!", 1)
    END IF

    DO i=1,nobs
       CALL NULLIFY(obs_scoords(i))
       CALL NULLIFY(obsy_ccoords(i))
    END DO
    DEALLOCATE(obs_masks, obs_coords, comp_coords, cosdec0, &
         residuals, obs_scoords, obsy_ccoords, &
         information_matrix_obs, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / MCMCRanging", &
            "Could not deallocate memory (20).", 1)
       RETURN
    END IF

  END SUBROUTINE MCMCRanging





  SUBROUTINE observationSampling(this, orb_arr)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this
    TYPE (Orbit), DIMENSION(:), INTENT(in) :: orb_arr    

    TYPE (StochasticOrbit) :: storb
    TYPE (Observations) :: obss
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: scoord_arr => NULL()
    TYPE (Time) :: t
    TYPE (Orbit), DIMENSION(:), ALLOCATABLE :: orb_arr_init
    TYPE (Orbit), DIMENSION(7) :: orb_arr_tmp
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=DYN_MODEL_LEN) :: dyn_model, dyn_model_
    CHARACTER(len=INTEGRATOR_LEN) :: orb_integrator
    CHARACTER(len=32) :: str
    REAL(bp), DIMENSION(:,:,:), POINTER :: cov_mat_obs => NULL(), &
         center_and_absbound_arr => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: orb_additional_perturbers => NULL(), &
         stddev_arr => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: mean_arr, coordinates_arr
    REAL(bp), DIMENSION(6) :: elements, elements_, coordinates, coordinates_
    REAL(bp) :: sigma_multiplier_rms, orb_integration_step, a_r, &
         pdv, chi2, dchi2, rchi2, ran, pdv_previous
    INTEGER :: i, j, k, err, nobs, info_verb_, iorb, ipreli
    LOGICAL :: first, accept

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "MCMC OBSERVATION SAMPLING"
       WRITE(stdout,"(1X)")
       WRITE(stdout,"(2X,A)") "Parameters:"
       WRITE(stdout,"(2X,A,1X,L1)") "outlier_rejection_prm", this%outlier_rejection_prm
       WRITE(stdout,"(2X,A,1X,F10.5)") "generat_multiplier_prm", this%generat_multiplier_prm
       WRITE(stdout,"(2X,A,1X,L1)") "generat_gaussian_deviates_prm", this%generat_gaussian_deviates_prm
       WRITE(stdout,"(2X,A,1X,I0)") "os_sampling_type_prm", this%os_sampling_type_prm
       WRITE(stdout,"(2X,A,1X,A)") "dyn_model_prm", TRIM(this%dyn_model_prm)
       WRITE(stdout,"(2X,A,1X,A)") "element_type_prm", TRIM(this%element_type_prm)
       WRITE(stdout,"(2X,A,1X,I0)") "os_norb_prm", this%os_norb_prm
       WRITE(stdout,"(2X,A,1X,I0)") "os_ntrial_prm", this%os_ntrial_prm
    END IF

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / observationSampling", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(orb_arr(1))) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / observationSampling", &
            "Preliminary orbit (#1) has not been initialized.", 1)
       RETURN
    END IF

    IF (this%outlier_rejection_prm) THEN
       IF (this%outlier_multiplier_prm < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / observationSampling", &
               "Outlier criterion has not been defined.", 1)
          RETURN
       END IF
    END IF

    IF (.NOT.ASSOCIATED(this%obs_masks_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / observationSampling", &
            "Observation mask is missing.", 1)
       RETURN
    END IF

    ! Observational information
    nobs = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "observationSampling", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    IF (this%generat_gaussian_deviates_prm) THEN
       cov_mat_obs => getCovarianceMatrices(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (10)", 1)
          DEALLOCATE(cov_mat_obs, stat=err)
          RETURN
       END IF
       ALLOCATE(mean_arr(nobs,6))
       mean_arr = 0.0_bp
    ELSE
       stddev_arr => getStandardDeviations(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (15)", 1)
          DEALLOCATE(cov_mat_obs, stat=err)
          RETURN
       END IF
       ALLOCATE(center_and_absbound_arr(nobs,6,2))
       center_and_absbound_arr = 0.0_bp
       DO i=1,nobs
          DO j=2,3
             center_and_absbound_arr(i,j,2) = this%generat_multiplier_prm * stddev_arr(i,j)
          END DO
       END DO
       DEALLOCATE(stddev_arr)
    END IF
    IF (info_verb >= 3 .OR. this%os_sampling_type_prm == 2) THEN
       scoord_arr => getObservationSCoords(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (20)", 1)
          DEALLOCATE(cov_mat_obs, stat=err)
          RETURN
       END IF
       ALLOCATE(coordinates_arr(nobs,6))
       DO i=1,nobs
          coordinates_arr(i,:) = getCoordinates(scoord_arr(i))
          CALL NULLIFY(scoord_arr(i))
       END DO
       DEALLOCATE(scoord_arr)
    END IF

    ! Orbital information
    ! Inversion epoch is equal to epoch of first preliminary orbit:
    t = getTime(orb_arr(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "observationSampling", &
            "TRACE BACK (25)", 1)
       RETURN
    END IF
    CALL getParameters(orb_arr(1), dyn_model=dyn_model)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "observationSampling", &
            "TRACE BACK (30)", 1)
       RETURN
    END IF
    IF (dyn_model /= this%dyn_model_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "observationSampling", &
            "Inconsistent propagation schemes: " // &
            "orb=" // TRIM(dyn_model) // " and storb=" // &
            TRIM(this%dyn_model_prm) // ".", 1)
       RETURN
    END IF
    dyn_model_ = dyn_model
    ! Initialize elements:
    frame = getFrame(orb_arr(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / " // &
            "observationSampling", &
            "TRACE BACK (35)", 1)
       RETURN
    END IF
    element_type = this%element_type_prm
    CALL locase(element_type, error)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / observationSampling", &
            "The element type string contains forbidden characters.", 1)
       RETURN
    END IF
    ALLOCATE(orb_arr_init(MAX(7,SIZE(orb_arr))))
    IF (SIZE(orb_arr) >= 7) THEN
       DO i=1,SIZE(orb_arr)
          orb_arr_init(i) = copy(orb_arr(i))
       END DO
    ELSE
       elements = getElements(orb_arr(1), element_type, frame=frame)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (45)", 1)
          DEALLOCATE(cov_mat_obs, stat=err)
          RETURN
       END IF
       orb_arr_init(1) = copy(orb_arr(1))
       DO i=1,6
          elements_ = elements
          elements_(i) = 1.01_bp*elements_(i)
          CALL NEW(orb_arr_init(i+1), elements_, element_type, frame, t)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "observationSampling", &
                  "TRACE BACK (50)", 1)
             DEALLOCATE(cov_mat_obs, stat=err)
             RETURN
          END IF
       END DO
       DO i=1,7
          CALL setParameters(orb_arr_init(i), &
               dyn_model=this%dyn_model_prm, &
               perturbers=this%perturbers_prm, &
               asteroid_perturbers=this%ast_perturbers_prm, &
               integration_step=this%integration_step_prm, &
               integrator=this%integrator_prm)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "observationSampling", &
                  "TRACE BACK (55)", 1)
             DEALLOCATE(cov_mat_obs, stat=err)
             DEALLOCATE(center_and_absbound_arr, stat=err)
             RETURN
          END IF
       END DO
    END IF
    DO i=1,SIZE(orb_arr_init)
       SELECT CASE (TRIM(element_type))
       CASE ("cartesian")
          CALL toCartesian(orb_arr_init(i), frame=frame)
       CASE ("keplerian")
          CALL toKeplerian(orb_arr_init(i))
       CASE default
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / observationSampling", &
               "Can not use elements of type: " // TRIM(element_type), 1)
          RETURN
       END SELECT
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (40)", 1)
          RETURN
       END IF
    END DO

    first = .TRUE.
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DO i=1,SIZE(this%orb_arr_cmp)
          CALL NULLIFY(this%orb_arr_cmp(i))
       END DO
       DEALLOCATE(this%orb_arr_cmp)
    END IF
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp)
    END IF
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp)
    END IF
    IF (ASSOCIATED(this%repetition_arr_cmp)) THEN
       DEALLOCATE(this%repetition_arr_cmp)
    END IF
    ALLOCATE(this%orb_arr_cmp(this%os_norb_prm), &
         this%rchi2_arr_cmp(this%os_norb_prm), &
         this%pdf_arr_cmp(this%os_norb_prm), &
         this%repetition_arr_cmp(this%os_norb_prm))
    this%repetition_arr_cmp = 0
    iorb = 0
    pdv_previous = 1
    DO i=0,this%os_ntrial_prm

       ! Delete working copies
       CALL NULLIFY(storb)
       DO j=1,7
          CALL NULLIFY(orb_arr_tmp(j))
       END DO

       ! Make working copy of the original observations and
       ! configuration for orbit inversion:
       storb = copy(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (65)", 1)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(center_and_absbound_arr, stat=err)
          RETURN
       END IF

       IF (.NOT.first) THEN
          ! Add noise to original observations:
          IF (this%generat_gaussian_deviates_prm) THEN
             CALL addMultinormalDeviates(storb%obss, mean_arr, this%generat_multiplier_prm**2 * cov_mat_obs)
          ELSE
             CALL addUniformDeviates(storb%obss, center_and_absbound_arr)
          END IF
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "observationSampling", &
                  "TRACE BACK (70)", 1)
             DEALLOCATE(cov_mat_obs, stat=err)
             DEALLOCATE(center_and_absbound_arr, stat=err)
             RETURN
          END IF
       END IF

       ! Make a working copy of the set of initial orbits
       IF (ipreli+7 > SIZE(orb_arr_init)) THEN
          ipreli = 0
       END IF
       DO j=1,7
          orb_arr_tmp(j) = copy(orb_arr_init(ipreli+j))
       END DO
       ipreli = ipreli + 7

       ! Run simplex on the set of modified observations
       info_verb_ = info_verb
       info_verb = info_verb - 1
       CALL simplexOrbits(storb, orb_arr_tmp)
       info_verb = info_verb_
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "WARNING: Error message from simplexOrbits", 1)
          !DEALLOCATE(cov_mat_obs, stat=err)
          !DEALLOCATE(center_and_absbound_arr, stat=err)
          !RETURN
          error = .FALSE.
          CYCLE
       END IF

       ! Compute chi2 on between the best-fitting simplex orbit and
       ! the original observations
       chi2 = getChi2(this, storb%orb_ml_cmp)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / " // &
               "observationSampling", &
               "TRACE BACK (80)", 1)
          DEALLOCATE(cov_mat_obs, stat=err)
          DEALLOCATE(center_and_absbound_arr, stat=err)
          RETURN
       END IF

       ! Compute the "reduced" chi2 by subtracting the number of
       ! measured angles (or number of measured coordinates in
       ! general)
       rchi2 = chi2 - COUNT(this%obs_masks_prm)

       ! Compute probability density value
       pdv = EXP(-0.5_bp*rchi2)

       ! If this is the first trial then the observations were the
       ! nominal ones and the computations were done only to set up
       ! the variables used in the MCMC comparison
       IF (first) THEN
          pdv_previous = pdv
          this%chi2_min_prm = chi2
          first = .FALSE.
          CYCLE
       END IF

       IF (this%mh_acceptance_prm) THEN

          ! Compute the pdv ratio between the previous accepted orbit
          ! and the current trial orbit and use the MCMC MH criterion
          ! to decide whether the trial orbit should be accepted or
          ! rejected
          a_r = pdv/MAX(TINY(a_r),pdv_previous)
          CALL randomNumber(ran)
          accept = .FALSE.
          IF (ran < a_r) THEN
             accept = .TRUE.
          END IF

       ELSE IF (this%dchi2_rejection_prm) THEN

          ! Compute dchi2 between the best fit orbit and the current
          ! trial orbit, and use that and the maximum allowed dchi2 to
          ! decide whether the trial orbit should be accepted or rejected
          dchi2 = chi2 - this%chi2_min_prm
          accept = .FALSE.
          IF (dchi2 < this%dchi2_prm) THEN
             accept = .TRUE.
          END IF

       ELSE

          accept = .TRUE.

       END IF

       ! Write cometary elements and other information for each trial
       ! orbit (incl whether it was rejected or accepted)
       IF (info_verb >= 3) THEN
          elements = getElements(storb%orb_ml_cmp, "cometary")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "observationSampling", &
                  "TRACE BACK (85)", 1)
             DEALLOCATE(cov_mat_obs, stat=err)
             RETURN
          END IF
          elements(3:5) = elements(3:5)/rad_deg
          IF (accept) THEN
             WRITE(stdout,"(A,1X,10(F15.8,1X),A,1X,I0)") &
                  "COMETARY", elements, pdv, chi2, a_r, ran, &
                  "ACCEPTED", iorb + 1
          ELSE
             WRITE(stdout,"(A,1X,10(F15.8,1X),A,1X,I0)") &
                  "COMETARY", elements, pdv, chi2, a_r, ran, &
                  "REJECTED", i - iorb
          END IF
       END IF

       ! Write offset (from original measurement) for each angle and,
       ! if using "dependence" sampling and trial orbit is accepted,
       ! update the offsets from the original position
       IF (info_verb >= 3 .OR. (accept .AND. this%os_sampling_type_prm == 2)) THEN
          scoord_arr => getObservationSCoords(storb%obss)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / " // &
                  "observationSampling", &
                  "TRACE BACK (90)", 1)
             DEALLOCATE(cov_mat_obs, stat=err)
             DEALLOCATE(center_and_absbound_arr, stat=err)
             RETURN
          END IF
          DO j=1,nobs
             coordinates = getCoordinates(scoord_arr(j))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / " // &
                     "observationSampling", &
                     "TRACE BACK (95)", 1)
                DEALLOCATE(cov_mat_obs, stat=err)
                RETURN
             END IF
             IF (accept .AND. this%os_sampling_type_prm == 2) THEN ! "dependence" sampling
                IF (this%generat_gaussian_deviates_prm) THEN
                   mean_arr(j,:) = coordinates - coordinates_arr(j,:)
                ELSE
                   center_and_absbound_arr(j,:,1) = coordinates - coordinates_arr(j,:)
                END IF
             END IF
             IF (info_verb >= 3) THEN
                WRITE(stdout,"(A,I0,2(1X,F10.5))") "NOISE IN ARCSEC FOR RA,DEC PAIR #", &
                     j, (coordinates(2:3)-coordinates_arr(j,2:3))/rad_asec
             END IF
             CALL NULLIFY(scoord_arr(j))
          END DO
          DEALLOCATE(scoord_arr)
       END IF

       ! Update the solution if the trial orbit is accepted
       IF (accept) THEN
          iorb = iorb + 1
          this%orb_arr_cmp(iorb) = copy(storb%orb_ml_cmp)
          this%rchi2_arr_cmp(iorb) = rchi2
          this%pdf_arr_cmp(iorb) = pdv
          pdv_previous = pdv
       END IF
       IF (iorb /= 0) THEN
          this%repetition_arr_cmp(iorb) = this%repetition_arr_cmp(iorb) + 1
       END IF

       IF (info_verb >= 2 .AND. MOD(iorb,100) == 0) THEN
          WRITE(stdout,"(2X,A,3(I0,1X))") "Nr of accepted orbits and trials: ", &
               iorb, i
       END IF

       ! Exit the loop when enough sample orbits have been found
       IF (iorb == this%os_norb_prm) THEN
          EXIT
       END IF

    END DO

    this%orb_arr_cmp => reallocate(this%orb_arr_cmp, iorb)
    this%rchi2_arr_cmp => reallocate(this%rchi2_arr_cmp, iorb)
    this%pdf_arr_cmp => reallocate(this%pdf_arr_cmp, iorb)
    this%repetition_arr_cmp => reallocate(this%repetition_arr_cmp, iorb)
    ! Compute pdf based on repetitions
    this%pdf_arr_cmp = REAL(this%repetition_arr_cmp,bp)/SUM(this%repetition_arr_cmp)
    DO i=1,7
       CALL NULLIFY(orb_arr_tmp(i))
       CALL NULLIFY(orb_arr_init(i))
    END DO
    DEALLOCATE(orb_arr_init, stat=err)
    DEALLOCATE(mean_arr, stat=err)
    DEALLOCATE(cov_mat_obs, stat=err)
    DEALLOCATE(center_and_absbound_arr, stat=err)
    CALL NULLIFY(storb)

  END SUBROUTINE observationSampling





  SUBROUTINE autoVolumeOfVariation(this, preliminary_orbit)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), INTENT(in)              :: preliminary_orbit
    TYPE (Orbit) :: orb
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type
    CHARACTER(len=64) :: frmt = "(2X,A,1X,6(F7.3,1X))"
    REAL(bp), DIMENSION(:,:), ALLOCATABLE:: element_arr, histo
    REAL(bp), DIMENSION(:), ALLOCATABLE::  elem_data
    REAL(bp), DIMENSION(6,2) :: scaling_cmp
    REAL(bp), DIMENSION(6) :: elements_1, elements_2
    REAL(bp) :: dchi2, ddchi2, xmin, xmax, mean1, stdev1, &
         tmp, mean2, stdev2
    INTEGER, DIMENSION(2) :: scaling
    INTEGER :: i, j, k, l, err, iiter, indx, norb, norb_final, &
         ntrial_final, imap_zero, nmax, nbin, &
         ibin_zero, imap, max_indx, nmap, scalar
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask1, mask2, mask3
    LOGICAL, DIMENSION(6,2) :: scaling_ready

    ! For automated volumeOfVariation:
    ! I) Iterate to find the unbiased interval(s) of variation for the mapping parameter(s)
    !    (iterate scaling parameter for all parameters?)
    !    1) Check for unpopulated mapping bins, increase mapping interval (in practise, psc1)
    !       if necessary and repeat, at maximum niter_max (=5?) times
    !    2) Use a smaller number of sample orbits (i.e., 500) in the iteration of
    !       the mapping interval?
    ! II) Iterate the maximum point of the rigorous p.d.f. until convergence
    !     (as in Ranging: d(dchi2)=abs(2.e0_bp*log(pdf_ml_final/pdf_ml)) < 2.0)
    !     Note: should inform the user if the global ls-solution does not correspond
    !           to the maximum point!

    iiter = -1

    orb = copy(preliminary_orbit)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    norb_final = this%vov_norb_prm
    ntrial_final = this%vov_ntrial_prm
    CALL setParameters(this, &
         vov_norb=this%vov_norb_iter_prm, &
         vov_ntrial=this%vov_ntrial_iter_prm)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

    ! Find mapping element:
    DO i=1,6
       IF (this%vov_mapping_mask_prm(i)) THEN
          indx = i
          EXIT
       END IF
    END DO
    ! Initialize scaling factors:
    scaling_cmp = this%vov_scaling_prm

    ddchi2 = HUGE(ddchi2)
    scaling_ready = .FALSE.
    DO

       iiter = iiter + 1
       this%vov_niter_cmp = iiter

       IF ((ALL(scaling_ready) .AND. &
            norb >= this%vov_norb_iter_prm .AND. & 
            ddchi2 < 2.0_bp) .OR. iiter == this%vov_niter_prm) THEN
          CALL setParameters(this, &
               vov_norb=norb_final, &
               vov_ntrial=ntrial_final)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
                  "TRACE BACK (15)", 1)
             RETURN
          END IF
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(2X)")
             WRITE(stdout,"(2X,A)")       "===================="
             WRITE(stdout,"(2X,A,1X,I0)") "STARTING FINAL ROUND"
             WRITE(stdout,"(2X,A)")       "===================="
             WRITE(stdout,"(2X)")
          END IF
       ELSE
          IF (info_verb >= 2) THEN
             WRITE(stdout,"(2X)")
             WRITE(stdout,"(2X,A)")       "===================="
             WRITE(stdout,"(2X,A,1X,I0)") "STARTING ITERATION",iiter
             WRITE(stdout,"(2X,A)")       "===================="
             WRITE(stdout,"(2X)")
          END IF
       END IF

       ! Starting from the 2nd iteration:
       IF (iiter > 0) THEN
          ! Use iterated scaling parameters
          this%vov_scaling_prm = scaling_cmp
          CALL NULLIFY(orb)
          ! Use existing ml orbit as preliminary orbit
          max_indx = MINLOC(this%rchi2_arr_cmp,dim=1)
          this%chi2_min_prm = MIN(this%chi2_min_prm,this%chi2_min_cmp)
          orb = copy(this%orb_arr_cmp(max_indx))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
          ! Remove nominal orbit
          CALL NULLIFY(this%orb_ml_cmp)
          IF (ASSOCIATED(this%cov_ml_cmp)) THEN
             DEALLOCATE(this%cov_ml_cmp, stat=err)
             IF (err /= 0) THEN
                CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
                     "Could not deallocate memory ()", 1)
                RETURN
             END IF
          END IF
       END IF

       CALL volumeOfVariation(this, orb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "TRACE BACK (25)", 1)
          RETURN       
       END IF
       norb = SIZE(this%orb_arr_cmp,dim=1)
       IF ((ALL(scaling_ready) .AND. norb >= norb_final) &
            .OR. iiter == this%vov_niter_prm) THEN
          CALL NULLIFY(orb)
          EXIT
       END IF
       dchi2 = MAXVAL(this%rchi2_arr_cmp) - MINVAL(this%rchi2_arr_cmp)
       ddchi2 = ABS(dchi2 - this%dchi2_prm)
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,F10.4)") "Delta dchi2: ", ddchi2
       END IF
       ALLOCATE(element_arr(norb,6), elem_data(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "Could not allocate memory (5).", 1)
          RETURN       
       END IF

       DO i=1,norb
          element_arr(i,1:6) = getElements(this%orb_arr_cmp(i), this%element_type_prm)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
                  "TRACE BACK (30)", 1)
             RETURN       
          END IF
       END DO

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A)") "------------"
          WRITE(stdout,"(2X,A)") "UPDATING ..."
          WRITE(stdout,"(2X,A)") "------------"       
       END IF

       scaling_cmp = this%vov_scaling_cmp
       imap_zero = NINT(this%vov_nmap_prm*this%vov_scaling_cmp(indx,1) / &
            SUM(this%vov_scaling_cmp(indx,:)))

       ! Histogram for mapping parameter
       xmin = MINVAL(this%vov_map_cmp(:,indx))
       xmax = MAXVAL(this%vov_map_cmp(:,indx))
       nmax = norb
       nbin = this%vov_nmap_prm
       ibin_zero = imap_zero
       IF (nbin > 1000) THEN
          nbin = 101
          ! Not true if unsymmetric scaling parameters
          ibin_zero = 51
       END IF
       ALLOCATE(histo(nbin,2), mask1(nbin), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "Could not allocate memory (10).", 1)
          RETURN       
       END IF
       histo = 0.0_bp
       elem_data = element_arr(:,indx)
       CALL histogram(elem_data(1:nmax), histo, xmin_in=xmin, xmax_in=xmax)
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,I0,A)") "Mapping element (", indx, "):"
          WRITE(stdout,"(2X,A,2(F15.7,1X))") "Histo min/max:", xmin, xmax
          WRITE(stdout,"(2X,A,2(F15.7,1X))") "Data min/max: ", &
               MINVAL(elem_data(1:nmax)), &
               MAXVAL(elem_data(1:nmax))
       END IF

       ! Assume Gaussian statistics and adjust the scaling parameter to 
       ! cover the 3-sigma range: 
       CALL moments(elem_data(1:nmax), mean=mean1, std_dev=stdev1, &
            errstr=errstr)
       IF (LEN_TRIM(errstr) > 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "Computation of moments failed: " // TRIM(errstr), 1)
          RETURN
       END IF
       tmp = 3.0_bp*stdev1/(this%vov_scaling_cmp(indx,1)*this%vov_map_cmp(1,6+indx))
       scaling_cmp(indx,1:2) = REAL(CEILING(tmp*this%vov_scaling_cmp(indx,1:2)),bp)
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,2(1X,F15.7))") "3-sigma_pdf vs. n-sigma_ls:", &
               3.0_bp*stdev1, &
               (this%vov_scaling_cmp(indx,1)*this%vov_map_cmp(1,6+indx))
       END IF

       ! Assume that scaling is ok, and check 
       ! whether this assumption holds:
       scaling_ready = .TRUE.
       ! Empty bins in the mapping parameter direction:
       mask1 = .FALSE.
       WHERE (histo(:,2) < 1.0_bp)
          mask1 = .TRUE.
       END WHERE
       IF (ANY(.NOT.mask1(1:CEILING(0.05_bp*nbin)))) THEN
          scaling_ready(indx,1:2) = .FALSE.
       ELSE IF (ANY(.NOT.mask1(nbin+1-CEILING(0.05_bp*nbin):nbin))) THEN
          scaling_ready(indx,1:2) = .FALSE.
       END IF
       DEALLOCATE(histo, mask1, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "Could not deallocate memory (5).", 1)
          RETURN       
       END IF

       nbin = 20
       nmap = SIZE(this%vov_map_cmp, dim=1)
       ALLOCATE(histo(nbin,2), mask1(norb), mask2(nbin), mask3(nbin), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "Could not allocate memory (15).", 1)
          RETURN       
       END IF
       DO i=1,6
          IF (i == indx) THEN
             CYCLE
          END IF
          scaling = 0.0_bp
          mask3 = .FALSE.
          DO j=1,nmap
             ! Select all points around a specific
             ! mapping point:
             mask1 = .FALSE.
             IF (j == 1) THEN
                xmin = this%vov_map_cmp(j,indx)
             ELSE
                xmin = 0.5_bp * (this%vov_map_cmp(j,indx) + &
                     this%vov_map_cmp(j-1,indx))
             END IF
             IF (j == nmap) THEN
                xmax = this%vov_map_cmp(j,indx)
             ELSE
                xmax = 0.5_bp * (this%vov_map_cmp(j,indx) + &
                     this%vov_map_cmp(j+1,indx))
             END IF
             WHERE(element_arr(:,indx) < xmax .AND. element_arr(:,indx) >= xmin)
                mask1 = .TRUE.
             END WHERE
             elem_data(1:COUNT(mask1)) = PACK(element_arr(:,i), mask1)
             xmin = this%vov_map_cmp(j,i) - &
                  this%vov_scaling_cmp(i,1)*this%vov_map_cmp(j,6+i)
             xmax = this%vov_map_cmp(j,i) + &
                  this%vov_scaling_cmp(i,2)*this%vov_map_cmp(j,6+i)
             histo = 0.0_bp
             CALL histogram(elem_data(1:COUNT(mask1)), histo, &
                  xmin_in=xmin, xmax_in=xmax)
             ! Which bins corresponding to this particular mapping
             ! element value do have solutions?
             mask2 = .TRUE.
             WHERE(histo(:,2) < 1.0_bp)
                mask2 = .FALSE.
             END WHERE
             ! Which bins covering the entire mapping space do have
             ! solutions?
             WHERE (mask2)
                mask3 = .TRUE.
             END WHERE
             IF (ANY(mask3(1:CEILING(0.1_bp*nbin))) .OR. &
                  ANY(mask3(nbin+1-CEILING(0.1_bp*nbin):nbin))) THEN
                ! The current interval needs to be widened
                EXIT
             END IF
             ! Find min and max of the topical orbital element
             ! around the mapping point (if orbits exist!):
             IF (COUNT(mask1) /= 0) THEN
                xmin = MINVAL(element_arr(:,i), mask1)
                xmax = MAXVAL(element_arr(:,i), mask1)
                ! Compare with used sampling intervals and find
                ! the smallest integer scaling factor that accepts all
                ! sample orbits from the previous round:
                scalar = CEILING(ABS((xmin - this%vov_map_cmp(j,i))/this%vov_map_cmp(j,6+i)))
                IF (COUNT(mask1) /= 0 .AND. scalar > scaling(1)) THEN 
                   scaling(1) = scalar
                   IF (info_verb >= 3) THEN
                      WRITE(stdout,"(2X,2(A,1X,I0,1X))") &
                           "Lower scaling factor for element", i, &
                           "is", scaling(1)
                   END IF
                END IF
                scalar = CEILING(ABS((xmax - this%vov_map_cmp(j,i))/this%vov_map_cmp(j,6+i)))
                IF (COUNT(mask1) /= 0 .AND. scalar > scaling(2)) THEN 
                   scaling(2) = scalar
                   IF (info_verb >= 3) THEN
                      WRITE(stdout,"(2X,2(A,1X,I0,1X))") &
                           "Upper scaling factor for element", i, &
                           "is", scaling(2)
                   END IF
                END IF
             END IF
          END DO
          IF (info_verb >= 2) THEN
             !WRITE(stdout,"(2X,A,1X,I0)")    "Nr of empty bins      =", &
             !     COUNT(mask1)
             WRITE(stdout,*) " Are first n bins empty for element ", i, &
                  "? ", .NOT.mask3(1:CEILING(0.1_bp*nbin))
             WRITE(stdout,*) " Are last n bins empty for element  ", i, &
                  "? ", .NOT.mask3(nbin+1-CEILING(0.1_bp*nbin):nbin)
          END IF
          IF (ANY(mask3(1:CEILING(0.1_bp*nbin))) .OR. &
               ANY(mask3(nbin+1-CEILING(0.1_bp*nbin):nbin))) THEN
             ! Lower end of the sampling region is too small -> 
             ! increase scaling parameter value by at most 50% or
             ! at least one unit:
             scaling_cmp(i,1) = scaling_cmp(i,1) + &
                  MAX(1.0_bp,REAL(FLOOR(0.3_bp*scaling_cmp(i,1)),bp))
             scaling_ready(i,1) = .FALSE.
             ! Upper end of the sampling region is too small -> 
             ! increase scaling parameter value by at most 50% or
             ! at least one unit:
             scaling_cmp(i,2) = scaling_cmp(i,2) + &
                  MAX(1.0_bp,REAL(FLOOR(0.3_bp*scaling_cmp(i,2)),bp))
             scaling_ready(i,2) = .FALSE.
          END IF
          IF (.NOT.ANY(mask3(1:CEILING(0.1_bp*nbin))) .AND. &
               CEILING(1.2*scaling(1)) < NINT(scaling_cmp(i,1))) THEN
             ! Sampling region may be unnecessary wide in the lower end -> 
             ! optimize scaling parameters:
             scaling_cmp(i,1) = CEILING(1.2*scaling(1))
          END IF
          IF (.NOT.ANY(mask3(nbin+1-CEILING(0.1_bp*nbin):nbin)) .AND. &
               CEILING(1.2*scaling(2)) < NINT(scaling_cmp(i,2))) THEN
             ! Sampling region may be unnecessary wide in the upper end -> 
             ! optimize scaling parameters:
             scaling_cmp(i,2) = CEILING(1.2*scaling(2))
          END IF
       END DO

       DEALLOCATE(element_arr, elem_data, mask1, mask2, mask3, &
            histo, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "Could not deallocate memory (15).", 1)
          RETURN       
       END IF

       IF (info_verb >= 2) THEN
          WRITE(stdout,TRIM(frmt)) "Adjusted scaling parameters (lo):", &
               scaling_cmp(:,1)
          WRITE(stdout,TRIM(frmt)) "Adjusted scaling parameters (up):", &
               scaling_cmp(:,2)
          WRITE(stdout,"(2X,A,6(1X,L1))") "Scaling ready?", &
               scaling_ready(:,1)
          WRITE(stdout,"(2X,A,6(1X,L1))") "Scaling ready?", &
               scaling_ready(:,2)
       END IF

    END DO

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "... ITERATIONS READY"       
       WRITE(stdout,TRIM(frmt)) "Final scaling parameters (lo):", &
            this%vov_scaling_cmp(:,1)
       WRITE(stdout,TRIM(frmt)) "Final scaling parameters (up):", &
            this%vov_scaling_cmp(:,2)
    END IF

    this%vov_scaling_ready_cmp = scaling_ready
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,6(1X,L1))") "Scaling ready?", scaling_ready(:,1)
       WRITE(stdout,"(2X,A,6(1X,L1))") "Scaling ready?", scaling_ready(:,2)
       dchi2 = MAXVAL(this%rchi2_arr_cmp) - MINVAL(this%rchi2_arr_cmp)
       ddchi2 = ABS(dchi2 - this%dchi2_prm)
       WRITE(stdout,"(2X,A,F10.4)") "Delta dchi2: ", ddchi2
       element_type = "keplerian"
       elements_1 = getElements(this%orb_ml_cmp, element_type)
       elements_1(3:6) = elements_1(3:6)/rad_deg
       elements_2 = getElements(this%orb_arr_cmp(indx), element_type)
       elements_2(3:6) = elements_2(3:6)/rad_deg
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
               "TRACE BACK (50)", 1)
          error = .FALSE.
          element_type = "cartesian"
          elements_1 = getElements(this%orb_ml_cmp, element_type)
          elements_2 = getElements(this%orb_arr_cmp(indx), element_type)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / autoVolumeOfVariation", &
                  "TRACE BACK (55)", 1)
             RETURN
          END IF
       END IF
       WRITE(stdout,"(2X,A)") "Global least-squares solution in " // &
            TRIM(element_type) // " orbital elements:"
       WRITE(stdout,"(2X,6(F15.10,1X))") elements_1
       WRITE(stdout,"(2X,A)") "Maximum likelihood solution in " // &
            TRIM(element_type) // " orbital elements:"
       WRITE(stdout,"(2X,6(F15.10,1X))") elements_2
    END IF

  END SUBROUTINE autoVolumeOfVariation





  SUBROUTINE volumeOfVariation(this, preliminary_orbit)

    ! VoV parameters: vov_nmap, vov_norb, vov_ntrial, vov_niter_max
    !                 vov_scaling_parameters(6), vov_automatic
    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), INTENT(in)              :: preliminary_orbit

    TYPE (Time) :: &
         t0, &
         t
    TYPE (Orbit), DIMENSION(:), POINTER :: &
         orb_arr => NULL()
    TYPE (Orbit), DIMENSION(this%vov_norb_prm) :: &
         orb_accepted
    TYPE (Orbit) :: &
         orb, &
         orb_global, &
         orb_local
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: &
         comp_scoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
         obs_scoords => NULL()
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoords => NULL()
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=DYN_MODEL_LEN) :: &
         storb_dyn_model, &
         orb_dyn_model
    CHARACTER(len=INTEGRATOR_LEN) :: &
         storb_integrator, &
         orb_integrator
    CHARACTER(len=64) :: &
         frmt = "(F20.15,1X)", &
         efrmt = "(E10.4,1X)"
    CHARACTER(len=64) :: &
         str
    REAL(bp), DIMENSION(:,:,:,:), POINTER :: &
         partials4 => NULL()
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         residuals3 => NULL(), &
         information_matrix_obs
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE :: &
         principal_axes
    REAL(bp), DIMENSION(:,:), POINTER :: &
         residuals2 => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         obs_coords, &
         jac_arr, &
         rms_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         debiasing_factor_map, &
         debiasing_factor_arr, &
         reg_apriori_arr, &
         rchi2_arr, &
         pdf_arr, &
         cosdec0, &
                                !chi2_local, &
         diff
    REAL(bp), DIMENSION(this%vov_nmap_prm,6) :: &
         elements_local_arr
    REAL(bp), DIMENSION(6,6) :: &
         covariance_global, &
         partial_covariance_global, &
         covariance_local, &
         information_matrix_global, &
         information_matrix_local, &
         correlation_matrix, &
         partial_inverse_stdev_global, &
         eigenvectors, jacobian_matrix
    REAL(bp), DIMENSION(6) :: &
         stdev_global, &
         partial_stdev_global, &
         elements_global, &
         elements_local, &
         elements, &
         partial_eigenvalues_global, &
         eigenvalues, &
         comp_coord, &
         storb_finite_diff, &
         orb_finite_diff, &
         ran_arr, &
         upper_limit, &
         lower_limit, &
         diff_elem
    REAL(bp) :: &
         apriori, &
         chi2_ml_global_ls, &
         partial_product_global, &
         product_local, &
         mapping_interval, &
         mapping_point, &
         chi2, &
         dchi2, &
         pdf, &
         variation, &
         storb_integration_step, &
         orb_integration_step,&
         jac_car_kep, jac_equ_kep, sma, tmp, &
         obs_, comp_
    INTEGER, DIMENSION(:), ALLOCATABLE :: &
         failed_flag, &
         imap_arr
    INTEGER, DIMENSION(6) :: &
         n0, n0_, &
         ind_arr
    INTEGER :: &
         i, &
         j, &
         k, &
         imap, &
         jmap, &
         itrial, &
         iorb, &
         mmap, &
         indx, &
         nrot, &
         err, &
         nobs, &
         nset, &
         nfailed, &
         norb, &
         naccepted, &
         ielem, &
         iobs, &
         imulti, &
         isigma, &
         imap_zero, &
         ires, &
         ipdf, &
         info_verb_, &
         err_verb_
    LOGICAL, DIMENSION(:,:), POINTER :: &
         mask_arr2 => NULL()
    LOGICAL, DIMENSION(6) :: &
         element_mask
    LOGICAL :: &
         parameters_agree, outlier_rejection_

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF
    IF (.NOT.exist(preliminary_orbit)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Preliminary orbit has not yet been initialized.", 1)
       RETURN
    END IF

    nobs = getNrOfObservations(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    IF (nobs < 2) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Less than two observations available.", 1)
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(this%res_accept_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Window for accepted residuals not set.", 1)
       RETURN
    END IF
    IF (this%vov_norb_prm < 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Required number of sample orbits not set.", 1)
       RETURN
    END IF
    CALL comparePropagationParameters(this, preliminary_orbit, &
         parameters_agree=parameters_agree)
    IF (error .OR. .NOT.parameters_agree) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (1)", 1)
       RETURN
    END IF

    IF (MOD(this%vov_nmap_prm,2) == 0) THEN
       mmap = this%vov_nmap_prm/2
       this%vov_nmap_prm = this%vov_nmap_prm + 1
    ELSE
       mmap = (this%vov_nmap_prm-1)/2
    END IF
    nobs = getNrOfObservations(this%obss)
    ALLOCATE( &
         principal_axes(this%vov_nmap_prm,6,6), &
         debiasing_factor_arr(this%vov_norb_prm), &
         debiasing_factor_map(this%vov_nmap_prm), &
         obs_coords(nobs,6), &
         cosdec0(nobs), &
         residuals3(this%vov_norb_prm,nobs,6), &
         mask_arr2(nobs,6), &
         failed_flag(8), &
         reg_apriori_arr(this%vov_norb_prm), &
         pdf_arr(this%vov_norb_prm), &
         rchi2_arr(this%vov_norb_prm), &
         jac_arr(this%vov_norb_prm,3), &
         rms_arr(this%vov_norb_prm,6), &
                                !chi2_local(this%vov_nmap_prm), &
         diff(this%vov_nmap_prm), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       RETURN
    END IF

    ! Set type of mapping elements:
    IF (LEN_TRIM(this%element_type_prm) == 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Element type missing.", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       RETURN
    ELSE       
       element_type = this%element_type_prm
    END IF

    obs_scoords => getObservationSCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (45)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       RETURN
    END IF
    DO i=1,nobs
       obs_coords(i,:) = getCoordinates(obs_scoords(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (50)",1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          RETURN
       END IF
       cosdec0(i) = COS(obs_coords(i,3))
    END DO
    obsy_ccoords => getObservatoryCCoords(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (55)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       RETURN
    END IF
    IF (this%chi2_min_prm < 0.0_bp) THEN
       !       first = .TRUE.
       IF (this%chi2_min_init_prm <= 0.0_bp) THEN
          this%chi2_min_prm = REAL(COUNT(this%obs_masks_prm),bp)
       ELSE
          this%chi2_min_prm = this%chi2_min_init_prm
       END IF
       !    ELSE
       !       first = .FALSE.
    END IF

    orb = copy(preliminary_orbit)
    frame = getFrame(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (70)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb)
       RETURN
    END IF
    IF (element_type == "keplerian") THEN
       CALL toKeplerian(orb)
    ELSE IF (element_type == "cartesian") THEN
       CALL toCartesian(orb, frame=frame)
    ELSE
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Cannot recognize element type: " // TRIM(element_type), 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb)
       RETURN       
    END IF
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (60)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb)
       RETURN
    END IF
    ! Epoch is bound to the epoch of the preliminary orbit
    t0 = getTime(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (65)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb)
       CALL NULLIFY(t0)
       RETURN
    END IF

    !!
    !! 1) INTERVAL FOR MAPPING PARAMETER VIA GLOBAL LEAST-SQUARES:
    !!

    ! Check number of mapping parameters; only one allowed for now!
    IF (COUNT(this%vov_mapping_mask_prm) /= 1) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Too many or no mapping parameters chosen.", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb)
       CALL NULLIFY(t0)
       RETURN
    END IF
    ! Mapping index
    DO i=1,6
       IF (this%vov_mapping_mask_prm(i)) THEN
          indx = i
          EXIT
       END IF
    END DO

    ! Only do least squares if (nominal) ml orbit and covariance matrix 
    ! do not yet exist
    IF (.NOT. exist(this%orb_ml_cmp) .AND. .NOT. ASSOCIATED(this%cov_ml_cmp)) THEN
       info_verb_ = info_verb
       IF (this%vov_niter_cmp > 0) THEN
          info_verb = info_verb - 1
       END IF
       CALL levenbergMarquardt(this, orb)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Global least-squares solution not found.", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          CALL NULLIFY(orb)
          CALL NULLIFY(t0)
          RETURN
       END IF
       IF (info_verb == 1 .AND. info_verb_ == 2) THEN
          t = getTime(this%orb_ml_cmp)
          str = getCalendarDateString(t, 'TDT')
          CALL NULLIFY(t)
          elements = getElements(this%orb_ml_cmp, "keplerian")
          IF (error) THEN
             error = .FALSE.
             elements = getElements(orb, "cartesian", frame="ecliptic")
             WRITE(stdout,"(2X,A)") "Cartesian ecliptic elements resulting " // &
                  "from global least squares:"
             WRITE(stdout,"(2X,6(1X,F20.13),1X,A)") elements, &
                  TRIM(str)
          ELSE
             WRITE(stdout,"(2X,A)") "Keplerian elements resulting " // &
                  "from global least squares:"
             WRITE(stdout,"(2X,6(1X,F20.13),1X,A)") elements(1:2), &
                  elements(3:6)/rad_deg, TRIM(str)
          END IF
       END IF
       info_verb = info_verb_
    END IF
    CALL NULLIFY(orb)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (90)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(t0)
       RETURN
    END IF
    ! LS orbit
    orb_global = getNominalOrbit(this)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (95)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    ! Orbital elements at the specified epoch:
    elements_global = getElements(orb_global, element_type)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (100)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    ! Correlation/Standard deviation matrix:
    covariance_global = getCovarianceMatrix(this, element_type, frame)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (105)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    DO i=1,6
       stdev_global(i) = SQRT(covariance_global(i,i))
    END DO
    ! Maximum likelihood point
    information_matrix_global = matinv(covariance_global, errstr, "Cholesky")
    IF (LEN_TRIM(errstr) > 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Error in matrix inversion: " // TRIM(errstr), 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    IF (this%regularization_prm) THEN
       ! Jeffrey's apriori
       apriori = SQRT(ABS(determinant(information_matrix_global, errstr)))
       IF (LEN_TRIM(errstr) > 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Error in computation of determinant: " // TRIM(errstr), 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(t0)
          RETURN
       END IF
    ELSE
       apriori = 1.0_bp
    END IF
    ! Residuals for the global fit:
    residuals2 => getResiduals(this, orb_global)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (130)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    information_matrix_obs => getBlockDiagInformationMatrix(this%obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (135)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    DO i=1,6
       n0(i) = COUNT(this%obs_masks_prm(:,i))
    END DO
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,6(1X,I0))") "Nr of included observations:", n0
    END IF

    ! Chi2 for the global fit
    chi2_ml_global_ls = chi_square(residuals2, information_matrix_obs, this%obs_masks_prm, errstr)
    IF (LEN_TRIM(errstr) > 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Computation of chi2 for the global solution failed: " // TRIM(errstr), 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF

    ! 
    ! S :: stdev (diagonal) matrix
    ! S^(-1) :: inverse stdev (diagonal) matrix
    ! 5x5 global covariance matrix:
    information_matrix_global(indx,:) = 0.0_bp
    information_matrix_global(:,indx) = 0.0_bp
    information_matrix_global(indx,indx) = 1.0_bp
    partial_covariance_global = matinv(information_matrix_global, errstr, "Cholesky")
    IF (LEN_TRIM(errstr) > 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Error in inversion of partial global information matrix: " // TRIM(errstr), 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    partial_stdev_global = 0.0_bp
    partial_inverse_stdev_global = 0.0_bp
    DO i=1,6
       IF (this%vov_mapping_mask_prm(i)) THEN
          CYCLE
       ELSE
          partial_stdev_global(i) = SQRT(partial_covariance_global(i,i))
          partial_inverse_stdev_global(i,i) = 1.0_bp/partial_stdev_global(i)
       END IF
    END DO
    ! C = S^(-1) Sigma S^(-1)
    correlation_matrix = MATMUL(MATMUL(partial_inverse_stdev_global, &
         partial_covariance_global), partial_inverse_stdev_global)
    ! 
    CALL eigen_decomposition_jacobi(correlation_matrix, &
         partial_eigenvalues_global, eigenvectors, nrot, errstr)
    IF (LEN_TRIM(errstr) > 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Error in eigen decomposition: " // TRIM(errstr), 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    partial_product_global = 1.0_bp
    DO i=1,6
       IF (this%vov_mapping_mask_prm(i)) THEN
          CYCLE
       END IF
       IF (partial_eigenvalues_global(i) < 0.0_bp) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Negative global eigenvalue:", 1)
          IF (err_verb >= 1) THEN
             WRITE(0,"(E15.7)") partial_eigenvalues_global
          END IF
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(t0)
          RETURN
       END IF
       partial_product_global = partial_product_global*ABS(partial_eigenvalues_global(i))
    END DO

    ! If initial estimate for chi2 corresponding to maximum likelihood
    ! solution remains unknown (indicated by a negative value), set it
    ! equal to the chi2 computed for the global least squares
    ! solution:
    IF (this%chi2_min_prm < 0.0_bp) THEN
       this%chi2_min_prm = chi2_ml_global_ls
    END IF
    IF (this%vov_niter_cmp < 1 .AND. info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "Global orbital elements " // &
            "(element type=" // TRIM(element_type) // &
            " and frame=" // TRIM(getFrame(orb_global)) // "):"
       t = getTime(orb_global)
       str = getCalendarDateString(t,"tdt")
       IF (element_type == "keplerian") THEN
          WRITE(stdout,"(2X,6"//TRIM(frmt)//",A,2X,2(A,1X))") &
               elements_global(1:2), elements_global(3:6)/rad_deg, &
               TRIM(str), TRIM(element_type), TRIM(getFrame(orb_global))
       ELSE
          WRITE(stdout,"(2X,6"//TRIM(frmt)//",A,2X,2(A,1X))") &
               elements_global, TRIM(str), TRIM(element_type), &
               TRIM(getFrame(orb_global))
       END IF
       CALL NULLIFY(t)
       WRITE(stdout,"(2X,A)") "Global standard deviations:"
       WRITE(stdout, "(2X)", advance="no")
       DO i=1,6
          IF (element_type == "keplerian" .AND. i>=3) THEN
             WRITE(stdout, TRIM(frmt),advance="no") &
                  SQRT(covariance_global(i,i))/rad_deg
          ELSE
             WRITE(stdout,TRIM(frmt),advance="no") &
                  SQRT(covariance_global(i,i))
          END IF
       END DO
       WRITE(stdout,*)
       WRITE(stdout,"(2X,A)") "Global correlation:"
       CALL matrix_print(correlation_matrix, stdout, errstr, frmt=frmt)
       WRITE(stdout,"(2X,A)") "Global 5x5 eigenvalues:"
       WRITE(stdout,"(2X,6"//TRIM(frmt)//")") partial_eigenvalues_global
       WRITE(stdout,"(2X,A)") "Global 5x5 eigenvectors:"
       CALL matrix_print(eigenvectors, stdout, errstr, frmt=frmt)
       WRITE(stdout,"(2X,A,1X,3(E10.4,1X))") "Global apriori:", apriori
       WRITE(stdout,"(2X,A,1X,3(E10.4,1X))") "Global chi2:", chi2
       WRITE(stdout,*)
    END IF

    ! -----------------------------------
    ! Precomputed map:
    ! Map plausible semimajor axis region
    ! -----------------------------------
    ! Note: 1) number of mapping points, nmap, used is an input
    !          parameter (case-dependent)
    !          (nmap => vov_nmap, typically 101 or 1001)
    !       2) mapping could be made asymmetric with respect to 
    !          the global solution, i.e., different psc1 for
    !          decreasing/increasing loop below
    !          (psc1 => vov_scaling_parameter(1))

    ! Save initial element mask:
    element_mask = this%ls_elem_mask_prm

    ! Save initial outlier rejection, do not use it here
    outlier_rejection_ = this%outlier_rejection_prm
    CALL setParameters(this, outlier_rejection=.FALSE.)

    ! Adjust scaling parameter for Keplerian elements if needed
    ! (requirements: a>0 and 0<e<1):
    IF (element_type == "keplerian") THEN
       IF (elements_global(1) - &
            this%vov_scaling_prm(1,1)*stdev_global(1) < &
            planetary_radii(11)) THEN
          WRITE(stderr,*) "Problem in semimajor axis interval:", &
               elements_global(1) - &
               this%vov_scaling_prm(1,1)*stdev_global(1)
          this%vov_scaling_prm(1,1) = (elements_global(1) - &
               planetary_radii(11))/stdev_global(1)
       END IF
       IF (elements_global(2) - &
            this%vov_scaling_prm(2,1)*stdev_global(2) <= &
            0.0_bp) THEN
          WRITE(stderr,*) "Problem in eccentricity interval:", &
               elements_global(2) - &
               this%vov_scaling_prm(2,1)*stdev_global(2)
          this%vov_scaling_prm(2,1) = (elements_global(2) - &
               EPSILON(elements_global(2)))/stdev_global(2)
       END IF
       IF (elements_global(2) + &
            this%vov_scaling_prm(2,2)*stdev_global(2) >= &
            1.0_bp) THEN
          WRITE(stderr,*) "Problem in eccentricity interval:", &
               elements_global(2) + &
               this%vov_scaling_prm(2,2)*stdev_global(2)
          this%vov_scaling_prm(2,2) = (1.0_bp - &
               EPSILON(elements_global(2) - elements_global(2))) / &
               stdev_global(2)
       END IF
    END IF

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "STARTING MAPPING ..."
       WRITE(stdout,"(2X,A,1X,6(F7.3,1X))") "Scaling parameters (lo):", this%vov_scaling_prm(:,1)
       WRITE(stdout,"(2X,A,1X,6(F7.3,1X))") "Scaling parameters (hi):", this%vov_scaling_prm(:,2)
    END IF


    !!
    !! 2) LOCAL MAPS VIA LOCAL LEAST SQUARES
    !!
    !!
    imap_zero = NINT(this%vov_nmap_prm*this%vov_scaling_prm(indx,1) / &
         SUM(this%vov_scaling_prm(indx,:)))
    mapping_interval = SUM(this%vov_scaling_prm(indx,:))*stdev_global(indx)
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A,1X,I0)") "Central mapping element:", imap_zero
       WRITE(stdout,"(2X,A,1X,F10.6)") "Mapping interval:", mapping_interval
    END IF
    ! Elements used for local least squares solutions:
    CALL setParameters(this, ls_element_mask=.NOT.this%vov_mapping_mask_prm)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (230)", 1)
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(t0)
       RETURN
    END IF
    ! Mapping parameters lower or equal to midpoint (<= global ls)
    DO imap=imap_zero,1,-1
       IF (imap == imap_zero) THEN
          ! First local solution -> start from the global solution:
          elements_local = elements_global
       ELSE
          ! Start from the nearest local solution:
          elements_local = elements_local_arr(imap+1,1:6)
       END IF
       elements_local(indx) = elements_global(indx) + &
            REAL(imap-imap_zero,bp)/(this%vov_nmap_prm-imap_zero) * &
            this%vov_scaling_prm(indx,1)*stdev_global(indx)
       CALL NEW(orb_local, elements_local, element_type, frame, t0)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (235)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       CALL setParameters(orb_local, &
            dyn_model=this%dyn_model_prm, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            integration_step=this%integration_step_prm, &
            integrator=this%integrator_prm)
       !, &
       !finite_diff=this%finite_diff_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (245)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       ! Local least-squares fits:
       info_verb_ = info_verb
       IF (this%vov_niter_cmp > 0) THEN
          info_verb = info_verb - 1
       END IF
       CALL levenbergMarquardt(this, orb_local)
       info_verb = info_verb_
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (250)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       CALL NULLIFY(orb_local)
       orb_local = copy(this%orb_ml_cmp)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (255)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       elements_local_arr(imap,1:6) = getElements(orb_local, element_type)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (260)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       ! 6x6 correlation matrix where elements on the lines
       ! corresponding to the mapping parameter are zero and the
       ! diagonal element is 1:
       covariance_local = getCovarianceMatrix(this, element_type, frame)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (265)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       ! C = S^(-1) Sigma S^(-1)
       correlation_matrix = MATMUL(MATMUL(partial_inverse_stdev_global, &
            covariance_local), partial_inverse_stdev_global)
       ! Find eigenvalues and eigenvectors:
       CALL eigen_decomposition_jacobi(correlation_matrix, &
            eigenvalues, eigenvectors, nrot, errstr)
       IF (LEN_TRIM(errstr) > 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Eigen decomposition failed: " // TRIM(errstr), 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       DO i=1,6
          principal_axes(imap,1:6,i) = SQRT(eigenvalues(i))*partial_stdev_global(1:6)*eigenvectors(1:6,i)
       END DO
       product_local = 1.0_bp
       DO i=1,6
          IF (this%vov_mapping_mask_prm(i)) THEN
             CYCLE
          END IF
          IF (eigenvalues(i) < 0.0_bp) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "Negative local eigenvalue:", 1)
             IF (err_verb >= 1) THEN
                WRITE(stderr,"(6(E15.7,1X))") eigenvalues
             END IF
          END IF
          product_local = product_local*ABS(eigenvalues(i))
       END DO
       debiasing_factor_map(imap) = SQRT(product_local/partial_product_global)
       CALL NULLIFY(orb_local)
       IF (info_verb >= 2) THEN
          IF (MOD(imap-1,100) == 0) THEN
             WRITE(stdout,"(2X,I0,1X,A)") imap_zero-imap, "local least squares computed..."
          END IF
!!$          IF (MOD(imap,100) == 0) THEN
!!$             WRITE(stdout,*) imap
!!$          END IF
       END IF
    END DO

    ! Mapping parameters larger than midpoint (> global ls)
    DO imap=imap_zero+1,this%vov_nmap_prm
       ! Start from the nearest local solution:
       elements_local = elements_local_arr(imap-1,1:6)
       elements_local(indx) = elements_global(indx) + &
            REAL(imap-imap_zero,bp)/(this%vov_nmap_prm-imap_zero) * &
            this%vov_scaling_prm(indx,2)*stdev_global(indx)
       CALL NEW(orb_local, elements_local, element_type, frame, t0)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (235)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       CALL setParameters(orb_local, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            dyn_model=this%dyn_model_prm, &
            integration_step=this%integration_step_prm, &
            integrator=this%integrator_prm)
       !, &
       !finite_diff=this%finite_diff_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (245)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       ! Local least-squares fits:
       info_verb_ = info_verb
       IF (this%vov_niter_cmp > 0) THEN
          info_verb = info_verb - 1
       END IF
       CALL levenbergMarquardt(this, orb_local)
       info_verb = info_verb_
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (250)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       CALL NULLIFY(orb_local)
       orb_local = copy(this%orb_ml_cmp)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (255)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       elements_local_arr(imap,1:6) = getElements(orb_local, element_type)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (260)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       ! 6x6 correlation matrix where elements on the lines
       ! corresponding to the mapping parameter are zero and the
       ! diagonal element is 1:
       covariance_local = getCovarianceMatrix(this, element_type, frame)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (265)", 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       ! C = S^(-1) Sigma S^(-1)
       correlation_matrix = MATMUL(MATMUL(partial_inverse_stdev_global, &
            covariance_local), partial_inverse_stdev_global)
       ! Find eigenvalues and eigenvectors:
       CALL eigen_decomposition_jacobi(correlation_matrix, &
            eigenvalues, eigenvectors, nrot, errstr)
       IF (LEN_TRIM(errstr) > 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Eigen decomposition failed: " // TRIM(errstr), 1)
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
       DO i=1,6
          principal_axes(imap,1:6,i) = SQRT(eigenvalues(i))*partial_stdev_global(1:6)*eigenvectors(1:6,i)
       END DO
       product_local = 1.0_bp
       DO i=1,6
          IF (this%vov_mapping_mask_prm(i)) THEN
             CYCLE
          END IF
          IF (eigenvalues(i) < 0.0_bp) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "Negative local eigenvalue:", 1)
             IF (err_verb >= 1) THEN
                WRITE(stderr,"(6(E15.7,1X))") eigenvalues
             END IF
          END IF
          product_local = product_local*ABS(eigenvalues(i))
       END DO
       debiasing_factor_map(imap) = SQRT(product_local/partial_product_global)
       CALL NULLIFY(orb_local)
       IF (info_verb >= 2) THEN
          IF (MOD(imap,100) == 0) THEN
             WRITE(stdout,"(2X,I0,1X,A)") imap, "local least squares computed..."
          END IF
       END IF
    END DO

    !!
    !! 3) MONTE CARLO SIMULATION USING MAPPED INTERVALS
    !!
    CALL setParameters(this, outlier_rejection=outlier_rejection_)

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "STARTING MONTE CARLO SAMPLING ..."
    END IF
    failed_flag = 0
    iorb = 0
    itrial = 0
    norb = this%vov_norb_prm
    IF (norb > norb_simult_max) THEN
       norb = norb_simult_max
    END IF
    ires = 0 ; ipdf = 0
    vov_main: DO WHILE (iorb < this%vov_norb_prm .AND. itrial < this%vov_ntrial_prm)

       itrial = itrial + norb
       naccepted = 0
       IF (.NOT.ASSOCIATED(orb_arr)) THEN
          ALLOCATE(orb_arr(norb), imap_arr(norb), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "Could not allocate memory (10).", 1)
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(principal_axes, stat=err)
             DEALLOCATE(debiasing_factor_arr, stat=err)
             DEALLOCATE(debiasing_factor_map, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(cosdec0, stat=err)
             DEALLOCATE(residuals3, stat=err)
             DEALLOCATE(mask_arr2, stat=err)
             DEALLOCATE(failed_flag, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(jac_arr, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms_arr, stat=err)
             DEALLOCATE(diff, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(residuals2, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(imap_arr, stat=err)
             CALL NULLIFY(orb_global)
             CALL NULLIFY(orb_local)
             CALL NULLIFY(t0)
             RETURN
          END IF
       END IF

       !!
       !! 4) ORBIT GENERATION
       !! 

       ! Uniform sampling for the mapping parameter over the predefined mapping inteval:

       ! NOTES:
       ! 1) Currently we use directly the precomputed map, that is, we choose the closest 
       !    mapping point of the local maximum-p.d.f.'s points and use the corresponding 
       !    local intervals of variation. Should this be done with interpolation as stated 
       !    in Muinonen et al.?
       i = 0
       DO

          CALL randomNumber(ran_arr)
          ! Generate random point along the mapping axis:
          mapping_point = elements_local_arr(1,indx) + &
               ran_arr(indx) * mapping_interval
          ! Find mapping point closest to the random point:
          imap = NINT(1+ABS(mapping_point - elements_local_arr(1,indx)) / &
               (mapping_interval/(this%vov_nmap_prm-1)))
          ! Uniform sampling for the remaining parameters using the (not
          ! interpolated!) precomputed map:
          DO ielem=1,6
             IF (this%vov_mapping_mask_prm(ielem)) THEN
                ! Use the generated mapping point
                elements(ielem) = mapping_point
             ELSE
                ! Use 5x5 covariances of the nearest mapping point
                variation = this%vov_scaling_prm(ielem,1) * &
                     SUM((2.0_bp*ran_arr - 1.0_bp)*principal_axes(imap,ielem,:))
                elements(ielem) = elements_local_arr(imap,ielem) + variation 
             END IF
          END DO
          IF (element_type == "keplerian") THEN
             IF (elements(2) < 0.0_bp) THEN
                ! Eccentricity out of bounds (or should non-elliptic orbits be accepted?).
                itrial = itrial + 1
                failed_flag(1) = failed_flag(1) + 1
                WRITE(*,*) "ran", ran_arr
                CYCLE
             END IF
             IF (elements(3) < 0.0_bp .OR. elements(3) > pi) THEN
                ! Inclination not defined.
                itrial = itrial + 1
                failed_flag(2) = failed_flag(2) + 1
                WRITE(*,*) "ran", ran_arr
                CYCLE
             END IF
             ! 0 <= angle < 2pi :
             elements(4:6) = MODULO(elements(4:6), two_pi)
          END IF
          i = i + 1
          CALL NULLIFY(orb_arr(i))
          CALL NEW(orb_arr(i), elements, element_type, frame, t0)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "TRACE BACK (315)", 1)
             error = .FALSE.
             CALL NULLIFY(orb_arr(i))
             itrial = itrial + 1
             i = i - 1
             CYCLE
          END IF
          CALL setParameters(orb_arr(i), &
               perturbers=this%perturbers_prm, &
               asteroid_perturbers=this%ast_perturbers_prm, &
               dyn_model=this%dyn_model_prm, &
               integration_step=this%integration_step_prm, &
               integrator=this%integrator_prm)
          !, &
          !finite_diff=this%finite_diff_prm)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "TRACE BACK (325)", 1)
             DO i=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(i))
             END DO
             DEALLOCATE(principal_axes, stat=err)
             DEALLOCATE(debiasing_factor_arr, stat=err)
             DEALLOCATE(debiasing_factor_map, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(cosdec0, stat=err)
             DEALLOCATE(residuals3, stat=err)
             DEALLOCATE(mask_arr2, stat=err)
             DEALLOCATE(failed_flag, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(jac_arr, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms_arr, stat=err)
             DEALLOCATE(diff, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(residuals2, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(imap_arr, stat=err)
             CALL NULLIFY(orb_global)
             CALL NULLIFY(orb_local)
             CALL NULLIFY(t0)
             RETURN
          END IF
          imap_arr(i) = imap
          IF (i == norb) THEN
             EXIT
          END IF
       END DO

       IF (ASSOCIATED(comp_scoords)) THEN
          DEALLOCATE(comp_scoords)
       END IF
       NULLIFY(comp_scoords)
       IF (ASSOCIATED(partials4)) THEN
          DEALLOCATE(partials4)
       END IF
       NULLIFY(partials4)

       !!
       !! 5) ACCEPTANCE / REJECTION OF GENERATED ORBIT
       !!
       CALL getEphemerides(orb_arr, obsy_ccoords, comp_scoords, &
            partials_arr=partials4)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (330)",1)
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(imap_arr, stat=err)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials4, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF

       DO i=1,norb
          IF (.NOT. boundOrbit(orb_arr(i), this%apriori_a_max_prm, sma)) THEN
             IF (err_verb >= 2) THEN
                WRITE(stderr,*) " PROBLEM: Unbound orbit!", sma
                CYCLE
             END IF
          END IF
          DO iobs=1,nobs
             ! Multiply RA partials with cosine of observed declination:
             partials4(i,2,:,iobs) = partials4(i,2,:,iobs)*cosdec0(iobs)
             ! Sky-plane residuals and chi-squares:
             comp_coord = getCoordinates(comp_scoords(i,iobs))
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                     "TRACE BACK (335)",1)
                DO j=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(j))
                END DO
                DEALLOCATE(principal_axes, stat=err)
                DEALLOCATE(debiasing_factor_arr, stat=err)
                DEALLOCATE(debiasing_factor_map, stat=err)
                DEALLOCATE(obs_coords, stat=err)
                DEALLOCATE(cosdec0, stat=err)
                DEALLOCATE(residuals3, stat=err)
                DEALLOCATE(mask_arr2, stat=err)
                DEALLOCATE(failed_flag, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(jac_arr, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms_arr, stat=err)
                DEALLOCATE(diff, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(residuals2, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(imap_arr, stat=err)
                DEALLOCATE(comp_scoords, stat=err)
                DEALLOCATE(partials4, stat=err)
                RETURN
             END IF
             residuals3(iorb+1,iobs,1:6) = obs_coords(iobs,1:6) - comp_coord(1:6)
             residuals3(iorb+1,iobs,2) = residuals3(iorb+1,iobs,2) * cosdec0(iobs)
             IF (ABS(residuals3(iorb+1,iobs,2)) > pi) THEN
                obs_ = obs_coords(iobs,2)
                comp_ = comp_coord(2)
                IF (obs_ < comp_) THEN
                   obs_ = obs_ + two_pi
                ELSE
                   comp_ = comp_ + two_pi
                END IF
                residuals3(iorb+1,iobs,2) = (obs_ - comp_) * cosdec0(iobs)
             END IF
             IF (info_verb >= 4) THEN
                WRITE(stdout,"(2X,A,3"//TRIM(frmt)//")") "  observed pos.", obs_coords(iobs,1:3)
                WRITE(stdout,"(2X,A,3"//TRIM(frmt)//")") "  computed pos.", comp_coord(1:3)
             END IF
          END DO
          mask_arr2 = .FALSE.
          WHERE (this%obs_masks_prm .AND. ABS(residuals3(iorb+1,:,:)) > this%res_accept_prm)
             mask_arr2 = .TRUE.
          END WHERE
          IF (info_verb >= 4) THEN
             DO iobs=1,nobs
                WRITE(stdout,"(2X,A,2"//TRIM(frmt)//")") "  O-C residuals (RA, Dec):", &
                     residuals3(iorb+1,iobs,2:3)/rad_asec
             END DO
             WRITE(stdout,"(2X,A,I0,A,I0)") &
                  " No of omitted obs/included obs: ", &
                  COUNT(mask_arr2),"/",n0(2)
          END IF
          IF (COUNT(mask_arr2) > 0) THEN
             ! Residuals are too large for at least one observation.
             failed_flag(3) = failed_flag(3) + 1
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A)") &
                     " Failed (residuals are too large)"
             END IF
             elements = getElements(orb_arr(i),element_type)
             IF (elements(indx) < elements_global(indx)) THEN
                ires = ires + 1
             END IF
             CYCLE
          END IF

          ! Compute chi2:
          chi2 = chi_square(residuals3(iorb+1,:,:), information_matrix_obs, this%obs_masks_prm, errstr)
          IF (LEN_TRIM(errstr) > 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "Computation of chi2 for a sampled orbit failed: " // TRIM(errstr), 1)
             DO j=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(j))
             END DO
             DEALLOCATE(principal_axes, stat=err)
             DEALLOCATE(debiasing_factor_arr, stat=err)
             DEALLOCATE(debiasing_factor_map, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(cosdec0, stat=err)
             DEALLOCATE(residuals3, stat=err)
             DEALLOCATE(mask_arr2, stat=err)
             DEALLOCATE(failed_flag, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(jac_arr, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms_arr, stat=err)
             DEALLOCATE(diff, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(residuals2, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(imap_arr, stat=err)
             DEALLOCATE(comp_scoords, stat=err)
             DEALLOCATE(partials4, stat=err)
             CALL NULLIFY(orb_global)
             CALL NULLIFY(orb_local)
             CALL NULLIFY(t0)
             RETURN
          END IF

          ! Compute dchi2 wrt best fit orbit
          dchi2 = chi2 - this%chi2_min_prm
          IF (this%dchi2_rejection_prm .AND. &
               dchi2 > this%dchi2_prm) THEN
             ! The dchi2 is used and its value is not acceptable.
             failed_flag(4) = failed_flag(4) + 1
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A)") &
                     "Failed (Delta chi2 value not acceptable)"
             END IF
             elements = getElements(orb_arr(i),element_type)
             IF (elements(indx) < elements_global(indx)) THEN
                ipdf = ipdf + 1
             END IF
             CYCLE
          END IF

          ! Sigma_elements^(-1) = A^T Sigma_obs^(-1) A, where A is the
          ! partial derivatives matrix of ephemerides wrt elements:
          information_matrix_local(:,:) = 0.0_bp
          DO j=1,nobs
             information_matrix_local = information_matrix_local + &
                  MATMUL(MATMUL(TRANSPOSE(partials4(i,1:6,1:6,j)), &
                  information_matrix_obs(j,1:6,1:6)), &
                  partials4(i,1:6,1:6,j))
          END DO
          IF (this%regularization_prm) THEN
             ! Jeffrey's apriori:
             apriori = SQRT(ABS(determinant(information_matrix_local, errstr)))
             IF (LEN_TRIM(errstr) > 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                     "Computation of determinant for local information matrix failed: " //TRIM(errstr), 1)
                DO j=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(j))
                END DO
                DEALLOCATE(principal_axes, stat=err)
                DEALLOCATE(debiasing_factor_arr, stat=err)
                DEALLOCATE(debiasing_factor_map, stat=err)
                DEALLOCATE(obs_coords, stat=err)
                DEALLOCATE(cosdec0, stat=err)
                DEALLOCATE(residuals3, stat=err)
                DEALLOCATE(mask_arr2, stat=err)
                DEALLOCATE(failed_flag, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(jac_arr, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms_arr, stat=err)
                DEALLOCATE(diff, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(residuals2, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(imap_arr, stat=err)
                DEALLOCATE(comp_scoords, stat=err)
                DEALLOCATE(partials4, stat=err)
                CALL NULLIFY(orb_global)
                CALL NULLIFY(orb_local)
                CALL NULLIFY(t0)
                RETURN
             END IF
          ELSE
             apriori = 1.0_bp
          END IF
          ! Probability density function:
          pdf = apriori*EXP(-0.5_bp*(chi2 - SUM(n0(1:6))))*debiasing_factor_map(imap_arr(i))
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(2X,A)") "Sample information matrix:"
             CALL matrix_print(information_matrix_local, stdout, errstr)
             WRITE(stdout,"(2X,A,1X,2"//TRIM(frmt)//")") "Sample chi2:", chi2, chi2-SUM(n0(1:6))
             WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Sample apriori:", apriori
             WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Sample pdf:", pdf
             WRITE(stdout,*)
          END IF

          ! Jaocobians between different elements
          CALL partialsCartesianWrtKeplerian(orb_arr(i), jacobian_matrix, "equatorial")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "TRACE BACK (356) ",4)
             RETURN
          ELSE
             jac_car_kep = ABS(determinant(jacobian_matrix, errstr))
             IF (LEN_TRIM(errstr) > 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                     "Computation of determinant for local " // &
                     "jacobian matrix failed: " // TRIM(errstr), 1)
                DO j=1,SIZE(orb_arr)
                   CALL NULLIFY(orb_arr(j))
                END DO
                DEALLOCATE(principal_axes, stat=err)
                DEALLOCATE(debiasing_factor_arr, stat=err)
                DEALLOCATE(debiasing_factor_map, stat=err)
                DEALLOCATE(obs_coords, stat=err)
                DEALLOCATE(cosdec0, stat=err)
                DEALLOCATE(residuals3, stat=err)
                DEALLOCATE(mask_arr2, stat=err)
                DEALLOCATE(failed_flag, stat=err)
                DEALLOCATE(reg_apriori_arr, stat=err)
                DEALLOCATE(pdf_arr, stat=err)
                DEALLOCATE(jac_arr, stat=err)
                DEALLOCATE(rchi2_arr, stat=err)
                DEALLOCATE(rms_arr, stat=err)
                DEALLOCATE(diff, stat=err)
                DEALLOCATE(obs_scoords, stat=err)
                DEALLOCATE(obsy_ccoords, stat=err)
                DEALLOCATE(residuals2, stat=err)
                DEALLOCATE(information_matrix_obs, stat=err)
                DEALLOCATE(orb_arr, stat=err)
                DEALLOCATE(imap_arr, stat=err)
                DEALLOCATE(comp_scoords, stat=err)
                DEALLOCATE(partials4, stat=err)
                CALL NULLIFY(orb_global)
                CALL NULLIFY(orb_local)
                CALL NULLIFY(t0)
                RETURN
             END IF
          END IF
          elements = getElements(orb_arr(i), "keplerian")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "TRACE BACK (357) ",4)
             !error = .FALSE.
             !jac_equ_kep = -1.0_bp
             RETURN
          ELSE
             jac_equ_kep = 0.5_bp*elements(2)*SIN(0.5_bp*elements(3)) / &
                  COS(0.5_bp*elements(3))**3
          END IF
          !err_verb = err_verb_
          naccepted = naccepted + 1
          iorb = iorb + 1
          orb_accepted(iorb) = copy(orb_arr(i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "TRACE BACK (355)", 1)
             DO j=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(j))
             END DO
             DEALLOCATE(principal_axes, stat=err)
             DEALLOCATE(debiasing_factor_arr, stat=err)
             DEALLOCATE(debiasing_factor_map, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(cosdec0, stat=err)
             DEALLOCATE(residuals3, stat=err)
             DEALLOCATE(mask_arr2, stat=err)
             DEALLOCATE(failed_flag, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(jac_arr, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms_arr, stat=err)
             DEALLOCATE(diff, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(residuals2, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(imap_arr, stat=err)
             DEALLOCATE(comp_scoords, stat=err)
             DEALLOCATE(partials4, stat=err)
             CALL NULLIFY(orb_global)
             CALL NULLIFY(orb_local)
             CALL NULLIFY(t0)
             RETURN
          END IF

          reg_apriori_arr(iorb) = apriori
          pdf_arr(iorb) = pdf
          debiasing_factor_arr(iorb) = debiasing_factor_map(imap_arr(i))
          rchi2_arr(iorb) = chi2 - SUM(n0(1:6))
          jac_arr(iorb,1) = -1.0_bp
          jac_arr(iorb,2) = jac_car_kep
          jac_arr(iorb,3) = jac_equ_kep
          n0_ = n0
          WHERE (n0_ == 0)
             n0_ = 1
          END WHERE
          rms_arr(iorb,:) = SQRT(SUM(residuals3(iorb,:,:)**2.0_bp,dim=1,mask=this%obs_masks_prm)/n0_)
          ! P.d.f.:
          ! - exponential part
          ! - a priori p.d.f. for invariance (optionally)
          IF (iorb == this%vov_norb_prm) THEN
             itrial = itrial - (norb - i)
             EXIT
          END IF
       END DO

       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,2(I0,1X))") "Nr of accepted orbits and trials: ",naccepted,norb
       ENDIF
       IF (naccepted == 0) THEN
          norb = norb*100
       ELSE
          norb = NINT((this%vov_norb_prm - iorb)*(1.2_bp*norb/naccepted))
       END IF
       IF (norb > norb_simult_max) THEN
          norb = norb_simult_max
       END IF
       IF (norb /= SIZE(orb_arr,dim=1)) THEN
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, imap_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / volumeOfVariation", &
                  "Could not deallocate memory (10).", 1)
             DEALLOCATE(principal_axes, stat=err)
             DEALLOCATE(debiasing_factor_arr, stat=err)
             DEALLOCATE(debiasing_factor_map, stat=err)
             DEALLOCATE(obs_coords, stat=err)
             DEALLOCATE(cosdec0, stat=err)
             DEALLOCATE(residuals3, stat=err)
             DEALLOCATE(mask_arr2, stat=err)
             DEALLOCATE(failed_flag, stat=err)
             DEALLOCATE(reg_apriori_arr, stat=err)
             DEALLOCATE(pdf_arr, stat=err)
             DEALLOCATE(jac_arr, stat=err)
             DEALLOCATE(rchi2_arr, stat=err)
             DEALLOCATE(rms_arr, stat=err)
             DEALLOCATE(diff, stat=err)
             DEALLOCATE(obs_scoords, stat=err)
             DEALLOCATE(obsy_ccoords, stat=err)
             DEALLOCATE(residuals2, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(orb_arr, stat=err)
             DEALLOCATE(imap_arr, stat=err)
             DEALLOCATE(comp_scoords, stat=err)
             DEALLOCATE(partials4, stat=err)
             CALL NULLIFY(orb_global)
             CALL NULLIFY(orb_local)
             CALL NULLIFY(t0)
             RETURN
          END IF
       END IF
       DEALLOCATE(comp_scoords, partials4, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Could not deallocate memory (15).", 1)
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(imap_arr, stat=err)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials4, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF

    END DO vov_main

    this%vov_ntrial_cmp = itrial
    this%vov_norb_cmp = iorb
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "Final number of orbits and the required trials:"
       WRITE(stdout,"(2X,2(I0,2X))") iorb, itrial
       WRITE(stdout,*) " Total failure percentage (1), and failure " // &
            "due to (2) eccentricity, (3) inclination, " // &
            "(4) residuals, and (5) pdf:"
       nfailed = SUM(failed_flag)
       nfailed = MAX(1,nfailed)
       WRITE(stdout,"(2X,5"//TRIM(frmt)//")") &
            100.0_lp*REAL(SUM(failed_flag),lp)/itrial, &
            100.0_lp*REAL(failed_flag(1),lp)/nfailed, &
            100.0_lp*REAL(failed_flag(2),lp)/nfailed, &
            100.0_lp*REAL(failed_flag(3),lp)/nfailed, &
            100.0_lp*REAL(failed_flag(4),lp)/nfailed
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DO i=1,SIZE(this%orb_arr_cmp)
          CALL NULLIFY(this%orb_arr_cmp(i))
       END DO
       DEALLOCATE(this%orb_arr_cmp, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "Could not deallocate memory (20).", 1)
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          DEALLOCATE(imap_arr, stat=err)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials4, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN       
       END IF
    END IF
    ALLOCATE(this%orb_arr_cmp(iorb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "Could not allocate memory (15).", 1)
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(imap_arr, stat=err)
       DEALLOCATE(comp_scoords, stat=err)
       DEALLOCATE(partials4, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(orb_local)
       CALL NULLIFY(t0)
       RETURN       
    END IF
    IF (iorb == 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "No sample orbit found.", 1)
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(orb_arr, stat=err)
       DEALLOCATE(imap_arr, stat=err)
       DEALLOCATE(comp_scoords, stat=err)
       DEALLOCATE(partials4, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(orb_local)
       CALL NULLIFY(t0)
       RETURN
    END IF

    ! Output: sample orbits, mapping orbits and local intervals of variation (for plotting)
    DO i=1,iorb
       this%orb_arr_cmp(i) = copy(orb_accepted(i))
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / volumeOfVariation", &
               "TRACE BACK (365)", 1)
          IF (ASSOCIATED(orb_arr)) THEN
             DO j=1,SIZE(orb_arr)
                CALL NULLIFY(orb_arr(j))
             END DO
             DEALLOCATE(orb_arr, stat=err)
          END IF
          DEALLOCATE(principal_axes, stat=err)
          DEALLOCATE(debiasing_factor_arr, stat=err)
          DEALLOCATE(debiasing_factor_map, stat=err)
          DEALLOCATE(obs_coords, stat=err)
          DEALLOCATE(cosdec0, stat=err)
          DEALLOCATE(residuals3, stat=err)
          DEALLOCATE(mask_arr2, stat=err)
          DEALLOCATE(failed_flag, stat=err)
          DEALLOCATE(reg_apriori_arr, stat=err)
          DEALLOCATE(pdf_arr, stat=err)
          DEALLOCATE(jac_arr, stat=err)
          DEALLOCATE(rchi2_arr, stat=err)
          DEALLOCATE(rms_arr, stat=err)
          DEALLOCATE(diff, stat=err)
          DEALLOCATE(obs_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(residuals2, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(imap_arr, stat=err)
          DEALLOCATE(comp_scoords, stat=err)
          DEALLOCATE(partials4, stat=err)
          CALL NULLIFY(orb_global)
          CALL NULLIFY(orb_local)
          CALL NULLIFY(t0)
          RETURN
       END IF
    END DO
    CALL propagate(this%orb_arr_cmp, t0)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (370)", 1)
       IF (ASSOCIATED(orb_arr)) THEN
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
       END IF
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(imap_arr, stat=err)
       DEALLOCATE(comp_scoords, stat=err)
       DEALLOCATE(partials4, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(orb_local)
       CALL NULLIFY(t0)
       RETURN
    END IF
    ! Maximum likelihood point:
    CALL NULLIFY(this%orb_ml_cmp)
    this%orb_ml_cmp = copy(orb_global)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / volumeOfVariation", &
            "TRACE BACK (375)", 1)
       IF (ASSOCIATED(orb_arr)) THEN
          DO i=1,SIZE(orb_arr)
             CALL NULLIFY(orb_arr(i))
          END DO
          DEALLOCATE(orb_arr, stat=err)
       END IF
       DEALLOCATE(principal_axes, stat=err)
       DEALLOCATE(debiasing_factor_arr, stat=err)
       DEALLOCATE(debiasing_factor_map, stat=err)
       DEALLOCATE(obs_coords, stat=err)
       DEALLOCATE(cosdec0, stat=err)
       DEALLOCATE(residuals3, stat=err)
       DEALLOCATE(mask_arr2, stat=err)
       DEALLOCATE(failed_flag, stat=err)
       DEALLOCATE(reg_apriori_arr, stat=err)
       DEALLOCATE(pdf_arr, stat=err)
       DEALLOCATE(jac_arr, stat=err)
       DEALLOCATE(rchi2_arr, stat=err)
       DEALLOCATE(rms_arr, stat=err)
       DEALLOCATE(diff, stat=err)
       DEALLOCATE(obs_scoords, stat=err)
       DEALLOCATE(obsy_ccoords, stat=err)
       DEALLOCATE(residuals2, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(imap_arr, stat=err)
       DEALLOCATE(comp_scoords, stat=err)
       DEALLOCATE(partials4, stat=err)
       CALL NULLIFY(orb_global)
       CALL NULLIFY(orb_local)
       CALL NULLIFY(t0)
       RETURN
    END IF
    CALL NULLIFY(orb_global)
    IF (.NOT.ASSOCIATED(this%cov_ml_cmp)) THEN
       ALLOCATE(this%cov_ml_cmp(6,6), stat=err)
    END IF
    this%cov_ml_cmp = covariance_global
    this%cov_type_prm = this%element_type_prm
    this%ls_elem_mask_prm = element_mask
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp, stat=err)
    END IF
    ALLOCATE(this%pdf_arr_cmp(iorb), stat=err)
    this%pdf_arr_cmp = pdf_arr(1:iorb)
    IF (ASSOCIATED(this%reg_apr_arr_cmp)) THEN
       DEALLOCATE(this%reg_apr_arr_cmp, stat=err)
    END IF
    ALLOCATE(this%reg_apr_arr_cmp(iorb), stat=err)
    this%reg_apr_arr_cmp = reg_apriori_arr(1:iorb)
    IF (ASSOCIATED(this%jac_arr_cmp)) THEN
       DEALLOCATE(this%jac_arr_cmp, stat=err)
    END IF
    ALLOCATE(this%jac_arr_cmp(iorb,3), stat=err)
    this%jac_arr_cmp = jac_arr(1:iorb,1:3)
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp, stat=err)
    END IF
    ALLOCATE(this%rchi2_arr_cmp(iorb), stat=err)
    this%rchi2_arr_cmp = rchi2_arr(1:iorb)
    this%chi2_min_cmp = MINVAL(this%rchi2_arr_cmp) + SUM(n0(1:6))
    IF (ASSOCIATED(this%rms_arr_cmp)) THEN
       DEALLOCATE(this%rms_arr_cmp, stat=err)
    END IF
    ALLOCATE(this%rms_arr_cmp(iorb,6), stat=err)
    this%rms_arr_cmp = rms_arr(1:iorb,1:6)
    IF (ASSOCIATED(this%vov_map_cmp)) THEN
       DEALLOCATE(this%vov_map_cmp, stat=err)
    END IF
    ALLOCATE(this%vov_map_cmp(this%vov_nmap_prm,12), stat=err)
    DO i=1,this%vov_nmap_prm
       this%vov_map_cmp(i,1:6) = elements_local_arr(i,1:6)
       DO j=1,6
          this%vov_map_cmp(i,6+j) = 0.0_bp
          IF (j == indx) THEN
             this%vov_map_cmp(i,6+j) = stdev_global(j)
             CYCLE
          END IF

          DO k=1,6
             IF (principal_axes(i,j,k) > 0.0_bp) THEN
                this%vov_map_cmp(i,6+j) = this%vov_map_cmp(i,6+j) + principal_axes(i,j,k)
             ELSE
                this%vov_map_cmp(i,6+j) = this%vov_map_cmp(i,6+j) - principal_axes(i,j,k)                
             END IF
          END DO
       END DO
    END DO
    this%vov_scaling_cmp = this%vov_scaling_prm

    IF (ASSOCIATED(orb_arr)) THEN
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(orb_arr, stat=err)
    END IF
    DO i=1,SIZE(orb_accepted)
       CALL NULLIFY(orb_accepted(i))
    END DO
    DEALLOCATE(principal_axes, stat=err)
    DEALLOCATE(debiasing_factor_arr, stat=err)
    DEALLOCATE(debiasing_factor_map, stat=err)
    DEALLOCATE(obs_coords, stat=err)
    DEALLOCATE(cosdec0, stat=err)
    DEALLOCATE(residuals3, stat=err)
    DEALLOCATE(mask_arr2, stat=err)
    DEALLOCATE(failed_flag, stat=err)
    DEALLOCATE(reg_apriori_arr, stat=err)
    DEALLOCATE(pdf_arr, stat=err)
    DEALLOCATE(jac_arr, stat=err)
    DEALLOCATE(rchi2_arr, stat=err)
    DEALLOCATE(rms_arr, stat=err)
    DEALLOCATE(diff, stat=err)
    DEALLOCATE(obs_scoords, stat=err)
    DEALLOCATE(obsy_ccoords, stat=err)
    DEALLOCATE(residuals2, stat=err)
    DEALLOCATE(information_matrix_obs, stat=err)
    DEALLOCATE(imap_arr, stat=err)
    DEALLOCATE(comp_scoords, stat=err)
    DEALLOCATE(partials4, stat=err)

  END SUBROUTINE volumeOfVariation





  SUBROUTINE virtualObservationMCMC(this, orb_arr_in)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this
    TYPE (Orbit), DIMENSION(:), INTENT(in) :: orb_arr_in

    TYPE (Time) :: &
         t0
    TYPE (Orbit), DIMENSION(:), POINTER :: &
         orb_arr => NULL()
    TYPE (Orbit) :: &
         orb
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=DYN_MODEL_LEN) :: &
         storb_dyn_model, &
         orb_dyn_model
    CHARACTER(len=INTEGRATOR_LEN) :: &
         storb_integrator, &
         orb_integrator
    CHARACTER(len=64) :: &
         frmt = "(F25.15,1X)", &
         efrmt = "(E10.4,1X)"
    CHARACTER(len=64) :: &
         str
    REAL(bp), DIMENSION(:,:), POINTER :: &
         elements_arr => NULL()
    REAL(bp), DIMENSION(3) :: &
         ran_arr
    REAL(bp), DIMENSION(6) :: &
         elements, elements_
    REAL(bp) :: &
         a, &
         a_r, &
         pdf, pdf_, &
         rchi2, rchi2_, &
         storb_integration_step, &
         q, &
         orb_integration_step
    INTEGER :: &
         i, &
         j, &
         err, &
         norb
    LOGICAL, DIMENSION(:,:), ALLOCATABLE :: &
         obs_masks_
    LOGICAL, DIMENSION(:), ALLOCATABLE :: &
         mask_arr
    LOGICAL, DIMENSION(6) :: &
         mask
    LOGICAL :: &
         accept, &
         burn_in = .FALSE., &
         first, &
         outlier_rejection_, &
         parameters_agree

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (getNrOfObservations(this%obss) < 3) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "Less than two observations available.", 1)
       RETURN
    END IF
    IF (.NOT.ASSOCIATED(this%obs_masks_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "Observation masks not set.", 1)
       RETURN
    END IF
    ALLOCATE(obs_masks_(SIZE(this%obs_masks_prm,dim=1),SIZE(this%obs_masks_prm,dim=2)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / " // &
            "observationSampling", &
            "Could not allocate memory.", 1)
       DEALLOCATE(obs_masks_, stat=err)
       RETURN
    END IF
    obs_masks_ = this%obs_masks_prm

    IF (.NOT.ASSOCIATED(this%res_accept_prm)) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "Window for accepted residuals not set.", 1)
       RETURN
    END IF
    IF (this%vomcmc_norb_prm < 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "Required number of sample orbits not set.", 1)
       RETURN
    END IF
!!$    CALL comparePropagationParameters(this, orb_arr(1), &
!!$         parameters_agree=parameters_agree)
!!$    IF (error .OR. .NOT.parameters_agree) THEN
!!$       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
!!$            "TRACE BACK (1)", 1)
!!$       RETURN
!!$    END IF

    frame = getFrame(orb_arr_in(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "TRACE BACK (70)", 1)
       RETURN
    END IF
    ! Inversion epoch is bound to the epoch of the first input orbit:
    t0 = getTime(orb_arr_in(1))
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "TRACE BACK (65)", 1)
       RETURN
    END IF

    !!
    !! 1) EXTRACT OR GENERATE ORBITAL ELEMENTS THAT WILL FORM THE
    !! DISCRETE PROPOSAL
    !!

    ! Generate orbits focusing on the relevant phase-space
    ! regime. Also set this%chi2_min_prm to the chi2 corresponding to
    ! the best fit orbit to the nominal set of obserations:
    IF (.FALSE.) THEN
       CALL observationSampling(this, orb_arr_in)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
               "TRACE BACK (66)", 1)       
          RETURN
       END IF
       ALLOCATE(orb_arr(SIZE(this%orb_arr_cmp)))
       DO i=1,SIZE(this%orb_arr_cmp)
          orb_arr(i) = copy(this%orb_arr_cmp(i))
       END DO
    ELSE
       ALLOCATE(orb_arr(SIZE(orb_arr_in)))
       DO i=1,SIZE(orb_arr_in)
          orb_arr(i) = copy(orb_arr_in(i))
       END DO
    END IF
    norb = SIZE(orb_arr)
    ALLOCATE(elements_arr(norb,6))
    DO i=1,norb
       elements_arr(i,:) = getElements(orb_arr(i), this%element_type_prm)
    END DO


    !!
    !! 2) MCMC USING DISCRETE PROPOSAL
    !!

    CALL setParameters(this, outlier_rejection=outlier_rejection_)
    ! Generate a cumulative distribution for the mapping parameter
    ! based on the local sampling volumes:
    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       DEALLOCATE(this%orb_arr_cmp)
    END IF
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp)
    END IF
    IF (ASSOCIATED(this%pdf_arr_cmp)) THEN
       DEALLOCATE(this%pdf_arr_cmp)
    END IF
    IF (ASSOCIATED(this%repetition_arr_cmp)) THEN
       DEALLOCATE(this%repetition_arr_cmp)
    END IF
    ALLOCATE( &
         this%orb_arr_cmp(this%vomcmc_norb_prm), &
         this%rchi2_arr_cmp(this%vomcmc_norb_prm), &
         this%pdf_arr_cmp(this%vomcmc_norb_prm), &
         this%repetition_arr_cmp(this%vomcmc_norb_prm))
    this%repetition_arr_cmp = 0

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "DMCMC SAMPLING ..."
    END IF
    this%vomcmc_norb_cmp = 0
    this%vomcmc_ntrial_cmp = -1
    DO WHILE (this%vomcmc_norb_cmp < this%vomcmc_norb_prm .AND. &
         this%vomcmc_ntrial_cmp < this%vomcmc_ntrial_prm)

       this%vomcmc_ntrial_cmp = this%vomcmc_ntrial_cmp + 1

       !!
       !! 3) ORBIT GENERATION
       !! 

       ! elements_ has the last accepted set of elements
       IF (this%vomcmc_ntrial_cmp == 0) THEN
          elements_ = elements_arr(1,:)
          elements = elements_
       ELSE
          CALL randomNumber(ran_arr)
          i = CEILING(ran_arr(1)*norb)
          j = CEILING(ran_arr(2)*norb)
          IF (i == j) THEN
             CYCLE
          END IF
          ! Use orbital-element differences as a discrete, symmetric
          ! proposal:
          elements = elements_ + (elements_arr(i,:)-elements_arr(j,:))
          ! Make sure keplerian and cometary orbital elements make
          ! sense:
          IF (this%element_type_prm == "keplerian" .OR. &
               this%element_type_prm == "cometary") THEN
             WHERE (elements(1:3) < 0.0_bp)
                elements(1:3) = EPSILON(elements(1))
             END WHERE
             IF (elements(3) > pi) THEN
                elements(3) = pi - EPSILON(elements(3))
             END IF
             elements(4:5) = MODULO(elements(4:5),two_pi)
             IF (this%element_type_prm == "keplerian") THEN
                IF (elements(2) >= 1.0_bp) THEN
                   elements(2) = 1.0_bp - EPSILON(elements(2))
                END IF
                elements(6) = MODULO(elements(6),two_pi)
             END IF
          END IF

       END IF

       CALL NULLIFY(orb)
       CALL NEW(orb, elements, this%element_type_prm, frame, t0)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
               "TRACE BACK (315)", 1)
          RETURN
       END IF
       CALL setParameters(orb, &
            perturbers=this%perturbers_prm, &
            asteroid_perturbers=this%ast_perturbers_prm, &
            dyn_model=this%dyn_model_prm, &
            integration_step=this%integration_step_prm, &
            integrator=this%integrator_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
               "TRACE BACK (325)", 1)
          RETURN
       END IF

       !!
       !! 4) ACCEPTANCE / REJECTION OF GENERATED ORBIT
       !!

       ! Informative priors
       IF (this%informative_apriori_prm .AND. this%vomcmc_ntrial_cmp > 0) THEN

          ! Semimajor axis:
          a = -1.0_bp
          IF (this%apriori_a_min_prm >= 0.0_bp .OR. &
               this%apriori_a_max_prm >= 0.0_bp) THEN
             a = getSemimajorAxis(orb)
             IF (this%apriori_a_min_prm >= 0.0_bp .AND. &
                  a < this%apriori_a_min_prm) THEN
                ! Semimajor axis too small
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (semimajor axis too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             IF (this%apriori_a_max_prm >= 0.0_bp .AND. &
                  a > this%apriori_a_max_prm) THEN
                ! Semimajor axis too large
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (semimajor axis too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
          ! Periapsis distance:
          IF (this%apriori_periapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_periapsis_max_prm >= 0.0_bp) THEN
             IF (a >= 0.0_bp) THEN
                CALL getPeriapsisDistance(orb, q, a)
             ELSE
                CALL getPeriapsisDistance(orb, q)
             END IF
             IF (error) THEN
                error = .FALSE.
                CYCLE
             END IF
             ! Periapsis distance too small:
             IF (this%apriori_periapsis_min_prm >= 0.0_bp .AND. &
                  q < this%apriori_periapsis_min_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (periapsis distance too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             ! Periapsis distance too large:
             IF (this%apriori_periapsis_max_prm >= 0.0_bp .AND. &
                  q > this%apriori_periapsis_max_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (periapsis distance too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
          ! Apoapsis distance:
          IF (this%apriori_apoapsis_min_prm >= 0.0_bp .OR. &
               this%apriori_apoapsis_max_prm >= 0.0_bp) THEN
             IF (a >= 0.0_bp) THEN
                CALL getApoapsisDistance(orb, Q, a)
             ELSE
                CALL getApoapsisDistance(orb, Q)
             END IF
             ! Apoapsis distance too small:
             IF (this%apriori_apoapsis_min_prm >= 0.0_bp .AND. &
                  Q < this%apriori_apoapsis_min_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (apoapsis distance too small: ", a, " au)"
                END IF
                CYCLE
             END IF
             ! Apoapsis distance too large:
             IF (this%apriori_apoapsis_max_prm >= 0.0_bp .AND. &
                  Q > this%apriori_apoapsis_max_prm) THEN
                IF (info_verb >= 4) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (apoapsis distance too large: ", a, " au)"
                END IF
                CYCLE
             END IF
          END IF
       END IF

       ! Compute "reduced" chi2:
       rchi2 = getChi2(this, orb) - COUNT(obs_masks_)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
               "TRACE BACK (330)", 1)
          RETURN
       END IF

       ! Probability density function value:
       pdf = EXP(-0.5_bp*(rchi2))
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,1X,2"//TRIM(frmt)//")") "Sample reduced chi2:", rchi2
          WRITE(stdout,"(2X,A,1X,1"//TRIM(efrmt)//")") "Sample pdf:", pdf
          WRITE(stdout,*)
       END IF

       ! Set initial values and cycle if this is the first loop
       IF (this%vomcmc_ntrial_cmp == 0) THEN
          rchi2_ = rchi2
          pdf_ = pdf
          CYCLE
       END IF

       ! Compute the pdv ratio between the previous accepted orbit and
       ! the current trial orbit and use the MH criterion to decide
       ! whether the trial orbit should be accepted or rejected
       a_r = pdf/pdf_

       IF (a_r > HUGE(a_r)) THEN
          WRITE(stdout,*) "huge avalue", MIN(rchi2,rchi2_)
          a_r = rchi2_/rchi2
          burn_in = .TRUE.
       END IF

       accept = .FALSE.
       IF (a_r > ran_arr(3)) THEN
          accept = .TRUE.
       END IF

       IF (accept) THEN
          IF (.NOT.burn_in) THEN
             this%vomcmc_norb_cmp = this%vomcmc_norb_cmp + 1
             this%orb_arr_cmp(this%vomcmc_norb_cmp) = copy(orb)
             this%rchi2_arr_cmp(this%vomcmc_norb_cmp) = rchi2
             this%pdf_arr_cmp(this%vomcmc_norb_cmp) = pdf
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
                     "TRACE BACK (355)", 1)
                RETURN
             END IF
          ELSE
             burn_in = .FALSE.
          END IF
          elements_ = elements
          rchi2_ = rchi2
          pdf_ = pdf
       END IF
       IF (this%vomcmc_norb_cmp > 0) THEN
          this%repetition_arr_cmp(this%vomcmc_norb_cmp) = &
               this%repetition_arr_cmp(this%vomcmc_norb_cmp) + 1
       END IF

       IF (info_verb >= 3) THEN
          IF (this%element_type_prm == "cometary") THEN
             elements(3:5) = elements(3:5)/rad_deg
          END IF
          IF (accept) THEN
             WRITE(stdout,"(A,1X,9(F15.8,1X),A,1X,I0)") &
                  TRIM(this%element_type_prm), &
                  elements, pdf, rchi2, a_r, &
                  "ACCEPTED", this%vomcmc_norb_cmp + 1
          ELSE
             WRITE(stdout,"(A,1X,9(F15.8,1X),A,1X,I0)") &
                  TRIM(this%element_type_prm), &
                  elements, pdf, rchi2, a_r, "REJECTED", &
                  this%vomcmc_ntrial_cmp - this%vomcmc_norb_cmp
          END IF
       END IF

       IF (info_verb >= 2 .AND. MOD(this%vomcmc_ntrial_cmp,500) == 0) THEN
          WRITE(stdout,"(2X,A,3(I0,1X))") "Nr of accepted orbits and trials: ", &
               this%vomcmc_norb_cmp, this%vomcmc_ntrial_cmp
       END IF

    END DO

    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2X,A)") "Final number of orbits and the required trials:"
       WRITE(stdout,"(2X,2(I0,2X))") this%vomcmc_norb_cmp, this%vomcmc_ntrial_cmp
    END IF
    IF (this%vomcmc_norb_cmp == 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / virtualObservationMCMC", &
            "No sample orbit found.", 1)
       RETURN
    END IF

    this%orb_arr_cmp => &
         reallocate(this%orb_arr_cmp, this%vomcmc_norb_cmp)
    this%rchi2_arr_cmp => &
         reallocate(this%rchi2_arr_cmp, this%vomcmc_norb_cmp)
    this%repetition_arr_cmp => &
         reallocate(this%repetition_arr_cmp, this%vomcmc_norb_cmp)
    this%pdf_arr_cmp => &
         reallocate(this%pdf_arr_cmp, this%vomcmc_norb_cmp)
    ! Compute pdf based on repetitions
    this%pdf_arr_cmp = REAL(this%repetition_arr_cmp,bp)/SUM(this%repetition_arr_cmp)
    CALL NULLIFY(this%orb_ml_cmp)
    i = MAXLOC(this%pdf_arr_cmp,dim=1)
    this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
    this%chi2_min_cmp = MINVAL(this%rchi2_arr_cmp) + COUNT(obs_masks_)

    CALL NULLIFY(orb)
    DO i=1,SIZE(orb_arr)
       CALL NULLIFY(orb_arr(i))
    END DO
    DEALLOCATE(orb_arr, elements_arr, stat=err)

  END SUBROUTINE virtualObservationMCMC





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

    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
         obs_scoords => NULL(), &
         ephemerides => NULL()
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoords => NULL()
    TYPE (Time) :: t, t_
    TYPE (Orbit) :: orb
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: &
         obsy_codes => NULL()
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type, element_type_ls
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=DYN_MODEL_LEN) :: dyn_model, dyn_model_
    CHARACTER(len=INTEGRATOR_LEN) :: orb_integrator
    CHARACTER(len=32) :: str
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         partials_arr => NULL(), &
         cov_mat_obs => NULL(), &
         inform_mat_obs_bd => NULL()
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE :: design_mat ! (data,obs,parameter)
    REAL(bp), DIMENSION(:,:), POINTER :: &
         stdev_arr_obs => NULL(), &
         elements_iter_arr => NULL(), &
         rmss_iter_arr => NULL(), &
         orb_additional_perturbers => NULL()
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
               finite_diff=orb_finite_diff, &
               additional_perturbers=orb_additional_perturbers)
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
               asteroid_perturbers=this%ast_perturbers_prm, &
               integration_step=orb_integration_step, &
               integrator=orb_integrator, &
               finite_diff=orb_finite_diff, &
               additional_perturbers=orb_additional_perturbers)
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
  !! Keplerian) and their covariances in the two-body or full
  !! many-body approaches using the Levenberg-Marquardt algorithm.
  !!
  !! Returns error.
  !!
  SUBROUTINE levenbergMarquardt_SO(this, preliminary_orbit)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), INTENT(in)              :: preliminary_orbit    

    INTEGER, PARAMETER :: niter_size = 3
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
         obs_scoords => NULL()
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoords => NULL()
    TYPE (Time) :: t, t_
    TYPE (Orbit) :: orb
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: &
         obsy_codes => NULL()
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type_
    CHARACTER(len=FRAME_LEN) :: frame_
    CHARACTER(len=DYN_MODEL_LEN) :: dyn_model_
    CHARACTER(len=INTEGRATOR_LEN) :: integrator_
    CHARACTER(len=32) :: str
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         information_matrix_measur => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: &
         stdev_arr_measur => NULL()
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
             WRITE(stdout,"(2X,I0,A,I0,A)") j, "th minor and ", i, &
                  "th major call to 'LevenbergMarquardt_private'" // &
                  " from loop in 'levenbergMarquardt_SO'."
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
             EXIT
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
       ELSE ! IF (.NOT.this%outlier_rejection_prm) THEN
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
         asteroid_perturbers=this%ast_perturbers_prm, &
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

    IF (rchi2 > this%ls_rchi2_acceptable_prm) THEN
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
      TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
           ephemerides => NULL()
      REAL(bp), DIMENSION(:,:,:), POINTER :: &
           partials_arr => NULL()
      REAL(bp), DIMENSION(SIZE(residuals,dim=1),SIZE(residuals,dim=2)) :: computed
      REAL(bp) :: chi2
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
           asteroid_perturbers=this%ast_perturbers_prm, &
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
               t_ = getTime(obsy_ccoords(j))
               WRITE(stdout,"(2X,A,2(F20.7,1X),A,1X,A)") &
                    " ", &
                    residuals(j,2:3)/rad_asec, " ", &
                    TRIM(obsy_codes(j))
            ELSE
               WRITE(stdout,"(2X,A,2(F15.7,1X),A,1X,A)") "(", &
                    residuals(j,2:3)/rad_asec, ")", TRIM(obsy_codes(j))
            END IF
         END IF
      END DO

      ! Compute chi2:
      chi2 = chi_square(residuals, information_matrix_measur, mask_measur, errstr)

      ! Compute reduced chi2:
      rchi2 = chi2 / REAL(COUNT(mask_measur)-COUNT(mask_param),bp)

      IF (info_verb >= 2) THEN
         WRITE(stdout,"(2X,A,2(1X,F15.7))") "RMS RA & Dec [arcsec]: ", &
              SQRT(SUM(residuals(:,2)**2,mask=mask_measur(:,2))/COUNT(mask_measur(:,2)))/rad_asec, &
              SQRT(SUM(residuals(:,3)**2,mask=mask_measur(:,3))/COUNT(mask_measur(:,3)))/rad_asec
         WRITE(stdout,"(2X,A,1X,F20.5)") "Chi2: ", &
              chi2
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

    REAL(bp), DIMENSION(:,:,:), POINTER :: jacobians => NULL()
    REAL(bp), DIMENSION(:), POINTER     :: pdf_arr => NULL()
    REAL(bp), DIMENSION(6,6)            :: jacobian, cov
    REAL(bp)                            :: det
    INTEGER                             :: i, err

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
       pdf_arr => getDiscretePDF(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / propagate", &
               "TRACE BACK (10)", 1)
          DEALLOCATE(jacobians, stat=err)
          RETURN
       END IF
       DO i=1,SIZE(this%orb_arr_cmp)
          det = determinant(jacobians(i,:,:), errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / propagate", &
                  "TRACE BACK (15) " // TRIM(errstr), 1)
             errstr = ""
             DEALLOCATE(jacobians, stat=err)
             RETURN
          END IF
          ! Changed 2008-12-14
          !this%pdf_arr_cmp(i) = pdf_arr(i) / det
          !this%jac_arr_cmp(i,1) = this%jac_arr_cmp(i,1) / det
          this%pdf_arr_cmp(i) = pdf_arr(i) * det
          this%jac_arr_cmp(i,1) = this%jac_arr_cmp(i,1) * det
       END DO
       DEALLOCATE(jacobians, pdf_arr, stat=err)
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
               "TRACE BACK (20)", 1)
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

    REAL(bp), DIMENSION(:,:), POINTER :: stdevs => NULL()
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
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
            "TRACE BACK 5", 1)
       RETURN       
    END IF
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
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / setAcceptanceWindow", &
            "TRACE BACK 10", 1)
       RETURN
    END IF
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

    REAL(bp), DIMENSION(:,:), POINTER :: stdevs => NULL()
    REAL(bp) :: mean
    INTEGER :: i, j, nobs, err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setGenerationWindow", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (this%generat_multiplier_prm < 0.0_bp) THEN
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
             CALL moments(this%res_arr_cmp(:,i,j), mean=mean, errstr=errstr)
             IF (LEN_TRIM(errstr) /= 0) THEN
                error = .TRUE.
                CALL errorMessage("StochasticOrbit / setGenerationWindow", &
                     "Could not compute moments. " // TRIM(errstr), 1)
                errstr = ""
                RETURN
             END IF
             this%sor_deviates_prm(i,j,1) = mean
             !this%sor_deviates_prm(i,j,2) = this%generat_multiplier_prm*stdev
          END DO
       END DO
    ELSE
       stdevs => getStandardDeviations(this%obss)
       this%sor_deviates_prm(1:nobs,1:6,1) = 0.0_bp          
       IF (PRESENT(offset)) THEN
          this%sor_deviates_prm(1,2:3,1) = offset(1:2)
          this%sor_deviates_prm(nobs,2:3,1) = offset(3:4)
       END IF
       this%sor_deviates_prm(1:nobs,1:6,2) = this%generat_multiplier_prm*stdevs(1:nobs,1:6)
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
  SUBROUTINE setParameters_SO(this, orb_ml, &
       dyn_model, perturbers, asteroid_perturbers, integration_step, integrator, &
       finite_diff, &
       t_inv, element_type, multiple_objects, outlier_rejection, &
       dchi2_rejection, dchi2_max, regularized_pdf, jacobians_pdf, &
       accept_multiplier, outlier_multiplier, &
       gaussian_pdf, chi2_min, chi2_min_init, &
       apriori_a_max, apriori_a_min, apriori_periapsis_max, apriori_periapsis_min, &
       apriori_apoapsis_max, apriori_apoapsis_min, apriori_rho_max, apriori_rho_min, &
       sor_norb, sor_norb_sw, sor_ntrial, sor_ntrial_sw, &
       sor_rho1_l, sor_rho1_u, sor_rho2_l, sor_rho2_u, &
       sor_random_obs_selection, sor_niter, generat_multiplier, &
       sor_generat_offset, sor_2point_method, sor_2point_method_sw, &
       sor_iterate_bounds, &
       vov_norb, vov_ntrial, vov_norb_iter, vov_ntrial_iter, &
       vov_nmap, vov_niter, vov_scaling, vov_mapping_mask, &
       vomcmc_norb, vomcmc_ntrial, vomcmc_norb_iter, vomcmc_ntrial_iter, &
       vomcmc_nmap, vomcmc_niter, vomcmc_scaling, vomcmc_mapping_mask, &
       ls_correction_factor, ls_niter_major_max, ls_niter_major_min, ls_niter_minor, &
       ls_element_mask, ls_rchi2_acceptable, &
       cos_nsigma, cos_norb, cos_ntrial, cos_gaussian, &
       smplx_tol, smplx_niter, smplx_force, smplx_similarity_tol, &
       os_norb, os_ntrial, os_sampling_type, generat_gaussian_deviates, &
       set_acceptance_window, center)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Time), INTENT(in), OPTIONAL :: t_inv
    TYPE (Orbit), INTENT(in), OPTIONAL :: orb_ml
    CHARACTER(len=*), INTENT(in), OPTIONAL :: &
         dyn_model, &
         integrator, &
         element_type, &
         sor_2point_method, &
         sor_2point_method_sw
    REAL(bp), DIMENSION(6,2), INTENT(in), OPTIONAL :: &
         vov_scaling, &
         vomcmc_scaling
    REAL(bp), DIMENSION(6), INTENT(in), OPTIONAL :: &
         finite_diff
    REAL(bp), DIMENSION(4), INTENT(in), OPTIONAL :: &
         sor_generat_offset
    REAL(bp), INTENT(in), OPTIONAL :: &
         integration_step, &
         accept_multiplier, &
         outlier_multiplier, &
         chi2_min, &
         chi2_min_init, &
         dchi2_max, &
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
         generat_multiplier, &
         ls_correction_factor, &
         ls_rchi2_acceptable, &
         cos_nsigma, &
         smplx_tol, &
         smplx_similarity_tol
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
         vomcmc_norb, &
         vomcmc_ntrial, &
         vomcmc_norb_iter, &
         vomcmc_ntrial_iter, &
         vomcmc_nmap, &
         vomcmc_niter, &
         ls_niter_major_max, &
         ls_niter_major_min, &
         ls_niter_minor, &
         cos_norb, &
         cos_ntrial, &
         smplx_niter, &
         os_norb, &
         os_ntrial, &
         os_sampling_type, &
         center
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL :: &
         perturbers, &
         sor_iterate_bounds, &
         vov_mapping_mask, &
         vomcmc_mapping_mask, &
         ls_element_mask
    LOGICAL, INTENT(in), OPTIONAL :: &
         multiple_objects, &
         regularized_pdf, &
         dchi2_rejection, &
         jacobians_pdf, &
         sor_random_obs_selection, &
         gaussian_pdf, &
         outlier_rejection, &
         cos_gaussian, &
         smplx_force, &
         generat_gaussian_deviates, &
         set_acceptance_window, &
         asteroid_perturbers

    CHARACTER(len=256) :: str
    INTEGER :: i, err

    IF (.NOT.this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / setParameters", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(orb_ml)) THEN
       this%orb_ml_prm=copy(orb_ml)
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
    IF (PRESENT(asteroid_perturbers)) THEN
       this%ast_perturbers_prm = asteroid_perturbers
       IF (ASSOCIATED(this%orb_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL setParameters(this%orb_arr_cmp(i), asteroid_perturbers=this%ast_perturbers_prm)
          END DO
       END IF
       IF (exist(this%orb_ml_cmp)) THEN
          CALL setParameters(this%orb_ml_cmp, asteroid_perturbers=this%ast_perturbers_prm)
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
    IF (PRESENT(dchi2_rejection)) THEN
       this%dchi2_rejection_prm = dchi2_rejection
    END IF
    IF (PRESENT(dchi2_max)) THEN
       this%dchi2_prm = dchi2_max
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
       IF (PRESENT(set_acceptance_window)) THEN
          ! If set_acceptance_window present, then set acceptance
          ! window only if set_acceptance_window=.true.
          IF (set_acceptance_window) THEN
             CALL setAcceptanceWindow(this)
          END IF
       ELSE
          CALL setAcceptanceWindow(this)          
       END IF
    END IF
    IF (PRESENT(chi2_min)) THEN
       this%chi2_min_prm = chi2_min
    END IF
    IF (PRESENT(chi2_min_init)) THEN
       this%chi2_min_init_prm = chi2_min_init
    END IF
    IF (PRESENT(center)) THEN
       this%center_prm = center
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
       IF (sor_norb < 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Number of sample orbits must be positive.", 1)
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
    IF (PRESENT(generat_multiplier) .AND. PRESENT(sor_generat_offset)) THEN
       IF (generat_multiplier <= 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Scaling of generation window is zero or negative.", 1)
          RETURN
       END IF
       this%generat_multiplier_prm = generat_multiplier
       CALL setGenerationWindow(this, sor_generat_offset)
    ELSE IF (PRESENT(generat_multiplier)) THEN
       IF (generat_multiplier <= 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / setParameters", &
               "Scaling of generation window is zero or negative.", 1)
          RETURN
       END IF
       this%generat_multiplier_prm = generat_multiplier
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
    IF (PRESENT(vomcmc_norb)) THEN
       this%vomcmc_norb_prm = vomcmc_norb
    END IF
    IF (PRESENT(vomcmc_ntrial)) THEN
       this%vomcmc_ntrial_prm = vomcmc_ntrial
    END IF
    IF (PRESENT(vomcmc_norb_iter)) THEN
       this%vomcmc_norb_iter_prm = vomcmc_norb_iter
    END IF
    IF (PRESENT(vomcmc_ntrial_iter)) THEN
       this%vomcmc_ntrial_iter_prm = vomcmc_ntrial_iter
    END IF
    IF (PRESENT(vomcmc_nmap)) THEN
       this%vomcmc_nmap_prm = vomcmc_nmap
    END IF
    IF (PRESENT(vomcmc_niter)) THEN
       this%vomcmc_niter_prm = vomcmc_niter
    END IF
    IF (PRESENT(vomcmc_scaling)) THEN
       this%vomcmc_scaling_prm = vomcmc_scaling
    END IF
    IF (PRESENT(vomcmc_mapping_mask)) THEN
       this%vomcmc_mapping_mask_prm = vomcmc_mapping_mask
    END IF
    IF (PRESENT(ls_correction_factor)) THEN
       this%ls_corr_fac_prm = ls_correction_factor
    END IF
    IF (PRESENT(ls_rchi2_acceptable)) THEN
       this%ls_rchi2_acceptable_prm = ls_rchi2_acceptable
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
    IF (PRESENT(smplx_similarity_tol)) THEN
       this%smplx_similarity_tol_prm = smplx_similarity_tol
    END IF
    IF (PRESENT(smplx_niter)) THEN
       this%smplx_niter_prm = smplx_niter
    END IF
    IF (PRESENT(smplx_force)) THEN
       this%smplx_force_prm = smplx_force
    END IF
    IF (PRESENT(os_norb)) THEN
       this%os_norb_prm = os_norb
    END IF
    IF (PRESENT(os_ntrial)) THEN
       this%os_ntrial_prm = os_ntrial
    END IF
    IF (PRESENT(os_sampling_type)) THEN
       this%os_sampling_type_prm = os_sampling_type
    END IF
    IF (PRESENT(generat_gaussian_deviates)) THEN
       this%generat_gaussian_deviates_prm = generat_gaussian_deviates
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

    INTEGER, DIMENSION(:,:), POINTER      :: idx_pair => NULL()
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
  !! Determines the new range bounds from the 3-sigma cutoff values of
  !! the a posteriori range probability density.
  !!
  !! Returns error.
  !!
  SUBROUTINE setRangeBounds_3sigma(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)    :: this
    REAL(bp), DIMENSION(2)                   :: mean, stdev
    REAL(bp), DIMENSION(:,:), ALLOCATABLE    :: histo
    REAL(bp), DIMENSION(:), ALLOCATABLE      :: topo_range
    REAL(bp), DIMENSION(2,2)                 :: rho_
    REAL(bp)                                 :: rhomin1, rhomax1, &
         rhomin2, rhomax2
    INTEGER                                  :: grid_num, err

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
            std_dev=stdev(1), errstr=errstr)
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
       ALLOCATE(topo_range(SIZE(this%sor_rho_arr_cmp,dim=1)), stat=err)
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
            errstr=errstr)
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
  SUBROUTINE simplexOrbits(this, orb_arr, mask, force_earth_impact_at_epoch)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: this
    TYPE (Orbit), DIMENSION(:), INTENT(inout) :: orb_arr
    LOGICAL, DIMENSION(6), INTENT(in), OPTIONAL :: mask
    LOGICAL, INTENT(in), OPTIONAL :: force_earth_impact_at_epoch

    REAL(bp), PARAMETER :: TINY = 1.0e-10_bp

    TYPE (Observatories) :: obsies
    TYPE (CartesianCoordinates) :: ccoord
    CHARACTER(len=FRAME_LEN) :: frame
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: p, p_init
    REAL(bp), DIMENSION(:), ALLOCATABLE :: y, psum, p_best
    REAL(bp), DIMENSION(6) :: coord_obj, coord_earth
    REAL(bp) :: y_best
    INTEGER(ibp) :: ihi, ilo, ndim ! Global variables (within this subroutine).  
    INTEGER :: err, i, info_verb_, nobs
    LOGICAL :: force_earth_impact_at_epoch_

    IF (PRESENT(mask)) THEN
       ndim = COUNT(mask)
    ELSE
       ndim = 6
    END IF
    IF (PRESENT(force_earth_impact_at_epoch)) THEN
       force_earth_impact_at_epoch_ = force_earth_impact_at_epoch
    ELSE
       force_earth_impact_at_epoch_ = .FALSE.
    END IF

    ALLOCATE(p(ndim+1,ndim), p_init(ndim+1,ndim), y(ndim+1), &
         psum(ndim), p_best(ndim))

    IF (force_earth_impact_at_epoch_) THEN
       CALL NEW(obsies)
       ccoord = getObservatoryCCoord(obsies, "500", this%t_inv_prm)
       CALL rotateToEcliptic(ccoord)
       coord_earth = getCoordinates(ccoord)
       CALL NULLIFY(ccoord)
       CALL NULLIFY(obsies)
    END IF

    IF (info_verb >= 3) THEN
       WRITE(stdout,"(2X,A)") "SIMPLEX OPTIMIZATION"
       WRITE(stdout,"(1X)")
       WRITE(stdout,"(2X,A)") "Parameters:"
       WRITE(stdout,"(2X,A,1X,F10.5)") "smplx_tol_prm", this%smplx_tol_prm
       WRITE(stdout,"(2X,A,1X,L1)") "smplx_force_prm", this%smplx_force_prm
       WRITE(stdout,"(2X,A,1X,F10.8)") "smplx_similarity_tol_prm", this%smplx_similarity_tol_prm
       WRITE(stdout,"(2X,A,1X,I0)") "smplx_niter_prm", this%smplx_niter_prm
    END IF

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
    IF (ASSOCIATED(this%rchi2_arr_cmp)) THEN
       DEALLOCATE(this%rchi2_arr_cmp, stat=err)
    END IF

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
    DO i=1,ndim+1
       p(i,1:6) = getElements(orb_arr(i), this%element_type_prm, &
            frame=frame)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / simplexOrbits", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       p_init(i,1:6) = p(i,1:6)
       info_verb_ = info_verb
       info_verb = info_verb_ - 1
       y(i) = getChi2(this, orb_arr(i))
       info_verb = info_verb_
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
       DO i=1,ndim+1
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
       DEALLOCATE(p, p_init, y, psum, p_best)
       RETURN
    END IF
    ALLOCATE(this%orb_arr_cmp(ndim+1), this%rchi2_arr_cmp(ndim+1), stat=err)
    DO i=1,ndim+1
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
            asteroid_perturbers=this%ast_perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / simplexOrbits", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       this%rchi2_arr_cmp(i) = y(i) - REAL(COUNT(this%obs_masks_prm),bp)
    END DO
    this%orb_ml_cmp = copy(this%orb_arr_cmp(1))
    DEALLOCATE(p, p_init, y, psum, p_best)

  CONTAINS 

    SUBROUTINE simplex_private 

      IMPLICIT NONE
      TYPE (Orbit) :: orb
      REAL(bp), DIMENSION(ndim) :: ptry_ref, ptry_exp, ptry_co
      REAL(bp) :: rtol, ytry_ref, ytry_exp, ytry_co, ytmp
      INTEGER(ibp) :: i, j, inhi
      LOGICAL :: contraction_successful

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
         ! reinitialize or exit if very similar values
         IF (2.0_bp*ABS(y(ihi)-y(ilo))/(ABS(y(ihi))+ABS(y(ilo))+TINY) &
              < this%smplx_similarity_tol_prm) THEN
            nobs = getNrOfObservations(this%obss)
            !WRITE(*,*) y_best, nobs,  y_best-2*nobs, this%smplx_tol_prm
            IF (this%smplx_force_prm .AND. &
                 y_best-2*nobs > this%smplx_tol_prm) THEN 
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
                       asteroid_perturbers=this%ast_perturbers_prm, &
                       integrator=this%integrator_prm, &
                       integration_step=this%integration_step_prm)
                  info_verb_ = info_verb
                  info_verb = info_verb_ - 1
                  y(i) = getChi2(this, orb_arr(i))
                  info_verb = info_verb_
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
            ELSE
               ! If returning, put best point and value in slot 1. 
               CALL swap(y(1),y(ilo))
               CALL swap(p(1,:),p(ilo,:)) 
               ilo = 1
               RETURN 
            END IF
         END IF
         IF (info_verb >= 3) THEN
            DO i=1,ndim+1
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
         IF (this%smplx_force_prm .AND. rtol < this%smplx_tol_prm) THEN 
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
                          asteroid_perturbers=this%ast_perturbers_prm, &
                          integrator=this%integrator_prm, &
                          integration_step=this%integration_step_prm)
                     IF (error) THEN
                        CALL errorMessage("StochasticOrbit / simplexOrbits / simplex_private", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     info_verb_ = info_verb
                     info_verb = info_verb_ - 1
                     y(i) = getChi2(this, orb)
                     info_verb = info_verb_
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
      INTEGER :: i

      ptry(:) = (1.0_bp+fac)*(psum(:)-p(ihi,:))/ndim - fac*p(ihi,:)
      ! Evaluate the function at the trial point.
      IF (this%element_type_prm == "keplerian") THEN
         IF (ptry(1) < this%apriori_a_min_prm .OR. &
              ptry(1) > this%apriori_a_max_prm .OR. &
              ptry(2) < 0.0_bp .OR. &
              ptry(3) < 0.0_bp .OR. &
              ptry(3) > pi) THEN
            !CALL randomNumber(ran)
            !ptry = p_init(1+NINT(6*ran),:)
            DO i=1,6
               ptry(i) = SUM(p_init(:,i))
            END DO
            ptry = ptry/7
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
           asteroid_perturbers=this%ast_perturbers_prm, &
           integrator=this%integrator_prm, &
           integration_step=this%integration_step_prm)
      IF (error) THEN
         CALL errorMessage("StochasticOrbit / simplexOrbits / simtry", &
              "TRACE BACK (45)", 1)
         RETURN
      END IF
      ! Get heliocentric cartesion coordinates if forcing Earth impact
      IF (force_earth_impact_at_epoch_) THEN
         coord_obj = getElements(orb, "cartesian", "ecliptic")
         IF (SQRT(SUM((coord_obj(1:3) - coord_earth(1:3))**2)) > planetary_radii(3) + 60.0_bp/km_au) THEN
            ytry = HUGE(ytry)
            RETURN
         END IF
      END IF
      info_verb_ = info_verb
      info_verb = info_verb_ - 1
      ytry = getChi2(this, orb)
      info_verb = info_verb_
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
  !!      - No random observation selection (note: is available but cannot be used because of the multi-propagation scheme)
  !!
  !! Returns error.
  !!
  SUBROUTINE statisticalRanging(this)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)                :: this

    TYPE (Orbit), DIMENSION(:), POINTER                  :: orb_arr_ => NULL()
    TYPE (Orbit), DIMENSION(:), POINTER                  :: orb_arr => NULL()
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER   :: obsies_ccoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: scoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER   :: obs_scoords => NULL()
    TYPE (CartesianCoordinates)                          :: obs_ccoord_topo1, &
         obs_ccoord_focus1, obs_ccoord_topo2, obs_ccoord_focus2
    TYPE (SphericalCoordinates)                          :: obs_scoord1, &
         obs_scoord2, obs_scoord_
    TYPE (Time)                                          :: t1, t2
    CHARACTER(len=DYN_MODEL_LEN)                         :: dyn_model
    CHARACTER(len=INTEGRATOR_LEN)                        :: integrator
    REAL(bp), DIMENSION(:,:,:,:), POINTER                :: partials_arr => NULL()
    REAL(bp), DIMENSION(:,:,:), POINTER                  :: cov_matrices => NULL(), &
         sphdev => NULL()
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE              :: residuals
    REAL(bp), DIMENSION(:,:,:), POINTER                  :: information_matrix_obs => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE                :: rms, &
         rho_distribution, jacobians
    REAL(bp), DIMENSION(:), POINTER                      :: &
         rho1 => NULL(), &
         rho2 => NULL(), &
         mjd_lt => NULL()
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
         apriori, jac_sph_inv, rho_comp1, rho_comp2, tdt, &
         rho_mid, a, ran, tmp, chi2, dchi2, integration_step, &
         trials_per_orb, jac_car_kep, jac_equ_kep, obs_arc, &
         rang, distance, q, ftol
    INTEGER, DIMENSION(:,:), POINTER                     :: &
         obs_pair_arr => NULL()
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

    IF (this%chi2_min_prm < 0.0_bp) THEN
       first = .TRUE.
       IF (this%chi2_min_init_prm <= 0.0_bp) THEN
          this%chi2_min_prm = REAL(COUNT(this%obs_masks_prm),bp)
       ELSE
          this%chi2_min_prm = this%chi2_min_init_prm
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

    IF (this%sor_random_obs_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / StatisticalRanging", &
            "Use of random observation pairs not accepted.", 1)
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

    ! Extract focuscentric equatorial observatory coordinates
    obsies_ccoords => getObservatoryCCoords(this%obss)
!!$! XXX
!!$obsies_ccoords(1)%position = 0.0_bp
!!$obsies_ccoords(1)%velocity = 0.0_bp
!!$obsies_ccoords(2)%position = 0.0_bp
!!$obsies_ccoords(2)%velocity = 0.0_bp
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
       IF (.NOT.this%dchi2_rejection_prm) THEN
          WRITE(stdout,"(2X,A)") "WARNING: NOT using dchi2 in acceptance!"
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
       IF (.NOT.this%dchi2_rejection_prm) THEN
          WRITE(stdout,"(2X,A)") "NOT using dchi2 in acceptance!"
       ENDIF
       IF (.NOT. this%regularization_prm) WRITE(stdout,"(2X,A)",advance="no") "NOT"
       WRITE(stdout,"(2X,A,L2)") "Using regularization by Jeffreys"
       WRITE(stdout,"(2X,A,A)") "Epoch                   : ", &
            TRIM(getCalendarDateString(this%t_inv_prm,"tdt"))
       WRITE(stdout,'(2X,A,E10.4)') "Initial minimum chi2     : ", &
            this%chi2_min_prm
       WRITE(stdout,"(2X,A,A,A)") "Using ",TRIM(dyn_model)," dynamical model"
       WRITE(stdout,"(2X,A,A,A)") "Using ", &
            TRIM(this%sor_2point_method_prm), " method"
       IF (this%sor_gaussian_pdf_prm) THEN
          WRITE(stdout,"(2X,A)") "Using Gaussian deviate for rho!"
       END IF
       WRITE(stdout,"(2X,A,4(1X,E10.4))") "Initial rho bounds (au) :", &
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
          CALL NULLIFY(obs_ccoord_focus1)
          CALL NULLIFY(obs_ccoord_topo2)
          CALL NULLIFY(obs_ccoord_focus2)
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
                     rho1(i+1), " au)"
             END IF
             CYCLE sor_orb_gen
          END IF
          IF (this%informative_apriori_prm) THEN
             IF (this%apriori_rho_min_prm >= 0.0_bp .AND. &
                  rho1(i+1) < this%apriori_rho_min_prm) THEN
                ! rho1 too small
                IF (info_verb >= 5) THEN
                   WRITE(stdout,"(2X,A,F13.7,A)") &
                        "Failed (rho1 too small: ", a, " au)"
                END IF
                CYCLE sor_orb_gen
             END IF
             IF (this%apriori_rho_max_prm >= 0.0_bp .AND. &
                  rho1(i+1) > this%apriori_rho_max_prm) THEN
                ! rho1 too large
                IF (info_verb >= 5) THEN
                   WRITE(stdout,"(2X,A,F10.7,A)") &
                        "Failed (rho1 too large: ", a, " au)"
                END IF
                CYCLE sor_orb_gen
             END IF
          END IF
          rho_mid = 0.5_bp*(this%sor_rho_prm(2,1) + this%sor_rho_prm(2,2))
          bounds(1,1) = rho1(i+1) + rho_mid 
          IF (bounds(1,1) <= planetary_radii(3)) THEN
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,F10.7,A)") &
                     "Failed (Mean of rho2 smaller than the Earth radius: ", &
                     bounds(1,1), " au)"
             END IF
             CYCLE sor_orb_gen
          END IF
          bounds(2:3,1:2) = this%sor_deviates_prm(obs_pair_arr(i+1,2),2:3,1:2)
          rho2(i+1) = 0.0_bp
          ! Make sure rho2 is larger than zero length units
          DO WHILE (ABS(rho2(i+1)) <= EPSILON(rho2(i+1)))
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
             IF (this%informative_apriori_prm) THEN
                IF (this%apriori_rho_min_prm >= 0.0_bp .AND. &
                     rho2(i+1) < this%apriori_rho_min_prm) THEN
                   ! rho2 too small
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F13.7,A)") &
                           "Failed (rho2 too small: ", a, " au)"
                   END IF
                   rho2(i+1) = 0.0_bp ! try again...
                END IF
                IF (this%apriori_rho_max_prm >= 0.0_bp .AND. &
                     rho2(i+1) > this%apriori_rho_max_prm) THEN
                   ! rho2 too large
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (rho2 too large: ", a, " au)"
                   END IF
                   rho2(i+1) = 0.0_bp ! try again...
                END IF
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

          ! Cartesian focuscentric coordinates:
          obs_ccoord_focus1 = &
               copy(obsies_ccoords(obs_pair_arr(i+1,1)) + obs_ccoord_topo1)
          obs_ccoord_focus2 = &
               copy(obsies_ccoords(obs_pair_arr(i+1,2)) + obs_ccoord_topo2)

          IF (info_verb >= 5) THEN
             t1 = getTime(obsies_ccoords(obs_pair_arr(i+1,1)))
             t2 = getTime(obs_ccoord_topo1)
             WRITE(stdout,"(2X,A,2(1X,F15.7))") &
                  "Focuscentric times 1:", &
                  getMJD(t1,"tdt"), getMJD(t2,"tdt")
             position1 = getPosition(obs_ccoord_focus1)
             position2 = getPosition(obsies_ccoords(obs_pair_arr(i+1,1)))
             WRITE(stdout,"(2X,A,1X,3(F15.10,1X))") &
                  "Focuscentric coordinates for observatory at t0", position2
             WRITE(stdout,"(2X,A,1X,3(F15.10,1X))") &
                  "Focuscentric coordinates generated for object at t0", position1
             position1 = getPosition(obs_ccoord_focus2)
             position2 = getPosition(obsies_ccoords(obs_pair_arr(i+1,2)))
             WRITE(stdout,"(2X,A,1X,3(F15.10,1X))") &
                  "Focuscentric coordinates for observatory at t1", position2
             WRITE(stdout,"(2X,A,1X,3(F15.10,1X))") &
                  "Focuscentric coordinates generated for object at t1", position1
             CALL NULLIFY(t1)
             CALL NULLIFY(t2)
          END IF

          ! ------------------------------------------------
          ! Cartesian ecliptic position and velocity from 
          ! two cartesian ecliptic coordinates. P-iteration.
          ! ------------------------------------------------

          CALL estimateLightTime(obs_ccoord_focus1, rho1(i+1))
          CALL estimateLightTime(obs_ccoord_focus2, rho2(i+1))
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
             CALL NULLIFY(obs_ccoord_focus1)
             CALL NULLIFY(obs_ccoord_focus2)
             RETURN
          END IF
          t1 = getTime(obs_ccoord_focus1)
          mjd_lt(i+1) = getMJD(t1,"TDT")
          CALL NULLIFY(t1)
          t1 = getTime(obs_scoords(obs_pair_arr(i+1,1)))

          IF (info_verb >= 5) THEN
             !             t2 = getTime(obs_ccoord_focus2)
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
          !err_verb = 0
          ! Find orbit candidate at the epoch of the first observation by
          ! using the chosen method to solve the 2-point boundary value
          ! problem:
          CALL NEW(orb_arr(i+1), obs_ccoord_focus1, obs_ccoord_focus2, &
               this%sor_2point_method_prm, this%apriori_a_max_prm, &
               ftol=ftol, perturbers=this%perturbers_prm, &
               asteroid_perturbers=this%ast_perturbers_prm, &
               integrator=this%integrator_prm, &
               integration_step=this%integration_step_prm, &
               center=this%center_prm)
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
                           "Failed (semimajor axis too small: ", a, " au)"
                   END IF
                   failed_flag(2) = failed_flag(2) + 1
                   CYCLE sor_orb_gen
                END IF
                IF (this%apriori_a_max_prm >= 0.0_bp .AND. &
                     a > this%apriori_a_max_prm) THEN
                   ! Semimajor axis too large
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (semimajor axis too large: ", a, " au)"
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
                   error = .FALSE.
                   CYCLE sor_orb_gen
                END IF
                ! Periapsis distance too small:
                IF (this%apriori_periapsis_min_prm >= 0.0_bp .AND. &
                     q < this%apriori_periapsis_min_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F13.7,A)") &
                           "Failed (periapsis distance too small: ", a, " au)"
                   END IF
                   failed_flag(4) = failed_flag(4) + 1
                   CYCLE sor_orb_gen
                END IF
                ! Periapsis distance too large:
                IF (this%apriori_periapsis_max_prm >= 0.0_bp .AND. &
                     q > this%apriori_periapsis_max_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (periapsis distance too large: ", a, " au)"
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
                           "Failed (apoapsis distance too small: ", a, " au)"
                   END IF
                   failed_flag(4) = failed_flag(4) + 1
                   CYCLE sor_orb_gen
                END IF
                ! Apoapsis distance too large:
                IF (this%apriori_apoapsis_max_prm >= 0.0_bp .AND. &
                     Q > this%apriori_apoapsis_max_prm) THEN
                   IF (info_verb >= 5) THEN
                      WRITE(stdout,"(2X,A,F10.7,A)") &
                           "Failed (apoapsis distance too large: ", a, " au)"
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
               asteroid_perturbers=this%ast_perturbers_prm, &
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
             CALL NULLIFY(obs_ccoord_focus1)
             CALL NULLIFY(obs_ccoord_focus2)
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
                CALL NULLIFY(obs_ccoord_focus1)
                CALL NULLIFY(obs_ccoord_focus2)
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
       CALL NULLIFY(obs_ccoord_focus1)
       CALL NULLIFY(obs_ccoord_focus2)

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
             !write(*,*) "RES", observed_coord(2:3)/rad_deg, computed_coord(2:3)/rad_deg, residuals(iorb+1,j,2:3)/rad_asec
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
                !write(*,*) ABS(residuals(iorb+1,:,2:3))
                !write(*,*) this%res_accept_prm(:,2:3)
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

          ! Compute dchi2 wrt best fit orbit
          dchi2 = chi2 - this%chi2_min_prm
          IF (this%dchi2_rejection_prm .AND. &
               dchi2 > this%dchi2_prm) THEN
             ! The dchi2 is used and its value is not acceptable.
             failed_flag(7) = failed_flag(7) + 1
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,1X,E10.5)") &
                     "Failed (dchi2 too large)", dchi2
             END IF
             CYCLE sor_orb_acc
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
                   WRITE(stderr,"(A)") "cosdec0_arr(obs_pair_arr(i,1)):"
                   WRITE(stderr,"(F15.10)") cosdec0_arr(obs_pair_arr(i,1))
                   WRITE(stderr,"(A)") "partials_arr(i,1:3,:,obs_pair_arr(i,1)) :"
                   CALL matrix_print(partials_arr(i,1:3,:,obs_pair_arr(i,1)) , stderr, errstr)
                   WRITE(stderr,"(A)") "cosdec0(nobs):"
                   WRITE(stderr,"(F15.10)") cosdec0_arr(obs_pair_arr(i,2))
                   WRITE(stderr,"(A)") "partials_arr(i,1:3,:,obs_pair_arr(i,2)):"
                   CALL matrix_print(partials_arr(i,1:3,:,obs_pair_arr(i,2)), stderr, errstr)
                   WRITE(stderr,"(A)") "jacobian_matrix:"
                   CALL matrix_print(jacobian_matrix, stderr, errstr)
                END IF
                errstr = ""
                CYCLE
             END IF

             ! Compute determinant of Jacobian between Cartesian and
             ! Keplerian orbital elements ("Cartesian Wrt
             ! Keplerian"). Temporarily reduce error messaging if
             ! hyperbolic orbits are accepted.
             IF (this%apriori_a_min_prm < 0.0_bp) THEN
                err_verb = err_verb - 1
             END IF
             CALL partialsCartesianWrtKeplerian(orb_arr(i), &
                  jacobian_matrix, "equatorial")
             IF (this%apriori_a_min_prm < 0.0_bp) THEN
                err_verb = err_verb + 1
             END IF
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
             elements = getElements(orb_arr(i), "cometary")
             jac_equ_kep = 0.5_bp*elements(2) * &
                  SIN(0.5_bp*elements(3)) / COS(0.5_bp*elements(3))**3
          ELSE
             jac_sph_inv = 1.0_bp
             jac_car_kep = 1.0_bp
             jac_equ_kep = 1.0_bp
          END IF

          ! Probability density function (note that the '- nobs' term
          ! is there for practical reasons):
          pdf_val = apriori*EXP(-0.5_bp*(chi2 - COUNT(this%obs_masks_prm)))/jac_sph_inv

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
          ! Nonlinear "reduced chi2" 
          ! (note that usually rchi2 = chi2/(number of measurements - number of fitted parameters))
          rchi2_arr(iorb) = chi2 - REAL(COUNT(this%obs_masks_prm),bp)

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
       norb = NINT((this%sor_norb_prm - iorb)*(1.0_bp*itrial/MAX(REAL(iorb,bp),1.0_bp)))
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,1X,I0,A,I0)") "Number of sample orbits accepted so far:", iorb, "/", this%sor_norb_prm
          WRITE(stdout,"(2X,A,1X,I0)") "Number of trial orbits generated so far:", itrial
          WRITE(stdout,"(2X,A,1X,I0)") "Number of trial orbits expected to be required for completion:", norb
       END IF
       IF (norb > norb_simult_max) THEN
          norb = norb_simult_max
       END IF
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,1X,I0)") "Number of trial orbits to be generated for next batch:", norb
          WRITE(stdout,*)
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
       i = MAXLOC(this%pdf_arr_cmp,dim=1)
       CALL NULLIFY(this%orb_ml_cmp)
       this%orb_ml_cmp        = copy(this%orb_arr_cmp(i)) 
       this%pdf_arr_cmp       = pdf_arr(1:iorb)
       this%reg_apr_arr_cmp   = reg_apriori_arr(1:iorb)
       this%jac_arr_cmp       = jacobians(1:iorb,1:3)
       this%rchi2_arr_cmp     = rchi2_arr(1:iorb)
       this%rms_arr_cmp       = rms(1:iorb,:)
       this%sor_rho_arr_cmp   = rho_distribution(1:iorb,:)
       this%res_arr_cmp       = residuals(1:iorb,:,:)
       this%chi2_min_cmp      = MINVAL(rchi2_arr(1:iorb)) + REAL(COUNT(this%obs_masks_prm),bp)
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
       WRITE(stdout,"(2X,A,2(1X,E10.4))") "Initial - final minimum chi2 values    :", &
            this%chi2_min_prm - this%chi2_min_cmp
       WRITE(stdout,"(2X,A,4(1X,E10.4))") "Min/max rhos, @ end    :" , &
            this%sor_rho_cmp(1,1:2), this%sor_rho_cmp(2,1:2)
       acc(2:3) = acc(2:3)/rad_asec
       WRITE(stdout,"(2X,A,3(1X,E10.4),A)") "Accur. (rho/R.A./Dec.) :", acc, " (au/as/as)"
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

    TYPE (Observations)                       :: obss, obss_next
    TYPE (Observation), DIMENSION(:), POINTER :: obs_arr => NULL()
    TYPE (Observation)                        :: obs_next
    TYPE (StochasticOrbit)                    :: storb
    CHARACTER(len=4)                          :: str1, str2
    REAL(bp), DIMENSION(4)                    :: rho
    REAL(bp)                                  :: chi2_min_, chi2_min_final, ddchi2
    INTEGER(ibp)                              :: nobs, nobs_max_, err, &
         i, j, k

    IF (.NOT. this%is_initialized_prm) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / stepwiseRanging", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    IF (ASSOCIATED(this%orb_arr_cmp)) THEN
       CALL constrainRangeDistributions(this, this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / stepwiseRanging", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       CALL setRangeBounds(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / stepwiseRanging", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
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
          IF (i /= nobs) THEN 
             CALL NEW(obss_next, obs_arr(k))
          END IF
          j = j + 1
       ELSE
          CALL addObservation(obss, obs_arr(k))
          IF (i /= nobs) THEN 
             CALL NEW(obss_next, obs_arr(j))
          END IF
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
            asteroid_perturbers=this%ast_perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm, &
            t_inv=this%t_inv_prm, &
            element_type=this%element_type_prm, &
            multiple_objects=this%multiple_obj_prm, &
            outlier_rejection=this%outlier_rejection_prm, &
            dchi2_rejection=this%dchi2_rejection_prm, &
            center=this%center_prm, &
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
            apriori_rho_max=this%apriori_rho_max_prm, &
            sor_norb=this%sor_norb_sw_prm, &
                                ! Number of trial orbits increases with number of observations:
                                !sor_ntrial=this%sor_ntrial_sw_prm*i, &
            sor_ntrial=this%sor_ntrial_sw_prm, &
            sor_rho1_l=rho(1), &
            sor_rho1_u=rho(2), &
            sor_rho2_l=rho(3), &
            sor_rho2_u=rho(4), &
            sor_niter=1, &
            generat_multiplier=this%generat_multiplier_prm, &
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
          storb%chi2_min_prm = chi2_min_
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
             CALL toString(i, str2, error)
             CALL errorMessage("StochasticOrbit / stepwiseRanging", &
                  TRIM(str1) // "sample orbits found when " // &
                  TRIM(str2) // " observations were included." , 1)
             error = .TRUE.
             DEALLOCATE(obs_arr, stat=err)
             CALL NULLIFY(obss)
             CALL NULLIFY(storb)
             RETURN
          END IF
          IF (exist(obss_next)) THEN
             CALL constrainRangeDistributions(storb, obss_next)
             CALL NULLIFY(obss_next)
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
          chi2_min_ = storb%chi2_min_cmp
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
            (this%sor_rho_histo_cmp > 0 .OR. &
            (this%dchi2_rejection_prm .AND. ddchi2 > 2.0_bp) .OR. &
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
          this%chi2_min_prm = chi2_min_
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
          ddchi2  = MAXVAL(this%rchi2_arr_cmp) - MINVAL(this%rchi2_arr_cmp) - this%dchi2_prm
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
          chi2_min_ = MIN(this%chi2_min_cmp,this%chi2_min_prm)
       END DO
    ELSE
       CALL NULLIFY(this)
       this = copy(storb)
    END IF
    CALL NULLIFY(storb)

  END SUBROUTINE stepwiseRanging





  !! *Description*:
  !!
  !! Optimizes the range distribution corresponding to the first and
  !! last observation of a set of observations which is a combination
  !! of this%obss and the additional observations supplied with this
  !! subroutine. For the estimation of the ranges the algorithm
  !! selects all orbits which reproduce the observations with
  !! acceptable residuals (limit is set by the 1-sigma astrometric
  !! uncertainty multiplied by this%accept_multiplier_prm) or, if i<10
  !! orbits with acceptable residuals are found, selects the 10-i
  !! orbits that produce the smallest rchi2 values.
  !!
  !! Also updates orb_ml
  !! Sets error=.TRUE. if an error occurs.
  !!
  SUBROUTINE constrainRangeDistributions2(this, obss)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this
    TYPE (Observations), INTENT(in) :: obss

    TYPE (Orbit), DIMENSION(:), POINTER :: orb_arr => NULL()
    TYPE (Observations) :: obss_
    TYPE (Observation) :: obs
    TYPE (CartesianCoordinates), DIMENSION(2) :: observers
    TYPE (SphericalCoordinates) , DIMENSION(:,:), POINTER :: &
         ephemerides => NULL()
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         residuals => NULL(), &
         information_matrix_obs => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: &
         stdevs => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         sor_rho_arr_cmp
    REAL(bp), DIMENSION(:), POINTER :: &
         dates_orig => NULL(), &
         dates_add => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: chi2_arr
    INTEGER :: err, i, j, norb, imin
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask, mask_

    ! Get residuals between predicted positions and additional
    ! observations:
    residuals => getResiduals(this, obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(residuals, stat=err)
       RETURN
    END IF

    ! Get astrometric uncertainty for additional observations:
    stdevs => getStandardDeviations(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(stdevs, stat=err)
       RETURN
    END IF

    ! Find out which orbits reproduce the additional observations
    ! within the set limits:
    norb = SIZE(residuals,dim=2)
    ALLOCATE(mask(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not allocate memory (5)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(stdevs, stat=err)
       DEALLOCATE(mask, stat=err)
       RETURN
    END IF

    mask = .TRUE.
    DO i=1,norb
       ! Note that RA,Dec is hardwired here
       IF (ANY(ABS(residuals(:,i,2:3)) > this%accept_multiplier_prm*stdevs(:,2:3))) THEN
          mask(i) = .FALSE.
       END IF
    END DO
    DEALLOCATE(stdevs, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not deallocate memory (5)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       RETURN
    END IF
    IF (info_verb >= 2 .AND. COUNT(mask) > 0) THEN
       WRITE(stdout,"(2X,2(A,1X),I0,1X,A)") "constrainRangeDistributions2:", &
            "RA,Dec residuals corresponding to the", COUNT(mask), &
            "orbits having residuals smaller than the acceptance window [asec]:"
       DO i=1,SIZE(mask)
          IF (info_verb >= 3 .AND. mask(i)) THEN
             DO j=1,SIZE(residuals,dim=1)
                WRITE(stdout,"(2X,2(A,1X,I0,1X),A,2(1X,F10.3))") &
                     "Orbit #", i, "& observation #", j, ":", residuals(j,i,2:3)/rad_asec
             END DO
          END IF
       END DO
    END IF

    ! Always compute chi2 
    information_matrix_obs => getBlockDiagInformationMatrix(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "TRACE BACK (15)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       RETURN
    END IF
    ALLOCATE(chi2_arr(norb), mask_(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not allocate memory (10)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(information_matrix_obs, stat=err)
       DEALLOCATE(chi2_arr, stat=err)
       DEALLOCATE(mask_, stat=err)
       RETURN
    END IF
    mask_ = mask
    DO i=1,norb
       chi2_arr(i) = chi_square(residuals(:,i,:), information_matrix_obs, errstr=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               TRIM(errstr), 1)             
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(chi2_arr, stat=err)
          DEALLOCATE(mask_, stat=err)
          RETURN
       END IF
    END DO
    ! Find best-fit orbit
    imin = MINLOC(chi2_arr, dim=1)

    ! If less than 10 orbits acceptably reproduce the additional
    ! observations, then select (10 - norb) orbits that have the
    ! smallest chi2 wrt the additional observations:
    IF (COUNT(mask) < 10) THEN
       j = COUNT(mask)
       DO WHILE (COUNT(mask) < 10 .AND. j < norb)
          j = j + 1
          i = MINLOC(chi2_arr,1,.NOT.mask)
          mask(i) = .TRUE.
       END DO
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A,1X,I0,1X,A)") "constrainRangeDistributions2:", &
               "RA,Dec residuals corresponding to the", 10-COUNT(mask_), &
               "additional orbits included [asec]:"
          DO i=1,SIZE(mask)
             IF (info_verb >= 3 .AND. mask(i) .AND. .NOT.mask_(i)) THEN
                DO j=1,SIZE(residuals,dim=1)
                   WRITE(stdout,"(2X,2(A,1X,I0,1X),A,2(1X,F10.3))") &
                        "Orbit #", i, "& observation #", j, ":", residuals(j,i,2:3)/rad_asec
                END DO
             END IF
          END DO
       END IF
       DEALLOCATE(information_matrix_obs, chi2_arr, mask_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "Could not deallocate memory (10)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(chi2_arr, stat=err)
          DEALLOCATE(mask_, stat=err)
          RETURN
       END IF
    END IF

    ! Get observation dates for original data and additional data
    IF (exist(this%obss)) THEN
       dates_orig => getDates(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "TRACE BACK (20)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          RETURN
       END IF
       dates_add => getDates(obss)       
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "TRACE BACK (25)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          RETURN
       END IF
    ELSE
       ALLOCATE(dates_orig(1), dates_add(1), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "Could not allocate memory (15)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          RETURN
       END IF
       dates_orig = 1.0_bp
       dates_add  = 0.0_bp
    END IF

    ! Re-calculate this%sor_rho_arr_cmp if...
    IF (& 
                                ! ...rhos don't exist:
         .NOT.ASSOCIATED(this%sor_rho_arr_cmp) .OR. & 
                                ! ...previous observations don't exist -> no way to find out if
                                ! update needed:
         .NOT.exist(this%obss) .OR. & 
                                ! ...additional data earlier than original data:
         MINVAL(dates_add) < MINVAL(dates_orig) .OR. & 
                                ! ...additional data later than original data
         MAXVAL(dates_add) > MAXVAL(dates_orig)) THEN 
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A)") "constrainRangeDistributions2:", &
               "Re-calculating the rho1 and rho2 distributions..."
       END IF
       IF (.NOT.exist(this%obss)) THEN
          obs = getObservation(obss,1)
          observers(1) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
          obs = getObservation(obss,getNrOfObservations(obss))
          observers(2) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
       ELSE
          obss_ = this%obss + obss
          obs = getObservation(obss_,1)
          observers(1) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
          obs = getObservation(obss_,getNrOfObservations(obss_))
          observers(2) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
          CALL NULLIFY(obss_)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "TRACE BACK (30)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          RETURN
       END IF
       orb_arr => getSampleOrbits(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "TRACE BACK (35)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err) 
          DEALLOCATE(orb_arr, stat=err)         
          RETURN
       END IF
       ! Save best-fit orbit (only if not already defined)
       IF (.NOT.exist(this%orb_ml_cmp)) THEN
          this%orb_ml_cmp = copy(orb_arr(imin))
       END IF

       CALL getEphemerides(orb_arr, observers, ephemerides)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "TRACE BACK (40)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          DEALLOCATE(orb_arr, stat=err)         
          DEALLOCATE(ephemerides, stat=err)         
          RETURN
       END IF
       IF (ASSOCIATED(this%sor_rho_arr_cmp)) THEN
          DEALLOCATE(this%sor_rho_arr_cmp, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
                  "Could not deallocate memory (15)", 1)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(mask, stat=err)
             DEALLOCATE(dates_orig, stat=err)
             DEALLOCATE(dates_add, stat=err)
             DEALLOCATE(orb_arr, stat=err)         
             DEALLOCATE(ephemerides, stat=err)         
             RETURN
          END IF
       END IF
       ALLOCATE(this%sor_rho_arr_cmp(norb,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "Could not allocate memory (20)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          DEALLOCATE(orb_arr, stat=err)         
          DEALLOCATE(ephemerides, stat=err)         
          RETURN
       END IF
       DO i=1,norb
          this%sor_rho_arr_cmp(i,1) = getDistance(ephemerides(i,1))
          this%sor_rho_arr_cmp(i,2) = getDistance(ephemerides(i,2))
          CALL NULLIFY(ephemerides(i,1))
          CALL NULLIFY(ephemerides(i,2))
          CALL NULLIFY(orb_arr(i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
                  "TRACE BACK (45)", 1)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(mask, stat=err)
             DEALLOCATE(dates_orig, stat=err)
             DEALLOCATE(dates_add, stat=err)
             DEALLOCATE(orb_arr, stat=err)         
             DEALLOCATE(ephemerides, stat=err)         
             RETURN
          END IF
       END DO
       CALL NULLIFY(observers(1))
       CALL NULLIFY(observers(2))
       DEALLOCATE(ephemerides, orb_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
               "Could not deallocate memory (20)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          DEALLOCATE(orb_arr, stat=err)         
          DEALLOCATE(ephemerides, stat=err)         
          RETURN
       END IF
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A)") "constrainRangeDistributions2:", &
               "Re-calculating the rho1 and rho2 distributions... done"
       END IF
    END IF
    DEALLOCATE(dates_orig, dates_add, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not deallocate memory (25)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(dates_orig, stat=err)
       DEALLOCATE(dates_add, stat=err)
       RETURN
    END IF

    ! Rewrite sor_rho_arr_cmp (N.B. The size of the 1st dimension of
    ! sor_rho_arr_cmp is no longer necessary equal to norb!)
    ALLOCATE(sor_rho_arr_cmp(COUNT(mask),2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not allocate memory (25)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(sor_rho_arr_cmp, stat=err)
       RETURN
    END IF
    IF (info_verb >= 3) THEN
       WRITE(stdout,"(2X,A,1X,A)") "constrainRangeDistributions2:", &
            "Constrained (rho1,rho2-rho1) distribution [au]: "
    END IF
    j = 0
    DO i=1,SIZE(this%sor_rho_arr_cmp,dim=1)
       IF (mask(i)) THEN
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(2X,2(1X,F11.7))") &
                  this%sor_rho_arr_cmp(i,1), &
                  this%sor_rho_arr_cmp(i,2)-this%sor_rho_arr_cmp(i,1)
          END IF
          j = j + 1
          sor_rho_arr_cmp(j,:) = this%sor_rho_arr_cmp(i,:)
       END IF
    END DO
    DEALLOCATE(this%sor_rho_arr_cmp, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not deallocate memory (30)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(sor_rho_arr_cmp, stat=err)
       RETURN
    END IF
    ALLOCATE(this%sor_rho_arr_cmp(COUNT(mask),2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not allocate memory (30)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       RETURN
    END IF
    this%sor_rho_arr_cmp = sor_rho_arr_cmp
    DEALLOCATE(sor_rho_arr_cmp, residuals, mask, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions2", &
            "Could not deallocate memory (35)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(sor_rho_arr_cmp, stat=err)
       RETURN
    END IF

  END SUBROUTINE constrainRangeDistributions2





  !! *Description*:
  !!
  !! Optimizes the range distribution corresponding to the first and
  !! last observation of a set of observations which is a combination
  !! of this%obss and the additional observations supplied with this
  !! subroutine. For the estimation of the ranges the algorithm
  !! selects all orbits which reproduce the observations with
  !! acceptable residuals (limit is set by the 1-sigma astrometric
  !! uncertainty multiplied by this%accept_multiplier_prm) or, if i<10
  !! orbits with acceptable residuals are found, selects the 10-i
  !! orbits that produce the smallest rchi2 values.
  !!
  !! Sets error=.TRUE. if an error occurs.
  !!
  SUBROUTINE constrainRangeDistributions(this, obss)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this
    TYPE (Observations), INTENT(in) :: obss

    TYPE (Orbit), DIMENSION(:), POINTER :: orb_arr => NULL()
    TYPE (Observations) :: obss_
    TYPE (Observation) :: obs
    TYPE (CartesianCoordinates), DIMENSION(2) :: observers
    TYPE (SphericalCoordinates) , DIMENSION(:,:), POINTER :: &
         ephemerides => NULL()
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         residuals => NULL(), &
         information_matrix_obs => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: &
         stdevs => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: sor_rho_arr_cmp
    REAL(bp), DIMENSION(:), POINTER :: &
         dates_orig => NULL(), &
         dates_add => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: chi2_arr
    INTEGER :: err, i, j, norb
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask, mask_

    ! Get residuals between predicted positions and additional
    ! observations:
    residuals => getResiduals(this, obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(residuals, stat=err)
       RETURN
    END IF

    ! Get astrometric uncertainty for additional observations:
    stdevs => getStandardDeviations(obss)
    IF (error) THEN
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "TRACE BACK (10)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(stdevs, stat=err)
       RETURN
    END IF

    ! Find out which orbits reproduce the additional observations
    ! within the set limits:
    norb = SIZE(residuals,dim=2)
    ALLOCATE(mask(norb), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not allocate memory (5)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(stdevs, stat=err)
       DEALLOCATE(mask, stat=err)
       RETURN
    END IF

    mask = .TRUE.
    DO i=1,norb
       ! Note that RA,Dec is hardwired here
       IF (ANY(ABS(residuals(:,i,2:3)) > this%accept_multiplier_prm*stdevs(:,2:3))) THEN
          mask(i) = .FALSE.
       END IF
    END DO
    DEALLOCATE(stdevs, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not deallocate memory (5)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       RETURN
    END IF
    IF (info_verb >= 2 .AND. COUNT(mask) > 0) THEN
       WRITE(stdout,"(2X,2(A,1X),I0,1X,A)") "constrainRangeDistributions:", &
            "RA,Dec residuals corresponding to the", COUNT(mask), &
            "orbits having residuals smaller than the acceptance window [asec]:"
       DO i=1,SIZE(mask)
          IF (info_verb >= 3 .AND. mask(i)) THEN
             DO j=1,SIZE(residuals,dim=1)
                WRITE(stdout,"(2X,2(A,1X,I0,1X),A,2(1X,F10.3))") &
                     "Orbit #", i, "& observation #", j, ":", residuals(j,i,2:3)/rad_asec
             END DO
          END IF
       END DO
    END IF

    ! If less than 10 orbits acceptably reproduce the additional
    ! observations, then select (10 - norb) orbits that have the
    ! smallest chi2 wrt the additional observations:
    IF (COUNT(mask) < 10) THEN
       information_matrix_obs => getBlockDiagInformationMatrix(obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "TRACE BACK (15)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          RETURN
       END IF
       ALLOCATE(chi2_arr(norb), mask_(norb), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "Could not allocate memory (10)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(chi2_arr, stat=err)
          DEALLOCATE(mask_, stat=err)
          RETURN
       END IF
       mask_ = mask
       DO i=1,norb
          chi2_arr(i) = chi_square(residuals(:,i,:), information_matrix_obs, errstr=errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
                  TRIM(errstr), 1)             
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(mask, stat=err)
             DEALLOCATE(information_matrix_obs, stat=err)
             DEALLOCATE(chi2_arr, stat=err)
             DEALLOCATE(mask_, stat=err)
             RETURN
          END IF
       END DO
       j = COUNT(mask)
       DO WHILE (COUNT(mask) < 10 .AND. j < norb)
          j = j + 1
          i = MINLOC(chi2_arr,1,.NOT.mask)
          mask(i) = .TRUE.
       END DO
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A,1X,I0,1X,A)") "constrainRangeDistributions:", &
               "RA,Dec residuals corresponding to the", 10-COUNT(mask_), &
               "additional orbits included [asec]:"
          DO i=1,SIZE(mask)
             IF (info_verb >= 3 .AND. mask(i) .AND. .NOT.mask_(i)) THEN
                DO j=1,SIZE(residuals,dim=1)
                   WRITE(stdout,"(2X,2(A,1X,I0,1X),A,2(1X,F10.3))") &
                        "Orbit #", i, "& observation #", j, ":", residuals(j,i,2:3)/rad_asec
                END DO
             END IF
          END DO
       END IF
       DEALLOCATE(information_matrix_obs, chi2_arr, mask_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "Could not deallocate memory (10)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(information_matrix_obs, stat=err)
          DEALLOCATE(chi2_arr, stat=err)
          DEALLOCATE(mask_, stat=err)
          RETURN
       END IF
    END IF

    ! Get observation dates for original data and additional data
    IF (exist(this%obss)) THEN
       dates_orig => getDates(this%obss)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "TRACE BACK (20)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          RETURN
       END IF
       dates_add => getDates(obss)       
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "TRACE BACK (25)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          RETURN
       END IF
    ELSE
       ALLOCATE(dates_orig(1), dates_add(1), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "Could not allocate memory (15)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          RETURN
       END IF
       dates_orig = 1.0_bp
       dates_add  = 0.0_bp
    END IF

    ! Re-calculate this%sor_rho_arr_cmp if...
    IF (& 
                                ! ...rhos don't exist:
         .NOT.ASSOCIATED(this%sor_rho_arr_cmp) .OR. & 
                                ! ...previous observations don't exist -> no way to find out if
                                ! update needed:
         .NOT.exist(this%obss) .OR. & 
                                ! ...additional data earlier than original data:
         MINVAL(dates_add) < MINVAL(dates_orig) .OR. & 
                                ! ...additional data later than original data
         MAXVAL(dates_add) > MAXVAL(dates_orig)) THEN 
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A)") "constrainRangeDistributions:", &
               "Re-calculating the rho1 and rho2 distributions..."
       END IF
       IF (.NOT.exist(this%obss)) THEN
          obs = getObservation(obss,1)
          observers(1) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
          obs = getObservation(obss,getNrOfObservations(obss))
          observers(2) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
       ELSE
          obss_ = this%obss + obss
          obs = getObservation(obss_,1)
          observers(1) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
          obs = getObservation(obss_,getNrOfObservations(obss_))
          observers(2) = getObservatoryCCoord(obs)
          CALL NULLIFY(obs)
          CALL NULLIFY(obss_)
       END IF
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "TRACE BACK (30)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          RETURN
       END IF
       orb_arr => getSampleOrbits(this)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "TRACE BACK (35)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err) 
          DEALLOCATE(orb_arr, stat=err)         
          RETURN
       END IF
       CALL getEphemerides(orb_arr, observers, ephemerides)
       IF (error) THEN
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "TRACE BACK (40)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          DEALLOCATE(orb_arr, stat=err)         
          DEALLOCATE(ephemerides, stat=err)         
          RETURN
       END IF
       IF (ASSOCIATED(this%sor_rho_arr_cmp)) THEN
          DEALLOCATE(this%sor_rho_arr_cmp, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
                  "Could not deallocate memory (15)", 1)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(mask, stat=err)
             DEALLOCATE(dates_orig, stat=err)
             DEALLOCATE(dates_add, stat=err)
             DEALLOCATE(orb_arr, stat=err)         
             DEALLOCATE(ephemerides, stat=err)         
             RETURN
          END IF
       END IF
       ALLOCATE(this%sor_rho_arr_cmp(norb,2), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "Could not allocate memory (20)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          DEALLOCATE(orb_arr, stat=err)         
          DEALLOCATE(ephemerides, stat=err)         
          RETURN
       END IF
       DO i=1,norb
          this%sor_rho_arr_cmp(i,1) = getDistance(ephemerides(i,1))
          this%sor_rho_arr_cmp(i,2) = getDistance(ephemerides(i,2))
          CALL NULLIFY(ephemerides(i,1))
          CALL NULLIFY(ephemerides(i,2))
          CALL NULLIFY(orb_arr(i))
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
                  "TRACE BACK (45)", 1)
             DEALLOCATE(residuals, stat=err)
             DEALLOCATE(mask, stat=err)
             DEALLOCATE(dates_orig, stat=err)
             DEALLOCATE(dates_add, stat=err)
             DEALLOCATE(orb_arr, stat=err)         
             DEALLOCATE(ephemerides, stat=err)         
             RETURN
          END IF
       END DO
       CALL NULLIFY(observers(1))
       CALL NULLIFY(observers(2))
       DEALLOCATE(ephemerides, orb_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
               "Could not deallocate memory (20)", 1)
          DEALLOCATE(residuals, stat=err)
          DEALLOCATE(mask, stat=err)
          DEALLOCATE(dates_orig, stat=err)
          DEALLOCATE(dates_add, stat=err)
          DEALLOCATE(orb_arr, stat=err)         
          DEALLOCATE(ephemerides, stat=err)         
          RETURN
       END IF
       IF (info_verb >= 2) THEN
          WRITE(stdout,"(2X,A,1X,A)") "constrainRangeDistributions:", &
               "Re-calculating the rho1 and rho2 distributions... done"
       END IF
    END IF
    DEALLOCATE(dates_orig, dates_add, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not deallocate memory (25)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(dates_orig, stat=err)
       DEALLOCATE(dates_add, stat=err)
       RETURN
    END IF

    ! Rewrite sor_rho_arr_cmp (N.B. The size of the 1st dimension of
    ! sor_rho_arr_cmp is no longer necessary equal to norb!)
    ALLOCATE(sor_rho_arr_cmp(COUNT(mask),2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not allocate memory (25)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(sor_rho_arr_cmp, stat=err)
       RETURN
    END IF
    IF (info_verb >= 3) THEN
       WRITE(stdout,"(2X,A,1X,A)") "constrainRangeDistributions:", &
            "Constrained (rho1,rho2-rho1) distribution [au]: "
    END IF
    j = 0
    DO i=1,SIZE(this%sor_rho_arr_cmp,dim=1)
       IF (mask(i)) THEN
          IF (info_verb >= 3) THEN
             WRITE(stdout,"(2X,2(1X,F11.7))") &
                  this%sor_rho_arr_cmp(i,1), &
                  this%sor_rho_arr_cmp(i,2)-this%sor_rho_arr_cmp(i,1)
          END IF
          j = j + 1
          sor_rho_arr_cmp(j,:) = this%sor_rho_arr_cmp(i,:)
       END IF
    END DO
    DEALLOCATE(this%sor_rho_arr_cmp, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not deallocate memory (30)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(sor_rho_arr_cmp, stat=err)
       RETURN
    END IF
    ALLOCATE(this%sor_rho_arr_cmp(COUNT(mask),2), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not allocate memory (30)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       RETURN
    END IF
    this%sor_rho_arr_cmp = sor_rho_arr_cmp
    DEALLOCATE(sor_rho_arr_cmp, residuals, mask, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("StochasticOrbit / constrainRangeDistributions", &
            "Could not deallocate memory (35)", 1)
       DEALLOCATE(residuals, stat=err)
       DEALLOCATE(mask, stat=err)
       DEALLOCATE(sor_rho_arr_cmp, stat=err)
       RETURN
    END IF

  END SUBROUTINE constrainRangeDistributions





  SUBROUTINE toCartesian_SO(this, frame)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)  :: this
    CHARACTER(len=*), INTENT(in) :: frame

    CHARACTER(len=FRAME_LEN) :: frame_
    REAL(bp), DIMENSION(:), POINTER :: pdf_arr => NULL()
    INTEGER :: i, err

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
       IF (exist(this%orb_ml_cmp) .AND. containsDiscretePDF(this)) THEN
          pdf_arr => getDiscretePDF(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCartesian", &
                  "TRACE BACK (5)", 1)
             RETURN
          END IF
          i = MAXLOC(pdf_arr,dim=1)
          DEALLOCATE(pdf_arr, stat=err)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.containsDiscretePDF(this)) THEN
          CALL toCartesian(this%orb_ml_cmp, frame_)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCartesian", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
       END IF
    ELSE IF (this%element_type_prm == "cometary" .OR. &
         this%element_type_prm == "keplerian") THEN
       IF (ASSOCIATED(this%orb_arr_cmp) .AND. &
            ASSOCIATED(this%pdf_arr_cmp) .AND. &
            ASSOCIATED(this%jac_arr_cmp)) THEN
          DO i=1,SIZE(this%orb_arr_cmp)
             CALL toCartesian(this%orb_arr_cmp(i), frame_)
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / toCartesian", &
                     "TRACE BACK (15)", 1)
                RETURN
             END IF
             !this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i)/this%jac_arr_cmp(i,2)
          END DO
          pdf_arr => getDiscretePDF(this, "cartesian")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCartesian", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
          this%pdf_arr_cmp = pdf_arr
          DEALLOCATE(pdf_arr)
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "cartesian", frame_)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCartesian", &
                  "TRACE BACK (25)", 1)
             RETURN
          END IF
          this%cov_type_prm = "cartesian"
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          pdf_arr => getDiscretePDF(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCartesian", &
                  "TRACE BACK (30)", 1)
             RETURN
          END IF
          i = MAXLOC(pdf_arr,dim=1)
          DEALLOCATE(pdf_arr, stat=err)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
          CALL toCartesian(this%orb_ml_cmp, frame_)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCartesian", &
                  "TRACE BACK (35)", 1)
             RETURN
          END IF
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

    REAL(bp), DIMENSION(:), POINTER :: pdf_arr => NULL()
    INTEGER :: i, err

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
          pdf_arr => getDiscretePDF(this, "cometary")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCometary", &
                  "TRACE BACK (5)", 1)
             RETURN
          END IF
          this%pdf_arr_cmp = pdf_arr
          DEALLOCATE(pdf_arr, stat=err)
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "cometary")
          this%cov_type_prm = "cometary"
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          pdf_arr => getDiscretePDF(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toCometary", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
          i = MAXLOC(pdf_arr,dim=1)
          DEALLOCATE(pdf_arr, stat=err)
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

    REAL(bp), DIMENSION(:), POINTER :: pdf_arr => NULL()
    INTEGER :: i, err

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
             IF (error) THEN
                CALL errorMessage("StochasticOrbit / toKeplerian", &
                     "TRACE BACK (5)", 1)
                RETURN
             END IF
             !this%pdf_arr_cmp(i) = this%pdf_arr_cmp(i)*this%jac_arr_cmp(i,2)
          END DO
          pdf_arr => getDiscretePDF(this, "keplerian")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toKeplerian", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
          this%pdf_arr_cmp = pdf_arr
          DEALLOCATE(pdf_arr, stat=err)
       END IF
       IF (ASSOCIATED(this%cov_ml_cmp)) THEN
          this%cov_ml_cmp = getCovarianceMatrix(this, "keplerian")
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toKeplerian", &
                  "TRACE BACK (15)", 1)
             RETURN
          END IF
          this%cov_type_prm = "keplerian"
       END IF
       IF (exist(this%orb_ml_cmp) .AND. ASSOCIATED(this%pdf_arr_cmp)) THEN
          pdf_arr => getDiscretePDF(this)
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toKeplerian", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
          i = MAXLOC(pdf_arr,dim=1)
          DEALLOCATE(pdf_arr, stat=err)
          CALL NULLIFY(this%orb_ml_cmp)
          this%orb_ml_cmp = copy(this%orb_arr_cmp(i))
       ELSE IF (exist(this%orb_ml_cmp) .AND. .NOT.ASSOCIATED(this%pdf_arr_cmp)) THEN
          CALL toKeplerian(this%orb_ml_cmp)          
          IF (error) THEN
             CALL errorMessage("StochasticOrbit / toKeplerian", &
                  "TRACE BACK (25)", 1)
             RETURN
          END IF
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
  !! Updates: 
  !!      - range bounds, currently from the 3-sigma cutoff values 
  !!         of the range probability density.
  !!      - Minimum chi2
  !!      - Optional outlier detection: recognizes outlier observations
  !!         and optionally flags them (obs_masks -> .false.).
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE updateRanging(this, automatic)

    IMPLICIT NONE
    TYPE(StochasticOrbit), INTENT(inout)  :: this
    LOGICAL, INTENT(in), OPTIONAL         :: automatic

    CHARACTER(len=128) :: str, frmt
    REAL(bp), DIMENSION(:,:), POINTER :: stdevs => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: ra_mean, dec_mean
    REAL(bp) :: ra_std, dec_std
    INTEGER, DIMENSION(:), ALLOCATABLE :: nr_array
    INTEGER :: i, nobs, err, nr_of_omitted
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks => NULL()
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
            std_dev=ra_std, errstr=errstr)
       IF (LEN_TRIM(errstr) /= 0) THEN
          error = .TRUE.
          CALL errorMessage("StochasticOrbit / updateRanging", &
               "Error in moment computation for RA " // &
               TRIM(errstr), 1)
          RETURN
       END IF
       CALL moments(this%res_arr_cmp(:,i,3), mean=dec_mean(i), &
            std_dev=dec_std, errstr=errstr)
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

    ! Update minimum chi2 value and generating windows
    this%chi2_min_prm = this%chi2_min_cmp

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
