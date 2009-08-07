!
! LSST Data Management System
! Copyright 2008, 2009 LSST Corporation.
!
! This product includes software developed by the
! LSST Project (http://www.lsst.org/).
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the LSST License Statement and
! the GNU General Public License along with this program.  If not,
! see <http://www.lsstcorp.org/LegalNotices/>.
!
! 
! OpenOrb library
! F. Pierfederici <fpierfed@gmail.com>
! 
! Description
! 
! Notes
! error is a boolean global and is defined somewhere in liboorb.so!
! 
! 
! 
! -------------------------------- IMPORTANT -----------------------------------
! All angles coming in and out are in radians. 
! All times coming in and out in MJD TAI
! -------------------------------- IMPORTANT -----------------------------------
! 
! 
! The way I find out the max value used for errorCodes is:
! $NF is the number of fields in the record and so it points to the last column.
! awk '/errorCode = / { if (x < $NF) x = $NF} END { print x }' *.f90

module oorb
  use Base_cl
  ! use File_cl
  ! use PhysicalParameters_cl
  use Time_cl
  use CartesianCoordinates_cl
  use SphericalCoordinates_cl
  use Observatories_cl
  use Orbit_cl
  use Observation_cl
  use Observations_cl
  use StochasticOrbit_cl
  use planetary_data
  use utilities
  use io

  implicit none
  save
  private
  type (Observatories), public            :: obsies
  ! Orbital elements that we support.
  character(len=11), dimension(4), public :: ORBITAL_ELEMENTS = (/             &
                                                                 "keplerian  ",&
                                                                 "delaunay   ",&
                                                                 "poincare   ",&
                                                                 "equinoctial" &
                                                                 /)
  ! One can change the internal timescale to something other than TAI but let's
  ! be quiet about it ;-)
  character(len=3), public                 :: internal_timescale = "TAI"
  
  ! Public API
  public ::                                 &
       ! These are the functions we provide here.
       init,                                &
       propagateOrbit,                      &
       ranging,                             &
       ! mcmc,                                &
       lsl,                                 &
       ephemeris,                           &
       classification,                      &
       computeMoid,                         &
       dumpRangingOrbits,                   &
       dumpLslOrbit,                        &
       dumpEphemeris,                       &
       readOpenOrbOrbit,                    &
       readObsFile,                         &
       obssFromCoords,                      &
       exportRangingOrbits,                 &
       exportLslOrbit,                      &
       calendarDateToMjd,                   &
       mjdConvert,                          &
       rangingOrbitsToStochasticOrbit,      &
       lslOrbitToStochasticOrbit,           &
       exportEphemeris
       
       

  ! Constants and global variables.


  ! Code!
contains
  subroutine init(ephemFileName, verbosity, errorCode)
    ! Initialize the OpenOrb module.
    !
    ! @param ephemFileName: full path of the JPL ephem file (usually 
    !        $OORB_DATA/JPL_ephemeris/de405.dat).
    ! @param verbosity: verbosity level for OpenOrb calls [0, 5] (default 0)
    ! @return errorCode: int error code. 0 = success, otherwise, failure.

    ! Input/Output variables.
    character(len=*), intent(IN)    :: ephemFileName
    integer, optional, intent(in)   :: verbosity
    integer, intent(OUT)            :: errorCode

    ! Variable declaration.
    logical                         :: error = .false.
    type (Time)                     :: t

    
    ! Defaults.
    error = .false.
    errorCode = 0
    
    ! Verbosity levels.
    if(present(verbosity) .and. verbosity >= 0 .and. verbosity <= 5) then 
        info_verb = verbosity
        err_verb = verbosity
    else
        info_verb = 0
        err_verb = 0
    end if

    ! Init global variable OORB_DATA.
    CALL setAccessToDataFiles()
    if(error) then
        ! CALL errorMessage("oorb / init", "TRACE BACK (1)", 1)
        errorCode = 1
        return
    end if

    ! Init JPL ephemeris.
    call JPL_ephemeris_init(error, filename=trim(ephemFileName))
    if(error) then
        ! Error Could not initialize planetary ephemerides.
        CALL errorMessage("oorb / init", "TRACE BACK (2)", 1)
        errorCode = 1
        return
    end if

    ! Read data from ET-UT.dat, TAI-UTC.dat, OBSCODE.dat
    call new(t)
    call nullify(t)
    call new(obsies)
    if(error) then
        ! Error Could not initialize observatory codes.
        errorCode = 45
        return
    end if
    return
  end subroutine init

  
  subroutine ranging(obs_in, &
       element_type, &
       epoch_mjd, &
       dyn_model, &
       ! From here below it is either config file params or output
       integration_step, &
       perturbers, &
       apriori_a_max, &
       apriori_a_min, &
       apriori_periapsis_max, &
       apriori_periapsis_min, &
       apriori_apoapsis_max, &
       apriori_apoapsis_min, &
       apriori_rho_max, &
       apriori_rho_min, &
       outlier_rejection, &
       outlier_multiplier, &
       sor_type_prm, &
       sor_2point_method, &
       sor_2point_method_sw, &
       sor_norb, &
       sor_norb_sw, &
       sor_ntrial, &
       sor_ntrial_sw, &
       sor_niter, &
       sor_rho_init, &
       sor_genwin_multiplier, &
       sor_genwin_offset, &
       accwin_multiplier, &
       gaussian_rho, &
       pdf_ml_init, &
       uniform, &
       regularized, &
       random_obs, &
       outOrbit, &
       errorCode)
    ! Orbital inversion using statistical orbital ranging, that is,
    ! without making any assumptions on the shape of the resulting
    ! orbital-element pdf.

    ! Input/Output variables.
    ! Observations (i.e. DiASources) to be used as input.
    TYPE (Observations), intent(inout)                  :: obs_in
    ! Element type used in the orbit computation (default keplerian)
    CHARACTER(len=*), optional, intent(in)              :: element_type
    ! Computation epoch TAI (MJD) 
    ! If left unspecified and doing inversion, the midnight closest to the 
    ! observational mid-epoch will be used.
    REAL(bp), optional, intent(in)                      :: epoch_mjd
    ! Dynamical model [ "2-body" | "n-body" ]. Defaults to "2-body".
    character(len=6), optional, intent(in)              :: dyn_model
    ! Integrator step length (in days)
    real(bp), intent(in)                                :: integration_step
    LOGICAL, DIMENSION(:), pointer, intent(in)          :: perturbers
    ! BAYESIAN (INFORMATIVE) A PRIORI PARAMETERS
    ! Upper limit for semimajor axis in AU
    real(bp), optional, intent(in)                      :: apriori_a_max
    ! Lower limit for semimajor axis in AU (default: r_Sun = 0.00465424 AU)
    real(bp), optional, intent(in)                      :: apriori_a_min
    ! Upper limit for perihelion distance in AU
    real(bp), optional, intent(in)                      :: apriori_periapsis_max
    ! Lower limit for perihelion distance in AU
    real(bp), optional, intent(in)                      :: apriori_periapsis_min
    ! Upper limit for aphelion distance in AU
    real(bp), optional, intent(in)                      :: apriori_apoapsis_max
    ! Lower limit for aphelion distance in AU
    real(bp), optional, intent(in)                      :: apriori_apoapsis_min
    ! Upper limit for rho in AU (default: not defined)
    real(bp), optional, intent(in)                      :: apriori_rho_max
    ! Lower limit for rho in AU (default: 10*r_Earth = 0.000425641 AU)
    real(bp), optional, intent(in)                      :: apriori_rho_min
    ! Toggle outlier rejection (default false)
    logical, optional, intent(in)                       :: outlier_rejection
    ! Outlier criterion (sigma multiplier) (default 3.0)
    real(bp), optional, intent(in)                      :: outlier_multiplier
    ! Type of ranging (1=basic, 2=automatic, 3=stepwise)
    integer, intent(in)                                 :: sor_type_prm
    ! Method used for 2-point boundary-value problem. 
    ! [ continued fraction | p-iteration | n-body amoeba ]
    character(len=*),intent(in)                         :: sor_2point_method
    ! Method for solving the two-point boundary-value problem during the 
    ! preliminary steps of stepwise Ranging
    ! [ continued fraction | p-iteration | n-body amoeba ]
    character(len=*),intent(in)                         :: sor_2point_method_sw
    ! Number of requested sample orbits
    integer, intent(in)                                 :: sor_norb
    ! Number of requested sample orbits for the preliminary steps of stepwise 
    ! ranging
    integer, intent(in)                                 :: sor_norb_sw
    ! Maximum number of trial orbits
    integer, intent(in)                                 :: sor_ntrial
    ! Maximum number of trial orbits for the preliminary steps of stepwise 
    ! ranging
    integer, intent(in)                                 :: sor_ntrial_sw
    ! Maximum number of iterations in automatic/stepwise Ranging
    integer, intent(in)                                 :: sor_niter
    ! [Lower topocentric range bound corresponding to the first date 
    !  (observer -> target) [AU],
    !  Upper topocentric range bound corresponding to the first date
    !  (observer -> target) [AU],
    !  Lower topocentric range bound of the second date relative to the 
    !  generated range of the first date [AU],
    !  Upper topocentric range bound of the second date relative to the 
    !  generated range of the first date [AU]]
    ! (observer -> target) [AU]
    REAL(bp), DIMENSION(4), intent(inout)               :: sor_rho_init
    ! Sigma multiplier for generation windows
    real(bp), intent(in)                                :: sor_genwin_multiplier
    ! Offset for generation windows RA1, Dec1, RA2, Dec2 [asec]
    REAL(bp), DIMENSION(4), intent(in)                  :: sor_genwin_offset
    ! Sigma multiplier for acceptance windows
    real(bp), intent(in)                                :: accwin_multiplier
    ! sor.rho.gauss
    logical, optional, intent(in)                       :: gaussian_rho
    ! pdf.init (default -1.0)
    real(bp), optional, intent(in)                      :: pdf_ml_init
    ! Toggle use of a uniform p.d.f. (default off)
    logical, optional, intent(in)                       :: uniform
    ! Toggle use of regularization (default off)
    logical, optional, intent(in)                       :: regularized
    ! Toggle use of random observation pair. Default is off, that is, use fixed
    ! pair (cronologically first and last) (default false)
    logical, optional, intent(in)                       :: random_obs
    ! Output Sthocastic orbit.
    TYPE (StochasticOrbit), intent(out)                 :: outOrbit
    ! Output error code
    integer, intent(out)                                :: errorCode

    ! Internal variables.
    ! Integrator (only used if dynamical model is different from 2-body).
    ! [ "bulirsch-stoer" ]
    character(len=INTEGRATOR_LEN)               :: integrator = "bulirsch-stoer"!
    CHARACTER(len=DESIGNATION_LEN)              :: id
    real(bp)                                    :: dt
    TYPE (Time)                                 :: t
    ! Element type to be used in computations.
    CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
    integer                                     :: j = 0
    ! Bounds to be iterated in automatical versions 
    ! (rho1_lo, rho1_hi, rho2_lo, rho2_hi)
    LOGICAL, DIMENSION(4)                       :: sor_iterate_bounds
    TYPE (Time)                                 :: epoch
    character(len=6)                            :: dynModel = "2-body"
    real(bp)                                    :: aprioriAMax = -1.0_bp
    real(bp)                                    :: aprioriAMin = 0.00465424_bp
    real(bp)                                    :: aprioriPeriapsisMax = -1.0_bp
    real(bp)                                    :: aprioriPeriapsisMin = -1.0_bp
    real(bp)                                    :: aprioriApoapsisMax = -1.0_bp
    real(bp)                                    :: aprioriApoapsisMin = -1.0_bp
    real(bp)                                    :: aprioriRhoMax = -1.0_bp
    real(bp)                                    :: aprioriRhoMin = -1.0_bp
    logical                                     :: outlierRejection = .false.
    real(bp)                                    :: outlierMultiplier = 3.0_bp
    logical                                     :: gaussianRho = .false.
    real(bp)                                    :: pdfMLInit = -1.0_bp
    logical                                     :: uniformPDF = .false.
    logical                                     :: regularizedPDF = .false.
    logical                                     :: randomObs = .false.
    
!    TYPE (SphericalCoordinates), DIMENSION(:), POINTER  :: scoords
!    REAL(8), DIMENSION(:,:,:), POINTER                  :: cov_matrices
!    integer                                             :: i = 0
!    type (Time)                                         :: scoords_t
!    CHARACTER(len=DESIGNATION_LEN)                      :: obsId = " "
!    type (Observation)                                  :: obs_in_i
!    type (Time)                                         :: obs_t
    
    
!    cov_matrices => getCovarianceMatrices(obs_in)
!    scoords => getObservationSCoords(obs_in)
!    obsId = getID(obs_in)
!    
!    write(*, *) "ID: ", obsId
!    DO i=1,size(scoords)
!        obs_in_i = getObservation(obs_in, i)
!        scoords_t = getTime(scoords(i))
!        if(error) then
!            write(*, *) "Error getting Time object from SphCoords."
!        end if
!        obs_t = getTime(obs_in_i)
!        
!        write(*, *) "mag: ", getMagnitude(obs_in_i)
!        write(*, *) "RA:  ", getRA(obs_in_i)
!        write(*, *) "Dec: ", getDec(obs_in_i)
!        write(*, *) "MJD: ", getMJD(obs_t, internal_timescale)
!        write(*, *) "RA:  ", getLongitude(scoords(i))
!        write(*, *) "Dec: ", getLatitude(scoords(i))
!        write(*, *) "MJD: ", getMJD(scoords_t, internal_timescale)
!        write(*, *) "COV: ", cov_matrices(i,:,:)
!        call nullify(scoords_t)
!    END DO
!    write(*, *) "element_type: ", element_type
!    write(*, *) "dyn_model: ", dyn_model
!    write(*, *) "integration_step: ", integration_step
!    write(*, *) "perturbers: ", perturbers
!    write(*, *) "apriori_a_max: ", apriori_a_max
!    write(*, *) "apriori_a_min: ", apriori_a_min
!    write(*, *) "apriori_rho_min: ", apriori_rho_min
!    write(*, *) "apriori_rho_max: ", apriori_rho_max
!    write(*, *) "outlier_rejection: ", outlier_rejection
!    write(*, *) "outlier_multiplier: ", outlier_multiplier
!    write(*, *) "sor_type_prm: ", sor_type_prm
!    write(*, *) "sor_2point_method: ", sor_2point_method
!    write(*, *) "sor_2point_method_sw: ", sor_2point_method_sw
!    write(*, *) "sor_norb: ", sor_norb
!    write(*, *) "sor_norb_sw: ", sor_norb_sw
!    write(*, *) "sor_ntrial: ", sor_ntrial
!    write(*, *) "sor_ntrial_sw: ", sor_ntrial_sw
!    write(*, *) "sor_niter: ", sor_niter
!    write(*, *) "sor_rho_init: ", sor_rho_init
!    write(*, *) "sor_genwin_multiplier: ", sor_genwin_multiplier
!    write(*, *) "sor_genwin_offset: ", sor_genwin_offset
!    write(*, *) "accwin_multiplier: ", accwin_multiplier
!    write(*, *) "gaussian_rho: ", gaussian_rho
!    write(*, *) "pdf_ml_init: ", pdf_ml_init
!    write(*, *) "uniform: ", uniform
!    write(*, *) "regularized: ", regularized
!    write(*, *) "random_obs: ", random_obs
    
    ! Init optional and output vars.
    errorCode = 0
    error = .false.
    if(present(element_type)) then
        elementType = element_type
    end if
    if(.not.present(epoch_mjd)) then
        call nullify(epoch)
    else
        call NEW(epoch, epoch_mjd, internal_timescale)
    end if
    if(present(dyn_model) .and. &
       dyn_model .ne. "2-body" .and. &
       dyn_model .ne. "n-body") then
        ! Error: unsupported dynamical model.
        errorCode = 21
        return
    elseif(dyn_model .eq. "n-body") then
        dynModel = "n-body"
    end if
    if(present(apriori_a_max)) then
        aprioriAMax = apriori_a_max
    end if
    if(present(apriori_a_min)) then
        aprioriAMin = apriori_a_min
    end if
    if(present(apriori_periapsis_max)) then
        aprioriPeriapsisMax = apriori_periapsis_max
    end if
    if(present(apriori_periapsis_min)) then
        aprioriPeriapsisMin = apriori_periapsis_min
    end if
    if(present(apriori_apoapsis_max)) then
        aprioriApoapsisMax = apriori_apoapsis_max
    end if
    if(present(apriori_apoapsis_min)) then
        aprioriApoapsisMin = apriori_apoapsis_min
    end if
    if(present(apriori_rho_max)) then
        aprioriRhoMax = apriori_rho_max
    end if
    if(present(apriori_rho_min)) then
        aprioriRhoMin = apriori_rho_min
    end if
    if(present(outlier_rejection)) then
        outlierRejection = outlier_rejection
    end if
    if(present(outlier_multiplier)) then
        outlierMultiplier = outlier_multiplier
    end if
    if(present(gaussian_rho)) then
        gaussianRho = gaussian_rho
    end if
    if(present(pdf_ml_init)) then
        pdfMLInit = pdf_ml_init
    end if
    if(present(uniform)) then
        uniformPDF = uniform
    end if
    if(present(regularized)) then
        regularizedPDF = regularized
    end if
    if(present(random_obs)) then
        randomObs = random_obs
    end if
    
    ! Init global error variables.
    error  = .FALSE.

    ! Init other variables.
    j = 0
    
    ! sor_iterate_bounds
    sor_iterate_bounds(1) = .true.
    sor_iterate_bounds(2) = .true.
    sor_iterate_bounds(3) = .true.
    sor_iterate_bounds(4) = .true.
    
    id = getID(obs_in)
    IF (error) THEN
       ! Error in getID()
       errorCode = 2
       return
    END IF

    dt = getObservationalTimespan(obs_in)
    IF (error) THEN
       ! Error in getObservationalTimeInterval()
       errorCode = 3
       return
    END IF
    
    ! write(*, *) "dt: ", dt
    

    IF (sor_rho_init(3) > HUGE(sor_rho_init(3))/2) THEN
       ! Initialize the rho2-rho1 range based on the observational timespan.
       IF (dt > 10) THEN
          sor_rho_init(3) = -0.5_bp
       ELSE
          sor_rho_init(3) = -0.05_bp*dt
       END IF
       sor_rho_init(4) = -sor_rho_init(3)
    END IF
    
    ! write(*, *) "sor_rho_init: ", sor_rho_init
    
    
    IF (.NOT.exist(epoch)) THEN
        call epochFromObservations(obs_in, t, errorCode)
        if(errorCode /= 0) then
            return
        end if
    ELSE
       CALL NULLIFY(t)
       t = copy(epoch)
    END IF
    CALL NEW(outOrbit, obs_in)
    IF (error) THEN
       ! Error in new(outOrbit)
       errorCode = 8
       return
    END IF
    CALL setParameters(outOrbit, dyn_model=dynModel, perturbers=perturbers, &
         integrator=integrator, integration_step=integration_step, &
         outlier_rejection=outlierRejection, &
         outlier_multiplier=outlierMultiplier, t_inv=t, &
         element_type=elementType, &
         regularized_pdf=regularizedPDF, uniform_pdf=uniformPDF, &
         pdf_ml=pdfMLInit, accept_multiplier=accwin_multiplier, &
         apriori_a_max=aprioriAMax, apriori_a_min=aprioriAMin, &
         apriori_periapsis_max=aprioriPeriapsisMax, &
         apriori_periapsis_min=aprioriPeriapsisMin, &
         apriori_apoapsis_max=aprioriApoapsisMax, &
         apriori_apoapsis_min=aprioriApoapsisMin, &
         apriori_rho_min=aprioriRhoMin, &
         sor_2point_method=sor_2point_method, &
         sor_2point_method_sw=sor_2point_method_sw, sor_norb=sor_norb, &
         sor_ntrial=sor_ntrial, sor_rho1_l=sor_rho_init(1), &
         sor_rho1_u=sor_rho_init(2), sor_rho2_l=sor_rho_init(3), &
         sor_rho2_u=sor_rho_init(4), sor_iterate_bounds=sor_iterate_bounds, &
         sor_random_obs_selection=randomObs, gaussian_pdf=gaussianRho, &
         sor_generat_multiplier=sor_genwin_multiplier, &
         sor_generat_offset=sor_genwin_offset)
    IF (error) THEN
       ! Error in setParameters()
       errorCode = 9
       return
    END IF
    
!    write(*, *) "dynModel: ", dynModel
!    write(*, *) "perturbers: ", perturbers
!    write(*, *) "integrator: ", integrator
!    write(*, *) "integration_step: ", integration_step
!    write(*, *) "outlierRejection: ", outlierRejection
!    write(*, *) "outlierMultiplier: ", outlierMultiplier
!    write(*, *) "t (MJD TAI): ", getMjd(t, internal_timescale)
!    write(*, *) "elementType: ", elementType
!    write(*, *) "regularizedPDF: ", regularizedPDF
!    write(*, *) "uniformPDF: ", uniformPDF
!    write(*, *) "pdfMLInit: ", pdfMLInit
!    write(*, *) "accwin_multiplier: ", accwin_multiplier
!    write(*, *) "aprioriAMax: ", aprioriAMax
!    write(*, *) "aprioriAMin: ", aprioriAMin
!    write(*, *) "aprioriPeriapsisMax: ", aprioriPeriapsisMax
!    write(*, *) "aprioriPeriapsisMin: ", aprioriPeriapsisMin
!    write(*, *) "aprioriApoapsisMax: ", aprioriApoapsisMax
!    write(*, *) "aprioriApoapsisMin: ", aprioriApoapsisMin
!    write(*, *) "aprioriRhoMin: ", aprioriRhoMin
!    write(*, *) "sor_2point_method: ", sor_2point_method
!    write(*, *) "sor_2point_method_sw: ", sor_2point_method_sw
!    write(*, *) "sor_norb: ", sor_norb
!    write(*, *) "sor_ntrial: ", sor_ntrial
!    write(*, *) "sor_rho_init(1): ", sor_rho_init(1)
!    write(*, *) "sor_rho_init(2): ", sor_rho_init(2)
!    write(*, *) "sor_rho_init(3): ", sor_rho_init(3)
!    write(*, *) "sor_rho_init(4): ", sor_rho_init(4)
!    write(*, *) "sor_iterate_bounds: ", sor_iterate_bounds
!    write(*, *) "randomObs: ", randomObs
!    write(*, *) "gaussianRho: ", gaussianRho
!    write(*, *) "sor_genwin_multiplier: ", sor_genwin_multiplier
!    write(*, *) "sor_genwin_offset: ", sor_genwin_offset
    
    

    SELECT CASE (sor_type_prm)
    CASE (1)
       CALL statisticalRanging(outOrbit)
    CASE (2)
       CALL setParameters(outOrbit, sor_niter=sor_niter)
       IF (error) THEN
          ! Error in setParameters()
          errorCode = 10
          return
       END IF
       CALL autoStatisticalRanging(outOrbit)
    CASE (3)
       CALL setParameters(outOrbit, &
            sor_norb_sw=sor_norb_sw, &
            sor_ntrial_sw=sor_ntrial_sw, &
            sor_niter=sor_niter)
       IF (error) THEN
          ! Error in setParameters()
          errorCode = 11
          return
       END IF
       CALL stepwiseRanging(outOrbit, nobs_max=-1)
    CASE default
       ! Error: unknown type of ranging
       errorCode = 12
       return
    END SELECT

    ! Did we get any error?
    IF (error) THEN
       ! Error
       errorCode = 13
       error  = .FALSE.
       CALL NULLIFY(outOrbit)
       return
    end if

    ! Cleanup everything
    return
  end subroutine ranging


  subroutine mcmc
    write (*,*)  "mcmc"
    return
  end subroutine mcmc


  subroutine lsl(obs_in, &
       element_type, &
       epoch_mjd, &
       dyn_model, &
       dyn_model_init, &
       ! From here below it is either config file params or output
       integration_step, &
       integration_step_init, &
       perturbers, &
       ls_correction_factor, &
       ls_element_mask, &
       ls_niter_major_max, &
       ls_niter_major_min, &
       ls_niter_minor, &
       outlier_rejection, &
       outlier_multiplier, &
       accwin_multiplier, &
       inStochOrbit, &
       outStochOrbit, &
       errorCode)
    ! Orbital inversion using least squares with linearized
    ! covariances, that is, fixing the resulting shape of
    ! the orbital-element pdf to a Gaussian.

    ! Input/Output variables.
    ! Observations (i.e. DiASources) to be used as input.
    TYPE (Observations), intent(inout)                  :: obs_in
    ! Element type used in the orbit computation (default keplerian)
    CHARACTER(len=*), optional, intent(in)              :: element_type
    ! Computation epoch TAI (MJD) 
    ! If left unspecified and doing inversion, the midnight closest to the 
    ! observational mid-epoch will be used.
    REAL(bp), optional, intent(in)                      :: epoch_mjd
    ! Dynamical model [ "2-body" | "n-body" ]. Defaults to "n-body".
    character(len=6), optional, intent(in)              :: dyn_model
    ! Dynamical model for initial orbit [ "2-body" | "n-body" ]. Defaults to 
    ! "2-body".
    character(len=6), optional, intent(in)              :: dyn_model_init
    ! Integrator step length (in days)
    real(bp), intent(in)                                :: integration_step
    ! Integrator step length (in days) for initial orbit
    real(bp), intent(in)                                :: integration_step_init
    ! Perturbers array (planets and moon).
    LOGICAL, DIMENSION(:), pointer, intent(in)          :: perturbers
    ! Correction factor for the iterative solution of the least-squares
    ! problem [0:1] (default 0.2)
    real(bp), optional, intent(in)                      :: ls_correction_factor
    ! Elements to included in the correction process (indicated with T). Fixed 
    ! elements should be marked with F. (default all true)
    LOGICAL, DIMENSION(6), optional, intent(in)         :: ls_element_mask
    ! Maximum number of major iterations (default 20)
    integer, optional, intent(in)                       :: ls_niter_major_max
    ! Minimum number of major iterations (default 2)
    integer, optional, intent(in)                       :: ls_niter_major_min
    ! Number of iterations per scheme (default 10)
    integer, optional, intent(in)                       :: ls_niter_minor
    ! Toggle outlier rejection (default false)
    logical, optional, intent(in)                       :: outlier_rejection
    ! Outlier criterion (sigma multiplier) (default 3.0)
    real(bp), optional, intent(in)                      :: outlier_multiplier
    ! Sigma multiplier for acceptance windows
    real(bp), intent(in)                                :: accwin_multiplier
    ! Input Sthocastic orbit from ranging.
    TYPE (StochasticOrbit), intent(in)                  :: inStochOrbit
    ! Output Stichastic Orbit.
    TYPE (StochasticOrbit), intent(out)                 :: outStochOrbit
    ! Output error code
    integer, intent(out)                                :: errorCode
    
    ! Internal variables.
    integer                                             :: nobs = 0
    integer                                             :: norb = 0
    TYPE (Orbit), DIMENSION(:), POINTER                 :: orb_arr
    type (Time)                                         :: t
    ! Internal variables.
    ! Integrator (only used if dynamical model is different from 2-body).
    ! [ "bulirsch-stoer" ]
    character(len=INTEGRATOR_LEN)               :: integrator = "bulirsch-stoer"
    character(len=INTEGRATOR_LEN)          :: integratorInit = "bulirsch-stoer"
    INTEGER                                     :: err = 0
    ! Element type to be used in computations.
    CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
    integer                                     :: j = 0
    TYPE (Time)                                 :: epoch
    character(len=6)                            :: dynModel = "n-body"
    character(len=6)                            :: dynModelInit = "2-body"
    logical                                     :: outlierRejection = .false.
    real(bp)                                    :: outlierMultiplier = 3.0_bp
    LOGICAL, DIMENSION(6)                       :: lsElementMask = .true.
    integer                                     :: lsNIterMajorMax = 20
    integer                                     :: lsNIterMajorMin = 2
    integer                                     :: lsNIterMinor = 10
    real(bp)                                    :: lsCorrectionFactor = 0.2_bp
    
    ! Init
    call nullify(t)
    errorCode = 0
    error = .false.
    if(present(element_type)) then
        elementType = element_type
    end if
    if(.not.present(epoch_mjd)) then
        call nullify(epoch)
    else
        call NEW(epoch, epoch_mjd, internal_timescale)
    end if
    if(present(dyn_model) .and. &
       dyn_model .ne. "2-body" .and. &
       dyn_model .ne. "n-body") then
        ! Error: unsupported dynamical model.
        errorCode = 21
        return
    elseif(dyn_model .eq. "2-body") then
        dynModel = "2-body"
    end if
    if(present(dyn_model_init) .and. &
       dyn_model_init .ne. "2-body" .and. &
       dyn_model_init .ne. "n-body") then
        ! Error: unsupported dynamical model.
        errorCode = 21
        return
    elseif(dyn_model_init .eq. "n-body") then
        dynModelInit = "n-body"
    end if
    if(present(ls_correction_factor)) then
        lsCorrectionFactor = ls_correction_factor
    end if
    if(present(ls_element_mask)) then
        lsElementMask = ls_element_mask
    end if
    if(present(ls_niter_major_max)) then
        lsNIterMajorMax = ls_niter_major_max
    end if
    if(present(ls_niter_major_min)) then
        lsNIterMajorMin = ls_niter_major_min
    end if
    if(present(ls_niter_minor)) then
        lsNIterMinor = ls_niter_minor
    end if
    if(present(outlier_rejection)) then
        outlierRejection = outlier_rejection
    end if
    if(present(outlier_multiplier)) then
        outlierMultiplier = outlier_multiplier
    end if
    
    ! Make sure that we have at least 4 observations!
    nobs = getNrOfObservations(obs_in)
    IF (error) THEN
        ! Error in getNrOfObservations()
        errorCode = 22
        return
    END IF
    IF (nobs < 4) THEN
        ! Error: Too few observations!
        errorCode = 23
        return
    END IF
    
    ! Get the sample orbits from ranging.
    IF (containsSampledPDF(inStochOrbit)) THEN
       orb_arr => getSampleOrbits(inStochOrbit)
       IF (error) THEN
          ! Error in getSampleOrbits()
          errorCode = 16
          return
       END IF
       norb = SIZE(orb_arr)
    ELSE
       ALLOCATE(orb_arr(1))
       orb_arr(1) = getNominalOrbit(inStochOrbit)
       IF (error) THEN
          ! Error in getNominalOrbit()
          errorCode = 24
          return
       END IF
       norb = 1
    END IF
    IF (norb == 0) THEN
       CALL errorMessage("oorb4mops / lsl", &
            "Initial orbit not available.", 1)
       STOP
    END IF
    
    IF (.NOT.exist(epoch)) THEN
        call epochFromObservations(obs_in, t, errorCode)
        if(errorCode /= 0) then
            return
        end if
    ELSE
       t = copy(epoch)
    END IF
    
    
    
!    write(*, *) "element_type: ", elementType
!    write(*, *) "dyn_model: ", dynModel
!    write(*, *) "dyn_model_init: ", dynModelInit
!    write(*, *) "integration_step: ", integration_step
!    write(*, *) "integration_step_init: ", integration_step_init
!    write(*, *) "perturbers: ", perturbers
!    write(*, *) "ls_correction_factor: ", lsCorrectionFactor
!    write(*, *) "ls_element_mask: ", lsElementMask
!    write(*, *) "ls_niter_major_max: ", lsNIterMajorMax
!    write(*, *) "ls_niter_major_min: ", lsNIterMajorMin
!    write(*, *) "ls_niter_minor: ", lsNiterMinor
!    write(*, *) "outlier_rejection: ", outlierRejection
!    write(*, *) "outlier_multiplier: ", outlierMultiplier
!    write(*, *) "accwin_multiplier: ", accwin_multiplier
    
    
    
    ! Create the output Stochastic Orbit.
    CALL NEW(outStochOrbit, obs_in)
    IF (error) THEN
       ! Error in createing a new Stochasic Orbit
       errorCode = 8
       return
    END IF
    CALL setParameters(outStochOrbit, &
                       dyn_model=dynModel, &
                       perturbers=perturbers, &
                       integrator=integrator, &
                       integration_step=integration_step, &
                       outlier_rejection=outlierRejection, &
                       outlier_multiplier=outlierMultiplier, &
                       t_inv=t, &
                       element_type=elementType, &
                       accept_multiplier=accwin_multiplier, &
                       ls_correction_factor=lsCorrectionFactor, &
                       ls_element_mask=lsElementMask, &
                       ls_niter_major_max=lsNIterMajorMax, &
                       ls_niter_major_min=lsNIterMajorMin, &
                       ls_niter_minor=lsNiterMinor)
    IF (error) THEN
       ! Error in setParameters()
       errorCode = 9
       return
    END IF
    
    DO j=1,SIZE(orb_arr,dim=1)
       CALL setParameters(orb_arr(j), &
            dyn_model=dynModelInit, &
            perturbers=perturbers, &
            integrator=integratorInit, &
            integration_step=integration_step_init)
       IF (error) THEN
           ! Error in setParameters()
           errorCode = 9
           return
       END IF
       
       CALL propagate(orb_arr(j), t)
       IF (error) THEN
          ! Error in propagate()
          errorCode = 25
          return
       END IF
       CALL setParameters(orb_arr(j), &
            dyn_model=dynModel, &
            perturbers=perturbers, &
            integrator=integrator, &
            integration_step=integration_step)
       IF (error) THEN
           ! Error in setParameters()
           errorCode = 9
           return
       END IF
       
       CALL leastSquares(outStochOrbit, orb_arr(j))
       IF (.NOT.error) THEN
          EXIT
       ELSE IF (j < norb) THEN
          error = .FALSE.
       END IF
    END DO
    IF (error) THEN
        ! Error: least squares fit failed.
        errorCode = 26
        return
    end if
    
    ! Cleanup and quit.
    DEALLOCATE(orb_arr, stat=err)
    IF (err /= 0) THEN
        ! Error in deallocating orbit array.
        errorCode = 27
    END IF
    return
  end subroutine lsl


subroutine propagateOrbit(storb, &
                          dyn_model, &
                          integration_step, &
                          perturbers, &
                          epoch_mjd, &
                          errorCode)
    ! Input/Output variables.
    type(StochasticOrbit), intent(inout)                :: storb
    ! Dynamical model [ "2-body" | "n-body" ]. Defaults to "2-body".
    character(len=6), optional, intent(in)              :: dyn_model
    ! Integrator step length (in days)
    real(bp), intent(in)                                :: integration_step
    LOGICAL, DIMENSION(:), pointer, intent(in)          :: perturbers
    ! Propagate the orbit to epoch_mjd (TAI)
    real(bp), intent(in)                                :: epoch_mjd
    ! Output error code
    integer, intent(out)                                :: errorCode
    
    ! Internal variables.
    ! Integrator (only used if dynamical model is different from 2-body).
    ! [ "bulirsch-stoer" ]
    character(len=INTEGRATOR_LEN)               :: integrator = "bulirsch-stoer"!
    type(Time)                                  :: epoch
    character(len=6)                            :: dynModel = "2-body"
    
    
    ! Init vars
    error = .false.
    errorCode = 0
    CALL NULLIFY(epoch)
    IF (epoch_mjd > 0.0_bp) THEN
        CALL NEW(epoch, epoch_mjd, internal_timescale)
        IF (error) THEN
           ! Error in creating a new epoch
           errorCode = 39
           return
        END IF
     END IF
    IF (.NOT.exist(epoch)) THEN
        ! Error creating a viable epoch.
        errorCode = 40
        return
    END IF
    if(present(dyn_model) .and. &
       dyn_model .ne. "2-body" .and. &
       dyn_model .ne. "n-body") then
        ! Error: unsupported dynamical model.
        errorCode = 21
        return
    elseif(dyn_model .eq. "n-body") then
        dynModel = "n-body"
    end if
    
    ! Set integration parameters
    CALL setParameters(storb, &
                       dyn_model=dynModel, &
                       perturbers=perturbers, &
                       integrator=integrator, &
                       integration_step=integration_step)
    IF (error) THEN
        ! Error in setParameters()
        errorCode = 37
        return
    END IF
    
    ! Compute topocentric ephemerides
    CALL propagate(storb, epoch)
    IF (error) THEN
        ! Error in propagate()
        errorCode = 41
        return
    END IF
    return
  end subroutine propagateOrbit


  subroutine ephemeris(storb, &
                       element_type, &
                       dyn_model, &
                       integration_step, &
                       perturbers, &
                       obsy_code, &
                       timespan, &
                       step, &
                       ephemerides, &
                       cov_arr, &
                       pdfs_arr, &
                       observers, &
                       errorCode)
    ! Input/Output variables.
    type(StochasticOrbit), intent(inout)                :: storb
    ! Element type used in the orbit computation (default keplerian)
    CHARACTER(len=*), optional, intent(in)              :: element_type
    ! Dynamical model [ "2-body" | "n-body" ]. Defaults to "2-body".
    character(len=6), optional, intent(in)              :: dyn_model
    ! Integrator step length (in days)
    real(bp), intent(inout)                             :: integration_step
    LOGICAL, DIMENSION(:), pointer, intent(in)          :: perturbers
    ! IAU/MPC designated observatory code.
    CHARACTER(len=OBSY_CODE_LEN), intent(in)            :: obsy_code
    ! Range in days for ephemeris production.
    real(bp), intent(in)                                :: timespan
    ! Step in days for ephemeris production.
    real(bp), intent(inout)                             :: step
    ! Output ephem array (intent(inout)).
    TYPE (SphericalCoordinates), DIMENSION(:,:),POINTER :: ephemerides
    ! Ephem uncertainty matrices (intent(inout)).
    REAL(bp), DIMENSION(:,:,:), POINTER                 :: cov_arr
    REAL(bp), DIMENSION(:,:), POINTER                   :: pdfs_arr
    ! observatory coordinates and ephem time array (intent(inout)).
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: observers
    ! Output error code
    integer, intent(out)                                :: errorCode
    
    ! Internal variables.
    ! Integrator (only used if dynamical model is different from 2-body).
    ! [ "bulirsch-stoer" ]
    character(len=INTEGRATOR_LEN)               :: integrator = "bulirsch-stoer"!
    TYPE (Time)                                 :: t
    ! Element type to be used in computations.
    CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
    integer                                     :: j = 0
    character(len=6)                            :: dynModel = "2-body"
    integer                                     :: nstep = 0
    real(bp)                                    :: mjd_tai
    
    
    ! Init vars
    error = .false.
    errorCode = 0
    IF (step <= 0.0_bp) THEN
        nstep = 1        
    ELSE
       step = SIGN(ABS(step),timespan)
       IF (ABS(timespan) > 10.0_bp*EPSILON(timespan) .AND. &
            ABS(timespan) < ABS(step)) THEN
          step = timespan
       END IF
       nstep = NINT(timespan/step) + 1
    END IF
    integration_step = MIN(ABS(step),integration_step)
    
    if(present(element_type)) then
        elementType = element_type
    end if
    if(present(dyn_model) .and. &
       dyn_model .ne. "2-body" .and. &
       dyn_model .ne. "n-body") then
        ! Error: unsupported dynamical model.
        errorCode = 21
        return
    elseif(dyn_model .eq. "n-body") then
        dynModel = "n-body"
    end if
    
    ! Use equatorial coordinates:
    CALL toCartesian(storb, "equatorial")
    
    ! Get orbit epoch at the observatory.
    t = getTime(storb)
    mjd_tai = getMJD(t, internal_timescale)
    CALL NULLIFY(t)
    ALLOCATE(observers(nstep))
    DO j=1,nstep
        CALL NEW(t, mjd_tai+(j-1)*step, internal_timescale)
        IF (error) THEN
            ! Error in creating a new Time object.
            errorCode = 35
            return
       END IF
       ! Compute heliocentric observatory coordinates
       observers(j) = getObservatoryCCoord(obsies, obsy_code, t)
       IF (error) THEN
            ! Error in getObservatoryCCoord()
            errorCode = 36
            return
       END IF
       CALL rotateToEquatorial(observers(j))
       CALL NULLIFY(t)
    END DO
    
    ! Set integration parameters
    CALL setParameters(storb, &
                       dyn_model=dynModel, &
                       perturbers=perturbers, &
                       integrator=integrator, &
                       integration_step=integration_step)
    IF (error) THEN
        ! Error in setParameters()
        errorCode = 37
        return
    END IF
    
    ! Compute topocentric ephemerides
    CALL getEphemerides(storb, &
                        observers, &
                        ephemerides, &
                        cov_arr=cov_arr, &
                        pdfs_arr=pdfs_arr)
    IF (error) THEN
        ! Error in getEphemerides()
        errorCode = 38
        return
    END IF
    return
  end subroutine ephemeris
  

  subroutine classification(storb, group_name_arr, weight_arr, errorCode)
    ! Compute probability for an object with the input
    ! orbital-element pdf to belong to a group of asteroids. Caveat:
    ! an, ap, M are not currently taken into account, that is,
    ! e.g. Jupiter Trojans not correctly accounted for.
    ! Input/Output variables.
    type(StochasticOrbit), intent(inout)                    :: storb
    ! Class names and probability arrays.
    REAL(bp), DIMENSION(:), POINTER                         :: weight_arr
    CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER   :: group_name_arr
    integer, intent(out)                                    :: errorCode
    
    
    ! Init
    error = .false.
    errorCode = 0
    
    IF (.NOT.containsSampledPDF(storb)) THEN
        ! Error: we need sampled uncertainty information required for this task.
        errorCode = 44
        return
    END IF
    
    ! Propagate input orbital-element pdf to Keplerian-element pdf:
    CALL toKeplerian(storb)
    ! Compute weights for each class that has been defined:
    CALL getGroupWeights(storb, weight_arr, group_name_arr)
    return
  end subroutine classification
    

  subroutine computeMoid(storb, moid, errorCode)
    ! Compute Earth MOID (w/o uncertainty) for the input orbit.
    ! Input/Output variables.
    type(StochasticOrbit), intent(in)                   :: storb
    real(bp), intent(out)                               :: moid
    integer, intent(out)                                :: errorCode
    
    ! Internal variables.
    real(bp), dimension(6)                              :: elements
    type(Orbit)                                         :: orb
    type(Orbit)                                         :: ref_orb
    type(Time)                                          :: epoch
    real(bp)                                            :: mjd_tai
    integer                                             :: err = 0
    REAL(bp), DIMENSION(:,:), POINTER                   :: planeph
    
    
    ! Init
    error = .false.
    errorCode = 0
    moid = -1.0_bp
    
    ! Get the nominal orbit.
    orb = getNominalOrbit(storb)
    elements = getElements(orb, "cometary")
    IF (elements(2) > 1.0_bp) THEN
       ! Error: MOID for hyperbolic orbits currently not computed...
       errorCode = 42
       moid = -1.0_bp
       return
    end if
    
    ! Epoch of the asteroid orbit:
    epoch = getTime(orb)
    mjd_tai = getMJD(epoch, internal_timescale)
    ! Earth's osculating elements for the epoch of the asteroid orbit:
    planeph => JPL_ephemeris(mjd_tai, 3, 11, error)
    CALL NEW(ref_orb, planeph(1,:), "cartesian", "equatorial", epoch)
    ! MOID between Earth and asteroid orbits:
    moid = getMOID(orb, ref_orb)
    IF (error) THEN
       ! Error in getMOID()
       moid = -1.0_bp
       error = .FALSE.
       errorCode = 43
    END IF
    
    ! Cleanup.
    CALL NULLIFY(epoch)
    CALL NULLIFY(ref_orb)
    call nullify(orb)
    DEALLOCATE(planeph, stat=err)
    return
  end subroutine computeMoid
  
  
  subroutine epochFromObservations(obs_in, t, errorCode)
    TYPE (Observations), intent(in)             :: obs_in
    TYPE (Time), intent(out)                    :: t
    integer, intent(OUT)                        :: errorCode
    
    ! Internal variables.
    real(bp)                                    :: dt = 0
    real(bp)                                    :: mjd = 0
    TYPE (Observation)                          :: obs
    
    
    ! init
    errorCode = 0
    error = .false.
    CALL NULLIFY(t)
    call nullify(obs)
    dt = getObservationalTimespan(obs_in)
    IF (error) THEN
       ! Error in getObservationalTimeInterval()
       errorCode = 3
       return
    END IF
    obs = getObservation(obs_in,1)
    IF (error) THEN
      ! Error in getObservation()
      errorCode = 4
      return
    END IF
    t = getTime(obs)
    IF (error) THEN
      ! Error in getTime()
      errorCode = 5
      return
    END IF
    CALL NULLIFY(obs)
    mjd = getMJD(t, internal_timescale)
    IF (error) THEN
      ! Error in getMJD()
      errorCode = 6
      return
    END IF
    CALL NULLIFY(t)
    mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
    CALL NEW(t, mjd, internal_timescale)
    IF (error) THEN
      ! Error in new(t)
      errorCode = 7
      return
    END IF
    
  end subroutine epochFromObservations


    subroutine dumpRangingOrbits(storb, element_type, errorCode)
        TYPE (StochasticOrbit), intent(in)          :: storb
        CHARACTER(len=*), optional, intent(in)      :: element_type
        integer, intent(out)                        :: errorCode
        
        ! Internal variables.
        REAL(bp), DIMENSION(2,2)                    :: sor_rho_cmp
        TYPE (Orbit), DIMENSION(:), POINTER         :: orb_arr_cmp
        REAL(bp), DIMENSION(:), POINTER             :: pdf_arr_cmp
        REAL(bp), DIMENSION(:), POINTER             :: rchi2_arr_cmp
        REAL(bp), DIMENSION(:), POINTER             :: reg_apr_arr_cmp
        REAL(bp), DIMENSION(:,:), POINTER           :: jac_arr_cmp
        integer                                     :: numOutOrbits
        REAL(bp), DIMENSION(6)                      :: elements
        CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
        TYPE (Time)                                 :: t0
        real(bp)                                    :: orbMjd
        integer                                     :: j = 0
        integer                                     :: err = 0
        character(len=85)                           :: fmt
        
        integer                                     :: year
        integer                                     :: month
        real(bp)                                    :: day
        character(len=9)                            :: id = "K08K0042V"
        
        
        
        ! Init
        errorCode = 0
        error = .false.
        fmt = "(A16,6(1X,E21.14),1X,I4,1X,I2,1X,F8.5,1X,A12,1X,E18.10,1X,E18.10,1X,E18.10,1X,E18.10)"
        if(present(element_type)) then
            elementType = element_type
        end if
        
        CALL getResults(storb, sor_rho_cmp=sor_rho_cmp)
        IF (error) THEN
           write(*, *) "Error in getResults()"
           errorCode = 15
           return
        end if
        ! Get ORBITAL-ELEMENT PDF
        orb_arr_cmp => getSampleOrbits(storb)
        IF (error) THEN
           write(*, *) "Error in getSampleOrbits()"
           errorCode = 16
           return
        END IF
        pdf_arr_cmp => getPDFValues(storb)
        IF (error) THEN
           write(*, *) "Error in getPDFValues()"
           errorCode = 17
           return
        END IF
        rchi2_arr_cmp => getReducedChi2Distribution(storb)
        IF (error) THEN
           write(*, *) "Error in getReducedChi2Distribution()"
           errorCode = 17
           return
        END IF
        CALL getResults(storb, &
             reg_apr_arr=reg_apr_arr_cmp, &
             jac_arr=jac_arr_cmp)
        IF (error) THEN
           write(*, *) "Error in getResults()"
           numOutOrbits = 0
           errorCode = 18
           return
        END IF
    
        ! This is where we wrere originally writing the .orb file(s) out.
        numOutOrbits = SIZE(orb_arr_cmp,dim=1)
        DO j=1,numOutOrbits
           elements = getElements(orb_arr_cmp(j), &
                                  elementType, &
                                  frame="ecliptic")
           IF (error) THEN
              write(*, *) "Error in fetching the orbit"
              errorCode = 19
              return
           END IF
           IF (elementType == "keplerian") THEN
              elements(3:6) = elements(3:6)/rad_deg
           ELSE IF (elementType == "delaunay") THEN
              elements(1:3) = elements(1:3)/rad_deg
           ELSE IF (elementType == "poincare") THEN
              elements(4) = elements(4)/rad_deg
           ELSE IF (elementType == "equinoctial") THEN
              elements(6) = elements(6)/rad_deg
           END IF
           
           t0 = getTime(orb_arr_cmp(j))
           CALL getCalendarDate(t0, "tdt", year, month, day)
           orbMjd = getMjd(t0, internal_timescale)
           call nullify(t0)
           IF (error) THEN
              write(*, *) "Error in getTime/getMjd"
              errorCode = 5
              return
           END IF
           
           
           WRITE(*, fmt) id, elements(1:6), year, month, day, elementType, &
                 pdf_arr_cmp(j), rchi2_arr_cmp(j), reg_apr_arr_cmp(j), &
                 jac_arr_cmp(j, 1:3)
           
!           write(*, *) "elements:             ", elements(1:6)
!           write(*, *) "epoch:                ", orbMjd
!           write(*, *) "Un-normalized p.d.f.: ", pdf_arr_cmp(j)
!           write(*, *) "Reduced chi2:         ", rchi2_arr_cmp(j)
!           write(*, *) "Regularized apr:      ", reg_apr_arr_cmp(j)
!           write(*, *) "Jacobian det:         ", jac_arr_cmp(j, 1:3)
        END DO
        ! WRITE RESIDUALS
    
        ! This is where we were originally writing out the .res file.
        ! CALL writeResiduals(storb, obs_in, getUnit(out_file))
    
        ! Cleanup everything
        DO j=1,SIZE(orb_arr_cmp)
           CALL NULLIFY(orb_arr_cmp(j))
        END DO
        DEALLOCATE(orb_arr_cmp, stat=err)
        DEALLOCATE(pdf_arr_cmp, stat=err)
        DEALLOCATE(rchi2_arr_cmp, stat=err)
        DEALLOCATE(reg_apr_arr_cmp, stat=err)
        DEALLOCATE(jac_arr_cmp, stat=err)
        IF (err /= 0) THEN
           write(*, *) "Error: Could not deallocate memory"
           errorCode = 19
           return
        END IF
        return
    end subroutine dumpRangingOrbits
    
    
    subroutine dumpLslOrbit(storb, element_type, errorCode)
        TYPE (StochasticOrbit), intent(in)          :: storb
        CHARACTER(len=*), optional, intent(in)      :: element_type
        integer, intent(out)                        :: errorCode
        
        ! Internal variables.
        REAL(bp), DIMENSION(6)                      :: elements
        CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
        TYPE (Time)                                 :: t0
        real(bp)                                    :: orbMjd
        integer                                     :: j = 0
        integer                                     :: k = 0
        TYPE (Orbit)                                :: orb
        REAL(bp), DIMENSION(6,6)                    :: cov
        REAL(bp), DIMENSION(6,6)                    :: corr
        REAL(bp), DIMENSION(6)                      :: sigmas
        
        
        
        ! Init
        errorCode = 0
        error = .false.
        if(present(element_type)) then
            elementType = element_type
        end if
        
        ! Get the nominal orbit out.
        orb = getNominalOrbit(storb)
        IF (error) THEN
           write(*, *) "Error in getNominalOrbit()"
           errorCode = 28
           return
        end if
        
        ! Get the orbital elements.
        elements = getElements(orb, elementType, frame="ecliptic")
        IF (error) THEN
           write(*, *) "Error in fetching the orbit"
           errorCode = 19
           return
        END IF
        IF (elementType == "keplerian") THEN
           elements(3:6) = elements(3:6)/rad_deg
        ELSE IF (elementType == "delaunay") THEN
           elements(1:3) = elements(1:3)/rad_deg
        ELSE IF (elementType == "poincare") THEN
           elements(4) = elements(4)/rad_deg
        ELSE IF (elementType == "equinoctial") THEN
           elements(6) = elements(6)/rad_deg
        END IF
        
        ! Get orbit epoch.
        t0 = getTime(orb)
        orbMjd = getMjd(t0, internal_timescale)
        call nullify(t0)
        IF (error) THEN
           write(*, *) "Error in getTime/getMjd"
           errorCode = 5
           return
        END IF
        
        ! Covariance matrix.
        cov = getCovarianceMatrix(storb, elementType, "ecliptic")
        IF (error) THEN
           write(*, *) "Error in getCovarianceMatrix"
           errorCode = 29
           return
        END IF
        DO j=1,6
           sigmas(j) = SQRT(cov(j, j))
           IF (elementType == "keplerian" .AND. j >= 3) THEN
              sigmas(j) = sigmas(j) / rad_deg
           END IF
        END DO
        
        ! Correlations.
        DO j=1,6
           DO k=1,6
              corr(j,k) = cov(j,k) / (sigmas(j)*sigmas(k))
           END DO
        END DO
        
        
        write(*, *) "elements:             ", elements(1:6)
        write(*, *) "epoch:                ", orbMjd
        write(*, *) "Covariance:           ", cov
        write(*, *) "Sigmas:               ", sigmas
        write(*, *) "Correlation:          ", corr
    
        ! Cleanup everything
        CALL NULLIFY(orb)
        return
    end subroutine dumpLslOrbit
    
        
    subroutine dumpEphemeris(storb, ephemerides, cov, pdfs, observers, &
                             obscode, errorCode)
        ! Input/Output variables.
        type(StochasticOrbit), intent(inout)                :: storb
        ! Output ephem array (intent(inout)).
        TYPE (SphericalCoordinates), DIMENSION(:,:),POINTER :: ephemerides
        ! Ephem uncertainty matrices (intent(inout)).
        REAL(bp), DIMENSION(:,:,:), POINTER                 :: cov
        REAL(bp), DIMENSION(:,:), POINTER                   :: pdfs
        ! observatory coordinates and ephem time array (intent(inout)).
        TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: observers
        ! Observatory code.
        character(len=*), intent(in)                        :: obscode
        ! Output error code
        integer, intent(out)                                :: errorCode
        
        ! Internal vars.
        TYPE (Time)                                         :: t
        integer                                             :: j = 0
        integer                                             :: k = 0
        integer                                             :: l = 0
        real(bp)                                            :: mjd_tai = 0.0_bp
        REAL(bp), DIMENSION(6)                              :: coordinates
        REAL(bp), DIMENSION(6)                              :: stdev
        REAL(bp), DIMENSION(6,6)                            :: corr
        
        
        ! Init
        error = .false.
        errorCode = 0
        
        
        DO j=1,SIZE(observers)
            t = getTime(observers(j))
            mjd_tai = getMJD(t, internal_timescale)
            CALL NULLIFY(t)
            IF (containsSampledPDF(storb)) THEN
               ! Input orbits correspond to one or more sampled pdfs.
               ! Loop over sampled orbits:
               DO k=1,SIZE(ephemerides, dim=1)
                  ! Make sure that the ephemeris is equatorial:
                  CALL rotateToEquatorial(ephemerides(k,j))
                  coordinates = getCoordinates(ephemerides(k,j))
                  WRITE(*, *) TRIM(obscode), &
                       mjd_tai, coordinates(1), &
                       coordinates(2:3)/rad_deg, coordinates(4), &
                       coordinates(5:6)/rad_deg, pdfs(k,j)
               END DO
            ELSE
               ! Input orbits correspond to one or more single-point estimates of the pdf.
               ! Make sure that the ephemeris is equatorial:
               CALL rotateToEquatorial(ephemerides(1,j))
               coordinates = getCoordinates(ephemerides(1,j))
               DO k=1,6
                  stdev(k) = SQRT(cov(k,k,j)) 
               END DO
               DO k=1,6
                  DO l=1,6
                     corr(k,l) = cov(k,l,j) / &
                          (stdev(k)*stdev(l))
                  END DO
               END DO
               WRITE(*, "(A18,28(1X,F18.10))\") TRIM(obscode), mjd_tai, &
                    coordinates(1), &
                    coordinates(2:3)/rad_deg, coordinates(4), &
                    coordinates(5:6)/rad_deg, stdev(1), &
                    stdev(2:3)/rad_deg, stdev(4), &
                    stdev(5:6)/rad_deg, corr(1,2:6), &
                    corr(2,3:6), corr(3,4:6), corr(4,5:6), corr(5,6)
            END IF
      END DO
    end subroutine dumpEphemeris
    
    
    subroutine readOpenOrbOrbit(fileName, storb, errorCode)
        character(len=*), intent(in)                :: fileName
        TYPE (StochasticOrbit), intent(out)         :: storb
        integer, intent(out)                        :: errorCode
        
        ! Internal vars.
        integer                                     :: norb = 0
        TYPE (File)                                 :: orb_in_file
        CHARACTER(len=1024), DIMENSION(4)           :: header
        REAL(bp), DIMENSION(:,:), POINTER           :: HG_arr_in
        integer                                     :: j = 0
        CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: id_arr_in
        CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: element_type_pdf_arr_in
        REAL(bp), DIMENSION(:,:), ALLOCATABLE       :: jac_arr_in
        REAL(bp), DIMENSION(:,:,:), POINTER         :: cov_arr_in
        REAL(bp), DIMENSION(:), ALLOCATABLE         :: rchi2_arr_in, &
                                                       pdf_arr_in, &
                                                       reg_apr_arr_in
        CHARACTER(len=ELEMENT_TYPE_LEN)             :: element_type_in
        CHARACTER(len=DESIGNATION_LEN)              :: id
        TYPE (Orbit), DIMENSION(:), POINTER         :: orb_arr_in
        
        ! Init vars.
        error = .false.
        errorCode = 0
        
        CALL NEW(orb_in_file, fileName)
        IF (error) THEN
           ! Error in creating input file instance
           errorCode = 30
           return
        END IF
        CALL setStatusOld(orb_in_file)
        CALL OPEN(orb_in_file)
        IF (error) THEN
           ! Error in opeining input file.
           errorCode = 31
           return
        END IF
        ! 4 header lines are not taken into account.
        norb = getNrOfLines(orb_in_file) - 4
        IF (error) THEN
           ! Error in getNrOfLines()
           errorCode = 32
           return
        END IF
        ALLOCATE(id_arr_in(norb), orb_arr_in(norb), &
             element_type_pdf_arr_in(norb), cov_arr_in(norb,6,6), &
             HG_arr_in(norb,2), pdf_arr_in(norb), &
             rchi2_arr_in(norb), jac_arr_in(norb,3), &
             reg_apr_arr_in(norb), stat=errorCode)
        IF (errorCode /= 0) THEN
           ! Error in memory allocation
           errorCode = 33
           return
        END IF
        jac_arr_in = -1.0_bp
        pdf_arr_in = -1.0_bp
        cov_arr_in = -1.0_bp
        id_arr_in = " "
        element_type_pdf_arr_in = " "
        header(1:4)(:) = " "
        FORALL (j=1:norb)
           HG_arr_in(j,1:2) = (/ 99.0_bp, 9.9_bp /)
        END FORALL
        DO j=1,norb
           CALL readOpenOrbOrbitFile(getUnit(orb_in_file), header, &
                element_type_in=element_type_in, id=id_arr_in(j), &
                orb=orb_arr_in(j), &
                element_type_pdf=element_type_pdf_arr_in(j), &
                cov=cov_arr_in(j,:,:), H=HG_arr_in(j,1), &
                G=HG_arr_in(j,2), pdf=pdf_arr_in(j), &
                rchi2=rchi2_arr_in(j), reg_apr=reg_apr_arr_in(j), &
                jac_sph_inv=jac_arr_in(j,1), &
                jac_car_kep=jac_arr_in(j,2), &
                jac_equ_kep=jac_arr_in(j,3))
           IF (error) THEN
              ! Error in readOpenOrbOrbitFile()
               errorCode = 34
               return
           END IF
        END DO
        CALL NULLIFY(orb_in_file)

        ! Initialize stochasticorbits if uncertainty information available:
        IF (norb > 1 .AND. &
            ALL(pdf_arr_in > 0.0_bp) .AND. &
            ALL(jac_arr_in > 0.0_bp)) THEN
           call nullify(storb)
           CALL NEW(storb, orb_arr_in, pdf_arr_in, &
                element_type_pdf_arr_in(1), jac_arr=jac_arr_in, &
                reg_apr_arr=reg_apr_arr_in, &
                rchi2_arr=rchi2_arr_in)
           id = id_arr_in(1)
           DEALLOCATE(id_arr_in)
           ALLOCATE(id_arr_in(1))
           id_arr_in(1) = id
        ELSE IF (norb == 1 .AND. ALL(cov_arr_in(:,1,1) > 0.0_bp)) THEN
            CALL NEW(storb, orb_arr_in(1), cov_arr_in(1,:,:), &
                     cov_type=element_type_in, element_type=element_type_in)
        ELSE
           ! Error in StochasticOrbit instantiation
           errorCode = 34
           return
        END IF
        
        ! Cleanup.
        DO j=1,norb
           CALL NULLIFY(orb_arr_in(j))
        END DO
        DEALLOCATE(orb_arr_in, pdf_arr_in, rchi2_arr_in, jac_arr_in, &
                   reg_apr_arr_in, element_type_pdf_arr_in)
        return
    end subroutine readOpenOrbOrbit


    subroutine rangingOrbitsToStochasticOrbit(orbits, norb, storb, errorCode)
        ! Flattened orbit array. Each elememnt is
        !  (id, elements(1:6), epoch_mjd, Un-normalized p.d.f., Reduced chi2, 
        !   Regularized apr, Jacobian det(1:3), element_type_index)
        integer, intent(in)                         :: norb
        real(bp), dimension(norb, 15), intent(in)   :: orbits
        TYPE (StochasticOrbit), intent(out)         :: storb
        integer, intent(out)                        :: errorCode
        
        ! Internal vars.
        REAL(bp), DIMENSION(:,:), POINTER           :: HG_arr_in
        integer                                     :: j = 0
        integer                                     :: i = 0
        CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: id_arr_in
        CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: element_type_pdf_arr_in
        REAL(bp), DIMENSION(:,:), ALLOCATABLE       :: jac_arr_in
        REAL(bp), DIMENSION(:,:,:), POINTER         :: cov_arr_in
        REAL(bp), DIMENSION(:), ALLOCATABLE         :: rchi2_arr_in, &
                                                       pdf_arr_in, &
                                                       reg_apr_arr_in
        CHARACTER(len=DESIGNATION_LEN)              :: id
        TYPE (Orbit), DIMENSION(:), POINTER         :: orb_arr_in
        REAL(bp), DIMENSION(15)                     :: correlation
        real(bp), dimension(6)                      :: stdev
        real(bp), dimension(6)                      :: elements
        type(Time)                                  :: t
        character(len=11)                           :: element_type
        integer                                     :: element_type_index = -1
        
        
        ! Init vars.
        error = .false.
        errorCode = 0
        
        ! Get the element type from the input flattened orbits (the are supposed 
        ! be all the same).
        element_type_index = orbits(1, 15)
        if(element_type_index .le. 0 .or.                               &
           element_type_index .gt. size(ORBITAL_ELEMENTS)) then
            ! Error: unsupported orbital elements.
            errorCode = 58
            return
        end if
        element_type = ORBITAL_ELEMENTS(element_type_index)
        
        ! Allocate memory for the orbit arrays.
        ALLOCATE(id_arr_in(norb), orb_arr_in(norb), &
             element_type_pdf_arr_in(norb), cov_arr_in(norb,6,6), &
             HG_arr_in(norb,2), pdf_arr_in(norb), &
             rchi2_arr_in(norb), jac_arr_in(norb,3), &
             reg_apr_arr_in(norb), stat=errorCode)
        IF (errorCode /= 0) THEN
           ! Error in memory allocation
           errorCode = 33
           return
        END IF
        jac_arr_in = -1.0_bp
        pdf_arr_in = -1.0_bp
        cov_arr_in = 0.0_bp
        id_arr_in = " "
        ! There are no H/G in the input orbits.
        FORALL (j=1:norb)
           HG_arr_in(j,1:2) = (/ 99.0_bp, 9.9_bp /)
        END FORALL
        stdev = -1.0_bp
        correlation = 2.0_bp
        element_type_pdf_arr_in = element_type
        
        ! Get each flattened orbit and create an Orbit instance.
        do j=1,norb
            ! Just to beat on a dead horse:
            ! orbits(j,1):      id
            ! orbits(j,2:7):    elements(1:6)
            ! orbits(j,8):      epoch_mjd
            ! orbits(j,9):      Un-normalized p.d.f.
            ! orbits(j,10):     Reduced chi2
            ! orbits(j,11):     Regularized apr
            ! orbits(j,12:14):  Jacobian det(1:3)
            ! orbits(j,15):     element_type_index
            ! Convert angles to radians, if needed.
            write(id_arr_in(j), fmt="(A,I10.10)") "TRK", IDINT(orbits(j, 1))
            
            elements(1:6) = orbits(j,2:7)
            IF (element_type == "keplerian") THEN
                elements(3:6) = elements(3:6) * rad_deg
            ELSE IF (element_type == "delaunay") THEN
                elements(1:3) = elements(1:3) * rad_deg
            ELSE IF (element_type == "poincare") THEN
                elements(4) = elements(4) * rad_deg
            ELSE IF (element_type == "equinoctial") THEN
                elements(6) = elements(6) * rad_deg
            END IF
            
            ! Create a Time instance.
            call new(t, orbits(j, 8), internal_timescale)
            if(error) then
                ! Error in creating a Time instance.
                errorCode = 57
                return
            end if
            
            ! Now create an Orbit instance.
            call new(orb_arr_in(j),                                 &
                     elements(1:6),                                 &
                     element_type,                                  &
                     "ecliptic",                                    &
                     copy(t))
            call nullify(t)
            
            ! Covariance: we fake it since we do not have it.
            cov_arr_in(j,:,:) = 0.0_bp
            DO i=1,6
                cov_arr_in(j,i,i) = -1.0_bp
            END DO
            
            ! Populate the uncertinty arrays.
            pdf_arr_in(j) = orbits(j, 9)
            rchi2_arr_in(j) = orbits(j, 10)
            reg_apr_arr_in(j) = orbits(j, 11)
            jac_arr_in(j,1:3) = orbits(j, 12:14)
        end do

        ! Initialize stochasticorbits if uncertainty information available:
        IF (norb > 1 .AND. &
            ALL(pdf_arr_in > 0.0_bp) .AND. &
            ALL(jac_arr_in > 0.0_bp)) THEN
           call nullify(storb)
           CALL NEW(storb, orb_arr_in, pdf_arr_in, &
                element_type_pdf_arr_in(1), jac_arr=jac_arr_in, &
                reg_apr_arr=reg_apr_arr_in, &
                rchi2_arr=rchi2_arr_in)
           id = id_arr_in(1)
           DEALLOCATE(id_arr_in)
           ALLOCATE(id_arr_in(1))
           id_arr_in(1) = id
        ELSE IF (norb == 1 .AND. ALL(cov_arr_in(:,1,1) > 0.0_bp)) THEN
            CALL NEW(storb, orb_arr_in(1), cov_arr_in(1,:,:), &
                     cov_type=element_type, element_type=element_type)
        ELSE
           ! Error in StochasticOrbit instantiation
           errorCode = 34
           return
        END IF
        
        ! Cleanup.
        DO j=1,norb
           CALL NULLIFY(orb_arr_in(j))
        END DO
        DEALLOCATE(orb_arr_in, pdf_arr_in, rchi2_arr_in, jac_arr_in, &
                   reg_apr_arr_in, element_type_pdf_arr_in)
        return
    end subroutine rangingOrbitsToStochasticOrbit
    
    
    subroutine lslOrbitToStochasticOrbit(in_orbit,                          &
                                         in_covariance,                     &
                                         storb,                             &
                                         errorCode)
        ! Input/Output variables.
        ! Input flattened orbit:
        !  (track_id, elements(1:6), epoch, H, G, element_type_index)
        ! FIXME: H and G are not currently computed.
        real(8),dimension(11), intent(in)                   :: in_orbit
        ! Uncertainty matrices:
        REAL(8), DIMENSION(6,6), intent(in)                 :: in_covariance
        type(StochasticOrbit), intent(out)                  :: storb
        integer, intent(out)                                :: errorCode
        
        ! Internal vars.
        REAL(bp), DIMENSION(:,:), POINTER           :: HG_arr_in
        integer                                     :: j = 0
        CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: id_arr_in
        CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: element_type_pdf_arr_in
        REAL(bp), DIMENSION(:,:), ALLOCATABLE       :: jac_arr_in
        REAL(bp), DIMENSION(:,:,:), POINTER         :: cov_arr_in
        REAL(bp), DIMENSION(:), ALLOCATABLE         :: rchi2_arr_in, &
                                                       pdf_arr_in, &
                                                       reg_apr_arr_in
        CHARACTER(len=DESIGNATION_LEN)              :: id
        TYPE (Orbit), DIMENSION(:), POINTER         :: orb_arr_in
        real(bp), dimension(6)                      :: elements
        type(Time)                                  :: t
        character(len=11)                           :: element_type
        integer                                     :: element_type_index = -1
        integer                                     :: norb = 1
        
        
        ! Init vars.
        error = .false.
        errorCode = 0
        ! We try to keep these two subroutines similar so that we can merge in 
        ! the future, hence the norb = 1 thing here.
        norb = 1
        
        ! Get the element type from the input flattened in_orbit.
        element_type_index = in_orbit(11)
        if(element_type_index .le. 0 .or.                               &
           element_type_index .gt. size(ORBITAL_ELEMENTS)) then
            ! Error: unsupported orbital elements.
            errorCode = 58
            return
        end if
        element_type = ORBITAL_ELEMENTS(element_type_index)
        
        ! Allocate memory for the orbit arrays.
        ALLOCATE(id_arr_in(norb), orb_arr_in(norb), &
             element_type_pdf_arr_in(norb), cov_arr_in(norb,6,6), &
             HG_arr_in(norb,2), pdf_arr_in(norb), &
             rchi2_arr_in(norb), jac_arr_in(norb,3), &
             reg_apr_arr_in(norb), stat=errorCode)
        IF (errorCode /= 0) THEN
           ! Error in memory allocation
           errorCode = 33
           return
        END IF
        jac_arr_in = -1.0_bp
        pdf_arr_in = -1.0_bp
        cov_arr_in = 0.0_bp
        id_arr_in = " "
        ! There are no H/G in the input in_orbit.
        FORALL (j=1:norb)
           HG_arr_in(j,1:2) = (/ 99.0_bp, 9.9_bp /)
        END FORALL
        element_type_pdf_arr_in = element_type
        
        ! Get each flattened orbit and create an Orbit instance.
        do j=1,norb
            ! Just to beat on a dead horse:
            ! in_orbit(1):      id
            ! in_orbit(2:7):    elements(1:6)
            ! in_orbit(8):      epoch_mjd
            ! in_orbit(9):      H
            ! in_orbit(10):     G
            ! in_orbit(j,11):   element_type_index
            ! Convert angles to radians, if needed.
            write(id_arr_in(j), fmt="(A,I10.10)") "TRK", IDINT(in_orbit(1))
            
            elements(1:6) = in_orbit(2:7)
            IF (element_type == "keplerian") THEN
                elements(3:6) = elements(3:6) * rad_deg
            ELSE IF (element_type == "delaunay") THEN
                elements(1:3) = elements(1:3) * rad_deg
            ELSE IF (element_type == "poincare") THEN
                elements(4) = elements(4) * rad_deg
            ELSE IF (element_type == "equinoctial") THEN
                elements(6) = elements(6) * rad_deg
            END IF
            
            ! Create a Time instance.
            call new(t, in_orbit(8), internal_timescale)
            if(error) then
                ! Error in creating a Time instance.
                errorCode = 57
                return
            end if
            
            ! Now create an Orbit instance.
            call new(orb_arr_in(j),                                 &
                     elements(1:6),                                 &
                     element_type,                                  &
                     "ecliptic",                                    &
                     copy(t))
            call nullify(t)
            
            ! Covariance. We do not need to re-compute it since we already have 
            ! it.
            cov_arr_in(j,:,:) = in_covariance(:,:)
        end do

        ! Initialize stochasticorbits if uncertainty information available:
        IF (norb > 1 .AND. &
            ALL(pdf_arr_in > 0.0_bp) .AND. &
            ALL(jac_arr_in > 0.0_bp)) THEN
           call nullify(storb)
           CALL NEW(storb, orb_arr_in, pdf_arr_in, &
                element_type_pdf_arr_in(1), jac_arr=jac_arr_in, &
                reg_apr_arr=reg_apr_arr_in, &
                rchi2_arr=rchi2_arr_in)
           id = id_arr_in(1)
           DEALLOCATE(id_arr_in)
           ALLOCATE(id_arr_in(1))
           id_arr_in(1) = id
        ELSE IF (norb == 1 .AND. ALL(cov_arr_in(:,1,1) > 0.0_bp)) THEN
            CALL NEW(storb, orb_arr_in(1), cov_arr_in(1,:,:), &
                     cov_type=element_type, element_type=element_type)
        ELSE
           ! Error in StochasticOrbit instantiation
           errorCode = 34
           return
        END IF
        
        ! Cleanup.
        DO j=1,norb
           CALL NULLIFY(orb_arr_in(j))
        END DO
        DEALLOCATE(orb_arr_in, pdf_arr_in, rchi2_arr_in, jac_arr_in, &
                   reg_apr_arr_in, element_type_pdf_arr_in)
        return
    end subroutine lslOrbitToStochasticOrbit


    
    subroutine readObsFile(fileName, stdev, obss, errorCode)
        ! Input/output vars
        character(len=*), intent(in)                            :: fileName
        ! This most likely comes from a config file. If not, remember that all
        ! angles are in radians and that stdev(2) = RAErr and stdev(3) = DecErr.
        REAL(bp), DIMENSION(6), optional, intent(in)            :: stdev
        ! Remember pointers need to be intent(inout), which is the default.
        TYPE (Observations), dimension(:), pointer              :: obss
        integer, intent(out)                                    :: errorCode
        
        ! Internal vars
        real(bp), dimension(6)                                  :: stdevArray
        type(File)                                              :: obs_file
        type(Observations)                                      :: rawObss
        
        
        ! Init
        error = .false.
        errorCode = 0
        stdevArray = 0.3_bp * rad_asec
        if(present(stdev)) then
            stdevArray = stdev
        end if
        
        ! Open the file.
        CALL NEW(obs_file, TRIM(fileName))
        CALL setStatusOld(obs_file)
        CALL OPEN(obs_file)
        if (error) then
            ! Error in opening observation file.
            errorCode = 46
            return
        end if
        CALL NEW(rawObss, obs_file, stdev=stdevArray)
        CALL NULLIFY(obs_file)
        
        ! Split the raw observations into sets.
        obss => getSeparatedSets(rawObss)
        CALL NULLIFY(rawObss)
        return
    end subroutine readObsFile


    subroutine obssFromCoords(trackId, n, coords, mjds, mags, filters, &
                              obscodes, obss, errorCode)
        ! obssFromCoords
        ! 
        ! Create a list of Observations instances from the input coordinate 
        ! array.
        ! @param trackId: object ID. IDs are used to group observations of the 
        !        same object together (equivalent to track IDs in PS-MOPS). IDs 
        !        are numbers.
        ! @param n: number of coordinates to process.
        ! @param coords: array of [(ra, ra_err, dec, dec_err), ] coordinates in
        !        radians.
        ! @param mjds: array of MJDs (TAI), one per coords element. These are
        !        the MJDs of the input observations/DiaSources/Detections.
        ! @param mags: array of observed magnitudes, one per coords element.
        !        mags are floats.
        ! @param filters: array of filter names, one per coords element.
        !        Filters are (usually 1 or 2-letter) strings
        ! @param obscodes: array of MPC observatory codes, one per coords 
        !        element. Observatory codes are (up to 4 letter) strings.
        ! 
        ! @param obss: array of Observations instances, one element per unique 
        !        input ID.
        ! @param errorCode: error code. Anything other than 0 signals an error.
        integer, intent(in)                                     :: trackId
        integer, intent(in)                                     :: n
        real(bp), dimension(n,4), intent(in)                    :: coords
        real(bp), dimension(n), intent(in)                      :: mjds
        real(bp), dimension(n), intent(in)                      :: mags
        character(len=2), dimension(n), intent(in)              :: filters
        character(len=OBSY_CODE_LEN), dimension(n),intent(in)   :: obscodes
        ! Pointers have to be intent(inout)!
        TYPE (Observations), DIMENSION(:), POINTER              :: obss
        integer, intent(out)                                    :: errorCode
        
        ! Internal variables.
        type(Observations)                                      :: rawObss
        LOGICAL, DIMENSION(6)                                   :: obs_mask
        integer                                                 :: i = 0
        type(Time)                                              :: t
        TYPE (SphericalCoordinates)                             :: obs_scoord
        REAL(bp), DIMENSION(6,6)                                :: covariance
        real(bp)                                                :: ra = 0.0_bp
        real(bp)                                                :: dec = 0.0_bp
        real(bp)                                                :: raErr=0.0_bp
        real(bp)                                                :: decErr=0.0_bp
        TYPE (Observatory)                                      :: obsy
        TYPE (CartesianCoordinates)                             :: obsy_ccoord
        type(Observation)                                       :: obs
        character(len=255)                                      :: designation
        
        
        ! Init
        error = .false.
        errorCode = 0
        designation = " "
        ! We only have RA and Dec and obs_mask is
        ! (radial distance, longitude, latitude + their time derivatives), so:
        obs_mask = (/ .FALSE., .TRUE., .TRUE., .FALSE., .FALSE., .FALSE. /)
        call new(rawObss)
        if(n .le. 0) then
            ! Now enough data!
            errorCode = 50
            return
        end if
        
        ! Loop through the input arrays and create Observations!
        DO i=1,n
            ! Extract RA/Dec etc from coords.
            ra = coords(i, 1)
            raErr = coords(i, 2)
            dec = coords(i, 3)
            decErr = coords(i, 4)
            
            ! Create a Time instance using the input mjd.
            call new(t, mjds(i), internal_timescale)
            if(error) then
                ! Error in creating Time instnce.
                errorCode = 51
                return
            end if
            
            ! Create SphericalCoordinates object containing R.A. Dec. epoch.
            CALL NEW(obs_scoord, ra, dec, t)
            IF (error) THEN
                ! Error in creating a neSphericalCoordinates instance.
                errorCode = 52
                return
            END IF
            
            ! Now set the covariance matrix.
            covariance = 0.0_bp
            covariance(2,2) = raErr**2.0_bp
            covariance(3,3) = decErr**2.0_bp
            
!             write(*, *) "cov: ", covariance
            
            ! Compute the heliocentric position of the observer at epoch t:
            obsy = getObservatory(obsies, TRIM(obscodes(i)))
            IF (error) THEN
                ! Error in deriving the observatory location!
                errorCode = 48
                return
            END IF
            obsy_ccoord = getObservatoryCCoord(obsies, obsy, t)
            
            ! Create observation object:
            CALL NULLIFY(obs)
            write(designation, fmt="(A,I10.10)") "TRK", trackId
            call new(obs, &
                     number=trackId, &
                     designation=trim(designation), &
                     discovery=.false., &
                     note1=" ", &
                     note2=" ", &
                     obs_scoord=obs_scoord, &
                     covariance=covariance, &
                     obs_mask=obs_mask, &
                     mag=mags(i), &
                     filter=trim(filters(i)), &
                     obsy=obsy, &
                     obsy_ccoord=obsy_ccoord)
            IF (error) THEN
                ! Error in creating a new Observation instance.
                errorCode = 49
                return
            END IF
            
            call addObservation(rawObss, obs, sort=.true.)
            if(error) then
                ! Error in adding a new observation.
                errorCode = 53
                return
            end if
            
            call nullify(obs)
            CALL NULLIFY(obs_scoord)
            CALL NULLIFY(t)
            CALL NULLIFY(obsy)
            CALL NULLIFY(obsy_ccoord)
        END DO
        
        ! Split the raw observatikons into one list per ID.
        obss => getSeparatedSets(rawObss)
        CALL NULLIFY(rawObss)
        return
    end subroutine obssFromCoords


    subroutine exportRangingOrbits(storb, element_type, track_id, outOrbits, &
                                   errorCode)
        TYPE (StochasticOrbit), intent(in)              :: storb
        CHARACTER(len=*), optional, intent(in)          :: element_type
        integer, intent(in)                             :: track_id
        ! outOrbits has the form ((id, elements(1:6), epoch, Un-normalized pdf,
        ! Reduced chi2, Regularized apr, Jacobian det(1:3), el_type_index), )
        real(bp),dimension(:,:),intent(out)             :: outOrbits
        integer, intent(out)                            :: errorCode
        
        ! Internal variables.
        REAL(bp), DIMENSION(2,2)                    :: sor_rho_cmp
        TYPE (Orbit), DIMENSION(:), POINTER         :: orb_arr_cmp
        REAL(bp), DIMENSION(:), POINTER             :: pdf_arr_cmp
        REAL(bp), DIMENSION(:), POINTER             :: rchi2_arr_cmp
        REAL(bp), DIMENSION(:), POINTER             :: reg_apr_arr_cmp
        REAL(bp), DIMENSION(:,:), POINTER           :: jac_arr_cmp
        integer                                     :: numOutOrbits
        REAL(bp), DIMENSION(6)                      :: elements
        CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
        TYPE (Time)                                 :: t0
        real(bp)                                    :: orbMjd
        integer                                     :: j = 0
        integer                                     :: err = 0
        integer                                     :: element_type_index = -1
        
        
        ! Init
        error = .false.
        errorCode = 0
        if(present(element_type)) then
            elementType = element_type
        end if
        ! Figure out the element_type_index
        do j=1,size(ORBITAL_ELEMENTS)
            if(elementType == ORBITAL_ELEMENTS(j)) then
                element_type_index = j
                exit
            end if
        end do
        if(element_type_index .le. 0) then
            ! Error: unsupported element type!
            errorCode = 58
            return
        end if
        
        CALL getResults(storb, sor_rho_cmp=sor_rho_cmp)
        IF (error) THEN
           ! write(*, *) "Error in getResults()"
           errorCode = 15
           return
        end if
        ! Get ORBITAL-ELEMENT PDF
        orb_arr_cmp => getSampleOrbits(storb)
        IF (error) THEN
           ! write(*, *) "Error in getSampleOrbits()"
           errorCode = 16
           return
        END IF
        pdf_arr_cmp => getPDFValues(storb)
        IF (error) THEN
           ! write(*, *) "Error in getPDFValues()"
           errorCode = 17
           return
        END IF
        rchi2_arr_cmp => getReducedChi2Distribution(storb)
        IF (error) THEN
           ! write(*, *) "Error in getReducedChi2Distribution()"
           errorCode = 17
           return
        END IF
        CALL getResults(storb, &
             reg_apr_arr=reg_apr_arr_cmp, &
             jac_arr=jac_arr_cmp)
        IF (error) THEN
           ! write(*, *) "Error in getResults()"
           numOutOrbits = 0
           errorCode = 18
           return
        END IF
    
        ! This is where we wrere originally writing the .orb file(s) out.
        numOutOrbits = SIZE(orb_arr_cmp,dim=1)
!        if(.not. allocated(outOrbits)) then
!            allocate(outOrbits(numOutOrbits, 14), stat=err)
!        end if
!        if(err /= 0) then
!            ! Error in allocating memory.
!            errorCode = 54
!            return
!        end if
        DO j=1,numOutOrbits
           elements = getElements(orb_arr_cmp(j), &
                                  elementType, &
                                  frame="ecliptic")
           IF (error) THEN
              ! write(*, *) "Error in fetching the orbit"
              errorCode = 19
              return
           END IF
           IF (elementType == "keplerian") THEN
              elements(3:6) = elements(3:6)/rad_deg
           ELSE IF (elementType == "delaunay") THEN
              elements(1:3) = elements(1:3)/rad_deg
           ELSE IF (elementType == "poincare") THEN
              elements(4) = elements(4)/rad_deg
           ELSE IF (elementType == "equinoctial") THEN
              elements(6) = elements(6)/rad_deg
           END IF
           
           t0 = getTime(orb_arr_cmp(j))
           orbMjd = getMjd(t0, internal_timescale)
           call nullify(t0)
           IF (error) THEN
              ! write(*, *) "Error in getTime/getMjd"
              errorCode = 5
              return
           END IF
           
           ! Add the values to outOrbits.
           ! outOrbits has the form (elements(1:6), epoch, Un-normalized p.d.f.,
           ! Reduced chi2, Regularized apr, Jacobian det(1:3))
           outOrbits(j, 1) = track_id
           outOrbits(j, 2:7) = elements(1:6)
           outOrbits(j, 8) = orbMjd
           outOrbits(j, 9) = pdf_arr_cmp(j)
           outOrbits(j, 10) = rchi2_arr_cmp(j)
           outOrbits(j, 11) = reg_apr_arr_cmp(j)
           outOrbits(j, 12:14) = jac_arr_cmp(j, 1:3)
           outOrbits(j, 15) = element_type_index
           
!           write(*, *) "elements:             ", elements(1:6)
!           write(*, *) "epoch:                ", orbMjd
!           write(*, *) "Un-normalized p.d.f.: ", pdf_arr_cmp(j)
!           write(*, *) "Reduced chi2:         ", rchi2_arr_cmp(j)
!           write(*, *) "Regularized apr:      ", reg_apr_arr_cmp(j)
!           write(*, *) "Jacobian det:         ", jac_arr_cmp(j, 1:3)
        END DO
        
        ! Cleanup everything
        DO j=1,numOutOrbits
           CALL NULLIFY(orb_arr_cmp(j))
        END DO
        DEALLOCATE(orb_arr_cmp, stat=err)
        DEALLOCATE(pdf_arr_cmp, stat=err)
        DEALLOCATE(rchi2_arr_cmp, stat=err)
        DEALLOCATE(reg_apr_arr_cmp, stat=err)
        DEALLOCATE(jac_arr_cmp, stat=err)
        IF (err /= 0) THEN
           ! write(*, *) "Error: Could not deallocate memory"
           errorCode = 19
           return
        END IF
        return
    end subroutine exportRangingOrbits


    subroutine exportLslOrbit(storb, element_type, track_id, outOrbit,      &
                              outCovariance, outSigmas, outCorrelation,     &
                              errorCode)
        TYPE (StochasticOrbit), intent(in)          :: storb
        CHARACTER(len=*), optional, intent(in)      :: element_type
        integer, intent(in)                         :: track_id
        ! outOrbit has the form (track_id, elements(1:6), epoch, H, G)
        real(bp), dimension(11), intent(out)        :: outOrbit
        REAL(bp), DIMENSION(6,6), intent(out)       :: outCovariance
        REAL(bp), DIMENSION(6)                      :: outSigmas
        REAL(bp), DIMENSION(6,6)                    :: outCorrelation
        integer, intent(out)                        :: errorCode
        
        ! Internal variables.
        REAL(bp), DIMENSION(6)                      :: elements
        CHARACTER(len=ELEMENT_TYPE_LEN)             :: elementType = "keplerian"
        TYPE (Time)                                 :: t0
        real(bp)                                    :: orbMjd
        integer                                     :: j = 0
        integer                                     :: k = 0
        TYPE (Orbit)                                :: orb
        integer                                     :: element_type_index = -1
        
        
        ! Init
        errorCode = 0
        error = .false.
        if(present(element_type)) then
            elementType = element_type
        end if
        ! Figure out the element_type_index
        do j=1,size(ORBITAL_ELEMENTS)
            if(elementType == ORBITAL_ELEMENTS(j)) then
                element_type_index = j
                exit
            end if
        end do
        if(element_type_index .le. 0) then
            ! Error: unsupported element type!
            errorCode = 58
            return
        end if
        
        ! Get the nominal orbit out.
        orb = getNominalOrbit(storb)
        IF (error) THEN
           ! write(*, *) "Error in getNominalOrbit()"
           errorCode = 28
           return
        end if
        
        ! Get the orbital elements.
        elements = getElements(orb, elementType, frame="ecliptic")
        IF (error) THEN
           ! write(*, *) "Error in fetching the orbit"
           errorCode = 19
           return
        END IF
        IF (elementType == "keplerian") THEN
           elements(3:6) = elements(3:6)/rad_deg
        ELSE IF (elementType == "delaunay") THEN
           elements(1:3) = elements(1:3)/rad_deg
        ELSE IF (elementType == "poincare") THEN
           elements(4) = elements(4)/rad_deg
        ELSE IF (elementType == "equinoctial") THEN
           elements(6) = elements(6)/rad_deg
        END IF
        
        ! Get orbit epoch.
        t0 = getTime(orb)
        orbMjd = getMjd(t0, internal_timescale)
        call nullify(t0)
        IF (error) THEN
           ! write(*, *) "Error in getTime/getMjd"
           errorCode = 5
           return
        END IF
        
        ! FIXME: compute H and G.
        outOrbit(1) = track_id
        outOrbit(2:7) = elements(1:6)
        outOrbit(8) = orbMjd
        outOrbit(9) = 99.0_bp                               ! FIXME: compute H
        outOrbit(10) = 9.9_bp                               ! FIXME: compute G
        outOrbit(11) = element_type_index
        
        ! Covariance matrix.
        outCovariance = getCovarianceMatrix(storb, elementType, "ecliptic")
        IF (error) THEN
           ! write(*, *) "Error in getCovarianceMatrix"
           errorCode = 29
           return
        END IF
        DO j=1,6
           outSigmas(j) = SQRT(outCovariance(j, j))
           IF (elementType == "keplerian" .AND. j >= 3) THEN
              outSigmas(j) = outSigmas(j) / rad_deg
           END IF
        END DO
        
        ! Correlations.
        DO j=1,6
           DO k=1,6
              outCorrelation(j,k) = outCovariance(j,k) /                &
                                    (outSigmas(j) * outSigmas(k))
           END DO
        END DO
        
!        write(*, *) "elements:             ", elements(1:6)
!        write(*, *) "epoch:                ", orbMjd
!        write(*, *) "Covariance:           ", outCovariance
!        write(*, *) "Sigmas:               ", outSigmas
!        write(*, *) "Correlation:          ", outCorrelation
    
        ! Cleanup everything
        CALL NULLIFY(orb)
        return
    end subroutine exportLslOrbit


    subroutine exportEphemeris(storb, ephemerides, cov, observers, &
                               outEphems, errorCode)
        ! Input/Output variables.
        type(StochasticOrbit), intent(inout)                :: storb
        ! Output ephem array (intent(inout)).
        TYPE (SphericalCoordinates), DIMENSION(:,:),POINTER :: ephemerides
        ! Ephem uncertainty matrices (intent(inout)).
        REAL(bp), DIMENSION(:,:,:), POINTER                 :: cov
        ! observatory coordinates and ephem time array (intent(inout)).
        TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: observers
        ! Output ephems:
        !   ((dist, ra, dec, mag, mjd, raErr, decErr, smia, smaa, pa), )
        real(bp), dimension(SIZE(observers),10),intent(out) :: outEphems
        ! Output error code
        integer, intent(out)                                :: errorCode
        
        ! Internal vars.
        TYPE (Time)                                         :: t
        integer                                             :: j = 0
        integer                                             :: k = 0
        integer                                             :: l = 0
        real(bp)                                            :: mjd = 0.0_bp
        REAL(bp), DIMENSION(6)                              :: coordinates
        REAL(bp), DIMENSION(6)                              :: stdev
        REAL(bp), DIMENSION(6,6)                            :: corr
        
        
        ! Init
        error = .false.
        errorCode = 0
        
        
        DO j=1,SIZE(observers)
            t = getTime(observers(j))
            mjd = getMJD(t, internal_timescale)
            CALL NULLIFY(t)
            IF (containsSampledPDF(storb)) THEN
               ! We do not support exporting ephems from ranging orbits yet.
               ! FIXME: support ephems from ranging orbits?
               errorCode = 59
               return
            end if

            ! Input orbits correspond to one or more single-point estimates of the pdf.
            ! Make sure that the ephemeris is equatorial:
            CALL rotateToEquatorial(ephemerides(1,j))
            coordinates = getCoordinates(ephemerides(1,j))
            DO k=1,6
               stdev(k) = SQRT(cov(k,k,j)) 
            END DO
            DO k=1,6
               DO l=1,6
                  corr(k,l) = cov(k,l,j) / &
                       (stdev(k)*stdev(l))
               END DO
            END DO
            
            ! Write the output ephem array.
            outEphems(j, 1) = coordinates(1)                ! distance
            outEphems(j, 2:3) = coordinates(2:3)/rad_deg    ! ra/dec
            ! FIXME: compute predicted magnitude!
            outEphems(j, 4) = 99.0_bp                       ! mag
            outEphems(j, 5) = mjd                           ! ephem mjd
            ! FIXME: compute positional uncertainties.
            outEphems(j, 6) = 99.0_bp                       ! raErr
            outEphems(j, 7) = 99.0_bp                       ! decErr
            outEphems(j, 8) = 99.0_bp                       ! semi-major axis
            outEphems(j, 9) = 99.0_bp                       ! semi-minor axis
            outEphems(j, 10) = 99.0_bp                      ! position angle
      END DO
    end subroutine exportEphemeris



    subroutine calendarDateToMjd(year, month, day, mjd, timescale, errorCode)
        integer, intent(in)                             :: year
        integer, intent(in)                             :: month
        ! day + time in fractional days.
        real(bp), intent(in)                            :: day
        CHARACTER(len=*), intent(in)                    :: timescale
        real(bp), intent(out)                           :: mjd
        integer, intent(out)                            :: errorCode
        
        ! Internal vars.
        type(Time)                                      :: t
        
        
        ! Init
        error = .false.
        errorCode = 0
        
        ! Create a Time instance.
        call new(t, year, month, day, trim(timescale))
        if(error) then
            ! Error in creating a Time object.
            errorCode = 55
            return
        end if
        
        ! Convert to MJD TAI
        mjd = getMjd(t, internal_timescale)
        if(error) then
            ! Error in converting to MJD.
            errorCode = 56
            return
        end if
        
        ! Cleanup.
        call nullify(t)
        return
    end subroutine calendarDateToMjd
    
    
    subroutine mjdConvert(mjdIn, timescaleIn, mjdOut, timescaleOut, errorCode)
        ! Input variables
        real(bp), intent(in)                            :: mjdIn
        CHARACTER(len=*), intent(in)                    :: timescaleIn
        real(bp), intent(out)                           :: mjdOut
        CHARACTER(len=*), intent(in)                    :: timescaleOut
        integer, intent(out)                            :: errorCode
        
        ! Internal vars.
        type(Time)                                      :: t
        
        
        ! Init
        error = .false.
        errorCode = 0
        mjdOut = -1.0_bp
        
        ! Create a Time instance.
        call new(t, mjdIn, trim(timescaleIn))
        if(error) then
            ! Error in creating a Time object.
            errorCode = 55
            return
        end if
        
        ! Convert to the desired timescale.
        mjdOut = getMjd(t, trim(timescaleOut))
        if(error) then
            ! Error in converting to MJD.
            errorCode = 56
            return
        end if
        
        ! Cleanup.
        call nullify(t)
        return
    end subroutine mjdConvert
end module oorb

































