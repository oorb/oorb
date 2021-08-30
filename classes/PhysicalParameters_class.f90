!====================================================================!
!                                                                    !
! Copyright 2002-2015,2016                                           !
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
!! Type and routines for various physical parameters such as absolute
!! magnitude H and slope parameter G and the new H,G1,G2 and H,G12
!! systems.
!!
!! @see StochasticOrbit_class 
!!
!! @author  MG
!! @version 2016-03-09
!!
MODULE PhysicalParameters_cl

  USE Base_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  USE CartesianCoordinates_cl
  USE Observations_cl
  USE StochasticOrbit_cl

  USE utilities
  USE planetary_data
  USE linal
  USE sort
  USE statistics

  IMPLICIT NONE

  PRIVATE :: new_PP
  PRIVATE :: new_PP_storb
  PRIVATE :: nullify_PP
  PRIVATE :: copy_PP
  PRIVATE :: exist_PP

  TYPE PhysicalParameters

     PRIVATE
     REAL(bp), DIMENSION(:,:), POINTER :: H0_arr => NULL()
     REAL(bp)                          :: H0_nominal
     REAL(bp)                          :: H0_unc
     REAL(bp), DIMENSION(:,:), POINTER :: G_arr => NULL()
     REAL(bp)                          :: G_nominal
     REAL(bp)                          :: G_unc
     REAL(bp), DIMENSION(:,:), POINTER :: m_arr => NULL()
     REAL(bp)                          :: m_nominal
     REAL(bp)                          :: m_unc     
     TYPE (StochasticOrbit)            :: storb
     LOGICAL                           :: is_initialized = .FALSE.

  END TYPE PhysicalParameters

  INTERFACE NEW
     MODULE PROCEDURE new_PP
     MODULE PROCEDURE new_PP_storb
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_PP
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_PP
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_PP
  END INTERFACE exist

CONTAINS


  !! *Description*:
  !!
  !! Initializes a PhysicalParameters object. 
  !!
  !! Returns error.
  !!
  SUBROUTINE new_PP(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(inout) :: this

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%H0_nominal = 99.9_bp
    this%G_nominal = -9.9_bp
    CALL NEW(this%storb)
    this%is_initialized  = .TRUE.

  END SUBROUTINE new_PP





  !! *Description*:
  !!
  !! Initializes a PhysicalParameters object. 
  !!
  !! Returns error.
  !!
  SUBROUTINE new_PP_storb(this, storb)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(inout) :: this
    TYPE (StochasticOrbit), INTENT(in) :: storb

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%H0_nominal = 99.9_bp
    this%G_nominal = -9.9_bp
    this%storb = copy(storb)
    this%is_initialized  = .TRUE.

  END SUBROUTINE new_PP_storb





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_PP(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(inout) :: this

    INTEGER :: err

    CALL NULLIFY(this%storb)
    this%H0_nominal = 99.9_bp
    this%H0_unc = 99.0_bp
    this%G_nominal = -9.9_bp
    IF (ASSOCIATED(this%H0_arr)) THEN
       DEALLOCATE(this%H0_arr, stat=err)
    END IF
    IF (ASSOCIATED(this%G_arr)) THEN
       DEALLOCATE(this%G_arr, stat=err)
    END IF
    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_PP





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  !! Returns error.
  !!
  FUNCTION copy_PP(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(in) :: this
    TYPE (PhysicalParameters)             :: copy_PP

    copy_PP%H0_nominal     = this%H0_nominal
    copy_PP%G_nominal      = this%G_nominal
    copy_PP%is_initialized = this%is_initialized

  END FUNCTION copy_PP





  !! *Description*:
  !!
  !! Returns the status of the object, i.e. whether it exists or not.
  !!
  LOGICAL FUNCTION exist_PP(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(in) :: this

    exist_PP = this%is_initialized

  END FUNCTION exist_PP



  !! *Description*:
  !!
  !! Compute blackbody spectrum in the wavelength domain. (Ref:
  !! Tähtitieteen perusteet (eng. Fundamentals of Astronomy), pp. 159)
  !!
  !! Input: 
  !!
  !!  T          :: temperature [K] 
  !!  lambda_min :: lower limit for wavelength [m]
  !!  lambda_max :: upper limit for wavelength [m]
  !!
  !! Output:
  !!
  !!  B          :: blackbody spectrum [W m^-2 m^-1 sterad^-1]
  !!
  !! Returns error.
  !!
  SUBROUTINE blackbodyFluxWavelength(T, lambda_min, lambda_max, B)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: T, lambda_min, lambda_max
    REAL(bp), DIMENSION(:), INTENT(out) :: B

    REAL(bp) :: c, c1, c2, c3, dlambda, lambda
    INTEGER :: i

    ! speed of light [m/s]
    c = sol * m_au / sec_day
    ! 2hc^2 [J m^2 s^-1]
    c1 = 2.0_bp * hPl * c**2
    ! hc [J m] 
    c2 = hPl * c
    ! kT [J]
    c3 = kB * T

    dlambda = (lambda_max - lambda_min)/SIZE(B)
    DO i=1,SIZE(B)
       lambda = lambda_min + (i - 0.5_bp) * dlambda
       B(i) = (c1/lambda**5) * 1.0_bp/(EXP(c2/(lambda*c3)) - 1.0_bp)
    END DO

  END SUBROUTINE blackbodyFluxWavelength




  !! *Description*:
  !!
  !! Compute blackbody spectrum in the frequency domain. (Ref:
  !! Tähtitieteen perusteet (eng. Fundamentals of Astronomy), pp. 159)
  !!
  !! Input: 
  !!
  !!  T          :: temperature [K] 
  !!  nu_min     :: lower limit for frequency [Hz]
  !!  nu_max     :: upper limit for frequency [Hz]
  !!
  !! Output:
  !!
  !!  B          :: blackbody spectrum [W m^-2 Hz^-1 sterad^-1]
  !!
  !! Returns error.
  !!
  SUBROUTINE blackbodyFluxFrequency(T, nu_min, nu_max, B)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: T, nu_min, nu_max
    REAL(bp), DIMENSION(:), INTENT(out) :: B

    REAL(bp) :: c, c1, c2, dnu, nu
    INTEGER :: i

    ! speed of light [m/s]
    c = sol * m_au / sec_day
    ! 2hc^-2 [J m^2 s^-1]
    c1 = 2.0_bp * hPl / c**2
    ! h/(kB T) [J m] 
    c2 = hPl / (kB * T)

    dnu = (nu_max - nu_min)/SIZE(B)
    DO i=1,SIZE(B)
       nu = nu_min + (i - 0.5_bp) * dnu
       B(i) = (c1 * nu**3) * 1.0_bp/(EXP(c2 * nu) - 1.0_bp)
    END DO

  END SUBROUTINE blackbodyFluxFrequency





  !! *Description*:
  !!
  !! Compute near-Earth-asteroid thermal model (NEATM) flux in the
  !! wavelength domain. (Ref: Michael Mommert's PhD thesis 2013)
  !!
  !! Input: 
  !!
  !!  T           :: temperature [K] 
  !!  lambda_min  :: lower limit for wavelength [m]
  !!  lambda_max  :: upper limit for wavelength [m]
  !!  phase_angle :: solar phase angle [rad]
  !!
  !! Output:
  !!
  !!  B           :: blackbody spectrum [W m^-2 m^-1 sterad^-1]
  !!
  !! Returns error.
  !!
  SUBROUTINE NEATMFluxWavelength(Tss, lambda_min, lambda_max, phase_angle, B)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: Tss, lambda_min, lambda_max, phase_angle
    REAL(bp), DIMENSION(:), INTENT(out) :: B

    REAL(bp) :: c, dlambda, lambda, theta, dtheta, phi, dphi, temp

    INTEGER :: i, j, k

    ! speed of light [m/s]
    c = sol * m_au / sec_day

    dtheta = (pi/2.0_bp)/10.0_bp
    dphi = pi/20.0_bp
    dlambda = (lambda_max - lambda_min)/SIZE(B)
    DO k=1,SIZE(B)
       lambda = lambda_min + (i - 0.5_bp) * dlambda
       B(k) = 0.0_bp
       DO i=1,10
          theta = dtheta * (i-0.5)
          DO j=1,20
             phi = -(pi/2.0_bp) + phase_angle + dphi * (j-0.5)
             IF (COS(phi)*COS(theta) >= 0.0_bp) THEN
                temp = Tss * COS(phi)**0.25_bp * COS(theta)**0.25_bp
                B(k) = B(k) + &
                     COS(theta)**2 * COS(phi - phase_angle) / &
                     (EXP((hPl*c)/(lambda*kB*temp)) - 1.0_bp) * &
                     dtheta * dphi
             END IF
             !write(*,*) theta/rad_deg, phi/rad_deg, flux
          END DO
       END DO
       ! Absolute flux at 1 au for 1-m body
       B(:) = (hPl * c**2 / lambda**5) * B(:)
    END DO

  END SUBROUTINE NEATMFluxWavelength




  !! *Description*:
  !!
  !! Compute fast-rotating [thermal] model (FRM) flux in the
  !! wavelength domain. (Ref: Michael Mommert's PhD thesis 2013)
  !!
  !! Input: 
  !!
  !!  T          :: temperature [K] 
  !!  lambda_min :: lower limit for wavelength [m]
  !!  lambda_max :: upper limit for wavelength [m]
  !!
  !! Output:
  !!
  !!  B          :: blackbody spectrum [W m^-2 m^-1 sterad^-1]
  !!
  !! Returns error.
  !!
  SUBROUTINE FRMFluxWavelength(Tss, lambda_min, lambda_max, B)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: Tss, lambda_min, lambda_max
    REAL(bp), DIMENSION(:), INTENT(out) :: B

    REAL(bp) :: c, dlambda, lambda, theta, dtheta, temp

    INTEGER :: i, k

    ! speed of light [m/s]
    c = sol * m_au / sec_day

    dtheta = (pi/2.0_bp)/10.0_bp
    dlambda = (lambda_max - lambda_min)/SIZE(B)
    DO k=1,SIZE(B)
       lambda = lambda_min + (i - 0.5_bp) * dlambda
       B(k) = 0.0_bp
       DO i=1,10
          theta = dtheta * (i-0.5)
          temp = Tss * COS(theta)**0.25_bp
          B(k) = B(k) + &
               COS(theta)**2 * SIN(theta) / &
               (EXP((hPl*c)/(lambda*kB*temp)) - 1.0_bp) * &
               dtheta
       END DO
       ! Absolute flux at 1 au for 1-m body
       B(:) = (hPl * c**2 / lambda**5) * B(:)
    END DO

  END SUBROUTINE FRMFluxWavelength






  !! *Description*:
  !!
  !! Compute crude estimate for H by calculating the exact value for
  !! each orbit and observation, and averaging over observations. That
  !! is, for single point estimates (e.g., from LSL) we calculate nobs
  !! H values and assume their average represents the overall H. For a
  !! sampled solution each orbit gets an H value which has been
  !! averaged over observations.
  !!
  !! NB: different bands (AKA filters) are not yet properly treated!
  !!
  !! Returns error.
  !!
  SUBROUTINE crudeHEstimate(this, obss, input_G)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(inout) :: this    
    TYPE (Observations), INTENT(in) :: obss
    REAL(bp), INTENT(in) :: input_G

    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoord_arr => NULL() 
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: &
         ephemerides_arr => NULL()
    TYPE (Orbit), DIMENSION(:,:), POINTER :: &
         orb_lt_corr_arr => NULL()
    TYPE (Orbit), DIMENSION(:), POINTER :: &
         orb_arr => NULL()
    CHARACTER(len=2), DIMENSION(:), POINTER :: &
         filter_arr => NULL()
    CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: &
         filter_arr_
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         cov_arr => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: &
         pdfs_arr => NULL(), &
         phase_angle_arr => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         helio_dist_arr, &
         topo_dist_arr
    REAL(bp), DIMENSION(:), POINTER :: &
         obs_mag_arr => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         H_arr, &
         obs_mag_arr_
    REAL(bp), DIMENSION(6) :: &
         coordinates
    INTEGER :: &
         err, i, j, nobs, norb, nmag

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(this%storb)) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "Orbital information not available.", 1)
       RETURN
    END IF

    IF (.NOT.exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "Photometric information not available.", 1)
       RETURN
    END IF

    ! Extract observational information
    obsy_ccoord_arr => getObservatoryCCoords(obss)
    nobs = SIZE(obsy_ccoord_arr)
    obs_mag_arr => getMagnitudes(obss)
    filter_arr => getFilters(obss)
    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Photometry extracted..."
    END IF

    ! Extract information of observations reporting brightness
    DO i=1,nobs
       IF (obs_mag_arr(i) > 90.0_bp) THEN
          CALL NULLIFY(obsy_ccoord_arr(i))
       END IF
    END DO
    obsy_ccoord_arr => reallocate(obsy_ccoord_arr)
    nmag = SIZE(obsy_ccoord_arr)
    ALLOCATE(obs_mag_arr_(nmag), filter_arr_(nmag), stat=err)
    IF (err /= 0) THEN
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "Could not allocate memory (5).", 1)
       RETURN
    END IF
    j = 0
    DO i=1,nobs
       IF (obs_mag_arr(i) < 90.0_bp) THEN
          j = j + 1
          obs_mag_arr_(j) = obs_mag_arr(i)
          filter_arr_(j) = filter_arr(i)
       END IF
    END DO
    DEALLOCATE(obs_mag_arr, filter_arr)
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2(I0,A))") nmag, " observations out of ", nobs, &
            " contain brightness information."
    END IF

    ! Compute heliocentric and topocentric distances, and phase angles
    IF (containsDiscretePDF(this%storb)) THEN
       orb_arr => getSampleOrbits(this%storb)
    ELSE
       ALLOCATE(orb_arr(1), stat=err)
       IF (err /= 0) THEN
          CALL errorMessage("PhysicalParameters / crudeHEstimate", &
               "Could not allocate memory (10).", 1)
          RETURN
       END IF
       orb_arr(1) = getNominalOrbit(this%storb)
    END IF
    IF (error) THEN
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "TRACE BACK (10).", 1)
       RETURN
    END IF
    CALL getEphemerides(orb_arr, obsy_ccoord_arr, &
         ephemerides_arr, this_lt_corr_arr=orb_lt_corr_arr)
    IF (error) THEN
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "TRACE BACK (15).", 1)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Ephemerides computed..."
    END IF
    norb = SIZE(ephemerides_arr,dim=1)
    ALLOCATE(helio_dist_arr(norb,nmag), topo_dist_arr(norb,nmag), &
         H_arr(nmag), stat=err)
    IF (err /= 0) THEN
       CALL errorMessage("PhysicalParameters / crudeHEstimate", &
            "Could not allocate memory (15).", 1)
       RETURN
    END IF
    DO i=1,nmag
       DO j=1,norb
          coordinates = getElements(orb_lt_corr_arr(j,i), "cartesian", "ecliptic")
          helio_dist_arr(j,i) = SQRT(SUM(coordinates(1:3)**2.0_bp))
          coordinates = getCoordinates(ephemerides_arr(j,i))
          topo_dist_arr(j,i) = coordinates(1)
       END DO
       IF (info_verb >= 3) THEN
          WRITE(stdout,*) coordinates(2:3)/rad_deg
       END IF
    END DO
    CALL getPhaseAngles(this%storb, obsy_ccoord_arr, &
         phase_angle_arr)
    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Phase angles and distances computed..."
    END IF

    ! Estimate H
    IF (containsDiscretePDF(this%storb)) THEN
       ALLOCATE(this%H0_arr(norb,2), this%G_arr(norb,2))
       this%G_arr(:,1) = input_G
       this%G_arr(:,2) = 0.0_bp
       DO i=1,norb
          DO j=1,nmag
             H_arr(j) = getExactH(obs_mag_arr_(j), helio_dist_arr(i,j), &
                  topo_dist_arr(i,j), phase_angle_arr(i,j), &
                  this%G_arr(i,1))
          END DO
          this%H0_arr(i,1) = SUM(H_arr)/nmag
          this%H0_arr(i,2) = 0.0_bp
          IF (info_verb >= 3) THEN
             WRITE(stdout,*) "H0_arr: ", this%H0_arr(i,:), &
                  "   G_arr: ", this%G_arr(i,:) 
          END IF
       END DO
       IF (info_verb >= 2) THEN
          WRITE(stdout,*) "H0_min: ", MINVAL(this%H0_arr(:,1)-this%H0_arr(:,2)), &
               "  H0_max: ", MAXVAL(this%H0_arr(:,1)+this%H0_arr(:,2))
       END IF
    ELSE
       this%G_nominal = input_G
       this%G_unc = 0.0_bp
       DO j=1,nmag
          H_arr(j) = getExactH(obs_mag_arr_(j), helio_dist_arr(1,j), &
               topo_dist_arr(1,j), phase_angle_arr(1,j), &
               this%G_nominal)
          !WRITE(*,*) j, H_arr(j)
       END DO
       this%H0_nominal = SUM(H_arr)/nmag
       this%H0_unc = 0.0_bp
       IF (info_verb >= 2) THEN
          WRITE(stdout,*) "H0: ", this%H0_nominal, "  dH0: ", this%H0_unc
       END IF
    END IF

    ! Deallocate memory
    DEALLOCATE(H_arr, stat=err)
    DEALLOCATE(obsy_ccoord_arr, stat=err)
    DEALLOCATE(obs_mag_arr_, stat=err)
    DEALLOCATE(filter_arr_, stat=err)
    DEALLOCATE(ephemerides_arr, stat=err)
    DEALLOCATE(orb_lt_corr_arr, stat=err)
    DEALLOCATE(orb_arr, stat=err)
    IF (ASSOCIATED(cov_arr)) THEN
       DEALLOCATE(cov_arr, stat=err)
    END IF
    IF (ASSOCIATED(pdfs_arr)) THEN
       DEALLOCATE(pdfs_arr, stat=err)
    END IF
    DEALLOCATE(phase_angle_arr, stat=err)
    DEALLOCATE(helio_dist_arr, stat=err)
    DEALLOCATE(topo_dist_arr, stat=err)

  END SUBROUTINE crudeHEstimate





  !! *Description*:
  !!
  !! Compute estimates for H and G as specified in Bowell et al. 1989.
  !!
  SUBROUTINE estimateHAndG(this, obss, input_G, input_delta_G)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(inout) :: this    
    TYPE (Observations), INTENT(in) :: obss
    REAL(bp), INTENT(in), OPTIONAL :: input_G, input_delta_G

    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
         obsy_ccoord_arr => NULL() 
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: &
         ephemerides_arr => NULL()  
    TYPE (Orbit), DIMENSION(:,:), POINTER :: &
         orb_lt_corr_arr => NULL() 
    CHARACTER(len=2), DIMENSION(:), POINTER :: &
         filter_arr => NULL()
    CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: &
         filter_arr_
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         cov_arr => NULL()
    REAL(bp), DIMENSION(:,:), POINTER :: &
         pdfs_arr => NULL(), &
         phase_angle_arr => NULL()
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         helio_dist_arr, &
         topo_dist_arr
    REAL(bp), DIMENSION(:), POINTER :: &
         obs_mag_arr => NULL(), &
         obs_mag_unc_arr => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         obs_mag_arr_, &
         obs_mag_unc_arr_
    REAL(bp), DIMENSION(6) :: &
         coordinates
    INTEGER :: &
         err, i, j, nmag, norb, nobs

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / estimateHAndG", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(this%storb)) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / estimateHAndG", &
            "Orbital information not available.", 1)
       RETURN
    END IF

    IF (.NOT.exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / estimateHAndG", &
            "Photometric information not available.", 1)
       RETURN
    END IF

    ! Extract observational information
    obsy_ccoord_arr => getObservatoryCCoords(obss)
    nobs = SIZE(obsy_ccoord_arr)
    obs_mag_arr => getMagnitudes(obss)
    obs_mag_unc_arr => getMagnitudeUncertainties(obss)
    filter_arr => getFilters(obss)
    IF (error) THEN
       CALL errorMessage("PhysicalParameters / estimateHAndG", &
            "TRACE BACK (5).", 1)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Photometry extracted..."
    END IF
    IF (info_verb >= 3) THEN
       WRITE(stdout,*) obs_mag_unc_arr(:)
    END IF

    ! Extract information of observations reporting brightness
    DO i=1,nobs
       IF (obs_mag_arr(i) > 90.0_bp) THEN
          CALL NULLIFY(obsy_ccoord_arr(i))
       END IF
    END DO
    obsy_ccoord_arr => reallocate(obsy_ccoord_arr)
    IF (.NOT.ASSOCIATED(obsy_ccoord_arr)) THEN
       CALL warningMessage("PhysicalParameters / estimateHAndG", &
            "0 observations contain brightness information -> " // &
            "estimation of H and G is impossible.", 1)
       RETURN
    END IF
    nmag = SIZE(obsy_ccoord_arr)
    ALLOCATE(obs_mag_arr_(nmag), filter_arr_(nmag), obs_mag_unc_arr_(nmag))
    j = 0
    DO i=1,nobs
       IF (obs_mag_arr(i) < 90.0_bp) THEN
          j = j + 1
          obs_mag_arr_(j) = obs_mag_arr(i)
          obs_mag_unc_arr_(j) = obs_mag_unc_arr(i)
          filter_arr_(j) = filter_arr(i)
       END IF
    END DO
    DEALLOCATE(obs_mag_arr, obs_mag_unc_arr, filter_arr)
    IF (info_verb >= 2) THEN
       WRITE(stdout,"(2(I0,A))") nmag, " observations out of ", nobs, &
            " contain brightness information."
    END IF

    ! Compute heliocentric and topocentric distances, and phase angles
    IF (containsDiscretePDF(this%storb)) THEN
       CALL getEphemerides(this%storb, obsy_ccoord_arr, &
            ephemerides_arr, pdf_arr=pdfs_arr, this_lt_corr_arr=orb_lt_corr_arr)
    ELSE
       CALL getEphemerides(this%storb, obsy_ccoord_arr, &
            ephemerides_arr, cov_arr=cov_arr, this_lt_corr_arr=orb_lt_corr_arr)       
    END IF
    IF (error) THEN
       CALL errorMessage("PhysicalParameters / estimateHAndG", &
            "TRACE BACK (10).", 1)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Ephemerides computed..."
    END IF
    norb = SIZE(ephemerides_arr,dim=1)
    ALLOCATE(helio_dist_arr(norb,nmag), topo_dist_arr(norb,nmag))
    DO i=1,nmag
       DO j=1,norb
          coordinates = getElements(orb_lt_corr_arr(j,i), "cartesian", "ecliptic")
          helio_dist_arr(j,i) = SQRT(SUM(coordinates(1:3)**2.0_bp))
          coordinates = getCoordinates(ephemerides_arr(j,i))
          topo_dist_arr(j,i) = coordinates(1)
          IF (error) THEN
             CALL errorMessage("PhysicalParameters / estimateHAndG", &
                  "TRACE BACK (15).", 1)
             RETURN
          END IF
       END DO
       IF (info_verb >= 3) THEN
          WRITE(stdout,*) coordinates(2:3)/rad_deg
       END IF
    END DO
    CALL getPhaseAngles(this%storb, obsy_ccoord_arr, &
         phase_angle_arr)
    IF (error) THEN
       CALL errorMessage("PhysicalParameters / estimateHAndG", &
            "TRACE BACK (20).", 1)
       RETURN
    END IF
    IF (info_verb >= 2) THEN
       WRITE(stdout,*) "Phase angles and distances computed..."
    END IF

    ! Estimate H and G
    IF (containsDiscretePDF(this%storb)) THEN
       ALLOCATE(this%H0_arr(norb,2), this%G_arr(norb,2))
       DO i=1,norb
          IF (PRESENT(input_G)) THEN
             CALL HG(obs_mag_arr_, obs_mag_unc_arr_, helio_dist_arr(i,:), &
                  topo_dist_arr(i,:), phase_angle_arr(i,:), &
                  this%H0_arr(i,1), this%H0_arr(i,2), this%G_arr(i,1), &
                  this%G_arr(i,2), input_G, input_delta_G)
          ELSE
             CALL HG(obs_mag_arr_, obs_mag_unc_arr_, helio_dist_arr(i,:), &
                  topo_dist_arr(i,:), phase_angle_arr(i,:), &
                  this%H0_arr(i,1), this%H0_arr(i,2), this%G_arr(i,1), &
                  this%G_arr(i,2))
          END IF
          IF (error) THEN
             CALL errorMessage("PhysicalParameters / estimateHAndG", &
                  "TRACE BACK (25).", 1)
             RETURN
          END IF
          IF (info_verb >= 3) THEN
             WRITE(stdout,*) "H0_arr: ", this%H0_arr(i,:), &
                  "   G_arr: ", this%G_arr(i,:) 
          END IF
       END DO
       IF (info_verb >= 2) THEN
          WRITE(stdout,*) "H0_min: ", MINVAL(this%H0_arr(:,1)-this%H0_arr(:,2)), &
               "  H0_max: ", MAXVAL(this%H0_arr(:,1)+this%H0_arr(:,2))
          WRITE(stdout,*) "G_min:  ", MINVAL(this%G_arr(:,1)-this%G_arr(:,2)), &
               "  G_max:  ", MAXVAL(this%G_arr(:,1)+this%G_arr(:,2))
       END IF
    ELSE
       IF (PRESENT(input_G)) THEN
          CALL HG(obs_mag_arr_, obs_mag_unc_arr_, helio_dist_arr(1,:), &
               topo_dist_arr(1,:), phase_angle_arr(1,:), &
               this%H0_nominal, this%H0_unc, this%G_nominal, &
               this%G_unc, input_G, input_delta_G)
       ELSE
          CALL HG(obs_mag_arr_, obs_mag_unc_arr_, helio_dist_arr(1,:), &
               topo_dist_arr(1,:), phase_angle_arr(1,:), &
               this%H0_nominal, this%H0_unc, this%G_nominal, &
               this%G_unc)
       END IF
       IF (error) THEN
          CALL errorMessage("PhysicalParameters / estimateHAndG", &
               "TRACE BACK (30).", 1)
          RETURN
       END IF
       IF (info_verb >= 2) THEN
          WRITE(stdout,*) "H0: ", this%H0_nominal, "  dH0: ", this%H0_unc
          WRITE(stdout,*) "G:  ", this%G_nominal,  "  dG:  ", this%G_unc
       END IF
    END IF

    ! Deallocate memory
    DEALLOCATE(obsy_ccoord_arr, stat=err)
    DEALLOCATE(obs_mag_arr_, stat=err)
    DEALLOCATE(obs_mag_unc_arr_, stat=err)
    DEALLOCATE(filter_arr_, stat=err)
    DEALLOCATE(ephemerides_arr, stat=err)
    DEALLOCATE(orb_lt_corr_arr, stat=err)
    IF (ASSOCIATED(cov_arr)) THEN
       DEALLOCATE(cov_arr, stat=err)
    END IF
    IF (ASSOCIATED(pdfs_arr)) THEN
       DEALLOCATE(pdfs_arr, stat=err)
    END IF
    DEALLOCATE(phase_angle_arr, stat=err)
    DEALLOCATE(helio_dist_arr, stat=err)
    DEALLOCATE(topo_dist_arr, stat=err)

  END SUBROUTINE estimateHAndG





  !! *Description*:*
  !!
  !! The apparent magnitude from the H,G magnitude system.
  !!
  REAL(bp) FUNCTION getApparentHGMagnitude(H, G, r, Delta, phase_angle)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: H, G, r, Delta, phase_angle

    REAL(bp), DIMENSION(2) :: phi

    phi = HGPhaseFunctions(phase_angle)
    getApparentHGMagnitude = H - 2.5_bp*LOG10((1.0_bp - G)*phi(1) + &
         G*phi(2)) + 5.0_bp*LOG10(r*Delta)

  END FUNCTION getApparentHGMagnitude





  !! *Description*:*
  !!
  !! Returns the apparent magnitude from the H,G1,G2 magnitude system
  !! (see Muinonen et al. 2010, Icarus, 209, 542-555 for details).
  !!
  REAL(bp) FUNCTION getApparentHG1G2Magnitude(H, G1, G2, r, Delta, phase_angle)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: H, G1, G2, r, Delta, phase_angle

    REAL(bp), DIMENSION(3) :: phi

    phi = HG1G2BasisFunctions(phase_angle)
    getApparentHG1G2Magnitude = H - 2.5_bp*LOG10(G1*phi(1) + G2*phi(2) + &
         (1.0_bp - G1 - G2)*phi(3)) + 5.0_bp*LOG10(r*Delta)

  END FUNCTION getApparentHG1G2Magnitude





  !! *Description*:*
  !!
  !! Returns the apparent magnitude from the H,G12 magnitude system
  !! (see Muinonen et al. 2010, Icarus, 209, 542-555 for details).
  !!
  REAL(bp) FUNCTION getApparentHG12Magnitude(H, G12, r, Delta, phase_angle)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: H, G12, r, Delta, phase_angle

    REAL(bp), DIMENSION(3) :: phi
    REAL(bp) :: G1, G2

    IF (G12 < 0.2_bp) THEN
       G1 = 0.7527_bp*G12 + 0.06164_bp
       G2 = -0.9612_bp*G12 + 0.6270_bp
    ELSE
       G1 = 0.9529_bp*G12 + 0.02162_bp
       G2 = -0.6125_bp*G12 + 0.5572_bp
    END IF

    phi = HG1G2BasisFunctions(phase_angle)

    getApparentHG12Magnitude = H - 2.5_bp*LOG10(G1*phi(1) + G2*phi(2) + &
         (1.0_bp - G1 - G2)*phi(3)) + 5.0_bp*LOG10(r*Delta)

  END FUNCTION getApparentHG12Magnitude





  !! *Description*:*
  !!
  !! The apparent cometary magnitude using theory reviewed in Cook et
  !! al. 2012, in prep.
  !!
  REAL(bp) FUNCTION getApparentCometaryMagnitude(H10, G, r, Delta, phase_angle, n)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: H10, G, r, Delta, phase_angle, n

    REAL(bp), DIMENSION(2) :: phi

    phi = HGPhaseFunctions(phase_angle)

    getApparentCometaryMagnitude = &
         H10 - &
         2.5_bp*LOG10((1.0_bp - G)*phi(1) + G*phi(2)) + &
         2.5_bp*(0.5_bp*n*LOG10(r**2) + LOG10(Delta**2))

  END FUNCTION getApparentCometaryMagnitude





  !! *Description*:*
  !!
  !! Returns the slope parameter G and its uncertainty.
  !!
  FUNCTION getG(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(in) :: this    
    REAL(bp), DIMENSION(2) :: getG

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / getG", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getG = (/ this%G_nominal, this%G_unc /)

  END FUNCTION getG





  !! *Description*:*
  !!
  !! Returns a distribution of slope parameters G and their
  !! uncertainties corresponding to an orbital-element distribution.
  !!
  FUNCTION getGDistribution(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(in) :: this    
    REAL(bp), DIMENSION(:,:), POINTER :: getGDistribution

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / getGDistribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ALLOCATE(getGDistribution(SIZE(this%G_arr,dim=1),2))
    getGDistribution = this%G_arr

  END FUNCTION getGDistribution





  !! *Description*:*
  !!
  !! Returns the H(alpha=0) magnitude and its uncertainty.
  !!
  FUNCTION getH0(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(in) :: this    
    REAL(bp), DIMENSION(2) :: getH0

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / getH0", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getH0 = (/ this%H0_nominal, this%H0_unc /)

  END FUNCTION getH0





  !! *Description*:*
  !!
  !! Returns a distribution of H(alpha=0) magnitudes and their
  !! uncertainties corresponding to an orbital-element distribution.
  !!
  FUNCTION getH0Distribution(this)

    IMPLICIT NONE
    TYPE (PhysicalParameters), INTENT(in) :: this    
    REAL(bp), DIMENSION(:,:), POINTER :: getH0Distribution

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / getH0Distribution", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.ASSOCIATED(this%H0_arr)) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / getH0Distribution", &
            "Object does not contain H0 distribution.", 1)
       RETURN
    END IF

    ALLOCATE(getH0Distribution(SIZE(this%H0_arr,dim=1),2))
    getH0Distribution = this%H0_arr

  END FUNCTION getH0Distribution





  !! *Description*:
  !!
  !! Compute estimates for H and G and their uncertainties as
  !! specified in Bowell et al. 1989.
  !!
  SUBROUTINE HG(obs_mag_arr, obs_mag_unc_arr, helio_dist_arr, &
       topo_dist_arr, phase_angle_arr, H0, delta_H0, G, delta_G, input_G, &
       input_delta_G)

    IMPLICIT NONE

    REAL(bp), DIMENSION(:), INTENT(in) :: obs_mag_arr, &
         obs_mag_unc_arr, helio_dist_arr, topo_dist_arr, phase_angle_arr
    REAL(bp), INTENT(out) :: H0, G, delta_H0, delta_G
    REAL(bp), INTENT(in), OPTIONAL :: input_G, input_delta_G

    REAL(bp), DIMENSION(SIZE(obs_mag_arr),2) :: phi_arr
    REAL(bp), DIMENSION(SIZE(obs_mag_arr)) :: reduced_mag_arr, &
         mag_resid_arr, delta_malpha_arr
    REAL(bp), DIMENSION(2,2) :: h
    REAL(bp), DIMENSION(2) :: gg, a
    REAL(bp) :: bigI, D, alpha0, sigma2_Halpha0, sigma2_beta, s2, &
         delta_Halpha0, delta_beta
    INTEGER :: i

    ! Note that geo_dist has been changed to topo_dist:
    reduced_mag_arr = obs_mag_arr - 5.0_bp*LOG10(helio_dist_arr*topo_dist_arr)

    DO i=1,SIZE(reduced_mag_arr)
       phi_arr(i,1:2) = HGPhaseFunctions(phase_angle_arr(i))
       !phi_arr(i,1:2) = simplifiedHGPhaseFunctions(phase_angle_arr(i))    
    END DO

    IF (ANY(obs_mag_unc_arr == 0.0_bp)) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / HG", &
            "One or more reported brightness uncertainties are equal to zero.", 1)
       RETURN       
    END IF

    IF (SIZE(obs_mag_arr) < 3) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / HG", &
            "Needs 3 or more brightness measurements to estimate uncertainty on H and/or G.", 1)
       RETURN       
    END IF

    IF (PRESENT(input_G)) THEN

       ! Compute nominal LS value for H(alpha=0) assuming fixed G

       G = input_G       
       delta_malpha_arr = -2.5_bp * LOG10((1.0_bp - G)*phi_arr(:,1) + G*phi_arr(:,2))
       H0 = SUM((reduced_mag_arr - delta_malpha_arr)/obs_mag_unc_arr**2.0_bp) / &
            SUM(1.0_bp/obs_mag_unc_arr**2.0_bp)

    ELSE

       ! Compute nominal LS values for both H(alpha=0) and G 

       h = 0.0_bp
       gg = 0.0_bp
       DO i=1,SIZE(reduced_mag_arr)
          bigI = 10**(-0.4_bp*reduced_mag_arr(i))
          h(1,1) = h(1,1) + phi_arr(i,1)**2.0_bp/(obs_mag_unc_arr(i)**2.0_bp*bigI**2.0_bp)
          h(2,2) = h(2,2) + phi_arr(i,2)**2.0_bp/(obs_mag_unc_arr(i)**2.0_bp*bigI**2.0_bp)
          h(1,2) = h(1,2) + phi_arr(i,1)*phi_arr(i,2)/(obs_mag_unc_arr(i)**2.0_bp*bigI**2.0_bp)
          gg = gg + phi_arr(i,1:2)/(obs_mag_unc_arr(i)**2.0_bp*bigI**2.0_bp)
       END DO
       D = h(1,1)*h(2,2) - h(1,2)**2.0_bp
       a(1) = (h(2,2)*gg(1) - h(1,2)*gg(2))/D
       a(2) = (h(1,1)*gg(2) - h(1,2)*gg(1))/D
       H0 = -2.5_bp*LOG10(SUM(a))
       ! Note that the definition of G is wrong in the AII chapter!
       G = a(2)/SUM(a)

    END IF


    IF (PRESENT(input_G) .AND. PRESENT(input_delta_G)) THEN

       ! Compute uncertainty for H(alpha=0) assuming fixed G and delta_G

       mag_resid_arr = reduced_mag_arr - H0 + &
            2.5_bp * LOG10((1.0_bp - G)*phi_arr(:,1) + G*phi_arr(:,2))
       alpha0 = SUM(phase_angle_arr/obs_mag_unc_arr**2.0_bp) / &
            SUM(1.0_bp/obs_mag_unc_arr**2.0_bp)
       sigma2_Halpha0 = 1.0_bp/SUM(1.0_bp/obs_mag_unc_arr**2.0_bp)
       s2 = 1.0_bp/(SIZE(obs_mag_arr) - 2.0_bp)*SUM((mag_resid_arr/obs_mag_unc_arr)**2.0_bp)
       delta_Halpha0 = SQRT(s2*sigma2_Halpha0)
       delta_G = input_delta_G
       delta_beta = delta_G*(0.0673_bp - 0.1132_bp*G + 0.0615_bp*G**2.0_bp)
       delta_H0 = SQRT(s2*sigma2_Halpha0 + (delta_beta*alpha0)**2.0_bp)       

    ELSE IF (.NOT.PRESENT(input_G) .AND. .NOT.PRESENT(input_delta_G)) THEN

       ! Compute uncertainties for H(alpha=0) and G

       mag_resid_arr = reduced_mag_arr - H0 + &
            2.5_bp * LOG10((1.0_bp - G)*phi_arr(:,1) + G*phi_arr(:,2))
       alpha0 = SUM(phase_angle_arr/obs_mag_unc_arr**2.0_bp) / &
            SUM(1.0_bp/obs_mag_unc_arr**2.0_bp)
       sigma2_Halpha0 = 1.0_bp/SUM(1.0_bp/obs_mag_unc_arr**2.0_bp)
       sigma2_beta = 1.0_bp / (SUM(phase_angle_arr**2.0_bp/obs_mag_unc_arr**2.0_bp) - &
            alpha0**2.0_bp*SUM(1.0_bp/obs_mag_unc_arr**2.0_bp))
       s2 = 1.0_bp/(SIZE(obs_mag_arr) - 2.0_bp)*SUM((mag_resid_arr/obs_mag_unc_arr)**2.0_bp)
       delta_Halpha0 = SQRT(s2*sigma2_Halpha0)
       delta_beta = SQRT(s2*sigma2_beta)
       delta_G = delta_beta/(0.0673_bp - 0.1132_bp*G + 0.0615_bp*G**2.0_bp)
       delta_H0 = SQRT(s2*sigma2_Halpha0 + s2*sigma2_beta*alpha0**2.0_bp)

    ELSE

       ! Uncertainties cannot be computed

       delta_H0 = -1.0_bp
       delta_G = -1.0_bp

    END IF

  END SUBROUTINE HG





  !! *Description*:
  !!
  !! Compute the exact value for H corresponding to an observed magnitude,
  !! a heliocentric distance, a topocentric distance, a phase angle, and
  !! a slope parameter.
  !!
  REAL(bp) FUNCTION getExactH(obs_mag, helio_dist, &
       topo_dist, phase_angle, G)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: obs_mag, &
         helio_dist, &
         topo_dist, &
         phase_angle, &
         G

    REAL(bp), DIMENSION(2) :: phi
    REAL(bp) :: reduced_mag, delta_malpha

    ! Note that geo_dist has been changed to topo_dist:
    reduced_mag = obs_mag - 5.0_bp*LOG10(helio_dist*topo_dist)
    phi(1:2) = HGPhaseFunctions(phase_angle)

    ! Compute nominal value for H(alpha=0) assuming fixed G
    delta_malpha = -2.5_bp * LOG10((1.0_bp - G)*phi(1) + G*phi(2))
    getExactH = reduced_mag - delta_malpha

  END FUNCTION getExactH





  FUNCTION HGPhaseFunctions(phase_angle)

    IMPLICIT NONE
    REAL(bp), DIMENSION(2) :: HGPhaseFunctions
    REAL(bp), INTENT(in) :: phase_angle

    REAL(bp), DIMENSION(2) :: a, b, c
    REAL(bp) :: w, phis, phil
    INTEGER :: i

    a(1) = 3.332_bp
    a(2) = 1.862_bp
    b(1) = 0.631_bp
    b(2) = 1.218_bp
    c(1) = 0.986_bp
    c(2) = 0.238_bp
    w = EXP(-90.56_bp*TAN(0.5_bp*phase_angle)**2)

    DO i=1,2
       phis = 1.0_bp - C(i)*SIN(phase_angle) / &
            (0.119_bp+1.341_bp*SIN(phase_angle)-0.754_bp*SIN(phase_angle)**2) 
       phil = EXP(-a(i)*TAN(0.5_bp*phase_angle)**b(i))
       HGPhaseFunctions(i) = w*phis + (1.0_bp-w)*phil
    END DO

  END FUNCTION HGPhaseFunctions





  !! *Description*:
  !!
  !! Returns the values of the basis functions phi1, phi2 and phi3 of
  !! the H,G1,G2 system at the given phase angle (see Muinonen et
  !! al. 2010, Icarus, 209, 542-555 for details).
  !!
  !! Based on implementation by H. Penttilä.
  !!
  FUNCTION HG1G2BasisFunctions(phase_angle)

    IMPLICIT NONE
    REAL(bp), DIMENSION(3) :: HG1G2BasisFunctions
    REAL(bp), INTENT(in) :: phase_angle

    ! The following splines for basis functions phi1, phi2 and phi3
    ! are defined in Muinonen et al. (2010, Icarus 209, 542-555) and
    ! the derivatives are computed using an algorithm based on the
    ! public sw produced by H. Penttilä. The knots12 and knots3 arrays
    ! contain phase angle at knot, value of basis function at knot,
    ! and derivative at knot for each spline representation
    ! (i=function, j=knot, k=(angle,value,derivative)).

    ! phi1 and phi2 
    REAL(bp), DIMENSION(2,6,3), PARAMETER :: knots12 = RESHAPE((/ &
         0.13089969389957470_bp, 0.13089969389957470_bp, &
         0.52359877559829882_bp, 0.52359877559829882_bp, &
         1.0471975511965976_bp, 1.0471975511965976_bp, &
         1.5707963267948966_bp, 1.5707963267948966_bp, &
         2.0943951023931953_bp, 2.0943951023931953_bp, &
         2.6179938779914944_bp, 2.6179938779914944_bp, &
         0.75000000000000000_bp, 0.92500000000000004_bp, &
         0.33486016000000002_bp, 0.62884169000000001_bp, &
         0.13410559999999999_bp, 0.31755495000000000_bp, &
         5.11047560000000012E-002_bp, 0.12716367000000001_bp, &
         2.14656870000000007E-002_bp, 2.23739030000000005E-002_bp, &
         3.63969890000000011E-003_bp, 1.65056890000000008E-004_bp, &
         -1.9098593171027440_bp, -0.57295779513082323_bp, &
         -0.55463432402954305_bp, -0.76705367127692980_bp, &
         -0.24404598605661942_bp, -0.45665789025666331_bp, &
         -9.49804380669390519E-002_bp, -0.28071808974438800_bp, &
         -2.14114236377018034E-002_bp, -0.11173256930106369_bp, &
         -9.13286120000000035E-002_bp, -8.65731380000000014E-008_bp &
         /), (/ 2,6,3 /))
    ! phi3
    REAL(bp), DIMENSION(9,3) :: knots3 = RESHAPE((/ &
         0.00000000000000000E+000_bp, 5.23598775598298812E-003_bp, &
         1.74532925199432955E-002_bp, 3.49065850398865909E-002_bp, &
         6.98131700797731819E-002_bp, 0.13962634015954636_bp, &
         0.20943951023931956_bp, 0.34906585039886590_bp, &
         0.52359877559829882_bp, 1.0000000000000000_bp, &
         0.83381185000000002_bp, 0.57735424000000002_bp, &
         0.42144772000000003_bp, 0.23174230000000001_bp, &
         0.10348178000000000_bp, 6.17334729999999970E-002_bp, &
         1.61070060000000001E-002_bp, 0.00000000000000000E+000_bp, &
         -0.10630096999999999_bp, -41.180439484750558_bp, &
         -10.366914741841031_bp, -7.5784614901018310_bp, &
         -3.6960950182737276_bp, -0.78605651603897575_bp, &
         -0.46527011787268774_bp, -0.20459544750797629_bp, &
         0.00000000000000000E+000_bp /), (/ 9, 3 /))

    REAL(bp) :: phase_angle_lo, phase_angle_hi, &
         basis_function_lo, basis_function_hi, &
         phase_angle_derivative_lo, phase_angle_derivative_hi, &
         Delta_phase_angle, Delta_basis_function, &
         a, b, t
    INTEGER :: i, j

    ! Basis functions 1 and 2
    IF (phase_angle < 7.5_bp*rad_deg) THEN
       DO i=1,2
          HG1G2BasisFunctions(i) = & 
               (phase_angle - knots12(i,1,1))*knots12(i,1,3) + &
               knots12(i,1,2)
       END DO
    ELSE
       DO i=1,2 ! loop over basis functions 1 and 2
          DO j=1,SIZE(knots12,dim=2) - 1 ! loop over phase-angle bins
             IF (knots12(i,j,1) <= phase_angle .AND. &
                  phase_angle <= knots12(i,j+1,1)) EXIT
          END DO
          phase_angle_lo = knots12(i,j,1)
          phase_angle_hi = knots12(i,j+1,1)
          basis_function_lo = knots12(i,j,2)
          basis_function_hi = knots12(i,j+1,2)
          phase_angle_derivative_lo = knots12(i,j,3)
          phase_angle_derivative_hi = knots12(i,j+1,3)
          Delta_phase_angle = phase_angle_hi - phase_angle_lo
          Delta_basis_function = basis_function_hi - basis_function_lo
          a = phase_angle_derivative_lo*Delta_phase_angle - &
               Delta_basis_function
          b = -phase_angle_derivative_hi*Delta_phase_angle + &
               Delta_basis_function
          t = (phase_angle - phase_angle_lo) / Delta_phase_angle
          HG1G2BasisFunctions(i) = (1-t)*basis_function_lo + &
               t*basis_function_hi + t*(1-t)*(a*(1-t)+b*t)
       END DO
    END IF

    ! Basis function 3
    IF (phase_angle > 30.0_bp*rad_deg) THEN
       HG1G2BasisFunctions(3) = 0.0_bp
    ELSE
       DO j=1,SIZE(knots3,dim=1) - 1 ! loop over phase-angle bins
          IF (knots3(j,1) <= phase_angle .AND. &
               phase_angle <= knots3(j+1,1)) EXIT
       END DO
       phase_angle_lo = knots3(j,1)
       phase_angle_hi = knots3(j+1,1)
       basis_function_lo = knots3(j,2)
       basis_function_hi = knots3(j+1,2)
       phase_angle_derivative_lo = knots3(j,3)
       phase_angle_derivative_hi = knots3(j+1,3)
       Delta_phase_angle = phase_angle_hi - phase_angle_lo
       Delta_basis_function = basis_function_hi - basis_function_lo
       a = phase_angle_derivative_lo*Delta_phase_angle - &
            Delta_basis_function
       b = -phase_angle_derivative_hi*Delta_phase_angle + &
            Delta_basis_function
       t = (phase_angle - phase_angle_lo) / Delta_phase_angle
       HG1G2BasisFunctions(3) = (1-t)*basis_function_lo + &
            t*basis_function_hi + t*(1-t)*(a*(1-t)+b*t)
    END IF

  END FUNCTION HG1G2BasisFunctions





  SUBROUTINE spline_coefs(xval, yval, deriv, d)

    IMPLICIT NONE
    REAL(bp), DIMENSION(:), INTENT(IN) :: xval, yval
    REAL(bp), DIMENSION(:), INTENT(out) :: d
    REAL(bp), DIMENSION(2), INTENT(IN) :: deriv

    INTEGER :: i, j, n
    REAL(bp) :: bet
    REAL(bp), DIMENSION(SIZE(xval,1)) :: a, b, c, r, gam, u

    n = SIZE(xval,1)

    ! Spline system matrix A coefficients

    a(1) = 0.0_bp
    DO i=2,n-1
       a(i) = 1.0_bp/(xval(i)-xval(i-1))
    END DO
    a(n) = 0.0_bp

    b(1) = 1.0_bp
    DO i=2,n-1
       b(i) = 2.0_bp/(xval(i)-xval(i-1)) + 2.0_bp/(xval(i+1)-xval(i))
    END DO
    b(n) = 1.0_bp

    c(1) = 0.0_bp
    DO i=2,n-1
       c(i) = 1.0_bp/(xval(i+1)-xval(i))
    END DO
    c(n) = 0.0_bp

    r(1) = deriv(1)
    DO i=2,n-1
       r(i) = 3.0_bp * ((yval(i)-yval(i-1))/((xval(i)-xval(i-1))**2) + &
            (yval(i+1)-yval(i))/((xval(i+1)-xval(i))**2))
    END DO
    r(n) = deriv(2)

    ! Tri-diagonal solver
    bet=b(1)
    u(1)=r(1)/bet
    DO i=2,n
       gam(i)=c(i-1)/bet
       bet=b(i)-a(i)*gam(i)
       u(i)=(r(i)-a(i)*u(i-1))/bet
    END DO

    DO i=n-1,1,-1
       u(i)=u(i)-gam(i+1)*u(i+1) ! derivative
    END DO

    DO i=1,n
       d(i) = u(i) 
       !write(*,*) xval(i), yval(i), u(i)
    END DO

  END SUBROUTINE spline_coefs





  FUNCTION simplifiedHGPhaseFunctions(phase_angle)

    IMPLICIT NONE
    REAL(bp), DIMENSION(2) :: simplifiedHGPhaseFunctions
    REAL(bp), INTENT(in) :: phase_angle
    REAL(bp), DIMENSION(2) :: a, b
    INTEGER :: i

    a(1) = 3.33_bp
    a(2) = 1.87_bp
    b(1) = 0.63_bp
    b(2) = 1.22_bp

    DO i=1,2
       simplifiedHGPhaseFunctions(i) = EXP(-a(i)*TAN(0.5_bp*phase_angle)**b(i))
    END DO

  END FUNCTION simplifiedHGPhaseFunctions




  !! MCMC mass estimation algorithm. See Siltala & Granvik (2020,2017)
  !! https://doi.org/10.1051/0004-6361/201935608
  !! https://doi.org/10.1016/j.icarus.2017.06.028
  !! Author: LS, MG
  SUBROUTINE massEstimation_MCMC(storb_arr, orb_arr, &
       proposal_density_masses, norb, iorb_init, itrial_init, nburn_arr, &
       estimated_masses, accepted_solutions, nominal_arr, &
       adaptation, delayed_rejection,out_fname, input_cov_matrix, mass_lock)

    IMPLICIT NONE

    TYPE (StochasticOrbit), DIMENSION(:), INTENT(inout) :: storb_arr
    TYPE (Orbit), DIMENSION(:), INTENT(out) :: nominal_arr
    TYPE (Orbit), DIMENSION(:), INTENT(inout) :: &
         orb_arr
    REAL(bp), DIMENSION(:), INTENT(in) :: proposal_density_masses
    REAL(bp) :: nsigma_mass, nsigma_orbit
    INTEGER, INTENT(in) :: norb
    INTEGER, INTENT(inout) :: iorb_init, itrial_init
    LOGICAL, INTENT(inout) :: mass_lock, delayed_rejection
    INTEGER, DIMENSION(2), OPTIONAL :: nburn_arr
    REAL(bp), DIMENSION(:,:), INTENT(out) :: estimated_masses
    REAL(bp), DIMENSION(:,:), INTENT(inout) :: accepted_solutions
    REAL(bp), DIMENSION(:,:), POINTER :: input_cov_matrix
    CHARACTER(len=256), INTENT(inout) :: &
         adaptation
    CHARACTER(len=FNAME_LEN), INTENT(inout ):: &
         out_fname

    TYPE (SparseArray) :: resids
    TYPE (Observations) :: obss
    TYPE (Observations), DIMENSION(:), ALLOCATABLE :: obs_arr
    TYPE (Orbit), DIMENSION(:,:), ALLOCATABLE :: orb_arr2
    TYPE (Orbit) :: &
         orb
    TYPE (Time) :: &
         epoch
    TYPE (File) :: mcmc_out_file, meanres_file, cov_file, res_file
    CHARACTER(len=DYN_MODEL_LEN) :: &
         dyn_model                                                 !! Dynamical model.
    CHARACTER(len=INTEGRATOR_LEN) :: &
         integrator                                                !! Integrator.
    CHARACTER(len=12) :: &
         str1, &
         str2
    CHARACTER(len=64) :: &
         efrmt = "(E11.4,1X)"
    REAL(bp), DIMENSION(:,:,:), POINTER :: information_matrix => NULL()

    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE :: &
         A_arr, covariance2
    REAL(bp), DIMENSION(:,:), POINTER :: &
         res
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         cov_matrix, A_arr2, ok_cov_matrix, previous_cov_matrix ,&
         additional_perturbers, &
         additional_perturbers_, &
         elements0_arr, &
         mean_resids, &
         elem_temp,&
         ya_temp, &
         deviates_matrix,&
         lambda_arr,&
         elements_arr!, &
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         chi2_arr, &
         chi2_arr_updated, &
         last_accepted_chi2_arr, &
         mean_arr, &
         old_mean, &
         another_temp,&
         ran_arr, &
         deviates, &
         last_proposal, &
         rchi2_arr, p2!, pdv_list
    REAL(hp), DIMENSION(:), ALLOCATABLE :: pdv_list
    REAL(bp), DIMENSION(8) :: temp ! LS uses this, will want to change sometime
    REAL(bp), DIMENSION(6,6) :: &
         covariance, &
         covariance_temp, &
         eigentestmatrix !remove later
    REAL(bp), DIMENSION(6) :: &
         elements, &
         p, &
         ran6, &
         eigentestd
    REAL(bp), DIMENSION(1,6) :: &
         mean1_temp,&
         mean2_temp
    REAL(bp), DIMENSION(2) :: &
         bounds
    REAL(bp) :: &
         integration_step, &
         mjd_tt, &
         peak, &
         probability_mass, &
         ran, &
         tmp, tt, &
         last_accepted_chi2
    REAL(hp) :: pdv, a_r, expo, last_accepted_expo
    TYPE (Time) :: t
    INTEGER, DIMENSION(:), ALLOCATABLE :: &
         indx_arr, &
         nobs_arr
    INTEGER :: &
         err, &
         i, &
         iorb, &
         itrial, &
         ntrial, &
         j, &
         k, &
         l, &
         nperturber, &
         nstorb, &
         nrot,&
         nchain, noutlier, &
         proposal_stage,&
         ncur
    LOGICAL, DIMENSION(:,:), POINTER :: &
         obs_masks
    LOGICAL, DIMENSION(10) :: &
         perturbers
    LOGICAL :: &
         accept, &
         last, &
         first, &
         burnin_done, &
         asteroid_perturbers, &
         chi2_compared

    ntrial = 10000000
    burnin_done = .FALSE.
    chi2_compared = .FALSE.
    nstorb = SIZE(storb_arr)
    nperturber = 0
    last_accepted_chi2 = 0.0_bp
    nchain = 1
    ncur = 1
    proposal_stage = 1

    DO i=1,SIZE(proposal_density_masses)
       IF (proposal_density_masses(i) > 0.0_bp) THEN
          nperturber = nperturber + 1
       END IF
    END DO

    ALLOCATE(additional_perturbers(nperturber,8), &
         additional_perturbers_(nperturber-1,8), &
         elements0_arr(nstorb,6), elements_arr(nstorb,6), &
         chi2_arr(nstorb),chi2_arr_updated(nstorb), rchi2_arr(nstorb), &
         nobs_arr(nstorb), obs_arr(nstorb), A_arr(nstorb,6,6), stat=err)
    ALLOCATE(pdv_list(norb))
    ALLOCATE(cov_matrix(6*nstorb+nperturber,6*nstorb+nperturber))
    ALLOCATE(ok_cov_matrix(6*nstorb+nperturber,6*nstorb+nperturber))
    ALLOCATE(previous_cov_matrix(6*nstorb+nperturber,6*nstorb+nperturber))
    ALLOCATE(mean_arr(6*nstorb+nperturber), old_mean(6*nstorb+nperturber))
    ALLOCATE(covariance2(nstorb,6,6))
    ALLOCATE(A_arr2(6*nstorb+nperturber,6*nstorb+nperturber))
    ALLOCATE(ya_temp(6*nstorb+nperturber,6*nstorb+nperturber))
    ALLOCATE(p2(6*nstorb+nperturber))
    ALLOCATE(another_temp(6*nstorb+nperturber))
    ALLOCATE(ran_arr(6*nstorb+nperturber))
    ALLOCATE(deviates(6*nstorb+nperturber),lambda_arr(1,6*nstorb+nperturber))
    ALLOCATE(deviates_matrix(1,6*nstorb+nperturber))
    ALLOCATE(elem_temp(1,6*nstorb+nperturber))
    ALLOCATE(last_proposal(nstorb*8+3))
    ALLOCATE(last_accepted_chi2_arr(nstorb))

    lambda_arr(1,:) = 2.4_bp**2/(6.0_bp*nstorb+nperturber) ! Use this with AM; default scale factor for all
    cov_matrix = 0.0_bp
    A_arr2 = 0.0_bp
    mean_arr = 0.0_bp
    nsigma_mass = 1
    nsigma_orbit = 1

    CALL NEW(mcmc_out_file, TRIM(out_fname))
    CALL OPEN(mcmc_out_file)

    CALL NEW(meanres_file, TRIM(out_fname)//trim('.mean.res'))
    CALL OPEN(meanres_file)

    CALL NEW(res_file, TRIM(out_fname)//trim('.res'))
    CALL OPEN(res_file)

    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / massEstimation_MCMC:", &
            "Could not allocate memory.", 1)
       RETURN
    END IF

    ! Extract integration parameters
    CALL getParameters(storb_arr(1), dyn_model=dyn_model, &
         perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, integrator=integrator, &
         integration_step=integration_step)

    ! extract observational and orbital information
    DO i=1,nstorb

       obs_arr(i) = getObservations(storb_arr(i)) ! Need to store these for later
       ! extract information on how many and which observations should
       ! be used in the analysis
       obs_masks => getObservationMasks(storb_arr(i))
       nobs_arr(i) = COUNT(obs_masks)
       DEALLOCATE(obs_masks)

       t = getTime(orb_arr(i))
       mjd_tt = getMJD(t,"TT")
       elements0_arr(i,:) = getElements(orb_arr(i), "cartesian", "equatorial")
       elements_arr(i,:) = elements0_arr(i,:)
       mean_arr(1+(i-1)*6:6+(i-1)*6) = elements_arr(i,:)
       ! prepare matrix needed for covariance sampling
       covariance = getCovarianceMatrix(storb_arr(i), "cartesian", "equatorial")
       DO j=1, 6
          WRITE(stderr, *) covariance(j, :)
       END DO
       WRITE(stderr, * ) "----------------------------"
       covariance2(i,:,:) = getCovarianceMatrix(storb_arr(i), "cartesian", "equatorial")
       ! Merge covariance matrices into one
       cov_matrix(1+(i-1)*6:6+(i-1)*6,1+(i-1)*6:6+(i-1)*6) = covariance(:,:) * lambda_arr(1,1)

       ! populate perturber array for integrations
       IF (i <= nperturber) THEN
          additional_perturbers(i,1:6) = elements0_arr(i,:)
          additional_perturbers(i,7) = getMJD(t, "TT")
          additional_perturbers(i,8) = proposal_density_masses(i)

       END IF

       CALL NULLIFY(t)
       CALL NULLIFY(obs_arr(i))
    END DO

    ! Add masses into cov matrix
    IF (.NOT. ASSOCIATED(input_cov_matrix)) THEN
      DO i=1, nperturber
         cov_matrix(6*nstorb+i, 6*nstorb+i) = proposal_density_masses(i)*1e-19_bp
      END DO
    END IF

    IF (ASSOCIATED(input_cov_matrix)) THEN
      cov_matrix(:,:) = input_cov_matrix(:,:)
    END IF

    ok_cov_matrix = cov_matrix
    A_arr2(:,:) = cov_matrix(:,:)

    DO i=1,nstorb*6+nperturber
       WRITE(stderr,*) cov_matrix(i,:)
    END DO

    WRITE(stderr, * ) "----------------------------"

    CALL cholesky_decomposition(A_arr2(:,:), p2, errstr)

    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("PhysicalParameters / massEstimation_MCMC:", &
            "Cholesky decomposition unsuccessful:", 1)
       WRITE(stderr,"(A)") TRIM(errstr)
       RETURN
    END IF

    DO j=1,6*nstorb+nperturber
       A_arr2(j,j) = p2(j)
    END DO

    iorb = iorb_init
    itrial = itrial_init
    first = .TRUE.

    DO
       IF ( PRESENT(nburn_arr) .EQV. .FALSE.) THEN
          burnin_done = .TRUE.
       END IF

       ncur = ncur + 1
       itrial = itrial + 1

       IF (PRESENT(nburn_arr)) THEN
          IF ((burnin_done .EQV. .FALSE.) .AND. (ncur > nburn_arr(nchain))) THEN
             burnin_done = .TRUE.
          END IF
       END IF

       DO i=1,nstorb
          ! Compute goodness of fit
          IF (i <= nperturber) THEN
             k = 0
             additional_perturbers_ = 0

             DO j=1,nperturber
                IF (j /= i) THEN
                   k = k + 1
                   additional_perturbers_(k,:) = additional_perturbers(j,:)
                END IF
             END DO

             CALL setParameters(orb_arr(i), &
                  dyn_model=dyn_model, &
                  perturbers=perturbers, &
                  asteroid_perturbers=asteroid_perturbers, &
                  integrator=integrator, &
                  integration_step=integration_step, &
                  mass=additional_perturbers(i,8))
          ELSE
             CALL setParameters(orb_arr(i), &
                  dyn_model=dyn_model, &
                  perturbers=perturbers, &
                  asteroid_perturbers=asteroid_perturbers, &
                  integrator=integrator, &
                  integration_step=integration_step)
          END IF
       END DO

       WRITE(stderr, *) "----------------------------------"
       chi2_arr = getChi2(storb_arr,orb_arr,residuals=resids)
       rchi2_arr = chi2_arr / (nobs_arr - nstorb*6 - nperturber)
       !       pdv = EXP(-0.5_bp*SUM(chi2_arr)/(SUM(nobs_arr)-nstorb*6-nperturber-1))
       !       expo = -0.5_bp*(SUM(chi2_arr)-SUM(nobs_arr)-nstorb*6-nperturber)
       expo = -0.5_hp*SUM(chi2_arr)
       !expo = -0.5_bp*(SUM(chi2_arr)/(2*SUM(nobs_arr)-nstorb*6-nperturber))
       pdv = EXP(expo)

       WRITE(stderr, *) "iorb:", iorb
       WRITE(stderr, *) "pdv:", pdv
       WRITE(stderr, *) "chi2:", SUM(chi2_arr)
       WRITE(stderr, *) "last accepted chi2:", last_accepted_chi2
       WRITE(stderr, *) expo

       IF (pdv == 0.0_bp) THEN
          pdv = TINY(pdv)
       END IF

       ran = 1.0_bp

       IF (first) THEN
          accept = .TRUE.
          last_accepted_chi2 = SUM(chi2_arr)
          last_accepted_chi2_arr = chi2_arr
          first = .FALSE.
       ELSE
          chi2_compared = .FALSE.
          a_r = 1.0_hp

          IF (expo - last_accepted_expo > 1000) THEN
             a_r = 1.0_hp
          ELSE IF  (expo - last_accepted_expo < -1000) THEN
             a_r = 0.0_hp
          ELSE
             a_r = EXP(expo - last_accepted_expo)
          END IF

          ! Alternative acceptance probability. See Siltala & Granvik (2021a)
          !DO j=1, nstorb
          !   a_r = a_r * EXP(-0.5_bp*(chi2_arr(j) - last_accepted_chi2_arr(j)))
          !END DO

          ! This catches NaN values (can occur with very bad proposals)
          ! Values of infinity may also occur but they are handled by the next
          ! if clause.
          IF (a_r /= a_r) THEN
             a_r = 0.0_hp
          END IF
          WRITE(stderr, *) "a_r:", a_r

          IF (a_r > 1.0_hp) THEN
             WRITE(stderr, *) "BETTER!"
             a_r = 1.0_hp
             accept = .TRUE.
          ELSE
             CALL randomNumber(ran)
             IF (ran < a_r) THEN
                accept = .TRUE.
             ELSE IF (NINT(accepted_solutions(8*nstorb+3,iorb)) == 100) THEN
                accept = .TRUE. !  Force accept new proposal if the chain is stuck for some reason.
             ELSE
                accept = .FALSE.
                WRITE(stderr, *) "NOT ACCEPTED!"
             END IF
          END IF
       END IF

       IF (accept .EQV. .TRUE. .AND. chi2_compared .EQV. .FALSE.) THEN

          CALL updateMeanResids(storb_arr,resids)

          IF (iorb == 0 .OR. MOD(iorb,500) == 0) THEN ! Writes the residuals every 500 accepted proposals.
             WRITE(getUnit(res_file), *) "#", iorb
             CALL writeResiduals_SO(storb_arr,resids,getUnit(res_file))
             WRITE(getUnit(meanres_file), *) "#", iorb
             CALL writeMeanResids(storb_arr,orb_arr,getUnit(meanres_file))
          END IF

          ! Outlier detection happens here. We also need to get
          ! updated chi2 values after.
          IF (((iorb /= 0) .AND. (MOD(iorb,500) == 0 )) .OR. iorb == 15) THEN
             CALL outlierDetection(storb_arr)

             IF (ASSOCIATED(resids%vectors(1)%elements)) THEN
                DEALLOCATE(resids%vectors)
             END IF

             chi2_arr = getChi2(storb_arr,orb_arr,residuals=resids)
             expo = -0.5_bp*SUM(chi2_arr)
             pdv = EXP(expo)
             noutlier = 0

             DO j=1,nstorb
                obs_arr(j) = getObservations(storb_arr(j)) ! Need to store these for later
                ! extract information on how many and which observations should
                ! be used in the analysis
                obs_masks => getObservationMasks(storb_arr(j))
                nobs_arr(j) = COUNT(obs_masks(:,2)) ! We meed to update the total
                noutlier = noutlier + SIZE(obs_masks(:,2)) - nobs_arr(j)
                DEALLOCATE(obs_masks)          ! amount of used observations after outlier rejection.
                ALLOCATE(mean_resids(getNrOfObservations(obs_arr(j)),6))
                CALL getMeanResids(storb_arr(j),mean_resids)
                information_matrix => getBlockDiagInformationMatrix(obs_arr(j))
                obs_masks => getObservationMasks(storb_arr(j))
                chi2_arr_updated(j) = chi_square(mean_resids, &
                     information_matrix, obs_masks, errstr)
                DEALLOCATE(obs_masks)
                DEALLOCATE(mean_resids)
                DEALLOCATE(information_matrix)
                CALL NULLIFY(obs_arr(j))
             END DO
          END IF

          iorb = iorb + 1
          first = .FALSE.
          last_accepted_expo = -0.5_bp*SUM(chi2_arr)
          last_accepted_chi2 = SUM(chi2_arr)
          last_accepted_chi2_arr = chi2_arr
          last_accepted_expo = expo
          WRITE(stderr, *) "ACCEPTED!"

          IF (delayed_rejection) THEN
             WRITE(stderr,*) "DR at stage", proposal_stage
          END IF

          IF (info_verb >= 3) THEN
             WRITE(stderr, *) "iorb: ", iorb
          END IF

          DO i=1,nstorb
             IF (i <= nperturber) THEN
                accepted_solutions(8*(i-1)+1:8*i,iorb) = &
                     (/ additional_perturbers(i,1:6), additional_perturbers(i,8), chi2_arr(i) /)
             ELSE
                accepted_solutions(8*(i-1)+1:8*i,iorb) = &
                     (/ elements_arr(i,1:6), proposal_density_masses(i), chi2_arr(i) /)
             END IF
          END DO

          accepted_solutions(8*nstorb+1:8*nstorb+3,iorb) = &
               (/ SUM(chi2_arr), REAL(pdv,8), 1.0_bp /)
          pdv_list(iorb) = pdv
       ELSE
          accepted_solutions(8*nstorb+3,iorb) = &
               accepted_solutions(8*nstorb+3,iorb) + 1.0_bp
       END IF

       DO j=1, nstorb
          DEALLOCATE(resids%vectors(j)%elements)
       END DO

       DEALLOCATE(resids%vectors)

       IF (info_verb >= 3) THEN
          IF (accept .EQV. .FALSE.) THEN
             WRITE(stderr, *) "NOT ACCEPTED!!!"
          END IF
       END IF

       IF (accept) THEN
          ! Print out the current proposal matrix every 500 accepted proposals.
          IF (MODULO(iorb,500) == 0) THEN 
             CALL NEW(cov_file, TRIM(out_fname)//trim('.cov'))
             CALL OPEN(cov_file)
             WRITE(getUnit(cov_file), *) "# iorb at", iorb, "and itrial at:", itrial, &
                  " and a total of", noutlier, "outliers. Prop. matrix follows:"
             DO i=1,6*nstorb+nperturber
                WRITE(getUnit(cov_file),*) cov_matrix(i,:)
             END DO
             CALL NULLIFY(cov_file)
             ok_cov_matrix = cov_matrix
          END IF

          IF ((.NOT. first) .and. (iorb > 1)) THEN ! This is how many repetitions there were. can't be printed until the next
            ! proposal gets accepted because we don't know the amount until then..
             WRITE(getUnit(mcmc_out_file), "(1(I5))", advance="yes") NINT(accepted_solutions(8*nstorb+3,iorb-1))
          END IF

          WRITE(getUnit(mcmc_out_file), "(1(I7,1X),1(I10,1X), 2(I5,1X), 1(F8.2,1X))", advance="NO") iorb, itrial, nperturber, &
                nstorb-nperturber, mjd_tt
          DO i=1,nstorb
             WRITE(getUnit(mcmc_out_file),"(A5,1X,1(3(F19.15,1X),3(F19.15,1X),3(E12.5,1X)))",advance="NO") &
                  getID(storb_arr(i)), &
                  accepted_solutions(8*(i-1)+1:8*i,iorb), rchi2_arr(i)!, & ! accepted: 6 fitted orbital elements, fitted mass, chi2
          END DO

          WRITE(getUnit(mcmc_out_file),"(3(F19.2,1X))", advance="NO") SUM(chi2_arr), &
                SUM(chi2_arr)/(SUM(nobs_arr)-6*nstorb-nperturber), EXP(-0.5*SUM(chi2_arr)/SUM(nobs_arr))
       END IF

       IF (info_verb >= 3) THEN
          WRITE(stderr, *) "-----------------------------------"
       END IF

       IF (iorb == norb) THEN
          EXIT
       END IF

       ! Generate new trial orbits half way through / if chain gets stuck
       old_mean = mean_arr
       IF (iorb == norb/2) THEN
          DO i=1,nperturber
             additional_perturbers(i,8) = proposal_density_masses(i)*2.0
          END DO
          DO i=1,nstorb
             ! When half-way through re-start chain from initial
             ! orbits
             elements_arr(i,:) = elements0_arr(i,:)
             first = .TRUE.
             IF (PRESENT(nburn_arr)) THEN
                burnin_done = .FALSE.
                ncur = 0
                nchain = 2
             END IF
          END DO
       END IF

       IF (iorb == 5000 .AND. accept .AND. adaptation == "gaswam") THEN
          adaptation = "ram" ! This block switches from aswam to ram.
          itrial = 1
          A_arr2(:,:) = lambda_arr(1,1)*A_arr2(:,:)
       END IF

       ! Adaptation occurs within this IF block.
       IF ((iorb - iorb_init) > 1 .AND. iorb .NE. norb/2) THEN
          IF (accept) THEN
            previous_cov_matrix = cov_matrix
          END IF
          ! ----------------------- RAM ADAPTATION -------------------

          IF (adaptation == "ram") THEN
             lambda_arr(1,:) = 1.0_bp
             deviates_matrix(1,:) = ran_arr(:)
             cov_matrix = MATMUL(MATMUL(A_arr2,(identity_matrix(6*nstorb+nperturber)+(itrial**(-0.5_bp) &
                  *(a_r - 0.237_bp) &
                  *  MATMUL(TRANSPOSE(deviates_matrix),deviates_matrix)/matnorm(deviates_matrix)**2))),TRANSPOSE(A_arr2))
          END IF

          ! ---------------------- ASWAM SCALE ADAPTATION -----------

          IF ((adaptation == "aswam" .OR. adaptation == "gaswam") .AND. (iorb - iorb_init) > 19) THEN
             lambda_arr(1,:) = EXP(LOG(lambda_arr(1,1))+itrial**(-0.5_bp)*(a_r-0.237_bp)) !"Optimal": 0.237
             WRITE(stderr, *) "lambda:", lambda_arr(1,1)

             IF (.NOT. lambda_arr(1,1) / lambda_arr(1,1) == 1) THEN ! Apparently this catches NaNs
               WRITE(stderr, *) "WARNING: NaN lambda detected! Resetting.."
               lambda_arr(1,:) = 0.237_bp
             END IF

          END IF

          ! ------------------ ASWAM/AM COVARIANCE MATRIX UPDATE ------

          IF (adaptation == "aswam" .OR. adaptation == "am" .OR. adaptation == "gaswam") THEN
             DO i=1, nstorb ! Orbital elements mean
                mean_arr(1+(i-1)*6:6+(i-1)*6) = ((itrial-1.0_bp)/itrial)*mean_arr(1+(i-1)*6:6+(i-1)*6) &
                     + 1.0_bp/itrial*accepted_solutions(8*(i-1)+1:8*(i-1)+6,iorb)
             END DO

             DO i=1, nperturber !Mean masses
                mean_arr(6*nstorb+i) = ((itrial-1.0_bp)/itrial)*mean_arr(6*nstorb+i) &
                     + 1.0_bp/itrial*accepted_solutions(8*(i-1)+7,iorb)
             END DO

             IF ((iorb - iorb_init) > 19 .AND. iorb .NE. norb/2) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stderr, *) "Adapting cov matrix.."
                   WRITE(stderr, *) "Old covariance matrix:"
                   DO j=1,6
                      WRITE(stderr, *) covariance2(i,j,:)
                   END DO
                END IF
                ya_temp = 0.0_bp
                another_temp = 0.0_bp
                DO j=1, iorb
                   DO i=1, nstorb
                      another_temp(1+(i-1)*6:6+(i-1)*6) =  &
                           accepted_solutions(8*(i-1)+1:8*(i-1)+6,j) - mean_arr(1+(i-1)*6:6+(i-1)*6)
                   END DO
                   DO i=1, nperturber ! masses
                      another_temp(6*nstorb+i) =  &
                           accepted_solutions(8*(i-1)+7,j) - mean_arr(6*nstorb+i)
                   END DO
                   elem_temp(1,:) = another_temp(:)
                   ya_temp = ya_temp + MATMUL(TRANSPOSE(elem_temp), elem_temp)
                END DO

                elem_temp(1,:) = another_temp(:)
                cov_matrix = 1.0_bp/(iorb - 1.0_bp) * &
                     ya_temp + 1e-23_bp*identity_matrix(6*nstorb+nperturber)
             END IF
          END IF

          IF (adaptation .NE. "none") THEN

             A_arr2(:,:) = cov_matrix(:,:)
             CALL cholesky_decomposition(A_arr2(:,:), p2, errstr)
             ! If Cholesky decomposition fails for some reason, fall back
             ! with an earlier proposal matrix.
             IF (LEN_TRIM(errstr) /= 0) THEN
                error = .TRUE.
                CALL errorMessage("PhysicalParameters / massEstimation_MCMC:", &
                     "Cholesky decomposition unsuccessful:", 1)
                WRITE(stderr,"(A)") TRIM(errstr)
                WRITE(stderr, *) "Previous cov. matrix follows:"
                DO j=1,6*nstorb+nperturber
                   WRITE(stderr, *) previous_cov_matrix(j,:)
                END DO
                WRITE(stderr, *) "Resuming from earlier matrix."
                A_arr2 = ok_cov_matrix
                CALL cholesky_decomposition(A_arr2, p2, errstr)
             END IF

             DO j=1,6*nstorb+nperturber
                A_arr2(j,j) = p2(j)
                A_arr2(j,j+1:6*nstorb+nperturber) = 0.0_bp
             END DO
          END IF
       END IF

       elem_temp(1,:) = p2(:)
       CALL randomGaussian(ran_arr)
       deviates = lambda_arr(1,:)*MATMUL(A_arr2(:,:),ran_arr)

       ! New coordinates = last accepted coordinates + deviates:
       ! Second stage proposal from Mira's symmetric delayed rejection algorithm.
       IF ((accept .EQV. .FALSE.) .AND. (delayed_rejection) .AND. (proposal_stage == 1)) THEN
          WRITE(0, *) "DEBUG: DR ON!"
          ! This is triggered when we want to do DR and the previous first-stage proposal was rejeceted.
          DO i=1,nstorb
             IF (i <= nperturber) THEN
                last_proposal(8*(i-1)+1:8*i) = &
                     (/ additional_perturbers(i,1:6), additional_perturbers(i,8), chi2_arr(i) /)
             ELSE
                last_proposal(8*(i-1)+1:8*i) = &
                     (/ elements_arr(i,1:6), proposal_density_masses(i), chi2_arr(i) /)
             END IF
          END DO
          last_proposal(8*nstorb+1:8*nstorb+3) = &
               (/ SUM(chi2_arr), REAL(pdv,8), 1.0_bp /)
          proposal_stage = 2
       ELSE
          ! Currently, we do 2-stage proposals with DR. If the second is also rejected, return to "zero-stage".
          proposal_stage = 1
          last_proposal(:) = accepted_solutions(:,iorb)
          !ELSE
       END IF

       DO i=1, nstorb
          elements_arr(i, :) = last_proposal(8*(i-1)+1:8*(i-1)+6) + deviates(1+(i-1)*6:6+(i-1)*6)
       END DO
       IF (iorb - iorb_init > 499) THEN 
         mass_lock = .FALSE. 
       END IF
       IF (mass_lock .EQV. .FALSE.) THEN
         DO i=1, nperturber
            additional_perturbers(i,8) = last_proposal(8*(i-1)+7) + deviates(6*nstorb+i)
            IF (additional_perturbers(i,8) < 0) THEN
               WRITE(stderr, *) "WARNING: NEGATIVE MASS FOUND"
            END IF
         END DO

         DO WHILE (ANY(additional_perturbers(:,8) < 0)) !Repeat if neg mass found---------------------
             CALL randomGaussian(ran_arr)
             deviates = lambda_arr(1,:)*MATMUL(A_arr2(:,:),ran_arr)
             DO i=1, nstorb
                elements_arr(i, :) = last_proposal(8*(i-1)+1:8*(i-1)+6) + deviates(1+(i-1)*6:6+(i-1)*6)
             END DO
             DO i=1, nperturber
                additional_perturbers(i,8) = last_proposal(8*(i-1)+7) + deviates(6*nstorb+i)
                IF (additional_perturbers(i,8) < 0) THEN
                   WRITE(stderr, *) "WARNING: NEGATIVE MASS FOUND"
                END IF
             END DO

         END DO
       END IF

       DO i=1, nstorb
          t = getTime(orb_arr(i))
          CALL NULLIFY(orb_arr(i))
          IF (i <= nperturber) THEN
             additional_perturbers(i,1:6) = elements_arr(i,1:6)
             CALL NEW(orb_arr(i), elements_arr(i,:), "cartesian", "equatorial", t, mass=additional_perturbers(i,8))
          ELSE
             CALL NEW(orb_arr(i), elements_arr(i,:), "cartesian", "equatorial", t)
          END IF
          CALL NULLIFY(t)
       END DO
    END DO

    ! ----------------------------- MCMC IS DONE -----------------------
    ! We have to print out the number of accepted proposals for the last one separately.
    WRITE(getUnit(mcmc_out_file), "(1(I5))", advance="yes") 1
    ! Post-fit computation of statistics for masses
    IF (last) THEN
       ALLOCATE(indx_arr(norb))
       WRITE(getUnit(mcmc_out_file), "(1(I5))", advance="yes") NINT(accepted_solutions(8*nstorb+3,iorb-1))
       DO i=1,nperturber
          ! Nominal:
          estimated_masses(i,1) = accepted_solutions(8*i-1,MAXLOC(accepted_solutions(8*nstorb+2,1:iorb),dim=1))
          ! N-sigma bounds:
          DO j=1,3
             IF (j == 1) THEN
                probability_mass = 0.67_bp
             ELSE IF (j == 2) THEN
                probability_mass = 0.95_bp
             ELSE IF (j == 3) THEN
                probability_mass = 0.9973_bp
             END IF
             CALL credible_region(accepted_solutions(8*nstorb+1,1:iorb), & ! pdf
                  probability_mass=probability_mass, &
                  indx_arr=indx_arr, &
                  errstr=errstr, &
                  repetition_arr=NINT(accepted_solutions(8*nstorb+3,1:iorb))) ! repetitions
             IF (LEN_TRIM(errstr) /= 0) THEN
                WRITE(stderr,"(2A)") "Error: ", TRIM(errstr)
                STOP
             END IF
             bounds = (/ HUGE(probability_mass), -HUGE(probability_mass) /)
             DO k=1,iorb
                IF (indx_arr(k) /= -1) THEN
                   IF (accepted_solutions(8*(i-1)+7,indx_arr(k)) < bounds(1)) THEN
                      bounds(1) = accepted_solutions(8*(i-1)+7,indx_arr(k))
                   ELSE IF (accepted_solutions(8*(i-1)+7,indx_arr(k)) > bounds(2)) THEN
                      bounds(2) = accepted_solutions(8*(i-1)+7,indx_arr(k))
                   END IF
                END IF
             END DO
             estimated_masses(i,2+2*(j-1):3+2*(j-1)) = bounds
             IF (j == 1) THEN
                WRITE(getUnit(mcmc_out_file),"(A,1X,E13.6)") "# Nominal: ", &
                     estimated_masses(i,1)
                WRITE(getUnit(mcmc_out_file),"(A,2(1X,E13.6))") "# 1-sigma: ", bounds
             ELSE IF (j == 2) THEN
                WRITE(getUnit(mcmc_out_file),"(A,2(1X,E13.6))") "# 2-sigma: ", bounds
             ELSE IF (j == 3) THEN
                WRITE(getUnit(mcmc_out_file),"(A,2(1X,E13.6))") "# 3-sigma: ", bounds
             END IF
          END DO
       END DO

       DO i=1,nstorb
          CALL NULLIFY(storb_arr(i))
          CALL NULLIFY(orb_arr(i))
       END DO
    END IF

    IF (last .EQV. .FALSE.) THEN
       i = 1
       DO i=1,nstorb
          t = getTime(orb_arr(i))
          temp = accepted_solutions(8*(i-1)+1:8*i, MAXLOC(accepted_solutions(8*nstorb+2, 1:iorb), dim=1))
          elements_arr(i, :) = temp(1:6)

          CALL NULLIFY(orb_arr(i))
          CALL NEW(orb_arr(i), elements_arr(i,:), "cartesian", "equatorial", t)
          t = getTime(orb_arr(i))
          nominal_arr(i) = copy(orb_arr(i))

       END DO

    END IF

    CALL NULLIFY(mcmc_out_file)
    CALL NULLIFY(meanres_file)
    CALL NULLIFY(res_file)
    DEALLOCATE(additional_perturbers, additional_perturbers_, &
         indx_arr, A_arr, elements_arr, &
         elements0_arr, chi2_arr, rchi2_arr, stat=err)

  END SUBROUTINE massEstimation_MCMC

  ! Asteroid mass estimation 'marching' algorithm.
  ! Note that this is an approximation at best and should not be used by itself for
  ! any serious work!
  ! See Siltala & Granvik (2017)
  ! https://doi.org/10.1016/j.icarus.2017.06.028
  ! Author: LS
  SUBROUTINE massEstimation_march(storb_arr, orb_arr, HG_arr, dyn_model, integrator, integration_step, &
       perturbers, asteroid_perturbers, mass, out_fname, resolution)
    IMPLICIT NONE
    TYPE (StochasticOrbit), DIMENSION(:), INTENT(inout) :: storb_arr
    TYPE (Orbit), DIMENSION(:), INTENT(inout)           :: orb_arr
    CHARACTER(len=FNAME_LEN), INTENT(in)                :: out_fname
    INTEGER, INTENT(INOUT)                              :: resolution
    REAL(bp), INTENT(in), DIMENSION(:,:)                :: HG_arr

    TYPE (Orbit), DIMENSION(:,:), POINTER               :: orb_arr2 => NULL()
    TYPE (Time)                                         :: t, t_prop
    TYPE (File)                                         :: march_out_file, res_file
    INTEGER, DIMENSION(:), ALLOCATABLE                  :: nobs_arr
    LOGICAL, DIMENSION(:,:), POINTER                    :: obs_masks
    REAL(bp), DIMENSION(:,:), POINTER                   :: stdev_arr_measur, residuals
    REAL(bp), DIMENSION(:,:), ALLOCATABLE               :: chi2_arr
    REAL(bp), DIMENSION(6)                              :: perturber_elems, rms6, elements
    REAL(bp), DIMENSION(1,8)                            :: additional_perturber
    REAL(bp), DIMENSION(2)                              :: chi2_arr2
    REAL(bp), DIMENSION(:), ALLOCATABLE                 :: marching_masses, chi2_sum_arr
    REAL(bp)                                            :: density, albedo, const, &
         rough_estimate, lower_mass_bound, upper_mass_bound, chi_sum, integration_step, &
         mass, avgstdev, rms1, rms2, best_chi, chi_perturber, step
    LOGICAL, DIMENSION(10)                              :: perturbers
    LOGICAL                                             :: print_data, asteroid_perturbers
    LOGICAL, DIMENSION(6)                               :: obs_masks2, obs_masks_true
    INTEGER                                             :: k, i, j, l, m, nobs, outliertot
    INTEGER, DIMENSION(6)                               :: n0
    CHARACTER(len=DYN_MODEL_LEN)                        :: dyn_model
    CHARACTER(len=INTEGRATOR_LEN)                       :: integrator
    REAL(bp), DIMENSION(:,:), ALLOCATABLE               :: residuals_, observed_coords, &
         computed_coords, elements_arr, stdev, mean
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: &
         obsy_ccoords => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER  :: &
         observed_scoords => NULL(), &
         computed_scoords => NULL()
    TYPE (Observations)                                 :: obss
    TYPE (SparseArray)                                  :: resids

    CALL NEW(march_out_file, TRIM(out_fname))
    CALL OPEN(march_out_file)

    CALL NEW(res_file, TRIM(out_fname)//trim('.res'))
    CALL OPEN(res_file)

    ALLOCATE(nobs_arr(size(orb_arr)))

    IF (resolution == 0) THEN
      resolution = 300 ! Default value if resolution is not given.
    END IF

    perturber_elems(:) = getElements(orb_arr(1), "cartesian", "equatorial")

    ! Rough mass estimate to sample in approximately right range:
    density = 2.5_bp ! g/cm3
    albedo = 0.15_bp
    const = 391222381.5_bp * pi * density * (1.0e12_bp / kg_solar) / (albedo*SQRT(albedo))
    rough_estimate = const*10**(-0.6_bp*HG_arr(1,1))
    WRITE(getUnit(march_out_file), *) "# Rough mass estimate for march:", rough_estimate
    lower_mass_bound = 0.1_bp * rough_estimate
    upper_mass_bound = 10.0_bp * rough_estimate

    ALLOCATE(marching_masses(resolution))
    ALLOCATE(chi2_arr(resolution,size(orb_arr)))
    ALLOCATE(chi2_sum_arr(resolution))

    ! Generates an array of evenly spaced masses.
    step = 1.0_bp/resolution * upper_mass_bound
    marching_masses(1) = lower_mass_bound

    DO i=1,resolution
      mass = lower_mass_bound + (i-1) * step
      marching_masses(i) = mass
    END DO

    DO i=1, resolution
       CALL setParameters(orb_arr(1), mass=marching_masses(i))
       chi2_arr(i,:) = getChi2(storb_arr, orb_arr, resids)
       chi2_sum_arr(i) = SUM(chi2_arr(i,:))
       WRITE(getUnit(march_out_file), *) marching_masses(i), chi2_arr(i,:), chi2_sum_arr(i), i
    END DO
    mass = marching_masses(MINLOC(chi2_sum_arr,dim=1))
    WRITE(getUnit(march_out_file), *) "# Best mass: ", mass
    ALLOCATE(stdev(2,2))
    ALLOCATE(mean(2,2))

    CALL writeResiduals_SO(storb_arr,orb_arr,getUnit(res_file))

  END SUBROUTINE massEstimation_march


END MODULE PhysicalParameters_cl
