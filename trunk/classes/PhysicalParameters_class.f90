!====================================================================!
!                                                                    !
! Copyright 2002,2003,2004,2005,2006,2007,2008,2009,2010,2011        !
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
!! magnitude H and slope parameter G.
!!
!! @see StochasticOrbit_class 
!!
!! @author  MG
!! @version 2011-01-26
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
  use statistics

  IMPLICIT NONE

  PRIVATE :: new_PP
  PRIVATE :: new_PP_storb
  PRIVATE :: nullify_PP
  PRIVATE :: copy_PP
  PRIVATE :: exist_PP

  TYPE PhysicalParameters

     PRIVATE
     REAL(bp), DIMENSION(:,:), POINTER :: H0_arr
     REAL(bp)                          :: H0_nominal
     REAL(bp)                          :: H0_unc
     REAL(bp), DIMENSION(:,:), POINTER :: G_arr
     REAL(bp)                          :: G_nominal
     REAL(bp)                          :: G_unc
     REAL(bp), DIMENSION(:,:), POINTER :: m_arr
     REAL(bp)                          :: m_nominal
     REAL(bp)                          :: m_unc     
     TYPE (StochasticOrbit)            :: storb
     LOGICAL                           :: is_initialized = .FALSE.

  END TYPE PhysicalParameters

  INTERFACE NEW
     MODULE PROCEDURE new_PP
     MODULE PROCEDURE new_PP_storb
  END INTERFACE

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_PP
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_PP
  END INTERFACE

  INTERFACE exist
     MODULE PROCEDURE exist_PP
  END INTERFACE

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

    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoord_arr 
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides_arr  
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr 
    TYPE (Orbit), DIMENSION(:), POINTER :: orb_arr 
    CHARACTER(len=2), DIMENSION(:), POINTER :: filter_arr
    CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: filter_arr_
    REAL(bp), DIMENSION(:,:,:), POINTER :: cov_arr
    REAL(bp), DIMENSION(:,:), POINTER :: pdfs_arr, phase_angle_arr
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: helio_dist_arr, topo_dist_arr
    REAL(bp), DIMENSION(:), POINTER :: obs_mag_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE :: H_arr, &
         obs_mag_arr_
    REAL(bp), DIMENSION(6) :: coordinates
    INTEGER :: err, i, j, nobs, norb, nmag

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
    IF (containsSampledPDF(this%storb)) THEN
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
    IF (containsSampledPDF(this%storb)) THEN
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
       END DO
       IF (info_verb >= 2) THEN
          WRITE(*,*) "H0_min: ", MINVAL(this%H0_arr(:,1)-this%H0_arr(:,2)), &
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
          WRITE(*,*) "H0: ", this%H0_nominal, "  dH0: ", this%H0_unc
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

    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoord_arr 
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides_arr  
    TYPE (Orbit), DIMENSION(:,:), POINTER :: orb_lt_corr_arr 
    CHARACTER(len=2), DIMENSION(:), POINTER :: filter_arr
    CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: filter_arr_
    REAL(bp), DIMENSION(:,:,:), POINTER :: cov_arr
    REAL(bp), DIMENSION(:,:), POINTER :: pdfs_arr, phase_angle_arr
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: helio_dist_arr, topo_dist_arr
    REAL(bp), DIMENSION(:), POINTER :: obs_mag_arr, obs_mag_unc_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE :: obs_mag_arr_, &
         obs_mag_unc_arr_
    REAL(bp), DIMENSION(6) :: coordinates
    INTEGER :: err, i, j, nmag, norb, nobs

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
    IF (containsSampledPDF(this%storb)) THEN
       CALL getEphemerides(this%storb, obsy_ccoord_arr, &
            ephemerides_arr, pdfs_arr=pdfs_arr, this_lt_corr_arr=orb_lt_corr_arr)
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

    ! Estimate H and G
    IF (containsSampledPDF(this%storb)) THEN
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
       END DO
       IF (info_verb >= 2) THEN
          WRITE(*,*) "H0_min: ", MINVAL(this%H0_arr(:,1)-this%H0_arr(:,2)), &
               "  H0_max: ", MAXVAL(this%H0_arr(:,1)+this%H0_arr(:,2))
          WRITE(*,*) "G_min:  ", MINVAL(this%G_arr(:,1)-this%G_arr(:,2)), &
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
       IF (info_verb >= 2) THEN
          WRITE(*,*) "H0: ", this%H0_nominal, "  dH0: ", this%H0_unc
          WRITE(*,*) "G:  ", this%G_nominal,  "  dG:  ", this%G_unc
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
  REAL(bp) FUNCTION getApparentMagnitude(H, G, r, Delta, phase_angle)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: H, G, r, Delta, phase_angle

    REAL(bp), DIMENSION(2) :: phi

    phi = HGPhaseFunctions(phase_angle)
    getApparentMagnitude = H - 2.5_bp*LOG10((1.0_bp - G)*phi(1) + &
         G*phi(2)) + 5.0_bp*LOG10(r*Delta)

  END FUNCTION getApparentMagnitude





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





  !! MCMC mass estimation algorithm
  SUBROUTINE massEstimation_mcmc(this, perturber_arr, proposal_density_masses, estimated_masses)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: this
    TYPE (StochasticOrbit), DIMENSION(:), INTENT(in) :: perturber_arr
    REAL(bp), DIMENSION(:), INTENT(in) :: proposal_density_masses
    REAL(bp), DIMENSION(:,:), POINTER :: estimated_masses

  END SUBROUTINE massEstimation_mcmc




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






END MODULE PhysicalParameters_cl
