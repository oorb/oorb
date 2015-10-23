!====================================================================!
!                                                                    !
! Copyright 2002-2014,2015                                           !
! Mikael Granvik, Jenni Virtanen, Karri Muinonen, Teemu Laakso,      !
! Dagmara Oszkiewicz, Grigori Fedorets                               !
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
!! Type and routines for orbits. Contains the most
!! fundamental routines used in orbital inversion. 
!!
!! @see StochasticOrbit_class 
!!
!! @author  MG, TL, KM, JV, GF
!! @version 2015-10-23
!!
MODULE Orbit_cl

  USE Base_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  USE CartesianCoordinates_cl

  USE utilities
  USE planetary_data
  USE integrators
  USE linal
  USE sort

  IMPLICIT NONE
  PRIVATE :: new_Orb
  PRIVATE :: new_Orb_cartesian
  PRIVATE :: new_Orb_spherical
  PRIVATE :: new_Orb_elements
  PRIVATE :: new_Orb_2point
  PRIVATE :: nullify_Orb
  PRIVATE :: copy_Orb
  PRIVATE :: exist_Orb
  PRIVATE :: GaussfgJacobian_Orb
  PRIVATE :: getCartesianElements
  PRIVATE :: getDelaunayElements
  PRIVATE :: getEphemeris_Orb_single
  PRIVATE :: getEphemeris_Orb_multiple
  PRIVATE :: getEphemerides_Orb_single
  PRIVATE :: getEphemerides_Orb_multiple
  PRIVATE :: getFrame_Orb
  PRIVATE :: getKeplerianElements
  PRIVATE :: getParameters_Orb
  PRIVATE :: getPhaseAngles_Orb
  PRIVATE :: getPoincareElements
  PRIVATE :: getPosition_Orb
  PRIVATE :: getPhaseAngle_Orb
  PRIVATE :: getSCoord_Orb
  PRIVATE :: getStumpffFunctions
  PRIVATE :: getTime_Orb
  PRIVATE :: getVelocity_Orb
  PRIVATE :: opposite_Orb
  PRIVATE :: propagate_Orb_single
  PRIVATE :: propagate_Orb_multiple
  PRIVATE :: reallocate_Orb
  PRIVATE :: reallocate_Orb_1
  PRIVATE :: reallocate_Orb_2
  PRIVATE :: rotateToEcliptic_Orb
  PRIVATE :: rotateToEquatorial_Orb
  PRIVATE :: estimateCosDf
  PRIVATE :: setParameters_Orb
  PRIVATE :: switchCenter_Orb
  !  PRIVATE :: solveKeplerEquation_stumpff
  !  PRIVATE :: solveKeplerEquation_newton
  PRIVATE :: toCartesian_Orb
  PRIVATE :: toCometary_Orb
  PRIVATE :: toKeplerian_Orb

  TYPE Orbit

     !     PRIVATE
     TYPE (Time)                         :: t
     CHARACTER(len=FRAME_LEN)            :: frame                = "equatorial"
     CHARACTER(len=ELEMENT_TYPE_LEN)     :: element_type         = "cartesian"
     REAL(bp), DIMENSION(6)              :: elements
     LOGICAL                             :: is_initialized       = .FALSE.
     ! Central body, default is Sun (for body id's, see module planetary_data)
     INTEGER                             :: center         = 11
     ! Parameters for propagation:
     CHARACTER(len=DYN_MODEL_LEN)        :: dyn_model_prm        = "2-body"
     CHARACTER(len=INTEGRATOR_LEN)       :: integrator_prm       = "bulirsch-stoer"
     REAL(bp), DIMENSION(:,:), POINTER   :: additional_perturbers => NULL() ! (car equ + mjdtt + mass)
     REAL(bp), DIMENSION(6)              :: finite_diff_prm      = -1.0_bp
     REAL(bp)                            :: integration_step_prm = 5.0_bp 
     LOGICAL, DIMENSION(10)              :: perturbers_prm       = .FALSE.
     ! Mass of the asteroid in solar masses (M_sol). A negative value
     ! indicates that this orbit is integrated as a test particle, ie,
     ! it has no effect on the orbits of other objects. Default is -1.
     REAL(bp)                            :: mass_prm             = -1.0_bp

  END TYPE Orbit

  INTERFACE NEW
     MODULE PROCEDURE new_Orb
     MODULE PROCEDURE new_Orb_cartesian
     MODULE PROCEDURE new_Orb_spherical
     MODULE PROCEDURE new_Orb_elements
     MODULE PROCEDURE new_Orb_2point
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_Orb
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_Orb
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_Orb
  END INTERFACE exist

  INTERFACE GaussfgJacobian
     MODULE PROCEDURE GaussfgJacobian_Orb
  END INTERFACE GaussfgJacobian

  INTERFACE getApoapsisDistance
     MODULE PROCEDURE getApoapsisDistance_Orb
  END INTERFACE getApoapsisDistance

  INTERFACE getCCoord
     MODULE PROCEDURE getCCoord_Orb
  END INTERFACE getCCoord

  INTERFACE getEphemeris
     MODULE PROCEDURE getEphemeris_Orb_single
     MODULE PROCEDURE getEphemeris_Orb_multiple
  END INTERFACE getEphemeris

  INTERFACE getEphemerides
     MODULE PROCEDURE getEphemerides_Orb_single
     MODULE PROCEDURE getEphemerides_Orb_multiple
  END INTERFACE getEphemerides

  INTERFACE getFrame
     MODULE PROCEDURE getFrame_Orb
  END INTERFACE getFrame

  INTERFACE getJacobiConstants
     MODULE PROCEDURE getJacobiConstants_Orb
  END INTERFACE getJacobiConstants

  INTERFACE getParameters
     MODULE PROCEDURE getParameters_Orb
  END INTERFACE getParameters

  INTERFACE getPeriapsisDistance
     MODULE PROCEDURE getPeriapsisDistance_Orb
  END INTERFACE getPeriapsisDistance

  INTERFACE getPhaseAngle
     MODULE PROCEDURE getPhaseAngle_Orb
  END INTERFACE getPhaseAngle

  INTERFACE getPhaseAngles
     MODULE PROCEDURE getPhaseAngles_Orb
  END INTERFACE getPhaseAngles

  INTERFACE getPosition
     MODULE PROCEDURE getPosition_Orb
  END INTERFACE getPosition

  INTERFACE getSCoord
     MODULE PROCEDURE getSCoord_Orb
  END INTERFACE getSCoord

  INTERFACE getSolarElongation
     MODULE PROCEDURE getSolarElongation_Orb
  END INTERFACE getSolarElongation

  INTERFACE getTime
     MODULE PROCEDURE getTime_Orb
  END INTERFACE getTime

  INTERFACE getTisserandsParameters
     MODULE PROCEDURE getTisserandsParameters_Orb
  END INTERFACE getTisserandsParameters

  INTERFACE getVelocity
     MODULE PROCEDURE getVelocity_Orb
  END INTERFACE getVelocity

  INTERFACE propagate
     MODULE PROCEDURE propagate_Orb_single
     MODULE PROCEDURE propagate_Orb_multiple
  END INTERFACE propagate

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_Orb
     MODULE PROCEDURE reallocate_Orb_1
     MODULE PROCEDURE reallocate_Orb_2
  END INTERFACE reallocate

  INTERFACE rotateToEquatorial
     MODULE PROCEDURE rotateToEquatorial_Orb
  END INTERFACE rotateToEquatorial

  INTERFACE rotateToEcliptic
     MODULE PROCEDURE rotateToEcliptic_Orb
  END INTERFACE rotateToEcliptic

  INTERFACE setParameters
     MODULE PROCEDURE setParameters_Orb
  END INTERFACE setParameters

  INTERFACE solveKeplerEquation
     MODULE PROCEDURE solveKeplerEquation_stumpff
     MODULE PROCEDURE solveKeplerEquation_newton
  END INTERFACE solveKeplerEquation

  INTERFACE switchCenter
     MODULE PROCEDURE switchCenter_Orb
  END INTERFACE switchCenter

  INTERFACE opposite
     MODULE PROCEDURE opposite_Orb
  END INTERFACE opposite

  INTERFACE toCartesian
     MODULE PROCEDURE toCartesian_Orb
  END INTERFACE toCartesian

  INTERFACE toCometary
     MODULE PROCEDURE toCometary_Orb
  END INTERFACE toCometary

  INTERFACE toKeplerian
     MODULE PROCEDURE toKeplerian_Orb
  END INTERFACE toKeplerian

CONTAINS





  !! *Description*:
  !!
  !! Initializes a Orbit-object. Position and velocity are
  !! zero-vectors. Coordinate frame is equatorial and propagation
  !! scheme is 2-body. Central body is the Sun.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout) :: this

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%elements        = 0.0_bp
    this%element_type    = "cartesian"
    this%frame           = "equatorial"
    CALL NULLIFY(this%t)
    this%center    = 11
    this%dyn_model_prm   = "2-body"
    this%finite_diff_prm = -1.0_bp
    this%is_initialized  = .TRUE.

  END SUBROUTINE new_Orb





  !! *Description*:
  !!
  !! Initializes a Orbit-object using a given
  !! CartesianCoordinates-object. Propagation scheme
  !! is 2-body and central body is the Sun.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Orb_cartesian(this, ccoord, center)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)             :: this
    TYPE (CartesianCoordinates), INTENT(in) :: ccoord
    INTEGER, INTENT(in), OPTIONAL           :: center

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    this%elements = getCoordinates(ccoord)
    IF (error) THEN
       CALL errorMessage("Orbit / new",&
            "TRACE BACK 1", 1)
       RETURN
    END IF
    this%element_type = "cartesian"
    this%frame = getFrame(ccoord)
    NULLIFY(this%additional_perturbers)
    this%t = getTime(ccoord)
    IF (PRESENT(center)) THEN
       this%center = center
    ELSE
       this%center = 11
    END IF
    this%dyn_model_prm = "2-body"
    this%finite_diff_prm = -1.0_bp
    this%is_initialized = .TRUE.

  END SUBROUTINE new_Orb_cartesian





  !! *Description*:
  !!
  !! Initializes a Orbit-object using given elements. Propagation
  !! scheme is 2-body. Optional argument can be used to define the
  !! central body other than the Sun.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Orb_elements(this, elements, element_type, frame, t, center, mass)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)        :: this
    REAL(bp), DIMENSION(6), INTENT(in) :: elements
    CHARACTER(len=*), INTENT(in)       :: element_type, frame
    TYPE (Time), INTENT(in)            :: t
    INTEGER, INTENT(in), OPTIONAL      :: center
    REAL(bp), INTENT(in), OPTIONAL     :: mass

    TYPE (SphericalCoordinates)        :: scoord
    TYPE (CartesianCoordinates)        :: ccoord
    CHARACTER(len=32)                  :: tmp

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(t)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Given epoch has not been initialized.", 1)
       RETURN
    END IF

    this%t = copy(t)
    this%dyn_model_prm = "2-body"
    this%finite_diff_prm = -1.0_bp
    tmp = element_type
    CALL locase(tmp, error)
    IF (error) THEN
       CALL errorMessage("Orbit / new", &
            "The element type string contains forbidden characters.", 1)
       RETURN
    END IF
    this%element_type = TRIM(tmp)

    IF (PRESENT(center)) THEN
       this%center = center
    ELSE
       this%center = 11
    END IF

    IF (PRESENT(mass)) THEN
       this%mass_prm = mass
    ELSE
       this%mass_prm = -1.0_bp
    END IF

    SELECT CASE (this%element_type)

    CASE ("cartesian")

       this%elements = elements
       tmp = frame
       CALL locase(tmp, error)
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "The frame string contains forbidden characters (5).", 1)
          RETURN
       END IF
       this%frame = TRIM(tmp)
       this%is_initialized = .TRUE.

    CASE ("spherical")

       tmp = frame
       CALL locase(tmp, error)
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "The frame string contains forbidden characters (10).", 1)
          RETURN
       END IF
       this%frame = TRIM(tmp)
       CALL NEW(scoord, elements, this%frame, this%t)
       CALL NEW(ccoord, scoord)
       this%elements = getCoordinates(ccoord)
       CALL NULLIFY(ccoord)
       CALL NULLIFY(scoord)
       IF (error) THEN
          CALL errorMessage("Orbit / new",&
               "TRACE BACK 1", 1)
          RETURN
       END IF
       this%element_type = "cartesian"
       this%is_initialized = .TRUE.

    CASE ("keplerian")

       ! Check soundness of Keplerian elements:
       IF (elements(1) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Semimajor Axis is negative.", 1)
          RETURN
       ELSE
          this%elements(1) = elements(1)
       END IF

       IF (elements(2) == 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is exactly 1 (not possible, too few digits used).", 1)
          RETURN
       ELSE IF (elements(2) < 0.0_bp .OR. elements(2) > 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is outside the range [0,1[.", 1)
          RETURN
       ELSE
          this%elements(2) = elements(2)
       END IF

       IF (elements(3) < 0.0_bp .OR. elements(3) > pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Inclination is outside the range [0,pi].", 1)
          RETURN
       ELSE
          this%elements(3) = elements(3)
       END IF

       IF (elements(4) < 0.0_bp .OR. elements(4) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Longitude of Ascending Node is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(4) = elements(4)
       END IF

       IF (elements(5) < 0.0_bp .OR. elements(5) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Argument of Perihelion is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(5) = elements(5)
       END IF

       IF (elements(6) < 0.0_bp .OR. elements(6) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Mean Anomaly is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(6) = elements(6)
       END IF

       this%element_type = "keplerian"
       this%frame = frame
       this%is_initialized = .TRUE.

    CASE ("cometary")

       ! Check soundness of cometary elements:
       IF (elements(1) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Perihelion distance is negative.", 1)
          RETURN
       ELSE
          this%elements(1) = elements(1)
       END IF

       IF (elements(2) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is negative.", 1)
          RETURN
       ELSE IF (elements(2) == 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is exactly 1 (not possible, too few digits used).", 1)
          RETURN
       ELSE
          this%elements(2) = elements(2)
       END IF

       IF (elements(3) < 0.0_bp .OR. elements(3) > pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Inclination is outside the range [0,pi].", 1)
          RETURN
       ELSE
          this%elements(3) = elements(3)
       END IF

       IF (elements(4) < 0.0_bp .OR. elements(4) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Longitude of Ascending Node is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(4) = elements(4)
       END IF

       IF (elements(5) < 0.0_bp .OR. elements(5) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Argument of Perihelion is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(5) = elements(5)
       END IF

       ! Time of perihelion can be anything:
       this%elements(6) = elements(6)

       this%element_type = "cometary"
       this%frame = frame
       this%is_initialized = .TRUE.

    CASE ("cometary_ta")

       !! Same as 'cometary', but time of perihelion changed to the
       !! true anomaly.

       ! Check soundness of cometary elements:
       IF (elements(1) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Perihelion distance is negative.", 1)
          RETURN
       ELSE
          this%elements(1) = elements(1)
       END IF

       IF (elements(2) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is negative.", 1)
          RETURN
       ELSE IF (elements(2) == 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is exactly 1 (not possible, too few digits used).", 1)
          RETURN
       ELSE
          this%elements(2) = elements(2)
       END IF

       IF (elements(3) < 0.0_bp .OR. elements(3) > pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Inclination is outside the range [0,pi].", 1)
          RETURN
       ELSE
          this%elements(3) = elements(3)
       END IF

       IF (elements(4) < 0.0_bp .OR. elements(4) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Longitude of Ascending Node is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(4) = elements(4)
       END IF

       IF (elements(5) < 0.0_bp .OR. elements(5) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Argument of Perihelion is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(5) = elements(5)
       END IF

       IF (elements(6) < 0.0_bp .OR. elements(6) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "True Anomaly is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(6) = elements(6)
       END IF

       this%element_type = "cometary_ta"
       this%frame = frame
       this%is_initialized = .TRUE.

    CASE ("cometary_ma")

       !! Same as 'cometary', but time of perihelion changed to the
       !! mean anomaly.

       ! Check soundness of cometary elements:
       IF (elements(1) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Perihelion distance is negative.", 1)
          RETURN
       ELSE
          this%elements(1) = elements(1)
       END IF

       IF (elements(2) < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is negative.", 1)
          RETURN
       ELSE IF (elements(2) == 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Eccentricity is exactly 1 (not possible, too few digits used).", 1)
          RETURN
       ELSE
          this%elements(2) = elements(2)
       END IF

       IF (elements(3) < 0.0_bp .OR. elements(3) > pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Inclination is outside the range [0,pi].", 1)
          RETURN
       ELSE
          this%elements(3) = elements(3)
       END IF

       IF (elements(4) < 0.0_bp .OR. elements(4) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Longitude of Ascending Node is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(4) = elements(4)
       END IF

       IF (elements(5) < 0.0_bp .OR. elements(5) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Argument of Perihelion is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(5) = elements(5)
       END IF

       IF (elements(6) < 0.0_bp .OR. elements(6) > two_pi) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Mean Anomaly is outside the range [0,two_pi].", 1)
          RETURN
       ELSE
          this%elements(6) = elements(6)
       END IF

       this%element_type = "cometary_ma"
       this%frame = frame
       this%is_initialized = .TRUE.

    CASE default

       error = .TRUE.
       CALL errorMessage("Orbit / new",&
            "Type of elements (" // TRIM(this%element_type) // &
            ") is incorrect.", 1)
       RETURN       

    END SELECT

    NULLIFY(this%additional_perturbers)


  END SUBROUTINE new_Orb_elements





  !! *Description*:
  !!
  !! Derives Cartesian orbital elements (this) from two positions
  !! (ccoord0 and ccoord1) using either (i) the p-iteration method
  !! (2-body) by Herrick and Liu, (ii) the continued fraction method
  !! by Hansen (2-body), or (iii) the amoeba method by Granvik and
  !! Muinonen (n-body). The frame of the output orbit is defined by
  !! the equal frames of the input positions. The Sun is the central
  !! body.
  !! 
  !! Returns error.
  !!
  RECURSIVE SUBROUTINE new_Orb_2point(this, ccoord0, ccoord1, &
       method, smamax, ftol, iter, perturbers, integrator, &
       integration_step, center)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                       :: this
    TYPE (CartesianCoordinates), INTENT(in)           :: ccoord0, ccoord1
    CHARACTER(len=*), INTENT(in)                      :: method
    REAL(bp), INTENT(in)                              :: smamax
    REAL(bp), INTENT(in), OPTIONAL                    :: ftol
    INTEGER(ibp), INTENT(out), OPTIONAL               :: iter
    LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL       :: perturbers
    CHARACTER(len=*), INTENT(in), OPTIONAL            :: integrator
    REAL(bp), INTENT(in), OPTIONAL                    :: integration_step
    INTEGER, INTENT(in), OPTIONAL                     :: center

    REAL(bp), PARAMETER :: dp_init = 1.0_bp !1.0_bp
    REAL(bp), PARAMETER :: p_init = 1.0e-7_bp
    REAL(bp), PARAMETER :: tolp = 1.0e-8_bp !1.0e-13_bp
    REAL(bp), PARAMETER :: told = 1.0e-16_bp !1.0e-16
    REAL(bp), PARAMETER :: tol = 1.0e-10_bp !1.0e-13_bp
    REAL(bp), PARAMETER :: toli = 1.0e-8_bp
    ! Require 100m accuracy (comparable to better than 0.5" astrometric
    ! accuracy when the topocentric distance is more than 40,000 km) in
    ! the n-body amoeba routine:
    REAL(bp), PARAMETER :: ftol_prm = 100.0_bp/m_au
    INTEGER(ibp), PARAMETER :: itmax_prm = 1000!500
    TYPE (Orbit), DIMENSION(4) :: orb_arr
    TYPE (Orbit) :: orb
    TYPE (Time) :: t0, t1
    TYPE (CartesianCoordinates) :: ccoord_
    CHARACTER(len=FRAME_LEN) :: frame
    CHARACTER(len=64) :: method_
    REAL(bp), DIMENSION(4,3) :: p_matrix
    REAL(bp), DIMENSION(4) :: y_vector
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp), DIMENSION(3) :: psum, pos0, vel0, pos1, orbit_normal, &
         ran
    REAL(bp) :: r0, r1, r01, cdf, sdf, diff1, diff2, p, p1, p2, dp, &
         cos_df, f, g, kappa, tau, h, cf, dcf, len_orbit_normal, ftol_
    INTEGER :: y, i, ntest, nsec, ntmax, nsmax, ihi, ilo, ndim, &
         err_verb_tmp, iter_, center_

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(ccoord0)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "First position is unavailable.", 1)
       RETURN
    END IF

    IF (.NOT.exist(ccoord1)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Second position is unavailable.", 1)
       RETURN
    END IF

    frame = getFrame(ccoord0)
    IF (frame /= getFrame(ccoord1)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Frames of input coordinates are different.", 1)
       RETURN
    END IF

    t0 = getTime(ccoord0)
    t1 = getTime(ccoord1)
    ccoord_ = copy(ccoord0)
    pos0 = getPosition(ccoord_)
    ! Initialize this object (at least) temporarily:
    CALL NEW(this, ccoord_)
    IF (error) THEN
       CALL errorMessage("Orbit / new", &
            "Could not initialize this object temporarily.", 1)
       RETURN
    END IF
    CALL NULLIFY(ccoord_)
    ccoord_ = copy(ccoord1)
    pos1 = getPosition(ccoord_)
    r0  = SQRT(DOT_PRODUCT(pos0,pos0))
    r1  = SQRT(DOT_PRODUCT(pos1,pos1))
    cdf = DOT_PRODUCT(pos0,pos1) / (r0*r1)
    IF (ABS(cdf) > 1.0_bp) THEN
       cdf = SIGN(1.0_bp, cdf)
    END IF

    IF (PRESENT(perturbers)) THEN
       this%perturbers_prm = perturbers
    ELSE
       this%perturbers_prm = .TRUE.
    END IF
    IF (PRESENT(integrator)) THEN
       this%integrator_prm = integrator
    ELSE
       this%integrator_prm = "bulirsch-stoer"
    END IF
    IF (PRESENT(integration_step)) THEN
       this%integration_step_prm = integration_step
    ELSE
       this%integration_step_prm = 1.0_bp
    END IF
    IF (PRESENT(center)) THEN
       this%center = center
       center_ = center
    ELSE
       this%center = 11
       center_ = 11
    END IF

    method_ = method
    CALL locase(method_, error)
    IF (error) THEN
       CALL errorMessage("Orbit / new", &
            "The method string contains forbidden characters.", 1)
       RETURN
    END IF
    SELECT CASE (method_)

    CASE ("p-iteration")

       !!
       !! p-iteration
       !!
       !! The orbital segment is not assumed to include the attracting 
       !! focus by default (y=1). However, if the computation is unsuccesful, 
       !! the attracting focus will be included in the orbital segment and
       !! the computations will be done again.
       !!
       !! Reference: J. M. A. Danby, "Fundamentals of celestial mechanics", 1992
       !!
       y = 1

       DO i=1,2

          sdf = y * SQRT(1.0_bp - cdf**2.0_bp)

          ! First bracketing:
          diff1 = 0.0_bp
          ntest = 1
          p2 = p_init
          dp = dp_init
          ! Determine maximum number of iteration steps needed to cover
          ! the phase space region defined by maximum allowed semimajor axis
          ! (maximum value for orbital parameter):
          ntmax = INT(smamax/dp)
          cos_df = estimateCosDf(this, ccoord_, p2, y)
          IF (error) THEN
             CALL errorMessage("Orbit / new", &
                  "TRACE BACK 5", 1)
             CALL NULLIFY(this)
             RETURN
          END IF
          diff2 = cdf - cos_df
          DO WHILE (diff1*diff2 >= 0.0_bp)
             ntest = ntest + 1
             IF (ntest > ntmax) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / new", &
                     "Too many steps in bracketing (1).", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             p1    = p2
             p2    = p1 + dp
             cos_df = estimateCosDf(this, ccoord_, p2, y)
             IF (error) THEN
                CALL errorMessage("Orbit / new", &
                     "TRACE BACK 10", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             diff1 = diff2
             diff2 = cdf - cos_df
          END DO

          ! Second bracketing:
          diff2 = diff1
          ntest = 1
          p2 = p1
          dp = dp_init/10.0_bp
          ! Determine maximum number of iteration steps needed to cover
          ! the phase space region defined by maximum allowed semimajor axis
          ! (maximum value for orbital parameter):
          ntmax = INT(smamax/dp)
          DO WHILE (diff1*diff2 >= 0.0_bp)
             ntest = ntest + 1
             IF (ntest > ntmax) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / new", &
                     "Too many steps in bracketing (2).", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             p1    = p2
             p2    = p1 + dp
             cos_df = estimateCosDf(this, ccoord_, p2, y)
             IF (error) THEN
                CALL errorMessage("Orbit / new", &
                     "TRACE BACK 15", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             diff1 = diff2
             diff2 = cdf - cos_df
          END DO

          ! Secant method:
          nsec = 0
          ! Determine maximum number of iteration steps needed to cover
          ! the phase space region defined by maximum allowed semimajor axis
          ! (maximum value for orbital parameter):
          !nsmax = INT(smamax/dp_init)
          nsmax = INT(smamax/dp)
          cos_df = estimateCosDf(this, ccoord_, p2, y)
          IF (error) THEN
             CALL errorMessage("Orbit / new", &
                  "TRACE BACK 20", 1)
             CALL NULLIFY(this)
             RETURN
          END IF
          diff2 = cdf - cos_df
          IF (ABS(diff2-diff1) < told) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / new", &
                  "Step size approaches (+/-)-infinity (1).", 1)
             CALL NULLIFY(this)
             RETURN
          END IF
          dp = -diff2*(p2-p1)/(diff2-diff1)
          DO WHILE (ABS(dp) >= tolp) 
             nsec = nsec+1
             IF (nsec > nsmax .OR. p2+dp <= 0.0_bp) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / new", &
                     "Too many steps in iteration.", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             p1 = p2
             p2 = p2 + dp
             cos_df = estimateCosDf(this, ccoord_, p2, y)
             IF (error) THEN
                CALL errorMessage("Orbit / new", &
                     "TRACE BACK 25", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             diff1 = diff2
             diff2 = cdf - cos_df
             IF (ABS(diff2-diff1) < told) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / new", &
                     "Step size approaches (+/-)-infinity (2).", 1)
                CALL NULLIFY(this)
                RETURN
             END IF
             dp = -diff2*(p2-p1)/(diff2-diff1)
          END DO
          p = p2
          f = (r1/p)*(cdf - 1.0_bp) + 1.0_bp
          g = ((r0*r1) / SQRT(planetary_mu(this%center)*p)) * sdf

          IF (g*y >= 0.0_bp) THEN
             ! Found correct orbital solution.
             EXIT
          ELSE
             ! Found incorrect orbital solution; include attracting focus
             ! in the orbital segment and redo the computation:
             y = -1 * y
          END IF

       END DO

       IF (g*y < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Could not find correct orbit.", 1)
          CALL NULLIFY(this)
          RETURN
       ELSE IF (ABS(g) < TINY(g)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Attempted division by zero (g-function).", 1)
          CALL NULLIFY(this)
          RETURN
       ELSE
          vel0 = (pos1 - f*pos0)/g
          this%elements(1:3) = pos0
          this%elements(4:6) = vel0
       END IF


    CASE ("continued fraction")

       !!
       !! continued fraction
       !!
       !! Derives Cartesian orbital elements (this) from two Cartesian
       !! positions (ccoord0 and ccoord1) using Hansen's continued
       !! fraction for the sector-triangle ratio (2-body).
       !!
       !! Reference: A. D. Dubyago, "The determination of orbits", 1961
       !!

       r01 = DOT_PRODUCT(pos0,pos1)
       tau = SQRT(planetary_mu(this%center))*(getMJD(t1,"TT") - getMJD(t0,"TT"))
       IF (ABS(tau) < tol) THEN
          ! t0 and t1 practically the same.
          CALL errorMessage("Orbit / new", &
               "Epochs t0 and t1 practically the same.", 1)
          RETURN
       END IF
       kappa = SQRT(2.0_bp*(r0*r1+r01))
       orbit_normal = cross_product(pos0,pos1)
       len_orbit_normal = SQRT(DOT_PRODUCT(orbit_normal,orbit_normal))
       IF (len_orbit_normal < toli) THEN
          ! r0 and r1 almost parallel.
          CALL errorMessage("Orbit / new", &
               "Vectors r0 and r1 almost parallel.", 1)
          RETURN
       END IF
       h = tau**2/(kappa**2*(kappa/3.0_bp+(r0+r1)/2.0_bp))
       CALL continuedFraction(h, tol, cf, dcf)
       IF (dcf > tol) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Failure in continued fraction", 1)
          RETURN
       END IF
       p = (len_orbit_normal*cf/tau)**2.0_bp
       f = (r1/p)*(cdf - 1.0_bp) + 1.0_bp
       sdf = SQRT(DOT_PRODUCT(orbit_normal,orbit_normal))/(r0*r1)
       g = ((r0*r1) / SQRT(planetary_mu(this%center)*p)) * sdf
       IF (ABS(g) < TINY(g)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / new", &
               "Attempted division by zero (g-function).", 1)
          CALL NULLIFY(this)
          RETURN
       END IF
       vel0 = (pos1 - f*pos0)/g
       this%elements(1:3) = pos0
       this%elements(4:6) = vel0


    CASE ("2-body amoeba")

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

       this%dyn_model_prm = "2-body"
       IF (PRESENT(ftol)) THEN
          ftol_ = ftol
       ELSE
          ftol_ = ftol_prm
       END IF
       elements = getCoordinates(ccoord0)
       CALL NEW(orb_arr(1), elements, "cartesian", frame, t0, center=center_)
       CALL setParameters(orb_arr(1), &
            dyn_model=this%dyn_model_prm, &
            perturbers=this%perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)
       y_vector(1) = distance_2b(orb_arr(1))
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
       IF (y_vector(1) < ftol_) THEN
          this = copy(orb)
          this%dyn_model_prm = "2-body"
          this%perturbers_prm = perturbers
          this%integrator_prm = integrator
          this%integration_step_prm = integration_step
          CALL NULLIFY(orb)
          IF (PRESENT(iter)) THEN
             iter = 0
          END IF
          RETURN
       END IF
       p_matrix(1,1:3) = elements(4:6)
       CALL NULLIFY(orb)
       DO i=2,4
          p_matrix(i,1:3) = elements(4:6) + 0.1_bp*elements(4:6)
          elements(4:6) =  p_matrix(i,1:3)
          CALL NEW(orb_arr(i), elements, "cartesian", frame, t0, center=center_)
       END DO
       DO i=2,4
          CALL setParameters(orb_arr(i), &
               dyn_model=this%dyn_model_prm, &
               perturbers=this%perturbers_prm, &
               integrator=this%integrator_prm, &
               integration_step=this%integration_step_prm)
          y_vector(i) = distance_2b(orb_arr(i))
          IF (error) THEN
             CALL errorMessage("Orbit / new", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
       END DO
       CALL amoeba_private_2b
       IF (PRESENT(iter)) THEN
          iter = iter_
       END IF
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
       ilo = iminloc(y_vector(:))
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,4(F20.15,1X))") "Velocity at " // &
               "epoch 0 and resulting distance from target " // &
               "position:", p_matrix(ilo,:), y_vector(ilo)
          WRITE(stdout,"(2X,A,1X,I0)") "Required number of " // &
               "iterations:", iter_
          WRITE(stdout,"(1X)")
       END IF
       this%elements(1:3) = elements(1:3)
       this%elements(4:6) = p_matrix(ilo,1:3)
       this%dyn_model_prm = "2-body"

       !write(*,*) "2-BODY AMOEBA ELEMENTS:", this%elements, center_, this%center

    CASE ("n-body amoeba")

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

       this%dyn_model_prm = "n-body"
       err_verb_tmp = err_verb
       err_verb = 0
       CALL NEW(orb, ccoord0, ccoord1, "continued fraction", smamax, center=center_)
       IF (error) THEN
          error = .FALSE.
          CALL NULLIFY(orb)
          CALL NEW(orb, ccoord0, ccoord1, "p-iteration", smamax, center=center_)          
       END IF
       err_verb = err_verb_tmp
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "Could not solve the 2-point boundary value " // & 
               "problem in the 2-body approximation.", 1)
          RETURN
       END IF
       IF (PRESENT(ftol)) THEN
          ftol_ = ftol
       ELSE
          ftol_ = ftol_prm
       END IF
       orb_arr(1) = copy(orb)
       CALL setParameters(orb_arr(1), &
            dyn_model="n-body", &
            perturbers=this%perturbers_prm, &
            integrator=this%integrator_prm, &
            integration_step=this%integration_step_prm)
       y_vector(1) = distance(orb_arr(1))
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
       IF (y_vector(1) < ftol_) THEN
          this = copy(orb)
          this%dyn_model_prm = "n-body"
          this%perturbers_prm = perturbers
          this%integrator_prm = integrator
          this%integration_step_prm = integration_step
          CALL NULLIFY(orb)
          IF (PRESENT(iter)) THEN
             iter = 0
          END IF
          RETURN
       END IF
       elements = getElements(orb, "cartesian", frame)
       p_matrix(1,1:3) = elements(4:6)
       CALL NULLIFY(orb)
       DO i=2,4
          CALL randomNumber(ran)
          p_matrix(i,1:3) = elements(4:6) + 0.01_bp*(2.0_bp*ran-1.0_bp)*elements(4:6)
          elements(4:6) =  p_matrix(i,1:3)
          CALL NEW(orb_arr(i), elements, "cartesian", frame, t0, center=center_)
       END DO
       DO i=2,4
          CALL setParameters(orb_arr(i), &
               dyn_model="n-body", &
               perturbers=this%perturbers_prm, &
               integrator=this%integrator_prm, &
               integration_step=this%integration_step_prm)
          y_vector(i) = distance(orb_arr(i))
          IF (error) THEN
             CALL errorMessage("Orbit / new", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
       END DO
       CALL amoeba_private
       IF (PRESENT(iter)) THEN
          iter = iter_
       END IF
       IF (error) THEN
          CALL errorMessage("Orbit / new", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
       ilo = iminloc(y_vector(:))
       IF (info_verb >= 3) THEN
          WRITE(stdout,"(2X,A,4(F20.15,1X))") "Velocity at " // &
               "epoch 0 and resulting distance from target " // &
               "position:", p_matrix(ilo,:), y_vector(ilo)
          WRITE(stdout,"(2X,A,1X,I0)") "Required number of " // &
               "iterations:", iter_
          WRITE(stdout,"(1X)")
       END IF
       this%elements(1:3) = elements(1:3)
       this%elements(4:6) = p_matrix(ilo,1:3)
       this%dyn_model_prm = "n-body"
       !write(*,*) "N-BODY AMOEBA ELEMENTS:", this%elements, center_, this%center

    CASE default 

       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Method chosen for solving the 2-point boundary value " // &
            "problem not available:" // TRIM(method) // ".", 1)
       RETURN

    END SELECT

    CALL NULLIFY(ccoord_)
    NULLIFY(this%additional_perturbers)


  CONTAINS 


    SUBROUTINE amoeba_private_2b 

      IMPLICIT NONE
      TYPE (Orbit) :: orb
      REAL(bp) :: ysave, ytry, ytmp
      INTEGER(ibp) :: i, ilo, inhi 

      ndim = SIZE(p_matrix,dim=2)
      IF (ndim /= SIZE(p_matrix,dim=1) - 1 .OR. ndim /= SIZE(y_vector) - 1) THEN
         error = .TRUE.
         CALL errorMessage("Orbit / new", &
              "Matrix and vector dimensions incompatible.", 1)
         RETURN
      END IF
      iter_ = 0 
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
         IF (info_verb >= 3) THEN
            DO i=1,SIZE(y_vector,dim=1)
               IF (i == ilo) THEN
                  WRITE(stdout,"('ILO:',1X,4(F20.15,1X))") p_matrix(i,:), y_vector(i) 
               ELSE IF (i == ihi) THEN
                  WRITE(stdout,"('IHI:',1X,4(F20.15,1X))") p_matrix(i,:), y_vector(i) 
               ELSE
                  WRITE(stdout,"(5X,4(F20.15,1X))") p_matrix(i,:), y_vector(i) 
               END IF
            END DO
            WRITE(stdout,"(1X)")
         END IF
         ! Compute the fractional range from highest to lowest and
         ! return if satisfactory.
         IF (y_vector(ilo) < ftol_) THEN 
            ! If returning, put best point and value in slot 1. 
            CALL swap(y_vector(1),y_vector(ilo))
            CALL swap(p_matrix(1,:),p_matrix(ilo,:)) 
            RETURN 
         END IF
         IF (iter_ >= ITMAX_PRM) THEN
            ! TMAX exceeded in amoeba
            error = .TRUE.
            CALL errorMessage("Orbit / new", &
                 "Maximum number of iterations exceeded.", 1)
            RETURN
         END IF
         ! Begin a new iteration. First extrapolate by a factor -1
         ! through the face of the simplex across from the high
         ! point, i.e., reflect the simplex from the high point.
         ytry = amotry_2b(-1.0_bp) 
         IF (error) THEN
            CALL errorMessage("Orbit / new", &
                 "TRACE BACK (30)", 1)
            RETURN
         END IF
         iter_ = iter_ + 1 
         IF (ytry <= y_vector(ilo)) THEN
            ! Gives a result better than the best point, so try an
            ! additional extrapolation by a factor of 2.  
            ytry = amotry_2b(2.0_bp)
            IF (error) THEN
               CALL errorMessage("Orbit / new", &
                    "TRACE BACK (35)", 1)
               RETURN
            END IF
            iter_ = iter_ + 1 
         ELSE IF (ytry >= y_vector(inhi)) THEN 
            ! The reflected point is worse than the second highest,
            ! so look for an intermediate lower point, i.e., do a
            ! one-dimensional contraction.
            ysave = y_vector(ihi) 
            ytry = amotry_2b(0.5_bp)
            IF (error) THEN
               CALL errorMessage("Orbit / new", &
                    "TRACE BACK (40)", 1)
               RETURN
            END IF
            iter_ = iter_ + 1
            IF (ytry >= ysave) THEN
               ! Can't seem to get rid of that high point. Better
               ! contract around the lowest (best) point.
               p_matrix(:,:) = 0.5_bp*(p_matrix(:,:)+SPREAD(p_matrix(ilo,:),1,SIZE(p_matrix,1)))
               DO i=1,ndim+1
                  IF (i /= ilo) THEN
                     elements(4:6) = p_matrix(i,:)
                     CALL NEW(orb, elements, "cartesian", frame, t0, center=center_)
                     IF (error) THEN
                        CALL errorMessage("Orbit / new", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     CALL setParameters(orb, &
                          dyn_model=this%dyn_model_prm, &
                          perturbers=this%perturbers_prm, &
                          integrator=this%integrator_prm, &
                          integration_step=this%integration_step_prm)
                     IF (error) THEN
                        CALL errorMessage("Orbit / new", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     y_vector(i) = distance_2b(orb)
                     IF (error) THEN
                        CALL errorMessage("Orbit / new", &
                             "TRACE BACK (50)", 1)
                        RETURN
                     END IF
                     CALL NULLIFY(orb)
                  END IF
               END DO
               iter_ = iter_ + ndim ! Keep track of function evaluations.
               psum(:) = SUM(p_matrix(:,:),dim=1)
            END IF
         END IF
      END DO ! Go back for the test of doneness and the next iteration.

    END SUBROUTINE amoeba_private_2b



    SUBROUTINE amoeba_private 

      IMPLICIT NONE
      TYPE (Orbit) :: orb
      REAL(bp) :: ysave, ytry, ytmp
      INTEGER(ibp) :: i, ilo, inhi 

      ndim = SIZE(p_matrix,dim=2)
      IF (ndim /= SIZE(p_matrix,dim=1) - 1 .OR. ndim /= SIZE(y_vector) - 1) THEN
         error = .TRUE.
         CALL errorMessage("Orbit / new", &
              "Matrix and vector dimensions incompatible.", 1)
         RETURN
      END IF
      iter_ = 0 
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
         IF (info_verb >= 3) THEN
            DO i=1,SIZE(y_vector,dim=1)
               IF (i == ilo) THEN
                  WRITE(stdout,"('ILO:',1X,4(F20.15,1X))") p_matrix(i,:), y_vector(i) 
               ELSE IF (i == ihi) THEN
                  WRITE(stdout,"('IHI:',1X,4(F20.15,1X))") p_matrix(i,:), y_vector(i) 
               ELSE
                  WRITE(stdout,"(5X,4(F20.15,1X))") p_matrix(i,:), y_vector(i) 
               END IF
            END DO
            WRITE(stdout,"(1X)")
         END IF
         ! Compute the fractional range from highest to lowest and
         ! return if satisfactory.
         IF (y_vector(ilo) < ftol_) THEN 
            ! If returning, put best point and value in slot 1. 
            CALL swap(y_vector(1),y_vector(ilo))
            CALL swap(p_matrix(1,:),p_matrix(ilo,:)) 
            RETURN 
         END IF
         IF (iter_ >= ITMAX_PRM) THEN
            ! TMAX exceeded in amoeba
            error = .TRUE.
            CALL errorMessage("Orbit / new", &
                 "Maximum number of iterations exceeded.", 1)
            RETURN
         END IF
         ! Begin a new iteration. First extrapolate by a factor -1
         ! through the face of the simplex across from the high
         ! point, i.e., reflect the simplex from the high point.
         ytry = amotry(-1.0_bp) 
         IF (error) THEN
            CALL errorMessage("Orbit / new", &
                 "TRACE BACK (30)", 1)
            RETURN
         END IF
         iter_ = iter_ + 1 
         IF (ytry <= y_vector(ilo)) THEN
            ! Gives a result better than the best point, so try an
            ! additional extrapolation by a factor of 2.  
            ytry = amotry(2.0_bp)
            IF (error) THEN
               CALL errorMessage("Orbit / new", &
                    "TRACE BACK (35)", 1)
               RETURN
            END IF
            iter_ = iter_ + 1 
         ELSE IF (ytry >= y_vector(inhi)) THEN 
            ! The reflected point is worse than the second highest,
            ! so look for an intermediate lower point, i.e., do a
            ! one-dimensional contraction.
            ysave = y_vector(ihi) 
            ytry = amotry(0.5_bp)
            IF (error) THEN
               CALL errorMessage("Orbit / new", &
                    "TRACE BACK (40)", 1)
               RETURN
            END IF
            iter_ = iter_ + 1
            IF (ytry >= ysave) THEN
               ! Can't seem to get rid of that high point. Better
               ! contract around the lowest (best) point.
               p_matrix(:,:) = 0.5_bp*(p_matrix(:,:)+SPREAD(p_matrix(ilo,:),1,SIZE(p_matrix,1)))
               DO i=1,ndim+1
                  IF (i /= ilo) THEN
                     elements(4:6) = p_matrix(i,:)
                     CALL NEW(orb, elements, "cartesian", frame, t0, center=center_)
                     IF (error) THEN
                        CALL errorMessage("Orbit / new", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     CALL setParameters(orb, &
                          dyn_model=this%dyn_model_prm, &
                          perturbers=this%perturbers_prm, &
                          integrator=this%integrator_prm, &
                          integration_step=this%integration_step_prm)
                     IF (error) THEN
                        CALL errorMessage("Orbit / new", &
                             "TRACE BACK (45)", 1)
                        RETURN
                     END IF
                     y_vector(i) = distance(orb)
                     IF (error) THEN
                        CALL errorMessage("Orbit / new", &
                             "TRACE BACK (50)", 1)
                        RETURN
                     END IF
                     CALL NULLIFY(orb)
                  END IF
               END DO
               iter_ = iter_ + ndim ! Keep track of function evaluations.
               psum(:) = SUM(p_matrix(:,:),dim=1)
            END IF
         END IF
      END DO ! Go back for the test of doneness and the next iteration.

    END SUBROUTINE amoeba_private



    !! *Description*:
    !!
    !! Extrapolates by a factor fac through the face of the
    !! simplex across from the high point, tries it, and replaces
    !! the high point if the new point is better.
    !!
    REAL(bp) FUNCTION amotry_2b(fac)

      IMPLICIT NONE
      REAL(bp), INTENT(IN) :: fac

      TYPE (Orbit) :: orb
      REAL(bp), DIMENSION(SIZE(p_matrix,2)) :: ptry 
      REAL(bp) :: fac1, fac2, ytry

      fac1 = (1.0_bp-fac)/ndim
      fac2 = fac1 - fac 
      ptry(:) = psum(:)*fac1-p_matrix(ihi,:)*fac2 
      ! Evaluate the function at the trial point.
      elements(4:6) = ptry
      CALL NEW(orb, elements, "cartesian", frame, t0, center=center_)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (55)", 1)
         RETURN
      END IF
      CALL setParameters(orb, &
           dyn_model=this%dyn_model_prm, &
           perturbers=this%perturbers_prm, &
           integrator=this%integrator_prm, &
           integration_step=this%integration_step_prm)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (45)", 1)
         RETURN
      END IF
      ytry = distance_2b(orb)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (60)", 1)
         RETURN
      END IF
      CALL NULLIFY(orb)
      IF (ytry < y_vector(ihi)) THEN 
         ! If it's better than the highest, then replace
         ! the highest.
         y_vector(ihi) = ytry
         psum(:) = psum(:) - p_matrix(ihi,:) + ptry(:)
         p_matrix(ihi,:) = ptry(:)
      END IF
      amotry_2b = ytry 

    END FUNCTION amotry_2b



    !! *Description*:
    !!
    !! Extrapolates by a factor fac through the face of the
    !! simplex across from the high point, tries it, and replaces
    !! the high point if the new point is better.
    !!
    REAL(bp) FUNCTION amotry(fac)

      IMPLICIT NONE
      REAL(bp), INTENT(IN) :: fac

      TYPE (Orbit) :: orb
      REAL(bp), DIMENSION(SIZE(p_matrix,2)) :: ptry 
      REAL(bp) :: fac1, fac2, ytry

      fac1 = (1.0_bp-fac)/ndim
      fac2 = fac1 - fac 
      ptry(:) = psum(:)*fac1-p_matrix(ihi,:)*fac2 
      ! Evaluate the function at the trial point.
      elements(4:6) = ptry
      CALL NEW(orb, elements, "cartesian", frame, t0, center=center_)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (55)", 1)
         RETURN
      END IF
      CALL setParameters(orb, &
           dyn_model=this%dyn_model_prm, &
           perturbers=this%perturbers_prm, &
           integrator=this%integrator_prm, &
           integration_step=this%integration_step_prm)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (45)", 1)
         RETURN
      END IF
      ytry = distance(orb)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (60)", 1)
         RETURN
      END IF
      CALL NULLIFY(orb)
      IF (ytry < y_vector(ihi)) THEN 
         ! If it's better than the highest, then replace
         ! the highest.
         y_vector(ihi) = ytry
         psum(:) = psum(:) - p_matrix(ihi,:) + ptry(:)
         p_matrix(ihi,:) = ptry(:)
      END IF
      amotry = ytry 

    END FUNCTION amotry



    REAL(bp) FUNCTION distance_2b(orb)

      IMPLICIT NONE
      TYPE (Orbit), INTENT(inout) :: orb 

      REAL(bp), DIMENSION(3) :: pos_orb

      CALL propagate(orb, t1)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (65)", 1)
         RETURN
      END IF
      pos_orb = getPosition(orb)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (75)", 1)
         RETURN
      END IF
      distance_2b = SQRT(SUM((pos1-pos_orb)**2))

    END FUNCTION distance_2b




    REAL(bp) FUNCTION distance(orb)

      IMPLICIT NONE
      TYPE (Orbit), INTENT(inout) :: orb 

      REAL(bp), DIMENSION(3) :: pos_orb1

      CALL propagate(orb, t1)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (65)", 1)
         RETURN
      END IF
      pos_orb1 = getPosition(orb)
      IF (error) THEN
         CALL errorMessage("Orbit / new", &
              "TRACE BACK (75)", 1)
         RETURN
      END IF
      distance = SQRT(SUM((pos1-pos_orb1)**2))

    END FUNCTION distance

  END SUBROUTINE new_Orb_2point





  !! *Description*:
  !!
  !! Initializes a Orbit-object using a given
  !! SphericalCoordinates-object. Coordinate frame s equal to the
  !! frame of the input frame and propagation scheme is 2-body.
  !! Central body is the Sun.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_Orb_spherical(this, scoord)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout) :: this
    TYPE (SphericalCoordinates), INTENT(in) :: scoord
    TYPE (CartesianCoordinates) :: ccoord

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    ! Transform to Cartesian coordinates:
    CALL NEW(ccoord, scoord)
    IF (error) THEN
       CALL errorMessage("Orbit / new",&
            "TRACE BACK 1", 1)
       RETURN
    END IF
    this%elements = getCoordinates(ccoord)
    IF (error) THEN
       CALL errorMessage("Orbit / new",&
            "TRACE BACK 2", 1)
       RETURN
    END IF
    this%element_type = "cartesian"
    this%frame = getFrame(ccoord)
    this%t = getTime(ccoord)
    CALL NULLIFY(ccoord)
    NULLIFY(this%additional_perturbers)
    this%dyn_model_prm = "2-body"
    this%finite_diff_prm = -1.0_bp
    this%is_initialized = .TRUE.
    this%center = 11

  END SUBROUTINE new_Orb_spherical





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout) :: this
    INTEGER :: err

    this%elements = 0.0_bp
    this%element_type = ""
    CALL NULLIFY(this%t)
    this%finite_diff_prm = -1.0_bp
    this%frame = ""
    this%mass_prm = -1.0_bp
    this%dyn_model_prm    = ""
    this%perturbers_prm = .FALSE.
    IF (ASSOCIATED(this%additional_perturbers)) THEN
       DEALLOCATE(this%additional_perturbers, stat=err)
    ELSE
       NULLIFY(this%additional_perturbers)
    END IF
    this%is_initialized = .FALSE.
    this%center   = 11

  END SUBROUTINE nullify_Orb





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  !! Returns error.
  !!
  FUNCTION copy_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    TYPE (Orbit)             :: copy_Orb
    INTEGER :: err

    copy_Orb%elements       = this%elements
    copy_Orb%element_type   = this%element_type
    copy_Orb%frame          = this%frame
    copy_Orb%mass_prm           = this%mass_prm
    copy_Orb%t              = copy(this%t)
    copy_Orb%center   = this%center
    copy_Orb%dyn_model_prm  = this%dyn_model_prm
    copy_Orb%perturbers_prm(:) = this%perturbers_prm(:)
    IF (ASSOCIATED(this%additional_perturbers)) THEN
       ALLOCATE(copy_Orb%additional_perturbers(SIZE(this%additional_perturbers,dim=1), &
            SIZE(this%additional_perturbers,dim=2)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / copy", &
               "Could not allocate memory.",1)
          RETURN
       END IF
       copy_Orb%additional_perturbers = this%additional_perturbers
    ELSE
       NULLIFY(copy_Orb%additional_perturbers)
    END IF
    copy_Orb%finite_diff_prm = this%finite_diff_prm
    copy_Orb%integrator_prm = this%integrator_prm
    copy_Orb%integration_step_prm = this%integration_step_prm
    copy_Orb%is_initialized = this%is_initialized

  END FUNCTION copy_Orb





  !! *Description*:
  !!
  !! Returns the status of the object, i.e. whether it exists or not.
  !!
  LOGICAL FUNCTION exist_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this

    exist_Orb = this%is_initialized

  END FUNCTION exist_Orb





  REAL(bp) FUNCTION aeidist(ptry,param,errstr)

    IMPLICIT NONE
    REAL(bp), DIMENSION(:), INTENT(in) :: ptry
    REAL(bp), DIMENSION(:), INTENT(in) :: param
    CHARACTER(len=*), INTENT(inout) :: errstr

    TYPE (Orbit) :: orb
    TYPE (Time) :: t
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp), DIMENSION(3) :: vel

    vel = ptry 
    elements = (/ param(1:3), vel /)
    CALL NEW(t, param(7), "TT")
    IF (error) THEN
       errstr = "error12"
       RETURN
    END IF
    CALL NEW(orb, elements, "cartesian", "ecliptic", t)
    IF (error) THEN
       errstr = "error13"
       RETURN
    END IF
    elements = getCometaryElements(orb, qei_only=.TRUE.)
    IF (error) THEN
       errstr = "error14"
       RETURN
    END IF
    CALL NULLIFY(orb)
    CALL NULLIFY(t)
    elements(1) = elements(1)/(1.0_bp-elements(2))
    aeidist = SQRT(SUM((ABS(elements(1:3)-param(4:6))/param(4:6))**2))

  END FUNCTION aeidist





  REAL(bp) FUNCTION xyzdist(ptry,param,errstr)

    IMPLICIT NONE
    REAL(bp), DIMENSION(:), INTENT(in) :: ptry
    REAL(bp), DIMENSION(:), INTENT(in) :: param
    CHARACTER(len=*), INTENT(inout) :: errstr

    TYPE (Orbit) :: orb
    TYPE (Time) :: t
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp), DIMENSION(3) :: nam

    nam = ptry 
    elements = (/ param(1:3), nam /)
    CALL NEW(t, param(7), "TT")
    IF (error) THEN
       errstr = "error12"
       RETURN
    END IF
    CALL NEW(orb, elements, "keplerian", "ecliptic", t)
    IF (error) THEN
       errstr = "error13"
       RETURN
    END IF
    elements = getCartesianElements(orb, "ecliptic")
    IF (error) THEN
       errstr = "error14"
       RETURN
    END IF
    CALL NULLIFY(orb)
    CALL NULLIFY(t)
    xyzdist = SQRT(SUM((elements(1:3)-param(4:6))**2))

  END FUNCTION xyzdist





  !! *Description*:
  !!
  !! Checks the orbit, bound/unbound.
  !!
  !! Returns error.
  !!
  LOGICAL FUNCTION boundOrbit(this, smamax, sma)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)     :: this
    REAL(bp), INTENT(in)            :: smamax

    REAL(bp), OPTIONAL, INTENT(out) :: sma
    REAL(bp), DIMENSION(6)          :: elements
    REAL(bp)                        :: a, r, alpha, tmp

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / boundOrbit", &
            "Object has not been initialized.", 1)
       RETURN
    END IF

    elements = getCartesianElements(this,"equatorial")
    IF (error) THEN
       CALL errorMessage("Orbit / boundOrbit", &
            "TRACE BACK 1",1)
       RETURN
    END IF

    r = SQRT(DOT_PRODUCT(elements(1:3),elements(1:3)))
    ! alpha = rv^2/mu
    alpha = r*DOT_PRODUCT(elements(4:6),elements(4:6))/planetary_mu(this%center)
    ! h = v^2/2 - mu/r and for elliptical orbits a = -mu/2h
    ! -> a = r / (2 - rv^2/mu)
    tmp = 2.0_bp - alpha
    IF (ABS(tmp) < 100.0_bp*EPSILON(tmp)) THEN
       boundOrbit = .FALSE.
       RETURN
    END IF
    a = r / tmp

    IF (PRESENT(sma)) THEN
       sma = a
    END IF

    ! Make sure the semimajor axis is within the boundary values.
    ! Lower bound:
    IF (a < planetary_radii(this%center)) THEN
       boundOrbit = .FALSE.
       RETURN
    END IF
    ! Upper bound:
    IF (a > smamax) THEN
       boundOrbit = .FALSE.
       RETURN
    END IF

    boundOrbit = .TRUE.

  END FUNCTION boundOrbit





  !! *Description*:
  !!
  !! The continued fraction for computing the orbit parameter
  !! Comments: Based on f77 routine by Karri Muinonen,
  !!           partly based on software written by Olof Hernius.  
  !! 
  !! Reference: "The determination of orbits" (1961), A. D. Dubyago.
  !!
  !! Returns error.
  !!
  SUBROUTINE continuedFraction(h, tol, sector_to_triangel_ratio, &
       dsector_to_triangel_ratio)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: h, tol
    REAL(bp), INTENT(out) :: sector_to_triangel_ratio, &
         dsector_to_triangel_ratio
    INTEGER, PARAMETER :: niter = 10000
    REAL (bp) :: tmp1, tmp2, d1, d2
    INTEGER :: i

    d1 = 10.0_bp/11.0_bp
    d2 = 11.0_bp/9.0_bp
    tmp1 = d2*h/(1.0_bp+d2*h)
    tmp2 = tmp1
    DO i=1,niter
       tmp1 = d2*h/(1.0_bp + tmp1)
       dsector_to_triangel_ratio = ABS(tmp1-tmp2)
       IF (dsector_to_triangel_ratio < tol) THEN
          sector_to_triangel_ratio = 1.0_bp + d1*tmp1
          RETURN
       END IF
       tmp2 = tmp1
    END DO

  END SUBROUTINE continuedFraction





  !! *Description*:
  !!
  !! Computes an estimate for the cosine of the true anomaly
  !! difference.
  !! 
  !! Returns error.
  !!
  REAL(bp) FUNCTION estimateCosDf(this, ccoord, p, y)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                :: this
    TYPE (CartesianCoordinates), INTENT(in) :: ccoord
    REAL(bp), INTENT(in)                    :: p
    INTEGER, INTENT(in)                     :: y
    TYPE (Orbit)                            :: orb
    TYPE (Time)                             :: t1
    REAL(bp), DIMENSION(3)                  :: pos0, pos1, pos, vel0
    REAL(bp)                                :: r0, r1, cdf, sdf, f, g

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / estimateCosDf", &
            "This object has not yet been initialized.", 1)
       RETURN       
    END IF

    IF (this%element_type /= "cartesian") THEN
       error = .TRUE.
       CALL errorMessage("Orbit / estimateCosDf", &
            "Orbital elements must be Cartesian.", 1)
       RETURN
    END IF

    IF (.NOT. exist(ccoord)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / estimateCosDf", &
            "Second position has not been initialized.", 1)
       RETURN       
    END IF

    pos0 = this%elements(1:3)
    pos1 = getPosition(ccoord)

    r0  = SQRT(DOT_PRODUCT(pos0,pos0))
    r1  = SQRT(DOT_PRODUCT(pos1,pos1))
    cdf = DOT_PRODUCT(pos0,pos1) / (r0*r1)
    sdf = y * SQRT(ABS(1.0_bp - cdf**2))

    f = (r1/p)*(cdf - 1.0_bp) + 1.0_bp
    g = ((r0*r1) / SQRT(planetary_mu(this%center)*p)) * sdf
    IF (ABS(g) < TINY(g)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / estimateCosDf", &
            "Attempted division by zero (g-function).", 1)
       RETURN
    END IF
    vel0 = (pos1 - f*pos0)/g
    orb = copy(this)
    orb%elements(4:6) = vel0
    t1 = getTime(ccoord)
    CALL propagate(orb, t1)
    IF (error) THEN
       CALL errorMessage("Orbit / estimateCosDf", &
            "TRACE BACK 4", 1)
       RETURN
    END IF
    pos = orb%elements(1:3)
    estimateCosDf = DOT_PRODUCT(pos0,pos) / (r0*SQRT(DOT_PRODUCT(pos,pos)))

    CALL NULLIFY(orb)

  END FUNCTION estimateCosDf





  !! *Description*:
  !!
  !! Calculates the Jacobian matrix of final coordinates wrt initial
  !! coordinates when Gauss' f- and g-functions are used.
  !!
  SUBROUTINE GaussfgJacobian_Orb(this, r0, u, alpha, stumpff_cs, s, &
       f, g, df, dg, pos, r1, jacobian)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)              :: this
    REAL(bp), INTENT(in)                  :: r0, u, alpha, s, f, g, df, dg, r1
    REAL(bp), DIMENSION(0:3), INTENT(in)  :: stumpff_cs
    REAL(bp), DIMENSION(3), INTENT(in)    :: pos
    REAL(bp), DIMENSION(6,6), INTENT(out) :: jacobian
    REAL(bp), PARAMETER                   :: tol = 1.0e-20_bp
    REAL(bp), DIMENSION(6,3)              :: tmp_S, pstumpff_cs
    REAL(bp), DIMENSION(6)                :: pr0, pdr0, palpha, pr1, tmp_A, &
         tmp_B, tmp_C, ps, tmp_e, pf, pg, pdf, pdg
    REAL(bp)                              :: dr0, tmp_d, mu_
    INTEGER                               :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / GaussfgJacobian", &
            "This object has not yet been initialized.", 1)
       RETURN       
    END IF

    IF (ABS(r0) < tol .OR. ABS(alpha) < tol .OR. &
         ABS(s) < tol .OR. ABS(r1) < tol) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / GaussfgJacobian", &
            "Preliminary criterions are not fulfilled.", 1)
       RETURN
    END IF

    ! Define mu parameter
    mu_ = planetary_mu(this%center)

    ! Partial derivatives of r0, dr0, and alpha.
    pr0(1:3)    = this%elements(1:3)/r0
    pr0(4:6)    = 0
    dr0         = u/r0
    pdr0(1:3)   =  -this%elements(1:3)*dr0/r0**2 + this%elements(4:6)/r0
    pdr0(4:6)   = this%elements(1:3)/r0
    palpha(1:3) = -2.0_bp*mu_/r0**3.0_bp*this%elements(1:3) 
    palpha(4:6) = -2.0_bp*this%elements(4:6)

    ! Partial derivatives of s.
    tmp_A = pr0*stumpff_cs(1) + (pr0*dr0+r0*pdr0)*stumpff_cs(2)
    tmp_B = (r0-mu_/alpha)*(stumpff_cs(1) - s*stumpff_cs(0)) + &
         (2.0_bp*stumpff_cs(2)-s*stumpff_cs(1))*u + 2.0_bp*mu_*stumpff_cs(3)
    tmp_C = (r0*stumpff_cs(1) + 2.0_bp*u*stumpff_cs(2) + &
         3.0_bp*mu_*stumpff_cs(3) - tmp_B)/s
    IF (ANY(ABS(tmp_C) < tol)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / GaussfgJacobian", &
            "Criterions are not fulfilled.", 1)
       RETURN
    END IF
    ps = (0.5_bp*palpha/alpha*tmp_B - tmp_A)/tmp_C

    ! Partial derivatives of Stumpff-functions 1-3 (0 not needed).
    tmp_d      = 1/(alpha*s**2.0_bp)
    tmp_e      = palpha*s**2.0_bp + 2.0_bp*s*alpha*ps
    tmp_S(:,1) = 0.5_bp*tmp_d*(s*stumpff_cs(0)-stumpff_cs(1))*tmp_e
    tmp_S(:,2) = tmp_d*(0.5_bp*s*stumpff_cs(1)-stumpff_cs(2))*tmp_e
    tmp_S(:,3) = tmp_d*(0.5_bp*tmp_d*s**2_bp*(stumpff_cs(1)-s*stumpff_cs(0)) &
         -stumpff_cs(3))*tmp_e

    pstumpff_cs(:,1) = ps*stumpff_cs(1)/s+tmp_S(:,1)
    pstumpff_cs(:,2) = 2.0_bp*ps*stumpff_cs(2)/s+tmp_S(:,2)
    pstumpff_cs(:,3) = 3.0_bp*ps*stumpff_cs(3)/s+tmp_S(:,3)

    ! Partial derivatives of Gauss" f- and g-functions
    ! and the Jacobian.
    pf  = mu_/r0*(pr0/r0*stumpff_cs(2)-pstumpff_cs(:,2))
    pg  = -mu_*pstumpff_cs(:,3)

    DO i=1,3
       jacobian(i,:)   = pf*this%elements(i) + pg*this%elements(3+i)
       jacobian(i,i)   = jacobian(i,i) + f
       jacobian(i,i+3) = jacobian(i,i+3) + g
    END DO
    DO i=1,6
       pr1(i) = DOT_PRODUCT(pos,jacobian(1:3,i))/r1
    END DO
    pdf = mu_/(r0*r1)*((pr0/r0+pr1/r1)*stumpff_cs(1) - pstumpff_cs(:,1))
    pdg = mu_/r1*(pr1/r1*stumpff_cs(2) - pstumpff_cs(:,2))
    DO i=4,6
       jacobian(i,:)   = pdf*this%elements(i-3) + pdg*this%elements(i)
       jacobian(i,i-3) = jacobian(i,i-3) + df
       jacobian(i,i)   = jacobian(i,i) + dg
    END DO

  END SUBROUTINE GaussfgJacobian_Orb





  !! *Description*:
  !!
  !! Returns the apoapsis distance Q for the input orbit. The
  !! semimajor axis 'a' is computed if it is not given as a
  !! parameter. Optionally, returns partial derivatives between Q and
  !! Keplerian elements.
  !!
  !! Returns error.
  !!
  SUBROUTINE getApoapsisDistance_Orb(this, Q, a, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                      :: this
    REAL(bp), INTENT(out)                         :: Q
    REAL(bp), INTENT(in), OPTIONAL                :: a
    REAL(bp), DIMENSION(6), INTENT(out), OPTIONAL :: partials ! Partials

    TYPE (Orbit) :: this_
    REAL(bp), DIMENSION(3) :: pos, vel
    REAL(bp) :: a_, e, r, e_sin_ea, e_cos_ea, alpha, gamma

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getApoapsisDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    SELECT CASE (this%element_type)

    CASE ("cometary")

       IF (this%elements(2) < 1.0_bp) THEN
          ! Elliptic orbit
          Q = this%elements(1) * &
               (1.0_bp + this%elements(2)) / &
               (1.0_bp - this%elements(2))
       ELSE
          ! Parabolic or hyperbolic orbit
          Q = HUGE(Q)
       END IF
       IF (PRESENT(partials)) THEN
          partials = 0.0_bp
          partials(1:2) = (/ 1.0_bp + this%elements(2), this%elements(1) / &
               (1.0_bp - this%elements(2)) /)
       END IF

    CASE ("keplerian")

       Q = this%elements(1)*(1.0_bp+this%elements(2))
       IF (PRESENT(partials)) THEN
          partials = 0.0_bp
          partials(1:2) = (/ 1.0_bp + this%elements(2), this%elements(1) /)
       END IF

    CASE ("cartesian")

       this_ = copy(this)
       CALL rotateToEcliptic(this_)
       pos = this_%elements(1:3)
       vel = this_%elements(4:6)
       r = SQRT(DOT_PRODUCT(pos,pos))
       ! alpha = rv^2/mu
       alpha = r*DOT_PRODUCT(vel,vel)/planetary_mu(this_%center)
       IF (PRESENT(a)) THEN
          a_ = a
       ELSE
          ! h = v^2/2 - mu/r and for elliptical orbits a = -mu/2h
          ! -> a = r / (2 - rv^2/mu)
          IF (ABS(2.0_bp - alpha) < 10.0_bp*EPSILON(alpha)) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getApoapsisDistance", &
                  "a is approaching infinity.", 1)
             RETURN
          END IF
          a_ = r / (2.0_bp - alpha)
          IF (a_ < 0.0_bp) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getApoapsisDistance", &
                  "Semimajor axis is negative.", 1)
             RETURN
          END IF
       END IF
       ! gamma = sqrt(mu*a)
       gamma = SQRT(planetary_mu(this_%center)*a_)
       IF (ABS(gamma) < TINY(gamma)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "Gamma is computationally too small.", 1)
          RETURN
       END IF
       ! r(1+e*cos(f))=a|1-e^2|, r*cos(f)=a*cos(E)-ae, and for elliptical 
       ! orbits a=r/(2-rv^2/mu) (see above) and |1-e^2|=1-e^2
       ! -> e*cos(E)=rv^2/mu-1=alpha-1
       e_cos_ea = alpha - 1.0_bp
       !
       e_sin_ea = DOT_PRODUCT(pos,vel)/gamma
       e = SQRT(e_cos_ea**2.0_bp + e_sin_ea**2.0_bp)
       IF (e > 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "Orbit is hyperbolic.", 1)
          RETURN
       ELSE IF (e < 10.0_bp*EPSILON(e)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getApoapsisDistance", &
               "Orbit is almost circular (1).", 1)
          WRITE(stderr,*) e
          RETURN
       END IF
       CALL NULLIFY(this_)
       Q = a_*(1.0_bp + e)

       IF (PRESENT(partials)) THEN
          partials = 0.0_bp
          partials(1:2) = (/ 1.0_bp + e, a_ /)
       END IF

    END SELECT

  END SUBROUTINE getApoapsisDistance_Orb





  !! *Description*:
  !!
  !! Returns Cartesian orbital elements.
  !!
  !!   - array(1) = x
  !!   - array(2) = y
  !!   - array(3) = z
  !!   - array(4) = dx/dt
  !!   - array(5) = dy/dt
  !!   - array(6) = dz/dt
  !!
  !!
  !! *Usage*:
  !!
  !! elements = getCartesianElements(myorbit)
  !! 
  !! Returns error.
  !!
  FUNCTION getCartesianElements(this, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)     :: this
    CHARACTER(len=*), INTENT(in) :: frame
    REAL(bp), DIMENSION(6)       :: getCartesianElements

    TYPE (Orbit) :: this_
    TYPE (Time) :: t
    CHARACTER(len=FRAME_LEN) :: frame_
    REAL(bp), DIMENSION(3,3) :: R
    REAL(bp), DIMENSION(6) :: celements
    REAL(bp) :: sea, cea, ea, dot_ea, b

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getCartesianElements", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    frame_ = frame
    CALL locase(frame_, error)
    IF (error) THEN
       CALL errorMessage("Orbit / getCartesianElements", &
            "The frame string contains forbidden characters.", 1)
       RETURN
    END IF
    IF (frame_ /= "ecliptic" .AND. &
         frame_ /= "equatorial") THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getCartesianElements", &
            "Frame " // TRIM(frame_) // " not recognized.", 1)
    END IF

    IF (this%element_type == "cartesian") THEN
       getCartesianElements(1:6) = this%elements(1:6)
       IF (frame_ == "equatorial" .AND. &
            this%frame /= "equatorial") THEN
          CALL rotateToEquatorial(getCartesianElements(1:6))
       ELSE IF (frame_ == "ecliptic" .AND. &
            this%frame /= "ecliptic") THEN
          CALL rotateToEcliptic(getCartesianElements(1:6))
       END IF
       RETURN
    END IF

    SELECT CASE (this%element_type)

    CASE ("cometary")

       ! Make transformation cometary -> cartesian at periapsis:
       celements(1:3) = (/ this%elements(1), 0.0_bp, 0.0_bp /)
       ! v^2 = mu * (2/r - 1/a) = mu * (1+e)/q
       celements(4:6) = (/ 0.0_bp, &
            SQRT(planetary_mu(this%center) * &
            (1.0_bp+this%elements(2))/this%elements(1)), 0.0_bp /)

       ! Orbital-plane Cartesian elements to ecliptical Cartesian
       ! elements:
       R = getTransformationMatrix(this)
       celements(1:3) = MATMUL(R,celements(1:3))
       celements(4:6) = MATMUL(R,celements(4:6))

       ! Propagate elements (must use 2b!) to the initial epoch
       CALL NEW(t, this%elements(6), "TT")
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       CALL NEW(this_, celements, "cartesian", "ecliptic", t)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL setParameters(this_, dyn_model="2-body")
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
       CALL propagate(this_, this%t)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
       getCartesianElements = getElements(this_, "cartesian", frame_)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (30)", 1)
          RETURN
       END IF
       CALL NULLIFY(this_)
       CALL NULLIFY(t)

    CASE ("cometary_ta", "cometary_ma")

       this_ = copy(this)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (35)", 1)
          RETURN
       END IF
       CALL toCometary(this_)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (40)", 1)
          RETURN
       END IF
       getCartesianElements = getElements(this_, "cartesian", frame_)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       CALL NULLIFY(this_)

    CASE ("keplerian")

       ! Compute needed quantities:
       CALL solveKeplerEquation(this, this%t, ea)
       IF (error) THEN
          CALL errorMessage("Orbit / getCartesianElements", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       cea = COS(ea)
       sea = SIN(ea)
       b = this%elements(1) * SQRT(1.0_bp - this%elements(2)**2.0_bp)
       dot_ea = SQRT(planetary_mu(this%center)/this%elements(1)**3.0_bp) / &
            (1.0_bp - this%elements(2)*cea)
       ! Keplerian elements to polar Cartesian elements:
       ! -positions:
       celements(1) = this%elements(1)*(cea - this%elements(2))
       celements(2) = b*sea
       celements(3) = 0.0_bp
       ! -velocities:
       celements(4) = -this%elements(1)*dot_ea*sea
       celements(5) = b*dot_ea*cea
       celements(6) = 0.0_bp

       ! Polar Cartesian elements to ecliptical Cartesian elements:
       R = getTransformationMatrix(this)
       celements(1:3) = MATMUL(R,celements(1:3))
       celements(4:6) = MATMUL(R,celements(4:6))

       IF (frame_ == "equatorial") THEN
          CALL rotateToEquatorial(celements)
       END IF
       getCartesianElements(1:6) = celements(1:6)

    CASE default

       error = .TRUE.
       CALL errorMessage("Orbit / getCartesianElements", &
            "Conversion from " // TRIM(this%element_type) // &
            " elements to cartesian elements has not yet been implemented.", 1)
       RETURN

    END SELECT

  END FUNCTION getCartesianElements





  !! *Description*:
  !!
  !! Returns cartesian heliocentric coordinates.
  !!
  !! Returns error.
  !!
  FUNCTION getCCoord_Orb(this, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)    :: this
    TYPE (CartesianCoordinates) :: getCCoord_Orb
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame
    CHARACTER(len=FRAME_LEN) :: frame_
    REAL(bp), DIMENSION(6) :: celements

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getCCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(frame)) THEN
       frame_ = frame
       CALL locase(frame_, error)
       IF (error) THEN
          CALL errorMessage("Orbit / getCCoord", &
               "The frame string contains forbidden characters.", 1)
          RETURN
       END IF
    ELSE IF (this%element_type == "cartesian") THEN
       frame_ = this%frame
    ELSE IF (.NOT.PRESENT(frame)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getCCoord", &
            "Frame missing when using Keplerian or cometary input elements.", 1)
       RETURN       
    END IF

    celements(1:6) = getCartesianElements(this, frame=frame_)
    IF (error) THEN
       CALL errorMessage("Orbit / getCCoord", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    CALL NEW(getCCoord_Orb, celements, frame_, copy(this%t))
    IF (error) THEN
       CALL errorMessage("Orbit / getCCoord", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF

  END FUNCTION getCCoord_Orb





  !! *Description*:
  !!
  !! Calculates the cometary orbital elements from the heliocentric
  !! (ecliptical or equatorial) Cartesian orbital elements, Keplerian
  !! elements, or modified cometary elements (true anomaly instead of
  !! time of perihelion) for the current epoch.
  !!
  !! Note that the time of perihelion of elliptic orbits is the
  !! closest one in time.
  !!
  !! Returns error.
  !! 
  FUNCTION getCometaryElements(this, qei_only)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    LOGICAL, INTENT(in), OPTIONAL :: qei_only

    REAL(bp), DIMENSION(6)   :: getCometaryElements
    TYPE (Orbit) :: this_
    !! The following four tolerance parameters have been tuned using
    !! PS S3M populations SI, SL, and SH (Earth impactors, long-period
    !! comets, hyperbolic objects). Nearly-parabolic orbits with
    !! eccentricities having offsets larger than 10^(-13) from unity
    !! should be correctly converted to cometary elements. The
    !! relative error of the conversion com->car->com is, e.g, < 10^-7
    !! for SI whereas for e~1-10^-13 SL orbits the error on q is up to
    !! 0.2%. Note that the latter is only true for e~1-10^-13 orbits,
    !! not for e~1+10^-13 orbits which are a lot more accurate!
    !! REAL(bp), PARAMETER :: tol1 = 1.0e-6_bp
    !! REAL(bp), PARAMETER :: tol2 = 1.0e-8_bp
    !! REAL(bp), PARAMETER :: tol3 = 1.0e-11_bp
    !! REAL(bp), PARAMETER :: tol4 = 1.0e-3_bp
    !! The tolerances have subsequently been updated using the S1b
    !! population and when generating initial orbits for the tidal
    !! disruption studies.
!!$    REAL(bp), PARAMETER :: tol1 = 1.0e-6_bp
!!$    REAL(bp), PARAMETER :: tol2 = 1.0e-2_bp
!!$    REAL(bp), PARAMETER :: tol3 = 1.0e-1_bp
!!$    REAL(bp), PARAMETER :: tol4 = 1.0e-2_bp
    REAL(bp), PARAMETER :: tol1 = 1.0e-6_bp
    REAL(bp), PARAMETER :: tol2 = 1.0e-6_bp
    REAL(bp), PARAMETER :: tol3 = 1.0e-9_bp
    REAL(bp), PARAMETER :: tol4 = 1.0e-3_bp
    INTEGER, PARAMETER :: max_iter = 100
    REAL(bp), DIMENSION(0:3) :: stumpff_c, stumpff_cs
    REAL(bp), DIMENSION(3) :: pos, vel, k, cos_angles, evec, fb, gb
    REAL(bp) :: r0, ru, alpha, a, e, i, an, ap, varpi, tmp1, tmp2, &
         div, q, tp, r, rp, rpp, xv, s, ds, x, dt, cosu, u0, mjd_tt, &
         p, ea, sin_ea, cos_ea, ma, mm
    INTEGER :: iiter

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getCometaryElements", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Return initial elements if error occurs:
    getCometaryElements(1:6) = this%elements(1:6)

    SELECT CASE (this%element_type)

    CASE ("cometary")

       ! Return immediately if elements are cometary
       RETURN

    CASE ("cometary_ta")

       this_ = copy(this)
       IF (this%elements(2) < 1.0_bp) THEN
          ! cos(ea) = (cos(f)+e)/(1+e*cos(f))
          cos_ea = (COS(this%elements(6)) + this%elements(2)) / &
               (1.0_bp + this%elements(2)*COS(this%elements(6)))
          ! sin(ea) = sqrt(1-e^2)*sin(f)/(1+e*cos(f))
          sin_ea = SQRT(1.0_bp - this%elements(2)**2.0_bp)*SIN(this%elements(6)) / &
               (1.0_bp + this%elements(2)*COS(this%elements(6)))
          ea = ATAN2(sin_ea,cos_ea)
          ma = ea - this%elements(2)*SIN(ea)
          ! Time of periapsis:
          mm = SQRT((this_%elements(1) / &
               (1.0_bp-this_%elements(2)))**3.0_bp / &
               planetary_mu(this_%center))
          mjd_tt = getMJD(this_%t,"TT")
          dt = ma * mm
          p = two_pi * mm
          ! Select the closest time of perihelion:
          IF (ma <= pi) THEN
             getCometaryElements(6) = mjd_tt - dt
          ELSE
             getCometaryElements(6) = mjd_tt + (p - dt)
          END IF
       ELSE
          error = .TRUE.
          CALL errorMessage("Orbit / getCometaryElements", &
               "Conversion from ta to tp not yet available for e>=1 orbits.", 1)
          RETURN
       END IF
       CALL NULLIFY(this_)

    CASE ("cometary_ma")

       this_ = copy(this)
       ! Time of periapsis:
       mm = SQRT(ABS(this_%elements(1) / &
            (1.0_bp-this_%elements(2)))**3.0_bp / &
            planetary_mu(this_%center))
       mjd_tt = getMJD(this_%t,"TT")
       dt = this%elements(6) * mm
       p = two_pi * mm
       ! Select the closest time of perihelion:
       IF (this%elements(6) <= pi) THEN
          getCometaryElements(6) = mjd_tt - dt
       ELSE
          getCometaryElements(6) = mjd_tt + (p - dt)
       END IF
       CALL NULLIFY(this_)

    CASE ("keplerian")

       this_ = copy(this)
       ! Distance of periapsis:
       getCometaryElements(1) = &
            this_%elements(1) * (1.0_bp - this_%elements(2))
       ! period p:
       mjd_tt = getMJD(this_%t,"TT")
       ! |a| is there to allow use with hyperbolic orbits
       mm = SQRT(ABS(this_%elements(1))**3.0_bp / planetary_mu(this_%center))
       dt = this_%elements(6) * mm
       p = two_pi * mm
       ! Remove full periods from dt:
       IF (dt > p) THEN
          dt = MODULO(dt,p)
       END IF
       ! Select the closest time of perihelion:
       IF (this_%elements(6) <= pi) THEN
          getCometaryElements(6) = mjd_tt - dt
       ELSE
          getCometaryElements(6) = mjd_tt + (p - dt)
       END IF
       CALL NULLIFY(this_)

    CASE ("cartesian")

       this_ = copy(this)

!!$       call toKeplerian(this_)
!!$       getCometaryElements = getElements(this_, "cometary")
!!$       call nullify(this_)
!!$       if (error) then
!!$          CALL errorMessage("Orbit / getCometaryElements", &
!!$               "TRACE BACK (10)", 1)
!!$          return
!!$       end if
!!$       
!!$
       CALL rotateToEcliptic(this_)

       ! Semimajor axis, eccentricity and mean anomaly:
       pos = this_%elements(1:3)
       r0  = SQRT(DOT_PRODUCT(pos,pos))
       vel = this_%elements(4:6)
       ! h = v^2/2 - mu/r and for elliptical orbits a = -mu/2h
       ! -> a = r / (2 - rv^2/mu)
       ! Semimajor axis (a<0 for e>1 orbits):
       a = r0 / (2.0_bp - r0*DOT_PRODUCT(vel,vel) / &
            planetary_mu(this_%center))
       ! Angular momentum:
       k = cross_product(pos,vel)
       ! Eccentricity vector:
       evec = cross_product(vel,k)/planetary_mu(this_%center) - pos/r0
       ! Eccentricity:
       e = SQRT(DOT_PRODUCT(evec,evec))
       ! Periapsis distance (note that if e>1 then a<0):
       q = a * (1.0_bp - e)

       ! Inclination and ascending node:
       ru = SQRT(DOT_PRODUCT(k,k))
       IF (ru < EPSILON(ru)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getCometaryElements", &   
               "Position and velocity almost parallel.", 1)
          RETURN
       END IF
       k = k/ru
       cos_angles(1) = k(3)
       IF (ABS(cos_angles(1)) > 1.0_bp) THEN
          cos_angles(1) = SIGN(1.0_bp, cos_angles(1))
       END IF
       i = ACOS(cos_angles(1))
       IF (PRESENT(qei_only)) THEN
          IF (qei_only) THEN
             ! Return (q,e,i) only
             CALL NULLIFY(this_)
             getCometaryElements = (/ q, e, i, -1.0_bp, -1.0_bp, -1.0_bp /)
             RETURN
          END IF
       END IF
       ! Notice computational abs!!
       ru = SQRT(ABS(1.0_bp - k(3)**2.0_bp))
       IF (ru < EPSILON(ru)) THEN
          CALL errorMessage("Orbit / getCometaryElements", &
               "Warning: Inclination (almost) zero, setting Longitude of Node to 0 deg.", 1)
          an = 0.0_bp          
       ELSE
          cos_angles(2) = -k(2)/ru
          IF (ABS(cos_angles(2)) > 1.0_bp) THEN
             cos_angles(2) = SIGN(1.0_bp, cos_angles(2))
          END IF
          an = ACOS(cos_angles(2))
          IF (k(1) < 0.0_bp) THEN
             an = two_pi - an
          END IF
          IF (an == two_pi) THEN
             an = 0.0_bp
          END IF
       END IF

       ! Argument of periapsis:
       div = 1.0_bp + k(3)
       fb(1) = 1.0_bp - k(1)**2.0_bp/div
       fb(2) = -k(1)*k(2)/div
       fb(3) = -k(1)
       gb = cross_product(k, fb)
       tmp1 = DOT_PRODUCT(evec, gb)
       tmp2 = DOT_PRODUCT(evec, fb)
       varpi = ATAN2(tmp1, tmp2)
       ap = varpi - an
       ap = MODULO(ap,two_pi)
       IF (ap == two_pi) THEN
          ap = 0.0_bp
       END IF

       ! Time of periapsis
       alpha = 2.0_bp*planetary_mu(this_%center)/r0 - &
            DOT_PRODUCT(vel,vel)
       xv = DOT_PRODUCT(pos,vel)
       IF (alpha <= 0.0_bp) THEN
          ! Hyperbolic orbit
          cosu = (1.0_bp+r0/a)/e
          IF (cosu <= -1.0_bp) THEN
             u0 = pi
          ELSE IF (cosu >= 1.0_bp) THEN
             u0 = 0.0_bp
          ELSE
             u0 = ACOS(cosu)
          END IF
          s = -SIGN(1.0_bp,xv)*u0/SQRT(-alpha)
       ELSE
          ! Elliptic orbit
          IF (ABS(xv) > tol1) THEN
             s = (q-r0)/xv
          ELSE
             s = 0.0_bp
          END IF
       END IF
       !write(*,*) 
       !write(*,*) 
       !write(*,*) q,e,i/rad_deg,an/rad_deg,ap/rad_deg
       !write(*,*) 

       DO iiter=1,max_iter
          x = s**2.0_bp * alpha
          CALL getStumpffFunctions(x, stumpff_c)
          IF (error) THEN
             CALL errorMessage("Orbit / getCometaryElements", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
          stumpff_cs(0) = stumpff_c(0)
          stumpff_cs(1) = stumpff_c(1) * s
          stumpff_cs(2) = stumpff_c(2) * s**2.0_bp
          stumpff_cs(3) = stumpff_c(3) * s**3.0_bp
          r = r0 * stumpff_cs(0) + xv * stumpff_cs(1) + &
               planetary_mu(this_%center) * stumpff_cs(2)
          rp = (-r0*alpha + planetary_mu(this_%center)) * &
               stumpff_cs(1) + xv * stumpff_cs(0)
          rpp = (-r0*alpha + planetary_mu(this_%center)) * &
               stumpff_cs(0) - xv * alpha * stumpff_cs(1)
          ds = -rp/rpp
          s = s + ds
          !write(*,*) s, x, q, r, rp, rpp, ds
          IF (ABS(ds) < tol2 .AND. ABS(rp) < tol3 .AND. &
               (ABS(r-q) < tol4 .OR. (e < 1.0_bp .AND. &
               ABS(r-q/(1.0_bp-e)*(1.0_bp+e)) < tol4))) THEN
             EXIT
          END IF
       END DO
       IF (ABS(ds) > tol2 .OR. ABS(rp) > tol3 .OR. &
            (e >= 1.0_bp .AND. ABS(r-q) > tol4) .OR. &
            (e < 1.0_bp .AND. ABS(r-q) > tol4 .AND. &
            ABS(r-q/(1.0_bp-e)*(1.0_bp+e)) > tol4)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getCometaryElements", &
               "Computation of time of perihelion did not converge.", 1)
          IF (err_verb >= 1) THEN
             WRITE(stderr,*) "abs(ds) > tol2: ", ABS(ds) > tol2, &
                  "(", ABS(ds), " > ", tol2, ") or"
             WRITE(stderr,*) "abs(rp) > tol3: ", ABS(rp) > tol3, &
                  "(", ABS(rp), " > ", tol3, ") or"
             WRITE(stderr,*) "e >= 1 and abs(r-q) > tol4: ", &
                  e >= 1.0_bp .AND. ABS(r-q) > tol4, &
                  "(", ABS(r-q), " > ", tol4, ")"
             WRITE(stderr,*) "e < 1 and abs(r-q) > tol4 and abs(r-Q) > tol4: ", &
                  e < 1.0_bp .AND. ABS(r-q) > tol4 .AND. ABS(r-q/(1.0_bp-e)*(1.0_bp+e)) > tol4, &
                  "(", ABS(r-q), " > ", tol4, " and ", ABS(r-q/(1.0_bp-e)*(1.0_bp+e)), " > ", tol4, ")"
          END IF
          RETURN          
       END IF
       dt = r0 * stumpff_cs(1) + xv * stumpff_cs(2) + &
            planetary_mu(this_%center)*stumpff_cs(3)
       mjd_tt = getMJD(this_%t,"TT")
       ! Choose the time of perihelion closest to the epoch in case of
       ! an elliptic orbit:
       IF (e < 1.0_bp) THEN
          ! period p:
          p = two_pi*SQRT(a**3.0_bp/planetary_mu(this_%center))
          ! Remove full periods from dt:
          IF (ABS(dt) > p) THEN
             dt = MODULO(dt,p)
          END IF
          ! If dt refers to the time of aphelion, subtract (if dt>0,
          ! otherwise add) half a period to get to perihelion (scale
          ! tolerance with e to make sure ~circular orbits are
          ! processed correctly):
          IF (ABS(r-q/(1-e)*(1+e)) < tol4*e) THEN
             dt = dt - SIGN(1.0_bp,dt) * 0.5_bp * p
          END IF
          ! Select the closest time of perihelion:
          IF (ABS(dt) < ABS(dt - SIGN(1.0_bp,dt)*p)) THEN
             tp = mjd_tt + dt
          ELSE
             tp = mjd_tt + dt - SIGN(1.0_bp,dt)*p
          END IF
       ELSE
          ! For hyperbolic orbits always:
          tp = mjd_tt + dt
       END IF
       CALL NULLIFY(this_)

       getCometaryElements = (/ q, e, i, an, ap, tp /)

    END SELECT

  END FUNCTION getCometaryElements





  !! *Description*:
  !!
  !! Returns Delaunay's elements calculated from Keplerian elements. The mass 
  !! of the target body is assumed to be zero compared to the mass of the Sun.
  !!
  !!   - de(1) = l      = Mean Motion (Mean Anomaly)
  !!   - de(2) = g      = Argument of Perihelion
  !!   - de(3) = &theta = Longitude of the Ascending Node
  !!   - de(4) = L      = sqrt(mu*(Semimajor Axis)) = related to the 
  !!                                                  two-body orbital energy
  !!   - de(5) = G      = L*sqrt(1-Eccentricity^2)  = magnitude of the orbital 
  !!                                                  angular momentum
  !!   - de(6) = &Theta = G*cos(Inclination)        = z-component of the orbital 
  !!                                                  angular momentum
  !!
  !! Returns error.
  !! 
  FUNCTION getDelaunayElements(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(6)   :: getDelaunayElements
    REAL(bp), DIMENSION(6)   :: kep

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getDelaunayElements", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    kep = getElements(this, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / getDelaunayElements", &
            "TRACE BACK", 1)
       RETURN
    END IF
    getDelaunayElements(1) = kep(6)
    getDelaunayElements(2) = kep(5)
    getDelaunayElements(3) = kep(4)
    getDelaunayElements(4) = SQRT(planetary_mu(this%center) * kep(1))
    getDelaunayElements(5) = getDelaunayElements(4) * SQRT(1.0_bp - &
         kep(2)**2.0_bp)
    getDelaunayElements(6) = getDelaunayElements(5) * &
         COS(kep(3))

  END FUNCTION getDelaunayElements





  FUNCTION getElements(this, element_type, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)        :: this
    CHARACTER(len=*), INTENT(in)    :: element_type
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame
    REAL(bp), DIMENSION(6)          :: getElements
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type_
    REAL(bp), DIMENSION(6)          :: elem

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getElements", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (LEN_TRIM(element_type) > ELEMENT_TYPE_LEN) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getElements", &
            "Element type indicator too long.", 1)
       RETURN       
    END IF

    element_type_ = element_type
    CALL locase(element_type_, error)
    IF (error) THEN
       CALL errorMessage("Orbit / getElements", &
            "The element type string ('" // TRIM(element_type_) &
            // "') contains forbidden characters.", 1)
       RETURN
    END IF
    SELECT CASE (element_type_)
    CASE ("cometary")
       getElements = getCometaryElements(this)
    CASE ("cometary_ta")
       getElements = getCometaryElements(this)
    CASE ("keplerian")
       getElements = getKeplerianElements(this)
    CASE ("cartesian")
       IF (PRESENT(frame)) THEN
          getElements = getCartesianElements(this, frame)
       ELSE IF (this%element_type == "cartesian") THEN
          getElements = getCartesianElements(this, this%frame)
       ELSE
          error = .TRUE.
          CALL errorMessage("Orbit / getElements", &
               "Frame must be given explicitly if Keplerian " // &
               "elements are to be converted to Cartesian elements.", 1)
          RETURN
       END IF
    CASE ("delaunay")
       getElements = getDelaunayElements(this)
    CASE ("poincare")
       getElements = getPoincareElements(this)
    CASE ("equinoctial")
       ! Ref. Milani 1999, Icarus 137, pp. 269--292
       elem = getKeplerianElements(this)
       getElements(1) = elem(1)
       getElements(2) = elem(2)*SIN(elem(4) + elem(5))
       getElements(3) = elem(2)*COS(elem(4) + elem(5))
       getElements(4) = TAN(elem(3)/2.0_bp)*SIN(elem(4))
       getElements(5) = TAN(elem(3)/2.0_bp)*COS(elem(4))
       getElements(6) = MODULO(SUM(elem(4:6)), two_pi)
    CASE default
       error = .TRUE.
       CALL errorMessage("Orbit / getElements", &
            "Unkown option: " // TRIM(element_type_) // ".", 1)
       RETURN       
    END SELECT
    IF (error) THEN
       CALL errorMessage("Orbit / getElements", &
            "TRACE BACK", 1)
       RETURN       
    END IF

  END FUNCTION getElements





  CHARACTER(len=ELEMENT_TYPE_LEN) FUNCTION getElementType(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getElementType", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getElementType = this%element_type

  END FUNCTION getElementType





  !! *Description*:
  !!
  !! Computes spherical topocentric equatorial light-time-corrected
  !! coordinates from the object's heliocentric Cartesian orbit and 
  !! the observers" heliocentric Cartesian coordinates. Both 
  !! Cartesian elements and coordinates can be given either as 
  !! equatorial or ecliptic. Optionally, partial derivatives of the
  !! spherical topocentric equatorial light-time-corrected coordinates
  !! wrt object's orbital elements are calculated. The output format
  !! of partial derivatives is (row;column;observer). Optionally,
  !! the light-time corrected orbit is returned.
  !! 
  !! Returns error.
  !!
  SUBROUTINE getEphemerides_Orb_single(this, observers, ephemerides, &
       lt_corr, partials_arr, this_lt_corr_arr)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                              :: this
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: observers
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER    :: ephemerides
    LOGICAL, INTENT(in), OPTIONAL                         :: lt_corr
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL         :: partials_arr
    TYPE (Orbit), DIMENSION(:), POINTER, OPTIONAL         :: this_lt_corr_arr

    TYPE (Orbit), DIMENSION(1)                           :: this_arr
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: ephemerides_arr => NULL()
    REAL(bp), DIMENSION(:,:,:,:), POINTER                :: partials_arr_ => NULL()
    TYPE (Orbit), DIMENSION(:,:), POINTER                :: this_lt_corr_arr_ => NULL()
    INTEGER                                              :: i, err, nobsy
    LOGICAL                                              :: lt_corr_

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (single)", &
            "Object has not yet been initialized.", 1)
       RETURN       
    END IF

    IF (PRESENT(lt_corr)) THEN
       lt_corr_ = lt_corr
    ELSE
       lt_corr_ = .TRUE.
    END IF
    nobsy = SIZE(observers,dim=1)
    this_arr(1) = copy(this)
    IF (PRESENT(partials_arr)) THEN 
       IF (PRESENT(this_lt_corr_arr)) THEN
          CALL getEphemerides(this_arr, observers, ephemerides_arr, &
               lt_corr=lt_corr_, partials_arr=partials_arr_, &
               this_lt_corr_arr=this_lt_corr_arr_)
       ELSE
          CALL getEphemerides(this_arr, observers, ephemerides_arr, &
               lt_corr=lt_corr_, partials_arr=partials_arr_)
       END IF
       IF (error) THEN
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "TRACE BACK (5)", 1)
          DEALLOCATE(partials_arr_, stat=err)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(this_lt_corr_arr_, stat=err)
          RETURN
       END IF
       ALLOCATE(partials_arr(6,6,nobsy), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "Could not allocate memory (5).", 1)
          DEALLOCATE(partials_arr_, stat=err)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(this_lt_corr_arr_, stat=err)
          RETURN       
       END IF
       partials_arr = partials_arr_(1,:,:,:)
       DEALLOCATE(partials_arr_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "Could not deallocate memory (5).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(this_lt_corr_arr_, stat=err)
          RETURN       
       END IF
    ELSE
       IF (PRESENT(this_lt_corr_arr)) THEN
          CALL getEphemerides(this_arr, observers, ephemerides_arr, &
               lt_corr=lt_corr_, &
               this_lt_corr_arr=this_lt_corr_arr_)
       ELSE
          CALL getEphemerides(this_arr, observers, ephemerides_arr, &
               lt_corr=lt_corr_)
       END IF
       IF (error) THEN
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "TRACE BACK (10)", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          DEALLOCATE(this_lt_corr_arr_, stat=err)
          RETURN
       END IF
    END IF
    IF (PRESENT(this_lt_corr_arr)) THEN
       ALLOCATE(this_lt_corr_arr(nobsy), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          DEALLOCATE(this_lt_corr_arr_, stat=err)
          DEALLOCATE(this_lt_corr_arr, stat=err)
          RETURN       
       END IF
       DO i=1,nobsy
          this_lt_corr_arr(i) = copy(this_lt_corr_arr_(1,i))
       END DO
       DEALLOCATE(this_lt_corr_arr_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          DEALLOCATE(this_lt_corr_arr, stat=err)
          RETURN       
       END IF
    END IF
    IF (ASSOCIATED(ephemerides)) THEN
       DEALLOCATE(ephemerides, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (single)", &
               "Could not allocate memory (12).", 1)
          DEALLOCATE(ephemerides_arr, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          RETURN       
       END IF
    END IF
    ALLOCATE(ephemerides(nobsy), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (single)", &
            "Could not allocate memory (15).", 1)
       DEALLOCATE(ephemerides_arr, stat=err)
       DEALLOCATE(ephemerides, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       RETURN       
    END IF
    DO i=1,nobsy
       ephemerides(i) = copy(ephemerides_arr(1,i))
    END DO
    DEALLOCATE(ephemerides_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (single)", &
            "Could not allocate memory (15).", 1)
       DEALLOCATE(ephemerides, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       RETURN       
    END IF
    CALL NULLIFY(this_arr(1))

  END SUBROUTINE getEphemerides_Orb_single





  !! *Description*:
  !!
  !! Computes spherical topocentric equatorial light-time-corrected
  !! coordinates from the object's heliocentric Cartesian orbit and 
  !! the observers" heliocentric Cartesian coordinates. Both 
  !! Cartesian elements and coordinates can be given either as 
  !! equatorial or ecliptic. Optionally, partial derivatives of the
  !! spherical topocentric equatorial light-time-corrected coordinates
  !! wrt object's orbital elements are calculated. The output format
  !! of partial derivatives is (orbit;row;column;observer). Optionally,
  !! the light-time corrected orbit is returned.
  !! 
  !! Returns error.
  !!
  SUBROUTINE getEphemerides_Orb_multiple(this_arr, observers, ephemerides, &
       lt_corr, partials_arr, this_lt_corr_arr)

    IMPLICIT NONE
    TYPE (Orbit), DIMENSION(:), INTENT(in)                :: this_arr
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: observers
    TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER  :: ephemerides
    LOGICAL, INTENT(in), OPTIONAL                         :: lt_corr
    REAL(bp), DIMENSION(:,:,:,:), POINTER, OPTIONAL       :: partials_arr
    TYPE (Orbit), DIMENSION(:,:), POINTER, OPTIONAL       :: this_lt_corr_arr

    TYPE (Orbit), DIMENSION(:), POINTER                 :: this_prop_arr_next => NULL()
    TYPE (Orbit), DIMENSION(:), ALLOCATABLE             :: this_prop_arr_current
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER  :: ephemerides_ => NULL()
    TYPE (Orbit), DIMENSION(:), POINTER                 :: this_lt_corr_arr_ => NULL()
    TYPE (Time)                                         :: t
    CHARACTER(len=FRAME_LEN), DIMENSION(:), ALLOCATABLE :: frames
    REAL(bp), DIMENSION(:,:,:), POINTER                 :: partials_arr_ => NULL(), &
         jacobian_prop_arr => NULL()
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE             :: jacobian_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE                 :: mjd_tt_arr
    REAL(bp)                                            :: mjd_tt_orb
    INTEGER, DIMENSION(:), ALLOCATABLE                  :: indx_arr, indx_smaller, indx_larger
    INTEGER                                             :: i, j, k, err, nthis, nobsy, indx
    LOGICAL                                             :: lt_corr_

    nthis = SIZE(this_arr)
    nobsy = SIZE(observers)
    ALLOCATE(this_prop_arr_current(nthis), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(this_prop_arr_current, stat=err)
       RETURN
    END IF
    DO i=1,nthis
       IF (.NOT.this_arr(i)%is_initialized) THEN
          error = .TRUE.
          DEALLOCATE(this_prop_arr_current, stat=err)
          CALL errorMessage("Orbit / getEphemerides (multiple)", &
               "All objects have not yet been initialized.", 1)
          RETURN       
       END IF
       this_prop_arr_current(i) = copy(this_arr(i))
    END DO
    ALLOCATE(ephemerides(nthis,nobsy), frames(nthis), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            "Could not allocate memory (10).", 1)
       DEALLOCATE(this_prop_arr_current, stat=err)
       DEALLOCATE(ephemerides, stat=err)
       DEALLOCATE(frames, stat=err)
       RETURN
    END IF
    IF (PRESENT(lt_corr)) THEN
       lt_corr_ = lt_corr
    ELSE
       lt_corr_ = .TRUE.
    END IF
    IF (PRESENT(partials_arr)) THEN
       ALLOCATE(partials_arr(nthis,6,6,SIZE(observers)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (multiple)", &
               "Could not allocate memory (15).", 1)
          DEALLOCATE(this_prop_arr_current, stat=err)
          DEALLOCATE(ephemerides, stat=err)
          DEALLOCATE(frames, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          RETURN
       END IF
       ALLOCATE(jacobian_arr(nthis,6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (multiple)", &
               "Could not allocate memory (20).", 1)
          DEALLOCATE(this_prop_arr_current, stat=err)
          DEALLOCATE(ephemerides, stat=err)
          DEALLOCATE(frames, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          DEALLOCATE(jacobian_arr, stat=err)
          RETURN
       END IF
       DO i=1,nthis
          jacobian_arr(i,:,:) = identity_matrix(6)
       END DO
    END IF
    IF (PRESENT(this_lt_corr_arr)) THEN
       ALLOCATE(this_lt_corr_arr(nthis,nobsy), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (multiple)", &
               "Could not allocate memory (25).", 1)
          DEALLOCATE(this_prop_arr_current, stat=err)
          DEALLOCATE(ephemerides, stat=err)
          DEALLOCATE(frames, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_arr, stat=err)
          DEALLOCATE(this_lt_corr_arr, stat=err)
          RETURN
       END IF
    END IF

    ! This section splits ephemeris dates to two different groups:
    ! 1) those earlier than the orbit epoch and 2) those later than
    ! the orbit epoch:
    mjd_tt_orb = getMJD(this_prop_arr_current(1)%t, "TT")
    ALLOCATE(mjd_tt_arr(nobsy), indx_arr(nobsy), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            "Could not allocate memory (30).", 1)
       DEALLOCATE(this_prop_arr_current, stat=err)
       DEALLOCATE(ephemerides, stat=err)
       DEALLOCATE(frames, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       DEALLOCATE(jacobian_arr, stat=err)
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       DEALLOCATE(mjd_tt_arr, stat=err)
       DEALLOCATE(indx_arr, stat=err)
       RETURN
    END IF
    DO i=1,nobsy
       t = getTime(observers(i))
       mjd_tt_arr(i) = getMJD(t, "TT")
       CALL NULLIFY(t)
    END DO
    CALL quicksort(mjd_tt_arr, indx_arr, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            TRIM(errstr), 1)
       RETURN
    END IF
    indx = findLocation(mjd_tt_orb, mjd_tt_arr, indx_arr)
    ! NB: Here we use tricks possible with Fortran. The arrays below 
    ! have indices from 1:indx-1 and indx:nobsy, that is, the second
    ! one does not generally start from 0 or 1
    ALLOCATE(indx_smaller(1:indx-1), &
         indx_larger(indx:SIZE(observers)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            "Could not allocate memory (35).", 1)
       DEALLOCATE(this_prop_arr_current, stat=err)
       DEALLOCATE(ephemerides, stat=err)
       DEALLOCATE(frames, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       DEALLOCATE(jacobian_arr, stat=err)
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       DEALLOCATE(mjd_tt_arr, stat=err)
       DEALLOCATE(indx_arr, stat=err)
       DEALLOCATE(indx_smaller, stat=err)
       DEALLOCATE(indx_larger, stat=err)
       RETURN
    END IF
    indx_smaller(1:indx-1) = indx_arr(1:indx-1)
    indx_larger(indx:) = indx_arr(indx:)
    DEALLOCATE(mjd_tt_arr, indx_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            "Could not deallocate memory (5).", 1)
       DEALLOCATE(this_prop_arr_current, stat=err)
       DEALLOCATE(ephemerides, stat=err)
       DEALLOCATE(frames, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       DEALLOCATE(jacobian_arr, stat=err)
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       DEALLOCATE(mjd_tt_arr, stat=err)
       DEALLOCATE(indx_arr, stat=err)
       DEALLOCATE(indx_smaller, stat=err)
       DEALLOCATE(indx_larger, stat=err)
       RETURN
    END IF

    DO i=1,nthis
       IF (this_prop_arr_current(i)%element_type == "cartesian" &
            .AND. this_prop_arr_current(i)%frame == "ecliptic") THEN
          frames(i) = this_prop_arr_current(i)%frame
          CALL rotateToEquatorial(this_prop_arr_current(i))
       END IF
    END DO
    DO i=1, nobsy
       IF (i <= SIZE(indx_smaller)) THEN
          ! First propagate backwards from the orbit epoch 
          j = indx_smaller(SIZE(indx_smaller)-i+1)
       ELSE IF (i == SIZE(indx_smaller) + 1) THEN
          ! Then propagate forwards from the orbit epoch
          DO j=1,nthis
             CALL NULLIFY(this_prop_arr_current(j))
             this_prop_arr_current(j) = copy(this_arr(j))
             IF (PRESENT(partials_arr)) THEN
                jacobian_arr(j,:,:) = identity_matrix(6)
             END IF
          END DO
          j = indx_larger(i)
       ELSE
          j = indx_larger(i)
       END IF
       IF (PRESENT(partials_arr)) THEN 
          IF (PRESENT(this_lt_corr_arr)) THEN
             CALL getEphemeris(this_prop_arr_current, &
                  observers(j), ephemerides_, &
                  lt_corr=lt_corr_, &
                  partials_arr=partials_arr_, &
                  jacobian_prop_arr=jacobian_prop_arr, &
                  this_lt_corr_arr=this_lt_corr_arr_, &
                  this_prop_arr=this_prop_arr_next)
             IF (error) THEN
                CALL errorMessage("Orbit / getEphemerides (multiple)", &
                     "TRACE BACK (5)", 1)
                DEALLOCATE(this_prop_arr_current, stat=err)
                DEALLOCATE(ephemerides, stat=err)
                DEALLOCATE(frames, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(jacobian_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr, stat=err)
                DEALLOCATE(indx_smaller, stat=err)
                DEALLOCATE(indx_larger, stat=err)
                DEALLOCATE(ephemerides_, stat=err)
                DEALLOCATE(partials_arr_, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr_, stat=err)
                DEALLOCATE(this_prop_arr_next, stat=err)
                RETURN
             END IF
             DO k=1,nthis
                this_lt_corr_arr(k,j) = copy(this_lt_corr_arr_(k))
                CALL NULLIFY(this_lt_corr_arr_(k))
             END DO
             DEALLOCATE(this_lt_corr_arr_, stat=err)
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / getEphemerides (multiple)", &
                     "Could not deallocate memory (10).", 1)
                DEALLOCATE(this_prop_arr_current, stat=err)
                DEALLOCATE(ephemerides, stat=err)
                DEALLOCATE(frames, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(jacobian_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr, stat=err)
                DEALLOCATE(indx_smaller, stat=err)
                DEALLOCATE(indx_larger, stat=err)
                DEALLOCATE(ephemerides_, stat=err)
                DEALLOCATE(partials_arr_, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr_, stat=err)
                DEALLOCATE(this_prop_arr_next, stat=err)
                RETURN
             END IF
          ELSE
             CALL getEphemeris(this_prop_arr_current, &
                  observers(j), ephemerides_, &
                  lt_corr=lt_corr_, &
                  partials_arr=partials_arr_, &
                  jacobian_prop_arr=jacobian_prop_arr, &
                  this_prop_arr=this_prop_arr_next)
             IF (error) THEN
                CALL errorMessage("Orbit / getEphemerides (multiple)", &
                     "TRACE BACK (10)", 1)
                DEALLOCATE(this_prop_arr_current, stat=err)
                DEALLOCATE(ephemerides, stat=err)
                DEALLOCATE(frames, stat=err)
                DEALLOCATE(partials_arr, stat=err)
                DEALLOCATE(jacobian_arr, stat=err)
                IF (PRESENT(this_lt_corr_arr)) THEN
                   DEALLOCATE(this_lt_corr_arr, stat=err)
                END IF
                DEALLOCATE(indx_smaller, stat=err)
                DEALLOCATE(indx_larger, stat=err)
                DEALLOCATE(ephemerides_, stat=err)
                DEALLOCATE(partials_arr_, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr_, stat=err)
                DEALLOCATE(this_prop_arr_next, stat=err)
                RETURN
             END IF
          END IF
          DO k=1,nthis
             partials_arr(k,:,:,j) = MATMUL(partials_arr_(k,:,:), jacobian_arr(k,:,:))
             jacobian_arr(k,:,:) = MATMUL(jacobian_prop_arr(k,:,:), jacobian_arr(k,:,:))
          END DO
          DEALLOCATE(partials_arr_, jacobian_prop_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getEphemerides (multiple)", &
                  "Could not deallocate memory (15).", 1)
             DEALLOCATE(this_prop_arr_current, stat=err)
             DEALLOCATE(ephemerides, stat=err)
             DEALLOCATE(frames, stat=err)
             DEALLOCATE(partials_arr, stat=err)
             DEALLOCATE(jacobian_arr, stat=err)
             IF (PRESENT(this_lt_corr_arr)) THEN
                DEALLOCATE(this_lt_corr_arr, stat=err)
             END IF
             DEALLOCATE(indx_smaller, stat=err)
             DEALLOCATE(indx_larger, stat=err)
             DEALLOCATE(ephemerides_, stat=err)
             DEALLOCATE(partials_arr_, stat=err)
             DEALLOCATE(jacobian_prop_arr, stat=err)
             DEALLOCATE(this_lt_corr_arr_, stat=err)
             DEALLOCATE(this_prop_arr_next, stat=err)
             RETURN
          END IF
       ELSE
          IF (PRESENT(this_lt_corr_arr)) THEN
             CALL getEphemeris(this_prop_arr_current, &
                  observers(j), ephemerides_, &
                  lt_corr=lt_corr_, &
                  this_lt_corr_arr=this_lt_corr_arr_, &
                  this_prop_arr=this_prop_arr_next)
             IF (error) THEN
                CALL errorMessage("Orbit / getEphemerides (multiple)", &
                     "TRACE BACK (15)", 1)
                DEALLOCATE(this_prop_arr_current, stat=err)
                DEALLOCATE(ephemerides, stat=err)
                DEALLOCATE(frames, stat=err)
                DEALLOCATE(jacobian_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr, stat=err)
                DEALLOCATE(indx_smaller, stat=err)
                DEALLOCATE(indx_larger, stat=err)
                DEALLOCATE(ephemerides_, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr_, stat=err)
                DEALLOCATE(this_prop_arr_next, stat=err)
                RETURN
             END IF
             DO k=1,nthis
                this_lt_corr_arr(k,j) = copy(this_lt_corr_arr_(k))
                CALL NULLIFY(this_lt_corr_arr_(k))
             END DO
             DEALLOCATE(this_lt_corr_arr_, stat=err)
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / getEphemerides (multiple)", &
                     "Could not deallocate memory (20).", 1)
                DEALLOCATE(this_prop_arr_current, stat=err)
                DEALLOCATE(ephemerides, stat=err)
                DEALLOCATE(frames, stat=err)
                DEALLOCATE(jacobian_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr, stat=err)
                DEALLOCATE(indx_smaller, stat=err)
                DEALLOCATE(indx_larger, stat=err)
                DEALLOCATE(ephemerides_, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr_, stat=err)
                DEALLOCATE(this_prop_arr_next, stat=err)
                RETURN
             END IF
          ELSE
             CALL getEphemeris(this_prop_arr_current, &
                  observers(j), ephemerides_, &
                  lt_corr=lt_corr_, &
                  this_prop_arr=this_prop_arr_next)
             IF (error) THEN
                CALL errorMessage("Orbit / getEphemerides (multiple)", &
                     "TRACE BACK (20)", 1)
                DEALLOCATE(this_prop_arr_current, stat=err)
                DEALLOCATE(ephemerides, stat=err)
                DEALLOCATE(frames, stat=err)
                DEALLOCATE(jacobian_arr, stat=err)
                IF (PRESENT(this_lt_corr_arr)) THEN
                   DEALLOCATE(this_lt_corr_arr, stat=err)
                END IF
                DEALLOCATE(indx_smaller, stat=err)
                DEALLOCATE(indx_larger, stat=err)
                DEALLOCATE(ephemerides_, stat=err)
                DEALLOCATE(jacobian_prop_arr, stat=err)
                DEALLOCATE(this_lt_corr_arr_, stat=err)
                DEALLOCATE(this_prop_arr_next, stat=err)
                RETURN
             END IF
          END IF
       END IF
       DO k=1,nthis
          ephemerides(k,j) = copy(ephemerides_(k))
          CALL NULLIFY(ephemerides_(k))
          CALL NULLIFY(this_prop_arr_current(k))
          this_prop_arr_current(k) = copy(this_prop_arr_next(k))
          CALL NULLIFY(this_prop_arr_next(k))
       END DO
       DEALLOCATE(ephemerides_, this_prop_arr_next, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (multiple)", &
               "Could not deallocate memory (25).", 1)
          DEALLOCATE(this_prop_arr_current, stat=err)
          DEALLOCATE(ephemerides, stat=err)
          DEALLOCATE(frames, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_arr, stat=err)
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          DEALLOCATE(indx_smaller, stat=err)
          DEALLOCATE(indx_larger, stat=err)
          DEALLOCATE(ephemerides_, stat=err)
          DEALLOCATE(partials_arr_, stat=err)
          DEALLOCATE(jacobian_prop_arr, stat=err)
          DEALLOCATE(this_lt_corr_arr_, stat=err)
          DEALLOCATE(this_prop_arr_next, stat=err)
          RETURN
       END IF
    END DO
    DO i=1,nthis
       IF (PRESENT(partials_arr) .AND. &
            this_prop_arr_current(i)%element_type == "cartesian" &
            .AND. frames(i) == "ecliptic") THEN
          DO j=1,SIZE(partials_arr,dim=4)
             DO k=1,6
                CALL rotateToEcliptic(partials_arr(i,k,:,j))
             END DO
          END DO
       END IF
    END DO

    DO j=1,SIZE(this_prop_arr_current)
       CALL NULLIFY(this_prop_arr_current(j))
    END DO
    DEALLOCATE(frames, this_prop_arr_current, indx_smaller, &
         indx_larger, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemerides (multiple)", &
            "Could not deallocate memory (30).", 1)
       DEALLOCATE(this_prop_arr_current, stat=err)
       DEALLOCATE(ephemerides, stat=err)
       DEALLOCATE(frames, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       DEALLOCATE(jacobian_arr, stat=err)
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       DEALLOCATE(indx_smaller, stat=err)
       DEALLOCATE(indx_larger, stat=err)
       DEALLOCATE(ephemerides_, stat=err)
       DEALLOCATE(partials_arr_, stat=err)
       DEALLOCATE(jacobian_prop_arr, stat=err)
       DEALLOCATE(this_lt_corr_arr_, stat=err)
       DEALLOCATE(this_prop_arr_next, stat=err)
       RETURN
    END IF
    IF (ALLOCATED(jacobian_arr)) THEN
       DEALLOCATE(jacobian_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemerides (multiple)", &
               "Could not deallocate memory (35).", 1)
          DEALLOCATE(ephemerides, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          RETURN
       END IF
    END IF

  END SUBROUTINE getEphemerides_Orb_multiple





  !! *Description*:
  !!
  !! Computes spherical topocentric equatorial light-time-corrected
  !! coordinates from the object's heliocentric Cartesian orbit and 
  !! the observer's heliocentric Cartesian coordinates. Both 
  !! Cartesian elements and coordinates can be given either as 
  !! equatorial or ecliptic. Optionally, partial derivatives of the
  !! spherical topocentric equatorial light-time-corrected coordinates
  !! wrt object's orbital elements are calculated. Optionally,
  !! the light-time corrected orbit is returned. 
  !!
  !! The given orbit is not propagated to the ephemeris epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE getEphemeris_Orb_single(this, observer, ephemeris, &
       lt_corr, partials, this_prop, jacobian_prop, this_lt_corr, &
       jacobian_lt_corr)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                     :: this
    TYPE (CartesianCoordinates), INTENT(in)         :: observer
    TYPE (SphericalCoordinates), INTENT(out)        :: ephemeris
    LOGICAL, INTENT(in), OPTIONAL                   :: lt_corr
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: partials
    TYPE (Orbit), INTENT(out), OPTIONAL             :: this_prop
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: jacobian_prop
    TYPE (Orbit), INTENT(out), OPTIONAL             :: this_lt_corr
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: jacobian_lt_corr

    TYPE (Orbit), DIMENSION(1)                         :: this_arr
    TYPE (Orbit), DIMENSION(:), POINTER                :: this_prop_arr => NULL()
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: ephemeris_arr => NULL()
    TYPE (Orbit), DIMENSION(:), POINTER                :: this_lt_corr_arr => NULL()
    REAL(bp), DIMENSION(:,:,:), POINTER                :: partials_arr => NULL()
    REAL(bp), DIMENSION(:,:,:), POINTER                :: jacobian_prop_arr => NULL(), &
         jacobian_lt_corr_arr => NULL()
    INTEGER :: i, err
    LOGICAL :: lt_corr_

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemeris (single)", &
            "This object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemeris (single)", &
            "'observer' has not been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(lt_corr)) THEN
       lt_corr_ = lt_corr
    ELSE
       lt_corr_ = .TRUE.
    END IF
    this_arr(1) = copy(this)
    IF (PRESENT(partials)) THEN 
       IF (PRESENT(jacobian_prop)) THEN
          IF (PRESENT(this_prop)) THEN
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             END IF
          ELSE
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             END IF
          END IF
       ELSE
          IF (PRESENT(this_prop)) THEN
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_prop_arr=this_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        this_prop_arr=this_prop_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_prop_arr=this_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_prop_arr=this_prop_arr)
                END IF
             END IF
          ELSE
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_lt_corr_arr=this_lt_corr_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr, &
                        this_lt_corr_arr=this_lt_corr_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        partials_arr=partials_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        partials_arr=partials_arr)
                END IF
             END IF
          END IF
       END IF
    ELSE
       IF (PRESENT(jacobian_prop)) THEN
          IF (PRESENT(this_prop)) THEN
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_prop_arr=this_prop_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             END IF
          ELSE
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        jacobian_prop_arr=jacobian_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        jacobian_prop_arr=jacobian_prop_arr)
                END IF
             END IF
          END IF
       ELSE
          IF (PRESENT(this_prop)) THEN
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_prop_arr=this_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr, &
                        this_prop_arr=this_prop_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_prop_arr=this_prop_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_prop_arr=this_prop_arr)
                END IF
             END IF
          ELSE
             IF (PRESENT(this_lt_corr)) THEN
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        this_lt_corr_arr=this_lt_corr_arr)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_, &
                        this_lt_corr_arr=this_lt_corr_arr)
                END IF
             ELSE
                IF (PRESENT(jacobian_lt_corr)) THEN
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
                        lt_corr=lt_corr_)
                ELSE
                   CALL getEphemeris(this_arr, observer, ephemeris_arr, &
                        lt_corr=lt_corr_)
                END IF
             END IF
          END IF
       END IF
    END IF
    CALL NULLIFY(this_arr(1))
    IF (error) THEN
       CALL errorMessage("Orbit / getEphemeris (single)", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(partials_arr, stat=err)
       DEALLOCATE(this_lt_corr_arr, stat=err)
       DEALLOCATE(this_prop_arr, stat=err)
       DEALLOCATE(jacobian_prop_arr, stat=err)
       DEALLOCATE(jacobian_lt_corr_arr, stat=err)
       DEALLOCATE(ephemeris_arr, stat=err)
       RETURN
    END IF
    IF (PRESENT(partials)) THEN 
       partials = partials_arr(1,:,:)
       DEALLOCATE(partials_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (single)", &
               "Could not deallocate memory (5)", 1)
          DEALLOCATE(this_lt_corr_arr, stat=err)
          DEALLOCATE(this_prop_arr, stat=err)
          DEALLOCATE(jacobian_prop_arr, stat=err)
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          DEALLOCATE(ephemeris_arr, stat=err)
          RETURN
       END IF
    END IF
    IF (PRESENT(this_lt_corr)) THEN
       this_lt_corr = copy(this_lt_corr_arr(1))
       DO i=1,SIZE(this_lt_corr_arr)
          CALL NULLIFY(this_lt_corr_arr(i))
       END DO
       DEALLOCATE(this_lt_corr_arr, stat=err)
    END IF
    IF (PRESENT(this_prop)) THEN
       this_prop = copy(this_prop_arr(1))
       DO i=1,SIZE(this_prop_arr)
          CALL NULLIFY(this_prop_arr(i))
       END DO
       DEALLOCATE(this_prop_arr, stat=err)
    END IF
    IF (PRESENT(jacobian_prop)) THEN
       jacobian_prop = jacobian_prop_arr(1,:,:)
       DEALLOCATE(jacobian_prop_arr, stat=err)
    END IF
    IF (PRESENT(jacobian_lt_corr)) THEN
       jacobian_lt_corr = jacobian_lt_corr_arr(1,:,:)
       DEALLOCATE(jacobian_lt_corr_arr, stat=err)
    END IF
    ephemeris = copy(ephemeris_arr(1))
    DEALLOCATE(ephemeris_arr, stat=err)

  END SUBROUTINE getEphemeris_Orb_single





  !! *Description*:
  !!
  !! Computes spherical topocentric equatorial light-time-corrected
  !! coordinates from the objects" heliocentric Cartesian orbits and 
  !! the observer's heliocentric Cartesian coordinates. Both 
  !! Cartesian elements and coordinates can be given either as 
  !! equatorial or ecliptic. Optionally, partial derivatives of the
  !! spherical topocentric equatorial light-time-corrected coordinates
  !! wrt object's orbital elements are calculated. Optionally,
  !! the light-time corrected orbit is returned.
  !!
  !! The input orbits are propagated to the ephemeris epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE getEphemeris_Orb_multiple(this_arr, observer, &
       ephemeris, lt_corr, partials_arr, this_prop_arr, &
       jacobian_prop_arr, this_lt_corr_arr, jacobian_lt_corr_arr)

    IMPLICIT NONE
    TYPE (Orbit), DIMENSION(:), INTENT(in)             :: this_arr
    TYPE (CartesianCoordinates), INTENT(in)            :: observer
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: ephemeris
    LOGICAL, INTENT(in), OPTIONAL                      :: lt_corr
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL      :: partials_arr
    TYPE (Orbit), DIMENSION(:), POINTER, OPTIONAL      :: this_prop_arr
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL      :: jacobian_prop_arr
    TYPE (Orbit), DIMENSION(:), POINTER, OPTIONAL      :: this_lt_corr_arr
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL      :: jacobian_lt_corr_arr

    TYPE (Orbit), DIMENSION(:), ALLOCATABLE :: this_arr_
    TYPE (Orbit)                            :: this_1, this_2
    TYPE (CartesianCoordinates)             :: observer_
    TYPE (Time)                             :: t_, t_observer
    CHARACTER(len=FRAME_LEN)                :: frame
    REAL(bp), DIMENSION(:,:,:), POINTER     :: jacobian_prop_arr_ => NULL()
    REAL(bp), DIMENSION(6,6)                :: jacobian, jacobian_lt_corr
    REAL(bp), DIMENSION(6,6)                :: scoord_partials
    REAL(bp), DIMENSION(6)                  :: observer_coordinates, elements
    INTEGER                                 :: i, nthis, err
    LOGICAL                                 :: lt_corr_

    nthis = SIZE(this_arr,dim=1)
    ALLOCATE(this_arr_(nthis), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemeris (multiple)", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(this_arr_, stat=err)
       RETURN
    END IF
    DO i=1,nthis
       IF (.NOT. this_arr(i)%is_initialized) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "All objects have not yet been initialized.", 1)
          DEALLOCATE(this_arr_, stat=err)
          RETURN
       END IF
       CALL NULLIFY(this_arr_(i))
       this_arr_(i) = copy(this_arr(i))
       ! Ensure that the coordinate system is heliocentric:
       IF (this_arr_(i)%center /= 11) THEN
          CALL switchCenter(this_arr_(i),11)
       END IF
    END DO
    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemeris (multiple)", &
            "'observer' has not been initialized.", 1)
       DEALLOCATE(this_arr_, stat=err)
       RETURN
    END IF
    observer_ = copy(observer)

    ALLOCATE(ephemeris(nthis), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemeris (multiple)", &
            "Could not allocate memory (5).", 1)
       DEALLOCATE(this_arr_, stat=err)
       DEALLOCATE(ephemeris, stat=err)
       CALL NULLIFY(observer_)
       RETURN
    END IF
    IF (PRESENT(lt_corr)) THEN
       lt_corr_ = lt_corr
    ELSE
       lt_corr_ = .TRUE.
    END IF
    IF (PRESENT(partials_arr)) THEN
       ALLOCATE(partials_arr(nthis,6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          DEALLOCATE(partials_arr, stat=err)
          CALL NULLIFY(observer_)
          RETURN
       END IF
    END IF
    IF (PRESENT(this_prop_arr)) THEN
       IF (ASSOCIATED(this_prop_arr)) THEN
          DEALLOCATE(this_prop_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getEphemeris (multiple)", &
                  "Could not deallocate memory (10).", 1)
             DEALLOCATE(this_arr_, stat=err)
             DEALLOCATE(ephemeris, stat=err)
             IF (PRESENT(partials_arr)) THEN
                DEALLOCATE(partials_arr, stat=err)
             END IF
             CALL NULLIFY(observer_)
             RETURN
          END IF
       END IF
       ALLOCATE(this_prop_arr(nthis), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          DEALLOCATE(this_prop_arr, stat=err)
          CALL NULLIFY(observer_)
          RETURN
       END IF
    END IF
    IF (PRESENT(jacobian_prop_arr)) THEN
       ALLOCATE(jacobian_prop_arr(nthis,6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "Could not allocate memory (10).", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          IF (PRESENT(this_prop_arr)) THEN
             DEALLOCATE(this_prop_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_prop_arr, stat=err)
          CALL NULLIFY(observer_)
          RETURN
       END IF
    END IF
    IF (PRESENT(this_lt_corr_arr)) THEN
       ALLOCATE(this_lt_corr_arr(nthis), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "Could not allocate memory (15).", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          IF (PRESENT(this_prop_arr)) THEN
             DEALLOCATE(this_prop_arr, stat=err)
          END IF
          DEALLOCATE(this_lt_corr_arr, stat=err)
          IF (PRESENT(jacobian_prop_arr)) THEN
             DEALLOCATE(jacobian_prop_arr, stat=err)
          END IF
          CALL NULLIFY(observer_)
          RETURN
       END IF
    END IF
    IF (PRESENT(jacobian_lt_corr_arr)) THEN
       ALLOCATE(jacobian_lt_corr_arr(nthis,6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "Could not allocate memory (15).", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          IF (PRESENT(this_prop_arr)) THEN
             DEALLOCATE(this_prop_arr, stat=err)
          END IF
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          IF (PRESENT(jacobian_prop_arr)) THEN
             DEALLOCATE(jacobian_prop_arr, stat=err)
          END IF
          CALL NULLIFY(observer_)
          RETURN
       END IF
       IF (.NOT.lt_corr_) THEN
          DO i=1,nthis
             jacobian_lt_corr_arr(i,:,:) = identity_matrix(6) 
          END DO
       END IF
    END IF

    ! Propagate to the observer epoch:
    IF (info_verb >= 4) THEN
       WRITE(stdout,"(1X,2(1X,A),I0,A)") &
            "Orbit / getEphemeris_multiple:", &
            "Propagate orbit array (size=", nthis, &
            ") to the observer epoch."
    END IF
    t_observer = getTime(observer_)
    IF (PRESENT(partials_arr) .OR. PRESENT(jacobian_prop_arr) &
         .OR. PRESENT(jacobian_lt_corr_arr)) THEN
       CALL propagate(this_arr_, t_observer, jacobian=jacobian_prop_arr_)
    ELSE
       CALL propagate(this_arr_, t_observer)
    END IF
    IF (error) THEN
       CALL errorMessage("Orbit / getEphemeris (multiple)", &
            "TRACE BACK (5)", 1)
       DEALLOCATE(this_arr_, stat=err)
       DEALLOCATE(ephemeris, stat=err)
       DEALLOCATE(jacobian_prop_arr_, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       IF (PRESENT(this_prop_arr)) THEN
          DEALLOCATE(this_prop_arr, stat=err)
       END IF
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       IF (PRESENT(jacobian_lt_corr_arr)) THEN
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
       END IF
       IF (PRESENT(jacobian_prop_arr)) THEN
          DEALLOCATE(jacobian_prop_arr, stat=err)
       END IF
       CALL NULLIFY(observer_)
       CALL NULLIFY(t_observer)
       RETURN
    END IF
    IF (PRESENT(jacobian_prop_arr)) THEN
       jacobian_prop_arr = jacobian_prop_arr_
    END IF

    ! Make the light-time correction (since the new epoch is orbit
    ! dependent, all orbits must be propagated individually):
    IF (info_verb >= 4 .AND. lt_corr_) THEN
       WRITE(stdout,"(1X,2(1X,A))") &
            "Orbit / getEphemeris_multiple:", &
            "Make light-time correction."
    END IF
    DO i=1,nthis
       IF (PRESENT(this_prop_arr)) THEN
          this_prop_arr(i) = copy(this_arr_(i))
       END IF
       this_1 = copy(this_arr_(i))
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,1X,5A)") &
               "Orbit / getEphemeris_multiple:", &
               "Orbital elements (" // &
               TRIM(this_1%element_type) // ", " // &
               TRIM(this_1%frame) // &
               ") at observation date without light-time correction:"
          elements = getElements(this_1, this_1%element_type, this_1%frame)
          IF (this_1%element_type == "keplerian") THEN
             WRITE(stdout,"(6(F14.10,1X))") elements(1:2), &
                  elements(3:6)/rad_deg
          ELSE IF (this_1%element_type == "cometary") THEN
             WRITE(stdout,"(6(F14.10,1X))") elements(1:2), &
                  elements(3:5)/rad_deg, elements(6)
          ELSE
             WRITE(stdout,"(6(F14.10,1X))") elements
          END IF
       END IF
       IF (lt_corr_) THEN
          ! Light-time correction:
          t_ = getLightTimeCorrectedTime(this_1, observer_)
          IF (error) THEN
             CALL errorMessage("Orbit / getEphemeris (multiple)", &
                  "TRACE BACK 10", 1)
             DEALLOCATE(this_arr_, stat=err)
             DEALLOCATE(ephemeris, stat=err)
             DEALLOCATE(jacobian_prop_arr_, stat=err)
             IF (PRESENT(partials_arr)) THEN
                DEALLOCATE(partials_arr, stat=err)
             END IF
             IF (PRESENT(this_prop_arr)) THEN
                DEALLOCATE(this_prop_arr, stat=err)
             END IF
             IF (PRESENT(this_lt_corr_arr)) THEN
                DEALLOCATE(this_lt_corr_arr, stat=err)
             END IF
             IF (PRESENT(jacobian_lt_corr_arr)) THEN
                DEALLOCATE(jacobian_lt_corr_arr, stat=err)
             END IF
             IF (PRESENT(jacobian_prop_arr)) THEN
                DEALLOCATE(jacobian_prop_arr, stat=err)
             END IF
             CALL NULLIFY(observer_)
             CALL NULLIFY(t_observer)
             CALL NULLIFY(this_1)
             CALL NULLIFY(t_)
             RETURN
          END IF
          IF (info_verb >= 5) THEN
             WRITE(stdout,"(1X,2(1X,A),1X,I0,A,I0,1X,A)") &
                  "Orbit / getEphemeris_multiple:", "Propagate orbit:", i, &
                  "/", nthis, "due to light-time correction."
          END IF
          IF (PRESENT(partials_arr) .OR. PRESENT(jacobian_lt_corr_arr)) THEN
             CALL propagate(this_1, t_, jacobian=jacobian_lt_corr)
             frame = this_1%frame
          ELSE
             CALL propagate(this_1, t_)
          END IF
          IF (error) THEN
             CALL errorMessage("Orbit / getEphemeris (multiple)", &
                  "TRACE BACK 15", 1)
             DEALLOCATE(this_arr_, stat=err)
             DEALLOCATE(ephemeris, stat=err)
             DEALLOCATE(jacobian_prop_arr_, stat=err)
             IF (PRESENT(partials_arr)) THEN
                DEALLOCATE(partials_arr, stat=err)
             END IF
             IF (PRESENT(this_prop_arr)) THEN
                DEALLOCATE(this_prop_arr, stat=err)
             END IF
             IF (PRESENT(this_lt_corr_arr)) THEN
                DEALLOCATE(this_lt_corr_arr, stat=err)
             END IF
             IF (PRESENT(jacobian_lt_corr_arr)) THEN
                DEALLOCATE(jacobian_lt_corr_arr, stat=err)
             END IF
             IF (PRESENT(jacobian_prop_arr)) THEN
                DEALLOCATE(jacobian_prop_arr, stat=err)
             END IF
             CALL NULLIFY(observer_)
             CALL NULLIFY(t_observer)
             CALL NULLIFY(this_1)
             CALL NULLIFY(t_)
             RETURN
          END IF
          IF (PRESENT(jacobian_lt_corr_arr)) THEN
             jacobian_lt_corr_arr(i,:,:) = jacobian_lt_corr
          END IF
       END IF
       this_2 = copy(this_1)
       IF (lt_corr_) THEN
          CALL NULLIFY(this_1%t)
          this_1%t = copy(t_observer)
       END IF
       ! Optionally, return orbit corresponding to
       ! the ephemeris coordinates:
       IF (PRESENT(this_lt_corr_arr)) THEN
          this_lt_corr_arr(i) = copy(this_1)
          IF (error) THEN
             CALL errorMessage("Orbit / getEphemeris (multiple)", &
                  "TRACE BACK 20", 1)
             DEALLOCATE(this_arr_, stat=err)
             DEALLOCATE(ephemeris, stat=err)
             DEALLOCATE(jacobian_prop_arr_, stat=err)
             IF (PRESENT(partials_arr)) THEN
                DEALLOCATE(partials_arr, stat=err)
             END IF
             IF (PRESENT(this_prop_arr)) THEN
                DEALLOCATE(this_prop_arr, stat=err)
             END IF
             IF (PRESENT(this_lt_corr_arr)) THEN
                DEALLOCATE(this_lt_corr_arr, stat=err)
             END IF
             IF (PRESENT(jacobian_lt_corr_arr)) THEN
                DEALLOCATE(jacobian_lt_corr_arr, stat=err)
             END IF
             IF (PRESENT(jacobian_prop_arr)) THEN
                DEALLOCATE(jacobian_prop_arr, stat=err)
             END IF
             CALL NULLIFY(observer_)
             CALL NULLIFY(t_observer)
             CALL NULLIFY(this_1)
             CALL NULLIFY(t_)
             CALL NULLIFY(this_2)
             RETURN
          END IF
       END IF
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,1X,5A)") &
               "Orbit / getEphemeris_multiple:", &
               "Orbital elements (" // &
               TRIM(this_1%element_type) // ", " // &
               TRIM(this_1%frame) // &
               ") at observation date with light-time correction:"
          elements = getElements(this_1, this_1%element_type, this_1%frame)
          IF (this_1%element_type == "keplerian") THEN
             WRITE(stdout,"(6(F14.10,1X))") elements(1:2), &
                  elements(3:6)/rad_deg
          ELSE IF (this_1%element_type == "cometary") THEN
             WRITE(stdout,"(6(F14.10,1X))") elements(1:2), &
                  elements(3:5)/rad_deg, elements(6)
          ELSE
             WRITE(stdout,"(6(F14.10,1X))") elements
          END IF
       END IF

       ! Transform to equatorial topocentric coordinates:
       CALL toCartesian(this_1, frame="equatorial")
       CALL rotateToEquatorial(observer_)
       observer_coordinates(1:6) = getCoordinates(observer_)
       this_1%elements(1:6) = this_1%elements(1:6) - &
            observer_coordinates(1:6)
       ephemeris(i) = getSCoord(this_1, frame="equatorial")
       IF (error) THEN
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "TRACE BACK 25", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          DEALLOCATE(jacobian_prop_arr_, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          IF (PRESENT(this_prop_arr)) THEN
             DEALLOCATE(this_prop_arr, stat=err)
          END IF
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          IF (PRESENT(jacobian_lt_corr_arr)) THEN
             DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          END IF
          IF (PRESENT(jacobian_prop_arr)) THEN
             DEALLOCATE(jacobian_prop_arr, stat=err)
          END IF
          CALL NULLIFY(observer_)
          CALL NULLIFY(t_observer)
          CALL NULLIFY(this_1)
          CALL NULLIFY(t_)
          CALL NULLIFY(this_2)
          RETURN
       END IF

       ! Partial derivatives:
       IF (PRESENT(partials_arr)) THEN

          ! Chain rule.
          jacobian = MATMUL(jacobian_lt_corr, jacobian_prop_arr_(i,:,:))

          SELECT CASE (this_arr_(i)%element_type)

          CASE ("cartesian")

             ! Topocentric, spherical, equatorial coordinates wrt Cartesian elements.

!!$             ! There may be something wrong with the jacobians since
!!$             ! the analytical and numerical approaches do not
!!$             ! match. Or it may be that the numerical solution is
!!$             ! inaccurate because of multiplying matrices with large
!!$             ! condition numbers.
!!$             write(stdout,*)
!!$             write(stdout,*) "Analytical method"
!!$             write(stdout,*)
!!$             CALL partialsSCoordWrtCartesian_d(this_1, scoord_partials)
!!$             call matrix_print(scoord_partials, stdout, errstr)
!!$             write(stdout,*)
!!$             write(stdout,*) cond_nr(scoord_partials, errstr)
!!$             write(stdout,*)
!!$             write(stdout,*)
!!$             write(stdout,*) "Inverse or numerical method"
!!$             write(stdout,*)
!!$             CALL partialsSCoordWrtCartesian_i(this_2, observer_, scoord_partials)
!!$             call matrix_print(scoord_partials, stdout, errstr)
!!$             write(stdout,*)
!!$             write(stdout,*) cond_nr(scoord_partials, errstr)
!!$             write(stdout,*)


             !CALL partialsSCoordWrtCartesian_i(this_2, observer_, scoord_partials)
             CALL partialsSCoordWrtCartesian_d(this_1, scoord_partials)

             IF (error) THEN
                CALL errorMessage("Orbit / getEphemeris (multiple)", &
                     "TRACE BACK 30", 1)
                DEALLOCATE(this_arr_, stat=err)
                DEALLOCATE(ephemeris, stat=err)
                DEALLOCATE(jacobian_prop_arr_, stat=err)
                IF (PRESENT(partials_arr)) THEN
                   DEALLOCATE(partials_arr, stat=err)
                END IF
                IF (PRESENT(this_prop_arr)) THEN
                   DEALLOCATE(this_prop_arr, stat=err)
                END IF
                IF (PRESENT(this_lt_corr_arr)) THEN
                   DEALLOCATE(this_lt_corr_arr, stat=err)
                END IF
                IF (PRESENT(jacobian_lt_corr_arr)) THEN
                   DEALLOCATE(jacobian_lt_corr_arr, stat=err)
                END IF
                IF (PRESENT(jacobian_prop_arr)) THEN
                   DEALLOCATE(jacobian_prop_arr, stat=err)
                END IF
                CALL NULLIFY(observer_)
                CALL NULLIFY(t_observer)
                CALL NULLIFY(this_1)
                CALL NULLIFY(t_)
                CALL NULLIFY(this_2)
                RETURN
             END IF

          CASE ("cometary")

             ! Topocentric equatorial coordinates wrt cometary orbital elements.
             CALL partialsSCoordWrtCometary(this_2, observer_, scoord_partials)
             IF (error) THEN
                CALL errorMessage("Orbit / getEphemeris (multiple)", &
                     "TRACE BACK 35", 1)
                DEALLOCATE(this_arr_, stat=err)
                DEALLOCATE(ephemeris, stat=err)
                DEALLOCATE(jacobian_prop_arr_, stat=err)
                IF (PRESENT(partials_arr)) THEN
                   DEALLOCATE(partials_arr, stat=err)
                END IF
                IF (PRESENT(this_prop_arr)) THEN
                   DEALLOCATE(this_prop_arr, stat=err)
                END IF
                IF (PRESENT(this_lt_corr_arr)) THEN
                   DEALLOCATE(this_lt_corr_arr, stat=err)
                END IF
                IF (PRESENT(jacobian_lt_corr_arr)) THEN
                   DEALLOCATE(jacobian_lt_corr_arr, stat=err)
                END IF
                IF (PRESENT(jacobian_prop_arr)) THEN
                   DEALLOCATE(jacobian_prop_arr, stat=err)
                END IF
                CALL NULLIFY(observer_)
                CALL NULLIFY(t_observer)
                CALL NULLIFY(this_1)
                CALL NULLIFY(t_)
                CALL NULLIFY(this_2)
                RETURN
             END IF

          CASE ("keplerian")

             ! Topocentric equatorial coordinates wrt Keplerian orbital elements.
             CALL partialsSCoordWrtKeplerian(this_2, observer_, scoord_partials)
             IF (error) THEN
                CALL errorMessage("Orbit / getEphemeris (multiple)", &
                     "TRACE BACK 35", 1)
                DEALLOCATE(this_arr_, stat=err)
                DEALLOCATE(ephemeris, stat=err)
                DEALLOCATE(jacobian_prop_arr_, stat=err)
                IF (PRESENT(partials_arr)) THEN
                   DEALLOCATE(partials_arr, stat=err)
                END IF
                IF (PRESENT(this_prop_arr)) THEN
                   DEALLOCATE(this_prop_arr, stat=err)
                END IF
                IF (PRESENT(this_lt_corr_arr)) THEN
                   DEALLOCATE(this_lt_corr_arr, stat=err)
                END IF
                IF (PRESENT(jacobian_lt_corr_arr)) THEN
                   DEALLOCATE(jacobian_lt_corr_arr, stat=err)
                END IF
                IF (PRESENT(jacobian_prop_arr)) THEN
                   DEALLOCATE(jacobian_prop_arr, stat=err)
                END IF
                CALL NULLIFY(observer_)
                CALL NULLIFY(t_observer)
                CALL NULLIFY(this_1)
                CALL NULLIFY(t_)
                CALL NULLIFY(this_2)
                RETURN
             END IF

          CASE default

             error = .TRUE.
             CALL errorMessage("Orbit / getEphemeris (multiple)", &
                  "Could not choose between element types" // &
                  "('cartesian', 'cometary', or 'keplerian'): " // &
                  TRIM(this_arr_(i)%element_type), 1)
             RETURN

          END SELECT

          ! Chain rule.
          partials_arr(i,:,:) = MATMUL(scoord_partials, jacobian)

       END IF

       CALL NULLIFY(this_1)
       CALL NULLIFY(this_2)
       CALL NULLIFY(t_)

    END DO

    CALL NULLIFY(observer_)
    CALL NULLIFY(t_observer)
    IF (ASSOCIATED(jacobian_prop_arr_)) THEN
       DEALLOCATE(jacobian_prop_arr_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getEphemeris (multiple)", &
               "Could not deallocate memory (5).", 1)
          DEALLOCATE(this_arr_, stat=err)
          DEALLOCATE(ephemeris, stat=err)
          IF (PRESENT(partials_arr)) THEN
             DEALLOCATE(partials_arr, stat=err)
          END IF
          IF (PRESENT(this_prop_arr)) THEN
             DEALLOCATE(this_prop_arr, stat=err)
          END IF
          IF (PRESENT(this_lt_corr_arr)) THEN
             DEALLOCATE(this_lt_corr_arr, stat=err)
          END IF
          IF (PRESENT(jacobian_lt_corr_arr)) THEN
             DEALLOCATE(jacobian_lt_corr_arr, stat=err)
          END IF
          IF (PRESENT(jacobian_prop_arr)) THEN
             DEALLOCATE(jacobian_prop_arr, stat=err)
          END IF
          RETURN
       END IF
    END IF
    DO i=1,SIZE(this_arr_)
       CALL NULLIFY(this_arr_(i))
    END DO
    DEALLOCATE(this_arr_, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getEphemeris (multiple)", &
            "Could not deallocate memory (15).", 1)
       DEALLOCATE(ephemeris, stat=err)
       IF (PRESENT(partials_arr)) THEN
          DEALLOCATE(partials_arr, stat=err)
       END IF
       IF (PRESENT(this_prop_arr)) THEN
          DEALLOCATE(this_prop_arr, stat=err)
       END IF
       IF (PRESENT(this_lt_corr_arr)) THEN
          DEALLOCATE(this_lt_corr_arr, stat=err)
       END IF
       IF (PRESENT(jacobian_lt_corr_arr)) THEN
          DEALLOCATE(jacobian_lt_corr_arr, stat=err)
       END IF
       IF (PRESENT(jacobian_prop_arr)) THEN
          DEALLOCATE(jacobian_prop_arr, stat=err)
       END IF
       RETURN
    END IF

  END SUBROUTINE getEphemeris_Orb_multiple





  !! *Description*:
  !!
  !! Returns the coordinate frame.
  !!
  !! Returns error.
  !!
  CHARACTER(len=FRAME_LEN) FUNCTION getFrame_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getFrame", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getFrame_Orb = TRIM(this%frame)

  END FUNCTION getFrame_Orb





  !! *Description*:
  !!
  !! Calculates the Keplerian orbital elements (for elliptical orbits)
  !! from the heliocentric (ecliptical or equatorial) Cartesian orbital
  !! elements or elliptic cometary elements for the current epoch. 
  !! Hyperbolic cometary elements lead to an error.
  !!
  !! Returns error.
  !! 
  FUNCTION getKeplerianElements(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(6)   :: getKeplerianElements

    TYPE (CartesianCoordinates) :: ccoord
    REAL(bp), DIMENSION(:,:), POINTER :: planeph => NULL()
    REAL(bp), DIMENSION(6)   :: elements, coordinates
    REAL(bp) :: mjd

    TYPE (Orbit):: this_
    TYPE (Time) :: t
    REAL(bp), DIMENSION(3) :: pos, vel, k, fb, gb, evec, cos_angles
    REAL(bp) :: r, ru, ea, e_cos_ea, e_sin_ea, cos_ea, sin_ea, &
         alpha, gamma, mjd_tt, a, e, i, an, ap, ma, div, tmp1, tmp2, &
         varpi

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getKeplerianElements", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    ! Return initial elements if error occurs:
    getKeplerianElements(1:6) = this%elements(1:6)

    SELECT CASE (TRIM(this%element_type))

    CASE ("keplerian")

       RETURN

    CASE ("cartesian")

       this_ = copy(this)
       CALL rotateToEcliptic(this_)

       ! Semimajor axis, eccentricity and mean anomaly:
       pos = this_%elements(1:3)
       r   = SQRT(DOT_PRODUCT(pos,pos))
       vel = this_%elements(4:6)
       ! Angular momentum:
       k = cross_product(pos,vel)
       ! Eccentricity vector:
       evec = cross_product(vel,k)/planetary_mu(this_%center) - pos/r

       ! alpha = rv^2/mu
       alpha = r*DOT_PRODUCT(vel,vel)/planetary_mu(this%center)
       ! h = v^2/2 - mu/r and for elliptical orbits a = -mu/2h
       ! -> a = r / (2 - rv^2/mu)
       IF (ABS(2.0_bp - alpha) < 10.0_bp*EPSILON(alpha)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getKeplerianElements", &
               "a is approaching infinity.", 1)
          RETURN
       END IF
       a = r / (2.0_bp - alpha)
       ! gamma = sqrt(mu*a) 
       ! Note that we use |a| below in order to allow for hyperbolic
       ! orbits.
       gamma = SQRT(planetary_mu(this%center)*ABS(a))
       IF (ABS(gamma) < TINY(gamma)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getKeplerianElements", &
               "Gamma is computationally too small.", 1)
          RETURN
       END IF
       IF (a > 0.0_bp) THEN ! elliptic motion

          ! r(1+e*cos(f))=a|1-e^2|, r*cos(f)=a*cos(E)-ae, and for elliptical 
          ! orbits a=r/(2-rv^2/mu) (see above) and |1-e^2|=1-e^2
          ! -> e*cos(E)=rv^2/mu-1=alpha-1
          e_cos_ea = alpha - 1.0_bp
          !
          e_sin_ea = DOT_PRODUCT(pos,vel)/gamma
          e = SQRT(e_cos_ea**2.0_bp + e_sin_ea**2.0_bp)
          IF (e == 1.0_bp) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getKeplerianElements", &
                  "Orbit is parabolic.", 1)
             RETURN
          ELSE IF (e < EPSILON(e)) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getKeplerianElements", &
                  "Orbit is (almost) circular.", 1)
             WRITE(stderr,*) e
             RETURN
          END IF
          cos_ea = e_cos_ea/e
          sin_ea = e_sin_ea/e
          IF (ABS(cos_ea) > 1.0_bp) cos_ea = SIGN(1.0_bp, cos_ea)
          ea = ACOS(cos_ea)
          IF (sin_ea < 0.0_bp) ea = two_pi - ea
          ma = ea - e_sin_ea
          IF (ma == two_pi) THEN
             ma = 0.0_bp
          END IF
          ma = MODULO(ma,two_pi)

       ELSE IF (a < 0.0_bp) THEN ! hyperbolic motion

          error = .TRUE.
          CALL errorMessage("Orbit / getKeplerianElements", &
               "Orbit is hyperbolic.", 1)
          RETURN

       END IF

       ! Inclination and ascending node:
       ru = SQRT(DOT_PRODUCT(k,k))
       IF (ru < EPSILON(ru)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getKeplerianElements", &   
               "Position and velocity almost parallel.", 1)
          RETURN
       END IF
       k = k/ru
       cos_angles(1) = k(3) ! = cos(i)
       IF (ABS(cos_angles(1)) > 1.0_bp) THEN
          cos_angles(1) = SIGN(1.0_bp, cos_angles(1))
       END IF
       i = ACOS(cos_angles(1))
       ! Notice computational abs!!
       ru = SQRT(ABS(1.0_bp - k(3)**2.0_bp)) ! = sin(i)
       IF (ru < EPSILON(ru)) THEN
          CALL errorMessage("Orbit / getKeplerianElements", &
               "Warning: Inclination (almost) zero, setting Longitude of Node to 0 deg.", 1)
          an = 0.0_bp
       ELSE
          cos_angles(2) = -k(2)/ru
          IF (ABS(cos_angles(2)) > 1.0_bp) THEN
             cos_angles(2) = SIGN(1.0_bp, cos_angles(2))
          END IF
          an = ACOS(cos_angles(2))
          IF (k(1) < 0.0_bp) THEN
             an = two_pi - an
          END IF
       END IF
       IF (an == two_pi) THEN
          an = 0.0_bp
       END IF

       ! Argument of periapsis:
       div = 1.0_bp + k(3)
       fb(1) = 1.0_bp - k(1)**2.0_bp/div
       fb(2) = -k(1)*k(2)/div
       fb(3) = -k(1)
       gb = cross_product(k, fb)
       tmp1 = DOT_PRODUCT(evec, gb)
       tmp2 = DOT_PRODUCT(evec, fb)
       varpi = ATAN2(tmp1, tmp2)
       ap = varpi - an
       ap = MODULO(ap,two_pi)
       IF (ap == two_pi) THEN
          ap = 0.0_bp
       END IF

       getKeplerianElements = (/ a, e, i, an, ap, ma /)
       CALL NULLIFY(this_)

    CASE ("cometary")

       IF (this%elements(2) >= 1.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getKeplerianElements", &
               "Hyperbolic orbit; cannot return Keplerian elements.", 1)
          RETURN
       END IF

       ! Semimajor axis:
       getKeplerianElements(1) = &
            this%elements(1) / (1.0_bp - this%elements(2))
       ! Mean anomaly:
       t = copy(this%t)
       mjd_tt = getMJD(t, "TT")
       CALL NULLIFY(t)
       getKeplerianElements(6) = &
            SQRT(planetary_mu(this%center) / &
            getKeplerianElements(1)**3.0_bp) * &
            (mjd_tt - this%elements(6))
       getKeplerianElements(6) = &
            MODULO(getKeplerianElements(6), two_pi)

    CASE ("cometary_ta", "cometary_ma")

       this_ = copy(this)
       IF (error) THEN
          CALL errorMessage("Orbit / getKeplerianElements", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       CALL toCometary(this_)
       IF (error) THEN
          CALL errorMessage("Orbit / getKeplerianElements", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       getKeplerianElements = getElements(this_, "keplerian")
       IF (error) THEN
          CALL errorMessage("Orbit / getKeplerianElements", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL NULLIFY(this_)

    CASE default

       error = .TRUE.
       CALL errorMessage("Orbit / getKeplerianElements", &
            "Conversion from " // TRIM(this%element_type) // &
            " elements to keplerian elements has not yet been implemented.", 1)
       RETURN

    END SELECT

  END FUNCTION getKeplerianElements





  !! *Description*:
  !!
  !! Calculates the k-vector.
  !!
  !! Returns error.
  !! 
  FUNCTION getKVector(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(3)   :: getKVector
    REAL(bp), DIMENSION(6)   :: elements

    IF (.NOT.this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getKVector", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this, "cartesian", "ecliptic")
    IF (error) THEN
       CALL errorMessage("Orbit / getKVector", &
            "TRACE BACK", 1)
       RETURN
    END IF

    getKVector = cross_product(elements(1:3),elements(4:6))

  END FUNCTION getKVector





  !! *Description*:
  !!
  !! Computes the light-time corrected epoch via linear extrapolation
  !! from the objects orbital elements (this) to the retarded,
  !! apparent position seen from the observer's coordinates
  !! (observer) at the epoch of the object's orbital elements
  !! (this).
  !!
  !! Returns error.
  !!
  FUNCTION getLightTimeCorrectedTime(this, observer)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                :: this
    TYPE (CartesianCoordinates), INTENT(in) :: observer
    TYPE (Time)                             :: getLightTimeCorrectedTime
    TYPE (Orbit)                            :: this_
    TYPE (CartesianCoordinates)             :: observer_
    REAL(bp), DIMENSION(6)                  :: observer_coord
    REAL(bp)                                :: lt, rr, rv, vv, mjd_tt

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getLightTimeCorrectedTime", &
            "Object object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getLightTimeCorrectedTime", &
            "Observer object has not yet been initialized.", 1)
       RETURN
    END IF

    this_ = copy(this)
    observer_ = copy(observer)
    mjd_tt = getMJD(this_%t, "TT")
    CALL toCartesian(this_, frame="equatorial")
    CALL rotateToEquatorial(observer_)
    observer_coord = getCoordinates(observer_)
    rr = DOT_PRODUCT((this_%elements(1:3) - observer_coord(1:3)), &
         (this_%elements(1:3) - observer_coord(1:3)))
    rv = DOT_PRODUCT((this_%elements(1:3) - observer_coord(1:3)), &
         this_%elements(4:6))
    vv = DOT_PRODUCT(this_%elements(4:6),this_%elements(4:6))
    IF (rv**2.0_bp+rr*(sol**2.0_bp-vv) < 0.0_bp) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getLightTimeCorrectedTime", &
            "Object is moving faster than light.", 1)
       RETURN
    END IF
    lt = (SQRT(rv**2.0_bp+rr*(sol**2.0_bp-vv))-rv) / (sol**2.0_bp-vv)
    CALL NULLIFY(this_)
    CALL NULLIFY(observer_)
    CALL NEW(getLightTimeCorrectedTime, mjd_tt-lt, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / getLightTimeCorrectedTime", &
            "TRACE BACK 1", 1)
       RETURN
    END IF

  END FUNCTION getLightTimeCorrectedTime





  !! *Description*:
  !!
  !! The Sitarski algorithm (Acta Astronomica, Vol.18, No.2, p.171)
  !! for computing the minimum distance between two Keplerian orbits
  !! (MOID = minimum orbital intersection distance).
  !!
  REAL(bp) FUNCTION getMOID(this, that)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this, that
    TYPE (Orbit) :: this_, that_
    TYPE (Time) :: t
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: ta2
    REAL(bp), DIMENSION(:), ALLOCATABLE :: ta1
    REAL(bp), DIMENSION(10,10) :: local_minima
    REAL(bp), DIMENSION(3,3) :: rot1, rot2
    REAL(bp), DIMENSION(2,2) :: rot
    REAL(bp), DIMENSION(6) :: elem1, elem2
    REAL(bp), DIMENSION(3) :: coord1, coord2, coord11, coord22, &
         rad1, rad2, dr, der2
    REAL(bp), DIMENSION(2) :: sine, cosine, f1, f2, ff
    REAL(bp) :: par1, par2, par11, par22, r1, r2, ta0, dta, q1, q2, &
         f, ea, sin_ea, cos_ea, ma
    INTEGER :: i, j, k, N, m1, m2, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getMOID", &
            "Object object has not yet been initialized (1).", 1)
       RETURN
    END IF

    IF (.NOT. that%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getMOID", &
            "Object object has not yet been initialized (2).", 1)
       RETURN
    END IF

    this_ = copy(this)
    that_ = copy(that)

    ! Propagate orbital elements of "that" to 
    ! the epoch of "this" if needed:
    IF (.NOT.equal(getTime(this_),getTime(that_))) THEN
       t = getTime(this_)
       CALL propagate(that_, t)
       IF (error) THEN
          CALL errorMessage("Orbit / getMOID", &
               "TRACE TRACK (5).", 1)
          RETURN
       END IF
    END IF

    ! Convert to Keplerian elements:
    CALL toKeplerian(this_)
    CALL toKeplerian(that_)
    IF (error) THEN
       CALL errorMessage("Orbit / getMOID", &
            "TRACE TRACK (10).", 1)
       RETURN
    END IF

    ! Rotation matrix of the first orbit (polar->ecliptic):
    rot1 = getTransformationMatrix(this_)
    ! Rotation matrix of the second orbit (--- " ---):
    rot2 = getTransformationMatrix(that_)

    rot(1,1) = rot1(1,1)*rot2(1,1) + rot1(2,1)*rot2(2,1) + rot1(3,1)*rot2(3,1)
    rot(1,2) = rot1(1,1)*rot2(1,2) + rot1(2,1)*rot2(2,2) + rot1(3,1)*rot2(3,2)
    rot(2,1) = rot1(1,2)*rot2(1,1) + rot1(2,2)*rot2(2,1) + rot1(3,2)*rot2(3,1)
    rot(2,2) = rot1(1,2)*rot2(1,2) + rot1(2,2)*rot2(2,2) + rot1(3,2)*rot2(3,2)

    ! Orbital parameters:
    elem1 = getElements(this_, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / getMOID", &
            "TRACE TRACK (15).", 1)
       RETURN
    END IF
    elem2 = getElements(that_, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / getMOID", &
            "TRACE TRACK (20).", 1)
       RETURN
    END IF
    IF (info_verb >= 3) THEN
       WRITE(stdout,"(A)") "Orbit #1:"
       WRITE(stdout,"(6(1X,F20.15))") elem1(1:2), elem1(3:6)/rad_deg
       WRITE(stdout,"(A)") "Orbit #2:"
       WRITE(stdout,"(6(1X,F20.15))") elem2(1:2), elem2(3:6)/rad_deg
    END IF

    par1 = (1-elem1(2)**2)*elem1(1)
    par2 = (1-elem2(2)**2)*elem2(1)
    q1 = elem1(1)*(1.0_bp - elem1(2))
    q2 = elem2(1)*(1.0_bp - elem2(2))

    k = 0
    N = 360
    local_minima(:,1) = HUGE(local_minima(:,1))

    !-----------------------------------------------------------
    ! Find minimum, increase true anomaly sampling if necessary:
    !-----------------------------------------------------------

    DO WHILE (N <= 360000 .AND. &
         (k == 0 .OR. MINVAL(local_minima(:,1)) > q2 + q1))

       ALLOCATE(ta1(N), ta2(2,N), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getMOID", &
               "Could not allocate memory.", 1)
          RETURN
       END IF

       dta = two_pi/N
       ta0 = 0.0_bp

       ! Moving along the first orbit with steps of dta:
       ! (N=360, dta=1 deg)
       DO i=1,N
          ta1(i) = ta0 + (i-1)*dta
          CALL SOLVE(elem1, ta1(i), elem2, rot, sine, cosine)
          ta2(1,i) = angle(cosine(1), sine(1))
          ta2(2,i) = angle(cosine(2), sine(2))
       END DO

       ! Sign-changes-intervals:
       ! Find intervals of true anomalies (F,f) where the second
       ! derivative of the function to be minimized (see Sitarski) 
       ! changes its sign.

       DO i=1,N
          DO j=1,2

             ! Set (F,f) intervals
             IF (i == 1) THEN
                f1(1) = ta1(N)
                f1(2) = ta1(i)
                f2(1) = ta2(j,N)
                f2(2) = ta2(j,i)
             ELSE
                f1(1) = ta1(i-1)
                f1(2) = ta1(i)
                f2(1) = ta2(j,i-1)
                f2(2) = ta2(j,i)
             END IF

             ! Check that upper interval limits are less than 0 deg,
             ! ie 359.9999... deg:
             IF (f1(2) == 0.0_bp) THEN
                f1(2) = two_pi - 100*EPSILON(f1(2))
             END IF
             IF (f2(2) == 0.0_bp) THEN
                f2(2) = two_pi - 100*EPSILON(f2(2))
             END IF

             IF (info_verb >= 3) THEN
                WRITE(stdout,"(A)") "(F,f) intervals:"
                WRITE(stdout,"(4(1X,F10.6))") f1/rad_deg, f2/rad_deg
             END IF

             par11 = (1.0_bp - elem1(2)**2)*elem1(1)
             r1 = par11/(1.0_bp + elem1(2)*COS(f1(1)))
             coord1(1) = r1*COS(f1(1))
             coord1(2) = r1*SIN(f1(1))
             coord1(3) = 0.0_bp

             par22 = (1.0_bp - elem2(2)**2)*elem2(1)
             r2 = par22/(1.0_bp + elem2(2)*COS(f2(1)))
             coord2(1) = r2*COS(f2(1))
             coord2(2) = r2*SIN(f2(1))
             coord2(3) = 0.0_bp

             m1 = NINT(SIGN(1.0_bp, derivative(coord1, coord2, r2, rot, elem2)))

             par11 = (1.0_bp - elem1(2)**2)*elem1(1)
             r1 = par11/(1.0_bp + elem1(2)*COS(f1(2)))
             coord1(1) = r1*COS(f1(2))
             coord1(2) = r1*SIN(f1(2))
             coord1(3) = 0.0_bp

             par22 = (1.0_bp - elem2(2)**2)*elem2(1)
             r2 = par22/(1.0_bp + elem2(2)*COS(f2(2)))
             coord2(1) = r2*COS(f2(2))
             coord2(2) = r2*SIN(f2(2))
             coord2(3) = 0.0_bp

             m2 = NINT(SIGN(1.0_bp, derivative(coord1, coord2, r2, rot, elem2)))

             IF (m1*m2 < 0) THEN

                ! Refine the initial interval by successive halving

                CALL HALF(j, elem1, f1, elem2, f2, rot, m1, m2, ff)

                par11 = (1.0_bp - elem1(2)**2)*elem1(1)
                r1 = par11/(1.0_bp + elem1(2)*COS(ff(1)))
                coord11(1) = r1*COS(ff(1))
                coord11(2) = r1*SIN(ff(1))
                coord11(3) = 0.0_bp

                par22 = (1.0_bp - elem2(2)**2)*elem2(1)
                r2 = par22/(1.0_bp + elem2(2)*COS(ff(2)))
                coord22(1) = r2*COS(ff(2))
                coord22(2) = r2*SIN(ff(2))
                coord22(3) = 0.0_bp

                ! Study the found extrema, compute second derivatives:
                der2(1) = (r1/par1)*(elem1(2)*r1**2/par1*(elem1(2)*r1+coord11(1)) + &
                     coord11(1)*(rot(1,1)*coord22(1)+rot(1,2)*coord22(2)) + &
                     coord11(2)*(rot(2,1)*coord22(1)+rot(2,2)*coord22(2)))

                der2(2) = (r2/par2)*(elem2(2)*r2**2/par2*(elem2(2)*r2+coord22(1)) + &
                     coord22(1)*(rot(1,1)*coord11(1)+rot(2,1)*coord11(2)) + &
                     coord22(2)*(rot(1,2)*coord11(1)+rot(2,2)*coord11(2)))

                der2(3) = (r1*r2/par1/par2)*((elem2(2)*r2+coord22(1))*(rot(2,2)* &
                     (elem1(2)*r1+coord11(1)) - rot(1,2)*coord11(2)) - &
                     coord22(2)*(rot(2,1)*(elem1(2)*r1+coord11(1)) - &
                     rot(1,1)*coord11(2)))

                ! Check for local minimum:
                IF (der2(1) > 0.0_bp .AND. der2(2) > 0.0_bp .AND.&
                     (der2(1)*der2(2)-der2(3)**2) > 0.0_bp) THEN
                   k = k + 1
                   IF (k > 10) THEN
                      error = .TRUE.
                      CALL errorMessage("Orbit / getMOID", &
                           "More than 10 minima in the MOID function.", 1)
                      RETURN
                   END IF

                   ! Compute distance in local minimum:
                   rad1 = MATMUL(rot1, coord11)
                   rad2 = MATMUL(rot2, coord22)
                   dr = rad1 - rad2
                   IF (i == 1) THEN
                      local_minima(k,1:6) = (/ &
                           SQRT(SUM(dr**2)), &
                           REAL(N,bp), &
                           ta1(N), &
                           ta1(i), &
                           ta2(j,N), &
                           ta2(j,i) /)
                   ELSE
                      local_minima(k,1:6) = (/ &
                           SQRT(SUM(dr**2)), &
                           REAL(N,bp), &
                           ta1(i-1), &
                           ta1(i), &
                           ta2(j,i-1), &
                           ta2(j,i) /)
                   END IF
                   local_minima(k,7:10) = (/ &
                        f1(1), &
                        f1(2), &
                        f2(1), &
                        f2(2) /)
                END IF
             END IF
          END DO
       END DO

       DEALLOCATE(ta1, ta2, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getMOID", &
               "Could not deallocate memory.", 1)
          DEALLOCATE(ta1, stat=err)
          DEALLOCATE(ta2, stat=err)
          RETURN
       END IF

       ! Increase the true anomaly sampling in case no minimum was found:
       N = 10*N

    END DO

    !----------------------
    ! Global minimum, MOID:
    !----------------------
    IF (k /= 0) THEN
       getMOID = MINVAL(local_minima(1:k,1),1)
       IF (info_verb >= 2) THEN
          i = MINLOC(local_minima(1:k,1),1)
          WRITE(stdout,"(10(F13.9,1X))") local_minima(i,1:2), &
               local_minima(i,3:10)/rad_deg
          f = 0.5_bp*SUM(local_minima(i,7:8))
          cos_ea = (elem1(2) + COS(f))/(1.0_bp + elem1(2)*COS(f))
          sin_ea = (SQRT(1.0_bp - elem1(2)**2)*SIN(f))/(1.0_bp + elem1(2)*COS(f))
          ea = ATAN2(sin_ea, cos_ea)
          ma = MODULO(ea - elem1(2)*SIN(ea), two_pi)
          WRITE(stdout,"(A,6(1X,F15.10))") &
               "Keplerian elements for the 1st orbit at the epoch:", &
               elem1(1:2), elem1(3:6)/rad_deg
          WRITE(stdout,"(A,1X,F15.10)") &
               "Mean anomaly M for the 1st orbit when MOID takes place:", &
               ma/rad_deg
          f = 0.5_bp*SUM(local_minima(i,9:10))
          cos_ea = (elem2(2) + COS(f))/(1.0_bp + elem2(2)*COS(f))
          sin_ea = (SQRT(1.0_bp - elem2(2)**2)*SIN(f))/(1.0_bp + elem2(2)*COS(f))
          ea = ATAN2(sin_ea, cos_ea)
          ma = MODULO(ea - elem2(2)*SIN(ea), two_pi)
          WRITE(stdout,"(A,6(1X,F15.10))") &
               "Keplerian elements for the 2nd orbit at the epoch of the 1st orbit:", &
               elem2(1:2), elem2(3:6)/rad_deg
          WRITE(stdout,"(A,1X,F15.10)") &
               "Mean anomaly M for the 2nd orbit when MOID takes place:", &
               ma/rad_deg
       END IF
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / getMOID", &
            "Minimum not found.", 1)       
       getMOID = HUGE(local_minima(1,1))
    END IF

    CALL NULLIFY(this_)
    CALL NULLIFY(that_)

  END FUNCTION getMOID





  SUBROUTINE solve(elem1, ta1, elem2, rot, v1, v2)

    ! See Sitarski, Eq. (10).

    IMPLICIT NONE
    REAL(bp), DIMENSION(2,2), INTENT(in) :: rot
    REAL(bp), DIMENSION(6), INTENT(in) :: elem1, elem2
    REAL(bp), DIMENSION(2), INTENT(out) :: v1, v2
    REAL(bp), INTENT(in) :: ta1

    COMPLEX(cbp), DIMENSION(2) :: z1, z2, zz
    COMPLEX(cbp) :: s, t, w
    REAL(bp), DIMENSION(3) :: coord1
    REAL(bp) :: r1, par1, par2

    par1 = (1.0_bp - elem1(2)**2)*elem1(1)
    r1 = par1/(1.0_bp + elem1(2)*COS(ta1))
    coord1(1) = r1*COS(ta1)
    coord1(2) = r1*SIN(ta1)
    coord1(3) = 0.0_bp

    par2 = (1-elem2(2)**2)*elem2(1)
    s = elem1(2)*r1*coord1(2)/par2
    t = rot(1,2)*coord1(2)-rot(2,2)*(elem1(2)*r1+coord1(1))
    w = elem2(2)*s+rot(1,1)*coord1(2)-rot(2,1)*(elem1(2)*r1+coord1(1))

    z1(1) = (-t*s+w*SQRT(t**2+w**2-s**2))/(t**2+w**2)
    z1(2) = (-t*s-w*SQRT(t**2+w**2-s**2))/(t**2+w**2)

    zz = CMPLX(z1)
    IF (AIMAG(zz(1)) == 0.0_bp) THEN
       v1(1) = REAL(z1(1),bp)
    END IF
    IF (AIMAG(zz(2)) == 0.0_bp) THEN
       v1(2) = REAL(z1(2),bp)
    END IF

    z2(1) = (-w*s-t*SQRT(t**2+w**2-s**2))/(t**2+w**2)
    z2(2) = (-w*s+t*SQRT(t**2+w**2-s**2))/(t**2+w**2)

    zz = CMPLX(z2)
    IF (AIMAG(zz(1)) == 0.0_bp) THEN
       v2(1) = REAL(z2(1),bp)
    END IF
    IF (AIMAG(zz(2)) == 0.0_bp) THEN
       v2(2) = REAL(z2(2),bp)
    END IF

  END SUBROUTINE solve





  SUBROUTINE half(i, elem1, f1, elem2, f2, rot, m1, m2, df)

    ! Halving of the sign-change-intervals.
    ! Generalize?

    IMPLICIT NONE
    REAL(bp), DIMENSION(2,2) :: rot
    REAL(bp), DIMENSION(6) :: elem1, elem2
    REAL(bp), DIMENSION(3) :: coord1, coord2
    REAL(bp), DIMENSION(2) :: df, f1, f2, sine, cosine
    REAL(bp) :: par1, par2, r1, r2, eps
    INTEGER :: i, m, m1, m2

    eps = 1.e-10_bp
    DO WHILE (ABS(f1(1)-f1(2)) > eps)
       df(1) = 0.5_bp*(f1(1) + f1(2))

       CALL solve(elem1, df(1), elem2, rot, sine, cosine)
       df(2) = angle(cosine(i),sine(i))

       par1 = (1.0_bp - elem1(2)**2)*elem1(1)
       r1 = par1/(1.0_bp + elem1(2)*COS(df(1)))
       coord1(1) = r1*COS(df(1))
       coord1(2) = r1*SIN(df(1))
       coord1(3) = 0.0_bp
       par2 = (1.0_bp - elem2(2)**2)*elem2(1)
       r2 = par2/(1.0_bp + elem2(2)*COS(df(2)))
       coord2(1) = r2*COS(df(2))
       coord2(2) = r2*SIN(df(2))
       coord2(3) = 0.0_bp

       m = NINT(SIGN(1.0_bp, derivative(coord1, coord2, r2, rot, elem2)))
       IF (m == m1) THEN
          f1(1) = df(1)
          f2(1) = df(2)
       ELSE IF (m == m2) THEN
          f1(2) = df(1)
          f2(2) = df(2)
       END IF
    END DO

  END SUBROUTINE half





  !! Partial derivative of the function to be minimized
  !! (df/dv, see Sitarski 1968).
  !!
  REAL(bp) FUNCTION derivative(coord1, coord2, r2, rot, elem2)

    IMPLICIT NONE
    REAL(bp), DIMENSION(2,2) :: rot
    REAL(bp), DIMENSION(6)   :: elem2
    REAL(bp), DIMENSION(3)   :: coord1, coord2
    REAL(bp)                 :: r2, par2

    par2 = (1.0_bp-elem2(2)**2.0_bp)*elem2(1)
    derivative = (r2/par2) * &
         (coord2(2)*(elem2(2)*r2 + rot(1,1)*coord1(1) + rot(2,1)*coord1(2)) - &
         (elem2(2)*r2 + coord2(1)) * (rot(1,2)*coord1(1) + rot(2,2)*coord1(2)))

  END FUNCTION derivative





  !! *Description*:
  !!
  !! Returns the semimajor axis for the input orbit.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getSemimajorAxis(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)       :: this

    TYPE (Orbit) :: this_
    REAL(bp), DIMENSION(3) :: pos, vel
    REAL(bp) :: alpha, r

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getSemimajorAxis", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    SELECT CASE (this%element_type)

    CASE ("cometary")

       getSemimajorAxis = this%elements(1) / &
            (1.0_bp - this%elements(2))

    CASE ("keplerian")

       getSemimajorAxis = this%elements(1)

    CASE ("cartesian")

       this_ = copy(this)
       CALL rotateToEcliptic(this_)
       pos = this_%elements(1:3)
       vel = this_%elements(4:6)
       r = SQRT(DOT_PRODUCT(pos,pos))
       ! alpha = rv^2/mu
       alpha = r*DOT_PRODUCT(vel,vel)/planetary_mu(this_%center)
       ! h = v^2/2 - mu/r and for elliptical orbits a = -mu/2h
       ! -> a = r / (2 - rv^2/mu)
       IF (ABS(2.0_bp - alpha) < 10.0_bp*EPSILON(alpha)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getSemimajorAxis", &
               "a is approaching infinity.", 1)
          RETURN
       END IF
       getSemimajorAxis = r / (2.0_bp - alpha)

    END SELECT

  END FUNCTION getSemimajorAxis





  !! *Description*:
  !!
  !!
  !! Returns error.
  !!
  SUBROUTINE getParameters_Orb(this, mass, dyn_model, &
       integration_step, integrator, finite_diff, &
       additional_perturbers)

    IMPLICIT NONE
    TYPE(Orbit), INTENT(in) :: this
    REAL(bp), INTENT(out), OPTIONAL :: mass
    CHARACTER(len=DYN_MODEL_LEN), INTENT(out), OPTIONAL :: dyn_model
    REAL(bp), INTENT(out), OPTIONAL :: integration_step
    CHARACTER(len=INTEGRATOR_LEN), INTENT(out), OPTIONAL :: integrator
    REAL(bp), DIMENSION(6), INTENT(out), OPTIONAL :: finite_diff
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL :: additional_perturbers

    INTEGER :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getParameters", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(mass)) THEN
       mass = this%mass_prm
    END IF

    IF (PRESENT(dyn_model)) THEN
       dyn_model = this%dyn_model_prm
    END IF
    IF (PRESENT(integration_step)) THEN
       integration_step = this%integration_step_prm
    END IF
    IF (PRESENT(integrator)) THEN
       integrator = this%integrator_prm
    END IF
    IF (PRESENT(finite_diff)) THEN
       finite_diff = this%finite_diff_prm
    END IF
    IF (PRESENT(additional_perturbers)) THEN
       !! Returns the current set of orbital elements for additional
       !! perturbers if available. The shape of the table is (1:n,1:8)
       !! where n is the number of perturbers. The first six elements on a
       !! row (n,1:6) contain the Cartesian equatorial orbital elements,
       !! (n,7) is the epoch given as MJD TT, and (n,8) is the mass of the
       !! perturber given as solar masses.
       IF (ASSOCIATED(this%additional_perturbers)) THEN
          ALLOCATE(additional_perturbers(SIZE(this%additional_perturbers,dim=1), &
               SIZE(this%additional_perturbers,dim=2)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getParameters", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
          additional_perturbers = this%additional_perturbers
       ELSE
          NULLIFY(additional_perturbers)
       END IF
    END IF

  END SUBROUTINE getParameters_Orb





  !! *Description*:
  !!
  !! Returns the periapsis distance q for the input orbit. The
  !! semimajor axis 'a' is computed if it is not given as a
  !! parameter. Optionally, returns partial derivatives between q and
  !! Keplerian elements.
  !!
  !! Returns error.
  !!
  SUBROUTINE getPeriapsisDistance_Orb(this, q, a, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                      :: this
    REAL(bp), INTENT(out)                         :: q
    REAL(bp), INTENT(in), OPTIONAL                :: a
    REAL(bp), DIMENSION(6), INTENT(out), OPTIONAL :: partials ! Partials

    TYPE (Orbit) :: this_
    REAL(bp), DIMENSION(3) :: pos, vel
    REAL(bp) :: a_, e, r, e_sin_ea, e_cos_ea, alpha, gamma

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPeriapsisDistance", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    SELECT CASE (this%element_type)

    CASE ("cometary")

       q = this%elements(1)
       IF (PRESENT(partials)) THEN
          partials = 0.0_bp
          partials(1) = 1.0_bp 
       END IF

    CASE ("keplerian")

       q = this%elements(1) * (1.0_bp - this%elements(2))
       IF (PRESENT(partials)) THEN
          partials = 0.0_bp
          partials(1:2) = (/ 1.0_bp - this%elements(2), -this%elements(1) /)
       END IF

    CASE ("cartesian")

       this_ = copy(this)
       CALL rotateToEcliptic(this_)
       pos = this_%elements(1:3)
       vel = this_%elements(4:6)
       r = SQRT(DOT_PRODUCT(pos,pos))
       ! alpha = rv^2/mu
       alpha = r*DOT_PRODUCT(vel,vel)/planetary_mu(this_%center)
       IF (PRESENT(a)) THEN
          a_ = a
       ELSE
          ! h = v^2/2 - mu/r and for elliptical orbits a = -mu/2h
          ! -> a = r / (2 - rv^2/mu)
          IF (ABS(2.0_bp - alpha) < 10.0_bp*EPSILON(alpha)) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / getPeriapsisDistance", &
                  "a is approaching infinity.", 1)
             RETURN
          END IF
          a_ = r / (2.0_bp - alpha)
       END IF
       ! gamma = sqrt(mu*a)
       gamma = SQRT(planetary_mu(this_%center)*a_)
       IF (ABS(gamma) < TINY(gamma)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "Gamma is computationally too small.", 1)
          RETURN
       END IF
       ! r(1+e*cos(f))=a|1-e^2|, r*cos(f)=a*cos(E)-ae, and for elliptical 
       ! orbits a=r/(2-rv^2/mu) (see above) and |1-e^2|=1-e^2
       ! -> e*cos(E)=rv^2/mu-1=alpha-1
       e_cos_ea = alpha - 1.0_bp
       !
       e_sin_ea = DOT_PRODUCT(pos,vel)/gamma
       e = SQRT(e_cos_ea**2.0_bp + e_sin_ea**2.0_bp)
       IF (e < 10.0_bp*EPSILON(e)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "Orbit is almost circular (1).", 1)
          WRITE(stderr,*) e
          RETURN
       END IF
       CALL NULLIFY(this_)
       IF ((1.0_bp-e)*a_ < 0.0_bp) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPeriapsisDistance", &
               "Signs of semimajor axis and eccentricity are incompatible.", 1)
          RETURN
       END IF
       q = a_*(1.0_bp - e)

       IF (PRESENT(partials)) THEN
          partials = 0.0_bp
          partials(1:2) = (/ 1.0_bp - e, -a_ /)
       END IF

    END SELECT

  END SUBROUTINE getPeriapsisDistance_Orb











  !! *Description*:
  !!
  !! Returns the previous Newtonian periapsis time computed in the
  !! 2-b approximation. Optionally, returns periapsis time closest to
  !! the given date.
  !!
  !! Returns error.
  !! 
  FUNCTION getPeriapsisTime(this, t)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)          :: this
    TYPE (Time), INTENT(in), OPTIONAL :: t
    TYPE (Time)                       :: getPeriapsisTime

    TYPE (Time) :: t_
    REAL(bp), DIMENSION(6) :: elements 
    REAL(bp) :: tt_mjd, mean_motion, tt_mjd_periapsis, &
         tt_mjd_periapsis_, period, dt_min

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPeriapsisTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this,"keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / getPeriapsisTime", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    t_ = copy(this%t)
    tt_mjd = getMJD(t_, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / getPeriapsisTime", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    CALL NULLIFY(t_)
    mean_motion = SQRT(planetary_mu(this%center)/elements(1)**3.0_bp)
    tt_mjd_periapsis = tt_mjd - elements(6)/mean_motion

    IF (PRESENT(t)) THEN
       IF (.NOT.exist(t)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPeriapsisTime", &
               "Date not initialized.", 1)
          RETURN
       END IF
       t_ = copy(t)
       tt_mjd = getMJD(t_, "TT")
       IF (error) THEN
          CALL errorMessage("Orbit / getPeriapsisTime", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL NULLIFY(t_)
       period = two_pi/mean_motion
       tt_mjd_periapsis_ = tt_mjd_periapsis
       dt_min = ABS(tt_mjd - tt_mjd_periapsis_)
       IF (tt_mjd_periapsis_ < tt_mjd) THEN
          DO WHILE (tt_mjd_periapsis_ < tt_mjd)
             tt_mjd_periapsis_ = tt_mjd_periapsis_ + period
             IF (ABS(tt_mjd - tt_mjd_periapsis_) < dt_min) THEN
                tt_mjd_periapsis = tt_mjd_periapsis_
                dt_min = ABS(tt_mjd - tt_mjd_periapsis_)
             END IF
          END DO
       ELSE
          DO WHILE (tt_mjd_periapsis_ > tt_mjd)
             tt_mjd_periapsis_ = tt_mjd_periapsis_ - period
             IF (ABS(tt_mjd - tt_mjd_periapsis_) < dt_min) THEN
                tt_mjd_periapsis = tt_mjd_periapsis_
                dt_min = ABS(tt_mjd - tt_mjd_periapsis_)
             END IF
          END DO
       END IF
    END IF

    CALL NEW(getPeriapsisTime, tt_mjd_periapsis, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / getPeriapsisTime", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF

  END FUNCTION getPeriapsisTime





  !! *Description*:
  !!
  !! Returns the phase angle for an object on this orbit at a given
  !! epoch as seen from a given observer. Optionally, returns partials
  !! of phase angle wrt the cartesian position.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getPhaseAngle_Orb(this, observer, phase_angle, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                   :: this
    TYPE (CartesianCoordinates), INTENT(in)       :: observer
    REAL(bp), INTENT(out)                         :: phase_angle
    REAL(bp), DIMENSION(6), INTENT(out), OPTIONAL :: partials

    TYPE (Orbit)                :: this_lt_corr
    TYPE (CartesianCoordinates) :: ephemeris_ccoord, observer_
    TYPE (SphericalCoordinates) :: ephemeris_scoord
    REAL(bp), DIMENSION(6,6)    :: jacobian, jacobian_lt_corr, jacobian_prop
    REAL(bp), DIMENSION(1,6)    :: partials_
    REAL(bp), DIMENSION(3)      :: observer_pos, ephemeris_pos, heliocentric_pos
    REAL(bp)                    :: cos_alpha, sqrt_

    IF (.NOT. this%is_initialized) THEN
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

    ! Topocentric light-time-corrected ephemeris:
    IF (PRESENT(partials)) THEN
       CALL getEphemeris(this, observer, ephemeris_scoord, &
            this_lt_corr=this_lt_corr, &
            jacobian_lt_corr=jacobian_lt_corr, &
            jacobian_prop=jacobian_prop)
       jacobian = MATMUL(jacobian_lt_corr, jacobian_prop)
    ELSE
       CALL getEphemeris(this, observer, ephemeris_scoord, &
            this_lt_corr=this_lt_corr)
    END IF
    IF (error) THEN
       CALL errorMessage("Orbit / getPhaseAngle", &
            "TRACE BACK 5", 1)
       RETURN
    END IF
    CALL NEW(ephemeris_ccoord, ephemeris_scoord)
    IF (error) THEN
       CALL errorMessage("Orbit / getPhaseAngle", &
            "TRACE BACK 10", 1)
       RETURN
    END IF

    IF (getFrame(this) == "equatorial") THEN
       CALL rotateToEquatorial(ephemeris_ccoord)
       CALL rotateToEquatorial(this_lt_corr)
    ELSE IF (getFrame(this) == "ecliptic") THEN
       CALL rotateToEcliptic(ephemeris_ccoord)
       CALL rotateToEcliptic(this_lt_corr)
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngle", &
            "Unknown frame.", 1)
       RETURN       
    END IF
    ephemeris_pos = getPosition(ephemeris_ccoord)
    heliocentric_pos = getPosition(this_lt_corr)

    ! Cosine of phase angle:
    cos_alpha = DOT_PRODUCT(heliocentric_pos,ephemeris_pos)/ &
         (SQRT(DOT_PRODUCT(heliocentric_pos,heliocentric_pos))* &
         SQRT(DOT_PRODUCT(ephemeris_pos,ephemeris_pos)))
    IF (ABS(cos_alpha) > 1.0_bp) THEN
       cos_alpha = SIGN(1.0_bp,cos_alpha)
    END IF

    ! Phase angle:
    IF (cos_alpha == -1.0_bp) THEN
       phase_angle = pi
    ELSE IF (cos_alpha == 1.0_bp) THEN
       phase_angle = 0.0_bp
    ELSE
       phase_angle = ACOS(cos_alpha)
    END IF

    ! Partials of phase angle wrt Cartesian position:
    IF (PRESENT(partials)) THEN
       observer_ = copy(observer)
       IF (getFrame(this) == "equatorial") THEN
          CALL rotateToEquatorial(observer_)
       ELSE IF (getFrame(this) == "ecliptic") THEN
          CALL rotateToEcliptic(observer_)
       END IF
       observer_pos = getPosition(observer_)
       CALL NULLIFY(observer_)
       partials_ = 0.0_bp
       sqrt_ = DOT_PRODUCT(heliocentric_pos,heliocentric_pos) * &
            DOT_PRODUCT(ephemeris_pos,ephemeris_pos) * (1.0_bp - &
            DOT_PRODUCT(heliocentric_pos,ephemeris_pos)**2.0_bp / &
            (DOT_PRODUCT(heliocentric_pos,heliocentric_pos) * &
            DOT_PRODUCT(ephemeris_pos,ephemeris_pos)))
       IF (sqrt_ < 10.0_bp*EPSILON(sqrt_)) THEN
          partials_(1,1:3) = 1.0E-07_bp
       ELSE
          sqrt_ = SQRT(sqrt_)
       END IF
       partials_(1,1:3) = (observer_pos - heliocentric_pos * (2.0_bp - &
            DOT_PRODUCT(heliocentric_pos,ephemeris_pos) * &
            (1.0_bp / DOT_PRODUCT(ephemeris_pos,ephemeris_pos) + &
            1.0_bp / DOT_PRODUCT(heliocentric_pos,heliocentric_pos)))) / &
            sqrt_
       partials_ = MATMUL(partials_,jacobian)
       partials(1:3) = partials_(1,1:3)
    END IF

  END SUBROUTINE getPhaseAngle_Orb





  !! *Description*:
  !!
  !! Returns the phase angles for an object on this orbit at a given
  !! epoch as seen from given observers.
  !!
  !! Returns error.
  !! 
  SUBROUTINE getPhaseAngles_Orb(this, observers, phase_angles)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                           :: this
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: observers
    REAL(bp), DIMENSION(:), POINTER                       :: phase_angles

    TYPE (Orbit), DIMENSION(:), POINTER :: this_lt_corr_arr => NULL()
    TYPE (CartesianCoordinates) :: ephemeris_ccoord, observer_
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
         ephemerides => NULL()
!!$    REAL(bp), DIMENSION(:,:,:), POINTER :: jacobian_lt_corr_arr, jacobian_prop_arr
!!$    REAL(bp), DIMENSION(6,6) :: jacobian
!!$    REAL(bp), DIMENSION(1,6) :: partials_
    REAL(bp), DIMENSION(3) :: observer_pos, ephemeris_pos, heliocentric_pos
    REAL(bp) :: cos_alpha, sqrt_
    INTEGER :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPhaseAngles", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    DO i=1,SIZE(observers)
       IF (.NOT.exist(observers(i))) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPhaseAngles", &
               "Observer object has not been initialized.", 1)
          RETURN
       END IF
    END DO

    ALLOCATE(phase_angles(SIZE(observers)))
    ! Topocentric light-time-corrected ephemeris:
!!$    IF (PRESENT(partials)) THEN
!!$       ALLOCATE(partials(SIZE(observers),6))
!!$       CALL getEphemerides(this, observers, ephemerides, &
!!$            this_lt_corr_arr=this_lt_corr_arr, &
!!$            jacobian_lt_corr_arr=jacobian_lt_corr_arr, &
!!$            jacobian_prop_arr=jacobian_prop_arr)
!!$    ELSE
    CALL getEphemerides(this, observers, ephemerides, &
         this_lt_corr_arr=this_lt_corr_arr)
!!$    END IF
    IF (error) THEN
       CALL errorMessage("Orbit / getPhaseAngles", &
            "TRACE BACK 5", 1)
       RETURN
    END IF

    DO i=1,SIZE(observers)

       CALL NEW(ephemeris_ccoord, ephemerides(i))
       IF (error) THEN
          CALL errorMessage("Orbit / getPhaseAngles", &
               "TRACE BACK 10", 1)
          RETURN
       END IF

       IF (getFrame(this) == "equatorial") THEN
          CALL rotateToEquatorial(ephemeris_ccoord)
          CALL rotateToEquatorial(this_lt_corr_arr(i))
       ELSE IF (getFrame(this) == "ecliptic") THEN
          CALL rotateToEcliptic(ephemeris_ccoord)
          CALL rotateToEcliptic(this_lt_corr_arr(i))
       ELSE
          error = .TRUE.
          CALL errorMessage("Orbit / getPhaseAngles", &
               "Unknown frame.", 1)
          RETURN       
       END IF
       ephemeris_pos = getPosition(ephemeris_ccoord)
       heliocentric_pos = getPosition(this_lt_corr_arr(i))
       CALL NULLIFY(ephemeris_ccoord)

       ! Cosine of phase angle:
       cos_alpha = DOT_PRODUCT(heliocentric_pos,ephemeris_pos)/ &
            (SQRT(DOT_PRODUCT(heliocentric_pos,heliocentric_pos))* &
            SQRT(DOT_PRODUCT(ephemeris_pos,ephemeris_pos)))
       IF (ABS(cos_alpha) > 1.0_bp) THEN
          cos_alpha = SIGN(1.0_bp,cos_alpha)
       END IF

       ! Phase angle:
       IF (cos_alpha == -1.0_bp) THEN
          phase_angles(i) = pi
       ELSE IF (cos_alpha == 1.0_bp) THEN
          phase_angles(i) = 0.0_bp
       ELSE
          phase_angles(i) = ACOS(cos_alpha)
       END IF

!!$       ! Partials of phase angle wrt Cartesian position:
!!$       IF (PRESENT(partials)) THEN
!!$          observer_ = copy(observers(i))
!!$          IF (getFrame(this) == "equatorial") THEN
!!$             CALL rotateToEquatorial(observer_)
!!$          ELSE IF (getFrame(this) == "ecliptic") THEN
!!$             CALL rotateToEcliptic(observer_)
!!$          END IF
!!$          observer_pos = getPosition(observer_)
!!$          CALL NULLIFY(observer_)
!!$          partials_ = 0.0_bp
!!$          sqrt_ = DOT_PRODUCT(heliocentric_pos,heliocentric_pos) * &
!!$               DOT_PRODUCT(ephemeris_pos,ephemeris_pos) * (1.0_bp - &
!!$               DOT_PRODUCT(heliocentric_pos,ephemeris_pos)**2.0_bp / &
!!$               (DOT_PRODUCT(heliocentric_pos,heliocentric_pos) * &
!!$               DOT_PRODUCT(ephemeris_pos,ephemeris_pos)))
!!$          IF (sqrt_ < 10.0_bp*EPSILON(sqrt_)) THEN
!!$             partials_(1,1:3) = 1.0E-07_bp
!!$          ELSE
!!$             sqrt_ = SQRT(sqrt_)
!!$          END IF
!!$          partials_(1,1:3) = (observer_pos - heliocentric_pos * (2.0_bp - &
!!$               DOT_PRODUCT(heliocentric_pos,ephemeris_pos) * &
!!$               (1.0_bp / DOT_PRODUCT(ephemeris_pos,ephemeris_pos) + &
!!$               1.0_bp / DOT_PRODUCT(heliocentric_pos,heliocentric_pos)))) / &
!!$               sqrt_
!!$          jacobian = MATMUL(jacobian_lt_corr_arr(i,:,:), jacobian_prop_arr(i,:,:))
!!$          partials_ = MATMUL(partials_,jacobian)
!!$          partials(i,1:3) = partials_(1,1:3)
!!$          partials(i,4:6) = 0.0_bp
!!$          CALL NULLIFY(observer_)
!!$       END IF

    END DO

!!$    IF (PRESENT(partials)) THEN
!!$       DEALLOCATE(this_lt_corr_arr)
!!$       DEALLOCATE(jacobian_lt_corr_arr)
!!$       DEALLOCATE(jacobian_prop_arr)
!!$    END IF

  END SUBROUTINE getPhaseAngles_Orb





  !! *Description*:
  !!
  !! Returns the previous Newtonian plane-crossing (corresponding to
  !! the date when the object is at the ascending node) time computed
  !! in the 2-b approximation. Optionally, returns plane-crossing time
  !! closest to the given date.
  !!
  !! Returns error.
  !! 
  FUNCTION getPlaneCrossingTime(this, t)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)          :: this
    TYPE (Time), INTENT(in), OPTIONAL :: t
    TYPE (Time)                       :: getPlaneCrossingTime

    TYPE (Time) :: t_
    REAL(bp), DIMENSION(6) :: elements 
    REAL(bp) :: tt_mjd, mean_motion, tt_mjd_planecrossing, &
         tt_mjd_planecrossing_, period, dt_min

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPlaneCrossingTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this,"keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / getPlaneCrossingTime", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF

    t_ = copy(this%t)
    tt_mjd = getMJD(t_, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / getPlaneCrossingTime", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    CALL NULLIFY(t_)
    mean_motion = SQRT(planetary_mu(this%center)/elements(1)**3.0_bp)
    tt_mjd_planecrossing = tt_mjd - MODULO(SUM(elements(5:6)), two_pi)/mean_motion

    IF (PRESENT(t)) THEN
       IF (.NOT.exist(t)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getPlaneCrossingTime", &
               "Date not initialized.", 1)
          RETURN
       END IF
       t_ = copy(t)
       tt_mjd = getMJD(t_, "TT")
       IF (error) THEN
          CALL errorMessage("Orbit / getPlaneCrossingTime", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL NULLIFY(t_)
       period = two_pi/mean_motion
       tt_mjd_planecrossing_ = tt_mjd_planecrossing
       dt_min = ABS(tt_mjd - tt_mjd_planecrossing_)
       IF (tt_mjd_planecrossing_ < tt_mjd) THEN
          DO WHILE (tt_mjd_planecrossing_ < tt_mjd)
             tt_mjd_planecrossing_ = tt_mjd_planecrossing_ + period
             IF (ABS(tt_mjd - tt_mjd_planecrossing_) < dt_min) THEN
                tt_mjd_planecrossing = tt_mjd_planecrossing_
                dt_min = ABS(tt_mjd - tt_mjd_planecrossing_)
             END IF
          END DO
       ELSE
          DO WHILE (tt_mjd_planecrossing_ > tt_mjd)
             tt_mjd_planecrossing_ = tt_mjd_planecrossing_ - period
             IF (ABS(tt_mjd - tt_mjd_planecrossing_) < dt_min) THEN
                tt_mjd_planecrossing = tt_mjd_planecrossing_
                dt_min = ABS(tt_mjd - tt_mjd_planecrossing_)
             END IF
          END DO
       END IF
    END IF

    CALL NEW(getPlaneCrossingTime, tt_mjd_planecrossing, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / getPlaneCrossingTime", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF

  END FUNCTION getPlaneCrossingTime





  !! *Description*:
  !!
  !! Returns Poincaré elements calculated from Delaunay's
  !! elements. The mass of the target body is assumed to be negligible
  !! compared to the mass of the Sun.
  !  (&radic = square root)
  !!
  !!   - pv(1) = &Lambda = L
  !!   - pv(2) = &xi     = &radic(2(L - G)) * cos(g + &theta)
  !!   - pv(3) = p       = &radic(2(G - &Theta)) * cos(&theta)
  !!   - pv(4) = &lambda = l + g + &theta
  !!   - pv(5) = &eta    = -&radic(2(L - G)) * sin(g + &theta)
  !!   - pv(6) = q       = -&radic(2(G - &Theta)) * sin(&theta)
  !!
  !! Returns error.
  !! 
  FUNCTION getPoincareElements(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(6)   :: getPoincareElements
    REAL(bp), DIMENSION(6)   :: de

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPoincareElements", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    de = getDelaunayElements(this)
    IF (error) THEN
       CALL errorMessage("Orbit / getPoincareElements", &
            "TRACE BACK", 1)
       RETURN
    END IF

    getPoincareElements(1) = de(4)
    getPoincareElements(2) = SQRT(2.0_bp*(de(4) - de(5))) * COS(de(2) + de(3))
    getPoincareElements(3) = SQRT(2.0_bp*(de(5) - de(6))) * COS(de(3))
    getPoincareElements(4) = de(1) + de(2) + de(3)
    ! Make sure that lambda=[0,2*pi]:
    getPoincareElements(4) = MODULO(getPoincareElements(4),two_pi)
    getPoincareElements(5) = -1.0_bp * SQRT(2.0_bp*(de(4) - de(5))) * SIN(de(2) + de(3))
    getPoincareElements(6) = -1.0_bp * SQRT(2.0_bp*(de(5) - de(6))) * SIN(de(3))

  END FUNCTION getPoincareElements





  !! *Description*:
  !!
  !! Returns the Cartesian position vector as AUs.
  !!
  !! Returns error.
  !!
  FUNCTION getPosition_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(3)   :: getPosition_Orb
    REAL(bp), DIMENSION(6)   :: celements

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getPosition", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    celements(1:6) = getCartesianElements(this, this%frame)
    IF (error) THEN
       CALL errorMessage("Orbit / getPosition", &
            "TRACE BACK", 1)
       RETURN
    END IF
    getPosition_Orb(1:3) = celements(1:3)

  END FUNCTION getPosition_Orb





  !! *Description*:
  !!
  !! Returns equatorial spherical elements.
  !!
  !! Returns error.
  !!
  FUNCTION getSCoord_Orb(this, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)    :: this
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame
    TYPE (SphericalCoordinates) :: getSCoord_Orb
    TYPE (CartesianCoordinates) :: ccoord
    CHARACTER(len=FRAME_LEN)    :: frame_

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getSCoord", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(frame)) THEN
       frame_ = frame
       CALL locase(frame_, error)
       IF (error) THEN
          CALL errorMessage("Orbit / getSCoord", &
               "The frame string contains forbidden characters.", 1)
          RETURN
       END IF
    ELSE IF (this%element_type == "cartesian") THEN
       frame_ = this%frame
    ELSE IF (this%element_type == "keplerian" .OR. &
         this%element_type == "cometary") THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getSCoord", &
            "Frame must be given explicitly if the " // &
            "input is Keplerian or cometary elements.", 1)
       RETURN
    END IF
    ! orbit -> cartesian coordinates
    ccoord = getCCoord(this, frame_) 
    IF (error) THEN
       CALL errorMessage("Orbit / getSCoord", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    ! cartesian coordinates -> spherical coordinates
    CALL NULLIFY(getSCoord_Orb)
    getSCoord_Orb = getSCoord(ccoord)
    IF (error) THEN
       CALL errorMessage("Orbit / getSCoord", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    CALL NULLIFY(ccoord)

  END FUNCTION getSCoord_Orb





  !! *Description*:
  !!
  !! Returns the solar elongation for an object on this orbit at a given
  !! epoch as seen from a given observer.
  !!
  !! Returns error.
  !! 
  REAL(bp) FUNCTION getSolarElongation_Orb(this, observer)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                   :: this
    TYPE (CartesianCoordinates), INTENT(in)       :: observer

    TYPE (CartesianCoordinates) :: observer_
    REAL(bp), DIMENSION(3)      :: asteroid2sun, asteroid2observer
    REAL(bp), DIMENSION(6)      :: asteroidelems

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getSolarElongation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT.exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getSolarElongation", &
            "Observer object has not been initialized.", 1)
       RETURN
    END IF

    ! Ecliptic position vector from asteroid to the Sun:
    asteroidelems = -1.0_bp*getElements(this, "cartesian", "ecliptic")
    asteroid2sun = asteroidelems(1:3)

    IF (error) THEN
       CALL errorMessage("Orbit / getSolarElongation", &
            "TRACE BACK (5)", 1)
       RETURN       
    END IF

    ! Ecliptic position vector from asteroid to the observer:
    observer_ = copy(observer)
    CALL rotateToEcliptic(observer_)
    asteroid2observer = asteroid2sun + getPosition(observer_)
    IF (error) THEN
       CALL errorMessage("Orbit / getSolarElongation", &
            "TRACE BACK (10)", 1)
       RETURN       
    END IF

    ! Compute the solar elongation:
    getSolarElongation_Orb = ACOS(DOT_PRODUCT(asteroid2sun,asteroid2observer) / &
         (SQRT(DOT_PRODUCT(asteroid2sun,asteroid2sun)) * &
         SQRT(DOT_PRODUCT(asteroid2observer,asteroid2observer))))

  END FUNCTION getSolarElongation_Orb





  !! *Description*:
  !!
  !! Returns the values of the Stumpff c0, c1, c2, and c3 
  !! functions.
  !!
  !! Returns error.
  !!
  !! Reference: 
  !! Danby, J. M.: Fundamentals of Celestial Mechanics, 2nd ed., 
  !!               rev. ed. Richmond, VA: Willmann-Bell, pp. 170-178, 1988
  !!
  SUBROUTINE getStumpffFunctions(x, stumpff_c)

    USE functions
    IMPLICIT NONE
    REAL(bp), INTENT(in)                  :: x
    REAL(bp), DIMENSION(0:3), INTENT(out) :: stumpff_c
    REAL(bp), PARAMETER                   :: xm = 0.1_bp 
    REAL(bp), PARAMETER                   :: xc = -1.0e6_bp ! -1.0e5
    REAL(bp)                              :: y
    INTEGER                               :: n, i

    IF (x < xc) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getStumpffFunctions", &
            "x < xc", 1)
       WRITE(stderr,*) x
       RETURN
    END IF

    y = x
    n = 0
    DO WHILE (.NOT. ABS(y) < xm)
       n = n + 1
       y = y/4.0_bp
    END DO

    stumpff_c(3) = stumpff(y,3)
    stumpff_c(2) = stumpff(y,2)
    stumpff_c(1) = stumpff(y,1,stumpff_c(3))
    stumpff_c(0) = stumpff(y,0,stumpff_c(2))

    DO i=n, 1, -1
       stumpff_c(3) = 0.25_bp * (stumpff_c(2) + &
            stumpff_c(0)*stumpff_c(3))
       stumpff_c(2) = 0.5_bp * stumpff_c(1)**2.0_bp
       stumpff_c(1) = stumpff_c(0) * stumpff_c(1)
       stumpff_c(0) = 2.0_bp * stumpff_c(0)**2.0_bp - 1.0_bp
    END DO

  END SUBROUTINE getStumpffFunctions





  !! *Description*:
  !!
  !! Returns the Time-object corresponding to this orbit.
  !!
  !! Returns error.
  !!
  FUNCTION getTime_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    TYPE (Time)                       :: getTime_Orb

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getTime", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getTime_Orb = copy(this%t)

  END FUNCTION getTime_Orb





  !! *Description*:
  !!
  !! Returns the Tisserand's parameters for this orbit wrt the 8
  !! planets, Pluto, and the Moon.
  !!
  !! Returns error.
  !!
  FUNCTION getTisserandsParameters_Orb(this) RESULT(tisserands_parameter)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(10) :: tisserands_parameter 

    TYPE (Orbit) :: orb
    TYPE (Time) :: t
    REAL(bp), DIMENSION(:,:), POINTER :: planeph => NULL()
    REAL(bp), DIMENSION(6) :: elements, elements_
    REAL(bp) :: mjd_tt
    INTEGER :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getTisserandsParameters", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this,"keplerian")
    t = getTime(this)
    mjd_tt = getMJD(t, "TT")
    planeph => JPL_ephemeris(mjd_tt, -10, 11, error)

    DO i=1,SIZE(planeph,dim=1)
       CALL NEW(orb, planeph(i,1:6), "cartesian", "equatorial", t) 
       elements_ = getElements(orb, "keplerian")
       CALL NULLIFY(orb)
       tisserands_parameter(i) = elements_(1)/elements(1) + &
            2.0_bp*SQRT(elements(1)/elements_(1)*(1.0_bp-elements(2)**2)) * &
            COS(elements(3)-elements_(3))
    END DO
    DEALLOCATE(planeph)
    CALL NULLIFY(t)

  END FUNCTION getTisserandsParameters_Orb





  !! *Description*:
  !!
  !! Returns the Jacobi constants, C, for this orbit wrt the 8
  !! planets, Pluto, and the Moon.
  !!
  !! Returns error.
  !!
  FUNCTION getJacobiConstants_Orb(this) RESULT(jacobi_constant)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(10) :: jacobi_constant 

    TYPE (Orbit) :: orb
    TYPE (Time) :: t
    REAL(bp), DIMENSION(:,:), POINTER :: planeph => NULL()
    REAL(bp), DIMENSION(6) :: elements, elements_, coordinates, &
         coordinates_
    REAL(bp) :: mjd_tt, d, mass_ratio
    INTEGER :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getJacobiConstants", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this,"keplerian")
    coordinates = getElements(this,"cartesian","ecliptic")
    t = getTime(this)
    mjd_tt = getMJD(t, "TT")
    planeph => JPL_ephemeris(mjd_tt, -10, 11, error)

    DO i=1,SIZE(planeph,dim=1)
       CALL NEW(orb, planeph(i,1:6), "cartesian", "equatorial", t) 
       IF (error) THEN
          RETURN
       END IF
       !elements_ = getElements(orb, "keplerian")
       coordinates_ = getElements(orb,"cartesian","ecliptic")
       d = SQRT(SUM(coordinates_(1:3)**2))
       mass_ratio = planetary_masses(i) / (1.0_bp+planetary_masses(i))
       CALL NULLIFY(orb)
       jacobi_constant(i) = -1.0_bp/(2.0_bp*elements(1)) - &
            SQRT(elements(1)*(1.0_bp-elements(2)**2)) * &
            COS(elements(3)-elements_(3)) - mass_ratio * &
            (d/SQRT(SUM((coordinates_(1:3)-coordinates(1:3))**2)) - &
            d/SQRT(SUM(coordinates(1:3)**2)))
    END DO
    DEALLOCATE(planeph)
    CALL NULLIFY(t)

  END FUNCTION getJacobiConstants_Orb





  !! *Description*:
  !!
  !! Computes spherical topocentric equatorial
  !! coordinates from the object's heliocentric Cartesian orbit and 
  !! the observers" heliocentric Cartesian coordinates. Both 
  !! Cartesian elements and coordinates can be given either as 
  !! equatorial or ecliptic. Optionally, partial derivatives of the
  !! spherical topocentric equatorial coordinates
  !! wrt object's heliocentric Cartesian orbital elements are
  !! calculated.
  !! 
  !! Returns error.
  !!
  SUBROUTINE getTopocentricSCoords(this, topocenters, scoords_tc, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                           :: this
    TYPE (CartesianCoordinates), DIMENSION(:), INTENT(in) :: topocenters
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER    :: scoords_tc
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL         :: partials
    INTEGER                                               :: i, err

    ALLOCATE(scoords_tc(SIZE(topocenters)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getTopocentricSCoords", &
            "Could not allocate pointer (1).", 1)
       DEALLOCATE(scoords_tc, stat=err)
       RETURN
    END IF
    IF (PRESENT(partials)) THEN
       ALLOCATE(partials(6,6,SIZE(topocenters)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / getTopocentricSCoords", &
               "Could not allocate pointer (2).", 1)
          DEALLOCATE(scoords_tc, stat=err)
          DEALLOCATE(partials, stat=err)
          RETURN
       END IF
    END IF

    CALL toCartesian(this, frame="equatorial")
    DO i=1, SIZE(topocenters,dim=1)
       !       if (present(partials)) then
       !          CALL getTopocentricSCoord(this, topocenters(i), scoords_tc(i), partials(:,:,i))
       !       else
       CALL getTopocentricSCoord(this, topocenters(i), scoords_tc(i))
       !       end if
       IF (error) THEN
          CALL errorMessage("Orbit / getTopocentricSCoords", &
               "TRACE BACK 10", 1)
          DEALLOCATE(scoords_tc, stat=err)
          IF (PRESENT(partials)) THEN
             DEALLOCATE(partials, stat=err)
             RETURN
          END IF
       END IF
    END DO

  END SUBROUTINE getTopocentricSCoords





  !! *Description*:
  !!
  !! Computes spherical topocentric equatorial
  !! coordinates from the object's heliocentric Cartesian orbit and 
  !! the observer's heliocentric Cartesian coordinates. Both 
  !! Cartesian elements and coordinates can be given either as 
  !! equatorial or ecliptic. Optionally, partial derivatives of the
  !! spherical topocentric equatorial coordinates
  !! wrt object's heliocentric Cartesian orbital elements are
  !! calculated. NOT APPLICABLE FOR THE TIME BEING: The given orbit is
  !! propagated to the coordinate epoch.
  !!
  !! Returns error.
  !!
  SUBROUTINE getTopocentricSCoord(this, topocenter, scoord_tc, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                     :: this
    TYPE (CartesianCoordinates), INTENT(in)         :: topocenter
    TYPE (SphericalCoordinates), INTENT(out)        :: scoord_tc
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: partials
    TYPE (Orbit)                                    :: this_
    TYPE (CartesianCoordinates)                     :: topocenter_
    TYPE (Time)                                     :: t_topocenter
    REAL(bp), DIMENSION(6,6)                        :: jacobian
    REAL(bp), DIMENSION(6,6)                        :: scoord_partials
    REAL(bp), DIMENSION(6)                          :: topocenter_coordinates
    INTEGER                                         :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getTopocentricSCoord", &
            "This object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(topocenter)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getTopocentricSCoord", &
            "'topocenter' has not been initialized.", 1)
       RETURN
    END IF

    this_ = copy(this)
    topocenter_ = copy(topocenter)
    t_topocenter = getTime(topocenter)
    IF (PRESENT(partials)) THEN 
       CALL propagate(this_, t_topocenter, jacobian=jacobian)
    ELSE
       CALL propagate(this_, t_topocenter)
    END IF
    IF (error) THEN
       CALL errorMessage("Orbit / getTopocentricSCoord", &
            "TRACE BACK 1", 1)
       RETURN
    END IF

    ! Transform to topocentric coordinates:
    !    this_ = copy(this)
    CALL toCartesian(this_, frame="equatorial")
    CALL rotateToEquatorial(topocenter_)
    topocenter_coordinates(1:6) = getCoordinates(topocenter_)
    this_%elements(1:6) = this_%elements(1:6) - &
         topocenter_coordinates(1:6)
    scoord_tc = getSCoord(this_)
    IF (error) THEN
       CALL errorMessage("Orbit / getTopocentricSCoord", &
            "TRACE BACK 6", 1)
       RETURN
    END IF

    ! Partial derivatives:
    IF (PRESENT(partials)) THEN

       SELECT CASE (this%element_type)

       CASE ("cartesian")

          ! Rotate each row of Jacobian (gradient) to ecliptical frame.
          DO i=1,6
             CALL rotateToEcliptic(jacobian(i,:))
          END DO
          ! Topocentric equatorial (ra,dec) wrt topocentric equatorial Cartesian. 
          CALL partialsSCoordWrtCartesian_d(this_, scoord_partials)
          IF (error) THEN
             CALL errorMessage("Orbit / getTopocentricSCoord", &
                  "TRACE BACK 8", 1)
             RETURN
          END IF

       CASE ("keplerian")

          ! Topocentric equatorial (ra,dec) wrt Keplerian orbital elements.
          CALL partialsSCoordWrtKeplerian(this_, topocenter_, scoord_partials)
          IF (error) THEN
             CALL errorMessage("Orbit / getTopocentricSCoord", &
                  "TRACE BACK 10", 1)
             RETURN
          END IF

       CASE default

          error = .TRUE.
          CALL errorMessage("Orbit / getTopocentricSCoord", &
               "Could not choose between element types" // &
               "('cartesian' or 'keplerian'): " // &
               TRIM(this%element_type), 1)
          RETURN

       END SELECT

       ! Chain rule.
       partials = MATMUL(scoord_partials, jacobian)

    END IF

  END SUBROUTINE getTopocentricSCoord





  !! *Description*:
  !!
  !! Returns the transformation matrix from polar coordinate
  !! system to the ecliptical one.
  !!
  !!
  !! Returns error.
  !! 
  FUNCTION getTransformationMatrix(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(3,3) :: getTransformationMatrix
    REAL(bp), DIMENSION(6)   :: elements
    REAL(bp), DIMENSION(3)   :: sin_angles, cos_angles

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getTransformationMatrix", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%element_type == "keplerian" .OR. &
         this%element_type == "cometary") THEN
       sin_angles = (/ SIN(this%elements(3)), SIN(this%elements(4)), SIN(this%elements(5)) /)
       cos_angles = (/ COS(this%elements(3)), COS(this%elements(4)), COS(this%elements(5)) /)
    ELSE
       elements = getElements(this, "keplerian")
       sin_angles = (/ SIN(elements(3)), SIN(elements(4)), SIN(elements(5)) /)
       cos_angles = (/ COS(elements(3)), COS(elements(4)), COS(elements(5)) /)       
    END IF

    getTransformationMatrix(1,1) = cos_angles(2)*cos_angles(3) - &
         sin_angles(2)*sin_angles(3)*cos_angles(1)
    getTransformationMatrix(1,2) = -(cos_angles(2)*sin_angles(3) + &
         sin_angles(2)*cos_angles(3)*cos_angles(1))
    getTransformationMatrix(1,3) = sin_angles(2)*sin_angles(1)

    getTransformationMatrix(2,1) = sin_angles(2)*cos_angles(3) + &
         cos_angles(2)*sin_angles(3)*cos_angles(1)
    getTransformationMatrix(2,2) = -(sin_angles(2)*sin_angles(3) - &
         cos_angles(2)*cos_angles(3)*cos_angles(1))
    getTransformationMatrix(2,3) = -cos_angles(2)*sin_angles(1)

    getTransformationMatrix(3,1) = sin_angles(3)*sin_angles(1)
    getTransformationMatrix(3,2) = cos_angles(3)*sin_angles(1)
    getTransformationMatrix(3,3) = cos_angles(1)

  END FUNCTION getTransformationMatrix





  !! *Description*:
  !!
  !! Returns the Cartesian velocity vector as AUs per day.
  !!
  !! Returns error.
  !!
  FUNCTION getVelocity_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    REAL(bp), DIMENSION(3)   :: getVelocity_Orb
    REAL(bp), DIMENSION(6)   :: celements

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / getVelocity", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    celements(1:6) = getCartesianElements(this, this%frame)
    IF (error) THEN
       CALL errorMessage("Orbit / getVelocity", &
            "TRACE BACK", 1)
       RETURN
    END IF
    getVelocity_Orb(1:3) = celements(4:6)

  END FUNCTION getVelocity_Orb





  !! *Description*:
  !!
  !! Computes the opposite Cartesian orbit, e.g., position and
  !! velocity multiplied by -1.
  !!
  !! Returns error. 
  !!
  FUNCTION opposite_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: this
    TYPE (Orbit)             :: opposite_Orb

    opposite_Orb = copy(this)
    IF (opposite_Orb%element_type /= "cartesian") THEN
       CALL toCartesian(opposite_Orb, frame=this%frame)
    END IF
    opposite_Orb%elements = -1.0_bp * this%elements

  END FUNCTION opposite_Orb





  !! *Description*:
  !!
  !! Returns the partial derivatives of Cartesian orbital elements
  !! wrt cometary orbital elements.
  !!
  !!       dx/dq      dx/de      dx/di      dx/dOmega      dx/domega      dx/dt
  !!
  !!       dy/dq      dy/de      dy/di      dy/dOmega      dy/domega      dy/dt
  !!
  !!       dz/dq      dz/de      dz/di      dz/dOmega      dz/domega      dz/dt
  !!
  !!      ddx/dq     ddx/de     ddx/di     ddx/dOmega     ddx/domega     ddx/dt
  !!
  !!      ddy/dq     ddy/de     ddy/di     ddy/dOmega     ddy/domega     ddy/dt
  !!
  !!      ddz/dq     ddz/de     ddz/di     ddz/dOmega     ddz/domega     ddz/dt
  !!
  SUBROUTINE partialsCartesianWrtCometary(this, partials, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)               :: this
    REAL(bp), DIMENSION(6,6), INTENT(out)  :: partials
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame

    CHARACTER(len=1024) :: errstr

    ! Cometary wrt Cartesian:
    IF (PRESENT(frame)) THEN
       CALL partialsCometaryWrtCartesian(this, partials, frame)
    ELSE
       CALL partialsCometaryWrtCartesian(this, partials)
    END IF
    IF (error) THEN
       CALL errorMessage("Orbit / partialsCartesianWrtCometary", &
            "TRACE BACK 5", 1)
       RETURN
    END IF
    errstr = ""
    partials = matinv(partials, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCartesianWrtCometary", &
            "From matinv in linal: " // TRIM(errstr), 1)
       RETURN
    END IF

  END SUBROUTINE partialsCartesianWrtCometary





  !! *Description*:
  !!
  !! Returns the partial derivatives of Cartesian orbital elements
  !! wrt Keplerian orbital elements.
  !!
  !!       dx/da      dx/de      dx/di      dx/dOmega      dx/domega      dx/dM
  !!
  !!       dy/da      dy/de      dy/di      dy/dOmega      dy/domega      dy/dM
  !!
  !!       dz/da      dz/de      dz/di      dz/dOmega      dz/domega      dz/dM
  !!
  !!      ddx/da     ddx/de     ddx/di     ddx/dOmega     ddx/domega     ddx/dM
  !!
  !!      ddy/da     ddy/de     ddy/di     ddy/dOmega     ddy/domega     ddy/dM
  !!
  !!      ddz/da     ddz/de     ddz/di     ddz/dOmega     ddz/domega     ddz/dM
  !!
  SUBROUTINE partialsCartesianWrtKeplerian(this, partials, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)               :: this
    REAL(bp), DIMENSION(6,6), INTENT(out)  :: partials
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame

    REAL(bp), DIMENSION(6,6) :: pecl
    REAL(bp), DIMENSION(4,6) :: ppolar
    REAL(bp), DIMENSION(3,3) :: RM
    REAL(bp), DIMENSION(6) :: kep_elements
    REAL(bp), DIMENSION(4) :: polar_coord
    REAL(bp), DIMENSION(3) :: vector3, sin_angles, cos_angles
    REAL(bp) :: tmp1, tmp2, tmp3, tmp4, b, mm
    INTEGER :: i
    CHARACTER(len=FRAME_LEN) :: frame_

    partials = 0.0_bp

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    kep_elements = getElements(this, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
            "TRACE BACK 1", 1)
       RETURN
    END IF

    ! Semiminor axis:
    b = kep_elements(1) * SQRT(1.0_bp - kep_elements(2)**2)

    ! Mean motion:
    mm = SQRT(planetary_mu(this%center))/SQRT(ABS(kep_elements(1))**3.0_bp)

    ! Sines and cosines of the inclination, the longitude of the
    ! ascending node, and the argument of periapsis:
    sin_angles(1) = SIN(kep_elements(3))
    sin_angles(2) = SIN(kep_elements(4))
    sin_angles(3) = SIN(kep_elements(5))
    cos_angles(1) = COS(kep_elements(3))
    cos_angles(2) = COS(kep_elements(4))
    cos_angles(3) = COS(kep_elements(5))

    CALL partialsPolarCoordWrtKeplerian(this, ppolar, polar_coord)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
            "TRACE BACK 10", 1)
       RETURN
    END IF

    ! Partials for Cartesian ecliptic: a, e, M
    DO i=1,3
       vector3(:) = 0.0_bp
       RM = getTransformationMatrix(this)
       IF (error) THEN
          CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
               "TRACE BACK 20", 1)
          RETURN
       END IF
       vector3(1:2) = ppolar(1:2,i)
       pecl(1:3,i) = MATMUL(RM,vector3)

       vector3(:) = 0.0_bp
       vector3(1:2) = ppolar(3:4,i)
       pecl(4:6,i) = MATMUL(RM,vector3)
    END DO
    pecl(:,6) = pecl(:,3)

    ! Partials for Cartesian ecliptic: i, Om, om
    tmp1 = -sin_angles(2)*cos_angles(3) - &
         cos_angles(2)*sin_angles(3)*cos_angles(1)

    tmp2 =  sin_angles(2)*sin_angles(3) - &
         cos_angles(2)*cos_angles(3)*cos_angles(1)

    tmp3 =  cos_angles(2)*cos_angles(3) - &
         sin_angles(2)*sin_angles(3)*cos_angles(1)

    tmp4 = -cos_angles(2)*sin_angles(3) - &
         sin_angles(2)*cos_angles(3)*cos_angles(1)

    ! Inclination:
    ! coordinates
    pecl(1,3) = polar_coord(1)*sin_angles(2)*sin_angles(3)*sin_angles(1) + &
         polar_coord(2)*sin_angles(2)*cos_angles(3)*sin_angles(1)
    pecl(2,3) =-polar_coord(1)*cos_angles(2)*sin_angles(3)*sin_angles(1) - &
         polar_coord(2)*cos_angles(2)*cos_angles(3)*sin_angles(1)
    pecl(3,3) = polar_coord(1)*sin_angles(3)*cos_angles(1) + &
         polar_coord(2)*cos_angles(3)*cos_angles(1)
    ! velocities
    pecl(4,3) = polar_coord(3)*sin_angles(2)*sin_angles(3)*sin_angles(1) + &
         polar_coord(4)*sin_angles(2)*cos_angles(3)*sin_angles(1)
    pecl(5,3) =-polar_coord(3)*cos_angles(2)*sin_angles(3)*sin_angles(1) - &
         polar_coord(4)*cos_angles(2)*cos_angles(3)*sin_angles(1)
    pecl(6,3) = polar_coord(3)*sin_angles(3)*cos_angles(1) + &
         polar_coord(4)*cos_angles(3)*cos_angles(1)

    ! Longitude of ascending node
    pecl(1,4) = polar_coord(1) * tmp1 + polar_coord(2) * tmp2
    pecl(2,4) = polar_coord(1) * tmp3 + polar_coord(2) * tmp4
    pecl(3,4) = 0.0_bp
    pecl(4,4) = polar_coord(3) * tmp1 + polar_coord(4) * tmp2
    pecl(5,4) = polar_coord(3) * tmp3 + polar_coord(4) * tmp4
    pecl(6,4) = 0.0_bp

    ! Argument of periapsis
    pecl(1,5) = polar_coord(1) * tmp4 - polar_coord(2) * tmp3
    pecl(2,5) =-polar_coord(1) * tmp2 + polar_coord(2) * tmp1
    pecl(3,5) = polar_coord(1)*cos_angles(3)*sin_angles(1) - &
         polar_coord(2)*sin_angles(3)*sin_angles(1)
    pecl(4,5) = polar_coord(3) * tmp4 - polar_coord(4) * tmp3
    pecl(5,5) =-polar_coord(3) * tmp2 + polar_coord(4) * tmp1
    pecl(6,5) = polar_coord(3)*cos_angles(3)*sin_angles(1) - &
         polar_coord(4)*sin_angles(3)*sin_angles(1)

    IF (PRESENT(frame)) THEN
       frame_ = frame
       CALL locase(frame_, error)
       IF (error) THEN
          CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
               "The frame string contains forbidden characters.", 1)
          RETURN
       END IF
    ELSE IF (this%element_type == "cartesian") THEN
       frame_ = this%frame
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
            "Coordinate frame undefined.", 1)
       RETURN
    END IF

    IF (frame_ == "equatorial") THEN
       DO i=1,6 ! i is element
          CALL rotateToEquatorial(pecl(:,i))
          partials(:,i) = pecl(:,i)
       END DO
    ELSE IF (frame_ == "ecliptic") THEN
       partials = pecl
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCartesianWrtKeplerian", &
            "Unknown frame:" // TRIM(frame_) // ".", 1)
       RETURN       
    END IF

  END SUBROUTINE partialsCartesianWrtKeplerian






  !! *Description*:
  !!
  !! Returns the partial derivatives of cometary orbital elements
  !! wrt Keplerian orbital elements.
  !!
  !!       dq/da      dq/de      dq/di      dq/dOmega      dq/domega      dq/dM
  !!
  !!       de/da      de/de      de/di      de/dOmega      de/domega      de/dM
  !!
  !!       di/da      di/de      di/di      di/dOmega      di/domega      di/dM
  !!
  !!   dOmega/da  dOmega/de  dOmega/di  dOmega/dOmega  dOmega/domega  dOmega/dM
  !!
  !!   domega/da  domega/de  domega/di  domega/dOmega  domega/domega  domega/dM
  !!
  !!       dt/da      dt/de      dt/di      dt/dOmega      dt/domega      dt/dM
  !!
  SUBROUTINE partialsCometaryWrtKeplerian(this, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)               :: this
    REAL(bp), DIMENSION(6,6), INTENT(out)  :: partials

    REAL(bp), DIMENSION(6) :: elements

    ! Keplerian elements:
    elements = getElements(this, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / partialsCometaryWrtKeplerian", &
            "TRACE BACK 5", 1)
       RETURN
    END IF
    partials = identity_matrix(6)
    ! 1 - e
    partials(1,1) = 1.0_bp - elements(2)
    ! a
    partials(1,2) = elements(1)
    ! Inverse of mean motion:
    partials(6,6) = SQRT(elements(1)**3.0_bp)/SQRT(planetary_mu(this%center))

  END SUBROUTINE partialsCometaryWrtKeplerian





  !! *Description*:
  !!
  !! Returns the partial derivatives of Keplerian orbital elements
  !! wrt Cartesian orbital elements.
  !!
  !! partials-matrix:
  !!
  !!       dq/dx     dq/dy     dq/dz     dq/ddx     dq/ddy     dq/ddz
  !!       de/dx     de/dy     de/dz     de/ddx     de/ddy     de/ddz
  !!       di/dx     di/dy     di/dz     di/ddx     di/ddy     di/ddz
  !!   dOmega/dx dOmega/dy dOmega/dz dOmega/ddx dOmega/ddy dOmega/ddz
  !!   domega/dx domega/dy domega/dz domega/ddx domega/ddy domega/ddz
  !!       dt/dx     dt/dy     dt/dz     dt/ddx     dt/ddy     dt/ddz
  !!
  SUBROUTINE partialsCometaryWrtCartesian(this, partials, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)               :: this
    REAL(bp), DIMENSION(6,6), INTENT(out)  :: partials
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame

    CHARACTER(len=FRAME_LEN) :: frame_
    REAL(bp), DIMENSION(6,6) :: com_kep, kep_car

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCometaryWrtCartesian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(frame)) THEN
       frame_ = frame
       CALL locase(frame_, error)
       IF (error) THEN
          CALL errorMessage("Orbit / partialsCometaryWrtCartesian", &
               "The frame string contains forbidden characters.", 1)
          RETURN
       END IF
    ELSE IF (this%element_type == "cartesian") THEN
       frame_ = this%frame
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCometaryWrtCartesian", &
            "Coordinate frame undefined.", 1)
       RETURN
    END IF
    IF (frame_ /= "ecliptic" .AND. frame_ /= "equatorial") THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsCometaryWrtCartesian", &
            "Unknown frame:" // TRIM(frame_) // ".", 1)
       RETURN       
    END IF

    CALL partialsKeplerianWrtCartesian(this, kep_car, frame_)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsCometaryWrtCartesian", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    CALL partialsCometaryWrtKeplerian(this, com_kep)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsCometaryWrtCartesian", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    partials = MATMUL(com_kep,kep_car)

  END SUBROUTINE partialsCometaryWrtCartesian





  !! *Description*:
  !!
  !! Returns the partial derivatives of Keplerian orbital elements
  !! wrt Cartesian orbital elements.
  !!
  !! partials-matrix:
  !!
  !!       da/dx     da/dy     da/dz     da/ddx     da/ddy     da/ddz
  !!       de/dx     de/dy     de/dz     de/ddx     de/ddy     de/ddz
  !!       di/dx     di/dy     di/dz     di/ddx     di/ddy     di/ddz
  !!   dOmega/dx dOmega/dy dOmega/dz dOmega/ddx dOmega/ddy dOmega/ddz
  !!   domega/dx domega/dy domega/dz domega/ddx domega/ddy domega/ddz
  !!       dM/dx     dM/dy     dM/dz     dM/ddx     dM/ddy     dM/ddz
  !!
  SUBROUTINE partialsKeplerianWrtCartesian(this, partials, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)               :: this
    REAL(bp), DIMENSION(6,6), INTENT(out)  :: partials
    CHARACTER(len=*), INTENT(in), OPTIONAL :: frame

    CHARACTER(len=FRAME_LEN) :: frame_
    REAL(bp), DIMENSION(6) :: kep_elements, car_elements, dr, dv, &
         da, db, dc, du, dea
    REAL(bp), DIMENSION(3) :: vector3, sin_angles, cos_angles
    REAL(bp) :: tmp4, tmp5, ea, cea, sea, b, r, v, rv
    INTEGER :: i

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(frame)) THEN
       frame_ = frame
       CALL locase(frame_, error)
       IF (error) THEN
          CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
               "The frame string contains forbidden characters.", 1)
          RETURN
       END IF
    ELSE IF (this%element_type == "cartesian") THEN
       frame_ = this%frame
    ELSE
       error = .TRUE.
       CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
            "Coordinate frame undefined.", 1)
       RETURN
    END IF
    IF (frame_ /= "ecliptic" .AND. frame_ /= "equatorial") THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
            "Unknown frame:" // TRIM(frame_) // ".", 1)
       RETURN       
    END IF

    ! (Ecliptic heliocentric) Keplerian elements:
    kep_elements = getElements(this, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
            "TRACE BACK 5", 1)
       RETURN
    END IF

    ! Ecliptic heliocentric Cartesian elements:
    car_elements = getElements(this, "cartesian", frame="ecliptic")
    IF (error) THEN
       CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
            "TRACE BACK 10", 1)
       RETURN
    END IF

    ! Semiminor axis:
    b = kep_elements(1) * SQRT(1.0_bp - kep_elements(2)**2)

    ! Sines and cosines of the inclination, the longitude of the
    ! ascending node, and the argument of periapsis:
    sin_angles(1) = SIN(kep_elements(3))
    sin_angles(2) = SIN(kep_elements(4))
    sin_angles(3) = SIN(kep_elements(5))
    cos_angles(1) = COS(kep_elements(3))
    cos_angles(2) = COS(kep_elements(4))
    cos_angles(3) = COS(kep_elements(5))

    ! Sine and cosine of eccentric anomaly:
    CALL solveKeplerEquation(this, this%t, ea)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsKeplerianWrtCartesian", &
            "TRACE BACK 15", 1)
       RETURN
    END IF
    cea = COS(ea)
    sea = SIN(ea)

    ! r
    r = SQRT(DOT_PRODUCT(car_elements(1:3),car_elements(1:3)))
    ! v
    v = SQRT(DOT_PRODUCT(car_elements(4:6),car_elements(4:6)))
    ! rv
    rv = DOT_PRODUCT(car_elements(1:3),car_elements(4:6))

    dr(1:3) = car_elements(1:3)/r
    dr(4:6) = 0.0_bp
    dv(1:3) = 0.0_bp
    dv(4:6) = car_elements(4:6)/v

    ! Semimajor axis
    partials(1,1:6) = kep_elements(1)*(1.0_bp + &
         kep_elements(1)*v**2/planetary_mu(this%center))*dr/r + &
         2.0_bp*kep_elements(1)**2*v*dv/planetary_mu(this%center)

    ! Eccentricity and mean anomaly
    da      = v*(v*dr+2.0_bp*r*dv)/planetary_mu(this%center)
    db(1:3) = (car_elements(4:6)-rv*partials(1,1:3)/(2.0_bp*kep_elements(1))) / &
         (SQRT(planetary_mu(this%center))*SQRT(kep_elements(1)))
    db(4:6) = (car_elements(1:3)-rv*partials(1,4:6)/(2.0_bp*kep_elements(1))) / &
         (SQRT(planetary_mu(this%center))*SQRT(kep_elements(1)))
    dea     = (-sea*da+cea*db)/kep_elements(2)

    partials(2,1:6) = cea*da+sea*db
    partials(6,1:6) = dea*(1.0_bp-kep_elements(2)*cea) - &
         sea*partials(2,1:6)

    ! Inclination and longitude of ascending node
    vector3 = cross_product(car_elements(1:3),car_elements(4:6))
    tmp4 = SQRT(DOT_PRODUCT(vector3,vector3))

    da = (/ 0.0_bp, car_elements(6), -car_elements(5), &
         0.0_bp, -car_elements(3), car_elements(2) /)
    db = (/ -car_elements(6), 0.0_bp, car_elements(4), &
         car_elements(3), 0.0_bp, -car_elements(1) /)
    dc = (/ car_elements(5), -car_elements(4), 0.0_bp, &
         -car_elements(2), car_elements(1), 0.0_bp /)

    du = (vector3(1)*da + vector3(2)*db + vector3(3)*dc)/tmp4
    da = da/tmp4 - vector3(1)*du/tmp4**2
    db = db/tmp4 - vector3(2)*du/tmp4**2
    dc = dc/tmp4 - vector3(3)*du/tmp4**2

    tmp5 = SQRT(1.0_bp-(vector3(3)/tmp4)**2)

    partials(3,1:6) = -dc/sin_angles(1)
    partials(4,1:6) = (cos_angles(2)*da + sin_angles(2)*db + &
         (cos_angles(2)*vector3(1) + &
         sin_angles(2)*vector3(2))*vector3(3)*dc/ &
         (tmp4*tmp5)**2)/tmp5

    ! Argument of periapsis
    da = partials(1,1:6)*(cea-kep_elements(2)) + &
         kep_elements(1)*(-sea*dea-partials(2,1:6))
    db = partials(1,1:6)*SQRT(1.0_bp-kep_elements(2)**2)*sea - &
         partials(2,1:6)*kep_elements(1)*kep_elements(2)*sea / &
         SQRT(1.0_bp-kep_elements(2)**2) + dea*b*cea
    dc = 0.0_bp
    dc(3) = 1.0_bp

    da = da*sin_angles(1) + &
         kep_elements(1)*(cea-kep_elements(2))*cos_angles(1)*partials(3,1:6)
    db = db*sin_angles(1) + b*sea*cos_angles(1)*partials(3,1:6)

    partials(5,1:6) = (dc-da*sin_angles(3)-db*cos_angles(3)) / &
         (kep_elements(1)*(cea-kep_elements(2))*sin_angles(1)*cos_angles(3) - &
         b*sea*sin_angles(1)*sin_angles(3))

    IF (frame_ == "equatorial") THEN
       DO i=1,6 ! i is element
          CALL rotateToEquatorial(partials(i,:))
       END DO
    END IF

  END SUBROUTINE partialsKeplerianWrtCartesian





  !! *Description*:
  !!
  !! Returns the partial derivatives of Keplerian orbital elements
  !! wrt cometary orbital elements.
  !!
  !!       da/dq      da/de      da/di      da/dOmega      da/domega      da/dt
  !!
  !!       de/dq      de/de      de/di      de/dOmega      de/domega      de/dt
  !!
  !!       di/dq      di/de      di/di      di/dOmega      di/domega      di/dt
  !!
  !!   dOmega/dq  dOmega/de  dOmega/di  dOmega/dOmega  dOmega/domega  dOmega/dt
  !!
  !!   domega/dq  domega/de  domega/di  domega/dOmega  domega/domega  domega/dt
  !!
  !!       dM/dq      dM/de      dM/di      dM/dOmega      dM/domega      dM/dt
  !!
  SUBROUTINE partialsKeplerianWrtCometary(this, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)               :: this
    REAL(bp), DIMENSION(6,6), INTENT(out)  :: partials

    CHARACTER(len=1024) :: errstr

    ! Cometary wrt Keplerian:
    CALL partialsCometaryWrtKeplerian(this, partials)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsKeplerianWrtCometary", &
            "TRACE BACK 5", 1)
       RETURN
    END IF
    errstr = ""
    partials = matinv(partials, errstr)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsKeplerianWrtCometary", &
            "From matinv in linal: " // TRIM(errstr), 1)
       RETURN
    END IF

  END SUBROUTINE partialsKeplerianWrtCometary





  !! *Description*:
  !!
  !! Returns the partial derivatives of polar coordinates
  !! wrt Keplerian orbital elements.
  !!
  !! partials-matrix:
  !!
  !!       dr/da      dr/de      dr/di      dr/dOmega      dr/domega      dr/dM
  !!
  !!   dalpha/da  dalpha/de  dalpha/di  dalpha/dOmega  dalpha/domega  dalpha/dM
  !!
  !!   ddelta/da  ddelta/de  ddelta/di  ddelta/dOmega  ddelta/domega  ddelta/dM
  !!
  !!      ddr/da     ddr/de     ddr/di     ddr/dOmega     ddr/domega     ddr/dM
  !!
  !!  ddalpha/da ddalpha/de ddalpha/di ddalpha/dOmega ddalpha/domega ddalpha/dM
  !!
  !!  dddelta/da dddelta/de dddelta/di dddelta/dOmega dddelta/domega dddelta/dM
  !!
  SUBROUTINE partialsPolarCoordWrtKeplerian(this, partials, polar_coord)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                      :: this
    REAL(bp), DIMENSION(4,6), INTENT(out)         :: partials
    REAL(bp), DIMENSION(4), INTENT(out), OPTIONAL :: polar_coord

    REAL(bp), DIMENSION(3) :: xpolar
    REAL(bp), DIMENSION(6) :: kep_elements
    REAL(bp) :: ea, cea, sea, r, b, mm

    partials = 0.0_bp

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsPolarCoordWrtKeplerian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    kep_elements = getElements(this, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / partialsPolarCoordWrtKeplerian", &
            "TRACE BACK 5", 1)
       RETURN
    END IF

    ! Semiminor axis:
    b = kep_elements(1) * SQRT(1.0_bp - kep_elements(2)**2)
    ! Mean motion:
    mm = SQRT(planetary_mu(this%center))/SQRT(ABS(kep_elements(1))**3.0_bp)

    CALL solveKeplerEquation(this, this%t, ea)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsPolarCoordWrtKeplerian", &
            "TRACE BACK 15", 1)
       RETURN
    END IF

    cea = COS(ea)
    sea = SIN(ea)
    r = kep_elements(1)*(1.0_bp - kep_elements(2)*cea)
    xpolar(1) = kep_elements(1)*(cea-kep_elements(2))
    xpolar(2) = b*sea
    xpolar(3) = 0.0_bp

    IF (PRESENT(polar_coord)) THEN
       polar_coord(1:2)=xpolar(1:2)
       polar_coord(3)=-mm*kep_elements(1)**2*sea/r
       polar_coord(4)=mm*kep_elements(1)*b*cea/r
    END IF

    ! NOTE: THE USE OF ARRAY ELEMENTS CAN BE MISLEADING,
    ! partials wrt M should be the 6th?
    ! Partials for polar coordinates: a, e, M (others -> 0.0)
    ! Position: 
    partials(1,1) = cea-kep_elements(2)
    partials(1,2) = -(kep_elements(1)*sea)**2/r - &
         kep_elements(1)
    partials(1,3) = -sea*kep_elements(1)**2/r

    partials(2,1) = b*sea/kep_elements(1)
    partials(2,2) = sea*(cea-kep_elements(2)) * &
         kep_elements(1)**3/(b*r)
    partials(2,3) = kep_elements(1)*b*cea/r

    ! Velocity:
    !partials(3,1) = 0.5_bp*mm*kep_elements(1) * &
    !     xpolar(2)/(r*b) + 1.5_bp*(ma-this%ma)*mm*kep_elements(1)**2*xpolar(1)/r**3
    partials(3,1) = 0.5_bp*mm*kep_elements(1) * &
         xpolar(2)/(r*b)
    partials(3,2) = -mm*kep_elements(1)**2*xpolar(2) * &
         ((kep_elements(1)+r)*xpolar(1) + &
         r*kep_elements(1)*kep_elements(2))/(b*r**3)
    partials(3,3) = -mm*kep_elements(1)**3*xpolar(1)/r**3

    partials(4,1) = -0.5_bp * mm * b * &
         (xpolar(1)+kep_elements(1)*kep_elements(2)) / &
         (kep_elements(1)*r) 
    !    partials(4,1) = -0.5_bp * mm * b * &
    !         (xpolar(1)+kep_elements(1)*kep_elements(2)) / &
    !         (kep_elements(1)*r) + 1.5_bp*(ma-this%ma)*mm*kep_elements(1)**2*xpolar(2)/r**3
    partials(4,2) = mm * kep_elements(1)**2 * b * &
         (xpolar(1) * (xpolar(1) + kep_elements(1) * &
         kep_elements(2)) / kep_elements(1)- &
         r**2*(kep_elements(1)-r+xpolar(2)**2/r)/b**2)/r**3
    partials(4,3) = -mm * kep_elements(1)**3 * xpolar(2) / &
         r**3

  END SUBROUTINE partialsPolarCoordWrtKeplerian





  !! *Description*:
  !!
  !! Returns the partial derivatives of spherical coordinates wrt
  !! Cartesian orbital elements. The partial derivatives are obtained
  !! indirectly by utilizing partial derivatives of Keplerian elements
  !! in the process.
  !!
  !! partials-matrix:
  !!
  !!       dr/dx      dr/dy      dr/dz      dr/ddx      dr/ddy      dr/ddz
  !!
  !!   dalpha/dx  dalpha/dy  dalpha/dz  dalpha/ddx  dalpha/ddy  dalpha/ddz
  !!
  !!   ddelta/dx  ddelta/dy  ddelta/dz  ddelta/ddx  ddelta/ddy  ddelta/ddz
  !!
  !!      ddr/dx     ddr/dy     ddr/dz     ddr/ddx     ddr/ddy     ddr/ddz
  !!
  !!  ddalpha/dx ddalpha/dy ddalpha/dz ddalpha/ddx ddalpha/ddy ddalpha/ddz
  !!
  !!  dddelta/dx dddelta/dy dddelta/dz dddelta/ddx dddelta/ddy dddelta/ddz
  !!
  SUBROUTINE partialsSCoordWrtCartesian_i(this, observer, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                :: this
    TYPE (CartesianCoordinates), INTENT(in) :: observer
    REAL(bp), DIMENSION(6,6), INTENT(out)   :: partials

    REAL(bp), DIMENSION(6,6) :: scoord_kep, kep_car

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtCartesian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtCartesian", &
            "Observer has not yet been initialized.", 1)
       RETURN
    END IF

    CALL partialsSCoordWrtKeplerian(this, observer, scoord_kep)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsSCoordWrtCartesian", &
            "TRACE BACK (5).", 1)
       RETURN
    END IF

    ! Partials between Keplerian and cometary elements
    CALL partialsKeplerianWrtCartesian(this, kep_car, this%frame)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsSCoordWrtCartesian", &
            "TRACE BACK (10).", 1)
       RETURN
    END IF
    partials = MATMUL(scoord_kep,kep_car)

!!$    call matrix_print(scoord_kep, stdout, errstr)
!!$    write(stdout,*)
!!$    write(stdout,*) cond_nr(scoord_kep, errstr)
!!$    write(stdout,*)
!!$    call matrix_print(kep_car, stdout, errstr)
!!$    write(stdout,*)
!!$    write(stdout,*) cond_nr(kep_car, errstr)
!!$    write(stdout,*)

  END SUBROUTINE partialsSCoordWrtCartesian_i





  !! *Description*:
  !!
  !! Returns the partial derivatives of spherical coordinates wrt
  !! Cartesian orbital elements. Note that the origin and frame for
  !! both coordinate systems must coincide.
  !!
  !! partials-matrix:
  !!
  !!       dr/dx      dr/dy      dr/dz      dr/ddx      dr/ddy      dr/ddz
  !!
  !!   dalpha/dx  dalpha/dy  dalpha/dz  dalpha/ddx  dalpha/ddy  dalpha/ddz
  !!
  !!   ddelta/dx  ddelta/dy  ddelta/dz  ddelta/ddx  ddelta/ddy  ddelta/ddz
  !!
  !!      ddr/dx     ddr/dy     ddr/dz     ddr/ddx     ddr/ddy     ddr/ddz
  !!
  !!  ddalpha/dx ddalpha/dy ddalpha/dz ddalpha/ddx ddalpha/ddy ddalpha/ddz
  !!
  !!  dddelta/dx dddelta/dy dddelta/dz dddelta/ddx dddelta/ddy dddelta/ddz
  !!
  SUBROUTINE partialsSCoordWrtCartesian_d(this, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)              :: this
    REAL(bp), DIMENSION(6,6), INTENT(out) :: partials

    REAL(bp), DIMENSION(3) :: pos, vel
    REAL(bp) :: x2, y2, z2, xy, xy2, x2y2, x2y2z2, inv_x2y2, &
         inv_sqrt_x2y2, inv_x2y2z2, inv_sqrt_x2y2z2, xdxydy, &
         sum3x23y2z2

    partials = 0.0_bp

    ! Common terms.
    pos = getPosition(this)
    vel = getVelocity(this)
    x2 = pos(1)**2.0_bp
    y2 = pos(2)**2.0_bp
    z2 = pos(3)**2.0_bp
    xy = SUM(pos(1:2))
    xy2 = xy**2.0_bp
    x2y2 = x2 + y2
    inv_x2y2 = 1.0_bp/x2y2
    inv_sqrt_x2y2 = SQRT(inv_x2y2)
    x2y2z2 = x2y2 + z2    
    inv_x2y2z2 = 1.0_bp/x2y2z2
    inv_sqrt_x2y2z2 = SQRT(inv_x2y2z2)
    xdxydy = pos(1)*vel(1) + pos(2)*vel(2)
    sum3x23y2z2 = 3.0_bp*x2 + 3.0_bp*y2 + z2

    ! First row: r
    partials(1,1) = pos(1)*inv_sqrt_x2y2z2
    partials(1,2) = pos(2)*inv_sqrt_x2y2z2
    partials(1,3) = pos(3)*inv_sqrt_x2y2z2

    ! Second row: ra
    partials(2,1) = -pos(2)*inv_x2y2
    partials(2,2) = pos(1)*inv_x2y2

    ! Third row: dec
    partials(3,1) = -pos(1)*pos(3)*inv_sqrt_x2y2*inv_x2y2z2
    partials(3,2) = -pos(2)*pos(3)*inv_sqrt_x2y2*inv_x2y2z2
    partials(3,3) = inv_x2y2z2/inv_sqrt_x2y2

    ! Fourth row: dr/dt
    partials(4,1) = (vel(1)*(y2 + z2) - &
         pos(1)*(pos(2)*vel(2) + pos(3)*vel(3))) * inv_x2y2z2
    partials(4,2) = (vel(2)*(x2 + z2) - &
         pos(2)*(pos(1)*vel(1) + pos(3)*vel(3))) * inv_x2y2z2
    partials(4,3) = (vel(3)*x2y2 - pos(3)*xdxydy) * inv_x2y2z2
    partials(4,4) = pos(1)*inv_sqrt_x2y2z2
    partials(4,5) = pos(2)*inv_sqrt_x2y2z2
    partials(4,6) = pos(3)*inv_sqrt_x2y2z2

    ! Fifth row: dra/dt
    partials(5,1) = (vel(1)*pos(2)*(2.0_bp*pos(1)+pos(2)) - &
         x2*vel(2)) / (x2*xy2)
    partials(5,2) = -(vel(1) + vel(2)) / xy2
    partials(5,4) = -pos(2) / (pos(1)*xy)
    partials(5,5) = 1.0_bp / xy

    ! Sixth row: ddec/dt
    partials(6,1) = ((2.0_bp*vel(3)*pos(1) - pos(3)*vel(1))*x2y2z2 * &
         x2y2 - pos(1)*sum3x23y2z2*(vel(3)*x2y2 - pos(3)*xdxydy)) * &
         inv_x2y2z2**2.0_bp * inv_sqrt_x2y2**3.0_bp
    partials(6,2) = ((2.0_bp*vel(3)*pos(2) - pos(3)*vel(2))*x2y2z2 * &
         x2y2 - pos(2)*sum3x23y2z2*(vel(3)*x2y2 - pos(3)*xdxydy)) * &
         inv_x2y2z2**2.0_bp * inv_sqrt_x2y2**3.0_bp
    partials(6,3) = (-xdxydy*(x2y2-z2) - 2.0_bp*pos(3)*vel(3)*x2y2)* &
         inv_x2y2z2**2.0_bp*inv_sqrt_x2y2
    partials(6,4) = -pos(1)*pos(3)*inv_x2y2z2*inv_sqrt_x2y2
    partials(6,5) = -pos(2)*pos(3)*inv_x2y2z2*inv_sqrt_x2y2
    partials(6,6) = inv_x2y2z2 * SQRT(x2y2)

  END SUBROUTINE partialsSCoordWrtCartesian_d





  !! *Description*:
  !!
  !! Returns the partial derivatives of spherical coordinates wrt
  !! cometary orbital elements.
  !!
  !! partials-matrix:
  !!
  !!       dr/dq      dr/de      dr/di      dr/dOmega      dr/domega      dr/dt
  !!
  !!   dalpha/dq  dalpha/de  dalpha/di  dalpha/dOmega  dalpha/domega  dalpha/dt
  !!
  !!   ddelta/dq  ddelta/de  ddelta/di  ddelta/dOmega  ddelta/domega  ddelta/dt
  !!
  !!      ddr/dq     ddr/de     ddr/di     ddr/dOmega     ddr/domega     ddr/dt
  !!
  !!  ddalpha/dq ddalpha/de ddalpha/di ddalpha/dOmega ddalpha/domega ddalpha/dt
  !!
  !!  dddelta/dq dddelta/de dddelta/di dddelta/dOmega dddelta/domega dddelta/dt
  !!
  SUBROUTINE partialsSCoordWrtCometary(this, observer, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                :: this
    TYPE (CartesianCoordinates), INTENT(in) :: observer
    REAL(bp), DIMENSION(6,6), INTENT(out)   :: partials

    REAL(bp), DIMENSION(6,6) :: scoord_kep, kep_com

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtCometary", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtCometary", &
            "Observer has not yet been initialized.", 1)
       RETURN
    END IF

    CALL partialsSCoordWrtKeplerian(this, observer, scoord_kep)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsSCoordWrtCometary", &
            "TRACE BACK (5).", 1)
       RETURN
    END IF

    ! Partials between Keplerian and cometary elements
    CALL partialsKeplerianWrtCometary(this, kep_com)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsSCoordWrtCometary", &
            "TRACE BACK (10).", 1)
       RETURN
    END IF
    partials = MATMUL(scoord_kep,kep_com)

  END SUBROUTINE partialsSCoordWrtCometary





  !! *Description*:
  !!
  !! Returns the partial derivatives of spherical coordinates wrt
  !! Keplerian orbital elements.
  !!
  !! partials-matrix:
  !!
  !!       dr/da      dr/de      dr/di      dr/dOmega      dr/domega      dr/dM
  !!
  !!   dalpha/da  dalpha/de  dalpha/di  dalpha/dOmega  dalpha/domega  dalpha/dM
  !!
  !!   ddelta/da  ddelta/de  ddelta/di  ddelta/dOmega  ddelta/domega  ddelta/dM
  !!
  !!      ddr/da     ddr/de     ddr/di     ddr/dOmega     ddr/domega     ddr/dM
  !!
  !!  ddalpha/da ddalpha/de ddalpha/di ddalpha/dOmega ddalpha/domega ddalpha/dM
  !!
  !!  dddelta/da dddelta/de dddelta/di dddelta/dOmega dddelta/domega dddelta/dM
  !!
  SUBROUTINE partialsSCoordWrtKeplerian(this, observer, partials)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in)                :: this
    TYPE (CartesianCoordinates), INTENT(in) :: observer
    REAL(bp), DIMENSION(6,6), INTENT(out)   :: partials

    CHARACTER(len=FRAME_LEN) :: frame
    REAL(bp), DIMENSION(6,6) :: pcar_equ
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp), DIMENSION(3) :: xequ
    REAL(bp) :: tmp1, tmp2, tmp3, tmp4, tmp5
    INTEGER :: i

    partials = 0.0_bp

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtKeplerian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(observer)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtKeplerian", &
            "Observer has not yet been initialized.", 1)
       RETURN
    END IF

    frame = getFrame(observer)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / partialsSCoordWrtKeplerian", &
            "TRACE BACK (30).", 1)
       RETURN
    END IF

    elements = getElements(this, "cartesian", frame)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsSCoordWrtKeplerian", &
            "TRACE BACK (20).", 1)
       RETURN
    END IF
    xequ = elements(1:3)

    ! Partials for Cartesian equatorial coordinates
    CALL partialsCartesianWrtKeplerian(this, pcar_equ, frame=frame)
    IF (error) THEN
       CALL errorMessage("Orbit / partialsSCoordWrtKeplerian", &
            "TRACE BACK (25).", 1)
       RETURN
    END IF

    ! Partials for spherical equatorial coordinates
    ! Topocentric coordinates:
    xequ = xequ - getPosition(observer)

    tmp1 = SQRT(DOT_PRODUCT(xequ,xequ))
    tmp2 = SQRT(xequ(1)**2+xequ(2)**2)/tmp1
    tmp3 = xequ(3)/tmp1
    tmp4 = xequ(1)/(tmp1*tmp2)
    tmp5 = xequ(2)/(tmp1*tmp2)

    DO i=1,6 ! i is element
       ! rho
       partials(1,i) = tmp2*tmp4*pcar_equ(1,i) + &
            tmp2*tmp5*pcar_equ(2,i) + tmp3*pcar_equ(3,i)
       ! RA
       partials(2,i) = (tmp4*pcar_equ(2,i) - &
            tmp5*pcar_equ(1,i))/(tmp1*tmp2)
       ! Dec
       partials(3,i) = (tmp2*pcar_equ(3,i) - &
            tmp3*tmp4*pcar_equ(1,i) - tmp3*tmp5*pcar_equ(2,i))/tmp1
    END DO

    ! Note that x/dx and dx/x are missing! the following numbers are
    ! there to allow, e.g., computation of the determinant of the
    ! matrix:
    partials(4,4)   = 1.0_bp
    partials(5,5)   = 1.0_bp
    partials(6,6)   = 1.0_bp

  END SUBROUTINE partialsSCoordWrtKeplerian





  !! *Description*:
  !!
  !! Propagates the Cartesian or Keplerian orbital elements to the
  !! given epoch using either the 2-body or the n-body model. The
  !! dynamical model as well as other propagation parameters can be
  !! set with the _setParameters_-routine. Optionally,
  !! the Jacobian matrix of final coordinates wrt. initial coordinates
  !! is calculated. Doesn"t alter the original orbit, if an error occurs.
  !!
  !! *2-body*propagation*:
  !! 
  !! In Keplerian elements, trivially computes a new mean anomaly.
  !!
  !! In Cartesian elements, transforms position and velocity at one
  !! epoch into those at another epoch by using f- and g-functions.
  !!
  !! *N-body*propagation*:
  !!
  !! If necessary, transforms input orbit to equatorial Cartesian
  !! (assumed heliocentric).
  !!
  !! Transforms position and velocity at one epoch into those at
  !! another epoch, using the Bulirsch-Stoer method, full force
  !! integration including nine planets and the Moon.
  !!
  !! Transforms output to match the element type of the input orbit.
  !!
  !! *Jacobian*:
  !! 
  !! Optionally, a Jacobian matrix (partial derivatives between final
  !! orbital elements and initial orbital elements) is
  !! computed. Depending on the elements used, the layout of the
  !! Jacobian matrix is either (Keplerian elements)
  !!
  !!  da1/da0     da1/de0     da1/di0     da1/dOmega0     da1/domega0     da1/dM0
  !!
  !!  de1/da0     de1/de0     de1/di0     de1/dOmega0     de1/domega0     de1/dM0
  !!
  !!  di1/da0     di1/de0     di1/di0     di1/dOmega0     di1/domega0     di1/dM0
  !!
  !!  dOmega1/da0 dOmega1/de0 dOmega1/di0 dOmega1/dOmega0 dOmega1/domega0 dOmega1/dM0
  !!
  !!  domega1/da0 domega1/de0 domega1/di0 domega1/dOmega0 domega1/domega0 domega1/dM0
  !!
  !!  dM1/da0     dM1/de0     dM1/di0     dM1/dOmega0     dM1/domega0     dM1/dM0
  !!
  !!
  !! or (Cartesian elements)
  !!  
  !!  dx1/dx0     dx1/dy0     dx1/dz0     dx1/ddx0        dx1/ddy0        dx1/ddz0
  !!
  !!  dy1/dx0     dy1/dy0     dy1/dz0     dy1/ddx0        dy1/ddy0        dy1/ddz0
  !!
  !!  dz1/dx0     dz1/dy0     dz1/dz0     dz1/ddx0        dz1/ddy0        dz1/ddz0
  !!
  !!  ddx1/dx0    ddx1/dy0    ddx1/dz0    ddx1/ddx0       ddx1/ddy0       ddx1/ddz0
  !!
  !!  ddy1/dx0    ddy1/dy0    ddy1/dz0    ddy1/ddx0       ddy1/ddy0       ddy1/ddz0
  !!
  !!  ddz1/dx0    ddz1/dy0    ddz1/dz0    ddz1/ddx0       ddz1/ddy0       ddz1/ddz0
  !!
  !!
  !! Returns error.
  !!
  !!
  SUBROUTINE propagate_Orb_single(this, t, jacobian, encounters)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)                     :: this
    TYPE (Time), INTENT(in)                         :: t
    REAL(bp), DIMENSION(6,6), INTENT(out), OPTIONAL :: jacobian
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL   :: encounters

    TYPE (Orbit), DIMENSION(1)                      :: this_
    REAL(bp), DIMENSION(:,:,:), POINTER             :: jacobian_ => NULL()
    INTEGER :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / propagate (single)", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    this_(1) = copy(this)
    IF (PRESENT(jacobian)) THEN
       IF (PRESENT(encounters)) THEN
          CALL propagate(this_, t, jacobian=jacobian_, encounters=encounters)
       ELSE
          CALL propagate(this_, t, jacobian=jacobian_)
       END IF
       IF (.NOT.error) THEN
          jacobian = jacobian_(1,:,:)
          DEALLOCATE(jacobian_, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (single)", &
                  "Could not deallocate memory.", 1)
             CALL NULLIFY(this_(1))
             RETURN
          END IF
       END IF
    ELSE
       IF (PRESENT(encounters)) THEN
          CALL propagate(this_, t, encounters=encounters)
       ELSE
          CALL propagate(this_, t)
       END IF
    END IF
    IF (error) THEN
       CALL errorMessage("Orbit / propagate (single)", &
            "TRACE BACK (single orbit)", 1)
       IF (PRESENT(jacobian)) THEN
          jacobian = identity_matrix(6)
          DEALLOCATE(jacobian_, stat=err)
       END IF
       CALL NULLIFY(this_(1))
       RETURN
    END IF
    CALL NULLIFY(this)
    this = copy(this_(1))
    CALL NULLIFY(this_(1))

  END SUBROUTINE propagate_Orb_single





  !! *Description*:
  !!
  !! Propagates the Cartesian or Keplerian orbital elements to the
  !! given epoch using either the 2-body or the n-body model. The
  !! dynamical model as well as other propagation parameters can be
  !! set with the _setParameters_-routine. Optionally,
  !! the Jacobian matrix of final coordinates wrt. initial coordinates
  !! is calculated. Doesn"t alter the original orbit, if an error occurs.
  !!
  !! *2-body*propagation*:
  !! 
  !! In Keplerian elements, trivially computes a new mean anomaly.
  !!
  !! In Cartesian elements, transforms position and velocity at one
  !! epoch into those at another epoch by using f- and g-functions.
  !!
  !! *N-body*propagation*:
  !!
  !! If necessary, transforms input orbit to equatorial Cartesian
  !! (assumed heliocentric).
  !!
  !! Transforms position and velocity at one epoch into those at
  !! another epoch, using the Bulirsch-Stoer method, full force
  !! integration including nine planets and the Moon.
  !!
  !! Transforms output to match the element type of the input orbit.
  !!
  !! *Jacobian*:
  !! 
  !! Optionally, a Jacobian matrix (partial derivatives between final
  !! orbital elements and initial orbital elements) is computed for
  !! each orbit. Depending on the elements used, the layout of the
  !! Jacobian matrix is either (Keplerian elements)
  !!
  !!  da1/da0     da1/de0     da1/di0     da1/dOmega0     da1/domega0     da1/dM0
  !!
  !!  de1/da0     de1/de0     de1/di0     de1/dOmega0     de1/domega0     de1/dM0
  !!
  !!  di1/da0     di1/de0     di1/di0     di1/dOmega0     di1/domega0     di1/dM0
  !!
  !!  dOmega1/da0 dOmega1/de0 dOmega1/di0 dOmega1/dOmega0 dOmega1/domega0 dOmega1/dM0
  !!
  !!  domega1/da0 domega1/de0 domega1/di0 domega1/dOmega0 domega1/domega0 domega1/dM0
  !!
  !!  dM1/da0     dM1/de0     dM1/di0     dM1/dOmega0     dM1/domega0     dM1/dM0
  !!
  !!
  !! or (Cartesian elements)
  !!  
  !!  dx1/dx0     dx1/dy0     dx1/dz0     dx1/ddx0        dx1/ddy0        dx1/ddz0
  !!
  !!  dy1/dx0     dy1/dy0     dy1/dz0     dy1/ddx0        dy1/ddy0        dy1/ddz0
  !!
  !!  dz1/dx0     dz1/dy0     dz1/dz0     dz1/ddx0        dz1/ddy0        dz1/ddz0
  !!
  !!  ddx1/dx0    ddx1/dy0    ddx1/dz0    ddx1/ddx0       ddx1/ddy0       ddx1/ddz0
  !!
  !!  ddy1/dx0    ddy1/dy0    ddy1/dz0    ddy1/ddx0       ddy1/ddy0       ddy1/ddz0
  !!
  !!  ddz1/dx0    ddz1/dy0    ddz1/dz0    ddz1/ddx0       ddz1/ddy0       ddz1/ddz0
  !!
  !!
  !! Returns error.
  !!
  !!
  SUBROUTINE propagate_Orb_multiple(this_arr, t, jacobian, encounters)

    IMPLICIT NONE
    TYPE (Orbit), DIMENSION(:), INTENT(inout)     :: this_arr
    TYPE (Time), INTENT(in)                       :: t
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL :: jacobian
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL :: encounters

    TYPE (Time) :: t_
    TYPE (Orbit) :: this_
    CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: &
         element_type_arr
    CHARACTER(len=FRAME_LEN) , DIMENSION(:), ALLOCATABLE :: frame_arr
    CHARACTER(len=DYN_MODEL_LEN) :: dyn_model
    CHARACTER(len=INTEGRATOR_LEN) :: integrator
    REAL(bp), DIMENSION(:,:,:), ALLOCATABLE :: partials0, jacobian_
    REAL(bp), DIMENSION(:,:), POINTER :: elm_arr => NULL()
    REAL(bp), DIMENSION(:), ALLOCATABLE :: masses
    REAL(bp), DIMENSION(6,6) :: partials1
    REAL(bp), DIMENSION(0:3) :: stumpff_cs, ffs
    REAL(bp), DIMENSION(3) :: pos, vel
    REAL(bp) :: mjd_tt, mjd_tt0, dt, r0, u, alpha, s, f, g, df, &
         dg, mean_motion, step, mu_, p
    INTEGER :: i, j, err, nthis, center, naddit, nmassive
    LOGICAL :: multiple_t0

    nthis = SIZE(this_arr,dim=1)
    IF (nthis == 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / propagate (multiple)", &
            "Orbit array is of zero length.", 1)
       RETURN
    END IF
    DO i=1,nthis
       IF (.NOT. this_arr(i)%is_initialized) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "All objects have not yet been initialized.", 1)
          WRITE(stderr,*) i
          RETURN
       END IF
    END DO

    multiple_t0 = .FALSE.
    DO i=2,nthis
       IF (.NOT.equal(this_arr(1)%t,this_arr(i)%t)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)",&
               "All orbits do not have an equal starting epoch.", 1)          
          RETURN
          multiple_t0 = .TRUE.
       END IF
       IF (error) THEN
          CALL errorMessage("Orbit / propagate (multiple)",&
               "TRACE BACK (5)", 1)
          RETURN
       END IF
    END DO

    IF (PRESENT(jacobian)) THEN
       ALLOCATE(jacobian(nthis,6,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "Could not allocate memory (5).", 1)
          RETURN
       END IF
    END IF

    IF (PRESENT(encounters)) THEN
       ALLOCATE(encounters(nthis,11,4), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "Could not allocate memory (10).", 1)
          RETURN
       END IF
    END IF

    t_ = copy(t)
    IF (equal(this_arr(1)%t,t_)) THEN
       IF (PRESENT(jacobian)) THEN
          DO i=1,nthis
             jacobian(i,:,:) = identity_matrix(6)
          END DO
       END IF
       CALL NULLIFY(t_)
       RETURN
    END IF

    dyn_model = this_arr(1)%dyn_model_prm
    IF (error) THEN
       CALL errorMessage("Orbit / propagate (multiple)",&
            "TRACE BACK (10)", 1)
       IF (PRESENT(jacobian)) THEN
          DEALLOCATE(jacobian, stat=err)
       END IF
       RETURN
    END IF
    DO i=2,nthis
       IF (dyn_model /= this_arr(i)%dyn_model_prm) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "All orbits do not share the same propagation scheme.", 1)
          IF (PRESENT(jacobian)) THEN
             DEALLOCATE(jacobian, stat=err)
          END IF
          RETURN
       END IF
    END DO

    mjd_tt0 = getMJD(this_arr(1)%t, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / propagate (multiple)",&
            "TRACE BACK (15)", 1)
       IF (PRESENT(jacobian)) THEN
          DEALLOCATE(jacobian, stat=err)
       END IF
       RETURN
    END IF
    mjd_tt = getMJD(t_, "TT")
    IF (error) THEN
       CALL errorMessage("Orbit / propagate (multiple)",&
            "TRACE BACK (20)", 1)
       IF (PRESENT(jacobian)) THEN
          DEALLOCATE(jacobian, stat=err)
       END IF
       RETURN
    END IF

    center = this_arr(1)%center
    DO i=1,nthis
       IF (center /= this_arr(i)%center) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "All orbits do not share the same central body.", 1)
          IF (PRESENT(jacobian)) THEN
             DEALLOCATE(jacobian, stat=err)
          END IF
          RETURN
       END IF
    END DO
    ! Define mu parameter
    mu_ = planetary_mu(center)

    ! Select dynamical model; either 2-body or n-body: 
    SELECT CASE (TRIM(dyn_model))

    CASE ("2-body")

       IF (info_verb >= 5) THEN
          WRITE(stdout,"(2X,A,1X,A)") &
               "Orbit / propagate_Orb_multiple:", &
               "Preparing for 2-body propagation..."
       END IF

       dt = mjd_tt - mjd_tt0
       IF (PRESENT(jacobian) .AND. ALL(this_arr(1)%finite_diff_prm > 0.0_bp)) THEN
          ALLOCATE(elm_arr(2,6), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Could not allocate memory (15).", 1)
             DEALLOCATE(jacobian, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             RETURN
          END IF
       END IF

       DO i=1,nthis

          SELECT CASE (TRIM(this_arr(i)%element_type))

          CASE ("cartesian")

             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,1X,A,I0,A)") &
                     "Orbit / propagate_Orb_multiple:", &
                     "Carrying out 2-body propagation for orbit #", i, "..."
             END IF
             r0      = SQRT(DOT_PRODUCT(this_arr(i)%elements(1:3),this_arr(i)%elements(1:3))) ! r_0
             u       = DOT_PRODUCT(this_arr(i)%elements(1:3),this_arr(i)%elements(4:6)) 
             alpha   = 2.0_bp*mu_/r0 - DOT_PRODUCT(this_arr(i)%elements(4:6),this_arr(i)%elements(4:6)) ! 
             CALL solveKeplerEquation(r0, u, alpha, dt, stumpff_cs, ffs, s, &
                  this_arr(i)%center)
             IF (error) THEN
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "2-body propagation using Stumpff-functions was unsuccessful.", 1)
                IF (PRESENT(jacobian)) THEN
                   DEALLOCATE(jacobian, stat=err)
                END IF
                DEALLOCATE(elm_arr, stat=err)
                RETURN
             END IF
             f   = 1.0_bp - (mu_/r0) * stumpff_cs(2)
             g   = dt - mu_*stumpff_cs(3)
             df  = -mu_/(ffs(1)*r0) * stumpff_cs(1)
             dg  = 1.0_bp - (mu_/ffs(1)) * stumpff_cs(2)
             pos(1:3) = f * this_arr(i)%elements(1:3) +  g * this_arr(i)%elements(4:6)
             vel(1:3) = df * this_arr(i)%elements(1:3) + dg * this_arr(i)%elements(4:6)
             IF (PRESENT(jacobian)) THEN
                CALL GaussfgJacobian(this_arr(i), r0, u, alpha, stumpff_cs, s, &
                     f, g, df, dg, pos, ffs(1), jacobian(i,:,:))
                IF (error) THEN
                   CALL errorMessage("Orbit / propagate (multiple)",&
                        "TRACE BACK (25)", 1)
                   DEALLOCATE(jacobian, stat=err)
                   DEALLOCATE(elm_arr, stat=err)
                   RETURN
                END IF
             END IF
             this_arr(i)%elements(1:3) = pos(1:3)
             this_arr(i)%elements(4:6) = vel(1:3)
             IF (info_verb >= 5) THEN
                WRITE(stdout,"(2X,A,1X,A,I0,A)") &
                     "Orbit / propagate_Orb_multiple:", &
                     "2-body propagation carried out for orbit #", i, "."
             END IF

          CASE ("cometary")

             IF (PRESENT(jacobian)) THEN
                error = .TRUE.
                CALL errorMessage("Orbit / propagate (multiple)",&
                     "Jacobians not available for cometary elements.", 1)
                RETURN
             END IF
             p = two_pi * SQRT((this_arr(i)%elements(1) / &
                  (1.0_bp-this_arr(i)%elements(2)))**3.0_bp / &
                  planetary_mu(this_arr(i)%center))
             ! Select time of perihelion closest to epoch:
             IF (ABS(MOD(dt,p)) <= 0.5*p) THEN
                this_arr(i)%elements(6) = this_arr(i)%elements(6) + INT(dt/p)*p
             ELSE
                this_arr(i)%elements(6) = this_arr(i)%elements(6) + INT(dt/p)*p + SIGN(p,dt)
             END IF

          CASE ("keplerian")

             IF (PRESENT(jacobian) .AND. &
                  .NOT.ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                jacobian(i,:,:) = identity_matrix(6)
                mean_motion = SQRT(mu_/this_arr(i)%elements(1)**3.0_bp)
                ! The only nonzero element (apart from the diagonal) is dM/da:
                jacobian(i,6,1) = -1.5_bp * mean_motion*dt/this_arr(i)%elements(1)
             ELSE IF (PRESENT(jacobian) .AND. ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                DO j=1,6
                   elm_arr(1,:) = this_arr(i)%elements
                   elm_arr(2,:) = this_arr(i)%elements
                   elm_arr(1,j) = elm_arr(1,j) + this_arr(i)%finite_diff_prm(j)
                   elm_arr(2,j) = elm_arr(2,j) - this_arr(i)%finite_diff_prm(j)
                   mean_motion = SQRT(mu_/elm_arr(1,1)**3.0_bp)
                   elm_arr(1,6) = MODULO(elm_arr(1,6) + mean_motion*dt, two_pi)
                   mean_motion = SQRT(mu_/elm_arr(2,1)**3.0_bp)
                   elm_arr(2,6) = MODULO(elm_arr(2,6) + mean_motion*dt, two_pi)
                   jacobian(i,:,j) = (elm_arr(1,:) - elm_arr(2,:)) / &
                        (2.0_bp*this_arr(i)%finite_diff_prm(j))
                END DO
             END IF
             mean_motion = SQRT(mu_/this_arr(i)%elements(1)**3.0_bp)
             this_arr(i)%elements(6) = MODULO(this_arr(i)%elements(6) + mean_motion*dt, two_pi)

          CASE default

             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Element type cannot be propagated: " // &
                  TRIM(this_arr(i)%element_type), 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             DEALLOCATE(elm_arr, stat=err)
             RETURN

          END SELECT

          this_arr(i)%t = copy(t_)

       END DO

       IF (ASSOCIATED(elm_arr)) THEN
          DEALLOCATE(elm_arr, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Could not deallocate memory (5).", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             RETURN
          END IF
       END IF

    CASE ("n-body")

       ! Check for inconcistensies:
       DO i=2,nthis
          IF (ANY(this_arr(1)%perturbers_prm .NEQV. this_arr(i)%perturbers_prm)) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Different orbits require different perturbers.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             RETURN
          END IF
       END DO
       integrator = this_arr(1)%integrator_prm
       step = this_arr(1)%integration_step_prm
       DO i=2,nthis
          IF (integrator /= this_arr(i)%integrator_prm) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "All orbits do not share the same integrator option.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             RETURN
          END IF
          IF ((ALL(this_arr(1)%finite_diff_prm > 0.0_bp) .AND. &
               .NOT.ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) .OR. &
               (.NOT.ALL(this_arr(1)%finite_diff_prm > 0.0_bp) .AND. &
               ALL(this_arr(i)%finite_diff_prm > 0.0_bp))) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "All orbits do not share the same jacobian approach.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             RETURN
          END IF
          IF (ABS(step-this_arr(i)%integration_step_prm) > &
               100.0_bp*EPSILON(step)) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "All orbits do not share the same integration step.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             RETURN
          END IF
       END DO
       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,1X,I0)") "Number of standard perturbers:", &
               COUNT(this_arr(1)%perturbers_prm)
       END IF

       IF (ASSOCIATED(this_arr(1)%additional_perturbers)) THEN
          naddit = SIZE(this_arr(1)%additional_perturbers,dim=1)
       ELSE
          naddit = 0
       END IF

       ALLOCATE(elm_arr(6,nthis+naddit), element_type_arr(nthis), frame_arr(nthis), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "Could not allocate memory (20).", 1)
          IF (PRESENT(jacobian)) THEN
             DEALLOCATE(jacobian, stat=err)
          END IF
          DEALLOCATE(element_type_arr, stat=err)
          DEALLOCATE(frame_arr, stat=err)
          DEALLOCATE(elm_arr, stat=err)
          RETURN
       END IF
       DO i=1,nthis
          elm_arr(1:6,i) = this_arr(i)%elements
       END DO
       DO i=1,naddit
          elm_arr(1:6,nthis+i) = this_arr(1)%additional_perturbers(i,1:6)
       END DO
       IF (PRESENT(jacobian) .AND. ALL(this_arr(1)%finite_diff_prm > 0.0_bp)) THEN
          elm_arr => reallocate(elm_arr, 6, 13*(nthis+naddit))
          ! Initialize for finite difference technique:
          DO i=1,nthis+naddit
             DO j=1,6
                ! Change one element at a time by adding and subtracting by a finite amount
                elm_arr(:,nthis+naddit+(i-1)*12+j) = elm_arr(:,i)
                elm_arr(:,nthis+naddit+(i-1)*12+j+6) = elm_arr(:,i)
                elm_arr(j,nthis+naddit+(i-1)*12+j) = elm_arr(j,nthis+naddit+(i-1)*12+j) + &
                     this_arr(i)%finite_diff_prm(j)
                elm_arr(j,nthis+naddit+(i-1)*12+j+6) = elm_arr(j,nthis+naddit+(i-1)*12+j+6) - &
                     this_arr(i)%finite_diff_prm(j)
             END DO
          END DO
       ELSE IF (PRESENT(jacobian) .AND. .NOT.ALL(this_arr(1)%finite_diff_prm > 0.0_bp)) THEN
          ALLOCATE(jacobian_(6,6,nthis+naddit), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Could not allocate memory (25).", 1)
             DEALLOCATE(jacobian, stat=err)
             DEALLOCATE(element_type_arr, stat=err)
             DEALLOCATE(frame_arr, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             DEALLOCATE(jacobian_, stat=err)
             RETURN
          END IF
       END IF

       ! Make copies of original element types and frames, and make
       ! sure equatorial Cartesian elements are used during integration:
       DO i=1,nthis
          element_type_arr(i) = "cartesian"
          IF (this_arr(i)%element_type /= "cartesian") THEN
             ! Non-Cartesian elements
             element_type_arr(i) = this_arr(i)%element_type
             this_ = copy(this_arr(i))
             this_%elements = elm_arr(:,i)
             elm_arr(:,i) = getElements(this_, "cartesian", frame="equatorial")
             IF (error) THEN
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "TRACE BACK (30)", 1)
                IF (PRESENT(jacobian)) THEN
                   DEALLOCATE(jacobian, stat=err)
                END IF
                DEALLOCATE(element_type_arr, stat=err)
                DEALLOCATE(frame_arr, stat=err)
                DEALLOCATE(elm_arr, stat=err)
                DEALLOCATE(jacobian_, stat=err)
                RETURN
             END IF
             IF (PRESENT(jacobian) .AND. ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                DO j=1,12
                   this_%elements = elm_arr(:,nthis+naddit+(i-1)*12+j)
                   elm_arr(:,nthis+naddit+(i-1)*12+j) = getElements(this_, "cartesian", frame="equatorial")
                   IF (error) THEN
                      CALL errorMessage("Orbit / propagate (multiple)", &
                           "TRACE BACK (35)", 1)
                      DEALLOCATE(jacobian, stat=err)
                      DEALLOCATE(element_type_arr, stat=err)
                      DEALLOCATE(frame_arr, stat=err)
                      DEALLOCATE(elm_arr, stat=err)
                      DEALLOCATE(jacobian_, stat=err)
                      RETURN
                   END IF
                END DO
             END IF
             CALL NULLIFY(this_)
             IF (PRESENT(jacobian) .AND. .NOT.ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                IF (.NOT.ALLOCATED(partials0)) THEN
                   ALLOCATE(partials0(nthis+naddit,6,6), stat=err)
                   IF (err /= 0) THEN
                      error = .TRUE.
                      CALL errorMessage("Orbit / propagate (multiple)", &
                           "Could not allocate memory (30).", 1)
                      DEALLOCATE(jacobian, stat=err)
                      DEALLOCATE(element_type_arr, stat=err)
                      DEALLOCATE(frame_arr, stat=err)
                      DEALLOCATE(elm_arr, stat=err)
                      DEALLOCATE(jacobian_, stat=err)
                      DEALLOCATE(partials0, stat=err)
                      RETURN
                   END IF
                END IF
                IF (this_arr(i)%element_type == "keplerian") THEN
                   CALL partialsCartesianWrtKeplerian(this_arr(i), partials0(i,:,:), frame="equatorial")
                ELSE IF (this_arr(i)%element_type == "cometary") THEN
                   CALL partialsCartesianWrtCometary(this_arr(i), partials0(i,:,:), frame="equatorial")                   
                END IF
                IF (error) THEN
                   CALL errorMessage("Orbit / propagate (multiple)", &
                        "TRACE BACK (40)", 1)
                   DEALLOCATE(jacobian, stat=err)
                   DEALLOCATE(element_type_arr, stat=err)
                   DEALLOCATE(frame_arr, stat=err)
                   DEALLOCATE(elm_arr, stat=err)
                   DEALLOCATE(jacobian_, stat=err)
                   DEALLOCATE(partials0, stat=err)
                   RETURN
                END IF
             END IF
          ELSE
             ! Cartesian elements
             frame_arr(i) = this_arr(i)%frame
             IF (this_arr(i)%frame /= "equatorial") THEN
                ! Non-equatorial frame (here assumed to be ecliptic):
                CALL rotateToEquatorial(elm_arr(1:6,i))
                IF (PRESENT(jacobian) .AND. ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                   DO j=1,12
                      CALL rotateToEquatorial(elm_arr(1:6,nthis+naddit+(i-1)*12+j))
                   END DO
                END IF
             END IF
          END IF
       END DO

       IF (info_verb >= 4) THEN
          WRITE(stdout,"(2X,A,1X,I0)") "Number of additional perturbers:", naddit
       END IF
       ALLOCATE(masses(nthis+naddit), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)",&
               "Could not allocate memory (35).", 1)          
          RETURN
       END IF
       DO i=1,nthis
          masses(i) = this_arr(i)%mass_prm
       END DO
       IF (naddit /= 0) THEN
          masses(nthis+1:nthis+naddit) = this_arr(1)%additional_perturbers(1:naddit,8)
       END IF

       SELECT CASE (integrator)

       CASE ("bulirsch-stoer")

          ! Check whether the additional perturbers share the same epoch
          ! with the orbits to be integrated and integrate them to the
          ! correct epoch if required:
          IF (ASSOCIATED(this_arr(1)%additional_perturbers) .AND. &
               SIZE(this_arr(1)%additional_perturbers,dim=1) > 0) THEN
             IF (naddit >= 2) THEN
                DO i=2,naddit
                   IF (ABS(this_arr(1)%additional_perturbers(i-1,7) - &
                        this_arr(1)%additional_perturbers(i,7)) > &
                        10.0*EPSILON(this_arr(1)%additional_perturbers(1,7))) THEN
                      error = .TRUE.
                      CALL errorMessage("Orbit / propagate (multiple)",&
                           "Additional perturbers do not have a common epoch.", 1)          
                      RETURN
                   END IF
                END DO
             END IF
             IF (ABS(mjd_tt0 - this_arr(1)%additional_perturbers(1,7)) > &
                  10.0*EPSILON(this_arr(1)%additional_perturbers(1,7))) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,2(A,1X,F15.7,1X),A)") &
                        "Integrating additional perturbers from epoch", &
                        this_arr(1)%additional_perturbers(1,7), &
                        "to starting epoch", mjd_tt0, "."
                   WRITE(stdout,"(2X,A)") &
                        "Initial elements (cartesian equatorial), epoch, and mass:"
                   DO i=1,naddit
                      WRITE(stdout,"(2X,I0,7(1X,F17.9),1X,E10.4)") &
                           i, elm_arr(:,SIZE(elm_arr,dim=2)-naddit+i), &
                           this_arr(1)%additional_perturbers(i,7:8)
                   END DO
                END IF
                CALL bulirsch_full_jpl(this_arr(1)%additional_perturbers(1,7), &
                     mjd_tt0, elm_arr(:,SIZE(elm_arr,dim=2)-naddit+1:), &
                     this_arr(1)%perturbers_prm, error, &
                     step=this_arr(1)%integration_step_prm, &
                     ncenter=center, &
                     masses=this_arr(1)%additional_perturbers(:,8), &
                     info_verb=info_verb)
             END IF
             IF (info_verb >= 3) THEN
                WRITE(stdout,"(2X,A)") &
                     "Orbital elements for additional perturbers (cartesian equatorial), " // &
                     "their epoch, and mass at starting epoch:"
                DO i=1,naddit
                   WRITE(stdout,"(2X,I0,7(1X,F17.9),1X,E10.4)") &
                        i, elm_arr(:,SIZE(elm_arr,dim=2)-naddit+i), &
                        mjd_tt0, this_arr(1)%additional_perturbers(i,8)
                END DO
             END IF
          END IF

          IF (PRESENT(jacobian) .AND. .NOT.ALL(this_arr(1)%finite_diff_prm > 0.0_bp)) THEN
             ! Variational equations technique:
             DO i=1,nthis+naddit
                jacobian_(:,:,i) = identity_matrix(6)
             END DO
             IF (PRESENT(encounters)) THEN
                CALL bulirsch_full_jpl(mjd_tt0, mjd_tt, elm_arr, &
                     this_arr(1)%perturbers_prm, &
                     error, step=this_arr(1)%integration_step_prm, &
                     jacobian=jacobian_, ncenter=center, &
                     encounters=encounters, masses=masses, &
                     info_verb=info_verb)
             ELSE
                CALL bulirsch_full_jpl(mjd_tt0, mjd_tt, elm_arr, &
                     this_arr(1)%perturbers_prm, &
                     error, step=this_arr(1)%integration_step_prm, &
                     jacobian=jacobian_, ncenter=center, &
                     masses=masses, info_verb=info_verb)
             END IF
             DO i=1,nthis+naddit
                jacobian(i,1:6,1:6) = jacobian_(1:6,1:6,i)
             END DO
          ELSE IF (.NOT.PRESENT(jacobian) .OR. &
               (PRESENT(jacobian) .AND. ALL(this_arr(1)%finite_diff_prm > 0.0_bp))) THEN
             ! No Jacobian requested or Jacobian requested through
             ! finite differences technique:
             IF (PRESENT(encounters)) THEN
                CALL bulirsch_full_jpl(mjd_tt0, mjd_tt, elm_arr, &
                     this_arr(1)%perturbers_prm, error, &
                     step=this_arr(1)%integration_step_prm, &
                     ncenter=center, &
                     encounters=encounters, &
                     masses=masses, &
                     info_verb=info_verb)
             ELSE
                CALL bulirsch_full_jpl(mjd_tt0, mjd_tt, elm_arr, &
                     this_arr(1)%perturbers_prm, error, &
                     step=this_arr(1)%integration_step_prm, &
                     ncenter=center, masses=masses, &
                     info_verb=info_verb)
             END IF
          END IF
          IF (error) THEN
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Bulirsch-Stoer integration was unsuccessful.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             DEALLOCATE(element_type_arr, stat=err)
             DEALLOCATE(frame_arr, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             DEALLOCATE(jacobian_, stat=err)
             DEALLOCATE(partials0, stat=err)
             RETURN
          END IF

       CASE ("gauss-radau")

          ! Check whether the additional perturbers share the same epoch
          ! with the orbits to be integrated and integrate them to the
          ! correct epoch if required:
          IF (ALLOCATED(masses)) THEN
             IF (naddit >= 2) THEN
                DO i=2,naddit
                   IF (ABS(this_arr(1)%additional_perturbers(i-1,7) - &
                        this_arr(1)%additional_perturbers(i,7)) > &
                        10.0*EPSILON(this_arr(1)%additional_perturbers(1,7))) THEN
                      error = .TRUE.
                      CALL errorMessage("Orbit / propagate (multiple)",&
                           "Additional perturbers do not have a common epoch.", 1)          
                      RETURN
                   END IF
                END DO
             END IF
             IF (ABS(mjd_tt0 - this_arr(1)%additional_perturbers(1,7)) > &
                  10.0*EPSILON(this_arr(1)%additional_perturbers(1,7))) THEN
                IF (info_verb >= 3) THEN
                   WRITE(stdout,"(2X,2(A,1X,F15.7,1X),A)") &
                        "Integrating additional perturbers from epoch", &
                        this_arr(1)%additional_perturbers(1,7), &
                        "to starting epoch", mjd_tt0, "."
                   WRITE(stdout,"(2X,A)") &
                        "Initial elements (cartesian equatorial), epoch, and mass:"
                   DO i=1,naddit
                      WRITE(stdout,"(2X,I0,7(1X,F17.9),1X,E10.4)") &
                           i, elm_arr(:,SIZE(elm_arr,dim=2)-naddit+i), &
                           this_arr(1)%additional_perturbers(i,7:8)
                   END DO
                END IF
                CALL gauss_radau_15_full_jpl( &
                     this_arr(1)%additional_perturbers(1,7), &
                     mjd_tt0, &
                     elm_arr(:,SIZE(elm_arr,dim=2)-naddit+1:), 12, &
                     2, this_arr(1)%perturbers_prm, error, &
                     step=this_arr(1)%integration_step_prm, &
                     ncenter=center, &
                     masses=masses)
             END IF
             IF (info_verb >= 3) THEN
                WRITE(stdout,"(2X,A)") &
                     "Orbital elements for additional perturbers (cartesian equatorial), " // &
                     "their epoch, and mass at starting epoch:"
                DO i=1,naddit
                   WRITE(stdout,"(2X,I0,7(1X,F17.9),1X,E10.4)") &
                        i, elm_arr(:,SIZE(elm_arr,dim=2)-naddit+i), &
                        mjd_tt0, masses
                END DO
             END IF
          END IF

          IF (PRESENT(jacobian) .AND. .NOT.ALL(this_arr(1)%finite_diff_prm > 0.0_bp)) THEN
             ! Jacobians through variational equations technique:
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)",&
                  "Variational equations technique not available for Gauss-Radau.", 1)          
             RETURN             
             DO i=1,nthis+naddit
                jacobian_(:,:,i) = identity_matrix(6)
             END DO
             IF (PRESENT(encounters)) THEN
                CALL gauss_radau_15_full_jpl(mjd_tt0, mjd_tt, &
                     elm_arr, 12, 2, this_arr(1)%perturbers_prm, &
                     error, jacobian=jacobian_, &
                     step=this_arr(1)%integration_step_prm, &
                     ncenter=center, encounters=encounters, &
                     masses=masses)
             ELSE
                CALL gauss_radau_15_full_jpl(mjd_tt0, mjd_tt, &
                     elm_arr, 12, 2, this_arr(1)%perturbers_prm, &
                     error, jacobian=jacobian_, &
                     step=this_arr(1)%integration_step_prm, &
                     ncenter=center, &
                     masses=masses)
             END IF
             DO i=1,nthis+naddit
                jacobian(i,1:6,1:6) = jacobian_(1:6,1:6,i)
             END DO
          ELSE IF (.NOT.PRESENT(jacobian) .OR. &
               (PRESENT(jacobian) .AND. ALL(this_arr(1)%finite_diff_prm > 0.0_bp))) THEN
             ! No Jacobian requested or Jacobian requested through
             ! finite differences technique:
             IF (PRESENT(encounters)) THEN
                CALL gauss_radau_15_full_jpl(mjd_tt0, mjd_tt, &
                     elm_arr, 12, 2, this_arr(1)%perturbers_prm, &
                     error, step=this_arr(1)%integration_step_prm, &
                     ncenter=center, encounters=encounters, &
                     masses=masses)
             ELSE
                CALL gauss_radau_15_full_jpl(mjd_tt0, mjd_tt, &
                     elm_arr, 12, 2, this_arr(1)%perturbers_prm, &
                     error, step=this_arr(1)%integration_step_prm, &
                     ncenter=center, &
                     masses=masses)
             END IF
          END IF
          IF (error) THEN
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Gauss-Radau integration was unsuccessful.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             DEALLOCATE(element_type_arr, stat=err)
             DEALLOCATE(frame_arr, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             DEALLOCATE(jacobian_, stat=err)
             DEALLOCATE(partials0, stat=err)
             RETURN
          END IF

       CASE default

          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "No such integrator:" // TRIM(this_arr(1)%integrator_prm), 1)
          IF (PRESENT(jacobian)) THEN
             DEALLOCATE(jacobian, stat=err)
          END IF
          DEALLOCATE(element_type_arr, stat=err)
          DEALLOCATE(frame_arr, stat=err)
          DEALLOCATE(elm_arr, stat=err)
          DEALLOCATE(jacobian_, stat=err)
          DEALLOCATE(partials0, stat=err)
          RETURN

       END SELECT

       IF (naddit > 0 .AND. info_verb >= 3) THEN
          WRITE(stdout,"(2X,A)") &
               "Orbital elements for additional perturbers (cartesian equatorial), " // &
               "their epoch, and mass at final epoch:"
          DO i=1,naddit
             WRITE(stdout,"(2X,I0,7(1X,F17.9),1X,E10.4)") &
                  i, elm_arr(:,SIZE(elm_arr,dim=2)-naddit+i), &
                  mjd_tt, this_arr(1)%additional_perturbers(i,8)
          END DO
       END IF

       DO i=1,nthis
          this_arr(i)%t = copy(t_)
          IF (element_type_arr(i) == "cartesian" .AND. frame_arr(i) == "equatorial") THEN
             this_arr(i)%elements = elm_arr(:,i)
          ELSE IF (element_type_arr(i) == "cartesian" .AND. frame_arr(i) == "ecliptic") THEN
             CALL rotateToEcliptic(elm_arr(1:6,i))
             this_arr(i)%elements = elm_arr(:,i)
             IF (PRESENT(jacobian) .AND. ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                DO j=1,12
                   CALL rotateToEcliptic(elm_arr(1:6,nthis+naddit+(i-1)*12+j))
                END DO
             END IF
          ELSE IF (element_type_arr(i) == "keplerian") THEN
             CALL NEW(this_, elm_arr(1:6,i), "cartesian", "equatorial", copy(t_))
             IF (error) THEN
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "TRACE BACK (45)", 1)
                IF (PRESENT(jacobian)) THEN
                   DEALLOCATE(jacobian, stat=err)
                END IF
                DEALLOCATE(element_type_arr, stat=err)
                DEALLOCATE(frame_arr, stat=err)
                DEALLOCATE(elm_arr, stat=err)
                DEALLOCATE(jacobian_, stat=err)
                DEALLOCATE(partials0, stat=err)
                RETURN
             END IF
             CALL toKeplerian(this_)
             IF (error) THEN
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "TRACE BACK (50)", 1)
                IF (PRESENT(jacobian)) THEN
                   DEALLOCATE(jacobian, stat=err)
                END IF
                DEALLOCATE(element_type_arr, stat=err)
                DEALLOCATE(frame_arr, stat=err)
                DEALLOCATE(elm_arr, stat=err)
                DEALLOCATE(jacobian_, stat=err)
                DEALLOCATE(partials0, stat=err)
                RETURN
             END IF
             this_arr(i)%elements = this_%elements
             CALL NULLIFY(this_)
             IF (PRESENT(jacobian) .AND. ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                DO j=1,12
                   CALL NEW(this_, elm_arr(1:6,nthis+naddit+(i-1)*12+j), "cartesian", "equatorial", copy(t_))
                   IF (error) THEN
                      CALL errorMessage("Orbit / propagate (multiple)", &
                           "TRACE BACK (55)", 1)
                      DEALLOCATE(jacobian, stat=err)
                      DEALLOCATE(element_type_arr, stat=err)
                      DEALLOCATE(frame_arr, stat=err)
                      DEALLOCATE(elm_arr, stat=err)
                      DEALLOCATE(jacobian_, stat=err)
                      DEALLOCATE(partials0, stat=err)
                      RETURN
                   END IF
                   elm_arr(1:6,nthis+naddit+(i-1)*12+j) = getElements(this_, "keplerian")
                   IF (error) THEN
                      CALL errorMessage("Orbit / propagate (multiple)", &
                           "TRACE BACK (60)", 1)
                      DEALLOCATE(jacobian, stat=err)
                      DEALLOCATE(element_type_arr, stat=err)
                      DEALLOCATE(frame_arr, stat=err)
                      DEALLOCATE(elm_arr, stat=err)
                      DEALLOCATE(jacobian_, stat=err)
                      DEALLOCATE(partials0, stat=err)
                      RETURN
                   END IF
                   CALL NULLIFY(this_)
                END DO
             END IF
          ELSE IF (element_type_arr(i) == "cometary") THEN
             CALL NEW(this_, elm_arr(1:6,i), "cartesian", "equatorial", copy(t_))
             IF (error) THEN
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "TRACE BACK (45)", 1)
                IF (PRESENT(jacobian)) THEN
                   DEALLOCATE(jacobian, stat=err)
                END IF
                DEALLOCATE(element_type_arr, stat=err)
                DEALLOCATE(frame_arr, stat=err)
                DEALLOCATE(elm_arr, stat=err)
                DEALLOCATE(jacobian_, stat=err)
                DEALLOCATE(partials0, stat=err)
                RETURN
             END IF
             CALL toCometary(this_)
             IF (error) THEN
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "TRACE BACK (50)", 1)
                IF (PRESENT(jacobian)) THEN
                   DEALLOCATE(jacobian, stat=err)
                END IF
                DEALLOCATE(element_type_arr, stat=err)
                DEALLOCATE(frame_arr, stat=err)
                DEALLOCATE(elm_arr, stat=err)
                DEALLOCATE(jacobian_, stat=err)
                DEALLOCATE(partials0, stat=err)
                RETURN
             END IF
             this_arr(i)%elements = this_%elements
             CALL NULLIFY(this_)
             IF (PRESENT(jacobian) .AND. ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                DO j=1,12
                   CALL NEW(this_, elm_arr(1:6,nthis+naddit+(i-1)*12+j), "cartesian", "equatorial", copy(t_))
                   IF (error) THEN
                      CALL errorMessage("Orbit / propagate (multiple)", &
                           "TRACE BACK (55)", 1)
                      DEALLOCATE(jacobian, stat=err)
                      DEALLOCATE(element_type_arr, stat=err)
                      DEALLOCATE(frame_arr, stat=err)
                      DEALLOCATE(elm_arr, stat=err)
                      DEALLOCATE(jacobian_, stat=err)
                      DEALLOCATE(partials0, stat=err)
                      RETURN
                   END IF
                   elm_arr(1:6,nthis+naddit+(i-1)*12+j) = getElements(this_, "cometary")
                   IF (error) THEN
                      CALL errorMessage("Orbit / propagate (multiple)", &
                           "TRACE BACK (60)", 1)
                      DEALLOCATE(jacobian, stat=err)
                      DEALLOCATE(element_type_arr, stat=err)
                      DEALLOCATE(frame_arr, stat=err)
                      DEALLOCATE(elm_arr, stat=err)
                      DEALLOCATE(jacobian_, stat=err)
                      DEALLOCATE(partials0, stat=err)
                      RETURN
                   END IF
                   CALL NULLIFY(this_)
                END DO
             END IF
          ELSE 
             error =.TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "'" // TRIM(element_type_arr(i)) // &
                  "' elements not supported.", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             DEALLOCATE(element_type_arr, stat=err)
             DEALLOCATE(frame_arr, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             DEALLOCATE(jacobian_, stat=err)
             DEALLOCATE(partials0, stat=err)
             RETURN
          END IF

          IF (PRESENT(jacobian)) THEN
             IF (ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                ! If fd approach, compute the Jacobian in the following way 
                ! for all element types and frames:
                DO j=1,6
                   jacobian(i,:,j) = (elm_arr(:,nthis+naddit+(i-1)*12+j) - &
                        elm_arr(:,nthis+naddit+(i-1)*12+j+6)) / &
                        (2.0_bp*this_arr(i)%finite_diff_prm(j))
                END DO
             ELSE IF ((element_type_arr(i) == "keplerian" .OR. &
                  element_type_arr(i) == "cometary") .AND. &
                  .NOT.ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                ! If element type is Keplerian and ve approach:
                IF (element_type_arr(i) == "keplerian") THEN
                   CALL partialsKeplerianWrtCartesian(this_arr(i), &
                        partials1, frame="equatorial")
                ELSE IF (element_type_arr(i) == "cometary") THEN
                   CALL partialsCometaryWrtCartesian(this_arr(i), &
                        partials1, frame="equatorial")                   
                END IF
                IF (error) THEN
                   CALL errorMessage("Orbit / propagate (multiple)", &
                        "TRACE BACK (65)", 1)
                   DEALLOCATE(jacobian, stat=err)
                   DEALLOCATE(element_type_arr, stat=err)
                   DEALLOCATE(frame_arr, stat=err)
                   DEALLOCATE(elm_arr, stat=err)
                   DEALLOCATE(jacobian_, stat=err)
                   DEALLOCATE(partials0, stat=err)
                   RETURN
                END IF
                jacobian(i,:,:) = MATMUL(MATMUL(partials1,jacobian(i,:,:)),partials0(i,:,:))
             ELSE IF (element_type_arr(i) == "cartesian" .AND. frame_arr(i) == "ecliptic" .AND. &
                  .NOT.ALL(this_arr(i)%finite_diff_prm > 0.0_bp)) THEN
                ! If element type is Cartesian, frame is ecliptical and ve approach:
                DO j=1,6
                   CALL rotateToEcliptic(jacobian(i,:,j))
                END DO
                DO j=1,6
                   CALL rotateToEcliptic(jacobian(i,j,:))
                END DO
             ELSE IF (element_type_arr(i) /= "cartesian" .AND. &
                  element_type_arr(i) /= "cometary" .AND. &
                  element_type_arr(i) /= "keplerian") THEN 
                error =.TRUE.
                CALL errorMessage("Orbit / propagate (multiple)", &
                     "Jacobian of '" // TRIM(element_type_arr(i)) // &
                     "' elements not supported.", 1)
                DEALLOCATE(element_type_arr, stat=err)
                DEALLOCATE(frame_arr, stat=err)
                DEALLOCATE(elm_arr, stat=err)
                DEALLOCATE(jacobian, stat=err)
                DEALLOCATE(jacobian_, stat=err)
                DEALLOCATE(partials0, stat=err)
                RETURN
             END IF
          END IF
          ! Update additional perturbers' orbital elements and epochs
          DO j=1,naddit
             this_arr(i)%additional_perturbers(j,1:6) = elm_arr(:,nthis+j)
             this_arr(i)%additional_perturbers(j,7) = mjd_tt
          END DO
       END DO

       IF (ALLOCATED(partials0)) THEN
          DEALLOCATE(partials0, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Could not deallocate memory (10).", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             DEALLOCATE(element_type_arr, stat=err)
             DEALLOCATE(frame_arr, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             DEALLOCATE(jacobian_, stat=err)
             RETURN
          END IF
       END IF
       IF (ALLOCATED(jacobian_)) THEN
          DEALLOCATE(jacobian_, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / propagate (multiple)", &
                  "Could not deallocate memory (15).", 1)
             IF (PRESENT(jacobian)) THEN
                DEALLOCATE(jacobian, stat=err)
             END IF
             DEALLOCATE(element_type_arr, stat=err)
             DEALLOCATE(frame_arr, stat=err)
             DEALLOCATE(elm_arr, stat=err)
             RETURN
          END IF
       END IF
       DEALLOCATE(elm_arr, element_type_arr, frame_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / propagate (multiple)", &
               "Could not deallocate memory (20).", 1)
          IF (PRESENT(jacobian)) THEN
             DEALLOCATE(jacobian, stat=err)
          END IF
          DEALLOCATE(element_type_arr, stat=err)
          DEALLOCATE(frame_arr, stat=err)
          DEALLOCATE(elm_arr, stat=err)
          RETURN
       END IF

    CASE default

       error = .TRUE.
       CALL errorMessage("Orbit / propagate (multiple)", &
            "Could not choose between 2-body and n-body: " // &
            TRIM(dyn_model), 1)
       IF (PRESENT(jacobian)) THEN
          DEALLOCATE(jacobian, stat=err)
       END IF
       DEALLOCATE(element_type_arr, stat=err)
       DEALLOCATE(frame_arr, stat=err)
       DEALLOCATE(elm_arr, stat=err)
       DEALLOCATE(jacobian_, stat=err)
       DEALLOCATE(partials0, stat=err)
       RETURN

    END SELECT

    CALL NULLIFY(this_)
    CALL NULLIFY(t_)

  END SUBROUTINE propagate_Orb_multiple





  !! *Description*:
  !!
  !! Reallocates a pointer array of Orbit-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myorbits => reallocate(myorbits,4)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_Orb_1(array,n)

    IMPLICIT NONE
    TYPE (Orbit), DIMENSION(:), POINTER :: reallocate_Orb_1, array
    INTEGER, INTENT(in)                 :: n
    INTEGER                             :: i, nold, err

    ALLOCATE(reallocate_Orb_1(n), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_Orb_1 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    DO i=1, MIN(n,nold)
       CALL NULLIFY(reallocate_Orb_1(i))
       reallocate_Orb_1(i) = copy(array(i))
    END DO
    DO i=1,nold
       CALL NULLIFY(array(i))
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_Orb_1





  !! *Description*:
  !!
  !! Reallocates a pointer array of Orbit-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! myorbits => reallocate(myorbits,4,2)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_Orb_2(array, n, m)

    IMPLICIT NONE
    TYPE (Orbit), DIMENSION(:,:), POINTER :: reallocate_Orb_2, array
    INTEGER, INTENT(in)                   :: n, m
    INTEGER                               :: i, j, nold, mold, err

    ALLOCATE(reallocate_Orb_2(n,m), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_Orb_2 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array,dim=1)
    mold = SIZE(array,dim=2)
    DO i=1, MIN(n,nold)
       DO j=1, MIN(m,mold)
          CALL NULLIFY(reallocate_Orb_2(i,j))
          reallocate_Orb_2(i,j) = copy(array(i,j))
       END DO
    END DO
    DO i=1,nold
       DO j=1,mold
          CALL NULLIFY(array(i,j))
       END DO
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_Orb_2





  !! *Description*:
  !!
  !! Reallocates a pointer array of Orbit-objects by removing
  !! non-initialized objects.
  !!
  !! *Usage*:
  !!
  !! myorbits => reallocate(myorbits)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_Orb(array)

    IMPLICIT NONE
    TYPE (Orbit), DIMENSION(:), POINTER :: reallocate_Orb, array
    INTEGER                             :: i, j, err

    IF (.NOT. ASSOCIATED(array)) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Parameter pointer is not associated.", 1)
       reallocate_Orb => NULL()
       RETURN
    END IF

    ALLOCATE(reallocate_Orb(SIZE(array,dim=1)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_Orb => NULL()
       RETURN
    END IF
    j = 0
    DO i=1, SIZE(array,dim=1)
       IF (exist(array(i))) THEN
          j = j + 1 
          reallocate_Orb(j) = copy(array(i))
       END IF
    END DO
    DO i=1,SIZE(array)
       CALL NULLIFY(array(i))
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF
    reallocate_Orb => reallocate(reallocate_Orb,j)

  END FUNCTION reallocate_Orb





  !! *Description*:
  !!
  !! Rotates coordinates to the ecliptical coordinate frame, if not
  !! already the case.
  !!
  SUBROUTINE rotateToEcliptic_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout) :: this

    IF (this%frame /= "ecliptic" .AND. this%is_initialized &
         .AND. this%element_type == "cartesian") THEN
       CALL rotateToEcliptic(this%elements)
       this%frame = "ecliptic"
    END IF

  END SUBROUTINE rotateToEcliptic_Orb





  !! *Description*:
  !!
  !! Rotates coordinates to the equatorial coordinate frame, if not
  !! already the case.
  !!
  SUBROUTINE rotateToEquatorial_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout) :: this

    IF (this%frame /= "equatorial" .AND. this%is_initialized &
         .AND. this%element_type == "cartesian") THEN
       CALL rotateToEquatorial(this%elements)
       this%frame = "equatorial"
    END IF

  END SUBROUTINE rotateToEquatorial_Orb





  !! *Description*:
  !!
  !! Use finite diffference approach instead of the default
  !! variational equations in n-body integrations. Note that also the
  !! setParameters routine must be called to force n-body integration
  !! to be used instead of the default 2-body approximations.
  !!
  !! Returns error.
  !!
  SUBROUTINE setParameters_Orb(this, &
       mass, &
       dyn_model, &
       integration_step, &
       integrator, &
       finite_diff, &
       perturbers, &
       additional_perturbers)

    IMPLICIT NONE
    TYPE(Orbit), INTENT(inout)                   :: this
    CHARACTER(len=*), INTENT(in), OPTIONAL       :: dyn_model
    CHARACTER(len=*), INTENT(in), OPTIONAL       :: integrator
    REAL(bp), DIMENSION(:,:), OPTIONAL           :: additional_perturbers ! (car equ + mjd_tt + mass)
    REAL(bp), DIMENSION(6), INTENT(in), OPTIONAL :: finite_diff
    REAL(bp), INTENT(in), OPTIONAL               :: integration_step
    REAL(bp), INTENT(in), OPTIONAL               :: mass
    LOGICAL, DIMENSION(10), OPTIONAL             :: perturbers

    CHARACTER(len=256) :: str
    INTEGER :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / setParameters", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(mass)) THEN
       this%mass_prm = mass
    END IF
    IF (PRESENT(dyn_model)) THEN
       IF (LEN_TRIM(dyn_model) > DYN_MODEL_LEN) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / setParameters", &
               "Parameter 'propagation' too long.", 1)
          RETURN
       END IF
       str = dyn_model
       CALL locase(str, error)
       IF (error) THEN
          CALL errorMessage("Orbit / setParameters", &
               "The dynamical model string contains forbidden characters.", 1)
          RETURN
       END IF
       IF (str /= "2-body" .AND. &
            str /= "n-body") THEN
          error = .TRUE.
          CALL errorMessage("Orbit / setParameters", &
               "Option " // TRIM(dyn_model) // " not available.", 1)
          RETURN
       END IF
       this%dyn_model_prm = TRIM(str)
    END IF
    IF (PRESENT(integration_step)) THEN
       this%integration_step_prm = integration_step
    END IF
    IF (PRESENT(integrator)) THEN
       str = integrator
       CALL locase(str, error)
       IF (error) THEN
          CALL errorMessage("Orbit / setParameters", &
               "The integrator string contains forbidden characters.", 1)
          RETURN
       END IF
       this%integrator_prm = TRIM(str)
    END IF
    IF (PRESENT(finite_diff)) THEN
       ! Only use the finite difference technique, if all finite
       ! differences are positive.
       IF (ALL(finite_diff > 0.0_bp)) THEN
          this%finite_diff_prm = finite_diff
       ELSE
          this%finite_diff_prm = -1.0_bp
       END IF
    END IF
    IF (PRESENT(perturbers)) THEN
       this%perturbers_prm = perturbers
    END IF
    IF (PRESENT(additional_perturbers)) THEN
       IF (SIZE(additional_perturbers,dim=2) < 8) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / setParameters", &
               "Too few parameters provided for additional perturbers. " // &
               "Minimum information includes orbital elements, epoch, and perturber mass.", 1)
          RETURN          
       END IF
       IF (.NOT.ASSOCIATED(this%additional_perturbers)) THEN
          ALLOCATE(this%additional_perturbers(SIZE(additional_perturbers,dim=1), &
               SIZE(additional_perturbers,dim=2)), stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / setParameters", &
                  "Could not allocate memory.", 1)
             RETURN
          END IF
       END IF
       this%additional_perturbers = additional_perturbers
    END IF

  END SUBROUTINE setParameters_Orb




  !! *Description*:
  !!
  !! Iterates Kepler's equation for input Keplerian orbital elements
  !! and epoch with accelerated Newton's method, and returns the 
  !! eccentric anomaly at the given time. If an error occurs, the
  !! subroutine returns a negative eccentric anomaly.
  !!
  !! References:
  !!
  !! J. M. A. Danby: "Fundamentals of Celestial mechanics"
  !! pp. 149-154
  !!
  !! Returns error.
  !!
  SUBROUTINE solveKeplerEquation_newton(this, t, ea)

    IMPLICIT NONE 
    TYPE(Orbit), INTENT(in) :: this
    TYPE(Time), INTENT(in)  :: t
    REAL(bp), INTENT(out)   :: ea
    INTEGER, PARAMETER      :: nmax = 10000
    TYPE(Time)              :: t_
    REAL(bp), PARAMETER     :: k = 0.85_bp
    REAL(bp), PARAMETER     :: tol = 1.0e-14_bp
    REAL(bp), DIMENSION(6)  :: elements
    REAL(bp)                :: mjd_tt0, mjd_tt, ma, sigma, x, dx, &
         esinx, ecosx, f, fp, fpp, fppp
    INTEGER                 :: i

    ea = -HUGE(ea)

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / solveKeplerEquation", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%element_type == "keplerian") THEN
       elements = this%elements
    ELSE
       elements = getElements(this, "keplerian")
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "TRACE BACK 5", 1)
          RETURN
       END IF
    END IF
    t_ = copy(this%t)
    mjd_tt0 = getMJD(t_, "TT")
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / solveKeplerEquation", &
            "TRACE BACK 10", 1)
       RETURN
    END IF
    CALL NULLIFY(t_)
    t_ = copy(t)
    mjd_tt = getMJD(t_, "TT")
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / solveKeplerEquation", &
            "TRACE BACK 15", 1)
       RETURN
    END IF
    CALL NULLIFY(t_)

    ! Find initial guess:
    ma = elements(6) + SQRT(planetary_mu(this%center)/elements(1)**3.0_bp)*(mjd_tt-mjd_tt0)
    ma = MODULO(ma,two_pi)
    sigma = SIGN(1.0_bp,SIN(ma))
    x = ma + sigma*k*elements(2)

    ! Solve Kepler's equation iteratively using Newton's accelerated method:
    i = 1
    esinx = elements(2)*SIN(x)
    f = x - esinx - ma
    DO WHILE (ABS(f) >= tol)
       IF (i > nmax) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "Uncertain convergence.", 1)
          RETURN
       END IF
       ecosx = elements(2)*COS(x)
       fp    = 1.0_bp - ecosx
       fpp   = esinx
       fppp  = ecosx
       dx    = -f/fp
       dx    = -f/(fp+0.5_bp*dx*fpp)
       dx    = -f/(fp+0.5_bp*dx*fpp+dx*dx*fppp/6.0_bp)
       x     = x + dx
       esinx = elements(2)*SIN(x)
       f     = x - esinx - ma
       i     = i + 1
    END DO

    ea = MODULO(x,two_pi)

  END SUBROUTINE solveKeplerEquation_newton





  !! *Description*:
  !!
  !! Solves Kepler's equation using universal variables and returns the values 
  !! of the scaled Stumpff functions c1, c2, and c3. 
  !!
  !! Returns error.
  !!
  !! Reference: 
  !! Danby, J. M.: Fundamentals of Celestial Mechanics, 2nd ed., 
  !!               rev. ed. Richmond, VA: Willmann-Bell, pp. 170-178, 1988
  !!
  SUBROUTINE solveKeplerEquation_stumpff(r0, u, alpha, dt, stumpff_cs, ffs, s, center)

    IMPLICIT NONE
    REAL(bp), INTENT(in)                  :: r0, u, alpha
    REAL(bp), INTENT(inout)               :: dt
    REAL(bp), DIMENSION(0:3), INTENT(out) :: stumpff_cs, ffs
    INTEGER, INTENT(in)                   :: center
    REAL(bp), INTENT(out)                 :: s
    INTEGER,  PARAMETER                   :: nnew = 40 ! 8
    INTEGER,  PARAMETER                   :: nlag = 100 ! 20
    ! Care must be taken that tols is not too small!
    REAL(bp), PARAMETER                   :: tols = 1.0e-11_bp ! e-10
    REAL(bp), PARAMETER                   :: toll = 1.0e+150_bp
    REAL(bp), DIMENSION(0:3)              :: stumpff_c
    REAL(bp)                              :: s_, ds, a, e, en, ec, es, tmp
    REAL(bp)                              :: ech, esh, x, y, sigma, dm, ln, mu_
    INTEGER                               :: n

    ! Define mu parameter
    mu_ = planetary_mu(center)

    ! For small ds, equal initial guesses for elliptic and hyperbolic motion;
    ! for large ds, separate initial guesses:
    IF (ABS(dt/r0) <= 0.2_bp) THEN
       IF (info_verb >= 4) THEN
          WRITE(stdout,*) "Small ds -> equal initial guess for e<1 and e>1 orbits..."
       END IF
       IF (ABS(r0**3.0_bp) < 10*EPSILON(r0)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "r0**3 equal to zero.", 1)
          RETURN
       ELSE
          s = dt/r0 - dt**2.0_bp * 0.5_bp*u/(r0**3.0_bp) ! r0 can be small !!!
       END IF
    ELSE
       a = mu_/alpha
       IF (ABS(a**3.0_bp) < 10**(-15)*EPSILON(tmp)) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "a**3 equal to zero.", 1)
          RETURN
       ELSE
          ! If alpha < 0, then mu < 0, and vice versa:
          en = SQRT(SIGN(mu_,alpha)/(a**3.0_bp))
       END IF
       IF (alpha > EPSILON(alpha)) THEN ! alpha positive -> elliptic motion
          IF (info_verb >= 4) THEN
             WRITE(stdout,*) "Large ds and alpha > 0 -> initial guess for e<1 orbit..."
          END IF
          ec = 1.0_bp - r0/a
          tmp = en*a**2.0_bp
          IF (ABS(tmp) < 10*EPSILON(tmp)) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / solveKeplerEquation", &
                  "1-r0/a equal to zero (5).", 1)
             RETURN
          ELSE
             es = u/tmp
          END IF
          e = SQRT(ec**2.0_bp + es**2.0_bp)
          dt = dt - INT(en*dt/two_pi) * two_pi/en
          y = en*dt - es
          sigma = SIGN(1.0_bp,es*COS(y)+ec*SIN(y))
          x = y + sigma*0.85_bp*e
          s = x/SQRT(alpha)
       ELSE IF (alpha < -EPSILON(alpha)) THEN ! alpha negative -> hyperbolic motion
          IF (info_verb >= 4) THEN
             WRITE(stdout,*) "Large ds and alpha < 0 -> initial guess for e>1 orbit..."
          END IF
          ech = 1.0_bp - r0/a
          tmp = SQRT(-a*mu_)
          esh = u/tmp
          tmp = ech**2.0_bp - esh**2.0_bp
          IF (tmp < 0.0_bp) THEN
             error = .TRUE.
             CALL errorMessage("Orbit / solveKeplerEquation", &
                  "ech**2.0_bp - esh**2.0_bp negative.", 1)
             RETURN
          ELSE
             e = SQRT(tmp)
          END IF
          dm = en*dt
          IF (dm >= 0.0_bp) THEN
             tmp = ech + esh
             s = LOG((2.0_bp*dm + 1.8_bp*e)/tmp)/SQRT(-alpha)
          ELSE
             tmp = ech - esh
             s = -LOG((-2.0_bp*dm + 1.8_bp*e)/tmp)/SQRT(-alpha)
          END IF
       ELSE
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "Large ds and alpha ~0 -> no initial guess for s...", 1)
          RETURN
       END IF
    END IF
    s_ = s

    ! Newton's method for solving s:
    n  = 0
    DO WHILE (n < nnew)
       n             = n + 1
       !write(*,*) "nnew", n, s, alpha, s**2.0_bp, s**2.0_bp * alpha
       x             = s**2.0_bp * alpha
       CALL getStumpffFunctions(x, stumpff_c)
       IF (error) THEN
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "TRACE BACK 5", 1)
          RETURN
       END IF
       stumpff_cs(0) = stumpff_c(0)
       stumpff_cs(1) = stumpff_c(1) * s
       stumpff_cs(2) = stumpff_c(2) * s**2.0_bp
       stumpff_cs(3) = stumpff_c(3) * s**3.0_bp
       ffs(0) = r0*stumpff_cs(1) + u*stumpff_cs(2) + mu_*stumpff_cs(3) - dt
       ffs(1) = r0*stumpff_cs(0) + u*stumpff_cs(1) + mu_*stumpff_cs(2)
       ffs(2) = (-r0*alpha+mu_)*stumpff_cs(1) + u*stumpff_cs(0)
       ffs(3) = (-r0*alpha+mu_)*stumpff_cs(0) - u*alpha*stumpff_cs(1)
       IF (ABS(ffs(1)) < 10*EPSILON(ffs(1))) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "Division by zero (40).", 1)
          RETURN
       ELSE
          ds = -ffs(0) / ffs(1)
       END IF
       tmp = ffs(1) + 0.5_bp*ds*ffs(2)
       IF (ABS(tmp) < 10*EPSILON(ffs(1))) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "Division by zero (45).", 1)
          RETURN
       ELSE
          ds = -ffs(0) / tmp
       END IF
       tmp = ffs(1) + 0.5_bp*ds*ffs(2) + ds**2.0_bp*ffs(3)/6.0_bp
       IF (ABS(tmp) < 10*EPSILON(ffs(1))) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "Division by zero (50).", 1)
          RETURN
       ELSE
          ds = -ffs(0) / tmp
       END IF
       s = s + ds
       ! Return if/when solution is found:
       IF (ABS(ds) < tols) THEN
          RETURN
       END IF
    END DO

    ! Laguerre's method for solving s:
    s  = s_
    ln = 5.0_bp
    n  = 0
    DO WHILE (n < nlag)
       n = n + 1
       !write(*,*) "nlag", n, s, r0, u, alpha, dt, center, stumpff_cs(0:3), ffs(0:2)
       !write(*,*) "nlag", n, s, stumpff_cs(0:3), ffs(0:2)
       x = s**2.0_bp * alpha
       CALL getStumpffFunctions(x, stumpff_c)
       IF (error) THEN
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "TRACE BACK 10", 1)
          RETURN
       END IF
       stumpff_cs(0) = stumpff_c(0)
       stumpff_cs(1) = stumpff_c(1) * s
       stumpff_cs(2) = stumpff_c(2) * s**2.0_bp
       stumpff_cs(3) = stumpff_c(3) * s**3.0_bp
       ffs(0)        = r0*stumpff_cs(1) + u*stumpff_cs(2) + mu_*stumpff_cs(3) - dt
       ffs(1)        = r0*stumpff_cs(0) + u*stumpff_cs(1) + mu_*stumpff_cs(2)
       ffs(2)        = (-r0*alpha+mu_)*stumpff_cs(1) + u*stumpff_cs(0)
       tmp = ffs(1)+SIGN(1.0_bp,ffs(1)) * &
            SQRT(ABS((ln-1.0_bp)**2.0_bp*ffs(1)**2.0_bp - &
            (ln-1.0_bp)*ln*ffs(0)*ffs(2)))
       IF (ABS(tmp) < 10*EPSILON(ffs(1))) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "Division by zero (55).", 1)
          RETURN
       ELSE IF (ABS(ffs(1)) > toll) THEN
          error = .TRUE.
          CALL errorMessage("Orbit / solveKeplerEquation", &
               "tmp blew up.", 1)
          RETURN
       ELSE
          ds = -ln*ffs(0) / tmp
       END IF
       s = s + ds
       ! Return if/when solution is found:
       IF (ABS(ds) < tols) THEN
          RETURN
       END IF
    END DO

    ! Both Newton's and Laguerre's methods were unsuccesful. Report error:
    error = .TRUE.
    CALL errorMessage("Orbit / solveKeplerEquation", &
         "Newton's and Laguerre's methods were unsuccessful.", 1)

  END SUBROUTINE solveKeplerEquation_stumpff




  SUBROUTINE switchCenter_Orb(this, center)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout) :: this
    INTEGER, INTENT(in) :: center

    ! Scratch
    TYPE (Time) ::  t
    TYPE(CartesianCoordinates) :: ccoord
    CHARACTER(len=ELEMENT_TYPE_LEN) :: element_type
    CHARACTER(len=FRAME_LEN) :: frame
    REAL(bp), DIMENSION(:,:), POINTER :: planeph => NULL()
    REAL(bp), DIMENSION(6) :: helioc_elements, plan_elements
    REAL(bp) :: mjd_tt

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / switchCenter", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (center == this%center) THEN
       RETURN
    END IF

    IF (center > 18) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / switchCenter", &
            "An input new center is not specified.", 1)
       RETURN
    END IF

    element_type = this%element_type
    frame = this%frame
    t = getTime(this)
    mjd_tt = getMJD(t, "TT")    

    ! Compute heliocentric equivalent of the old elements
    IF (this%center /= 11) THEN
       ! old elements not heliocentric
       planeph => JPL_ephemeris(mjd_tt, this%center, 11, error)
       CALL NEW(ccoord, planeph(1,1:6), "equatorial", t)
       CALL rotateToEcliptic(ccoord)
       plan_elements = getCoordinates(ccoord)
       CALL NULLIFY(ccoord)
       DEALLOCATE(planeph)   
       helioc_elements = plan_elements + getElements(this, "cartesian", "ecliptic")
    ELSE IF (this%center == 11) THEN
       ! old elements heliocentric
       helioc_elements = getElements(this, "cartesian", "ecliptic")
    END IF

    IF (center == 11) THEN
       plan_elements = 0.0_bp
    ELSE
       ! coordinates for the new center
       planeph => JPL_ephemeris(mjd_tt, center, 11, error)
       CALL NEW(ccoord, planeph(1,1:6), "equatorial", t)
       CALL rotateToEcliptic(ccoord)
       plan_elements = getCoordinates(ccoord)
       CALL NULLIFY(ccoord)
       DEALLOCATE(planeph)   
    END IF

    ! change from heliocentric equivalent of the old elements to the
    ! elements of the center
    this%center = center
    this%frame = "ecliptic"
    this%element_type = "cartesian"
    this%elements = helioc_elements - plan_elements

    IF (frame == "equatorial") THEN
       CALL rotateToEquatorial(this)       
    END IF

    ! Return to original element type
    IF (element_type == "keplerian") THEN
       CALL toKeplerian(this)
    ELSE IF (element_type == "cometary") THEN
       CALL toCometary(this)
    END IF

  END SUBROUTINE switchCenter_Orb





  LOGICAL FUNCTION checkCapture_Orb(this, forCenter) 

    IMPLICIT NONE
    ! Input interface
    TYPE (Orbit), INTENT(inout) :: this
    INTEGER, INTENT(in) :: forCenter

    ! PRIVATE
    TYPE (Orbit) :: this_
    REAL(bp) :: specific_energy
    ! BEGIN

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / checkCapture", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF
    this_ = copy(this)

    IF (forCenter > 18) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / checkCapture", &
            "An input center for which the check is performed is not specified.", 1)
       RETURN
    END IF

    ! If cartesian, switch to cometary
    IF (this_%element_type /= "cartesian") THEN
       CALL toCartesian(this_, this_%frame)
    END IF

    ! Input orbit
    ! Check center, if not same as forcenter, switchCenter to the right one
    IF (this_%center /= forCenter) THEN
       !WRITE(*,*) 'Switching center...'
       CALL switchCenter(this_, forCenter)
    END IF

    ! Negative specific energy implies a bound orbit
    specific_energy = 0.5_bp * &
         DOT_PRODUCT(this_%elements(4:6),this_%elements(4:6)) - &
         planetary_mu(this_%center) / &
         SQRT(DOT_PRODUCT(this_%elements(4:6),this_%elements(4:6)))
    IF (specific_energy < 0.0_bp) THEN
       checkCapture_Orb = .TRUE.
    ELSE
       checkCapture_Orb = .FALSE.        
    END IF

    CALL NULLIFY(this_)

  END FUNCTION checkCapture_Orb





  SUBROUTINE toCartesian_Orb(this, frame)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)  :: this
    CHARACTER(len=*), INTENT(in) :: frame
    CHARACTER(len=FRAME_LEN)     :: frame_
    REAL(bp), DIMENSION(6)       :: elements

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / toCartesian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    frame_ = frame
    CALL locase(frame_, error)
    IF (error) THEN
       CALL errorMessage("Orbit / toCartesian", &
            "The frame string contains forbidden characters.", 1)
       RETURN
    END IF
    IF (frame_ /= "equatorial" .AND. &
         frame_ /= "ecliptic") THEN
       error = .TRUE.
       CALL errorMessage("Orbit / toCartesian", &
            "Frame " // TRIM(frame_) // " not recognized.", 1)
    END IF

    IF (info_verb >= 4) THEN
       WRITE(stdout,"(2X,A,1X,A)") &
            "Orbit / toCartesian:", &
            "Conversion to " // TRIM(frame_) // &
            " Cartesian elements. Initial " // TRIM(this%frame) // &
            " " // TRIM(this%element_type) // " elements:"
       IF (this%element_type == "cartesian") THEN
          WRITE(stdout,"(6(F22.15))") this%elements
       ELSE IF (this%element_type == "cometary") THEN
          WRITE(stdout,"(6(F22.15))") this%elements(1:2), &
               this%elements(3:5)/rad_deg, this%elements(6)
       ELSE IF (this%element_type == "keplerian") THEN
          WRITE(stdout,"(6(F22.15))") this%elements(1:2), &
               this%elements(3:6)/rad_deg
       END IF
    END IF

    IF (this%element_type /= "cartesian") THEN
       elements = getElements(this, "cartesian", frame_)
       IF (error) THEN
          CALL errorMessage("Orbit / toCartesian", &
               "TRACE BACK (5).", 1)
          RETURN
       END IF
       this%elements = elements
       this%element_type = "cartesian"
       this%frame = frame_
    ELSE
       IF (frame_ == "equatorial") THEN
          CALL rotateToEquatorial(this)
       ELSE
          CALL rotateToEcliptic(this)
       END IF
    END IF

    IF (info_verb >= 4) THEN
       WRITE(stdout,"(2X,A,1X,A)") &
            "Orbit / toCartesian:", &
            "Final " // TRIM(this%frame) // &
            " " // TRIM(this%element_type) // " elements:"
       WRITE(stdout,"(6(F22.15))") this%elements
    END IF

  END SUBROUTINE toCartesian_Orb





  SUBROUTINE toCometary_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)  :: this
    REAL(bp), DIMENSION(6)       :: elements

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / toCometary", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this, "cometary")
    IF (error) THEN
       CALL errorMessage("Orbit / toCometary", &
            "TRACE BACK (5).", 1)
       RETURN
    END IF
    this%elements = elements
    this%element_type = "cometary"
    this%frame = "ecliptic"

  END SUBROUTINE toCometary_Orb





  SUBROUTINE toKeplerian_Orb(this)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(inout)  :: this
    REAL(bp), DIMENSION(6)       :: elements

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Orbit / toKeplerian", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    elements = getElements(this, "keplerian")
    IF (error) THEN
       CALL errorMessage("Orbit / toKeplerian", &
            "TRACE BACK (5).", 1)
       RETURN
    END IF
    this%elements = elements
    this%element_type = "keplerian"
    this%frame = "ecliptic"

  END SUBROUTINE toKeplerian_Orb





END MODULE Orbit_cl
