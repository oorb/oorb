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
!! Contains integrators and force routines.
!!
!! @author  TL, MG, JV
!! @version 2009-06-12
!!
MODULE integrators

  USE planetary_data
  USE linal

  IMPLICIT NONE
  PRIVATE
  !! High precicion real type for integration variables.
  INTEGER, PARAMETER :: prec = rprec8
  !! Value of pi (32 digits).
  REAL(prec), PARAMETER :: pi = 3.14159265358979323846264338327950_prec
  !! Gaussian gravitational constant.
  REAL(prec), PARAMETER :: ggc = 0.01720209895_prec
  !! Constant of gravitation.
  REAL(prec), PARAMETER :: gc = ggc**2
  !! Speed of light (AU/day).
  REAL(prec), PARAMETER :: c = 173.14463272_prec
  REAL(prec), PARAMETER :: ic2 = 1.0_prec/c**2
  !! Decimal precision. 
  REAL(prec), PARAMETER :: rstep_tol = 10*EPSILON(pi)
  !! Decimal precision for extrapolation (or step convergence). 
  REAL(prec), PARAMETER :: bs_extrapol_tol = 10*EPSILON(pi)
  !! Default step length (1 days):
  REAL(prec), PARAMETER :: step_def = 1.0_prec
  !! Sequence of substeps (Bulirsch and Stoer).
  INTEGER, PARAMETER :: nseq = 27
  INTEGER, DIMENSION(nseq), PARAMETER :: &
       seq = (/ 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, &
       256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6144, 8192, &
       12288, 16384, 24576 /) ! n_i = 2 * n_(i-2)

  LOGICAL, PARAMETER :: relativity = .TRUE.
  ! Default central body for the dynamical system
  INTEGER, PARAMETER :: central_body_prm = 11
  INTEGER :: central_body

  ! NOTE: central_body should be an input argument in all
  ! the public routines, since from thereon it is currently
  ! passed on as an global module parameter.
  PUBLIC :: bulirsch_full_jpl

  INTERFACE ratf_extrapolation
     MODULE PROCEDURE ratf_extrapolation_vec, ratf_extrapolation_mat, &
          ratf_extrapolation_vec_n, ratf_extrapolation_mat_n 
  END INTERFACE

  INTERFACE polf_extrapolation
     MODULE PROCEDURE polf_extrapolation_vec, polf_extrapolation_mat, &
          polf_extrapolation_vec_n, polf_extrapolation_mat_n
  END INTERFACE

CONTAINS





  !! Description: 
  !!   Bulirsch-Stoer scheme.
  !!   Rational function extrapolation w/ modified midpoint. 
  !!
  !! References:
  !!   [1] Press et al. 1989, Numerical Recipes. 
  !!   [2] Stoer & Bulirsch 1980, Introduction to Num. Anal.
  !!
  !! Usage:
  !!   CALL bulirsch_full_jpl(mjd_tdt0, mjd_tdt1, &
  !!                          celements, error, jacobian, step)
  !!
  !! mjd_tdt0         integration start (MJD TT)
  !! mjd_tdt1         integration stop (MJD TT)
  !! celements        initial coordinates for massless particles
  !! error            returns true, if something fails 
  !! jacobian         jacobian matrix (coordinates wrt initial coordinates)
  !! step             step size
  !! ncenter          number of solar-system object to use as center (default=Sun)
  !! encounters       table containing the closest distances to, or earliest
  !!                  time of impact with solar-system objects
  !! addit_mcelements initial coordinates for additional perturbing bodies to be
  !!                  integrated simultaneously
  !! addit_masses     masses for additional perturbing bodies
  !!
  SUBROUTINE bulirsch_full_jpl(mjd_tdt0, mjd_tdt1, celements, &
       perturbers, error, jacobian, step, ncenter, encounters, &
       addit_celements, addit_masses)

    REAL(prec), INTENT(in)                                :: mjd_tdt0, mjd_tdt1
    REAL(prec), DIMENSION(:,:), INTENT(inout)             :: celements
    LOGICAL, DIMENSION(:), INTENT(in)                     :: perturbers
    LOGICAL, INTENT(inout)                                :: error
    REAL(prec), DIMENSION(:,:,:), INTENT(inout), OPTIONAL :: jacobian
    REAL(prec), INTENT(in), OPTIONAL                      :: step
    INTEGER,  INTENT(in), OPTIONAL                        :: ncenter
    REAL(prec), DIMENSION(:,:,:), INTENT(out), OPTIONAL   :: encounters
    REAL(prec), DIMENSION(:,:), INTENT(inout), OPTIONAL   :: addit_celements
    REAL(prec), DIMENSION(:), INTENT(in), OPTIONAL        :: addit_masses

    REAL(prec), DIMENSION(:,:,:), ALLOCATABLE :: pws, encounters_
    REAL(prec), DIMENSION(:,:), ALLOCATABLE :: ws
    REAL(prec) :: mjd_tdt, tmp, istep, rstep
    INTEGER    :: k, m, n, total, err

    IF (PRESENT(addit_celements) .AND. .NOT.PRESENT(addit_masses) .OR. &
         .NOT.PRESENT(addit_celements) .AND. PRESENT(addit_masses)) THEN
       error = .TRUE.
       WRITE(0,"(A)") "bulirsch_full_jpl: Both addit_celements and addit_masses must be specified."
       RETURN       
    END IF

    ALLOCATE(ws(SIZE(celements,dim=1),SIZE(celements,dim=2)), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       WRITE(0,"(A)") "bulirsch_full_jpl: Could not allocate memory (5)."
       RETURN
    END IF

    IF (PRESENT(jacobian)) THEN
       ALLOCATE(pws(SIZE(jacobian,dim=1),SIZE(jacobian,dim=2), &
            SIZE(jacobian,dim=3)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "bulirsch_full_jpl: Could not allocate memory (10)."
          DEALLOCATE(ws, pws, stat=err)
          RETURN
       END IF
    END IF

    ! Reset the global module parameter
    IF (PRESENT(ncenter)) THEN
       central_body = ncenter
    ELSE 
       central_body = central_body_prm
    END IF

    ! Initialize the encounter table
    IF (PRESENT(encounters)) THEN
       IF (SIZE(encounters,dim=1) < SIZE(celements,dim=2)) THEN
          error = .TRUE.
          WRITE(0,"(A)") "bulirsch_full_jpl: 'encounters' array too small (1)."
          RETURN
       END IF
       IF (SIZE(encounters,dim=2) < 11) THEN
          error = .TRUE.
          WRITE(0,"(A)") "bulirsch_full_jpl: 'encounters' array too small (2)."
          RETURN
       END IF
       IF (SIZE(encounters,dim=3) < 4) THEN
          error = .TRUE.
          WRITE(0,"(A)") "bulirsch_full_jpl: 'encounters' array too small (3)."
          RETURN
       END IF
       encounters = HUGE(encounters)
       encounters(:,:,2) = 3
       ALLOCATE(encounters_(SIZE(encounters,dim=1),SIZE(encounters,dim=2),SIZE(encounters,dim=3)))
    END IF

    ! Set integration parameters.
    tmp = mjd_tdt1 - mjd_tdt0
    IF (PRESENT(step)) THEN
       istep = SIGN(step,tmp)
    ELSE
       ! Use default step length:
       istep = SIGN(step_def,tmp)
    END IF
    rstep = SIGN(MOD(tmp,istep),tmp)
    total = ABS(INT(tmp/istep))

    ! Integration loop
    k       = 1
    mjd_tdt = mjd_tdt0
    ws      = celements
    IF (PRESENT(jacobian)) THEN
       pws = jacobian
       DO WHILE (k <= total)
          IF (PRESENT(encounters)) THEN
             CALL step_bulirsch_full_jpl(mjd_tdt, istep, perturbers, ws, &
                  ws, error, pws, pws, encounters=encounters_)
             ! Log closest non-impacting encounter during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
             ! Log earliest time of impact during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) < 1.1_prec .AND. encounters_(m,n,2) < 1.1_prec .AND. &
                  encounters_(m,n,1) < encounters(m,n,1))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
          ELSE
             CALL step_bulirsch_full_jpl(mjd_tdt, istep, perturbers, ws, &
                  ws, error, pws, pws)
          END IF
          IF (error) THEN
             DEALLOCATE(ws, pws, stat=err)
             WRITE(0,"(A)") "bulirsch_full_jpl: TRACE BACK (5)"
             RETURN
          END IF
          mjd_tdt = mjd_tdt0 + k * istep
          k       = k + 1
       END DO
       IF (ABS(rstep) > rstep_tol) THEN
          IF (PRESENT(encounters)) THEN
             CALL step_bulirsch_full_jpl(mjd_tdt, rstep, perturbers, ws, &
                  ws, error, pws, pws, encounters=encounters_)
          ELSE
             CALL step_bulirsch_full_jpl(mjd_tdt, rstep, perturbers, ws, &
                  ws, error, pws, pws)
          END IF
       ELSE
          IF (PRESENT(encounters)) THEN
             CALL step_midpoint_full_jpl(mjd_tdt, rstep, 10, perturbers, &
                  ws, ws, error, pws, pws, encounters=encounters_)
          ELSE
             CALL step_midpoint_full_jpl(mjd_tdt, rstep, 10, perturbers, &
                  ws, ws, error, pws, pws)             
          END IF
       END IF
       IF (PRESENT(encounters)) THEN
          ! Log closest non-impacting encounter during the integration step
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
             encounters(m,n,:) = encounters_(m,n,:)
          END FORALL
          ! Log earliest time of impact during the integration step
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters(m,n,2) < 1.1_prec .AND. encounters_(m,n,2) < 1.1_prec .AND. &
               encounters_(m,n,1) < encounters(m,n,1))
             encounters(m,n,:) = encounters_(m,n,:)
          END FORALL
       END IF
       celements = ws
       jacobian  = pws
       DEALLOCATE(pws, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "bulirsch_full_jpl: Could not deallocate memory (5)."
          RETURN
       END IF
    ELSE
       DO WHILE (k <= total)
          IF (PRESENT(encounters)) THEN
             CALL step_bulirsch_full_jpl(mjd_tdt, istep, perturbers, ws, &
                  ws, error, encounters=encounters_)
             ! Log closest non-impacting encounter during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
             ! Log earliest time of impact during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) < 1.1_prec .AND. encounters_(m,n,2) < 1.1_prec .AND. &
                  encounters_(m,n,1) < encounters(m,n,1))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
          ELSE
             CALL step_bulirsch_full_jpl(mjd_tdt, istep, perturbers, ws, &
                  ws, error)
          END IF
          IF (error) THEN
             DEALLOCATE(ws, stat=err)
             WRITE(0,"(A)") "bulirsch_full_jpl: TRACE BACK (10)"
             RETURN
          END IF
          mjd_tdt = mjd_tdt0 + k * istep
          k       = k + 1
       END DO
       IF (ABS(rstep) > rstep_tol) THEN
          IF (PRESENT(encounters)) THEN
             CALL step_bulirsch_full_jpl(mjd_tdt, rstep, perturbers, ws, &
                  ws, error, encounters=encounters_)
          ELSE
             CALL step_bulirsch_full_jpl(mjd_tdt, rstep, perturbers, ws, &
                  ws, error)
          END IF
       ELSE
          IF (PRESENT(encounters)) THEN
             CALL step_midpoint_full_jpl(mjd_tdt, rstep, 10, perturbers, &
                  ws, ws, error, encounters=encounters_)
          ELSE
             CALL step_midpoint_full_jpl(mjd_tdt, rstep, 10, perturbers, &
                  ws, ws, error)
          END IF
       END IF
       IF (PRESENT(encounters)) THEN
          ! Log closest non-impacting encounter during the integration step
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
             encounters(m,n,:) = encounters_(m,n,:)
          END FORALL
          ! Log earliest time of impact during the integration step
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters(m,n,2) < 1.1_prec .AND. encounters_(m,n,2) < 1.1_prec .AND. &
               encounters_(m,n,1) < encounters(m,n,1))
             encounters(m,n,:) = encounters_(m,n,:)
          END FORALL
       END IF
       celements = ws
    END IF
    DEALLOCATE(ws, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       WRITE(0,"(A)") "bulirsch_full_jpl: Could not deallocate memory (10)."
       RETURN
    END IF
    IF (PRESENT(encounters)) THEN
       DEALLOCATE(encounters_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "bulirsch_full_jpl: Could not deallocate memory (15)."
          RETURN
       END IF
    END IF

  END SUBROUTINE bulirsch_full_jpl





  !! One integration step.
  !!
  !!  mjd_tdt  modified Julian date (for ephemerides)
  !!  H        step size
  !!  ws0      initial coordinates for the massless bodies
  !!  ws1      final coordinates for the massless bodies
  !!  error    returns true, if something fails
  !!  pws0     initial values of the partial derivatives
  !!  pws1     final values of the partial derivatives
  !!  
  SUBROUTINE step_bulirsch_full_jpl(mjd_tdt, H, perturbers, ws0, ws1, &
       error, pws0, pws1, encounters)

    REAL(prec), INTENT(in)                                :: mjd_tdt, H
    LOGICAL, DIMENSION(:), INTENT(in)                     :: perturbers
    REAL(prec), DIMENSION(:,:), INTENT(inout)             :: ws0, ws1
    LOGICAL, INTENT(inout)                                :: error
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(inout) :: pws0, pws1
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(out)   :: encounters

    REAL(prec), DIMENSION(:,:,:,:,:), ALLOCATABLE :: pddif
    REAL(prec), DIMENSION(:,:,:,:), ALLOCATABLE :: pwa0, pwa1, ddif
    REAL(prec), DIMENSION(:,:,:), ALLOCATABLE :: wa0, wa1, pwst, encounters_
    REAL(prec), DIMENSION(:,:), ALLOCATABLE :: wst
    REAL(prec), DIMENSION(:), ALLOCATABLE :: hseq
    INTEGER, DIMENSION(:), ALLOCATABLE :: ws_index, pws_index
    INTEGER :: ws_final, pws_final, i, k, l, m, n, NS, nseq, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: ws_converged, pws_converged

    ALLOCATE(hseq(SIZE(seq)), &
         wa0(6,SIZE(seq),SIZE(ws0,dim=2)), &
         wa1(6,SIZE(seq),SIZE(ws0,dim=2)), &
         ddif(6,SIZE(seq),SIZE(seq),SIZE(ws0,dim=2)), &
         wst(6,SIZE(ws0,dim=2)), &
         ws_converged(SIZE(ws0,dim=2)), &
         ws_index(SIZE(ws0,dim=2)), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       WRITE(0,"(A)") "step_bulirsch_full_jpl: Could not allocate memory (5)."
       DEALLOCATE(hseq, wa0, wa1, pwa0, pwa1, ddif, pddif, wst, &
            pwst, ws_converged, pws_converged, ws_index, pws_index, &
            stat=err)
       RETURN
    END IF
    ws_index = 0
    ws_converged = .FALSE.

    IF (PRESENT(encounters)) THEN
       ALLOCATE(encounters_(SIZE(encounters,dim=1),SIZE(encounters,dim=2),SIZE(encounters,dim=3)))
       encounters = HUGE(encounters)
       encounters(:,:,2) = 3
    END IF

    IF (PRESENT(pws0) .AND. PRESENT(pws1)) THEN
       ALLOCATE(pwa0(6,SIZE(seq),6,SIZE(ws0,dim=2)), &
            pwa1(6,SIZE(seq),6,SIZE(ws0,dim=2)), &
            pddif(6,SIZE(seq),SIZE(seq),6,SIZE(ws0,dim=2)), &
            pwst(6,6,SIZE(ws0,dim=2)), &
            pws_converged(SIZE(ws0,dim=2)), &
            pws_index(SIZE(ws0,dim=2)), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "step_bulirsch_full_jpl: Could not allocate memory (10)."
          DEALLOCATE(hseq, wa0, wa1, pwa0, pwa1, ddif, pddif, wst, &
               pwst, ws_converged, pws_converged, ws_index, &
               pws_index, stat=err)
          RETURN
       END IF
       pws_index = 0
       pws_converged = .FALSE.
    END IF

    ws_final = 0
    pws_final = 0
    nseq = SIZE(seq)
    NS = SIZE(ws0,dim=2)

    DO i=1, nseq

       hseq(i) = (H / seq(i)) ** 2

       ! Integration with seq(i) substeps.
       IF (PRESENT(pws0) .AND. PRESENT(pws1)) THEN
          IF (PRESENT(encounters)) THEN
             CALL step_midpoint_full_jpl(mjd_tdt, H, seq(i), perturbers, &
                  ws0, wst, error, pws0, pwst, encounters=encounters_)
             ! Log closest non-impacting encounter during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
             ! Log earliest time of impact during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) < 1.1_prec .AND. encounters_(m,n,2) < 1.1_prec .AND. &
                  encounters_(m,n,1) < encounters(m,n,1))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
          ELSE
             CALL step_midpoint_full_jpl(mjd_tdt, H, seq(i), perturbers, &
                  ws0, wst, error, pws0, pwst)
          END IF
          DO k=1, NS
             wa0(:,i,k) = wst(:,k)
             DO l=1, 6
                pwa0(:,i,l,k) = pwst(:,l,k)
             END DO
          END DO
       ELSE
          IF (PRESENT(encounters)) THEN
             CALL step_midpoint_full_jpl(mjd_tdt, H, seq(i), perturbers, &
                  ws0, wst, error, encounters=encounters_)
             ! Log closest non-impacting encounter during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
             ! Log earliest time of impact during the integration step
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) < 1.1_prec .AND. encounters_(m,n,2) < 1.1_prec .AND. &
                  encounters_(m,n,1) < encounters(m,n,1))
                encounters(m,n,:) = encounters_(m,n,:)
             END FORALL
          ELSE
             CALL step_midpoint_full_jpl(mjd_tdt, H, seq(i), perturbers, &
                  ws0, wst, error)
          END IF
          DO k=1, NS
             wa0(:,i,k) = wst(:,k)
          END DO
       END IF
       IF (error) THEN
          WRITE(0,"(A)") "step_bulirsch_full_jpl: TRACE BACK (5)."
          DEALLOCATE(hseq, wa0, wa1, pwa0, pwa1, ddif, pddif, wst, &
               pwst, ws_converged, pws_converged, ws_index, pws_index, &
               stat=err)
          RETURN
       END IF

       ! Fill in the extrapolation table and check for convergence.
       IF (ws_final == 0) THEN
          CALL polf_extrapolation(i, hseq, wa0, wa1, ddif, ws_converged, error)
          !          CALL ratf_extrapolation(i, hseq, wa0, wa1, ddif, ws_converged, error)
          IF (error) THEN
             WRITE(0,"(A)") "step_bulirsch_full_jpl: TRACE BACK (10)."
             DEALLOCATE(hseq, wa0, wa1, pwa0, pwa1, ddif, pddif, wst, &
                  pwst, ws_converged, pws_converged, ws_index, pws_index, &
                  stat=err)
             RETURN
          END IF
          WHERE (ws_converged .AND. ws_index == 0) 
             ws_index = i
          END WHERE
          IF (ALL(ws_index > 0)) THEN
             ws_final = i
          END IF
       END IF
       IF (PRESENT(pws0) .AND. PRESENT(pws1)) THEN
          IF (pws_final == 0) THEN
             CALL polf_extrapolation(i, hseq, pwa0, pwa1, pddif, pws_converged, error)
             !             CALL ratf_extrapolation(i, hseq, pwa0, pwa1, pddif, pws_converged, error)
             IF (error) THEN
                WRITE(0,"(A)") "step_bulirsch_full_jpl: TRACE BACK (15)."
                DEALLOCATE(hseq, wa0, wa1, pwa0, pwa1, ddif, pddif, wst, &
                     pwst, ws_converged, pws_converged, ws_index, pws_index, &
                     stat=err)
                RETURN
             END IF
             WHERE (pws_converged .AND. pws_index == 0)
                pws_index = i
             END WHERE
             IF (ALL(pws_index > 0)) THEN
                pws_final = i
             END IF
          END IF
       ELSE
          pws_final = 1
       END IF
       IF (ws_final /= 0 .AND. pws_final /= 0) THEN
          EXIT
       END IF

    END DO

    ! If there were problems with convergence.
    IF (ws_final == 0) THEN
       ws_final = nseq
       DO k=1, NS
          IF (ws_index(k) == 0) THEN
             ! Error report
             ws_index(k) = nseq
          END IF
       END DO
    END IF
    IF (pws_final == 0) THEN
       pws_final = nseq
       DO k=1, NS
          IF (pws_index(k) == 0) THEN
             ! Error report
             pws_index(k) = nseq
          END IF
       END DO
    END IF

    IF (PRESENT(pws0) .AND. PRESENT(pws1)) THEN
       DO k=1, NS
          ws1(:,k) = wa1(:,ws_index(k),k)
          DO l=1, 6
             pws1(:,l,k) = pwa1(:,pws_index(k),l,k)
          END DO
       END DO
    ELSE
       DO k=1, NS
          ws1(:,k) = wa1(:,ws_index(k),k)
       END DO
    END IF

    DEALLOCATE(hseq, wa0, wa1, ddif, wst, ws_converged, ws_index, &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       WRITE(0,"(A)") "step_bulirsch_full_jpl: Could not deallocate memory (5)."
       RETURN
    END IF
    IF (PRESENT(pws0) .AND. PRESENT(pws1)) THEN
       DEALLOCATE(pwa0, pwa1, pddif, pwst, pws_converged, pws_index, &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "step_bulirsch_full_jpl: Could not deallocate memory (10)."
          RETURN
       END IF
    END IF
    IF (PRESENT(encounters)) THEN
       DEALLOCATE(encounters_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "step_bulirsch_full_jpl: Could not deallocate memory (15)."
          RETURN
       END IF
    END IF

  END SUBROUTINE step_bulirsch_full_jpl






  !! Description: 
  !!   Modified midpoint step (2nd order).
  !! 
  !! References: 
  !!   [1] Press et al. 1989, Numerical Recipes.
  !!
  !! Usage:
  !!   CALL step_midpoint_full_jpl(mjd_tdt, h, nsteps, ws0, ws1, 
  !!                               error, pws0, pws1)
  !!
  !! One integration step.
  !!
  !!  mjd_tdt  modified Julian date (for ephemerides)
  !!  h        (big) step size
  !!  nsteps   number of substeps
  !!  ws0      initial coordinates for the massless bodies
  !!  ws1      final coordinates for the massless bodies
  !!  error    returns true, if something fails                        
  !!  pws0     initial values of the partial derivatives
  !!  pws1     final values of the partial derivatives
  !!
  SUBROUTINE step_midpoint_full_jpl(mjd_tdt, h, nsteps, perturbers, &
       ws0, ws1, error, pws0, pws1, encounters)

    REAL(prec), INTENT(in)                                :: mjd_tdt, h
    INTEGER, INTENT(in)                                   :: nsteps
    LOGICAL, DIMENSION(:), INTENT(in)                     :: perturbers
    REAL(prec), DIMENSION(:,:), INTENT(inout)             :: ws0, ws1
    LOGICAL, INTENT(inout)                                :: error
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(inout) :: pws0, pws1
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(out)   :: encounters

    INTEGER               :: err, i, k, m, n, NS
    INTEGER, DIMENSION(2) :: sw
    REAL(prec)            :: dt

    REAL(prec), DIMENSION(6,SIZE(ws0,dim=2),2)   :: q
    REAL(prec), DIMENSION(6,SIZE(ws0,dim=2))     :: qd
    REAL(prec), DIMENSION(6,6,SIZE(ws0,dim=2),2) :: pq
    REAL(prec), DIMENSION(6,6,SIZE(ws0,dim=2))   :: pqd
    REAL(prec), DIMENSION(:,:,:), ALLOCATABLE    :: encounters_

    REAL(prec), DIMENSION(6,6) :: pqt

    NS = SIZE(ws0,dim=2)

    ! Substep size.
    dt = h / nsteps

    ! Algorithm needs information from two previous points.
    ! When calculating at k+1, sw(1) points to variables at k-1,
    ! and sw(2) to variables at k.   
    sw(1)         = 1
    sw(2)         = 2
    q(:,:,sw(1))  = ws0
    IF (PRESENT(encounters)) THEN
       ALLOCATE(encounters_(SIZE(encounters,dim=1), &
            SIZE(encounters,dim=2), SIZE(encounters,dim=3)), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "step_midpoint_full_jpl: Could not allocate memory (5)."
          RETURN
       END IF
       encounters = 0.0_prec
       encounters(:,:,2) = 3
    END IF

    ! If partial derivatives are needed.
    IF (PRESENT(pws0) .AND. PRESENT(pws1)) THEN
       pq(:,:,:,sw(1)) = pws0
       IF (PRESENT(encounters)) THEN
          CALL interact_full_jpl(q(:,:,sw(1)), mjd_tdt, perturbers, &
               COUNT(perturbers), qd, error, pqd, encounters=encounters_)
          ! Initialize log
          encounters = encounters_
          encounters(:,:,4) = dt
       ELSE
          CALL interact_full_jpl(q(:,:,sw(1)), mjd_tdt, perturbers, &
               COUNT(perturbers), qd, error, pqd)
       END IF
       IF (error) THEN
          RETURN
       END IF
       q(:,:,sw(2)) = q(:,:,sw(1)) + dt * qd
       DO i=1, NS
          pqt = MATMUL(pqd(:,:,i),pq(:,:,i,sw(1)))
          pq(:,:,i,sw(2)) = pq(:,:,i,sw(1)) + dt * pqt
       END DO
       DO k=2, nsteps     
          IF (PRESENT(encounters)) THEN
             CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + (k - 1) * dt, &
                  perturbers, COUNT(perturbers), qd, error, pqd, encounters=encounters_)
             ! Log closest non-impacting encounter during the integration substep
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
                encounters(m,n,1:3) = encounters_(m,n,1:3)
             END FORALL
             ! Log earliest time of impact during the integration substep
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters_(m,n,2) < 1.1_prec .AND. encounters_(m,n,1) < encounters(m,n,1))
                encounters(m,n,1:3) = encounters_(m,n,1:3)
             END FORALL
          ELSE
             CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + (k - 1) * dt, &
                  perturbers, COUNT(perturbers), qd, error, pqd)
          END IF
          IF (error) THEN
             RETURN
          END IF
          q(:,:,sw(1)) = q(:,:,sw(1)) + 2.0_prec * dt * qd
          DO i = 1, NS
             pqt = MATMUL(pqd(:,:,i),pq(:,:,i,sw(2)))
             pq(:,:,i,sw(1)) = pq(:,:,i,sw(1)) + 2.0_prec * dt * pqt
          END DO
          sw(1:2) = sw(2:1:-1)
       END DO
       IF (PRESENT(encounters)) THEN
          CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + h, perturbers, &
               COUNT(perturbers), qd, error, pqd, encounters=encounters_)
          ! Log closest non-impacting encounter during the integration substep
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
             encounters(m,n,1:3) = encounters_(m,n,1:3)
          END FORALL
          ! Log earliest time of impact during the integration substep
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters_(m,n,2) < 1.1_prec .AND. encounters_(m,n,1) < encounters(m,n,1))
             encounters(m,n,1:3) = encounters_(m,n,1:3)
          END FORALL
       ELSE
          CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + h, perturbers, &
               COUNT(perturbers), qd, error, pqd)
       END IF
       IF (error) THEN
          RETURN
       END IF
       q(:,:,sw(1)) = q(:,:,sw(1)) + dt * qd
       DO i=1, NS
          pqt = MATMUL(pqd(:,:,i),pq(:,:,i,sw(2)))
          pq(:,:,i,sw(1)) = pq(:,:,i,sw(1)) + dt * pqt
       END DO
       ws1  = (q(:,:,sw(1)) + q(:,:,sw(2))) / 2.0_prec
       pws1 = (pq(:,:,:,sw(1)) + pq(:,:,:,sw(2))) / 2.0_prec
    ELSE
       ! Plain integration.
       IF (PRESENT(encounters)) THEN
          CALL interact_full_jpl(q(:,:,sw(1)), mjd_tdt, perturbers, &
               COUNT(perturbers), qd, error, encounters=encounters_)
          ! Initialize log
          encounters(:,:,1:3) = encounters_(:,:,1:3)
          encounters(:,:,4) = dt
       ELSE
          CALL interact_full_jpl(q(:,:,sw(1)), mjd_tdt, perturbers, &
               COUNT(perturbers), qd, error)
       END IF
       IF (error) THEN
          RETURN
       END IF
       q(:,:,sw(2)) = q(:,:,sw(1)) + dt * qd
       DO k=2, nsteps     
          IF (PRESENT(encounters)) THEN
             CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + (k - 1) * dt, &
                  perturbers, COUNT(perturbers), qd, error, encounters=encounters_)
             ! Log closest non-impacting encounter during the integration substep
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
                encounters(m,n,1:3) = encounters_(m,n,1:3)
             END FORALL
             ! Log earliest time of impact during the integration substep
             FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
                  encounters_(m,n,2) < 1.1_prec .AND. encounters_(m,n,1) < encounters(m,n,1))
                encounters(m,n,1:3) = encounters_(m,n,1:3)
             END FORALL
          ELSE
             CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + (k - 1) * dt, &
                  perturbers, COUNT(perturbers), qd, error)
          END IF
          IF (error) THEN
             RETURN
          END IF
          q(:,:,sw(1)) = q(:,:,sw(1)) + 2.0_prec * dt * qd
          sw(1:2)      = sw(2:1:-1)
       END DO
       IF (PRESENT(encounters)) THEN
          CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + h, perturbers, &
               COUNT(perturbers), qd, error, encounters=encounters_)
          ! Log closest non-impacting encounter during the integration substep
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters(m,n,2) > 1.1_prec .AND. encounters_(m,n,3) < encounters(m,n,3))
             encounters(m,n,1:3) = encounters_(m,n,1:3)
          END FORALL
          ! Log earliest time of impact during the integration substep
          FORALL (m=1:SIZE(encounters,dim=1), n=1:SIZE(encounters,dim=2), &
               encounters_(m,n,2) < 1.1_prec .AND. encounters_(m,n,1) < encounters(m,n,1))
             encounters(m,n,1:3) = encounters_(m,n,1:3)
          END FORALL
       ELSE
          CALL interact_full_jpl(q(:,:,sw(2)), mjd_tdt + h, perturbers, &
               COUNT(perturbers), qd, error)
       END IF
       IF (error) THEN
          RETURN
       END IF
       q(:,:,sw(1)) = q(:,:,sw(1)) + dt * qd
       ws1 = (q(:,:,sw(1)) + q(:,:,sw(2))) / 2.0_prec
    END IF

    IF (PRESENT(encounters)) THEN
       DEALLOCATE(encounters_, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          WRITE(0,"(A)") "step_midpoint_full_jpl: Could not deallocate memory (5)."
          RETURN
       END IF
    END IF


  END SUBROUTINE step_midpoint_full_jpl





  !! *Description*: 
  !!
  !!   Rational and polynomial function extrapolation for
  !!   (Gragg-)Bulirsch-Stoer integration.
  !!
  !!   These routines are applicable to single vectors or single matrices
  !!   or to groups of vectors or matrices. 
  !!   Extrapolation is always towards zero. 
  !!   Extrapolation converges when all of the elements of a vector/matrix
  !!   are converged. For groups of N vectors/matrices, convergence is checked
  !!   individually N times.
  !!
  !! References:
  !!   [1] Press et al. 1989, Numerical Recipes. 
  !!   [2] Stoer & Bulirsch 1980, Introduction to Num. Anal.
  !!
  !! Usage:
  !!   CALL ratf_extrapolation(z, h, w0, w1, ddif, converged)
  !!   CALL polf_extrapolation(z, h, w0, w1, ddif, converged)
  !!
  !! Fill one diagonal in rational function extrapolation table.
  !! Applicable to vectors.
  !!
  !!  z          index to process in extrapolation table  
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h 
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences (see below) 
  !!  converged  flag for convergence
  !!
  SUBROUTINE ratf_extrapolation_vec(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                         :: z
    REAL(prec), DIMENSION(:), INTENT(in)        :: h
    REAL(prec), DIMENSION(:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:), INTENT(inout) :: ddif
    LOGICAL, INTENT(out)                        :: converged
    LOGICAL, INTENT(inout)                      :: error

    REAL(prec), DIMENSION(:), ALLOCATABLE :: cdif, tmp1, tmp2
    REAL(prec) :: delta
    INTEGER :: i, j, k, err

    ! Fill in the extrapolation table: (see [1] pp. 80-86 & 563-568)
    ! Here, tables are indexed as (i,j), where the rows i represent
    ! integrations with seq(i) substeps and columns j extrapolations
    ! with j extrapolation points (cf. [1]: j=m+1 in (3.2.8)).
    ! Table ddif is filled like this:
    !   w0(1)   ddif(1,2) ddif(1,3)     ...
    !   w0(2)   ddif(2,2)         0     ...
    !   w0(3)           0         0     ...
    !     ...
    ! Extrapolated value to stepsize=0 is w0(3)+ddif(2,2)+ddif(1,3).

    converged = .FALSE.
    w1(:,z) = w0(:,z)

    IF (z == 1) THEN
       RETURN
    END IF

    ALLOCATE(cdif(SIZE(w0,dim=1)), tmp1(SIZE(w0,dim=1)), &
         tmp2(SIZE(w0,dim=1)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(cdif, tmp1, tmp2, stat=err)
       RETURN
    END IF

    ddif(:,z-1,1) = w0(:,z-1)
    cdif          = w0(:,z)

    DO i=z-1, 1, -1
       j    = z - i + 1
       tmp1 = h(i) / h(z) * ddif(:,i,j-1)
       tmp2 = (tmp1 - cdif)
       DO k=1, 6
          ! Check for division by zero.
          IF (tmp2(k) /= 0) THEN
             tmp2(k)     = (cdif(k) - ddif(k,i,j-1)) / tmp2(k)
             ddif(k,i,j) = cdif(k) * tmp2(k)
             cdif(k)     = tmp1(k) * tmp2(k)
             ! Add difference towards final extrapolation.
             w1(k,z)     = w1(k,z) + ddif(k,i,j)
          ELSE
             ddif(k,i,j) = 0.0_prec
             cdif(k)     = 0.0_prec                
          END IF
       END DO
    END DO

    ! Check if extrapolation is converged.
    delta = MAXVAL(ABS(ddif(:,1,z)))
    IF (delta < bs_extrapol_tol) THEN
       converged = .TRUE.
    END IF

    DEALLOCATE(cdif, tmp1, tmp2, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE ratf_extrapolation_vec





  !! Fill one diagonal in polynomial function extrapolation table.
  !! Applicable to vectors.
  !!
  !!  z          index to process in extrapolation table  
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h 
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences 
  !!  converged  flag for convergence
  !!
  SUBROUTINE polf_extrapolation_vec(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                         :: z
    REAL(prec), DIMENSION(:), INTENT(in)        :: h
    REAL(prec), DIMENSION(:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:), INTENT(inout) :: ddif
    LOGICAL, INTENT(out)                        :: converged
    LOGICAL, INTENT(inout)                      :: error

    REAL(prec)                            :: delta
    REAL(prec), DIMENSION(:), ALLOCATABLE :: cdif, tmp
    INTEGER                               :: i, j, err

    converged = .FALSE.
    w1(:,z) = w0(:,z)

    IF (z == 1) THEN
       RETURN
    END IF

    ALLOCATE(cdif(SIZE(w0,dim=1)), tmp(SIZE(w0,dim=1)), &
         stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(cdif, tmp, stat=err)
       RETURN
    END IF

    ddif(:,z-1,1) = w0(:,z-1)
    cdif          = w0(:,z)

    DO i=z-1, 1, -1
       j           = z - i + 1
       tmp         = (cdif - ddif(:,i,j-1)) / (h(i) - h(z))
       ddif(:,i,j) = h(z) * tmp
       cdif        = h(i) * tmp
       w1(:,z)     = w1(:,z) + ddif(:,i,j)
    END DO

    ! Check if extrapolation is converged.
    delta = MAXVAL(ABS(ddif(:,1,z)))
    IF (delta < bs_extrapol_tol) THEN
       converged = .TRUE.
    END IF

    DEALLOCATE(cdif, tmp, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE polf_extrapolation_vec





  !! Fill one diagonal in rational function extrapolation table.
  !! Applicable to matrices.
  !!
  !!  z          index to process in extrapolation table  
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h 
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences (see below)
  !!  converged  flag for convergence
  !!
  SUBROUTINE ratf_extrapolation_mat(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                           :: z
    REAL(prec), DIMENSION(:), INTENT(in)          :: h
    REAL(prec), DIMENSION(:,:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:,:), INTENT(inout) :: ddif
    LOGICAL, INTENT(out)                          :: converged
    LOGICAL, INTENT(inout)                        :: error

    INTEGER                            :: i, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: column_ok

    ! This routine calls ratf_extrapolation_vec (for vectors) for each
    ! column of the extrapolated matrix. Therefore, the last index of 
    ! w0, w1, and ddif is the column index of the matrix. Note this 
    ! when accessing the corresponding actual arguments.

    converged = .FALSE.

    ALLOCATE(column_ok(SIZE(w0,dim=3)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(column_ok, stat=err)
       RETURN
    END IF

    DO i=1, SIZE(w0,dim=3)
       CALL ratf_extrapolation_vec(z, h, w0(:,:,i), w1(:,:,i), &
            ddif(:,:,:,i), column_ok(i), error)
       IF (error) THEN
          DEALLOCATE(column_ok, stat=err)
          RETURN
       END IF
    END DO
    IF (ALL(column_ok)) THEN
       converged = .TRUE.
    END IF

    DEALLOCATE(column_ok, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE ratf_extrapolation_mat





  !! Fill one diagonal in polynomial function extrapolation table.
  !! Applicable to matrices.
  !!
  !!  z          index to process in extrapolation table
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences
  !!  converged  flag for convergence
  !! 
  SUBROUTINE polf_extrapolation_mat(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                           :: z
    REAL(prec), DIMENSION(:), INTENT(in)          :: h
    REAL(prec), DIMENSION(:,:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:,:), INTENT(inout) :: ddif
    LOGICAL, INTENT(out)                          :: converged
    LOGICAL, INTENT(inout)                        :: error

    INTEGER                            :: i, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: column_ok

    converged = .FALSE.

    ALLOCATE(column_ok(SIZE(w0,dim=3)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(column_ok, stat=err)
       RETURN
    END IF

    DO i=1, SIZE(w0,dim=3) 
       CALL polf_extrapolation_vec(z, h, w0(:,:,i), w1(:,:,i), & 
            ddif(:,:,:,i), column_ok(i), error)
       IF (error) THEN
          DEALLOCATE(column_ok, stat=err)
          RETURN
       END IF
    END DO
    IF (ALL(column_ok)) THEN 
       converged = .TRUE. 
    END IF

    DEALLOCATE(column_ok, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE polf_extrapolation_mat





  !! Fill one diagonal in rational function extrapolation table.
  !! Applicable to N individual vectors.
  !!
  !!  z          index to process in extrapolation table
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences (see below)
  !!  converged  flags for convergence
  !!
  SUBROUTINE ratf_extrapolation_vec_n(z, h, w0, w1, ddif, converged, error)

    INTEGER, INTENT(in)                           :: z 
    REAL(prec), DIMENSION(:), INTENT(in)          :: h
    REAL(prec), DIMENSION(:,:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:,:), INTENT(inout) :: ddif
    LOGICAL, DIMENSION(:), INTENT(inout)          :: converged
    LOGICAL, INTENT(inout)                        :: error

    ! The last index of w0, w1, and ddif is now N, the number of vectors.
    ! "converged"-flags are used to check which vectors need further
    ! extrapolation.

    INTEGER :: i

    DO i=1, SIZE(w0,dim=3)
       IF (converged(i)) CYCLE
       CALL ratf_extrapolation_vec(z, h, w0(:,:,i), w1(:,:,i), & 
            ddif(:,:,:,i), converged(i), error)
       IF (error) THEN
          RETURN
       END IF
    END DO

  END SUBROUTINE ratf_extrapolation_vec_n





  !! Fill one diagonal in polynomial function extrapolation table.
  !! Applicable to N individual vectors.
  !!
  !!  z          index to process in extrapolation table
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences (see below)
  !!  converged  flags for convergence
  !!
  SUBROUTINE polf_extrapolation_vec_n(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                           :: z 
    REAL(prec), DIMENSION(:), INTENT(in)          :: h
    REAL(prec), DIMENSION(:,:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:,:), INTENT(inout) :: ddif
    LOGICAL, DIMENSION(:), INTENT(inout)          :: converged
    LOGICAL, INTENT(inout)                        :: error

    INTEGER :: i

    DO i=1, SIZE(w0,dim=3)
       IF (converged(i)) CYCLE
       CALL polf_extrapolation_vec(z, h, w0(:,:,i), w1(:,:,i), & 
            ddif(:,:,:,i), converged(i), error)
       IF (error) THEN
          RETURN
       END IF
    END DO

  END SUBROUTINE polf_extrapolation_vec_n





  !! Fill one diagonal in rational function extrapolation table.
  !! Applicable to N individual matrices.
  !!
  !!  z          index to process in extrapolation table
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences (see below)
  !!  converged  flags for convergence
  !!
  SUBROUTINE ratf_extrapolation_mat_n(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                             :: z 
    REAL(prec), DIMENSION(:), INTENT(in)            :: h
    REAL(prec), DIMENSION(:,:,:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:,:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:,:,:), INTENT(inout) :: ddif
    LOGICAL, DIMENSION(:), INTENT(inout)            :: converged
    LOGICAL, INTENT(inout)                          :: error

    ! The last index of w0, w1, and ddif is N, the number of vectors.
    ! The index before the last one is the column index of the matrix.
    ! Only matrices that are not already converged are extrapolated further.

    INTEGER :: i, j, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: column_ok

    ALLOCATE(column_ok(SIZE(w0,dim=3)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(column_ok, stat=err)
       RETURN
    END IF

    DO i=1, SIZE(w0,dim=4)
       IF (converged(i)) CYCLE
       DO j=1, SIZE(w0,dim=3)
          CALL ratf_extrapolation_vec(z, h, w0(:,:,j,i), w1(:,:,j,i), & 
               ddif(:,:,:,j,i), column_ok(j), error)
          IF (error) THEN
             DEALLOCATE(column_ok, stat=err)
             RETURN
          END IF
       END DO
       IF (ALL(column_ok)) THEN 
          converged(i) = .TRUE. 
       END IF
    END DO

    DEALLOCATE(column_ok, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE ratf_extrapolation_mat_n





  !! Fill one diagonal in polynomial function extrapolation table.
  !! Applicable to N individual matrices.
  !!
  !!  z          index to process in extrapolation table
  !!  h          sequence of step sizes used
  !!  w0         integrations corresponding the step sizes h
  !!  w1         extrapolated values with z integrations
  !!  ddif       extrapolation table with differences (see below)
  !!  converged  flags for convergence
  !!
  SUBROUTINE polf_extrapolation_mat_n(z, h, w0, w1, ddif, converged, &
       error)

    INTEGER, INTENT(in)                             :: z 
    REAL(prec), DIMENSION(:), INTENT(in)            :: h
    REAL(prec), DIMENSION(:,:,:,:), INTENT(in)      :: w0
    REAL(prec), DIMENSION(:,:,:,:), INTENT(out)     :: w1
    REAL(prec), DIMENSION(:,:,:,:,:), INTENT(inout) :: ddif
    LOGICAL, DIMENSION(:), INTENT(inout)            :: converged
    LOGICAL, INTENT(inout)                          :: error

    INTEGER :: i, j, err
    LOGICAL, DIMENSION(:), ALLOCATABLE :: column_ok

    ALLOCATE(column_ok(SIZE(w0,dim=3)), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(column_ok, stat=err)
       RETURN
    END IF

    DO i=1, SIZE(w0,dim=4)
       IF (converged(i)) CYCLE
       DO j=1, SIZE(w0,dim=3)
          CALL polf_extrapolation_vec(z, h, w0(:,:,j,i), w1(:,:,j,i), & 
               ddif(:,:,:,j,i), column_ok(j), error)
          IF (error) THEN
             DEALLOCATE(column_ok, stat=err)
             RETURN
          END IF
       END DO
       IF (ALL(column_ok)) THEN 
          converged(i) = .TRUE. 
       END IF
    END DO

    DEALLOCATE(column_ok, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE polf_extrapolation_mat_n





  !! Description: 
  !!   Evaluation of the full Newtonian force function for several 
  !!   massless bodies. Positions of the massive bodies are read 
  !!   from JPL ephemerides. A relativistic term due to the Sun is included.
  !!   Optional argument triggers evaluation of the partial derivatives
  !!   of the force function wrt Cartesian coordinates.
  !!
  !! References:
  !!   [1] Karttunen, Taivaanmekaniikka
  !!
  !! Usage:
  !!   CALL interact_full_jpl(ws, mjd_tdt, wds, error, pwds)
  !!
  !! Interaction.
  !!
  !!  ws          heliocentric, ecliptical, Cartesian coordinates for the
  !!                 massless bodies
  !!  mjd_tdt     modified Julian date (for ephemerides)
  !!  wds         evaluated force function
  !!  error       true, if reading from ephemerides fails
  !!  pwds        evaluated partial derivatives of the force function
  !!  encounters  encounters(i,j) where i=massless particle,
  !!              j=perturber, encounters(i,j)=1 or 2, where 1=inside
  !!              Roche limit and 2=within planetary radius
  !!              (~collision)
  !!
  SUBROUTINE interact_full_jpl(ws, mjd_tdt, perturbers, N, wds, error, pwds, encounters)

    REAL(prec), DIMENSION(:,:), INTENT(in)              :: ws
    REAL(prec), INTENT(in)                              :: mjd_tdt
    LOGICAL, DIMENSION(:), INTENT(in)                   :: perturbers
    INTEGER, INTENT(in)                                 :: N ! = objects with masses available
    REAL(prec), DIMENSION(:,:), INTENT(out)             :: wds
    LOGICAL, INTENT(inout)                              :: error
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: pwds 
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: encounters

    ! Number of massless bodies.
    INTEGER            :: NS

    ! Coordinates for the massive bodies.
    REAL(prec), DIMENSION(:,:), POINTER :: wc

    ! Various distances.
    REAL(prec)                 :: r2s, ir3s, ir5s
    REAL(prec), DIMENSION(N)   :: r2c, ir3c
    REAL(prec), DIMENSION(3,N) :: drs 
    REAL(prec), DIMENSION(N)   :: r2d, ir3d, ir5d
    REAL(prec)                 :: ir4s, v2s, us, ir6s, dist

    ! Utility variables.
    REAL(prec), DIMENSION(3,3) :: A, B, P1, P2 
    INTEGER :: i, j, k, l, err

    ! Get positions of massive bodies (-10 = 9 planets + Moon).
    !wc => JPL_ephemeris(mjd_tdt, perturbers(:), 11, error)
    wc => JPL_ephemeris(mjd_tdt, -10, 11, error)
    IF (error) THEN 
       DEALLOCATE(wc, stat=err)
       RETURN 
    END IF

    ! Useful quantities. 
    r2c = 0.0_prec ; ir3c = 0.0_prec
    DO i=1,N
       IF (perturbers(i)) THEN
          r2c(i)  = DOT_PRODUCT(wc(i,1:3), wc(i,1:3)) 
          ir3c(i) = 1.0_prec / (r2c(i) * SQRT(r2c(i))) 
       END IF
    END DO

    ! Number of massless bodies.
    NS = SIZE(ws,dim=2)

    ! Loop over massless bodies.
    DO i=1,NS

       r2s  = DOT_PRODUCT(ws(1:3,i), ws(1:3,i)) ! ws(6,NS) 
       ir3s = 1.0_prec / (r2s * SQRT(r2s)) 
       drs = 0.0_prec ; r2d = 0.0_prec ; ir3d = 0.0_prec
       ! Log impacts and distances to solar-system objects
       DO j=1,N
          IF (perturbers(j)) THEN
             drs(1:3,j) = wc(j,1:3) - ws(1:3,i) 
             r2d(j)     = DOT_PRODUCT(drs(1:3,j), drs(1:3,j)) 
             ir3d(j)    = 1.0_prec / (r2d(j) * SQRT(r2d(j)))
          END IF
          IF (PRESENT(encounters)) THEN
             dist = SQRT(r2d(j))
             IF (dist < planetary_radii(j)) THEN
                encounters(i,j,1) = mjd_tdt
                encounters(i,j,2) = 1          
                encounters(i,j,3) = dist          
             ELSE
                encounters(i,j,1) = mjd_tdt
                encounters(i,j,2) = 2
                encounters(i,j,3) = dist
             END IF
          END IF
       END DO
       IF (PRESENT(encounters)) THEN
          ! Distance to Sun
          dist = SQRT(DOT_PRODUCT(ws(1:3,i), ws(1:3,i)))
          IF (dist < planetary_radii(11)) THEN
             encounters(i,11,1) = mjd_tdt
             encounters(i,11,2) = 1
             encounters(i,11,3) = dist
          ELSE
             encounters(i,11,1) = mjd_tdt
             encounters(i,11,2) = 2
             encounters(i,11,3) = dist
          END IF
       END IF

       ! Initialize force function.
       wds(1:3,i) = ws(4:6,i)
       wds(4:6,i) = 0.0_prec

       ! Non-integrable part of the interaction.
       DO j=1,N
          IF (perturbers(j)) THEN
             wds(4:6,i) = wds(4:6,i) + &
                  planetary_masses(j)* (drs(1:3,j) * ir3d(j) - wc(j,1:3) * ir3c(j))
          END IF
       END DO
       wds(4:6,i) = gc * wds(4:6,i)

       ! Integrable part of the interaction.
       wds(4:6,i) = wds(4:6,i) - gc * ws(1:3,i) * ir3s

       ! Relativistic term from the Sun (isotropic coordinates).
       ! (Sitarski (1983) AcA 33)
       ir4s = 1.0_prec / (r2s ** 2)
       v2s  = DOT_PRODUCT(ws(4:6,i), ws(4:6,i))
       us   = DOT_PRODUCT(ws(1:3,i), ws(4:6,i))
       ! Optionally, no relativistic term (see global module parameter)
       IF (relativity) THEN
          wds(4:6,i) = wds(4:6,i) &
               + gc * ((4.0_prec * gc * ir4s - v2s * ir3s) * ws(1:3,i) &
               + 4.0_prec * us * ir3s * ws(4:6,i)) * ic2
       END IF

       ! Optional part.
       IF (PRESENT(pwds)) THEN

          ! Some useful quantities.
          ir5s = ir3s / r2s
          ir5d = ir3d / r2d
          ir6s = ir4s / r2s

          ! Non-integrable part of the interaction
          A = 0.0_prec
          DO j = 1, N
             IF (perturbers(j)) THEN
                DO k = 1, 3
                   DO l = 1, 3
                      A(k,l) = A(k,l) + 3.0_prec * planetary_masses(j) &
                           * drs(k,j) * drs(l,j) * ir5d(j)
                   END DO
                   A(k,k) = A(k,k) - planetary_masses(j) * ir3d(j)
                END DO
             END IF
          END DO
          A = gc * A

          ! Integrable part of the interaction
          DO k = 1, 3
             DO l = 1, 3 
                B(k,l) = 3.0_prec * ws(k,i) * ws(l,i) * ir5s
             END DO
             B(k,k) = B(k,k) - ir3s
          END DO
          B = gc * B

          ! Relativistic term from the Sun
          IF (relativity) THEN
             DO k = 1, 3
                DO l = 1, 3 
                   P1(k,l) = (3.0_prec * v2s * ir5s - &
                        16.0_prec * gc * ir6s) * ws(k,i) * ws(l,i) - &
                        12.0_prec * us * ir5s * ws(k+3,i) * ws(l,i) + &
                        4.0_prec * ir3s * ws(k+3,i) * ws(l+3,i)
                END DO
                P1(k,k) = P1(k,k) + 4.0_prec *gc * ir4s - v2s * ir3s
             END DO
             P1 = gc * ic2 * P1

             DO k = 1, 3
                DO l = 1, 3 
                   P2(k,l) = 2.0_prec * ws(k+3,i) * ws(l,i) - ws(k,i) * ws(l+3,i)
                END DO
                P2(k,k) = P2(k,k) + 2.0_prec * us
             END DO
             P2 = 2.0_prec * gc * ic2 * ir3s * P2
          ELSE
             P1 = 0.0_prec
             P2 = 0.0_prec
          END IF

          ! Partial derivatives of the force function.
          ! Put the parts together.
          pwds(:,:,i) = 0.0_prec
          DO k = 1, 3
             pwds(k,k+3,i) = 1.0_prec
          END DO
          pwds(4:6,1:3,i) = A + B + P1
          pwds(4:6,4:6,i) = P2

       END IF

    END DO

    DEALLOCATE(wc, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE interact_full_jpl






  !! Description: 
  !!   Evaluation of the full Newtonian force function for several 
  !!   massless bodies. Positions of the massive bodies are read 
  !!   from JPL ephemerides. A relativistic term due to the Sun is included.
  !!   Optional argument triggers evaluation of the partial derivatives
  !!   of the force function wrt Cartesian coordinates.
  !!
  !! References:
  !!   [1] Karttunen, Taivaanmekaniikka
  !!
  !! Usage:
  !!   CALL interact_full_jpl(ws, mjd_tdt, wds, error, pwds)
  !!
  !! Interaction.
  !!
  !!  ws       planetocentric, ecliptical (are not they equatorial??), Cartesian coordinates
  !!                 for the massless bodies
  !!  mjd_tdt  modified Julian date (for ephemerides)
  !!  wds      evaluated force function
  !!  error    true, if reading from ephemerides fails
  !!  pwds     evaluated partial derivatives of the force function
  !!
  SUBROUTINE interact_full_jpl_center(ws, mjd_tdt, wds, error, pwds)

    REAL(prec), DIMENSION(:,:), INTENT(in)              :: ws
    REAL(prec), INTENT(in)                              :: mjd_tdt
    REAL(prec), DIMENSION(:,:), INTENT(out)             :: wds
    LOGICAL, INTENT(inout)                              :: error
    REAL(prec), DIMENSION(:,:,:), OPTIONAL, INTENT(out) :: pwds 

    ! Number of massive bodies.
    INTEGER, PARAMETER :: N = 11
    ! Number of massless bodies.
    INTEGER            :: NS

    ! Coordinates for the massive bodies.
    REAL(prec), DIMENSION(:,:), POINTER :: wc_sun
    REAL(prec), DIMENSION(:,:), ALLOCATABLE :: wc
    !REAL(prec), DIMENSION(10,6) :: wc

    ! Planetary masses.
    REAL(prec), DIMENSION(N) :: m

    ! Various distances.
    REAL(prec)                 :: r2s, ir3s, ir5s
    REAL(prec), DIMENSION(N)   :: r2c, ir3c
    REAL(prec), DIMENSION(3,N) :: drs 
    REAL(prec), DIMENSION(N)   :: r2d, ir3d, ir5d
    REAL(prec)                 :: ir4s, v2s, us, ir6s

    ! Utility variables.
    REAL(prec), DIMENSION(3,3) :: A, B, P1, P2 
    INTEGER :: i, j, k, l, err


    IF (central_body /= 11 .AND. relativity) THEN
       error = .TRUE.
       ! Relativity term can currently only be used in heliocentric system.
       RETURN
    END IF

    ! Get heliocentric positions of massive bodies (-10 = 9 planets + Moon).
    wc_sun => JPL_ephemeris(mjd_tdt, -10, 11, error)
    IF (error) THEN 
       DEALLOCATE(wc_sun, stat=err)
       RETURN 
    END IF

    ALLOCATE(wc(SIZE(wc_sun,dim=1)+1,SIZE(wc_sun,dim=2)),stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       DEALLOCATE(wc_sun, wc, stat=err)
       RETURN
    END IF

    ! Change of origin
    IF (central_body /= 11) THEN
       ! Sun
       wc(N,1:3) = -1.0_prec * wc_sun(central_body, 1:3)
       ! Other massive bodies
       DO i = 1, N-1
          wc(i,1:3) = wc_sun(i,1:3) - wc_sun(central_body, 1:3)
       END DO
    ELSE
       wc(1:N-1,1:3) = wc_sun(1:N-1,1:3)
    END IF

    ! Planetary masses + Sun.
    m(1:N) = planetary_masses(1:N)

    ! Useful quantities. 
    DO i = 1, N
       r2c(i)  = DOT_PRODUCT(wc(i,1:3), wc(i,1:3)) 
       ir3c(i) = 1.0_prec / (r2c(i) * SQRT(r2c(i))) 
    END DO

    ! Number of massless bodies.
    NS = SIZE(ws,dim=2)

    ! Loop over massless bodies.
    DO i = 1, NS

       r2s  = DOT_PRODUCT(ws(1:3,i), ws(1:3,i)) ! ws(6,NS) 
       ir3s = 1.0_prec / (r2s * SQRT(r2s)) 
       DO j = 1, N 
          drs(1:3,j) = wc(j,1:3) - ws(1:3,i) 
          r2d(j)     = DOT_PRODUCT(drs(1:3,j), drs(1:3,j)) 
          ir3d(j)    = 1.0_prec / (r2d(j) * SQRT(r2d(j))) 
       END DO

       ! Initialize force function.
       wds(1:3,i) = ws(4:6,i)
       wds(4:6,i) = 0.0_prec

       ! Non-integrable part of the interaction.
       DO j = 1, N
          IF (j == central_body) CYCLE
          wds(4:6,i) = wds(4:6,i) + &
               m(j)* (drs(1:3,j) * ir3d(j) - wc(j,1:3) * ir3c(j))
          !if (j == 3) print*,m(j)* (drs(1:3,j) * ir3d(j) - wc(j,1:3) * ir3c(j))
       END DO
       wds(4:6,i) = gc * wds(4:6,i)

       ! Integrable part of the interaction.
       !wds(4:6,i) = wds(4:6,i) - gc * ws(1:3,i) * ir3s
       wds(4:6,i) = wds(4:6,i) - planetary_mu(central_body) * ws(1:3,i) * ir3s

       ! Relativistic term from the Sun (isotropic coordinates).
       ! CHECK THIS FOR NON-HELIOCENTRIC FRAME!
       ! (Sitarski (1983) AcA 33)
       ir4s = 1.0_prec / (r2s ** 2)
       v2s  = DOT_PRODUCT(ws(4:6,i), ws(4:6,i))
       us   = DOT_PRODUCT(ws(1:3,i), ws(4:6,i))
       ! Optionally, no relativistic term (see global module parameter)
       IF (relativity) THEN
          wds(4:6,i) = wds(4:6,i) &
               + gc * ((4.0_prec * gc * ir4s - v2s * ir3s) * ws(1:3,i) &
               + 4.0_prec * us * ir3s * ws(4:6,i)) * ic2
       END IF

       ! Optional part.
       IF (PRESENT(pwds)) THEN

          ! Some useful quantities.
          ir5s = ir3s / r2s
          ir5d = ir3d / r2d
          ir6s = ir4s / r2s

          ! Non-integrable part of the interaction
          A = 0.0_prec
          DO j = 1, N
             IF (j == central_body) CYCLE
             DO k = 1, 3
                DO l = 1, 3
                   A(k,l) = A(k,l) + 3.0_prec * m(j) &
                        * drs(k,j) * drs(l,j) * ir5d(j)
                END DO
                A(k,k) = A(k,k) - m(j) * ir3d(j)
             END DO
          END DO
          A = gc * A

          ! Integrable part of the interaction
          DO k = 1, 3
             DO l = 1, 3 
                B(k,l) = 3.0_prec * ws(k,i) * ws(l,i) * ir5s
             END DO
             B(k,k) = B(k,k) - ir3s
          END DO
          B = gc * B

          ! Relativistic term from the Sun
          IF (relativity) THEN
             DO k = 1, 3
                DO l = 1, 3 
                   P1(k,l) = (3.0_prec * v2s * ir5s - &
                        16.0_prec * gc * ir6s) * ws(k,i) * ws(l,i) - &
                        12.0_prec * us * ir5s * ws(k+3,i) * ws(l,i) + &
                        4.0_prec * ir3s * ws(k+3,i) * ws(l+3,i)
                END DO
                P1(k,k) = P1(k,k) + 4.0_prec *gc * ir4s - v2s * ir3s
             END DO
             P1 = gc * ic2 * P1

             DO k = 1, 3
                DO l = 1, 3 
                   P2(k,l) = 2.0_prec * ws(k+3,i) * ws(l,i) - ws(k,i) * ws(l+3,i)
                END DO
                P2(k,k) = P2(k,k) + 2.0_prec * us
             END DO
             P2 = 2.0_prec * gc * ic2 * ir3s * P2
          ELSE
             P1 = 0.0_prec
             P2 = 0.0_prec
          END IF

          ! Partial derivatives of the force function.
          ! Put the parts together.
          pwds(:,:,i) = 0.0_prec
          DO k = 1, 3
             pwds(k,k+3,i) = 1.0_prec
          END DO
          pwds(4:6,1:3,i) = A + B + P1
          pwds(4:6,4:6,i) = P2

       END IF

    END DO

    DEALLOCATE(wc, wc_sun, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       RETURN
    END IF

  END SUBROUTINE interact_full_jpl_center


END MODULE integrators




