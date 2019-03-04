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
!! *Description*:
!!
!!   Asteroid identification software for linking over apparitions.
!!
!!   More details can be found in these published works and, please,
!!   do also cite them if you make use of the code in your work:
!!
!!     Granvik, M. & Muinonen, K. (2005), 'Asteroid identification at
!!     discovery', Icarus 179, 109-127.
!!
!!     Granvik, M. & Muinonen, K. (2008), 'Asteroid identification
!!     over apparitions', Icarus 198, 130-137.
!!
!!     Granvik, M. (2007), 'Asteroid identification using statistical
!!     orbital inversion methods', PhD thesis, Department of
!!     Astronomy, Faculty of Science, University of Helsinki.
!!
!!   @author  MG 
!!   @version 2019-02-20
!!
PROGRAM orbit_linking

  USE Base_cl
  USE File_cl
  USE Time_cl
  USE Orbit_cl
  USE StochasticOrbit_cl
  USE io
  USE cl_options
  USE utilities
  USE data_structures
  USE sort

#ifdef __INTEL_COMPILER
  USE IFPORT
#endif

  IMPLICIT NONE
  TYPE (rb_tree_i8_i4arr), POINTER :: address_tree
  TYPE (rb_tree_ch32_8r8), POINTER :: tree_2bmc
  TYPE (rb_tree_ch32), POINTER :: &
       pos3links_tree, &
       tree_c2bmc, &
       tree_e2bmc, &
       tree_c2blsl, &
       tree_e2blsl, &
       tree_cnblsl, &
       tree_enblsl, &
       tree_2blsl, &
       tree_nblsl
  TYPE (rb_tree_node_i8_i4arr), POINTER :: address_tree_node1, &
       address_tree_node2
  TYPE (rb_tree_node_ch32), POINTER :: tree_node_ch32, &
       tree_node_ch32_1, tree_node_ch32_2
  TYPE (rb_tree_node_ch32_8r8), POINTER :: tree_node_ch32_8r8, &
       tree_node_ch32_8r8_1, tree_node_ch32_8r8_2
  TYPE (File) :: &
       logfile, &
       obsfile, &
       orbidfile, &
       orbfile, &
       tmpfile, &
       addressfile, &
       lnksfile
  TYPE (Time), DIMENSION(3) :: &
       t_arr
  TYPE (Time) :: &
       epoch, &
       t0, &
       t, &
       t_inv
  TYPE (SphericalCoordinates) :: &
       scoord
  TYPE (Observation), DIMENSION(:), POINTER :: &
       obs_arr_all_first
  TYPE (Observation) :: &
       obs
  TYPE (Observations), DIMENSION(:), POINTER :: &
       obss_sep_arr, &
       obss_sep_arr_
  TYPE (Observations) :: &
       obss, &
       obss_, &
       obss_1, &
       obss_2, &
       obss_all, &
       obss_all_, &
       obss_selected, &
       obss_all_first, &
       obss_nb
  TYPE (Orbit), DIMENSION(:), POINTER :: &
       orb_arr
  TYPE (Orbit) :: &
       orb_mc, &
       orb, &
       orb_init
  TYPE (StochasticOrbit) :: &
       storb
  CHARACTER(len=FNAME_LEN) :: &
       fname, &
       obsfname, &
       obsallfname, &
       addressfname, &
       lnksfname
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:,:), POINTER :: &
       correct_linkages
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: &
       id_arr, &
       obs_set_ids
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), ALLOCATABLE :: &
       id_arr_prm, &
       id_arr_
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(3) :: &
       id_arr3
  CHARACTER(len=DESIGNATION_LEN) :: &
       id1, id2, id3, des1
  CHARACTER(len=ELEMENT_TYPE_LEN) :: &
       element_type_in_prm, &
       element_type_nblsl
  CHARACTER(len=DYN_MODEL_LEN) :: &
       dyn_model_step1, &
       dyn_model_step2, &
       dyn_model_step3
  CHARACTER(len=INTEGRATOR_LEN) :: &
       integrator
  CHARACTER(len=1024), DIMENSION(4) :: &
       header
  CHARACTER(len=64) :: &
       orbid, &
       comparison_variable_type, &
       obs_type, &
       key3, &
       str, &
       pidstr
  REAL(bp), DIMENSION(:,:), POINTER :: &
       residuals
  REAL(bp), DIMENSION(6,2) :: &
       bounds
  REAL(bp), DIMENSION(:), ALLOCATABLE :: &
       mjds0
  REAL(bp), DIMENSION(8) :: &
       rarr8
  REAL(bp), DIMENSION(6) :: &
       elements, &
       elements_, &
       delements, &
       width, &
       rms, &
       rans, &
       stdevs
  REAL(bp), DIMENSION(3) :: &
       pos
  REAL(bp) :: &
       integration_step, &
       computed_upper_bound, &
       mjd0, &
       tt_mjd, &
       res_max_mc, &
       rms_max_2blsl, &
       rms_max_nblsl, &
       res_max_addlink, &
       ls_correction_factor, &
       dt_min, & 
       dt1, dt2, &
       rms_val, &
       rms_2blsl_step3_min, &
       rms_2blsl_step3_max
  INTEGER(ihp) :: indx
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       sets
  INTEGER, DIMENSION(6) :: box_nrs
  INTEGER, DIMENSION(3) :: i4arr
  INTEGER :: &
       filter_start, &
       filter_stop, &
       ifilter, &
       err, &
       i, &
       j, &
       k, &
       l, &
       m, &
       n, &
       ii, &
       jj, &
       kk, &
       nobj, &
       norb, &
       norb_mc, &
       iorb, &
       lu, &
       nlinks, &
       nlinks_, &
       nlink_1, nlink_2, &
       nnolinks, &
       nclink, &
       nelink, &
       naddress, &
       nset, &
       nset_max, &
       loop_start, &
       loop_stop, &
       pid, &
       info_verb_, &
       ls_niter_major_min, &
       ls_niter_major_max, &
       ls_niter_minor, &
       nr
  LOGICAL, DIMENSION(10) :: &
       perturbers
  LOGICAL, DIMENSION(6) :: &
       ls_element_mask, &
       elm
  LOGICAL :: &
       common_epoch, &
       at_least_one_correct, &
       insert_header
  INTEGER :: ret

  IF (get_cl_option("--help", .FALSE.)) THEN
     WRITE(stdout,*) "orbit_linking options:"
     WRITE(stdout,*) ""
     WRITE(stdout,*) "  --help"
     WRITE(stdout,*) "  --logfile"
     WRITE(stdout,*) "  --filter-start"
     WRITE(stdout,*) "  --filter-stop"
     WRITE(stdout,*) "  --orbid"
     WRITE(stdout,*) "  --addressfile"
     WRITE(stdout,*) "  --obsall"
     WRITE(stdout,*) "  --2blnks"
     WRITE(stdout,*) "  --loop-start"
     WRITE(stdout,*) "  --loop-stop"
     WRITE(stdout,*) "  --2borb"
     WRITE(stdout,*) "  --nbobs"
     WRITE(stdout,*) "  --nborb"
     WRITE(stdout,*) "  --nblsl"
     WRITE(stdout,*) ""
     STOP
  END IF

  ! Get process ID
  pid = getpid()
  CALL toString(pid, pidstr, error)
  IF (error) THEN
     CALL errorMessage("orbit_linking", &
          "TRACE BACK (1)", 1)
     STOP
  END IF

  ! PARAMETERS
  simulated_observations = .FALSE.
  !simulated_observations = .TRUE.

  stdevs = -1.0_bp*rad_asec
  fname = get_cl_option("--logfile=", " ")
  IF (LEN_TRIM(fname) == 0) THEN
     lu = stdout
  ELSE
     CALL NEW(logfile, TRIM(fname))
     IF (error) THEN
        CALL errorMessage("orbit_linking", &
             "TRACE BACK (5)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     CALL setPositionAppend(logfile)
     IF (error) THEN
        CALL errorMessage("orbit_linking", &
             "TRACE BACK (10)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     CALL OPEN(logfile)
     IF (error) THEN
        CALL errorMessage("orbit_linking", &
             "TRACE BACK (10)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     lu = getUnit(logfile)
     IF (error) THEN
        CALL errorMessage("orbit_linking", &
             "TRACE BACK (15)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
  END IF
  obs_type = "mpc3"
  common_epoch = .FALSE.
  dyn_model_step1 = "2-body"
  dyn_model_step2 = "2-body"
  norb_mc = 100000 ! NEO
  !norb_mc = 10000 ! MBO
  res_max_mc = 4.0_bp*rad_deg ! 4 deg during tests incl in submitted MS
  rms_max_2blsl = 100.0_bp*rad_asec !100.0_bp*rad_asec
  element_type_nblsl = "keplerian"
  rms_max_nblsl = 1.5_bp*rad_asec
  rms_2blsl_step3_min = 0.0_bp ! as
  rms_2blsl_step3_max = 10.0_bp ! as
  res_max_addlink = 5.0_bp*rad_amin
  dt_min = 0.0_bp
  dyn_model_step3 = "n-body"
  perturbers = .TRUE.
  integrator = "bulirsch-stoer"
  integration_step = 10.0_bp
  info_verb = 1
  err_verb = 1
  CALL NEW(epoch, 2000, 1, 1.0_bp, "TT")
  IF (error) THEN
     CALL errorMessage("orbit_linking", &
          "TRACE BACK (20)", 1)
     CALL NULLIFY(logfile)
     STOP
  END IF

!!$  ! Keplerian
!!$  comparison_variable_type = "keplerian"
!!$  elm = .TRUE.
!!$  elm(6) = .FALSE.
!!$  bounds(:,1) = (/   0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp /)
!!$  bounds(:,2) = (/ 500.0_bp, 1.0_bp, pi, two_pi, two_pi, two_pi /)
!!$  width = (/ 0.1_bp, 0.1_bp, 0.5_bp*rad_deg, 5.0_bp*rad_deg, &
!!$       45.0_bp*rad_deg, 10.0_bp*rad_deg /)

  ! Modified Keplerian
  comparison_variable_type = "modified keplerian"
  elm = .TRUE.
  bounds(:,1) = (/   0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 0.0_bp, 46066.0_bp /)
  bounds(:,2) = (/  50.0_bp, 1.0_bp, pi, two_pi, two_pi, 57023.0_bp /)
  width = (/ 0.1_bp, 0.1_bp, 0.5_bp*rad_deg, 10.0_bp*rad_deg, &
       45.0_bp*rad_deg, 180.0_bp /)

!!$  ! Equinoctial
!!$  comparison_variable_type = "equinoctial"
!!$  elm = .TRUE.
!!$  bounds(:,1) = (/   0.0_bp, -1.0_bp, -1.0_bp, -100.0_bp, -100.0_bp, 46066.0_bp /)
!!$  bounds(:,2) = (/  10.0_bp*rad_deg,  1.0_bp,  1.0_bp,  100.0_bp,  100.0_bp, 57023.0_bp /)
!!$  width =         (/ 0.01_bp*rad_deg,  0.05_bp, 0.05_bp,   0.05_bp,   0.05_bp, 30.0_bp /)

  box_nrs = CEILING((bounds(:,2) - bounds(:,1))/width)
  DO i=1,SIZE(box_nrs)
     computed_upper_bound = bounds(i,1) + box_nrs(i)*width(i)
     IF (ABS(computed_upper_bound - bounds(i,2)) > &
          1000.0_bp*EPSILON(computed_upper_bound)) THEN
        box_nrs(i) = box_nrs(i) - 1
     END IF
  END DO
  ! Which filters are to be used?
  filter_start = get_cl_option("--filter-start=", 1)
  filter_stop = get_cl_option("--filter-stop=", 4)

  ! Write logfile
  WRITE(lu,"(1X,A,L1)") "Simulated observations... ", simulated_observations
  WRITE(lu,"(1X,A,6(F10.5,1X))") "Assumed observational RA and Dec uncertainties [asec]:", stdevs(2:3)/rad_asec
  WRITE(lu,"(1X,2(A,1X,I0,1X))") "Using filters", filter_start, "through", filter_stop
  WRITE(lu,"(1X,A,1X,F7.3)") "Maximum O-C residual in 2bMC [deg]:", res_max_mc/rad_deg
  WRITE(lu,"(1X,A,1X,I0)") "Number of trial orbits in 2bMC:", norb_mc
  WRITE(lu,"(1X,A,1X,F7.3)") "Maximum O-C residual rms in 2bLSL [arcsec]:", rms_max_2blsl/rad_asec
  WRITE(lu,"(1X,A,1X,F7.3)") "Maximum O-C residual rms in nbLSL [arcsec]:", rms_max_nblsl/rad_asec
  WRITE(lu,"(1X,A,1X,F7.3)") "Maximum O-C residual in search for additional linkages [arcmin]:", res_max_addlink/rad_amin

  DO ifilter=filter_start,filter_stop

     SELECT CASE (ifilter)

     CASE (1)

        !!
        !! PHASE 1: FIND CANDIDATE LINKAGES (VIA ADDRESS COMPARISON OR FROM FILE)
        !!

        WRITE(lu,*)
        WRITE(lu,"(1X,A)") "Using the following parameters for the address comparison:"
        WRITE(lu,"(1X,A,A)") " - comparison variable : ", TRIM(ADJUSTL(comparison_variable_type))
        WRITE(lu,"(1X,A,6(6X,L1,4X))") " - elements to be used : ", elm
        IF (comparison_variable_type == "keplerian") THEN
           WRITE(lu,"(1X,A,6(F10.4,1X))") " - lower bounds        : ", bounds(1:2,1), bounds(3:6,1)/rad_deg 
           WRITE(lu,"(1X,A,6(F10.4,1X))") " - upper bounds        : ", bounds(1:2,2), bounds(3:6,2)/rad_deg 
           WRITE(lu,"(1X,A,6(F10.4,1X))") " - bin widths          : ", width(1:2), width(3:6)/rad_deg
           WRITE(lu,"(1X,A,6(I10,1X))")   " - bin numbers         : ", box_nrs
        ELSE IF (comparison_variable_type == "modified keplerian") THEN
           WRITE(lu,"(1X,A,6(F13.4,1X))") " - lower bounds        : ", bounds(1:2,1), bounds(3:5,1)/rad_deg, bounds(6,1)
           WRITE(lu,"(1X,A,6(F13.4,1X))") " - upper bounds        : ", bounds(1:2,2), bounds(3:5,2)/rad_deg, bounds(6,2)
           WRITE(lu,"(1X,A,6(F13.4,1X))") " - bin widths          : ", width(1:2), width(3:5)/rad_deg, width(6)
           WRITE(lu,"(1X,A,6(I13,1X))")   " - bin numbers         : ", box_nrs
        ELSE
           WRITE(lu,*) " - lower bounds        : ", bounds(:,1) 
           WRITE(lu,*) " - upper bounds        : ", bounds(:,2) 
           WRITE(lu,*) " - bin widths          : ", width 
           WRITE(lu,*) " - bin numbers         : ", box_nrs
        END IF
        IF (common_epoch) THEN
           WRITE(lu,"(1X,A,A)") " - common epoch        : ", getCalendarDateString(epoch, "TT")
        END IF
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1005)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        WRITE(lu,*)
        indx = 1_ihp
        DO i=1,6
           IF (.NOT.elm(i)) THEN
              CYCLE
           END IF
           indx = indx * INT(box_nrs(i),ihp)
        END DO
        WRITE(lu,*)
        WRITE(lu,"(1X,A,1X,I0)") "Number of addresses in the binned phase space is", indx
        WRITE(lu,*)

        ! READ FILE CONTAINING ORBIT FILENAMES:
        fname = get_cl_option("--orbid=", " ")
        IF (LEN_TRIM(fname) == 0) THEN
           error = .TRUE.
           CALL errorMessage("orbit_linking", &
                "Orbit ID file missing.", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NEW(orbidfile, TRIM(fname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(orbidfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(orbidfile)  
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1020)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        nset_max = getNrOfLines(orbidfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1025)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        nset = 0
        DO
           READ(getUnit(orbidfile),"(A)",iostat=err) orbid
           IF (err > 0) THEN
              error = .TRUE.
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1030)", 1)
              CALL NULLIFY(logfile)
              STOP
           ELSE IF (err < 0) THEN
              EXIT
           END IF
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " ID of set read."
           END IF
           ret = system("gunzip -c " // TRIM(orbid) // ".orb.gz > " // TRIM(orbid) // ".orb")
           CALL NEW(orbfile,TRIM(orbid) // ".orb")
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1035)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setStatusOld(orbfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1040)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(orbfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1045)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Orbit file opened."
           END IF
           norb = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1050)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ALLOCATE(id_arr_prm(norb), &
                orb_arr(norb), &
                stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not allocate memory (1005).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           id_arr_prm = " "
           header(1:4)(:) = " "
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Reading orbit file."
           END IF
           DO i=1,norb
              CALL readOpenOrbOrbitFile(getUnit(orbfile), header, &
                   element_type_in=element_type_in_prm, &
                   id=id_arr_prm(i), &
                   orb=orb_arr(i))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not read orbit file.", 1)
                 WRITE(lu,*) i
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL setParameters(orb_arr(i), &
                   dyn_model=dyn_model_step1, &
                   perturbers=perturbers, &
                   integrator=integrator, &
                   integration_step=integration_step)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1055)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL NULLIFY(orbfile)
           ret = system("rm -f " // TRIM(orbid) // ".orb")
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Orbit file read."
           END IF

           ! Propagate to the comparison epoch if necessary
           t = getTime(orb_arr(1))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1060)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (.NOT.equal(t,epoch) .AND. common_epoch) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Propagation needed."
              CALL propagate(orb_arr, epoch)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1065)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              IF (info_verb >= 1) THEN
                 WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Propagation ready."
              END IF
           END IF

           ! Compute addresses of the orbits and put them into the
           ! comparison table:
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Compute addresses."
           END IF
           nset = nset + 1
           norb = SIZE(orb_arr, dim=1)

           IF (.NOT.ASSOCIATED(id_arr)) THEN
              ALLOCATE(id_arr(nset_max), stat=err)
              IF (err /= 0) THEN
                 error = .TRUE.
                 CALL errorMessage("orbit_linking", &
                      'Could not allocate memory (1010).',1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              id_arr = " "
              CALL init(address_tree)
           END IF

!!!$omp parallel do &
!!!$omp private(i, elements, t, tt_mjd, error, scoord, orb, pos, err_verb, indx)
           orbloop:DO iorb=1,norb
              SELECT CASE (comparison_variable_type)
              CASE ("modified keplerian", "keplerian", "delaunay", &
                   "poincare", "equinoctial")
                 err_verb = 0
              END SELECT
              SELECT CASE (comparison_variable_type)
              CASE ("keplerian", "delaunay", "poincare", "equinoctial", "cartesian")
                 elements = getElements(orb_arr(iorb), comparison_variable_type)
                 elements(1) = SQRT(planetary_mu(11)/elements(1)**3.0_bp)
                 t = getPlaneCrossingTime(orb_arr(iorb), epoch)
                 tt_mjd = getMJD(t, "TT")
                 CALL NULLIFY(t)
                 elements(6) = tt_mjd
              CASE ("modified keplerian")
                 elements = getElements(orb_arr(iorb), "keplerian")
                 ! Temp solution for NEO analysis
                 IF (elements(1) > 5.5_bp) THEN
                    CYCLE orbloop
                 END IF
                 IF (elements(3) == pi) THEN
                    elements(3) = 0.0_bp
                 END IF
                 IF (elements(4) == two_pi) THEN
                    elements(4) = 0.0_bp
                 END IF
                 IF (elements(5) == two_pi) THEN
                    elements(5) = 0.0_bp
                 END IF
                 IF (elements(6) == two_pi) THEN
                    elements(6) = 0.0_bp
                 END IF
                 t = getPlaneCrossingTime(orb_arr(iorb), epoch)
                 tt_mjd = getMJD(t, "TT")
                 CALL NULLIFY(t)
                 elements(6) = tt_mjd
                 IF (error) THEN
                    error = .FALSE.
                    CYCLE orbloop
                 END IF
              CASE ("k-vector")
                 elements = 0.0_bp
                 elements(1:3) = getKVector(orb_arr(iorb))
                 scoord = getSCoord(orb_arr(iorb))
                 elements(4:6) = getPosition(scoord)
                 CALL NULLIFY(scoord)
              CASE ("sphpos")
                 orb = copy(orb_arr(iorb))
                 CALL setParameters(orb, &
                      dyn_model="2-body")
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (1070)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 DO i=1,3
                    CALL propagate(orb, t_arr(i))
                    IF (error) THEN
                       error = .FALSE.
                       CYCLE orbloop
                    END IF
                    scoord = getSCoord(orb)
                    pos = getPosition(scoord)
                    elements((i-1)*2+1:(i-1)*2+2) = pos(2:3)
                    CALL NULLIFY(scoord)
                 END DO
              END SELECT
              SELECT CASE (comparison_variable_type)
              CASE ("modified keplerian", "keplerian", "delaunay", &
                   "poincare", "equinoctial")
                 err_verb = 1
              END SELECT
              IF (error) THEN
                 SELECT CASE (comparison_variable_type)
                 CASE ("modified keplerian", "keplerian", "delaunay", &
                      "poincare", "equinoctial")
                    error = .FALSE.
                    CYCLE
                 END SELECT
                 CALL errorMessage("orbit_linking", &
                      'TRACE BACK (1075)',1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              indx = arrayToInteger(elements, elm, width, box_nrs, &
                   bounds, error)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      'Could not transform array to integer.', 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              IF (indx == 0_ihp) THEN
                 WRITE(lu,*) elements(1:2), elements(3:6)/rad_deg
              END IF
!!!$omp critical
              CALL insert_tree_node(address_tree, indx, nset)
!!!!!$omp end critical
           END DO orbloop
!!!!!$omp end parallel do
           id_arr(nset) = id_arr_prm(1)
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // &
                   " Computation of addresses ready."
           END IF

           DEALLOCATE(id_arr_prm, &
                orb_arr, &
                stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (1005).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
        END DO
        id_arr => reallocate(id_arr, nset)

        ret = system("ps ux | grep orbit_linking")

        key3 = " "
        m = 0
        addressfname = get_cl_option("--addressfile=","address.out")
        CALL NEW(addressfile, TRIM(addressfname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1080)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(addressfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (1085)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        IF (simulated_observations) THEN
           CALL init(pos3links_tree)
        END IF
        address_tree_node1 => minimum(address_tree, address_tree%root)
        DO WHILE (.NOT.ASSOCIATED(address_tree_node1, address_tree%nil))
           m = m + SIZE(address_tree_node1%data_nodes)
           DO k=1,SIZE(address_tree_node1%data_nodes)
              nset = SIZE(address_tree_node1%data_nodes(k)%data_array)
              WRITE(getUnit(addressfile),"(I15,1X,I10)",advance="NO") &
                   address_tree_node1%data_nodes(k)%key, nset
              DO i=1,nset
                 WRITE(getUnit(addressfile),"(1X,A16)",advance="NO") &
                      TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(i)))
              END DO
              WRITE(getUnit(addressfile),"(1X)")
              IF (simulated_observations) THEN
                 DO i=1,nset-2
                    DO j=i+1,nset-1
                       DO l=j+1,nset
                          id_arr3(1) = TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(i)))
                          id_arr3(2) = TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(j)))
                          id_arr3(3) = TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(l)))
                          CALL insertionSort(id_arr3)
                          key3 = TRIM(id_arr3(1)) // " " // &
                               TRIM(id_arr3(2)) // " " // &
                               TRIM(id_arr3(3))
                          tree_node_ch32 => search(pos3links_tree, key3)
                          IF (ASSOCIATED(tree_node_ch32,pos3links_tree%nil)) THEN
                             CALL insert_tree_node(pos3links_tree, key3)
                          END IF
                       END DO
                    END DO
                 END DO
              END IF
              DEALLOCATE(address_tree_node1%data_nodes(k)%data_array, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (1010).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           address_tree_node2 => address_tree_node1
           address_tree_node1 => successor(address_tree, address_tree_node1)
           address_tree_node2 => delete_tree_node(address_tree, address_tree_node2)
           DEALLOCATE(address_tree_node2%data_nodes, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (1015).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           DEALLOCATE(address_tree_node2, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (1020).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
        END DO
        CALL delete_tree(address_tree)
        CALL NULLIFY(addressfile)
        WRITE(lu,*) m, " addresses used."

        IF (simulated_observations) THEN
           obsfname = get_cl_option("--obsall=", " ")
           IF (LEN_TRIM(obsfname) /= 0) THEN
              CALL NEW(obsfile, TRIM(obsfname))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1090)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL setStatusOld(obsfile)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1095)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              WRITE(lu,"(A,A)") "Open observation file ", TRIM(obsfname)
              CALL OPEN(obsfile)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1100)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              WRITE(lu,"(A,A)") "Read observation file ", TRIM(obsfname)
              IF (ANY(stdevs < 0.0_bp)) THEN
                 CALL NEW(obss, obsfile)
              ELSE
                 CALL NEW(obss, obsfile, stdev=stdevs)
              END IF
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1105)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NULLIFY(obsfile)
           ELSE
              WRITE(0,*) "Could not read observation file ", TRIM(obsfname)
              CALL NULLIFY(logfile)
              STOP
           END IF
           WRITE(lu,"(A)") "Observation file read."
           obss_sep_arr => getSeparatedSets(obss)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (1110)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           nset = SIZE(obss_sep_arr,dim=1)
           WRITE(lu,"(1X,A,I0,A)") "Altogether ", &
                nset, " observation sets."
           CALL NULLIFY(obss)
           WRITE(lu,"(1X,A)") "Getting designations and searching " // &
                "all correct linkages between the sets."
           l = 0
           DO i=1,nset-2
              des1 = getDesignation(obss_sep_arr(i))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1115)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              id1 = getID(obss_sep_arr(i))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (1120)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              DO j=i+1,nset-1
                 IF (des1 == getDesignation(obss_sep_arr(j))) THEN
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (1125)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    id2 = getID(obss_sep_arr(j))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (1130)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    DO k=j+1,nset
                       IF (des1 == getDesignation(obss_sep_arr(k))) THEN
                          IF (error) THEN
                             CALL errorMessage("orbit_linking", &
                                  "TRACE BACK (1135)", 1)
                             CALL NULLIFY(logfile)
                             STOP
                          END IF
                          id3 = getID(obss_sep_arr(k))
                          IF (error) THEN
                             CALL errorMessage("orbit_linking", &
                                  "TRACE BACK (1140)", 1)
                             CALL NULLIFY(logfile)
                             STOP
                          END IF
                          l = l + 1
                          IF (.NOT.ASSOCIATED(correct_linkages)) THEN
                             ALLOCATE(correct_linkages(l,3), stat=err)
                             IF (err /= 0) THEN
                                CALL errorMessage("orbit_linking", &
                                     "Could not allocate memory (1015).", 1)
                                CALL NULLIFY(logfile)
                                STOP
                             END IF
                          ELSE IF (l > SIZE(correct_linkages,dim=1)) THEN
                             correct_linkages => reallocate(correct_linkages, 2*l, 3)
                          END IF
                          id_arr3(1) = TRIM(id1)
                          id_arr3(2) = TRIM(id2)
                          id_arr3(3) = TRIM(id3)
                          CALL insertionSort(id_arr3)
                          correct_linkages(l,1) = TRIM(id_arr3(1))
                          correct_linkages(l,2) = TRIM(id_arr3(2))
                          correct_linkages(l,3) = TRIM(id_arr3(3))
                       END IF
                    END DO
                 END IF
              END DO
              CALL NULLIFY(obss_sep_arr(i))
           END DO
           CALL NULLIFY(obss_sep_arr(nset-1))
           CALL NULLIFY(obss_sep_arr(nset))        
           DEALLOCATE(obss_sep_arr, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (1025).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           WRITE(lu,"(1X,A,I0,A)") "Altogether ", l, " correct 3-linkages to be detected."
           correct_linkages => reallocate(correct_linkages, l, 3)
           j = 0
           tree_node_ch32 => minimum(pos3links_tree, pos3links_tree%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32, pos3links_tree%nil))
              tree_node_ch32 => successor(pos3links_tree, tree_node_ch32)
              j = j + 1
           END DO
           WRITE(lu,*) "3-links found: ", j
           j = 0
           DO i=1,l
              key3 = TRIM(correct_linkages(i,1)) // " " // &
                   TRIM(correct_linkages(i,2)) // " " // &
                   TRIM(correct_linkages(i,3))
              tree_node_ch32 => search(pos3links_tree, key3)
              IF (ASSOCIATED(tree_node_ch32,pos3links_tree%nil)) THEN
                 WRITE(lu,*) "Missing: ", TRIM(key3)
                 j = j + 1
              END IF
           END DO
           WRITE(lu,*) "Missing correct 3-links: ", j
           tree_node_ch32_1 => minimum(pos3links_tree, pos3links_tree%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, pos3links_tree%nil))
              tree_node_ch32_2 => tree_node_ch32_1
              tree_node_ch32_1 => successor(pos3links_tree, tree_node_ch32_1)
              tree_node_ch32_2 => delete_tree_node(pos3links_tree, tree_node_ch32_2)
              DEALLOCATE(tree_node_ch32_2, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (1030).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL delete_tree(pos3links_tree)
           DEALLOCATE(correct_linkages, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (1035).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF

        END IF
        DEALLOCATE(id_arr, stat=err)


     CASE (2)

        !!
        !! PHASE 2: FIND COMMON 2B ORBIT
        !!

        obsfname = get_cl_option("--obsall=", " ")
        IF (LEN_TRIM(obsfname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Observation file missing.", 1)
           CALL NULLIFY(logfile)
           STOP     
        END IF
        addressfname = get_cl_option("--addressfile="," ")
        IF (LEN_TRIM(addressfname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Address file missing.", 1)
           CALL NULLIFY(logfile)
           STOP     
        END IF
        CALL NEW(addressfile, TRIM(addressfname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2005)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(addressfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(addressfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        IF (simulated_observations) THEN
           CALL init(tree_c2bmc)
           CALL init(tree_e2bmc)
           CALL init(tree_c2blsl)
           CALL init(tree_e2blsl)
        ELSE
           CALL init(tree_2blsl)
        END IF
        nlinks = 0
        fname = get_cl_option("--2blnks-found=", " ")
        IF (LEN_TRIM(fname) /= 0) THEN
           IF (simulated_observations) THEN
              CALL errorMessage("orbit_linking", &
                   "Using the --2blnks-found option for simulated " // &
                   "data may produce erroneous results.", 1)
           END IF
           CALL NEW(tmpfile, TRIM(fname))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2020)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(tmpfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2025)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           nlinks = getNrOfLines(tmpfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2030)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (nlinks /= 0) THEN
              WRITE(lu,"(3(1X,A),1X,I0,1X,A)") "1st 2-b links file", &
                   TRIM(fname), "contains", nlinks, &
                   "lines. Reading the links."
              DO i=1,nlinks
                 READ(getUnit(tmpfile),*) id1, rms_val, id_arr3
                 key3 = TRIM(id_arr3(1)) // " " // TRIM(id_arr3(2)) // " " // TRIM(id_arr3(3))
                 CALL insert_tree_node(tree_2blsl, TRIM(key3))
                 WRITE(*,*) i, id1, rms_val, TRIM(key3)
              END DO
           END IF
           CALL NULLIFY(tmpfile)
        END IF
        fname = get_cl_option("--2blnks=", " ")
        IF (LEN_TRIM(fname) /= 0) THEN
           IF (simulated_observations) THEN
              CALL errorMessage("orbit_linking", &
                   "Using the --2blnks option for simulated " // &
                   "data may produce erroneous results.", 1)
           END IF
           CALL NEW(tmpfile, TRIM(fname))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2035)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(tmpfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2040)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           nlinks_ = getNrOfLines(tmpfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2045)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (nlinks_ /= 0) THEN
              WRITE(lu,"(3(1X,A),1X,I0,1X,A)") "2nd 2-b links file", &
                   TRIM(fname), "contains", nlinks_, &
                   "lines. Reading the links."
              DO i=1,nlinks_
                 READ(getUnit(tmpfile),*) id1, rms_val, id_arr3
                 key3 = TRIM(id_arr3(1)) // " " // TRIM(id_arr3(2)) // " " // TRIM(id_arr3(3))
                 CALL insert_tree_node(tree_2blsl, TRIM(key3))
                 WRITE(*,*) i, id1, rms_val, TRIM(key3)
              END DO
              nlinks = nlinks + nlinks_
           END IF
           CALL NULLIFY(tmpfile)
        END IF
        ! Check whether a header is needed for the output orbit file:
        fname = get_cl_option("--2borb=", "2b_rms_acc.orb")
        CALL NEW(tmpfile,TRIM(fname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2050)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setPositionAppend(tmpfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2060)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(tmpfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2055)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        i = getNrOfLines(tmpfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (2065)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        IF (i == 0) THEN
           insert_header = .TRUE.
        ELSE
           insert_header = .FALSE.
        END IF
        CALL NULLIFY(tmpfile)

        WRITE(lu,*)
        WRITE(lu,"(1X,A)") "Using the following parameters for extracting orbital elements from the addresses:"
        WRITE(lu,"(1X,A,A)") " - comparison variable : ", TRIM(ADJUSTL(comparison_variable_type))
        WRITE(lu,"(1X,A,6(6X,L1,4X))") " - elements to be used : ", elm
        IF (comparison_variable_type == "keplerian") THEN
           WRITE(lu,"(1X,A,6(F10.4,1X))") " - lower bounds        : ", bounds(1:2,1), bounds(3:6,1)/rad_deg 
           WRITE(lu,"(1X,A,6(F10.4,1X))") " - upper bounds        : ", bounds(1:2,2), bounds(3:6,2)/rad_deg 
           WRITE(lu,"(1X,A,6(F10.4,1X))") " - bin widths          : ", width(1:2), width(3:6)/rad_deg
           WRITE(lu,"(1X,A,6(I10,1X))")   " - bin numbers         : ", box_nrs
        ELSE IF (comparison_variable_type == "modified keplerian") THEN
           WRITE(lu,"(1X,A,6(F13.4,1X))") " - lower bounds        : ", bounds(1:2,1), bounds(3:5,1)/rad_deg, bounds(6,1)
           WRITE(lu,"(1X,A,6(F13.4,1X))") " - upper bounds        : ", bounds(1:2,2), bounds(3:5,2)/rad_deg, bounds(6,2)
           WRITE(lu,"(1X,A,6(F13.4,1X))") " - bin widths          : ", width(1:2), width(3:5)/rad_deg, width(6)
           WRITE(lu,"(1X,A,6(I13,1X))")   " - bin numbers         : ", box_nrs
        ELSE
           WRITE(lu,*) " - lower bounds        : ", bounds(:,1) 
           WRITE(lu,*) " - upper bounds        : ", bounds(:,2) 
           WRITE(lu,*) " - bin widths          : ", width 
           WRITE(lu,*) " - bin numbers         : ", box_nrs
        END IF
        WRITE(lu,*)
        naddress = 0
        loop_start = get_cl_option("--loop-start=", 1)
        loop_stop = get_cl_option("--loop-stop=", HUGE(loop_stop))

        ! Scan addresses one by one
        addressloop: DO
           naddress = naddress + 1
           WRITE(stdout,*) naddress
           IF (naddress > loop_stop) THEN
              EXIT addressloop
           END IF

           ! Read address
           READ(getUnit(addressfile),"(I15,1X,I10)",advance="NO",iostat=err) indx, nset
           IF (err > 0) THEN
              CALL NULLIFY(logfile)
              STOP "Error while reading address tree file."
           ELSE IF (err < 0) THEN
              EXIT addressloop
           END IF
           IF (nset == 1) THEN
              READ(getUnit(addressfile),"(1(1X,A16))") id1
              CYCLE addressloop
           ELSE IF (nset == 2) THEN
              READ(getUnit(addressfile),"(2(1X,A16))") id1, id2
              CYCLE addressloop
           END IF
           CALL toString(nset, str, error)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2020)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ALLOCATE(id_arr_(nset), stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not allocate memory (2005)", 1)
              CALL NULLIFY(logfile)
              STOP        
           END IF
           READ(getUnit(addressfile),"(" // TRIM(str) // "(1X,A16))") id_arr_
           IF (naddress < loop_start) THEN
              DEALLOCATE(id_arr_, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2001).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CYCLE addressloop
           END IF
           ! Address read

           IF (indx < 1_ihp) THEN
              DEALLOCATE(id_arr_, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2001).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CYCLE addressloop
           END IF
           WRITE(lu,"(1X,2(A,I0))") "Analyzing address #", naddress, ": ", indx
           ALLOCATE(mjds0(nset), stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not allocate memory (2006)", 1)
              CALL NULLIFY(logfile)
              STOP        
           END IF

           ! If only one possible 3-link available, check whether it
           ! is already found:
           IF (nset == 3) THEN
              CALL insertionSort(id_arr_)
              key3 = TRIM(ADJUSTL(id_arr_(1))) // " " // &
                   TRIM(ADJUSTL(id_arr_(2))) // " " // &
                   TRIM(ADJUSTL(id_arr_(3)))
              IF (simulated_observations) THEN
                 tree_node_ch32_1 => search(tree_c2blsl, key3)
                 tree_node_ch32_2 => search(tree_e2blsl, key3)
                 IF (.NOT.ASSOCIATED(tree_node_ch32_1,tree_c2blsl%nil) .OR. &
                      .NOT.ASSOCIATED(tree_node_ch32_2,tree_e2blsl%nil)) THEN
                    DEALLOCATE(id_arr_, mjds0, stat=err)
                    IF (err /= 0) THEN
                       CALL errorMessage("orbit_linking", &
                            "Could not deallocate memory (2005).", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CYCLE addressloop
                 END IF
              ELSE
                 tree_node_ch32 => search(tree_2blsl, key3)
                 IF (.NOT.ASSOCIATED(tree_node_ch32,tree_2blsl%nil)) THEN
                    DEALLOCATE(id_arr_, mjds0, stat=err)
                    IF (err /= 0) THEN
                       CALL errorMessage("orbit_linking", &
                            "Could not deallocate memory (2010).", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CYCLE addressloop
                 END IF
              END IF
           END IF

           ! Read the observations refering to the address:
           ret = system("rm -f orbit_linking_tmp" // &
                TRIM(pidstr) // "." // TRIM(obs_type))
           DO i=1,SIZE(id_arr_)
              ret = system("grep " // TRIM(id_arr_(i)) // &
                   " " // TRIM(obsfname) // &
                   " >> orbit_linking_tmp" // TRIM(pidstr) // &
                   "." // TRIM(obs_type))
           END DO
           CALL NEW(obsfile, "orbit_linking_tmp" // &
                TRIM(pidstr) // "." // TRIM(obs_type))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2025)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setStatusOld(obsfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2030)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(obsfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2035)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (ANY(stdevs < 0.0_bp)) THEN
              CALL NEW(obss_all, obsfile)
           ELSE
              CALL NEW(obss_all, obsfile, stdev=stdevs)
           END IF
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2040)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obsfile)
           DEALLOCATE(id_arr_, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (2015).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ! Observations read.

           ! Extract different sets as well as a single observation
           ! and the epoch of each set:
           obss_sep_arr => getSeparatedSets(obss_all)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2045)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NEW(obss_selected)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2050)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obs)
           DO i=1,SIZE(obss_sep_arr)
              obs = getObservation(obss_sep_arr(i),1)
              CALL addObservation(obss_selected, obs)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2055)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              t = getTime(obs)
              mjds0(i) = getMJD(t, "TT")
              CALL NULLIFY(t)
              CALL NULLIFY(obs)
           END DO
           ! Extraction done.

           ! Check that observation sets from at least three different
           ! apparitions exist in the same address:
           DO i=1,SIZE(mjds0)-1
              k = 0
              DO j=i+1,SIZE(mjds0)
                 IF (ABS(mjds0(i)-mjds0(j)) < dt_min) THEN
                    k = k + 1
                 END IF
              END DO
              IF (nset - k < 3) THEN
                 CALL NULLIFY(obss_all)
                 CALL NULLIFY(obss_selected)
                 DO j=1,SIZE(obss_sep_arr)
                    CALL NULLIFY(obss_sep_arr(j))
                 END DO
                 DEALLOCATE(obss_sep_arr, mjds0, stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("orbit_linking", &
                         "Could not deallocate memory (2020).", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CYCLE addressloop
              END IF
           END DO
           ! Check done.

           obs = getObservation(obss_all,1)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2075)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF

           ! Compute orbital-element epoch:
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2080)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           mjd0 = getMJD(t, "TT")
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2085)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(t)
           CALL NULLIFY(obs)
           mjd0 = mjd0 + 0.5_bp*getObservationalTimespan(obss_all)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2090)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NEW(t0, mjd0, "TT")
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2095)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obss_all)
           ! Orbital-element epoch computed.

           ! Compute orbital-elements and their bounds
           elements = integerToArray(indx, elm, width, box_nrs, bounds, error)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (2070)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ! Move elements from mid-point of interval to the beginning:
           elements = elements - 0.5_bp*width
           delements = width
           ! For elements that were not included in binning, use the
           ! upper and lower limits for acceptable values:
           DO i=1,6
              IF (elm(i)) THEN
                 CYCLE
              END IF
              elements(i) = bounds(i,1)
              delements(i) = bounds(i,2) - bounds(i,1)
           END DO
           IF (info_verb >= 2) THEN
              IF (comparison_variable_type == "keplerian") THEN
                 WRITE(lu,"(6(1X,F10.5))") elements(1:2), elements(3:6)/rad_deg
              ELSE IF (comparison_variable_type == "modified keplerian") THEN
                 WRITE(lu,"(6(1X,F10.5))") elements(1:2), elements(3:5)/rad_deg, elements(6)
              ELSE
                 WRITE(lu,"(A)") "Not an option to write out elements of this type..."
              END IF
           END IF
           ! Orbital-elements and their bounds computed.

           ! Monte Carlo (MC) sampling:
           nlink_1 = 0
           CALL init(tree_2bmc)
           DO i=1,norb_mc
              ! Generate a random orbit within the computed bounds:
              CALL randomNumber(rans)
              elements_ = elements + rans*delements
              ! Change plane-crossing time to mean anomaly, if needed:
              IF (comparison_variable_type == "modified keplerian") THEN
                 elements_(6) = MODULO( &
                      SQRT(planetary_mu(11)/elements_(1)**3.0_bp) * &
                      (mjd0 - elements_(6)) - elements_(5), two_pi)
              END IF
              CALL NEW(orb_mc, elements_, "keplerian", "ecliptic", t0)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2105)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              ! random orbit generated.

              ! Compute residuals between the observed positions and
              ! the MC orbit:
              CALL NEW(storb, obss_selected)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2100)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              residuals => getResiduals(storb, orb_mc)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2110)", 1)
                 error = .FALSE.
                 CALL NULLIFY(storb)
                 CALL NULLIFY(orb_mc)
                 DEALLOCATE(residuals, stat=err)
                 CYCLE                 
              END IF
              ! Residuals computed.

              ! Check how many observation sets have acceptable
              ! residuals:
              nset = 0
              DO j=1,SIZE(residuals,dim=1)
                 IF (ALL(ABS(residuals(j,2:3)) < res_max_mc)) THEN
                    nset = nset + 1
                 END IF
              END DO
              IF (nset < 3) THEN
                 CALL NULLIFY(storb)
                 CALL NULLIFY(orb_mc)
                 DEALLOCATE(residuals, stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("orbit_linking", &
                         "Could not deallocate memory (2025).", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CYCLE
              END IF
              ALLOCATE(sets(nset), stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not allocate memory (2010).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              k = 0
              DO j=1,SIZE(residuals,dim=1)
                 IF (ALL(ABS(residuals(j,2:3)) < res_max_mc)) THEN
                    k = k + 1
                    sets(k) = j
                 END IF
              END DO
              ! Residual check performed.

              ! Scan linkages
              DO ii=1,nset-2
                 DO jj=ii+1,nset-1
                    DO kk=jj+1,nset
                       IF (ABS(mjds0(sets(ii))-mjds0(sets(jj))) < dt_min .OR. &
                            ABS(mjds0(sets(ii))-mjds0(sets(kk))) < dt_min .OR. &
                            ABS(mjds0(sets(jj))-mjds0(sets(kk))) < dt_min) THEN
                          CYCLE
                       END IF
                       id_arr3(1) = TRIM(getID(obss_sep_arr(sets(ii))))
                       id_arr3(2) = TRIM(getID(obss_sep_arr(sets(jj))))
                       id_arr3(3) = TRIM(getID(obss_sep_arr(sets(kk))))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (2120)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL insertionSort(id_arr3)
                       key3 = TRIM(id_arr3(1)) // " " // &
                            TRIM(id_arr3(2)) // " " // &
                            TRIM(id_arr3(3))
                       IF (simulated_observations) THEN
                          tree_node_ch32_1 => search(tree_c2blsl, key3)
                          tree_node_ch32_2 => search(tree_e2blsl, key3)
                          IF (.NOT.ASSOCIATED(tree_node_ch32_1,tree_c2blsl%nil) .OR. &
                               .NOT.ASSOCIATED(tree_node_ch32_2,tree_e2blsl%nil)) THEN
                             CYCLE
                          END IF
                          IF (getDesignation(obss_sep_arr(sets(ii))) == &
                               getDesignation(obss_sep_arr(sets(jj))) .AND. &
                               getDesignation(obss_sep_arr(sets(ii))) == &
                               getDesignation(obss_sep_arr(sets(kk)))) THEN
                             IF (error) THEN
                                CALL errorMessage("orbit_linking", &
                                     "TRACE BACK (2125)", 1)
                                CALL NULLIFY(logfile)
                                STOP
                             END IF
                             tree_node_ch32 => search(tree_c2bmc, key3)
                             IF (ASSOCIATED(tree_node_ch32,tree_c2bmc%nil)) THEN
                                CALL insert_tree_node(tree_c2bmc, key3)
                             END IF
                          ELSE
                             tree_node_ch32 => search(tree_e2bmc, key3)
                             IF (ASSOCIATED(tree_node_ch32,tree_e2bmc%nil)) THEN
                                CALL insert_tree_node(tree_e2bmc, key3)
                             END IF
                             ! Do not process erroneous 2bMC linkages further
                             !CYCLE
                          END IF
                       ELSE
                          tree_node_ch32 => search(tree_2blsl, key3)
                          IF (.NOT.ASSOCIATED(tree_node_ch32,tree_2blsl%nil)) THEN
                             CYCLE
                          END IF
                       END IF
                       rms_val = SQRT(SUM(residuals(sets(ii),2:3)**2) + &
                            SUM(residuals(sets(jj),2:3)**2) + &
                            SUM(residuals(sets(kk),2:3)**2)/6.0_bp)
                       rarr8 = (/ elements_, mjd0, rms_val /)
                       tree_node_ch32_8r8 => search(tree_2bmc, key3)
                       IF (ASSOCIATED(tree_node_ch32_8r8,tree_2bmc%nil)) THEN
                          i4arr = (/ sets(ii), sets(jj), sets(kk) /)
                          CALL insert_tree_node(tree_2bmc, key3, rarr8, i4arr)
                          nlink_1 = nlink_1 + 1
                       ELSE IF (tree_node_ch32_8r8%r8arr(8) > rms_val) THEN
                          tree_node_ch32_8r8%r8arr = rarr8
                       END IF
                    END DO
                 END DO
              END DO
              CALL NULLIFY(storb)
              CALL NULLIFY(orb_mc)
              DEALLOCATE(residuals, sets, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2030).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL NULLIFY(t0)
           CALL NULLIFY(obss_selected)
           DEALLOCATE(mjds0, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (2040).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ! Linkages scanned.

           ! Compute LSL:
           nlink_2 = 0
           i = 0
           tree_node_ch32_8r8_1 => minimum(tree_2bmc, tree_2bmc%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_8r8_1,tree_2bmc%nil))
              nlink_2 = nlink_2 + 1
              i = i + 1
              WRITE(lu,"(1X,I15,1(1X,I5),7(1X,F11.5))") naddress, i, &
                   tree_node_ch32_8r8_1%r8arr(1:2), &
                   tree_node_ch32_8r8_1%r8arr(3:6)/rad_deg, &
                   tree_node_ch32_8r8_1%r8arr(8)/rad_asec
              obss_ = obss_sep_arr(tree_node_ch32_8r8_1%i4arr(1)) + &
                   obss_sep_arr(tree_node_ch32_8r8_1%i4arr(2))
              obss = obss_ + obss_sep_arr(tree_node_ch32_8r8_1%i4arr(3))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2130)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NULLIFY(obss_)
              IF (info_verb >= 3) THEN
                 CALL writeObservationFile(obss, lu, TRIM(obs_type))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (2135)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
              END IF
              obs = getObservation(obss,1)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2140)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              t = getTime(obs)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2145)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              mjd0 = getMJD(t, "TT")
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2150)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NULLIFY(t)
              CALL NULLIFY(obs)
              mjd0 = mjd0 + 0.5_bp*getObservationalTimespan(obss)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2155)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NEW(t_inv, mjd0, "TT")
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2160)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NEW(storb, obss)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2165)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NULLIFY(obss)
              CALL setParameters(storb, &
                   dyn_model=dyn_model_step2, &
                   perturbers=perturbers, &
                   integrator=integrator, &
                   integration_step=integration_step, &
                   outlier_rejection=.FALSE., &
                   outlier_multiplier=3.0_bp, &
                   t_inv=t_inv, &
                   element_type="keplerian", &
                   multiple_objects=.TRUE., &
                   accept_multiplier=3.0_bp, &
                   ls_niter_major_max=15, &
                   ls_niter_major_min=2, &
                   ls_niter_minor=50)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2170)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NEW(t0, tree_node_ch32_8r8_1%r8arr(7), "TT")
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2095)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NEW(orb, tree_node_ch32_8r8_1%r8arr(1:6), "keplerian", "ecliptic", t0)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2175)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NULLIFY(t0)
              CALL setParameters(orb, &
                   dyn_model=dyn_model_step2, &
                   perturbers=perturbers, &
                   integrator=integrator, &
                   integration_step=integration_step)                 
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (2180)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL propagate(orb, t_inv)
              IF (error) THEN
                 error = .FALSE.
                 CALL NULLIFY(storb)
                 CALL NULLIFY(orb)
                 CYCLE
              END IF
              CALL NULLIFY(t_inv)
              ! First least squares incorporating three
              ! observation sets and using the 2-b model
              DO j=1,5
                 ls_element_mask = .TRUE.
                 ls_correction_factor=0.1_bp
                 SELECT CASE (j)
                 CASE (2)
                    ls_element_mask(2:6) = .FALSE.
                    ls_correction_factor=0.1_bp
                 CASE (3)
                    ls_element_mask(3:6) = .FALSE.
                    ls_correction_factor=0.1_bp
                 CASE (4)
                    ls_element_mask(3:4) = .FALSE.
                    ls_correction_factor=0.1_bp
                 END SELECT
                 CALL setParameters(storb, &
                      ls_element_mask=ls_element_mask, &
                      ls_correction_factor=ls_correction_factor)
                 IF (error .AND. j /= 5) THEN
                    error = .FALSE.
                    CYCLE
                 END IF
                 CALL leastSquares(storb, orb)
                 IF (error .AND. j /= 5) THEN
                    error = .FALSE.
                 ELSE IF (.NOT.error .AND. j > 1 .AND. j < 5) THEN
                    CALL NULLIFY(orb)
                    orb = getNominalOrbit(storb)
                 ELSE IF (.NOT.error .AND. j == 1) THEN
                    CALL NULLIFY(orb)
                    orb = getNominalOrbit(storb)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2182)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    rms = getRMS(storb, orb)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2183)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    IF (SQRT(SUM(rms(2:3)**2)) < rms_max_2blsl) THEN
                       EXIT
                    END IF
                 END IF
              END DO
              CALL NULLIFY(orb)
              IF (error) THEN
                 error = .FALSE.
              ELSE
                 orb = getNominalOrbit(storb)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (2185)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 rms = getRMS(storb, orb)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (2190)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF

                 IF (SQRT(SUM(rms(2:3)**2)) < rms_max_2blsl) THEN

                    ! 2BLSL FILTER PASSED

                    nlinks = nlinks + 1
                    CALL toString(nlinks, id1, error)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2195)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    DO WHILE (LEN_TRIM(id1) < 7)
                       id1 = "0" // TRIM(id1)
                    END DO
                    IF (simulated_observations) THEN
                       IF (getDesignation(obss_sep_arr(tree_node_ch32_8r8_1%i4arr(1))) == &
                            getDesignation(obss_sep_arr(tree_node_ch32_8r8_1%i4arr(2))) .AND. &
                            getDesignation(obss_sep_arr(tree_node_ch32_8r8_1%i4arr(1))) == &
                            getDesignation(obss_sep_arr(tree_node_ch32_8r8_1%i4arr(3)))) THEN
                          tree_node_ch32 => search(tree_c2blsl, tree_node_ch32_8r8_1%key)
                          IF (ASSOCIATED(tree_node_ch32,tree_c2blsl%nil)) THEN
                             CALL insert_tree_node(tree_c2blsl, tree_node_ch32_8r8_1%key)
                          END IF
                       ELSE
                          tree_node_ch32 => search(tree_e2blsl, tree_node_ch32_8r8_1%key)
                          IF (ASSOCIATED(tree_node_ch32,tree_e2blsl%nil)) THEN
                             CALL insert_tree_node(tree_e2blsl, tree_node_ch32_8r8_1%key)
                          END IF
                       END IF
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (2200)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                    ELSE
                       tree_node_ch32 => search(tree_2blsl, tree_node_ch32_8r8_1%key)
                       IF (ASSOCIATED(tree_node_ch32,tree_2blsl%nil)) THEN
                          CALL insert_tree_node(tree_2blsl, tree_node_ch32_8r8_1%key)
                       END IF
                    END IF
                    fname = get_cl_option("--2blnks=", "2blsl.lnks")
                    CALL NEW(tmpfile, TRIM(fname))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (3005)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL setPositionAppend(tmpfile)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (3015)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL OPEN(tmpfile)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (3015)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    WRITE(getUnit(tmpfile),"(A,1X,F0.2,1X,A)") &
                         TRIM(id1), &
                         SQRT(SUM(rms(2:3)**2))/rad_asec, &
                         TRIM(tree_node_ch32_8r8_1%key)
                    CALL NULLIFY(tmpfile)                                
                    fname = get_cl_option("--2borb=", "2b_rms_acc.orb")
                    CALL NEW(tmpfile,TRIM(fname))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2225)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL setPositionAppend(tmpfile)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2230)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL OPEN(tmpfile)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2235)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL writeOpenOrbOrbitFile(getUnit(tmpfile), insert_header, &
                         "keplerian", TRIM(id1), orb, "keplerian", &
                         getCovarianceMatrix(storb, "keplerian"))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (2240)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    IF (insert_header) THEN
                       insert_header = .FALSE.
                    END IF
                    CALL NULLIFY(tmpfile)
                    IF (info_verb >= 2) THEN
                       obss = getObservations(storb)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (2245)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL toString(nlinks, str, error)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (2246)", 1)
                          CALL NULLIFY(logfile)
                          STOP                          
                       END IF
                       CALL writeObservationFile(obss, lu, TRIM(obs_type), number=TRIM(str))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (2247)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL NULLIFY(obss)
                    END IF
                 END IF
                 CALL NULLIFY(orb)
              END IF
              CALL NULLIFY(storb)
              CALL NULLIFY(orb)
              tree_node_ch32_8r8_2 => tree_node_ch32_8r8_1
              tree_node_ch32_8r8_1 => successor(tree_2bmc, tree_node_ch32_8r8_1)
              tree_node_ch32_8r8_2 => delete_tree_node(tree_2bmc, tree_node_ch32_8r8_2)
              DEALLOCATE(tree_node_ch32_8r8_2)
           END DO
           CALL delete_tree(tree_2bmc)
           DO i=1,SIZE(obss_sep_arr)
              CALL NULLIFY(obss_sep_arr(i))
           END DO
           DEALLOCATE(obss_sep_arr, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (2040).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (nlink_1 /= nlink_2) THEN
              WRITE(*,*) " nlink1, nlink2: ", nlink_1, nlink_2
              STOP
           END IF

           IF (simulated_observations) THEN
              j = 0
              tree_node_ch32 => minimum(tree_c2bmc, tree_c2bmc%root)
              DO WHILE (.NOT.ASSOCIATED(tree_node_ch32,tree_c2bmc%nil))
                 tree_node_ch32 => successor(tree_c2bmc, tree_node_ch32)
                 j = j + 1
              END DO
              WRITE(lu,"(1X,A,1X,I0)") "Number of correct 3-links found in 2BMC:", j
              j = 0
              tree_node_ch32 => minimum(tree_e2bmc, tree_e2bmc%root)
              DO WHILE (.NOT.ASSOCIATED(tree_node_ch32,tree_e2bmc%nil))
                 tree_node_ch32 => successor(tree_e2bmc, tree_node_ch32)
                 j = j + 1
              END DO
              WRITE(lu,"(1X,A,1X,I0)") "Number of erroneous 3-links found in 2BMC:", j
              j = 0
              tree_node_ch32 => minimum(tree_c2blsl, tree_c2blsl%root)
              DO WHILE (.NOT.ASSOCIATED(tree_node_ch32,tree_c2blsl%nil))
                 tree_node_ch32 => successor(tree_c2blsl, tree_node_ch32)
                 j = j + 1
              END DO
              WRITE(lu,"(1X,A,1X,I0)") "Number of correct 3-links found in 2BLSL:", j
              j = 0
              tree_node_ch32 => minimum(tree_e2blsl, tree_e2blsl%root)
              DO WHILE (.NOT.ASSOCIATED(tree_node_ch32,tree_e2blsl%nil))
                 tree_node_ch32 => successor(tree_e2blsl, tree_node_ch32)
                 j = j + 1
              END DO
              WRITE(lu,"(1X,A,1X,I0)") "Number of erroneous 3-links found in 2BLSL:", j
           ELSE
              WRITE(lu,"(1X,A,1X,I0)") "Number of 3-links found in 2BLSL:", nlinks
           END IF
        END DO addressloop
        CALL NULLIFY(addressfile)
        ! Delete and deallocate trees:
        IF (simulated_observations) THEN
           tree_node_ch32_1 => minimum(tree_c2bmc, tree_c2bmc%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, tree_c2bmc%nil))
              tree_node_ch32_2 => tree_node_ch32_1
              tree_node_ch32_1 => successor(tree_c2bmc, tree_node_ch32_1)
              tree_node_ch32_2 => delete_tree_node(tree_c2bmc, tree_node_ch32_2)
              DEALLOCATE(tree_node_ch32_2, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2045).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL delete_tree(tree_c2bmc)
           tree_node_ch32_1 => minimum(tree_e2bmc, tree_e2bmc%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, tree_e2bmc%nil))
              tree_node_ch32_2 => tree_node_ch32_1
              tree_node_ch32_1 => successor(tree_e2bmc, tree_node_ch32_1)
              tree_node_ch32_2 => delete_tree_node(tree_e2bmc, tree_node_ch32_2)
              DEALLOCATE(tree_node_ch32_2, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2050).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL delete_tree(tree_e2bmc)
           tree_node_ch32_1 => minimum(tree_c2blsl, tree_c2blsl%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, tree_c2blsl%nil))
              tree_node_ch32_2 => tree_node_ch32_1
              tree_node_ch32_1 => successor(tree_c2blsl, tree_node_ch32_1)
              tree_node_ch32_2 => delete_tree_node(tree_c2blsl, tree_node_ch32_2)
              DEALLOCATE(tree_node_ch32_2, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2055).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL delete_tree(tree_c2blsl)
           tree_node_ch32_1 => minimum(tree_e2blsl, tree_e2blsl%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, tree_e2blsl%nil))
              tree_node_ch32_2 => tree_node_ch32_1
              tree_node_ch32_1 => successor(tree_e2blsl, tree_node_ch32_1)
              tree_node_ch32_2 => delete_tree_node(tree_e2blsl, tree_node_ch32_2)
              DEALLOCATE(tree_node_ch32_2, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2060).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL delete_tree(tree_e2blsl)
        ELSE
           tree_node_ch32_1 => minimum(tree_2blsl, tree_2blsl%root)
           DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, tree_2blsl%nil))
              tree_node_ch32_2 => tree_node_ch32_1
              tree_node_ch32_1 => successor(tree_2blsl, tree_node_ch32_1)
              tree_node_ch32_2 => delete_tree_node(tree_2blsl, tree_node_ch32_2)
              DEALLOCATE(tree_node_ch32_2, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory (2070).", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END DO
           CALL delete_tree(tree_2blsl)
        END IF
        ret = system("rm -f orbit_linking_tmp" // TRIM(pidstr) // &
             "." // TRIM(obs_type))


     CASE (3)

        !! PHASE 3: FIND COMMON NB ORBIT

        nnolinks = 0
        nclink = 0
        nelink = 0
        ! Obtain name of observation file
        obsfname = get_cl_option("--obsall=", " ")
        IF (LEN_TRIM(obsfname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Observation file missing.", 1)
           CALL NULLIFY(logfile)
           STOP     
        END IF
        ! Read input 2b orbits
        fname = get_cl_option("--2borb=", " ")
        IF (LEN_TRIM(fname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Input orbits missing.", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NEW(orbfile, TRIM(fname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3005)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(orbfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(orbfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        norb = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3020)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        ALLOCATE(id_arr_prm(norb), orb_arr(norb), stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("orbit_linking", &
                "Could not allocate memory (3005).", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        id_arr_prm = " "
        header(1:4)(:) = " "
        DO i=1,norb
           CALL readOpenOrbOrbitFile(getUnit(orbfile), header, &
                element_type_in=element_type_in_prm, &
                id=id_arr_prm(i), orb=orb_arr(i))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not read orbit file.", 1)
              IF (id_arr_prm(i)(1:1) == "#") THEN
                 CALL errorMessage("orbit_linking", &
                      "Check orbit file for extra headers.", 1)
              END IF
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (element_type_nblsl == "keplerian") THEN
              CALL toKeplerian(orb_arr(i))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3025)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END IF
        END DO
        CALL NULLIFY(orbfile)
        ! Find out how many linkages have already been found:
        fname = get_cl_option("--nborb=", "nb_rms_acc.orb")
        CALL NEW(orbfile, TRIM(fname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3030)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(orbfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3035)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        nlinks = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3040)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        IF (nlinks < 0) THEN
           nlinks = 0
        END IF
        CALL NULLIFY(orbfile)
        lnksfname = get_cl_option("--2blnks="," ")
        IF (LEN_TRIM(lnksfname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Linkages file missing.", 1)
           CALL NULLIFY(logfile)
           STOP     
        END IF
        CALL NEW(lnksfile, TRIM(lnksfname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3045)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(lnksfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3050)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(lnksfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (3055)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        loop_start = get_cl_option("--loop-start=", 1)
        loop_stop = get_cl_option("--loop-stop=", norb)
        DO i=1,loop_stop
           !READ(getUnit(lnksfile),"(A,1X,A,1X,A)",iostat=err) id_arr3
           READ(getUnit(lnksfile),*,iostat=err) id1, rms_val, id_arr3
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not read linkages file.", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF

           ! 
           IF (TRIM(id1) /= TRIM(id_arr_prm(i))) THEN
              CALL errorMessage("orbit_linking", &
                   "IDs in 2blnks file and 2b orbit file differ.", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF

           ! Skip candidate 2b links if starting from later
           IF (i < loop_start) THEN
              CYCLE
           END IF
           WRITE(lu,"(1X,A,I0,2A)") "Analyzing 2-b orbit #", i, ": ", TRIM(id_arr_prm(i))

           ! Skip 2b lsl solutions that do not meet the rms criteria
           ! for this filter:
           IF (rms_val < rms_2blsl_step3_min .OR. &
                rms_val > rms_2blsl_step3_max) THEN
              CYCLE
           END IF

           ! Extract observations:
           ret = system("rm -f orbit_linking_tmp" // &
                TRIM(pidstr) // "." // TRIM(obs_type))
           DO j=1,3
              ret = system("grep " // TRIM(id_arr3(j)) // &
                   " " // TRIM(obsfname) // &
                   " >> orbit_linking_tmp" // TRIM(pidstr) // &
                   "." // TRIM(obs_type))
           END DO
           CALL NEW(obsfile, "orbit_linking_tmp" // &
                TRIM(pidstr) // "." // TRIM(obs_type))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3060)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setStatusOld(obsfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3065)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(obsfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3070)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (ANY(stdevs < 0.0_bp)) THEN
              CALL NEW(obss, obsfile)
           ELSE
              CALL NEW(obss, obsfile, stdev=stdevs)
           END IF
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3075)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obsfile)
           IF (simulated_observations) THEN
              CALL setNumber(obss, 0)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3080)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              obs_set_ids => getObjects(obss)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3085)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              nobj = SIZE(obs_set_ids)
              DEALLOCATE(obs_set_ids, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "Could not deallocate memory.", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END IF
           CALL toInt(TRIM(id_arr_prm(i)), nr, error)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3090)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setNumber(obss, nr)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3095)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ! Observations extracted.

           key3 = TRIM(id_arr3(1)) // " " // &
                TRIM(id_arr3(2)) // " " // &
                TRIM(id_arr3(3))
           CALL NEW(storb, obss)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3100)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setParameters(storb, &
                dyn_model=dyn_model_step3, &
                perturbers=perturbers, &
                integrator=integrator, &
                integration_step=integration_step, &
                outlier_rejection=.FALSE., &
                outlier_multiplier=3.0_bp, &
                element_type="keplerian", &
                multiple_objects=.TRUE., &
                accept_multiplier=3.0_bp, &
                ls_niter_major_max=10, &
                ls_niter_major_min=2, &
                ls_niter_minor=20)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3105)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           orb = copy(orb_arr(i))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3110)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setParameters(orb, &
                dyn_model=dyn_model_step3, &
                perturbers=perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (3115)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           DO j=1,3
              ls_element_mask = .TRUE.
              ls_correction_factor = 1.0_bp
              SELECT CASE (j)
              CASE (1)
                 ls_niter_major_max =  10
                 ls_niter_major_min =   2
                 ls_niter_minor =      20
              CASE (2)
                 ls_element_mask(3:6) = .FALSE.
                 ls_correction_factor = 0.1_bp
                 ls_niter_major_max =   5
                 ls_niter_major_min =   2
                 ls_niter_minor =      10
              CASE (3)
                 ls_niter_major_max =   5
                 ls_niter_major_min =   2
                 ls_niter_minor =      10                 
              END SELECT
              CALL setParameters(storb, &
                   ls_element_mask=ls_element_mask, &
                   ls_correction_factor=ls_correction_factor, &
                   ls_niter_major_max=ls_niter_major_max, &
                   ls_niter_major_min=ls_niter_major_min, &
                   ls_niter_minor=ls_niter_minor)
              IF (error .AND. j /= 3) THEN
                 error = .FALSE.
                 CYCLE
              END IF
              ! Second least squares incorporating three
              ! observation sets and using the full n-b model
              info_verb_ = info_verb
              info_verb = 3
              CALL leastSquares(storb, orb)
              info_verb = info_verb_
              IF (.NOT.error .AND. j == 1) THEN
                 EXIT
              ELSE IF (error .AND. j == 1) THEN
                 error = .FALSE.
              ELSE IF (j == 2) THEN
                 error = .FALSE.
                 CALL NULLIFY(orb)
                 orb = getNominalOrbit(storb)                 
              END IF
           END DO
           CALL NULLIFY(orb)
           IF (error) THEN
              error = .FALSE.
              CALL NEW(tmpfile,"nb_inv_failed.lnks")
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3120)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL setPositionAppend(tmpfile)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3125)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL OPEN(tmpfile)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3130)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              WRITE(getUnit(tmpfile),"(A)",iostat=err) TRIM(key3)
              IF (err /= 0) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3135)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL NULLIFY(tmpfile)
           ELSE
              orb = getNominalOrbit(storb)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3140)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              rms = getRMS(storb, orb)
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (3145)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              IF (SQRT(SUM(rms(2:3)**2)) < rms_max_nblsl) THEN
                 ! RMS ACCEPTABLE
                 nlinks = nlinks + 1
                 IF (simulated_observations) THEN
                    IF (nobj == 1) THEN
                       nclink = nclink + 1
                    ELSE
                       nelink = nelink + 1
                    END IF
                 END IF
                 WRITE(lu,"(1X,2(A,1X),2(F15.5,1X))") &
                      "Possible 3-link:", TRIM(key3), rms(2:3)/rad_asec

                 ! WRITE OBSERVATION FILE
                 fname = get_cl_option("--nbobs=", "nb_rms_acc." // TRIM(obs_type))
                 CALL NEW(tmpfile,TRIM(fname))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3150)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setPositionAppend(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3155)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL OPEN(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3160)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL writeObservationFile(obss, getUnit(tmpfile), TRIM(obs_type))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3165)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 WRITE(getUnit(tmpfile),*)
                 CALL NULLIFY(tmpfile)

                 ! WRITE ORBIT FILE
                 fname = get_cl_option("--nborb=", "nb_rms_acc.orb")
                 CALL NEW(tmpfile,TRIM(fname))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3170)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setPositionAppend(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3175)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL OPEN(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3180)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL writeOpenOrbOrbitFile(getUnit(tmpfile), nlinks==1, &
                      "keplerian", TRIM(id1), orb, "keplerian", &
                      getCovarianceMatrix(storb, "keplerian"))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3185)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL NULLIFY(tmpfile)

                 ! WRITE LSL FILE INCL. ASTROMETRY, RESIDUALS,
                 ! ELEMENTS, UNCERTAINTIES, CORRELATIONS:
                 fname = get_cl_option("--nblsl=", "nb_rms_acc.lsl")
                 CALL NEW(tmpfile,TRIM(fname))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3190)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setPositionAppend(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3195)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL OPEN(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3200)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL writeObservationFile(obss, getUnit(tmpfile), &
                      TRIM(obs_type)) 
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3205)", 1)
                    WRITE(getUnit(tmpfile),"(A)") &
                         "Could not write observations for object " // TRIM(id1) 
                    error = .FALSE.
                 END IF
                 WRITE(getUnit(tmpfile),"(A)") "#"
                 CALL writeNominalSolution(storb, obss, "keplerian", getUnit(tmpfile))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3210)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL NULLIFY(tmpfile)

              ELSE
                 ! RMS TOO LARGE
                 nnolinks = nnolinks + 1
                 CALL NEW(tmpfile,"nb_rms_failed.lnks")
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3215)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setPositionAppend(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3220)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL OPEN(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3225)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 WRITE(getUnit(tmpfile),"(A)",iostat=err) TRIM(key3)
                 IF (err /= 0) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3230)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL NULLIFY(tmpfile)
                 CALL NEW(tmpfile,"nb_rms_failed.orb")
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3235)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setPositionAppend(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3240)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL OPEN(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3245)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL writeopenOrbOrbitFile(getUnit(tmpfile), nnolinks==1, &
                      "keplerian", TRIM(id1), orb, "keplerian", &
                      getCovarianceMatrix(storb, "keplerian"))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (3250)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL NULLIFY(tmpfile)
              END IF
           END IF
           CALL NULLIFY(obss)
           CALL NULLIFY(storb)
           CALL NULLIFY(orb)
           IF (simulated_observations) THEN
              WRITE(lu,"(1X,A,1X,I0)") "Number of correct 3-links found in NBLSL:", nclink
              WRITE(lu,"(1X,A,1X,I0)") "Number of erroneous 3-links found in NBLSL:", nelink
           ELSE
              WRITE(lu,"(1X,A,1X,I0)") "Number of 3-links found in NBLSL:", nlinks
           END IF
           CALL NULLIFY(orb_arr(i))
        END DO
        CALL NULLIFY(lnksfile)
        DEALLOCATE(id_arr_prm, orb_arr)





     CASE (4)

        ! PHASE 4: FIND ADDITIONAL LINKAGES

        ls_element_mask = .TRUE.
        nlinks = 0
        nelink = 0
        ! Read input astrometry
        obsfname = get_cl_option("--nbobs=", " ")
        IF (LEN_TRIM(obsfname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Observation file for accepted n-body linkages missing.", 1)
           CALL NULLIFY(logfile)
           STOP           
        END IF
        CALL NEW(obsfile, TRIM(obsfname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4005)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(obsfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(obsfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NEW(obss_nb, obsfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4020)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NULLIFY(obsfile)
        obss_sep_arr => getSeparatedSets(obss_nb)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4025)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NULLIFY(obss_nb)
        ! Read input nb orbits
        fname = get_cl_option("--nborb=", " ")
        IF (LEN_TRIM(fname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Orbit file for accepted n-body linkages missing.", 1)
           CALL NULLIFY(logfile)
           STOP           
        END IF
        CALL NEW(orbfile, TRIM(fname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4030)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(orbfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4035)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(orbfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4040)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        norb = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4045)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        ALLOCATE(id_arr_prm(norb), orb_arr(norb), stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("orbit_linking", &
                "Could not allocate memory (20).", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        id_arr_prm = " "
        header(1:4)(:) = " "
        DO i=1,norb
           CALL readOpenOrbOrbitFile(getUnit(orbfile), header, &
                element_type_in=element_type_in_prm, &
                id=id_arr_prm(i), orb=orb_arr(i))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not read orbit file.", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (element_type_nblsl == "keplerian") THEN
              CALL toKeplerian(orb_arr(i))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (4050)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END IF
        END DO
        CALL NULLIFY(orbfile)
        obsallfname = get_cl_option("--obsall=", " ")
        IF (LEN_TRIM(obsallfname) == 0) THEN
           CALL errorMessage("orbit_linking", &
                "Complete observation file missing.", 1)
           CALL NULLIFY(logfile)
           STOP           
        END IF
        fname = get_cl_option("--obsfirst=", " ")
        IF (LEN_TRIM(fname) /= 0) THEN
           CALL CPU_TIME(dt1)
           CALL NEW(obsfile, TRIM(fname))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4055)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setStatusOld(obsfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4060)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           WRITE(lu,"(1X,2A)") "Open observation file ", TRIM(fname)
           CALL OPEN(obsfile)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4065)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           WRITE(lu,"(1X,2A)") "Read observation file ", TRIM(fname)
           IF (ANY(stdevs < 0.0_bp)) THEN
              CALL NEW(obss_all_first, obsfile)
           ELSE
              CALL NEW(obss_all_first, obsfile, stdev=stdevs)
           END IF
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4070)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obsfile)
           WRITE(lu,"(1X,3A)") "Observation file ", TRIM(fname), " read"
           CALL CPU_TIME(dt2)
           WRITE(lu,*) 'Time used for reading input astrometry: ', TRIM(secToHMS(dt2-dt1, error))
        ELSE
           CALL errorMessage("orbit_linking", &
                "File containing the first observation of each set is missing.", 1)
           CALL NULLIFY(logfile)
           STOP     
        END IF
        obs_arr_all_first => getObservations(obss_all_first)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (4075)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        DO i=1,SIZE(obss_sep_arr)
           CALL NULLIFY(orb_init)
           id1 = getID(obss_sep_arr(i))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4080)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           WRITE(lu,"(1X,A,I0,2(1X,A))") "Analyzing linkage #", i, ":", TRIM(id1)
           IF (id1 == id_arr_prm(i)) THEN
              orb_init = copy(orb_arr(i))
           ELSE
              DO j=1,norb
                 IF (id1 == id_arr_prm(j)) THEN
                    orb_init = copy(orb_arr(j))
                    EXIT
                 END IF
              END DO
           END IF
           CALL setParameters(orb_init, &
                dyn_model=dyn_model_step3, &
                perturbers=perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4085)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NEW(storb, obss_all_first)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4090)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setParameters(storb, &
                dyn_model=dyn_model_step3, &
                perturbers=perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4095)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(orb)
           orb = copy(orb_init)
           residuals => getResiduals(storb, orb)
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4100)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (info_verb >= 3) THEN
              CALL writeObservationFile(obss_all_first, lu, TRIM(obs_type))
              IF (error) THEN
                 CALL errorMessage("orbit_linking", &
                      "TRACE BACK (4105)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
           END IF
           CALL NULLIFY(storb)
           CALL NULLIFY(orb)
           id_arr => getDesignations(obss_sep_arr(i))
           IF (error) THEN
              CALL errorMessage("orbit_linking", &
                   "TRACE BACK (4110)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           DO j=1,SIZE(residuals,dim=1)
              IF (ALL(ABS(residuals(j,2:3)) < res_max_addlink)) THEN
                 id2 = getID(obs_arr_all_first(j))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4115)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 DO WHILE (LEN_TRIM(id2) < 7)
                    id2 = "0" // TRIM(id2)
                 END DO
                 IF (ANY(id_arr == TRIM(id2))) THEN
                    CYCLE
                 END IF
                 IF (info_verb >= 1) THEN
                    WRITE(lu,"(1X,A,2(1X,F20.5))") TRIM(id2), residuals(j,2:3)/rad_amin
                 END IF
                 ret = system("grep " // TRIM(id2) // &
                      " " // TRIM(obsallfname) // &
                      " > orbit_linking_tmp" // TRIM(pidstr) // &
                      "." // TRIM(obs_type))
                 CALL NEW(obsfile, "orbit_linking_tmp" // &
                      TRIM(pidstr) // "." // TRIM(obs_type))
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4120)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setStatusOld(obsfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4125)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL OPEN(obsfile)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4130)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 IF (ANY(stdevs < 0.0_bp)) THEN
                    CALL NEW(obss_2, obsfile)
                 ELSE
                    CALL NEW(obss_2, obsfile, stdev=stdevs)
                 END IF
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4135)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL NULLIFY(obsfile)
                 IF (getNrOfObservations(obss_2) == 0) THEN
                    CALL errorMessage("orbit_linking", &
                         "Additional observations not found in the " // &
                         "complete astrometry file.", 1)
                    CALL NULLIFY(logfile)
                    STOP                    
                 END IF
                 obss_1 = obss_sep_arr(i) + obss_2
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4140)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 IF (info_verb >= 1) THEN
                    CALL writeObservationFile(obss_1, lu, TRIM(obs_type))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4145)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                 END IF
                 CALL NEW(storb, obss_1)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4150)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 CALL setParameters(storb, &
                      dyn_model=dyn_model_step3, &
                      perturbers=perturbers, &
                      integrator=integrator, &
                      integration_step=integration_step, &
                      outlier_rejection=.FALSE., &
                      outlier_multiplier=3.0_bp, &
                      element_type="keplerian", &
                      multiple_objects=.TRUE., &
                      accept_multiplier=3.0_bp, &
                      ls_element_mask=ls_element_mask, &
                      ls_correction_factor=1.0_bp, &
                      ls_niter_major_max=10, &
                      ls_niter_major_min=2, &
                      ls_niter_minor=10)
                 IF (error) THEN
                    CALL errorMessage("orbit_linking", &
                         "TRACE BACK (4155)", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
                 ! Third least squares incorporating three
                 ! observation sets and using the full n-b model
                 orb = copy(orb_init)
                 info_verb_ = info_verb
                 info_verb = 3
                 CALL leastSquares(storb, orb)
                 info_verb = info_verb_
                 CALL NULLIFY(orb)
                 IF (error) THEN
                    error = .FALSE.
                    CALL NEW(tmpfile,"add_nb_inv_failed." // TRIM(obs_type))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4160)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL setPositionAppend(tmpfile)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4165)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL OPEN(tmpfile)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4170)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    CALL writeObservationFile(obss_1, getUnit(tmpfile), TRIM(obs_type))
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4175)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    WRITE(getUnit(tmpfile),*)
                    CALL NULLIFY(tmpfile)
                 ELSE
                    orb = getNominalOrbit(storb)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4180)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    rms = getRMS(storb, orb)
                    IF (error) THEN
                       CALL errorMessage("orbit_linking", &
                            "TRACE BACK (4185)", 1)
                       CALL NULLIFY(logfile)
                       STOP
                    END IF
                    IF (SQRT(SUM(rms(2:3)**2)) < 3.0_bp*rad_asec) THEN

                       nlinks = nlinks + 1
                       CALL setNumber(obss_1, nlinks)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4188)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL toString(nlinks, id1, error)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4190)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       DO WHILE (LEN_TRIM(id1) < 7)
                          id1 = "0" // TRIM(id1)
                       END DO
                       WRITE(lu,"(3(1X,A,1X),2(F15.5,1X))") &
                            "Possible link:", &
                            TRIM(id1), &
                            TRIM(id2), &
                            rms(2:3)/rad_asec

                       ! WRITE OBSERVATION FILE
                       fname = get_cl_option("--addnbobs=", "add_nb_rms_acc." // TRIM(obs_type))
                       CALL NEW(tmpfile,TRIM(fname))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4195)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL setPositionAppend(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4200)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL OPEN(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4205)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL writeObservationFile(obss_1, &
                            getUnit(tmpfile), TRIM(obs_type))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4210)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       WRITE(getUnit(tmpfile),*)
                       CALL NULLIFY(tmpfile)

                       ! WRITE ORBIT FILE
                       fname = get_cl_option("--addnborb=", "add_nb_rms_acc.orb")
                       CALL NEW(tmpfile,TRIM(fname))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4215)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL setPositionAppend(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4220)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL OPEN(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4225)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL writeOpenOrbOrbitFile(getUnit(tmpfile), nlinks==1, &
                            "keplerian", TRIM(id1), orb, "keplerian", &
                            getCovarianceMatrix(storb, "keplerian"))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4230)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL NULLIFY(tmpfile)

                       ! WRITE LSL FILE INCL. ASTROMETRY, RESIDUALS,
                       ! ELEMENTS, UNCERTAINTIES, CORRELATIONS:
                       fname = get_cl_option("--addnblsl=", "add_nb_rms_acc.lsl")
                       CALL NEW(tmpfile,TRIM(fname))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4235)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL setPositionAppend(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4240)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL OPEN(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4245)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL writeObservationFile(obss_1, getUnit(tmpfile), &
                            TRIM(obs_type))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4250)", 1)
                          WRITE(getUnit(tmpfile),"(A)") &
                               "Could not write observations for object " // TRIM(id1) 
                          error = .FALSE.
                       END IF
                       WRITE(getUnit(tmpfile),"(A)") "#"
                       CALL writeNominalSolution(storb, obss_1, "cartesian", getUnit(tmpfile))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4255)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL NULLIFY(tmpfile)

                    ELSE

                       nelink = nelink + 1
                       CALL setNumber(obss_1, getNumber(obss_sep_arr(i)))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4193)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF

                       ! WRITE OBSERVATION FILE
                       fname = get_cl_option("--faddnbobs=", "add_nb_rms_failed." // TRIM(obs_type))
                       CALL NEW(tmpfile,TRIM(fname))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4195)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL setPositionAppend(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4200)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL OPEN(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4205)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL writeObservationFile(obss_1, &
                            getUnit(tmpfile), TRIM(obs_type))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4210)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       WRITE(getUnit(tmpfile),*)
                       CALL NULLIFY(tmpfile)

                       ! WRITE ORBIT FILE
                       fname = get_cl_option("--faddnborb=", "add_nb_rms_failed.orb")
                       CALL NEW(tmpfile,TRIM(fname))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4215)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL setPositionAppend(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4220)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL OPEN(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4225)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL writeOpenOrbOrbitFile(getUnit(tmpfile), nelink==1, &
                            "keplerian", TRIM(id1), orb, "keplerian", &
                            getCovarianceMatrix(storb, "keplerian"))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4230)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL NULLIFY(tmpfile)

                       ! WRITE LSL FILE INCL. ASTROMETRY, RESIDUALS,
                       ! ELEMENTS, UNCERTAINTIES, CORRELATIONS:
                       fname = get_cl_option("--faddnblsl=", "add_nb_rms_failed.lsl")
                       CALL NEW(tmpfile,TRIM(fname))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4235)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL setPositionAppend(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4240)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL OPEN(tmpfile)
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4245)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL writeObservationFile(obss_1, getUnit(tmpfile), &
                            TRIM(obs_type)) 
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4250)", 1)
                          WRITE(getUnit(tmpfile),"(A)") &
                               "Could not write observations for object " // TRIM(id1) 
                          error = .FALSE.
                       END IF
                       WRITE(getUnit(tmpfile),"(A)") "#"
                       CALL writeNominalSolution(storb, obss_1, "cartesian", getUnit(tmpfile))
                       IF (error) THEN
                          CALL errorMessage("orbit_linking", &
                               "TRACE BACK (4255)", 1)
                          CALL NULLIFY(logfile)
                          STOP
                       END IF
                       CALL NULLIFY(tmpfile)

                    END IF
                 END IF
                 CALL NULLIFY(obss_2)
              END IF
              CALL NULLIFY(orb)
              CALL NULLIFY(storb)
              CALL NULLIFY(obss_1)
           END DO
           DEALLOCATE(residuals, id_arr, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("orbit_linking", &
                   "Could not deallocate memory (65).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obss_sep_arr(i))
           CALL NULLIFY(orb_arr(i))
        END DO
        DO i=1,SIZE(obs_arr_all_first)
           CALL NULLIFY(obs_arr_all_first(i))
        END DO
        DEALLOCATE(obss_sep_arr, orb_arr, id_arr_prm, &
             obs_arr_all_first)
        ret = system("rm -f orbit_linking_tmp" // TRIM(pidstr) // &
             "." // TRIM(obs_type))

     CASE (5)

        ! PHASE 5: FIND 2-LINKAGES BY SPLITTING 3-LINKAGES
        !          UNACCEPTABLE IN TERMS OF THEIR RMS

        ls_element_mask = .TRUE.
        nlinks = 0
        ! Read input astrometry
        obsfname = get_cl_option("--nbobs=", " ")
        IF (LEN_TRIM(obsfname) == 0) THEN
           obsfname = "nb_rms_acc." // TRIM(obs_type)
        END IF
        CALL NEW(obsfile, TRIM(obsfname))
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (5005)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(obsfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (5010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(obsfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (5015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NEW(obss_nb, obsfile)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (5020)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NULLIFY(obsfile)
        obss_sep_arr => getSeparatedSets(obss_nb)
        IF (error) THEN
           CALL errorMessage("orbit_linking", &
                "TRACE BACK (5025)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NULLIFY(obss_nb)
        DO i=1,SIZE(obss_sep_arr)
           CALL setNumber(obss_sep_arr(i),0)
           obss_sep_arr_ => getSeparatedSets(obss_sep_arr(i))
           k = 0
           DO j=1,SIZE(obss_sep_arr_)
              IF (getObservationalTimespan(obss_sep_arr_(j)) > 0.5_bp) THEN
                 k = k + 1
              END IF
           END DO
           IF (k <= 1) THEN

           END IF
        END DO


     END SELECT

  END DO

  ! Close logfile if used:
  fname = get_cl_option("--logfile=", " ")
  IF (LEN_TRIM(fname) /= 0) THEN
     CALL NULLIFY(logfile)
  END IF

END PROGRAM orbit_linking
