!====================================================================!
!                                                                    !
! Copyright 2002-2017,2018                                           !
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
!!   Asteroid identification software for linking within apparitions.
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
!!   @version 2018-12-10
!!
PROGRAM ephemeris_linking

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
  TYPE (rb_tree_r16_i4arr), POINTER :: &
       address_tree
  TYPE (rb_tree_ch32), POINTER :: &
       trial_linkage_tree
  TYPE (rb_tree_node_r16_i4arr), POINTER :: &
       address_tree_node1, &
       address_tree_node2
  TYPE (rb_tree_node_ch32), POINTER :: &
       tree_node_ch32, &
       tree_node_ch32_1, &
       tree_node_ch32_2
  TYPE (rb_tree_node_ch32_8r8), POINTER :: &
       tree_node_ch32_8r8, &
       tree_node_ch32_8r8_1, &
       tree_node_ch32_8r8_2
  TYPE (File) :: &
       optfile, &
       obsfile, &
       obsfile_, &
       orbidfile, &
       orbfile, &
       addressfile, &
       triallinkfile, &
       apriorifile, &
       logfile, &
       outfile, &
       resfile, &
       tmpfile
  TYPE (Time), DIMENSION(:), ALLOCATABLE :: &
       epoch_arr
  TYPE (Time) :: &
       epoch, &
       t0, &
       t, &
       t_inv
  TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
       scoord_arr
  TYPE (SphericalCoordinates) :: &
       scoord
  TYPE (CartesianCoordinates), DIMENSION(:), ALLOCATABLE :: &
       observer_arr, &
       topocenter_arr
  TYPE (CartesianCoordinates) :: &
       observer, &
       ccoord_first, &
       ccoord_last
  TYPE (Observatory) :: &
       obsy
  TYPE (Observatories) :: &
       obsies
  TYPE (Observation), DIMENSION(:), POINTER :: &
       obs_arr
  TYPE (Observation) :: &
       obs_first, &
       obs_last
  TYPE (Observations), DIMENSION(:), POINTER :: &
       obss_arr, &
       obss_sep_arr
  TYPE (Observations) :: &
       obss, &
       obss_
  TYPE (Orbit), DIMENSION(:), POINTER :: &
       orb_arr, &
       orb_arr_, &
       orb_arr_prm, &
       orb_arr_cmp, &
       orb_lt_corr_arr, &
       orb_arr_1, &
       orb_arr_2
  TYPE (Orbit) :: &
       orb
  TYPE (StochasticOrbit) :: &
       storb
  CHARACTER(len=FNAME_LEN), DIMENSION(:), POINTER :: &
       obsfnames, &
       outfnames
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:,:), allocatable :: &
       correct_linkages, &
       correct_linkages_tmp, &
       linkages, &
       linkages_tmp
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: &
       designation_arr
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), ALLOCATABLE :: &
       id_arr_prm, &
       id_arr_tmp, &
       id_arr
  CHARACTER(len=1024), DIMENSION(4) :: &
       header
  CHARACTER(len=32), DIMENSION(6) :: &
       element_str_arr, &
       stdev_str_arr
  CHARACTER(len=32), DIMENSION(5) :: &
       corr_str_arr
  CHARACTER(len=FNAME_LEN) :: &
       fname, &
       optfname, &
       orbidfname, &
       orbfname, &
       apriorifname, &
       addressfname, &
       triallinkfname, &
       obsfname
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: &
       ids
  CHARACTER(len=DESIGNATION_LEN) :: &
       id, id1, id2, des1, des2, id_prm
  CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), POINTER :: &
       element_type_inv_arr_prm, &
       element_type_inv_arr_cmp
  CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:,:), POINTER :: &
       obsy_code_arr
  CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: &
       obsy_code_arr_prm
  CHARACTER(len=ELEMENT_TYPE_LEN) :: &
       element_type_comp_prm, &
       element_type_in_prm, &
       element_type_out_prm, &
       element_type_inv, &
       comparison_variable_type
  CHARACTER(len=DYN_MODEL_LEN) :: &
       dyn_model
  CHARACTER(len=INTEGRATOR_LEN) :: &
       integrator
  CHARACTER(len=64) :: &
       orbid, &
       str, &
       method, &
       key2, &
       pidstr, &
       obs_type, &
       lstr
  REAL(bp), DIMENSION(:,:,:), POINTER :: &
       cov_arr, &
       cov_arr_prm, &
       cov_arr_cmp
  REAL(bp), DIMENSION(:,:), POINTER :: &
       pdf_arr2, &
       pdf_arr2_cmp, &
       residuals, &
       jac_arr_prm, &
       jac_arr_cmp, &
       jac_arr, &
       vov_map, &
       stdev_arr, &
       apriori_arr
  REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
       sigma_arr, &
       HG_arr_prm, &
       coord_arr, &
       element_arr, &
       elm_arr_1, &
       elm_arr_2, &
       bounds
  REAL(bp), DIMENSION(6,6) :: &
       cov, &
       cov_lt_corr, &
       corr
  REAL(bp), DIMENSION(6,2) :: &
       vov_scaling
  REAL(bp), DIMENSION(:), POINTER :: &
       pdf_arr, &
       pdf_arr_1, &
       pdf_arr_2, &
       pdf_arr_prm, &
       pdf_arr1_cmp, &
       chi2_ndof_arr_prm, &
       chi2_ndof_arr_cmp, &
       reg_apr_arr, &
       reg_apr_arr_prm, &
       reg_apr_arr_cmp, &
       eph_dt_since_last_obs, &
       pdf_lt_corr_arr
  REAL(bp), DIMENSION(:), ALLOCATABLE :: &
       dt_arr, &
       dist_arr, &
       coords, &
       width
  REAL(bp), DIMENSION(6) :: &
       elements, &
       stdevs, &
       coord, &
       lower_limit, &
       upper_limit, &
       mean_arr
  REAL(bp), DIMENSION(4) :: &
       sor_rho_init, &
       rho
  REAL(bp), DIMENSION(3) :: &
       pos, &
       d_min
  REAL(bp) :: &
       integration_step, &
       ls_correction_factor, &
       mjd, &
       dt, &
       day, &
       accwin_multiplier, &
       sor_genwin_multiplier, &
       stdev, &
       d, &
       cos_d, &
       c1, &
       c2, &
       c3, &
       c4, &
       c5, &
       c6, &
       chi2, &
       chi2_, &
       computed_upper_bound
  REAL(hp), DIMENSION(:,:), POINTER :: &
       indx_arr
  REAL(hp), DIMENSION(:), POINTER :: &
       trial_indxs
  INTEGER, DIMENSION(:), POINTER :: &
       nind_arr
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       indx_arr_, &
       box_nrs
  INTEGER, DIMENSION(3) :: &
       j_min, &
       k_min
  INTEGER :: &
       task, &
       err, &
       i, &
       j, &
       k, &
       l, &
       m, &
       n, &
       norb, &
       ntrial, &
       norb_sw, &
       ntrial_sw, &
       ninit, &
       nobsy, &
       nobs, &
       year, &
       month, &
       sor_niter, &
       nepoch, &
       nsimult, &
       var, &
       lu, &
       j_max, &
       k_max, &
       norb_old, &
       iter, &
       iindx, &
       nbin, &
       ifilter, &
       filter_start, &
       filter_stop, &
       ntriallink, &
       nlink, &
       loop_start, &
       loop_stop, &
       pid, &
       i2phase
  LOGICAL, DIMENSION(:,:), POINTER :: &
       obs_mask
  LOGICAL, DIMENSION(:), POINTER :: &
       eph_lt_correction_arr, &
       detected_correct_linkages
  LOGICAL, DIMENSION(:), ALLOCATABLE :: &
       elm
  LOGICAL, DIMENSION(10) :: &
       perturbers
  LOGICAL, DIMENSION(6) :: &
       ls_element_mask, &
       vov_mapping_mask
  LOGICAL :: &
       simulation, &
       obsfnames_from_file, &
       plot_results, &
       plot_open, &
       multiple_ids, &
       outlier_rejection_prm, &
       linkage, &
       common_epoch, &
       apriori, &
       found

  REAL(bp) :: dt1, dt2
  REAL(hp) :: indx
  INTEGER :: nset, nset_max, iset, iset_, norb_, iorb
  INTEGER :: ret

!!$  ! Set path to configuration file:
!!$  ! First, try the environment variable:
!!$  CALL getenv("OORB_CONF", conf_fname)
!!$  IF (LEN_TRIM(conf_fname) == 0) THEN
!!$     ! Second, if environment variable not defined use default:
!!$     conf_fname = "./oorb.conf"
!!$  END IF
!!$  ! Third, (if specified) the command-line option overrides the previous:
!!$  conf_fname = get_cl_option("--conf=",conf_fname)
!!$  CALL NEW(conf_file, conf_fname)
!!$  CALL setActionRead(conf_file)
!!$  CALL setStatusOld(conf_file)
!!$  CALL OPEN(conf_file)
!!$  IF (error) THEN
!!$     CALL errorMessage("oorb", &
!!$          "Configuration file missing or error while opening it.", 1)
!!$     STOP
!!$  END IF
!!$  CALL readConfigurationFile(conf_file, &
!!$       planetary_ephemeris_fname=planetary_ephemeris_fname, &
!!$       err_verb=err_verb, &
!!$       info_verb=info_verb, &
!!$       element_type_comp=element_type_comp_prm, &
!!$       element_type_out=element_type_out_prm, &
!!$       obs_stdev_arr=obs_stdev_arr_prm, &
!!$       observation_format_out=observation_format_out, &
!!$       orbit_format_out=orbit_format_out, &
!!$       outlier_rejection=outlier_rejection_prm, &
!!$       outlier_multiplier=outlier_multiplier_prm, &
!!$       dchi2_rejection=dchi2_rejection, &
!!$       dchi2_max=dchi2_max, &
!!$       plot_open=plot_open, &
!!$       plot_results=plot_results, &
!!$       dyn_model=dyn_model, &
!!$       perturbers=perturbers, &
!!$       integrator=integrator, &
!!$       integration_step=integration_step, &
!!$       relativity=relativity, &
!!$       simint=simint, &
!!$       pp_H_estimation=pp_H_estimation, &
!!$       pp_G=pp_G, &
!!$       pp_G_unc=pp_G_unc, &
!!$       asteroid_perturbers=asteroid_perturbers  )
!!$  IF (error) THEN
!!$     CALL errorMessage("oorb", &
!!$          "TRACE BACK (15)", 1)
!!$     STOP
!!$  END IF
!!$  IF (.NOT.ALL(obs_stdev_arr_prm < 0.0_bp)) THEN
!!$     WHERE (obs_stdev_arr_prm < 0.0_bp)
!!$        obs_stdev_arr_prm = 0.0_bp
!!$     END WHERE
!!$  END IF
!!$  CALL set_relativity(relativity)

  ! Set path to data files:
  CALL setAccessToDataFiles()
!!$  IF (LEN_TRIM(planetary_ephemeris_fname) == 0) THEN
!!$     planetary_ephemeris_fname = TRIM(EPH_FNAME)
!!$  END IF
!!$  CALL JPL_ephemeris_init(error, &
!!$       filename=TRIM(OORB_DATA_DIR) // "/" // TRIM(planetary_ephemeris_fname)) 
!!$  IF (error) THEN
!!$     CALL errorMessage("oorb", &
!!$          "Could not initialize planetary ephemerides using the " // &
!!$          TRIM(OORB_DATA_DIR) // "/" // TRIM(planetary_ephemeris_fname) // " file.", 1)
!!$     STOP
!!$  END IF

!!$  ! Set path to Gnuplot scripts using environment variable:
!!$  CALL getenv("OORB_GNUPLOT_SCRIPTS_DIR", gnuplot_scripts_dir)

  ! Get process ID
  pid = getpid()
  CALL toString(pid, pidstr, error)
  IF (error) THEN
     CALL errorMessage("ephemeris_linking", &
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
        CALL errorMessage("ephemeris_linking", &
             "TRACE BACK (5)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     CALL setPositionAppend(logfile)
     IF (error) THEN
        CALL errorMessage("ephemeris_linking", &
             "TRACE BACK (10)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     CALL OPEN(logfile)
     IF (error) THEN
        CALL errorMessage("ephemeris_linking", &
             "TRACE BACK (10)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     lu = getUnit(logfile)
     IF (error) THEN
        CALL errorMessage("ephemeris_linking", &
             "TRACE BACK (15)", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
  END IF
  obs_type = "mpc3"
  dyn_model = "2-body"
  perturbers = .TRUE.
  integrator = "bulirsch-stoer"
  integration_step = 10.0_bp
  info_verb = 1
  err_verb = 1
  CALL NEW(epoch, 2007, 1, 1.0_bp, "TT")
  IF (error) THEN
     CALL errorMessage("ephemeris_linking", &
          "TRACE BACK (20)", 1)
     CALL NULLIFY(logfile)
     STOP
  END IF
  ALLOCATE(epoch_arr(3))
  ALLOCATE(topocenter_arr(SIZE(epoch_arr)), &
       coords(3*SIZE(epoch_arr)), &
       elm(3*SIZE(epoch_arr)), &
       bounds(3*SIZE(epoch_arr),2), &
       width(3*SIZE(epoch_arr)), &
       box_nrs(3*SIZE(epoch_arr)))
  elm = .TRUE.
  CALL NEW(epoch_arr(1), 2004, 1, 24.0_bp, "TT")
  CALL NEW(epoch_arr(2), 2004, 1, 26.0_bp, "TT")
  CALL NEW(epoch_arr(3), 2004, 1, 28.0_bp, "TT")
  IF (error) THEN
     CALL errorMessage("ephemeris_linking", &
          "TRACE BACK (50)", 1)
     STOP
  END IF
  CALL NEW(obsies)
  DO i=1,SIZE(topocenter_arr)
     topocenter_arr(i) = getObservatoryCCoord(obsies, "500", epoch_arr(i))
     elm((i-1)*3+1) = .FALSE.
     bounds((i-1)*3+1:(i-1)*3+3,1) = (/   0.0_bp, 0.0_bp, -pi/2.0_bp /)
     bounds((i-1)*3+1:(i-1)*3+3,2) = (/ 100.0_bp, two_pi,  pi/2.0_bp /)
     width((i-1)*3+1:(i-1)*3+3)    = (/   0.1_bp, 3.0_bp*rad_amin, 3.0_bp*rad_amin /)
  END DO
  CALL NULLIFY(obsies)
  box_nrs = CEILING((bounds(:,2) - bounds(:,1))/width)
  DO i=1,SIZE(box_nrs)
     computed_upper_bound = bounds(i,1) + box_nrs(i)*width(i)
     IF (ABS(computed_upper_bound - bounds(i,2)) > &
          1000.0_bp*EPSILON(computed_upper_bound)) THEN
        box_nrs(i) = box_nrs(i) - 1
     END IF
  END DO
  indx = 1.0_hp
  DO i=1,SIZE(box_nrs)
     indx = indx*REAL(box_nrs(i),hp)
  END DO
  WRITE(lu,*) "Number of addresses in the binned phase space is ", indx
  c1 = 1.25_bp ! = 5/4
  c2 = 2.0_bp
  c3 = 2.0_bp
  c4 = 1.0E-01_bp
  c5 = 1.0E-01_bp
  ! Which filters are to be used?
  filter_start = get_cl_option("--filter_start=", 1)
  filter_stop = get_cl_option("--filter_stop=", 4)

  ! Write logfile

  WRITE(lu,"(1X,A,L1)") "Simulated observations... ", simulated_observations
  WRITE(lu,"(1X,A,6(F10.5,1X))") "Assumed observational RA and Dec uncertainties [asec]:", &
       stdevs(2:3)/rad_asec
  WRITE(lu,"(1X,2(A,1X,I0,1X))") "Using filters", filter_start, "through", filter_stop
  WRITE(lu,"(1X,A,1X,I0)") "Number of trial orbits:", norb


  ! For simulations, compute 2-linkages missing
  IF (simulated_observations) THEN
     obsfname = get_cl_option("--obs-in=", " ")
     IF (LEN_TRIM(obsfname) /= 0) THEN
        CALL NEW(obsfile, TRIM(obsfname))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1090)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(obsfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1095)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        WRITE(lu,"(A,A)") "Open observation file ", TRIM(obsfname)
        CALL OPEN(obsfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
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
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1105)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NULLIFY(obsfile)
     ELSE
        WRITE(stderr,*) "Could not read observation file ", TRIM(obsfname)
        CALL NULLIFY(logfile)
        STOP
     END IF
     WRITE(lu,"(A)") "Observation file read."
     obss_sep_arr => getSeparatedSets(obss)
     IF (error) THEN
        CALL errorMessage("ephemeris_linking", &
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
     k = 0
     DO i=1,nset-1
        des1 = getDesignation(obss_sep_arr(i))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1115)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        id1 = getID(obss_sep_arr(i))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1120)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        DO j=i+1,nset
           IF (des1 == getDesignation(obss_sep_arr(j))) THEN
              IF (error) THEN
                 CALL errorMessage("ephemeris_linking", &
                      "TRACE BACK (1125)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              id2 = getID(obss_sep_arr(j))
              IF (error) THEN
                 CALL errorMessage("ephemeris_linking", &
                      "TRACE BACK (1130)", 1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              k = k + 1
              IF (.NOT.AllOCATED(correct_linkages)) THEN
                 ALLOCATE(correct_linkages(k,2), stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "Could not allocate memory (1015).", 1)
                    CALL NULLIFY(logfile)
                    STOP
                 END IF
              ELSE IF (k > SIZE(correct_linkages,dim=1)) THEN
                 allocate(correct_linkages_tmp(size(correct_linkages,dim=1),2))
                 correct_linkages_tmp = correct_linkages
                 deallocate(correct_linkages)
                 allocate(correct_linkages(2*k,2))
                 correct_linkages = correct_linkages_tmp
                 deallocate(correct_linkages_tmp)
              END IF
              IF (TRIM(id1) < TRIM(id2)) THEN
                 correct_linkages(k,1:2) = (/ id1, id2 /)
              ELSE
                 correct_linkages(k,1:2) = (/ id2, id1 /)
              END IF
           END IF
        END DO
        CALL NULLIFY(obss_sep_arr(i))
     END DO
     CALL NULLIFY(obss_sep_arr(nset-1))
     CALL NULLIFY(obss_sep_arr(nset))        
     DEALLOCATE(obss_sep_arr, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("ephemeris_linking", &
             "Could not deallocate memory (1025).", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
     WRITE(lu,"(1X,A,I0,A)") "Altogether ", k, " correct 2-linkages to be detected."
     allocate(correct_linkages_tmp(k,2))
     correct_linkages_tmp = correct_linkages(1:k,:)
     deallocate(correct_linkages)
     allocate(correct_linkages(k,2))
     correct_linkages = correct_linkages_tmp
     deallocate(correct_linkages_tmp)
  END IF



  DO ifilter=filter_start,filter_stop

     SELECT CASE (ifilter)

     CASE (1)

        !!
        !! PHASE 1: FIND CANDIDATE LINKAGES (VIA ADDRESS COMPARISON OR FROM FILE)
        !!

        WRITE(lu,*)
        WRITE(lu,"(1X,A)") "Using the following parameters for the address comparison:"
        WRITE(lu,"(1X,A,6(6X,L1,4X))") " - parameters to be used : ", elm
        WRITE(lu,*) " - lower bounds        : ", bounds(1,1), bounds(2:3,1)/rad_deg 
        WRITE(lu,*) " - upper bounds        : ", bounds(1,2), bounds(2:3,2)/rad_deg 
        WRITE(lu,*) " - bin widths          : ", width(1), width(2:3)/rad_deg 
        WRITE(lu,*) " - bin numbers         : ", box_nrs
        DO i=1,SIZE(epoch_arr)
           WRITE(lu,"(1X,A,I0,2A)") " - common epoch #", i,"      : ", &
                getCalendarDateString(epoch_arr(i), "TT")
        END DO
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
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
        WRITE(lu,"(1X,A,1X,F32.0)") "Number of addresses in the binned phase space is", indx
        WRITE(lu,*)

        ! READ FILE CONTAINING ORBIT FILENAMES:
        fname = get_cl_option("--orb-fname-in=", " ")
        IF (LEN_TRIM(fname) == 0) THEN
           error = .TRUE.
           CALL errorMessage("ephemeris_linking", &
                "Orbit ID file missing.", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL NEW(orbidfile, TRIM(fname))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(orbidfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(orbidfile)  
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1020)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        nset_max = getNrOfLines(orbidfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1025)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        nset = 0
        DO
           READ(getUnit(orbidfile),"(A)",iostat=err) orbid
           IF (err > 0) THEN
              error = .TRUE.
              CALL errorMessage("ephemeris_linking", &
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
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (1035)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setStatusOld(orbfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (1040)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(orbfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (1045)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Orbit file opened."
           END IF
           norb = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (1050)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           ALLOCATE(id_arr_prm(norb), &
                orb_arr(norb), &
                stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("ephemeris_linking", &
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
                 CALL errorMessage("ephemeris_linking", &
                      "Could not read orbit file.", 1)
                 WRITE(lu,*) i
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              CALL setParameters(orb_arr(i), &
                   dyn_model=dyn_model, &
                   perturbers=perturbers, &
                   integrator=integrator, &
                   integration_step=integration_step)
              IF (error) THEN
                 CALL errorMessage("ephemeris_linking", &
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

           ! Compute addresses of the orbits and put them into the
           ! comparison table:
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // " Compute addresses."
           END IF
           nset = nset + 1
           norb = SIZE(orb_arr, dim=1)
           
           IF (.NOT.Allocated(id_arr_tmp)) THEN
              ALLOCATE(id_arr_tmp(nset_max), stat=err)
              IF (err /= 0) THEN
                 error = .TRUE.
                 CALL errorMessage("ephemeris_linking", &
                      'Could not allocate memory (1010).',1)
                 CALL NULLIFY(logfile)
                 STOP
              END IF
              id_arr_tmp = " "
              CALL init(address_tree)
           END IF

           CALL NEW(tmpfile,TRIM(orbid) // ".eph")
           CALL OPEN(tmpfile)
           orbloop:DO iorb=1,norb
              CALL getTopocentricSCoords(orb_arr(iorb), topocenter_arr, scoord_arr)
              IF (error) THEN
                 error = .FALSE.
                 CALL errorMessage('identification_mod / scoordComparison' , &
                      'TRACE BACK (200)', 1)
                 CYCLE
              END IF
              DO k=1,SIZE(scoord_arr,dim=1)
                 coords((k-1)*3+1:(k-1)*3+3) = getPosition(scoord_arr(k))
              END DO
              WRITE(getUnit(tmpfile),*) coords(2:3)/rad_deg, coords(5:6)/rad_deg, coords(8:9)/rad_deg
              DEALLOCATE(scoord_arr, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("ephemeris_linking", &
                      "Could not deallocate memory (15).", 1)
                 STOP
              END IF
              indx = arrayToReal(coords, elm, width, box_nrs, &
                   bounds, error)
              IF (error) THEN
                 CALL errorMessage("ephemeris_linking", &
                      'Could not transform array to integer.', 1)
                 STOP
              END IF
              CALL insert_tree_node(address_tree, indx, nset)
           END DO orbloop
           CALL NULLIFY(tmpfile)

           id_arr_tmp(nset) = id_arr_prm(1)
           IF (info_verb >= 1) THEN
              WRITE(lu,*) "[" // TRIM(orbid) // "] " // &
                   " Computation of addresses ready."
           END IF

           DEALLOCATE(id_arr_prm, &
                orb_arr, &
                stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("ephemeris_linking", &
                   "Could not deallocate memory (1005).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
        END DO
        allocate(id_arr(nset))
        id_arr = id_arr_tmp(1:nset)
        deallocate(id_arr_tmp)

        ret = system("ps ux | grep ephemeris_linking")

        ! Output addresses and corresponding sets, and build
        ! trial_linkage_tree:
        key2 = " "
        m = 0
        addressfname = get_cl_option("--address-out=", "address.out")
        CALL NEW(addressfile, TRIM(addressfname))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1080)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(addressfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1085)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL init(trial_linkage_tree)
        address_tree_node1 => minimum(address_tree, address_tree%root)
        DO WHILE (.NOT.ASSOCIATED(address_tree_node1, address_tree%nil))
           m = m + SIZE(address_tree_node1%data_nodes)
           DO k=1,SIZE(address_tree_node1%data_nodes)
              nset = SIZE(address_tree_node1%data_nodes(k)%data_array)
              WRITE(getUnit(addressfile),"(F32.0,1X,I10)",advance="NO") &
                   address_tree_node1%data_nodes(k)%key, nset
              DO i=1,nset
                 WRITE(getUnit(addressfile),"(1X,A16)",advance="NO") &
                      TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(i)))
              END DO
              WRITE(getUnit(addressfile),"(1X)")
              DO i=1,nset-1
                 DO j=i+1,nset
                    IF (TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(i))) < &
                         TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(j)))) THEN
                       key2 = TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(i))) // &
                            " " // TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(j)))
                    ELSE
                       key2 = TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(j))) // &
                            " " // TRIM(id_arr(address_tree_node1%data_nodes(k)%data_array(i)))
                    END IF
                    tree_node_ch32 => search(trial_linkage_tree, key2)
                    IF (ASSOCIATED(tree_node_ch32,trial_linkage_tree%nil)) THEN
                       CALL insert_tree_node(trial_linkage_tree, key2)
                    END IF
                 END DO
              END DO
              DEALLOCATE(address_tree_node1%data_nodes(k)%data_array, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("ephemeris_linking", &
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
              CALL errorMessage("ephemeris_linking", &
                   "Could not deallocate memory (1015).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           DEALLOCATE(address_tree_node2, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("ephemeris_linking", &
                   "Could not deallocate memory (1020).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
        END DO
        CALL delete_tree(address_tree)
        CALL NULLIFY(addressfile)
        WRITE(lu,*) m, " addresses used."

        ret = system("ps ux | grep ephemeris_linking")

        ! Output 2-linkages found
        triallinkfname = get_cl_option("--triallink-out=","triallink.out")
        CALL NEW(triallinkfile, TRIM(triallinkfname))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1080)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(triallinkfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (1085)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        i = 0
        tree_node_ch32 => minimum(trial_linkage_tree, trial_linkage_tree%root)
        DO WHILE (.NOT.ASSOCIATED(tree_node_ch32, trial_linkage_tree%nil))
           i = i + 1
           WRITE(getUnit(triallinkfile),"(A)") TRIM(tree_node_ch32%key)  
           tree_node_ch32 => successor(trial_linkage_tree, tree_node_ch32)
        END DO
        CALL NULLIFY(triallinkfile)
        WRITE(lu,*) "2-links found: ", i

        ! For simulations, compute 2-linkages missing
        IF (simulated_observations) THEN
           j = 0
           DO i=1,SIZE(correct_linkages,dim=1)
              key2 = TRIM(correct_linkages(i,1)) // " " // &
                   TRIM(correct_linkages(i,2))
              tree_node_ch32 => search(trial_linkage_tree, key2)
              IF (ASSOCIATED(tree_node_ch32,trial_linkage_tree%nil)) THEN
                 WRITE(lu,*) "Missing: ", TRIM(key2)
                 j = j + 1
              END IF
           END DO
           WRITE(lu,*) "Missing correct 2-links: ", j
        END IF

        ! Deallocate trial_linkage_tree
        tree_node_ch32_1 => minimum(trial_linkage_tree, trial_linkage_tree%root)
        DO WHILE (.NOT.ASSOCIATED(tree_node_ch32_1, trial_linkage_tree%nil))
           tree_node_ch32_2 => tree_node_ch32_1
           tree_node_ch32_1 => successor(trial_linkage_tree, tree_node_ch32_1)
           tree_node_ch32_2 => delete_tree_node(trial_linkage_tree, tree_node_ch32_2)
           DEALLOCATE(tree_node_ch32_2, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("ephemeris_linking", &
                   "Could not deallocate memory (1030).", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
        END DO
        CALL delete_tree(trial_linkage_tree)

     CASE (2)

        !!
        !! PHASE 2: FIND LINKING 2B ORBIT
        !!

        ! Open file with trial 2-linkages:
        triallinkfname = get_cl_option("--triallink-in=", triallinkfname)
        IF (LEN_TRIM(triallinkfname) == 0) THEN
           CALL errorMessage("ephemeris_linking", &
                "File with trial linkages not specified.", 1)
           CALL NULLIFY(logfile)
           STOP           
        END IF
        CALL NEW(triallinkfile, TRIM(triallinkfname))
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (2005)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL setStatusOld(triallinkfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (2010)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        CALL OPEN(triallinkfile)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (2015)", 1)
           CALL NULLIFY(logfile)
           STOP
        END IF
        ! Read name of observation file
        obsfname = get_cl_option("--obs-in=", " ")
        IF (LEN_TRIM(obsfname) == 0) THEN
           CALL errorMessage("ephemeris_linking", &
                "Observation file missing.", 1)
           CALL NULLIFY(logfile)
           STOP     
        END IF

        ntriallink = 0
        nlink = 0
        loop_start = get_cl_option("--loop-start=", 1)
        loop_stop = get_cl_option("--loop-stop=", HUGE(loop_stop))

        ! Get process ID
        pid = getpid()
        CALL toString(pid, pidstr, error)
        IF (error) THEN
           CALL errorMessage("ephemeris_linking", &
                "TRACE BACK (2020)", 1)
           STOP
        END IF

        ! Scan trial linkages one by one
        triallinkloop: DO
           ntriallink = ntriallink + 1
           IF (ntriallink > loop_stop) THEN
              EXIT triallinkloop
           END IF

           ! Read address
           READ(getUnit(triallinkfile),*,iostat=err) id1, id2
           IF (err > 0) THEN
              CALL NULLIFY(logfile)
              STOP "Error while reading address tree file."
           ELSE IF (err < 0) THEN
              EXIT triallinkloop
           END IF
           IF (ntriallink < loop_start) THEN
              CYCLE triallinkloop
           END IF

           ! Trial linkage read
           linkage = .FALSE.
           method = " "

           ret = system("rm -f ephemeris_linking_tmp" // &
                TRIM(pidstr) // "." // TRIM(obs_type))
           ret = system("grep " // TRIM(id1) // &
                " " // TRIM(obsfname) // &
                " > ephemeris_linking_tmp" // TRIM(pidstr) // &
                "." // TRIM(obs_type) // &
                " ; grep " // TRIM(id2) // &
                " " // TRIM(obsfname) // &
                " >> ephemeris_linking_tmp" // TRIM(pidstr) // &
                "." // TRIM(obs_type))
           CALL NEW(obsfile, "ephemeris_linking_tmp" // &
                TRIM(pidstr) // "." // TRIM(obs_type))
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2025)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setStatusOld(obsfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2030)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL OPEN(obsfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2035)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NEW(obss, obsfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2040)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL NULLIFY(obsfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2045)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           CALL setNumber(obss,1)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2050)", 1)
              CALL NULLIFY(logfile)
              STOP
           END IF
           obs_first = getObservation(obss, 1)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2055)", 1)
              STOP
           END IF
           ccoord_first = getObservatoryCCoord(obs_first)
           nobs = getNrOfObservations(obss)
           obs_last = getObservation(obss, nobs)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2060)", 1)
              STOP
           END IF
           ccoord_last = getObservatoryCCoord(obs_last)
           coord = getCoordinates(ccoord_first) - getCoordinates(ccoord_last)
           rho(1:2) = (/ 0.0_bp, 100.0_bp /)
           dt = getObservationalTimespan(obss)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2065)", 1)
              STOP
           END IF
           ! Maximum bounds according to different observatory coordinates +
           ! additional width due to motion of target (0.05 AU/d ~ 87 km/s)
           rho(4) = SQRT(DOT_PRODUCT(coord(1:3),coord(1:3))) + dt*0.05_bp
           rho(3) = -rho(4)
           t = getTime(obs_first)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2070)", 1)
              STOP
           END IF
           mjd = getMJD(t, "TT")
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2075)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t_inv, mjd, "TT")
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (2080)", 1)
              STOP
           END IF
           CALL NULLIFY(obs_first)
           CALL NULLIFY(obs_last)

           DO i2phase=1,1!3

              SELECT CASE (i2phase)

              CASE (1)

                 ! Stepwise Ranging
                 CALL NULLIFY(storb)
                 CALL NEW(storb, obss)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (2085)", 1)
                    STOP
                 END IF
                 CALL setParameters(storb, &
                      dyn_model=dyn_model, &
                      integrator=integrator, &
                      integration_step=integration_step, &
                      t_inv=t_inv, &
                      element_type="cartesian", &
                      regularized_pdf=.FALSE., &
                      jacobians_pdf=.FALSE., &
                      accept_multiplier=4.0_bp, &
                      generat_multiplier=4.0_bp, &
                      sor_rho1_l=rho(1), &
                      sor_rho1_u=rho(2), &
                      sor_rho2_l=rho(3), &
                      sor_rho2_u=rho(4), &
                      sor_norb=50, &
                      sor_ntrial=10000, &
                      sor_niter=1, &
                      sor_norb_sw=100, &
                      sor_ntrial_sw=10000, &
                      sor_2point_method="continued fraction")
                 IF (error) THEN
                    STOP
                 END IF
                 info_verb = 0
                 CALL stepwiseRanging(storb, nobs_max=-1)
                 info_verb = 1
                 IF (error) THEN
                    error = .FALSE.
                 ELSE
                    linkage = .TRUE.
                    method = "RANGING"
                    orb_arr => getSampleOrbits(storb)
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (2090)", 1)
                       STOP
                    END IF
                    chi2_ = HUGE(chi2_)
                    DO j=1,SIZE(orb_arr)
                       chi2 = getChi2(storb, orb_arr(j))
                       IF (error) THEN
                          CALL errorMessage("ephemeris_linking", &
                               "TRACE BACK (2095)", 1)
                          STOP
                       END IF
                       IF (chi2 < chi2_) THEN
                          chi2_ = chi2
                          k = j
                       END IF
                    END DO
                    orb = copy(orb_arr(k))
                    DEALLOCATE(orb_arr, stat=err)
                    IF (err /= 0) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "Could not deallocate memory (57).", 1)
                       STOP
                    END IF
                 END IF

              CASE (2)

                 ! Similar orbits
                 ret = system("gunzip -c " // TRIM(id1) &
                      // ".orb.gz > " // TRIM(id1) // &
                      ".orb")
                 CALL NEW(orbfile, TRIM(id1) // ".orb")
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (390)", 1)
                    STOP
                 END IF
                 CALL setStatusOld(orbfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (395)", 1)
                    STOP
                 END IF
                 CALL OPEN(orbfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (400)", 1)
                    STOP
                 END IF
                 WRITE(lu,*) "[" // TRIM(id1) // "] " // " Orbit file opened."
                 norb = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (405)", 1)
                    STOP
                 END IF
                 ALLOCATE(orb_arr_1(norb), &
                      stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "Could not allocate memory (20).", 1)
                    STOP
                 END IF
                 header(1:4)(:) = " "
                 WRITE(lu,*) "[" // TRIM(id1) // "] " // " Reading orbit file."
                 DO j=1,norb
                    CALL readOpenOrbOrbitFile(getUnit(orbfile), &
                         header, &
                         element_type_in=element_type_in_prm, &
                         id=id_prm, &
                         orb=orb_arr_1(j))
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "Could not read orbit file.", 1)
                       WRITE(stderr,*) j
                       STOP
                    END IF
                    CALL setParameters(orb_arr_1(j), dyn_model=dyn_model, &
                         integrator=integrator, integration_step=integration_step)
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (410)", 1)
                       STOP
                    END IF
                 END DO
                 CALL NULLIFY(orbfile)
                 ret = system("rm -f " // TRIM(id1) // ".orb")
                 WRITE(lu,*) "[" // TRIM(id1) // "] " // " Orbit file read."

                 ret = system("gunzip -c " // TRIM(id2) &
                      // ".orb.gz > " // TRIM(id2) // &
                      ".orb")
                 CALL NEW(orbfile,TRIM(id2) // ".orb")
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (415)", 1)
                    STOP
                 END IF
                 CALL setStatusOld(orbfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (420)", 1)
                    STOP
                 END IF
                 CALL OPEN(orbfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (425)", 1)
                    STOP
                 END IF
                 WRITE(lu,*) "[" // TRIM(id2) // "] " // " Orbit file opened."
                 norb = getNrOfLines(orbfile) - 4 ! 4 header lines are not taken into count
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (430)", 1)
                    STOP
                 END IF
                 ALLOCATE(orb_arr_2(norb), stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "Could not allocate memory (20).", 1)
                    STOP
                 END IF
                 header(1:4)(:) = " "
                 WRITE(lu,*) "[" // TRIM(id2) // "] " // " Reading orbit file."
                 DO j=1,norb
                    CALL readOpenOrbOrbitFile(getUnit(orbfile), &
                         header, &
                         element_type_in=element_type_in_prm, &
                         id=id_prm, &
                         orb=orb_arr_2(j))
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "Could not read orbit file.", 1)
                       WRITE(stderr,*) j
                       STOP
                    END IF
                    CALL setParameters(orb_arr_2(j), dyn_model=dyn_model, &
                         integrator=integrator, integration_step=integration_step)
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (435)", 1)
                       STOP
                    END IF
                 END DO
                 CALL NULLIFY(orbfile)
                 ret = system("rm -f " // TRIM(id2) // ".orb")
                 WRITE(lu,*) "[" // TRIM(id2) // "] " // " Orbit file read."

                 ALLOCATE(elm_arr_1(SIZE(orb_arr_1,dim=1),7), &
                      elm_arr_2(SIZE(orb_arr_2,dim=1),7))
                 err_verb = 0
                 k = 1
                 DO j=1,SIZE(orb_arr_1,dim=1)
                    elm_arr_1(k,1:6) = getElements(orb_arr_1(j), "keplerian", "ecliptic")
                    IF (error) THEN
                       error = .FALSE.
                    ELSE
                       elm_arr_1(k,7) = MOD(SUM(elm_arr_1(k,4:5)),two_pi)
                       k = k + 1
                    END IF
                 END DO
                 j_max = k - 1
                 k = 1
                 DO j=1,SIZE(orb_arr_2,dim=1)
                    elm_arr_2(k,1:6) = getElements(orb_arr_2(j), "keplerian", "ecliptic")
                    IF (error) THEN
                       error = .FALSE.
                    ELSE
                       elm_arr_2(k,7) = MOD(SUM(elm_arr_2(k,4:5)),two_pi)
                       k = k + 1
                    END IF
                 END DO
                 k_max = k - 1
                 WRITE(lu,*) "[" // TRIM(id2) // "] " // " Orbits transformed."
                 err_verb = 1
                 d_min = HUGE(d_min)
                 DO j=1,j_max
                    DO k=1,k_max
                       ! Different choices for the D parameter
                       ! From Nesvorny et al. 2002, Nature 417
                       !d = SQRT(mu/elm_arr_1(j,1)) * SQRT( &
                       !     c1*((elm_arr_1(j,1)-elm_arr_2(k,1))/elm_arr_1(j,1))**2 + &
                       !     c2*(elm_arr_1(j,2)-elm_arr_2(k,2))**2 + &
                       !     c3*(SIN(elm_arr_1(j,3))-SIN(elm_arr_2(k,3)))**2)
                       ! Modified Nesvorny et al.
                       d = SQRT(planetary_mu(3)/elm_arr_1(j,1)) * SQRT( &
                            c1*((elm_arr_1(j,1)-elm_arr_2(k,1))/elm_arr_1(j,1))**2 + &
                            c2*(elm_arr_1(j,2)-elm_arr_2(k,2))**2 + &
                            c3*(SIN(elm_arr_1(j,3))-SIN(elm_arr_2(k,3)))**2 + &
                            c4*(SIN(elm_arr_1(j,7))-SIN(elm_arr_2(k,7)))**2)! + &
                       !    c5*(sin(elm_arr_1(j,5))-sin(elm_arr_2(k,5)))**2)
                       IF (d < d_min(1)) THEN
                          d_min(1) = d
                          j_min(1) = j
                          k_min(1) = k
                       ELSE IF (d < d_min(2)) THEN
                          d_min(2) = d
                          j_min(2) = j
                          k_min(2) = k
                       ELSE IF (d < d_min(3)) THEN
                          d_min(3) = d
                          j_min(3) = j
                          k_min(3) = k
                       END IF
                    END DO
                 END DO
                 IF (info_verb >= 2) THEN
                    DO j=1,3
                       WRITE(lu,*) d_min(j)
                       WRITE(lu,*) elm_arr_1(j_min(j),1:2), elm_arr_1(j_min(j),3:6)/rad_deg
                       WRITE(lu,*) elm_arr_2(k_min(j),1:2), elm_arr_2(k_min(j),3:6)/rad_deg
                    END DO
                 END IF
                 elements = 0.0_bp
                 ALLOCATE(orb_arr(7), orb_arr_(7))
                 DO j=1,3
                    t = getTime(orb_arr_1(1))
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (440)", 1)
                       STOP
                    END IF
                    CALL NEW(orb_arr((j-1)*2+1), elm_arr_1(j_min(j),:), "keplerian", "ecliptic", t)
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (445)", 1)
                       STOP
                    END IF
                    CALL NULLIFY(t)
                    t = getTime(orb_arr_2(1))
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (450)", 1)
                       STOP
                    END IF
                    CALL NEW(orb_arr((j-1)*2+2), elm_arr_2(k_min(j),:), "keplerian", "ecliptic", t)
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (455)", 1)
                       STOP
                    END IF
                    CALL NULLIFY(t)
                 END DO
                 t = getTime(orb_arr_2(1))
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (460)", 1)
                    STOP
                 END IF
                 elements = elm_arr_1(j_min(1),:) + elm_arr_2(k_min(1),:)
                 elements = elements/2.0_bp
                 CALL NEW(orb_arr(7), elements, "keplerian", "ecliptic", t)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (465)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(t)
                 DEALLOCATE(orb_arr_1, orb_arr_2, elm_arr_1, elm_arr_2, stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "Could not allocate memory (5).", 1)        
                    STOP
                 END IF
                 DO j=1,SIZE(orb_arr,dim=1)
                    CALL setParameters(orb_arr(j), dyn_model=dyn_model, &
                         integrator=integrator, integration_step=integration_step)
                    CALL propagate(orb_arr(j), t_inv)
                    IF (error) THEN
                       CALL errorMessage("ephemeris_linking", &
                            "TRACE BACK (470)", 1)
                       STOP
                    END IF
                    orb_arr_(j) = copy(orb_arr(j))
                 END DO
                 stdev_arr => getStandardDeviations(obss)

                 ! Amoeba
                 CALL NULLIFY(storb)
                 CALL NEW(storb, obss)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (475)", 1)
                    STOP
                 END IF
                 CALL setParameters(storb, &
                      dyn_model=dyn_model, &
                      integrator=integrator, &
                      integration_step=integration_step, &
                      t_inv=t_inv, &
                      element_type="cartesian", &
                      smplx_niter=1000, &
                      smplx_tol=2*nobs*2.0_bp)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (480)", 1)
                    STOP
                 END IF
                 !                 iter = 1000
                 !                 CALL amoebaOrbits(storb, orb_arr, ftol=2*nobs*2.0_bp, iter=iter)
                 CALL simplexOrbits(storb, orb_arr)
                 IF ((error .AND. iter >= 1000) .OR. .NOT.error) THEN
                    error = .FALSE.
                    IF (info_verb >= 2) THEN
                       WRITE(lu,*) "Number of function evaluations in amoebaOrbits: ", iter
                    END IF
                    j = 0
                    DO WHILE (.NOT.linkage)
                       j = j + 1
                       IF (j > 7) THEN
                          EXIT
                       END IF
                       residuals => getResiduals(storb, orb_arr(j))
                       IF (error) THEN
                          error = .FALSE.
                          DEALLOCATE(residuals, stat=err)
                          CALL errorMessage("ephemeris_linking", &
                               "TRACE BACK (485)", 1)
                          CYCLE
                       END IF
                       IF (ALL(ABS(residuals(:,2:3)) < 4.0_bp*stdev_arr(:,2:3))) THEN
                          linkage = .TRUE.
                          method = "AMOEBA"
                          orb = copy(orb_arr(j))
                       END IF
                       DEALLOCATE(residuals, stat=err)
                       IF (err /= 0) THEN
                          CALL errorMessage("ephemeris_linking", &
                               "Could not deallocate memory (45).", 1)
                          STOP
                       END IF
                    END DO
                 ELSE IF (error .AND. iter < 1000) THEN
                    error = .FALSE.
                    DO j=1,7
                       CALL NULLIFY(orb_arr(j))
                       orb_arr(j) = copy(orb_arr_(j))
                    END DO
                 END IF

              CASE (3)

                 ! Least squares
                 CALL NULLIFY(storb)
                 CALL NEW(storb, obss)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (490)", 1)
                    STOP
                 END IF
                 ls_element_mask = .TRUE.
                 CALL setParameters(storb, &
                      dyn_model=dyn_model, &
                      integrator=integrator, &
                      integration_step=integration_step, &
                      t_inv=t_inv, &
                      element_type="cartesian", &
                      accept_multiplier=4.0_bp, &
                      ls_correction_factor=1.0_bp, &
                      ls_element_mask=ls_element_mask)
                 DO j=1,7
                    CALL leastSquares(storb, orb_arr(j))
                    IF (error) THEN
                       error = .FALSE.
                    ELSE
                       orb = getNominalOrbit(storb)
                       residuals => getResiduals(storb, orb)
                       IF (ALL(ABS(residuals(:,2:3)) < 4.0_bp*stdev_arr(:,2:3))) THEN
                          linkage = .TRUE.
                          method = "LS"
                          IF (info_verb >= 2) THEN
                             CALL writeNominalSolution(storb, obss, "cartesian", stdout)
                          END IF
                       END IF
                       DEALLOCATE(residuals, stat=err)
                       IF (err /= 0) THEN
                          CALL errorMessage("ephemeris_linking", &
                               "Could not deallocate memory (50).", 1)
                          STOP
                       END IF
                       IF (linkage) THEN
                          EXIT
                       END IF
                       CALL NULLIFY(orb)
                    END IF
                 END DO
                 DEALLOCATE(orb_arr, orb_arr_, stdev_arr, stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "Could not deallocate memory (55).", 1)
                    STOP
                 END IF

              END SELECT

              IF (linkage) THEN
                 nlink = nlink + 1
                 IF (.NOT.allocated(linkages)) THEN
                    ALLOCATE(linkages(nlink,2), stat=err)
                    IF (err /= 0) THEN
                       error = .TRUE.
                       CALL errorMessage("ephemeris_linking", &
                            'Could not allocate memory (10).',1)
                       STOP
                    END IF
                 ELSE IF (nlink > SIZE(linkages,dim=1)) THEN
                    allocate(linkages_tmp(size(linkages,dim=1),2))
                    linkages_tmp = linkages
                    deallocate(linkages)
                    allocate(linkages(2*nlink,2))
                    linkages = linkages_tmp
                    deallocate(linkages_tmp)
                 END IF

                 linkages(nlink,:) = (/ id1, id2 /)
                 CALL NEW(tmpfile,"lnk.out")
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (500)", 1)
                    STOP
                 END IF
                 CALL setPositionAppend(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (505)", 1)
                    STOP
                 END IF
                 CALL OPEN(tmpfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (510)", 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmpfile),"(3(A,1X))") TRIM(id1), "<-->", TRIM(id2)
                 CALL  NULLIFY(tmpfile)

                 ! Write linked observation sets:
                 CALL NEW(obsfile,"lnk.mpc2")
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (515)", 1)
                    STOP
                 END IF
                 CALL setPositionAppend(obsfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (520)", 1)
                    STOP
                 END IF
                 CALL OPEN(obsfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (525)", 1)
                    STOP
                 END IF
                 CALL toString(l, lstr, error)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (527)", 1)
                    STOP
                 END IF
                 CALL writeObservationFile(obss, getUnit(obsfile), "mpc2", number=lstr)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (530)", 1)
                    STOP
                 END IF
                 CALL  NULLIFY(obsfile)

                 ! Write orbit that link the observation sets:
                 CALL NEW(orbfile,"lnk.orb")
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (535)", 1)
                    STOP
                 END IF
                 CALL setPositionAppend(orbfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (540)", 1)
                    STOP
                 END IF
                 CALL OPEN(orbfile)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (545)", 1)
                    STOP
                 END IF
                 CALL toString(nlink, str, error)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (550)", 1)           
                    STOP
                 END IF
                 DO WHILE (LEN_TRIM(str) < 7)
                    str = "0" // TRIM(str)
                 END DO
                 CALL writeopenOrbOrbitFile(getUnit(orbfile), &
                      print_header=nlink==1, &
                      element_type_out="cartesian", &
                      id=TRIM(str), &
                      orb=orb)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (555)", 1)
                    STOP
                 END IF
                 CALL  NULLIFY(orbfile)
                 CALL setNumber(obss,0)
                 designation_arr => getObjects(obss)
                 IF (error) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "TRACE BACK (560)", 1)
                    STOP
                 END IF
                 str = " "
                 IF (SIZE(designation_arr) == 1) THEN
                    str = "CORRECT"
                 ELSE
                    str = "ERRONEOUS"
                 END IF
                 DEALLOCATE(designation_arr, stat=err)
                 IF (err /= 0) THEN
                    CALL errorMessage("ephemeris_linking", &
                         "Could not deallocate memory (60).", 1)
                    STOP
                 END IF
                 WRITE(lu,*) "LINKAGE: ", TRIM(id1), " <--> ", &
                      TRIM(id2), " ", TRIM(method), " ", TRIM(str)
                 EXIT

              END IF

           END DO

           IF (.NOT.linkage) THEN
              WRITE(lu,*) "NO LINKAGE: ", TRIM(id1), " <--> ", TRIM(id2)
           END IF
           CALL NULLIFY(obss)
           CALL NULLIFY(t_inv)
           CALL NULLIFY(storb)
           CALL NULLIFY(orb)

        END DO triallinkloop

        ret = system("rm -f ephemeris_linking_tmp" // &
             TRIM(pidstr) // "." // TRIM(obs_type))
        IF (nlink /= 0) THEN
           allocate(linkages_tmp(nlink,2))
           linkages_tmp = linkages(1:nlink,:)
           deallocate(linkages)
           allocate(linkages(nlink,2))
           linkages = linkages_tmp
           deallocate(linkages_tmp)
           WRITE(lu,*) "Found ", nlink, " linkages."
        ELSE
           WRITE(lu,*) "Found zero linkages."
        END IF

        IF (simulated_observations) THEN

           ! Find all detected correct linkages (indicated by a common designation)
           ALLOCATE(detected_correct_linkages(SIZE(correct_linkages,dim=1)))
           detected_correct_linkages = .FALSE.
           DO i=1,SIZE(linkages,dim=1)
              DO j=1,SIZE(correct_linkages,dim=1)
                 IF ((linkages(i,1) == correct_linkages(j,1) &
                      .AND. linkages(i,2) == correct_linkages(j,2)) &
                      .OR. &
                      (linkages(i,2) == correct_linkages(j,1) &
                      .AND. linkages(i,1) == correct_linkages(j,2))) THEN
                    detected_correct_linkages(j) = .TRUE.
                    EXIT
                 END IF
              END DO
           END DO

           ! Write all correct linkages that were detected to a file:
           CALL NEW(tmpfile,"det_cor_lnk.out")
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (565)", 1)
              STOP
           END IF
           CALL OPEN(tmpfile)
           IF (error) THEN
              CALL errorMessage("ephemeris_linking", &
                   "TRACE BACK (570)", 1)
              STOP
           END IF
           DO i=1,SIZE(correct_linkages,dim=1)
              IF (detected_correct_linkages(i)) THEN
                 WRITE(getUnit(tmpfile),"(3(1X,A))") TRIM(correct_linkages(i,1)), &
                      "<-->", TRIM(correct_linkages(i,2))
              END IF
           END DO
           CALL  NULLIFY(tmpfile)

           WRITE(lu,"(1X,I0,1X,A,1X,I0,1X,A,1X,F10.8,A)") &
                SIZE(linkages,dim=1), &
                "linkages found out of", &
                SIZE(id_arr,dim=1)*(SIZE(id_arr,dim=1)-1), &
                "possibilities (acceptance ratio:", &
                SIZE(linkages,dim=1) / &
                REAL(SIZE(id_arr,dim=1)*(SIZE(id_arr,dim=1)-1)), &
                ")."
           WRITE(lu,"(1X,I0,1X,A)") &
                COUNT(detected_correct_linkages), &
                "correct linkages were detected."
           WRITE(lu,"(1X,I0,1X,A)") &
                SIZE(detected_correct_linkages) - &
                COUNT(detected_correct_linkages), &
                "correct linkages were not detected."

           ! Write all correct linkages that weren't detected to a file:
           IF (SIZE(detected_correct_linkages) - &
                COUNT(detected_correct_linkages) > 0) THEN
              CALL NEW(tmpfile,"not_det_cor_lnk.out")
              IF (error) THEN
                 CALL errorMessage("ephemeris_linking", &
                      "TRACE BACK (575)", 1)
                 STOP
              END IF
              CALL OPEN(tmpfile)
              IF (error) THEN
                 CALL errorMessage("ephemeris_linking", &
                      "TRACE BACK (580)", 1)
                 STOP
              END IF
              DO i=1,SIZE(detected_correct_linkages,dim=1)
                 IF (.NOT.detected_correct_linkages(i)) THEN
                    WRITE(getUnit(tmpfile),"(3(1X,A))") TRIM(correct_linkages(i,1)), &
                         "<-->", TRIM(correct_linkages(i,2))
                 END IF
              END DO
              CALL NULLIFY(tmpfile)
           END IF

           DEALLOCATE(detected_correct_linkages, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("ephemeris_linking", &
                   "Could not deallocate memory (65).", 1)
              STOP
           END IF

        END IF

        CALL NULLIFY(triallinkfile)

     END SELECT

  END DO

  IF (AllOCATED(id_arr)) THEN
     DEALLOCATE(id_arr, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("ephemeris_linking", &
             "Could not deallocate memory (70).", 1)
        STOP
     END IF
  END IF

  IF (AllocATED(linkages)) THEN
     DEALLOCATE(linkages, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("ephemeris_linking", &
             "Could not deallocate memory (70).", 1)
        STOP
     END IF
  END IF

  IF (AllOCATED(correct_linkages)) THEN
     DEALLOCATE(correct_linkages, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("ephemeris_linking", &
             "Could not deallocate memory (1035).", 1)
        CALL NULLIFY(logfile)
        STOP
     END IF
  END IF

END PROGRAM ephemeris_linking



