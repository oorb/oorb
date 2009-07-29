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
!! *Description*:
!!
!! Main program for various tasks that include orbit computation.
!!
!! @author  MG
!! @version 2009-07-28
!!
PROGRAM oorb

  USE Base_cl
  USE File_cl
  USE Time_cl
  USE CartesianCoordinates_cl
  USE SphericalCoordinates_cl
  USE Observatories_cl
  USE Orbit_cl
  USE Observation_cl
  USE Observations_cl
  USE StochasticOrbit_cl
  USE PhysicalParameters_cl
  USE cl_options
  USE planetary_data  
  USE io

  IMPLICIT NONE
  TYPE (PhysicalParameters) :: &
       physparam
  TYPE (StochasticOrbit), DIMENSION(:), ALLOCATABLE :: &
       storb_arr_in
  TYPE (StochasticOrbit) :: &
       storb
  TYPE (Orbit), DIMENSION(:), POINTER :: &
       orb_arr_in
  TYPE (Orbit), DIMENSION(:), POINTER :: &
       orb_arr, &
       orb_arr_cmp, &
       orb_lt_corr_arr
  TYPE (Orbit) :: &
       orb, &
       ref_orb
  TYPE (Observations), DIMENSION(:), POINTER :: &
       obss_sep
  TYPE (Observations) :: &
       obss_in  
  TYPE (Observation) :: &
       obs
  TYPE (Observatories) :: &
       obsies
  TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: &
       ephemerides_arr
  TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
       ephemerides
  TYPE (SphericalCoordinates) :: &
       scoord
  TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
       observers
  TYPE (CartesianCoordinates) :: &
       ccoord, &
       obsy_ccoord
  TYPE (Time) :: &
       epoch, &
       epoch0, &
       t
  TYPE (File) :: &
       conf_file, &                                                 !! Configuration file.
       obs_file, &                                                  !! Generic observation file.
       orb_in_file, &                                               !! Input orbit file.
       orb_out_file, &                                              !! Output orbit file.
       out_file, &                                                  !! Generic output file.
       tmp_file                                                     !! Generic temporary file.
  CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: &
       element_type_pdf_arr_in                                      !! Element type of input orbital-element PDF.
  CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: &
       obsy_code_arr                                                !! IAU/MPC designated observatory code.
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: &
       group_name_arr, &
       id_arr_in
  CHARACTER(len=32), DIMENSION(:), ALLOCATABLE :: &
       str_arr
  CHARACTER(len=32), DIMENSION(6) :: &
       element_str_arr, &
       stdev_str_arr
  CHARACTER(len=32), DIMENSION(5) :: &
       corr_str_arr
  CHARACTER(len=1024), DIMENSION(4) :: &
       header                                                       !! Generic header.
  CHARACTER(len=FNAME_LEN) :: &
       conf_fname, &                                                !! Path to configuration file (incl. fname). 
       gnuplot_scripts_dir, &                                       !! Path to Gnuplot scripts directory. 
       obs_fname, &                                                 !! Path to observation file (incl. fname).
       orb_in_fname, &                                              !! Path to input orbit file (incl. fname).
       orb_out_fname, &                                             !! Path to output orbit file (incl. fname).
       out_fname, &                                                 !! Path to generic output file (incl. fname).
       tmp_fname
  CHARACTER(len=ELEMENT_TYPE_LEN) :: &
       element_type_comp_prm, &                                     !! Element type to be used in computations.
       element_type_in, &                                           !! Element type of input orbit(s).
       element_type_out_prm                                         !! Element type of output orbit(s).
  CHARACTER(len=DYN_MODEL_LEN) :: &
       dyn_model, &                                                 !! Dynamical model.
       dyn_model_init                                               !! Dynamical model used to propagate init orbit to epoch in LSL.
  CHARACTER(len=INTEGRATOR_LEN) :: &
       integrator, &                                                !! Integrator.
       integrator_init                                              !! Integrator used to propagate init orbit to epoch in LSL.
  CHARACTER(len=256) :: &
       frmt, &
       str, &
       suffix, &
       task
  CHARACTER(len=64) :: &
       sor_2point_method, &                                         !! Method used for 2-point boundary-value problem. 
       sor_2point_method_sw                                         !! Method used for 2-point boundary-value problem in stepwise Ranging.
  CHARACTER(len=DESIGNATION_LEN) :: &
       id                                                           !! Number or temporary designation of asteroid.
  CHARACTER(len=8) :: &
       observation_format_out, &                                    !! Format of output observations
       orbit_format_out                                             !! Format of output orbits
  CHARACTER(len=OBSY_CODE_LEN) :: &
       obsy_code                                                    !! IAU/MPC designated observatory code.
  REAL(bp), DIMENSION(:,:,:), POINTER :: &
       cov_arr, &
       cov_arr_in, &                                                !! Input array of covariance matrices.
       encounters
  REAL(bp), DIMENSION(:,:), POINTER :: &
       apoapsis_distance_pdf, &
       HG_arr_in, &
       jac_arr_cmp, &
       pdfs_arr, &
       periapsis_distance_pdf, &
       planeph
  REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
       elements_arr, &
       jac_arr_in, &
       temp_arr
  REAL(bp), DIMENSION(6,6) :: &
       corr, &
       cov
  REAL(bp), DIMENSION(2,2) :: &
       sor_rho_cmp
  REAL(bp), DIMENSION(:), POINTER :: &
       rchi2_arr_cmp, &
       pdf_arr_cmp, &
       pdf_arr_in, &
       reg_apr_arr_cmp, &
       weight_arr
  REAL(bp), DIMENSION(:), ALLOCATABLE :: &
       arc_arr, &
       rchi2_arr_in, &
       reg_apr_arr_in, &
       real_arr
  REAL(bp), DIMENSION(6) :: &
       comp_coord, &
       coordinates, &
       elements, &
       h_ecl_car_coord_obj, &
       h_ecl_car_coord_obsy, &
       obs_stdev_arr_prm, &
       stdev_arr
  REAL(bp), DIMENSION(4) :: &
       sor_genwin_offset, &
       sor_rho_init
  REAL(bp), DIMENSION(3) :: &
       geoc_obsy, &
       obsy_moon, &
       obsy_obj, &
       obsy_pos, &
       obsy_sun, &
       pos, &
       pos_opp, &
       sun_moon, &
       vec3
  REAL(bp) :: &
       Delta, &
       H_value, &
       G_value, &
       accwin_multiplier, &
       apriori_a_max, &
       apriori_a_min, &
       apriori_apoapsis_max, &
       apriori_apoapsis_min, &
       apriori_periapsis_max, &
       apriori_periapsis_min, &
       apriori_rho_min, &
       cos_obj_phase, &
       day0, day1, &
       dDelta, &
       ddec, &
       dec, &
       dra, &
       dt, dt_fulfill_night, &
       ephemeris_r2, &
       H_max, &
       hdist, &
       heliocentric_r2, &
       hlat, &
       hlon, &
       hoclat, &
       hoclon, &
       i_min, i_max, &
       integration_step, &
       integration_step_init, &
       ls_correction_factor, &
       lunar_alt, lunar_alt_max, &
       lunar_elongation, lunar_elongation_min, &
       lunar_phase, lunar_phase_min, lunar_phase_max, &
       mjd, mjd_tai, mjd_tt, mjd_utc, &
       moid, &
       obj_alt, obj_alt_min, &
       obj_phase, &
       obj_vmag, obj_vmag_max, &
       observer_r2, &
       obsy_moon_r2, &
       opplat, &
       opplon, &
       outlier_multiplier_prm, &
       pdf_ml_init, &
       pp_G, pp_G_unc, &
       ra, &
       sec, &
       solar_elongation, solar_elon_min, solar_elon_max, &
       solar_alt, solar_alt_max, &
       sor_genwin_multiplier, &
       stdev, &
       step, &
       sun_moon_r2, &
       timespan, &
       tlat, &
       tlon, &
       toclat, &
       toclon
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       indx_arr, &
       int_arr
  INTEGER :: &
       err, &
       err_verb_, &
       i, &
       iday, &
       indx, &
       istep, &
       j, &
       k, &
       l, &
       ls_niter_major_max, &
       ls_niter_major_min, &
       ls_niter_minor, &
       lu, &
       lu_orb_out, &
       min, minstep, &
       month, month0, month1, &
       nobj, &
       nobs, &
       norb, &
       noutlier, &
       nstep, &
       sor_niter, &
       sor_norb, &
       sor_norb_sw, &
       sor_ntrial, &
       sor_ntrial_sw, &
       sor_type_prm, &
       year, year0, year1
  LOGICAL, DIMENSION(:,:), POINTER :: &
       obs_masks
  LOGICAL, DIMENSION(:), POINTER :: &
       perturbers
  LOGICAL, DIMENSION(6) :: &
       ls_element_mask
  LOGICAL, DIMENSION(4) :: &
       sor_iterate_bounds
  LOGICAL :: &
       compress, &
       first, &
       gaussian_rho, &
       separately, &                        !! Output orbit(s)/ephemerides/etc separately for each object
       outlier_rejection_prm, &
       plot_open, &
       plot_results, &
       pp_H_estimation, &
       random_obs, &
       regularized, &
       uniform, &
       write_residuals

  ! Defaults:
  task = " "
  compress = .FALSE.
  orbit_format_out = "des"
  element_type_comp_prm = "keplerian"
  element_type_out_prm = "keplerian"
  err_verb = 1
  info_verb = 1
  gnuplot_scripts_dir = "."
  obs_stdev_arr_prm = -1.0_bp
  observation_format_out = "mpc3"
  obsy_code = "500"
  orb_in_fname = " "
  orb_out_fname = " "
  separately = .FALSE.
  outlier_multiplier_prm = 3.0_bp
  outlier_rejection_prm = .FALSE.
  plot_open = .FALSE.
  plot_results = .FALSE.
  pp_H_estimation = .FALSE.
  pp_G = 99.9_bp
  pp_G_unc = 99.9_bp

  IF (get_cl_option("--version",.FALSE.)) THEN
     WRITE(stdout,"(A)") ""
     WRITE(stdout,"(A)") "OpenOrb v0.9.3"
     WRITE(stdout,"(A)") "Copyright 2009 Mikael Granvik, Jenni Virtanen, Karri Muinonen,"
     WRITE(stdout,"(A)") "               Teemu Laakso, Dagmara Oszkiewicz"
     WRITE(stdout,"(A)") ""
     WRITE(stdout,"(A)") "OpenOrb comes with NO WARRANTY, to the extent permitted by law."
     WRITE(stdout,"(A)") "You may redistribute copies of OpenOrb under the terms of the "
     WRITE(stdout,"(A)") "GNU General Public License. For more information about these "
     WRITE(stdout,"(A)") "matters, see the file named COPYING or visit "
     WRITE(stdout,"(A)") "<http://www.gnu.org/licenses/>."
     WRITE(stdout,"(A)") ""
  END IF

  task = get_cl_option("--task=",task)
  IF (LEN_TRIM(task) == 0 .OR. get_cl_option("--help",.FALSE.)) THEN
     WRITE(stdout,"(A)") "Usage:" 
     WRITE(stdout,"(A)") "oorb [ --version | --help | --task=TASK ] [ options ]"
     WRITE(stdout,"(A)") " or"
     WRITE(stdout,"(A)") "oorb [ --version | --help | --task=TASK ] [ options ]"
     WRITE(stdout,"(A)") ""
     WRITE(stdout,"(A)") ""
     STOP
  END IF

  ! Set path to data files:
  CALL setAccessToDataFiles()
  CALL JPL_ephemeris_init(error, &
       filename=TRIM(OORB_DATA_DIR) // "/" // TRIM(EPH_FNAME)) 
  IF (error) THEN
     CALL errorMessage("oorb", &
          "Could not initialize planetary ephemerides using the " // &
          TRIM(OORB_DATA_DIR) // "/" // TRIM(EPH_FNAME) // " file.", 1)
     STOP
  END IF

  ! Set path to Gnuplot scripts using environment variable:
  CALL getenv("OORB_GNUPLOT_SCRIPTS_DIR", gnuplot_scripts_dir)

  ! Set path to configuration file:
  ! First, try the environment variable:
  CALL getenv("OORB_CONF", conf_fname)
  IF (LEN_TRIM(conf_fname) == 0) THEN
     ! Second, if environment variable not defined use default:
     conf_fname = "./oorb.conf"
  END IF
  ! Third, (if specified) the command-line option overrides the previous:
  conf_fname = get_cl_option("--conf=",conf_fname)
  CALL NEW(conf_file, conf_fname)
  CALL setStatusOld(conf_file)
  CALL OPEN(conf_file)
  IF (error) THEN
     CALL errorMessage("oorb", &
          "Configuration file missing or error while opening it.", 1)
     STOP
  END IF
  CALL readConfigurationFile(conf_file, &
       info_verb=info_verb, &
       err_verb=err_verb, &
       plot_results=plot_results, &
       plot_open=plot_open, &
       obs_stdev_arr=obs_stdev_arr_prm, &
       outlier_rejection=outlier_rejection_prm, &
       outlier_multiplier=outlier_multiplier_prm, &
       element_type_comp=element_type_comp_prm, &
       element_type_out=element_type_out_prm, &
       observation_format_out=observation_format_out, &
       orbit_format_out=orbit_format_out, &
       pp_H_estimation=pp_H_estimation, &
       pp_G=pp_G, &
       pp_G_unc=pp_G_unc)
  IF (error) THEN
     CALL errorMessage("oorb", &
          "TRACE BACK (15)", 1)
     STOP
  END IF
  IF (.NOT.ALL(obs_stdev_arr_prm < 0.0_bp)) THEN
     WHERE (obs_stdev_arr_prm < 0.0_bp)
        obs_stdev_arr_prm = 0.0_bp
     END WHERE
  END IF

  ! Read observation file if given:
  obs_fname = get_cl_option("--obs-in="," ")
  IF (LEN_TRIM(obs_fname) /= 0) THEN
     CALL NEW(obs_file, TRIM(obs_fname))
     CALL setStatusOld(obs_file)
     CALL OPEN(obs_file)
     IF (error) THEN
        CALL errorMessage("oorb", &
             "TRACE BACK (30)", 1)
        STOP
     END IF
     IF (ANY(obs_stdev_arr_prm < 0.0_bp)) THEN
        CALL NEW(obss_in, obs_file)
     ELSE
        CALL NEW(obss_in, obs_file, stdev=obs_stdev_arr_prm)
     END IF
     IF (error) THEN
        CALL errorMessage("oorb", &
             "TRACE BACK (35)", 1)
        STOP
     END IF
     CALL NULLIFY(obs_file)
     indx = INDEX(obs_fname,".",back=.TRUE.)
     out_fname = obs_fname(1:indx-1)
  ELSE IF (task == "tompc3" .OR. &
       task == "tompc" .OR. &
       task == "ranging" .OR. &
       task == "lsl") THEN
     CALL errorMessage("oorb", &
          "Path to observation file must be supplied using " // &
          "the '--obs-in=FILE' option.", 1)
     STOP
  END IF

  ! Read orbit file if given:
  orb_in_fname = get_cl_option("--orb-in=",orb_in_fname)
  IF (LEN_TRIM(orb_in_fname) /= 0) THEN

     ! Open orbit file
     CALL NEW(orb_in_file,TRIM(orb_in_fname))
     IF (error) THEN
        CALL errorMessage("oorb", &
             "TRACE BACK (40)", 1)
        STOP
     END IF
     CALL setStatusOld(orb_in_file)
     CALL OPEN(orb_in_file)
     IF (error) THEN
        CALL errorMessage("oorb", &
             "TRACE BACK (45)", 1)
        STOP
     END IF

     ! Decide which routine should read the input 
     ! file by checking the suffix:
     indx = INDEX(orb_in_fname, ".", back=.TRUE.)
     IF (indx == 0) THEN
        error = .TRUE.
        CALL errorMessage("oorb", &
             "Suffix (either .orb or .des) required in orbit file name.",1)
        STOP
     END IF
     suffix = orb_in_fname(indx+1:LEN_TRIM(orb_in_fname))
     CALL locase(suffix, error)
     IF (error) THEN
        CALL errorMessage("oorb", &
             "The orbit file suffix contains forbidden characters.", 1)
        STOP
     END IF

     SELECT CASE (suffix)
     CASE ("orb")
        ! OpenOrb orbit file

        norb = getNrOfLines(orb_in_file) - 4 ! 4 header lines are not taken into count
        IF (error) THEN
           CALL errorMessage("oorb", &
                "TRACE BACK (50)", 1)
           STOP
        END IF
        ALLOCATE(id_arr_in(norb), orb_arr_in(norb), &
             element_type_pdf_arr_in(norb), cov_arr_in(norb,6,6), &
             HG_arr_in(norb,4), pdf_arr_in(norb), &
             rchi2_arr_in(norb), jac_arr_in(norb,3), &
             reg_apr_arr_in(norb), stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb", &
                "Could not allocate memory (20).", 1)
           STOP
        END IF
        jac_arr_in = -1.0_bp
        pdf_arr_in = -1.0_bp
        cov_arr_in = -1.0_bp
        id_arr_in = " "
        element_type_pdf_arr_in = " "
        header(1:4)(:) = " "
        FORALL (i=1:norb)
           HG_arr_in(i,1:4) = 99.9_bp
        END FORALL
        DO i=1,norb
           CALL readOpenOrbOrbitFile(getUnit(orb_in_file), header, &
                element_type_in=element_type_in, id=id_arr_in(i), &
                orb=orb_arr_in(i), &
                element_type_pdf=element_type_pdf_arr_in(i), &
                cov=cov_arr_in(i,:,:), H=HG_arr_in(i,1), &
                G=HG_arr_in(i,3), pdf=pdf_arr_in(i), &
                rchi2=rchi2_arr_in(i), reg_apr=reg_apr_arr_in(i), &
                jac_sph_inv=jac_arr_in(i,1), &
                jac_car_kep=jac_arr_in(i,2), &
                jac_equ_kep=jac_arr_in(i,3))
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "Could not read orbit file.", 1)
              STOP
           END IF
        END DO
        CALL NULLIFY(orb_in_file)

        ! Calculate the number of different objects in the orbit file:
        ALLOCATE(indx_arr(SIZE(id_arr_in)))
        CALL quickSort(id_arr_in, indx_arr, error)
        nobj = 1
        DO i=1,SIZE(id_arr_in)-1
           IF (id_arr_in(indx_arr(i)) /= id_arr_in(indx_arr(i+1))) THEN
              nobj = nobj + 1
           END IF
        END DO

        ! Initialize stochasticorbits if uncertainty information available:
        IF (nobj == 1 .AND. norb > 1 .AND. &
             ALL(pdf_arr_in > 0.0_bp) .AND. &
             ALL(jac_arr_in > 0.0_bp)) THEN
           ALLOCATE(storb_arr_in(1))
           CALL NEW(storb_arr_in(1), orb_arr_in, pdf_arr_in, &
                element_type_pdf_arr_in(1), jac_arr=jac_arr_in, &
                reg_apr_arr=reg_apr_arr_in, &
                rchi2_arr=rchi2_arr_in)
           id = id_arr_in(1)
           DEALLOCATE(id_arr_in)
           ALLOCATE(id_arr_in(1))
           id_arr_in(1) = id
        ELSE IF (nobj == norb .AND. &
             ALL(cov_arr_in(:,1,1) > 0.0_bp)) THEN
           ALLOCATE(storb_arr_in(norb))
           DO i=1,norb
              CALL NEW(storb_arr_in(i), orb_arr_in(i), cov_arr_in(i,:,:), &
                   cov_type=element_type_in, element_type=element_type_in)            
              !WRITE(*,*) "cond nr", cond_nr(cov_arr_in(i,:,:), error)
           END DO
        ELSE IF (nobj > 1 .AND. nobj /= norb) THEN
           CALL errorMessage("oorb", &
                "Input of more than one sampled orbital-element pdf not supported.", 1)
           STOP
        END IF

     CASE ("des")
        ! Data exchange format

        norb = getNrOfLines(orb_in_file) 
        IF (error) THEN
           CALL errorMessage("oorb", &
                "TRACE BACK (55)", 1)
           STOP
        END IF
        ALLOCATE(id_arr_in(norb), orb_arr_in(norb), &
             HG_arr_in(norb,4), stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb", &
                "Could not allocate memory (20).", 1)
           STOP
        END IF
        id_arr_in = " "
        header = " "
        CALL readDESOrbitFile(getUnit(orb_in_file), norb, header(1), &
             id_arr_in, orb_arr_in, HG_arr_in(:,1))
        IF (error) THEN
           CALL errorMessage("oorb", &
                "Could not read orbit file.", 1)
           STOP
        END IF
        CALL NULLIFY(orb_in_file)
        HG_arr_in(:,3) = 0.15_bp
        HG_arr_in(:,2) = 99.9_bp
        HG_arr_in(:,4) = 99.9_bp
        id_arr_in => reallocate(id_arr_in, norb)
        orb_arr_in => reallocate(orb_arr_in, norb)
        HG_arr_in => reallocate(HG_arr_in, norb, 4)

     CASE default

        ! Not a valid suffix.
        CALL errorMessage("oorb", &
             "'." // TRIM(suffix) // &
             "' is not a valid suffix for orbit files. " // &
             "Use either '.orb' or '.des'.", 1)
        STOP

     END SELECT

     ! Convert to the desired element type if needed:
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (element_type_in == "cartesian" .AND. &
                element_type_comp_prm == "keplerian") THEN
              CALL toKeplerian(storb_arr_in(i))
           ELSE IF (element_type_in == "keplerian" .AND. &
                element_type_comp_prm == "cartesian") THEN
              CALL toCartesian(storb_arr_in(i), "ecliptic")
           ELSE IF (element_type_comp_prm /= "cartesian" .AND. &
                element_type_comp_prm /= "keplerian") THEN
              CALL errorMessage("oorb", &
                   "No such option: " // TRIM(element_type_comp_prm), 1)
              STOP
           END IF
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "TRACE BACK (60)", 1)
              STOP
           END IF
        END DO
        ! Delete 'raw' orbit data
        DO i=1,norb
           CALL NULLIFY(orb_arr_in(i))
        END DO
        DEALLOCATE(orb_arr_in, pdf_arr_in, rchi2_arr_in, jac_arr_in, &
             reg_apr_arr_in, element_type_pdf_arr_in)
     ELSE
        DO i=1,norb
           IF (element_type_comp_prm == "keplerian") THEN
              CALL toKeplerian(orb_arr_in(i))
           ELSE IF (element_type_comp_prm == "cartesian") THEN
              CALL toCartesian(orb_arr_in(i), "ecliptic")
           ELSE IF (element_type_comp_prm == "cometary") THEN
              CALL toCometary(orb_arr_in(i))
           END IF
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
        END DO
     END IF

  ELSE IF (task == "lsl" .OR. &
       task == "ephemeris" .OR. &
       task == "propagation") THEN
     CALL errorMessage("oorb", &
          "Path to orbit file must be supplied using " // &
          "the '--orb-in=FILE' option.", 1)
     STOP
  END IF

  ! Choose logical unit for output orbit file (all orbits in same file)
  orb_out_fname = get_cl_option("--orb-out=",orb_out_fname)
  IF (LEN_TRIM(orb_out_fname) == 0) THEN
     lu_orb_out = stdout
  ELSE
     CALL NEW(orb_out_file,TRIM(orb_out_fname))
     CALL OPEN(orb_out_file)
     IF (error) THEN
        CALL errorMessage("oorb", &
             "TRACE BACK (70)", 1)
        STOP
     END IF
     lu_orb_out = getUnit(orb_out_file)
  END IF

  ! Output orbit(s)/ephemerides/etc written into separate files for
  ! each object (overrides --orb-out)
  separately = get_cl_option("--separately",separately)

  SELECT CASE (task)

  CASE ("none")

     CONTINUE

  CASE ("toorbpdf")

     ! Convert orbits + uncertainties
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (containsSampledPDF(storb_arr_in(i))) THEN
              orb_arr_in => getSampleOrbits(storb_arr_in(i))
              pdf_arr_in => getPDFValues(storb_arr_in(i), element_type_out_prm)
              !pdf_arr_in = pdf_arr_in/sum(pdf_arr_in)
              DO j=1,SIZE(orb_arr_in)
                 CALL writeOpenOrbOrbitFile(lu_orb_out, i==1.AND.j==1, &
                      element_type_out_prm, id_arr_in(i), &
                      orb_arr_in(j), pdf=pdf_arr_in(j), &
                      element_type_pdf=element_type_out_prm, &
                      H=HG_arr_in(i,1), &
                      G=HG_arr_in(i,3))
                 CALL NULLIFY(orb_arr_in(j))
              END DO
              DEALLOCATE(orb_arr_in, pdf_arr_in)
           ELSE
              orb = getNominalOrbit(storb_arr_in(i))
              cov = getCovarianceMatrix(storb_arr_in(i), element_type_out_prm)
              CALL writeOpenOrbOrbitFile(lu_orb_out, i==1, &
                   element_type_out_prm, id_arr_in(i), &
                   orb=orb, cov=cov, H=HG_arr_in(i,1), &
                   G=HG_arr_in(i,3))
           END IF
        END DO
     END IF

  CASE ("toorb")

     norb = HUGE(norb)
     norb = get_cl_option("--norb=",norb)
     H_max = get_cl_option("--H-max=",HUGE(H_max))

     ! Convert orbits
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)        
           orb_arr_in => getSampleOrbits(storb_arr_in(i))
           k = 0
           DO j=1,SIZE(orb_arr_in)
              IF (get_cl_option("--only-prograde",.FALSE.)) THEN
                 elements = getElements(orb_arr_in(j), "keplerian")
                 IF (elements(3) > pi/2) THEN
                    ! retrograde orbit -> skip it
                    CYCLE
                 END IF
              END IF
              IF (get_cl_option("--H-max=",.FALSE.)) THEN
                 IF (HG_arr_in(i,1) > H_max) THEN
                    ! too large H magnitude -> skip it
                    CYCLE
                 END IF
              END IF
              k = k + 1
              CALL writeOpenOrbOrbitFile(lu_orb_out, i==1 &
                   .AND. (j==1 .OR. k==1), element_type_out_prm, id_arr_in(i), &
                   orb_arr_in(j), H=HG_arr_in(i,1), G=HG_arr_in(i,3))
              IF (k ==  norb) THEN
                 EXIT
              END IF
           END DO
        END DO
     ELSE IF (ASSOCIATED(orb_arr_in)) THEN
        j = 0
        DO i=1,SIZE(orb_arr_in)
           IF (get_cl_option("--only-prograde",.FALSE.)) THEN
              elements = getElements(orb_arr_in(i), "keplerian")
              IF (elements(3) > pi/2) THEN
                 ! retrograde orbit -> skip it
                 CYCLE
              END IF
           END IF
           IF (get_cl_option("--H-max=",.FALSE.)) THEN
              IF (HG_arr_in(i,1) > H_max) THEN
                 ! too large H magnitude -> skip it
                 CYCLE
              END IF
           END IF
           j = j + 1
           CALL writeOpenOrbOrbitFile(lu_orb_out, (i==1 .OR. j==1), &
                element_type_out_prm, id_arr_in(i), orb_arr_in(i), &
                H=HG_arr_in(i,1), G=HG_arr_in(i,3))
           IF (j == norb) THEN
              EXIT
           END IF
        END DO
     END IF

  CASE ("orbitstodes")

     H_max = get_cl_option("--H-max=",HUGE(H_max))
     i_min = get_cl_option("--i-min=",0.0_bp)
     i_max = get_cl_option("--i-max=",HUGE(i_max))

     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (containsSampledPDF(storb_arr_in(i))) THEN
              orb_arr_in => getSampleOrbits(storb_arr_in(i))
              DO j=1,SIZE(orb_arr_in)
                 IF (get_cl_option("--H-max=",.FALSE.)) THEN
                    IF (HG_arr_in(i,1) > H_max) THEN
                       ! too large H magnitude -> skip it
                       CYCLE
                    END IF
                 END IF
                 IF (get_cl_option("--only-neos",.FALSE.) .OR. &
                      get_cl_option("--i-min=",.FALSE.) .OR. &
                      get_cl_option("--i-max=",.FALSE.)) THEN
                    elements = getElements(orb_arr_in(j),"cometary")
                 END IF
                 IF (get_cl_option("--only-neos",.FALSE.)) THEN
                    IF (elements(1) > 1.3_bp) THEN
                       ! too large q -> skip it
                       CYCLE
                    END IF
                 END IF
                 IF (get_cl_option("--i-min=",.FALSE.)) THEN
                    IF (elements(3)/rad_deg < i_min) THEN
                       ! inclination too small -> skip it
                       CYCLE
                    END IF
                 END IF
                 IF (get_cl_option("--i-max=",.FALSE.)) THEN
                    IF (elements(3)/rad_deg > i_max) THEN
                       ! inclination too large -> skip it
                       CYCLE
                    END IF
                 END IF
                 CALL writeDESOrbitFile(lu_orb_out, i==1 .AND. j==1, "cometary", &
                      id_arr_in(i), orb_arr_in(j), HG_arr_in(i,1))
                 IF (error) THEN
                    CALL errorMessage("oorb", &
                         "DES output failed at orbit:", 1)
                    WRITE(stderr,*) j
                    STOP
                 END IF
                 CALL NULLIFY(orb_arr_in(j))
              END DO
              DEALLOCATE(orb_arr_in)
           ELSE
              IF (HG_arr_in(i,1) > H_max) THEN
                 ! too large H magnitude -> skip it
                 CYCLE
              END IF
              orb = getNominalOrbit(storb_arr_in(i))
              IF (get_cl_option("--only-neos",.FALSE.)) THEN
                 elements = getElements(orb,"cometary")
                 IF (elements(1) > 1.3_bp) THEN
                    ! too large q -> skip it
                    CYCLE
                 END IF
              END IF
              CALL writeDESOrbitFile(lu_orb_out, i==1, "cometary", &
                   id_arr_in(i), orb, HG_arr_in(i,1))              
              CALL NULLIFY(orb)
           END IF
        END DO
        STOP
     END IF
     IF (ASSOCIATED(orb_arr_in)) THEN
        DO i=1,SIZE(orb_arr_in)
           IF (get_cl_option("--H-max=",.FALSE.)) THEN
              IF (HG_arr_in(i,1) > H_max) THEN
                 ! too large H magnitude -> skip it
                 CYCLE
              END IF
           END IF
           IF (get_cl_option("--only-neos",.FALSE.)) THEN
              elements = getElements(orb_arr_in(i),"cometary")
              IF (elements(1) > 1.3_bp) THEN
                 ! too large q -> skip it
                 CYCLE
              END IF
           END IF
           CALL writeDESOrbitFile(lu_orb_out, i==1, "cometary", &
                id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1))
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "DES output failed at orbit:", 1)
              WRITE(stderr,*) i
              STOP
           END IF
        END DO
     END IF

  CASE ("observationstodes")

     ! Output input observations in DES format.

     obs_fname = " "
     obs_fname = get_cl_option("--obs-out=",obs_fname)
     IF (LEN_TRIM(obs_fname) == 0) THEN
        lu = stdout
     ELSE
        CALL NEW(obs_file,TRIM(obs_fname))
        CALL OPEN(obs_file)
        IF (error) THEN
           CALL errorMessage("oorb / tompc3", &
                "TRACE BACK", 1)
           STOP
        END IF
        lu = getUnit(obs_file)
     END IF
     CALL writeObservationFile(obss_in, lu, "des")
     IF (LEN_TRIM(obs_fname) /= 0) THEN
        CALL NULLIFY(obs_file)
     END IF

  CASE ("tompc3")

     ! Output input observations in new MPC format.

     obs_fname = " "
     obs_fname = get_cl_option("--obs-out=",obs_fname)
     IF (LEN_TRIM(obs_fname) == 0) THEN
        lu = stdout
     ELSE
        CALL NEW(obs_file,TRIM(obs_fname))
        CALL OPEN(obs_file)
        IF (error) THEN
           CALL errorMessage("oorb / tompc3", &
                "TRACE BACK", 1)
           STOP
        END IF
        lu = getUnit(obs_file)
     END IF
     CALL writeObservationFile(obss_in, lu, "mpc3")
     IF (LEN_TRIM(obs_fname) /= 0) THEN
        CALL NULLIFY(obs_file)
     END IF

  CASE ("tompc")

     ! Output input observations in current MPC format.

     obs_fname = " "
     obs_fname = get_cl_option("--obs-out=",obs_fname)
     IF (LEN_TRIM(obs_fname) == 0) THEN
        lu = stdout
     ELSE
        CALL NEW(obs_file,TRIM(obs_fname))
        CALL OPEN(obs_file)
        IF (error) THEN
           CALL errorMessage("oorb / tompc", &
                "TRACE BACK", 1)
           STOP
        END IF
        lu = getUnit(obs_file)
     END IF
     CALL writeObservationFile(obss_in, lu, "mpc")
     IF (LEN_TRIM(obs_fname) /= 0) THEN
        CALL NULLIFY(obs_file)
     END IF

  CASE ("astorbtoorb")

     ! Convert Ted's astorb file to the .orb format:

     orb_in_fname = " "
     orb_in_fname = get_cl_option("--astorb=",orb_in_fname)
     IF (LEN_TRIM(orb_in_fname) == 0) THEN
        CALL errorMessage("oorb / astorbtoorb", &
             "ASTORB file not given.", 1)
        STOP        
     END IF
     CALL NEW(orb_in_file, TRIM(orb_in_fname))
     CALL setStatusOld(orb_in_file)
     CALL OPEN(orb_in_file)
     IF (error) THEN
        CALL errorMessage("oorb / astorbtoorb", &
             "TRACE BACK (5)", 1)
        STOP        
     END IF
     norb = getNrOfLines(orb_in_file)
     ! Format for astorb.dat:
     frmt = '(A6,1X,A18,1X,A15,1X,F5.2,1X,F5.2,1X,A4,1X,A5,1X,A4,' // &
          '1X,6I4,1X,2I5,1X,I4,2I2.2,3(1X,F10.6),F10.6,1X,F10.8,1X,' // &
          'F12.8,1X,I4,2I2.2,1X,F7.2,1X,F8.2,1X,I4,2I2,' // &
          '3(1X,F7.2,1X,I4,2I2))'
     i = 0
     ALLOCATE(str_arr(6), real_arr(5), int_arr(23))
     str_arr(1)(1:LEN(str_arr(1))) = " "
     str_arr(2)(1:LEN(str_arr(2))) = " "
     DO
        str_arr(1)(1:6) = " "
        str_arr(2)(1:18) = " "
        READ(getUnit(orb_in_file),TRIM(frmt),iostat=err) &
             str_arr(1)(1:6), str_arr(2)(1:18), str_arr(3)(1:15), &
             H_value, G_value, str_arr(4)(1:4), str_arr(5)(1:5), &
             str_arr(6)(1:4), int_arr(1:8), year, month, iday, &
             elements(6:1:-1), int_arr(9:11), real_arr(1:2), &
             int_arr(12:14), real_arr(3), int_arr(15:17), &
             real_arr(4), int_arr(18:20), real_arr(5), int_arr(21:23)
        IF (err > 0) THEN
           CALL errorMessage("oorb / astorbtoorb", &
                'Error while reading line from file.',1)
           STOP
        ELSE IF (err < 0) THEN
           ! Input file read.
           EXIT
        END IF
!!$        ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!!$        ! Discard some orbits depending on the following requirement:
!!$        IF (int_arr(7) < 90 .OR. & ! orbital arc in days
!!$                                !             elements(1)*(1.0_bp-elements(2)) < 1.666_bp .OR. &
!!$                                !             elements(1) <= 2.0_bp .or. elements(1) >= 3.5_bp) then ! .or. &
!!$             elements(1)*(1.0_bp-elements(2)) < 1.3_bp .OR. &
!!$             elements(1) > 4.5_bp .OR. &
!!$             elements(3) > 60.0_bp .OR. &
!!$             H_value > 15.5_bp) THEN
!!$           CYCLE
!!$        END IF
!!$        ! End requirement
!!$        ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        IF (LEN_TRIM(str_arr(2)) /= 0) THEN
           CALL encodeMPC3Designation(str_arr(2))
           IF (error) THEN
              CALL errorMessage("oorb / astorbtoorb", &
                   "TRACE BACK (10)",1)
              STOP
           END IF
        END IF
        DO j=1,6
           IF (IACHAR(str_arr(1)(j:j)) == 0) THEN
              str_arr(1)(j:j) = CHAR(32)
           END IF
        END DO
        DO j=1,18
           IF (IACHAR(str_arr(2)(j:j)) == 0) THEN
              str_arr(2)(j:j) = CHAR(32)
           END IF
        END DO
        i = i + 1
        CALL NEW(epoch, year, month, REAL(iday,bp), "TT")
        elements(3:6) = elements(3:6)*rad_deg
        CALL NEW(orb, elements, "keplerian", "ecliptic", epoch)
        IF (error) THEN
           CALL errorMessage("oorb / astorbtoorb", &
                "TRACE BACK (15)",1)
           STOP
        END IF
        CALL NULLIFY(epoch)
        id(1:LEN(id)) = " "
        IF (LEN_TRIM(str_arr(1)) /= 0) THEN
           id = TRIM(ADJUSTL(str_arr(1)))
           DO WHILE (LEN_TRIM(id) < 7)
              id = '0' // TRIM(id)
           END DO
        ELSE
           id = TRIM(ADJUSTL(str_arr(2)))
        END IF
        IF (orbit_format_out == "des") THEN
           CALL writeDESOrbitFile(lu_orb_out, i==1, "cometary", &
                id, orb, H_value, 1, 6, &
                -1.0_bp, "OPENORB")
        ELSE IF (orbit_format_out == "orb") THEN
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id), orb=orb, &
                H=H_value, G=G_value)
        END IF
        CALL NULLIFY(orb)
     END DO
     CALL NULLIFY(orb_in_file)

  CASE ("mpcorbtoorb")

     ! Convert MPC's MPCORB.DAT file to the .orb format:

     orb_in_fname = " "
     orb_in_fname = get_cl_option("--mpcorb=",orb_in_fname)
     IF (LEN_TRIM(orb_in_fname) == 0) THEN
        CALL errorMessage("oorb / mpcorbtoorb", &
             "MPCORB file not given.", 1)
        STOP        
     END IF
     CALL NEW(orb_in_file, TRIM(orb_in_fname))
     CALL setStatusOld(orb_in_file)
     CALL OPEN(orb_in_file)
     IF (error) THEN
        CALL errorMessage("oorb / mpcorbtoorb", &
             "TRACE BACK (5)", 1)
        STOP        
     END IF
     norb = getNrOfLines(orb_in_file)
     DEALLOCATE(id_arr_in, orb_arr_in, HG_arr_in, arc_arr, stat=err)
     ALLOCATE(id_arr_in(norb), orb_arr_in(norb), HG_arr_in(norb,2), arc_arr(norb))
     CALL readMPCOrbitFile(getUnit(orb_in_file), norb, id_arr_in, orb_arr_in, HG_arr_in, arc_arr)
     IF (error) THEN
        CALL errorMessage("oorb / mpcorbtoorb", &
             "TRACE BACK (10)", 1)
        STOP        
     END IF
     CALL NULLIFY(orb_in_file)

     DO i=1,norb
        elements = getElements(orb_arr_in(i), "keplerian")
        IF (orbit_format_out == "des") THEN
           CALL writeDESOrbitFile(lu_orb_out, i==1, "cometary", &
                id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), 1, 6, &
                -1.0_bp, "OPENORB")
        ELSE IF (orbit_format_out == "orb") THEN
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
                H=HG_arr_in(i,1), G=HG_arr_in(i,2))
        END IF
        IF (error) THEN
           CALL errorMessage("oorb / mpcorbtoorb", &
                "TRACE BACK (15)", 1)
           STOP        
        END IF
     END DO


  CASE ("ranging")

     ! Orbital inversion using statistical orbital ranging, that is,
     ! without making any assumptions on the shape of the resulting
     ! orbital-element pdf.

     CALL NULLIFY(epoch)
     dyn_model = " "
     integrator = " "
     integration_step = -1.0_bp
     apriori_a_max = -1.0_bp
     apriori_a_min = -1.0_bp
     apriori_periapsis_max = -1.0_bp
     apriori_periapsis_min = -1.0_bp
     apriori_apoapsis_max = -1.0_bp
     apriori_apoapsis_min = -1.0_bp
     apriori_rho_min = -1.0_bp
     sor_type_prm = -1
     sor_2point_method = " "
     sor_2point_method_sw = " "
     sor_norb = -1
     sor_norb_sw = -1
     sor_ntrial = -1
     sor_ntrial_sw = -1
     sor_niter = -1
     sor_rho_init = HUGE(sor_rho_init)
     sor_genwin_multiplier = -1.0_bp
     sor_genwin_offset = -1.0_bp
     sor_iterate_bounds = .TRUE.
     accwin_multiplier = -1.0_bp
     gaussian_rho = .FALSE.
     regularized = .FALSE.
     uniform = .FALSE.
     write_residuals = .FALSE.
     pdf_ml_init = -1.0_bp
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          integrator=integrator, &
          integration_step=integration_step, &
          apriori_a_max=apriori_a_max, &
          apriori_a_min=apriori_a_min, &
          apriori_periapsis_max=apriori_periapsis_max, &
          apriori_periapsis_min=apriori_periapsis_min, &
          apriori_apoapsis_max=apriori_apoapsis_max, &
          apriori_apoapsis_min=apriori_apoapsis_min, &
          apriori_rho_min=apriori_rho_min, &
          sor_type=sor_type_prm, sor_2point_method=sor_2point_method, &
          sor_2point_method_sw=sor_2point_method_sw, &
          sor_norb=sor_norb, sor_norb_sw=sor_norb_sw, &
          sor_ntrial=sor_ntrial, sor_ntrial_sw=sor_ntrial_sw, &
          sor_niter=sor_niter, sor_rho_init=sor_rho_init, &
          sor_genwin_multiplier=sor_genwin_multiplier, &
          sor_genwin_offset=sor_genwin_offset, &
          sor_iterate_bounds=sor_iterate_bounds, &
          accwin_multiplier=accwin_multiplier, &
          uniform_pdf=uniform, regularized_pdf=regularized, &
          sor_random_obs=random_obs, sor_rho_gauss=gaussian_rho, &
          write_residuals=write_residuals, &
          pdf_ml=pdf_ml_init)

     IF (error) THEN
        CALL errorMessage("oorb / ranging", &
             "TRACE BACK (5)", 1)
        STOP
     END IF

     obss_sep => getSeparatedSets(obss_in)
     IF (error) THEN
        CALL errorMessage("oorb / ranging", &
             "TRACE BACK (10)", 1)
        STOP
     END IF
     CALL NULLIFY(obss_in)
     DO i=1,SIZE(obss_sep)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (25)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(2(1X,A))") "Object:", TRIM(id)
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (sor_rho_init(3) > HUGE(sor_rho_init(3))/2) THEN
           ! Initialize the rho2-rho1 range based on the observational timespan.
           IF (dt > 10) THEN
              sor_rho_init(3) = -0.5_bp
           ELSE
              sor_rho_init(3) = -0.05_bp*dt
           END IF
           sor_rho_init(4) = -sor_rho_init(3)
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (35)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (40)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t, mjd, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL NEW(storb, obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (55)", 1)
           STOP
        END IF
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             regularized_pdf = regularized, &
             uniform_pdf = uniform, &
             pdf_ml=pdf_ml_init, &
             accept_multiplier=accwin_multiplier, &
             apriori_a_max=apriori_a_max, apriori_a_min=apriori_a_min, &
             apriori_periapsis_max=apriori_periapsis_max, &
             apriori_periapsis_min=apriori_periapsis_min, &
             apriori_apoapsis_max=apriori_apoapsis_max, &
             apriori_apoapsis_min=apriori_apoapsis_min, &
             apriori_rho_min=apriori_rho_min, &
             sor_2point_method=sor_2point_method, &
             sor_2point_method_sw=sor_2point_method_sw, &
             sor_norb=sor_norb, sor_ntrial=sor_ntrial, &
             sor_rho1_l=sor_rho_init(1), sor_rho1_u=sor_rho_init(2), &
             sor_rho2_l=sor_rho_init(3), sor_rho2_u=sor_rho_init(4), &
             sor_iterate_bounds=sor_iterate_bounds, &
             sor_random_obs_selection=.FALSE., &
             gaussian_pdf=gaussian_rho, &
             sor_generat_multiplier=sor_genwin_multiplier, &
             sor_generat_offset=sor_genwin_offset)
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (60)", 1)
           STOP
        END IF

        SELECT CASE (sor_type_prm)
        CASE (1)
           CALL statisticalRanging(storb)
        CASE (2)
           CALL setParameters(storb, &
                sor_niter=sor_niter)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (85)", 1)
              STOP
           END IF
           CALL autoStatisticalRanging(storb)
        CASE (3)
           CALL setParameters(storb, &
                sor_norb_sw=sor_norb_sw, &
                sor_ntrial_sw=sor_ntrial_sw, &
                sor_niter=sor_niter)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (90)", 1)
              STOP
           END IF
           CALL stepwiseRanging(storb, nobs_max=-1)
        CASE default
           CALL errorMessage("oorb / ranging", &
                "Unknown type of ranging:", 1)
           IF (err_verb >= 1) THEN
              WRITE(stderr,"(2X,I0)") sor_type_prm
           END IF
           STOP              
        END SELECT

        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (115)", 1)
           error  = .FALSE.
           CALL NEW(out_file, "problematic_observation_sets" // "." // &
                TRIM(observation_format_out))
           CALL setPositionAppend(out_file)
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (130)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (135)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           CALL NULLIFY(storb)
        ELSE
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Ranging inversion for object", &
                   TRIM(id), "is ready."
           END IF

           IF (pp_H_estimation) THEN
              CALL NEW(physparam, storb)
              IF (pp_G > 99.0_bp) THEN
                 CALL estimateHAndG(physparam, obss_sep(i))
              ELSE
                 CALL estimateHAndG(physparam, obss_sep(i), &
                      input_G=pp_G, input_delta_G=pp_G_unc)
              END IF
              HG_arr_in => getH0Distribution(physparam)
              ALLOCATE(temp_arr(SIZE(HG_arr_in,dim=1),4))
              temp_arr(:,1:2) = HG_arr_in(:,1:2)
              DEALLOCATE(HG_arr_in)
              HG_arr_in => getGDistribution(physparam)
              temp_arr(:,3:4) = HG_arr_in(:,1:2)
              ALLOCATE(HG_arr_in(SIZE(temp_arr,dim=1),4))
              HG_arr_in = temp_arr
              DEALLOCATE(temp_arr)              
              CALL NULLIFY(physparam)
           END IF

           CALL NEW(out_file, TRIM(id) // ".sor")
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (150)", 1)
              STOP
           END IF
           ! WRITE OBSERVATIONS:
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (155)", 1)
              STOP
           END IF
           WRITE(getUnit(out_file),"(A)") "#" 
           ! WRITE RANGING PARAMETERS
           CALL writeSORResults(storb, obss_sep(i), getUnit(out_file))
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (160)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file) 
           CALL getResults(storb, sor_rho_cmp=sor_rho_cmp)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (170)", 1)
              STOP
           END IF
           ! WRITE ORBITAL-ELEMENT PDF
           orb_arr_cmp => getSampleOrbits(storb)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (175)", 1)
              STOP
           END IF
           pdf_arr_cmp => getPDFValues(storb)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (180)", 1)
              STOP
           END IF
           rchi2_arr_cmp => getReducedChi2Distribution(storb)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (185)", 1)
              STOP
           END IF
           CALL getResults(storb, &
                reg_apr_arr=reg_apr_arr_cmp, &
                jac_arr=jac_arr_cmp)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (190)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NEW(out_file, TRIM(id) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (200)", 1)
                 STOP
              END IF
              lu_orb_out = getUnit(out_file)
           END IF
           DO j=1,SIZE(orb_arr_cmp,dim=1)
              IF (orbit_format_out == "orb") THEN
                 CALL writeOpenOrbOrbitFile(lu_orb_out, &
                      print_header=j==1, &
                      element_type_out=element_type_out_prm, &
                      id=id, &
                      orb=orb_arr_cmp(j), &
                      element_type_pdf=element_type_comp_prm, &
                      pdf=pdf_arr_cmp(j), &
                      rchi2=rchi2_arr_cmp(j), &
                      reg_apr=reg_apr_arr_cmp(j), &
                      jac_sph_inv=jac_arr_cmp(j,1), &
                      jac_car_kep=jac_arr_cmp(j,2), &
                      jac_equ_kep=jac_arr_cmp(j,3), &
                      H=HG_arr_in(j,1), &
                      G=HG_arr_in(j,3))
              ELSE IF (orbit_format_out == "des") THEN
                 CALL errorMessage("oorb / ranging ", &
                      "DES format not yet supported for Ranging output.", 1)
                 STOP                 
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (205)", 1)
                 STOP
              END IF
           END DO
           IF (separately) THEN
              CALL NULLIFY(out_file)
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(id) // ".orb")
              END IF
           END IF
           ! WRITE RESIDUALS
           IF (write_residuals) THEN
              CALL NEW(out_file, TRIM(out_fname) // ".res")
              CALL setPositionAppend(out_file)
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (225)", 1)
                 STOP
              END IF
              CALL writeResiduals(storb, obss_sep(i), getUnit(out_file))
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (230)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
           END IF
           IF (plot_results) THEN
              CALL toString(dt, str, error, frmt="(F10.2)")
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (235)", 1)
                 STOP
              END IF
              str = TRIM(id) // "_"// TRIM(str)
              CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
                   "_sor_residual_stamps.eps")
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (240)", 1)
                 STOP
              END IF
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_sor_residual_stamps.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_sor_residual_stamps.eps* &")
              END IF
              ALLOCATE(elements_arr(SIZE(orb_arr_cmp,dim=1),7), stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / ranging", &
                      "Could not allocate memory (3).", 1)
                 STOP
              END IF
              CALL NEW(tmp_file, TRIM(str)// "_sor_orbits.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (250)", 1)
                 STOP
              END IF
              DO j=1,SIZE(orb_arr_cmp,dim=1)
                 IF (element_type_comp_prm == "cartesian") THEN
                    CALL rotateToEcliptic(orb_arr_cmp(j))
                 END IF
                 elements_arr(j,1:6) = getElements(orb_arr_cmp(j), element_type_comp_prm)
                 IF (error) THEN
                    CALL errorMessage("oorb / ranging", &
                         "TRACE BACK (255)", 1)
                    STOP
                 END IF
                 IF (element_type_comp_prm == "keplerian") THEN
                    elements_arr(j,3:6) = elements_arr(j,3:6)/rad_deg
                 END IF
                 elements_arr(j,7) = pdf_arr_cmp(j)
                 t = getTime(orb_arr_cmp(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / ranging ", &
                         "TRACE BACK (260)", 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmp_file),"(7(E23.15,1X),A)") &
                      elements_arr(j,1:6), &
                      pdf_arr_cmp(j), &
                      getCalendarDateString(t,"tdt")
                 IF (error) THEN
                    CALL errorMessage("oorb / ranging ", &
                         "TRACE BACK (265)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(t)
              END DO
              CALL NULLIFY(tmp_file)
              CALL NEW(tmp_file, TRIM(str) // &
                   "_sor_sample_standard_deviations.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (275)", 1)
                 STOP
              END IF
              WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                   advance="no") &
                   getObservationalTimespan(obss_sep(i))
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (280)", 1)
                 STOP
              END IF
              DO j=1,6
                 CALL moments(elements_arr(:,j), &
                      pdf=pdf_arr_cmp, std_dev=stdev, error=error)
                 WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                      advance="no") stdev
              END DO
              WRITE(getUnit(tmp_file),*)
              CALL NULLIFY(tmp_file)
              DEALLOCATE(elements_arr, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / ranging", &
                      "Could not deallocate memory (5).", 1)
                 STOP
              END IF
              ! Make plot using gnuplot:
              CALL system("cp " // TRIM(str) // &
                   "_sor_orbits.out sor_orbits.out")
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_car.gp")
              ELSE
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_kep.gp")
              END IF
              CALL system("cp sor_results.eps " // TRIM(str) // &
                   "_sor_" // TRIM(element_type_comp_prm) // &
                   "_results.eps")
              CALL system("rm -f sor_orbits.out sor_results.eps " // & 
                   TRIM(str) // "_sor_orbits.out " // TRIM(str) // &
                   "_sor_sample_standard_deviations.out")
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_sor_" // &
                      TRIM(element_type_comp_prm) // "_results.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_sor_" // &
                      TRIM(element_type_comp_prm) // &
                      "_results.eps* &")
              END IF
           END IF
           DO j=1,SIZE(orb_arr_cmp)
              CALL NULLIFY(orb_arr_cmp(j))
           END DO
           DEALLOCATE(orb_arr_cmp, stat=err)
           DEALLOCATE(pdf_arr_cmp, stat=err)
           DEALLOCATE(rchi2_arr_cmp, stat=err)
           DEALLOCATE(reg_apr_arr_cmp, stat=err)
           DEALLOCATE(jac_arr_cmp, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("oorb / ranging", &
                   "Could not deallocate memory (10).", 1)
              STOP
           END IF
           CALL NULLIFY(obss_sep(i))
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Object", &
                   TRIM(id), "successfully processed."
           END IF
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb >= 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
        CALL NULLIFY(obss_sep(i))
     END DO
     DEALLOCATE(obss_sep, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("oorb / ranging", &
             "Could not deallocate memory (15).", 1)
        STOP
     END IF


  CASE ("lsl")

     ! Orbital inversion using least squares with linearized
     ! covariances, that is, fixing the resulting shape of
     ! the orbital-element pdf to a Gaussian.

     CALL NULLIFY(epoch)
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          integrator=integrator, &
          integration_step=integration_step, &
          dyn_model_init=dyn_model_init, &
          integrator_init=integrator_init, &
          integration_step_init=integration_step_init, &
          ls_correction_factor=ls_correction_factor, &
          ls_element_mask=ls_element_mask, &
          ls_niter_major_max=ls_niter_major_max, &
          ls_niter_major_min=ls_niter_major_min, &
          ls_niter_minor=ls_niter_minor, &
          sor_norb=sor_norb, &
          sor_ntrial=sor_ntrial, &
          sor_norb_sw=sor_norb_sw, &
          sor_ntrial_sw=sor_ntrial_sw, &
          sor_niter=sor_niter, &
          sor_rho_init=sor_rho_init, &
          sor_genwin_multiplier=sor_genwin_multiplier, &
          accwin_multiplier=accwin_multiplier)
     IF (error) THEN
        CALL errorMessage("oorb / lsl", &
             "TRACE BACK (5)", 1)
        STOP
     END IF

     obss_sep => getSeparatedSets(obss_in)
     IF (error) THEN
        CALL errorMessage("oorb / lsl", &
             "TRACE BACK (20)", 1)
        STOP
     END IF
     CALL NULLIFY(obss_in)
     ! Print header before printing first orbit:
     first = .TRUE.
     DO i=1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (25)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")."
        END IF
        nobs = getNrOfObservations(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (nobs < 4) THEN
           CALL errorMessage("oorb / lsl", &
                "Too few observations:", 1)
           WRITE(stderr,*) "ID: ", TRIM(id), "  and number of observations: ", nobs
           CYCLE
        END IF
        norb = 0
        IF (ALLOCATED(storb_arr_in)) THEN
           DO j=1,SIZE(id_arr_in,dim=1)
              IF (id_arr_in(j) == id) THEN
                 IF (containsSampledPDF(storb_arr_in(j))) THEN
                    orb_arr => getSampleOrbits(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / lsl", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = SIZE(orb_arr)
                    EXIT
                 ELSE
                    ALLOCATE(orb_arr(1))
                    orb_arr(1) = getNominalOrbit(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / lsl", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = 1
                    EXIT                    
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / lsl", &
                   "Initial orbit not available.", 1)
              STOP
           END IF
        ELSE IF (ASSOCIATED(orb_arr_in)) THEN
           DO j=1,SIZE(id_arr_in,dim=1)
              IF (id_arr_in(j) == id) THEN
                 norb = norb + 1
                 IF (.NOT.ASSOCIATED(orb_arr)) THEN
                    orb_arr => reallocate(orb_arr,2*norb)
                 ELSE IF (SIZE(orb_arr) < norb) THEN
                    orb_arr => reallocate(orb_arr,2*norb)
                 END IF
                 orb_arr(norb) = copy(orb_arr_in(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / lsl", &
                   "Initial orbit not available.", 1)
              STOP
           ELSE
              orb_arr => reallocate(orb_arr,norb)
           END IF
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (70)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL NEW(storb, obss_sep(i))        
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (75)", 1)
           STOP
        END IF
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             accept_multiplier=accwin_multiplier, &
             ls_correction_factor=ls_correction_factor, &
             ls_element_mask=ls_element_mask, &
             ls_niter_major_max=ls_niter_major_max, &
             ls_niter_major_min=ls_niter_major_min, &
             ls_niter_minor=ls_niter_minor)
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (80)", 1)
           STOP
        END IF
        DO j=1,SIZE(orb_arr,dim=1)
           CALL setParameters(orb_arr(j), &
                dyn_model=dyn_model_init, &
                perturbers=perturbers, &
                integrator=integrator_init, &
                integration_step=integration_step_init)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (85)", 1)
              STOP
           END IF
           CALL propagate(orb_arr(j), t)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (90)", 1)
              STOP
           END IF
           CALL setParameters(orb_arr(j), &
                dyn_model=dyn_model, &
                perturbers=perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (95)", 1)
              STOP
           END IF
           CALL leastSquares(storb, orb_arr(j))
           IF (.NOT.error) THEN
              EXIT
           ELSE IF (j < norb) THEN
              error = .FALSE.
           END IF
        END DO
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "Least squares failed:", 1)
           IF (err_verb >= 1) THEN
              WRITE(stderr,"(3A,1X,I0)") "ID: ", TRIM(id), &
                   " and number of observations: ", nobs
           END IF
           error  = .FALSE.
           CALL NEW(out_file, "problematic_observation_sets" // "." // &
                TRIM(observation_format_out))
           CALL setPositionAppend(out_file)
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (120)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (125)", 1)
              STOP
           END IF
        ELSE
           obs_masks => getObservationMasks(storb)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (130)", 1)
              STOP
           END IF
           noutlier = 0
           DO j=1,SIZE(obs_masks,dim=1)
              IF (ALL(.NOT.obs_masks(j,:))) THEN
                 noutlier = noutlier + 1
              END IF
           END DO
           IF (noutlier > SIZE(obs_masks,dim=1)/5) THEN
              ! In case of too many outliers (>20% of obs), throw the
              ! observation set in a separate bin:
              CALL NEW(out_file, "observation_sets_with_many_outliers." // TRIM(observation_format_out))
              CALL setPositionAppend(out_file)
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (145)", 1)
                 STOP
              END IF
              WRITE(getUnit(out_file),"(1X)")
              WRITE(getUnit(out_file),"(A,I0,A)",advance="no") "# ", &
                   noutlier, " outliers: "
              DO j=1,SIZE(obs_masks,dim=1)
                 IF (ALL(.NOT.obs_masks(j,:))) THEN
                    WRITE(getUnit(out_file),"(A)",advance="no") "*"
                 ELSE
                    WRITE(getUnit(out_file),"(A)",advance="no") "-"
                 END IF
              END DO
              WRITE(getUnit(out_file),"(1X)")
              CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                   TRIM(observation_format_out))
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (150)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (155)", 1)
                 STOP
              END IF
              CYCLE
           END IF
           DEALLOCATE(obs_masks, stat=err)

           IF (pp_H_estimation) THEN
              CALL NEW(physparam, storb)
              IF (pp_G > 99.0_bp) THEN
                 CALL estimateHAndG(physparam, obss_sep(i))
              ELSE
                 CALL estimateHAndG(physparam, obss_sep(i), &
                      input_G=pp_G, input_delta_G=pp_G_unc)
              END IF
              HG_arr_in(i,1:2) = getH0(physparam)
              HG_arr_in(i,3:4) = getG(physparam)
           END IF


           CALL NEW(out_file, TRIM(out_fname) // ".ls")
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (160)", 1)
              STOP
           END IF
           CALL setPositionAppend(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (165)", 1)
              STOP
           END IF
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (170)", 1)
              STOP
           END IF

           ! WRITE OBSERVATIONS:
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (175)", 1)
              WRITE(getUnit(out_file),"(A)") &
                   "Could not write observations for object " // TRIM(id) 
              error = .FALSE.
           END IF
           WRITE(getUnit(out_file),"(A)") "#"

           ! WRITE LS RESULTS
           CALL writeNominalSolution(storb, obss_sep(i), &
                element_type_out_prm, getUnit(out_file))
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (180)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)

           ! ORBITAL ELEMENTS:
           CALL NULLIFY(orb)
           orb = getNominalOrbit(storb)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (185)", 1)
              STOP
           END IF
           cov = getCovarianceMatrix(storb, element_type_out_prm, "ecliptic")
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (190)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NEW(out_file, TRIM(id) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (200)", 1)
                 STOP
              END IF
              lu_orb_out = getUnit(out_file)
           END IF
           IF (orbit_format_out == "orb") THEN
              CALL writeOpenOrbOrbitFile(lu_orb_out, &
                   print_header=first, &
                   element_type_out=element_type_out_prm, &
                   id=id, &
                   orb=orb, &
                   cov=cov, &
                   H=HG_arr_in(i,1), &
                   G=HG_arr_in(i,3))
           ELSE IF (orbit_format_out == "des") THEN
              CALL errorMessage("oorb / lsl", &
                   "DES format not yet supported for LSL output.", 1)
              STOP                 
           END IF
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (195)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NULLIFY(out_file)
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(id) // ".orb")
              END IF
           END IF
           first = .FALSE.

           IF (plot_results) THEN
              ! Prepare data for plotting
              CALL toString(dt, str, error, frmt="(F10.2)")
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (200)", 1)
                 STOP
              END IF
              str = TRIM(id) // "_"// TRIM(str)
              elements = getElements(orb, element_type_comp_prm, "ecliptic")
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (205)", 1)
                 STOP
              END IF
              t = getTime(orb)
              CALL NULLIFY(orb)
              cov = getCovarianceMatrix(storb, element_type_comp_prm, "ecliptic")
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (207)", 1)
                 STOP
              END IF
              CALL NEW(tmp_file, TRIM(str)// &
                   "_ls_nominal_orbit_stdevs_corrs.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (215)", 1)
                 STOP
              END IF
              IF (element_type_comp_prm == "keplerian") THEN
                 elements(3:6) = elements(3:6)/rad_deg
              END IF
              WRITE(getUnit(tmp_file), "(6(E22.15,1X))", advance="no") &
                   elements(1:6)
              DO j=1,6
                 CALL toString(elements(j), element_str_arr(j), error, frmt="(F15.9)")
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (220)", 1)
                    STOP
                 END IF
                 stdev = SQRT(cov(j,j))
                 IF (element_type_comp_prm == "keplerian" .AND. j>=3) THEN
                    stdev = stdev/rad_deg
                 END IF
                 WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                      advance="no") stdev
                 CALL toString(stdev, stdev_str_arr(j), error, frmt="(F22.16)")
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (225)", 1)
                    STOP
                 END IF
              END DO
              stdev = SQRT(cov(1,1))
              DO j=2,6
                 WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                      advance="no") cov(1,j)/(stdev*SQRT(cov(j,j)))
                 CALL toString(cov(1,j)/(stdev*SQRT(cov(j,j))), &
                      corr_str_arr(j-1), error, frmt="(F15.9)")
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (230)", 1)
                    STOP
                 END IF
              END DO
              WRITE(getUnit(tmp_file), "(F22.15,1X,A)") &
                   getObservationalTimespan(obss_sep(i)), &
                   getCalendarDateString(t,"tdt")
              IF (error) THEN
                 CALL errorMessage("oorb / lsl", &
                      "TRACE BACK (235)", 1)
                 STOP
              END IF
              CALL NULLIFY(tmp_file)
              ! Make plot using gnuplot:
              CALL system("cp " // TRIM(str) // &
                   "_ls_nominal_orbit_stdevs_corrs.out " // &
                   "ls_nominal_orbit_stdevs_corrs.out")
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL system("cp " // TRIM(gnuplot_scripts_dir) // "/ls_plot.gp " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'x [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'y [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(2)) // &
                      "+" // TRIM(stdev_str_arr(2)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(1)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+ 3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // &
                      TRIM(element_str_arr(2)) // &
                      "+3*" // TRIM(stdev_str_arr(2)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(1)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:2 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.5,0.66 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'x [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'z [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(3)) // &
                      "+" // TRIM(stdev_str_arr(3)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(2)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(3)) // &
                      "+3*" // TRIM(stdev_str_arr(3)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(2)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:3 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.0,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'x [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'dx/dt [AU/d]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(4)) // &
                      "+" // TRIM(stdev_str_arr(4)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(3)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(4)) // &
                      "+3*" // TRIM(stdev_str_arr(4)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(3)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:4 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'x [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'dy/dt [AU/d]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(5)) // &
                      "+" // TRIM(stdev_str_arr(5)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(4)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(5)) // &
                      "+3*" // TRIM(stdev_str_arr(5)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(4)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:5 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.0,0.0 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'x [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'dz/dt [AU/d]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(6)) // &
                      "+" // TRIM(stdev_str_arr(6)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(5)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(6)) // &
                      "+3*" // TRIM(stdev_str_arr(6)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(5)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:6 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo unset parametric  >> " // TRIM(str) // ".gp")
                 CALL system("echo unset multiplot >> " // TRIM(str) // ".gp")
                 CALL system("echo reset >> " // TRIM(str) // ".gp")
                 CALL system("echo set terminal x11 >> " // TRIM(str) // ".gp")
                 CALL system("gnuplot " //  TRIM(str) // ".gp")
              ELSE IF (element_type_comp_prm == "keplerian") THEN
                 CALL system("cp " // TRIM(gnuplot_scripts_dir) // "/ls_plot.gp " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'a [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'e\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(2)) // &
                      "+" // TRIM(stdev_str_arr(2)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(1)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+ 3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // &
                      TRIM(element_str_arr(2)) // &
                      "+3*" // TRIM(stdev_str_arr(2)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(1)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:2 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.5,0.66 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'a [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'i [deg]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(3)) // &
                      "+" // TRIM(stdev_str_arr(3)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(2)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(3)) // &
                      "+3*" // TRIM(stdev_str_arr(3)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(2)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:3 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.0,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'a [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'\{/Symbol O\} [deg]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(4)) // &
                      "+" // TRIM(stdev_str_arr(4)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(3)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(4)) // &
                      "+3*" // TRIM(stdev_str_arr(4)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(3)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:4 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'a [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'\{/Symbol o\} [deg]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(5)) // &
                      "+" // TRIM(stdev_str_arr(5)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(4)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(5)) // &
                      "+3*" // TRIM(stdev_str_arr(5)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(4)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:5 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo set size 0.5,0.33 >> " // TRIM(str) // ".gp")
                 CALL system("echo set origin 0.0,0.0 >> " // TRIM(str) // ".gp")
                 CALL system("echo set xlabel \'a [AU]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set ylabel \'M [deg]\' >> " // TRIM(str) // ".gp")
                 CALL system("echo set parametric >> " // TRIM(str) // ".gp")
                 CALL system("echo plot " // TRIM(element_str_arr(1)) // "+" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(6)) // &
                      "+" // TRIM(stdev_str_arr(6)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(5)) // "*pi/2\), " // &
                      TRIM(element_str_arr(1)) // "+3*" // &
                      TRIM(stdev_str_arr(1)) // "*cos\(t\)," // TRIM(element_str_arr(6)) // &
                      "+3*" // TRIM(stdev_str_arr(6)) // "*sin\(t+" // &
                      TRIM(corr_str_arr(5)) // "*pi/2\), " // &
                      "\'ls_nominal_orbit_stdevs_corrs.out\' " // &
                      "using 1:6 with points pt 3 ps 3.0 >> " // &
                      TRIM(str) // ".gp")
                 CALL system("echo unset parametric  >> " // TRIM(str) // ".gp")
                 CALL system("echo unset multiplot >> " // TRIM(str) // ".gp")
                 CALL system("echo reset >> " // TRIM(str) // ".gp")
                 CALL system("echo set terminal x11 >> " // TRIM(str) // ".gp")
                 CALL system("gnuplot " //  TRIM(str) // ".gp")
              END IF
              CALL system("cp ls_results.eps " // TRIM(str) // &
                   "_ls_" // TRIM(element_type_comp_prm) // &
                   "_results.eps")
              CALL system("rm -f " // &
                   "ls_nominal_orbit_stdevs_corrs.out " // &
                   "ls_results.eps " // TRIM(str) // &
                   "_ls_nominal_orbit_stdevs_corrs.out " // &
                   TRIM(str) // ".gp")
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_ls_" // &
                      TRIM(element_type_comp_prm) // "_results.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_ls_" // &
                      TRIM(element_type_comp_prm) // &
                      "_results.eps* &")
              END IF
           END IF
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Object", &
                   TRIM(id), "successfully processed."
           END IF
        END IF
        DEALLOCATE(orb_arr, stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb / lsl", &
                "Could not deallocate memory (10)", 1)
           STOP
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb > 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
        CALL NULLIFY(obss_sep(i))
     END DO
     DEALLOCATE(obss_sep, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("oorb / lsl", &
             "Could not deallocate memory (15)", 1)
        STOP
     END IF

  CASE ("propagation")

     CALL readConfigurationFile(conf_file, &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          integrator=integrator, &
          integration_step=integration_step)

     CALL NULLIFY(epoch)
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / propagation", &
                "TRACE BACK (5)", 1)
           STOP
        END IF
     ELSE IF (get_cl_option("--epoch-mjd-utc=", .FALSE.)) THEN
        mjd_utc = get_cl_option("--epoch-mjd-utc=", 0.0_bp)
        CALL NEW(epoch, mjd_utc, "UTC")
        IF (error) THEN
           CALL errorMessage("oorb / propagation", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF
     IF (.NOT.exist(epoch) .AND. &
          .NOT.get_cl_option("--delta-epoch-mjd=", .FALSE.) .AND. &
          .NOT.exist(obss_in)) THEN
        CALL errorMessage("oorb / propagation", &
             "No epoch specified (5).", 1)
        STOP
     END IF

     first = .TRUE.
     IF (ALLOCATED(storb_arr_in)) THEN
        ! Input orbits contain uncertainty information.

        DO i=1,SIZE(storb_arr_in)
           ! Set integration parameters
           CALL setParameters(storb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (15)", 1)
              STOP
           END IF

           IF (.NOT.exist(epoch) .AND. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) THEN
              dt = get_cl_option("--delta-epoch-mjd=", 0.0_bp)
              epoch = getTime(storb_arr_in(i))
              mjd_tt = getMJD(epoch, "TT")
              CALL NULLIFY(epoch)
              CALL NEW(epoch, mjd_tt + dt, "TT")
           ELSE IF (.NOT.exist(epoch) .AND. & 
                .NOT.get_cl_option("--delta-epoch-mjd=", .FALSE.) .AND. &
                exist(obss_in)) THEN
              dt = get_cl_option("--delta-epoch-mjd=", 0.0_bp)
              epoch = getTime(storb_arr_in(i))
              mjd_tt = getMJD(epoch, "TT")
              CALL NULLIFY(epoch)
              CALL NEW(epoch, mjd_tt + dt, "TT")
           END IF

           IF (info_verb >= 2) THEN
              epoch0 = getTime(storb_arr_in(i))
           END IF

           ! Propagate orbital-element pdf from one epoch (=input) to another:
           CALL propagate(storb_arr_in(i), epoch, encounters=encounters)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (20)", 1)
              STOP
           END IF

           IF (info_verb >= 2) THEN
              CALL getCalendarDate(epoch0, "TT", year0, month0, day0)
              CALL getCalendarDate(epoch, "TT", year1, month1, day1)
              CALL NULLIFY(epoch0)
              WRITE(stderr,'(A)') ""
              WRITE(stderr,'(A,I0,A,I0,A,F0.5,A,I0,A,I0,F0.5,A)') &
                   "Planetary encounters by object " // TRIM(id_arr_in(i)) // &
                   " between ", year0, "-", month0, "-", day0, " and ", &
                   year1, "-", month1, "-", day1, ":"
              DO j=1,SIZE(encounters,dim=1)
                 WRITE(stderr,'(A,I0,A)') "Orbit #", j, ":"
                 DO k=1,SIZE(encounters,dim=2)
                    IF (encounters(j,k,2) < 1.1_bp) THEN
                       WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F12.6,1X,A,1X,F0.6,1X,A)') &
                            "!! I-M-P-A-C-T !! with", planetary_locations(k), &
                            "at MJD", encounters(j,k,1), "TT given a stepsize of", &
                            encounters(j,k,4), "days."
                    ELSE
                       WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F12.6,1X,A,1X,F11.8,1X,A,1X,F0.6,1X,A)') &
                            "Closest encounter with", planetary_locations(k), &
                            "at MJD", encounters(j,k,1), &
                            "TT at a distance of", encounters(j,k,3), &
                            "AU given a stepsize of", encounters(j,k,4), "days."
                    END IF
                 END DO
              END DO
           END IF
           DEALLOCATE(encounters, stat=err)

           IF (containsSampledPDF(storb_arr_in(1))) THEN
              ! Sampled orbital-element pdf:
              orb_arr_cmp => getSampleOrbits(storb_arr_in(i))
              IF (error) THEN
                 CALL errorMessage("oorb / propagation ", &
                      "TRACE BACK (25)", 1)
                 STOP
              END IF
              pdf_arr_cmp => getPDFValues(storb_arr_in(i))
              IF (error) THEN
                 CALL errorMessage("oorb / propagation ", &
                      "TRACE BACK (30)", 1)
                 STOP
              END IF
              rchi2_arr_cmp => getReducedChi2Distribution(storb_arr_in(i))
              IF (error) THEN
                 CALL errorMessage("oorb / propagation ", &
                      "TRACE BACK (35)", 1)
                 STOP
              END IF
              CALL getResults(storb_arr_in(i), &
                   reg_apr_arr=reg_apr_arr_cmp, &
                   jac_arr=jac_arr_cmp)
              IF (error) THEN
                 CALL errorMessage("oorb / propagation ", &
                      "TRACE BACK (40)", 1)
                 STOP
              END IF
              IF (separately) THEN
                 CALL NEW(out_file, TRIM(id_arr_in(i)) // ".orb")
                 CALL OPEN(out_file)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (45)", 1)
                    STOP
                 END IF
                 lu_orb_out = getUnit(out_file)
              END IF
              DO j=1,SIZE(orb_arr_cmp,dim=1)
                 IF (orbit_format_out == "orb") THEN
                    CALL writeOpenOrbOrbitFile(lu_orb_out, &
                         print_header=first.OR.separately, &
                         element_type_out=element_type_out_prm, &
                         id=id_arr_in(i), &
                         orb=orb_arr_cmp(j), &
                         element_type_pdf=element_type_comp_prm, &
                         pdf=pdf_arr_cmp(j), &
                         rchi2=rchi2_arr_cmp(j), &
                         reg_apr=reg_apr_arr_cmp(j), &
                         jac_sph_inv=jac_arr_cmp(j,1), &
                         jac_car_kep=jac_arr_cmp(j,2), &
                         jac_equ_kep=jac_arr_cmp(j,3), &
                         H=HG_arr_in(i,1), &
                         G=HG_arr_in(i,3))
                 ELSE IF (orbit_format_out == "des") THEN
                    CALL errorMessage("oorb / propagation", &
                         "DES format not yet supported for propagation of sampled pdfs.", 1)
                    STOP                 
                 END IF
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation ", &
                         "TRACE BACK (50)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(orb_arr_cmp(j))
                 first = .FALSE.
              END DO
              IF (separately) THEN
                 CALL NULLIFY(out_file)
                 IF (compress) THEN
                    CALL system("gzip -f " // TRIM(id_arr_in(i)) // ".orb")
                 END IF
              END IF
              DEALLOCATE(orb_arr_cmp, pdf_arr_cmp, rchi2_arr_cmp, &
                   reg_apr_arr_cmp, jac_arr_cmp)
           ELSE
              ! LSL orbit:
              orb = getNominalOrbit(storb_arr_in(i))
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (55)", 1)
                 STOP
              END IF
              cov = getCovarianceMatrix(storb_arr_in(i), element_type_out_prm, "ecliptic")
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (60)", 1)
                 STOP
              END IF
              IF (separately) THEN
                 CALL NEW(out_file, TRIM(id_arr_in(i)) // ".orb")
                 CALL OPEN(out_file)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (65)", 1)
                    STOP
                 END IF
                 lu_orb_out = getUnit(out_file)
              END IF
              IF (orbit_format_out == "orb") THEN
                 CALL writeOpenOrbOrbitFile(lu_orb_out, &
                      print_header=first.OR.separately, &
                      element_type_out=element_type_out_prm, &
                      id=id_arr_in(i), &
                      orb=orb, &
                      cov=cov, &
                      H=HG_arr_in(i,1), &
                      G=HG_arr_in(i,3))
              ELSE IF (orbit_format_out == "des") THEN
                 CALL errorMessage("oorb / propagation", &
                      "DES format not yet supported for propagation of covariance matrices.", 1)
                 STOP                 
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (70)", 1)
                 STOP
              END IF
              IF (separately) THEN
                 CALL NULLIFY(out_file)
                 IF (compress) THEN
                    CALL system("gzip -f " // TRIM(id_arr_in(i)) // ".orb")
                 END IF
              END IF
              CALL NULLIFY(orb)
              first = .FALSE.
           END IF
        END DO

     ELSE

        IF (exist(obss_in)) THEN
           obss_sep => getSeparatedSets(obss_in)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (10)", 1)
              STOP
           END IF
           CALL NULLIFY(obss_in)
        END IF

        dt = get_cl_option("--delta-epoch-mjd=", 0.0_bp)
        DO i=1,SIZE(orb_arr_in)
           CALL setParameters(orb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (75)", 1)
              STOP
           END IF
           IF (.NOT.exist(epoch) .AND. &
                .NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.)) .AND. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) THEN
              epoch = getTime(orb_arr_in(i))
              mjd_tt = getMJD(epoch, "TT")
              CALL NULLIFY(epoch)
              CALL NEW(epoch, mjd_tt + dt, "TT")
           ELSE IF (.NOT.exist(epoch) .AND. &
                .NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.) .OR. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) .AND. &
                ASSOCIATED(obss_sep)) THEN
              DO j=1,SIZE(obss_sep)
                 IF (id_arr_in(i) == getID(obss_sep(j))) THEN
                    EXIT
                 END IF
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (25)", 1)
                    STOP
                 END IF
              END DO
              IF (get_cl_option("--last-observation-date", .FALSE.)) THEN
                 obs = getObservation(obss_sep(j),getNrOfObservations(obss_sep(j)))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
                 epoch = getTime(obs)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (40)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(obs)                 
              ELSE ! observational mid-date
                 dt = getObservationalTimespan(obss_sep(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (30)", 1)
                    STOP
                 END IF
                 obs = getObservation(obss_sep(j),1)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
                 epoch = getTime(obs)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (40)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(obs)
                 mjd = getMJD(epoch, "tt")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (45)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(epoch)
                 mjd = mjd + dt/2.0_bp
                 CALL NEW(epoch, mjd, "tt")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (50)", 1)
                    STOP
                 END IF
              END IF
           ELSE IF (.NOT.exist(epoch)) THEN
              CALL errorMessage("oorb / propagation", &
                   "No epoch specified (10).", 1)
              STOP              
           END IF
           IF (info_verb >= 2) THEN
              epoch0 = getTime(orb_arr_in(i))
           END IF
           CALL propagate(orb_arr_in(i), epoch, encounters=encounters)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (80)", 1)
              STOP
           END IF
           IF (info_verb >= 2) THEN
              CALL getCalendarDate(epoch0, "TT", year0, month0, day0)
              CALL getCalendarDate(epoch, "TT", year1, month1, day1)
              CALL NULLIFY(epoch0)
              WRITE(stderr,'(A)') ""
              WRITE(stderr,'(A,I0,A,I0,A,I0,A,F0.5,A,I0,A,I0,A,F0.5,A)') &
                   "Planetary encounters by object " // TRIM(id_arr_in(i)) // &
                   " (orbit #", i, ") between ", year0, "-", month0, &
                   "-", day0, " and ", year1, "-", month1, "-", day1, ":"
              DO j=1,SIZE(encounters,dim=2)
                 IF (encounters(1,j,2) < 1.1_bp) THEN
                    WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F12.6,1X,A,1X,F0.6,1X,A)') &
                         "!! I-M-P-A-C-T !! with", planetary_locations(j), &
                         "at MJD", encounters(1,j,1), "TT given a stepsize of", &
                         encounters(1,j,4), "days."
                 ELSE
                    WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F12.6,1X,A,1X,F11.8,1X,A,1X,F0.6,1X,A)') &
                         "Closest encounter with", planetary_locations(j), &
                         "at MJD", encounters(1,j,1), "TT at a distance of", &
                         encounters(1,j,3), "AU given a stepsize of", encounters(1,j,4), "days."
                 END IF
              END DO
           END IF
           DEALLOCATE(encounters, stat=err)
           IF (.NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.)) .AND. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) THEN
              CALL NULLIFY(epoch)
           END IF
        END DO
        DO i=1,SIZE(orb_arr_in)
           IF (orbit_format_out == "des") THEN
              CALL writeDESOrbitFile(lu_orb_out, i==1, "cometary", &
                   id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), 1, 6, &
                   -1.0_bp, "OPENORB")
           ELSE IF (orbit_format_out == "orb") THEN
              CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                   element_type_out=element_type_out_prm, &
                   id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
                   H=HG_arr_in(i,1), G=HG_arr_in(i,3))
           END IF
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (85) " // TRIM(id_arr_in(i)), 1)
              STOP              
           END IF
        END DO

     END IF

  CASE ("ephemeris")

     CALL readConfigurationFile(conf_file, &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          integrator=integrator, &
          integration_step=integration_step)

     ! Input observatory code
     obsy_code = get_cl_option("--code=", obsy_code)

     ! Input evolutionary timespan [days]
     timespan = get_cl_option("--timespan=", 0.0_bp)

     ! Input time step [days]
     step = get_cl_option("--step=", 1.0_bp)
     IF (step == 0.0_bp) THEN
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

     CALL NEW(obsies)
     IF (error) THEN
        CALL errorMessage('oorb / ephemeris', &
             'TRACE BACK (5)',1)
        STOP
     END IF

     IF (ALLOCATED(storb_arr_in)) THEN
        ! Input orbits contain uncertainty information.

        DO i=1,SIZE(storb_arr_in)

           ! Use equatorial coordinates:
           CALL toCartesian(storb_arr_in(i), "equatorial")

           IF (exist(obss_in)) THEN
              observers => getObservatoryCCoords(obss_in)
              DO j=1,SIZE(observers)
                 CALL rotateToEquatorial(observers(j))
              END DO
              obsy_code_arr => getObservatoryCodes(obss_in)
           ELSE
              t = getTime(storb_arr_in(i))
              mjd_tt = getMJD(t, "TT")
              CALL NULLIFY(t)
              ALLOCATE(observers(nstep), obsy_code_arr(nstep))
              DO j=1,nstep
                 CALL NEW(t, mjd_tt+(j-1)*step, "TT")
                 IF (error) THEN
                    CALL errorMessage("oorb / ephemeris", &
                         "TRACE BACK (10)", 1)
                    STOP
                 END IF
                 obsy_code_arr(j) = obsy_code
                 ! Compute heliocentric observatory coordinates
                 observers(j) = getObservatoryCCoord(obsies, obsy_code_arr(j), t)
                 IF (error) THEN
                    CALL errorMessage('oorb / ephemeris', &
                         'TRACE BACK (15)',1)
                    STOP
                 END IF
                 CALL rotateToEquatorial(observers(j))
                 CALL NULLIFY(t)
              END DO
           END IF

           ! Set integration parameters
           CALL setParameters(storb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / ephemeris", &
                   "TRACE BACK (20)", 1)
              STOP
           END IF

           ! Compute topocentric ephemerides
           CALL getEphemerides(storb_arr_in(i), observers, ephemerides_arr, &
                cov_arr=cov_arr, pdfs_arr=pdfs_arr)
           IF (error) THEN
              CALL errorMessage('oorb / ephemeris', &
                   'TRACE BACK (25)',1)
              STOP
           END IF

           ! Prepare output file
           IF (separately) THEN
              CALL NEW(tmp_file, TRIM(id_arr_in(i)) // ".eph")
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (25)',1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (30)',1)
                 STOP
              END IF
              lu = getUnit(tmp_file)
           ELSE
              lu = stdout
           END IF
           IF (separately .OR. i == 1) THEN
              IF (containsSampledPDF(storb_arr_in(1))) THEN
                 WRITE(lu,"(A,A11,5X,9(1X,A18))") "#", "Designation  ", &
                      "Observatory_code ", "MJD_UTC  ", "Delta  ", "RA  ", &
                      "Dec  ", "dDelta/dt  ", "dRA/dt  ", "dDec/dt  ", "PDF_value  "
              ELSE
                 WRITE(lu,"(A,A11,5X,29(1X,A18))") "#", "Designation  ", &
                      "Observatory_code ", "MJD_UTC  ", "Delta  ", "RA  ", "Dec  ", &
                      "dDelta/dt  ", "dRA/dt  ", "dDec/dt  ", "Delta_unc  ", &
                      "RA_unc  ", "Dec_unc  ", "dDelta/dt_unc  ", "dRA/dt_unc  ", &
                      "dDec/dt_unc  ", "cor(Delta,RA)  ", "cor(Delta,Dec)  ", &
                      "cor(Delta,dDelta)  ", "cor(Delta,dRA)  ", &
                      "cor(Delta,dDec)  ", "cor(RA,Dec)  ", &
                      "cor(RA,dDelta)  ", "cor(RA,dRA)  ", &
                      "cor(RA,dDec)  ", "cor(Dec,dDelta)  ", &
                      "cor(Dec,dRA)  ", "cor(Dec,dDec)  ", &
                      "cor(dDelta,dRA)  ", "cor(dDelta,dDec)  ", &
                      "cor(dRA,dDec)  "
              END IF
           END IF

           ! Loop over time steps:
           DO j=1,SIZE(observers)
              t = getTime(observers(j))
              mjd_utc = getMJD(t, "UTC")
              CALL NULLIFY(t)
              IF (containsSampledPDF(storb_arr_in(i))) THEN
                 ! Input orbits correspond to one or more sampled pdfs.
                 ! Loop over sampled orbits:
                 DO k=1,SIZE(ephemerides_arr,dim=1)
                    ! Make sure that the ephemeris is equatorial:
                    CALL rotateToEquatorial(ephemerides_arr(k,j))
                    coordinates = getCoordinates(ephemerides_arr(k,j))
                    WRITE(lu,"(2(A18,1X),7(F18.10,1X),E18.10)") &
                         TRIM(id_arr_in(i)), TRIM(obsy_code_arr(j)), &
                         mjd_utc, coordinates(1), &
                         coordinates(2:3)/rad_deg, coordinates(4), &
                         coordinates(5:6)/rad_deg, pdfs_arr(k,j)
                    CALL NULLIFY(ephemerides_arr(k,j))
                 END DO
              ELSE
                 ! Input orbits correspond to one or more single-point estimates of the pdf.
                 ! Make sure that the ephemeris is equatorial:
                 CALL rotateToEquatorial(ephemerides_arr(1,j))
                 coordinates = getCoordinates(ephemerides_arr(1,j))
                 DO k=1,6
                    stdev_arr(k) = SQRT(cov_arr(k,k,j)) 
                 END DO
                 DO k=1,6
                    DO l=1,6
                       corr(k,l) = cov_arr(k,l,j) / &
                            (stdev_arr(k)*stdev_arr(l))
                    END DO
                 END DO
                 WRITE(lu,"(2(A18,1X),28(F18.10,1X))") TRIM(id_arr_in(i)), &
                      TRIM(obsy_code_arr(j)), mjd_utc, coordinates(1), &
                      coordinates(2:3)/rad_deg, coordinates(4), &
                      coordinates(5:6)/rad_deg, stdev_arr(1), &
                      stdev_arr(2:3)/rad_deg, stdev_arr(4), &
                      stdev_arr(5:6)/rad_deg, corr(1,2:6), &
                      corr(2,3:6), corr(3,4:6), corr(4,5:6), corr(5,6)
                 CALL NULLIFY(ephemerides_arr(1,j))
              END IF
           END DO
           DEALLOCATE(ephemerides_arr)
           IF (ASSOCIATED(pdfs_arr)) THEN
              DEALLOCATE(pdfs_arr)
           END IF
           IF (ASSOCIATED(cov_arr)) THEN
              DEALLOCATE(cov_arr)
           END IF
           CALL NULLIFY(storb_arr_in(i))
           IF (separately) THEN
              CALL NULLIFY(tmp_file)
           END IF
        END DO
        DEALLOCATE(storb_arr_in, observers)

     ELSE


        DO i=1,norb

           IF (exist(obss_in)) THEN
              observers => getObservatoryCCoords(obss_in)
              DO j=1,SIZE(observers)
                 CALL rotateToEquatorial(observers(j))
              END DO
              obsy_code_arr => getObservatoryCodes(obss_in)
           ELSE
              t = getTime(orb_arr_in(i))
              mjd_tt = getMJD(t, "TT")
              ALLOCATE(observers(nstep), obsy_code_arr(nstep))
              DO j=1,nstep
                 obsy_code_arr(j) = obsy_code
                 ! Compute heliocentric observatory coordinates
                 observers(j) = getObservatoryCCoord(obsies, obsy_code_arr(j), t)
                 IF (error) THEN
                    CALL errorMessage('oorb / ephemeris', &
                         'TRACE BACK (10)',1)
                    STOP
                 END IF
                 CALL rotateToEquatorial(observers(j))
                 CALL NULLIFY(t)
                 CALL NEW(t, mjd_tt+j*step, "TT")
                 IF (error) THEN
                    CALL errorMessage("oorb / ephemeris", &
                         "TRACE BACK (15)", 1)
                    STOP
                 END IF
              END DO
           END IF

           ! Set integration parameters
           CALL setParameters(orb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / ephemeris", &
                   "TRACE BACK (15)", 1)
              STOP
           END IF

           ! Compute topocentric ephemerides
           CALL getEphemerides(orb_arr_in(i), observers, ephemerides, &
                this_lt_corr_arr=orb_lt_corr_arr)
           IF (error) THEN
              CALL errorMessage('oorb / ephemeris', &
                   'TRACE BACK (20)',1)
              STOP
           END IF

           IF (separately) THEN
              CALL NEW(tmp_file, TRIM(id_arr_in(i)) // ".eph")
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (25)',1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (30)',1)
                 STOP
              END IF
              lu = getUnit(tmp_file)
           ELSE
              lu = stdout
           END IF

           IF (separately .OR. i == 1) THEN
              WRITE(lu,'(A,A11,1X,A,1X,35(A18,1X))') "#", &
                   "Designation", "Code", "MJD_UTC/UT1", "Delta", &
                   "RA", "Dec", "dDelta/dt", "dRA/dt", "dDec/dt", &
                   "VMag", "Alt", "Phase", "LunarElon", "LunarAlt", &
                   "LunarPhase", "SolarElon", "SolarAlt", "r", &
                   "HLon", "HLat", "TLon", "TLat", "TOCLon", &
                   "TOCLat", "HOCLon", "HOCLat", "TOppLon", &
                   "TOppLat", "HEclObj_X", "HEclObj_Y", "HEclObj_Z", &
                   "HEclObj_dX/dt", "HEclObj_dY/dt", & 
                   "HEclObj_dZ/dt", "HEclObsy_X", "HEclObsy_Y", &
                   "HEclObsy_Z" 
           END IF

           DO j=1,SIZE(observers)

              t = getTime(observers(j))
              mjd_tt = getMJD(t, "TT")
              err_verb_ = err_verb
              err_verb = 0
              mjd_utc = getMJD(t, "UTC")
              err_verb = err_verb_
              IF (error) THEN
                 error = .FALSE.
                 mjd_utc = getMJD(t, "UT1")
                 IF (error) THEN
                    CALL errorMessage('oorb / ephemeris', &
                         'TRACE BACK (30)',1)
                    STOP
                 END IF
              END IF

              ! Extract topocentic equatorial coordinates
              CALL rotateToEquatorial(ephemerides(j))        
              comp_coord = getCoordinates(ephemerides(j))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (35)',1)
                 STOP
              END IF
              Delta = comp_coord(1)
              ra = comp_coord(2)
              dec = comp_coord(3)

              ! Extract instantaneous topocentric equatorial sky
              ! motions (for coordinate motions, divide sky dra by
              ! cos(dec), ie remove cos(comp_coord(3)))
              dDelta = comp_coord(4)
              dra = comp_coord(5)*COS(comp_coord(3))
              ddec = comp_coord(6)

              ! Extract topocentric ecliptic lon and lat
              CALL rotateToEcliptic(ephemerides(j))        
              comp_coord = getCoordinates(ephemerides(j))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (40)',1)
                 STOP
              END IF
              tlon = comp_coord(2)
              tlat = comp_coord(3)

              ! Compute topocentric opposition coordinates (actually,
              ! get heliocentric observatory coordinates)
              scoord = getSCoord(observers(j))
              CALL rotateToEcliptic(scoord)
              pos_opp = getPosition(scoord)
              opplon = pos_opp(2)
              opplat = pos_opp(3)
              CALL NULLIFY(scoord)

              ! Compute opposition-centered topocentric ecliptic coordinates
              toclon = tlon - opplon
              toclat = tlat - opplat
              IF (toclon > pi) THEN
                 toclon = toclon - two_pi
              ELSE IF (toclon < -pi) THEN
                 toclon = toclon + two_pi
              END IF

              ! Extract heliocentric ecliptic lon and lat
              ccoord = getCCoord(orb_lt_corr_arr(j), "ecliptic")
              scoord = getSCoord(ccoord)
              comp_coord = getCoordinates(scoord)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (45)',1)
                 STOP
              END IF
              hlon = comp_coord(2)
              hlat = comp_coord(3)
              CALL NULLIFY(ccoord)
              CALL NULLIFY(scoord)

              ! Compute opposition-centered topocentric ecliptic coordinates
              hoclon = hlon - opplon
              hoclat = hlat - opplat
              IF (hoclon > pi) THEN
                 hoclon = hoclon - two_pi
              ELSE IF (hoclon < -pi) THEN
                 hoclon = hoclon + two_pi
              END IF

              ! Compute phase angle
              CALL NEW(ccoord, ephemerides(j))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (50)',1)
                 STOP
              END IF
              CALL rotateToEquatorial(ccoord)
              obsy_obj = getPosition(ccoord)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (55)',1)
                 STOP
              END IF
              ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
              CALL toCartesian(orb_lt_corr_arr(j), frame='equatorial')
              pos = getPosition(orb_lt_corr_arr(j))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (60)',1)
                 STOP
              END IF
              heliocentric_r2 = DOT_PRODUCT(pos,pos)
              obsy_pos = getPosition(observers(j))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (65)',1)
                 STOP
              END IF
              observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)
              cos_obj_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
                   observer_r2) / (SQRT(heliocentric_r2) * &
                   SQRT(ephemeris_r2))
              obj_phase = ACOS(cos_obj_phase)

              ! Compute apparent brightness:
              ! Input slope parameter
              G_value = get_cl_option("--G=", HG_arr_in(i,3))
              obj_vmag = getApparentMagnitude(H=HG_arr_in(i,1), &
                   G=G_value, r=SQRT(heliocentric_r2), &
                   Delta=Delta, phase_angle=obj_phase)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (70)',1)
                 STOP
              END IF

              ! Compute (approximate) altitude of the target
              obsy_ccoord = getGeocentricObservatoryCCoord(obsies, obsy_code_arr(j), t)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (75)',1)
                 STOP
              END IF
              CALL rotateToEquatorial(obsy_ccoord)
              geoc_obsy = getPosition(obsy_ccoord)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (80)',1)
                 STOP
              END IF
              vec3 = cross_product(geoc_obsy,obsy_obj)
              obj_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_obj))

              ! Compute (approximate) altitude of the Sun
              ! Position of the geocenter as seen from the Sun:
              planeph => JPL_ephemeris(mjd_tt, 3, 11, error)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (85)',1)
                 STOP
              END IF
              ! Position of the Sun as seen from the observatory:
              obsy_sun = -(planeph(1,1:3) + geoc_obsy)
              DEALLOCATE(planeph)
              vec3 = cross_product(geoc_obsy,obsy_sun)
              solar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_sun))
              ! Compute the solar elongation:
              vec3 = cross_product(obsy_obj,obsy_sun)
              solar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_sun))

              ! Compute phase of the Moon:
              ! Position of the geocenter as seen from the Moon:
              planeph => JPL_ephemeris(mjd_tt, 3, 10, error)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (90)',1)
                 STOP
              END IF
              ! Position of the Moon as seen from the observatory:
              obsy_moon = -(planeph(1,1:3) + geoc_obsy)
              DEALLOCATE(planeph)
              ! Position of the Sun as seen from the Moon:
              planeph => JPL_ephemeris(mjd_tt, 10, 11, error)
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (95)',1)
                 STOP
              END IF
              sun_moon = planeph(1,1:3)
              DEALLOCATE(planeph)
              obsy_moon_r2 = DOT_PRODUCT(obsy_moon,obsy_moon)
              sun_moon_r2 = DOT_PRODUCT(sun_moon,sun_moon)
              cos_obj_phase = 0.5_bp * (sun_moon_r2 + obsy_moon_r2 - &
                   observer_r2) / (SQRT(sun_moon_r2) * &
                   SQRT(obsy_moon_r2))
              lunar_phase = (pi-ACOS(cos_obj_phase))/pi
              ! Compute (approximate) distance between the target and the Moon:
              vec3 = cross_product(obsy_obj,obsy_moon)
              lunar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_moon))
              ! Compute (approximate) altitude of the Moon:
              vec3 = cross_product(geoc_obsy,obsy_moon)
              lunar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_moon))

              ! Extract heliocentric distance
              hdist = SQRT(heliocentric_r2)

              ! Extract heliocentric ecliptic cartesian coordinates for the object
              h_ecl_car_coord_obj = getElements(orb_lt_corr_arr(j), "cartesian", "ecliptic")
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (90)',1)
                 STOP
              END IF

              ! Extract heliocentric ecliptic cartesian coordinates for the observer
              CALL rotateToEcliptic(observers(j))
              h_ecl_car_coord_obsy = getCoordinates(observers(j))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (95)',1)
                 STOP
              END IF

              CALL NULLIFY(observers(j))
              CALL NULLIFY(ephemerides(j))
              CALL NULLIFY(orb_lt_corr_arr(j))
              CALL NULLIFY(ccoord)
              CALL NULLIFY(obsy_ccoord)

              ra = ra/rad_deg
              dec = dec/rad_deg
              dDelta = dDelta
              dra = dra/rad_deg
              ddec = ddec/rad_deg
              obj_alt = obj_alt/rad_deg
              obj_phase = obj_phase/rad_deg
              lunar_elongation = lunar_elongation/rad_deg
              lunar_alt = lunar_alt/rad_deg
              solar_elongation = solar_elongation/rad_deg
              solar_alt = solar_alt/rad_deg
              tlon = tlon/rad_deg
              tlat = tlat/rad_deg
              toclon = toclon/rad_deg
              toclat = toclat/rad_deg
              hlon = hlon/rad_deg
              hlat = hlat/rad_deg
              hoclon = hoclon/rad_deg
              hoclat = hoclat/rad_deg
              opplon = opplon/rad_deg
              opplat = opplat/rad_deg

              WRITE(lu,'(2(A,1X),35(F18.10,1X))') &
                   id_arr_in(i), TRIM(obsy_code_arr(j)), mjd_utc, Delta, &
                   ra, dec, dDelta, dra, ddec, obj_vmag, obj_alt, &
                   obj_phase, lunar_elongation, lunar_alt, &
                   lunar_phase, solar_elongation, solar_alt, hdist, &
                   hlon, hlat, tlon, tlat, toclon, toclat, hoclon, &
                   hoclat, opplon, opplat, h_ecl_car_coord_obj, &
                   h_ecl_car_coord_obsy(1:3)

           END DO

           IF (separately) THEN
              CALL NULLIFY(tmp_file)
           END IF
           DEALLOCATE(observers, ephemerides, orb_lt_corr_arr)

        END DO

     END IF


  CASE ("classification")

     ! Compute probability for an object with the input
     ! orbital-element pdf to belong to a group of asteroids. Caveat:
     ! an, ap, M are not currently taken into account, that is,
     ! e.g. Jupiter Trojans not correctly accounted for.

     IF (.NOT.ALLOCATED(storb_arr_in)) THEN
        CALL errorMessage("oorb / classification", &
             "Input orbits do not contain uncertainty " // &
             "information required for this task.", 1)        
        STOP
     END IF

     ! Write header
     WRITE(stdout,"(A2,3X,A12,9X,A6,7X,A12)") "# ", &
          "Designation", "Group", "Probability"

     ! Compute probabilities for each object
     DO i=1,SIZE(storb_arr_in)
        IF (.NOT.containsSampledPDF(storb_arr_in(i))) THEN
           CALL errorMessage("oorb / classification", &
                "Input orbits do not contain sampled " // &
                "uncertainty information required for this task.", 1)        
           STOP
        END IF
        ! Propagate input orbital-element pdf to Keplerian-element pdf:
        CALL toKeplerian(storb_arr_in(i))
        ! Compute weights for each class that has been defined:
        CALL getGroupWeights(storb_arr_in(i), weight_arr, group_name_arr)
        ! Write output
        DO j=1,SIZE(weight_arr)
           WRITE(stdout,"(3X,A16,1X,A18,1X,F12.7)") TRIM(id_arr_in(i)), &
                group_name_arr(j), weight_arr(j)
        END DO
        WRITE(stdout,*)
        DEALLOCATE(group_name_arr, weight_arr)
        CALL NULLIFY(storb_arr_in(i))
     END DO
     DEALLOCATE(storb_arr_in, id_arr_in)

  CASE ("classification_apriori")

     ! Compute probability for an object with the input
     ! orbital-element pdf to belong to a group of asteroids. Caveat:
     ! an, ap, M are not currently taken into account, that is,
     ! e.g. Jupiter Trojans not correctly accounted for.

     IF (.NOT.ALLOCATED(storb_arr_in)) THEN
        CALL errorMessage("oorb / classification", &
             "Input orbits do not contain uncertainty " // &
             "information required for this task.", 1)        
        STOP
     END IF

     ! Write header
     WRITE(stdout,"(A2,3X,A12,9X,A6,7X,A12)") "# ", &
          "Designation", "Group", "Probability"

     ! Compute probabilities for each object
     DO i=1,SIZE(storb_arr_in)
        IF (.NOT.containsSampledPDF(storb_arr_in(i))) THEN
           CALL errorMessage("oorb / classification", &
                "Input orbits do not contain sampled " // &
                "uncertainty information required for this task.", 1)        
           STOP
        END IF
        ! Propagate input orbital-element pdf to Keplerian-element pdf:
        CALL toKeplerian(storb_arr_in(i))
        ! Compute weights for each class that has been defined:
        CALL getGroupWeights(storb_arr_in(i), weight_arr, group_name_arr)
        weight_arr(5) = weight_arr(5)*35.0_bp
        weight_arr(6:SIZE(weight_arr)-1) = weight_arr(6:SIZE(weight_arr)-1)*400.0_bp
        weight_arr = weight_arr/SUM(weight_arr)
        ! Write output
        DO j=1,SIZE(weight_arr)
           WRITE(stdout,"(3X,A16,1X,A18,1X,F12.7)") TRIM(id_arr_in(i)), &
                group_name_arr(j), weight_arr(j)
        END DO
        WRITE(stdout,*)
        DEALLOCATE(group_name_arr, weight_arr)
        CALL NULLIFY(storb_arr_in(i))
     END DO
     DEALLOCATE(storb_arr_in, id_arr_in)


  CASE("apoapsis_distance")

     ! Compute apoapsis distance(s) for input orbit(s)
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (containsSampledPDF(storb_arr_in(i))) THEN
              CALL getApoapsisDistance(storb_arr_in(i), apoapsis_distance_pdf)
              DO j=1,SIZE(apoapsis_distance_pdf,dim=1)
                 WRITE(stdout,*) apoapsis_distance_pdf(j,:) 
              END DO
              DEALLOCATE(apoapsis_distance_pdf)
           ELSE
              CALL errorMessage("oorb / apoapsis_distance", &
                   "Input orbits do not contain sampled " // &
                   "uncertainty information currently required for this task.", 1)        
              STOP
              !              orb = getNominalOrbit(storb_arr_in(i))
              !              cov = getCovarianceMatrix(storb_arr_in(i), element_type_out_prm)
              !              CALL writeOpenOrbOrbitFile(lu_orb_out, i==1, &
              !                   element_type_out_prm, id_arr_in(i), &
              !                   orb=orb, cov=cov, H=HG_arr_in(i,1), &
              !                   G=HG_arr_in(i,2))
           END IF
        END DO
     ELSE
        CALL errorMessage("oorb / apoapsis_distance", &
             "Input orbits do not contain sampled " // &
             "uncertainty information currently required for this task.", 1)        
        STOP
     END IF


  CASE("periapsis_distance")

     ! Compute periapsis distance(s) for input orbit(s)
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (containsSampledPDF(storb_arr_in(i))) THEN
              CALL getPeriapsisDistance(storb_arr_in(i), periapsis_distance_pdf)
              DO j=1,SIZE(periapsis_distance_pdf,dim=1)
                 WRITE(stdout,*) periapsis_distance_pdf(j,:) 
              END DO
              DEALLOCATE(periapsis_distance_pdf)
           ELSE
              CALL errorMessage("oorb / periapsis_distance", &
                   "Input orbits do not contain sampled " // &
                   "uncertainty information currently required for this task.", 1)        
              STOP
              !              orb = getNominalOrbit(storb_arr_in(i))
              !              cov = getCovarianceMatrix(storb_arr_in(i), element_type_out_prm)
              !              CALL writeOpenOrbOrbitFile(lu_orb_out, i==1, &
              !                   element_type_out_prm, id_arr_in(i), &
              !                   orb=orb, cov=cov, H=HG_arr_in(i,1), &
              !                   G=HG_arr_in(i,2))
           END IF
        END DO
     ELSE
        CALL errorMessage("oorb / periapsis_distance", &
             "Input orbits do not contain sampled " // &
             "uncertainty information currently required for this task.", 1)        
        STOP
     END IF


  CASE ("moid")

     ! Compute Earth MOID (w/o uncertainty) for input orbits and write
     ! a DES orbit file with the results.

     DO i=1,SIZE(orb_arr_in)
        elements = getElements(orb_arr_in(i), "cometary")
        IF (elements(2) > 1.0_bp) THEN
           ! MOID for hyperbolic orbits currently not computed...
           moid = -1.0_bp
        ELSE
           err_verb_ = err_verb 
           err_verb = 0
           ! Epoch of the asteroid orbit:
           epoch = getTime(orb_arr_in(i))
           mjd_tt = getMJD(epoch, "TT")
           ! Earth's osculating elements for the epoch of the asteroid orbit:
           planeph => JPL_ephemeris(mjd_tt, 3, 11, error)
           CALL NEW(ref_orb, planeph(1,:), "cartesian", "equatorial", epoch)
           ! MOID between Earth and asteroid orbits:
           moid = getMOID(orb_arr_in(i), ref_orb)
           CALL NULLIFY(epoch)
           CALL NULLIFY(ref_orb)
           DEALLOCATE(planeph, stat=err)
           err_verb = err_verb_
           IF (error) THEN
              moid = -1.0_bp
              error = .FALSE.
           END IF
        END IF
        IF (orbit_format_out == "des") THEN
           CALL writeDESOrbitFile(lu_orb_out, i==1, "cometary", &
                id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), 1, 6, &
                moid, "OPENORB")
        ELSE IF (orbit_format_out == "orb") THEN
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
                H=HG_arr_in(i,1), G=HG_arr_in(i,3))
        END IF
        IF (error) THEN
           WRITE(stderr,*) i
           error = .FALSE.
        END IF
     END DO
     CALL NULLIFY(orb_out_file)


  CASE ("obsplanner")

     CALL readConfigurationFile(conf_file, &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          integrator=integrator, &
          integration_step=integration_step)

     ! Input observatory code
     obsy_code = get_cl_option("--code=", obsy_code)

     ! Input evolutionary timespan [days]
     timespan = get_cl_option("--timespan=", 0.0_bp)

     ! Input time step [days]
     step = get_cl_option("--step=", 1.0_bp)
     IF (step == 0.0_bp) THEN
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

     ! Input how long the requirements have to be fulfilled each night
     dt_fulfill_night = get_cl_option("--timespan-fulfill-night=", 0.05_bp)
     minstep = MAX(1,CEILING(dt_fulfill_night/step))
     ALLOCATE(temp_arr(minstep,35))

     ! Input minimum altitude for object (wrt the horizon)
     obj_alt_min = get_cl_option("--object-altitude-min=", 30.0_bp)
     obj_alt_min = obj_alt_min*rad_deg

     ! Input maximum apparent magnitude (minimum brightness) for object
     obj_vmag_max = get_cl_option("--object-magnitude-max=", 23.5_bp)

     ! Input maximum altitude for Sun
     solar_alt_max = get_cl_option("--solar-altitude-max=", -18.0_bp)
     solar_alt_max = solar_alt_max*rad_deg

     ! Input maximum solar elongation
     solar_elon_max = get_cl_option("--solar-elongation-max=", 180.0_bp)
     solar_elon_max = solar_elon_max*rad_deg

     ! Input minimum solar elongation
     solar_elon_min = get_cl_option("--solar-elongation-min=", 45.0_bp)
     solar_elon_min = solar_elon_min*rad_deg

     ! Input maximum altitude for Moon
     lunar_alt_max = get_cl_option("--lunar-altitude-max=", 90.0_bp)
     lunar_alt_max = lunar_alt_max*rad_deg

     ! Input minimum distance to Moon
     lunar_elongation_min = get_cl_option("--lunar-elongation-min=", 45.0_bp)
     lunar_elongation_min = lunar_elongation_min*rad_deg

     ! Input minimum phase of Moon
     lunar_phase_min = get_cl_option("--lunar-phase-min=", 0.25_bp)

     ! Input maximum phase of Moon
     lunar_phase_max = get_cl_option("--lunar-phase-max=", 0.75_bp)

     CALL NEW(obsies)
     IF (error) THEN
        CALL errorMessage('oorb / obsplanner', &
             'TRACE BACK (5)',1)
        STOP
     END IF

     ! Write header

     IF (ALLOCATED(storb_arr_in)) THEN

        CALL errorMessage("oorb / obsplanner", &
             "storb option not available yet", 1)
        STOP

     ELSE

        DO i=1,norb
           t = getTime(orb_arr_in(i))
           mjd_tt = getMJD(t, "TT")
           ALLOCATE(observers(nstep))
           DO j=1,nstep
              ! Compute heliocentric observatory coordinates
              observers(j) = getObservatoryCCoord(obsies, obsy_code, t)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (10)',1)
                 STOP
              END IF
              CALL rotateToEquatorial(observers(j))
              CALL NULLIFY(t)
              CALL NEW(t, mjd_tt+j*step, "TT")
              IF (error) THEN
                 CALL errorMessage("oorb / obsplanner", &
                      "TRACE BACK (15)", 1)
                 STOP
              END IF
           END DO

           ! Set integration parameters
           CALL setParameters(orb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / obsplanner", &
                   "TRACE BACK (15)", 1)
              STOP
           END IF

           ! Compute topocentric ephemerides
           CALL getEphemerides(orb_arr_in(i), observers, ephemerides, &
                this_lt_corr_arr=orb_lt_corr_arr)
           IF (error) THEN
              CALL errorMessage('oorb / obsplanner', &
                   'TRACE BACK (20)',1)
              STOP
           END IF

           IF (separately) THEN
              CALL NEW(tmp_file, TRIM(id_arr_in(i)) // ".eph")
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (25)',1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (30)',1)
                 STOP
              END IF
              lu = getUnit(tmp_file)
           ELSE
              lu = stdout
           END IF

           IF (separately .OR. i == 1) THEN
              WRITE(lu,'(A3,1X,A11,1X,A,1X,35(A18,1X))') "#RN", &
                   "Designation", "Code", "MJD_UTC/UT1", "Delta", &
                   "RA", "Dec", "dDelta/dt", "dRA/dt", "dDec/dt", &
                   "VMag", "Alt", "Phase", "LunarElon", "LunarAlt", &
                   "LunarPhase", "SolarElon", "SolarAlt", "r", &
                   "HLon", "HLat", "TLon", "TLat", "TOCLon", &
                   "TOCLat", "HOCLon", "HOCLat", "TOppLon", &
                   "TOppLat", "HEclObj X", "HEclObj Y", "HEclObj Z", &
                   "HEclObj dX/dt", "HEclObj dY/dt", & 
                   "HEclObj dZ/dt", "HEclObsy X", "HEclObsy Y", &
                   "HEclObsy Z" 
           END IF

           istep = 0
           DO j=1,nstep

              istep = istep + 1
              t = getTime(observers(j))
              mjd_tt = getMJD(t, "TT")
              err_verb_ = err_verb
              err_verb = 0
              mjd_utc = getMJD(t, "UTC")
              err_verb = err_verb_
              IF (error) THEN
                 error = .FALSE.
                 mjd_utc = getMJD(t, "UT1")
                 IF (error) THEN
                    CALL errorMessage('oorb / obsplanner', &
                         'TRACE BACK (30)',1)
                    STOP
                 END IF
              END IF

              ! Extract topocentic equatorial coordinates
              CALL rotateToEquatorial(ephemerides(j))        
              comp_coord = getCoordinates(ephemerides(j))
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (35)',1)
                 STOP
              END IF
              Delta = comp_coord(1)
              ra = comp_coord(2)
              dec = comp_coord(3)

              ! Extract instantaneous topocentric equatorial coordinate velocities
              dDelta = comp_coord(4)
              dra = comp_coord(5)*COS(comp_coord(3))
              ddec = comp_coord(6)

              ! Extract topocentric ecliptic lon and lat
              CALL rotateToEcliptic(ephemerides(j))        
              comp_coord = getCoordinates(ephemerides(j))
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (40)',1)
                 STOP
              END IF
              tlon = comp_coord(2)
              tlat = comp_coord(3)

              ! Compute topocentric opposition coordinates (actually,
              ! get heliocentric observatory coordinates)
              scoord = getSCoord(observers(j))
              CALL rotateToEcliptic(scoord)
              pos_opp = getPosition(scoord)
              opplon = pos_opp(2)
              opplat = pos_opp(3)
              CALL NULLIFY(scoord)

              ! Compute opposition-centered topocentric ecliptic coordinates
              toclon = tlon - opplon
              toclat = tlat - opplat
              IF (toclon > pi) THEN
                 toclon = toclon - two_pi
              ELSE IF (toclon < -pi) THEN
                 toclon = toclon + two_pi
              END IF

              ! Extract heliocentric ecliptic lon and lat
              ccoord = getCCoord(orb_lt_corr_arr(j), "ecliptic")
              scoord = getSCoord(ccoord)
              comp_coord = getCoordinates(scoord)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (45)',1)
                 STOP
              END IF
              hlon = comp_coord(2)
              hlat = comp_coord(3)
              CALL NULLIFY(ccoord)
              CALL NULLIFY(scoord)

              ! Compute opposition-centered topocentric ecliptic coordinates
              hoclon = hlon - opplon
              hoclat = hlat - opplat
              IF (hoclon > pi) THEN
                 hoclon = hoclon - two_pi
              ELSE IF (hoclon < -pi) THEN
                 hoclon = hoclon + two_pi
              END IF

              ! Compute phase angle
              CALL NEW(ccoord, ephemerides(j))
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (50)',1)
                 STOP
              END IF
              CALL rotateToEquatorial(ccoord)
              obsy_obj = getPosition(ccoord)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (55)',1)
                 STOP
              END IF
              ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
              CALL toCartesian(orb_lt_corr_arr(j), frame='equatorial')
              pos = getPosition(orb_lt_corr_arr(j))
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (60)',1)
                 STOP
              END IF
              heliocentric_r2 = DOT_PRODUCT(pos,pos)
              obsy_pos = getPosition(observers(j))
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (65)',1)
                 STOP
              END IF
              observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)
              cos_obj_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
                   observer_r2) / (SQRT(heliocentric_r2) * &
                   SQRT(ephemeris_r2))
              obj_phase = ACOS(cos_obj_phase)

              ! Compute apparent brightness:
              ! Input slope parameter
              G_value = get_cl_option("--G=", HG_arr_in(i,3))
              obj_vmag = getApparentMagnitude(H=HG_arr_in(i,1), &
                   G=G_value, r=SQRT(heliocentric_r2), &
                   Delta=Delta, phase_angle=obj_phase)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (70)',1)
                 STOP
              END IF
              IF (obj_vmag > obj_vmag_max) THEN
                 istep = 0
                 CYCLE
              END IF

              ! Compute (approximate) altitude of the target
              obsy_ccoord = getGeocentricObservatoryCCoord(obsies, obsy_code, t)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (75)',1)
                 STOP
              END IF
              CALL rotateToEquatorial(obsy_ccoord)
              geoc_obsy = getPosition(obsy_ccoord)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (80)',1)
                 STOP
              END IF
              vec3 = cross_product(geoc_obsy,obsy_obj)
              obj_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_obj))
              IF (obj_alt < obj_alt_min) THEN
                 istep = 0
                 CYCLE
              END IF

              ! Compute (approximate) altitude of the Sun
              ! Position of the geocenter as seen from the Sun:
              planeph => JPL_ephemeris(mjd_tt, 3, 11, error)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (85)',1)
                 STOP
              END IF
              ! Position of the Sun as seen from the observatory:
              obsy_sun = -(planeph(1,1:3) + geoc_obsy)
              DEALLOCATE(planeph)
              vec3 = cross_product(geoc_obsy,obsy_sun)
              solar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_sun))
              IF (solar_alt > solar_alt_max) THEN
                 istep = 0
                 CYCLE
              END IF
              ! Compute the solar elongation:
              vec3 = cross_product(obsy_obj,obsy_sun)
              solar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_sun))
              IF (solar_elongation < solar_elon_min) THEN
                 istep = 0
                 CYCLE
              ELSE IF (solar_elongation > solar_elon_max) THEN
                 istep = 0
                 CYCLE
              END IF

              ! Compute phase of the Moon:
              ! Position of the geocenter as seen from the Moon:
              planeph => JPL_ephemeris(mjd_tt, 3, 10, error)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (90)',1)
                 STOP
              END IF
              ! Position of the Moon as seen from the observatory:
              obsy_moon = -(planeph(1,1:3) + geoc_obsy)
              DEALLOCATE(planeph)
              ! Position of the Sun as seen from the Moon:
              planeph => JPL_ephemeris(mjd_tt, 10, 11, error)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (95)',1)
                 STOP
              END IF
              sun_moon = planeph(1,1:3)
              DEALLOCATE(planeph)
              obsy_moon_r2 = DOT_PRODUCT(obsy_moon,obsy_moon)
              sun_moon_r2 = DOT_PRODUCT(sun_moon,sun_moon)
              cos_obj_phase = 0.5_bp * (sun_moon_r2 + obsy_moon_r2 - &
                   observer_r2) / (SQRT(sun_moon_r2) * &
                   SQRT(obsy_moon_r2))
              lunar_phase = (pi-ACOS(cos_obj_phase))/pi
              IF (lunar_phase < lunar_phase_min) THEN
                 istep = 0
                 CYCLE
              ELSE IF (lunar_phase > lunar_phase_max) THEN
                 istep = 0
                 CYCLE
              END IF
              ! Compute (approximate) distance between the target and the Moon:
              vec3 = cross_product(obsy_obj,obsy_moon)
              lunar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_moon))
              IF (lunar_elongation < lunar_elongation_min) THEN
                 istep = 0
                 CYCLE
              END IF
              ! Compute (approximate) altitude of the Moon:
              vec3 = cross_product(geoc_obsy,obsy_moon)
              lunar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_moon))
              IF (lunar_alt > lunar_alt_max) THEN
                 istep = 0
                 CYCLE
              END IF

              ! Extract heliocentric distance
              hdist = SQRT(heliocentric_r2)

              ! Extract heliocentric ecliptic cartesian coordinates for the object
              h_ecl_car_coord_obj = getElements(orb_lt_corr_arr(j), "cartesian", "ecliptic")
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (90)',1)
                 STOP
              END IF

              ! Extract heliocentric ecliptic cartesian coordinates for the observer
              CALL rotateToEcliptic(observers(j))
              h_ecl_car_coord_obsy = getCoordinates(observers(j))
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (95)',1)
                 STOP
              END IF

              CALL NULLIFY(observers(j))
              CALL NULLIFY(ephemerides(j))
              CALL NULLIFY(orb_lt_corr_arr(j))
              CALL NULLIFY(ccoord)
              CALL NULLIFY(obsy_ccoord)

              ra = ra/rad_deg
              dec = dec/rad_deg
              dDelta = dDelta
              dra = dra/rad_deg
              ddec = ddec/rad_deg
              obj_alt = obj_alt/rad_deg
              obj_phase = obj_phase/rad_deg
              lunar_elongation = lunar_elongation/rad_deg
              lunar_alt = lunar_alt/rad_deg
              solar_elongation = solar_elongation/rad_deg
              solar_alt = solar_alt/rad_deg
              tlon = tlon/rad_deg
              tlat = tlat/rad_deg
              toclon = toclon/rad_deg
              toclat = toclat/rad_deg
              hlon = hlon/rad_deg
              hlat = hlat/rad_deg
              hoclon = hoclon/rad_deg
              hoclat = hoclat/rad_deg
              opplon = opplon/rad_deg
              opplat = opplat/rad_deg

              IF (istep <= minstep) THEN
                 temp_arr(istep,:) = &
                      (/ mjd_utc, Delta, &
                      ra, dec, dDelta, dra, ddec, obj_vmag, obj_alt, &
                      obj_phase, lunar_elongation, lunar_alt, &
                      lunar_phase, solar_elongation, solar_alt, hdist, &
                      hlon, hlat, tlon, tlat, toclon, toclat, hoclon, &
                      hoclat, opplon, opplat, h_ecl_car_coord_obj, &
                      h_ecl_car_coord_obsy(1:3) /)
              END IF

              IF (istep == minstep) THEN
                 DO k=1,istep
                    WRITE(lu,'(I0,1X,2(A,1X),35(F18.10,1X))') &
                         i, id_arr_in(i), TRIM(obsy_code), temp_arr(k,:)
                 END DO
              ELSE IF (istep > minstep) THEN
                 WRITE(lu,'(I0,1X,2(A,1X),35(F18.10,1X))') &
                      i, id_arr_in(i), TRIM(obsy_code), mjd_utc, Delta, &
                      ra, dec, dDelta, dra, ddec, obj_vmag, obj_alt, &
                      obj_phase, lunar_elongation, lunar_alt, &
                      lunar_phase, solar_elongation, solar_alt, hdist, &
                      hlon, hlat, tlon, tlat, toclon, toclat, hoclon, &
                      hoclat, opplon, opplat, h_ecl_car_coord_obj, &
                      h_ecl_car_coord_obsy(1:3)
              END IF

           END DO

           IF (separately) THEN
              CALL NULLIFY(tmp_file)
           END IF
           DEALLOCATE(observers, ephemerides, orb_lt_corr_arr)

        END DO

     END IF

     DEALLOCATE(temp_arr)


  CASE ("tai2utc")

     mjd_tai = get_cl_option("--epoch-mjd-tai=", 0.0_bp)
     CALL NEW(epoch, mjd_tai, "TAI")
     IF (error) THEN
        CALL errorMessage("oorb / tai2utc", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     mjd_utc = getMJD(epoch, "UTC")
     IF (error) THEN
        CALL errorMessage("oorb / tai2utc", &
             "TRACE BACK (10)", 1)
        STOP
     END IF
     WRITE(stdout,'(F15.8)') mjd_utc
     CALL NULLIFY(epoch)


  CASE ("utc2tai")

     mjd_utc = get_cl_option("--epoch-mjd-utc=", 0.0_bp)
     CALL NEW(epoch, mjd_utc, "UTC")
     IF (error) THEN
        CALL errorMessage("oorb / utc2tai", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     mjd_tai = getMJD(epoch, "TAI")
     IF (error) THEN
        CALL errorMessage("oorb / utc2tai", &
             "TRACE BACK (10)", 1)
        STOP
     END IF
     WRITE(stdout,'(F15.8)') mjd_tai
     CALL NULLIFY(epoch)


  CASE ("degreestosexagesimal")

     ! Input RA in degrees
     IF (get_cl_option("--ra-degrees=", .FALSE.)) THEN
        ra = get_cl_option("--ra-degrees=", -1.0_bp)
        IF (ra < 0.0_bp) THEN
           CALL errorMessage("oorb / degreestosexagesimal", &
                "Negative degrees not possible for RA.", 1)        
           STOP
        END IF
        ra = ra*rad_deg
        CALL radiansToHMS(ra, i, j, sec)
        WRITE(stdout,"(A,1X,2(I0,1X),F15.12)") "RA = ", i, j, sec
     END IF

     ! Input Dec in degrees
     IF (get_cl_option("--dec-degrees=", .FALSE.)) THEN     
        dec = get_cl_option("--dec-degrees=", -1.0_bp)
        IF (dec < -90.0_bp .OR. dec > 90.0_bp) THEN
           CALL errorMessage("oorb / degreestosexagesimal", &
                "Declination not in the interval [-90 deg,90 deg].", 1)        
           STOP
        END IF
        dec = dec*rad_deg
        CALL radiansToDAMAS(dec, i, j, sec)
        WRITE(stdout,"(A,1X,2(I0,1X),F15.12)") "Dec = ", i, j, sec
     END IF

  CASE ("encode_designation")

     IF (.NOT.(get_cl_option("--mpc", .FALSE.) .OR. &
          get_cl_option("--mpc3", .FALSE.))) THEN     
        CALL errorMessage("oorb / encode_designation", &
             "Either the '--mpc' or the '--mpc3' option has to be specified.", 1)                
        STOP
     ELSE IF (get_cl_option("--mpc", .FALSE.) .AND. &
          get_cl_option("--mpc3", .FALSE.)) THEN     
        CALL errorMessage("oorb / encode_designation", &
             "Both '--mpc' and '--mpc3' options cannot be used simultaneously.", 1)                
        STOP
     END IF

     IF (get_cl_option("--data-in=", .FALSE.)) THEN     
        tmp_fname = get_cl_option("--data-in="," ")
        CALL NEW(tmp_file, TRIM(tmp_fname))
        CALL setStatusOld(tmp_file)
        CALL OPEN(tmp_file)
        IF (error) THEN
           CALL errorMessage("oorb / encode_designation", &
                "TRACE BACK", 1)        
           STOP
        END IF
        DO 
           READ(getUnit(tmp_file),"(A)",iostat=err) id
           IF (err < 0) THEN
              EXIT
           ELSE IF (err > 0) THEN
              CALL errorMessage("oorb / encode_designation", &
                   "Could not read designation from file.", 1)        
              STOP
           END IF
           IF (get_cl_option("--mpc", .FALSE.)) THEN
              CALL encodeMPCDesignation(id)
           ELSE IF (get_cl_option("--mpc3", .FALSE.)) THEN
              CALL encodeMPC3Designation(id)
           END IF
           WRITE(stdout,"(A)") TRIM(id)
        END DO
        CALL NULLIFY(tmp_file)
     ELSE
        CALL errorMessage("oorb / encode_designation", &
             "No input file given (use '--data-in=' option).", 1)                
     END IF


  CASE ("decode_designation")

     IF (get_cl_option("--data-in=", .FALSE.)) THEN     
        tmp_fname = get_cl_option("--data-in="," ")
        CALL NEW(tmp_file, TRIM(tmp_fname))
        CALL setStatusOld(tmp_file)
        CALL OPEN(tmp_file)
        IF (error) THEN
           CALL errorMessage("oorb / decode_designation", &
                "TRACE BACK", 1)        
           STOP
        END IF
        DO 
           READ(getUnit(tmp_file),"(A)",iostat=err) id
           IF (err < 0) THEN
              EXIT
           ELSE IF (err > 0) THEN
              CALL errorMessage("oorb / decode_designation", &
                   "Could not read designation from file.", 1)        
              STOP
           END IF
           IF (LEN_TRIM(id) > 7) THEN
              CALL decodeMPC3Designation(id)
           ELSE
              CALL decodeMPCDesignation(id)
           END IF
           WRITE(stdout,"(A)") TRIM(id)
        END DO
        CALL NULLIFY(tmp_file)
     END IF

  CASE default

     IF (LEN_TRIM(task) == 0) THEN
        CALL errorMessage("oorb", &
             "No task specified.", 1)
     ELSE
        CALL errorMessage("oorb", &
             "No task such as '" // TRIM(task) // "' available.", 1)
     END IF
     STOP

  END SELECT

  ! Close output orbit file
  CALL NULLIFY(orb_out_file)

  ! Deallocate memory
  IF (ASSOCIATED(id_arr_in)) THEN
     DEALLOCATE(id_arr_in)
  END IF
  IF (ASSOCIATED(orb_arr_in)) THEN
     DO i=1,SIZE(orb_arr_in)
        CALL NULLIFY(orb_arr_in(i))
     END DO
     DEALLOCATE(orb_arr_in)
  END IF
  IF (ALLOCATED(storb_arr_in)) THEN
     DO i=1,SIZE(storb_arr_in)
        CALL NULLIFY(storb_arr_in(i))
     END DO
     DEALLOCATE(storb_arr_in)
  END IF
  IF (ASSOCIATED(HG_arr_in)) THEN
     DEALLOCATE(HG_arr_in)
  END IF
  CALL JPL_ephemeris_nullify()
  CALL nullifyTime()
  IF (ASSOCIATED(perturbers)) THEN
     DEALLOCATE(perturbers, stat=err)
  END IF

END PROGRAM oorb
