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
!! Main program for various tasks that include orbit computation.
!!
!! @author  MG
!! @version 2018-06-12
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
  USE statistics
  USE io

  IMPLICIT NONE
  INCLUDE "version.h"

  TYPE (PhysicalParameters) :: &
       physparam
  TYPE (StochasticOrbit), DIMENSION(:), ALLOCATABLE :: &
       storb_arr_in
  TYPE (StochasticOrbit) :: &
       storb
  TYPE (Orbit), DIMENSION(:,:), POINTER :: &
       orb_lt_corr_arr2 => NULL()
  TYPE (Orbit), DIMENSION(:), POINTER :: &
       orb_arr_in => NULL()
  TYPE (Orbit), DIMENSION(:), POINTER :: &
       orb_arr => NULL(), &
       orb_arr_ => NULL(), &       
       orb_arr_cmp => NULL(), &
       orb_lt_corr_arr => NULL()
  TYPE (Orbit) :: &
       orb, &
       ref_orb
  TYPE (Observations), DIMENSION(:), POINTER :: &
       obss_sep => NULL()
  TYPE (Observations) :: &
       obss, &
       obss_in, &
       obss_in_temp
  TYPE (Observation) :: &
       obs
  TYPE (Observation), DIMENSION(:), POINTER :: &
       obs_arr => NULL()
  TYPE (Observatories) :: &
       obsies
  TYPE (SphericalCoordinates), DIMENSION(:,:), POINTER :: &
       ephemerides_arr => NULL()
  TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: &
       ephemerides => NULL()
  TYPE (SphericalCoordinates) :: &
       ephemeris, &
       scoord
  TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: &
       observers => NULL()
  TYPE (CartesianCoordinates) :: &
       ccoord, &
       obsy_ccoord, &
       sun_ccoord
  TYPE (Time) :: &
       epoch, &
       epoch0, &
       epoch1, &
       t
  TYPE (File) :: &
       conf_file, &                                                 !! Configuration file.
       obs_file, &                                                  !! Generic observation file.
       orb_in_file, &                                               !! Input orbit file.
       orb_out_file, &                                              !! Output orbit file.
       out_file, &                                                  !! Generic output file.
       tmp_file, &                                                     !! Generic temporary file.
       tmp_file2
  CHARACTER(len=ELEMENT_TYPE_LEN), DIMENSION(:), ALLOCATABLE :: &
       element_type_pdf_arr_in                                      !! Element type of input orbital-element PDF.
  CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), POINTER :: &
       obsy_code_arr => NULL()                                      !! IAU/MPC designated observatory code.
  CHARACTER(len=DESIGNATION_LEN), DIMENSION(:), POINTER :: &
       group_name_arr => NULL(), &
       id_arr_storb_in => NULL(), &
       id_arr_in => NULL(), &
       id_arr => NULL()
  CHARACTER(len=2), DIMENSION(:), POINTER :: &
       filters => NULL()
  CHARACTER(len=1), DIMENSION(:), ALLOCATABLE :: &
       strlen1
  CHARACTER(len=32), DIMENSION(:), ALLOCATABLE :: &
       str_arr
  CHARACTER(len=32), DIMENSION(6) :: &
       element_str_arr, &
       stdev_str_arr
  CHARACTER(len=32), DIMENSION(5) :: &
       corr_str_arr
  CHARACTER(len=1024), DIMENSION(4) :: &
       header                                                       !! Generic header.
  CHARACTER(len=4096) :: &
       str1
  CHARACTER(len=FNAME_LEN) :: &
       conf_fname, &                                                !! Path to configuration file (incl. fname). 
       planetary_ephemeris_fname, &
       gnuplot_scripts_dir, &                                       !! Path to Gnuplot scripts directory. 
       obs_fname, &                                                 !! Path to observation file (incl. fname).
       orb_in_fname, &                                              !! Path to input orbit file (incl. fname).
       orb_out_fname, &                                             !! Path to output orbit file (incl. fname).
       out_fname, &                                                 !! Path to generic output file (incl. fname).
       tmp_fname, &
       buffer, &
       line
  CHARACTER(len=FNAME_LEN) , DIMENSION(:), ALLOCATABLE :: &
       fnames
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
  CHARACTER(len=FRAME_LEN) :: &
       frame, &
       frame_
  CHARACTER(len=256) :: &
       flavor, &
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
       cov_arr => NULL(), &
       cov_arr_in => NULL(), &                                                !! Input array of covariance matrices.
       encounters => NULL(), &
       HG_arr_storb_in => NULL()
  REAL(bp), DIMENSION(:,:), POINTER :: &
       apoapsis_distance_pdf => NULL(), &
       HG_arr_in => NULL(), &
       jac_arr_cmp => NULL(), &
       pdfs_arr => NULL(), &
       periapsis_distance_pdf => NULL(), &
       planeph => NULL(), &
       vov_map => NULL(), &
       vomcmc_map => NULL(), &
       sor_rho_arr => NULL()
  REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
       elements_arr, &
       ephem_, &
       hist, &
       jac_arr_in, &
       mean_arr, &
       temp_arr, &
       histo1, histo2
  REAL(bp), DIMENSION(6,6) :: &
       corr, &
       cov, &
       eigenvectors
  REAL(bp), DIMENSION(6,2) :: &
       vov_scaling_cmp, &
       vov_scaling_prm, &
       vomcmc_scaling_cmp, &
       vomcmc_scaling_prm
  REAL(bp), DIMENSION(2,2) :: &
       sor_rho_cmp
  REAL(bp), DIMENSION(:), POINTER :: &
       mags => NULL(), &
       pdf_arr => NULL(), &
       pdf_arr_cmp => NULL(), &
       pdf_arr_in => NULL(), &
       rchi2_arr_cmp => NULL(), &
       reg_apr_arr_cmp => NULL(), &
       weight_arr => NULL()
  REAL(bp), DIMENSION(:), ALLOCATABLE :: &
       arc_arr, &
       rchi2_arr_in, &
       reg_apr_arr_in, &
       real_arr
  REAL(bp), DIMENSION(10) :: &
       tisserands_parameters,  &
       jacobi_constants
  REAL(bp), DIMENSION(7) :: &
       ran7
  REAL(bp), DIMENSION(6) :: &
       comp_coord, &
       coordinates, &
       eigenvalues, &
       elements, &
       h_ecl_car_coord_obj, &
       h_ecl_car_coord_obsy, &
       lower_limit, &
       upper_limit, &
       mean, &
       obs_stdev_arr_prm, &
       ran6, &
       stdev_arr
  REAL(bp), DIMENSION(4) :: &
       sor_genwin_offset, &
       sor_rho_init
  REAL(bp), DIMENSION(3) :: &
       dgeoc_obsy, &
       geoc_obsy, &
       moon_obsy, &
       moon_sun, &
       obsy_moon, &
       obsy_obj, &
       obsy_pos, &
       obsy_sun, &
       obsy_vel, &
       pos, &
       pos_opp, &
       pos_sun, &
       sun_moon, &
       pos_ast, &
       vec3, &
       vel
  REAL(bp), DIMENSION(2) :: &
       bounds
  REAL(bp) :: &
       b1, b2, &
       Delta, diameter, &
       H_max, H_value, H10_value, H, Hmin, Hmax, &
       G_value, G12, &
       accwin_multiplier, amin, amax, &
       angular_distance_dec, angular_distance_ra, &
       apoapsis_distance, apriori_a_max, &
       apriori_a_min, apriori_apoapsis_max, apriori_apoapsis_min, &
       apriori_periapsis_max,  apriori_periapsis_min, &
       apriori_rho_min, apriori_rho_max, &
       caldate, chi2_min_init, cos_nsigma, cos_obj_phase, &
       day0, day1, dDelta, ddec, dec, dra, dt, dt_, dt_fulfill_night, dchi2_max, &
       emin, emax, ephemeris_r2, &
       geometric_albedo, &
       hdist, heliocentric_r2, hlat, hlon, hoclat, hoclon, &
       i_min, i_max, imin, imax, integration_step, integration_step_init, &
       ls_correction_factor, ls_rchi2_acceptable, lunar_alt, lunar_alt_max, &
       lunar_elongation, lunar_elongation_min, lunar_phase, &
       lunar_phase_min, lunar_phase_max, &
       mag, mjd, mjd0, mjd1, mjd_tai, mjd_tt, mjd_utc, mjd_utc0, moid, &
       moon_obsy_r2, moon_sun_r2, &
       n, &
       obj_alt, obj_alt_min, obj_phase, obj_vmag, obj_vmag_cometary, obj_vmag_max, &
       observer_r2, obsy_moon_r2, opplat, opplon, outlier_multiplier_prm, &
       output_interval, &
       pa, periapsis_distance, peak, pp_G, pp_G_unc, probability_mass, &
       ra, &
       sec, smplx_tol, smplx_similarity_tol, &
       solar_elongation, solar_elon_min, solar_elon_max, &
       solar_alt, solar_alt_max, solelon_max, solelon_min, generat_multiplier, stdev, &
       step, sun_moon_r2, sunlat, sunlon, &
       timespan, tlat, tlon, toclat, toclon, tsclat, tsclon, &
       ta_s, ta_c, fak, ecc_anom, true_anom, chi2, chi2_min
  INTEGER, DIMENSION(:), POINTER :: &
       repetition_arr_cmp => NULL(), &
       repetition_arr_in => NULL()
  INTEGER, DIMENSION(:), ALLOCATABLE :: &
       indx_arr, &
       int_arr
  INTEGER :: &
       cos_norb, cos_ntrial, &
       err, &
       err_verb_, &
       i, &
       iday, &
       indx, indx_max, indx_min, &
       iobj, &
       iorb, &
       istep, &
       j, &
       k, k_max, &
       l, &
       m, &
       ls_niter_major_max, &
       ls_niter_major_min, &
       ls_niter_minor, &
       lu, &
       lu_orb_out, &
       min, minstep, &
       month, month0, month1, &
       center, nhist, nobj, nobs, norb, noutlier, nstep, &
       os_norb, os_ntrial, os_sampling_type, &
       smplx_niter, &
       sor_niter, sor_norb, sor_norb_sw, sor_ntrial, sor_ntrial_sw, &
       sor_type_prm, &
       vov_type, vov_type_prm, vov_norb, vov_ntrial, vov_niter, &
       vov_norb_iter, vov_ntrial_iter, vov_nmap, &
       vomcmc_type, vomcmc_type_prm, vomcmc_norb, vomcmc_ntrial, vomcmc_niter, &
       vomcmc_norb_iter, vomcmc_ntrial_iter, vomcmc_nmap, &
       year, year0, year1, &
       loc, nfile
  LOGICAL, DIMENSION(:,:), POINTER :: &
       obs_masks => NULL()
  LOGICAL, DIMENSION(:), POINTER :: &
       perturbers => NULL()
  LOGICAL, DIMENSION(6) :: &
       ls_element_mask, &
       vov_mapping_mask, &
       vomcmc_mapping_mask
  LOGICAL, DIMENSION(4) :: &
       sor_iterate_bounds
  LOGICAL :: &
       reading, &
       compress, &
       cos_gaussian, &
       first, &
       force_earth_impact_at_epoch, &
       gaussian_rho, &
       generat_gaussian_deviates, &
       interval, &
       mjd_epoch, &
       noise, &
       outlier_rejection_prm, &
       plot_open, &
       plot_results, &
       relativity, &
       pp_H_estimation, &
       radians, &
       random_obs, &
       regularized, &
       separately, separately_, & !! Output orbit(s)/ephemerides/etc separately for each object
       smplx_force, &
       dchi2_rejection, &
       write_residuals, &
       asteroid_perturbers
  ! Defaults:
  error = .FALSE.
  task = " "
  compress = .FALSE.
  orbit_format_out = "des"
  element_type_comp_prm = "keplerian"
  element_type_out_prm = "keplerian"
  planetary_ephemeris_fname = "de430.dat"
  err_verb = 1
  info_verb = 1
  frame = "ecliptic"
  gnuplot_scripts_dir = "."
  obs_stdev_arr_prm = -1.0_bp
  observation_format_out = "des"
  obsy_code = "500"
  dyn_model = " "
  integrator = " "
  integration_step = -1.0_bp
  orb_in_fname = " "
  orb_out_fname = " "
  separately = .FALSE.
  simint = 1
  outlier_multiplier_prm = 3.0_bp
  outlier_rejection_prm = .FALSE.
  plot_open = .FALSE.
  plot_results = .FALSE.
  radians = .FALSE.
  relativity = .TRUE.
  pp_H_estimation = .FALSE.
  pp_G = 99.9_bp
  pp_G_unc = 99.9_bp
  mjd_epoch = .TRUE.
  smplx_force = .FALSE.
  write_residuals = .FALSE.

  IF (get_cl_option("--version",.FALSE.)) THEN
     WRITE(stdout,"(A)") ""
     WRITE(stdout,"(A)") "OpenOrb v" // VERSION
     WRITE(stdout,"(A)") "Copyright 2011 Mikael Granvik, Jenni Virtanen, Karri Muinonen,"
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
     WRITE(stdout, "(A)") ""
     STOP
  END IF

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
  CALL setActionRead(conf_file)
  CALL setStatusOld(conf_file)
  CALL OPEN(conf_file)
  IF (error) THEN
     CALL errorMessage("oorb", &
          "Configuration file missing or error while opening it.", 1)
     STOP
  END IF
  CALL readConfigurationFile(conf_file, &
       planetary_ephemeris_fname=planetary_ephemeris_fname, &
       err_verb=err_verb, &
       info_verb=info_verb, &
       element_type_comp=element_type_comp_prm, &
       element_type_out=element_type_out_prm, &
       obs_stdev_arr=obs_stdev_arr_prm, &
       observation_format_out=observation_format_out, &
       orbit_format_out=orbit_format_out, &
       outlier_rejection=outlier_rejection_prm, &
       outlier_multiplier=outlier_multiplier_prm, &
       dchi2_rejection=dchi2_rejection, &
       dchi2_max=dchi2_max, &
       plot_open=plot_open, &
       plot_results=plot_results, &
       dyn_model=dyn_model, &
       perturbers=perturbers, &
       integrator=integrator, &
       integration_step=integration_step, &
       relativity=relativity, &
       simint=simint, &
       pp_H_estimation=pp_H_estimation, &
       pp_G=pp_G, &
       pp_G_unc=pp_G_unc, &
       asteroid_perturbers=asteroid_perturbers  )
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
  CALL set_relativity(relativity)

  ! Set path to data files:
  CALL setAccessToDataFiles()
  IF (LEN_TRIM(planetary_ephemeris_fname) == 0) THEN
     planetary_ephemeris_fname = TRIM(EPH_FNAME)
  END IF
  CALL JPL_ephemeris_init(error, &
       filename=TRIM(OORB_DATA_DIR) // "/" // TRIM(planetary_ephemeris_fname)) 
  IF (error) THEN
     CALL errorMessage("oorb", &
          "Could not initialize planetary ephemerides using the " // &
          TRIM(OORB_DATA_DIR) // "/" // TRIM(planetary_ephemeris_fname) // " file.", 1)
     STOP
  END IF

  ! Set path to Gnuplot scripts using environment variable:
  CALL getenv("OORB_GNUPLOT_SCRIPTS_DIR", gnuplot_scripts_dir)

  ! Read observation file if given:
  obs_fname = get_cl_option("--obs-in="," ")
  IF (LEN_TRIM(obs_fname) /= 0) THEN
     loc = 0
     nfile = 1         ! Need to start at 1 because the while loop below
     buffer = obs_fname ! won't take the last one into account.
     reading = .TRUE.
     DO WHILE (reading .EQV. .TRUE.)      ! Check amount of commas (ie amount of obs files)
        loc = INDEX(buffer, ",")
        IF (loc == 0) THEN
           reading = .FALSE.
           CYCLE
        END IF
        buffer = buffer(loc+1:)
        nfile = nfile + 1
     END DO
     ! Now that we know the amount of files, get each individual name
     buffer = obs_fname
     ALLOCATE(fnames(nfile))
     DO i=1, nfile-1
        loc = INDEX(buffer, ",")
        fnames(i) = buffer(1:loc-1)
        buffer = buffer(loc+1:)
     END DO
     fnames(nfile) = buffer
     ! Now that we have the filenames, read them
     DO i=1, SIZE(fnames)
        CALL NEW(obs_file, TRIM(fnames(i)))
        CALL setActionRead(obs_file)
        CALL setStatusOld(obs_file)
        CALL OPEN(obs_file)
        IF (error) THEN
           CALL errorMessage("oorb", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (i == 1) THEN
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
        ELSE
           IF (ANY(obs_stdev_arr_prm < 0.0_bp)) THEN
              CALL NEW(obss_in_temp, obs_file)
           ELSE
              CALL NEW(obss_in_temp, obs_file, stdev=obs_stdev_arr_prm)
           END IF
           obss_in = obss_in + obss_in_temp
           CALL NULLIFY(obss_in_temp)
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "TRACE BACK (35)", 1)
              STOP
           END IF
        END IF

        CALL NULLIFY(obs_file)
        ! Now we come up with the output filenames.
        indx = INDEX(obs_fname,".",back=.TRUE.)
        out_fname = obs_fname(1:indx-1)
        reading = .TRUE.
        DO WHILE (reading .EQV. .TRUE.) ! Remove commas from filename
           indx = INDEX(out_fname,",",back=.TRUE.)
           IF (indx == 0) THEN
              reading = .FALSE.
           ELSE
              out_fname(indx:indx) = "_"
           END IF
        END DO
        reading = .TRUE.
        j = 1
        DO WHILE (reading .EQV. .TRUE.) ! Remove backslashes from filename

           indx = INDEX(out_fname,"/",back=.TRUE.)
           IF (indx == 0) THEN
              reading = .FALSE.

           ELSE
              out_fname(indx:indx) = "_"
           END IF
        END DO
     END DO
     obss_sep => getSeparatedSets(obss_in)
     IF (error) THEN
        CALL errorMessage("oorb", &
             "TRACE BACK (37)", 1)
        STOP
     END IF
  ELSE IF (task == "2mpc3" .OR. &
       task == "2mpc" .OR. &
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
     CALL setActionRead(orb_in_file)
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
             reg_apr_arr_in(norb), repetition_arr_in(norb), stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb", &
                "Could not allocate memory (20).", 1)
           STOP
        END IF
        jac_arr_in = -1.0_bp
        pdf_arr_in = -1.0_bp
        cov_arr_in = -1.0_bp
        repetition_arr_in = 0
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
                jac_equ_kep=jac_arr_in(i,3), &
                repetitions=repetition_arr_in(i))
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "Could not read orbit file.", 1)
              STOP
           END IF
           !DO WHILE (LEN_TRIM(id_arr_in(i)) < 7)
           !   id_arr_in(i) = "0" // TRIM(id_arr_in(i))
           !END DO
        END DO
        CALL NULLIFY(orb_in_file)

        ! Calculate the number of different objects in the orbit file:
        ALLOCATE(indx_arr(SIZE(id_arr_in)))
        CALL quickSort(id_arr_in, indx_arr, errstr)
        IF (LEN_TRIM(errstr) /= 0) THEN
           CALL errorMessage("oorb", &
                "Could not sort orbit id's. " // TRIM(errstr), 1)
           STOP
        END IF
        nobj = 1
        DO i=1,SIZE(id_arr_in)-1
           IF (id_arr_in(indx_arr(i)) /= id_arr_in(indx_arr(i+1))) THEN
              nobj = nobj + 1
           END IF
        END DO
        DEALLOCATE(indx_arr)

        IF (info_verb >= 2) THEN
           WRITE(stdout,"(A,1X,I0,1X,A)") &
                "Input orbit file contains orbits for ", nobj, &
                "different objects."
           WRITE(stdout,"(A,1X,2A)") & 
                "The type of input orbital elements is", & 
                TRIM(element_type_in), "."
           WRITE(stdout,"(A,1X,2A)") & 
                "The type of input pdf values is", & 
                TRIM(element_type_pdf_arr_in(1)), "."
        END IF

        ! Initialize stochasticorbits if uncertainty information available:
        IF (norb > nobj .AND. ALL(pdf_arr_in /= -1.0_bp)) THEN
           ! Sampled PDF available
           ALLOCATE(storb_arr_in(nobj))
           i = 0
           j = 1
           ! generate storbs for sets 1...n-1
           DO k=2,SIZE(id_arr_in)
              IF (id_arr_in(k-1) /= id_arr_in(k)) THEN
                 i = i + 1
                 IF (ASSOCIATED(obss_sep)) THEN
                    IF (getID(obss_sep(1)) == id_arr_in(k-1)) THEN
                       CALL NEW(storb_arr_in(i), orb_arr_in(j:k-1), pdf_arr_in(j:k-1), &
                            element_type_pdf_arr_in(1), jac_arr=jac_arr_in(j:k-1,:), &
                            reg_apr_arr=reg_apr_arr_in(j:k-1), &
                            rchi2_arr=rchi2_arr_in(j:k-1), &
                            repetition_arr=repetition_arr_in(j:k-1), &
                            obss=obss_sep(1))
                    END IF
                 END IF
                 IF (.NOT.exist(storb_arr_in(i))) THEN
                    CALL NEW(storb_arr_in(i), orb_arr_in(j:k-1), pdf_arr_in(j:k-1), &
                         element_type_pdf_arr_in(1), jac_arr=jac_arr_in(j:k-1,:), &
                         reg_apr_arr=reg_apr_arr_in(j:k-1), &
                         rchi2_arr=rchi2_arr_in(j:k-1), &
                         repetition_arr=repetition_arr_in(j:k-1))
                 END IF
                 j = k
              END IF
           END DO
           ! generate storb for set n (note k-1 above)
           i = i + 1
           IF (ASSOCIATED(obss_sep)) THEN
              IF (getID(obss_sep(1)) == id_arr_in(j)) THEN
                 CALL NEW(storb_arr_in(i), orb_arr_in(j:), pdf_arr_in(j:), &
                      element_type_pdf_arr_in(1), jac_arr=jac_arr_in(j:,:), &
                      reg_apr_arr=reg_apr_arr_in(j:), &
                      rchi2_arr=rchi2_arr_in(j:), &
                      repetition_arr=repetition_arr_in(j:), &
                      obss=obss_sep(1))
              END IF
           END IF
           IF (.NOT.exist(storb_arr_in(i))) THEN
              CALL NEW(storb_arr_in(i), orb_arr_in(j:), pdf_arr_in(j:), &
                   element_type_pdf_arr_in(1), jac_arr=jac_arr_in(j:,:), &
                   reg_apr_arr=reg_apr_arr_in(j:), &
                   rchi2_arr=rchi2_arr_in(j:), &
                   repetition_arr=repetition_arr_in(j:))
           END IF
           ALLOCATE(id_arr(nobj), HG_arr_storb_in(nobj,norb,4))
           id_arr = ""
           id_arr(1) = id_arr_in(1)
           HG_arr_storb_in = 99.9_bp 
           HG_arr_storb_in(1,1,1:4) = HG_arr_in(1,1:4)
           j = 1
           k = 1
           k_max = 1
           DO i=2,SIZE(id_arr_in)
              k = k + 1
              IF (ALL(id_arr(1:j) /= id_arr_in(i))) THEN
                 j = j + 1
                 id_arr(j) = id_arr_in(i)
                 k = 1
                 HG_arr_storb_in(j,k,1:4) = HG_arr_in(i,1:4)
              ELSE
                 HG_arr_storb_in(j,k,1:4) = HG_arr_in(i,1:4)
              END IF
              IF (k > k_max) THEN
                 k_max = k
              END IF
           END DO
           ALLOCATE(id_arr_storb_in(nobj))
           id_arr_storb_in = id_arr
           HG_arr_storb_in => reallocate(HG_arr_storb_in, nobj, k_max, 4)
           DEALLOCATE(id_arr, id_arr_in, HG_arr_in)
        ELSE IF (nobj == norb .AND. &
             ALL(cov_arr_in(:,1,1) >= 0.0_bp)) THEN
           ! Covariance matrix/ces available
           ALLOCATE(storb_arr_in(nobj))
           i = 0
           j = 1
           ! generate storbs for sets 1...n-1
           DO k=2,SIZE(id_arr_in)
              IF (id_arr_in(k-1) /= id_arr_in(k)) THEN
                 i = i + 1
                 IF (ASSOCIATED(obss_sep)) THEN
                    IF (getID(obss_sep(1)) == id_arr_in(k-1)) THEN
                       CALL NEW(storb_arr_in(i), orb_arr_in(i), cov_arr_in(i,:,:), &
                            cov_type=element_type_in, element_type=element_type_in, &
                            obss=obss_sep(1))
                    END IF
                 END IF
                 IF (.NOT.exist(storb_arr_in(i))) THEN
                    CALL NEW(storb_arr_in(i), orb_arr_in(i), cov_arr_in(i,:,:), &
                         cov_type=element_type_in, element_type=element_type_in, &
                         obss=obss_sep(1))
                 END IF
                 j = k
              END IF
           END DO
           ! generate storb for set n (note k-1 above)
           i = i + 1
           IF (ASSOCIATED(obss_sep)) THEN
              IF (getID(obss_sep(1)) == id_arr_in(j)) THEN
                 CALL NEW(storb_arr_in(i), orb_arr_in(i), cov_arr_in(i,:,:), &
                      cov_type=element_type_in, element_type=element_type_in, &
                      obss=obss_sep(1))
              END IF
           END IF
           IF (.NOT.exist(storb_arr_in(i))) THEN
              CALL NEW(storb_arr_in(i), orb_arr_in(i), cov_arr_in(i,:,:), &
                   cov_type=element_type_in, element_type=element_type_in)
           END IF
           ALLOCATE(id_arr(nobj), HG_arr_storb_in(nobj,norb,4))
           id_arr = ""
           id_arr(1) = id_arr_in(1)
           HG_arr_storb_in = 99.9_bp 
           HG_arr_storb_in(1,1,1:4) = HG_arr_in(1,1:4)
           j = 1
           k = 1
           k_max = 1
           DO i=2,SIZE(id_arr_in)
              k = k + 1
              IF (ALL(id_arr(1:j) /= id_arr_in(i))) THEN
                 j = j + 1
                 id_arr(j) = id_arr_in(i)
                 k = 1
                 HG_arr_storb_in(j,k,1:4) = HG_arr_in(i,1:4)
              ELSE
                 HG_arr_storb_in(j,k,1:4) = HG_arr_in(i,1:4)
              END IF
              IF (k > k_max) THEN
                 k_max = k
              END IF
           END DO
           ALLOCATE(id_arr_storb_in(nobj))
           id_arr_storb_in = id_arr
           HG_arr_storb_in => reallocate(HG_arr_storb_in, nobj, k_max, 4)
           DEALLOCATE(id_arr, id_arr_in, HG_arr_in)
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
!!$        DO i=1,SIZE(id_arr_in)
!!$           DO WHILE (LEN_TRIM(id_arr_in(i)) < 7)
!!$              id_arr_in(i) = "0" // TRIM(id_arr_in(i))
!!$           END DO
!!$        END DO

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
           IF (element_type_in /= "cartesian" .AND. &
                element_type_comp_prm == "cartesian") THEN
              CALL toCartesian(storb_arr_in(i), "ecliptic")
           ELSE IF (element_type_in /= "cometary" .AND. &
                element_type_comp_prm == "cometary") THEN
              CALL toCometary(storb_arr_in(i))
           ELSE IF (element_type_in /= "keplerian" .AND. &
                element_type_comp_prm == "keplerian") THEN
              CALL toKeplerian(storb_arr_in(i))
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
        DEALLOCATE(orb_arr_in, stat=err)
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
     DEALLOCATE(pdf_arr_in, rchi2_arr_in, jac_arr_in, &
          repetition_arr_in, reg_apr_arr_in, &
          element_type_pdf_arr_in, stat=err)

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

  ! OpenOrb orbit files written use mjd instead of cal date:
  mjd_epoch = get_cl_option("--oorb-mjd", mjd_epoch)

  center = get_cl_option("--center=",11)
  IF (center < 0 .OR. center > 13) THEN
     CALL errorMessage("oorb", &
          "New center must be given with the --center=CENTER " // &
          "option where 1 <= CENTER <= 13. Default is the Sun (11).", 1)
     STOP
  END IF

  frame = get_cl_option("--frame=", "ecliptic")

  ! Initialize random number generator using current working
  ! directory and system clock
  IF (.NOT.get_cl_option("--fixran",.FALSE.)) THEN
     CALL GETCWD(str1)
     ALLOCATE(strlen1(1:LEN_TRIM(str1)))
     DO i=1,LEN_TRIM(str1)
        strlen1(i) = str1(i:i)
     END DO
     CALL SYSTEM_CLOCK(i)
     CALL initializeRandomNumberGenerator(-1*(i+SUM(ICHAR(strlen1))))
     IF (info_verb >= 2) THEN
        WRITE(stdout,"(A,I0)") "Random number generator initialized with ", &
             -1*(i+SUM(ICHAR(strlen1)))
     END IF
  END IF

  SELECT CASE (task)

  CASE ("none")

     CONTINUE

  CASE ("2oorbpdf", "toorbpdf")

     first = .TRUE.
     tmp_fname = get_cl_option("--id-in=","")
     IF (LEN_TRIM(tmp_fname) /= 0) THEN
        ! Read and sort ids
        CALL NEW(tmp_file, TRIM(tmp_fname))
        CALL setStatusOld(tmp_file)
        CALL OPEN(tmp_file)
        IF (error) THEN
           STOP
        END IF
        j = getNrOfLines(tmp_file)
        ALLOCATE(id_arr(j))
        DO i=1,j
           READ(getUnit(tmp_file),*) id_arr(i)
        END DO
        CALL quickSort(id_arr,errstr)
     END IF

     ! Convert orbits + uncertainties
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (get_cl_option("--id-in=",.FALSE.)) THEN
              j = binarySearch(id_arr_storb_in(i),id_arr,errstr)
              IF (j < 0 .OR. LEN_TRIM(errstr) /= 0) THEN
                 CYCLE
              END IF
           END IF
           IF (containsDiscretePDF(storb_arr_in(i))) THEN
              orb_arr_in => getSampleOrbits(storb_arr_in(i))
              pdf_arr_in => getDiscretePDF(storb_arr_in(i), element_type_out_prm)
              !pdf_arr_in = pdf_arr_in/sum(pdf_arr_in)

              DO j=1,SIZE(orb_arr_in)
                 CALL writeOpenOrbOrbitFile(lu_orb_out, first, &
                      element_type_out_prm, id_arr_storb_in(i), &
                      orb_arr_in(j), pdf=pdf_arr_in(j), &
                      element_type_pdf=element_type_out_prm, &
                      H=HG_arr_storb_in(i,j,1), &
                      G=HG_arr_storb_in(i,j,3), &
                      mjd=mjd_epoch, frame=frame)
                 first = .FALSE.
                 CALL NULLIFY(orb_arr_in(j))
              END DO
              DEALLOCATE(orb_arr_in, pdf_arr_in)
           ELSE
              orb = getNominalOrbit(storb_arr_in(i))
              cov = getCovarianceMatrix(storb_arr_in(i), &
                   element_type_out_prm, frame=frame)
              CALL writeOpenOrbOrbitFile(lu_orb_out, first, &
                   element_type_out_prm, id_arr_storb_in(i), &
                   orb=orb, cov=cov, H=HG_arr_storb_in(i,j,1), &
                   G=HG_arr_storb_in(i,j,3), mjd=mjd_epoch, &
                   frame=frame)
              first = .FALSE.
           END IF
        END DO
     ELSE
        CALL errorMessage("oorb / 2oorbpdf", &
             "Orbital-element PDFs were not detected in the input.", 1)
     END IF

     IF (ASSOCIATED(id_arr)) THEN
        DEALLOCATE(id_arr, stat=err)
     END IF

  CASE ("2oorb", "toorb")

     norb = HUGE(norb)
     norb = get_cl_option("--norb=",norb)
     H_max = get_cl_option("--H-max=",HUGE(H_max))
     first = .TRUE.

     tmp_fname = get_cl_option("--id-in=","")
     IF (LEN_TRIM(tmp_fname) /= 0) THEN
        ! Read and sort ids
        CALL NEW(tmp_file, TRIM(tmp_fname))
        CALL setStatusOld(tmp_file)
        CALL OPEN(tmp_file)
        IF (error) THEN
           STOP
        END IF
        j = getNrOfLines(tmp_file)
        ALLOCATE(id_arr(j))
        DO i=1,j
           READ(getUnit(tmp_file),*) id_arr(i)
        END DO
        CALL quickSort(id_arr,errstr)
     END IF

     ! Convert orbits
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)        
           IF (get_cl_option("--id-in=",.FALSE.)) THEN
              j = binarySearch(id_arr_storb_in(i),id_arr,errstr)
              IF (j < 0 .OR. LEN_TRIM(errstr) /= 0) THEN
                 CYCLE
              END IF
           END IF
           IF (containsDiscretePDF(storb_arr_in(i))) THEN
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
                    IF (HG_arr_storb_in(i,j,1) > H_max) THEN
                       ! too large H magnitude -> skip it
                       CYCLE
                    END IF
                 END IF
                 k = k + 1
                 CALL writeOpenOrbOrbitFile(lu_orb_out, first, &
                      element_type_out_prm, id_arr_storb_in(i), &
                      orb_arr_in(j), H=HG_arr_storb_in(i,j,1), &
                      G=HG_arr_storb_in(i,j,3), mjd=mjd_epoch, &
                      frame=frame)
                 first = .FALSE.
                 IF (k == norb) THEN
                    EXIT
                 END IF
              END DO
           ELSE
              orb = getNominalOrbit(storb_arr_in(i))
              CALL writeOpenOrbOrbitFile(lu_orb_out, first, &
                   element_type_out_prm, id_arr_storb_in(i), &
                   orb=orb, H=HG_arr_storb_in(i,1,1), &
                   G=HG_arr_storb_in(i,1,3), mjd=mjd_epoch, &
                   frame=frame)
              first = .FALSE.
           END IF
        END DO
     ELSE IF (ASSOCIATED(orb_arr_in)) THEN
        j = 0
        DO i=1,SIZE(orb_arr_in)
           IF (get_cl_option("--id-in=",.FALSE.)) THEN
              k = binarySearch(id_arr_in(i),id_arr,errstr)
              IF (k < 0 .OR. LEN_TRIM(errstr) /= 0) THEN
                 CYCLE
              END IF
           END IF
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
           CALL writeOpenOrbOrbitFile(lu_orb_out, first, &
                element_type_out_prm, id_arr_in(i), orb_arr_in(i), &
                H=HG_arr_in(i,1), G=HG_arr_in(i,3), mjd=mjd_epoch, &
                frame=frame)
           first = .FALSE.
           IF (j == norb) THEN
              EXIT
           END IF
        END DO
     END IF
     IF (ASSOCIATED(id_arr)) THEN
        DEALLOCATE(id_arr, stat=err)
     END IF

  CASE ("orb2des", "orbitstodes")

     H_max = get_cl_option("--H-max=",HUGE(H_max))
     i_min = get_cl_option("--i-min=",0.0_bp)
     i_max = get_cl_option("--i-max=",HUGE(i_max))
     first = .TRUE.

     tmp_fname = get_cl_option("--id-in=","")
     IF (LEN_TRIM(tmp_fname) /= 0) THEN
        ! Read and sort ids
        CALL NEW(tmp_file, TRIM(tmp_fname))
        CALL setStatusOld(tmp_file)
        CALL OPEN(tmp_file)
        IF (error) THEN
           STOP
        END IF
        j = getNrOfLines(tmp_file)
        ALLOCATE(id_arr(j))
        DO i=1,j
           READ(getUnit(tmp_file),*) id_arr(i)
        END DO
        CALL quickSort(id_arr,errstr)
     END IF

     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           IF (get_cl_option("--id-in=",.FALSE.)) THEN
              j = binarySearch(id_arr_storb_in(i),id_arr,errstr)
              IF (j < 0 .OR. LEN_TRIM(errstr) /= 0) THEN
                 CYCLE
              END IF
           END IF
           IF (containsDiscretePDF(storb_arr_in(i))) THEN
              orb_arr_in => getSampleOrbits(storb_arr_in(i))
              DO j=1,SIZE(orb_arr_in)
                 IF (get_cl_option("--H-max=",.FALSE.)) THEN
                    IF (HG_arr_storb_in(i,j,1) > H_max) THEN
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
                 CALL writeDESOrbitFile(lu_orb_out, first, &
                      element_type_out_prm, id_arr_storb_in(i), &
                      orb_arr_in(j), HG_arr_storb_in(i,j,1), &
                      frame=frame, center=center)
                 IF (error) THEN
                    CALL errorMessage("oorb", &
                         "DES output failed at orbit:", 1)
                    WRITE(stderr,*) j
                    STOP
                 END IF
                 first = .FALSE.
                 CALL NULLIFY(orb_arr_in(j))
              END DO
              DEALLOCATE(orb_arr_in)
           ELSE
              IF (HG_arr_storb_in(i,1,1) > H_max) THEN
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
              CALL writeDESOrbitFile(lu_orb_out, first, &
                   element_type_out_prm, id_arr_storb_in(i), orb, &
                   HG_arr_storb_in(i,1,1), frame=frame, &
                   center=center)
              first = .FALSE.
              CALL NULLIFY(orb)
           END IF
        END DO
        STOP
     END IF
     IF (ASSOCIATED(orb_arr_in)) THEN
        DO i=1,SIZE(orb_arr_in)
           IF (get_cl_option("--id-in=",.FALSE.)) THEN
              j = binarySearch(id_arr_in(i),id_arr,errstr)
              IF (j < 0 .OR. LEN_TRIM(errstr) /= 0) THEN
                 CYCLE
              END IF
           END IF
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
           CALL writeDESOrbitFile(lu_orb_out, first, element_type_out_prm, &
                id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), frame=frame, &
                center=center)
           IF (error) THEN
              CALL errorMessage("oorb", &
                   "DES output failed at orbit:", 1)
              WRITE(stderr,*) i
              STOP
           END IF
           first = .FALSE.
        END DO
     END IF
     IF (ASSOCIATED(id_arr)) THEN
        DEALLOCATE(id_arr, stat=err)
     END IF

  CASE ("obs2des", "observationstodes")

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

  CASE ("2mpc3", "tompc3")

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

  CASE ("2mpc", "tompc")

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

  CASE ("astorb")

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
     iorb = 0
     ALLOCATE(str_arr(6), real_arr(5), int_arr(23))
     str_arr(1)(1:LEN(str_arr(1))) = " "
     str_arr(2)(1:LEN(str_arr(2))) = " "
     DO
        line(:) = " "
        str_arr(1)(1:6) = " "
        str_arr(2)(1:18) = " "
        iorb = iorb + 1
        READ(getUnit(orb_in_file),"(A)",iostat=err) line
        IF (err > 0) THEN
           CALL errorMessage("oorb / astorbtoorb", &
                'Error while reading line:',1)
           WRITE(stderr,"(I0)") iorb
           STOP
        ELSE IF (err < 0) THEN
           ! Input file read.
           EXIT
        END IF
        READ(line,TRIM(frmt),iostat=err) &
             str_arr(1)(1:6), str_arr(2)(1:18), str_arr(3)(1:15), &
             H_value, G_value, str_arr(4)(1:4), str_arr(5)(1:5), &
             str_arr(6)(1:4), int_arr(1:8), year, month, iday, &
             elements(6:1:-1), int_arr(9:11), real_arr(1:2), &
             int_arr(12:14), real_arr(3), int_arr(15:17), &
             real_arr(4), int_arr(18:20), real_arr(5), int_arr(21:23)
        IF (err > 0) THEN
           CALL errorMessage("oorb / astorbtoorb", &
                'Error while converting line:',1)
           WRITE(stderr,"(I0)") iorb
           STOP
        END IF
        READ(line(43:47),*,iostat=err) H_value
        IF (err > 0) THEN
           CALL errorMessage("oorb / astorbtoorb", &
                'Error while converting H str to real on line:',1)
           WRITE(stderr,"(I0)") iorb
           STOP
        END IF
!!$        ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!!$        ! Discard some orbits depending on the following requirement:
!!$        IF (int_arr(7) < 30 .OR. & ! orbital arc in days
!!$                                !elements(1)*(1.0_bp-elements(2)) < 1.666_bp .OR. &
!!$                                !elements(1) <= 2.0_bp .or. elements(1) >= 3.5_bp) then ! .or. &
!!$             elements(1)*(1.0_bp-elements(2)) < 1.3_bp .OR. &
!!$             elements(1) > 4.5_bp .OR. &
!!$                                !elements(3) > 60.0_bp .OR. &
!!$             H_value > 15.5_bp) then
!!$             int_arr(9) >= 2008) THEN
!!$           CYCLE
!!$        END IF
!!$        ! End requirement
!!$        ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        IF (LEN_TRIM(str_arr(2)) /= 0) THEN
           !! 9-char alphanumeric encoded designation
           !CALL encodeMPC3Designation(str_arr(2))
           !! 7-char alphanumeric encoded designation (current MPC format)
           CALL encodeMPCDesignation(str_arr(2))
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
           !! 7-char numeric number
           !DO WHILE (LEN_TRIM(id) < 7)
           !   id = '0' // TRIM(id)
           !END DO
           !! 5-char alphanumeric number (current MPC format)
           IF (LEN_TRIM(id) == 6) THEN
              CALL toInt(id(1:2), i, error)
              id(2:5) = id(3:6)
              id(6:6) = " "
              id(1:1) = mpc_conv_table(i)
           ELSE IF (LEN_TRIM(id) > 6) THEN
              CALL errorMessage("oorb / astorbtoorb", &
                   "Number (" // TRIM(id) // ") too large -> cannot encode.",1)
              STOP             
           END IF
           DO WHILE (LEN_TRIM(id) < 5)
              id = '0' // TRIM(id)
           END DO
        ELSE
           id = TRIM(ADJUSTL(str_arr(2)))
        END IF
        SELECT CASE (TRIM(orbit_format_out))
        CASE ("des")
           CALL writeDESOrbitFile(lu_orb_out, i==1, &
                element_type_out_prm, id, orb, H_value, 1, 6, &
                REAL(int_arr(7),bp), "OPENORB", frame=frame, &
                center=center)
        CASE ("orb")
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id), orb=orb, H=H_value, G=G_value, &
                mjd=mjd_epoch, frame=frame)
        CASE default
           CALL errorMessage("oorb / astorbtoorb", &
                "Orbit format " // TRIM(orbit_format_out) // &
                " not supported.",1)
           STOP           
        END SELECT
        CALL NULLIFY(orb)
     END DO
     CALL NULLIFY(orb_in_file)

  CASE ("mpcorb")

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
!!$        elements = getElements(orb_arr_in(i), "keplerian")
!!$        ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!!$        ! Discard some orbits depending on the following requirement:
!!$        IF (elements(1) < 2.1_bp .or. elements(1) > 3.3_bp .or. & 
!!$             elements(2) > 0.35_bp .or. elements(3) > 26.0_bp*rad_deg) THEN ! not-MBO
!!$        IF (.not.(elements(1)*(1.0_bp-elements(2)) <= 1.3_bp .and. &
!!$             elements(1)*(1.0_bp-elements(2)) > 1.017_bp)) THEN ! not-Amor
!!$        IF (.not.(elements(1)*(1.0_bp-elements(2)) <= 1.017_bp .and. &
!!$             elements(1) > 1.0_bp)) THEN ! not-Apollo
!!$        IF (.not.(elements(1)*(1.0_bp+elements(2)) > 0.983_bp .and. &
!!$             elements(1) <= 1.0_bp)) THEN ! not-Aten
!!$           CYCLE
!!$        END IF
!!$        ! End requirement
!!$        ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        SELECT CASE (TRIM(orbit_format_out))
        CASE ("des")
           CALL writeDESOrbitFile(lu_orb_out, i==1, element_type_out_prm, &
                id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), 1, 6, &
                arc_arr(i), "OPENORB", frame=frame, center=center)
        CASE ("orb")
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
                H=HG_arr_in(i,1), G=HG_arr_in(i,2), &
                mjd=mjd_epoch, frame=frame)
        CASE default
           CALL errorMessage("oorb / mpcorbtoorb", &
                "Orbit format " // TRIM(orbit_format_out) // &
                " not supported.",1)
           STOP           
        END SELECT
        IF (error) THEN
           CALL errorMessage("oorb / mpcorbtoorb", &
                "TRACE BACK (15)", 1)
           STOP        
        END IF
     END DO

  CASE ("2planetocentric")

     IF (ALLOCATED(storb_arr_in)) THEN
        IF (SIZE(storb_arr_in) == 1) THEN
           orb_arr_in => getSampleOrbits(storb_arr_in(1))
           ALLOCATE(id_arr_in(SIZE(orb_arr_in)), HG_arr_in(SIZE(orb_arr_in),4))
           id_arr_in = id_arr_storb_in(1)
           HG_arr_in = HG_arr_storb_in(1,:,:)
        ELSE
           CALL errorMessage("oorb / 2planetocentric", &
                "Sampled pdf can only be used for one (1) object at a time.", 1)
           STOP
        END IF
     END IF

     DO i=1,SIZE(orb_arr_in)
        t = getTime(orb_arr_in(i))
        mjd = getMJD(t, "TT")
        ! Get Sun's coordinates at the epoch as seen from the planet:
        planeph => JPL_ephemeris(mjd, 11, center, error)
        CALL NEW(ccoord, planeph(1,1:6), "equatorial", t)
        DEALLOCATE(planeph)
        CALL rotateToEcliptic(ccoord)
        coordinates = getCoordinates(ccoord)
        CALL NULLIFY(ccoord)
        elements = getElements(orb_arr_in(i), "cartesian", "ecliptic")
        ! New center -> asteroid = new_center -> Sun + Sun -> asteroid
        CALL NEW(orb, coordinates + elements, "cartesian", "ecliptic", t, center=center)
        SELECT CASE (TRIM(orbit_format_out))
        CASE ("des")
           CALL writeDESOrbitFile(lu_orb_out, i==1, element_type_out_prm, &
                id_arr_in(i), orb, HG_arr_in(i,1), 1, 6, &
                -1.0_bp, "OPENORB", frame=frame, center=center)
        CASE ("orb")
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id_arr_in(i)), orb=orb, &
                H=HG_arr_in(i,1), G=HG_arr_in(i,2), &
                mjd=mjd_epoch, frame=frame)
        CASE default
           CALL errorMessage("oorb / planetocentricorbits", &
                "Orbit format " // TRIM(orbit_format_out) // &
                " not supported.",1)
           STOP           
        END SELECT
        IF (error) THEN
           CALL errorMessage("oorb / planetocentricorbits", &
                "TRACE BACK (15)", 1)
           STOP        
        END IF
        CALL NULLIFY(orb)        
        CALL NULLIFY(t)
     END DO


  CASE ("obs2ecl")

     ! Convert observed RA,Dec to topocentric ecliptic coordinates.

     DO i=1,getNrOfObservations(obss_in)
        obs = getObservation(obss_in,i)
        ! Get number or designation
        id = getID(obs)
        ! Get observation date
        t = getTime(obs)
        mjd_utc = getMJD(t,"UTC")
        CALL NULLIFY(t)
        ! Get apparent magnitude
        mag = getMagnitude(obs)
        ! Compute topocentric opposition coordinates starting from the
        ! heliocentric observatory coordinates.
        ccoord = getObservatoryCCoord(obs)
        scoord = getSCoord(ccoord)
        CALL NULLIFY(ccoord)
        CALL rotateToEcliptic(scoord)
        pos_opp = getPosition(scoord)
        CALL NULLIFY(scoord)
        scoord = getObservationSCoord(obs)
        CALL rotateToEcliptic(scoord)
        pos_ast = getPosition(scoord)
        CALL NULLIFY(scoord)
        ! Compute opposition-centered topocentric ecliptic coordinates
        toclon = pos_ast(2) - pos_opp(2)
        toclat = pos_ast(3) - pos_opp(3)
        IF (toclon > pi) THEN
           toclon = toclon - two_pi
        ELSE IF (toclon < -pi) THEN
           toclon = toclon + two_pi
        END IF
        WRITE(stdout,"(A,1X,I0,1X,F7.5,2(1X,F12.7),1X,F6.2)") &
             TRIM(id), FLOOR(mjd_utc), mjd_utc-FLOOR(mjd_utc), &
             toclon/rad_deg, toclat/rad_deg, mag
        CALL NULLIFY(obs)
     END DO





  CASE ("ranging")

     ! MCMC

     ! Orbital inversion using various flavors of statistical
     ! [orbital] ranging, that is, sampling in spherical topocentric
     ! coordinates without making any assumptions on the shape of the
     ! resulting orbital-element pdf.

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
     generat_multiplier = -1.0_bp
     sor_genwin_offset = -1.0_bp
     sor_iterate_bounds = .TRUE.
     accwin_multiplier = -1.0_bp
     gaussian_rho = .FALSE.
     regularized = .FALSE.
     write_residuals = .FALSE.
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          asteroid_perturbers=asteroid_perturbers, &
          integrator=integrator, &
          integration_step=integration_step, &
          dyn_model_init=dyn_model_init, &
          integrator_init=integrator_init, &
          integration_step_init=integration_step_init, &
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
          generat_multiplier=generat_multiplier, &
          sor_genwin_offset=sor_genwin_offset, &
          sor_iterate_bounds=sor_iterate_bounds, &
          accwin_multiplier=accwin_multiplier, &
          regularized_pdf=regularized, &
          sor_random_obs=random_obs, sor_rho_gauss=gaussian_rho, &
          write_residuals=write_residuals, &
          ls_correction_factor=ls_correction_factor, &
          ls_element_mask=ls_element_mask, &
          ls_niter_major_max=ls_niter_major_max, &
          ls_niter_major_min=ls_niter_major_min, &
          ls_niter_minor=ls_niter_minor)
     IF (error) THEN
        CALL errorMessage("oorb / ranging", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        CALL NULLIFY(epoch)
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     ! Print header before printing first orbit:
     first = .TRUE.

     DO i=1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (15)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")."
        END IF
        nobs = getNrOfObservations(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (20)", 1)
           STOP
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(4X,A,1X,F8.2,1X,A)") &
                "Observational arc: ", dt, "days"
        END IF

        ! Initialize the rho2-rho1 range based on the observational timespan.
        IF (sor_rho_init(3) > HUGE(sor_rho_init(3))/2) THEN
           IF (dt > 10) THEN
              sor_rho_init(3) = -0.5_bp
           ELSE
              sor_rho_init(3) = -0.05_bp*dt
           END IF
           sor_rho_init(4) = -sor_rho_init(3)
        END IF

        ! Initialize StochasticOrbit
        IF (ASSOCIATED(id_arr_storb_in) .AND. ALLOCATED(storb_arr_in)) THEN
           ! ...with an estimate for the orbit + observations:
           DO j=1,SIZE(id_arr_storb_in)
              IF (TRIM(id) == TRIM(id_arr_storb_in(j))) THEN
                 EXIT
              END IF
           END DO
           IF (j <= SIZE(id_arr_storb_in)) THEN
              storb = copy(storb_arr_in(j))
           END IF
        ELSE IF (ASSOCIATED(id_arr_in) .AND. ASSOCIATED(orb_arr_in)) THEN
           ! ...with a sample of orbits + observations:
           norb = 0
           DO j=1,SIZE(id_arr_in)
              IF (TRIM(id) == TRIM(id_arr_in(j))) THEN
                 norb = norb + 1
                 orb_arr => reallocate(orb_arr, norb)
                 orb_arr(norb) = copy(orb_arr_in(j))
              END IF
           END DO
           ALLOCATE(pdf_arr(norb), stat=err)
           pdf_arr = 1.0_bp
           pdf_arr = pdf_arr/SUM(pdf_arr)
           CALL NEW(storb, orb_arr=orb_arr, pdf_arr=pdf_arr, &
                element_type=getElementType(orb_arr(1)), &
                obss=obss_sep(i))
        ELSE
           ! ...with just observations:
           CALL NEW(storb, obss_sep(i))        
        END IF
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (65)", 1)
           STOP
        END IF

        ! Determine inversion epoch
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (55)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           IF (get_cl_option("--exact-mid-epoch", .FALSE.)) THEN
              mjd = mjd+dt/2.0_bp
           ELSE
              mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           END IF
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (60)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF

        ! Set parameters for inversion
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             dchi2_rejection = dchi2_rejection, &
             dchi2_max = dchi2_max, &
             regularized_pdf = regularized, &
             accept_multiplier=accwin_multiplier, &
             apriori_a_max=apriori_a_max, apriori_a_min=apriori_a_min, &
             apriori_periapsis_max=apriori_periapsis_max, &
             apriori_periapsis_min=apriori_periapsis_min, &
             apriori_apoapsis_max=apriori_apoapsis_max, &
             apriori_apoapsis_min=apriori_apoapsis_min, &
             apriori_rho_min=apriori_rho_min, &
             sor_2point_method=sor_2point_method, &
             sor_2point_method_sw=sor_2point_method_sw, &
             sor_norb_sw=sor_norb_sw, sor_ntrial_sw=sor_ntrial_sw, &
             sor_norb=sor_norb, sor_ntrial=sor_ntrial, &
             sor_rho1_l=sor_rho_init(1), sor_rho1_u=sor_rho_init(2), &
             sor_rho2_l=sor_rho_init(3), sor_rho2_u=sor_rho_init(4), &
             sor_niter=sor_niter, sor_iterate_bounds=sor_iterate_bounds, &
             sor_random_obs_selection=.FALSE., &
             gaussian_pdf=gaussian_rho, &
             generat_multiplier=generat_multiplier, &
             sor_generat_offset=sor_genwin_offset)
        IF (error) THEN
           CALL errorMessage("oorb / ranging", &
                "TRACE BACK (70)", 1)
           STOP
        END IF

        ! These are here to provide partial backwards compatibility
        IF (sor_type_prm == 2) THEN
           flavor = "mc"
        ELSE IF (sor_type_prm == 3) THEN
           flavor = "stepwise-mc"
        ELSE
           flavor = "mc"
        END IF

        flavor = get_cl_option("--flavor=",flavor)

        SELECT CASE (TRIM(flavor))

        CASE ("mc")

           CALL setParameters(storb, &
                sor_niter=sor_niter)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (85)", 1)
              STOP
           END IF
           CALL autoStatisticalRanging(storb)

        CASE ("stepwise-mc")

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

        CASE ("random-walk")

           ! Orbital inversion using random-walk ranging (Muinonen et
           ! al. 2015), that is, unlike in MCMC ranging, the
           ! acceptance value is based not on the quotient of proposal
           ! p.d.f.:s, but on the comparison of difference of proposal
           ! and previous chi^2 values to the predefined boundary
           ! value, usually 20.1 (corresponding to 3*sigma).

           CALL autoRandomWalkRanging(storb)

        CASE ("mcmc")

           CALL autoMCMCRanging(storb)

        CASE default

           CALL errorMessage("oorb / ranging", &
                "Unknown flavor of ranging:", 1)
           IF (err_verb >= 1) THEN
              WRITE(stderr,"(2X,A)") TRIM(flavor)
           END IF
           STOP              

        END SELECT

        IF (error) THEN

           CALL errorMessage("oorb / ranging", &
                "Ranging failed:", 1)
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
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (90)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (95)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (100)", 1)
              STOP
           END IF

        ELSE

           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Ranging for object", &
                   TRIM(id), "is ready."
           END IF

           IF (pp_H_estimation) THEN
              IF (info_verb >= 2) THEN
                 WRITE(stdout,"((1X,A))") "Start estimation of H and G."
              END IF
              CALL NEW(physparam, storb)
              IF (pp_G > 99.0_bp) THEN
                 CALL estimateHAndG(physparam, obss_sep(i))
              ELSE
                 CALL estimateHAndG(physparam, obss_sep(i), &
                      input_G=pp_G, input_delta_G=pp_G_unc)
              END IF
              HG_arr_in => getH0Distribution(physparam)
              IF (.NOT.error) THEN
                 ALLOCATE(temp_arr(SIZE(HG_arr_in,dim=1),4))
                 temp_arr(:,1:2) = HG_arr_in(:,1:2)
                 DEALLOCATE(HG_arr_in)
                 HG_arr_in => getGDistribution(physparam)
                 temp_arr(:,3:4) = HG_arr_in(:,1:2)
                 DEALLOCATE(HG_arr_in)
                 ALLOCATE(HG_arr_in(SIZE(temp_arr,dim=1),4))
                 HG_arr_in(:,:) = temp_arr
                 DEALLOCATE(temp_arr)
              ELSE
                 error = .FALSE.
                 NULLIFY(HG_arr_in)
              END IF
              CALL NULLIFY(physparam)
              IF (info_verb >= 2) THEN
                 WRITE(stdout,"((1X,A))") "End estimation of H and G."
              END IF
           END IF

           CALL toString(dt, str, error, frmt="(F10.2)")
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (165)", 1)
              STOP
           END IF
           str = TRIM(id) // "_"// TRIM(str)
           ! WRITE RANGING OUTPUT FILE WITH
           CALL NEW(out_file, TRIM(id) // ".sor")
           !           CALL NEW(out_file, TRIM(str) // ".sor")
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (105)", 1)
              STOP
           END IF
           ! 1) OBSERVATIONS AND
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (110)", 1)
              STOP
           END IF
           WRITE(getUnit(out_file),"(A)") "#" 
           ! 2) RANGING PARAMETERS.
           CALL writeSORResults(storb, obss_sep(i), getUnit(out_file))
           IF (error) THEN
              CALL errorMessage("oorb / ranging", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)

           ! WRITE ORBITAL-ELEMENT PDF
           orb_arr_cmp => getSampleOrbits(storb)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (125)", 1)
              STOP
           END IF
           pdf_arr_cmp => getDiscretePDF(storb, element_type_comp_prm)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (130)", 1)
              STOP
           END IF
           rchi2_arr_cmp => getReducedChi2Distribution(storb)
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (135)", 1)
              STOP
           END IF
           SELECT CASE (TRIM(flavor))
           CASE ("mc", "stepwise-mc")
              CALL getResults(storb, &
                   reg_apr_arr=reg_apr_arr_cmp, &
                   jac_arr=jac_arr_cmp)
           CASE ("mcmc", "random-walk")
              CALL getResults(storb, &
                   repetition_arr_cmp=repetition_arr_cmp)
           END SELECT
           IF (error) THEN
              CALL errorMessage("oorb / ranging ", &
                   "TRACE BACK (140)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NEW(out_file, TRIM(id) // ".orb")
              !              CALL NEW(out_file, TRIM(str) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (145)", 1)
                 STOP
              END IF
              lu_orb_out = getUnit(out_file)
           END IF
           DO j=1,SIZE(orb_arr_cmp,dim=1)
              IF (orbit_format_out == "orb") THEN
                 IF (ASSOCIATED(HG_arr_in)) THEN
                    SELECT CASE (TRIM(flavor))
                    CASE ("mc", "stepwise-mc")
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
                            G=HG_arr_in(j,3), &
                            mjd=mjd_epoch)
                    CASE ("mcmc", "random-walk")
                       CALL writeOpenOrbOrbitFile(lu_orb_out, &
                            print_header=j==1, &
                            element_type_out=element_type_out_prm, &
                            id=id, &
                            orb=orb_arr_cmp(j), &
                            element_type_pdf=element_type_comp_prm, &
                            pdf=pdf_arr_cmp(j), &
                            rchi2=rchi2_arr_cmp(j), &
                            repetitions=repetition_arr_cmp(j), &
                            H=HG_arr_in(j,1), &
                            G=HG_arr_in(j,3), &
                            mjd=mjd_epoch)
                    END SELECT
                 ELSE
                    SELECT CASE (TRIM(flavor))
                    CASE ("mc", "stepwise-mc")
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
                            mjd=mjd_epoch)
                    CASE ("mcmc", "random-walk")
                       CALL writeOpenOrbOrbitFile(lu_orb_out, &
                            print_header=j==1, &
                            element_type_out=element_type_out_prm, &
                            id=id, &
                            orb=orb_arr_cmp(j), &
                            element_type_pdf=element_type_comp_prm, &
                            pdf=pdf_arr_cmp(j), &
                            rchi2=rchi2_arr_cmp(j), &
                            repetitions=repetition_arr_cmp(j), &
                            mjd=mjd_epoch)
                    END SELECT
                 END IF
              ELSE IF (orbit_format_out == "des") THEN
                 CALL errorMessage("oorb / ranging ", &
                      "DES format not yet supported for ranging output.", 1)
                 STOP                 
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (150)", 1)
                 STOP
              END IF
           END DO
           IF (separately) THEN
              CALL NULLIFY(out_file)
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(id) // ".orb")
              END IF
           END IF

           ! WRITE RESIDUALS FILE
           IF (write_residuals) THEN
              CALL NEW(out_file, TRIM(out_fname) // ".res")
              CALL setPositionAppend(out_file)
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (155)", 1)
                 STOP
              END IF
              CALL writeResiduals(storb, obss_sep(i), getUnit(out_file))
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (160)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
           END IF

           ! MAKE PLOTS
           IF (plot_results) THEN
              CALL toString(dt, str, error, frmt="(F10.2)")
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (165)", 1)
                 STOP
              END IF
              str = TRIM(id) // "_"// TRIM(str)
              CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
                   "_ranging_residual_stamps.eps")
              !CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
              !     "_" // TRIM(flavor) //"_residual_stamps.eps")
              IF (error) THEN
                 CALL errorMessage("oorb / ranging", &
                      "TRACE BACK (170)", 1)
                 STOP
              END IF
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_ranging_residual_stamps.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_ranging_residual_stamps.eps* &")
              END IF
              ALLOCATE(elements_arr(SIZE(orb_arr_cmp,dim=1),7), stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / ranging", &
                      "Could not allocate memory (3).", 1)
                 STOP
              END IF

              ! Topocentric ranges
              CALL NEW(tmp_file, "sor_histo.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (175)", 1)
                 STOP
              END IF
              CALL getResults(storb, sor_rho_arr_cmp=sor_rho_arr)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (201)", 1)
                 STOP
              END IF
              ALLOCATE(histo1(50,2), histo2(50,2), stat=err)
              CALL histogram(sor_rho_arr(:,1), histo1)
              CALL histogram(sor_rho_arr(:,1), histo2, pdf=pdf_arr_cmp)
              DO j=1,50
                 WRITE(getUnit(tmp_file),*) histo1(j,1:2), histo2(j,2)
              END DO
              DEALLOCATE(histo1,histo2, stat=err)
              CALL NULLIFY(tmp_file)

              CALL NEW(tmp_file, TRIM(str)// "_ranging_orbits.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (175)", 1)
                 STOP
              END IF
              CALL NEW(tmp_file2, TRIM(str)// "_ranging_ranges.out")
              CALL OPEN(tmp_file2)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (176)", 1)
                 STOP
              END IF
              DO j=1,SIZE(orb_arr_cmp,dim=1)
                 IF (element_type_comp_prm == "cartesian") THEN
                    CALL rotateToEcliptic(orb_arr_cmp(j))
                 END IF
                 elements_arr(j,1:6) = getElements(orb_arr_cmp(j), element_type_comp_prm)
                 IF (error) THEN
                    CALL errorMessage("oorb / ranging", &
                         "TRACE BACK (180)", 1)
                    STOP
                 END IF
                 IF (element_type_comp_prm == "keplerian") THEN
                    elements_arr(j,3:6) = elements_arr(j,3:6)/rad_deg
                 END IF
                 elements_arr(j,7) = pdf_arr_cmp(j)
                 t = getTime(orb_arr_cmp(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / ranging ", &
                         "TRACE BACK (185)", 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmp_file),*) &
                      elements_arr(j,1:6), &
                      pdf_arr_cmp(j), &
                      getCalendarDateString(t,"tdt")
                 IF (error) THEN
                    CALL errorMessage("oorb / ranging ", &
                         "TRACE BACK (190)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(t)
                 WRITE(getUnit(tmp_file2),*) sor_rho_arr(j,1:2)
              END DO
              CALL NULLIFY(tmp_file)
              CALL NULLIFY(tmp_file2)
              CALL NEW(tmp_file, TRIM(str) // &
                   "_ranging_sample_standard_deviations.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (195)", 1)
                 STOP
              END IF
              WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                   advance="no") &
                   getObservationalTimespan(obss_sep(i))
              IF (error) THEN
                 CALL errorMessage("oorb / ranging ", &
                      "TRACE BACK (200)", 1)
                 STOP
              END IF
              DO j=1,6
                 SELECT CASE (TRIM(flavor))
                 CASE ("mc", "random-walk")
                    CALL moments(elements_arr(:,j), &
                         pdf=pdf_arr_cmp, std_dev=stdev, errstr=errstr)
                 CASE ("mcmc")
                    CALL moments(elements_arr(:,j), &
                         std_dev=stdev, errstr=errstr)
                 END SELECT
                 WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                      advance="no") stdev
                 stdev_arr(j)=stdev
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
                   "_ranging_orbits.out sor_orbits.out")
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_car.gp")
              ELSE
                 IF (stdev_arr(1) < 1.0_bp) THEN
                    CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_kep_nolog.gp")
                 ELSE
                    CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_kep.gp")
                 END IF
              END IF
              CALL system("cp sor_results.eps " // TRIM(str) // &
                   "_ranging_" // TRIM(element_type_comp_prm) // &
                   "_results.eps")
              CALL system("rm -f sor_orbits.out sor_results.eps " // & 
                   TRIM(str) // "_ranging_orbits.out ")
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_ranging_" // &
                      TRIM(element_type_comp_prm) // "_results.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_ranging_" // &
                      TRIM(element_type_comp_prm) // &
                      "_results.eps* &")
              END IF
              !              CALL system("cp " // TRIM(str) // &
              !                   "_ranging_ranges.out sor_ranges.out")
              CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_range.gp")
              CALL system("cp sor_ranges.eps " // TRIM(str) // &
                   "_ranging_ranges.eps")
              CALL system("rm -f sor_histo.out sor_ranges.eps")
           END IF

           ! DEALLOCATE MEMORY
           DO j=1,SIZE(orb_arr_cmp)
              CALL NULLIFY(orb_arr_cmp(j))
           END DO
           DEALLOCATE(orb_arr_cmp, stat=err)
           DEALLOCATE(pdf_arr_cmp, stat=err)
           DEALLOCATE(rchi2_arr_cmp, stat=err)
           DEALLOCATE(repetition_arr_cmp, stat=err)
           DEALLOCATE(reg_apr_arr_cmp, stat=err)
           DEALLOCATE(jac_arr_cmp, stat=err)
           DEALLOCATE(sor_rho_arr, stat=err)

           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Object", &
                   TRIM(id), "successfully processed."
           END IF
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb > 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
     END DO




  CASE ("simplex")

     ! Orbit optimization using the downhill simplex method.

     force_earth_impact_at_epoch = &
          get_cl_option("--force-earth-impact-at-epoch", .FALSE.)

     ALLOCATE(orb_arr_(7))

     CALL NULLIFY(epoch)
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model_init=dyn_model_init, &
          integrator_init=integrator_init, &
          integration_step_init=integration_step_init, &
          accwin_multiplier=accwin_multiplier, &
          smplx_niter=smplx_niter, &
          smplx_tol=smplx_tol, &
          smplx_force=smplx_force)
     IF (error) THEN
        CALL errorMessage("oorb / simplex", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        CALL NULLIFY(epoch)
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     DEALLOCATE(HG_arr_storb_in, stat=err)
     ALLOCATE(HG_arr_storb_in(SIZE(obss_sep),7,4))
     HG_arr_storb_in = 99.9_bp

     ! Print header before printing first orbit:
     first = .TRUE.
     DO i=1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "TRACE BACK (25)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")."
        END IF
        nobs = getNrOfObservations(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (nobs < 3) THEN
           CALL errorMessage("oorb / simplex", &
                "Too few observations:", 1)
           WRITE(stderr,*) "ID: ", TRIM(id), "  and number of observations: ", nobs
           CYCLE
        END IF
        norb = 0
        IF (ALLOCATED(storb_arr_in)) THEN
           DO j=1,SIZE(id_arr_storb_in,dim=1)
              IF (id_arr_storb_in(j) == id) THEN
                 IF (containsDiscretePDF(storb_arr_in(j))) THEN
                    orb_arr => getSampleOrbits(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / simplex", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = SIZE(orb_arr)
                    EXIT
                 ELSE
                    ALLOCATE(orb_arr(1))
                    orb_arr(1) = getNominalOrbit(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / simplex", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = 1
                    EXIT                    
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / simplex", &
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
                    CALL errorMessage("oorb / simplex", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / simplex", &
                   "Initial orbit not available.", 1)
              STOP
           ELSE
              orb_arr => reallocate(orb_arr,norb)
           END IF
        END IF
        IF (norb < 7) THEN
           CALL errorMessage("oorb / simplex", &
                "Not enough initial orbits available - at least 7 needed.", 1)
           STOP
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (70)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL NEW(storb, obss_sep(i))        
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "TRACE BACK (75)", 1)
           STOP
        END IF
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             accept_multiplier=accwin_multiplier, &
             smplx_niter=smplx_niter, &
             smplx_tol=smplx_tol, &
             smplx_force=smplx_force)
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "TRACE BACK (80)", 1)
           STOP
        END IF
        DO j=1,SIZE(orb_arr,dim=1),7
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(A,1X,I0,1X))") "Orbits", j, "to", &
                   j+6, "out of", SIZE(orb_arr,dim=1)
           END IF
           DO k=1,7
              CALL NULLIFY(orb_arr_(k))
              orb_arr_(k) = copy(orb_arr(j+k-1))
              CALL setParameters(orb_arr_(k), &
                   dyn_model=dyn_model_init, &
                   perturbers=perturbers, &
                   asteroid_perturbers=asteroid_perturbers, &
                   integrator=integrator_init, &
                   integration_step=integration_step_init)
              IF (error) THEN
                 CALL errorMessage("oorb / simplex", &
                      "TRACE BACK (85)", 1)
                 STOP
              END IF
              CALL propagate(orb_arr_(k), t)
              IF (error) THEN
                 CALL errorMessage("oorb / simplex", &
                      "TRACE BACK (90)", 1)
                 STOP
              END IF
              CALL setParameters(orb_arr_(k), &
                   dyn_model=dyn_model, &
                   perturbers=perturbers, &
                   asteroid_perturbers=asteroid_perturbers, &
                   integrator=integrator, &
                   integration_step=integration_step)
              IF (error) THEN
                 CALL errorMessage("oorb / simplex", &
                      "TRACE BACK (95)", 1)
                 STOP
              END IF
           END DO
           CALL simplexOrbits(storb, orb_arr_, &
                force_earth_impact_at_epoch=force_earth_impact_at_epoch)
           IF (.NOT.error) THEN
              EXIT
           ELSE IF (j + 13 < norb) THEN
              error = .FALSE.
           ELSE
              EXIT
           END IF
        END DO
        IF (error) THEN
           CALL errorMessage("oorb / simplex", &
                "Simplex failed:", 1)
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
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (120)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (125)", 1)
              STOP
           END IF
        ELSE
           obs_masks => getObservationMasks(storb)
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
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
                 CALL errorMessage("oorb / simplex", &
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
                 CALL errorMessage("oorb / simplex", &
                      "TRACE BACK (150)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / simplex", &
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
              HG_arr_in => getH0Distribution(physparam)
              HG_arr_storb_in(i,:,1:2) = HG_arr_in
              DEALLOCATE(HG_arr_in)
              HG_arr_in => getGDistribution(physparam)
              HG_arr_storb_in(i,:,3:4) = HG_arr_in
              DEALLOCATE(HG_arr_in)
              CALL NULLIFY(physparam)
           END IF

           CALL NEW(out_file, TRIM(out_fname) // ".smplx")
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (160)", 1)
              STOP
           END IF
           CALL setPositionAppend(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (165)", 1)
              STOP
           END IF
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (170)", 1)
              STOP
           END IF

           ! WRITE OBSERVATIONS:
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / simplex", &
                   "TRACE BACK (175)", 1)
              WRITE(getUnit(out_file),"(A)") &
                   "Could not write observations for object " // TRIM(id) 
              error = .FALSE.
           END IF
           WRITE(getUnit(out_file),"(A)") "#"

           CALL NULLIFY(out_file)

           ! WRITE ORBITS
           orb_arr_cmp => getSampleOrbits(storb)
           IF (error) THEN
              CALL errorMessage("oorb / simplex ", &
                   "TRACE BACK (175)", 1)
              STOP
           END IF
           rchi2_arr_cmp => getReducedChi2Distribution(storb)
           IF (error) THEN
              CALL errorMessage("oorb / simplex ", &
                   "TRACE BACK (185)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NEW(out_file, TRIM(id) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / simplex", &
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
                      rchi2=rchi2_arr_cmp(j), &
                      H=HG_arr_storb_in(i,j,1), &
                      G=HG_arr_storb_in(i,j,3), &
                      mjd=mjd_epoch)
              ELSE IF (orbit_format_out == "des") THEN
                 CALL errorMessage("oorb / simplex ", &
                      "DES format not yet supported for Simplex output.", 1)
                 STOP                 
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / simplex ", &
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
                 CALL errorMessage("oorb / simplex", &
                      "TRACE BACK (225)", 1)
                 STOP
              END IF
              CALL writeResiduals(storb, obss_sep(i), getUnit(out_file), compute=.TRUE.)
              IF (error) THEN
                 CALL errorMessage("oorb / simplex", &
                      "TRACE BACK (230)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
           END IF

           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Object", &
                   TRIM(id), "successfully processed."
           END IF
        END IF
        DEALLOCATE(orb_arr, stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb / simplex", &
                "Could not deallocate memory (10)", 1)
           STOP
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb > 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
     END DO
     DEALLOCATE(orb_arr_, stat=err)
     IF (err /= 0) THEN
        CALL errorMessage("oorb / simplex", &
             "Could not deallocate memory (15)", 1)
        STOP
     END IF




  CASE ("observation_sampling")

     ! Orbit inversion using a combination of MCMC sampling in
     ! observation space and optimization by simplex

     accwin_multiplier = -1.0_bp
     chi2_min_init = -1.0_bp
     CALL NULLIFY(epoch)
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model_init=dyn_model_init, &
          integrator_init=integrator_init, &
          integration_step_init=integration_step_init, &
          generat_multiplier=generat_multiplier, &
          generat_gaussian_deviates=generat_gaussian_deviates, &
          accwin_multiplier=accwin_multiplier, &
          smplx_niter=smplx_niter, &
          smplx_similarity_tol=smplx_similarity_tol, &
          smplx_tol=smplx_tol, &
          smplx_force=smplx_force, &
          os_norb=os_norb, &
          os_ntrial=os_ntrial, &
          os_sampling_type=os_sampling_type, &
          chi2_min=chi2_min_init)
     IF (error) THEN
        CALL errorMessage("oorb / observation_sampling", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        CALL NULLIFY(epoch)
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     DEALLOCATE(HG_arr_storb_in, stat=err)
     ALLOCATE(HG_arr_storb_in(SIZE(obss_sep),os_norb,4))
     HG_arr_storb_in = 99.9_bp
     ! Print header before printing first orbit:
     first = .TRUE.
     DO i=1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "TRACE BACK (25)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")."
        END IF
        nobs = getNrOfObservations(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        norb = 0
        IF (ALLOCATED(storb_arr_in)) THEN
           DO j=1,SIZE(id_arr_storb_in,dim=1)
              IF (id_arr_storb_in(j) == id) THEN
                 IF (containsDiscretePDF(storb_arr_in(j))) THEN
                    orb_arr => getSampleOrbits(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / observation_sampling", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = SIZE(orb_arr)
                    EXIT
                 ELSE
                    ALLOCATE(orb_arr(1))
                    orb_arr(1) = getNominalOrbit(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / observation_sampling", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = 1
                    EXIT                    
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / observation_sampling", &
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
                    CALL errorMessage("oorb / observation_sampling", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "Initial orbit not available.", 1)
              STOP
           ELSE
              orb_arr => reallocate(orb_arr,norb)
           END IF
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           IF (get_cl_option("--exact-mid-epoch", .FALSE.)) THEN
              mjd = mjd+dt/2.0_bp
           ELSE
              mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           END IF
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (70)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL NEW(storb, obss_sep(i))        
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "TRACE BACK (75)", 1)
           STOP
        END IF
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             dchi2_rejection = dchi2_rejection, &
             dchi2_max = dchi2_max, &
             chi2_min=chi2_min_init, &
             generat_multiplier=generat_multiplier, &
             generat_gaussian_deviates=generat_gaussian_deviates, &
             accept_multiplier=accwin_multiplier, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             smplx_niter=smplx_niter, &
             smplx_force=smplx_force, &
             smplx_tol=smplx_tol, &
             smplx_similarity_tol=smplx_similarity_tol, &
             os_norb=os_norb, &
             os_ntrial=os_ntrial, &
             os_sampling_type=os_sampling_type)
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "TRACE BACK (80)", 1)
           STOP
        END IF
        DO k=1,SIZE(orb_arr)
           CALL setParameters(orb_arr(k), &
                dyn_model=dyn_model_init, &
                perturbers=perturbers, &
                asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator_init, &
                integration_step=integration_step_init)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (85)", 1)
              STOP
           END IF
           CALL propagate(orb_arr(k), t)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (90)", 1)
              STOP
           END IF
           CALL setParameters(orb_arr(k), &
                dyn_model=dyn_model, &
                perturbers=perturbers, &
                asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (95)", 1)
              STOP
           END IF
        END DO
        CALL observationSampling(storb, orb_arr)
        IF (error) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "Observation_sampling failed:", 1)
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
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (120)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (125)", 1)
              STOP
           END IF
        ELSE
           obs_masks => getObservationMasks(storb)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
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
                 CALL errorMessage("oorb / observation_sampling", &
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
                 CALL errorMessage("oorb / observation_sampling", &
                      "TRACE BACK (150)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling", &
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
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling", &
                      "TRACE BACK (156)", 1)
                 STOP
              END IF
              HG_arr_in => getH0Distribution(physparam)
              HG_arr_storb_in(i,:,1:2) = HG_arr_in
              DEALLOCATE(HG_arr_in)
              HG_arr_in => getGDistribution(physparam)
              HG_arr_storb_in(i,:,3:4) = HG_arr_in
              DEALLOCATE(HG_arr_in)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling", &
                      "TRACE BACK (157)", 1)
                 STOP
              END IF
              CALL NULLIFY(physparam)
           END IF

           CALL NEW(out_file, TRIM(out_fname) // ".os")
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (160)", 1)
              STOP
           END IF
           CALL setPositionAppend(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (165)", 1)
              STOP
           END IF
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (170)", 1)
              STOP
           END IF

           ! WRITE OBSERVATIONS:
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (175)", 1)
              WRITE(getUnit(out_file),"(A)") &
                   "Could not write observations for object " // TRIM(id) 
              error = .FALSE.
           END IF
           WRITE(getUnit(out_file),"(A)") "#"

           CALL NULLIFY(out_file)

           ! WRITE ORBITS
           orb_arr_cmp => getSampleOrbits(storb)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling ", &
                   "TRACE BACK (175)", 1)
              STOP
           END IF
           rchi2_arr_cmp => getReducedChi2Distribution(storb)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling ", &
                   "TRACE BACK (185)", 1)
              STOP
           END IF
           pdf_arr_cmp => getDiscretePDF(storb, element_type_comp_prm)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "TRACE BACK (186)", 1)
              STOP
           END IF
           CALL getResults(storb, repetition_arr_cmp=repetition_arr_cmp)
           IF (error) THEN
              CALL errorMessage("oorb / observation_sampling ", &
                   "TRACE BACK (187)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NEW(out_file, TRIM(id) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling", &
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
                      rchi2=rchi2_arr_cmp(j), &
                      element_type_pdf=element_type_comp_prm, &
                      pdf=pdf_arr_cmp(j), &
                      H=HG_arr_storb_in(i,j,1), &
                      G=HG_arr_storb_in(i,j,3), &
                      mjd=mjd_epoch, &
                      repetitions=repetition_arr_cmp(j))
              ELSE IF (orbit_format_out == "des") THEN
                 CALL errorMessage("oorb / observation_sampling ", &
                      "DES format not yet supported for observation_sampling output.", 1)
                 STOP                 
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling ", &
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
                 CALL errorMessage("oorb / observation_sampling", &
                      "TRACE BACK (225)", 1)
                 STOP
              END IF
              CALL writeResiduals(storb, obss_sep(i), getUnit(out_file), compute=.TRUE.)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling", &
                      "TRACE BACK (230)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
           END IF

           IF (plot_results) THEN
              CALL toString(dt, str, error, frmt="(F10.2)")
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling ", &
                      "TRACE BACK (235)", 1)
                 STOP
              END IF
              str = TRIM(id) // "_"// TRIM(str)
              CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
                   "_os_residual_stamps.eps", compute=.TRUE.)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling", &
                      "TRACE BACK (240)", 1)
                 STOP
              END IF
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_os_residual_stamps.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_os_residual_stamps.eps* &")
              END IF
              ALLOCATE(elements_arr(SIZE(orb_arr_cmp,dim=1),7), stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / observation_sampling", &
                      "Could not allocate memory (3).", 1)
                 STOP
              END IF
              CALL NEW(tmp_file, TRIM(str)// "_os_orbits.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling ", &
                      "TRACE BACK (250)", 1)
                 STOP
              END IF
              DO j=1,SIZE(orb_arr_cmp,dim=1)
                 IF (element_type_comp_prm == "cartesian") THEN
                    CALL rotateToEcliptic(orb_arr_cmp(j))
                 END IF
                 elements_arr(j,1:6) = getElements(orb_arr_cmp(j), element_type_comp_prm)
                 IF (error) THEN
                    CALL errorMessage("oorb / observation_sampling", &
                         "TRACE BACK (255)", 1)
                    STOP
                 END IF
                 IF (element_type_comp_prm == "keplerian") THEN
                    elements_arr(j,3:6) = elements_arr(j,3:6)/rad_deg
                 END IF
                 elements_arr(j,7) = pdf_arr_cmp(j)
                 t = getTime(orb_arr_cmp(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / observation_sampling ", &
                         "TRACE BACK (260)", 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmp_file),"(7(E23.15,1X),A)") &
                      elements_arr(j,1:6), &
                      pdf_arr_cmp(j), &
                      getCalendarDateString(t,"tdt")
                 IF (error) THEN
                    CALL errorMessage("oorb / observation_sampling ", &
                         "TRACE BACK (265)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(t)
              END DO
              CALL NULLIFY(tmp_file)
              CALL NEW(tmp_file, TRIM(str) // &
                   "_os_sample_standard_deviations.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling ", &
                      "TRACE BACK (275)", 1)
                 STOP
              END IF
              WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                   advance="no") &
                   getObservationalTimespan(obss_sep(i))
              IF (error) THEN
                 CALL errorMessage("oorb / observation_sampling ", &
                      "TRACE BACK (280)", 1)
                 STOP
              END IF
              DO j=1,6
                 CALL moments(elements_arr(:,j), &
                      pdf=pdf_arr_cmp, std_dev=stdev, errstr=errstr)
                 IF (LEN_TRIM(errstr) /= 0) THEN
                    CALL errorMessage("oorb / observation_sampling", &
                         "Could not compute moments. " // TRIM(errstr), 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                      advance="no") stdev
              END DO
              WRITE(getUnit(tmp_file),*)
              CALL NULLIFY(tmp_file)
              DEALLOCATE(elements_arr, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / observation_sampling", &
                      "Could not deallocate memory (5).", 1)
                 STOP
              END IF
              ! Make plot using gnuplot:
              CALL system("cp " // TRIM(str) // &
                   "_os_orbits.out sor_orbits.out")
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_car.gp")
              ELSE
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_kep.gp")
              END IF
              CALL system("cp sor_results.eps " // TRIM(str) // &
                   "_os_" // TRIM(element_type_comp_prm) // &
                   "_results.eps")
              CALL system("rm -f sor_orbits.out sor_results.eps " // & 
                   TRIM(str) // "_os_orbits.out " // TRIM(str) // &
                   "_os_sample_standard_deviations.out")
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_os_" // &
                      TRIM(element_type_comp_prm) // "_results.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_os_" // &
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
           DEALLOCATE(repetition_arr_cmp, stat=err)           
           IF (err /= 0) THEN
              CALL errorMessage("oorb / observation_sampling", &
                   "Could not deallocate memory (10).", 1)
              STOP
           END IF
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Object", &
                   TRIM(id), "successfully processed."
           END IF
        END IF
        DEALLOCATE(orb_arr, stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb / observation_sampling", &
                "Could not deallocate memory (10)", 1)
           STOP
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb > 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
     END DO


  CASE ("vov")

     ! Orbit inversion using the 6D phase-space volume-of-variation
     ! sampling.

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
          ls_rchi2_acceptable=ls_rchi2_acceptable, &
          vov_type=vov_type_prm, &
          vov_norb=vov_norb, &
          vov_ntrial=vov_ntrial, &
          vov_niter=vov_niter, &
          vov_norb_iter=vov_norb_iter, &
          vov_ntrial_iter=vov_ntrial_iter, &
          vov_nmap=vov_nmap, &
          vov_mapping_mask=vov_mapping_mask, &
          vov_scaling=vov_scaling_prm, &
          accwin_multiplier=accwin_multiplier)
     IF (error) THEN
        CALL errorMessage("oorb / vov", &
             "TRACE BACK (2)", 1)
        STOP
     END IF
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        CALL NULLIFY(epoch)
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / vov", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     DO i=j+1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / vov", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")"
        END IF
        norb = 0
        norb = 0
        IF (ALLOCATED(storb_arr_in)) THEN
           DO j=1,SIZE(id_arr_storb_in,dim=1)
              IF (id_arr_storb_in(j) == id) THEN
                 IF (containsDiscretePDF(storb_arr_in(j))) THEN
                    orb_arr => getSampleOrbits(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / vov", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = SIZE(orb_arr)
                    EXIT
                 ELSE
                    ALLOCATE(orb_arr(1))
                    orb_arr(1) = getNominalOrbit(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / vov", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = 1
                    EXIT                    
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / vov", &
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
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / vov", &
                   "Initial orbit not available.", 1)
              STOP
           ELSE
              orb_arr => reallocate(orb_arr,norb)
           END IF
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / vov", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (70)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL propagate(orb_arr, t)
        IF (error) THEN
           CALL errorMessage("oorb / vov", &
                "TRACE BACK (65)", 1)
           STOP
        END IF
        DO k=1,norb
           CALL setParameters(orb_arr(k), &
                dyn_model=dyn_model, &
                perturbers=perturbers, &
                asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (35)", 1)
              STOP
           END IF
        END DO
        CALL NEW(storb, obss_sep(i))        
        IF (error) THEN
           CALL errorMessage("oorb / vov", &
                "TRACE BACK (70)", 1)
           STOP
        END IF
        vov_k: DO k=1,MIN(10,norb)
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(1X,A,1X,I0,2(A,1X))") "Trying", &
                   k, ". orbit for object", TRIM(id)
           END IF
           CALL setParameters(storb, &
                dyn_model=dyn_model, &
                perturbers=perturbers, &
                asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, &
                integration_step=integration_step, &
                outlier_rejection=outlier_rejection_prm, &
                outlier_multiplier=outlier_multiplier_prm, &
                element_type=element_type_comp_prm, &
                accept_multiplier=accwin_multiplier, &
                vov_norb=vov_norb, &
                vov_ntrial=vov_ntrial, &
                vov_niter=vov_niter, &
                vov_norb_iter=vov_norb_iter, &
                vov_ntrial_iter=vov_ntrial_iter, &
                vov_nmap=vov_nmap, &
                vov_mapping_mask=vov_mapping_mask, &
                vov_scaling=vov_scaling_prm, &
                ls_correction_factor=ls_correction_factor, &
                ls_element_mask=ls_element_mask, &
                ls_rchi2_acceptable=ls_rchi2_acceptable)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (75)", 1)
              STOP
           END IF
           SELECT CASE (vov_type_prm)
           CASE (1)
              CALL volumeOfVariation(storb, orb_arr(k))
           CASE (2)
              CALL autoVolumeOfVariation(storb, orb_arr(k))
           CASE default
              CALL errorMessage("oorb / vov", &
                   "Unknown type of VoV:", 1)
              IF (err_verb >= 1) THEN
                 WRITE(stderr,"(2X,I0)") vov_type_prm
              END IF
              STOP
           END SELECT
           IF (error .OR. getNrOfSampleOrbits(storb) < INT(0.5*vov_norb)) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (80)", 1)
              error  = .FALSE.
              CALL NEW(out_file, "problematic_observation_sets." // &
                   TRIM(observation_format_out))
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (85)", 1)
                 STOP
              END IF
              CALL setPositionAppend(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (90)", 1)
                 STOP
              END IF
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (95)", 1)
                 STOP
              END IF
              CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                   TRIM(observation_format_out))
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (100)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (105)", 1)
                 STOP
              END IF
              CYCLE vov_k
           END IF
           ! In case of too many outliers and the 2-b
           ! approximation, throw the observation set in a separate
           ! bin:
           IF (dyn_model == "2-body") THEN
              obs_masks => getObservationMasks(storb)
              noutlier = 0
              DO l=1,SIZE(obs_masks,dim=1)
                 IF (ALL(.NOT.obs_masks(l,:))) THEN
                    noutlier = noutlier + 1
                 END IF
              END DO
              IF (noutlier > 0) THEN
                 CALL NEW(out_file, "outliers." // TRIM(observation_format_out))
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (85)", 1)
                    STOP
                 END IF
                 CALL setPositionAppend(out_file)
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (90)", 1)
                    STOP
                 END IF
                 CALL OPEN(out_file)
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (95)", 1)
                    STOP
                 END IF
                 WRITE(getUnit(out_file),"(1X)")
                 WRITE(getUnit(out_file),"(A,I0,A)",advance="no") "# ", &
                      noutlier, " outliers: "
                 DO l=1,SIZE(obs_masks,dim=1)
                    IF (ALL(.NOT.obs_masks(l,:))) THEN
                       WRITE(getUnit(out_file),"(A)",advance="no") "*"
                    ELSE
                       WRITE(getUnit(out_file),"(A)",advance="no") "-"
                    END IF
                 END DO
                 WRITE(getUnit(out_file),"(1X)")
                 CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                      TRIM(observation_format_out))
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (100)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(out_file)
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (105)", 1)
                    STOP
                 END IF
              END IF
              DEALLOCATE(obs_masks, stat=err)
           END IF
           CALL NEW(out_file, TRIM(id) // ".vov")
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (10)", 1)
              STOP
           END IF
           ! WRITE OBSERVATIONS:
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (85)", 1)
              STOP
           END IF
           WRITE(getUnit(out_file),"(A)") "#"
           CALL writeVOVResults(storb, obss_sep(i), getUnit(out_file))
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK ()", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)

           ! SAVE NOMINAL ORBIT INFORMATION TO A "LS" OUTPUT FILE?
           CALL NEW(out_file, TRIM(id) // ".nominal")
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (10)", 1)
              STOP
           END IF
           CALL writeNominalSolution(storb, obss_sep(i), &
                element_type_out_prm, getUnit(out_file))
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK ()", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           orb_arr_cmp => getSampleOrbits(storb)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (120)", 1)
              STOP
           END IF
           pdf_arr_cmp => getDiscretePDF(storb)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (125)", 1)
              STOP
           END IF
           rchi2_arr_cmp => getReducedChi2Distribution(storb)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (130)", 1)
              STOP
           END IF
           CALL getResults(storb, &
                reg_apr_arr=reg_apr_arr_cmp, &
                jac_arr=jac_arr_cmp)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (135)", 1)
              STOP
           END IF
           CALL NEW(orb_out_file, TRIM(id) // ".orb")
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (15)", 1)
              STOP
           END IF
           CALL OPEN(orb_out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vov", &
                   "TRACE BACK (20)", 1)
              STOP
           END IF
           DO l=1,SIZE(orb_arr_cmp,dim=1)
              CALL writeOpenOrbOrbitFile(getUnit(orb_out_file), &
                   print_header=l==1.AND.i==1, &
                   element_type_out=element_type_out_prm, &
                   id=id, &
                   orb=orb_arr_cmp(l), &
                   element_type_pdf=element_type_comp_prm, &
                   pdf=pdf_arr_cmp(l), &
                   rchi2=rchi2_arr_cmp(l), &
                   reg_apr=reg_apr_arr_cmp(l), &
                   jac_car_kep=jac_arr_cmp(l,2), &
                   jac_equ_kep=jac_arr_cmp(l,3))
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (140)", 1)
                 STOP
              END IF
           END DO
           CALL NULLIFY(orb_out_file)
           IF (plot_results) THEN
              CALL toString(dt, str, error, frmt="(F10.2)")
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (145)", 1)
                 STOP
              END IF
              str = TRIM(id) // "_"// TRIM(str)
              CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
                   "_vov_residual_stamps.ps",  compute=.TRUE.)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (142)", 1)
                 STOP
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_vov_residual_stamps.ps &")
              END IF
              CALL getResults(storb, &
                   vov_map_cmp=vov_map, &
                   vov_scaling_cmp=vov_scaling_cmp)
              !vov_mapping_mask_prm=vov_mapping_mask)
              CALL getParameters(storb, &
                   vov_mapping_mask=vov_mapping_mask)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (135)", 1)
                 STOP
              END IF
              ALLOCATE(elements_arr(SIZE(orb_arr_cmp,dim=1),7), stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / vov", &
                      "Could not allocate memory (10)", 1)
                 STOP
              END IF
              CALL NEW(tmp_file, TRIM(str)// "_vov_orbits.out")
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (150)", 1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (155)", 1)
                 STOP
              END IF
              DO l=1,SIZE(orb_arr_cmp,dim=1)
                 IF (element_type_comp_prm == "cartesian") THEN
                    CALL rotateToEcliptic(orb_arr_cmp(l))
                 END IF
                 elements_arr(l,1:6) = getElements(orb_arr_cmp(l), element_type_comp_prm)
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (160)", 1)
                    STOP
                 END IF
                 elements_arr(l,7) = pdf_arr_cmp(l)
                 t = getTime(orb_arr_cmp(l))
                 IF (error) THEN
                    CALL errorMessage("oorb / vov", &
                         "TRACE BACK (165)", 1)
                    STOP
                 END IF
                 IF (element_type_comp_prm == "keplerian") THEN
                    WRITE(getUnit(tmp_file),*) &
                         elements_arr(l,1:2), &
                         elements_arr(l,3:6)/rad_deg, &
                         pdf_arr_cmp(l), &
                         getCalendarDateString(t,"tdt")
                 ELSE
                    WRITE(getUnit(tmp_file),*) &
                         elements_arr(l,1:6), &
                         pdf_arr_cmp(l), &
                         getCalendarDateString(t,"tdt")
                 END IF
                 CALL NULLIFY(t)
              END DO
              CALL NULLIFY(tmp_file)
              CALL NULLIFY(orb)
              orb = getNominalOrbit(storb)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (170)", 1)
                 STOP
              END IF
              elements = getElements(orb, element_type_comp_prm)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (175)", 1)
                 STOP
              END IF
              t = getTime(orb)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (180)", 1)
                 STOP
              END IF
              CALL NULLIFY(orb)
              CALL NEW(tmp_file, TRIM(str)// "_vov_nominal_orbit.out")
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (185)", 1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (190)", 1)
                 STOP
              END IF
              IF (element_type_comp_prm == "keplerian") THEN
                 WRITE(getUnit(tmp_file),*) elements(1:2), &
                      elements(3:6)/rad_deg, &
                      getCalendarDateString(t,"tdt")
              ELSE
                 WRITE(getUnit(tmp_file),*) elements(1:6), &
                      getCalendarDateString(t,"tdt")
              END IF
              CALL NULLIFY(tmp_file)
              CALL NULLIFY(t)
              vov_nmap = SIZE(vov_map,dim=1)
              CALL NEW(tmp_file, TRIM(str)// "_vov_sampling_grid.out")
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (195)", 1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (200)", 1)
                 STOP
              END IF
              DO l=1,vov_nmap,vov_nmap/10
                 IF (element_type_comp_prm == "keplerian") THEN
                    WRITE(getUnit(tmp_file), "(6(F22.15,1X))", &
                         advance="no") vov_map(l,1:2), &
                         vov_map(l,3:6)/rad_deg
                 ELSE
                    WRITE(getUnit(tmp_file), "(6(F22.15,1X))", &
                         advance="no") vov_map(l,1:6)
                 END IF
                 lower_limit = vov_map(l,1:6) - &
                      vov_scaling_cmp(:,1)*vov_map(l,7:12)
                 upper_limit = vov_map(l,1:6) + &
                      vov_scaling_cmp(:,2)*vov_map(l,7:12)
                 DO m=1,6
                    IF (vov_mapping_mask(m)) THEN
                       CYCLE
                    END IF
                    IF (element_type_comp_prm == "keplerian" .AND. m>=3) THEN
                       WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                            advance="no") lower_limit(m)/rad_deg
                       WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                            advance="no") upper_limit(m)/rad_deg
                    ELSE
                       WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                            advance="no") lower_limit(m)
                       WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                            advance="no") upper_limit(m)
                    END IF
                 END DO
                 WRITE(getUnit(tmp_file),*)
              END DO
              CALL NULLIFY(tmp_file)
              CALL NEW(tmp_file, TRIM(str) // &
                   "_vov_sample_standard_deviations.out")
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (205)", 1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (210)", 1)
                 STOP
              END IF
              WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                   advance="no") &
                   getObservationalTimespan(obss_sep(i))
              IF (error) THEN
                 CALL errorMessage("oorb / vov", &
                      "TRACE BACK (215)", 1)
                 STOP
              END IF
              DO l=1,6
                 CALL moments(elements_arr(:,l), &
                      pdf=pdf_arr_cmp, std_dev=stdev, errstr=errstr)
                 IF (LEN_TRIM(errstr) /= 0) THEN
                    CALL errorMessage("oorb / vov", &
                         "Could not compute moments: " // TRIM(errstr), 1)
                    STOP
                 END IF
                 IF (element_type_comp_prm == "keplerian" .AND. l>=3) THEN
                    WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                         advance="no")  stdev/rad_deg
                 ELSE
                    WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                         advance="no") stdev
                 END IF
              END DO
              WRITE(getUnit(tmp_file),*)
              CALL NULLIFY(tmp_file)
              DEALLOCATE(elements_arr, vov_map, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / vov", &
                      "Could not deallocate memory (10)", 1)
                 STOP
              END IF
              ! Make plot using gnuplot:
              CALL system("cp " // TRIM(str) // &
                   "_vov_sampling_grid.out vov_sampling_grid.out")
              CALL system("cp " // TRIM(str) // &
                   "_vov_orbits.out vov_orbits.out")
              CALL system("cp " // TRIM(str) // &
                   "_vov_nominal_orbit.out vov_nominal_orbit.out")
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/vov_plot_car.gp")
              ELSE
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/vov_plot_kep.gp")
              END IF
              CALL system("cp vov_results.ps " // TRIM(str) // &
                   "_vov_" // TRIM(element_type_comp_prm) // &
                   "_results.ps")
              CALL system("rm -f vov_sampling_grid.out vov_orbits.out " // &
                   "vov_nominal_orbit.out vov_results.ps")
              CALL system("rm -f " // TRIM(str) // &
                   "_vov_sampling_grid.out " // TRIM(str) &
                   // "_vov_orbits.out " // TRIM(str) // &
                   "_vov_nominal_orbit.out " // TRIM(str) // &
                   "_vov_sample_standard_deviations.out")
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_vov_" // &
                      TRIM(element_type_comp_prm) // &
                      "_results.ps &")
              END IF
           END IF
           DEALLOCATE(orb_arr_cmp, pdf_arr_cmp, rchi2_arr_cmp, &
                reg_apr_arr_cmp, jac_arr_cmp, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("oorb / vov", &
                   "Could not deallocate memory (10).", 1)
              STOP
           END IF
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Object", &
                   TRIM(id), "successfully processed."
           END IF
           ! No more trials needed:
           EXIT
        END DO vov_k
        DEALLOCATE(orb_arr, stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb / vov", &
                "Could not deallocate memory (15)", 1)
           STOP
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb > 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
     END DO

  CASE ("vomcmc")

     ! Orbit inversion using the virtual-observation MCMC sampling.

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
          ls_rchi2_acceptable=ls_rchi2_acceptable, &
          vomcmc_type=vomcmc_type_prm, &
          vomcmc_norb=vomcmc_norb, &
          vomcmc_ntrial=vomcmc_ntrial, &
          vomcmc_niter=vomcmc_niter, &
          vomcmc_norb_iter=vomcmc_norb_iter, &
          vomcmc_ntrial_iter=vomcmc_ntrial_iter, &
          vomcmc_nmap=vomcmc_nmap, &
          vomcmc_mapping_mask=vomcmc_mapping_mask, &
          vomcmc_scaling=vomcmc_scaling_prm, &
          accwin_multiplier=accwin_multiplier, &
          generat_multiplier=generat_multiplier, &
          generat_gaussian_deviates=generat_gaussian_deviates, &
          os_norb=os_norb, &
          os_ntrial=os_ntrial, &
          os_sampling_type=os_sampling_type, &
          smplx_niter=smplx_niter, &
          smplx_similarity_tol=smplx_similarity_tol)
     IF (error) THEN
        CALL errorMessage("oorb / vomcmc", &
             "TRACE BACK (2)", 1)
        STOP
     END IF

     DO i=1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")"
        END IF
        norb = 0
        IF (ALLOCATED(storb_arr_in)) THEN
           DO j=1,SIZE(id_arr_storb_in,dim=1)
              IF (id_arr_storb_in(j) == id) THEN
                 IF (containsDiscretePDF(storb_arr_in(j))) THEN
                    orb_arr => getSampleOrbits(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / vomcmc", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = SIZE(orb_arr)
                    EXIT
                 ELSE
                    ALLOCATE(orb_arr(1))
                    orb_arr(1) = getNominalOrbit(storb_arr_in(j))
                    IF (error) THEN
                       CALL errorMessage("oorb / vomcmc", &
                            "TRACE BACK (35)", 1)
                       STOP
                    END IF
                    norb = 1
                    EXIT                    
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "Initial orbit not available (1).", 1)
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
                    CALL errorMessage("oorb / vomcmc", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
              END IF
           END DO
           IF (norb == 0) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "Initial orbit not available (2).", 1)
              STOP
           ELSE
              orb_arr => reallocate(orb_arr,norb)
           END IF
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (70)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL propagate(orb_arr, t)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (65)", 1)
           STOP
        END IF
        DO k=1,norb
           CALL setParameters(orb_arr(k), &
                dyn_model=dyn_model, &
                perturbers=perturbers, &
                asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (35)", 1)
              STOP
           END IF
        END DO
        CALL NEW(storb, obss_sep(i))        
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (70)", 1)
           STOP
        END IF
        CALL setParameters(storb, &
             t_inv=t, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             dchi2_rejection = dchi2_rejection, &
             dchi2_max = dchi2_max, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             element_type=element_type_comp_prm, &
             accept_multiplier=accwin_multiplier, &
             vomcmc_norb=vomcmc_norb, &
             vomcmc_ntrial=vomcmc_ntrial, &
             vomcmc_niter=vomcmc_niter, &
             vomcmc_norb_iter=vomcmc_norb_iter, &
             vomcmc_ntrial_iter=vomcmc_ntrial_iter, &
             vomcmc_nmap=vomcmc_nmap, &
             vomcmc_mapping_mask=vomcmc_mapping_mask, &
             vomcmc_scaling=vomcmc_scaling_prm, &
             generat_multiplier=generat_multiplier, &
             generat_gaussian_deviates=generat_gaussian_deviates, &
             os_norb=os_norb, &
             os_ntrial=os_ntrial, &
             os_sampling_type=os_sampling_type, &
             smplx_niter=smplx_niter, &
             smplx_force=.FALSE., &
             smplx_similarity_tol=smplx_similarity_tol)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (75)", 1)
           STOP
        END IF
        SELECT CASE (vomcmc_type_prm)
        CASE (1,2)
           CALL virtualObservationMCMC(storb, orb_arr)
        CASE default
           CALL errorMessage("oorb / vomcmc", &
                "Unknown type of Vomcmc:", 1)
           IF (err_verb >= 1) THEN
              WRITE(stderr,"(2X,I0)") vomcmc_type_prm
           END IF
           STOP
        END SELECT
        IF (error) THEN! .OR. getNrOfSampleOrbits(storb) < INT(0.5*vomcmc_norb)) THEN
           WRITE(stderr,*) error, getNrOfSampleOrbits(storb), INT(0.5*vomcmc_norb)
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (80)", 1)
           error  = .FALSE.
           CALL NEW(out_file, "problematic_observation_sets." // &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (85)", 1)
              STOP
           END IF
           CALL setPositionAppend(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (90)", 1)
              STOP
           END IF
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (95)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (100)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (105)", 1)
              STOP
           END IF
           CALL NULLIFY(storb)
           CYCLE
        END IF
        ! In case of too many outliers and the 2-b
        ! approximation, throw the observation set in a separate
        ! bin:
        IF (dyn_model == "2-body") THEN
           obs_masks => getObservationMasks(storb)
           noutlier = 0
           DO l=1,SIZE(obs_masks,dim=1)
              IF (ALL(.NOT.obs_masks(l,:))) THEN
                 noutlier = noutlier + 1
              END IF
           END DO
           IF (noutlier > 0) THEN
              CALL NEW(out_file, "outliers." // TRIM(observation_format_out))
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (85)", 1)
                 STOP
              END IF
              CALL setPositionAppend(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (90)", 1)
                 STOP
              END IF
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (95)", 1)
                 STOP
              END IF
              WRITE(getUnit(out_file),"(1X)")
              WRITE(getUnit(out_file),"(A,I0,A)",advance="no") "# ", &
                   noutlier, " outliers: "
              DO l=1,SIZE(obs_masks,dim=1)
                 IF (ALL(.NOT.obs_masks(l,:))) THEN
                    WRITE(getUnit(out_file),"(A)",advance="no") "*"
                 ELSE
                    WRITE(getUnit(out_file),"(A)",advance="no") "-"
                 END IF
              END DO
              WRITE(getUnit(out_file),"(1X)")
              CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                   TRIM(observation_format_out))
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (100)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (105)", 1)
                 STOP
              END IF
           END IF
           DEALLOCATE(obs_masks, stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "Could not deallocate memory (5)", 1)
              STOP
           END IF
        END IF
        CALL NEW(out_file, TRIM(id) // ".vomcmc")
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (115)", 1)
           STOP
        END IF
        CALL OPEN(out_file)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
        ! WRITE OBSERVATIONS:
        CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
             TRIM(observation_format_out)) 
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (85)", 1)
           STOP
        END IF
        WRITE(getUnit(out_file),"(A)") "#"
        CALL writeVOMCMCResults(storb, obss_sep(i), getUnit(out_file))
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK ()", 1)
           STOP
        END IF
        CALL NULLIFY(out_file)
        ! WRITE SAMPLE ORBITS TO OUTPUT FILE
        orb_arr_cmp => getSampleOrbits(storb)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (120)", 1)
           STOP
        END IF
        pdf_arr_cmp => getDiscretePDF(storb, element_type_comp_prm)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (125)", 1)
           STOP
        END IF
        rchi2_arr_cmp => getReducedChi2Distribution(storb)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (130)", 1)
           STOP
        END IF
        CALL getResults(storb, repetition_arr_cmp=repetition_arr_cmp)
        IF (error) THEN
           CALL errorMessage("oorb / vomcmc", &
                "TRACE BACK (135)", 1)
           STOP
        END IF
        IF (separately) THEN
           CALL NEW(out_file, TRIM(id) // ".orb")
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (200)", 1)
              STOP
           END IF
           lu_orb_out = getUnit(out_file)
        END IF
        DO l=1,SIZE(orb_arr_cmp,dim=1)
           CALL writeOpenOrbOrbitFile(lu_orb_out, &
                print_header=l==1.AND.i==1, &
                element_type_out=element_type_out_prm, &
                id=id, &
                orb=orb_arr_cmp(l), &
                element_type_pdf=element_type_comp_prm, &
                pdf=pdf_arr_cmp(l), &
                rchi2=rchi2_arr_cmp(l), &
                repetitions=repetition_arr_cmp(l), &
                mjd=mjd_epoch)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (140)", 1)
              STOP
           END IF
        END DO
        IF (separately) THEN
           CALL NULLIFY(out_file)
           IF (compress) THEN
              CALL system("gzip -f " // TRIM(id) // ".orb")
           END IF
        END IF
        IF (plot_results) THEN
           CALL toString(dt, str, error, frmt="(F10.2)")
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (145)", 1)
              STOP
           END IF
           str = TRIM(id) // "_"// TRIM(str)
           CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
                "_vomcmc_residual_stamps.eps",  compute=.TRUE.)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (142)", 1)
              STOP
           END IF
           IF (plot_open) THEN
              CALL system("gv " // TRIM(str) // "_vomcmc_residual_stamps.eps &")
           END IF
!!$           CALL getResults(storb, &
!!$                vomcmc_map_cmp=vomcmc_map, &
!!$                vomcmc_scaling_cmp=vomcmc_scaling_cmp)
           !vomcmc_mapping_mask_prm=vomcmc_mapping_mask)
           CALL getParameters(storb, &
                vomcmc_mapping_mask=vomcmc_mapping_mask)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (135)", 1)
              STOP
           END IF
           ALLOCATE(elements_arr(SIZE(orb_arr_cmp,dim=1),7), stat=err)
           IF (err /= 0) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "Could not allocate memory (10)", 1)
              STOP
           END IF
           CALL NEW(tmp_file, TRIM(str)// "_vomcmc_orbits.out")
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (150)", 1)
              STOP
           END IF
           CALL OPEN(tmp_file)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (155)", 1)
              STOP
           END IF
           DO l=1,SIZE(orb_arr_cmp,dim=1)
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL rotateToEcliptic(orb_arr_cmp(l))
              END IF
              elements_arr(l,1:6) = getElements(orb_arr_cmp(l), element_type_comp_prm)
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (160)", 1)
                 STOP
              END IF
              elements_arr(l,7) = pdf_arr_cmp(l)
              t = getTime(orb_arr_cmp(l))
              IF (error) THEN
                 CALL errorMessage("oorb / vomcmc", &
                      "TRACE BACK (165)", 1)
                 STOP
              END IF
              IF (element_type_comp_prm == "keplerian") THEN
                 WRITE(getUnit(tmp_file),*) &
                      elements_arr(l,1:2), &
                      elements_arr(l,3:6)/rad_deg, &
                      pdf_arr_cmp(l), &
                      getCalendarDateString(t,"tdt")
              ELSE
                 WRITE(getUnit(tmp_file),*) &
                      elements_arr(l,1:6), &
                      pdf_arr_cmp(l), &
                      getCalendarDateString(t,"tdt")
              END IF
              CALL NULLIFY(t)
           END DO
           CALL NULLIFY(tmp_file)
           CALL NULLIFY(orb)
           orb = getNominalOrbit(storb)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (170)", 1)
              STOP
           END IF
           elements = getElements(orb, element_type_comp_prm)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (175)", 1)
              STOP
           END IF
           t = getTime(orb)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (180)", 1)
              STOP
           END IF
           CALL NULLIFY(orb)
           CALL NEW(tmp_file, TRIM(str)// "_vomcmc_nominal_orbit.out")
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (185)", 1)
              STOP
           END IF
           CALL OPEN(tmp_file)
           IF (error) THEN
              CALL errorMessage("oorb / vomcmc", &
                   "TRACE BACK (190)", 1)
              STOP
           END IF
           IF (element_type_comp_prm == "keplerian") THEN
              WRITE(getUnit(tmp_file),*) elements(1:2), &
                   elements(3:6)/rad_deg, &
                   getCalendarDateString(t,"tdt")
           ELSE
              WRITE(getUnit(tmp_file),*) elements(1:6), &
                   getCalendarDateString(t,"tdt")
           END IF
           CALL NULLIFY(tmp_file)
           CALL NULLIFY(t)
           !CALL system("cp " // TRIM(str) // &
           !     "_vomcmc_orbits.out vomcmc_orbits.out")
           CALL system("cp " // TRIM(str) // &
                "_vomcmc_orbits.out sor_orbits.out")
           !           IF (element_type_comp_prm == "cartesian") THEN
           !              CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/vomcmc_plot_car.gp")
           !           ELSE
           !              CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/vomcmc_plot_kep.gp")
           !           END IF
           !           CALL system("cp vomcmc_results.ps " // TRIM(str) // &
           !                "_vomcmc_" // TRIM(element_type_comp_prm) // &
           !                "_results.ps")
           IF (element_type_comp_prm == "cartesian") THEN
              CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_car.gp")
           ELSE
              CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/sor_plot_kep.gp")
           END IF
           CALL system("mv sor_results.eps " // TRIM(str) // &
                "_vomcmc_" // TRIM(element_type_comp_prm) // &
                "_results.eps")
           CALL system("rm -f vomcmc_sampling_grid.out vomcmc_orbits.out " // &
                "vomcmc_nominal_orbit.out vomcmc_results.ps sor_orbits.out")
           CALL system("rm -f " // TRIM(str) // &
                "_vomcmc_sampling_grid.out " // TRIM(str) &
                // "_vomcmc_orbits.out " // TRIM(str) // &
                "_vomcmc_nominal_orbit.out " // TRIM(str) // &
                "_vomcmc_sample_standard_deviations.out")
           IF (plot_open) THEN
              CALL system("gv " // TRIM(str) // "_vomcmc_" // &
                   TRIM(element_type_comp_prm) // &
                   "_results.ps &")
           END IF
        END IF
        DEALLOCATE(orb_arr_cmp, pdf_arr_cmp, rchi2_arr_cmp, &
             repetition_arr_cmp, stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb / vomcmc", &
                "Could not deallocate memory (15).", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(3(1X,A))") "Object", &
                TRIM(id), "successfully processed."
        END IF
        DEALLOCATE(orb_arr, stat=err)
        IF (err /= 0) THEN
           CALL errorMessage("oorb / vomcmc", &
                "Could not deallocate memory (20)", 1)
           STOP
        END IF
        CALL NULLIFY(storb)
        CALL NULLIFY(orb)
        IF (info_verb > 2) THEN
           WRITE(stdout,*)
           WRITE(stdout,*)
        END IF
     END DO


  CASE ("lsl")

     ! Orbital inversion using least squares with linearized
     ! covariances, that is, fixing the resulting shape of the
     ! orbital-element pdf to a multidimensional Gaussian.

     CALL NULLIFY(epoch)
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model_init=dyn_model_init, &
          integrator_init=integrator_init, &
          integration_step_init=integration_step_init, &
          ls_correction_factor=ls_correction_factor, &
          ls_rchi2_acceptable=ls_rchi2_acceptable, &
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
          generat_multiplier=generat_multiplier, &
          accwin_multiplier=accwin_multiplier)
     IF (error) THEN
        CALL errorMessage("oorb / lsl", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        CALL NULLIFY(epoch)
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     DEALLOCATE(HG_arr_storb_in, stat=err)
     ALLOCATE(HG_arr_storb_in(SIZE(obss_sep),1,4))
     HG_arr_storb_in = 99.9_bp
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
           DO j=1,SIZE(id_arr_storb_in,dim=1)
              IF (id_arr_storb_in(j) == id) THEN
                 IF (containsDiscretePDF(storb_arr_in(j))) THEN
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
        CALL NULLIFY(storb)
        CALL NEW(storb, obss_sep(i))        
        IF (error) THEN
           CALL errorMessage("oorb / lsl", &
                "TRACE BACK (75)", 1)
           STOP
        END IF
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             accept_multiplier=accwin_multiplier, &
             ls_correction_factor=ls_correction_factor, &
             ls_rchi2_acceptable=ls_rchi2_acceptable, &
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
           IF (info_verb >= 2) THEN
              WRITE(stdout,"(A,I0)") "Orbit #", j
           END IF
           CALL setParameters(orb_arr(j), &
                dyn_model=dyn_model_init, &
                perturbers=perturbers, &
                asteroid_perturbers=asteroid_perturbers, &
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
                asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, &
                integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / lsl", &
                   "TRACE BACK (95)", 1)
              STOP
           END IF
           CALL levenbergMarquardt(storb, orb_arr(j))
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
           DO k=1,SIZE(obs_masks,dim=1)
              IF (ALL(.NOT.obs_masks(k,:))) THEN
                 noutlier = noutlier + 1
              END IF
           END DO
           IF (noutlier > SIZE(obs_masks,dim=1)/5) THEN
              CALL warningMessage("oorb / lsl", &
                   "More than 20% of the observations have been" // &
                   "discarded as outliers.", 1)
              ! In case of too many outliers (>20% of obs), try new
              ! initial orbit or declare error.
              IF (j == SIZE(orb_arr,dim=1)) THEN
                 ! If no more initial orbits left, throw the
                 ! observation set in a separate bin and declare
                 ! error.
                 CALL errorMessage("oorb / lsl", &
                      "No more initial orbits left", 1)
                 CALL NEW(out_file, "observation_sets_with_many_outliers." // &
                      TRIM(observation_format_out))
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
                 DEALLOCATE(obs_masks, stat=err)
                 STOP
              ELSE
                 CALL warningMessage("oorb / lsl", &
                      "Will try another initial orbit.", 2)
                 DEALLOCATE(obs_masks, stat=err)
                 CYCLE
              END IF
           END IF

           IF (pp_H_estimation) THEN
              CALL NEW(physparam, storb)
              IF (pp_G > 99.0_bp) THEN
                 CALL estimateHAndG(physparam, obss_sep(i))
              ELSE
                 CALL estimateHAndG(physparam, obss_sep(i), &
                      input_G=pp_G, input_delta_G=pp_G_unc)
              END IF
              HG_arr_storb_in(i,1,1:2) = getH0(physparam)
              HG_arr_storb_in(i,1,3:4) = getG(physparam)
              CALL NULLIFY(physparam)
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
           SELECT CASE (TRIM(orbit_format_out))
           CASE ("des")
              CALL errorMessage("oorb / lsl", &
                   "DES format not yet supported for LSL output.", 1)
              STOP                 
           CASE ("orb")
              CALL writeOpenOrbOrbitFile(lu_orb_out, &
                   print_header=first, &
                   element_type_out=element_type_out_prm, &
                   id=id, &
                   orb=orb, &
                   cov=cov, &
                   H=HG_arr_storb_in(i,1,1), &
                   G=HG_arr_storb_in(i,1,3), &
                   mjd=mjd_epoch)
           CASE default
              CALL errorMessage("oorb / lsl", &
                   "Orbit format " // TRIM(orbit_format_out) // &
                   " not supported.",1)
              STOP           
           END SELECT
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
                 CALL toString(elements(j), element_str_arr(j), error, frmt="(E22.15)")
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (220)", 1)
                    STOP
                 END IF
                 stdev = SQRT(cov(j,j))
                 IF (element_type_comp_prm == "keplerian" .AND. j>=3) THEN
                    stdev = stdev/rad_deg
                 END IF
                 WRITE(getUnit(tmp_file), "(E22.15,1X)", &
                      advance="no") stdev
                 CALL toString(stdev, stdev_str_arr(j), error, frmt="(E22.15)")
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (225)", 1)
                    STOP
                 END IF
              END DO
              stdev = SQRT(cov(1,1))
              DO j=2,6
                 WRITE(getUnit(tmp_file), "(E22.15,1X)", &
                      advance="no") cov(1,j)/(stdev*SQRT(cov(j,j)))
                 CALL toString(cov(1,j)/(stdev*SQRT(cov(j,j))), &
                      corr_str_arr(j-1), error, frmt="(E22.15)")
                 IF (error) THEN
                    CALL errorMessage("oorb / lsl", &
                         "TRACE BACK (230)", 1)
                    STOP
                 END IF
              END DO
              WRITE(getUnit(tmp_file), "(E22.15,1X,A)") &
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
     END DO


  CASE ("covariance_sampling")

     ! Guided orbital inversion by randomly sampling an error
     ! hyperellipsoid defined by the covariance matrix of the global
     ! least-squares solution.

     CALL NULLIFY(epoch)
     CALL readConfigurationFile(conf_file, &
          t0=epoch, &
          dyn_model_init=dyn_model_init, &
          integrator_init=integrator_init, &
          integration_step_init=integration_step_init, &
          accwin_multiplier=accwin_multiplier, &
          cos_gaussian=cos_gaussian, &
          cos_nsigma=cos_nsigma, &
          cos_norb=cos_norb, &
          cos_ntrial=cos_ntrial)
     IF (error) THEN
        CALL errorMessage("oorb / covariance_sampling", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        CALL NULLIFY(epoch)
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     DEALLOCATE(HG_arr_storb_in, stat=err)
     ALLOCATE(HG_arr_storb_in(SIZE(obss_sep),cos_norb,4))
     HG_arr_storb_in = 99.9_bp
     ! Print header before printing first orbit:
     first = .TRUE.
     DO i=1,SIZE(obss_sep,dim=1)
        id = getID(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (25)", 1)
           STOP
        END IF
        IF (info_verb >= 2) THEN
           WRITE(stdout,"(1X,I0,3(A),1X,I0,A,I0)") i, &
                ". observation set (", TRIM(id), ")."
        END IF
        nobs = getNrOfObservations(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (30)", 1)
           STOP
        END IF
        IF (ALLOCATED(storb_arr_in)) THEN
           DO j=1,SIZE(id_arr_storb_in,dim=1)
              IF (id_arr_storb_in(j) == id) THEN
                 IF (containsDiscretePDF(storb_arr_in(j))) THEN
                    CALL errorMessage("oorb / covariance_sampling", &
                         "TRACE BACK (35)", 1)
                    STOP
                 ELSE
                    EXIT                    
                 END IF
              END IF
           END DO
           IF (j > SIZE(id_arr_storb_in,dim=1)) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "Initial orbit not available.", 1)
              STOP
           END IF
        ELSE
           CALL errorMessage("oorb / covariance_sampling", &
                "Initial orbit + covariance matrix not available.", 1)
           STOP
        END IF
        orb = getNominalOrbit(storb_arr_in(j))
        frame_ = getFrame(orb)
        cov = getCovarianceMatrix(storb_arr_in(j),element_type_comp_prm,frame_)
        CALL NEW(storb, orb, cov, element_type_comp_prm, &
             element_type_comp_prm, obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (75)", 1)
           STOP
        END IF
        dt = getObservationalTimespan(obss_sep(i))
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (40)", 1)
           STOP
        END IF
        IF (.NOT.exist(epoch)) THEN
           CALL NULLIFY(t)
           obs = getObservation(obss_sep(i),1)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (45)", 1)
              STOP
           END IF
           t = getTime(obs)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (50)", 1)
              STOP
           END IF
           CALL NULLIFY(obs)
           mjd = getMJD(t, "tt")
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (65)", 1)
              STOP
           END IF
           CALL NULLIFY(t)
           mjd = REAL(NINT(mjd+dt/2.0_bp),bp)
           CALL NEW(t, mjd, "tt")   
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (70)", 1)
              STOP
           END IF
        ELSE
           CALL NULLIFY(t)
           t = copy(epoch)
        END IF
        CALL setParameters(storb, &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step, &
             outlier_rejection=outlier_rejection_prm, &
             outlier_multiplier=outlier_multiplier_prm, &
             t_inv=t, &
             element_type=element_type_comp_prm, &
             accept_multiplier=accwin_multiplier, &
             cos_gaussian=cos_gaussian, &
             cos_nsigma=cos_nsigma, &
             cos_norb=cos_norb, &
             cos_ntrial=cos_ntrial)
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (80)", 1)
           STOP
        END IF
        CALL propagate(storb, t)
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "TRACE BACK (90)", 1)
           STOP
        END IF
        CALL covarianceSampling(storb)
        IF (error) THEN
           CALL errorMessage("oorb / covariance_sampling", &
                "Covariance sampling failed:", 1)
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
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (115)", 1)
              STOP
           END IF
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out))
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (120)", 1)
              STOP
           END IF
           CALL NULLIFY(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (125)", 1)
              STOP
           END IF

        ELSE

           IF (info_verb >= 2) THEN
              WRITE(stdout,"(3(1X,A))") "Covariance sampling for object", &
                   TRIM(id), "is ready."
           END IF

           IF (pp_H_estimation) THEN
              IF (info_verb >= 2) THEN
                 WRITE(stdout,"(1X,A)") "Computing HG parameters..."
              END IF
              CALL NEW(physparam, storb)
              IF (pp_G > 99.0_bp) THEN
                 CALL estimateHAndG(physparam, obss_sep(i))
              ELSE
                 CALL estimateHAndG(physparam, obss_sep(i), &
                      input_G=pp_G, input_delta_G=pp_G_unc)
              END IF
              HG_arr_in => getH0Distribution(physparam)
              HG_arr_storb_in(i,:,1:2) = HG_arr_in
              DEALLOCATE(HG_arr_in)
              HG_arr_in => getGDistribution(physparam)
              HG_arr_storb_in(i,:,3:4) = HG_arr_in
              DEALLOCATE(HG_arr_in)
              CALL NULLIFY(physparam)
              IF (info_verb >= 2) THEN
                 WRITE(stdout,"(1X,A)") "Computing HG parameters... done"
              END IF

           END IF

           CALL NEW(out_file, TRIM(id) // ".cos")
           CALL OPEN(out_file)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (150)", 1)
              STOP
           END IF
           ! WRITE OBSERVATIONS:
           CALL writeObservationFile(obss_sep(i), getUnit(out_file), &
                TRIM(observation_format_out)) 
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (155)", 1)
              STOP
           END IF
           WRITE(getUnit(out_file),"(A)") "#" 
!!$           ! WRITE COVARIANCE_SAMPLING PARAMETERS
!!$           CALL writeCOSResults(storb, obss_sep(i), getUnit(out_file))
!!$           IF (error) THEN
!!$              CALL errorMessage("oorb / covariance_sampling", &
!!$                   "TRACE BACK (160)", 1)
!!$              STOP
!!$           END IF
           CALL NULLIFY(out_file) 
!!$           CALL getResults(storb, sor_rho_cmp=sor_rho_cmp)
!!$           IF (error) THEN
!!$              CALL errorMessage("oorb / covariance_sampling", &
!!$                   "TRACE BACK (170)", 1)
!!$              STOP
!!$           END IF
           ! WRITE ORBITAL-ELEMENT PDF
           orb_arr_cmp => getSampleOrbits(storb)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling ", &
                   "TRACE BACK (175)", 1)
              STOP
           END IF
           pdf_arr_cmp => getDiscretePDF(storb)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling ", &
                   "TRACE BACK (180)", 1)
              STOP
           END IF
           rchi2_arr_cmp => getReducedChi2Distribution(storb)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling ", &
                   "TRACE BACK (185)", 1)
              STOP
           END IF
           CALL getResults(storb, &
                reg_apr_arr=reg_apr_arr_cmp, &
                jac_arr=jac_arr_cmp)
           IF (error) THEN
              CALL errorMessage("oorb / covariance_sampling", &
                   "TRACE BACK (190)", 1)
              STOP
           END IF
           IF (separately) THEN
              CALL NEW(out_file, TRIM(id) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling", &
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
                      H=HG_arr_storb_in(i,j,1), &
                      G=HG_arr_storb_in(i,j,3), &
                      mjd=mjd_epoch)
              ELSE
                 CALL errorMessage("oorb / covariance_sampling", &
                      TRIM(orbit_format_out) // &
                      " format not yet supported for Covariance_sampling output.", 1)
                 STOP                 
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling ", &
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
                 CALL errorMessage("oorb / covariance_sampling", &
                      "TRACE BACK (225)", 1)
                 STOP
              END IF
              CALL writeResiduals(storb, obss_sep(i), getUnit(out_file))
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling", &
                      "TRACE BACK (230)", 1)
                 STOP
              END IF
              CALL NULLIFY(out_file)
           END IF
           IF (plot_results) THEN
              CALL toString(dt, str, error, frmt="(F10.2)")
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling ", &
                      "TRACE BACK (235)", 1)
                 STOP
              END IF
              str = TRIM(id) // "_"// TRIM(str)
              CALL makeResidualStamps(storb, obss_sep(i), TRIM(str) // &
                   "_cos_residual_stamps.eps")
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling", &
                      "TRACE BACK (240)", 1)
                 STOP
              END IF
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_cos_residual_stamps.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_cos_residual_stamps.eps* &")
              END IF
              ALLOCATE(elements_arr(SIZE(orb_arr_cmp,dim=1),7), stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / covariance_sampling", &
                      "Could not allocate memory (3).", 1)
                 STOP
              END IF
              CALL NEW(tmp_file, TRIM(str)// "_cos_orbits.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling ", &
                      "TRACE BACK (250)", 1)
                 STOP
              END IF
              DO j=1,SIZE(orb_arr_cmp,dim=1)
                 IF (element_type_comp_prm == "cartesian") THEN
                    CALL rotateToEcliptic(orb_arr_cmp(j))
                 END IF
                 elements_arr(j,1:6) = getElements(orb_arr_cmp(j), element_type_comp_prm)
                 IF (error) THEN
                    CALL errorMessage("oorb / covariance_sampling", &
                         "TRACE BACK (255)", 1)
                    STOP
                 END IF
                 IF (element_type_comp_prm == "keplerian") THEN
                    elements_arr(j,3:6) = elements_arr(j,3:6)/rad_deg
                 END IF
                 elements_arr(j,7) = pdf_arr_cmp(j)
                 t = getTime(orb_arr_cmp(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / covariance_sampling ", &
                         "TRACE BACK (260)", 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmp_file),"(7(E23.15,1X),A)") &
                      elements_arr(j,1:6), &
                      pdf_arr_cmp(j), &
                      getCalendarDateString(t,"tdt")
                 IF (error) THEN
                    CALL errorMessage("oorb / covariance_sampling ", &
                         "TRACE BACK (265)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(t)
              END DO
              CALL NULLIFY(tmp_file)
              CALL NEW(tmp_file, TRIM(str) // &
                   "_cos_sample_standard_deviations.out")
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling ", &
                      "TRACE BACK (275)", 1)
                 STOP
              END IF
              WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                   advance="no") &
                   getObservationalTimespan(obss_sep(i))
              IF (error) THEN
                 CALL errorMessage("oorb / covariance_sampling ", &
                      "TRACE BACK (280)", 1)
                 STOP
              END IF
              DO j=1,6
                 CALL moments(elements_arr(:,j), &
                      pdf=pdf_arr_cmp, std_dev=stdev, errstr=errstr)
                 IF (LEN_TRIM(errstr) /= 0) THEN
                    CALL errorMessage("oorb / covariance_sampling", &
                         "Could not compute moments. " // TRIM(errstr), 1)
                    STOP
                 END IF
                 WRITE(getUnit(tmp_file), "(F22.15,1X)", &
                      advance="no") stdev
              END DO
              WRITE(getUnit(tmp_file),*)
              CALL NULLIFY(tmp_file)
              DEALLOCATE(elements_arr, stat=err)
              IF (err /= 0) THEN
                 CALL errorMessage("oorb / covariance_sampling", &
                      "Could not deallocate memory (5).", 1)
                 STOP
              END IF
              ! Make plot using gnuplot:
              CALL system("cp " // TRIM(str) // &
                   "_cos_orbits.out cos_orbits.out")
              IF (element_type_comp_prm == "cartesian") THEN
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/cos_plot_car.gp")
              ELSE
                 CALL system("gnuplot " // TRIM(gnuplot_scripts_dir) // "/cos_plot_kep.gp")
              END IF
              CALL system("cp cos_results.eps " // TRIM(str) // &
                   "_cos_" // TRIM(element_type_comp_prm) // &
                   "_results.eps")
              CALL system("rm -f cos_orbits.out cos_results.eps " // & 
                   TRIM(str) // "_cos_orbits.out " // TRIM(str) // &
                   "_cos_sample_standard_deviations.out")
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(str) // "_cos_" // &
                      TRIM(element_type_comp_prm) // "_results.eps")
              END IF
              IF (plot_open) THEN
                 CALL system("gv " // TRIM(str) // "_cos_" // &
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
              CALL errorMessage("oorb / covariance_sampling", &
                   "Could not deallocate memory (10).", 1)
              STOP
           END IF
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
        IF (ASSOCIATED(HG_arr_in)) THEN
           DEALLOCATE(HG_arr_in, stat=err)
        END IF
     END DO


  CASE ("residual-stamps")

     CALL readConfigurationFile(conf_file, &
          accwin_multiplier=accwin_multiplier)
     IF (error) THEN
        CALL errorMessage("oorb / residual-stamps", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     id = getID(obss_sep(1))
     dt = getObservationalTimespan(obss_sep(1))
     CALL toString(dt, str, error, frmt="(F10.2)")
     IF (error) THEN
        CALL errorMessage("oorb / residual-stamps", &
             "TRACE BACK (15)", 1)
        STOP
     END IF
     str = TRIM(id) // "_"// TRIM(str)
     CALL setParameters(storb_arr_in(1), &
          dyn_model=dyn_model, &
          perturbers=perturbers, &
          asteroid_perturbers=asteroid_perturbers, &
          integrator=integrator, &
          integration_step=integration_step, &
          outlier_rejection=outlier_rejection_prm, &
          outlier_multiplier=outlier_multiplier_prm, &
          accept_multiplier=accwin_multiplier)
     CALL makeResidualStamps(storb_arr_in(1), obss_sep(1), TRIM(str) // &
          "_residual_stamps.eps",  compute=.TRUE.)
     IF (error) THEN
        CALL errorMessage("oorb / residual-stamps", &
             "TRACE BACK (20)", 1)
        STOP
     END IF

  CASE ("propagation")

     first = .TRUE.

     CALL NULLIFY(epoch1)
     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        ! New epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
        CALL NEW(epoch1, mjd_tt, "TT")
        IF (error) THEN
           CALL errorMessage("oorb / propagation", &
                "TRACE BACK (5)", 1)
           STOP
        END IF
     ELSE IF (get_cl_option("--epoch-mjd-utc=", .FALSE.)) THEN
        mjd_utc = get_cl_option("--epoch-mjd-utc=", 0.0_bp)
        ! New epoch given as MJD UTC
        CALL NEW(epoch1, mjd_utc, "UTC")
        IF (error) THEN
           CALL errorMessage("oorb / propagation", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
     END IF

     IF (.NOT.exist(epoch1) .AND. &
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

           separately_ = separately

           IF (.NOT.exist(epoch1) .AND. &
                .NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.)) .AND. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) THEN
              ! New epoch relative to current epoch of orbit 
              dt = get_cl_option("--delta-epoch-mjd=", 0.0_bp)
              epoch1 = getTime(storb_arr_in(i))
              mjd_tt = getMJD(epoch1, "TT")
              CALL NULLIFY(epoch1)
              CALL NEW(epoch1, mjd_tt + dt, "TT")
           ELSE IF (.NOT.exist(epoch1) .AND. &
                .NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.) .OR. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) .AND. &
                ASSOCIATED(obss_sep)) THEN
              ! New epoch relative to observations...
              DO j=1,SIZE(obss_sep)
                 IF (id_arr_storb_in(i) == getID(obss_sep(j))) THEN
                    EXIT
                 END IF
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (20)", 1)
                    STOP
                 END IF
              END DO
              IF (get_cl_option("--last-observation-date", .FALSE.)) THEN
                 ! New epoch equal to last observation date
                 obs = getObservation(obss_sep(j),getNrOfObservations(obss_sep(j)))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (25)", 1)
                    STOP
                 END IF
                 epoch1 = getTime(obs)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (30)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(obs)                 
              ELSE
                 ! New epoch equal to observational mid-date
                 dt = getObservationalTimespan(obss_sep(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (35)", 1)
                    STOP
                 END IF
                 obs = getObservation(obss_sep(j),1)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (40)", 1)
                    STOP
                 END IF
                 epoch1 = getTime(obs)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (45)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(obs)
                 mjd = getMJD(epoch1, "tt")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (50)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(epoch1)
                 mjd = mjd + dt/2.0_bp
                 CALL NEW(epoch1, mjd, "tt")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (55)", 1)
                    STOP
                 END IF
              END IF
           ELSE IF (.NOT.exist(epoch1)) THEN
              CALL errorMessage("oorb / propagation", &
                   "No epoch specified (10).", 1)
              STOP              
           END IF
           epoch0 = getTime(storb_arr_in(i))
           mjd0 = getMJD(epoch0, "TT")
           mjd = mjd0
           mjd1 = getMJD(epoch1, "TT")

           ! How frequently should the elements be reported? Default is to
           ! only report them for the last epoch (output_interval < 0).
           output_interval = get_cl_option("--output-interval-days=", -1.0_bp)
           IF (output_interval < 0.0_bp) THEN
              interval = .FALSE.
           ELSE
              interval = .TRUE.
           END IF
           IF (interval .AND. output_interval >= integration_step) THEN
              output_interval = SIGN(output_interval,mjd1-mjd0)
           ELSE IF (interval .AND. output_interval < integration_step) THEN
              integration_step = output_interval
              output_interval = SIGN(output_interval,mjd1-mjd0)
           END IF

           ! Set integration parameters
           CALL setParameters(storb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (15)", 1)
              STOP
           END IF

           IF (separately) THEN
              CALL NEW(out_file, TRIM(id_arr_storb_in(i)) // ".orb")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (85)", 1)
                 STOP
              END IF
              lu_orb_out = getUnit(out_file)
           END IF

           ! Loop over the integration interval and output
           ! intermediate results during each step if requested
           integration_interval_storb: DO 
              IF ((interval .AND. output_interval < 0.0_bp .AND. &
                   mjd + output_interval <= mjd1) .OR. &
                   (interval .AND. output_interval > 0.0_bp .AND. &
                   mjd + output_interval >= mjd1) .OR. &
                   .NOT.interval) THEN
                 mjd = mjd1
              ELSE
                 mjd = mjd + output_interval
              END IF
              CALL NEW(epoch, mjd, "TT")

              ! Propagate orbital-element pdf from one epoch (=input) to another:
              IF (info_verb >= 2 .AND. dyn_model /= "2-body") THEN
                 CALL propagate(storb_arr_in(i), epoch, encounters=encounters)
              ELSE
                 CALL propagate(storb_arr_in(i), epoch)
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (60)", 1)
                 STOP
              END IF

              IF (info_verb >= 2 .AND. dyn_model /= "2-body") THEN
                 CALL getCalendarDate(epoch0, "TT", year0, month0, day0)
                 CALL getCalendarDate(epoch, "TT", year1, month1, day1)
                 WRITE(stderr,'(A)') ""
                 WRITE(stderr,'(A,I0,A,I0,A,F0.5,A,I0,A,I0,A,F0.5,A)') &
                      "Planetary encounters by object " // TRIM(id_arr_storb_in(i)) // &
                      " between ", year0, "-", month0, "-", day0, " and ", &
                      year1, "-", month1, "-", day1, ":"
                 DO j=1,SIZE(encounters,dim=1)
                    WRITE(stderr,'(A,I0,A)') "Orbit #", j, ":"
                    DO k=1,SIZE(encounters,dim=2)
                       IF (encounters(j,k,2) < 1.1_bp) THEN
                          WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F15.6,1X,A,1X,F0.6,1X,A)') &
                               "!! I-M-P-A-C-T !! with", planetary_locations(k), &
                               "at MJD", encounters(j,k,1), "TT given a stepsize of", &
                               encounters(j,k,4), "days."
                       ELSE
                          WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F15.6,1X,A,1X,F11.8,1X,A,1X,F0.6,1X,A)') &
                               "Closest encounter with", planetary_locations(k), &
                               "at MJD", encounters(j,k,1), &
                               "TT at a distance of", encounters(j,k,3), &
                               "AU given a stepsize of", encounters(j,k,4), "days."
                       END IF
                    END DO
                 END DO
              END IF

              IF (get_cl_option("--discard-impactors", .FALSE.)) THEN
                 DO j=1,SIZE(encounters,dim=2)
                    IF (encounters(1,j,2) < 1.1_bp) THEN 
                       DEALLOCATE(encounters, stat=err)
                       CALL NULLIFY(epoch0)
                       CALL NULLIFY(epoch1)
                       CALL NULLIFY(epoch)
                       EXIT integration_interval_storb
                    END IF
                 END DO
              END IF

              DEALLOCATE(encounters, stat=err)
              CALL NULLIFY(epoch0)
              IF (info_verb >= 2 .AND. dyn_model /= "2-body") THEN
                 epoch0 = copy(epoch)
              END IF
              CALL NULLIFY(epoch)

              IF (containsDiscretePDF(storb_arr_in(i))) THEN
                 ! Sampled orbital-element pdf:
                 orb_arr_cmp => getSampleOrbits(storb_arr_in(i))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation ", &
                         "TRACE BACK (65)", 1)
                    STOP
                 END IF
                 pdf_arr_cmp => getDiscretePDF(storb_arr_in(i))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation ", &
                         "TRACE BACK (70)", 1)
                    STOP
                 END IF
                 rchi2_arr_cmp => getReducedChi2Distribution(storb_arr_in(i))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation ", &
                         "TRACE BACK (75)", 1)
                    STOP
                 END IF
                 CALL getResults(storb_arr_in(i), &
                      reg_apr_arr=reg_apr_arr_cmp, &
                      jac_arr=jac_arr_cmp)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation ", &
                         "TRACE BACK (80)", 1)
                    STOP
                 END IF
                 DO j=1,SIZE(orb_arr_cmp,dim=1)
                    IF (orbit_format_out == "orb") THEN
                       CALL writeOpenOrbOrbitFile(lu_orb_out, &
                            print_header=first.OR.separately_, &
                            element_type_out=element_type_out_prm, &
                            id=id_arr_storb_in(i), &
                            orb=orb_arr_cmp(j), &
                            element_type_pdf=element_type_comp_prm, &
                            pdf=pdf_arr_cmp(j), &
                            rchi2=rchi2_arr_cmp(j), &
                            reg_apr=reg_apr_arr_cmp(j), &
                            jac_sph_inv=jac_arr_cmp(j,1), &
                            jac_car_kep=jac_arr_cmp(j,2), &
                            jac_equ_kep=jac_arr_cmp(j,3), &
                            H=HG_arr_storb_in(i,j,1), &
                            G=HG_arr_storb_in(i,j,3), &
                            mjd=mjd_epoch)
                    ELSE IF (orbit_format_out == "des") THEN
                       CALL errorMessage("oorb / propagation", &
                            "DES format not yet supported for propagation of sampled pdfs.", 1)
                       STOP                 
                    END IF
                    IF (error) THEN
                       CALL errorMessage("oorb / propagation ", &
                            "TRACE BACK (90)", 1)
                       STOP
                    END IF
                    CALL NULLIFY(orb_arr_cmp(j))
                    first = .FALSE.
                    separately_ = .FALSE.
                 END DO
                 DEALLOCATE(orb_arr_cmp, pdf_arr_cmp, rchi2_arr_cmp, &
                      reg_apr_arr_cmp, jac_arr_cmp)
              ELSE
                 ! LSL orbit:
                 orb = getNominalOrbit(storb_arr_in(i))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (95)", 1)
                    STOP
                 END IF
                 cov = getCovarianceMatrix(storb_arr_in(i), element_type_out_prm, "ecliptic")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (100)", 1)
                    STOP
                 END IF
                 IF (orbit_format_out == "orb") THEN
                    CALL writeOpenOrbOrbitFile(lu_orb_out, &
                         print_header=first.OR.separately_, &
                         element_type_out=element_type_out_prm, &
                         id=id_arr_storb_in(i), &
                         orb=orb, &
                         cov=cov, &
                         H=HG_arr_storb_in(i,1,1), &
                         G=HG_arr_storb_in(i,1,3), &
                         mjd=mjd_epoch)
                 ELSE IF (orbit_format_out == "des") THEN
                    CALL errorMessage("oorb / propagation", &
                         "DES format not yet supported for propagation of covariance matrices.", 1)
                    STOP                 
                 END IF
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (110)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(orb)
                 first = .FALSE.
                 separately_ = .FALSE.
              END IF

              IF ((interval .AND. ABS(mjd-mjd1) < EPSILON(mjd)) .OR. .NOT.interval) THEN
                 EXIT
              END IF

           END DO integration_interval_storb

           IF (separately) THEN
              CALL NULLIFY(out_file)
              IF (compress) THEN
                 CALL system("gzip -f " // TRIM(id_arr_storb_in(i)) // ".orb")
              END IF
           END IF

           IF (.NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.))) THEN
              CALL NULLIFY(epoch1)
           END IF

        END DO

     ELSE

        first = .TRUE.
        dt = get_cl_option("--delta-epoch-mjd=", 0.0_bp)

        IF (info_verb >= 2) THEN
           WRITE(stderr,"(A,1X,I0,1X,A)") "Integrating", &
                SIZE(orb_arr_in), "particles."
        END IF
        DO i=1,SIZE(orb_arr_in)

           IF (separately) THEN
              CALL NEW(out_file, TRIM(id_arr_in(i)) // ".des")
              CALL OPEN(out_file)
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (85)", 1)
                 STOP
              END IF
              lu_orb_out = getUnit(out_file)
           END IF

           IF (.NOT.exist(epoch1) .AND. &
                .NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.)) .AND. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) THEN
              ! New epoch relative to current epoch of orbit 
              dt = get_cl_option("--delta-epoch-mjd=", 0.0_bp)
              epoch1 = getTime(orb_arr_in(i))
              mjd_tt = getMJD(epoch1, "TT")
              CALL NULLIFY(epoch1)
              CALL NEW(epoch1, mjd_tt + dt, "TT")
           ELSE IF (.NOT.exist(epoch1) .AND. &
                .NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.) .OR. & 
                get_cl_option("--delta-epoch-mjd=", .FALSE.)) .AND. &
                ASSOCIATED(obss_sep)) THEN
              ! New epoch relative to observations...
              DO j=1,SIZE(obss_sep)
                 IF (id_arr_in(i) == getID(obss_sep(j))) THEN
                    EXIT
                 END IF
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (125)", 1)
                    STOP
                 END IF
              END DO
              IF (get_cl_option("--last-observation-date", .FALSE.)) THEN
                 ! New epoch equal to last observation date
                 obs = getObservation(obss_sep(j),getNrOfObservations(obss_sep(j)))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (130)", 1)
                    STOP
                 END IF
                 epoch1 = getTime(obs)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (135)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(obs)                 
              ELSE
                 ! New epoch equal to observational mid-date
                 dt = getObservationalTimespan(obss_sep(j))
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (140)", 1)
                    STOP
                 END IF
                 obs = getObservation(obss_sep(j),1)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (145)", 1)
                    STOP
                 END IF
                 epoch1 = getTime(obs)
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (150)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(obs)
                 mjd = getMJD(epoch1, "tt")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (155)", 1)
                    STOP
                 END IF
                 CALL NULLIFY(epoch1)
                 mjd = mjd + dt/2.0_bp
                 CALL NEW(epoch1, mjd, "tt")
                 IF (error) THEN
                    CALL errorMessage("oorb / propagation", &
                         "TRACE BACK (160)", 1)
                    STOP
                 END IF
              END IF
           ELSE IF (.NOT.exist(epoch1)) THEN
              CALL errorMessage("oorb / propagation", &
                   "No epoch specified (15).", 1)
              STOP              
           END IF
           epoch0 = getTime(orb_arr_in(i))
           mjd0 = getMJD(epoch0, "TT")
           mjd = mjd0
           mjd1 = getMJD(epoch1, "TT")

           ! How frequently should the elements be reported? Default is to
           ! only report them for the last epoch (output_interval < 0).
           output_interval = get_cl_option("--output-interval-days=", -1.0_bp)
           IF (output_interval < 0.0_bp) THEN
              interval = .FALSE.
           ELSE
              interval = .TRUE.
           END IF
           IF (interval .AND. output_interval >= integration_step) THEN
              output_interval = SIGN(output_interval,mjd1-mjd0)
           ELSE IF (interval .AND. output_interval < integration_step) THEN
              integration_step = output_interval
              output_interval = SIGN(output_interval,mjd1-mjd0)
           END IF

           ! Set integration parameters
           CALL setParameters(orb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / propagation", &
                   "TRACE BACK (120)", 1)
              STOP
           END IF

           ! Loop over the integration interval and output
           ! intermediate results during each step if requested
           integration_interval_orb: DO 
              IF ((interval .AND. output_interval < 0.0_bp .AND. &
                   mjd + output_interval <= mjd1) .OR. &
                   (interval .AND. output_interval > 0.0_bp .AND. &
                   mjd + output_interval >= mjd1) .OR. &
                   .NOT.interval) THEN
                 mjd = mjd1
              ELSE
                 mjd = mjd + output_interval
              END IF
              CALL NEW(epoch, mjd, "TT")
              IF ((info_verb >= 2 .OR. &
                   get_cl_option("--discard-impactors", .FALSE.)) .AND. &
                   dyn_model /= "2-body") THEN
                 CALL propagate(orb_arr_in(i), epoch, encounters=encounters)
              ELSE
                 CALL propagate(orb_arr_in(i), epoch)
              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "TRACE BACK (165)", 1)
                 STOP
              END IF
              IF (info_verb >= 2 .AND. dyn_model /= "2-body") THEN
                 CALL getCalendarDate(epoch0, "TT", year0, month0, day0)
                 CALL getCalendarDate(epoch, "TT", year1, month1, day1)
                 WRITE(stderr,'(A)') ""
                 WRITE(stderr,'(A,I0,A,I0,A,I0,A,F0.5,A,I0,A,I0,A,F0.5,A)') &
                      "Planetary encounters by object " // TRIM(id_arr_in(i)) // &
                      " (orbit #", i, ") between ", year0, "-", month0, &
                      "-", day0, " and ", year1, "-", month1, "-", day1, ":"
                 DO j=1,SIZE(encounters,dim=2)
                    IF (encounters(1,j,2) < 1.1_bp) THEN 
                       ! encounters(1,j,2) == 1 == impact
                       WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F15.6,1X,A,1X,F0.6,1X,A)') &
                            "!! I-M-P-A-C-T !! with", planetary_locations(j), &
                            "at MJD", encounters(1,j,1), "TT given a stepsize of", &
                            encounters(1,j,4), "days."
                    ELSE
                       WRITE(stderr,'(A22,1X,A7,1X,A6,1X,F15.6,1X,A,1X,F11.8,1X,A,1X,F0.6,1X,A)') &
                            "Closest encounter with", planetary_locations(j), &
                            "at MJD", encounters(1,j,1), "TT at a distance of", &
                            encounters(1,j,3), "AU given a stepsize of", encounters(1,j,4), "days."
                    END IF
                 END DO
              END IF
              IF (get_cl_option("--discard-impactors", .FALSE.)) THEN
                 DO j=1,SIZE(encounters,dim=2)
                    IF (encounters(1,j,2) < 1.1_bp) THEN 
                       DEALLOCATE(encounters, stat=err)
                       CALL NULLIFY(epoch0)
                       CALL NULLIFY(epoch1)
                       CALL NULLIFY(epoch)
                       EXIT integration_interval_orb
                    END IF
                 END DO
              END IF
              DEALLOCATE(encounters, stat=err)
              CALL NULLIFY(epoch0)
              IF (info_verb >= 2 .AND. dyn_model /= "2-body") THEN
                 epoch0 = copy(epoch)
              END IF
              CALL NULLIFY(epoch)

              IF (orbit_format_out == "des") THEN
                 CALL writeDESOrbitFile(lu_orb_out, & 
                      first, element_type_out_prm, &
                      id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), &
                      1, 6, -1.0_bp, "OPENORB", frame=frame, center=center)
              ELSE IF (orbit_format_out == "orb") THEN
                 CALL writeOpenOrbOrbitFile(lu_orb_out, &
                      print_header=first, &
                      element_type_out=element_type_out_prm, &
                      id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
                      H=HG_arr_in(i,1), G=HG_arr_in(i,3), &
                      mjd=mjd_epoch, frame=frame)
              END IF

!!$              IF (orbit_format_out == "des") THEN
!!$                 CALL writeDESOrbitFile(lu_orb_out, & 
!!$                      first, element_type_out_prm, &
!!$                      id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), &
!!$                      1, 6, -1.0_bp, "OPENORB")
!!$              ELSE IF (orbit_format_out == "orb") THEN
!!$                 CALL writeOpenOrbOrbitFile(lu_orb_out, &
!!$                      print_header=first, &
!!$                      element_type_out=element_type_out_prm, &
!!$                      id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
!!$                      H=HG_arr_in(i,1), G=HG_arr_in(i,3), &
!!$                      mjd=mjd_epoch)
!!$              END IF
              IF (error) THEN
                 CALL errorMessage("oorb / propagation", &
                      "Unable to propagate orbit for " // TRIM(id_arr_in(i)), 1)
                 STOP              
              END IF
              first = .FALSE.
              IF ((interval .AND. ABS(mjd-mjd1) < EPSILON(mjd)) .OR. .NOT.interval) THEN
                 EXIT
              END IF

           END DO integration_interval_orb

           IF (.NOT.(get_cl_option("--epoch-mjd-tt=", .FALSE.) .OR. &
                get_cl_option("--epoch-mjd-utc=", .FALSE.))) THEN
              CALL NULLIFY(epoch1)
           END IF

           IF (separately) THEN
              CALL NULLIFY(out_file)
           END IF

        END DO

     END IF


  CASE ("credible_region")

     probability_mass = get_cl_option("--probability_mass=", 0.9973_bp)

     IF (ALLOCATED(storb_arr_in)) THEN
        ! Input orbits contain uncertainty information.

        first = .TRUE.
        DO i=1,SIZE(storb_arr_in)

           orb_arr => getSampleOrbits(storb_arr_in(i), probability_mass)

           DO j=1,SIZE(orb_arr)

              IF (orbit_format_out == "des") THEN
                 CALL writeDESOrbitFile(lu_orb_out, & 
                      first, element_type_out_prm, &
                      id_arr_storb_in(i), orb_arr(j), HG_arr_storb_in(i,j,1), &
                      1, 6, -1.0_bp, "OPENORB", frame=frame, center=center)
              ELSE IF (orbit_format_out == "orb") THEN
                 CALL writeOpenOrbOrbitFile(lu_orb_out, &
                      print_header=first, &
                      element_type_out=element_type_out_prm, &
                      id=TRIM(id_arr_storb_in(i)), orb=orb_arr(j), &
                      H=HG_arr_storb_in(i,j,1), G=HG_arr_storb_in(i,j,3), &
                      mjd=mjd_epoch, frame=frame)
              END IF
              first = .FALSE.
              CALL NULLIFY(orb_arr(j))
           END DO

           DEALLOCATE(orb_arr)

        END DO

     ELSE

        CALL errorMessage("oorb / credible_region", &
             "Cannot compute credible region due to lack of uncertainty information", 1)
        STOP

     END IF

  CASE ("ephemeris")

     ! Input observatory code
     obsy_code = get_cl_option("--code=", obsy_code)

     ! Input start date [MJD in UTC]
     mjd_utc0 = get_cl_option("--start-mjd-utc=", -HUGE(mjd_utc0))

     ! Input evolutionary timespan [days]
     timespan = get_cl_option("--timespan=", 0.0_bp)

     ! Input time step [days]
     step = get_cl_option("--step=", integration_step)
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

     radians = get_cl_option("--radians", radians)

     G12 = get_cl_option("--G12=", -HUGE(G12))

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
                perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / ephemeris", &
                   "TRACE BACK (20)", 1)
              STOP
           END IF

           ! Compute topocentric ephemerides
           CALL getEphemerides(storb_arr_in(i), observers, ephemerides_arr, &
                cov_arr=cov_arr, pdf_arr=pdfs_arr)
           IF (error) THEN
              CALL errorMessage('oorb / ephemeris', &
                   'TRACE BACK (25)',1)
              STOP
           END IF

           ! Prepare output file
           IF (separately) THEN
              CALL NEW(tmp_file, TRIM(id_arr_storb_in(i)) // ".eph")
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
              IF (containsDiscretePDF(storb_arr_in(1)) .AND. &
                   .NOT.get_cl_option("--single-point-estimate", .FALSE.)) THEN
                 WRITE(lu,"(A,A11,5X,9(1X,A18))") "#", "Designation", &
                      "Observatory_code", "MJD_UTC", "Delta", "RA", &
                      "Dec", "dDelta/dt", "dRA/dt", "dDec/dt", "PDF"
              ELSE
                 WRITE(lu,"(A,A11,5X,30(1X,A18))") "#", "Designation", &
                      "Observatory_code ", "MJD_UTC", "Delta", "RA", "Dec", &
                      "dDelta/dt", "dRA/dt", "dDec/dt", "Delta_unc", &
                      "RA_unc", "Dec_unc", "dDelta/dt_unc", "dRA/dt_unc", &
                      "dDec/dt_unc", "cor(Delta,RA)", "cor(Delta,Dec)", &
                      "cor(Delta,dDelta)", "cor(Delta,dRA)", &
                      "cor(Delta,dDec)", "cor(RA,Dec)", &
                      "cor(RA,dDelta)", "cor(RA,dRA)", &
                      "cor(RA,dDec)", "cor(Dec,dDelta)", &
                      "cor(Dec,dRA)", "cor(Dec,dDec)", &
                      "cor(dDelta,dRA)", "cor(dDelta,dDec)", &
                      "cor(dRA,dDec)", "CurrentEphUnc"
              END IF
           END IF

           ! Loop over time steps:
           DO j=1,SIZE(observers)
              t = getTime(observers(j))
              mjd_utc = getMJD(t, "UTC")
              CALL NULLIFY(t)
              IF (containsDiscretePDF(storb_arr_in(i)) .AND. &
                   .NOT.get_cl_option("--single-point-estimate", .FALSE.)) THEN
                 ! Input orbits correspond to one or more sampled pdfs.
                 ! Loop over sampled orbits:
                 DO k=1,SIZE(ephemerides_arr,dim=1)
                    ! Make sure that the ephemeris is equatorial:
                    CALL rotateToEquatorial(ephemerides_arr(k,j))
                    coordinates = getCoordinates(ephemerides_arr(k,j))
                    WRITE(lu,"(2(A18,1X),7(F18.10,1X),E18.10)") &
                         TRIM(id_arr_storb_in(i)), TRIM(obsy_code_arr(j)), &
                         mjd_utc, coordinates(1), &
                         coordinates(2:3)/rad_deg, coordinates(4), &
                         coordinates(5:6)/rad_deg, pdfs_arr(k,j)
                    CALL NULLIFY(ephemerides_arr(k,j))
                 END DO
              ELSE
                 ! Input orbits correspond to one or more single-point
                 ! estimates of the orbital-element pdf, or the user
                 ! requests a single-point estimate for the ephemeris.
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
                 WRITE(lu,"(2(A18,1X),29(F18.10,1X))") &
                      TRIM(id_arr_storb_in(i)), &
                      TRIM(obsy_code_arr(j)), mjd_utc, coordinates(1), &
                      coordinates(2:3)/rad_deg, coordinates(4), &
                      coordinates(5:6)/rad_deg, stdev_arr(1), &
                      stdev_arr(2:3)/rad_deg, stdev_arr(4), &
                      stdev_arr(5:6)/rad_deg, corr(1,2:6), &
                      corr(2,3:6), corr(3,4:6), corr(4,5:6), corr(5,6), &
                      SQRT(2.3_bp*(cov_arr(2,2,j)*COS(coordinates(3))**2+cov_arr(3,3,j)))/rad_deg
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

        IF (ASSOCIATED(obss_sep)) THEN
           id_arr => getObjects(obss_in)
        END IF

        DO i=1,norb

           IF (ASSOCIATED(obss_sep)) THEN
              k = -1
              DO j=1,SIZE(id_arr)
                 IF (TRIM(id_arr_in(i)) == TRIM(id_arr(j))) THEN
                    k = j
                    EXIT
                 END IF
              END DO
              IF (k < 0) THEN
                 IF (err_verb >= 1) THEN
                    WRITE(stderr,*) "WARNING: Observations for object " // &
                         TRIM(id_arr_in(i)) // " not found..."
                 END IF
                 CYCLE
              END IF
              observers => getObservatoryCCoords(obss_sep(k))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (35)',1)
                 STOP
              END IF
              IF (info_verb >= 3) THEN

                 WRITE(stdout,"(18(1X,A23))") &
                      "GEO2OBSY_X_[AU]", &
                      "GEO2OBSY_Y_[AU]", &
                      "GEO2OBSY_Z_[AU]", &
                      "GEO2OBSY_dX/dt_[AU/d]", &
                      "GEO2OBSY_dY/dt_[AU/d]", &
                      "GEO2OBSY_dZ/dt_[AU/d]", &
                      "HELIO2GEO_X_[AU]", &
                      "HELIO2GEO_Y_[AU]", &
                      "HELIO2GEO_Z_[AU]", &
                      "HELIO2GEO_dX/dt_[AU/d]", &
                      "HELIO2GEO_dY/dt_[AU/d]", &
                      "HELIO2GEO_dZ/dt_[AU/d]", &
                      "HELIO2OBSY_X_[AU]", &
                      "HELIO2OBSY_Y_[AU]", &
                      "HELIO2OBSY_Z_[AU]", &
                      "HELIO2OBSY_dX/dt_[AU/d]", &
                      "HELIO2OBSY_dY/dt_[AU/d]", &
                      "HELIO2OBSY_dZ/dt_[AU/d]"

                 DO j=1,SIZE(observers)
                    !CALL rotateToEcliptic(observers(j))
                    CALL rotateToEquatorial(observers(j))
                    obsy_pos = getPosition(observers(j))
                    obsy_vel = getVelocity(observers(j))
                    t = getTime(observers(j))
                    ccoord = getObservatoryCCoord(obsies,"500",t)
                    !CALL rotateToEcliptic(ccoord)
                    CALL rotateToEquatorial(ccoord)
                    pos = getPosition(ccoord)
                    vel = getVelocity(ccoord)
                    geoc_obsy = obsy_pos - pos
                    dgeoc_obsy = obsy_vel - vel
                    WRITE(stdout,"(18(1X,F23.16))") &
                         geoc_obsy, dgeoc_obsy, &
                         pos, vel, obsy_pos, obsy_vel
                 END DO
              END IF
              DO j=1,SIZE(observers)
                 CALL rotateToEquatorial(observers(j))
              END DO
              obsy_code_arr => getObservatoryCodes(obss_sep(k))
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (40)',1)
                 STOP
              END IF
           ELSE
              IF (get_cl_option("--start-mjd-utc",.FALSE.)) THEN
                 CALL NEW(t, mjd_utc0, "UTC")
                 mjd_tt = getMJD(t, "TT")
              ELSE
                 t = getTime(orb_arr_in(i))
                 mjd_tt = getMJD(t, "TT")
              END IF
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
                 CALL NEW(t, mjd_tt+(j-1)*step, "TT")
                 IF (error) THEN
                    CALL errorMessage("oorb / ephemeris", &
                         "TRACE BACK (15)", 1)
                    STOP
                 END IF
              END DO
           END IF

           ! Set integration parameters
           CALL setParameters(orb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, asteroid_perturbers=asteroid_perturbers,&
                integrator=integrator, integration_step=integration_step)
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
              WRITE(lu,'(A,A11,1X,A,1X,38(A18,1X))') "#", &
                   "Designation", "Code", "MJD_UTC/UT1", "Delta", &
                   "RA", "Dec", "dDelta/dt", "dRA/dt", "dDec/dt", &
                   "VMag", "Alt", "PhaseAngle", "LunarElon", "LunarAlt", &
                   "LunarPhase", "SolarElon", "SolarAlt", "r", &
                   "HLon", "HLat", "TLon", "TLat", "TOCLon", &
                   "TOCLat", "HOCLon", "HOCLat", "TOppLon", &
                   "TOppLat", "HEclObj_X", "HEclObj_Y", "HEclObj_Z", &
                   "HEclObj_dX/dt", "HEclObj_dY/dt", & 
                   "HEclObj_dZ/dt", "HEclObsy_X", "HEclObsy_Y", &
                   "HEclObsy_Z", "EccAnom", "TrueAnom", "PosAngle"
           END IF

           DO j=1,SIZE(ephemerides)

              t = getTime(ephemerides(j))
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

              ! Compute position angle for direction of motion
              pa = ATAN2(dra,ddec)
              IF (pa < 0.0_bp) THEN
                 pa = two_pi + pa
              END IF

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

              ! Compute opposition-centered heliocentric ecliptic coordinates
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
                   observer_r2) / SQRT(heliocentric_r2 * ephemeris_r2)
              obj_phase = ACOS(cos_obj_phase)

              ! Compute apparent brightness:
              ! Input absolute magnitude
              IF (get_cl_option("--H=", .FALSE.)) THEN
                 H_value = get_cl_option("--H=", 99.9_bp)
              ELSE
                 H_value = HG_arr_in(i,1)
              END IF
              ! Input slope parameter
              IF (get_cl_option("--G12", .FALSE.)) THEN
                 obj_vmag = getApparentHG12Magnitude(H=H_value, &
                      G12=G12, r=SQRT(heliocentric_r2), &
                      Delta=Delta, phase_angle=obj_phase)
              ELSE
                 G_value = get_cl_option("--G=", HG_arr_in(i,3))
                 obj_vmag = getApparentHGMagnitude(H=H_value, &
                      G=G_value, r=SQRT(heliocentric_r2), &
                      Delta=Delta, phase_angle=obj_phase)
              END IF
              IF (error) THEN
                 CALL errorMessage('oorb / ephemeris', &
                      'TRACE BACK (70)',1)
                 STOP
              END IF

              IF (INDEX(obsy_code_arr(j),"-") == 0) THEN

                 ! Parameters relevant for an Earth-based observer (for
                 ! now, at least!):

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
                 ! Position of the Moon as seen from the Sun:
                 planeph => JPL_ephemeris(mjd_tt, 10, 11, error)
                 IF (error) THEN
                    CALL errorMessage('oorb / ephemeris', &
                         'TRACE BACK (95)',1)
                    STOP
                 END IF
                 sun_moon = planeph(1,1:3)
                 DEALLOCATE(planeph)
                 ! Angle between Sun and Moon as seen from the observatory:
                 obsy_moon = sun_moon - obsy_pos
                 obsy_moon_r2 = DOT_PRODUCT(obsy_moon,obsy_moon)
                 cos_obj_phase = DOT_PRODUCT(obsy_moon,-obsy_pos) / &
                      (SQRT(observer_r2) * SQRT(obsy_moon_r2))
                 lunar_phase = (1.0_bp-cos_obj_phase)/2.0_bp

                 ! Compute (approximate) distance between the target and the Moon:
                 vec3 = cross_product(obsy_obj,obsy_moon)
                 lunar_elongation = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_moon))

                 ! Compute (approximate) altitude of the Moon:
                 vec3 = cross_product(geoc_obsy,obsy_moon)
                 lunar_alt = pi/2.0_bp - ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(geoc_obsy,obsy_moon))

              ELSE

                 obj_alt = pi/2.0_bp
                 solar_elongation = -pi
                 lunar_phase = -1.0_bp
                 lunar_elongation = -pi
                 lunar_alt = pi/2.0_bp

              END IF

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

              elements = getElements(orb_arr_in(i), "cometary")
              IF (elements(2) < 1.0_bp) THEN
                 ! Compute eccentric and true anomalies if orbit is elliptic:
                 CALL solveKeplerEquation(orb_arr_in(i), orb_arr_in(i)%t, ecc_anom)
                 ta_s = SIN(ecc_anom)
                 ta_c = COS(ecc_anom)
                 fak = SQRT(1 - elements(2) * elements(2))
                 true_anom = MODULO(ATAN2(fak * ta_s, ta_c - elements(2)),two_pi)
              ELSE
                 ecc_anom = -99.0_bp
                 true_anom = -99.0_bp
              END IF

              IF (.NOT.radians) THEN 
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
                 ecc_anom = ecc_anom/rad_deg
                 true_anom = true_anom/rad_deg
                 pa = pa/rad_deg
              END IF

              WRITE(lu,'(2(A,1X),38(F18.10,1X))') &
                   id_arr_in(i), TRIM(obsy_code_arr(j)), mjd_utc, Delta, &
                   ra, dec, dDelta, dra, ddec, obj_vmag, obj_alt, &
                   obj_phase, lunar_elongation, lunar_alt, &
                   lunar_phase, solar_elongation, solar_alt, hdist, &
                   hlon, hlat, tlon, tlat, toclon, toclat, hoclon, &
                   hoclat, opplon, opplat, h_ecl_car_coord_obj, &
                   h_ecl_car_coord_obsy(1:3), ecc_anom, true_anom, &
                   pa

              CALL NULLIFY(observers(j))
              CALL NULLIFY(ephemerides(j))
              CALL NULLIFY(orb_lt_corr_arr(j))
              CALL NULLIFY(ccoord)
              CALL NULLIFY(obsy_ccoord)


           END DO

           IF (separately) THEN
              CALL NULLIFY(tmp_file)
           END IF
           DEALLOCATE(observers, ephemerides, orb_lt_corr_arr, &
                obsy_code_arr, stat=err)

        END DO

     END IF



  CASE ("phasecurve")

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
        CALL errorMessage('oorb / phasecurve', &
             'TRACE BACK (5)',1)
        STOP
     END IF


     IF (.NOT.ASSOCIATED(obss_sep)) THEN
        CALL errorMessage('oorb / phasecurve', &
             'Observations must be provided.',1)
        STOP
     END IF

     id_arr => getObjects(obss_in)

     DO i=1,norb

        k = -1
        DO j=1,SIZE(id_arr)
           IF (TRIM(id_arr_in(i)) == TRIM(id_arr(j))) THEN
              k = j
              EXIT
           END IF
        END DO
        IF (k < 0) THEN
           IF (err_verb >= 1) THEN
              WRITE(stderr,*) "WARNING: Observations for object " // &
                   TRIM(id_arr_in(i)) // " not found..."
           END IF
           CYCLE
        END IF
        observers => getObservatoryCCoords(obss_sep(k))
        IF (error) THEN
           CALL errorMessage('oorb / phasecurve', &
                'TRACE BACK (35)',1)
           STOP
        END IF
        nobs = SIZE(observers)
        DO j=1,nobs
           CALL rotateToEquatorial(observers(j))
        END DO
        obsy_code_arr => getObservatoryCodes(obss_sep(k))
        IF (error) THEN
           CALL errorMessage('oorb / phasecurve', &
                'TRACE BACK (40)',1)
           STOP
        END IF
        mags => getMagnitudes(obss_sep(k))
        IF (error) THEN
           CALL errorMessage('oorb / phasecurve', &
                'TRACE BACK (45)',1)
           STOP
        END IF
        filters => getFilters(obss_sep(k))
        IF (error) THEN
           CALL errorMessage('oorb / phasecurve', &
                'TRACE BACK (50)',1)
           STOP
        END IF

        ! Set integration parameters
        CALL setParameters(orb_arr_in(i), dyn_model=dyn_model, &
             perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, integration_step=integration_step)
        IF (error) THEN
           CALL errorMessage("oorb / phasecurve", &
                "TRACE BACK (15)", 1)
           STOP
        END IF

        ! Compute topocentric ephemerides
        CALL getEphemerides(orb_arr_in(i), observers, ephemerides, &
             this_lt_corr_arr=orb_lt_corr_arr)
        IF (error) THEN
           CALL errorMessage('oorb / phasecurve', &
                'TRACE BACK (20)',1)
           STOP
        END IF

        IF (separately) THEN
           CALL NEW(tmp_file, TRIM(id_arr_in(i)) // ".pc")
           IF (error) THEN
              CALL errorMessage('oorb / phasecurve', &
                   'TRACE BACK (25)',1)
              STOP
           END IF
           CALL OPEN(tmp_file)
           IF (error) THEN
              CALL errorMessage('oorb / phasecurve', &
                   'TRACE BACK (30)',1)
              STOP
           END IF
           lu = getUnit(tmp_file)
        ELSE
           lu = stdout
        END IF

        ! Calculate number of observations for which no reported
        ! magnitudes are available:
        k = 0
        DO j=1,nobs
           IF (mags(j) > 90.0_bp) THEN
              k = k + 1
           END IF
        END DO

        ! Header for the output (file)
        WRITE(lu,'("#",A,1X,I0)') TRIM(id_arr_in(i)), nobs-k

        DO j=1,nobs

           IF (mags(j) > 90.0_bp) THEN
              ! No reported magnitudes available
              CYCLE
           END IF

           ! Compute phase angle
           CALL NEW(ccoord, ephemerides(j))
           IF (error) THEN
              CALL errorMessage('oorb / phasecurve', &
                   'TRACE BACK (50)',1)
              STOP
           END IF
           CALL rotateToEquatorial(ccoord)
           obsy_obj = getPosition(ccoord)
           IF (error) THEN
              CALL errorMessage('oorb / phasecurve', &
                   'TRACE BACK (55)',1)
              STOP
           END IF
           ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
           CALL toCartesian(orb_lt_corr_arr(j), frame='equatorial')
           pos = getPosition(orb_lt_corr_arr(j))
           IF (error) THEN
              CALL errorMessage('oorb / phasecurve', &
                   'TRACE BACK (60)',1)
              STOP
           END IF
           heliocentric_r2 = DOT_PRODUCT(pos,pos)
           obsy_pos = getPosition(observers(j))
           IF (error) THEN
              CALL errorMessage('oorb / phasecurve', &
                   'TRACE BACK (65)',1)
              STOP
           END IF
           observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)
           cos_obj_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
                observer_r2) / SQRT(heliocentric_r2 * ephemeris_r2)
           obj_phase = ACOS(cos_obj_phase)

           ! Compute reduced magnitude:
           mag = mags(j) - 5.0_bp*LOG10(SQRT(heliocentric_r2*ephemeris_r2))

           ! Filter transformations to V:
           IF (obsy_code_arr(j) == "F51") THEN

              ! PS1SC internal analysis by A. Fitzsimmons 08-11-2011.
              ! Average results for (S+C).
              SELECT CASE (filters(j))
              CASE ("g")
                 obj_vmag = mag - 0.28_bp
              CASE ("r")
                 obj_vmag = mag + 0.23_bp
              CASE ("i")
                 obj_vmag = mag + 0.39_bp
              CASE ("z")
                 obj_vmag = mag + 0.37_bp
              CASE ("y")
                 obj_vmag = mag + 0.36_bp
              CASE ("w")
                 obj_vmag = mag + 0.16_bp
              CASE default
                 IF (error) THEN
                    CALL errorMessage('oorb / phasecurve', &
                         'Unknown filter for F51: ' // TRIM(filters(j)),1)
                    STOP
                 END IF
              END SELECT

           ELSE

              ! MPC transformations from G. Williams 05-11-2013.
              SELECT CASE (filters(j))
              CASE ("")
                 ! Values are stored as V-"Band" 
                 obj_vmag = mag - 0.8_bp
              CASE ("U")
                 obj_vmag = mag - 1.3_bp
              CASE ("B")
                 obj_vmag = mag - 0.8_bp
              CASE ("g")
                 ! [V-g = -0.60(B-V)+0.12, Jester et al. ,2005]
                 obj_vmag = mag - 0.35_bp
              CASE ("V")
                 obj_vmag = mag
              CASE ("r")
                 ! [V-r = 0.42(B-V)-0.11, Jester et al. ,2005]
                 obj_vmag = mag + 0.14_bp
              CASE ("R")
                 obj_vmag = mag + 0.4_bp
              CASE ("C")
                 obj_vmag = mag + 0.4_bp
              CASE ("W")
                 obj_vmag = mag + 0.4_bp
              CASE ("i")
                 ! [V-i = 0.91(R-I)+0.42(B-V)-0.31, Jester et al. ,2005]
                 obj_vmag = mag + 0.32_bp
              CASE ("z")
                 ! [V-z = 1.72(R-I)+0.42(B-V)-0.52, Jester et al. ,2005]
                 obj_vmag = mag + 0.26_bp
              CASE ("I")
                 obj_vmag = mag + 0.8_bp
              CASE ("J")
                 ! Was +1.4 
                 obj_vmag = mag + 1.2_bp
              CASE ("w")
                 obj_vmag = mag - 0.13_bp
              CASE ("y")
                 obj_vmag = mag + 0.32_bp
              CASE ("L")
                 obj_vmag = mag + 0.2_bp
              CASE default
                 IF (error) THEN
                    CALL errorMessage('oorb / phasecurve', &
                         'Unknown filter for ' // TRIM(obsy_code_arr(j)) // &
                         ': ' // TRIM(filters(j)),1)
                    STOP
                 END IF
              END SELECT

           END IF
!!$
!!$           obj_vmag = mag

           obj_phase = obj_phase/rad_deg
           t = getTime(ephemerides(j))
           err_verb_ = err_verb
           err_verb = 0
           mjd_utc = getMJD(t,"UTC") 
           IF (error) THEN
              error = .FALSE.
              mjd_utc = getMJD(t,"UT1") 
           END IF
           err_verb = err_verb_

           WRITE(lu,'(2(F13.8,1X),A,1X,F12.5)') obj_phase, obj_vmag, &
                TRIM(obsy_code_arr(j)), mjd_utc

           CALL NULLIFY(t)
           CALL NULLIFY(observers(j))
           CALL NULLIFY(ephemerides(j))
           CALL NULLIFY(orb_lt_corr_arr(j))
           CALL NULLIFY(ccoord)
           CALL NULLIFY(obsy_ccoord)


        END DO

        IF (separately) THEN
           CALL NULLIFY(tmp_file)
        END IF
        DEALLOCATE(observers, ephemerides, orb_lt_corr_arr, &
             obsy_code_arr, mags, filters, stat=err)

     END DO

  CASE ("synthetic_phasecurve")

     H = get_cl_option("--H=", 0.0_bp)
     G12 = get_cl_option("--G12=", 0.0_bp)

     DO i=0,150*5

        WRITE(stdout,*) i/5.0_bp, getApparentHG12Magnitude(H=H, &
             G12=G12, r=1.0_bp, Delta=1.0_bp, &
             phase_angle=i/5.0_bp*rad_deg)

     END DO


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
        IF (.NOT.containsDiscretePDF(storb_arr_in(i))) THEN
           CALL errorMessage("oorb / classification", &
                "Input orbits do not contain sampled " // &
                "uncertainty information required for this task.", 1)        
           STOP
        END IF
        ! Propagate input orbital-element pdf to cometary pdf:
        !CALL toCometary(storb_arr_in(i))
        ! Compute weights for each class that has been defined:
        CALL getGroupWeights(storb_arr_in(i), weight_arr, group_name_arr)
        ! Write output
        DO j=1,SIZE(weight_arr)
           WRITE(stdout,"(3X,A16,1X,A18,1X,F12.7)") TRIM(id_arr_storb_in(i)), &
                group_name_arr(j), weight_arr(j)
        END DO
        WRITE(stdout,*)
        DEALLOCATE(group_name_arr, weight_arr)
        CALL NULLIFY(storb_arr_in(i))
     END DO
     DEALLOCATE(storb_arr_in, id_arr_storb_in)

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
        IF (.NOT.containsDiscretePDF(storb_arr_in(i))) THEN
           CALL errorMessage("oorb / classification", &
                "Input orbits do not contain sampled " // &
                "uncertainty information required for this task.", 1)        
           STOP
        END IF
        ! Propagate input orbital-element pdf to cometary pdf:
        !CALL toCometary(storb_arr_in(i))
        ! Compute weights for each class that has been defined:
        CALL getGroupWeights(storb_arr_in(i), weight_arr, group_name_arr)
        weight_arr(5) = weight_arr(5)*35.0_bp
        !weight_arr(6:SIZE(weight_arr)-1) = weight_arr(6:SIZE(weight_arr)-1)*400.0_bp
        weight_arr(6:10) = weight_arr(6:10)*500.0_bp
        weight_arr = weight_arr/SUM(weight_arr)
        ! Write output
        DO j=1,SIZE(weight_arr)
           WRITE(stdout,"(3X,A16,1X,A18,1X,F12.7)") TRIM(id_arr_storb_in(i)), &
                group_name_arr(j), weight_arr(j)
        END DO
        WRITE(stdout,*)
        DEALLOCATE(group_name_arr, weight_arr)
        CALL NULLIFY(storb_arr_in(i))
     END DO
     DEALLOCATE(storb_arr_in, id_arr_storb_in)

  CASE ("ephemeris_planets")

     IF (get_cl_option("--epoch-mjd-tt=", .FALSE.)) THEN
        ! Epoch given as MJD TT
        mjd_tt = get_cl_option("--epoch-mjd-tt=", 0.0_bp)
     ELSE IF (get_cl_option("--epoch-mjd-utc=", .FALSE.)) THEN
        mjd_utc = get_cl_option("--epoch-mjd-utc=", 0.0_bp)
        ! Epoch given as MJD UTC
        CALL NEW(epoch, mjd_utc, "UTC")
        IF (error) THEN
           CALL errorMessage("oorb / propagation", &
                "TRACE BACK (10)", 1)
           STOP
        END IF
        mjd_tt = getMJD(epoch, "TT")
        CALL NULLIFY(epoch)
     END IF
     CALL NEW(epoch, mjd_tt, "TT")

     planeph => JPL_ephemeris(mjd_tt, -10, 11, error)
     coordinates = 0.0_bp
     elements = 0.0_bp
     i = 11
     WRITE(stdout,"(A,17(1X,E25.18))") &
          TRIM(planetary_locations(i)), 2400000.5_bp+mjd_tt, &
          coordinates, elements(1:2), elements(3:6)/rad_deg, &
          planetary_masses(i), &
          planetary_densities(i)/kgm3_smau3/1000.0_bp, &
          planetary_mu(i), 0.0_bp     
     DO i=1,SIZE(planeph,dim=1)
        CALL NEW(ccoord, planeph(i,:), "equatorial", epoch)
        CALL rotateToEcliptic(ccoord)
        coordinates = getCoordinates(ccoord)
        CALL NEW(orb,ccoord)
        elements = getElements(orb, "keplerian")
        CALL NULLIFY(ccoord)
        CALL NULLIFY(orb)
        WRITE(stdout,"(A,17(1X,E25.18))") &
             TRIM(planetary_locations(i)), 2400000.5_bp+mjd_tt, &
             coordinates, elements(1:2), elements(3:6)/rad_deg, &
             planetary_masses(i), &
             planetary_densities(i)/kgm3_smau3/1000.0_bp, &
             planetary_mu(i), Hill_radius(planetary_masses(i), &
             planetary_masses(11), elements(1), elements(2))
     END DO
     DEALLOCATE(planeph)

  CASE ("tisserand_parameters")

     ! Compute Tisserand's parameter(s) for input orbit(s)
     IF (ALLOCATED(storb_arr_in)) THEN
        WRITE(stdout,*) "--task=tisserand_parameters is not yet implemented for orbit PDFs."
     ELSE
        DO i=1,SIZE(orb_arr_in)
           tisserands_parameters = getTisserandsParameters(orb_arr_in(i))
           WRITE(stdout,*) TRIM(id_arr_in(i)), tisserands_parameters 
        END DO
     END IF


  CASE ("jacobi_constants")

     ! Compute Jacobi constants for input orbit(s)
     IF (ALLOCATED(storb_arr_in)) THEN
        WRITE(stdout,*) "--task=jacobi_constants is not yet implemented for orbit PDFs."
     ELSE
        DO i=1,SIZE(orb_arr_in)
           jacobi_constants = getJacobiConstants(orb_arr_in(i))
           WRITE(stdout,"(A,10(1X,F10.5))") TRIM(id_arr_in(i)), jacobi_constants
        END DO
     END IF


  CASE ("apoapsis_distance")

     ! Compute apoapsis distance(s) for input orbit(s)
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           CALL getApoapsisDistance(storb_arr_in(i), apoapsis_distance_pdf)
           IF (containsDiscretePDF(storb_arr_in(i))) THEN
              DO j=1,SIZE(apoapsis_distance_pdf,dim=1)
                 WRITE(stdout,*) apoapsis_distance_pdf(j,:) 
              END DO
           ELSE
              WRITE(stdout,*) apoapsis_distance_pdf(1,:) 
           END IF
           DEALLOCATE(apoapsis_distance_pdf)
        END DO
     ELSE
        DO i=1,SIZE(orb_arr_in)
           CALL getApoapsisDistance(orb_arr_in(i), apoapsis_distance)
           WRITE(stdout,*) apoapsis_distance 
        END DO
     END IF


  CASE ("periapsis_distance")

     ! Compute periapsis distance(s) for input orbit(s)
     IF (ALLOCATED(storb_arr_in)) THEN
        DO i=1,SIZE(storb_arr_in)
           CALL getPeriapsisDistance(storb_arr_in(i), periapsis_distance_pdf)
           IF (containsDiscretePDF(storb_arr_in(i))) THEN
              DO j=1,SIZE(periapsis_distance_pdf,dim=1)
                 WRITE(stdout,*) periapsis_distance_pdf(j,:) 
              END DO
           ELSE
              WRITE(stdout,*) periapsis_distance_pdf(1,:) 
           END IF
           DEALLOCATE(periapsis_distance_pdf)
        END DO
     ELSE
        DO i=1,SIZE(orb_arr_in)
           CALL getPeriapsisDistance(orb_arr_in(i), periapsis_distance)
           WRITE(stdout,*) periapsis_distance 
        END DO
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
           CALL writeDESOrbitFile(lu_orb_out, i==1, element_type_out_prm, &
                id_arr_in(i), orb_arr_in(i), HG_arr_in(i,1), 1, 6, &
                moid, "OPENORB", frame=frame, center=center)
        ELSE IF (orbit_format_out == "orb") THEN
           CALL writeOpenOrbOrbitFile(lu_orb_out, print_header=i==1, &
                element_type_out=element_type_out_prm, &
                id=TRIM(id_arr_in(i)), orb=orb_arr_in(i), &
                H=HG_arr_in(i,1), G=HG_arr_in(i,3), &
                mjd=mjd_epoch, frame=frame)
        END IF
        IF (error) THEN
           WRITE(stderr,*) i
           error = .FALSE.
        END IF
     END DO
     CALL NULLIFY(orb_out_file)



  CASE ("obsplanner")

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
     dt_fulfill_night = get_cl_option("--timespan-fulfill-night=", 0.05_bp) ! 0.05
     minstep = MAX(1,CEILING(dt_fulfill_night/step))
     ALLOCATE(temp_arr(minstep,37))

     ! Input minimum altitude for object (wrt the horizon)
     obj_alt_min = get_cl_option("--object-altitude-min=", -90.0_bp) ! 30
     obj_alt_min = obj_alt_min*rad_deg

     ! Input maximum apparent magnitude (minimum brightness) for object
     obj_vmag_max = get_cl_option("--object-magnitude-max=", HUGE(obj_vmag_max))

     ! Input maximum altitude for Sun
     solar_alt_max = get_cl_option("--solar-altitude-max=", 90.0_bp) ! -18
     solar_alt_max = solar_alt_max*rad_deg

     ! Input maximum solar elongation
     solar_elon_max = get_cl_option("--solar-elongation-max=", 180.0_bp) ! 180
     solar_elon_max = solar_elon_max*rad_deg

     ! Input minimum solar elongation
     solar_elon_min = get_cl_option("--solar-elongation-min=", 0.0_bp) ! 45
     solar_elon_min = solar_elon_min*rad_deg

     ! Input maximum altitude for Moon
     lunar_alt_max = get_cl_option("--lunar-altitude-max=", 90.0_bp) ! 90
     lunar_alt_max = lunar_alt_max*rad_deg

     ! Input minimum distance to Moon
     lunar_elongation_min = get_cl_option("--lunar-elongation-min=", 0.0_bp) ! 45
     lunar_elongation_min = lunar_elongation_min*rad_deg

     ! Input minimum phase of Moon
     lunar_phase_min = get_cl_option("--lunar-phase-min=", 0.0_bp) ! 0

     ! Input maximum phase of Moon
     lunar_phase_max = get_cl_option("--lunar-phase-max=", 1.0_bp) ! 1

     IF (.NOT.get_cl_option("--cometary-magnitude", .FALSE.)) THEN
        obj_vmag_cometary = -99.0_bp
     END IF

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
                perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, integration_step=integration_step)
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

           ! Size or absolute magnitude of object
           H_value = 99.9_bp
           H10_value = 99.9_bp
           b1 = 99.9_bp
           b2 = 99.9_bp
           IF (get_cl_option("--diameter-in-km=", .FALSE.)) THEN
              diameter = get_cl_option("--diameter-in-km=", 0.0_bp)
              ! Input albedo
              geometric_albedo = get_cl_option("--geometric-albedo=", -1.0_bp)
              IF (geometric_albedo < 0.0_bp) THEN
                 CALL errorMessage("oorb / obsplanner", &
                      "Providing the geometric albedo is " // &
                      "mandatory when providing diameter " // &
                      "of object.", 1)
                 STOP                 
              END IF
              H_value = (LOG10(diameter) + &
                   0.5_bp*LOG10(geometric_albedo) - &
                   3.1236_bp) / (-0.2_bp)
              IF (get_cl_option("--cometary-magnitude", .FALSE.)) THEN
                 b1 = get_cl_option("--b1=", b1)
                 b2 = get_cl_option("--b2=", b2)
                 IF (b1 > 99.0_bp .OR. b2 > 99.0_bp) THEN
                    CALL errorMessage("oorb / obsplanner", &
                         "Providing b1 and b2 is " // &
                         "mandatory when requesting" // &
                         "cometary magnitudes.", 1)
                    STOP                 
                 END IF
                 H10_value = (LOG10(diameter) - b2) / b1
              END IF
           ELSE IF (get_cl_option("--H=", .FALSE.)) THEN
              H_value = get_cl_option("--H=", 99.9_bp)
           ELSE
              H_value = HG_arr_in(i,1)
           END IF

           ! Input slope parameter
           G_value = get_cl_option("--G=", HG_arr_in(i,3))

           IF (separately .OR. i == 1) THEN
              WRITE(lu,'(A3,1X,A11,1X,A,1X,43(A18,1X))') "#RN", &
                   "Designation", "Code", "MJD_UTC", "Delta", &
                   "RA", "Dec", "dDelta/dt", "dRA/dt", "dDec/dt", &
                   "VMag", "Alt", "PhaseAngle", "LunarElon", "LunarAlt", &
                   "LunarPhase", "SolarElon", "SolarAlt", "r", &
                   "HLon", "HLat", "TLon", "TLat", "TOCLon", &
                   "TOCLat", "HOCLon", "HOCLat", "TOppLon", &
                   "TOppLat", "HEclObj_X", "HEclObj_Y", "HEclObj_Z", &
                   "HEclObj_dX/dt", "HEclObj_dY/dt", & 
                   "HEclObj_dZ/dt", "HEclObsy_X", "HEclObsy_Y", &
                   "HEclObsy_Z", "PosAngle", "CalendarDate_UTC", &
                   "H", "G", "Diameter", "p_V", "VMag_cometary", "H10"
           END IF

           istep = 0
           DO j=1,nstep

              istep = istep + 1
              t = getTime(observers(j))
              mjd_tt = getMJD(t, "TT")
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (29)',1)
                 STOP
              END IF
              mjd_utc = getMJD(t, "UTC")
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (30)',1)
                 STOP
              END IF
              CALL getCalendarDate(t, "UTC", caldate)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (31)',1)
                 STOP
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

              ! Compute position angle for direction of motion
              pa = ATAN2(dra,ddec)
              IF (pa < 0.0_bp) THEN
                 pa = two_pi + pa
              END IF

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
              obj_vmag = getApparentHGMagnitude(H=H_value, &
                   G=G_value, r=SQRT(heliocentric_r2), &
                   Delta=Delta, phase_angle=obj_phase)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (70)',1)
                 STOP
              END IF

              ! Add effect of cometary brightening in the inner solar
              ! system (see Cook et al. 2012, in prep).
              IF (get_cl_option("--cometary-magnitude", .FALSE.)) THEN
                 elements = getElements(orb_lt_corr_arr(j), "cometary")
                 IF (error) THEN
                    CALL errorMessage('oorb / obsplanner', &
                         'TRACE BACK (70)',1)
                    STOP
                 END IF
                 IF (mjd_tt < elements(6)) THEN
                    n = 5.0_bp
                 ELSE
                    n = 3.5_bp
                 END IF
                 obj_vmag_cometary = getApparentCometaryMagnitude(&
                      H10=H10_value, G=G_value, &
                      r=SQRT(heliocentric_r2), Delta=Delta, &
                      phase_angle=obj_phase, n=n)
                 IF (obj_vmag_cometary < obj_vmag) THEN
                    obj_vmag = obj_vmag_cometary
                 END IF
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
              ! Position of the Moon as seen from the Sun:
              planeph => JPL_ephemeris(mjd_tt, 10, 11, error)
              IF (error) THEN
                 CALL errorMessage('oorb / obsplanner', &
                      'TRACE BACK (87)',1)
                 STOP
              END IF
              sun_moon = planeph(1,1:3)
              DEALLOCATE(planeph)
              ! Angle between Sun and Moon as seen from the observatory:
              obsy_moon = sun_moon - obsy_pos
              obsy_moon_r2 = DOT_PRODUCT(obsy_moon,obsy_moon)
              cos_obj_phase = DOT_PRODUCT(obsy_moon,-obsy_pos) / &
                   (SQRT(observer_r2) * SQRT(obsy_moon_r2))
              lunar_phase = (1.0_bp-cos_obj_phase)/2.0_bp
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
              pa = pa/rad_deg

              IF (istep <= minstep) THEN
                 temp_arr(istep,:) = &
                      (/ mjd_utc, Delta, &
                      ra, dec, dDelta, dra, ddec, obj_vmag, obj_alt, &
                      obj_phase, lunar_elongation, lunar_alt, &
                      lunar_phase, solar_elongation, solar_alt, hdist, &
                      hlon, hlat, tlon, tlat, toclon, toclat, hoclon, &
                      hoclat, opplon, opplat, h_ecl_car_coord_obj, &
                      h_ecl_car_coord_obsy(1:3), pa, caldate /)
              END IF

              IF (istep == minstep) THEN
                 DO k=1,istep
                    WRITE(lu,'(I0,1X,2(A,1X),36(F18.6,1X),F18.1,6(1X,F18.6))') &
                         i, id_arr_in(i), TRIM(obsy_code), &
                         temp_arr(k,:), H_value, G_value, &
                         diameter, geometric_albedo, &
                         obj_vmag_cometary, H10_value
                 END DO
              ELSE IF (istep > minstep) THEN
                 WRITE(lu,'(I0,1X,2(A,1X),36(F18.6,1X),F18.1,6(1X,F18.6))') &
                      i, id_arr_in(i), TRIM(obsy_code), mjd_utc, Delta, &
                      ra, dec, dDelta, dra, ddec, obj_vmag, obj_alt, &
                      obj_phase, lunar_elongation, lunar_alt, &
                      lunar_phase, solar_elongation, solar_alt, &
                      hdist, hlon, hlat, tlon, tlat, toclon, toclat, &
                      hoclon, hoclat, opplon, opplat, &
                      h_ecl_car_coord_obj, &
                      h_ecl_car_coord_obsy(1:3), pa, caldate, H_value, &
                      G_value, diameter, geometric_albedo, &
                      obj_vmag_cometary, H10_value
              END IF

           END DO

           IF (separately) THEN
              CALL NULLIFY(tmp_file)
           END IF
           DEALLOCATE(observers, ephemerides, orb_lt_corr_arr)

        END DO

     END IF

     DEALLOCATE(temp_arr)
     CALL NULLIFY(obsies)



  CASE ("fou")

     !! Produces a table with the quantities required for the
     !! computation of the Figure of Urgency proposed for scheduling
     !! follow-up observations with NEOSSat. The three key quantities
     !! are the minimum solar elongation, the minimum apparent
     !! brightness (max mag), and the 3-sigma (~equivalent) ephemeris
     !! uncertainty.

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
        CALL errorMessage('oorb / fou', &
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
                    CALL errorMessage("oorb / fou", &
                         "TRACE BACK (10)", 1)
                    STOP
                 END IF
                 obsy_code_arr(j) = obsy_code
                 ! Compute heliocentric observatory coordinates
                 observers(j) = getObservatoryCCoord(obsies, obsy_code_arr(j), t)
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (15)',1)
                    STOP
                 END IF
                 CALL rotateToEquatorial(observers(j))
                 CALL NULLIFY(t)
              END DO
           END IF

           ! Set integration parameters
           CALL setParameters(storb_arr_in(i), dyn_model=dyn_model, &
                perturbers=perturbers, asteroid_perturbers=asteroid_perturbers, &
                integrator=integrator, integration_step=integration_step)
           IF (error) THEN
              CALL errorMessage("oorb / fou", &
                   "TRACE BACK (20)", 1)
              STOP
           END IF

           ! Compute topocentric ephemerides
           CALL getEphemerides(storb_arr_in(i), observers, ephemerides_arr, &
                cov_arr=cov_arr, pdf_arr=pdfs_arr, this_lt_corr_arr=orb_lt_corr_arr2)
           IF (error) THEN
              CALL errorMessage('oorb / fou', &
                   'TRACE BACK (25)',1)
              STOP
           END IF

           IF (separately) THEN
              CALL NEW(tmp_file, TRIM(id_arr_storb_in(i)) // ".fou")
              IF (error) THEN
                 CALL errorMessage('oorb / fou', &
                      'TRACE BACK (30)',1)
                 STOP
              END IF
              CALL OPEN(tmp_file)
              IF (error) THEN
                 CALL errorMessage('oorb / fou', &
                      'TRACE BACK (35)',1)
                 STOP
              END IF
              lu = getUnit(tmp_file)
           ELSE
              lu = stdout
           END IF

           IF (separately .OR. i == 1) THEN
              WRITE(lu,'(9(A,2X))') "#JD", "SolElon_min", &
                   "SolElon_max", "VMag_max", "EELon", "ELat", &
                   "EphWidth", "RA", "Dec"
           END IF

           ! Loop over time steps:
           DO j=1,nstep

              t = getTime(observers(j))
              mjd_tt = getMJD(t, "TT")
              mjd_utc = getMJD(t, "UTC")

              obsy_pos = getPosition(observers(j))
              IF (error) THEN
                 CALL errorMessage('oorb / fou', &
                      'TRACE BACK (40)',1)
                 STOP
              END IF
              observer_r2 = DOT_PRODUCT(obsy_pos,obsy_pos)

              obsy_ccoord = getObservatoryCCoord(obsies, obsy_code_arr(j), t)
              IF (error) THEN
                 CALL errorMessage('oorb / fou', &
                      'TRACE BACK (45)',1)
                 STOP
              END IF
              CALL rotateToEcliptic(obsy_ccoord)
              sun_ccoord = opposite(obsy_ccoord)
              CALL NULLIFY(obsy_ccoord)
              obsy_sun = getPosition(sun_ccoord)

              ! Compute topocentric coordinates for the Sun (actually,
              ! get heliocentric observatory coordinates)
              scoord = getSCoord(sun_ccoord)
              CALL NULLIFY(sun_ccoord)
              CALL rotateToEcliptic(scoord)
              pos_sun = getPosition(scoord)
              sunlon = pos_sun(2)
              sunlat = pos_sun(3)
              CALL NULLIFY(scoord)

              IF (containsDiscretePDF(storb_arr_in(i))) THEN

                 ! Input orbits correspond to one or more sampled pdfs.
                 ! Loop over sampled orbits:
                 norb = SIZE(ephemerides_arr,dim=1)
                 ALLOCATE(temp_arr(norb,6))
                 DO k=1,norb

                    ! RA & DEC
                    CALL rotateToEquatorial(ephemerides_arr(k,j))        
                    comp_coord = getCoordinates(ephemerides_arr(k,j))
                    IF (error) THEN
                       CALL errorMessage('oorb / fou', &
                            'TRACE BACK (50)',1)
                       STOP
                    END IF
                    Delta = comp_coord(1)
                    temp_arr(k,5:6) = comp_coord(2:3)

                    ! Compute topocentric ecliptic coordinates
                    ! relative to the Sun
                    CALL rotateToEcliptic(ephemerides_arr(k,j))        
                    comp_coord = getCoordinates(ephemerides_arr(k,j))
                    IF (error) THEN
                       CALL errorMessage('oorb / fou', &
                            'TRACE BACK (55)',1)
                       STOP
                    END IF
                    tsclon = comp_coord(2) - sunlon
                    tsclat = comp_coord(3) - sunlat
                    IF (tsclon > pi) THEN
                       tsclon = tsclon - two_pi
                    ELSE IF (tsclon < -pi) THEN
                       tsclon = tsclon + two_pi
                    END IF
                    temp_arr(k,1:2) = (/ tsclon, tsclat /)

                    ! APPARENT BRIGHTNESS:
                    ! Compute phase angle
                    CALL NEW(ccoord, ephemerides_arr(k,j))
                    IF (error) THEN
                       CALL errorMessage('oorb / fou', &
                            'TRACE BACK (60)',1)
                       STOP
                    END IF
                    CALL rotateToEcliptic(ccoord)
                    obsy_obj = getPosition(ccoord)
                    IF (error) THEN
                       CALL errorMessage('oorb / fou', &
                            'TRACE BACK (65)',1)
                       STOP
                    END IF
                    ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
                    CALL toCartesian(orb_lt_corr_arr2(k,j), frame='ecliptic')
                    pos = getPosition(orb_lt_corr_arr2(k,j))
                    IF (error) THEN
                       CALL errorMessage('oorb / fou', &
                            'TRACE BACK (70)',1)
                       STOP
                    END IF
                    heliocentric_r2 = DOT_PRODUCT(pos,pos)
                    cos_obj_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
                         observer_r2) / (SQRT(heliocentric_r2) * &
                         SQRT(ephemeris_r2))
                    obj_phase = ACOS(cos_obj_phase)
                    ! Input absolute magnitude
                    H_value = get_cl_option("--H=", HG_arr_storb_in(i,k,1))
                    ! Input slope parameter
                    G_value = get_cl_option("--G=", HG_arr_storb_in(i,k,3))
                    temp_arr(k,3) = getApparentHGMagnitude(H=H_value, &
                         G=G_value, r=SQRT(heliocentric_r2), &
                         Delta=Delta, phase_angle=obj_phase)
                    IF (error) THEN
                       CALL errorMessage('oorb / fou', &
                            'TRACE BACK (75)',1)
                       STOP
                    END IF

                    ! SOLAR ELONGATION
                    vec3 = cross_product(obsy_obj,obsy_sun)
                    temp_arr(k,4) = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_sun))

                    CALL NULLIFY(ccoord)
                    CALL NULLIFY(ephemerides_arr(k,j))
                    CALL NULLIFY(orb_lt_corr_arr2(k,j))

                 END DO

                 CALL NULLIFY(observers(j))

                 ! WIDTH OF EPHEMERIS UNCERTAINTY
                 ! First, try in measuring largest width in RA:
                 nhist = 360
                 ALLOCATE(hist(nhist,2))
                 CALL histogram(temp_arr(:,5), hist, 0.0_bp, two_pi)
                 DO k=1,nhist
                    IF (hist(k,2) == 0.0_bp) THEN
                       EXIT
                    END IF
                 END DO
                 DEALLOCATE(hist)
                 ALLOCATE(ephem_(norb,2))
                 WHERE (temp_arr(1:norb,5) < two_pi*((k-1)/nhist))
                    ephem_(1:norb,1) = temp_arr(1:norb,5) + two_pi*(1 - (k-1)/nhist)
                    ephem_(1:norb,2) = temp_arr(1:norb,6)
                 END WHERE
                 WHERE (temp_arr(1:norb,5) > two_pi*((k-1)/nhist))
                    ephem_(1:norb,1) = temp_arr(1:norb,5) - two_pi*((k-1)/nhist)
                    ephem_(1:norb,2) = temp_arr(1:norb,6)
                 END WHERE
                 indx_min = MINLOC(ephem_(1:norb,1),dim=1)
                 indx_max = MAXLOC(ephem_(1:norb,1),dim=1)
                 angular_distance_ra = angularDistance(ephem_(indx_min,1), ephem_(indx_min,2), &
                      ephem_(indx_max,1), ephem_(indx_max,2))
                 DEALLOCATE(ephem_)
                 ! Second, try in measuring largest width in Dec:
                 indx_min = MINLOC(temp_arr(1:norb,6),dim=1)
                 indx_max = MAXLOC(temp_arr(1:norb,6),dim=1)
                 angular_distance_dec = angularDistance(temp_arr(indx_min,5), temp_arr(indx_min,6), &
                      temp_arr(indx_max,5), temp_arr(indx_max,6))

                 WRITE(lu,"(F13.5,8(2X,F9.4))") &
                      2400000.5_bp + mjd_utc, &
                      MINVAL(temp_arr(:,4))/rad_deg, & 
                      MAXVAL(temp_arr(:,4))/rad_deg, & 
                      MAXVAL(temp_arr(:,3)), & 
                      SUM(temp_arr(:,1))/REAL(SIZE(temp_arr,dim=1),bp)/rad_deg, &
                      SUM(temp_arr(:,2))/REAL(SIZE(temp_arr,dim=1),bp)/rad_deg, &
                      MAX(angular_distance_ra,angular_distance_dec)/rad_deg, &
                      SUM(temp_arr(:,5))/REAL(SIZE(temp_arr,dim=1),bp)/rad_deg, &
                      SUM(temp_arr(:,6))/REAL(SIZE(temp_arr,dim=1),bp)/rad_deg

                 DEALLOCATE(temp_arr)

              ELSE

                 ! Input orbits correspond to one or more single-point estimates of the pdf.
                 ALLOCATE(temp_arr(1,6))

                 ! RA & DEC
                 CALL rotateToEquatorial(ephemerides_arr(1,j))        
                 comp_coord = getCoordinates(ephemerides_arr(1,j))
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (80)',1)
                    STOP
                 END IF
                 Delta = comp_coord(1)
                 temp_arr(1,5:6) = comp_coord(2:3)

                 ! Compute topocentric ecliptic coordinates
                 ! relative to the Sun
                 CALL rotateToEcliptic(ephemerides_arr(1,j))        
                 comp_coord = getCoordinates(ephemerides_arr(1,j))
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (85)',1)
                    STOP
                 END IF
                 tsclon = comp_coord(2) - sunlon
                 tsclat = comp_coord(3) - sunlat
                 IF (tsclon > pi) THEN
                    tsclon = tsclon - two_pi
                 ELSE IF (tsclon < -pi) THEN
                    tsclon = tsclon + two_pi
                 END IF
                 temp_arr(1,1:2) = (/ tsclon, tsclat /)

                 ! APPARENT BRIGHTNESS:
                 ! Compute phase angle
                 CALL NEW(ccoord, ephemerides_arr(1,j))
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (90)',1)
                    STOP
                 END IF
                 CALL rotateToEcliptic(ccoord)
                 obsy_obj = getPosition(ccoord)
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (95)',1)
                    STOP
                 END IF
                 ephemeris_r2 = DOT_PRODUCT(obsy_obj,obsy_obj)
                 CALL toCartesian(orb_lt_corr_arr2(1,j), frame='ecliptic')
                 pos = getPosition(orb_lt_corr_arr2(1,j))
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (100)',1)
                    STOP
                 END IF
                 heliocentric_r2 = DOT_PRODUCT(pos,pos)
                 cos_obj_phase = 0.5_bp * (heliocentric_r2 + ephemeris_r2 - &
                      observer_r2) / (SQRT(heliocentric_r2) * &
                      SQRT(ephemeris_r2))
                 obj_phase = ACOS(cos_obj_phase)
                 ! Input absolute magnitude
                 H_value = get_cl_option("--H=", HG_arr_storb_in(i,1,1))
                 ! Input slope parameter
                 G_value = get_cl_option("--G=", HG_arr_storb_in(i,1,3))
                 temp_arr(1,3) = getApparentHGMagnitude(H=H_value, &
                      G=G_value, r=SQRT(heliocentric_r2), &
                      Delta=Delta, phase_angle=obj_phase)
                 IF (error) THEN
                    CALL errorMessage('oorb / fou', &
                         'TRACE BACK (105)',1)
                    STOP
                 END IF

                 ! SOLAR ELONGATION
                 vec3 = cross_product(obsy_obj,obsy_sun)
                 temp_arr(1,4) = ATAN2(SQRT(SUM(vec3**2)),DOT_PRODUCT(obsy_obj,obsy_sun))

                 CALL NULLIFY(ccoord)
                 CALL NULLIFY(ephemerides_arr(1,j))
                 CALL NULLIFY(orb_lt_corr_arr2(1,j))

                 CALL NULLIFY(observers(j))

                 ! WIDTH OF EPHEMERIS UNCERTAINTY
                 ! Make sure that the ephemeris is equatorial:
                 DO k=1,6
                    stdev_arr(k) = SQRT(cov_arr(k,k,j)) 
                 END DO

                 ! First, try measuring largest width in RA:
                 angular_distance_ra = 3.0_bp*2.0_bp*stdev_arr(2)
                 ! Second, try measuring largest width in Dec:
                 angular_distance_dec = 3.0_bp*2.0_bp*stdev_arr(3)

                 ! MG 20100202:
                 ! NOTE THAT THE MIN,MAX VALUES FOR SOLELON ARE FOR
                 ! NOW USING A TEMPORARY SOLUTION WHICH PROVIDES AN
                 ! APPROXIMATION OF THE UNCERTAINTY. ASSUMPTION IS THE
                 ! ELONGATION UNCERTAINTY IS MOSTLY DEPENDENT ON THE
                 ! RA UNCERTAINTY AND LESS ON THE DEC UNCERTAINTY.
                 ! MG 20100414:
                 ! NEED TO MAKE SURE THAT THE MIN/MAX FOR SOLELON ARE
                 ! REASONABLE (IE BETWEEN 0 AND 180 DEG AT ALL TIMES).

                 solelon_min = (temp_arr(1,4)-angular_distance_ra/2.0_bp)/rad_deg
                 IF (solelon_min < 0.0_bp) THEN
                    solelon_min = 0.0_bp
                 END IF
                 solelon_max = (temp_arr(1,4)+angular_distance_ra/2.0_bp)/rad_deg
                 IF (solelon_max > 180.0_bp) THEN
                    solelon_max = 180.0_bp
                 END IF

                 WRITE(lu,"(F13.5,8(2X,F9.4))") &
                      2400000.5_bp + mjd_utc, &
                      solelon_min, & 
                      solelon_max, & 
                      temp_arr(1,3), & 
                      temp_arr(1,1)/rad_deg, &
                      temp_arr(1,2)/rad_deg, &
                      MAX(angular_distance_ra,angular_distance_dec)/rad_deg, &
                      temp_arr(:,5:6)/rad_deg

                 DEALLOCATE(temp_arr)

              END IF

              CALL NULLIFY(t)

           END DO

           DEALLOCATE(observers, obsy_code_arr)
           DEALLOCATE(ephemerides_arr, orb_lt_corr_arr2)
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
        DEALLOCATE(storb_arr_in)

     ELSE

        CALL errorMessage("oorb / fou", &
             "Uncertainty information not available for input orbits.",1)
        STOP

     END IF

  CASE ("obs_timespan")

     !! Returns the total observational timespan of the input
     !! observations.

     WRITE(stdout,"(F20.6)") getObservationalTimespan(obss_in)

  CASE ("obs_angular_arc")

     !! Returns the total observational angular arc of the input
     !! observations (first to last).

     WRITE(stdout,"(F14.10)") getObservationalAngularArc(obss_in)/rad_deg

  CASE ("uncertainty")

     !! Fractional (dx/x) and absolute uncertainty (dx) on a,e,i;
     !! where dx either covers 68.27% of the total probability mass or
     !! is equal to the 1-sigma limits.

     DO i=1,SIZE(storb_arr_in)
        IF (containsDiscretePDF(storb_arr_in(i))) THEN
           orb_arr_in => getSampleOrbits(storb_arr_in(i), probability_mass=0.6827_bp)
           ALLOCATE(elements_arr(SIZE(orb_arr_in),6))
           DO j=1,SIZE(orb_arr_in)
              elements_arr(j,:) = getElements(orb_arr_in(j), "keplerian")
              elements_arr(j,3:6) = elements_arr(j,3:6)/rad_deg
              CALL NULLIFY(orb_arr_in(j))
           END DO
           pdf_arr_in => getDiscretePDF(storb_arr_in(i), "keplerian")
           j = MAXLOC(pdf_arr_in,dim=1)
           orb = getSampleOrbit(storb_arr_in(i),j)
           elements = getElements(orb, "keplerian")
           CALL NULLIFY(orb)
           elements(3) = elements(3)/rad_deg
           DO j=1,3
              bounds(1) = MINVAL(elements_arr(:,j))
              bounds(2) = MAXVAL(elements_arr(:,j))
              WRITE(stdout,'(2(F12.6,1X))',advance="NO") &
                   (bounds(2)-bounds(1))/elements(j), &
                   bounds(2)-bounds(1)
           END DO
           WRITE(stdout,*)
           DEALLOCATE(orb_arr_in, elements_arr, pdf_arr_in)
        ELSE
           orb = getNominalOrbit(storb_arr_in(i))
           elements = getElements(orb, "keplerian")
           elements(3:6) = elements(3:6)/rad_deg
           CALL NULLIFY(orb)
           cov = getCovarianceMatrix(storb_arr_in(i), "keplerian")
           cov(:,3:) = cov(:,3:)/rad_deg
           cov(3:,:) = cov(3:,:)/rad_deg
           DO j=1,3
              WRITE(stdout,'(2(F12.6,1X))',advance="NO") &
                   (1.0_bp*2.0_bp*SQRT(cov(j,j)))/elements(j), &
                   1.0_bp*2.0_bp*SQRT(cov(j,j))
           END DO
           WRITE(stdout,*)
        END IF
     END DO




  CASE ("uncertainty.old")

     !! Fractional (dx/x) and absolute uncertainty (dx) on a,e,i;
     !! where dx either covers 68.27% of the total probability mass or
     !! is equal to the 1-sigma limits.

     DO i=1,SIZE(storb_arr_in)
        IF (containsDiscretePDF(storb_arr_in(i))) THEN
           orb_arr_in => getSampleOrbits(storb_arr_in(i))
           ALLOCATE(elements_arr(SIZE(orb_arr_in),6))
           DO j=1,SIZE(orb_arr_in)
              elements_arr(j,:) = getElements(orb_arr_in(j), "keplerian")
              elements_arr(j,3:6) = elements_arr(j,3:6)/rad_deg
              CALL NULLIFY(orb_arr_in(j))
           END DO
           pdf_arr_in => getDiscretePDF(storb_arr_in(i), "keplerian")
           DO j=1,3
              CALL confidence_limits(elements_arr(:,j), pdf_arr_in, &
                   probability_mass=0.6827_bp, peak=peak, bounds=bounds, &
                   errstr=errstr)
              WRITE(stdout,'(2(F12.6,1X))',advance="NO") &
                   (bounds(2)-bounds(1))/peak, &
                   bounds(2)-bounds(1)
           END DO
           WRITE(stdout,*)
           DEALLOCATE(orb_arr_in, elements_arr, pdf_arr_in)
        ELSE
           orb = getNominalOrbit(storb_arr_in(i))
           elements = getElements(orb, "keplerian")
           elements(3:6) = elements(3:6)/rad_deg
           CALL NULLIFY(orb)
           cov = getCovarianceMatrix(storb_arr_in(i), "keplerian")
           cov(:,3:) = cov(:,3:)/rad_deg
           cov(3:,:) = cov(3:,:)/rad_deg
           DO j=1,3
              WRITE(stdout,'(2(F12.6,1X))',advance="NO") &
                   (1.0_bp*2.0_bp*SQRT(cov(j,j)))/elements(j), &
                   1.0_bp*2.0_bp*SQRT(cov(j,j))
           END DO
           WRITE(stdout,*)
        END IF
     END DO





  CASE ("afrac")

     !! Fractional uncertainty on the semimajor axis; da/a where da
     !! either covers 99.73% of the total probability mass or is equal
     !! to the 3-sigma limits.

     DO i=1,SIZE(storb_arr_in)
        IF (containsDiscretePDF(storb_arr_in(i))) THEN
           orb_arr_in => getSampleOrbits(storb_arr_in(i))
           ALLOCATE(elements_arr(SIZE(orb_arr_in),6))
           DO j=1,SIZE(orb_arr_in)
              elements_arr(j,:) = getElements(orb_arr_in(j), "keplerian")
              CALL NULLIFY(orb_arr_in(j))
           END DO
           pdf_arr_in => getDiscretePDF(storb_arr_in(i), "keplerian")
           CALL confidence_limits(elements_arr(:,1), pdf_arr_in, &
                probability_mass=0.9973_bp, peak=peak, bounds=bounds, &
                errstr=errstr)
           WRITE(stdout,'(F12.6)') (bounds(2)-bounds(1))/peak
           DEALLOCATE(orb_arr_in, elements_arr, pdf_arr_in)
        ELSE
           orb = getNominalOrbit(storb_arr_in(i))
           elements = getElements(orb, "keplerian")
           CALL NULLIFY(orb)
           cov = getCovarianceMatrix(storb_arr_in(i), "keplerian")
           WRITE(stdout,'(F12.6)') (3.0_bp*2.0_bp*SQRT(cov(1,1)))/elements(1)           
        END IF
     END DO

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
     WRITE(stdout,'(F20.12)') mjd_utc
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
     WRITE(stdout,'(F20.12)') mjd_tai
     CALL NULLIFY(epoch)


  CASE ("utc2tt")

     mjd_utc = get_cl_option("--epoch-mjd-utc=", 0.0_bp)
     CALL NEW(epoch, mjd_utc, "UTC")
     IF (error) THEN
        CALL errorMessage("oorb / utc2tt", &
             "TRACE BACK (5)", 1)
        STOP
     END IF
     mjd_tt = getMJD(epoch, "TT")
     IF (error) THEN
        CALL errorMessage("oorb / utc2tt", &
             "TRACE BACK (10)", 1)
        STOP
     END IF
     WRITE(stdout,'(F20.12)') mjd_tt
     CALL NULLIFY(epoch)


  CASE ("deg2radec", "degreestosexagesimal")

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
          get_cl_option("--9char-mpc", .FALSE.))) THEN     
        CALL errorMessage("oorb / encode_designation", &
             "Either the '--mpc' or the '--9char-mpc' option has to be specified.", 1)                
        STOP
     ELSE IF (get_cl_option("--mpc", .FALSE.) .AND. &
          get_cl_option("--9char-mpc", .FALSE.)) THEN     
        CALL errorMessage("oorb / encode_designation", &
             "Both '--mpc' and '--9char-mpc' options cannot be used simultaneously.", 1)                
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
           ELSE IF (get_cl_option("--9char-mpc", .FALSE.)) THEN
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
           id = ""
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

  CASE ("synthetic_astrometry")

     noise = .NOT.get_cl_option("--no-noise", .FALSE.)

     IF (ALLOCATED(storb_arr_in)) THEN
        ALLOCATE(orb_arr_in(SIZE(storb_arr_in)), id_arr_in(SIZE(storb_arr_in)))
        DO i=1,SIZE(storb_arr_in)
           orb_arr_in(i) = getNominalOrbit(storb_arr_in(i))
           id_arr_in(i) = id_arr_storb_in(i)
        END DO
     END IF

     norb = SIZE(orb_arr_in)
     nobj = SIZE(obss_sep)
     DO i=1,norb

        CALL setParameters(orb_arr_in(i), &
             dyn_model=dyn_model, &
             perturbers=perturbers, &
             asteroid_perturbers=asteroid_perturbers, &
             integrator=integrator, &
             integration_step=integration_step)

        ! find astrometry for this orbit
        iobj = -1
        DO j=1,nobj
           IF (TRIM(id_arr_in(i)) == TRIM(getID(obss_sep(j)))) THEN
              iobj = j
           END IF
        END DO
        IF (iobj < 0 ) THEN
           WRITE(stderr,*) "astrometry not available for " // TRIM(id_arr_in(i))  
           STOP 
        END IF

        observers => getObservatoryCCoords(obss_sep(iobj))
        obs_arr => getObservations(obss_sep(iobj))
        mean = 0.0_bp
        CALL getEphemerides(orb_arr_in(i), observers, ephemerides)
        IF (error) THEN
           IF (error) THEN
              CALL errorMessage('oorb / synthetic_astrometry', &
                   'TRACE BACK (15)',1)
              STOP
           END IF
        END IF
        DO j=1,SIZE(ephemerides)
           CALL rotateToEquatorial(ephemerides(j))
           CALL setObservationSCoord(obs_arr(j), ephemerides(j))
           IF (noise) THEN
              cov = getCovarianceMatrix(obs_arr(j))
              CALL addMultinormalDeviate(obs_arr(j), mean, cov, combined_covariance=.FALSE.)
           END IF
           IF (error) THEN
              CALL errorMessage('oorb / synthetic_astrometry', &
                   'TRACE BACK (20)',1)
              STOP
           END IF
        END DO
        CALL NULLIFY(orb)

        CALL NEW(obss, obs_arr)
        CALL writeObservationFile(obss, stdout, &
             TRIM(observation_format_out))
        CALL NULLIFY(obss)

        DO j=1,SIZE(observers)
           CALL NULLIFY(obs_arr(j))
           CALL NULLIFY(observers(j))
           CALL NULLIFY(ephemerides(j))
        END DO
        DEALLOCATE(observers, obs_arr, ephemerides)

     END DO


  CASE ("add_noise")

     mean = 0.0_bp
     cov = 0.0_bp
     cov(2,2) = get_cl_option("--sigma-ra=", cov(2,2))
     cov(3,3) = get_cl_option("--sigma-dec=", cov(3,3))
     cov(2:3,2:3) = (cov(2:3,2:3)*rad_asec)**2
     DO i=1,SIZE(obss_sep)
        ALLOCATE(mean_arr(getNrOfObservations(obss_sep(i)),6), &
             cov_arr(getNrOfObservations(obss_sep(i)),6,6))
        DO j=1,SIZE(mean_arr,dim=1)
           mean_arr(j,:) = mean
           cov_arr(j,:,:) = cov
        END DO
        CALL addMultinormalDeviates(obss_sep(i), mean_arr, cov_arr)
        CALL writeObservationFile(obss_sep(i), stdout, &
             TRIM(observation_format_out))
        DEALLOCATE(mean_arr, cov_arr)
     END DO


  CASE ("synthetic_orbits")

     norb = get_cl_option("--norb=", 10000)
     amin = get_cl_option("--a-min=", 2.5_bp)
     amax = get_cl_option("--a-max=", 2.6_bp)
     emin = get_cl_option("--e-min=", 0.56_bp)
     emax = get_cl_option("--e-max=", 0.6_bp)
     imin = get_cl_option("--i-min=", 0.0_bp)
     imax = get_cl_option("--i-max=", 4.0_bp)
     imin = imin*rad_deg
     imax = imax*rad_deg
     Hmin = get_cl_option("--H-min=", 17.5_bp)
     Hmax = get_cl_option("--H-max=", 18.0_bp)
     DO i=1,norb
        CALL randomNumber(ran7)
        elements(1) = amin + (amax - amin)*ran7(1)
        elements(2) = emin + (emax - emin)*ran7(2)
        elements(3) = imin + (imax - imin)*ran7(3)
        elements(4) = two_pi*ran7(4)
        elements(5) = two_pi*ran7(5)
        elements(6) = two_pi*ran7(6)
        H = Hmin + (Hmax - Hmin)*ran7(7)
        elements(3:6) = elements(3:6)/rad_deg
        WRITE(*,*) i, "KEP", elements, H, "56000.0 1 6 -1 OpenOrb"
     END DO

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
  IF (ASSOCIATED(id_arr_storb_in)) THEN
     DEALLOCATE(id_arr_storb_in)
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
  IF (exist(obss_in)) THEN
     CALL NULLIFY(obss_in)
  END IF
  IF (ASSOCIATED(obss_sep)) THEN
     DO i=1,SIZE(obss_sep)
        CALL NULLIFY(obss_sep(i))
     END DO
     DEALLOCATE(obss_sep)
  END IF
  IF (ASSOCIATED(HG_arr_in)) THEN
     DEALLOCATE(HG_arr_in)
  END IF
  CALL JPL_ephemeris_nullify()
  CALL nullifyTime()
  IF (ASSOCIATED(perturbers)) THEN
     DEALLOCATE(perturbers, stat=err)
  END IF
  DEALLOCATE(element_type_pdf_arr_in, stat=err)
  DEALLOCATE(cov_arr_in, stat=err)
  DEALLOCATE(pdf_arr_in, stat=err)
  DEALLOCATE(rchi2_arr_in, stat=err)
  DEALLOCATE(jac_arr_in, stat=err)
  DEALLOCATE(reg_apr_arr_in, stat=err)
  DEALLOCATE(HG_arr_storb_in, stat=err)

END PROGRAM oorb
