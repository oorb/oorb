!====================================================================!
!                                                                    !
! Copyright 2002,2003,2004,2005,2006,2007,2008,2009,2010             !
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
!! *Module*description*:
!!
!! Contains routines for IO of options, orbits and residuals. To be
!! called from main programs.
!!
!! @author  MG, JV
!! @version 2010-07-13
!!
MODULE io

  USE Base_cl
  USE File_cl
  USE Time_cl
  USE SphericalCoordinates_cl
  USE Orbit_cl
  USE StochasticOrbit_cl

  USE utilities
  IMPLICIT NONE

CONTAINS




  !! *Description*:
  !!
  !! Decodes an encoded MPC type date.
  !!
  !! Example:
  !!
  !! K01AM -> 2001 11 22
  !!
  !! Returns error.
  !!
  SUBROUTINE decodeMPCDate(encoded_date, year, month, day)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: encoded_date
    INTEGER, INTENT(out) :: year, month
    REAL(bp), INTENT(out) :: day

    REAL(bp) :: tmp
    INTEGER :: i

    CALL removeLeadingBlanks(encoded_date)
    CALL toInt(encoded_date(2:3), year, error)
    IF (encoded_date(1:1) == "J") THEN
       year = year + 1900
    ELSE IF (encoded_date(1:1) == "K") THEN
       year = year + 2000
    END IF
    DO i=0,SIZE(mpc_conv_table)-1
       IF (mpc_conv_table(i) == encoded_date(4:4)) THEN
          EXIT
       END IF
    END DO
    month = i
    DO i=0,SIZE(mpc_conv_table)-1
       IF (mpc_conv_table(i) == encoded_date(5:5)) THEN
          EXIT
       END IF
    END DO
    day = REAL(i,bp)
    IF (LEN_TRIM(encoded_date) > 5) THEN
       CALL toReal("0." // TRIM(encoded_date(6:)), tmp, error)
       day = day + tmp
    END IF

  END SUBROUTINE decodeMPCDate





  SUBROUTINE makeResidualStamps(storb, obss, fname, compute, residuals)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: storb
    TYPE (Observations), INTENT(in) :: obss
    CHARACTER(len=*), INTENT(in) :: fname
    LOGICAL, INTENT(in), OPTIONAL :: compute
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL :: residuals

    TYPE (Orbit), DIMENSION(:), POINTER  :: orb_arr
    TYPE (File) :: tmpfile
    TYPE (Time) :: t
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: observed_scoords, &
         computed_scoords
    TYPE (Observation), DIMENSION(:), POINTER :: obs_arr
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=32), DIMENSION(:), ALLOCATABLE :: dates
    CHARACTER(len=DESIGNATION_LEN) :: id
    CHARACTER(len=64) :: str, str1, str2
    REAL(bp), DIMENSION(:,:,:), POINTER :: residuals_
    REAL(bp), DIMENSION(:,:), POINTER :: res_accept_prm
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: observed_coords, &
         computed_coords
    REAL(bp) :: day, xsize, ysize
    INTEGER :: i, j, k, l, nobs, err, nrow, nrow_max, ncolumn_max, page, &
         year, month
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks
    LOGICAL :: compute_

    ! NOTE: info on obs_mask is not output, should it be?.

    IF (.NOT. exist(storb)) THEN
       error = .TRUE.
       CALL errorMessage("io / makeResidualStamps", &
            "StochasticOrbit object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("io / makeResidualStamps", &
            "Observations object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(compute)) THEN
       compute_ = compute
    ELSE
       compute_ = .FALSE.
    END IF

    IF (.NOT. compute_) THEN
       residuals_ => getResidualDistribution(storb)
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "Residuals are not available.", 1)
          RETURN
       END IF
    ELSE
       ! Compute residuals
       observed_scoords => getObservationSCoords(obss)
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       obsy_ccoords => getObservatoryCCoords(obss)
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       nobs = SIZE(observed_scoords,dim=1)
       orb_arr => getSampleOrbits(storb)
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       ALLOCATE(residuals_(SIZE(orb_arr,dim=1),nobs,6), &
            observed_coords(nobs,6), computed_coords(nobs,6), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / makeResidualStamps", &
               "Could not allocate memory (5).", 1)
          RETURN
       END IF
       DO j=1,SIZE(orb_arr, dim=1)
          CALL getEphemerides(orb_arr(j), obsy_ccoords, &
               computed_scoords)
          IF (error) THEN
             CALL errorMessage("io / makeResidualStamps", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
          observed_coords = 0.0_bp
          computed_coords = 0.0_bp
          DO i=1,nobs
             observed_coords(i,:) = getCoordinates(observed_scoords(i))
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (25)", 1)
                RETURN
             END IF
             computed_coords(i,:) = getCoordinates(computed_scoords(i))
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (30)", 1)
                RETURN
             END IF
          END DO
          residuals_(j,1:nobs,1:6) = observed_coords(1:nobs,1:6) - &
               computed_coords(1:nobs,1:6)        
          residuals_(j,1:nobs,2) = residuals_(j,1:nobs,2) * &
               COS(observed_coords(1:nobs,3))
          DEALLOCATE(computed_scoords, stat=err)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / makeResidualStamps", &
                  "Could not deallocate memory (5).", 1)
             RETURN
          END IF
       END DO
       DO i=1,SIZE(orb_arr)
          CALL NULLIFY(orb_arr(i))
       END DO
       DEALLOCATE(observed_coords, computed_coords, observed_scoords, &
            obsy_ccoords, orb_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / makeResidualStamps", &
               "Could not deallocate memory (10).", 1)
          DEALLOCATE(observed_coords, stat=err)
          DEALLOCATE(computed_coords, stat=err)
          DEALLOCATE(observed_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          RETURN
       END IF
    END IF

    obs_arr => getObservations(obss)
    IF (error) THEN
       CALL errorMessage("io / makeResidualStamps", &
            "TRACE BACK (35)", 1)
       RETURN
    END IF
    nobs = SIZE(obs_arr, dim=1)

    ALLOCATE(dates(nobs), codes(nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / makeResidualStamps", &
            "Could not allocate memory (10).", 1)
       RETURN
    END IF

    DO i=1,nobs
       t = getTime(obs_arr(i))
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (40)", 1)
          RETURN
       END IF
       CALL getCalendarDate(t, "TT", year, month, day)
       IF (year >= 1972) THEN
          dates(i) = getCalendarDateString(t, "UTC")
       ELSE
          dates(i) = getCalendarDateString(t, "UT1")
       END IF
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       codes(i) = getCode(obs_arr(i))
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (50)", 1)
          RETURN
       END IF
    END DO

    id = getID(obss)
    IF (error) THEN
       CALL errorMessage("io / makeResidualStamps", &
            "TRACE BACK (55)", 1)
       RETURN
    END IF
    CALL toString(nobs, str, error)
    IF (error) THEN
       CALL errorMessage("io / makeResidualStamps", &
            "TRACE BACK (60)", 1)
       RETURN
    END IF

    CALL getParameters(storb, res_accept=res_accept_prm)
    IF (error) THEN
       CALL errorMessage("io / makeResidualStamps", &
            "TRACE BACK (65)", 1)
       RETURN
    END IF
    obs_masks => getObservationMasks(storb)
    IF (error) THEN
       CALL errorMessage("io / makeResidualStamps", &
            "TRACE BACK (70)", 1)
       RETURN
    END IF

    nrow_max = 4
    ncolumn_max = 4
    nrow = CEILING(nobs/REAL(ncolumn_max))

    str = " "
    str = TRIM(id) // "_" // TRIM(str) // "resid_stamp.gp"
    CALL system("echo set terminal postscript eps enhanced color linewidth 1.0 8.0 > " // TRIM(str))
    CALL system("echo set out \'residual_stamps.eps\' >> " // TRIM(str))
    CALL system("echo unset key >> " // TRIM(str))
    CALL system("echo set pointsize 0.2 >> " // TRIM(str))
    IF (nrow >= nrow_max) THEN
       xsize = 1.0
       ysize = 1.0
       CALL system("echo set size 1.0,1.0 >> " // TRIM(str))
    ELSE IF (nrow < nrow_max) THEN
       IF (nobs > ncolumn_max) THEN
          xsize = 1.0
          ysize = 1.0/nrow_max*nrow
          CALL toString(ysize, str1, error)
          IF (error) THEN
             CALL errorMessage("io / makeResidualStamps", &
                  "TRACE BACK (80)", 1)
             RETURN
          END IF
          CALL system("echo set size 1.0," // TRIM(str1) // " >> " // TRIM(str))
       ELSE
          xsize = 1.0/ncolumn_max*nobs
          ysize = 1.0/nrow_max*nrow
          CALL toString(xsize, str1, error)
          IF (error) THEN
             CALL errorMessage("io / makeResidualStamps", &
                  "TRACE BACK (80)", 1)
             RETURN
          END IF
          CALL toString(ysize, str2, error)
          IF (error) THEN
             CALL errorMessage("io / makeResidualStamps", &
                  "TRACE BACK (80)", 1)
             RETURN
          END IF
          CALL system("echo set size " // TRIM(str1) // "," // &
               TRIM(str2) // " >> " // TRIM(str))
       END IF
    END IF
    CALL system("echo set origin 0.0,0.0 >> " // TRIM(str))
    page = 1

    DO
       !       CALL system('echo set multiplot title \"'  // TRIM(id) // &
       !            '\" >> ' // TRIM(str))
       CALL system('echo set multiplot >> ' // TRIM(str))
       CALL toString(1.0/ncolumn_max, str1, error)
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (75)", 1)
          RETURN
       END IF
       CALL toString(1.0/nrow_max, str2, error)
       IF (error) THEN
          CALL errorMessage("io / makeResidualStamps", &
               "TRACE BACK (80)", 1)
          RETURN
       END IF
       CALL system("echo set size " // TRIM(str1) // &
            "," // TRIM(str2) // " >> " // TRIM(str))
       DO i=1,nrow
          DO j=1,ncolumn_max
             k = (page-1)*nrow_max*ncolumn_max + (i-1)*ncolumn_max + j
             IF (k > MIN(nobs,page*nrow_max*ncolumn_max)) THEN
                ! All residuals plotted for this page or overall
                EXIT
             END IF
             IF (nobs >= ncolumn_max) THEN
                CALL toString((j-1)*(xsize/ncolumn_max), str1, error)
             ELSE
                CALL toString((j-1)*(xsize/nobs), str1, error)                
             END IF
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (85)", 1)
                RETURN
             END IF
             CALL toString(ysize-i*(1.0/nrow_max), str2, error)
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (90)", 1)
                RETURN
             END IF
             CALL system("echo set origin " // TRIM(str1) // &
                  "," // TRIM(str2) // " >> " // TRIM(str))
             CALL system("echo set size square >> " // TRIM(str))
             CALL system("echo set xzeroaxis >> " // TRIM(str))
             CALL system("echo set yzeroaxis >> " // TRIM(str))
             CALL toString(res_accept_prm(k,2)/rad_asec, &
                  str2, error)
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (95)", 1)
                RETURN
             END IF
             str1 = "-" // TRIM(str2)
             CALL system("echo set xrange [" // TRIM(str1) // &
                  ":" // TRIM(str2) // "] >> " // TRIM(str))          
             CALL toString(res_accept_prm(k,3)/rad_asec, &
                  str2, error)
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (100)", 1)
                RETURN
             END IF
             str1 = "-" // TRIM(str2)
             CALL system("echo set yrange [" // TRIM(str1) // &
                  ":" // TRIM(str2) // "] >> " // TRIM(str))
             !             CALL system('echo set title \"' // TRIM(id) // "\\" // "n" // &
             CALL system('echo set title \"' // TRIM(codes(k)) // " " &
                  // TRIM(dates(k)) // '\" >> ' // TRIM(str))
             CALL system("echo set xlabel \'{/Symbol D}RA [asec]\' >> " // &
                  TRIM(str))
             CALL system("echo set ylabel \'{/Symbol D}Dec [asec]\' >> " // &
                  TRIM(str))
             CALL toString(k, str1, error)
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (105)", 1)
                RETURN
             END IF
             str1 = TRIM(str) // "_tmp" // TRIM(str1)
             CALL NEW(tmpfile, TRIM(str1))
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (110)", 1)
                RETURN
             END IF
             CALL OPEN(tmpfile)
             IF (error) THEN
                CALL errorMessage("io / makeResidualStamps", &
                     "TRACE BACK (115)", 1)
                RETURN
             END IF
             DO l=1,SIZE(residuals_,dim=1)
                WRITE(getUnit(tmpfile),*) residuals_(l,k,2:3)/rad_asec
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / makeResidualStamps", &
                        "Write error (15).", 1)
                   RETURN
                END IF
             END DO
             CALL NULLIFY(tmpfile)
             IF (ALL(obs_masks(k,2:3))) THEN
                CALL system("echo plot \'" // TRIM(str1) // &
                     "\' using 1:2 0 >> " // TRIM(str))
             ELSE
                CALL system("echo plot \'" // TRIM(str1) // &
                     "\' using 1:2 1 >> " // TRIM(str))
             END IF
          END DO
       END DO
       CALL system("echo unset multiplot >> " // TRIM(str))
       IF (nobs > page*nrow_max*ncolumn_max) THEN
          ! Turn to next page if there is stamps to plot
          page = page + 1
       ELSE
          ! All stamps have been plotted
          EXIT
       END IF
    END DO
    CALL system("gnuplot " // TRIM(str))
    CALL system("mv residual_stamps.eps " // TRIM(fname))
    CALL system("rm -f " // TRIM(str) // " " // TRIM(str) // "_tmp*")

    IF (PRESENT(residuals)) THEN
       IF (ASSOCIATED(residuals)) THEN
          DEALLOCATE(residuals, stat=err)
       END IF
       ALLOCATE(residuals(SIZE(residuals_,dim=1),nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / makeResidualStamps", &
               "Could not allocate memory (15).", 1)
          RETURN
       END IF
       residuals = residuals_
    END IF

    DEALLOCATE(obs_arr, residuals_, dates, codes, res_accept_prm, &
         obs_masks, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / makeResidualStamps", &
            "Could not deallocate memory (15).", 1)
       DEALLOCATE(obs_arr, stat=err)
       DEALLOCATE(residuals_, stat=err)
       DEALLOCATE(dates, stat=err)
       DEALLOCATE(codes, stat=err)
       DEALLOCATE(res_accept_prm, stat=err)
       DEALLOCATE(obs_masks, stat=err)
       RETURN
    END IF

  END SUBROUTINE makeResidualStamps





  SUBROUTINE readConfigurationFile(config_file, &
       info_verb, err_verb, planetary_ephemeris_fname, &
       task, &
       t0, &
       multiple_ids, &
       obs_stdev_arr, &
       element_type_comp, &
       element_type_out, &
       observation_format_out, orbit_format_out, &
       plot_results, &
       plot_open, &
       dyn_model, perturbers, integrator, integration_step, &
       dyn_model_init, integrator_init, integration_step_init, &
       accwin_multiplier, &
       uniform_pdf, regularized_pdf, &
       pdf_ml, & 
       outlier_rejection, outlier_multiplier, &
       apriori_a_min, apriori_a_max, apriori_periapsis_min, &
       apriori_periapsis_max, apriori_apoapsis_min, &
       apriori_apoapsis_max, apriori_rho_min, apriori_rho_max, &
       sor_type, sor_2point_method, sor_2point_method_sw, sor_norb, &
       sor_ntrial, sor_norb_sw, sor_ntrial_sw, sor_niter, &
       sor_random_obs, sor_rho_init, sor_genwin_multiplier, &
       sor_genwin_offset, sor_rho_gauss, sor_rho_init2, &
       sor_rho_gauss2, sor_iterate_bounds, &
       vov_type, vov_norb, vov_ntrial, vov_niter, vov_norb_iter, &
       vov_ntrial_iter, vov_nmap, vov_mapping_mask, vov_scaling, &
       ls_type, ls_element_mask, ls_correction_factor, ls_rchi2_max, &
       ls_niter_major_max, ls_niter_major_min, ls_niter_minor, &
       cos_nsigma, cos_norb, cos_ntrial, cos_gaussian, &
       smplx_tol, smplx_niter, &
       pp_H_estimation, pp_G, pp_G_unc, &
       eph_lt_correction, eph_dt_since_last_obs, eph_obsy_code, &
       eph_date, &
       masked_obs, write_residuals, &
       sor_rhofname)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: &
         config_file
    TYPE (Time), DIMENSION(:), POINTER, OPTIONAL :: &
         eph_date
    TYPE (Time), INTENT(inout), OPTIONAL :: &
         t0
    CHARACTER(len=*), DIMENSION(:), POINTER, OPTIONAL :: &
         eph_obsy_code
    CHARACTER(len=*), INTENT(inout), OPTIONAL :: &
         sor_rhofname, &
         element_type_comp, &
         element_type_out, &
         dyn_model, &
         dyn_model_init, &
         integrator, &
         integrator_init, &
         sor_2point_method, &
         sor_2point_method_sw, &
         observation_format_out, &
         orbit_format_out, &
         planetary_ephemeris_fname
    REAL(bp), DIMENSION(6,2), INTENT(inout), OPTIONAL :: &
         vov_scaling
    REAL(bp), DIMENSION(:), POINTER, OPTIONAL :: &
         eph_dt_since_last_obs
    REAL(bp), DIMENSION(6), INTENT(inout), OPTIONAL :: &
         obs_stdev_arr
    REAL(bp), DIMENSION(4), INTENT(inout), OPTIONAL :: &
         sor_rho_init, &
         sor_rho_init2, &
         sor_genwin_offset
    REAL(bp), INTENT(inout), OPTIONAL :: &
         accwin_multiplier, &
         outlier_multiplier, &
         sor_genwin_multiplier, &
         ls_correction_factor, &
         ls_rchi2_max, &
         cos_nsigma, &
         smplx_tol, &
         integration_step, &
         integration_step_init, &
         pdf_ml, &
         apriori_a_min, &
         apriori_a_max, &
         apriori_periapsis_min, &
         apriori_periapsis_max, &
         apriori_apoapsis_min, &
         apriori_apoapsis_max, &
         apriori_rho_min, &
         apriori_rho_max, &
         pp_G, pp_G_unc
    INTEGER, INTENT(inout), OPTIONAL :: &
         info_verb, &
         err_verb, &
         task, &
         sor_type, &
         sor_norb, &
         sor_ntrial, &
         sor_norb_sw, &
         sor_ntrial_sw, &
         sor_niter, &
         vov_type, &
         vov_norb, &
         vov_ntrial, &
         vov_niter, &
         vov_norb_iter, &
         vov_ntrial_iter, &
         vov_nmap, &
         ls_type, &
         ls_niter_major_max, &
         ls_niter_major_min, &
         ls_niter_minor, &
         cos_norb, &
         cos_ntrial, &
         smplx_niter
    LOGICAL, DIMENSION(:), POINTER, OPTIONAL :: &
         eph_lt_correction, &
         perturbers
    LOGICAL, DIMENSION(6), INTENT(inout), OPTIONAL :: &
         ls_element_mask, &
         vov_mapping_mask
    LOGICAL, DIMENSION(4), INTENT(inout), OPTIONAL :: &
         sor_iterate_bounds
    LOGICAL, INTENT(inout), OPTIONAL :: &
         plot_results, &
         plot_open, &
         multiple_ids, &
         uniform_pdf, &
         regularized_pdf, &
         masked_obs, &
         outlier_rejection, &
         sor_random_obs, &
         sor_rho_gauss, sor_rho_gauss2, &
         cos_gaussian, &
         write_residuals, &
         pp_H_estimation

    CHARACTER(len=256) :: line, par_id, par_val
    REAL(bp) :: day, jd, mjd
    INTEGER :: err, i, indx, indx1, indx2, year, month
    LOGICAL :: last

    REWIND(getUnit(config_file))
    DO
       READ(getUnit(config_file),"(A)", iostat=err) line
       IF (err > 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readConfigurationFile", &
               "Error while reading option file.", 1)
          RETURN
       ELSE IF (err < 0) THEN ! end-of-file
          EXIT
       END IF
       IF (LEN_TRIM(line) == 0 .OR. line(1:1) == "#") THEN
          CYCLE
       END IF
       indx = INDEX(line,":")
       par_id = line(1:indx-1)
       par_val = line(indx+1:LEN_TRIM(line))
       CALL removeLeadingBlanks(par_val)
       indx = INDEX(par_id,".",back=.TRUE.)
       IF (IACHAR(par_id(indx+1:indx+1)) >= 48 .AND. &
            IACHAR(par_id(indx+1:indx+1)) <= 57) THEN
          CALL toInt(TRIM(par_id(indx+1:)), i, error)          
          par_id = par_id(1:indx-1)
       ELSE
          i = 0
       END IF
       SELECT CASE(TRIM(par_id))

          ! GENERAL PARAMETERS
       CASE ("verbose.info")
          IF (PRESENT(info_verb)) THEN
             CALL toInt(TRIM(par_val), info_verb, error)
          END IF
       CASE ("verbose.error")
          IF (PRESENT(err_verb)) THEN
             CALL toInt(TRIM(par_val), err_verb, error)
          END IF
       CASE ("planetary_ephemeris_fname")
          IF (PRESENT(planetary_ephemeris_fname)) THEN
             planetary_ephemeris_fname = TRIM(par_val)
          END IF
       CASE ("task")
          IF (PRESENT(task)) THEN
             CALL toInt(TRIM(par_val), task, error)
          END IF
       CASE ("write.residuals")
          IF (PRESENT(write_residuals)) THEN
             write_residuals = .TRUE.
          END IF
       CASE ("plot.results")
          IF (PRESENT(plot_results)) THEN
             plot_results = .TRUE.
          END IF
       CASE ("plot.open")
          IF (PRESENT(plot_open)) THEN
             plot_open = .TRUE.
          END IF
       CASE ("element_type_out")
          IF (PRESENT(element_type_out)) THEN
             CALL locase(par_val, error)
             IF (error) THEN
                CALL errorMessage("io / readConfigurationFile", &
                     "The type of output orbital elements " // &
                     "string contains forbidden characters.", 1)
             END IF
             IF (par_val /= "cartesian" .AND. &
                  par_val /= "keplerian" .AND. &
                  par_val /= "cometary") THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Type of output orbital elements " // &
                     "is not recognized: " // TRIM(par_val) &
                     // ".", 1)
             END IF
             element_type_out = TRIM(par_val)
          END IF
       CASE ("observation.format.out")
          IF (PRESENT(observation_format_out)) THEN
             observation_format_out = TRIM(par_val)
          END IF
       CASE ("orbit.format.out")
          IF (PRESENT(orbit_format_out)) THEN
             orbit_format_out = TRIM(par_val)
          END IF

          ! INPUT OBSERVATION PARAMETERS:
       CASE ("sor.rhofname")
          IF (PRESENT(sor_rhofname)) THEN
             sor_rhofname = TRIM(par_val)
          END IF
       CASE ("stdev.ra")
          IF (PRESENT(obs_stdev_arr)) THEN
             CALL toReal(TRIM(par_val), obs_stdev_arr(2), error)
             obs_stdev_arr(2) = obs_stdev_arr(2)*rad_asec
          END IF
       CASE ("stdev.dec")
          IF (PRESENT(obs_stdev_arr)) THEN
             CALL toReal(TRIM(par_val), obs_stdev_arr(3), error)
             obs_stdev_arr(3) = obs_stdev_arr(3)*rad_asec
          END IF
       CASE ("multiple_ids")
          IF (PRESENT(multiple_ids)) THEN
             multiple_ids = .TRUE.
          END IF


          ! GENERAL INVERSION PARAMETERS
       CASE ("element_type_comp")
          IF (PRESENT(element_type_comp)) THEN
             CALL locase(par_val, error)
             IF (error) THEN
                CALL errorMessage("io / readConfigurationFile", &
                     "The type of computation orbital " // &
                     "elements string contains forbidden characters.", 1)
             END IF
             IF (par_val /= "cartesian" .AND. &
                  par_val /= "keplerian" .AND. &
                  par_val /= "cometary") THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Type of orbital elements to be used in computations " // &
                     "is not recognized: " // TRIM(par_val) &
                     // ".", 1)
             END IF
             element_type_comp = TRIM(par_val)
          END IF
       CASE ("epoch.cal")
          IF (PRESENT(t0)) THEN
             indx1 = INDEX(TRIM(par_val),"/")
             indx2 = INDEX(TRIM(par_val),"/",back=.TRUE.)
             CALL toInt(par_val(1:indx1-1), year, error)
             CALL toInt(par_val(indx1+1:indx2-1), month, error)
             CALL toReal(par_val(indx2+1:), day, error)
             IF (exist(t0)) THEN
                CALL NULLIFY(t0)
             END IF
             CALL NEW(t0, year, month, day, "tdt")
          END IF
       CASE ("epoch.mjd")
          IF (PRESENT(t0)) THEN
             CALL toReal(TRIM(par_val), mjd, error)
             IF (exist(t0)) THEN
                CALL NULLIFY(t0)
             END IF
             CALL NEW(t0, mjd, "tdt")
          END IF
       CASE ("epoch.jd")
          IF (PRESENT(t0)) THEN
             CALL toReal(TRIM(par_val), jd, error)
             IF (exist(t0)) THEN
                CALL NULLIFY(t0)
             END IF
             CALL NEW(t0, jd-2400000.5_bp, "tdt")
          END IF
       CASE ("outlier_rejection")
          IF (PRESENT(outlier_rejection)) THEN
             outlier_rejection = .TRUE.
          END IF
       CASE ("outlier.multiplier")
          IF (PRESENT(outlier_multiplier)) THEN
             CALL toReal(TRIM(par_val), outlier_multiplier, error)
          END IF
       CASE ("accwin.multiplier")
          IF (PRESENT(accwin_multiplier)) THEN
             CALL toReal(TRIM(par_val), accwin_multiplier, error)
          END IF


          ! PROPAGATION PARAMETERS:
       CASE ("dynamical_model")
          IF (PRESENT(dyn_model)) THEN
             dyn_model = TRIM(par_val)
          END IF
       CASE ("perturber.Mercury")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(1) = .TRUE.
                CASE ("f", "F")
                   perturbers(1) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Venus")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(2) = .TRUE.
                CASE ("f", "F")
                   perturbers(2) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Earth")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(3) = .TRUE.
                CASE ("f", "F")
                   perturbers(3) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Mars")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(4) = .TRUE.
                CASE ("f", "F")
                   perturbers(4) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Jupiter")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(5) = .TRUE.
                CASE ("f", "F")
                   perturbers(5) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Saturn")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(6) = .TRUE.
                CASE ("f", "F")
                   perturbers(6) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Uranus")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(7) = .TRUE.
                CASE ("f", "F")
                   perturbers(7) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Neptune")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(8) = .TRUE.
                CASE ("f", "F")
                   perturbers(8) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Pluto")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(9) = .TRUE.
                CASE ("f", "F")
                   perturbers(9) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("perturber.Moon")
          IF (PRESENT(perturbers)) THEN
             IF (.NOT.ASSOCIATED(perturbers)) THEN
                ALLOCATE(perturbers(10), stat=err)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not allocate perturbers array.", 1)
                END IF
             END IF
             IF (.NOT.error) THEN
                SELECT CASE (ADJUSTL(par_val))
                CASE ("t", "T")
                   perturbers(10) = .TRUE.
                CASE ("f", "F")
                   perturbers(10) = .FALSE.
                CASE default
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Cannot understand logical value: " // &
                        TRIM(ADJUSTL(par_val)) // ".", 1)
                END SELECT
             END IF
          END IF
       CASE ("integrator")
          IF (PRESENT(integrator)) THEN
             integrator = TRIM(par_val)
          END IF
       CASE ("integration_step")
          IF (PRESENT(integration_step)) THEN
             CALL toReal(TRIM(par_val), integration_step, error)
             IF (integration_step <= 0.0_bp) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Integration step must be larger than zero.", 1)
             END IF
          END IF
       CASE ("dynamical_model_init")
          IF (PRESENT(dyn_model_init)) THEN
             dyn_model_init = TRIM(par_val)
          END IF
       CASE ("integrator_init")
          IF (PRESENT(integrator_init)) THEN
             integrator_init = TRIM(par_val)
          END IF
       CASE ("integration_step_init")
          IF (PRESENT(integration_step_init)) THEN
             CALL toReal(TRIM(par_val), integration_step_init, error)
             IF (integration_step_init <= 0.0_bp) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Integration step must be larger than zero.", 1)
             END IF
          END IF


          ! STATISTICAL PARAMETERS
       CASE ("uni.pdf")
          IF (PRESENT(uniform_pdf)) THEN
             uniform_pdf = .TRUE.
          END IF
       CASE ("reg.pdf")
          IF (PRESENT(regularized_pdf)) THEN
             regularized_pdf = .TRUE.
          END IF
       CASE ("obs.mask")
          IF (PRESENT(masked_obs)) THEN
             masked_obs = .TRUE.
          END IF
       CASE ("pdf.init")
          IF (PRESENT(pdf_ml)) THEN
             CALL toReal(par_val, pdf_ml, error)
          END IF


          ! BAYESIAN (INFORMATIVE) A PRIORI PARAMETERS
       CASE ("apriori.a.min")
          IF (PRESENT(apriori_a_min)) THEN
             CALL toReal(par_val, apriori_a_min, error)
          END IF
       CASE ("apriori.a.max")
          IF (PRESENT(apriori_a_max)) THEN
             CALL toReal(par_val, apriori_a_max, error)
          END IF
       CASE ("apriori.q.min")
          IF (PRESENT(apriori_periapsis_min)) THEN
             CALL toReal(par_val, apriori_periapsis_min, error)
          END IF
       CASE ("apriori.q.max")
          IF (PRESENT(apriori_periapsis_max)) THEN
             CALL toReal(par_val, apriori_periapsis_max, error)
          END IF
       CASE ("apriori.Q.min")
          IF (PRESENT(apriori_apoapsis_min)) THEN
             CALL toReal(par_val, apriori_apoapsis_min, error)
          END IF
       CASE ("apriori.Q.max")
          IF (PRESENT(apriori_apoapsis_max)) THEN
             CALL toReal(par_val, apriori_apoapsis_max, error)
          END IF
       CASE ("apriori.rho.min")
          IF (PRESENT(apriori_rho_min)) THEN
             CALL toReal(par_val, apriori_rho_min, error)
          END IF
       CASE ("apriori.rho.max")
          IF (PRESENT(apriori_rho_max)) THEN
             CALL toReal(par_val, apriori_rho_max, error)
          END IF


          ! STATISTICAL ORBITAL RANGING PARAMETERS:
       CASE ("sor.type")
          IF (PRESENT(sor_type)) THEN
             CALL toInt(TRIM(par_val), sor_type, error)
          END IF
       CASE ("sor.two_point_method")
          IF (PRESENT(sor_2point_method)) THEN
             sor_2point_method = TRIM(par_val)
          END IF
       CASE ("sor.norb")
          IF (PRESENT(sor_norb)) THEN
             CALL toInt(TRIM(par_val), sor_norb, error)
          END IF
       CASE ("sor.ntrial")
          IF (PRESENT(sor_ntrial)) THEN
             CALL toInt(TRIM(par_val), sor_ntrial, error)
          END IF
       CASE ("sor.two_point_method_sw")
          IF (PRESENT(sor_2point_method_sw)) THEN
             sor_2point_method_sw = TRIM(par_val)
          END IF
       CASE ("sor.norb.sw")
          IF (PRESENT(sor_norb_sw)) THEN
             CALL toInt(TRIM(par_val), sor_norb_sw, error)
          END IF
       CASE ("sor.ntrial.sw")
          IF (PRESENT(sor_ntrial_sw)) THEN
             CALL toInt(TRIM(par_val), sor_ntrial_sw, error)
          END IF
       CASE ("sor.niter")
          IF (PRESENT(sor_niter)) THEN
             CALL toInt(TRIM(par_val), sor_niter, error)
          END IF
       CASE ("sor.ran.obs")
          IF (PRESENT(sor_random_obs)) THEN
             sor_random_obs = .TRUE.
          END IF
       CASE ("sor.rho.gauss")
          IF (PRESENT(sor_rho_gauss)) THEN
             sor_rho_gauss = .TRUE.
          END IF
       CASE ("sor.rho.gauss2")
          IF (PRESENT(sor_rho_gauss2)) THEN
             sor_rho_gauss2 = .TRUE.
          END IF
       CASE ("sor.rho11.init")
          IF (PRESENT(sor_rho_init)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init(1), error)
          END IF
       CASE ("sor.rho12.init")
          IF (PRESENT(sor_rho_init)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init(2), error)
          END IF
       CASE ("sor.rho21.init")
          IF (PRESENT(sor_rho_init)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init(3), error)
          END IF
       CASE ("sor.rho22.init")
          IF (PRESENT(sor_rho_init)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init(4), error)
          END IF
       CASE ("sor.genwin.multiplier")
          IF (PRESENT(sor_genwin_multiplier)) THEN
             CALL toReal(TRIM(par_val), sor_genwin_multiplier, error)
          END IF
       CASE ("sor.genwin.offset")
          IF (PRESENT(sor_genwin_offset)) THEN
             indx1 = 1
             last = .FALSE.
             CALL removeLeadingBlanks(par_val)
             i = 0
             DO WHILE (.NOT.last)
                ! Find the location of the next space: 
                indx2 = INDEX(TRIM(par_val(indx1:))," ")
                IF (indx2 == 0) THEN
                   ! No more spaces left -> last parameter
                   last = .TRUE.
                   indx2 = LEN_TRIM(par_val(indx1:))
                END IF
                i = i + 1
                IF (i > 4) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Too many input values for sor.genwin.offset; only 4 allowed!.", 1)
                   RETURN
                END IF
                CALL toReal(TRIM(par_val(indx1:indx1+indx2-1)), sor_genwin_offset(i), error)
                ! Move starting index to the next parameter 
                indx1 = indx1 + indx2
                CALL removeLeadingBlanks(par_val(indx1:))
             END DO
             sor_genwin_offset = sor_genwin_offset*rad_asec
          END IF
       CASE ("sor.rho11.init2")
          IF (PRESENT(sor_rho_init2)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init2(1), error)
          END IF
       CASE ("sor.rho12.init2")
          IF (PRESENT(sor_rho_init2)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init2(2), error)
          END IF
       CASE ("sor.rho21.init2")
          IF (PRESENT(sor_rho_init2)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init2(3), error)
          END IF
       CASE ("sor.rho22.init2")
          IF (PRESENT(sor_rho_init2)) THEN
             CALL toReal(TRIM(par_val), sor_rho_init2(4), error)
          END IF
       CASE ("sor.iterate_bounds")
          IF (PRESENT(sor_iterate_bounds)) THEN
             indx1 = 1
             last = .FALSE.
             CALL removeLeadingBlanks(par_val)
             i = 0
             DO WHILE (.NOT.last)
                ! Find the location of the next space: 
                indx2 = INDEX(TRIM(par_val(indx1:))," ")
                IF (indx2 == 0) THEN
                   ! No more spaces left -> last parameter
                   last = .TRUE.
                   indx2 = LEN_TRIM(par_val(indx1:))
                END IF
                i = i + 1
                IF (i > 4) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Too many input values for sor.iterate_bounds; only 4 allowed!.", 1)
                   RETURN
                END IF
                READ(par_val(indx1:indx1+indx2-1),"(L1)") sor_iterate_bounds(i)
                ! Move starting index to the next parameter 
                indx1 = indx1 + indx2
                CALL removeLeadingBlanks(par_val(indx1:))
             END DO
          END IF

          ! VOLUME-OF-VARIATION PARAMETERS:
       CASE ("vov.type")
          IF (PRESENT(vov_type)) THEN
             CALL toInt(TRIM(par_val), vov_type, error)
          END IF
       CASE ("vov.norb")
          IF (PRESENT(vov_norb)) THEN
             CALL toInt(TRIM(par_val), vov_norb, error)
          END IF
       CASE ("vov.ntrial")
          IF (PRESENT(vov_ntrial)) THEN
             CALL toInt(TRIM(par_val), vov_ntrial, error)
          END IF
       CASE ("vov.niter")
          IF (PRESENT(vov_niter)) THEN
             CALL toInt(TRIM(par_val), vov_niter, error)
          END IF
       CASE ("vov.norb_iter")
          IF (PRESENT(vov_norb_iter)) THEN
             CALL toInt(TRIM(par_val), vov_norb_iter, error)
          END IF
       CASE ("vov.ntrial_iter")
          IF (PRESENT(vov_ntrial_iter)) THEN
             CALL toInt(TRIM(par_val), vov_ntrial_iter, error)
          END IF
       CASE ("vov.nmap")
          IF (PRESENT(vov_nmap)) THEN
             CALL toInt(TRIM(par_val), vov_nmap, error)
          END IF
       CASE ("vov.mapping_mask")
          IF (PRESENT(vov_mapping_mask)) THEN
             vov_mapping_mask = .TRUE.
             READ(par_val, *, iostat=err) vov_mapping_mask
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Could not read parameter value (5).", 1)
             END IF
          END IF
       CASE ("vov.scaling.lo")
          IF (PRESENT(vov_scaling)) THEN
             READ(par_val, *, iostat=err) vov_scaling(:,1)
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Could not read parameter value (10).", 1)
             END IF
          END IF
       CASE ("vov.scaling.hi")
          IF (PRESENT(vov_scaling)) THEN
             READ(par_val, *, iostat=err) vov_scaling(:,2)
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Could not read parameter value (15).", 1)
             END IF
          END IF


          ! LEAST SQUARES PARAMETERS:
       CASE ("ls.type")
          IF (PRESENT(ls_type)) THEN
             CALL toInt(TRIM(par_val), ls_type, error)
          END IF
       CASE ("ls.correction_factor")
          IF (PRESENT(ls_correction_factor)) THEN
             CALL toReal(TRIM(par_val), ls_correction_factor, error)
          END IF
       CASE ("ls.rchi2.max")
          IF (PRESENT(ls_rchi2_max)) THEN
             CALL toReal(TRIM(par_val), ls_rchi2_max, error)
          END IF
       CASE ("ls.element_mask")
          IF (PRESENT(ls_element_mask)) THEN
             ls_element_mask = .TRUE.
             READ(par_val, *, iostat=err) ls_element_mask
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Could not read parameter value (20).", 1)
             END IF
          END IF
       CASE ("ls.niter_major.max")
          IF (PRESENT(ls_niter_major_max)) THEN
             CALL toInt(TRIM(par_val), ls_niter_major_max, error)
          END IF
       CASE ("ls.niter_major.min")
          IF (PRESENT(ls_niter_major_min)) THEN
             CALL toInt(TRIM(par_val), ls_niter_major_min, error)
          END IF
       CASE ("ls.niter_minor")
          IF (PRESENT(ls_niter_minor)) THEN
             CALL toInt(TRIM(par_val), ls_niter_minor, error)
          END IF


          ! COVARIANCE SAMPLING PARAMETERS:
       CASE ("cos.nsigma")
          IF (PRESENT(cos_nsigma)) THEN
             CALL toReal(TRIM(par_val), cos_nsigma, error)
          END IF
       CASE ("cos.norb")
          IF (PRESENT(cos_norb)) THEN
             CALL toInt(TRIM(par_val), cos_norb, error)
          END IF
       CASE ("cos.ntrial")
          IF (PRESENT(cos_ntrial)) THEN
             CALL toInt(TRIM(par_val), cos_ntrial, error)
          END IF
       CASE ("cos.gaussian")
          IF (PRESENT(cos_gaussian)) THEN
             READ(par_val, *, iostat=err) cos_gaussian
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Could not read parameter value (21).", 1)
             END IF
          END IF


          ! SIMPLEX OPTIMIZATION PARAMETERS:
       CASE ("smplx.tol")
          IF (PRESENT(smplx_tol)) THEN
             CALL toReal(TRIM(par_val), smplx_tol, error)
          END IF
       CASE ("smplx.niter")
          IF (PRESENT(smplx_niter)) THEN
             CALL toInt(TRIM(par_val), smplx_niter, error)
          END IF


          ! PHYSICAL PARAMETERS:
       CASE ("pp.H_estimation")
          IF (PRESENT(pp_H_estimation)) THEN
             READ(par_val, *, iostat=err) pp_H_estimation
             IF (err /= 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "Could not read parameter value (25).", 1)
             END IF
          END IF
       CASE ("pp.G")
          IF (PRESENT(pp_G)) THEN
             CALL toReal(TRIM(par_val), pp_G, error)
          END IF
       CASE ("pp.G_unc")
          IF (PRESENT(pp_G_unc)) THEN
             CALL toReal(TRIM(par_val), pp_G_unc, error)
          END IF


          ! PREDICTION PARAMETERS:
       CASE ("eph.observatory.code")
          IF (PRESENT(eph_obsy_code)) THEN
             IF (LEN_TRIM(par_val) == 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "No codes specified for option 'eph.observatory.code'.", 1)
                RETURN
             END IF
             i = 0
             indx1 = 1
             last = .FALSE.
             CALL removeLeadingBlanks(par_val)
             DO WHILE (.NOT.last)
                ! Find the location of the next space: 
                indx2 = INDEX(TRIM(par_val(indx1:))," ")
                IF (indx2 == 0) THEN
                   ! No more spaces left -> last parameter
                   last = .TRUE.
                   indx2 = LEN_TRIM(par_val(indx1:))
                END IF
                i = i + 1
                IF (.NOT.ASSOCIATED(eph_obsy_code)) THEN
                   ALLOCATE(eph_obsy_code(i), stat=err)
                   IF (err /= 0) THEN
                      error = .TRUE.
                      CALL errorMessage("io / readConfigurationFile", &
                           "Could not allocate memory.", 1)
                      RETURN
                   END IF
                ELSE IF (SIZE(eph_obsy_code,dim=1) < i) THEN
                   eph_obsy_code => reallocate(eph_obsy_code, 2*i)
                END IF
                eph_obsy_code(i) = TRIM(par_val(indx1:indx1+indx2-1))
                ! Move starting index to the next parameter 
                indx1 = indx1 + indx2
                CALL removeLeadingBlanks(par_val(indx1:))
             END DO
             eph_obsy_code => reallocate(eph_obsy_code, i)
          END IF

       CASE ("eph.date")
          IF (PRESENT(eph_date)) THEN
             IF (LEN_TRIM(par_val) == 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "No dates specified for option 'eph.date'.", 1)
                RETURN
             END IF
             i = 0
             indx1 = 1
             last = .FALSE.
             CALL removeLeadingBlanks(par_val)
             DO WHILE (.NOT.last)
                ! Find the location of the next space: 
                indx2 = INDEX(TRIM(par_val(indx1:))," ")
                IF (indx2 == 0) THEN
                   ! No more spaces left -> last parameter
                   last = .TRUE.
                   indx2 = LEN_TRIM(par_val(indx1:))
                END IF
                i = i + 1
                IF (.NOT.ASSOCIATED(eph_date)) THEN
                   ALLOCATE(eph_date(i), stat=err)
                   IF (err /= 0) THEN
                      error = .TRUE.
                      CALL errorMessage("io / readConfigurationFile", &
                           "Could not allocate memory.", 1)
                      RETURN
                   END IF
                ELSE IF (SIZE(eph_date,dim=1) < i) THEN
                   eph_date => reallocate(eph_date, 2*i)
                END IF
                CALL toInt(TRIM(par_val(indx1:indx1+3)), year, error)
                CALL toInt(TRIM(par_val(indx1+4:indx1+5)), month, error)
                CALL toReal(TRIM(par_val(indx1+6:indx1+indx2-1)), day, error)
                IF (year >= 1972) THEN
                   CALL NEW(eph_date(i), year, month, day, "UTC")
                ELSE
                   CALL NEW(eph_date(i), year, month, day, "UT1")
                END IF
                ! Move starting index to the next parameter 
                indx1 = indx1 + indx2
                CALL removeLeadingBlanks(par_val(indx1:))
             END DO
             eph_date => reallocate(eph_date, i)
          END IF

       CASE ("eph.dt_since_last_obs")
          IF (PRESENT(eph_dt_since_last_obs)) THEN
             IF (LEN_TRIM(par_val) == 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "No dt's specified for option 'eph.dt_since_last_obs'.", 1)
                RETURN
             END IF
             i = 0
             indx1 = 1
             last = .FALSE.
             CALL removeLeadingBlanks(par_val)
             DO WHILE (.NOT.last)
                ! Find the location of the next space: 
                indx2 = INDEX(TRIM(par_val(indx1:))," ")
                IF (indx2 == 0) THEN
                   ! No more spaces left -> last parameter
                   last = .TRUE.
                   indx2 = LEN_TRIM(par_val(indx1:))
                END IF
                i = i + 1
                IF (.NOT.ASSOCIATED(eph_dt_since_last_obs)) THEN
                   ALLOCATE(eph_dt_since_last_obs(i), stat=err)
                   IF (err /= 0) THEN
                      error = .TRUE.
                      CALL errorMessage("io / readConfigurationFile", &
                           "Could not allocate memory.", 1)
                      RETURN
                   END IF
                ELSE IF (SIZE(eph_dt_since_last_obs,dim=1) < i) THEN
                   eph_dt_since_last_obs => reallocate(eph_dt_since_last_obs, 2*i)
                END IF
                CALL toReal(TRIM(par_val(indx1:indx1+indx2-1)), eph_dt_since_last_obs(i), error)
                ! Move starting index to the next parameter 
                indx1 = indx1 + indx2
                CALL removeLeadingBlanks(par_val(indx1:))
             END DO
             eph_dt_since_last_obs => reallocate(eph_dt_since_last_obs, i)
          END IF

       CASE ("eph.lt_correction")
          IF (PRESENT(eph_lt_correction)) THEN
             ! Check for missing parameters
             IF (LEN_TRIM(par_val) == 0) THEN
                error = .TRUE.
                CALL errorMessage("io / readConfigurationFile", &
                     "No values specified for option 'eph.lt_correction'.", 1)
                RETURN
             END IF
             i = 0
             indx1 = 1
             last = .FALSE.
             CALL removeLeadingBlanks(par_val)
             DO WHILE (.NOT.last)
                ! Find the location of the next space: 
                indx2 = INDEX(TRIM(par_val(indx1:))," ")
                IF (indx2 == 0) THEN
                   ! No more spaces left -> last parameter
                   last = .TRUE.
                   indx2 = LEN_TRIM(par_val(indx1:))
                END IF
                i = i + 1
                IF (.NOT.ASSOCIATED(eph_lt_correction)) THEN
                   ALLOCATE(eph_lt_correction(i), stat=err)
                   IF (err /= 0) THEN
                      error = .TRUE.
                      CALL errorMessage("io / readConfigurationFile", &
                           "Could not allocate memory.", 1)
                      RETURN
                   END IF
                ELSE IF (SIZE(eph_lt_correction,dim=1) < i) THEN
                   eph_lt_correction => reallocate(eph_lt_correction, 2*i)
                END IF
                READ(par_val(indx1:indx1+indx2-1), *, iostat=err) eph_lt_correction(i)
                IF (err /= 0) THEN
                   error = .TRUE.
                   CALL errorMessage("io / readConfigurationFile", &
                        "Could not read parameter value (30).", 1)
                   RETURN
                END IF
                ! Move starting index to the next parameter 
                indx1 = indx1 + indx2
                CALL removeLeadingBlanks(par_val(indx1:))
             END DO
             eph_lt_correction => reallocate(eph_lt_correction, i)
          END IF

       END SELECT
       IF (error) THEN
          CALL errorMessage("io / readConfigurationFile", &
               "Problems with the following line:", 1)
          IF (err_verb >= 1) THEN
             WRITE(stderr,"(A)") TRIM(line)
          END IF
          RETURN
       END IF
    END DO

  END SUBROUTINE readConfigurationFile





  SUBROUTINE readDESOrbitFile(lu, norb, header, id_arr, orb_arr, H_arr)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: lu
    INTEGER, INTENT(out) :: norb
    CHARACTER(len=*), INTENT(out) :: header
    CHARACTER(len=*), DIMENSION(:), INTENT(out) :: id_arr
    TYPE (Orbit), DIMENSION(:), INTENT(out) :: orb_arr
    REAL(bp), DIMENSION(:), INTENT(out) :: H_arr
    TYPE (Time) :: epoch
    CHARACTER(len=16) :: compcode, str
    CHARACTER(len=3) :: frmt
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp) :: mjd_epoch, moid
    INTEGER :: err, indx_, npar

    IF (LEN(header) < 128) THEN
       error = .TRUE.
       CALL errorMessage("io / readDESOrbitFile", &
            "Length of header variable must be at least 128 characters.", 1)
       RETURN       
    END IF

    ! Read header if not yet known:
    READ(lu, "(A)", iostat=err) header
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / readDESOrbitFile", &
            "Could not read 1st header line of orbit file.", 1)
       RETURN
    END IF
    IF (.NOT.(header(1:2) == "!!" .OR. header(1:1) == "#")) THEN
       REWIND(lu)
       header = " "
    END IF

    norb = 0
    DO WHILE (norb + 1 <= SIZE(id_arr))
       ! Read parameters found on a single line
       READ(lu, *, iostat=err) id_arr(norb+1), frmt, elements, &
            H_arr(norb+1), mjd_epoch, indx_, npar, moid, compcode
       IF (err > 0) THEN
          CALL toString(norb+1, str, error)
          error = .TRUE.
          CALL errorMessage("io / readDESOrbitFile", &
               "Problem reading line #" // TRIM(str) // " of the data file.", 1)
          RETURN
       ELSE IF (err < 0) THEN
          EXIT
       ELSE IF (id_arr(norb+1)(1:1) == "#") THEN
          CYCLE
       END IF
       CALL removeLeadingBlanks(id_arr(norb+1))
       CALL NEW(epoch, mjd_epoch, "TT")
       IF (error) THEN
          CALL errorMessage("io / readDESOrbitFile", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       IF (frmt == "COM") THEN
          elements(3:5) = elements(3:5)*rad_deg
          CALL NEW(orb_arr(norb+1), elements, "cometary", "ecliptic", epoch)
       ELSE IF (frmt == "COT") THEN
          elements(3:6) = elements(3:6)*rad_deg
          CALL NEW(orb_arr(norb+1), elements, "cometary_ta", "ecliptic", epoch)
       ELSE IF (frmt == "KEP") THEN
          elements(3:6) = elements(3:6)*rad_deg
          CALL NEW(orb_arr(norb+1), elements, "keplerian", "ecliptic", epoch)
       ELSE IF (frmt == "CAR") THEN
          CALL NEW(orb_arr(norb+1), elements, "cartesian", "ecliptic", epoch)
       ELSE
          error = .TRUE.
          CALL errorMessage("io / readDESOrbitFile", &
               "No such option available: " // TRIM(frmt), 1)
          RETURN          
       END IF
       IF (error) THEN
          CALL errorMessage("io / readDESOrbitFile", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       CALL NULLIFY(epoch)
       norb = norb + 1
    END DO

  END SUBROUTINE readDESOrbitFile





  SUBROUTINE readMPCOrbitFile(lu, norb, id_arr, orb_arr, HG_arr, arc_arr)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: lu
    INTEGER, INTENT(out) :: norb
    CHARACTER(len=*), DIMENSION(:), INTENT(out) :: id_arr
    TYPE (Orbit), DIMENSION(:), INTENT(out) :: orb_arr
    REAL(bp), DIMENSION(:,:), INTENT(out) :: HG_arr
    REAL(bp), DIMENSION(:), INTENT(out) :: arc_arr

    TYPE (Time) :: epoch
    CHARACTER(len=256) :: line
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp) :: day, n
    INTEGER :: err, y1, y2, year, month, nlines

    ! Jump over header:
    line = ""
    nlines = 0
    DO WHILE (line(1:5) /= "-----")
       READ(lu, "(A)", iostat=err) line
       IF (err > 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readMPCOrbitFile", &
               "Could not jump over header line of MPC orbit file.", 1)
          RETURN
       ELSE IF (err < 0 .AND. nlines > 0) THEN
          ! Probably no header present, let's try rewinding to
          ! beginning of file and read orbits:
          REWIND(lu)
          EXIT
       END IF
       nlines = nlines + 1
    END DO

    norb = 0
    DO WHILE (norb + 1 <= SIZE(id_arr))
       ! Read parameters found on a single line
       READ(lu, "(A)", iostat=err) line
       IF (err > 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readMPCOrbitFile", &
               "Read error.", 1)
          RETURN
       ELSE IF (err < 0) THEN
          EXIT
       ELSE IF (line(1:1) == "#") THEN
          CYCLE
       END IF
       norb = norb + 1
       ! Designation or number
       id_arr(norb) = line(1:7)
       ! H
       IF (LEN_TRIM(line(9:13)) == 0) THEN
          HG_arr(norb,1) = 99.9_bp
       ELSE
          CALL toReal(line(9:13), HG_arr(norb,1), error)
       END IF
       ! G
       IF (LEN_TRIM(line(15:19)) == 0) THEN
          HG_arr(norb,2) = 9.9_bp
       ELSE
          CALL toReal(line(15:19), HG_arr(norb,2), error)
       END IF
       IF (error) THEN
          CALL errorMessage("io / readMPCOrbitFile", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       ! Epoch
       CALL decodeMPCDate(line(21:25), year, month, day)
       CALL NEW(epoch, year, month, day, "TT")
       IF (error) THEN
          CALL errorMessage("io / readMPCOrbitFile", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       ! M
       CALL toReal(line(27:35), elements(6), error)
       ! ap
       CALL toReal(line(38:46), elements(5), error)
       ! an
       CALL toReal(line(49:57), elements(4), error)
       ! i
       CALL toReal(line(60:68), elements(3), error)
       ! e
       CALL toReal(line(71:79), elements(2), error)
       ! n
       CALL toReal(line(81:91), n, error)       
       ! a
       CALL toReal(line(93:103), elements(1), error)
       elements(3:6) = elements(3:6)*rad_deg
       CALL NEW(orb_arr(norb), elements, "keplerian", "ecliptic", epoch)
       IF (error) THEN
          CALL errorMessage("io / readMPCOrbitFile", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL NULLIFY(epoch)
       IF (line(132:132) == "-") THEN
          CALL toInt(line(128:131), y1, error)
          CALL toInt(line(133:136), y2, error)
          arc_arr(norb) = REAL(y2-y1)*day_year
       ELSE
          CALL toReal(line(128:131), arc_arr(norb), error)          
       END IF
       IF (error) THEN
          CALL errorMessage("io / readMPCOrbitFile", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
    END DO

  END SUBROUTINE readMPCOrbitFile





  SUBROUTINE readOpenOrbOrbitFile(lu, header, element_type_in, id, orb, &
       element_type_pdf, cov, pdf, rchi2, reg_apr, jac_sph_inv, &
       jac_car_kep, jac_equ_kep, H, G, rho1, rho2)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: lu
    CHARACTER(len=*), DIMENSION(4), INTENT(inout) :: header
    CHARACTER(len=*), INTENT(inout) :: element_type_in
    CHARACTER(len=*), INTENT(inout) :: id
    TYPE (Orbit), INTENT(inout) :: orb
    CHARACTER(len=*), INTENT(inout), OPTIONAL :: element_type_pdf
    REAL(bp), DIMENSION(:,:), INTENT(inout), OPTIONAL :: cov
    REAL(bp), INTENT(inout), OPTIONAL :: pdf, rchi2, &
         reg_apr, jac_sph_inv, jac_car_kep, jac_equ_kep, H, G, &
         rho1, rho2
    TYPE (Time) :: t
    CHARACTER(len=128) :: str
    REAL(bp), DIMENSION(15) :: correlation
    REAL(bp), DIMENSION(6) :: elements, stdev
    REAL(bp) :: day, r, mjd
    INTEGER :: year, month, err, i

    correlation = 2.0_bp

    IF (LEN(header(1)(:)) < 1024) THEN
       error = .TRUE.
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "Length of header variable must be at least 1024.", 1)
       RETURN       
    END IF

    ! Read header if not yet known:
    IF (LEN_TRIM(header(1)) == 0) THEN
       READ(lu, "(A)", iostat=err) header(1)(:)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Could not read 1st header line of orbit file.", 1)
          RETURN
       END IF
       READ(lu, "(A)", iostat=err) header(2)(:)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Could not read 2nd header line of orbit file.", 1)
          RETURN
       END IF
       READ(lu, "(A)", iostat=err) header(3)(:)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Could not read 3rd header line of orbit file.", 1)
          RETURN
       END IF
       READ(lu, "(A)", iostat=err) header(4)(:)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Could not read 4th header line of orbit file.", 1)
          RETURN
       END IF
    END IF

    ! Read id, elements, and epoch
    READ(lu,"(A16,6(1X,E21.14))", advance="no", iostat=err) id, &
         elements
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "Read error (5).", 1)
       RETURN
    END IF
    CALL removeLeadingBlanks(id)
    DO i=1,LEN(id)
       IF (IACHAR(id(i:i)) == 0) THEN
          id(i:i) = CHAR(32)
       END IF
    END DO
    IF (INDEX(header(4),"-0008-") /= 0) THEN
       READ(lu,"(1X,I4,1X,I2,1X,F8.5)", advance="no", iostat=err) year, &
            month, day
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (6).", 1)
          RETURN
       END IF
       CALL NEW(t, year, month, day, "TT")
    ELSE IF (INDEX(header(4),"-0074-") /= 0) THEN
       READ(lu,"(1X,F16.10)", advance="no", iostat=err) mjd
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (7).", 1)
          RETURN
       END IF
       CALL NEW(t, mjd, "TT")
    ELSE
       error = .TRUE.
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "Epoch not given or given in wrong format.", 1)
       RETURN       
    END IF
    IF (error) THEN
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    IF (INDEX(header(4),"-0002-") /= 0 .AND. &
         INDEX(header(4),"-0003-") /= 0 .AND. &
         INDEX(header(4),"-0004-") /= 0 .AND. &
         INDEX(header(4),"-0005-") /= 0 .AND. &
         INDEX(header(4),"-0006-") /= 0 .AND. &
         INDEX(header(4),"-0007-") /= 0) THEN
       element_type_in = "keplerian"
       elements(3:6) = elements(3:6)*rad_deg
    ELSE IF (INDEX(header(4),"-0039-") /= 0 .AND. &
         INDEX(header(4),"-0040-") /= 0 .AND. &
         INDEX(header(4),"-0041-") /= 0 .AND. &
         INDEX(header(4),"-0042-") /= 0 .AND. &
         INDEX(header(4),"-0043-") /= 0 .AND. &
         INDEX(header(4),"-0044-") /= 0) THEN
       element_type_in = "cartesian"
    ELSE
       error = .TRUE.
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "Orbital elements are unknown to the software.", 1)       
       RETURN
    END IF
    CALL NEW(orb, elements, element_type_in, "ecliptic", copy(t))
    IF (error) THEN
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    CALL NULLIFY(t)
    IF (INDEX(header(4),"-0038-") /= 0) THEN
       str = " "
       READ(lu,"(1X,A12)", advance="no", iostat=err) str(1:12)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (10).", 1)
          RETURN
       END IF
       IF (PRESENT(element_type_pdf)) THEN
          element_type_pdf = str(1:12)
       END IF
    END IF
    stdev = -1.0_bp
    IF (INDEX(header(4),"-0009-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) stdev(1)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (10).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0010-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) stdev(2)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (15).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0011-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) stdev(3)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (20).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0012-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) stdev(4)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (25).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0013-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) stdev(5)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (30).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0014-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) stdev(6)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (35).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0015-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(1)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (40).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0016-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(2)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (45).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0017-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(3)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (50).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0018-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(4)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (55).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0019-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(5)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (60).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0020-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(6)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (65).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0021-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(7)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (70).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0022-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(8)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (75).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0023-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(9)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (80).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0024-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(10)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (85).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0025-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(11)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (90).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0026-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(12)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (95).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0027-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(13)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (100).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0028-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(14)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (105).", 1)
          RETURN
       END IF
    END IF
    IF (INDEX(header(4),"-0029-") /= 0) THEN
       READ(lu,"(1X,E15.7)", advance="no", iostat=err) correlation(15)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (110).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(cov)) THEN
       IF (INDEX(header(4),"-0009-") /= 0) THEN
          ! Standard deviations appear to be present, so build the
          ! covariance using the available information.
          IF (element_type_in == "keplerian") THEN
             stdev(3:6) = stdev(3:6)*rad_deg
          END IF
          DO i=1,6
             cov(i,i) = stdev(i)**2
          END DO
          cov(1,2:6) = stdev(1)*stdev(2:6)*correlation(1:5)
          cov(2:6,1) = cov(1,2:6)
          cov(2,3:6) = stdev(2)*stdev(3:6)*correlation(6:9)
          cov(3:6,2) = cov(2,3:6)
          cov(3,4:6) = stdev(3)*stdev(4:6)*correlation(10:12)
          cov(4:6,3) = cov(3,4:6)
          cov(4,5:6) = stdev(4)*stdev(5:6)*correlation(13:14)
          cov(5:6,4) = cov(4,5:6)
          cov(5,6) = stdev(5)*stdev(6)*correlation(15)
          cov(6,5) = cov(5,6)
       END IF
    END IF
    IF (INDEX(header(4),"-0030-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (115).", 1)
          RETURN
       END IF
       IF (PRESENT(pdf)) THEN
          pdf = r
       END IF
    END IF
    IF (INDEX(header(4),"-0031-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (120).", 1)
          RETURN
       END IF
       IF (PRESENT(rchi2)) THEN
          rchi2 = r
       END IF
    END IF
    IF (INDEX(header(4),"-0032-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (125).", 1)
          RETURN
       END IF
       IF (PRESENT(reg_apr)) THEN
          reg_apr = r
       END IF
    END IF
    IF (INDEX(header(4),"-0033-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (130).", 1)
          RETURN
       END IF
       IF (PRESENT(jac_sph_inv)) THEN
          jac_sph_inv = r
       END IF
    END IF
    IF (INDEX(header(4),"-0034-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (135).", 1)
          RETURN
       END IF
       IF (PRESENT(jac_car_kep)) THEN
          jac_car_kep = r
       END IF
    END IF
    IF (INDEX(header(4),"-0035-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (140).", 1)
          RETURN
       END IF
       IF (PRESENT(jac_equ_kep)) THEN
          jac_equ_kep = r
       END IF
    END IF
    IF (INDEX(header(4),"-0036-") /= 0) THEN
       !READ(lu,"(1X,F9.6)", advance="no", iostat=err) r
       READ(lu, "(1X,A9)", advance="no", iostat=err) str(1:9)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (145).", 1)
          RETURN
       END IF
       IF (PRESENT(H)) THEN
          CALL toReal(str(1:9), H, error)
          IF (error) THEN
             CALL errorMessage("io / readOpenOrbOrbitFile", &
                  "Conversion error (145).", 1)
             RETURN
          END IF
          !          H = r
       END IF
    END IF
    IF (INDEX(header(4),"-0037-") /= 0) THEN
       READ(lu,"(1X,F9.6)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (150).", 1)
          RETURN
       END IF
       IF (PRESENT(G)) THEN
          G = r
       END IF
    END IF
    IF (INDEX(header(4),"-0070-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (150).", 1)
          RETURN
       END IF
       IF (PRESENT(rho1)) THEN
          rho1 = r
       END IF
    END IF
    IF (INDEX(header(4),"-0071-") /= 0) THEN
       READ(lu,"(1X,E18.10)", advance="no", iostat=err) r
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / readOpenOrbOrbitFile", &
               "Read error (155).", 1)
          RETURN
       END IF
       IF (PRESENT(rho2)) THEN
          rho2 = r
       END IF
    END IF
    READ(lu,"(1X)",iostat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / readOpenOrbOrbitFile", &
            "Read error (160).", 1)
       RETURN
    END IF

  END SUBROUTINE readOpenOrbOrbitFile





  SUBROUTINE writeDESOrbitFile(lu, print_header, element_type_out, &
       id, orb, H, indx, npar, moid, compcode)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: lu
    LOGICAL, INTENT(in) :: print_header
    CHARACTER(len=*), INTENT(in) :: element_type_out
    CHARACTER(len=*), INTENT(in) :: id
    TYPE (Orbit), INTENT(in) :: orb
    REAL(bp), INTENT(in) :: H
    INTEGER, INTENT(in), OPTIONAL :: indx, npar
    REAL(bp), INTENT(in), OPTIONAL :: moid
    CHARACTER(len=*), INTENT(in), OPTIONAL :: compcode

    TYPE (Time) :: epoch
    CHARACTER(len=32) :: compcode_
    CHARACTER(len=3) :: frmt
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp) :: mjd_tt, moid_
    INTEGER :: err, indx_, npar_

    elements = getElements(orb, element_type_out, "ecliptic")
    IF (error) THEN
       CALL errorMessage("io / writeDESOrbitFile", &
            "TRACE BACK", 1)
       RETURN       
    END IF
    IF (element_type_out == "cometary") THEN
       frmt = "COM"
       elements(3:5) = elements(3:5)/rad_deg
    ELSE IF (element_type_out == "keplerian") THEN
       frmt = "KEP"
       elements(3:6) = elements(3:6)/rad_deg
    ELSE IF (element_type_out == "cartesian") THEN
       frmt = "CAR"
    ELSE
       error = .TRUE.
       CALL errorMessage("io / writeDESOrbitFile", &
            "Element type '" // TRIM(element_type_out) // &
            "' not yet supported (5).", 1)
       RETURN
    END IF
    epoch = getTime(orb)
    mjd_tt = getMJD(epoch, "TT")
    CALL NULLIFY(epoch)
    IF (PRESENT(indx)) THEN
       indx_ = indx
    ELSE
       indx_ = 1
    END IF
    IF (PRESENT(npar)) THEN
       npar_ = npar
    ELSE
       npar_ = 6
    END IF
    IF (PRESENT(moid)) THEN
       moid_ = moid
    ELSE
       moid_ = -1.0_bp
    END IF
    IF (PRESENT(compcode)) THEN
       compcode_ = compcode
    ELSE
       compcode_ = "OPENORB"
    END IF
    IF (print_header) THEN
       SELECT CASE (element_type_out)
       CASE ("cometary")
          WRITE(lu,"(A)") "!!OID FORMAT q e i node argperi t_p H t_0 INDEX N_PAR MOID COMPCODE"
       CASE ("keplerian")
          WRITE(lu,"(A)") "!!OID FORMAT a e i node argperi M H t_0 INDEX N_PAR MOID COMPCODE"
       CASE ("cartesian")
          WRITE(lu,"(A)") "!!OID FORMAT x y z dx/dt dy/dt dz/dt H t_0 INDEX N_PAR MOID COMPCODE"
       CASE default
          error = .TRUE.
          CALL errorMessage("io / writeDESOrbitFile", &
               "Element type '" // TRIM(element_type_out) // &
               "' not yet supported (10).", 1)
          RETURN
       END SELECT
    END IF
    WRITE(lu,"(2(A,1X),8(E22.15,1X),2(I0,1X),E22.15,1X,A)",iostat=err) &
         TRIM(id), TRIM(frmt), elements, H, mjd_tt, indx_, npar_, &
         moid_, TRIM(compcode_)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeDESOrbitFile", &
            "Write error (5).", 1)
       RETURN
    END IF

  END SUBROUTINE writeDESOrbitFile





  SUBROUTINE writeNominalSolution(storb, obss, element_type, lu)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)            :: storb
    TYPE (Observations), INTENT(in)                  :: obss
    CHARACTER(len=*), INTENT(in)                     :: element_type
    INTEGER, INTENT(in)                              :: lu

    TYPE (Orbit) :: orb
    TYPE (Time) :: t
    CHARACTER(len=DESIGNATION_LEN) :: &
         id
    CHARACTER(len=ELEMENT_TYPE_LEN) :: &
         element_type_
    CHARACTER(len=DYN_MODEL_LEN) :: &    
         dyn_model
    CHARACTER(len=INTEGRATOR_LEN) :: &
         integrator
    CHARACTER(len=64) :: &
         str
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         inform_mat_obs_bd
    REAL(bp), DIMENSION(:,:), POINTER :: &
         residuals
    REAL(bp), DIMENSION(6,6)   :: cov, corr
    REAL(bp), DIMENSION(6)     :: elements, sigmas
    REAL(bp) :: obsarc, rchi2
    INTEGER :: k, l, err
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks
    LOGICAL, DIMENSION(6) :: ls_element_mask

    IF (.NOT. exist(storb)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeNominalSolution", &
            "StochasticOrbit object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeNominalSolution", &
            "Observations object has not yet been initialized.", 1)
       RETURN
    END IF

    orb = getNominalOrbit(storb)
    IF (error) THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (105)", 1)
       RETURN
    END IF
    t = getTime(orb)
    IF (error) THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (110)", 1)
       RETURN
    END IF

    element_type_ = TRIM(element_type)
    str = "#CAR  "
    elements = getElements(orb, TRIM(element_type_))
    IF (error .AND. TRIM(element_type_) == "cartesian") THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (111)", 1)
       RETURN       
    ELSE IF (error .AND. TRIM(element_type_) == "keplerian") THEN
       error = .FALSE.
       element_type_ = "cartesian"
       elements = getElements(orb, TRIM(element_type_), "ecliptic")
    ELSE IF (.NOT.error .AND. TRIM(element_type_) == "keplerian") THEN
       elements(3:6) = elements(3:6)/rad_deg
       str = "#KEP  "
    END IF
    id = getID(obss)
!!$    obsarc = getObservationalTimespan(obss)
!!$    WRITE(lu,"(A,F14.4)") "#Observational arc = ",obsarc

    CALL getParameters(orb, dyn_model=dyn_model, integrator=integrator)
    WRITE(lu,"(A,A)") "#DYNAMICAL MODEL = ", TRIM(dyn_model)
    IF (TRIM(dyn_model) == 'n-body') THEN
       WRITE(lu,"(A,A)") "#INTEGRATOR = ", TRIM(integrator)
    END IF
    WRITE(lu,"(A)") "#"

    IF (element_type_ == "keplerian") THEN
       WRITE(lu,"(A6,2X,A7,1X,7(3X,A12,2X))") str(1:6), id(1:7), & 
            "a [AU]", "e", "i [deg]", "node [deg]", "ap [deg]", & 
            "M [deg]", "Epoch"
    ELSE
       WRITE(lu,"(A6,2X,A7,1X,7(3X,A12,2X))") str(1:6), id(1:7), & 
            "x [AU]", "y [AU]", "z [AU]", "dx/dt [AU/d]", & 
            "dy/dt [AU/d]", "dz/dt [AU/d]", "Epoch" 
    END IF
    WRITE(lu,"(A6,2X,A7,1X,6(F16.12,1X),A,'TT')") &
         str(1:6), id(1:7), elements, getCalendarDateString(t,"tt")
    IF (error) THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (115)", 1)
       RETURN
    END IF
    CALL NULLIFY(t)
    WRITE(lu,"(A)") "#"
    ! STANDARD DEVIATIONS:
    WRITE(lu,"(A6,2X,A7,1X)",advance="no") "#STDEV", id(1:7)
    IF (TRIM(element_type_) == "keplerian") THEN
       cov = getCovarianceMatrix(storb, TRIM(element_type_))
    ELSE
       cov = getCovarianceMatrix(storb, TRIM(element_type_), "ecliptic")
    END IF
    IF (error) THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (116)", 1)
       RETURN
    END IF
    DO l=1,6
       sigmas(l) = SQRT(cov(l,l))
       IF (TRIM(element_type_) == "keplerian" .AND. l >= 3) THEN
          WRITE(lu,"(F16.12,1X)",advance="no") sigmas(l)/rad_deg
       ELSE
          WRITE(lu,"(F16.12,1X)",advance="no") sigmas(l)
       END IF
    END DO
    WRITE(lu,*)
    WRITE(lu,"(A)") "#"
    ! CORRELATIONS:
    DO l=1,6
       DO k=1,6
          corr(l,k) = cov(l,k) / &
               (sigmas(l)*sigmas(k))
       END DO
    END DO
    WRITE(lu,"(A6,2X,A7,6(1X,F16.12))") &
         "#CORR ", id(1:7), corr(1,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,F16.12))") &
         "#CORR ", id(1:7), corr(2,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,F16.12))") &
         "#CORR ", id(1:7), corr(3,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,F16.12))") &
         "#CORR ", id(1:7), corr(4,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,F16.12))") &
         "#CORR ", id(1:7), corr(5,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,F16.12))") &
         "#CORR ", id(1:7), corr(6,1:6)
    WRITE(lu,"(A)") "#"
    ! COVARIANCE:
    IF (TRIM(element_type_) == "keplerian") THEN
       cov(3:6,:) = cov(3:6,:)/rad_deg
       cov(:,3:6) = cov(:,3:6)/rad_deg
    END IF
    WRITE(lu,"(A6,2X,A7,6(1X,E16.8))") &
         "#COV  ", id(1:7), cov(1,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,E16.8))") &
         "#COV  ", id(1:7), cov(2,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,E16.8))") &
         "#COV  ", id(1:7), cov(3,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,E16.8))") &
         "#COV  ", id(1:7), cov(4,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,E16.8))") &
         "#COV  ", id(1:7), cov(5,1:6)
    WRITE(lu,"(A6,2X,A7,6(1X,E16.8))") &
         "#COV  ", id(1:7), cov(6,1:6)
    WRITE(lu,"(A)") "#"
    ! WRITE RESIDUAL BLOCK:
    obs_masks => getObservationMasks(storb)
    IF (error) THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (120)", 1)
       RETURN
    END IF
    WRITE(str,"(A6,2X,A7)") "#RES  ", id(1:7)
    CALL writeResidualBlock(orb, obss, obs_masks, &
         TRIM(str), lu, residuals)
    IF (error) THEN
       CALL errorMessage("io / writeNominalSolution", &
            "TRACE BACK (125)", 1)
       RETURN
    END IF
    WRITE(lu,"(A)") "#"
    ! WRITE RMS:
    WRITE(lu,"(A6,2X,A7,3X,2(F11.6,1X))") "#RMS  ", id(1:7), &
         SQRT(SUM(residuals(:,2)**2.0_bp,obs_masks(:,2) &
         .AND. obs_masks(:,3))/COUNT(obs_masks(:,2) &
         .AND. obs_masks(:,3)))/rad_asec, &
         SQRT(SUM(residuals(:,3)**2.0_bp,obs_masks(:,2) &
         .AND. obs_masks(:,3))/COUNT(obs_masks(:,2) &
         .AND. obs_masks(:,3)))/rad_asec
    WRITE(lu,"(A)") "#"
    ! WRITE REDUCED CHI2:
    CALL getParameters(storb, ls_element_mask=ls_element_mask)
    inform_mat_obs_bd => getBlockDiagInformationMatrix(obss)
    rchi2 = chi_square(residuals, inform_mat_obs_bd, obs_masks, errstr) / &
         REAL(COUNT(obs_masks)-COUNT(ls_element_mask),bp)
    IF (LEN_TRIM(errstr) /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeNominalSolution", &
            "Could not compute chi2. " // TRIM(errstr), 1)
       errstr = ""
       RETURN
    END IF
    WRITE(lu,"(A6,2X,A7,1X,F14.6)") "#RCHI2", id(1:7), rchi2
    WRITE(lu,"(A)") "#"
    DEALLOCATE(inform_mat_obs_bd, stat=err)
    ! WRITE OBSERVATIONAL TIMESPAN:
    obsarc = getObservationalTimespan(obss)
    WRITE(lu,"(A7,1X,A7,1X,F14.4)") "#OBSARC", id(1:7), obsarc
    WRITE(lu,*)

    DEALLOCATE(obs_masks, stat=err)
    DEALLOCATE(residuals, stat=err)
    CALL NULLIFY(orb)

  END SUBROUTINE writeNominalSolution





  SUBROUTINE writeOpenOrbOrbitFile(lu, print_header, element_type_out, &
       id, orb, element_type_pdf, cov, pdf, rchi2, reg_apr, &
       jac_sph_inv, jac_car_kep, jac_equ_kep, H, G, rho1, rho2, &
       rms, npdf, mjd)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: lu
    LOGICAL, INTENT(in) :: print_header
    CHARACTER(len=*), INTENT(in) :: element_type_out
    CHARACTER(len=*), INTENT(in) :: id
    TYPE (Orbit), INTENT(in) :: orb
    CHARACTER(len=*), INTENT(in), OPTIONAL :: element_type_pdf
    REAL(bp), DIMENSION(:,:), INTENT(in), OPTIONAL :: cov
    REAL(bp), INTENT(in), OPTIONAL :: pdf, rchi2, reg_apr, &
         jac_sph_inv, jac_car_kep, jac_equ_kep, H, G, rho1, rho2, &
         rms, npdf
    LOGICAL, INTENT(in), OPTIONAL :: mjd

    TYPE (Time) :: t
    CHARACTER(len=1024), DIMENSION(4) :: header
    REAL(bp), DIMENSION(6) :: elements, stdev
    REAL(bp) :: day, jac, mjd_tt
    INTEGER :: year, month, err, indx, i, j

    ! Define format for output and construct header:
    IF (print_header) THEN
       indx = 1
       header(1:4)(1:LEN(header(1))) = " "
       header(1)(indx:indx+16) = "#    Number     "
       header(2)(indx:indx+16) = "#      or       "
       header(3)(indx:indx+16) = "#  designation  "
       header(4)(indx:indx+16) = "#-----0001-----<"
       indx = indx + 16
!!$       header(1)(17:32) = "   Periapsis    "
!!$       header(2)(17:32) = "   Distance q   "
!!$       header(3)(17:32) = "      [AU]      "
!!$       header(4)(17:32) = ">-----0002-----<"
       IF (element_type_out == "keplerian") THEN
          header(1)(indx:indx+22) = "   Semimajor axis a   "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [AU]         "
          header(4)(indx:indx+22) = ">--------0002--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "    Eccentricity e    "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">---------0003-------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "     Inclination i    "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0004--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "     Longitude of     "
          header(2)(indx:indx+22) = "   Ascending Node O   "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0005--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "      Argument of     "
          header(2)(indx:indx+22) = "      Periapsis o     "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0006--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "     Mean Anomaly M   "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0007--------<"
          indx = indx + 22
       ELSE IF (element_type_out == "cartesian") THEN
          header(1)(indx:indx+22) = "      Ecliptic x      "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [AU]         "
          header(4)(indx:indx+22) = ">--------0039--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "      Ecliptic y      "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [AU]         "
          header(4)(indx:indx+22) = ">--------0040--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "      Ecliptic z      "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [AU]         "
          header(4)(indx:indx+22) = ">--------0041--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "    Ecliptic dx/dt    "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "        [AU/d]        "
          header(4)(indx:indx+22) = ">--------0042--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "    Ecliptic dy/dt    "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "        [AU/d]        "
          header(4)(indx:indx+22) = ">--------0043--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "    Ecliptic dz/dt    "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "        [AU/d]        "
          header(4)(indx:indx+22) = ">--------0044--------<"
          indx = indx + 22
       ELSE IF (element_type_out == "delaunay") THEN
          !!   l      = Mean Motion (Mean Anomaly)
          !!   g      = Argument of Periapsis
          !!   &theta = Longitude of the Ascending Node
          !!   L      = sqrt(mu*(Semimajor Axis)) = related to the 
          !!                                        two-body orbital energy
          !!   G      = L*sqrt(1-Eccentricity^2)  = magnitude of the orbital 
          !!                                        angular momentum
          !!   &Theta = G*cos(Inclination)        = z-component of the orbital 
          !!                                        angular momentum
          header(1)(indx:indx+22) = "           l          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0045--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           g          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0046--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "         theta        "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0047--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           L          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0048--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           G          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0049--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "         Theta        "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0050--------<"
          indx = indx + 22
       ELSE IF (element_type_out == "poincare") THEN
          !! Output Poincaré elements calculated from Delaunay's
          !! elements. The mass of the target body is assumed to be negligible
          !! compared to the mass of the Sun.
          !  (&radic = square root)
          !!
          !! &Lambda = L
          !! &xi     = &radic(2(L - G)) * cos(g + &theta)
          !! p       = &radic(2(G - &Theta)) * cos(&theta)
          !! &lambda = l + g + &theta
          !! &eta    = -&radic(2(L - G)) * sin(g + &theta)
          !! q       = -&radic(2(G - &Theta)) * sin(&theta)
          header(1)(indx:indx+22) = "        Lambda        "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0051--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "          xi          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0052--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           p          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0053--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "        lambda        "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0054--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "          eta         "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0055--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           q          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0056--------<"
          indx = indx + 22
       ELSE IF (element_type_out == "equinoctial") THEN
          !! Output equinoctial elements.
          header(1)(indx:indx+22) = "   Semimajor axis a   "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "         [AU]         "
          header(4)(indx:indx+22) = ">--------0057--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           h          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0058--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           k          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0059--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           p          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0060--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "           q          "
          header(2)(indx:indx+22) = "                      "
          header(3)(indx:indx+22) = "                      "
          header(4)(indx:indx+22) = ">--------0061--------<"
          indx = indx + 22
          header(1)(indx:indx+22) = "    Mean longitude    "
          header(2)(indx:indx+22) = "        lambda        "
          header(3)(indx:indx+22) = "         [deg]        "
          header(4)(indx:indx+22) = ">--------0062--------<"
          indx = indx + 22
       ELSE
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Output orbital elements are unknown to the software.", 1)       
          RETURN          
       END IF
       IF (PRESENT(mjd)) THEN
          IF (mjd) THEN
             header(1)(indx:indx+17) = "    Epoch (TT)   "
             header(2)(indx:indx+17) = "       MJD       "
             header(3)(indx:indx+17) = "                 "
             header(4)(indx:indx+17) = ">------0074-----<"
             indx = indx + 17             
          ELSE
             header(1)(indx:indx+17) = "    Epoch (TT)   "
             header(2)(indx:indx+17) = " YYYY MM DD.ddddd"
             header(3)(indx:indx+17) = "                 "
             header(4)(indx:indx+17) = ">------0008-----<"
             indx = indx + 17
          END IF
       ELSE
          header(1)(indx:indx+17) = "    Epoch (TT)   "
          header(2)(indx:indx+17) = " YYYY MM DD.ddddd"
          header(3)(indx:indx+17) = "                 "
          header(4)(indx:indx+17) = ">------0008-----<"
          indx = indx + 17
       END IF
       IF (PRESENT(element_type_pdf)) THEN
          header(1)(indx:indx+13) = "  Inversion  "
          header(2)(indx:indx+13) = "   orbital   "
          header(3)(indx:indx+13) = "  elements   "
          header(4)(indx:indx+13) = ">----0038---<"
          indx = indx + 13
       END IF
       IF (PRESENT(cov)) THEN
          header(1)(indx:indx+16) = "    sigma e1    "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0009-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "    sigma e2    "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0011-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "    sigma e3    "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0010-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "    sigma e4    "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0012-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "    sigma e5    "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0013-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "    sigma e6    "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0014-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e1,e2)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0015-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e1,e3)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0016-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e1,e4)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0017-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e1,e5)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0018-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e1,e6)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0019-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e2,e3)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0020-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e2,e4)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0021-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e2,e5)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0022-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e2,e6)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0023-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e3,e4)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0024-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e3,e5)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0025-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e3,e6)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0026-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e4,e5)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0027-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e4,e6)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0028-----<"
          indx = indx + 16
          header(1)(indx:indx+16) = "   cor(e5,e6)   "
          header(2)(indx:indx+16) = "                "
          header(3)(indx:indx+16) = "                "
          header(4)(indx:indx+16) = ">-----0029-----<"
          indx = indx + 16
       END IF
       IF (PRESENT(pdf)) THEN
          header(1)(indx:indx+19) = "    Unnormalized   "
          header(2)(indx:indx+19) = "       p.d.f       "
          header(3)(indx:indx+19) = "                   "
          header(4)(indx:indx+19) = ">-------0030------<"
          indx = indx + 19
       END IF
       IF (PRESENT(rchi2)) THEN
          header(1)(indx:indx+19) = "    Reduced chi2   "
          header(2)(indx:indx+19) = "                   "
          header(3)(indx:indx+19) = "                   "
          header(4)(indx:indx+19) = ">-------0031------<"
          indx = indx + 19
       END IF
       IF (PRESENT(reg_apr)) THEN
          header(1)(indx:indx+19) = "    Regularizing   "
          header(2)(indx:indx+19) = "       apriori     "
          header(3)(indx:indx+19) = "                   "
          header(4)(indx:indx+19) = ">-------0032------<"
          indx = indx + 19
       END IF
       IF (PRESENT(jac_sph_inv)) THEN
          header(1)(indx:indx+19) = "    Determinant    "
          header(2)(indx:indx+19) = "    of Jacobian    "
          header(3)(indx:indx+19) = "    (sph/inv.elm.) "
          header(4)(indx:indx+19) = ">-------0033------<"
          indx = indx + 19
       END IF
       IF (PRESENT(jac_car_kep)) THEN
          header(1)(indx:indx+19) = "    Determinant    "
          header(2)(indx:indx+19) = "    of Jacobian    "
          header(3)(indx:indx+19) = "     (car/kep)     "
          header(4)(indx:indx+19) = ">-------0034------<"
          indx = indx + 19
       END IF
       IF (PRESENT(jac_equ_kep)) THEN
          header(1)(indx:indx+19) = "    Determinant    "
          header(2)(indx:indx+19) = "    of Jacobian    "
          header(3)(indx:indx+19) = "     (equ/kep)     "
          header(4)(indx:indx+19) = ">-------0035------<"
          indx = indx + 19
       END IF
       IF (PRESENT(H)) THEN
          IF (H < 99.0_bp) THEN
             header(1)(indx:indx+10) = " Absolute "
             Header(2)(indx:indx+10) = " magnitude"
             header(3)(indx:indx+10) = "     H    "
             header(4)(indx:indx+10) = ">---0036-<"
             indx = indx + 10
          END IF
       END IF
       IF (PRESENT(G)) THEN
          IF (G < 9.0_bp) THEN
             header(1)(indx:indx+10) = "  Slope   "
             header(2)(indx:indx+10) = " parameter"
             header(3)(indx:indx+10) = "     G    "
             header(4)(indx:indx+10) = ">---0037-<"
             indx = indx + 10
          END IF
       END IF
       IF (PRESENT(rho1)) THEN
          header(1)(indx:indx+19) = "    Distance at    "
          header(2)(indx:indx+19) = "    first epoch    "
          header(3)(indx:indx+19) = "    in Ranging     "
          header(4)(indx:indx+19) = ">-------0070------<"
          indx = indx + 19
       END IF
       IF (PRESENT(rho2)) THEN
          header(1)(indx:indx+19) = "    Distance at    "
          header(2)(indx:indx+19) = "    second epoch   "
          header(3)(indx:indx+19) = "    in Ranging     "
          header(4)(indx:indx+19) = ">-------0071------<"
          indx = indx + 19
       END IF
       IF (PRESENT(rms)) THEN
          header(1)(indx:indx+13) = "  O-C resid. "
          header(2)(indx:indx+13) = "     RMS     "
          header(3)(indx:indx+13) = "   [arcsec]  "
          header(4)(indx:indx+13) = ">----0072---<"
          indx = indx + 13
       END IF
       IF (PRESENT(npdf)) THEN
          header(1)(indx:indx+19) = "     Normalized    "
          header(2)(indx:indx+19) = "       p.d.f.      "
          header(3)(indx:indx+19) = "                   "
          header(4)(indx:indx+19) = ">-------0073------<"
          indx = indx + 19
       END IF

       WRITE(lu, "(A)", iostat=err) TRIM(header(1))
       WRITE(lu, "(A)", iostat=err) TRIM(header(2))
       WRITE(lu, "(A)", iostat=err) TRIM(header(3))
       WRITE(lu, "(A)", iostat=err) TRIM(header(4))
    END IF
    elements = getElements(orb, element_type_out, frame="ecliptic")
    IF (error) THEN
       CALL errorMessage("io / writeOpenOrbOrbitFile", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
!!$    ! a -> q
!!$    elements(1) = elements(1)*abs(1-elements(2))
    IF (element_type_out == "keplerian") THEN
       elements(3:6) = elements(3:6)/rad_deg
    ELSE IF (element_type_out == "delaunay") THEN
       elements(1:3) = elements(1:3)/rad_deg
    ELSE IF (element_type_out == "poincare") THEN
       elements(4) = elements(4)/rad_deg
    ELSE IF (element_type_out == "equinoctial") THEN
       elements(6) = elements(6)/rad_deg
    END IF
    t = getTime(orb)
    IF (error) THEN
       CALL errorMessage("io / writeOpenOrbOrbitFile", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    IF (PRESENT(mjd)) THEN
       IF (mjd) THEN
          mjd_tt = getMJD(t, "TT")
          WRITE(lu, "(A16,6(1X,E21.14),1X,F16.8)", &
               advance="no", iostat=err) id, elements(1:6), mjd_tt
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / writeOpenOrbOrbitFile", &
                  "Write error (4).", 1)
             RETURN
          END IF
          CALL NULLIFY(t)
       END IF
    END IF
    IF (exist(t)) THEN ! MJD not requested...
       CALL getCalendarDate(t, "tdt", year, month, day)
       IF (error) THEN
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL NULLIFY(t)
       WRITE(lu, "(A16,6(1X,E21.14),1X,I4,1X,I2,1X,F8.5)", &
            advance="no", iostat=err) id, elements(1:6), year, month, &
            day
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (5).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(element_type_pdf)) THEN
       IF (LEN_TRIM(element_type_pdf) > 12) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Type description of orbital elements used in inversion too long.", 1)
          RETURN
       END IF
       WRITE(lu, "(1X,A12)", advance="no", iostat=err) element_type_pdf
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (20).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(cov)) THEN
       stdev = (/ (SQRT(cov(i,i)), i=1,6) /)
       IF (element_type_out == "keplerian") THEN
          WRITE(lu, "(6(1X,E15.7))", advance="no", iostat=err) &
               stdev(1:2), stdev(3:6)/rad_deg
       ELSE
          WRITE(lu, "(6(1X,E15.7))", advance="no", iostat=err) &
               stdev(1:6)
       END IF
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (10).", 1)
          RETURN
       END IF
       DO i=1,5
          DO j=i+1,6
             IF (stdev(i) == 0.0_bp .OR. stdev(j) == 0.0_bp) THEN
                WRITE(lu, "(1X,E15.7)", advance="no", iostat=err) 0.0_bp
             ELSE                
                WRITE(lu, "(1X,E15.7)", advance="no", iostat=err) cov(i,j)/(stdev(i)*stdev(j))
             END IF
          END DO
       END DO
    END IF
    IF (PRESENT(pdf)) THEN
       jac = 1.0_bp
       IF (PRESENT(element_type_pdf)) THEN
          IF (element_type_pdf == "cartesian" .AND. &
               element_type_out == "keplerian") THEN
             IF (PRESENT(jac_car_kep)) THEN
                jac = jac_car_kep
             ELSE
                error = .TRUE.
                CALL errorMessage("io / writeOpenOrbOrbitFile", &
                     "Jacobian between Cartesian and Keplerian " // &
                     "elements needed, but missimg.", 1)
                RETURN
             END IF
          ELSE IF (element_type_pdf == "keplerian" .AND. &
               element_type_out == "cartesian") THEN
             IF (PRESENT(jac_car_kep)) THEN
                jac = 1.0_bp/jac_car_kep
             ELSE
                error = .TRUE.
                CALL errorMessage("io / writeOpenOrbOrbitFile", &
                     "Jacobian between Cartesian and Keplerian " // &
                     "elements needed, but missimg.", 1)
                RETURN
             END IF
          END IF
          WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) pdf*jac
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / writeOpenOrbOrbitFile", &
                  "Write error (15).", 1)
             RETURN
          END IF
       ELSE          
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Missing type of orbital elements corresponding to p.d.f.", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(rchi2)) THEN
       WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) rchi2
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (20).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(reg_apr)) THEN
       WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) reg_apr
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (25).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(jac_sph_inv)) THEN
       WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) jac_sph_inv
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (30).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(jac_car_kep)) THEN
       WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) jac_car_kep
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (35).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(jac_equ_kep)) THEN
       WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) jac_equ_kep
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (40).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(H)) THEN
       IF (H < 99.0_bp) THEN     
          WRITE(lu, "(1X,F9.5)", advance="no", iostat=err) H
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / writeOpenOrbOrbitFile", &
                  "Write error (45).", 1)
             RETURN
          END IF
       END IF
    END IF
    IF (PRESENT(G)) THEN
       IF (G < 9.0_bp) THEN
          WRITE(lu, "(1X,F9.6)", advance="no", iostat=err) G
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / writeOpenOrbOrbitFile", &
                  "Write error (50).", 1)
             RETURN
          END IF
       END IF
    END IF
    IF (PRESENT(rho1)) THEN
       WRITE(lu, '(1X,E18.10)', advance='no', iostat=err) rho1
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (55).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(rho2)) THEN
       WRITE(lu, '(1X,E18.10)', advance='no', iostat=err) rho2
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (60).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(rms)) THEN
       WRITE(lu, '(1X,E12.4)', advance='no', iostat=err) rms/rad_asec
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (60).", 1)
          RETURN
       END IF
    END IF
    IF (PRESENT(npdf)) THEN
       WRITE(lu, "(1X,E18.10)", advance="no", iostat=err) npdf
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeOpenOrbOrbitFile", &
               "Write error (60).", 1)
          RETURN
       END IF
    END IF
    WRITE(lu,"(A)") ""

  END SUBROUTINE writeOpenOrbOrbitFile





  SUBROUTINE writeProbabilities(storb, lu, apriori)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)          :: storb
    INTEGER, INTENT(in)                            :: lu
    REAL(bp), DIMENSION(:,:), INTENT(in), OPTIONAL :: apriori

    CHARACTER(len=16), DIMENSION(:), POINTER :: group_arr
    REAL(bp), DIMENSION(:), POINTER :: probability_arr, apriori_pdf
    INTEGER :: err

    IF (PRESENT(apriori)) THEN
       CALL getAPrioriWeights(storb, apriori, apriori_pdf)
       IF (error) THEN
          CALL errorMessage("io / writeProbabilities", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       CALL getGroupWeights(storb, probability_arr, group_arr, apriori_pdf)
    ELSE
       CALL getGroupWeights(storb, probability_arr, group_arr)
    END IF
    IF (error) THEN
       CALL errorMessage("io / writeProbabilities", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    WRITE(lu,"('#')")
    WRITE(lu,"('#',3X,A)") &
         "Probabilities of belonging to one of the following groups"
    IF (PRESENT(apriori)) THEN
       WRITE(lu,"('#',3X,A)") &
            "(in %, using a priori information):" 
    ELSE
       WRITE(lu,"('#',3X,A)") &
            "(in %, without using a priori information):" 
    END IF
    probability_arr = 100.0_bp*probability_arr
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(1), probability_arr(1)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(2), probability_arr(2)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(3), probability_arr(3)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(4), probability_arr(4)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(5), probability_arr(5)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(6), probability_arr(6)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(7), probability_arr(7)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(8), probability_arr(8)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(9), probability_arr(9)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(10), probability_arr(10)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(11), probability_arr(11)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(12), probability_arr(12)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(13), probability_arr(13)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(14), probability_arr(14)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(15), probability_arr(15)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(16), probability_arr(16)
    WRITE(lu,"('#',5X,A16,2X,'=',2X,F9.4)") group_arr(17), probability_arr(17)
    WRITE(lu,"('#')")
    WRITE(lu,"('#',3X,A,2X,'=',2X,F9.5)") &
         "PHA probability (without prior knowledge and size requirement)", &
         100.0_bp*getPHAProbability(storb)
    IF (error) THEN
       CALL errorMessage("io / writeProbabilities", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    DEALLOCATE(probability_arr, group_arr, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeProbabilities", &
            "TRACE BACK (20)", 1)
       RETURN
    END IF

  END SUBROUTINE writeProbabilities





  SUBROUTINE writeResidualBlock(orb, obss, obs_mask, str, lu, residuals)

    IMPLICIT NONE
    TYPE (Orbit), INTENT(in) :: orb
    TYPE (Observations), INTENT(in) :: obss
    CHARACTER(len=*) :: str
    INTEGER, INTENT(in) :: lu
    LOGICAL, DIMENSION(:,:), INTENT(in) :: obs_mask
    REAL(bp), DIMENSION(:,:), POINTER, OPTIONAL :: residuals

    TYPE (Time) :: t
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: observed_scoords, &
         computed_scoords
    TYPE (Observation), DIMENSION(:), POINTER :: obs_arr
    CHARACTER(len=64), DIMENSION(:,:), ALLOCATABLE :: residual_block
    CHARACTER(len=64), DIMENSION(:), ALLOCATABLE :: residual_arr 
    CHARACTER(len=2) :: month_str, day_str
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: residuals_, observed_coords, &
         computed_coords
    REAL(bp) :: day
    INTEGER :: i, j, nobs, year, month, err

    ! Compute residuals
    observed_scoords => getObservationSCoords(obss)
    IF (error) THEN
       CALL errorMessage("io / writeResidualBlock", &
            "TRACE BACK (5)", 1)
       RETURN
    END IF
    obsy_ccoords => getObservatoryCCoords(obss)
    IF (error) THEN
       CALL errorMessage("io / writeResidualBlock", &
            "TRACE BACK (10)", 1)
       RETURN
    END IF
    nobs = SIZE(observed_scoords,dim=1)
    ALLOCATE(residuals_(nobs,6), observed_coords(nobs,6), &
         computed_coords(nobs,6), residual_arr(CEILING(nobs/3.0)*3), &
         residual_block(CEILING(nobs/3.0),3), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResidualBlock", &
            "Could not allocate memory (5).", 1)
       RETURN
    END IF
    CALL getEphemerides(orb, obsy_ccoords, computed_scoords)
    IF (error) THEN
       CALL errorMessage("io / writeResidualBlock", &
            "TRACE BACK (15)", 1)
       RETURN
    END IF
    observed_coords = 0.0_bp
    computed_coords = 0.0_bp
    DO i=1,nobs
       observed_coords(i,:) = getCoordinates(observed_scoords(i))
       IF (error) THEN
          CALL errorMessage("io / writeResidualBlock", &
               "TRACE BACK (20)", 1)
          RETURN
       END IF
       computed_coords(i,:) = getCoordinates(computed_scoords(i))
       IF (error) THEN
          CALL errorMessage("io / writeResidualBlock", &
               "TRACE BACK (25)", 1)
          RETURN
       END IF
    END DO
    residuals_(1:nobs,1:6) = observed_coords(1:nobs,1:6) - &
         computed_coords(1:nobs,1:6)        
    residuals_(1:nobs,2) = residuals_(1:nobs,2) * &
         COS(observed_coords(1:nobs,3))
    obs_arr => getObservations(obss)
    IF (error) THEN
       CALL errorMessage("io / writeResidualBlock", &
            "TRACE BACK (30)", 1)
       RETURN
    END IF
    residual_arr = " "
    DO i=1,nobs
       t = getTime(obs_arr(i))
       IF (error) THEN
          CALL errorMessage("io / writeResidualBlock", &
               "TRACE BACK (35)", 1)
          RETURN
       END IF
       CALL getCalendarDate(t, "TT", year, month, day)
       IF (year >= 1972) THEN
          CALL getCalendarDate(t, "UTC", year, month, day)
       ELSE
          CALL getCalendarDate(t, "UT1", year, month, day)
       END IF
       IF (error) THEN
          CALL errorMessage("io / writeResidualBlock", &
               "TRACE BACK (40)", 1)
          RETURN
       END IF
       CALL NULLIFY(t)
       CALL toString(month, month_str, error)
       IF (error) THEN
          CALL errorMessage("io / writeResidualBlock", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       IF (LEN_TRIM(month_str) == 1) THEN
          month_str = "0" // TRIM(month_str)
       END IF
       CALL toString(FLOOR(day), day_str, error)
       IF (error) THEN
          CALL errorMessage("io / writeResidualBlock", &
               "TRACE BACK (50)", 1)
          RETURN
       END IF
       IF (LEN_TRIM(day_str) == 1) THEN
          day_str = "0" // TRIM(day_str)
       END IF
       IF (ALL(obs_mask(i,2:3))) THEN
          WRITE(residual_arr(i), &
               "(I4,2(A2),2X,A4,2X,2(F9.5,1X))", iostat=err) &
               year, month_str, day_str, &
               getCode(obs_arr(i)), &
               residuals_(i,2:3)/rad_asec
       ELSE
          WRITE(residual_arr(i), &
               "(I4,2(A2),2X,A4,1X,A1,2(F9.5,1X),A1)",iostat=err) &
               year, month_str, day_str, &
               getCode(obs_arr(i)), "(", &
               residuals_(i,2:3)/rad_asec, ")"
       END IF
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResidualBlock", &
               "Write error (5).", 1)
          RETURN
       END IF
    END DO
    residual_block = " "
    residual_block = RESHAPE(residual_arr, &
         SHAPE(residual_block))
    DO i=1,SIZE(residual_block,dim=1)
       WRITE(lu,"(A,2X)",advance="no",iostat=err) TRIM(str)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResidualBlock", &
               "Write error (10).", 1)
          RETURN
       END IF
       DO j=1,SIZE(residual_block,dim=2)
          WRITE(lu,"(A,2X)",advance="no",iostat=err) &
               residual_block(i,j)(1:37)
          IF (err /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / writeResidualBlock", &
                  "Write error (15).", 1)
             RETURN
          END IF
       END DO
       WRITE(lu,*,iostat=err)
    END DO

    IF (PRESENT(residuals)) THEN
       ALLOCATE(residuals(nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResidualBlock", &
               "Could not allocate memory (10).", 1)
          RETURN
       END IF
       residuals = residuals_
    END IF
    DEALLOCATE(observed_coords, stat=err)
    DEALLOCATE(computed_coords, stat=err)
    DEALLOCATE(observed_scoords, stat=err)
    DEALLOCATE(obsy_ccoords, stat=err)
    DEALLOCATE(computed_scoords, stat=err) 
    DEALLOCATE(obs_arr, stat=err)
    DEALLOCATE(residual_arr, stat=err)
    DEALLOCATE(residual_block, stat=err)
    DEALLOCATE(residuals_, stat=err)

  END SUBROUTINE writeResidualBlock





  SUBROUTINE writeResiduals(storb, obss, lu, residuals, compute)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(in) :: storb
    TYPE (Observations), INTENT(in) :: obss
    INTEGER, INTENT(in) :: lu
    LOGICAL, INTENT(in), OPTIONAL :: compute
    REAL(bp), DIMENSION(:,:,:), POINTER, OPTIONAL :: residuals

    TYPE (Orbit), DIMENSION(:), POINTER  :: orb_arr
    TYPE (Time) :: t
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER :: obsy_ccoords
    TYPE (SphericalCoordinates), DIMENSION(:), POINTER :: observed_scoords, &
         computed_scoords
    TYPE (Observation), DIMENSION(:), POINTER :: obs_arr
    CHARACTER(len=OBSY_CODE_LEN), DIMENSION(:), ALLOCATABLE :: codes
    CHARACTER(len=DESIGNATION_LEN) :: id
    CHARACTER(len=10) :: str
    REAL(bp), DIMENSION(:,:,:), POINTER :: residuals_
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: observed_coords, &
         computed_coords
    REAL(bp), DIMENSION(:), ALLOCATABLE :: jds
    REAL(bp) :: day
    INTEGER :: i, j, nobs, err, year, month
    LOGICAL :: compute_

    ! NOTE: info on obs_mask is not output, should it be?.

    IF (.NOT. exist(storb)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResiduals", &
            "StochasticOrbit object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResiduals", &
            "Observations object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(compute)) THEN
       compute_ = compute
    ELSE
       compute_ = .FALSE.
    END IF

    IF (.NOT. compute_) THEN
       residuals_ => getResidualDistribution(storb)
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "Residuals are not available.", 1)
          RETURN
       END IF
    ELSE
       ! Compute residuals
       observed_scoords => getObservationSCoords(obss)
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       obsy_ccoords => getObservatoryCCoords(obss)
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "TRACE BACK (10)", 1)
          RETURN
       END IF
       nobs = SIZE(observed_scoords,dim=1)

       orb_arr => getSampleOrbits(storb)
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "TRACE BACK (15)", 1)
          RETURN
       END IF

       ALLOCATE(residuals_(SIZE(orb_arr,dim=1),nobs,6), &
            observed_coords(nobs,6), computed_coords(nobs,6), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResiduals", &
               "Could not allocate memory (5).", 1)
          RETURN
       END IF

       DO j=1,SIZE(orb_arr, dim=1)
          IF (ASSOCIATED(computed_scoords)) THEN
             NULLIFY(computed_scoords)
          END IF
          CALL getEphemerides(orb_arr(j), obsy_ccoords, &
               computed_scoords)
          IF (error) THEN
             CALL errorMessage("io / writeResiduals", &
                  "TRACE BACK (20)", 1)
             RETURN
          END IF
          observed_coords = 0.0_bp
          computed_coords = 0.0_bp
          DO i=1,nobs
             observed_coords(i,:) = getCoordinates(observed_scoords(i))
             IF (error) THEN
                CALL errorMessage("io / writeResiduals", &
                     "TRACE BACK (25)", 1)
                RETURN
             END IF
             computed_coords(i,:) = getCoordinates(computed_scoords(i))
             IF (error) THEN
                CALL errorMessage("io / writeResiduals", &
                     "TRACE BACK (30)", 1)
                RETURN
             END IF
          END DO
          residuals_(j,1:nobs,1:6) = observed_coords(1:nobs,1:6) - &
               computed_coords(1:nobs,1:6)        
          residuals_(j,1:nobs,2) = residuals_(j,1:nobs,2) * &
               COS(observed_coords(1:nobs,3))
       END DO

       DEALLOCATE(observed_coords, computed_coords, observed_scoords, &
            obsy_ccoords, computed_scoords, orb_arr, stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResiduals", &
               "Could not deallocate memory (5).", 1)
          DEALLOCATE(observed_coords, stat=err)
          DEALLOCATE(computed_coords, stat=err)
          DEALLOCATE(observed_scoords, stat=err)
          DEALLOCATE(obsy_ccoords, stat=err)
          DEALLOCATE(computed_scoords, stat=err)
          DEALLOCATE(orb_arr, stat=err)
          RETURN
       END IF

    END IF

    obs_arr => getObservations(obss)
    IF (error) THEN
       CALL errorMessage("io / writeResiduals", &
            "TRACE BACK (35)", 1)
       RETURN
    END IF
    nobs = SIZE(obs_arr, dim=1)

    ALLOCATE(jds(nobs), codes(nobs), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResiduals", &
            "Could not allocate memory (10).", 1)
       RETURN
    END IF

    DO i=1,nobs
       t = getTime(obs_arr(i))
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "TRACE BACK (40)", 1)
          RETURN
       END IF
       CALL getCalendarDate(t, "TT", year, month, day)
       IF (year >= 1972) THEN
          jds(i) = getJD(t, "UTC")
       ELSE
          jds(i) = getJD(t, "UT1")
       END IF
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "TRACE BACK (45)", 1)
          RETURN
       END IF
       codes(i) = getCode(obs_arr(i))
       IF (error) THEN
          CALL errorMessage("io / writeResiduals", &
               "TRACE BACK (50)", 1)
          RETURN
       END IF
    END DO

    id = getID(obss)
    IF (error) THEN
       CALL errorMessage("io / writeResiduals", &
            "TRACE BACK (52)", 1)
       RETURN
    END IF
    CALL toString(nobs, str, error)
    IF (error) THEN
       CALL errorMessage("io / writeResiduals", &
            "TRACE BACK (55)", 1)
       RETURN
    END IF

    WRITE(lu, "('#',4X,'ID',4X,"//TRIM(str)//"(7X,F16.8,5X))", iostat=err) jds
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResiduals", &
            "Write error (5).", 1)
       RETURN
    END IF
    WRITE(lu, "('#',10X,"//TRIM(str)//"(11X,A4,13X))", iostat=err) codes
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResiduals", &
            "Write error (10).", 1)
       RETURN
    END IF

    DO i=1,SIZE(residuals_,dim=1)
       WRITE(lu, "(A10,1X,"//TRIM(str)//"(2(E13.6,1X)))", iostat=err) &
            TRIM(id), residuals_(i,:,2:3)/rad_asec
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResiduals", &
               "Write error (15).", 1)
          RETURN
       END IF
    END DO

    IF (PRESENT(residuals)) THEN
       ALLOCATE(residuals(SIZE(residuals_,dim=1),nobs,6), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeResiduals", &
               "Could not allocate memory (15).", 1)
          RETURN
       END IF
       residuals = residuals_
    END IF

    DEALLOCATE(obs_arr, residuals_, jds, codes, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("io / writeResiduals", &
            "Could not deallocate memory (15).", 1)
       DEALLOCATE(obs_arr, stat=err)
       DEALLOCATE(residuals_, stat=err)
       DEALLOCATE(jds, stat=err)
       DEALLOCATE(codes, stat=err)
       RETURN
    END IF


  END SUBROUTINE writeResiduals





  SUBROUTINE writeSORResults(storb, obss, lu, impact)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout)            :: storb
    TYPE (Observations), INTENT(in)                  :: obss
    INTEGER, INTENT(in)                              :: lu
    LOGICAL, INTENT(in), OPTIONAL                    :: impact

    TYPE (Orbit), DIMENSION(:), POINTER :: orb_arr_cmp
    TYPE (Time) :: t
    CHARACTER(len=ELEMENT_TYPE_LEN) :: &
         element_type_prm
    CHARACTER(len=DYN_MODEL_LEN) :: &    
         dyn_model_prm
    CHARACTER(len=INTEGRATOR_LEN) :: &
         integrator
    CHARACTER(len=1024) :: frmt
    CHARACTER(len=145) :: str1, str2, str3
    CHARACTER(len=64) :: sor_2point_method
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         res_arr_cmp, sor_deviates_prm
    REAL(bp), DIMENSION(:,:), POINTER :: &
         sor_rho_arr_cmp, res_accept_prm, obs_stdev_arr, jac_arr_cmp, &
         rms_arr_cmp
    REAL(bp), DIMENSION(:), POINTER :: &
         pdf_arr_cmp, reg_apr_arr_cmp, rchi2_arr_cmp
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         res_arr, elements_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         ellipse_fac
    REAL(bp), DIMENSION(2,6,2) :: &
         conf_limits 
    REAL(bp), DIMENSION(2,2) :: &
         sor_rho_prm, &
         sor_rho_cmp 
    REAL(bp), DIMENSION(6) :: elements
    REAL(bp) :: &
         sor_generat_multiplier_prm, accept_multiplier_prm, obsarc, &
         prob_mass, pdf_ml_prm, apriori_a_min_prm, apriori_a_max_prm, &
         apriori_periapsis_min_prm, apriori_periapsis_max_prm, &
         apriori_apoapsis_min_prm, apriori_apoapsis_max_prm, &
         apriori_rho_min_prm, apriori_rho_max_prm
    INTEGER, DIMENSION(:,:), POINTER :: &
         sor_pair_arr_prm
    INTEGER, DIMENSION(2) :: &
         npoints
    INTEGER :: & 
         sor_norb_cmp, sor_ntrial_cmp, &
         sor_norb_prm, sor_ntrial_prm, sor_rho_histo_cmp
    INTEGER :: &
         nobs, i, k, err, indx_ml, nra, ndec
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks
    LOGICAL, DIMENSION(:), ALLOCATABLE :: mask_arr
    LOGICAL :: impact_, sor_random_obs_prm, regularization_prm, uniform_pdf_prm

    IF (.NOT. exist(storb)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "StochasticOrbit object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "Observations object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (PRESENT(impact)) THEN
       impact_ = impact
    ELSE
       impact_ = .FALSE.
    END IF

    IF (impact_) THEN
       WRITE(lu,"(A,3X,A,A)", advance="no") "#","IMPACT PROBABILITY (MC) FOR ", getID(obss)
    ELSE
       WRITE(lu,"(A,3X,A,A)", advance="no") "#","STATISTICAL ORBITAL RANGING FOR ", getID(obss)
    END IF

    nobs = getNumberOfObservations(obss)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "TRACE BACK 20", 1)
       RETURN
    END IF
    obsarc = getObservationalTimespan(obss)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "TRACE BACK 20", 1)
       RETURN
    END IF
    obs_stdev_arr => getStandardDeviations(obss)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "TRACE BACK 20", 1)
       RETURN
    END IF
    obs_stdev_arr = obs_stdev_arr/rad_asec
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "TRACE BACK 20", 1)
       RETURN
    END IF
    obs_masks => getObservationMasks(storb)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "TRACE BACK 10", 1)
       RETURN
    END IF

    nra=COUNT(obs_masks(:,2))
    ndec=COUNT(obs_masks(:,3))

    frmt = "(/'#',3X,'Number of initial observations = ',I14/" // &
         "'#',3X,'Number of R.A. included        = ',I14/" // &
         "'#',3X,'Number of Dec. included        = ',I14/" // &
         "'#',3X,'Observational time arc         = ',F14.4,3X,'days'/" // &
         "'#')"
    WRITE(lu,TRIM(frmt)) nobs, nra, ndec, obsarc

    CALL getParameters(storb, &
         dyn_model=dyn_model_prm, integrator=integrator, &
         element_type = element_type_prm, &
         uniform_pdf = uniform_pdf_prm, regularized_pdf = regularization_prm, &
         accept_multiplier = accept_multiplier_prm, &
         res_accept = res_accept_prm, &
         pdf_ml_prm = pdf_ml_prm, &
         prob_mass = prob_mass, &
         apriori_a_min = apriori_a_min_prm, &
         apriori_a_max = apriori_a_max_prm, &
         apriori_periapsis_min = apriori_periapsis_min_prm, &
         apriori_periapsis_max = apriori_periapsis_max_prm, &
         apriori_apoapsis_min = apriori_apoapsis_min_prm, &
         apriori_apoapsis_max = apriori_apoapsis_max_prm, &
         apriori_rho_min = apriori_rho_min_prm, &
         apriori_rho_max = apriori_rho_max_prm, &
         sor_2point_method = sor_2point_method, &
         sor_norb = sor_norb_prm, sor_ntrial = sor_ntrial_prm, &
         sor_rho1_l=sor_rho_prm(1,1), sor_rho1_u=sor_rho_prm(1,2), &
         sor_rho2_l=sor_rho_prm(2,1), sor_rho2_u=sor_rho_prm(2,2), &
         sor_random_obs_selection=sor_random_obs_prm, & 
         sor_generat_multiplier = sor_generat_multiplier_prm, &
         sor_deviates = sor_deviates_prm)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeSORResults", &
            "TRACE BACK 15", 1)
       RETURN
    END IF

    IF (.NOT.impact_) THEN

       orb_arr_cmp => getSampleOrbits(storb)
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "TRACE BACK 20", 1)
          RETURN
       END IF
       t = getTime(orb_arr_cmp(1))

       pdf_arr_cmp => getPDFValues(storb)
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "TRACE BACK 20", 1)
          RETURN
       END IF
       rchi2_arr_cmp => getReducedChi2Distribution(storb)
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "TRACE BACK 20", 1)
          RETURN
       END IF
       res_arr_cmp => getResidualDistribution(storb)
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "TRACE BACK 20", 1)
          RETURN
       END IF
       rms_arr_cmp => getRMSDistribution(storb)
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "TRACE BACK 20", 1)
          RETURN
       END IF
       CALL getResults(storb, reg_apr_arr=reg_apr_arr_cmp, &
            jac_arr=jac_arr_cmp, &
            sor_norb_cmp=sor_norb_cmp, sor_ntrial_cmp=sor_ntrial_cmp, &
            sor_rho_cmp=sor_rho_cmp, sor_rho_histo_cmp=sor_rho_histo_cmp)
       IF (error) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "TRACE BACK 25", 1)
          RETURN
       END IF
       indx_ml = MAXLOC(pdf_arr_cmp,dim=1)
       ALLOCATE(elements_arr(SIZE(orb_arr_cmp),6))
       str1 = "KEP elements     = "
       str2 = "Confid. limit lo = "
       str3 = "Confid. limit hi = "
       DO i=1,SIZE(orb_arr_cmp)
          elements_arr(i,1:6) =  getElements(orb_arr_cmp(i), "keplerian")
          IF (error) THEN
             EXIT
          END IF
       END DO
       IF (error) THEN
          error = .FALSE.
          str1 = "CAR elements     = "
          str2 = "Confid. limit lo = "
          str3 = "Confid. limit hi = "
          DO i=1,SIZE(orb_arr_cmp)
             elements_arr(i,1:6) =  getElements(orb_arr_cmp(i), "cartesian", "ecliptic")
             IF (error) THEN
                CALL errorMessage("io / writeSORResults", &
                     "TRACE BACK 26", 1)
                RETURN
             END IF
          END DO
       ELSE
          elements_arr(:,3:6) = elements_arr(:,3:6)/rad_deg
       END IF
       ! Compute 1-sigma-equivalent and 3-sigma-equivalent bounds for elements 
       DO i=1,6
          CALL confidence_limits(elements_arr(:,i), pdf_arr_cmp, &
               probability_mass=0.6827_bp, peak=elements(i), bounds=conf_limits(1,i,:), &
               error=errstr)
          CALL confidence_limits(elements_arr(:,i), pdf_arr_cmp, &
               probability_mass=0.9973_bp, peak=elements(i), bounds=conf_limits(2,i,:), &
               error=errstr)
          IF (LEN_TRIM(errstr) /= 0) THEN
             error = .TRUE.
             CALL errorMessage("io / writeSORResults", &
                  "Could not compute confidence limits. " // TRIM(errstr), 1)
             errstr = ""
             RETURN
          END IF
       END DO
       IF (ANY(elements-elements_arr(indx_ml,:) /= 0.0_bp)) THEN
          WRITE(stdout,*) elements
          WRITE(stdout,*) elements_arr(indx_ml,:)
       END IF
       DEALLOCATE(elements_arr, stat=err)

       frmt = "('#',3X,'ORBITAL-ELEMENT PDF' /" // &
            "'#',3X,' Epoch            = ',A,' = ',F13.5,' TDT'/" // &
            "'#',3X,' Maximum likelihood (ML) orbit' /" // &
            "'#',3X,'  ',A19,6(F15.10,1X)/" // &
            "'#',3X,' 68.27% credible intervals' /" // &
            "'#',3X,'  ',A19,6(F15.10,1X)/" // &
            "'#',3X,'  ',A19,6(F15.10,1X)/" // &
            "'#',3X,' 99.73% credible intervals' /" // &
            "'#',3X,'  ',A19,6(F15.10,1X)/" // &
            "'#',3X,'  ',A19,6(F15.10,1X)/" // &
            "'#',3X,'  ML value        =',E16.6/" // &
            "'#',3X,'  ML reduced chi2 =',E16.6/" // &
            "'#',3X,'  ML rms          =',E16.6,3X,'arcsec'/" // &            
            "'#',3X,' Apriori pdf, min =',E16.6/" // &
            "'#',3X,'              max =',E16.6/" // &
            "'#',3X,' Jacobian,    min =',E16.6/" // &
            "'#',3X,'              max =',E16.6/" // &
            "'#',3X,' Rms,         min =',E16.6,3X,'arcsec'/" // &
            "'#',3X,'              max =',E16.6,3X,'arcsec'/" // &
            "'#')"
       WRITE(lu,TRIM(frmt)) getCalendarDateString(t,"tdt"), getJD(t,"tdt"), &
            str1(1:19), elements, str2(1:19), conf_limits(1,:,1), str3(1:19), &
            conf_limits(1,:,2), str2(1:19), conf_limits(2,:,1), str3(1:19), &
            conf_limits(2,:,2), pdf_arr_cmp(indx_ml),rchi2_arr_cmp(indx_ml)+nra+ndec,&
            SQRT(0.5*(rms_arr_cmp(indx_ml,2)**2+rms_arr_cmp(indx_ml,3)**2))/rad_asec,&
            MINVAL(reg_apr_arr_cmp),MAXVAL(reg_apr_arr_cmp),&
            MINVAL(jac_arr_cmp(:,1)),MAXVAL(jac_arr_cmp(:,1)),&
            MINVAL(SQRT(0.5*(rms_arr_cmp(:,2)**2+rms_arr_cmp(:,3)**2)))/rad_asec,&
            MAXVAL(SQRT(0.5*(rms_arr_cmp(:,2)**2+rms_arr_cmp(:,3)**2)))/rad_asec

       frmt = "('#',3X,'COMPUTATIONAL PARAMETERS'/" // &
            "'#',3X,'Element set      = ',A/" // &
            "'#',3X,'Two-point method = ',A/" // &
            "'#',3X,'Dynamical model  = ',A/" // &
            "'#',3X,'Regularization   = ',L2/" // &
            "'#',3X,'Uniform PDF      = ',L2/" // &
            "'#')"
       WRITE(lu,TRIM(frmt)) element_type_prm, &
            TRIM(sor_2point_method), &
            TRIM(dyn_model_prm), &
            regularization_prm, &
            uniform_pdf_prm


       frmt = "('#',3X,'BAYESIAN A PRIORI INFORMATION'/" // &
            "'#',3X,'Semimajor axis (AU), min      = ',F15.9/" // &
            "'#',3X,'                     max      = ',F15.9/" // &
            "'#',3X,'Periapsis distance (AU),  min = ',F15.9/" // &
            "'#',3X,'                          max = ',F15.9/" // &
            "'#',3X,'Apoapsis distance (AU),   min = ',F15.9/" // &
            "'#',3X,'                          max = ',F15.9/" // &
            "'#',3X,'Topocentric range (AU),   min = ',F15.9/" // &
            "'#',3X,'                          max = ',F15.9/" // &
            "'#')"
       WRITE(lu,TRIM(frmt)) apriori_a_min_prm, &
            apriori_a_max_prm, &
            apriori_periapsis_min_prm, &
            apriori_periapsis_max_prm, &
            apriori_apoapsis_min_prm, &
            apriori_apoapsis_max_prm, &
            apriori_rho_min_prm, &
            apriori_rho_max_prm
    END IF

    frmt = "('#',3X,'Final   number of sample orbits  = ',I14/" // &
         "'#',3X,'Initial number of sample orbits  = ',I14/" // &
         "'#',3X,'Final   number of trials         = ',I14/" // &
         "'#',3X,'Initial number of trials         = ',I14/" // &
         "'#')"
    WRITE(lu,TRIM(frmt)) sor_norb_cmp, &
         sor_norb_prm, &
         sor_ntrial_cmp,&
         sor_ntrial_prm

    frmt = "('#',3X,'RESIDUALS AND PDF'/" // &
         "'#',3X,'R.A.*cos Dec. std        (min)   = ',E14.4,3X,'arcsec'/" // &
         "'#',3X,'Dec. std                 (min)   = ',E14.4,3X,'arcsec'/" // &
         "'#',3X,'Acceptance, residuals'/" // &
         "'#',3X,'  sigma multiplier               = ',E14.4/" // &
         "'#',3X,'  window for 1st R.A.            = ',E14.4,3X,'arcsec'/" // &
         "'#',3X,'  window for 1st Dec.            = ',E14.4,3X,'arcsec'/" // &
         "'#',3X,'Acceptance, PDF                  ',3X,/" // &
         "'#',3X,'  reference ML value             = ',E14.4/" // &
         "'#',3X,'  chi-square difference          = ',E14.4/" // &
         "'#',3X,'  corresponding relative bound   = ',E14.4/" // &
         "'#',3X,'  Delta dchi-square              = ',E14.4/" // &
         "'#')"
    WRITE(lu,TRIM(frmt)) & 
         MINVAL(obs_stdev_arr(:,2:3),dim=1), &
         accept_multiplier_prm, &
         res_accept_prm(1,2:3)/rad_asec,&
         pdf_ml_prm, &
         prob_mass, &
         EXP(-0.5_bp*prob_mass), &
         ABS(2.0_bp*LOG(pdf_arr_cmp(indx_ml)/pdf_ml_prm))

    IF (impact_) THEN

       mask_arr(:) = .FALSE.
       WHERE(sor_rho_arr_cmp(:,1) <= planetary_radii(3))
          mask_arr = .TRUE.
       END WHERE
       frmt = "('#',3X,'Impact probability  = ',E14.4/)"
       WRITE(lu,TRIM(frmt)) SUM(pdf_arr_cmp,mask_arr)/SUM(pdf_arr_cmp)

    ELSE

       sor_pair_arr_prm => getObservationPairs(storb)
       IF (.NOT. sor_random_obs_prm) THEN
          CALL toString(sor_pair_arr_prm(1,1), str1, error)
          IF (error) THEN
             error = .TRUE.
             CALL errorMessage("io / writeSORResults", &
                  "TRACE BACK 20", 1)
             RETURN
          END IF
          CALL toString(sor_pair_arr_prm(1,2), str2, error)
          IF (error) THEN
             error = .TRUE.
             CALL errorMessage("io / writeSORResults", &
                  "TRACE BACK 20", 1)
             RETURN
          END IF
       ELSE
          str1 = "<random>"
          str2 = "<random>"
       END IF

       ALLOCATE(mask_arr(sor_norb_cmp), &
            res_arr(2,sor_norb_cmp), &
            ellipse_fac(sor_norb_cmp), &
            stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("io / writeSORResults", &
               "Could not allocate memory.", 1)
          RETURN       
       END IF

       DO k=1,2
          res_arr(1,:) = res_arr_cmp(:,sor_pair_arr_prm(1,k),2)
          res_arr(2,:) = res_arr_cmp(:,sor_pair_arr_prm(1,k),3)       
          ellipse_fac = &
               (res_arr(1,:)/res_accept_prm(sor_pair_arr_prm(1,k),2))**2 + &
               (res_arr(2,:)/res_accept_prm(sor_pair_arr_prm(1,k),3))**2
          mask_arr(:) = .FALSE.
          WHERE (ellipse_fac > 1.0_bp)
             mask_arr = .TRUE.
          endwhere
          npoints(k) = COUNT(mask_arr)
       END DO

       sor_deviates_prm = sor_deviates_prm/rad_asec

       frmt = "('#',3X,'GENERATION WINDOWS FOR RHO, R.A., AND DEC.' /" // &
            "'#',3X,'  Id. number of 1st observation  = ',A14/" // &
            "'#',3X,'  Id. number of 2nd observation  = ',A14/" // &
            "'#',3X,'  Bound for rho1;      lower     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'                       upper     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'  Bound for rho2-rho1; lower     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'                       upper     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'  sigma multiplier               = ',E14.4/" // &
            "'#',3X,'  window shift for 1st R.A.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'                   1st Dec.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'  window width for 1st R.A.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'                   1st Dec.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'  window shift for 2nd R.A.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'                   2nd Dec.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'  window width for 2nd R.A.      = ',E14.4,3X,'arcsec'/" // &
            "'#',3X,'                   2nd Dec.      = ',E14.4,3X,'arcsec'/" // &
            "'#')"
       WRITE(lu,TRIM(frmt)) ADJUSTR(str1(1:14)), &
            ADJUSTR(str2(1:14)), &
            sor_rho_prm(1,1:2), &
            sor_rho_prm(2,1:2), &
            sor_generat_multiplier_prm, &
            sor_deviates_prm(sor_pair_arr_prm(1,1),2:3,1),&
            sor_deviates_prm(sor_pair_arr_prm(1,1),2:3,2),&
            sor_deviates_prm(sor_pair_arr_prm(1,2),2:3,1),&
            sor_deviates_prm(sor_pair_arr_prm(1,2),2:3,2)

       frmt = "('#',3X,'COMPUTED VALUES FOR RHO, R.A., AND DEC.'/" // &
            "'#',3X,'Computed, rho' /" // &
            "'#',3X,'  Bound for rho1,      lower     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'                       upper     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'  Bound for rho2-rho1, lower     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'                       upper     = ',E14.4,3X,'AU'/" // &
            "'#',3X,'  Histogram flag                 = ',I14/" // &
            "'#',3X,'Computed, R.A. and Dec. residuals' /" // &
            "'#',3X,'  Fraction outside ref. ellipse, 1st obs =',E12.4/" // &
            "'#',3X,'  Fraction outside ref. ellipse, 2nd obs =',E12.4/" // &
            "'#')"
       WRITE(lu,TRIM(frmt)) sor_rho_cmp(1,1:2), &
            sor_rho_cmp(2,1:2), &
            sor_rho_histo_cmp, & 
            npoints(1)/REAL(sor_norb_cmp), &
            npoints(2)/REAL(sor_norb_cmp)

    END IF

    DO i=1,SIZE(orb_arr_cmp)
       CALL NULLIFY(orb_arr_cmp(i))
    END DO
    DEALLOCATE(orb_arr_cmp, stat=err)
    DEALLOCATE(res_arr_cmp, stat=err)
    DEALLOCATE(rms_arr_cmp, stat=err)
    DEALLOCATE(sor_deviates_prm, stat=err)
    DEALLOCATE(sor_rho_arr_cmp, stat=err)
    DEALLOCATE(res_accept_prm, stat=err)
    DEALLOCATE(jac_arr_cmp, stat=err)
    DEALLOCATE(pdf_arr_cmp, stat=err)
    DEALLOCATE(reg_apr_arr_cmp, stat=err)
    DEALLOCATE(rchi2_arr_cmp, stat=err)
    DEALLOCATE(res_arr, stat=err)
    DEALLOCATE(ellipse_fac, stat=err)
    DEALLOCATE(sor_pair_arr_prm, stat=err)
    DEALLOCATE(obs_stdev_arr, stat=err)
    DEALLOCATE(obs_masks, stat=err)
    DEALLOCATE(mask_arr, stat=err)

  END SUBROUTINE writeSORResults





  SUBROUTINE writeVOVResults(storb, obss, lu)

    IMPLICIT NONE
    TYPE (StochasticOrbit), INTENT(inout) :: storb
    TYPE (Observations), INTENT(in)       :: obss
    INTEGER, INTENT(in)                   :: lu
    TYPE (Orbit), DIMENSION(:), POINTER   :: orb_arr_cmp
    TYPE (Time) :: t
    REAL(bp), DIMENSION(:,:,:), POINTER :: &
         res_arr_cmp
    REAL(bp), DIMENSION(:,:), POINTER :: &
         res_accept_prm, stdevs, jac_arr_cmp, rms_arr_cmp
    REAL(bp), DIMENSION(:), POINTER :: &
         pdf_arr_cmp, reg_apr_arr_cmp, rchi2_arr_cmp
    REAL(bp), DIMENSION(:,:), ALLOCATABLE :: &
         res_arr
    REAL(bp), DIMENSION(:), ALLOCATABLE :: &
         ellipse_fac
    REAL(bp), DIMENSION(6,2) :: vov_scaling_prm
    REAL(bp), DIMENSION(6) :: elements, elem_stdevs
    REAL(bp) :: &
         accept_multiplier_prm, obsarc, &
         prob_mass, pdf_ml_prm
    INTEGER :: & 
         vov_norb_cmp, vov_ntrial_cmp, &
         vov_norb_prm, vov_ntrial_prm, vov_nmap_prm 
    INTEGER :: &
         nobs, i, err, indx_ml, nra, ndec
    LOGICAL, DIMENSION(:,:), POINTER :: obs_masks
    LOGICAL, DIMENSION(6,2) :: &
         vov_scaling_ready_cmp
    LOGICAL, DIMENSION(6) :: &
         vov_mapping_mask_prm
    LOGICAL :: regularization_prm, uniform_pdf_prm
    CHARACTER(len=ELEMENT_TYPE_LEN) :: &
         element_type_prm
    CHARACTER(len=DYN_MODEL_LEN) :: &    
         dyn_model_prm
    CHARACTER(len=INTEGRATOR_LEN) :: &
         integrator
    CHARACTER(len=20) :: str1, str2

    IF (.NOT. exist(storb)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVOVResults", &
            "StochasticOrbit object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. exist(obss)) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVOVResults", &
            "Observations object has not yet been initialized.", 1)
       RETURN
    END IF

    WRITE(lu,"(A,3X,A,A)", advance="no") "#", &
         "VOLUME-OF-VARIATION SAMPLING FOR ", getID(obss)

    nobs = getNumberOfObservations(obss)
    obsarc = getObservationalTimespan(obss)
    stdevs => getStandardDeviations(obss)
    stdevs = stdevs/rad_asec
    obs_masks => getObservationMasks(storb)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVOVResults", &
            "TRACE BACK", 1)
       RETURN
    END IF

    nra=COUNT(obs_masks(:,2))
    ndec=COUNT(obs_masks(:,3))

    WRITE(lu,200) nobs, nra, ndec, obsarc
200 FORMAT(/"#",3X,"Number of initial observations = ",I14/ &
         "#",3X,"Number of R.A. included        = ",I14/ &
         "#",3X,"Number of Dec. included        = ",I14/ &
         "#",3X,"Observational time arc         = ",F14.4,3X,"days")

    WRITE(lu,"(A)") "#"

    CALL getParameters(storb, dyn_model=dyn_model_prm, integrator=integrator, &
         element_type = element_type_prm, &
         uniform_pdf = uniform_pdf_prm, regularized_pdf = regularization_prm, &
         accept_multiplier = accept_multiplier_prm, &
         res_accept = res_accept_prm, &
         pdf_ml_prm = pdf_ml_prm, &
         prob_mass = prob_mass, &
         vov_norb = vov_norb_prm, vov_ntrial = vov_ntrial_prm, &
         vov_nmap = vov_nmap_prm, vov_scaling = vov_scaling_prm, &
         vov_mapping_mask=vov_mapping_mask_prm)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVOVResults", &
            "TRACE BACK", 1)
       RETURN
    END IF

    orb_arr_cmp => getSampleOrbits(storb)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVOVResults", &
            "TRACE BACK", 1)
       RETURN
    END IF
    t = getTime(orb_arr_cmp(1))

    pdf_arr_cmp => getPDFValues(storb)
    rchi2_arr_cmp => getReducedChi2Distribution(storb)
    rms_arr_cmp => getRMSDistribution(storb)
    CALL getResults(storb, reg_apr_arr=reg_apr_arr_cmp, &
         jac_arr=jac_arr_cmp, &
         vov_norb_cmp=vov_norb_cmp, &
         vov_ntrial_cmp=vov_ntrial_cmp, &
         vov_scaling_ready_cmp=vov_scaling_ready_cmp)
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVOVResults", &
            "TRACE BACK", 1)
       RETURN
    END IF

    indx_ml = MAXLOC(pdf_arr_cmp,dim=1)

    elements = getElements(orb_arr_cmp(indx_ml), "keplerian")
    IF (error) THEN
       error = .FALSE.
       elements = getElements(orb_arr_cmp(indx_ml), "cartesian")
       str1 = "CAR elements    = "
       str2 = "CAR stdevs      = "
    ELSE
       elements(3:6) = elements(3:6)/rad_deg
       str1 = "KEP elements    = "
       str2 = "KEP stdevs      = "
    END IF

    elem_stdevs = getStandardDeviations(storb, "keplerian")
    IF (error) THEN
       error = .TRUE.
       CALL errorMessage("io / writeVoVResults", &
            "TRACE BACK 30", 1)
       RETURN
    END IF
    elem_stdevs(3:6) = elem_stdevs(3:6)/rad_deg

    WRITE(lu,700) getCalendarDateString(t,"tdt"), getJD(t,"tdt"), &
         str1(1:18), elements, str2(1:18), elem_stdevs, &
         pdf_arr_cmp(indx_ml),rchi2_arr_cmp(indx_ml)+nra+ndec,&
         SQRT(0.5*(rms_arr_cmp(indx_ml,2)**2+rms_arr_cmp(indx_ml,3)**2))/rad_asec,&
         MINVAL(reg_apr_arr_cmp),MAXVAL(reg_apr_arr_cmp),&
         MINVAL(jac_arr_cmp(:,1)),MAXVAL(jac_arr_cmp(:,1)),&
         MINVAL(SQRT(0.5*(rms_arr_cmp(:,2)**2+rms_arr_cmp(:,3)**2)))/rad_asec,&
         MAXVAL(SQRT(0.5*(rms_arr_cmp(:,2)**2+rms_arr_cmp(:,3)**2)))/rad_asec

700 FORMAT("#",3X,"ORBITAL-ELEMENT PDF" /&
         "#",3X," Epoch            = ",A," = ",F13.5," TDT"/&
         "#",3X," Maximum likelihood (ML) orbit" /&
         "#",3X,"  ",A18,6(F15.10,1X)/&
         "#",3X,"  ",A18,6(F15.10,1X)/&
         "#",3X,"  ML value        =",E16.6/ &
         "#",3X,"  ML reduced chi2 =",E16.6/ &
         "#",3X,"  ML rms          =",E16.6,3X,"arcsec"/ &
         "#",3X," Apriori pdf, min =",E16.6/ &
         "#",3X,"              max =",E16.6/ &
         "#",3X," Jacobian,    min =",E16.6/ &
         "#",3X,"              max =",E16.6/ &
         "#",3X," Rms,         min =",E16.6,3X,"arcsec"/ &
         "#",3X,"              max =",E16.6,3X,"arcsec")
    WRITE(lu,"(A)") "#"

    WRITE(lu,220) element_type_prm, &
         TRIM(dyn_model_prm), &
         regularization_prm, uniform_pdf_prm
220 FORMAT("#",3X,"COMPUTATIONAL PARAMETERS"/&
         "#",3X,"Element set      = ",A/ &
         "#",3X,"Dynamical model  = ",A/ &
         "#",3X,"Regularization   = ",L2/&
         "#",3X,"Uniform PDF      = ",L2)

    WRITE(lu,"(A)") "#"

    WRITE(lu,250) vov_norb_cmp, vov_norb_prm, vov_ntrial_cmp,&
         vov_ntrial_prm
250 FORMAT("#",3X,"Final   number of sample orbits  = ",I14/ &
         "#",3X,"Initial number of sample orbits  = ",I14/ &
         "#",3X,"Final   number of trials         = ",I14/ &
         "#",3X,"Initial number of trials         = ",I14)

    WRITE(lu,"(A)") "#"

    WRITE(lu,300) & 
         MINVAL(stdevs(:,2:3),dim=1),&
         accept_multiplier_prm,res_accept_prm(1,2:3)/rad_asec,&
         pdf_ml_prm,prob_mass,EXP(-0.5_bp*prob_mass)

300 FORMAT(&
         "#",3X,"R.A.*cos Dec. std        (min)   = ",E14.4,3X,"arcsec"/ &
         "#",3X,"Dec. std                 (min)   = ",E14.4,3X,"arcsec"/ &
         "#",3X,"Acceptance, residuals"/  &
         "#",3X,"  sigma multiplier               = ",E14.4/ &
         "#",3X,"  window for 1st R.A.            = ",E14.4,3X,"arcsec"/ &
         "#",3X,"  window for 1st Dec.            = ",E14.4,3X,"arcsec"/ &
         "#",3X,"Acceptance, PDF                  ",3X,/&
         "#",3X,"  reference ML value             = ",E14.4/ &
         "#",3X,"  chi-square difference          = ",E14.4/ &
         "#",3X,"  corresponding relative bound   = ",E14.4)

    WRITE(lu,"(A)") "#"

    WRITE(lu,400) vov_mapping_mask_prm, &
         vov_scaling_prm(:,1), &
         vov_scaling_prm(:,2), &
         vov_scaling_ready_cmp(:,1), &
         vov_scaling_ready_cmp(:,2)
400 FORMAT("#",3X,"Mapping mask         = ",6(3X,L2,1X)/ &
         "#",3X,"Scaling factors, lower = ",6(F5.1,1X)/ &
         "#",3X,"                 upper = ",6(F5.1,1X)/ &
         "#",3X,"Scaling ready,   lower = ",6(L3,3X)/ &
         "#",3X,"                 upper = ",6(L3,3X))

!!$    ALLOCATE(mask_arr(sor_norb_cmp),res_arr(2,sor_norb_cmp),&
!!$         ellipse_fac(sor_norb_cmp),stat=err)
!!$    IF (err /= 0) THEN
!!$       error = .TRUE.
!!$       CALL errorMessage("io / writeVOVResults", &
!!$            "Could not allocate memory.", 1)
!!$       RETURN       
!!$    END IF
!!$
!!$    DO k=1,2
!!$       res_arr(1,:) = res_arr_cmp(:,sor_pair_arr_prm(1,k),2)
!!$       res_arr(2,:) = res_arr_cmp(:,sor_pair_arr_prm(1,k),3)       
!!$       ellipse_fac = &
!!$            (res_arr(1,:)/res_accept_prm(sor_pair_arr_prm(1,k),2))**2 + &
!!$            (res_arr(2,:)/res_accept_prm(sor_pair_arr_prm(1,k),3))**2
!!$       mask_arr(:) = .FALSE.
!!$       WHERE (ellipse_fac > 1.0_bp)
!!$          mask_arr = .TRUE.
!!$       endwhere
!!$       npoints(k) = COUNT(mask_arr)
!!$    END DO
!!$
!!$    sor_deviates_prm=sor_deviates_prm/rad_asec
!!$
!!$    WRITE(lu,400) ADJUSTR(str1(1:14)), ADJUSTR(str2(1:14)), &
!!$         sor_generat_multiplier_prm, &
!!$         sor_deviates_prm(sor_pair_arr_prm(1,1),2:3,1),&
!!$         sor_deviates_prm(sor_pair_arr_prm(1,1),2:3,2),&
!!$         sor_deviates_prm(sor_pair_arr_prm(1,2),2:3,1),&
!!$         sor_deviates_prm(sor_pair_arr_prm(1,2),2:3,2)
!!$
!!$400 FORMAT("#",3X,"Generated, R.A. and Dec. "/ &
!!$         "#",3X,"  Id. number of 1st observation  = ",A14/ &
!!$         "#",3X,"  Id. number of 2nd observation  = ",A14/ &
!!$         "#",3X,"  sigma multiplier               = ",E14.4/ &
!!$         "#",3X,"  window shift for 1st R.A.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"                   1st Dec.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"  window width for 1st R.A.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"                   1st Dec.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"  window shift for 2nd R.A.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"                   2nd Dec.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"  window width for 2nd R.A.      = ",E14.4,3X,"arcsec"/ &
!!$         "#",3X,"                   2nd Dec.      = ",E14.4,3X,"arcsec")
!!$    WRITE(lu,"(A)") "#"
!!$
!!$    WRITE(lu,550) npoints(1)/REAL(sor_norb_cmp),npoints(2)/REAL(sor_norb_cmp)
!!$550 FORMAT("#",3X,"Computed, R.A. and Dec. residuals" /&
!!$         "#",3X,"  Fraction outside ref. ellipse, 1st obs =",E12.4/ &
!!$         "#",3X,"  Fraction outside ref. ellipse, 2nd obs =",E12.4)
    DO i=1,SIZE(orb_arr_cmp)
       CALL NULLIFY(orb_arr_cmp(i))
    END DO
    DEALLOCATE(orb_arr_cmp, stat=err)
    DEALLOCATE(res_arr_cmp, stat=err)
    DEALLOCATE(rms_arr_cmp, stat=err)
    DEALLOCATE(res_accept_prm, stat=err)
    DEALLOCATE(stdevs, stat=err)
    DEALLOCATE(jac_arr_cmp, stat=err)
    DEALLOCATE(pdf_arr_cmp, stat=err)
    DEALLOCATE(reg_apr_arr_cmp, stat=err)
    DEALLOCATE(obs_masks, stat=err)
    DEALLOCATE(rchi2_arr_cmp, stat=err)
    DEALLOCATE(res_arr, stat=err)
    DEALLOCATE(ellipse_fac, stat=err)

  END SUBROUTINE writeVOVResults





END MODULE io

