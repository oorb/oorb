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
! Implementation notes:
!   1. both oorb_ranging_fast and oorb_lsl_fast have a useless loop in them:
!      either clean them up or make them useful!
! 
! 

! Orbital elements that we support.
! ORBITAL_ELEMENTS = (/ "keplerian",  "delaunay", "poincare", "equinoctial" /)





subroutine oorb_init(ephemFileName, verbosity, error_code)
    ! Initialize the OpenOrb module.
    !
    ! @param ephemFileName: full path of the JPL ephem file (usually 
    !        $OORB_DATA/JPL_ephemeris/de405.dat).
    ! @param verbosity: verbosity level for OpenOrb calls [0, 5] (default 0)
    ! @return error_code: int error code. 0 = success, otherwise, failure.
    use oorb, only: init
    ! Input/Output variables.
    character(len=*), intent(IN)    :: ephemFileName
    integer, intent(in)             :: verbosity
    integer, intent(OUT)            :: error_code

    
    call init(trim(ephemFileName), verbosity, error_code)
    return
end subroutine oorb_init



subroutine oorb_ranging_fast(obscodes,                                  &
                             filters,                                   &
                             track_id,                                  &
                             num_coords,                                &
                             coords,                                    &
                             mjds,                                      &
                             mags,                                      &
                             element_type,                              &
                             num_orbits,                                &
                             out_orbits,                                &
                             error_code)
    use oorb
    use Observations_cl
    use StochasticOrbit_cl
    ! Orbital inversion using statistical orbital ranging, that is,
    ! without making any assumptions on the shape of the resulting
    ! orbital-element pdf.
    
    ! All angles coming in are in decimal degrees. All times in MJD TAI.

    ! Input/Output variables.
    ! Observations (i.e. DiASources) to be used as input. It is assumed that 
    ! they all belong to the same track.
    ! We have to use separate arrays since they are of different type.
    integer, intent(in)                                 :: track_id
    integer, intent(in)                                 :: num_coords
    real(8), dimension(num_coords), intent(in)          :: mjds
    real(8), dimension(num_coords, 4), intent(in)       :: coords
    real(8), dimension(num_coords), intent(in)          :: mags
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=1), dimension(num_coords), intent(in) :: filters
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=4), dimension(num_coords), intent(in) :: obscodes
    ! Element type used in the orbit computation (default keplerian)
    CHARACTER(len=11), intent(in)                       :: element_type
    ! How may orbits do we want out?
    integer, intent(in)                                 :: num_orbits
    ! Output flattened orbits.
    ! out_orbits has the form (id, (elements(1:6), epoch, Un-normalized p.d.f.,
    ! Reduced chi2, Regularized apr, Jacobian det(1:3), element_type_index), )
    real(8),dimension(num_orbits,15),intent(out)        :: out_orbits
    ! Output error code
    integer, intent(out)                                :: error_code
    
    
    ! Internal variables.
    character(len=6)                                    :: dyn_model
    real(8)                                             :: integration_step
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    real(8)                                             :: apriori_a_max
    real(8)                                             :: apriori_a_min
    real(8)                                             :: apriori_rho_max
    real(8)                                             :: apriori_rho_min
    logical                                             :: outlier_rejection
    real(8)                                             :: outlier_multiplier
    integer                                             :: sor_type_prm
    character(len=80)                                   :: sor_2point_method
    character(len=80)                                   :: sor_2point_method_sw
    integer                                             :: sor_norb
    integer                                             :: sor_norb_sw
    integer                                             :: sor_ntrial
    integer                                             :: sor_ntrial_sw
    integer                                             :: sor_niter
    REAL(8), DIMENSION(4)                               :: sor_rho_init
    real(8)                                             :: sor_genwin_multiplier
    REAL(8), DIMENSION(4)                               :: sor_genwin_offset
    real(8)                                             :: accwin_multiplier
    logical                                             :: gaussian_rho
    real(8)                                             :: pdf_ml_init
    logical                                             :: uniform
    logical                                             :: regularized
    logical                                             :: random_obs
    TYPE (Observations), DIMENSION(:), POINTER          :: obss
    TYPE (StochasticOrbit)                              :: storb
    integer                                             :: j
    integer                                             :: element_type_index
    
!    REAL(8), DIMENSION(:,:,:), POINTER                  :: cov_matrices
!    integer                                             :: nobs = 0
!    integer                                             :: i = 0
!    TYPE (SphericalCoordinates), DIMENSION(:), POINTER  :: scoords
    
!    do j=1,num_coords
!        write(*, *) "obscodes: ", obscodes(j)
!    end do
    
    
    
!    write(*, *) "obscodes(1) ", obscodes(1)
!    write(*, *) "obscodes(2) ", obscodes(2)
!    write(*, *) "obscodes(3) ", obscodes(3)
!    write(*, *) "filters(1) ", filters(1)
!    write(*, *) "filters(2) ", filters(2)
!    write(*, *) "filters(3) ", filters(3)
!    write(*, *) "track_id ", track_id
!    write(*, *) "num_coords ", num_coords
!    write(*, *) "coords(1) ", coords(1, :)
!    write(*, *) "coords(2) ", coords(2, :)
!    write(*, *) "coords(3) ", coords(3, :)
!    write(*, *) "mjds ", mjds
!    write(*, *) "mags ", mags
!    write(*, *) "element_type '", element_type, "'"
!    write(*, *) "num_orbits ", num_orbits
    
    ! Init (these values come from best practices and default OpenOrb values).
    error_code = 0
    dyn_model = "2-body"
    integration_step = 5.0_8
    apriori_a_max = 1000.0_8
    apriori_a_min = 4.65423999999999959E-003_8
    apriori_rho_max = -1.0_8
    apriori_rho_min = -1.0_8
    outlier_rejection = .false.
    outlier_multiplier = 4.0_8
    sor_type_prm = 2
    sor_2point_method = "continued fraction"
    sor_2point_method_sw = "continued fraction"
    sor_norb = num_orbits
    sor_norb_sw = sor_norb / 10
    sor_ntrial = sor_norb * 2000
    sor_ntrial_sw = sor_norb_sw * 400
    sor_niter = 3
    sor_rho_init = huge(sor_rho_init)
    sor_rho_init(1) = 0.0_8
    sor_rho_init(2) = 100.0_8
    sor_genwin_multiplier = 3.0_8
    sor_genwin_offset = 0.0_8
    accwin_multiplier = 3.0_8
    gaussian_rho = .false.
    pdf_ml_init = -1.0_8
    uniform = .false.
    regularized = .true.
    random_obs = .false.
    obss => null()
    call nullify(storb)
    j = 0
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Figure out the element_type_index
    element_type_index = -1
    do j=1,size(ORBITAL_ELEMENTS)
        if(trim(element_type) == trim(ORBITAL_ELEMENTS(j))) then
            element_type_index = j
            exit
!        else
!            write(*, *) "'" // trim(element_type) // "' =/ '" // trim(ORBITAL_ELEMENTS(j)) // "'"
!            write(*, *) " ", len(trim(element_type)), " =/ ", len(trim(ORBITAL_ELEMENTS(j)))
        end if
    end do
    if(element_type_index .le. 0) then
        ! Error: unsupported element type!
        error_code = 58
        return
    end if
    
    ! Convert the input observation array into a list of Observations instances.
    ! Input coordinates (and their uncertainties) are in decimal degrees.
    ! We need them in radians.
    call obssFromCoords(track_id, num_coords, coords*rad_deg, mjds, mags,   &
                        filters, obscodes, obss, error_code)
    if(error_code /= 0) then
        return
    end if
    
!    do j=1,size(obss)
!        cov_matrices => getCovarianceMatrices(obss(j))
!        nobs = getNrOfObservations(obss(j))
!        DO i=1,nobs
!            write(*, *) "COV: ", cov_matrices(i,:,:)
!        END DO
!    end do
    
    ! Now that we have a set of Observations, pass them to the ranging routine.
    ! obss has already been separated into tracks, pass one track at the time.
    ! We clearly expect one track only, so this loop is here only so that we can
    ! support multiple track processing later on. But it could just as well be
    ! eliminated...
    do j=1,size(obss)
        call nullify(storb)
        
!        scoords => getObservationSCoords(obss(j))
!        do i=1,size(scoords)
!            write(*, *) "RA:  ", getLongitude(scoords(i))
!            write(*, *) "Dec: ", getLatitude(scoords(i))
!        end do
        
        
        call ranging(obs_in=obss(j), element_type=trim(element_type), &
                     dyn_model=trim(dyn_model), &
                     integration_step=integration_step, &
                     perturbers=perturbers, apriori_a_max=apriori_a_max, &
                     apriori_a_min=apriori_a_min, &
                     apriori_rho_max=apriori_rho_max, &
                     apriori_rho_min=apriori_rho_min, &
                     outlier_rejection=outlier_rejection, &
                     outlier_multiplier=outlier_multiplier, &
                     sor_type_prm=sor_type_prm, &
                     sor_2point_method=trim(sor_2point_method), &
                     sor_2point_method_sw=trim(sor_2point_method_sw), &
                     sor_norb=sor_norb, sor_norb_sw=sor_norb_sw, &
                     sor_ntrial=sor_ntrial, sor_ntrial_sw=sor_ntrial_sw, &
                     sor_niter=sor_niter, sor_rho_init=sor_rho_init, &
                     sor_genwin_multiplier=sor_genwin_multiplier, &
                     sor_genwin_offset=sor_genwin_offset, &
                     accwin_multiplier=accwin_multiplier, &
                     gaussian_rho=gaussian_rho, &
                     pdf_ml_init=pdf_ml_init, uniform=uniform, &
                     regularized=regularized, random_obs=random_obs, &
                     outOrbit=storb, errorCode=error_code)
        if(error_code /= 0) then
            return
        end if
        
        ! Dump the storb orbit.
        ! call dumpRangingOrbits(storb, errorCode=error_code)
        
        ! Now flatten the ranging orbits.
        call exportRangingOrbits(storb, element_type, track_id, out_orbits, &
                                 error_code)
        call nullify(storb)
        if(error_code /= 0) then
            return
        end if
    end do
    
    deallocate(perturbers)
    return
end subroutine oorb_ranging_fast



subroutine oorb_ranging(obscodes,                                  &
                        filters,                                   &
                        track_id,                                  &
                        num_coords,                                &
                        coords,                                    &
                        mjds,                                      &
                        mags,                                      &
                        element_type,                              &
                        num_orbits,                                &
                        out_orbits,                                &
                        error_code)
    use oorb
    use Observations_cl
    use StochasticOrbit_cl
    ! Orbital inversion using statistical orbital ranging, that is,
    ! without making any assumptions on the shape of the resulting
    ! orbital-element pdf.
    
    ! All angles coming in are in decimal degrees. All times in MJD TAI.

    ! Input/Output variables.
    ! Observations (i.e. DiASources) to be used as input. It is assumed that 
    ! they all belong to the same track.
    ! We have to use separate arrays since they are of different type.
    integer, intent(in)                                 :: track_id
    integer, intent(in)                                 :: num_coords
    real(8), dimension(num_coords), intent(in)          :: mjds
    real(8), dimension(num_coords, 4), intent(in)       :: coords
    real(8), dimension(num_coords), intent(in)          :: mags
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=1), dimension(num_coords), intent(in) :: filters
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=4), dimension(num_coords), intent(in) :: obscodes
    ! Element type used in the orbit computation (default keplerian)
    CHARACTER(len=11), intent(in)                       :: element_type
    ! How may orbits do we want out?
    integer, intent(in)                                 :: num_orbits
    ! Output flattened orbits.
    ! out_orbits has the form (id, (elements(1:6), epoch, Un-normalized p.d.f.,
    ! Reduced chi2, Regularized apr, Jacobian det(1:3), element_type_index), )
    real(8),dimension(num_orbits,15),intent(out)        :: out_orbits
    ! Output error code
    integer, intent(out)                                :: error_code
    
    
    ! Internal variables.
    character(len=6)                                    :: dyn_model
    real(8)                                             :: integration_step
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    real(8)                                             :: apriori_a_max
    real(8)                                             :: apriori_a_min
    real(8)                                             :: apriori_rho_max
    real(8)                                             :: apriori_rho_min
    logical                                             :: outlier_rejection
    real(8)                                             :: outlier_multiplier
    integer                                             :: sor_type_prm
    character(len=80)                                   :: sor_2point_method
    character(len=80)                                   :: sor_2point_method_sw
    integer                                             :: sor_norb
    integer                                             :: sor_norb_sw
    integer                                             :: sor_ntrial
    integer                                             :: sor_ntrial_sw
    integer                                             :: sor_niter
    REAL(8), DIMENSION(4)                               :: sor_rho_init
    real(8)                                             :: sor_genwin_multiplier
    REAL(8), DIMENSION(4)                               :: sor_genwin_offset
    real(8)                                             :: accwin_multiplier
    logical                                             :: gaussian_rho
    real(8)                                             :: pdf_ml_init
    logical                                             :: uniform
    logical                                             :: regularized
    logical                                             :: random_obs
    TYPE (Observations), DIMENSION(:), POINTER          :: obss
    TYPE (StochasticOrbit)                              :: storb
    integer                                             :: j
    integer                                             :: element_type_index
    
!    REAL(8), DIMENSION(:,:,:), POINTER                  :: cov_matrices
!    integer                                             :: nobs = 0
!    integer                                             :: i = 0
!    TYPE (SphericalCoordinates), DIMENSION(:), POINTER  :: scoords
    
!    do j=1,num_coords
!        write(*, *) "obscodes: ", obscodes(j)
!    end do
    
    
    
!    write(*, *) "obscodes(1) ", obscodes(1)
!    write(*, *) "obscodes(2) ", obscodes(2)
!    write(*, *) "obscodes(3) ", obscodes(3)
!    write(*, *) "filters(1) ", filters(1)
!    write(*, *) "filters(2) ", filters(2)
!    write(*, *) "filters(3) ", filters(3)
!    write(*, *) "track_id ", track_id
!    write(*, *) "num_coords ", num_coords
!    write(*, *) "coords(1) ", coords(1, :)
!    write(*, *) "coords(2) ", coords(2, :)
!    write(*, *) "coords(3) ", coords(3, :)
!    write(*, *) "mjds ", mjds
!    write(*, *) "mags ", mags
!    write(*, *) "element_type '", element_type, "'"
!    write(*, *) "num_orbits ", num_orbits
    
    ! Init (these values come from best practices and default OpenOrb values).
    error_code = 0
    dyn_model = "n-body"
    integration_step = 5.0_8
    apriori_a_max = 1000.0_8
    apriori_a_min = 4.65423999999999959E-003_8
    apriori_rho_max = -1.0_8
    apriori_rho_min = -1.0_8
    outlier_rejection = .false.
    outlier_multiplier = 4.0_8
    sor_type_prm = 2
    sor_2point_method = "continued fraction"
    sor_2point_method_sw = "continued fraction"
    sor_norb = num_orbits
    sor_norb_sw = sor_norb / 10
    sor_ntrial = sor_norb * 2000
    sor_ntrial_sw = sor_norb_sw * 400
    sor_niter = 3
    sor_rho_init = huge(sor_rho_init)
    sor_rho_init(1) = 0.0_8
    sor_rho_init(2) = 100.0_8
    sor_genwin_multiplier = 3.0_8
    sor_genwin_offset = 0.0_8
    accwin_multiplier = 3.0_8
    gaussian_rho = .false.
    pdf_ml_init = -1.0_8
    uniform = .false.
    regularized = .true.
    random_obs = .false.
    obss => null()
    call nullify(storb)
    j = 0
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Figure out the element_type_index
    element_type_index = -1
    do j=1,size(ORBITAL_ELEMENTS)
        if(trim(element_type) == trim(ORBITAL_ELEMENTS(j))) then
            element_type_index = j
            exit
!        else
!            write(*, *) "'" // trim(element_type) // "' =/ '" // trim(ORBITAL_ELEMENTS(j)) // "'"
!            write(*, *) " ", len(trim(element_type)), " =/ ", len(trim(ORBITAL_ELEMENTS(j)))
        end if
    end do
    if(element_type_index .le. 0) then
        ! Error: unsupported element type!
        error_code = 58
        return
    end if
    
    ! Convert the input observation array into a list of Observations instances.
    ! Input coordinates (and their uncertainties) are in decimal degrees.
    ! We need them in radians.
    call obssFromCoords(track_id, num_coords, coords*rad_deg, mjds, mags,   &
                        filters, obscodes, obss, error_code)
    if(error_code /= 0) then
        return
    end if
    
!    do j=1,size(obss)
!        cov_matrices => getCovarianceMatrices(obss(j))
!        nobs = getNrOfObservations(obss(j))
!        DO i=1,nobs
!            write(*, *) "COV: ", cov_matrices(i,:,:)
!        END DO
!    end do
    
    ! Now that we have a set of Observations, pass them to the ranging routine.
    ! obss has already been separated into tracks, pass one track at the time.
    ! We clearly expect one track only, so this loop is here only so that we can
    ! support multiple track processing later on. But it could just as well be
    ! eliminated...
    do j=1,size(obss)
        call nullify(storb)
        
!        scoords => getObservationSCoords(obss(j))
!        do i=1,size(scoords)
!            write(*, *) "RA:  ", getLongitude(scoords(i))
!            write(*, *) "Dec: ", getLatitude(scoords(i))
!        end do
        
        
        call ranging(obs_in=obss(j), element_type=trim(element_type), &
                     dyn_model=trim(dyn_model), &
                     integration_step=integration_step, &
                     perturbers=perturbers, apriori_a_max=apriori_a_max, &
                     apriori_a_min=apriori_a_min, &
                     apriori_rho_max=apriori_rho_max, &
                     apriori_rho_min=apriori_rho_min, &
                     outlier_rejection=outlier_rejection, &
                     outlier_multiplier=outlier_multiplier, &
                     sor_type_prm=sor_type_prm, &
                     sor_2point_method=trim(sor_2point_method), &
                     sor_2point_method_sw=trim(sor_2point_method_sw), &
                     sor_norb=sor_norb, sor_norb_sw=sor_norb_sw, &
                     sor_ntrial=sor_ntrial, sor_ntrial_sw=sor_ntrial_sw, &
                     sor_niter=sor_niter, sor_rho_init=sor_rho_init, &
                     sor_genwin_multiplier=sor_genwin_multiplier, &
                     sor_genwin_offset=sor_genwin_offset, &
                     accwin_multiplier=accwin_multiplier, &
                     gaussian_rho=gaussian_rho, &
                     pdf_ml_init=pdf_ml_init, uniform=uniform, &
                     regularized=regularized, random_obs=random_obs, &
                     outOrbit=storb, errorCode=error_code)
        if(error_code /= 0) then
            return
        end if
        
        ! Dump the storb orbit.
        ! call dumpRangingOrbits(storb, errorCode=error_code)
        
        ! Now flatten the ranging orbits.
        call exportRangingOrbits(storb, element_type, track_id, out_orbits, &
                                 error_code)
        call nullify(storb)
        if(error_code /= 0) then
            return
        end if
    end do
    
    deallocate(perturbers)
    return
end subroutine oorb_ranging



subroutine oorb_lsl_fast(obscodes,                                  &
                         filters,                                   &
                         track_id,                                  &
                         num_coords,                                &
                         coords,                                    &
                         mjds,                                      &
                         mags,                                      &
                         in_orbits,                                 &
                         num_orbits,                                &
                         out_orbit,                                 &
                         out_covariance,                            &
                         out_sigmas,                                &
                         out_correlation,                           &
                         error_code)
    use oorb
    use StochasticOrbit_cl
    ! Orbital inversion using statistical orbital ranging, that is,
    ! without making any assumptions on the shape of the resulting
    ! orbital-element pdf.
    
    ! All angles coming in are in decimal degrees. All times in MJD TAI.

    ! Input/Output variables.
    ! Observations (i.e. DiASources) to be used as input. It is assumed that 
    ! they all belong to the same track.
    ! We have to use separate arrays since they are of different type.
    integer, intent(in)                                 :: track_id
    integer, intent(in)                                 :: num_coords
    real(8), dimension(num_coords), intent(in)          :: mjds
    real(8), dimension(num_coords, 4), intent(in)       :: coords
    real(8), dimension(num_coords), intent(in)          :: mags
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=1), dimension(num_coords), intent(in) :: filters
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=4), dimension(num_coords), intent(in) :: obscodes
    integer, intent(in)                                 :: num_orbits
    ! Input flattened orbits.
    ! in_orbits has the form ((id, elements(1:6), epoch, Un-normalized p.d.f.,
    ! Reduced chi2, Regularized apr, Jacobian det(1:3), el_type_index), )
    real(8),dimension(num_orbits,15),intent(in)         :: in_orbits
    ! Output flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, el_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(out)                  :: out_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(out)                :: out_covariance
    REAL(8), DIMENSION(6), intent(out)                  :: out_sigmas
    REAL(8), DIMENSION(6,6), intent(out)                :: out_correlation
    ! Output error code
    integer, intent(out)                                :: error_code
    
    ! Internal variables.
    character(len=6)                                    :: dyn_model
    character(len=6)                                    :: dyn_model_init
    real(8)                                             :: integration_step
    real(8)                                             :: integration_step_init
    real(8)                                             :: ls_correction_factor
    LOGICAL, DIMENSION(6)                               :: ls_element_mask
    integer                                             :: ls_niter_major_max
    integer                                             :: ls_niter_major_min
    integer                                             :: ls_niter_minor
    logical                                             :: outlier_rejection
    real(8)                                             :: outlier_multiplier
    real(8)                                             :: accwin_multiplier
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    TYPE (Observations), DIMENSION(:), POINTER          :: obss
    TYPE (StochasticOrbit)                              :: storb
    TYPE (StochasticOrbit)                              :: in_storb
    integer                                             :: j
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    
    
    
    
!    do j=1,num_coords
!        write(*, *) "obscodes(", j, ") ", obscodes(j)
!    end do
!    do j=1,num_coords
!        write(*, *) "filters(", j, ") ", filters(j)
!    end do
!    write(*, *) "track_id ", track_id
!    write(*, *) "num_coords ", num_coords
!    do j=1,num_coords
!        write(*, *) "coords(", j, ") ", coords(j,:)
!    end do
!    write(*, *) "mjds ", mjds
!    write(*, *) "mags ", mags
!    write(*, *) "num_orbits ", num_orbits
    
    ! Init variables (from OpenOrb config file and best practices).
    error_code = 0
    dyn_model = "2-body"
    dyn_model_init = "2-body"
    integration_step = 5.0_8
    integration_step_init = 5.0_8
    ls_correction_factor = 0.2_8
    ls_element_mask = .true.
    ls_niter_major_max = 20
    ls_niter_major_min = 2
    ls_niter_minor = 10
    outlier_rejection = .false.
    outlier_multiplier = 3.0_8
    accwin_multiplier = 3.0_8
    
    ! Quick sanity check.
    if(size(in_orbits) .le. 0) then
        ! Nothing to do.
        return
    end if
    
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Get the element type from the input flattened orbits (the are supposed to 
    ! be all the same).
    element_type_index = in_orbits(1, 15)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    
    ! Convert the input observation array into a list of Observations instances.
    ! Input coordinates (and their uncertainties) are in decimal degrees.
    ! We need them in radians.
    call obssFromCoords(track_id, num_coords, coords*rad_deg, mjds, mags,   &
                        filters, obscodes, obss, error_code)
    if(error_code /= 0) then
        return
    end if
    
!    do j=1,size(obss)
!        cov_matrices => getCovarianceMatrices(obss(j))
!        nobs = getNrOfObservations(obss(j))
!        DO i=1,nobs
!            write(*, *) "COV: ", cov_matrices(i,:,:)
!        END DO
!    end do
    
    ! Convert the input orbits (which come from the ranging routine) into a 
    ! StochasticOrbit instance.
    call nullify(in_storb)
    call rangingOrbitsToStochasticOrbit(in_orbits,                          &
                                        num_orbits,                         &
                                        in_storb,                           &
                                        error_code)
    if(error_code /= 0) then
        return
    end if
    
    
    ! Now that we have a set of Observations, pass them to the lsl routine.
    ! obss has already been separated into tracks, pass one track at the time.
    do j=1,size(obss)
        call nullify(storb)
        
        ! Do the LSL thing!
        call lsl(obs_in=obss(j),                                            &
                 element_type=element_type,                                 &
                 dyn_model=dyn_model,                                       &
                 dyn_model_init=dyn_model_init,                             &
                 integration_step=integration_step,                         &
                 integration_step_init=integration_step_init,               &
                 perturbers=perturbers,                                     &
                 ls_correction_factor=ls_correction_factor,                 &
                 ls_element_mask=ls_element_mask,                           &
                 ls_niter_major_max=ls_niter_major_max,                     &
                 ls_niter_major_min=ls_niter_major_min,                     &
                 ls_niter_minor=ls_niter_minor,                             &
                 outlier_rejection=outlier_rejection,                       &
                 outlier_multiplier=outlier_multiplier,                     &
                 accwin_multiplier=accwin_multiplier,                       &
                 inStochOrbit=in_storb,                                     &
                 outStochOrbit=storb,                                       &
                 errorCode=error_code)
        if(error_code /= 0) then
            return
        end if
        
        ! Now flatten the lsl orbit.
        call exportLslOrbit(storb, element_type, track_id, out_orbit,       &
                            out_covariance, out_sigmas, out_correlation,    &
                            error_code)
        call nullify(storb)
        if(error_code /= 0) then
            return
        end if
    end do
    deallocate(perturbers)
end subroutine oorb_lsl_fast



subroutine oorb_lsl(obscodes,                                  &
                    filters,                                   &
                    track_id,                                  &
                    num_coords,                                &
                    coords,                                    &
                    mjds,                                      &
                    mags,                                      &
                    in_orbits,                                 &
                    num_orbits,                                &
                    out_orbit,                                 &
                    out_covariance,                            &
                    out_sigmas,                                &
                    out_correlation,                           &
                    error_code)
    use oorb
    use StochasticOrbit_cl
    ! Orbital inversion using statistical orbital ranging, that is,
    ! without making any assumptions on the shape of the resulting
    ! orbital-element pdf.
    
    ! All angles coming in are in decimal degrees. All times in MJD TAI.

    ! Input/Output variables.
    ! Observations (i.e. DiASources) to be used as input. It is assumed that 
    ! they all belong to the same track.
    ! We have to use separate arrays since they are of different type.
    integer, intent(in)                                 :: track_id
    integer, intent(in)                                 :: num_coords
    real(8), dimension(num_coords), intent(in)          :: mjds
    real(8), dimension(num_coords, 4), intent(in)       :: coords
    real(8), dimension(num_coords), intent(in)          :: mags
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=1), dimension(num_coords), intent(in) :: filters
    ! FIXME: f2py bug that only passes the first char of element of an array of
    ! strings.
    character(len=4), dimension(num_coords), intent(in) :: obscodes
    integer, intent(in)                                 :: num_orbits
    ! Input flattened orbits.
    ! in_orbits has the form ((id, elements(1:6), epoch, Un-normalized p.d.f.,
    ! Reduced chi2, Regularized apr, Jacobian det(1:3), el_type_index), )
    real(8),dimension(num_orbits,15),intent(in)         :: in_orbits
    ! Output flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, el_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(out)                  :: out_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(out)                :: out_covariance
    REAL(8), DIMENSION(6), intent(out)                  :: out_sigmas
    REAL(8), DIMENSION(6,6), intent(out)                :: out_correlation
    ! Output error code
    integer, intent(out)                                :: error_code
    
    ! Internal variables.
    character(len=6)                                    :: dyn_model
    character(len=6)                                    :: dyn_model_init
    real(8)                                             :: integration_step
    real(8)                                             :: integration_step_init
    real(8)                                             :: ls_correction_factor
    LOGICAL, DIMENSION(6)                               :: ls_element_mask
    integer                                             :: ls_niter_major_max
    integer                                             :: ls_niter_major_min
    integer                                             :: ls_niter_minor
    logical                                             :: outlier_rejection
    real(8)                                             :: outlier_multiplier
    real(8)                                             :: accwin_multiplier
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    TYPE (Observations), DIMENSION(:), POINTER          :: obss
    TYPE (StochasticOrbit)                              :: storb
    TYPE (StochasticOrbit)                              :: in_storb
    integer                                             :: j
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    
    
    
    
!    do j=1,num_coords
!        write(*, *) "obscodes(", j, ") ", obscodes(j)
!    end do
!    do j=1,num_coords
!        write(*, *) "filters(", j, ") ", filters(j)
!    end do
!    write(*, *) "track_id ", track_id
!    write(*, *) "num_coords ", num_coords
!    do j=1,num_coords
!        write(*, *) "coords(", j, ") ", coords(j,:)
!    end do
!    write(*, *) "mjds ", mjds
!    write(*, *) "mags ", mags
!    write(*, *) "num_orbits ", num_orbits
    
    ! Init variables (from OpenOrb config file and best practices).
    error_code = 0
    dyn_model = "n-body"
    dyn_model_init = "2-body"
    integration_step = 5.0_8
    integration_step_init = 5.0_8
    ls_correction_factor = 0.2_8
    ls_element_mask = .true.
    ls_niter_major_max = 20
    ls_niter_major_min = 2
    ls_niter_minor = 10
    outlier_rejection = .false.
    outlier_multiplier = 3.0_8
    accwin_multiplier = 3.0_8
    
    ! Quick sanity check.
    if(size(in_orbits) .le. 0) then
        ! Nothing to do.
        return
    end if
    
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Get the element type from the input flattened orbits (the are supposed to 
    ! be all the same).
    element_type_index = in_orbits(1, 15)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    
    ! Convert the input observation array into a list of Observations instances.
    ! Input coordinates (and their uncertainties) are in decimal degrees.
    ! We need them in radians.
    call obssFromCoords(track_id, num_coords, coords*rad_deg, mjds, mags,   &
                        filters, obscodes, obss, error_code)
    if(error_code /= 0) then
        return
    end if
    
!    do j=1,size(obss)
!        cov_matrices => getCovarianceMatrices(obss(j))
!        nobs = getNrOfObservations(obss(j))
!        DO i=1,nobs
!            write(*, *) "COV: ", cov_matrices(i,:,:)
!        END DO
!    end do
    
    ! Convert the input orbits (which come from the ranging routine) into a 
    ! StochasticOrbit instance.
    call nullify(in_storb)
    call rangingOrbitsToStochasticOrbit(in_orbits,                          &
                                        num_orbits,                         &
                                        in_storb,                           &
                                        error_code)
    if(error_code /= 0) then
        return
    end if
    
    
    ! Now that we have a set of Observations, pass them to the lsl routine.
    ! obss has already been separated into tracks, pass one track at the time.
    do j=1,size(obss)
        call nullify(storb)
        
        ! Do the LSL thing!
        call lsl(obs_in=obss(j),                                            &
                 element_type=element_type,                                 &
                 dyn_model=dyn_model,                                       &
                 dyn_model_init=dyn_model_init,                             &
                 integration_step=integration_step,                         &
                 integration_step_init=integration_step_init,               &
                 perturbers=perturbers,                                     &
                 ls_correction_factor=ls_correction_factor,                 &
                 ls_element_mask=ls_element_mask,                           &
                 ls_niter_major_max=ls_niter_major_max,                     &
                 ls_niter_major_min=ls_niter_major_min,                     &
                 ls_niter_minor=ls_niter_minor,                             &
                 outlier_rejection=outlier_rejection,                       &
                 outlier_multiplier=outlier_multiplier,                     &
                 accwin_multiplier=accwin_multiplier,                       &
                 inStochOrbit=in_storb,                                     &
                 outStochOrbit=storb,                                       &
                 errorCode=error_code)
        if(error_code /= 0) then
            return
        end if
        
        ! Now flatten the lsl orbit.
        call exportLslOrbit(storb, element_type, track_id, out_orbit,       &
                            out_covariance, out_sigmas, out_correlation,    &
                            error_code)
        call nullify(storb)
        if(error_code /= 0) then
            return
        end if
    end do
    deallocate(perturbers)
end subroutine oorb_lsl



subroutine oorb_propagate_orbit_fast(in_orbit, in_covariance, in_mjd,       &
                                     out_orbit, out_covariance, out_sigmas, &
                                     out_correlation, error_code)
    use oorb
    use StochasticOrbit_cl
    
    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(in)                   :: in_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(in)                 :: in_covariance
    ! The epoch (MJD TAI) we want to propagate to.
    real(8), intent(in)                                 :: in_mjd
    ! Output flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(out)                  :: out_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(out)                :: out_covariance
    REAL(8), DIMENSION(6), intent(out)                  :: out_sigmas
    REAL(8), DIMENSION(6,6), intent(out)                :: out_correlation
    ! Output error code
    integer, intent(out)                                :: error_code
    
    ! Internal variables.
    TYPE (StochasticOrbit)                              :: storb
    character(len=6)                                    :: dyn_model
    real(8)                                             :: integration_step
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    integer                                             :: track_id
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    
    
    ! Init
    error_code = 0
    dyn_model = "2-body"
    integration_step = 5.0_8
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Get the element type from the input flattened orbit.
    element_type_index = in_orbit(11)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    
    ! Create a StochastisOrbit instance form the flattened input orbit.
    call lslOrbitToStochasticOrbit(in_orbit,                                &
                                   in_covariance,                           &
                                   storb,                                   &
                                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now call the orbit propagation routine.
    call propagateOrbit(storb, dyn_model, integration_step, perturbers,  &
                        in_mjd, error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! And now re-flatten the StochasticOrbit.
    track_id = IDINT(in_orbit(1))
    call exportLslOrbit(storb, element_type, track_id, out_orbit,       &
                        out_covariance, out_sigmas, out_correlation,    &
                        error_code)
    if(error_code /= 0) then
        return
    end if
    call nullify(storb)
    deallocate(perturbers)
    return
end subroutine oorb_propagate_orbit_fast



subroutine oorb_propagate_orbit(in_orbit, in_covariance, in_mjd, out_orbit, &
                                out_covariance, out_sigmas,                 &
                                out_correlation, error_code)
    use oorb
    use StochasticOrbit_cl
    
    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(in)                   :: in_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(in)                 :: in_covariance
    ! The epoch (MJD TAI) we want to propagate to.
    real(8), intent(in)                                 :: in_mjd
    ! Output flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(out)                  :: out_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(out)                :: out_covariance
    REAL(8), DIMENSION(6), intent(out)                  :: out_sigmas
    REAL(8), DIMENSION(6,6), intent(out)                :: out_correlation
    ! Output error code
    integer, intent(out)                                :: error_code
    
    ! Internal variables.
    TYPE (StochasticOrbit)                              :: storb
    character(len=6)                                    :: dyn_model
    real(8)                                             :: integration_step
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    integer                                             :: track_id
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    
    
    ! Init
    error_code = 0
    dyn_model = "n-body"
    integration_step = 5.0_8
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Get the element type from the input flattened orbit.
    element_type_index = in_orbit(11)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    
    ! Create a StochastisOrbit instance form the flattened input orbit.
    call lslOrbitToStochasticOrbit(in_orbit,                                &
                                   in_covariance,                           &
                                   storb,                                   &
                                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now call the orbit propagation routine.
    call propagateOrbit(storb, dyn_model, integration_step, perturbers,  &
                        in_mjd, error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! And now re-flatten the StochasticOrbit.
    track_id = IDINT(in_orbit(1))
    call exportLslOrbit(storb, element_type, track_id, out_orbit,       &
                        out_covariance, out_sigmas, out_correlation,    &
                        error_code)
    if(error_code /= 0) then
        return
    end if
    call nullify(storb)
    deallocate(perturbers)
    return
end subroutine oorb_propagate_orbit



subroutine oorb_ephemeris_fast(in_orbit,                                &
                               in_covariance,                           &
                               in_obscode,                              &
                               in_num_ephems,                           &
                               in_step,                                 &
                               out_ephems,                              &
                               error_code)
    use oorb
    use StochasticOrbit_cl
    
    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(in)                   :: in_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(in)                 :: in_covariance
    character(len=4), intent(in)                        :: in_obscode
    ! Compute epehemeris from the orbit epoch to epoch+in_step*in_num_ephems
    ! Number of ephemeris to compute.
    integer, intent(in)                                 :: in_num_ephems
    ! Ephemeris step in fractional days.
    real(8), intent(in)                                 :: in_step
    ! Output ephemeris
    ! out_ephems = ((dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa), )
    real(8), dimension(in_num_ephems,10), intent(out)   :: out_ephems
    ! Output error code
    integer, intent(out)                                :: error_code

    ! Internal variables.
    TYPE (StochasticOrbit)                              :: storb
    character(len=6)                                    :: dyn_model
    real(8)                                             :: integration_step
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    TYPE (SphericalCoordinates), DIMENSION(:,:),POINTER :: ephem_arr
    REAL(8), DIMENSION(:,:,:), POINTER                  :: cov_arr
    REAL(8), DIMENSION(:,:), POINTER                    :: pdfs_arr
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: observers
    real(8)                                             :: timespan
    real(8)                                             :: step
    
    
    ! Init
    error_code = 0
    dyn_model = "2-body"
    integration_step = 5.0_8
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Get the element type from the input flattened orbit.
    element_type_index = in_orbit(11)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    ! Compute timespan given step and num_ephems.
    timespan = in_step * in_num_ephems
    step = in_step          ! Might be modified by the code below.
    
    ! Create a StochastisOrbit instance form the flattened input orbit.
    call lslOrbitToStochasticOrbit(in_orbit,                                &
                                   in_covariance,                           &
                                   storb,                                   &
                                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now call the orbit propagation routine.
    call ephemeris(storb,                                                   &
                   element_type,                                            &
                   dyn_model,                                               &
                   integration_step,                                        &
                   perturbers,                                              &
                   in_obscode,                                              &
                   timespan,                                                &
                   step,                                                    &
                   ephem_arr,                                               &
                   cov_arr,                                                 &
                   pdfs_arr,                                                &
                   observers,                                               &
                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now Export the ephem_arr to a flar array.
    call exportEphemeris(storb, ephem_arr, cov_arr, observers, out_ephems,  &
                         error_code)
    if(error_code /= 0) then
        return
    end if
    return
end subroutine oorb_ephemeris_fast



subroutine oorb_ephemeris(in_orbit,                                         &
                          in_covariance,                                    &
                          in_obscode,                                       &
                          in_num_ephems,                                    &
                          in_step,                                          &
                          out_ephems,                                       &
                          error_code)
    use oorb
    use StochasticOrbit_cl
    
    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(in)                   :: in_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(in)                 :: in_covariance
    character(len=4), intent(in)                        :: in_obscode
    ! Compute epehemeris from the orbit epoch to epoch+in_step*in_num_ephems
    ! Number of ephemeris to compute.
    integer, intent(in)                                 :: in_num_ephems
    ! Ephemeris step in fractional days.
    real(8), intent(in)                                 :: in_step
    ! Output ephemeris
    ! out_ephems = ((dist, ra, dec, mag, mjd, raErr, decErr, smaa, smia, pa), )
    real(8), dimension(in_num_ephems,10), intent(out)   :: out_ephems
    ! Output error code
    integer, intent(out)                                :: error_code

    ! Internal variables.
    TYPE (StochasticOrbit)                              :: storb
    character(len=6)                                    :: dyn_model
    real(8)                                             :: integration_step
    LOGICAL, DIMENSION(:), pointer                      :: perturbers
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    TYPE (SphericalCoordinates), DIMENSION(:,:),POINTER :: ephem_arr
    REAL(8), DIMENSION(:,:,:), POINTER                  :: cov_arr
    REAL(8), DIMENSION(:,:), POINTER                    :: pdfs_arr
    TYPE (CartesianCoordinates), DIMENSION(:), POINTER  :: observers
    real(8)                                             :: timespan
    real(8)                                             :: step
    
    
    ! Init
    error_code = 0
    dyn_model = "n-body"
    integration_step = 5.0_8
    allocate(perturbers(10), stat=error_code)
    if(error_code /= 0) then
        ! Error in allocating memory!
        error_code = 57
        return
    end if
    perturbers = .true.
    
    ! Get the element type from the input flattened orbit.
    element_type_index = in_orbit(11)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    ! Compute timespan given step and num_ephems.
    timespan = in_step * in_num_ephems
    step = in_step          ! Might be modified by the code below.
    
    ! Create a StochastisOrbit instance form the flattened input orbit.
    call lslOrbitToStochasticOrbit(in_orbit,                                &
                                   in_covariance,                           &
                                   storb,                                   &
                                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now call the orbit propagation routine.
    call ephemeris(storb,                                                   &
                   element_type,                                            &
                   dyn_model,                                               &
                   integration_step,                                        &
                   perturbers,                                              &
                   in_obscode,                                              &
                   timespan,                                                &
                   step,                                                    &
                   ephem_arr,                                               &
                   cov_arr,                                                 &
                   pdfs_arr,                                                &
                   observers,                                               &
                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now Export the ephem_arr to a flar array.
    call exportEphemeris(storb, ephem_arr, cov_arr, observers, out_ephems,  &
                         error_code)
    if(error_code /= 0) then
        return
    end if
    return
end subroutine oorb_ephemeris



subroutine oorb_moid(in_orbit, in_covariance, moid, error_code)
    use oorb
    use StochasticOrbit_cl
    
    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    real(8),dimension(11), intent(in)                   :: in_orbit
    ! Uncertainty matrices:
    REAL(8), DIMENSION(6,6), intent(in)                 :: in_covariance
    ! Output MOID
    real(8), intent(out)                                :: moid
    ! Output error code
    integer, intent(out)                                :: error_code
    
    ! Internal vars
    TYPE (StochasticOrbit)                              :: storb
    character(len=11)                                   :: element_type
    integer                                             :: element_type_index
    
    
    ! Init
    error_code = 0
    
    ! Get the element type from the input flattened orbit.
    element_type_index = in_orbit(11)
    if(element_type_index .le. 0 .or.                                       &
       element_type_index .gt. size(ORBITAL_ELEMENTS)) then
        ! Error: unsupported orbital elements.
        error_code = 58
        return
    end if
    element_type = ORBITAL_ELEMENTS(element_type_index)
    
    ! Create a StochastisOrbit instance form the flattened input orbit.
    call lslOrbitToStochasticOrbit(in_orbit,                                &
                                   in_covariance,                           &
                                   storb,                                   &
                                   error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now call the orbit propagation routine.
    call computeMoid(storb, moid, error_code)
    if(error_code /= 0) then
        return
    end if
    return
end subroutine oorb_moid



subroutine priv_classification(in_orbits, num_orbits, weights, error_code)
    use oorb
    use StochasticOrbit_cl
    
    ! Input/Output variables.
    ! Input flattened orbit:
    !  (track_id, elements(1:6), epoch, H, G, element_type_index)
    ! FIXME: H and G are not currently computed.
    integer, intent(in)                                 :: num_orbits
    ! Input flattened orbits.
    ! in_orbits has the form ((id, elements(1:6), epoch, Un-normalized p.d.f.,
    ! Reduced chi2, Regularized apr, Jacobian det(1:3), el_type_index), )
    real(8),dimension(num_orbits,15),intent(in)         :: in_orbits
    ! Class names and probability arrays.
    ! CHARACTER(len=16), DIMENSION(17), intent(out)       :: classes
    REAL(8), DIMENSION(17), intent(out)                 :: weights
    ! Output error code
    integer, intent(out)                                :: error_code
    
    ! Internal vars
    TYPE (StochasticOrbit)                              :: storb
    integer                                             :: j = 0
    REAL(bp), DIMENSION(:), POINTER                     :: weight_arr
    CHARACTER(len=16), DIMENSION(:), POINTER            :: group_name_arr
    
    
    ! Init
    error_code = 0
    
    ! Create a StochastisOrbit instance form the flattened input orbits.
    call nullify(storb)
    call rangingOrbitsToStochasticOrbit(in_orbits,                          &
                                        num_orbits,                         &
                                        storb,                              &
                                        error_code)
    if(error_code /= 0) then
        return
    end if
    
    ! Now call the orbit propagation routine.
    call classification(storb, group_name_arr, weight_arr, error_code)
    if(error_code /= 0) then
        return
    end if
    call nullify(storb)
    
    ! Copy values from the pointers to the regular arrays.
    do j=1,min(size(weights), size(weight_arr))
!        classes(j) = group_name_arr(j)
        weights(j) = weight_arr(j)
!        write(*, *) "classes(j): ", classes(j)
!        write(*, *) "weights(j): ", weights(j)
    end do
    deallocate(group_name_arr)
    deallocate(weight_arr)
    return
end subroutine priv_classification



subroutine oorb_calendardate_to_mjd(year, month, day, mjd, timescale, &
                                    error_code)
    use oorb, only: calendarDateToMjd
    
    ! Input: year, month, fractional days in whatever timescale they are 
    ! (usually UTC or UT1). Output MJD TAI.
    
    integer, intent(in)                             :: year
    integer, intent(in)                             :: month
    ! day + time in fractional days.
    real(8), intent(in)                             :: day
    CHARACTER(len=*), intent(in)                    :: timescale
    real(8), intent(out)                            :: mjd
    integer, intent(out)                            :: error_code
    
    
    call calendarDateToMjd(year, month, day, mjd, trim(timescale), error_code)
    return
end subroutine oorb_calendardate_to_mjd



subroutine oorb_mjdutc_to_mjdtai(mjd_utc, mjd_tai, error_code)
    use oorb, only: mjdConvert
    
    ! Input:  MJD UTC
    ! Output: MJD TAI
    real(8), intent(in)                             :: mjd_utc
    real(8), intent(out)                            :: mjd_tai
    integer, intent(out)                            :: error_code
    
    
    call mjdConvert(mjd_utc, "UTC", mjd_tai, "TAI", error_code)
    return
end subroutine oorb_mjdutc_to_mjdtai



subroutine oorb_mjdtai_to_mjdutc(mjd_tai, mjd_utc, error_code)
    use oorb, only: mjdConvert
    
    ! Input:  MJD UTC
    ! Output: MJD TAI
    real(8), intent(in)                             :: mjd_tai
    real(8), intent(out)                            :: mjd_utc
    integer, intent(out)                            :: error_code
    
    
    call mjdConvert(mjd_tai, "TAI", mjd_utc, "UTC", error_code)
    return
end subroutine oorb_mjdtai_to_mjdutc



subroutine oorb_mjdtt_to_mjdtai(mjd_tt, mjd_tai, error_code)
    use oorb, only: mjdConvert
    
    ! Input:  MJD TT = MJD TAI + 32.184 seconds
    ! Output: MJD TAI
    real(8), intent(in)                             :: mjd_tt
    real(8), intent(out)                            :: mjd_tai
    integer, intent(out)                            :: error_code
    
    
    call mjdConvert(mjd_tt, "TT", mjd_tai, "TAI", error_code)
    return
end subroutine oorb_mjdtt_to_mjdtai



subroutine oorb_mjdtai_to_mjdtt(mjd_tai, mjd_tt, error_code)
    use oorb, only: mjdConvert
    
    ! Input:  MJD TAI
    ! Output: MJD TT = MJD TAI + 32.184 seconds
    real(8), intent(in)                             :: mjd_tai
    real(8), intent(out)                            :: mjd_tt
    integer, intent(out)                            :: error_code
    
    
    call mjdConvert(mjd_tai, "TAI", mjd_tt, "TT", error_code)
    return
end subroutine oorb_mjdtai_to_mjdtt



subroutine sexagesimal_hours_to_decimal_degrees(h, m, s, deg, error_code)
    integer, intent(in)                             :: h
    integer, intent(in)                             :: m
    real(8), intent(in)                             :: s
    real(8), intent(out)                            :: deg
    integer, intent(out)                            :: error_code
    
    
    error_code = 0
    call sexagesimal_degrees_to_decimal_degrees(h, m, s, deg, error_code)
    deg = deg * 15.0_8
    return
end subroutine sexagesimal_hours_to_decimal_degrees



subroutine sexagesimal_degrees_to_decimal_degrees(d, m, s, deg, error_code)
    integer, intent(in)                             :: d
    integer, intent(in)                             :: m
    real(8), intent(in)                             :: s
    real(8), intent(out)                            :: deg
    integer, intent(out)                            :: error_code

    
    error_code = 0
    if(d .ge. 0.0_8) then
        deg = d + m / 60.0_8 + s / 3600.0_8
    else
        deg = d - m / 60.0_8 - s / 3600.0_8
    end if
    return
end subroutine sexagesimal_degrees_to_decimal_degrees




















