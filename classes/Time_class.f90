!====================================================================!
!                                                                    !
! Copyright 2002-2014,2015                                           !
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
!! Type and routines for time.
!!
!! @see CartesianCoordinates_class
!! @see SphericalCoordinates_class
!! @see Orbit_class
!! 
!! @author  MG, JV
!! @version 2015-10-23
!!
MODULE Time_cl

  USE Base_cl
  USE File_cl

  USE utilities

  IMPLICIT NONE
  ! Error limit for iterations:
  REAL(bp), PARAMETER, PRIVATE                 :: epst = 1.0e-10_bp
  ! Max number of iterations:
  INTEGER, PARAMETER, PRIVATE                  :: nitmax = 5
  INTEGER, PARAMETER, PRIVATE                  :: nitutc = 20
  REAL(bp), DIMENSION(:), POINTER, PRIVATE     :: tv => NULL()
  REAL(bp), DIMENSION(:), POINTER, PRIVATE     :: dtv => NULL()
  INTEGER, DIMENSION(:), POINTER, PRIVATE      :: mjdv => NULL()
  INTEGER, DIMENSION(:), POINTER, PRIVATE      :: idv => NULL()
  INTEGER, PRIVATE                             :: jp
  INTEGER, PRIVATE                             :: ipos
  INTEGER, PRIVATE                             :: taiut_size
  INTEGER, PRIVATE                             :: etut_size
  LOGICAL, PRIVATE                             :: first = .TRUE.

  PRIVATE :: new_T
  PRIVATE :: new_T_cd_long
  PRIVATE :: new_T_cd_short
  PRIVATE :: new_T_mjd
  PRIVATE :: nullify_T  
  PRIVATE :: copy_T  
  PRIVATE :: exist_T  
  PRIVATE :: equal_T  
  PRIVATE :: deltaAT
  PRIVATE :: deltaT
  PRIVATE :: getCalendarDate_long
  PRIVATE :: getCalendarDate_short
  PRIVATE :: getCurrentTime_values
  PRIVATE :: getMJD_t
  PRIVATE :: getMJD_cd
  PRIVATE :: reallocate_T_1
  PRIVATE :: timescaleConversion
  PRIVATE :: toNormalForm

  TYPE Time
     PRIVATE
     REAL(bp) :: tdt            =  0.0_bp
     REAL(bp) :: utc            =  0.0_bp
     REAL(bp) :: tai            =  0.0_bp
     REAL(bp) :: ut1            =  0.0_bp
     LOGICAL  :: is_initialized = .FALSE.
  END TYPE Time

  INTERFACE NEW
     MODULE PROCEDURE new_T
     MODULE PROCEDURE new_T_cd_long
     MODULE PROCEDURE new_T_cd_short
     MODULE PROCEDURE new_T_mjd
     MODULE PROCEDURE new_T_MPC
  END INTERFACE NEW

  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_T
  END INTERFACE NULLIFY

  INTERFACE copy
     MODULE PROCEDURE copy_T
  END INTERFACE copy

  INTERFACE exist
     MODULE PROCEDURE exist_T
  END INTERFACE exist

  INTERFACE equal
     MODULE PROCEDURE equal_T
  END INTERFACE equal

  INTERFACE getCalendarDate
     MODULE PROCEDURE getCalendarDate_long
     MODULE PROCEDURE getCalendarDate_short
  END INTERFACE getCalendarDate

  INTERFACE getCurrentTime
     MODULE PROCEDURE getCurrentTime_values
  END INTERFACE getCurrentTime

  INTERFACE getMJD
     MODULE PROCEDURE getMJD_t
     MODULE PROCEDURE getMJD_cd
  END INTERFACE getMJD

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_T_1
  END INTERFACE reallocate

CONTAINS




  !! *Desription*:
  !!
  !! Initializes a new object using default values.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_T(this)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout) :: this
    TYPE (Time) :: t

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (first) THEN
       CALL NEW(t, 0.0_bp, "tdt")
       CALL NULLIFY(t)
    END IF

    this%tdt = -1.0_bp
    this%ut1 = -1.0_bp
    this%utc = -1.0_bp
    this%tai = -1.0_bp

    this%is_initialized = .TRUE.

  END SUBROUTINE new_T





  !! *Desription*:
  !!
  !! Initializes a new object using values given as a calendar 
  !! date.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_T_cd_long(this, year, month, day, hour, min, sec, timescale)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    INTEGER, INTENT(in)          :: year, month, day, hour, min
    REAL(bp), INTENT(in)         :: sec
    CHARACTER(len=*), INTENT(in) :: timescale
    REAL(bp)                     :: mjd, tmp

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    tmp = REAL(day,bp) + (REAL(hour,bp) + REAL(min,bp)/60.0_bp + sec/3600.0_bp)/24.0_bp
    mjd = getMJD(year, month, tmp)
    CALL NEW(this, mjd, timescale)
    IF (error) THEN
       CALL errorMessage("Time / new", &
            "TRACE BACK", 1)
       RETURN
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE new_T_cd_long





  !! *Desription*:
  !!
  !! Initializes a new object using values given as a calendar 
  !! date.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_T_cd_short(this, year, month, day, timescale)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    INTEGER, INTENT(in)          :: year, month
    REAL(bp), INTENT(in)         :: day
    CHARACTER(len=*), INTENT(in) :: timescale
    REAL(bp)                     :: mjd

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF
    mjd = getMJD(year, month, day)
    CALL NEW(this, mjd, timescale)
    IF (error) THEN
       CALL errorMessage("Time / new", &
            "TRACE BACK", 1)
       RETURN
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE new_T_cd_short





  !! *Desription*:
  !!
  !! Initializes a new object using values 
  !! given as a MPC packed date.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_T_MPC(this, packed_date, timescale)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: packed_date
    CHARACTER(len=*), INTENT(in) :: timescale
    CHARACTER(len=31), PARAMETER :: coding = &
         "123456789ABCDEFGHIJKLMNOPQRSTUV"
    INTEGER                      :: year, month
    REAL(bp)                     :: day
    INTEGER                      :: tmp, n

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    tmp = INDEX(coding, packed_date(1:1))
    CALL toInt(packed_date(2:3), year, error)
    IF (error) THEN
       CALL errorMessage("Time / new", &
            "Could not transform string to integer (1).", 1)
       RETURN
    END IF
    year = tmp*100 + year
    month = INDEX(coding, packed_date(4:4))
    day = REAL(INDEX(coding, packed_date(5:5)),bp)
    n = LEN_TRIM(packed_date)
    IF (n > 5) THEN
       CALL toInt(packed_date(6:n), tmp, error)
       IF (error) THEN
          CALL errorMessage("Time / new", &
               "Could not transform string to integer (2).", 1)
          RETURN
       END IF
       day = day + tmp/10.0_bp**(n-5)
    END IF

    CALL NEW(this, year, month, day, timescale)
    IF (error) THEN
       CALL errorMessage("Time / new", &
            "TRACE BACK", 1)
       RETURN
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE new_T_MPC





  !! *Description*:
  !!
  !! Initializes a new object using values given as a modified 
  !! Julian date.
  !!
  !! Returns error.
  !!
  SUBROUTINE new_T_mjd(this, mjd, timescale)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    REAL(bp), INTENT(in)         :: mjd
    CHARACTER(len=*), INTENT(in) :: timescale
    TYPE (File)                  :: datafile
    CHARACTER(len=10)            :: record
    REAL(bp)                     :: dt, iii
    INTEGER                      :: err, year, month, day, &
         taiut_max, etut_max

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (first) THEN

       ! Load table of TAI-UTC as a function of the day:
       CALL NEW(datafile, TRIM(OORB_DATA_DIR) // &
            "/" // TRIM(TAIUTC_FNAME))
       CALL setActionRead(datafile)
       CALL setStatusOld(datafile)
       CALL OPEN(datafile)
       IF (error) THEN
          CALL errorMessage("Time / new", &
               "TRACE BACK 3", 1)
          RETURN
       END IF
       taiut_max = getNrOfLines(datafile)
       ALLOCATE(idv(taiut_max), mjdv(taiut_max), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Time / new", &
               "Could not allocate arrays for TAI-UT.", 1)
          RETURN
       END IF
       DO
          READ(getUnit(datafile), "(A)", iostat=err) record
          IF (err <  0) THEN
             error = .TRUE.
             CALL errorMessage("Time / new", &
                  "End of file " // TRIM(TAIUTC_FNAME), 1)
             RETURN
          ELSE IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Time / new", &
                  "Could not read file " // TRIM(TAIUTC_FNAME), 1)
             RETURN
          ELSE
             IF (record == "----------") EXIT
          END IF
       END DO
       taiut_size = 0
       DO
          READ(getUnit(datafile),*,iostat=err) day, month, year, iii
          IF (err <  0) THEN
             EXIT
          ELSE IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Time / new", &
                  "Could not read file " // TRIM(TAIUTC_FNAME), 1)
             RETURN
          ELSE
             taiut_size = taiut_size + 1
             idv(taiut_size) = iii
             mjdv(taiut_size) = NINT(getMJD(year,month,day*1.0_bp))
             IF (taiut_size > 1) THEN
                IF (mjdv(taiut_size) <= mjdv(taiut_size-1)) THEN
                   error  = .TRUE.
                   CALL errorMessage("Time / new", &
                        "File " // TRIM(TAIUTC_FNAME) // " is not sorted.", 1)
                   RETURN
                END IF
             END IF
          END IF
       END DO
       jp = 2
       CALL NULLIFY(datafile)
       IF (error) THEN
          CALL errorMessage("Time / new", &
               "TRACE BACK 4", 1)
          RETURN
       END IF
       IF (taiut_size < 2) THEN
          error = .TRUE.
          CALL errorMessage("Time / new", &
               "File " // TRIM(TAIUTC_FNAME) // &
               " contains less than two data records.", 1)
          RETURN
       END IF

       ! Load table of ET-UT as a function of UT:
       CALL NEW(datafile, TRIM(OORB_DATA_DIR) // "/" &
            // TRIM(ETUT_FNAME))
       CALL setActionRead(datafile)
       CALL setStatusOld(datafile)
       CALL OPEN(datafile)
       IF (error) THEN
          CALL errorMessage("Time / new", &
               "TRACE BACK 7", 1)
          RETURN
       END IF
       etut_max = getNrOfLines(datafile)
       ALLOCATE(tv(etut_max), dtv(etut_max), stat=err)
       IF (err /= 0) THEN
          error = .TRUE.
          CALL errorMessage("Time / new", &
               "Could not allocate arrays for ET-UT.", 1)
          RETURN
       END IF
       DO
          READ(getUnit(datafile), "(A)", iostat=err) record
          IF (err <  0) THEN
             error = .TRUE.
             CALL errorMessage("Time / new", &
                  "End of file" // TRIM(ETUT_FNAME), 1)
             RETURN
          ELSE IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Time / new", &
                  "Could not read file " // TRIM(ETUT_FNAME), 1)
             RETURN
          ELSE
             IF (record == "----------") EXIT
          END IF
       END DO
       etut_size = 0
       DO
          READ(getUnit(datafile),*,iostat=err) day, month, year, dt
          IF (err <  0) THEN
             EXIT
          ELSE IF (err > 0) THEN
             error = .TRUE.
             CALL errorMessage("Time / new", &
                  "Could not read file " // TRIM(ETUT_FNAME), 1)
             RETURN
          ELSE
             etut_size = etut_size + 1
             tv(etut_size) = getMJD(year, month, day*1.0_bp)
             dtv(etut_size) = dt
          END IF
       END DO
       ipos = 1
       CALL NULLIFY(datafile)
       IF (etut_size < 2) THEN
          error = .TRUE.
          CALL errorMessage("Time / new", &
               "File " // TRIM(ETUT_FNAME) // &
               " contains less than two data records.", 1)
          RETURN
       END IF

       first = .FALSE.

    END IF

    SELECT CASE (TRIM(timescale))
    CASE ("ut1", "UT1")
       this%ut1 = mjd
       CALL timescaleConversion(mjd, "UT1", this%tdt, "TDT")
       this%utc = -1.0_bp
       this%tai = -1.0_bp
    CASE ("tai", "TAI")
       this%tai = mjd
       CALL timescaleConversion(mjd, "TAI", this%tdt, "TDT")
       this%ut1 = -1.0_bp
       this%utc = -1.0_bp
    CASE ("utc", "UTC")
       this%utc = mjd
       CALL timescaleConversion(mjd, "UTC", this%tdt, "TDT")
       this%tai = -1.0_bp
       this%ut1 = -1.0_bp
    CASE ("tdt", "TDT", "et", "ET", "tt", "TT")
       this%tdt  = mjd
       this%utc = -1.0_bp
       this%tai = -1.0_bp
       this%ut1 = -1.0_bp
    CASE ("tcb", "TCB")
       ! This assumes TT == TDB (error should be <2 millisecs)
       this%tdt = mjd - 1.550505e-8_bp * (mjd - 43144.0_bp)
       this%utc = -1.0_bp
       this%tai = -1.0_bp
       this%ut1 = -1.0_bp
    CASE default
       error = .TRUE.
       CALL errorMessage("Time / new", &
            "Type of time is either missing or erroneus.", 1)
       RETURN
    END SELECT
    IF (error) THEN
       CALL errorMessage("Time / new", &
            "TRACE BACK 10", 1)
       RETURN
    END IF
    this%is_initialized = .TRUE.

  END SUBROUTINE new_T_mjd





  !! *Description*:
  !!
  !! Nullifies this object.
  !!
  SUBROUTINE nullify_T(this)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout) :: this

    this%tdt            = -1.0_bp
    this%utc            = -1.0_bp
    this%tai            = -1.0_bp
    this%ut1            = -1.0_bp
    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_T





  !! *Description*:
  !!
  !! Deallocates the memory required for storing time-difference tables.
  !!
  SUBROUTINE nullifyTime()

    IMPLICIT NONE
    INTEGER :: err

    IF (ASSOCIATED(tv)) THEN
       DEALLOCATE(tv, stat=err)
    END IF
    IF (ASSOCIATED(dtv)) THEN
       DEALLOCATE(dtv, stat=err)
    END IF
    IF (ASSOCIATED(mjdv)) THEN
       DEALLOCATE(mjdv, stat=err)
    END IF
    IF (ASSOCIATED(idv)) THEN
       DEALLOCATE(idv, stat=err)
    END IF

  END SUBROUTINE nullifyTime





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_T(this)

    IMPLICIT NONE
    TYPE (Time), INTENT(in) :: this
    TYPE (Time)             :: copy_T

    copy_T%tdt            = this%tdt
    copy_T%utc            = this%utc
    copy_T%ut1            = this%ut1
    copy_T%tai            = this%tai
    copy_T%is_initialized = this%is_initialized

  END FUNCTION copy_T





  !! *Description*:
  !!
  !! Returns the status of this object, i.e. whether
  !! it exists or not.
  !!
  LOGICAL FUNCTION exist_T(this)

    IMPLICIT NONE
    TYPE (Time), INTENT(in) :: this

    exist_T = this%is_initialized

  END FUNCTION exist_T





  !! *Description*:
  !!
  !! Returns difference DAT = TAI - UTC as a function of UTC.
  !!
  !! INPUT:    mjdc    - Modified Julian Day (UTC)
  !!
  !! OUTPUT:   deltaAT - TAI-UTC given as an integer number
  !!           of seconds
  !!
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  !! Returns error.
  !!
  INTEGER FUNCTION deltaAT(mjdc)

    IMPLICIT NONE
    INTEGER, INTENT(in)  :: mjdc
    INTEGER              :: i

    ! Trying to use previous value
    IF(mjdc >= mjdv(jp-1) .AND. mjdc < mjdv(jp)) THEN
       deltaAT = idv(jp-1)
       RETURN
    END IF

    ! Selecting the records of the table before and after the date
    ! supplied
    IF (mjdc < mjdv(1)) THEN
       error = .TRUE.
       CALL errorMessage("Time / deltaAT", &
            "TJM too small.", 1)
       RETURN
    ELSE IF (mjdc >= mjdv(taiut_size)) THEN
       error = .TRUE.
       CALL errorMessage("Time / deltaAT", &
            "TJM too large.", 1)
       RETURN
    ELSE
       DO i=2, taiut_size
          IF (mjdc < mjdv(i)) THEN
             jp = i
             deltaAT = idv(jp-1)
             RETURN
          END IF
       END DO
    END IF

  END FUNCTION deltaAT





  !! Description:
  !!
  !! Returns delta T = ET - UT.
  !!
  !! INPUT:    tjm       -  Modified Julian Day (UT1)
  !!
  !! OUTPUT:   deltaT    -  ET - UT1 (in seconds)
  !!
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION deltaT(tjm)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: tjm
    REAL(bp)             :: c1, c2

    IF (.NOT.ASSOCIATED(dtv)) THEN
       error = .TRUE.
       CALL errorMessage("Time / deltaT", &
            "dtv has not been initialized.", 1)
       RETURN
    END IF

    ! Trying to use previous value
    IF (tjm >= tv(ipos) .AND. tjm <= tv(ipos+1)) THEN
       c1 = (tv(ipos+1) - tjm) / (tv(ipos+1) - tv(ipos))
       c2 = 1 - c1
       deltaT = c1*dtv(ipos) + c2*dtv(ipos+1)
       RETURN
    END IF
    ! Selecting the records of the table before and after the date
    ! supplied
    IF (tjm < tv(1) .OR. tjm > tv(etut_size)) THEN
       error = .TRUE.
       CALL errorMessage("Time / deltaT", &
            "TJM out of range", 1)
       RETURN
    ELSE
       DO ipos=1,etut_size-1
          IF (tjm >= tv(ipos) .AND. tjm <= tv(ipos+1)) THEN
             c1 = (tv(ipos+1) - tjm) / (tv(ipos+1) - tv(ipos))
             c2 = 1 - c1
             deltaT = c1*dtv(ipos) + c2*dtv(ipos+1)
             RETURN
          END IF
       END DO
       error = .TRUE.
       CALL errorMessage("Time / deltaT", &
            "Internal error.", 1)
       RETURN
    END IF

  END FUNCTION deltaT





  !! *Description*:
  !!
  !! Returns .true., if the absolute difference of the
  !! objects is less than epsilon, and .false. otherwise.
  !!
  LOGICAL FUNCTION equal_T(this, that)

    IMPLICIT NONE
    TYPE (Time), INTENT(in) :: this, that

    IF (ABS(this%tdt-that%tdt) < EPSILON(this%tdt)) THEN
       equal_T = .TRUE.
    ELSE
       equal_T = .FALSE.
    END IF

  END FUNCTION equal_T





  !! *Description*:
  !!
  !! Returns
  !! 
  !!  - TDT (Terrestrial Dynamical Time; 
  !!          = ET (Ephemeris Time) 
  !!          = TT (Terrestrial Time)), or
  !!  - TAI (Atomic Time, french Temps Atomique International), or
  !!  - UTC (Coordinated Universal Time), or
  !!  - UT1 (UT + pole variation)
  !!
  !! as a Calendar Date. Default timescale is UTC.
  !!
  !! Returns error.
  !! 
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  SUBROUTINE getCalendarDate_long(this, timescale, year, month, day, hour, min, sec)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: timescale
    INTEGER, INTENT(out)         :: year, month, day, hour, min
    REAL(bp), INTENT(out)        :: sec
    REAL(bp)                     :: aa, mjd, tmp
    INTEGER                      :: a, b, c, d, e, f

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / getCalendarDate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd = getMJD(this, timescale)
    IF (error) THEN
       CALL errorMessage("Time / getCalendarDate", &
            "TRACE BACK", 1)
       RETURN
    END IF

    aa = mjd + 2400001.0_bp
    a = INT(aa)
    IF (a < 2299161) THEN
       c = a + 1524
    ELSE
       b = FLOOR((a-1867216.25_bp)/36524.25_bp)
       c = a + b - FLOOR(b/4.0_bp) + 1525
    END IF
    d = FLOOR((c-122.1_bp)/365.25_bp)
    e = FLOOR(365.25_bp*d)
    f = FLOOR((c-e)/30.6001_bp)
    tmp = c - e - FLOOR(30.6001_bp*f) + (aa - a)
    ! date:
    day = FLOOR(tmp)
    month = f - 1 - 12*FLOOR(f/14.0_bp)
    year = d - 4715 - FLOOR((7+month)/10.0_bp)
    ! time:
    tmp = 24.0_bp*(tmp - day)
    hour = FLOOR(tmp)
    tmp = 60.0_bp*(tmp - hour)
    min = FLOOR(tmp)
    sec = 60.0_bp*(tmp - min)

  END SUBROUTINE getCalendarDate_long





  !! *Description*:
  !!
  !! Returns
  !! 
  !!  - TDT (Terrestrial Dynamical Time; 
  !!          = ET (Ephemeris Time) 
  !!          = TT (Terrestrial Time)), or
  !!  - TAI (Atomic Time, french Temps Atomique International), or
  !!  - UTC (Coordinated Universal Time), or
  !!  - UT1 (UT + pole variation)
  !!
  !! as a Calendar Date. Default timescale is UTC.
  !!
  !! Returns error.
  !! 
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  SUBROUTINE getCalendarDate_short(this, timescale, year, month, day)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: timescale
    INTEGER, INTENT(out)         :: year, month
    REAL(bp), INTENT(out)        :: day
    REAL(bp)                     :: aa, mjd
    INTEGER                      :: a, b, c, d, e, f

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / getCalendarDate", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd = getMJD(this, timescale)
    IF (error) THEN
       CALL errorMessage("Time / getCalendarDate", &
            "TRACE BACK", 1)
       RETURN
    END IF

    aa = mjd + 2400001.0_bp
    a = INT(aa)
    IF (a < 2299161) THEN
       c = a + 1524
    ELSE
       b = FLOOR((a-1867216.25_bp)/36524.25_bp)
       c = a + b - FLOOR(b/4.0_bp) + 1525
    END IF
    d = FLOOR((c-122.1_bp)/365.25_bp)
    e = FLOOR(365.25_bp*d)
    f = FLOOR((c-e)/30.6001_bp)
    day = c - e - FLOOR(30.6001_bp*f) + (aa - a)
    month = f - 1 - 12*FLOOR(f/14.0_bp)
    year = d - 4715 - FLOOR((7+month)/10.0_bp)

  END SUBROUTINE getCalendarDate_short





  !! *Description*:
  !!
  !! Returns the calendar date as a string (_yyyy_mm_dd.dddddd_).
  !!
  !! Returns error.
  !!
  CHARACTER(len=17) FUNCTION getCalendarDateString(this, timescale) RESULT(str)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: timescale
    REAL(bp)                     :: day, sec
    INTEGER                      :: err, year, month, day_, hour, min

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / getCalendarDateString", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    CALL getCalendarDate(this, timescale, year, month, day_, hour, min, sec)
    IF (error) THEN
       CALL errorMessage("Time / getCalendarDateString", &
            "TRACE BACK", 1)
       RETURN
    END IF

    day = day_ + (hour + min/60.0_bp + sec/3600.0_bp)/24.0_bp
    WRITE(str, "(I4,1X,I2,1X,F9.6)", iostat=err) year, month, day
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Time / getCalendarDateString", &
            "Could not write output string.", 1)
       str = "*****ERROR*****"
       RETURN
    END IF

  END FUNCTION getCalendarDateString






  !! *Description*:
  !!
  !! This routine returns the current date and time of any timezone
  !! given the difference between the needed timezone and the timezone
  !! of the computer (e.g., for a Helsinki based computer timezone
  !! should be set to -2 to get Greenwich time).
  !!
  !! Programming interface:
  !!
  !! integer getCurrentTime(1) is day
  !! integer getCurrentTime(2) is month (Jan=1,...,Dec=12)
  !! integer getCurrentTime(3) is year
  !! integer getCurrentTime(4) is hour
  !! integer getCurrentTime(5) is minute
  !! integer getCurrentTime(6) is second
  !!
  !! Returns error.
  !!
  FUNCTION getCurrentTime_values(timezone)

    IMPLICIT NONE
    INTEGER, DIMENSION(6) :: getCurrentTime_values
    REAL(bp), INTENT(in)  :: timezone
    TYPE (Time)           :: t
    REAL(kind=bp)         :: mjd, s, day
    INTEGER, DIMENSION(8) :: date_time

    IF (first) THEN
       CALL NEW(t, 54000.0_bp, "tdt")
       CALL NULLIFY(t)
    END IF

    CALL DATE_AND_TIME(values=date_time)
    day = REAL(date_time(3),bp) + (REAL(date_time(5),bp) + &
         REAL(date_time(6),bp)/60.0_bp + &
         REAL(date_time(7),bp)/3600.0_bp)/hour_day
    mjd = getMJD(date_time(1), date_time(2), day) + timezone/24.0_bp
    CALL NEW(t, mjd, "utc")
    IF (error) THEN
       CALL errorMessage("Time / getCurrentTime", "TRACE BACK (5)", 1)
       RETURN
    END IF
    CALL getCalendarDate(t, "utc", date_time(1), date_time(2), &
         date_time(3), date_time(5), date_time(6), s)
    IF (error) THEN
       CALL errorMessage("Time / getCurrentTime", "TRACE BACK (10)", 1)
       RETURN
    END IF
    date_time(7) = NINT(s)
    CALL NULLIFY(t)

    ! If seconds rounds to 60 then add half a second to the Julian
    ! date and recalculate the calendar date:
    IF (date_time(7) == 60) THEN
       mjd = mjd + 0.5d0/86400.0d0
       CALL NEW(t, mjd, "utc")
       IF (error) THEN
          CALL errorMessage("Time / getCurrentTime", "TRACE BACK (15)", 1)
          RETURN
       END IF
       CALL getCalendarDate(t, "utc", date_time(1), date_time(2), &
            date_time(3), date_time(5), date_time(6), s)
       IF (error) THEN
          CALL errorMessage("Time / getCurrentTime", "TRACE BACK (20)", 1)
          RETURN
       END IF
       date_time(7) = NINT(s)
       CALL NULLIFY(t)
    END IF
    getCurrentTime_values(1:3) = date_time(1:3)
    getCurrentTime_values(4:6) = date_time(5:7)

  END FUNCTION getCurrentTime_values





  !!
  !! Greenwich Mean Sidereal Time as a function of UT1        
  !!
  !! INPUT:    TJM       -  Modified Julian Time (UT1)
  !!
  !! OUTPUT:   GMST      -  Greenwich Mean Sidereal Time referred
  !!                        to the mean equinox of date (rad)
  !!
  !! Based on routines by Mario Carpino (carpino@brera.mi.astro.it)
  !!
  REAL(bp) FUNCTION getGMST(this)

    IMPLICIT NONE

    TYPE(Time), INTENT(inout) :: this
    INTEGER                :: imjd,i
    REAL(bp)               :: t0,mjd,gmst0,gmst,h
    REAL(bp), PARAMETER    :: c0  = 24110.54841_bp, c1 = 8640184.812866_bp,&
         c2  = 9.3104e-2_bp,   c3 =-6.2e-6_bp, &
         rap = 1.00273790934_bp

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / getGMST", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd = getMJD(this, "UT1")
    ! Sidereal time at 0h UT1
    imjd = mjd
    t0 = (imjd-51544.5_bp)/36525.0_bp
    gmst0 = ((c3*t0+c2)*t0+c1)*t0 + c0
    gmst0 = gmst0*two_pi/86400.0_bp

    ! Increment in GMST from 0h
    h = (mjd-imjd)*two_pi
    gmst = gmst0 + h*rap
    i = gmst/two_pi
    IF (gmst < 0.0_bp) THEN
       i = i - 1
    END IF
    gmst = gmst - i*two_pi

    getGMST = gmst

  END FUNCTION getGMST





  !! *Description*:
  !!
  !! Returns
  !! 
  !!  - TDT (Terrestrial Dynamical Time; 
  !!          = ET (Ephemeris Time) 
  !!          = TT (Terrestrial Time)), or
  !!  - TAI (Atomic Time, french Temps Atomique International), or
  !!  - UTC (Coordinated Universal Time), or
  !!  - UT1 (UT + pole variation)
  !! 
  !! of this Time object as a Julian Date.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getJD(this, timescale) RESULT(jd)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: timescale
    REAL(bp)                     :: mjd

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / getJD", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    mjd = getMJD(this, timescale)
    IF (error) THEN
       CALL errorMessage("Time / getJD", &
            "TRACE BACK", 1)
       RETURN
    END IF
    jd = mjd + 2400000.5_bp 

  END FUNCTION getJD





  !! *Description*:
  !!
  !! Returns
  !! 
  !!  - TDT (Terrestrial Dynamical Time; 
  !!          = ET (Ephemeris Time) 
  !!          = TT (Terrestrial Time)), or
  !!  - TAI (Atomic Time, french Temps Atomique International), or
  !!  - UTC (Coordinated Universal Time), or
  !!  - UT1 (UT + pole variation)
  !! 
  !! of this Time object as a Modified Julian Date.
  !!
  !! Returns error.
  !!
  REAL(bp) FUNCTION getMJD_t(this, timescale) RESULT(mjd)

    IMPLICIT NONE
    TYPE (Time), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: timescale

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("Time / getMJD", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    SELECT CASE (TRIM(timescale))
    CASE ("ut1", "UT1")
       IF (this%ut1 < 0.0_bp) THEN
          CALL timescaleConversion(this%tdt, "TDT", mjd, "UT1")
          IF (error) THEN
             CALL errorMessage("Time / getMJD", &
                  "TRACE BACK (5)", 1)
             RETURN
          END IF
          this%ut1 = mjd
       ELSE
          mjd = this%ut1
       END IF
    CASE ("tai", "TAI")
       IF (this%tai < 0.0_bp) THEN
          CALL timescaleConversion(this%tdt, "TDT", mjd, "TAI")
          IF (error) THEN
             CALL errorMessage("Time / getMJD", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
          this%tai = mjd
       ELSE
          mjd = this%tai
       END IF
    CASE ("utc", "UTC")
       IF (this%utc < 0.0_bp) THEN
          CALL timescaleConversion(this%tdt, "TDT", mjd, "UTC")
          IF (error) THEN
             CALL errorMessage("Time / getMJD", &
                  "TRACE BACK (15)", 1)
             RETURN
          END IF
          this%utc = mjd
       ELSE
          mjd = this%utc
       END IF
    CASE ("tdt", "TDT", "et", "ET", "tt", "TT")
       mjd = this%tdt
    CASE default
       error = .TRUE.
       CALL errorMessage("Time / getMJD", &
            "Timescale " // TRIM(timescale) // " is erroneus.", 1)
       RETURN
    END SELECT

  END FUNCTION getMJD_t





  !! *Description*:
  !!
  !! Computes the Modified Julian Date from Calendar Date.
  !!
  !! INPUT:    year   -  Year (e.g.: 1987)
  !!           month  -  Month of the year ( 1 <= month <= 12 )
  !!           day    -  Day of the month ( 1.0 <= day < 32.0 )
  !!
  !! OUTPUT:   getMJD -  Modified Julian Day MJD = JD - 2,400,000.5
  !!
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  REAL(bp) FUNCTION getMJD_cd(year, month, day)

    IMPLICIT NONE
    INTEGER, INTENT(in)  :: month, year
    REAL(bp), INTENT(in) :: day
    INTEGER              :: k1, k2, ib, month_, year_

    year_ = year
    month_ = month
    DO WHILE (month_ > 12)
       month_ = month_ - 12
       year_ = year_ + 1
    END DO
    DO WHILE (month_ <= 2)
       year_ = year_ - 1
       month_ = month_ + 12
    END DO
    IF (year_ > 1582) THEN
       ib = year_/400 - year_/100
    ELSE
       ib = -2
       IF (year_ == 1582) THEN
          IF (month_ > 10) THEN
             ib = year_/400 - year_/100
          ELSE IF (month_ == 10 .AND. day >= 15) THEN
             ib = year_/400 - year_/100
          END IF
       END IF
    END IF
    k1 = 365.25_bp*year_
    k2 = 30.6001_bp*(month_+1)
    getMJD_cd = k1 + k2 + ib - 679004 + day

  END FUNCTION getMJD_cd





  !! *Description*:
  !!
  !! Reallocates a pointer array of Time-objects and
  !! copies the data from the old array to the new array (if it fits).
  !!
  !! *Usage*:
  !!
  !! mytimes => reallocate(mytimes,4)
  !!
  !! Returns error.
  !!
  FUNCTION reallocate_T_1(array,n)

    IMPLICIT NONE
    TYPE (Time), DIMENSION(:), POINTER :: reallocate_T_1, array
    INTEGER, INTENT(in)                 :: n
    INTEGER                             :: i, nold, err

    ALLOCATE(reallocate_T_1(n), stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Time / reallocate", &
            "Could not allocate memory.", 1)
       reallocate_T_1 => NULL()
       RETURN
    END IF
    IF (.NOT. ASSOCIATED(array)) RETURN
    nold = SIZE(array)
    DO i=1, MIN(n,nold)
       reallocate_T_1(i) = copy(array(i))
    END DO
    DEALLOCATE(array, stat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("Time / reallocate", &
            "Could not deallocate memory.", 1)
       RETURN
    END IF

  END FUNCTION reallocate_T_1





  !! Description:
  !!
  !! Conversion of a modified Julian date from one timescale 
  !! to another.
  !!
  !! INPUT:    mjd1      -  Modified Julian Day
  !!           scale1    -  Input time scale
  !!           scale2    -  Output (required) time scale
  !!
  !! OUTPUT:   mjd2      -  Modified Julian Day
  !!
  !! Supported time scales:
  !!        UT1
  !!        TAI
  !!        UTC
  !!        TDT = TT = ET
  !!
  !! Returns error.
  !!
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  SUBROUTINE timescaleConversion(mjd1, scale1, mjd2, scale2)

    IMPLICIT NONE
    REAL(bp), INTENT(in)         :: mjd1
    CHARACTER(len=*), INTENT(in) :: scale1, scale2
    REAL(bp), INTENT(out)        :: mjd2
    CHARACTER(len=3)             :: eqsc, eqsc2, scale
    REAL(bp)                     :: mjd2_, dt, mjdt, sec2r, diff, &
         sect, dat, sec1, sec2
    INTEGER                      :: loops, nit, mjd2_real, mjdt_int, mjd2_int

    sec1  = 86400.0_bp*(mjd1 - INT(mjd1))
    mjd2_int = INT(mjd1)
    sec2 = sec1
    ! Current timescale (in which mjd2_int,sec2 are given)
    scale = scale1
    eqsc  = scale1
    eqsc2 = scale2

    loops = 0
    DO
       CALL toNormalForm(mjd2_int, sec2, eqsc)
       IF (error) THEN
          CALL errorMessage("Time / timescaleConversion", &
               "TRACE BACK (5)", 1)
          RETURN
       END IF
       ! Required timescale has been reached
       IF (eqsc == eqsc2) THEN
          mjd2 = mjd2_int + sec2/86400.0_bp
          RETURN
       END IF
       ! Check on infinite loops
       IF (loops > 6) THEN
          error = .TRUE.
          CALL errorMessage("Time / timescaleConversion", &
               "Too many loops.", 1)
          RETURN
       END IF
       loops = loops + 1
       ! Transformations are performed according to the following path:
       !
       !                UT1 -- TDT -- TAI -- UTC
       !
       IF(eqsc == "UT1") THEN
          ! Conversion UT1 --> TDT
          mjd2_ = mjd2_int + sec2/86400.0_bp
          dt = deltaT(mjd2_)
          sec2 = sec2 + dt
          eqsc = "TDT"
       ELSE IF (eqsc == "TDT") THEN
          IF (eqsc2 == "UT1") THEN
             ! Conversion TDT --> UT1 (iterative method)
             !   a) computation of DT = TDT - UT1 using (mjd2_int,sec2) (TDT) as an
             !      approximate value of UT1
             mjd2_ = mjd2_int + sec2/86400.0_bp
             dt = deltaT(mjd2_)
             !   b) subtract DT from (mjd2_int,sec2), finding a first approximation
             !      for UT1
             mjdt_int = mjd2_int
             sect = sec2 - dt
             !   start iterations
             nit = 0
             DO
                CALL toNormalForm(mjdt_int, sect, "UT1")
                IF (error) THEN
                   CALL errorMessage("Time / timescaleConversion", &
                        "TRACE BACK (10)", 1)
                   RETURN
                END IF
                nit = nit + 1
                IF (nit > nitmax) THEN
                   error = .TRUE.
                   CALL errorMessage("Time / timescaleConversion", &
                        "Abnormal end.",1)
                   RETURN
                ENDIF
                !   c) try to find the starting value of TDT from the approximate
                !      value of UT1
                mjdt = mjdt_int + sect/86400.0_bp
                dt = deltaT(mjdt)
                mjd2_real = mjdt_int
                sec2r = sect + dt
                CALL toNormalForm(mjd2_real, sec2r, "TDT")
                IF (error) THEN
                   CALL errorMessage("Time / timescaleConversion", &
                        "TRACE BACK (15)", 1)
                   RETURN
                END IF
                !   d) computation of error and correction of the approximate value
                diff = (mjd2_real-mjd2_int)*86400.0_bp + sec2r - sec2
                IF(ABS(diff) > epst) THEN
                   sect = sect - diff
                ELSE
                   EXIT
                END IF
             END DO
             mjd2_int = mjdt_int
             sec2 = sect
             eqsc = "UT1"
          ELSE
             ! Conversion TDT --> TAI
             sec2 = sec2 - etmtai
             eqsc = "TAI"
          END IF
       ELSE IF (eqsc == "TAI") THEN
          IF (eqsc2 == "UTC") THEN
             ! Conversion TAI --> UTC (iterative method)
             !   a) computation of DAT = TAI - UTC using (mjd2_int,sec2) (TAI) as an
             !      approximate value of UTC
             mjdt_int = mjd2_int
             sect = sec2
             dat = deltaAT(mjdt_int)
             IF (error) THEN
                CALL errorMessage("Time / timescaleConversion", &
                     "TRACE BACK (20)", 1)
                RETURN
             END IF
             !   b) subtract DAT from (mjd2_int,sec2), finding a first approximation
             !      for UTC
             sect = sec2 - dat
             !   start iterations
             nit = 0
             DO
                CALL toNormalForm(mjdt_int, sect, "UTC")
                IF (error) THEN
                   CALL errorMessage("Time / timescaleConversion", &
                        "TRACE BACK (25)", 1)
                   RETURN
                END IF
                nit = nit + 1
                IF (nit > nitmax) THEN
                   error = .TRUE.
                   CALL errorMessage("Time / timescaleConversion", &
                        "Abnormal end.",1)
                   RETURN
                END IF
                !   c) try to find the starting value of TAI from the approximate
                !      value of UTC
                mjd2_real = mjdt_int
                sec2r = sect + deltaAT(mjdt_int)
                IF (error) THEN
                   CALL errorMessage("Time / timescaleConversion", &
                        "TRACE BACK (30)", 1)
                   RETURN
                END IF
                CALL toNormalForm(mjd2_real, sec2r, "TAI")
                IF (error) THEN
                   CALL errorMessage("Time / timescaleConversion", &
                        "TRACE BACK (35)", 1)
                   RETURN
                END IF
                !   d) computation of error and correction of the approximate value
                diff = (mjd2_real-mjd2_int)*86400.0_bp + sec2r - sec2 + &
                     deltaAT(mjd2_real) - deltaAT(mjd2_int)
                IF (error) THEN
                   CALL errorMessage("Time / timescaleConversion", &
                        "TRACE BACK (40)", 1)
                   RETURN
                END IF
                IF (ABS(diff) > epst) THEN
                   sect = sect - diff
                ELSE
                   EXIT   
                END IF
             END DO
             mjd2_int = mjdt_int
             sec2 = sect
             eqsc = "UTC"
          ELSE
             ! Conversion TAI --> TDT
             sec2 = sec2 + etmtai
             eqsc = "TDT"
          END IF
       ELSE IF (eqsc == "UTC") THEN
          ! Conversione UTC --> TAI
          sec2 = sec2 + deltaAT(mjd2_int)
          IF (error) THEN
             CALL errorMessage("Time / timescaleConversion", &
                  "TRACE BACK (45)", 1)
             RETURN
          END IF
          eqsc = "TAI"
       ELSE
          error = .TRUE.
          CALL errorMessage("Time / timescaleConversion", &
               "Abnormal end.",1)
          RETURN
       END IF

    END DO

  END SUBROUTINE timescaleConversion





  !! Description:
  !!
  !! Reduction of time to normal form.
  !!
  !! INPUT:    mjd       -  Modified Julian Day (integer part)
  !!           sec       -  Seconds within day
  !!           scale     -  Time scale
  !!
  !! OUTPUT:   mjd and sec are normalized, namely sec is reduced within
  !!           the limits 0 <= sec < lod, and mjd is changed accordingly;
  !!           lod (length of the day) is usually 86400 s, but can be
  !!           different from that value in the case of UTC, due to
  !!           leap seconds.
  !!
  !! Returns error.
  !! 
  !! Acknowledgements: Based on Mario Carpino's f77 routine.
  !!
  SUBROUTINE toNormalForm(mjd, sec, scale)

    IMPLICIT NONE
    INTEGER, INTENT(inout)       :: mjd
    REAL(bp), INTENT(inout)      :: sec
    CHARACTER(len=*), INTENT(in) :: scale
    INTEGER                      :: nit, idur, isec, k
    REAL(bp)                     :: fsec

    IF (scale == "UTC") THEN
       ! Non-trivial case: UTC (the duration of the day can be different
       ! from 86400 s)
       nit=0
       DO
          nit = nit + 1
          IF (nit > nitutc) THEN
             error = .TRUE.
             CALL errorMessage("Time / toNormalForm", &
                  "Abnormal end.", 1)
             RETURN
          END IF
          IF (sec < 0.0_bp) THEN
             ! Duration in seconds of the previous day
             idur = 86400 + deltaAT(mjd) - deltaAT(mjd-1)
             IF (error) THEN
                CALL errorMessage("Time / toNormalForm", &
                     "TRACE BACK (5)", 1)
                RETURN
             END IF
             sec = sec + idur
             mjd = mjd - 1
          ELSE
             EXIT  
          END IF
       END DO
       ! Decomposition of SEC into integer part (ISEC) + fraction (FSEC),
       ! where 0 <= FSEC < 1
       isec = sec
       fsec = sec - isec
       ! Duration in seconds of today (MJD)
       DO
          idur = 86400 + deltaAT(mjd+1) - deltaAT(mjd)
          IF (error) THEN
             CALL errorMessage("Time / toNormalForm", &
                  "TRACE BACK (10)", 1)
             RETURN
          END IF
          ! Renormalization of time
          IF (isec >= idur) THEN
             isec = isec - idur
             mjd = mjd + 1
          ELSE
             EXIT
          END IF
       END DO
       sec = isec + fsec

    ELSE

       ! Trivial case: the duration of the day is always 86400 s
       ! Also this case requires iterations, due to rounding-off problems.
       ! EXAMPLE: Let us suppose that the starting values are MJD=48000,
       ! SEC=-1.d-14. The result of the first iteration is then:
       !    SEC --> SEC+86400 = 86400 EXACTLY (due to rounding off)
       !    MJD --> MJD-1 = 47999
       ! Therefore, a second iteration is required, giving:
       !    SEC --> SEC-86400 = 0
       !    MJD --> MJD+1 = 48000
       nit = 0
       DO
          nit = nit + 1
          IF (nit > nitmax) THEN
             error = .TRUE.
             CALL errorMessage("Time / toNormalForm", &
                  "Abnormal end.", 1)
             RETURN
          END IF
          k = sec/86400.0_bp
          IF (sec < 0.0_bp) k = k - 1
          IF (k == 0) RETURN
          mjd = mjd + k
          sec = sec - k*86400.0_bp
       END DO
    END IF

  END SUBROUTINE toNormalForm





END MODULE Time_cl





