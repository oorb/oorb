!====================================================================!
!                                                                    !
! Copyright 2002-2015,2016                                           !
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
!! This is actually not a class in a strict sense, but can anyway be
!! thought of as a sort of a class that all other classes should
!! inherit. It contains the most fundamental parameters, types and
!! routines, that should be available for all classes.
!!
!! IAU units are used (if not otherwise stated):
!! - mass = solar mass
!! - distance = astronomical unit (AU)
!! - time = day
!! - angle = rad
!!
!! @author  MG, JV, TL
!! @version 2016-03-23
!!
MODULE Base_cl

  USE utilities
  USE parameters
  IMPLICIT NONE
  PRIVATE :: calendarDateToJulianDate
  PRIVATE :: coordinatedUniversalTime
  PRIVATE :: FNAME_LEN             !! defined in module parameters

  INTEGER, PARAMETER :: stderr = 0 !! Standard error logical unit. 
  !!                                  Shouldn't be connected to a file.
  INTEGER, PARAMETER :: stdin  = 5 !! Standard input logical unit 
  !!                                  (refers usually to the keyboard).
  INTEGER, PARAMETER :: stdout = 6 !! Standard output logical unit 
  !!                                  (refers usually to the screen).
  INTEGER, PARAMETER :: debug  = 7 !! Debugging logical unit. 
  INTEGER, PARAMETER :: cbp =  SELECTED_REAL_KIND(p=12) !! Base precision kind for type complex.
  INTEGER, PARAMETER :: hp  =  SELECTED_REAL_KIND(p=18) !! High precision kind for type real.
  INTEGER, PARAMETER :: bp  =  SELECTED_REAL_KIND(p=12) !! Base precision kind for type real.
  INTEGER, PARAMETER :: lp  =  SELECTED_REAL_KIND(p=12)  !! Low precision kind for type real. p=6
  INTEGER, PARAMETER :: ihp =  SELECTED_INT_KIND(12)    !! High precision kind for type integer.
  INTEGER, PARAMETER :: ibp =  SELECTED_INT_KIND(7)     !! Base precision kind for type integer.
  INTEGER, PARAMETER :: ilp =  SELECTED_INT_KIND(4)     !! Low precision kind for type integer.
  REAL(bp) :: timezone = 10.0_bp  !! Difference between UT and local time [h]. Default: 10 (HST).
  CHARACTER(len=3), DIMENSION(12), PARAMETER :: month_abbr = (/ &
       "Jan", &
       "Feb", &
       "Mar", &
       "Apr", &
       "May", &
       "Jun", &
       "Jul", &
       "Aug", &
       "Sep", &
       "Oct", &
       "Nov", &
       "Dec"   /)
  ! Lengths of months in days:
  INTEGER, DIMENSION(12), PARAMETER :: monthlen = (/ 31, 28, 31, &
       30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  REAL(bp), PARAMETER :: pi       = 3.14159265358979324_bp  !! Pi.
  REAL(bp), PARAMETER :: two_pi   = 2.0_bp*pi               !! 2*Pi.
  REAL(bp), PARAMETER :: rad_hour = pi/12.0_bp              !! Radians per hour.
  REAL(bp), PARAMETER :: rad_min  = rad_hour/60.0_bp        !! Radians per minute.
  REAL(bp), PARAMETER :: rad_sec  = rad_min/60.0_bp         !! Radians per second.
  REAL(bp), PARAMETER :: rad_deg  = pi/180.0_bp             !! Radians per degree.
  REAL(bp), PARAMETER :: rad_amin = rad_deg/60.0_bp         !! Radians per arcminute.
  REAL(bp), PARAMETER :: rad_asec = rad_amin/60.0_bp        !! Radians per arcsecond.
  REAL(bp), PARAMETER :: km_au    = 1.495978707e8_bp       !! Kilometres per astronomical unit (IAU 2012, HORIZONS).
  REAL(bp), PARAMETER :: m_au     = 1.495978707e11_bp      !! Metres per astronomical unit (IAU 2012, HORIZONS).
  REAL(bp), PARAMETER :: day_year = 365.256363004_bp        !! Days per sidereal year (J2000.0).
  REAL(bp), PARAMETER :: hour_day = 24.0_bp                 !! Hours per day.
  REAL(bp), PARAMETER :: min_day  = hour_day*60.0_bp        !! Minutes per day.
  REAL(bp), PARAMETER :: sec_day  = min_day*60.0_bp         !! Seconds per day.
  REAL(bp), PARAMETER :: sec_year = sec_day*day_year        !! Seconds per year.
  REAL(bp), PARAMETER :: kg_solar = 1.9892e30_bp            !! Kilograms per solar mass (IAU 1976).  
  REAL(bp), PARAMETER :: kB = 1.38e-23_bp                   !! Boltzmans constant kB [J K^-1]
  REAL(bp), PARAMETER :: hPl = 6.63e-34_bp                  !! Planck's constant h [J s]
  REAL(bp), PARAMETER :: sol      = 299792.458_bp/km_au*sec_day !! Speed of light (AU/d).
  !REAL(bp), PARAMETER :: eps = 23.439280833_bp*rad_deg     !! Obliquity of ecliptic (J2000.0) 
  !!                                                           Gaia 2006 (84381.41100 arcsec).
  REAL(bp), PARAMETER :: eps = 23.43929111111111_bp*rad_deg !! Obliquity of ecliptic (J2000.0;
  !!                                                           84381.448 arcsec), IAU 1976
  !!                                                           (used by HORIZONS, MPC, OrbFit)
  REAL(bp), PARAMETER :: r_earth  = 6378.140_bp/km_au       !! Earth equatorial radius (not mean!).
  REAL(bp), PARAMETER :: smamax   = 500.0_bp                !! Maximum for semimajor axis (AU).
  REAL(bp), PARAMETER :: rmoon    = 2.57e-3_bp              !! Earth-Moon mean distance (AU).
  REAL(bp), PARAMETER :: kgm3_smau3 = (km_au)**3/kg_solar   !! density conversion
  INTEGER, PARAMETER :: DIR_LEN = 256
  INTEGER, PARAMETER :: OBS_RECORD_LEN = 256
  INTEGER, PARAMETER :: ELEMENT_TYPE_LEN = 16
  INTEGER, PARAMETER :: FRAME_LEN = 16
  INTEGER, PARAMETER :: DYN_MODEL_LEN = 16
  INTEGER, PARAMETER :: INTEGRATOR_LEN = 32
  INTEGER, PARAMETER :: DESIGNATION_LEN = 16
  INTEGER, PARAMETER :: OBSY_CODE_LEN = 4
  CHARACTER(len=FNAME_LEN) :: OORB_DATA_DIR
  INCLUDE "../prefix.h"

  CHARACTER(len=FNAME_LEN), PARAMETER :: EPH_FNAME = "de430.dat"
  ! OBS CODES
  CHARACTER(len=FNAME_LEN), PARAMETER :: CODE_FNAME = "OBSCODE.dat"
  ! Default name of the file containing the difference between ET and UT:
  CHARACTER(len=FNAME_LEN), PARAMETER :: ETUT_FNAME   = "ET-UT.dat"
  ! Default name of the file containing the difference between TAI and UTC:
  CHARACTER(len=FNAME_LEN), PARAMETER :: TAIUTC_FNAME = "TAI-UTC.dat"
  ! ET-TAI in seconds:
  REAL(bp), PARAMETER :: ETMTAI = 32.184_bp

  ! MPC conversion table
  CHARACTER(len=1), DIMENSION(0:61), PARAMETER :: mpc_conv_table = &
       (/ "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", &
       "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", &
       "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", &
       "Z", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", &
       "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", &
       "x", "y", "z" /)
  ! MPC conversion string
  CHARACTER(len=62), PARAMETER :: mpc_conv_str = &
       "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

  ! Verbosity of error messages
  INTEGER :: err_verb = 1
  ! Verbosity of information messages
  INTEGER :: info_verb = 1 
  ! Number of massless particles integrated simultaneously
  INTEGER :: simint = 1
  ! 
  LOGICAL :: simulated_observations = .FALSE.
  ! Character string error flag
  CHARACTER(len=1024) :: errstr = ""
  ! Logical error flag
  LOGICAL :: error = .FALSE.

  TYPE :: Vector
     REAL(bp), DIMENSION(:,:), POINTER :: elements
  END TYPE Vector

  TYPE :: SparseArray
     TYPE(Vector), DIMENSION(:), POINTER :: vectors
  END TYPE SparseArray

CONTAINS




  !! *Description*:
  !!
  !! Given cosinus and sinus of an angle, returns the angle in radians.
  !!
  REAL(bp) FUNCTION angle(cos_angle, sin_angle)

    IMPLICIT NONE
    REAL(bp) :: cos_angle, sin_angle

    ! Safety check (problems with accuracy)
    IF (ABS(cos_angle) > 1.0_bp) cos_angle = SIGN(1.0_bp,cos_angle)*1.0_bp
    IF (ABS(sin_angle) > 1.0_bp) sin_angle = SIGN(1.0_bp,sin_angle)*1.0_bp

    IF (SIGN(1.0_bp,cos_angle) == SIGN(1.0_bp,sin_angle)) THEN
       IF (cos_angle >= 0.0_bp .OR. sin_angle == 0.0_bp) THEN
          angle = ACOS(cos_angle)
       ELSE
          angle = two_pi - ACOS(cos_angle)
       ENDIF
    ELSE IF (SIGN(1.0_bp,cos_angle) /= SIGN(1.0_bp,sin_angle)) THEN
       IF (SIGN(1.0_bp,cos_angle) < 0.0_bp) THEN
          angle = ACOS(cos_angle)
       ELSE
          angle = two_pi-ACOS(cos_angle)
       END IF
    END IF

  END FUNCTION angle





  !! Description:
  !!
  !! This routine converts calendar date and time to the corresponding
  !! Julian date. The algorithm from Van Flandern and Pulkkinen,
  !! Ap.J. Suppl. vol 41, page 392.
  !!
  SUBROUTINE calendarDateToJulianDate(year, month, day, h, min, sec, jd)

    IMPLICIT NONE
    INTEGER, INTENT(in)   :: year, month, day, h, min
    REAL(bp), INTENT(in)  :: sec
    REAL(bp), INTENT(out) :: jd
    INTEGER               :: iaux

    iaux = -7*(year+(month+9)/12)/4 - &
         3*((year+(month-9)/7)/100+1)/4 + 275*month/9
    jd = iaux + day + 367.0_bp*year + 1721028.5_bp + &
         (h + min/60.0_bp + sec/3600.0_bp)/24.0_bp

  END SUBROUTINE calendarDateToJulianDate





  !! *Description*:
  !!
  !! Transforms degrees, arcminutes, and/or -seconds to radians.
  !!
  SUBROUTINE DAMAStoRadians(deg, arcmin, arcsec, rad)

    IMPLICIT NONE
    INTEGER, INTENT(in)   :: deg, arcmin
    REAL(bp), INTENT(in)  :: arcsec
    REAL(bp), INTENT(out) :: rad

    rad = (deg + arcmin/60.0_bp + arcsec/3600.0_bp)*pi/180.0_bp

  END SUBROUTINE DAMAStoRadians





  !! Description:
  !!
  !! Converts fractional days to hours, minutes, and seconds. Note
  !! that "d" must be positive!
  !!
  SUBROUTINE dayToHMS(d, h, m, s)

    IMPLICIT NONE
    REAL(bp), INTENT(in)  :: d
    REAL(bp), INTENT(out) :: s
    INTEGER, INTENT(out)  :: h, m
    REAL(bp)              :: fh, fm

    fh = 24.0_bp*d
    fh = fh + EPSILON(fh)
    h = INT(fh)
    fm = 60.0_bp*(fh - REAL(h, bp))
    fm = fm + EPSILON(fm)
    m = INT(fm)
    s = 60.0_bp*(fm - REAL(m, bp))

  END SUBROUTINE dayToHMS





  !! *Description*:
  !! 
  !! This routines reports errors in a standard fashion. The first
  !! character string is the name of the routine from which this
  !! routine is invoked. It is nice to know from where this message
  !! originates. The other character string is the error
  !! description. The error message is written to the standard error
  !! logical unit, and can be presented to the end-user.
  !!
  !! The standard array of I/O errors can occur but they are ignored
  !! because this is an error handling routine and it would be easy
  !! to get into an error loop.
  !!
  !! *Usage*:
  !!
  !! <pre>
  !! use Base_cl !(contains global logical variable "error")
  !! .
  !! .
  !! .
  !! subroutine thisroutine(...)
  !! .
  !! .
  !! .
  !! if (error) then
  !! call errorMessage("thisroutine","Experienced this error.",1)
  !! return
  !! end if
  !! </pre>
  !!
  SUBROUTINE errorMessage(routine_name, msg_str, vrbs)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)    :: routine_name, msg_str
    INTEGER, INTENT(in)             :: vrbs
    CHARACTER(len=3), DIMENSION(12) :: month_abbrev
    CHARACTER(len=256)              :: form
    INTEGER, DIMENSION(3)           :: date, time
    INTEGER                         :: err

    IF (err_verb >= vrbs) THEN
       month_abbrev = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", &
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)
       CALL coordinatedUniversalTime(date,time)
       form = "('***ERROR*** ', i2, ' ', a3, ' ', i4, ' ', i2.2, " // &
            "':', i2.2, ':', i2.2, 'UTC (', a, ') ', a)"
       WRITE(stderr, TRIM(form), iostat=err) date(1), month_abbrev(date(2)), &
            date(3), time(1), time(2), time(3), TRIM(routine_name), &
            TRIM(msg_str)
       IF (err /= 0) WRITE(stderr,*) "Could not write error message!" 
    END IF

  END SUBROUTINE errorMessage





  !! *Description*:
  !! 
  !! This routines reports warnings in a standard fashion. The first
  !! character string is the name of the routine from which this
  !! routine is invoked. It is nice to know from where this message
  !! originates. The other character string is the warning
  !! description. The warning message is written to the standard error
  !! logical unit, and can be presented to the end-user.
  !!
  !! The standard array of I/O errors can occur but they are ignored
  !! because this is an error handling routine and it would be easy
  !! to get into an error loop.
  !!
  !! *Usage*:
  !!
  !! <pre>
  !! use Base_cl !(contains global logical variable "error")
  !! .
  !! .
  !! .
  !! subroutine thisroutine(...)
  !! .
  !! .
  !! .
  !! if (error) then
  !! call warningMessage("thisroutine","Experienced this event.",1)
  !! return
  !! end if
  !! </pre>
  !!
  SUBROUTINE warningMessage(routine_name, msg_str, vrbs)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)    :: routine_name, msg_str
    INTEGER, INTENT(in)             :: vrbs
    CHARACTER(len=3), DIMENSION(12) :: month_abbrev
    CHARACTER(len=256)              :: form
    INTEGER, DIMENSION(3)           :: date, time
    INTEGER                         :: err

    IF (err_verb >= vrbs) THEN
       month_abbrev = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", &
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)
       CALL coordinatedUniversalTime(date,time)
       form = "('***WARNING*** ', i2, ' ', a3, ' ', i4, ' ', i2.2, " // &
            "':', i2.2, ':', i2.2, 'UTC (', a, ') ', a)"
       WRITE(stderr, TRIM(form), iostat=err) date(1), month_abbrev(date(2)), &
            date(3), time(1), time(2), time(3), TRIM(routine_name), &
            TRIM(msg_str)
       IF (err /= 0) WRITE(stderr,*) "Could not write warning message!" 
    END IF

  END SUBROUTINE warningMessage





  !! *Description*:
  !!
  !! Transforms hours, minutes, and/or seconds to radians.
  !!
  SUBROUTINE HMStoRadians(hour, min, sec, rad)

    IMPLICIT NONE
    INTEGER, INTENT(in)   :: hour, min
    REAL(bp), INTENT(in)  :: sec
    REAL(bp), INTENT(out) :: rad

    rad = (hour + min/60.0_bp + sec/3600.0_bp)*two_pi/24.0_bp

  END SUBROUTINE HMStoRadians





  !! *Description*:
  !! 
  !! This routines prints information in a standard fashion. The first
  !! character string is the name of the routine from which this
  !! routine is invoked. It is nice to know from where this message
  !! originates. The other character string is the message. The info
  !! message is written to the standard error logical unit, and can be
  !! presented to the end-user.
  !!
  !! *Usage*:
  !!
  !! <pre>
  !! use Base_cl
  !! .
  !! .
  !! .
  !! call infoMessage("thisroutine", "Experienced this thing.", stdout, 1)
  !! return
  !! end if
  !! </pre>
  !!
  SUBROUTINE infoMessage(routine_name, msg_str, lu, vrbs)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in)    :: routine_name, msg_str
    INTEGER, INTENT(in)             :: lu, vrbs
    CHARACTER(len=3), DIMENSION(12) :: month_abbrev
    CHARACTER(len=256)              :: form
    INTEGER, DIMENSION(3)           :: date, time
    INTEGER                         :: err

    IF (info_verb >= vrbs) THEN
       month_abbrev = (/ "Jan", "Feb", "Mar", "Apr", "May", "Jun", &
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec" /)
       CALL coordinatedUniversalTime(date,time)
       form = "('***INFO*** ', i2, ' ', a3, ' ', i4, ' ', i2.2, " // &
            "':', i2.2, ':', i2.2, 'UTC (', a, ') ', a)"
       WRITE(lu,TRIM(form),iostat=err) date(1), month_abbrev(date(2)), &
            date(3), time(1), time(2), time(3), TRIM(routine_name), &
            TRIM(msg_str)
       IF (err /= 0) WRITE(stderr,*) "Could not write message!" 
    END IF

  END SUBROUTINE infoMessage





  !! *Description*:
  !!
  !! This routine converts the Julian date to calendar date and time.
  !! The algorithm from Q.J.R. Astr. Soc. 1984, vol 25) page 53.
  !!
  !! Note:
  !!
  !! For 100 year dates which are not leap years
  !! (ie 1900, 2100, etc) the conversion at mar 1 is wrong.
  !!       eg    2415078.5  --->  1900,2,28  (correct)
  !!             2415079.5  --->  1900,2,29  (should be 1900,3,1 !)
  !!             2415080.5  --->  1900,3,2   (correct)
  !!
  SUBROUTINE julianDateToCalendarDate(jd, year, month, day, hour, min, sec)

    IMPLICIT NONE
    REAL(bp), INTENT(in)  :: jd
    REAL(bp), INTENT(out) :: sec
    INTEGER, INTENT(out)  :: year, month, day, hour, min
    REAL(bp)              :: fjd, jul, fday
    INTEGER               :: ijul, ig, idp

    jul = jd
    fjd = jul - INT(jul)
    IF(fjd >= 0.5_bp) THEN
       ijul = jul + 1.0_bp
       fday = fjd - 0.5_bp
    ELSE
       ijul = jul
       fday = fjd + 0.5_bp
    ENDIF
    ig = INT(INT((jul-4479.5_bp)/36524.25)*0.75_bp + 0.5_bp) - 37
    jul = ijul + ig
    year = INT(jul/365.25_bp) - 4712
    idp = INT(MOD(jul-59.25_bp, 365.25_bp))
    month = 1 + MOD(INT((idp+0.5_bp)/30.6_bp)+2,12)
    day = INT(MOD(idp+0.5_bp, 30.6_bp)) + 1
    CALL dayToHMS(fday,hour,min,sec)

  END SUBROUTINE julianDateToCalendarDate





  !! Mean obliquity of ecliptic (J2000.0?). See Astronomical Almanac
  !! 1987, B18.
  !!
  !!  INPUT:    mjd   -   Modified Julian Date (TT)
  !!
  !!  OUTPUT:   in radians
  !!
  REAL(bp) FUNCTION meanObliquity(mjd_tt)

    IMPLICIT NONE

    REAL(bp), INTENT(in) :: mjd_tt
    REAL(bp)             :: t0
    REAL(bp), PARAMETER  :: ob0 = &
         (23.0_bp*3600.0_bp+26.0_bp*60.0_bp+21.448_bp)*rad_asec, & ! IAU value - 0.02"
         ob1 = -46.8150_bp  * rad_asec, & ! -46.815_bp
         ob2 = -0.00059_bp  * rad_asec, & ! -0.0006_bp
         ob3 =  0.001813_bp * rad_asec    ! 0.00181_bp

    t0     = ( mjd_tt - 51544.5_bp ) / 36525.0_bp
    meanObliquity = (( ob3 * t0 + ob2 ) * t0 + ob1 ) * t0 + ob0

  END FUNCTION meanObliquity





  !! *Description*:
  !!
  !! Transforms radians to degrees, arcminutes, and -seconds.
  !!
  SUBROUTINE radiansToDAMAS(rad, deg, arcmin, arcsec)

    IMPLICIT NONE
    REAL(bp), INTENT(in)  :: rad
    INTEGER, INTENT(out)  :: deg, arcmin
    REAL(bp), INTENT(out) :: arcsec

    deg    = INT(ABS(rad)/rad_deg)
    arcmin = INT(60.0_bp*(ABS(rad)/rad_deg-deg))
    arcsec = 3600.0_bp*(ABS(rad)/rad_deg-(deg+arcmin/60.0_bp))
    deg    = INT(rad/rad_deg)

  END SUBROUTINE radiansToDAMAS





  !! *Description*:
  !!
  !! Transforms radians to hours, minutes, and -seconds.
  !!
  SUBROUTINE radiansToHMS(rad, hour, min, sec)

    IMPLICIT NONE
    REAL(bp), INTENT(in)  :: rad
    INTEGER, INTENT(out)  :: hour, min
    REAL(bp), INTENT(out) :: sec

    hour = INT(rad/(rad_deg*15.0_bp))
    min  = INT(60.0_bp*(rad/(rad_deg*15.0_bp)-hour))
    sec  = 3600.0_bp*(rad/(rad_deg*15.0_bp)-(hour+min/60.0_bp))

  END SUBROUTINE radiansToHMS





  !! *Description*:
  !!
  !! Rotation matrix around i'th axis: If X are "old" coordinates and
  !! X' are "new" coordinates (referred to a frame which is rotated by
  !! an angle alpha around an axis in direct sense), then X' = R X,
  !! where R is the rotation matrix.
  !!
  !! Returns error.
  !!
  FUNCTION rotationMatrix(alpha, axis)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: alpha
    INTEGER, INTENT(in) :: axis
    REAL(bp), DIMENSION(3,3) :: rotationMatrix 
    INTEGER :: i1, i2, i3

    IF (axis < 1 .OR. axis > 3) THEN
       error = .TRUE.
       CALL errorMessage("Base / rotationMatrix", &
            "Axis out of range. It has to be 1, 2, or 3.", 1)
       RETURN
    END IF

    i1 = axis
    IF(i1 > 3) THEN
       i1 = i1 - 3
    END IF
    i2 = i1 + 1
    IF(i2 > 3) THEN
       i2 = i2 - 3
    END IF
    i3 = i2 + 1
    IF(i3 > 3) THEN
       i3 = i3 - 3
    END IF

    rotationMatrix(i1,i1) = 1.0_bp
    rotationMatrix(i1,i2) = 0.0_bp
    rotationMatrix(i1,i3) = 0.0_bp

    rotationMatrix(i2,i1) = 0.0_bp
    rotationMatrix(i2,i2) = COS(alpha)
    rotationMatrix(i2,i3) = SIN(alpha)

    rotationMatrix(i3,i1) = 0.0_bp
    rotationMatrix(i3,i2) = -SIN(alpha)
    rotationMatrix(i3,i3) = COS(alpha)

  END FUNCTION rotationMatrix


  FUNCTION resolveDirectory(subdir, envvar) RESULT(s2)
    CHARACTER(*), INTENT(IN) :: subdir, envvar
    CHARACTER(FNAME_LEN) :: s2

    ! If overriden by an environmental variable, prefer that
    CALL getenv(envvar, s2)
    IF (LEN_TRIM(s2) /= 0) THEN
       RETURN
    END IF

    ! Otherwise, return <PREFIX>/<subdir>
    s2 = TRIM(PREFIX) // "/" // subdir

  END FUNCTION resolveDirectory


  SUBROUTINE setAccessToDataFiles()

    IMPLICIT NONE
    OORB_DATA_DIR = resolveDirectory("share/oorb", "OORB_DATA")

  END SUBROUTINE setAccessToDataFiles




  !! *Description*:
  !!
  !! Sets the global error verbose variable to a given value.
  !! The value can be any integer between 1 and 5. By setting
  !! it to 1, only the most critical error messages are shown.
  !! Number 5 corresponds to all error messages.   
  !!
  SUBROUTINE setErrorVerbose(error)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: error

    SELECT CASE (error)
    CASE (:-1)
       err_verb = 0
    CASE (0:5) 
       err_verb = error
    CASE (6:)
       err_verb = 5
    CASE default
       err_verb = 1
    END SELECT

  END SUBROUTINE setErrorVerbose





  !! *Description*:
  !!
  !! Sets the global information verbose variable to a given value.
  !! The value can be any integer between 1 and 5. By setting
  !! it to 1, only the most critical informative messages are shown.
  !! Number 5 corresponds to all informative messages.   
  !!
  SUBROUTINE setInfoVerbose(info)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: info

    SELECT CASE (info)
    CASE (:-1)
       info_verb = 0
    CASE (0:5) 
       info_verb = info
    CASE (6:)
       info_verb = 5
    CASE default
       info_verb = 1
    END SELECT

  END SUBROUTINE setInfoVerbose





  !! *Description*:
  !!
  !! If simulated observations (no light time correction, and
  !! epoch given in TDT) are used, set "simulated" to _true_.
  !! The default value is _false_.
  !!
  SUBROUTINE setSimulatedObservations(simulated)

    IMPLICIT NONE
    LOGICAL, INTENT(in) :: simulated

    simulated_observations = simulated

  END SUBROUTINE setSimulatedObservations





  !! *Description*:
  !!
  !! Returns the angular distance in radians between two points on the
  !! sky (that is, a unit sphere).
  !!
  REAL(bp) FUNCTION angularDistance(ra1, dec1, ra2, dec2)

    IMPLICIT NONE
    REAL(bp), INTENT(in) :: ra1, dec1, ra2, dec2
    REAL(bp) :: cos_distance

    cos_distance = COS(ra2-ra1)*COS(dec1)*COS(dec2) + &
         SIN(dec1)*SIN(dec2)
    IF (ABS(cos_distance) > 1.0_bp) THEN
       cos_distance = SIGN(1.0_bp, cos_distance)
    END IF
    angularDistance = ACOS(cos_distance)

  END FUNCTION angularDistance





  !! *Description*:
  !!
  !! Transforms MPC3 type designation to MPC type.
  !!
  !! Example:
  !!
  !! K04A0001S -> K04A01S
  !!
  !! Returns error.
  !!
  SUBROUTINE MPC3DesToMPCDes(designation)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: designation
    INTEGER :: i

    CALL removeLeadingBlanks(designation)
    CALL toInt(designation(6:7), i, error)
    IF (error) THEN
       CALL errorMessage("Base / MPC3DesToMPCDes", &
            "Could not convert string to integer.", 1)
       RETURN
    END IF
    designation(1:7) = designation(1:4) // mpc_conv_table(i) // &
         designation(8:9)
    designation(8:) = " "

  END SUBROUTINE MPC3DesToMPCDes





  !! *Description*:
  !!
  !! Transforms MPC type designation to MPC3 type.
  !!
  !! Example:
  !!
  !! K04A01S -> K04A0001S
  !!
  !! Returns error.
  !!
  SUBROUTINE MPCDesToMPC3Des(designation)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: designation
    INTEGER :: i

    CALL removeLeadingBlanks(designation)
    IF (.NOT.(designation(1:1) == "I" .OR. &
         designation(1:1) == "J" .OR. &
         designation(1:1) == "K")) THEN
       RETURN
    END IF
    DO i=0,SIZE(mpc_conv_table,dim=1)-1
       IF (mpc_conv_table(i) == designation(5:5)) THEN
          EXIT
       END IF
    END DO
    designation(5:5) = "0"
    designation(8:9) = designation(6:7)
    IF (i <= 9) THEN
       designation(6:6) = "0"
       CALL toString(i, designation(7:7), error)
    ELSE IF (i > 9) THEN
       CALL toString(i, designation(6:7), error)
    END IF

  END SUBROUTINE MPCDesToMPC3Des





  !! *Description*:
  !!
  !! Decodes an encoded MPC type designation.
  !!
  !! Example:
  !!
  !! K04A01S -> 2004 AS1
  !!
  !! Returns error.
  !!
  SUBROUTINE decodeMPCDesignation(designation)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: designation
    CHARACTER(len=DESIGNATION_LEN) :: designation_
    CHARACTER(len=2) :: str
    INTEGER :: i

    CALL removeLeadingBlanks(designation)
    i = IACHAR(designation(1:1))
    IF (INDEX(TRIM(designation)," ") /= 0 .OR. &
         (i >= 48 .AND.  i <= 57)) THEN
       ! Designation already seems to be decoded, because it 
       ! contains empty space or starts with a number.
       CONTINUE
    ELSE IF (LEN_TRIM(designation) == 5) THEN
       i = INDEX(mpc_conv_str, designation(1:1)) - 1
       CALL toString(i, str, error)
       designation = TRIM(str) // designation(2:5)
    ELSE IF (designation(3:3) == "S") THEN
       ! Designation of type T3S1234:
       designation_ = ""
       designation_(1:4) = designation(4:7)
       designation_(5:5) = " "
       designation_(6:8) = designation(1:1) // "-" // &
            designation(2:2)
       designation = " "
       designation = TRIM(designation_)
    ELSE
       ! Designation of type K01A12B:
       designation_ = " "
       DO i=0,SIZE(mpc_conv_table,dim=1)-1
          IF (mpc_conv_table(i) == designation(1:1)) THEN
             EXIT
          END IF
       END DO
       CALL toString(i, designation_(1:2), error)
       IF (error) THEN
          CALL errorMessage("Base / decodeMPCDesignation", &
               "Could not convert integer to string (5).", 1)       
          WRITE(stderr,*) TRIM(designation)
          RETURN
       END IF
       designation_(3:4) = designation(2:3)
       designation_(5:5) = " "
       designation_(6:6) = designation(4:4)
       designation_(7:7) = designation(7:7)
       DO i=0,SIZE(mpc_conv_table,dim=1)-1
          IF (mpc_conv_table(i) == designation(5:5)) THEN
             EXIT
          END IF
       END DO
       CALL toString(i, designation_(8:9), error)
       IF (error) THEN
          CALL errorMessage("Base / decodeMPCDesignation", &
               "Could not convert integer to string (10).", 1)       
          WRITE(stderr,*) TRIM(designation)
          RETURN
       END IF
       IF (designation_(8:9) == "00") THEN
          designation_(8:) = " "
       ELSE IF (designation_(8:8) == "0") THEN
          designation_(8:8) = designation_(9:9)
       END IF
       IF (LEN_TRIM(designation_) >= 8 .OR. designation(6:6) /= "0") THEN
          designation_ = TRIM(designation_) // designation(6:6)
       END IF
       designation = TRIM(designation_)
    END IF
    ! Remove leading zeros
    DO WHILE (designation(1:1) == "0")
       designation_ = ""
       designation_(1:) = designation(2:LEN_TRIM(designation))
       designation = ""
       designation = TRIM(designation_)
    END DO

  END SUBROUTINE decodeMPCDesignation





  !! *Description*:
  !!
  !! Decodes an encoded MPC type designation.
  !!
  !! Example:
  !!
  !! K04A0001S -> 2004 AS1
  !!
  !! Returns error.
  !!
  SUBROUTINE decodeMPC3Designation(designation)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: designation
    CHARACTER(len=DESIGNATION_LEN) :: designation_
    INTEGER :: i

    CALL removeLeadingBlanks(designation)
    i = IACHAR(designation(1:1))
    IF (INDEX(TRIM(designation)," ") /= 0 .OR. &
         (i >= 48 .AND.  i <= 57)) THEN
       ! Designation already seems to be decoded, because it 
       ! contains empty space or starts with a number.
       RETURN
    END IF

    IF (designation(3:3) == "S") THEN
       ! Designation of type T3S1234:
       designation_(1:4) = designation(4:7)
       designation_(5:5) = " "
       designation_(6:8) = designation(1:1) // "-" // &
            designation(2:2)
       designation = " "
       designation = TRIM(designation_)
    ELSE
       ! Designation of type K01A0012B:
       designation_ = " "
       DO i=0,SIZE(mpc_conv_table,dim=1)-1
          IF (mpc_conv_table(i) == designation(1:1)) THEN
             EXIT
          END IF
       END DO
       CALL toString(i, designation_(1:2), error)
       IF (error) THEN
          CALL errorMessage("Base / decodeMPCDesignation", &
               "Could not convert integer to string (5).", 1)       
          WRITE(stderr,*) TRIM(designation)
          RETURN
       END IF
       designation_(3:4) = designation(2:3)
       designation_(5:5) = " "
       designation_(6:6) = designation(4:4)
       designation_(7:7) = designation(9:9)
       i = 5
       DO WHILE (designation(i:i) == "0")
          i = i + 1
          IF (i > 8) THEN
             EXIT
          END IF
       END DO
       IF (i <= 8) THEN
          designation_ = TRIM(designation_) // designation(i:8)
       END IF
       designation = TRIM(designation_)
    END IF

  END SUBROUTINE decodeMPC3Designation





  !! *Description*:
  !!
  !! Encodes a decoded MPC type designation.
  !!
  !! Example:
  !!
  !! 2004 AS1 -> K04A01S 
  !!
  !! Returns error.
  !!
  SUBROUTINE encodeMPCDesignation(designation)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: designation
    CHARACTER(len=DESIGNATION_LEN) :: designation_
    INTEGER :: i

    CALL removeLeadingBlanks(designation)
    i = IACHAR(designation(1:1))
    IF (INDEX(TRIM(designation)," ") == 0 .OR. &
         i < 48 .OR.  i > 57) THEN
       ! Designation already seems to be encoded, because it has no
       ! empty space or starts with a non-number.
       RETURN
    END IF

    IF (INDEX(TRIM(designation),"-") /= 0) THEN
       ! Designation of type 1234 A-B:
       designation_(1:1) = designation(6:6)
       designation_(2:2) = designation(8:8)
       designation_(3:7) = "S" // designation(1:4)
       designation = designation_(1:7)
       designation(8:) = " "
    ELSE
       ! Designation of type 2001 AB12:
       CALL toInt(designation(1:2), i, error)
       IF (error) THEN
          CALL errorMessage("Base / encodeMPCDesignation", &
               "Could not convert string to integer (5).", 1)       
          WRITE(stderr,*) TRIM(designation)
          RETURN
       END IF
       designation_ = " "
       designation_(1:3) = mpc_conv_table(i) // designation(3:4)
       designation_(4:4) = designation(6:6)
       IF (LEN_TRIM(designation) == 7) THEN
          designation_(5:6) = "00" // designation(8:8)
       ELSE IF (LEN_TRIM(designation) == 8) THEN
          designation_(5:6) = "0" // designation(8:8)
       ELSE IF (LEN_TRIM(designation) == 9) THEN
          designation_(5:6) = designation(8:9)       
       ELSE
          CALL toInt(designation(8:9), i, error)
          IF (error) THEN
             CALL errorMessage("Base / encodeMPCDesignation", &
                  "Could not convert string to integer (10).", 1)       
             WRITE(stderr,*) TRIM(designation)
             RETURN
          END IF
          designation_(5:6) = mpc_conv_table(i) // designation(10:10)
       END IF
       designation_(7:7) = designation(7:7)
       designation = " "
       designation = TRIM(designation_)
    END IF

  END SUBROUTINE encodeMPCDesignation





  !! *Description*:
  !!
  !! Encodes a number using the MPC scheme.
  !!
  !! Example:
  !!
  !! 100001 -> A0001 
  !!
  !! Returns error.
  !!
  SUBROUTINE encodeMPCNumber(number, length)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: number
    INTEGER, INTENT(in) :: length
    CHARACTER(len=DESIGNATION_LEN) :: number_
    INTEGER :: i

    CALL removeLeadingBlanks(number)
    i = IACHAR(number(1:1))
    IF (i < 48 .OR.  i > 57 .OR. LEN_TRIM(number) <= length) THEN
       ! Number either doesn't need to be coded or already seems to be
       ! encoded, because it starts with a non-number.
       RETURN
    END IF

    ! Number of type 100001:
    CALL toInt(number(1:2), i, error)
    IF (error) THEN
       CALL errorMessage("Base / encodeMPCNumber", &
            "Could not convert string to integer (5).", 1)       
       WRITE(stderr,*) TRIM(number)
       RETURN
    END IF
    number_ = " "
    number_(1:) = mpc_conv_table(i) // TRIM(number(3:))
    number = " "
    number = TRIM(number_)

  END SUBROUTINE encodeMPCNumber





  !! *Description*:
  !!
  !! Encodes a decoded MPC type designation to the MPC3 format.
  !!
  !! Example:
  !!
  !! 2004 AS1 -> K04A0001S 
  !!
  !! Returns error.
  !!
  SUBROUTINE encodeMPC3Designation(designation)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(inout) :: designation
    CHARACTER(len=DESIGNATION_LEN) :: designation_
    INTEGER :: i

    CALL removeLeadingBlanks(designation)
    i = IACHAR(designation(1:1))
    IF (INDEX(TRIM(designation)," ") == 0 .OR. &
         i < 48 .OR.  i > 57) THEN
       ! Designation already seems to be encoded, because it has no
       ! empty space or starts with a non-number.
       RETURN
    END IF

    IF (INDEX(TRIM(designation),"-") /= 0) THEN
       ! Designation of type 1234 A-B:
       designation_(1:1) = designation(6:6)
       designation_(2:2) = designation(8:8)
       designation_(3:9) = "S00" // designation(1:4)
       designation = designation_(1:9)
       designation(10:) = " "
    ELSE
       ! Designation of type 2001 AB12:
       CALL toInt(designation(1:2), i, error)
       IF (error) THEN
          CALL errorMessage("Base / encodeMPC3Designation", &
               "Could not convert string to integer (5).", 1)       
          WRITE(stderr,*) TRIM(designation)
          RETURN
       END IF
       designation_ = TRIM(ADJUSTL(designation(8:)))
       DO WHILE (LEN_TRIM(designation_) < 4) 
          designation_ = "0" // TRIM(designation_)
       END DO
       designation_ = mpc_conv_table(i) // designation(3:4) // &
            designation(6:6) // TRIM(designation_)
       designation_(9:9) = designation(7:7)
       designation = " "
       designation = TRIM(designation_)
    END IF

  END SUBROUTINE encodeMPC3Designation





  !! *Description*:
  !!
  !! This routine returns the date and time in Greenwich England.
  !! date(1) is day, date(2) is month (Jan=1,...,Dec=12), date(3) is
  !! year, t(1) is hour, t(2) is minute, and t(3) is second.
  !!
  SUBROUTINE coordinatedUniversalTime(date, time)

    IMPLICIT NONE
    INTEGER, DIMENSION(3), INTENT(out) :: date, time
    REAL(bp)                           :: jd, s
    INTEGER, DIMENSION(8)              :: date_time

    CALL DATE_AND_TIME(values=date_time)
    CALL calendarDateToJulianDate(date_time(1), date_time(2), date_time(3), &
         date_time(5), date_time(6), REAL(date_time(7), bp), jd)
    ! Change jd to get correct UTC:
    jd = jd - date_time(4)/(60.0_bp*24.0_bp)
    CALL julianDateToCalendarDate(jd, date_time(1), date_time(2), date_time(3), &
         date_time(5), date_time(6), s)
    date_time(7) = NINT(s)
    ! If seconds rounds to 60 then add half a second to the Julian date
    ! and recalculate the calendar date:
    IF (date_time(7) == 60) THEN
       jd = jd + 0.5_bp/86400.0_bp
       CALL julianDateToCalendarDate(jd, date_time(1), date_time(2), date_time(3), &
            date_time(5), date_time(6), s)
       date_time(7) = NINT(s)
    ENDIF
    date(1) = date_time(3)
    date(2) = date_time(2)
    date(3) = date_time(1)
    time(1) = date_time(5)
    time(2) = date_time(6)
    time(3) = date_time(7)

  END SUBROUTINE coordinatedUniversalTime





END MODULE Base_cl
