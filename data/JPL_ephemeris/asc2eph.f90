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
!! ASC2EPH creates a direct access JPL Planetary Ephemeris file from
!! one or more ascii text files.
!!
!! This program, 'asc2eph', requires (via standard input) an ascii
!! header file ('header.XXX'), followed by one or more ascii ephemeris 
!! data files ('ascSYYYY.XXX').  All files must have the same ephemeris
!! number, XXX.  Further, the data files must be consecutive in time
!! with no gaps between them. 
!!
!! By default, the output ephemeris will span the same interval as the input
!! text file(s).  If you are interested in only a portion of data, set the
!! below T1 and T2 to the begin and end times of the span you desire.  T1
!! and T2 must be specified in  Julian Ephemeris Days (ET).
!!
!! A sample sequence of files might be:
!!
!!   header.405  asc+1920.405 asc+1940.405 asc+1960.405 asc+1980.405  
!!
!! (The data files for DE200 and DE405 contain 20 years each; for DE406, 100 years)
!!
!! @author  MG
!! @version 2015-11-04
!!
PROGRAM asc2eph

  USE cl_options
  IMPLICIT NONE

  !                               *** NOTE ***
  !
  ! The units in which the length of a direct access record is specified 
  ! are PROCESSOR DEPENDENT. The parameter NRECL, the number of units per word, 
  ! controls the length of a record in the direct access ephemeris.
  !
  !*****  The user must choose one of the following statements  *****
  !
  INTEGER, PARAMETER :: NRECL = 4
  !INTEGER, PARAMETER :: NRECL = 1
  !
  !******************************************************************
  !
  ! This parameter is the 'kind'-attribute of the real value to be read
  ! from the ascii files. The value eight (8) corresponds to double
  ! precision numbers, which are used by JPL.
  INTEGER, PARAMETER :: dp = 8 
  !
  ! This parameter is the 'kind'-attribute of the binary real value, to 
  ! which the ascii values are converted. The value eight (8) corresponds 
  ! to bp (base precision), which are used by the JPL_ephemeris routine.
  INTEGER, PARAMETER :: bp = 8 

  CHARACTER(len=6), DIMENSION(400) :: cnam
  CHARACTER(len=6), DIMENSION(14,3) :: ttl
  CHARACTER(len=12) :: header
  CHARACTER(len=3) :: eph_type

  REAL(kind=bp) :: au, emrat, t1, t2
  REAL(kind=bp), DIMENSION(400) :: cval
  REAL(kind=bp), DIMENSION(3) :: ss
  REAL(kind=bp), DIMENSION(3000) :: db
  REAL(kind=bp) :: db2z = 0.0_bp
  REAL(kind=dp), DIMENSION(400) :: cval_dp
  REAL(kind=dp), DIMENSION(3) :: ss_dp
  REAL(kind=dp), DIMENSION(3000) :: db_dp

  INTEGER, DIMENSION(3,12) :: ipt
  INTEGER, DIMENSION(3) :: lpt
  INTEGER :: i, j, k
  INTEGER :: in, out, ksize, irecsz
  INTEGER :: n, ncon, nrout, ncoeff, nrw, numde

  LOGICAL :: first = .TRUE.


  !By default, the output ephemeris will span the same interval as the
  !input ascii data file(s).  The user may reset these to other JED's.
  t1 = -HUGE(t1)
  t2 = HUGE(t2)

  eph_type = get_cl_option("--eph-type=", "405")

  ! Write a fingerprint to the screen.
  WRITE(*,*) " JPL ASCII-TO-DIRECT-I/O program. Last modified 4-November-2015."
  WRITE(*,*) 
  WRITE(*,*) " ASSUMING TYPE OF INPUT ASCII EPHEMERIS FILE IS DE" // TRIM(eph_type)

  IF (nrecl == 0) THEN
     STOP "*** ERROR: User did not set NRECL ***"
  END IF

  ! Read the size and number of main ephemeris records.
  READ(*,"(6X,I6)") ksize
  WRITE(*,'("  KSIZE = ",I0)') ksize

  irecsz = nrecl * ksize

  ! Now for the alphameric heading records (GROUP 1010)
  CALL nxtgrp(header)
  IF (header /= "GROUP   1010") CALL errprt(1010, "not header")
  READ(*,"(14A6)") ttl
  WRITE(*,"(/(2X,14A6))") ttl


  ! Read start, end and record span  (GROUP 1030)
  CALL nxtgrp(header)
  IF (header /= "GROUP   1030") CALL errprt(1030, " not header")
  READ(*,"(3D12.0)") ss_dp
  ss = REAL(ss_dp,bp)


  ! Read number of constants and names of constants (GROUP 1040/4).
  CALL nxtgrp(header)
  IF (header /= "GROUP   1040 ") CALL errprt(1040, " not header")
  READ(*,"(I6)") n
  READ(*,"(10A8)") cnam(1:n)
  ncon = n


  ! Read number of values and values (GROUP 1041/4)
  CALL nxtgrp(header)
  IF (header /= "GROUP   1041") CALL errprt(1041, " not header")
  READ(*,"(I6)") n
  READ(*,"(3D26.18)") cval_dp(1:n)
  cval = REAL(cval_dp,bp)
  DO  i = 1, n
     IF (cnam(i) .EQ. "AU    ")  au    = cval(i)
     IF (cnam(i) .EQ. "EMRAT ")  emrat = cval(i)
     IF (cnam(i) .EQ. "DENUM ")  numde = cval(i)
  END DO
  WRITE(*,"(500(/2(A8,D24.16)))") (cnam(i),cval(i),i=1,n)


  ! Read pointers needed by interp (GROUP 1050)
  CALL nxtgrp(header)
  IF (header /= "GROUP   1050") CALL errprt(1050, " not header")
  READ(*,"(13I6)") ((ipt(i,j),j=1,12),lpt(i),i=1,3)
  WRITE(*,"(/(3I5))") ipt, lpt


  ! Open direct-access output file ('deXXX.dat')

  OPEN(unit=12, file="de" // TRIM(eph_type) // ".dat", &
       access="DIRECT", form="UNFORMATTED", recl=irecsz, status="NEW")

  ! Read and write the ephemeris data records (GROUP 1070).
  CALL nxtgrp(header)
  IF (header /= "GROUP   1070") CALL errprt(1070," not header")
  nrout = 0
  in = 0
  out = 0

  nrw = 0
  DO WHILE (nrw == 0)
     READ(*,"(2I6)") nrw, ncoeff
  END DO

  READ(*,*,iostat=in) db_dp(1:ncoeff)
  db = REAL(db_dp,bp)

  DO WHILE (in == 0 .AND. db(2) < t2)

     IF (2*ncoeff /= ksize) THEN
        CALL errprt(ncoeff, " 2*ncoeff not equal to ksize")
     END IF

     ! Skip this data block if the end of the interval is less
     ! than the specified start time or if the it does not begin
     ! where the previous block ended.
     IF ((db(2) >= t1) .AND. (db(1) >= db2z)) THEN

        IF (first) THEN
           ! Don't worry about the intervals overlapping
           ! or abutting if this is the first applicable
           ! interval.
           db2z  = db(1)
           first = .FALSE.
        END IF

        IF (db(1) /= db2z) THEN
           ! Beginning of current interval is past the end
           ! of the previous one.
           CALL errprt (nrw, " Records do not overlap or abut")
        END IF

        db2z = db(2)
        nrout = nrout + 1

        WRITE(12, rec=nrout+2, iostat=out) db(1:ncoeff)
        IF (out /= 0) CALL errprt(nrout, "th record not written because of error")

        ! Save this block's starting date, its interval span, and its end 
        ! date.
        IF (nrout == 1) THEN
           ss(1) = db(1)
           ss(3) = db(2) - db(1)
        END IF
        ss(2) = db(2)

        ! Update the user as to our progress every 10th block.
        IF (MOD(nrout,10) == 1) THEN
           IF (db(1) >= t1) THEN
              WRITE(*,'(I6," ephemeris records written. Last jed = ", F12.2)') nrout, db(2)
           ELSE
              WRITE(*,*) " Searching for first requested record..."
           END IF
        END IF
     END IF

     !READ(*,"(2I6,3000(/3D26.18))",iostat=in) nrw, ncoeff, db_dp(1:ncoeff)
     READ(*,"(2I6)",iostat=in) nrw, ncoeff
     IF (in == 0) THEN
        READ(*,*,iostat=in) db_dp(1:ncoeff)
        db = REAL(db_dp,bp)
     END IF

  END DO

  WRITE(*,'(/I6," ephemeris records written. Last jed = ", F12.2)') nrout, db(2)


  ! Write header records onto output file.
  nrout = 1
  WRITE(12, rec=1, iostat=out) ttl, cnam, ss, ncon, au, emrat, ipt, numde, lpt
  IF (out /= 0) CALL errprt(nrout, "th record not written because of error")

  nrout = 2
  WRITE(12, rec=2, iostat=out) cval
  IF (out /= 0) CALL errprt(nrout, "th record not written because of error")

  ! We're through.  Wrap it up.
  CLOSE (12)
  STOP "ASCII to binary conversion successful."



CONTAINS



  SUBROUTINE errprt(i, msg)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: msg
    INTEGER, INTENT(in) :: i

    WRITE(*,'("ERROR #",I8,2X,A50)') i, msg
    STOP " ERROR "

  END SUBROUTINE errprt



  SUBROUTINE nxtgrp(header)

    IMPLICIT NONE
    CHARACTER(len=*), INTENT(out) :: header
    CHARACTER(len=12) :: blank

    !Start with nothing.
    header = " "

    !The next non-blank line we encounter is a header record.
    !The group header and data are seperated by a blank line.
    DO WHILE (header == " ")
       READ(*,"(A)") header
    END DO

    !Found the header. Read the blank line so we can get at the data.
    IF (header /= "GROUP   1070") THEN
       READ(*,"(A)") blank
    END IF

  END SUBROUTINE nxtgrp



END PROGRAM asc2eph
