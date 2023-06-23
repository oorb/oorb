MODULE stats_err

  USE parameters
  IMPLICIT NONE

  PUBLIC :: stats_err_init
  PUBLIC :: StatsErr

  PRIVATE
  LOGICAL :: first = .TRUE.
  
  TYPE StatsErr
     SEQUENCE
     CHARACTER(len = 3) :: obs
     INTEGER            :: start, end
     INTEGER            :: low_mag, upper_mag
     REAL               :: ra_rms, dec_rms
  END TYPE StatsErr

  CHARACTER(len=FNAME_LEN), PARAMETER :: STATS_FNAME = 'STATS-ERR.dat'
  TYPE (StatsErr), dimension (:), allocatable :: stats_errs
  
CONTAINS
  
  SUBROUTINE stats_err_init(error, filename)

    IMPLICIT NONE
    LOGICAL, INTENT(inout)                   :: error
    CHARACTER(len = *), OPTIONAL, INTENT(in) :: filename
    CHARACTER(len = FNAME_LEN)               :: fname, OORB_DATA_DIR    

    INTEGER :: i, fid = 11, err = 0, lines = 0
    CHARACTER(len=256) :: tmp

    IF (.NOT.first) THEN
       RETURN
    END IF
    
    OPEN(fid, file = STATS_FNAME, status = "old")

    READ(fid, *) ! skip 2 lines of headers
    READ(fid, *)    
    DO WHILE (err == 0)
       lines = lines + 1
       READ(fid, *, iostat = err) tmp
    END DO
    lines = lines - 1
    write(*, '(A,I0)') "total number of lines: ", lines

    ALLOCATE(stats_errs(lines))
    
    REWIND(fid)
    READ(fid, *)
    READ(fid, *)

    DO i = 1, lines
       READ(fid, *) stats_errs(i)%obs
    END DO
    
  END SUBROUTINE stats_err_init
  
END MODULE stats_err

