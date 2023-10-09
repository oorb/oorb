MODULE stats_err

  USE Time_cl
  IMPLICIT NONE

  PUBLIC :: stats_errors  
  PUBLIC :: stats_errors_init
  PUBLIC :: stats_err_file
  PUBLIC :: stats_errs       
  PUBLIC :: stats_errs_size
  PUBLIC :: StatsErr

  PRIVATE
  LOGICAL :: first_time = .FALSE.
  LOGICAL :: stats_err_file = .FALSE.
  
  TYPE StatsErr
     SEQUENCE
     CHARACTER(len = 3) :: obs
     INTEGER            :: start, end
     INTEGER            :: low_mag, upper_mag
     REAL               :: ra_rms, dec_rms
  END TYPE StatsErr

  CHARACTER(len=FNAME_LEN), PARAMETER :: STATS_FNAME = 'STATS-ERR.dat'
  TYPE (StatsErr), DIMENSION (:), ALLOCATABLE :: stats_errs
  INTEGER :: stats_errs_size
  
CONTAINS

  FUNCTION stats_errors(obs, mag, tim)

    IMPLICIT NONE

    CHARACTER(len = 3), INTENT(in)    :: obs
    INTEGER, INTENT(in)               :: mag
    TYPE (Time), INTENT(in), OPTIONAL :: tim

    TYPE (StatsErr), DIMENSION(:), POINTER :: stats_errors
    INTEGER  :: i, j, err
    TYPE (Time) :: tim_    
    REAL(bp) :: jdt

    if (PRESENT(tim)) THEN
       tim_ = tim
       jdt = getJD(tim_, "UTC")
    END IF
    
    ! todo: make the redundant looping more efficient later--correctness first
    j = 0
    DO i = 1, stats_errs_size
       IF (stats_errs(i)%obs == obs .AND. stats_errs(i)%low_mag <= mag .AND. mag < stats_errs(i)%upper_mag) THEN
          IF (PRESENT(tim)) THEN
             IF (real(stats_errs(i)%start) <= jdt .AND. jdt < real(stats_errs(i)%end)) then
                j = j + 1
             END IF
          ELSE
             j = j + 1
          END IF
       END IF
    END DO

    if (j > 0) THEN
       ALLOCATE(stats_errors(1:j), STAT = err) ! *todo: test for error later
       j = 0
       DO i = 1, stats_errs_size
          IF (stats_errs(i)%obs == obs .AND. stats_errs(i)%low_mag <= mag .AND. mag < stats_errs(i)%upper_mag) THEN
             IF (PRESENT(tim)) THEN
                IF (real(stats_errs(i)%start) <= jdt .AND. jdt < real(stats_errs(i)%end)) then
                   j = j + 1
                   stats_errors(j) = stats_errs(i)                   
                END IF
             ELSE
                j = j + 1
                stats_errors(j) = stats_errs(i)
             END IF
          END IF
       END DO
    ELSE
       if (stats_err_file) THEN
          ! get the last row for "obs"--it's sorted to begin with
          ALLOCATE(stats_errors(1), STAT = err) ! *todo: test for error later
          DO i = 1, stats_errs_size
             IF (stats_errs(i)%obs == obs) THEN
                stats_errors(1) = stats_errs(i)
             END IF
          END DO
       ELSE
          ALLOCATE(stats_errors(0), STAT = err) ! *todo: test for error later          
       END IF
    END IF

  END FUNCTION stats_errors

  
  SUBROUTINE stats_errors_init(error, filename)

    IMPLICIT NONE
    LOGICAL, INTENT(inout)                   :: error
    CHARACTER(len = *), OPTIONAL, INTENT(in) :: filename
    CHARACTER(len = FNAME_LEN)               :: fname, OORB_DATA_DIR    

    INTEGER :: i, fid = 11, err = 0, lines = 0
    CHARACTER(len=256) :: ctmp
    INTEGER :: itmp
    REAL :: rtmp    
    

    IF (first_time) THEN
       RETURN
    END IF

    if (PRESENT(filename) .AND. LEN_TRIM(filename) <= FNAME_LEN) THEN
       fname = TRIM(filename)
    ELSE
       CALL getenv("OORB_DATA", OORB_DATA_DIR)
       IF (LEN_TRIM(OORB_DATA_DIR) == 0) THEN
          OORB_DATA_DIR = "."
       END IF
       fname = TRIM(OORB_DATA_DIR) // "/" // TRIM(STATS_FNAME)
    END IF
    
    OPEN(fid, file = fname, status = "old")

    READ(fid, *) ! skip the first 2 lines (headers)
    READ(fid, *)    
    DO WHILE (err == 0)
       lines = lines + 1
       READ(fid, *, iostat = err) ctmp
    END DO
    lines = lines - 1
    write(*, '(A,I0)') "total number of lines: ", lines

    ALLOCATE(stats_errs(lines))
    
    REWIND(fid)
    READ(fid, *)
    READ(fid, *)
    DO i = 1, lines
       READ(fid, *) &
            stats_errs(i)%obs, &
            stats_errs(i)%start, &
            stats_errs(i)%end, &
            stats_errs(i)%low_mag, &
            stats_errs(i)%upper_mag, &
            ctmp, &
            ctmp, &
            stats_errs(i)%ra_rms, &
            stats_errs(i)%dec_rms, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            itmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp, &
            rtmp
       IF (stats_errs(i)%obs == "Y28") EXIT ! some fields missing after obs="Y28"
    END DO

    stats_errs_size = SIZE(stats_errs) - 4 ! last 4 lines are missing data columns

    ! PRINT *, "stats_errs_size=", stats_errs_size
    
    ! DO i = 1, stats_errs_size
    !    print *, &
    !         stats_errs(i)%obs, &
    !         stats_errs(i)%start, &
    !         stats_errs(i)%end, &                        
    !         stats_errs(i)%low_mag, &
    !         stats_errs(i)%upper_mag, &
    !         stats_errs(i)%ra_rms, &
    !         stats_errs(i)%dec_rms
    ! END DO
    
    first_time = .TRUE.
  END SUBROUTINE stats_errors_init
  
END MODULE stats_err

