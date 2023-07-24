PROGRAM stats


  USE Base_cl
  USE Time_cl
  USE cl_options
  USE utilities
  USE stats_err
  
  IMPLICIT NONE

  TYPE (StatsErr), DIMENSION(:), POINTER :: ses
  INTEGER :: i
  LOGICAL :: err = .FALSE.

  TYPE (Time) :: tim
  REAL(bp) :: tmp
  tmp = getMJD(tim, "UTC")
  
  CALL stats_errors_init(err)

  PRINT *, "err: ", err
  PRINT *, "size: ", stats_errs_size

  ses => stats_errors("Y00", 1)
  PRINT *, "size-of-ses: ", SIZE(ses)
  DO i = 1, SIZE(ses)
     print *, "stats_error", ses(i)
  END DO
  PRINT *, "size-of-ses: ", SIZE(ses)
  DEALLOCATE(ses)
  NULLIFY(ses)
  
END PROGRAM stats
