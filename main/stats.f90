PROGRAM stats

  USE cl_options
  USE Base_cl
  USE Time_cl
  USE stats_err
  
  IMPLICIT NONE

  TYPE (StatsErr), DIMENSION(:), POINTER :: ses
  INTEGER :: i
  LOGICAL :: err = .TRUE.
  CALL stats_errors_init(err)

  PRINT *, "err: ", err
  PRINT *, "size: ", SIZE(stats_errs)

  ses => stats_errors(1)
  PRINT *, "size-of-ses: ", SIZE(ses)
  DO i = 1, SIZE(ses)
     print *, "stats_error", ses(i)
  END DO
  PRINT *, "size-of-ses: ", SIZE(ses)
  DEALLOCATE(ses)
  NULLIFY(ses)
  
END PROGRAM stats
