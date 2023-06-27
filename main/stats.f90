PROGRAM stats

  USE cl_options
  USE Base_cl
  USE Time_cl
  USE stats_err
  
  IMPLICIT NONE

  LOGICAL :: err = .TRUE.
  
  PRINT *, "Hello, World!"

  CALL stats_err_init(err)
  
END PROGRAM stats
