!====================================================================!
!                                                                    !
! Copyright 2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012   !
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
!! Defines a File-object and the basic file handling routines like open, 
!! close, read (character string), write (character string) etc. 
!!
!! *Example*:
!!
!! <pre>
!!program myprog
!!
!! use File_cl
!! implicit none
!! type (File) :: f1, f2
!! character(len=32) :: str
!!
!! call new(f1, "myfile1.txt")
!! if (error) stop "Error 1"
!!
!! call open(f1)
!! if (error) stop "Error 2"
!!
!! frmt = "(A20)"
!! call readString(f1, trim(frmt), str)
!! if (error) stop "Error 3"
!!
!! call nullify(f1)
!! if (error) stop "Error 4"
!!
!! call new(f2, "myfile2.txt")
!! if (error) stop "Error 5"
!!
!! call open(f2)
!! if (error) stop "Error 6"
!!
!! frmt = "('Copy from myfile1.txt: ',A)"
!! call writeString(f2, trim(frmt), str)
!! if (error) stop "Error 7"
!!
!! call nullify(f2)
!! if (error) stop "Error 8"
!!
!!end program myprog
!! </pre>
!!
!! @author  MG
!! @version 2012-01-17
!! 
!! @see Observations_class
!! @see Time_class
!!
MODULE File_cl

  USE Base_cl
  USE Unit_cl

  USE utilities

  IMPLICIT NONE
  PRIVATE :: new_F_scratch
  PRIVATE :: new_F_name
  PRIVATE :: nullify_F
  PRIVATE :: copy_F
  PRIVATE :: writeString_unformatted
  PRIVATE :: writeString_formatted
  PRIVATE :: readString_unformatted
  PRIVATE :: readString_formatted
  PRIVATE :: getUnit_F

  TYPE File
     PRIVATE
     TYPE (Unit) :: lu                   ! Logical unit.
     CHARACTER(len=FNAME_LEN) :: fname   ! File name.
     CHARACTER(len=16) :: open_access    ! File access (SEQUENTIAL/DIRECT).
     CHARACTER(len=16) :: open_action    ! Action (READWRITE/READ/WRITE).
     CHARACTER(len=16) :: open_form      ! Form (FORMATTED/UNFORMATTED/BINARY).
     CHARACTER(len=16) :: open_position  ! Position (ASIS/APPEND/REWIND).
     CHARACTER(len=16) :: open_status    ! File status (UNKOWN/SCRATCH/OLD/NEW).
     INTEGER :: open_recl                ! Record length in bytes (max length for sequential files).
     LOGICAL :: opened         = .FALSE. ! File opened? (true=yes, false=no) 
     LOGICAL :: is_initialized = .FALSE. ! Object initialized? (true=yes, false=no)
  END TYPE File

  !! Initializes a File-object. Can be 
  !! used with or without a filename.
  INTERFACE NEW
     MODULE PROCEDURE new_F_scratch
     MODULE PROCEDURE new_F_name
  END INTERFACE NEW

  !! Nullifies a File-object, e.g. closes the file.
  INTERFACE NULLIFY
     MODULE PROCEDURE nullify_F
  END INTERFACE NULLIFY

  !! Makes a copy of a File-object.
  INTERFACE copy
     MODULE PROCEDURE copy_F
  END INTERFACE copy

  !! Returns the status of the object, i.e. whether
  !! it exists or not.
  INTERFACE exist
     MODULE PROCEDURE exist_F
  END INTERFACE exist

  !! Writes a formatted or list-directed 
  !! character string to a file.
  INTERFACE writeString
     MODULE PROCEDURE writeString_unformatted
     MODULE PROCEDURE writeString_formatted
  END INTERFACE writeString

  !! Reads a formatted or list-directed 
  !! character string from a file.
  INTERFACE readString
     MODULE PROCEDURE readString_unformatted
     MODULE PROCEDURE readString_formatted
  END INTERFACE readString

  !! Returns logical unit of a file object.
  INTERFACE getUnit
     MODULE PROCEDURE getUnit_F
  END INTERFACE getUnit





CONTAINS





  !! *Description*:
  !!
  !! Initializes a simple File instance, i.e. a scratchfile with 
  !! sequential access, and both read and write permissions.
  !!
  !! Returns error.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call new(myfile)
  !!
  SUBROUTINE new_F_scratch(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / new", &
            "Object has already been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / new", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%fname          = ""
    this%open_access    = "sequential"
    this%open_action    = "readwrite"
    this%open_form      = "formatted"
    this%open_position  = "asis"
    this%open_recl      = -1
    this%open_status    = "scratch"
    this%opened         = .FALSE.
    this%is_initialized = .TRUE.

  END SUBROUTINE new_F_scratch





  !! *Description*:
  !!
  !! Initializes a new File instance. It's a named file (either old
  !! or new) with sequential access, and both read and write permissions.
  !! Returns _error_ =
  !!     - _false_, if everything is ok.
  !!     - _true_, if the file name is too long (more than 64 letters)
  !! 
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call new(myfile, "myfile.dat")
  !!
  SUBROUTINE new_F_name(this, file_name)

    IMPLICIT NONE
    TYPE (File), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: file_name

    CALL NEW(this)
    IF (error) THEN
       CALL errorMessage("File / new", &
            "TRACE BACK", 1)
       RETURN
    END IF
    this%is_initialized = .FALSE.
    this%open_status = "unknown"
    IF (LEN_TRIM(file_name) > FNAME_LEN) THEN
       error = .TRUE.
       CALL errorMessage("File / new", &
            "Filename (" // TRIM(file_name) // &
            ") too long; adjust parameters.", 1)
       RETURN
    ELSE
       this%fname    = TRIM(file_name)
    END IF

    this%is_initialized = .TRUE.

  END SUBROUTINE new_F_name





  !! *Description*:
  !!
  !! Closes a file. Returns _error_ =
  !!     - _false_, if the file is properly closed, or if it 
  !!                has not been opened.
  !!     - _true_,  if an error occurs during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call closeFile(myfile) 
  !!
  SUBROUTINE nullify_F(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (this%opened) THEN
       CALL NULLIFY(this%lu)
       IF (error) THEN
          CALL errorMessage("File / nullify", &
               "TRACE BACK", 1)
          RETURN
       END IF
    END IF
    this%opened         = .FALSE.
    this%is_initialized = .FALSE.

  END SUBROUTINE nullify_F





  !! *Description*:
  !!
  !! Returns a copy of this object.
  !!
  FUNCTION copy_F(this)

    IMPLICIT NONE
    TYPE (File) :: this
    TYPE (File) :: copy_F

    copy_F%lu             = copy(this%lu)
    copy_F%fname          = this%fname
    copy_F%open_status    = this%open_status
    copy_F%open_access    = this%open_access
    copy_F%open_action    = this%open_action
    copy_F%open_form      = this%open_form
    copy_F%open_position  = this%open_position
    copy_F%open_recl      = this%open_recl
    copy_F%opened         = this%opened
    copy_F%is_initialized = this%is_initialized

  END FUNCTION copy_F





  !! *Description*:
  !!
  !! Returns the status of the object, i.e. whether
  !! it exists or not.
  !!
  LOGICAL FUNCTION exist_F(this)

    IMPLICIT NONE
    TYPE (File) :: this

    exist_F = this%is_initialized

  END FUNCTION exist_F





  !! *Description*:
  !!
  !! Returns the name of the file connected to this object. Returns
  !! _error_ =
  !!     - _false_, no errors occur.
  !!     - _true_, if an error occurs during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! character(len=FNAME_LEN) :: fname
  !!
  !! ...
  !!
  !! fname = getFileName(myfile)
  !!
  CHARACTER(len=FNAME_LEN) FUNCTION getFileName(this)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / getFileName", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    getFileName = TRIM(this%fname)

  END FUNCTION getFileName





  !! *Description*:
  !!
  !! Counts the number of columns in a file based on the first row and
  !! assuming an empty space as the column delimiter. Returns _error_ =
  !!     - _false_, no errors occur.
  !!     - _true_, if an error occurs during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! integer :: nr_of_columns
  !!
  !! ...
  !!
  !! nr_of_columns = getNrOfColumns(myfile, error) 
  !!
  INTEGER FUNCTION getNrOfColumns(this)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: this
    CHARACTER(len=4096)     :: line, line_
    INTEGER                 :: i, err, indx

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfColumns", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfColumns", &
            "File has not yet been opened.", 1)
       RETURN
    END IF

    REWIND(getUnit(this), iostat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfColumns", &
            "Error while rewinding the file.", 1)
       RETURN
    END IF

    line = " "
    READ(getUnit(this), "(A)", iostat=err) line
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfColumns", &
            "Error while reading first line of file " // TRIM(this%fname), 1)
       RETURN
    END IF

    REWIND(getUnit(this), iostat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfColumns", &
            "Error while rewinding the file.", 1)
       RETURN
    END IF

    i = 0
    DO WHILE (LEN_TRIM(line) /= 0)
       Line_ = " "
       IF (line(1:1) /= " ") THEN
          i = i + 1
          indx = INDEX(line," ")
          line_ = line(indx:)
          line = " "
          line = line_
       ELSE
          line_ = line(2:)
          line = " "
          line = line_          
       END IF
    END DO

    getNrOfColumns = i

  END FUNCTION getNrOfColumns





  !! *Description*:
  !!
  !! Counts the number of lines in a file. Returns _error_ =
  !!     - _false_, no errors occur.
  !!     - _true_, if an error occurs during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! integer :: nr_of_lines
  !!
  !! ...
  !!
  !! nr_of_lines = getNrOfLines(myfile, error) 
  !!
  INTEGER FUNCTION getNrOfLines(this)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: this
    CHARACTER(len=5)        :: line
    INTEGER                 :: i, err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfLines", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfLines", &
            "File has not yet been opened.", 1)
       RETURN
    END IF

    REWIND(getUnit(this), iostat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfLines", &
            "Error while rewinding the file.", 1)
       RETURN
    END IF

    i = 0
    DO
       READ(getUnit(this), "(a)", iostat=err) line
       IF (err == 0) THEN
          i = i + 1
       ELSE IF (err < 0) THEN
          ! It's the end-of-file: 
          EXIT
       ELSE
          error = .TRUE.
          CALL errorMessage("File / getNrOfLines", &
               "Error while counting lines in file.", 1)
          RETURN
       END IF
    END DO

    REWIND(getUnit(this), iostat=err)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / getNrOfLines", &
            "Error while rewinding the file.", 1)
       RETURN
    END IF

    getNrOfLines = i

  END FUNCTION getNrOfLines





  !! *Description*:
  !!
  !! Returns
  !!  - record length, if the file has been opened for direct access.
  !!  - maximum record length, if the file has been opened for
  !!    sequential access.
  !!  - 0, if the file has not been opened, or it does not exist.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! integer :: reclen
  !!
  !! ...
  !!
  !! recl = getRecordLength(myfile) 
  !! 
  INTEGER FUNCTION getRecordLength(this)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: this
    INTEGER                 :: record_length, err

    IF (.NOT. this%is_initialized) THEN
       CALL errorMessage("File / getRecordLength", &
            "Object has not yet been initialized.", 1)
       getRecordLength = 0
       RETURN
    END IF

    IF (LEN_TRIM(this%fname) /= 0) THEN 
       INQUIRE(file=TRIM(this%fname), recl=record_length, iostat=err)
       IF (err == 0) THEN
          getRecordLength = record_length
       ELSE
          getRecordLength = 0
       END IF
    ELSE
       getRecordLength = 0
    END IF

  END FUNCTION getRecordLength





  !! *Description*:
  !!
  !! Returns
  !!     - logical unit number, if the file has been opened.
  !!     - -1, if the file has not been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! integer :: logical_unit
  !!
  !! ...
  !!
  !! logical_unit = getUnit(myfile) 
  !!
  INTEGER FUNCTION getUnit_F(this)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: this

    IF (.NOT. this%is_initialized) THEN
       getUnit_F = -1
    ELSE IF (.NOT. this%opened) THEN
       getUnit_F = -1
    ELSE
       getUnit_F = getUnit(this%lu)
    END IF

  END FUNCTION getUnit_F





  !! *Description*:
  !!
  !! Returns
  !!     - _true_, if the file has been opened.
  !!     - _false_, if the file has not been opened
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! logical :: open
  !!
  !! ...
  !!
  !! open = isOpen(myfile) 
  !!
  LOGICAL FUNCTION isOpen(this)

    IMPLICIT NONE
    TYPE (File), INTENT(in) :: this

    isOpen = this%opened

  END FUNCTION isOpen





  !! *Description*:
  !!
  !! Opens a file. Open specifiers can be altered with the 
  !! _set_ routines (see below). Returns _error_ =
  !!  - _false_, if the file is opened, or if it already has been
  !!             opened.
  !!  - _true_,  if an error occurs during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call openFile(myfile) 
  !!
  SUBROUTINE OPEN(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this
    CHARACTER(len=8)           :: str
    INTEGER                    :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / open", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       CALL NEW(this%lu)
       IF (error) THEN
          CALL errorMessage("File / open", &
               "TRACE BACK", 1)
          RETURN
       END IF
       IF (this%open_recl == -1) THEN
          IF (this%open_position == "none") THEN
             OPEN(unit=getUnit(this%lu), file=TRIM(this%fname), form=TRIM(this%open_form), &
                  access=TRIM(this%open_access), status=TRIM(this%open_status), &
                  action=TRIM(this%open_action), iostat=err)
          ELSE
             OPEN(unit=getUnit(this%lu), file=TRIM(this%fname), form=TRIM(this%open_form), &
                  access=TRIM(this%open_access), status=TRIM(this%open_status), &
                  position=TRIM(this%open_position), action=TRIM(this%open_action), iostat=err)
          END IF
       ELSE
          IF (this%open_position == "none") THEN
             OPEN(unit=getUnit(this%lu), file=TRIM(this%fname), form=TRIM(this%open_form), &
                  access=TRIM(this%open_access), status=TRIM(this%open_status), &
                  recl=this%open_recl, action=TRIM(this%open_action), iostat=err)
          ELSE
             OPEN(unit=getUnit(this%lu), file=TRIM(this%fname), form=TRIM(this%open_form), &
                  access=TRIM(this%open_access), status=TRIM(this%open_status), &
                  position=TRIM(this%open_position), recl=this%open_recl, &
                  action=TRIM(this%open_action), iostat=err)
          END IF
       END IF
       IF (err == 0) THEN
          this%opened = .TRUE.
       ELSE
          this%opened = .FALSE.
          CALL NULLIFY(this%lu)
          CALL toString(err, str, error)
          error = .TRUE.
          CALL errorMessage("File / open", &
               "Could not open file (" // TRIM(this%fname) // &
               "). Error code " // TRIM(str) // ".", 1)
          RETURN     
       END IF
    END IF

  END SUBROUTINE OPEN





  !! *Description*:
  !!
  !! Reads a formatted (_len_ = *) character string (_len_ = *) from 
  !! file and moves to the following line. Returns _error_ = 
  !!     - _false_, if the file has been opened and 
  !!       no errors occured during the procedure.
  !!     - _true_, if the file has not been opened, or
  !!       an error occured during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! character(len=20) :: mystring
  !!
  !! character(len=20) :: myformat = "(a30)"
  !!
  !! ...
  !!
  !! call readString(myfile, trim(myformat), mystring) 
  !!
  SUBROUTINE readString_formatted(this, frmt, str)

    IMPLICIT NONE
    TYPE (File), INTENT(in)       :: this
    CHARACTER(len=*), INTENT(in)  :: frmt
    CHARACTER(len=*), INTENT(out) :: str
    INTEGER                       :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / readString", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / readString", &
            "File has not yet been opened.", 1)
       RETURN
    END IF

    READ(getUnit(this%lu), fmt=TRIM(frmt), iostat=err) str
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / readString", &
            "Could not read string from file.", 1)
       RETURN
    END IF

  END SUBROUTINE readString_formatted





  !! *Description*:
  !!
  !! Reads a singel, list-directed character string (_len_ = *) 
  !! from file and moves to the following line. Returns _error_ = 
  !!     - _false_, if the file has been opened and 
  !!       no errors occured during the procedure.
  !!     - _true_, if the file has not been opened, or
  !!       an error occured during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! character(len=20) :: mystring
  !!
  !! ...
  !!
  !! call readString(myfile, mystring) 
  !!
  SUBROUTINE readString_unformatted(this, str)

    IMPLICIT NONE
    TYPE (File), INTENT(in)       :: this
    CHARACTER(len=*), INTENT(out) :: str
    INTEGER                       :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / readString", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / readString", &
            "File has not yet been opened.", 1)
       RETURN
    END IF

    READ(getUnit(this%lu), *, iostat=err) str
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / readString", &
            "Could not read string from file.", 1)
       RETURN
    END IF

  END SUBROUTINE readString_unformatted





  !! *Description*:
  !!
  !! Allows unformatted direct access with a record length of 
  !! _record___length_ bytes for this File-instance. Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! integer :: reclen = 2
  !!
  !! ...
  !!
  !! call setAccessDirect(myfile, reclen) 
  !!
  SUBROUTINE setAccessDirect(this, record_length)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this
    INTEGER, INTENT(in)        :: record_length

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setAccessDirect", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setAccessDirect", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_access = "direct"
    this%open_form   = "unformatted"
    this%open_recl   = record_length
    this%open_position = "none"

  END SUBROUTINE setAccessDirect





  !! *Description*:
  !!
  !! Allows formatted sequential access for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setAccessSequential(myfile) 
  !!
  SUBROUTINE setAccessSequential(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setAccessSequential", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setAccessSequential", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_access = "sequential"
    this%open_form   = "formatted"
    this%open_recl   = -1


  END SUBROUTINE setAccessSequential





  !! *Description*:
  !!
  !! Sets flag _only_reading_allowed_ for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setActionRead(myfile) 
  !!
  SUBROUTINE setActionRead(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setActionRead", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setActionRead", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_action = "read"

  END SUBROUTINE setActionRead





  !! *Description*:
  !!
  !! Sets flag _both_reading_and_writing_allowed_
  !! for this File-instance. Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setActionReadWrite(myfile) 
  !!
  SUBROUTINE setActionReadWrite(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setActionReadWrite", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setActionReadWrite", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_action = "readwrite"

  END SUBROUTINE setActionReadWrite





  !! *Description*:
  !!
  !! Sets flag _only_writing_allowed_ for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setActionWrite(myfile) 
  !!
  SUBROUTINE setActionWrite(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setActionWrite", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setActionWrite", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_action = "write"

  END SUBROUTINE setActionWrite





  !! *Description*:
  !!
  !! Sets position flag _append_ for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setPositionAppend(myfile) 
  !!
  SUBROUTINE setPositionAppend(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setPositionAppend", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setPositionAppend", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_position = "append"
    this%open_access   = "sequential"

  END SUBROUTINE setPositionAppend





  !! *Description*:
  !!
  !! Sets position flag _asis_ (default) for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setPositionAsis(myfile) 
  !!
  SUBROUTINE setPositionAsis(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setPositionAsis", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setPositionAsis", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_position = "asis"
    this%open_access   = "sequential"

  END SUBROUTINE setPositionAsis





  !! *Description*:
  !!
  !! Sets position flag _rewind_ for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setPositionRewind(myfile)
  !!
  SUBROUTINE setPositionRewind(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setPositionRewind", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setPositionRewind", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_position = "rewind"
    this%open_access   = "sequential"

  END SUBROUTINE setPositionRewind





  !! *Description*:
  !!
  !! Sets flag _file_must_not_exist for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setStatusNew(myfile) 
  !!
  SUBROUTINE setStatusNew(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setStatusNew", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setStatusNew", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_status = "new"

  END SUBROUTINE setStatusNew





  !! *Description*:
  !!
  !! Sets flag _file_must_exist_ for this File-instance. 
  !! Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setStatusOld(myfile) 
  !!
  SUBROUTINE setStatusOld(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setStatusOld", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setStatusOld", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_status = "old"

  END SUBROUTINE setStatusOld





  !! *Description*:
  !!
  !! Sets flag _this_file_replaces_an_existing_file_ 
  !! for this File-instance. Returns _error_ =
  !!     - _false_, if the file has not yet been opened.
  !!     - _true_, if the file has already been opened.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! ...
  !!
  !! call setStatusReplace(myfile) 
  !!
  SUBROUTINE setStatusReplace(this)

    IMPLICIT NONE
    TYPE (File), INTENT(inout) :: this

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / setStatusReplace", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / setStatusReplace", &
            "File has already been opened.", 1)
       RETURN
    END IF

    this%open_status = "replace"

  END SUBROUTINE setStatusReplace





  !! *Description*:
  !!
  !! Writes a formatted (_len_ = *) string (_len_ = *) to 
  !! file and moves to the following line. Returns _error_ =
  !!     - _false_, if the file has been opened and 
  !!       no errors occured during the procedure.
  !!     - _true_, if the file has not been opened, or
  !!       an error occured during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! character(len=20) :: mystring = "Hello!"
  !!
  !! character(len=20) :: myformat = "(a10)"
  !!
  !! ...
  !!
  !! call writeString(myfile, myformat, mystring) 
  !!
  SUBROUTINE writeString_formatted(this, frmt, str)

    IMPLICIT NONE
    TYPE (File), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: frmt
    CHARACTER(len=*), INTENT(in) :: str
    INTEGER                      :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / writeString", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / writeString", &
            "File has not yet been opened.", 1)
       RETURN
    END IF

    WRITE(getUnit(this%lu), TRIM(frmt), iostat=err) TRIM(str)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / writeString", &
            "Could not write string to file.", 1)
       RETURN
    END IF

  END SUBROUTINE writeString_formatted





  !! *Description*:
  !!
  !! Writes a string (_len_ = *) to file and moves 
  !! to the following line. Returns _error_ = 
  !!     - _false_, if the file has been opened and 
  !!       no errors occured during the procedure.
  !!     - _true_, if the file has not been opened, or
  !!       an error occured during the procedure.
  !!
  !! *Usage*:
  !!
  !! type (File) :: myfile
  !!
  !! character(len=20) :: mystring = "Hello!"
  !!
  !! ...
  !!
  !! call writeString(myfile, mystring) 
  !!
  SUBROUTINE writeString_unformatted(this, str)

    IMPLICIT NONE
    TYPE (File), INTENT(inout)   :: this
    CHARACTER(len=*), INTENT(in) :: str
    INTEGER                      :: err

    IF (.NOT. this%is_initialized) THEN
       error = .TRUE.
       CALL errorMessage("File / writeString", &
            "Object has not yet been initialized.", 1)
       RETURN
    END IF

    IF (.NOT. this%opened) THEN
       error = .TRUE.
       CALL errorMessage("File / writeString", &
            "File has not yet been opened.", 1)
       RETURN
    END IF

    WRITE(getUnit(this%lu), *, iostat=err) TRIM(str)
    IF (err /= 0) THEN
       error = .TRUE.
       CALL errorMessage("File / writeString", &
            "Could not write string to file.", 1)
       RETURN       
    END IF

  END SUBROUTINE writeString_unformatted





END MODULE File_cl





