!====================================================================!
!                                                                    !
! Copyright 2002-2022,2023                                           !
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
!! *Module*description*:
!!
!! Parameters for modules.
!!
!! @author  MG
!! @version 2023-04-05
!!
MODULE parameters

  IMPLICIT NONE
  INTEGER, PARAMETER :: iprec1 = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: iprec4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: iprec8 = SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: iprec16 = SELECTED_INT_KIND(36)
  INTEGER, PARAMETER :: rprec4 = SELECTED_REAL_KIND(p=6)
  INTEGER, PARAMETER :: rprec8 = SELECTED_REAL_KIND(p=15)
  INTEGER, PARAMETER :: rprec16 = SELECTED_REAL_KIND(p=15)

  INTEGER, PARAMETER :: FNAME_LEN = 512

END MODULE parameters
