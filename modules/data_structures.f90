!====================================================================!
!                                                                    !
! Copyright 2002-2017,2018                                           !
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
!! Tools for building and manipulating data strucures such as trees
!! and lists.
!!
!! @author  MG
!! @version 2008-06-06
!!
MODULE data_structures

  USE parameters
  USE utilities
  IMPLICIT NONE
  INTEGER(iprec4), PARAMETER :: r16arr_size_max = 1270 !1270
  INTEGER(iprec4), PARAMETER :: i8arr_size_max = 1270 !1270

  TYPE rb_tree_ch32
     TYPE (rb_tree_node_ch32), POINTER :: root => NULL()
     TYPE (rb_tree_node_ch32), POINTER :: nil => NULL()     
  END TYPE rb_tree_ch32

  TYPE rb_tree_ch32_8r8
     TYPE (rb_tree_node_ch32_8r8), POINTER :: root => NULL()
     TYPE (rb_tree_node_ch32_8r8), POINTER :: nil => NULL()     
  END TYPE rb_tree_ch32_8r8

  TYPE rb_tree_r16_i4
     TYPE (rb_tree_node_r16_i4), POINTER :: root => NULL()
     TYPE (rb_tree_node_r16_i4), POINTER :: nil => NULL()     
  END TYPE rb_tree_r16_i4

  TYPE rb_tree_r16_i4arr
     TYPE (rb_tree_node_r16_i4arr), POINTER :: root => NULL()
     TYPE (rb_tree_node_r16_i4arr), POINTER :: nil => NULL()     
  END TYPE rb_tree_r16_i4arr

  TYPE rb_tree_i8_i4arr
     TYPE (rb_tree_node_i8_i4arr), POINTER :: root => NULL()
     TYPE (rb_tree_node_i8_i4arr), POINTER :: nil => NULL()     
  END TYPE rb_tree_i8_i4arr

  TYPE rb_tree_i8_i4
     TYPE (rb_tree_node_i8_i4), POINTER :: root => NULL()
     TYPE (rb_tree_node_i8_i4), POINTER :: nil => NULL()     
  END TYPE rb_tree_i8_i4

  TYPE rb_tree_i8_ch16arr
     TYPE (rb_tree_node_i8_ch16arr), POINTER :: root => NULL()
     TYPE (rb_tree_node_i8_ch16arr), POINTER :: nil => NULL()     
  END TYPE rb_tree_i8_ch16arr

  TYPE rb_tree_node_i8_ch16arr
     TYPE (rb_tree_node_i8_ch16arr), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_i8_ch16arr), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_i8_ch16arr), POINTER :: parent => NULL()
     LOGICAL :: red
     INTEGER(iprec8) :: key
     CHARACTER(len=16), DIMENSION(:), POINTER :: data_array
  END TYPE rb_tree_node_i8_ch16arr

  TYPE rb_tree_node_i8_i4
     TYPE (rb_tree_node_i8_i4), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_i8_i4), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_i8_i4), POINTER :: parent => NULL()
     LOGICAL :: red
     INTEGER(iprec8) :: key
     TYPE (list_node_i4), POINTER :: data_list => NULL()
  END TYPE rb_tree_node_i8_i4

  TYPE rb_tree_node_r16_i4arr
     TYPE (rb_tree_node_r16_i4arr), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_r16_i4arr), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_r16_i4arr), POINTER :: parent => NULL()
     LOGICAL :: red
     REAL(rprec16) :: lkey
     REAL(rprec16) :: ukey
     !INTEGER(iprec1), DIMENSION(:), POINTER :: indx
     TYPE (r16_i4arr_node), DIMENSION(:), POINTER :: data_nodes
  END TYPE rb_tree_node_r16_i4arr

  TYPE r16_i4arr_node
     REAL(rprec16) :: key
     INTEGER(iprec4), DIMENSION(:), POINTER :: data_array
  END TYPE r16_i4arr_node

  TYPE rb_tree_node_i8_i4arr
     TYPE (rb_tree_node_i8_i4arr), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_i8_i4arr), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_i8_i4arr), POINTER :: parent => NULL()
     LOGICAL :: red
     INTEGER(iprec8) :: lkey
     INTEGER(iprec8) :: ukey
     !INTEGER(iprec1), DIMENSION(:), POINTER :: indx
     TYPE (i8_i4arr_node), DIMENSION(:), POINTER :: data_nodes
  END TYPE rb_tree_node_i8_i4arr

  TYPE i8_i4arr_node
     INTEGER(iprec8) :: key
     INTEGER(iprec4), DIMENSION(:), POINTER :: data_array
  END TYPE i8_i4arr_node

  TYPE rb_tree_node_r16_i4
     TYPE (rb_tree_node_r16_i4), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_r16_i4), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_r16_i4), POINTER :: parent => NULL()
     LOGICAL :: red
     REAL(rprec16) :: key
     TYPE (list_node_i4), POINTER :: data_list => NULL()
  END TYPE rb_tree_node_r16_i4

  TYPE rb_tree_node_ch32
     TYPE (rb_tree_node_ch32), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_ch32), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_ch32), POINTER :: parent => NULL()
     LOGICAL :: red
     CHARACTER(len=32) :: key
  END TYPE rb_tree_node_ch32

  TYPE rb_tree_node_ch32_8r8
     TYPE (rb_tree_node_ch32_8r8), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_ch32_8r8), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_ch32_8r8), POINTER :: parent => NULL()
     LOGICAL :: red
     CHARACTER(len=32) :: key
     REAL(rprec8), DIMENSION(8) :: r8arr
     INTEGER(iprec4), DIMENSION(3) :: i4arr
  END TYPE rb_tree_node_ch32_8r8

  TYPE rb_tree_i8
     TYPE (rb_tree_node_i8), POINTER :: root => NULL()
     TYPE (rb_tree_node_i8), POINTER :: nil => NULL()     
  END TYPE rb_tree_i8

  TYPE rb_tree_node_i8
     TYPE (rb_tree_node_i8), POINTER :: lchild => NULL()
     TYPE (rb_tree_node_i8), POINTER :: rchild => NULL()
     TYPE (rb_tree_node_i8), POINTER :: parent => NULL()
     LOGICAL :: red
     INTEGER(iprec8) :: key
  END TYPE rb_tree_node_i8

  TYPE binary_tree_i8
     TYPE (binary_tree_node_i8), POINTER :: root => NULL()
     TYPE (binary_tree_node_i8), POINTER :: nil => NULL()     
  END TYPE binary_tree_i8

  TYPE binary_tree_node_i8
     TYPE (binary_tree_node_i8), POINTER :: lchild => NULL()
     TYPE (binary_tree_node_i8), POINTER :: rchild => NULL()
     TYPE (binary_tree_node_i8), POINTER :: parent => NULL()
     INTEGER(iprec8) :: key
  END TYPE binary_tree_node_i8

  TYPE list_i4
     TYPE (list_node_i4), POINTER :: nil => NULL()
  END TYPE list_i4

  TYPE list_node_i4
     INTEGER(iprec4) :: key
     TYPE (list_node_i4), POINTER :: next => NULL()
  END TYPE list_node_i4

  TYPE list_i8
     INTEGER(iprec8) :: key
     TYPE (list_i8), POINTER :: next => NULL()
  END TYPE list_i8

  TYPE dbl_lnkd_list_i4
     TYPE (dbl_lnkd_list_node_i4), POINTER :: nil => NULL()
  END TYPE dbl_lnkd_list_i4

  TYPE dbl_lnkd_list_i8
     TYPE (dbl_lnkd_list_node_i8), POINTER :: nil => NULL()
  END TYPE dbl_lnkd_list_i8

  TYPE dbl_lnkd_list_node_i4
     INTEGER(iprec4) :: key
     TYPE (dbl_lnkd_list_node_i4), POINTER :: previous => NULL()
     TYPE (dbl_lnkd_list_node_i4), POINTER :: next => NULL()
  END TYPE dbl_lnkd_list_node_i4

  TYPE dbl_lnkd_list_node_i8
     INTEGER(iprec8) :: key
     TYPE (dbl_lnkd_list_node_i8), POINTER :: previous => NULL()
     TYPE (dbl_lnkd_list_node_i8), POINTER :: next => NULL()
  END TYPE dbl_lnkd_list_node_i8

  TYPE lnkd_list_ch1024_i4
     TYPE (lnkd_list_node_ch1024_i4), POINTER :: nil => NULL()
     TYPE (lnkd_list_node_ch1024_i4), POINTER :: head => NULL()
     TYPE (lnkd_list_node_ch1024_i4), POINTER :: tail => NULL()
  END TYPE lnkd_list_ch1024_i4

  TYPE lnkd_list_node_ch1024_i4
     CHARACTER(len=1024) :: str
     INTEGER(iprec4) :: key
     TYPE (lnkd_list_node_ch1024_i4), POINTER :: next => NULL()
  END TYPE lnkd_list_node_ch1024_i4

  TYPE dbl_lnkd_list_ch1024_i4
     TYPE (dbl_lnkd_list_node_ch1024_i4), POINTER :: nil => NULL()
     TYPE (dbl_lnkd_list_node_ch1024_i4), POINTER :: head => NULL()
     TYPE (dbl_lnkd_list_node_ch1024_i4), POINTER :: tail => NULL()
  END TYPE dbl_lnkd_list_ch1024_i4

  TYPE dbl_lnkd_list_node_ch1024_i4
     CHARACTER(len=1024) :: str
     INTEGER(iprec4) :: key
     TYPE (dbl_lnkd_list_node_ch1024_i4), POINTER :: previous => NULL()
     TYPE (dbl_lnkd_list_node_ch1024_i4), POINTER :: next => NULL()
  END TYPE dbl_lnkd_list_node_ch1024_i4

  INTERFACE init
     MODULE PROCEDURE init_list_i4
     MODULE PROCEDURE init_dbl_lnkd_list_i4
     MODULE PROCEDURE init_dbl_lnkd_list_i8
     MODULE PROCEDURE init_lnkd_list_ch1024_i4
     MODULE PROCEDURE init_rb_tree_ch32
     MODULE PROCEDURE init_rb_tree_ch32_8r8
     MODULE PROCEDURE init_rb_tree_r16_i4
     MODULE PROCEDURE init_rb_tree_i8_i4
     MODULE PROCEDURE init_rb_tree_i8_ch16arr
     MODULE PROCEDURE init_rb_tree_r16_i4arr
     MODULE PROCEDURE init_rb_tree_i8_i4arr
     MODULE PROCEDURE init_rb_tree_i8
     MODULE PROCEDURE init_binary_tree_i8
  END INTERFACE init

  INTERFACE insert_list_node
     MODULE PROCEDURE insert_list_node_i4_list
     MODULE PROCEDURE insert_list_node_i4_head
     MODULE PROCEDURE insert_list_node_ch1024_i4
     MODULE PROCEDURE insert_list_node_dbl_i4
     MODULE PROCEDURE insert_list_node_dbl_i8
  END INTERFACE insert_list_node

  INTERFACE delete_list_node
     MODULE PROCEDURE delete_list_node_i4_list
     MODULE PROCEDURE delete_list_node_i4_head
     MODULE PROCEDURE delete_list_node_dbl_i4
     MODULE PROCEDURE delete_list_node_dbl_i8
     MODULE PROCEDURE delete_list_node_ch1024_i4
  END INTERFACE delete_list_node

  INTERFACE delete_list
     MODULE PROCEDURE delete_list_lnkd_i4
     MODULE PROCEDURE delete_list_dbl_lnkd_i4
     MODULE PROCEDURE delete_list_dbl_lnkd_i8
  END INTERFACE delete_list

  INTERFACE delete_tree
     MODULE PROCEDURE delete_tree_rb_ch32
     MODULE PROCEDURE delete_tree_rb_ch32_8r8
     MODULE PROCEDURE delete_tree_rb_r16_i4
     MODULE PROCEDURE delete_tree_rb_i8_i4
     MODULE PROCEDURE delete_tree_rb_i8_ch16arr
     MODULE PROCEDURE delete_tree_rb_r16_i4arr
     MODULE PROCEDURE delete_tree_rb_i8_i4arr
     MODULE PROCEDURE delete_tree_rb_i8
     MODULE PROCEDURE delete_tree_binary_i8
  END INTERFACE delete_tree

  INTERFACE minimum
     MODULE PROCEDURE minimum_rb_tree_node_ch32
     MODULE PROCEDURE minimum_rb_tree_node_ch32_8r8
     MODULE PROCEDURE minimum_rb_tree_node_r16_i4
     MODULE PROCEDURE minimum_rb_tree_node_i8_i4
     MODULE PROCEDURE minimum_rb_tree_node_i8_ch16arr
     MODULE PROCEDURE minimum_rb_tree_node_r16_i4arr
     MODULE PROCEDURE minimum_rb_tree_node_i8_i4arr
     MODULE PROCEDURE minimum_rb_tree_node_i8
  END INTERFACE minimum

  INTERFACE maximum
     MODULE PROCEDURE maximum_rb_tree_node_r16_i4
     MODULE PROCEDURE maximum_rb_tree_node_i8_i4
     MODULE PROCEDURE maximum_rb_tree_node_i8
  END INTERFACE maximum

  INTERFACE successor
     MODULE PROCEDURE succ_rb_tree_node_ch32
     MODULE PROCEDURE succ_rb_tree_node_ch32_8r8
     MODULE PROCEDURE succ_list_node_ch1024_i4
     MODULE PROCEDURE succ_rb_tree_node_r16_i4
     MODULE PROCEDURE succ_rb_tree_node_i8_i4
     MODULE PROCEDURE succ_rb_tree_node_i8_ch16arr
     MODULE PROCEDURE succ_rb_tree_node_r16_i4arr
     MODULE PROCEDURE succ_rb_tree_node_i8_i4arr
     MODULE PROCEDURE succ_rb_tree_node_i8
  END INTERFACE successor

  INTERFACE preorder_tree_walk
     MODULE PROCEDURE preorder_tree_walk_rb_i8_i4
     MODULE PROCEDURE preorder_tree_walk_rb_i8
  END INTERFACE preorder_tree_walk

  INTERFACE inorder_tree_walk
     MODULE PROCEDURE inorder_tree_walk_rb_i8_i4
     MODULE PROCEDURE inorder_tree_walk_rb_i8
  END INTERFACE inorder_tree_walk

  INTERFACE postorder_tree_walk
     MODULE PROCEDURE postorder_tree_walk_rb_i8_i4
     MODULE PROCEDURE postorder_tree_walk_rb_i8
  END INTERFACE postorder_tree_walk

  INTERFACE leftrotate
     MODULE PROCEDURE leftrotate_rb_tree_ch32
     MODULE PROCEDURE leftrotate_rb_tree_ch32_8r8
     MODULE PROCEDURE leftrotate_rb_tree_r16_i4
     MODULE PROCEDURE leftrotate_rb_tree_i8_i4
     MODULE PROCEDURE leftrotate_rb_tree_i8_ch16arr
     MODULE PROCEDURE leftrotate_rb_tree_r16_i4arr
     MODULE PROCEDURE leftrotate_rb_tree_i8_i4arr
     MODULE PROCEDURE leftrotate_rb_tree_i8
  END INTERFACE leftrotate

  INTERFACE rightrotate
     MODULE PROCEDURE rightrotate_rb_tree_ch32
     MODULE PROCEDURE rightrotate_rb_tree_ch32_8r8
     MODULE PROCEDURE rightrotate_rb_tree_r16_i4
     MODULE PROCEDURE rightrotate_rb_tree_i8_i4
     MODULE PROCEDURE rightrotate_rb_tree_i8_ch16arr
     MODULE PROCEDURE rightrotate_rb_tree_r16_i4arr
     MODULE PROCEDURE rightrotate_rb_tree_i8_i4arr
     MODULE PROCEDURE rightrotate_rb_tree_i8
  END INTERFACE rightrotate

  INTERFACE insert_tree_node
     MODULE PROCEDURE insert_tree_node_rb_ch32
     MODULE PROCEDURE insert_tree_node_rb_ch32_8r8
     MODULE PROCEDURE insert_tree_node_rb_r16_i4
     MODULE PROCEDURE insert_tree_node_rb_i8_i4
     MODULE PROCEDURE insert_tree_node_rb_i8_ch16arr
     MODULE PROCEDURE insert_tree_node_rb_r16_i4arr
     MODULE PROCEDURE insert_tree_node_rb_i8_i4arr
     MODULE PROCEDURE insert_tree_node_rb_i8
     MODULE PROCEDURE insert_tree_node_binary_i8
  END INTERFACE insert_tree_node

  INTERFACE rb_insert_fixup
     MODULE PROCEDURE rb_insert_fixup_ch32
     MODULE PROCEDURE rb_insert_fixup_ch32_8r8
     MODULE PROCEDURE rb_insert_fixup_r16_i4
     MODULE PROCEDURE rb_insert_fixup_i8_i4
     MODULE PROCEDURE rb_insert_fixup_i8_ch16arr
     MODULE PROCEDURE rb_insert_fixup_r16_i4arr
     MODULE PROCEDURE rb_insert_fixup_i8_i4arr
     MODULE PROCEDURE rb_insert_fixup_i8
  END INTERFACE rb_insert_fixup

  INTERFACE delete_tree_node
     MODULE PROCEDURE delete_tree_node_rb_ch32
     MODULE PROCEDURE delete_tree_node_rb_ch32_8r8
     MODULE PROCEDURE delete_tree_node_rb_r16_i4
     MODULE PROCEDURE delete_tree_node_rb_i8_i4
     MODULE PROCEDURE delete_tree_node_rb_i8_ch16arr
     MODULE PROCEDURE delete_tree_node_rb_r16_i4arr
     MODULE PROCEDURE delete_tree_node_rb_i8_i4arr
     MODULE PROCEDURE delete_tree_node_rb_i8
  END INTERFACE delete_tree_node

  INTERFACE rb_delete_fixup
     MODULE PROCEDURE rb_delete_fixup_ch32
     MODULE PROCEDURE rb_delete_fixup_ch32_8r8
     MODULE PROCEDURE rb_delete_fixup_r16_i4
     MODULE PROCEDURE rb_delete_fixup_i8_i4
     MODULE PROCEDURE rb_delete_fixup_i8_ch16arr
     MODULE PROCEDURE rb_delete_fixup_r16_i4arr
     MODULE PROCEDURE rb_delete_fixup_i8_i4arr
     MODULE PROCEDURE rb_delete_fixup_i8
  END INTERFACE rb_delete_fixup

  INTERFACE search
     MODULE PROCEDURE search_dbl_lnkd_list_i4
     MODULE PROCEDURE search_dbl_lnkd_list_i8
     MODULE PROCEDURE search_lnkd_list_ch1024_i4
     MODULE PROCEDURE search_rb_tree_ch32
     MODULE PROCEDURE search_rb_tree_ch32_8r8
     MODULE PROCEDURE search_rb_tree_r16_i4
     MODULE PROCEDURE search_rb_tree_i8_i4
     MODULE PROCEDURE search_rb_tree_i8_ch16arr
     !MODULE PROCEDURE search_rb_tree_r16_i4arr
     MODULE PROCEDURE search_rb_tree_i8_i4arr
     MODULE PROCEDURE search_rb_tree_i8
  END INTERFACE search

  INTERFACE reallocate
     MODULE PROCEDURE reallocate_r16_i4arr_node_1
     MODULE PROCEDURE reallocate_i8_i4arr_node_1
  END INTERFACE reallocate


CONTAINS


  SUBROUTINE init_list_i4(list)

    IMPLICIT NONE
    TYPE (list_i4), POINTER :: list

    ALLOCATE(list)
    ALLOCATE(list%nil)
    list%nil%next => list%nil

  END SUBROUTINE init_list_i4





  SUBROUTINE init_dbl_lnkd_list_i4(list)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i4), POINTER :: list

    ALLOCATE(list)
    ALLOCATE(list%nil)
    list%nil%next => list%nil
    list%nil%previous => list%nil

  END SUBROUTINE init_dbl_lnkd_list_i4





  SUBROUTINE init_dbl_lnkd_list_i8(list)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i8), POINTER :: list

    ALLOCATE(list)
    ALLOCATE(list%nil)
    list%nil%next => list%nil
    list%nil%previous => list%nil

  END SUBROUTINE init_dbl_lnkd_list_i8





  SUBROUTINE delete_list_lnkd_i4(list)

    IMPLICIT NONE
    TYPE (list_i4), POINTER :: list

    list%nil%next => NULL()
    DEALLOCATE(list%nil)
    DEALLOCATE(list)

  END SUBROUTINE delete_list_lnkd_i4





  SUBROUTINE delete_list_dbl_lnkd_i4(list)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i4), POINTER :: list

    list%nil%next => NULL()
    list%nil%previous => NULL()
    DEALLOCATE(list%nil)
    DEALLOCATE(list)

  END SUBROUTINE delete_list_dbl_lnkd_i4





  SUBROUTINE delete_list_dbl_lnkd_i8(list)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i8), POINTER :: list

    list%nil%next => NULL()
    list%nil%previous => NULL()
    DEALLOCATE(list%nil)
    DEALLOCATE(list)

  END SUBROUTINE delete_list_dbl_lnkd_i8





  SUBROUTINE init_lnkd_list_ch1024_i4(list)

    IMPLICIT NONE
    TYPE (lnkd_list_ch1024_i4), POINTER :: list

    ALLOCATE(list)
    ALLOCATE(list%nil)
    list%nil%next => list%nil
    list%head => list%nil
    list%tail => list%nil

  END SUBROUTINE init_lnkd_list_ch1024_i4





  SUBROUTINE insert_list_node_i4_list(list, key)

    IMPLICIT NONE
    TYPE (list_i4), POINTER :: list
    INTEGER(iprec4) :: key
    TYPE (list_node_i4), POINTER :: new_node, tmp_node

    ALLOCATE(new_node)
    new_node%next => list%nil%next
    list%nil%next => new_node
    new_node%key = key

  END SUBROUTINE insert_list_node_i4_list





  SUBROUTINE insert_list_node_i4_head(head, key)

    IMPLICIT NONE
    TYPE (list_node_i4), POINTER :: head
    INTEGER(iprec4) :: key
    TYPE (list_node_i4), POINTER :: new_node, tmp_node

    ALLOCATE(new_node)
    new_node%next => head
    head => new_node
    new_node%key = key

  END SUBROUTINE insert_list_node_i4_head





  SUBROUTINE insert_list_node_dbl_i4(list, key)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i4), POINTER :: list
    INTEGER(iprec4) :: key
    TYPE (dbl_lnkd_list_node_i4), POINTER :: new_node

    ALLOCATE(new_node)
    new_node%next => list%nil%previous
    list%nil%next%previous => new_node
    list%nil%next => new_node
    new_node%previous => list%nil
    new_node%key = key

  END SUBROUTINE insert_list_node_dbl_i4





  SUBROUTINE insert_list_node_dbl_i8(list, key)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i8), POINTER :: list
    INTEGER(iprec8) :: key
    TYPE (dbl_lnkd_list_node_i8), POINTER :: new_node

    ALLOCATE(new_node)
    new_node%next => list%nil%previous
    list%nil%next%previous => new_node
    list%nil%next => new_node
    new_node%previous => list%nil
    new_node%key = key

  END SUBROUTINE insert_list_node_dbl_i8





!!$  SUBROUTINE insert_list_node_ch1024_i4(list, key)
!!$
!!$    IMPLICIT NONE
!!$    TYPE (lnkd_list_ch1024_i4), POINTER :: list
!!$    INTEGER(iprec4), INTENT(in) :: key
!!$    TYPE (lnkd_list_node_ch1024_i4), POINTER :: new_node
!!$
!!$    ALLOCATE(new_node)
!!$    new_node%next => list%nil
!!$    new_node%key = key
!!$    !new_node%str = str    
!!$    IF (ASSOCIATED(list%tail,list%nil) .AND. &
!!$         ASSOCIATED(list%head,list%nil)) THEN
!!$       list%head => new_node
!!$       list%tail => new_node
!!$    ELSE
!!$       list%tail%next => new_node
!!$       list%tail => list%tail%next
!!$    END IF
!!$
!!$  END SUBROUTINE insert_list_node_ch1024_i4





  SUBROUTINE insert_list_node_ch1024_i4(list, key, str)

    IMPLICIT NONE
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    CHARACTER(len=*), INTENT(in) :: str
    INTEGER(iprec4), INTENT(in) :: key
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: new_node
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: list_node

    list_node => list%head
    DO WHILE (list_node%key /= key .AND. &
         .NOT.ASSOCIATED(list_node,list%nil))
       list_node => list_node%next       
    END DO
    IF (.NOT.ASSOCIATED(list_node,list%nil)) THEN
       list_node%str = TRIM(list_node%str) // TRIM(str)
    ELSE
       ALLOCATE(new_node)
       new_node%next => list%nil
       new_node%key = key
       new_node%str = str    
       IF (ASSOCIATED(list%tail,list%nil) .AND. &
            ASSOCIATED(list%head,list%nil)) THEN
          list%head => new_node
          list%tail => new_node
       ELSE
          list%tail%next => new_node
          list%tail => list%tail%next
       END IF
    END IF

  END SUBROUTINE insert_list_node_ch1024_i4





  SUBROUTINE delete_list_node_i4_list(list, node)

    IMPLICIT NONE
    TYPE (list_i4), POINTER :: list
    TYPE (list_node_i4), POINTER :: node
    TYPE (list_node_i4), POINTER :: node_

    node_ => list%nil
    DO WHILE (.NOT.ASSOCIATED(node,node_%next))
       node_ => node_%next
    END DO
    node_%next => node%next
    DEALLOCATE(node)

  END SUBROUTINE delete_list_node_i4_list




  SUBROUTINE delete_list_node_i4_head(head, node)

    IMPLICIT NONE
    TYPE (list_node_i4), POINTER :: head
    TYPE (list_node_i4), POINTER :: node
    TYPE (list_node_i4), POINTER :: node_

    IF (ASSOCIATED(head,node)) THEN
       head => node%next
    ELSE
       node_ => head
       DO WHILE (.NOT.ASSOCIATED(node,node_%next))
          node_ => node_%next
       END DO
       node_%next => node%next
    END IF
    DEALLOCATE(node)

  END SUBROUTINE delete_list_node_i4_head




  SUBROUTINE delete_list_node_dbl_i4(node)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_node_i4), POINTER :: node

    node%previous%next => node%next
    node%next%previous => node%previous
    DEALLOCATE(node)

  END SUBROUTINE delete_list_node_dbl_i4




  SUBROUTINE delete_list_node_dbl_i8(node)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_node_i8), POINTER :: node

    node%previous%next => node%next
    node%next%previous => node%previous
    DEALLOCATE(node)

  END SUBROUTINE delete_list_node_dbl_i8




  SUBROUTINE delete_list_node_ch1024_i4(list)

    IMPLICIT NONE
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: x

    x => list%tail
    list%tail => list%nil
    DEALLOCATE(x)

  END SUBROUTINE delete_list_node_ch1024_i4




  SUBROUTINE init_rb_tree_ch32(tree)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_ch32




  SUBROUTINE init_rb_tree_ch32_8r8(tree)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_ch32_8r8




  SUBROUTINE init_rb_tree_r16_i4(tree)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_r16_i4




  SUBROUTINE init_rb_tree_i8_i4(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_i8_i4




  SUBROUTINE init_rb_tree_i8_ch16arr(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_i8_ch16arr




  SUBROUTINE init_rb_tree_r16_i4arr(tree)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_r16_i4arr




  SUBROUTINE init_rb_tree_i8_i4arr(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_i8_i4arr




  SUBROUTINE init_rb_tree_i8(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%nil%red = .FALSE.
    tree%root => tree%nil

  END SUBROUTINE init_rb_tree_i8




  SUBROUTINE init_binary_tree_i8(tree)

    IMPLICIT NONE
    TYPE (binary_tree_i8), POINTER :: tree

    ALLOCATE(tree)
    ALLOCATE(tree%nil)
    tree%nil%lchild => tree%nil
    tree%nil%rchild => tree%nil
    tree%nil%parent => tree%nil
    tree%root => tree%nil

  END SUBROUTINE init_binary_tree_i8




  SUBROUTINE delete_tree_rb_ch32(tree)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_ch32




  SUBROUTINE delete_tree_rb_ch32_8r8(tree)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_ch32_8r8




  SUBROUTINE delete_tree_rb_r16_i4(tree)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_r16_i4




  SUBROUTINE delete_tree_rb_i8_i4(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_i8_i4




  SUBROUTINE delete_tree_rb_i8_ch16arr(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_i8_ch16arr




  SUBROUTINE delete_tree_rb_r16_i4arr(tree)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_r16_i4arr




  SUBROUTINE delete_tree_rb_i8_i4arr(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_i8_i4arr




  SUBROUTINE delete_tree_rb_i8(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_rb_i8




  SUBROUTINE delete_tree_binary_i8(tree)

    IMPLICIT NONE
    TYPE (binary_tree_i8), POINTER :: tree
    INTEGER :: err

    DEALLOCATE(tree%nil, stat=err)
    DEALLOCATE(tree, stat=err)

  END SUBROUTINE delete_tree_binary_i8




  FUNCTION minimum_rb_tree_node_ch32(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: x
    TYPE (rb_tree_node_ch32), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_ch32




  FUNCTION minimum_rb_tree_node_ch32_8r8(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x
    TYPE (rb_tree_node_ch32_8r8), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_ch32_8r8




  FUNCTION minimum_rb_tree_node_r16_i4(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (rb_tree_node_r16_i4), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_r16_i4




  FUNCTION minimum_rb_tree_node_i8_i4(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (rb_tree_node_i8_i4), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_i8_i4




  FUNCTION minimum_rb_tree_node_i8_ch16arr(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_i8_ch16arr




  FUNCTION minimum_rb_tree_node_r16_i4arr(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x
    TYPE (rb_tree_node_r16_i4arr), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_r16_i4arr




  FUNCTION minimum_rb_tree_node_i8_i4arr(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x
    TYPE (rb_tree_node_i8_i4arr), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_i8_i4arr




  FUNCTION minimum_rb_tree_node_i8(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x
    TYPE (rb_tree_node_i8), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%lchild))
       y => y%lchild
    END DO

  END FUNCTION minimum_rb_tree_node_i8




  FUNCTION maximum_rb_tree_node_r16_i4(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (rb_tree_node_r16_i4), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%rchild))
       y => y%rchild
    END DO

  END FUNCTION maximum_rb_tree_node_r16_i4




  FUNCTION maximum_rb_tree_node_i8_i4(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (rb_tree_node_i8_i4), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%rchild))
       y => y%rchild
    END DO

  END FUNCTION maximum_rb_tree_node_i8_i4




  FUNCTION maximum_rb_tree_node_i8(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x
    TYPE (rb_tree_node_i8), POINTER :: y

    y => x
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y%rchild))
       y => y%rchild
    END DO

  END FUNCTION maximum_rb_tree_node_i8




  FUNCTION succ_list_node_ch1024_i4(list, x) RESULT(y)

    IMPLICIT NONE
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: x
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: y

    y => x%next

  END FUNCTION succ_list_node_ch1024_i4




  FUNCTION succ_rb_tree_node_ch32(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: x
    TYPE (rb_tree_node_ch32), POINTER :: y
    TYPE (rb_tree_node_ch32), POINTER :: xx

    xx => x
    IF (.NOT.ASSOCIATED(tree%nil,xx%rchild)) THEN
       y => minimum(tree, xx%rchild)
       RETURN
    END IF
    y => xx%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(xx,y%rchild))
       xx => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_ch32





  FUNCTION succ_rb_tree_node_ch32_8r8(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x
    TYPE (rb_tree_node_ch32_8r8), POINTER :: y
    TYPE (rb_tree_node_ch32_8r8), POINTER :: xx

    xx => x
    IF (.NOT.ASSOCIATED(tree%nil,xx%rchild)) THEN
       y => minimum(tree, xx%rchild)
       RETURN
    END IF
    y => xx%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(xx,y%rchild))
       xx => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_ch32_8r8





  FUNCTION succ_rb_tree_node_r16_i4(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (rb_tree_node_r16_i4), POINTER :: y

    IF (.NOT.ASSOCIATED(tree%nil,x%rchild)) THEN
       y => minimum(tree, x%rchild)
       RETURN
    END IF
    y => x%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(x,y%rchild))
       x => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_r16_i4





  FUNCTION succ_rb_tree_node_i8_i4(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (rb_tree_node_i8_i4), POINTER :: y

    IF (.NOT.ASSOCIATED(tree%nil,x%rchild)) THEN
       y => minimum(tree, x%rchild)
       RETURN
    END IF
    y => x%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(x,y%rchild))
       x => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_i8_i4





  FUNCTION succ_rb_tree_node_i8_ch16arr(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: y

    IF (.NOT.ASSOCIATED(tree%nil,x%rchild)) THEN
       y => minimum(tree, x%rchild)
       RETURN
    END IF
    y => x%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(x,y%rchild))
       x => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_i8_ch16arr





  FUNCTION succ_rb_tree_node_r16_i4arr(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x
    TYPE (rb_tree_node_r16_i4arr), POINTER :: y

    IF (.NOT.ASSOCIATED(tree%nil,x%rchild)) THEN
       y => minimum(tree, x%rchild)
       RETURN
    END IF
    y => x%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(x,y%rchild))
       x => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_r16_i4arr





  FUNCTION succ_rb_tree_node_i8_i4arr(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x
    TYPE (rb_tree_node_i8_i4arr), POINTER :: y

    IF (.NOT.ASSOCIATED(tree%nil,x%rchild)) THEN
       y => minimum(tree, x%rchild)
       RETURN
    END IF
    y => x%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(x,y%rchild))
       x => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_i8_i4arr





  FUNCTION succ_rb_tree_node_i8(tree, x) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x
    TYPE (rb_tree_node_i8), POINTER :: y

    IF (.NOT.ASSOCIATED(tree%nil,x%rchild)) THEN
       y => minimum(tree, x%rchild)
       RETURN
    END IF
    y => x%parent
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y) .AND. ASSOCIATED(x,y%rchild))
       x => y
       y => y%parent
    END DO

  END FUNCTION succ_rb_tree_node_i8





  SUBROUTINE leftrotate_rb_tree_ch32(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: x
    TYPE (rb_tree_node_ch32), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_ch32





  SUBROUTINE leftrotate_rb_tree_ch32_8r8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x
    TYPE (rb_tree_node_ch32_8r8), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_ch32_8r8





  SUBROUTINE leftrotate_rb_tree_r16_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (rb_tree_node_r16_i4), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_r16_i4





  SUBROUTINE leftrotate_rb_tree_i8_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (rb_tree_node_i8_i4), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_i8_i4





  SUBROUTINE leftrotate_rb_tree_i8_ch16arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_i8_ch16arr





  SUBROUTINE leftrotate_rb_tree_r16_i4arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x
    TYPE (rb_tree_node_r16_i4arr), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_r16_i4arr





  SUBROUTINE leftrotate_rb_tree_i8_i4arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x
    TYPE (rb_tree_node_i8_i4arr), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_i8_i4arr





  SUBROUTINE leftrotate_rb_tree_i8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x
    TYPE (rb_tree_node_i8), POINTER :: y, w

    ! Set y and w
    y => x%rchild
    w => x
    ! Turn y's left subtree into x's right subtree
    x%rchild => y%lchild
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       y%lchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%lchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE leftrotate_rb_tree_i8





  SUBROUTINE rightrotate_rb_tree_ch32(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: x
    TYPE (rb_tree_node_ch32), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_ch32




  SUBROUTINE rightrotate_rb_tree_ch32_8r8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x
    TYPE (rb_tree_node_ch32_8r8), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_ch32_8r8




  SUBROUTINE rightrotate_rb_tree_r16_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (rb_tree_node_r16_i4), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_r16_i4




  SUBROUTINE rightrotate_rb_tree_i8_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (rb_tree_node_i8_i4), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_i8_i4




  SUBROUTINE rightrotate_rb_tree_i8_ch16arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_i8_ch16arr




  SUBROUTINE rightrotate_rb_tree_r16_i4arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x
    TYPE (rb_tree_node_r16_i4arr), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_r16_i4arr




  SUBROUTINE rightrotate_rb_tree_i8_i4arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x
    TYPE (rb_tree_node_i8_i4arr), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_i8_i4arr




  SUBROUTINE rightrotate_rb_tree_i8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x
    TYPE (rb_tree_node_i8), POINTER :: y, w

    ! Set y and w
    y => x%lchild
    w => x
    ! Turn y's right subtree into x's left subtree
    x%lchild => y%rchild
    IF (.NOT.ASSOCIATED(tree%nil, y%rchild)) THEN
       y%rchild%parent => x
    END IF
    ! Link x's parent to y
    y%parent => x%parent
    IF (ASSOCIATED(tree%nil,x%parent)) THEN
       tree%root => y ! root[tree] <- y
    ELSE
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          x%parent%lchild => y
       ELSE
          x%parent%rchild => y
       END IF
    END IF
    y%rchild => w
    w%parent => y

    ! Deallocate memory 
    y => NULL()
    w => NULL()

  END SUBROUTINE rightrotate_rb_tree_i8





  SUBROUTINE insert_tree_node_rb_ch32(tree, key)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    CHARACTER(len=*), INTENT(in) :: key
    TYPE (rb_tree_node_ch32), POINTER :: x, y, z

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key == x%key) THEN
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    ALLOCATE(z)
    z%key = key
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil
    z%red = .TRUE.
    CALL rb_insert_fixup(tree, z)

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_ch32





  SUBROUTINE insert_tree_node_rb_ch32_8r8(tree, key, r8arr, i4arr)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    CHARACTER(len=*), INTENT(in) :: key
    REAL(rprec8), DIMENSION(8), INTENT(in) :: r8arr
    INTEGER(iprec4), DIMENSION(3), INTENT(in) :: i4arr
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x, y, z

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key == x%key) THEN
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    ALLOCATE(z)
    z%key = key
    z%r8arr = r8arr
    z%i4arr = i4arr
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil
    z%red = .TRUE.
    CALL rb_insert_fixup(tree, z)

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_ch32_8r8





  SUBROUTINE insert_tree_node_rb_i8_i4(tree, key, data_value)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    INTEGER(iprec4), INTENT(in), OPTIONAL :: data_value
    TYPE (rb_tree_node_i8_i4), POINTER :: x, y, z

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key == x%key) THEN
          IF (PRESENT(data_value)) THEN
             CALL insert_list_node(x%data_list, data_value)
          END IF
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    ALLOCATE(z)
    z%key = key
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil
    z%red = .TRUE.
    IF (PRESENT(data_value)) THEN
       CALL insert_list_node(z%data_list, data_value)
    END IF
    CALL rb_insert_fixup(tree, z)

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_i8_i4





  SUBROUTINE insert_tree_node_rb_i8_ch16arr(tree, key, data_array)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    CHARACTER(len=16), DIMENSION(:), INTENT(in), OPTIONAL :: data_array
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x, y, z

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key == x%key) THEN
          IF (PRESENT(data_array)) THEN
             ALLOCATE(x%data_array(SIZE(data_array)))
             x%data_array(:) = data_array(:)
          END IF
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    ALLOCATE(z)
    z%key = key
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil
    z%red = .TRUE.
    IF (PRESENT(data_array)) THEN
       ALLOCATE(z%data_array(SIZE(data_array)))
       z%data_array(:) = data_array(:)
    END IF
    CALL rb_insert_fixup(tree, z)

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_i8_ch16arr





  SUBROUTINE insert_tree_node_rb_r16_i4(tree, key, data_value)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    REAL(rprec16), INTENT(in) :: key
    INTEGER(iprec4), INTENT(in), OPTIONAL :: data_value
    TYPE (rb_tree_node_r16_i4), POINTER :: x, y, z

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (ABS(key - x%key) < 100.0_rprec16*EPSILON(key)) THEN
          IF (PRESENT(data_value)) THEN
             CALL insert_list_node(x%data_list, data_value)
          END IF
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    ALLOCATE(z)
    z%key = key
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil
    z%red = .TRUE.
    IF (PRESENT(data_value)) THEN
       CALL insert_list_node(z%data_list, data_value)
    END IF
    CALL rb_insert_fixup(tree, z)

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_r16_i4





  RECURSIVE SUBROUTINE insert_tree_node_rb_r16_i4arr(tree, key, data_value)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    REAL(rprec16), INTENT(in) :: key
    INTEGER(iprec4), INTENT(in) :: data_value
    TYPE (r16_i4arr_node), DIMENSION(:), POINTER :: data_nodes
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x, y, z, w
    !INTEGER(iprec1), DIMENSION(:), POINTER :: indx
    REAL(rprec16) :: mean_key
    INTEGER :: i, j, k, nnodes, nkey

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key >= x%lkey .AND. key <= x%ukey) THEN
          DO i=1,SIZE(x%data_nodes)
             IF (key == x%data_nodes(i)%key) THEN
                ! Found key in tree. Insert data value if it's not
                ! the latest inserted item:
                IF (x%data_nodes(i)%data_array(SIZE(x%data_nodes(i)%data_array)) /= data_value) THEN
                   x%data_nodes(i)%data_array => &
                        reallocate(x%data_nodes(i)%data_array, &
                        SIZE(x%data_nodes(i)%data_array) + 1)
                   x%data_nodes(i)%data_array(SIZE(x%data_nodes(i)%data_array)) = data_value
                END IF
                RETURN
             END IF
          END DO
          IF (SIZE(x%data_nodes) < r16arr_size_max) THEN
             nnodes = SIZE(x%data_nodes) + 1
             x%data_nodes => reallocate(x%data_nodes, nnodes)
!!$             x%indx => reallocate(x%indx, nnodes)
!!$             IF (key < x%data_nodes(x%indx(1))%key) THEN
!!$                x%indx(2:nnodes) = x%indx(1:nnodes-1)
!!$                x%indx(1) = nnodes
!!$             ELSE IF (key > x%data_nodes(x%indx(nnodes-1))%key) THEN
!!$                x%indx(nnodes) = nnodes                
!!$             ELSE
!!$                DO i=1,nnodes-1
!!$                   IF (x%data_nodes(x%indx(i))%key < key .AND. &
!!$                        key < x%data_nodes(x%indx(i+1))%key) THEN
!!$                      x%indx(i+2:nnodes) = x%indx(i+1:nnodes-1)
!!$                      x%indx(i+1) = nnodes
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$             END IF
             x%data_nodes(nnodes)%key = key
             ALLOCATE(x%data_nodes(nnodes)%data_array(1))
             x%data_nodes(nnodes)%data_array(1) = data_value
          ELSE
             mean_key = (x%data_nodes(1)%key + &
                  x%data_nodes(r16arr_size_max/4)%key + &
                  x%data_nodes(r16arr_size_max/2)%key + &
                  x%data_nodes(3*(r16arr_size_max/4))%key + &
                  x%data_nodes(r16arr_size_max)%key)/5
             nkey = 0
             DO i=1,r16arr_size_max
                IF (x%data_nodes(i)%key <= mean_key) THEN
                   nkey = nkey + 1
                END IF
             END DO
             ALLOCATE(z)
             ALLOCATE(z%data_nodes(nkey), &
                                !z%indx(r16arr_size_max/2), &
                                !indx(r16arr_size_max - r16arr_size_max/2), &
                  data_nodes(r16arr_size_max - nkey))
!!$             ! Copy first half to new node
!!$             DO i=1,r16arr_size_max/2
!!$                z%indx(i) = INT(i,iprec1)
!!$                z%data_nodes(i)%key = x%data_nodes(x%indx(i))%key
!!$                z%data_nodes(i)%data_array => x%data_nodes(x%indx(i))%data_array
!!$             END DO
!!$             ! Copy rest to temporary arrays 
!!$             DO i=1,r16arr_size_max - r16arr_size_max/2
!!$                indx(i) = INT(i,iprec1)
!!$                data_nodes(i)%key = x%data_nodes(x%indx(r16arr_size_max/2 + i))%key
!!$                data_nodes(i)%data_array => x%data_nodes(x%indx(r16arr_size_max/2 + i))%data_array
!!$             END DO
             ! Copy data from one set of nodes to two new sets of nodes
             j = 0
             k = 0
             DO i=1,r16arr_size_max
                IF (x%data_nodes(i)%key <= mean_key) THEN
                   j = j + 1
                   z%data_nodes(j)%key = x%data_nodes(i)%key
                   z%data_nodes(j)%data_array => x%data_nodes(i)%data_array
                ELSE
                   k = k + 1
                   data_nodes(k)%key = x%data_nodes(i)%key
                   data_nodes(k)%data_array => x%data_nodes(i)%data_array
                END IF
             END DO
             ! Deallocate arrays of the current node 
             DEALLOCATE(x%data_nodes)!, x%indx)
             ! Initialize the arrays of the current node with the temporary arrays
             x%data_nodes => data_nodes
             !x%indx => indx
             ! Adjust lower and upper keys to match the new structure
             z%lkey = x%lkey
!!$             x%lkey = x%data_nodes(1)%key
!!$             z%ukey = x%lkey - 1_iprec8
             x%lkey = mean_key + 1_iprec8
             z%ukey = mean_key
             ! Search for the correct position of the new node
             w => tree%root
             y => tree%nil
             DO WHILE (.NOT.ASSOCIATED(tree%nil,w))
                y => w
                IF (z%lkey < w%lkey) THEN
                   w => w%lchild
                ELSE
                   w => w%rchild
                END IF
             END DO
             z%parent => y
             IF (z%lkey < y%lkey) THEN
                y%lchild => z
             ELSE
                y%rchild => z
             END IF
             z%lchild => tree%nil 
             z%rchild => tree%nil
             z%red = .TRUE.
             CALL rb_insert_fixup(tree, z)
             ! Insert the new key and its first value
             CALL insert_tree_node(tree, key, data_value)
          END IF
          RETURN
       ELSE IF (key < x%lkey) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    IF (ASSOCIATED(tree%nil,y)) THEN
       ! Make new root node
       ALLOCATE(z)
       ALLOCATE(z%data_nodes(1))!, z%indx(1))
       z%data_nodes(1)%key = key
       !z%indx(1) = 1_iprec1
       ALLOCATE(z%data_nodes(1)%data_array(1))
       z%data_nodes(1)%data_array(1) = data_value
       z%parent => y
       tree%root => z
       z%lkey = 0_iprec8
       z%ukey = HUGE(key)
       z%lchild => tree%nil 
       z%rchild => tree%nil
       z%red = .TRUE.
       CALL rb_insert_fixup(tree, z)
    ELSE
       WRITE(0,*) 'error in insert_tree_node_rb_r16_i4arr'
       WRITE(0,'(I20)') HUGE(key)
       WRITE(0,'(2(I20,1X))') tree%root%lkey, tree%root%ukey
       WRITE(0,'(4(I20,1X))') tree%root%lchild%lkey, &
            tree%root%lchild%ukey, &
            tree%root%rchild%lkey, &
            tree%root%rchild%ukey
       WRITE(0,'(8(I20,1X))') tree%root%lchild%lchild%lkey, &
            tree%root%lchild%lchild%ukey, &
            tree%root%lchild%rchild%lkey, &
            tree%root%lchild%rchild%ukey, &
            tree%root%rchild%lchild%lkey, &
            tree%root%rchild%lchild%ukey, &
            tree%root%rchild%rchild%lkey, &
            tree%root%rchild%rchild%ukey
       STOP
    END IF

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_r16_i4arr





  RECURSIVE SUBROUTINE insert_tree_node_rb_i8_i4arr(tree, key, data_value)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    INTEGER(iprec4), INTENT(in) :: data_value
    TYPE (i8_i4arr_node), DIMENSION(:), POINTER :: data_nodes
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x, y, z, w
    !INTEGER(iprec1), DIMENSION(:), POINTER :: indx
    INTEGER(iprec8) :: mean_key
    INTEGER :: i, j, k, nnodes, nkey

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key >= x%lkey .AND. key <= x%ukey) THEN
          DO i=1,SIZE(x%data_nodes)
             IF (key == x%data_nodes(i)%key) THEN
                ! Found key in tree. Insert data value if it's not
                ! the latest inserted item:
                IF (x%data_nodes(i)%data_array(SIZE(x%data_nodes(i)%data_array)) /= data_value) THEN
                   x%data_nodes(i)%data_array => &
                        reallocate(x%data_nodes(i)%data_array, &
                        SIZE(x%data_nodes(i)%data_array) + 1)
                   x%data_nodes(i)%data_array(SIZE(x%data_nodes(i)%data_array)) = data_value
                END IF
                RETURN
             END IF
          END DO
          IF (SIZE(x%data_nodes) < i8arr_size_max) THEN
             nnodes = SIZE(x%data_nodes) + 1
             x%data_nodes => reallocate(x%data_nodes, nnodes)
!!$             x%indx => reallocate(x%indx, nnodes)
!!$             IF (key < x%data_nodes(x%indx(1))%key) THEN
!!$                x%indx(2:nnodes) = x%indx(1:nnodes-1)
!!$                x%indx(1) = nnodes
!!$             ELSE IF (key > x%data_nodes(x%indx(nnodes-1))%key) THEN
!!$                x%indx(nnodes) = nnodes                
!!$             ELSE
!!$                DO i=1,nnodes-1
!!$                   IF (x%data_nodes(x%indx(i))%key < key .AND. &
!!$                        key < x%data_nodes(x%indx(i+1))%key) THEN
!!$                      x%indx(i+2:nnodes) = x%indx(i+1:nnodes-1)
!!$                      x%indx(i+1) = nnodes
!!$                      EXIT
!!$                   END IF
!!$                END DO
!!$             END IF
             x%data_nodes(nnodes)%key = key
             ALLOCATE(x%data_nodes(nnodes)%data_array(1))
             x%data_nodes(nnodes)%data_array(1) = data_value
          ELSE
             mean_key = (x%data_nodes(1)%key + &
                  x%data_nodes(i8arr_size_max/4)%key + &
                  x%data_nodes(i8arr_size_max/2)%key + &
                  x%data_nodes(3*(i8arr_size_max/4))%key + &
                  x%data_nodes(i8arr_size_max)%key)/5
             nkey = 0
             DO i=1,i8arr_size_max
                IF (x%data_nodes(i)%key <= mean_key) THEN
                   nkey = nkey + 1
                END IF
             END DO
             ALLOCATE(z)
             ALLOCATE(z%data_nodes(nkey), &
                                !z%indx(i8arr_size_max/2), &
                                !indx(i8arr_size_max - i8arr_size_max/2), &
                  data_nodes(i8arr_size_max - nkey))
!!$             ! Copy first half to new node
!!$             DO i=1,i8arr_size_max/2
!!$                z%indx(i) = INT(i,iprec1)
!!$                z%data_nodes(i)%key = x%data_nodes(x%indx(i))%key
!!$                z%data_nodes(i)%data_array => x%data_nodes(x%indx(i))%data_array
!!$             END DO
!!$             ! Copy rest to temporary arrays 
!!$             DO i=1,i8arr_size_max - i8arr_size_max/2
!!$                indx(i) = INT(i,iprec1)
!!$                data_nodes(i)%key = x%data_nodes(x%indx(i8arr_size_max/2 + i))%key
!!$                data_nodes(i)%data_array => x%data_nodes(x%indx(i8arr_size_max/2 + i))%data_array
!!$             END DO
             ! Copy data from one set of nodes to two new sets of nodes
             j = 0
             k = 0
             DO i=1,i8arr_size_max
                IF (x%data_nodes(i)%key <= mean_key) THEN
                   j = j + 1
                   z%data_nodes(j)%key = x%data_nodes(i)%key
                   z%data_nodes(j)%data_array => x%data_nodes(i)%data_array
                ELSE
                   k = k + 1
                   data_nodes(k)%key = x%data_nodes(i)%key
                   data_nodes(k)%data_array => x%data_nodes(i)%data_array
                END IF
             END DO
             ! Deallocate arrays of the current node 
             DEALLOCATE(x%data_nodes)!, x%indx)
             ! Initialize the arrays of the current node with the temporary arrays
             x%data_nodes => data_nodes
             !x%indx => indx
             ! Adjust lower and upper keys to match the new structure
             z%lkey = x%lkey
!!$             x%lkey = x%data_nodes(1)%key
!!$             z%ukey = x%lkey - 1_iprec8
             x%lkey = mean_key + 1_iprec8
             z%ukey = mean_key
             ! Search for the correct position of the new node
             w => tree%root
             y => tree%nil
             DO WHILE (.NOT.ASSOCIATED(tree%nil,w))
                y => w
                IF (z%lkey < w%lkey) THEN
                   w => w%lchild
                ELSE
                   w => w%rchild
                END IF
             END DO
             z%parent => y
             IF (z%lkey < y%lkey) THEN
                y%lchild => z
             ELSE
                y%rchild => z
             END IF
             z%lchild => tree%nil 
             z%rchild => tree%nil
             z%red = .TRUE.
             CALL rb_insert_fixup(tree, z)
             ! Insert the new key and its first value
             CALL insert_tree_node(tree, key, data_value)
          END IF
          RETURN
       ELSE IF (key < x%lkey) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    IF (ASSOCIATED(tree%nil,y)) THEN
       ! Make new root node
       ALLOCATE(z)
       ALLOCATE(z%data_nodes(1))!, z%indx(1))
       z%data_nodes(1)%key = key
       !z%indx(1) = 1_iprec1
       ALLOCATE(z%data_nodes(1)%data_array(1))
       z%data_nodes(1)%data_array(1) = data_value
       z%parent => y
       tree%root => z
       z%lkey = 0_iprec8
       z%ukey = HUGE(key)
       z%lchild => tree%nil 
       z%rchild => tree%nil
       z%red = .TRUE.
       CALL rb_insert_fixup(tree, z)
    ELSE
       WRITE(0,*) 'error in insert_tree_node_rb_i8_i4arr'
       WRITE(0,'(I20)') HUGE(key)
       WRITE(0,'(2(I20,1X))') tree%root%lkey, tree%root%ukey
       WRITE(0,'(4(I20,1X))') tree%root%lchild%lkey, &
            tree%root%lchild%ukey, &
            tree%root%rchild%lkey, &
            tree%root%rchild%ukey
       WRITE(0,'(8(I20,1X))') tree%root%lchild%lchild%lkey, &
            tree%root%lchild%lchild%ukey, &
            tree%root%lchild%rchild%lkey, &
            tree%root%lchild%rchild%ukey, &
            tree%root%rchild%lchild%lkey, &
            tree%root%rchild%lchild%ukey, &
            tree%root%rchild%rchild%lkey, &
            tree%root%rchild%rchild%ukey
       STOP
    END IF

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_i8_i4arr





  SUBROUTINE insert_tree_node_rb_i8(tree, key)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (rb_tree_node_i8), POINTER :: x, y, z

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (key == x%key) THEN
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    ALLOCATE(z)
    z%key = key
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil
    z%red = .TRUE.
    CALL rb_insert_fixup(tree, z)

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_rb_i8





  SUBROUTINE insert_tree_node_binary_i8(tree, key)

    IMPLICIT NONE
    TYPE (binary_tree_i8), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (binary_tree_node_i8), POINTER :: x, y, z

    ALLOCATE(z)
    z%key = key

    x => tree%root
    y => tree%nil
    DO WHILE (.NOT.ASSOCIATED(tree%nil,x))
       y => x
       IF (z%key == x%key) THEN
          RETURN
       ELSE IF (key < x%key) THEN
          x => x%lchild
       ELSE
          x => x%rchild
       END IF
    END DO
    z%parent => y
    IF (ASSOCIATED(tree%nil,y)) THEN
       tree%root => z
    ELSE
       IF (z%key < y%key) THEN
          y%lchild => z
       ELSE
          y%rchild => z
       END IF
    END IF
    z%lchild => tree%nil 
    z%rchild => tree%nil

    ! Deallocate memory 
    x => NULL()
    y => NULL()
    z => NULL()

  END SUBROUTINE insert_tree_node_binary_i8





  SUBROUTINE rb_insert_fixup_ch32(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: z
    TYPE (rb_tree_node_ch32), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_ch32





  SUBROUTINE rb_insert_fixup_ch32_8r8(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: z
    TYPE (rb_tree_node_ch32_8r8), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_ch32_8r8





  SUBROUTINE rb_insert_fixup_r16_i4(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: z
    TYPE (rb_tree_node_r16_i4), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_r16_i4





  SUBROUTINE rb_insert_fixup_i8_i4(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: z
    TYPE (rb_tree_node_i8_i4), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_i8_i4





  SUBROUTINE rb_insert_fixup_i8_ch16arr(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: z
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_i8_ch16arr





  SUBROUTINE rb_insert_fixup_r16_i4arr(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: z
    TYPE (rb_tree_node_r16_i4arr), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_r16_i4arr





  SUBROUTINE rb_insert_fixup_i8_i4arr(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: z
    TYPE (rb_tree_node_i8_i4arr), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_i8_i4arr





  SUBROUTINE rb_insert_fixup_i8(tree, z)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: z
    TYPE (rb_tree_node_i8), POINTER :: y, w

    DO WHILE (z%parent%red)
       IF (ASSOCIATED(z%parent,z%parent%parent%lchild)) THEN
          y => z%parent%parent%rchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%rchild)) THEN
                z => z%parent
                CALL leftrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL rightrotate(tree, w)
          END IF
       ELSE
          y => z%parent%parent%lchild
          IF (y%red) THEN
             z%parent%red = .FALSE.
             y%red = .FALSE.
             z%parent%parent%red = .TRUE.
             z => z%parent%parent
          ELSE
             IF (ASSOCIATED(z,z%parent%lchild)) THEN
                z => z%parent
                CALL rightrotate(tree, z)
             END IF
             z%parent%red = .FALSE.
             z%parent%parent%red = .TRUE.
             w => z%parent%parent
             CALL leftrotate(tree, w)
          END IF
       END IF

       ! Deallocate memory 
       y => NULL()
       w => NULL()

    END DO
    tree%root%red = .FALSE.

  END SUBROUTINE rb_insert_fixup_i8





  FUNCTION delete_tree_node_rb_ch32(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: z, y
    TYPE (rb_tree_node_ch32), POINTER :: x

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%key = y%key
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_ch32




  FUNCTION delete_tree_node_rb_ch32_8r8(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: z, y
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%key = y%key
       z%r8arr = y%r8arr
       z%i4arr = y%i4arr
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_ch32_8r8




  FUNCTION delete_tree_node_rb_r16_i4(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: z, y
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (list_node_i4), POINTER :: w

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%key = y%key
       DO
          w => z%data_list%next
          IF (ASSOCIATED(w,z%data_list)) THEN
             EXIT
          END IF
          CALL delete_list_node(z%data_list,w)
       END DO
       DEALLOCATE(z%data_list)
       z%data_list => y%data_list
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_r16_i4




  FUNCTION delete_tree_node_rb_i8_i4(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: z, y
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (list_node_i4), POINTER :: w

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%key = y%key
       DO
          w => z%data_list%next
          IF (ASSOCIATED(w,z%data_list)) THEN
             EXIT
          END IF
          CALL delete_list_node(z%data_list,w)
       END DO
       DEALLOCATE(z%data_list)
       z%data_list => y%data_list
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_i8_i4




  FUNCTION delete_tree_node_rb_i8_ch16arr(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: z, y
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x
    INTEGER :: err

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%key = y%key
       IF (ASSOCIATED(z%data_array)) THEN
          DEALLOCATE(z%data_array)
       END IF
       z%data_array => y%data_array
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_i8_ch16arr




  FUNCTION delete_tree_node_rb_r16_i4arr(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: z, y
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x
    INTEGER :: i

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%lkey = y%lkey
       z%ukey = y%ukey
       IF (ASSOCIATED(z%data_nodes)) THEN
          DO i=1,SIZE(z%data_nodes)
             IF (ASSOCIATED(z%data_nodes(i)%data_array)) THEN
                DEALLOCATE(z%data_nodes(i)%data_array)
             END IF
          END DO
          DEALLOCATE(z%data_nodes)
       END IF
       z%data_nodes => y%data_nodes
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_r16_i4arr




  FUNCTION delete_tree_node_rb_i8_i4arr(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: z, y
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x
    INTEGER :: i

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%lkey = y%lkey
       z%ukey = y%ukey
       IF (ASSOCIATED(z%data_nodes)) THEN
          DO i=1,SIZE(z%data_nodes)
             IF (ASSOCIATED(z%data_nodes(i)%data_array)) THEN
                DEALLOCATE(z%data_nodes(i)%data_array)
             END IF
          END DO
          DEALLOCATE(z%data_nodes)
       END IF
       z%data_nodes => y%data_nodes
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_i8_i4arr




  FUNCTION delete_tree_node_rb_i8(tree, z) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: z, y
    TYPE (rb_tree_node_i8), POINTER :: x

    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
       y => z
    ELSE
       y => successor(tree, z)
    END IF
    IF (.NOT.ASSOCIATED(tree%nil,y%lchild)) THEN
       x => y%lchild
    ELSE
       x => y%rchild
    END IF
    x%parent => y%parent
    IF (ASSOCIATED(tree%nil,y%parent)) THEN
       tree%root => x
    ELSE
       IF (ASSOCIATED(y,y%parent%lchild)) THEN
          y%parent%lchild => x
       ELSE
          y%parent%rchild => x
       END IF
    END IF
    IF (.NOT.ASSOCIATED(y,z)) THEN
       z%key = y%key
    END IF
    IF (.NOT.y%red) THEN
       CALL rb_delete_fixup(tree, x)
    END IF

    ! Deallocate memory 
    x => NULL()

  END FUNCTION delete_tree_node_rb_i8




  !! *Description*:
  !!
  !! This algorithm was presented during the "Tietorakenteet" course
  !! organized in Spring 2006 at the Department of Computer Science
  !! at the University of Helsinki. The results are identical with the
  !! version by Cormen et al. (Introduction to algorithms, 2003, MIT
  !! Press).
  !! 
!!$  FUNCTION delete_tree_node_rb_i8(tree, z) RESULT(y)
!!$
!!$    IMPLICIT NONE
!!$    TYPE (rb_tree_i8), POINTER :: tree
!!$    TYPE (rb_tree_node_i8), POINTER :: z, y
!!$    TYPE (rb_tree_node_i8), POINTER :: x, w
!!$
!!$    IF (ASSOCIATED(z%lchild,tree%nil) .OR. ASSOCIATED(z%rchild,tree%nil)) THEN
!!$       if (.not.associated(z%lchild,tree%nil)) then
!!$          x => z%lchild
!!$       else
!!$          x => z%rchild
!!$       end if
!!$       w => z%parent
!!$       if (associated(w,tree%nil)) then
!!$          tree%root => x
!!$       else
!!$          if (associated(z,w%lchild)) then
!!$             w%lchild => x
!!$          else
!!$             w%rchild => x
!!$          end if
!!$       end if
!!$       x%parent => w
!!$       if (.not.z%red) then
!!$          CALL rb_delete_fixup(tree, x)
!!$       end if
!!$       y => z
!!$    else
!!$       y => successor(tree, z)
!!$       x => y%rchild
!!$       w => y%parent
!!$       IF (ASSOCIATED(y,w%lchild)) THEN
!!$          w%lchild => x
!!$       else
!!$          w%rchild => x
!!$       end IF
!!$       x%parent => w
!!$       z%key = y%key
!!$       z%data_list => y%data_list
!!$       IF (.NOT.y%red) THEN
!!$          CALL rb_delete_fixup(tree, x)
!!$       END IF
!!$    end IF
!!$
!!$  END FUNCTION delete_tree_node_rb_i8





  SUBROUTINE rb_delete_fixup_ch32(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    TYPE (rb_tree_node_ch32), POINTER :: x
    TYPE (rb_tree_node_ch32), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_ch32





  SUBROUTINE rb_delete_fixup_ch32_8r8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    TYPE (rb_tree_node_ch32_8r8), POINTER :: x
    TYPE (rb_tree_node_ch32_8r8), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_ch32_8r8





  SUBROUTINE rb_delete_fixup_r16_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    TYPE (rb_tree_node_r16_i4), POINTER :: x
    TYPE (rb_tree_node_r16_i4), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_r16_i4





  SUBROUTINE rb_delete_fixup_i8_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    TYPE (rb_tree_node_i8_i4), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_i8_i4





  SUBROUTINE rb_delete_fixup_i8_ch16arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: x
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_i8_ch16arr





  SUBROUTINE rb_delete_fixup_r16_i4arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    TYPE (rb_tree_node_r16_i4arr), POINTER :: x
    TYPE (rb_tree_node_r16_i4arr), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_r16_i4arr





  SUBROUTINE rb_delete_fixup_i8_i4arr(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    TYPE (rb_tree_node_i8_i4arr), POINTER :: x
    TYPE (rb_tree_node_i8_i4arr), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_i8_i4arr





  SUBROUTINE rb_delete_fixup_i8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x
    TYPE (rb_tree_node_i8), POINTER :: w

    DO WHILE (.NOT.ASSOCIATED(x,tree%root) .AND. .NOT.x%red)
       IF (ASSOCIATED(x,x%parent%lchild)) THEN
          w => x%parent%rchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL leftrotate(tree,x%parent)
             w => x%parent%rchild
          END IF
          IF (.NOT.w%lchild%red .AND. .NOT.w%rchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%rchild%red) THEN
                w%lchild%red = .FALSE.
                w%red = .TRUE.
                CALL rightrotate(tree, w)
                w => x%parent%rchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%rchild%red = .FALSE.
             CALL leftrotate(tree, x%parent)
             x => tree%root
          END IF
       ELSE
          w => x%parent%lchild
          IF (w%red) THEN
             w%red = .FALSE.
             x%parent%red = .TRUE.
             CALL rightrotate(tree,x%parent)
             w => x%parent%lchild
          END IF
          IF (.NOT.w%rchild%red .AND. .NOT.w%lchild%red) THEN
             w%red = .TRUE.
             x => x%parent
          ELSE
             IF (.NOT.w%lchild%red) THEN
                w%rchild%red = .FALSE.
                w%red = .TRUE.
                CALL leftrotate(tree, w)
                w => x%parent%lchild
             END IF
             w%red = x%parent%red
             x%parent%red = .FALSE.
             w%lchild%red = .FALSE.
             CALL rightrotate(tree, x%parent)
             x => tree%root
          END IF
       END IF
       ! Deallocate memory 
       w => NULL()
    END DO
    x%red = .FALSE.

  END SUBROUTINE rb_delete_fixup_i8





  RECURSIVE SUBROUTINE preorder_tree_walk_rb_i8_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x

    IF (.NOT.ASSOCIATED(tree%nil,x)) THEN
       WRITE(*,*) x%key
       CALL preorder_tree_walk(tree, x%lchild)
       CALL preorder_tree_walk(tree, x%rchild)
    END IF

  END SUBROUTINE preorder_tree_walk_rb_i8_i4




  RECURSIVE SUBROUTINE preorder_tree_walk_rb_i8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x

    IF (.NOT.ASSOCIATED(tree%nil,x)) THEN
       WRITE(*,*) x%key
       CALL preorder_tree_walk(tree, x%lchild)
       CALL preorder_tree_walk(tree, x%rchild)
    END IF

  END SUBROUTINE preorder_tree_walk_rb_i8




  RECURSIVE SUBROUTINE inorder_tree_walk_rb_i8_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x

    IF (.NOT.ASSOCIATED(tree%nil,x)) THEN
       CALL inorder_tree_walk(tree, x%lchild)
       WRITE(*,*) x%key, x%red
       CALL inorder_tree_walk(tree, x%rchild)
    END IF

  END SUBROUTINE inorder_tree_walk_rb_i8_i4




  RECURSIVE SUBROUTINE inorder_tree_walk_rb_i8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x

    IF (.NOT.ASSOCIATED(tree%nil,x)) THEN
       CALL inorder_tree_walk(tree, x%lchild)
       WRITE(*,*) x%key, x%red
       CALL inorder_tree_walk(tree, x%rchild)
    END IF

  END SUBROUTINE inorder_tree_walk_rb_i8




  RECURSIVE SUBROUTINE postorder_tree_walk_rb_i8_i4(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (rb_tree_node_i8_i4), POINTER :: x

    IF (.NOT.ASSOCIATED(tree%nil,x)) THEN
       CALL postorder_tree_walk(tree, x%lchild)
       CALL postorder_tree_walk(tree, x%rchild)
       WRITE(*,*) x%key
    END IF

  END SUBROUTINE postorder_tree_walk_rb_i8_i4




  RECURSIVE SUBROUTINE postorder_tree_walk_rb_i8(tree, x)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    TYPE (rb_tree_node_i8), POINTER :: x

    IF (.NOT.ASSOCIATED(tree%nil,x)) THEN
       CALL postorder_tree_walk(tree, x%lchild)
       CALL postorder_tree_walk(tree, x%rchild)
       WRITE(*,*) x%key
    END IF

  END SUBROUTINE postorder_tree_walk_rb_i8




  SUBROUTINE print_list(list)

    IMPLICIT NONE
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: list_node

    list_node => list%head
    DO WHILE (.NOT.ASSOCIATED(list_node,list%nil))
       WRITE(*,*) "key: ", list_node%key, " str: ", TRIM(ADJUSTL(list_node%str))
       list_node => successor(list, list_node)
    END DO

  END SUBROUTINE print_list




  SUBROUTINE print_tree(tree)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: list_node

    CALL init(list)
    CALL buildlist(tree, list, tree%root, 1_iprec4)
    CALL print_list(list)
    !do while (list%head

  END SUBROUTINE print_tree




  RECURSIVE SUBROUTINE buildlist(tree, list, x, key)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    TYPE (rb_tree_node_i8_i4), POINTER :: x
    INTEGER(iprec4), INTENT(in) :: key
    CHARACTER(len=1024) :: str
    LOGICAL :: error

    ! Same idea as in, for example, preordertreewalk, i.e., all levels are
    ! scanned from left to right. Add string to an element of a linked
    ! circular list where the key is containing the level information of
    ! the tree. 

    IF (.NOT.ASSOCIATED(x,tree%nil)) THEN
       ! Write correct string according to the number of children
       CALL toString(x%key, str, error)
       IF (x%red) THEN
          str = TRIM(str) // "R"
       ELSE
          str = TRIM(str) // "B"
       END IF
       IF (.NOT.ASSOCIATED(x%lchild,tree%nil) .AND. &
            .NOT.ASSOCIATED(x%rchild, tree%nil)) THEN
          str =  " _" // TRIM(ADJUSTL(str)) // "_ "
       ELSE IF (.NOT.ASSOCIATED(x%lchild,tree%nil) .AND. &
            ASSOCIATED(x%rchild, tree%nil)) THEN
          str = " _" // TRIM(ADJUSTL(str)) // " "
       ELSE IF (ASSOCIATED(x%lchild,tree%nil) .AND. &
            .NOT.ASSOCIATED(x%rchild, tree%nil)) THEN
          str = " " // TRIM(ADJUSTL(str)) // "_ "
       ELSE IF (ASSOCIATED(x%lchild,tree%nil) .AND. &
            ASSOCIATED(x%rchild, tree%nil)) THEN
          str = " " // TRIM(ADJUSTL(str)) // " "
       END IF

       CALL insert_list_node(list, key, str)

       ! When stepping down to the childrens level the level is
       ! increased by one:
       CALL buildlist(tree, list, x%lchild, key+1_iprec4)
       CALL buildlist(tree, list, x%rchild, key+1_iprec4)
    END IF

  END SUBROUTINE buildlist





  FUNCTION search_dbl_lnkd_list_i4(list, key) RESULT(node)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i4), POINTER :: list
    INTEGER(iprec4), INTENT(in) :: key
    TYPE (dbl_lnkd_list_node_i4), POINTER :: node

    node => list%nil%next
    DO WHILE (node%key /= key .AND. &
         .NOT.ASSOCIATED(node,list%nil))
       node => node%next       
    END DO

  END FUNCTION search_dbl_lnkd_list_i4




  FUNCTION search_dbl_lnkd_list_i8(list, key) RESULT(node)

    IMPLICIT NONE
    TYPE (dbl_lnkd_list_i8), POINTER :: list
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (dbl_lnkd_list_node_i8), POINTER :: node

    node => list%nil%next
    DO WHILE (node%key /= key .AND. &
         .NOT.ASSOCIATED(node,list%nil))
       node => node%next       
    END DO

  END FUNCTION search_dbl_lnkd_list_i8




  FUNCTION search_lnkd_list_ch1024_i4(list, key) RESULT(node)

    IMPLICIT NONE
    TYPE (lnkd_list_ch1024_i4), POINTER :: list
    INTEGER, INTENT(in) :: key
    TYPE (lnkd_list_node_ch1024_i4), POINTER :: node

    node => list%head
    DO WHILE (node%key /= key .AND. &
         .NOT.ASSOCIATED(node,list%nil))
       node => node%next       
    END DO

  END FUNCTION search_lnkd_list_ch1024_i4





  FUNCTION search_rb_tree_ch32(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32), POINTER :: tree
    CHARACTER(len=*), INTENT(in) :: key
    TYPE (rb_tree_node_ch32), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key == y%key) THEN
          EXIT
       ELSE IF (key < y%key) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_ch32





  FUNCTION search_rb_tree_ch32_8r8(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_ch32_8r8), POINTER :: tree
    CHARACTER(len=*), INTENT(in) :: key
    TYPE (rb_tree_node_ch32_8r8), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key == y%key) THEN
          EXIT
       ELSE IF (key < y%key) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_ch32_8r8





  FUNCTION search_rb_tree_r16_i4(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4), POINTER :: tree
    REAL(rprec16), INTENT(in) :: key
    TYPE (rb_tree_node_r16_i4), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (ABS(key - y%key) < 100.0_rprec16*EPSILON(key)) THEN
          EXIT
       ELSE IF (key < y%key) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_r16_i4




  FUNCTION search_rb_tree_i8_i4(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (rb_tree_node_i8_i4), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key == y%key) THEN
          EXIT
       ELSE IF (key < y%key) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_i8_i4




  FUNCTION search_rb_tree_i8_ch16arr(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_ch16arr), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (rb_tree_node_i8_ch16arr), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key == y%key) THEN
          EXIT
       ELSE IF (key < y%key) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_i8_ch16arr




  FUNCTION search_rb_tree_r16_i4arr(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_r16_i4arr), POINTER :: tree
    REAL(rprec16), INTENT(in) :: key
    TYPE (rb_tree_node_r16_i4arr), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key >= y%lkey .AND. key <= y%ukey) THEN
          EXIT
       ELSE IF (key < y%lkey) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_r16_i4arr




  FUNCTION search_rb_tree_i8_i4arr(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8_i4arr), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (rb_tree_node_i8_i4arr), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key >= y%lkey .AND. key <= y%ukey) THEN
          EXIT
       ELSE IF (key < y%lkey) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_i8_i4arr




  FUNCTION search_rb_tree_i8(tree, key) RESULT(y)

    IMPLICIT NONE
    TYPE (rb_tree_i8), POINTER :: tree
    INTEGER(iprec8), INTENT(in) :: key
    TYPE (rb_tree_node_i8), POINTER :: y

    y => tree%root
    DO WHILE (.NOT.ASSOCIATED(tree%nil,y))
       IF (key == y%key) THEN
          EXIT
       ELSE IF (key < y%key) THEN
          y => y%lchild
       ELSE ! key > y%key
          y => y%rchild
       END IF
    END DO

  END FUNCTION search_rb_tree_i8





  FUNCTION reallocate_r16_i4arr_node_1(array, n)

    IMPLICIT NONE
    TYPE (r16_i4arr_node), DIMENSION(:), POINTER :: array
    TYPE (r16_i4arr_node), DIMENSION(:), POINTER :: reallocate_r16_i4arr_node_1
    INTEGER, INTENT(in) :: n
    INTEGER :: nold, i

    ALLOCATE(reallocate_r16_i4arr_node_1(n))
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array, dim=1)
       DO i=1,MIN(n,nold)
          reallocate_r16_i4arr_node_1(i)%key = array(i)%key          
          reallocate_r16_i4arr_node_1(i)%data_array => array(i)%data_array
       END DO
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_r16_i4arr_node_1



  FUNCTION reallocate_i8_i4arr_node_1(array, n)

    IMPLICIT NONE
    TYPE (i8_i4arr_node), DIMENSION(:), POINTER :: array
    TYPE (i8_i4arr_node), DIMENSION(:), POINTER :: reallocate_i8_i4arr_node_1
    INTEGER, INTENT(in) :: n
    INTEGER :: nold, i

    ALLOCATE(reallocate_i8_i4arr_node_1(n))
    IF (ASSOCIATED(array)) THEN
       nold = SIZE(array, dim=1)
       DO i=1,MIN(n,nold)
          reallocate_i8_i4arr_node_1(i)%key = array(i)%key          
          reallocate_i8_i4arr_node_1(i)%data_array => array(i)%data_array
       END DO
       DEALLOCATE(array)
    END IF

  END FUNCTION reallocate_i8_i4arr_node_1



END MODULE data_structures
