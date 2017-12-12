! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd
! Smaller adjustments for personal usage by Nikola Mirkov

module qsort_c_module

implicit none

  interface Qsort
    module procedure QsortC
    module procedure QsortCInt
  end interface

  interface QsortIndx
    module procedure QsortCIndx
    module procedure QsortCIntIndx
  end interface  

private

public :: Qsort, QsortIndx

contains

recursive subroutine QsortC(A)
  real, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition


recursive subroutine QsortCInt(A)
  integer, intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call PartitionInt(A, iq)
     call QsortCInt(A(:iq-1))
     call QsortCInt(A(iq:))
  endif
end subroutine QsortCInt

subroutine PartitionInt(A, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  ! integer :: temp
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        ! temp = A(i)
        ! A(i) = A(j)
        ! A(j) = temp
        call swapInt( A(i), A(j) )
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine PartitionInt


recursive subroutine QsortCIndx(A,ind)
  real, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: ind
  integer :: iq

  if(size(A) > 1) then
     call PartitionIndx(A, ind, iq)
     call QsortCIndx(A(:iq-1), ind(:iq-1))
     call QsortCIndx(A(iq:), ind(iq:))
  endif
end subroutine QsortCIndx

subroutine PartitionIndx(A, ind, marker)
  real, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: ind
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
         
        call swapInt( ind(i), ind(j) )

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine PartitionIndx


recursive subroutine QsortCIntIndx(A,ind)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: ind
  integer :: iq

  if(size(A) > 1) then
     call PartitionIntIndx(A, ind, iq)
     call QsortCIntIndx(A(:iq-1), ind(:iq-1))
     call QsortCIntIndx(A(iq:), ind(iq:))
  endif
end subroutine QsortCIntIndx

subroutine PartitionIntIndx(A, ind, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: ind
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
         
        call swapInt( ind(i), ind(j) )

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine PartitionIntIndx


subroutine swap(a,b)
real, intent(inout) :: a,b
real :: temp
  temp = a
  a = b
  b = temp
end subroutine swap

subroutine swapInt(a,b)
integer, intent(inout) :: a,b
integer :: temp
  temp = a
  a = b
  b = temp
end subroutine swapInt

end module qsort_c_module
