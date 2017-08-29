! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
!
! Made F conformant by Walt Brainerd
! Modified for use in Cappuccino code by Nikola Mirkov (nmirkov@vinca.rs)
!

module qsort_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A, A1,A2,A3,A4, B1,B2,B3,B4)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: A1,A2,A3,A4
  integer, intent(in out), dimension(:) :: B1,B2,B3,B4
  integer :: iq
  if(size(A) > 1) then
     call Partition(A, iq, A1,A2,A3,A4, B1,B2,B3,B4)
     call QsortC( A(:iq-1), A1(:iq-1),A2(:iq-1),A3(:iq-1),A4(:iq-1), B1(:iq-1),B2(:iq-1),B3(:iq-1),B4(:iq-1) )
     call QsortC( A(iq:), A1(iq:),A2(iq:),A3(iq:),A4(iq:), B1(iq:),B2(iq:),B3(iq:),B4(iq:) )
  endif
end subroutine QsortC

subroutine Partition(A, marker, A1,A2,A3,A4, B1,B2,B3,B4)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: A1,A2,A3,A4
  integer, intent(in out), dimension(:) :: B1,B2,B3,B4
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
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
        temp = A(i)
        A(i) = A(j)
        A(j) = temp

        ! do the companion matrix:

        temp = A1(i)
        A1(i) = A1(j)
        A1(j) = temp

        temp = A2(i)
        A2(i) = A2(j)
        A2(j) = temp

        temp = A3(i)
        A3(i) = A3(j)
        A3(j) = temp

        temp = A4(i)
        A4(i) = A4(j)
        A4(j) = temp

        ! do the other companion matrix:

        temp = B1(i)
        B1(i) = B1(j)
        B1(j) = temp

        temp = B2(i)
        B2(i) = B2(j)
        B2(j) = temp

        temp = B3(i)
        B3(i) = B3(j)
        B3(j) = temp

        temp = B4(i)
        B4(i) = B4(j)
        B4(j) = temp

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition


recursive subroutine QsortC2(A, A1, B1,B2,B3,B4)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: A1
  integer, intent(in out), dimension(:) :: B1,B2,B3,B4
  integer :: iq
  if(size(A) > 1) then
     call Partition2(A, iq, A1, B1,B2,B3,B4)
     call QsortC2( A(:iq-1), A1(:iq-1), B1(:iq-1),B2(:iq-1),B3(:iq-1),B4(:iq-1) )
     call QsortC2( A(iq:), A1(iq:), B1(iq:),B2(iq:),B3(iq:),B4(iq:) )
  endif
end subroutine QsortC2

subroutine Partition2(A, marker, A1, B1,B2,B3,B4)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: A1
  integer, intent(in out), dimension(:) :: B1,B2,B3,B4
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
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
        temp = A(i)
        A(i) = A(j)
        A(j) = temp

        ! do the companion matrix:

        temp = A1(i)
        A1(i) = A1(j)
        A1(j) = temp

        ! do the other companion matrix:

        temp = B1(i)
        B1(i) = B1(j)
        B1(j) = temp

        temp = B2(i)
        B2(i) = B2(j)
        B2(j) = temp

        temp = B3(i)
        B3(i) = B3(j)
        B3(j) = temp

        temp = B4(i)
        B4(i) = B4(j)
        B4(j) = temp

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition2


recursive subroutine QsortC3(A, B1,B2,B3,B4)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B1,B2,B3,B4
  integer :: iq
  if(size(A) > 1) then
     call Partition3(A, iq, B1,B2,B3,B4)
     call QsortC3( A(:iq-1), B1(:iq-1),B2(:iq-1),B3(:iq-1),B4(:iq-1) )
     call QsortC3( A(iq:), B1(iq:),B2(iq:),B3(iq:),B4(iq:) )
  endif
end subroutine QsortC3

subroutine Partition3(A, marker, B1,B2,B3,B4)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B1,B2,B3,B4
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
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
        temp = A(i)
        A(i) = A(j)
        A(j) = temp

        ! do the companion matrix:

        temp = B1(i)
        B1(i) = B1(j)
        B1(j) = temp

        temp = B2(i)
        B2(i) = B2(j)
        B2(j) = temp

        temp = B3(i)
        B3(i) = B3(j)
        B3(j) = temp

        temp = B4(i)
        B4(i) = B4(j)
        B4(j) = temp

     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition3

end module qsort_c_module
