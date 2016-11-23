!***********************************************************************
!
 subroutine writehistory
!
!***********************************************************************
!  Writes values at each monitoring point, at each time-step for 
!  transient simulation, to a specific monitor file.
!
!***********************************************************************
!
  use types
  use parameters
  use variables
  use title_mod
  use k_epsilon_std

  implicit none
!
!***********************************************************************
!
   integer :: inp,imon

  if(ltransient) then
  do imon=1,mpoints
    read(89,*) inp
    write(91+imon,'(2x,1p7e14.5,2x)') time,u(inp),v(inp),w(inp),te(inp),ed(inp)
  end do
  rewind 89
  end if

end subroutine
