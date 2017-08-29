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
  use parameters, only: ltransient,time,mpoints
  use variables, only: u,v,w,te,ed

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
