!     # Create and save a frame for animation
      write(timechar,'(i5)') itime

      Open(Unit=87,File=Trim(Out_Folder_Path)//'/tecplot-vel-'//timechar//'.txt') 
      Rewind 87

      do i=1,numCells
      write(87,'(a,es11.4,1x,es11.4,1x,es11.4,a)') '(',u(i),v(i),w(i),')'
      enddo

      Close(87)


      ! call execute_command_line("tec360 -b -p makro-film.mcr")
      ! call execute_command_line("mv untitled100.png "//timechar//".png")
