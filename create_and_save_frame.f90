!     # Create and save a frame for animation
      Open(Unit=87,File=Trim(Out_Folder_Path)//'/tecplot-vel.plt') 
      Rewind 87
      Write(87,*) 'Title     = " "'
      Write(87,*) 'Variables = "X"'
      Write(87,*) '"Y"'
      Write(87,*) '"Z"'
      Write(87,*) '"U"'
      Write(87,*) '"V"'
      Write(87,*) '"W"'
      Write(87,*) 'Zone T=" "'
      ! Write(87,*) 'I=',Ni, ' ,J=',Nj, ' ,K=',Nk,', F=Point'
      ! Do k=1,nk; do i=1,ni; do j=1,nj;
      ! Inp=Lk(K)+Li(I)+J
      ! Write(87,*) Xc(Inp),Yc(Inp),Zc(Inp),U(Inp),V(Inp),W(Inp)
      ! Enddo; Enddo; Enddo 
      Close(87)

      write(timechar,'(i5)') itime
      call execute_command_line("tec360 -b -p makro-film.mcr")
      call execute_command_line("mv untitled100.png "//timechar//".png")
