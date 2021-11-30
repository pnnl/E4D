module forward_seq

  !_________________________________________________________________
  !! Author: Tim Johnson
  !!
  !! email: tj@pnnl.gov
  !!
  !! Purpose: This module computes forward ert simulation runs given
  !! a list of conductivity distributions (sigma files)
  !!
  !! Revision History: created January 26, 2017
  !!
  !__________________________________________________________________

  use vars
  use input
  !use output
  use master
  implicit none

contains
  !__________________________________________________________________
  subroutine exec_mforward
    implicit none
    integer :: i, nt
    character*80 :: str

    !setup a mpi communication group for the slaves only (see module: master)
    call build_sgroupm

    !loop over the conductivity files
    nt = size(tl_cfils(:,1))
    do i=1,nt
       sigfile = tl_cfils(i,1)
       call read_conductivity
       call send_sigma
       if(i_flag) call send_sigmai

       
       if(i.eq.1) then
          call setup_forward
          call nreport(65)
          call get_time
          call treport(1)
       
          if(i_flag) call send_command(102)
          call nreport(66)
          call get_time
          call treport(1)
       end if

       !build the coupling matrix (command 3 for slaves)
       !if this is a complex conductivity forward run then 
       !setup the complex conductivity coupling matrix also
       call send_command(3)
       if(i_flag) call send_command(103)

        !execute a forward run
       !call nreport(67)
       !call nreport(68)
       if(i_flag) then
          !call send_command(5)
          call run_forward
          call nreport(71)
          !call send_command(5)
          call run_forwardi
       else
          !call send_command(5)
          call run_forward
       end if

       call get_dpred 
       write(str,*) tl_cfils(i,2)
  
       if (ms_flag==1) then
		  call print_ms_dpd_mf(str)
		  call mwrite_ms_pots(i)
		  call mbuild_ms_srv(i)
	   else
		  call print_dpd_mf(str)
		  call mwrite_pots(i)
		  call mbuild_srv(i)
	   end if

       

    end do
  end subroutine exec_mforward
  !__________________________________________________________________



   !_______________________________________________________________________________________
  subroutine mbuild_srv(indx)
    implicit none
    character*80 :: fname=""
    integer :: i,indx

    write(fname,"(A,A)") trim(tl_cfils(indx,2)),'.srv'
    call nreport(69)
    open(12,file=fname,status='replace',action='write')
    write(12,*) ne
    do i=1,ne
       write(12,"(I10,3F15.5,I10)") i,e_pos(i,1)+xorig,e_pos(i,2)+yorig,e_pos(i,3)+zorig,int(e_pos(i,4))
    end do
    write(12,*)
    write(12,*) nm
    if(i_flag) then
       do i=1,nm
          write(12,"(I8,4I10,4G15.5)") i,s_conf(i,1:4),dpred(i),0.05*abs(dpred(i)),dpredi(i),0.05*abs(dpredi(i))
       end do
    else
       do i=1,nm
          write(12,"(I8,4I10,2G15.5)") i,s_conf(i,1:4),dpred(i),0.05*abs(dpred(i))+0.01
       end do
    end if
    close(12)
    
  end subroutine mbuild_srv
  !_______________________________________________________________________________________

   !_______________________________________________________________________________________
  subroutine mbuild_ms_srv(indx)
    implicit none
    character*80 :: fname=""
    integer :: i,indx

    write(fname,"(A,A)") trim(tl_cfils(indx,2)),'.srv'
    call nreport(69)
    open(12,file=fname,status='replace',action='write')
    write(12,*) ne
    do i=1,ne
       write(12,"(I10,3F15.5,I10)") i,e_pos(i,1)+xorig,e_pos(i,2)+yorig,e_pos(i,3)+zorig,int(e_pos(i,4))
    end do
    write(12,*)
    write(12,*) nm
    if(i_flag) then
       do i=1,nm
          write(12,"(I8,6I10,6G15.5)") i,ms_conf(i,1:6),ms_currents(i,1:2),dpred(i),0.05*abs(dpred(i)),dpredi(i),0.05*abs(dpredi(i))
       end do
    else
       do i=1,nm
          write(12,"(I8,6I10,4G15.5)") i,ms_conf(i,1:6),ms_currents(i,1:2),dpred(i),0.05*abs(dpred(i))+0.01
       end do
    end if
    close(12)
    
  end subroutine mbuild_ms_srv
  !_______________________________________________________________________________________




  !_________________________________________________________________
  subroutine mwrite_pots(indx)
    implicit none
    integer :: dp_flag,pot_flag,npot,o_opt,ist,jflag,indx
    logical :: fcheck
    character*80 :: dp_file
    character*20 :: fname,jformat
    integer :: i,a,b,j,emin,emax,ra,rb
    integer, dimension(2) :: spack
    real, dimension(nnodes) :: pa,pb,rp,cp
    integer ::  status(MPI_STATUS_SIZE)

   
    inquire(file=trim(outfile),exist=fcheck); if(.not.fcheck) goto 10
  
    open(15,file=outfile,status='old',action='read')
    read(15,*,IOSTAT=ist) dp_flag; if(ist.ne.0) goto 11
    read(15,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
    read(15,*,IOSTAT=ist) npot   ; if(ist.ne.0) goto 13

    if(npot>0) then
       allocate(ipot(npot,2))
       do i=1,npot
          read(15,*,IOSTAT=ist) ipot(i,1); if(ist.ne.0) goto 14
       end do
    else
       close(15)
       return
    end if
    close(15)

    call nreport(21)
    do i=1,npot
       pa=0
       pb=0
       if(ipot(i,1)>nm) goto 100
       a=s_conf(ipot(i,1),1)
       b=s_conf(ipot(i,1),2)
       do j=1,n_rank-1
          emin=eind(j,1); emax=eind(j,2)
          if((emin .le. a) .and. (emax .ge. a)) ra = j	
          if((emin .le. b) .and. (emax .ge. b)) rb = j 
       end do
    
       if(a .ne. 0) then
          spack(1) = ra
          spack(2) = a
          call send_commando(23)
          call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
          call MPI_RECV(pa,nnodes,MPI_REAL,ra,0,E4D_COMM,status,ierr)
       end if
       
       if(b .ne. 0) then
          spack(1) = rb
          spack(2) = b
          call send_commando(23)
          call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
          call MPI_RECV(pb,nnodes,MPI_REAL,rb,0,E4D_COMM,status,ierr)
       end if
       
       do j=1,nnodes
          rp(j)=pa(j)-pb(j)
       end do


       if(i_flag) then
          pa=0
          pb=0
          if(a .ne. 0) then
             spack(1) = ra
             spack(2) = a
             call send_commando(123)
             call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
             call MPI_RECV(pa,nnodes,MPI_REAL,ra,0,E4D_COMM,status,ierr)
          end if
         
          if(b .ne. 0) then
             spack(1) = rb
             spack(2) = b
             call send_commando(123)
             call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
             call MPI_RECV(pb,nnodes,MPI_REAL,rb,0,E4D_COMM,status,ierr)
          end if
          
          do j=1,nnodes
             cp(j) = pa(j)-pb(j)
          end do
    
       end if

       write(fname,"(A,I0)") "pot_"//trim(tl_cfils(indx,2))//".",ipot(i,1)
       open(27,file=fname,status='replace',action='write')
       if(i_flag) then
          write(27,*) nnodes, 2, ipot(i,1)
          do j=1,nnodes
             write(27,*) rp(j),cp(j)
          end do
       else
          write(27,*) nnodes, 1, ipot(i,1)
          do j=1,nnodes
             write(27,*) rp(j)
          end do
       end if
       close(27)

       
100    continue
    end do
    deallocate(ipot)
    return

   
10  continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' Cannot find the output options file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' Cannot find the output options file: ',trim(outfile)

      return
      
11    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' The was a problem reading the first line in the output file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' There was a problem reading the first line in the output file: ',trim(outfile)
      return

12    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) 'There was a problem reading the predicted data file name in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) 'The was a problem reading the predicted data file in: ',trim(outfile)
      return

13    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      return

14    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      close(51)
      write(*,*)
      write(*,*) ' There was a problem reading potential field index: ',i,' in: ',trim(outfile)
       deallocate(ipot)
      return



  end subroutine mwrite_pots
  !___________________________________________________________


  !_________________________________________________________________
  subroutine mwrite_ms_pots(indx)
    implicit none
    integer :: dp_flag,pot_flag,npot,o_opt,ist,jflag,indx
    logical :: fcheck
    character*80 :: dp_file
    character*20 :: fname,jformat
    integer :: i,a1,a2,b1,b2,j,emin,emax,ra1,ra2,rb1,rb2
    integer, dimension(2) :: spack
    real, dimension(nnodes) :: pa1,pa2,pb1,pb2,rp,cp
    integer ::  status(MPI_STATUS_SIZE)

   
    inquire(file=trim(outfile),exist=fcheck); if(.not.fcheck) goto 10
  
    open(15,file=outfile,status='old',action='read')
    read(15,*,IOSTAT=ist) dp_flag; if(ist.ne.0) goto 11
    read(15,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
    read(15,*,IOSTAT=ist) npot   ; if(ist.ne.0) goto 13

    if(npot>0) then
       allocate(ipot(npot,2))
       do i=1,npot
          read(15,*,IOSTAT=ist) ipot(i,1); if(ist.ne.0) goto 14
       end do
    else
       close(15)
       return
    end if
    close(15)

    call nreport(21)
    do i=1,npot
       pa1=0
       pb1=0
       
       pa2=0
       pb2=0
       
       if(ipot(i,1)>nm) goto 100
       a1=ms_conf(ipot(i,1),1)
       b1=ms_conf(ipot(i,1),2)

       a2=ms_conf(ipot(i,1),3)
       b2=ms_conf(ipot(i,1),4)
       
       do j=1,n_rank-1
          emin=eind(j,1); emax=eind(j,2)
          if((emin .le. a1) .and. (emax .ge. a1)) ra1 = j	
          if((emin .le. b1) .and. (emax .ge. b1)) rb1 = j 
          
          if((emin .le. a2) .and. (emax .ge. a2)) ra2 = j	
          if((emin .le. b2) .and. (emax .ge. b2)) rb2 = j 

       end do
    
       if(a1 .ne. 0) then
          spack(1) = ra1
          spack(2) = a1
          call send_commando(23)
          call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
          call MPI_RECV(pa1,nnodes,MPI_REAL,ra1,0,E4D_COMM,status,ierr)
       end if
       
       if(b1 .ne. 0) then
          spack(1) = rb1
          spack(2) = b1
          call send_commando(23)
          call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
          call MPI_RECV(pb1,nnodes,MPI_REAL,rb1,0,E4D_COMM,status,ierr)
       end if

       if(a2 .ne. 0) then
          spack(1) = ra2
          spack(2) = a2
          call send_commando(23)
          call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
          call MPI_RECV(pa2,nnodes,MPI_REAL,ra2,0,E4D_COMM,status,ierr)
       end if
       
       if(b2 .ne. 0) then
          spack(1) = rb2
          spack(2) = b2
          call send_commando(23)
          call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
          call MPI_RECV(pb2,nnodes,MPI_REAL,rb2,0,E4D_COMM,status,ierr)
       end if

       
       do j=1,nnodes
          rp(j)=pa1(j)-pb1(j) + pa2(j)-pb2(j)
       end do


       if(i_flag) then
          pa1=0
          pb1=0
          
          pa2=0
          pb2=0
          
          if(a1 .ne. 0) then
             spack(1) = ra1
             spack(2) = a1
             call send_commando(123)
             call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
             call MPI_RECV(pa1,nnodes,MPI_REAL,ra1,0,E4D_COMM,status,ierr)
          end if
         
          if(b1 .ne. 0) then
             spack(1) = rb1
             spack(2) = b1
             call send_commando(123)
             call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
             call MPI_RECV(pb1,nnodes,MPI_REAL,rb1,0,E4D_COMM,status,ierr)
          end if

          if(a2 .ne. 0) then
             spack(1) = ra2
             spack(2) = a2
             call send_commando(123)
             call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
             call MPI_RECV(pa2,nnodes,MPI_REAL,ra2,0,E4D_COMM,status,ierr)
          end if
         
          if(b2 .ne. 0) then
             spack(1) = rb2
             spack(2) = b2
             call send_commando(123)
             call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
             call MPI_RECV(pb2,nnodes,MPI_REAL,rb2,0,E4D_COMM,status,ierr)
          end if

          
          do j=1,nnodes
             cp(j) = pa1(j)-pb1(j) + pa2(j)-pb2(j)
          end do
    
       end if

       write(fname,"(A,I0)") "pot_"//trim(tl_cfils(indx,2))//".",ipot(i,1)
       open(27,file=fname,status='replace',action='write')
       if(i_flag) then
          write(27,*) nnodes, 2, ipot(i,1)
          do j=1,nnodes
             write(27,*) rp(j),cp(j)
          end do
       else
          write(27,*) nnodes, 1, ipot(i,1)
          do j=1,nnodes
             write(27,*) rp(j)
          end do
       end if
       close(27)

       
100    continue
    end do
    deallocate(ipot)
    return

   
10  continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' Cannot find the output options file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' Cannot find the output options file: ',trim(outfile)

      return
      
11    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' The was a problem reading the first line in the output file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' There was a problem reading the first line in the output file: ',trim(outfile)
      return

12    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) 'There was a problem reading the predicted data file name in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) 'The was a problem reading the predicted data file in: ',trim(outfile)
      return

13    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      return

14    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      close(51)
      write(*,*)
      write(*,*) ' There was a problem reading potential field index: ',i,' in: ',trim(outfile)
       deallocate(ipot)
      return



  end subroutine mwrite_ms_pots
  !___________________________________________________________


  !__________________________________________________________________________
  subroutine print_dpd_mf(dp_fout_name)
    implicit none
    character*80, intent(in) :: dp_fout_name
    integer :: i
    write(*,*) " WRITING e4d_",trim(adjustl(dp_fout_name))//".dpd"
    open(15,file="e4d_"//trim(adjustl(dp_fout_name))//".dpd",status='replace',action='write')
    write(15,*) nm
    if(i_flag) then
       do i=1,nm
          write(15,1002) i,s_conf(i,1:4),dobs(i),dpred(i),dobsi(i),dpredi(i)
       end do
       
    else
       do i=1,nm
          write(15,1001) i,s_conf(i,1:4),dobs(i),dpred(i)
       end do
    end if
    close(15)
    return
1001 format(1I8,4I8,2g15.6)
1002 format(1I8,4I8,4g15.6)
  end subroutine print_dpd_mf
  !__________________________________________________________________________


  !__________________________________________________________________________
  subroutine print_ms_dpd_mf(dp_fout_name)
    implicit none
    character*80, intent(in) :: dp_fout_name
    integer :: i
    write(*,*) " WRITING e4d_",trim(adjustl(dp_fout_name))//".dpd"
    open(15,file="e4d_"//trim(adjustl(dp_fout_name))//".dpd",status='replace',action='write')
    write(15,*) nm
    if(i_flag) then
       do i=1,nm
          write(15,1002) i,ms_conf(i,1:6),dobs(i),dpred(i),dobsi(i),dpredi(i)
       end do
       
    else
       do i=1,nm
          write(15,1001) i,ms_conf(i,1:6),dobs(i),dpred(i)
       end do
    end if
    close(15)
    return
1001 format(1I8,6I8,2g15.6)
1002 format(1I8,6I8,4g15.6)
  end subroutine print_ms_dpd_mf
  !__________________________________________________________________________

end module forward_seq
