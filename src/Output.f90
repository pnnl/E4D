module output
  
  use vars
  use report
  use reorder_mesh

  integer, dimension(:,:), allocatable :: ipot
  real, dimension(:,:), allocatable :: pot
  integer :: njrows_out
  integer, dimension(:), allocatable :: jrows_out
  logical :: print_jaco_rows = .false.
  
contains
 
  !_______________________________________________________________________________________
  subroutine output_dpred
    implicit none
    
    integer :: conflag,i,j,dp_flag,o_opt,opt,ist
    logical :: gs_flag,exst
    character*80 :: fname,dp_file
    real, dimension(:,:), allocatable :: pot
    real, dimension(:), allocatable :: fpot
 
    inquire(file=trim(outfile),exist=exst); if(.not.exst) goto 10
    open(15,file=trim(outfile),status='old',action='read')
    read(15,*,IOSTAT=ist) dp_flag;  if(ist.ne.0) goto 11
    read(15,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
    close(15)
    
    if(dp_flag==1) then

       !if(o_opt==1) then
       open(15,file=dp_file,status='replace',action='write')
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
    end if
    return
!!$       elseif(o_opt==2) then 
!!$          write(dp_file,"(A,A,I0)") trim(dp_file),'_',iter
!!$          open(15,file=dp_file,status='replace',action='write')
!!$          write(15,*) nm
!!$          if(opt==2) then
!!$             do i=1,nm
!!$                write(15,1002) i,s_conf(i,1:4),dobs(i),dpred(i),dobsi(i),dpredi(i)
!!$             end do
!!$          else
!!$             do i=1,nm
!!$                write(15,1001) i,s_conf(i,1:4),dobs(i),dpred(i)
!!$             end do
!!$          end if
!!$          close(15)
!!$       end if
!!$
!!$    end if

10  continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: Cannot find the output options file: ',trim(outfile)
      close(51)
      write(*, *) ' E4D: Cannot find the output options file: ',trim(outfile)
      return
      
11    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: The was a problem reading the first line in the output file: ',trim(outfile)
      write(51,*) ' E4D: aborting'
      close(51)
      write(*, *) ' E4D: The was a problem reading the first line in the output file: ',trim(outfile)
      write(*, *) ' E4D: aborting'
      return

12    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: The was a problem reading the predicted data file in: ',trim(outfile)
      write(51,*) ' E4D: aborting'
      close(51)
      write(*, *) ' E4D: The was a problem reading the predicted data file in: ',trim(outfile)
      write(*, *) ' E4D: aborting'
      return

13    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: The was a problem reading the number of potential fields to write in: ',trim(outfile)
      write(51,*) ' E4D: aborting'
      close(51)
      write(*, *) ' E4D: The was a problem reading the number of potential fields to write in: ',trim(outfile)
      write(*, *) ' E4D: aborting'
      return

14    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: The was a problem reading potential index: ',i,' in: ',trim(outfile)
      write(51,*) ' E4D: aborting'
      close(51)
      write(*, *) ' E4D: The was a problem reading potential index: ',i,' in: ',trim(outfile)
      write(*, *) ' E4D: aborting'

      return

   
1001 format(1I8,4I8,2g15.6)
1002 format(1I8,4I8,4g15.6)

  end subroutine output_dpred
  !_______________________________________________________________________________________

  !_______________________________________________________________________________________
  subroutine build_srv
    implicit none
    character*80 :: fname=""
    integer :: i

    write(fname,"(A,A)") trim(sigfile),'.srv'
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
    
  end subroutine build_srv
  !_______________________________________________________________________________________
  !_______________________________________________________________________________________
  subroutine write_sigma
    implicit none
    integer :: i
    character(20) :: fname=""
    
    write(fname,"(A,I0)") "sigma.",iter
    open(12,file=fname,status='replace',action='write')
    
    if(allocated(element_map)) then
       if(i_flag) then
          write(12,*) size(element_map), 2, chi2
          do i=1,size(element_map)
             if(element_map(i) .ne. 0) then
                write(12,*) sigma(element_map(i)), real(sigmai(element_map(i)))
             else
                write(12,*) -999,-999
             end if
          end do
       else
          write(12,*) size(element_map), 1, chi2
           do i=1,size(element_map)
             if(element_map(i) .ne. 0) then
                write(12,*) sigma(element_map(i))
             else
                write(12,*) -999
             end if
          end do
       end if
       
    else

       if(i_flag) then
          write(12,*) nelem, 2, chi2
       else
          write(12,*) nelem, 1, chi2
       end if
       
       if(i_flag) then
          do i=1,nelem
             write(12,*) sigma(i),real(sigmai(i))
          end do
       else
          do i=1,nelem
             write(12,*) sigma(i)
          end do
       end if

    end if
    close(12)
  end subroutine write_sigma
  !_______________________________________________________________________________________
  !_______________________________________________________________________________________
  subroutine write_sigmai
    implicit none
    integer :: i
    character(20) :: fname=""
    
    write(fname,"(A,I0)") "sigmai.",iter
    open(12,file=fname,status='replace',action='write')

    if(allocated(element_map)) then
       write(12,*) size(element_map), 2
       do i=1,size(element_map)
          if(element_map(i) .ne. 0) then
             write(12,*) sigma(element_map(i)), real(sigmai(element_map(i)))
          else
             write(12,*) -999,-999
          end if
       end do
    else
       write(12,*) nelem, 2
       do i=1,nelem
          write(12,*) sigma(i), sigmai(i)
       end do
    end if
       close(12)
  end subroutine write_sigmai
  !_______________________________________________________________________________________

  !_______________________________________________________________________________________
  subroutine write_sigma_rttl(fname)
    implicit none
    character*40 :: fname, oname
    integer :: i,np
    
   
    !!find the location of the . in the file name
    np=40 
    do i=1,np
       if(fname(i:i) == '.') then
          np=i
          exit
      end if
    end do

    oname=fname
    oname(np:np+3) = '.sig'
    open(12,file=trim(oname)//'.part',status='replace',action='write')
    if(i_flag) then
       write(12,*) nelem, 2, chi2
    else
       write(12,*) nelem, 1, chi2
    end if
    
    if(i_flag) then
       do i=1,nelem
          write(12,*) sigma(i),real(sigmai(i))
       end do
    else 
       do i=1,nelem
          write(12,*) sigma(i)
       end do
    end if
    close(12)
    call system('mv '//trim(oname)//'.part '//trim(oname))

  end subroutine write_sigma_rttl
  !_______________________________________________________________________________________
  subroutine write_sigma_tl(tm)
    implicit none
    integer :: i
    real :: tm
    character(8) ::ts
    character(20) :: fname=""
    write(ts,'(f8.3)') tm
    write(fname,"(A6,A8)") "tl_sigma_",adjustl(ts)
    open(12,file=fname,status='replace',action='write')
    
    if(allocated(element_map)) then
       if(i_flag) then
          write(12,*) size(element_map), 2, chi2
          do i=1,size(element_map)
             if(element_map(i) .ne. 0) then
                write(12,*) sigma(element_map(i)), real(sigmai(element_map(i)))
             else
                write(12,*) -999,-999
             end if
          end do
       else
          write(12,*) size(element_map), 1, chi2
           do i=1,size(element_map)
             if(element_map(i) .ne. 0) then
                write(12,*) sigma(element_map(i))
             else
                write(12,*) -999
             end if
          end do
       end if
    else
       if(i_flag) then
          write(12,*) nelem, 2, chi2
       else
          write(12,*) nelem, 1, chi2
       end if
       
       if(i_flag) then
          do i=1,nelem
             write(12,*) sigma(i),real(sigmai(i))
          end do
       else
          do i=1,nelem
             write(12,*) sigma(i)
          end do
       end if
    end if
    close(12)
  end subroutine write_sigma_tl
  !_______________________________________________________________________________________



  !_______________________________________________________________________________________
  subroutine write_pots
    implicit none
    integer :: dp_flag,pot_flag,npot,o_opt,ist,jflag
    logical :: fcheck
    character*80 :: dp_file
    character*20 :: fname,jformat
    integer :: i,a,b,j,emin,emax,ra,rb
    integer, dimension(2) :: spack
    real, dimension(nnodes) :: pa,pb,rp,cp
    integer ::  status(MPI_STATUS_SIZE)

   
    inquire(file=trim(outfile),exist=fcheck); if(.not.fcheck) goto 10
  
    call nreport(21)
    open(15,file=outfile,status='old',action='read')
    read(15,*,IOSTAT=ist) dp_flag; if(ist.ne.0) goto 11
    read(15,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
    read(15,*,IOSTAT=ist) npot   ; if(ist.ne.0) goto 13

    if(npot>0) then
       allocate(ipot(npot,2))
       do i=1,npot
          read(15,*,IOSTAT=ist) ipot(i,1); if(ist.ne.0) goto 14
       end do
    end if
  
    !!read the jacobian output flag
    read(15,*,IOSTAT=ist) jflag
    if(ist .ne.0) then
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) ' There was a problem Jacobian matrix output option: ',i,' in: ',trim(outfile)
       write(51,*) ' Not printing the Jacobian matrix.'
       close(51)
       write(*,*)
       write(*, *) ' There was a problem Jacobian matrix output option: ',i,' in: ',trim(outfile)
       write(*, *) ' Not printing the Jacobian matrix.'
       goto 9
    end if
    read(15,*,IOSTAT=ist) jformat
    if(ist .ne.0) then
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) ' There was a problem Jacobian output format option: ',i,' in: ',trim(outfile)
       write(51,*) ' Printing in binary format'
       close(51)
       write(*,*)
       write(*, *) ' There was a problem Jacobian matrix output option: ',i,' in: ',trim(outfile)
       write(*, *) ' Printing in binary format'
      
    end if
  
    if(jflag==1) then
       jaco_out_opt = .true.
       jaco_ascii_opt = .true.
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) ' Printing Jacobian matrix in '
       write(*, *) ' Printing Jacobian matrix in '

       if(trim(jformat)=='ASCII'.or.trim(jformat)=='ascii') then
          jaco_ascii_opt = .true.
          write(51,*) ' ascii format'
          write(*, *) ' ascii format'
       else
          jaco_ascii_opt = .false.
          write(51,*) ' binary format'
          write(*, *) ' binary format'
       end if
       close(51)
    else
       jaco_out_opt = .false.
        open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) ' Not printing Jacobian matrix'
       write(*, *) ' Not printing Jacobian matrix'
       close(51)
    end if


9   continue
    close(15)
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

       write(fname,"(A,I0)") "potential.",ipot(i,1)
       open(27,file=fname,status='replace',action='write')
       
       if(allocated(node_map)) then

          if(i_flag) then
             write(27,*) size(node_map), 2, ipot(i,1)
             do j=1,size(node_map)
                if(node_map(i) .ne. 0) then
                   write(27,*) rp(node_map(j)),cp(node_map(j))
                else
                   write(27,*) -999, -999
                end if
             end do
          else
             write(27,*) size(node_map), 1, ipot(i,1)
             do j=1,size(node_map)
                if(node_map(i) .ne. 0) then
                   write(27,*) rp(node_map(j))
                else
                   write(27,*) -999
                end if
             end do
          end if

       else
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

       end if
       close(27)

       
100    continue
    end do
    return

   
10  continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: Cannot find the output options file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*, *) ' E4D: Cannot find the output options file: ',trim(outfile)
      return
      
11    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: The was a problem reading the first line in the output file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*, *) ' E4D: There was a problem reading the first line in the output file: ',trim(outfile)
      return

12    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: There was a problem reading the predicted data file name in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*, *) ' E4D: The was a problem reading the predicted data file in: ',trim(outfile)
      return

13    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: There was a problem reading the number of potential fields to write in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*, *) ' E4D: There was a problem reading the number of potential fields to write in: ',trim(outfile)
      return

14    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      close(51)
      write(*,*)
      write(*, *) ' E4D: There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      return



      

  end subroutine write_pots
  !_______________________________________________________________________________________

  !_______________________________________________________________________________________
  subroutine write_sigiter
    implicit none
    character*40 :: fiter
    integer :: i
   
 
    write(fiter,"(A,I0)") "si.",iter
    open(23,file=trim(fiter),status='replace',action='write')
    if(i_flag) then
       write(23,*) nelem, 2, chi2
       do i=1,nelem
          write(23,*) sigma(i),sigmai(i)
       end do

    else
       write(23,*) nelem, 1, chi2
       do i=1,nelem
          write(23,*) sigma(i)
       end do

    end if
    close(23)

  end subroutine write_sigiter
  !_______________________________________________________________________________________

   !_______________________________________________________________________________________
  subroutine check_jrows_out
    implicit none
    integer :: dp_flag,npot,o_opt,ist,nfres,junk
    logical :: fcheck
    character*80 :: dp_file
    character*20 :: fname
    integer :: i,a,j,smin,smax,ra
    integer, dimension(2) :: spack
    real, dimension(nnodes) :: pa
    integer ::  status(MPI_STATUS_SIZE)

    
    inquire(file=trim(outfile),exist=fcheck); if(.not.fcheck) goto 10
    
    !call nreport_fmm(21)
    open(15,file=outfile,status='old',action='read')
    read(15,*,IOSTAT=ist) dp_flag; if(ist.ne.0) goto 11
    read(15,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
    read(15,*,IOSTAT=ist) npot   ; if(ist.ne.0) goto 13
    
    if(npot>0) then
       do i=1,npot
          read(15,*,IOSTAT=ist) junk; if(ist.ne.0) goto 14
       end do
    end if
    
    !this is the flag to print JTJ
    read(15,*,IOSTAT=ist) junk; if(ist.ne.0) goto 17
 
    read(15,*,IOSTAT=ist) njrows_out;  if(ist.ne.0) goto 15
    
    if(njrows_out>0) then
       allocate(jrows_out(njrows_out))
       do i=1,njrows_out
          read(15,*,IOSTAT=ist) jrows_out(i); if(ist.ne.0) goto 16
       end do
       print_jaco_rows = .true.
    end if
    
    close(15)
   
    
    
    return
    
10  continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: Cannot find the output options file: ',trim(outfile)
      close(51)
      write(*, *) 
      write(*, *) ' E4D: Cannot find the output options file: ',trim(outfile)
      return
      
11    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: The was a problem reading the first line in the output file: ',trim(outfile)
      close(51)
      write(*, *) 
      write(*, *) ' E4D: There was a problem reading the first line in the output file: ',trim(outfile)
      return

12    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: There was a problem reading the predicted data file name in: ',trim(outfile)
      close(51)
      write(*, *) 
      write(*, *) ' E4D: The was a problem reading the predicted data file in: ',trim(outfile)
      return

13    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' E4D: There was a problem reading the number of potential fields to write in: ',trim(outfile)
      close(51)
      write(*, *) 
      write(*, *) ' E4D: There was a problem reading the number of potential fields to write in: ',trim(outfile)
      return

14    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      close(51)
      write(*,*)
      write(*, *) ' E4D: There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      return 
      
15    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: There was a problem reading the number of rows of JTJ matrix outputs in: ',trim(outfile)
      close(51)
      write(*,*)
      write(*, *) ' E4D: There was a problem reading the number of rows of JTJ matrix outputs in: ',trim(outfile)
      return

16    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: There was a problem reading the row index of JTJ matrix for measurement: ',i,' in: ',trim(outfile)
      close(51)
      write(*, *)
      write(*, *) ' E4D: There was a problem reading the row index of JTJ matrix for measurement: ',i,' in: ',trim(outfile)
      return

17    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' E4D: There was a problem reading JTJ output option in: ',trim(outfile)
      write(51,*) ' E4D: Not printing Jacobian rows.'
      close(51)
      write(*, *)
      write(*, *) ' E4D: There was a problem reading JTJ output option in: ',trim(outfile)
      write(*, *) ' E4D: Not printing JTJ rows.'
      return
            
    end subroutine check_jrows_out
  !_______________________________________________________________________________________


  !_____________________________________________________________________
  subroutine send_commando(com)
    !!Send a general command to the slaves
    integer :: com
    
    call MPI_BCAST(com,1,MPI_INTEGER,0,E4D_COMM,ierr)
    
  end subroutine send_commando
  !____________________________________________________________________

 
 
end module output
 
