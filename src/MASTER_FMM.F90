module master_fmm

  use vars
  use fmm_vars
  use report_fmm
  use forward_fmm
  use input_fmm
  use output_fmm
  use master
contains

  !____________________________________________________________________
  subroutine send_info_fmm
    implicit none

    call send_command_fmm(7)
    call MPI_BCAST(nm_fmm, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
    call MPI_BCAST(s_conf_fmm, 2*nm_fmm,MPI_INTEGER,0,FMM_COMM,ierr)
    call MPI_BCAST(fresnel,1,MPI_LOGICAL,0,FMM_COMM,ierr)

  end subroutine send_info_fmm

  !_____________________________________________________________________
  !subroutine send_command_fmm(com)
  !  !!Send a general command to the slaves
  !  integer :: com
 !
  !     call MPI_BCAST(com,1,MPI_INTEGER,0,FMM_COMM,ierr)

 ! end subroutine send_command_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_dists_fmm
    !read inputs and distribute the run info to the slave
    implicit none
    integer :: neven,nextra,ce,i,j,k
    integer :: npos,nneg,np_ranks,nn_ranks,tn_rank
    integer, dimension(n_rank_fmm-1) :: tmp

    
    !!divide the slave processors given sources
    tns=ns
    tn_rank=n_rank_fmm
   
    if((n_rank_fmm-1)>tns) then
       tn_rank=n_rank_fmm
       n_rank_fmm=tns+1
    end if

!   nneg = tns 
    nn_ranks=n_rank_fmm-1
    
   
    if(allocated(sind)) deallocate(sind)
    allocate(sind(tn_rank-1,2))
    sind=0
  
    if(tns > 0) then
       neven = tns/nn_ranks
       nextra = tns - neven*nn_ranks
       if(neven>0) then
          ce = 1
          do i=1,nn_ranks-nextra
             sind(i,1) = ce
             sind(i,2) = ce+neven-1
             ce=ce+neven
          end do
          
          if(nextra>0) then
             do i=nn_ranks+1-nextra,nn_ranks
                sind(i,1)=ce
                sind(i,2)=ce+neven
                ce=ce+neven+1
             end do
          end if
       else
          
          do i=1,nn_ranks
             sind(i,1) = i
             sind(i,2) = i
          end do
       end if

    end if
    
    n_rank_fmm = tn_rank
    
    !!Build jaco assignments 
    if(allocated(jind)) deallocate(jind)
    allocate(jind(n_rank_fmm-1,2))
    jind = 0
    tmp = 0

    if(fresnel) then
       if(nm_fmm > 0) then
          neven = nm_fmm/nn_ranks
          nextra = nm_fmm - neven*nn_ranks
          if(neven>0) then
             ce = 1
             do i=1,nn_ranks-nextra
                jind(i,1) = ce
                jind(i,2) = ce+neven-1
                ce=ce+neven
             end do
             
             if(nextra>0) then
                do i=nn_ranks+1-nextra,nn_ranks
                   jind(i,1)=ce
                   jind(i,2)=ce+neven
                   ce=ce+neven+1
                end do
             end if
          else
          
             do i=1,nn_ranks
                jind(i,1) = i
                jind(i,2) = i
             end do
          end if
          
       end if
    else
       do i=1,nm_fmm
          do j=1,n_rank_fmm - 1
             if(s_conf_fmm(i,1).ge.sind(j,1) .and. s_conf_fmm(i,1).le.sind(j,2)) then
                tmp(j)=tmp(j)+1
                exit
             end if
          end do
       end do
       jind(1,1) = 1
       jind(1,2) = tmp(1)
       do i=2,n_rank_fmm -1
          jind(i,1) = jind(i-1,2) + 1
          jind(i,2) = jind(i-1,2) + tmp(i)
       end do
    end if

    call MPI_BCAST(tns,1,MPI_INTEGER , 0, FMM_COMM, ierr )
    call MPI_BCAST(s_pos,3*tns,MPI_REAL, 0, FMM_COMM, ierr)
    if(.not.fresnel) then
       call MPI_BCAST(nrc,1,MPI_INTEGER , 0, FMM_COMM, ierr )
       call MPI_BCAST(rc_pos,3*nrc,MPI_REAL , 0, FMM_COMM, ierr )
    end if
    call MPI_BCAST(jind,2*(n_rank_fmm-1),MPI_INTEGER , 0, FMM_COMM, ierr )
    call MPI_BCAST(sind,2*(n_rank_fmm-1),MPI_INTEGER , 0, FMM_COMM, ierr )
    if(fresnel) then
       call MPI_BCAST(frq,tns,MPI_real , 0, FMM_COMM, ierr )
    end if
    call nreport_fmm(64)
  end subroutine send_dists_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine setup_forward_fmm
    implicit none

    integer :: npre,i,nchr,dim,bflag,itmp
    real :: jnk,jnk1 
    integer :: status(MPI_STATUS_SIZE)
    logical :: stat
    
   
    !read the node info
    call nreport_fmm(60)
    call read_nodes_fmm(stat)
    if(.not.stat) call crash_exit_fmm
    call nreport_fmm(61)
   
    !!read the element info 
    call read_elements_fmm(stat)
    if(.not.stat) call crash_exit_fmm
    call nreport_fmm(62)
   
    if(nspd .ne. nelem ) then
       call nreport_fmm(57)
       call crash_exit_fmm
    end if

    !!read the neighbors
    call read_neighbors_fmm(stat)
    if(.not.stat) call crash_exit_fmm
    
    !command slaves to do the setup
    call send_command_fmm(2)

    !send everything to slave
    call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, FMM_COMM, ierr )
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, FMM_COMM, ierr )
    call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, FMM_COMM, ierr )
    call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,FMM_COMM,ierr)
    !call MPI_BCAST(nfaces, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
    call MPI_BCAST(nbrs, nelem*4,MPI_INTEGER,0,FMM_COMM,ierr)
    call MPI_BCAST(zones, nelem,MPI_INTEGER,0,FMM_COMM,ierr)
    
    !get/check the source node positions
    call mget_source_nodes(stat)
    if(.not. stat) call crash_exit_fmm
    
    !send the assignments
    call send_dists_fmm 

    !build and send the use_ele vector
    call build_use_ele
    !write(*,*) use_ele
    !write(*,*) nelem
    call MPI_BCAST(use_ele,nelem,MPI_LOGICAL,0,FMM_COMM,ierr)

    !instruct the slave to build the first ring map
    call send_command_fmm(27)

  end subroutine setup_forward_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine run_forward_fmm
    implicit none
    integer :: i,ierr,nmin,nmax,ntot,itm
    real :: rtm,ts,tc,lrep,tmax,tmin,ttot
    real, dimension(2) :: pck
    integer ::  status(MPI_STATUS_SIZE)
    
    call send_command_fmm(6)

    return

    call cpu_time(ts)
    tmin=1e9
    tmax=0
    ttot=0
    
    nmin=1e9
    nmax=0
    ntot=0
    lrep=ts
    do i=1,tns
       call MPI_RECV(pck,2,MPI_REAL,MPI_ANY_SOURCE,1,FMM_COMM,status,ierr)
   
       rtm=pck(1)
       itm=int(pck(2))
       if(rtm<tmin) tmin = rtm
       if(rtm>tmax) tmax = rtm
       ttot = ttot+rtm

       if(itm<nmin) nmin=itm
       if(itm>nmax) nmax=itm
       ntot = ntot+itm

       call cpu_time(tc)
       if((tc-lrep)>30) then
          write(*,"(A15,F5.2,A11,F6.3,A8)") "  FORWARD RUN: ",100*real(i)/real(tns),"% DONE IN :",(tc-ts)/60," MINUTES"
          write(*,*) "   MIN / MAX / AVE  RUN TIMES (sec.)",tmin,"/",tmax,"/",ttot/i 
          write(*,*) "   MIN / MAX / AVE  ITERATIONS      ",nmin,"/",nmax,"/",ntot/i
          lrep=tc
       end if
      
    end do
    call cpu_time(tc)
    write(*,"(A27,F6.3,A8)") "  DONE WITH FORWARD RUN IN:",(tc-ts)/60," minutes"
    write(*,*) 
  end subroutine run_forward_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine read_nodes_fmm(st)
       implicit none
       logical :: st
       logical :: exst
       integer :: nchr,npre
       integer :: ist,bflag,dim,jnk,jnk1
       integer :: i
       logical, dimension(:), allocatable :: nbi

       st = .true.
       if(allocated(nodes)) return
     
       !get the meshfile prefix
       nchr=len_trim(mshfile)
       do i=1,nchr
          if(mshfile(i:i)=='.') then
             npre=i+1;
             exit
          end if
       end do
          
       inquire(file=mshfile(1:npre)//".node",exist=exst)
       if(.not.exst) goto 10

       open(10,file=mshfile(1:npre)//".node",status="old",action="read")
       read(10,*,IOSTAT=ist) nnodes,dim,jnk,bflag
       if(ist .ne. 0) goto 11
       
       allocate(nodes(nnodes,3),nbounds(nnodes))
       do i=1,nnodes
          read(10,*,IOSTAT=ist) jnk,nodes(i,1:3),jnk1,nbounds(i)
          if(ist .ne. 0) goto 12
       end do
       close(10)

       allocate(nbi(nnodes))
       nbi=.false.
       do i=1,nnodes
          if(nbounds(i)<0) then
             nbi(abs(nbounds(i)))=.true.
          end if
       end do
       n_met = 0
       do i=1,nnodes
          if(nbi(i)) n_met=n_met+1
       end do

       return

10     continue
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*)
       write(51,*) ' Cannot find the node file : ',mshfile(1:npre)//'.node'
       write(51,*) ' Aborting ...'
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the node file : ',mshfile(1:npre)//'.node'
       write(*,*) ' Aborting ...'
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*)
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the node file: ',mshfile(1:npre)//'.node'
       write(51,*) ' Aborting ...'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the node file: ',mshfile(1:npre)//'.node'
       st=.false.
       return

12     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the node file: ',mshfile(1:npre)//'.node'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the node file: ',mshfile(1:npre)//'.node'
       st=.false.
       return

13     continue
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) ' Node number ',i,' has a negative boundary flag'
       write(51,*) ' which indicates a infinite conductivity boundary.'
       write(51,*) ' Infinite conductivity boundaries are not '
       write(51,*) ' implemented in this version of fmm.'
       write(51,*) ' Aborting ...'
       close(51) 
       write(*,*) 
       write(*,*) ' Node number ',i,' has a negative boundary flag'
       write(*,*) ' which indicates a infinite conductivity boundary.'
       write(*,*) ' Infinite conductivity boundaries are not '
       write(*,*) ' implemented in this version of fmm.'
       write(*,*) ' Aborting ...'
       st = .false.
       return
     end subroutine read_nodes_fmm
     !__________________________________________________________________________

     !__________________________________________________________________________
     subroutine read_elements_fmm(st)
       implicit none
       logical :: st
       logical :: exst
       integer :: i,nchr,npre
       integer :: ist,dim,jnk,nzn
       logical, dimension(:), allocatable :: zne
 
       if(allocated(elements)) return

       st = .true.
       !get the meshfile prefix
       nchr=len_trim(mshfile)
       do i=1,nchr
          if(mshfile(i:i)=='.') then
             npre=i+1;
             exit
          end if
       end do
          
       inquire(file=mshfile(1:npre)//".ele",exist=exst)
       if(.not.exst) goto 10

       open(10,file=mshfile(1:npre)//".ele",status="old",action="read")
       read(10,*,IOSTAT=ist) nelem,dim,jnk
       if(ist .ne. 0) goto 11

       allocate(elements(nelem,4),zones(nelem))
       do i=1,nelem
          read(10,*,IOSTAT=ist) jnk,elements(i,1:4),zones(i)
          if(ist .ne. 0) goto 12
       end do
       close(10)
       
       allocate(zne(nelem))
       zne=.false.
       do i=1,nelem
          if(zones(i) .le. 0) then
             nrz = i
             call nreport_fmm(63)
             st=.false.
             return
          end if
          zne(zones(i))=.true.
       end do
       nrz=0
       do i=1,nelem
          if(zne(i)) nrz=nrz+1
       end do
       deallocate(zne)
       return

10     continue
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' Cannot find the element file : ',mshfile(1:npre)//'.ele'
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the ele file : ',mshfile(1:npre)//'.ele'
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       close(51)
       write(*,*)
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       st=.false.
       return

12     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       st=.false.
       return

     end subroutine read_elements_fmm
     !__________________________________________________________________________

     !__________________________________________________________________________
     subroutine read_faces_fmm(st)
     implicit none
     logical :: st
     logical :: exst
     integer :: npre,nchr,i,jnk,ist

     st = .true.
     !get the meshfile prefix
     nchr=len_trim(mshfile)
     do i=1,nchr
        if(mshfile(i:i)=='.') then
           npre=i+1;
           exit
        end if
     end do
     
     inquire(file=mshfile(1:npre)//".face",exist=exst)
     if(.not.exst) goto 10

     open(10,file=mshfile(1:npre)//".face",status="old",action="read")
     
     read(10,*,IOSTAT=ist) nfaces,jnk
     if(ist .ne. 0) goto 11

     allocate(faces(nfaces,4))
     do i=1,nfaces
        read(10,*,IOSTAT=ist) jnk,faces(i,1:4)
        if(ist .ne. 0) goto 12
     end do
     close(10)
     return

10   continue
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' Cannot find the face file : ',mshfile(1:npre)//'.face'
       write(51,*) " Aborting ..."
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the face file : ',mshfile(1:npre)//'.face'
       write(*,*) " Aborting ..."
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write')
       write(51,*) 
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the face file: ',mshfile(1:npre)//'.face'
       write(*,*) " Aborting ..."
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the face file: ',mshfile(1:npre)//'.face'
       write(*,*) " Aborting"
       st=.false.
       return

12     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the face file: ',mshfile(1:npre)//'.face'
       close(51)
       write(*,*)
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the face file: ',mshfile(1:npre)//'.face'
       st=.false.
       return

     end subroutine read_faces_fmm
   !__________________________________________________________________________
   
  !__________________________________________________________________________
  subroutine read_neighbors_fmm(st)
    implicit none

    logical :: st
    logical :: exst
    integer :: i,nchr,npre
    integer :: ist,dim,jnk,nzn,nneigh
    logical, dimension(:), allocatable :: zne
    st = .true.
  
    !get the meshfile prefix
    nchr=len_trim(mshfile)
    do i=1,nchr
       if(mshfile(i:i)=='.') then
          npre=i+1;
          exit
       end if
    end do
    
    inquire(file=mshfile(1:npre)//".neigh",exist=exst)
    if(.not.exst) goto 10

    open(10,file=mshfile(1:npre)//".neigh",status="old",action="read")
    read(10,*,IOSTAT=ist) nneigh,jnk
    if(ist .ne. 0 .or. nneigh .ne. nelem) goto 11

    allocate(nbrs(nelem,4))
    do i=1,nelem
       read(10,*,IOSTAT=ist) jnk,nbrs(i,1:4)
       if(ist .ne. 0) goto 12
    end do
    close(10)
       
    return

10     continue
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' Cannot find the neighbor file : ',mshfile(1:npre)//'.neigh'
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the neighbor file : ',mshfile(1:npre)//'.neigh'
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the neighbor file: ',mshfile(1:npre)//'.neigh'
       close(51)
       write(*,*)
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the neighbor file: ',mshfile(1:npre)//'.neigh'
       st=.false.
       return

12     continue
       close(10)
       open(51,file='fmm.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the neighbor file: ',mshfile(1:npre)//'.neigh'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the neighbor file: ',mshfile(1:npre)//'.neigh'
       st=.false.
       return
    
     end subroutine read_neighbors_fmm
  !__________________________________________________________________________
  !____________________________________________________________________
   subroutine send_slowness
    implicit none

    call send_command_fmm(4)
    call MPI_BCAST(nspd, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
    call MPI_BCAST(speed, nspd,MPI_REAL,0,FMM_COMM,ierr)

   end subroutine send_slowness
  !____________________________________________________________________
  !____________________________________________________________________
  !get travel time
  subroutine get_ttpred
    implicit none
    integer :: opt
    integer :: i,j
    integer :: nadd
    integer :: nbuff
    integer, dimension(nm_fmm*2) :: ibuff
    real, dimension(nm_fmm) :: rbuff
    integer ::  status(MPI_STATUS_SIZE)
  
    !instruct slave to assemble and send the prediceted data
    call send_command_fmm(8)
 
    !allocate ttpred if not already done and zero
    if(.not. allocated(ttpred)) then
       allocate(ttpred(nm_fmm),dpred(nm_fmm))
    end if
    ttpred = 0.0
    rbuff = 0.0
    ibuff = 1
    do i=1,n_rank_fmm-1
       call MPI_RECV(nbuff,1,MPI_INTEGER,i,0,FMM_COMM,status,ierr)
       call MPI_RECV(ibuff(1:nbuff),nbuff,MPI_INTEGER,i,0,FMM_COMM,status,ierr)
       call MPI_RECV(rbuff(1:nbuff),nbuff,MPI_REAL,i,0,FMM_COMM,status,ierr)
       
       do j=1,nbuff
          ttpred(ibuff(j))=ttpred(ibuff(j))+rbuff(j)
       end do
          
    end do
    dpred = ttpred
    call output_ttpred
  end subroutine get_ttpred
  !____________________________________________________________________
  !____________________________________________________________________
   subroutine build_sgroupm_fmm
    implicit none
    integer :: orig_group,new_group
    integer , dimension(n_rank_fmm-1) :: sranks
    integer :: i
    call send_command_fmm(22)

    !setup a new mpi group including only the slaves for forward runs
    !note this is only done in master because all processors in 
    !E4D_COMM must participate
    call MPI_COMM_GROUP(FMM_COMM,orig_group,ierr)
    do i=1,n_rank_fmm-1
       sranks(i)=i
    end do
    call MPI_GROUP_INCL(orig_group,n_rank_fmm-1,sranks,new_group,ierr)
    call MPI_COMM_CREATE(FMM_COMM,new_group,SCOMM_FMM,ierr)

  end subroutine build_sgroupm_fmm
  !____________________________________________________________________
  
!!$  !____________________________________________________________________
!!$  subroutine setup_inv_fmm
!!$    implicit none
!!$    integer :: status(MPI_STATUS_SIZE)
!!$    integer :: i,n
!!$    call send_command_fmm(10)
!!$    
!!$    do i=1,n_rank_fmm-1
!!$         call MPI_RECV(n,1,MPI_INTEGER,i,0,FMM_COMM,status,ierr)
!!$         if(i .eq. 1) then
!!$            jind(i,1) = 1
!!$            jind(i,2) = n
!!$         else
!!$            jind(i,1) = jind(i-1,2) + 1
!!$            jind(i,2) = jind(i-1,2) + n
!!$         end if
!!$     end do
!!$     
!!$    return
!!$  end subroutine setup_inv_fmm
!!$  !____________________________________________________________________

  !____________________________________________________________________
  subroutine make_jaco_fmm
    implicit none
    integer :: status(MPI_STATUS_SIZE)
    integer :: i,n,ntot
    integer, dimension(nm_fmm) :: failed,tmp

    call send_command_fmm(10)
    if(.not.fresnel) then
       ntot = 0
       do i=1,n_rank_fmm-1
          
          call MPI_RECV(n,1,MPI_INTEGER,i,0,FMM_COMM,status,ierr)
          call MPI_RECV(tmp(1:n),n,MPI_INTEGER,i,i,FMM_COMM,status,ierr)
          failed(ntot+1:ntot+n) = tmp(1:n)
          ntot=ntot+n
          
       end do
       open(67,file='fmm.log',status='old',action='write',position='append')  
       write(67,"(A3,I5.5,A53)") '   ',ntot," Rays failed to trace back to the source "
       close(67)
       return
       
       do i=1,ntot
          write(67,"(A10,I5.5)") "         ",failed(i) 
       end do
       write(67,*)
       close(67)

    else


    end if

  end subroutine make_jaco_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine print_sens_fmm
    call send_command(52)  
  end subroutine print_sens_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine build_use_ele
    implicit none
    integer :: i,j
    allocate(use_ele(nelem))
    use_ele = .true.

    !!if the number of zones used in the simulation is not set to zeros
    !!the find the zones we're not using
    if(nzf .ne. 0) then
       use_ele = .false.
       do i=1,nelem
          do j=1,nzf
             if(zones(i) .eq. zsims(j)) then
                use_ele(i) = .true.
                exit
             end if
          end do
       end do
    end if
  end subroutine build_use_ele
  !____________________________________________________________________

end module master_fmm

