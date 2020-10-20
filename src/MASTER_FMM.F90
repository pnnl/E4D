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
    write(*,"(A32,F6.3,A8)") "  FMM: DONE WITH FORWARD RUN IN:",(tc-ts)/60," minutes"
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
    real :: ts,tc    

    call cpu_time(ts)
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

    call cpu_time(tc)
    write(*,"(A29,F5.1,A8)") " FMM: DONE BUILDING JACO IN:",(tc-ts)/60," minutes"

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

   !****************************************************************************
   !* validate joint inversion
   !****************************************************************************
   subroutine validate_jointInv_fmm
     implicit none

     ! local variables
     integer :: tag,lenchar,mnchar1,mnchar2
     integer ::  status(MPI_STATUS_SIZE)
     integer ::  jnk1, jnk2, jnk3, jnk4     
     integer :: tjnk1,tjnk2,tjnk3,tjnk4
     character(len=len(mshfile)) :: tmpmshfile

     
     ! STart
     if (.not.cgmin_flag(1).and..not.cgmin_flag(2)) then
        open(51,file='fmm.log',status='old',action='write',position='append')                     
        write(51,*) '======================================================================='
        write(51,*) '                  FMM: RUNNING IN SEPARATE INVERSION MODE ... '
        write(51,*) '======================================================================='
        close(51)
        write(*, *) '======================================================================='
        write(*, *) '                  FMM: RUNNING IN SEPARATE INVERSION MODE ... '
        write(*, *) '======================================================================='
     elseif (cgmin_flag(1).and..not.cgmin_flag(2)) then
        open(51,file='fmm.log',status='old',action='write',position='append')
        write(51,*) '=============================== WARNING ==============================='
        write(51,*) ' E4D inverse option file doesnot request joint inversion but'
        write(51,*) ' FMM inverse option file ',trim(invfile),' requests joint inversion!'
        write(51,*) ' NO JOINT INVERSION WILL BE PERFORMED.'
        write(51,*) '                  FMM: RUNNING IN SEPARATE INVERSION MODE ... '
        write(51,*) '======================================================================='
        close(51)
        write(*, *) '=============================== WARNING ==============================='
        write(*, *) ' E4D inverse option file doesnot request joint inversion but'
        write(*, *) ' FMM inverse option file ',trim(invfile),' requests joint inversion!'
        write(*, *) ' NO JOINT INVERSION WILL BE PERFORMED.'
        write(*, *) '                  FMM: RUNNING IN SEPARATE INVERSION MODE ... '
        write(*, *) '======================================================================='

        !reset cgmin_flag
        cgmin_flag = .false.
        
     elseif (.not.cgmin_flag(1).and.cgmin_flag(2)) then
        open(51,file='fmm.log',status='old',action='write',position='append')
        write(51,*) '=============================== WARNING ==============================='
        write(51,*) ' E4D inverse option file requests joint inversion but FMM inverse'
        write(51,*) ' option file ',trim(invfile),' doesnot request joint inversion!'
        write(51,*) ' NO JOINT INVERSION WILL BE PERFORMED.'
        write(51,*) '                  FMM: RUNNING IN SEPARATE INVERSION MODE ... '
        write(51,*) '======================================================================='
        close(51)
        write(*, *) '=============================== WARNING ==============================='
        write(*, *) ' E4D inverse option file requests joint inversion but FMM inverse'
        write(*, *) ' option file ',trim(invfile),' doesnot request joint inversion!'
        write(*, *) ' NO JOINT INVERSION WILL BE PERFORMED.'
        write(*, *) '                  FMM: RUNNING IN SEPARATE INVERSION MODE ... '
        write(*, *) '======================================================================='

        !reset cgmin_flag
        cgmin_flag = .false.
        
     elseif (cgmin_flag(1).and.cgmin_flag(2)) then
        open(51,file='fmm.log',status='old',action='write',position='append')                     
        write(51,*) '======================================================================='
        write(51,*) '                  FMM: RUNNING IN JOINT INVERSION MODE ... '
        write(51,*) '======================================================================='
        close(51)
        write(*, *) '======================================================================='
        write(*, *) '                  FMM: RUNNING IN JOINT INVERSION MODE ... '
        write(*, *) '======================================================================='

        tag = 0
        lenchar = len(mshfile)
        call MPI_SEND(mshfile,   lenchar,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_RECV(tmpmshfile,lenchar,MPI_CHARACTER,0,tag,MPI_COMM_WORLD,status,ierr)
     
        ! Check if both physics have same mesh files
        mnchar1=index(mshfile,'.')
        mnchar2=index(tmpmshfile,'.')

        if (mshfile(1:mnchar1) .ne. tmpmshfile(1:mnchar2)) then
           open(51,file='fmm.log',status='old',action='write',position='append')
           write(51,*) '------------------------------- WARNING -------------------------------'
           write(51,*) ' FMM mesh file name ',trim(mshfile),' is different from '
           write(51,*) ' E4D mesh file name ',trim(tmpmshfile)
           write(51,*) ' Checking if they both have the same number of nodes, elements, and faces.'
           write(51,*) '-----------------------------------------------------------------------'
           close(51)
           write(*, *) '------------------------------- WARNING -------------------------------'
           write(*, *) ' FMM mesh file name ',trim(mshfile),' is different from '
           write(*, *) ' E4D mesh file name ',trim(tmpmshfile)
           write(*, *) ' Checking if they both have the same number of nodes, elements, and faces.'
           write(*, *) '-----------------------------------------------------------------------'

           ! Checking number of nodes
           open(331,file=   mshfile(1:mnchar1+1)//".node",status="old",action="read")
           open(332,file=tmpmshfile(1:mnchar2+1)//".node",status="old",action="read")

           read(331,*)  jnk1, jnk2, jnk3, jnk4
           read(332,*) tjnk1,tjnk2,tjnk3,tjnk4

           close(331); close(332)
          
           if (jnk1 .ne. tjnk1) then
              open(51,file='fmm.log',status='old',action='write',position='append')
              write(51,*) '-------------------------------- ERROR --------------------------------'
              write(51,*) ' FMM node file ',mshfile(1:mnchar1+1)//".node",' has different number '
              write(51,*) ' of nodes than E4D node file ',tmpmshfile(1:mnchar2+1)//".node"
              write(51,*) ' FMM node file has ',jnk1,' nodes but E4D has ',tjnk1,' nodes.'
              write(51,*) ' FMM node file should have the same number of nodes as E4D node file  '
              write(51,*) ' for the cross-gradient joint inversion.'
              write(51,*) ' Aborting ...'
              write(51,*) '-----------------------------------------------------------------------'
              close(51)
              write(*, *) '-------------------------------- ERROR --------------------------------'
              write(*, *) ' FMM node file ',mshfile(1:mnchar1+1)//".node",' has different number '
              write(*, *) ' of nodes than E4D node file ',tmpmshfile(1:mnchar2+1)//".node"
              write(*, *) ' FMM node file has ',jnk1,' nodes but E4D has ',tjnk1,' nodes.'
              write(*, *) ' FMM node file should have the same number of nodes as E4D node file  '
              write(*, *) ' for the cross-gradient joint inversion.'
              write(*, *) ' Aborting ...'
              write(*, *) '-----------------------------------------------------------------------'
              call crash_exit_fmm
           else
              open(51,file='fmm.log',status='old',action='write',position='append')              
              write(51,*) ' FMM: Passed nodes - same number of nodes.'
              close(51)
              write(*, *) ' FMM: Passed nodes - same number of nodes.'
           endif

           ! checking number of elements
           open(331,file=   mshfile(1:mnchar1+1)//".ele",status="old",action="read")
           open(332,file=tmpmshfile(1:mnchar2+1)//".ele",status="old",action="read")

           read(331,*)  jnk1, jnk2, jnk3
           read(332,*) tjnk1,tjnk2,tjnk3

           close(331); close(332)
        
           if (jnk1 .ne. tjnk1) then
              open(51,file='fmm.log',status='old',action='write',position='append')
              write(51,*) '-------------------------------- ERROR --------------------------------'
              write(51,*) ' FMM element file ',mshfile(1:mnchar1+1)//".ele",' has different number '
              write(51,*) ' of elements than E4D element file ',tmpmshfile(1:mnchar2+1)//".ele"
              write(51,*) ' FMM element file has ',jnk1,' elements but E4D has ',tjnk1,' elements.'
              write(51,*) ' FMM element file should have the same number of elements as E4D element'
              write(51,*) ' file for the cross-gradient joint inversion.'           
              write(51,*) ' Aborting ...'
              write(51,*) '-----------------------------------------------------------------------'
              close(51)
              write(*, *) '-------------------------------- ERROR --------------------------------'
              write(*, *) ' FMM element file ',mshfile(1:mnchar1+1)//".ele",' has different number '
              write(*, *) ' of elements than E4D element file ',tmpmshfile(1:mnchar2+1)//".ele"
              write(*, *) ' FMM element file has ',jnk1,' elements but E4D has ',tjnk1,' elements.'
              write(*, *) ' FMM element file should have the same number of elements as E4D element'
              write(*, *) ' file for the cross-gradient joint inversion.'           
              write(*, *) ' Aborting ...'
              write(*, *) '-----------------------------------------------------------------------'
              call crash_exit_fmm
           else
              open(51,file='fmm.log',status='old',action='write',position='append')              
              write(51,*) ' FMM: Passed elements - same number of elements.'
              close(51)
              write(*, *) ' FMM: Passed elements - same number of elements.'
           endif

           open(331,file=   mshfile(1:mnchar1+1)//".face",status="old",action="read")
           open(332,file=tmpmshfile(1:mnchar2+1)//".face",status="old",action="read")

           read(331,*)  jnk1, jnk2
           read(332,*) tjnk1,tjnk2

           close(331); close(332)
        
           if (jnk1 .ne. tjnk1) then
              open(51,file='fmm.log',status='old',action='write',position='append')
              write(51,*) '-------------------------------- ERROR --------------------------------'
              write(51,*) ' FMM face file ',mshfile(1:mnchar1+1)//".face",' has different number '
              write(51,*) ' of faces than E4D face file ',tmpmshfile(1:mnchar2+1)//".face"
              write(51,*) ' FMM face file has ',jnk1,' faces but E4D has ',tjnk1,' faces.'
              write(51,*) ' FMM face file should have the same number of faces as E4D face file'
              write(51,*) ' for the cross-gradient joint inversion.'           
              write(51,*) ' Aborting ...'
              write(51,*) '-----------------------------------------------------------------------'
              close(51)
              write(*, *) '-------------------------------- ERROR --------------------------------'
              write(*, *) ' FMM face file ',mshfile(1:mnchar1+1)//".face",' has different number '
              write(*, *) ' of faces than E4D face file ',tmpmshfile(1:mnchar2+1)//".face"
              write(*, *) ' FMM face file has ',jnk1,' faces but E4D has ',tjnk1,' faces.'
              write(*, *) ' FMM face file should have the same number of faces as E4D face file'
              write(*, *) ' for the cross-gradient joint inversion.'  
              write(*, *) ' Aborting ...'
              write(*, *) '-----------------------------------------------------------------------'
              call crash_exit_fmm
           else
              open(51,file='fmm.log',status='old',action='write',position='append')              
              write(51,*) ' FMM: Passed faces - same number of faces.'
              close(51)
              write(*, *) ' FMM: Passed faces - same number of faces.'
           endif                  
           
        endif ! mesh files are different

     endif ! cgmin(1) and cgmin(2) -> true
             
   end subroutine validate_jointInv_fmm
     

end module master_fmm

