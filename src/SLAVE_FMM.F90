module slave_fmm

  use vars
  use fmm_vars
  use forward_fmm
  use assemble_fmm
  use jaco_fmm
  use slave
  
  contains
   
    !__________________________________________________________________
    subroutine go_slave_fmm
      implicit none
      integer :: command
    
      im_fmm = .true.
100   continue
      !Recieve a command from master
      call MPI_BCAST(command,1,MPI_INTEGER,0,FMM_COMM,ierr)

      select case(command)

      case(-1) 
         call go_slave
         goto 100

      !return to main
      case(0)
         return

      !recieve elements   
      case(1)
         call receive_dists_fmm
         goto 100
      
      case(2)
         call setup_frun_fmm
         goto 100

      case(3) 
!	 call cpu_time(Cstart)     
!         call build_A
!         call cpu_time(Cend)
!         my_abt = Cend-Cstart
         goto 100

      case(4) 
        call receive_slowness
        goto 100

     case(5) 
!        call cpu_time(Cstart)
!        call build_ksp
!        call cpu_time(Cend)
!        my_kspt = Cend-Cstart
        goto 100

     case(6) 
        call cpu_time(Cstart)
        call forward_run_fmm
        call cpu_time(Cend)
        my_frt_fmm = Cend-Cstart
        goto 100

     case(7) 
        call receive_info_fmm
        goto 100
        
     case(8) 
        call send_ttpred 
        goto 100

     case(10) 
        if(fresnel) then
           call build_jaco_fresnel
        else
           call build_jaco_raytrace  !see JACOBIAN_FMM.F90
        end if
        goto 100

     case(11)
!        call cpu_time(Cstart)
!        call build_jaco
!        call cpu_time(Cend)
!        my_jbt = Cend-Cstart
        goto 100
        
     case(12) 
!        call weight_jaco
        goto 100

     case(13)
!        call compute_pmatvec2
        goto 100
        
     case(14) 
!        call compute_pmatvec1
        goto 100

     case(15) 
        call MPI_SEND(my_frt,1,MPI_REAL,0,0,FMM_COMM,ierr)
        goto 100
        
     case(16)
        call MPI_SEND(my_abt,1,MPI_REAL,0,0,FMM_COMM,ierr)
        goto 100

     case(17) 
        call MPI_SEND(my_kspt,1,MPI_REAL,0,0,FMM_COMM,ierr)
        goto 100
        
     case(18) 
        call MPI_SEND(my_jbt,1,MPI_REAL,0,0,FMM_COMM,ierr)
        goto 100
        
     case(19)
!        call receive_rrseq
        goto 100

     case(20) 
!        call start_inv
        goto 100
     
     case(21)
!        call end_inv
        goto 100

     case(22)
        call build_sgroups_fmm
        goto 100

     case(23) 
        call send_tt
        goto 100

     case(24)
        write(*,*) "FMM_Sending"
        call send_sens
        goto 100

     case(25)
!        call send_Jcol
        goto 100

     case(26) 
        call receive_dists_fmm
        goto 100
        
     case(27)
        call build_rings !see FORWARD_FMM.F90
        goto 100

     case(51)
        call get_J_on_off
        goto 100

     case(52)
        call print_mysens
        goto 100
        
     case DEFAULT 
        goto 100
        
        
     end select
     
    end subroutine go_slave_fmm
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine receive_dists_fmm
      implicit none
      
      integer :: i
      if(allocated(jind)) deallocate(jind)
      if(allocated(sind)) deallocate(sind)
      if(allocated(s_pos)) deallocate(s_pos) 
      allocate(sind(n_rank_fmm-1,2))
      allocate(jind(n_rank_fmm-1,2))

      call MPI_BCAST(tns,1,MPI_INTEGER , 0, FMM_COMM, ierr )
      ns=tns
      allocate(s_pos(tns,3))
      call MPI_BCAST(s_pos,3*tns,MPI_REAL, 0, FMM_COMM, ierr)
      if(.not.fresnel) then
         call MPI_BCAST(nrc,1,MPI_INTEGER , 0, FMM_COMM, ierr )
         allocate(rc_pos(nrc,3))
         call MPI_BCAST(rc_pos,3*nrc,MPI_REAL , 0, FMM_COMM, ierr )
      end if
      call MPI_BCAST(jind,2*(n_rank_fmm-1),MPI_INTEGER , 0, FMM_COMM, ierr )
      call MPI_BCAST(sind,2*(n_rank_fmm-1),MPI_INTEGER , 0, FMM_COMM, ierr )
      if(fresnel) then
         allocate(frq(tns))
         call MPI_BCAST(frq,tns,MPI_real , 0, FMM_COMM, ierr )
      end if

      my_ns = sind(my_rank_fmm,2)-sind(my_rank_fmm,1)+1
      
    end subroutine receive_dists_fmm
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine receive_info_fmm
      implicit none
      integer :: i
      
      if(allocated(s_conf_fmm)) deallocate(s_conf_fmm)
      call MPI_BCAST(nm_fmm, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
      allocate(s_conf_fmm(nm_fmm,2))
      call MPI_BCAST(s_conf_fmm, 2*nm_fmm,MPI_INTEGER,0,FMM_COMM,ierr)
      call MPI_BCAST(fresnel,1,MPI_LOGICAL,0,FMM_COMM,ierr)

    end subroutine receive_info_fmm
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine setup_frun_fmm
      implicit none

      integer :: status(MPI_STATUS_SIZE)
      integer :: neven, nextra, ce, i,itmp
      integer :: beg,end
      logical :: eqf = .false.
      logical :: stat
      real, dimension(:), allocatable :: dist
      integer, dimension(1) :: tmp
      real :: mx,my,mz

      !receive nodes from master
      call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, FMM_COMM, ierr )
      allocate(nodes(nnodes,3),nbounds(nnodes))
      call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, FMM_COMM, ierr )
      call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, FMM_COMM, ierr )
      
    
      !receive elements from master
      call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
      allocate(elements(nelem,4),nbrs(nelem,4))
      call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,FMM_COMM,ierr)
      call MPI_BCAST(nbrs, nelem*4,MPI_INTEGER,0,FMM_COMM,ierr)

      !call MPI_BCAST(nfaces, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
      !allocate(faces(nfaces,4))
      allocate(zones(nelem))
      !call MPI_BCAST(faces, nfaces*4,MPI_INTEGER,0,FMM_COMM,ierr)
      call MPI_BCAST(zones, nelem,MPI_INTEGER,0,FMM_COMM,ierr)
     
      !receive the assignments
      call receive_dists_fmm
      
      !get the source and receiver nodes
      call mget_source_nodes(stat)

      !allocate and receive the use_ele vector
      allocate(use_ele(nelem))
      call MPI_BCAST(use_ele,nelem,MPI_LOGICAL,0,FMM_COMM,ierr)
      
      !allocate travel time array
      allocate(ttimes(nnodes,my_ns))
      ttimes = 0.0
 
      return

    end subroutine setup_frun_fmm
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine build_sgroups_fmm
      implicit none
      
      integer :: orig_group,new_group
      integer, dimension(n_rank_fmm-1) :: sranks
      integer :: i,j,k
      integer :: mysrank,nsrank
      call MPI_COMM_GROUP(FMM_COMM,orig_group,ierr)
      do i=1,n_rank_fmm-1
         sranks(i)=i
      end do
      call MPI_GROUP_INCL(orig_group,n_rank_fmm-1,sranks,new_group,ierr)
      call MPI_COMM_CREATE(FMM_COMM,new_group,SCOMM_FMM,ierr)
      call MPI_COMM_RANK(SCOMM_FMM,mysrank,ierr)
      call MPI_COMM_SIZE(SCOMM_FMM, nsrank, ierr)
      if(mysrank .ne. my_rank_fmm-1) then
         write(*,*) 'WARNING: my_srank is not equal to my_rank-1',my_rank_fmm
      end if

      
    end subroutine build_sgroups_fmm
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine receive_slowness
      implicit none
      
      call MPI_BCAST(nspd, 1, MPI_INTEGER, 0,FMM_COMM,ierr)
      if(.not.allocated(velocity)) allocate(velocity(nspd))
      call MPI_BCAST(velocity, nspd,MPI_REAL,0,FMM_COMM,ierr)

    end subroutine receive_slowness
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine send_ttpred
      implicit none
      integer :: i

      call fmm_assemble_data
      call MPI_SEND(nmy_drows,1,MPI_INTEGER,0,0,FMM_COMM, ierr)
      call MPI_SEND(my_drows,nmy_drows,MPI_INTEGER,0,0,FMM_COMM,ierr)
      call MPI_SEND(my_dvals,nmy_drows,MPI_REAL,0,0,FMM_COMM,ierr)
    
    end subroutine send_ttpred
    !__________________________________________________________________
 
    !__________________________________________________________________
     subroutine send_tt
       implicit none
       integer, dimension(2) :: spack
       integer :: es

       call MPI_BCAST(spack,2,MPI_INTEGER,0,FMM_COMM,ierr)

       if(spack(1)==my_rank_fmm) then
          es=spack(2)-sind(my_rank_fmm,1) + 1
          call MPI_SEND(ttimes(:,es),nnodes,MPI_REAL,0,0,FMM_COMM,ierr)
       end if

     end subroutine send_tt
     !__________________________________________________________________

     !__________________________________________________________________
     subroutine print_mysens
       implicit none
       integer :: i,j
       character*15 :: fname
       !if(allocated(Jaco)) then
          do i = 1,my_ns
             write(fname,'(A12,I1)') 'sensitivity.',sind(my_rank_fmm,1)+i-1
             write(*,*) 'printing sensitivities to ',fname,sum(Jaco(i,:))
             open(13,file=fname,status='replace',action='write')
             write(13,*) nelem,1
             do j=1,nelem
                write(13,*) Jaco(i,j)
             end do
             close(13)
          end do
       !end if
     end subroutine print_mysens
     !__________________________________________________________________
     
    

   
end module slave_fmm

