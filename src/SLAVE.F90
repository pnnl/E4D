module slave

  use vars
  use fmm_vars
  use build_amap
  use forward
  use assemble
  use jacob
 
#ifdef imimode
  use imi
#endif
#ifdef resmode
  use dd_opt_slave
#endif
  
  contains
   
    !__________________________________________________________________
    subroutine go_slave
      implicit none
      integer :: command
      integer :: COMM
      if(im_fmm) then
         COMM = FMM_COMM
      else
         COMM = E4D_COMM
      end if
     
100   continue 
      !Recieve a command from master
      call MPI_BCAST(command,1,MPI_INTEGER,0,COMM,ierr)

      select case(command)
      !return to main
      case(0)
         return

      !recieve elements   
      case(1)
         call receive_dists
         goto 100
      
      case(2)
         call setup_frun
         goto 100

      case(3) 
	 call cpu_time(Cstart)     
         call build_A
         call cpu_time(Cend)
         my_abt = Cend-Cstart
         goto 100

      case(4) 
        call receive_sigma
        goto 100

     case(5) 
        call cpu_time(Cstart)
        call build_ksp
        call cpu_time(Cend)
        my_kspt = Cend-Cstart
        goto 100

     case(6) 
        call cpu_time(Cstart)
        call forward_run
        call cpu_time(Cend)
        my_frt = Cend-Cstart
        goto 100

     case(7) 
        call receive_info
        goto 100
        
     case(8) 
        call send_dpred(1) 
        goto 100

     case(10) 
        call report_nind
        goto 100

     case(11)
        call cpu_time(Cstart)
        call build_jaco
        call cpu_time(Cend)
        my_jbt = Cend-Cstart
        goto 100
        
     case(12) 
        call weight_jaco
        goto 100

     case(13)
        call compute_pmatvec2
        goto 100
        
     case(14) 
        call compute_pmatvec1
        goto 100

     case(15) 
        call MPI_SEND(my_frt,1,MPI_REAL,0,0,E4D_COMM,ierr)
        goto 100
        
     case(16)
        call MPI_SEND(my_abt,1,MPI_REAL,0,0,E4D_COMM,ierr)
        goto 100

     case(17) 
        call MPI_SEND(my_kspt,1,MPI_REAL,0,0,E4D_COMM,ierr)
        goto 100
        
     case(18) 
        call MPI_SEND(my_jbt,1,MPI_REAL,0,0,E4D_COMM,ierr)
        goto 100
        
     case(19)
        call receive_rrseq
        goto 100

     case(20) 
        call start_inv
        goto 100
     
     case(21)
        call end_inv
        goto 100

     case(22)
        call build_sgroups
        goto 100

     case(23) 
        call send_pot
        goto 100

     case(24)
        call send_sens
        goto 100

     case(25)
        call send_Jcol
        goto 100

     case(26) 
        call receive_dists
        goto 100
        
     case(27) 
        call rec_excon
        goto 100

     case(50)
        call write_myjaco
        goto 100

     case(51) 
        call get_J_on_off
        goto 100

     case(52)
        call send_full_jaco
        goto 100
   
     case(102) 
        call setup_fruni
        goto 100
        
     case(103) 
        call cpu_time(Cstart)     
        call build_Ai
        call cpu_time(Cend)
        my_abt = Cend-Cstart
        goto 100

     case(104) 
        call receive_sigmai
        goto 100
        
     case(106) 
        call cpu_time(Cstart)
        call forward_runi
        call cpu_time(Cend)
        my_frt = Cend-Cstart
        goto 100
        
     case(108)
        call send_dpred(2) 
        goto 100

     case(111) 
        call cpu_time(Cstart)
        call cpu_time(Cend)
        my_jbt = Cend-Cstart
        goto 100

     case(112) 
        goto 100
  
     case(113)
        call compute_pmatvec2i
        goto 100

     case(114)
        call compute_pmatvec1i
        goto 100

     case(120) 
        call start_invi
        goto 100

     case(121) 
        call end_invi
        goto 100

     case(123) 
        call send_poti
        goto 100

     case(124)
        call slave_print_jaco_rows
        goto 100
        
     case(213) 
        call compute_pmatvec2_dbl
        goto 100
        
     case(214) 
        call compute_pmatvec1_dbl
        goto 100
 
     case(215) 
        invi = .false.
        goto 100

     case(216) 
        invi = .true.
        goto 100

     case(1000) 
        call recieve_dd_opts
        goto 100

     case(1001)
        call ddopt_gather_poles(1)
        goto 100

     case(1002)
        call build_dd_st_en
        goto 100
        
    case(1003)
        call ddopt_gather_poles(0)
        goto 100

     case(1004)
        call ddopt_get_next_set
        goto 100

     case(1005)
        call send_next_best
        goto 100

     case (1006)
        call receive_new_jinds
        goto 100

     case(1007) 
        call rebuild_dobs_dpred_Wd
        goto 100

     case(1008)
        call ddopt_receive_poles
        goto 100
  
     case(1009)
        call sbuild_inverse_dist
        goto 100
     
     case DEFAULT 
        goto 100
        
        
     end select
     
    end subroutine go_slave
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine receive_dists
      implicit none
      
      integer :: i
      if(allocated(jind)) deallocate(jind)
      if(allocated(eind)) deallocate(eind)
      if(allocated(e_pos)) deallocate(e_pos)
      allocate(eind(n_rank-1,2))
      allocate(jind(n_rank-1,2))

      call MPI_BCAST(tne, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
      allocate(e_pos(tne,4))
      call MPI_BCAST(e_pos,4*tne,MPI_REAL,0,E4D_COMM,ierr)
      call MPI_BCAST(jind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
      call MPI_BCAST(eind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )

      my_ne = eind(my_rank,2)-eind(my_rank,1)+1
      
    end subroutine receive_dists
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine receive_info
      implicit none
      integer :: i
      
      if(allocated(s_conf)) deallocate(s_conf)
      if(allocated(ms_conf)) deallocate(ms_conf)
      if(allocated(ms_currents)) deallocate(ms_currents)
      call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
      call MPI_BCAST(ms_flag, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
      
      if (ms_flag==1) then
		allocate(ms_conf(nm,6))
		allocate(ms_currents(nm,2))
		call MPI_BCAST(ms_conf, 6*nm,MPI_INTEGER,0,E4D_COMM,ierr)
		call MPI_BCAST(ms_currents, 2*nm,MPI_REAL,0,E4D_COMM,ierr)
	  else
	    allocate(s_conf(nm,4))
	    call MPI_BCAST(s_conf, 4*nm,MPI_INTEGER,0,E4D_COMM,ierr)
      end if
  
    end subroutine receive_info
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine rec_excon
      implicit none
     
      call MPI_BCAST(nec, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
      if(nec .ne. 0) then
         allocate(ex_cols(nec))  
         call MPI_BCAST(ex_cols,nec,MPI_INTEGER,0,E4D_COMM,ierr)
      end if
    
    end subroutine rec_excon
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine setup_frun
      implicit none

      integer :: status(MPI_STATUS_SIZE)
      integer :: neven, nextra, ce, i,itmp
      integer :: beg,end
      logical :: eqf = .false.
      real, dimension(:), allocatable :: dist
      integer, dimension(1) :: tmp
      real :: mx,my,mz

      !receive nodes from master
      call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
      allocate(nodes(nnodes,3),nbounds(nnodes))
      call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
      call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, E4D_COMM, ierr )
      
    
      !receive elements from master
      call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
      allocate(elements(nelem,4))
      call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)

      call MPI_BCAST(nfaces, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
      allocate(faces(nfaces,4))
      allocate(zones(nelem))
      call MPI_BCAST(faces, nfaces*4,MPI_INTEGER,0,E4D_COMM,ierr)
      call MPI_BCAST(zones, nelem,MPI_INTEGER,0,E4D_COMM,ierr)
     
      !receive the tank flag, choose zero potential node if true
      call MPI_BCAST(itmp,1,MPI_INTEGER,0,E4D_COMM,ierr)
      if(itmp==1) then
         tank_flag = .true.
!!$         allocate(dist(nnodes))
!!$         mx=sum(nodes(:,1))/nnodes;!mx=0
!!$         my=sum(nodes(:,2))/nnodes;!my=0
!!$         mz=sum(nodes(:,3))/nnodes;!mz=0
!!$         
!!$         dist=sqrt((mx-nodes(:,1))**2 + (my-nodes(:,2))**2 + (mz-nodes(:,3))**2)
!!$         !tmp = maxloc(dist)
!!$         tmp=minloc(dist)
!!$         !i_zpot=tmp(1)
!!$         deallocate(dist)
      end if
      

     
      !receive the assignments
      call receive_dists
      
     
      !build A_map and delA.......................................
      call build_delA
     
       !!Initialize the Petsc A matrix
      call MatCreateSeqAIJ(PETSC_COMM_SELF,nnodes,nnodes,d_nz,d_nnz,A,perr)
      call MatSetFromOptions(A,perr)
    
        
      
      !Get the electrode ownership indexes 
      call get_electrode_nodes
      call MatGetType(A,tp,perr)
    
  
      !Set up the source and solution vectors
      call VecCreate(PETSC_COMM_SELF,X,perr)
      call VecSetSizes(X,nnodes,PETSC_DECIDE,perr)
      call VecSetFromOptions(X,perr)

      !allocate and setup the pole solution vector
      call VecCreate(PETSC_COMM_SELF,psol,perr)
      call VecSetSizes(psol,nnodes,PETSC_DECIDE,perr)
      call VecSetFromOptions(psol,perr)
      allocate(poles(nnodes,my_ne))
      poles = 0
     
      call VecCreate(PETSC_COMM_SELF,B,perr)
      call VecSetSizes(B,nnodes,PETSC_DECIDE,perr)
      call VecSetFromOptions(B,perr)
    
      if(.not.eqf) neq = 0
          
#ifdef imimode
      !A_map and S_map given the coupling info so we can deallocate the 
      !nodes and elements
      if(minval(nbounds)<0) then
          call setup_infrastructure  !see FORWARD.F90
      end if
#endif
     
      deallocate(nodes)
      deallocate(elements)
      deallocate(faces)
      deallocate(zones)
    end subroutine setup_frun
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine setup_fruni
      implicit none
      call MatCreateSeqAIJ(PETSC_COMM_SELF,nnodes,nnodes,d_nz,d_nnz,Ai,perr)
      call MatSetFromOptions(Ai,perr)
      allocate(polesi(nnodes,my_ne))
      polesi=0
      allocate(Jdp(nnodes))
      Jdp=0
      call set_iflg(.true.)
    end subroutine setup_fruni
    !__________________________________________________________________

    
    !__________________________________________________________________
    subroutine report_nind
      !send the node assignments set by PETSC to master
      call MPI_SEND(nind(my_rank,1:2),2,MPI_INTEGER,0,0,E4D_COMM,ierr)
    end subroutine report_nind
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine build_sgroups
      implicit none
      
      integer :: orig_group,new_group
      integer, dimension(n_rank-1) :: sranks
      integer :: i,j,k
      call MPI_COMM_GROUP(E4D_COMM,orig_group,ierr)
      do i=1,n_rank-1
         sranks(i)=i
      end do
      call MPI_GROUP_INCL(orig_group,n_rank-1,sranks,new_group,ierr)
      call MPI_COMM_CREATE(E4D_COMM,new_group,SCOMM,ierr)
      call MPI_COMM_RANK(SCOMM,my_srank,ierr)
      call MPI_COMM_SIZE(SCOMM, n_srank, ierr)
      if(my_srank .ne. my_rank-1) then
         write(*,*) 'WARNING: my_srank is not equal to my_rank-1',my_rank
      end if

      
    end subroutine build_sgroups
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine build_A
      implicit none

      integer, dimension(nnodes) :: ncolss
      integer, dimension(50) :: colss
      integer, dimension(neq) :: llink,flink

      integer :: i,lrnl,lrnu,row,col,rbv,cbv,j,bnum,enum
      logical :: ilow,iup, eflag,check

      integer :: rw   
      integer :: ncls 
      integer :: cls(50) 
      
      !check to see if I'm computing a infrastructure source or not
      enum=eind(my_rank,1) 
      bnum=nbounds(e_nods(enum)) 
     
      !zero A
      call MatZeroEntries(A,perr)
      
      !if(.not. tank_flag) then
         do i=1,10*nelem
            row=rows(A_map(i))
            col=cols(A_map(i))
            rbv=nbounds(row)
            cbv=nbounds(col)
            
            eflag = .true.
            if(bnum<0 .and. (rbv==bnum .or. cbv==bnum)) eflag = .false.
         
            !lower triangle
            if(((rbv .eq. 2) .or. (cbv .eq. 2)) .and. .not. tank_flag) then
               !one or both nodes are on an outerboundary so set to zero for bc's
               val(1) = 0
#ifdef imimode
            elseif( (rbv<0 .or. cbv < 0) .and. eflag ) then
               !one or both nodes are one an ifrastructure boundary (and this is for an electrode pole)
               val(1) = 0
#endif
            else
               
               val(1) = sigma(S_map(i))*delA(i)
            end if
            
            prn(1) = row-1
            pcn(1) = col-1
            
            call MatSetValues(A,1,prn,1,pcn,val,ADD_VALUES,perr)          
            !upper triangle
            if(row .ne. col) then
               call MatSetValues(A,1,pcn,1,prn,val,ADD_VALUES,perr)
            end if
            
         end do
        
            
         !!fill in the diagonal for the zero potential bc's
         do i=1,nnodes
            !I own this node
          
            if((nbounds(i).eq.2) .and. .not. tank_flag) then
               !this node is on a boundary, set val to 1.0 for bc
               prn(1)=i-1;
               pcn(1)=i-1;
               val(1)=1;
               call MatSetValues(A,1,prn,1,pcn,val,ADD_VALUES,perr)
#ifdef imimode
            elseif(nbounds(i) < 0 .and. nbounds(i) .ne. bnum ) then
               !one or both nodes are on an ifrastructure boundary (and this is for and electrode pole)
               prn(1)=i-1;
               pcn(1)=i-1;
               val(1)=1;
               call MatSetValues(A,1,prn,1,pcn,val,ADD_VALUES,perr)
#endif
            end if
         end do

      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr)
   
    end subroutine build_A
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine build_Ai
      implicit none

      integer :: i,lrnl,lrnu,row,col,rbv,cbv,j,bnum,enum
      logical :: ilow,iup, eflag,check

      integer :: rw   
      integer :: ncls 
      integer :: cls(50) 
      
      !check to see if I'm computing a infrastructure source or not
      enum=eind(my_rank,1) 
      bnum=nbounds(e_nods(enum)) 
     
      !zero Ai
      call MatZeroEntries(Ai,perr)

      !if(.not. tank_flag) then
      do i=1,10*nelem
         row=rows(A_map(i))
         col=cols(A_map(i))
         rbv=nbounds(row)
         cbv=nbounds(col)
         
         eflag = .true.
         if(bnum<0 .and. (rbv==bnum .or. cbv==bnum)) eflag = .false.
         
         !lower triangle
         if(((rbv .eq. 2) .or. (cbv .eq. 2)) .and. .not. tank_flag) then
            !one or both nodes are on an outerboundary so set to zero for bc's
            val(1) = 0
#ifdef imimode
!         elseif( (rbv<0 .or. cbv < 0) .and. eflag ) then
!            !one or both nodes are one an ifrastructure boundary (and this is for an electrode pole)
!           val(1) = 0
#endif
         else
            
            val(1) = sigmai(S_map(i))*delA(i)
         end if
         
         prn(1) = row-1
         pcn(1) = col-1
         
         call MatSetValues(Ai,1,prn,1,pcn,val,ADD_VALUES,perr)          
         !upper triangle
         if(row .ne. col) then
            call MatSetValues(Ai,1,pcn,1,prn,val,ADD_VALUES,perr)
         end if
         
      end do
      
      
      !!fill in the diagonal for the zero potential bc's
      do i=1,nnodes
         !I own this node
         
         if((nbounds(i).eq.2) .and. .not. tank_flag) then
            !this node is on a boundary, set val to 1.0 for bc
            prn(1)=i-1;
            pcn(1)=i-1;
            val(1)=1;
            call MatSetValues(Ai,1,prn,1,pcn,val,ADD_VALUES,perr)
#ifdef imimode
!         elseif(nbounds(i) < 0 .and. nbounds(i) .ne. bnum) then
!            !one or both nodes are on an ifrastructure boundary (and this is for and electrode pole)
!            prn(1)=i-1;
!            pcn(1)=i-1;
!            val(1)=1;
!            call MatSetValues(Ai,1,prn,1,pcn,val,ADD_VALUES,perr)
#endif
         end if
      end do
      
      call MatAssemblyBegin(Ai,MAT_FINAL_ASSEMBLY,perr)
      call MatAssemblyEnd(Ai,MAT_FINAL_ASSEMBLY,perr)
  
    end subroutine build_Ai
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine receive_sigma
      implicit none
      
      call MPI_BCAST(nsig, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
      if(.not.allocated(sigma)) allocate(sigma(nsig))
      call MPI_BCAST(sigma, nsig,MPI_REAL,0,E4D_COMM,ierr)

    end subroutine receive_sigma
    !__________________________________________________________________

    !__________________________________________________________________
    subroutine receive_sigmai

      implicit none      
      integer :: i
      if(.not.allocated(sigmai)) allocate(sigmai(nsig))
      call MPI_BCAST(sigmai, nsig,MPI_REAL,0,E4D_COMM,ierr)
     
    end subroutine receive_sigmai
    !__________________________________________________________________

  
    !__________________________________________________________________
    subroutine send_dpred(flg)
      implicit none
      integer :: flg,i
      integer :: COMM

       if(im_fmm) then
          COMM = FMM_COMM
       else
          COMM = E4D_COMM
       end if
      
       !write(*,*) !get rid of this later 
      call assemble_data(flg)
      call MPI_SEND(nmy_drows,1,MPI_INTEGER,0,0,COMM, ierr)
      call MPI_SEND(my_drows,nmy_drows,MPI_INTEGER,0,0,COMM,ierr)
      call MPI_SEND(my_dvals,nmy_drows,MPI_REAL,0,0,COMM,ierr)
    
    end subroutine send_dpred
    !__________________________________________________________________


    !__________________________________________________________________
    subroutine weight_jaco
      implicit none
       
       integer :: ierr,i,j,row
       real, dimension(:), allocatable :: Wdt
       integer :: COMM
       
       if(im_fmm) then
          COMM = FMM_COMM
          nm = nm_fmm
       else
          COMM = E4D_COMM
       end if
       allocate(Wdt(nm))

       call MPI_BCAST(Wdt,nm,MPI_REAL,0,COMM,ierr)
 
       if(im_fmm) then
          do i=1,jind(my_rank_fmm,2) - jind(my_rank_fmm,1) + 1
             row = jind(my_rank_fmm,1) + i-1
             Jaco(i,:) = (Wdt(row))*Jaco(i,:)
          end do
       else
          do i=1,jind(my_rank,2) - jind(my_rank,1) + 1
             row = jind(my_rank,1) + i-1
             Jaco(i,:) = (Wdt(row))*Jaco(i,:)
          end do
       end if
       deallocate(Wdt)
       return
     
    end subroutine weight_jaco
    !__________________________________________________________________

  !__________________________________________________________________
    subroutine weight_jacoi
      implicit none
       
       integer :: ierr,i,j,row
       real, dimension(nm) :: Wdt
       integer :: COMM

      
       call MPI_BCAST(Wdt,nm,MPI_REAL,0,E4D_COMM,ierr)
       
       do i=1,jind(my_rank,2) - jind(my_rank,1) + 1
          row = jind(my_rank,1) + i-1
          Jacoi(i,:) = dble(Wdt(row))*Jacoi(i,:)
       end do
  
       return
     
    end subroutine weight_jacoi
    !__________________________________________________________________


    !__________________________________________________________________
     subroutine compute_pmatvec2

       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: m,i,j,ind1,ind2,tag,ierr
       real :: temp
       !real*8, dimension(:), allocatable :: x,sol
       integer :: COMM

       if(im_fmm) then
          COMM = FMM_COMM
          nm=nm_fmm
       else
          COMM = E4D_COMM
       end if

       tag = 0
       !call MPI_BCAST(m,1,MPI_INTEGER,0,E4D_COMM,ierr)
       
       !allocate(x(m))
       call MPI_BCAST(X2,nm,MPI_REAL,0,COMM,ierr)
      
       !allocate(sol(nelem))
       sol2 = 0
       if(im_fmm) then
          ind1 = jind(my_rank_fmm,1)
          ind2 = jind(my_rank_fmm,2)
       else
          ind1 = jind(my_rank,1)
          ind2 = jind(my_rank,2)
       end if
       
       !!do M'x
       do i=1,nelem
          sol2(i) = (dot_product((Jaco(:,i)),(X2(ind1:ind2))))
       end do
      
       call MPI_REDUCE(sol2,sol2,nelem,MPI_REAL,MPI_SUM,0,COMM,ierr)
      
       return

       !!Return the solution to master
       call MPI_SEND(sol2,nelem,MPI_REAL,0,tag,COMM, ierr)
      
     end subroutine compute_pmatvec2
     !__________________________________________________________________

     !__________________________________________________________________
     subroutine compute_pmatvec2_dbl

       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: m,i,j,ind1,ind2,tag,ierr
       real :: temp
       real*8, dimension(nm) :: X2d
       real*8, dimension(nelem) ::sol2d
       integer :: COMM
       
       if(im_fmm) then
          COMM = FMM_COMM
       else
          COMM = E4D_COMM
       end if
       tag = 0
       !call MPI_BCAST(m,1,MPI_INTEGER,0,E4D_COMM,ierr)

       call MPI_BCAST(X2d,nm,MPI_DOUBLE_PRECISION,0,COMM,ierr)
       
       !allocate(sol(nelem))
       sol2d = 0
       ind1 = jind(my_rank,1)
       ind2 = jind(my_rank,2)
    
       !!do M'x
       do i=1,nelem
          sol2d(i) = dot_product(dble(Jaco(:,i)),X2d(ind1:ind2))
       end do
      
       call MPI_REDUCE(sol2d,sol2d,nelem,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr)
    
       return

       !!Return the solution to master
       !!call MPI_SEND(sol2,nelem,MPI_REAL,0,tag,E4D_COMM, ierr)
      
     end subroutine compute_pmatvec2_dbl
     !__________________________________________________________________

    !__________________________________________________________________
     subroutine start_inv
       implicit none   
       if(im_fmm) nm=nm_fmm
       allocate(X1(nelem),sol1(nj_rows),X2(nm),sol2(nelem))
     end subroutine start_inv
    !__________________________________________________________________

     !__________________________________________________________________
     subroutine start_invi
       implicit none   
       allocate(X1i(nelem),sol1i(nj_rows),X2i(nm),sol2i(nelem))
     end subroutine start_invi
    !__________________________________________________________________
    
     !__________________________________________________________________
     subroutine send_Jcol
       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: n,i,j,tag,ierr,nrows,col
       real, dimension(1) :: dummy
       integer, dimension(1) :: recvcnts, displs
       !real*8, dimension(:), allocatable :: x,sol
       integer :: COMM
       recvcnts=0
       displs=0

       if(im_fmm) then
          COMM = FMM_COMM
       else
          COMM = E4D_COMM
       end if
       tag = 0
   
       call MPI_BCAST(col,1,MPI_REAL,0,COMM,ierr)
       
       do i=1,nj_rows
          sol1(i) = Jaco(i,col)
       end do
       
       call MPI_GATHERV(sol1,nj_rows,MPI_REAL,dummy,recvcnts,displs,MPI_REAL,0,COMM,ierr)
       return
     end subroutine send_Jcol
    !__________________________________________________________________


    !__________________________________________________________________
     subroutine compute_pmatvec1
       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: n,i,j,tag,ierr,nrows
       real, dimension(1) :: dummy
       integer, dimension(1) :: recvcnts, displs
       !real*8, dimension(:), allocatable :: x,sol
       integer :: COMM
       recvcnts=0
       displs=0

       if(im_fmm) then
          COMM = FMM_COMM
          nm=nm_fmm
       else
          COMM = E4D_COMM
       end if

       tag = 0
       
       !call MPI_BCAST(n,1,MPI_INTEGER,0,E4D_COMM,ierr)
       !allocate(x(n))
       call MPI_BCAST(X1,nelem,MPI_REAL,0,COMM,ierr)
       
       !nrows = jind(my_rank,2) - jind(my_rank,1) + 1
       !allocate(sol(nrows))
       
       !!compute Mx
       do i=1,nj_rows
          sol1(i) = (dot_product((Jaco(i,:)),(X1)))
       end do
       
       call MPI_GATHERV(sol1,nj_rows,MPI_REAL,dummy,recvcnts,displs,MPI_REAL,0,COMM,ierr)
       return
       
       !!Return the solution to master
       call MPI_SEND(sol1,nj_rows,MPI_REAL,0,tag,E4D_COMM, ierr)
      
     end subroutine compute_pmatvec1
    !__________________________________________________________________

 !__________________________________________________________________
     subroutine compute_pmatvec1_dbl
       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: n,i,j,tag,ierr,nrows
       real, dimension(1) :: dummy
       integer, dimension(1) :: idummy=0
       real*8, dimension(nelem) :: X1d
       integer :: COMM

       if(im_fmm) then
          COMM = FMM_COMM
       else
          COMM = E4D_COMM
       end if
       
       tag = 0
       
       !call MPI_BCAST(n,1,MPI_INTEGER,0,E4D_COMM,ierr)
       !allocate(x(n))
       
       call MPI_BCAST(X1d,nelem,MPI_DOUBLE_PRECISION,0,COMM,ierr)
       
       !nrows = jind(my_rank,2) - jind(my_rank,1) + 1
       !allocate(sol(nrows))
       
       !!compute Mx
       do i=1,nj_rows
          sol1(i) = real(dot_product(dble(Jaco(i,:)),X1d))
       end do
      
       call MPI_GATHERV(sol1,nj_rows,MPI_REAL,dummy,idummy,idummy,MPI_REAL,0,COMM,ierr)
      
       return
       
       !!Return the solution to master
       call MPI_SEND(sol1,nj_rows,MPI_REAL,0,tag,COMM, ierr)
      
     end subroutine compute_pmatvec1_dbl
    !__________________________________________________________________



     !__________________________________________________________________
     subroutine compute_pmatvec2i

       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: m,i,j,ind1,ind2,tag,ierr
       !real*8, dimension(:), allocatable :: x,sol
       
       tag = 0
       !call MPI_BCAST(m,1,MPI_INTEGER,0,E4D_COMM,ierr)

       !allocate(x(m))
       call MPI_BCAST(X2i,nm,MPI_DOUBLE_PRECISION,0,E4D_COMM,ierr)
 
       !allocate(sol(nelem))
       sol2i = 0
       ind1 = jind(my_rank,1)
       ind2 = jind(my_rank,2)
        
       !!do M'x
       do i=1,nelem
          sol2i(i) = dot_product(Jacoi(:,i),X2i(ind1:ind2))
       end do
       
       call MPI_REDUCE(sol2i,sol2i,nelem,MPI_DOUBLE_PRECISION,MPI_SUM,0,E4D_COMM,ierr)
       return

       !!Return the solution to master
       call MPI_SEND(sol2,nelem,MPI_REAL,0,tag,E4D_COMM, ierr)
      
     end subroutine compute_pmatvec2i
     !__________________________________________________________________

     !__________________________________________________________________
     subroutine compute_pmatvec1i
       !!receive the vector x from master and multiply it by my part of the 
       !!Jacobian transpose, then return the results
       implicit none
       integer :: n,i,j,tag,ierr,nrows
       real*8, dimension(1) :: dummy
       integer, dimension(1) :: idummy=0
       !real*8, dimension(:), allocatable :: x,sol
       
       tag = 0
       
       !call MPI_BCAST(n,1,MPI_INTEGER,0,E4D_COMM,ierr)
       !allocate(x(n))
       call MPI_BCAST(X1i,nelem,MPI_DOUBLE_PRECISION,0,E4D_COMM,ierr)
       
       !nrows = jind(my_rank,2) - jind(my_rank,1) + 1
       !allocate(sol(nrows))
       
       !!compute Mx
       do i=1,nj_rows
          sol1i(i) = dot_product(Jacoi(i,:),X1i)
       end do
       
       call MPI_GATHERV(sol1i,nj_rows,MPI_DOUBLE_PRECISION,dummy,idummy,idummy,MPI_DOUBLE_PRECISION,0,E4D_COMM,ierr)
       return
       
       !!Return the solution to master
       call MPI_SEND(sol1,nj_rows,MPI_REAL,0,tag,E4D_COMM, ierr)
      
     end subroutine compute_pmatvec1i
    !__________________________________________________________________

     !_________________________________________________________________
     subroutine get_J_on_off
       implicit none
       integer :: COMM

       if(im_fmm) then
          COMM = FMM_COMM
       else
          COMM = E4D_COMM
       end if

       if(.not. allocated(J_on_off)) allocate(J_on_off(nelem))
       J_on_off = .true.
      
       call MPI_BCAST(J_on_off,nelem,MPI_LOGICAL,0,COMM,ierr)
     
     end subroutine get_J_on_off
     !_________________________________________________________________


    !__________________________________________________________________
     subroutine end_inv
       implicit none
       deallocate(X1,sol1,X2,sol2)
     end subroutine end_inv
    !__________________________________________________________________

     !__________________________________________________________________
     subroutine end_invi
       implicit none
       deallocate(X1i,sol1i,X2i,sol2i)
     end subroutine end_invi
    !__________________________________________________________________


    !__________________________________________________________________
     subroutine write_myjaco
       implicit none
       integer :: i,j
       do i=1,nj_rows
          do j=1,nelem
             write(my_rank+100,*) jaco(i,j)
          end do
       end do
     end subroutine write_myjaco
    !__________________________________________________________________

    !__________________________________________________________________
     subroutine receive_rrseq
       implicit none
       if(allocated(rr_seq)) deallocate(rr_seq)
       allocate(rr_seq(nm))
       call MPI_BCAST(rr_seq,nm,MPI_INTEGER,0,E4D_COMM,ierr)
       call MPI_BCAST(pcheck,4,MPI_LOGICAL,0,E4D_COMM,ierr)

     end subroutine receive_rrseq
    !__________________________________________________________________

    !__________________________________________________________________
     subroutine send_pot
       implicit none
       integer, dimension(2) :: spack
       integer :: es

       call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)

       if(spack(1)==my_rank) then
          es=spack(2)-eind(my_rank,1) + 1
          call MPI_SEND(poles(:,es),nnodes,MPI_REAL,0,0,E4D_COMM,ierr)
       end if

     end subroutine send_pot
     !__________________________________________________________________
     
     !__________________________________________________________________
     subroutine send_poti
       implicit none
       integer, dimension(2) :: spack
       integer :: es
       
       call MPI_BCAST(spack,2,MPI_INTEGER,0,E4D_COMM,ierr)
       
       if(spack(1)==my_rank) then
          es=spack(2)-eind(my_rank,1) + 1
          call MPI_SEND(real(polesi(:,es)),nnodes,MPI_REAL,0,0,E4D_COMM,ierr)
       end if
       
     end subroutine send_poti
     !__________________________________________________________________
     
    !__________________________________________________________________
     subroutine send_sens
       implicit none
       integer :: njaco,i,tag,ierr
       real, dimension(:), allocatable :: sens
       integer :: COMM
       COMM = E4D_COMM
       njaco = jind(my_rank,2)-jind(my_rank,1)+1
       if(im_fmm) then
          COMM = FMM_COMM
          njaco = jind(my_rank_fmm,2)-jind(my_rank_fmm,1)+1
       end if
       
       allocate(sens(nelem))
       sens=0
       tag = 0
       do i=1,nelem
          sens(i)=real(dot_product(Jaco(1:njaco,i),Jaco(1:njaco,i)))
          !sens(i)=sum(Jaco(1:njaco,i))
       end do
       
       call MPI_SEND(sens,nelem,MPI_REAL,0,tag,COMM,ierr)
     
       deallocate(sens)
       
     end subroutine send_sens



    !__________________________________________________________________

    !__________________________________________________________________
     subroutine send_full_jaco
       implicit none
       integer :: njaco,i,tag,ierr
       njaco = jind(my_rank,2)-jind(my_rank,1)+1
       call MPI_SEND(Jaco,njaco*nelem,MPI_REAL,0,my_rank,E4D_COMM,ierr)
     end subroutine send_full_jaco
    !__________________________________________________________________


     !__________________________________________________________________
     subroutine slave_print_jaco_rows
       implicit none
       integer :: i,njrows,indx,j
       integer, dimension(:), allocatable :: ifres
       character(len=10) :: sindx
       
       call MPI_BCAST(njrows,1,MPI_INTEGER,0,E4D_COMM,ierr)
      
       allocate(ifres(njrows))
       call MPI_BCAST(ifres,njrows,MPI_INTEGER,0,E4D_COMM,ierr)
    
   
       do i=1,njrows
          if( (ifres(i) >= jind(my_rank,1)) .and. (ifres(i)) <= jind(my_rank,2)) then
             indx = ifres(i)-jind(my_rank,1)+1
             write(sindx,'(I6.6)') ifres(i)

             open(13,file='jacobian_row_'//trim(sindx)//'.txt',status='replace',action='write')
             write(13,*) nelem,1
             do j=1,nelem
                write(13,*) Jaco(indx,j)
             end do
             close(13)
          end if
       end do
       deallocate(ifres)
       
     end subroutine slave_print_jaco_rows
     !__________________________________________________________________
end module slave
