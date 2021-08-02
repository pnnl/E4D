module forward
 
  use vars
  use build_amap
#ifdef imimode
  use imi
#endif

  implicit none   
  logical :: i_flg = .false.

  
contains

 !____________________________________________________________________
  subroutine get_electrode_nodes
    implicit none
    real, dimension(nnodes) :: dist
    integer :: i
    integer, dimension(1) :: indb
    logical :: infra_node = .false.
    real :: mx,my,mz
   
    if(.not.allocated(e_nods)) then
       allocate (e_nods(tne))
    end if
  
    do i=1,tne
       dist=sqrt( (e_pos(i,1) - nodes(:,1))**2 + &
            (e_pos(i,2) - nodes(:,2))**2 + &
            (e_pos(i,3) - nodes(:,3))**2)
       
       indb = minloc(dist)

       e_nods(i) = indb(1)	   
    end	do

    !Compute the FE weighting functions to determine the current dist
    !at the 4 nodes of the element containing the point electrode
    !This will enable sources to be place off-node (tcj 5/1/2019)
    call compute_source_currents
    return

1001   format(" WARNING: ELECTRODE ",I5," HAS BEEN MOVED ",F10.4," TO THE NEAREST NODE")   
  end subroutine get_electrode_nodes
  !____________________________________________________________________
  
  !__________________________________________________________________
  subroutine build_ksp
    
#ifdef petsc_7
implicit none
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscksp.h90"
#else
#include "petsc/finclude/petscksp.h"
use petscksp
implicit none
#endif
    real*8 :: rtol = 1e-6 
    real*8 :: atol = 1e-12
    real*8 :: dtol = 500
    integer :: maxints = 10000
    !KSPSetTolerances(KSP ksp,double rtol,double atol,double dtol,int maxits);
    
    !Set up the KSP context
    call KSPCreate(PETSC_COMM_SELF,KS,perr)
    !call KSPSetOperators(KS,A,A,SAME_PRECONDITIONER,perr)
    call KSPSetOperators(KS,A,A,perr)
    
    call KSPGetPC(KS,P,perr)
    !call KSPSetType(KS,KSPGMRES,perr) !use default
    !call KSPGMRESSetRestart(KS,1000,perr);
    !call KSPGetTolerances(KS,rtol,atol,dtol,maxints,perr)
    call KSPSetTolerances(KS,rtol,atol,dtol,maxints,perr)
    !call KSPSetTolerances(KS,PETS_DEFAULT,PETS_DEFAULT,PETS_DEFAULT,PETS_DEFAULT,perr)
    call KSPSetFromOptions(KS,perr)
   
  end subroutine build_ksp
  !__________________________________________________________________

  !____________________________________________________________________
  subroutine forward_run
#ifdef petsc_7
implicit none
#else
#include "petsc/finclude/petscvec.h"
use petscvec
implicit none
#endif

    integer :: i,m,n,niter,j,enum,bnum
    !integer, dimension(1) :: eindx
    PetscInt, dimension(1) :: eindx
    real, dimension(2) :: pck
    PetscScalar :: val
    PetscInt :: p_int 
    real :: tstart, tend
    character*40 :: fl
    
    p_int = 1
    do i=1,my_ne
       call cpu_time(tstart)
       
       call VecGetArrayF90(psol,vloc,ierr)
       vloc(1:nnodes)=dble(poles(:,i))
       call VecRestoreArrayF90(psol,vloc,ierr)
       
       enum=eind(my_rank,1)+i-1

       val=0.0
       call VecSet(B,val,perr)

       if(i_flg) then
          !call Add_Jpp(i)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
       end if

       eindx(1)=e_nods(enum)
       val=1.0
 
#ifdef imimode
       if(nbounds(eindx(1))<0) then
          !build_A_inf build both A and B for non-point sources
          call build_A_inf(nbounds(eindx(1)))      
       else
         
          !call VecSetValues(B,p_int,eindx(1)-1,val,ADD_VALUES,perr)
          do j=1,4
             eindx(1) = source_nodes(enum,j)
             val = source_currents(enum,j)
             call VecSetValues(B,p_int,eindx(1)-1,val,ADD_VALUES,perr)
          end do
          if(tank_flag) then
             eindx(1)=e_nods(tne)
             val = -1.0
             call VecSetValues(B,p_int,eindx(1)-1,val,ADD_VALUES,perr)
          end if
       end if
#else
       p_int = 1
       !call VecSetValues(B,p_int,eindx(1)-1,val,ADD_VALUES,perr)
       do j=1,4
          eindx(1) = source_nodes(enum,j) 
          val = source_currents(enum,j)
          call VecSetValues(B,p_int,eindx(1)-1,val,ADD_VALUES,perr)
       end do
       if(tank_flag) then
          eindx(1)=e_nods(tne)
          val = -1.0
          call VecSetValues(B,p_int,eindx(1)-1,val,ADD_VALUES,perr)
       end if
#endif
       
       call VecAssemblyBegin(B,perr)
       call VecAssemblyEnd(B,perr) 
 
       call KSPSolve(KS,B,psol,perr)
       
       call VecGetArrayF90(psol,vloc,ierr)
       poles(:,i)= real(vloc(1:nnodes))
      
      
       call VecRestoreArrayF90(psol,vloc,ierr)

       call KSPGetIterationNumber(KS,niter,perr)
       call cpu_time(tend)
       pck(1)=tend-tstart
       pck(2)=real(niter)
       
       call MPI_SEND(pck,2,MPI_REAL,0,1,E4D_COMM,ierr)

    end do

    if(first_sol) then
       call KSPSetInitialGuessNonzero(KS,PETSC_TRUE,perr)
       first_sol=.false.
    end if
    
    call KSPDestroy(KS,perr)

#ifdef imimode    
    if(minval(nbounds)<0) then
       !if infrastracture inclusions exist, scale and 
       !sum solutions to honor physical flux conditions
       call infrastructure_adjust
    end if
    call check_fluxes
#endif

  end subroutine forward_run
  !____________________________________________________________________


  !____________________________________________________________________
  subroutine forward_runi
    implicit none
    integer :: i,m,n,niter,j,enum,row,col,rbv,cbv
    integer, dimension(1) :: eindx
    real, dimension(2) :: pck
    PetscScalar :: val
    real :: tstart, tend
   

    !!Get the ghost pole solution if this is a tank sim
  
    do i=1,my_ne
       call cpu_time(tstart) 
     
       !build B
       call VecGetArrayF90(psol,vloc,ierr)
       vloc(1:nnodes)=-dble(poles(:,i))
       call VecRestoreArrayF90(psol,vloc,ierr)
       call MatMult(Ai,psol,B,ierr)
       !call MatMult(A,psol,B,ierr)
       call VecAssemblyBegin(B,perr)
       call VecAssemblyEnd(B,perr)
       
#ifdef imimode    
    if(minval(nbounds)<0) then
       !we need to build A without the infinite conductivity
       !boundaries for the complex solution
       call build_A_common
    end if
#endif
      

       if(my_rank==1 .and. i==1) then
          call VecGetArrayF90(B,vloc,ierr)
          open(13,file='Jpp.pot',status='replace',action='write');
          write(13,*) nnodes, ' 1'
          do j=1,nnodes
             write(13,*) vloc(j)
          end do
          close(13)
          call VecRestoreArrayF90(B,vloc,ierr)
       end if;
  
       !set the solution vector psol
       call VecGetArrayF90(psol,vloc,ierr)
       vloc(1:nnodes)=dble(polesi(:,i))
       call VecRestoreArrayF90(psol,vloc,ierr)
         
       call KSPSolve(KS,B,psol,perr)
       
       call VecGetArrayF90(psol,vloc,ierr)
       polesi(:,i)= vloc(1:nnodes)
       !polesi(:,i)= vloc(1:nnodes)/poles(:,i)
       
       call VecRestoreArrayF90(psol,vloc,ierr)
       
       call KSPGetIterationNumber(KS,niter,perr)
       call cpu_time(tend)
       pck(1)=tend-tstart
       pck(2)=real(niter)
       call MPI_SEND(pck,2,MPI_REAL,0,1,E4D_COMM,ierr)
       
    end do
 
    call KSPDestroy(KS,perr)
    

  end subroutine forward_runi
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine Add_Jpp(i)
    implicit none
    integer :: i
    call VecGetArrayF90(psol,vloc,ierr)
    vloc(1:nnodes)=dble(polesi(:,i))
    call VecRestoreArrayF90(psol,vloc,ierr)
    call MatMult(Ai,psol,B,ierr)
  end subroutine Add_Jpp
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine write_Jpp(i)
    implicit none
    integer :: i,j
    PetscScalar :: val
    val=0.0
    call VecSet(B,val,perr)
    call VecGetArrayF90(psol,vloc,ierr)
    vloc(1:nnodes)=dble(polesi(:,i))
    call VecRestoreArrayF90(psol,vloc,ierr)
    call MatMult(Ai,psol,B,ierr)

    call VecAssemblyBegin(B,perr)
    call VecAssemblyEnd(B,perr) 

    call VecGetArrayF90(B,vloc,ierr)
    open(13,file='Jpp.pot',status='replace',action='write');
    write(13,*) nnodes,' 1 '
    do j=1,nnodes
       write(13,*) vloc(j)
    end do
    write(*,*) my_rank,sum((abs(vloc)))
    close(13)
    call VecRestoreArrayF90(B,vloc,ierr)
  end subroutine write_Jpp
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine set_iflg(set)
    implicit none
    logical :: set
    i_flg=set
  end subroutine set_iflg
  !____________________________________________________________________
  
  !__________________________________________________________________
    subroutine build_A_common
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
            
         end if
      end do
      
      call MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY,perr)
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,perr)
      
    end subroutine build_A_common
    !__________________________________________________________________


  
 !____________________________________________________________________
  subroutine mget_electrode_nodes(passed)
    implicit none
    logical :: passed
    real, dimension(nnodes) :: dist
    integer :: i,ierr
    integer, dimension(1) :: indb
    logical :: infra_node = .false.
   
    passed = .true.

    if(.not.allocated(e_nods)) then
       allocate (e_nods(tne))
    end if
  
    do i=1,tne
       dist=sqrt( (e_pos(i,1) - nodes(:,1))**2 + &
            (e_pos(i,2) - nodes(:,2))**2 + &
            (e_pos(i,3) - nodes(:,3))**2)
       
       indb = minloc(dist)
       
       if(my_rank==0) then
          if(dist(indb(1)) > 0.001) then
             !!Print a notice if an electrode does not lie directly on a node
             write(*,1001) i,dist(indb)
          end if
       end if
       e_nods(i) = indb(1)

       if(nbounds(e_nods(i)) .lt. 0 ) infra_node = .true.
       if(infra_node .and. (nbounds(e_nods(i)) .ne. int(e_pos(i,4)))) goto 100
       if(nbounds(e_nods(i)) .ge.0 .and. infra_node) goto 101
       
	   
    end	do

    return

 100    continue
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,"(A,I7.7)") "  For electrode number ",i
       write(51,*) " the survey file specifies and infrastructure boundary flag of ",int(e_pos(i,4))
       write(51,*) " but the electrode lies on node number ",e_nods(i)
       write(51,*) " which has a boundary flag of ",nbounds(e_nods(i))
       write(51,*) " The infrastructure boundary flag must match the node boundary flag."
       write(51,*) " Aborting ..."
       close(51)
       write(*,*) 
       write(*,"(A,I7.7)") "  For electrode number ",i
       write(*,*) " the survey file specifies and infrastructure boundary flag of ",int(e_pos(i,4))
       write(*,*) " but the electrode lies on node number ",e_nods(i)
       write(*,*) " which has a boundary flag of ",nbounds(e_nods(i))
       write(*,*) " The infrastructure boundary flag must match the node boundary flag."
       write(*,*) " Aborting ..."
       passed=.false.
       return
       
101    continue
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,"(A,I7.7,A)") "  Electrode ",i," is on a node with a positive boundary flag"
       write(51,*) " but one or more previous electrodes have been placed"
       write(51,*) " on a node with a negative boundary flag (i.e. an infrastructure node)."
       write(51,*) " Infrastructure electrodes must be listed after non-infrastructure "
       write(51,*) " in the survey file."
       write(51,*) " Aborting ..."
       close(51)
       write(51,*) 
       write(*,"(A,I7.7,A)") "  Electrode ",i," is on a node with a positive boundary flag"
       write(*,*) " but one or more previous electrodes have been placed"
       write(*,*) " on a node with a negative boundary flag (i.e. an infrastructure node)."
       write(*,*) " Infrastructure electrodes must be listed after non-infrastructure "
       write(*,*) " in the survey file."
       write(*,*) " Aborting ..."
       passed=.false.
       return
       
       
1001   format(" NOTICE: ELECTRODE ",I5," IS ",E10.4," DISTANCE UNITS FROM THE NEAREST NODE")   
  end subroutine mget_electrode_nodes
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine compute_source_currents
    !computes the source current distribution across the nodes of the
    !elements that contain the electrodes owned by this slave process
    implicit none
    integer :: i,j,k,near_node,enum,ecnt,fnd,ii
    real*8 :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,px,py,pz
    integer*8, dimension(4) :: INDX
    real*8 :: D0,D(4)
    real*8, dimension(4,4) :: A,A0
    real*8, dimension(4) :: X
    logical, dimension(nelem) :: check_ele
    
    allocate(source_nodes(tne,4),source_currents(tne,4))
    source_nodes = 0
    source_currents = 0.0
    
    do i=1,my_ne
       A0(:,4) = 1.0
       enum = eind(my_rank,1)+i-1
          
       !near_node is the nearest node to this electrode
       near_node = e_nods(enum)

       !find all of the elements that contain this node
       check_ele = .false.
       do j=1,4
          ecnt = COUNT(elements(:,j)==near_node)
          fnd = 0
          do k = 1,nelem
             if(elements(k,j) == near_node) then 
                check_ele(k) = .true.
                fnd=fnd+1
                if(fnd == ecnt) exit
             end if
          end do
       end do

       !loop over the relevant elements to determine which
       !element contains the source position
       ecnt = COUNT(check_ele)
       do ii = 1,nelem
          if(check_ele(ii)) then
             do j=1,4
                do k=1,3
                   A0(j,k)=nodes(elements(ii,j),k)
                end do
             end do
             A=A0
             D0 = det4(A)

             do j=1,4
                A=A0
                A(j,1:3)=[e_pos(enum,1),e_pos(enum,2),e_pos(enum,3)]
                D(j) = det4(A)

             end do
             if(D0<0) then
                if(D(1)<0.and.D(2)<0.and.D(3)<0.and.D(4)<0) then
                   exit
                end if
             elseif(D0>0) then
                if(D(1)>0.and.D(2)>0.and.D(3)>0.and.D(4)>0) then
                   exit
                end if
             end if
           
          end if
       end do
       if(ii>nelem) then
          !if we're here then the source lies on a node
          source_nodes(enum,:)=0
          source_nodes(enum,1)= near_node
          source_currents(enum,:)= 0.0
          source_currents(enum,1) = 1.0
       else
          source_nodes(enum,:) = elements(ii,:)
          A0(1,:)=1.0
          X(1) = 1.0
          do j=2,4
             do k=1,4
                A0(j,k) = nodes(elements(ii,k),j-1)
             end do
             X(j) = e_pos(enum,j-1)
          end do
          call MIGS(A0,4,A,INDX)
       
          do j=1,4
             source_currents(enum,j) = dot_product(A(j,:),X)
          end do
          
       end if
       !do j=1,4
       !   write(*,*) my_rank,i,source_nodes(i,j),source_currents(i,j)
       !end do
    end do

    !distribute the source currents and source nodes
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_nodes,tne*4,MPI_INTEGER,MPI_SUM,SCOMM,ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE,source_currents,tne*4,MPI_REAL,MPI_SUM,SCOMM,ierr)
 
  end subroutine compute_source_currents
  !____________________________________________________________________
end module forward
 
