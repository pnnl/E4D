module joint_invert

  use vars
  use fmm_vars
  use master
  use mod_con
  use input
  implicit none

 
contains
  
  !_____________________________________________________________________________________
  subroutine joint_pcgls
    implicit none
    
    integer :: i,j,ncon,k,indefinite,unstable,info,rbi,nbr
    integer :: offset,ind1,ind2,ii
    real :: et,t1,t2,cs,ce
    real :: norms0,xmax,normx,resNE,resNE_old,gammaNE
    real :: alpha,gamma,delta,norms,gamma1,gbeta,beta4cg
    real, dimension(:), allocatable :: q,r,b,u,v,r1
    real, dimension(nelem) :: s,p,xpar,pinv,h,s1
    real, dimension(nm) :: wdobs
    logical :: infoII,comm_flag=.true.
    
    call cpu_time(Cstart)
    cs=Cstart

    call send_command(20)
    if(im_fmm) then
       write(*,*) "FMM: STARTING INVERSION"
    else
       write(*,*) "E4D: STARTING INVERSION"
    end if

    if(im_fmm .and. iter.eq.1) then
       do i=1,nrz
          write(*,*) 'FMM smetric: ',smetric(i,:),', zwts: ',zwts(i,:)
       end do
    elseif(iter.eq.1) then
       do i=1,nrz
          write(*,*) 'E4D smetric: ',smetric(i,:),', zwts: ',zwts(i,:)
       end do
    end if

    !ncon=ccount !maxval(wrows)   
    xpar = 0
    

    if(im_fmm) then
       sigma = log(sigma)
       speed = sqrt(speed)
       nm=nm_fmm
       allocate(q(nm+ccount),r(nm+ccount),b(nm+ccount),r1(nm+ccount))
       b = 0
    else
       allocate(q(nm+ccount),r(nm+ccount),b(nm+ccount),r1(nm+ccount))
       b = 0
       if(invi) then
          sigmai = log(sigmai)
          sigma = log(sigma)
       else
          sigma = log(sigma)
          refsig=log(refsig)
          if(allocated(prefsig)) prefsig=log(prefsig)
       end if
       speed = sqrt(speed)
    end if
 
    if(invi) then
       do i=1,nm
          b(i) = (Wdi(i)*Wd_cull(i)*(dobsi(i)-dpredi(i)))
          wdobs(i) = Wdi(i)*Wd_cull(i)*dobsi(i)
       end do
    else
       do i=1,nm
          b(i) = (Wd(i)*Wd_cull(i)*(dobs(i)-dpred(i)))
          wdobs(i) = Wd(i)*Wd_cull(i)*dobs(i)
       end do
    end if

    if(wopt) then
       call cpu_time(ce)
       write(*,"(A,g10.4,A)") "  SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
       call send_data_noise
    end if

    !!Build the model part of b
    call cpu_time(ce)
    write(*,"(A,g10.4,A)") "  CONSTRUCTING CONSTRAINT RESIDUALS AT: ",ce-cs," seconds"
    do i=1,ccount
       
       if(invi .and. Wm(i) .eq. 0) then 
          goto 10
       elseif(Wm(i).eq.0) then
          goto 10
       end if

       rbi=rblock(i,3)
       select case(smetric(rbi,2))
       case(1) 
          if(im_fmm) then
               b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1))) 
          else
             if(invi) then
                !b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1))) 
                b(nm+i)=Wm(i)*(sigma(rblock(i,1))-sigmai(rblock(i,1)) - sigma(rblock(i,2)) + sigmai(rblock(i,2))) 
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1))) 
             end if
          end if
       case(2)
          if(im_fmm) then
             b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1))) 
          else
             if(invi) then
                !b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1)))
                b(nm+i)=Wm(i)*(sigma(rblock(i,1))-sigmai(rblock(i,1)) - sigma(rblock(i,2)) + sigmai(rblock(i,2))) 
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1))) 
             end if
          end if
       case(3)
          select case(smetric(rbi,3))
          case(0)
             if(im_fmm) then
                b(nm+i)=Wm(i)*(C_targ(rbi) - speed(rblock(i,1)))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*((C_targ(rbi))+ sigma(rblock(i,1)) - sigmai(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*(C_targ(rbi) - sigma(rblock(i,1)))
                end if
             end if
          case(1)
             if(im_fmm) then
                b(nm+i)=Wm(i)*(refsig(rblock(i,1)) - speed(rblock(i,1)))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(refsig(rblock(i,1))+ sigma(rblock(i,1)) - sigmai(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*(refsig(rblock(i,1)) - sigma(rblock(i,1)))
                end if
             end if
          case(2)
             if(im_fmm) then
                b(nm+i)=Wm(i)*(prefsig(rblock(i,1)) - speed(rblock(i,1)))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(prefsig(rblock(i,1)) + sigma(rblock(i,1)) - sigmai(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*(prefsig(rblock(i,1)) - sigma(rblock(i,1)))
                end if
             end if
          end select
          
          
       case(4)
          select case(smetric(rbi,3))
          case(0)
             if(im_fmm) then
                b(nm+i)=Wm(i)*(C_targ(rbi) - speed(rblock(i,1)))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(C_targ(rbi) + sigma(rblock(i,1)) - sigmai(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*(C_targ(rbi) - sigma(rblock(i,1)))
                end if
             end if
          case(1)
             if(im_fmm) then
                b(nm+i)=Wm(i)*(refsig(rblock(i,1)) - speed(rblock(i,1)))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(refsig(rblock(i,1)) - sigmai(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*(refsig(rblock(i,1)) - sigma(rblock(i,1)))
                end if
             end if
          case(2)
             if(im_fmm) then
                b(nm+i)=Wm(i)*(prefsig(rblock(i,1)) - speed(rblock(i,1)))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(prefsig(rblock(i,1)) - sigmai(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*(prefsig(rblock(i,1)) - sigma(rblock(i,1)))
                end if
             end if
          end select
          
       case(5)
          if(im_fmm) then
             b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1))) 
          else
             if(invi) then
                b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1))) 
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1))) 
             end if
          end if

       case(6)
          if(im_fmm) then
             b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1))) 
          else
             if(invi) then
                b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1)))   
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1))) 
             end if
          end if

       case(7)
          select case(smetric(rbi,3))
          case(0)
             if(im_fmm) then
                b(nm+i)=Wm(i)*( (C_targ(rbi) - speed(rblock(i,2))) - (C_targ(rbi) - speed(rblock(i,1))))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*( (C_targ(rbi) - sigmai(rblock(i,2))) - (C_targ(rbi) - sigmai(rblock(i,1))))
                else
                   b(nm+i)=Wm(i)*( (C_targ(rbi) - sigma(rblock(i,2))) - (C_targ(rbi) - sigma(rblock(i,1))))
                end if
             end if
          case(1)
             if(im_fmm) then
                b(nm+i)=Wm(i)*((speed(rblock(i,2))-refsig(rblock(i,2))) - &
                     (speed(rblock(i,1))-refsig(rblock(i,1))))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*((sigmai(rblock(i,2))-refsig(rblock(i,2))) - &
                        (sigmai(rblock(i,1))-refsig(rblock(i,1))))
                else
                   b(nm+i)=Wm(i)*((sigma(rblock(i,2))-refsig(rblock(i,2))) -&
                        (sigma(rblock(i,1))-refsig(rblock(i,1))))
                end if
             end if
          case(2)
             if(im_fmm) then
                b(nm+i)=Wm(i)*((speed(rblock(i,2))-prefsig(rblock(i,2))) - &
                     (speed(rblock(i,1))-prefsig(rblock(i,1))))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*((sigmai(rblock(i,2))-prefsig(rblock(i,2))) - &
                        (sigmai(rblock(i,1))-prefsig(rblock(i,1))))
                else
                   b(nm+i)=Wm(i)*((sigma(rblock(i,2))-prefsig(rblock(i,2))) - &
                        (sigma(rblock(i,1))-prefsig(rblock(i,1))))
                end if
             end if
          end select

      
       case(8)
          select case(smetric(rbi,3))
          case(0)
             if(im_fmm) then
                b(nm+i)=Wm(i)*( (C_targ(rbi) - speed(rblock(i,2))) &
                     - (C_targ(rbi) - speed(rblock(i,1))))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*( (C_targ(rbi) - sigmai(rblock(i,2))) &
                        - (C_targ(rbi) - sigmai(rblock(i,1))))
                else
                   b(nm+i)=Wm(i)*( (C_targ(rbi) - sigma(rblock(i,2))) &
                        - (C_targ(rbi) - sigma(rblock(i,1))))
                end if
             end if
          case(1)
             if(im_fmm) then
                b(nm+i)=Wm(i)*((speed(rblock(i,2))-refsig(rblock(i,2))) - &
                     (speed(rblock(i,1))-refsig(rblock(i,1))))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(sigma(rblock(i,1))-sigmai(rblock(i,1)) - & 
                        sigma(rblock(i,2)) + sigmai(rblock(i,2))+refsig(rblock(i,2))-refsig(rblock(i,1))) 
                else
                   b(nm+i)=Wm(i)*((sigma(rblock(i,2))-refsig(rblock(i,2))) - &
                        (sigma(rblock(i,1))-refsig(rblock(i,1))))
                end if
             end if
          case(2)
             if(im_fmm) then
                b(nm+i)=Wm(i)*((speed(rblock(i,2))-prefsig(rblock(i,2))) - &
                     (speed(rblock(i,1))-prefsig(rblock(i,1))))
             else
                if(invi) then
                   b(nm+i)=Wm(i)*(sigma(rblock(i,1))-sigmai(rblock(i,1)) - & 
                        sigma(rblock(i,2)) + sigmai(rblock(i,2))+prefsig(rblock(i,2))-prefsig(rblock(i,1)))
                else
                   b(nm+i)=Wm(i)*((sigma(rblock(i,2))-prefsig(rblock(i,2))) - &
                        (sigma(rblock(i,1))-prefsig(rblock(i,1))))
                end if
             end if
          end select

       case(9)
          if(im_fmm) then
             b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1)))
          else
             if(invi) then
                b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1))) 
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1)))
             end if
          end if
          
       case(10)
          if(im_fmm) then
             b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1))) 
          else
             if(invi) then
                b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1))) 
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1))) 
             end if
          end if
       case(11)
          if(im_fmm) then
             b(nm+i)=Wm(i)*(speed(rblock(i,2)) - speed(rblock(i,1))) 
          else
             if(invi) then
                b(nm+i)=Wm(i)*(sigmai(rblock(i,2)) - sigmai(rblock(i,1))) 
             else
                b(nm+i)=Wm(i)*(sigma(rblock(i,2)) - sigma(rblock(i,1))) 
             end if
          end if

       case(12)      
       	  if(im_fmm) then
             b(nm+i) = -Wm(i)*cg_wts(rblock(i,2),1)*speed(rblock(i,1))
             do j=1,4
                nbr = neighbors(rblock(i,1),j)
                if(nbr>0) then
                   b(nm+i) = b(nm+i) - Wm(i)*cg_wts(rblock(i,2),j+1)*speed(nbr)
                end if
             end do

          else
             b(nm+i) = -Wm(i)*cg_wts(rblock(i,2),1)*sigma(rblock(i,1))
             do j=1,4
                nbr = neighbors(rblock(i,1),j)
                if(nbr>0) then
                   b(nm+i) = b(nm+i) - Wm(i)*cg_wts(rblock(i,2),j+1)*sigma(nbr)
                end if
             end do

          end if
          
       case DEFAULT
       end select
10     continue
    end do
   
    b(nm+1:nm+ccount) = (sqrt(beta))*b(nm+1:nm+ccount)

    !At this point we need to compute B and v -> SHOULD COMPUTE ONLY WHEN beta4cg is non-zero
    call build_BM
    call build_v
    allocate(u(nnz_cgrad*(nphys-1)*nphys/10))
    allocate(v(nnz_cgrad*(nphys-1)*nphys/10))

    beta4cg = 0
    if (iter.eq.1) tmp_phi_0 = dot_product(b(1:nm),b(1:nm))    
    !if(dot_product(b(1:nm),b(1:nm))<0.06*dot_product(wdobs,wdobs)) then
    if(dot_product(b(1:nm),b(1:nm))<1.0*tmp_phi_0) then
       beta4cg = sqrt(betalist(1))*sqrt(betalist(2))
    end if
    if(cgmin_flag(2)) then
       call MPI_ALLREDUCE(MPI_IN_PLACE,beta4cg,1,MPI_REAL,MPI_PROD,M_COMM,ierr)
       beta4cg = sqrt(beta4cg)
    end if

    ! DBG
    !!beta4cg = sqrt(betalist(1))*sqrt(betalist(2))
    !beta4cg = 1.0 !beta4cg*2 !!for test
    
    if(im_fmm) write(*,*) 'FMM beta4cg:',iter,beta4cg
    if(.not. im_fmm) write(*,*) 'E4D beta4cg',iter,beta4cg

    ! if (im_fmm) then
    !    if (iter.eq.1) then
    !       open(359,file='Wm_fmm.txt',status='replace',action='write')
    !       write(359,*) ccount
    !       do i=1,ccount
    !          write(359,*) Wm(i)
    !       end do
    !       close(359)

    !       open(359,file='beta_fmm.txt',status='replace',action='write')
    !    else
    !       open(359,file='beta_fmm.txt',status='old',action='write',position='append')
    !    end if
    !    write(359,*) iter,beta,beta4cg,betalist(1)
    ! else
    !    if (iter.eq.1) then
    !       open(359,file='Wm_e4d.txt',status='replace',action='write')
    !       write(359,*) ccount
    !       do i=1,ccount
    !          write(359,*) Wm(i)
    !       end do
    !       close(359)

    !       open(359,file='beta_e4d.txt',status='replace',action='write')
    !    else
    !       open(359,file='beta_e4d.txt',status='old',action='write',position='append')
    !    end if
    !    write(359,*) iter,beta,beta4cg,betalist(1)
    ! end if
    ! close(359)

    do i=1,nphys-1
       offset=(block4cgrad(i,3)-i)*ncgrad*3
       ind1=1+(i-1)*ncgrad*3
       ind2=i*ncgrad*3
       !v(ind1+offset:ind2+offset) = -beta4cg*tnod(ind1:ind2)        ! YZ
       v(ind1+offset:ind2+offset) = -sqrt(beta4cg)*tnod(ind1:ind2)   ! PJ
    end do

    !Compute preconditioner
    call build_precondII(pinv,beta4cg)

    !Setup is done. The joint CGLS algorithm starts here
    r = b

    !!For test -> WHY NEED gammaNE in printing?
    r1(1:nm) = wdobs
    r1(nm+1:nm+ccount) = 0
    s1 = pmatvec2II(r1)
    gammaNE = dot_product(s1,pinv*s1)
    if (comm_flag .and. cgmin_flag(2)) then
       call MPI_ALLREDUCE(MPI_IN_PLACE,gammaNE,1,MPI_REAL,MPI_SUM,M_COMM,ierr)
    end if    
    
    ! Get vector s = pmatvec2II(r) + sqrt(beta4cg)*transpose(B)*v 
    s = pmatvec2II(r) + bmatvec2(dble(v),dble(beta4cg))
    h = pinv * s
    p = h
    
    gamma = dot_product(s,h)
    if (comm_flag .and. cgmin_flag(2)) then
       call MPI_ALLREDUCE(MPI_IN_PLACE,gamma,1,MPI_REAL,MPI_SUM,M_COMM,ierr)
    end if
    
    norms0 = sqrt(gamma)
    xmax = 0
    normx = 0
    k = 0
    info = 0
    infoII = .false.
    indefinite = 0
    unstable = 0
    resNE = 0
    call cpu_time(ce)
    write(*,"(A,g10.4,A)") "  STARTING JOINT INNER ITERATIONS AT: ",ce-cs," seconds"

    do i=1,max_initer

       !if( info .ne. 0) then
       if (infoII) exit

       !this is q1 or q2 depending on the physics in Yue's notation
       q = pmatvec1II(p)

       !compute u = sqrt(beta4cg)*B*p where B is the cross grad matrix (computed prior to inner iterations)
       u = bmatvec1(dble(p),dble(beta4cg))

       !add the u vectors for each type of physics to get u = sqrt(beta4cg)*B1*p1 + sqrt(beta4cg)*B2*p2
       if(comm_flag .and. cgmin_flag(2)) then
          call MPI_ALLREDUCE(MPI_IN_PLACE,u,nnz_cgrad*(nphys-1)*nphys/10,MPI_REAL,MPI_SUM,M_COMM,ierr)
       end if

       !compute my part of delta and then add the delta from the other physics using
       !MPI_ALLREDUCE again, then add u*trans(u)
       delta = dot_product(q,q) 
       if (comm_flag .and. cgmin_flag(2)) then          
          call MPI_ALLREDUCE(MPI_IN_PLACE,delta,1,MPI_REAL,MPI_SUM,M_COMM,ierr)
       end if
           
       delta = delta + dot_product(u,u)
       
       if (delta <= 0) indefinite = 1
       if (delta == 0) delta = epsilon(delta)

       alpha = gamma/delta
       
       xpar = xpar + alpha*p
       r = r - alpha*q
       !update v
       v = v - alpha*u

       ! s = pmatvec2II(r) + sqrt(beta4cg)*transpose(B)*v
       s = pmatvec2II(r) + bmatvec2(dble(v),dble(beta4cg))
       h = pinv * s
       
       norms = sqrt(dot_product(s,h))
       gamma1 = gamma
       gamma = norms*norms
       !add the component of gamma from the other physics
       if (comm_flag .and. cgmin_flag(2)) then          
          call MPI_ALLREDUCE(MPI_IN_PLACE,gamma,1,MPI_REAL,MPI_SUM,M_COMM,ierr)
       end if
  
       norms = sqrt(gamma)

       gbeta = gamma/gamma1
       p = h + gbeta*p

!print*,norm2(r),gamma       
       
       normx = sqrt(dot_product(xpar,xpar))
       if (xmax < normx) xmax = normx
       
       if( (norms <= norms0*initer_conv) .or. (normx*initer_conv >= 1)) then
          info = 1;
          infoII = .true.
       end if
       
       resNE_old = resNE
       resNE = norms/norms0
       call cpu_time(ce)
       if(im_fmm) then
          write(*,"(A,I3.3,A,f10.3,A,f10.4,f10.4)") "  FMM: ITERATION ",i," FINISHED AT: ",&
               (ce-cs)/60," minutes  ",resNE,norms/sqrt(gammaNE)
       else
           write(*,"(A,I3.3,A,f10.3,A,f10.4,f10.4)") "  E4D: ITERATION ",i," FINISHED AT: ",&
               (ce-cs)/60," minutes  ",resNE,norms/sqrt(gammaNE)
       end if

       ! Convergance criteria
       if( abs((resNE_old-resNe)/resNE_old) < delta_initer .and. i > min_initer) then
          !exit
          infoII = .true.
       end if
       if(comm_flag .and. cgmin_flag(2)) then
          call MPI_ALLREDUCE(MPI_IN_PLACE,infoII,1,MPI_LOGICAL,MPI_LAND,M_COMM,ierr)
       end if
     
    end do
    
    if(im_fmm) then
       sigma = exp(sigma)
       speed = speed**2
       sig_up = real(xpar)
    else
       if(invi) then
          sigmai = exp(sigmai)
          sigma = exp(sigma)
       else
          sigma = exp(sigma)
          refsig=exp(refsig)
          if(allocated(prefsig)) prefsig=exp(prefsig)
       end if
       !!map xpar to sigup
       sig_up = real(xpar)
       speed = speed**2
    end if
 100 continue
  
    call send_command(21)
    call cpu_time(Cend)
    etm=Cend-Cstart
 
    !deallocate(u,v)
    return
  end subroutine joint_pcgls
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine build_precondII(pinv,beta4cg)
    implicit none
    real, intent(in) :: beta4cg
    real, dimension(nelem), intent(out) :: pinv
    real*8, dimension(5)   :: diag
    real*8, dimension(5,5) :: atild
    real*8, dimension(5,5) :: atildi
    integer, dimension(5) :: indx
    integer :: i,j,k,iphys,incr_b,ind1,ind2,nbr

    call build_precond
    pinv=0
    incr_b=0
    do iphys=1,nphys-1
       do i=1,nelem
          if(use4cgrad(i) .and. beta4cg>0) then
             ! Compute Full Matrix -> No NEED
             !atild = 0
             !do j=1,5
             !   ind1=(j-1)*3+incr_b               
             !   do k=1,5
             !      ind2=(k-1)*3+incr_b
             !      atild(j,k)=dot_product(BM(ind1+1:ind1+3),BM(ind2+1:ind2+3))
             !   end do
             !end do

             !pinv(i) = pinv(i)+atild(1,1)
             !do j=1,4
             !   !pinv(neighbors(i,j)) = pinv(neighbors(i,j))+atildi(j+1,j+1)
             !   pinv(neighbors(i,j)) = pinv(i)+atild(j+1,j+1)
             !end do

             !! OR -> Uses only diagonal computations 
             !! Get only diagonal elements
             !!diag(1) = dot_product(BM(1+incr_b:3+incr_b),BM(1+incr_b:3+incr_b))
             pinv(i) = pinv(i)+dot_product(BM(1+incr_b:3+incr_b),BM(1+incr_b:3+incr_b))
             do j=1,4
                nbr = neighbors(i,j)
                !diag(j+1) = dot_product(BM(1+3*j+incr_b:3+3*j+incr_b),BM(1+3*j+incr_b:3+3*j+incr_b))
                pinv(nbr) = pinv(nbr)+dot_product(BM(1+3*j+incr_b:3+3*j+incr_b),BM(1+3*j+incr_b:3+3*j+incr_b))
             enddo
             
             !!atild(1,1)=dble(Precond(i))+dble(beta)+dble(beta4cg)*atild(1,1)
             !!do j=1,4
             !!   atild(j+1,j+1)=dble(Precond(neighbors(i,j)))+dble(beta)+&
             !!        dble(beta4cg)*atild(j+1,j+1)
             !!end do

             incr_b=incr_b+15

             !!call MIGS(atild,5,atildi,indx)
             !!pinv(i) = pinv(i)+atildi(1,1)
             !!pinv(i) = pinv(i)+1/(dble(Precond(i))+dble(beta)+atild(1,1))
             !pinv(i) = pinv(i)+atild(1,1)
             !do j=1,4
                !pinv(neighbors(i,j)) = pinv(neighbors(i,j))+atildi(j+1,j+1)
             !   pinv(neighbors(i,j)) = pinv(i)+atild(j+1,j+1)
             !end do
          !!else
          !!   pinv(i) = pinv(i)+1/(dble(Precond(i))+dble(beta))
          end if
       end do
    end do
    
    do i = 1,nelem
       !if (im_fmm.and.beta4cg>0) print*,i,pinv(i)       
       !!pinv(i) = 1/(dble(Precond(i))+dble(beta))
       pinv(i) = 1/(dble(Precond(i))+dble(beta)+dble(beta4cg)*pinv(i))
    end do

    !call MPI_BARRIER(M_COMM,ierr)    
    !stop
    
  end subroutine build_precondII
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  function bmatvec1(x,beta4cg)
    implicit none
    real*8 :: beta4cg
    real*8, dimension(nelem) :: x
    real*8, dimension(nnz_cgrad*(nphys-1)*nphys/10) :: bmatvec1
    integer :: i,j,nbr,iphys,incr,incr_b,offset

    ! DBG
    !if (im_fmm) then
    !   open(104,file='x2.txt',status='replace',action='write')
    !   open(105,file='B2x2.txt',status='replace',action='write')
    !else
    !   open(104,file='x1.txt',status='replace',action='write')
    !   open(105,file='B1x1.txt',status='replace',action='write')
    !endif
    !do i=1,nelem
    !   write(104,*) x(i)
    !enddo
    !close(104)
    
    bmatvec1=0
    incr=0
    incr_b=0
    do iphys=1,nphys-1
       offset=(block4cgrad(iphys,3)-iphys)*ncgrad*3
       incr=incr+offset
       do i=1,nelem
          if(use4cgrad(i)) then
             bmatvec1(1+incr)=BM(1+incr_b)*x(i)
             bmatvec1(2+incr)=BM(2+incr_b)*x(i)
             bmatvec1(3+incr)=BM(3+incr_b)*x(i)
!if (.not.im_fmm) print*,i,BM(1+incr_b),BM(2+incr_b),BM(3+incr_b)            
             do j=1,4
                nbr=neighbors(i,j)
                bmatvec1(1+incr)=bmatvec1(1+incr)+&
                     BM(1+j*3+incr_b)*x(nbr)
                bmatvec1(2+incr)=bmatvec1(2+incr)+&
                     BM(2+j*3+incr_b)*x(nbr)
                bmatvec1(3+incr)=bmatvec1(3+incr)+&
                     BM(3+j*3+incr_b)*x(nbr)
!if (.not.im_fmm) print*,nbr,BM(1+j*3+incr_b), BM(2+j*3+incr_b), BM(3+j*3+incr_b)                 
             end do
!call MPI_BARRIER(M_COMM,ierr)
!if (incr.eq.6) stop
             incr=incr+3
             incr_b=incr_b+15
          end if
       end do
       incr=incr-offset
    end do

!do i=1,3*ncgrad
!   write(105,*)i-1, bmatvec1(i)
!enddo
!close(105)
!call MPI_BARRIER(M_COMM,ierr)
!stop    
    
    !bmatvec1 = beta4cg*bmatvec1 !!is this correct
    bmatvec1 = sqrt(beta4cg)*bmatvec1 !!is this correct

  end function bmatvec1
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  function bmatvec2(x,beta4cg)
     implicit none
     real*8 :: beta4cg
     real*8, dimension(nnz_cgrad*(nphys-1)*nphys/10) :: x
     real*8, dimension(nelem) :: bmatvec2
     integer :: i,j,nbr,iphys,incr,incr_b,offset


    ! DBG
    !if (im_fmm) then
    !   open(104,file='y2.txt',status='replace',action='write')
    !   open(105,file='B2y2.txt',status='replace',action='write')
    !else
    !   open(104,file='y1.txt',status='replace',action='write')
    !   open(105,file='B1y1.txt',status='replace',action='write')
    !endif
    !do i=1,3*ncgrad
    !   write(104,*) x(i)
    !enddo
    !close(104)     
     
     bmatvec2=0
     incr=0
     incr_b=0
     do iphys=1,nphys-1
        offset=(block4cgrad(iphys,3)-iphys)*ncgrad*3
        incr=incr+offset
        do i=1,nelem
           if(use4cgrad(i)) then
              bmatvec2(i)=bmatvec2(i)+dot_product(BM(1+incr_b:3+incr_b),&
                   x(1+incr:3+incr))
              do j=1,4
                 nbr=neighbors(i,j)
                 bmatvec2(nbr)=bmatvec2(nbr)+dot_product(BM(1+j*3+incr_b:&
                      3+j*3+incr_b),x(1+incr:3+incr))
              end do

              incr=incr+3
              incr_b=incr_b+15
           end if
        end do
        incr=incr-offset
     end do

!DBG
!do i=1,nelem
!   write(105,*) i-1, bmatvec2(i)
!enddo
!close(105)
!call MPI_BARRIER(M_COMM,ierr)
!if (beta4cg.ne.0.0) stop 
     
     
     !bmatvec2 = beta4cg*bmatvec2 !!is this correct
     bmatvec2 = sqrt(beta4cg)*bmatvec2 !!is this correct
     
  end function bmatvec2
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  function pmatvec1II(x)
    
    implicit none
    real, dimension(nelem) :: x
    real, dimension(nm+ccount) :: pmatvec1II
    integer :: i,rbi,nbr,j,k
   
    pmatvec1II = 0
   
    
    call do_pmatvec1(nelem,nm,pmatvec1II(1:nm),x(1:nelem))
   
    do i=1,ccount
       if(Wm(i).eq.0) goto 10
       
       rbi=rblock(i,3)
       if(smetric(rbi,2) .lt. 12) then
          if(smetric(rbi,2).eq.3 .or. smetric(rbi,2).eq.4 ) then
             pmatvec1II(nm+i)=pmatvec1II(nm+i) + sqrt(beta)*Wm(i)*x(rblock(i,1))
          else
             pmatvec1II(nm+i)=pmatvec1II(nm+i) + sqrt(beta)*Wm(i)*(x(rblock(i,1))-x(rblock(i,2)))
          end if
       end if

!print*,"kl",my_rank,rbi,smetric(rbi,2) 
!       if(smetric(rbi,2).eq.12) then      
!          pmatvec1II(nm+i) = pmatvec1II(nm+i) + sqrt(beta)*Wm(i)*x(rblock(i,1))*cg_wts(rblock(i,2),1)
!          do j=1,4
!             nbr = neighbors(rblock(i,1),j)
!             if(nbr>0) then
!                pmatvec1II(nm+i) = pmatvec1II(nm+i) + &
!                     sqrt(beta)*Wm(i)*x(nbr)*cg_wts(rblock(i,2),j+1)
!             end if
!          end do
!          print*,"cg_wt:",cg_wts(rblock(i,2),j+1)
!       end if

10     continue
    end do
   
!stop
    
  end function pmatvec1II
  !_____________________________________________________________________________________ 


  !_____________________________________________________________________________________

  function pmatvec2II(r)
 
    !!computes A'x    
    implicit none 
    integer :: rbi,i,j,nbr
    real, dimension(nm+ccount) :: r
    real, dimension(nelem) :: pmatvec2II
  
   
    !!do the full matrix (Jacobian) part of the multiplication in parallel
    pmatvec2II = 0
    call do_pmatvec2(nelem,nm,pmatvec2II(1:nelem),r(1:nm))
  
    do i=1,ccount
      
       if(Wm(i).eq.0) goto 10
       rbi=rblock(i,3)

       if(smetric(rbi,2) .lt. 12) then
          if(smetric(rbi,2).eq.3 .or. smetric(rbi,2).eq.4 ) then
             pmatvec2II(rblock(i,1)) = pmatvec2II(rblock(i,1))+sqrt(beta)*r(nm+i)*Wm(i)
          else
             pmatvec2II(rblock(i,1)) = pmatvec2II(rblock(i,1))+sqrt(beta)*r(nm+i)*Wm(i)
             pmatvec2II(rblock(i,2)) = pmatvec2II(rblock(i,2))-sqrt(beta)*r(nm+i)*Wm(i)
          end if
       end if

!       if(smetric(rbi,2).eq.12) then
!          pmatvec2II(rblock(i,1)) = pmatvec2II(rblock(i,1)) + &
!               sqrt(beta)*Wm(i)*r(nm+i)*cg_wts(rblock(i,2),1)
!          do j=1,4
!             nbr = neighbors(rblock(i,1),j)
!             if(nbr .gt. 0) then
!                pmatvec2II(nbr) = pmatvec2II(nbr) + &
!                     sqrt(beta)*Wm(i)*r(nm+i)*cg_wts(rblock(i,2),j+1)
!             end if
!          end do
!       end if

10     continue
    end do
    
  end function pmatvec2II
  !_____________________________________________________________________________________


  !_____________________________________________________________________________________
  function pmatvec2(nm,nelem,ncon,ccount,wrows,wcols,Wm,r,beta,npar,par_map,nec,ex_cols,ex_vals)
 
    !!computes A'x    
    implicit none
    integer :: i,nm,nelem,ncon,ccount,npar,nec
    real :: beta
    integer, dimension(ccount) :: wrows,wcols
    integer, dimension(nm) :: par_map
    real, dimension(ccount) :: Wm
    real, dimension(nm+ncon) :: r
    integer, dimension(nec) :: ex_cols
    real, dimension(nec,3) :: ex_vals
    real, dimension(nelem) :: pmatvec
    real, dimension(npar) :: pmatvec2
    
    !!do the full matrix (Jacobian) part of the multiplication in parallel
    pmatvec = 0
    pmatvec2 = 0
 
    call do_pmatvec2(nelem,nm,pmatvec,r(1:nm))

    do i=1,nelem
       pmatvec2(par_map(i))= pmatvec2(par_map(i)) + pmatvec(i)
    end do

    !!do the sparse matrix part
    !!Wm=Wm*sqrt(beta)
    do i=1,ccount
       pmatvec2(wcols(i)) = pmatvec2(wcols(i)) + (sqrt(beta)* r(nm+wrows(i))*Wm(i)) 
    end do
    do i=1,nec
       pmatvec2(ex_cols(i)) = pmatvec2(ex_cols(i)) + (sqrt(beta)*r(nm+ncon+i)*ex_vals(i,1)*ex_vals(i,2))
    end do
  end function pmatvec2
  !_____________________________________________________________________________________

 !_____________________________________________________________________________________
  function pmatvec2_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,r,beta,npar,par_map,nec,ex_cols,ex_vals)
 
    !!computes A'x    
    implicit none
    integer :: i,nm,nelem,ncon,ccount,npar,nec
    real :: beta
    integer, dimension(ccount) :: wrows,wcols
    integer, dimension(nm) :: par_map
    real, dimension(ccount) :: Wm
    real, dimension(nm+ncon) :: r
    integer, dimension(nec) :: ex_cols
    real, dimension(nec,3) :: ex_vals
    real*8, dimension(nelem) :: pmatvec
    real, dimension(npar) :: pmatvec2_dbl
    
    !!do the full matrix (Jacobian) part of the multiplication in parallel
    pmatvec = 0
    pmatvec2_dbl = 0
 
    call do_pmatvec2_dbl(nelem,nm,pmatvec,r(1:nm))

    do i=1,nelem
       pmatvec2_dbl(par_map(i))= pmatvec2_dbl(par_map(i)) + real(pmatvec(i))
    end do

    !!do the sparse matrix part
    !!Wm=Wm*sqrt(beta)
    do i=1,ccount
       pmatvec2_dbl(wcols(i)) = pmatvec2_dbl(wcols(i)) + (sqrt(beta)* r(nm+wrows(i))*Wm(i)) 
    end do
    do i=1,nec
       pmatvec2_dbl(ex_cols(i)) = pmatvec2_dbl(ex_cols(i)) + (sqrt(beta)*r(nm+ncon+i)*ex_vals(i,1)*ex_vals(i,2))
    end do
   
  end function pmatvec2_dbl
  !_____________________________________________________________________________________



  !_____________________________________________________________________________________
  function pmatvec1(m,n,ncon,cc,Wm_row,Wm_col,Wm,x,reg,npar,par_map,nec,exc,exv)
    
    !!computes the matrix vector multiplication using
    !!the full Jacobian (J) and the sparse weighting matrix
    !!Wm
    implicit none
    integer :: m,n,cc,ncon,npar,nec,i
    integer, dimension(cc) :: Wm_col,Wm_row
    real :: reg
    integer, dimension(npar) :: par_map
    integer, dimension(nec) :: exc
    real, dimension(nec,3) :: exv
    real, dimension(cc) :: Wm
    real, dimension(npar) :: x
    real, dimension(n) :: xtmp
    real, dimension(m+ncon+nec) :: pmatvec1
    
    !!compute the full matrix mutliplication
    !!matvec1(1:m) = matmul(J,x)
    
    do i=1,n
       xtmp(i)=x(par_map(i))
    end do
    pmatvec1 = 0
    
    
    call do_pmatvec1(n,m,pmatvec1(1:m),xtmp)
    
    !!do the sparse part
    !!Wm=Wm*sqrt(reg)
    do i=1,cc
       pmatvec1(m+Wm_row(i)) =  pmatvec1(m+Wm_row(i)) + (sqrt(reg)*Wm(i)*x(Wm_col(i)))
    end do
    do i=1,nec
       pmatvec1(m+ncon+i) = pmatvec1(m+ncon+i) + (sqrt(reg)*exv(i,1)*exv(i,2)*x(exc(i)))
    end do
    
  end function pmatvec1
  !_____________________________________________________________________________________


 !_____________________________________________________________________________________
  function pmatvec1_dbl(m,n,ncon,cc,Wm_row,Wm_col,Wm,x,reg,npar,par_map,nec,exc,exv)
    
    !!computes the matrix vector multiplication using
    !!the full Jacobian (J) and the sparse weighting matrix
    !!Wm
    implicit none
    integer :: m,n,cc,ncon,npar,nec,i
    integer, dimension(cc) :: Wm_col,Wm_row
    real :: reg
    integer, dimension(npar) :: par_map
    integer, dimension(nec) :: exc
    real, dimension(nec,3) :: exv
    real, dimension(cc) :: Wm
    real, dimension(npar) :: x
    real, dimension(n) :: xtmp
    real, dimension(m+ncon+nec) :: pmatvec1_dbl
    
    !!compute the full matrix mutliplication
    !!matvec1(1:m) = matmul(J,x)
    
    do i=1,n
       xtmp(i)=x(par_map(i))
    end do
    pmatvec1_dbl = 0
    
    
    call do_pmatvec1_dbl(n,m,pmatvec1_dbl(1:m),dble(xtmp))
    
    !!do the sparse part
    !!Wm=Wm*sqrt(reg)
    do i=1,cc
       pmatvec1_dbl(m+Wm_row(i)) =  pmatvec1_dbl(m+Wm_row(i)) + (sqrt(reg)*Wm(i)*x(Wm_col(i)))
    end do
    do i=1,nec
       pmatvec1_dbl(m+ncon+i) = pmatvec1_dbl(m+ncon+i) + (sqrt(reg)*exv(i,1)*exv(i,2)*x(exc(i)))
    end do
    
  end function pmatvec1_dbl
  !_____________________________________________________________________________________







 !_____________________________________________________________________________________
  subroutine pcglsi
    implicit none
    
    integer :: i,ncon,k,indefinite,unstable,info
    real*8 :: et,t1,t2,cs,ce
    real*8 :: norms0,xmax,normx,resNE,resNE_old
    real*8 :: alpha,gamma,delta,norms,gamma1,gbeta
    real*8, dimension(:), allocatable :: q,r,b
    real*8, dimension(npar) :: s,p,xpar
    
    call cpu_time(Cstart)
    cs=Cstart
    call send_command(120)
    write(*,*) "STARTING COMPLEX COMPONENT INVERSION"

    ncon=maxval(wrows)  
    allocate(q(nm+ncon+nec),r(nm+ncon+nec),b(nm+ncon+nec))
    
    xpar = 0
    b = 0
    sigmai = log(sigmai)
    !refsig=log(refsig)
    !if(tl_flag .and. tl_ly) sig_not=log(sig_not)
 
    do i=1,nelem
       sigma_par(par_map(i))=sigmai(i)
       !rsigma_par(par_map(i))=refsig(i)
    end do
   
    !!Add the data weighting
    !if(tl_flag .and. tl_ly) then
       !do i=1,nm
          !b(i) = Wd(i)*Wd_cull(i)*((dobs(i)-dobs_not(i)) - (dpred(i)-dpred_not(i)))
       !end do
    !else
       do i=1,nm
          b(i) = dble(Wdi(i)*Wd_cull(i)*(dobsi(i)-dpredi(i)))
       end do
       
    !end if
    
    if(wopt) then
       call cpu_time(ce)
       write(*,*) "SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
       call send_data_noisei
       call cpu_time(ce)
       write(*,*) "DONE SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
    end if
    
   
    !!Add the model weighting 
    do i=1,ccount
       if(reg_opt(zones(wcols(i)))<3) then
          if(tl_flag .and. tl_ly) then
             b(nm+wrows(i)) = b(nm+wrows(i)) + dble(( sig_not(wcols(i)) - sigma_par(wcols(i)) )*Wm(i))
          else
             b(nm+wrows(i)) = b(nm+wrows(i)) - dble(sigma_par(wcols(i))*Wm(i))
          end if
       end if
       if(reg_opt(zones(wcols(i)))==3) then
          b(nm+wrows(i)) = b(nm+wrows(i)) + dble((rsigma_par(wcols(i))-sigma_par(wcols(i)))*Wm(i))
       end if
    end do

    !!add the external constraints
    do i=1,nec
       b(nm+ncon+i) = dble((ex_vals(i,2))*((ex_vals(i,3)-sigma_par(ex_cols(i)))))
    end do
    
    b(nm+1:nm+ncon+nec) = dble(sqrt(beta))*b(nm+1:nm+ncon+nec)

    r = b
    
    !!s = matmul(transpose(J),b)
    s = pmatvec2i(nm,nelem,ncon,ccount,wrows,wcols,Wm,r(1:nm+ncon),beta,npar,par_map,nec,ex_cols,ex_vals)  
 
    p = s
    gamma = dot_product(s,s)
  
    norms0 = sqrt(gamma)
    xmax = 0
    normx = 0
    k = 0
    info = 0
    indefinite = 0
    unstable = 0
    resNE = 0
    call cpu_time(ce)
    write(*,*) "STARTING INNER ITERATIONS AT: ",ce-cs," seconds"
   
    do i=1,max_initer
       !write(*,*) "ITERATION ",i
       if( info .ne. 0) then
          exit
       end if

       !call cpu_time(ce); write(*,*)"      executing pmatvec1 at:",(ce-cs)/60,' minutes'
       q = pmatvec1i(nm,nelem,ncon,ccount,wrows,wcols,Wm,p,beta,npar,par_map,nec,ex_cols,ex_vals) 
       !call cpu_time(ce); write(*,*)"      done with pmatvec1 at:",(ce-cs)/60,' minutes'

       delta = dot_product(q,q) 
       
       if(delta <= 0) then
          indefinite = 1
       end if
       
       if(delta == 0) then
          delta = epsilon(delta)
       end if
       
       alpha = gamma/delta
       
       xpar = xpar + alpha*p
       r = r - alpha*q

       !call cpu_time(ce); write(*,*)"      executing pmatvec2 at:",(ce-cs)/60,' minutes'
       s = pmatvec2i(nm,nelem,ncon,ccount,wrows,wcols,Wm,r,beta,npar,par_map,nec,ex_cols,ex_vals)
       !call cpu_time(ce); write(*,*)"      done with pmatvec2 at:",(ce-cs)/60,' minutes'

       norms = sqrt(dot_product(s,s))
       gamma1 = gamma
       gamma = norms*norms       
       gbeta = gamma/gamma1
       p = s + gbeta*p
      
       normx = sqrt(dot_product(xpar,xpar))
       if(xmax < normx) then
          xmax = normx
       end if
       
       if( (norms <= norms0*initer_conv) .or. (normx*initer_conv >= 1)) then
          info = 1;
       end if
       
       resNE_old = resNE
       resNE = norms/norms0
       call cpu_time(ce)
       write(*,*) "ITERATION ",i," FINISHED AT: ",(ce-cs)/60," minutes",resNE
       !write(*,'(I5,3F10.3)') i,resNE,abs((resNE_old-resNe)/resNE_old)

       if( abs((resNE_old-resNe)/resNE_old) < delta_initer .and. i > min_initer) then
          exit
       end if
    end do
    
    sigmai = exp(sigmai)
    !refsig=exp(refsig)
    if(tl_flag .and. tl_ly) sig_not=exp(sig_not)
 
    !!map xpar to sigup
    do i=1,nelem
       sig_up(i)=real(xpar(par_map(i)))
    end do
    write(*,*) maxval(sig_up),minval(sig_up)
    call send_command(121)
    call cpu_time(Cend)
    etm=Cend-Cstart
  
  end subroutine pcglsi
  !_____________________________________________________________________________________


  !_____________________________________________________________________________________
  function pmatvec2i(nm,nelem,ncon,ccount,wrows,wcols,Wm,r,beta,npar,par_map,nec,ex_cols,ex_vals)
    !!computes A'x
    
    implicit none
    integer :: i,nm,nelem,ncon,ccount,npar,nec
    real :: beta
    integer, dimension(ccount) :: wrows,wcols
    integer, dimension(nm) :: par_map
    real, dimension(ccount) :: Wm
    real*8, dimension(nm+ncon) :: r
    integer, dimension(nec) :: ex_cols
    real, dimension(nec,3) :: ex_vals
    real*8, dimension(nelem) :: pmatveci
    real*8, dimension(npar) :: pmatvec2i
    
    !!do the full matrix (Jacobian) part of the multiplication in parallel
    pmatveci = 0
    pmatvec2i = 0
    call do_pmatvec2i(nelem,nm,pmatveci,r(1:nm))
    do i=1,nelem
       pmatvec2i(par_map(i))= pmatvec2i(par_map(i)) + pmatveci(i)
    end do
    
    !!do the sparse matrix part
    !!Wm=Wm*sqrt(beta)
    do i=1,ccount
       pmatvec2i(wcols(i)) = pmatvec2i(wcols(i)) + dble(sqrt(beta))* r(nm+wrows(i))*dble(Wm(i)) 
    end do
    do i=1,nec
       pmatvec2i(ex_cols(i)) = pmatvec2i(ex_cols(i)) + dble(sqrt(beta))*r(nm+ncon+i)*dble(ex_vals(i,1)*ex_vals(i,2))
    end do
  end function pmatvec2i
  !_____________________________________________________________________________________


  !_____________________________________________________________________________________
  function pmatvec1i(m,n,ncon,cc,Wm_row,Wm_col,Wm,x,reg,npar,par_map,nec,exc,exv)
    
    !!computes the matrix vector multiplication using
    !!the full Jacobian (J) and the sparse weighting matrix
    !!Wm
    implicit none
    integer :: m,n,cc,ncon,npar,nec,i
    integer, dimension(cc) :: Wm_col,Wm_row
    real :: reg
    integer, dimension(npar) :: par_map
    integer, dimension(nec) :: exc
    real, dimension(nec,3) :: exv
    real, dimension(cc) :: Wm
    real*8, dimension(npar) :: x
    real*8, dimension(n) :: xtmp
    real*8, dimension(m+ncon+nec) :: pmatvec1i
    
    !!compute the full matrix mutliplication
    !!matvec1(1:m) = matmul(J,x)
    
    do i=1,n
       xtmp(i)=x(par_map(i))
    end do
    pmatvec1i = 0
    
    
    call do_pmatvec1i(n,m,pmatvec1i(1:m),xtmp)
    
    !!do the sparse part
    !!Wm=Wm*sqrt(reg)
    do i=1,cc
       pmatvec1i(m+Wm_row(i)) =  pmatvec1i(m+Wm_row(i)) + dble(sqrt(reg)*Wm(i)*x(Wm_col(i)))
    end do
    do i=1,nec
       pmatvec1i(m+ncon+i) = pmatvec1i(m+ncon+i) + dble(sqrt(reg)*exv(i,1)*exv(i,2)*x(exc(i)))
    end do
    
  end function pmatvec1i
  !_____________________________________________________________________________________



  
  
end module joint_invert
  
