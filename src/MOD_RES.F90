module mod_res
 
  use master
  use vars
  use forward
  use input
  use output
  use report
  use invert
  use mod_con
  use build_amap

  implicit none
  real,  dimension(:), allocatable :: RD,tk,qk,vk,dk,Rres,Rbase,Rtest,Wd_all,Wd_alli
  real, dimension(:), allocatable :: rg_rank,dist,Sbase,Stest,sg_rank,vlm
  integer :: k,ngrp,nbase,nmt,nR_inds,n_add,ntot
  integer, dimension(:),allocatable :: recvcounts,displs,mgrp,np_grp,R_inds,ig_rank
  integer, dimension(:,:), allocatable :: s_conf_all
  real, dimension(:,:), allocatable :: R_pts,mids
  logical, dimension(:),allocatable :: sused
  logical :: bflag
  logical :: setgrp0 = .false.
  logical :: cpsf_flag = .false.
  logical :: use_R = .false.
contains
  !____________________________________________________________________________________
  subroutine opt_srv
    implicit none
    integer :: t,i,j,itr,n_added
    integer, dimension(1) :: ind
    real :: trank
   
   
    call setup_res
    call get_volume
   
    allocate(Rres(nR_inds),Rbase(nR_inds),Rtest(nR_inds))
    allocate(rg_rank(ngrp),Sbase(nR_inds),Stest(nR_inds))
    allocate(sused(ngrp))
    allocate(sg_rank(ngrp))
    allocate(dist(nelem))
    sused=.false. 
    itr=0
    call write_survey(itr)
    open(23,file='opt_sequence.txt',status='replace',action='write')
    write(23,*) "ITERATION, GROUP_ADDED, QUALITY_FACTOR"
    close(23)    
    
10  continue
    itr=itr+1
    !!get the baseline survey res values
    bflag=.true.
    write(*,*)
    write(*,*) ">>>UPDATING RANKING FOR BASELINE SURVEY AT ITER:",itr," <<<"
    WD_cull = 1.0
    call comp_psf
    write(*,*)
    bflag=.false.
    
    rg_rank(t)=0
    sg_rank(t)=0

    do t=1,ngrp
       if(sused(t)) goto 15
       !!rebuild the survey with this group
       deallocate(s_conf,Wd,Wd_cull)
       nm=nbase+np_grp(t)
    
       allocate(s_conf(nm,4),Wd(nm),Wd_cull(nm))
       Wd_cull=1.0
       j=0
       do i=1,nmt
          if(mgrp(i)==0) then
             j=j+1
             s_conf(j,1:4)=s_conf_all(i,1:4)
             Wd(j)=Wd_all(i)
          end if
       end do
       do i=1,nmt
          if(mgrp(i)==t .or. sused(mgrp(i)) .and. (mgrp(i).ne.0)) then
             j=j+1
             s_conf(j,1:4)=s_conf_all(i,1:4)
             Wd(j)=Wd_all(i)
          end if
       end do
       
       !!send the new Jacobian assignments
       call send_info
       call send_command(26)
       call send_dists
       call build_rrseq
       call send_rrseq

       !instruct slaves to build the jacobian matrix
       call treport(0)
       !call send_command(11)
       call mjaco
       !call get_jtimes
       !call treport(3)
       write(*,*)
       write(*,*) ">>>UPDATING RANKING FOR MEASURMENT GROUP ",t," <<<"
       call comp_psf
     
      
       write(*,*) "  CONTROL POINT VALUES FOR GROUP",t," ITERATION ",itr
       do i=1,nR_inds
          write(*,*) "  control point: ",i," PSF: ",Sbase(i)/Stest(i)," R: ",Rtest(i)/Rbase(i)
       end do

       rg_rank(t)=sum(Rtest/(Rbase+1e-8))/nR_inds
       sg_rank(t)=sum(Sbase/(Stest+1e-9))/nR_inds
       write(*,*) "  RANKING FOR GROUP: ",t
       write(*,*) "  PSF: ",sg_rank(t)," R: ",rg_rank(t)
15     continue
    end do
    
    write(*,*)
    write(*,*) "  >>>RANKING SUMMARY<<<"
    do t=1,ngrp
       if(.not.sused(t)) then
          write(*,*) "  meas. group ",t," PSF: ",sg_rank(t)," R: ",rg_rank(t)
       end if
    end do

    n_added=0
    do t=1,ngrp
       if(sused(t)) goto 20
       if(use_R) then
          ind=maxloc(rg_rank)
       else
          ind=maxloc(sg_rank)
       end if
       sused(ind(1))=.true.
       n_added=n_added+np_grp(t)

       if(use_R) then
          trank = rg_rank(ind(1))
       else
          trank = sg_rank(ind(1))
       end if
       write(*,*)
       write(*,*) 'ADDING MEASUREMENT GROUP: ',ind(1),' WITH RANKING / #MEAS:',trank,np_grp(t)
       open(23,file='opt_sequence.txt',status='old',action='write',position='append')
       write(23,*) itr,ind(1),trank
       close(23)

       
       sg_rank(ind(1))=-sg_rank(ind(1))
       rg_rank(ind(1))=-rg_rank(ind(1))
       if(n_added>=n_add) goto 30
       
20     continue
    end do

30  continue
    nm=0
    do i=1,nmt
       if(mgrp(i)==0 .or. sused(mgrp(i))) nm=nm+1
    end do
    nbase=nm
    deallocate(s_conf,Wd,Wd_cull)
    allocate(s_conf(nm,4),Wd(nm),Wd_cull(nm))
    Wd_cull=1.0
    j=0
    do i=1,nmt
       if(mgrp(i)==0 .or. sused(mgrp(i))) then
          j=j+1
          s_conf(j,1:4)=s_conf_all(i,1:4)
          Wd(j)=Wd_all(i)
       end if
    end do
    if(j.ne.nm) write(*,*)"WAAAAAAAAIIIIIIIIIT"
    write(*,*)
    write(*,*) "___________________________________________________________________"
    write(*,*) "  STARTING NEXT OPTIMIZATION ITERATION WITH ",nm," MEASUREMENTS"

    call write_survey(itr)
    if(nm>=ntot) goto 40
    !!send the new Jacobian assignments
    call send_info
    call send_command(26)
    call send_dists
    call build_rrseq
    call send_rrseq
    
    !instruct slaves to build the jacobian matrix
    !call treport(0)
    !call send_command(11)
    call mjaco
    !call get_jtimes
    !call treport(3)
    goto 10

40  continue
    write(51,*) "SURVEY OPTIMIZATION COMPLETE"
    return
  end subroutine opt_srv
  !____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine gen_psf
    implicit none
    character*80 :: ofl
    integer :: i
    
    setgrp0 = .true.
    call setup_res
    call get_volume
    
    allocate(Rres(nR_inds),Rbase(nR_inds),Rtest(nR_inds))
    allocate(rg_rank(ngrp),Sbase(nR_inds),Stest(nR_inds))
    allocate(sused(ngrp))
    allocate(sg_rank(ngrp))
    allocate(dist(nelem))
    
    !!get the baseline survey res values
    bflag=.true.
    write(*,*)
    write(*,*) ">>> COMPUTING POINT SPREAD FUNCTIONS <<<"
    WD_cull = 1.0
    cpsf_flag = .true.
    call comp_psf
    cpsf_flag = .false. 
    setgrp0 = .false.

    write(ofl,"(A,A4)") trim(efile),".spd"
    write(*,*) ">>> WRITING POINT SPREADS TO: ",trim(ofl)," <<<"
    open(23,file=trim(ofl),action='write',status='replace')
    write(23,*) "ELEMENT NUMBER,X_POS,Y_POS,Z_POS,SPREAD,DIAG_R"
    do i=1,nR_inds
       write(23,*) R_inds(i),mids(i,1),mids(i,2),mids(i,3),Sbase(i),Rbase(i)
    end do
    close(23)
    
     
    return

  end subroutine gen_psf
  !_____________________________________________________________________________________
  

  !____________________________________________________________________________________
  subroutine compute_res
    implicit none
    integer :: i

    setgrp0 = .true. 
    call setup_res
    call get_volume
    call compute_mids
    allocate(Rbase(nR_inds),Sbase(nR_inds))
    allocate(dist(nelem))

    write(*,*) ">>> COMPUTING TRUE RESOLUTION AT CONTROL POINTS <<<"
    bflag = .true.
    Wd_cull = 1.0
    call comp_psf

    write(*,*) ">>> ESTIMATING RESOLUTION MATRIX DIAGONAL <<<"
    call rpcgls
   
    open(23,file='R_pred.txt',action='write',status='replace')
    write(23,*) nelem
    do i=1,nelem
       write(23,*) RD(i)
    end do
    close(23)
    
  end subroutine compute_res
  !____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine setup_res
    implicit none 
    integer :: i
    
    !read the survey info
    call get_survey

    !Send things the slaves need
    call send_info
    
    !send the conductivity vector to the slaves
    call send_sigma
    if(i_flag) call send_sigmai

    !setup the forward coupling matrix mappings
    call setup_forward
    call get_cpts
    call get_time
    call treport(1)
    if(i_flag) call send_command(102)
    call get_time
    call treport(1)
   
   
    !build the coupling matrix (command 3 for slaves)
    call send_command(3)
    if(i_flag) call send_command(103)

    !set the iteration number
    iter = 0
    
    !execute a forward run
    call send_command(5)
    call run_forward
    if(i_flag) then
       call send_command(5)
       call run_forwardi
    end if
   
    !Read the inverse options and build the reg matrix
    !while slaves are solving
    call get_inv_opts
    call nreport(2)
    call build_Wm
    call nreport(3)
    call build_rrseq
    call treport(12)
    call build_sgroupm
 
    !get the forward run times from slaves and report
    call treport(0)
    call get_abtimes
    call treport(7)
    call get_frtimes
    call treport(6)

    !send the Jacobian build sequence to the slave
    call send_rrseq
 
    !instruct slaves to build the jacobian matrix
    call treport(0)
    !call send_command(11)
    call mjaco
    call get_jtimes
    call treport(3)

    !allocate the update vector if not done
    !call alloc_sigup
   
    !!info used in vk_matmul
    allocate(recvcounts(n_rank),displs(n_rank))

   
  end subroutine setup_res
  !_____________________________________________________________________________________


  !_____________________________________________________________________________________
  subroutine rpcgls
    implicit none
    
    integer :: i,ncon,j,k,indefinite,unstable,info,ind1,ind2,nrow
    real :: et,t1,t2,cs,ce
    real :: norms0,xmax,normx,resNE,resNE_old
    real :: alpha,gamma,delta,norms,gamma1,gbeta
    real, dimension(:), allocatable :: q,r,b
    real, dimension(npar) :: s,p,xpar
    
    call cpu_time(Cstart)
    cs=Cstart
    call send_command(20)
    write(*,*) "STARTING RESOLUTION MATRIX BUILD"
    if(allocated(dk)) deallocate(dk)
    allocate(dk(nm))
    allocate(vk(nelem),tk(nelem),qk(nelem),RD(nelem))
    tk=0
    qk=0
    dk=0
    RD=0

    ncon=maxval(wrows)  
    allocate(q(nm+ncon+nec),r(nm+ncon+nec),b(nm+ncon+nec))
    
    xpar = 0
    b = 0
    sigma = log(sigma)
    refsig=log(refsig)
    if(tl_flag .and. tl_ly) sig_not=log(sig_not)
 
    do i=1,nelem
       sigma_par(par_map(i))=sigma(i)
       rsigma_par(par_map(i))=refsig(i)
    end do
    
    !data noise already sent by  comp_psf
    !if(wopt) then
    !   call cpu_time(ce)
    !   write(*,*) "SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
    !   call send_data_noise
    !   call cpu_time(ce)
    !   write(*,*) "DONE SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
    !end if
    
    
    !!Add the model weighting 
    do i=1,ccount
       if(reg_opt(zones(wcols(i)))<3) then
          if(tl_flag .and. tl_ly) then
             b(nm+wrows(i)) = b(nm+wrows(i)) + ( sig_not(wcols(i)) - sigma_par(wcols(i)) )*Wm(i)
          else
             b(nm+wrows(i)) = b(nm+wrows(i)) - sigma_par(wcols(i))*Wm(i)
          end if
       end if
       if(reg_opt(zones(wcols(i)))==3) then
          b(nm+wrows(i)) = b(nm+wrows(i)) + ((rsigma_par(wcols(i))-sigma_par(wcols(i)))*Wm(i))
       end if
    end do

    !!add the external constraints
    do i=1,nec
       b(nm+ncon+i) = (ex_vals(i,2))*((ex_vals(i,3)-sigma_par(ex_cols(i))))
    end do
    
    b(nm+1:nm+ncon+nec) = (sqrt(beta))*b(nm+1:nm+ncon+nec)

    recvcounts=0
    displs=0
    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       displs(i+1)=ind1-1
       recvcounts(i+1)=nrow
    end do

    call init_random_seed()
    do k=1,50000

       xpar=0
       do j=1,nelem
          vk(j)=random_normal()
       end do
       call vk_matmul
       b(1:nm)=dk

       r = b

       !!s = matmul(transpose(J),b)
       s = pmatvec2_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,r(1:nm+ncon),beta,npar,par_map,nec,ex_cols,ex_vals)  
       
       p = s
       gamma = dot_product(s,s)
       
       norms0 = sqrt(gamma)
       xmax = 0
       normx = 0
       
       info = 0
       indefinite = 0
       unstable = 0
       resNE = 0
       call cpu_time(ce)
       !write(*,*) "STARTING ITERATION ",k, " AT: ",ce-cs," seconds"
       do i=1,max_initer
          !write(*,*) "ITERATION ",i
          if( info .ne. 0) then
             exit
          end if
          
          !call cpu_time(ce); write(*,*)"      executing pmatvec1 at:",(ce-cs)/60,' minutes'
          q = pmatvec1_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,p,beta,npar,par_map,nec,ex_cols,ex_vals) 
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
          s = pmatvec2_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,r,beta,npar,par_map,nec,ex_cols,ex_vals)
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
          !call cpu_time(ce)
          !write(*,*) "ITERATION ",i," FINISHED AT: ",(ce-cs)/60," minutes",resNE
          !write(*,'(I5,3F10.3)') i,resNE,abs((resNE_old-resNe)/resNE_old)
          
          if( abs((resNE_old-resNe)/resNE_old) < delta_initer .and. i > min_initer) then
             exit
          end if
       end do
       
       do j=1,nelem
          tk(j)=tk(j)+real(xpar(par_map(j)))*vk(j)
          qk(j)=qk(j)+vk(j)*vk(j)
          RD(j)=tk(j)/qk(j)
          !write(*,*) real(xpar(par_map(j))),vk(j),tk(j),qk(j),RD(j)
       end do
       write(*,*) "____________________________________________________________"
       do i=1,nR_inds
          write(*,*) k,i,RD(R_inds(i)),Rbase(i),(Rbase(i)-RD(R_inds(i)))/Rbase(i)
       end do    

    end do
    !do j=1,nelem
    !   write(10,*) tk(j),qk(j),RD(j)
    !end do

    sigma = exp(sigma)
    refsig=exp(refsig)
    if(tl_flag .and. tl_ly) sig_not=exp(sig_not)
 
 
    call send_command(21)
    call cpu_time(Cend)
    etm=Cend-Cstart
  
  end subroutine rpcgls
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine comp_psf
    implicit none
    
    integer :: i,ncon,j,k,indefinite,unstable,info,ind1,ind2,nrow
    real :: et,t1,t2,cs,ce,C
    real :: norms0,xmax,normx,resNE,resNE_old
    real :: alpha,gamma,delta,norms,gamma1,gbeta
    real, dimension(:), allocatable :: q,r,b
    real*8, dimension(:), allocatable :: sup
    real, dimension(npar) :: s,p,xpar
    character*20 :: ofl
    
    if(.not.allocated(Rres)) allocate(Rres(nR_inds))
    if(allocated(dk)) deallocate(dk)
    allocate(dk(nm))
    allocate(sup(nelem))
    call cpu_time(Cstart)
    cs=Cstart
    call send_command(20)
   

    ncon=maxval(wrows)  
    allocate(q(nm+ncon+nec),r(nm+ncon+nec),b(nm+ncon+nec))
    
    xpar = 0
    b = 0
    sigma = log(sigma)
    refsig=log(refsig)
    if(tl_flag .and. tl_ly) sig_not=log(sig_not)
 
    do i=1,nelem
       sigma_par(par_map(i))=sigma(i)
       rsigma_par(par_map(i))=refsig(i)
    end do
    
    !if(wopt) then
       !call cpu_time(ce)
       !write(*,*) "SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
       call send_data_noise
       !call cpu_time(ce)
       !write(*,*) "DONE SENDING DATA NOISE TO SLAVES AT: ",ce-cs," seconds"
    !end if
    
    
    !!Add the model weighting 
    do i=1,ccount
       if(reg_opt(zones(wcols(i)))<3) then
          if(tl_flag .and. tl_ly) then
             b(nm+wrows(i)) = b(nm+wrows(i)) + ( sig_not(wcols(i)) - sigma_par(wcols(i)) )*Wm(i)
          else
             b(nm+wrows(i)) = b(nm+wrows(i)) - sigma_par(wcols(i))*Wm(i)
          end if
       end if
       if(reg_opt(zones(wcols(i)))==3) then
          b(nm+wrows(i)) = b(nm+wrows(i)) + ((rsigma_par(wcols(i))-sigma_par(wcols(i)))*Wm(i))
       end if
    end do

    !!add the external constraints
    do i=1,nec
       b(nm+ncon+i) = (ex_vals(i,2))*((ex_vals(i,3)-sigma_par(ex_cols(i))))
    end do
    
    b(nm+1:nm+ncon+nec) = (sqrt(beta))*b(nm+1:nm+ncon+nec)

   

    !!info used in vk_matmul

    !allocate(recvcounts(n_rank),displs(n_rank))
    recvcounts=0
    displs=0
    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       displs(i+1)=ind1-1
       recvcounts(i+1)=nrow
    end do

    !open(21,file='res.txt',action='write',status='replace')

    do k=1,nR_inds

       xpar=0
       call get_Jcol(R_inds(k))
       b(1:nm)=dk
   

       r = b

       !!s = matmul(transpose(J),b)
       s = pmatvec2_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,r(1:nm+ncon),beta,npar,par_map,nec,ex_cols,ex_vals)  
       
       p = s
       gamma = dot_product(s,s)
       
       norms0 = sqrt(gamma)
       xmax = 0
       normx = 0
       
       info = 0
       indefinite = 0
       unstable = 0
       resNE = 0
       call cpu_time(ce)
       if(bflag) write(*,*) "  Computing baseline PSF at optimization control point: ",k
       if(.not.bflag) write(*,*) "  Computing PSF at optimization control point: ",k
       do i=1,max_initer
          !write(*,*) "ITERATION ",i
          if( info .ne. 0) then
             exit
          end if
          
          !call cpu_time(ce); write(*,*)"      executing pmatvec1 at:",(ce-cs)/60,' minutes'
          q = pmatvec1_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,p,beta,npar,par_map,nec,ex_cols,ex_vals) 
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
          s = pmatvec2_dbl(nm,nelem,ncon,ccount,wrows,wcols,Wm,r,beta,npar,par_map,nec,ex_cols,ex_vals)
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
          !write(*,*) "ITERATION ",i," FINISHED AT: ",(ce-cs)/60," minutes",resNE
          !write(*,'(I5,3F10.3)') i,resNE,abs((resNE_old-resNe)/resNE_old)
          
          if( abs((resNE_old-resNe)/resNE_old) < delta_initer .and. i > min_initer) then
             exit
          end if
       end do
       

       do i=1,nelem
          sup(i)=xpar(par_map(i))
       end do
      
       dist=((mids(:,1)-mids(R_inds(k),1))**2 + &
            (mids(:,2)-mids(R_inds(k),2))**2 + &
            (mids(:,3)-mids(R_inds(k),3))**2)

       if(bflag) then
          Sbase(k)=0
     
          Sbase(k)= real((dot_product(dble(sqrt(dist))*abs(sup),dble(vlm)) + &
               (abs(sup(R_inds(k)))-1)*dble(vlm(R_inds(k))))/(1e-8+dot_product(abs(sup),dble(vlm))))
                         
          Rbase(k)=real(abs(sup(R_inds(k))))
  
       else
          Stest(k)=0         
          Stest(k)= real((dot_product(dble(sqrt(dist))*abs(sup),dble(vlm)) + &
               (abs(sup(R_inds(k)))-1)*dble(vlm(R_inds(k))))/(1e-8+dot_product(abs(sup),dble(vlm))))

          Rtest(k)= real(abs(sup(R_inds(k))))
         
       end if 
       
    end do
    

    sigma = exp(sigma)
    refsig=exp(refsig)
    if(tl_flag .and. tl_ly) sig_not=exp(sig_not)
 
 
    call send_command(21)
    call cpu_time(Cend)
    etm=Cend-Cstart
   
  end subroutine comp_psf
  !_____________________________________________________________________________________

 
  !_____________________________________________________________________________________
  subroutine write_psf(pt)
    implicit none
    integer :: pt,i
    character*20 :: fname
    real :: C

    write(*,*) 'WRITING PSF: ',pt
    write(fname,'(A4,I0,A4)') "psf_",k,".txt"
    write(*,*) fname
    open(21,file=trim(fname),action='write',status='replace')
    write(21,*) nelem
    C=nelem*(sum(vlm*sqrt(sig_up**2)))
    !C=nelem*(sum(sig_up**2))

    do i=1,nelem
       write(21,*)  vlm(i)*sqrt(dist(i))*sqrt(sig_up(i)**2)/C
       !write(21,*)  sqrt(dist(i))*sqrt(sig_up(i)**2)/(C*vlm(i))
    end do
    close(21)
    return
  end subroutine write_psf
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine vk_matmul
    implicit none
    integer :: ierr
    real :: dummy

    call send_command(14)
    call MPI_BCAST(vk,nelem,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_GATHERV(dummy,0,MPI_REAL,dk,recvcounts,displs,MPI_REAL,0,E4D_COMM,ierr)
    return

  end subroutine vk_matmul
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine get_Jcol(k)

    implicit none
    integer :: k
    real :: dummy

    call send_command(25)
    call MPI_BCAST(k,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_GATHERV(dummy,0,MPI_REAL,dk,recvcounts,displs,MPI_REAL,0,E4D_COMM,ierr)
    return

  end subroutine get_Jcol
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine get_survey
    implicit none
    integer :: i,j,junk,ropt
    real, dimension(4) :: etmp
   

    open(51,file='e4d.log',status='old',action='write',position='append')
    open(10,file=efile,status='old',action='read')   
    read(10,*,IOSTAT=ios) ne
    if(ios .ne. 0) goto 1007
    write(51,*) "NUMBER OF ELECTRODE LOCATIONS = ",ne
    
    
    allocate(e_pos(ne,4))
    do i=1,ne
       read(10,*,IOSTAT=ios) junk,etmp
       if(ios .ne. 0) goto 1008
       if(junk>ne) then
          write(51,*) 'ELECTRODE ',junk,' IS GREATER THAN THE NUMBER OF ELLECTRODES ',ne
          close(51)
          call crash_exit
       end if
       e_pos(junk,1:4)=etmp
    end do
    close(51)
    
    !!translate electrodes 
    open(51,file='e4d.log',status='old',action='write',position='append')
    if(mode > 1) then
       open(21,file='trans.txt',status='old')
       read(21,*,IOSTAT=ios) xorig,yorig,zorig
       if(ios .ne. 0) goto 1009
       close(21)
       write(51,*) 'TRANSLATING ELECTRODES'
       e_pos(:,1) = e_pos(:,1)-xorig
       e_pos(:,2) = e_pos(:,2)-yorig
       e_pos(:,3) = e_pos(:,3)-zorig
    end if
    close(51)
    
    !!Read in the survey
    open(51,file='e4d.log',status='old',action='write',position='append')
    read(10,*,IOSTAT=ios) nm,n_add,ntot,ropt
    if(ios .ne. 0) goto 1010
    if(ntot>=nm) goto 1015
    if(ropt == 2) then
       use_R = .true.
    end if
    write(51,*) "NUMBER OF LISTED MEASUREMENTS = ",nm
    write(51,*) "NUMBER OF MEASUREMENTS TO ADD PER ITERATION = ",n_add
    write(51,*) "TARGET NUMBER OF MEASUREMENTS = ",ntot
    
    allocate(dobs(nm),s_conf_all(nm,4),Wd_all(nm),mgrp(nm))
    if(i_flag) then
       allocate(dobsi(nm),Wd_alli(nm))
       do i=1,nm
          read(10,*,IOSTAT=ios) junk,s_conf_all(i,1:4),dobs(i),Wd_all(i),dobsi(i),Wd_alli(i),mgrp(i)
          if(ios .ne. 0) goto 1011
          if(Wd_all(i) <= 0 .or. Wd_alli(i) <= 0) then
             write(*,*) "NEGATIVE OR ZERO STD DEVIATION AT DATUM ",i
             write(*,*) "SETTING DEVIATION TO VERY LARGE "
             Wd_all(i)=1e15
          end if
          do j=1,4
             if(s_conf_all(i,j)>ne .or. s_conf_all(i,j)<0) then
                write(*,*) 'MEASURMENT ',i,'IS OUT OF ELECTRODE RANGE ... ABORTING'
                call crash_exit
             end if
          end do
          Wd_all(i)= 1/Wd_all(i)
          Wd_alli(i) = 1/Wd_alli(i)
       end do
       close(10)
       close(51)
       
    else
       do i=1,nm
          read(10,*,IOSTAT=ios) junk,s_conf_all(i,1:4),dobs(i),Wd_all(i),mgrp(i)
          if(ios .ne. 0) goto 1011
          if(Wd_all(i) <= 0) then
             write(*,*) "NEGATIVE OR ZERO STD DEVIATION AT DATUM ",i
             write(*,*) "SETTING DEVIATION TO VERY LARGE "
             Wd_all(i)=1e15
          end if
          do j=1,4
             if(s_conf_all(i,j)>ne .or. s_conf_all(i,j)<0) then
                write(*,*) 'MEASURMENT ',i,'IS OUT OF ELECTRODE RANGE ... ABORTING'
                call crash_exit
             end if
          end do
          Wd_all(i)= 1/Wd_all(i)
          
       end do
       close(51)  
    end if

    open(51,file='e4d.log',status='old',action='write',position='append')
    ngrp=maxval(mgrp)
    write(51,*) "THERE ARE ",ngrp," MEASUREMENT GROUPS"
 
   allocate(np_grp(ngrp))
   nbase=0
   np_grp=0
   do i=1,nm
      if(mgrp(i)==0) then
         nbase=nbase+1;
      else
         np_grp(mgrp(i)) = np_grp(mgrp(i))+1
      end if
   end do
   write(51,*) "THERE ARE ",nbase, " BASE MEASUREMENTS"
   write(51,*) "THERE ARE ",sum(np_grp), "CANDIDATE MEASUREMENTS"
   write(51,*) "THE LARGEST GROUP HAS",maxval(np_grp)," MEASUREMENTS"
   

   if(setgrp0) then
      nmt=nm
      nbase=nm
      
      allocate(s_conf(nm,4),Wd(nm))
      s_conf=s_conf_all
      Wd=Wd_all
   else
      nmt=nm
      nm=nbase
      allocate(s_conf(nm,4),Wd(nm))
      j=0
      do i=1,nmt
         if(mgrp(i)==0) then
            j=j+1
            s_conf(j,1:4)=s_conf_all(i,1:4)
            Wd(j)=Wd_all(i)
         end if
      end do
   end if

   read(10,*,IOSTAT=ios) nR_inds
   if(ios .ne.0) goto 1012
   if(nR_inds<1 .or. n_add<1) goto 1013
   
   !open(51,file='e4d.log',status='old',action='write',position='append')
   write(51,*) "THERE ARE ",nR_inds," RESOLUTION CONTROL POINTS"

   allocate(R_pts(nR_inds,3),R_inds(nR_inds))
   do i=1,nR_inds
      read(10,*,IOSTAT=ios) junk,R_pts(i,1),R_pts(i,2),R_pts(i,3)
      if(ios .ne. 0) goto 1014
   end do
   close(10)

   R_pts(:,1)=R_pts(:,1)-xorig
   R_pts(:,2)=R_pts(:,2)-yorig
   R_pts(:,3)=R_pts(:,3)-zorig
   close(51)
   return


1007 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading the number of electrodes in file: ",efile
    close(51)
    call crash_exit
    return
    
1008 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading the electrode position on line: ",i," of file: ",efile
    close(51)
    call crash_exit
    return

1009 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading trans.txt"
    close(51)
    call crash_exit
    return

1010 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading the number of measurement from file: ",efile
    close(51)
    call crash_exit
    return

1011 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading measurment: ",i," in file: ",efile
    close(51)
    call crash_exit
    return

1012 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading the number of resolution control points in file: ",efile
    close(51)
    call crash_exit
    return

1013 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There must be 1 or more resolution control points in file: ",efile
    write(51,*) "I read ",nR_inds
    close(51)
    call crash_exit
    return

1014 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "There was a problem reading resolution control point ",i," in file: ",efile
    close(51)
    call crash_exit
    return

1015 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    write(51,*) "THE NUMBER OF TARGET MEASUREMENTS: ",ntot
    write(51,*) "IS GREATER THAN THE NUMBER OF CANDIDATES: ",nm
    close(51)
    call crash_exit
    return


  end subroutine get_survey
  !_____________________________________________________________________________________

  !______________________________________________________________________________________
  subroutine get_cpts
    implicit none
    integer :: i,ii
    real, dimension(nelem) :: dist
    integer, dimension(1) :: indx
  
  
    call compute_mids

    do i=1,nR_inds
       dist=sqrt( (mids(:,1)-R_pts(i,1))**2 +(mids(:,2)-R_pts(i,2))**2 + (mids(:,3)-R_pts(i,3))**2)
      indx=minloc(dist)
      R_inds(i)=indx(1)
    end do
    
    open(50,file='psf_control_pts.txt',status='replace',action='write')
    do i=1,nR_inds
       write(50,*) i,mids(R_inds(i),1)+xorig,mids(R_inds(i),2)+yorig,mids(R_inds(i),3)+zorig
    end do
    close(50)

    deallocate(R_pts)

  end subroutine get_cpts
  !______________________________________________________________________________________

  !______________________________________________________________________________________
  subroutine compute_mids
    implicit none
    integer :: i

    if(allocated(mids)) deallocate(mids)
    allocate(mids(nelem,3))

    do i=1,nelem
       mids(i,1) = 0.25*sum(nodes(elements(i,1:4),1))
       mids(i,2) = 0.25*sum(nodes(elements(i,1:4),2))
       mids(i,3) = 0.25*sum(nodes(elements(i,1:4),3))
    end do
       
  end subroutine compute_mids
  !______________________________________________________________________________________

  !______________________________________________________________________________________
  function random_normal() result(fn_val)
    
    implicit none
    real :: fn_val
    real :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472
    real :: r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
    real :: half = 0.5
    
    !call RANDOM_NUMBER(fn_val)
    !return
	 !     Generate P = (u,v) uniform in rectangle enclosing acceptance region
    do
       call RANDOM_NUMBER(u)
       call RANDOM_NUMBER(v)
       v = 1.7156 * (v - half)
       
       !     Evaluate the quadratic form
       x = u - s
       y = abs(v) - t
       q = x**2 + y*(a*y - b*x)
       
       !     Accept P if inside inner ellipse
       if (q < r1) exit
       !     Reject P if outside outer ellipse
       if (q > r2) cycle
       !     Reject P if outside acceptance region
       if (v**2 < -4.0*log(u)*u**2) exit
    end do
    
    !     Return ratio of P's coordinates as the normal deviate
    fn_val = v/u
    return
    
  end function random_normal
  !_________________________________________________________________________________
  

  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    CALL SYSTEM_CLOCK(COUNT=clock)
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed
  

  !_________________________________________________________________________________
  subroutine write_survey(it)
    implicit none
    integer :: it,i,j
    character*40 :: sfile
    
    write(sfile,'(A,I0,A)') 'opt_survey_',it,'.srv'
    open(25,file=trim(sfile),status='replace',action='write')
    write(25,'(I0)') ne
    do i=1,ne
       write(25,'(I5,3f10.3,I5)') i,e_pos(i,1)+xorig,e_pos(i,2)+yorig, &
                                    e_pos(i,3)+zorig,int(e_pos(i,4))
    end do
    write(25,*)
   
    write(25,'(I0)') nm
    j=0
    do i=1,nmt
       if(mgrp(i)==0 .or. sused(mgrp(i))) then
          j=j+1
          write(25,'(7I10,I10,f10.3)') j,s_conf_all(i,1:4),-999,-999,mgrp(i),-sg_rank(mgrp(i))
       end if
    end do
    close(25)
  end subroutine write_survey
  !_________________________________________________________________________________

  !_________________________________________________________________________________
  subroutine get_volume
    implicit none
    integer :: kk
    real :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
 
    allocate(vlm(nelem))
    do kk=1,nelem
       
       !!get the 4 nodes for this element 
       x1 = nodes(elements(kk,1),1)
       y1 = nodes(elements(kk,1),2)
       z1 = nodes(elements(kk,1),3)
       
       x2 = nodes(elements(kk,2),1)
       y2 = nodes(elements(kk,2),2)
       z2 = nodes(elements(kk,2),3)
       
       x3 = nodes(elements(kk,3),1)
       y3 = nodes(elements(kk,3),2)
       z3 = nodes(elements(kk,3),3) 
       
       x4 = nodes(elements(kk,4),1)
       y4 = nodes(elements(kk,4),2)
       z4 = nodes(elements(kk,4),3)
       
      
       !!get the volume of this element
       !call get_vol (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,vlm(kk)) 
    end do
  end subroutine get_volume
  !_________________________________________________________________________________


end module mod_res
