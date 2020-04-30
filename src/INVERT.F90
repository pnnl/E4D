module invert

  use vars
  use fmm_vars
  use master
  use mod_con
  use input
  implicit none

 
contains
  
  !_____________________________________________________________________________________
  subroutine pcgls
    implicit none
    
    integer :: i,j,ncon,k,indefinite,unstable,info,rbi,nbr
    real :: et,t1,t2,cs,ce
    real :: norms0,xmax,normx,resNE,resNE_old
    real :: alpha,gamma,delta,norms,gamma1,gbeta
    real, dimension(:), allocatable :: q,r,b
    real, dimension(nelem) :: s,p,xpar
    
    call cpu_time(Cstart)
    cs=Cstart

    call send_command(20)
    if(im_fmm) then
       write(*,*) "FMM: STARTING INVERSION"
    else
       write(*,*) "E4D: STARTING INVERSION"
    end if

    !ncon=ccount !maxval(wrows)   
    xpar = 0
    

    if(im_fmm) then
       !speed is slowness squared for the forward fmm computations
       !change speed to slowness for the inverse update
       speed = sqrt(speed)
       nm=nm_fmm
       allocate(q(nm+ccount),r(nm+ccount),b(nm+ccount))
       b = 0
    else
       allocate(q(nm+ccount),r(nm+ccount),b(nm+ccount))
       b = 0
       if(invi) then
          sigmai = log(sigmai)
          sigma = log(sigma)
       else
          sigma = log(sigma)
          refsig=log(refsig)
          if(allocated(prefsig)) prefsig=log(prefsig)
       end if
    end if
 
    if(invi) then
       do i=1,nm
          b(i) = (Wdi(i)*Wd_cull(i)*(dobsi(i)-dpredi(i)))
       end do
    else
       do i=1,nm
          b(i) = (Wd(i)*Wd_cull(i)*(dobs(i)-dpred(i)))
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
    r = b

    !!s = matmul(transpose(J),b)
    !!s = pmatvec2(nm,nelem,ncon,ccount,wrows,wcols,Wm,r(1:nm+ncon),beta,npar,par_map,nec,ex_cols,ex_vals
    s = pmatvec2II(r)  
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
    write(*,"(A,g10.4,A)") "  STARTING INNER ITERATIONS AT: ",ce-cs," seconds"
    do i=1,max_initer

       if( info .ne. 0) then
          exit
       end if
  
       q = pmatvec1II(p)
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
       s = pmatvec2II(r)
 
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
       if(im_fmm) then
          write(*,"(A,I3.3,A,f10.3,A,f10.4)") "  FMM: ITERATION ",i," FINISHED AT: ",&
               (ce-cs)/60," minutes  ",resNE
       else
           write(*,"(A,I3.3,A,f10.3,A,f10.4)") "  E4D: ITERATION ",i," FINISHED AT: ",&
               (ce-cs)/60," minutes  ",resNE
        end if
       
       if( abs((resNE_old-resNe)/resNE_old) < delta_initer .and. i > min_initer) then
          exit
       end if
     
    end do
     
    if(im_fmm) then
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
    end if
 100 continue
  
    call send_command(21)
    call cpu_time(Cend)
    etm=Cend-Cstart
 
    return
  end subroutine pcgls
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

       if(smetric(rbi,2).eq.12) then      
          pmatvec1II(nm+i) = pmatvec1II(nm+i) + sqrt(beta)*Wm(i)*x(rblock(i,1))*cg_wts(rblock(i,2),1)
          do j=1,4
             nbr = neighbors(rblock(i,1),j)
             if(nbr>0) then
                pmatvec1II(nm+i) = pmatvec1II(nm+i) + &
                     sqrt(beta)*Wm(i)*x(nbr)*cg_wts(rblock(i,2),j+1)
             end if
          end do    
       end if

10     continue
    end do
   
    
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

       if(smetric(rbi,2).eq.12) then
          pmatvec2II(rblock(i,1)) = pmatvec2II(rblock(i,1)) + &
               sqrt(beta)*Wm(i)*r(nm+i)*cg_wts(rblock(i,2),1)
          do j=1,4
             nbr = neighbors(rblock(i,1),j)
             if(nbr .gt. 0) then
                pmatvec2II(nbr) = pmatvec2II(nbr) + &
                     sqrt(beta)*Wm(i)*r(nm+i)*cg_wts(rblock(i,2),j+1)
             end if
          end do
       end if

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



  
  
end module invert
  
