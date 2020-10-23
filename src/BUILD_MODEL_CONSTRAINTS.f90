module mod_con

  use vars
  use fmm_vars
  use input
  use reorder_mesh
 

  implicit none
  integer, dimension(:), allocatable :: wrows,wcols,par_map,master_ranks
  real, dimension(:), allocatable :: Wm,sigma_par,rsigma_par,sigi,v,Wmw
  real*8, dimension(:), allocatable :: WA,BM,tnod,wts4cg
  real, dimension(3) :: CV,CK,CJ
  integer :: ccount,npar,ji_count,cur_el,ncgrad,nnz_cgrad,nphys
  integer, dimension(:,:), allocatable :: neighbors,block4cgrad
  logical :: itst=.true.
  logical :: iitst=.false.
  logical, dimension(:),allocatable :: use4cgrad
contains

  !_____________________________________________________________________
  subroutine build_WmII
    implicit none
    integer :: n_homo,i,j,k,nzmax,n_homo_elem,count,row,col,i1,i2,lzone
    integer :: znn,nbr,nbr2,x1,x2,x3,z1,z2,z3,ii,nzn,nelem_tmp,cnt,jnk
    integer, dimension(:), allocatable :: homo_par
    integer, dimension(4) :: neighbors_temp          
    real :: midxi,midyi,midzi,midxn,midyn,midzn,rx,ry,rz,r,x,rho,eps
    logical :: ffound, extfil
    character*80 :: extfile
    
    if(invi .and. itst) then
       itst = .false.
       iitst = .true.
       if(allocated(Wm)) deallocate(Wm)

    elseif(iitst .and. .not. invi) then
       itst = .true.
       iitst = .false.
       if(allocated(Wm)) deallocate(Wm)
    end if
   
    if(allocated(Wm)) then
       do i=1,ccount
          call comp_Wm(i,Wm(i))
       end do
       return
    else
      
       if(allocated(neighbors)) deallocate(neighbors)
       if(allocated(rblock)) deallocate(rblock)
       
       !!Check to see if there are NN regularized zones. If so, read neighbors
       !if(sum(smetric(:,2))>nrz) then          
       allocate(neighbors(nelem,4))
       do i=1,80
          if(mshfile(i:i)=='.') then
             open(21,file=mshfile(1:i+1)//".neigh",status='old',action='read')
                         
             if(allocated(element_map)) then
                !if element_map is allocated then there are inactive elements 
                !and we need to re_map the element neighbors
                read(21,*) nelem_tmp
                cnt = 0
                do j=1,nelem_tmp
                   read(21,*) jnk,neighbors_temp(1:4)
                   if(element_map(j) .ne. 0) then
                      cnt=cnt+1
                      do k=1,4
                         if(neighbors_temp(k)>0) then
                            neighbors(cnt,k) = element_map(neighbors_temp(k))
                         else
                            neighbors(cnt,k) = neighbors_temp(k)
                         end if
                      end do
                   end if
                end do
                close(21)

             else
                read(21,*) 
                do j=1,nelem
                   read(21,*) k,neighbors(j,1:4)
                end do
                close(21)
             end if
             exit
          end if
       end do
       nzn = maxval(zones)
     
       !count the number of constraint equations
       row=0
   
       do i=1,nelem
          do j=1,nrz
             !!12 and 13 are joint inversion constraints, skip them for now
             if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) goto 10
             if(smetric(j,1)==nzn+1 .and. i==1) then
                !this constraint block is an external constraint block
                !read the number of constraint equations but only do this once
                open(11,file='external_files.txt',status='old',action='read',IOSTAT=ios)
                if(ios .ne. 0) then
                   write(*,*) 'Did not find the external constraint list file external_confiles.txt'
                   write(*,*) 'Not implementing external constraint'
                   close(11)
                else
                   !loop over the list file 
                   ffound = .false.
                   do while(.not.ffound)
                      read(11,*,IOSTAT=ios) ii,extfile
                      if(ii==j) then
                         ffound = .true.
                         close(11)
                      end if
                   end do
                   if(ffound) then
                      open(11,file=extfile,status='old',action='read')
                      read(11,*) ii
                      ccount=ccount+ii
                      !write(*,*) 'adding ',ii,' constraints for block ',j
                   end if
                end if
             else
                if(smetric(j,2) > 0 ) then
                   if(smetric(j,1)==zones(i)) then
                      if(smetric(j,2)==3 .or. smetric(j,2)==4 ) then
                         ccount=ccount+1
                         
                      else
                         do k=1,4
                            nbr=neighbors(i,k)
                            
                            if(zones(nbr) .ne. zones(i)) then
                               do ii=2,zone_links(j,1)+1
                                  if(zone_links(j,ii)==zones(nbr)) then
                                     ccount=ccount+1
                                  end if
                               end do
                            else
                               ccount=ccount+1
                            end if
                            
                         end do
                      end if
                   end if
                end if
             end if
10           continue
          end do
       end do

       !check to see if this is a joint inversion and there are
       !cross-gradient constraints specified. If so, setup the
       !cross-gradient inversion
       if(cgmin_flag(1) .and. cgmin_flag(2)) then
          do j = 1,nrz
             if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) then
                if(.not. allocated(WA)) then
                   call setup_crossgrad
                end if
                exit
             end if
          end do
       end if
       
!!$       !count the constraints for the joint inversion if applicable
!!$       if(cgmin_flag(1)) then
!!$          ji_count = 0
!!$          do i=1,nelem
!!$             do j=1,nrz
!!$                if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) then
!!$                   if(smetric(j,1) .eq. zones(i)) then
!!$                      !one constraint for each component of the 
!!$                      !cross-gradient vector
!!$                      ji_count=ji_count+3
!!$                      ccount=ccount+3
!!$                   end if
!!$                end if
!!$             end do
!!$          end do
!!$          allocate(cg_wts(ji_count,5))
!!$          cg_wts = 0
!!$          cur_el = 0
!!$       end if
 
       allocate(Wm(ccount),rblock(ccount,3))
       Wm = 0.0
       ccount=0
       ji_count=0
       row=0
       col=0
    
       rblock=0
       do i=1,nelem
          do j=1,nrz
             if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) goto 20
             
             if(smetric(j,1)==nzn+1 .and. i==1) then
                !this constraint block is an external constraint block
                !read the constraints only once
                open(11,file='external_files.txt',status='old',action='read',IOSTAT=ios)
                if(ios .ne. 0) then   
                   close(11)
                else
                   ffound = .false.
                   do while(.not.ffound)
                      read(11,*,IOSTAT=ios) ii,extfile
                      if(ii==j) then
                         ffound = .true.
                         close(11)
                         !write(*,*) 'Adding external constraints for block: ',j
                      end if
                   end do
                   if(ffound) then
                      open(11,file=extfile,status='old',action='read')
                      read(11,*) ii
                      do k=1,ii
                         ccount=ccount+1
                         read(11,*) rblock(ccount,1:2)
                         rblock(ccount,3) = j
                      end do
                      close(11)
                      write(*,*) ' ADDING EXTERNAL CONSTRAINTS FROM FILE: ',trim(extfile),' FOR REG BLOCK: ',j
                   end if
                end if
		
             else      
           
                if(smetric(j,2)>0) then
                   if(smetric(j,1)==zones(i)) then
                      if(smetric(j,2)==3 .or. smetric(j,2)==4) then
                         ccount=ccount+1
                         rblock(ccount,1)=i
                         rblock(ccount,3)=j
                         
                      else
                         do k=1,4
                            nbr=neighbors(i,k)
                            
                            if(zones(nbr) .ne. zones(i)) then
                               do ii=2,zone_links(j,1)+1
                                  if(zone_links(j,ii)==zones(nbr) ) then
                                     ccount=ccount+1
                                     rblock(ccount,1)=i
                                     rblock(ccount,2)=nbr
                                     rblock(ccount,3)=j
                                  end if
                               end do
                            else
                               if(smetric(j,2)<9 .or. smetric(j,2)>10) then
                                  ccount=ccount+1
                                  rblock(ccount,1)=i
                                  rblock(ccount,2)=nbr
                                  rblock(ccount,3)=j
                               end if
                            end if
                            
                         end do
                      end if
                      
                   end if
                end if
             end if
20           continue
          end do
       end do
    end if
    
  
    
    !count the joint inversion constraints if necessary
!!$    if(cgmin_flag(1)) then
!!$       do i=1,nelem
!!$          do j=1,nrz
!!$             if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) then
!!$                if(smetric(j,1) .eq. zones(i)) then
!!$                   !one constraint for each component of the 
!!$                   !cross-gradient vector
!!$                   do k=1,3
!!$                      ccount=ccount+1
!!$                      ji_count = ji_count+1
!!$                      rblock(ccount,1) = i
!!$                      rblock(ccount,2) = ji_count
!!$                      rblock(ccount,3) = j
!!$                   end do
!!$                end if
!!$             end if
!!$          end do
!!$       end do
!!$    end if
    
    !set the on/off flag for element estimation
    if(.not.allocated(J_on_off)) then
       allocate(J_on_off(nelem))
       J_on_off = .true.
       do i=1,nrz

          if(smetric(i,2)==0 .and. smetric(i,1).ne.nzn+1) then
             do j=1,nelem
                if(smetric(i,1)==zones(j)) then
                   J_on_off(j) = .false.
                end if
             end do
          end if

          if(smetric(i,2)==0 .and. smetric(i,1)==nzn+1) then
             !these are specified in and external constraint file
             open(11,file='external_files.tmp',status='old',action='read')
             ffound = .false.
             do while(.not.ffound)
                read(11,*,IOSTAT=ios) ii,extfile
                if(ii==i) then
                   ffound = .true.
                   close(11)
                end if
             end do
             if(ffound) then
                open(11,file=extfile,status='old',action='read')
                read(11,*) ii
                do k=1,ii
                   read(11,*) j
                   J_on_off(j) = .false.
                end do
                close(11)
             end if
          end if

       end do
    end if
   
    do i=1,ccount    
       call comp_Wm(i,Wm(i))
    end do
 
  end subroutine build_WmII
  !_____________________________________________________________________
  

  !_____________________________________________________________________
  subroutine comp_Wm(indx,wm)
    implicit none
    integer, intent(in) :: indx
    real,intent(out) :: wm
    integer :: rbi,eli,eln,i
    real :: X, mn,sd
    real :: midxi,midyi,midzi,midxn,midyn,midzn,rx,ry,rz,r,rho,eps
    
    
    rbi=abs(rblock(indx,3))
    select case (smetric(rbi,2))
 
    case(1) 
       if(im_fmm) then
           X=(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2))))
       else
         if(invi) then
            X=(log(phase(rblock(indx,1)))-log(phase(rblock(indx,2))))
         else
            X=(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2))))
         end if
       end if
       
    case(2)
       if(im_fmm) then
           X=abs(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2))))
       else
         if(invi) then
            !X=abs(log(sigmai(rblock(indx,1)))-log(sigmai(rblock(indx,2))))
             X=abs(log(phase(rblock(indx,1))) - log(phase(rblock(indx,2))))
         else
             X=abs(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2))))
         end if
       end if
       
    case(3)
       select case(smetric(rbi,3))
       case(0)
          if(im_fmm) then
             X=sqrt(velocity(rblock(indx,1)))-(C_targ(rbi))
          else
             if(invi) then
                X=log(phase(rblock(indx,1)))-(C_targ(rbi))
             else
                X=log(sigma(rblock(indx,1)))-(C_targ(rbi))
             end if
          end if
       case(1)
          if(im_fmm) then
             X=sqrt(velocity(rblock(indx,1)))-(refsig(rblock(indx,1)))
          else
             if(invi) then
                X=log(phase(rblock(indx,1)))-log(refsig(rblock(indx,1)))
             else
                X=log(sigma(rblock(indx,1)))-log(refsig(rblock(indx,1)))
             end if
          end if   
       case(2)
          if(im_fmm) then
             X=sqrt(velocity(rblock(indx,1)))-(prefsig(rblock(indx,1)))
	  else
             if(invi) then
                X=log(phase(rblock(indx,1)))-log(prefsig(rblock(indx,1)))
             else
                X=log(sigma(rblock(indx,1)))-log(prefsig(rblock(indx,1)))
             end if
          end if
       end select
      
    case(4)
       select case(smetric(rbi,3))
       case(0)
       	  if(im_fmm) then
       	     X=abs(sqrt(velocity(rblock(indx,1)))-(C_targ(rbi)))
       	  else
             if(invi) then
                X=abs(log(phase(rblock(indx,1)))-(C_targ(rbi)))
             else 
                X=abs(log(sigma(rblock(indx,1)))-(C_targ(rbi)))
             end if
          end if
       case(1)
          if(im_fmm) then
             X=abs(sqrt(velocity(rblock(indx,1)))-(refsig(rblock(indx,1))))
          else
             if(invi) then
                X=abs(log(phase(rblock(indx,1)))-log(refsig(rblock(indx,1))))
             else
                X=abs(log(sigma(rblock(indx,1)))-log(refsig(rblock(indx,1))))
             end if
          end if
       case(2)
          if(im_fmm) then
             X=abs(sqrt(velocity(rblock(indx,1)))-(prefsig(rblock(indx,1))))
          else
             if(invi) then
                X=abs(log(phase(rblock(indx,1)))-log(prefsig(rblock(indx,1))))
             else
                X=abs(log(sigma(rblock(indx,1)))-log(prefsig(rblock(indx,1))))
             end if
          end if   
       end select
     
    case(5)
       if(im_fmm) then
          X=(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2)))) 
       else
          if(invi) then
             X=(log(phase(rblock(indx,1)))-log(phase(rblock(indx,2)))) 
          else
             X=(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2)))) 
          end if
       end if
       eli=rblock(indx,1)
       eln=rblock(indx,2)
       midxi=0.25*sum(nodes(elements(eli,1:4),1))
       midyi=0.25*sum(nodes(elements(eli,1:4),2))
       midzi=0.25*sum(nodes(elements(eli,1:4),3))  
       midxn=0.25*sum(nodes(elements(eln,1:4),1))
       midyn=0.25*sum(nodes(elements(eln,1:4),2))
       midzn=0.25*sum(nodes(elements(eln,1:4),3))
       !rx=((midxi-midxn)**2)
       !ry=((midyi-midyn)**2)
       !rz=((midzi-midzn)**2)
       !r=(rx+ry+rz)             
       rx=((midxi-midxn))
       ry=((midyi-midyn))
       rz=((midzi-midzn))
       r=sqrt((rx**2+ry**2+rz**2))
       rx = rx/r
       ry = ry/r
       rz = rz/r
       
    case(6)
       if(im_fmm) then
          X=abs(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2))))  
       else
          if(invi) then
             X=abs(log(phase(rblock(indx,1)))-log(phase(rblock(indx,2))))
          else
             X=abs(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2))))
          end if
       end if
       eli=rblock(indx,1)
       eln=rblock(indx,2)
       midxi=0.25*sum(nodes(elements(eli,1:4),1))
       midyi=0.25*sum(nodes(elements(eli,1:4),2))
       midzi=0.25*sum(nodes(elements(eli,1:4),3))  
       midxn=0.25*sum(nodes(elements(eln,1:4),1))
       midyn=0.25*sum(nodes(elements(eln,1:4),2))
       midzn=0.25*sum(nodes(elements(eln,1:4),3))
              
       rx=((midxi-midxn))
       ry=((midyi-midyn))
       rz=((midzi-midzn))
       r=sqrt((rx**2+ry**2+rz**2))
       rx = rx/r
       ry = ry/r
       rz = rz/r

       !!extra for fracture testing
       !rx=abs((midxi-midxn))
       !ry=abs((midyi-midyn))
       !rz=abs((midzi-midzn))
       !rx = rx/r
       !ry = ry/r
       !rz = rz/r
       
       !r=sqrt(midxi**2 + midyi**2 + midzi**2)
       !midxi=midxi/r
       !midyi=midyi/r
       !midzi=midzi/r
       !r=rx*midxi + ry*midyi + rz*midzi
       
    case(7)
       select case(smetric(rbi,3))
       case(0)
          if(im_fmm) then
             X=(sqrt(velocity(rblock(indx,1)))-C_targ(rbi)) - &
               (sqrt(velocity(rblock(indx,2)))-C_targ(rbi))
          else
             if(invi) then
                X=(log(phase(rblock(indx,1)))-C_targ(rbi)) - &
                     (log(phase(rblock(indx,2)))-C_targ(rbi))
             else
                X=(log(sigma(rblock(indx,1)))-C_targ(rbi)) - &
                     (log(sigma(rblock(indx,2)))-C_targ(rbi))
             end if
          end if
       case(1)
          if(im_fmm) then
              X=(sqrt(velocity(rblock(indx,1)))-(refsig(rblock(indx,1)))) - &
                (sqrt(velocity(rblock(indx,2)))-(refsig(rblock(indx,2)))) 
          else
             if(invi) then
                X=(log(phase(rblock(indx,1)))-log(refsig(rblock(indx,1)))) - &
                     (log(phase(rblock(indx,2)))-log(refsig(rblock(indx,2))))
              
             else
                X=(log(sigma(rblock(indx,1)))-log(refsig(rblock(indx,1)))) - &
                     (log(sigma(rblock(indx,2)))-log(refsig(rblock(indx,2))))
             end if
          end if
      case(2)
         if(im_fmm) then
            X=(sqrt(velocity(rblock(indx,1)))-(prefsig(rblock(indx,1)))) - &
              (sqrt(velocity(rblock(indx,2)))-(prefsig(rblock(indx,2))))         
         else
            if(invi) then
               X=(log(phase(rblock(indx,1)))-log(prefsig(rblock(indx,1)))) - &
                  (log(phase(rblock(indx,2)))-log(prefsig(rblock(indx,2))))
            else
               X=(log(sigma(rblock(indx,1)))-log(prefsig(rblock(indx,1)))) - &
                 (log(sigma(rblock(indx,2)))-log(prefsig(rblock(indx,2))))
            end if
         end if
         
      end select

    case(8)
       select case(smetric(rbi,3))
       case(0)
          if(im_fmm) then
             X=abs((sqrt(velocity(rblock(indx,1)))-C_targ(rbi)) - &
                   (sqrt(velocity(rblock(indx,2)))-C_targ(rbi)))
          else
             if(invi) then
                X=abs((log(phase(rblock(indx,1)))-C_targ(rbi)) - &
                     (log(phase(rblock(indx,2)))-C_targ(rbi)))
             else
                X=abs((log(sigma(rblock(indx,1)))-C_targ(rbi)) - &
                     (log(sigma(rblock(indx,2)))-C_targ(rbi)))
             end if
          end if
       case(1)
          if(im_fmm) then
             X=abs((sqrt(velocity(rblock(indx,1)))-(refsig(rblock(indx,1)))) - &
                   (sqrt(velocity(rblock(indx,2)))-(refsig(rblock(indx,2)))))
          else
             if(invi) then
                X=abs((log(phase(rblock(indx,1)))-log(refsig(rblock(indx,1)))) - &
                     (log(phase(rblock(indx,2)))-log(refsig(rblock(indx,2)))))
             else
                X=abs((log(sigma(rblock(indx,1)))-log(refsig(rblock(indx,1)))) - &
                     (log(sigma(rblock(indx,2)))-log(refsig(rblock(indx,2)))))
             end if
          end if
       case(2)
          if(im_fmm) then
             X=abs((sqrt(velocity(rblock(indx,1)))-(prefsig(rblock(indx,1)))) - &
                  (sqrt(velocity(rblock(indx,2))) -(prefsig(rblock(indx,2)))))
          else
             if(invi) then
                X=abs((log(phase(rblock(indx,1)))-log(prefsig(rblock(indx,1)))) - &
                     (log(phase(rblock(indx,2)))-log(prefsig(rblock(indx,2)))))
             else
                X=abs((log(sigma(rblock(indx,1)))-log(prefsig(rblock(indx,1)))) - &
                     (log(sigma(rblock(indx,2)))-log(prefsig(rblock(indx,2)))))
             end if
          end if
       end select

    case(9)
       if(im_fmm) then
          X=(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2))))
       else
         if(invi) then
            X=(log(phase(rblock(indx,1)))-log(phase(rblock(indx,2))))
         else
            X=(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2))))
         end if
       end if
       
    case(10)
       if(im_fmm) then
          X=abs(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2))))
       else
          if(invi) then
             X=abs(log(phase(rblock(indx,1)))-log(phase(rblock(indx,2))))
          else
             X=abs(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2))))
          end if
       end if
       
    case(11) !test case for horizontal radial regularization TCJ 11/04/15
       if(im_fmm) then
           X=abs(sqrt(velocity(rblock(indx,1)))-sqrt(velocity(rblock(indx,2))))
       else
          if(invi) then
             X=abs(log(phase(rblock(indx,1)))-log(phase(rblock(indx,2))))
          else
             X=abs(log(sigma(rblock(indx,1)))-log(sigma(rblock(indx,2))))
          end if
       end if
       eli=rblock(indx,1)
       eln=rblock(indx,2)
       
       midxi=0.25*sum(nodes(elements(eli,1:4),1))!-zwts(rbi,2)
       midyi=0.25*sum(nodes(elements(eli,1:4),2))!-zwts(rbi,3)
       !midzi=0.25*sum(nodes(elements(eli,1:4),3))-zwts(rbi,4)  
       midxn=0.25*sum(nodes(elements(eln,1:4),1))!-zwts(rbi,2)
       midyn=0.25*sum(nodes(elements(eln,1:4),2))!-zwts(rbi,3)
       !midzn=0.25*sum(nodes(elements(eln,1:4),3))-zwts(rbi,4)
       !get the normal from neighbor to me
       midxn=midxn-midxi
       midyn=midyn-midyi
       r=sqrt(midxn**2 + midyn**2)
       if(r==0) then
          r=1
       else   
          midxn=midxn/r
          midyn=midyn/r

          r=sqrt(midxi**2 + midyi**2);
          if(r==0) then
             r=1
          else
             midxi=midxi/r
             midyi=midyi/r
             r=(midxn*midxi + midyn*midyi)**2
          end if
       end if
 
    case(12)
       !build the joint inversion inversion constraint
       !for element rblock(indx,1). Note rblock(indx,2)
       !is the row of cg_wts that holds this constraint.
       call build_cgwt(rblock(indx,1),rblock(indx,2))
       if(im_fmm) then
          X=sqrt(velocity(rblock(indx,1)))*cg_wts(rblock(indx,2),1)
          do i=1,4
             eln=neighbors(rblock(indx,1),i)
             X=X+sqrt(velocity(eln))*cg_wts(rblock(indx,2),i+1)
          end do
       else
          X=log(sigma(rblock(indx,1)))*cg_wts(rblock(indx,2),1)
          do i=1,4
             eln=neighbors(rblock(indx,1),i)
             X=X+log(sigma(eln))*cg_wts(rblock(indx,2),i+1)
          end do
       end if	
       X=abs(X)

    case DEFAULT
       
    end select
    
    mn=Fw_parm(rbi,1)
    sd=Fw_parm(rbi,2)
    
    select case (Fw_type(rbi))
       
    case(1)
       wm = .5* (1- erf( (X-mn)/sqrt(2*sd**2)))
    case(2)
       wm = .5* (1+ erf( (X-mn)/sqrt(2*sd**2)))
       
    case(3)
       wm = 1-exp(-((X-mn)**2)/(2*sd**2))
        
    case(4)
       wm = exp(-((X-mn)**2)/(2*sd**2))

    case(5)
       if((X-mn)<0) then
          wm=sd**(-2)
       else
          wm=sd**(2) * ((X-mn)**2 + sd**2)**(-2)
       end if

    case(6)
       if((X-mn)>0) then
          wm=sd**(-2)
       else
          wm=sd**(2) * ((X-mn)**2 + sd**2)**(-2)
       end if
      
    case DEFAULT
    end select
    
    if(smetric(rbi,2)==5 .or. smetric(rbi,2)==6) then
       wm = (1 - abs( rx*zwts(rbi,2) + ry*zwts(rbi,3) + rz*zwts(rbi,4)))**2
      
    end if
    if(smetric(rbi,2)==11) then
       wm=r*wm
    end if
    wm=wm*zwts(rbi,1)

    !! test adjustment with face area
    !wm=wm*Wmw(indx)	    

  end subroutine comp_Wm
  !_____________________________________________________________________

  !_____________________________________________________________________
  subroutine build_BM
    implicit none
    integer :: i,j,nbr,iphys,incr_b,ind1,ind2
    real*8, dimension(3) :: grad
    !TEMP-DBG
    !integer :: row
    !integer, dimension(5) :: itemp

    ! ! DBG
    ! if (im_fmm) then
    !    open(103,file='B2.txt',status='replace',action='write')
    ! else
    !    open(103,file='B1.txt',status='replace',action='write')
    ! endif
    ! row = 0
    
    
    incr_b=0
    do iphys=1,nphys-1
       select case (block4cgrad(iphys,2))

       case(1) ! B2
          ncgrad=0
          do i=1,nelem
             if(use4cgrad(i)) then
                ncgrad=ncgrad+1
                grad(1)=WA(1+(ncgrad-1)*15)*dble(sigma(i))
                grad(2)=WA(2+(ncgrad-1)*15)*dble(sigma(i))
                grad(3)=WA(3+(ncgrad-1)*15)*dble(sigma(i))
                ! DGB: Store col index
                !itemp(1) = i
                do j=1,4
                   nbr=neighbors(i,j)
                   ! DBG: store col index
                   !itemp(j+1) = nbr
                   
                   grad(1)=grad(1)+WA(1+j*3+(ncgrad-1)*15)*&
                       dble(sigma(nbr))
                   grad(2)=grad(2)+WA(2+j*3+(ncgrad-1)*15)*&
                        dble(sigma(nbr))
                   grad(3)=grad(3)+WA(3+j*3+(ncgrad-1)*15)*&
                        dble(sigma(nbr))
                end do
                
!if (iter.eq.25.and.incr_b.eq.15*98990) print*,"ielem, grad_m1:",i,grad
    
                do j=0,4
                   BM(1+j*3+incr_b)=grad(2)*WA(3+j*3+(ncgrad-1)*15)-&
                                    grad(3)*WA(2+j*3+(ncgrad-1)*15)
                   BM(2+j*3+incr_b)=grad(3)*WA(1+j*3+(ncgrad-1)*15)-&
                                    grad(1)*WA(3+j*3+(ncgrad-1)*15)
                   BM(3+j*3+incr_b)=grad(1)*WA(2+j*3+(ncgrad-1)*15)-&
                                    grad(2)*WA(1+j*3+(ncgrad-1)*15)
                end do

                ! ! DBG: To write sparse B matrix
                ! do j=1,3
                !    row = row + 1
                !    write(103,*) row,itemp(1),BM( 1+(j-1)+incr_b)
                !    write(103,*) row,itemp(2),BM( 4+(j-1)+incr_b)
                !    write(103,*) row,itemp(3),BM( 7+(j-1)+incr_b)
                !    write(103,*) row,itemp(4),BM(10+(j-1)+incr_b)
                !    write(103,*) row,itemp(5),BM(13+(j-1)+incr_b)                   
                ! enddo
                
                incr_b=incr_b+15
             end if
          end do

       case(2) ! B1
          ncgrad=0
          do i=1,nelem
             if(use4cgrad(i)) then
                ncgrad=ncgrad+1
                grad(1)=WA(1+(ncgrad-1)*15)*dble(velocity(i))
                grad(2)=WA(2+(ncgrad-1)*15)*dble(velocity(i))
                grad(3)=WA(3+(ncgrad-1)*15)*dble(velocity(i))
                ! DGB: Store col index
                !itemp(1) = i
                do j=1,4
                   nbr=neighbors(i,j)
                   ! DBG: store col index
                   !itemp(j+1) = nbr
                   
                   grad(1)=grad(1)+WA(1+j*3+(ncgrad-1)*15)*&
                       dble(velocity(nbr))
                   grad(2)=grad(2)+WA(2+j*3+(ncgrad-1)*15)*&
                        dble(velocity(nbr))
                   grad(3)=grad(3)+WA(3+j*3+(ncgrad-1)*15)*&
                        dble(velocity(nbr))
                end do

!if (iter.eq.25.and.incr_b.eq.15*98990) print*,"ielem, grad_m2:",i,grad                
               
                do j=0,4
                   BM(1+j*3+incr_b)=grad(3)*WA(2+j*3+(ncgrad-1)*15)-&
                                    grad(2)*WA(3+j*3+(ncgrad-1)*15)
                   BM(2+j*3+incr_b)=grad(1)*WA(3+j*3+(ncgrad-1)*15)-&      ! CAN BE VECTORIZED
                                    grad(3)*WA(1+j*3+(ncgrad-1)*15)
                   BM(3+j*3+incr_b)=grad(2)*WA(1+j*3+(ncgrad-1)*15)-&
                                    grad(1)*WA(2+j*3+(ncgrad-1)*15)
                end do

                ! ! DBG: To write sparse B matrix
                ! do j=1,3
                !    row = row + 1
                !    write(103,*) row,itemp(1),BM( 1+(j-1)+incr_b)
                !    write(103,*) row,itemp(2),BM( 4+(j-1)+incr_b)
                !    write(103,*) row,itemp(3),BM( 7+(j-1)+incr_b)
                !    write(103,*) row,itemp(4),BM(10+(j-1)+incr_b)
                !    write(103,*) row,itemp(5),BM(13+(j-1)+incr_b)                   
                ! enddo                
                
                incr_b=incr_b+15
             end if
          end do

       end select

       !if(block4cgrad(iphys,1)>block4cgrad(iphys,2)) then
       !   ind1 = 1+(iphys-1)*nnz_cgrad
       !   ind2 = iphys*nnz_cgrad
       !   BM(ind1:ind2) = -BM(ind1:ind2)
       !end if
    end do

    ! DBG
    !print*,ncgrad,nelem
    !close(103)    
    
    do iphys=1,nphys-1
       do i=1,ncgrad
          do j=1,15
             BM(j+(i-1)*15+(iphys-1)*15*ncgrad) = &
                  BM(j+(i-1)*15+(iphys-1)*15*ncgrad)*wts4cg(i)
          end do
       end do
    end do
  end subroutine build_BM
  !_____________________________________________________________________

  !_____________________________________________________________________
  subroutine build_V
    implicit none
    integer :: i,j,nbr,iphys,incr,incr_b

    incr=0
    incr_b=0
    do iphys=1,nphys-1
       select case (block4cgrad(iphys,1))

       case(1)
          do i=1,nelem
             if(use4cgrad(i)) then
                tnod(1+incr)=BM(1+incr_b)*dble(sigma(i))
                tnod(2+incr)=BM(2+incr_b)*dble(sigma(i))
                tnod(3+incr)=BM(3+incr_b)*dble(sigma(i))

!if (iter.eq.25.and.incr.eq.3*98990) then
!   print*,"Sigma Center:",sigma(i)                
!   print*,'i:',0.25*dble(sum(nodes(elements(i,1:4),1))),0.25*dble(sum(nodes(elements(i,1:4),2))), &
!        0.25*dble(sum(nodes(elements(i,1:4),3)))
!endif

                do j=1,4
                   nbr=neighbors(i,j)
                   tnod(1+incr)=tnod(1+incr)+BM(1+j*3+incr_b)*&
                        dble(sigma(nbr))
                   tnod(2+incr)=tnod(2+incr)+BM(2+j*3+incr_b)*&
                        dble(sigma(nbr))
                   tnod(3+incr)=tnod(3+incr)+BM(3+j*3+incr_b)*&
                        dble(sigma(nbr))

!if (iter.eq.25.and.incr.eq.3*98990) then                   
!   print*,"Sigma Nbrs: ",nbr,j,sigma(nbr)
!   print*,0.25*dble(sum(nodes(elements(nbr,1:4),1))),0.25*dble(sum(nodes(elements(nbr,1:4),2))), &
!        0.25*dble(sum(nodes(elements(nbr,1:4),3)))
!endif

                end do

!if (iter.eq.25.and.incr.eq.3*98990) then                   
!   print*,"tau_sigma",tnod(1+incr),tnod(2+incr),tnod(3+incr)
!endif

                incr=incr+3
                incr_b=incr_b+15
                
             end if
          end do

       case(2)
          do i=1,nelem
             if(use4cgrad(i)) then
                tnod(1+incr)=BM(1+incr_b)*dble(velocity(i))
                tnod(2+incr)=BM(2+incr_b)*dble(velocity(i))
                tnod(3+incr)=BM(3+incr_b)*dble(velocity(i))

!if (iter.eq.25.and.incr.eq.3*98990) then                
!   print*,"Velocity Center:",velocity(i)
!   print*,"i:",0.25*dble(sum(nodes(elements(i,1:4),1))),0.25*dble(sum(nodes(elements(i,1:4),2))), &
!        0.25*dble(sum(nodes(elements(i,1:4),3)))                
!endif

                do j=1,4
                   nbr=neighbors(i,j)
                   tnod(1+incr)=tnod(1+incr)+BM(1+j*3+incr_b)*&
                        dble(velocity(nbr))
                   tnod(2+incr)=tnod(2+incr)+BM(2+j*3+incr_b)*&
                        dble(velocity(nbr))
                   tnod(3+incr)=tnod(3+incr)+BM(3+j*3+incr_b)*&
                        dble(velocity(nbr))

!if (iter.eq.25.and.incr.eq.3*98990) then                   
!   print*,"Velocity Nbrs: ",nbr,j,velocity(nbr)
!   print*,0.25*dble(sum(nodes(elements(nbr,1:4),1))),0.25*dble(sum(nodes(elements(nbr,1:4),2))), &
!       0.25*dble(sum(nodes(elements(nbr,1:4),3)))
!endif

                end do

!if (iter.eq.25.and.incr.eq.3*98990) then                
!   print*,"tau_velocity",tnod(1+incr),tnod(2+incr),tnod(3+incr)
!endif

                incr=incr+3
                incr_b=incr_b+15
                
             end if
          end do

       end select
    end do
  end subroutine build_V
  !_____________________________________________________________________

  
  !_____________________________________________________________________
  subroutine setup_crossgrad
    !count the number of non-zeros in the cross gradient constraint matrices
    !structural metric 12 indicates cross gradient constraints
    implicit none
    integer :: i,j,k,ii,nbr,iphys,jphys,incr,incr_b,ind1,ind2
    logical :: found

    if(im_fmm) then
       write(*,*) ' FMM: MASTER SETTING UP CROSS-GRADIENT CONSTRAINTS'
    else
       write(*,*) ' E4D: MASTER SETTING UP CROSS-GRADIENT CONSTRAINTS'
    end if

    !!modify nphys and master_ranks if needed
    nphys = 2
    allocate(master_ranks(nphys))
    master_ranks(1) = 0
    master_ranks(2) = master_proc_fmm

    !!block4cgrad holds the block row index of BM
    !!block4cgrad(1:nphys-1,1) current phys_id
    !!block4cgrad(1:nphys-1,2) coupled phys_id
    !!block4cgrad(1:nphys-1,3) block row index
    allocate(block4cgrad(nphys-1,3))
    incr=0
    incr_b=0
    do iphys=1,nphys-1
       do jphys=iphys+1,nphys
          incr=incr+1
          if(master_ranks(iphys).eq.my_rank) then
             incr_b=incr_b+1
             block4cgrad(incr_b,1)=iphys
             block4cgrad(incr_b,2)=jphys
             block4cgrad(incr_b,3)=incr
          elseif(master_ranks(jphys).eq. my_rank) then
             incr_b=incr_b+1
             block4cgrad(incr_b,1)=jphys
             block4cgrad(incr_b,2)=iphys
             block4cgrad(incr_b,3)=incr
          end if
       end do
       if(master_ranks(iphys).eq.my_rank) then
          exit
       end if
    end do

    !if(im_fmm) then
    !   write(*,*) 'FMM: ',block4cgrad
    !else
    !   write(*,*) 'E4D: ',block4cgrad
    !end if

    !!use4cgrad determines which elements will have cross gradient
    !!constraints
    allocate(use4cgrad(nelem))
    use4cgrad = .false.

    !!ncgrad is the number of elements that will have cross gradient
    !!constraints
    ncgrad = 0

    !!nnz_cgrad is the number of nonzeros in the gradient operator WA
    nnz_cgrad = 0
!print*,zone_links        
    do i=1,nelem
       do j=1,nrz
          !!12 and 13 are joint inversion constraints
          if(smetric(j,2).eq.12 .or. smetric(j,2).eq.13) then
             if(smetric(j,1)==zones(i)) then
                do k = 1,4
                   nbr = neighbors(i,k)
!if (nbr<1) print*,i,j,nbr
                   if(nbr<1) goto 10 ! skip element which has boundary neighbor
                   if(zones(nbr) .ne. zones(i)) then  ! neighbor element is in different zone
                      found = .false.
                      do ii=2,zone_links(j,1)+1
                         if(zone_links(j,ii) .eq. zones(nbr)) then
                            found=.true.
                            exit
                         end if
                      end do
                      if(.not. found) goto 10
                   end if
                end do
                !if we made it to hear we add 15 to ccount
                !5 non-zero elements for each dimension
                nnz_cgrad = nnz_cgrad+15
                use4cgrad(i) = .true.
                ncgrad = ncgrad+1
             end if
          end if
10        continue
       end do
    end do

    !!WA holds the flattened 3D gradient operator for each
    !!element that has cross gradient constraints
    !!the third dimension references the element, the second
    !!dimension reference the neighbors, and the first 
    !!dimension reference the principle direction (x,y,or z)
    allocate(WA(nnz_cgrad))
    ncgrad=0
    do i=1,nelem
       if(use4cgrad(i)) then
          ncgrad = ncgrad+1
          ind1 = 1+(ncgrad-1)*15
          ind2 = 15+(ncgrad-1)*15
          call compute_WA(i,WA(ind1:ind2))
!do k=ind1,ind2
!   if (WA(k).eq.0) print*,i,ind1,ind2,WA(k)
!enddo
       end if
    end do

    !!wts4cg holds the relative weight for each element
    !!that has cross gradient constraints
    allocate(wts4cg(ncgrad))
    ncgrad=0
    do i=1,nelem
       if(use4cgrad(i)) then
          ncgrad=ncgrad+1
          do j=1,nrz
             if(smetric(j,1)==zones(i)) then
                wts4cg(ncgrad)=dble(zwts(j,1))
                !!if(im_fmm) then
                !!   write(*,*) zones(i),wts4cg(ncgrad)
                !!end if
                exit
             end if
          end do
       end if
    end do

    allocate(BM(15*ncgrad*(nphys-1)))
    allocate(tnod(3*ncgrad*(nphys-1)))
  
  end subroutine setup_crossgrad
  !_____________________________________________________________________

  !_____________________________________________________________________
  subroutine build_evol
    !build the vector of element volumes
    implicit none
    integer :: i,j,k
    real*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
    real*8, dimension(4,4) :: a
    
    if(allocated(evol)) deallocate(evol)
    allocate(evol(nelem))

    do k=1,nelem
       x1 = dble(nodes(elements(k,1),1))
       y1 = dble(nodes(elements(k,1),2))
       z1 = dble(nodes(elements(k,1),3))
       
       x2 = dble(nodes(elements(k,2),1))
       y2 = dble(nodes(elements(k,2),2))
       z2 = dble(nodes(elements(k,2),3))
       
       x3 = dble(nodes(elements(k,3),1))
       y3 = dble(nodes(elements(k,3),2))
       z3 = dble(nodes(elements(k,3),3)) 
       
       x4 = dble(nodes(elements(k,4),1))
       y4 = dble(nodes(elements(k,4),2))
       z4 = dble(nodes(elements(k,4),3))

       a(1:4,1) = (/ x1, x2, x3, x4 /)
       a(1:4,2) = (/ y1, y2, y3, y4 /)
       a(1:4,3) = (/ z1, z2, z3, z4 /)
       a(1:4,4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00 /)
      
       evol(k) = abs ( det4(a) ) / 6.0E+00
    end do
    
  end subroutine build_evol
  !_____________________________________________________________________
  !___________________________________________________________________________
  function det4 ( a1 )

    !computes the determinant of a 4x4 matrix
    implicit none    
    real*8 :: a1(4,4)
    real*8 :: a(4,4)
    real*8 :: det4
    integer ::i
    real*8 :: c1,c2,c3,c4
    
    a = dble(a1)
    c1=a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) 
    
    c2=a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))
    
    c3=a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    c4=a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    det4 = real(c1 - c2 + c3 - c4)
    
    return
  end function det4
  !__________________________________________________________________________
  

  subroutine build_cgwt(elm,rw)
    implicit none
    integer, intent(in) :: elm, rw
    integer :: comp,ti,i,j,nbr
    real*8, dimension(4) :: A
    real*8, dimension(4,3) :: N
    real*8, dimension(3) :: u,v,x
    real*8 :: vol
    real :: avesig,C1,C2
    logical :: hcheck
    real, dimension(5,4) :: AP
    real*8, dimension(4,4) :: ATA, XX
    real*8, dimension(4,5) :: G
    integer, dimension(4) :: INDX
    real, dimension(4) :: etvs,ftvs

    !if this node is on the boundary don't use cross gradients
    do i=1,4   
       if(neighbors(elm,i).le.0) then
          cg_wts(rw,:) = 0
          return
       end if
    end do

    !determine which component we're computing
    ti=mod(rw,3)
    if(ti.eq.0) then
       comp = 3
    else
       comp = ti
    end if

    !update CV,CK, and CJ if we're starting a new element
    if(cur_el .ne. elm) then   

       !get the midpoints
       AP(:,4) = 1
       do j=1,3
          AP(1,j) = 0.25*sum(nodes(elements(elm,:),j))
       end do
       do i=1,4
          nbr = neighbors(elm,i)
          do j=1,3
             AP(i+1,j) = 0.25*sum(nodes(elements(nbr,:),j))
          end do
 
       end do
      
       ATA = matmul(transpose(AP),AP)
     
       call MIGS(ATA,4,XX,INDX)
       G=matmul(XX,transpose(AP))
      
 
       cur_el = elm
       

    end if


    !Compute weights for the x-component of the cross gradient
   if(comp.eq.1) then
       C1 = 0
       C2 = 0
       if(im_fmm) then
          C1 = G(3,1)*log(sigma(elm))
          C2 = G(2,1)*log(sigma(elm))
          do i=2,5
             nbr = neighbors(elm,i-1)
             C1 = C1+G(3,i)*log(sigma(nbr))
             C2 = C2+G(2,i)*log(sigma(nbr))
          end do

       else
          C1 = G(3,1)*sqrt(velocity(elm))
          C2 = G(2,1)*sqrt(velocity(elm))
          do i=2,5
             nbr = neighbors(elm,i-1)
             C1 = C1+G(3,i)*sqrt(velocity(nbr))
             C2 = C2+G(2,i)*sqrt(velocity(nbr))
          end do  
       end if

       do i=1,5
          cg_wts(rw,i) = C1*G(2,i)-C2*G(3,i)
       end do
       if(im_fmm) cg_wts(rw,:) = -cg_wts(rw,:)

    end if


    !Compute weights for the y-component of the cross gradient
    if(comp.eq.2) then
       C1 = 0
       C2 = 0
       if(im_fmm) then
          C1 = G(1,1)*log(sigma(elm))
          C2 = G(3,1)*log(sigma(elm))
          do i=2,5
             nbr = neighbors(elm,i-1)
             C1 = C1+G(1,i)*log(sigma(nbr))
             C2 = C2+G(3,i)*log(sigma(nbr))
          end do

       else
          C1 = G(1,1)*sqrt(velocity(elm))
          C2 = G(3,1)*sqrt(velocity(elm))
          do i=2,5
             nbr = neighbors(elm,i-1)
             C1 = C1+G(1,i)*sqrt(velocity(nbr))
             C2 = C2+G(3,i)*sqrt(velocity(nbr))
          end do  
       end if

       do i=1,5
          cg_wts(rw,i) = C1*G(3,i)-C2*G(1,i)
       end do
       if(im_fmm) cg_wts(rw,:) = -cg_wts(rw,:)

    end if
          

    !Compute weights for the z-component of the cross gradient
    if(comp.eq.3) then
       C1 = 0
       C2 = 0
       if(im_fmm) then
          C1 = G(2,1)*log(sigma(elm))
          C2 = G(1,1)*log(sigma(elm))
          do i=2,5
             nbr = neighbors(elm,i-1)
             C1 = C1+G(2,i)*log(sigma(nbr))
             C2 = C2+G(1,i)*log(sigma(nbr))
          end do

       else
          C1 = G(2,1)*sqrt(velocity(elm))
          C2 = G(1,1)*sqrt(velocity(elm))
          do i=2,5
             nbr = neighbors(elm,i-1)
             C1 = C1+G(2,i)*sqrt(velocity(nbr))
             C2 = C2+G(1,i)*sqrt(velocity(nbr))
          end do  
       end if

       do i=1,5
          cg_wts(rw,i) = C1*G(1,i)-C2*G(2,i)
       end do
       if(im_fmm) cg_wts(rw,:) = -cg_wts(rw,:)

    end if
    return

  end subroutine build_cgwt
    !_____________________________________________________________________

  !_____________________________________________________________________
  subroutine get_areas_normal_vol(el,A,N,vol)
    implicit none
    integer, intent(in) :: el
    real*8, dimension(4), intent(out) :: A
    real*8, dimension(4,3), intent(out) :: N
    real*8, intent(out) :: vol
    real*8, dimension(3) :: u,v,x,elmid,pmid
    real*8 :: vlen
    integer :: i
    
    do i=1,3
       elmid(i) = .25*sum(nodes(elements(el,:),i))
    end do

    !face 1,2,3
    u = nodes(elements(el,2),:)-nodes(elements(el,1),:)
    v = nodes(elements(el,3),:)-nodes(elements(el,1),:)
    x=cross_prod(u,v)
    vlen = sqrt(dot_product(x,x))
    A(1)=.5*vlen
    N(1,:) = x/vlen
    
    !make sure the normal is pointing outward
    do i=1,3
       pmid(i) = (nodes(elements(el,1),i)+nodes(elements(el,2),i) + &
            nodes(elements(el,3),i))/3
    end do
    u=pmid-elmid
    if(dot_product(u,N(1,:))<0) N(1,:)=-N(1,:)
    
    !get the volume while we're here
    x = nodes(elements(el,4),:)-nodes(elements(el,1),:)
    !vol = sqrt(dot_product(u,cross_prod(v,x)))/6
    vol = (dot_product(u,cross_prod(v,x)))/6
    
    !face 1,2,4
    u = nodes(elements(el,2),:)-nodes(elements(el,1),:)
    v = nodes(elements(el,4),:)-nodes(elements(el,1),:)
    x=cross_prod(u,v)
    vlen = sqrt(dot_product(x,x))
    A(2)=.5*vlen
    N(2,:) = x/vlen
    
    !make sure the normal is pointing outward
    do i=1,3
       pmid(i) = (nodes(elements(el,1),i)+nodes(elements(el,2),i) + &
            nodes(elements(el,4),i))/3
    end do
    u=pmid-elmid
    if(dot_product(u,N(2,:))<0) N(2,:)=-N(2,:)
    
    !face 1,3,4
    u = nodes(elements(el,3),:)-nodes(elements(el,1),:)
    v = nodes(elements(el,4),:)-nodes(elements(el,1),:)
    x=cross_prod(u,v)
    vlen = sqrt(dot_product(x,x))
    A(3)=.5*vlen
    N(3,:) = x/vlen
    
    !make sure the normal is pointing outward
    do i=1,3
       pmid(i) = (nodes(elements(el,1),i)+nodes(elements(el,3),i) + &
            nodes(elements(el,4),i))/3
    end do
    u=pmid-elmid
    if(dot_product(u,N(3,:))<0) N(3,:)=-N(3,:)
    
    !face 2,3,4
    u = nodes(elements(el,3),:)-nodes(elements(el,2),:)
    v = nodes(elements(el,4),:)-nodes(elements(el,2),:)
    x=cross_prod(u,v)
    vlen = sqrt(dot_product(x,x))
    A(4)=.5*vlen
    N(4,:) = x/vlen
    
    !make sure the normal is pointing outward
    do i=1,3
       pmid(i) = (nodes(elements(el,2),i)+nodes(elements(el,3),i) + &
            nodes(elements(el,4),i))/3
    end do
    u=pmid-elmid
    if(dot_product(u,N(4,:))<0) N(4,:)=-N(4,:)
 

   
  end subroutine get_areas_normal_vol
  !_____________________________________________________________________
  
  !_____________________________________________________________________
  subroutine comp_Wmw
  implicit none
  !computes the area of the face between neighbor tets and assigns to Wmw
  integer :: i,j,k,icnt
  integer, dimension(3) :: n_indexes
  real, dimension(3) :: u,v,cp
  
  do i = 1,ccount
      !compute only if this is a spatial constraing (i.e. has a neighbor specified)
      if(rblock(i,2) .ne. 0) then
      	  icnt=0
          do j=1,4
              do k=1,4
                 if(elements(rblock(i,1),j)==elements(rblock(i,2),k)) then
                     icnt = icnt+1
                     n_indexes(icnt) = elements(rblock(i,1),j)
                     exit
                 end if
              end do
              if(icnt==3) then
                  exit
              end if
          end do
          if(icnt .ge. 3) then
              u = nodes(n_indexes(2),:) - nodes(n_indexes(1),:)
              v = nodes(n_indexes(3),:) - nodes(n_indexes(1),:)
              cp = 0.5*cross_prod(dble(u),dble(v))
              Wmw(i) = 0.5*sqrt(cp(1)**2 +cp(2)**2 + cp(3)**2)
          end if
      end if
  end do
  !adjust Wmw so that the average value is 1.0
  Wmw = Wmw*ccount/sum(Wmw)
  
  end subroutine comp_Wmw
  !_____________________________________________________________________
  !_____________________________________________________________________
  function cross_prod(u,v)
    real*8, dimension(3) :: cross_prod,u,v
    cross_prod(1) = u(2)*v(3) - u(3)*v(2)
    cross_prod(2) = u(3)*v(1) - u(1)*v(3)
    cross_prod(3) = u(1)*v(2) - u(2)*v(1)
  end function cross_prod
  !______________________________________________________________________

	






























  !_____________________________________________________________________
  subroutine build_Wm

    implicit none
    integer :: n_homo,i,j,k,nzmax,n_homo_elem,count,row,col,i1,i2,lzone
    integer :: znn,nbr,nbr2,x1,x2,x3,z1,z2,z3
    integer, dimension(:), allocatable :: homo_par
    integer, dimension(:,:), allocatable :: neighbors          
    real :: midxi,midyi,midzi,midxn,midyn,midzn,rx,ry,rz,r,x,rho,eps

    
    if(gs_flag) then
       if(allocated(Wm)) then 
          return
       end if
       allocate(Wm(nelem),wrows(nelem),wcols(nelem))
       ccount=nelem
       do i=1,nelem
          wrows(i)=i
          wcols(i)=i
          Wm(i)=1
       end do
       return
    end if


    !!If Wm is already constructed then there is nothing to do
    if(allocated(Wm)) then
       !return   
       deallocate(par_map)
       deallocate(sigma_par,rsigma_par)
       deallocate(wrows,wcols,Wm)       
    end if
    
    !!check to see if there are heterogeneous zones. if so, read in the neighbors
    if(sum(reg_opt)>0) then          
       allocate(neighbors(nelem,4))
       do i=1,80
          if(mshfile(i:i)=='.') then
             open(21,file=mshfile(1:i+1)//".neigh",status='old',action='read');
             read(21,*) 
             do j=1,nelem
                read(21,*) k,neighbors(j,1:4)
             end do
             goto 25
          end if
       end do
    end if
25  continue
    close(21)
  
    
    
    !!send a warning if there is a discrepancy be the number of zones 
    !!and the number of zones which are regularized
    nzmax=int(maxval(zones))

    if(nzmax .ne. nrz) then
       if(.not.gs_flag) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) '!!!WARNING!!! THERE ARE ',nzmax,' ZONES IN THE ELEMENTS FILE AND ',nrz,' ZONES REGULARIZED'
          close(51)
       else
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) '!!!WARNING!!! THERE AR',nzmax,' ZONES AND ',nrz, ' ARE REGULARIZED'
          write(51,*) 'MAKES SURE THAT ZONES ',nrz+1,' THROUGH ',nzmax,' HAVE SEMIVARIOGRAM CONSTRAINTS'
          write(51,*) 'THIS IS NOT CHECKED INTERNALLY'
          close(51)
       end if

    end if
    
    

    !!determine how many parameters we will estimate 
    n_homo = 0
    allocate(homo_par(nrz))
    do i=1,nrz
       if(reg_opt(i)==0 ) then
          n_homo=n_homo+1
          homo_par(n_homo)=i
       end if
    end do
    
    !!count the number of elements in homogeneous zones
    n_homo_elem=0
    do i=1,nelem
       if(reg_opt(zones(i))==0 ) then
          n_homo_elem = n_homo_elem+1
       end if
    end do
   
    !!allocate the new parameter vector
    npar=nelem-n_homo_elem+n_homo
    allocate(sigma_par(npar),rsigma_par(npar))
  
    
    !!build the vector mapping the conductivities to the inversion parameters
    allocate(par_map(nelem))
    count=0
    do i=1,nelem
       if(reg_opt(zones(i))==0 ) then
          !!this is an element in a homogeneous zone ... find which zone
          do j=1,n_homo
             if(homo_par(j)==zones(i)) then
                par_map(i)=j
                if(invi) then
                   sigma_par(j)=sigmai(i)
                else
                   sigma_par(j)=sigma(i)
                end if
             end if
          end do
       else
          count=count+1
          par_map(i) = n_homo+count
          if(invi) then
             sigma_par(par_map(i)) = sigmai(i)
             !rsigma_par(par_map(i))=refsigi(i)
          else
             sigma_par(par_map(i)) = sigma(i)
             rsigma_par(par_map(i))= refsig(i)
          end if
       end if
    
    end do
   
    !!now map sigma_par back to sigma to make sure they're consistent
    if(invi) then
       do i=1,nelem
          sigmai(i) = sigma_par(par_map(i))
       end do
    else
       do i=1,nelem
          sigma(i)=sigma_par(par_map(i))
          refsig(i)=rsigma_par(par_map(i))
       end do
    end if
    
    !!Now we're ready to build the sparse regularization matrix between inversion
    !!parameters. We start with the homogeneous zones and then move to the
    !!heterogeneous zones. Since we don't know how many elements are in the matrix
    !!apriori, we'll do a false run to count, allocate the matrix, and then do the 
    !!computations to fill in the matrix.
    
    ccount = 0
    !!Homogeneous zone regularization count
    do i=1,nrz
       !!check to see if this is a homogeneous zone and is linked to other zones
       if(reg_opt(i)==0  .and. zone_links(i,1)>0) then
          do j=1,zone_links(i,1)
             if(reg_opt(zone_links(i,j+1))==0) ccount=ccount+2
          end do
       end if      
    end do
 
    
    !!smallest model regularization count
    do i=1,nelem
       if(reg_opt(zones(i)) .ge. 1 .and.  reg_opt(zones(i)) .le. 4 ) then
          !!this element belongs to a zone with smallest model smoothing
          do j=1,4
             nbr=neighbors(i,j)                           !!nbr = element index of the neighbor
             if(nbr<=0) goto 102              		  !!there is no neighbor on face j so skip this neighbor
                                                          !!or we have already covered this pair in nbr<i.
             znn = zones(nbr)                             !!zone of the neighbor          
             if(zones(i)==znn) goto 101                   !!elements i and nbr are in the same zone
             
             do k=1,zone_links(zones(i),1)
                if(zone_links(zones(i),k+1)==znn) goto 101          !!elements i and nbr are linked neighbors
             end do
             goto 102
             
101          continue
             ccount=ccount+2
             
102          continue
          end do
       end if
    end do

    
    !!if there are homogeneous zone, external constraints are not allowed
    !!check to make sure we don't have both and allocate 
    if(allocated(ex_vals) .and. n_homo>0) then
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) "Homogeneous zones are specified. External constraints will be ignored"
       close(51)
       nec=0
    end if
    allocate(Wm(ccount),wrows(ccount),wcols(ccount))
  
       
    
    !!build the regularization matrix
    ccount=0
    row=0
    col=0
    !!build the homogeneous zone regularization
    do i=1,nrz
       !!check to see if this is a homogeneous zone and is linked to another zones
       if(reg_opt(i)==0 .and. zone_links(i,1)>0) then
          do j=1,zone_links(i,1)
             
             !!find which parameters hold the linked zones
             lzone = zone_links(i,j+1)
             if(reg_opt(lzone)==0) then
                
                do k=1,n_homo
                   if(homo_par(k)==lzone) i2=k 
                   if(homo_par(k)==i) i1 = k
                end do
                
                row = row+1; col = i1; ccount=ccount+1
                wrows(ccount)=row; wcols(ccount)=col; Wm(ccount)=zwts(i,1)
                
                col = i2; ccount=ccount+1
                wrows(ccount)=row; wcols(ccount)=col; Wm(ccount)=-zwts(i,1)
               
             end if
          end do
       end if     
    end do
    
    
    do i=1,nelem
       !!smallest model regularization
       if(reg_opt(zones(i)) .ge. 1 .and. reg_opt(zones(i)) .le. 4 ) then
          
          
          do j=1,4
             nbr=neighbors(i,j)                        !!nbr = element index of the neighbor
             if(nbr<=0) goto 111            !!there is no neighbor on face j so skip this neighbor
                                                       !!or we have already covered this pair
             znn = zones(nbr)                          !!zone of the neighbor          
             if(zones(i)==znn) goto 110                !!elements i and nbr are in the same zone
             
             do k=1,zone_links(zones(i),1)
                if(zone_links(zones(i),k+1)==znn)  goto 110   !!elements i and nbr are linked neighbors
             end do
             goto 111
             
110          continue
          
             if(reg_opt(zones(i))==1 .or. reg_opt(zones(i))==3 ) then              !!smallest model smoothing
                row=row+1; col=par_map(i); ccount=ccount+1
                wrows(ccount)=row;  wcols(ccount)=col;  Wm(ccount)=zwts(zones(i),1)
                
                col=par_map(nbr); ccount=ccount+1
                wrows(ccount)=row;  wcols(ccount)=col;  Wm(ccount)=-zwts(zones(i),1) 
             end if
             
            
             if(reg_opt(zones(i))==4) then
                eps=zwts(zones(i),2)
                x=log(sigma_par(par_map(i)))*zwts(zones(i),1) - log(sigma_par(par_map(nbr)))*zwts(zones(i),1)
                rho = (eps**2)/(x**2 + eps**2)
                !write(*,*) eps,x,rho
                row=row+1; col=par_map(i); ccount=ccount+1
                wrows(ccount)=row;  wcols(ccount)=col;  Wm(ccount)=rho*zwts(zones(i),1)
                !write(*,*) wrows(ccount),wcols(ccount), Wm(ccount)
                col=par_map(nbr); ccount=ccount+1
                wrows(ccount)=row;  wcols(ccount)=col;  Wm(ccount)=-rho*zwts(zones(i),1) 
                !write(*,*) wrows(ccount),wcols(ccount), Wm(ccount)
                
             end if

             if(reg_opt(zones(i))==2) then              !!first order weighted smoothing

                midxi=0.25*sum(nodes(elements(i,1:4),1))
                midyi=0.25*sum(nodes(elements(i,1:4),2))
                midzi=0.25*sum(nodes(elements(i,1:4),3))  
                midxn=0.25*sum(nodes(elements(nbr,1:4),1))
                midyn=0.25*sum(nodes(elements(nbr,1:4),2))
                midzn=0.25*sum(nodes(elements(nbr,1:4),3))
                rx=(midxi-midxn)**2
                ry=(midyi-midyn)**2
                rz=(midzi-midzn)**2
                r=rx+ry+rz

               
                row=row+1; col=par_map(i); ccount=ccount+1
                wrows(ccount)=row;  wcols(ccount)=col;  
                Wm(ccount)= zwts(zones(i),2)*rx/r + zwts(zones(i),3)*ry/r + zwts(zones(i),4)*rz/r
                
                col=par_map(nbr); ccount=ccount+1
                wrows(ccount)=row;  wcols(ccount)=col;  Wm(ccount)=-Wm(ccount-1)!zwts(zones(i),1) 

                
                !row = row+1; col=par_map(i); ccount=ccount+1
                !wrows(ccount)=row; wcols(ccount)=col;
                !Wm(ccount) =1! zwts(zones(i),2)*rx/r + zwts(zones(i),3)*ry/r + zwts(zones(i),4)*rz/r
                
                !col=par_map(nbr); ccount=count+1
                !wrows(ccount)=row; wcols(ccount)=col;
                !Wm(ccount)=-1!Wm(ccount-1)
                
             end if

111          continue
          end do

       end if
    end do

    return

1001 continue
    
    allocate(wrows(nelem),wcols(nelem),par_map(nelem))
    allocate(Wm(nelem),sigma_par(nelem),rsigma_par(nelem))
    ccount=nelem
    npar=nelem
    do i=1,nelem
       wrows(i)=i
       wcols(i)=i
       Wm(i)=1
       par_map(i)=i
       if(invi) then
          sigma_par(i)=sigmai(i)
       else
          sigma_par(i)=sigma(i)
       end if
    end do
 
   !!now map sigma_par back to sigma to make sure they're consistent
    if(invi) then
       do i=1,nelem
          sigmai(i) = sigma_par(par_map(i))
       end do
    else
       do i=1,nelem
          sigma(i) = sigma_par(par_map(i))
       end do
    end if
 

  end subroutine build_Wm
  !_____________________________________________________________________

  !_____________________________________________________________________
  subroutine compute_WA(el,grd)
    implicit none
    integer, intent(in) :: el
    real*8, dimension(3,5), intent(out) :: grd

    real*8, dimension(3,3) :: atild
    real*8, dimension(3,3) :: atildi
    real*8, dimension(4,3) :: dist
    real*8, dimension(3)   :: elmid
    real*8, dimension(3,5) :: temp_a
    real*8, dimension(4,5) :: temp_b
    integer, dimension(4,3):: iface
    integer, dimension(4)  :: indx
    integer :: j,k,nbr

    elmid(1)=0.25*dble(sum(nodes(elements(el,1:4),1)))
    elmid(2)=0.25*dble(sum(nodes(elements(el,1:4),2)))
    elmid(3)=0.25*dble(sum(nodes(elements(el,1:4),3)))

    ! Get C-matrix
    do j=1,4
       nbr = neighbors(el,j)
       do k=1,3
          dist(j,k) = 0.25*dble(sum(nodes(elements(nbr,1:4),k))) - elmid(k)
          !print*,"Nbr",j,"elmid",k,0.25*dble(sum(nodes(elements(nbr,1:4),k)))
          !print*,j,k,dist(j,k)
       end do
    end do

    !build the shape function -> C^TC
    atild = 0
    do j=1,3
       do k=1,3
          atild(j,k)=dot_product(dist(1:4,j),dist(1:4,k))
          !print*,j,k,atild(j,k)
       end do
    end do
    
    !!invert atild to get atildi
    !!note MIGS routine is located in mat_inv.f
    call MIGS(atild,3,atildi,indx)
    !print*,atildi

    ! L-matrix
    temp_b(1,:) = (/-1,1,0,0,0/)
    temp_b(2,:) = (/-1,0,1,0,0/)
    temp_b(3,:) = (/-1,0,0,1,0/)
    temp_b(4,:) = (/-1,0,0,0,1/)
    do j=1,5
       do k=1,3
          ! C^T*L matrix
          temp_a(k,j) = dot_product(dist(1:4,k),temp_b(1:4,j))
       end do
    end do

    do j=1,5
       do k=1,3
          grd(k,j) = dot_product(atildi(k,1:3),temp_a(1:3,j))
          !print*,k,j,grd(k,j)
       end do
    end do

  end subroutine compute_WA
  !_____________________________________________________________________
  
end module mod_con
