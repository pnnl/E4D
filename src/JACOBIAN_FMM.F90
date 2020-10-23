Module jaco_fmm
  
  use vars
  use fmm_vars
  use build_amap
  use input_fmm
  
  logical, dimension(:), allocatable :: eligible
  
contains
  
  !_____________________________________________________________________________
  subroutine build_jaco_fresnel
    implicit none
    integer :: my_srank_fmm, n_srank_fmm,ierr,src,myindx
    integer :: i,j,mn,a,m,jrow,crank,el
    real :: tsr,tel,wisum
    real*8 :: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,evol
    real, dimension(nnodes) :: tsum
    real, dimension(:,:), allocatable :: att
    logical, dimension(nnodes) :: fzon
    real, dimension(nelem) :: wi

    !allocate Jaco
    nj_rows=jind(my_rank_fmm,2)-jind(my_rank_fmm,1)+1
    if(.not. allocated(Jaco)) then
       allocate(Jaco(nj_rows,nelem))
    end if
    Jaco = 0  

    allocate(att(nnodes,tns))
    do i=1,tns
       do crank = 1,n_rank_fmm-1
          if(i.ge.sind(crank,1) .and. i.le.sind(crank,2)) then
             if(my_rank_fmm .eq. crank) then
                att(:,i) = ttimes(:,i-sind(crank,1)+1)
             end if
             call MPI_BCAST(att(:,i),nnodes,MPI_REAL,crank-1,SCOMM_FMM,ierr)
             exit
          end if
       end do
    end do
    
    jrow = 0
    do mn=jind(my_rank_fmm,1),jind(my_rank_fmm,2)
       fzon = .false.
       a = s_conf_fmm(mn,1)
       m = s_conf_fmm(mn,2)
       tsr = att(s_nods(m),a)
       tsum = att(:,a)+att(:,m)
       fzon = (tsum - tsr) .le. .5/frq(a)
       wi = 0
       do el=1,nelem
          i=0
          do j=1,4
             if(fzon(elements(el,j))) i=i+1
          end do
          if(i.lt.3) goto 11

          tel = 0
          do j=1,4
             tel = tel+tsum(elements(el,j))
          end do
          tel = .25*tel

          x1 = dble(nodes(elements(el,1),1))
          y1 = dble(nodes(elements(el,1),2))
          z1 = dble(nodes(elements(el,1),3))
          
          x2 = dble(nodes(elements(el,2),1))
          y2 = dble(nodes(elements(el,2),2))
          z2 = dble(nodes(elements(el,2),3))
          
          x3 = dble(nodes(elements(el,3),1))
          y3 = dble(nodes(elements(el,3),2))
          z3 = dble(nodes(elements(el,3),3)) 
          
          x4 = dble(nodes(elements(el,4),1))
          y4 = dble(nodes(elements(el,4),2))
          z4 = dble(nodes(elements(el,4),3))
          
          call get_vol (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,evol)
          wi(el) = evol*sin(1-2*frq(a)*(tel-tsr))
          if(wi(el)<0) wi(el)=0
11           continue
       end do
       jrow=jrow+1
       !Jaco(jrow,el) = tsr*wi/sum(wi)
       if(sum(wi)>0) wi = wi/sum(wi)
       do el=1,nelem
          if(velocity(el).gt.0) then
             Jaco(jrow,el) = wi(el)*tsr/sqrt(velocity(el))
          end if
       end do

!!$       if(my_rank_fmm .eq.6) then
!!$          open(20,file='jaco_test.txt',status='replace',action='write')
!!$          write(20,*) nelem,' 1'
!!$          do i=1,nelem
!!$             write(20,*) Jaco(jrow,i)
!!$          end do
!!$          close(20)
!!$          return
!!$       end if
    end do
 
    deallocate(att)
    if(my_rank_fmm==1) then
       open(13,file='jaco_fmm.txt',status='replace',action='write')
       write(13,*) nelem,' 1'
       do i=1,nelem
          write(13,*) Jaco(1,i)
       end do
    end if
    close(13)
    !Jaco = 0
  end subroutine build_jaco_fresnel
  !_____________________________________________________________________________

  !_____________________________________________________________________________
  subroutine build_jaco_raytrace
    implicit none
    logical :: check_faces = .true.
    logical :: upcheck,fcfound
    logical :: near_source
    
    integer :: c_row,i,j,edge,fc,nd,rec_nod,rrow,fc_exit,obsn,old_el
    integer :: cur_el,my_src,src_nod,st_el,nbr,k,l,cnt,m,ofc,step,ierr,nfailed
    
    real*8, dimension(4,4) :: atild              !pre-shape linear shape function matrix 
    real*8, dimension(4,4) :: atildi             !shape function matrix
    real*8, dimension(4,4) :: temp_a             !temporary placeholder
    
    integer, dimension(4) :: indx,tindx,oindx    !computation space for MIGS routine
    integer, dimension(:,:),allocatable :: A_tmp !temporary accounting matrix
    integer, dimension(:), allocatable :: failed
    
    real*8, dimension(4,3) :: t_n
    real*8 :: evol,rt                            !nodal volume
    
    real*8 :: A_kij                                !forward matrix coefficient sum
    
    real*8, dimension(4) :: local_tt                 !travel time on four tet nodes
    real*8, dimension(3) :: pcur,pup,p1,p2,p3,p4,u,grad,ograd,opcur,avegrad
    real*8, dimension(3) :: v1,v2,v3,v4,Ur,Vr,Vs,pip,n
    real*8, dimension(3,3) :: Us
    integer, dimension(3) ::   edge_rot,n1,n2,n3
    integer, dimension(4,4) :: edge_seq
    real*8 :: d, sens, rad, radtest, rlim, ttcur, ttold
    character*40 :: str

    !allocate Jaco
    nj_rows=jind(my_rank_fmm,2)-jind(my_rank_fmm,1)+1
    if(.not. allocated(Jaco)) then
       allocate(Jaco(nj_rows,nelem))
    end if
    Jaco = 0
    
    if(.not. allocated(eligible)) allocate(eligible(nelem))
    allocate(failed(nelem))
    eligible = .true.

    edge_seq(1,1:4) = (/ 1,2,3,1 /)
    edge_seq(2,1:4) = (/ 2,3,4,2 /)
    edge_seq(3,1:4) = (/ 3,4,1,3 /)
    edge_seq(4,1:4) = (/ 4,1,2,4 /)
    
    !loop over the survey and do the ray tracing for sources I own
    c_row = 0
    step = 0
    nfailed = 0
   
    do obsn=1,nm_fmm
 
       if(obsn .ge. jind(my_rank_fmm,1) .and. obsn .le. jind(my_rank_fmm,2)) then
      
          my_src = s_conf_fmm(obsn,1) - sind(my_rank_fmm,1) + 1
          eligible = .true.
          
          !f(s_conf_fmm(i,1).ge.sind(my_rank_fmm,1) .and. s_conf_fmm(i,1).le.sind(my_rank_fmm,2)) then
          c_row = c_row + 1
          src_nod = s_nods(s_conf_fmm(obsn,1))
          
          !set the starting node
          rec_nod = r_nods(s_conf_fmm(obsn,2))
          pcur = nodes(rec_nod,:)
          upcheck = .false.
          v1=pcur - nodes(src_nod,:)
          rad=sqrt(dot_product(v1,v1))
          rlim = 0.2*rad
          ttold = ttimes(rec_nod,my_src)
          
          !loop over the elements that use this node and
          !find the one that contains this ray
          avegrad = 0
          do rrow = ring_map(rec_nod-1)+1,ring_map(rec_nod)
             fcfound = .false.
             cur_el = rings(rrow)
             !compute the tt gradient
             call compute_grad(cur_el, my_src,grad)
             grad=-grad
             avegrad = avegrad+grad
             
             !find the opposite face
             do fc=1,4
                do j=1,3
                   indx(j)=elements(cur_el,edge_seq(fc,j))
                end do
                if(indx(1).ne.rec_nod .and. indx(2).ne.rec_nod .and. indx(3).ne.rec_nod) exit
             end do
             
             !determine if the ray passes through the opposite face and compute the crossing point
             !call check_cross(indx(1:3),dble(nodes(rec_nod,:)),grad,fcfound,pcur,ierr)
             call check_cross(indx(1:3),pcur,grad,fcfound,pup,ierr)
             !fcfound is true then we've found the starting element and can exit
             if(fcfound) exit
             
          end do
          
          !if we haven't found the element then use the average gradient, which should guarantee
          !an exit face
          if(.not. fcfound) then
             avegrad = avegrad/real(rrow)
             do rrow = ring_map(rec_nod-1)+1,ring_map(rec_nod)
                fcfound = .false.
                cur_el = rings(rrow)
                
                !find the opposite face
                do fc=1,4
                   do j=1,3
                      indx(j)=elements(cur_el,edge_seq(fc,j))
                   end do
                   if(indx(1).ne.rec_nod .and. indx(2).ne.rec_nod .and. indx(3).ne.rec_nod) exit
                end do
                
                !determine if the ray passes through the opposite face and compute the crossing point
                call check_cross(indx(1:3),pcur,avegrad,fcfound,pup,ierr)
                
                !fcfound is true then we've found the starting element and can exit
                if(fcfound) exit
                
             end do
          end if
          
          if(.not.fcfound) then          
             !the starting element was not found
             write(*,*) my_rank_fmm,"Couldn't find the starting element for measurement",obsn
             nfailed=nfailed+1
             failed(nfailed) = obsn 
             goto 1000
          end if
         
          !raytrace back to the source
100       continue
          
          !make sure we're getting closer to the source
 
          call interp_tt(pcur,cur_el,my_src,ttcur)
 

          if(real(ttcur) .gt. real(ttold)) then
             !write(*,*) 'Time increasing at: ',ttold,ttcur
             !write(*,*) my_rank_fmm, 'Done with measurement: ',obsn
             goto 1000
          endif
          ttold=ttcur
 
          if(elements(cur_el,1).eq.src_nod .or. elements(cur_el,2).eq.src_nod .or. &
               elements(cur_el,3).eq.src_nod .or. elements(cur_el,4).eq.src_nod) then
             !if we're in the source element then exit the ray tracing loop

             v1=pcur-nodes(src_nod,:)
             Jaco(c_row,cur_el) = sqrt(dot_product(v1,v1))
         
             !write(*,*) my_rank," all done"
             goto 1000
          end if
 
          !loop over the neighbors to find the one that shares this face
          do i = 1,4
 
             nbr = nbrs(cur_el,i)
             if(nbr .lt. 1) cycle
 
              
             cnt=0
             do j=1,4
 
                do k=1,3
                   if(elements(nbr,j).eq.indx(k)) then
                      cnt=cnt+1
                      tindx(cnt)=elements(nbr,j)
                      exit
                   end if
                end do
               
        
                if(cnt ==3) goto 101
             end do
          end do
          
          write(*,*) my_rank_fmm," Couldn't find my neighbor"
          nfailed=nfailed+1
          failed(nfailed) = obsn 
          goto 1000
          
101       continue  
         

          old_el = cur_el
          cur_el = nbr
          oindx = indx
          indx = tindx
          ofc=fc
          opcur = pcur
          ograd = grad
          
          !find the face for this element
          do fc=1,4
             cnt = 0
             do j=1,3
                
        
                do k=1,3
                   if(elements(cur_el,edge_seq(fc,j)).eq.indx(k)) then
                      cnt=cnt+1
                      exit
                   end if
                end do
                
             
                if(cnt .eq.3) goto 102
             end do
          end do
          
          write(*,*) my_rank_fmm," Couldn't find the neighbor face"
          nfailed=nfailed+1
          failed(nfailed) = obsn 
          goto 1000
          
102       continue
          
         
          !compute the gradient for this element
        
       
          call compute_grad(cur_el, my_src,grad)
          
            
          grad=-grad
     
          !loop over the faces and see if the ray crosses this face 
          do i=1,4
             if(i.ne.fc) then
                
       
                do j=1,3
                   indx(j)=elements(cur_el,edge_seq(i,j))
                end do
              
  

           
                call check_cross(indx(1:3),pcur,grad,fcfound,pup,ierr)
  
          
                if(fcfound) then
                   
           
                   v1=pup-pcur
                   Jaco(c_row,cur_el) = sqrt(dot_product(v1,v1))
                   pcur = pup
                   eligible(cur_el) = .false.
                   
               
                   goto 100
                end if
                
             end if
          end do
          
          !this part is a contingency in case the ray bounces back
          !back into the previous element, based on the gradient in
          !in this element. If so, we use the gradient from the previous
          !element to compute the ray through this one. 
          grad = ograd
          do i=1,4
             if(i.ne.fc) then

           
                do j=1,3
                   indx(j)=elements(cur_el,edge_seq(i,j))
                end do
              
                    
                call check_cross(indx(1:3),pcur,grad,fcfound,pup,ierr)
          
                    
                if(fcfound) then
                  
                   v1=pup-pcur
                   Jaco(c_row,cur_el) = sqrt(dot_product(v1,v1))
                   pcur = pup
                   eligible(cur_el) = .false.
                   goto 100
                end if
                
             end if
          end do
          
          
          
          
          !write(*,*) my_rank_fmm, "couldn't find the exit for: ",obsn
          nfailed=nfailed+1
          failed(nfailed) = obsn 
          
         
1000      continue
         
       end if
    end do
  
    !report any failed traces
    call MPI_SEND(nfailed,1,MPI_INTEGER,0,0,FMM_COMM, ierr)
    call MPI_SEND(failed(1:nfailed),nfailed,MPI_INTEGER,0,my_rank_fmm,FMM_COMM, ierr)
    deallocate(failed)
    return
    !This part will print all of the ray traces in sigma file format
    !if the return statement above is commented out
    c_row = 0
    do obsn=1,nm_fmm
       
       if(obsn .ge. jind(my_rank_fmm,1) .and. obsn .le. jind(my_rank_fmm,2)) then
          my_src = s_conf_fmm(obsn,1) - sind(my_rank_fmm,1) + 1
          c_row = c_row + 1
          write(str,"(A4,I3.3,A1,I2.2,A1,I2.2,A4)") 'obs_',obsn,'_',s_conf_fmm(obsn,1),'_',s_conf_fmm(obsn,2),'.txt'
          open(23,file=trim(str), status='replace',action='write')
          write(23,*) nelem," 1"
          do i=1,nelem
             write(23,*) Jaco(c_row,i)
          end do
          close(23)
       end if
    end do
 
    
end subroutine build_jaco_raytrace
!_____________________________________________________________________________

!_____________________________________________________________________________
subroutine trep(spot)
  implicit none
  integer :: spot
  if(my_rank_fmm .eq. 1) write(21,*),spot
 if(my_rank_fmm .eq. 2) write(22,*),spot
if(my_rank_fmm .eq. 3) write(23,*),spot
if(my_rank_fmm .eq. 4) write(24,*),spot
if(my_rank_fmm .eq. 5) write(25,*),spot
if(my_rank_fmm .eq. 6) write(26,*),spot
if(my_rank_fmm .eq. 7) write(27,*),spot
end subroutine trep
!_____________________________________________________________________________
!_____________________________________________________________________________
subroutine interp_tt(pt,el,src,tt)
  implicit none
  real*8, dimension(3), intent(in) :: pt
  integer, intent(in) :: src,el
  real*8, intent(out) :: tt
  
  real*8, dimension(4,4) :: atild              !pre-shape linear shape function matrix 
  real*8, dimension(4,4) :: atildi             !shape function matrix
  real*8, dimension(4,4) :: temp_a             !temporary placeholder
  real*8, dimension(4) :: local_tt,lambda      !travel time on four tet nodes
  real*8, dimension(4,3) :: t_n
  integer, dimension(4) :: indx        
  integer :: j,k
  
  do j=1,4
     t_n(j,:) = dble(nodes(elements(el,j),:))
  end do
  
  !build the shape function
  atild(1,:) = 1.0
  do j=2,4
     do k=1,4
        atild(j,k)=t_n(k,j-1)
     end do
  end do
  temp_a=atild
    
  !!invert temp_a to get atildi
  !!note MIGS routine is located in mat_inv.f  
  call MIGS(temp_a,4,atildi,indx)
    
  !get the travel times for each node in this element
  do j=1,4   
     local_tt(j) = ttimes(elements(el,j),src)
     lambda(j) = atildi(j,1) + dot_product(atildi(j,2:4),pt)
  end do
    
  tt = dot_product(local_tt,lambda)
  
  
end subroutine interp_tt
!_____________________________________________________________________________

!_____________________________________________________________________________
subroutine check_cross(npts,rpos,grd,fnd,cp,code)
  implicit none
  integer, dimension(3), intent(in) :: npts
  real*8, dimension(3), intent(in) :: grd,rpos
  logical, intent(inout) :: fnd
  real*8, dimension(3), intent(out) :: cp
  integer, intent(out) :: code

  real*8, dimension(3) :: p0, r1, r2, r3, n
  real*8, dimension(3) :: u, v, w,vxw,vxu,uxw,uxv
  real*8 :: d,dp1,dp2,denom,r,t 
  integer :: i

!!$  if(my_rank_fmm ==1) then
!!$     do i=1,3
!!$        write(*,*) nodes(npts(i),:),1.5
!!$     end do
!!$     write(*,*) rpos
!!$  end if
  ! this is the wrong face by default
  fnd = .false.

  !find the point in space where the ray crosses the plane
  !find the normal to the plane
  r1 = dble(nodes(npts(2),:)-nodes(npts(1),:))
  r2 = dble(nodes(npts(3),:)-nodes(npts(2),:))
  cp = cross_prod(r1,r2)
  n = cp/dot_product(cp,cp)
  p0 = nodes(npts(1),:)
  
  denom = dot_product(grd,n)
  !if(my_rank_fmm .eq. 1) write(*,*) 
  !if(my_rank_fmm .eq.1) write(*,*) 'denom ',denom
  if(denom .eq. 0) then
     if(abs(denom)<1-6) then
        code = 11
     else
        code = 1
     end if
     !if(my_rank_fmm .eq. 1) write(*,*) 'denom',denom
     return
  end if

  d = dot_product((p0-rpos),n)/denom
 !if(my_rank_fmm .eq.1) write(*,*) 'd ',d
  if(d<0.0) then
     if(abs(d)<1e-6) then
        code = 22
     else
        code = 2
     end if
     !if(my_rank_fmm .eq. 1) write(*,*) 'd less zero',d
     return                !the crossing point is in the wrong direction
  end if
  !cp is now the point where the ray crosses the plane
  cp = rpos + d*grd

  !check to see if cp is within the triangle formed by this face
  w = cp - nodes(npts(1),:)
  u = nodes(npts(2),:) - nodes(npts(1),:)
  v = nodes(npts(3),:) - nodes(npts(1),:)
  vxw = cross_prod(v,w)
  vxu = cross_prod(v,u)
  uxw = cross_prod(u,w)
  uxv = cross_prod(u,v)

  dp1 = dot_product(vxw,vxu)
 !if(my_rank_fmm .eq.1) write(*,*) 'dp1 ',dp1
  if(dp1 < 0.0) then
     if(abs(dp1)<1e-6) then
        code = 33
     else
        code = 3
     end if
     !if(my_rank_fmm .eq. 1) write(*,*) 'dp1 less zero',dp1
     return
  end if
  
  dp2 = dot_product(uxw,uxv)
  !if(my_rank_fmm .eq.1) write(*,*) 'dp2 ',dp2
  if(dp2 < 0.0) then
     if(abs(dp2)<0) then
        code = 44
     else
        code = 4
     end if
     !if(my_rank_fmm .eq. 1) write(*,*) 'dp2 less zero',dp2
     return
  end if

  denom = sqrt(dot_product(vxu,vxu))
  r = sqrt(dot_product(vxw,vxw))/denom
  !  if(my_rank_fmm .eq.1) write(*,*) 'r ',r
  if(r.ge.1.0) then
     if(abs(r-1.0)<1e-6) then
        code = 55
     else
        code = 5
     end if
     !if(my_rank_fmm .eq. 1) write(*,*) 'r > 1',r
     return
  end if
  t = sqrt(dot_product(uxw,uxw))/denom
 !if(my_rank_fmm .eq.1) write(*,*) 't ',t
  if(t .ge. 1.0 .or. r+t .ge. 1.0) then
     if(abs(t-1)<1e-6 .or. abs(r+t-1)<1e-6) then
        code = 66
     else
        code = 6
     end if
     !if(my_rank_fmm .eq. 1) write(*,*) 't or r+t > 1',t,r
     return
  end if
  fnd = .true.
!!$  if(my_rank_fmm==1) then
!!$     do i=1,3
!!$        write(20,*) nodes(npts(i),:),2
!!$     end do
!!$     write(20,*) cp,3
!!$  end if
  return

end subroutine check_cross
!_____________________________________________________________________________
!_____________________________________________________________________________
subroutine compute_grad(el,src,grd)
  implicit none
  integer, intent(in) :: el,src
  real*8, dimension(3), intent(out) :: grd
  
  real*8, dimension(4,4) :: atild              !pre-shape linear shape function matrix 
  real*8, dimension(4,4) :: atildi             !shape function matrix
  real*8, dimension(4,4) :: temp_a             !temporary placeholder
  real*8, dimension(4) :: local_tt                 !travel time on four tet nodes
  real*8, dimension(4,3) :: t_n
  integer, dimension(4) :: indx        
  integer :: j,k

  do j=1,4
     t_n(j,:) = dble(nodes(elements(el,j),:))
  end do

  !build the shape function
  atild(1,:) = 1.0
  do j=2,4
     do k=1,4
        atild(j,k)=t_n(k,j-1)
     end do
  end do
  temp_a=atild
  
  !!invert temp_a to get atildi
  !!note MIGS routine is located in mat_inv.f  
  call MIGS(temp_a,4,atildi,indx)
  
  !get the travel times for each node in this element
  do j=1,4   
     local_tt(j) = ttimes(elements(el,j),src)
  end do
  
  !compute the travel time gradient though this element Ur
  do j=1,3
     grd(j) = dot_product(atildi(:,j+1),local_tt)
  end do
  
end subroutine compute_grad
!______________________________________________________________________________

  !_____________________________________________________________________________
  subroutine build_jaco_fmm
    implicit none
    logical :: check_faces = .true.
    logical :: upcheck,fcfound
    logical :: near_source
    logical, dimension(:), allocatable :: eligible
    
    integer :: c_row,i,j,edge,fc,nd,st_nod,rrow,fc_exit
    integer :: cur_el,my_src,src_nod,st_el,nbr,k,l,cnt,m,ofc,step
    
    real*8, dimension(4,4) :: atild              !pre-shape linear shape function matrix 
    real*8, dimension(4,4) :: atildi             !shape function matrix
    real*8, dimension(4,4) :: temp_a             !temporary placeholder
    
    integer, dimension(4) :: indx                !computation space for MIGS routine
    integer, dimension(:,:),allocatable :: A_tmp !temporary accounting matrix
    
    real*8, dimension(4,3) :: t_n
    real*8 :: evol                                 !nodal volume
    
    real*8 :: A_kij                                !forward matrix coefficient sum
    
    real*8, dimension(4) :: local_tt                 !travel time on four tet nodes
    real*8, dimension(3) :: pcur,pup,p1,p2,p3,p4,u
    real*8, dimension(3) :: v1,v2,v3,v4,Ur,Vr,Vs,pip,n
    real*8, dimension(3,3) :: Us
    integer, dimension(3) ::   edge_rot,n1,n2,n3
    integer, dimension(4,4) :: edge_seq
    real*8 :: d, sens
    
    !allocate Jaco
    nj_rows=jind(my_rank_fmm,2)-jind(my_rank_fmm,1)+1
    if(.not. allocated(Jaco)) allocate(Jaco(nj_rows,nelem))

    allocate(eligible(nelem))

    edge_seq(1,1:4) = (/ 1,2,3,1 /)
    edge_seq(2,1:4) = (/ 1,2,4,1 /)
    edge_seq(3,1:4) = (/ 1,4,3,1 /)
    edge_seq(4,1:4) = (/ 2,3,4,2 /)
    
    !loop over the survey and do the ray tracing for sources I own
    c_row = 0
    step = 0
    do i=1,nm_fmm
       do my_src=1,my_ns

          if(s_conf_fmm(i,1) .eq. sind(my_rank_fmm,1)-my_src+1)then
             eligible = .true.

             !f(s_conf_fmm(i,1).ge.sind(my_rank_fmm,1) .and. s_conf_fmm(i,1).le.sind(my_rank_fmm,2)) then
             c_row = c_row + 1
             src_nod = s_nods(s_conf_fmm(i,1))
             
             !set the starting node
             st_nod = r_nods(s_conf_fmm(i,2))
             pcur = nodes(st_nod,:)
             upcheck = .false.
             
             !loop over the elements that use this node and
             !find the one that contains this ray
             do rrow = ring_map(st_nod-1)+1,ring_map(st_nod)
             	
                cur_el = rings(rrow)
 	        do j=1,4
                   t_n(j,:) = dble(nodes(elements(cur_el,j),:))
                end do
          
                !!build the linear shape function matrix for this element
 	        atild(1,:) = 1.0
                do j=2,4
                   do k=1,4
                      atild(j,k)=t_n(k,j-1)
                   end do
                end do
 	        temp_a=atild
 	       
                
                !!invert temp_a to get atildi
                !!note MIGS routine is located in mat_inv.f  
                call MIGS(temp_a,4,atildi,indx)
                
                !get the travel times for each node in this element
                do j=1,4   
                   local_tt(j) = ttimes(elements(cur_el,j),my_src)
                end do
                
                !compute the travel time gradient though this element Ur
                do j=1,3
                   Ur(j) = dot_product(atildi(:,j+1),local_tt)
                end do
               
                
                !Ur = Ur/sqrt(dot_product(Ur,Ur))
                Vr = cross_prod(Ur,pcur)
                
                !test each face to see if the ray passes through
                do fc = 1,4
                   !test the direction of rotation between the gradient vector
                   !and each edge of this face
                   edge_rot = 0
                   do edge = 1,3
                      
                      Us(edge,:) = t_n(edge_seq(fc,edge+1),:) - t_n(edge_seq(fc,edge),:)
                      Vs = cross_prod(Us(edge,:),t_n(edge_seq(fc,edge),:))
                      pip(edge)= dot_product(Ur,Vs) + dot_product(Us(edge,:),Vr)    
                      if(pip(edge)>0) edge_rot(edge)=1
                      if(pip(edge)<0) edge_rot(edge)=-1
                      
                   end do
                   !if(my_rank==2) write(*,*) pip
                   if(sum(edge_rot)==3 .or. sum(edge_rot)==-3) then
                      !compute the normal to face fc
                      n=cross_prod(Us(2,:),Us(1,:))
                      n=n/sqrt(dot_product(n,n))
                      
                      !compute the distance from the current point to the face
                      if(dot_product(n,Ur) .ne. 0) then
                         d = dot_product(t_n(edge_seq(fc,1),:)-pcur,n)/dot_product(Ur,n)
                         sens=sqrt(dot_product(d*Ur,d*Ur))
                         !if(my_rank == 2) write(*,*) d,sens,'......................'
                      end if
                     
                      if(d<0 .and. sens>1e-6) then
                         
                         pup = pcur + d*Ur
                         upcheck = .true.
                         if(check_faces) then
                            v1=cross_prod(Us(1,:),pup-t_n(edge_seq(fc,1),:))
                            v2=cross_prod(Us(2,:),pup-t_n(edge_seq(fc,2),:))
                            v3=cross_prod(Us(3,:),pup-t_n(edge_seq(fc,3),:))
                            v1=v1/sqrt(dot_product(v1,v1))
                            v2=v2/sqrt(dot_product(v2,v2))
                            v3=v3/sqrt(dot_product(v3,v3))
                            
                            !write(*,*) my_rank_fmm,sens
                            !write(*,*) pup
                            if( abs(dot_product(v1,v2)-1) >1e-6  .or. &
                                 abs(dot_product(v2,v3)-1)>1e-6 .or. &
                                 abs(dot_product(v3,v1)-1)>1e-6) then
                               write(*,*) 'Error locating the right point'
                               upcheck = .false.
                            end if
                         end if
                      
                         goto 100
                      end if
                   end if
                end do
		
             end do
          
            
100          continue   !we're sent here if the ray was initialized
             
             step=step+1
             if(.not. upcheck) then
             	write(*,*) my_rank_fmm, "Unable to start ray tracing for: ", s_conf_fmm(i,:)
                write(*,*) "Measurement: ", s_conf_fmm(i,:)," will be ignored"
                goto 200
             end if
           
             !set the sensitivity for the current element
             Jaco(c_row,cur_el) = sens
             eligible(cur_el) = .false.
             pcur = pup
             st_el = cur_el
             
             !here we continue the ray tracing until we're in an element that contains
             !the source node
             do rrow = ring_map(src_nod-1)+1,ring_map(src_nod)
                if(st_el .eq. rings(rrow)) then
                   write(*,*) 'finishing this ray'
                   goto 110
                end if
             end do
             
             !find the element that shares this face
             !n1 holds the nodes of the face that contains the pcur 
           
             n1 = elements(st_el,edge_seq(fc,1:3))
             fcfound = .false.
             do j=1,4 !loop over neighbors
                nbr = nbrs(st_el,j)  !jth neigbor to downwind element
                do m=1,4 !loop over faces
                   n2 = elements(nbr,edge_seq(m,1:3)) !face m of neighbor
                   !check to see if n2 has the same nodes n1
                   cnt = 0
                   do k=1,3
                      if(n2(k) .eq. n1(1) .or. n2(k) .eq. n1(2) .or. n2(k) .eq. n1(3)) then
                         cnt=cnt+1
                      end if
                   end do
                   if(cnt .eq. 3) then
                      fc = m
                      cur_el = nbr
                      fcfound = .true.
                      !if(.not.eligible(nbr)) then
                      !   do k=1,3
                      !      v1(1,:) = .25*sum(nodes(elements(st_el,k),1))
!!$                      if(step .eq. 40) then
!!$                         if(my_rank_fmm==1) then
!!$                           
!!$                               
!!$                              
!!$                               do k=1,4
!!$                                  write(21,*) nodes(elements(st_el,k),:),1
!!$                                  write(21,*) nodes(elements(nbr,k),:),2 
!!$                               end do
!!$                               
!!$                               write(21,*) nodes(n1(1),:),3
!!$                               write(21,*) nodes(n1(2),:),3
!!$                               write(21,*) nodes(n1(3),:),3
!!$                               write(21,*) nodes(n2(1),:),4
!!$                               write(21,*) nodes(n2(2),:),4
!!$                               write(21,*) nodes(n2(3),:),4
!!$                               write(21,*) pcur,6
!!$                               write(*,*) fc,nbr,step
!!$                               return
!!$                            
!!$
!!$                            end if
!!$                         end if
                      goto 101
                   end if
               end do
             end do
             write(*,*) my_rank_fmm," couldn't find the right element"
101 continue
             !if(my_rank .eq. 4) write(*,*) st_el,cur_el
             if(.not.fcfound) then
                write(*,*) my_rank_fmm,"couldn't find the face"
                return
             else
                !write(*,*) n1
                !write(*,*) n2
             end if
              
             !compute the gradient in element cur_el______________________________
             do j=1,4
                t_n(j,:) = dble(nodes(elements(cur_el,j),:))
             end do
             
             atild(1,:) = 1.0
             do j=2,4
                do k=1,4
                   atild(j,k)=t_n(k,j-1)
                end do
             end do
             
             temp_a=atild
             call MIGS(temp_a,4,atildi,indx)
             
             do j=1,4   
                local_tt(j) = ttimes(elements(cur_el,j),my_src)
             end do
             
             do j=1,3
                Ur(j) = dot_product(atildi(:,j+1),local_tt)
             end do
             Vr = cross_prod(Ur,pcur)
             !_______________________________________________________________________
           
            
       
             !loop over the other 3 faces to determine which the ray exits
             fcfound = .false.
             do j=1,4
                if(j .ne. fc) then
                   edge_rot = 0
                   do edge = 1,3   
                      Us(edge,:) = t_n(edge_seq(j,edge+1),:) - t_n(edge_seq(j,edge),:)
                      Vs = cross_prod(Us(edge,:),t_n(edge_seq(j,edge),:))
                      pip(edge)= dot_product(Ur,Vs) + dot_product(Us(edge,:),Vr)    
                      if(pip(edge)>0) edge_rot(edge)=1
                      if(pip(edge)<0) edge_rot(edge)=-1   
                   end do 
                   if(my_rank_fmm == 1) write(*,*) step,pip
                   if(sum(edge_rot)==3 .or. sum(edge_rot) == -3) then
                      fcfound = .true.
                      fc_exit= j
                      
                      do m=1,3
                         u(m)=pip(m)/sum(pip)
                      end do
                      v1= u(1)*t_n(edge_seq(j,1),:)+u(2)*t_n(edge_seq(j,2),:)+u(3)*t_n(edge_seq(j,3),:) 
                      exit
                   end if
                end if
             end do
             if(.not.fcfound) then
                write(*,*) my_rank_fmm,"couldn't find the exit face"
                return
             end if

             n=cross_prod(Us(2,:),Us(1,:))
             n=n/sqrt(dot_product(n,n))

             !compute the distance from the current point to the face
             if(dot_product(n,Ur) .ne. 0) then
                d = dot_product(t_n(edge_seq(fc_exit,1),:)-pcur,n)/dot_product(Ur,n)
                sens=sqrt(dot_product(d*Ur,d*Ur))
             end if
             !write(*,*) my_rank_fmm, d, dot_product(n,Ur),sens
             fc = fc_exit
             pup = pcur + d*Ur       
             
          
             
             
             if(my_rank_fmm == 1) then
                !write(20,*)
                !write(20,*) pup
                !write(20,*) v1; !-Ur/sqrt(dot_product(Ur,Ur))
                !write(20,*) v2; !(pup-pcur)/sqrt(dot_product(pup-pcur,pup-pcur))
                !write(20,*) dot_product(Ur,n)
                write(20,*) t_n(edge_seq(j,1),:),step,1
                write(20,*) t_n(edge_seq(j,2),:),step,1
                write(20,*) t_n(edge_seq(j,3),:),step,1
                write(20,*) pup,2,step
                write(20,*) 
                !write(*,*) step,cur_el,st_el
                !write(*,*) d,Ur
                !write(*,*)
             end if
             goto 100

             !Trace the ray through element st_el
             !upcheck = .false.
!!$             
!!$        
!!$             do j=1,4
!!$                cnt = 0
!!$                cur_el = nbrs(st_el,j) 
!!$                do m=1,3
!!$                   k=edge_seq(fc,m)
!!$                   do l=1,4
!!$                      if(elements(cur_el,l) .eq. elements(st_el,k)) then
!!$                         cnt=cnt+1
!!$                         exit
!!$                      end if
!!$                   end do    
!!$                end do
!!$                if(cnt==3) exit
!!$             end do
!!$
!!$             !find the face on the current element that holds the point
!!$             ofc=fc
!!$             n1=elements(st_el,edge_seq(ofc,1:3))
!!$             fc = 0
!!$             do j=1,4
!!$                cnt = 0
!!$                n2=elements(cur_el,edge_seq(j,1:3))
!!$                do m=1,3
!!$                   if(n2(m).eq.n1(1) .or. n2(m).eq.n1(2) .or. n2(m).eq.n1(3)) cnt=cnt+1
!!$                end do
!!$                if(cnt==3) then
!!$                   fc = j
!!$                   exit
!!$                end if
!!$             end do
!!$             if(fc.eq.0) then
!!$                write(*,*) my_rank, "couldn't find the common face"
!!$                return
!!$             end if
!!$
!!$             !compute the gradient in element cur_el______________________________
!!$             do j=1,4
!!$                t_n(j,:) = dble(nodes(elements(cur_el,j),:))
!!$             end do
!!$
!!$             atild(1,:) = 1.0
!!$             do j=2,4
!!$                do k=1,4
!!$                   atild(j,k)=t_n(k,j-1)
!!$                end do
!!$             end do
!!$
!!$             temp_a=atild
!!$             call MIGS(temp_a,4,atildi,indx)
!!$                   
!!$             do j=1,4   
!!$                local_tt(j) = ttimes(elements(cur_el,j),my_src)
!!$             end do
!!$            
!!$             do j=1,3
!!$                Ur(j) = dot_product(atildi(:,j+1),local_tt)
!!$             end do
!!$             Vr = cross_prod(Ur,pcur)
!!$             !_______________________________________________________________________
!!$             
!!$             !loop over the other 3 faces to determine which the ray exits
!!$             fcfound = .false.
!!$             do j=1,4
!!$                if(j .ne. fc) then
!!$                   edge_rot = 0
!!$                   do edge = 1,3   
!!$                      Us(edge,:) = t_n(edge_seq(j,edge+1),:) - t_n(edge_seq(j,edge),:)
!!$                      Vs = cross_prod(Us(edge,:),t_n(edge_seq(j,edge),:))
!!$                      pip(edge)= dot_product(Ur,Vs) + dot_product(Us(edge,:),Vr)    
!!$                      if(pip(edge)>0) edge_rot(edge)=1
!!$                      if(pip(edge)<0) edge_rot(edge)=-1   
!!$                   end do 
!!$                end if
!!$                if(my_rank .eq. 4) write(*,*) pip
!!$                if(sum(edge_rot)==3 .or. sum(edge_rot) == -3) then
!!$                   fcfound = .true.
!!$                   fc= j
!!$                   exit
!!$                end if
!!$             end do
!!$             if(.not.fcfound) then
!!$                write(*,*) my_rank_fmm,"couldn't find the exit face"
!!$                return
!!$             end if
!!$             
!!$             n=cross_prod(Us(2,:),Us(1,:))
!!$             n=n/sqrt(dot_product(n,n))
!!$                         
!!$             !compute the distance from the current point to the face
!!$             if(dot_product(n,Ur) .ne. 0) then
!!$                d = dot_product(t_n(edge_seq(fc_exit,1),:)-pcur,n)/dot_product(Ur,n)
!!$                sens=sqrt(dot_product(d*Ur,d*Ur))
!!$             end if
!!$             pup = pcur + d*Ur
!!$          
!!$
!!$             if(my_rank_fmm==4) then
!!$                write(*,*) st_el,cur_el,sens
!!$                write(*,*) pup
!!$             end if
!!$             goto 100
!!$             
!!$             
!!$
!!$
!!$
!!$
!!$
!!$
!!$








!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$
!!$             do m=1,4
!!$                cur_el = nbrs(st_el,m)
!!$                if(eligible(cur_el)) then
!!$                   if(my_rank==4) then
!!$                      write(*,*) st_el,cur_el
!!$                      if(.not.eligible(cur_el)) then
!!$                         cur_el = st_el
!!$                         !call find_me(pcur,Ur,cur_el)
!!$                         !write(*,*) st_el,cur_el
!!$                      end if
!!$                      
!!$                      !do m=1,4
!!$                      !write(*,*) nodes(elements(st_el,m),:),1
!!$                      !end do
!!$                   end if
!!$                   
!!$                   
!!$                   !comput the tt gradient
!!$                   do j=1,4
!!$                      t_n(j,:) = dble(nodes(elements(cur_el,j),:))
!!$                      !if(my_rank == 4) write(*,*) t_n(j,:),2
!!$                   end do
!!$                   !if(my_rank == 4) write(*,*) pcur,3
!!$                   
!!$                   !!build the linear shape function matrix for this element
!!$                   atild(1,:) = 1.0
!!$                   do j=2,4
!!$                      do k=1,4
!!$                         atild(j,k)=t_n(k,j-1)
!!$                      end do
!!$                   end do
!!$                   temp_a=atild       
!!$                   
!!$                   !!invert temp_a to get atildi
!!$                   !!note MIGS routine is located in mat_inv.f  
!!$                   call MIGS(temp_a,4,atildi,indx)
!!$                   
!!$                   !get the travel times for each node in this element
!!$                   do j=1,4   
!!$                      local_tt(j) = ttimes(elements(cur_el,j),my_src)
!!$                   end do
!!$                   
!!$                   !compute the travel time gradient though this element Ur
!!$                   do j=1,3
!!$                      Ur(j) = dot_product(atildi(:,j+1),local_tt)
!!$                   end do
!!$                   Vr = cross_prod(Ur,pcur)
!!$                   
!!$                   !test each face to see if the ray passes through
!!$                   
!!$                   do fc = 1,4
!!$                      !test the direction of rotation between the gradient vector
!!$                      !and each edge of this face
!!$                      edge_rot = 0
!!$                      do edge = 1,3
!!$                         
!!$                         Us(edge,:) = t_n(edge_seq(fc,edge+1),:) - t_n(edge_seq(fc,edge),:)
!!$                         Vs = cross_prod(Us(edge,:),t_n(edge_seq(fc,edge),:))
!!$                         pip(edge)= dot_product(Ur,Vs) + dot_product(Us(edge,:),Vr)    
!!$                         if(pip(edge)>0) edge_rot(edge)=1
!!$                         if(pip(edge)<0) edge_rot(edge)=-1
!!$                         
!!$                      end do
!!$                      !if(my_rank_fmm==4) write(*,*) pip                 
!!$                      if(sum(edge_rot)==3 .or. sum(edge_rot)==-3) then
!!$                         !compute the normal to face fc
!!$                         n=cross_prod(Us(2,:),Us(1,:))
!!$                         n=n/sqrt(dot_product(n,n))
!!$                         
!!$                         !compute the distance from the current point to the face
!!$                         if(dot_product(n,Ur) .ne. 0) then
!!$                            d = dot_product(t_n(edge_seq(fc,1),:)-pcur,n)/dot_product(Ur,n)
!!$                            sens=sqrt(dot_product(d*Ur,d*Ur))
!!$                            
!!$                         end if
!!$                         !if(my_rank_fmm ==4) write(*,*) Ur,d,sens
!!$                         if(d<0 .and. sens>1e-12) then
!!$                            if(my_rank_fmm ==4) write(*,*) Ur,d,sens
!!$                            pup = pcur + d*Ur
!!$                            upcheck = .true.
!!$                            if(check_faces) then
!!$                               v1=cross_prod(Us(1,:),pup-t_n(edge_seq(fc,1),:))
!!$                               v2=cross_prod(Us(2,:),pup-t_n(edge_seq(fc,2),:))
!!$                               v3=cross_prod(Us(3,:),pup-t_n(edge_seq(fc,3),:))
!!$                               v1=v1/sqrt(dot_product(v1,v1))
!!$                               v2=v2/sqrt(dot_product(v2,v2))
!!$                               v3=v3/sqrt(dot_product(v3,v3))
!!$                               
!!$                               if( abs(dot_product(v1,v2)-1) >1e-6  .or. &
!!$                                    abs(dot_product(v2,v3)-1)>1e-6 .or. &
!!$                                    abs(dot_product(v3,v1)-1)>1e-6) then
!!$                                  write(*,*) 'Error locating the right point'
!!$                                  upcheck = .false.
!!$                               end if
!!$                            end if
!!$                            
!!$                            if(my_rank == 4) write(*,*) pup,4
!!$                            if(upcheck) goto 100 
!!$                         end if
!!$                      end if
!!$                      
!!$                   end do
!!$                end if
!!$             end do
!!$             goto 300
             
          end if
          
110       continue   ! we're sent here if we're at an element that contains the source
         
200       continue   
       end do
       
    end do
300 continue
  
end subroutine build_jaco_fmm
!_____________________________________________________________________________

!_____________________________________________________________________________
subroutine find_me(pt,grd,el)
  implicit none
  real, dimension(3), intent(in) :: pt,grd
  real, dimension(3) :: nrm
  integer :: el
  integer, dimension(1) :: indx,indx2
  real, dimension(nnodes) :: dist
  real,dimension(:,:), allocatable :: mids
  integer :: ntet,i
  
  !find the nearest node
  dist = sqrt((nodes(:,1)-pt(1))**2 + (nodes(:,2)-pt(2))**2 + (nodes(:,3)-pt(3))**2)  
  indx = minloc(dist)

  !count the number of tets using this node
  ntet = ring_map(indx(1))-ring_map(indx(1)-1)

  !compute the midpoints of the tets
  allocate(mids(ntet,3))
  do i=1,ntet
     el=rings(ring_map(indx(1)-1)+i)
     mids(i,1) = .25*sum(nodes(elements(el,:),1))
     mids(i,2) = .25*sum(nodes(elements(el,:),2))
     mids(i,3) = .25*sum(nodes(elements(el,:),3))
  end do
  
  !now build vectors from pt to each of the midpoints
  do i=1,ntet
    mids(i,:) = pt-mids(i,:)
    dist(i) = sqrt(dot_product(mids(i,:),mids(i,:)))  
    mids(i,:) = mids(i,:)/dist(i)
  end do

  !compute the normal vector for the grad
  nrm=-grd/sqrt(dot_product(grd,grd))
  do i=1,ntet
     dist(i)=dot_product(nrm,mids(i,:))
  end do
  indx2 = maxloc(dist(1:ntet))
  el=rings(ring_map(indx(1)-1)+indx2(1))

  deallocate(mids)
  return
end subroutine find_me
!_____________________________________________________________________________

!____________________________________________________________________
  !compute the cross product between two vectors
  function cross_prod(u,v)

    implicit none
    real*8, dimension(3), intent(in) :: u,v
    real*8, dimension(3) :: cross_prod
    cross_prod(1) = u(2)*v(3) - u(3)*v(2)
    cross_prod(2) = u(3)*v(1) - u(1)*v(3)
    cross_prod(3) = u(1)*v(2) - u(2)*v(1)
  
  end function cross_prod
!____________________________________________________________________


end module jaco_fmm
