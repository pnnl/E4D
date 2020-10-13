module forward_fmm
 
  use fmm_vars
  use vars
  
  implicit none   
  
contains
  !____________________________________________________________________
  subroutine forward_run_fmm
    implicit none
    integer :: cur_src,s,i,j,k,snod,dpos,ndi,nring,cur_ele,tnode,nct
    integer, dimension(4) :: cur_nd_ups
    logical :: rec_done_flag
    real :: pd
    integer :: heap_map_tmp, heap_tmp
!    real, dimension(1,1) :: tt_test
    do s=1,my_ns
       cur_src = sind(my_rank_fmm,1)+s-1
       snod = s_nods(cur_src)
       call init_heap(snod)
       
       
       rec_done_flag = .false.
       !here we go
      
      
       !keep going if any receivers are not upstream
       !do while(.not.rec_done_flag)
       do while(heap_end.gt.0)
          !move the first node in the heap to upstream
          num_heap_done = num_heap_done + 1
          dpos = nnodes - num_heap_done + 1
!          heap(dpos) = heap(1)
!          heap_map(heap(1)) = dpos
! to avoid overlap due to bottom up (upstream) and top down (close) of heap
          heap_tmp = heap(1)
          heap_map_tmp = dpos

          pd = real(num_heap_done)/real(nnodes)
          !write(*,*) pd*1000
          
          ndi = heap(1)
          upstream(ndi) = .true.
          clos(ndi) = .false.
          heap(1) = heap(heap_end)
          heap_map(heap(1)) = 1
          heap_end = heap_end - 1
          call sift_down
          heap(dpos) = heap_tmp
          heap_map(heap_tmp) = heap_map_tmp

          if(ndi .ne. 1) then
!             nring = ring_map(ndi) - ring_map(ndi-1) - 1
             nring = ring_map(ndi) - ring_map(ndi-1)
          else
             nring = ring_map(1)
          end if
          
          do i=1,nring
             cur_ele = rings(ring_map(ndi-1)+i)
             cur_nd_ups(1:4) = 0
             do k = 1,4
                tnode = elements(cur_ele,k)
                if(upstream(tnode)) cur_nd_ups(k) = 1
             end do
             do k=1,4
              tnode = elements(cur_ele,k)
              if(.not. upstream(tnode) .and. sum(cur_nd_ups) < 4) then
                 
                 !set the location and travel time of this node
                 !in position 4 
                 XP(4,1:3) = nodes(tnode,1:3)
!                 LTT(4) = TT(tnode)
                 nct = 0

                 !get the positions and travel times of the other nodes
                 do j=1,4
                    if(j.ne.k) then
                       nct=nct+1
                       XP(nct,1:3) = nodes(elements(cur_ele,j),1:3)
                       LTT(nct) = TT(elements(cur_ele,j))
                    end if
                 end do
                 M(1,1) = speed(cur_ele); M(2,2)=speed(cur_ele); M(3,3)=speed(cur_ele)
                 call local_tt(xp(1,:),xp(2,:),xp(3,:),xp(4,:),ltt(1),ltt(2),ltt(3),M,tt_test(1,1))
     
                 !if the test travel time is less than the current, update TT and heap
                 if(TT(tnode) .gt. tt_test(1,1)) then
                    TT(tnode) = tt_test(1,1)
                    
                    !if tnode isn't close yet then add it to the heap
                    if(.not.clos(tnode)) then
                       clos(tnode) = .true.
                       heap_end = heap_end+1
                       heap(heap_end) = tnode
                       heap_map(tnode) = heap_end
                    end if
 
                    call bubble_up(tnode)
              
                 end if
              end if
           end do
        end do
        !check to see if all the receiver positions are upstream
        call check_recs(rec_done_flag)


!     do i=1,nnodes
!          j=nnodes-i+1
!          ttimes(heap(j),s) = TT(heap(j))        
!     end do
     end do
     ttimes(:,s) = tt(:)
!check if there are nodes not reached
!     do cur_ele = 1,nelem 
!        i = 0
!        do k = 1,4
!           tnode = elements(cur_ele,k)
!           if(upstream(tnode)) i = i+1
!        enddo
!        if(i /= 4) print *,'cur_ele',cur_ele,i
!     enddo
          
  end do
 
!  if(my_rank_fmm==1) then
!     do i=1,nnodes
!        j=nnodes-i+1
!        write(20,*) i, j, TT(heap(j)), heap(j)
!     end do
!  end if
  return



  end subroutine forward_run_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine check_recs(rdf)
    implicit none
    logical :: rdf
    integer :: i
    rdf = .true.
    do i=1,nrc
       if(.not. upstream(r_nods(i))) then
          rdf = .false.
          exit
       end if
    end do
  end subroutine check_recs
  !____________________________________________________________________
  !____________________________________________________________________
  subroutine init_heap(snode)
    implicit none
    integer, intent(in) :: snode
    integer :: nring,i,j,cur_ele,tnode,nd,ndi,nct,k,kx
    integer, dimension(4) :: cur_nd_ups
    real, dimension(3,1) :: e_vec,t_vec
!    real, dimension(1,1) :: tt_test 
    real :: tmp

    if(.not. allocated(heap)) allocate(heap(nnodes))
    heap = 0

    if(.not. allocated(TT)) allocate(TT(nnodes)) 
    TT=tt_max

    if(.not. allocated(heap_map)) allocate(heap_map(nnodes))
    heap_map = 0

    if(.not. allocated(upstream)) allocate(upstream(nnodes))
    upstream = .false.

    if(.not. allocated(clos)) allocate(clos(nnodes))
    clos = .false.

    heap_end = 0
!
    num_heap_done = 0
    
    !loop over all the elements that have snode as a vertex
    !and compute the travel time to node connected to snode

    !nring is the number of elements that have snode as a vertex
    if(snode .ne. 1) then
!       nring = ring_map(snode) - ring_map(snode-1) - 1
       nring = ring_map(snode) - ring_map(snode-1)
    else
       nring = ring_map(1)
    end if
    TT(snode) = 0
    upstream(snode) = .true.
    heap(nnodes) = snode
    heap_map(snode) = nnodes
    num_heap_done = num_heap_done+1

    M=0
    do i=1,nring
       cur_ele = rings(ring_map(snode-1)+i)
       do j=1,4
          tnode = elements(cur_ele,j)
          if(.not. upstream(tnode)) then
             
             !compute the test travel time to this tnode
             e_vec(1:3,1) = nodes(tnode,1:3) - nodes(snode,1:3)
             M(1,1) = speed(cur_ele); M(2,2)=speed(cur_ele); M(3,3)=speed(cur_ele)
             tt_test = sqrt(MATMUL(TRANSPOSE(e_vec),MATMUL(M,e_vec)))
             
             !if the test travel time is less than the current, update TT and heap
             if(TT(tnode)>tt_test(1,1)) then
                TT(tnode) = tt_test(1,1)

                !if tnode isn't close yet then add it to the heap
                if(.not.clos(tnode)) then
                   clos(tnode) = .true.
                   heap_end = heap_end+1
                   heap(heap_end) = tnode
                   heap_map(tnode) = heap_end
                end if
             end if
          end if
       end do
    end do

    !put the heap in min heap order
    call heapify

    !move the elements of the heap to upstream in order
    do while(heap_end>0)
       heap(nnodes-num_heap_done) = heap(1)
       heap_map(heap(1)) = nnodes-num_heap_done
       num_heap_done = num_heap_done+1
       clos(heap(1)) = .false.
       upstream(heap(1)) = .true.
       if(heap_end .eq. 1) then
         heap_end = heap_end-1
         exit
       endif
       heap_map(heap(heap_end)) = 1
       heap(1) = heap(heap_end)
       heap_end = heap_end-1
       call sift_down
    end do
    
    !now update the heap with the close nodes
    do nd=1,num_heap_done 
       ndi = heap(nnodes-nd+1)
       if(ndi .ne. 1) then
!          nring = ring_map(ndi) - ring_map(ndi-1) - 1
          nring = ring_map(ndi) - ring_map(ndi-1)
       else
          nring = ring_map(1)
       end if
       
       do i=1,nring
           cur_ele = rings(ring_map(ndi-1)+i)
           cur_nd_ups(1:4) = 0
           do k = 1,4
              tnode = elements(cur_ele,k)
              if(upstream(tnode)) cur_nd_ups(k) = 1
           end do
           do k=1,4
              tnode = elements(cur_ele,k)
              if(.not. upstream(tnode) .and. sum(cur_nd_ups) < 4) then
                 !set the location and travel time of this node
                 !in position 4 
                 XP(4,1:3) = nodes(tnode,1:3)
!                 LTT(4) = TT(tnode)
                 nct = 0

                 !get the positions and travel times of the other node
                 do j=1,4
                    if(j.ne.k) then
                       nct=nct+1
                       XP(nct,1:3) = nodes(elements(cur_ele,j),1:3)
                       LTT(nct) = TT(elements(cur_ele,j))
                    end if
                 end do
                 M(1,1) = speed(cur_ele); M(2,2)=speed(cur_ele); M(3,3)=speed(cur_ele)
                 call local_tt(xp(1,:),xp(2,:),xp(3,:),xp(4,:),ltt(1),ltt(2),ltt(3),M,tt_test(1,1))
     
                 !if the test travel time is less than the current, update TT and heap
                 if(TT(tnode)>tt_test(1,1)) then
                    TT(tnode) = tt_test(1,1)
                    
                    !if tnode isn't close yet then add it to the heap
                    if(.not.clos(tnode)) then
                       clos(tnode) = .true.
                       heap_end = heap_end+1
                       heap(heap_end) = tnode
                       heap_map(tnode) = heap_end
                    end if
                 end if
              end if
           end do
        end do
     end do
     
     !get the heap in order
     call heapify
  
  end subroutine init_heap
  !____________________________________________________________________

  !____________________________________________________________________
! local solver 
  recursive subroutine local_tt(x1,x2,x3,x4,t1,t2,t3,M,l_tt)
    implicit none
    integer :: i
    real, dimension(3) :: x1,x2,x3,x4
    real, dimension(3,3), intent(in) :: M
    real, intent(out) :: l_tt
    real :: t1,t2,t3
    real, dimension(3) :: vec,tsplit
    real, dimension(3,1) :: e41,e42,e43,e54,nvec,e13,e23,e34
    real, dimension(1,1):: tmp_
!    real, dimension(1,1) :: tt_test
    real :: ltt(3),xp(4,3), tt(3), ttmp(3), x5(3), lttx(3)
    real :: m11(1,1), m21(1,1), m31(1,1), m22(1,1), m32(1,1),m33(1,1)
    real :: ttf, sa, tmp
    real :: pi
    real :: lam1(2),lam2(2),lam3(2)
    real :: x1f(3), x2f(3), x3f(3)
    real :: big_solid
    real :: t5, t13, t23, c1, c2, c3, c4, c5, c6, c7, a, b, c, d, t1f, t2f
    real :: t41, t42, t43
! if any of the 3 points haven't been touched yet then
!    the smallest travel time is the one that's direct from
!    the upstream point
    if(t1 >= 1e12 .or. t2 >= 1e12 .or. t3 >= 1e12) then
      e41(:,1) = x4(:) - x1(:)
      tmp_ = MATMUL(TRANSPOSE(e41(:,:)),MATMUL(M(:,:),e41(:,:)))
      ttf = sqrt(tmp_(1,1))
      t41 = t1 + ttf
      e42(:,1) = x4(:) - x2(:)
      tmp_ = MATMUL(TRANSPOSE(e42(:,:)),MATMUL(M(:,:),e42(:,:)))
      ttf = sqrt(tmp_(1,1))
      t42 = t2 + ttf
      e43(:,1) = x4(:) - x3(:)
      tmp_ = MATMUL(TRANSPOSE(e43(:,:)),MATMUL(M(:,:),e43(:,:)))
      ttf = sqrt(tmp_(1,1))
      t43 = t3 + ttf
      tt_test(1,1) = min(min(t41,t42),t43)
      ttf = tt_test(1,1)
      l_tt = ttf
      return
    endif
! if the front is normal to the x1,x2,x3 plane
    if(t1 == t2 .and. t1 == t3) then
      vec = cross_prod(x3-x1,x2-x1)
      e41(:,1) = x4(:) - x1(:)
      nvec(:,1) = vec/NORM2(vec)
      e54 = MATMUL(nvec(:,:),MATMUL(TRANSPOSE(e41(:,:)),nvec(:,:)))
      tmp_ = MATMUL(TRANSPOSE(e54(:,:)),MATMUL(M(:,:),e54(:,:)))
! shouldn't t1 be added to?
      ttf = sqrt(tmp_(1,1)) + t1
!      ttf = sqrt(tmp_(1,1))
      tt_test(1,1) = ttf  
      l_tt = ttf
    endif

    if(t1 == t3) then
      tmp = t2
      t2 = t1
      t1 = tmp
      ttmp = x2
      x2 = x1
      x1 = ttmp
    endif

!    check the solid angle at x4 to see if we need to divide and conquer
    sa = solid_angle(x1,x2,x3,x4,4)
    pi = 4.0*ATAN(1.0)
    if(sa > pi/2) then
!        these 4 lines need to be changed to use barycentric coords
!        disp('big solid angle')
      big_solid = big_solid + 1
      x5(1) = (x1(1)+x2(1)+x3(1))/3.0
      x5(2) = (x1(2)+x2(2)+x3(2))/3.0
      x5(3) = (x1(3)+x2(3)+x3(3))/3.0
      t5 = (t1 + t2+ t3)/3

      call local_tt(x1,x2,x5,x4,t1,t2,t5,M,l_tt)
      tsplit(1) = l_tt
      call local_tt(x2,x3,x5,x4,t2,t3,t5,M,l_tt)
      tsplit(2) = l_tt
      call local_tt(x3,x1,x5,x4,t3,t1,t5,M,l_tt)
      tsplit(3) = l_tt

      ttf = minval(tsplit)
      tt_test(1,1) = ttf  
      l_tt = ttf
      return
    end if
    e13(:,1) = x3(:) - x1(:)
    e23(:,1) = x3(:) - x2(:)
    e34(:,1) = x4(:) - x3(:)
    M11 = MATMUL(TRANSPOSE(e13(:,:)),MATMUL(M,e13(:,:)))
    M21 = MATMUL(TRANSPOSE(e23(:,:)),MATMUL(M,e13(:,:)))
    M31 = MATMUL(TRANSPOSE(e34(:,:)),MATMUL(M,e13(:,:)))
    M22 = MATMUL(TRANSPOSE(e23(:,:)),MATMUL(M,e23(:,:)))
    M32 = MATMUL(TRANSPOSE(e23(:,:)),MATMUL(M,e34(:,:)))
    M33 = MATMUL(TRANSPOSE(e34(:,:)),MATMUL(M,e34(:,:)))

    t13=t3-t1
    t23=t3-t2

    C1 = (t13*M21(1,1)-t23*M22(1,1))/(t23*M21(1,1) - t13*M11(1,1))
    C2 = (t13*M31(1,1)-t23*M32(1,1))/(t23*M21(1,1) - t13*M11(1,1))

    C3 = M11(1,1)*C1 + 2*M21(1,1)*C1 + M22(1,1)
    C4 = 2*M11(1,1)*C1*C2 + 2*M21(1,1)*C2 + 2*M31(1,1)*C1+2*M32(1,1)
    C5 = M11(1,1)*C2*C2 + 2*M31(1,1)*C2 + M33(1,1)

    C6 = (C1*M11(1,1)+M21(1,1))/(t13)
    C7 = (C2*M11(1,1)+M31(1,1))/(t13)

    a = C3 - C6*C6
    b = C4 - 2*C6*C7
    c = C5 - C7*C7
    
    lttx(:) = 1e12
    d = b*b - 4*a*c
    if(d >= 0) then
      lam2(1) = (-b + sqrt(d))/(2*a)
      lam2(2) = (-b - sqrt(d))/(2*a)
      lam1(1) = lam2(1)*C1+C2
      lam1(2) = lam2(2)*C1+C2
      do i=1,2
         if(lam1(i)>=0 .and. lam1(i)<=1 .and. lam2(i)>=0 .and. lam2(i)<=1) then
             lam3(i)=1-lam1(i)-lam2(i)
             x5 = x1*lam1(i) + x2*lam2(i) + x3*lam3(i)
             e54(:,1) = x4-x5
             tmp_ = MATMUL(TRANSPOSE(e54(:,:)), MATMUL(M(:,:),e54(:,:)))
             lttx(i)=t1*lam1(i) + t2*lam2(i) + t3*lam3(i) +sqrt(tmp_(1,1))
          else
             ttmp(1)=face_edge_tt(x1,x2,x4,t1,t2,M)
             ttmp(2)=face_edge_tt(x2,x3,x4,t2,t3,M)
             ttmp(3)=face_edge_tt(x3,x1,x4,t3,t1,M)
             lttx(i)=minval(ttmp)
          endif
      end do
      ttf=minval(lttx)
      tt_test(1,1) = ttf  
      l_tt = ttf
      return
    else
!       face and edge solutions go here
      x1f=x1
      x2f=x2
      x3f=x4
      t1f=t1
      t2f=t2
      ttmp(1)=face_edge_tt(x1,x2,x4,t1,t2,M)
      x1f=x2
      x2f=x3
      t1f=t2
      t2f=t3
      ttmp(2)=face_edge_tt(x2,x3,x4,t2,t3,M)
      x1f=x3
      x2f=x1
      t1f=t3
      t2f=t1
      ttmp(3)=face_edge_tt(x3,x1,x4,t3,t1,M)
      ttf=minval(ttmp)
      tt_test(1,1) = ttf  
      l_tt = ttf
    end if
    
  end subroutine local_tt
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine bubble_up(tnd)
    implicit none
    integer, intent(in) :: tnd
    integer :: tmp,parent,child

    if(heap_map(tnd) .gt. 1) then
       
       child = heap_map(tnd)
       parent = child/2

       do while(TT(heap(parent)) .gt.  TT(heap(child)))
     
          tmp=heap(child)
          heap_map(heap(child)) = parent
          heap_map(heap(parent))= child 
          heap(child) = heap(parent)
          heap(parent) = tmp
          child=parent
          parent=parent/2
          
          if(parent<1) then
             exit
          end if

       end do
    end if
   
  end subroutine bubble_up
  !____________________________________________________________________
  !____________________________________________________________________
  subroutine sift_down
    implicit none
    integer :: parent, child,tmp

    parent = 1
    child = 2*parent
    if(child+1 .le. heap_end) then
       if(TT(heap(child)) .gt. TT(heap(child+1))) then
          child = child+1
       end if
    end if

    if(child .le. heap_end) then

       do while(TT(heap(parent)) .gt. TT(heap(child)))

          tmp=heap(child)
          heap_map(heap(child)) = parent
          heap_map(heap(parent)) = child
          heap(child) = heap(parent)
          heap(parent) = tmp
          parent = child
          child = 2*parent

          if(child .gt. heap_end) then
             exit
          end if

          if(child+1 .le. heap_end) then
             if(TT(heap(child)) .gt. TT(heap(child+1))) then
                child = child + 1
             end if
          end if

       end do

    end if
  end subroutine sift_down
  !____________________________________________________________________
  !____________________________________________________________________
  subroutine heapify
    implicit none
    integer :: start,parent, child,bottom,i,tmp
 
    do i=1,(heap_end-2)/2+1
       start = (heap_end-2)/2+2-i
       parent = start
       bottom = heap_end

       do while(2*parent .le. bottom)
          
          child = 2*parent
          
          if(child+1 .le. bottom) then
             if(TT(heap(child)) .gt. TT(heap(child+1))) then
                child = child+1
             end if
          end if

          if(TT(heap(parent)) > TT(heap(child))) then
             tmp = heap(child)
             heap_map(heap(child)) = parent
             heap_map(heap(parent)) = child
             heap(child) = heap(parent)
             heap(parent) = tmp
             parent = child
          else
             exit
          end if
       end do
    end do
    
  end subroutine heapify
  !____________________________________________________________________

  !____________________________________________________________________

  !____________________________________________________________________
  subroutine mget_source_nodes(passed)
    implicit none
    logical :: passed
    real, dimension(nnodes) :: dist
    integer :: i
    integer, dimension(1) :: indb
    logical :: infra_node = .false.
   
    passed = .true.

    if(allocated(s_nods)) deallocate(s_nods)
    if(allocated(r_nods)) deallocate(r_nods)
    allocate(s_nods(tns))
    allocate(r_nods(nrc)) 
    
    !find the source node indexes   
    do i=1,tns
       dist=sqrt( (s_pos(i,1) - nodes(:,1))**2 + &
            (s_pos(i,2) - nodes(:,2))**2 + &
            (s_pos(i,3) - nodes(:,3))**2)
       
       indb = minloc(dist)
       
       if(my_rank_fmm==0) then
          if(dist(indb(1)) > 0.01) then
             !!There should always be a node at all
             !!source locations. If not, something is
             !!wrong
             write(*,1001) i,dist(indb)
          end if
       end if
       s_nods(i) = indb(1)
       
    end	do
    
    !find the receiver node indexes
    
    do i=1,nrc
       dist=sqrt( (rc_pos(i,1) - nodes(:,1))**2 + &
            (rc_pos(i,2) - nodes(:,2))**2 + &
            (rc_pos(i,3) - nodes(:,3))**2)
       
       indb = minloc(dist)
       
       if(my_rank_fmm==0) then
          if(dist(indb(1)) > 0.01) then
             !!There should always be a node at all
             !!source locations. If not, something is
             !!wrong
             write(*,1002) i,dist(indb)
          end if
       end if
       r_nods(i) = indb(1)

    end	do

       
1001   format(" WARNING: SOURCE ",I5," HAS BEEN MOVED ",F10.4," TO THE NEAREST NODE")
1002   format(" WARNING: RECEIVER ",I5," HAS BEEN MOVED ",F10.4," TO THE NEAREST NODE")
    
  end subroutine mget_source_nodes
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine build_rings
    implicit none
    integer :: i,j,row,ierr,div,nstart,nend,ii,tag
    integer, dimension(nnodes) :: tmap
    integer :: status(MPI_STATUS_SIZE)

    !divide the elements of fmm slave processes
    div=floor(real(nelem)/real(n_rank_fmm-1))
  
    if(my_rank_fmm == n_rank_fmm-1) then
       nstart=(my_rank_fmm-1)*div+1
       nend = nelem
    else
       nstart=(my_rank_fmm-1)*div+1
       nend=my_rank_fmm*div
    end if
    
    !allocate ring map
    if (allocated(ring_map)) deallocate(ring_map)

    allocate(ring_map(nnodes))
    
    !build ring map
    !at this position i of ring map holds the number of
    !elements that use node i
    ring_map = 0
    do i=nstart,nend
       if(use_ele(i)) then
          do j=1,4
             row=elements(i,j);
             ring_map(row)=ring_map(row)+1;
          end do
       end if
    end do
    
    !this part initializes the temporary vector tmap so that each 
    !process put their element indices in the right place within rings
    !rings is a list of elements that use each node, listed consecutively
    !by noe
    tmap = 0
    do i=2,n_rank_fmm-1
       if(my_rank_fmm .eq. i-1) then
          call MPI_SEND(tmap+ring_map,nnodes,MPI_INTEGER,i,0,FMM_COMM,ierr)
       end if
       if(my_rank_fmm .eq. i) then
          call MPI_RECV(tmap,nnodes,MPI_INTEGER,i-1,0,FMM_COMM,status,ierr)
       end if
    end do

    !sum up ring_map
    call MPI_ALLREDUCE(MPI_IN_PLACE,ring_map,nnodes,MPI_INTEGER,MPI_SUM,SCOMM_FMM,ierr)
    
    !here ring_map is changed so that ring_map(i) holds the position withing rings 
    !that marks the end of the element list for node(i).
    do i=1,nnodes
       ii=nnodes-i+1
       ring_map(ii) = sum(ring_map(1:ii))
    end do
   
    
    !allocate and populate rings
    !tmap keeps track of where we are in the rings array
    if(allocated(rings)) deallocate(rings)
    allocate(rings(ring_map(nnodes)))
    
    rings = 0
    do i=nstart,nend
       if(use_ele(i)) then
          do j=1,4
             row=elements(i,j)
             tmap(row)=tmap(row)+1
             if(row.ne.1) then
                rings(ring_map(row-1)+tmap(row))=i
             else
                rings(tmap(row))=i
             end if
          end do
       end if
    end do
 
    !combine the contributions from each element t complete rings
    call MPI_ALLREDUCE(MPI_IN_PLACE,rings,ring_map(nnodes),MPI_INTEGER,MPI_SUM,SCOMM_FMM,ierr)

  end subroutine build_rings
  !____________________________________________________________________

  !____________________________________________________________________
! compute the cross product between two vectors
  function cross_prod(u,v)
!
    implicit none
    real, dimension(3), intent(in) :: u,v
    real, dimension(3) :: cross_prod
    cross_prod(1) = u(2)*v(3) - u(3)*v(2)
    cross_prod(2) = u(3)*v(1) - u(1)*v(3)
    cross_prod(3) = u(1)*v(2) - u(2)*v(1)
  
  end function cross_prod
  !____________________________________________________________________

  !____________________________________________________________________
! compute solid angle
  function solid_angle(x1,x2,x3,x4,nd)

    implicit none
    real, dimension(3), intent(in) :: x1,x2,x3,x4
    integer, intent(in) :: nd
    real, dimension(3) :: v1,v2,v3,di
    real, dimension(3,1) :: nv1,nv2,nv3
    real :: tmp(1,1)
    real :: pi, solid_angle
 
    if(nd == 1) then
      v1(:) = cross_prod(x4(:)-x1(:), x2(:)-x1(:))
      nv1(:,1) = v1(:)/NORM2(v1)
 
      v2(:) = cross_prod(x3(:)-x1(:), x4(:)-x1(:))
      nv2(:,1) = v2(:)/NORM2(v2)

      v3(:) = cross_prod(x2(:)-x1(:), x3(:)-x1(:))
      nv3(:,1) = v3(:)/NORM2(v3)

    elseif(nd ==2) then
      v1(:) = cross_prod(x1(:)-x2(:), x4(:)-x2(:))
      nv1(:,1) = v1(:)/NORM2(v1)
 
      v2(:) = cross_prod(x4(:)-x2(:), x3(:)-x2(:))
      nv2(:,1) = v2(:)/NORM2(v2)

      v3(:) = cross_prod(x3(:)-x2(:), x1(:)-x2(:))
      nv3(:,1) = v3(:)/NORM2(v3)
    
    elseif(nd ==3) then
      v1(:) = cross_prod(x2(:)-x3(:), x4(:)-x3(:))
      nv1(:,1) = v1(:)/NORM2(v1)
 
      v2(:) = cross_prod(x4(:)-x3(:), x1(:)-x3(:))
      nv2(:,1) = v2(:)/NORM2(v2)

      v3(:) = cross_prod(x1(:)-x3(:), x2(:)-x3(:))
      nv3(:,1) = v3(:)/NORM2(v3)
    
    elseif(nd ==4) then
      v1(:) = cross_prod(x2(:)-x4(:), x1(:)-x4(:))
      nv1(:,1) = v1(:)/NORM2(v1)
 
      v2(:) = cross_prod(x3(:)-x4(:), x2(:)-x4(:))
      nv2(:,1) = v2(:)/NORM2(v2)

      v3(:) = cross_prod(x1(:)-x4(:), x3(:)-x4(:))
      nv3(:,1) = v3(:)/NORM2(v3)
   
    else
      print *,'Check nd ....'
      solid_angle = 0.0
      return
    endif
    
    tmp = ACOS(-MATMUL(TRANSPOSE(nv1(:,:)),nv2(:,:)))
    di(1) = tmp(1,1)
    tmp = ACOS(-MATMUL(TRANSPOSE(nv2(:,:)),nv3(:,:)))
    di(2) = tmp(1,1)
    tmp = ACOS(-MATMUL(TRANSPOSE(nv3(:,:)),nv1(:,:)))
    di(3) = tmp(1,1)
    
    pi = 4.0*ATAN(1.0)
    solid_angle = SUM(di)-pi

     
  end function solid_angle
  
  !____________________________________________________________________

  !____________________________________________________________________
! compute norm_2
  function NORM2(v)
    implicit none
    real, dimension(3), intent(in) :: v
    real :: norm2
    norm2 = sqrt(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
  end function NORM2
  !____________________________________________________________________

  !____________________________________________________________________
! compute travel time on face and edge
  function face_edge_tt(x1f,x2f,x3f,t1f,t2f,M)
   
    implicit none
    real, dimension(3), intent(in) :: x1f,x2f,x3f
    real, dimension(3,3), intent(in) :: M
    real, dimension(3,1) :: e12, e13, e31, e32
    real :: ax(1,1),bx(1,1),cx(1,1)
    real :: t1f,t2f, t1f2,a,b,c, c1,c2,c3, lam1,lam2
    real :: t3_1,t3_2, del_tt(1,1)
    real :: face_edge_tt
    integer :: check1, check2

    e12(:,1) = x1f(:)-x2f(:)
    e13(:,1) = x1f(:)-x3f(:)
    t1f2 = t2f-t1f
    ax = MATMUL(TRANSPOSE(e12(:,:)),MATMUL(M(:,:),e12(:,:)))
    bx = 2.0*MATMUL(TRANSPOSE(e12(:,:)),MATMUL(M(:,:),e13(:,:)))
    cx = MATMUL(TRANSPOSE(e13(:,:)),MATMUL(M(:,:),e13(:,:)))
    
    a = ax(1,1)
    b = bx(1,1)
    c = cx(1,1)
    c1 = a*(t1f2*t1f2) - a*a
    c2 = b*a-(t1f2*t1f2)*b
    c3 = c*(t1f2*t1f2) - 0.25*b*b
   
    check1 = 0
    check2 = 0
  
    if(c2*c2-4*c1*c3 >= 0 .and. c1 /= 0) then
!        if we're here then the normal to the wavefront that
!        passes through x3f is in this triangle
       lam1 = (-c2 + sqrt(c2*c2-4*c1*c3))/(2.0*c1)
       lam2 = (-c2 - sqrt(c2*c2-4*c1*c3))/(2.0*c1)
       t3_1 = 1.e30
       t3_2 = 1.e30
       if(lam1 >= 0.0 .and. lam1 <= 1.0) then
         check1 = 1
         del_tt = sqrt(MATMUL(TRANSPOSE((e13-lam1*e12)),MATMUL(M,(e13-lam1*e12))))
         t3_1 = t1f + lam1*t1f2 + del_tt(1,1)
       endif

       if(lam2 >= 0.0 .and. lam2 <= 1.0) then
         check2 = 1
         del_tt = sqrt(MATMUL(TRANSPOSE((e13-lam2*e12)),MATMUL(M,(e13-lam2*e12))))
         t3_2 = t1f + lam2*t1f2 + del_tt(1,1)
       endif

       if(check1 + check2 > 0) then
         face_edge_tt = min(t3_1,t3_2)
       else
!            if we're here then compute the minimum travel time along
!            the edges
         e31(:,1) = x3f-x1f
         e32(:,1) = x3f-x2f
         del_tt = sqrt(MATMUL(TRANSPOSE(e31(:,:)),MATMUL(M,e31(:,:))))
         t3_1 = del_tt(1,1)+t1f
         del_tt = sqrt(MATMUL(TRANSPOSE(e32(:,:)),MATMUL(M,e32(:,:))))
         t3_2 = del_tt(1,1)+t2f
         face_edge_tt = min(t3_1,t3_2)
       endif
    else
!        if we're here then compute the minimum travel time along
!        the edges
       e31(:,1) = x3f-x1f
       e32(:,1) = x3f-x2f
       del_tt = sqrt(MATMUL(TRANSPOSE(e31(:,:)),MATMUL(M,e31(:,:))))
       t3_1 = del_tt(1,1)+t1f
       del_tt = sqrt(MATMUL(TRANSPOSE(e32(:,:)),MATMUL(M,e32(:,:))))
       t3_2 = del_tt(1,1)+t2f
       face_edge_tt = min(t3_1,t3_2)
    endif

    if(face_edge_tt < t1f .or. face_edge_tt < t2f) then
       face_edge_tt = 1.e30
    end if

  end function face_edge_tt

end module forward_fmm
 
