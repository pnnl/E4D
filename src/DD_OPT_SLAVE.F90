module dd_opt_slave
  !_________________________________________________________________
  !! Author: Tim Johnson
  !!
  !! email: tj@pnnl.gov
  !!
  !! Purpose: This module performs data-deficiency ERT survey
  !!          optimization: routines specific to slaves
  !!
  !! Revision History: Last modified 9/23/2015
  !!
  !__________________________________________________________________

  use vars
  use build_amap
  implicit none


  integer :: src_cnt, min_group_indx
  integer :: dd_gsz,dd_ngrps, dd_ad
  integer*8 :: my_tcnt
  integer*8, dimension(:), allocatable :: src_cnt_v,dd_st,dd_en,test_tcnt
  integer*8, dimension(:,:), allocatable :: best_tcnt
  logical*8, dimension(:), allocatable :: snn_flag
  real :: dd_alpha, dd_beta, dd_thresh, dr, min_group_dr, dr_sum
  real, dimension(:,:), allocatable :: epots_t, epots
  integer, dimension(:,:,:), allocatable :: best_groups
  integer, dimension(:,:), allocatable :: test_group
  real, dimension(:,:), allocatable :: best_dr
  real, dimension(:), allocatable :: test_dr, group_sums
  logical :: sens_allocated = .false.
  real, dimension(:), allocatable :: sens
  
contains
  !__________________________________________________________________
  subroutine recieve_dd_opts
    implicit none
    call MPI_BCAST(dd_alpha,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_beta,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_thresh,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_gsz,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_ngrps,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_ad,1,MPI_INTEGER,0,E4D_COMM,ierr)
    res_flag = .true.
  end subroutine recieve_dd_opts
  !__________________________________________________________________
  
  !__________________________________________________________________
  subroutine send_next_best
    implicit none
    integer :: i,ierr,tag,j
    integer, dimension(:,:), allocatable :: rank_group
    
    tag = 0
    do i=1,dd_ngrps
       call MPI_SEND(best_groups(i,:,:),dd_gsz*4,MPI_INTEGER,0,tag,E4D_COMM,ierr)
       call MPI_SEND(best_dr(:,i),dd_gsz,MPI_INTEGER,0,tag,E4D_COMM,ierr)      
    end do
    
    allocate(rank_group(dd_gsz*dd_ad,2))
    call MPI_BCAST(rank_group,dd_gsz*dd_ad*2,MPI_INTEGER,0,E4D_COMM,ierr)
    do i=1,dd_gsz*dd_ad
       if(rank_group(i,1)==my_rank) then
          do j=1,dd_gsz
             snn_flag(best_tcnt(j,rank_group(i,2)))=.false.
          end do
       end if
    end do
    deallocate(rank_group)
    
  end subroutine send_next_best
  !__________________________________________________________________
  !__________________________________________________________________
  subroutine ddopt_get_next_set
    implicit none
 
    integer :: a,b,m,n,i,j
    integer*8 :: tcnt
    !real :: min_group_dr

    if(.not. allocated(best_groups)) then
       allocate(best_groups(dd_ngrps,dd_gsz,4))
       allocate(test_group(dd_gsz,4))
       allocate(best_dr(dd_gsz,dd_ngrps))
       allocate(test_dr(dd_gsz))
       allocate(best_tcnt(dd_gsz,dd_ngrps))
       allocate(test_tcnt(dd_gsz))
       allocate(snn_flag(my_tcnt))
       snn_flag = .true.
    end if
    
    best_groups = 0
    test_group = 0
    best_dr = 0
    test_dr = 0
    best_tcnt = 0
    test_tcnt = 0
    tcnt = 0
    min_group_indx =1
    min_group_dr = 0
    src_cnt= 0
    !...................
    do a=1,tne-3
       b=a+1
       src_cnt=src_cnt+1
       test_dr = 0
       test_group = 0

       
       if(src_cnt.ge.dd_st(my_rank) .and. src_cnt.le.dd_en(my_rank)) then
          do m=b+1,tne-1
             do n=m+1,tne
                tcnt=tcnt+1
                if(snn_flag(tcnt)) then
                   call get_dr(a,b,m,n,dr,tcnt)
                   do i=1,dd_gsz   
                      if(dr > test_dr(i)) then
                         test_dr(i)=dr 
                         test_tcnt(i) = tcnt
                         test_group(i,1:4) = [a,b,m,n]
                         exit
                      end if
                   end do
                end if
             end do
          end do
          
          dr_sum = sum(test_dr)
          if(min_group_dr < dr_sum) then
             call swap_group
          end if
       end if
       
    end do
    
    !...................
    do a=1,tne-3
       b=a+2
       src_cnt=src_cnt+1
       test_dr = 0
       test_tcnt = 0
       test_group = 0
       
       if(src_cnt.ge.dd_st(my_rank) .and. src_cnt.le.dd_en(my_rank)) then
          
          do m=b+1,tne-1  
             do n=m+1,tne
                tcnt=tcnt+1
                if(snn_flag(tcnt)) then
                   call get_dr(a,b,m,n,dr,tcnt)
                   do i=1,dd_gsz   
                      if(dr > test_dr(i)) then
                         test_dr(i)=dr 
                         test_tcnt(i) = tcnt
                         test_group(i,1:4) = [a,b,m,n]
                         exit
                      end if
                   end do
                end if
             end do
          end do
          
          do m=a+1,b-1
             do n=b+1,tne
                tcnt=tcnt+1
                snn_flag(tcnt)=.false. !!interleaved measurement
                if(snn_flag(tcnt)) then
                   call get_dr(a,b,m,n,dr,tcnt)
                   do i=1,dd_gsz   
                      if(dr > test_dr(i)) then
                         test_dr(i)=dr 
                         test_tcnt(i)=tcnt
                         test_group(i,1:4) = [a,b,m,n]
                         exit
                      end if
                   end do
                end if
             end do
          end do
          
          dr_sum = sum(test_dr)
          if(min_group_dr < dr_sum) then
             call swap_group
          end if
       end if
    end do
    
    !..................
    do a=1,tne-3
       do b=a+3,tne-2
          src_cnt = src_cnt+1
          test_dr = 0
          test_tcnt = 0
          test_group = 0
          
          if(src_cnt.ge.dd_st(my_rank) .and. src_cnt.le.dd_en(my_rank)) then
             do m=b+1,tne-1
                do n=m+1,tne  
                   tcnt=tcnt+1
                   if(snn_flag(tcnt)) then
                      call get_dr(a,b,m,n,dr,tcnt)
                      do i=1,dd_gsz   
                         if(dr > test_dr(i)) then
                            test_dr(i)=dr 
                            test_tcnt(i) = tcnt
                            test_group(i,1:4) = [a,b,m,n]
                            exit
                         end if
                      end do
                   end if
                end do
             end do
             
             do m=a+1,b-1
                do n=b+1,tne  
                   tcnt=tcnt+1
                   snn_flag(tcnt)=.false. !!interleaved measurement
                   if(snn_flag(tcnt)) then
                      call get_dr(a,b,m,n,dr,tcnt)
                      do i=1,dd_gsz   
                         if(dr > test_dr(i)) then
                            test_dr(i)=dr 
                            test_tcnt(i) = tcnt
                            test_group(i,1:4) = [a,b,m,n]
                            exit
                         end if
                      end do
                   end if
                end do
             end do
             
             do m=a+1,b-2
                do n=m+1,b-1
                   tcnt=tcnt+1
                   if(snn_flag(tcnt)) then
                      call get_dr(a,b,m,n,dr,tcnt)
                      do i=1,dd_gsz   
                         if(dr > test_dr(i)) then
                            test_dr(i)=dr 
                            test_tcnt(i)=tcnt
                            test_group(i,1:4) = [a,b,m,n]
                            exit
                         end if
                      end do
                   end if
                end do
             end do
             dr_sum = sum(test_dr)
             if(min_group_dr < dr_sum) then
                call swap_group
             end if

          end if
       end do
    end do
    
    !...................
    do a=1,tne-3
       b=tne-1
       src_cnt=src_cnt+1
       test_dr = 0
       test_group = 0
       
       if(src_cnt.ge.dd_st(my_rank) .and. src_cnt.le.dd_en(my_rank)) then
          do m=a+1,b-1   
             do n=b+1,tne
                tcnt=tcnt+1
                if(snn_flag(tcnt)) then
                   call get_dr(a,b,m,n,dr,tcnt)
                   do i=1,dd_gsz   
                      if(dr > test_dr(i)) then
                         test_dr(i)=dr
                         test_tcnt(i) = tcnt
                         test_group(i,1:4) = [a,b,m,n]
                         exit
                      end if
                   end do
                end if
             end do
          end do
          
          do m=a+1,b-2    
             do n=m+1,b-1
                tcnt=tcnt+1
                snn_flag(tcnt)=.false. !!interleaved measurement
                if(snn_flag(tcnt)) then
                   call get_dr(a,b,m,n,dr,tcnt)
                   do i=1,dd_gsz   
                      if(dr > test_dr(i)) then
                         test_dr(i)=dr
                         test_tcnt(i)=tcnt
                         test_group(i,1:4) = [a,b,m,n]
                         exit
                      end if
                   end do
                end if
             end do
          end do
          dr_sum = sum(test_dr)
          if(min_group_dr < dr_sum) then
             call swap_group
          end if

       end if
       
    end do
    
    !...............
    do a=1,tne-3
       b=tne
       src_cnt=src_cnt+1
       test_dr =0
       test_group = 0
       if(src_cnt.ge.dd_st(my_rank) .and. src_cnt.le.dd_en(my_rank)) then
          do m=a+1,b-2  
             do n=m+1,b-1
                tcnt=tcnt+1
                if(snn_flag(tcnt)) then
                   call get_dr(a,b,m,n,dr,tcnt)
                   do i=1,dd_gsz   
                      if(dr > test_dr(i)) then
                         test_dr(i)=dr
                         test_tcnt(i) = tcnt
                         test_group(i,1:4) = [a,b,m,n]
                         exit
                      end if
                   end do
                end if
             end do
          end do
          
          dr_sum = sum(test_dr)
          if(min_group_dr < dr_sum) then
             call swap_group
          end if

       end if
       
    end do
  
100 continue
   
    return
  
    if(my_rank == 6) then
       do i=1,dd_ngrps
          write(*,*) i,'........................'
          do j=1,dd_gsz
             write(*,*) best_groups(i,j,:),best_dr(j,i),best_tcnt(j,i)
          end do
          write(*,*) 
       end do
    end if
  end subroutine ddopt_get_next_set
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine swap_group
    implicit none
    integer :: j

    !replace the minimum
    best_groups(min_group_indx,:,:) = test_group
    best_dr(:,min_group_indx) = test_dr
    best_tcnt(:,min_group_indx) = test_tcnt
    
    !!find the new minimum 
    min_group_indx = 1
    min_group_dr = sum(best_dr(:,min_group_indx))
    do j=2,dd_ngrps
       dr = sum(best_dr(:,j))
       if(dr<min_group_dr) then
          min_group_dr = dr
          min_group_indx = j
       end if
    end do
    
   
  end subroutine swap_group
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine ddopt_gather_poles(opt)
    implicit none
    integer :: opt
    integer :: i,j,cnt,ecnt
    real, dimension(:), allocatable :: tpol
    character*40 :: string

    if(.not.allocated(epots)) allocate(epots_t(tne,tne),epots(tne,tne))
    allocate(tpol(tne))
    epots=0
    ecnt=0
    do i=1,n_rank-1
       cnt=0
       do j=eind(i,1),eind(i,2)
          cnt=cnt+1
          ecnt=ecnt+1
          tpol = 0
          if(i==my_rank) then
             tpol=poles(e_nods,cnt)
          end if
          call MPI_BCAST(tpol,tne,MPI_REAL,i-1,SCOMM,ierr)
          epots(ecnt,1:tne) = tpol
       end do
    end do
    if(opt==1) then
       epots_t=0
       epots_t=epots
    end if
    return
    if(opt==1 .and. my_rank==1) then
       do i=1,tne-1
          do j=i+1,tne  
             write(*,*) i,j,epots_t(i,j),epots(i,j)
          end do
       end do
    end if
    deallocate(tpol)
 
  end subroutine ddopt_gather_poles
  !__________________________________________________________________
  
  !__________________________________________________________________
  subroutine ddopt_receive_poles
    implicit none
    integer :: opt
    integer :: i,j,cnt,ecnt
    real, dimension(:), allocatable :: tpol
    character*40 :: string

    if(.not.allocated(epots)) allocate(epots_t(tne,tne),epots(tne,tne))
    allocate(tpol(tne))
    epots=0
    epots_t=0
    tpol=0
    ecnt=0
    do i=1,n_rank-1
       cnt=0
       do j=eind(i,1),eind(i,2)
          cnt=cnt+1
          ecnt=ecnt+1
          call MPI_BCAST(tpol,tne,MPI_REAL,0,E4D_COMM,ierr)
          epots_t(ecnt,1:tne) = tpol
    
       end do
    end do
    
    deallocate(tpol)
 
  end subroutine ddopt_receive_poles
  !__________________________________________________________________
  

  
  !__________________________________________________________________
  subroutine get_dr(a,b,m,n,val,tct)
    implicit none
    integer*8 :: tct
    integer, intent(in) :: a,b,m,n
    integer :: i
    real :: val, Ptest,Ptrue,Pma,Pmb,Pna,Pnb,sdev,numer,denom
    real, dimension(nelem) :: G
    Pma = epots_t(a,m)
    Pmb = epots_t(b,m)
    Pna = epots_t(a,n)
    Pnb = epots_t(b,n)
    Ptrue = Pma-Pmb-Pna+Pnb
    sdev = dd_alpha*abs(Ptrue) + dd_beta
    
    !make sure this measurement is above the noise threshold
    if(abs(Ptrue)/sdev < dd_thresh) then
       val = 0
       snn_flag(tct)=.false.
       return
    end if
    
    Pma = epots(a,m)
    Pmb = epots(b,m)
    Pna = epots(a,n)
    Pnb = epots(b,n)
    Ptest = Pma-Pmb-Pna+Pnb

    if(.not.sens_allocated) then
       allocate(sens(my_tcnt))
       sens=0
       sens_allocated = .true.
    end if
 

    if(sens(tct)==0) then
       G=inv_dist(a,:)*inv_dist(m,:)+inv_dist(n,:)*inv_dist(b,:)-inv_dist(a,:)*inv_dist(n,:)-inv_dist(m,:)*inv_dist(b,:)
       G=G*evol
       sens(tct)=sqrt(dot_product(G,G))
    end if
    
    val = abs((Ptest-Ptrue)/sdev)*sens(tct)
 
  end subroutine get_dr
  !__________________________________________________________________

  !____________________________________________________________________
  subroutine build_dd_st_en
    !this routine counts the number of unique current source pairs and
    !the number of unique potenital measurements associated with each pair
    !it then divides the number of measurements to investigate evenly
    !between the slave processors
    !dd_st and dd_en are respectively the first and last source pairs
    !in the list that the corresponding slave process is responsible
    !for investigating
    
    implicit none
    integer*8 :: i,j,tcnt,ctarg,targ,sm2,sm1,d1,d2
    integer :: a,b,m,n

    if(allocated(dd_st)) deallocate(dd_st)
    if(allocated(dd_en)) deallocate(dd_en)
    
    src_cnt = 0
    tcnt = 0
    src_cnt=4*(tne-3)
    do a=1,tne-3
       do b=a+3,tne-2
          src_cnt=src_cnt+1
       end do
    end do
    allocate(src_cnt_v(src_cnt))
    src_cnt_v = 0
    
    !left end
    src_cnt=0
    do a=1,tne-3
       b=a+1
       src_cnt=src_cnt+1
       do m=b+1,tne-1
          do n=m+1,tne
             src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
             tcnt=tcnt+1
          end do
       end do
    end do
    
    do a=1,tne-3
       b=a+2
       src_cnt=src_cnt+1
       do m=b+1,tne-1
          do n=m+1,tne
             src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
             tcnt=tcnt+1
          end do
       end do
       
       do m=a+1,b-1
          do n=b+1,tne
             src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
             tcnt=tcnt+1
          end do
       end do
    end do
    
    !center
    do a=1,tne-3
       do b=a+3,tne-2
          src_cnt = src_cnt+1
          
          do m=b+1,tne-1
             do n=m+1,tne
                src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
                tcnt=tcnt+1
             end do
          end do
          
          do m=a+1,b-1
             do n=b+1,tne
                src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
                tcnt=tcnt+1
             end do
          end do
          
          do m=a+1,b-2
             do n=m+1,b-1
                src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
                tcnt=tcnt+1
             end do
          end do
          
       end do
    end do
    
    !right end
    do a=1,tne-3
       b=tne-1
       src_cnt = src_cnt+1
       
       do m=a+1,b-1
          do n=b+1,tne
             src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
             tcnt=tcnt+1
          end do
       end do
       
       do m=a+1,b-2
          do n=m+1,b-1
             src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
             tcnt=tcnt+1
          end do
       end do
    end do
    
    do a=1,tne-3
       b=tne
       src_cnt = src_cnt+1
       do m=a+1,b-2
          do n=m+1,b-1
             src_cnt_v(src_cnt) = src_cnt_v(src_cnt)+1
             tcnt=tcnt+1
          end do
       end do
    end do
    
    targ = tcnt/(n_rank-1)
    
    if(allocated(dd_st)) then
       deallocate(dd_st)
       deallocate(dd_en)
    end if
    allocate(dd_st(n_rank-1),dd_en(n_rank-1))
    
    dd_st(1)=1
    do i=1,n_rank-2
       ctarg=i*targ
       
       do j=dd_st(1),src_cnt
          sm2 = sum(src_cnt_v(1:j))
          
          if(sm2>ctarg) then
             sm1=sum(src_cnt_v(1:j-1))
             d1=abs(ctarg-sm1)
             d2=abs(ctarg-sm2)
             
             if(d2>d1) then
                dd_en(i)=j
             else
                dd_en(i)=j-1
             end if
             dd_st(i+1)=dd_en(i)+1
             
             exit
          end if
       end do
       
    end do
    
    dd_en(n_rank-1) = src_cnt
    call get_my_tcnt
    
  end subroutine build_dd_st_en
  !____________________________________________________________________


  !____________________________________________________________________
  subroutine receive_new_jinds
    implicit none
    
    call MPI_BCAST(jind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    
  end subroutine receive_new_jinds
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine rebuild_dobs_dpred_Wd
    implicit none
    integer :: i,a,b,m,n
    real :: Pma,Pmb,Pna,Pnb

    if(allocated(dobs)) deallocate(dobs)
    if(allocated(dpred)) deallocate(dpred)
    if(allocated(Wd)) deallocate(Wd)
    allocate(dobs(nm),dpred(nm),Wd(nm))
    
    do i=1,nm
       a=s_conf(i,1)
       b=s_conf(i,2)
       m=s_conf(i,3)
       n=s_conf(i,4)

       Pma = epots(a,m)
       Pmb = epots(b,m)
       Pna = epots(a,n)
       Pnb = epots(b,n)
       dpred(i) = Pma-Pmb-Pna+Pnb

       Pma = epots_t(a,m)
       Pmb = epots_t(b,m)
       Pna = epots_t(a,n)
       Pnb = epots_t(b,n)
       dobs(i) = Pma-Pmb-Pna+Pnb
       
       Wd(i) = 1/(dd_alpha*abs(dobs(i)) + dd_beta)
    end do

    !update from rank 1
    call MPI_BCAST(dobs,nm,MPI_REAL , 1, E4D_COMM, ierr )
    call MPI_BCAST(dpred,nm,MPI_REAL , 1, E4D_COMM, ierr )
    call MPI_BCAST(Wd,nm,MPI_REAL , 1, E4D_COMM, ierr )
    
  end subroutine rebuild_dobs_dpred_Wd
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_my_tcnt
    implicit none
    integer :: a,b,m,n
    integer*8 :: s_cnt
    
    my_tcnt = 0
    s_cnt = 0
    !...................
    do a=1,tne-3
       b=a+1
       s_cnt=s_cnt+1
       if(s_cnt.ge.dd_st(my_rank) .and. s_cnt.le.dd_en(my_rank)) then
          do m=b+1,tne-1
             do n=m+1,tne
                my_tcnt = my_tcnt+1               
             end do
          end do
       end if     
    end do
    
    !...................
    do a=1,tne-3
       b=a+2
       s_cnt=s_cnt+1
       if(s_cnt.ge.dd_st(my_rank) .and. s_cnt.le.dd_en(my_rank)) then
          
          do m=b+1,tne-1  
             do n=m+1,tne
                my_tcnt=my_tcnt+1         
             end do
          end do
          
          do m=a+1,b-1
             do n=b+1,tne
                my_tcnt=my_tcnt+1
             end do
          end do
          
       end if
    end do
    
    !..................
    do a=1,tne-3
       do b=a+3,tne-2
          s_cnt = s_cnt+1
          if(s_cnt.ge.dd_st(my_rank) .and. s_cnt.le.dd_en(my_rank)) then
             
             do m=b+1,tne-1
                do n=m+1,tne  
                   my_tcnt=my_tcnt+1
                end do
             end do
             
             do m=a+1,b-1
                do n=b+1,tne  
                   my_tcnt=my_tcnt+1 
                end do
             end do
             
             do m=a+1,b-2
                do n=m+1,b-1
                   my_tcnt=my_tcnt+1 
                end do
             end do
             
          end if
       end do
    end do
    
    !...................
    do a=1,tne-3
       b=tne-1
       s_cnt=s_cnt+1    
       if(s_cnt.ge.dd_st(my_rank) .and. s_cnt.le.dd_en(my_rank)) then
          
          do m=a+1,b-1   
             do n=b+1,tne
                my_tcnt=my_tcnt+1
             end do
          end do
          
          do m=a+1,b-2    
             do n=m+1,b-1
                my_tcnt=my_tcnt+1
             end do
          end do
       end if
    end do
    
    !...............
    do a=1,tne-3
       b=tne
       s_cnt=s_cnt+1
       if(s_cnt.ge.dd_st(my_rank) .and. s_cnt.le.dd_en(my_rank)) then

          do m=a+1,b-2  
             do n=m+1,b-1
                my_tcnt=my_tcnt+1
             end do
          end do
    
       end if   
    end do
  
  end subroutine get_my_tcnt
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine sbuild_inverse_dist
    !Build the inverse distance matrix for my electrodes and then 
    !broadcast 
    implicit none
    integer :: i,j
    real*8, dimension(4) :: x,y,z
    real*8 :: vol
    real, dimension(:), allocatable :: dist
    real, dimension(:,:), allocatable :: mids
    

    if(allocated(nodes)) deallocate(nodes)
    if(allocated(elements)) deallocate(elements)
    allocate(nodes(nnodes,3))
    allocate(elements(nelem,4))
    allocate(dist(nelem),mids(nelem,3))
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)
    
    allocate(inv_dist(tne,nelem),evol(nelem))
    inv_dist = 0
    evol = 0

    !get the midpoints of the elements
    do i=1,nelem
       do j=1,3
          mids(i,j) = 0.25*sum(nodes(elements(i,:),j))
       end do
       x = nodes(elements(i,:),1)
       y = nodes(elements(i,:),2)
       z = nodes(elements(i,:),3)
       call get_vol( x(1),y(1),z(1),x(2),y(2),z(2),x(3),y(3),z(3),x(4),y(4),z(4),vol)
       evol(i) = real(vol)
    end do
    
    do i=eind(my_rank,1),eind(my_rank,2)
       dist=sqrt( (mids(:,1)-e_pos(i,1))**2 + (mids(:,2)-e_pos(i,2))**2 + (mids(:,3)-e_pos(i,3))**2)
       inv_dist(i,:) = dist**(-1)
    end do
    
    call MPI_ALLREDUCE(MPI_IN_PLACE,inv_dist,tne*nelem,MPI_REAL,MPI_SUM,SCOMM,ierr)
    
    deallocate(dist)
    deallocate(mids)
  end subroutine sbuild_inverse_dist
  !____________________________________________________________________

!!$  !__________________________________________________________________
!!$  subroutine ddopt_get_high_snn
!!$    !this routine checks the signal to noise threshold of each
!!$    !measurement for which this processor is responsible. If the
!!$    !snr is below the user supplied threshold, then the flag for
!!$    !this measurement is marked as false, meaning it will not
!!$    !be consider in the optimization
!!$    
!!$    implicit none
!!$    integer*8 :: i,j,tcnt
!!$    integer :: a,b,m,n
!!$    logical :: cflag
!!$    
!!$    if(.not. allocated(dd_st)) then
!!$       call build_dd_st_en
!!$    end if
!!$    return
!!$  
!!$    allocate(snn_flag((tne-3)*(tne-2)*(tne-1)*tne/8+1))
!!$    snn_flag = .false.
!!$
!!$    src_cnt = 0
!!$    tcnt = 0
!!$   
!!$    !...................
!!$    do a=1,tne-3
!!$       b=a+1
!!$       src_cnt=src_cnt+1
!!$       cflag = .false.
!!$       !see if this source pair is mine
!!$       if(src_cnt>=dd_st(my_rank) .and. src_cnt<=dd_en(my_rank)) then
!!$          cflag = .true.
!!$       end if
!!$       do m=b+1,tne-1
!!$          do n=m+1,tne
!!$             tcnt=tcnt+1
!!$             if(cflag) then
!!$                call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$             end if
!!$          end do
!!$       end do
!!$      
!!$    end do
!!$
!!$    !...................
!!$    do a=1,tne-3
!!$       b=a+2
!!$       src_cnt=src_cnt+1
!!$       cflag = .false.
!!$
!!$       if(src_cnt>=dd_st(my_rank) .and. src_cnt<=dd_en(my_rank)) then
!!$          cflag = .true.
!!$       end if
!!$
!!$       do m=b+1,tne-1
!!$          do n=m+1,tne
!!$             tcnt=tcnt+1
!!$             if(cflag) then
!!$                call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$             end if
!!$          end do
!!$       end do
!!$
!!$       do m=a+1,b-1
!!$          do n=b+1,tne
!!$             tcnt=tcnt+1
!!$             if(cflag) then
!!$                call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$             end if
!!$          end do
!!$       end do
!!$ 
!!$    end do
!!$
!!$    !..................
!!$    do a=1,tne-3
!!$       do b=a+3,tne-2
!!$          src_cnt = src_cnt+1
!!$          cflag = .false.
!!$
!!$          if(src_cnt>=dd_st(my_rank) .and. src_cnt<=dd_en(my_rank)) then
!!$             cflag = .true.
!!$          end if
!!$
!!$          do m=b+1,tne-1
!!$             do n=m+1,tne
!!$                tcnt=tcnt+1
!!$                if(cflag) then
!!$                   call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$                end if
!!$             end do
!!$          end do
!!$
!!$          do m=a+1,b-1
!!$             do n=b+1,tne
!!$                tcnt=tcnt+1
!!$                if(cflag) then
!!$                   call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$                end if
!!$             end do
!!$          end do
!!$
!!$          do m=a+1,b-2
!!$             do n=m+1,b-1
!!$                tcnt=tcnt+1
!!$                if(cflag) then
!!$                   call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$                end if
!!$             end do
!!$          end do
!!$
!!$       end do
!!$    end do
!!$
!!$    !...................
!!$    do a=1,tne-3
!!$       b=tne-1
!!$       src_cnt=src_cnt+1
!!$       cflag = .false.
!!$       if(src_cnt>=dd_st(my_rank) .and. src_cnt<=dd_en(my_rank)) then
!!$          cflag = .true.
!!$       end if
!!$
!!$       do m=a+1,b-1
!!$          do n=b+1,tne
!!$             tcnt=tcnt+1
!!$             if(cflag) then
!!$                call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$             end if
!!$          end do
!!$       end do
!!$       
!!$       do m=a+1,b-2
!!$          do n=m+1,b-1
!!$             tcnt=tcnt+1
!!$             if(cflag) then
!!$                call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$             end if
!!$          end do
!!$       end do
!!$      
!!$    end do
!!$    
!!$    !...............
!!$    do a=1,tne-3
!!$       b=tne
!!$       src_cnt=src_cnt+1
!!$       cflag = .false.
!!$       if(src_cnt>=dd_st(my_rank) .and. src_cnt<=dd_en(my_rank)) then
!!$          cflag = .true.
!!$       end if
!!$
!!$       do m=a+1,b-2
!!$          do n=m+1,b-1
!!$             tcnt=tcnt+1
!!$             if(cflag) then
!!$                call check_thresh(a,b,m,n,snn_flag(tcnt))
!!$             end if
!!$          end do
!!$       end do
!!$       
!!$    end do
!!$    return
!!$    j=0
!!$    do i=1,tcnt
!!$       if(snn_flag(i)) j=j+1
!!$    end do
!!$    write(*,*) my_rank,tcnt,j
!!$  end subroutine ddopt_get_high_snn
!!$  !__________________________________________________________________
!!$
!!$  !__________________________________________________________________
!!$  subroutine check_thresh(a,b,m,n,chk)
!!$    implicit none
!!$ 
!!$    integer :: a,b,m,n
!!$    logical :: chk
!!$    real :: v_over_i
!!$    real :: Pma,Pmb,Pna,Pnb
!!$    real :: sdev
!!$  
!!$    Pma = epots_t(a,m)
!!$    Pmb = epots_t(b,m)
!!$    Pna = epots_t(a,n)
!!$    Pnb = epots_t(b,n)
!!$    v_over_i = abs((Pma-Pmb)-(Pna-Pnb))
!!$    sdev = dd_alpha*v_over_i + dd_beta
!!$    
!!$    chk = .false.
!!$    if(v_over_i/sdev >= dd_thresh) then
!!$       chk = .true.
!!$    end if
!!$    
!!$    !write(*,*) my_rank
!!$    !write(*,*) v_over_i,sdev,dd_thresh,chk
!!$  end subroutine check_thresh
!!$  !__________________________________________________________________
end module dd_opt_slave
