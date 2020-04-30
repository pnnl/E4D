module dd_opt
  !_________________________________________________________________
  !! Author: Tim Johnson
  !!
  !! email: tj@pnnl.gov
  !!
  !! Purpose: This module performs data-deficiency ERT survey
  !!          optimization
  !!
  !! Revision History: Last modified 9/23/2015
  !!
  !__________________________________________________________________

  use vars
  use master
  use slave
  use obj
  use invert
  implicit none

  logical :: poles_provided =.false.
  integer :: dd_gsize, riter                            !number of data points per source
  integer :: dd_ngroups                            !number of sources to report
  integer :: dd_add                              !number of source to add per update
  integer :: npm                                 !number of pole measurements in pole file
  integer :: cnt                                 !for error reporting
  real :: dd_alph                                !noise model alpha
  real :: dd_bet                                 !noise model beta
  real :: dd_thrsh                               !threshold for accept
  real :: dd_sig_start
  

  !integer*8, dimension(:), allocatable :: tcnt_keep
  integer, dimension(:,:), allocatable :: add_srv
  character*40 :: dd_sigfile, pole_file
 
 
contains
  !__________________________________________________________________
  subroutine data_def_opt
    implicit none
    logical :: stat
    integer :: i

    open(25,file='dr_convergence.txt',status='replace',action='write')
    write(25,*) "Survey Number, Redundancy Ratio"
    close(25)
    !! read the optimization options
    call read_dd_opts(stat)
    if(.not.stat) return
    
    !! send options to slaves
    call send_dd_opts
    
    !!send the elements and nodes again and instructe the slaves
    !!to build their inverse distance matrices and volume vectors
    call build_inv_dist
    
    !!send the true pole solutions if they are provided,
    !!otherwise instruct the slaves to gather them from
    !!the forward run.
    if(poles_provided) then
       call send_poles
    else
       call send_command(1001)
    end if

    !!instruct the slaves to determine their abmn assignments
    call send_command(1002)  
   
    !!send the starting sigma value and instruct the slave
    !!to compute the poles
    !sigma=dd_sig_start
    sigfile = dd_sigfile
    call read_conductivity
    call send_sigma
    call send_command(3)
    call send_command(5)
    call nreport(67)
    call nreport(68)
    call run_forward
 

    riter=0
101 continue
    iter = 0
    riter=riter+1
    !!instruct the slaves to gather the estimated pole solutions
    call send_command(1003)
    
    !!instruct the slave to find least redundant measurements
    call send_command(1004)
   
    !!gather the optimal measurements from each slave
    call get_next_best
 
    !!send the next set of measurements
    call build_send_new_survey
    
    !!update the jacobian build assignments
    call build_send_jinds
 
    !!instruct the slaves to rebuild the dobs,dpred, and Wm vectors
    call update_data_vecs

    !write the new model and survey
    call write_survey(riter)
    
    !!update the inverted model
    invi = .false.
    call get_inv_optsII 
    call nreport(2)                                  !see module: report
    call build_WmII                                  !see module: mod_con
    call send_J_on_off
    call build_rrseq

    !send the Jacobian build sequence to the slaves
    call send_rrseq                                     !see module: master
    !check the data fit
    call check_convergence                              !see module: obj
  
    !report the convergence
    call nreport(1)
    call treport(-1)
 
    !update the model
    call dd_invert

    !write the solution for this opt step
    call write_sigma_dd(riter)
 goto 101  
  
  end subroutine data_def_opt
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine dd_invert
    implicit none
    
  
    do while (.not. con_flag)
       iter = iter+1
       call nreport(67)
       
       !instruct slaves to build the jacobian matrix
       call treport(0)
       call nreport(72)
       call mjaco
       call get_jtimes
       call treport(3)
       
       !allocate the update vector if necessary
       call alloc_sigup
       
       !build the update vector
       !do the inversion
       call pcgls    
       call treport(4)
       !update the conductivity
       call update_sigma
       !send the conductivity to the slaves
       call send_sigma
       if(i_flag) then
          call update_sigi
          call send_sigmai
       end if
       call send_command(3)
       call write_sigma
       
       !output the conductivity for this iteration
       !call write_sigma
       call send_command(5)
       !execute a forward run
       call run_forward 
       !build the new constraint equations
       call build_WmII
       call get_abtimes
       call treport(7)
       call get_frtimes
       call treport(6)
       
       call get_time; etm1=etm
       !assemble the simulated data
       call get_dpred
       call get_time; etm=etm-etm1
       call treport(2)
       !check for convergence
       call check_convergence
       
       call nreport(1)
       !see if we need to reduce beta
       call check_beta
       
       call treport(-1)
       call nreport(73)
       
    end do
    !adjust the solution to the target chi-squared 
    !if the solution was updated from the starting model
    !!!!!!!this part is not adjusting correctly so return here
    !!!!!!!for now TCJ 11/21/17
    return

   if(iter>0) then
     iter=iter+1 
     call nreport(4)
     call nreport(56)
     call sadjust(1) 
     call send_sigma
     if(i_flag) then
        call update_sigi
        call send_sigmai
     end if
     call send_command(3)
     call write_sigma
     call send_command(5)
     call run_forward 
     
     call get_abtimes
     call treport(7)
     call get_ksptimes
     call treport(8)
     call get_frtimes
     call treport(6)
  
     call get_time; etm1=etm
     call get_dpred
     call get_time; etm=etm-etm1
     call treport(2)
  
 
     call check_convergence
     call nreport(55)
  end if
  
  end subroutine dd_invert
  !__________________________________________________________________
  !__________________________________________________________________
  subroutine update_data_vecs
    implicit none
    integer :: i
    !instruct the slaves to rebuild the observed and predicted 
    !data vectors and Wd
    call send_command(1007)
    
    !receive the updates from slave rank 1
    if(allocated(dobs)) deallocate(dobs)
    if(allocated(dpred)) deallocate(dpred)
    if(allocated(Wd)) deallocate(Wd)
    allocate(dobs(nm),dpred(nm),Wd(nm))

    call MPI_BCAST(dobs,nm,MPI_REAL , 1, E4D_COMM, ierr )
    call MPI_BCAST(dpred,nm,MPI_REAL , 1, E4D_COMM, ierr )
    call MPI_BCAST(Wd,nm,MPI_REAL , 1, E4D_COMM, ierr )
   
  end subroutine update_data_vecs
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine build_send_new_survey
    implicit none 
    integer :: i,cnt
    integer, dimension(nm,4) :: tsurv
    logical :: st
 
    tsurv = s_conf
    deallocate(s_conf)
    if(nm+dd_gsize*dd_add<tne) then
       !there are more processors than measurements, so fill the survey
       !with redundant measurements
       call check_ddinp(14,st)
       allocate(s_conf(tne,4))
       cnt=0
       do while(nm+cnt<tne)
          do i=1,dd_gsize*dd_add
             cnt=cnt+1
             s_conf(nm+cnt,1:4) =  add_srv(i,1:4)
             if(nm+cnt.eq.tne) then
                exit
             end if
          end do
       end do
       nm=nm+cnt
    else
       allocate(s_conf(nm+dd_gsize*dd_add,4))
       s_conf(1:nm,1:4) = tsurv
       do i=1,dd_gsize*dd_add
          s_conf(nm+i,1:4) = add_srv(i,1:4)
       end do
       nm=nm+dd_gsize*dd_add
    end if
    call send_info

  end subroutine build_send_new_survey
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine build_send_jinds
    implicit none
    integer :: neven,nextra,ce, i

    call send_command(1006)
    !!Build jaco assignments
    jind = 0
    if(nm < (n_rank-1)) then
       do i=1,nm
          jind(i,1:2)=i
       end do
    else
       neven = nm/(n_rank-1)
       nextra = nm - neven*(n_rank-1)
       
       jind=0
       if(neven > 0) then
          ce = 1
          do i=1,n_rank-1-nextra
             jind(i,1) = ce
             jind(i,2) = ce+neven-1
             ce = ce+neven
          end do
          
          if(nextra > 0) then
             do i=n_rank-nextra,n_rank-1    
                jind(i,1) = ce
                jind(i,2) = ce+neven
                ce=ce+neven+1
             end do
          end if
       end if
    end if
   
    !!send assignments
    call MPI_BCAST(jind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    call nreport(64)
  end subroutine build_send_jinds
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine get_next_best
    implicit none
 
    integer, dimension(:,:), allocatable :: test_grp
    integer, dimension(:,:), allocatable :: rank_group
    real, dimension(:), allocatable :: test_dr,dr_sums
    real :: dr_sum,dr_min
    integer :: dr_mindex
    integer :: r,i,j
    integer :: tag, ierr , status(MPI_STATUS_SIZE)

    call send_command(1005)
    if (.not. allocated(add_srv)) then
       allocate(add_srv(dd_gsize*dd_add,4))
    end if
    allocate(rank_group(dd_gsize*dd_add,2))
    allocate(test_grp(dd_gsize,4))
    allocate(test_dr(dd_gsize))
    allocate(dr_sums(dd_add))
    
    add_srv = 0
    test_grp = 0
    test_dr = 0
    tag = 0
    dr_mindex = 1
    dr_min = 0
    dr_sums = 0
    rank_group = 0
    
    do r = 1,n_rank -1
       do i=1,dd_ngroups
          call MPI_RECV(test_grp,dd_gsize*4, MPI_INTEGER,r,tag, E4D_COMM, status, ierr)
          call MPI_RECV(test_dr,dd_gsize, MPI_REAL,r,tag,E4D_COMM, status, ierr)
              
          dr_sum = sum(test_dr)
       
          if(dr_sum>dr_min) then
             add_srv((dr_mindex-1)*dd_gsize+1:dr_mindex*dd_gsize,1:4)=test_grp
             rank_group((dr_mindex-1)*dd_gsize+1:dr_mindex*dd_gsize,1)=r
             rank_group((dr_mindex-1)*dd_gsize+1:dr_mindex*dd_gsize,2)=i
             dr_sums(dr_mindex) = dr_sum

             dr_min=dr_sums(1)
             dr_mindex = 1
             do j=1,dd_add
                if(dr_sums(j)<dr_min) then
                   dr_min=dr_sums(j)
                   dr_mindex = j
                end if
             end do
          end if
       end do
    end do
    open(25,file='dr_convergence.txt',status='old',action='write',position='append')
    write(25,*) riter,sum(sqrt(dr_sums))/(dd_gsize*dd_ngroups)
    close(25)

    call MPI_BCAST(rank_group,dd_gsize*dd_add*2,MPI_INTEGER,0,E4D_COMM,ierr)
    deallocate(test_grp,test_dr,dr_sums,rank_group)

  end subroutine get_next_best
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine read_dd_opts(stat)
    implicit none
    integer :: ist
    logical :: exst, stat
    call check_ddinp(0,stat); if(.not.stat) return
    if(stat) then
       open(13,file='e4d_ddopt.inp',status='old',action='read')
 
       read(13,*,IOSTAT=ist) dd_alph, dd_bet 
       if(ist .ne. 0) then
          call check_ddinp(1,stat); 
          return
       end if
       if(dd_alph<=0 .or. dd_bet<=0) then
          call check_ddinp(2,stat) 
          return
       end if

       read(13,*,IOSTAT=ist) dd_thrsh
       if(ist .ne. 0) then
          call check_ddinp(3,stat)
          return
       end if
       if(dd_thrsh<=0) then
          call check_ddinp(4,stat)
          return
       end if
    
       read(13,*,IOSTAT=ist) dd_sigfile
       if(ist .ne. 0) then
          call check_ddinp(5,stat)
          return
       end if
       !if(dd_sig_start<=0) then
       !   call check_ddinp(6,stat)
       !   return
       !end if
       
       read(13,*) dd_gsize, dd_ngroups, dd_add
       if(ist .ne. 0) then
          call check_ddinp(7,stat)
          return
       end if
         
       read(13,*,IOSTAT=ist) pole_file
       if(ist .ne. 0) then
          call check_ddinp(8,stat)
          return
       end if
       call check_pole_file(ist)
       if(ist .ne. 0) then
          stat = .false.
          return
       end if
       
       close(13)
    end if
    
  end subroutine read_dd_opts
  !__________________________________________________________________

  !__________________________________________________________________
  subroutine send_dd_opts
    implicit none
    call send_command(1000)
    call MPI_BCAST(dd_alph,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_bet,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_thrsh,1,MPI_REAL,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_gsize,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_ngroups,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(dd_add,1,MPI_INTEGER,0,E4D_COMM,ierr)

  end subroutine send_dd_opts
  !__________________________________________________________________


  !__________________________________________________________________
  subroutine check_ddinp(opt,st)
    implicit none
    integer :: opt
    logical :: exst,st

    select case(opt)
       case(0)
          inquire(file='e4d_ddopt.inp',exist=st)
          if(.not.st) then
             open(51,file='e4d.log',status='old',action='write',position='append')
             write(51,*) "Cannot find e4d_ddopt.inp."
             write(51,*) "aborting ..."
             close(51)
             
             write(*,*) "Cannot find e4d_ddoppt.inp."
             write(*,*) "aborting ..."
             st=.false. 
          end if
          return

       case(1)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading dd_alph or dd_bet on line 1"
          write(51,*) "of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)
           
          write(*,*) "There was a problem reading dd_alph or dd_bet on line 1"
          write(*,*) "of e4d_ddopt.inp"
          write(*,*) "aborting ..."
          st=.false. 
          return

       case(2)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "dd_alph and dd_bet must be greater than 0"
          write(51,*) "you entered",dd_alph,dd_bet," on line 1 of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)
          write(*,*) "dd_alph and dd_bet must be greater than 0"
          write(*,*) "you entered",dd_alph,dd_bet," on line 1 of e4d_ddopt.inp"
          write(*,*) "aborting ..."
          st=.false.
          return
        
       case(3)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading dd_thrsh on line 2"
          write(51,*) "of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)
           
          write(*,*) "There was a problem reading dd_thrsh on line 1"
          write(*,*) "of e4d_ddopt.inp"
          write(*,*) "aborting ..."
          st=.false. 
          return

       case(4)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "dd_thrsh must be greater than 0"
          write(51,*) "you entered",dd_thrsh," on line 2 of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)
          write(*,*) "dd_thrsh must be greater than 0"
          write(*,*) "you entered",dd_thrsh," on line 1 of e4d_ddopt.inp"
          write(*,*) "aborting ..."
          st=.false.
          return

       case(5)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading dd_sigfile on line 3"
          write(51,*) "of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)
           
          write(*,*) "There was a problem reading dd_sigfile on line 3"
          write(*,*) "of e4d_ddopt.inp"
          write(*,*) "aborting ..."
          st=.false. 
          return

       case(6)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "dd_sig_start must be greater than 0"
          write(51,*) "you entered",dd_sig_start," on line 3 of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)
          write(*,*) "dd_sig_start must be greater than 0"
          write(*,*) "you entered",dd_sig_start," on line 3 of e4d_ddopt.inp"
          write(*,*) "aborting ..."
          st=.false.
          return

       case(7)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "there was a problem reading the options on line 4 of"
          write(51,*) "of e4d_ddopt.inp"
          write(51,*) "aborting ..."
          close(51)

       case(8)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading the pole file name on line 5"
          write(51,*) "of e4d_ddopt.inp"
          write(51,*) "Survey will be optimzed for conductivity specified on line 4"
          write(51,*) "of e4d.inp."
          close(51)

       case(9)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading the number of pole measurements"
          write(51,*) "in the pole file specified in e4d_ddopt.inp, which is: ",trim(pole_file)
          write(51,*) "Instead of using the pole measurment in: ",trim(pole_file)," ,"
          write(51,*) "the survey will be optimzed for conductivity specified on line 4"
          write(51,*) "of e4d.inp."
          close(51)

       case(10)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "The number of pole measurements specified in ",trim(pole_file)
          write(51,*) "is ",npm," . The number of electrodes specified in ",trim(efile)
          write(51,*) "is ",ne, ". The number of pole measurements should be (ne*(ne-1))/2,"
          write(51,*) "which is: ", ne*(ne-1)/2
          write(51,*) "aborting ..."
          close(51)
          
      case(11)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading pole measurement ",cnt
          write(51,*) "in the pole file ",trim(pole_file)
          write(51,*) "aborting ..."
          close(51)

       case(12)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "The current or potential electrode for pole ",cnt
          write(51,*) "in the pole file ",trim(pole_file)," is listed out of order."
          write(51,*) "If (a) is the current electrode, (m) is the potential electrode, "
          write(51,*) "and ne is the number of electrodes, then the poles potentials "
          write(51,*) "should be listed in the order specified by the following "
          write(51,*) "psuedocode: "
          write(51,*) "count = 0"
          write(51,*) "for a = 1,ne-1 "
          write(51,*) "  for m=a+1,ne "
          write(51,*) "    count=count+1 "
          write(51,*) "    write to file: count,a,m,potential "
          write(51,*) "  end"
          write(51,*) "end"
          write(51,*) "aborting ..."
          close(51)

       case(13)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) "______________________________________________________"
          write(51,*) "Using pole measurements specified in ",trim(pole_file)
          write(51,*) "for optimization"
          close(51)

       case(14)
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) "WARNING: There are more electrodes than meausurements"                                         
          write(51,*) "Adding redundant measurements to equal the number"
          write(51,*) "of electrodes for this iteration"
          close(51)
          write(*,*)
          write(*,*) "WARNING: There are more electrodes than meausurements"                                         
          write(*,*) "Adding redundant measurements to equal the number"
          write(*,*) "of electrodes for this iteration"
         
          
          
    case DEFAULT
    end select
    
  end subroutine check_ddinp
  !__________________________________________________________________

  !_________________________________________________________________________________
  subroutine write_survey(it)
    implicit none
    integer :: it,i,j
    character*40 :: sfile
    
    write(sfile,'(A,I0,A)') 'opt_survey_',it,'.srv'
    open(25,file=trim(sfile),status='replace',action='write')
    write(25,'(I0)') ne
    do i=1,ne
       write(25,'(I5,3f20.3,I5)') i,e_pos(i,1)+xorig,e_pos(i,2)+yorig, &
                                    e_pos(i,3)+zorig,int(e_pos(i,4))
    end do
    write(25,*)
   
    write(25,'(I0)') nm
   
    do i=1,nm
       write(25,'(I10,4I6,E15.5,E15.5)') i,s_conf(i,1:4),dobs(i),1/Wd(i)
    end do   
    close(25)

!!$    write(sfile,'(A,I0,A)') 'opt_survey_',it,'.sig'
!!$    open(25,file=trim(sfile),status='replace',action='write')
!!$    write(25,*) nelem,1
!!$    do i=1,nelem
!!$       write(25,*) sigma(i)
!!$    end do
!!$    close(25)
    
  end subroutine write_survey
  !_________________________________________________________________________________

  !_________________________________________________________________________________
  subroutine write_sigma_dd(it)
    implicit none
    integer :: it,i,j
    character*40 :: sfile
    
    write(sfile,'(A,I0,A)') 'opt_survey_',it,'.sig'
    open(25,file=trim(sfile),status='replace',action='write')
    write(25,*) nelem,1
    do i=1,nelem
       write(25,*) sigma(i)
    end do
    close(25)
  end subroutine write_sigma_dd
  !_________________________________________________________________________________

  !_________________________________________________________________________________
  subroutine check_pole_file(st)
    implicit none
    integer, intent(out) :: st
    logical :: ext
    integer :: ist
    integer :: i, crw, ic, ip, a, m
    real :: tmp_pole

    poles_provided = .false.
    inquire(file=pole_file,exist=ext)
    
    if(.not.ext) then
       st = 0
       return
    end if

    open(21,file=pole_file,status='old',action='read')
    read(21,*,IOSTAT=ist) npm
    if(ist .ne. 0) then
       call check_ddinp(9,ext)
       st=0
       return
    end if
    
    if(npm .ne. ne*(ne-1)/2) then
       call check_ddinp(10,ext)
       st = -1
       return
    end if
    
    if(allocated(poles)) deallocate(poles)
    allocate(poles(1,npm))
    
    do a=1,ne-1
       do m=a+1, ne
          cnt=cnt+1
          read(21,*, IOSTAT=ist) crw,ic,ip,tmp_pole
          if(ist .ne. 0) then
             call check_ddinp(11,ext)
             st=-1
             return
          end if
          if(ic .ne. a .or. ip .ne. m) then
             call check_ddinp(12,ext)
             st=-1
             return
          else
             poles(1,cnt)=tmp_pole
          end if
       end do
    end do
    close(21)
    call check_ddinp(13,ext)
    poles_provided = .true.
  end subroutine check_pole_file
  !_________________________________________________________________________________

  !_________________________________________________________________________________
  subroutine send_poles
    implicit none
    real, dimension(:), allocatable :: tpol
    integer :: i,j,a,m,pcnt
    
    call send_command(1008)
    
    if(allocated(tpol)) deallocate(tpol)
    allocate(tpol(ne))

    do i=1,n_rank -1
       do j=eind(i,1),eind(i,2)
          tpol = 0
          pcnt = 0
          do a=1,ne-1
             do m=a+1,ne
                pcnt=pcnt+1
                if(a.eq.j) then
                   tpol(m) = poles(1,pcnt)
                end if
                if(m.eq.j) then
                   tpol(a) = poles(1,pcnt)
                end if
             end do
          end do
        
          call MPI_BCAST(tpol,tne,MPI_REAL,0,E4D_COMM,ierr)
       end do
    end do
    

  end subroutine send_poles
  !_________________________________________________________________________________
end module dd_opt
