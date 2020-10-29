module distributor
  !!This module divides the processes and sets corresponding variables 
  !!depending on command line input
  !!Yilin Fang and Tim Johnson, 12/18/2015
  use vars
  use fmm_vars
  implicit none
  integer :: nargs, iarg, stat
  integer :: my_color, my_key, my_comm, my_group, my_grank, my_commsize
  integer :: my_mcolor, my_mkey, my_mcomm, my_mgroup, my_mgrank, my_mcommsize
  integer :: my_wrank, n_wrank
  character*4 :: arg,tstr
  logical :: petsc_initialized, simulate_e4d = .true.


contains
  
  subroutine distribute
    implicit none
    
    !call print_banner

    !Get command line argument to see if we're running FMM
    !or E4D or both. 
    nargs = command_argument_count()
    simulate_fmm = .false.
    num_fmm_procs = 0
    iarg = 0
    !do while(iarg < nargs)

    iarg = iarg + 1
    call get_command_argument(iarg,arg)
    
    if( arg == '-fmm' ) then
       iarg = iarg + 1
       call get_command_argument(iarg,tstr)
       read( tstr,*,iostat=stat ) num_fmm_procs
       call dist_report(1,stat)
       if(stat == -1) then
          call dist_abort
       end if
    end if
   
    !if there are fmm processors then set the fmm flag
    !to true and make some checks
    if(num_fmm_procs > 0) then
       
       simulate_fmm = .true.
       
       !FMM requires at least two processes
       !if(num_fmm_procs < 2) then
       !   call dist_report(2,-1)
       !   call dist_abort
       !end if
       
       !If all the processes are assigned to FMM
       !then turn off the e4d flag
       if(n_rank-num_fmm_procs == 0) then
          simulate_e4d = .false.
       end if
       
       !!if there are processes left over make
       !!sure there are at least 2 for e4d
       !if(n_rank-num_fmm_procs < 2) then
       !   call dist_report(3,-1)
       !   call dist_abort
       !end if
       
       !!if there aren't enough processors then abort
       if(n_rank - num_fmm_procs < 0) then
          call dist_report(4,-1)
          call dist_abort
       end if
       
    end if
    !end do
 
    !split the processes
    master_proc_fmm = n_rank - num_fmm_procs
   
    if (simulate_fmm) then
       !divide the processes
       if( my_rank >= master_proc_fmm ) then
          !!I'm a fmm process
          my_color = 1
          my_key = my_rank - master_proc_fmm
          
       else
          !I'm an e4d process
          my_color = 0
          my_key = my_rank
       endif
       
       call MPI_COMM_SPLIT( MPI_COMM_WORLD, my_color, my_key, my_comm,ierr )
       if( my_rank >= master_proc_fmm ) then
          call MPI_COMM_DUP( my_comm, FMM_COMM, ierr )
       else
          call MPI_COMM_DUP( my_comm, E4D_COMM, ierr )
       endif
       call MPI_COMM_GROUP( my_comm, my_group, ierr )
       call MPI_COMM_RANK( my_comm, my_grank, ierr )
       call MPI_COMM_SIZE( my_comm, my_commsize, ierr )
       if(my_rank >= master_proc_fmm) then
          my_rank_fmm = my_grank
          n_rank_fmm = my_commsize
       else
          n_rank = my_commsize
       endif
       !Build a communicator for the masters
       if(my_rank == 0 .or. my_rank == master_proc_fmm) then
          my_mcolor = 2
          if(my_rank == 0) then
             my_mkey = 0
          else
             my_mkey = 1
          end if
       end if
       call MPI_COMM_SPLIT( MPI_COMM_WORLD, my_mcolor, my_mkey, my_mcomm,ierr )
       call MPI_COMM_DUP( my_mcomm, M_COMM, ierr )
       call MPI_COMM_GROUP( my_mcomm, my_mgroup, ierr )
       call MPI_COMM_RANK( my_mcomm, my_mgrank, ierr )
       call MPI_COMM_SIZE( my_mcomm, my_mcommsize, ierr )

       call dist_report(5,-1)
    else
       
       !!only running e4d
       call MPI_COMM_DUP( mpi_comm_world, E4D_COMM, ierr )
       call MPI_COMM_GROUP( E4D_COMM, my_group, ierr )
       call MPI_COMM_RANK( E4D_COMM, my_grank, ierr )
       call MPI_COMM_SIZE( E4D_COMM, my_commsize, ierr )
    end if
    
   
  end subroutine distribute
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine dist_report(which,stat)
    implicit none

    integer :: which
    integer :: stat

    select case(which)
       case(1)
          if(stat.ne.0) then
             if(my_rank==0) then
                write(*,*) " The FFM flag -fmm was specified, but there was"
                write(*,*) " an error reading the number of processes to use"
                write(*,*) " for FMM. Please specify the number of FFM "
                write(*,*) " processes to use after the -fmm flag."
                write(*,*) " Aborting ..."
             end if
             stat = -1
          end if

       case(2)
          if(my_rank==0) then
             write(*,*) " The -fmm flag was specified, but only ", num_fmm_procs
             write(*,*) " processes were specified on the command line."
             write(*,*) " FMM requires at least 2 processes."
             write(*,*) " Aborting ..."
          end if

       case(3)
          if(my_rank == 0) then
             write(*,*) " Number of processes = ",n_rank
             write(*,*) " Number of FMM processes = ",num_fmm_procs
             write(*,*) " Number of E4D processes = ",n_rank - num_fmm_procs
             write(*,*) " E4D requires at least 2 processes."
             write(*,*) " Aborting ..."
          end if

       case(4)
          if(my_rank==0) then
             write(*,*) " Number of processes = ",n_rank
             write(*,*) " Number of FMM processes requested = ",num_fmm_procs
             write(*,*) " There aren't enough processors."
             write(*,*) " Aborting ..."
          end if
          
       case(5)
          if(my_rank_fmm == 0 .and. my_rank >= master_proc_fmm) then
             write(*,*) "Running FMM with: ",n_rank_fmm," processes"
          end if
          if(my_rank==0 .and. my_rank < master_proc_fmm) then
             write(*,*) "Running E4D with: ",n_rank," processes"
          end if
      
       case default
          return

    end select

  end subroutine dist_report
  !_____________________________________________________________________________________

  !_____________________________________________________________________________________
  subroutine dist_abort
    implicit none
    
    call PetscInitialized(petsc_initialized, ierr)
    if( petsc_initialized ) then
       call PetscFinalize(ierr)
    endif
    stop
  end subroutine dist_abort
  !_____________________________________________________________________________________


  !_____________________________________________________________________________________
  subroutine print_banner
    
    !This subroutine initializes e4d.log
    implicit none
    character*8 :: date
    character*10 :: time,year,month,day,hour,minute,second
    character*5 :: zone
    integer, dimension(10) :: values
    
    call DATE_AND_TIME(date,time,zone,values)
    select case(values(2))
    case(1)
       write(month,*) 'January'
    case(2)
       write(month,*) 'February'
    case(3)
       write(month,*) 'March'
    case(4)
       write(month,*) 'April'
    case(5)
       write(month,*) 'May'
    case(6)
       write(month,*) 'June'
    case(7)
       write(month,*) 'July'
    case(8)
       write(month,*) 'August'
    case(9)
       write(month,*) 'September'
    case(10)
       write(month,*) 'October'
    case(11)
       write(month,*) 'November'
    case(12)
       write(month,*) 'December'
    end select
    
    write(year,"(I4)")   values(1)
    write(day,"(I2.2)")    values(3)
    write(hour,"(I2.2)")   values(5)
    write(minute,"(I2.2)") values(6)
    write(second,"(I2.2)") values(7)
    write(*,*)
    write(*,*)"************************ WELCOME TO E4D ************************"
    write(*,*)"Copyright © 2014, Battelle Memorial Institute"
    write(*,*) "All rights reserved."
    write(*,*)"Current date: ",trim(month)," ",trim(day),", ",trim(year)
    write(*,*)"Current time:  ",trim(hour),":",trim(minute),":",trim(second)
    write(*,"(A,I7.7,A)")" Running on ",n_rank," processing cores"
    write(*,*)"Please refer to e4d.log for further logging information ..."
    write(*,*)"****************************************************************"

  end subroutine print_banner
  !_____________________________________________________________________________________
end module distributor
