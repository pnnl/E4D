module fmm_main

  use vars
  use buildmesh
  use fmm_vars
  use master_fmm
  use master
  use slave_fmm
  use report_fmm
  use report
  use input_fmm
  use input
  use mod_con
  use obj
  use invert
  use joint_invert
  
contains

  subroutine fmm
    !this is the main routine for the fmm processes
    implicit none

    if(my_rank_fmm > 0) then
       !slave nodes for fmm enter and the slave routine
       !and exit fmm here
       call go_slave_fmm
       return
    end if

    !NOTE: Only the fmm master processes sees the rest of this code 
    !call start_timer
  
    !initialize the log file fmm.log
    call init_log_fmm
 
    !set the iteration number
    !iter = 0
 
    
    !Read the inputs ... fmm_read_inp is in the module
    call read_input_fmm

    !fmm requires at least 2 processors except in mesh build modes
    if (mode_fmm > 1 .and. n_rank < 2) then
       call nreport_fmm(59)
       call PetscFinalize(perr)
       stop
    end if

    if(mode_fmm == 1) then
      !build the mesh and exit
      if(.not. simulate_e4d) then 
         call build_tetgen4(i_flag)
      end if
      call send_command(0)
      !call get_time
      !call treport(11)
      call PetscFinalize(perr)
      stop
    end if
  
    !if the number of sources is smaller than the number of slaves 
    !print an error message and exit
    if(n_rank_fmm-1 .gt. ns) then
      call send_command_fmm(0)
      call nreport_fmm(58)
      call PetscFinalize(perr)
      stop
    end if
    !send the slaves the survey information: master_fmm
    call send_info_fmm

    !send the slaves the slowness information: master_fmm
    call send_slowness

    !build a communicator for the fmm slaves only: master_fmm
    call build_sgroupm_fmm
  
    !setup the forward runs: master_fmm
    call setup_forward_fmm
    call nreport_fmm(65)
    !call get_time_fmm
    !call treport_fmm(1)
    ! clean up and exit
   
    if(mode_fmm > 1) then
       call nreport_fmm(2)
       call run_forward_fmm
    end if
     
    if(mode_fmm==3) then
       
       call get_inv_optsII
       if(simulate_e4d) call sync_joint                 !see master
       
       ! validate joint inversion run
       call validate_jointInv_fmm
       
       if(cgmin_flag(1) .and. cgmin_flag(2)) then
       	  write(*,*) " FMM: waiting to trade with E4D"
       	  call get_other_dists
          call sync_zwts
          call sync_beta
       !!else
       !!   betalist(2) = beta
       end if
     
       !call nreport(2)                                  !see module: report
       call build_WmII                                  !see module: mod_con
       call send_J_on_off       
   end if
   
    !assemble the simulated data
    call get_ttpred
  
    if(mode_fmm == 2) then
     !output any requested travel time distributions
	
      call write_tt

      !build a synthetic survey file based on the simulated data                                  
      call fmm_build_srv     

      call send_command_fmm(0)
      call nreport_fmm(70)
      return
    endif

    !check the data fit
    
    call check_convergence                              !see module: obj
    call nreport(1)
   
    !the outer iterations start here
    iter = 0
    do while(.not. con_flag) 
       call alloc_sigup
       iter = iter+1
       call nreport(67)
       
       !compute the Jacobian
       call nreport_fmm(72)
       call make_jaco_fmm
    
       
       !instruct the slave to go into the e4d slave subroutine from
       !slave_fmm
       call send_command(-1)
       !do the inversion
       if(cgmin_flag(1)) then 
           call joint_pcgls
        else
           call pcgls
        end if
        
       !instruct slaves to leave the e4d slave subroutine
       call send_command(0) 

       !update the slowness field
    
       call update_sigma
     

       !send the updated field to the slave
       call send_slowness
   
       !write the solution to file
       call write_velocity

       !update the travel times (i.e. run fmm)
       call run_forward_fmm
  
       !get conductivity and send slowness to E4D if
       !this is a joint inversion
       if(cgmin_flag(1) .and. cgmin_flag(2)) then
       	  write(*,*) " FMM: waiting to trade with E4D"
       	  call get_other_dists
       end if
       
       !update the constraints
       call build_WmII

       !assemble the simulated data
       call get_ttpred
       
    
       !check for convergence
        call check_convergence
        call nreport(1)
        
     
        !see if we need to reduce beta
        call check_beta

        if(cgmin_flag(1) .and. cgmin_flag(2)) then
           !call sync_convergence
           cgmin_flag(1) = .not. con_flag
           call sync_joint
           call sync_beta
        !!else
        !!   betalist(2) = beta
        !!   write(*,*) betalist
        end if

        call nreport(73)
     end do

     
     

     !check to see if fresnel volume outputs are requested.
     call check_fresnel_output
    
     !print the sum of squared sensitivities if required
     if(mode_fmm == 3 .and. jprnt) then
        call make_jaco_fmm
        call record_sens
     end if     

     ! clean up and exit
     call send_command_fmm(0)
     return     

    
  end subroutine fmm

  
end module fmm_main
