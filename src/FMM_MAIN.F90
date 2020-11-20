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
       call nreport_fmm(67)
       
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

        call nreport_fmm(73)
     end do

     ! The baseline inversion is done. If this is a time-lapse
     ! inversion, invert the time lapse data files.     
     if (tl_ly) call time_lapse_1_fmm
          
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

  !==========================================================================================
  ! Time-Lapse Inversion
  ! Author: Piyoosh Jaysaval
  ! email : piyoosh.jaysaval@pnnl.gov
  !
  ! Descriptions: This subroutine executes the time-lapse inversion for seismic tomogrpahy
  ! Adopted from time_lapse_1 for ERT
  !==========================================================================================

  subroutine time_lapse_1_fmm
    implicit none

    !allocate the update vector if needed
    if(.not.allocated(Jaco)) then
       call alloc_sigup
    end if

    call nreport_fmm(50)

    !set the time lapse flag to true
    tl_flag = .true.

    ! Real time -> ntl=1e6
    if(rt_flag) ntl=1e6

    !!iterate over the survey files
    do i_tl=1,ntl

       call nreport_fmm(51)
       
       !reset the beta value
       if(r_last .and. conv_opt .ne. 2) then
          beta = beta/beta_red                    ! ASK??
       else
          beta=beta_s
       end if

       !set the starting model to the reference model (i.e. baseline model) if specified.
       !Otherwise the previous solution will be used as starting model.
       !(i.e. velocity is currently the previous solution)
       !NB. velocity=1/velocity^2 & refsig = 1/velocity
       if(.not.r_last) then         ! start model is the baseline sol in ref model 
          velocity = refsig*refsig
          !!velocity = 1/4.0                                          ! TMP    <------------ REMOVE      
       end if       

       !set the previous solution 
       if(allocated(prefsig)) prefsig = sqrt(velocity)       ! use last/ref velocity model as a ref model if pref
       con_flag = .false.
       iter = 0

       !get the observed data for this time step
       if(rt_flag) then
          stop 1500
          !call get_dobs_rt_fmm(i_tl)
       else
          call get_dobs_tl_fmm(i_tl)
       end if       

       !!--------- TMP ----------------
       !!====================================
       !!call send_slowness
       !!call run_forward_fmm
       !!call build_WmII
       !!call get_ttpred
       !!====================================
       !!------------------------------
       
       !check convergence
       call check_convergence           

       !write the current solution to file
       call write_veliter
       call nreport(1)

       !!these are the outer iterations for this data set       
       !if(rt_flag) con_flag = .false.
       do while(.not. con_flag)

          iter = iter+1

          ! print iteration number
          call nreport_fmm(67)  

          !compute the Jacobian          
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

          !send the updated velocity/slowness to the slave
          call send_slowness

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
          
          !write the solution to file
          call write_veliter

          if(cgmin_flag(1) .and. cgmin_flag(2)) then             
             !call sync_convergence             
             cgmin_flag(1) = .not. con_flag             
             call sync_joint             
             call sync_beta             
          end if

          ! print end iteration
          call nreport_fmm(73)
          
       enddo

       if(rt_flag) then         
          !call write_velocity_rttl(tl_dfils(1))          
       else          
          call write_velocity_tl(tlt(i_tl))          
       end if
              
    enddo

    return
    
  end subroutine time_lapse_1_fmm  
  
  
end module fmm_main
