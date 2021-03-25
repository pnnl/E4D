program main
!!____________________________________________________________________________________________  
!!  E4D
!!  Copyright  -A) 2014, Battelle Memorial Institute
!!  All rights reserved.
!!
!!  1 . Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any
!!  person or entity lawfully obtaining a copy of this software and associated documentation 
!!  files (hereinafter “the Software”) to redistribute and use the Software in source and 
!!  binary forms, with or without modification.  Such person or entity may use, copy, modify, 
!!  merge, publish, distribute, sublicense, and/or sell copies of the Software, and may permit 
!!  others to do so, subject to the following conditions:
!! 
!!   -Redistributions of source code must retain the above copyright notice, this list of 
!!    conditions and the following disclaimers. 
!!   -Redistributions in binary form must reproduce the above copyright notice, this list 
!!    of conditions and the following disclaimer in the documentation and/or other materials 
!!    provided with the distribution. 
!!   -Other than as used herein, neither the name Battelle Memorial Institute or Battelle may
!!    be used in any form whatsoever without the express written consent of Battelle.  
!!   -Redistributions of the software in any form, and publications based on work performed 
!!    using the software should include the following citation as a reference:
!!    [Cite published manuscript.  If this does not apply, delete this bulleted item].
!!
!!  2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
!!  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
!!  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
!!  BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!!  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!!  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!!  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
!!  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
!!  POSSIBILITY OF SUCH DAMAGE.
!!
!!  DISCLAIMER
!!  The Software was produced by Battelle under Contract No. DE-AC05-76RL01830 with the 
!!  Department of Energy.  The U.S. Government is granted for itself and others acting on its 
!!  behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, 
!!  prepare derivative works, distribute copies to the public, perform publicly and display 
!!  publicly, and to permit others to do so.  The specific term of the license can be identified
!!  by inquiry made to Battelle or DOE.  Neither the United States nor the United States 
!!  Department of Energy, nor any of their employees, makes any warranty, express or implied, 
!!  or assumes any legal liability or responsibility for the accuracy, completeness or 
!!  usefulness of any data, apparatus, product or process disclosed, or represents that its use 
!!  would not infringe privately owned rights.
!!______________________________________________________________________________________________
!!______________________________________________________________________________________________
!! Author: Tim Johnson
!!
!! email: e4d@pnnl.gov
!!
!! Purpose: This is the main program for E4D
!!     
!! Input variables: see individual subroutines
!!
!! Output variables: see individual subroutines
!!______________________________________________________________________________________________


  
  use vars                                
  use buildmesh
  use master
  use slave
  use input
  use mod_con
  use obj
  use report
  use invert
  use joint_invert
  use v_analytic
#ifdef resmode
  use dd_opt
#endif
  use fmm_main
  use fmm_vars
  use distributor
!#ifdef fseq
  use forward_seq
!#endif

  implicit none
!#include "finclude/petscsys.h"
  
 
 
  !Initialize MPI and PETSC
  call PetscInitialize(PETSC_NULL_CHARACTER,perr)

  !get my_rank and number of processors
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_rank, ierr)

  !Each process has a world rank and a local rank. 
  !Set the world rank
  my_wrank = my_rank
  n_wrank = n_rank

  !!read command line arguments and distribute the processors
  !!accordingly
  call distribute

  !!Enter and leave the FMM side of the code here
  if(simulate_fmm .and. my_wrank >= master_proc_fmm) then
     im_fmm = .true.
     call fmm
     call dist_abort
  end if
 
  !start the slave processes
   if(my_rank>0 ) then
     !Enter the slave subroutines (see SLAVE.F90)
     call go_slave

     !if we're here then the master process has commanded this
     !slave to clean up and exit
     call dist_abort
     
   end if


  !NOTE: Only the e4d master processes sees the rest of this code 
  !start_timer is located in the module: master
   call start_timer
  
  !initialize the log file e4d.log
   call init_log
 
  !set the iteration number
   iter = 0
  
  !Read the inputs ... read_inp is in the module: input
   call read_input
  

  !Populate the mesh with the average apparent conductivity if required
   if(ave_sig) call get_ave_sig

  !e4d requires at least 2 processors except in analytic and mesh build modes
   if(mode > 1 .and. n_rank < 2 .and. .not. analytic) then
     call nreport(59)
     call PetscFinalize(perr)
     stop
   end if

   if(analytic) then
     !compute the analytic solution and exit (see module: v_analytic)
     call compute_analytic
     call send_command(0)
     call PetscFinalize(perr)
     stop
   end if
  
   if(mode==0) then
     !! test read on input files specified in e4d.inp
     ! in read_inp
     call check_files
     if (cfg_flag) then
         ! in buildmesh
         call build_tetgen4(i_flag)    
     else
         ! in master
         call check_meshfiles
         if (inv_opt_flag) then
            ! in read_inp
			call get_inv_optsII
		 end if
     end if
     call testrun_report()
     
     call send_command(0)
     call PetscFinalize(perr)
     stop
   end if

   if(mode==1) then
     !build the mesh and exit
     if(.not. simulate_fmm) then 
        call build_tetgen4(i_flag)
     end if
     call send_command(0)
     call get_time
     call treport(11)
     call PetscFinalize(perr)
     stop
   end if

  !if the number of slaves is greater than the number of electrodes 
  !print an error message and exit
   if(n_rank-1 .gt. ne) then
     call send_command(0)
     call nreport(58)
     call PetscFinalize(perr)
     stop
   end if
  
  !Send the slaves some info they  need (see module: master)
   call send_info
  
  !do a multi-forward run if called for
   if(multi_forward) then
!#ifdef fseq
      call exec_mforward
!#endif
      call send_command(0)
      call PetscFinalize(perr)
      stop
   end if

  !send the conductivity vector to the slaves (see module: master)
  !if this is a complex conductivity computation then 
  !also send the complex conductivity
   call send_sigma
   if(i_flag) call send_sigmai

  !setup a mpi communication group for the slaves only (see module: master)
   call build_sgroupm
  
  !setup the forward runs: (see module: master)
   call setup_forward
   call nreport(65)
   call get_time
   call treport(1)

 
  !if this is a complex conductivity run then get the 
  !complex conductivity forward comps ready (see module: master)
   if(i_flag) call send_command(102)
   call nreport(66)
   call get_time
   call treport(1)
  
 
 
  !build the coupling matrix (command 3 for slaves)
  !if this is a complex conductivity forward run then 
  !setup the complex conductivity coupling matrix also
   call send_command(3)
   if(i_flag) call send_command(103)
 
  !execute a forward run
   call nreport(67)
   call nreport(68)
   if(i_flag) then
     call run_forward
     call nreport(71)
     call run_forwardi
     !do another forward run after the imaginary solution to include
     !the real current arising from the imaginary potential
     call run_forward
   else
     call run_forward
   end if

   if(res_flag) then
     !survey optimization subroutines
#ifdef resmode
     call data_def_opt
#endif
     call send_command(0)
     call PetscFinalize(perr)
     stop
   end if
  
 

  !Read the inverse options and build the constraint matrix
  !while slaves are solving 
   if(mode==3) then
     !if(i_flag) then
     !   invi = .true.
     !   call get_inv_optsII                           !see module: input
     !end if
     invi = .false.
     call get_inv_optsII 

     !if fmm is running sycronize with fmm_master to determine
     !what needs to be communicated for a joint inversion
     if(simulate_fmm) call sync_joint

     ! validate joint inversion run
     call validate_jointInv

     ! save the cgmin_flag values
     cgmin_flag_start = cgmin_flag     
     
     !send and/or get parmeters needed for cross gradient constraints     
     if(cgmin_flag(1) .and. cgmin_flag(2)) then
        write(*,*) " E4D: waiting to trade with FMM"
        call get_other_dists
        call sync_zwts
        call sync_beta
     !!else
     !!   betalist(1) = beta
     end if
     
     !call nreport(2)                                  !see module: report
     call build_WmII                                  !see module: mod_con
     call send_J_on_off
     call build_rrseq                                 !see module: master

   end if

  !get the forward run times from slaves and report
   call treport(0)                                     !see module: report
   call get_abtimes                                    !see module: master
   call treport(7)               
   call get_frtimes                                    !see module: master
   call treport(6)
 
  !Assemble the simulated data 
   call get_dpred 

  !if this is just a forward run then stop here
   if(mode==2) then 
     !output any requested potential distributions
     call write_pots                                  !see module: output
     
      !build a synthetic survey file based on the 
     !simulated data                                  
     call build_srv                                   !see module: output
    
     !check for Jacobian output
     if(jaco_out_opt .and. .not. im_fmm) then
        call build_rrseq
        call send_rrseq
        call mjaco
        !call print_jaco                               !see module master
        !check to see if specific sensitivity row outputs are requested.
        call check_jaco_row_output        
     end if
     
     !instruct slaves to clean up and exit
     !clean up and exit
     call send_command(0)
     call get_time
     call treport(11)
     call PetscFinalize(perr)
     call nreport(70)
     stop

   end if

  !send the Jacobian build sequence to the slaves
   call send_rrseq                                     !see module: master
  !check the data fit
   call check_convergence                              !see module: obj

  !report the convergence
   call nreport(1)
   call treport(-1)

  !get the phase if this is an sip inversion  
   if(i_flag) call get_phase
   
  !the inverse iterations start here
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
     if(ls_beta) then
        !the beta line search option is currently broken
        call beta_line_search
     else
        !do the inversion
        if(cgmin_flag(1)) then
           call joint_pcgls
        else
           call pcgls
        end if
        call treport(4)

        !update the conductivity
        call update_sigma

        !send the conductivity to the slaves
        call send_sigma
        if(i_flag) then
           call update_sigi
           call send_sigmai
           call send_command(103)
        end if

        call send_command(3)
        !output the conductivity for this iteration
        call write_sigma
        call send_command(5)
        !execute a forward run
        call run_forward
    
        
        !get and send the slowness update if this is a joint inversion
        if(cgmin_flag(1) .and. cgmin_flag(2)) then
           write(*,*) " E4D: waiting to trade with FMM"
           call get_other_dists
        end if
     	
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

        if(cgmin_flag(1) .and. cgmin_flag(2)) then
           !call sync_convergence
           cgmin_flag(1) = .not. con_flag
           call sync_joint
           call sync_beta
        !!else
        !!   betalist(1) = beta
        !!   write(*,*) betalist
        end if        
   
      end if

      call treport(-1)
      call nreport(73)

   end do
  
  !adjust the solution to the target chi-squared 
  !if the solution was updated from the starting model
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

  !invert the complex data if this is a sip inversion
   if(mode == 3 .and. i_flag) then
     call complex_inv
   end if
  
  !The baseline inversion is done. If this is a time-lapse
  !inversion, invert the time lapse data files.
   if(tl_ly) call time_lapse_1 
  
  !print the sum of squared sensitivities if required
   if(mode == 3 .and. jprnt .and. .not. i_flag) then
     call mjaco
     call record_sens
   end if

   !check to see if specific measurement sensitivity outputs are requested.
   call check_jaco_row_output
  
  !clean up and exit
   call send_command(0)
   call get_time
   call treport(11)
   call dist_abort
   

contains
  
  !____________________________________________________________________
  subroutine complex_inv
    !________________________________________________________________________
    !! Author: Tim Johnson
    !!
    !! email: tj@pnnl.gov
    !!
    !! Description: This subroutine inverts the complex conductivity data
    !!          
    !_________________________________________________________________________
    implicit none

    !set the complex inversion flag
    invi = .true.

    !set invi true on slaves
    call send_command(216)

    !read the phase inversion options
    call get_inv_optsII
    call update_sigi
    call nreport(2)
    if(allocated(J_on_off)) deallocate(J_on_off)
    call build_WmII
    !call send_J_on_off

    call send_sigmai
    iter = 0
   
    !check complex data convergence
    call check_convergence
    call nreport(1)
    call treport(-1)

    !allocate the update vector if necessary
    call alloc_sigup

    do while (.not. con_flag) 
       iter = iter+1
       call mjaco
       call get_jtimes
       call treport(3)
       call pcgls    
       call treport(4)
       call update_sigma
       call get_phase
       call send_sigmai
       call send_command(103)              
       call write_sigmai
       call send_command(5)
       call run_forwardi
       call build_WmII
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
       
       call nreport(1)
       call check_beta
       
    end do
    return
  end subroutine complex_inv
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine time_lapse_1 
    !________________________________________________________________________
    !! Author: Tim Johnson
    !!
    !! email: tj@pnnl.gov
    !!
    !! Description: This subroutine executes the time-lapse inversion
    !!          
    !_________________________________________________________________________
    implicit none
    
    !allocate the update vector if needed
    if(.not.allocated(Jaco)) then
       call alloc_sigup
    end if
    
    call treport(-1) 
    call get_jtimes
    call treport(3)
    call nreport(50)
    
    !set the time lapse flag to true
    tl_flag = .true.
    
    
    !!iterate over the survey files
    if(rt_flag) ntl=1e6
    do i_tl=1,ntl

       if(i_flag) then
           call build_WmII
		   
		   !set the complex inversion flag to false to be reset during complex inversion in complex_inv
		   invi = .false.

		   !set invi false on slaves
		   call send_command(215)

		   !read the conductivity inversion options
		   call get_inv_optsII
       end if
 
       call nreport(51)
       
       !reset the beta value
       if(r_last .and. conv_opt .ne. 2) then
          beta = beta/beta_red
       else
          beta=beta_s
       end if
             
       !set the starting model to the reference model 
       !if specified. Otherwise the previous solution
       !will be used (i.e. sigma is currently the previous
       !solution)
       if(.not.r_last) then  ! use baseline sol as start 
         sigma=refsig          
       end if
       
       !set the previous solution 
       if(allocated(prefsig)) prefsig=sigma
       con_flag = .false.
       iter = 0
       
       !get the observed data for this time step
       if(rt_flag) then
          call get_dobs_rt(i_tl)
       else
          call get_dobs_tl(i_tl)
       end if

       ! reset cgmin_flag and betalist
       cgmin_flag = cgmin_flag_start
       if(cgmin_flag(1) .and. cgmin_flag(2)) then                 
          write(*,*) " E4D: waiting to trade with FMM"          
          call get_other_dists          
          call sync_zwts          
          call sync_beta
       endif              

       !check convergence
       call check_convergence
       !write the current solution to file
       call write_sigiter
       call nreport(1)      
       
       !!these are the outer iterations for this data set
       !if(rt_flag) con_flag = .false.
       do while (.not. con_flag)
          
          iter = iter+1

          ! print iteration number
          call nreport(67)
          
          call treport(0)
          call mjaco
          call get_jtimes
          call treport(3)

          !do the inversion
          if(cgmin_flag(1)) then
             call joint_pcgls             
          else             
             call pcgls
          endif
          
          call treport(4)
          call update_sigma
          call send_sigma
          
          if(i_flag) then
            call update_sigi
            call send_sigmai
            call send_command(103)
          end if
          
          
          call send_command(3)             
          call send_command(5)
          call run_forward          

          !get and send the slowness update if this is a joint inversion
          if(cgmin_flag(1) .and. cgmin_flag(2)) then
             write(*,*) " E4D: waiting to trade with FMM"
             call get_other_dists
          end if
          
          call build_WmII
          
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
          call nreport(1)
          call check_beta                    
          call write_sigiter

          if(cgmin_flag(1) .and. cgmin_flag(2)) then             
             !call sync_convergence             
             cgmin_flag(1) = .not. con_flag             
             call sync_joint             
             call sync_beta             
          end if

          ! print end iteration
          call nreport(73)
          
       end do
               
       if(iter>0) then          
       
          call nreport(52)          
          call nreport(53)          

          !!final solution adjustment
          call sadjust(2)
          call send_sigma
       
          if(i_flag) then
            call update_sigi
            call send_sigmai
          end if
       
          call send_command(3)              
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
       
          !call map_sigma_par
          call check_convergence
          call nreport(55)
       end if
          
      !invert the complex data if this is a sip inversion
      if(i_flag) then
        call complex_inv
      end if
   
      if(rt_flag) then         
         call write_sigma_rttl(tl_dfils(1))         
      else         
         call write_sigma_tl(tlt(i_tl))         
      end if      
       
   end do   

   return
   
 end subroutine time_lapse_1
 
!____________________________________________________________________

!____________________________________________________________________
  subroutine beta_line_search
    
    !________________________________________________________________________
    !! Author: Tim Johnson
    !!
    !! email: tj@pnnl.gov
    !!
    !! Description: This subroutine conducts a line search on beta. It is 
    !!              currently not working, and is not available for use. 
    !_________________________________________________________________________
    implicit none
    
    real, dimension(nelem) :: stemp,snew
    real, dimension(nm) :: tdpd
    real, dimension(5) :: btest,x2,dnorm,mnorm,obj,Y
!!$    real, dimension(5,3) :: A
!!$    real, dimension(3,3) :: ATA,INDX,ATAI
!!$    real, dimension(3) :: ATY
    real :: cx2,cdnorm,cmnorm,cobj,min_obj,chi_last,chi_min
    integer :: i,min_ibet
    
  
    wopt = .true.

    !if(iopt==1) then
    tdpd = dpred
    stemp=sigma
    !elseif(iopt==2) then
    !   tdpd = dpredi
    !   stemp=sigi
    !end if

    
 
    cx2 = chi2
    cdnorm = phi_data
    cmnorm = phi_model
    cobj = phi_tot
    min_obj=cobj
    min_ibet=0
    wopt = .true.
    chi_min=chi2

    btest(1)= 2*beta
    btest(2)=beta
    btest(3)=beta/2
    btest(4)=beta/5
    btest(5)=beta/10

    call nreport(100)
    call treport(9)

    do i=1,5 
       beta= btest(i)
       call treport(10)
       call pcgls
       call treport(4)
       wopt = .false.
       call update_sigma
       call send_sigma
       call send_command(3)
       call send_command(5)
       call run_forward
       !call send_command(6)
       
       call get_abtimes
       call treport(7)
       call get_ksptimes
       call treport(8)
       call get_frtimes
       call treport(6)
       
       call map_sigma_par
       
       
       call get_time; etm1=etm
       call get_dpred
       call get_time; etm=etm-etm1
       call treport(2)

       call check_convergence
       call nreport(101)
  

       x2(i) = chi2
       dnorm(i) = phi_data
       mnorm(i) = phi_model
       obj(i) = phi_tot
       if(chi2 < cx2) then 
          min_ibet=i
          snew=sigma
          cx2=chi2
       end if
       
       if(i>1) then
          if(x2(i-1) .le.x2(i)) then
             goto 100
          end if
       end if
       if(con_flag) then
          goto 1000
       end if
       sigma=stemp
       dpred=tdpd
    end do

100 continue

    if(min_ibet > 0 ) then
       beta=btest(min_ibet)
       !if(iopt==1) then
          sigma=snew
       !elseif(iopt==2) then
       !   sigi=snew
       !end if
    else
       call nreport(102)
    end if    

    call send_sigma    
    call send_command(3)
    call send_command(5)
    call run_forward
    !call send_command(6)
    call build_Wm

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
    call nreport(1)
    call write_sigma
    return

1000 continue
    beta=2*beta
    call send_command(11)
    call get_jtimes
    call treport(3)
    wopt=.true.

1001 continue
    call nreport(6)

    iter=iter+1
    call pcgls
    call treport(4)

    wopt=.false.
    call update_sigma
    call send_sigma
    call send_command(3)
    call send_command(5)
    call run_forward
    !call send_command(6)

    call get_abtimes
    call treport(7)
    call get_ksptimes
    call treport(8)
    call get_frtimes
    call treport(6)

    call write_sigma
    call map_sigma_par
    
    call get_time; etm1=etm
    call get_dpred
    call get_time; etm=etm-etm1
    call treport(2)

    call check_convergence
    call nreport(1)
    if( abs(1-chi2/norm_chi2) < 0.1) return
    if(chi2<norm_chi2) then
       beta=2*beta
       goto 1001
    else
       beta=beta/2
       goto 1001
    end if
   
 
  end subroutine beta_line_search
  !____________________________________________________________________


end program main
