module report

  use obj
  use vars
  use mod_con
  use input
  use buildmesh
 
 
  contains

    !____________________________________________________________________________________
    subroutine init_log
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
      !write(*,"(A,I7.7,A)")" (1 master process, and ",n_rank-1," slave processes)"
      write(*,*)"Please refer to e4d.log for further logging information ..."
      write(*,*)"****************************************************************"

      open(51,file='e4d.log',status='replace',action='write',IOSTAT=ios)
      if(ios.ne.0) then
         write(*,*)
         write(*,*) " THERE WAS A PROBLEM CREATING THE LOG FILE e4d.log"
         write(*,*) " Aborting ..."
         call crash_exit
         return
      end if
      write(51,*)"************************ WELCOME TO E4D ************************"
      write(51,*) "Copyright © 2014, Battelle Memorial Institute"
      write(51,*) "All rights reserved."
      write(51,*)"Current date: ",trim(month)," ",trim(day),", ",trim(year)
      write(51,*)"Current time:  ",trim(hour),":",trim(minute),":",trim(second)
      write(51,"(A,I7.7,A)")" Running on ",n_rank," processing cores"
      !write(51,"(A,I7.7,A)")" (1 master process, and ",n_rank-1," slave processes)"
      write(51,*)"****************************************************************"

      close(51)

      

    end subroutine init_log
    !____________________________________________________________________________________

    !____________________________________________________________________________________
    subroutine nreport(tag)
      implicit none
      integer :: tag,i,npre,nchr,i1,i2,i3
      character*20 :: fstring

      if(mode==1) return;
      if(tag==1) then
         if(im_fmm) then
            open(67,file='fmm.log',status='old',position='append')
         else
             open(67,file='e4d.log',status='old',position='append')
         end if
         write(67,*)
         if(iter==0) then
            write(67,*) "********** CONVERGENCE STATISTICS AT STARTING MODEL *******************"
         else
            write(67,"(A,I3.3,A)") " ********** CONVERGENCE STATISTICS AFTER INVERSE UPDATE # ",iter," **********"
         end if

         if (cgmin_flag(1)) then ! joint inversion
            write(67,"(5A15)") "    Phi_dat    ","    Phi_Mod    "," Phi_Mod/Beta  ","     Phi_CG    ","    Phi_Tot    "            
            write(67,"(5g15.5)") phi_data,phi_model,phi_model/beta,phi_cg,phi_tot            
         else ! separate inversion
            write(67,"(4A15)") "    Phi_dat    ","    Phi_Mod    "," Phi_Mod/Beta  ","    Phi_Tot    "            
            write(67,"(4g15.5)") phi_data,phi_model,phi_model/beta,phi_tot            
         endif
         

         write(67,*)
         !if(sum(imod_vec)>0) then
         write(67,"(A,I10.10)")" Total Number of Constraint Eqs: ",sum(inum_vec)
         write(67,"(A6,A18,A21,A25)") "Block","# Constraints","% of Total Error","Error per Constraint"
         do i=1,inrb(1)
            if(imod_vec(i)==0) then
               write(67,"(I6,I18,g21.5,g25.5)") i,inum_vec(i),0.0,0.0
            else
               write(67,"(I6,I18,g21.5,g25.5)") i,inum_vec(i),imod_vec(i)/sum(imod_vec) * 100, &
                    imod_vec(i)/inum_vec(i)
            end if
         end do
         write(67,*)
         !end if
         write(67,"(A,g10.4,A,g10.5)") " Chi2 is currently ",chi2,"   Target value is ",norm_chi2    
         write(67,"(A,g10.4)") " Mean error is ",errm
         write(67,"(A,g10.4)") " RMS error is ",sqrt(chi2) !erms
         write(67,"(A,I7.7)") " Number of data culled is ",ncull
         if(cull_flag==1) then
            write(67,*) "Prior to data culling ..."
            !write(67,"(A,g15.5)") "   Phi_dat is ",phi_data0
            write(67,"(A,g10.4)") "    Chi2 is ", chi20
            write(67,"(A,g10.4)") "    Mean error is ",errm0
            write(67,"(A,g10.4)") "    RMS error is ",sqrt(chi20) !erms0
         end if
         
         write(67,*)"***********************************************************************"
         close(67)
      end if

      ! 2
      if(tag==2) then
         
         write(*,*) " E4D: EXECUTING FORWARD RUN ..."
         !open(67,file='e4d.log',status='old',position='append')
         !write(67,*) "Master building model weighting matrix"
         !close(67)

      end if


      ! 3
      if(tag==3) then
         open(67,file='e4d.log',status='old',action='write',position='append')
         write(67,"(A12,I8,A19,I8,A9)") " Estimating ",npar," parameters out of ",nelem," elements"
         close(67)
      end if

      ! 4
      if(tag==4) then
         open(67,file='e4d.log',status='old',action='write',position='append')
         write(67,*) "SOLUTION CONVERGED"
         close(67)
      end if

      if(tag==6) then
         open(67,file='e4d.log',status='old',action='write',position='append')
         write(67,*) "__________ATTEMPTING FINAL SOLUTION ADJUSTMENT__________"
         close(67)
      end if

      if(tag==20) then
          open(67,file='e4d.log',status='old',action='write',position='append')
          write(67,*) "CANNOT FIND OUTPUT FILE: ",trim(outfile)
          write(67,*) "SKIPPING POTENTIAL OUTPUT"
          close(67)
       end if

       if(tag==21) then
          open(67,file='e4d.log',status='old',action='write',position='append')
          write(67,*) 
          write(67,*) " WRITING POTENTIAL FILES"
          close(67)
          write(*,*) " WRITING POTENTIAL FILES"
       end if


       if(tag==50) then
          open(67,file='e4d.log',status='old',action='write',position='append')
          write(67,*)
          write(67,*)
          write(67,*) "======================================================================="
          write(67,*) "++++++++++++++++++++ STARTING TIME-LAPSE_INVERSIONS +++++++++++++++++++"
          write(67,*) "======================================================================="
          close(67)
       end if

       if(tag==51) then
          open(67,file='e4d.log',status='old',action='write',position='append')
          write(67,*)
          write(67,*)
          write(67,*)
          write(67,*) "------------------------------------    SURVEY TIME: ",tlt(i_tl)
          write(67,*)
          write(67,*)
          write(67,*)
          close(67)
       end if

       if(tag==52) then
          open(67,file='e4d.log',status='old',action='write',position='append')   
          write(67,*) "_____________SOLUTION CONVERGED FOR TIME: ",tlt(i_tl)
          close(67)
       end if

       if(tag==53) then
          open(67,file='e4d.log',status='old',action='write',position='append')  
          write(67,*) "_____________ADJUSTING FINAL SOLUTION FOR TIME: ",tlt(i_tl)
          close(67)
       end if

       if(tag==54) then
          open(67,file='e4d.log',status='old',action='write',position='append')  
          write(67,*) "THE CHI-SQUARED VALUES ARE NOT ORDERED PROPERLY ... SKIPPING ADJUSTMENT"
          close(67)
       end if
       
       if(tag==55) then
         
         open(67,file='e4d.log',status='old',position='append')
         write(67,*)
         if(iter==0) then
            write(67,*) "********** DATA FIT AT STARTING MODEL *******************"
         else
            write(67,*) "********** FINAL SOLUTION ADJUSTMENT **************"
         end if
     
         write(67,"(4A15)") "Phi_dat","Phi_Mod","Phi_Mod/Beta","Phi_Tot"
         write(67,"(4F15.3)") phi_data,phi_model,phi_model/beta,phi_tot
         write(67,*) "Chi2 is currently",chi2,"Target value is ",norm_chi2    
         write(67,*) "Mean error is ",errm
         write(67,*) "RMS error is ",erms
         write(67,*) "Number of data culled is ",ncull
         if(cull_flag==1) then
            write(67,*) "PRIOR TO DATA CULLING"
            write(67,*) "Phi_dat is ",phi_data0
            write(67,*) "Chi2 is ", chi20
            write(67,*) "Mean error is ",errm0
            write(67,*) "RMS error is ",erms0
         end if
         
         write(67,*) "********************************************************"
         close(67)
      end if
       
      if(tag==56) then
         open(67,file='e4d.log',status='old',action='write',position='append')  
         write(67,*) "_____________ADJUSTING FINAL SOLUTION FOR TIME_________________"
         close(67)
      end if
      
      if(tag==57) then
         !get the meshfile prefix
         nchr=len_trim(mshfile)
         do i=1,nchr
            if(mshfile(i:i)=='.') then
               npre=i+1;
               exit
            end if
         end do

         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " The number of conductivity values is not equal to the number of mesh elements"
         write(67,*) " Number of elements in ",trim(mshfile(1:npre))//".ele is: ",nelem
         write(67,*) " Number of conductivity values in ",trim(sigfile)," is: ",nsig
         write(67,*) " Aborting ..."
         close(67)

         write(*,*)
         write(*,*) " The number of conductivity values is not equal to the number of mesh elements"
         write(*,*) " Number of elements in ",trim(mshfile(1:npre))//".ele is: ",nelem
         write(*,*) " Number of conductivity values in ",trim(sigfile)," is: ",nsig
         write(*,*) " Aborting ..."
      
         
      end if

      if(tag==58) then
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " The number of processors minus 1 must not exceed the number electrodes"
         write(67,*) " Number of electrodes specified in the survey file ",trim(efile)//" is: ",ne
         write(67,*) " Number of processors minus 1 is: ",n_rank - 1
         write(67,*) " To use more processors, add dummy electrodes to the survey file."
         write(67,*) " Aborting ..."
         close(67)
 
         write(*,*)
         write(*,*) " The number of processors minus 1 must not exceed the number electrodes"
         write(*,*) " Number of electrodes specified in the survey file ",trim(efile)//" is: ",ne
         write(*,*) " Number of processors minus 1 is: ",n_rank - 1
         write(*,*) " To use more processors, add dummy electrodes to the survey file."
         write(*,*) " Aborting ..."  
      end if

      if(tag==59) then
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " E4D requires at least 2 processors to run in modes greater than 1"
         write(67,"(A,I3.3)") "  Running in mode: ",mode
         write(67,"(A,I7.7)") "  Number of processors is: ",n_rank
         write(67,*) " Aborting ..."
         write(67,*)
         close(67)
 
         write(*,*)
         write(*,*) " E4D requires at least 2 processors to run in modes greater than 1"
         write(*,"(A,I3.3)") "  Running in mode: ",mode
         write(*,"(A,I7.7)") "  Number of processors is: ",n_rank
         write(*,*) " Aborting ..."
      end if

      if(tag==60) then
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " MESH SUMMARY"
         close(67)
      end if

      if(tag==61) then
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,"(A,I10.10)") "  Number of nodes:                  ",nnodes
         write(67,"(A,I10.10)") "  Number of metallic boundaries:    ",n_met
         
         close(67)
      end if

      if(tag==62) then
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,"(A,I10.10)") "  Number of elements:               ",nelem
         write(67,"(A,I10.10)") "  Number of zones:                  ",nrz        
         close(67)
      end if

      if(tag == 63) then
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,"(A,I10.10,A)") "  Element ",nrz, " has a zone label that is"
         write(67,*) " not greater than zero. All zone labels must"
         write(67,*) " be integers that are greater than zero or -999 for inactive elements."
         write(67,*) " Aborting ..."        
         close(67)
         write(*,*)
         write(*,"(A,I10.10,A)") "  Element ",nrz, " has a zone label that is"
         write(*,*) " not greater than zero. All zone labels must"
         write(*,*) " be integers that are greater than zero or -999 for inactive elements."
         write(*,*) " Aborting ..."   
         close(67)
      end if

      if(tag==64) then
         i1=0
         i2=0
         do i=1,n_rank-1
            if(eind(i,2)-eind(i,1)+1 > i1) i1 = eind(i,2)-eind(i,1)+1
            if(jind(i,2)-jind(i,1)+1 > i2) i2 = jind(i,2)-jind(i,1)+1
         end do
         open(67,file='e4d.log',status='old',action='write',position='append') 
         write(67,*) 
         write(67,*) " LOAD SUMMARY"
         write(67,"(A,I10.10)") "  Max # of pole soln's per core:    ",i1
         write(67,"(A,I10.10)") "  Max # of data rows per core:      ",i2
         i1=1e9
         i2=1e9
         do i=1,n_rank-1
            if(eind(i,2)-eind(i,1)+1 < i1) i1 = eind(i,2)-eind(i,1)+1
            if(jind(i,2)-jind(i,1)+1 < i2) i2 = jind(i,2)-jind(i,1)+1
         end do 
         write(67,"(A,I10.10)") "  Min # of pole soln's per core:    ",i1
         write(67,"(A,I10.10)") "  Min # of data rows per core:      ",i2
          
         close(67)
         
      end if

      if(tag==65) then
         write(*,*)
         write(*,*) " E4D: BUILDING FORWARD MATRIX MAPPING VECTORS "
      end if

      if(tag==66) then 
         write(*,*) " E4D: BUILDING FORWARD COUPLING MATRIX "
      end if

      if(tag==67) then
         write(*,*)             
         write(*,"(A,I3.3,A)") " ------------------------- E4D: ITERATION ",iter," --------------------------"
      end if

      if(tag==68) then
         write(*,*) " E4D: EXECUTING FORWARD RUN ..."
      end if

      if(tag==69) then
          open(67,file='e4d.log',status='old',action='write',position='append')
          write(67,*) 
          write(67,*) " WRITING SIMULATED SURVEY FILE: ",trim(sigfile),".srv"
          close(67)
          write(*,*) " WRITING SIMULATED SURVEY FILE: ",trim(sigfile),".srv"
       end if

       if(tag==70) then
          open(67,file='e4d.log',status='old',action='write',position='append')
          write(67,*) 
          write(67,*) " Forward run complete."
          close(67)
          write(*,*) " Forward run complete."
       end if

       if(tag==71) then
          write(*,*) " E4D: EXECUTING FORWARD IRUN ..."
       end if

       if(tag==72) then
          write(*,*) " E4D: BUILDING JACOBIAN MATRIX "
       end if

       if(tag==73) then
          write(*,"(A,I3.3,A)") " --------------------- E4D: END ITERATION ",iter," --------------------------"
       end if

      ! 100
      if(tag==100) then
         open(67,file='e4d.log',status='old',action='write',position='append')
         write(67,*) 
         write(67,*) "EXECUTING LINE SEARCH ON BETA"
         write(67,"(A12,A12,A12,A12,A12)") "BETA","PHI_DATA","PHI_MODEL","CHI_SQ","RMS_ERR"
         close(67)
      end if

      ! 101
      if(tag==101) then
         open(67,file='e4d.log',status='old',action='write',position='append')
         write(67,"(F12.3,F12.0,F12.0,F12.3,F12.3)") beta,phi_data,phi_model,chi2,erms
         close(67)
      end if
      
      ! 102
      if(tag==102) then
         open(67,file='e4d.log',status='old',action='write',position='append')
         write(67,*) "!!!!  UNABLE TO DECREASE OBJECTIVE FUNCTION   !!!!"
         close(67)
      end if
      
 
      

      
      
    end subroutine nreport
    !______________________________________________________________________________________

    !______________________________________________________________________________________
    subroutine treport(tag)
      implicit none
      integer :: tag
      if(mode==1) return;
      if(tag == -1) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         write(12,*)"____________________________________________________________________________"
         close(12)
      end if

      if(tag==0) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         write(12,*)"____________________________________________________________________________"
         write(12,*) "ITERATION: ",iter 
         close(12)
      end if

      if(tag==1) then
         open(12,file='run_time.txt',status='replace',action='write')
         write(12,*) "               F3D TIMING LOG    "
         write(12,*)
         if(etm>60.0) then
            write(12,*) " Setup time: ",etm/60," minutes"
         else
            write(12,*) " Setup time: ",etm," seconds"
         end if
         write(12,*) 
         close(12)
      end if

      if(tag==2) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(etm > 60.0) then
            write(12,*) 'Data assembly time         :',etm/60,' minutes'
         else
            write(12,*) 'Data assembly time         :',etm,' seconds'
         end if
         close(12)
      end if

      if(tag==3) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(jmin > 60.0) then   
            write(12,*) 'Minimum Jacobian build time: ',jmin/60,' minutes on slave :',sjmin
         else
            write(12,*) 'Minimum Jacobian build time: ',jmin,' seconds on slave :',sjmin
         end if

         if(jmax > 60.0) then
            write(12,*) 'Maximum Jacobian build time: ',jmax/60,' minutes on slave :',sjmax
         else
            write(12,*) 'Maximum Jacobian build time: ',jmax,' seconds on slave :',sjmax
         end if
         close(12)
      
      end if

      if(tag==4) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(etm > 60.0) then
            write(12,*) 'Parallel inversion time    :',etm/60,' minutes'
         else
            write(12,*) 'Parallel inversion time    :',etm,' seconds'
         end if
         close(12)
      end if

      if(tag==5) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         write(12,*) 'A build time at iteration       ',iter,' : ',etm1,' seconds'
         close(12)
      end if

      if(tag==6) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(rt_min > 60.0) then
            write(12,*) 'Minimum forward run time   : ',rt_min/60,' minutes on slave :',smin
         else
            write(12,*) 'Minimum forward run time   : ',rt_min,' seconds on slave :',smin
         end if

         if(rt_max > 60.0) then
            write(12,*) 'Maximum forward run time   : ',rt_max/60,' minutes on slave :',smax
         else
            write(12,*) 'Maximum forward run time   : ',rt_max,' seconds on slave :',smax
         end if
         close(12)
      end if
      
      if(tag==7) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(abmin > 60.0) then
            write(12,*) 'Minimum A-matrix build time: ',abmin/60,' minutes on slave :',sabmin
         else
            write(12,*) 'Minimum A-matrix build time: ',abmin,' seconds on slave :',sabmin
         end if
         
         if(abmax > 60.0) then
            write(12,*) 'Maximum A-matrix build time: ',abmax/60,' minutes on slave :',sabmax
         else
            write(12,*) 'Maximum A-matrix build time: ',abmax,' seconds on slave :',sabmax
         end if
         close(12)
      end if

      if(tag==8) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(kspmin > 60.0) then
            write(12,*) 'Minimum KSP build time     : ',kspmin/60,' minutes on slave :',skmin
         else
            write(12,*) 'Minimum KSP build time     : ',kspmin,' seconds on slave :',skmin
         end if
         
         if(kspmax > 60.0) then
            write(12,*) 'Maximum KSP build time     : ',kspmax/60,' minutes on slave :',skmax
         else
            write(12,*) 'Maximum KSP build time     : ',kspmax,' seconds on slave :',skmax
         end if
         close(12)
      end if

      
      if(tag==9) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         write(12,*) "EXECUTING LINE SEARCH FOR BETA"
         close(12)
      end if
      
      if(tag==10) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         write(12,*) "TESTING BETA AT: ",beta
         close(12)
      end if

      if(tag==11) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         write(12,*)
         write(12,*)"____________________________________________________________________________"
         if(etm > 60.0) then
            write(12,*) 'TOTAL RUN TIME             : ',etm/60,' minutes'
         else
            write(12,*) 'TOTAL RUN TIME             : ',etm,' seconds'
         end if
         close(12)
      end if


      if(tag==12) then
         open(12,file='run_time.txt',status='old',action='write',position='append')
         if(sbt >= 60.0) then
            write(12,*) ' Jaco. sequence build time  :',sbt/60,'minutes'
         else
            write(12,*) ' Jaco. sequence build time  :',sbt,'seconds'
         end if
         close(12)
      end if



    end subroutine treport
    
    
    subroutine testrun_report
		open(51,file='e4d.log',status='old',action='write',position='append')
		write(51,*) ""
		write(51,*) "******************* TEST RUN RESULTS *******************"
		write(51,*) ""
		if (cfg_flag) then
			write(51,*) "  The .cfg file entered in e4d.inp passed validity checks."
			write(51,*) ""
			write(51,*) "********************************************************"
			close(51)
			
			write(*,*) "******************* TEST RUN RESULTS *******************"
			write(*,*) ""
			write(*,*) "  The .cfg file entered in e4d.inp passed validity checks."
			write(*,*) ""
			write(*,*) "********************************************************"
		else
			write(51,*) "  The files entered in e4d.inp passed validity checks."
			write(51,*) ""
			write(51,*) "********************************************************"
			close(51)

			write(*,*) "******************* TEST RUN RESULTS *******************"
			write(*,*) ""
			write(*,*) "  The files entered in e4d.inp pass validity checks."
			write(*,*) ""
			write(*,*) "********************************************************"
			
		end if
		write(*,*) ""
		call crash_exit
		close(51)

    end subroutine testrun_report
    
    !______________________________________________________________________________________
end module report
