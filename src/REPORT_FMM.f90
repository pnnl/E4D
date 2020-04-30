module report_fmm

  use fmm_vars
  use vars
  use input_fmm
  
  contains

    !____________________________________________________________________________________
    subroutine init_log_fmm
      !This subroutine initializes fmm.log
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
      write(*,*)"************************ WELCOME TO FMM ************************"
      write(*,*)"Copyright © 2014, Battelle Memorial Institute"
      write(*,*) "All rights reserved."
      write(*,*)"Current date: ",trim(month)," ",trim(day),", ",trim(year)
      write(*,*)"Current time:  ",trim(hour),":",trim(minute),":",trim(second)
      write(*,"(A,I7.7,A)")" Running on ",n_rank," processing cores"
      !write(*,"(A,I7.7,A)")" (1 master process, and ",n_rank-1," slave processes)"
      write(*,*)"Please refer to fmm.log for further logging information ..."
      write(*,*)"****************************************************************"

      open(51,file='fmm.log',status='replace',action='write',IOSTAT=ios)
      if(ios.ne.0) then
         write(*,*)
         write(*,*) " THERE WAS A PROBLEM CREATING THE LOG FILE fmm.log"
         write(*,*) " Aborting ..."
         call crash_exit_fmm
         return
      end if
      write(51,*)"************************ WELCOME TO FMM ************************"
      write(51,*) "Copyright © 2014, Battelle Memorial Institute"
      write(51,*) "All rights reserved."
      write(51,*)"Current date: ",trim(month)," ",trim(day),", ",trim(year)
      write(51,*)"Current time:  ",trim(hour),":",trim(minute),":",trim(second)
      write(51,"(A,I7.7,A)")" Running on ",n_rank," processing cores"
      !write(51,"(A,I7.7,A)")" (1 master process, and ",n_rank-1," slave processes)"
      write(51,*)"****************************************************************"

      close(51)

      

    end subroutine init_log_fmm
    !____________________________________________________________________________________

    !____________________________________________________________________________________
    subroutine nreport_fmm(tag)
      implicit none
      integer :: tag,i,npre,nchr,i1,i2,i3
      character*20 :: fstring

      if(mode==1) return;
      if(tag==1) then
         
         open(67,file='fmm.log',status='old',position='append')
         write(67,*)
         
         write(67,*)"***********************************************************************"
         close(67)
      end if

      ! 2
      if(tag==2) then   
         write(*,*) " FMM: EXECUTING FORWARD RUN"
      end if


      ! 3
      if(tag==3) then
      end if

      ! 4
      if(tag==4) then
         open(67,file='fmm.log',status='old',action='write',position='append')
         write(67,*) "SOLUTION CONVERGED"
         close(67)
      end if

      if(tag==6) then
         open(67,file='fmm.log',status='old',action='write',position='append')
         write(67,*) "__________ATTEMPTING FINAL SOLUTION ADJUSTMENT__________"
         close(67)
      end if

      if(tag==20) then
          open(67,file='fmm.log',status='old',action='write',position='append')
          write(67,*) "CANNOT FIND OUTPUT FILE: ",trim(outfile)
          write(67,*) "SKIPPING POTENTIAL OUTPUT"
          close(67)
       end if

       if(tag==21) then
          open(67,file='fmm.log',status='old',action='write',position='append')
          write(67,*) 
          write(67,*) " WRITING TRAVEL TIME FILES"
          close(67)
          write(*,*) " WRITING TRAVEL TIME FILES"
       end if



       if(tag==51) then
          open(67,file='fmm.log',status='old',action='write',position='append')
          write(67,*)
          write(67,*)
          write(67,*)
          write(67,*) "_____________SURVEY TIME: ",tlt_fmm(i_tl_fmm)
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

         open(67,file='fmm.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " The number of slowness values is not equal to the number of mesh elements"
         write(67,*) " Number of elements in ",trim(mshfile(1:npre))//".ele is: ",nelem
         write(67,*) " Number of slowness values in ",trim(spdfile)," is: ",nspd
         write(67,*) " Aborting ..."
         close(67)

         write(*,*)
         write(*,*) " The number of slowness values is not equal to the number of mesh elements"
         write(*,*) " Number of elements in ",trim(mshfile(1:npre))//".ele is: ",nelem
         write(*,*) " Number of slowness values in ",trim(spdfile)," is: ",nspd
         write(*,*) " Aborting ..."
      
         
      end if

      if(tag==59) then
         open(67,file='fmm.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " FMM requires at least 2 processors to run in modes greater than 1"
         write(67,"(A,I3.3)") "  Running in mode: ",mode
         write(67,"(A,I7.7)") "  Number of processors is: ",n_rank
         write(67,*) " Aborting ..."
         write(67,*)
         close(67)
 
         write(*,*)
         write(*,*) " FMM requires at least 2 processors to run in modes greater than 1"
         write(*,"(A,I3.3)") "  Running in mode: ",mode
         write(*,"(A,I7.7)") "  Number of processors is: ",n_rank
         write(*,*) " Aborting ..."
      end if

      if(tag==60) then
         open(67,file='fmm.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,*) " MESH SUMMARY"
         close(67)
      end if

      if(tag==61) then
         open(67,file='fmm.log',status='old',action='write',position='append') 
         write(67,"(A,I10.10)") "  Number of nodes:                  ",nnodes
         close(67)
      end if

      if(tag==62) then
         open(67,file='fmm.log',status='old',action='write',position='append') 
         write(67,"(A,I10.10)") "  Number of elements:               ",nelem
         write(67,"(A,I10.10)") "  Number of zones:                  ",nrz        
         close(67)
      end if

      if(tag == 63) then
         open(67,file='fmm.log',status='old',action='write',position='append') 
         write(67,*)
         write(67,"(A,I10.10,A)") "  Element ",nrz, " has a zone label that is"
         write(67,*) " not greater than zero. All zone labels must"
         write(67,*) " be integers that are greater than zero."
         write(67,*) " Aborting ..."        
         close(67)
         write(*,*)
         write(*,"(A,I10.10,A)") "  Element ",nrz, " has a zone label that is"
         write(*,*) " not greater than zero. All zone labels must"
         write(*,*) " be integers that are greater than zero."
         write(*,*) " Aborting ..."   
         close(67)
      end if


      if(tag==68) then
         write(*,*) " FMM: EXECUTING FORWARD RUN ..."
      end if

      if(tag==69) then
          open(67,file='fmm.log',status='old',action='write',position='append')
          write(67,*) 
          write(67,*) " WRITING SIMULATED SURVEY FILE: ",trim(spdfile),".srv"
          close(67)
          write(*,*) " WRITING SIMULATED SURVEY FILE: ",trim(spdfile),".srv"
       end if

       if(tag==70) then
          open(67,file='fmm.log',status='old',action='write',position='append')
          write(67,*) 
          write(67,*) " Forward run complete."
          close(67)
          write(*,*) " Forward run complete."
       end if

       if(tag==71) then
          write(*,*) " EXECUTING FORWARD IRUN ..."
       end if

       if(tag==72) then
          write(*,*) " BUILDING JACOBIAN MATRIX "
       end if

       if(tag==73) then
          write(*,"(A,I3.3,A)") "---------------------- END ITERATION ",iter," ----------------------"
       end if

       if(tag==74) then
          open(67,file='fmm.log',status='old',action='write',position='append')
          write(67,*)
          write(67,*)"-----------------------------------------------------------------"
          write(67,*) "Computing the FMM Jacobian Matrix "
          close(67)
       end if

      ! 100
      
    end subroutine nreport_fmm
    !______________________________________________________________________________________

end module report_fmm
