module input_fmm

  use fmm_vars
  use vars
  use input
  implicit none

  !Variables needed only by the master process are declared here
  !integer :: ios                                           !!io status
  !integer :: mnchar                                        !!used to record number of character in a string
  !real :: xorig,yorig,zorig                                !!x,y,z mesh translation values
  !integer :: i_tl_fmm                                       !!current time lapse file index
  !real, dimension(:) , allocatable :: tlt_fmm               !!time lapse data time markers
  integer, dimension(:,:), allocatable :: nbrs             

contains

  !____________________________________________________________________________________________
  subroutine read_input_fmm
    implicit none
    
    character*40 :: smode
    integer :: nchar,junk,i,check,j
    logical :: exst
    logical :: wdwarn = .false.


    call check_inp_fmm(0,junk)
    open(10,file='fmm.inp',status='old',action='read') 
    read(10,*,IOSTAT=ios) smode; call check_inp_fmm(101,junk)

    ! check for alpha characters and make them all lower case
    smode=lcase(smode)
    select case(trim(adjustl(smode)))
           case('fmm1')
              mode_fmm = 1; call check_inp_fmm(1,junk)
           case('fmm2')
              mode_fmm = 2; call check_inp_fmm(1,junk)
           case('fmm3')
              mode_fmm = 3; call check_inp_fmm(1,junk)
           case('fmm4')
              mode_fmm = 4; call check_inp_fmm(1,junk)
           case default
              read(smode,*,IOSTAT=ios) mode_fmm; call check_inp_fmm(1,junk)
    end select
    
    !read mesh file
    read(10,*,IOSTAT=ios) mshfile;         call check_inp_fmm(2,junk)
       
    !if mode is > 1 then read the zones to be used in the simulation
    if(mode_fmm > 1) then
       backspace(10)
       read(10,*,IOSTAT=ios) smode, nzf;  call check_inp_fmm(54,junk)
       if(ios .ne. 0) then
          nzf = 0
          backspace(10)
       else
          backspace(10)
          allocate(zsims(nzf))
          read(10,*,IOSTAT=ios) smode,nzf,zsims(1:nzf); call check_inp_fmm(55,junk)
          if(ios.ne.0) then
             nzf = 0
             backspace(10)
             deallocate(zsims)
          end if
       end if

       ! read survey file       
       read(10,*,IOSTAT=ios) tfile;           call check_inp_fmm(3,junk)       
       ! read velocity file       
       read(10,*,IOSTAT=ios) spdfile;         call check_inp_fmm(4,junk)       
       read(10,*,IOSTAT=ios) outfile_fmm;     call check_inp_fmm(5,junk)              
    end if    
    
    if(mode_fmm == 3 .or. mode_fmm == 4) then
       read(10,*,IOSTAT=ios) invfile;         call check_inp_fmm(6,junk)
       read(10,*,IOSTAT=ios) refmod_file;     call check_inp_fmm(7,junk)
    end if

    ! Time-lapse
    if(mode_fmm == 4) tl_ly = .true.
        
    !!read the time lapse data list file
    if (tl_ly) then
       mode_fmm = 3
       check = 0

       read(10,*,IOSTAT=ios) tl_file,check;   call check_inp_fmm(8,check)
       if(check == 2) r_last = .true.
      
       if(rt_flag) then
          ntl=1 
          allocate(tl_dfils(ntl),tlt(ntl))   
       else
          open(21,file=trim(tl_file),status='old',action='read',IOSTAT=ios)  
          read(21,*,IOSTAT=ios) ntl; call check_inp_fmm(9,junk)
          allocate(tl_dfils(ntl),tlt(ntl))
          
          do i=1,ntl
             read(21,*,IOSTAT=ios) tl_dfils(i), tlt(i); call check_inp_fmm(10,i)
             inquire(file=trim(tl_dfils(i)),exist=exst)
             if(.not. exst) call check_inp_fmm(11,i)
          end do
          close(21)
       end if
    endif    
    
    close(10)

    !!Determine if the mesh file is a .cfg file or if meshfiles are provided
    mnchar = 0
    do i=1,40
       if(mshfile(i:i) == '.') then
          mnchar = i
          exit
       end if
    end do
    if(mnchar == 0) call check_inp_fmm(21,0)

    !!Check for compatibility between mshfile and mode
    !!Check for config file
    if (mode_fmm == 1 ) then
       if(mshfile(mnchar+1:mnchar+3) == "cfg") then
          inquire(file=trim(mshfile),exist=exst)
          if(.not.exst) call check_inp_fmm(21,1)
       else if(mshfile(mnchar+1:mnchar+3).ne."cfg") then
          call check_inp_fmm(21,0)
       end if

       ! Both E4D and FMM can't build mesh file together so
       if (simulate_e4d) call check_inp_fmm(58,junk)
       
       ! return for mesh generation mode       
       return       
    end if     
    
    !!check mesh files
    if (mode_fmm > 1)then 
       call check_inp_fmm(121,mnchar)
       call check_inp_fmm(122,mnchar)
    end if
    
    !!Allocate/read the source positions and survey configuration
    if(mode_fmm>1) then
       !if we're inverting then check to see if we're doing fresnel volume inversion
       !or not
       !call set_fresnel
       call read_survey_fmm
       call translate_source
    end if

    !!read the velocity file
    call read_velocity

    ! examine min and max Fresnel volume
    !if (fresnel) call check_fresnel

  end subroutine read_input_fmm
  !___________________________________________________________________________________


  !____________________________________________________________________________________________
  subroutine check_inp_fmm(spot,indx)
    implicit none
    integer :: spot,indx
    logical :: exst, mcheck

    select case(spot)

    case(0)
       inquire(file='fmm.inp',exist=exst)
       if(.not. exst) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) "FMM: Cannot find the primary input file fmm.inp: aborting"
          write(*, *) "FMM: Cannot find the primary input file fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       end if

    case(101)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) "FMM: There was a problem reading the mode in fmm.inp: aborting"
          write(*, *) "FMM: There was a problem reading the mode in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       end if

    case(1)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) " FMM: There was a problem reading the mode in fmm.inp."
          write(51,*) " FMM: Aborting ..."
          write(*, *) " FMM: There was a problem reading the mode in fmm.inp."
          write(*, *) " FMM: Aborting ..."
          close(51)
          call crash_exit_fmm
       end if
       mcheck = .false.
       select case(mode_fmm)
       case(0) 
       case(1)
          mcheck = .true.
       case(2) 
          mcheck = .true.
       case(3)
          mcheck = .true.
       case(4)
          mcheck = .true.
       end select

       if(.not.mcheck) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,"(A33,I3,A19)") "The mode selected in fmm.inp is ",mode_fmm,", which is invalid."
          write(51,*)"Valid run modes include: < see below > "
          write(51,*)  "---------------------------------------"
          write(51,*)  " -FUNCTION-                      -MODE-"
          write(51,*)  "---------------------------------------"
          write(51,*)  " FMM Mesh Generation          1 or FMM1"
          write(51,*)  " FMM Forward                  2 or FMM2"
          write(51,*)  " FMM Inversion                3 or FMM3"
          write(51,*)  " FMM Time-lapse Inversion     4 or FMM4"
          write(51,*)  "---------------------------------------"
          write(51,*)  " Aborting ..."
          write(*,"(A33,I3,A19)") "The mode selected in fmm.inp is ",mode_fmm,", which is invalid."
          write(*, *)"Valid run modes include: < see below > "
          write(*, *)  "---------------------------------------"
          write(*, *)  " -FUNCTION-                      -MODE-"
          write(*, *)  "---------------------------------------"
          write(*, *)  " FMM Mesh Generation          1 or FMM1"
          write(*, *)  " FMM Forward                  2 or FMM2"
          write(*, *)  " FMM Inversion                3 or FMM3"
          write(*, *)  " FMM Time-lapse Inversion     4 or FMM4"
          write(*, *)  "---------------------------------------"
          write(*, *)  " Aborting ..."
          close(51)
          call crash_exit_fmm
         
       else
          
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) 
             select case(mode_fmm)
             case(0) 
             case(1) 
             case(2) 
                write(51,*) "***************** RUNNING IN FMM FORWARD MODE *******************"
                write(51,"(A,I3.3)") "  Mode:                             ",mode_fmm
             case(3)
             	write(51,*) "**************** RUNNING IN FMM INVERSION MODE ******************"
                write(51,"(A,I3.3)") "  Mode:                             ",mode_fmm
             case(4)                
             	write(51,*) "********** RUNNING IN FMM Time-lapse INVERSION MODE *************"
             	write(51,"(A,I3.3)") "  Mode:                             ",mode_fmm
             end select
          close(51)
       end if

    case(2)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then

             open(51,file='fmm.log',status='old',action='write',position='append')
             write(51,*) "  FMM: There was a problem reading the mesh file name in fmm.inp"
             write(51,*) "  FMM: Aborting ..."
             close(51)
             write(*, *) "  FMM: There was a problem reading the mesh file name in fmm.inp"
             write(*, *) "  FMM: Aborting ..."
          call crash_exit_fmm

       else
             write(51,*) " Mesh file:                        ",trim(mshfile)
       end if
       close(51)

    case(3)

       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the survey file name in fmm.inp: aborting"
          write(*, *) "FMM: There was a problem reading the survey file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Survey configuration file:        ",trim(tfile)
       end if
       close(51)

    case(4)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the velocity file name in fmm.inp: aborting"
          write(*, *) "FMM: There was a problem reading the velocity file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Velocity file:                       ",trim(spdfile)
       end if
       close(51)

    case(5)
      
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the output options file name in fmm.inp: aborting"
          write(*, *) "FMM: There was a problem reading the output options file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Output options file:              ",trim(outfile_fmm)
       end if
       close(51)

    case(6)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the inverse options file name in fmm.inp: aborting"
          write(*, *) "FMM: There was a problem reading the inverse options file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Inverse options file:             ",trim(invfile)
       endif
       close(51)

    case(7)
       open(51,file='fmm.log',status='old',action='write',position='append')      
       if (ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the reference model file name in fmm.inp: aborting"
          close(51)
          write(*, *) "FMM: There was a problem reading the reference model file name in fmm.inp: aborting"
          call crash_exit_fmm          
       else
          write(51,*) " Reference model file:             ",trim(refmod_file)
       end if        
       close(51) 

    case(8)
       open(51,file='fmm.log',status='old',action='write',position='append')       
       if (ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the time-lapse survey file name "
          write(51,*) "FMM: and/or reference model update option in fmm.inp: Aborting ..."
          close(51)
          write(*, *) "FMM: There was a problem reading the time-lapse survey file name "
          write(*, *) "FMM: and/or reference model update option in fmm.inp: Aborting ..."
          call crash_exit_fmm          
       else
          write(51,*) " Time-lapse survey file:           ",trim(tl_file)
       end if

       if (indx==2) then
          write(51,*) " Reference model update option:    previous solution"
       else
          write(51,*) " Reference model update option:    baseline solution"
       end if

       inquire(file=trim(tl_file),exist=exst)       
       if (.not.exst .and. .not. rt_flag) then
          write(*, *) 
          write(51,*) "FMM: Cannot find the time lapse survey list file: ",trim(tl_file)," ... aborting"
          write(*, *) "FMM: Cannot find the time lapse survey list file: ",trim(tl_file)," ... aborting"
          close(51)
          call crash_exit_fmm
       end if
       
       close(51) 

    case(9)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if (ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading the number of survey files in ",trim(tl_file)," :aborting"
          write(*, *) "FMM: There was a problem reading the number of survey files in ",trim(tl_file)," :aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*)
          write(51,*) " Time-lapse survey list file:      ",trim(tl_file)
          write(51,*) " Num. of time-lapse survey files:  ",ntl            
       end if       
       close(51)
              
    case(10)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if (ios .ne. 0) then
          write(51,*) "FMM: There was a problem reading time-lapse file at time ",indx
          write(51,*) "FMM: in the time lapse survey file: ",trim(tl_file)
          write(51,*) "FMM: Aborting ... "
          write(*, *) "FMM: There was a problem reading time-lapse file at time ",indx
          write(*, *) "FMM: in the time lapse survey file: ",trim(tl_file)
          write(*, *) "FMM: Aborting ... "
          close(51)
          call crash_exit_fmm
       else         
          if (indx==1) then
             write(51,*) " Index     Time-Lapse Survey File       Time Stamp"            
          end if          
          write(51,"(I6,A28,g17.5)"),indx,trim(tl_dfils(indx)),tlt(indx)          
       end if       
      
       close(51)       
       
    case(11)
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*) "FMM: There was a problem finding the time-lapse file: ",trim(tl_dfils(indx))
       write(51,*) "FMM: aborting ..."
       write(*, *) "FMM: There was a problem finding the time-lapse file: ",trim(tl_dfils(indx))
       write(*, *) "FMM: aborting ..."
       close(51)
       call crash_exit_fmm
              
    case(12)

       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) 
          write(51,*) " FMM: There was a problem reading the number of sources in the survey file ",trim(tfile)
          write(51,*) " FMM: Aborting ..."
          write(*,*) 
          write(*,*) " FMM: There was a problem reading the number of sources in the survey file ",trim(tfile)
          write(*,*) " FMM: Aborting ..."
          
          close(51)
          call crash_exit_fmm
       else
          write(51,"(A,I7.7)") "  Number of sources:                ",ns
       end if
       close(51)

    case(13)      
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,"(A,I7.7)") "  There was a problem reading source number ",indx
          write(51,*) " FMM: in the survey file file: ",trim(tfile)
          write(51,*) " FMM: Aborting ..."
          if(fresnel) then
             write(51,*) "Running in fresnel mode ... be sure to include the positive"
             write(51,*) "frequency value in the last column"
          end if
          close(51)
          write(*,*)
          write(*,"(A,I7.7)") " There was a problem reading source number ",indx
          write(*,*) "FMM: in the survey file file: ",trim(tfile)
          write(*,*) "FMM: Aborting ..."
          if(fresnel) then
             write(*,*) "Running in fresnel mode ... be sure to include the positive"
             write(*,*) "frequency value in the last column"
          end if
          call crash_exit_fmm
       end if

       if (fresnel) then
          if ((frq(indx).lt.0.01) .or. (frq(indx).gt.100.0)) then
             open(51,file='fmm.log',status='old',action='write',position='append')
             write(51,*)
             write(51,*) "Running in fresnel mode. The dominant frequency should be within 0.01-100."
             write(51,*) "The frequency for source number",indx," in survey file: ",trim(tfile)," is: ",frq(indx)
             write(51,*) "Aborting ..."
             close(51)
             write(*, *) "Running in fresnel mode. The dominant frequency should be within 0.01-100."
             write(*, *) "The frequency for source number",indx,"in survey file: ",trim(tfile)," is:",frq(indx)
             write(*,*) "Aborting ..."
             call crash_exit_fmm
          endif
       endif
       

    case(14)
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) " FMM: The source index specified for source: ",indx
       write(51,*) " FMM: is greater than the total number of sources."
       write(51,*) " FMM: Aborting ..."
       write(*, *) 
       write(*, *) " FMM: The source index specified for source: ",indx
       write(*, *) " FMM: is greater than the total number of sources."
       write(*, *) " FMM: Aborting ..."
     
       close(51)
       
       call crash_exit_fmm

    case(112)

       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) 
          write(51,*) " FMM: There was a problem reading the number of receivers in the survey file ",trim(tfile)
          write(51,*) " FMM: Aborting ..."
          write(*,*) 
          write(*,*) " FMM: There was a problem reading the number of receivers in the survey file ",trim(tfile)
          write(*,*) " FMM: Aborting ..."
          
          close(51)
          call crash_exit_fmm
       else
          write(51,"(A,I7.7)") "  Number of receivers:             ",nrc
       end if
       close(51)

    case(113)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,"(A,I7.7)") "  There was a problem reading receiver number ",indx
          write(51,*) " in the survey file file: ",trim(tfile)
          write(51,*) " Aborting ..."
          close(51)
          write(*,*)
          write(*,"(A,I7.7)") " There was a problem reading receiver number ",indx
          write(*,*) "in the survey file file: ",trim(tfile)
          write(*,*) "Aborting ..."
          call crash_exit_fmm
       end if

    case(114)
       open(51,file='fmm.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) " FMM: The receiver index specified for receiver: ",indx
       write(51,*) " FMM: is greater than the total number of receivers."
       write(51,*) " FMM: Aborting ..."
       write(*, *) 
       write(*, *) " FMM: The receiver index specified for receiver: ",indx
       write(*, *) " FMM: is greater than the total number of receivers."
       write(*, *) " FMM: Aborting ..."
     
       close(51)
       
       call crash_exit_fmm
    case(15)
       inquire(file=mshfile(1:mnchar)//'trn',exist=exst)
       if(.not. exst) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " FMM: Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn' 
          write(51,*) " FMM: Aborting..."
          write(*, *)
          write(*, *) " FMM: Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn' 
          write(*, *) " FMM: Aborting..."
          close(51)
          call crash_exit_fmm
       end if

    case(16)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " FMM: There was a problem reading the mesh "
          write(51,*) " FMM: translation numbers in: ",trim(mshfile(1:mnchar))//'trn' 
          write(51,*) " FMM: Aborting ... "
          close(51)
          write(*, *)
          write(*, *) " FMM: There was a problem reading the mesh "
          write(*, *) " FMM: translation numbers in: ",trim(mshfile(1:mnchar))//'trn' 
          write(*, *) " FMM: Aborting ... "
          close(51)
          call crash_exit_fmm
       end if

    case(17)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " FMM: There was a problem reading the number of measurements" 
          write(51,*) " FMM: in the survey file: ",trim(tfile)
          close(51)
          write(*,*)
          write(*,*) " FMM: There was a problem reading the number of measurements" 
          write(*,*) " FMM: in the survey file: ",trim(tfile)
          close(51)
          call crash_exit_fmm
       elseif(nm_fmm .le. 0) then
          write(51,*) " FMM: The number of measurements is not positive: ",indx
          write(51,*) " FMM: Aborting ..."
          write(*,*)
          write(*,*) " FMM: The number of measurements is not positive: ",indx
          write(*,*) " FMM: Aborting ..."
          close(51)
          call crash_exit_fmm
       else
           write(51,"(A,I7.7)") "  Number of measurements:           ",nm_fmm
       end if
       close(51)
       
    case(18)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
             write(51,*)
             write(51,"(A,I8.8)")" There was a problem reading measurement number ",indx
             write(51,*) " in the survey file: ",trim(tfile)
             write(51,*) " aborting ..."
             write(*,*)
             write(*,"(A,I8.8)") "  There was a problem reading measurement number ",indx
             write(*,*) " in the survey file: ",trim(tfile)
             write(*,*) " aborting ..."
          close(51)
          call crash_exit_fmm
       end if

   case(19)
      open(51,file='fmm.log',status='old',action='write',position='append')
      write(51,*)
      write(51,"(A,I8.8)") "  !!! WARNING: MEASUREMENT ",indx
      write(51,*) " !!! AND POSSIBLY OTHERS SPECIFY A NEGATIVE OR ZERO STANDARD DEVIATION"
      write(51,*) " !!! IN THE SURVEY FILE ",trim(tfile)
      write(51,*) " !!! SETTING TO LARGE STANDARD DEVIATION"
      close(51)

      write(*,*)
      write(*,"(A,I8.8)") " !!! WARNING: MEASUREMENT ",indx
      write(*,*) "!!! AND POSSIBLY OTHERS SPECIFY A NEGATIVE OR ZERO STANDARD DEVIATION"
      write(*,*) "!!! IN THE SURVEY FILE ",trim(tfile)
      write(*,*) "!!! SETTING TO LARGE STANDARD DEVIATION" 
      close(51)
      return

   case(20)
      open(51,file='fmm.log',status='old',action='write',position='append')
      write(51,*)
      write(51,"(A,I8.8,A)") "  Measurement ",indx," uses sources/receivers that are out of range." 
      write(51,"(A,I8.8,A)") "  There are ",ns," sources."
      write(51,"(A,I8.8,A,4I10.8)") "  A B for measurement ",indx," reads ",s_conf_fmm(indx,:)
      write(51,*) " Aborting ..."
      write(51,*)
      close(51)

      write(*,*)
      write(*,"(A,I8.8,A)") "  Measurement ",indx," uses source that are out of range." 
      write(*,"(A,I8.8,A)") "  There are ",ns," sources."
      write(*,"(A,I8.8,A,4I10.8)") "  A B for measurement ",indx," reads ",s_conf_fmm(indx,:)
      write(*,*) " Aborting ..."
      write(*,*)
      call crash_exit_fmm
      
 
   case(21)
      open(51,file='fmm.log',status='old',action='write',position='append')
      if(indx==0) then
         write(51,*)
         write(51,*) " FMM: In mode = 1 or FMM1 you must provide a mesh configuration (*.cfg) file "
         write(51,*) " FMM: In all other modes you must provide a node or element file name"
         write(51,*) " FMM: with .*.node or .*.ele file name extension where * is a single letter. "
         write(51,*) " FMM: You provided: ",trim(mshfile)
         write(51,*) " FMM: Aborting ..."
         write(*, *)
         write(*, *) " FMM: In mode = 1 or FMM1 you must provide a mesh configuration (*.cfg) file "
         write(*, *) " FMM: In all other modes you must provide a node or element file name"
         write(*, *) " FMM: with .*.node or .*.ele file name extension where * is a single letter. "
         write(*, *) " FMM: You provided: ",trim(mshfile)
         write(*, *) " FMM: Aborting ..."
      else if(indx==1) then
         write(51,*) 
         write(51,*) " FMM: Cannot find the mesh configuration file: ",trim(mshfile)
         write(51,*) " FMM: Aborting ..."
         write(*, *) 
         write(*, *) " FMM: Cannot find the mesh configuration file: ",trim(mshfile)
         write(*, *) " FMM: Aborting ..."
      end if
      close(51)
      call crash_exit_fmm
      
   case(22)
      if(ios .ne. 0) then
         open(51,file='fmm.log',status='old',action='write',position='append')
         write(51,*) "FMM: There was a problem reading the number of velocities in"
         write(51,*) "FMM: in the velocity file: ",trim(spdfile)
         write(51,*) "FMM: aborting."
         close(51)
         call crash_exit_fmm
      end if


   case(23)
      if(indx == 0) then
         open(51,file='fmm.log',status='old',action='write',position='append')
         if(ios == 0) then
            write(51,*)
            write(51,*) " VELOCITY FILE SUMMARY "
            write(51,"(A,I10.10)") "  Number of velocity values:    ",nspd
            close(51)
         else
            write(51,*) 
            write(51,*) " FMM: There was a problem reading the number of "
            write(51,*) " FMM: velocity values in ",trim(spdfile)
            write(51,*) " FMM: Aborting ..."
            close(51)
            write(*,*) 
            write(*,*) " FMM: There was a problem reading the number of "
            write(*,*) " FMM: velocity values in ",trim(spdfile)
            write(*,*) " FMM: Aborting ..."
            close(51)
            call crash_exit_fmm
         end if
      end if
      if(ios .ne. 0) then
         open(51,file='fmm.log',status='old',action='write',position='append')
         write(51,*)
         write(51,"(A,I10.10)") "  There was a problem reading velocity number ",indx
         write(51,*) " in the velocity file: ",trim(spdfile),"."
         write(51,*) " Aborting ..."
         close(51)
         write(*,*)
         write(*,"(A,I10.10)") "  There was a problem reading source number ",indx
         write(*,*) " in the source file: ",trim(spdfile),"."
         write(*,*) " Aborting ..."
         close(51)
         call crash_exit_fmm
      end if

   case(24)
      open(51,file='fmm.log',status='old',action='write',position='append')
      if(indx == 0) then
         write(51,*)
         write(51,*) " FMM: Can't find the survey file: ",trim(tfile)
         write(51,*) " FMM: aborting."
         write(*,*)
         write(*,*) " FMM: Can't find the survey file: ",trim(tfile)
         write(*,*) " FMM: aborting."
         close(51)
         call crash_exit_fmm
      else         
         write(51,*) 
         write(51,*) " SURVEY FILE SUMMARY"
      end if
      close(51)
      
   case(25)     
      open(51,file='fmm.log',status='old',action='write',position='append')      
      write(51,*) 
      write(51,*) " FMM: Can't find the velocity file: ",trim(spdfile)
      write(51,*) " FMM: Aborting ..."
      close(51)
      write(*,*) 
      write(*,*) " FMM: Can't find the velocity file: ",trim(spdfile)
      write(*,*) " FMM: Aborting ..."
      close(51)
      call crash_exit_fmm

   case(52)      
      if(ios.ne.0 .or. indx .ne. 0 .or. indx .ne. 1) then
         open(51,file='fmm.log',status='old',action='write',position='append')
         write(51,*) "FMM: There was a problem reading the first line of the "
         write(51,*) "FMM: fmm survey file ",trim(tfile)
         write(51,*) "FMM: In fmm mode, the first line of the survey file"
         write(51,*) "FMM: should contain two integers: the number of source positions "
         write(51,*) "FMM: and the fresnel volume flag"
         write(51,*) "FMM: 0 for ray-based and 1 for fresnel volume"
         write(51,*) "FMM: Using ray-based by default "
         write(51,*) "FMM: Aborting ..."
         close(51)
         write(*,*)
         write(*,*) "FMM: There was a problem reading the first line of the "
         write(*,*) "FMM: fmm survey file ",trim(tfile)
         write(*,*) "FMM: In fmm mode, the first line of the survey file"
         write(*,*) "FMM: should contain two integers: the number of source positions "
         write(*,*) "FMM: and the fresnel volume flag"
         write(*,*) "FMM: 0 for ray-based and 1 for fresnel volume"
         write(*,*) "FMM: Aborting ..."
         call crash_exit_fmm
      end if
                
   case(53)      
      open(51,file='fmm.log',status='old',action='write',position='append')      
      write(51,*) "FMM: The frequency for source number: ",indx      
      write(51,*) "FMM: Is less than or equal to zero. Specified frequencies"      
      write(51,*) "FMM: must be positive."      
      write(51,*) "FMM: Aborting..."      
      close(51)
      write(*,*) "FMM: The frequency for source number: ",indx
      write(*,*) "FMM: Is less than or equal to zero. Specified frequencies"
      write(*,*) "FMM: must be positive."
      write(*,*) "FMM: Aborting..."
      call crash_exit_fmm

   case(54)      
      if(ios .ne. 0) then         
         open(51,file='fmm.log',status='old',action='write',position='append')         
         write(51,*)         
         write(51,*) '------------------------------- WARNING -------------------------------'         
         write(51,*) "There was a problem reading the number of zones to include"         
         write(51,*) "in the forward travel time simulation after the mesh file name."         
         write(51,*) "Using all zones."         
         write(51,*) '-----------------------------------------------------------------------'         
         write(51,*)         
         write(*, *)         
         write(*, *) '------------------------------- WARNING -------------------------------'        
         write(*, *) "There was a problem reading the number of zones to include"         
         write(*, *) "in the forward travel time simulation after the mesh file name."         
         write(*, *) "Using all zones."         
         write(*, *) '-----------------------------------------------------------------------'         
         write(*, *)         
      end if      

   case(55)      
      if(ios .ne. 0) then
         open(51,file='fmm.log',status='old',action='write',position='append')
         write(51,*)
         write(51,*) '------------------------------- WARNING -------------------------------'  
         write(51,*) "There was a problem reading which zones to include"
         write(51,*) "in the forward travel time simulation."
         write(51,*) "Using all zones" 
         write(51,*) '-----------------------------------------------------------------------'
         write(51,*)
         write(*, *)
         write(*, *) '------------------------------- WARNING -------------------------------'              
         write(*, *) "There was a problem reading which zones to include"
         write(*, *) "in the forward travel time simulation."
         write(*, *) "Using all zones"
         write(*, *) '-----------------------------------------------------------------------'
         write(*, *)
      end if      

   case(56)
      open(51,file='fmm.log',status='old',action='write',position='append')
      write(51,*) '-------------------------------- ERROR --------------------------------'
      write(51,*) "There was a problem reading the source position and/or the dominant"
      write(51,'(2(A,I8.5))') " frequency for source number: ",indx-1,"   or ",indx
      write(51,*) "Please make sure the x-, y-, and z-positions and frequency for "
      write(51,'(2(A,I8.5))') " for source number: ",indx-1,"   or ", indx
      write(51,*) "are specified correctly in fmm survey file: ",trim(tfile)
      write(51,*) "Aborting..."
      write(51,*) '-----------------------------------------------------------------------'
      close(51)
      write(*, *) '-------------------------------- ERROR --------------------------------'
      write(*, *) "There was a problem reading the source position and/or the dominant"
      write(*, '(2(A,I8.5))') " frequency for source number: ",indx-1,"   or ",indx
      write(*, *) "Please make sure the x-, y-, and z-positions and frequency for "
      write(*, '(2(A,I8.5))') " for source number: ",indx-1,"   or ", indx
      write(*, *) "are specified correctly in fmm survey file: ",trim(tfile)
      write(*, *) "Aborting..."      
      write(*, *) '-----------------------------------------------------------------------'      
      call crash_exit_fmm      
             
   case(57)      
      open(51,file='fmm.log',status='old',action='write',position='append')      
      write(51,*) '-------------------------------- ERROR --------------------------------'      
      write(51,*) "There was a problem reading the receiver position and/or data"      
      write(51,'(2(A,I8.5))') " for receiver number: ",indx-1,"   or ",indx      
      write(51,*) "Please make sure the x-, y-, and z-positions and/or data for "      
      write(51,'(2(A,I8.5))') " for receiver number: ",indx-1,"   or ", indx      
      write(51,*) "are specified correctly in fmm survey file: ",trim(tfile)      
      if (indx.eq.1) write(51,*) "NB: Problem could be with the last source position/frequency."      
      write(51,*) "Aborting..."      
      write(51,*) '-----------------------------------------------------------------------'      
      close(51)
      
      write(*, *) '-------------------------------- ERROR --------------------------------'      
      write(*, *) "There was a problem reading the receiver position and/or data"      
      write(*, '(2(A,I8.5))') " for receiver number: ",indx-1,"   or ",indx      
      write(*, *) "Please make sure the x-, y-, and z-positions and/or data for "      
      write(*, '(2(A,I8.5))') " for receiver number: ",indx-1,"   or ", indx     
      write(*, *) "are specified correctly in fmm survey file: ",trim(tfile)      
      if (indx.eq.1) write(*, *) "NB: Problem could be with the last source position/frequency."      
      write(*, *) "Aborting..."      
      write(*, *) '-----------------------------------------------------------------------'      
      call crash_exit_fmm      

   case(58)
      open(51,file='fmm.log',status='old',action='write',position='append')
      write(51,*)
      write(51,*) "FMM: E4D and FMM cannot build mesh files together."
      write(51,*) "FMM: Please build mesh files either using E4D or FMM."
      write(51,*) "FMM: Aborting ..."
      write(*, *)
      write(*, *) "FMM: E4D and FMM cannot build mesh files together."
      write(*, *) "FMM: Please build mesh files either using E4D or FMM."
      write(*, *) "FMM: Aborting ..."
      close(51)
      call crash_exit_fmm

   case(121)
      open(51,file='fmm.log',status='old',action='write',position='append')
      if(mshfile(mnchar+2:mnchar+6) == ".node") then
         inquire(file=trim(mshfile),exist=exst)
         if(.not.exst) then
            write(51,*)
            write(*,*)
            write(51,*) "FMM: Cannot find the specified mesh node file: ",trim(mshfile)
            write(*, *) "FMM: Cannot find the specified mesh node file: ",trim(mshfile)
            close(51)
            call crash_exit_fmm
         end if
      elseif(mshfile(mnchar+2:mnchar+5) == ".ele") then
         inquire(file=trim(mshfile),exist=exst)
         if(.not.exst) then
            write(51,*)
            write(*,*)
            write(51,*) "FMM: Cannot find the specified mesh element file: ",trim(mshfile)
            write(*, *) "FMM: Cannot find the specified mesh element file: ",trim(mshfile)
            close(51)
            call crash_exit_fmm
         end if         
      else
         write(51,*)
         write(*,*)
         write(51,*) "FMM: If mode > 1 you must provide the name of the mesh with extension"
         write(51,*) "FMM: node file (.*.node) or mesh element file (.*.ele) where * is a single letter."
         write(51,*) "FMM: You provided: ",trim(mshfile)
         write(*, *) "FMM: If mode > 1 you must provide the name of the mesh with extension"
         write(*, *) "FMM: node file (.*.node) or mesh element file (.*.ele)  where * is a single letter."
         write(*, *) "FMM: You provided: ",trim(mshfile)
         close(51)
         call crash_exit_fmm
      end if      
      close(51)

   case(122)
      open(51,file='fmm.log',status='old',action='write',position='append')
      inquire(file=trim(mshfile(1:mnchar))//'1.node',exist=exst)
      if (.not.exst) then
         write(51,*) "FMM: Cannot find the mesh node file: ",trim(mshfile(1:mnchar))//'1.node'
         write(*, *) "FMM: Cannot find the mesh node file: ",trim(mshfile(1:mnchar))//'1.node'
         close(51)
         call crash_exit_fmm
      endif

      inquire(file=trim(mshfile(1:mnchar))//'1.ele',exist=exst)
      if (.not.exst) then
         write(51,*) "FMM: Cannot find the mesh element file: ",trim(mshfile(1:mnchar))//'1.ele'
         write(*, *) "FMM: Cannot find the mesh element file: ",trim(mshfile(1:mnchar))//'1.ele'
         close(51)
         call crash_exit_fmm
      endif

      inquire(file=trim(mshfile(1:mnchar))//'1.neigh',exist=exst)
      if (.not.exst) then
         write(51,*) "FMM: Cannot find the mesh neighbor file: ",trim(mshfile(1:mnchar))//'1.neigh'
         write(*, *) "FMM: Cannot find the mesh neighbor file: ",trim(mshfile(1:mnchar))//'1.neigh'
         close(51)
         call crash_exit_fmm
      endif

      inquire(file=trim(mshfile(1:mnchar))//'1.face',exist=exst)
      if (.not.exst) then
         write(51,*) "FMM: Cannot find the mesh face file: ",trim(mshfile(1:mnchar))//'1.face'
         write(*, *) "FMM: Cannot find the mesh face file: ",trim(mshfile(1:mnchar))//'1.face'
         close(51)
         call crash_exit_fmm
      endif

      inquire(file=trim(mshfile(1:mnchar))//'trn',exist=exst)
      if (.not.exst) then
         write(51,*) "FMM: Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn'
         write(*, *) "FMM: Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn'
         close(51)
         call crash_exit_fmm
      endif
      
   case DEFAULT

   end select
   
  end subroutine check_inp_fmm
  !_________________________________________________________________________
  
  !_________________________________________________________________________
  subroutine set_fresnel
    implicit none
    integer :: i1,i2

    open(10,file=invfile,status='old',action='read')   
    read(10,*,IOSTAT=ios) i1, i2 ; call check_inp_fmm(52,i2)
    fresnel = .false.
    if(ios.ne.0 .and. i2 .eq. 1) then
       fresnel = .true.
    end if

  end subroutine set_fresnel
  !_________________________________________________________________________

  !________________________________________________________________________
  subroutine read_survey_fmm
    implicit none
    logical :: exst
    integer :: i,j,junk,frflag
    real, dimension(3) :: etmp
    logical :: wdwarn = .true.

   
    inquire(file=trim(trim(tfile)),exist=exst)
    if(.not. exst) then
       call check_inp_fmm(24,0)
    else
       call check_inp_fmm(24,1)
    end if

    open(10,file=tfile,status='old',action='read')  

    !read in a ray based inversion file
   
    read(10,*,IOSTAT=ios) ns,junk;     call check_inp_fmm(12,junk)
    
    fresnel = .false.
    if(junk .eq. 1) fresnel = .true.
    if(.not. fresnel) then
       ! read source locations    
       allocate(s_pos(ns,3))       
       do i=1,ns
          read(10,*,IOSTAT=ios) junk,etmp; call check_inp_fmm(13,i)         
          if(junk>ns) call check_inp_fmm(14,i)
          s_pos(junk,1:3)=etmp
       end do
       ! read receiver locations    
       read(10,*,IOSTAT=ios) nrc;     call check_inp_fmm(112,junk)
       
       allocate(rc_pos(nrc,3))
       do i=1,nrc
          read(10,*,IOSTAT=ios) junk,etmp; call check_inp_fmm(113,i)
          if(junk>nrc) call check_inp_fmm(114,i)
          rc_pos(junk,1:3)=etmp
       end do
    else

       allocate(s_pos(ns,3),frq(ns))  
       
       do i=1,ns
          read(10,*,IOSTAT=ios) junk,etmp,frq(i); call check_inp_fmm(13,i)
          if (junk.ne.i)  call check_inp_fmm(56,i)
          if(frq(i).le.0) call check_inp_fmm(53,i)
          if(junk>ns) call check_inp_fmm(14,i)
          s_pos(junk,1:3)=etmp
       end do
       nrc = 0
    
    end if

    !!Read in the survey
    read(10,*,IOSTAT=ios) nm_fmm;     call check_inp_fmm(17,nm_fmm)
    allocate(dobs_fmm(nm_fmm),s_conf_fmm(nm_fmm,2),Wd_fmm(nm_fmm))
       do i=1,nm_fmm
          read(10,*,IOSTAT=ios) junk,s_conf_fmm(i,1:2),dobs_fmm(i),Wd_fmm(i); call check_inp_fmm(18,i)
          if (junk.ne.i) call check_inp_fmm(57,i)
          if(Wd_fmm(i) <= 0) then
             Wd_fmm(i)=1e15
             if( wdwarn) call check_inp_fmm(19,i)
             wdwarn = .false.
          endif

          Wd_fmm(i) = 1/Wd_fmm(i)
                    
          if(.not.fresnel) then
             do j=1,2
                !             if(s_conf_fmm(i,j)>ns .or. s_conf_fmm(i,j)<0) call check_inp_fmm(20,i)
                if(s_conf_fmm(i,1)>ns .or. s_conf_fmm(i,2) > nrc .or. s_conf_fmm(i,j)<0) call check_inp_fmm(20,i)
             end do
          else
             if(s_conf_fmm(i,1)>ns .or. s_conf_fmm(i,2)>ns) call check_inp_fmm(20,i)
          end if
       end do
       close(10)

       nm=nm_fmm
       allocate(dobs(nm),Wd(nm))
       dobs=dobs_fmm
       Wd=Wd_fmm
      
  end subroutine read_survey_fmm
  !_________________________________________________________________________

  !=========================================================================
  ! Subroutine to read time-lapse FMM data
  !=========================================================================
  subroutine get_dobs_tl_fmm(ind)
    implicit none

    integer, intent(in) :: ind

    ! local variables
    integer :: ns_test,i,j,nmt,junk,nrc_test
    real , dimension(3) :: etp
    real :: frq_test
    real :: perr=0.001
    real :: dobst,Wdt
    integer, dimension(2) :: ab

    open(11,file=trim(tl_dfils(ind)),status='old',action='read')
    read(11,*,IOSTAT=ios) ns_test,junk
    if(ios .ne. 0) goto 99
    if(ns_test .ne. ns) goto 100
    if (junk.ne.1.and.fresnel) goto 103

    if(.not. fresnel) then
       ! read source locations
       do i=1,ns
          read(11,*,IOSTAT=ios) junk,etp
          if(ios .ne. 0) goto 104
          if(abs(etp(1)-s_pos(i,1)-xorig)>perr .or. abs(etp(2)-s_pos(i,2)-yorig)>perr .or. abs(etp(3)-s_pos(i,3)-zorig)>perr) goto 101
       enddo

       ! read receiver locations    
       read(11,*,IOSTAT=ios) nrc_test
       if(nrc_test .ne. nrc) goto 105
       
       do i=1,nrc          
          read(11,*,IOSTAT=ios) junk,etp
          if(ios .ne. 0) goto 106
          if(abs(etp(1)-rc_pos(i,1)-xorig)>perr .or. abs(etp(2)-rc_pos(i,2)-yorig)>perr .or. abs(etp(3)-rc_pos(i,3)-zorig)>perr) goto 107
       enddo       
    else
       ! read source locations       
       do i=1,ns
          read(11,*,IOSTAT=ios) junk,etp,frq_test
          if(ios .ne. 0) goto 104
          if(abs(etp(1)-s_pos(i,1)-xorig)>perr .or. abs(etp(2)-s_pos(i,2)-yorig)>perr .or. abs(etp(3)-s_pos(i,3)-zorig)>perr) goto 101
          if (frq_test.ne.frq(i)) goto 108
       enddo       
    endif

    ! read the survey
    read(11,*,IOSTAT=ios) nmt  
    if(ios .ne. 0) goto 109
    if(nmt .ne. nm_fmm) goto 102
    
    do i=1,nm_fmm
       read(11,*,IOSTAT=ios) junk,ab,dobst,Wdt
       if(ios .ne. 0) goto 110
       if(Wdt.le.0) Wdt=1e9

       do j=1,2
          if(ab(j) .ne. s_conf_fmm(i,j)) goto 111
       enddo

       dobs_fmm(i) = dobst
       Wd_fmm(i)   = 1/Wdt
       
    enddo
            
    close(11)

    dobs=dobs_fmm
    Wd=Wd_fmm

    return
        
 99 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "There was an error reading number of source or fresnel option in: ",trim(tl_dfils(ind))
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

100 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The number of sources in: ",trim(tl_dfils(ind))," is: ",ns_test
    write(51,*) "The number of sources in the base survey is: ",ns
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

101 continue 
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The source position in: ",trim(tl_dfils(ind))," for source: ",i
    write(51,*) "is different than the same baseline source position."
    write(51,*) trim(tl_dfils(ind))," gives: ",etp
    write(51,*) "The baseline file gives: ",s_pos(i,1)-xorig,s_pos(i,2)-yorig,s_pos(i,3)-zorig
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

102 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The number of measurements in: ",trim(tl_dfils(ind))," is: ",nmt
    write(51,*) "The number of measurements in the base survey is: ",nm_fmm
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

103 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The fresnel option in: ",trim(tl_dfils(ind))," is not 1"
    write(51,*) "The fresnel option in the base survey is 1"
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

104 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "There was an error reading ",i,"th source position in: ",trim(tl_dfils(ind))
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

105 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The number of receivers in: ",trim(tl_dfils(ind))," is: ",nrc_test
    write(51,*) "The number of receivers in the base survey is: ",nrc
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

106 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "There was an error reading ",i,"th receiver position in: ",trim(tl_dfils(ind))
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

107 continue 
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The receiver position in: ",trim(tl_dfils(ind))," for receiver: ",i
    write(51,*) "is different than the same baseline receiver position."
    write(51,*) trim(tl_dfils(ind))," gives: ",etp
    write(51,*) "The baseline file gives: ",rc_pos(i,1)-xorig,rc_pos(i,2)-yorig,rc_pos(i,3)-zorig
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

108 continue 
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "The dominant frequency for source: ",i," in: ",trim(tl_dfils(ind))," is: ",frq_test
    write(51,*) "The dominant frequency for source: ",i," in the baseline survey is: ",frq(i)
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

109 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "There was an error reading number of measurements in: ",trim(tl_dfils(ind))
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

110 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "There was an error reading ",i,"th measurement in: ",trim(tl_dfils(ind))
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return

111 continue
    write(*,*) "FMM: INPUT ERROR: See fmm.log"  
    open(51,file='fmm.log',status='old',action='write',position='append')
    write(51,*) "A B for measurement: ",i," in: ",trim(tl_dfils(ind))," is: ",ab
    write(51,*) "The baseline A B is: ",s_conf_fmm(i,:)
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit_fmm
    return    

  end subroutine get_dobs_tl_fmm  
  !_________________________________________________________________________
  subroutine translate_source
    implicit none
    integer :: junk
    integer :: i
    mnchar = 0
    do i=1,40
       if(mshfile(i:i) == '.') then
          mnchar = i
          exit
       end if
    end do

    call check_inp_fmm(15,junk)
    open(21,file=mshfile(1:mnchar)//'trn',status='old')
    read(21,*,IOSTAT=ios) xorig,yorig,zorig; call check_inp_fmm(16,junk)
    close(21) 
    s_pos(:,1) = s_pos(:,1)-xorig
    s_pos(:,2) = s_pos(:,2)-yorig
    s_pos(:,3) = s_pos(:,3)-zorig
    if(.not. fresnel) then
       rc_pos(:,1) = rc_pos(:,1)-xorig
       rc_pos(:,2) = rc_pos(:,2)-yorig
       rc_pos(:,3) = rc_pos(:,3)-zorig
    end if
    
  end subroutine translate_source
  !_________________________________________________________________________

  !_________________________________________________________________________
  subroutine read_velocity
    implicit none
    integer :: i,junk,npre,nchr
    logical :: exst
    real :: tspd
    
    if(mode_fmm .ne. 1) then
       inquire(file=trim(trim(spdfile)),exist=exst)
       if(.not.exst) then
          read(spdfile,*,IOSTAT=ios) tspd
          if(ios .ne. 0 ) then
             call check_inp_fmm(25,junk)
          else
             nchr=len_trim(mshfile)
             do i=1,nchr
                if(mshfile(i:i)=='.') then
                   npre=i+1;
                   exit
                end if
             end do
             inquire(file=mshfile(1:npre)//".ele",exist=exst)
             if(.not.exst) then
                open(51,file='fmm.log',status='old',action='write')
                write(51,*)
                write(51,*) ' Cannot find the element file : ',mshfile(1:npre)//'.ele'
                close(51)
                write(*,*)
                write(*,*) ' Cannot find the ele file : ',mshfile(1:npre)//'.ele'
                close(51)
                call crash_exit_fmm
             else
                open(10,file=mshfile(1:npre)//".ele",status='old',action='read')
                read(10,*) nspd
                close(10)
                allocate(velocity(nspd))
                !slowness is 1/velocity**2 for the forward computations
                ! NB: velocity is defined here as slowness^2
                velocity=tspd**(-2)
                return
             end if
          end if
       end if
       open(10,file=spdfile,status='old',action='read')
       read(10,*,IOSTAT=ios) nspd; call check_inp_fmm(23,0) 
       
       allocate(velocity(nspd))
       
       do i=1,nspd
          read(10,*,IOSTAT=ios) velocity(i); call check_inp_fmm(23,i)
          !convert velocity to 1/(velocity^2) for travel time computations
          ! NB: velocity is defined here as slowness^2
          velocity(i) = velocity(i)**(-2)
       end do
       close(10)
    end if
       
  end subroutine read_velocity
  !_________________________________________________________________________
  subroutine check_fresnel
    implicit none

    integer :: i
    real    :: avgVel,length_sr
    real    :: FresVol,minFresVol,maxFresVol

    ! Start
    minFresVol = 1e15
    maxFresVol = 0.
    
    ! get average velocity
    avgVel = sum(1/sqrt(velocity))/nspd

    do i=1,nm_fmm
       ! get Fresnel volume
       length_sr = norm2( s_pos(s_conf_fmm(i,1),:) - s_pos(s_conf_fmm(i,2),:) )
       FresVol   = sqrt(avgVel*length_sr/frq(s_conf_fmm(i,1)))
       if (FresVol<minFresVol) minFresVol = FresVol
       if (FresVol>maxFresVol) maxFresVol = FresVol    
    enddo

    !open(51,file='fmm.log',status='old',action='write')
    !write(51,*)
    !write(51,*) " ======================================================================="
    !write(51,*) " FMM: Minimum Fresnel Volume: ",minFresVol
    !write(51,*) " FMM: Maximum Fresnel Volume: ",maxFresVol
    !write(51,*) " ======================================================================="
    !write(51,*)
    ! close(51)
    ! write(*, *)
    ! write(*, *) " ======================================================================="
    ! write(*, *) "  FMM: Minimum Fresnel Volume: ",minFresVol
    ! write(*, *) "  FMM: Maximum Fresnel Volume: ",maxFresVol
    ! write(*, *) " ======================================================================="
    ! write(*, *)
        
  end subroutine check_fresnel  
  !_________________________________________________________________________
  subroutine crash_exit_fmm
    implicit none
    integer, parameter :: errcode=100
    
    !call MPI_BCAST(0,1,MPI_INTEGER,0,FMM_COMM,ierr)
    call MPI_ABORT(MPI_COMM_WORLD,errcode,ierr)
    call PetscFinalize(perr)    
    stop
  end subroutine crash_exit_fmm
  !__________________________________________________________________________

  
end module input_fmm
