module input_fmm

  use fmm_vars
  use vars
  use input
  implicit none

  !Variables needed only by the master process are declared here
  !integer :: ios                                           !!io status
  !integer :: mnchar                                        !!used to record number of character in a string
  !real :: xorig,yorig,zorig                                !!x,y,z mesh translation values
  integer :: i_tl_fmm                                       !!current time lapse file index
  real, dimension(:) , allocatable :: tlt_fmm               !!time lapse data time markers
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
    read(smode,*,IOSTAT=ios) mode_fmm; call check_inp_fmm(1,junk)

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
    end if    
    
    ! read survey file  
    read(10,*,IOSTAT=ios) tfile;           call check_inp_fmm(3,junk)
    ! read velocity file
    read(10,*,IOSTAT=ios) spdfile;         call check_inp_fmm(4,junk)
    read(10,*,IOSTAT=ios) outfile_fmm;     call check_inp_fmm(5,junk)

    if(mode_fmm == 3) then
       read(10,*,IOSTAT=ios) invfile;      call check_inp(6,junk)
       read(10,*,IOSTAT=ios) refmod_file;  call check_inp(7,junk)
    end if
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


    !!check mesh files
    if (mode_fmm > 1)then 
       call check_inp_fmm(121,mnchar)
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
          write(51,*) "Cannot find the primary input file fmm.inp: aborting"
          write(*, *) "Cannot find the primary input file fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       end if

    case(101)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading the mode in fmm.inp: aborting"
          write(*, *) "There was a problem reading the mode in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       end if

    case(1)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) " There was a problem reading the mode in fmm.inp."
          write(51,*) " Aborting ..."
          write(*,*) " There was a problem reading the mode in fmm.inp."
          write(*,*) " Aborting ..."
          close(51)
          call crash_exit_fmm
       end if
       mcheck = .false.
       select case(mode_fmm)
       case(0) 
       case(1) 
       case(2) 
          mcheck = .true.
       case(3)
          mcheck = .true.
       end select

       if(.not.mcheck) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,"(A33,I3,A19)") "The mode selected in fmm.inp is ",mode_fmm,", which is invalid."
          write(51,*)"Valid run modes include: < see below > "
          write(51,*)  "---------------------------------------"
          write(51,*)  " -FUNCTION-                      -MODE-"
          write(51,*)  "---------------------------------------"
          write(51,*)  " FMM Forward                       2"
          write(51,*)  " FMM Inversion                     3"
          write(51,*)  "---------------------------------------"
          write(51,*)  " Aborting ..."
          write(*,"(A33,I3,A19)") "The mode selected in fmm.inp is ",mode_fmm,", which is invalid."
          write(*, *)"Valid run modes include: < see below > "
          write(*, *)  "---------------------------------------"
          write(*, *)  " -FUNCTION-                      -MODE-"
          write(*, *)  "---------------------------------------"
          write(*, *)  " FMM Forward                       2"
          write(*, *)  " FMM Inversion                     3"
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
                write(51,*) "***************** RUNNING IN FMM FORWARD MODE ******************"
                write(51,"(A,I3.3)") "  Mode:                             ",mode_fmm
             case(3)
             	write(51,*) "***************** RUNNING IN FMM INVERSE MODE *******************"
             	write(51,"(A,I3.3)") "  Mode:                             ",mode_fmm
             end select
          close(51)
       end if

    case(2)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then

             open(51,file='fmm.log',status='old',action='write',position='append')
             write(51,*) "  ERROR: There was a problem reading the mesh file name in fmm.inp"
             write(51,*) "  Aborting ..."
             close(51)
             write(*,*) "  ERROR: There was a problem reading the mesh file name in fmm.inp"
             write(*,*) "  Aborting ..."
          call crash_exit_fmm

       else
             write(51,*) " Mesh file:                        ",trim(mshfile)
       end if
       close(51)

    case(3)

       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the survey file name in fmm.inp: aborting"
          write(*,*) "There was a problem reading the survey file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Survey configuration file:        ",trim(tfile)
       end if
       close(51)

    case(4)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the velocity file name in fmm.inp: aborting"
          write(*,*) "There was a problem reading the velocity file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Velocity file:                       ",trim(spdfile)
       end if
       close(51)

    case(5)
      
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the output options file name in fmm.inp: aborting"
          write(*,*) "There was a problem reading the output options file name in fmm.inp: aborting"
          close(51)
          call crash_exit_fmm
       else
          write(51,*) " Output options file:              ",trim(outfile_fmm)
       end if
       close(51)

    case(6)

    case(7)

    case(8)


    case(9)

    case(10)
    case(11)
      
    case(12)

       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) 
          write(51,*) " There was a problem reading the number of sources in the survey file ",trim(tfile)
          write(51,*) " Aborting ..."
          write(*,*) 
          write(*,*) " There was a problem reading the number of sources in the survey file ",trim(tfile)
          write(*,*) " Aborting ..."
          
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
          write(51,*) " in the survey file file: ",trim(tfile)
          write(51,*) " Aborting ..."
          if(fresnel) then
             write(51,*) "Running in fresnel mode ... be sure to include the positive"
             write(51,*) "frequency value in the last column"
          end if
          close(51)
          write(*,*)
          write(*,"(A,I7.7)") " There was a problem reading source number ",indx
          write(*,*) "in the survey file file: ",trim(tfile)
          write(*,*) "Aborting ..."
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
       write(51,*) " The source index specified for source: ",indx
       write(51,*) " is greater than the total number of sources."
       write(51,*) " Aborting ..."
       write(*,*) 
       write(*,*) " The source index specified for source: ",indx
       write(*,*) " is greater than the total number of sources."
       write(*,*) " Aborting ..."
     
       close(51)
       
       call crash_exit_fmm

    case(112)

       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*) 
          write(51,*) " There was a problem reading the number of receivers in the survey file ",trim(tfile)
          write(51,*) " Aborting ..."
          write(*,*) 
          write(*,*) " There was a problem reading the number of receivers in the survey file ",trim(tfile)
          write(*,*) " Aborting ..."
          
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
       write(51,*) " The receiver index specified for receiver: ",indx
       write(51,*) " is greater than the total number of receivers."
       write(51,*) " Aborting ..."
       write(*,*) 
       write(*,*) " The receiver index specified for receiver: ",indx
       write(*,*) " is greater than the total number of receivers."
       write(*,*) " Aborting ..."
     
       close(51)
       
       call crash_exit_fmm
    case(15)
       inquire(file=mshfile(1:mnchar)//'trn',exist=exst)
       if(.not. exst) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn' 
          write(51,*) " Aborting..."
          write(*,*)
          write(*,*) " Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn' 
          write(*,*) " Aborting..."
          close(51)
          call crash_exit_fmm
       end if

    case(16)
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " There was a problem reading the mesh "
          write(51,*) " translation numbers in: ",trim(mshfile(1:mnchar))//'trn' 
          write(51,*) " Aborting ... "
          close(51)
          write(*,*)
          write(*,*) " There was a problem reading the mesh "
          write(*,*) " translation numbers in: ",trim(mshfile(1:mnchar))//'trn' 
          write(*,*) " Aborting ... "
          close(51)
          call crash_exit_fmm
       end if

    case(17)
       open(51,file='fmm.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='fmm.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " There was a problem reading the number of measurements" 
          write(51,*) " in the survey file: ",trim(tfile)
          close(51)
          write(*,*)
          write(*,*) " There was a problem reading the number of measurements" 
          write(*,*) " in the survey file: ",trim(tfile)
          close(51)
          call crash_exit_fmm
       elseif(nm_fmm .le. 0) then
          write(51,*) " The number of measurements is not positive: ",indx
          write(51,*) " Aborting ..."
          write(*,*)
          write(*,*) " The number of measurements is not positive: ",indx
          write(*,*) " Aborting ..."
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
         write(51,*) " In mode > 1 you must provide a node or element file name"
         write(51,*) " with .*.node or .*.ele file name extension where * is a single letter. "
         write(51,*) " You provided: ",trim(mshfile)
         write(51,*) " Aborting ..."
         write(*,*)
         write(*,*) " In mode > 1 you must provide a node or element file name"
         write(*,*) " with .*.node or .*.ele file name extension where * is a single letter. "
         write(*,*) " You provided: ",trim(mshfile)
         write(*,*) " Aborting ..."
      else if(indx==1) then
         write(51,*) 
         write(51,*) " Cannot find the mesh configuration file: ",trim(mshfile)
         write(51,*) " Aborting ..."
         write(*,*) 
         write(*,*) " Cannot find the mesh configuration file: ",trim(mshfile)
         write(*,*) " Aborting ..."
      end if
      close(51)
      call crash_exit_fmm
      
   case(22)
      if(ios .ne. 0) then
         open(51,file='fmm.log',status='old',action='write',position='append')
         write(51,*) "There was a problem reading the number of velocities in"
         write(51,*) "in the velocity file: ",trim(spdfile)
         write(51,*) "aborting."
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
            write(51,*) " There was a problem reading the number of "
            write(51,*) " velocity values in ",trim(spdfile)
            write(51,*) " Aborting ..."
            close(51)
            write(*,*) 
            write(*,*) " There was a problem reading the number of "
            write(*,*) " velocity values in ",trim(spdfile)
            write(*,*) " Aborting ..."
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
        write(51,*) " Can't find the survey file: ",trim(tfile)
        write(51,*) " aborting."
        write(*,*)
        write(*,*) " Can't find the survey file: ",trim(tfile)
        write(*,*) " aborting."
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
     write(51,*) " Can't find the velocity file: ",trim(spdfile)
     write(51,*) " Aborting ..."
     close(51)
     write(*,*) 
     write(*,*) " Can't find the velocity file: ",trim(spdfile)
     write(*,*) " Aborting ..."
     close(51)
     call crash_exit_fmm
      

     case(121)
        open(51,file='fmm.log',status='old',action='write',position='append')
        if(mshfile(mnchar+2:mnchar+6) == ".node") then
           inquire(file=trim(mshfile),exist=exst)
           if(.not.exst) then
              write(51,*)
              write(*,*)
              write(51,*) " Cannot find the specified mesh node file: ",trim(mshfile)
              write(*,*) " Cannot find the specified mesh node file: ",trim(mshfile)
              close(51)
              call crash_exit_fmm
           end if
        elseif(mshfile(mnchar+2:mnchar+5) == ".ele") then
           inquire(file=trim(mshfile),exist=exst)
           if(.not.exst) then
              write(51,*)
              write(*,*)
              write(51,*) " Cannot find the specified mesh element file: ",trim(mshfile)
              write(*,*) " Cannot find the specified mesh element file: ",trim(mshfile)
              close(51)
              call crash_exit_fmm
           end if
        else
           write(51,*)
           write(*,*)
           write(51,*) " If mode > 1 you must provide the name of the mesh with extension"
           write(51,*) " node file (.*.node) or mesh element file (.*.ele) where * is a single letter."
           write(51,*) " You provided: ",trim(mshfile)
           write(*,*) " If mode > 1 you must provide the name of the mesh with extension"
           write(*,*) " node file (.*.node) or mesh element file (.*.ele)  where * is a single letter."
           write(*,*) " You provided: ",trim(mshfile)
           close(51)
           call crash_exit_fmm
        end if
        close(51)
      
      case(52)
         if(ios.ne.0 .or. indx .ne. 0 .or. indx .ne. 1) then
            open(51,file='fmm.log',status='old',action='write',position='append')
            write(51,*) "There was a problem reading the first line of the "
            write(51,*) "fmm survey file ",trim(tfile)
            write(51,*) "In fmm mode, the first line of the survey file"
            write(51,*) "should contain two integers: the number of source positions "
            write(51,*) "and the fresnel volume flag"
            write(51,*) "0 for ray-based and 1 for fresnel volume"
            write(51,*) "Using ray-based by default "
            write(51,*) "Aborting ..."
            close(51)
            write(*,*)
            write(*,*) "There was a problem reading the first line of the "
            write(*,*) "fmm survey file ",trim(tfile)
            write(*,*) "In fmm mode, the first line of the survey file"
            write(*,*) "should contain two integers: the number of source positions "
            write(*,*) "and the fresnel volume flag"
            write(*,*) "0 for ray-based and 1 for fresnel volume"
            write(*,*) "Aborting ..."
            call crash_exit_fmm
         end if

         case(53)
             open(51,file='fmm.log',status='old',action='write',position='append')
             write(51,*) "The frequency for source number: ",indx
             write(51,*) "Is less than or equal to zero. Specified frequencies"
             write(51,*) "must be positive."
             write(51,*) "Aborting..."
             close(51)
             write(*,*) "The frequency for source number: ",indx
             write(*,*) "Is less than or equal to zero. Specified frequencies"
             write(*,*) "must be positive."
             write(*,*) "Aborting..."
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
          else
             Wd_fmm(i) = 1/Wd_fmm(i)
          end if
          
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
  !_________________________________________________________________________
  subroutine crash_exit_fmm
    
    call MPI_BCAST(0,1,MPI_INTEGER,0,FMM_COMM,ierr)
    call PetscFinalize(perr)
    stop
  end subroutine crash_exit_fmm
  !__________________________________________________________________________

  
end module input_fmm
