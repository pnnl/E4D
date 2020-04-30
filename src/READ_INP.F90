module input
  !_________________________________________________________________
  !! Author: Tim Johnson
  !!
  !! email: tj@pnnl.gov
  !!
  !! Purpose: This module contains subroutines to read user inputs and
  !! check their validity as possible. This module is used only by the   
  !! master process, and variable declared herein are used only by the 
  !! master process.
  !!
  !! Revision History: Last modified to include documentation by TJ 
  !! on 06/05/2014.
  !!
  !__________________________________________________________________

  use vars
  use distributor
  
  implicit none

  !Variables needed only by the master process are declared here

  integer :: up_opt                                        !!inversion update option
  integer :: nlsp                                          !!number of line search scaling paramters
  integer :: cull_flag                                     !!flag to ignore data outliers
  integer :: i_tl                                          !!current time lapse file index
  integer :: mnchar                                        !!used to record number of character in a string
  integer :: ntl                                           !!number of time-lapse data files
  
  real :: wtx,wty,wtz                                      !!x,y,z anisotropic weighting factors
  real :: red_fac                                          !!beta reduction factor
  real :: xorig,yorig,zorig                                !!x,y,z mesh translation values
  real :: cull_dev                                         !!data culling cutoff (standard deviation)

  character*40 :: refmod_file                              !!reference model file name
  character*40 :: gs_file                                  !!geostatistical constraint file name
  character*40 :: tl_file                                  !!time-lapse survey list file
  character*40, dimension(:), allocatable :: tl_dfils      !!time-lapse survey file names
  character*40, dimension(:,:), allocatable :: tl_cfils    !!conductivity file list
  character, pointer :: cptr

  real :: beta, beta_red, beta_s                           !!starting reg wt. reg.reduction factor
  real :: norm_chi2                                        !!target chi-squared value
  real :: del_obj                                          !!minimum change in objective function
  real :: delta_initer                                     !!minimum change at inner iteration
  real :: max_sig,min_sig                                  !!max/min sigma value
  real :: initer_conv                                      !!inner iteration convergence value

  integer, dimension(:), allocatable :: reg_opt            !!zonal regularization options
  integer, dimension(:,:), allocatable :: smetric          !!model structure zone and metric option
  integer, dimension(:), allocatable :: Fw_type            !!Weighting function type
  integer, dimension(:,:), allocatable :: zone_links       !!inter-zone regularization links
  integer, dimension(:,:), allocatable :: rblock           !!flag for which reg block to use
  integer, dimension(:), allocatable :: node_map           !!for mesh reording with inactive cells
  integer, dimension(:), allocatable :: element_map        !!for mesh reording with inactive cells
 
  integer :: ios                                           !!io status
  logical :: gs_flag = .false.                             !!semivariogram constraints flag
  logical :: tl_ly = .false.                               !!flag to execute Labrecque-Yang time-lapse
  logical :: i_flag = .false.                              !!flag for complex conductivity inversion
  logical :: r_last = .false.                              !!flag for which baseline to ref in t_lapse mode
  logical :: rt_flag                                       !!flag for real time imaging
  logical :: multi_forward = .false.                       !!true for multi conductivity field forward run
  logical :: cfg_flag = .false.                            !!flag for whether a file is a .cfg or not for testing of files
  logical :: inv_opt_flag = .true.                        !!flag for whether a inverse options file is used for testing of files
  logical :: refmod_flag = .false.                         !!flag for whether a reference model is used within inverse options for testing of files
  logical :: hemstiching = .false.                         !!flag indicating hemstiching in inversion
  logical :: diverging  = .false.                          !!inversion divergence flag

  real, dimension(:,:), allocatable :: zwts                !!zonal regularization weights
  real, dimension(:,:) , allocatable :: ex_vals            !!log of external constraints
  real, dimension(:) , allocatable :: lsp                  !!line search scaling values
  real, dimension(:) , allocatable :: Wd_cull              !!Culled data weights            
  real, dimension(:) , allocatable :: gwts                 !!geocon weights 
  real, dimension(:) , allocatable :: tlt                  !!time lapse data time markers
  real, dimension(:,:) , allocatable :: Fw_parm            !!Weighting function parameters
  real, dimension(:) , allocatable :: C_targ               !!constraint targets
  real, dimension(:) , allocatable :: phase                !!constraint targets

contains

  !____________________________________________________________________________________________
  subroutine read_input
    implicit none
    
    character*40 :: smode
    character*40 :: chk_sigfile                              !!temp variable to check conductivity file value in e4d.inp

    integer :: nchar,junk,i,check,j
    logical :: exst
    logical :: wdwarn = .false.


    call check_inp(0,junk)
    open(10,file='e4d.inp',status='old',action='read') 
    read(10,*,IOSTAT=ios) smode; call check_inp(101,junk)
    
    ! check for alpha characters and make them all lower case
    smode=lcase(smode)
    
    if(trim(adjustl(smode))=='analytic' .or. trim(adjustl(smode))=='ertanalytic') then
       analytic = .true.
       mode = 2; call check_inp(1,junk)
    else    
		 select case(trim(adjustl(smode)))
			   case('ert0') 
				  mode=0; call check_inp(1,junk)
			   case('ert1') 
				  mode=1; call check_inp(1,junk)
			   case('ert2') 
				  mode=2; call check_inp(1,junk)
			   case('ert3') 
				  mode=3; call check_inp(1,junk)
			   case('ert4') 
				  mode=4; call check_inp(1,junk)
			   case('ert5') 
				  mode=5; call check_inp(1,junk)
		
			   case('sip0') 
				  mode=0; call check_inp(1,junk)
				  i_flag = .true.
			   case('sip1') 
				  mode=21; call check_inp(1,junk)
			   case('sip2') 
				  mode=22; call check_inp(1,junk)
			   case('sip3') 
				  mode=23; call check_inp(1,junk)
			   case('sip4') 
				  mode=24; call check_inp(1,junk)
		
			   case('erttank0') 
				  mode=0; call check_inp(1,junk)
				  tank_flag = .true.
			   case('erttank1') 
				  mode=31; call check_inp(1,junk)
			   case('erttank2') 
				  mode=32; call check_inp(1,junk)
			   case('erttank3') 
				  mode=33; call check_inp(1,junk)
			   case('erttank4')
				  mode=34; call check_inp(1,junk)
			   case('erttank5')
				  mode=35; call check_inp(1,junk)
		
			   case('siptank0') 
				  mode=0; call check_inp(1,junk)
				  tank_flag = .true.
				  i_flag=.true.
			   case('siptank1') 
				  mode=41; call check_inp(1,junk)
			   case('siptank2') 
				  mode=42; call check_inp(1,junk)
			   case('siptank3') 
				  mode=43; call check_inp(1,junk)
			   case('siptank4') 
				  mode=44; call check_inp(1,junk)	   
			   case DEFAULT
				  read(smode,*,IOSTAT=ios) mode; call check_inp(1,junk)
		   end select
    end if
    
    !! if running in mode zero attempt the read the input files, report, and exit
    if(mode == 0) then
       !call test_read
       return
    end if
    
    read(10,*,IOSTAT=ios) mshfile;  call check_inp(2,junk)
  
    if(mode>1 .and. mode .ne. 21 .and. mode .ne. 31 .and. mode .ne.41) then
       read(10,*,IOSTAT=ios) efile;    call check_inp(3,junk)
       read(10,*,IOSTAT=ios) sigfile;  call check_inp(4,junk)
       
       ! check for alpha characters and make them all lower case
       chk_sigfile=sigfile
       chk_sigfile=lcase(chk_sigfile)
       if(trim(chk_sigfile)=="average") then
          ave_sig = .true.
       end if
       read(10,*,IOSTAT=ios) outfile;             call check_inp(5,junk)
    end if
    if(analytic) return

    if(mode==5 .or. mode==35) then
       rt_flag = .true.
       if(mode==5) mode = 4
       if(mode==35) mode = 34
    end if

    if(mode == 3 .or. mode == 4 &
         .or. mode == 33 &
         .or. mode == 34 &
         .or. mode .ge. 100) then
       
       read(10,*,IOSTAT=ios) invfile;      call check_inp(6,junk)
       read(10,*,IOSTAT=ios) refmod_file;  call check_inp(7,junk)      
    
    elseif(mode == 23 .or. mode == 43 .or. mode==24 .or. mode==44) then

       read(10,*,IOSTAT=ios) invfile, iinvfile;      call check_inp(6,junk)
       read(10,*,IOSTAT=ios) refmod_file;           call check_inp(7,junk)     

    end if
    

    if(mode==4 ) tl_ly = .true.

    if(mode==21) then
       mode=1
       i_flag=.true.
    elseif(mode==22) then
       mode=2
       i_flag = .true.
    elseif(mode==23 ) then
       mode=3
       i_flag = .true.
    elseif(mode==24 ) then
       i_flag = .true.
       tl_ly = .true.
    end if
    
    if( (mode.ge.31 .and. mode.le.34) .or. (mode .ge. 41 .and. mode .le. 44) ) then 
       tank_flag = .true.
       select case(mode)
          case(31)
             mode = 1
          case(32)
             mode = 2
          case(33)
             mode = 3
          case(34)
             tl_ly = .true.
          case(41)
             mode = 1
             i_flag = .true.
          case(42)
             mode = 2
             i_flag = .true.
          case(43)
             mode = 3
             i_flag = .true.
          case(44)
             i_flag = .true.
             tl_ly = .true.
       end select
   
    end if

    if(mode .ge.  100) then
       res_flag = .true.
       if(mode == 101) opt_flag = .true. 
       if(mode == 102) psf_flag = .true. 
    end if
      

    !!read the time lapse data list file
    if(tl_ly) then

       mode = 3
       check = 0

       read(10,*,IOSTAT=ios) tl_file,check;   call check_inp(8,check)
       if(check==2) r_last = .true.
      
       if(rt_flag) then
          ntl=1 
          allocate(tl_dfils(ntl),tlt(ntl))   
       else

          open(21,file=trim(tl_file),status='old',action='read',IOSTAT=ios)  
          read(21,*,IOSTAT=ios) ntl; call check_inp(9,junk)
          allocate(tl_dfils(ntl),tlt(ntl))
          
          do i=1,ntl
             read(21,*,IOSTAT=ios) tl_dfils(i), tlt(i); call check_inp(10,i)
             inquire(file=trim(tl_dfils(i)),exist=exst)
             if(.not. exst) call check_inp(11,i)
          end do
          close(21)
       end if

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
    if(mnchar == 0) call check_inp(21,0)

    !!Check for compatibility between mshfile and mode
    !!Check for config file
    if (mode==1) then
       if(mshfile(mnchar+1:mnchar+3) == "cfg") then
          inquire(file=trim(mshfile),exist=exst)
          if(.not.exst) call check_inp(21,1)
       else if(mshfile(mnchar+1:mnchar+3).ne."cfg") then
          call check_inp(21,0)
       end if
    end if
    
    if(mode == 1) return

    !!check mesh files
    if (mode > 1)then 
       call check_inp(121,mnchar)
    end if
    
    !!Allocate/read the electrode positions and survey configuration
    if(mode>1) then
       call read_survey
       call translate_electrodes
    end if
    
    !!if we're running in a forward mode, check to see if a list
    !!file was provided. If so, check the list file to make sure 
    !!all of the files exist and set the multi_forward flag to 
    !!true. 
    if(.not. ave_sig) then
       if(mode==2 .or. mode==22 .or. mode==32 .or. mode==42 .or. mode==52) then
          multi_forward = .false.
          call check_for_list
       end if
    end if

    !!read the conductivity file
    if(.not. ave_sig .and. .not. multi_forward) then
       call read_conductivity
    end if

  end subroutine read_input
  !___________________________________________________________________________________

  !__________________________________________________________________________________
  subroutine get_inv_optsII  
    !_________________________________________________________________
    !! Author: Tim Johnson
    !
    !! email: tj@pnnl.gov
    !
    !! Purpose: Reads the inversion options and checks their validity as possible   
    !!  
    !!
    !! USE variables accessed or modified:
    !! vars:
    !! nrz :                  number of regularization blocks
    !!
    !!
    !! input:
    !! invfile :             inversion options file name
    !__________________________________________________________________

    implicit none
    integer :: i,j,k                                      !!counters
    integer :: max_nlink                                  !!maximum number of zone
                                                          !!links in any regularization block
    integer :: izn                                        !!zone number,
    integer :: nstmp,junk                                 !!temporary placeholders 
    integer :: nzn                                        !!number of zones defined in mesh

    logical :: exst                                       !!for file exist check
    logical :: EXTRNL=.false.                           !!flag for external constraint file

    character*80 :: fnam,str1,str2                        !!input strings
    real :: rdum,rsig,isig                                !!placeholder

    integer, dimension(:), allocatable :: itmp,zflag
    
   
    !!Read the inversion options file and check for errors
    call check_inv_opts(0,junk)
    if(invi) then
       open(10,file=iinvfile,status='old',action='read')
    else
       open(10,file=invfile,status='old',action='read')
    end if

    
    if(allocated(smetric))    deallocate(smetric)
    if(allocated(Fw_type))    deallocate(Fw_type)
    if(allocated(Fw_parm))    deallocate(Fw_parm)
    if(allocated(C_targ))     deallocate(C_targ)
    if(allocated(zwts))       deallocate(zwts)
    if(allocated(itmp))       deallocate(itmp)
    if(allocated(zone_links)) deallocate(zone_links)
    if(allocated(zflag))      deallocate(zflag)
    if(allocated(lsp))        deallocate(Wd_cull)
    if(allocated(Wd_cull))    deallocate(Wd_cull)
  
    nrz=0
    read(10,*,IOSTAT=ios) nrz ;              call check_inv_opts(1,junk) 

    allocate(smetric(nrz,3),Fw_type(nrz),Fw_parm(nrz,2),C_targ(nrz),zwts(nrz,4))
    allocate(itmp(nrz))

    max_nlink=0
    do i=1,nrz
       !!Do and initial check on each constraint block specified
       !!in the inverse options file 

       !check for an external constraint file or zone number to constrain for this block
       EXTRNL=.false.
       read(10,*,IOSTAT=ios) str1;            call check_inv_opts(2,i)
       if(trim(str1).ne.'EXTERNAL'.and.trim(str1).ne.'External'.and.trim(str1).ne.'external') then
          read(str1,*,IOSTAT=ios) izn 
          call check_inv_opts(2,i)
       else 
          EXTRNL=.true.
       end if

       read(10,*,IOSTAT=ios) izn ;            call check_inv_opts(3,i)      
       read(10,*,IOSTAT=ios) izn ;            call check_inv_opts(4,i)

       if(EXTRNL) then
          read(10,*,IOSTAT=ios) str1;         call check_inv_opts(33,i)
       else
          read(10,*,IOSTAT=ios) itmp(i);         call check_inv_opts(5,i)
          if(itmp(i)>max_nlink) then
             max_nlink=itmp(i);     
          end if
       end if
       read(10,*,IOSTAT=ios) str1;            call check_inv_opts(6,i)
       read(10,*,IOSTAT=ios) rdum;            call check_inv_opts(7,i)
    end do

    !reset the external file list
    open(13,file='external_files.tmp',status='replace',action='write') 
    close(13)

    rewind(10)
    allocate(zone_links(nrz,max_nlink+1))
    zone_links = 0
    nzn=int(maxval(zones))
    zwts = 0

    read(10,*,IOSTAT=ios) nrz;                call check_inv_opts(8,junk)
    smetric=0
    C_targ=0
   
    do i=1,nrz
       EXTRNL = .false.
       call check_inv_opts(9,i)

       !check for external
       read(10,*,IOSTAT=ios) str1;      
       if(trim(str1).ne.'EXTRNL'.and.trim(str1).ne.'External'.and.trim(str1).ne.'external') then
          read(str1,*,IOSTAT=ios) smetric(i,1)
          call check_inv_opts(10,i)
       else 
          !!smetric(i,1) is the zone this regularization block operates one
          !!for external constraints that zone is specified as nzn+1 to distinguish
          !!from any zone that exists in the mesh
          EXTRNL=.true.
          smetric(i,1) = nzn+1
          call check_inv_opts(10,-1)
       end if

       read(10,*) smetric(i,2),zwts(i,2:4) ;  call check_inv_opts(11,i)

       if(smetric(i,2) .eq. 12 .or. smetric(i,2) .eq.13) then
          !structural metric 12 and 13 specify cross gradient
          !joint inversion for velocity and conductivity so 
          !both e4d and fmm must be running
          if(.not.simulate_e4d .or. .not.simulate_fmm) then
             call check_inv_opts(32,smetric(i,2))
          end if
          cgmin_flag(1) = .true.
       end if

       read(10,*) Fw_type(i),Fw_parm(i,1:2);  call check_inv_opts(12,i)
       if(EXTRNL) then
          read(10,*,IOSTAT=ios) str1;         call check_inv_opts(33,i)
          inquire(file=trim(str1),EXIST=exst)
          if(.not.exst) then
             call check_inv_opts(34,i)
          else
             tmpstr=str1
             call check_inv_opts(35,i)
             open(11,file='external_files.tmp',status='old',action='write',position='append')
             write(11,*) i,trim(str1)
             close(11)
          end if     
       else
          read(10,*) zone_links(i,1:itmp(i)+1);  call check_inv_opts(13,i)
       end if
   
       read(10,*,IOSTAT=ios) str1
       if(str1=="ref" .or. str1=="REF" .or. str1=="Ref") then
          smetric(i,3)=1;                     call check_inv_opts(14,i)
       elseif(str1=="pref" .or. str1=="PREF" .or. str1=="Pref") then
          if(.not. tl_ly) call check_inv_opts(30,i)
          smetric(i,3)=2;                     call check_inv_opts(14,i)
          if(.not.allocated(prefsig)) then
             allocate(prefsig(nelem))
             if(invi) then
                prefsig = phase
             else
                prefsig = sigma
             end if
          end if
       else
          read(str1,*,IOSTAT=ios) C_targ(i);  call check_inv_opts(14,i)
          if(smetric(i,2).eq.3 .or. smetric(i,2).eq.4 .or. smetric(i,2).eq.7 .or. smetric(i,2).eq.8) then
             if(C_targ(i) .le. 0)                call check_inv_opts(31,i)
             C_targ(i) = log(C_targ(i))
          end if
       end if
                                           
       read(10,*,IOSTAT=ios) zwts(i,1)  ;     call check_inv_opts(15,i)
                                              call check_inv_opts(16,i)

       rdum=sqrt(zwts(i,2)**2 + zwts(i,3)**2 + zwts(i,4)**2)
       zwts(i,2:4)=zwts(i,2:4)/rdum

    end do
    
    allocate(zflag(nzn))
    zflag=0
    do i=1,nrz
       if(smetric(i,1).le.nzn) then
          zflag(smetric(i,1))=1
       end if
    end do
    do i=1,nzn
       if(zflag(i)==0) then
          call check_inv_opts(17,i)
       end if
    end do
    
    !!check for reference model
    if(allocated(refsig)) deallocate(refsig)
    allocate(refsig(nelem))
    refsig=1.0
    do i=1,nrz
       if(smetric(i,3)==1) then
          refmod_flag=.true. 
          call check_inv_opts(18,i)
          open(11,file=trim(refmod_file),status='old',action='read')

          if(allocated(element_map)) then
             k=0
             read(11,*,IOSTAT=ios) nstmp
             if(invi) then
            
                do j=1,nstmp
                   read(11,*,IOSTAT=ios) rsig,isig
                   if(ios.ne.0) then
                      call check_inv_opts(20,j)
                   end if
                   if(element_map(j) .ne. 0) then
                      k=k+1
                      refsig(element_map(j))= atan(isig/rsig)   
                   end if
                end do
                
             else
          
                do j=1,nstmp
                   read(11,*,IOSTAT=ios) rsig
                   if(ios.ne.0) then
                      call check_inv_opts(20,j)
                   end if
                   if(element_map(j) .ne. 0 ) then
                      k=k+1
                      refsig(element_map(j)) = rsig
                   end if
                end do
                
             end if
            
             call check_inv_opts(19,k)
          else
             read(11,*,IOSTAT=ios) nstmp ;      
             if(invi) then
                k=0
                do j=1,nstmp
                   read(11,*,IOSTAT=ios) rsig,isig
                   if(ios.ne.0) then
                      call check_inv_opts(20,j)
                   end if
                   refsig(j)= atan(isig/rsig)                
                end do
                
             else
                
                do j=1,nstmp
                   read(11,*,IOSTAT=ios) refsig(j)
                   if(ios.ne.0) then
                      call check_inv_opts(20,j)
                   end if
                end do
                
             end if
             if(im_fmm) refsig = 1/refsig
          end if
          close(11)
          exit
       end if
    end do
    
    if (.not. refmod_flag) then 
       call check_inv_opts(36,junk)
    end if
    
    read(10,*,IOSTAT=ios) beta, del_obj,beta_red;     call check_inv_opts(21,junk)  
    beta_s = beta
    conv_opt = 1
  
    read(10,*,IOSTAT=ios) norm_chi2 ;                 call check_inv_opts(22,junk)
    read(10,*,IOSTAT=ios) min_initer, max_initer ;    call check_inv_opts(23,junk)
    delta_initer=0.001
    initer_conv = 1e-4;

    read(10,*,IOSTAT=ios) min_sig, max_sig ;          call check_inv_opts(24,junk)

    read(10,*,IOSTAT=ios) up_opt  ;                   call check_inv_opts(25,junk)
    if(up_opt == 1) then
       read(10,*,IOSTAT=ios) nlsp  ;                  call check_inv_opts(26,junk)
       nlsp=nlsp+1
       allocate(lsp(nlsp))
       read(10,*) lsp(2:nlsp)     ;                   call check_inv_opts(27,junk)
       lsp(1) = 0
    end if
   
    !default to no line search for now until further testing: tcj-Aug 7, 2014
    select case (up_opt)
       case(1)
          up_opt = 2
       case(3)
          conv_opt = 2
          up_opt = 2
       case DEFAULT
          up_opt = 2
       end select
   
    
    read(10,*,IOSTAT=ios) cull_flag, cull_dev ;       call check_inv_opts(28,junk)
    allocate(Wd_cull(nm))
    Wd_cull = 1.0
    nec = 0
    
    !call check_inv_opts(29,junk)

    close(10)
    return

  end subroutine get_inv_optsII
  !___________________________________________________________________________________  

 
  !__________________________________________________________________________________
  subroutine get_inv_opts
    !!reads in the inversion regularization and weighting options
    implicit none
    integer :: i,j,k,max_nlink,izn,nstmp
    logical :: exst
    character*40 :: fnam
    
    integer, dimension(:), allocatable :: itmp
    open(51,file='e4d.log',status='old',action='write',position='append')
    inquire(file=trim(invfile),exist=exst)
    if(.not.exst) then
       write(51,1001) trim(invfile)
       close(51)
       stop
    end if
    
    open(10,file=invfile,status='old',action='read')
    read(10,*) nrz
    write(51,*)
    write(51,*) "NUMBER OF REGULARIZED ZONES", nrz
    close(51)
    
    open(51,file='e4d.log',status='old',action='write',position='append')
    
    allocate(reg_opt(nrz),zwts(nrz,4),itmp(nrz))
    reg_opt = -1

    !!find the maximum number of zone links for a given zone
    max_nlink=0
    do i=1,nrz
       read(10,*) izn
       if(izn>nrz) then
          write(51,*) 'REGULARIZATION HAS BEEN SPECIFIED FOR ZONE ',izn
          write(51,*) 'BUT ONLY ',nrz,' ZONES ARE SPECIFIED IN ',trim(invfile)
          stop
       end if
       read(10,*) reg_opt(izn)       
       read(10,*) itmp(izn)
       if(itmp(izn)>max_nlink) then
          max_nlink=itmp(izn)
       end if
       read(10,*) zwts(izn,1)
    end do

    !!check for geostat constraints
    read(10,*) i
    if(i==1) then
       gs_flag = .true.
    end if

    rewind(10)
    allocate(zone_links(nrz,max_nlink+1))

    close(51)

    open(51,file='e4d.log',status='old',action='write',position='append')
    zwts = 0
    read(10,*) nrz
    write(51,*) '            ------------------------------------------'
    do i=1,nrz

       read(10,*) izn
       if(reg_opt(izn)==2 .or. reg_opt(izn)==4) then
          read(10,*) reg_opt(izn),zwts(izn,2:4)
       else
          read(10,*) reg_opt(izn)
       end if
       
       read(10,*) zone_links(izn,1:itmp(izn)+1)
       read(10,*) zwts(izn,1)
       write(51,*) '            zone: ',izn
       write(51,*) '            regularization option: ',reg_opt(izn)
       if(reg_opt(izn)==2) then
          write(51,*) '            x y z weights: ',zwts(izn,2:4)
       end if
       write(51,*) '            linked to zones: ',zone_links(izn,2:zone_links(izn,1)+1)
       write(51,*) '            zone weight: ',zwts(izn,1)
       write(51,*) '            ------------------------------------------'
       
       if(zone_links(izn,1)>nrz-1) then 
          write(51,*) "TOO MANY ZONE LINKS SPECIFIED"
          write(51,*) 'THERE ARE ONLY ',nrz-1, 'ZONES TO LINK WITH (see "linked to zones" above)'
          stop
       end if
       do j=2,zone_links(izn,1)+1
          if(zone_links(izn,j)==izn) then
             write(*,*) 'A ZONE MAY NOT BE LINKED WITH ITSELF (see "linked to zones" above)'
             stop
          end if
       end do
       
    end do
    close(51)

    
    if(.not.gs_flag) then
       read(10,*) i
    else
       read(10,*) i,gs_file
       !call get_geo_opts(nrz,reg_opt,gs_file)
    end if
    
    !!check for reference model
    do i=1,nrz
       if(reg_opt(i)==3) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "AT LEAST ONE ZONE USES TRANSIENT REGULARIZATION"
          write(51,*) "READING REFERENCE MODEL ",refmod_file
      
          open(11,file=trim(refmod_file),status='old',action='read')
          read(11,*) nstmp
          if(nstmp .ne. nsig) then
             write(51,*) "!!!!!WARNING!!!!!"
             write(51,*) "THE NUMBER OF ELEMENTS IN THE STARTING MODEL "
             write(51,*) "AND THE REFERENCE MODEL ARE DIFFERENT ",nsig,nstmp
             
          end if
          allocate(refsig(nsig))
          do j=1,nstmp
             read(11,*) refsig(j)
          end do
          close(11)
          close(51)
          goto 100 
       else
          allocate(refsig(nsig))  
          refsig=1.0
          goto 100
       end if
    end do

100 continue
    
    open(51,file='e4d.log',status='old',action='write',position='append')
    read(10,*) beta, beta_red
    beta_s = beta
    if(beta < 0) then
       write(51,*) "The beta value in the inverse options file must be greater than 0."
       write(51,*) "The value you provided is, ",beta
       stop
    end if
    if(beta_red < 0 .or. beta_red > 1) then
       write(51,*) "The beta reduction parameter in the inverse options file should"
       write(51,*) "be in the range 0 < beta_red < 1. The value you provided is ",beta_red
       beta_red = 0.5
       write(51,*) "beta_red has been changed to the default value of 0.5"
    end if
    write(51,*) "BETA AND BETA_RED ARE ",beta,beta_red
    close(51)

    read(10,*) conv_opt, norm_chi2, del_obj 
    
    if(conv_opt == 0) conv_opt = 1
    if(norm_chi2 == 0) norm_chi2 = 1.0
    if(del_obj == 0) del_obj = 0.05
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "CONVERGENCE OPTION = ",conv_opt
    write(51,*) "TARGET CHI-SQUARED = ",norm_chi2
    write(51,*) "MINIMUM OBJECTIVE FUNCTION REDUCTION = ",del_obj
    close(51)

    read(10,*) min_initer,max_initer,initer_conv,delta_initer
    if(min_initer == 0) min_initer = 20 
    if(max_initer == 0) max_initer = 100
    if(delta_initer == 0) delta_initer = 0.05
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "MINIMUM INNER ITERATIONS = ",min_initer
    write(51,*) "MAXIMUM INNER ITERATIONS = ",max_initer
    close(51)

    read(10,*) max_sig,min_sig
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "MIN/MAX CONDUCTIIVITY = ", min_sig,max_sig
 

    read(10,*) up_opt  
    if(up_opt == 1) then
       read(10,*) red_fac
       read(10,*) nlsp
       if(nlsp < 1) then
          write(51,*) "Line search updating has been specified in the inverse options file."
          write(51,*) "However, you have specified ",nlsp," line search points."
          stop
       end if
       nlsp=nlsp+1
       allocate(lsp(nlsp))
       read(10,*) lsp(2:nlsp)
       lsp(1) = 0
       write(51,*) "LINE SEARCH SPECIFIED"
       
    else if(up_opt == 2) then
       read(10,*) red_fac
       !read(10,*) js_perc !not used
       !read(10,*) js_perc !not used
       
       write(51,*) "NO LINE SEARCH SPECIFIED WITH UPDATE VECTOR SCALING OF ",red_fac
       
    else
       write(51,*) "The update option parameter up_opt in the inversion options"
       write(51,*) "file must be either 1 or 2. The value you provided is ", up_opt
       stop
       
    end if
    close(51)

    !!read in the cull options
    read(10,*) cull_flag, cull_dev
    allocate(Wd_cull(nm))
    Wd_cull = 1.0
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "CULL DATA OPTION = ",cull_flag," CULL STD. DEVIATION = ",cull_dev
    write(51,*) "____________________________________________________________________"
    write(51,*) 
    close(51)

    !!check for external constraints file and read if it exists
    nec=0
    if(gs_flag) then
       read(10,*) fnam
       inquire(file=trim(fnam),exist=exst)
       if(exst) then
         
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "READING EXTERNAL CONSTRAINTS FILE (GEOSTAT): ",trim(fnam)
          write(51,*) "____________________________________________________________________"
          write(51,*) 
          close(51)
          open(11,file=fnam,status='old',action='read')
          read(11,*) nec
          allocate(ex_cols(nec),ex_vals(nec,3),gwts(nec))
          do i=1,nec
             read(11,*) ex_cols(i),ex_vals(i,1:3),gwts(i)
          end do
          close(11)
       else  
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "CANNOT FIND EXTERNAL CONSTRAINTS FILE: ",trim(fnam)
          write(51,*) "EXTERNAL CONSTRAINTS WILL NOT BE IMPLEMENTED"
          write(51,*) "____________________________________________________________________"
          write(51,*) 
          close(51)
       end if

    else
       read(10,*) fnam
       inquire(file=trim(fnam),exist=exst)
       if(exst) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "READING EXTERNAL CONSTRAINTS FILE (REGULARIZED): ",trim(fnam)
          write(51,*) "____________________________________________________________________"
          write(51,*) 
          close(51)
          open(11,file=fnam,status='old',action='read')
          read(11,*) nec
          allocate(ex_cols(nec),ex_vals(nec,3),gwts(nec))
          do i=1,nec
             read(11,*) ex_cols(i),ex_vals(i,1:3)
          end do
          close(11)
          write(*,*) my_rank,nec
       else  
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "CANNOT FIND EXTERNAL CONSTRAINTS FILE: ",trim(fnam)
          write(51,*) "EXTERNAL CONSTRAINTS WILL NOT BE IMPLEMENTED"
          write(51,*) "____________________________________________________________________"
          write(51,*) 
          close(51)
       end if
       
    end if
    close(10)
    
1001 format(' CANNOT FIND THE INVERSION OPTIONS FILE ',A)
  end subroutine get_inv_opts
  !_________________________________________________________________________

  !_________________________________________________________________________
  subroutine get_dobs_rt(ind)
    implicit none
    integer :: ind
    character*80 :: cfil
    integer :: ios,npre
    logical :: exst
 
    npre = len(trim(tl_file))
    inquire(file=trim(efile), exist=exst)
    if(exst) then
         call system('mv '//trim(efile)//' '//trim(efile)//'.proc') 
    end if

    write(*,*) 'Checking for new survey files with the prefix: ',trim(tl_file)
    do while(.true.) 
       call system('/bin/ls -1 *.srv > slist.txt');
       open(13,file='slist.txt',status='old',action='read')
       read(13,*,IOSTAT=ios) tl_dfils(1)      
       cfil = tl_dfils(1)

       if(ios .eq. 0) then
          inquire(file=trim(cfil), exist=exst)
	  write(*,*) exst
	  write(*,*) cfil(1:npre)
	  write(*,*) tl_file(1:npre)
	  write(*,*) trim(efile),trim(cfil)
	  write(*,*) 
          !if(exst .and. (trim(efile) .ne. trim(cfil))) then
	  if(exst) then	
             if(cfil(1:npre) .eq. tl_file(1:npre)) then 
                call get_dobs_tl(1)
                call system('mv '//trim(cfil)//' '//trim(cfil)//'.proc')
                write(*,*) "FOUND A NEW SURVEY FILE: ",trim(cfil)
                open(51,file='e4d.log',status='old',action='write',position='append')
                write(51,*) "FOUND A NEW SURVEY FILE: ",trim(cfil)
                write(51,*) 
                close(51)
                close(13)
                return
             else
                write(*,*) "SKIPPING SURVEY FILE: ",trim(cfil)
                write(*,*) "BECAUSE "//cfil(1:npre)//" not equal to "//tl_file(1:npre)
                open(51,file='e4d.log',status='old',action='write',position='append')
                write(51,*) "SKIPPING SURVEY FILE: ",trim(cfil)
                write(51,*) "BECAUSE "//cfil(1:npre)//" not equal to "//tl_file(1:npre)
                write(51,*) 
                close(51)
                call system('mv '//trim(cfil)//' '//trim(cfil)//'.skip')         
             end if
          end if
       else
          call sleep(1)
       end if
       close(13)

    end do
  end subroutine get_dobs_rt
  !_________________________________________________________________________
  !_________________________________________________________________________
  subroutine get_dobs_tl(ind)
    implicit none
    integer :: ind,ne_test,sflg,jnk,i,nmt,ios,j
    real :: dobst,Wdt,dobsti,Wdti
    real , dimension(3) :: etp
    integer, dimension(4) :: abmn
    real :: perr=0.001

    open(11,file=trim(tl_dfils(ind)),status='old',action='read')      
    read(11,*,IOSTAT=ios) ne_test
    if(ios .ne. 0) goto 99
    if(ne_test .ne. ne) goto 100
    do i=1,ne
       read(11,*,IOSTAT=ios) jnk,etp,sflg
       if(ios .ne. 0) goto 99
       if(abs(etp(1)-e_pos(i,1)-xorig)>perr .or. abs(etp(2)-e_pos(i,2)-yorig)>perr .or. abs(etp(3)-e_pos(i,3)-zorig)>perr) goto 101
    end do

    read(11,*,IOSTAT=ios) nmt
 
  
    if(ios .ne. 0) goto 99
    if(nmt .ne. nm) goto 102
    do i=1,nm
       if (i_flag) then
          read(11,*,IOSTAT=ios) jnk,abmn,dobst,Wdt,dobsti,Wdti
       else
          read(11,*,IOSTAT=ios) jnk,abmn,dobst,Wdt
       end if
       if(ios .ne. 0) goto 99
       if(Wdt.le.0) Wdt=1e9
       if(Wdti.le.0) Wdti=1e9
       do j=1,4
          if(abmn(j) .ne. s_conf(i,j)) goto 103
          dobs(i)=dobst
          Wd(i)=1/Wdt   
          
          if (i_flag) then
			dobsi(i)=dobsti
			Wdi(i)=1/Wdti
          end if 
       end do
    end do
    close(11)

    return
99 continue
    write(*,*) "INPUT ERROR: See e4d.log"
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "There was an error reading: ",trim(tl_dfils(ind))
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit
    return

100 continue
    write(*,*) "INPUT ERROR: See e4d.log"  
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "The number of electrodes in: ",trim(tl_dfils(ind))," is: ",ne_test
    write(51,*) "The number of electrodes in the base survey is: ",ne
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit
    return

101 continue 
    write(*,*) "INPUT ERROR: See e4d.log"  
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "The electrode position in: ",trim(tl_dfils(ind))," for electrode: ",i
    write(51,*) "is different than the same baseline electrode position."
    write(51,*) trim(tl_dfils(ind))," gives: ",etp
    write(51,*) "The baseline file gives: ",e_pos(i,1)-xorig,e_pos(i,2)-yorig,e_pos(i,3)-zorig
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit
    return

102 continue
    write(*,*) "INPUT ERROR: See e4d.log"  
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "The number of measurements in: ",trim(tl_dfils(ind))," is: ",nmt
    write(51,*) "The number of measurements in the base survey is: ",nm
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit
    return

103 continue
    write(*,*) "INPUT ERROR: See e4d.log"  
    open(51,file='e4d.log',status='old',action='write',position='append')
    write(51,*) "ABMN for meas.: ",i," in: ",trim(tl_dfils(ind))," is: ",abmn
    write(51,*) "The baseline abmn is: ",abmn
    write(51,*) "Aborting ...."
    close(51)
    call crash_exit
    return

  end subroutine get_dobs_tl
  !_________________________________________________________________________


     

  !____________________________________________________________________________________________
  subroutine check_inp(spot,indx)
    implicit none
    integer :: spot,indx
    logical :: exst, mcheck

    select case(spot)

    case(0)
       inquire(file='e4d.inp',exist=exst)
       if(.not. exst) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "Cannot find the primary input file e4d.inp: aborting"
          close(51)
          call crash_exit
       end if

    case(101)
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) "There was a problem reading the mode in e4d.inp: aborting"
          close(51)
          call crash_exit
       end if

    case(1)
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) " There was a problem reading the mode in e4d.inp."
          write(51,*) " Aborting ..."
          write(*,*) " There was a problem reading the mode in e4d.inp."
          write(*,*) " Aborting ..."
          close(51)
          call crash_exit
       end if
       mcheck = .false.
       select case(mode)
       case(0) 
          mcheck = .true.
       case(1) 
          mcheck = .true.
       case(2) 
          mcheck = .true.
       case(3) 
          mcheck = .true.
       case(4) 
          mcheck = .true.
       case(5)
          mcheck = .true.
#ifdef sip
       case(21) 
          mcheck = .true.
       case(22) 
          mcheck = .true.
       case(23) 
          mcheck = .true.
       case(24) 
          mcheck = .true.
#endif
       case(31) 
          mcheck = .true.
       case(32) 
          mcheck = .true.
       case(33) 
          mcheck = .true.
       case(34)
          mcheck = .true.
       case(35)
          mcheck = .true.
#ifdef sip
       case(41) 
          mcheck = .true.
       case(42) 
          mcheck = .true.
       case(43) 
          mcheck = .true.
       case(44) 
          mcheck = .true.

#endif
#ifdef resmode
       case(101)
       mcheck = .true.
#endif   

       end select

       if(.not.mcheck) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,"(A33,I3,A19)") "The mode selected in e4d.inp is ",mode,", which is invalid."
          write(51,*)"Valid run modes include: ..."
          write(51,*)  "  FUNCION                        MODE "
          write(51,*)  " Input Check                       0"
          write(51,*)  " Mesh Generation (ERT)             1"
          write(51,*)  " ERT Forward                       2"
          write(51,*)  " ERT Inverse Static                3"
          write(51,*)  " ERT Inverse Time-lapse            4"
          write(51,*)  " Real Time Inverse Time-lapse      5"
#ifdef sip
          write(51,*)  " Mesh Generation (SIP)            21"
          write(51,*)  " SIP Forward                      22"
          write(51,*)  " SIP Inverse Static               23"
          write(51,*)  " SIP Inverse Multi-frequency      24"
#endif
          write(51,*)  " Tank ERT Forward                 32"
          write(51,*)  " Tank ERT Inverse Static          33"
          write(51,*)  " Tank ERT Inverse Time-lapse      34"
          write(51,*) 
#ifdef sip
          write(51,*)  " Tank SIP Mesh Generation         41"
          write(51,*)  " Tank SIP Forward                 42"
          write(51,*)  " Tank SIP Inverse Static          43"
          write(51,*)  " Tank SIP Inverse Multi-frequency 44"
          write(51,*)
#endif
#ifdef resmode
          write(51,*)  " Data Deficiency Survey Opt.     101"
#endif 
          write(51,*)  "Aborting ..."
          close(51)

          write(*,"(A33,I3,A19)") "The mode selected in e4d.inp is ",mode,", which is invalid."
          write(*,*)"Valid run modes include: ..."
          write(*,*)  "  FUNCION                        MODE "
          write(*,*)  " Input Check                       0"
          write(*,*)  " Mesh Generation (ERT)             1"
          write(*,*)  " ERT Forward                       2"
          write(*,*)  " ERT Inverse Static                3"
          write(*,*)  " ERT Inverse Time-lapse            4"
#ifdef sip          
          write(*,*)  " Mesh Generation (SIP)            21"
          write(*,*)  " SIP Forward                      22"
          write(*,*)  " SIP Inverse Static               23"
          write(*,*)  " SIP Inverse Multi-frequency      24"
#endif
          write(*,*)  " Tank ERT Forward                 32"
          write(*,*)  " Tank ERT Inverse Static          33"
          write(*,*)  " Tank ERT Inverse Time-lapse      34"
          write(*,*) 
#ifdef sip
          write(*,*)  " Tank SIP Mesh Generation         41"
          write(*,*)  " Tank SIP Forward                 42"
          write(*,*)  " Tank SIP Inverse Static          43"
          write(*,*)  " Tank SIP Inverse Multi-frequency 44"
          write(*,*)
#endif
#ifdef resmode
          write(51,*) " Data Deficiency Survey Opt.     101"
#endif 
          write(*,*)  " Aborting ..."
          
          call crash_exit
         
       else
          
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) 
          if(analytic) then
             write(51,*) "*************** RUNNING IN ANALYTIC FORWARD MODE ***************"
          else
             select case(mode)
             case(0) 
                write(51,*) "****************** RUNNING IN INPUT TEST MODE ******************"
             case(1) 
                write(51,*) "************* RUNNING IN MESH GENERATION MODE (ERT) ************"
             case(2) 
                if(analytic) then
                   write(51,*) "************* RUNNING IN ERT FORWARD ANALYTIC MODE *************"
                else
                   write(51,*) "***************** RUNNING IN ERT FORWARD MODE ******************"
                end if
             case(3) 
                write(51,*) "**************** RUNNING IN ERT INVERSION MODE *****************"
             case(4) 
                write(51,*) "*********** RUNNING IN TIME-LAPSE ERT INVERSION MODE ***********"
             case(5)
                write(51,*) "****** RUNNING IN REAL-TIME TIME-LAPSE ERT INVERSION MODE ******"
             case(21) 
                write(51,*) "************* RUNNING IN MESH GENERATION MODE (SIP) ************"
             case(22) 
                write(51,*) "***************** RUNNING IN SIP FORWARD MODE ******************"
             case(23) 
                write(51,*) "**************** RUNNING IN SIP INVERSION MODE *****************"
             case(24) 
                write(51,*) "**************** RUNNING IN SIP MULTI-FREQUENCY INVERSION MODE *****************"
             case(31) 
                write(51,*) "************* RUNNING IN MESH GENERATION MODE (TANK) ***********"
             case(32) 
                write(51,*) "*************** RUNNING IN ERT FORWARD MODE (TANK) *************"
             case(33)       
                write(51,*) "************** RUNNING IN ERT INVERSION MODE (TANK) ************"
             case(34)
                write(51,*) "******** RUNNING IN TIME-LAPSE ERT INVERSION MODE (TANK) *******"
             case(35)
                write(51,*) "*** RUNNING IN REAL-TIME TIME-LAPSE ERT INVERSION MODE (TANK) **"
             case(41)
                write(51,*) "********** RUNNING IN SIP MESH GENERATION MODE (TANK) **********"
             case(42)
                write(51,*) "************** RUNNING IN SIP FORWARD MODE (TANK) **************"
             case(43)
                write(51,*) "************* RUNNING IN SIP INVERSION MODE (TANK) *************"
             case(44) 
                write(51,*) "**************** RUNNING IN SIP MULTI-FREQUENCY INVERSION MODE (TANK) *****************"
             case(101)
                write(51,*) "*******RUNNING IN DATA DEFICIENCY SURVEY OPTIMIZATION MODE******"
              end select
             write(51,"(A,I3.3)") "  Mode:                             ",mode
          end if
          close(51)
       end if

    case(2)
       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then

          if(mode == 1 .or. mode == 21 .or. mode == 31) then
             write(51,*) "  ERROR: There was a problem reading the mesh config. file name in e4d.inp"
             write(51,*) "  Aborting ..."
             close(51)
             write(*,*) "  ERROR: There was a problem reading the mesh config. file name in e4d.inp"
             write(*,*) "  Aborting ..."
          else
             open(51,file='e4d.log',status='old',action='write',position='append')
             write(51,*) "  ERROR: There was a problem reading the mesh file name in e4d.inp"
             write(51,*) "  Aborting ..."
             close(51)
             write(*,*) "  ERROR: There was a problem reading the mesh file name in e4d.inp"
             write(*,*) "  Aborting ..."
          end if
          call crash_exit

       else
          if(mode == 1 .or. mode == 21 .or. mode == 31) then
             write(51,*) " Mesh configuration file:          ",trim(mshfile)
          else
             write(51,*) " Mesh file:                        ",trim(mshfile)
          end if
       end if
       close(51)

    case(3)

       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the survey file name in e4d.inp: aborting"
          write(*,*) "There was a problem reading the survey file name in e4d.inp: aborting"
          close(51)
          call crash_exit
       else
          write(51,*) " Survey configuration file:        ",trim(efile)
       end if
       close(51)

    case(4)
       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the conductivity file name in e4d.inp: aborting"
          write(*,*) "There was a problem reading the conductivity file name in e4d.inp: aborting"
          close(51)
          call crash_exit
       else
          if(mode == 2 .or. mode == 22 .or. mode == 32) then
             write(51,*) " Conductivity file:                ",trim(sigfile)
          else
             write(51,*) " Starting conductivity file:       ",trim(sigfile)
          end if
       end if
       close(51)

    case(5)
      
       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the output options file name in e4d.inp: aborting"
          write(*,*) "There was a problem reading the output options file name in e4d.inp: aborting"
          close(51)
          call crash_exit
       else
          write(51,*) " Output options file:              ",trim(outfile)
       end if
       close(51)

    case(6)
    	if(im_fmm) then
    	   open(51,file='fmm.log',status='old',action='write',position='append')
    	else	
       	   open(51,file='e4d.log',status='old',action='write',position='append')
       	end if
       	
       if(ios .ne. 0) then
         if(im_fmm) then
            write(51,*) "There was a problem reading the inverse options file name in fmm.inp: aborting"
            write(*,*) "There was a problem reading the inverse options file name in fmm.inp: aborting"
            close(51)
            call crash_exit
         else
            if(mode == 22 .or. mode == 23) then
                write(51,*) "There was a problem reading the two SIP inverse option file names in e4d.inp: aborting"
                write(*,*) "There was a problem reading the two SIP inverse option file names in e4d.inp: aborting"
             else if (mode==0) then
                if (i_flag) then
					write(51,*) "No ERT or SIP inverse options file name specified in e4d.inp."
					write(*,*) "No ERT or SIP inverse options file name specified in e4d.inp."
				else
					write(51,*) "No ERT inverse options file name specified in e4d.inp."
					write(*,*) "No ERT inverse options file name specified in e4d.inp."				
				end if
                inv_opt_flag=.false.
                return
             else
                write(51,*) "There was a problem reading the inverse options file name in e4d.inp: aborting"
                write(*,*) "There was a problem reading the inverse options file name in e4d.inp: aborting"
             end if
             close(51)
             call crash_exit
         end if
       
       else
         
         if(mode == 22 .or. mode == 23) then
             write(51,*) " Inverse options file (amp.):      ",trim(invfile)
             write(51,*) " Inverse options file (phase):     ",trim(iinvfile)
          else
             write(51,*) " Inverse options file:             ",trim(invfile)
          end if
          
       end if
       close(51)

    case(7)
        if(im_fmm) then
    	   open(51,file='fmm.log',status='old',action='write',position='append')
    	else	
       	   open(51,file='e4d.log',status='old',action='write',position='append')
       	end if
       
       if(ios .ne. 0) then
          if (mode==0) then
			  write(51,*) "No reference model file name specified in e4d.inp."
			  write(*,*) "No reference model file name specified in e4d.inp."
		  else
			  write(51,*) "There was a problem reading the reference model file name in e4d.inp: aborting"
			  write(*,*) "There was a problem reading the reference model file name in e4d.inp: aborting"
			  call crash_exit
		  end if
          close(51)
          
       else
          write(51,*) " Reference model file:             ",trim(refmod_file)
       end if
       close(51)


    case(8)

       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          if (mode==0) then
			  write(51,*) "There was no time-lapse survey file name specified and/or "
			  write(51,*) "reference model update option in e4d.inp."
			  close(51)
			  write(*,*) "There was no time-lapse survey file name specified and/or "
			  write(*,*) "reference model update option in e4d.inp."
			  return
          else 
			  write(51,*) "There was a problem reading the time-lapse survey file name"
			  write(51,*) "and/or reference model update option in e4d.inp: aborting"
			  close(51)
			  write(*,*) "There was a problem reading the time-lapse survey file name"
			  write(*,*) "and/or reference model update option in e4d.inp: aborting"
			  call crash_exit  
		  end if
          
       else
     	  tl_ly = .true.

          if (i_flag) then
            write(51,*) " Multiple frequency survey list file:      ",trim(tl_file)
          else
			write(51,*) " Time-lapse survey list file:      ",trim(tl_file)
		  end if
          if(indx==2) then
             write(51,*) " Reference model update option:    previous solution"
          else
             write(51,*) " Reference model update option:    baseline solution"
          end if
       end if
       
       inquire(file=trim(tl_file),exist=exst)
       if(.not.exst .and. .not. rt_flag) then
          write(*,*) 
          write(51,*) " Cannot find the time lapse survey list file: ",trim(tl_file)," ... aborting"
          write(*,*) " Cannot find the time lapse survey list file: ",trim(tl_file)," ... aborting"
          close(51)
          call crash_exit
       end if
       close(51)


    case(9)

       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading the number of survey files in ",trim(tl_file),"  :aborting"
          write(*,*) "There was a problem reading the number of survey files in ",trim(tl_file),"  :aborting"
          close(51)
          call crash_exit
       else
          write(51,*)
          if (i_flag) then
            write(51,*) " Frequency survey list file:      ",trim(tl_file)
            write(51,*) " Num. of frequency survey files:  ",ntl
          else
            write(51,*) " Time-lapse survey list file:      ",trim(tl_file)
			write(51,*) " Num. of time-lapse survey files:  ",ntl
		  end if
       end if 
       close(51)
       

    case(10)
       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          write(51,*) "There was a problem reading time-lapse file at time ",indx
          write(51,*) "in the time lapse survey file: ",trim(tl_file)
          write(*,*) "There was a problem reading time-lapse file at time ",indx
          write(*,*) "in the time lapse survey file: ",trim(tl_file)
          close(51)
          call crash_exit
       else
         
         if(indx==1) then
            if (i_flag) then
				write(51,*) " Index     Frequency Survey File       Data Stamp"
            else
				write(51,*) " Index     Time-Lapse Survey File       Time Stamp"
			end if
         end if
         write(51,"(I6,A28,g17.5)"),indx,trim(tl_dfils(indx)),tlt(indx)
       end if
      close(51)


    case(11)
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) " There was a problem finding the time-lapse file: ",trim(tl_dfils(indx))
       write(51,*) " aborting ..."
       write(*,*) " There was a problem finding the time-lapse file: ",trim(tl_dfils(indx))
       write(*,*) " aborting ..."
       close(51)
       call crash_exit
      
    case(12)

       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*) 
          write(51,*) " There was a problem reading the number of electrodes in the survey file ",trim(efile)
          write(51,*) " Aborting ..."
          write(*,*) 
          write(*,*) " There was a problem reading the number of electrodes in the survey file ",trim(efile)
          write(*,*) " Aborting ..."
          
          close(51)
          call crash_exit
       else
          write(51,"(A,I7.7)") "  Number of electrodes:             ",ne
       end if
       close(51)

    case(13)
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*)
          write(51,"(A,I7.7)") "  There was a problem reading electrode number ",indx
          write(51,*) " in the survey file file: ",trim(efile)
          write(51,*) " Aborting ..."
          close(51)
          write(*,*)
          write(*,"(A,I7.7)") " There was a problem reading electrode number ",indx
          write(*,*) "in the survey file file: ",trim(efile)
          write(*,*) "Aborting ..."
          call crash_exit
       end if

    case(14)
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) " The electrode index specified for electrode: ",indx
       write(51,*) " is greater than the total number of electrodes."
       write(51,*) " Aborting ..."
       write(*,*) 
       write(*,*) " The electrode index specified for electrode: ",indx
       write(*,*) " is greater than the total number of electrodes."
       write(*,*) " Aborting ..."
     
       close(51)
       
       call crash_exit

    case(15)
       inquire(file=mshfile(1:mnchar)//'trn',exist=exst)
       if(.not. exst) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn' 
          write(51,*) " Aborting..."
          write(*,*)
          write(*,*) " Cannot find the mesh translation file: ",trim(mshfile(1:mnchar))//'trn' 
          write(*,*) " Aborting..."
          close(51)
          call crash_exit
       end if

    case(16)
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
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
          call crash_exit
       end if

    case(17)
       open(51,file='e4d.log',status='old',action='write',position='append')
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          write(51,*)
          write(51,*) " There was a problem reading the number of measurements" 
          write(51,*) " in the survey file: ",trim(efile)
          close(51)
          write(*,*)
          write(*,*) " There was a problem reading the number of measurements" 
          write(*,*) " in the survey file: ",trim(efile)
          close(51)
          call crash_exit
       elseif(nm .le. 0 .and. .not. res_flag) then
          write(51,*) " The number of measurements is not positive: ",indx
          write(51,*) " Aborting ..."
          write(*,*)
          write(*,*) " The number of measurements is not positive: ",indx
          write(*,*) " Aborting ..."
          close(51)
          call crash_exit
       else
           write(51,"(A,I7.7)") "  Number of measurements:           ",nm
       end if
       close(51)
       
    case(18)
       if(ios .ne. 0) then
          open(51,file='e4d.log',status='old',action='write',position='append')
          if(i_flag) then
             write(51,*)
             write(51,"(A,I8.8)") "  There was a problem reading measurement number ",indx
             write(51,*) " in the survey file: ",trim(efile)
             write(51,*) " (hint: make sure both the real and complex transfer"
             write(51,*) " resistances and standard deviations are listed) "
             write(51,*) " aborting ..."
             write(*,*)
             write(*,"(A,I8.8)") "  There was a problem reading measurement number ",indx
             write(*,*) " in the survey file: ",trim(efile)
             write(*,*) " (hint: make sure both the real and complex transfer"
             write(*,*) " resistances and standard deviations are listed) "
             write(*,*) " aborting ..."
          else
             write(51,*)
             write(51,"(A,I8.8)")" There was a problem reading measurement number ",indx
             write(51,*) " in the survey file: ",trim(efile)
             write(51,*) " aborting ..."
             write(*,*)
             write(*,"(A,I8.8)") "  There was a problem reading measurement number ",indx
             write(*,*) " in the survey file: ",trim(efile)
             write(*,*) " aborting ..."
          end if
          close(51)
          call crash_exit
       end if

   case(19)
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*)
      write(51,"(A,I8.8)") "  !!! WARNING: MEASUREMENT ",indx
      write(51,*) " !!! AND POSSIBLY OTHERS SPECIFY A NEGATIVE OR ZERO STANDARD DEVIATION"
      write(51,*) " !!! IN THE SURVEY FILE ",trim(efile)
      write(51,*) " !!! SETTING TO LARGE STANDARD DEVIATION"
      close(51)

      write(*,*)
      write(*,"(A,I8.8)") " !!! WARNING: MEASUREMENT ",indx
      write(*,*) "!!! AND POSSIBLY OTHERS SPECIFY A NEGATIVE OR ZERO STANDARD DEVIATION"
      write(*,*) "!!! IN THE SURVEY FILE ",trim(efile)
      write(*,*) "!!! SETTING TO LARGE STANDARD DEVIATION" 
      close(51)
      return

   case(20)
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*)
      write(51,"(A,I8.8,A)") "  Measurement ",indx," uses electrodes that are out of range." 
      write(51,"(A,I8.8,A)") "  There are ",ne," electrodes."
      write(51,"(A,I8.8,A,4I10.8)") "  A B M N for measurement ",indx," reads ",s_conf(indx,:)
      write(51,*) " Aborting ..."
      write(51,*)
      close(51)

      write(*,*)
      write(*,"(A,I8.8,A)") "  Measurement ",indx," uses electrodes that are out of range." 
      write(*,"(A,I8.8,A)") "  There are ",ne," electrodes."
      write(*,"(A,I8.8,A,4I10.8)") "  A B M N for measurement ",indx," reads ",s_conf(indx,:)
      write(*,*) " Aborting ..."
      write(*,*)
      call crash_exit
      
   case(21)
      open(51,file='e4d.log',status='old',action='write',position='append')
      if(indx==0) then
         write(51,*)
         write(51,*) " In mode 1 you must provide a mesh configuration (*.cfg) file "
         write(51,*) " In all other modes you must provide a node or element file name"
         write(51,*) " You provided: ",trim(mshfile)
         write(51,*) " Aborting ..."
         write(*,*)
         write(*,*) " In mode 1 you must provide a mesh configuration (*.cfg) file "
         write(*,*) " In all other modes you must provide a node or element file name"
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
      call crash_exit
      
   case(22)
      if(ios .ne. 0) then
         open(51,file='e4d.log',status='old',action='write',position='append')
         write(51,*) "There was a problem reading the number of conductivities in"
         write(51,*) "in the conductivity file: ",trim(sigfile)
         write(51,*) "aborting."
         close(51)
         call crash_exit
      end if


   case(23)
      if(indx == 0) then
         open(51,file='e4d.log',status='old',action='write',position='append')
         if(ios == 0) then
            write(51,*)
            if(i_flag) then
               write(51,*) " COMPLEX CONDUCTIVITY FILE SUMMARY "
               write(51,"(A,I10.10)") "  Number of complex cond. values:   ",nsig
            else
               write(51,*) " CONDUCTIVITY FILE SUMMARY "
               write(51,"(A,I10.10)") "  Number of conductivity values:    ",nsig
            end if
            close(51)
         else
            write(51,*) 
            write(51,*) " There was a problem reading the number of "
            write(51,*) " conductivity values in ",trim(sigfile)
            write(51,*) " Aborting ..."
            close(51)
            write(*,*) 
            write(*,*) " There was a problem reading the number of "
            write(*,*) " conductivity values in ",trim(sigfile)
            write(*,*) " Aborting ..."
            close(51)
            call crash_exit
         end if
      end if
      if(ios .ne. 0) then
         open(51,file='e4d.log',status='old',action='write',position='append')
         write(51,*)
         write(51,"(A,I10.10)") "  There was a problem reading conductivity number ",indx
         write(51,*) " in the conductivity file: ",trim(sigfile),"."
         if(i_flag) write(51,*) " (hint: make sure both amplitude and phase are listed for each element"
         write(51,*) " Aborting ..."
         close(51)
         write(*,*)
         write(*,"(A,I10.10)") "  There was a problem reading conductivity number ",indx
         write(*,*) " in the conductivity file: ",trim(sigfile),"."
         if(i_flag) write(*,*) " (hint: make sure both amplitude and phase are listed for each element"
         write(*,*) " Aborting ..."
         close(51)
         call crash_exit
      end if

  case(24)
     open(51,file='e4d.log',status='old',action='write',position='append')
     if(indx == 0) then
        write(51,*)
        write(51,*) " Can't find the survey file: ",trim(efile)
        write(51,*) " aborting."
        write(*,*)
        write(*,*) " Can't find the survey file: ",trim(efile)
        write(*,*) " aborting."
        close(51)
        call crash_exit
     else
        write(51,*) 
        write(51,*) " SURVEY FILE SUMMARY"
     end if
    close(51)
      
  case(25)
     open(51,file='e4d.log',status='old',action='write',position='append')
     write(51,*) 
     write(51,*) " Can't find the conductivity file: ",trim(sigfile)
     write(51,*) " Aborting ..."
     close(51)
     write(*,*) 
     write(*,*) " Can't find the conductivity file: ",trim(sigfile)
     write(*,*) " Aborting ..."
     close(51)
     call crash_exit
      
  case(26)
     open(51,file='e4d.log',status='old',action='write',position='append')
     if(i_flag) then
        write(51,*)
        write(51,*) " All real and complex conductivities must be positive"
        write(51,*) " The conductivities on line: ",indx," are ",sigma(indx),sigmai(indx)
        write(51,*) " See conductivity file: ",trim(sigfile)
        write(51,*) " Aborting ..."
        
        write(*,*) 
        write(*,*) " All real and complex conductivities must be positive"
        write(*,*) " The conductivities on line: ",indx," are ",sigma(indx),sigmai(indx)
        write(*,*) " See conductivity file: ",trim(sigfile)
        write(*,*) " Aborting ..."
       
     else
        write(51,*)
        write(51,*) " All conductivities must be positive"
        write(51,*) " The conductivity on line: ",indx," is ",sigma(indx)
        write(51,*) " See conductivity file: ",trim(sigfile)
        write(51,*) " Aborting ..."
        write(51,*) 
        write(*,*)
        write(*,*) " All conductivities must be positive"
        write(*,*) " The conductivity on line: ",indx," is ",sigma(indx)
        write(*,*) " See conductivity file: ",trim(sigfile)
        write(*,*) " Aborting ..."
        write(*,*) 
        
     end if

     close(51)
     call crash_exit

     case(27)
        open(51,file='e4d.log',status='old',action='write',position='append')
     
        write(51,*)
        write(51,*) " In tank mode the last electrode specified in the survey file must"
        write(51,*) " be the ghost node (see E4D User Guide)"
        write(51,*) " The ghost node cannot be used as an electrode in a survey."
        write(51,*) " Measurement: ",indx," uses the ghost node as an electrode."
        write(51,*) " Aborting ..."
        write(51,*)
       
           
        write(*,*)
        write(*,*) " In tank mode the last electrode specified in the survey file must"
        write(*,*) " be the ghost node (see E4D User Guide)"
        write(*,*) " The ghost node cannot be used as an electrode in a survey."
        write(*,*) " Measurement: ",indx," uses the ghost node as an electrode."
        write(*,*) " Aborting ..."
        write(*,*)

        close(51)
        call crash_exit

     case(28)
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) 
        write(51,*) "Can't find the conductivity list file ",trim(sigfile)," specified in e4d.inp."
        write(51,*) "aborting ..."
        write(*,*) "Can't find the conductivity list file ",trim(sigfile)," specified in e4d.inp."
        write(*,*) "aborting ..."
        close(51)
        call crash_exit

     case(29)
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*) 
        write(51,*) "There was a problem reading the number of files listed "
        write(51,*) "on the first lint of ",trim(sigfile),"."
        write(51,*) "aborting ..."
        write(*,*) "There was a problem reading the number of files listed "
        write(*,*) "on the first lint of ",trim(sigfile),"."
        write(*,*) "aborting ..."

        close(51)
        call crash_exit

     case(30)
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*)
        write(51,*) "There was a problem reading conductivity file number ",indx
        write(51,*) "in the conductivity list file: ",trim(sigfile)
        write(51,*) "Remember to include the output file name in the second column"
        write(51,*) "aborting ..."
       
        write(*,*) "There was a problem reading conductivity file number ",indx
        write(*,*) "in the conductivity list file: ",trim(sigfile)
        write(51,*) "Remember to include the output file name in the second column"
        write(*,*) "aborting ..."
        
        close(51)
        call crash_exit
        
     case(31)
        open(51,file='e4d.log',status='old',action='write',position='append')
        write(51,*)
        write(51,*) "Cannot find find file ",trim(tl_cfils(indx,1))," listed in"
        write(51,*) "the conductivity list file: ",trim(sigfile)
        write(51,*) "aborting ..."
        write(*,*)
        write(*,*) "Cannot find find file ",trim(tl_cfils(indx,1))," listed in"
        write(*,*) "the conductivity list file: ",trim(sigfile)
        write(*,*) "aborting ..."
        close(51)
        call crash_exit

     case(121)
        open(51,file='e4d.log',status='old',action='write',position='append')
        if(mshfile(mnchar+2:mnchar+6) == ".node") then
           inquire(file=trim(mshfile),exist=exst)
           if(.not.exst) then
              write(51,*)
              write(*,*)
              write(51,*) " Cannot find the specified mesh node file: ",trim(mshfile)
              write(*,*) " Cannot find the specified mesh node file: ",trim(mshfile)
              close(51)
              call crash_exit
           end if
        elseif(mshfile(mnchar+2:mnchar+5) == ".ele") then
           inquire(file=trim(mshfile),exist=exst)
           if(.not.exst) then
              write(51,*)
              write(*,*)
              write(51,*) " Cannot find the specified mesh element file: ",trim(mshfile)
              write(*,*) " Cannot find the specified mesh element file: ",trim(mshfile)
              close(51)
              call crash_exit
           end if
        else
           write(51,*)
           write(*,*)
           write(51,*) " If mode > 1 you must provide the name of the mesh"
           write(51,*) " node file (*.node) or mesh element file (*.ele) ."
           write(51,*) " You provided: ",trim(mshfile)
           write(*,*) " If mode > 1 you must provide the name of the mesh"
           write(*,*) " node file (*.node) or mesh element file (*.ele) ."
           write(*,*) " You provided: ",trim(mshfile)
           close(51)
           call crash_exit
        end if
        close(51)
      
      
    case DEFAULT

    end select
  end subroutine check_inp
  !_________________________________________________________________________

  !_________________________________________________________________________
  subroutine check_inv_opts(spot,indx)
    implicit none
    integer :: spot,indx,ii,j,nzn,i
    logical :: exst, rpt
    character*40 :: str
    
    if(iter>0) return
   
    select case(spot)

    case(0)

       if(invi) then 
          inquire(file=trim(iinvfile),exist=exst)
       else
          inquire(file=trim(invfile),exist=exst)
       end if
       if(im_fmm) then
          open(51,file='fmm.log',status='old',action='write',position='append')
       else
          open(51,file='e4d.log',status='old',action='write',position='append')
       end if
       
       if(.not. exst) then

          write(51,*) 
          if(invi) then
             write(51,*) " Cannot find the phase inverse options file: ",trim(iinvfile)
             write(*,*) " Cannot find the phase inverse options file: ",trim(iinvfile)
          else
             write(51,*) " Cannot find the inverse options file: ",trim(invfile)
             write(*,*) " Cannot find the inverse options file: ",trim(invfile)
          end if
          write(51,*) " Aborting"
          write(*,*) " Aborting"
          close(51)
          call crash_exit

       else
 
          if(iter == 0) then
             write(51,*)
             if(invi) then
                write(51,*) " PHASE INVERSION OPTIONS "
                write(51,*) " Phase inversion options file:     ",trim(iinvfile)
             else
                if(i_flag) then
                   write(51,*) " AMPLITUDE INVERSION OPTIONS "
                   write(51,*) " Amplitude inversion options file: ",trim(invfile)
                else
                   write(51,*) " INVERSION OPTIONS "
                   write(51,*) " Inversion options file:           ",trim(invfile)
                end if
             end if
          end if
       end if
       close(51)

       case(1)
       	  if(im_fmm) then
       	     open(51,file='fmm.log',status='old',action='write',position='append')
       	  else
             open(51,file='e4d.log',status='old',action='write',position='append')
          end if
          
          if(ios .ne. 0) then
             write(51,*) 
             write(51,*) " There was a problem reading the number of constraint"
             if(invi) then
                write(51,*) " blocks in the inverse options file: ",trim(iinvfile)
             else
                write(51,*) " blocks in the inverse options file: ",trim(invfile)
             end if
             write(51,*) " Aborting"
             close(51)
             write(*,*) 
             write(*,*) " There was a problem reading the number of constraint"
             if(invi) then
                write(51,*) " blocks in the inverse options file: ",trim(iinvfile)
             else
                write(51,*) " blocks in the inverse options file: ",trim(invfile)
             end if
             write(*,*) " Aborting"
             close(51)
             call crash_exit

          elseif(nrz<0) then
             write(51,*) 
             write(51,*) " The number of constraint blocks is: ",nrz
             write(51,*) " The number of constraint blocks must be > 0"
             write(51,*) " Aborting"
             close(51)
             write(*,*) " The number of constraint blocks is: ",nrz
             write(*,*) " The number of constraint blocks must be > 0"
             write(*,*) " Aborting"
             close(51)
             call crash_exit

          else
             !write(51,*) " Inversion options file:           ",trim(invfile)
             !write(51,*) "   Number of constraint blocks:    ",nrz
   
       end if
          close(51)

       case(2)
      
          if(ios .ne. 0) then
             if(im_fmm) then
	       open(51,file='fmm.log',status='old',action='write',position='append')
	     else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
           
             write(51,*) 
             write(51,*) " There was a problem reading the zone number"
             write(51,*) " for constraint block ",indx,"."
             if(invi) then
                write(51,*) " Check the inverse options file: ",trim(iinvfile)
             else
                write(51,*) " Check the inverse options file: ",trim(invfile)
             end if
             write(51,*) " Aborting"
             close(51)

             write(*,*) 
             write(*,*) " There was a problem reading the zone number"
             write(*,*) " for constraint block ",indx,"."
             if(invi) then
                write(*,*) " Check the inverse options file: ",trim(iinvfile)
             else
                write(*,*) " Check the inverse options file: ",trim(invfile)
             end if
             write(*,*) " Aborting"
             close(51)
             call crash_exit

          end if

          case(3)
             if(ios .ne. 0) then
                if(im_fmm) then
	           open(51,file='fmm.log',status='old',action='write',position='append')
	        else
	           open(51,file='e4d.log',status='old',action='write',position='append')
                end if
                
                write(51,*) 
                write(51,*) " There was a problem reading either the structure metric"
                write(51,*) " or the spatial weights for constraint block ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)

                write(*,*) 
                write(*,*) " There was a problem reading either the structure metric"
                write(*,*) " or the spatial weights for constraint block ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting"

                call crash_exit
             end if

          case(4)
             if(ios .ne. 0) then
                if(im_fmm) then
		   open(51,file='fmm.log',status='old',action='write',position='append')
		else
		   open(51,file='e4d.log',status='old',action='write',position='append')
                end if
                
                write(51,*) 
                write(51,*) " There was a problem reading either the reweighting function"
                write(51,*) " or the reweighting function parameters for constraint block ",indx
                if(invi) then
                   write(51,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(51,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)

                write(*,*) 
                write(*,*) " There was a problem reading either the reweighting function"
                write(*,*) " or the reweighting function parameters for constraint block ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting"
                call crash_exit
             end if

          case(5)
             if(ios .ne. 0) then
                if(im_fmm) then
		   open(51,file='fmm.log',status='old',action='write',position='append')
		else
		   open(51,file='e4d.log',status='old',action='write',position='append')
                end if
                
                write(51,*) 
                write(51,*) " There was a problem reading the zone links"
                write(51,*) " for constraint block ",indx
                if(invi) then
                   write(51,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(51,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)

                write(*,*) 
                write(*,*) " There was a problem reading the zone links"
                write(*,*) " for constraint block ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting"

                call crash_exit
             end if

          case(6)
             if(ios .ne. 0) then
                if(im_fmm) then
		   open(51,file='fmm.log',status='old',action='write',position='append')
		else
		   open(51,file='e4d.log',status='old',action='write',position='append')
                end if
                write(51,*) 
                write(51,*) " There was a problem reading the reference"
                write(51,*) " for constraint block ",indx
                if(invi) then
                   write(51,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(51,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)

                write(*,*) 
                write(*,*) " There was a problem reading the reference"
                write(*,*) " for constraint block ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting"         
                call crash_exit
             end if

          case(7)
             if(ios .ne. 0) then
                if(im_fmm) then
		   open(51,file='fmm.log',status='old',action='write',position='append')
		else
		   open(51,file='e4d.log',status='old',action='write',position='append')
                end if
                write(51,*) 
                write(51,*) " There was a problem reading the relative weight"
                write(51,*) " for constraint block ",indx
                if(invi) then
                   write(51,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(51,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)

                write(*,*) 
                write(*,*) " There was a problem reading the relative weight"
                write(*,*) " for constraint block ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting"
                close(51)
                call crash_exit
             end if
          
          case(8)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             write(51,"(A,I3.3)") "  Number of constraint blocks:      ",nrz
             close(51)
             
          case(9)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             write(51,*) 
             write(51,"(A30,I10)") "CONSTRAINT BLOCK: ",indx
             close(51)
          
          case(10)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             if(indx==-1) then
                write(51,"(A47)")"Zone Number:          EXTERNAL"
             else
                write(51,"(A30,I10)") "Zone Number: ",smetric(indx,1)
             end if
             close(51)

          case(11)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if  
             if(ios .ne. 0 ) then
                write(51,*) 
                write(51,*) " There was a problem reading the directional weights"
                write(51,*) " for constraint block: ",indx
                if(invi) then
                   write(51,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(51,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)

                write(*,*) 
                write(*,*) " There was a problem reading the directional weights"
                write(*,*) " for constraint block: ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting"
                close(51)   
                call crash_exit

             else
                write(51,"(A30,I10)") "Structure Metric: ",smetric(indx,2)
                if(smetric(indx,2) < 0 .or. smetric(indx,2) > 13) then
                   write(51,*) " Invalid structural metric "
                   write(51,*) " Aborting"
                   close(51)
                   call crash_exit
                end if
                if(smetric(indx,2)==5 .or. smetric(indx,2)==6 .or. smetric(indx,2).ge.12) then
                   write(51,"(A30,3G10.3)") "X,Y,Z Weights: ",zwts(indx,2),zwts(indx,3),zwts(indx,4)
                end if
                close(51)
             end if
           
          case(12)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             if(ios .ne. 0) then
                write(51,*) 
                write(51,*) " There was a problem reading reweighting function parameters"
                write(51,*) " for constraint ",indx
                if(invi) then
                   write(51,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(51,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(51,*) " Aborting"
                close(51)
                
                write(*,*) 
                write(*,*) " There was a problem reading reweighting function parameters"
                write(*,*) " for constraint ",indx
                if(invi) then
                   write(*,*) " Check the inverse options file: ",trim(iinvfile)
                else
                   write(*,*) " Check the inverse options file: ",trim(invfile)
                end if
                write(*,*) " Aborting" 
                call crash_exit

             else
                write(51,"(A30,I10)") "Reweighting Function: ",Fw_type(indx)
                if(Fw_type(indx)<1 .or. Fw_type(indx)>6) then
                   write(51,*) " Invalid reweighting function"
                   write(51,*) " Aborting"
                   close(51)
                   call crash_exit
                end if
                if(Fw_parm(indx,2)<0.0) then
                   write(51,"(A30,8x,2g10.4)") "Reweighting Mean, St.Dev.: ",Fw_parm(indx,1),Fw_parm(indx,2)
                   write(51,*) " Reweighting function deviation must be greater than zero"
                   write(51,*) " Aborting ..."
                   close(51) 
                   call crash_exit
                end if
                write(51,"(A30,8x,2g10.4)") "Reweighting Mean, St.Dev.: ",Fw_parm(indx,1),Fw_parm(indx,2)
                close(51)
             end if

          case(13)
               
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             if(ios .ne. 0) then
                write(51,*) 
                write(51,*) " There was a problem reading the zone links for"
                write(51,*) " constraint block :",indx
                close(51)
                call crash_exit
             else
                write(51,"(A30,I10)") "Number of Zone Links: ",zone_links(indx,1)
                if(zone_links(indx,1)<0) then
                   write(51,*) " The number of zone links must be greater than or equal to zero"
                   write(51,*) " Aborting ..."
                   close(51) 
                   call crash_exit
                end if
                
                nzn=maxval(zones)
                do ii=1,zone_links(indx,1)
                   write(51,"(A30,I10)") "Linked to Zone: ",zone_links(indx,ii+1)
                   if(zone_links(indx,ii+1)>nzn) then
                      write(51,*) "  !!!WARNING: this zone link is out of range"
                   end if
                end do
             end if
             close(51)

          case(14)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             if(ios .ne. 0) then
                write(51,*) 
                write(51,*) " There was a problem reading the reference value"
                write(51,*) " for constraint block: ",indx
                write(51,*) " Aborting ..."
                close(51)
                call crash_exit
             else
                ii=smetric(indx,2)
                if(ii==3 .or. ii==4 .or. ii==7 .or. ii==8) then
                   if(smetric(indx,3)==1) then
                      write(51,"(A39,A)") "Reference model file:          ",trim(refmod_file)
                   elseif(smetric(indx,3)==2) then
                      write(51,"(A30,A26)") "Reference Value: ","         Previous Solution."
                   else
                      write(51,"(A30,g10.4)") "Reference Value: ",C_targ(indx)
                   end if
                end if
             end if
             close(51)

          case(15)
 
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             if(ios .ne. 0) then
                write(51,*) 
                write(51,*) " There was a problem reading the relative weight"
                write(51,*) " for constraint block ",indx
                write(51,*) " Aborting ..."
                close(51)
                call crash_exit
             else
                write(51,"(A30,8x,g10.4)") "Relative Weight: ",zwts(indx,1)
                if(zwts(indx,1)<0) then
                   write(51,*) " The relative weight for each block must be greater than or equal to zero"
                   write(51,*) " Aborting ..."
                   close(51)
                   call crash_exit
                end if
             end if
             close(51)

          case(16)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
             nzn=int(maxval(zones))
             do j=1,zone_links(indx,1)
                if(smetric(indx,1)==zone_links(indx,j+1)) then
                   write(51,"(A6,I5,A21,I5,A1)") " Zone ",smetric(indx,1)," is linked with zone ",zone_links(indx,j+1),"."
                   write(51,*)," A zone may not be linked with itself"
                   write(51,*) " Aborting"
                   close(51)
                   call crash_exit
                   return
                elseif(zone_links(indx,j+1).le.0 .and. zone_links(indx,j+1).gt.nzn) then
                   write(51,*) " A zone link is out of range for constraint block: ",indx
                   write(51,*) " Aborting"
                   close(51)
                   call crash_exit
                   return
                end if
             end do
             close(51)

          case(17)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
            write(51,*)
            write(51,"(A,I5,A)") "  !!!WARNING: ZONE ",indx," IS NOT CONSTRAINED"
            close(51)
             
         case(18)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
            write(51,*)
            inquire(file=trim(refmod_file),exist=exst)
            if(.not.exst) then
               write(51,*) " Cannot find the reference model file: ",trim(refmod_file)
               write(51,*) " Aborting ..."
               
               write(*,*) " Cannot find the reference model file: ",trim(refmod_file)
               close(51)
               call crash_exit
            else
               !write(51,*) "READING REFERENCE MODEL ",trim(refmod_file)
            end if
            close(51) 

         case(19)
            if(ios .ne. 0) then
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               write(51,*) 
               write(51,*) " There was a problem reading the number of values"
               write(51,*) " in the reference model file: ",trim(refmod_file)
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit
            elseif(indx .ne. nelem) then
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               write(51,*) 
               write(51,*) " The number of values specified in: ",trim(refmod_file)
               write(51,*) " is not equal to the number of values in: ",trim(sigfile)
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit
            end if

         case(20)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
            write(51,*) 
            write(51,*) " There was a problem reading line: ",indx
            write(51,*) " in the reference model file: ",trim(refmod_file)
            write(51,*) " Aborting ..."
            close(51)
            call crash_exit

         case(21)
            
            if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
            else
	       open(51,file='e4d.log',status='old',action='write',position='append')
            end if
            
            if(ios .ne. 0) then
               write(51,*) 
               write(51,*) " There was a problem reading the global starting constraint weight"
               write(51,*) " (beta) line in the inversion options file: ",trim(invfile)
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit
            else
               write(51,*)
               write(51,"(A30,8x,g10.4)") "Starting beta value: ",beta
               write(51,"(A30,9x,g10.4)") "Minimum objective reduction: ",del_obj
               write(51,"(A30,9x,g10.4)") "Beta reduction parameter: ",beta_red
               
            end if
     
            if(beta < 0) then
               write(51,*) 
               write(51,*) " The starting global constraint weight (beta) you provided"
               if(invi) then
                  write(51,*) " in ",trim(iinvfile)," is ",beta
               else
                  write(51,*) " in ",trim(invfile)," is ",beta
               end if
               write(51,*) " beta must be greater than zero."
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit 

            elseif(beta_red < 0 .or. beta_red > 1) then
               write(51,*) 
               write(51,*) " The beta reduction value (beta_red) you provided"
               if(invi) then
                  write(51,*) " in ",trim(iinvfile)," is ",beta_red
               else
                  write(51,*) " in ",trim(invfile)," is ",beta_red
               end if
               write(51,*) " beta_red must be greater than zero and less than one"
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit
 
            elseif(del_obj < 0 .or. del_obj > 1) then
               write(51,*) " The minimum objective function change ratio must be between zero and one"
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit
            end if 
            close(51)


         case(22)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
            if(ios .ne. 0) then
               write(51,*) 
               write(51,*) " There was a problem reading the target chi-squared value"
               if(invi) then
                  write(51,*) " in: ",trim(iinvfile)
               else
                  write(51,*) " in: ",trim(invfile)
               end if
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit

            elseif(norm_chi2 <= 0) then
               write(*,*)
               write(51,*) " The target chi-squared value you provided in ",trim(invfile)
               write(51,*) " is less than zero. The target chi-squared value must be"
               write(51,*) " greater than zero."
               write(51,*) " Aborting ..."
               close(51)
               call crash_exit

            else
               !write(51,*) 
               write(51,"(A30,8x,g10.4)") "Target chi-squared value: ",norm_chi2
            end if
            close(51)
         
 
            case(23)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               if(ios .ne. 0) then
                  write(51,*)
                  write(51,*) " There was a problem reading the maximum number of inner"
                  if(invi) then
                     write(51,*) " iterations in the inverse options file: ",trim(iinvfile)
                  else
                     write(51,*) " iterations in the inverse options file: ",trim(invfile)
                  end if
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
               else
                  !write(51,*)
                  write(51,"(A30,I10)") "Minimum # inner iterations: ",min_initer
                  write(51,"(A30,I10)") "Maximum # inner iterations: ",max_initer
                  !if(min_initer<30) then 
                  !   write(51,*) " WARNING"
                  !   write(51,*) " The recommended minimum number of inner iterations is 30"
                  !end if
               end if
               
               if(min_initer<0) then
                  write(51,*)
                  write(51,*) " The minimum number of inner iterations must be greater than zero"
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit;
               elseif(min_initer>max_initer) then
                  write(51,*) 
                  write(51,*) " The maximum number of inner iterations must be greater than the"
                  write(51,*) " minimum number of iterations"
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
               end if
               close(51)
               
            case(24)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               if(ios .ne. 0) then
                  write(51,*)
                  write(51,*) " There was a problem reading the max and min"
                  write(51,*) " conductivities in the inverse options file: ",trim(invfile)
                  write(51,*) " Aborting"
                  close(51)
                  call crash_exit
               else
                  !write(51,*) 
                  if(invi) then
                     write(51,"(A30,9x,g10.4)") "Maximum phase: ",max_sig
                     write(51,"(A30,9x,g10.4)") "Minimum phase: ",min_sig
                  else
                     write(51,"(A30,9x,g10.4)") "Maximum conductivity: ",max_sig
                     write(51,"(A30,9x,g10.4)") "Minimum conductivity: ",min_sig
                  end if
               end if

               if(min_sig < 0 .or. max_sig<0) then
                  write(51,*)
                  if(invi) then
                     write(51,*) " The phase limits must be greater than zero"
                     write(51,*) " Aborting ..."
                  else
                     write(51,*) " The conductivity limits must be greater than zero"
                     write(51,*) " Aborting ..."
                  end if
                  close(51)
                  call crash_exit

               elseif(max_sig < min_sig) then
                  write(51,*)
                  if(invi) then
                     write(51,*) " The maximum phase must be greater than the minimum phase"
                     write(51,*) " Aborting ..."
                  else
                     write(51,*) " The maximum conductivity must be greater than the minimum conductivy"
                     write(51,*) " Aborting ..."
                  end if
                  close(51)
                  call crash_exit
               end if
               close(51)


            case(25)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               if(ios .ne. 0) then
                  write(51,*)
                  write(51,*) " There was a problem reading the update option"
                  if(invi) then
                      write(51,*) " in the inverse options file: ",trim(iinvfile)
                  else
                     write(51,*) " in the inverse options file: ",trim(invfile)
                  end if
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
                

               elseif(up_opt==1) then
                  !write(51,*) 
                  write(51,"(A30,A)") "Update option: ","Line search"
               elseif(up_opt==2) then
                  !write(51,*) 
                  write(51,"(A30,9x,A)") "Update option: ","No line search"
               elseif(up_opt==3) then
                  write(51,"(A30,A)") "Update option: ", "No beta cooling"
               else
                  write(51,*) 
                  write(51,*) " The update option specified in the inverse options"
                  if(invi) then
                     write(51,*) " file: ",trim(iinvfile)," is ",up_opt
                  else
                     write(51,*) " file: ",trim(invfile)," is ",up_opt
                  end if
                  write(51,*) " Valid values are 1 for a line search and "
                  write(51,*) " 2 for no line search. "
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
               end if
               close(51)

            case(26)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               if(ios .ne. 0) then
                  write(51,*)
                  write(51,*) " There was a problem reading the number of line search "
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
                  close(51)
               else
                  write(51,*) 
                  write(51,"(A30,I10)") "# line search scalings: ",nlsp
               endif

               if(nlsp<1) then
                  write(51,*) 
                  write(51,*) " The number of line search scaling factors specified is ",nlsp
                  write(51,*) " There must be at least one line search scaling factor"
                  write(51,*) " Aborting ... "
                  close(51)
                  call crash_exit
               end if
               close(51)
               
            case(27)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               if(ios .ne. 0) then
                  write(51,*)
                  write(51,*) " There was a problem reading the line search scaling "
                  write(51,*) " factors."
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
               else
                  write(51,*) " Line search scaling factors: ",lsp(2:nlsp)
               end if
               
               do i=2,nlsp
                  if(lsp(i)<0) then
                     write(51,*)
                     write(51,*) " Line search scaling parameters must be greater than zero"
                     write(51,*) " Aborting ..."
                     close(51)
                     call crash_exit
                  end if
               end do
               do i=2,nlsp-1
                  if(lsp(i+1)<=lsp(i)) then
                     write(51,*)
                     write(51,*) " Each consequtive line search scaling parameter must be greater"
                     write(51,*) " than the last"
                     write(51,*) " Aborting ..."
                     close(51)
                     call crash_exit
                  end if
               end do
               close(51)

            case(28)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               if(ios .ne. 0) then
                  write(51,*)
                  write(51,*) " There was a problem reading the data culling options "
                  if(invi) then
                     write(51,*) " in the inverse options file: ",trim(iinvfile)
                  else
                     write(51,*) " in the inverse options file: ",trim(invfile)
                  end if
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit

               elseif((cull_flag .ne. 1) .and. (cull_flag .ne. 0)) then
                  write(51,*)
                  write(51,*) " The data cull flag provided is ",cull_flag
                  write(51,*) " Valid values are 0 for no data culling, or 1"
                  write(51,*) " for data culling. "
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit

               elseif(cull_flag == 1 .and. cull_dev <= 0) then
                  write(51,*)
                  write(51,*) " The data cull deviation provided is ",cull_dev
                  write(51,*) " Valid values are greater than zero"
                  write(51,*) " Aborting ..."
                  close(51)
                  call crash_exit
               else
                  if(cull_flag==0) then
                     write(51,"(A30)") "No data culling"
                  elseif(cull_flag==1) then
                     write(51,"(A30,8x,g10.4)") "Data culling st. deviation: ",cull_dev
                  end if
                  close(51)
               end if

            case(29)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               write(51,*)
               write(51,*)"****************************************************************"
               write(51,*)


            case(30)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               write(51,*) 
               write(51,*) " You specified a previous solution reference model"
               write(51,*) " for constraint block: ",indx
               write(51,*) " Previous reference models are only allowed for time"
               write(51,*) " lapse inversions "
               write(51,*) " Aborting ..."
               close(51)

			   write(*,*) 
               write(*,*) " You specified a previous solution reference model"
               write(*,*) " for constraint block: ",indx
               write(*,*) " Previous reference models are only allowed for time"
               write(*,*) " lapse inversions "
               
               call crash_exit
                
            case(31)
             if(im_fmm) then
               open(51,file='fmm.log',status='old',action='write',position='append')
             else
	       open(51,file='e4d.log',status='old',action='write',position='append')
             end if
               write(51,*) 
               if(invi) then
                  write(51,*) "The phase reference value for constraint block ",indx
               else
                  write(51,*) " The conductivity reference value for constraint block ",indx
               end if
               write(51,*) " is not greater than zero."
               write(51,*) " Aborting ..."
               close(51)

               write(*,*) 
               if(invi) then
                  write(*,*) "The phase reference value for constraint block ",indx
               else
                  write(*,*) " The conductivity reference value for constraint block ",indx
               end if
               write(*,*) " is not greater than zero."
               write(*,*) " Aborting ..."
               call crash_exit
                 
            case(32) 
               if(im_fmm) then
                  open(51,file='fmm.log',status='old',action='write',position='append')                                   
               else
                  open(51,file='e4d.log',status='old',action='write',position='append')
               end if
               write(51,*) "Structural metrics 12 and 13 are for joint inversion"
               write(51,*) "of travel time and resistivity data, and require both"
               write(51,*) "E4D and FMM to be running in inverse mode"
               if(.not. simulate_e4d) write(51,*) "E4D is not running ... Aborting"
               if(.not. simulate_fmm) write(51,*) "FMM is not running ... Aborting"
               close(51)
               write(*,*) "Joint inversion constraints have been specified, but not"
               write(*,*) "all of the required simulators are running. See log files."
               write(*,*) "Aborting ..."
               call crash_exit
               
            case(33)
               
               if(ios .ne. 0) then
                  if(im_fmm) then
                     open(52,file='fmm.log',status='old',action='write',position='append')     
                  else
                     open(52,file='e4d.log',status='old',action='write',position='append')
                  end if
                  write(52,*) "There was a problem reading the external constraint "
                  write(52,*) "file for constraint block ",indx
                  write(52,*) "Aborting ... "
                  write(*,*) "There was a problem reading the external constraint "
                  write(*,*) "file for constraint block ",indx
                  write(*,*) "Aborting ... "
                  close(52)
                  call crash_exit
               end if

            case(34)
               
               if(im_fmm) then
                  open(52,file='fmm.log',status='old',action='write',position='append')    
               else
                  open(52,file='e4d.log',status='old',action='write',position='append')
               end if
               write(*,*)
               write(52,*) "Cannot find the external constraint file specified "
               write(52,*) "for constraint block ",indx
               write(52,*) "Aborting ..."
               write(*,*)
               write(*,*) "Cannot find the external constraint file specified "
               write(*,*) "for constraint block ",indx
               write(*,*) "Aborting ..."
               close(52)
               call crash_exit
           

            case(35)
              
               if(im_fmm) then
                  open(52,file='fmm.log',status='old',action='write',position='append')    
               else
                  open(52,file='e4d.log',status='old',action='write',position='append')
               end if
               write(52,*) "            Constraint File:          "//trim(tmpstr)
               close(52)         


            case(36)
               
               if (mode==0) then
				   if(im_fmm) then
					  open(52,file='fmm.log',status='old',action='write',position='append')    
				   else
					  open(52,file='e4d.log',status='old',action='write',position='append')
				   end if
				   write(*,*)
				   write(52,*) "Reference model not used in inversion options file."
				   write(52,*) "The validity of this entry in e4d.inp was not checked."
				   write(*,*)
				   write(*,*) "Reference model not used in inversion options file."
				   write(*,*) "The validity of this entry in e4d.inp was not checked."

				   close(52)
				end if

          case default
    end select
  end subroutine check_inv_opts
  !_________________________________________________________________________  

  !_________________________________________________________________________
  subroutine read_survey
    implicit none
    logical :: exst
    integer :: i,j,junk
    real, dimension(4) :: etmp
    logical :: wdwarn = .true.

   
    inquire(file=trim(trim(efile)),exist=exst)
    if(.not. exst) then
       call check_inp(24,0)
    else
       call check_inp(24,1)
    end if

    open(10,file=efile,status='old',action='read')   
    read(10,*,IOSTAT=ios) ne;     call check_inp(12,junk)
    
    allocate(e_pos(ne,4))
    do i=1,ne
       read(10,*,IOSTAT=ios) junk,etmp; call check_inp(13,i)
       if(junk>ne) call check_inp(14,i)
       e_pos(junk,1:4)=etmp
    end do
    

    !!Read in the survey
    read(10,*,IOSTAT=ios) nm;     call check_inp(17,nm)
    allocate(dobs(nm),s_conf(nm,4),Wd(nm))
    if(i_flag) then 

       allocate(dobsi(nm),Wdi(nm))

       do i=1,nm         
          read(10,*,IOSTAT=ios) junk,s_conf(i,1:4),dobs(i),Wd(i),dobsi(i),Wdi(i); call check_inp(18,i)
          if(Wd(i) <= 0 .or. Wdi(i) <= 0) then
             Wd(i) = 1e15
             Wdi(i) = 1e15
             if(wdwarn) call check_inp(19,i)
             wdwarn = .false.
          end if
          
          do j=1,4
             if(s_conf(i,j)>ne .or. s_conf(i,j)<0) call check_inp(20,i)
             if(tank_flag .and. s_conf(i,j)==ne)   call check_inp(27,i)
          end do
          
          Wd(i)= 1/Wd(i)
          Wdi(i) = 1/Wdi(i)          
       end do
       close(10)
   
    else
       do i=1,nm
          read(10,*,IOSTAT=ios) junk,s_conf(i,1:4),dobs(i),Wd(i); call check_inp(18,i)
          
          if(Wd(i) <= 0) then
             Wd(i)=1e15
             if( wdwarn) call check_inp(19,i)
             wdwarn = .false.
          end if
          
          do j=1,4
             if(s_conf(i,j)>ne .or. s_conf(i,j)<0) call check_inp(20,i)
             if(tank_flag .and. s_conf(i,j)==ne)   call check_inp(27,i)
          end do
          Wd(i)= 1/Wd(i)
          
       end do
       close(10)
    end if
  end subroutine read_survey
  !_________________________________________________________________________

  !_________________________________________________________________________
  subroutine translate_electrodes
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

    call check_inp(15,junk)
    open(21,file=mshfile(1:mnchar)//'trn',status='old')
    read(21,*,IOSTAT=ios) xorig,yorig,zorig; call check_inp(16,junk)
    close(21) 
    e_pos(:,1) = e_pos(:,1)-xorig
    e_pos(:,2) = e_pos(:,2)-yorig
    e_pos(:,3) = e_pos(:,3)-zorig
    
  end subroutine translate_electrodes
  !_________________________________________________________________________

  !_________________________________________________________________________
  subroutine read_conductivity
    implicit none
    integer :: i,junk,npre,nchr
    logical :: exst
    real :: tsig
    if(allocated(sigma)) deallocate(sigma)
    if(mode .ne. 1) then
       inquire(file=trim(trim(sigfile)),exist=exst)
       if(.not.exst) then
          read(sigfile,*,IOSTAT=ios) tsig
          if(ios .ne. 0 ) then
             call check_inp(25,junk)
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
                open(51,file='e4d.log',status='old',action='write')
                write(51,*)
                write(51,*) ' Cannot find the element file : ',mshfile(1:npre)//'.ele'
                close(51)
                write(*,*)
                write(*,*) ' Cannot find the ele file : ',mshfile(1:npre)//'.ele'
                close(51)
                call crash_exit
             else
                open(10,file=mshfile(1:npre)//".ele",status='old',action='read')
                read(10,*) nsig
                close(10)
                allocate(sigma(nsig))
                sigma=tsig
                return
             end if
          end if
       end if
       open(10,file=sigfile,status='old',action='read')
       read(10,*,IOSTAT=ios) nsig; call check_inp(23,0) 
       
       allocate(sigma(nsig))
       
       if(i_flag) then
          if(allocated(sigmai)) deallocate(sigmai)
          allocate(sigmai(nsig))
          do i=1,nsig
             read(10,*,IOSTAT=ios) sigma(i),sigmai(i); call check_inp(23,i)
             if(sigma(i)<=0 .or. sigmai(i)<=0) call check_inp(26,i)
          end do
       else
          do i=1,nsig
             read(10,*,IOSTAT=ios) sigma(i); call check_inp(23,i)
             if(sigma(i)<=0) call check_inp(26,i)
          end do
          close(10)
       end if
    end if
       
  end subroutine read_conductivity
  !_________________________________________________________________________

  !_________________________________________________________________________
  subroutine check_for_list
    implicit none
    logical :: exst
    integer :: slen
    integer :: ios,i
    
    
    !!see if the file ends with .lst
    slen = len_trim(sigfile)
    if(sigfile(slen-3:slen) == ".lst") then
       !check to see if the file exists
       inquire(file=trim(trim(sigfile)),exist=exst)
       if(.not.exst) then 
          call check_inp(28,1)
       end if
       open(10,file=trim(sigfile),status='old',action='read')
       
       read(10,*,IOSTAT=ios) slen
       if(ios .ne. 0) then
          close(10)
          call check_inp(29,1)
       end if
       
       !!allocate an array of strings to hold the file names
       !!and read the file names
       allocate(tl_cfils(slen,2))
       do i=1,slen
          read(10,*,IOSTAT=ios) tl_cfils(i,1:2)
          if(ios .ne. 0) then
             close(10)
             call check_inp(30,i)
          end if
        
       end do
       close(10)

       !!check to make sure each file exists
       do i=1,slen
          inquire(file=trim(tl_cfils(i,1)),exist=exst)
          if(.not.exst) call check_inp(31,i)
       end do
       multi_forward = .true.
    end if
    
  end subroutine check_for_list
  !_________________________________________________________________________
  !_________________________________________________________________________
  subroutine crash_exit
    
    call MPI_BCAST(0,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call PetscFinalize(perr)
    stop
  end subroutine crash_exit
  !__________________________________________________________________________

  
  subroutine check_files
    implicit none
    
    integer, dimension(:,:), allocatable :: ipot
    character, dimension(:,:), allocatable :: output_str
    
    character*40 :: chk_sigfile
    character*80 :: dp_file
     
    integer :: i,junk,npot,ist,dp_flag,check
    logical :: exst,st
    
    junk=1
    read(10,*,IOSTAT=ios) mshfile
     
     !!Determine if the mesh file is a .cfg file or if meshfiles are provided
    mnchar = 0
    do i=1,40
       if(mshfile(i:i) == '.') then
          mnchar = i
          exit
       end if
    end do
    
    if(mshfile(mnchar+1:mnchar+3) == "cfg") then
       cfg_flag = .true. 
    end if
          
    inquire(file=trim(mshfile),exist=exst)
    open(51,file='e4d.log',status='old',action='write',position='append')
    if(.not.exst) then
	  if(cfg_flag) then
	     goto 6
	  else
		 goto 7
	  end if
   else
	  if(cfg_flag) then
		 write(51,*) " Mesh configuration file:          ",trim(mshfile)		 
	  else
		 write(51,*) " Mesh file:                        ",trim(mshfile)
	  end if
   end if
   close(51)
   
   ! if not a cfg file then check survey file, conductivity file, output options file
   if(.not.cfg_flag) then
     	read(10,*,IOSTAT=ios) efile;    call check_inp(3,junk)
        read(10,*,IOSTAT=ios) sigfile;  call check_inp(4,junk)
       
        chk_sigfile=sigfile
        chk_sigfile=lcase(chk_sigfile)
        if(trim(chk_sigfile)=="average") then
           ave_sig = .true.
        end if
        read(10,*,IOSTAT=ios) outfile;             call check_inp(5,junk)
        
        ! check for inversion options file entry in file (not complex inversion file)
        if (i_flag) then
			read(10,*,IOSTAT=ios) invfile, iinvfile;      call check_inp(6,junk)
			read(10,*,IOSTAT=ios) refmod_file;           call check_inp(7,junk)     
		else
			read(10,*,IOSTAT=ios) invfile;  call check_inp(6,junk)        
			read(10,*,IOSTAT=ios) refmod_file;  call check_inp(7,junk)  
		end if

        ! time-lapse optons
        check=0
        read(10,*,IOSTAT=ios) tl_file,check;   call check_inp(8,check) 
        
        
        ! check the survey file
        call read_survey
		
		
		!! check the conductivity file or list of files
		if(.not. ave_sig) then
			  multi_forward = .false.
			  call check_for_list
		end if

		if(.not. ave_sig .and. .not. multi_forward) then
		   call read_conductivity
		end if  
		
		!! check the output options file
		open(15,file=outfile,status='old',action='read')
		read(15,*,IOSTAT=ist) dp_flag; if(ist.ne.0) goto 11
		read(15,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
		read(15,*,IOSTAT=ist) npot   ; if(ist.ne.0) goto 13

		if(npot>0) then
		   allocate(ipot(npot,2))
		   do i=1,npot
			  read(15,*,IOSTAT=ist) ipot(i,1); if(ist.ne.0) goto 14
		   end do
		end if
	    close(15)
	    
	    
		!! check time-lapse options
	   if (tl_ly) then
		   open(21,file=trim(tl_file),status='old',action='read',IOSTAT=ios)  
		   read(21,*,IOSTAT=ios) ntl; call check_inp(9,junk)
		   allocate(tl_dfils(ntl),tlt(ntl))
	   
		   do i=1,ntl
			  read(21,*,IOSTAT=ios) tl_dfils(i), tlt(i); call check_inp(10,i)
			  inquire(file=trim(tl_dfils(i)),exist=exst)
			  if(.not. exst) call check_inp(11,i)
		   end do
	   end if
	   
       close(21)
	end if	
	
	close(51)
    close(10)
	    		
	
5	return

6	continue
	write(51,*) "  ERROR: Cannot find the mesh config. file name specified in e4d.inp"
	write(51,*) "  Aborting ..."
	close(51)
	write(*,*) "  ERROR: Cannot find the mesh config. file name specified in e4d.inp"
	write(*,*) "  Aborting ..."
	call crash_exit	
	return

7	continue
	write(51,*) "  ERROR: Cannot find the mesh file name specified in e4d.inp"
	write(51,*) "  Aborting ..."
	close(51)
	write(*,*) "  ERROR: Cannot find the mesh file name specified in e4d.inp"
	write(*,*) "  Aborting ..."
	call crash_exit	
	return
	      
11    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' The was a problem reading the first line in the output file: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' There was a problem reading the first line in the output file: ',trim(outfile)
      return

12    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) 'There was a problem reading the predicted data file name in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) 'The was a problem reading the predicted data file in: ',trim(outfile)
      return

13    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      close(51)
      write(*,*) 
      write(*,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      return

14    continue
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      close(51)
      write(*,*)
      write(*,*) ' There was a problem reading potential field index: ',i,' in: ',trim(outfile)
      return

  end subroutine check_files
  !__________________________________________________________________________


  character*40 function lcase(sinput)
	implicit none
	character*40, intent(inout) :: sinput
	integer :: slen, M, i
	
    slen = len(trim(sinput))
    if (slen > 1) then
        do i=1,slen 
			M=ichar(sinput(i:i))
			if( M .ge. 65 .AND. M .le. 90 ) then
				M = M + 32
				sinput(i:i) = char(M)
			end if
        end do
    end if 
    lcase=sinput
    return 
  end function lcase
  

end module input
