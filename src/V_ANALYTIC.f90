module v_analytic

  use vars
  use input
  use master
  use output

  implicit none
  real :: se,s_ave
  integer :: dp_flag,n_pot
  character*80 :: dp_file
  integer, dimension(:), allocatable :: poti
  real, dimension(:), allocatable :: an_pot
contains
  
    !______________________________________________________________________________
    subroutine compute_analytic
      implicit none
      logical :: stat
      integer :: i
     
      n_pot=0
      call get_out_opts(stat)
      if(.not.stat) return
      
      !call write_an_header
    
      if(dp_flag .eq.1 .or. n_pot>1) then
 
         call read_nodes(stat)
         if(.not.stat) call crash_exit
      
         call read_survey
         call translate_electrodes
         call read_conductivity

         call get_level
         
         if(dp_flag==1) then
            call build_an_dpred
            if(dp_flag==1) then
               call write_an_dpred
               call write_an_srv
            end if
         end if

         if(n_pot>0) then
            do i=1,n_pot
               if(poti(i) .le. nm) then
                  call build_pot_an(poti(i))
                  call write_pot_an(poti(i))
               end if
            end do
         end if
      end if

      if(allocated(poti)) deallocate(poti)
      if(allocated(an_pot)) deallocate(an_pot)

      open(51,file='e4d.log',status='old',action='write',position="append")
      write(51,*)
      write(51,*) "Computation complete."
      close(51)
      return
    end subroutine compute_analytic
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_ave_sig
      implicit none
      logical :: stat
      integer :: i,cc
      real :: gfam,gfan,gfbm,gfbn,ts,ssig,swt,gft,phz
      real :: mxsig,mnsig,mxphz,mnphz,phzi
      real, parameter :: pi = 3.14159265359

      call read_nodes(stat)
      if(.not.stat) call crash_exit
    
      call read_elements(stat)
      if(.not.stat) call crash_exit

      if(.not.allocated(sigma)) then
         allocate(sigma(nelem))
         nsig=nelem
      end if

      if(i_flag .and. .not. allocated(sigmai)) allocate(sigmai(nelem))

      call get_level 

      ssig=0
      swt=0
      mnsig=1e30
      mxsig=0
      
      do i=1,nm
         call get_gf(gfam,s_conf(i,1),s_conf(i,3))
         call get_gf(gfan,s_conf(i,1),s_conf(i,4))
         call get_gf(gfbm,s_conf(i,2),s_conf(i,3))
         call get_gf(gfbn,s_conf(i,2),s_conf(i,4))
         gft = (gfam- gfan - (gfbm - gfbn))/(4*pi)
        
         ts=gft/dobs(i)
         if(ts>0) then
            if(ts>mxsig) mxsig=ts
            if(ts<mnsig) mnsig=ts
            ssig=ssig+ts/Wd(i)
            swt=swt+1/Wd(i)
         end if
        
      end do
      mean_sig=ssig/swt
      sigma=mean_sig


      if(i_flag) then
         phz=0
         mnphz=1e30
         mxphz=0
         cc=0
         do i=1,nm
            phzi= -dobsi(i)
           
            if(phzi>0) then
               if(phzi>mxphz) mxphz=phzi
               if(phzi<mnphz) mnphz=phzi
               phz=phz+phzi
               cc=cc+1
            end if
         end do
         mean_phase=phz/cc
         sigmai=tan(mean_phase)*sigma
     
      end if
      
      call write_sigma

      open(51,file='e4d.log',status='old',action='write',position='append'); 
      write(51,"(A36,g10.5)") "  Minimum apparent conductivity :   ",mnsig
      write(51,"(A36,g10.5)") "  Maximum apparent conductivity :   ",mxsig
      write(51,"(A36,g10.5)") "  Using conductivity            :   ",mean_sig
      close(51)
      
      if(i_flag) then
         open(51,file='e4d.log',status='old',action='write',position='append'); 
         write(51,"(A36,g10.5)") "  Minimum apparent phase        :   ",mnphz
         write(51,"(A36,g10.5)") "  Maximum apparent phase        :   ",mxphz
         write(51,"(A36,g10.5)") "  Using phase                   :   ",mean_phase
         close(51)
      end if
   
      
    end subroutine get_ave_sig
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine write_pot_an(id)
      implicit none
      integer :: id
      integer :: i
      character*80 :: fstr
      
      write(fstr,"(A,I0)") 'an_potential.',id
      
      open(51,file='e4d.log',status='old',action='write',position="append")
      write(51,*) "  Writing analytic potential field ",trim(fstr)
      close(51)

      open(10,file=trim(fstr),status='replace',action='write')
      write(10,*) nnodes,1
      do i=1,nnodes
         write(10,*) an_pot(i)
      end do
      close(10)

    end subroutine write_pot_an
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine build_pot_an(id)
      implicit none
      integer :: id
      integer :: i,ei
      real :: ixp,iyp,izp,ixn,iyn,izn
      real, dimension(nnodes) :: r
      real, parameter :: pi = 3.14159265359
      
      if(.not. allocated(an_pot)) allocate(an_pot(nnodes))
  
      if(s_conf(id,1)>0) then
         ixp=e_pos(s_conf(id,1),1)
         iyp=e_pos(s_conf(id,1),2)
         izp=e_pos(s_conf(id,1),3)
      end if
      if(s_conf(id,2)>0) then
         ixn=e_pos(s_conf(id,2),1)
         iyn=e_pos(s_conf(id,2),2)
         izn=e_pos(s_conf(id,2),3)
      end if
      an_pot=0
      if(s_conf(id,1)>0) then
         r=sqrt( (ixp-nodes(:,1))**2 + (iyp-nodes(:,2))**2 + (izp-nodes(:,3))**2)
         r=r+1e-15
         an_pot = an_pot + 1/(4*pi*r*s_ave);

         izp=2*se-izp
         r=sqrt( (ixp-nodes(:,1))**2 + (iyp-nodes(:,2))**2 + (izp-nodes(:,3))**2)
         r=r+1e-15
         an_pot = an_pot + 1/(4*pi*r*s_ave);
      end if

      if(s_conf(id,2)>0) then
         r=sqrt( (ixn-nodes(:,1))**2 + (iyn-nodes(:,2))**2 + (izn-nodes(:,3))**2)
         r=r+1e-15


         an_pot = an_pot - 1/(4*pi*r*s_ave);

         izn=2*se-izn
         r=sqrt( (ixn-nodes(:,1))**2 + (iyn-nodes(:,2))**2 + (izn-nodes(:,3))**2)
         r=r+1e-15
         an_pot = an_pot - 1/(4*pi*r*s_ave);
      end if

    end subroutine build_pot_an
    !______________________________________________________________________________
    !______________________________________________________________________________
    subroutine write_an_dpred
      implicit none
      integer :: i
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) "  Writing analytic measurements to ",trim(dp_file)
      write(*,*) "  Writing analytic measurements to ",trim(dp_file)
      close(51)

      open(10,file=trim(dp_file),status='replace',action='write')
      write(10,*) nm
      do i=1,nm
         write(10,*) i,s_conf(i,1),s_conf(i,2),s_conf(i,3),s_conf(i,4),dobs(i),dpred(i)
      end do
      close(10)
    end subroutine write_an_dpred
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine write_an_srv
      implicit none
      integer :: i
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) 
      write(51,*) "  Writing analytic survey to ",trim(sigfile),".srv"
      write(*,*) "  Writing analytic survey to ",trim(sigfile),".srv"
      
      close(51)

      open(10,file=trim(sigfile)//".srv",status='replace',action='write')
      write(10,*) ne
      do i=1,ne
         write(10,"(I7,4F10.3,I5)") i,e_pos(i,1)+xorig,e_pos(i,2)+yorig,e_pos(i,3)+zorig,e_pos(i,4)
      end do
      write(10,*)
      write(10,*) nm
      do i=1,nm
         write(10,"(5I7,2g15.6)") i,s_conf(i,1),s_conf(i,2),s_conf(i,3),s_conf(i,4),dpred(i),0.05*abs(dpred(i))
      end do
      close(10)
    end subroutine write_an_srv
    !______________________________________________________________________________


    !______________________________________________________________________________
    subroutine build_an_dpred
      implicit none
      real :: gfam,gfan,gfbm,gfbn,gft
      integer :: i
      real, parameter :: pi = 3.14159265359
      
     
      if(allocated(dpred)) deallocate(dpred)
      allocate(dpred(nm))
      
      do i=1,nm
         call get_gf(gfam,s_conf(i,1),s_conf(i,3))
         call get_gf(gfan,s_conf(i,1),s_conf(i,4))
         call get_gf(gfbm,s_conf(i,2),s_conf(i,3))
         call get_gf(gfbn,s_conf(i,2),s_conf(i,4))

         gft = (gfam- gfan - (gfbm - gfbn))/(4*pi)
         dpred(i) = gft/s_ave
      end do
      
      
    end subroutine build_an_dpred
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_gf(gf,ei,ep)
      implicit none
      real :: gf
      integer :: ei,ep
      real :: xi,yi,zi,zii,xp,yp,zp,r,ri
     
      gf=0;
      if(ei==0 .or. ep==0) return
  
      xi=e_pos(ei,1); yi=e_pos(ei,2); zi=e_pos(ei,3); zii=2*se-e_pos(ei,3)
      
      xp=e_pos(ep,1); yp=e_pos(ep,2); zp=e_pos(ep,3)
      
      r=sqrt(  (xi-xp)**2 + (yi-yp)**2 + (zi-zp)**2);
      ri=sqrt( (xi-xp)**2 + (yi-yp)**2 + (zii-zp)**2);
      
      gf = 1/r + 1/ri;
      
    end subroutine get_gf
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_level
      implicit none
      integer :: i,cc
      logical :: flat = .true.
      real :: mx,mn
      cc = 0
      se = 0

      mx=nodes(1,3);
      mn=nodes(1,3);
      
      do i=2,nnodes
         if(nbounds(i)==1) then
            if(nodes(i,3)>mx) mx=nodes(i,3)
            if(nodes(i,3)<mn) mn=nodes(i,3)
         end if
      end do
      !se=.5*(mx+mn)
      se = mx
     
      open(51,file='e4d.log',status='old',action='write',position='append');
      write(51,"(A,I10.10)") "  Number of nodes:                  ",nnodes
      write(51,*) 
      write(51,*) " Maximum surface node elevation:   ",mx+zorig
      write(51,*) " Minimum surface node elevation:   ",mn+zorig
      write(51,*) " Using surface elevation of    :   ",se+zorig
      if(mx .ne. se .or. mn .ne. se) then
         write(51,*) "  !!! WARNING: It appears there may be some surface variability."
         write(51,*) "  !!! The analytic solution requires a flat surface."
      end if
      close(51)
      
      if(ave_sig) return

      mx=sigma(1);
      mn=sigma(1);
      
      do i=2,nelem
            if(sigma(i)>mx) mx=sigma(i)
            if(sigma(i)<mn) mn=sigma(i)
      end do
      s_ave=.5*(mx+mn)

      open(51,file='e4d.log',status='old',action='write',position='append');
      write(51,*) 
      write(51,*) " Maximum conductivity:             ",mx
      write(51,*) " Minimum conductivity:             ",mn
      write(51,*) " Using conductivity  :             ",s_ave
      if(mx .ne. s_ave .or. mn .ne. s_ave) then
         write(51,*) "  !!!WARNING: it appears there may be some conductive variability"
         write(51,*) "  !!!The analytic solution requires homogeneous conductivity"
      end if
      close(51)

    end subroutine get_level      
    !______________________________________________________________________________
  
    !______________________________________________________________________________
    subroutine write_an_header
      implicit none
           
      open(51,file='e4d.log',status='old',action='write')
      write(51,*) "RUNNING IN ANALYTIC FORWARD SOLUTION MODE"
      write(51,*) "MESH FILE = ",mshfile
      write(51,*) "SURVEY FILE = ",efile
      write(51,*) "CONDUCTIVITY FILE = ",sigfile
      write(51,*) "OUTPUT OPTIONS FILE = ",outfile
      close(51);
    end subroutine write_an_header
    !______________________________________________________________________________

    !______________________________________________________________________________
    subroutine get_out_opts(chck)
      implicit none
      logical :: chck
      logical :: exst
      integer :: ist,i

      chck=.true.
      inquire(file=trim(outfile),exist=exst)
      if(.not. exst) goto 10

      open(10,file=trim(outfile),status='old',action='read')
      read(10,*,IOSTAT=ist) dp_flag; if(ist.ne.0) goto 11
      read(10,*,IOSTAT=ist) dp_file; if(ist.ne.0) goto 12
      read(10,*,IOSTAT=ist) n_pot;   if(ist.ne.0) goto 13
      if(n_pot>0) then
         allocate(poti(n_pot))
         do i=1,n_pot
            read(10,*,IOSTAT=ist) poti(i); if(ist.ne.0) goto 14
         end do
      end if
      close(10)
      
      return

10    continue
      chck = .false.
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' Cannot find the output options file: ',trim(outfile)
      write(51,*) ' Aborting ...'
      close(51)
      write(*,*) ' Cannot find the output options file: ',trim(outfile)
      write(*,*) ' Aborting ...'
      return 
      
11    continue
      chck = .false.
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' The was a problem reading the first line in the output file: ',trim(outfile)
      write(51,*) ' Aborting ...'
      close(51)
      write(*,*) ' The was a problem reading the first line in the output file: ',trim(outfile)
      write(*,*) ' Aborting ...'
      return

12    continue
      chck = .false.
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' The was a problem reading the predicted data file in: ',trim(outfile)
      write(51,*) ' Aborting ...'
      close(51)
      write(*,*) ' The was a problem reading the predicted data file in: ',trim(outfile)
      write(*,*) ' Aborting'
      return

13    continue
      chck = .false.
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' There was a problem reading the number of potential fields to write in: ',trim(outfile)
      write(51,*) ' Aborting ...'
      close(51)
      write(*,*) ' The was a problem reading the number of potential fields to write in: ',trim(outfile)
      write(*,*) ' Aborting ...'
      return

14    continue
      chck = .false.
      open(51,file='e4d.log',status='old',action='write',position='append')
      write(51,*) ' The was a problem reading potential index: ',i,' in: ',trim(outfile)
      write(51,*) ' Aborting ...'
      close(51)
      write(*,*) ' The was a problem reading potential index: ',i,' in: ',trim(outfile)
      write(*,*) ' Aborting ...'

      return


    end subroutine get_out_opts
    !______________________________________________________________________________

end module v_analytic
