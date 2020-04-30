module exo_subs

  
  use subroutines
  implicit none
  include 'exodusII.inc'
  
  

  integer :: idexo,ierr
  
  contains

    !_______________________________________________________________
    subroutine init_exo
      implicit none

      !open the exodus file
      write(*,*) "Initializing ",trim(com_line(4))
      idexo=excre(trim(com_line(4)),EXCLOB,4,4,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error opening exodus file ",ierr
         stop
      end if
      
      !!initialize the file
      call expini(idexo,trim(com_line(4)),dim,nnods,nele,nzones,nnsets,0,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error initializing exodus file ",ierr
      end if
      
    end subroutine init_exo
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine open_exo
      implicit none
      real :: vers
      real :: dummy,tlast
      character*80 :: dbname
      character :: cdummy
      integer :: numess

      write(*,*) "Opening ",trim(com_line(4))
      idexo=exopen(trim(com_line(4)),EXWRIT,4,4,vers,ierr)
      if(ierr.ne.0) then
         write(*,*) 'Error opening ',trim(com_line(4))
         stop
      end if
      
      call exgini(idexo,dbname,dim,nnods,nele,nzones,nnsets,numess,ierr)

      call exinq(idexo,EXTIMS,n_times,dummy,cdummy,ierr)
      if(ierr.ne. 0) then
         write(*,*) "Error getting number of time steps"
         stop
      end if

      call exgtim(idexo,n_times,tlast,ierr)
      if(ierr.ne. 0) then
         write(*,*) "Error getting last time step value"
         stop
      end if
      write(*,*) "The original number of time steps is ",n_times
      write(*,*) "The last time recorded was ",tlast
      
    end subroutine open_exo
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine record_mesh
      implicit none
      integer :: ierr,i,j,z,count,nc
      character*(MXSTLN) :: coord_names(3)
      integer, dimension(nzones) :: idelb,numelb,numlnk,numatr,mkmaps
      character*(MXSTLN), dimension(nzones) :: namelb
      character*(MXSTLN), dimension(4) :: nameat
      integer, dimension(:,:), allocatable :: link
      real, dimension(:,:), allocatable :: atrib
      integer, dimension(:), allocatable :: ivec
      
      integer :: mzone
      
      coord_names(1) = "EASTING (m)"
      coord_names(2) = "NORTHING (m)"
      coord_names(3) = "ELEVATION (m)"

      !!record the node positions
      write(*,*) "Recording node locations"
      call expcor(idexo,nods(:,1),nods(:,2),nods(:,3),ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error recording the nodal coordinates ",ierr 
         stop
      end if
      
     
      call expcon(idexo,coord_names,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error recording coordinate names ",ierr
         stop
      end if

      !!record the node sets
      do i=1,nnsets
         write(*,*) 'Recording node set: ',nset_ids(i)
         !build the node id list for this set
         allocate(ivec(nset_sizes(i)))
         nc=0
         do j=1,nnods
            if(n_type(j,2) == nset_ids(i)) then
               nc=nc+1
               ivec(nc)=j
               if(nc==nset_sizes(i)) exit
            end if
         end do
       
         !call expnp(idexo,nset_ids(i),nset_sizes(i),nset_sizes(i),ierr)
         call expnp(idexo,nset_ids(i),nset_sizes(i),0,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording the nodes set parameters for set: ",i
         end if

         call expns(idexo,nset_ids(i),ivec,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording the nodes set for set: ",i
         end if
         
         deallocate(ivec)
      end do


      !!record the element block info
      mzone=minval(ele(:,5))
      do i=1,nzones
         idelb(i)=mzone+i-1
         namelb(i)='TETRAHEDRON'
         numlnk(i)=4   !4 connections per element
         numatr(i)=4  !4 neighbors per element
      end do
      numelb=0
      mkmaps=0
      do i=1,nele
         do j=1,nzones
            if(ele(i,5)==idelb(j)) then
               numelb(j)=numelb(j)+1
               goto 10
            end if
         end do
10       continue
      end do     

      call expclb(idexo,idelb,namelb,numelb,numlnk,numatr,mkmaps,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error recording element block info",ierr
         stop
      end if
      
      !record element blocks
      do z=1,nzones
         write(*,*) "Recording connections for zone ",z
         if(allocated(link)) deallocate(link)
         allocate(link(numlnk(z),numelb(z)))
         count=0
         do i=1,nele
            if(ele(i,5)==z) then
               count=count+1
               do j=1,numlnk(z)
                  link(j,count)=ele(i,j)
               end do
            end if
         end do
         call expelc(idexo,z,link,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording element connections",z,ierr
            stop
         end if
      end do
      if(allocated(link)) deallocate(link)

      !!update the file
      call exupda(idexo,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error updating exodus file",ierr
         stop
      end if

      return


      !!write the element attributes sig1,sig2,sig3,neigh1,neigh2,neigh3,neigh4
      !nameat(1)= 'sigma_x'
      !nameat(2)= 'sigma_y'
      !nameat(3)= 'sigma_z'
      nameat(1)= 'neighbor_1'
      nameat(2)= 'neighbor_2'
      nameat(3)= 'neighbor_3'
      nameat(4)= 'neighbor_4'

      do z=1,nzones
         write(*,*) "Recording attributes for zone ",z
         if(allocated(atrib)) deallocate(atrib)
         allocate(atrib(numatr(z),numelb(z)))
         count=0
         do i=1,nele
            if(ele(i,5)==z) then
               count=count+1
               !do j=1,3
               !   atrib(j,count) = sig(i,j)
               !end do
               do j=1,4
                  atrib(j,count)= neigh(i,j)
               end do
            end if
         end do
         call expeat(idexo,z,atrib,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording element attributes",z,ierr
            stop
         end if

         !call expean(idexo,z,nameat,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording attribute name",z,ierr
            stop
         end if
      end do
      if(allocated(atrib)) deallocate(atrib) 

      !!update the file at this point
      call exupda(idexo,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error updating exodus file",ierr
         stop
      end if
     
    end subroutine record_mesh
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine enter_data
      !!enters the data specified on the command line
      implicit none
      character*80 :: cfil
      integer :: ndat,ioserr,ii,i,jnk,ncol

      if(trim(com_line(1))=='-f' .or. trim(com_line(1))=='-file' .or. &
         trim(com_line(1))=='-af' .or. trim(com_line(1)) =='-afile') then
       
         nfil = 1
         cfil = com_line(3)
         open(25,file=trim(com_line(3)),status='old',action='read')
         read(25,*,IOSTAT=ioserr) ndat,ncol
         if(ioserr .ne. 0) goto 1001
      
         if(ndat==nnods .and. (ncol == 1 .or. ncol == 2)) then
           
            if(allocated(pot)) deallocate(pot)
            allocate(pot(nnods,ncol))
            if(ncol==1) then
               do i=1,nnods
                  read(25,*,IOSTAT=ioserr) pot(i,1)
                  if(ioserr .ne. 0) goto 1002
               end do
               close(25)
            elseif(ncol==2) then
               do i=1,nnods
                  read(25,*,IOSTAT=ioserr) pot(i,1:2)
                  if(ioserr .ne. 0) goto 1002
               end do
               close(25)
            end if

            !if(.not.allocated(nods)) call get_nodes
            call record_potential(nfil+n_times,ncol)
            deallocate(pot)

         elseif(ndat==nele .and. (ncol == 1 .or. ncol == 2)) then
            
            if(allocated(sig)) deallocate(sig)
            allocate(sig(nele,ncol))
            if(ncol==1) then
               do i=1,nele
                  read(25,*,IOSTAT=ioserr) sig(i,1)
                  if(ioserr .ne. 0) goto 1001
               end do
               close(25)
            elseif(ncol==2) then
               do i=1,nele
                  read(25,*,IOSTAT=ioserr) sig(i,1:2)
                  if(ioserr .ne. 0) goto 1001
               end do
               close(25)
            end if
            call record_sigma_1(nfil+n_times,ncol)

         end if

         
      elseif(trim(com_line(1))=='-l' .or. trim(com_line(1))=='-list' .or. &
          trim(com_line(1))=='-al' .or. trim(com_line(1))=='-alist') then
        
         open(25,file=trim(com_line(3)),status='old',action='read')
         read(25,*,IOSTAT=ioserr) nfil
         do ii=1,nfil
            read(25,*) cfil,time_stamp
            open(26,file=cfil,status='old',action='read')
            read(26,*) ndat,ncol

            if(ndat==nnods .and. (ncol==1 .or. ncol==2)) then
               if(allocated(pot)) deallocate(pot)
               allocate(pot(nnods,ncol))
               if(ncol==1) then
                  do i=1,nnods
                     read(26,*,IOSTAT=ioserr) pot(i,1)
                     if(ioserr .ne. 0) goto 1002
                  end do
               elseif(ncol==2) then
                  do i=1,nnods
                     read(26,*,IOSTAT=ioserr) pot(i,1:2)
                     if(ioserr .ne. 0) goto 1002
                  end do
               end if
               close(26)
               !if(.not.allocated(nods)) call get_nodes
               call record_potential(nfil+n_times,ncol)
               !call record_potential(25)
               close(26)
               
            elseif(ndat==nele .and. (ncol==1 .or. ncol==2)) then
               if(allocated(sig)) deallocate(sig)
               allocate(sig(nele,ncol))
               if(ncol==1) then
                  do i=1,nele
                     read(26,*,IOSTAT=ioserr) sig(i,1)
                     if(ioserr .ne. 0) goto 1001
                  end do
               elseif(ncol==2) then
                  do i=1,nele
                     read(26,*,IOSTAT=ioserr) sig(i,1:2)
                     if(ioserr .ne. 0) goto 1001
                  end do
               end if
               close(26)
           
               call record_sigma_1(ii+n_times,ncol)
            end if

         end do
      end if

      return

      !!update the file at this point
      call exupda(idexo,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error updating exodus file",ierr
         stop
      end if

1001  write(*,*) "THERE WAS A PROBLEM READING CONDUCTIVITY VALUE: ",i,nele
      write(*,*) "IN FILE: ",cfil
      stop

1002  write(*,*) "THERE WAS A PROBLEM READING POTENIAL VALUE: ",i,nnods
      write(*,*) "IN FILE: ",cfil
      stop


    end subroutine enter_data
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine record_potential(fn,rc)
      implicit none
      integer :: fn,rc
      integer :: ierr,i,j,nej,k,nct
      character*(MXSTLN) :: names(3)
      integer, dimension(:,:), allocatable :: isevok
      real, dimension(nnods) :: zsig
      
      if(fn==1) then
         !!record the number of node variables that will be written
         !!which is 1 in this case
         call expvp(idexo,"e",2,ierr)
         call expvp(idexo,"n",2,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording the number of node variables"
            write(*,*) "in record_pot"
            stop
         end if
         
         !!record the variable name
         names(1) = trim("real_conductivity")
         names(2) = trim("complex_conductivity")
         call expvan(idexo,"e",2,names(1:2),ierr)

         names(1) = trim("real_potential")
         names(2) = trim("complex_potenial")
         call expvan(idexo,"n",2,names,ierr)

         if(ierr .ne. 0) then
            write(*,*) "Error recording the variable names"
            write(*,*) "in record_pot"
            stop
         end if
      end if
        


      !!record the time step 
      call exptim(idexo,fn,time_stamp,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error recording the time_step"
         write(*,*) "in record_sigma"
         stop
      end if 
      
!!$      do i=1,nnsets    
!!$         nct=0
!!$         do j=1,nnods
!!$            if(n_type(j,2)==nset_ids(i)) then
!!$               nct=nct+1
!!$               zsig(nct)=sig(j)
!!$               if(nct==nset_sizes(i)) exit
!!$            end if
!!$         end do
!!$         write(*,*) idexo,fn,i,nset_ids(i),nset_sizes(i),maxval(zsig(1:nset_sizes(i)))
!!$         call expnsv(idexo,fn,i,nset_ids(i),nset_sizes(i),zsig(1:nset_sizes(i)),ierr)
!!$      end do

      !!record the potential
      if(rc==1) then
         call expnv(idexo,fn,1,nnods,pot(:,1),ierr)
      elseif(rc==2) then
         call expnv(idexo,fn,1,nnods,pot(:,1),ierr)
         call expnv(idexo,fn,2,nnods,pot(:,2),ierr)
      end if
      if(ierr .ne. 0) then
         write(*,*) "Error recording the potential"
         stop
      end if
      call exupda(idexo,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error updating exodus file",ierr
         stop
      end if
    end subroutine record_potential
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine record_sigma_1(fn,rc)
      !!records conductivity values
      implicit none
      integer :: i,j,k,nej,fn,rc
      integer, dimension(:,:), allocatable :: isevok
      real, dimension(nele) :: zsig
      character*(MXSTLN) :: names(3)
     
      if(fn==1) then

         !!enter the number of element variables (1)
         call expvp(idexo,"e",2,ierr)
         call expvp(idexo,"n",2,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording the number of element variables"
            write(*,*) "in record_sigma_1"
            stop
         end if
         
         !record the variable name

       names(1) = trim("real_conductivity")
         names(2) = trim("complex_conductivity")
         call expvan(idexo,"e",2,names(1:2),ierr)
  
         names(1) = trim("real_potential")
         names(2) = trim("complex_potenial")
         call expvan(idexo,"n",2,names,ierr)

         if(ierr .ne. 0) then
            write(*,*) "Error recording the variable name"
            write(*,*) "in record_sigma_1"
            stop
         end if
      end if

      !!record the time step 
      call exptim(idexo,fn,time_stamp,ierr)
      if(ierr .ne. 0) then
         write(*,*) "Error recording the time_step"
         write(*,*) "in record_sigma"
         stop
      end if
 
      if(fn==1) then
         !!write the element variable truth table
         !!allocate(isevok(nzones,3))
         allocate(isevok(nzones,2))
         isevok=1
         !!call expvtt(idexo,nzones,3,isevok,ierr)
         call expvtt(idexo,nzones,2,isevok,ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording the element truth table"
            write(*,*) "in record_sigma"
            stop
         end if
      end if

      !!read and record sigma
      i=1
      write(*,*)
      write(*,*) "RECORDING CONDUCTIVITY FOR TIME ",time_stamp
      do j=1,nzones
         write(*,*) 'recording conductivity for zone ',j
         nej=0
         !!count the number of elements in this zone
         do k=1,nele
            if(ele(k,5)==j) then
               nej=nej+1
               zsig(nej)=sig(k,1)
            end if
         end do
         call expev(idexo,fn,i,j,nej,zsig(1:nej),ierr)
         if(ierr .ne. 0) then
            write(*,*) "Error recording sigma",i,j,nej
            write(*,*) "in record_sigma_1"
            stop
         end if
      end do

      if(rc==2) then
         i=2
         write(*,*)
         write(*,*) "RECORDING COMPLEX CONDUCTIVITY FOR TIME ",time_stamp
         do j=1,nzones
            write(*,*) 'recording conductivity for zone ',j
            nej=0
            !!count the number of elements in this zone
            do k=1,nele
               if(ele(k,5)==j) then
                  nej=nej+1
                  zsig(nej)=sig(k,2)
               end if
            end do
            call expev(idexo,fn,i,j,nej,zsig(1:nej),ierr)
            if(ierr .ne. 0) then
               write(*,*) "Error recording sigma",i,j,nej
               write(*,*) "in record_sigma_1"
               stop
            end if
         end do
         
      end if

      
    end subroutine record_sigma_1
    !_______________________________________________________________


    !_______________________________________________________________
!!$    subroutine record_sigma
!!$      implicit none
!!$      integer :: ierr,i,j,nej,k
!!$      character*(MXSTLN) :: names(3)
!!$      integer, dimension(:,:), allocatable :: isevok
!!$      real, dimension(nele) :: zsig
!!$
!!$      !!record the number of element variables that will be written
!!$      !!which is 3 in this case, 1 for each direction
!!$      !!call expvp(idexo,"e",3,ierr)
!!$      call expvp(idexo,"e",1,ierr)
!!$      if(ierr .ne. 0) then
!!$         write(*,*) "Error recording the number of element variables"
!!$         write(*,*) "in record_sigma"
!!$         stop
!!$      end if
!!$
!!$      !!record the variable names
!!$      names(1)="sigma"
!!$      names(2)="sigma_y"
!!$      names(3)="sigma_z"
!!$      !!call expvan(idexo,"e",3,names,ierr)
!!$      call expvan(idexo,"e",1,names(1),ierr)
!!$      if(ierr .ne. 0) then
!!$         write(*,*) "Error recording the variable names"
!!$         write(*,*) "in record_sigma"
!!$         stop
!!$      end if
!!$
!!$      !!record the time step 
!!$      call exptim(idexo,1,0.0,ierr)
!!$      if(ierr .ne. 0) then
!!$         write(*,*) "Error recording the time_step"
!!$         write(*,*) "in record_sigma"
!!$         stop
!!$      end if
!!$
!!$      !!write the element variable truth table
!!$      !!allocate(isevok(nzones,3))
!!$      allocate(isevok(nzones,rc))
!!$      isevok=1
!!$      !!call expvtt(idexo,nzones,3,isevok,ierr)
!!$      call expvtt(idexo,nzones,2,isevok,ierr)
!!$      if(ierr .ne. 0) then
!!$         write(*,*) "Error recording the element truth table"
!!$         write(*,*) "in record_sigma"
!!$         stop
!!$      end if
!!$
!!$      !!write the sigma values
!!$      !do i=1,3
!!$         i=1
!!$         do j=1,nzones
!!$            write(*,*) 'recording sigma ',i,' for zone ',j
!!$            nej=0
!!$            !!count the number of elements in this zone
!!$            do k=1,nele
!!$               if(ele(k,5)==j) then
!!$                  nej=nej+1
!!$                  zsig(nej)=sig(k,1)
!!$               end if
!!$            end do
!!$            call expev(idexo,1,i,j,nej,zsig(1:nej),ierr)
!!$            if(ierr .ne. 0) then
!!$               write(*,*) "Error recording sigma",i,j,nej
!!$               write(*,*) "in record_sigma"
!!$               stop
!!$            end if
!!$         end do
!!$      !end do
!!$      
!!$       !!update the file at this point
!!$      call exupda(idexo,ierr)
!!$      if(ierr .ne. 0) then
!!$         write(*,*) "Error updating exodus file",ierr
!!$         stop
!!$      end if
!!$
!!$
!!$    end subroutine record_sigma
    !_______________________________________________________________


    !_______________________________________________________________
    subroutine close_exo
      implicit none
      integer :: ierr
      call exclos(idexo,ierr)
      write(*,*) ierr
    end subroutine close_exo
    !_______________________________________________________________
end module exo_subs
