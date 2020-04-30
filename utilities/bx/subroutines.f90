module subroutines

  implicit none
  character*80, dimension(5) :: com_line
  integer, dimension(5) :: com_lengths
  character*80 :: tet_prefix
  integer :: nnods,nele,nneigh,nface,dim,nsig,nzones,nnsets,nzvals,i,nfil
  integer :: n_times,nchar
  integer , dimension(:,:), allocatable :: ele,neigh,face,n_type
  real, dimension(:,:), allocatable :: nods
  real, dimension(:,:), allocatable :: sig,pot
  integer, dimension(:), allocatable :: zvals,nset_ids,nset_sizes
  real :: tx,ty,tz,time_stamp
  logical :: newexo = .true.
  character*80 :: tfname

  contains
    
    !_______________________________________________________________
    subroutine get_args
      !get the tetgen output file prefix from the command line
      implicit none
      logical :: lnode,lele,lneigh,lsig,lface,fabort,ltrans
      integer :: i, status
   
      n_times = 0
      do i=1,4
         call get_command_argument(i,com_line(i),com_lengths(i),status)
        
         if(status > 0) then
            write(*,"(A32,I4)") "FAILED TO READ ARGUMENT",i
            call print_help
            stop
         elseif(status<0) then
            write(*,"(A9,I4,A30)") "ARGUMENT ",i," MUST BE 80 CHARACTERS OR LESS"
            call print_help
            stop
         end if
         
         if(trim(com_line(i))=='-help'.or.trim(com_line(i))=='--help') then
            call print_help
            stop
         end if

         select case(i)
            case(1)
               call check_option(status)
               if(status<0) then
                  write(*,*) 'OPTION ',trim(com_line(i)),' IS INVALID'
                  call print_help
                  stop
               end if
         
            case(2)
               select case(trim(com_line(1)))
                  case('-al')
                  case('-alist')
                  case('-af')
                  case('-afile')
                  case default
                     call check_mesh(status)
                     if(status .ne. 0) stop
                  end select

            case(3)
               call check_input(status)
               if(status .ne. 0) stop

            case(4)
               call check_output(status)
               if(status .ne. 0) stop
           
         end select
      end do

      if(trim(com_line(1)) == "-f" .or. (com_line(1)) == "-file" .or. &
           trim(com_line(1)) == "-af" .or. trim(com_line(1)) == "-afile") then
         
         call get_command_argument(5,com_line(5),com_lengths(5),status)
   
         if(status .ne. 0) then
            write(*,*) "THERE WAS A PROBLEM READING THE TIMESTAMP"
            stop
         end if
         
         read(com_line(5),'(F10.0)') time_stamp

      end if

      return
    
    end subroutine get_args
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine print_help
      implicit none

      write(*,*) '------------------- bx3d command line arguements ------------------------'
      write(*,*)
      write(*,*) ' bx3d <opt> <mesh_file_prefix> <input> <output_filename>'
      write(*,*)
      write(*,*) '  opt ='
      write(*,*) '       -l or -list for a list of new files to enter into new exodus file'
      write(*,*) '       -al or -alist for a list of new files to add into the existing input exodus file'
      write(*,*) '       -f or -file to build a single file into a new exodus file'
      write(*,*) '       -af or -afile to add a single file to an existing exodus file'
      write(*,*) 
      write(*,*) '  mesh_file_prefix = the prefix of the mesh files'
      write(*,*) '  input = the list for file name to inserted into the output file'
      write(*,*) '  output_filename = the name of the output exodus file'
     
    end subroutine print_help
    !_______________________________________________________________
    !_______________________________________________________________
    subroutine check_output(status)
      !checks to make sure output file exists if the add option is 
      !chosen
      implicit none
      integer :: status
      logical :: exst

      status = 0
      if(trim(com_line(1)) == '-al' .or. trim(com_line(1))=='-alist' &
         .or. trim(com_line(1))== '-af' .or. trim(com_line(1))=='-afile') then
         
         inquire(file=trim(com_line(4)),exist=exst)
         if(.not. exst) then
            write(*,*) "CANNOT FIND THE EXISTING OUTPUT FILE: ",trim(com_line(4))
            write(*,*) "BUILDING NEW EXODUS FILE"
            status=0
         else
            newexo = .false.
         end if
      else
         write(*,*) "BUILDING NEW EXODUS FILE"
      end if
      return
    end subroutine check_output
    !_______________________________________________________________
    !_______________________________________________________________
    subroutine check_input(status)
      !checks the input file or input file list
      implicit none
      character*80 :: test_string, cfil
      integer :: status
      logical :: exst
      integer :: nn,nele,ntest,iostatus,flen,ppos,nfil
      real :: tstamp

      status = 0
      !determine if the input exists
      inquire(file=trim(com_line(3)),exist=exst)
      if(.not. exst) then
         write(*,*) 'CANNOT FIND THE INPUT FILE: ',trim(com_line(3))
         status = -1
         return
      end if
      
      !get the number of nodes and elements
      open(25,file=trim(com_line(2))//'.node',status='old',action='read')
      read(25,*) nn
      close(25)

      open(25,file=trim(com_line(2))//'.ele',status='old',action='read')
      read(25,*) nele
      close(25)
      
      ppos=index(trim(com_line(3)),'.',.true.)
      flen=len_trim(com_line(3))
      
      !if this is a file, make sure the number of values are correct
      if(trim(com_line(1)) == '-f' .or. trim(com_line(1))=='-file' &
         .or. trim(com_line(1))== '-af' .or. trim(com_line(1))=='-afile') then
         open(25,file=trim(com_line(3)),status='old',action='read')
         read(25,*,IOSTAT=iostatus) ntest
         close(25)
         if(iostatus .ne. 0) goto 100  
         if(ntest .eq. nn .or. ntest .eq. nele) then
            status = 0
            return
         else
            write(*,1001) trim(com_line(3))
            write(*,*) ntest,nn,nele
            status = -1
            return
         end if
       
      end if
     
      !if this is a list, make sure each file in the list exists and
      !had the correct number of values 
      if(trim(com_line(1)) == '-l' .or. trim(com_line(1))=='-list' &
           .or. trim(com_line(1))== '-al' .or. trim(com_line(1))=='-al') then
         
         open(25,file=trim(com_line(3)),status='old',action='read')
         read(25,*,IOSTAT=iostatus) nfil
         if(iostatus .ne. 0) goto 101
       
         do i=1,nfil
            read(25,*,IOSTAT=iostatus) cfil, tstamp
            if(iostatus .ne. 0) goto 102
            
            inquire(file=trim(cfil),exist=exst)
            if(.not.exst) goto 103
            open(26,file=trim(cfil),status='old',action='read')
            read(26,*,IOSTAT=iostatus) ntest
            if(iostatus .ne. 0) goto 104
            
            if(ntest .ne. nn .and. ntest .ne. nele) then
               write(*,*) i,ntest,nn,nele
               write(*,1001) trim(cfil)
               status=-1
               return
            end if
            close(26)

         end do
         close(25)
      end if

      status = 0
      return

100   if(iostatus>0) then
         write(*,*) 'THERE WAS A PROBLEM READING THE NUMBER OF VALUES IN: ',trim(com_line(3))
         status=-1
         return
      else
          write(*,*) 'REACHED END OF FILE READING THE NUMBER OF VALUES IN: ',trim(com_line(3))
         status=-1
         return
      end if

101   if(iostatus>0) then
         write(*,*) 'THERE WAS A PROBLEM READING THE NUMBER OF FILES IN: ',trim(com_line(3))
         status=-1
         return
      else
          write(*,*) 'REACHED END OF FILE READING THE NUMBER OF FILES IN: ',trim(com_line(3))
         status=-1
         return
      end if

102   if(iostatus>0) then
         write(*,*) 'THERE WAS A PROBLEM READING THE FILE NAME AND TIME STAMP IN',trim(com_line(3))
         write(*,*) 'FOR FILE NUMBER: ',i
         status=-1
         return
      end if
      
103   write(*,*) 'COULD NOT FIND: ',trim(cfil),' LISTED IN: ',trim(com_line(3)),' ON LINE: ',i+1
      status = -1
      return

104   if(iostatus>0) then
         write(*,*) 'THERE WAS A PROBLEM READING THE NUMBER OF VALUES IN: ',trim(cfil)
         status=-1
         return
      else
          write(*,*) 'REACHED END OF FILE READING THE NUMBER OF VALUES IN: ',trim(cfil)
         status=-1
         return
      end if


1001 FORMAT("THE NUMBER OF VALUES LISTED IN: ",A," DOES NOT MATCH THE NUMBER OF NODES OR ELEMENTS");
    end subroutine check_input
    !_______________________________________________________________    
    !_______________________________________________________________
    subroutine check_mesh(status)
      !checks to see if mesh files exist
      implicit none
      integer :: status
      logical :: exst
      integer :: i
      character*80 :: fname
      status=-1;

      write(tfname,"(A80)") com_line(2)
      nchar=90
      do i=1,90
         if(tfname(i:i)== '.') then
            nchar = i-1
            exit
         end if
      end do
    
      do i=1,5
         select case(i)
            case(1)
               inquire(file=trim(com_line(2))//'.node',exist=exst)
               if(.not.exst)  write(fname,*) trim(com_line(2))//'.node'
            case(2)
               inquire(file=trim(com_line(2))//'.ele',exist=exst)
               if(.not.exst)  write(fname,*) trim(com_line(2))//'.ele'
            case(3)
               inquire(file=trim(com_line(2))//'.face',exist=exst)
               if(.not.exst) write(fname,*) trim(com_line(2))//'.face'
            case(4)
               inquire(file=trim(com_line(2))//'.neigh',exist=exst)
               if(.not.exst) write(fname,*) trim(com_line(2))//'.neigh'
            case(5)      
               inquire(file=tfname(1:nchar)//'.trn',exist=exst)
               if(.not.exst) write(fname,*) trim(tfname)//'.trn'
            end select
         
            if(.not.exst) then      
               write(*,1001) trim(fname) 
               return
            end if
      end do

      status=0
      return
      
1001 FORMAT("CANNOT FIND THE MESH FILE:",A)
    end subroutine check_mesh
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine check_option(status)
      !checks to make sure the command line option is valid
      implicit none
      integer :: status
      status = 0
      select case(trim(com_line(1)))
         case('-l')
            return
         case('-list')
            return
         case('-al')
            return
         case('-alist')
            return   
         case('-f')
            return
         case('-file')
            return
         case('-af')
            return
         case('-afile')
            return
         case default
            status = -1
            return
         end select
    end subroutine check_option
    !_______________________________________________________________
    
    !_______________________________________________________________
    subroutine load_tetgen
      !LOAD THE MESH FILES
      implicit none
      integer :: j1,j2,i,np,nat

      !!allocate and read the node file
      open(10,file=trim(com_line(2))//".node",status='old',action='read')
      read(10,*) nnods,dim 
      allocate(nods(nnods,dim),n_type(nnods,2))
      write(*,*) "Reading ",trim(com_line(2))//".node"
      do i=1,nnods
         read(10,*) j1,nods(i,1:3),n_type(i,1:2)
      end do

      !nnsets=maxval(n_type)
      nnsets=1
      close(10)
     
      !!translate the nodes
      open(10,file=tfname(1:nchar)//'.trn',status='old')
      read(10,*) tx,ty,tz
      close(10)
      nods(:,1)=nods(:,1)+tx
      nods(:,2)=nods(:,2)+ty
      nods(:,3)=nods(:,3)+tz

      
      
      !!allocate and read the elements
      open(10,file=trim(com_line(2))//".ele",status='old',action='read')
      read(10,*) nele,np,nat
      allocate(ele(nele,np+nat))
      write(*,*) "Reading ",trim(com_line(2))//".ele"
      do i=1,nele
         read(10,*) j1,ele(i,1:np+nat)
      end do
      close(10)
      nzones=maxval(ele(:,np+nat))-minval(ele(:,np+nat))+1

      !!allocate and read the neighbors file
      open(10,file=trim(com_line(2))//".neigh",status='old',action='read')
      read(10,*) nneigh,np
      allocate(neigh(nneigh,np))
      write(*,*) "Reading ",trim(com_line(2))//".neigh"
      do i=1,nneigh
         read(10,*) j1,neigh(i,1:4)
      end do
      close(10)

      !allocate and read the faces file
      open(10,file=trim(com_line(2))//".face",status='old',action='read')
      read(10,*) nface
      allocate(face(nface,4))
      write(*,*) "Reading ",trim(com_line(2))//".face"
      do i=1,nface
         read(10,*) face(i,1:4)
      end do
      close(10)

     call get_nnode_sets
     
     write(*,*) "THE NUMBER OF NODES IS:         ",nnods
     write(*,*) "THE NUMBER OF NODE SETS IS:     ",nnsets
     write(*,*) "THE NUMBER OF ELEMENTS IS:      ",nele
     write(*,*) "THE NUMBER OF ELEMENT ZONES IS: ",nzones
    end subroutine load_tetgen
    !_______________________________________________________________
    
    !_______________________________________________________________
    subroutine get_nnode_sets
      !counts the number of node sets
      implicit none
      integer :: i,j
      integer, dimension(10000) :: nsets
      nsets(1)=n_type(1,2)
      nnsets=1
      do i=2,nnods
         do j=1,nnsets
            if(nsets(j)==n_type(i,2)) goto 100
         end do
        nnsets=nnsets+1
        nsets(nnsets)=n_type(i,2)
100     continue
      end do
      
      allocate(nset_ids(nnsets),nset_sizes(nnsets))
      nset_sizes = 0
      nset_ids(1:nnsets) = nsets(1:nnsets)
      do i=1,nnods
         do j=1,nnsets
            if(n_type(i,2) == nset_ids(j)) then
               nset_sizes(j)=nset_sizes(j)+1
               goto 200
            end if
         end do
200      continue
      end do
      
      !do i=1,nnsets
      !   write(*,*) i,nset_ids(i),nset_sizes(i)
      !end do
    end subroutine get_nnode_sets
    !_______________________________________________________________

    !_______________________________________________________________
    subroutine get_elements
      implicit none
      integer :: i,j1,tnele,np,nat

      !!allocate and read the elements
      open(10,file=trim(com_line(2))//".ele",status='old',action='read')
      read(10,*) tnele,np,nat
      if(tnele .ne. nele) then
         write(*,*) 'THE NUMBER OF ELEMENTS IN ',trim(com_line(2))//".ele",tnele
         write(*,*) 'IS NOT EQUAL TO THE NUMBER OF ELEMENTS IN ',trim(com_line(4)),nele
         stop
      end if
      allocate(ele(tnele,np+nat))
      write(*,*) "Reading ",trim(com_line(2))//".ele"
      do i=1,nele
         read(10,*) j1,ele(i,1:np+nat)
      end do
      close(10)
      nzones=maxval(ele(:,np+nat))-minval(ele(:,np+nat))+1
    end subroutine get_elements
    !_______________________________________________________________
    

    !_______________________________________________________________
    subroutine get_nodes
      implicit none
      integer :: i,j1, tnnods,ndm,nat1,nat2
      open(10,file=trim(com_line(2))//".node",status='old',action='read')

      read(10,*) tnnods,ndm,nat1,nat2
      if(tnnods .ne.nnods) then
         write(*,*) 'THE NUMBER OF NODES IN ',trim(com_line(2))//".node",tnnods
         write(*,*) 'IS NOT EQUAL TO THE NUMBER OF NODES IN ',trim(com_line(4)),nnods
         stop
      end if
      
      allocate(nods(nnods,dim),n_type(nnods,2))
      write(*,*) "Reading ",trim(com_line(2))//".node"
      do i=1,nnods
         read(10,*) j1,nods(i,1:3),n_type(i,1:2)
      end do

      !nnsets=maxval(n_type)
      nnsets=1
      close(10)
     
      !!translate the nodes
      open(10,file=tfname(1:nchar)//'.trn',status='old')
      read(10,*) tx,ty,tz
      close(10)
      nods(:,1)=nods(:,1)+tx
      nods(:,2)=nods(:,2)+ty
      nods(:,3)=nods(:,3)+tz

      call get_nnode_sets

    end subroutine get_nodes
    !_______________________________________________________________
end module subroutines
