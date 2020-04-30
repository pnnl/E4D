program main
  
  use subroutines
  use exo_subs
  implicit none

  
  !!get and check the arguments from the command line
  call get_args
  
  
  if(newexo) then
     !! build a new exodus file
     call load_tetgen
     call init_exo
     call record_mesh

     !!record the files specied on the command line
     call enter_data

  else
     !!open the old exodus file
     call open_exo
     call get_elements
     call enter_data
  end if
  
  

  !!open and initialize the exodos file
  !call init_exo

  !!record the mesh nodes and elements
  !call record_mesh

  !!record the conductivity values
  !call record_sigma

  !!close the file
  call close_exo
  
  
!!$  character*80 :: fname="testfile.exo"
!!$  
!!$
!!$  integer :: icompws,iows,ierr,idex
!!$  
!!$  icompws = 4
!!$  iows = 4
!!$
!!$  idex = excre(fname,EXCLOB,icompws,iows,ierr)
!!$  write(*,*) idex,ierr
!!$  
!!$  call exclos(idex,ierr)
!!$  write(*,*) ierr


end program main
