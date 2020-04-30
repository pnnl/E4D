module reorder_mesh
  !this module includes subroutines to reorder the mesh in order
  !to remove inactive elements (zone flag of -999)

  use vars
  use input
  implicit none
  
 
contains
!____________________________________________________________________
  subroutine remove_inactive_elements
    !reorders the mesh to remove inactive elements and creates the 
    !corresponding mapping files between meshes

    implicit none
    integer :: i,j,nnew_nodes,nnew_elements,cpos,nnew_faces
    integer, dimension(:,:), allocatable :: temp_elements
    real, dimension(:,:), allocatable :: temp_nodes
    logical, dimension(:), allocatable :: inc_nodes, inc_elements, inc_faces
    
    write(*,*) 
    write(*,*) "RE-ORDING MESH TO REMOVE INACTIVE ELEMENTS"

    
    !!allocate and populate a logical array indicating which nodes are used
    allocate(inc_nodes(nnodes))
    inc_nodes = .false.
    do i=1,nelem
       if(zones(i) .ne. -999) then
          do j=1,4
             inc_nodes(elements(i,j))=.true.
          end do
       end if  
    end do
    
    !!construct the map from old node indices to new ones
    allocate(node_map(nnodes))
    node_map = 0
    cpos=0
    do i=1,nnodes
       if(inc_nodes(i)) then
          cpos = cpos+1
          node_map(i) = cpos
       end if
    end do
    nnew_nodes = cpos
 
    !!map the elements
    allocate(element_map(nelem))
    element_map = 0
    cpos = 0
    do i=1,nelem
       if(zones(i) .ne. -999) then
          cpos = cpos+1
          element_map(i) = cpos
          do j=1,4
             elements(cpos,j) = node_map(elements(i,j))
          end do
       end if
    end do
    nnew_elements = cpos

    !!rebuild the nodes and node flags
    allocate(temp_nodes(nnew_nodes,3),temp_elements(nnew_nodes,1))
    do i=1,nnodes
       if(node_map(i) .ne. 0) then
          temp_nodes(node_map(i),:) = nodes(i,:)
          temp_elements(node_map(i),1) = nbounds(i)
       end if
    end do
    deallocate(nodes,nbounds)
    allocate(nodes(nnew_nodes,3),nbounds(nnew_nodes))
    nodes = temp_nodes
    nbounds = temp_elements(:,1)
    nnodes = nnew_nodes
    deallocate(temp_nodes,temp_elements)

    !!rebuild the elements
    allocate(temp_elements(nnew_elements,4))
    temp_elements = elements(1:nnew_elements,:)
    deallocate(elements)
    allocate(elements(nnew_elements,4))
    elements = temp_elements  

    temp_elements=0

    !!rebuild the zones
    do i=1,nelem
       if(element_map(i) .ne. 0) then
          temp_elements(element_map(i),1) = zones(i)
       end if
    end do
    deallocate(zones)
    allocate(zones(nnew_elements))
    zones = temp_elements(:,1)
    

    !!rebuild the faces
    allocate(inc_faces(nfaces))
    inc_faces = .false.
    nnew_faces = 0
    do i=1,nfaces
       if(inc_nodes(faces(i,1)) .and. inc_nodes(faces(i,2)) .and. inc_nodes(faces(i,3))) then
          nnew_faces = nnew_faces+1
          temp_elements(nnew_faces,4) = faces(i,4)
          do j=1,3
             temp_elements(nnew_faces,j) = faces(node_map(i),j)
          end do
       end if
    end do
    deallocate(faces)
    allocate(faces(nnew_faces,4))
    faces=temp_elements(1:nnew_faces,:)
    nfaces = nnew_faces
    deallocate(temp_elements)

    !rebuild the conductivity file
    if(allocated(sigma)) then
       allocate(temp_nodes(nnew_elements,1))
       do i=1,nelem
          if(element_map(i) .ne. 0) then
             temp_nodes(element_map(i),1) = sigma(i)
          end if
       end do
       deallocate(sigma)
       allocate(sigma(nnew_elements))
       sigma = temp_nodes(1:nnew_elements,1)
       deallocate(temp_nodes)
    end if

    if(allocated(sigmai)) then
       allocate(temp_nodes(nnew_elements,1))
       do i=1,nelem
          if(element_map(i) .ne. 0) then
             temp_nodes(element_map(i),1) = sigmai(i)
          end if
       end do
       deallocate(sigmai)
       allocate(sigmai(nnew_elements))
       sigmai = temp_nodes(1:nnew_elements,1)
       deallocate(temp_nodes)
    end if
    nelem = nnew_elements

!!$    open(10,file='temp.1.node',status='replace',action='write')
!!$    write(10,*) nnodes,3,1,1
!!$    do i=1,nnodes
!!$       write(10,*) i,nodes(i,:),1,nbounds(i)
!!$    end do
!!$    close(10)
!!$
!!$    open(10,file='temp.1.ele',status='replace',action='write')
!!$    write(10,*) nelem,4,1
!!$    do i=1,nelem
!!$       write(10,*) i,elements(i,:),zones(i)
!!$    end do
!!$    close(10)
!!$
!!$    open(10,file='temp.sig',status='replace',action='write')
!!$    write(10,*) nelem,1
!!$    do i=1,nelem
!!$       write(10,*) sigma(i)
!!$    end do
!!$    close(10)
!!$
!!$    open(10,file='temp.1.face',status='replace',action='write')
!!$    write(10,*) nfaces,1
!!$    do i=1,nfaces
!!$       write(10,*) i,faces(i,:)
!!$    end do
!!$    close(10)


    
    
  end subroutine remove_inactive_elements
!____________________________________________________________________

end module reorder_mesh
