module build_amap
!!____________________________________________________________________________________________  
!!  E4D
!!  Copyright © 2014, Battelle Memorial Institute
!!  All rights reserved.
!!
!!  1 . Battelle Memorial Institute (hereinafter Battelle) hereby grants permission to any
!!  person or entity lawfully obtaining a copy of this software and associated documentation 
!!  files (hereinafter “the Software”) to redistribute and use the Software in source and 
!!  binary forms, with or without modification.  Such person or entity may use, copy, modify, 
!!  merge, publish, distribute, sublicense, and/or sell copies of the Software, and may permit 
!!  others to do so, subject to the following conditions:
!! 
!!   -Redistributions of source code must retain the above copyright notice, this list of 
!!    conditions and the following disclaimers. 
!!   -Redistributions in binary form must reproduce the above copyright notice, this list 
!!    of conditions and the following disclaimer in the documentation and/or other materials 
!!    provided with the distribution. 
!!   -Other than as used herein, neither the name Battelle Memorial Institute or Battelle may
!!    be used in any form whatsoever without the express written consent of Battelle.  
!!   -Redistributions of the software in any form, and publications based on work performed 
!!    using the software should include the following citation as a reference:
!!    [Cite published manuscript.  If this does not apply, delete this bulleted item].
!!
!!  2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
!!  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
!!  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
!!  BATTELLE OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
!!  OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!!  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!!  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
!!  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
!!  POSSIBILITY OF SUCH DAMAGE.
!!
!!  DISCLAIMER
!!  The Software was produced by Battelle under Contract No. DE-AC05-76RL01830 with the 
!!  Department of Energy.  The U.S. Government is granted for itself and others acting on its 
!!  behalf a nonexclusive, paid-up, irrevocable worldwide license in this data to reproduce, 
!!  prepare derivative works, distribute copies to the public, perform publicly and display 
!!  publicly, and to permit others to do so.  The specific term of the license can be identified
!!  by inquiry made to Battelle or DOE.  Neither the United States nor the United States 
!!  Department of Energy, nor any of their employees, makes any warranty, express or implied, 
!!  or assumes any legal liability or responsibility for the accuracy, completeness or 
!!  usefulness of any data, apparatus, product or process disclosed, or represents that its use 
!!  would not infringe privately owned rights.
!!______________________________________________________________________________________________
!!______________________________________________________________________________________________
!! Author: Tim Johnson
!!
!! email: e4d@pnnl.gov
!!
!! Purpose: The build_amap module builds the mapping vectors and scaling
!!          factors necessary to efficiently reconstruct the forward 
!!          coupling and inverse Jacobian matrices. 
!!
!! Input variables: see individual subroutines
!!
!! Output variables: see individual subroutines
!!____________________________________________________________________

  use vars

contains

  !__________________________________________________________________________
  subroutine build_delA
 
    !________________________________________________________________________
    !! Author: Tim Johnson
    !!
    !! email: tj@pnnl.gov
    !!
    !! Purpose: The build_delA subroutine builds the mapping vectors and scaling
    !!          factors necessary to efficiently reconstruct the forward 
    !!          coupling matrix and the Jacobian matrix. 
    !! USE variables accessed or modified:
    !! vars module
    !! rows  : holds the row index of each element in delA
    !! cols  : holds the column index of each element in delA
    !! delA  : holds the forward matrix coefficient scaling factors
    !! nvals : number of scaling factors in delA
    !! A_map : forward matrix reconstruction mapping vector
    !! nelem : number of mesh elements
    !! elements : holds the node indexes for each element
    !! S_map    : conductivity mapping vector
    !! d_nnz : petsc preallocation vector
    !_________________________________________________________________________
   
    implicit none
    integer :: nnds                              !number of nodes in mesh
    integer :: i,j,k,l                           !counters
    integer :: mcount                            !current index in mapping vectors
    integer :: rn,cn                             !row and column numbers
    integer :: ncl                               !number of non-zero columns on this row
  
    real*8, dimension(4,4) :: atild              !pre-shape linear shape function matrix 
    real*8, dimension(4,4) :: atildi             !shape function matrix
    real*8, dimension(4,4) :: temp_a             !temporary placeholder

    integer, dimension(4) :: indx                !computation space for MIGS routine

    integer, dimension(:,:),allocatable :: A_tmp !temporary accounting matrix

    real*8 :: x1,x2,x3,x4                          !node x-positions for the current element
    real*8 :: y1,y2,y3,y4                          !nodal y-positions for the current element
    real*8 :: z1,z2,z3,z4                          !nodal z-positions for the current element
    real*8 :: evol                                 !nodal volume

    real*8 :: A_kij                                !forward matrix coefficient sum

    integer, dimension(:), allocatable :: t1vec     !temporary placeholder
    integer, dimension(:), allocatable :: t2vec     !temporary placeholder
  
 
   
    !allocate the petsc matrix preallocation vectors
    allocate(d_nnz(nnodes))
    d_nnz=0
    

    !!Allocate the mapping vectors. There are at most 10 contributions 
    !!to the forward coupling matrix per element
    allocate(rows(10*nelem),cols(10*nelem))
    allocate(A_map(10*nelem),S_map(10*nelem))
    allocate(delA(10*nelem))
    
 
    !!Count the number of node pairs (i.e. coupling matrix elements)
    !!and reallocate the mapping and coefficient vectors
    rows(1)=elements(1,1)
    cols(1)=elements(1,1)
    nvals=1
    mcount = 0
    A_map = 0
    nnds=maxval(elements(:,:))
    
    !!Allocate a temporary matrix to aid in building the 
    !!mapping vectors
    allocate(A_tmp(nnds,80))
    A_tmp=0
    A_tmp(rows(1),1) = 1
    A_tmp(rows(1),2) = nvals
    A_tmp(rows(1),3) = cols(1)
 
    !!loop over the elements   
    do k=1,nelem    
       
       !!loop over each node in this element
       do i=1,4
         
          !!rn is a row index of the coupling matrix
          rn=elements(k,i)    
          
          !!loop over each node in this element
          do j=1,4
             
             !!cn is a column index of the coupling matrix
             cn=elements(k,j)
        
          

             !!check to see if this pair (rn,cn) is in the lower
             !!triangle. If not, go to the next pair
             if (cn <= rn) then
                mcount=mcount+1

                !!loop over each pair found thus far (nvals pairs)
                !!and determine if the pair (rn,cn) is already 
                !!represented
                if(A_tmp(rn,1) == 0) then
                   nvals = nvals + 1
                   A_tmp(rn,1) = 1
                   A_tmp(rn,2) = nvals
                   A_tmp(rn,3) = cn
                   rows(nvals) = rn
                   cols(nvals) = cn
                   A_map(mcount) = nvals
                           
                else
                   
                   do l=1,A_tmp(rn,1)
                      if(cn == A_tmp(rn,2*l+1)) then
                         A_map(mcount) = A_tmp(rn,2*l)
                         goto 11
                      end if
                   end do
                   
                   !!if were here then no column index was found
                   !!for this row, so we'll add one
                   nvals = nvals+1
                   ncl = A_tmp(rn,1) + 1
                   A_tmp(rn,1) = ncl
                   A_tmp(rn,2*ncl) = nvals
                   A_tmp(rn,2*ncl+1) = cn
                   rows(nvals) = rn
                   cols(nvals) = cn
                   A_map(mcount) = nvals
                
11                 continue
     
                end if

             end if

          end do
          
       end do
       
    end do
    
    !build the petsc allocation vector
    do i=1,nnodes
       !upper part
       d_nnz(i)=A_tmp(i,1)-1;
       do j=1,A_tmp(i,1)
          cn=A_tmp(i,2*j+1)
          d_nnz(cn)=d_nnz(cn)+1
       end do
    end do   
    deallocate(A_tmp)
   
    !!now reallocate for the correct number of matrix values
    allocate(t1vec(nvals),t2vec(nvals))
    t1vec(1:nvals)=rows(1:nvals)
    t2vec(1:nvals)=cols(1:nvals)
    deallocate(rows,cols)
    allocate(rows(nvals),cols(nvals))
    rows=t1vec
    cols=t2vec
    deallocate(t1vec,t2vec)
    allocate(trows(nvals),tcols(nvals))
    
    !!loop over the elements and build delA
    mcount = 0
    delA=0
  
    do k=1,nelem
   
       !!get the 4 nodes for this element 
       x1 = dble(nodes(elements(k,1),1))
       y1 = dble(nodes(elements(k,1),2))
       z1 = dble(nodes(elements(k,1),3))
       
       x2 = dble(nodes(elements(k,2),1))
       y2 = dble(nodes(elements(k,2),2))
       z2 = dble(nodes(elements(k,2),3))
       
       x3 = dble(nodes(elements(k,3),1))
       y3 = dble(nodes(elements(k,3),2))
       z3 = dble(nodes(elements(k,3),3)) 
       
       x4 = dble(nodes(elements(k,4),1))
       y4 = dble(nodes(elements(k,4),2))
       z4 = dble(nodes(elements(k,4),3))
       
       !!get the volume of this element
       call get_vol (x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,evol)  
   
       if(evol==0) goto 12
       !!build the linear shape function matrix for this element
       atild(1,:) = (/1,1,1,1/)
       atild(2,:) = (/x1,x2,x3,x4/)
       atild(3,:) = (/y1,y2,y3,y4/)
       atild(4,:) = (/z1,z2,z3,z4/)
       temp_a=atild
            
       !!invert temp_a to get atildi
       !!note MIGS routine is located in mat_inv.f  
       call MIGS(temp_a,4,atildi,indx)
       
      
12     continue
       !!use atildi to build the coupling coefficients
       !!for this element      
       
       !!loop over each node in this element
       do i=1,4
          !!rn is a row index of the coupling matrix
          rn=elements(k,i)
          
          !!loop over each node in this element
          do j=1,4
             
             !!cn is a column index of the coupling matrix
             cn=elements(k,j)
             
             !!check to see if this pair (rn,cn) is in the lower
             !!triangle. If not, go to the next pair
             if (cn <= rn) then
                mcount=mcount+1
                !!Compute the coupling coefficient for this element
                !!and node pair
                A_kij=0
               
                do l=2,4
                   A_kij=A_kij+(atildi(i,l)*atildi(j,l))
                end do
          
                delA(mcount)=(A_kij*evol)
                S_map(mcount)=k
               
             end if
             
21           continue
          end do
          
       end do
       
    end do  
  

    
  end subroutine build_delA
  !___________________________________________________________________________
  
 
  
  !___________________________________________________________________________
  subroutine get_vol( x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, volume )

    !!computes the volume of a tetrahedron

    implicit none
    
    real*8 a(4,4)
    real*8 r4_det
    real*8 volume
    real*8 x1
    real*8 x2
    real*8 x3
    real*8 x4
    real*8 y1
    real*8 y2
    real*8 y3
    real*8 y4
    real*8 z1
    real*8 z2
    real*8 z3
    real*8 z4
    
    a(1:4,1) = (/ x1, x2, x3, x4 /)
    a(1:4,2) = (/ y1, y2, y3, y4 /)
    a(1:4,3) = (/ z1, z2, z3, z4 /)
    a(1:4,4) = (/ 1.0E+00, 1.0E+00, 1.0E+00, 1.0E+00 /)
    r4_det = det4(a)
    volume = abs ( r4_det ) / 6.0E+00
    
    return
  end subroutine get_vol
  !___________________________________________________________________________
    
  !___________________________________________________________________________
  function det4 ( a1 )

    !computes the determinant of a 4x4 matrix
    implicit none    
    real*8 :: a1(4,4)
    real*8 :: a(4,4)
    real*8 :: det4
    integer ::i
    real*8 :: c1,c2,c3,c4
    
    a = dble(a1)
    c1=a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) 
    
    c2=a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))
    
    c3=a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    c4=a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*&
         &(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
    
    det4 = real(c1 - c2 + c3 - c4)
    
    return
  end function det4
  !__________________________________________________________________________
  
 
end module build_amap
