module assemble_fmm

!!____________________________________________________________________________________________  
!!
!! Description: The assemble module contains subroutines for collecting
!! the travel time from each slave processors that are necessary to   
!! to reconstruct a simulated FMM survey 
!!
!!____________________________________________________________________

  use fmm_vars
  use vars
  implicit none

contains

  !___________________________________________________________________
  subroutine fmm_assemble_data

    implicit none
    integer :: i                      !indexing variable
    integer :: a                    !current source indexes
    integer :: m                    !receiver index
    integer :: s1,s2                  !first and last source stored by this process
    integer :: indx                   !current index in my_drows and my_dvals
    integer :: p                      !solution index
    
   
    !! set the pole solution ownership bounds for this process
    s1=sind(my_rank_fmm,1)                
    s2=sind(my_rank_fmm,2) 
    
    !! find the number of measurements using a source owned by this process
    !! and allocate my_drows, my_dvals
!    if(allocated(my_drows).and.res_flag) then
!       deallocate(my_drows)
!       deallocate(my_dvals)
!    end if
    if(.not.allocated(my_drows)) then
       nmy_drows=0
       do i=1,nm_fmm
          a=s_conf_fmm(i,1)
          if( a>=s1 .and. a<=s2 ) then
             nmy_drows=nmy_drows+1
          end if
       end do
       allocate(my_drows(nmy_drows),my_dvals(nmy_drows))
    end if
    
    indx=0
    my_dvals=0
    
   
       !! Assemble the part of the simulated data owned by this process
       do i=1,nm_fmm
          !! Set a,b,m and n for this measurement
          a = s_conf_fmm(i,1)
          m = s_conf_fmm(i,2)
         
          if(a>=s1 .and. a<=s2) then
             !! If this process owns a source for this
             !! measurement, then compute the part of the measurement depending on
             !! on this source
             indx=indx+1
             my_drows(indx)=i
             
             do p=s1,s2
                !! Loop of the sources owned by this process, find the ones used by this
                !! measurement, and assign the simulated travel time at the
                !! receiver (m)
                if(p==a) then
                   if(fresnel) then
                      if(m.ne.0) my_dvals(indx)= real(ttimes(s_nods(m),a-s1+1))+.25/frq(a)
                   else
                      if(m.ne.0) my_dvals(indx)= real(ttimes(r_nods(m),a-s1+1))
                   end if
                end if

             end do
          end if
       end do
    
  end subroutine fmm_assemble_data
  !___________________________________________________________

end module assemble_fmm
