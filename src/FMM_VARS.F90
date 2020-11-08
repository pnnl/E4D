module fmm_vars

implicit none
!#include "include/finclude/petscsys.h"
!#include "include/finclude/petscvec.h"
!#include "include/finclude/petscvec.h90"
!#include "include/finclude/petscmat.h"
!#include "include/finclude/petscmat.h90"
!#include "include/finclude/petscviewer.h"
!#include "include/finclude/petscviewer.h90"  
!#include "include/finclude/petscksp.h"
!#include "include/finclude/petscksp.h90"
  real, parameter :: tt_max = 1e12
  logical :: simulate_fmm = .false.
  logical :: fresnel = .false.                                !!flag for fresnel volume sensitivities
  logical :: fresnel_out = .false.                            !!flage to output fresnel volumes
  logical, dimension(:), allocatable :: clos
  logical, dimension(:), allocatable :: upstream
  logical, dimension(:), allocatable :: use_ele
  
  
  !FILES
  character*40 :: mshfile_fmm                                 !!file containing the mesh options
  character*40 :: tfile                                       !!survey configuration file
  character*40 :: spdfile                                     !!slowness file
  character*40 :: outfile_fmm                                 !!output options file
  
  !Integers
  integer :: num_fmm_procs                                    !number of processes for fmm simulation
  integer :: FMM_COMM                                         !FMM communicator
  integer :: SCOMM_FMM                                        !!communicator for slaves only
  integer :: heap_end
  integer :: num_heap_done
  integer :: nzf                                              !!number of zones used in travel time simulation

  integer :: master_proc_fmm                                  !root of fmm processes
  integer :: my_rank_fmm                                      !!my mpi rank
  integer :: n_rank_fmm                                       !!number of processes
  integer :: my_ns                                            !!number of soures I'm assigned
  real :: my_frt_fmm                                          !!my forward run time
  real, dimension(:,:), allocatable :: s_pos                  !!source positions
  real, dimension(:,:), allocatable :: rc_pos                 !!receiver positions
  real, dimension(:), allocatable :: dobs_fmm,Wd_fmm          !!observed data and weighting vector
  real, dimension(:), allocatable :: velocity                    !!element velocity
  real, dimension(:), allocatable :: TT
  real, dimension(:,:), allocatable :: ttimes
  real, dimension(:), allocatable :: ttpred                   !!predicted data vector
  real, dimension(:), allocatable :: frq                      !!wave frequencies for fresnel sensitivities

  integer :: mode_fmm                                         !!run mode
  integer :: ns,tns                                           !!num my sources, total num sources
  integer :: nrc                                              !!num of receivers
  integer :: nspd                                             !!number of elements, velocity values
  integer :: nm_fmm                                           !!number of measurements
  integer, dimension(:,:), allocatable :: s_conf_fmm          !!survey configuration
  integer, dimension(:), allocatable :: s_nods                !!indices of source nodes
  integer, dimension(:), allocatable :: r_nods                !!indices of receiver nodes
  integer, dimension(:,:), allocatable :: sind                !!source assignments 
  integer, dimension(:,:),allocatable :: tind
  integer, dimension(:), allocatable :: ring_map      
  integer, dimension(:), allocatable :: rings
  integer, dimension(:), allocatable :: heap
  integer, dimension(:), allocatable :: heap_map
  integer, dimension(:), allocatable :: zsims                 !!zone numbers to include in travel time sims
  
  real, dimension(3,3) :: M
  real, dimension(4,3) :: XP
  real, dimension(4) :: LTT
  real, dimension(1,1) :: tt_test
  
end module fmm_vars
