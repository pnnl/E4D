module vars

#ifdef petsc_7
implicit none
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscmat.h90"
#include "petsc/finclude/petscviewer.h"
#include "petsc/finclude/petscviewer.h90"  
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscksp.h90"
#else
#include "petsc/finclude/petscksp.h"
use petscksp
implicit none
#endif


  !FILES
  character*40 :: mshfile                                  !!file containing the mesh options
  character*40 :: efile                                    !!survey configuration file
  character*40 :: sigfile                                  !!bulk conductivity file
  character*40 :: outfile                                  !!output options file
  character*40 :: invfile                                  !!inversion options file
  character*40 :: iinvfile                                  !!inversion options file for phase components
  character*10 :: com_inp                                  !!command line input

  character*40 :: tmpstr                                   !!shareable string
  
  !Integers
  integer :: my_rank                                       !!my mpi rank
  integer :: ierr                                          !!generall error
  integer :: n_rank                                        !!number of processes
  integer :: secflag                                       !!singularity removal flag
  integer :: frun_flag                                     !!solver flag
  integer :: mode                                          !!run mode
  integer :: ne,tne                                        !!num my electrodes, total num electrodes
  integer :: nsig                                          !!number of elements, sigma values
  integer :: nm                                            !!number of measurements
  integer :: msh_opt                                       !!mesh option flag
  integer :: nrz                                           !!number of regularization zones,total #zones
  integer :: conv_opt                                      !!convergence option
  integer :: min_initer,max_initer                         !!min/max inner iterations
  integer :: bound_flag                                    !!boundary options flag
  integer :: nnodes                                        !!number of nodes
  integer :: nelem                                         !!number of elements (same as nsig)
  integer :: nfaces                                        !!number of faces
  integer :: nvals                                         !!number of non-zeros in coupling mat
  integer :: n_srank                                       !!number slaves in communicator SCOMM
  integer :: SCOMM                                         !!communicator for slaves only
  integer :: my_srank                                      !!slave rank in communictor SCOMM
  integer :: n_mynod                                       !!number of nodes I'm assigned
  integer :: my_ne                                         !!number of electrodes I'm assigned
  integer :: nmy_drows                                     !!number of data in my assembly vector
  integer :: iopt = 1                                      !!flag for real or complex inversion
  integer :: iter                                          !!outer inverse iteration
  integer :: nj_rows                                       !!number of jacobian rows I own
  integer :: nsig_max                                      !!number of conductivities at max set value
  integer :: nsig_min                                      !!number of conductivities at min set value
  integer :: smin,smax                                     !!slave indices for min and max run times
  integer :: sabmin,sabmax,skmax,skmin,sjmax,sjmin
  integer :: neq                                           !!number of equipotential constraints
  !integer :: i_zpot                                        !!locate of zero potential node for no flow boundaries
  integer :: LCOMM                                         !!Task lead communicator
  integer :: M_COMM                                        !!Master communicator (for joint inversion)
  integer :: n_ogrps,ppg                                   !!number of opt groups, procs per group
  integer :: my_ogrp,my_orank                              !!my opt group, my opt rank
  integer :: nec                                           !!number of external constraints
  integer :: n_inc                                         !!number of negative boundaries in the mesh
  integer :: n_met                                         !!number of metallic inclosures in the mesh

  !Reals
  real :: Cstart,Cend,etm,etm1,etm2,Cbeg,sbt               !!timing variables
  real :: rt_min,rt_max,abmin,abmax                        !!max and min run times, A build times
  real :: kspmax,kspmin,jmax,jmin
  real :: my_frt,my_abt                                    !!my forward run time,my A build times
  real :: my_kspt                                          !!KSP solver context build time
  real :: my_jbt
  real :: mean_sig                                         !!average apparent conductivity
  real :: mean_phase                                       !!average apparent phase
  real, dimension(10) :: betalist                          !!list of reg wts (for joint inversion)
  real :: tmp_phi_0                                        !! USED in JOINT-INVERT -> MUST BE REMOVED
  !Logicals
  logical, dimension(:), allocatable :: J_on_off           !!which elements to estimate
  logical :: cflag = .true.
  logical :: first_sol = .true.
  logical :: wopt = .true.                                 !!flag to weight data and Jacobian
  logical :: ls_beta = .false.                             !!flag to do line search on beta
  logical, dimension(4) :: pcheck                          !!flag for pole a,b,m,n's
  logical :: peq=.true.
  logical :: jprnt=.true.                                  !!flag to print diag(JTJ)
  logical :: tl_flag=.false.
  logical :: invi=.false.                                  !!indicates whether currently inverting complex
  logical :: tank_flag = .false.                           !!indicates if all boundaries are no-flow
  logical :: res_flag = .false.                            !!flag to build resolution matrix
  logical :: opt_flag = .false.                            !!flag to optimize survey
  logical :: psf_flag = .false.                            !!flag to write point spread functions
  logical :: analytic = .false.
  logical :: ave_sig = .false.                             !!flag to use average apparent conductivity as starting model 
  logical :: im_fmm = .false.                              !!true if I'm a node running on the fmm side
  logical, dimension(10) :: cgmin_flag = .false.           !!cross grad joint inversion on/off flags
  logical, dimension(10) :: cgmin_flag_start = .false.     !!cross grad joint inversion on/off flags
  logical :: jaco_out_opt = .false.                        !!flag for outputting the jacobian matrix in mode 2
  logical :: jaco_ascii_opt = .false.                      !!jaco output format flag
  
  !Integer Arrays
  integer, dimension(:,:), allocatable :: s_conf           !!abmn survey configuration
  integer, dimension(:,:), allocatable :: eind             !!electrode assignments
  integer, dimension(:,:),allocatable :: jind              !!measurement assignments (jaco rows)  
  integer, dimension(:,:),allocatable :: tind
  integer, dimension(:,:), allocatable :: elements         !!elements connections
  integer, dimension(:,:), allocatable :: faces            !!face connections
  integer, dimension(:,:), allocatable :: inc                !!table of boundary numbers for infrastructure inclusions

  integer, dimension(:), allocatable :: rows,cols
  integer, dimension(:), allocatable :: trows,tcols
  integer, dimension(:), allocatable :: AI_pc,AJ_pc
  integer, dimension(:), allocatable :: A_map              !!coupling matrix mapping vector
  integer, dimension(:), allocatable :: S_map              !!Sigma mapping vector
  integer, dimension(:,:), allocatable :: nind             !!node assignments
  integer, dimension(:), allocatable :: nneq               !!number of nodes in each equipotential boundary
  integer, dimension(:), allocatable :: leq                !!last node in each equipot boundary
  integer, dimension(:), allocatable :: ex_cols_s            !!column numbers for external constraints
  integer, dimension(:), allocatable :: e_nods             !!indices of electrode nodes
  integer, dimension(:), allocatable :: my_drows           !!rows of my data assemble vector
  integer, dimension(:), allocatable :: nbounds,zones      !!node boundaries and element zones
  integer, dimension(:), allocatable :: rr_seq             !!sequence for build jacobian matrix
  integer, dimension(:), allocatable :: OCOMM,SOCOMM       !!optimization group communicators
  integer, dimension(:), allocatable :: ex_cols            !!column numbers for external constraints
  integer, dimension(:,:), allocatable :: source_nodes     !!nodes indexes of elements that contain current sources

  !Real Arrays
  real, dimension(:,:), allocatable :: e_pos,eptmp         !!electrode positions and surface flag
  real, dimension(:), allocatable :: dobs,Wd               !!observed data and weighting vector
  real, dimension(:), allocatable :: dobsi,Wdi             !!observed complex data and weighting vector
  real, dimension(:), allocatable :: Wd_not
  real, dimension(:), allocatable :: dobs_not              !!background observed data
  real, dimension(:), allocatable :: dpred,dpred_not       !!predicted data vector, background pred
  real, dimension(:), allocatable :: dpredi
  real, dimension(:), allocatable :: sigma,refsig,prefsig  !!element conductivities
  real, dimension(:), allocatable :: sigmai                !!imaginary element conductivities
  real, dimension(:), allocatable :: sig_not               !!background conductivity
  real, dimension(:,:), allocatable :: nodes               !!node positions
  real, dimension(:), allocatable :: my_dvals              !!values in my data assembly vector
  real, dimension(:), allocatable :: sig_up                !!conductivity update vector
  real, dimension(:), allocatable :: evol                  !!vector to store element volumes
  real, dimension(:), allocatable :: Precond               !!vector to store preconditioner
  real, dimension(:,:), allocatable :: Jaco                !!jacobian matrix
  real*8, dimension(:,:), allocatable :: Jacoi             !!complex component jacobian matrix
  real, dimension(:,:), allocatable :: poles               !!pole solutions
  real*8, dimension(:,:), allocatable :: polesi            !!complex component of pole solutions
  real, dimension(:), allocatable :: X1,sol1               !!pmatvec1 vectors
  real, dimension(:), allocatable :: X2,sol2               !!pmatvec2 vectors
  real*8, dimension(:), allocatable :: X1i,sol1i           !!pmatvec1i vectors
  real*8, dimension(:), allocatable :: X2i,sol2i           !!pmatvec2i vectors
  real, dimension(:), allocatable :: d_res                 !!diagonal of the resolution matrix
  real*8, dimension(:), allocatable :: Jdp                 !!real current source due to complex cond
  real, dimension(:,:), allocatable :: cg_wts              !!weights for cross gradient joint inversion
  real, dimension(:,:), allocatable :: inv_dist            !!inverse distance matrix for survey optimization
  real, dimension(:,:), allocatable :: source_currents     !!fraction of current on each node in an element containing an electrode
  
  !PETSC 
  PetscInt, dimension(:), allocatable :: d_nnz              !!petsc preallocation vector (diag blocks)
  PetscReal, dimension(:), allocatable :: delA

  Mat :: A,Ai,Ad
  PetscErrorCode :: perr
  MatType :: tp
  PetscInt :: prn(1),pcn(1)
  PetscReal :: val(1)
  PetscInt :: d_nz,o_nz
  Vec :: psol
  Vec :: X,Xd
  Vec :: B,Bd
  KSP :: KS,KSi,Kd
  PC :: P,Pd
  logical :: nzero_flag=.true.
  PetscScalar, pointer :: vloc(:)
  PetscScalar, pointer :: vlocd(:)
  PetscScalar, pointer :: aloc(:,:)
  integer :: E4D_COMM                                         !!communicator for E4D only
  

end module vars
