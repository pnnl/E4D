module master
!!____________________________________________________________________________________________  
!!  E4D
!!  Copyright  © 2014, Battelle Memorial Institute
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
!! Description: This module contains subroutines used by the master process to control program
!! execution
!! 
!!____________________________________________________________________

  use vars
  use reorder_mesh
  use forward
  use input
  use output
  use report

  implicit none
  
contains

  !_____________________________________________________________________
  subroutine send_command(com)
    !!Send a general command to the slaves
    integer :: com
    integer :: COMM
    if(im_fmm) then
       COMM = FMM_COMM
    else
       COMM = E4D_COMM
    end if
    call MPI_BCAST(com,1,MPI_INTEGER,0,COMM,ierr)
   
  end subroutine send_command
  !____________________________________________________________________
  !_____________________________________________________________________
  subroutine send_command_fmm(com)
    !!Send a general command to the slaves
    integer :: com
    
    call MPI_BCAST(com,1,MPI_INTEGER,0,FMM_COMM,ierr)
    
  end subroutine send_command_fmm
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_dists
    !read inputs and distribute the run info to the slave
    implicit none
    integer :: neven,nextra,ce,i,j,k
    integer :: npos,nneg,np_ranks,nn_ranks,tn_rank

    
    !!divide the processors between point and non-point electrodes
    tne=ne
    npos=0
    nneg=0
    tn_rank=n_rank

    if((n_rank-1)>tne) then
       tn_rank=n_rank
       n_rank=tne+1
    end if

    do i=1,tne
       if(nbounds(e_nods(i))>=0) then
          npos=npos+1
       else 
          nneg=nneg+1
       end if
    end do
  
    nn_ranks=int(real(n_rank-1)*real(nneg)/real(tne))
    np_ranks=(n_rank-1)-nn_ranks
   
    if(np_ranks == 0 .and. npos>0) then
       np_ranks=1
       nn_ranks=nn_ranks-1
    end if
    if(nn_ranks == 0 .and. nneg>0) then
       nn_ranks=1
       np_ranks=np_ranks-1
    end if


    if(allocated(eind)) deallocate(eind)
    allocate(eind(tn_rank-1,2))
    eind=0
    
    if(npos>0) then
       neven = npos/np_ranks
       nextra = npos - neven*np_ranks
       if(neven>0) then
          ce = 1
          do i=1,np_ranks-nextra
             eind(i,1) = ce
             eind(i,2) = ce+neven-1
             ce=ce+neven
          end do
          
          if(nextra>0) then
             do i=np_ranks+1-nextra,np_ranks
                eind(i,1)=ce
                eind(i,2)=ce+neven
                ce=ce+neven+1
             end do
          end if
       else
          
          do i=1,np_ranks
             eind(i,1) = i
             eind(i,2) = i
          end do
       end if
    end if

    if(nneg > 0) then
       neven = nneg/nn_ranks
       nextra = nneg - neven*nn_ranks
       if(neven>0) then
          ce = 1
          do i=1,nn_ranks-nextra
             eind(i+np_ranks,1) = ce
             eind(i+np_ranks,2) = ce+neven-1
             ce=ce+neven
          end do
          
          if(nextra>0) then
             do i=nn_ranks+1-nextra,nn_ranks
                eind(i+np_ranks,1)=ce
                eind(i+np_ranks,2)=ce+neven
                ce=ce+neven+1
             end do
          end if
       else
          
          do i=1,nn_ranks
             eind(i+np_ranks,1) = i
             eind(i+np_ranks,2) = i
          end do
       end if

       eind(np_ranks+1:n_rank-1,:) = eind(np_ranks+1:n_rank-1,:) + npos
    end if
    
    n_rank = tn_rank
      
    !!Build jaco assignments
    if(allocated(jind)) deallocate(jind)
    allocate(jind(n_rank-1,2))
    jind = 0
    if(nm < (n_rank-1)) then
       do i=1,nm
          jind(i,1:2)=i
       end do
    else
       neven = nm/(n_rank-1)
       nextra = nm - neven*(n_rank-1)
       
       jind=0
       if(neven > 0) then
          ce = 1
          do i=1,n_rank-1-nextra
             jind(i,1) = ce
             jind(i,2) = ce+neven-1
             ce = ce+neven
          end do
          
          if(nextra > 0) then
             do i=n_rank-nextra,n_rank-1    
                jind(i,1) = ce
                jind(i,2) = ce+neven
                ce=ce+neven+1
             end do
          end if
       end if
    end if
   
    !!send assignments
    call MPI_BCAST(tne, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(e_pos,4*tne ,MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(jind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(eind,2*(n_rank-1),MPI_INTEGER , 0, E4D_COMM, ierr )
   
    call nreport(64)
  end subroutine send_dists
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine setup_forward
    implicit none

    integer :: npre,i,nchr,dim,bflag,ns,itmp
    real :: jnk,jnk1 
    integer :: status(MPI_STATUS_SIZE)
    logical :: stat
    
   
    !read the node info
    call nreport(60)
    call read_nodes(stat)
    if(.not.stat) call crash_exit
    call nreport(61)
   
    !!read the element info 
    call read_elements(stat)
    if(.not.stat) call crash_exit
    call nreport(62)
   
    if(nsig .ne. nelem ) then
       call nreport(57)
       call crash_exit
    end if

    !!read the face info
    call read_faces(stat)
    if(.not.stat) call crash_exit

    !!check to see if we have any inactive elements and if so remove them
    if (any(zones == -999)) then
       call remove_inactive_elements !!see reorder mesh module
    end if

    !command slaves to do the setup
    call send_command(2)

    !send everything to slave
    call MPI_BCAST(nnodes, 1, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
    call MPI_BCAST(nbounds, nnodes, MPI_INTEGER , 0, E4D_COMM, ierr )
    call MPI_BCAST(nelem, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(nfaces, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(faces, nfaces*4,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(zones, nelem,MPI_INTEGER,0,E4D_COMM,ierr)
  
    !!send the tank_flag
    itmp=0
    if(tank_flag) itmp=1
    call MPI_BCAST(itmp,1,MPI_INTEGER,0,E4D_COMM,ierr)
   
    !get electrode nodes
    tne=ne
    call mget_electrode_nodes(stat)
    if(.not. stat) call crash_exit
   
    !send the assignments
    call send_dists   
     
  end subroutine setup_forward
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine run_forward
    implicit none
    integer :: i,ierr,nmin,nmax,ntot,itm
    real :: rtm,ts,tc,lrep,tmax,tmin,ttot
    real, dimension(2) :: pck
    integer ::  status(MPI_STATUS_SIZE)

    call send_command(5)
    call send_command(6)
    call cpu_time(ts)
    tmin=1e9
    tmax=0
    ttot=0
    
    nmin=1e9
    nmax=0
    ntot=0
    lrep=ts
    do i=1,tne
       call MPI_RECV(pck,2,MPI_REAL,MPI_ANY_SOURCE,1,E4D_COMM,status,ierr)
   
       rtm=pck(1)
       itm=int(pck(2))
       if(rtm<tmin) tmin = rtm
       if(rtm>tmax) tmax = rtm
       ttot = ttot+rtm

       if(itm<nmin) nmin=itm
       if(itm>nmax) nmax=itm
       ntot = ntot+itm

       call cpu_time(tc)
       if((tc-lrep)>30) then
          write(*,"(A15,F5.2,A11,F6.3,A8)") "  FORWARD RUN: ",100*real(i)/real(tne),"% DONE IN :",(tc-ts)/60," MINUTES"
          write(*,*) "   MIN / MAX / AVE  RUN TIMES (sec.)",tmin,"/",tmax,"/",ttot/i 
          write(*,*) "   MIN / MAX / AVE  ITERATIONS      ",nmin,"/",nmax,"/",ntot/i
          lrep=tc
       end if
      
    end do
    call cpu_time(tc)
    write(*,"(A32,F6.3,A8)") "  E4D: DONE WITH FORWARD RUN IN:",(tc-ts)/60," minutes"
    write(*,*) 
  end subroutine run_forward
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine run_forwardi
    implicit none
    integer :: i,ierr,nmin,nmax,ntot,itm,ict
    real :: rtm,ts,tc,lrep,tmax,tmin,ttot,mxcc,cc,mxerr,err
    real, dimension(2) :: pck
    integer ::  status(MPI_STATUS_SIZE)
    real, dimension(nnodes) :: dpt,dpti
    logical :: converged 
    

    tmin=1e9
    tmax=0
    ttot=0
    
    nmin=1e9
    nmax=0
    ntot=0
    lrep=ts
    ict=0
    converged = .false.
    do while(.not. converged)
       
       !update the counter
       ict = ict+1

       !solve the complex part
       call send_command(5)
       call send_command(106)
       
       call cpu_time(ts)
  
       do i=1,tne
          call MPI_RECV(pck,2,MPI_REAL,MPI_ANY_SOURCE,1,E4D_COMM,status,ierr)
          
          rtm=pck(1)
          itm=int(pck(2))
          if(rtm<tmin) tmin = rtm
          if(rtm>tmax) tmax = rtm
          ttot = ttot+rtm
          
          if(itm<nmin) nmin=itm
          if(itm>nmax) nmax=itm
          ntot = ntot+itm
          
          call cpu_time(tc)
          if((tc-lrep)>30) then
             write(*,"(A15,F5.2,A11,F6.3,A8)") "  FORWARD IRUN: ",100*real(i)/real(tne),"% DONE IN :",(tc-ts)/60," MINUTES"
             write(*,*) "   MIN / MAX / AVE  IRUN TIMES (sec.)",tmin,"/",tmax,"/",ttot/i 
             write(*,*) "   MIN / MAX / AVE  ITERATIONS      ",nmin,"/",nmax,"/",ntot/i
             lrep=tc
          end if
          
       end do
       call cpu_time(tc)
       write(*,"(A27,F6.3,A8)") "  E4D: DONE WITH FORWARD IRUN IN:",(tc-ts)/60," minutes"
     
       return

       mxcc = 0
       do i=1,nm
          cc = sigmai(i)/sigma(i)
          if(cc > mxcc) mxcc = cc
       end do
       
       !the phase threshold is hard coded here to 200 mrad before triggering an error check
       if(mxcc>=.2) then
          !get the simulated data
          call get_dpred

          if(ict > 1) then
             converged = .true.

             !if any errors are greater than 2% then we haven't converged 
             do i=1,nm
                if(abs((dpred(i)-dpt(i))/dpt(i)) > 0.02) then
                   converged = .false.
                   exit
                elseif(abs((dpredi(i)-dpti(i))/dpti(i)) > 0.02) then
                   converged = .false.
                   exit
                end if
             end do
          end if
          if(.not. converged) then
             dpt=dpred
             dpti=dpredi 
          end if
       else
          converged = .true.
       end if
       
    end do
    
  end subroutine run_forwardi
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine get_frtimes
    implicit none
    integer :: i
    real, dimension(n_rank-1) :: run_times
    integer ::  status(MPI_STATUS_SIZE)


    call send_command(15)
    rt_max = 0
    rt_min = 1e15
    smax=0
    smin=0
    do i=1,n_rank-1
       call MPI_RECV(run_times(i),1,MPI_REAL,i,0,E4D_COMM,status,ierr)
       if(run_times(i)>rt_max) then
          rt_max = run_times(i)
          smax=i
       end if
       if(run_times(i)<rt_min) then
          rt_min=run_times(i)
          smin=i
       end if
    end do

  end subroutine get_frtimes
  !____________________________________________________________________

  
  !____________________________________________________________________

  subroutine get_abtimes
    implicit none
    integer :: i
    real, dimension(n_rank-1) :: run_times
    integer ::  status(MPI_STATUS_SIZE)


    call send_command(16)
    abmax = 0
    abmin = 1e15
    sabmax=0
    sabmin=0
    do i=1,n_rank-1
       call MPI_RECV(run_times(i),1,MPI_REAL,i,0,E4D_COMM,status,ierr)
       if(run_times(i)>abmax) then
          abmax = run_times(i)
          sabmax=i
       end if
       if(run_times(i)<abmin) then
          abmin=run_times(i)
          sabmin=i
       end if
    end do

  end subroutine get_abtimes
  !____________________________________________________________________

  !____________________________________________________________________

  subroutine get_ksptimes
    implicit none
    integer :: i
    real, dimension(n_rank-1) :: run_times
    integer ::  status(MPI_STATUS_SIZE)


    call send_command(17)
    kspmax = 0
    kspmin = 1e15
    skmax=0
    skmin=0
    do i=1,n_rank-1
       call MPI_RECV(run_times(i),1,MPI_REAL,i,0,E4D_COMM,status,ierr)
       if(run_times(i)>kspmax) then
          kspmax = run_times(i)
          skmax=i
       end if
       if(run_times(i)<kspmin) then
          kspmin=run_times(i)
          skmin=i
       end if
    end do

  end subroutine get_ksptimes
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_jtimes
    implicit none
    integer :: i
    real, dimension(n_rank-1) :: run_times
    integer ::  status(MPI_STATUS_SIZE)


    call send_command(18)
    jmax = 0
    jmin = 1e15
    sjmax=0
    sjmin=0
    do i=1,n_rank-1
       call MPI_RECV(run_times(i),1,MPI_REAL,i,0,E4D_COMM,status,ierr)
       if(run_times(i)>jmax) then
          jmax = run_times(i)
          sjmax=i
       end if
       if(run_times(i)<jmin) then
          jmin=run_times(i)
          sjmin=i
       end if
    end do

  end subroutine get_jtimes
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine synch_nind
    implicit none
    integer :: i
    integer ::  status(MPI_STATUS_SIZE)
    !THIS SUBROUTINE IS NOT BEING USED (TCJ 7/10/2012)
    !receive and store the PETSC determined node assignments from the slaves
    call send_command(10)
    if(.not. allocated(nind)) then
       allocate(nind(n_rank-1,2))
    end if
    do i=1,n_rank-1
       call MPI_RECV(nind(i,1:2),2,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
    end do

  end subroutine synch_nind
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine build_sgroupm
    implicit none
    integer :: orig_group,new_group
    integer , dimension(n_rank-1) :: sranks
    integer :: i
    call send_command(22)

    !setup a new mpi group including only the slaves for forward runs
    !note this is only done in master because all processors in 
    !E4D_COMM must participate
    call MPI_COMM_GROUP(E4D_COMM,orig_group,ierr)
    do i=1,n_rank-1
       sranks(i)=i
    end do
    call MPI_GROUP_INCL(orig_group,n_rank-1,sranks,new_group,ierr)
    call MPI_COMM_CREATE(E4D_COMM,new_group,SCOMM,ierr)
    
  end subroutine build_sgroupm
  !____________________________________________________________________
  !____________________________________________________________________
  subroutine mjaco
    implicit none
    integer, dimension(n_rank-1) :: ndone
    integer :: di,i,mown,emin,emax,count,ierr,ii
    real :: ts,tc,lrep,lrep2
    integer ::  status(MPI_STATUS_SIZE)
    
    call cpu_time(ts)
    call send_command(11)
    lrep=ts
    lrep2=ts
    ndone=0
    do i=1,nm
     
       call MPI_RECV(mown,1,MPI_INTEGER,MPI_ANY_SOURCE,1,E4D_COMM,status,ierr)
       ndone(mown)=ndone(mown)+1;
       call cpu_time(tc)
       if((tc-lrep)>30) then
          !write(*,"(A15,F4.1,A11,F5.1,A8)") "BUILDING JACO: ",100*real(i)/real(nm),"% DONE IN :",(tc-ts)/60," minutes" 
          lrep=tc
       end if
       if((tc-lrep2)/60 > 5) then
          write(*,*) "JACO BUILD ELAPSED TIME: ",(tc-ts)/60," minutes"
          write(*,*) "REPORTING PROGRESS TO FILE: jaco_build_status.txt"
          open(25,file='jaco_build_status.txt',status='replace',action='write');
          write(25,*)"Elapsed Time: ",(tc-ts)/60," minutes"
          write(25,*) i,' out of ',nm,' rows complete: ',100*real(i)/real(nm),' % of total'
          write(25,*) 'Rows finished per slave listed below:'
          do ii=1,n_rank-1
             write(25,*) 'Slave ',ii, "has finished ",ndone(ii)," out of ",jind(ii,2)-jind(ii,1)+1,' rows '
          end do
          close(25)
          lrep2=tc
       end if
    end do
    call cpu_time(tc)
    write(*,"(A29,F5.1,A8)") " E4D: DONE BUILDING JACO IN:",(tc-ts)/60," minutes"

  end subroutine mjaco
  !____________________________________________________________________

 !____________________________________________________________________
  subroutine mjacoi
    implicit none
    integer, dimension(n_rank-1) :: ndone
    integer :: di,i,mown,emin,emax,count,ierr,ii
    real :: ts,tc,lrep,lrep2
    integer ::  status(MPI_STATUS_SIZE)
    
    call cpu_time(ts)
    call send_command(111)
    lrep=ts
    lrep2=ts
    ndone=0
    do i=1,nm
     
       call MPI_RECV(mown,1,MPI_INTEGER,MPI_ANY_SOURCE,1,E4D_COMM,status,ierr)
       ndone(mown)=ndone(mown)+1;
       call cpu_time(tc)
       if((tc-lrep)>30) then
          write(*,"(A15,F4.1,A11,F5.1,A8)") "BUILDING IJACO: ",100*real(i)/real(nm),"% DONE IN :",(tc-ts)/60," minutes" 
          lrep=tc
       end if
       if((tc-lrep2)/60 > 5) then
          write(*,*) "IJACO BUILD ELAPSED TIME: ",(tc-ts)/60," minutes"
          write(*,*) "REPORTING PROGRESS TO FILE: ijaco_build_status.txt"
          open(25,file='ijaco_build_status.txt',status='replace',action='write');
          write(25,*)"Elapsed Time: ",(tc-ts)/60," minutes"
          write(25,*) i,' out of ',nm,' rows complete: ',100*real(i)/real(nm),' % of total'
          write(25,*) 'Rows finished per slave listed below:'
          do ii=1,n_rank-1
             write(25,*) 'Slave ',ii, "has finished ",ndone(ii)," out of ",jind(ii,2)-jind(ii,1)+1,' rows '
          end do
          close(25)
          lrep2=tc
       end if
    end do
    call cpu_time(tc)
    write(*,"(A23,F5.1,A8)") "DONE BUILDING IJACO IN:",(tc-ts)/60," minutes"


  end subroutine mjacoi
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_sigma
    implicit none
 
    call send_command(4)
    call MPI_BCAST(nsig, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(sigma, nsig,MPI_REAL,0,E4D_COMM,ierr)
 
  end subroutine send_sigma
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_sigmai
    implicit none
 
    call send_command(104)
    call MPI_BCAST(sigmai, nsig,MPI_REAL,0,E4D_COMM,ierr)
 
  end subroutine send_sigmai
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine send_info
    implicit none

    call send_command(7)
    call MPI_BCAST(nm, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    call MPI_BCAST(s_conf, 4*nm,MPI_INTEGER,0,E4D_COMM,ierr)

  end subroutine send_info
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_excon
    implicit none
    call send_command(27)
    call MPI_BCAST(nec, 1, MPI_INTEGER, 0,E4D_COMM,ierr)
    if(nec .ne. 0) then
       call MPI_BCAST(ex_cols,nec,MPI_INTEGER,0,E4D_COMM,ierr)
    end if  
  end subroutine send_excon
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine get_dpred
    implicit none
    integer :: opt
    integer :: i,j
    integer :: nadd
    integer :: nbuff
    integer, dimension(nm*2) :: ibuff
    real, dimension(nm) :: rbuff
    integer ::  status(MPI_STATUS_SIZE)
  
    !if(opt==1) then
    !instruct slave to assemble and send the prediceted data
    call send_command(8)
 
    !allocate dpred if not already done and zero
    if(.not. allocated(dpred)) then
       allocate(dpred(nm))
    end if
    dpred = 0
    rbuff = 0
    ibuff = 1
    do i=1,n_rank-1
       call MPI_RECV(nbuff,1,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
       call MPI_RECV(ibuff(1:nbuff),nbuff,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
       call MPI_RECV(rbuff(1:nbuff),nbuff,MPI_REAL,i,0,E4D_COMM,status,ierr)
       
       do j=1,nbuff
          dpred(ibuff(j))=dpred(ibuff(j))+rbuff(j)
       end do
          
    end do
    
    !end if
   
    if(i_flag) then
       !instruct slave to assemble and send the predceted data
       call send_command(108)
     
       !allocate dpred if not already done and zero
       if(.not. allocated(dpredi)) then
          allocate(dpredi(nm))
       end if
   
       dpredi = 0
       rbuff = 0
       ibuff = 1
       do i=1,n_rank-1 
          call MPI_RECV(nbuff,1,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
          call MPI_RECV(ibuff(1:nbuff),nbuff,MPI_INTEGER,i,0,E4D_COMM,status,ierr)
          call MPI_RECV(rbuff(1:nbuff),nbuff,MPI_REAL,i,0,E4D_COMM,status,ierr)
          
          do j=1,nbuff
             dpredi(ibuff(j))=dpredi(ibuff(j))+rbuff(j)
          end do
        
       end do
  
       dpredi = atan(dpredi/dpred)
    
    end if

    !dpredi = dpredi/dpred
    !!output the predicted data to file
    
    call output_dpred
  end subroutine get_dpred
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine send_data_noise
       
    !Instructs slaves to add data noise to the Jacobian
    implicit none
    integer :: ierr,iopt
    logical :: im_fmm_flag
    integer :: COMM
    if(im_fmm) then
       COMM = FMM_COMM
       nm = nm_fmm
    else
       COMM = E4D_COMM
    end if
    call  send_command(12)
    if(invi) then
       !call MPI_BCAST(Wdi*Wd_cull,nm,MPI_REAL,0,E4D_COMM,ierr)
       call MPI_BCAST(Wdi*Wd_cull/dpred,nm,MPI_REAL,0,COMM,ierr)
    else
       call MPI_BCAST(Wd*Wd_cull,nm,MPI_REAL,0,COMM,ierr)
    end if
  
    return
  end subroutine send_data_noise
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine send_data_noisei
       
    !Instructs slaves to add data noise to the Jacobian
    implicit none
    integer :: ierr,iopt
    
    call  send_command(112)  
    call MPI_BCAST(Wdi*Wd_cull,nm,MPI_REAL,0,E4D_COMM,ierr)
    
    return
  end subroutine send_data_noisei
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine do_pmatvec2(n,m,mv2,x)
    
    !!Computes M'x in parallel. Used in PCGLS1
    implicit none
    integer :: n,m,ierr,ierr1,ierr2,i,tag,iend
    real, dimension(n) :: mv2,tsol
    real, dimension(m) :: x
    integer ::  status(MPI_STATUS_SIZE)
    integer :: COMM

    tag = 0 

    if(im_fmm) then
       COMM = FMM_COMM 
    else
       COMM = E4D_COMM
    end if
    call send_command(13)
    !!send command 8 to slaves and then broadcast x
    call MPI_BCAST(x(1:m),m,MPI_REAL,0,COMM,ierr1)
   
    call MPI_REDUCE(MPI_IN_PLACE,mv2,n,MPI_REAL,MPI_SUM,0,COMM,ierr2)
 
    return

    !!now receive each part of the solution from the nodes
    if(im_fmm) then
       iend=n_rank-1
    else
       iend=n_rank_fmm-1
    end if

    do i=1,iend
       call MPI_RECV(tsol, n, MPI_REAL,i, tag, COMM, status, ierr2)      
       mv2=mv2+tsol
       tsol = 0
    end do
    write(*,*) dot_product(mv2,mv2)
    stop
  end subroutine do_pmatvec2
  !____________________________________________________________________

    !____________________________________________________________________
  subroutine do_pmatvec2_dbl(n,m,mv2,x)
    
    !!Computes M'x in parallel. Used in PCGLS1
    implicit none
    integer :: n,m,ierr,ierr1,ierr2,i,tag
    real*8, dimension(n) :: mv2,tsol
    real, dimension(m) :: x
    integer ::  status(MPI_STATUS_SIZE)
    integer :: COMM

    tag = 0

    if(im_fmm) then
       COMM = FMM_COMM
    else
       COMM = E4D_COMM
    end if

    !!send command 8 to slaves and then broadcast x
    call send_command(213)
    !call MPI_BCAST(m,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(dble(x),m,MPI_DOUBLE_PRECISION,0,COMM,ierr1)
  
   
    mv2 = 0
    tsol = 0
    !call MPI_REDUCE(MPI_IN_PLACE,mv2,n,MPI_REAL,MPI_SUM,0,E4D_COMM,ierr2)

    call MPI_REDUCE(MPI_IN_PLACE,mv2,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM,ierr2)
    return

    !!now receive each part of the solution from the nodes
    do i=1,n_rank-1
       call MPI_RECV(tsol, n, MPI_DOUBLE_PRECISION,i, tag, COMM, status, ierr2)      
       mv2=mv2+tsol
       tsol = 0
    end do
    write(*,*) dot_product(mv2,mv2)
    stop
  end subroutine do_pmatvec2_dbl
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine do_pmatvec1(n,m,mv1,x)
    
    !!Computes Mx in parallel. Used in PCGLS1
    implicit none
    integer n,m,ierr,ierr1,ierr2,i,tag,ind1,ind2,nrow,iend
    real, dimension(m) :: mv1
    real, dimension(n) :: x
    real :: dummy
    integer, dimension(n_rank) :: recvcounts,displs
    integer ::  status(MPI_STATUS_SIZE)
    integer :: COMM

    tag = 0

    if(im_fmm) then
       COMM = FMM_COMM
    else
       COMM = E4D_COMM
    end if
    call send_command(14)
 
    call MPI_BCAST(x,n,MPI_REAL,0,COMM,ierr1)
   
    recvcounts=0
    displs=0

    if(im_fmm) then
       iend = n_rank_fmm-1
    else
       iend = n_rank-1
    end if

    do i=1,iend
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       displs(i+1)=ind1-1
       recvcounts(i+1)=nrow
    end do
    
    !!Recieve the results from the slaves
    mv1 = 0 
    call MPI_GATHERV(dummy,0,MPI_REAL,mv1,recvcounts,displs,MPI_REAL,0,COMM,ierr2)
  
    return
    

    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       call MPI_RECV(mv1(ind1:ind2),nrow, MPI_REAL,i, tag, COMM, status, ierr2)
    end do
  
  end subroutine do_pmatvec1
  !____________________________________________________________________

 !____________________________________________________________________
  subroutine do_pmatvec1_dbl(n,m,mv1,x)
    
    !!Computes Mx in parallel. Used in PCGLS1
    implicit none
    integer n,m,ierr,ierr1,ierr2,i,tag,ind1,ind2,nrow
    real, dimension(m) :: mv1
    real*8, dimension(n) :: x
    real :: dummy
    integer, dimension(n_rank) :: recvcounts,displs
    integer ::  status(MPI_STATUS_SIZE)
    integer :: COMM

    tag = 0
    if(im_fmm) then
       COMM = FMM_COMM
    else
       COMM = E4D_COMM
    end if
    
    call send_command(214)
    !call MPI_BCAST(n,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,COMM,ierr1)
    
    recvcounts=0
    displs=0
    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       displs(i+1)=ind1-1
       recvcounts(i+1)=nrow
    end do
    
    !!Recieve the results from the slaves
    mv1 = 0
    call MPI_GATHERV(dummy,0,MPI_REAL,mv1,recvcounts,displs,MPI_REAL,0,COMM,ierr2)
    return
    

    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       call MPI_RECV(mv1(ind1:ind2),nrow, MPI_REAL,i, tag, COMM, status, ierr2)
    end do
  
  end subroutine do_pmatvec1_dbl
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine do_pmatvec2i(n,m,mv2,x)
    
    !!Computes M'x in parallel. Used in PCGLS1
    implicit none
    integer :: n,m,ierr,ierr1,ierr2,i,tag
    real*8, dimension(n) :: mv2,tsol
    real*8, dimension(m) :: x
    integer ::  status(MPI_STATUS_SIZE)
    integer :: COMM

    tag = 0
    
    !!send command 113 to slaves and then broadcast x
    call send_command(113)
    !call MPI_BCAST(m,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(x,m,MPI_DOUBLE_PRECISION,0,E4D_COMM,ierr1)
    
    mv2 = 0
    tsol = 0
    call MPI_REDUCE(MPI_IN_PLACE,mv2,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,E4D_COMM,ierr2)
    return
    
    !!now receive each part of the solution from the nodes
    do i=1,n_rank-1
       call MPI_RECV(tsol, n, MPI_REAL,i, tag, E4D_COMM, status, ierr2)      
       mv2=mv2+tsol
       tsol = 0
    end do
  end subroutine do_pmatvec2i
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine do_pmatvec1i(n,m,mv1,x)
    
    !!Computes Mx in parallel. Used in PCGLS1
    implicit none
    integer n,m,ierr,ierr1,ierr2,i,tag,ind1,ind2,nrow
    real*8, dimension(m) :: mv1
    real*8, dimension(n) :: x
    real*8 :: dummy
    integer, dimension(n_rank) :: recvcounts,displs
    integer ::  status(MPI_STATUS_SIZE)
    integer :: COMM

    tag = 0

    call send_command(114)
    !call MPI_BCAST(n,1,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(x,n,MPI_DOUBLE_PRECISION,0,E4D_COMM,ierr1)
    
    recvcounts=0
    displs=0
    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       displs(i+1)=ind1-1
       recvcounts(i+1)=nrow
    end do
    
    !!Recieve the results from the slaves
    mv1 = 0
    call MPI_GATHERV(dummy,0,MPI_DOUBLE_PRECISION,mv1,recvcounts,displs,MPI_DOUBLE_PRECISION,0,E4D_COMM,ierr2)
    return
    

    do i=1,n_rank-1
       ind1=jind(i,1)
       ind2=jind(i,2)
       nrow = ind2-ind1+1
       call MPI_RECV(mv1(ind1:ind2),nrow, MPI_REAL,i, tag, E4D_COMM, status, ierr2)
    end do
  
  end subroutine do_pmatvec1i
  !____________________________________________________________________

  !____________________________________________________________________
  subroutine alloc_sigup
    if(.not. allocated(sig_up)) then
       allocate(sig_up(nelem))
    end if
  end subroutine alloc_sigup
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine update_sigma
    implicit none
    integer :: ii
    nsig_max=0
    nsig_min=0

    if(im_fmm) then
       !at this point velocity is 1/(velocity**2) (i.e slowness**2)
 
       do ii=1,nelem
          !change velocity to slowness and update (note sig_up is the update delta
          !in terms of slowness
          velocity(ii) = sqrt(velocity(ii))+sig_up(ii)
          if((1/velocity(ii)) < min_sig) velocity(ii) = 1/min_sig
          if((1/velocity(ii)) > max_sig) velocity(ii) = 1/max_sig
       end do

       !change velocity back to slowness squared for the forward computations
       velocity = velocity**2
     
    else
       
       if (hemstiching) then
          sig_up = 0.5*sig_up
       end if
       
       if(invi) then      
          do ii=1,nelem
             sigmai(ii) = exp(log(sigmai(ii)) + sig_up(ii))
             if(atan(sigmai(ii)/sigma(ii)) > max_sig) then
                sigmai(ii) = tan(max_sig)*sigma(ii)
                nsig_max=nsig_max+1
             end if
             if(atan(sigmai(ii)/sigma(ii)) < min_sig) then
                sigmai(ii) = tan(min_sig)*sigma(ii)
                nsig_min=nsig_min+1
             end if
             
          end do
       
       else
          
          do ii=1,nelem
             sigma(ii) = exp(log(sigma(ii)) + sig_up(ii))
             if(sigma(ii) > max_sig) then
                sigma(ii) = max_sig
                nsig_max=nsig_max+1
             end if
             if(sigma(ii) < min_sig) then
                sigma(ii) = min_sig
                nsig_min=nsig_min+1
             end if
          end do   
       end if

    end if
  end subroutine update_sigma
  !____________________________________________________________________

   !____________________________________________________________________
  subroutine update_sigmai
    implicit none
    integer :: ii
    nsig_max=0
    nsig_min=0
    do ii=1,nelem
       sigmai(ii) = exp(log(sigmai(ii)) + sig_up(ii))
    end do
  end subroutine update_sigmai
  !____________________________________________________________________

  

  !____________________________________________________________________
  subroutine map_sigma_par
    implicit none
    integer :: i
    do i=1,nelem
       sigma_par(par_map(i))=sigma(i)
    end do
  end subroutine map_sigma_par
  !____________________________________________________________________
  
  !____________________________________________________________________
  subroutine map_sigma_pari
    implicit none
    integer :: i
    do i=1,nelem
       sigma_par(par_map(i))=sigmai(i)
    end do
  end subroutine map_sigma_pari
  !____________________________________________________________________
  

  !__________________________________________________________________________
  subroutine start_timer
    implicit none
    call cpu_time(Cbeg)
  end subroutine start_timer
  !__________________________________________________________________________


  !__________________________________________________________________________
  subroutine get_time
    implicit none
    call cpu_time(Cend)
    etm = (Cend-Cbeg)
  end subroutine get_time
  !__________________________________________________________________________

  !__________________________________________________________________________
  subroutine msolve_forward
    implicit none
    integer :: ierr2
    integer ::  status(MPI_STATUS_SIZE)   
    
    call send_command(6)
    !call MPI_RECV(etm2, 1, MPI_REAL,1, 0, E4D_COMM, status, ierr2)
  
  end subroutine msolve_forward
  !__________________________________________________________________________

  !__________________________________________________________________________
  subroutine write_jaco
    implicit none
    integer :: i,j
    call send_command(50)
    
    do i=1, ccount
       write(99,*) wrows(i),wcols(i),Wm(i)
    end do

  end subroutine write_jaco
  !__________________________________________________________________________

  !__________________________________________________________________________
  subroutine print_jaco
    !this subroutine prints the full jacobian matrix to an ascii or text file
    implicit none
    integer :: i,j,k,nrows
    real, dimension(:,:), allocatable :: tjac
    integer :: status(MPI_STATUS_SIZE)

    call send_command(52)
    if(jaco_ascii_opt) then
       open(13,file='jacobian.txt',action='write',status='replace')
       write(13,*) nm,nelem
    else
       open(13,file='jacobian.bin',action='write',status='replace',form='unformatted')
       write(13) nm,nelem
    end if
    do i=1,n_rank-1
        nrows = jind(i,2)-jind(i,1)+1
        allocate(tjac(nrows,nelem))
        call MPI_RECV(tjac,nrows*nelem,MPI_REAL,i,i,E4D_COMM,status,ierr)
        if(jaco_ascii_opt) then
           do j=1,nrows
              write(13,*) tjac(j,:)
           end do
        else
           do j=1,nrows
              !do k=1,nelem
                 write(13) tjac(j,:)
              !end do
           end do
        end if
        deallocate(tjac)
    end do
    close(13)
  end subroutine print_jaco
  !__________________________________________________________________________

  !__________________________________________________________________________
  subroutine build_rrseq
    implicit none
    integer :: i,j,ii,jj,count,a,b,m,n,emin,emax,irow
    real :: ts,te
    integer, dimension(nm,5) :: players
    integer, dimension(5) :: lastgo,checkgo
    integer, dimension(4) :: last_conf
    logical, dimension(nm) :: dn
    logical :: found_one

    dn=.false.
    if(allocated(rr_seq)) deallocate(rr_seq)
    allocate(rr_seq(nm))
    call cpu_time(ts)
    open(27,file='jaco_build_sequence.txt',action='write',status='replace')
    

    do ii=1,nm

       a = s_conf(ii,1)
       b = s_conf(ii,2)
       m = s_conf(ii,3)
       n = s_conf(ii,4)

       do i=1,n_rank-1
          emin=eind(i,1); emax=eind(i,2)
          if((emin .le. a) .and. (emax .ge. a)) lastgo(1) = i	
          if((emin .le. b) .and. (emax .ge. b)) lastgo(2) = i 
          if((emin .le. m) .and. (emax .ge. m)) lastgo(3) = i
          if((emin .le. n) .and. (emax .ge. n)) lastgo(4) = i
          if((jind(i,1) .le. ii) .and. (jind(i,2) .ge. ii)) lastgo(5) = i
       end do
       players(ii,:) = lastgo
    end do
    

    pcheck=.true.
    do i=2,nm
       do j=1,4
          if(s_conf(i,j).ne.s_conf(i-1,j)) then
             pcheck(j)=.false.
          end if
       end do
    end do
  
    !_________________________________________________________________________
    !test implementation
    irow=0
    count=0
    do ii=1,nm
       irow=irow+1
       do i=1,n_rank-1
          if(irow .le. (jind(i,2)-jind(i,1)+1)) then
             count=count+1
             rr_seq(count)=jind(i,1)+irow-1
             write(27,*) count,rr_seq(count),players(jind(i,1)+irow-1,:)
          end if
       end do
       if(count == nm) goto 5
    end do
5   continue
    close(27)
    call cpu_time(te)
    sbt=te-ts
    
    return
    !_________________________________________________________________________

    count=1
    rr_seq(1)=1
    dn(1)=.true.
    lastgo=players(1,:)
    last_conf=s_conf(1,:)

  
10  continue
    found_one = .false.
    do i=1,nm

       if(dn(i)) goto 100
       checkgo=players(i,:)
     
       do ii=1,5
          do jj=1,5
             if(lastgo(5)==checkgo(5)) goto 100
             if(lastgo(ii)==checkgo(jj)) then
                if(.not.pcheck(ii) .and. .not. pcheck(jj)) then
                   if(last_conf(ii) .ne. s_conf(i,jj) .or. ii.ne.jj) goto 100
                end if
             end if
          end do
       end do
    
       count=count+1
       rr_seq(count)=i
       dn(i)=.true.
       found_one=.true.
       lastgo = players(i,:)
       last_conf = s_conf(i,:)
       write(27,*) i,count,players(i,:)
       100 continue
    end do

    if(count==nm) return
    if(.not. found_one) then
       do i=1,nm
          if(.not.dn(i)) then
             count=count+1
             rr_seq(count)=i
             dn(i) = .true.
             lastgo = players(i,:)
             write(27,*) i,count,players(i,:)
             goto 200
          end if
       end do
200    continue
    end if

    if(count.ne.nm) goto 10

    !final check 
    dn=.false.
    do i=1,nm
       dn(rr_seq(i)) = .true.
    end do
    do i=1,nm
       if(.not.dn(i)) write(*,*) "MEASUREMENT ",i," NOT IN rr_seq"
    end do

    call cpu_time(te)
    sbt=te-ts
    close(27)

  end subroutine build_rrseq
  !__________________________________________________________________________

  !__________________________________________________________________________
  subroutine send_rrseq
    implicit none

    call send_command(19)
    call MPI_BCAST(rr_seq,nm,MPI_INTEGER,0,E4D_COMM,ierr)
    call MPI_BCAST(pcheck,4,MPI_LOGICAL,0,E4D_COMM,ierr)
    
  end subroutine send_rrseq
  !__________________________________________________________________________

  !__________________________________________________________________________
     subroutine record_sens
       implicit none
       integer :: ierr,i
       integer :: status(MPI_STATUS_SIZE)
       !real :: xor,yor,zor
       real, dimension(:),allocatable :: sens
       real, dimension(:),allocatable :: jtmp
       
       call send_command(24)
       if (im_fmm) then
          write(*,*) ' FMM: COMPUTING J_TRANS_J'
       else
          write(*,*) ' E4D: COMPUTING J_TRANS_J'
       endif
       
       allocate(sens(nelem),jtmp(nelem))
       sens=0
       jtmp=0

       if(im_fmm) then
          do i=1,n_rank_fmm-1             
             call MPI_RECV(jtmp,nelem,MPI_REAL,i,0,FMM_COMM,status,ierr)
             sens=sens+jtmp             
          end do
       else
          do i=1,n_rank-1             
             call MPI_RECV(jtmp,nelem,MPI_REAL,i,0,E4D_COMM,status,ierr)
             sens=sens+jtmp             
          end do
       end if

       deallocate(jtmp)

       if (im_fmm) then
          open(10,file='sensitivity_fmm.txt',action='write',status='replace')
       else
          open(10,file='sensitivity.txt',action='write',status='replace')
       endif
       
       write(10,*) nelem,' 1'
       do i=1,nelem
          write(10,*) sens(i)
       end do
       close(10)

       deallocate(sens)
       
     end subroutine record_sens
     !__________________________________________________________________________
     !__________________________________________________________________________
     subroutine sadjust(opt)
       implicit none
       character*40 :: fname
       integer :: opt, jnk,i
       real, dimension(nelem) :: osig
       real :: owt,wt

       if(opt==1) then
          write(fname,"(A,I0)") "sigma.",iter-2
       elseif(opt==2) then
          write(fname,"(A,I0)") "si.",iter-1
       else
          return
       end if
       
       if(opt==1) then
          if(iter .eq. 2 .and. ave_sig) then
             osig = mean_sig
          else
             if(iter .eq. 2) then
                open(23,file=trim(sigfile),status='old',action='read');
             else
                open(23,file=trim(fname),status='old',action='read')
             end if
             read(23,*) i
             if(i .ne. nelem) write(*,*) "SOMETHING BAD IS HAPPENING IN OZ"
             do i=1,nelem
                read(23,*) osig(i)
             end do
             close(23)
          end if
       end if

       if(opt == 2) then
          open(23,file=trim(fname),status='old',action='read')
          read(23,*) i
          if(i .ne. nelem) write(*,*) "SOMETHING BAD IS HAPPENING IN OZ"
          do i=1,nelem
             read(23,*) osig(i)
          end do
          close(23)
       end if
  
       if(o_chi2 <= norm_chi2 .or. norm_chi2 <= chi2) then
          call nreport(54)
          return
       end if
       
       owt=1-(o_chi2-norm_chi2)/(o_chi2-chi2)
       wt= 1-(norm_chi2-chi2)/(o_chi2-chi2)
       sigma=wt*sigma + owt*osig
    
     end subroutine sadjust
     !__________________________________________________________________________
 
     !__________________________________________________________________________
     subroutine read_nodes(st)
       implicit none
       logical :: st
       logical :: exst
       integer :: nchr,npre
       integer :: ist,bflag,dim,jnk,jnk1
       integer :: i
       logical, dimension(:), allocatable :: nbi

       st = .true.
       if(allocated(nodes)) return
     
       !get the meshfile prefix
       nchr=len_trim(mshfile)
       do i=1,nchr
          if(mshfile(i:i)=='.') then
             npre=i+1;
             exit
          end if
       end do
          
       inquire(file=mshfile(1:npre)//".node",exist=exst)
       if(.not.exst) goto 10

       open(10,file=mshfile(1:npre)//".node",status="old",action="read")
       read(10,*,IOSTAT=ist) nnodes,dim,jnk,bflag
       if(ist .ne. 0) goto 11
       
       allocate(nodes(nnodes,3),nbounds(nnodes))
       do i=1,nnodes
          read(10,*,IOSTAT=ist) jnk,nodes(i,1:3),jnk1,nbounds(i)
          if(ist .ne. 0) goto 12
       end do
       close(10)

#ifndef imimode 		      
          do i=1,nnodes
             if(nbounds(i)<0) then
                goto 13
             end if
          end do
#endif

       allocate(nbi(nnodes))
       nbi=.false.
       do i=1,nnodes
          if(nbounds(i)<0) then
             nbi(abs(nbounds(i)))=.true.
          end if
       end do
       n_met = 0
       do i=1,nnodes
          if(nbi(i)) n_met=n_met+1
       end do

       return

10     continue
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*)
       write(51,*) ' Cannot find the node file : ',mshfile(1:npre)//'.node'
       write(51,*) ' Aborting ...'
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the node file : ',mshfile(1:npre)//'.node'
       write(*,*) ' Aborting ...'
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*)
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the node file: ',mshfile(1:npre)//'.node'
       write(51,*) ' Aborting ...'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the node file: ',mshfile(1:npre)//'.node'
       st=.false.
       return

12     continue
       close(10)
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the node file: ',mshfile(1:npre)//'.node'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the node file: ',mshfile(1:npre)//'.node'
       st=.false.
       return

13     continue
       open(51,file='e4d.log',status='old',action='write',position='append')
       write(51,*) 
       write(51,*) ' Node number ',i,' has a negative boundary flag'
       write(51,*) ' which indicates a infinite conductivity boundary.'
       write(51,*) ' Infinite conductivity boundaries are not '
       write(51,*) ' implemented in this version of e4d.'
       write(51,*) ' Aborting ...'
       close(51) 
       write(*,*) 
       write(*,*) ' Node number ',i,' has a negative boundary flag'
       write(*,*) ' which indicates a infinite conductivity boundary.'
       write(*,*) ' Infinite conductivity boundaries are not '
       write(*,*) ' implemented in this version of e4d.'
       write(*,*) ' Aborting ...'
       st = .false.
       return
     end subroutine read_nodes
     !__________________________________________________________________________

     !__________________________________________________________________________
     subroutine read_elements(st)
       implicit none
       logical :: st
       logical :: exst
       integer :: i,nchr,npre
       integer :: ist,dim,jnk,nzn
       logical, dimension(:), allocatable :: zne
 
       if(allocated(elements)) return

       st = .true.
       !get the meshfile prefix
       nchr=len_trim(mshfile)
       do i=1,nchr
          if(mshfile(i:i)=='.') then
             npre=i+1;
             exit
          end if
       end do
          
       inquire(file=mshfile(1:npre)//".ele",exist=exst)
       if(.not.exst) goto 10

       open(10,file=mshfile(1:npre)//".ele",status="old",action="read")
       read(10,*,IOSTAT=ist) nelem,dim,jnk
       if(ist .ne. 0) goto 11

       allocate(elements(nelem,4),zones(nelem))
       do i=1,nelem
          read(10,*,IOSTAT=ist) jnk,elements(i,1:4),zones(i)
          if(ist .ne. 0) goto 12
       end do
       close(10)
       
       allocate(zne(nelem))
       zne=.false.
       do i=1,nelem
          if(zones(i) .ne. -999) then
             if(zones(i) .le. 0 ) then
                nrz = i
                call nreport(63)
                st=.false.
                return
             end if
             zne(zones(i))=.true.
          end if
       end do
       nrz=0
       do i=1,nelem
          if(zne(i)) nrz=nrz+1
       end do
       deallocate(zne)
       return

10     continue
       open(51,file='e4d.log',status='old',action='write')
       write(51,*)
       write(51,*) ' Cannot find the element file : ',mshfile(1:npre)//'.ele'
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the ele file : ',mshfile(1:npre)//'.ele'
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='e4d.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       close(51)
       write(*,*)
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       st=.false.
       return

12     continue
       close(10)
       open(51,file='e4d.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the element file: ',mshfile(1:npre)//'.ele'
       st=.false.
       return

     end subroutine read_elements
     !__________________________________________________________________________

     !__________________________________________________________________________
     subroutine read_faces(st)
     implicit none
     logical :: st
     logical :: exst
     integer :: npre,nchr,i,jnk,ist

     st = .true.
     !get the meshfile prefix
     nchr=len_trim(mshfile)
     do i=1,nchr
        if(mshfile(i:i)=='.') then
           npre=i+1;
           exit
        end if
     end do
     
     inquire(file=mshfile(1:npre)//".face",exist=exst)
     if(.not.exst) goto 10

     open(10,file=mshfile(1:npre)//".face",status="old",action="read")
     
     read(10,*,IOSTAT=ist) nfaces,jnk
     if(ist .ne. 0) goto 11

     allocate(faces(nfaces,4))
     do i=1,nfaces
        read(10,*,IOSTAT=ist) jnk,faces(i,1:4)
        if(ist .ne. 0) goto 12
     end do
     close(10)
     return

10   continue
       open(51,file='e4d.log',status='old',action='write')
       write(51,*)
       write(51,*) ' Cannot find the face file : ',mshfile(1:npre)//'.face'
       write(51,*) " Aborting ..."
       close(51)
       write(*,*)
       write(*,*) ' Cannot find the face file : ',mshfile(1:npre)//'.face'
       write(*,*) " Aborting ..."
       st=.false.
       return
       
11     continue
       close(10)
       open(51,file='e4d.log',status='old',action='write')
       write(51,*) 
       write(51,*) ' There was a problem reading the first line'
       write(51,*) ' of the face file: ',mshfile(1:npre)//'.face'
       write(*,*) " Aborting ..."
       close(51)
       write(*,*) 
       write(*,*) ' There was a problem reading the first line'
       write(*,*) ' of the face file: ',mshfile(1:npre)//'.face'
       write(*,*) " Aborting"
       st=.false.
       return

12     continue
       close(10)
       open(51,file='e4d.log',status='old',action='write')
       write(51,*)
       write(51,*) ' There was a problem reading line: ',i
       write(51,*) ' of the face file: ',mshfile(1:npre)//'.face'
       close(51)
       write(*,*)
       write(*,*) ' There was a problem reading line :',i
       write(*,*) ' of the face file: ',mshfile(1:npre)//'.face'
       st=.false.
       return

   end subroutine read_faces
     !__________________________________________________________________________

   !____________________________________________________________________________
   subroutine get_phase
     implicit none

     if(.not.allocated(phase)) allocate(phase(nelem))
     phase=atan(sigmai/sigma)

   end subroutine get_phase
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine update_sigi
     implicit none

     sigmai = tan(phase)*sigma
   end subroutine update_sigi
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine send_J_on_off
     implicit none
     integer :: COMM
     if(im_fmm) then
        COMM = FMM_COMM
     else
        COMM = E4D_COMM
     end if
     
     call send_command(51)
     call MPI_BCAST(J_on_off,nelem,MPI_LOGICAL,0,COMM,ierr)
    
   end subroutine send_J_on_off
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine sync_joint
     !determine if the other groups are running joint inversions
     implicit none
     integer :: ierr,tag
     integer ::  status(MPI_STATUS_SIZE)
     logical, dimension(10) :: tflags
 
     tag = 0
     if(im_fmm) then
        call MPI_SEND(cgmin_flag,10,MPI_LOGICAL,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_RECV(tflags,10,MPI_LOGICAL,0,tag,MPI_COMM_WORLD,status,ierr)
     else
        call MPI_RECV(tflags,10,MPI_LOGICAL,master_proc_fmm,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(cgmin_flag,10,MPI_LOGICAL,master_proc_fmm,tag,MPI_COMM_WORLD,ierr)
     end if

     !if cgmin_flag(2) is true then the other group is doing a cross grad joint inversion
     !if(tflags(1)) cgmin_flag(2) = .true.
     cgmin_flag(2) = tflags(1)
     
   end subroutine sync_joint
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine sync_zwts
     implicit none
     integer :: ierr,tag,i,j,cnt
     integer :: status(MPI_STATUS_SIZE)
     real, dimension(:,:), allocatable :: twts,tvals

     cnt = 0
     do i=1,nrz
        if(smetric(i,2).eq.12 .or. smetric(i,2).eq.13) then
           cnt = cnt+1
        end if
     end do

     allocate(twts(cnt,6),tvals(cnt,6))
     cnt = 0
     do i=1,nrz
        if(smetric(i,2).eq.12 .or. smetric(i,2).eq.13) then
           cnt = cnt+1
           twts(cnt,1) = smetric(i,1)
           twts(cnt,2) = smetric(i,2)
           twts(cnt,3) = zwts(i,1)
           twts(cnt,4) = zwts(i,2)
           twts(cnt,5) = zwts(i,3)
           twts(cnt,6) = zwts(i,4)
        end if
     end do

     tag = 0
     if(im_fmm) then
        call MPI_SEND(twts,cnt*6,MPI_REAL,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_RECV(tvals,cnt*6,MPI_REAL,0,tag,MPI_COMM_WORLD,status,ierr)
     else
        call MPI_RECV(tvals,cnt*6,MPI_REAL,master_proc_fmm,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(twts,cnt*6,MPI_REAL,master_proc_fmm,tag,MPI_COMM_WORLD,ierr)
     end if

     do i=1,nrz
        if(smetric(i,2).eq.12 .or. smetric(i,2).eq.13) then
           do j=1,cnt
              if(tvals(j,1).eq.smetric(i,1) .and. tvals(j,2).eq.smetric(i,2)) then
                 zwts(i,1) = sqrt(zwts(i,1)*tvals(j,3))
                 zwts(i,2) = sqrt(zwts(i,2)*tvals(j,4))
                 zwts(i,3) = sqrt(zwts(i,3)*tvals(j,5))
                 zwts(i,4) = sqrt(zwts(i,4)*tvals(j,6))
                 exit
              end if
           end do
        end if
     end do

     deallocate(twts,tvals)
   end subroutine sync_zwts
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine sync_beta
     implicit none
     integer :: ierr,tag
     integer :: status(MPI_STATUS_SIZE)
     real, dimension(10) :: tvals

     betalist = 0
     betalist(1) = beta
     tag = 0
     if(im_fmm) then
        call MPI_SEND(betalist,10,MPI_REAL,0,tag,MPI_COMM_WORLD,ierr)
        call MPI_RECV(tvals,10,MPI_REAL,0,tag,MPI_COMM_WORLD,status,ierr)
     else
        call MPI_RECV(tvals,10,MPI_REAL,master_proc_fmm,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(betalist,10,MPI_REAL,master_proc_fmm,tag,MPI_COMM_WORLD,ierr)
     end if
     !betalist(2) holds the regularization weight of the other group
     betalist(2) = tvals(1)
        
   end subroutine sync_beta
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine sync_convergence
     call MPI_ALLREDUCE(MPI_IN_PLACE,con_flag,1,MPI_LOGICAL,MPI_LAND,M_COMM,ierr)
   end subroutine sync_convergence
   !____________________________________________________________________________  
   
   !____________________________________________________________________________
   subroutine get_other_dists
     !this subroutine gathers the element distributions from the other
     !forward models base base on the flag values in cgmin_flag
     implicit none
     integer :: ierr,tag
     integer ::  status(MPI_STATUS_SIZE)

     tag = 0
     if(cgmin_flag(1)) then
        if(.not. allocated(sigma)) allocate(sigma(nelem))
        if(.not. allocated(velocity)) allocate(velocity(nelem))
     end if
     if(cgmin_flag(2)) then
        call MPI_BCAST(sigma,nelem,MPI_REAL,0,M_COMM,ierr)
        call MPI_BCAST(velocity,nelem,MPI_REAL,1,M_COMM,ierr)
     end if     
     return
     
     if(im_fmm) then
      
        if(cgmin_flag(1)) then
           if(allocated(sigma)) deallocate(sigma)
           allocate(sigma(nelem))
           call MPI_RECV(sigma,nelem,MPI_REAL,0,tag,MPI_COMM_WORLD,status,ierr)
           !write(*,*) "FMM got sigma",maxval(sigma),minval(sigma)
        end if

        if(cgmin_flag(2)) then
           !send velocity to the E4D master
           call MPI_SEND(velocity,nelem,MPI_REAL,0,tag,MPI_COMM_WORLD,ierr)
        end if

     else
       
        if(cgmin_flag(2)) then
           !send sigma to the FMM master
           call MPI_SEND(sigma,nelem,MPI_REAL,master_proc_fmm,tag,MPI_COMM_WORLD,ierr)
        end if

        if(cgmin_flag(1)) then
           !get velocity from the FMM master
           if(allocated(velocity)) deallocate(velocity)
           allocate(velocity(nelem))
           call MPI_RECV(velocity,nelem,MPI_REAL,master_proc_fmm,tag,MPI_COMM_WORLD,status,ierr)
           !write(*,*) "E4D got velocity",maxval(velocity),minval(velocity)
        end if

     end if
 
   end subroutine get_other_dists
   !____________________________________________________________________________

   !____________________________________________________________________________
   subroutine build_inv_dist
     !send the slave the nodes and elements and instruct them to build and share
     !the inverse distance matrix and volume vectors
     call send_command(1009)
     call MPI_BCAST(nodes, nnodes*3, MPI_REAL , 0, E4D_COMM, ierr )
     call MPI_BCAST(elements, nelem*4,MPI_INTEGER,0,E4D_COMM,ierr)
     
     
   end subroutine build_inv_dist
   !____________________________________________________________________________
   !____________________________________________________________________________
   subroutine build_precond
     implicit none
     integer :: i,j,n_rank_comm
     real, dimension(:), allocatable :: sens
     integer :: COMM,status(MPI_STATUS_SIZE)

     if(.not. allocated(Precond)) then
        allocate(Precond(nelem))
     end if

     if(im_fmm) then
        COMM = FMM_COMM
        n_rank_comm = n_rank_fmm
     else
        COMM = E4D_COMM
        n_rank_comm = n_rank
     end if
     allocate(sens(nelem))
     sens = 0.

     call send_command(24)
     Precond = 0
     do i=1,n_rank_comm-1
        call MPI_RECV(sens,nelem,MPI_REAL,i,0,COMM,status,ierr)
        !do j=1,nelem
        !   Precond(j) = Precond(j)+sens(j)
        !end do
        Precond = Precond + sens
     end do
     
     deallocate(sens)
   end subroutine build_precond
   !____________________________________________________________________________
   
   !____________________________________________________________________________
   subroutine check_meshfiles
        logical :: st  
          
		call read_nodes(st)
		if (st) then
			call read_faces(st)
		end if
		if (st) then
		    call read_elements(st)
		end if
     
     
   end subroutine check_meshfiles
   !____________________________________________________________________________

   !*********************************************************************************
   !* validate joint inversion
   !*********************************************************************************
   subroutine validate_jointInv
     implicit none

     ! local variables
     integer :: tag,lenchar,mnchar1,mnchar2
     integer ::  status(MPI_STATUS_SIZE)
     integer ::  jnk1, jnk2, jnk3, jnk4
     integer :: tjnk1,tjnk2,tjnk3,tjnk4
     character(len=len(mshfile)) :: tmpmshfile


     ! START---

     if (.not.cgmin_flag(1).and..not.cgmin_flag(2)) then
        open(51,file='e4d.log',status='old',action='write',position='append')                     
        write(51,*) '======================================================================='
        write(51,*) '                  E4D: RUNNING IN SEPARATE INVERSION MODE ... '
        write(51,*) '======================================================================='
        close(51)
        write(*, *) '======================================================================='
        write(*, *) '                  E4D: RUNNING IN SEPARATE INVERSION MODE ... '
        write(*, *) '======================================================================='      
     elseif (cgmin_flag(1).and..not.cgmin_flag(2)) then
        open(51,file='e4d.log',status='old',action='write',position='append')
        !write(51,*) '=============================== WARNING ==============================='
        !write(51,*) ' E4D inverse option file ',trim(invfile),' requests joint inversion'
        !write(51,*) ' but FMM inverse option file does not request joint inversion!     '
        !write(51,*) ' NO JOINT INVERSION WILL BE PERFORMED.'
        !write(51,*) '                  E4D: RUNNING IN SEPARATE INVERSION MODE ... '
        !write(51,*) '======================================================================='

        write(51,*) '================================ ERROR ================================'
        write(51,*) ' E4D inverse option file ',trim(invfile),' requests joint inversion'
        write(51,*) ' but FMM inverse option file does not request joint inversion!     '
        write(51,*) ' NO INVERSION WILL BE PERFORMED...'        
        write(51,*) ' Please make sure BOTH inverse option files have joint inversion '
        write(51,*) ' constraint using structural metric code (s_met) = 12 or 13 for'
        write(51,*) ' running in joint inversion mode. Otherwise, no inverse option file'
        write(51,*) ' should have s_met = 12 or 13 to rather run in seperate inversion mode.'
        write(51,*) ' E4D: Aborting ... '
        write(51,*) '======================================================================='
        close(51)
        !write(*, *) '=============================== WARNING ==============================='
        !write(*, *) ' E4D inverse option file ',trim(invfile),' requests joint inversion'
        !write(*, *) ' but FMM inverse option file does not request joint inversion!     '
        !write(*, *) ' NO JOINT INVERSION WILL BE PERFORMED.'
        !write(*, *) '                  E4D: RUNNING IN SEPARATE INVERSION MODE ... '
        !write(*, *) '======================================================================='        

        write(*, *) '================================ ERROR ================================'
        write(*, *) ' E4D inverse option file ',trim(invfile),' requests joint inversion'
        write(*, *) ' but FMM inverse option file does not request joint inversion!     '
        write(*, *) ' NO INVERSION WILL BE PERFORMED...'        
        write(*, *) ' Please make sure BOTH inverse option files have joint inversion '
        write(*, *) ' constraint using structural metric code (s_met) = 12 or 13 for'
        write(*, *) ' running in joint inversion mode. Otherwise, no inverse option file'
        write(*, *) ' should have s_met = 12 or 13 to rather run in seperate inversion mode.'
        write(*, *) ' E4D: Aborting ... '
        write(*, *) '======================================================================='
        call crash_exit
        
        !reset cgmin_flag
        cgmin_flag = .false.        
        
     elseif (.not.cgmin_flag(1).and.cgmin_flag(2)) then
        open(51,file='e4d.log',status='old',action='write',position='append')
        ! write(51,*) '=============================== WARNING ==============================='
        ! write(51,*) ' E4D inverse option file ',trim(invfile),' does not request joint'
        ! write(51,*) ' inversion but FMM inverse option file requests joint inversion. '
        ! write(51,*) ' NO JOINT INVERSION WILL BE PERFORMED.'
        ! write(51,*) '                  E4D: RUNNING IN SEPARATE INVERSION MODE ... '
        ! write(51,*) '======================================================================='

        write(51,*) '================================ ERROR ================================'
        write(51,*) ' E4D inverse option file ',trim(invfile),' does not request joint'
        write(51,*) ' inversion but FMM inverse option file requests joint inversion. '
        write(51,*) ' NO INVERSION WILL BE PERFORMED...'
        write(51,*) ' Please make sure BOTH inverse option files have joint inversion '
        write(51,*) ' constraint using structural metric code (s_met) = 12 or 13 for'
        write(51,*) ' running in joint inversion mode. Otherwise, no inverse option file'
        write(51,*) ' should have s_met = 12 or 13 to rather run in seperate inversion mode.'
        write(51,*) ' E4D: Aborting ... '
        write(51,*) '======================================================================='
        close(51)

        ! write(*, *) '=============================== WARNING ==============================='
        ! write(*, *) ' E4D inverse option file ',trim(invfile),' does not request joint'
        ! write(*, *) ' inversion but FMM inverse option file requests joint inversion. '
        ! write(*, *) ' NO JOINT INVERSION WILL BE PERFORMED.'
        ! write(*, *) '                  E4D: RUNNING IN SEPARATE INVERSION MODE ... '
        ! write(*, *) '======================================================================='

        write(*, *) '================================ ERROR ================================'
        write(*, *) ' E4D inverse option file ',trim(invfile),' does not request joint'
        write(*, *) ' inversion but FMM inverse option file requests joint inversion. '
        write(*, *) ' NO INVERSION WILL BE PERFORMED...'
        write(*, *) ' Please make sure BOTH inverse option files have joint inversion '
        write(*, *) ' constraint using structural metric code (s_met) = 12 or 13 for'
        write(*, *) ' running in joint inversion mode. Otherwise, no inverse option file'
        write(*, *) ' should have s_met = 12 or 13 to rather run in seperate inversion mode.'
        write(*, *) ' E4D: Aborting ... '
        write(*, *) '======================================================================='
        call crash_exit
        
        !reset cgmin_flag
        cgmin_flag = .false.
        
     elseif (cgmin_flag(1).and.cgmin_flag(2)) then
        open(51,file='e4d.log',status='old',action='write',position='append')                     
        write(51,*) '======================================================================='
        write(51,*) '                  E4D: RUNNING IN JOINT INVERSION MODE ... '
        write(51,*) '======================================================================='
        close(51)
        write(*, *) '======================================================================='
        write(*, *) '                  E4D: RUNNING IN JOINT INVERSION MODE ... '
        write(*, *) '======================================================================='

        tag = 0
        lenchar = len(mshfile)
        call MPI_RECV(tmpmshfile,lenchar,MPI_CHARACTER,master_proc_fmm,tag,MPI_COMM_WORLD,status,ierr)
        call MPI_SEND(mshfile,   lenchar,MPI_CHARACTER,master_proc_fmm,tag,MPI_COMM_WORLD,ierr)

        ! Check if both physics have same mesh files
        mnchar1=index(mshfile,'.')
        mnchar2=index(tmpmshfile,'.')

        if (mshfile(1:mnchar1) .ne. tmpmshfile(1:mnchar2)) then
           open(51,file='e4d.log',status='old',action='write',position='append')
           write(51,*) '------------------------------- WARNING -------------------------------'
           write(51,*) ' E4D mesh file name ',trim(mshfile),' is different from '
           write(51,*) ' FMM mesh file name ',trim(tmpmshfile)
           write(51,*) ' Checking if they both have the same number of nodes, elements, and faces.'
           write(51,*) '-----------------------------------------------------------------------'
           close(51)
           write(*, *) '------------------------------- WARNING -------------------------------'
           write(*, *) ' E4D mesh file name ',trim(mshfile),' is different from '
           write(*, *) ' FMM mesh file name ',trim(tmpmshfile)
           write(*, *) ' Checking if they both have the same number of nodes, elements, and faces.'
           write(*, *) '-----------------------------------------------------------------------'

           ! Checking number of nodes
           open(331,file=   mshfile(1:mnchar1+1)//".node",status="old",action="read")
           open(332,file=tmpmshfile(1:mnchar2+1)//".node",status="old",action="read")

           read(331,*)  jnk1, jnk2, jnk3, jnk4
           read(332,*) tjnk1,tjnk2,tjnk3,tjnk4

           close(331); close(332)

           if (jnk1 .ne. tjnk1) then
              open(51,file='e4d.log',status='old',action='write',position='append')
              write(51,*) '-------------------------------- ERROR --------------------------------'
              write(51,*) ' E4D node file ',mshfile(1:mnchar1+1)//".node",' has different number '
              write(51,*) ' of nodes than FMM node file ',tmpmshfile(1:mnchar2+1)//".node"
              write(51,*) ' E4D node file has ',jnk1,' nodes but FMM has ',tjnk1,' nodes.'
              write(51,*) ' E4D node file should have the same number of nodes as FMM node file  '
              write(51,*) ' for the cross-gradient joint inversion.'
              write(51,*) ' Aborting ...'
              write(51,*) '-----------------------------------------------------------------------'
              close(51)
              write(*, *) '-------------------------------- ERROR --------------------------------'
              write(*, *) ' E4D node file ',mshfile(1:mnchar1+1)//".node",' has different number '
              write(*, *) ' of nodes than FMM node file ',tmpmshfile(1:mnchar2+1)//".node"
              write(*, *) ' E4D node file has ',jnk1,' nodes but FMM has ',tjnk1,' nodes.'
              write(*, *) ' E4D node file should have the same number of nodes as FMM node file  '
              write(*, *) ' for the cross-gradient joint inversion.'
              write(*, *) ' Aborting ...'
              write(*, *) '-----------------------------------------------------------------------'
              call crash_exit
           else
              open(51,file='e4d.log',status='old',action='write',position='append')              
              write(51,*) ' E4D: Passed nodes - same number of nodes.'
              close(51)
              write(*, *) ' E4D: Passed nodes - same number of nodes.'   
           endif

           ! checking number of elements
           open(331,file=   mshfile(1:mnchar1+1)//".ele",status="old",action="read")           
           open(332,file=tmpmshfile(1:mnchar2+1)//".ele",status="old",action="read")           

           read(331,*)  jnk1, jnk2, jnk3
           read(332,*) tjnk1,tjnk2,tjnk3

           close(331); close(332)
        
           if (jnk1 .ne. tjnk1) then
              open(51,file='e4d.log',status='old',action='write',position='append')
              write(51,*) ' -------------------------------- ERROR --------------------------------'
              write(51,*) ' E4D element file ',mshfile(1:mnchar1+1)//".ele",' has different number '
              write(51,*) ' of elements than FMM element file ',tmpmshfile(1:mnchar2+1)//".ele"
              write(51,*) ' E4D element file has ',jnk1,' elements but FMM has ',tjnk1,' elements.'
              write(51,*) ' E4D element file should have the same number of elements as FMM element'
              write(51,*) ' file for the cross-gradient joint inversion.'           
              write(51,*) ' Aborting ...'
              write(51,*) '-----------------------------------------------------------------------'
              close(51)
              write(*, *) '-------------------------------- ERROR --------------------------------'
              write(*, *) ' E4D element file ',mshfile(1:mnchar1+1)//".ele",' has different number '
              write(*, *) ' of elements than FMM element file ',tmpmshfile(1:mnchar2+1)//".ele"
              write(*, *) ' E4D element file has ',jnk1,' elements but FMM has ',tjnk1,' elements.'
              write(*, *) ' E4D element file should have the same number of elements as FMM element'
              write(*, *) ' file for the cross-gradient joint inversion.'           
              write(*, *) ' Aborting ...'
              write(*, *) '-----------------------------------------------------------------------'
              call crash_exit
           else
              open(51,file='e4d.log',status='old',action='write',position='append')              
              write(51,*) ' E4D: Passed elements - same number of elements.'
              close(51)
              write(*, *) ' E4D: Passed elements - same number of elements.'          
           endif
           
           ! Checking number of faces
           open(331,file=   mshfile(1:mnchar1+1)//".face",status="old",action="read")
           open(332,file=tmpmshfile(1:mnchar2+1)//".face",status="old",action="read")
 
           read(331,*)  jnk1, jnk2
           read(332,*) tjnk1,tjnk2

           close(331); close(332)
        
           if (jnk1 .ne. tjnk1) then
              open(51,file='e4d.log',status='old',action='write',position='append')
              write(51,*) '-------------------------------- ERROR --------------------------------'
              write(51,*) ' E4D face file ',mshfile(1:mnchar1+1)//".face",' has different number '
              write(51,*) ' of faces than FMM face file ',tmpmshfile(1:mnchar2+1)//".face"
              write(51,*) ' E4D face file has ',jnk1,' faces but FMM has ',tjnk1,' faces.'
              write(51,*) ' E4D face file should have the same number of faces as FMM face file'
              write(51,*) ' for the cross-gradient joint inversion.'           
              write(51,*) ' Aborting ...'
              write(51,*) '-----------------------------------------------------------------------'
              close(51)
              write(*, *) '-------------------------------- ERROR --------------------------------'
              write(*, *) ' E4D face file ',mshfile(1:mnchar1+1)//".face",' has different number '
              write(*, *) ' of faces than FMM face file ',tmpmshfile(1:mnchar2+1)//".face"
              write(*, *) ' E4D face file has ',jnk1,' faces but FMM has ',tjnk1,' faces.'
              write(*, *) ' E4D face file should have the same number of faces as FMM face file'
              write(*, *) ' for the cross-gradient joint inversion.'  
              write(*, *) ' Aborting ...'
              write(*, *) '-----------------------------------------------------------------------'
              call crash_exit
           else
              open(51,file='e4d.log',status='old',action='write',position='append')              
              write(51,*) ' E4D: Passed faces - same number of faces.'
              close(51)
              write(*, *) ' E4D: Passed faces - same number of faces.'
           endif           
        endif ! mesh file name are different
        
     endif ! cgmin(1) & cgmin(2) -> true
              
   end subroutine validate_jointInv
   
end module master
