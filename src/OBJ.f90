module obj

  use input
  use vars
  use mod_con
  
  implicit none
  real :: phi_tot,phi_data,phi_model,chi2,PHI_D,PHI_M,PLAM,phi_cg
  real :: PHI_T = 0,errm,erms,errm0,erms0,chi20,phi_data0,o_chi2
  real , dimension(15) :: phi=0
  integer :: ncull
  logical :: con_flag = .false.
  integer, dimension(1) :: inrb
  real, dimension(:),allocatable :: imod_vec
  integer, dimension(:), allocatable :: inum_vec
  logical, dimension(4) :: check_chidec = .false.
  
  contains
   
    !______________________________________________________________________
    subroutine check_convergence
      implicit none
      real :: lambda
      logical :: exst

      con_flag = .true.
      !!get the objective function value 
      !!call eval_obj
      if(im_fmm) then
         nm = nm_fmm
         call eval_objII_fmm
      else
         call eval_objII
      end if
      
      !!Check the chi2 criteria
      con_flag = .false.
      if((chi2 <= norm_chi2)) then
         con_flag = .true.
         return
      end if
      

    end subroutine check_convergence
    !______________________________________________________________________

    !______________________________________________________________________
    subroutine eval_objII
      implicit none
      
      real, dimension(:),allocatable :: mod_vec, sig_up_par
      real, dimension(nm) :: dat_vec
      integer :: nrow,i,rbi,j,k,n,sizeData
      real :: err_mean,err_sdev
	  
  
      o_chi2=chi2
  
      if(invi) then
         dat_vec=Wdi*(dobsi-dpredi)
      else
         dat_vec = Wd*(dobs-dpred)
      end if
      phi_data0 = dot_product(dat_vec,dat_vec)
      chi20 = phi_data0/(nm)
     
      Wd_cull = 1
      ncull = 0
      if(cull_flag == 1) then
         !!calculate the error mean and standard deviation
         err_mean = sum(dat_vec)/real(nm)
         err_sdev = sqrt(dot_product( (dat_vec-err_mean),(dat_vec-err_mean))/real(nm))
         
         !!set the data culling vector
         do i=1,nm
            !irls data culling to improve stability 10/7/2018
            !Wd_cull(i) = 0.5*(1-erf((abs(dat_vec(i))-(err_mean+cull_dev*err_sdev))/sqrt(2*err_sdev**2)))

            !write(*,*) Wd_cull(i)
            if(dat_vec(i) < (err_mean - cull_dev*err_sdev) .or. dat_vec(i) > err_mean+cull_dev*err_sdev) then
               Wd_cull(i)=0
               ncull=ncull+1
            end if
         end do 
         ncull = nm-int(sum(Wd_cull))
         
         
      end if
      
      if(invi) then
         dat_vec=Wdi*Wd_cull*(dobsi-dpredi)
      else
         dat_vec = Wd*Wd_cull*(dobs-dpred)
      end if
      
      phi_data = dot_product(dat_vec,dat_vec)
      chi2 = phi_data/(nm-ncull)
   
      if(invi) then
         do i=1,nm
            dat_vec(i)=(dobsi(i)-dpredi(i))*Wdi(i) - 1 !/dobsi(i)
         end do
      else
         do i=1,nm
            dat_vec(i)=(dobs(i)-dpred(i))*Wd(i) -1     !/dobs(i)
         end do
      end if
   
      errm0= sum(dat_vec/nm)
      erms0= sqrt(dot_product(dat_vec,dat_vec)/nm)
      errm = sum(dat_vec*Wd_cull/(nm-ncull))
      erms=sqrt(dot_product(dat_vec,dat_vec)/(nm-ncull))
        
   
      if(mode == 3 .or. res_flag) then
         inrb(1)=maxval(rblock(:,3))
         if(allocated(imod_vec)) deallocate(imod_vec,inum_vec)
         allocate(imod_vec(inrb(1)),inum_vec(inrb(1)))
         allocate(mod_vec(ccount))
         mod_vec = 0
         imod_vec = 0
         inum_vec = 0
      
         do i=1,ccount
            rbi=rblock(i,3)
            inum_vec(rbi) = inum_vec(rbi) + 1
            select case(smetric(rbi,2))
            case(1)
               if(invi) then
                  mod_vec(i)=Wm(i)*(log(phase(rblock(i,1)))-log(phase(rblock(i,2))))
               else
                  mod_vec(i)=Wm(i)*(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))
               end if

            case(2)
               if(invi) then
                  !mod_vec(i)=Wm(i)*abs(log(sigmai(rblock(i,1)))-log(sigmai(rblock(i,2))))
                  mod_vec(i)=Wm(i)*abs(log(phase(rblock(i,1))) - log(phase(rblock(i,2))))
               else
                  mod_vec(i)=Wm(i)*abs(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))
               end if

            case(3)
           
               select case(smetric(rbi,3))
               case(0)
                  if(invi) then
                     mod_vec(i) = Wm(i)*(log(phase(rblock(i,1)))-(C_targ(rbi)))
                  else
                     mod_vec(i) = Wm(i)*(log(sigma(rblock(i,1)))-(C_targ(rbi)))
                  end if
                  
               case(1)
                  if(invi) then
                     mod_vec(i) = Wm(i)*(log(phase(rblock(i,1)))-log(refsig(rblock(i,1))))
                  else
                     mod_vec(i) = Wm(i)*(log(sigma(rblock(i,1)))-log(refsig(rblock(i,1))))
                  end if

               case(2)
                  if(invi) then
                     mod_vec(i) = Wm(i)*(log(phase(rblock(i,1)))-log(prefsig(rblock(i,1))))
                  else
                     mod_vec(i) = Wm(i)*(log(sigma(rblock(i,1)))-log(prefsig(rblock(i,1))))
                  end if

               end select

            case(4)
               select case(smetric(rbi,3))
               case(0)
                  if(invi) then
                     mod_vec(i) = Wm(i)*abs(log(phase(rblock(i,1)))-(C_targ(rbi)))
                  else
                     mod_vec(i) = Wm(i)*abs(log(sigma(rblock(i,1)))-C_targ(rbi))
                  end if

               case(1)
                  if(invi) then
                     mod_vec(i) = Wm(i)*abs(log(phase(rblock(i,1)))-log(refsig(rblock(i,1))))
                  else
                     mod_vec(i) = Wm(i)*abs(log(sigma(rblock(i,1)))-log(refsig(rblock(i,1))))
                  end if

               case(2)
                  if(invi) then
                     mod_vec(i) = Wm(i)*abs(log(phase(rblock(i,1)))-log(prefsig(rblock(i,1))))
                  else
                     mod_vec(i) = Wm(i)*abs(log(sigma(rblock(i,1)))-log(prefsig(rblock(i,1))))
                  end if

               end select

              
         
            case(5)
               if(invi) then
                  mod_vec(i)=Wm(i)*(log(phase(rblock(i,1)))-log(phase(rblock(i,2))))
               else
                  mod_vec(i)=Wm(i)*(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))  
               end if

            case(6)
               if(invi) then
                  mod_vec(i)=Wm(i)*abs(log(phase(rblock(i,1)))-log(phase(rblock(i,2))))
               else
                  mod_vec(i)=Wm(i)*abs(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))
               end if

            case(7)
               select case(smetric(rbi,3))
               case(0)
                  if(invi) then
                     mod_vec(i)=Wm(i)*( (C_targ(rbi) - log(phase(rblock(i,2)))) - &
                          (C_targ(rbi) - log(phase(rblock(i,1)))))
                  else
                     mod_vec(i)=Wm(i)*( (C_targ(rbi) - log(sigma(rblock(i,2)))) - &
                          (C_targ(rbi) - log(sigma(rblock(i,1)))))
                  end if

               case(1)
                  if(invi) then
                     mod_vec(i)=Wm(i)*(log(phase(rblock(i,2)))-log(refsig(rblock(i,2))) - &
                          (log(phase(rblock(i,1)))-log(refsig(rblock(i,1)))))
                  else
                     mod_vec(i)=Wm(i)*(log(sigma(rblock(i,2)))-log(refsig(rblock(i,2))) - &
                          (log(sigma(rblock(i,1)))-log(refsig(rblock(i,1)))))
                  end if

               case(2)
                  if(invi) then
                     mod_vec(i)=Wm(i)*(log(phase(rblock(i,2)))-log(prefsig(rblock(i,2))) - &
                          (log(phase(rblock(i,1)))-log(prefsig(rblock(i,1)))))
                  else
                     mod_vec(i)=Wm(i)*(log(sigma(rblock(i,2)))-log(prefsig(rblock(i,2))) - &
                          (log(sigma(rblock(i,1)))-log(prefsig(rblock(i,1)))))
                  end if

               end select


            case(8)
               select case(smetric(rbi,3))
               case(0)
                
                  if(invi) then
                     mod_vec(i)=Wm(i)*( (C_targ(rbi) - log(phase(rblock(i,2)))) - &
                          (C_targ(rbi) - log(phase(rblock(i,1)))))
                  else
                     mod_vec(i)=Wm(i)*( (C_targ(rbi) - log(sigma(rblock(i,2)))) - &
                          (C_targ(rbi) - log(sigma(rblock(i,1)))))
                  end if
                  
               case(1)
               
                  if(invi) then
                     mod_vec(i)=Wm(i)*(log(phase(rblock(i,2)))-log(refsig(rblock(i,2))) - &
                          (log(phase(rblock(i,1)))-log(refsig(rblock(i,1)))))
                  else
                     mod_vec(i)=Wm(i)*(log(sigma(rblock(i,2)))-log(refsig(rblock(i,2))) - &
                          (log(sigma(rblock(i,1)))-log(refsig(rblock(i,1)))))
                  end if
               
                  
               case(2)
                  
                  if(invi) then
                     mod_vec(i)=Wm(i)*(log(phase(rblock(i,2)))-log(prefsig(rblock(i,2))) - &
                          (log(phase(rblock(i,1)))-log(prefsig(rblock(i,1)))))
                  else
                     mod_vec(i)=Wm(i)*(log(sigma(rblock(i,2)))-log(prefsig(rblock(i,2))) - &
                          (log(sigma(rblock(i,1)))-log(prefsig(rblock(i,1)))))
                  end if
                  
               end select
               
           case(9)
              if(invi) then
                 mod_vec(i)=Wm(i)*(log(phase(rblock(i,1)))-log(phase(rblock(i,2))))
              else
                 mod_vec(i)=Wm(i)*(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))
              end if

            case(10)
               if(invi) then
                  mod_vec(i)=Wm(i)*abs(log(phase(rblock(i,1)))-log(phase(rblock(i,2))))
               else
                  mod_vec(i)=Wm(i)*abs(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))
               end if

            case(11)
               if(invi) then
                  mod_vec(i)=Wm(i)*abs(log(phase(rblock(i,1)))-log(phase(rblock(i,2))))
               else
                  mod_vec(i)=Wm(i)*abs(log(sigma(rblock(i,1)))-log(sigma(rblock(i,2))))
               end if

            case(12)
               mod_vec(i)=Wm(i)*log(sigma(rblock(i,1)))*cg_wts(rblock(i,2),1)
               do j=1,4
                  k=neighbors(rblock(i,1),j)
                  if(k>0) then
                     mod_vec(i)=mod_vec(i)+Wm(i)*log(sigma(k))*cg_wts(rblock(i,2),j+1)
                  end if
               end do
                           
         
            case DEFAULT
            end select
            imod_vec(rbi) = imod_vec(rbi) +  abs(mod_vec(i))

         end do
      end if
     
      phi_model = beta * dot_product(mod_vec,mod_vec)    
      phi_tot = phi_data + phi_model
      deallocate(mod_vec)

      if (cgmin_flag(1)) phi_cg = dot_product(tnod,tnod)
      
      if(iter==0) then
          PHI_T = phi_tot
          PHI_D = phi_data
          PHI_M = phi_model
      end if

      !!check for hemstiching and divergence
      if(iter>0 .and. iter<5) then
         hemstiching = .false.
         if(chi2>o_chi2) then
            check_chidec(iter) = .false.
         else
            check_chidec(iter) = .true.
         end if
      elseif (iter>4) then
         do i=1,3
            check_chidec(i) = check_chidec(i+1)
         end do
         check_chidec(4) = .true.
         if(chi2 > o_chi2) then
            check_chidec(4) = .false.
         end if

         if((check_chidec(2) .neqv. check_chidec(1)) .and. &
              (check_chidec(3).neqv.check_chidec(2)) .and. &
              (check_chidec(4).neqv.check_chidec(3))) then
            hemstiching = .true.
            
         end if
      end if

      if(iter<3) then
         diverging = .false.
      else
         if(.not.check_chidec(1) .and. &
              .not.check_chidec(2) .and. &
              .not.check_chidec(3)) then
            diverging = .true.
         end if
      end if
      
    end subroutine eval_objII 
    !______________________________________________________________________
    
    
    !______________________________________________________________________
     subroutine eval_objII_fmm
          implicit none
          
          real, dimension(:),allocatable :: mod_vec, sig_up_par
          real, dimension(nm) :: dat_vec
          integer :: nrow,i,rbi,j,k
          real :: err_mean,err_sdev
     

          o_chi2=chi2
          
          if(invi) then
             dat_vec=Wdi*(dobsi-dpredi)
          else
             dat_vec = Wd*(dobs-dpred)
          end if
          
          phi_data0 = dot_product(dat_vec,dat_vec)
          chi20 = phi_data0/(nm)
        
          Wd_cull = 1
          ncull = 0
          if(cull_flag == 1) then
             !!calculate the error mean and standard deviation
             err_mean = sum(dat_vec)/real(nm)
             err_sdev = sqrt(dot_product( (dat_vec-err_mean),(dat_vec-err_mean))/real(nm))
             
             !!set the data culling vector
             do i=1,nm
                if(dat_vec(i) < (err_mean - cull_dev*err_sdev) .or. dat_vec(i) > err_mean+cull_dev*err_sdev) then
                   Wd_cull(i)=0
                   ncull=ncull+1
                end if
             end do
          end if
         
          if(invi) then
             dat_vec=Wdi*Wd_cull*(dobsi-dpredi)
          else
             dat_vec = Wd*Wd_cull*(dobs-dpred)
          end if
          
          phi_data = dot_product(dat_vec,dat_vec)
          chi2 = phi_data/(nm-ncull)
         
          if(invi) then
             do i=1,nm
                dat_vec(i)=(dobsi(i)-dpredi(i))*Wdi(i) - 1 !/dobsi(i)
             end do
          else
             do i=1,nm
                dat_vec(i)=(dobs(i)-dpred(i))*Wd(i) -1     !/dobs(i)
             end do
          end if
       
          errm0= sum(dat_vec/nm)
          erms0= sqrt(dot_product(dat_vec,dat_vec)/nm)
          errm = sum(dat_vec*Wd_cull/(nm-ncull))
          erms=sqrt(dot_product(dat_vec,dat_vec)/(nm-ncull))
            
          if(mode_fmm == 3 .or. res_flag) then
             inrb(1)=maxval(rblock(:,3))
             if(allocated(imod_vec)) deallocate(imod_vec,inum_vec)
             allocate(imod_vec(inrb(1)),inum_vec(inrb(1)))
             allocate(mod_vec(ccount))
             mod_vec = 0
             imod_vec = 0
             inum_vec = 0
          
             do i=1,ccount
                rbi=rblock(i,3)
                inum_vec(rbi) = inum_vec(rbi) + 1
                select case(smetric(rbi,2))
                case(1)
                      mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))
                   
                case(2)
                      mod_vec(i)=Wm(i)*abs(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))
                 
                case(3)
                   select case(smetric(rbi,3))
                   case(0)
                         mod_vec(i) = Wm(i)*(sqrt(speed(rblock(i,1)))-(C_targ(rbi)))
                                       
                   case(1)
                         mod_vec(i) = Wm(i)*(sqrt(speed(rblock(i,1)))-(refsig(rblock(i,1))))
                        
                   case(2)
                         mod_vec(i) = Wm(i)*(sqrt(speed(rblock(i,1)))-(prefsig(rblock(i,1))))
                   
                   end select
    
                case(4)
                   select case(smetric(rbi,3))
                   case(0)
                         mod_vec(i) = Wm(i)*abs(sqrt(speed(rblock(i,1)))-C_targ(rbi))
                      
    
                   case(1)
                         mod_vec(i) = Wm(i)*abs(sqrt(speed(rblock(i,1)))-(refsig(rblock(i,1))))                     
    
                   case(2)
                         mod_vec(i) = Wm(i)*abs(sqrt(speed(rblock(i,1)))-(prefsig(rblock(i,1))))
                   
                   end select
    
                  
             
                case(5)
                      mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))  
                  
                case(6)
                      mod_vec(i)=Wm(i)*abs(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))
                   
                case(7)
                   select case(smetric(rbi,3))
                   case(0)
                         mod_vec(i)=Wm(i)*( (C_targ(rbi) - sqrt(speed(rblock(i,2)))) - &
                              (C_targ(rbi) - sqrt(speed(rblock(i,1)))))
                  
                   case(1)
                         mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,2)))-(refsig(rblock(i,2))) - &
                              (sqrt(speed(rblock(i,1)))-(refsig(rblock(i,1)))))
                     
                   case(2)
                         mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,2)))-(prefsig(rblock(i,2))) - &
                              (sqrt(speed(rblock(i,1)))-(prefsig(rblock(i,1)))))
    
                   end select
    
    
                case(8)
                   select case(smetric(rbi,3))
                   case(0)
                         mod_vec(i)=Wm(i)*( (C_targ(rbi) - sqrt(speed(rblock(i,2)))) - &
                              (C_targ(rbi) - sqrt(speed(rblock(i,1)))))
                      
                      
                   case(1)
                         mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,2)))-(refsig(rblock(i,2))) - &
                              (sqrt(speed(rblock(i,1)))-(refsig(rblock(i,1)))))
                       
                      
                   case(2)
                         mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,2)))-(prefsig(rblock(i,2))) - &
                              (sqrt(speed(rblock(i,1)))-(prefsig(rblock(i,1)))))
                      
                   end select
    
               case(9)
                     mod_vec(i)=Wm(i)*(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))
                  
    
                case(10)
                      mod_vec(i)=Wm(i)*abs(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))
                   
                case(11)
                      mod_vec(i)=Wm(i)*abs(sqrt(speed(rblock(i,1)))-sqrt(speed(rblock(i,2))))
                   
                case(12)
                   mod_vec(i)=Wm(i)*sqrt(speed(rblock(i,1)))*cg_wts(rblock(i,2),1)
                   do j=1,4
                      k=neighbors(rblock(i,1),j)
                      if(k>0) then
                         mod_vec(i)=mod_vec(i)+Wm(i)*sqrt(speed(k))*cg_wts(rblock(i,2),j+1)
                      end if
                   end do
                  
            
                case DEFAULT
                end select
                imod_vec(rbi) = imod_vec(rbi) +  abs(mod_vec(i))
    
             end do
          end if
         
          phi_model = beta * dot_product(mod_vec,mod_vec)    
          phi_tot = phi_data + phi_model
          deallocate(mod_vec)

          if (cgmin_flag(1)) phi_cg = dot_product(tnod,tnod)

          if(iter==0) then
              PHI_T = phi_tot
              PHI_D = phi_data
              PHI_M = phi_model
           end if
          
        end subroutine eval_objII_fmm
    !______________________________________________________________________

    !______________________________________________________________________
    subroutine eval_obj
      implicit none
      
      real, dimension(:),allocatable :: mod_vec, sig_up_par
      real, dimension(nm) :: dat_vec
      integer :: nrow,i
      real :: err_mean,err_sdev
 
      o_chi2=chi2
   
      sigma_par=log(sigma_par)
      if(allocated(rsigma_par))  rsigma_par=log(rsigma_par)
   
 

      nrow = maxval(wrows)      
  
      !errm =  sum(abs(1 - dpred/dobs))/nm
    

      if(invi) then
         dat_vec=Wdi*(dobsi-dpredi)
      else
         dat_vec = Wd*(dobs-dpred)
      end if
      phi_data0 = dot_product(dat_vec,dat_vec)
      chi20 = phi_data0/(nm)

      Wd_cull = 1
      ncull = 0
      if(cull_flag == 1) then
         !!calculate the error mean and standard deviation
         err_mean = sum(dat_vec)/real(nm)
         err_sdev = sqrt(dot_product( (dat_vec-err_mean),(dat_vec-err_mean))/real(nm))
         
         !!set the data culling vector
         do i=1,nm
            if(dat_vec(i) < (err_mean - cull_dev*err_sdev) .or. dat_vec(i) > err_mean+cull_dev*err_sdev) then
               Wd_cull(i)=0
               ncull=ncull+1
            end if
         end do
      end if
      
      if(invi) then
         dat_vec=Wdi*Wd_cull*(dobsi-dpredi)
      else
         dat_vec = Wd*Wd_cull*(dobs-dpred)
      end if
      
      phi_data = dot_product(dat_vec,dat_vec)
      chi2 = phi_data/(nm-ncull)
      !!chi2 = sum(abs(dat_vec))/nm


      if(invi) then
         do i=1,nm
            dat_vec(i)=(dobsi(i)-dpredi(i))*Wdi(i) - 1 !/dobsi(i)
         end do
      else
         do i=1,nm
            dat_vec(i)=(dobs(i)-dpred(i))*Wd(i) -1     !/dobs(i)
         end do
      end if
      
      errm0= sum(dat_vec/nm)
      erms0= sqrt(dot_product(dat_vec,dat_vec)/nm)
      errm = sum(dat_vec*Wd_cull/(nm-ncull))
      erms=sqrt(dot_product(dat_vec,dat_vec)/(nm-ncull))
      
   

      if(mode == 3) then
         allocate(mod_vec(nrow))
         mod_vec = 0

         do i=1,ccount
            if(reg_opt(zones(wcols(i)))<3 .or. reg_opt(zones(wcols(i)))==4 ) then
               mod_vec(wrows(i)) = mod_vec(wrows(i)) - Wm(i)*sigma_par(wcols(i))
            end if
            if(reg_opt(zones(wcols(i)))==3) then
               mod_vec(wrows(i)) = mod_vec(wrows(i)) + (rsigma_par(wcols(i))-sigma_par(wcols(i)))*Wm(i)
 
            end if
         end do
      end if
     
      phi_model = beta * dot_product(mod_vec,mod_vec)    
      phi_tot = phi_data + phi_model
      deallocate(mod_vec)
      
      sigma_par=exp(sigma_par)
      if(allocated(rsigma_par)) rsigma_par=exp(rsigma_par)
   
      
      if(iter==0) then
          PHI_T = phi_tot
          PHI_D = phi_data
          PHI_M = phi_model
       end if
      
    end subroutine eval_obj
    !______________________________________________________________________

     !____________________________________________________________________
  subroutine check_beta
    implicit none
    real :: lambda
    character*40 :: fnam
    if(im_fmm) then
       write(fnam,"(A7)") "fmm.log"
    else
       write(fnam,"(A7)") "e4d.log"
    end if

    !!Check to see if we need to change beta yet
    open(67,file=trim(fnam),status='old',position='append')
    write(67,"(A,g10.4)") "Decrease in objective function is: ",( (PHI_T - phi_tot)/PHI_T)
    close(67)
    if(abs( (PHI_T - phi_tot)/PHI_T) <= del_obj) then
       lambda = 1-(phi_data - norm_chi2*nm)/(norm_chi2*nm)
       !if(beta_red > lambda) lambda = beta_red
       lambda=beta_red
       !!if we're here then we've found the solution at this beta value
       !!if conv_opt = 2 then we're done
       if(conv_opt == 2) then
          con_flag = .true.
          return
       end if
       open(67,file=trim(fnam),status='old',position='append')
       write(67,"(A,g10.4)") "Solution did not converge at beta = ",beta
       write(67,"(A,g10.4,A,g10.4,A,g10.4)") "Decreasing beta by ",lambda," from ",beta," to ",lambda*beta
       close(67)
       beta=lambda*beta
       phi_model = lambda*phi_model
    else 
       !write(*,*) "Decrease in objective function is sufficient to continue at beta = ",beta
       open(67,file=trim(fnam),status='old',position='append')
       write(67,"(A,g10.4)") "Decrease in objective function is sufficient to continue at beta = ",beta
       close(67)
    end if
    
    
    !!Update the objective function
    PHI_T = phi_tot
    PHI_D = phi_data
    PHI_M = phi_model

  end subroutine check_beta
  !____________________________________________________________________

 
end module obj
