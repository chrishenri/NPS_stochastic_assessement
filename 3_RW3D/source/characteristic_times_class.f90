
!********************************************************************************************
!   CHARACTERISTIC TIMES CLASS
!********************************************************************************************
  module ctimes_class

  implicit none

  private
  
  public :: ctimes_cl
  public ::                                                  &

            !calculate_characteristic_times_                , &
            update_characteristic_times_                   , &
            set_characteristic_times_                      , &
            !calculate_adv_characteristic_time_             , &
            !calculate_rx_characteristic_times_             , &
            assignment(=)                                  , &
            print_characteristic_times_                    , &
            initialize_characteristic_times_               , &
            reinitialize_characteristic_times     

  type ctimes_cl
	  real*8,  pointer  :: adv        => null()  !characteristic time for advection 
      real*8,  pointer  :: disp       => null()  !characteristic time for dispersion 
      real*8,  pointer  :: kreact     => null()  !characteristic time for kinetic reactions
      real*8,  pointer  :: decay      => null()  !characteristic time for decay network
      real*8,  pointer  :: mass_trans => null()  !characteristic time for multirate mass transfer
  end type

  interface assignment (=)
       module procedure set_characteristic_times_ ; end interface

  contains

    subroutine initialize_characteristic_times_ (tc,plume,geo,advection,dispersion,reaction)
         use global_variables, only: UNEST
         use array_class
         use advection_class
         use geometry_class
         use dispersion_class
         use reaction_class
         use plume_class
         use list_class
         use to_solve
         implicit none

		 type(ctimes_cl),            intent(inout) :: tc 
		 type(plume_cl),             intent(in)    :: plume
		 type(geometry_cl),          intent(in)    :: geo
		 type(advection_cl),         intent(in)    :: advection
		 type(dispersion_cl),        intent(in)    :: dispersion
		 type(reaction_cl),          intent(in)    :: reaction	   
		 integer                                   :: ip,irxn,ip1,ispe,izone
		 real*8                                    :: vp(3),rpt,rpt1,poro,dd,dmax,Dxx,Dyy,Dzz,kf,dx(3),kfmax
		 real*8                                    :: tcadv,tcdisp,tcreact,tcadvmin,tcdispmin,tcreactmin
		 real*8                                    :: ds,dm,dmTh,dmTV,al,ath,atv,DL,DTH,DTV,vv,D(3)
         logical, save                             :: FromVeloField =.TRUE.
                    
                     tcadvmin   = UNEST
                     tcdispmin  = UNEST
                     tcreactmin = UNEST
                     
                     dx(1) = average_array_ (geo%dx)
                     dx(2) = average_array_ (geo%dy)
                     dx(3) = average_array_ (geo%dz)		    
                     rpt = 1.d0
		             if (kineticACTION) then
		                 if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
		                    rpt   = average_array_ (reaction%sorption%LinearSorption%R(1))
		                    do ispe=2,plume%nspecie
		                        rpt1 = average_array_ (reaction%sorption%LinearSorption%R(ispe))
		                        if (rpt < rpt1) rpt = rpt1
		                    end do
		                 else
		                 rpt = 1.d0
		                 end if
		             end if      	             
      	             if (advection%action) then
                        vp(1) = average_array_ (advection%qx)
		                vp(2) = average_array_ (advection%qy)
		                vp(3) = average_array_ (advection%qz)
		                poro  = average_array_ (advection%poro)
                        vp = dabs(vp)/rpt/poro
                        tcadvmin = minval(dx/vp,mask=vp>0.)
                     end if
                     if (dispersion%action) then
 	                    dm   = average_array_ (dispersion%dm)
	                    dmTH = average_array_ (dispersion%dmTH)
	                    dmTV = average_array_ (dispersion%dmTV)
	                    if (advection%action) then 	           
                           vv = dsqrt(vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3))
 	                       al  = average_array_ (dispersion%aL)  
	                       ath = average_array_ (dispersion%aTH) 
	                       atv = average_array_ (dispersion%aTV)
	                       D(1) = ( al  * vv + dm   ) / rpt
                           D(2) = ( ath * vv + dmTH ) / rpt
                           D(3) = ( atv * vv + dmTV ) / rpt 
                        else
                           D(1) = ( dm   ) / rpt
                           D(2) = ( dmTH ) / rpt
                           D(3) = ( dmTV ) / rpt
                        end if    
	                    tcdispmin = minval(dx*dx/D,mask=D>0.)
		             end if
                     if (kineticACTION) then
                        kfmax = average_array_ (reaction%kinetic%kf(1))
                        tcreactmin = 1.d0/kfmax
                        do irxn=2,reaction%kinetic%nreact
                           kf = average_array_ (reaction%kinetic%kf(irxn))
                           tcreact = 1.d0/kf
 		                   if (ip ==1) then
 		                      tcreactmin = tcreact 
 		                   elseif (tcreact < tcreactmin ) then
 		                      tcreactmin = tcreact
 		                   end if                
                        end do
                     end if
!    !    mass transfer


         if (advection%action) then
              if (.not.associated(tc%adv)) allocate (tc%adv)
              tc % adv  = tcadvmin
         end if
         if (dispersion%action) then
              if (.not.associated(tc%disp)) allocate (tc%disp)
              tc % disp  = tcdispmin
         end if
         if (kineticACTION) then
              if (.not.associated(tc%kreact)) allocate (tc%kreact)
              tc % kreact  = tcreactmin
         end if

    
    end subroutine
!*************************************************************************************
    subroutine print_characteristic_times_ (tc,fname)
	  use gslib, only: open_fname,generate_unit
	  use global_variables, only: UNEST
	  implicit none
	  type(ctimes_cl),             intent(in) :: tc
	  character(len=*), optional,  intent(in) :: fname
	  integer                                 :: iunit
      logical                                 :: connected,exists
      real*8                                  :: tc1,tc2,tc3

       tc1 = UNEST; tc2=UNEST; tc3=UNEST
       if(associated(tc%adv))    tc1 = tc%adv 
       if(associated(tc%disp))   tc2 = tc%disp 
       if(associated(tc%kreact)) tc3 = tc%kreact        
       if (tc1==UNEST.and. tc2==UNEST.and.tc3==UNEST) return
       if(.not.present(fname)) then
           iunit=6
       else
	       inquire(file=fname,exist=exists,opened=connected)
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if  	              
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown',access='append')
           end if
           write(iunit,'(3(3x,a7,x,g13.6))') 'TC_ADV =',tc1,'TC_DISP =',tc2,'TC_REACT =',tc3
       end if
    end subroutine

    subroutine set_characteristic_times_ (tc,val)
         implicit none
		 type(ctimes_cl),            intent(inout) :: tc
		 real*8,                     intent(in)    :: val
         
         if (associated(tc%adv))    tc % adv    = val
         if (associated(tc%disp))   tc % disp   = val
         if (associated(tc%kreact)) tc % kreact = val
     
    end subroutine


!*************************************************************************************
    subroutine update_characteristic_times_ (particle,geo,advection,dispersion,reaction,tc)
         use global_variables, only: UNEST,calcul_time_method,nspedecay,idspedecay
         use particle_class
         use advection_class
         use dispersion_class
         use reaction_class
         use geometry_class
         use array_class
         use to_solve
         use heterogeneity_flags
         use loops_particles, only: nmove
         implicit none

		 type(ctimes_cl),            intent(inout) :: tc 
		 type(particle_cl),          intent(in)    :: particle
		 type(advection_cl),         intent(in)    :: advection
		 type(dispersion_cl),        intent(in)    :: dispersion
		 type(reaction_cl),          intent(in)    :: reaction
		 type(geometry_cl),          intent(in)    :: geo
		 real*8                                    :: vp(3),poro,dd,dmax,Dxx,Dyy,Dzz,kf,dx(3)
		 real*8                                    :: tcadv,tcdisp,tcreact,tcadvmin,tcdispmin,tcreactmin,kfmax
		 integer                                   :: irxn,i
 		 real*8                                    :: ds,dm,dmTh,dmTV,al,ath,atv,DL,DTH,DTV !rpt,vv
		 real*8                                    :: tcmt,tcmtmin,tck,tckmin
		 real*8                                    :: R,k,alpha,beta,alphabeta(2), test
         real*8, save                              :: vv,rpt

         select case (calcul_time_method)

         !.............
         case('CONSTANT_TIME','CONSTANT_DT')
            return

         !.............
         case('OPTIMUM_DT','CONSTANT_MOVE','CONSTANT_CU')

            dx(1) = value_array_ (geo%dx, particle%cell%num(1), 1, 1)
            dx(2) = value_array_ (geo%dy, 1, particle%cell%num(2), 1)
            dx(3) = value_array_ (geo%dz, 1, 1, particle%cell%num(3))
            !...calculate advective characteristic time at particle position
            if (advection%action) then
            !if (.not.vel_homogeneous .OR. .not.poro_homogeneous .OR. .not.sorption_homogeneous .OR. nmove==1) then
		        vp(1) = particle % vel % qpL(1)
		        vp(2) = particle % vel % qpL(2)
		        vp(3) = particle % vel % qpL(3)
	            rpt   = particle % prop % rp
		        poro  = particle % prop % poro
                vp = dabs(vp)/rpt/poro
                vv = dsqrt(vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3))
                ds = dsqrt(vp(1)*dx(1)*dx(1)/vv+vp(2)*dx(2)*dx(2)/vv+vp(3)*dx(3)*dx(3)/vv)
                tcadvmin  = ds/vv
                              
                !...check if update the characteristic times (if smaller than previsous one)
                if (.not.associated(tc%adv)) then
                    allocate (tc%adv)
                    tc = UNEST
                end if  
                if (tc%adv  > tcadvmin .or. tc%adv == UNEST) tc%adv  = tcadvmin
            !end if
            end if
            
            !...
            if (calcul_time_method=='OPTIMUM_DT') then
            
            if (dispersion%action) then
            !if (.not.disp_homogeneous .OR. .not.vel_homogeneous .OR. .not.sorption_homogeneous .OR. nmove==1) then
	            dm   = particle % prop % dm
	            dmTH = particle % prop % dmTH
	            dmTV = particle % prop % dmTV
	            if (advection%action) then
 	                al  = particle % prop % aL  
	                ath = particle % prop % aTH 
	                atv = particle % prop % aTV
	                DL  = ( al  * vv + dm   ) / rpt
                    DTH = ( ath * vv + dmTH ) / rpt
                    DTV = ( atv * vv + dmTV ) / rpt 
	                tcdispmin = ds*ds/max(DL,DTH,DTV)
                else
                    rpt   = particle % prop % rp
                    DL  = ( dm   ) / rpt
                    DTH = ( dmTH ) / rpt
                    DTV = ( dmTV ) / rpt
                    tcdispmin = min(dx(1)*dx(1),dx(2)*dx(2),dx(3)*dx(3))/max(DL,DTH,DTV)                  
                end if
                !...check if update the characteristic times (if smaller than previsous one)
                if (.not.associated(tc%disp)) then
                    allocate (tc%disp)
                    tc = UNEST
                end if
                if (tc%disp  > tcdispmin .or. tc%disp == UNEST) tc%disp  = tcdispmin
            !end if
            end if
            
            if (kineticACTION) then
            !if (.not.kinetic_homogeneous .OR. nmove==1) then
                kfmax = value_array_ (reaction%kinetic%kf(1),particle%cell%num)
                do irxn=2,reaction%kinetic%nreact
                    kf = value_array_ (reaction%kinetic%kf(irxn),particle%cell%num)
                    if (kf > kfmax) kfmax = kf
                end do
                tcreactmin = 1.d0/kfmax
                !...check if update the characteristic times (if smaller than previsous one)
                if (.not.associated(tc%kreact)) then
                    allocate (tc%kreact)
                    tc = UNEST
                end if
                if (tc%kreact  > tcreactmin .or. tc%kreact == UNEST) tc%kreact  = tcreactmin
            end if
            !end if
            
            end if

         !.............
         case('CONSTANT_DAMT')
           !if (.not.mass_trans_homogeneous .OR. .not.poro_homogeneous .OR. .not.sorption_homogeneous .OR. nmove==1) then
            alpha = particle%zone%alpha(particle%specie%num,particle%zone%num)
            beta  = particle%zone%beta(particle%specie%num,particle%zone%num)
            tcmt  = 1/(alpha*(1+beta))
            !...check if update the characteristic times (if smaller than previsous one)
            if (.not.associated(tc%mass_trans)) then
                allocate (tc%mass_trans)
                tc = UNEST
            end if
            if (tc%mass_trans  > tcmt .or. tc%mass_trans == UNEST) tc%mass_trans  = tcmt
           !end if
            
         !.............
         case('CONSTANT_DADECAY')
           !if (.not.decay_homogeneous .OR. .not.sorption_homogeneous .OR. nmove==1) then
            R = 1.d0
            do i=1,nspedecay
                if (idspedecay(i)==particle%specie%num) exit
            end do
            if (particle%zone%num==0) then 
                k = particle%decay%k(i)
                if (sorptionACTION .AND. sorptionTYPE=='LINEAR') R = particle%sorption%linear_sorp%R(particle%specie%num)
            elseif (particle%zone%num==0>0) then 
                k = particle%decay%kim(i,particle%zone%num)
                if (sorptionACTION .AND. sorptionTYPE=='LINEAR') R = particle%sorption%linear_sorp%Rim(particle%specie%num,particle%zone%num)
            end if
            tck = R/k
            !...check if update the characteristic times (if smaller than previsous one)
            if (.not.associated(tc%decay)) then
                allocate (tc%decay)
                tc = UNEST
            end if
            if (tc%decay  > tcmt .or. tc%decay == UNEST) tc%decay  = tck
           !end if
           
         end select

    end subroutine



!*************************************************************************************
!      Maximum velocity and dispersion in the plume
!*************************************************************************************
!    subroutine calculate_characteristic_times_ (tc,plume,geo,advection,dispersion,reaction) 
!         use global_variables, only: UNEST
!         use particle_class
!         use array_class
!         use advection_class
!         use geometry_class
!         use dispersion_class
!         use reaction_class
!         use particle_class
!         use plume_class
!         use list_class
!         use to_solve
!         implicit none
!
!		 type(ctimes_cl),            intent(inout) :: tc 
!		 type(plume_cl),             intent(in)    :: plume
!		 type(geometry_cl),          intent(in)    :: geo
!		 type(advection_cl),         intent(in)    :: advection
!		 type(dispersion_cl),        intent(in)    :: dispersion
!		 type(reaction_cl),          intent(in)    :: reaction
!		 type(particle_cl)                         :: particle
!         type(partID_cl), pointer                  :: part	   
!		 integer                                   :: ip,irxn,ip1,ispe,izone
!		 real*8                                    :: vp(3),rpt,rpt1,poro,dd,dmax,Dxx,Dyy,Dzz,kf,dx(3),kfmax
!		 real*8                                    :: tcadv,tcdisp,tcreact,tcadvmin,tcdispmin,tcreactmin
!		 real*8                                    :: ds,dm,dmTh,dmTV,al,ath,atv,DL,DTH,DTV,vv,D(3)
!         logical, save                             :: FromVeloField =.TRUE.
!
!         
!         if (.not.associated(tc%adv)) then 
!                  allocate (tc%adv)
!                  tc % adv  = UNEST; endif
!         if (.not.associated(tc%disp)) then
!                  allocate (tc%disp)
!                  tc % disp  = UNEST; endif
!         if (.not.associated(tc%kreact)) then
!                  allocate (tc%kreact)
!                  tc % kreact  = UNEST;endif
!
!         if (tc%disp/=UNEST.or.tc%adv/=UNEST.or.tc%kreact/=UNEST) return
!
!
!         if(plume%np <= 0) then
!                     dx(1) = average_array_ (geo%dx)
!                     dx(2) = average_array_ (geo%dy)
!                     dx(3) = average_array_ (geo%dz)		    
!                     rpt = 1.d0
!		             if (kineticACTION) then
!		                 if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
!		                    rpt   = average_array_ (reaction%sorption%LinearSorption%R(1))
!		                    do ispe=2,plume%nspecie
!		                        rpt1 = average_array_ (reaction%sorption%LinearSorption%R(ispe))
!		                        if (rpt < rpt1) rpt = rpt1
!		                    end do
!		                 else
!		                 rpt = 1.d0
!		                 end if
!		             end if      	             
!      	             if (advection%action) then
!                        vp(1) = average_array_ (advection%qx)
!		                vp(2) = average_array_ (advection%qy)
!		                vp(3) = average_array_ (advection%qz)
!		                poro  = average_array_ (advection%poro)
!                        vp = dabs(vp)/rpt/poro
!                        tcadvmin = minval(dx/vp,mask=vp>0.)
!                     end if
!                     if (dispersion%action) then
! 	                    dm   = average_array_ (dispersion%dm)
!	                    dmTH = average_array_ (dispersion%dmTH)
!	                    dmTV = average_array_ (dispersion%dmTV)
!	                    if (advection%action) then 	           
!                           vv = dsqrt(vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3))
! 	                       al  = average_array_ (dispersion%aL)  
!	                       ath = average_array_ (dispersion%aTH) 
!	                       atv = average_array_ (dispersion%aTV)
!	                       D(1) = ( al  * vv + dm   ) / rpt
!                           D(2) = ( ath * vv + dmTH ) / rpt
!                           D(3) = ( atv * vv + dmTV ) / rpt 
!                        else
!                           D(1) = ( dm   ) / rpt
!                           D(2) = ( dmTH ) / rpt
!                           D(3) = ( dmTV ) / rpt
!                        end if    
!	                    tcdispmin = minval(dx*dx/D,mask=D>0.)
!		             end if
!                     if (kineticACTION) then
!                        kfmax = average_array_ (reaction%kinetic%kf(1))
!                        tcreactmin = 1.d0/kfmax
!                        do irxn=2,reaction%kinetic%nreact
!                           kf = average_array_ (reaction%kinetic%kf(irxn))
!                           tcreact = 1.d0/kf
! 		                   if (ip ==1) then
! 		                      tcreactmin = tcreact 
! 		                   elseif (tcreact < tcreactmin ) then
! 		                      tcreactmin = tcreact
! 		                   end if                
!                        end do
!                     end if
!         else
!
!		 ip = 0
!		 
!         species: do ispe  = 1, plume%nspecie
!           zones: do izone = 0, plume%nzone-1
!        
!            part => plume%species(ispe)%zone(izone)%head
! 
!               particles:  do        ! move all particles one time step
!               
!                     if (.not.associated(part)) cycle zones
!
!                     ip = ip + 1
!                     call from_plumeparticle_to_particle_   ( particle, plume%time, part, izone, ispe)  !create particle
!		             call update_cell_location_particle_       ( particle, geo )
!	                 call update_properties_particle_          ( particle, geo, advection, dispersion, reaction   ) !update properties particle poro,rpt,aL,aTH,aTV
!                     call update_velocity_particle_            ( particle, geo, advection, dispersion )             !update velocities qL,qT,qnode,qfaces
!		             call update_dispersion_nodes_particle_    ( particle, geo, dispersion ) 
!                     dx(1) = value_array_ (geo%dx, particle%cell%num(1), 1, 1)
!                     dx(2) = value_array_ (geo%dy, 1, particle%cell%num(2), 1)
!                     dx(3) = value_array_ (geo%dz, 1, 1, particle%cell%num(3))		    
!     	             if (advection%action) then
!                        vp(1) = value_array_ (advection%qx,particle%cell%num)
!		                vp(2) = value_array_ (advection%qy,particle%cell%num)
!		                vp(3) = value_array_ (advection%qz,particle%cell%num)
!		                rpt   = particle % prop % rp
!		                poro  = particle % prop % poro
!                        vp = dabs(vp)/rpt/poro
!                        tcadv = minval(dx/vp,mask=vp>0.)
!                        if (ip ==1) then
!                             tcadvmin  = tcadv
! 		                elseif (tcadv  < tcadvmin  ) then
! 		                     tcadvmin  = tcadv
!                        end if
!                     end if
!                     if (dispersion%action) then
! 	                    dm   = particle % prop % dm
!	                    dmTH = particle % prop % dmTH
!	                    dmTV = particle % prop % dmTV
!	                    if (advection%action) then 	           
!                           vv = dsqrt(vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3))
! 	                       al  = particle % prop % aL  
!	                       ath = particle % prop % aTH 
!	                       atv = particle % prop % aTV
!	                       D(1) = ( al  * vv + dm   ) / rpt
!                           D(2) = ( ath * vv + dmTH ) / rpt
!                           D(3) = ( atv * vv + dmTV ) / rpt 
!                           tcdisp = minval(dx*dx/D,mask=D>0.)
!                        else
!                           rpt  = particle % prop % rp
!                           D(1) = ( dm   ) / rpt
!                           D(2) = ( dmTH ) / rpt
!                           D(3) = ( dmTV ) / rpt
!                           tcdisp = minval(dx*dx/D,mask=D>0.)                  
!                        end if    
!                        if (ip==1) then
!                            tcdispmin = tcdisp
!                        elseif(tcdisp < tcdispmin ) then
!                            tcdispmin = tcdisp
!                        end if
!		             end if
!                     if (kineticACTION) then
!                        kfmax = value_array_ (reaction%kinetic%kf(1),particle%cell%num)
!                        tcreactmin = 1.d0/kfmax
!                        do irxn=2,reaction%kinetic%nreact
!                           kf = value_array_ (reaction%kinetic%kf(irxn),particle%cell%num)
!                           tcreact = 1.d0/kf
! 		                   if (ip ==1) then
! 		                      tcreactmin = tcreact 
! 		                   elseif (tcreact < tcreactmin ) then
! 		                      tcreactmin = tcreact
! 		                   end if                
!                        end do
!                     end if
!		            
!		            part => part%next	 
!			 
!		     end do particles
!		 end do zones
!		 end do species
!        
!         end if
!!        end if
!
!         if (advection%action) then
!              if (.not.associated(tc%adv)) allocate (tc%adv)
!              tc % adv  = tcadvmin
!         end if
!         if (dispersion%action) then
!              if (.not.associated(tc%disp)) allocate (tc%disp)
!              tc % disp  = tcdispmin
!         end if
!         if (kineticACTION) then
!              if (.not.associated(tc%kreact)) allocate (tc%kreact)
!              tc % kreact  = tcreactmin
!         end if
!
!    end subroutine


!*************************************************************************************
!      Maximum velocity and dispersion in the plume
!*************************************************************************************
!    subroutine calculate_adv_characteristic_time_ (tc,plume,geo,advection,dispersion,reaction) 
!         use global_variables, only: UNEST
!         use particle_class
!         use array_class
!         use advection_class
!         use geometry_class
!         use dispersion_class
!         use reaction_class
!         use particle_class
!         use plume_class
!         use list_class
!         use to_solve
!         implicit none
!
!		 type(ctimes_cl),            intent(inout) :: tc 
!		 type(plume_cl),             intent(in)    :: plume
!		 type(geometry_cl),          intent(in)    :: geo
!		 type(advection_cl),         intent(in)    :: advection
!		 type(dispersion_cl),        intent(in)    :: dispersion
!		 type(reaction_cl),          intent(in)    :: reaction
!		 type(particle_cl)                         :: particle
!         type(partID_cl), pointer                  :: part	   
!		 integer                                   :: ip,irxn,ip1,ispe,izone
!		 real*8                                    :: vp(3),rpt,rpt1,poro,dd,dmax,Dxx,Dyy,Dzz,kf,dx(3),kfmax
!		 real*8                                    :: tcadv,tcdisp,tcreact,tcadvmin,tcdispmin,tcreactmin
!		 real*8                                    :: ds,dm,dmTh,dmTV,al,ath,atv,DL,DTH,DTV,vv
!         logical, save                             :: FromVeloField =.TRUE.
!
!         
!         if (.not.associated(tc%adv)) then 
!                  allocate (tc%adv)
!                  tc % adv  = UNEST; endif
!         if (.not.associated(tc%disp)) then
!                  allocate (tc%disp)
!                  tc % disp  = UNEST; endif
!         if (.not.associated(tc%kreact)) then
!                  allocate (tc%kreact)
!                  tc % kreact  = UNEST;endif
!
!         if (tc%adv/=UNEST) return 
!
!		 ip = 0
!
!         species: do ispe  = 1, plume%nspecie
!           zones: do izone = 0, plume%nzone-1
!        
!            part => plume%species(ispe)%zone(izone)%head
! 
!               particles:  do        ! move all particles one time step
!               
!                     if (.not.associated(part)) cycle zones
!
!                     ip = ip + 1
!                     call from_plumeparticle_to_particle_   ( particle, plume%time, part, izone, ispe)  !create particle
!		             call update_cell_location_particle_       ( particle, geo )
!	                 call update_properties_particle_          ( particle, geo, advection, dispersion, reaction   ) !update properties particle poro,rpt,aL,aTH,aTV
!                     call update_velocity_particle_            ( particle, geo, advection, dispersion )             !update velocities qL,qT,qnode,qfaces
!		             call update_dispersion_nodes_particle_    ( particle, geo, dispersion ) 
!                     dx(1) = value_array_ (geo%dx, particle%cell%num(1), 1, 1)
!                     dx(2) = value_array_ (geo%dy, 1, particle%cell%num(2), 1)
!                     dx(3) = value_array_ (geo%dz, 1, 1, particle%cell%num(3))		    
!     	             if (advection%action) then
!                        vp(1) = value_array_ (advection%qx,particle%cell%num)
!		                vp(2) = value_array_ (advection%qy,particle%cell%num)
!		                vp(3) = value_array_ (advection%qz,particle%cell%num)
!		                rpt   = particle % prop % rp
!		                poro  = particle % prop % poro
!                        vp = dabs(vp)/rpt/poro
!                        vv = dsqrt(vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3))
!                        ds = dsqrt(vp(1)*dx(1)*dx(1)/vv+vp(2)*dx(2)*dx(2)/vv+vp(3)*dx(3)*dx(3)/vv)
!                        tcadv = ds/vv
!                        if (ip ==1) then
!                             tcadvmin  = tcadv
! 		                elseif (tcadv  < tcadvmin  ) then
! 		                     tcadvmin  = tcadv
!                        end if
!                     end if
!                                 
!		            part => part%next	 
!			 
!		     end do particles
!		end do zones
!		end do species
!
!        if (advection%action) then
!             if (.not.associated(tc%adv)) allocate (tc%adv)
!             tc % adv  = tcadvmin
!        end if
!
!    end subroutine


!*************************************************************************************
!      Maximum velocity and dispersion in the plume
!*************************************************************************************
!    subroutine calculate_rx_characteristic_times_ (tc,plume,geo,reaction,mass_trans) 
!         use global_variables, only: UNEST,calcul_time_method,nspedecay
!         use particle_class
!         use array_class
!         use geometry_class
!         use reaction_class
!         use mass_trans_class
!         use particle_class
!         use plume_class
!         use list_class
!         use to_solve
!         implicit none
!
!		 type(ctimes_cl),       intent(inout) :: tc
!		 type(plume_cl),        intent(in)    :: plume
!		 type(geometry_cl),     intent(in)    :: geo
!		 type(reaction_cl),     intent(in)    :: reaction
!		 type(mass_trans_cl),   intent(in)    :: mass_trans
!		 type(particle_cl)                    :: particle
!         type(partID_cl),       pointer       :: part
!		 integer                              :: ip,ispe,izone,i
!		 real*8                               :: tcmt,tcmtmin,tck,tckmin
!		 real*8                               :: R,k,alpha,beta,alphabeta(2)
!
!         if (.not.associated(tc%decay)) then 
!            allocate (tc%decay)
!            tc % decay  = UNEST ;endif
!         if (.not.associated(tc%mass_trans)) then
!            allocate (tc%mass_trans)
!            tc % mass_trans  = UNEST ;endif
!
!         if (tc%decay/=UNEST .OR. tc%mass_trans/=UNEST) return
!
!         ip = 0
!         species: do ispe  = 1, plume%nspecie
!           zones: do izone = 0, plume%nzone-1
!
!                part => plume%species(ispe)%zone(izone)%head
!
!                particles:  do        ! move all particles one time step
!
!                     if (.not.associated(part)) cycle zones
!                     ip = ip + 1
!                     call from_plumeparticle_to_particle_( particle, plume%time, part, izone, ispe )  !create particle
!		             call update_cell_location_particle_ ( particle, geo )
!	                 call update_sorption_               ( particle, geo, reaction )    !update retardation
!                     call update_decay_                  ( particle, geo, reaction )    !update decay parameters
!                     call update_mass_trans_             ( particle, geo, mass_trans )  !update mass transfer parameters
!
!                     select case (calcul_time_method)
!
!                        case ('CONSTANT_DAMT')
!                            if (.not.mass_transACTION) stop 'Cannot use CONSTANT_DAMT'
!                            alpha   = particle%zone%alpha(ispe,izone)
!                            beta    = particle%zone%beta(ispe,izone)
!                            tcmt = 1/(alpha*(1+beta))
!                            if (ip==1) then
!                                tcmtmin = tcmt
!                            elseif (tcmt < tcmtmin) then
!                                tcmtmin = tcmt
!                            end if
!
!                        case ('CONSTANT_DADECAY')
!                            if (.not.decayACTION) stop 'Cannot use CONSTANT_DADECAY'
!                            R = 1.d0
!                            if (izone==0) then 
!                                k = particle % decay % k(ispe)
!                                if (sorptionACTION .AND. sorptionTYPE=='LINEAR') R = particle % sorption % linear_sorp % R(ispe)
!                            elseif (izone>0) then 
!                                k = particle % decay % kim(ispe,izone)
!                                if (sorptionACTION .AND. sorptionTYPE=='LINEAR') R = particle % sorption % linear_sorp % Rim(ispe,izone)
!                            end if
!                            tck = R/k
!                            if (ip==1) then
!                                tckmin = tck
!                            elseif (tck < tckmin) then
!                                tckmin = tck
!                            end if
!
!                    end select
!                     
!		            part => part%next	 
!			 
!		        end do particles
!         end do zones
!         end do species
!
!        select case (calcul_time_method)
!            case ('CONSTANT_DAMT')
!                if (.not.associated(tc%mass_trans)) allocate (tc%mass_trans)
!                tc % mass_trans = tcmtmin
!            case ('CONSTANT_DADECAY')
!                if (.not.associated(tc%decay))      allocate (tc%decay)
!                tc % decay = tckmin
!        end select
!
!    end subroutine


!*************************************************************************************
    subroutine reinitialize_characteristic_times ( tc )
    use global_variables, only: UNEST
    use heterogeneity_flags
    implicit none

    type(ctimes_cl),       intent(inout) :: tc
        
        ! reinitialize the characteristic times only under heterogeneous conditions
        if ( .not.vel_homogeneous  .OR. .not.poro_homogeneous .OR. .not.sorption_homogeneous )      tc%adv = UNEST
        if ( .not.disp_homogeneous .OR. .not.vel_homogeneous  .OR. .not.sorption_homogeneous )      tc%disp = UNEST
        if ( .not.kinetic_homogeneous )                                                             tc%kreact = UNEST
        if ( .not.decay_homogeneous .OR. .not.sorption_homogeneous )                                tc%decay = UNEST
        if ( .not.mass_trans_homogeneous .OR. .not.poro_homogeneous .OR. .not.sorption_homogeneous) tc%mass_trans = UNEST
        
    end subroutine

!*************************************************************************************
!*************************************************************************************

  end module
