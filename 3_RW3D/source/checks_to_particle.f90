 module checks_to_particle

 public

 contains


!*************************************************************************************    
!       check all particles are within the domain initially
!*************************************************************************************
       subroutine check_plume_inside_domain (plume,geo)
	   use list_class
	   use array_class
	   use plume_class
       use geometry_class

	   implicit none
       
	   type(geometry_cl), intent(in) :: geo
	   type(plume_cl),    intent(in) :: plume
       type(partID_cl),     pointer  :: part	   
	   
	   integer             :: ispe,izone
	   real*8              :: Lx,Ly,Lz	   
       logical, save       :: evaluated =.FALSE.

        if (evaluated) return
        if (plume%np <=0) return

        Lx = get_xmesh_ (geo,geo%nx+1)
		Ly = get_ymesh_ (geo,geo%ny+1)
		Lz = get_zmesh_ (geo,geo%nz+1)
		
        species: do ispe  = 1, plume%nspecie
          zones: do izone = 0, plume%nzone-1
        
            part => plume%species(ispe)%zone(izone)%head
 
               particles:  do        ! move all particles one time step
               
                     if (.not.associated(part)) cycle zones

		             if (part%xp > Lx .or. part%xp < 0.) then
		                write(*,'(x,i8,x,3(g15.6,x))') part%id,part%xp,part%yp,part%zp
		                stop '>> error: particle outside x-domain at t=0'
	                 else if (part%yp > Ly .or. part%yp < 0.) then
		                write(*,'(x,i8,x,3(g15.6,x))') part%id,part%xp,part%yp,part%zp
		                stop '>> error: particle outside y-domain at t=0'
		             else if (part%zp > Lz .or. part%zp < 0.) then
		                write(*,'(x,i8,x,3(g15.6,x))') part%id,part%xp,part%yp,part%zp
		                stop '>> error: particle outside z-domain at t=0'
		             end if
		             
		part => part%next	 
			 
		end do particles
		end do zones
		end do species

        evaluated = .TRUE.

        end subroutine


!**************************************************************************************************
!      check boundary effects
!**************************************************************************************************
       subroutine check_boundary (particle,geo)
       use code_options, only: ibeff => idebug,SaveMemo
	   use global_variables, only: fname_exit,npOut
	   use gslib, only: open_fname
	   use array_class
	   use particle_class
	   use geometry_class

	   implicit none

	   type(particle_cl), intent(inout) :: particle
	   type(geometry_cl), intent(in)    :: geo
	   real*8                           :: Lx,Ly,Lz
	   integer                          :: nx,ny,nz
	   integer                          :: ib(3,2)
	   real*8                           :: xp,yp,zp,tp ! when bouncing correct position
	   real*8                           :: xpold(3)
	   integer                          :: nbounce
	   integer                          :: cellold(3)
       logical                          :: reflection, domextens
	   integer                          :: iunit,ip,ispecie,izone

      Lx = get_xmesh_ (geo,geo%nx+1)
      Ly = get_ymesh_ (geo,geo%ny+1)
      Lz = get_zmesh_ (geo,geo%nz+1)
     
	  ib = geo%ib

      xp = particle%position%xp(1)
      yp = particle%position%xp(2)
	  zp = particle%position%xp(3)
	  
	  ip = particle%id
	  tp = particle%position%tp
	  
	  ispecie = particle%specie%num
	  izone   = particle%zone%num

      
      !     particle sent to other side of the domain (domain extension)
      if (ib(1,1) == 2) then
      if (xp < 0.d0) then
	        xp = Lx + xp
            domextens = .TRUE.
			particle%control%npdomextens = particle%control%npdomextens + 1
	        if (ibeff > 0)write(*,*) 'particle mirroring at ib(1,1)'
	     end if
      end if

      if (ib(2,1) == 2) then
      if (yp < 0.d0) then
	        yp = Ly + yp
            domextens = .TRUE.
			particle%control%npdomextens = particle%control%npdomextens + 1
	        if (ibeff > 0)write(*,*) 'particle mirroring at ib(2,1)'
	     end if
      end if

      if (ib(3,1) == 2) then
      if (zp < 0.d0) then
	        zp = Lz + zp
            domextens = .TRUE.
			particle%control%npdomextens = particle%control%npdomextens + 1
	        if (ibeff > 0)write(*,*) 'particle mirroring at ib(3,1)'
	     end if
      end if

      if (ib(1,2) == 2) then
      if (xp > Lx) then
	        xp = xp - Lx
            domextens = .TRUE.
			particle%control%npdomextens = particle%control%npdomextens + 1
	        if (ibeff > 0)write(*,*) 'particle mirroring at ib(1,2)'
	     end if
      end if
 
      if (ib(2,2) == 2) then
      if (yp > Ly) then
	        yp = yp - Ly
            domextens = .TRUE.
			particle%control%npdomextens = particle%control%npdomextens + 1
	        if (ibeff > 0)write(*,*) 'particle mirroring at ib(2,2)'
	     end if
      end if
      
      if (ib(3,2) == 2) then
      if (zp > Lz) then
	        zp = zp - Lz
            domextens = .TRUE.
			particle%control%npdomextens = particle%control%npdomextens + 1
	        if (ibeff > 0)write(*,*) 'particle mirroring at ib(3,2)'
	     end if
      end if
      
!     particle reflecting at the boundaries

      reflection = .FALSE.

	  if (ib(1,1) == 1) then
         if (xp < 0.d0) then
	        xp = dabs(xp)
            reflection = .TRUE.
			particle%control%npbounce = particle%control%npbounce + 1
	        if (ibeff > 0)write(*,*) 'particle bouncing at ib(1,1)'
	     end if
	  end if

	  if (ib(2,1) == 1) then
         if (yp < 0.d0) then
	        yp = dabs(yp)
			reflection = .TRUE.
            particle%control%npbounce = particle%control%npbounce + 1
	        if (ibeff > 0)write(*,*) 'particle bouncing at ib(2,1)'
	     end if
	  end if 

	  if (ib(3,1) == 1) then
         if (zp < 0.d0) then
	        zp = dabs(zp)
			reflection = .TRUE.
            particle%control%npbounce = particle%control%npbounce + 1
	        if (ibeff > 0)write(*,*) 'particle bouncing at ib(3,1)'
	     end if
	  end if

	  if (ib(1,2) == 1) then
         if (xp > Lx ) then
	        xp = 2.d0*Lx-xp
			reflection = .TRUE.
            particle%control%npbounce = particle%control%npbounce + 1
 	        if (ibeff > 0)write(*,*) 'particle bouncing at ib(1,2)'
         end if
	  end if

	  if (ib(2,2) == 1) then
         if (yp > Ly ) then
	        yp = 2.d0*Ly-yp
			reflection = .TRUE.
            particle%control%npbounce = particle%control%npbounce + 1
 	        if (ibeff > 0)write(*,*) 'particle bouncing at ib(2,2)'
         end if
	  end if

	  if (ib(3,2) == 1) then
         if (zp > Lz ) then
	        zp = 2.d0*Lz-zp
			reflection = .TRUE.
            particle%control%npbounce = particle%control%npbounce + 1
 	        if (ibeff > 0)write(*,*) 'particle bouncing at ib(3,2)'
         end if
	  end if

!......correct particle position when bouncing

      if (reflection.or.domextens) then

         particle % position % xp(1) = xp
         particle % position % xp(2) = yp
	     particle % position % xp(3) = zp

		 call update_cell_location_particle_ ( particle, geo )   !update cell location and local coordinates, save switchcell-flag               

      end if

!      particle out of lateral boundaries

	  if (ib(1,1) == 0) then
 	     if(xp <= (0.d0)) then
           if (ibeff > 0) write(*,*) 'particle out of boundary ib(1,1)'
           if (ibeff > 0) write(*,*)  xp,yp,zp
		   if (ibeff >= 0) then 
		     call open_fname (fname_exit,iunit)
			 write(iunit,1) xp,yp,zp,ispecie,izone,ip,tp
           1 format(3(g15.6,x),3(i4,x),g15.6)
           end if		
		   particle%control%remove = .TRUE.
		   npOut = npOut + 1
		 end if
	  end if

	  if (ib(1,2) == 0) then
         if(xp >= Lx) then
           if (ibeff > 0 ) write(*,*) 'particle out of boundary ib(1,2)'
           if (ibeff > 0 ) write(*,*)  xp,yp,zp
		   if (ibeff >= 0) then 
		     call open_fname (fname_exit,iunit)
			 write(iunit,1) xp,yp,zp,ispecie,izone,ip,tp
		   end if
           particle%control%remove = .TRUE.
           npOut = npOut + 1
         end if
	  end if

	  if (ib(2,1) == 0) then
	     if(yp <= (0.d0)) then
           if (ibeff > 0 ) write(*,*) 'particle out of boundary ib(2,1)'
           if (ibeff > 0 ) write(*,*)  xp,yp,zp
		   if (ibeff >= 0) then
		     call open_fname (fname_exit,iunit)
			 write(iunit,1) xp,yp,zp,ispecie,izone,ip,tp
		   end if
		   particle%control%remove = .TRUE.
		   npOut = npOut + 1
         end if
	  end if

	  if (ib(2,2) == 0) then
         if(yp >= Ly) then
           if (ibeff > 0 ) write(*,*) 'particle out of boundary ib(2,2)'
           if (ibeff > 0 ) write(*,*)  xp,yp,zp
		   if (ibeff >= 0) then
		     call open_fname (fname_exit,iunit)
			 write(iunit,1) xp,yp,zp,ispecie,izone,ip,tp
		   end if
		   particle%control%remove = .TRUE.
		   npOut = npOut + 1
         end if
	  end if

	  if (ib(3,1) == 0) then
	     if(zp <= (0.d0)) then
           if (ibeff > 0 ) write(*,*) 'particle out of boundary ib(3,1)'
           if (ibeff > 0 ) write(*,*)  xp,yp,zp
		   if (ibeff >= 0) then 
		     call open_fname (fname_exit,iunit)
			 write(iunit,1) xp,yp,zp,ispecie,izone,ip,tp
		   end if
		   particle%control%remove = .TRUE.
		   npOut = npOut + 1
         end if
	  end if

	  if (ib(3,2) == 0) then
         if(zp >= Lz) then
           if (ibeff > 0 ) write(*,*) 'particle out of boundary ib(3,2)'
           if (ibeff > 0 ) write(*,*)  xp,yp,zp
		   if (ibeff >= 0) then 
		     call open_fname (fname_exit,iunit)
			 write(iunit,1) xp,yp,zp,ispecie,izone,ip,tp
		   end if
		   particle%control%remove = .TRUE.
		   npOut = npOut + 1
         end if
	  end if

 end subroutine

!***************************************************************************************************
! CHECK PARTICLE CROSSING A WELL
!***************************************************************************************************
  subroutine check_well_arrival (particle,well,reaction,recirculation)
    use gslib, only: root2
	use particle_class
	use well_vect_class
	use breakthru_class
	use reaction_class
	use recirculation_class
	use global_variables
	use code_options, only: idebug,SaveMemo
	use loops_particles, only: ip
	implicit none
    type(well_vect_cl), intent(inout)  :: well
    type(particle_cl),  intent(inout)  :: particle
	type(reaction_cl),  intent(in)     :: reaction
	type(recirculation_cl), intent(in) :: recirculation
	real*8                             :: dt
	real*8                             :: xp1,yp1,xp2,yp2,zp1,zp2,zp,xstep,ystep,tp1, zwelltop, zwellbot
	integer                            :: iwell,nwell,npbtc,spe
    real*8                             :: xwell,ywell,rwell
    real*8                             :: a,b,c,lambda,d1,d2,time,mp,DaI
    integer                            :: id

       if (.not. associated(well%num)) return

	   nwell = nwell_ (well)
	   
	   if ( nwell  <= 0 ) return 

       xp1    = particle%position%xpold(1)
	   yp1    = particle%position%xpold(2)
	   xp2    = particle%position%xp(1)
	   yp2    = particle%position%xp(2)
       zp2    = particle%position%xp(3)
	   mp     = particle%prop%mp
       spe    = particle%specie%num
       id     = particle%id

	   xstep  = xp2-xp1
	   ystep  = yp2-yp1
       
       !if (xstep == 0.d0 .AND. ystep == 0.d0 .AND. imove == 1) then
       !    particle%control%remove = .TRUE.
       !    return
       !end if
       
       tp1 = particle%position%tpold
	   dt  = particle%position%tp - tp1
       
	   do iwell=1,nwell
          
          if (spe>nspecie) cycle
        
		  xwell = well%num(iwell)%xwell
		  ywell = well%num(iwell)%ywell
		  rwell = well%num(iwell)%rwell
          
          zwelltop = well%num(iwell)%zwtop
          zwellbot = well%num(iwell)%zwbot
          
          if (zp2 > zwellbot .and. zp2 < zwelltop) then !particle elevation in screening elevation
          
              d1 = dsqrt((xwell-xp1)**2+(ywell-yp1)**2)
              d2 = dsqrt((xwell-xp2)**2+(ywell-yp2)**2)
	   
	          if( (d1 <= rwell .and. d2 >= rwell ) .or. (d1 >= rwell .and. d2 <= rwell) ) then

	               a = xstep * xstep + ystep * ystep 
	               b = 2.d0 * ( xstep * (xp1-xwell) + ystep * (yp1-ywell) )
	               c = ( (xp1-xwell)**2 ) + ( (yp1-ywell)**2 )-(rwell**2)
	               lambda = root2(a,b,c)
	               time   = tp1 + dt * lambda

               
			       if (.not.associated(well%num(iwell)%zwtop)) then
				  
				      ! For fully penetrating wells
                 		        
			            !if (well%num(iwell)%SaveBTC) call add_particle_to_btc_ (well%num(iwell)%btc(spe),id,time,mp)
                        if (well%num(iwell)%SaveBTC) call add_particle_to_btc_ (well%num(iwell)%btc(spe),id,time,mp,xp2,yp2,zp2)     
	                    if ( well%num(iwell)%out == 1 ) particle%control%remove = .TRUE.
                        particle%control%passwell = iwell
	                    npWell = npWell + 1
			       else 
			      
				       ! For partially penetrating wells  		   
				      zp1 = particle%position%xpold(3)
	                  zp2 = particle%position%xp(3)
			          zp = zp1 + (zp2-zp1)/dt * (time-tp1)
                      if (zp >= well%num(iwell)%zwbot .and.  zp <= well%num(iwell)%zwtop) then        
			            !if (well%num(iwell)%SaveBTC) call add_particle_to_btc_ (well%num(iwell)%btc(spe),id,time,mp)
                        if (well%num(iwell)%SaveBTC) call add_particle_to_btc_ (well%num(iwell)%btc(spe),id,time,mp,xp2,yp2,zp2)
	                    if ( well%num(iwell)%out == 1 ) particle%control%remove = .TRUE.
	                    particle%control%passwell = iwell 
	                    npWell = npWell + 1
			          end if
			       end if			   
			   
              end if
              
          end if
       
	  end do

  end subroutine

 
!***************************************************************************************************
! CHECK PARTICLE IN CELL WHERE A WELL IS LOCATED (TO CONSTRAIN DT)
!*************************************************************************************************** 
  subroutine check_well_cell (particle,well)
    use particle_class
    use well_vect_class
    use global_variables, only: particle_in_well_cell
    implicit none
    type(particle_cl),  intent(inout)   :: particle
    type(well_vect_cl), intent(in)      :: well
    
    logical             :: test(3)
    integer             :: iwll, nwll, partcol, partrow, partlay
   
    !reinitialize the flag  
    particle%control%part_in_wellcell = .FALSE.
    particle%control%wellID_in_partcell = 0
    
    partcol = particle%cell%num(1)
    partrow = particle%cell%num(2)
    partlay = particle%cell%num(3)
    nwll = size(well%num)
    
    do iwll=1,nwll
        test = .FALSE.
        if (partcol == well%num(iwll)%col) test(1) = .TRUE.
        if (partrow == well%num(iwll)%row) test(2) = .TRUE.
        if (partlay <= well%num(iwll)%toplay .AND. partlay >= well%num(iwll)%botlay) test(3) = .TRUE.
        !update the flag if needed
        if (test(1)==.TRUE. .AND. test(2)==.TRUE. .AND. test(3)==.TRUE.) then 
            particle_in_well_cell = .TRUE. 
            particle%control%part_in_wellcell = .TRUE.
            particle%control%wellID_in_partcell = iwll
            exit
        end if
    end do
    
  end subroutine
  
!****************************************************************************************************
! CHECK PARTICLE ARRIVAL AT A PLANE
!****************************************************************************************************
  subroutine check_plane_arrival (particle,plane, reaction)
     use gslib, only: upper_case, rotate_coord, vectorial_product
	 use particle_class
	 use plane_vect_class
	 use breakthru_class
	 use reaction_class
	 use global_variables
	 use code_options, only: idebug,SaveMemo,iwcshotpl,ipldisp
	 use loops_particles, only: ip
	 implicit none
     type(plane_vect_cl), intent(inout) :: plane
	 type(particle_cl),   intent(inout) :: particle
	 type(reaction_cl),   intent(in)    :: reaction
	 real*8                             :: dt
	 integer                            :: nplane,iplane,npbtc,spe
     real*8                             :: xp1,yp1,zp1,tp1,xp2,yp2,zp2,xp,yp,zp,xp0,yp0,zp0
	 real*8                             :: xstep,ystep,zstep,dis,time,mp
	 real*8                             :: A,B,C,D,dist_ptTOpl_1,dist_ptTOpl_2,factor
	 real*8                             :: xpart(3),l(3),m(3),n(3),DaI
	 integer                            :: id

      if (associated(plane%num)) then

      nplane = nplane_ (plane)  
	  
	  if (nplane <= 0) return
	  
      xp0 = particle%position%xpref(1)
	  yp0 = particle%position%xpref(2)
	  zp0 = particle%position%xpref(3)

      xp1 = particle%position%xpold(1)
	  yp1 = particle%position%xpold(2)
	  zp1 = particle%position%xpold(3)

	  xp2 = particle%position%xp(1)
	  yp2 = particle%position%xp(2)
	  zp2 = particle%position%xp(3)

      mp  = particle%prop%mp
      spe = particle%specie%num
      
      id  = particle%id

	  xstep = xp2-xp1
	  ystep = yp2-yp1
	  zstep = zp2-zp1

	  tp1 = particle%position%tpold 

	  dt =  particle%position%tp - tp1   

	  do iplane=1,nplane

           if (spe>nspecie) cycle

           if ( upper_case(plane%num(iplane)%tipo) == 'YY' ) then
	            dis=plane%num(iplane)%xplane
                if( (xp1 <= dis .and. xp2 >= dis) .or. (xp1 >= dis .and. xp2 <= dis) ) then
					 time = tp1 + dabs(dis-xp1)/dabs(xstep) * dt
					 xp   = plane%num(iplane)%xplane
					 yp   = yp1 + dabs(dis-xp1)/dabs(xstep) * ystep
					 zp   = zp1 + dabs(dis-xp1)/dabs(xstep) * zstep
					 if (ipldisp > 0) then
					     xp   = xp - xp0    ! calculate total displacement particle
					     yp   = yp - yp0
					     zp   = zp - zp0
                     end if
					 particle%control%passplane = iplane
                     call add_particle_to_btc_ (plane%num(iplane)%btc(spe),id,time,mp,xp,yp,zp)   
					 if ( plane%num(iplane)%out /= 0 ) then
	                      particle%control%remove = .TRUE.
	                      npPlane = npPlane + 1
	                end if
				end if

	       else if ( upper_case(plane%num(iplane)%tipo) == 'XX' ) then
                dis=plane%num(iplane)%yplane
		        if( (yp1 <= dis .and. yp2 >= dis).or. (yp1 >= dis .and. yp2 <= dis) ) then
					 time = tp1 + dabs(dis-yp1)/dabs(ystep) * dt
                     xp   = xp1 + dabs(dis-yp1)/dabs(ystep) * xstep
					 yp   = plane%num(iplane)%yplane
					 zp   = zp1 + dabs(dis-yp1)/dabs(ystep) * zstep
					 if (ipldisp > 0) then
					 	 xp   = xp - xp0                         !calculate total displacement particle
					     yp   = yp - yp0
					     zp   = zp - zp0
					 end if
					 particle%control%passplane = iplane
                     call add_particle_to_btc_ (plane%num(iplane)%btc(spe),id,time,mp,xp,yp,zp) 
					 if ( plane%num(iplane)%out /= 0 )  then
	                      particle%control%remove = .TRUE.
	                      npPlane = npPlane + 1
	                end if
                end if

	       else if ( upper_case(plane%num(iplane)%tipo) == 'ZZ' ) then
                dis=plane%num(iplane)%zplane
		        if( (zp1 <= dis .and. zp2 >= dis).or. (zp1 >= dis .and. zp2 <= dis) ) then
					 time = tp1 + dabs(dis-zp1)/dabs(zstep) * dt
                     xp   = xp1 + dabs(dis-zp1)/dabs(zstep) * xstep
					 yp   = yp1 + dabs(dis-zp1)/dabs(zstep) * ystep
					 zp   = plane%num(iplane)%zplane
					 if (ipldisp > 0) then
					 	 xp   = xp - xp0                         !calculate total displacement particle
					     yp   = yp - yp0
					     zp   = zp - zp0
					 end if
					 particle%control%passplane = iplane
                     call add_particle_to_btc_ (plane%num(iplane)%btc(spe),id,time,mp,xp,yp,zp) 
					 if ( plane%num(iplane)%out /= 0 )  then
	                      particle%control%remove = .TRUE.
	                      npPlane = npPlane + 1
	                end if
                end if
                
                
           else if ( upper_case(plane%num(iplane)%tipo) == 'GL' ) then
		        A = plane%num(iplane)%A 
				B = plane%num(iplane)%B
				C = plane%num(iplane)%C
				D = plane%num(iplane)%D
		        dist_ptTOpl_1 = (A*xp1 + B*yp1 + C*zp1 + D) / dsqrt(A*A+B*B+C*C) !min distance from point 1 to plane  
				dist_ptTOpl_2 = (A*xp2 + B*yp2 + C*zp2 + D) / dsqrt(A*A+B*B+C*C) !min distance from point 1 to plane
				if ( dist_ptTOpl_1/dist_ptTOpl_2 < 0.) then  ! particle crosses the plane
					 npbtc  = npbtc+1
					 factor = dabs(A*xp1 + B*yp1 + C*zp1 + D) / dabs(A*xstep + B*ystep + C*zstep)			 
					 time = tp1 + factor * dt                    
					 xp   = xp1 + factor * xstep
					 yp   = yp1 + factor * ystep
					 zp   = zp1 + factor * zstep
					 if (ipldisp > 0) then				 
					    xpart(1)   = xp - xp0                 !......calculate total displacement particle
					    xpart(2)   = yp - yp0
					    xpart(3)   = zp - zp0                 !......rotate coordenates
					 else
					    xpart(1)   = xp
						xpart(2)   = yp
						xpart(3)   = zp
					 end if
					 !...............................calculate directional cosines
                         l(1) = A/dsqrt(A*A+B*B+C*C)
		                 l(2) = B/dsqrt(A*A+B*B+C*C)
		                 l(3) = C/dsqrt(A*A+B*B+C*C)
						 if ( l(1) /= 0.d0 .or. l(2) /= 0.d0) then
						      m(1) = -l(2)
						      m(2) = +l(1)
						      m(3) =  0.d0  !horizontal plane
						      n = vectorial_product (l,m)
						 else if ( l(1) == 0.d0 .and. l(2) == 0.d0 ) then
						      m(1) =  l(3)
							  m(2) =  0.d0
							  m(3) = -l(1)
							  n = vectorial_product (l,m)
                         end if
					!.................................rotates coordenates increment displacement of a particle
					 if(ipldisp>0) xpart = rotate_coord (xpart,l,m,n) 
					     
					 particle%control%passplane = iplane
					 call add_particle_to_btc_ (plane%num(iplane)%btc(spe),id,time,mp,xpart(1),xpart(2),xpart(3))  
					 if ( plane%num(iplane)%out /= 0 )  then
	                      particle%control%remove = .TRUE.
	                      npPlane = npPlane + 1
	                end if
			    end if

		   end if

	  end do

	  end if

  end subroutine


!****************************************************************************************************
! CHECK PARTICLE LOCATED IN A BLOCK
!****************************************************************************************************
  subroutine check_block_arrival (particle,obs_block,reaction)
     !use gslib, only: upper_case, rotate_coord, vectorial_product
	 use particle_class
	 use block_vect_class
	 use breakthru_class
	 use reaction_class
	 !use global_variables
	 !use code_options, only: idebug,SaveMemo,iwcshotpl,ipldisp
	 !use loops_particles, only: ip
     implicit none
     type(block_vect_cl), intent(inout) :: obs_block
	 type(particle_cl),   intent(inout) :: particle
	 type(reaction_cl),   intent(in)    :: reaction
     integer                            :: nblock,spe,id,iblock
     real*8                             :: xp0,yp0,zp0,xp1,yp1,zp1,xp2,yp2,zp2,mp
     real*8                             :: xb1,xb2,yb1,yb2,zb1,zb2,tp1
     
     if (associated(obs_block%num)) then

     nblock = nblock_ (obs_block)  
	  
	 if (nblock <= 0) return
     
     xp0 = particle%position%xpref(1)
	 yp0 = particle%position%xpref(2)
	 zp0 = particle%position%xpref(3)
     
     xp1 = particle%position%xpold(1)
     yp1 = particle%position%xpold(2)
     zp1 = particle%position%xpold(3)

	 xp2 = particle%position%xp(1)
	 yp2 = particle%position%xp(2)
	 zp2 = particle%position%xp(3)

     mp  = particle%prop%mp
     spe = particle%specie%num
      
     id  = particle%id
     
     tp1 = particle%position%tpold 
     !dt =  particle%position%tp - tp1   
     
     do iblock=1,nblock
         
         xb1 = obs_block%num(iblock)%x1
         xb2 = obs_block%num(iblock)%x2
         yb1 = obs_block%num(iblock)%y1
         yb2 = obs_block%num(iblock)%y2
         zb1 = obs_block%num(iblock)%z1
         zb2 = obs_block%num(iblock)%z2
         
         !check if particle enters a block
         if ( xp2 > xb1 .and. yp2 > yb1 .and. zp2 > zb1 .and. xp2 < xb2 .and. yp2 < yb2 .and. zp2 < zb2 ) then !particle inside the block
             particle%control%inblock = iblock
             call add_particle_to_btc_ (obs_block%num(iblock)%btc(spe),id,tp1,mp,xp2,yp2,zp2)  !no interpolation of time and particle location
             obs_block%num(iblock)%mass(spe) = obs_block%num(iblock)%mass(spe) + mp
         end if
         
         !check if particle leaves a block
         if ( xp1 > xb1 .and. yp1 > yb1 .and. zp1 > zb1 .and. xp1 < xb2 .and. yp1 < yb2 .and. zp1 < zb2 ) then !particle was in the block
            if (xp2 < xb1 .or. yp2 < yb1 .or. zp2 < zb1 .or. xp2 > xb2 .or. yp2 > yb2 .or. zp2 > zb2) then
                particle%control%inblock = 0
                obs_block%num(iblock)%mass(spe) = obs_block%num(iblock)%mass(spe) - mp
            end if
         end if
         
     end do
     
     end if
     
  end subroutine
  
  
!***********************************************************************************************
! CHECK INACTIVE CELL (IF INACTIVE MAKE PARTICLE TO BOUNCE)
!***********************************************************************************************
  subroutine check_inactive_cell ( particle, geo )
     use loops_particles
	 use gslib, only: open_fname 
	 use global_variables, only: nmaxref,fdbg
	 use gslib, only: sortem
	 use particle_class
	 use geometry_class
	 implicit none
	   type(geometry_cl), intent(in)    :: geo
	   type(particle_cl), intent(inout) :: particle
	   integer                          :: ncell,icell,jcell,kcell,i,iselect
	   integer                          :: ip1,jp1,kp1,ip2,jp2,kp2,iref,unit
	   integer                          :: imin,jmin,kmin,imax,jmax,kmax
	   real*8                           :: dd,cosdir(3),xpnew(3),xpold(3),xint(3)
	   real*8                           :: x1,x2,y1,y2,z1,z2
	   real*8                           :: lambda(6),dx(3),lbdamin
	   logical                          :: Reflection
      
      if (.not. particle%control%switchcell ) return
	  if (.not. associated(geo%InactCell) )   return

!..... start big loop to account for more than one reflection at a time 

      loop_reflect: do iref=1,nmaxref

      xpnew = particle%position%xp
	  xpold = particle%position%xpold
     
	  ip1 = particle % cell % numold(1)
	  jp1 = particle % cell % numold(2)
      kp1 = particle % cell % numold(3)
     
	  ip2 = particle % cell % num(1)
	  jp2 = particle % cell % num(2)
      kp2 = particle % cell % num(3)

      imin = min(ip1,ip2)
	  jmin = min(jp1,jp2)
	  kmin = min(kp1,kp2)

	  if (imin <= 0) imin = 1
	  if (jmin <= 0) jmin = 1
	  if (kmin <= 0) kmin = 1

      imax = max(ip1,ip2)
	  jmax = max(jp1,jp2)
	  kmax = max(kp1,kp2)    

	  if (imax > geo%nx) imax = geo%nx
	  if (jmax > geo%ny) jmax = geo%ny
	  if (kmax > geo%nz) kmax = geo%nz     

      lambda  = -1.d0
	  lbdamin = 1.d+99
      Reflection = .FALSE.

	  do kcell = kmin, kmax
	  do jcell = jmin, jmax
	  do icell = imin, imax
		 
		 if ( icell == ip1 .and. jcell == jp1 .and. kcell == kp1) cycle 
		     
		 if ( geo % InactCell(icell,jcell,kcell) == 0 ) then

	          dd = dsqrt( ( xpnew(1) - xpold(1) ) **2 +  &
	                      ( xpnew(2) - xpold(2) ) **2 +  &
			              ( xpnew(3) - xpold(3) ) **2 )
			   
	          cosdir(1) =  ( xpnew(1) - xpold(1) ) / dd
	          cosdir(2) =  ( xpnew(2) - xpold(2) ) / dd
	          cosdir(3) =  ( xpnew(3) - xpold(3) ) / dd

	          ! calculate lambda's

              x1 = get_xmesh_ ( geo, icell ) 
	          y1 = get_ymesh_ ( geo, jcell )
	          z1 = get_zmesh_ ( geo, kcell )
	          x2 = get_xmesh_ ( geo, icell + 1 )
	          y2 = get_ymesh_ ( geo, jcell + 1 )
	          z2 = get_zmesh_ ( geo, kcell + 1 )

	          if (cosdir(1) /= 0.d0) lambda(1) = ( x1 - xpold(1) ) / cosdir(1)
	          if (cosdir(1) /= 0.d0) lambda(2) = ( x2 - xpold(1) ) / cosdir(1)
	
	          if (cosdir(2) /= 0.d0) lambda(3) = ( y1 - xpold(2) ) / cosdir(2)
	          if (cosdir(2) /= 0.d0) lambda(4) = ( y2 - xpold(2) ) / cosdir(2) 

	          if (cosdir(3) /= 0.d0) lambda(5) = ( z1 - xpold(3) ) / cosdir(3)
	          if (cosdir(3) /= 0.d0) lambda(6) = ( z2 - xpold(3) ) / cosdir(3) 
                 
              do i=1,6             

				 xint = xpold + cosdir * lambda(i)
                
				 if ( xint(1) >= x1 .and. xint(1) <= x2 .and. &
	                  xint(2) >= y1 .and. xint(2) <= y2 .and. &
		              xint(3) >= z1 .and. xint(3) <= z2 .and. &
					  lambda(i) < lbdamin .and. lambda(i) > 0.d0 ) then
		      		      Reflection = .TRUE.
						  lbdamin = lambda(i)	 
					      iselect = i
	             end if
              end do
        
		 end if

	  end do
	  end do
	  end do

	  if (Reflection) then
	  		
				 xint = xpold + cosdir * lbdamin

                      select case (iselect)

		              case(1)
			             particle%position%xp(1) = xpnew(1)-2.d0*(xpnew(1)-x1)
		              case(2)
			             particle%position%xp(1) = xpnew(1)-2.d0*(xpnew(1)-x2)
		              case(3)
			             particle%position%xp(2) = xpnew(2)-2.d0*(xpnew(2)-y1)
		              case(4)
			             particle%position%xp(2) = xpnew(2)-2.d0*(xpnew(2)-y2)
		              case(5)
			             particle%position%xp(3) = xpnew(3)-2.d0*(xpnew(3)-z1)
		              case(6)
			             particle%position%xp(3) = xpnew(3)-2.d0*(xpnew(3)-z2)
		              
					  end select                

				 call update_cell_location_particle_  ( particle, geo )   !update cell location and local coordinates, save switchcell-flag               
			
				 particle % position % xpold = xpold

				 particle % cell % numold(1) = ip1
				 particle % cell % numold(2) = jp1
				 particle % cell % numold(3) = kp1

      else if (.not. reflection ) then

         exit loop_reflect

      end if

  end do loop_reflect    

!.... in case particle keeps reflecting then remove the particle 

  if (iref == nmaxref) then

      particle % control % remove = .TRUE.
      call open_fname (fdbg,unit)
	  write(unit,*)
	  write(unit,*) ' **** WARNING: particle removed because too many &
	                           reflections in one time step' 

	  write(*   ,*) ' **** WARNING: particle removed because too many &
	                           reflections in one time step' 
      write(unit,*) ' PARTICLE:',ip, ' INJECTION: ',kinj,' NMOVE: ',nmove 
	  call print_position_particle_ (particle,fdbg)
	  write(unit,*) 

  end if

  end subroutine



 end module