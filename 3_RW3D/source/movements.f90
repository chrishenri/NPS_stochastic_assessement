 module particle_movements

 public

 contains


!*************************************************************************************
!      Move one step
!*************************************************************************************
    function move_part_by_advection_dispersion ( advection, dispersion, particle, geo, dt ) result ( dxp )
       use particle_class
	   use advection_class
	   use dispersion_class
	   use geometry_class
	   use heterogeneity_flags
	   use loops_particles, only: nmove,ip
       use to_solve,        only: decayACTION, decayTYPE
	   implicit none

	   type(advection_cl), intent (in)      :: advection
	   type(dispersion_cl),intent (in)      :: dispersion
	   type(particle_cl),  intent (in)      :: particle
	   type(geometry_cl),  intent (in)      :: geo
	   real*8,             intent (in)      :: dt

	   real*8              :: dxp(3)                            !particle displacement in one time step
	   real*8              :: dxp_AD(3),dxp_GD(3),dxp_BR(3)     !particle displacements due to ADVECTION, GRADIENT, DISPERSION
	   
	   if (decayACTION .AND. decayTYPE == 'SERIAL_REACTION_MOMENTS') then
	        dxp_AD = move_serial_mom_adv  ( advection, dispersion, particle, geo, dt )
            dxp_BR = move_serial_mom_bro  ( advection, geo, dispersion, particle, dt )

            dxp = dxp_AD + dxp_BR

       else
       	    dxp_AD = move_one_step_advective ( advection,  particle, geo, dt )    !calculate advective step
            dxp_GD = move_one_step_Gradient  ( dispersion, particle, geo, dt )    !calculate gradient step
            dxp_BR = move_one_step_Brownian  ( dispersion, particle, dt )         !calculate dispersion step

            dxp    = dxp_AD + dxp_GD + dxp_BR                   !increment particle position
       end if

	end function

!*************************************************************************************
!      1 Pure Advective Movement of 1 particle (without corrected velocity)
!*************************************************************************************
       function move_one_step_advective (advection,particle,geo,dt) result (xstep)
	   use particle_class
	   use advection_class
	   use geometry_class
	   use array_class
	   use global_variables, only: method_advection
	   use to_solve
	   implicit none

       type(advection_cl), intent (in) :: advection 
	   type(particle_cl),  intent (in) :: particle
	   type(geometry_cl),  intent (in) :: geo 
	   real*8,             intent (in) :: dt
	   real*8                          :: xstep(3)
	   real*8                          :: qx,qy,qz,poro,rpt,btot
       real*8                          :: dx(3),vface(3,2),vp(3),A(3)
	   integer                         :: i
	         
       xstep = 0.d0

       if (.not. advection % action) return !only moves if advection package is on

       if (particle % zone % num > 0) return !only moves particles in the Mobile Phase
       
       if (particle % phase % num > 0) return  !only moves aqueous phase

       if ( dt <= 0.d0 ) return 

       select case (method_advection)

       case ('EXPONENTIAL')

         dx(1) = value_array_ (geo%dx, particle%cell%num(1), 1, 1)
         dx(2) = value_array_ (geo%dy, 1, particle%cell%num(2), 1)
         dx(3) = value_array_ (geo%dz, 1, 1, particle%cell%num(3))

	     vface(1,1) = particle % cell % qxface(1) / particle % prop % poro
	     vface(1,2) = particle % cell % qxface(2) / particle % prop % poro
	     vface(2,1) = particle % cell % qyface(1) / particle % prop % poro
	     vface(2,2) = particle % cell % qyface(2) / particle % prop % poro
	     vface(3,1) = particle % cell % qzface(1) / particle % prop % poro
	     vface(3,2) = particle % cell % qzface(2) / particle % prop % poro    

	     vp(1) = particle % vel % qpL(1) / particle % prop % poro
	     vp(2) = particle % vel % qpL(2) / particle % prop % poro
	     vp(3) = particle % vel % qpL(3) / particle % prop % poro
	               
	     A(1) = ( vface(1,2) - vface(1,1) ) / dx(1)
	     A(2) = ( vface(2,2) - vface(2,1) ) / dx(2)
	     A(3) = ( vface(3,2) - vface(3,1) ) / dx(3)
       
         if (sorptionACTION .AND. sorptionTYPE == 'LINEAR') then
            if (decayACTION .AND. decayTYPE == 'SERIAL')    then
                rpt = Rap(particle)
            else
                rpt = particle % prop % rp
            end if
         else
            rpt = 1.d0
         end if

	     do i=1,3
            if ( A(i) /= 0.d0) then
			      xstep(i) = vp(i)/A(i)/rpt*(dexp(A(i)*dt)-1.d0)
			else
			      xstep(i) = vp(i)/rpt*dt
		    endif			      
         end do

	   case ('EULERIAN')

         qx   = particle % vel  % qpL(1) 
         qy   = particle % vel  % qpL(2) 
	     qz   = particle % vel  % qpL(3)
	     poro = particle % prop % poro
	     
         if (sorptionACTION .AND. sorptionTYPE == 'LINEAR') then
            if (decayACTION .AND. decayTYPE == 'SERIAL')    then
                rpt = Rap(particle)
            else
                rpt = particle % prop % rp
            end if
         else
            rpt = 1.d0
         end if

	     xstep(1) = qx / poro / rpt * dt
	     xstep(2) = qy / poro / rpt * dt
	     xstep(3) = qz / poro / rpt * dt

	   case default

         qx   = particle % vel  % qpL(1) 
         qy   = particle % vel  % qpL(2) 
	     qz   = particle % vel  % qpL(3)
	     poro = particle % prop % poro

         if (sorptionACTION .AND. sorptionTYPE == 'LINEAR') then
            if (decayACTION .AND. decayTYPE == 'SERIAL')    then
                rpt = Rap(particle)
            else
                rpt = particle % prop % rp
            end if
         else
            rpt = 1.d0
         end if

	     xstep(1) = qx / poro / rpt * dt
	     xstep(2) = qy / poro / rpt * dt
	     xstep(3) = qz / poro / rpt * dt

       end select


	   end function



!*************************************************************************************
!      1 Brownian Movement of 1 particle
!*************************************************************************************
   function move_one_step_Brownian (dispersion,particle,dt) result (move)
       use dispersion_class
	   use particle_class
	   use gslib, only: random_normal
	   use to_solve
	   implicit none
	   
       type(dispersion_cl), intent(in) :: dispersion 
       type(particle_cl),   intent(in) :: particle
	   real*8,              intent(in) :: dt
	   real*8                          :: move(3)
       real*8 :: poro,qx,qy,qz,al,ath,atv,dm,rpt,dmTH,dmTV
	   real*8 :: vx,vy,vz,v,vh,z1x,z1y,z1z,z2x,z2y,z3x,z3y,z3z,z1,z2,z3
	   real*8 :: v2,vh2 

       move = 0.d0

       if ( .not.dispersion%action ) return

       if (particle%zone%num > 0) return !only moves particles in the Mobile Phase
       
       if (particle % phase % num > 0) return  !only moves aqueous phase

       if ( dt <= 0.d0 ) return 

       poro  = particle % prop % poro 
	   qx   = particle % vel  % qpT(1)
	   qy   = particle % vel  % qpT(2)
	   qz   = particle % vel  % qpT(3)

	   dm     = particle % prop % dm
	   dmTH   = particle % prop % dmTH
	   dmTV   = particle % prop % dmTV
	   	   
	   if (sorptionACTION .AND. sorptionTYPE == 'LINEAR') then
          if (decayACTION .AND. decayTYPE == 'SERIAL')    then
              rpt = Rap(particle)
          else
              rpt = particle % prop % rp
          end if
       else
          rpt = 1.d0
	   end if

       vx = qx / poro
	   vy = qy / poro
	   vz = qz / poro

       v2 = vx * vx + vy * vy + vz * vz

       if (v2 == 0.d0) then; v = 0.d0
	                   else; v = dsqrt(vx * vx + vy * vy + vz * vz); endif

       vh2 = vx * vx + vy * vy

       if (vh2 == 0.d0) then; vh = 1.d-21
	                    else; vh = dsqrt( vx * vx + vy * vy); endif 


!------- special cases:
       
	   if ( v == 0.d0 .and. dm /=0.d0 ) then   !if no velocity allow molecular diffusion

             z1 = random_normal ()
             z2 = random_normal ()
             z3 = random_normal ()

             move(1) = z1 * dsqrt(2.d0*(dm)/rpt*dt)     
             move(2) = z2 * dsqrt(2.d0*(dmTH)/rpt*dt)
             move(3) = z3 * dsqrt(2.d0*(dmTV)/rpt*dt)       
	         return
	   
	   else if ( particle%zone%num > 0 .and. dm /=0.d0 ) then !if immobile region allow molecular diffusion

             z1 = random_normal ()
             z2 = random_normal ()
             z3 = random_normal ()

             move(1) = z1 * dsqrt(2.d0*(dm)/rpt*dt)     
             move(2) = z2 * dsqrt(2.d0*(dmTH)/rpt*dt)
             move(3) = z3 * dsqrt(2.d0*(dmTV)/rpt*dt)       
	         return

	   else if ( v  == 0.d0 .and. dm ==0.d0 ) then
	        
			 return
	   
	   else if ( particle%zone%num > 0 .and. dm ==0.d0 ) then
	   
	         return

	   end if

!--------- General Case:
!
	   al   = particle % prop % aL
	   ath  = particle % prop % aTH
	   atv  = particle % prop % aTV

!------ generation random numbers: std. normal dist. n(0,1)

        z1 = random_normal ()
        z2 = random_normal ()
        z3 = random_normal ()

!------- calculates projections rvs=[z1,z2,z3] to coordinates x,y,z

        z1x = vx/v
        z1y = vy/v
        z1z = vz/v

        z2x =  vy*vh/v/v +  vz*vz*vy/v/v/vh
        z2y = -vx*vh/v/v - vz*vz*vx/v/v/vh
!       z2z = 0.

        z3x = -vz*vx/v/vh
        z3y = -vz*vy/v/vh
        z3z =  vh/v         


!------- calculates incremental movement

        move(1) = z1 * dsqrt(2.d0*(al *v+dm)/rpt*dt)   * z1x +			&
                  z2 * dsqrt(2.d0*(ath*v+dmTH)/rpt*dt) * z2x +			&
                  z3 * dsqrt(2.d0*(atv*v+dmTV)/rpt*dt) * z3x

        move(2) = z1 * dsqrt(2.d0*(al *v+dm)/rpt*dt)   * z1y +			&
                  z2 * dsqrt(2.d0*(ath*v+dmTH)/rpt*dt) * z2y +			&
                  z3 * dsqrt(2.d0*(atv*v+dmTV)/rpt*dt) * z3y

        move(3) = z1 * dsqrt(2.d0*(al *v+dm)/rpt*dt)   * z1z +			&
                  z3 * dsqrt(2.d0*(atv*v+dmTV)/rpt*dt) * z3z      


   end function
   
!*******************************************************************************************************************
!      1 Advective Movement of 1 particle for corrected velocity 
!            corrected velocity = divergence(Dispersion*poro)
!*******************************************************************************************************************
   function move_one_step_Gradient (dispersion,particle,geo,dt) result(move)
        use dispersion_class
		use array_class
		use particle_class
		use geometry_class
		use to_solve
		implicit none 

        type(dispersion_cl), intent(in) :: dispersion
        type(particle_cl),   intent(in) :: particle
		type(geometry_cl),   intent(in) :: geo
		real*8,              intent(in) :: dt
  		real*8                          :: move(3)
	    real*8                          :: Dxxnode(2,2,2),Dyynode(2,2,2),Dzznode(2,2,2)
		real*8                          :: Dxynode(2,2,2),Dxznode(2,2,2),Dyznode(2,2,2)
		real*8                          :: GDxxx,GDxyx,GDxzx,GDyyy,GDxyy,GDyzy,GDzzz,GDxzz,GDyzz
		real*8                          :: qx,qy,qz
		real*8                          :: dx,dy,dz,xx,yy,zz,rpt,poro 

       move    = 0.d0

       if (.not. dispersion % action) return !If not dispersion package on do not move particle

       if (particle % zone % num > 0) return !If Immobile phase do not move particle due to gradient of dispersion

       if (particle % phase % num > 0) return  !only moves aqueous phase
       
       if ( dt <= 0.d0 ) return 

       Dxxnode = particle%cell%Dxxnode
	   Dyynode = particle%cell%Dyynode
	   Dzznode = particle%cell%Dzznode
       Dxynode = particle%cell%Dxynode
	   Dxznode = particle%cell%Dxznode
	   Dyznode = particle%cell%Dyznode

       dx = value_array_ (geo%dx,particle%cell%num(1),1,1)
       dy = value_array_ (geo%dy,1,particle%cell%num(2),1)
       dz = value_array_ (geo%dz,1,1,particle%cell%num(3))

	   xx = particle % cell % coord(1)
	   yy = particle % cell % coord(2)
	   zz = particle % cell % coord(3)


      GDxxx = zz        * (1.d0-yy) * ((Dxxnode(2,1,2)-Dxxnode(1,1,2))/dx) +	   &
              (1.d0-zz) * (1.d0-yy) * ((Dxxnode(2,1,1)-Dxxnode(1,1,1))/dx) +	   &
              zz        *       yy  * ((Dxxnode(2,2,2)-Dxxnode(1,2,2))/dx) +	   &
              (1.d0-zz) *       yy  * ((Dxxnode(2,2,1)-Dxxnode(1,2,1))/dx)	   
     
      GDxyx = zz        * (1.d0-yy) * ((Dxynode(2,1,2)-Dxynode(1,1,2))/dx) +	   &
              (1.d0-zz) * (1.d0-yy) * ((Dxynode(2,1,1)-Dxynode(1,1,1))/dx) +	   &
              zz        *       yy  * ((Dxynode(2,2,2)-Dxynode(1,2,2))/dx) +	   &
              (1.d0-zz) *       yy  * ((Dxynode(2,2,1)-Dxynode(1,2,1))/dx)
     
      GDxzx = zz        * (1.d0-yy) * ((Dxznode(2,1,2)-Dxznode(1,1,2))/dx) +	   &
              (1.d0-zz) * (1.d0-yy) * ((Dxznode(2,1,1)-Dxznode(1,1,1))/dx) +	   &
              zz        *       yy  * ((Dxznode(2,2,2)-Dxznode(1,2,2))/dx) +	   &
              (1.d0-zz) *       yy  * ((Dxznode(2,2,1)-Dxznode(1,2,1))/dx)
 
      GDyyy = zz        * (1.d0-xx) * ((Dyynode(1,2,2)-Dyynode(1,1,2))/dy) +	   &
              (1.d0-zz) * (1.d0-xx) * ((Dyynode(1,2,1)-Dyynode(1,1,1))/dy) +	   &
              xx        *       zz  * ((Dyynode(2,2,2)-Dyynode(2,1,2))/dy) + 	   &
              (1.d0-zz) *       xx  * ((Dyynode(2,2,1)-Dyynode(2,1,1))/dy)
  
      GDxyy = zz        * (1.d0-xx) * ((Dxynode(1,2,2)-Dxynode(1,1,2))/dy) +	   &
              (1.d0-zz) * (1.d0-xx) * ((Dxynode(1,2,1)-Dxynode(1,1,1))/dy) +	   &
              xx        *       zz  * ((Dxynode(2,2,2)-Dxynode(2,1,2))/dy) +	   &
              (1.d0-zz) *       xx  * ((Dxynode(2,2,1)-Dxynode(2,1,1))/dy)
     
      GDyzy = zz        * (1.d0-xx) * ((Dyznode(1,2,2)-Dyznode(1,1,2))/dy) +	   &
              (1.d0-zz) * (1.d0-xx) * ((Dyznode(1,2,1)-Dyznode(1,1,1))/dy) +	   &
                    xx  *       zz  * ((Dyznode(2,2,2)-Dyznode(2,1,2))/dy)    +	   &
              (1.d0-zz) *       xx  * ((Dyznode(2,2,1)-Dyznode(2,1,1))/dy)
     
      GDzzz = (1.d0-xx) * (1.d0-yy) * ((Dzznode(1,1,2)-Dzznode(1,1,1))/dz) +	   &
                    xx  * (1.d0-yy) * ((Dzznode(2,1,2)-Dzznode(2,1,1))/dz) +	   &
                    yy  * (1.d0-xx) * ((Dzznode(1,2,2)-Dzznode(1,2,1))/dz) +	   &
                    xx  *       yy  * ((Dzznode(2,2,2)-Dzznode(2,2,1))/dz)
  
      GDxzz = (1.d0-xx) * (1.d0-yy) * ((Dxznode(1,1,2)-Dxznode(1,1,1))/dz) +	   &
                    xx  * (1.d0-yy) * ((Dxznode(2,1,2)-Dxznode(2,1,1))/dz) +	   &
                    yy  * (1.d0-xx) * ((Dxznode(1,2,2)-Dxznode(1,2,1))/dz) +	   &
                    xx  *       yy  * ((Dxznode(2,2,2)-Dxznode(2,2,1))/dz)
     
      GDyzz = (1.d0-xx) * (1.d0-yy) * ((Dyznode(1,1,2)-Dyznode(1,1,1))/dz) +	   &
                    xx  * (1.d0-yy) * ((Dyznode(2,1,2)-Dyznode(2,1,1))/dz) +	   &
                    yy  * (1.d0-xx) * ((Dyznode(1,2,2)-Dyznode(1,2,1))/dz) +	   &
                    xx  *       yy  * ((Dyznode(2,2,2)-Dyznode(2,2,1))/dz)

	  qx = GDxxx + GDxyy + GDxzz
      qy = GDxyx + GDyyy + GDyzz
	  qz = GDxzx + GDyzy + GDzzz

	  if (sorptionACTION .AND. sorptionTYPE == 'LINEAR') then
          if (decayACTION .AND. decayTYPE == 'SERIAL')   then
              rpt = Rap(particle)
          else
              rpt = particle % prop % rp
          end if
      else
          rpt = 1.d0
	  end if

	  poro   = particle % prop % poro

      move(1) = qx / poro / rpt * dt
	  move(2) = qy / poro / rpt * dt
	  move(3) = qz / poro / rpt * dt


      end function


!*************************************************************************************
!      Motion by advection of a reactive particle
!*************************************************************************************
    function move_serial_mom_adv ( advection, dispersion, particle, geo, dt ) result (xstep)
	   use particle_class
	   use advection_class
	   use dispersion_class
	   use geometry_class
	   use loops_particles,  only: nmove,ip
	   use global_variables, only: nspecie
       use heterogeneity_flags
	   implicit none

	   type(advection_cl), intent (in)      :: advection 
	   type(dispersion_cl),intent (in)      :: dispersion
	   type(particle_cl),  intent (in)      :: particle
	   type(geometry_cl),  intent (in)      :: geo 
	   real*8,             intent (in)      :: dt
	   real*8                               :: xstep(3)

	   real*8                               :: qp(3), dx1(3), dx2(3)
       real*8,             save             :: dt_previous !=> -1.d0 
       logical                              :: heterogeneous

       if ( particle%specie%num > nspecie )    return
       heterogeneous = (.not.decay_homogeneous).or.(.not.poro_homogeneous)

       !*******************************************************************
       ! ADVECTIVE MOTION

       if ( (.not.heterogeneous .and. nmove == 1 .and. ip ==1 )     .or.  &
	        ( heterogeneous     .and. particle%control%switchcell ) .or.  &
	        ( dt /= dt_previous ) ) then 

          !* Get qp = q+grad(poro*D): adv. and disp. motion for dt=1
          dx1 = move_one_step_advective (advection,particle,geo,1.d0)
          dx2 = move_one_step_Gradient  (dispersion,particle,geo,1.d0)
          qp = (dx1 + dx2) * particle % prop % rp * particle % prop % poro

          !* Get first norm. spatial moment
          particle % decay % serial_mom % Areac = get_A_reaction(particle,dt,qp)

       end if

       xstep(1) = particle % decay % serial_mom % Areac (particle%specie%num,particle%specie%prenum,1)
       xstep(2) = particle % decay % serial_mom % Areac (particle%specie%num,particle%specie%prenum,2)
       xstep(3) = particle % decay % serial_mom % Areac (particle%specie%num,particle%specie%prenum,3)

       dt_previous = dt !save dt value

    end function 

!*************************************************************************************
!      Motion by dispersion of a reactive particle
!*************************************************************************************
    function move_serial_mom_bro  ( advection, geo, dispersion, particle, dt ) result (move)
       use advection_class
       use geometry_class
       use dispersion_class
	   use particle_class
	   use gslib, only: random_normal
	   implicit none

       type(dispersion_cl), intent(in)    :: dispersion 
       type(particle_cl)  , intent(in)    :: particle
       type(advection_cl) , intent(in)    :: advection
       type(geometry_cl)  , intent(in)    :: geo
	   real*8,              intent(in)    :: dt
	   real*8                             :: move(3)
       real*8               :: poro,qx,qy,qz,al,ath,atv,dm,rpt
	   real*8               :: vx,vy,vz,v,vh,z1x,z1y,z1z,z2x,z2y,z3x,z3y,z3z,z1,z2,z3
	   real*8               :: v2,vh2
	   real*8               :: dx1(3), dx2(3), qp(3)

       move = 0.d0

       if ( .not.dispersion%action ) return
       if ( dt <= 0.d0 ) return

       poro  = particle % prop % poro
       dm    = particle % prop % dm
       qx    = particle % vel  % qpT(1)
	   qy    = particle % vel  % qpT(2)
	   qz    = particle % vel  % qpT(3)

!------- velocities module
       vx = qx / poro
	   vy = qy / poro
	   vz = qz / poro

       v2 = vx * vx + vy * vy + vz * vz
       if (v2 == 0.d0)  then; v = 0.d0
	                    else; v = dsqrt(vx * vx + vy * vy + vz * vz) ; endif

       vh2 = vx * vx + vy * vy
       if (vh2 == 0.d0) then; vh = 1.d-21
	                    else; vh = dsqrt( vx * vx + vy * vy)         ; endif

!------- calculates projections rvs=[z1,z2,z3] to coordinates x,y,z
       z1x = vx/v
       z1y = vy/v
       z1z = vz/v

       z2x =  vy*vh/v/v +  vz*vz*vy/v/v/vh
       z2y = -vx*vh/v/v - vz*vz*vx/v/v/vh
!      z2z = 0.

       z3x = -vz*vx/v/vh
       z3y = -vz*vy/v/vh
       z3z =  vh/v

!------- get the brownian move
       dx1 = move_one_step_advective (advection ,particle,geo,1.d0)
       dx2 = move_one_step_Gradient  (dispersion,particle,geo,1.d0)
       qp  = (dx1 + dx2) * particle % prop % rp * particle % prop % poro

       move = Brown_mov_reactive( particle,dt,v,z1x,z1y,z1z,z2x,z2y,z3x,z3y,z3z,qp )

    end function



    !*************************************************************************************
    ! Subroutine to get the matrix of the center of mass position for a reactive system

    function get_A_reaction (particle,dt,qp) result(Ar)
       use particle_class
       use heterogeneity_flags
       use global_variables, only: TOL, nspecie
       use loops_particles , only: nmove,ip
       implicit none
       type(particle_cl),         intent (in)           :: particle
       real*8,                    intent (in)           :: dt,qp(3)
       logical                                          :: heterogeneous
       integer                                          :: a1,a2,i1,i2,i3,j1,nspe,i,l
       real*8                                           :: Ar   (nspecie, nspecie, 3)
       real*8                                           :: Reff (nspecie, nspecie)
       real*8                                           :: poro, sum
       real*8                                           :: sum1,SR
       real*8, dimension (nspecie,nspecie)              :: D, S, Sinv

       poro        = particle % prop % poro
       D           = 0.d0
       Ar          = 0.d0
       Reff        = 0.d0
       do i=1,nspecie
          D (i,i) = - particle% decay % k(i) / particle % sorption % linear_sorp % R(i)
       end do
       S    = particle % decay % serial_mom % S
       Sinv = particle % decay % serial_mom % Sinv

       !** build the normalized 1st moment matrix for a reactive system : Areac
       do a1=1,nspe; do a2=1,nspecie
          if (a2.le.a1) then
          !** GET EFFECTIVE RETARDATION COEFF.
          !** approx. if P is smaller than fixed tolerance (in commons)
          if (particle % decay % serial_mom % P(a1,a2)<TOL) then
            sum = 0.d0
            do l=a2,a1
                sum = sum + 1/ particle % sorption % linear_sorp % R(l)
            end do
            Reff(a1,a2) = ((a1-a2)+1)/sum

          !** Reff without approx.
          else
            sum1 = 0
            do i1=1,nspe ; do i2=1,nspe ; do i3=1,nspe
                SR = S(a1,i1)*Sinv(i1,i2)*S(i2,i3)*Sinv(i3,a2)/particle % sorption % linear_sorp % R(i2)
                if (SR == 0) cycle
                sum1 = sum1 + SR * Fpr(nspe,i1,i3,D,dt)
            end do; end do; end do
            Reff(a1,a2) = particle % decay % serial_mom % P(a1,a2) / sum1
          end if

          !** get the "advective" motion in the 3 direction
          if (particle % decay % serial_mom % P(a1,a2) /= 0.d0) then
             Ar(a1,a2,1) = qp(1)*dt / (poro * Reff(a1,a2))
             Ar(a1,a2,2) = qp(2)*dt / (poro * Reff(a1,a2))
             Ar(a1,a2,3) = qp(3)*dt / (poro * Reff(a1,a2))
          end if
          end if
       end do; end do

       !** save effective retardation matrices
       particle % decay % serial_mom % Reff  = Reff

    end function

    !*************************************************************************************
    ! Brownian Movement of a reactive particle

    function Brown_mov_reactive ( particle,dt,v,z1x,z1y,z1z,z2x,z2y,z3x,z3y,z3z,qp ) result ( move )
       use particle_class
       use heterogeneity_flags
       use global_variables, only: TOL, nspecie
       use loops_particles , only: nmove,ip
       use gslib, only: random_normal
       implicit none
       type(particle_cl),         intent (in)     :: particle
       real*8,                    intent (in)     :: dt,z1x,z1y,z1z,z2x,z2y,z3x,z3y,z3z,v,qp(3)
       real*8,                    save            :: dt_previous = -1.d0
       real*8                                     :: move(3)
       integer                                    :: i,j,a1,a2,a3,b1,b2,b3,b4,b5
       real*8                                     :: poro,dm,al,ath,atv,qL,qTH,qTV,z1,z2,z3
       logical                                    :: heterogeneous

       real*8                                    :: G, SR, sum2, W, H
       real*8, dimension (nspecie,nspecie,3)     :: B
       real*8, dimension (nspecie,nspecie)       :: D,RR
       real*8, dimension (nspecie,nspecie)       :: S,Sinv,P,Reff
       real*8, dimension (nspecie,nspecie,3)     :: Deff

       heterogeneous = (.not.decay_homogeneous).or.(.not.poro_homogeneous)

       if (particle%specie%num <= nspecie) then                !move only if the particle specie is < num of specie
            z1          = random_normal ()
            z2          = random_normal ()
            z3          = random_normal ()

            if ( (.not.heterogeneous .and. nmove == 1 .and. ip ==1 )   .or.  &
	           ( heterogeneous     .and. particle%control%switchcell ) .or.  &
	           ( dt /= dt_previous ) ) then                                             !calculate move only if conditions changed

            !** get coeff. and matrices
            dm          = particle % prop % dm
            poro        = particle % prop % poro
            al          = particle % prop % aL
			ath         = particle % prop % aTH
	        atv         = particle % prop % aTV
            S           = particle % decay % serial_mom % S
            Sinv        = particle % decay % serial_mom % Sinv
            Reff        = particle % decay % serial_mom % Reff

            D  = 0.d0
            RR = 0.d0
            do i=1,nspecie
               D (i,i) = -particle % decay % k(i) / particle % sorption % linear_sorp % R(i)
               RR(i,i) =  particle % sorption % linear_sorp % R(i)
            end do

            ! get the second spat. moment matrix B
            B    = 0.d0
            Deff = 0.d0
            do j=1,nspecie  ; do i=1,nspecie
                !** GET EFFECTIVE DISPERSION COEFF.
                if (i.le.j) then
                !** approx. if P is smaller than fixed tolerance (in commons)
                if (particle % decay % serial_mom % P(j,i)<TOL) then
                    Deff(j,i,1) = (al*v+dm)
                    Deff(j,i,2) = (ath*v+dm)
                    Deff(j,i,3) = (atv*v+dm)

                !** Deff without approx.
                else
                !** sum to calculate the function G
                    sum2=0 
                    do b1=1,nspecie ; do b2=1,nspecie ; do b3=1,nspecie ; do b4=1,nspecie ; do b5=1,nspecie
                        SR = (S(j,b1)*Sinv(b1,b2)*S(b2,b3)*Sinv(b3,b4)*S(b4,b5)*Sinv(b5,i)) / (RR(b2,b2)*RR(b4,b4))
                        if (SR == 0)  cycle

                        !** function Harp => H
                        if (b3 == b5) then
                            !** function War
                            if (b1 == b5) then
                                W = exp(D(b1,b1)*dt)/2
                            else
                                W = ( exp(D(b5,b5)*dt) * ((D(b5,b5)-D(b1,b1))*dt-1) + exp(D(b1,b1)*dt) ) / ( dt**2 * ((D(b1,b1)-D(b5,b5))**2) )
                            end if
                            H = W
                        else
                            H = (Fpr(nspecie,b1,b3,D,dt)-Fpr(nspecie,b1,b5,D,dt)) / (dt*(D(b3,b3)-D(b5,b5)))
                        end if
                        sum2 = sum2 + SR * H
                    end do; end do; end do; end do; end do
                    G = particle % decay % serial_mom % P(j,i) / sum2

                    Deff(j,i,1) = (al*v+dm)  + (Reff(j,i)/G - 1.d0/(2.d0*Reff(j,i)))*(qp(1)**2.d0/poro**2)*dt
                    Deff(j,i,2) = (ath*v+dm) + (Reff(j,i)/G - 1.d0/(2.d0*Reff(j,i)))*(qp(2)**2.d0/poro**2)*dt
                    Deff(j,i,3) = (atv*v+dm) + (Reff(j,i)/G - 1.d0/(2.d0*Reff(j,i)))*(qp(3)**2.d0/poro**2)*dt

                end if

                B(j,i,1) = 2.d0*Deff(j,i,1)/Reff(j,i)*dt
                B(j,i,2) = 2.d0*Deff(j,i,2)/Reff(j,i)*dt
                B(j,i,3) = 2.d0*Deff(j,i,3)/Reff(j,i)*dt

                end if

            end do; end do

            !save the normalized second moment matrix
            particle % decay % serial_mom % B2 = B

            end if          !end if conditions changed

            move(1) = z1 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,1)) * z1x +			&
                      z2 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,2)) * z2x +			&
                      z3 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,3)) * z3x

            move(2) = z1 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,1)) * z1y +			&
                      z2 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,2)) * z2y +			&
                      z3 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,3)) * z3y

            move(3) = z1 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,1)) * z1z +			&
                      z3 * sqrt(particle % decay % serial_mom % B2(particle%specie%num,particle%specie%prenum,3)) * z3z

       end if               !end if particle specie > nb of specie
       dt_previous = dt

    end function

    !*************************************************************************************
    ! Function Fpr
    function Fpr(nspecie,p,r,D,dt) result(x)
        implicit none
        integer          , intent (in)  :: p,r,nspecie
        real*8          , intent (in)   :: dt
        real*8          , intent (in)   :: D(nspecie,nspecie)
        real*8                          :: x

        if (p==r) then
            x = exp(D(r,r)*dt)
        else
            x = (exp(D(p,p)*dt)-exp(D(r,r)*dt))/(dt*(D(p,p)-D(r,r)))
        end if
    end function

    !*************************************************************************************
    ! Get an approximation of the effectie retardation (harmonic mean)
    function Rap(particle) result(Rapro)
        use particle_class
        use global_variables, only : nspecie
        type(particle_cl),    intent(in) :: particle
        integer                          :: prespe, spe, k
        real*8                           :: sum, Rapro

        j = particle%specie%prenum
        i = particle%specie%num
        if (i==j .OR. i==nspecie) then
            Rapro = particle % sorption % linear_sorp % R(j)
        else
            sum = 0
            do k=j,i
                sum = sum + 1 / particle % sorption % linear_sorp % R(j)
            end do
                Rapro = ((i-j)+1)/sum
        end if
        
    end function


!*************************************************************************************************
!     Time Step Computation
!     *********************
!
!    Courant Number................................: Cu = v * dt / dx
!    Peclet  Number................................: Pe = v * Dx / Dd
!
!    Characteristic time for Diffusion+Dispersion..: tD = Dx**2 / Dd
!    Characteristic time for Advection.............: tA = Dx / v
!
!    Time Step based on Courant = Cu * tA
!    Time Step based on Peclet  = Pe * tD
!
!*************************************************************************************************   
   function calculate_plume_time_step (tc,plume,geo,advection,dispersion,reaction,mass_trans,StressTime) result(dt)
      use global_variables, only: UNEST,DtStep,courant,peclet,calcul_time_method,DaKINETIC,DaDECAY,DaMMT,StartTimeInjection,tsim,particle_in_well_cell,courant_wll,peclet_wll,tshot,AlwaysPrintShots,ntshot
      use code_options,     only: iwcshot
	  use loops_particles,  only: nmove
	  use constants, only: pi
	  use geometry_class
	  use advection_class
	  use dispersion_class
	  use reaction_class
	  use mass_trans_class
	  use plume_class
	  use ctimes_class
      use StressTime_class
      use velocity_times
	  implicit none

	  type(geometry_cl),           intent(in)    :: geo
	  type(advection_cl),          intent(in)    :: advection
	  type(dispersion_cl),         intent(in)    :: dispersion
	  type(reaction_cl),           intent(in)    :: reaction
	  type(mass_trans_cl),         intent(in)    :: mass_trans
	  type(plume_cl),              intent(in)    :: plume
      type(StressTime_cl),         intent(inout) :: StressTime 
	  type(ctimes_cl),             intent(inout) :: tc

	  real*8                       :: dt
	  real*8, save                 :: dtold = UNEST
	  real*8                       :: dtstress
      integer                      :: loc
      integer, save                :: itshot = 1


          !...estimate dt based on method selected
          
          dtold = dabs(DtStep)

          if (nmove ==1)  call initialize_characteristic_times_ (tc,plume,geo,advection,dispersion,reaction)           
          
          select case (calcul_time_method)

	         case ('CONSTANT_TIME','CONSTANT_DT') 

	                 dt = dabs(DtStep)

             case ('CONSTANT_MOVE','CONSTANT_CU')

                     !call calculate_adv_characteristic_time_ ( tc, plume, geo, advection, dispersion, reaction ) !calculate characteristic times of plume particles if flow changes
			         dt = plume_constant_move_time_step ( courant, tc%adv )
                     !if (particle_in_well_cell==.FALSE.) dt = plume_constant_move_time_step ( courant, tc%adv )
                     !if (particle_in_well_cell==.TRUE.)  dt = plume_constant_move_time_step ( courant_wll, tc%adv )

			 case ('CONSTANT_DAMT','CONSTANT_DADECAY')

                     !call calculate_rx_characteristic_times_ ( tc, plume, geo, reaction, mass_trans ) !calculate reactive (decay network or/and mass transfer) characteristic times of plume particles
                     dt = calculate_dt_based_on_rx_characteristic_times_ ( tc, DaDECAY, DaMMT )

             case ('OPTIMUM_DT')

                     !call calculate_characteristic_times_ ( tc, plume, geo, advection, dispersion, reaction ) !calculate characteristic times of plume particles if flow changes
                     dt = calculate_dt_based_on_characteristic_times_ ( tc, courant, peclet, DaKINETIC )
                     !if (particle_in_well_cell==.FALSE.) dt = calculate_dt_based_on_characteristic_times_ ( tc, courant, peclet, DaKINETIC )
                     !if (particle_in_well_cell==.TRUE.)  dt = calculate_dt_based_on_characteristic_times_ ( tc, courant_wll, peclet_wll, DaKINETIC )

			 case default

		             stop ' Selected Time Step Method NOT valid'

          end select

          !...if time is negative, take the old time step
          
          if (dt <=0.d0) dt = dtold

          !...the first time restericted also to DtStep

          if (nmove ==1) dt = min(dabs(DtStep),dt)

          !...check that time step is not larger than the transport stress time interval

          !call find_loc_StressTime_ (StressTime, plume%time)
          call find_loc_StressTime_ (StressTime, velotime%time_flow)
          dtstress = StressTime%time(StressTime%loc+1)-StressTime%time(StressTime%loc)          
          dt = min(dt,dtstress)
          
          !...modify time step if greater than transport stress time

          if (velotime%time_flow < StressTime%time(StressTime%nt) .and. velotime%time_flow + dt > StressTime%time(StressTime%loc+1)) then
               dt = StressTime%time(StressTime%loc+1)-velotime%time_flow 
          end if
                  
        !...modify time step if time greater than next plume snapshot
          
        if (iwcshot>0) then
            if (.not.AlwaysPrintShots) then
                if (velotime%time_flow + dt > tshot(itshot)) then
                    if (itshot<ntshot) dt = tshot(itshot)-velotime%time_flow
                    if (itshot<ntshot) itshot = itshot + 1
                end if
            end if
        end if
        
          !...end of simulation

          if (plume%time + dt > tsim) dt = tsim - plume%time 

          !...save dt
          
          dtold = dt 
       
   end function

!**************************************************************************************
!  calculate dt based on Dave's approach
!**************************************************************************************
   function calculate_dt_based_on_characteristic_times_ (tc,Cu,Pe,Da) result (dt)
     use global_variables, only: UNEST
     use ctimes_class
     implicit none
       type(ctimes_cl), intent(in) :: tc
       real*8,          intent(in) :: Cu,Pe,Da 
       real*8                      :: tcdisp,tcadv,tcreact
       real*8                      :: dt
       real*8                      :: a(3)
       logical                     :: Active(3)
       integer                     :: n
       
      a = UNEST      
      Active = .FALSE.
      
      if (associated(tc%adv)    .and. tc%adv /= UNEST ) then
                Active(1) = .TRUE.
                     a(1) = Cu*tc%adv
      end if
      if (associated(tc%disp)   .and. tc%disp /= UNEST ) then
                 Active(2) = .TRUE.
                      a(2) = Pe*tc%disp 
      end if
      if (associated(tc%kreact) .and. tc%kreact /= UNEST ) then
                 Active(3) = .TRUE.
                      a(3) = Da*tc%kreact
      end if

      n=count(Active)

      if (count(Active)>0) then
         dt = minval(a,mask=Active)
      else
         dt = UNEST
         !stop ' characteristic time unestimated'
      end if

   end function

!**************************************************************************************
!  calculate dt based on Dave's approach
!**************************************************************************************
   function calculate_dt_based_on_rx_characteristic_times_ (tc, DaDECAY, DaMMT) result (dt)
      use global_variables, only: calcul_time_method,UNEST
      use ctimes_class
      implicit none

      type(ctimes_cl), intent(in) :: tc
      real*8,          intent(in) :: DaDECAY,DaMMT
      real*8                      :: dt
      
      select case (calcul_time_method)
      case ('CONSTANT_DAMT')
        dt = tc%mass_trans*DaMMT
      case ('CONSTANT_DADECAY')
        dt = tc%decay*DaDECAY
      case default
        stop ' characteristic time unestimated'
      end select

   end function
   
!**************************************************************************************
!  calculate constant move time step based on the courant number
!**************************************************************************************
   function plume_constant_move_time_step (courant,tcadv) result (dt)
	  implicit none
	     real*8, intent(in) :: courant
	     real*8, intent(in) :: tcadv
	     real*8             :: dt

	     dt    = courant * tcadv
   
   end function

!**************************************************************************************
!**************************************************************************************

 end module