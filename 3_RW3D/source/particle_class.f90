
!***************************************************************************************
!  PARTICLE CLASS
!***************************************************************************************

  module particle_class
	implicit none

    private
	public :: particle_cl        !class
	public ::                                            &   !methods
              allocate_particle_                       , &
              add_move_to_particle_                    , &
			  print_position_particle_                 , &
			  print_position_particle_tecplot_         , &
			  update_cell_location_particle_           , &
			  update_properties_particle_              , &
			  update_velocity_particle_                , &
			  update_dispersion_nodes_particle_        , &
			  delete_particle_                         , &
			  assign_to_particle_                      , &
			  velocity_centroid_                       , &
			  velocity_particle_                       , &
			  check_mobility_particle_                 , &
			  print_particle_passtime_                 , &
			  set_particle_to_zero_                    , &
			  print_block_residence_time_particle_     , &
			  update_block_region_location_particle_   , &
			  update_block_Ti_Xpi_particle_            , &
			  calc_module_Darcy_velocity_particle_     , &
			  from_plumeparticle_to_particle_          , &
			  from_particle_to_plumeparticle_          , &
    
			  update_mass_trans_                       , &
              update_decay_                            , &
              update_sorption_


    type position_cl
	    real*8  :: tp        ! actual time
		real*8  :: tpold     ! previous time
		real*8  :: xpref(3)  ! xpref = position at t=0
	    real*8  :: xp(3)     ! position in actual time
		real*8  :: xpold(3)  ! position in previous time
	end type

	type cell_cl
	    real*8  :: coord(3)    !local coordinates of the particle in the cell, origin is the lower corner of the cell
		integer :: num(3)      !column, row and layer of the cell where particle is
		integer :: numold(3)   !col,row,lay previous                                       
		real*8  :: qxnode(2,2,2),qynode(2,2,2),qznode(2,2,2)    !darcy velocities at the corner of the cell
		real*8  :: qxface(2),qyface(2),qzface(2)                !darcy velocities at the cell faces
	    real*8  :: Dxxnode(2,2,2),Dyynode(2,2,2),Dzznode(2,2,2) !dispersion tensor at the corners of the cells 
		real*8  :: Dxynode(2,2,2),Dxznode(2,2,2),Dyznode(2,2,2)
	end type

	type control_cl
	    integer          :: passplane    !index of plane passed last
		integer          :: passwell     !index of well passed last
		logical          :: remove       !specifies if particle is removed
		logical          :: switchcell   !specifies if particle switches cell
		logical          :: stuck        !specifies if particle gets stuck somewhere
		integer          :: npbounce     !specifies the number of bounces of a particle
        integer          :: npdomextens  !specifies the number of particle entendening the domain (for semi-infinite domain)
        logical          :: part_in_wellcell !specifies if particle located in a well cell
        integer          :: wellID_in_partcell !ID of the well located in the same cell than the particle (if so)
        integer          :: inblock      !index of plane passed last
	end type

    type zone_cl
		integer, pointer                 :: num    => null() !Mobile/Immobile zone number (number zero corresponds to mobile zone)
		integer, pointer                 :: prenum => null() !Mobile/Immobile zone number (number zero corresponds to mobile zone)
        real*8,  pointer, dimension(:,:) :: alpha  => null() !alpha parameter of mass transfer for the particle 
		real*8,  pointer, dimension(:,:) :: beta   => null() !beta parameter of mass transfer for the particle
		real*8,  pointer, dimension(:)   :: DaI    => null() !array with Damkohler numbers
	end type    

    type specie_cl
		integer, pointer                 :: num  => null() !Number of the specie
		integer, pointer                 :: prenum  => null()   !previous number of the species
	end type

    type phase_cl
        integer   :: num = 0
    end type

	type properties_cl      !properties at the particle positon
	    real*8              :: aL,aTH,aTV,dm,dmTH,dmTV,poro,mp,rp=1 
	end type

	type velocity_cl
	    real*8  :: qpL(3)   !particle darcy velocity obtained from linear interpolation 
	    real*8  :: qpT(3)   !particle darcy velocity obtained from trilinear interpolation 
	end type

    type block_cl       !Residence time of a particle within a block region (not a cell) 
        logical :: switchblock
        integer :: numold(3)    !block number in previous time
		integer :: num(3)       !block number (column, row, layer)
		real*8  :: Ti           !Arrival time Time of the particle at the block region
	    real*8  :: Xpi(3)       !Arrival poisition of the particle at the block region
	end type

    !** Mass transfer
    type multirate_cl
        real*8, pointer     :: alphap(:) => null()
        real*8, pointer     :: poroim(:) => null()
    end type

    type diffusion_cl
        real*8, pointer     :: alphap    => null()
        real*8, pointer     :: poroim    => null()
    end type

    type power_cl
        real*8, pointer     :: btot      => null() 
        real*8, pointer     :: Amin      => null()
        real*8, pointer     :: Amax      => null()
        real*8, pointer     :: power     => null()
    end type

    type lognormal_cl
        real*8, pointer     :: btot      => null()
        real*8, pointer     :: meanlog   => null()
        real*8, pointer     :: stdvlog   => null()
    end type

    type composite_cl
        type(multirate_cl)  :: multirateComp
        type(diffusion_cl)  :: sphericalComp
        type(diffusion_cl)  :: layeredComp
        type(diffusion_cl)  :: cylindricalComp
    end type

    type mass_transfer_cl
        type(multirate_cl)  :: multirateMT
        type(diffusion_cl)  :: diffusionMT
        type(power_cl)      :: powerMT
        type(lognormal_cl)  :: lognormalMT
        type(composite_cl)  :: compositeMT
    end type

    !** Decay
    type serial_mom_cl
        real*8, pointer     :: P(:,:)       => null()
        real*8, pointer     :: S(:,:)       => null()
        real*8, pointer     :: Sinv(:,:)    => null()
        real*8, pointer     :: Areac(:,:,:) => null()
        real*8, pointer     :: B2(:,:,:)    => null()
        real*8, pointer     :: Reff(:,:)    => null()
    end type

    type decay_network_cl
        real*8, pointer     :: k(:)     => null()
        real*8, pointer     :: kim(:,:) => null()
        real*8, pointer     :: y(:,:)   => null()
        type(serial_mom_cl) :: serial_mom
    end type

    !** Sorption
    type linear_sorption_cl
        real*8, pointer     :: R(:)     => null()
        real*8, pointer     :: Rim(:,:) => null()
    end type

    type CHTM_cl
        real*8, pointer     :: bd      => null()
        real*8, pointer     :: kf      => null()
        real*8, pointer     :: kb      => null()
    end type

    type sorption_cl
        type(linear_sorption_cl)    :: linear_sorp
        type(CHTM_cl)               :: CHTM
    end type

    ! Particle object
	type particle_cl
	    integer                     :: id
	    type(position_cl)           :: position
		type(cell_cl)               :: cell
		type(block_cl)              :: block
		type(control_cl)            :: control
		type(properties_cl)         :: prop
		type(velocity_cl)           :: vel
		type(zone_cl)               :: zone
		type(specie_cl)             :: specie
		type(phase_cl)              :: phase
		type(mass_transfer_cl)      :: mass_transfer
		type(sorption_cl)           :: sorption
		type(decay_network_cl)      :: decay
	end type

    interface update_properties_particle_
        module procedure update_properties_particle_1,update_properties_particle_2; end interface !,update_properties_particle_3; end interface

    interface print_particle_passtime_
	    module procedure print_particle_passtime_,print_particle_passtime_xyz_; end interface

    interface allocate_particle_
	    module procedure alloc_particle,allocate_particle_min; end interface

	contains


    subroutine delete_particle_ (this)
	   implicit none
	   type(particle_cl), intent(inout) :: this

		 if (associated(this%specie%num))        deallocate(this%specie%num)

		 if (associated(this%zone%num))          deallocate(this%zone%num)
		 if (associated(this%zone%prenum))       deallocate(this%zone%prenum)
		 if (associated(this%zone%alpha))        deallocate(this%zone%alpha)
		 if (associated(this%zone%beta))         deallocate(this%zone%beta)
		 
		 !delete mass transfer parameteres
         !* multirate
		 if (associated(this%mass_transfer%multirateMT%alphap))     deallocate(this%mass_transfer%multirateMT%alphap)
         if (associated(this%mass_transfer%multirateMT%poroim))     deallocate(this%mass_transfer%multirateMT%poroim)
         !* diffusion model
         if (associated(this%mass_transfer%diffusionMT%alphap))     deallocate(this%mass_transfer%diffusionMT%alphap)
         if (associated(this%mass_transfer%diffusionMT%poroim))     deallocate(this%mass_transfer%diffusionMT%poroim)
         !* power model
         if (associated(this%mass_transfer%powerMT%btot))           deallocate(this%mass_transfer%powerMT%btot)
         if (associated(this%mass_transfer%powerMT%Amin))           deallocate(this%mass_transfer%powerMT%Amin)
         if (associated(this%mass_transfer%powerMT%Amax))           deallocate(this%mass_transfer%powerMT%Amax)
         if (associated(this%mass_transfer%powerMT%power))          deallocate(this%mass_transfer%powerMT%power)
         !* LogNormal model
         if (associated(this%mass_transfer%lognormalMT%btot))       deallocate(this%mass_transfer%lognormalMT%btot)
         if (associated(this%mass_transfer%lognormalMT%meanlog))    deallocate(this%mass_transfer%lognormalMT%meanlog)
         if (associated(this%mass_transfer%lognormalMT%stdvlog))    deallocate(this%mass_transfer%lognormalMT%stdvlog)
         !* Composite model

    end subroutine

   !************************************
    subroutine add_move_to_particle_ (this,dxp,dt)
       implicit none
	   type(particle_cl), intent(inout) :: this
	   real*8,            intent(in)    :: dxp(3)
	   real*8,            intent(in)    :: dt
	   real*8                           :: xp1(3),tp1
	       xp1 = this % position % xp
		   tp1 = this % position % tp
	       this % position % xp    = xp1 + dxp
	       this % position % tp    = tp1 + dt
		   this % position % xpold = xp1
		   this % position % tpold = tp1
    end subroutine

   !************************************
    subroutine assign_to_particle_ (this,xp,mp,zone,specie) 
        implicit none
		type(particle_cl), intent(inout) :: this
		real*8,            intent(in)    :: xp(3)
		real*8,            intent(in)    :: mp
		integer,           intent(in)    :: zone,specie
		     this % position % xp  = xp
			 this % prop     % mp  = mp
			 this % zone     % num = zone
			 this % specie   % num = specie
    end subroutine

   !************************************
	subroutine from_plumeparticle_to_particle_   (this,time,part,izone,ispecie)
	   use list_class
	   implicit none
	   type(particle_cl),        intent(inout) :: this
       real*8,                   intent(in)    :: time
       type(partID_cl), pointer                :: part
       integer,                  intent(in)    :: izone,ispecie
	     this % id                  = part % id
	     this % position % xpref(1) = part % xp
		 this % position % xpref(2) = part % yp
		 this % position % xpref(3) = part % zp
	     this % position % xp(1)    = part % xp
		 this % position % xp(2)    = part % yp
		 this % position % xp(3)    = part % zp
		 this % prop     % mp       = part % mp
		 this % position % tp       = time
		 this % position % xpold    = this % position % xp  
	     if (.not.associated(this % zone % num))   allocate( this % zone % num)
	     if (.not.associated(this % specie % num)) allocate( this % specie % num)
	     if (.not.associated(this % zone % prenum))   allocate( this % zone % prenum)
	     if (.not.associated(this % specie % prenum)) allocate( this % specie % prenum)
	     this % zone     % num      = izone
	     this % specie   % num      = ispecie
	     this % zone     % prenum   = izone
	     this % specie   % prenum   = ispecie
		 this % prop     % rp       = part % rp
		 this % block    % Xpi      = this % position % xp
		 this % control  % remove   = .FALSE.
         this % control  % stuck    = .FALSE.
		 this % control  % passplane = 0 !initialize the flag
         this % control  % passwell  = 0 !initialize the flag       
         this % control  % part_in_wellcell = .FALSE.
         this % control  % wellID_in_partcell = 0
         this % control  % inblock = 0 !initialize the flag
	end subroutine

   !************************************
	subroutine from_particle_to_plumeparticle_ (particle,plumepart)
      use list_class
      implicit none
	   type(particle_cl), intent(in) :: particle
       type(partID_cl), pointer      :: plumepart
         plumepart % id = particle % id
         plumepart % xp = particle % position % xp(1)
         plumepart % yp = particle % position % xp(2)
         plumepart % zp = particle % position % xp(3)
         plumepart % mp = particle % prop % mp
         plumepart % rp = particle % prop % rp
    end subroutine

   !************************************
    subroutine print_position_particle_ (this,fname, nmove)
       use gslib, only: open_fname
	   use code_options, only: iwpath
	   use loops_particles, only: ip
	   use global_variables, only: pathfreq,pathpart
	   implicit none
	   type(particle_cl),          intent(in) :: this
	   character(len=*),           intent(in) :: fname
	   integer,          optional, intent(in) :: nmove
       integer                                :: i,iunit
	   integer, save                          :: step = 0
	   integer, save                          :: iter = 0
	   real*8                                 :: qqp
           if ( iwpath == 0 ) return
           if ( pathpart > 0 .and. pathpart /= ip ) return
		   if (step == 0) then
		      call open_fname(fname,iunit)
			  iter = iter +1
			  qqp = (this%vel%qpL(1)**2)+(this%vel%qpL(2)**2)+(this%vel%qpL(3)**2)
              write(iunit,1) this%id,(this%position%xp(i),i=1,3),this%prop%rp,this%zone%num,this%specie%num,this%position%tp,qqp 
		   end if   
		   step = step + 1
		   if (step == pathfreq) step = 0 
		   
         1 format(i10,x4(g15.6,x),2(i5,x),g15.6,x,g20.11)  	
    end subroutine

   !************************************
    subroutine print_position_particle_tecplot_ (fname,it,time,xishot,yishot,zishot,mass,rishot,zone,specie)
          use gslib, only: open_fname,generate_unit
		  use loops_particles, only: ip
		  use array_class
		  character(len=*), intent(in)  :: fname
		  real*8,           intent(in)  :: time,xishot,yishot,zishot,mass,rishot
		  integer,          intent(in)  :: it,zone,specie
		  integer                       :: iunit
	      logical                       :: connected
	       inquire(file=fname,opened=connected)
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown')
			   write(iunit,*) 
			   write(iunit,1) 'SNAPSHOT'      , &
			                  '    TIME    ', &
			                  '     XP     ', &
							  '     YP     ', &
							  '     ZP     ', &
							  '    MASS    ', &
							  '   RETARD   ', &
							  ' ZONE'       , &
							  'SPECIE'      , &
							  '  PART'
             1 format(a9,6(a20,x),a6,x,a8,x,a6)
	           write(iunit,*)
		   end if  
		        
           write(iunit,3) it,time,xishot,yishot,zishot,mass,rishot,zone,specie,ip 
	    3  format(i8,x,e20.11,x,5(e20.11,x),2(i6,x),i10)	 
	  end subroutine




!***********************************************************************************************
!     locate cell where particle moves
!***********************************************************************************************
      subroutine update_cell_location_particle_ (this,geo) !only used for update_cell_location
	  use array_class
	  use geometry_class
	  use gslib, only: locate
	  implicit none
	  type(particle_cl), intent(inout) :: this
	  type(geometry_cl), intent(in)    :: geo
	  integer                          :: nx,ny,nz,i,j,k,iprev,jprev,kprev
	  real*8                           :: dx,dy,dz,x1,y1,z1,x2,y2,z2,xp,yp,zp

             iprev = this%cell%num(1)
	         jprev = this%cell%num(2)
	         kprev = this%cell%num(3)

!............find cell location: new cell location is (i,j,k)

             nx = length_array_ (geo%dx)
			 ny = length_array_ (geo%dy)
			 nz = length_array_ (geo%dz)

             xp = this%position%xp(1)
			 yp = this%position%xp(2)
			 zp = this%position%xp(3)

			 if (nx == 1) then
                dx = value_array_ (geo%dx)
				i  = dint(xp/dx) + 1
                x1 = dx*(i-1)
			 else
				call locate(geo%xmesh%values,nx+1,1,nx+1,xp,i)
			    x1 = geo%xmesh%values(i,1,1)
             end if

			 if (ny == 1) then
                dy = value_array_ (geo%dy)
				j  = dint(yp/dy) + 1
                y1 = dy*(j-1)
			 else
				call locate(geo%ymesh%values,ny+1,1,ny+1,yp,j)
				y1 = geo%ymesh%values(1,j,1)
             end if

			 if (nz == 1) then
                dz = value_array_ (geo%dz)
				k  = dint(zp/dz) + 1
                z1 = dz*(k-1)
			 else
				call locate(geo%zmesh%values,nz+1,1,nz+1,zp,k)
			    z1 = geo%zmesh%values(1,1,k)
             end if     

!.............check if particle switches cell

	         if (i == iprev .and. j == jprev .and. k == kprev  ) then
	             this%control%switchcell = .FALSE.
	         else
	             this%control%switchcell = .TRUE.
	         end if

             this%cell%numold(1) = iprev
			 this%cell%numold(2) = jprev
			 this%cell%numold(3) = kprev

             this%cell%num(1) = i
	         this%cell%num(2) = j
	         this%cell%num(3) = k

!............calculate local coordinates of the particle within the cell
!if (i>nx) then
!   dx = 1
!end if
             dx = dx_ (geo,i,1,1)
             dy = dy_ (geo,1,j,1)
             dz = dz_ (geo,1,1,k)

             this%cell%coord(1) = (xp-x1)/dx
             this%cell%coord(2) = (yp-y1)/dy
             this%cell%coord(3) = (zp-z1)/dz

	  end subroutine


!***********************************************************************************************
!     locate block region where particle moves
!***********************************************************************************************
      subroutine update_block_region_location_particle_ (this,geo) !only used for upscaling
	  use array_class
	  use geometry_class
	  use gslib, only: locate
	  use code_options, only: irestime
	  implicit none
	  type(particle_cl), intent(inout) :: this
	  type(geometry_cl), intent(in)    :: geo
	  integer                          :: nx,ny,nz,i,j,k,iprev,jprev,kprev
	  real*8                           :: dx,dy,dz,x1,y1,z1,x2,y2,z2,xp,yp,zp

      if (.not.irestime ) return

             iprev = this%block%num(1)
	         jprev = this%block%num(2)
	         kprev = this%block%num(3)

!............find cell location: new cell location is (i,j,k)

             nx = length_array_ (geo%dx)
			 ny = length_array_ (geo%dy)
			 nz = length_array_ (geo%dz)

             xp = this%position%xp(1)
			 yp = this%position%xp(2)
			 zp = this%position%xp(3)

			 if (nx == 1) then
                dx = value_array_ (geo%dx)
				i  = dint(xp/dx) + 1
                x1 = dx*(i-1)
			 else
				call locate(geo%xmesh%values,nx+1,1,nx+1,xp,i)
			    x1 = geo%xmesh%values(i,1,1)
             end if

			 if (ny == 1) then
                dy = value_array_ (geo%dy)
				j  = dint(yp/dy) + 1
                y1 = dy*(j-1)
			 else
				call locate(geo%ymesh%values,ny+1,1,ny+1,yp,j)
				y1 = geo%ymesh%values(1,j,1)
             end if

			 if (nz == 1) then
                dz = value_array_ (geo%dz)
				k  = dint(zp/dz) + 1
                z1 = dz*(k-1)
			 else
				call locate(geo%zmesh%values,nz+1,1,nz+1,zp,k)
			    z1 = geo%zmesh%values(1,1,k)
             end if     

!.............check if particle switches block

	      if (i == iprev .and. j == jprev .and. k == kprev  ) then
	         this%block%switchblock = .FALSE.
	      else
	         this%block%switchblock = .TRUE.
             this%block%numold(1) = iprev
			 this%block%numold(2) = jprev
			 this%block%numold(3) = kprev

             this%block%num(1) = i
	         this%block%num(2) = j
	         this%block%num(3) = k	
	      end if



	  end subroutine

!***********************************************************************************************
!     initialize entrance time in a block region (upscaling)
!***********************************************************************************************
    subroutine update_block_Ti_Xpi_particle_ (this) !only used for upscaling
	  use array_class
	  use geometry_class
	  use code_options, only: irestime
	  type(particle_cl), intent(inout) :: this
         if (.not.irestime) return
	     if (.not.this%block%switchblock) return
	     if (this%block%switchblock) then	  
	          this%block%Ti  = this%position%tp
		      this%block%Xpi = this%position%xp
	     end if
    end subroutine

!***********************************************************************************************
!     print residence time when the particle leaves the block
!***********************************************************************************************
    subroutine print_block_residence_time_particle_ (this,fname,geoblock) !only used for upscaling
	  use array_class
	  use geometry_class
	  use loops_particles, only: ip
	  use gslib, only: generate_unit
	  use code_options, only: irestime
	  type(particle_cl), intent(inout) :: this
	  type(geometry_cl), intent(in)    :: geoblock
	  character(len=*),  intent(in)    :: fname
	  integer                          :: iunit
	  logical                          :: connected

      if (.not.irestime) return
	  if (.not.this%block%switchblock) return
	  if (this%block%switchblock) then
      	   inquire(file=fname,opened=connected)
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(200)
	           open(iunit,file=fname,status='unknown')
			   write(iunit,1) 'BLOCK PROPERTIES: ',' nx = ',geoblock%nx, &
			                                       ' ny = ',geoblock%ny, &
												   ' nz = ',geoblock%nz 
			 1 format(A18,3(a6,i6))
			   write(iunit,*) ' 11'
			   write(iunit,*) 'PARTICLE NUM  '
			   write(iunit,*) 'COLUMN  '
			   write(iunit,*) 'ROW '
			   write(iunit,*) 'LAYER '
			   write(iunit,*) 'RESIDENCE TIME'
			   write(iunit,*) 'BLOCK_DXp '
			   write(iunit,*) 'BLOCK_DYp '
               write(iunit,*) 'BLOCK_DZp '
			   write(iunit,*) 'Xpi'
			   write(iunit,*) 'Ypi'
			   write(iunit,*) 'Zpi'
		   end if  
       
	  write(iunit,2) ip, this%block%numold(1), this%block%numold(2),this%block%numold(3), &
	           this%position%tp    -  this%block%Ti      , &
               this%position%xp(1) -  this%block%Xpi(1)  , &
			   this%position%xp(2) -  this%block%Xpi(2)  , &
			   this%position%xp(3) -  this%block%Xpi(3)  , &
               this%block%Xpi(1),   &
			   this%block%Xpi(2),   &
			   this%block%Xpi(3)
      2 format(4(i7,x),7(x,g20.11))      
	  end if
    end subroutine

!***************************************************************************************************
!   two subroutines to update properties of a particle: dispersivities, diffusion coeff, retardation
!***************************************************************************************************
    subroutine update_properties_particle_1 (this,geo,advection,dispersion,reaction)
       use global_variables
       use array_class
	   use advection_class
	   use dispersion_class
	   use reaction_class
	   use geometry_class
	   use heterogeneity_flags
	   !use loops_particles, only: nmove
	   use to_solve
	   implicit none
	   
	   type(particle_cl),   intent(inout) :: this
	   type(advection_cl),  intent(in)    :: advection
	   type(dispersion_cl), intent(in)    :: dispersion
	   type(reaction_cl),   intent(in)    :: reaction
	   type(geometry_cl),   intent(in)    :: geo
	   real*8                             :: kd,bd,poro,fm,Rm
	   integer                            :: i

       !update porosity
       !if (poro_homogeneous == .FALSE. .or. nmove == 1) then
       this % prop % poro = update_property_block_discrete ( this,geo,advection%poro )!; end if

       !update dispersion particle parameters
       !if (disp_homogeneous == .FALSE. .or. nmove == 1) then
       if (dispersion%action) then
           if(.not.SpeciesDispersionDependent) then
	        this % prop % aL     = update_property ( this,geo,dispersion % aL    ) 
	        this % prop % aTH    = update_property ( this,geo,dispersion % aTH   )
	        this % prop % aTV    = update_property ( this,geo,dispersion % aTV   )
	        this % prop % dm     = update_property ( this,geo,dispersion % dm    )
	        this % prop % dmTH   = update_property ( this,geo,dispersion % dmTH  )
	        this % prop % dmTV   = update_property ( this,geo,dispersion % dmTV  )
	        else
	           if (this%specie%num <= nspe_aq) then
	             this % prop % aL     = update_property ( this,geo,dispersion % aL    ) * dispersion % MultA ( this%specie%num ) 
	             this % prop % aTH    = update_property ( this,geo,dispersion % aTH   ) * dispersion % MultA ( this%specie%num ) 
	             this % prop % aTV    = update_property ( this,geo,dispersion % aTV   ) * dispersion % MultA ( this%specie%num ) 
	             this % prop % dm     = update_property ( this,geo,dispersion % dm    ) * dispersion % MultD ( this%specie%num ) 
	             this % prop % dmTH   = update_property ( this,geo,dispersion % dmTH  ) * dispersion % MultD ( this%specie%num ) 
	             this % prop % dmTV   = update_property ( this,geo,dispersion % dmTV  ) * dispersion % MultD ( this%specie%num ) 
	            end if
	        end if
       end if !; end if

       !update sorption particle parameters
       if (sorptionACTION.AND. sorptionTYPE == 'LINEAR') then
            this % prop % rp = update_property_block_discrete (this,geo,reaction%sorption%LinearSorption%R(this%specie%num))
       else
            this % prop % rp = 1.0d0
       end if

	end subroutine


    subroutine update_properties_particle_2 (this,geo,advection,reaction)
       use advection_class
	   use geometry_class
	   use array_class
	   use reaction_class
	   use heterogeneity_flags
	  !use loops_particles, only: nmove
	   use to_solve
	   implicit none

	   type(particle_cl),   intent(inout) :: this
	   type(reaction_cl),   intent(in)    :: reaction
	   type(geometry_cl),   intent(in)    :: geo
	   type(advection_cl),  intent(in)    :: advection
	   real*8                             :: kd,poro,bd,fm,Rm
	   integer                            :: i

       !update porosity
       !if (poro_homogeneous == .FALSE. .or. nmove == 1) then
       this % prop % poro = update_property_block_discrete ( this,geo,advection%poro ) !; end if

       !update sorption particle parameters
       if (sorptionACTION.AND. sorptionTYPE == 'LINEAR') then
          this % prop % rp = update_property_block_discrete (this,geo,reaction%sorption%LinearSorption%R(this%specie%num)) 
       else
          this % prop % rp = 1.0d0
       end if

    end subroutine

!*******************************************************************************************
!   function to calculate the value of a property variable such as aL, R, poro at
!   the particle position using TRILINEAR INTERPOLATION with adjoint cells
!*******************************************************************************************
    function update_property (this,geo,array,coord) result (value)
       use array_class
	   use geometry_class
	   use gslib, only: interpolation_trilinear
	   implicit none
	   type(particle_cl), intent(inout) :: this
	   type(geometry_cl), intent(in)    :: geo   !geometry
	   type(array_cl),    intent(in)    :: array !array of the property that is evaluated
	   real*8, optional,  intent(in)    :: coord(3)
	   real*8                           :: value
	   integer                          :: n
             n = length_array_ (array)
	         if (n == 1) then
                       value = value_array_ (array,1,1,1)
                       return;  end if
             if ( present(coord) ) then
	                   value = interpolate_array_ (geo,array,this%cell%num,coord)
	              else
                       value = interpolate_array_ (geo,array,this%cell%num,this%cell%coord);  end if
	end function

!*******************************************************************************************
!   function to calculate the value of a property variable such as R and poro at
!   the particle position using grid-block value
!*******************************************************************************************
    function update_property_block_discrete (this,geo,array) result (value)
       use array_class
	   use geometry_class
	   implicit none
	   type(particle_cl), intent(inout) :: this
	   type(geometry_cl), intent(in)    :: geo   !geometry
	   type(array_cl),    intent(in)    :: array !array of the property that is evaluated
	   real*8                           :: value
	   integer                          :: n
             n = length_array_ (array)
	         if (n == 1) then
                       value = value_array_ (array,1,1,1)
			 else
			           value = value_array_ (array,this%cell%num)
			 end if
	end function

!  ***********************************************************************************************************************
    subroutine update_velocity_particle_ (this,geo,advection,dispersion,well)
	   use geometry_class
	   use advection_class
	   use dispersion_class
	   use heterogeneity_flags
       use well_vect_class
	   !use loops_particles, only: nmove
	   implicit none
	   type(particle_cl),             intent(inout) :: this
	   type(geometry_cl),             intent(in)    :: geo
	   type(advection_cl),            intent(in)    :: advection
	   type(dispersion_cl), optional, intent(in)    :: dispersion 
       type(well_vect_cl),  optional, intent(in)    :: well
       
         !if (vel_homogeneous== .TRUE. .and. disp_homogeneous== .TRUE. .and. nmove /= 1) return
         !if (vel_homogeneous== .FALSE. .or. nmove == 1) then
		    if (advection%action) then
                if (this%control%switchcell) call vel_faces_ (advection,this%cell%num,this%cell%qxface,this%cell%qyface,this%cell%qzface) 
                if (this%control%part_in_wellcell) then 
                    this%vel%qpL = velo_in_wellcell_ (well,geo,this)
                else
			        this%vel%qpL = vel_point_linear (this)
                end if   
            end if
         !if (disp_homogeneous== .FALSE. .or. nmove == 1) then
		    if (present(dispersion)) then
		        if (dispersion%action) then
		            if (this%control%switchcell) call vel_nodes_ (advection,geo,this%cell%num,this%cell%qxnode,this%cell%qynode,this%cell%qznode)
                    this%vel%qpT = vel_point_trilinear (this); end if; end if !; end if
	
	end subroutine


!**********************************************************************************************
!      calculate velocity with linear interpolation at particle position
!**********************************************************************************************
      function vel_point_linear (this) result (q)
      use gslib, only: interpolation_linear
	  implicit none

	  type(particle_cl), intent(in) :: this
      real*8                        :: q(3)

      q(1) = interpolation_linear (this%cell%qxface,this%cell%coord(1))
      q(2) = interpolation_linear (this%cell%qyface,this%cell%coord(2))
      q(3) = interpolation_linear (this%cell%qzface,this%cell%coord(3))

      end function
      
!**********************************************************************************************
!      calculate velocity with C. Zheng's (1993) interpolation scheme if particle located in a well cell
!**********************************************************************************************
      function velo_in_wellcell_ (well,geo,this) result (q)
      use well_vect_class
      use geometry_class
      use array_class,      only: value_array_
	  implicit none
      
      type(well_vect_cl), intent(in)    :: well
	  type(particle_cl),  intent(in)    :: this
      type(geometry_cl),  intent(in)    :: geo
      real*8                :: q(3)
      real*8                :: pi, Qw
      real*8                :: dz, anis, xp(3), qxface(2), qyface(2), qzface(2)
      integer               :: iwll, ilay
      
      pi = acos(-1.0_8)
      
      dz = value_array_ (geo%dz, 1, 1, this%cell%num(3))
      xp = this%position%xp
      qxface = this%cell%qxface
      qyface = this%cell%qyface
      qzface = this%cell%qzface
      
      anis = 1.d0
      iwll = this%control%wellID_in_partcell
      ilay = this%cell%num(3)
      Qw   = well%num(iwll)%Qw(ilay)
      
      q(1) = (Qw*sqrt(anis))/(2*pi*dz) * (xp(1)-well%num(iwll)%xwell)/((xp(1)-well%num(iwll)%xwell)**2/anis + (xp(2)-well%num(iwll)%ywell)**2) + (qxface(1)+qxface(2))/2
      q(2) = (Qw*sqrt(anis))/(2*pi*dz) * (xp(2)-well%num(iwll)%ywell)/((xp(1)-well%num(iwll)%xwell)**2/anis + (xp(2)-well%num(iwll)%ywell)**2) + (qyface(2)+qyface(1))/2
      q(3) = ((qzface(1)-qzface(2))/dz) * (xp(3)-dz/2) + qzface(1)
   
      end function
      
!*********************************************************************************************
!      calculate velocity with trilinear interpolation at particle position
!*********************************************************************************************
      function vel_point_trilinear (this) result (q)
	  use gslib, only: interpolation_trilinear
	  implicit none
	  type(particle_cl), intent(in)  :: this
	  real*8                         :: q(3)
      real*8                         :: xx,yy,zz

	  xx = this%cell%coord(1)
	  yy = this%cell%coord(2)
	  zz = this%cell%coord(3)

      q(1) = interpolation_trilinear (this%cell%qxnode,xx,yy,zz)
      q(2) = interpolation_trilinear (this%cell%qynode,xx,yy,zz)
      q(3) = interpolation_trilinear (this%cell%qznode,xx,yy,zz)

      end function
      
!**********************************************************************************************
!       update dispersion at nodes of the particle cell
!**********************************************************************************************
   subroutine update_dispersion_nodes_particle_ (this,geo,dispersion)
      use geometry_class
	  use dispersion_class
	  use heterogeneity_flags
	  !use loops_particles, only: nmove
	  implicit none
	  type(particle_cl),   intent(inout) :: this
	  type(geometry_cl),   intent(in)    :: geo
	  type(dispersion_cl), intent(in)    :: dispersion
	  real*8                             :: Dxxnode(2,2,2),Dyynode(2,2,2),Dzznode(2,2,2)
	  real*8                             :: Dxynode(2,2,2),Dxznode(2,2,2),Dyznode(2,2,2)
	  real*8                             :: qxnode(2,2,2),qynode(2,2,2),qznode(2,2,2)
      real*8                             :: coord(3)
      real*8                             :: al,ath,atv,dm,rpt,poro 
	  integer                            :: i,j,k


	   if (.not. this%control%switchcell) return
	   if (.not. dispersion%action) return
       !if (vel_homogeneous== .TRUE. .and. disp_homogeneous== .TRUE. .and. nmove /= 1) return
       
	   Dxxnode = 0.d0
	   Dyynode = 0.d0
	   Dzznode = 0.d0
       Dxynode = 0.d0
	   Dxznode = 0.d0
	   Dyznode = 0.d0

	   qxnode = this % cell % qxnode
       qynode = this % cell % qynode
	   qznode = this % cell % qznode 

	   dm     = this % prop % dm

            do i=1,2
	          do j=1,2
	            do k=1,2

                   coord(1) = dfloat(i) - 1.d0 
				   coord(2) = dfloat(i) - 1.d0
				   coord(3) = dfloat(i) - 1.d0

 	               al     = update_property (this,geo,dispersion%aL,coord)
	               ath    = update_property (this,geo,dispersion%aTH,coord)
	               atv    = update_property (this,geo,dispersion%aTV,coord)

                   call disp_tensor( qxnode(i,j,k) ,qynode(i,j,k) ,qznode(i,j,k) ,al,ath,atv,dm,   &
	                                 Dxxnode(i,j,k),Dyynode(i,j,k),Dzznode(i,j,k),              &
							         Dxynode(i,j,k),Dxznode(i,j,k),Dyznode(i,j,k) )
	            end do
	          end do
	        end do


   	   this%cell%Dxxnode = Dxxnode 
	   this%cell%Dyynode = Dyynode
	   this%cell%Dzznode = Dzznode
       this%cell%Dxynode = Dxynode
	   this%cell%Dxznode = Dxznode
	   this%cell%Dyznode = Dyznode

   end subroutine

!*************************************************************************************************************************
!     subroutine to calculate dispersion tensor at a point
!*************************************************************************************************************************

      subroutine disp_tensor (qx,qy,qz,al,ath,atv,dm,Dxx,Dyy,Dzz,Dxy,Dxz,Dyz)
	  implicit none

      real*8, intent(in)  :: qx,qy,qz,al,ath,atv,dm
	  real*8, intent(out) :: Dxx,Dyy,Dzz,Dxy,Dxz,Dyz
      real*8 :: uu

            uu = dsqrt(qx*qx + qy*qy + qz*qz)

            if(uu.eq.0.d0) uu = 10.d-21
        
            Dxx = al  * qx * qx / uu  +				   &
                  ath * qy * qy / uu  +		           &
                  atv * qz * qz / uu  + dm

            Dyy = ath * qx * qx / uu + 			       &
                  al  * qy * qy / uu +				   &
                  atv * qz * qz / uu + dm

            Dzz = al  * qz * qz / uu +				   &
                  atv * qx * qx / uu +				   &
                  atv * qy * qy / uu + dm
     
            Dxy = (al-ath) * qx * qy / uu
            Dxz = (al-atv) * qx * qz / uu
            Dyz = (al-atv) * qy * qz / uu

      end subroutine disp_tensor

      function velocity_centroid_ (this) result (v)
	     use gslib, only: norm
		 implicit none
		 type(particle_cl), intent(in) :: this
		 real*8                        :: v,vc(3),vx1,vx2,vy1,vy2,vz1,vz2
		     vx1 = this % cell % qxface(1) / this % prop % poro
	         vx2 = this % cell % qxface(2) / this % prop % poro
	         vy1 = this % cell % qyface(1) / this % prop % poro
	         vy2 = this % cell % qyface(2) / this % prop % poro
	         vz1 = this % cell % qzface(1) / this % prop % poro
	         vz2 = this % cell % qzface(2) / this % prop % poro    
              vc(1) = 0.5d0 * ( vx1 + vx2 )
			  vc(2) = 0.5d0 * ( vy1 + vy2 )
			  vc(3) = 0.5d0 * ( vz1 + vz2 )
              v = norm (vc)

	  end function

      function velocity_particle_ (this) result (v)
	     use gslib, only: norm
          implicit none
		  type(particle_cl), intent(in) :: this
          real*8                        :: v,vp(3)
			  vp(1) = this % vel % qpL(1) / this % prop % poro
	          vp(2) = this % vel % qpL(2) / this % prop % poro
	          vp(3) = this % vel % qpL(3) / this % prop % poro
                  v = norm (vp)
	  end function

!******************************************************************************************************************
!   allocate particle
!******************************************************************************************************************
	subroutine alloc_particle (this,mass_trans,nwell,nplane)
	   use mass_trans_class
	   use global_variables
	   use to_solve
	   implicit none

	   type(particle_cl), intent(inout) :: this
	   type(mass_trans_cl), intent(in)  :: mass_trans
	   integer, optional, intent(in)    :: nwell, nplane
	   integer                          :: nzones
	     this % position % xpref       = 0.d0
		 this % position % xp          = 0.d0
		 this % position % xpold       = 0.d0
		 this % position % tp          = 0.d0
		 this % position % tpold       = 0.d0
		 this % cell     % coord       = 0.d0
		 this % cell     % num         = 0
		 this % cell     % numold      = 0
		 this % cell     % qxnode      = 0.d0
		 this % cell     % qynode      = 0.d0
		 this % cell     % qznode      = 0.d0
		 this % cell     % Dxxnode     = 0.d0
		 this % cell     % Dyynode     = 0.d0
		 this % cell     % Dzznode     = 0.d0
		 this % cell     % Dxynode     = 0.d0
		 this % cell     % Dxznode     = 0.d0
		 this % cell     % Dyznode     = 0.d0
		 this % cell     % qxface      = 0.d0
		 this % cell     % qyface      = 0.d0
		 this % cell     % qzface      = 0.d0
		 this % control  % remove      = .FALSE.
		 this % control  % switchcell  = .TRUE.
		 this % control  % stuck       = .FALSE.
		 this % control  % npbounce    = 0
         this % control  % npdomextens = 0
         this%control%passplane = 0
		 this%control%passwell  = 0
		 this % vel  % qpL  = 0.d0
		 this % vel  % qpT  = 0.d0
		 this % prop % aL   = 0.d0
		 this % prop % aTH  = 0.d0
		 this % prop % aTV  = 0.d0
		 this % prop % dm   = 0.d0
		 this % prop % poro = 1.d0
		 this % prop % mp   = 1.d0
         this%control%inblock = 0

		 allocate (this % zone % num);      this % zone   % num    = 0
		 allocate (this % zone % prenum);   this % zone   % prenum = 0
		 allocate (this % specie % num);    this % specie % num    = 1
		 allocate (this % specie % prenum); this % specie % prenum = 1

         !** Mass Transfer
         if (mass_transACTION) then
            select case (mass_transTYPE)
            case ('MULTIRATE')
                allocate(this % mass_transfer % multirateMT % alphap(nzoneim))
                allocate(this % mass_transfer % multirateMT % poroim(nzoneim))

            case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                allocate(this % mass_transfer % diffusionMT % alphap)
                allocate(this % mass_transfer % diffusionMT % poroim)

            case ('POWER_LAW')
                allocate(this % mass_transfer % powerMT % btot)
                allocate(this % mass_transfer % powerMT % Amin)
                allocate(this % mass_transfer % powerMT % Amax)
                allocate(this % mass_transfer % powerMT % power)

            case ('LOGNORMAL_LAW')
                allocate(this % mass_transfer % lognormalMT % btot)
                allocate(this % mass_transfer % lognormalMT % meanlog)
                allocate(this % mass_transfer % lognormalMT % stdvlog)

            case ('COMPOSITE_MEDIA')
                !multirate
                allocate(this % mass_transfer % compositeMT % multirateComp % alphap(mass_trans%composite%nmrate))
                allocate(this % mass_transfer % compositeMT % multirateComp % poroim(mass_trans%composite%nmrate))
                !spherical
                allocate(this % mass_transfer % compositeMT % sphericalComp % alphap)
                allocate(this % mass_transfer % compositeMT % sphericalComp % poroim)
                !cylindrical
                allocate(this % mass_transfer % compositeMT % cylindricalComp % alphap)
                allocate(this % mass_transfer % compositeMT % cylindricalComp % poroim)
                !layered
                allocate(this % mass_transfer % compositeMT % layeredComp % alphap)
                allocate(this % mass_transfer % compositeMT % layeredComp % poroim)
            end select

            allocate(this % zone % alpha(nspecie,nzoneim))
            allocate(this % zone % beta (nspecie,nzoneim))

         end if

         !** Sorption
         if (sorptionACTION) then
            select case (sorptionTYPE)
            case ('LINEAR')
                allocate(this % sorption % linear_sorp % R(nspecie))
                if (mass_transACTION) then
                    select case (mass_transTYPE)
                    case ('MULTIRATE')
                        allocate(this % sorption % linear_sorp % Rim(nspecie,nzoneim))
                    case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        allocate(this % sorption % linear_sorp % Rim(nspecie,1))
                    case ('POWER_LAW')
                        allocate(this % sorption % linear_sorp % Rim(nspecie,1))
                    case ('LOGNORMAL_LAW')
                        allocate(this % sorption % linear_sorp % Rim(nspecie,1))
                    case ('COMPOSITE_MEDIA')
                        allocate(this % sorption % linear_sorp % Rim(nspecie,nzoneim))
                    end select
                end if
            case ('CHTM')
                allocate(this % sorption % CHTM % bd)
                allocate(this % sorption % CHTM % kf)
                allocate(this % sorption % CHTM % kb)
            end select
         end if

         !** Decay
         if (decayACTION) then
            select case (decayTYPE)
            case ('SERIAL','GENERIC')
                allocate(this % decay % k(nspedecay))
                allocate(this % decay % y(nspedecay,nspedecay))
                if (mass_transACTION) then
                    select case (mass_transTYPE)
                    case ('MULTIRATE')
                        allocate(this % decay % kim(nspedecay,nzoneim))
                    case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        allocate(this % decay % kim(nspedecay,1))
                    case ('POWER_LAW')
                        allocate(this % decay % kim(nspedecay,1))
                    case ('LOGNORMAL_LAW')
                        allocate(this % decay % kim(nspedecay,1))
                    case ('COMPOSITE_MEDIA')
                        allocate(this % decay % kim(nspedecay,nzoneim))
                    end select
                end if
            case ('SERIAL_MOMENTS')
                allocate(this % decay % k(nspedecay))
                allocate(this % decay % y(nspedecay,nspedecay))
                allocate(this % decay % serial_mom % P(nspedecay,nspedecay))
                allocate(this % decay % serial_mom % S(nspedecay,nspedecay))
                allocate(this % decay % serial_mom % Sinv(nspedecay,nspedecay))
                allocate(this % decay % serial_mom % Areac(nspedecay,nspedecay,3))
                allocate(this % decay % serial_mom % B2(nspedecay,nspedecay,3))
                allocate(this % decay % serial_mom % Reff(nspedecay,nspedecay))
                if (mass_transACTION) then
                    stop 'CANNOT SOLVE SERIAL MOMENT AND MASS TRANSFER'
                end if
            end select
         end if

        this % block % switchblock = .TRUE.
        this % block % numold(3)   = 0   
		this % block % num(3)      = 0      
		this % block % Ti          = 0.d0
		this % block % Xpi         = 0.d0
	end subroutine

	subroutine allocate_particle_min (this)
	   use global_variables
	   implicit none
	   type(particle_cl), intent(inout) :: this
	   integer                          :: nzones
		 this % position % xp          = 0.d0
		 this % position % tp          = 0.d0
		 this % cell     % coord       = 0.d0
		 this % cell     % num         = 0
		 this % cell     % numold      = 0
		 this % cell     % qxnode      = 0.d0
		 this % cell     % qynode      = 0.d0
		 this % cell     % qznode      = 0.d0
		 this % cell     % Dxxnode     = 0.d0
		 this % cell     % Dyynode     = 0.d0
		 this % cell     % Dzznode     = 0.d0
		 this % cell     % Dxynode     = 0.d0
		 this % cell     % Dxznode     = 0.d0
		 this % cell     % Dyznode     = 0.d0
		 this % cell     % qxface      = 0.d0
		 this % cell     % qyface      = 0.d0
		 this % cell     % qzface      = 0.d0
		 this % control  % remove      = .FALSE.
		 this % control  % switchcell  = .TRUE.
		 this % control  % stuck       = .FALSE.
		 this % control  % npbounce    = 0
         this % control  % npdomextens = 0
		 this % vel  % qpL  = 0.d0
		 this % vel  % qpT  = 0.d0
		 this % prop % aL   = 0.d0
		 this % prop % aTH  = 0.d0
		 this % prop % aTV  = 0.d0
		 this % prop % dm   = 0.d0
		 this % prop % poro = 1.d0
		 this % prop % mp   = 1.d0
		 this % prop % rp    = 1.d0
		 if (.not.associated(this % zone   % num)) allocate (this % zone % num)
		 if (.not.associated(this % specie % num)) allocate (this % specie % num)
		 this % zone % num   = 0
		 this % specie % num = 1
         this % block % switchblock = .TRUE.
         this % block % numold   = 0   
		 this % block % num      = 0      
	end subroutine

!****************************************************************************************
!      set particle to zero
!****************************************************************************************
    subroutine set_particle_to_zero_ (this)
	   use global_variables
	   implicit none
	   type(particle_cl), intent(inout) :: this
	   integer                          :: nzones
	     this % position % xpref       = 0.d0
		 this % position % xp          = 0.d0
		 this % position % xpold       = 0.d0
		 this % position % tp          = 0.d0
		 this % position % tpold       = 0.d0
		 this % cell     % coord       = 0.d0
		 this % cell     % num         = 0
		 this % cell     % numold      = 0
		 this % cell     % qxnode      = 0.d0
		 this % cell     % qynode      = 0.d0
		 this % cell     % qznode      = 0.d0
		 this % cell     % Dxxnode     = 0.d0
		 this % cell     % Dyynode     = 0.d0
		 this % cell     % Dzznode     = 0.d0
		 this % cell     % Dxynode     = 0.d0
		 this % cell     % Dxznode     = 0.d0
		 this % cell     % Dyznode     = 0.d0
		 this % cell     % qxface      = 0.d0
		 this % cell     % qyface      = 0.d0
		 this % cell     % qzface      = 0.d0
		 this % control  % remove      = .FALSE.
		 this % control  % switchcell  = .TRUE.
		 this % control  % stuck       = .FALSE.
		 this % control  % npbounce    = 0
         this % control  % npdomextens = 0
		 this % control  % passplane   = .FALSE.
         this % control  % passwell    = .FALSE.
		 this % vel  % qpL  = 0.d0
		 this % vel  % qpT  = 0.d0
		 this % prop % aL   = 0.d0
		 this % prop % aTH  = 0.d0
		 this % prop % aTV  = 0.d0
		 this % prop % dm   = 0.d0
		 this % prop % poro = 1.d0
		 this % prop % mp   = 1.d0
		 this % prop % mp   = 1.d0
         this % prop % rp   = 1.d0
         this % control  % inblock   = 0
		 if (associated(this % zone % num))          this % zone % num      = 0
		 if (associated(this % zone % prenum))       this % zone % prenum   = 0
		 if (associated(this % specie % num))        this % specie % num    = 1
		 if (associated(this % specie % prenum))     this % specie % prenum = 1
		 if (associated(this % zone % alpha))        this % zone % alpha    = 0.d0
	     if (associated(this % zone % beta ))        this % zone % beta     = 0.d0

         !Decay
         if (associated(this%decay%k))                this%decay%k          = 0.d0
         if (associated(this%decay%y))                this%decay%y          = 0.d0
         if (associated(this%decay%kim))              this%decay%kim        = 0.d0
         if (associated(this%decay%serial_mom%P))     this%decay%serial_mom%P     = 0.d0
         if (associated(this%decay%serial_mom%S))     this%decay%serial_mom%S     = 0.d0
         if (associated(this%decay%serial_mom%Sinv))  this%decay%serial_mom%Sinv  = 0.d0
         if (associated(this%decay%serial_mom%Areac)) this%decay%serial_mom%Areac = 0.d0
         if (associated(this%decay%serial_mom%B2))    this%decay%serial_mom%B2    = 0.d0
         if (associated(this%decay%serial_mom%Reff))  this%decay%serial_mom%Reff    = 0.d0
         
         !Mass Transfer
         !* multirate
		 if (associated(this%mass_transfer%multirateMT%alphap))     this%mass_transfer%multirateMT%alphap = 0.d0
         if (associated(this%mass_transfer%multirateMT%poroim))     this%mass_transfer%multirateMT%poroim = 0.d0
         !* diffusion model
         if (associated(this%mass_transfer%diffusionMT%alphap))     this%mass_transfer%diffusionMT%alphap = 0.d0
         if (associated(this%mass_transfer%diffusionMT%poroim))     this%mass_transfer%diffusionMT%poroim = 0.d0
         !* power model
         if (associated(this%mass_transfer%powerMT%btot))           this%mass_transfer%powerMT%btot       = 0.d0
         if (associated(this%mass_transfer%powerMT%Amin))           this%mass_transfer%powerMT%Amin       = 0.d0
         if (associated(this%mass_transfer%powerMT%Amax))           this%mass_transfer%powerMT%Amax       = 0.d0
         if (associated(this%mass_transfer%powerMT%power))          this%mass_transfer%powerMT%power      = 0.d0
         !* LogNormal model
         if (associated(this%mass_transfer%lognormalMT%btot))       this%mass_transfer%lognormalMT%btot    = 0.d0
         if (associated(this%mass_transfer%lognormalMT%meanlog))    this%mass_transfer%lognormalMT%meanlog = 0.d0
         if (associated(this%mass_transfer%lognormalMT%stdvlog))    this%mass_transfer%lognormalMT%stdvlog = 0.d0
         !* Composite model

         !Sorption
         if (associated(this%sorption%linear_sorp%R))        this%sorption%linear_sorp%R   = 0.d0
         if (associated(this%sorption%linear_sorp%Rim))      this%sorption%linear_sorp%Rim = 0.d0
		 if (associated(this%sorption%CHTM%bd))              this%sorption%CHTM%bd   = 0.d0
	     if (associated(this%sorption%CHTM%kf))              this%sorption%CHTM%kf   = 0.d0
	     if (associated(this%sorption%CHTM%kb))              this%sorption%CHTM%kb   = 0.d0

         this % block % switchblock = .TRUE.
         this % block % numold(3)   = 0   
		 this % block % num(3)      = 0      
		 this % block % Ti          = 0.d0
         this % block % Xpi         = 0.d0
	end subroutine

!****************************************************************************************
!      check if particle can move by advection or molecular diffusion 
!****************************************************************************************
    subroutine check_mobility_particle_ (this, advection, dispersion)
	   use advection_class
	   use dispersion_class
	   use array_class
	   use gslib, only: open_fname
	   use global_variables, only: fdbg,npStuck
	   use loops_particles,  only: ip,nmove
	   implicit none
	     type(particle_cl),   intent(inout) :: this
         type(advection_cl),  intent(in)    :: advection
		 type(dispersion_cl), intent(in)    :: dispersion
		 real*8                             :: Dm,DmTH,DmTV
		 integer                            :: unit  
		   if (dispersion%action) then
		      Dm   = value_array_ (dispersion%dm,this%cell%num)
		      DmTH = value_array_ (dispersion%dmTH,this%cell%num)
		      DmTV = value_array_ (dispersion%dmTV,this%cell%num) 
		   else
		      Dm   = 0.d0
		      DmTH = 0.d0
		      DmTV = 0.d0
		   end if
		   if ( this%vel%qpL(1) == 0.d0 .and. &
		        this%vel%qpL(2) == 0.d0 .and. &
				this%vel%qpL(3) == 0.d0 .and. &
				Dm == 0.d0 .and. DmTH == 0.d0 .and. DmTV == 0.d0 ) then  !particle cannot move by any mechanism
				this%control%stuck = .TRUE.
				npStuck = npStuck + 1
				write(*,1), ' Stuck: ', ' Part = ',ip,                &
				                        ' Move = ',nmove,             &
				                        ' X = ',this%position%xp(1),  &
				                        ' Y = ',this%position%xp(2),  &
							  		    ' Z = ',this%position%xp(3)
				call open_fname (fdbg,unit)
				write(unit,1), ' Stuck: ', ' Part = ',ip,                &
				                           ' Move = ',nmove,             &
				                           ' X = ',this%position%xp(1),  &
				                           ' Y = ',this%position%xp(2),  &
							  		       ' Z = ',this%position%xp(3)
	 1          format (a8,2(a8,i7),3(a5,g10.3))
	            close(unit)
		   end if 
	end subroutine

!*************************************************************************************
!   print the passage time of a particle to a file
!*************************************************************************************
    subroutine print_particle_passtime_ (this,fname,string,isurface,ip,time)
       use gslib, only: generate_unit
	   use loops_particles, only: kinj,nmove
	   implicit none
	   type(particle_cl), intent(inout) :: this
	   character(len=*),  intent(in)    :: fname,string
	   integer,           intent(in)    :: isurface,ip
	   real*8,            intent(in)    :: time
	   integer                          :: iunit
	   logical                          :: connected
	       inquire(file=fname,opened=connected)
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
			   iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown')
	           write(iunit,*) 'TITLE="TECPLOT input file" '
	           write(iunit,*) 'VARIABLE="CONTROL" "PART" "Tp" "MASS" "INJ" "MOVE" '
		   end if  
		   write(iunit,2) string,isurface,ip,time,this%prop%mp,kinj,nmove
         2 format(a8,x,i3,3x,i8,6x,x,2(e15.8,2x),2(i7,8x))
	end subroutine

    subroutine print_particle_passtime_xyz_ (this,fname,string,isurface,ip,time,xp,yp,zp)
       use gslib, only: generate_unit
	   implicit none
	   type(particle_cl), intent(inout) :: this
	   character(len=*),  intent(in)    :: fname,string
	   integer,           intent(in)    :: isurface,ip
	   real*8,            intent(in)    :: time
	   real*8,            intent(in)    :: xp,yp,zp
	   integer                          :: iunit
	   logical                          :: connected
	       inquire(file=fname,opened=connected)
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
			   iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown')
			   write(iunit,1) 
			   write(iunit,1) ' CONTROL NUM   ', &
			                  '   PARTICLE    ', &
			                  'ARRIVAL TIME   ', &
							  '    DXp        ', &
							  '    DYp        ', &
							  '    DZp        ', &
							  ' PART MASS     ', &
							  '    ZONE       ', &
							  '   SPECIE      ' 
             1 format(9(a15,x))
	           write(iunit,1)
		   end if  	        
		   write(iunit,2) string,isurface,ip,time,xp,yp,zp,this%prop%mp,this%zone%num,this%specie%num
         2 format(a8,x,i3,3x,i8,6x,x,5(e15.8,2x),2(i7,8x))

	end subroutine


!***************************************************************************
!    update properties for mass transfer
!***************************************************************************
    subroutine update_mass_trans_ ( this, geo, mass_trans )
    use geometry_class
    use mass_trans_class
    use loops_particles,  only: nmove
    use global_variables, only: nzoneim
    use heterogeneity_flags
    implicit none
    
    type(particle_cl),   intent(inout)  :: this
    type(mass_trans_cl), intent(in)     :: mass_trans
    type(geometry_cl),   intent(in)     :: geo
    logical                             :: heterogeneous
    integer                             :: izone

    if (.not. mass_trans % action) return
    heterogeneous = (.not.mass_trans_homogeneous).or.(.not.poro_homogeneous)
    if ( (.not.heterogeneous .and. nmove == 1 ) .or.   & 
	     ( heterogeneous     .and. this%control%switchcell ) ) then

	   select case (mass_trans % type_mass_trans)
	   case ('MULTIRATE')
            do izone=1,nzoneim
			this % mass_transfer % multirateMT % poroim(izone) = update_property_block_discrete (this,geo,mass_trans % multirate % poro_imm(izone))
			this % mass_transfer % multirateMT % alphap(izone) = update_property_block_discrete (this,geo,mass_trans % multirate % alphap(izone))
            enddo

	   case ('SPHERICAL_DIFFUSION')
			this % mass_transfer % diffusionMT % poroim = update_property_block_discrete (this,geo,mass_trans % sphericalDiff % poro_imm)
			this % mass_transfer % diffusionMT % alphap = update_property_block_discrete (this,geo,mass_trans % sphericalDiff % alphap)
	   
	   case ('LAYERED_DIFFUSION')
			this % mass_transfer % diffusionMT % poroim = update_property_block_discrete (this,geo,mass_trans % layeredDiff % poro_imm)
			this % mass_transfer % diffusionMT % alphap = update_property_block_discrete (this,geo,mass_trans % layeredDiff % alphap)
	   
	   case ('CYLINDRICAL_DIFFUSION')
			this % mass_transfer % diffusionMT % poroim = update_property_block_discrete (this,geo,mass_trans % cylindricalDiff % poro_imm)
			this % mass_transfer % diffusionMT % alphap = update_property_block_discrete (this,geo,mass_trans % cylindricalDiff % alphap)
	   
	   case ('POWER_LAW')
			this % mass_transfer % powerMT % btot  = update_property_block_discrete (this,geo,mass_trans % power % btot)
			this % mass_transfer % powerMT % Amin  = update_property_block_discrete (this,geo,mass_trans % power % Amin)
			this % mass_transfer % powerMT % Amax  = update_property_block_discrete (this,geo,mass_trans % power % Amax)
			this % mass_transfer % powerMT % power = update_property_block_discrete (this,geo,mass_trans % power % power)

	   case ('LOGNORMAL_LAW')
			this % mass_transfer % lognormalMT % btot     = update_property_block_discrete (this,geo,mass_trans % LogNormal % btot)
			this % mass_transfer % lognormalMT % meanlog  = update_property_block_discrete (this,geo,mass_trans % LogNormal % mean)
			this % mass_transfer % lognormalMT % stdvlog  = update_property_block_discrete (this,geo,mass_trans % LogNormal % stdv)

	   end select
       
       call get_alpha_beta_ ( this )
       
    end if
    
    end subroutine

!   GET MASS TRANSFER COEFFICIENT MATRICES ALPHA AND BETA
!***************************************************************************************************
    subroutine get_alpha_beta_ (particle)
        use constants, only: pi
        use array_class
        use gslib, only: BesselRoot
        use global_variables, only: nspecie, nzoneim
        use to_solve
        implicit none

        type(particle_cl),  intent(inout)       :: particle

	    integer                                 :: ispe,izone,i
        real*8                                  :: alpha(nspecie,nzoneim),beta(nspecie,nzoneim)
        real*8                                  :: poroim(nzoneim),alphap(nzoneim),Rim(nspecie,nzoneim),a(nzoneim)
        real*8                                  :: R(nspecie),alpha_sph(nspecie),beta_sph(nspecie)
	    real*8                                  :: sum_1,sum_2,Amax,Amin,power,btot,Aharm,meanlog,stdvlog,mean,stdv
	    real*8, allocatable                     :: ro(:)

        R   = 1.d0
        Rim = 1.d0
    !**********************************************************
    !Define alpha and beta
        select case (mass_transTYPE)

    !**********************************************************
    !Mass Transfer model = multirate

        case ('MULTIRATE')

             !Linear sorption parameters
	         do ispe=1,nspecie
	            if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
                   R(ispe)     = particle % sorption % linear_sorp % R(ispe)
                   do izone=1,nzoneim
                     Rim(ispe,izone) = particle % sorption % linear_sorp % Rim(ispe,izone)
                   end do
	            end if
             end do
             
             !Alpha & Beta
	         do izone=1,nzoneim
	            alphap(izone)   = particle % mass_transfer % multirateMT % alphap(izone)
	            poroim(izone)   = particle % mass_transfer % multirateMT % poroim(izone)
	         do ispe=1,nspecie
	            beta(ispe,izone)  = (poroim(izone)*Rim(ispe,izone))/(particle % prop % poro*R(ispe))
	            alpha(ispe,izone) = alphap(izone)/Rim(ispe,izone)
	         enddo; enddo

    !**********************************************************
    !Mass Transfer model = spherical diffusion
        case ('SPHERICAL_DIFFUSION')

             !Linear sorption parameters
	         do ispe=1,nspecie
	            if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
                   R(ispe)     = particle % sorption % linear_sorp % R(ispe)
                   do izone=1,nzoneim
                     Rim(ispe,izone) = particle % sorption % linear_sorp % Rim(ispe,izone)
                   end do
	            end if
             end do
             
             !Alpha & Beta
	         do ispe=1,nspecie
	            alpha_sph(ispe) =  particle % mass_transfer % diffusionMT % alphap / Rim(ispe,1)
	            beta_sph(ispe)  = (particle % mass_transfer % diffusionMT % poroim * Rim(ispe,1)) / (particle % prop % poro * R(ispe))
	         enddo

	         allocate (ro(nzoneim-1))
    !        .......first Terms in Multirate Series
	         do ispe=1,nspecie
	         do izone=1,nzoneim-1
		        ro(izone)         = dble(izone) * pi
		        alpha(ispe,izone) = alpha_sph(ispe) * (ro(izone)**2)
	            beta(ispe,izone)  = beta_sph(ispe)  * 6.d0 / (ro(izone)**2)
	         enddo; enddo

    !        .......final Term in Multirate Series
	         sum_1 = sum ( 1.d0 / (ro**2) )
	         sum_2 = sum ( 1.d0 / (ro**4) )
	         sum_1 = 1.d0 -  6.d0 * sum_1
	         sum_2 = 1.d0 - 90.d0 * sum_2
	         do ispe=1,nspecie
                alpha(ispe,nzoneim) = alpha_sph(ispe) * sum_1 / sum_2 * 15.d0
		        beta(ispe,nzoneim)  = beta_sph(ispe)  * sum_1; enddo
	         deallocate (ro)

    !**********************************************************
    !Mass Transfer model = layered diffusion
        case ('LAYERED_DIFFUSION')

             !Linear sorption parameters
	         do ispe=1,nspecie
	            if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
                   R(ispe)     = particle % sorption % linear_sorp % R(ispe)
                   do izone=1,nzoneim
                     Rim(ispe,izone) = particle % sorption % linear_sorp % Rim(ispe,izone)
                   end do
	            end if
             end do
             
	         do ispe=1,nspecie
	            alpha_sph(ispe) =  particle % mass_transfer % diffusionMT % alphap / Rim(ispe,1)
	            beta_sph(ispe)  = (particle % mass_transfer % diffusionMT % poroim * Rim(ispe,1)) / (particle % prop % poro * R(ispe))
	         enddo

	         allocate (ro(nzoneim-1))
    !        .......first Terms in Multirate Series
	         do ispe=1,nspecie
	         do izone=1,nzoneim-1
		        ro(izone)    = dble(2*izone-1) * pi
		        alpha(ispe,izone) = alpha_sph(ispe) * (ro(izone)**2) / 4.d0
	            beta(ispe,izone)  = beta_sph(ispe)  * 8.d0 / (ro(izone)**2)
	         enddo; enddo

    !        .......final Term in Multirate Series
	         sum_1 = sum ( 1.d0 / (ro**2) )
	         sum_2 = sum ( 1.d0 / (ro**4) )
	         sum_1 = 1.d0 -  8.d0 * sum_1
	         sum_2 = 1.d0 - 96.d0 * sum_2
	         do ispe=1,nspecie
                alpha(ispe,nzoneim) = alpha_sph(ispe) * sum_1 / sum_2 * 3.d0
		        beta(ispe,nzoneim)  = beta_sph(ispe)  * sum_1; enddo
	         deallocate (ro)

    !**********************************************************
    !Mass Transfer model = cylindral diffusion
        case ('CYLINDRICAL_DIFFUSION')

             !Linear sorption parameters
	         do ispe=1,nspecie
	            if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
                   R(ispe)     = particle % sorption % linear_sorp % R(ispe)
                   do izone=1,nzoneim
                     Rim(ispe,izone) = particle % sorption % linear_sorp % Rim(ispe,izone)
                   end do
	            end if
             end do
             
	         do ispe=1,nspecie
	            alpha_sph(ispe) =  particle % mass_transfer % diffusionMT % alphap / Rim(ispe,1)
	            beta_sph(ispe)  = (particle % mass_transfer % diffusionMT % poroim * Rim(ispe,1)) / (particle % prop % poro * R(ispe))
	         enddo

             allocate (ro(nzoneim-1))
    !        .......first Terms in Multirate Series
             !ro = square of Bessel roots 
             do izone=1,nzoneim-1
                ro(izone) = BesselRoot(izone)
                ro(izone) = ro(izone) * ro(izone)
             end do

             do ispe=1,nspecie
             do izone=1,nzoneim-1
                alpha(ispe,izone) = alpha_sph(ispe) * (ro(izone)**2)
	            beta (ispe,izone) = beta_sph(ispe)  * 4.d0 / (ro(izone)**2)
	         enddo; enddo

    !        .......final Term in Multirate Series
             sum_1 = sum ( 1.d0/(ro**2) )
             sum_2 = sum ( 1.d0/(ro**4) )
             sum_1 = 1.d0 -  4.d0 * sum_1
             sum_2 = 1.d0 - 32.d0 * sum_2
             do ispe=1,nspecie
                alpha(ispe,nzoneim) = alpha_sph(ispe) * sum_1 / sum_2 * 8.d0
                beta(ispe,nzoneim)  = beta_sph(ispe)  * sum_1; enddo
             deallocate ( ro )

    !**********************************************************
    !Mass Transfer model = power law
        case ('POWER_LAW')

             !Linear sorption parameters
	         do ispe=1,nspecie
	            if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
                   R(ispe)     = particle % sorption % linear_sorp % R(ispe)
                   do izone=1,nzoneim
                     Rim(ispe,izone) = particle % sorption % linear_sorp % Rim(ispe,izone)
                   end do
	            end if
             end do
             
             Amax  = particle % mass_transfer % powerMT % Amax
             Amin  = particle % mass_transfer % powerMT % Amin;  if (Amax < Amin) stop ' Input Data Error Amax < Amin: truncated power mass transfer '
             power = particle % mass_transfer % powerMT % power; if (power > 4)   stop  ' truncated power mass transfer model : exponent > 4 '
             btot  = particle % mass_transfer % powerMT % btot

    !        .......first Terms in Multirate Series
             do izone = 1,nzoneim
                a(izone) = Amin + (Amax-Amin)/ (dble(nzoneim-1)**(4-power)) * (dble(izone-1)**(4-power))
             end do
             do ispe=1,nspecie
             do izone = 1,nzoneim-1
                alpha(ispe,izone) = 0.5d0*(a(izone)+a(izone+1))
                if (power/=2.d0) then
	               beta(ispe,izone) = btot*(power-2.d0)/((Amax**(power-2.d0))-(Amin**(power-2.d0)))*(alpha(ispe,izone)**(power-3.d0))*(a(izone+1)-a(izone))
                else if (power==2.d0) then
		           beta(ispe,izone)  = btot/dlog(Amax/Amin)/alpha(ispe,izone)*(a(izone+1)-a(izone))
		        endif
             end do; end do

    !        .......final Term in Multirate Series
	         !Armonic mean
	         if     (power==2.d0) then; Aharm = dlog((Amax/Amin))*Amin*(Amax/Amin)/((Amax/Amin)-1.d0)
	         else if(power==3.d0) then; Aharm = Amax/dlog((Amax/Amin))
	         else;                      Aharm = Amin*((power-3)*((Amax/Amin)**(power-2.d0))-1.d0)/((power-2)*((Amax/Amin)**(power-3.d0))-1.d0) 
	         end if

	         do ispe=1,nspecie
	            if (nzoneim == 1) then
		           alpha (ispe,nzoneim) = Aharm
		           beta  (ispe,nzoneim) = btot
	            else
		           alpha(ispe,nzoneim) = Amax
		           beta(ispe,nzoneim)  = 0.d0
		           sum_1=0.d0
		           sum_2=0.d0
		           do izone=1,nzoneim
		              sum_1 = sum_1+beta(ispe,izone)
		              sum_2 = sum_2+(beta(ispe,izone)/alpha(ispe,izone)); enddo
                   alpha(ispe,nzoneim) = Aharm*(1.d0-sum_1/btot)/(1.d0-Aharm/btot*sum_2)
	               beta(ispe,nzoneim)    = btot-sum_1
	            end if
	         end do

    !**********************************************************
    !Mass Transfer model = lognormal law
        case ('LOGNORMAL_LAW')

             !Linear sorption parameters
	         do ispe=1,nspecie
	            if (sorptionACTION .and. sorptionTYPE == 'LINEAR') then
                   R(ispe)     = particle % sorption % linear_sorp % R(ispe)
                   do izone=1,nzoneim
                     Rim(ispe,izone) = particle % sorption % linear_sorp % Rim(ispe,izone)
                   end do
	            end if
             end do
             
             btot    = particle % mass_transfer % lognormalMT % btot
             meanlog = particle % mass_transfer % lognormalMT % meanlog
             stdvlog = particle % mass_transfer % lognormalMT % stdvlog

    !        .......first Terms in Multirate Series
             mean = dexp(meanlog+0.5d0*(stdvlog**2))
             stdv = dexp(2.d0*meanlog+stdvlog**2)*(dexp(stdvlog**2)-1.d0)
             stdv = dsqrt(stdv) 
             Amin = mean-3.d0*stdv
             Amax = mean+3.d0*stdv 
             do izone = 1,nzoneim
                a(izone) = Amin + (Amax - Amin)/ dble(nzoneim-1) * dble(izone-1)
             end do
             a = dexp(a)
             do ispe=1,nspecie
             do izone=1,nzoneim-1
                alpha(ispe,izone) = 0.5d0*(a(izone)+a(izone+1))
                beta(ispe,izone)  = btot/(dsqrt(2.d0*pi)*stdvlog*alpha(ispe,izone))*dexp(-((dlog(alpha(ispe,izone))-meanlog)**2)/(2.d0*(stdvlog**2)))*(a(izone+1)-a(izone))
             end do; end do

    !        .......final Term in Multirate Series
             Aharm = dexp(meanlog-0.5d0*(stdvlog**2))
             do ispe=1,nspecie
                if (izone == 1) then
	               alpha(ispe,nzoneim) = Aharm
	               beta (ispe,nzoneim) = btot
                else
		            alpha(ispe,nzoneim) = Amax
		            beta (ispe,nzoneim) = 0.d0
		            do izone=1,nzoneim
		               sum_1 = sum_1+beta(ispe,izone)
		               sum_2 = sum_2+(beta(ispe,izone)/alpha(ispe,izone)); enddo
                    alpha(ispe,nzoneim) = Aharm*(1.d0-sum_1/btot)/(1.d0-Aharm/btot*sum_2)
		            beta(ispe,izone)  = btot-sum_1
		        end if
             end do

    !**********************************************************
        end select
        particle % zone % alpha = alpha
        particle % zone % beta  = beta
        
    end subroutine get_alpha_beta_


!***************************************************************************
!    update properties for decay network
!***************************************************************************
    subroutine update_sorption_ ( this, geo, reaction )
    use reaction_class
    use geometry_class
    use loops_particles,  only: nmove
    use global_variables, only: nzoneim, nspecie
    use to_solve,         only: mass_transACTION, mass_transTYPE
    use heterogeneity_flags
    implicit none

    type(particle_cl),   intent(inout)  :: this
    type(reaction_cl),   intent(in)     :: reaction
    type(geometry_cl),   intent(in)     :: geo
    logical                             :: heterogeneous
    integer                             :: ispe, izone, izoneim

    if (.not. reaction % sorption % action) return
    heterogeneous = (.not.sorption_homogeneous)
    if ( (.not.heterogeneous .and. nmove == 1 ) .or.   & 
	     ( heterogeneous     .and. this%control%switchcell ) ) then

	     select case (reaction % sorption % type_sorption)
	     case ('LINEAR')
	          do ispe=1,nspecie
	     	     this % sorption % linear_sorp % R(ispe)   = update_property_block_discrete (this, geo, reaction % sorption % LinearSorption % R(ispe))
		         if (mass_transACTION .AND. nzoneim>0) then
	                select case (mass_transTYPE)
	                case ('MULTIRATE')
			            do izoneim=1, nzoneim
                        this % sorption % linear_sorp % Rim(ispe,izoneim)   = update_property_block_discrete (this, geo, reaction % sorption % LinearSorption % Rim(ispe,izoneim)); enddo
	                case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        this % sorption % linear_sorp % Rim(ispe,1)   = update_property_block_discrete (this, geo, reaction % sorption % LinearSorption % Rim(ispe,1))
	                case ('POWER_LAW')
                        this % sorption % linear_sorp % Rim(ispe,1)   = update_property_block_discrete (this, geo, reaction % sorption % LinearSorption % Rim(ispe,1))
	                case ('LOGNORMAL_LAW')
                        this % sorption % linear_sorp % Rim(ispe,1)   = update_property_block_discrete (this, geo, reaction % sorption % LinearSorption % Rim(ispe,1))
	                case ('COMPOSITE_MEDIA')
			            do izoneim=1, nzoneim
                        this % sorption % linear_sorp % Rim(ispe,izoneim)   = update_property_block_discrete (this, geo, reaction % sorption % LinearSorption % Rim(ispe,izoneim)); enddo
                    end select
	             end if
		      end do

	     case ('CHTM')
	          this % sorption % CHTM % bd   = update_property_block_discrete (this, geo, reaction % sorption % CHTM % bd)
		      this % sorption % CHTM % kf   = update_property_block_discrete (this, geo, reaction % sorption % CHTM % kf)
		      this % sorption % CHTM % kb   = update_property_block_discrete (this, geo, reaction % sorption % CHTM % kb)
		      if (mass_transACTION .AND. nzoneim>0) then
                 stop 'CANNOT SOLVE MASS TRANSFER AND SERIAL RX WITH HIGHER MOMENTS'
	          end if


	     end select

    end if
    
    end subroutine

!***************************************************************************
!    update properties for decay network
!***************************************************************************
    subroutine update_decay_ ( this, geo, reaction )
    use reaction_class
    use geometry_class
    use loops_particles,  only: nmove
    use global_variables, only: nzoneim, nspedecay
    use to_solve,         only: mass_transACTION, mass_transTYPE
    use heterogeneity_flags
    implicit none

    type(particle_cl),   intent(inout)  :: this
    type(reaction_cl),   intent(in)     :: reaction
    type(geometry_cl),   intent(in)     :: geo
    logical                             :: heterogeneous
    integer                             :: ispe, jspe, izone, izoneim

    if (.not. reaction % decay % action) return
    heterogeneous = (.not.decay_homogeneous)
    if ( (.not.heterogeneous .and. nmove == 1 ) .or.   & 
	     ( heterogeneous     .and. this%control%switchcell ) ) then

	     select case (reaction % decay % type_decay)
	     case ('SERIAL')
	          do ispe=1,nspedecay
	     	     this % decay % k(ispe)   = update_property_block_discrete (this, geo, reaction % decay % serial % k(ispe))
		         if (ispe>1) this % decay % y(ispe,ispe-1)   = update_property_block_discrete (this, geo, reaction % decay % serial % y(ispe,ispe-1))
		         if (mass_transACTION .AND. nzoneim>0) then
	                select case (mass_transTYPE)
	                case ('MULTIRATE')
			            do izoneim=1, nzoneim
                        this % decay % kim(ispe,izoneim) = update_property_block_discrete (this, geo, reaction % decay % serial % kim(ispe,izoneim)); enddo
	                case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        this % decay % kim(ispe,1) = update_property_block_discrete (this, geo, reaction % decay % serial % kim(ispe,1))
	                case ('POWER_LAW')
                        this % decay % kim(ispe,1) = update_property_block_discrete (this, geo, reaction % decay % serial % kim(ispe,1))
	                case ('LOGNORMAL_LAW')
                        this % decay % kim(ispe,1) = update_property_block_discrete (this, geo, reaction % decay % serial % kim(ispe,1))
	                case ('COMPOSITE_MEDIA')
			            do izoneim=1, nzoneim
                        this % decay % kim(ispe,izoneim) = update_property_block_discrete (this, geo, reaction % decay % serial % kim(ispe,izoneim)); enddo
                    end select
	             end if
		      end do

	     case ('SERIAL_MOMENTS')
	          do ispe=1,nspedecay
	          	 this % decay % k(ispe)   = update_property_block_discrete (this, geo, reaction % decay % serial_mom % k(ispe))
		         if (ispe>1) this % decay % y(ispe,ispe-1)   = update_property_block_discrete (this, geo, reaction % decay % serial_mom % y(ispe,ispe-1))
		         if (mass_transACTION .AND. nzoneim>0) then
                    stop 'CANNOT SOLVE MASS TRANSFER AND SERIAL RX WITH HIGHER MOMENTS'
	             end if
		      end do

	     case ('GENERIC')
	          do ispe=1,nspedecay
	          	 this % decay % k(ispe)   = update_property_block_discrete (this, geo, reaction % decay % generic_net % k(ispe))
	          	 do jspe=1,nspedecay
		            this % decay % y(ispe,jspe)   = update_property_block_discrete (this, geo, reaction % decay % generic_net % y(ispe,jspe))
		         end do
		         if (mass_transACTION .AND. nzoneim>0) then
	                select case (mass_transTYPE)
	                case ('MULTIRATE')
			            do izoneim=1, nzoneim
                        this % decay % kim(ispe,izoneim) = update_property_block_discrete (this, geo, reaction % decay % generic_net % kim(ispe,izoneim)); enddo
	                case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        this % decay % kim(ispe,1) = update_property_block_discrete (this, geo, reaction % decay % generic_net % kim(ispe,1))
	                case ('POWER_LAW')
                        this % decay % kim(ispe,1) = update_property_block_discrete (this, geo, reaction % decay % generic_net % kim(ispe,1))
	                case ('LOGNORMAL_LAW')
                        this % decay % kim(ispe,1) = update_property_block_discrete (this, geo, reaction % decay % generic_net % kim(ispe,1))
	                case ('COMPOSITE_MEDIA')
			            do izoneim=1, nzoneim
                        this % decay % kim(ispe,izoneim) = update_property_block_discrete (this, geo, reaction % decay % generic_net % kim(ispe,izoneim)); enddo
                    end select
	             end if
		      end do

	     end select

    end if

    end subroutine

!*****************************************************************************************
  function calc_module_Darcy_velocity_particle_ (this) result (val)
    implicit none
	type(particle_cl), intent(in) :: this
	real*8                        :: val
	  val = dsqrt(this%vel%qpl(1)**2 +  this%vel%qpl(2)**2 + this%vel%qpl(3)**2)
  end function

!*****************************************************************************************
!*****************************************************************************************
  end module 
