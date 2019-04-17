

  module advection_class
  use array_class
  implicit none

  private
  public :: advection_cl   !class
  public ::                                   &
            delete_advection_               , &
			vel_nodes_                      , &
			vel_faces_                      , &
			read_flux_from_mf2k_            , &
			print_velocity_property_        , &
			inquire_homogeneity_vel_        , &
			inquire_homogeneity_poro_       , &
			read_fluxes_                    , &
			calculate_advection_timeshots_  , &
			update_advection_timeshot_      , &
			print_velocity_timeshots_       , &
            print_velocity_field_TECPLOT_   , &
            initialize_advection_timeshot_         

  type advection_cl
      logical           :: action
	  type(array_cl)    :: qx,qy,qz,poro
	  integer           :: nt = 1              !total number of velocity snapshots
	  integer           :: it = 1              !timeshot interval associated with qx,qy,qz
	  integer           :: kper = 1            !stress period index
	  integer           :: kstp = 1            !timestep index (inside stress period)                 
	  real*8, pointer   :: time(:) => null()   !time discretization associated with velocity snapshots
      logical           :: switch = .TRUE.     !flag to switch velocity
      logical           :: loop_per = .FALSE.  !flag to loop stress period until final simulation time
   end type

    interface read_flux_from_mf2k_
	     module procedure read_flux_from_mf2k_1,read_flux_from_mf2k_2; end interface

	 
  contains

   subroutine print_velocity_field_TECPLOT_ (advection,geo,fname)
       use gslib, only: open_fname
       use geometry_class
       implicit none
       type(advection_cl), intent(in) :: advection
       type(geometry_cl),  intent(in) :: geo 
       character(len=*),   intent(in) :: fname
       integer                        :: iunit,it,klay,jrow,icol,nx,ny,nz
       integer                        :: cell(3),  k
       real*8                         :: qxval,qyval,qzval,qmodval,xc(3) 
       real*8                         :: qxface(2),qyface(2),qzface(2),time
          call open_fname(fname,iunit)
          if (.not.associated(advection%time)) return
          if (.not.advection%switch) return
           nx = geo%nx
           ny = geo%ny
           nz = geo%nz
           time = 0.5d0 * ( advection%time(advection%it) + advection%time(advection%it-1) )
		   write(iunit,1)   'ZONE T="',advection%time(advection%it),'"',' I=',nx,' J=',ny,' K=',nz,' F= POINT '
		 1 format (a8,x,g15.6,a2,3(a3,i4),a10) 
           do klay=1,nz
		     do jrow=1,ny
		       do icol=1,nx
			   cell(1)=icol
			   cell(2)=jrow
			   cell(3)=klay
               call vel_faces_ (advection,cell,qxface,qyface,qzface) 
               xc = get_cell_centroid_ (geo,cell)
			   qxval = 0.5d0*(qxface(1)+qxface(2))
			   qyval = 0.5d0*(qyface(1)+qyface(2))
			   qzval = 0.5d0*(qzface(1)+qzface(2))
			   qmodval = dsqrt(qxval*qxval+qyval*qyval+qzval*qzval)
!               if (qmodval > 0.d0) then
!			     qmodval = log(qmodval)
!			   else 
!			     qmodval = 0.
!			   end if
			   write(iunit,'(3(g15.6,x),4(g20.12,x))') xc(1),xc(2),xc(3),qxval,qyval,qzval,qmodval
			   end do
			 end do
		   end do 

   end subroutine


   subroutine print_velocity_timeshots_ (fname,advection)
       use gslib, only: open_fname
       implicit none
       type(advection_cl), intent(in) :: advection
       character(len=*),   intent(in) :: fname
       integer                        :: iunit,it
          call open_fname(fname,iunit)
          if (.not.associated(advection%time)) return
          write(iunit,*)
          write(iunit,*) ' velocity time intervals:'
          write(iunit,*) 
          do it =1,advection%nt
             write(iunit,*) it,advection%time(it-1),advection%time(it)
          end do
          write(iunit,*)
          close(iunit)
   end subroutine


   subroutine calculate_advection_timeshots_ (this)
       use velocity_times
       implicit none    
       type(advection_cl), intent(inout) :: this
       integer                           :: nt,iper,nper,ntstep,i,it
       real*8                            :: tmult,tlen,t0,dtf       
         nper = velotime%nper
         !get nt: total number of time steps
         nt = 0
         do iper=1,nper
           nt = nt + velotime%period(iper)%NSTP
         end do
         this%nt = nt
         
         !get time discretization
         if(.not. associated(this%time)) allocate(this%time(0:nt))
         this%time(0)=0.d0              
         t0 = 0.d0
         it = 0
         do iper=1,nper
            tmult  = velotime%period(iper)%MULT
            tlen   = velotime%period(iper)%PERLEN
            ntstep = velotime%period(iper)%NSTP
            
            if (tmult > 1) then
              it = it + 1
              dtf = tlen*((tmult-1.d0)/((tmult**ntstep)-1.d0))
              this%time(it)=t0+dtf
	          do i=2,ntstep
	             it = it + 1
	             dtf=dtf*tmult
	             this%time(it)=this%time(it-1)+dtf
              end do
              
	        else if (tmult <= 1.) then
	          it = it + 1
	          dtf = tlen/float(ntstep)
	          this%time(it)=t0+dtf
		      do i=2,ntstep
		         it = it + 1
		         this%time(it)=this%time(it-1)+dtf
		      end do
            end if
            
	        t0 = this%time(it)
	     end do   
    end subroutine



    subroutine initialize_advection_timeshot_ (this,time)
       use velocity_times
       implicit none
       type(advection_cl), intent(inout) :: this
       real*8,               intent(in)  :: time
       integer                           :: nper,itf,iper,istp,it
       integer                           :: kper,kstp
        
         it   = this%it
        
         if (time >= this%time(it-1) .and. time < this%time(it)) return

         if (time >= this%time(this%nt)) stop 'Could not find a velocity for this time regime'

         ! estimate new time index, period(kper), timestep (kstp)         
         
         kper = this%kper
         kstp = this%kstp
         
         this%switch = .TRUE.
         
         itf = it
         loop: do iper=kper,velotime%nper
                  if (time < this%time(itf)) exit loop
                  do istp=kstp,velotime%period(iper)%NSTP
                     itf = itf+1
                     if (time >= this%time(itf-1)) then
                        if (istp < velotime%period(iper)%NSTP) exit loop
                        exit
                     end if
                  end do
               end do loop

         this%it   = itf
         this%kper = iper
         this%kstp = istp
    
    end subroutine



    subroutine update_advection_timeshot_ (this,time)
       use velocity_times
       !use StressTime_class
       use global_variables, only: mf2k 
       implicit none
       type(advection_cl),   intent(inout) :: this
       real*8,               intent(in)    :: time
       integer                           :: nper,itf,iper,istp,it
       integer                           :: kper,kstp
       
       real*8                            :: time_flow
       integer                           :: iloop
       
         
         velotime%restart_flux_from_mf2k = .FALSE.
         
         ! if loop cbb file, modify the time to identify iper
         if (time >= this%time(this%nt) .and. velotime%loop_per) then
             iloop = ceiling((time+1E-5)/this%time(this%nt))                   !get cbb file loop index
             velotime%time_flow = time - this%time(this%nt)*(iloop-1)   !get time to consider in cbb file
             if (iloop /= velotime%iloop_per) then                      !new loop index: restart MF2K velocity file
                 velotime%restart_flux_from_mf2k = .TRUE.
                 velotime%iloop_per = iloop
                 rewind(mf2k%unitCBC)   !restart CBC file
             end if
         else
             velotime%time_flow = time
         end if
         
         it   = this%it
         ! check if effective time within incremental times; if yes, nothing to update
         if (velotime%time_flow >= this%time(it-1) .and. velotime%time_flow < this%time(it)) then
               this%switch = .FALSE.
               return
         end if
         
         ! keep last stress period (SP) flow condition if time larger than last SP time and cbb not looped
         if (velotime%time_flow >= this%time(this%nt) .and. velotime%loop_per==.FALSE.) then 
               this%switch = .FALSE.
               return
         end if
         
         ! reinitiate time index, period (kper), timestep (kstp) if just rewind the budget file
         if (velotime%restart_flux_from_mf2k == .TRUE.) then
            this%kper = 1 !reinitiate stress period
            this%kstp = 1 !reinitiate sp time step
            this%it   = 1
            it   = 1
            istp = 1
         end if 
         
         ! change of velocity, so estimate new time index, period (kper), timestep (kstp)
         kper = this%kper
         kstp = this%kstp
         
         this%switch = .TRUE.
         
         itf = it
         loop: do iper=kper,velotime%nper
                  if (velotime%time_flow < this%time(itf)) exit loop
                  do istp=kstp,velotime%period(iper)%NSTP
                     itf = itf+1
                     if (velotime%time_flow >= this%time(itf-1)) then
                        if (istp < velotime%period(iper)%NSTP) exit loop
                        exit
                     end if
                  end do
               end do loop

         this%it   = itf
         this%kper = iper
         this%kstp = istp
         
         
         
         !if (time >= this%time(it-1) .and. time < this%time(it)) then
         !      this%switch = .FALSE.
         !      return
         !end if
         !
         !if (time >= this%time(this%nt)) then
         !      this%switch = .FALSE.
         !      return
         !end if
         !
         !! change of velocity, so estimate new time index, period(kper), timestep (kstp)         
         !
         !kper = this%kper
         !kstp = this%kstp
         !
         !this%switch = .TRUE.
         !
         !itf = it
         !loop: do iper=kper,velotime%nper
         !         if (time < this%time(itf)) exit loop
         !         do istp=kstp,velotime%period(iper)%NSTP
         !            itf = itf+1
         !            if (time >= this%time(itf-1)) then
         !               if (istp < velotime%period(iper)%NSTP) exit loop
         !               exit
         !            end if
         !         end do
         !      end do loop
         !
         !this%it   = itf
         !this%kper = iper
         !this%kstp = istp
    
    end subroutine




  subroutine read_fluxes_ (advection,geo) 
     use geometry_class
  	 use global_variables,      only: dataqx,dataqy,dataqz,mf2k
  	 use heterogeneity_flags,   only: vel_homogeneous
     !use code_options,          only: simul_mode
     implicit none
	 type(advection_cl), intent(inout)   :: advection
	 type(geometry_cl),  intent(in)      :: geo
     

     if (.not.advection%action) return 
     if (.not.advection%switch) return
!
!         Read flow from MODFLOW cell-to-cell budget (not compact):
!
	      if (dataqx%flag == 2) then 
		      !call read_flux_from_mf2k_ (dataqx%file,advection,geo)
		      call read_flux_from_mf2k_ (mf2k,advection,geo)
			  if (dataqx%const /= 1.d0) advection%qx%values = dataqx%const * advection%qx%values
              if (dataqy%const /= 1.d0) advection%qy%values = dataqy%const * advection%qy%values
              if (dataqz%const /= 1.d0) advection%qz%values = dataqz%const * advection%qz%values
		  else 
!
!         Read flow from GSLIB file:
!
		      advection%qx = read_array_ (dataqx%file,dataqx%ivar,dataqx%const,dataqx%flag,geo%nx+1,geo%ny,geo%nz)
			  advection%qy = read_array_ (dataqy%file,dataqy%ivar,dataqy%const,dataqy%flag,geo%nx,geo%ny+1,geo%nz)
			  advection%qz = read_array_ (dataqz%file,dataqz%ivar,dataqz%const,dataqz%flag,geo%nx,geo%ny,geo%nz+1)
		  
          end if
          
          !if (simul_mode==2) then !if backward simulation, inverse velocities
          !    advection%qx%values = -1*advection%qx%values
          !    advection%qy%values = -1*advection%qy%values
          !    advection%qz%values = -1*advection%qz%values
          !end if
          
		  vel_homogeneous   = inquire_homogeneity_vel_ (advection)

  end subroutine

  function inquire_homogeneity_vel_ (this) result (flag)
      use array_class
	  implicit none
	  type(advection_cl), intent(in)   :: this
	  logical                          :: flag
	  integer                          :: i
	        flag = .TRUE.
            if(associated(this % qx % values )  .and. length_array_ (this % qx )  > 1)  flag = .FALSE.
			if(associated(this % qy % values )  .and. length_array_ (this % qy )  > 1)  flag = .FALSE.
			if(associated(this % qz % values )  .and. length_array_ (this % qz )  > 1)  flag = .FALSE.  
   end function

  function inquire_homogeneity_poro_ (this) result (flag)
      use array_class
	  implicit none
	  type(advection_cl), intent(in)  :: this
	  logical                         :: flag
	  integer                         :: i
	        flag = .TRUE.
            if(associated(this % poro % values )  .and. length_array_ (this % poro)   > 1)  flag = .FALSE.
   end function


 subroutine delete_advection_ (this) ! destructor
    use array_class
	implicit none	 
    type(advection_cl), intent(inout) :: this
	   this%action = .FALSE.
	   call delete_array_ (this%qx)
	   call delete_array_ (this%qy)
	   call delete_array_ (this%qz)
	   call delete_array_ (this%poro)
 end subroutine

 subroutine print_velocity_property_ (this,fname) 
    use array_class
	use gslib, only: open_fname
	implicit none
	type(advection_cl), intent(in) :: this
	character(len=*),   intent(in) :: fname
    type(array_cl)                 :: a
	real*8                         :: vmin,vmax,va
	integer                        :: iunit
    call open_fname (fname,iunit)
	write(iunit,*)
	write(iunit,*) ' Properties of the Velocity Field:'
    write(iunit,*)
!   ...x direction
    a    = equal_to_ (this%qx)
	if(count(a%values/=0.d0) == 0) then
	   vmax = 0.d0
	   va   = 0.d0
	   vmin = 0.d0
	   else
	   a%values = dabs(a%values)
	   vmin = minval(a%values,mask=a%values>0.d0)
	   vmax = maxval(a%values,mask=a%values>0.d0)
       va   = sum(a%values,mask=a%values>0.d0)/dble(count(a%values>0.d0))
	end if
    write(iunit,1) ' qx_min: ',vmin,' qx_mean: ',va,' qx_max: ',vmax
!   ...y direction
    a    = equal_to_ (this%qy)
	if(count(a%values/=0.d0) == 0) then
	   vmax = 0.d0
	   va   = 0.d0
	   vmin = 0.d0
	   else
	   a%values = dabs(a%values)
	   vmin = minval(a%values,mask=a%values>0.d0)
	   vmax = maxval(a%values,mask=a%values>0.d0)
       va   = sum(a%values,mask=a%values>0.d0)/dble(count(a%values>0.d0))
	end if
    write(iunit,1) ' qy_min: ',vmin,' qy_mean: ',va,' qy_max: ',vmax
!   ...z direction
    a    = equal_to_ (this%qz)
	if(count(a%values/=0.d0) == 0) then
	   vmax = 0.d0
	   va   = 0.d0
	   vmin = 0.d0
	   else
	   a%values = dabs(a%values)
	   vmin = minval(a%values,mask=a%values>0.d0)
	   vmax = maxval(a%values,mask=a%values>0.d0)
       va   = sum(a%values,mask=a%values>0.d0)/dble(count(a%values>0.d0))
	end if
    write(iunit,1) ' qz_min: ',vmin,' qz_mean: ',va,' qz_max: ',vmax

 1  format (3(a10,x,g15.6,x))

    close(iunit)

 end subroutine 


!*****************************************************************************************
!    calculates x,y,z velocity components at each node 
!*****************************************************************************************       
    subroutine vel_nodes_ (this,geo,cell,qxnode,qynode,qznode)
    use array_class
	use geometry_class
	implicit none

	type (geometry_cl),  intent(in)  :: geo 
	type (advection_cl), intent(in)  :: this
	integer,             intent(in)  :: cell(3)
	real*8,              intent(out) :: qxnode(2,2,2),qynode(2,2,2),qznode(2,2,2)
	real*8                           :: wx,wy,wz,sum
    integer                          :: m,n,i,j,k,icell,jcell,kcell,nx,ny,nz,loci,locj,lock

    loci = cell(1)
	locj = cell(2)
	lock = cell(3)

	nx = geo%nx 
	ny = geo%ny
	nz = geo%nz

! clear velocity arrays

      qxnode = 0.d0
	  qynode = 0.d0
	  qznode = 0.d0
        
! calculates velocity at nodes of the cell

    do i=1,2
	  do j=1,2
         do k=1,2

            sum=0.d0
            do m=0,1
			  do n=0,1
			      icell = loci+i-1
				  jcell = locj+j-1-m
				  kcell = lock+k-1-n
			      if (icell < 1 .or. icell > nx+1 .or. &
				      jcell < 1 .or. jcell > ny   .or. &
				      kcell < 1 .or. kcell > nz ) cycle
					  wx = weightx (geo,icell,jcell,kcell)
                      qxnode(i,j,k) = qxnode(i,j,k) + wx * value_array_ (this%qx,icell,jcell,kcell)
                      sum = sum + wx
			  end do
			end do
			if (sum > 0.d0) qxnode(i,j,k) = qxnode(i,j,k) / sum
            
            sum=0.d0
			do m=0,1
			  do n=0,1
			      icell = loci+i-1-m
				  jcell = locj+j-1
				  kcell = lock+k-1-n
			      if (icell < 1 .or. icell > nx   .or. &
			          jcell < 1 .or. jcell > ny+1 .or. &
				      kcell < 1 .or. kcell > nz ) cycle
					  wy = weighty (geo,icell,jcell,kcell)
                      qynode(i,j,k) = qynode(i,j,k) + wy * value_array_ (this%qy,icell,jcell,kcell)
                      sum = sum + wy
			  end do
			end do
			if (sum > 0.d0) qynode(i,j,k) = qynode(i,j,k) / sum

            sum=0.d0
            do m=0,1
			  do n=0,1
			      icell = loci+i-1-m
				  jcell = locj+j-1-n
				  kcell = lock+k-1
			      if (icell < 1 .or. icell > nx   .or. &
			          jcell < 1 .or. jcell > ny   .or. &
				      kcell < 1 .or. kcell > nz+1 ) cycle
                      wz = weightz (geo,icell,jcell,kcell)
                      qznode(i,j,k) = qznode(i,j,k) + wz * value_array_ (this%qz,icell,jcell,kcell)
					  sum = sum + wz
			  end do
			end do
			if (sum > 0.d0) qznode(i,j,k) = qznode(i,j,k) / sum
          
	        end do
	   end do
	end do


   end subroutine


   function weightx (geo,i,j,k) result (value)
	  use geometry_class
      use array_class
	  implicit none
	  type(geometry_cl),  intent(in) :: geo
	  integer,            intent(in) :: i,j,k
	  real*8                         :: dy,dz,value

      dy = 0.5d0*value_array_ (geo%dy,1,j,1)
	  dz = 0.5d0*value_array_ (geo%dz,1,1,k)

      value = 1.d0 / dsqrt(dy*dy+dz*dz)

   end function

   function weighty (geo,i,j,k) result (value)
      use array_class
	  use geometry_class
      implicit none
	  type(geometry_cl),  intent(in) :: geo
	  integer,            intent(in) :: i,j,k
	  real*8                         :: dx,dz,value

      dx = 0.5d0*value_array_ (geo%dx,i,1,1)
	  dz = 0.5d0*value_array_ (geo%dz,1,1,k)

      value = 1.d0 / dsqrt(dx*dx+dz*dz)

   end function

      function weightz (geo,i,j,k) result (value)
      use array_class
	  use geometry_class
      implicit none
	  type(geometry_cl),  intent(in) :: geo
	  integer,            intent(in) :: i,j,k
	  real*8                         :: dx,dy,value

      dx = 0.5d0*value_array_ (geo%dx,i,1,1)
	  dy = 0.5d0*value_array_ (geo%dy,1,j,1)

      value = 1.d0 / dsqrt(dy*dy+dx*dx)

   end function


!**********************************************************************************************
!      calculate velocity at faces
!**********************************************************************************************
      subroutine vel_faces_ (this,cell,qxface,qyface,qzface)
      use array_class
	  implicit none
	  type (advection_cl), intent(in)  :: this 
	  real*8,              intent(out) :: qxface(2),qyface(2),qzface(2)
	  integer,             intent(in)  :: cell(3)
      integer                          :: loci,locj,lock
      real*8						   :: qxijk,qyijk,qzijk

      loci = cell(1)
	  locj = cell(2)
	  lock = cell(3)
      
      qxijk = value_array_ (this%qx,loci,locj,lock)
      qyijk = value_array_ (this%qy,loci,locj,lock)
      qzijk = value_array_ (this%qz,loci,locj,lock)
      
      !if (qxijk>=0) then
      !  qxface(1) = qxijk
      !  qxface(2) = value_array_ (this%qx,loci+1,locj,lock)
      !else
      !  qxface(1) = value_array_ (this%qx,loci-1,locj,lock)
      !  qxface(2) = qxijk
      !end if
      
      !if (qyijk>=0) then
      !  qyface(1) = qyijk
      !  qyface(2) = value_array_ (this%qy,loci,locj+1,lock)
      !else
      !  qyface(1) = value_array_ (this%qy,loci,locj-1,lock)
      !  qyface(2) = qyijk
      !end if
      
      !if (qzijk>=0) then
      !  qzface(1) = qzijk
      !  qzface(2) = value_array_ (this%qz,loci,locj,lock+1)
      !else
      !  qzface(1) = value_array_ (this%qz,loci,locj,lock-1)
      !  qzface(2) = qzijk
      !end if
      
	  qxface(1) = value_array_ (this%qx,loci,locj,lock)
	  qxface(2) = value_array_ (this%qx,loci+1,locj,lock)
      
	  qyface(1) = value_array_ (this%qy,loci,locj,lock)
	  qyface(2) = value_array_ (this%qy,loci,locj+1,lock)
   
	  qzface(1) = value_array_ (this%qz,loci,locj,lock)
	  qzface(2) = value_array_ (this%qz,loci,locj,lock+1)


      end subroutine

!********************************************************************************************
!   SUBROUTINE TO READ FACE FLUXES FROM MODFLOW
!********************************************************************************************
! This program reads the output binary file from modflow cell-by-cell flow terms
! and calculates the darcy velocity field.
!
! - qx,qy,qz are darcy velocities of grid interfaces
!
! - watch-out: axes from modflow are not equal to axes in RW3D
!
! - how to go from (i,j,k) to modflow coordinates (jmod,imod,kmod):
!
!     jmod = i
!     imod = nrow-j+1
!     kmod = nlay-k+1
!
!*********************************************************************************************
  subroutine read_flux_from_mf2k_1 (fname,advection,geo)
	 use array_class
	 use geometry_class
	 use gslib,only: generate_unit,open_fname
	 use global_variables, only: fdbg
	 use loops_particles, only: nmove
     use velocity_times
	 implicit none
     character(len=*),   intent(in)    :: fname
	 type(geometry_cl),  intent(in)    :: geo
	 type(advection_cl), intent(inout) :: advection
     character(len=16)                 :: text
	 real, allocatable                 :: buff(:,:,:)
	 real*8                            :: dx,dy,dz
     integer*4                         :: kstp,kper,ncol,nrow,nlay
     integer                           :: iunit,ierror,ioerr,imod,jmod,kmod,i,j,k,iufdbg,alloc_err
	 logical                           :: exists,KeepReading
	 real*8                            :: Qsto
	 integer                           :: per,stp
	 
     per = advection%kper
     stp = advection%kstp

! inquire files and open them

     call open_fname (fname,iunit,'unknown','sequential','binary')
     if(nmove==1) call open_fname (fdbg,iufdbg,'unknown','append','formatted')
     
	 if(nmove==1) write(iufdbg,*)
	 write(*,*)
     write(*,*) 'reading flux from modflow',per,stp
     if(nmove==1) write(iufdbg,*) 'reading flux from modflow',per,stp


    KeepReading = .TRUE.

    do
    
        if (.not.KeepReading) return

        ! read identification line for constant head cells

        read(iunit) kstp,kper,text,ncol,nrow,nlay
    
        ! check if binary file is opened correctely and corresponds to the problem at hand

        if (ncol /= geo%nx .or. nrow /= geo%ny .or. nlay /= geo%nz) then
	       ! try to open with "unformatted" form
            close(iunit)
            call open_fname (fname,iunit,'unknown','sequential','unformatted')
            read(iunit) kstp,kper,text,ncol,nrow,nlay
		    if (ncol /= geo%nx .or. nrow /= geo%ny .or. nlay /= geo%nz) then
		        stop '*** Fatal error: modflow binary file with different mesh geometry'
		    end if
	    end if

        if(nmove==1) write(iufdbg,*) kstp,kper,text,ncol,nrow,nlay

        !if (kstp /= stp .or. kper /= per) stop 'velocity failure: STEP/PERIOD not correct'

        if (kstp == stp .and. kper == per) then
            KeepReading = .FALSE. 
        elseif (kper > per  .or. (kper == per .and. kstp > stp) ) then
            KeepReading = .TRUE.
        else
            stop 'velocity failure: STEP/PERIOD not correct'
        end if
    
       
        ! allocate arrays

        allocate (buff(ncol,nrow,nlay),stat=alloc_err)
	    if (alloc_err /=0) stop '>> Allocating error reading velocities from modflow' 
	    buff = 0.
	
	    if (.not.associated(advection%qx%values)) advection%qx = make_array_ (0.d0,ncol+1,nrow,  nlay  )
	    if (.not.associated(advection%qy%values)) advection%qy = make_array_ (0.d0,ncol,  nrow+1,nlay  )
	    if (.not.associated(advection%qz%values)) advection%qz = make_array_ (0.d0,ncol,  nrow,  nlay+1)

        ! initialize arrays

        advection%qx%values = 0.d0
        advection%qy%values = 0.d0
        advection%qz%values = 0.d0

        ! read storage if stress period is transient

        if( trim(adjustl(text)) == 'STORAGE') then
  
            read(iunit) buff
  
        !        do kmod=1,nlay
        !           do imod=1,nrow
        !             do jmod=1,ncol
        !               Qsto= buff(jmod,imod,kmod)
        !             end do
        !           end do
        !        end do

        ! read identification line for flow thru right face 
  
            read(iunit,iostat=ioerr,err=11,end=11) kstp,kper,text,ncol,nrow,nlay
            if(nmove==1) write(iufdbg,*) kstp,kper,text,ncol,nrow,nlay
  
        else

            if(nmove==1) write(iufdbg,*) ' skip storage'

        end if

        ! read flow from constant head cells

        if ( trim(adjustl(text)) == 'CONSTANT HEAD') then   
	
	        read(iunit) buff

            ! qx constant heads:

            do kmod=1,nlay
	        do imod=2,nrow-1
	             j  = nrow-imod+1
                 k  = nlay-kmod+1
	             dy = value_array_ (geo%dy,1,j,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qx%values(1,j,k)      = + buff(1,imod,kmod)/dy/dz
                 advection%qx%values(ncol+1,j,k) = - buff(ncol,imod,kmod)/dy/dz
	        end do
	        end do

	           imod=1
            do kmod=1,nlay
	             j  = nrow-imod+1
                 k  = nlay-kmod+1
	             dy = value_array_ (geo%dy,1,j,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qx%values(1,j,k)      = + 0.5d0 * buff(1,imod,kmod)/dy/dz
                 advection%qx%values(ncol+1,j,k) = - 0.5d0 * buff(ncol,imod,kmod)/dy/dz
	        end do

	           imod=nrow
            do kmod=1,nlay
	             j  = nrow-imod+1
                 k  = nlay-kmod+1
	             dy = value_array_ (geo%dy,1,j,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qx%values(1,j,k)      = + 0.5d0 * buff(1,imod,kmod)/dy/dz
                 advection%qx%values(ncol+1,j,k) = - 0.5d0 * buff(ncol,imod,kmod)/dy/dz
	        end do
	
            ! qy constant heads:

            do kmod=1,nlay
	        do jmod=2,ncol-1
                 i  = jmod
                 k  = nlay-kmod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qy%values(i,nrow+1,k) = - buff(jmod,1,kmod)/dx/dz
		         advection%qy%values(i,1,k)      = + buff(jmod,nrow,kmod)/dx/dz
	        end do
	        end do

	           jmod=1
            do kmod=1,nlay
                 i  = jmod
                 k  = nlay-kmod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qy%values(i,nrow+1,k) = - 0.5d0 * buff(jmod,1,kmod)/dx/dz
		         advection%qy%values(i,1,k)      = + 0.5d0 * buff(jmod,nrow,kmod)/dx/dz
	        end do
	
	           jmod=ncol
            do kmod=1,nlay
                 i  = jmod
                 k  = nlay-kmod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qy%values(i,nrow+1,k) = - 0.5d0 * buff(jmod,1,kmod)/dx/dz
		         advection%qy%values(i,1,k)      = + 0.5d0 * buff(jmod,nrow,kmod)/dx/dz
	        end do

            ! qz constant heads:
	
	        do imod=2,nrow-1
	        do jmod=2,ncol-1
                 i  = jmod
	             j  = nrow-imod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dy = value_array_ (geo%dy,1,j,1)
	             advection%qz%values(i,j,nlay+1) = - buff(jmod,imod,1)/dx/dy
		         advection%qz%values(i,j,1)      = + buff(jmod,imod,nlay)/dx/dy
	        end do
	        end do

        else
  
            if(nmove==1) write(iufdbg,*) ' skip flow thru constant heads'
            backspace(iunit)
    
        end if

        ! read identification line for flow thru right face 

        read(iunit,iostat=ioerr,err=11,end=11) kstp,kper,text,ncol,nrow,nlay

        ! read flow thru right faces

        if ( trim(adjustl(text)) == 'FLOW RIGHT FACE') then   

        if(nmove==1) write(iufdbg,*) kstp,kper,text,ncol,nrow,nlay

            buff = 0.
	        read(iunit) buff

            ! assign buffer to qx

            do kmod=1,nlay
	        do imod=1,nrow
	        do jmod=1,ncol-1
            i=jmod
	        j=nrow-imod+1
            k=nlay-kmod+1
	        dy = value_array_ (geo%dy,1,j,1)
	        dz = value_array_ (geo%dz,1,1,k)
	        advection%qx%values(i+1,j,k) = advection%qx%values(i+1,j,k) + buff(jmod,imod,kmod)/dy/dz
	        end do
	        end do
	        end do

        else

            if(nmove==1) write(iufdbg,*) ' skip flow thru right faces'
            backspace(iunit)
  
        end if

        11 if (ioerr.ne.0 .and. nmove==1) write(iufdbg,*) ' skip flow thru right faces'

        ! read identification line for flow thru front face

        read(iunit,iostat=ioerr,err=12,end=12) kstp,kper,text,ncol,nrow,nlay

        ! read flow thru front faces

        if ( trim(adjustl(text)) == 'FLOW FRONT FACE') then   

	        if(nmove==1) write(iufdbg,*) kstp,kper,text,ncol,nrow,nlay

            buff = 0.

	        read(iunit) buff

            ! assign buffer to qy

            do kmod=1,nlay
	        do imod=1,nrow  
	        do jmod=1,ncol
            i=jmod
	        j=nrow-imod
            k=nlay-kmod+1
	        dx = value_array_ (geo%dx,i,1,1)
	        dz = value_array_ (geo%dz,1,1,k)
	        advection%qy%values(i,j+1,k) = advection%qy%values(i,j+1,k) - buff(jmod,imod,kmod)/dx/dz
	        end do
	        end do
	        end do

        else

            if(nmove==1) write(iufdbg,*) ' skip reading flow thru front faces in modflow'
            backspace(iunit)

        end if

        12 if (ioerr.ne.0) write(iufdbg,*) ' skip flow thru front faces'

        ! read identification line for flow thru lower face

	    read(iunit,iostat=ioerr,err=14,end=14) kstp,kper,text,ncol,nrow,nlay

        ! read flow thru lower faces

        if ( trim(adjustl(text)) == 'FLOW LOWER FACE') then   

	    if(nmove==1) write(iufdbg,*) kstp,kper,text,ncol,nrow,nlay

            buff = 0.
	        read(iunit,iostat=ioerr,err=14,end=14) buff

            ! assign buffer to qz

            do kmod=1,nlay-1
	        do imod=1,nrow
	        do jmod=1,ncol
                i=jmod
	            j=nrow-imod+1
                k=nlay-kmod
	            dx = value_array_ (geo%dx,i,1,1)
	            dy = value_array_ (geo%dy,1,j,1)
	            advection%qz%values(i,j,k+1)=advection%qz%values(i,j,k+1)-buff(jmod,imod,kmod)/dx/dy
	        end do
	        end do
	        end do

        else

            if(nmove==1) write(iufdbg,*) ' skip reading flow thru lower faces in modflow'
            backspace(iunit)

        endif

        14 continue
    
	    if (ioerr.ne.0 .and. nmove==1) write(iufdbg,*) ' skip flow thru lower faces'


        ! read wells: border wells used as prescribed flux condition

	    read(iunit,iostat=ioerr,err=15,end=15) kstp,kper,text,ncol,nrow,nlay

        if ( trim(adjustl(text)) == 'WELLS') then   

	    if(nmove==1) write(iufdbg,*) kstp,kper,text,ncol,nrow,nlay

            buff = 0.
	        read(iunit,iostat=ioerr,err=15,end=15) buff

            ! qx border wells:

            do kmod=1,nlay
	        do imod=2,nrow-1
	             j  = nrow-imod+1
                 k  = nlay-kmod+1
	             dy = value_array_ (geo%dy,1,j,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qx%values(1,j,k)      = advection%qx%values(1,j,k)      + buff(1,imod,kmod)/dy/dz
                 advection%qx%values(ncol+1,j,k) = advection%qx%values(ncol+1,j,k) - buff(ncol,imod,kmod)/dy/dz
	        end do
	        end do

	           imod=1	
            do kmod=1,nlay
	             j  = nrow-imod+1
                 k  = nlay-kmod+1
	             dy = value_array_ (geo%dy,1,j,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qx%values(1,j,k)      = advection%qx%values(1,j,k)      + 0.5d0 * buff(1,imod,kmod)/dy/dz
                 advection%qx%values(ncol+1,j,k) = advection%qx%values(ncol+1,j,k) - 0.5d0 * buff(ncol,imod,kmod)/dy/dz
	        end do

	           imod=nrow
            do kmod=1,nlay
	             j  = nrow-imod+1
                 k  = nlay-kmod+1
	             dy = value_array_ (geo%dy,1,j,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qx%values(1,j,k)      = advection%qx%values(1,j,k)      + 0.5d0 * buff(1,imod,kmod)/dy/dz
                 advection%qx%values(ncol+1,j,k) = advection%qx%values(ncol+1,j,k) - 0.5d0 * buff(ncol,imod,kmod)/dy/dz
	        end do

            ! qy border wells:

            do kmod=1,nlay
	        do jmod=2,ncol-1
                 i  = jmod
                 k  = nlay-kmod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qy%values(i,nrow+1,k) = advection%qy%values(i,nrow+1,k) - buff(jmod,1,kmod)/dx/dz
		         advection%qy%values(i,1,k)      = advection%qy%values(i,1,k)      + buff(jmod,nrow,kmod)/dx/dz
	        end do
	        end do

	           jmod=1
            do kmod=1,nlay
                 i  = jmod
                 k  = nlay-kmod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qy%values(i,nrow+1,k) = advection%qy%values(i,nrow+1,k) - 0.5d0 * buff(jmod,1,kmod)/dx/dz
		         advection%qy%values(i,1,k)      = advection%qy%values(i,1,k)      + 0.5d0 * buff(jmod,nrow,kmod)/dx/dz
	        end do

	           jmod=ncol  
            do kmod=1,nlay

                 i  = jmod
                 k  = nlay-kmod+1
	             dx = value_array_ (geo%dx,i,1,1)
	             dz = value_array_ (geo%dz,1,1,k)
	             advection%qy%values(i,nrow+1,k) = advection%qy%values(i,nrow+1,k) - 0.5d0 * buff(jmod,1,kmod)/dx/dz
		         advection%qy%values(i,1,k)      = advection%qy%values(i,1,k)      + 0.5d0 * buff(jmod,nrow,kmod)/dx/dz
	        end do

            ! qz constant heads:

            ! needs improvements

        else

            if(nmove==1) write(iufdbg,*) ' skip reading wells in modflow'
            backspace(iunit)

        endif

        15 continue
    
	    if (ioerr.ne.0 .and. nmove==1) write(iufdbg,*) ' skip reading wells'

        deallocate (buff)

	    write(*,*)

        !close(iunit) !never close the unit for transient conditions 
	    if(nmove==1) close(iufdbg)
    
    end do
 
  end subroutine


!********************************************************************************************
!   SUBROUTINE TO READ FACE FLUXES FROM MODFLOW
!********************************************************************************************
! This program reads the output binary file from modflow cell-by-cell flow terms
! and calculates the darcy velocity field.
!
! - qx,qy,qz are darcy velocities of grid interfaces
!
! - watch-out: axes from modflow are not equal to axes in RW3D
!
! - to go from (i,j,k) to modflow coordinates (jmod,imod,kmod):
!
!     jmod = i
!     imod = nrow-j+1
!     kmod = nlay-k+1
!
!*********************************************************************************************
  subroutine read_flux_from_mf2k_2 (mf2k,advection,geo)
	 use array_class
	 use geometry_class
	 use gslib,only: generate_unit,open_fname, upper_case
	 use global_variables, only: fdbg, files_nam
	 use loops_particles, only: nmove
	 use mf2k_class
     use code_options, only: idebug
	 implicit none
	 type(mf2k_cl),      intent(in)    :: mf2k
	 type(geometry_cl),  intent(in)    :: geo
	 type(advection_cl), intent(inout) :: advection
     !character(len=16)                 :: text,ctmp
	 real*8                            :: dx,dy,dz
     !integer*4                         :: kstp,kper,ncol,nrow,nlay
     integer*4                         :: ncol,nrow,nlay
     
     integer                           :: iunit,ierror,ioerr,imod,jmod,kmod,i,j,k,iufdbg,alloc_err
	 logical                           :: exists
	 real*8                            :: Qsto
	 real*8                            :: dzero
	 !integer                           :: nc,nr,nl
	 integer                           :: inbud,iprec
     real*8                            :: deltd,pertimd,totimd,vald(20)
     real                              :: delt,pertim,totim,val(20)
     character(len=16)                 :: text1,text2,ctmp
     integer                           :: icode,nodes
     integer                           :: nlist,n,icell
     integer                           :: itype,nval
     integer, save, pointer            :: ibuff(:,:,:)
     real,    save, pointer            :: buff(:,:,:)
     real*8,  save, pointer            :: buffd(:,:,:)
     logical                           :: KeepReading,assignflow
     
     integer*4,         save           :: kstp,kper,nc,nr,nl
     character(len=16), save           :: text
     logical,           save           :: ReadHeader = .TRUE.
     
     real*8                            :: qxface(2),qyface(2),qzface(2)
	 integer                           :: cell(3),icol,jrow,klay,iunitVEL,nx,ny,nz
     real*8                            :: xc(3),qxval,qyval,qzval,qmodval

!....initialize

     dzero = 0.d0

     kper = advection%kper
     kstp = advection%kstp
     
     ncol = mf2k%nx
     nrow = mf2k%ny
     nlay = mf2k%nz
     
     if (.not.associated(buff))  allocate (buff(ncol,nrow,nlay))
     if (.not.associated(buffd)) allocate (buffd(ncol,nrow,nlay))
     if (.not.associated(ibuff)) allocate (ibuff(ncol,nrow,nlay))
     
     inbud = mf2k%unitCBC
     iprec = mf2k%preCBC

     assignflow = .FALSE.
     
     !ReadHeader = .TRUE.
!....read one time step


    loop: do
      
      if (ReadHeader == .TRUE.) then
        read(inbud,end=900,err=1000) kstp,kper,text,nc,nr,nl
      end if
      
      ! check if right stress period and time step is correct
      ! if still not there, exit and continue with the same velocity field
      if (kper > advection%kper  .or. (kper == advection%kper .and. kstp > advection%kstp) ) then
          !backspace(inbud) !pb with backspace for a binary file
          ReadHeader = .FALSE.
          exit
      end if
      ReadHeader = .TRUE.
      
      if(text.eq.'   CONSTANT HEAD') then
         write(*,12) 'Reading Budget for Period:', KPER, ' Timestep:', KSTP 
12       format(a30,i10,a12,i10) !,' rows',i10,' columns')
      end if
      
      
      if (kper == advection%kper  .and. kstp == advection%kstp ) assignflow =.TRUE.

      ! start reading one internal flow term
      
      itype=0
      if(nl.lt.0) then
         if(iprec.eq.1) then
           read(inbud) itype,delt,pertim,totim
         else
           read(inbud) itype,deltd,pertimd,totimd
         end if
         nval=1
         if(itype.eq.5) then
            read(inbud) nval
            if(nval.gt.1) then
               do n=2,nval
                 read(inbud) ctmp
               end do
            end if
         end if
         if(itype.eq. 2 .or. itype.eq.5) read(inbud) nlist
      end if

!------------------------------------------------------------------------------------------------------------------------------------------------------
!   read the budget term data under the following conditions:
!------------------------------------------------------------------------------------------------------------------------------------------------------
! if itype = 0 or 1, read a 3d array of values.
! if itype = 2 or 5, read a list of cells and their associated values.
! if itype = 3, read a 2d layer indicator array followed by a 2d array of ! values to be assigned to the layer indicated by the layer indicator array.
! if itype = 4 read a 2d array of values associated with layer 1.
!------------------------------------------------------------------------------------------------------------------------------------------------------
      if(itype.eq.0 .or. itype.eq.1) then
!  full 3-d array
         if(iprec.eq.1) then
           read(inbud) buff
           !buffd=buff
         else
           read(inbud) buffd
         end if
      else if(itype.eq.3) then
!  1-layer array with layer indicator array
         buffd=dzero
         read(inbud) ((ibuff(j,i,1),j=1,ncol),i=1,nrow)
         if(iprec.eq.1) then
           read(inbud) ((buff(j,i,1),j=1,ncol),i=1,nrow)
           do 265 i=1,nrow
           do 265 j=1,ncol
           buffd(j,i,1)=buff(j,i,1)
265        continue
         else
           read(inbud) ((buffd(j,i,1),j=1,ncol),i=1,nrow)
         end if
         do 270 i=1,nrow
         do 270 j=1,ncol
         if(ibuff(j,i,1).ne.1) then
            buffd(j,i,ibuff(j,i,1))=buffd(j,i,1)
            buffd(j,i,1)=dzero
         end if
270      continue
      else if(itype.eq.4) then
!  1-layer array that defines layer 1
         if(iprec.eq.1) then
           read(inbud) ((buff(j,i,1),j=1,ncol),i=1,nrow)
           do 275 i=1,nrow
           do 275 j=1,ncol
           buffd(j,i,1)=buff(j,i,1)
275        continue
         else
           read(inbud) ((buffd(j,i,1),j=1,ncol),i=1,nrow)
         end if
         if(nlay.gt.1) then
            do 280 k=2,nlay
            do 280 i=1,nrow
            do 280 j=1,ncol
            buffd(j,i,k)=dzero
280         continue
         end if
      else if(nlist.gt.0) then
!  list -- read only if the values need to be skipped.
         do 300 n=1,nlist
         if(iprec.eq.1) then
           read(inbud) icell,(val(i),i=1,nval)
!            k= (icell-1)/nrc + 1
!            i= ( (icell - (k-1)*nrc)-1 )/ncol +1
!            j= icell - (k-1)*nrc - (i-1)*ncol
         else
           read(inbud) icell,(vald(i),i=1,nval)
!            k= (icell-1)/nrc + 1
!            i= ( (icell - (k-1)*nrc)-1 )/ncol +1
!            j= icell - (k-1)*nrc - (i-1)*ncol
         end if
300      continue
      end if
!
!-----------------------------------------------------------------------
!     associate internal flow terms to darcy fluxes
!-----------------------------------------------------------------------      
      if (assignflow) then
      
      text = upper_case (text)
      
      if(text.eq.'   CONSTANT HEAD') then
      
         ! do nothing
      
      elseif(text.eq.'FLOW RIGHT FACE ') then
            
            ! assign buffer to qx

            advection%qx%values = 0.d0

            do kmod=1,nlay
	         do imod=1,nrow
	           do jmod=1,ncol-1
                   i=jmod
	               j=nrow-imod+1
                   k=nlay-kmod+1
	               dy = value_array_ (geo%dy,1,j,1)
	               dz = value_array_ (geo%dz,1,1,k)
	               if (iprec == 1) advection%qx%values(i+1,j,k) = buff(jmod,imod,kmod)/dy/dz
	               if (iprec /= 1) advection%qx%values(i+1,j,k) = buffd(jmod,imod,kmod)/dy/dz
	           end do
	         end do
	        end do     
      
      elseif(text.eq.'FLOW FRONT FACE ') then

            ! assign buffer to qy
     
            advection%qy%values = 0.d0
                        
            do kmod=1,nlay
	         do imod=1,nrow  
	          do jmod=1,ncol
                   i=jmod
	               j=nrow-imod
                   k=nlay-kmod+1
	               dx = value_array_ (geo%dx,i,1,1)
	               dz = value_array_ (geo%dz,1,1,k)
	               if (iprec == 1) advection%qy%values(i,j+1,k) = - buff(jmod,imod,kmod)/dx/dz
	               if (iprec /= 1) advection%qy%values(i,j+1,k) = - buffd(jmod,imod,kmod)/dx/dz
	          end do
	         end do
	        end do

      
      elseif(text.eq.'FLOW LOWER FACE ') then

            ! assign buffer to qz

            advection%qz%values = 0.d0

            do kmod=1,nlay-1
	         do imod=1,nrow
	          do jmod=1,ncol
                   i=jmod
	               j=nrow-imod+1
                   k=nlay-kmod
	               dx = value_array_ (geo%dx,i,1,1)
	               dy = value_array_ (geo%dy,1,j,1)
	               if (iprec == 1) advection%qz%values(i,j,k+1) = - buff(jmod,imod,kmod)/dx/dy
	               if (iprec /= 1) advection%qz%values(i,j,k+1) = - buffd(jmod,imod,kmod)/dx/dy
	          end do
	         end do
	        end do

      end if

      end if

      
    if (advection%action.and. idebug >=2) then
	    iunitVEL = generate_unit(100)
	    open(iunitVEL,file=trim(files_nam(14)))                 
        write(iunitVEL,*)   'ZONE I=',ncol,' J=',nrow,' K=',nlay,' F="POINT" '
	    do klay=1,nlay
		    do jrow=1,nrow
		    do icol=1,ncol
		    cell(1)=icol
		    cell(2)=jrow
		    cell(3)=klay
            call vel_faces_ (advection,cell,qxface,qyface,qzface) 
            xc = get_cell_centroid_ (geo,cell)
		    qxval = 0.5d0*(qxface(1)+qxface(2))
		    qyval = 0.5d0*(qyface(1)+qyface(2))
		    qzval = 0.5d0*(qzface(1)+qzface(2))
		    qmodval = dsqrt(qxval*qxval+qyval*qyval+qzval*qzval)
            if (qmodval > 0.d0) then
			    qmodval = log(qmodval)
		    else 
			    qmodval = 0.
		    end if
		    write(iunitVEL,'(3(g15.6,x),4(g20.12,x))') xc(1),xc(2),xc(3),qxval,qyval,qzval,qmodval
		    end do
		    end do
	    end do
	    close(iunitVEL)   
    end if   
        
    end do loop    

    
 return
 
 900 return  
1000 stop 'Error reading heading of fluxes'

     
     
  end subroutine

      



  end module

