
!********************************************************************************************
!   SYSTEM CLASS
!********************************************************************************************
  module sysinfo_class
  
  use plane_vect_class
  use well_vect_class

  implicit none

  private
  
  public :: sysinfo_cl
  public ::                                                  &
              assignment (=)                               , &
  			  update_system_particle_info_                 , &
			  print_system_info_                           , &
			  update_time_in_system_info_                  , &
			  update_move_in_system_info_                  , &
			  assign_exit_particle_                        , &
			  update_np_in_control_                        , &
			  count_exit_particle_                         , &
			  start_planes_in_system_                   

  type trans_prob_cl
      integer, pointer  :: count
	  real*8,  pointer  :: minim(:,:) => null() 
      real*8,  pointer  :: mean(:,:)  => null()
      real*8,  pointer  :: maxim(:,:) => null()
  end type

  type sysinfo_cl
      integer                       :: npbounce,npexceed,minmovesAll,maxmovesAll
	  integer                       :: minmovesOne,maxmovesOne         !minimum/maximum moves in one cell
      integer                       :: cellminmoves(3),cellmaxmoves(3) !cell minimum/maximum moves
	  real*8                        :: dtmin,dtmax                     !minim/maximum time step
	  integer                       :: celldtmin(3),celldtmax(3)       !cell minimum/maximum time step
	  integer                       :: npout                           !number particles exited
	  integer                       :: npstuck                         !number particles stuck
	  real*8, pointer               :: xpout(:),ypout(:),zpout(:)      !position particle exiting system
	  type(plane_vect_cl), pointer  :: planeDaI => null()              !plane vector with Damkohler information
	  type(well_vect_cl),  pointer  :: wellDaI  => null()              !well vector with Damkohler information
      integer, pointer              :: np_plane(:) => null() !number particles arriving at planes
	  integer, pointer              :: np_well(:)  => null() !number particles arriving at well
  end type

  interface assignment (=)
       module procedure is_equal_to_integer, is_equal_to_real; end interface

  contains

  subroutine assign_exit_particle_ (this,xp,yp,zp) 
      implicit none
	  real*8, intent(in)              :: xp,yp,zp
	  type(sysinfo_cl), intent(inout) :: this
	  real*8, allocatable             :: x(:),y(:),z(:)
	  integer                         :: n
		if ( this%npout == 0 .or. .not.associated(this%xpout)) then
		     allocate (this%xpout(1),this%ypout(1),this%zpout(1))
			 this%xpout(1)=xp
			 this%ypout(1)=yp
			 this%zpout(1)=zp
			 this%npout=1
		else
		     allocate (x(this%npout),y(this%npout),z(this%npout))
			 x = this%xpout
			 y = this%ypout
			 z = this%zpout
			 deallocate (this%xpout,this%ypout,this%zpout)
			 this%npout = this%npout + 1
			 n = this%npout
			 allocate (this%xpout(n),this%ypout(n),this%zpout(n))
			 this%xpout(1:n-1) = x(1:n-1)
			 this%ypout(1:n-1) = y(1:n-1)
			 this%zpout(1:n-1) = z(1:n-1)
			 this%xpout(n:n)   = xp
			 this%ypout(n:n)   = yp
			 this%zpout(n:n)   = zp
        end if
  end subroutine


  subroutine count_exit_particle_ (this) 
      implicit none
	  type(sysinfo_cl), intent(inout) :: this
			 this%npout = this%npout + 1
  end subroutine

  subroutine start_planes_in_system_ (this,nplane,nwell)
      implicit none
	  integer, intent(in)             :: nplane,nwell
	  type(sysinfo_cl), intent(inout) :: this
	  	  if (.not.associated(this%np_plane)) then
		       allocate(this%np_plane(nplane))
               this%np_plane = 0
		  end if
		  if (.not.associated(this%np_well )) then
		       allocate(this%np_well(nwell))
               this%np_well = 0
		  end if
  end subroutine


  subroutine is_equal_to_integer (this,value)
    use global_variables
	implicit none
	type(sysinfo_cl), intent(inout) :: this
	integer,          intent(in)    :: value      
	   if (value == 0.d0) then
	      this%npbounce    = 0
	      this%npexceed    = 0
	      this%minmovesOne = nmaxmove
	      this%maxmovesOne = 0  
	      this%minmovesAll = nmaxmove
	      this%maxmovesAll = 0 
	      this%dtmin       = 1.0d21
	      this%dtmax       = 0.d0
	      this%npout       = 0 
		  this%npstuck     = 0
		  if (associated(this%np_plane)) this%np_plane = 0
		  if (associated(this%np_well )) this%np_well  = 0
       else
	      this%npbounce    = value
	      this%npexceed    = value
	      this%minmovesOne = value
	      this%maxmovesOne = value  
	      this%minmovesAll = value
	      this%maxmovesAll = value
	      this%dtmin       = value
	      this%dtmax       = value  
	      this%npout       = value
		  this%npstuck     = value
		  if (associated(this%np_plane)) this%np_plane = value
		  if (associated(this%np_well )) this%np_well  = value
	   end if  
  end subroutine

  subroutine is_equal_to_real (this,value)
    use global_variables
	implicit none
	type(sysinfo_cl), intent(inout) :: this
	real,             intent(in)    :: value      
	   if (value == 0.) then
	       this%npbounce = 0
	       this%npexceed = 0
	       this%minmovesOne = nmaxmove
	       this%maxmovesOne = 0  
	       this%minmovesAll = nmaxmove
	       this%maxmovesAll = 0  
	       this%dtmin       = 1.0d21
	       this%dtmax       = 0.d0
	       this%npout       = 0
		   this%npstuck     = 0 
		   if (associated(this%np_plane)) this%np_plane = 0
		   if (associated(this%np_well )) this%np_well  = 0

       else
	       this%npbounce = value
	       this%npexceed = value
	       this%minmovesOne = value
	       this%maxmovesOne = value  
	       this%minmovesAll = value
	       this%maxmovesAll = value
	       this%dtmin = value
	       this%dtmax = value 
	       this%npout = value
		   this%npstuck = value 
		   if (associated(this%np_plane)) this%np_plane = 0
		   if (associated(this%np_well )) this%np_well  = 0
	   end if  
  end subroutine

  subroutine print_system_info_ (this,plane,well,fname)
	 use gslib, only: open_fname
	 use loops_particles
	 use plane_vect_class
	 use well_vect_class
     implicit none
     type(sysinfo_cl),    intent(in) :: this
     type(plane_vect_cl), intent(in) :: plane
	 type(well_vect_cl),  intent(in) :: well
	 character(len=*),    intent(in) :: fname
	 integer                         :: iunit,nwell,nplane,i,nzone,j
	 integer                         :: zero = 0

      call open_fname(fname,iunit)
 
!      write(iunit,'(70("_"),/)')
!	  write(iunit,*) 'Number of injections: ',kinj 
!	  write(iunit,'(70("_"))') 
	  
	  if (associated(this%np_plane).or.associated(this%np_well)) then
	    write(iunit,*)
	    write(iunit,*) 'Number of Particles Arriving at Control Surfaces:'
        write(iunit,*)
	    nplane = size(this%np_plane)
	    do i=1,nplane
		  if (.not.associated(this%np_plane)) then
		  write(iunit,4) 'PLANE',i,'..........: ',zero
		  else
	      write(iunit,4) 'PLANE',i,'..........: ',this%np_plane(i)
		  end if
	    end do           
  	    nwell = size(this%np_well)
	    do i=1,nwell
          if (.not.associated(this%np_well )) then
	      write(iunit,4) 'WELL ',i,'..........: ',zero
		  else
		  write(iunit,4) 'WELL ',i,'..........: ',this%np_well(i)
		  end if
	    end do
      end if
	  
	  write(iunit,1)
	  write(iunit,1) 'No. Bounces ................: ',this%npbounce
	  write(iunit,1) 'No. Part. Exceed Max. Moves.: ',this%npexceed
	  write(iunit,1) 'Total Min Moves ............: ',this%minmovesAll
	  write(iunit,1) 'Total Max Moves ............: ',this%maxmovesAll
!	  write(iunit,2) 'No. Min Moves One Cell......: ',this%minmovesOne,'    cell min moves (i,j,k)..: ',(this%cellminmoves(i),i=1,3)
!	  write(iunit,2) 'No. Max Moves One Cell......: ',this%maxmovesOne,'    cell max moves (i,j,k)..: ',(this%cellmaxmoves(i),i=1,3)
!	  write(iunit,3) 'Min Time Step...............: ',this%dtmin,      '    cell dt min (i,j,k).....: ',(this%celldtmin(i),i=1,3)
!	  write(iunit,3) 'Max Time Step...............: ',this%dtmax,      '    cell dt max (i,j,k).....: ',(this%celldtmax(i),i=1,3)
      write(iunit,*)
	  write(iunit,1) 'No. Particles Exit System...: ',this%npout
	  write(iunit,1) 'No. Particles Stuck ........: ',this%npstuck
	  write(iunit,*)

 1    format (5x,a30,i13)
 2    format (5x,a30,i13,a30,3(i6,x))
 3    format (5x,a30,g13.6,a30,3(i6,x))     
 4    format (5x,a5,i13,a12,i13)

	  if (this%npout > 0) then
	  if (associated(this%xpout)) then
	      write(iunit,*)
	      write(iunit,7) 'Xp_Exit','Yp_Exit','Zp_Exit'
	      write(iunit,*)
	         do i=1,this%npout
	            write(iunit,8) this%xpout(i),this%ypout(i),this%zpout(i) 
	         end do
		  write(iunit,*)
	  end if
	  end if

 7    format(3(a12,4x))
 8    format(3(g15.6,x))
     
	  close(iunit)

  end subroutine

  subroutine update_system_particle_info_ (this, particle, nmove )
     use particle_class
	 use global_variables
	 implicit none
	 type(sysinfo_cl),    intent(inout) :: this
	 type(particle_cl),   intent(in)    :: particle
	 integer,             intent(in)    :: nmove
	   if (nmove == nmaxmove) this%npexceed = this%npexceed + 1 ! count particle exceeds number of moves
	   if (this%minmovesAll > nmove) this%minmovesAll = nmove
	   if (this%maxmovesAll < nmove) this%maxmovesAll = nmove
       this%npbounce = this%npbounce + particle%control%npbounce ! count particle bounces 
	   if (particle%control%stuck) this%npstuck = this%npstuck + 1  
  end subroutine



    subroutine update_time_in_system_info_ ( this, particle, dt )
	    use particle_class
		implicit none
		type(sysinfo_cl),  intent(inout) :: this
		type(particle_cl), intent(in)    :: particle
		real*8,            intent(in)    :: dt
	    if ( this % dtmin > dt) then
		    this % dtmin = dt
			this % celldtmin = particle % cell % num; endif
	    if ( this % dtmax < dt) then
		    this % dtmax = dt
			this % celldtmax = particle % cell % num; endif
	end subroutine

    subroutine update_move_in_system_info_  ( this , particle, nmove )
	    use particle_class
		implicit none
		type(sysinfo_cl),  intent(inout) :: this 
		type(particle_cl), intent(in)    :: particle
		integer,           intent(in)    :: nmove
		integer,  save                   :: count
		if ( particle%control%switchcell ) then
		       if (this%minmovesOne > count .and. count /= 0 ) then
			           this%minmovesOne  = count
					   this%cellminmoves = particle%cell%num; endif
			   if (this%maxmovesOne < count .and. count /= 0) then
			           this%maxmovesOne = count
					   this%cellmaxmoves = particle%cell%num; endif
		    count = 0
			return
		else
		    count = count + 1
		end if
    end subroutine

    subroutine update_np_in_control_ (this,str,i)
	   use gslib, only: upper_case
	   implicit none
		type(sysinfo_cl),  intent(inout) :: this
		character(len=*),  intent(in)    :: str
		integer                          :: i  !plane or well number
		select case (trim(adjustl(upper_case(str))))
		   case('PLANE')
              if (.not.associated(this%np_plane)) then
			       allocate (this%np_plane(i))
                   this%np_plane(i) = 0
			  end if
              this%np_plane(i) = this%np_plane(i) + 1   
		   case('WELL')
              if (.not.associated(this%np_well)) then
			       allocate (this%np_well(i))
                   this%np_well(i) = 0
			  end if
              this%np_well(i) = this%np_well(i) + 1   
		end select
	end subroutine


  end module
