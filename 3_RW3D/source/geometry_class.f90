module geometry_class
    use array_class
	implicit none
	
	private
	public :: geometry_cl
	public ::                        & ! methods
              assign_geometry_     , &
			  dx_                  , &
			  dy_                  , &
			  dz_                  , &
			  locate_cell_         , &
			  delete_geometry_     , &
			  mesh_coordinates_    , &
			  get_xmesh_           , &
			  get_ymesh_           , &
			  get_zmesh_           , &
			  interpolate_array_   , &
			  read_geometry_       , &
			  read_boundary_       , &
			  get_cell_centroid_   , &
			  read_inactive_cells_ , &
			  length_cell_                

    type geometry_cl            
		 integer                 ::  nx,ny,nz          !number of cells
	     type(array_cl)          ::  dx,dy,dz          !contains the size of the mesh in x,y,z direction
		 type(array_cl)          ::  xmesh,ymesh,zmesh !interface position in x,y,z direction
		 integer,        pointer ::  ib(:,:) => null() !boundary information
		 integer,        pointer ::  InactCell(:,:,:) => null() !Defines inactive cells in the system
	end type

    interface get_cell_centroid_
	     module procedure get_cell_centroid1_, get_cell_centroid2_; end interface


	contains


    subroutine read_inactive_cells_ (this,fname,ivar,const,flag,nx,ny,nz) 
       use gslib, only: generate_unit
		 implicit none
         integer,           intent(in)    :: nx,ny,nz,ivar,flag
	     real*8,            intent(in)    :: const
         character(len=*),  intent(in)    :: fname
		 type(geometry_cl), intent(inout) :: this
		 integer                          :: ncell,iunit,i,j,k,nvar,icell,val,jcol
		 integer, allocatable             :: aline(:)
		 logical                          :: exists

         if (flag <= 0) then
		     return
         elseif (flag == 1) then !read in GSLIB format	      
		     inquire (file=fname,exist=exists)
	         if(.not.exists) then 
	           write(*,*) '>> file not found: ',trim(fname)
		       stop
	         end if
			 iunit = generate_unit(100)
             open(iunit,file=fname)
             !...allocate memory
             if(.not.associated(this%InactCell)) allocate (this%InactCell(nx,ny,nz))
			 this%InactCell = 1
             !...read inactive cells
             read(iunit,*,err=10,end=11)
	         read(iunit,*,err=10,end=11) nvar
	         do i=1,nvar
	            read(iunit,*,err=10,end=11)
	         end do
             allocate ( aline(nvar) )
             do k=1,nz
	           do j=1,ny
		         do i=1,nx
		         read(iunit,*,err=10,end=11) (aline(jcol),jcol=1,nvar)
				 this%InactCell(i,j,k) = aline(ivar) * const
		         end do
		       end do
	         end do
			 deallocate(aline)
			 close(iunit)
             return
         elseif (flag == 2) then !read in MODFLOW format
		     inquire (file=fname,exist=exists)
	         if(.not.exists) then 
	           write(*,*) '>> file not found: ',trim(fname)
		       stop
	         end if
			 iunit = generate_unit(100)
             open(iunit,file=fname)
             !...allocate memory
             if(.not.associated(this%InactCell)) allocate (this%InactCell(nx,ny,nz))
             this%InactCell = 1
             !...read inactive cells
             do k=nz,1,-1
	           do j=ny,1,-1
		         read(iunit,*,err=10,end=11) (this%InactCell(i,j,k),i=1,nx)
				 if (const /= 1.d0) then
				   do i=1,nx
				    this%InactCell(i,j,k) = this%InactCell(i,j,k) * const
		           end do
				 end if
			   end do
	         end do
			 close(iunit)
             return
         elseif (flag == 3) then !read in short format
		     inquire (file=fname,exist=exists)
	         if(.not.exists) then 
	           write(*,*) '>> file not found: ',trim(fname)
		       stop
	         end if
			 iunit = generate_unit(100)
             open(iunit,file=fname)
             !read file
			 read(iunit,*) ncell
             if(.not.associated(this%InactCell)) allocate (this%InactCell(nx,ny,nz))
			 this%InactCell = 1
             do icell=1,ncell
                  read(iunit,*) i,j,k
			      this%InactCell(i,j,k) = 0
             end do
			 close(iunit)
             return
         end if
             10   write(*,*) '>> error reading: ',trim(fname)
             11   write(*,*) '>> end-of-file reading: ',trim(fname)
             stop
	end subroutine


    subroutine read_geometry_ (this,unit) 
         use array_class
		 implicit none
		 integer,           intent(in)    :: unit
		 type(geometry_cl), intent(inout) :: this
		 integer                          :: nx,ny,nz,ivar,flag
		 character(len=200)               :: file
		 real*8                           :: const
           read(unit,*) nx,ny,nz;  this%nx=nx; this%ny=ny;  this%nz=nz
           read(unit,*) file,const,ivar,flag;  this%dx   = read_array_ (file,ivar,const,flag,nx,1,1)
           read(unit,*) file,const,ivar,flag;  this%dy   = read_array_ (file,ivar,const,flag,1,ny,1)
           read(unit,*) file,const,ivar,flag;  this%dz   = read_array_ (file,ivar,const,flag,1,1,nz)
    end subroutine

    subroutine read_boundary_ (this,unit) 
		 implicit none
		 integer,           intent(in)    :: unit
		 type(geometry_cl), intent(inout) :: this
		   if (.not.associated(this%ib)) allocate (this%ib(3,2))
	       read(unit,*) this%ib(1,1),this%ib(1,2),this%ib(2,1),this%ib(2,2),this%ib(3,1),this%ib(3,2)
    end subroutine

    subroutine assign_geometry_ (name,n1,n2,n3,d1,d2,d3,i) !constructor
	     use array_class
		 implicit none
		 type(geometry_cl), intent(inout)  :: name
		 integer, intent(in)               :: n1,n2,n3
	     real*8,  intent(in)               :: d1,d2,d3
		 integer, intent(in)               :: i(3,2)
		 character(len=1)                  :: fname = ' '
            name % nx = n1
			name % ny = n2
			name % nz = n3
            name % dx = read_array_ (fname,1,d1,0,1,1,1)
            name % dy = read_array_ (fname,1,d1,0,1,1,1)
			name % dz = read_array_ (fname,1,d1,0,1,1,1)
			name % ib = i
    end subroutine

    subroutine delete_geometry_ (name)  !destructor
         use array_class
		 implicit none
		 type(geometry_cl), intent(inout) :: name
		    call delete_array_ (name%dx)
			call delete_array_ (name%dy)
			call delete_array_ (name%dz)
			call delete_array_ (name%xmesh)
			call delete_array_ (name%ymesh)
			call delete_array_ (name%zmesh)
	end subroutine

    function length_cell_ (name,i,j,k) result (value) !accessor for dx
	     implicit none
		 type(geometry_cl), intent(in) :: name
		 integer,           intent(in) :: i,j,k
		 real*8                        :: value(3)
		 integer                       :: n
		    n = size(name%dx%values)
			if (n >  1) value(1) = name%dx%values(i,j,k)
			if (n == 1) value(1) = name%dx%values(1,1,1)
			n = size(name%dy%values)
			if (n >  1) value(2) = name%dy%values(i,j,k)
			if (n == 1) value(2) = name%dy%values(1,1,1)
		    n = size(name%dz%values)
			if (n >  1) value(3) = name%dz%values(i,j,k)
			if (n == 1) value(3) = name%dz%values(1,1,1)
	end function

    function dx_ (name,i,j,k) result (value) !accessor for dx
	     implicit none
		 type(geometry_cl), intent(in) :: name
		 integer,           intent(in) :: i,j,k
		 real*8                        :: value
		 integer                       :: n
		    n = size(name%dx%values)
			if (n >  1) value = name%dx%values(i,j,k)
			if (n == 1) value = name%dx%values(1,1,1)
	end function
		   
    function dy_ (name,i,j,k) result (value) !accessor for dy
	     implicit none
		 type(geometry_cl), intent(in) :: name
		 integer,           intent(in) :: i,j,k
		 real*8                        :: value
		 integer                       :: n
		    n = size(name%dy%values)
			if (n >  1) value = name%dy%values(i,j,k)
			if (n == 1) value = name%dy%values(1,1,1)
	end function
    
	function dz_ (name,i,j,k) result (value) !accessor for dz
	     implicit none
		 type(geometry_cl), intent(in) :: name
		 integer,           intent(in) :: i,j,k
		 real*8                        :: value
		 integer                       :: n
		    n = size(name%dz%values)
			if (n >  1) value = name%dz%values(i,j,k)
			if (n == 1) value = name%dz%values(1,1,1)
	end function

!***********************************************************************************************
!     locate cell 
!***********************************************************************************************
      subroutine locate_cell_ (this,xp,cell) 
	  use array_class
	  use gslib, only:locate
	  implicit none
	  type(geometry_cl), intent(in)    :: this
	  real*8,            intent(in)    :: xp(3)
	  integer,           intent(out)   :: cell(3)
	  real*8                           :: dx,dy,dz
	  real*8, allocatable              :: a(:)
      integer                          :: nx,ny,nz

      nx = length_array_ (this%dx,1)
	  ny = length_array_ (this%dy,2)
	  nz = length_array_ (this%dz,3)

      if ( nx == 1 ) then
             dx = value_array_ (this%dx)
             cell(1) = dint(xp(1) / dx) + 1
      else
	         allocate (a(nx)); a = this%dx%values(:,1,1)
	         call locate (a,nx,1,nx,xp(1),cell(1))
      end if

      if ( ny == 1 ) then
             dy = value_array_ (this%dy)
	         cell(2) = dint(xp(2) / dy) + 1
      else
	         allocate (a(ny)); a = this%dy%values(1,:,1)
	         call locate (a,ny,1,ny,xp(2),cell(2))
      end if

      if ( nz == 1 ) then
             dz = value_array_ (this%dz)
	         cell(3) = dint(xp(3) / dz) + 1
      else
	         allocate (a(nz)); a = this%dz%values(:,1,1)
	         call locate (a,nz,1,nz,xp(3),cell(3))
      end if


	  end subroutine
!*********************************************************************************************

      subroutine mesh_coordinates_ (name)
	     use array_class
		 implicit none
	     type(geometry_cl), intent(inout) :: name
		 integer                          :: n,i

		   n = name%nx
		   if (n>=1) then
		       name%xmesh = make_array_ (0.d0,n+1,1,1)
		       do i=1,n
		          name%xmesh%values(i+1,1,1) = name%xmesh%values(i,1,1) + value_array_ (name%dx,i,1,1)
		       end do
		   end if

		   n = name%ny
		   if (n>=1) then
		       name%ymesh = make_array_ (0.d0,1,n+1,1)
		       !do i=2,n
               do i=1,n
		          name%ymesh%values(1,i+1,1) = name%ymesh%values(1,i,1) + value_array_ (name%dy,1,i,1)
		       end do
		   end if

		   n = name%nz
		   if (n>=1) then
		       name%zmesh = make_array_ (0.d0,1,1,n+1)
		       !do i=2,n
               do i=1,n
		          name%zmesh%values(1,1,i+1) = name%zmesh%values(1,1,i) + value_array_ (name%dz,1,1,i)
		       end do
		   end if

	  end subroutine

      function get_xmesh_ (name,i) result (xmesh)
	     implicit none
		 type(geometry_cl), intent(in) :: name
	     integer,           intent(in) :: i
		 integer                       :: n
		 real*8                        :: xmesh
                n = length_array_ (name%dx)
                if (n == 1) then
				   xmesh = dfloat(i-1) * name%dx%values(1,1,1)
                elseif (n > 1) then
				   xmesh = name%xmesh%values(i,1,1)
				end if
	  end function


	  function get_cell_centroid1_ (name,cell) result (x)
	     implicit none
		 type(geometry_cl), intent(in)    :: name
		 integer,           intent(in)    :: cell(3)
		 real*8                           :: x(3)
		 integer                          :: n
                n = length_array_ (name%dx)
                if (n == 1) then
				   x(1) = (dfloat(cell(1)-1)+0.5d0) * name%dx%values(1,1,1)
                elseif (n > 1) then
				   x(1) = sum(name%dx%values(1:cell(1)-1,1,1)) + 0.5d0 * name%dx%values(cell(1),1,1)
				end if

				n = length_array_ (name%dy)
                if (n == 1) then
                   x(2) = (dfloat(cell(2)-1)+0.5d0) * name%dy%values(1,1,1)
                elseif (n > 1) then
				   x(2) = sum(name%dy%values(1,1:cell(2)-1,1)) + 0.5d0 * name%dy%values(1,cell(2),1)
				end if
                
				n = length_array_ (name%dz)
                if (n == 1) then
				   x(3) = (dfloat(cell(3)-1)+0.5d0) * name%dz%values(1,1,1)
                elseif (n > 1) then
				   x(3) = sum(name%dz%values(1,1,1:cell(3)-1)) + 0.5d0 * name%dz%values(1,1,cell(3))
				end if

      end function

	  function get_cell_centroid2_ (name,i,j,k) result (x)
	     implicit none
		 type(geometry_cl), intent(in)    :: name
		 integer,           intent(in)    :: i,j,k
		 real*8                           :: x(3)
		 integer                          :: n
                n = length_array_ (name%dx)
                if (n == 1) then
				   x(1) = (dfloat(i-1)+0.5d0) * name%dx%values(1,1,1)
                elseif (n > 1) then
				   x(1) = sum(name%dx%values(1:i-1,1,1)) + 0.5d0 * name%dx%values(i,1,1)
				end if

				n = length_array_ (name%dy)
                if (n == 1) then
                   x(2) = (dfloat(j-1)+0.5d0) * name%dy%values(1,1,1)
                elseif (n > 1) then
				   x(2) = sum(name%dy%values(1,1:j-1,1)) + 0.5d0 * name%dy%values(1,j,1)
				end if
                
				n = length_array_ (name%dz)
                if (n == 1) then
				   x(3) = (dfloat(k-1)+0.5d0) * name%dz%values(1,1,1)
                elseif (n > 1) then
				   x(3) = sum(name%dz%values(1,1,1:k-1)) + 0.5d0 * name%dz%values(1,1,k)
				end if

      end function



      function get_ymesh_ (name,j) result (ymesh)
	     implicit none
		 type(geometry_cl), intent(in) :: name
	     integer,           intent(in) :: j
		 integer                       :: n
		 real*8                        :: ymesh
                n = length_array_ (name%dy)
                if (n == 1) then
				   ymesh = dfloat(j-1) * name%dy%values(1,1,1)
                elseif (n > 1) then
				   ymesh = name%ymesh%values(1,j,1)
				end if
	  end function

      function get_zmesh_ (name,k) result (zmesh)
	     implicit none
		 type(geometry_cl), intent(in) :: name
	     integer,           intent(in) :: k
		 integer                       :: n
		 real*8                        :: zmesh
                n = length_array_ (name%dz)
                if (n == 1) then
				   zmesh = dfloat(k-1) * name%dz%values(1,1,1)
                elseif (n > 1) then
				   zmesh = name%zmesh%values(1,1,k)
				end if
	  end function


!*******************************************************************************************
!   TRILINEAR INTERPOLATION of the array properties at a given position
!*******************************************************************************************
!         
!            _________________
!          3| 131 | 231 | 331 |  
!           |_____|_____|_____|    
!          2| 121 | 222 | 321 |  
!           |_____|_____|_____|    
!          1| 111 | 211 | 311 |  
!           |_____|_____|_____|    
!              1     2     3 
!
!   ----------------------------------------------------------------------------------------
    function interpolate_array_ (this,array,cell,coord) result (value)
	   use array_class
	   use gslib, only: interpolation_trilinear
	   implicit none
	   type(geometry_cl), intent(in)    :: this           !geometry
	   type(array_cl),    intent(in)    :: array          !array of the property that is evaluated
	   real*8,            intent(in)    :: coord(3)       !local coordinates of the position in the cell
	   integer,           intent(in)    :: cell(3)
	   real*8                           :: value
	   real*8                           :: xx,yy,zz      !local coordinates of the position in the cell
	   real*8                           :: chi,nu,teta   !local coordinates of the position with respect adjoint cells
	   real*8                           :: dx(3),dy(3),dz(3)
	   real*8                           :: p(3,3,3),q(2,2,2)
	   integer                          :: col,row,lay,i,j,k,n
       logical                          :: inside

       n = length_array_ (array)
	   if (n == 1) then
          value = value_array_ (array,1,1,1)
          return
	   end if

       xx = coord(1)
       yy = coord(2)
       zz = coord(3)

       col = cell(1)
       row = cell(2)
	   lay = cell(3)

       do i=-1,1,1
	      dx(i+2) = value_array_geo (this%dx,this,col+i,1,1)
	      dy(i+2) = value_array_geo (this%dy,this,1,row+i,1)
	      dz(i+2) = value_array_geo (this%dz,this,1,1,lay+i)
	   end do

       do i=-1,1,1
	     do j=-1,1,1
		   do k=-1,1,1
           p(i+2,j+2,k+2) = value_array_geo (array,this,col+i,row+j,lay+k)
           end do
		 end do
	   end do
       
!........bottom layer
	      
		  if (xx >= 0.5d0 .and. yy < 0.5d0 .and. zz < 0.5d0) then
		      
			  chi  = ( xx*dx(2) - 0.5d0*dx(2) ) / ( 0.5d0*dx(2) + 0.5d0*dx(3) )
			  nu   = ( yy*dy(2) + 0.5d0*dy(1) ) / ( 0.5d0*dy(1) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) + 0.5d0*dz(1) ) / ( 0.5d0*dz(1) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p( 2+i-1, 1+j-1, 1+k-1 )
				  end do
				end do
			  end do

		  else if (xx >= 0.5d0 .and. yy >= 0.5d0 .and. zz < 0.5d0) then

			  chi  = ( xx*dx(2) - 0.5d0*dx(2) ) / ( 0.5d0*dx(2) + 0.5d0*dx(3) )
			  nu   = ( yy*dy(2) - 0.5d0*dy(2) ) / ( 0.5d0*dy(3) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) + 0.5d0*dz(1) ) / ( 0.5d0*dz(1) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p(2+i-1,2+j-1,1+k-1)
				  end do
				end do
			  end do

		  else if (xx < 0.5d0 .and. yy >= 0.5d0 .and. zz < 0.5d0) then

			  chi  = ( xx*dx(2) + 0.5d0*dx(1) ) / ( 0.5d0*dx(1) + 0.5d0*dx(2) )
			  nu   = ( yy*dy(2) - 0.5d0*dy(2) ) / ( 0.5d0*dy(3) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) + 0.5d0*dz(1) ) / ( 0.5d0*dz(1) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p(1+i-1,2+j-1,1+k-1)
				  end do
				end do
			  end do

		  else if (xx < 0.5d0 .and. yy < 0.5d0 .and. zz < 0.5d0) then
         
			  chi  = ( xx*dx(2) + 0.5d0*dx(1) ) / ( 0.5d0*dx(1) + 0.5d0*dx(2) )
			  nu   = ( yy*dy(2) + 0.5d0*dy(1) ) / ( 0.5d0*dy(1) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) + 0.5d0*dz(1) ) / ( 0.5d0*dz(1) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p(1+i-1,1+j-1,1+k-1)
				  end do
				end do
			  end do

!........upper layer

	      
		  else if (xx >= 0.5d0 .and. yy < 0.5d0 .and. zz >= 0.5d0) then
		      
			  chi  = ( xx*dx(2) - 0.5d0*dx(2) ) / ( 0.5d0*dx(2) + 0.5d0*dx(3) )
			  nu   = ( yy*dy(2) + 0.5d0*dy(1) ) / ( 0.5d0*dy(1) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) - 0.5d0*dz(2) ) / ( 0.5d0*dz(3) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p( 2+i-1, 1+j-1, 2+k-1 )
				  end do
				end do
			  end do

		  else if (xx >= 0.5d0 .and. yy >= 0.5d0 .and. zz >= 0.5d0) then

			  chi  = ( xx*dx(2) - 0.5d0*dx(2) ) / ( 0.5d0*dx(2) + 0.5d0*dx(3) )
			  nu   = ( yy*dy(2) - 0.5d0*dy(2) ) / ( 0.5d0*dy(3) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) - 0.5d0*dz(2) ) / ( 0.5d0*dz(3) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p(2+i-1,2+j-1,2+k-1)
				  end do
				end do
			  end do

		  else if (xx < 0.5d0 .and. yy >= 0.5d0 .and. zz >= 0.5d0) then

			  chi  = ( xx*dx(2) + 0.5d0*dx(1) ) / ( 0.5d0*dx(1) + 0.5d0*dx(2) )
			  nu   = ( yy*dy(2) - 0.5d0*dy(2) ) / ( 0.5d0*dy(3) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) - 0.5d0*dz(2) ) / ( 0.5d0*dz(3) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p(1+i-1,2+j-1,2+k-1)
				  end do
				end do
			  end do

		  else if (xx < 0.5d0 .and. yy < 0.5d0 .and. zz >= 0.5d0) then
         
			  chi  = ( xx*dx(2) + 0.5d0*dx(1) ) / ( 0.5d0*dx(1) + 0.5d0*dx(2) )
			  nu   = ( yy*dy(2) + 0.5d0*dy(1) ) / ( 0.5d0*dy(1) + 0.5d0*dy(2) )
			  teta = ( zz*dz(2) - 0.5d0*dz(2) ) / ( 0.5d0*dz(3) + 0.5d0*dz(2) )

              do i=1,2
			    do j=1,2
				  do k=1,2
			      q(i,j,k) = p(1+i-1,1+j-1,2+k-1)
				  end do
				end do
			  end do
		 
		  end if

    call check_boundary (q,chi,nu,teta)

    value = interpolation_trilinear (q,chi,nu,teta)

	end function

    subroutine check_boundary (q,xx,yy,zz) !subroutine only used for function update_property
        implicit none
        real*8, intent(in)    :: q(2,2,2)
		real*8, intent(inout) :: xx,yy,zz

		       if (sum(q(:,:,1)) == 0.d0) zz = 1.d0
		       if (sum(q(:,:,2)) == 0.d0) zz = 0.d0
		  
		       if (sum(q(1,:,:)) == 0.d0) xx = 1.d0
		       if (sum(q(2,:,:)) == 0.d0) xx = 0.d0
		  
		       if (sum(q(:,1,:)) == 0.d0) yy = 1.d0
		       if (sum(q(:,2,:)) == 0.d0) yy = 0.d0          

	end subroutine


    function value_array_geo (array,geo,loci,locj,lock) result (value)  ! function only used for function update_property
         use array_class
		 implicit none
		 type (geometry_cl), intent(in) :: geo
		 type (array_cl),    intent(in) :: array
		 integer,            intent(in) :: loci,locj,lock
         real*8                         :: value
		 integer                        :: n	       			   
			   if (loci > geo%nx .or. loci < 1) then
			        value = 0.d0
					return
			   else if (locj > geo%ny .or. locj < 1) then
			        value = 0.d0
					return
			   else if (lock > geo%nz .or. lock < 1) then
			        value = 0.d0
					return
			   else
			          value = value_array_ (array,loci,locj,lock)
               end if
	end function 

!  ***********************************************************************************************************************


end module


