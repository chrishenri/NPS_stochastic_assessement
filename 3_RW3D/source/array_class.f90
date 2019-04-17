!*************************************************************************************************
!   ARRAY CLASS
!*************************************************************************************************
    module array_class
    implicit none

	private
	public :: array_cl                ! class
	public ::                       & ! methods
              assign_array_       , & 
			  read_array_         , &
			  delete_array_       , &
			  print_array_        , &
			  list_array_         , &
			  length_array_       , &
			  value_array_        , &
			  make_array_         , &
			  equal_to_           , &
			  print_array_normal_ , &
			  average_array_      , &
              list_array_1D_     		  
			  
	type array_cl
		 real*8, pointer, dimension(:,:,:) :: values => null()
	end type array_cl


	interface assign_array_
	     module procedure assign_array,make_array ;	end interface
	interface read_array_
	     module procedure read_array_gslib ; end interface
    interface delete_array_
	     module procedure delete_array ; end interface
    interface print_array_
	     module procedure print_array_gslib ; end interface
    interface list_array_
	     module procedure list_array ; end interface
    interface length_array_
	     module procedure length_array ; end interface
    interface value_array_
	     module procedure value_array1, value_array2; end interface
    interface make_array_
	     module procedure make_array; end interface

	contains

	function assign_array (datos) result (name)  !array constructor
	     real*8, intent(in) :: datos(:,:,:)
         integer            :: size_nx,size_ny,size_nz
		 integer            :: ok
		 type(array_cl)     :: name
		      size_nx = size (datos,dim=1)
			  size_ny = size (datos,dim=2)
			  size_nz = size (datos,dim=3)
			  allocate (name%values(size_nx,size_ny,size_nz),stat=ok)
			  if (ok /= 0) stop '>> array not deallocated in delete_array'
			  name % values = datos
	end function

	function average_array_ (name) result (avg)  !array constructor
		 type(array_cl), intent(in) :: name
         integer                    :: n
         real*8                     :: avg
		      avg = sum(dabs(name%values),mask=dabs(name%values)>0.d0 )
		      n   = count(mask=dabs(name%values)>0.d0 )
              if (n>0) then 
                 avg = avg / dfloat(n)
              else
                 avg = 0.d0
              end if
	end function


    function make_array (const,size_nx,size_ny,size_nz) result (name) !optional constructor
         implicit none
		 integer, optional, intent(in) :: size_nx,size_ny,size_nz
		 real*8,            intent(in) :: const
		 type(array_cl)                :: name
		 integer                       :: ok
		   if (present(size_nx) .and. present(size_ny) .and. present(size_nz)) then
			  allocate (name%values(size_nx,size_ny,size_nz),stat=ok)
  			  if (ok /= 0) stop '>> array not allocated in make_array'
			  name % values = const
		   else
			  allocate (name%values(1,1,1),stat=ok)
			  if (ok /= 0) stop '>> array not allocated in make_array'
			  name % values = const
           end if
     end function
 
     function equal_to_ (a) result (this)
         implicit none
		 type(array_cl)                 :: this
		 type(array_cl), intent(in)     :: a
	     integer                        :: nx,ny,nz
         if (associated(this%values)) call delete_array (this)
           nx = length_array (a,1)
	       ny = length_array (a,2)
	       nz = length_array (a,3)
	       this = make_array (0.d0,nx,ny,nz)
		   this%values = a%values
	 end function

                              
    function read_array_gslib (file_name,ivar,const,flag,nx,ny,nz) result (name)
    use gslib, only: generate_unit
	implicit none
      integer,         intent(in) :: nx,ny,nz,ivar,flag
	  real*8,          intent(in) :: const
      character(len=*),intent(in) :: file_name
      type(array_cl)              :: name
      integer                     :: nvar,i,j,k,jcol,iunit
      real*8, allocatable         :: aline(:)
      logical                     :: exists
	  integer                     :: nline,iline
	  integer, save               :: npriorline 
	  real*8                      :: val
          if (flag == 0 ) then              ! parameter is constant and is not necessary 
             name = make_array (const)
			 return	     
          else if (flag == 1) then          ! read from GSLIB file and start from first line
             name = make_array (0.d0,nx,ny,nz)
	         inquire (file=file_name,exist=exists)
	         if(.not.exists) then 
	           write(*,*) '>> file not found: ',trim(file_name)
		       stop
	         end if
			 iunit = generate_unit(100)
             open(iunit,file=file_name)
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
		         name % values (i,j,k) = aline(ivar) * const
		         end do
		       end do
	         end do
             close(iunit)
             return
           else if (flag == 3) then        ! read from GSLIB but start from last line
			 name = make_array (0.d0,nx,ny,nz)
	         inquire (file=file_name,exist=exists)
	         if(.not.exists) then 
	           write(*,*) '>> file not found: ',trim(file_name)
		       stop
	         end if
			 iunit = generate_unit(100)
             open(iunit,file=file_name)
             read(iunit,*,err=10,end=11)
	         read(iunit,*,err=10,end=11) nvar
	         do i=1,nvar
	            read(iunit,*,err=10,end=11)
	         end do
             allocate ( aline(nvar) )
			 nline = npriorline
			 do iline=1,nline
			     read(iunit,*)
			 end do
             do k=1,nz
	           do j=1,ny
		         do i=1,nx
		         read(iunit,*,err=10,end=11) (aline(jcol),jcol=1,nvar)
				 nline = nline + 1
		         name % values (i,j,k) = aline(ivar) * const
		         end do
		       end do
	         end do
			 npriorline = nline
             close(iunit)
             return 
           else if (flag == 2) then        ! read from modflow type of file
             name = make_array (0.d0,nx,ny,nz)
	         inquire (file=file_name,exist=exists)
	         if(.not.exists) then 
	           write(*,*) '>> file not found: ',trim(file_name)
		       stop
	         end if
			 iunit = generate_unit(100)
             open(iunit,file=file_name)
             read(iunit,*,err=10,end=11)
             read(iunit,*,err=10,end=11)
             do k=nz,1,-1
	           do j=ny,1,-1
		         do i=1,nx
		         read(iunit,*,err=10,end=11) val 
		         name % values (i,j,k) = val * const
		         end do
		       end do
	         end do
             close(iunit)
             return
              
             
           end if
             10   write(*,*) '>> error reading: ',trim(file_name)
             11   write(*,*) '>> end-of-file reading: ',trim(file_name)
             stop 'NOT  a normal termination'
	  end function read_array_gslib

      subroutine delete_array (name)   ! destructor
	      type (array_cl), intent (inout) :: name
		  integer                         :: ok
		      if (associated(name % values)) then
		      deallocate (name % values, stat = ok); end if
			  !if (ok /= 0) stop '>> array not deallocated in delete_array'
	  end subroutine

      subroutine print_array_gslib (name,fname,nombre)
          use gslib, only: open_fname
		  implicit none
		  character(len=*), optional, intent(in)     :: fname,nombre
		  type (array_cl),            intent(in)     :: name
          logical                                    :: exists
		  integer                                    :: iunit,i,j,k,nx,ny,nz
             if (.not.associated(name%values)) return
             call open_fname (fname,iunit)			    
		     write(iunit,*) ' print variable in gslib format'
		     write(iunit,*) ' 1'
		         if (present(nombre))then
			        write(iunit,*) nombre
			     else
			        write(iunit,*) ' variable'
			     end if
             nx = size(name%values,dim=1)
			 ny = size(name%values,dim=2)
			 nz = size(name%values,dim=3)
			 do k=1,nz
	           do j=1,ny
		         do i=1,nx
		         write(iunit,*) name % values (i,j,k) 
		         end do
		       end do
	         end do          
			 close(iunit)
	  end subroutine

      subroutine print_array_normal_ (name,fname,nombre) !print in GSLIB without Access end-of-file
          use gslib, only: open_fname_normal
		  implicit none
		  character(len=*), optional, intent(in)     :: fname,nombre
		  type (array_cl),            intent(in)     :: name
          logical                                    :: exists
		  integer                                    :: iunit,i,j,k,nx,ny,nz
             if (.not.associated(name%values)) return
             call open_fname_normal (fname,iunit)			    
		     write(iunit,*) ' print variable in gslib format'
		     write(iunit,*) ' 1'
		         if (present(nombre))then
			        write(iunit,*) nombre
			     else
			        write(iunit,*) ' variable'
			     end if
             nx = size(name%values,dim=1)
			 ny = size(name%values,dim=2)
			 nz = size(name%values,dim=3)
			 do k=1,nz
	           do j=1,ny
		         do i=1,nx
		         write(iunit,*) name % values (i,j,k) 
		         end do
		       end do
	         end do          
			 close(iunit)
	  end subroutine


	  subroutine list_array (name,fname,nombre)
          use gslib, only: open_fname
		  implicit none
		  character(len=*), intent(in)           :: nombre
		  character(len=*), intent(in)           :: fname
		  type (array_cl),  intent(in)           :: name
		  integer                                :: iunit 
		  integer                                :: i,j,k,size_nx,size_ny,size_nz
		     size_nx = size(name%values,dim=1)
			 size_ny = size(name%values,dim=2)
			 size_nz = size(name%values,dim=3) 
             call open_fname(fname,iunit)
             
             if (size(name%values) == 1) then
			       write(iunit,1) trim(nombre),'  ',name%values(1,1,1)
	            1  format(a,a,g15.6)
             else
			      write(iunit,1) trim(nombre),'  '
                  write(iunit,*)
			      do k=size_nz,1,-1
			         write(iunit,2)
			         do j=size_ny,1,-1 
			            write(iunit,2) (name%values(i,j,k),i=1,size_nx)
			         end do			    
			      end do
			      write(iunit,2) 
		       2  format(<size_nx>(g15.6,x))
			 end if
			 close(iunit)
       end subroutine

	  subroutine list_array_1D_ (name,fname,nombre)
          use gslib, only: open_fname
		  implicit none
		  character(len=*), intent(in)           :: nombre
		  character(len=*), intent(in)           :: fname
		  type (array_cl),  intent(in)           :: name
		  integer                                :: iunit 
		  integer                                :: i,j,k,size_nx,size_ny,size_nz
		     size_nx = size(name%values,dim=1)
			 size_ny = size(name%values,dim=2)
			 size_nz = size(name%values,dim=3) 
             call open_fname(fname,iunit)
             
             if (size(name%values) == 1) then
			       write(iunit,1) trim(nombre),'  ',name%values(1,1,1)
	            1  format(a,a,g15.6)
             else
			      write(iunit,1) trim(nombre),'  '
                  write(iunit,*)
			      write(iunit,*) name%values
			      write(iunit,*) 
			 end if
			 close(iunit)
       end subroutine


      function length_array (name,index) result (n)
	     integer, optional, intent(in) :: index
		 type (array_cl), intent(in)   :: name
         integer                       :: n
		    if (present(index)) then
			n = size (name%values,dim=index)
			else
			n = size (name%values)
			end if
	  end function

      function value_array1 (name,loci,locj,lock) result (value)        ! accessor member
         integer, optional, intent(in) :: loci,locj,lock
		 type (array_cl), intent(in)   :: name
         real*8                        :: value
		 integer                       :: n
         if (present(loci) .and. present(locj) .and. present(lock)) then
			    n  = 0
				n  = length_array (name)
		        if (n > 1) then
			          value = name % values(loci,locj,lock)
			    else if ( n==0 ) then
				      value = 0.d0
				else
			          value = name % values(1,1,1)
			    end if 
		 else
		              value = name % values(1,1,1)
		 end if
	  end function 

      function value_array2 (name,cell) result (value)        ! accessor member
         integer, intent(in)           :: cell(3)
		 type (array_cl), intent(in)   :: name
         real*8                        :: value
		 integer                       :: n
		        n  = length_array (name)
		        if (n > 1) then
			          value = name % values(cell(1),cell(2),cell(3))
			    else
			          value = name % values(1,1,1)
			    end if 
	  end function 


	end module
