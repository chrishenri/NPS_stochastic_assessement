
!*****************************************************************************************************
!   PLANE CLASS
!*****************************************************************************************************   
    module block_class
	use breakthru_class
	implicit none

	private
	public :: block_cl     ! class
	public ::                               & ! methods
              assign_block_               , & 
              print_block_                , &
			  read_assign_block_          !, &
			  !print_block_breakthru_        

	type block_cl
		 real*8                          :: x1=0.d0,x2=0.d0,y1=0.d0,y2=0.d0,z1=0.d0,z2=0.d0      !block corner
         real*8, pointer                 :: mass(:)                                              !True = particle is removed when crossing 
         type(breakthru_cl), pointer     :: btc(:)                                               !Number of particles
	end type block_cl

    !interface print_block_
	   !  module procedure print_block_1_,print_block_2_; end interface

    contains


    subroutine assign_block_ (name,x1,x2,y1,y2,z1,z2,curve)
         use breakthru_class
		 implicit none
		 type(block_cl),                intent(inout) :: name
		 real*8,                        intent(in)    :: x1,x2,y1,y2,z1,z2
         type(breakthru_cl),  optional, intent(in)    :: curve
         name % x1 = x1 
         name % x2 = x2 
         name % y1 = y1 
         name % y2 = y2 
         name % z1 = z1 
         name % z2 = z2 
         if (.not.present(curve)) then
            return
		 else
			name%btc = curve
		 end if
	end subroutine  


    subroutine read_assign_block_ (name,iunit)
	     use gslib, only: ifcharacter
		 implicit none
		 type(block_cl),  intent(inout) :: name
		 integer,         intent(in)    :: iunit
		 logical                        :: exists
		 real*8                         :: x1,x2,y1,y2,z1,z2
		 integer                        :: flag
         
         inquire (unit=iunit,exist=exists)
		 if(.not.exists) then
			print *, 'file does not exist during read_assign_plane'
			stop
		 end if
	     read(iunit,*) x1,y1,z1,x2,y2,z2
		 call assign_block_ (name,x1,x2,y1,y2,z1,z2)
		      
    end subroutine
    
    
    subroutine print_block_ (name,fname,ispe)
         use gslib, only: open_fname, upper_case
		 implicit none
		 type(block_cl),   intent(inout)         :: name
         character(len=*), intent(in)            :: fname
		 integer,          intent(in)            :: ispe
         
         
    end subroutine
    
 !   subroutine print_block_1_ (name,fname)
 !        use gslib, only: open_fname, upper_case
	!	 implicit none
	!	 type(plane_cl), intent(inout)         :: name
	!	 character(len=*),optional, intent(in) :: fname
	!	 integer                               :: iunit          
	!	   call open_fname(fname,iunit)
	!	   name%tipo = upper_case(name%tipo)
	!	   if(name%tipo == 'XX') then
 !               write(iunit,1) 'ZONE T= "Yplane: ',name%yplane,' Type: ',name%tipo,' OUT: ',name%out, '"'
	!	   else if (name%tipo == 'YY') then
 !               write(iunit,1) 'ZONE T= "Xplane: ',name%xplane,' Type: ',name%tipo,' OUT: ',name%out, '"'
	!	   else if (name%tipo == 'ZZ') then
 !               write(iunit,1) 'ZONE T= "Zplane: ',name%zplane,' Type: ',name%tipo,' OUT: ',name%out, '"'
	!	   else if (name%tipo == 'GL') then
 !               write(iunit,2) 'ZONE T= "Ax+By+Cz+D=0 (A,B,C,D): ',name%A,name%B,name%C,name%D,' Type: ',name%tipo,' OUT: ',name%out,'"'
	!	   end if
	!	 2 format(a33,4(g15.6,x),a7,a2,x,a6,i2,a2)
 !        1 format((a17,g15.6,5x),a7,a2,x,a6,i2,a2)                                                       
	!	   close(iunit)
	!end subroutine
 !
 !   subroutine print_block_2_ (name,fname,ispe)
 !        use gslib, only: open_fname, upper_case
	!	 implicit none
	!	 type(plane_cl), intent(inout)         :: name
	!	 integer,        intent(in)            :: ispe
	!	 character(len=*),optional, intent(in) :: fname
	!	 integer                               :: iunit          
	!	   call open_fname(fname,iunit)
	!	   name%tipo = upper_case(name%tipo)
	!	   if(name%tipo == 'XX') then
 !               write(iunit,1) 'ZONE T= "Yplane: ',name%yplane,' Specie: ',ispe, '"'
	!	   else if (name%tipo == 'YY') then
 !               write(iunit,1) 'ZONE T= "Xplane: ',name%xplane,' Specie: ',ispe, '"'
	!	   else if (name%tipo == 'ZZ') then
 !               write(iunit,1) 'ZONE T= "Zplane: ',name%zplane,' Specie: ',ispe, '"'
	!	   else if (name%tipo == 'GL') then
 !               write(iunit,2) 'ZONE T= "General): ',name%A,name%B,name%C,name%D,' Specie: ',ispe, '"'
	!	   end if
	!	 2 format(a19,4(f7.1,x),a9,i2,a2)
 !        1 format((a17,f7.1,x),a9,i2,a2)                                                       
	!	   close(iunit)
	!end subroutine
 !
 !  subroutine print_block_breakthru_ (i,btc,fname)
 !   use gslib, only: open_fname
	!use breakthru_class
	!implicit none
	!type(breakthru_cl),         intent(in)    :: btc
 !   integer,                    intent(in)    :: i
	!character(len=*), optional, intent(in)    :: fname
	!integer                                   :: n,iunit
	!  if (btc%np<5) return
	!     call open_fname(fname,iunit)
	!	 write(iunit,2)'zone plane="',i,'"'
	!  2  format(a12,i5,a2) 
	!	 call print_breakthru_ (btc,fname)
 !  end subroutine
    

    end module
    
    
    
    
!******************************************************************************
!  VECTOR OF PLANES  
!******************************************************************************
 module block_vect_class
 use block_class

 implicit none

 public    !everything is public (inheritance of plane_class)
 public :: block_vect_cl                        !class
 public ::                                    & !methods
           alloc_block_vect_                !, &
           !print_obs_block
		   !print_block_dispersivity_vect_   , &
		   !print_block_breakthru_vect_      , &
		   !nblock_                          , & 
		   !print_block_moments_vect_      

 type block_vect_cl
      type(block_cl),pointer :: num(:) => null()
 end type

 !
 !interface print_block_moments_vect_
 !     module procedure print_block_moments_vect_,print_block_moments_vect_1_; end interface

 contains

 
 subroutine alloc_block_vect_ (name,n,nspecie)
    use block_class
    use breakthru_class
	implicit none
    type(block_vect_cl), intent(inout) :: name
    integer,             intent(in)    :: n
    integer,             intent(in)    :: nspecie
	integer                            :: i,ispe 
	      allocate(name%num(n))
		  do i=1,n
		    allocate (name%num(i)%mass(nspecie))
			do ispe=1,nspecie
			   name%num(i)%mass(ispe) = 0.d0
               call initialize_btc_ (name%num(i)%btc(ispe))
			end do
		  end do	  
 end subroutine
 
 
 !subroutine print_plane_breakthru_vect_ (name,fname)
 !   use gslib, only: open_fname
	!use breakthru_class
	!use plane_class
	!use global_variables, only: nspecie
	!implicit none
	!type(plane_vect_cl),        intent(inout) :: name
	!character(len=*), optional, intent(in)    :: fname
	!integer                                   :: i,n,iunit,np,ispe
	!  if (.not.associated(name%num)) return
	!  n = size(name%num)
	!     call open_fname(fname,iunit)
	!  do i=1,n
 !     do ispe=1,nspecie
	!	 np = name%num(i)%btc(ispe)%np
	!	 if (np > 0) write(iunit,2)'zone plane="',i,'"'
	!  2  format(a12,i5,a2) 
	!     close(iunit)
	!	 call print_breakthru_ (name%num(i)%btc(ispe),fname)
	!   end do
	!   end do
 !end subroutine

 
 function nblock_ (name) result (n)
    implicit none
	type(block_vect_cl), intent(in) :: name
    integer                         :: n
	if (associated(name%num)) then
		n = size(name%num)
	else
		n = 0
	end if
 end function
 
 
 end module