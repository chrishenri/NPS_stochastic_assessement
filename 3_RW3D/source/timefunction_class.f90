

module timefunction_class
   
	implicit none

	private
	public :: timefunction_cl     ! class
	public ::                                             & ! methods
              read_timefunction_                          , &
			  evaluate_timefunction_interpolate_          , &
			  evaluate_timefunction_upwind_               , &
              print_time_function_                 
			                
    type timefunction_cl    
         logical           :: action
         character(len=50) :: name
         integer           :: nt
         real*8            :: scale 
         real*8,  pointer  :: time(:)  => null()  !initial time
         real*8,  pointer  :: val(:)   => null()  !final time      
    end type      

	interface evaluate_timefunction_upwind_
	     module procedure evaluate_timefunction_upwind_1,evaluate_timefunction_upwind_2 ;	end interface


    
    contains

   subroutine print_time_function_ (fname,func)
       use gslib, only: open_fname
       implicit none
       type(timefunction_cl), intent(in) :: func
       character(len=*),      intent(in) :: fname
       integer                           :: iunit,it
          call open_fname(fname,iunit)
          if (.not.associated(func%time)) return
          write(iunit,*)
          write(iunit,*) 'time function ..........: '
          write(iunit,*)
          do it=1,func%nt
          write(iunit,*) it, func%time(it),func%val(it)*func%scale
          end do
          write(iunit,*)
          close(iunit)
   end subroutine
    
    subroutine read_timefunction_ (fname,scale,flag,func) !initialize the empty list
       use gslib, only: open_fname_normal
       implicit none
       character(len=*),      intent(in)    :: fname
       real*8,                intent(in)    :: scale
       integer,               intent(in)    :: flag
       type(timefunction_cl), intent(inout) :: func
       integer                              :: iunit,nt,i
       real*8                               :: tt
       
       if ( flag == 0 ) then
          func%action = .FALSE.
          func%nt = 1
          func%scale = scale
          return
       end if
       
       if (flag == 1) then
          func%action =.TRUE.
          call open_fname_normal(fname,iunit)
          read(iunit,*) !heading
          read(iunit,*) func%name          
          read(iunit,*) nt
          func%nt = nt
          func%scale = scale
          allocate(func%time(nt),func%val(nt))
          do i=1,nt
            read(iunit,*) func%time(i),func%val(i)
          end do
          close(iunit)
          return
      end if
       
    end subroutine

    function evaluate_timefunction_interpolate_ (func,time) result(val)
       use gslib, only: locate
       implicit none
       type(timefunction_cl), intent(in) :: func
       real*8,                intent(in) :: time
       integer                           :: it,nt
       real*8                            :: val
       nt = func%nt
       if (.not.associated(func%time)) then
          val = func%scale
          return
       end if
       if (time < func%time(1).or. time > func%time(nt)) then
          val = 0.d0
          return
       end if
       call locate(func%time,nt,1,nt,time,it)
       val = func%val(it) + (func%val(it+1)-func%val(it))/(func%time(it+1)-func%time(it))*(time-func%time(it)) 
       val = val * func%scale 
    end function

   function evaluate_timefunction_upwind_1 (func,time) result(val)
       use gslib, only: locate
       implicit none
       type(timefunction_cl), intent(in) :: func
       real*8,                intent(in) :: time
       integer                           :: it,nt, test
       real*8                            :: val
       nt = func%nt
       if (.not.associated(func%time)) then
          val = func%scale
          return
       end if
       if (time < func%time(1).or. time > func%time(nt)) then
          val = 0.d0
          return
       else if (time == func%time(1)) then
          val = func%val(1)
          return
       end if
       call locate(func%time,nt,1,nt,time,it)
       val = func%val(it) 
       val = val * func%scale 
    end function

   function evaluate_timefunction_upwind_2 (func,it) result(val)
       implicit none
       type(timefunction_cl), intent(in) :: func
       integer,               intent(in) :: it
       integer                           :: nt, test
       real*8                            :: val
       nt = func%nt
       if (.not.associated(func%time)) then
          val = func%scale
          return
       end if
       if (it<=nt) then 
          val = func%val(it)
          val = val * func%scale
       else if(it>nt) then
          val = func%val(nt) ! continue by projection
          val = val * func%scale
       end if 
    end function

end module timefunction_class