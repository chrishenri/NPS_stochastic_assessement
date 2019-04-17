


module StressTime_class
   
	implicit none

	private
	public :: StressTime_cl     ! class
	public ::                                             & ! methods
              calculate_StressTimes_                      , &
              find_loc_StressTime_                        , &
              print_transp_stress_times_
			                
    type StressTime_cl    
         integer           :: nt
         integer           :: loc 
         real*8,  pointer  :: time(:)  => null()  !initial time
    end type      

    
    contains

   subroutine print_transp_stress_times_ (fname,StressTime)
       use gslib, only: open_fname
       implicit none
       type(StressTime_cl)               :: StressTime
       character(len=*),   intent(in) :: fname
       integer                        :: iunit,it
          call open_fname(fname,iunit)
          if (.not.associated(StressTime%time)) return
          write(iunit,*)
          write(iunit,*) ' Transport Stress times: '
          write(iunit,*)
          write(iunit,*)  StressTime%time
          write(iunit,*)
          close(iunit)
   end subroutine

    
    function calculate_StressTimes_ (source,advection,tsim) result (StressTime)
       use source_vect_class
       use advection_class
       use velocity_times
       use gslib, only: remove_repeated_values, sortem
       implicit none
       type(source_vect_cl),  intent(in) :: source
       type(advection_cl),    intent(in) :: advection
       type(StressTime_cl)               :: StressTime
       real*8                            :: tsim
       integer                           :: isum,kinj,ninj,nt,nst,itime
       real*8                            :: v1,v2
       real*8, allocatable               :: time(:)

      ninj=size(source%num)
      isum = 0
      do kinj=1,ninj
        if (trim(adjustl(source%num(kinj)%typeinj)) == 'DIRAC') then
             isum = isum + 1
        elseif (source%num(kinj)%typeinj == 'GENERAL') then
          nt = source%num(kinj)%timefunct%nt
          if (nt<1) cycle
          v1  = source%num(kinj)%timefunct%val(1)
          if (v1 /= 0.d0) then
             isum = isum + 1
          end if  
          do itime=2,nt
             v2=source%num(kinj)%timefunct%val(itime)
             if (v2 == v1) then
                continue
             elseif (v2 /= v1 ) then
                isum = isum + 1
             end if
          v1 = v2
          end do
        end if
      end do

!     allocate memory

      nst = isum

      allocate (time(nst+advection%nt+1))
      
      time = 0.d0 !initialize

!     read stress times due to change in injection
 
      isum = 0
      do kinj=1,ninj
        if (trim(adjustl(source%num(kinj)%typeinj)) == 'DIRAC') then
             isum = isum + 1
             time(isum) = source%num(kinj)%TimeStartInj          
        elseif (source%num(kinj)%typeinj == 'GENERAL') then
          nt = source%num(kinj)%timefunct%nt
          if (nt<1) cycle
          v1  = source%num(kinj)%timefunct%val(1)
          if (v1 /= 0.d0) then
             isum = isum + 1
             time(isum) = source%num(kinj)%timefunct%time(1)
          end if  
          do itime=2,nt
             v2=source%num(kinj)%timefunct%val(itime)
             if (v2 == v1) then
                continue
             elseif (v2 /= v1 ) then
                isum = isum + 1
                time(isum) = source%num(kinj)%timefunct%time(itime)
             end if
          v1 = v2
          end do
        end if
      end do

!     add advection times

      time(nst+1:nst+advection%nt+1)=advection%time(0:advection%nt)

!     Last time is the simulation time if stress period are not looped
     
      if (.not.velotime%loop_per) then
          if (time(nst+advection%nt+1)<tsim) time(nst+advection%nt+1) = tsim
      end if

!     remove repeated Stress Times

      call remove_repeated_values (time)

      nst = size(time)

      call sortem (1,nst,time,0,time,time,time,time,time,time,time)

      allocate (StressTime%time(nst))

      StressTime%time = time
      StressTime%nt   = nst
      StressTime%loc  = 1

      deallocate(time)


    end function

    subroutine find_loc_StressTime_ (StressTime,time)
       use velocity_times
       implicit none
       type(StressTime_cl), intent(inout) :: StressTime
       real*8,              intent(in)    :: time
       integer                            :: i    
       
       if (velotime%restart_flux_from_mf2k == .TRUE.) StressTime%loc = 1

       do i=StressTime%loc,StressTime%nt-1
           if (StressTime%time(i)<=time.and.StressTime%time(i+1)>time) then
           StressTime%loc = i 
           return
           end if
       end do
       StressTime%loc = StressTime%nt

    end subroutine


end module StressTime_class