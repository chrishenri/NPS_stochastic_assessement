!**********************************************************************
!  calculates elapsed run times
!**********************************************************************

   module cal_elapsed_time

   contains


!  **********************************************************************************
!  subroutine get start time
!  **********************************************************************************

      subroutine get_start_time (ibdt,fname)
      use gslib, only:open_fname
      
	  implicit none
      
	  character(len=*),optional, intent(in)     :: fname
      character(len=10)    :: chdate, chtime, chzone
      integer, intent(out) :: ibdt(8)
	  integer :: i,unit
!     ------------------------------------------------------------------
!
!     get current date and time, assign to ibdt, and write to screen
!
      if (present(fname)) then
        call open_fname (fname,unit)
	    call date_and_time(chdate,chtime,chzone,ibdt)
        write(unit,*)
		write(unit,2) (ibdt(i),i=1,3),(ibdt(i),i=5,7)
	  end if
      
	  call date_and_time(chdate,chtime,chzone,ibdt)
      write(*,2) (ibdt(i),i=1,3),(ibdt(i),i=5,7)


 2    format(1x,'run start date and time (yyyy/mm/dd hh:mm:ss): ', &
      i4,'/',i2.2,'/',i2.2,1x,i2,':',i2.2,':',i2.2,/)

      end subroutine

!  **********************************************************************************
!  subroutine get end time and calculate elapsed time
!  **********************************************************************************
      subroutine get_end_time (ibdt,fname)
      use gslib, only:open_fname
      
	  implicit none
!
      character(len=10) chdate, chtime, chzone
	  character(len=*),optional, intent(in)     :: fname
      integer, intent(in) :: ibdt(8)
	  integer :: iedt(8), idpm(12)
	  integer :: nspd,i,leap,ndays,nhours,nmins,nrsecs,ibd,ied,mb,me,nm,mc, &
	             m,nsecs,msecs
	  integer :: unit

      real :: rsecs,elsec

      data idpm /31,28,31,30,31,30,31,31,30,31,30,31/  ! days per month
      data nspd /86400/                                ! seconds per day
!     ------------------------------------------------------------------
 1000 format(1x,'run end date and time (yyyy/mm/dd hh:mm:ss): ',i4,'/',i2.2,'/',i2.2,1x,i2,':',i2.2,':',i2.2)
 1010 format(1x,'elapsed run time: ',i3,' days, ',i2,' hours, ',i2,' minutes, ',i2,' seconds',/)
 1020 format(1x,'elapsed run time: ',i2,' hours, ',i2,' minutes, ',i2,' seconds',/)
 1030 format(1x,'elapsed run time: ',i2,' minutes, ',i2,'.',i3.3,' seconds',/)
 1040 format(1x,'elapsed run time: ',i2,'.',i3.3,' seconds',/)

!     get current date and time, assign to iedt, and write to screen
      if (present(fname)) then
        call open_fname (fname,unit)
        write(unit,*)
        write(unit,*)  
        call date_and_time(chdate,chtime,chzone,iedt)
        write(unit,1000) (iedt(i),i=1,3),(iedt(i),i=5,7)
	  end if

      call date_and_time(chdate,chtime,chzone,iedt)
      write(*,1000) (iedt(i),i=1,3),(iedt(i),i=5,7)
!
!     calculate elapsed time in days and seconds
      ndays=0
      leap=0
      if (mod(iedt(1),4).eq.0) leap = 1
      ibd = ibdt(3)            ! begin day
      ied = iedt(3)            ! end day
!     find days
      if (ibdt(2).ne.iedt(2)) then
!       months differ
        mb = ibdt(2)             ! begin month
        me = iedt(2)             ! end month
        nm = me-mb+1             ! number of months to look at
        if (mb.gt.me) nm = nm+12
        mc=mb-1
        do 10 m=1,nm
          mc=mc+1                ! mc is current month
          if (mc.eq.13) mc = 1
          if (mc.eq.mb) then
            ndays = ndays+idpm(mc)-ibd
            if (mc.eq.2) ndays = ndays + leap
          elseif (mc.eq.me) then
            ndays = ndays+ied
          else
            ndays = ndays+idpm(mc)
            if (mc.eq.2) ndays = ndays + leap
          endif
   10   continue
      elseif (ibd.lt.ied) then
!       start and end in same month, only account for days
        ndays = ied-ibd
      endif
      elsec=ndays*nspd
!
!     add or subtract seconds
      elsec = elsec+(iedt(5)-ibdt(5))*3600.0
      elsec = elsec+(iedt(6)-ibdt(6))*60.0
      elsec = elsec+(iedt(7)-ibdt(7))
      elsec = elsec+(iedt(8)-ibdt(8))*0.001
!
!     convert seconds to days, hours, minutes, and seconds
      ndays = elsec/nspd
      rsecs = mod(elsec,86400.0)
      nhours = rsecs/3600.0
      rsecs = mod(rsecs,3600.0)
      nmins = rsecs/60.0
      rsecs = mod(rsecs,60.0)
      nsecs = rsecs
      rsecs = mod(rsecs,1.0)
      msecs = nint(rsecs*1000.0)
      nrsecs = nsecs
      if (rsecs.ge.0.5) nrsecs=nrsecs+1
!
!     write elapsed time to screen

      if (present(fname)) then

        if (ndays.gt.0) then
          write(unit,1010) ndays,nhours,nmins,nrsecs
        elseif (nhours.gt.0) then
          write(unit,1020) nhours,nmins,nrsecs
        elseif (nmins.gt.0) then
          write(unit,1030) nmins,nsecs,msecs
        else
          write(unit,1040) nsecs,msecs
        endif
      end if

      if (ndays.gt.0) then
          write(*,1010) ndays,nhours,nmins,nrsecs
      elseif (nhours.gt.0) then
          write(*,1020) nhours,nmins,nrsecs
      elseif (nmins.gt.0) then
          write(*,1030) nmins,nsecs,msecs
      else
          write(*,1040) nsecs,msecs
      endif

      end subroutine get_end_time



   end module cal_elapsed_time