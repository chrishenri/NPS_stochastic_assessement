  module histogram_class

   implicit none

   private
   public :: histo_cl            !class
   public ::                             & !methods
             assign_histo_PDF_         , &
			 assign_histo_CDF_         , &
			 print_histo_              , &
			 assign_NUE_To_histo_PDF_  , &
			 assign_Ndata_To_histo_PDF_, &
			 set_original_histo_values_

    type histo_cl
       ! Parameters CDF
	   integer,          pointer :: inc => null()    !Used for CDF: increment for printing
       ! Parameters PDF
	   integer,          pointer :: ngrid  => null() !Number bins for PDF
	   integer,          pointer :: ndata  => null()
	   integer,          pointer :: NUE    => null() !Derivative Order
	   character(len=8), pointer :: Kernel => null() !Kernel type for PDF
	   real*8,           pointer :: bw   => null()
	   real*8,           pointer :: dt   => null()   !Bandwidth of bins, bin separation interval 
	   real*8,           pointer :: tmin => null()
	   real*8,           pointer :: tmax => null()   !Time interval to evaluate PDF
	   real*8,           pointer :: bwBTC   => null()
	   real*8,           pointer :: tminBTC => null()
	   real*8,           pointer :: tmaxBTC => null()   !Time interval to evaluate PDF
	end type 

    contains

    subroutine assign_NUE_To_histo_PDF_ (name,NUE)
	     use gslib, only: upper_case
		 implicit none
		 type(histo_cl),     intent(inout) :: name
		 integer,            intent(in)    :: NUE
		   if (.not.associated ( name%NUE  )) allocate ( name%NUE  )
		   name % NUE   = NUE
	end subroutine

    subroutine assign_Ndata_To_histo_PDF_ (name,ndata)
	     use gslib, only: upper_case
		 implicit none
		 type(histo_cl),     intent(inout) :: name
		 integer,            intent(in)    :: ndata
		   if (.not.associated ( name%ndata  )) allocate ( name%ndata  )
		   name % ndata   = ndata
	end subroutine

    subroutine set_original_histo_values_ (histo)
        implicit none
		type(histo_cl),     intent(inout) :: histo
        histo%tminBTC = histo%tmin
        histo%tmaxBTC = histo%tmax
        histo%bwBTC   = histo%bw
    end subroutine

    subroutine assign_histo_PDF_ (name,ngrid,Kernel,bw,tmin,tmax,dt)
	     use gslib, only: upper_case
		 implicit none
		 type(histo_cl),     intent(inout) :: name
		 integer,            intent(in)    :: ngrid
		 real*8,             intent(in)    :: bw,tmax,tmin
		 real*8, optional,   intent(in)    :: dt
		 character(len=*),   intent(in)    :: Kernel
		   if (.not.associated ( name%ngrid  )) allocate ( name%ngrid  )
		   if (.not.associated ( name%Kernel )) allocate ( name%Kernel )
		   if (.not.associated ( name%bw     )) allocate ( name%bw     )
		   if (.not.associated ( name%dt     )) allocate ( name%dt     )
		   if (.not.associated ( name%tmin   )) allocate ( name%tmin   )
		   if (.not.associated ( name%tmax   )) allocate ( name%tmax   )
		   if (.not.associated ( name%bwBTC     )) allocate ( name%bwBTC     )
		   if (.not.associated ( name%tminBTC   )) allocate ( name%tminBTC   )
		   if (.not.associated ( name%tmaxBTC   )) allocate ( name%tmaxBTC   )	   
		   name % ngrid   = ngrid
		   name % Kernel  = upper_case (Kernel)
		   name % bw      = bw
		   name % tmin    = tmin
		   name % tmax    = tmax
		   name % bwBTC      = bw
		   name % tminBTC    = tmin
		   name % tmaxBTC    = tmax
		   if (present(dt))      name % dt = dt
		   if (.not.present(dt)) name % dt = 0.d0
	end subroutine

    subroutine assign_histo_CDF_ (name,increment)
	     implicit none
		 type(histo_cl),     intent(inout) :: name
		 integer,            intent(in)    :: increment
  		   if (.not.associated(name%inc)) allocate ( name%inc ) 
		   name % inc     = increment
	end subroutine

	subroutine print_histo_ (this,fname)
         use gslib, only: open_fname
		 implicit none
		 type(histo_cl),             intent(in) :: this
		 character(len=*), optional, intent(in) :: fname
		 integer                                :: unit

         if (.not.present(fname)) then; unit = 6
		                          else; call open_fname (fname,unit); endif

		 write(unit,*)
		 if (associated(this % inc) ) then
		    write(unit,*) ' Printing Increment for CDF ..........:',this%inc
		 end if
		 write(unit,*)
         if (associated(this % ndata)) then
			write(unit,*) ' Number of data values................: ',this%ndata; end if
         if (associated(this % NUE)) then
			write(unit,*) ' Derivative Order PDF ................: ',this%NUE; end if
		 if (associated(this % ngrid)) then
            write(unit,*) ' Number of bins for PDF ..............: ',this%ngrid
			write(unit,*) ' Kernel Type for PDF .................: ',this%Kernel
			write(unit,*) ' Bandwidth for PDF ...................: ',this%bwBTC
			write(unit,*) ' Separation of Grid-Points of PDF ....: ',this%dt
			write(unit,*) ' Minimum Time.........................: ',this%tminBTC
			write(unit,*) ' Maximum Time.........................: ',this%tmaxBTC
		 end if
		 write(unit,*)		  
             
		  

	end subroutine

  end module