
  module dispersion_class
  use array_class
  implicit none

  private
  public :: dispersion_cl  !class
  public ::                              &
            delete_dispersion_         , &
			inquire_homogeneity_disp_  , &
			alloc_dispersion_

  type dispersion_cl
      logical                 :: action
	  type(array_cl), pointer :: aL  => null()
	  type(array_cl), pointer :: aTH => null()
	  type(array_cl), pointer :: aTV => null()
	  type(array_cl), pointer :: dm  => null()
	  type(array_cl), pointer :: dmTH => null()
	  type(array_cl), pointer :: dmTV => null()
	  real*8,         pointer :: MultA(:) => null()  !multiplicator for all dispersivities (species dependent)
	  real*8,         pointer :: MultD(:) => null()  !multiplicator for molecular diffusion (species dependent)
  end type

  contains

  subroutine alloc_dispersion_ (this)
     use global_variables, only: nspe_aq
     implicit none
	 type(dispersion_cl), intent(inout) :: this
         if (.not. associated (this % aL))  allocate (this % aL)
		 if (.not. associated (this % aTH)) allocate (this % aTH)
		 if (.not. associated (this % aTV)) allocate (this % aTV)
		 if (.not. associated (this % dm))  allocate (this % dm)
		 if (.not. associated (this % dmTH))   allocate (this % dmTH)
		 if (.not. associated (this % dmTV))   allocate (this % dmTV)
		 if (.not. associated (this % MultA))  allocate (this % MultA(nspe_aq))
		 if (.not. associated (this % MultD))  allocate (this % MultD(nspe_aq))
  end subroutine

  subroutine delete_dispersion_ (this)
     use array_class
     implicit none
	 type(dispersion_cl), intent(inout) :: this
	    this%action=.FALSE.
         if (associated (this % aL))  then
		      call delete_array_ (this%aL)
		      deallocate (this % aL); endif
		 if (associated (this % aTH)) then
		      call delete_array_ (this%aTH)
		      deallocate (this % aTH); endif
		 if (associated (this % aTV)) then
		      call delete_array_ (this%aTV)
			  deallocate (this % aTV); endif		 
		 if (associated (this % dm)) then
	  	      call delete_array_ (this%dm)	      
			  deallocate (this % dm); endif
		 if (associated (this % dmTH)) then
	  	      call delete_array_ (this%dmTH)	      
			  deallocate (this % dmTH); endif
		 if (associated (this % dmTV)) then
	  	      call delete_array_ (this%dmTV)	      
			  deallocate (this % dmTV); endif
		 if (associated (this % MultA)) deallocate (this % MultA)
		 if (associated (this % MultD)) deallocate (this % MultD)

  end subroutine

  function inquire_homogeneity_disp_ (this) result (flag)
      use array_class
	  implicit none
	  type(dispersion_cl), intent(in)  :: this
	  logical                          :: flag
	  integer                          :: i
	        flag = .TRUE.
            if(associated(this % aL  % values )  .and. length_array_ (this % aL)   > 1)  flag = .FALSE.
			if(associated(this % aTH % values )  .and. length_array_ (this % aTH)  > 1)  flag = .FALSE.
			if(associated(this % aTV % values )  .and. length_array_ (this % aTV)  > 1)  flag = .FALSE.  
            if(associated(this % dm  % values )  .and. length_array_ (this % dm)   > 1)  flag = .FALSE.
            if(associated(this % dmTH  % values )  .and. length_array_ (this % dmTH)   > 1)  flag = .FALSE.
            if(associated(this % dmTV  % values )  .and. length_array_ (this % dmTV)   > 1)  flag = .FALSE.
   end function



  end module
