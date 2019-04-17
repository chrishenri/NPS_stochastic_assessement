!*********************************************************************************************************************
!   Reaction class
!*********************************************************************************************************************
  module reaction_class   ! class defining equilibrium reaction (it is only implemented linear sorption)
  use sorption_class
  use decay_class
  use kinetic_reactions_class
  implicit none

  private
  public :: reaction_cl   !classes
  public ::                                    &       !methods
        read_reaction_                     ,   &
        delete_reaction_                   ,   &
		list_reaction_                     ,   &
		inquire_homogeneity_react_

        type reaction_cl
            logical                               :: action         !flag to distinguish if reactions are being used or not
            integer,                    pointer   :: nspecie        => null()
	        type(sorption_rx_cl),       pointer   :: sorption       => null()
	        type(decay_cl),             pointer   :: decay          => null()
	        type(kinetic_reactions_cl), pointer   :: kinetic        => null()
        end type

  contains

  !**************************** methods
  !************************************
  function read_reaction_ (unit,geo) result (this) 
      use geometry_class
      use sorption_class
      use decay_class
      use kinetic_reactions_class
      use gslib,            only: upper_case 
      use global_variables, only: nspecie
      implicit none

      integer,           intent(in)     :: unit
      type(geometry_cl), intent(in)     :: geo
      type(reaction_cl)                 :: this
      type(decay_cl)                    :: decay
      type(sorption_rx_cl)              :: sorption
      type(kinetic_reactions_cl)        :: kinetic

      allocate (this % nspecie)
      this % nspecie = nspecie
      
      ! Read sorption parameters
      read(unit,*); read(unit,*)
      allocate (this % sorption)
      this % sorption = read_sorption_ (unit,geo)

      ! Read first-order decay network parameters
      read(unit,*); read(unit,*); read(unit,*)
      allocate (this % decay)
      this % decay    = read_decay_    (unit,geo)

      ! Read kinetic reactions parameters
      read(unit,*); read(unit,*); read(unit,*)
      allocate (this % kinetic)
      this % kinetic  = read_kinetic_  (unit,geo)

      if (this % sorption % action) this % action = .TRUE.
      if (this % decay    % action) this % action = .TRUE.
      if (this % kinetic  % action) this % action = .TRUE.
      
  end function

  !************************************
  subroutine delete_reaction_ (this)
      use sorption_class
      use decay_class
      use kinetic_reactions_class
      implicit none

	  type(reaction_cl) :: this

	  if(associated(this%sorption))     call delete_sorption_    (this%sorption)
      if(associated(this%decay))        call delete_decay_       (this%decay)
      if(associated(this%kinetic))      call delete_kinetic_     (this%kinetic)

  end subroutine

  !************************************
  subroutine list_reaction_ (this,fname)
      use sorption_class
      use decay_class
      use kinetic_reactions_class
      implicit none

	  character(len=*), intent(in)                :: fname      
	  type(reaction_cl)                           :: this

	  if(associated(this%sorption))     call list_sorption_    (this%sorption,fname)
      if(associated(this%decay))        call list_decay_       (this%decay,fname)
      if(associated(this%kinetic))      call list_kinetic_     (this%kinetic,fname)

  end subroutine

  !************************************
  function inquire_homogeneity_react_ (this) result (flag)
      use sorption_class
      use decay_class
      use kinetic_reactions_class
	  implicit none

      type(reaction_cl), intent(in) :: this
	  !type(sorption_cl)             :: sorption
	  !type(decay_cl)                :: decay
	  !type(kinetic_reactions_cl)    :: kinetic
	  logical                       :: flag
	  flag = .TRUE.

	  flag = inquire_homogeneity_sorption_   (this%sorption)
	  flag = inquire_homogeneity_decay_      (this%decay)
	  flag = inquire_homogeneity_kinetic_    (this%kinetic)

  end function

!*********************************************************************************************************************
  end module


!*********************************************************************************************************************
!*********************************************************************************************************************
