!********************************************************************************************************************
!      multirate mass transfer model
!********************************************************************************************************************
  module multirate_class   ! class defining mass transfer processes
  use array_class
  implicit none

  private
  public :: multirate_mt_cl                          !classes
  public ::                                  &    !methods
            read_multirate_,                 &
            delete_multirate_,               &
            list_multirate_,                 &
            inquire_homogeneity_multirate_

  type multirate_mt_cl
      integer, pointer                  :: nzones => null()      !number of immobile zones
      type(array_cl), pointer           :: poro_imm(:) => null() !array with porosities immobile zones
      type(array_cl), pointer           :: alphap(:)   => null() !array with mass transfer coefficients for immobile zones
  end type

  contains

  !**************************** methods
  !************************************
  function read_multirate_ (mass_trans,unit,geo) result (this) 
      use array_class
      use geometry_class
      use global_variables, only: nzoneim
      implicit none

      integer,           intent(in)     :: unit
      logical,           intent(in)     :: mass_trans
      type(geometry_cl), intent(in)     :: geo
      character(len=200)                :: file
      real*8                            :: const
      integer                           :: ivar,flag,i
      type(multirate_mt_cl)                :: this

      read(unit,*) nzoneim
      if (mass_trans) then
         allocate (this%nzones)
         allocate (this%poro_imm(nzoneim))
         allocate (this%alphap(nzoneim))
         this%nzones = nzoneim ; endif

      do i=1,nzoneim
         read(unit,*) file,const,ivar,flag
         if (mass_trans)  then
             this%poro_imm(i)  = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); end if
      end do
      
      do i=1,nzoneim
         read(unit,*) file,const,ivar,flag
         if (mass_trans)  then
             this%alphap(i)    = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); end if
      end do

      
      
  end function  
  
  !************************************
  subroutine delete_multirate_ (this)
     use array_class
     implicit none

     type(multirate_mt_cl), intent(inout)          :: this

     if(associated(this % nzones))    deallocate(this % nzones)
     if(associated(this % poro_imm))  deallocate(this % poro_imm)
     if(associated(this % alphap))    deallocate(this % alphap)

  end subroutine

  !************************************
  subroutine list_multirate_ (this,fname)
      use gslib, only:open_fname
      use global_variables, only:fdbg
      use array_class
      use gslib, only: integer_to_char
      implicit none

      character(len=*), intent(in)              :: fname
      type (multirate_mt_cl), intent (inout)       :: this
      integer                                   :: i
      character(len=6)                          :: num
      integer                                   :: unit,nzone

      nzone=0
      if (associated( this%nzones ))  nzone = this%nzones
      call open_fname (fdbg,unit)
      write(unit,*)
      write(unit,*) ' MULTIRATE MASS TRANSFER MODEL'
      write(unit,*)
      write(unit,*) ' Number of Immobile Zones............: ',nzone
      write(unit,*)
      close(unit)

      if(associated(this%nzones)) then
          do i=1,this%nzones
            num = integer_to_char (i)
            if(associated(this % poro_imm)) call list_array_ (this%poro_imm(i),fname,'Porosity Immobile Zone '   //trim(num)//'............................................:')
            if(associated(this % alphap))   call list_array_ (this%alphap(i), fname,'Mass Transfer Coefficient '//trim(num)//'.........................................:')
          end do
      end if

  end subroutine

  !************************************
  function inquire_homogeneity_multirate_ (this) result (flag)
      use array_class
      implicit none

      type(multirate_mt_cl), intent(in)        :: this
      logical                               :: flag
      integer                               :: i

      flag = .TRUE.
      if(associated(this % nzones)) then
        do i=1,this%nzones
           if(associated(this % poro_imm) .and. length_array_ (this % poro_imm(i)) > 1)  flag = .FALSE.
           if(associated(this % alphap)   .and. length_array_ (this % alphap(i))   > 1)  flag = .FALSE.
        end do
      end if

   end function

  end module

!********************************************************************************************************************
!      spherical diffusion mass transfer model
!********************************************************************************************************************
 module spherical_diffusion_class   ! class defining mass transfer processes
 use array_class
  implicit none

  private
  public :: spherical_diffusion_mt_cl                        !classes
  public ::                                       &       !methods
            delete_spherical_diffusion_         , &
			list_spherical_diffusion_           , &
			read_spherical_diffusion_           , &
			inquire_homogeneity_spherical_
			 

  type spherical_diffusion_mt_cl
	  integer,        pointer   :: nzones      => null()   !number of immobile zones
	  type(array_cl), pointer   :: poro_imm    => null()   !array with porosities immobile zones
	  type(array_cl), pointer   :: alphap      => null()   !normalized difussion rate coefficient alpha = Dp / a^2
	                                                       ! Dp: is the effective pore diffusion coefficient
														   ! Da: (Harvey and Gorelick notation, 1995) is Dp/Rimm 
														   !  a: distance from the center to the immobile zone 
  end type

  contains

  !**************************** methods
  !************************************
  function read_spherical_diffusion_ (reaction,unit,geo) result (this) 
      use array_class
	  use geometry_class
	  use global_variables, only: nzoneim
	  implicit none
	  
      integer,           intent(in) :: unit
	  logical,           intent(in) :: reaction    
	  type(geometry_cl), intent(in) :: geo
	  character(len=200)            :: file
	  real*8                        :: const
	  integer                       :: ivar,flag
	  type(spherical_diffusion_mt_cl)  :: this
	  
         read(unit,*) nzoneim
		      if (reaction) then
			       allocate (this%nzones)
				   this%nzones=nzoneim; endif
         read(unit,*) file,const,ivar,flag
	          if (reaction)  then
			        allocate (this%poro_imm)
					this%poro_imm = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
			  if (reaction)  then
			        allocate (this%alphap)
			        this%alphap = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif	 
  end function

  !************************************
  subroutine delete_spherical_diffusion_ (this)   !destructor
     use array_class
	 implicit none
     integer                                     :: nzone
	 type(spherical_diffusion_mt_cl), intent(inout) :: this
		if (associated(this%nzones))   deallocate (this%nzones)
		if (associated(this%poro_imm)) deallocate (this%poro_imm)
		if (associated(this%alphap))   deallocate (this%alphap)
  end subroutine
  
  !************************************
  subroutine list_spherical_diffusion_ (this,fname)   !printing mass transfer parameters in a list
	  use gslib, only:open_fname
	  use global_variables, only:fdbg
	  use array_class
	  implicit none
	  character(len=*), intent(in)                :: fname
	  type(spherical_diffusion_mt_cl), intent(inout) :: this
	  integer                                     :: nzone,unit
	             nzone=0
				 if (associated( this%nzones ))  nzone  = this%nzones
                 call open_fname (fdbg,unit)
                 write(unit,*)
				 write(unit,*) ' SPHERICAL DIFFUSION MASS TRANSFER MODEL'
                 write(unit,*)
				 write(unit,*) ' Number of Immobile Zones............: ',nzone
				 write(unit,*)
				 close (unit)

                 if (associated(this%poro_imm)) call list_array_ (this%poro_imm,fname, 'Poro Immobile............:')
				 if (associated(this%alphap))   call list_array_ (this%alphap,fname,  'Mass Transf Coeff........:')
  end subroutine

  !************************************
  function inquire_homogeneity_spherical_ (this) result (flag)
      use array_class
	  implicit none
	  type(spherical_diffusion_mt_cl), intent(in) :: this
	  logical                                  :: flag
	        flag = .TRUE.
            if(associated(this % poro_imm) .and. length_array_ (this % poro_imm) > 1)  flag = .FALSE.
			if(associated(this % alphap)   .and. length_array_ (this % alphap)   > 1)  flag = .FALSE.
   end function

  end module

!********************************************************************************************************************
!      layered diffusion mass transfer model
!********************************************************************************************************************
 module layered_diffusion_class   ! class defining mass transfer processes
 use array_class
  implicit none

  private
  public :: layered_diffusion_mt_cl                        !classes
  public ::                                     &       !methods
            delete_layered_diffusion_         , &
			list_layered_diffusion_           , &
			read_layered_diffusion_           , &
			inquire_homogeneity_layered_


  type layered_diffusion_mt_cl
	  integer,        pointer   :: nzones   => null()   !number of immobile zones
	  type(array_cl), pointer   :: poro_imm => null()   !porosity immobile zone
	  type(array_cl), pointer   :: alphap   => null()   !normalized difussion rate coefficient alpha = Dp / a^2
	                                                    ! Dp: is the effective pore diffusion coefficient
														! Da: (Harvey and Gorelick notation, 1995) is Dp/Rimm
														!  a: distance from the center to the immobile zone 
  end type

  contains

  !**************************** methods
  !************************************
  function read_layered_diffusion_ (reaction,unit,geo) result (this) 
      use array_class
	  use geometry_class
	  use global_variables, only: nzoneim
	  implicit none

      integer,           intent(in) :: unit
	  logical,           intent(in) :: reaction
	  type(geometry_cl), intent(in) :: geo    
	  character(len=200)            :: file
	  real*8                        :: const
	  integer                       :: ivar,flag
	  type(layered_diffusion_mt_cl)    :: this

         read(unit,*) nzoneim
		      if (reaction) then
			       allocate (this%nzones)
				   this%nzones=nzoneim; endif
         read(unit,*) file,const,ivar,flag
	          if (reaction)  then
			        allocate (this%poro_imm)
					this%poro_imm = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
			  if (reaction)  then
			        allocate (this%alphap)
			        this%alphap = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif	 
  end function

  !************************************
  subroutine delete_layered_diffusion_ (this)   !destructor
     use array_class
	 implicit none

	 type(layered_diffusion_mt_cl), intent(inout)   :: this

     if (associated(this%nzones))   deallocate (this%nzones)
	 if (associated(this%poro_imm)) deallocate (this%poro_imm)
	 if (associated(this%alphap))   deallocate (this%alphap)
		
  end subroutine

  !************************************
  subroutine list_layered_diffusion_ (this,fname)   !printing mass transfer parameters in a list
	  use gslib, only:open_fname
	  use global_variables, only:fdbg
	  use array_class
	  implicit none
	  character(len=*), intent(in)              :: fname
	  type(layered_diffusion_mt_cl), intent(inout) :: this
	  integer                                   :: unit,nzone
	  
	  nzone=0
	  if (associated( this%nzones ))  nzone  = this%nzones
	  call open_fname (fdbg,unit)
	  write(unit,*)
	  write(unit,*) ' LAYERED DIFFUSION MASS TRANSFER MODEL'
	  write(unit,*)
	  write(unit,*) ' Number of Immobile Zones............: ',nzone
	  write(unit,*)
	  close (unit)

	  if (associated(this%poro_imm)) call list_array_ (this%poro_imm,fname, 'Poro Immobile............:')
	  if (associated(this%alphap))   call list_array_ (this%alphap,fname,  'Mass Transf Coeff........:')

  end subroutine

  !************************************
  function inquire_homogeneity_layered_ (this) result (flag)
      use array_class
	  implicit none
	  type(layered_diffusion_mt_cl), intent(in) :: this
	  logical                                :: flag
	  
	  flag = .TRUE.
      if(associated(this % poro_imm) .and. length_array_ (this % poro_imm) > 1)  flag = .FALSE.
	  if(associated(this % alphap)   .and. length_array_ (this % alphap)   > 1)  flag = .FALSE.

   end function

  end module


!********************************************************************************************************************
!      Cylindrical diffusion mass transfer model
!********************************************************************************************************************
 module cylindrical_diffusion_class   ! class defining mass transfer processes
 use array_class
 implicit none

  private
  public :: cylindrical_diffusion_mt_cl                        !classes
  public ::                                         &       !methods
            delete_cylindrical_diffusion_         , &
			list_cylindrical_diffusion_           , &
			read_cylindrical_diffusion_           , &
			inquire_homogeneity_cylindrical_
			 

  type cylindrical_diffusion_mt_cl
	  integer,        pointer   :: nzones   => null()   !number of immobile zones
	  type(array_cl), pointer   :: poro_imm => null()   !porosity immobile zone
	  type(array_cl), pointer   :: alphap   => null()   !normalized difussion rate coefficient alpha = Dp / a^2
	                                                    ! Dp: is the effective pore diffusion coefficient
														! Da: (Harvey and Gorelick notation, 1995) is Dp/Rimm
														!  a: distance from the center to the immobile zone 
  end type

  contains

  !**************************** methods
  !************************************
  function read_cylindrical_diffusion_ (reaction,unit,geo) result (this) 
      use array_class
	  use geometry_class
	  use global_variables, only: nzoneim
	  implicit none

      integer,           intent(in)     :: unit
	  logical,           intent(in)     :: reaction 
	  type(geometry_cl), intent(in)     :: geo   
	  character(len=200)                :: file
	  real*8                            :: const
	  integer                           :: ivar,flag
	  type(cylindrical_diffusion_mt_cl)    :: this
	  
         read(unit,*) nzoneim
		      if (reaction) then
			       allocate (this%nzones)
				   this%nzones=nzoneim; endif
         read(unit,*) file,const,ivar,flag
	          if (reaction)  then
			        allocate (this%poro_imm)
					this%poro_imm = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
			  if (reaction)  then
			        allocate (this%alphap)
			        this%alphap = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif

  end function

  !************************************
  subroutine delete_cylindrical_diffusion_ (this)   !destructor
     use array_class
	 implicit none
     integer                                       :: nzone
	 type(cylindrical_diffusion_mt_cl), intent(inout) :: this
		if (associated(this%nzones))   deallocate (this%nzones)
		if (associated(this%poro_imm)) deallocate (this%poro_imm)
		if (associated(this%alphap))   deallocate (this%alphap)
  end subroutine

  !************************************
  subroutine list_cylindrical_diffusion_ (this,fname)   !printing mass transfer parameters in a list
	  use gslib, only: open_fname
	  use global_variables, only: fdbg
	  use array_class
	  implicit none
	  character(len=*), intent(in)                  :: fname
	  type(cylindrical_diffusion_mt_cl), intent(inout) :: this
	  integer                                       :: nzone,unit

	  nzone=0
	  if (associated( this%nzones ))  nzone  = this%nzones
	  call open_fname (fdbg,unit)
	  write(unit,*)
	  write(unit,*) ' CYLINDRICAL DIFFUSION MASS TRANSFER MODEL'
	  write(unit,*)
	  write(unit,*) ' Number of Immobile Zones............: ',nzone
	  write(unit,*)
	  close (unit)

	  if (associated(this%poro_imm)) call list_array_ (this%poro_imm,fname, 'Poro Immobile............:')
	  if (associated(this%alphap))   call list_array_ (this%alphap,fname,  'Mass Transf Coeff........:')

  end subroutine

  !************************************
  function inquire_homogeneity_cylindrical_ (this) result (flag)
      use array_class
	  implicit none
	  type(cylindrical_diffusion_mt_cl), intent(in) :: this
	  logical                                  :: flag
	  integer                                  :: i
	        flag = .TRUE.
            if(associated(this % poro_imm) .and. length_array_ (this % poro_imm) > 1)  flag = .FALSE.
			if(associated(this % alphap)   .and. length_array_ (this % alphap)   > 1)  flag = .FALSE.
   end function

  end module

!********************************************************************************************************************
!      power truncated mass transfer model
!********************************************************************************************************************
 module power_distribution_class   ! class defining mass transfer processes
 use array_class
  implicit none

  private
  public :: power_distribution_mt_cl                        !classes
  public ::                                      &       !methods
            delete_power_distribution_         , &
            list_power_distribution_           , &
            read_power_distribution_           , &
            inquire_homogeneity_power_


  type power_distribution_mt_cl
      integer,        pointer   :: nzones   => null()   !number of immobile zones
      type(array_cl), pointer   :: btot     => null()   !total capacity coefficient
      type(array_cl), pointer   :: Amax     => null()   !maximum mass transfer coeff
      type(array_cl), pointer   :: Amin     => null()   !minimum mass transfer coeff
      type(array_cl), pointer   :: power    => null()   !exponent
  end type

  contains

  !**************************** methods
  !************************************
  function read_power_distribution_ (reaction,unit,geo) result (this)
      use array_class
      use geometry_class
      use global_variables, only: nzoneim
      implicit none
      integer,           intent(in) :: unit
      logical,           intent(in) :: reaction
      type(geometry_cl), intent(in) :: geo
      character(len=200)            :: file
      real*8                        :: const
      integer                       :: ivar,flag
      type(power_distribution_mt_cl)   :: this
         read(unit,*) nzoneim
              if (reaction) then
                   allocate (this%nzones)
                   this%nzones=nzoneim; endif
         read(unit,*) file,const,ivar,flag
              if (reaction)  then
                    allocate (this%btot)
                    this%btot = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
              if (reaction)  then
                    allocate (this%Amin)
                    this%Amin = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
              if (reaction)  then
                    allocate (this%Amax)
                    this%Amax = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
              if (reaction)  then
                    allocate (this%power)
                    this%power = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
  end function

  !************************************
  subroutine delete_power_distribution_ (this)   !destructor
     use array_class
     implicit none
     integer                                     :: nzone
     type(power_distribution_mt_cl), intent(inout)  :: this
        if (associated(this%nzones))     deallocate (this%nzones)
        if (associated(this%btot))       deallocate (this%btot)
        if (associated(this%Amin))       deallocate (this%Amin)
        if (associated(this%Amax))       deallocate (this%Amax)
        if (associated(this%power))      deallocate (this%power)
  end subroutine

  !************************************
  subroutine list_power_distribution_ (this,fname)   !printing mass transfer parameters in a list
      use gslib, only:open_fname
      use global_variables, only:fdbg
      use array_class
      implicit none
      character(len=*),            intent(in)    :: fname
      type(power_distribution_mt_cl), intent(inout) :: this
      integer                                    :: unit,nzone
                 nzone=0
                 if (associated( this%nzones ))  nzone  = this%nzones
                 call open_fname (fdbg,unit)
                 write(unit,*)
                 write(unit,*) ' TRUNCATED POWER DISTRIBUTION OF MASS TRANSFER COEFFICIENTS'
                 write(unit,*)
                 write(unit,*) ' Number of Immobile Zones............: ',nzone
                 write(unit,*)
                 close (unit)

                 if (associated(this%btot))     call list_array_ (this%btot, fname, 'Total Capacity................:')
                 if (associated(this%Amin))     call list_array_ (this%Amin, fname, 'Mimim Mass Transfer Coeff.....:')
                 if (associated(this%Amax))     call list_array_ (this%Amax, fname, 'Maxim Mass Transfer Coeff.....:')
                 if (associated(this%power))    call list_array_ (this%power,fname, 'Power Coefficient.............:')
  end subroutine

  !************************************
  function inquire_homogeneity_power_ (this) result (flag)
      use array_class
      implicit none
      type(power_distribution_mt_cl), intent(in)  :: this
      logical                                  :: flag
      integer                                  :: i
            flag = .TRUE.
            if(associated(this % btot)       .and. length_array_ (this % btot)     > 1)  flag = .FALSE.
            if(associated(this % Amin)       .and. length_array_ (this % Amin)     > 1)  flag = .FALSE.
            if(associated(this % Amax)       .and. length_array_ (this % Amax)     > 1)  flag = .FALSE.
            if(associated(this % power)      .and. length_array_ (this % power)    > 1)  flag = .FALSE.
   end function


  end module

!********************************************************************************************************************
!      lognormal distribution of mass transfer coefficient
!********************************************************************************************************************
 module lognormal_distribution_class   ! class defining mass transfer processes
 use array_class
  implicit none

  private
  public :: lognormal_distribution_mt_cl                        !classes
  public ::                                          &       !methods
            delete_lognormal_distribution_         , &
			list_lognormal_distribution_           , &
			read_lognormal_distribution_           , &
			inquire_homogeneity_lognormal_
			 

  type lognormal_distribution_mt_cl
	  integer,        pointer   :: nzones   => null()   !number of immobile zones
	  type(array_cl), pointer   :: btot     => null()   !total capacity coefficient
	  type(array_cl), pointer   :: mean     => null()   !mean of the natural log of mass transfer coeff
	  type(array_cl), pointer   :: stdv     => null()   !standard deviation of the natural log of mass transfer coeff
  end type

  contains

  !**************************** methods
  !************************************
  function read_lognormal_distribution_ (reaction,unit,geo) result (this) 
      use array_class
	  use geometry_class
	  use global_variables, only: nzoneim
	  implicit none
      integer,           intent(in)     :: unit
	  logical,           intent(in)     :: reaction    
	  type(geometry_cl), intent(in)     :: geo
	  character(len=200)                :: file
	  real*8                            :: const
	  integer                           :: ivar,flag,nzones
	  type(lognormal_distribution_mt_cl)   :: this
         read(unit,*) nzoneim
		      if (reaction) then
			       allocate (this%nzones)
				   this%nzones=nzoneim; endif
         read(unit,*) file,const,ivar,flag
	          if (reaction)  then
			        allocate (this%btot)
					this%btot = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
	          if (reaction)  then
			        allocate (this%mean) 
					this%mean = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
         read(unit,*) file,const,ivar,flag
	          if (reaction)  then
			        allocate (this%stdv)
					this%stdv = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); endif
  end function

  !************************************
  subroutine delete_lognormal_distribution_ (this)   !destructor
     use array_class
	 implicit none
     integer                                     :: nzone
	 type(lognormal_distribution_mt_cl), intent(inout)  :: this
		if (associated(this%nzones))     deallocate (this%nzones)
		if (associated(this%btot))       deallocate (this%btot)
		if (associated(this%mean))       deallocate (this%mean)
		if (associated(this%stdv))       deallocate (this%stdv)
  end subroutine

  !************************************
  subroutine list_lognormal_distribution_ (this,fname)   !printing mass transfer parameters in a list
	  use gslib, only:open_fname
	  use global_variables, only:fdbg
	  use array_class
	  implicit none
	  character(len=*),                intent(in)    :: fname
	  type(lognormal_distribution_mt_cl), intent(inout) :: this
	  integer                                        :: unit,nzone
	             nzone=0
				 if (associated( this%nzones ))  nzone  = this%nzones
                 call open_fname (fdbg,unit)
                 write(unit,*)
				 write(unit,*) ' LOGNORMAL DISTRIBUTION OF MASS TRANSFER COEFFICIENTS'
                 write(unit,*)
				 write(unit,*) ' Number of Immobile Zones............: ',nzone
				 write(unit,*)
				 close (unit)

                 if (associated(this%btot))         call list_array_ (this%btot, fname, 'Total Capacity.................:')
		       	 if (associated(this%mean))         call list_array_ (this%mean, fname, 'Mean ln Alpha .................:')
                 if (associated(this%stdv))         call list_array_ (this%stdv, fname, 'Standard Deviation ln Alpha....:')
  end subroutine

  !************************************
  function inquire_homogeneity_lognormal_ (this) result (flag)
      use array_class
	  implicit none
	  type(lognormal_distribution_mt_cl), intent(in)  :: this
	  logical                                  :: flag
	  integer                                  :: i
	        flag = .TRUE.
            if(associated(this % btot)       .and. length_array_ (this % btot)     > 1)  flag = .FALSE.
			if(associated(this % mean)       .and. length_array_ (this % mean)     > 1)  flag = .FALSE.  
            if(associated(this % stdv)       .and. length_array_ (this % stdv)     > 1)  flag = .FALSE.
   end function

  end module

!********************************************************************************************************************
!      Composite Media Mass Transfer Model
!********************************************************************************************************************
 module composite_media_class   ! class defining mass transfer processes
  use multirate_class
  use spherical_diffusion_class
  use layered_diffusion_class
  use cylindrical_diffusion_class
  implicit none

  private
  public :: composite_media_mt_cl                            !classes
  public ::                                         &     !methods
            delete_composite_media_               , &
			list_composite_media_                 , &
			read_composite_media_                 , &
			inquire_homogeneity_composite_       		 

  type composite_media_mt_cl
      integer,                        pointer :: nzones !total number of immobile zones
	  integer,                        pointer :: nmrate,nsph,ncyl,nlay  !number of multirate,spherical,cylindrical, and layered media
      real*8,                         pointer :: Fmrate(:) => null()    !volume fraction for multirate
	  real*8,                         pointer :: Fsph(:)   => null()    !volume fraction for spherical
	  real*8,                         pointer :: Fcyl(:)   => null()    !volume fraction for spherical
	  real*8,                         pointer :: Flay(:)   => null()    !volume fraction for spherical
	  type(multirate_mt_cl),             pointer :: Multirate(:)            => null()
	  type(spherical_diffusion_mt_cl),   pointer :: SphericalDiffusion(:)   => null()
	  type(layered_diffusion_mt_cl),     pointer :: LayeredDiffusion(:)     => null()
	  type(cylindrical_diffusion_mt_cl), pointer :: CylindricalDiffusion(:) => null()
  end type

  contains

  !**************************** methods
  !************************************
  function read_composite_media_ (reaction,unit,geo) result (this) 
      use geometry_class
	  use multirate_class
      use spherical_diffusion_class
      use layered_diffusion_class
      use cylindrical_diffusion_class
      use global_variables, only: nzoneim
	  implicit none
      integer,           intent(in) :: unit
	  logical,           intent(in) :: reaction    
	  type(geometry_cl), intent(in) :: geo
	  integer                       :: i,nmrate,nsph,ncyl,nlay
	  real*8, pointer               :: Fmrate(:),Fsph(:),Fcyl(:),Flay(:)
	  real*8                        :: wg = 0.d0
	  type(composite_media_mt_cl)      :: this

         !read number of zone for multirate, sph., cyl. and lay. diffusion models
         read(unit,*) nmrate,nsph,ncyl,nlay
		 if (reaction) then
			allocate (this%nmrate,this%nsph,this%ncyl,this%nlay)
			this%nmrate = nmrate
			this%nsph   = nsph
			this%ncyl   = ncyl
			this%nlay   = nlay; endif

         !read fraction of each zones for each mass transfer models (don't read if number of zone is 0)
		 allocate (Fmrate(nmrate),Fsph(nsph),Fcyl(ncyl),Flay(nlay))
		 if(nmrate>0) read(unit,*) ( Fmrate(i), i=1,nmrate )
		 if(nsph>0)   read(unit,*) ( Fsph(i),   i=1,nsph )
		 if(ncyl>0)   read(unit,*) ( Fcyl(i),   i=1,ncyl )
		 if(nlay>0)   read(unit,*) ( Flay(i),   i=1,nlay )
         wg = sum(Fmrate) + sum(Fsph) + sum(Fcyl) + sum(Flay)
		 Fmrate = Fmrate/wg
		 Fsph   = Fsph/wg
		 Fcyl   = Fcyl/wg
		 Flay   = Flay/wg

		 if (reaction) then
			allocate (this%Fmrate(nmrate),this%Fsph(nsph),this%Fcyl(ncyl),this%Flay(nlay))
			this%Fmrate = Fmrate
			this%Fsph   = Fsph
			this%Fcyl   = Fcyl
			this%Flay   = Flay; endif
		 deallocate (Fmrate,Fsph,Fcyl,Flay)

		 !read parameters for multirate MT model
		 if (nmrate > 0) then
            allocate (this%Multirate(nmrate))
            do i=1,nmrate
		        this%Multirate(i) = read_multirate_ (reaction,unit,geo); end do
		 endif
		 !read parameters for spherical diffusion MT model
		 if (nsph > 0) then
            allocate (this%SphericalDiffusion(nsph))
            do i=1,nsph
		        this%SphericalDiffusion(i) = read_spherical_diffusion_ (reaction,unit,geo); end do
		 endif
		 !read parameters for cylindrical diffusion MT model
		 if (ncyl > 0) then
            allocate (this%CylindricalDiffusion(ncyl))
            do i=1,ncyl
		        this%CylindricalDiffusion(i) = read_cylindrical_diffusion_ (reaction,unit,geo); end do
		 endif
		 !read parameters for layered MT model
		 if (nlay > 0) then
            allocate (this%LayeredDiffusion(nlay))
            do i=1,nlay
		        this%LayeredDiffusion(i) = read_layered_diffusion_ (reaction,unit,geo); end do
		 endif
		 !get total number of immobile zones 
		 if (nmrate+nsph+ncyl+nlay > 0) then
		    allocate (this%nzones)
			this%nzones = nzone_composite_ (this)
			nzoneim = this%nzones
		 end if
  end function  

  !************************************
  subroutine delete_composite_media_ (this)   !destructor
  use multirate_class
  use spherical_diffusion_class
  use layered_diffusion_class
  use cylindrical_diffusion_class
	 implicit none
     integer                                 :: i
	 type(composite_media_mt_cl), intent(inout) :: this
                 if (associated(this%Multirate))  then
				     do i=1,this%nmrate
					         call delete_multirate_  ( this%Multirate(i) ); enddo
					 deallocate ( this%Multirate )
				 end if
		       	 if (associated(this%SphericalDiffusion)) then
				     do i=1,this%nsph 
					         call delete_spherical_diffusion_   ( this%SphericalDiffusion(i) ); enddo
				     deallocate ( this%SphericalDiffusion )
				 end if
                 if (associated(this%LayeredDiffusion)) then
				     do i=1,this%nlay
					        call delete_layered_diffusion_     ( this%LayeredDiffusion(i) ); enddo
				     deallocate ( this%LayeredDiffusion )
				 end if
                 if (associated(this%CylindricalDiffusion)) then
				     do i=1,this%ncyl
					        call delete_cylindrical_diffusion_ (this%CylindricalDiffusion(i)); enddo
				     deallocate ( this%CylindricalDiffusion )
				 end if
				 if (associated( this%nzones ))  deallocate (this%nzones)
				 if (associated( this%nmrate ))  deallocate (this%nmrate)
				 if (associated( this%nsph ))    deallocate (this%nsph)
				 if (associated( this%ncyl ))    deallocate (this%ncyl)
				 if (associated( this%nlay ))    deallocate (this%nlay)
				 if (associated( this%Fmrate ))  deallocate (this%Fmrate)
				 if (associated( this%Fsph ))    deallocate (this%Fsph)
				 if (associated( this%Fcyl ))    deallocate (this%Fcyl)
				 if (associated( this%Flay ))    deallocate (this%Flay)
  end subroutine

  !************************************
  subroutine list_composite_media_ (this,fname)   !printing mass transfer parameters in a list
  use gslib, only: open_fname
  use global_variables, only: fdbg
  use multirate_class
  use spherical_diffusion_class
  use layered_diffusion_class
  use cylindrical_diffusion_class
	  implicit none
	  character(len=*), intent(in)            :: fname
	  type(composite_media_mt_cl), intent(inout) :: this
	  integer                                 :: i,nzone,nmrate,nsph,ncyl,nlay,unit

                 nzone=0; nmrate=0; nsph=0; ncyl=0; nlay=0
				 if (associated( this%nzones ))  nzone  = this%nzones
				 if (associated( this%nmrate ))  nmrate = this%nmrate
				 if (associated( this%nsph ))    nsph   = this%nsph
				 if (associated( this%ncyl ))    ncyl   = this%ncyl
				 if (associated( this%nlay ))    nlay   = this%nlay

                 call open_fname (fdbg,unit)
                	 
				 write(unit,*)
				 write(unit,*) ' Number of Immobile Zones............: ',nzone
				 write(unit,*) ' Number of Multirate Geometries......: ',nmrate
				 write(unit,*) ' Number of Spherical Geometries......: ',nsph
				 write(unit,*) ' Number of Cylindrical Geometries....: ',ncyl
				 write(unit,*) ' Number of Layered Geometries........: ',nlay
				 write(unit,*)

                 if (associated(this%Fmrate)) then
				     do i=1,this%nmrate
					    write(unit,*) ' Volume Fraction of Multirate  ',i,'........: ', this%Fmrate(i); enddo
				 end if
				 if (associated(this%Fsph)) then
				     do i=1,this%nsph
					    write(unit,*) ' Volume Fraction of Spherical  ',i,'........: ', this%Fsph(i); enddo
				 end if
				 if (associated(this%Fcyl)) then
				     do i=1,this%ncyl
					    write(unit,*) ' Volume Fraction of Cylindrical',i,'........: ', this%Fcyl(i); enddo
				 end if
				 if (associated(this%Flay)) then
				     do i=1,this%nlay
					    write(unit,*) ' Volume Fraction of Layered    ',i,'........: ', this%Flay(i); enddo
				 end if

				 close (unit)

                 if (associated(this%Multirate))  then
				     do i=1,this%nmrate
					         call list_multirate_  ( this%Multirate(i), fname ); enddo
				 end if
		       	 if (associated(this%SphericalDiffusion)) then
				     do i=1,this%nsph 
					         call list_spherical_diffusion_   ( this%SphericalDiffusion(i), fname ); enddo
				 end if
                 if (associated(this%LayeredDiffusion)) then
				     do i=1,this%nlay
					        call list_layered_diffusion_     ( this%LayeredDiffusion(i), fname ); enddo
				 end if
                 if (associated(this%CylindricalDiffusion)) then
				     do i=1,this%ncyl
					        call list_cylindrical_diffusion_ (this%CylindricalDiffusion(i),fname); enddo
				 end if
  end subroutine

  !************************************
  function inquire_homogeneity_composite_ (this) result (flag)
      use multirate_class
      use spherical_diffusion_class
      use layered_diffusion_class
      use cylindrical_diffusion_class
	  implicit none
	  type(composite_media_mt_cl), intent(in)     :: this
	  logical                                  :: flag
	  integer                                  :: i
	        flag = .TRUE.
            if (associated(this%Multirate))  then
				     do i=1,this%nmrate
					         flag = inquire_homogeneity_multirate_  ( this%Multirate(i) ); enddo
			end if
		    if (associated(this%SphericalDiffusion)) then
				     do i=1,this%nsph 
					         flag = inquire_homogeneity_spherical_  ( this%SphericalDiffusion(i) ); enddo
			end if
            if (associated(this%LayeredDiffusion)) then
				     do i=1,this%nlay
					        flag = inquire_homogeneity_layered_     ( this%LayeredDiffusion(i) ); enddo
			end if
            if (associated(this%CylindricalDiffusion)) then
				     do i=1,this%ncyl
					        flag = inquire_homogeneity_cylindrical_ (this%CylindricalDiffusion(i)); enddo
			end if
     end function

  !************************************
  function nzone_composite_ (this) result (nzone)
     implicit none
	 type (composite_media_mt_cl), intent(in) :: this
	 integer                               :: nzone,i
                 nzone = 0
				 if (associated(this%Multirate))  then
				     do i=1,this%nmrate
							 nzone   = nzone + this%Multirate(i)%nzones
					 enddo
				 end if
		       	 if (associated(this%SphericalDiffusion)) then
				     do i=1,this%nsph 
					         nzone   = nzone + this%SphericalDiffusion(i)%nzones; enddo
				 end if
                 if (associated(this%LayeredDiffusion)) then
				     do i=1,this%nlay
					       nzone   = nzone + this%LayeredDiffusion(i)%nzones; enddo
				 end if
                 if (associated(this%CylindricalDiffusion)) then
				     do i=1,this%ncyl
					        nzone   = nzone + this%CylindricalDiffusion(i)%nzones; enddo
				 end if
  end function

  end module

!*********************************************************************************************************************
!*********************************************************************************************************************
!   Mass Transfer class
!*********************************************************************************************************************
  module mass_trans_class
  use multirate_class
  use spherical_diffusion_class
  use layered_diffusion_class
  use cylindrical_diffusion_class
  use composite_media_class
  use power_distribution_class
  use lognormal_distribution_class

  implicit none

  private
  public :: mass_trans_cl               !classes
  public ::                                     &       !methods
            read_mass_trans_                ,   &
            list_mass_trans_                ,   &
            delete_mass_trans_              ,   &
            inquire_homogeneity_mass_trans_

       type mass_trans_cl
           logical                                     :: action         !flag to distinguish if sorption is being used or not
           character(len=30)                           :: type_mass_trans
           type(multirate_mt_cl),               pointer   :: multirate            => null()
           type(spherical_diffusion_mt_cl),     pointer   :: sphericalDiff        => null()
           type(layered_diffusion_mt_cl),       pointer   :: layeredDiff          => null()
           type(cylindrical_diffusion_mt_cl),   pointer   :: cylindricalDiff      => null()
           type(power_distribution_mt_cl),      pointer   :: power                => null()
           type(lognormal_distribution_mt_cl),  pointer   :: logNormal            => null()
           type(composite_media_mt_cl),         pointer   :: composite            => null()
       end type
  
  contains

  !**************************** methods
  !************************************
  function read_mass_trans_ (unit,geo) result (this) 
      use geometry_class
      use multirate_class
      use spherical_diffusion_class
      use layered_diffusion_class
      use cylindrical_diffusion_class
      use composite_media_class
      use power_distribution_class
      use lognormal_distribution_class
      use to_solve,             only: mass_transACTION, mass_transTYPE
      use gslib,                only: upper_case 
      use global_variables,     only: nzoneim
      implicit none

      integer,           intent(in) :: unit
      type(geometry_cl), intent(in) :: geo
      type(mass_trans_cl)           :: this

      read(unit,*) this % action
      mass_transACTION  = this % action

      read(unit,*) this % type_mass_trans
      this % type_mass_trans = upper_case (this % type_mass_trans)
      mass_transTYPE    = this % type_mass_trans

      select case (this % type_mass_trans)

            case ('MULTIRATE')
                 allocate (this % multirate)
                 this % multirate           = read_multirate_ (this % action,unit,geo)
            case ('SPHERICAL_DIFFUSION')
                 allocate (this % sphericalDiff)
                 this % sphericalDiff       = read_spherical_diffusion_ (this % action,unit,geo)
            case ('LAYERED_DIFFUSION')
                 allocate (this % layeredDiff)
                 this % layeredDiff         = read_layered_diffusion_ (this % action,unit,geo)
            case ('CYLINDRICAL_DIFFUSION')
                 allocate (this % cylindricalDiff)
                 this % cylindricalDiff     = read_cylindrical_diffusion_ (this % action,unit,geo)
            case ('POWER_LAW')
                 allocate (this % power)
                 this % power               = read_power_distribution_ (this % action,unit,geo)
            case ('LOGNORMAL_LAW')
                 allocate (this % logNormal)
                 this % LogNormal           = read_lognormal_distribution_ (this % action,unit,geo)
            case ('COMPOSITE_MEDIA')
                 allocate (this % composite)
                 this % Composite           = read_composite_media_ (this % action,unit,geo)

      case default
            if (this % action == .TRUE. ) then 
                stop 'Mass Transfer type not reconnized'
            end if

      end select

  end function

  !************************************
  subroutine delete_mass_trans_ (this)
     use multirate_class
     use spherical_diffusion_class
     use layered_diffusion_class
     use cylindrical_diffusion_class
     use composite_media_class
     use power_distribution_class
     use lognormal_distribution_class
     implicit none

     type(mass_trans_cl) :: this

     if(associated(this%multiRate))            call delete_multirate_              (this%multiRate)
     if(associated(this%sphericalDiff))        call delete_spherical_diffusion_    (this%sphericalDiff)
     if(associated(this%layeredDiff))          call delete_layered_diffusion_      (this%layeredDiff)
     if(associated(this%cylindricalDiff))      call delete_cylindrical_diffusion_  (this%cylindricalDiff)
     if(associated(this%power))                call delete_power_distribution_     (this%power)
     if(associated(this%composite))            call delete_composite_media_        (this%composite)
     if(associated(this%logNormal))            call delete_lognormal_distribution_ (this%logNormal)

  end subroutine


  !************************************
  subroutine list_mass_trans_ (this,fname)
      use multirate_class
      use spherical_diffusion_class
      use layered_diffusion_class
      use cylindrical_diffusion_class
      use composite_media_class
      use power_distribution_class
      use lognormal_distribution_class
      implicit none

      character(len=*), intent(in)                :: fname
      type(mass_trans_cl)                         :: this

      if(associated(this%multiRate))            call list_multirate_              (this%multiRate,fname)
      if(associated(this%sphericalDiff))        call list_spherical_diffusion_    (this%sphericalDiff,fname)
      if(associated(this%layeredDiff))          call list_layered_diffusion_      (this%layeredDiff,fname)
      if(associated(this%cylindricalDiff))      call list_cylindrical_diffusion_  (this%cylindricalDiff,fname)
      if(associated(this%power))                call list_power_distribution_     (this%power,fname)
      if(associated(this%composite))            call list_composite_media_        (this%composite,fname)
      if(associated(this%logNormal))            call list_lognormal_distribution_ (this%logNormal,fname)

  end subroutine

  !************************************
  function inquire_homogeneity_mass_trans_   (this) result (flag)
      use multirate_class
      use spherical_diffusion_class
      use layered_diffusion_class
      use cylindrical_diffusion_class
      use composite_media_class
      use power_distribution_class
      use lognormal_distribution_class
      implicit none

      type(mass_trans_cl)                         :: this
      logical                                     :: flag
      flag = .TRUE.

      select case (this % type_mass_trans)
            case ('MULTIRATE')
                 flag = inquire_homogeneity_multirate_ (this%multiRate)
            case ('SPHERICAL_DIFFUSION')
                 flag = inquire_homogeneity_spherical_ (this%sphericalDiff)
            case ('LAYERED_DIFFUSION')
                 flag = inquire_homogeneity_layered_ (this%layeredDiff)
            case ('CYLINDRICAL_DIFFUSION')
                 flag = inquire_homogeneity_cylindrical_ (this%cylindricalDiff)
            case ('POWER_LAW')
                 flag = inquire_homogeneity_power_ (this%power)
            case ('LOGNORMAL_LAW')
                 flag = inquire_homogeneity_lognormal_ (this%logNormal)
            case ('COMPOSITE_MEDIA')
                 flag = inquire_homogeneity_composite_ (this%composite)
      end select

  end function

  !************************************
  end module