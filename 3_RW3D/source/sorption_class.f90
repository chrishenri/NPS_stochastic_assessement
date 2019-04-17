!********************************************************************************************************************
!      linear sorption model
!********************************************************************************************************************
  module linear_sorption_class   ! class defining equilibrium reaction (it is only implemented linear sorption)
  use array_class
  use global_variables, only: nspecie
  implicit none

  private
  public :: linear_sorption_rx_cl                !classes
  public ::                           &       !methods
            read_linear_sorption_,    &
            delete_linear_sorption_,  &
            list_linear_sorption_,    &
            inquire_homogeneity_linear_sorption_

  type linear_sorption_rx_cl
	  type(array_cl), pointer             :: R(:)       => null() !retardation factor
	  type(array_cl), pointer             :: Rim(:,:)   => null() !retardation factor
  end type

  contains

  !**************************** methods
  !************************************
  function read_linear_sorption_ (sorption,unit,geo) result (this) 
      use array_class
	  use geometry_class
      use global_variables, only: nspecie, nzoneim
      use to_solve,         only: sorptionACTION, mass_transACTION, mass_transTYPE
	  implicit none

      integer,           intent(in) :: unit
	  logical,           intent(in) :: sorption    
	  type(geometry_cl), intent(in) :: geo
	  
	  character(len=200)            :: file
	  real*8                        :: const
	  integer                       :: ivar,flag,ispe, izoneim
	  type(linear_sorption_rx_cl)      :: this

         allocate(this%R(nspecie))
         do ispe=1,nspecie
            read(unit,*) file,const,ivar,flag
            if (sorption)  this%R(ispe) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
         end do

         select case (mass_transTYPE)
         case ('MULTIRATE')
            allocate(this%Rim(nspecie,nzoneim))
            do ispe=1,nspecie; do izoneim=1,nzoneim
                read(unit,*) file,const,ivar,flag
                if (mass_transACTION .AND. nzoneim>0)  this%Rim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
            end do; end do

         case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
            allocate(this%Rim(nspecie,1))
            do ispe=1,nspecie
                read(unit,*) file,const,ivar,flag
                if (mass_transACTION .AND. nzoneim>0)  this%Rim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                end do

         case ('POWER_LAW')
            allocate(this%Rim(nspecie,1))
            do ispe=1,nspecie
                read(unit,*) file,const,ivar,flag
                if (mass_transACTION .AND. nzoneim>0)  this%Rim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                end do

         case ('LOGNORMAL_LAW')
            allocate(this%Rim(nspecie,1))
            do ispe=1,nspecie
                read(unit,*) file,const,ivar,flag
                if (mass_transACTION .AND. nzoneim>0)  this%Rim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                end do

         case ('COMPOSITE_MEDIA')
            stop 'CANNOT SOLVE SORPTION AND COMPOSITE MEDIA'
                !allocate(this%kim(nspecie,nzoneim))
                !do ispe=1,nspecie; do izoneim=1,nzoneim
                !    this%Rim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                !end do; end do
         end select

  end function

  !************************************
  subroutine delete_linear_sorption_ (this) !destructor
     use array_class
	 implicit none

	 type(linear_sorption_rx_cl), intent(inout) :: this

	 if (associated(this%R))    deallocate (this%R)
	 if (associated(this%Rim))  deallocate (this%Rim)

  end subroutine

  !************************************
  subroutine list_linear_sorption_ (this,fname)   !printing mass transfer parameters in a list
	  use gslib, only: open_fname
	  use global_variables, only: fdbg, nspecie, nzoneim
	  use to_solve,         only: mass_transACTION, mass_transTYPE
	  use array_class
	  implicit none

	  character(len=*), intent(in)                   :: fname
	  type (linear_sorption_rx_cl), intent (inout)      :: this
	  integer                                        :: unit, ispe, izoneim

      call open_fname (fdbg,unit)
      write(unit,*)
	  write(unit,*) ' LINEAR SORPTION MOBILE ZONE'
	  close (unit)

      if (associated(this%R)) then
         do ispe=1,nspecie
            call list_array_ (this%R(ispe),fname,'Retardation.............:'); end do; end if

      if (mass_transACTION .AND. nzoneim>0) then
         select case (mass_transTYPE)
         case ('MULTIRATE')
              if (associated(this%Rim)) then
                 do ispe=1,nspecie; do izoneim=1,nzoneim
                    call list_array_ (this%Rim(ispe,izoneim),fname,'Rim ................:'); enddo; enddo; endif
         case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
              if (associated(this%Rim)) then
                 do ispe=1,nspecie
                    call list_array_ (this%Rim(ispe,1),fname,'Rim ................:'); enddo; endif
         case ('POWER_LAW')
              if (associated(this%Rim)) then
                 do ispe=1,nspecie
                    call list_array_ (this%Rim(ispe,1),fname,'Rim ................:'); enddo; endif
         case ('LOGNORMAL_LAW')
              if (associated(this%Rim)) then
                 do ispe=1,nspecie
                    call list_array_ (this%Rim(ispe,1),fname,'Rim ................:'); enddo; endif
         case ('COMPOSITE_MEDIA')
            stop 'CANNOT SOLVE SORPTION AND COMPOSITE MEDIA'
         end select
      end if
  end subroutine

  !************************************
  function inquire_homogeneity_linear_sorption_ (this) result (flag)
      use array_class
      use global_variables, only: nspecie, nzoneim
      use to_solve,         only: mass_transACTION, mass_transTYPE
	  implicit none

	  type(linear_sorption_rx_cl), intent(in)   :: this
	  logical                                :: flag
	  integer                                :: ispe, izoneim
	  flag = .TRUE.

      do ispe=1,nspecie
         if(associated(this % R) .and. length_array_ (this % R(ispe)) > 1)  flag = .FALSE.; end do

      if (mass_transACTION .AND. nzoneim>0) then
	     select case (mass_transTYPE)
	     case ('MULTIRATE')
			  do ispe=1,nspecie; do izoneim=1, nzoneim
                 if(associated(this % Rim) .and. length_array_ (this % Rim(ispe,izoneim)) > 1) flag = .FALSE.;enddo; enddo
	     case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
			  do ispe=1,nspecie
                 if(associated(this % Rim) .and. length_array_ (this % Rim(ispe,1)) > 1) flag = .FALSE.;enddo
	     case ('POWER_LAW')
			  do ispe=1,nspecie
                 if(associated(this % Rim) .and. length_array_ (this % Rim(ispe,1)) > 1) flag = .FALSE.;enddo
	     case ('LOGNORMAL_LAW')
			  do ispe=1,nspecie
                 if(associated(this % Rim) .and. length_array_ (this % Rim(ispe,1)) > 1) flag = .FALSE.;enddo
	     case ('COMPOSITE_MEDIA')
			  stop 'CANNOT SOLVE LINEAR SORPTION AND COMPOSITE MEDIA (for the moment)'
         end select
      end if

   end function

  end module


!********************************************************************************************************************
!      CHTM
!********************************************************************************************************************
  module CHTM_class   ! class defining equilibrium reaction (it is only implemented linear sorption)
  use array_class
  use global_variables, only: nspecie
  implicit none

  private
  public :: CHTM_rx_cl                !classes
  public ::                           &       !methods
            read_CHTM_,                 &
            delete_CHTM_,               &
            list_CHTM_,                 &
            inquire_homogeneity_CHTM_

  type CHTM_rx_cl
	  type(array_cl), pointer             :: bd    => null() !bulk density
	  type(array_cl), pointer             :: kf    => null() !forward rate
	  type(array_cl), pointer             :: kb    => null() !backward rate
  end type

  contains

  !**************************** methods
  !************************************
  function read_CHTM_ (sorption,unit,geo) result (this) 
      use array_class
	  use geometry_class
      use global_variables, only: nspecie, nzoneim
      use to_solve,         only: mass_transACTION
	  implicit none

      integer,           intent(in) :: unit
	  logical,           intent(in) :: sorption    
	  type(geometry_cl), intent(in) :: geo

	  character(len=200)            :: file
	  real*8                        :: const
	  integer                       :: ivar,flag
	  type(CHTM_rx_cl)                 :: this

      read(unit,*) file,const,ivar,flag
      if (mass_transACTION .OR. nspecie>1) then
         stop 'CANNOT SOLVE CHTM AND MASS TRANSFER'
      end if

	  if (sorption)	 then 
         allocate (this%bd, this%kf, this%kb)
		 this%bd = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
		 this%kf = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
		 this%kb = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
	  end if

  end function

  !************************************
  subroutine delete_CHTM_ (this) !destructor
     use array_class
	 implicit none

	 type(CHTM_rx_cl), intent(inout) :: this

	 if (associated(this%kb))    deallocate (this%kb)
	 if (associated(this%kf))    deallocate (this%kf)
	 if (associated(this%kb))    deallocate (this%kb)

  end subroutine

  !************************************
  subroutine list_CHTM_ (this,fname)   !printing mass transfer parameters in a list
	  use gslib, only: open_fname
	  use global_variables, only: fdbg, nzoneim
	  use to_solve,         only: mass_transACTION
	  use array_class
	  implicit none

	  character(len=*), intent(in)                   :: fname
	  type (CHTM_rx_cl), intent (inout)                 :: this
	  integer                                        :: unit

      call open_fname (fdbg,unit)
      write(unit,*)
	  write(unit,*) ' LINEAR SORPTION MOBILE ZONE'
	  close (unit)

      if (associated(this%bd)) then
         call list_array_ (this%bd,fname,'Retardation.............:'); end if

      if (associated(this%kf)) then
         call list_array_ (this%kf,fname,'Retardation.............:'); end if

      if (associated(this%kb)) then
         call list_array_ (this%kb,fname,'Retardation.............:'); end if

      if (mass_transACTION .AND. nzoneim>0) then
          stop 'CANNOT SOLVE CHTM AND MASS TRANSFER'
      end if
  end subroutine

  !************************************
  function inquire_homogeneity_CHTM_ (this) result (flag)
      use array_class
      use global_variables, only: nzoneim
      use to_solve,         only: mass_transACTION
	  implicit none

	  type(CHTM_rx_cl), intent(in)              :: this
	  logical                                :: flag
	  flag = .TRUE.

      if(associated(this % bd) .and. length_array_ (this % bd) > 1)  flag = .FALSE.
      if(associated(this % kf) .and. length_array_ (this % kf) > 1)  flag = .FALSE.
      if(associated(this % kb) .and. length_array_ (this % kb) > 1)  flag = .FALSE.

      if (mass_transACTION .AND. nzoneim>0) then
          stop 'CANNOT SOLVE CHTM AND MASS TRANSFER'
      end if

   end function

  end module

!*********************************************************************************************************************
!   Sorption class
!*********************************************************************************************************************
  module sorption_class
  use linear_sorption_class
  use CHTM_class
  implicit none

  private
  public :: sorption_rx_cl   !classes
  public ::                                    &       !methods
  			read_sorption_                 ,   &
			list_sorption_                 ,   &
			delete_sorption_               ,   &
			inquire_homogeneity_sorption_

       type sorption_rx_cl
           character(len=30)                      :: type_sorption
	       logical                                :: action         !flag to distinguish if sorption is being used or not 
	       type(linear_sorption_rx_cl),  pointer     :: LinearSorption => null()
	       type(CHTM_rx_cl),  pointer                :: CHTM => null()
       end type
  
  contains

  !**************************** methods
  !************************************
  function read_sorption_ (unit,geo) result (this) 
      use geometry_class
      use linear_sorption_class
      use to_solve,         only: sorptionACTION, sorptionTYPE
      use gslib,            only: upper_case 
      implicit none

      integer,           intent(in)         :: unit
      type(geometry_cl), intent(in)         :: geo
      type(sorption_rx_cl)                  :: this

      read(unit,*) this % action
      sorptionACTION    = this % action
      
      read(unit,*) this % type_sorption
      this % type_sorption = upper_case (this % type_sorption)
      sorptionTYPE    = this % type_sorption
      
      select case (this % type_sorption)

      case ('LINEAR')
            allocate (this % LinearSorption)
            this % LinearSorption = read_linear_sorption_ (this % action,unit,geo)
      case ('CHTM')
            allocate (this % CHTM)
            this % CHTM = read_CHTM_ (this % action,unit,geo)

      case default
            if (this % action == .TRUE. ) then 
                stop 'Sorption type not reconnized'
            end if

      end select

  end function

  !************************************
  subroutine list_sorption_ (this,fname)
      use linear_sorption_class
      implicit none
	  character(len=*), intent(in)                :: fname      
	  type(sorption_rx_cl)                        :: this

	  if(associated(this%LinearSorption))         call list_linear_sorption_ (this%LinearSorption,fname)
      if(associated(this%CHTM))                   call list_CHTM_ (this%CHTM,fname)

  end subroutine

  !************************************
  subroutine delete_sorption_ (this)
     use linear_sorption_class
	 implicit none
	 type(sorption_rx_cl)                         :: this

	 if(associated(this%LinearSorption))          call delete_linear_sorption_ (this%LinearSorption)
	 if(associated(this%CHTM))                    call delete_CHTM_ (this%CHTM)

  end subroutine

  !************************************
  function inquire_homogeneity_sorption_   (this) result (flag)
     use linear_sorption_class
     use heterogeneity_flags
	 implicit none
     type(sorption_rx_cl)                        :: this
     logical                                     :: flag
	 flag = .TRUE.

	 if(associated(this%LinearSorption))   flag = inquire_homogeneity_linear_sorption_ (this%LinearSorption)
	 if(associated(this%CHTM))             flag = inquire_homogeneity_CHTM_ (this%CHTM)
     sorption_homogeneous = flag
  end function

  !************************************
  end module