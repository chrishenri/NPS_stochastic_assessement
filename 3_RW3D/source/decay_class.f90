!********************************************************************************************************************
!      Serial reaction
!********************************************************************************************************************
  module serial_class
  use array_class

  implicit none

  private
  public :: serial_cl                           !classes
  public ::                                  &  !methods
            read_serial_reaction_         ,  &
            delete_serial_reaction_       ,  &
			list_serial_reaction_         ,  &
			inquire_homogeneity_serial_

  type serial_cl
	  type(array_cl), pointer :: k(:)     => null() !reaction rate in mobile zone
	  type(array_cl), pointer :: kim(:,:) => null() !reaction rate in immobile zone
	  type(array_cl), pointer :: y(:,:)   => null() !yield coefficient
  end type

  contains

  !**************************** methods
  !************************************
  function read_serial_reaction_ (unit,geo) result (this) !reader
      use array_class
	  use geometry_class
	  use global_variables, only: nzoneim, nspedecay
	  use to_solve
	  implicit none

      integer,           intent(in)     :: unit
	  type(geometry_cl), intent(in)     :: geo
	  character(len=200)                :: file
	  real*8                            :: const
	  integer                           :: ivar,flag,ispe,izoneim
	  type(serial_cl)                   :: this


      allocate (this%k(nspedecay),this%y(nspedecay,nspedecay))
      do ispe=1,nspedecay
         read(unit,*) file,const,ivar,flag
         if (decayACTION) this%k(ispe) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
      end do

      do ispe=1,nspedecay
         read(unit,*) file,const,ivar,flag
         if (decayACTION .AND. ispe>1) then
            this%y(ispe,ispe-1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
         end if
      end do

      select case (mass_transTYPE)
      case ('MULTIRATE')
         allocate(this%kim(nspedecay,nzoneim))
         do ispe=1,nspedecay; do izoneim=1,nzoneim
            read(unit,*) file,const,ivar,flag
            if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
         end do; end do

      case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
         allocate(this%kim(nspedecay,1))
         do ispe=1,nspedecay
            read(unit,*) file,const,ivar,flag
            if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
         end do

      case ('POWER_LAW')
         allocate(this%kim(nspedecay,1))
         do ispe=1,nspedecay
            read(unit,*) file,const,ivar,flag
            if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
         end do

      case ('LOGNORMAL_LAW')
         allocate(this%kim(nspedecay,1))
         do ispe=1,nspedecay
            read(unit,*) file,const,ivar,flag
            if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
         end do

      case ('COMPOSITE_MEDIA')
         stop 'CANNOT SOLVE NETWORK REACTION AND COMPOSITE MEDIA (for the moment)'
            !allocate(this%kim(nspedecay,nzoneim))
            !do ispe=1,nspedecay; do izoneim=1,nzoneim
            !    this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
            !end do; end do

      end select

  end function

  !************************************
  subroutine delete_serial_reaction_ (this) !destructor
     use array_class
	 implicit none
	 type(serial_cl), intent(inout) :: this
		if (associated(this%k))     deallocate (this%k)
		if (associated(this%kim))   deallocate (this%kim)
		if (associated(this%y))     deallocate (this%y)
  end subroutine
  
  !************************************
  subroutine list_serial_reaction_ (this,fname)   !printer
  	  use gslib, only: open_fname
	  use global_variables, only: fdbg, nzoneim, nspedecay
	  use to_solve,         only: mass_transACTION, mass_transTYPE
	  use array_class
	  implicit none

	  character(len=*), intent(in)          :: fname
	  type (serial_cl), intent (inout)      :: this
	  integer                               :: unit, ispe, izoneim

                 call open_fname (fdbg,unit)
                 write(unit,*)
				 write(unit,*) ' SERIAL REACTION METHOD'
				 write(unit,*)
				 close (unit)

                 if (associated(this%k)) then
                    do ispe=1,nspedecay
                        call list_array_ (this%k(ispe),fname,'k coeff ................:'); enddo; endif

                 if (associated(this%y)) then
                    do ispe=2,nspedecay
                        call list_array_ (this%y(ispe,ispe-1),fname,'y coeff ................:'); enddo; endif

                 if (mass_transACTION .AND. nzoneim>0) then
                    select case (mass_transTYPE)
                    case ('MULTIRATE')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay; do izoneim=1,nzoneim
                                call list_array_ (this%kim(ispe,izoneim),fname,'kim coeff ................:'); enddo; enddo; endif
                    case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay
                                call list_array_ (this%kim(ispe,1),fname,'kim coeff ................:'); enddo; endif
                    case ('POWER_LAW')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay
                                call list_array_ (this%kim(ispe,1),fname,'kim coeff ................:'); enddo; endif
                    case ('LOGNORMAL_LAW')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay
                                call list_array_ (this%kim(ispe,1),fname,'kim coeff ................:'); enddo; endif
                    case ('COMPOSITE_MEDIA')
                        stop 'CANNOT SOLVE NETWORK REACTION AND COMPOSITE MEDIA'
                !allocate(this%kim(nspedecay,nzoneim))
                !do ispe=1,nspedecay; do izoneim=1,nzoneim
                !    this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                !end do; end do
                    end select
                 end if

  end subroutine


  !************************************
  function inquire_homogeneity_serial_ (this) result (flag)
      use array_class
      use global_variables, only: nzoneim, nspedecay
      use to_solve,         only: mass_transACTION, mass_transTYPE
	  implicit none
	  type(serial_cl), intent(in)        :: this
	  logical                            :: flag
	  integer                            :: i, ispe, izoneim
	        flag = .TRUE.
	        do ispe=1,nspedecay
                if(associated(this % k)   .and. length_array_ (this % k(ispe)) > 1) flag = .FALSE.;enddo

            do ispe=2,nspedecay
			    if(associated(this % y)   .and. length_array_ (this % y(ispe,ispe-1)) > 1) flag = .FALSE.;enddo
                
	        if (mass_transACTION .AND. nzoneim>0) then
	            select case (mass_transTYPE)
	            case ('MULTIRATE')
			        do ispe=1,nspedecay; do izoneim=1, nzoneim
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,izoneim)) > 1) flag = .FALSE.;enddo; enddo
	            case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
			        do ispe=1,nspedecay
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,1)) > 1) flag = .FALSE.;enddo
	            case ('POWER_LAW')
			        do ispe=1,nspedecay
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,1)) > 1) flag = .FALSE.;enddo
	            case ('LOGNORMAL_LAW')
			        do ispe=1,nspedecay
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,1)) > 1) flag = .FALSE.;enddo
	            case ('COMPOSITE_MEDIA')
			        stop 'CANNOT SOLVE NETWORK REACTION AND COMPOSITE MEDIA'
                !allocate(this%kim(nspedecay,nzoneim))
                !do ispe=1,nspedecay; do izoneim=1,nzoneim
                !    this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                !end do; end do
                    end select
	        end if
   end function

  !************************************
  end module serial_class


!********************************************************************************************************************
!      Serial reaction - motion given by moments
!********************************************************************************************************************
  module serial_mom_class
  use array_class

  implicit none

  private
  public :: serial_mom_cl                !classes
  public ::                             &       !methods
            read_serial_reac_mom_         ,  &
            delete_serial_reac_mom_       ,  &
			list_serial_reac_mom_         ,  &
			inquire_homogeneity_serial_mom_

  type serial_mom_cl
      integer, pointer :: nzones                    => null()
	  type(array_cl), pointer :: k(:)               => null() !reaction rate in mobile zone
	  type(array_cl), pointer :: y(:,:)             => null() !yield coefficient
	  real*8, pointer         :: S(:,:)             => null() !transformation matrix / eigenvectors
	  real*8, pointer         :: Sinv(:,:)          => null() !inverse of S
	  real*8, pointer         :: Reff(:,:)          => null() !effective retardations matrix
	  real*8, pointer         :: P(:,:)             => null() !state transition probability matrix
	  real*8, pointer         :: Areac(:,:,:)       => null() !center of mass position 3D matrix
	  real*8, pointer         :: B2(:,:,:)          => null() !spread 3D matrix
  end type

  contains

  !**************************** methods
  !************************************
  subroutine delete_serial_reac_mom_ (this) !destructor
     use array_class
	 implicit none
	 type(serial_mom_cl), intent(inout) :: this
		if (associated(this%k))         deallocate (this%k)
		if (associated(this%y))         deallocate (this%y)
		if (associated(this%S))         deallocate (this%S)
		if (associated(this%Sinv))      deallocate (this%Sinv)
		if (associated(this%Reff))      deallocate (this%Reff)
		if (associated(this%P))         deallocate (this%P)
		if (associated(this%Areac))     deallocate (this%Areac)
		if (associated(this%B2))        deallocate (this%B2)
  end subroutine

  !************************************
  function read_serial_reac_mom_ (unit,geo) result (this) !reader
      use array_class
	  use geometry_class
	  use global_variables, only: nzoneim, nspedecay
	  use to_solve
	  implicit none

      integer,           intent(in) :: unit
	  type(geometry_cl), intent(in) :: geo
	  character(len=200)            :: file
	  real*8                        :: const
	  integer                       :: ivar,flag,ispe
	  type(serial_mom_cl)           :: this

         if (decayACTION) then
            allocate(this%k(nspedecay),this%y(nspedecay,nspedecay),                  & 
                     this%S(nspedecay,nspedecay),this%Sinv(nspedecay,nspedecay),this%Reff(nspedecay,nspedecay),     &
                     this%P(nspedecay,nspedecay),this%Areac(nspedecay,nspedecay,3),this%B2(nspedecay,nspedecay,3))
         endif

         do ispe=1,nspedecay
           read(unit,*) file,const,ivar,flag
	       if (decayACTION)	this%k(ispe) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); end do

	     do ispe=1,nspedecay
           read(unit,*) file,const,ivar,flag
           if (decayACTION .AND. ispe>1)	this%y(ispe,ispe-1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz); end do

	     if (decayACTION .AND. mass_transACTION .AND. nzoneim>0) then
	        stop 'CANNOT SOLVE MASS TRANSFER AND SERIAL RX WITH HIGHER MOMENTS'
	     end if
  end function

  !************************************
  subroutine list_serial_reac_mom_ (this,fname)   !printer
  	  use gslib, only: open_fname
	  use global_variables, only: fdbg, nspedecay
	  use array_class
	  implicit none
	  character(len=*), intent(in)                   :: fname
	  type (serial_mom_cl), intent (inout)           :: this
	  integer                                        :: unit, ispe

                 call open_fname (fdbg,unit)
                 write(unit,*)
				 write(unit,*) ' SERIAL REACTION METHOD'
				 write(unit,*)
				 close(unit)

                 if (associated(this%k)) then
                    do ispe=1,nspedecay
                        call list_array_ (this%k(ispe),fname,'k coeff ................:'); enddo; endif
                 if (associated(this%y)) then
                    do ispe=2,nspedecay
                        call list_array_ (this%y(ispe,ispe-1),fname,'y coeff ................:'); enddo; endif

  end subroutine

  !************************************
  function inquire_homogeneity_serial_mom_ (this) result (flag)
      use array_class
      use global_variables, only : nspedecay
	  implicit none
	  type(serial_mom_cl), intent(in)    :: this
	  logical                            :: flag
	  integer                            :: i, ispe
	        flag = .TRUE.
	        do ispe=1,nspedecay
                if(associated(this % k)   .and. length_array_ (this % k(ispe)) > 1) flag = .FALSE.;enddo
            do ispe=2,nspedecay
			    if(associated(this % y)   .and. length_array_ (this % y(ispe,ispe-1)) > 1) flag = .FALSE.;enddo

   end function


  end module serial_mom_class


!********************************************************************************************************************
!      Generic network reaction
!********************************************************************************************************************
  module generic_net_class
  use array_class

  implicit none

  private
  public :: generic_net_cl                    !classes
  public ::                              &    !methods
            read_generic_net_         ,  &
            delete_generic_net_       ,  &
			list_generic_net_         ,  &
			inquire_homogeneity_generic_net_

  type generic_net_cl
	  type(array_cl), pointer :: k(:)     => null() !reaction rate
	  type(array_cl), pointer :: kim(:,:) => null() !reaction rate
	  type(array_cl), pointer :: y(:,:)   => null() !yield coefficient
  end type

  contains

  !**************************** methods
  !************************************
  function read_generic_net_ (unit,geo) result (this) !reader
      use array_class
	  use geometry_class
	  use global_variables, only: nspedecay, nzoneim
	  use to_solve !,         only: mass_transACTION, mass_transTYPE 
	  implicit none

      integer,           intent(in) :: unit  
	  type(geometry_cl), intent(in) :: geo
	  character(len=200)            :: file
	  real*8                        :: const
	  integer                       :: ivar,flag,ispe,jspe,izoneim
	  type(generic_net_cl)          :: this

      if (decayACTION) allocate (this%k(nspedecay),this%y(nspedecay,nspedecay))

      do ispe=1,nspedecay
            read(unit,*) file,const,ivar,flag
	        if (decayACTION) this%k(ispe) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
	  end do

      do ispe=1,nspedecay
      do jspe=1,nspedecay
            read(unit,*) file,const,ivar,flag
	        if (decayACTION) this%y(ispe,jspe) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
	  end do; end do

	  select case (mass_transTYPE)
	  case ('MULTIRATE')
            allocate(this%kim(nspedecay,nzoneim))
            do ispe=1,nspedecay; do izoneim=1,nzoneim
                read(unit,*) file,const,ivar,flag
                if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
            end do; end do
            
	  case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
            allocate(this%kim(nspedecay,1))
            do ispe=1,nspedecay
                read(unit,*) file,const,ivar,flag
                if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
            end do
            
	  case ('POWER_LAW')
            allocate(this%kim(nspedecay,1))
            do ispe=1,nspedecay
                read(unit,*) file,const,ivar,flag
                if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
            end do
            
	  case ('LOGNORMAL_LAW')
            allocate(this%kim(nspedecay,1))
            do ispe=1,nspedecay
                read(unit,*) file,const,ivar,flag
                if (decayACTION .AND. mass_transACTION .AND. nzoneim>0)  this%kim(ispe,1) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
            end do
            
	  case ('COMPOSITE_MEDIA')
            stop 'CANNOT SOLVE NETWORK REACTION AND COMPOSITE MEDIA'
                
	  end select

  end function

  !************************************
  subroutine delete_generic_net_ (this) !destructor
     use array_class
	 implicit none
	 type(generic_net_cl), intent(inout) :: this
		if (associated(this%k))     deallocate (this%k)
		if (associated(this%kim))   deallocate (this%kim)
		if (associated(this%y))     deallocate (this%y)
  end subroutine
  
  !************************************
  subroutine list_generic_net_ (this,fname)   !printer
  	  use gslib, only: open_fname
	  use global_variables, only: fdbg, nspedecay, nzoneim
	  use to_solve,         only: mass_transACTION, mass_transTYPE 
	  use array_class
	  implicit none
	  character(len=*), intent(in)                   :: fname
	  type (generic_net_cl), intent (inout)          :: this
	  integer                                        :: unit, ispe, jspe, izoneim

                 call open_fname (fdbg,unit)
                 write(unit,*)
				 write(unit,*) ' SERIAL REACTION METHOD'
				 write(unit,*)
				 close (unit)

                 if (associated(this%k)) then
                    do ispe=1,nspedecay
                        call list_array_ (this%k(ispe),fname,'k coeff ................:'); enddo; endif
                 if (associated(this%y)) then
                    do ispe=1,nspedecay; do jspe=1,nspedecay
                        call list_array_ (this%y(ispe,jspe),fname,'y coeff ................:'); enddo; enddo; endif

                 if (mass_transACTION .AND. nzoneim>0) then
                    select case (mass_transTYPE)
                    case ('MULTIRATE')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay; do izoneim=1,nzoneim
                                call list_array_ (this%kim(ispe,izoneim),fname,'kim coeff ................:'); enddo; enddo; endif
                    case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay
                                call list_array_ (this%kim(ispe,1),fname,'kim coeff ................:'); enddo; endif
                    case ('POWER_LAW')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay
                                call list_array_ (this%kim(ispe,1),fname,'kim coeff ................:'); enddo; endif
                    case ('LOGNORMAL_LAW')
                        if (associated(this%kim)) then
                            do ispe=1,nspedecay
                                call list_array_ (this%kim(ispe,1),fname,'kim coeff ................:'); enddo; endif
                    case ('COMPOSITE_MEDIA')
                        stop 'CANNOT SOLVE NETWORK REACTION AND COMPOSITE MEDIA'
                !allocate(this%kim(nspedecay,nzoneim))
                !do ispe=1,nspedecay; do izoneim=1,nzoneim
                !    this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                !end do; end do
                    end select
                 end if

  end subroutine

  !************************************
  function inquire_homogeneity_generic_net_ (this) result (flag)
      use array_class
      use global_variables, only: nspedecay, nzoneim
      use to_solve,         only: mass_transACTION, mass_transTYPE
	  implicit none
	  type(generic_net_cl), intent(in)   :: this
	  logical                            :: flag
	  integer                            :: i, ispe, jspe, izoneim
	        flag = .TRUE.
	        do ispe=1,nspedecay
                if(associated(this % k)   .and. length_array_ (this % k(ispe)) > 1) flag = .FALSE.;enddo
                
            do ispe=1,nspedecay; do jspe=1,nspedecay
			    if(associated(this % y)   .and. length_array_ (this % y(ispe,jspe)) > 1) flag = .FALSE.;enddo; enddo
			    
	        if (mass_transACTION .AND. nzoneim>0) then
	            select case (mass_transTYPE)
	            case ('MULTIRATE')
			        do ispe=1,nspedecay; do izoneim=1, nzoneim
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,izoneim)) > 1) flag = .FALSE.;enddo; enddo
	            case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION')
			        do ispe=1,nspedecay
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,1)) > 1) flag = .FALSE.;enddo
	            case ('POWER_LAW')
			        do ispe=1,nspedecay
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,1)) > 1) flag = .FALSE.;enddo
	            case ('LOGNORMAL_LAW')
			        do ispe=1,nspedecay
                        if(associated(this % kim) .and. length_array_ (this % kim(ispe,1)) > 1) flag = .FALSE.;enddo
	            case ('COMPOSITE_MEDIA')
			        stop 'CANNOT SOLVE NETWORK REACTION AND COMPOSITE MEDIA'
                !allocate(this%kim(nspedecay,nzoneim))
                !do ispe=1,nspedecay; do izoneim=1,nzoneim
                !    this%kim(ispe,izoneim) = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
                !end do; end do
                    end select
	        end if
   end function

  !************************************
  end module generic_net_class

!*********************************************************************************************************************
!   FIRST-ORDER DECAY CLASS
!*********************************************************************************************************************
  module decay_class
  use serial_class
  use serial_mom_class
  use generic_net_class
  implicit none

  private
  public :: decay_cl            !classes
  public ::                                    &       !methods
            read_decay_                    ,   &
			list_decay_                    ,   &
			delete_decay_                  ,   &
			inquire_homogeneity_decay_

       type decay_cl
	       logical                               :: action         !flag to distinguish if sorption is being used or not
	       character(len=30)                     :: type_decay
	       integer                               :: nspedecay = 0
	       character(len=20), pointer            :: namespedecay(:) => null()
	       integer, pointer                      :: idspedecay(:)   => null()
	       type(serial_cl),          pointer     :: serial          => null()
	       type(serial_mom_cl),      pointer     :: serial_mom      => null()
	       type(generic_net_cl),     pointer     :: generic_net     => null()
       end type
  
  contains

  !**************************** methods
  !************************************
  function read_decay_ (unit,geo) result (this) 
      use geometry_class
      use serial_class
      use serial_mom_class
      use global_variables, only: nspecie, namespecie, nspedecay, idspedecay
      use to_solve,         only: decayACTION, decayTYPE
      use gslib,            only: upper_case 
      implicit none

      integer,           intent(in) :: unit
      type(geometry_cl), intent(in) :: geo
      type(decay_cl)                :: this
	  integer                       :: i, j

      read(unit,*) this % action
      decayACTION   = this % action

      ! read species and keep track of species indexes
      read(10,*) this % nspedecay
      nspedecay = this % nspedecay
      allocate(this % namespedecay (nspedecay), this % idspedecay(nspedecay))
      read(10,*) (this % namespedecay(i), i=1, nspedecay)
	  do i=1, nspedecay
	  	 this % namespedecay(i)=upper_case(this % namespedecay(i))
	     do j=1,nspecie
	        if (this % namespedecay(i) == namespecie(j)) then
	            this % idspedecay(i) = j
	            exit
	        end if
	     end do
	  end do
	  allocate(idspedecay(nspedecay))
      idspedecay = this % idspedecay

      ! read type of decay network and associated parameters
      read(unit,*) this % type_decay
      this % type_decay = upper_case (this % type_decay)
      decayTYPE     = this % type_decay

      select case (this % type_decay)
      case ('SERIAL')
            allocate (this % serial)
            this % serial = read_serial_reaction_ (unit,geo)
      case ('SERIAL_MOMENTS')
            allocate (this % serial_mom)
            this % serial_mom = read_serial_reac_mom_ (unit,geo)
      case ('GENERIC')
            allocate (this % generic_net)
            this % generic_net = read_generic_net_ (unit,geo)
      case default
            if (this % action == .TRUE. ) then 
                print *, 'Decay network type not reconnized'
                stop
            end if
      end select

  end function

  !************************************
  subroutine list_decay_ (this,fname)
      use serial_class
      implicit none
	  character(len=*), intent(in)               :: fname      
	  type(decay_cl)                             :: this

	  if(associated(this%serial))         call list_serial_reaction_ (this%serial,fname)
	  if(associated(this%serial_mom))     call list_serial_reac_mom_ (this%serial_mom,fname)
	  if(associated(this%generic_net))    call list_generic_net_ (this%generic_net,fname)

  end subroutine

  !************************************
  subroutine delete_decay_ (this)
     use linear_sorption_class
	 implicit none
	 type(decay_cl)                               :: this

	 if(associated(this%serial))         call delete_serial_reaction_ (this%serial)
	 if(associated(this%serial_mom))     call delete_serial_reac_mom_ (this%serial_mom)
	 if(associated(this%generic_net))    call delete_generic_net_ (this%generic_net)

  end subroutine

  !************************************
  function inquire_homogeneity_decay_ (this)   result (flag)
     use linear_sorption_class
     use heterogeneity_flags
	 implicit none
     type(decay_cl)                               :: this
     logical                                      :: flag
	 
	 flag = .TRUE.
	 if(associated(this%serial))         flag = inquire_homogeneity_serial_ (this%serial)
	 if(associated(this%serial_mom))     flag = inquire_homogeneity_serial_mom_ (this%serial_mom)
	 if(associated(this%generic_net))    flag = inquire_homogeneity_generic_net_ (this%generic_net)
     decay_homogeneous = flag
  end function

  !************************************
  end module
