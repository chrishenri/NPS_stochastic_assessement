!********************************************************************************************************************
!   Kinetic Reaction Class
!
!   Dave's approach
!
!********************************************************************************************************************
  module kinetic_reactions_class   ! class defining equilibrium reaction (it is only implemented linear sorption)
  use array_class
  use kdtree2_precision_module
  use kdtree2_module

  implicit none

  private
  public :: kinetic_reactions_cl                    !classes
  public ::                                 &       !methods
            read_kinetic_                ,  &
            delete_kinetic_              ,  &
            alloc_kinetic_               ,  &
            list_kinetic_                ,  &
            inquire_homogeneity_kinetic_ ,  &
            read_one_kinetic_reaction_            

  type rxn_cl
      integer                        :: nreactants           !total number of the species
      integer                        :: nproducts            !total number of the products
      integer, pointer, dimension(:) :: reactants => null()  !number of the species being the reactants of the reaction 
      integer, pointer, dimension(:) :: products  => null()  !number of the species being the products of the reaction
  end type

  type kinetic_reactions_cl 
      logical                     :: action 
      integer,            pointer :: nreact  => null()                  !number of chemical reactions
      type(rxn_cl),       pointer, dimension(:) :: rxn => null()        !chemical reaction
	  type(array_cl),     pointer, dimension(:) :: Kf => null()         !reaction rates coefficient (one for each kinetic reaction)
  end type

  contains

  !**************************** methods
  !************************************
  function read_kinetic_ (unit,geo) result(this)
      use array_class
	  use geometry_class
	  use global_variables
	  use to_solve
	  implicit none

      integer,           intent(in)     :: unit
	  type(geometry_cl), intent(in)     :: geo
	  character(len=200)                :: file
	  real*8                            :: const
	  integer                           :: ivar,flag,ispe,izoneim,nreact,irxn 
	  type(kinetic_reactions_cl)        :: this

      read(unit,*) this%action
      kineticACTION = this%action

	  read(unit,*) nreact

      select case (this%action)
	  case (.TRUE.)
	     call alloc_kinetic_ (this,nreact)
         do irxn=1,nreact
            call read_one_kinetic_reaction_ (this,irxn,unit)
         end do
         do irxn=1,nreact
            read(unit,*) file,const,ivar,flag      
	        this%Kf(irxn)  = read_array_ (file,ivar,const,flag,geo%nx,geo%ny,geo%nz)
	     end do
      case (.FALSE.)
         do irxn=1,nreact
            read(unit,*) !reactions
         end do
         do irxn=1,nreact
            read(unit,*) !read kf     
	     end do    
	  end select
	  
  end function
  
  
  !************************************
  subroutine read_one_kinetic_reaction_ (this,irxn,iunit)
     use gslib, only: generate_unit,count_words,locate_blanks,upper_case
     use global_variables, only: nspecie,namespecie
     implicit none

	 type(kinetic_reactions_cl), intent(inout) :: this
	 integer,                    intent(in)    :: irxn
     integer,                    intent(in)    :: iunit
     character (len=20), allocatable           :: words(:)
     logical                                   :: exists,match
     character (len=100)                       :: string
     integer                                   :: iword,nwords,nreactants,nproducts,ireactant,iproduct
     integer                                   :: locarrow,locend,locstart,isp
     integer, allocatable                      :: locblanks(:)     

     inquire (unit=iunit,exist=exists)
	 if(.not.exists) then
	     print *, 'file does not exist during read_assign_plane'
		 stop 'NOT a normal termination' 
	 end if
     read(iunit,'(a100)',err=10,end=11) string
     nwords = count_words (string)
     allocate(words(nwords))
     allocate(locblanks(1+nwords))
	 locblanks = locate_blanks(string)	  
	 do iword=1,nwords
		  words(iword)=string(locblanks(iword):locblanks(iword+1))
		  words(iword)=upper_case(words(iword))
		  if (trim(adjustl(words(iword)))== '-->') locarrow = iword
		  if (trim(adjustl(words(iword)))== '[')   locstart = iword
		  if (trim(adjustl(words(iword)))== ']')   locend   = iword
	 end do
     nreactants = (locarrow - locstart)/2
     nproducts  = (locend -locarrow)/2
     this%rxn(irxn)%nreactants = nreactants
     this%rxn(irxn)%nproducts  = nproducts
     if (nreactants > 0) allocate (this%rxn(irxn)%reactants(nreactants))
     if (nproducts  > 0) allocate (this%rxn(irxn)%products (nproducts ))
     if (nreactants > 0) then
        ireactant = 0
        do iword=locstart+1,locarrow,2
           ireactant = ireactant + 1
           match = .FALSE.
           do isp=1,nspecie
              if (trim(adjustl(words(iword)))== trim(adjustl(namespecie(isp)))) then
                  this%rxn(irxn)%reactants(ireactant) = isp
                  match = .TRUE.
              end if
           end do
           if (.not.match) stop '>> COULD NOT MATCH REACTANTS WITH SPECIES'
        end do
     end if
     if (nproducts>0) then
        iproduct  = 0
        do iword=locarrow+1,locend,2
          iproduct = iproduct + 1
          match = .FALSE.
          do isp=1,nspecie
              if (trim(adjustl(words(iword)))== trim(adjustl(namespecie(isp)))) then
                 this%rxn(irxn)%products(iproduct) = isp
                 match = .TRUE.
              end if       
          end do
          if (.not.match) write(*,*) '>> COULD NOT MATCH PRODUCTS WITH SPECIES'       
        end do
     end if
     return
     10   write(*,*) '>> error reading input file'
     11   write(*,*) '>> end-of-file reading input file'
     stop 'NOT a normal termination'   
        
  end subroutine


  !************************************
  subroutine delete_kinetic_ (this) ! destructor
      use array_class
	  implicit none	 
      type(kinetic_reactions_cl), intent(inout) :: this
      if (associated(this%Kf))     deallocate (this%Kf)
	  !this%action = .FALSE.
	  !call delete_array_ (this%Kf)
  end subroutine


  !************************************
  subroutine list_kinetic_ (this,fname)   !printer
  	  use gslib, only: open_fname
	  use global_variables, only: fdbg
	  use array_class
	  implicit none

	  character(len=*), intent(in)                  :: fname
	  type (kinetic_reactions_cl), intent (inout)   :: this
	  integer                                       :: unit, irxn

      call open_fname (fdbg,unit)
      write(unit,*)
      write(unit,*) ' KINETIC REACTION METHOD'
      write(unit,*)
      close (unit)

      if (associated(this%Kf)) then
        do irxn=1,this%nreact
            call list_array_ (this%Kf(irxn),fname,'Kf coeff ................:'); enddo; endif

  end subroutine


  !************************************
  subroutine alloc_kinetic_ (this,nreact)
     use global_variables, only: nspecie
     implicit none

	 type(kinetic_reactions_cl), intent(inout) :: this
     integer,                    intent(in)    :: nreact	 
         if (.not. associated (this % nreact)) then
               allocate (this % nreact)
               this%nreact = nreact
         end if
         if (.not. associated (this % rxn))        allocate (this % rxn(nreact))
         if (.not. associated (this % Kf))         allocate (this % Kf(nreact))
  end subroutine



  !************************************
  function inquire_homogeneity_kinetic_ (this) result (flag)
      use array_class
      use global_variables,    only: nspecie
      use heterogeneity_flags, only: kinetic_homogeneous
	  implicit none
	  type(kinetic_reactions_cl), intent(in):: this
	  logical                               :: flag
	  integer                               :: irx
	        flag = .TRUE.
            if(associated(this % nreact)) then
	        do irx=1,this%nreact
                if(associated(this % Kf)   .and. length_array_ (this % Kf(irx)) > 1) flag = .FALSE.
            enddo
            end if
            kinetic_homogeneous = flag

   end function

  !************************************
  
  end module kinetic_reactions_class



