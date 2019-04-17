!*************************************************************************************
!  GLOBAL VARIABLES: PROGRAM OPTIONS AND PACKAGES FLAGS 
!*************************************************************************************
  module global_variables
     use mf2k_class
     implicit none
     save
	 
	 type advectionDataFile_cl
        character(len=80) :: file
        real*8            :: const
        integer           :: ivar,flag
     end type

	 integer,            parameter      :: nfnam = 18, nmaxmove = 10000000, iwscreen = 1
	 !character(len=200), parameter      :: fdbg = 'rw3d.dbg'
     character(len=200)                 :: fdbg
     character(len=200), parameter      :: version='6.1',date='22 SEPT 2016'
 	 character(len=200)                 :: files_nam(nfnam)
	 !character(len=18),  parameter      :: fname_exit = 'rw3d_exit_part.dat' 
     character(len=200)                 :: fname_exit 
     character(len=200), allocatable    :: fname_well(:),fname_plane(:),fname_plume(:)
     real*8,             allocatable    :: tshot(:)
	 real*8                             :: courant,peclet,DaKINETIC,DaDECAY,DaMMT    !Grid courant, peclet number, Damkholer numbers
     real*8                             :: courant_wll=0.05,peclet_wll=0.05
     real*8                             :: DtStep            !Time step
	 real*8,             parameter      :: DaImin = 0.1      !Below this cut-off, particle moves as no mass transfer
	 real*8,             parameter      :: DaImax = 100.d0   !Above this cut-off, particle move as linear sorption
	 real*8,             parameter      :: UNEST  = -999.d0  !Value assigned to an unestimated variable
	 real*8                             :: TOL    = 1.0D-4   !Lower limit value of prob. to calculate spatial moments
     real*8,             parameter      :: EPS = 1.0d-15
	 character(len=200)                 :: calcul_time_method
	 character(len=15)                  :: method_advection
	 integer                            :: seed(2)
	 data                                  seed /123456789, 987654321/
	 integer*4                          :: seedzone = 123456789
	 integer*4                          :: seedbimo = 123456789
	 integer                            :: pathfreq,pathpart !frequency of printing particle path, and particle number
	 integer                            :: PostProcess_plume, PostProcess_btc
	 integer,           parameter       :: nmaxref = 20 !number of reflections in one time step (in case of inactive cells)
	 integer                            :: ndim = 3
	 logical                            :: ActiveDim(3) = .FALSE.
	 logical                            :: AlwaysPrintShots
	 logical                            :: SpeciesDispersionDependent =.FALSE.
	 real*8                             :: StartTimeInjection,EndTimeInjection
	 type(advectionDataFile_cl)         :: dataqx,dataqy,dataqz
	 real*8                             :: tsim
	 integer                            :: MinNumPart = 10
	 integer                            :: npInj   = 0
	 integer                            :: npWell  = 0
	 integer                            :: npPlane = 0
	 integer                            :: npOut   = 0
	 integer                            :: npStuck = 0
	 integer                            :: nspecie = 1
	 integer                            :: nspe_aq = 1
	 integer                            :: nspe_min = 0
	 integer                            :: nspedecay = 0
	 integer                            :: nzone   = 1
	 integer                            :: nzoneim = 0
	 integer, allocatable               :: np(:,:)
     real*8                             :: mpCUTOFF = 1.d-12
     integer                            :: ntshot = 0
     character(len=20),  pointer, dimension(:) :: namespecie  => null()     !name of species
     integer,            pointer, dimension(:) :: phasespecie => null()     !mobility of reactants and products (0=>mobile, 1=> immobile)
 	 integer,            pointer, dimension(:) :: idspedecay  => null()
	 type(mf2k_cl)                      :: mf2k
     real*8, allocatable                :: StressTime(:)
     logical                            :: particle_in_well_cell = .FALSE.
  end module
  
!*************************************************************************************
!  MODFLOW VARIABLES
!*************************************************************************************

  module velocity_times
     implicit none
     save
     
	 type period_cl
	     real*8  :: PERLEN
	     integer :: NSTP
	     real*8  :: MULT
	 end type period_cl
	
	type velotime_cl
		 integer                  :: NPER = 1       !total number of stress periods in velocity fields
		 type(period_cl), pointer :: PERIOD(:)      !stress period
         logical                  :: loop_per       !loop stress periods until tsim
         integer                  :: iloop_per      !loop index
         real*8                   :: time_flow
         logical                  :: restart_flux_from_mf2k = .FALSE.     !True if new loop index
    end type velotime_cl
    
    type(velotime_cl)    :: velotime
     
  end module


!*************************************************************************************
!  PROGRAM LOOPS COUNTS
!*************************************************************************************
  module loops_particles
     implicit none
     save    
	 integer :: ip,kinj,nmove
  end module


!*************************************************************************************
!  PROGRAM OPTIONS
!*************************************************************************************
  module code_options
     implicit none
     save    
	 integer :: idebug
	 integer :: ixmom,irmom,itmom,iwbtc,iwcbtc,iwDbtc,iwcshot,idilut,iwpath,ibcbtc,conv
     integer :: ixmompl,iwcshotpl,ipldisp
	 integer :: ipReStart = 1  ! particle from which restart simulation when Save Memory Mode
	 logical :: irestime
	 integer, parameter :: iprofile = 0
	 logical, parameter :: postprocess = .FALSE. 
     logical :: SaveMemo = .TRUE.
	 logical :: OC_plume = .FALSE.
	 logical :: OC_btc   = .FALSE.
     !integer :: simul_mode
  end module
  
  
!*************************************************************************************
!  CONSTANTS 
!*************************************************************************************
  module constants
     implicit none
	 save
     real*8, parameter :: pi = 3.14159265359
  end module 
  
  
!*************************************************************************************
!  PROBLEM DEFINITION
!*************************************************************************************
  module heterogeneity_flags
     implicit none
     save    
	 logical            :: vel_homogeneous         = .TRUE.   !darcy velocity
	 logical            :: poro_homogeneous        = .TRUE.   !porosity
	 logical            :: disp_homogeneous        = .TRUE.   !dispersivity
	 logical            :: mass_trans_homogeneous  = .TRUE.   !mass transfer
	 logical            :: react_homogeneous       = .TRUE.   !reaction 
	 logical            :: sorption_homogeneous    = .TRUE.   !sorption
	 logical            :: decay_homogeneous       = .TRUE.   !decay
	 logical            :: kinetic_homogeneous     = .TRUE.   !kinetic reaction
  end module
  
    !***************************
  module to_solve
     implicit none
     save    
     logical            :: mass_transACTION   = .FALSE.
     character(len=30)  :: mass_transTYPE     = ''
	 logical            :: sorptionACTION     = .FALSE.
	 character(len=30)  :: sorptionTYPE       = ''
	 logical            :: decayACTION        = .FALSE.
	 character(len=30)  :: decayTYPE          = ''
	 logical            :: kineticACTION      = .FALSE.
	 logical            :: recirculationACTION = .FALSE.
	 
  end module


!*************************************************************************************
!*************************************************************************************