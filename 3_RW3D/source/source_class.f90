

!******************************************************************************
!  SOURCE CLASS 
!******************************************************************************
 module source_class
   use timefunction_class
   implicit none

   private
   public :: source_cl             !class
   public ::                             & ! methods
                print_source_          , &
				delete_source_         , &
				read_source_           , &
                find_loc_TimeFunction_      
				  
   type parameters_cl
	   character(len=200), pointer :: file => null() 
	   real*8,  pointer :: xinj   => null(), yinj   => null(), zinj   => null()
	   real*8,  pointer :: zbot   => null(), ztop   => null(), rcyr   => null()
	   real*8,  pointer :: rcp    => null(), xdist  => null(), width  => null(), height => null()
       real*8,  pointer :: xinj_1 => null(), yinj_1 => null(), zinj_1 => null(), xinj_2 => null()
	   real*8,  pointer :: yinj_2 => null(), zinj_2 => null()
	   integer, pointer :: idwn   => null(), jdwn   => null(), kdwn   => null()
	   integer, pointer :: iup    => null(), jup    => null(), kup    => null()
	   real*8,  pointer :: np11x  => null(), np11y  => null(), np11z  => null()
	   real*8,  pointer :: const
	   integer, pointer :: ix(:)  => null(), iy(:)  => null(), iz(:)  => null()
       integer, pointer :: npblock(:) => null()
       real*8,  pointer :: mpblock(:) => null()
       character(len=100), pointer :: dist_part(:) => null()
   end type 
   
   
   type source_cl
	   character(len=200)     :: name 
       character(len=10)      :: TypeInj
       real*8                 :: TimeStartInj
       real*8                 :: TimeStopInj
       integer                :: np 
	   real*8                 :: pmass
       real*8                 :: pmassINI 
	   integer                :: zone
	   integer                :: specie
       type(parameters_cl)    :: par
       type(timefunction_cl)  :: timefunct
       integer                :: loc   !interval of time in timefunction in which the simulation is
       integer                :: freq
   end type
   

   interface delete_source_
       module procedure source_null; end interface

   
   contains



   subroutine find_loc_TimeFunction_ (source,time)
       implicit none
       type(source_cl),     intent(inout) :: source
       real*8,              intent(in)    :: time
       integer                            :: i       
       do i=source%loc,source%timefunct%nt-1
         if (source%timefunct%time(i)<=time.and.source%timefunct%time(i+1)>time) then
            source%loc = i 
            return
         end if
       end do
       source%loc = source%timefunct%nt
    end subroutine


!   --------------------------------------------------------------------------------
!   destructor
!	--------------------------------------------------------------------------------
	subroutine source_null (this)
    implicit none
	
	type(source_cl), intent(inout) :: this
	
	   this % name      = ' '
	   this % np        = 0
	   this % pmass     = 0.d0
       this % pmassINI  = 0.d0
	   this % zone      = 0
	   this % specie    = 0
	   this % freq      = 1
	   this % TimeStartInj = 0.d0
	   this % TimeStopInj  = 0.d0
       this % loc       = 1
	   
	   nullify(this % par % xinj)
	   nullify(this % par % yinj)
	   nullify(this % par % zinj)
	   nullify(this % par % zbot)
	   nullify(this % par % ztop)
	   nullify(this % par % rcyr)
	   nullify(this % par % rcp )
	   nullify(this % par % xdist)
	   nullify(this % par % width )
	   nullify(this % par % height)
	   nullify(this % par % xinj_1)
	   nullify(this % par % yinj_1)
	   nullify(this % par % zinj_1)
	   nullify(this % par % xinj_2)
	   nullify(this % par % yinj_2)
	   nullify(this % par % zinj_2)
	   nullify(this % par % idwn  )
	   nullify(this % par % jdwn  )
	   nullify(this % par % kdwn  )
	   nullify(this % par % iup   )
	   nullify(this % par % jup   )
	   nullify(this % par % kup   )
	   nullify(this % par % np11x )
	   nullify(this % par % np11y )
	   nullify(this % par % np11z )
	   nullify(this % par % file  )
	   nullify(this % par % ix    )
	   nullify(this % par % iy    )
	   nullify(this % par % iz    )	
       
    end subroutine

!   --------------------------------------------------------------------------------------
!   reader
!	-------------------------------------------------------------------------------------- 			     

    subroutine read_source_ (this,namein,TypeInjIn,fname)
	   use gslib, only: upper_case, open_fname, open_fname_normal
	   use global_variables, only: UNEST
	   implicit none
	   type(source_cl),          intent(inout) :: this
	   character(len=*),         intent(in) :: namein,TypeInjIn
	   character(len=*),         intent(in) :: fname
	   character(len=len(trim(adjustl(namein))))    :: name
	   character(len=len(trim(adjustl(TypeInjIn)))) :: name2
	   character(len=100)                   :: string
	   real*8                               :: xinj,yinj,zinj,zbot,ztop,rcyr,rcp,xdist,width,height
	   real*8                               :: xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2,pmass,totmass
	   integer                              :: zone,specie
       integer                              :: idwn,jdwn,kdwn,iup,jup,kup,np
	   real*8                               :: np11x,np11y,np11z
       integer                              :: unit,unit2
	   logical                              :: existeix
       character(len=100)                   :: file
       real*8                               :: const
       integer                              :: flag,freq
       logical                              :: ReadFromFile
       integer                              :: it
       
       integer,allocatable                  :: ix(:),iy(:),iz(:),npblock(:)
       real*8, allocatable                  :: mpblock(:)
       character(len=100)                   :: dist(3)

           call source_null(this)

           call open_fname (fname,unit)

           name  = trim(adjustl(namein))
           name2 = trim(adjustl(TypeInjIn))

		   name  = upper_case (name)
		   name2 = upper_case (name2)

		   this % name = name
		   this % TypeInj = name2	   
		   
		   np    = UNEST
		   
           ReadFromFile = .FALSE.
		   
		   select case (trim(adjustl(name)))

		    case ('BLOCK' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,iup,jup,kup !,np11x,np11y,np11z 
	               allocate(this % par % idwn); this % par % idwn  = idwn
	               allocate(this % par % jdwn); this % par % jdwn  = jdwn
	               allocate(this % par % kdwn); this % par % kdwn  = kdwn
                   allocate(this % par % iup);  this % par % iup   = iup
                   allocate(this % par % jup);  this % par % jup   = jup
                   allocate(this % par % kup);  this % par % kup   = kup
                   !allocate(this % np11x);  this % np11x = np11x
                   !allocate(this % np11y);  this % np11y = np11y
                   !allocate(this % np11z);  this % np11z = np11z
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
                   
		    case ('BLOCK_FLUX_WEIGHTED' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,iup,jup,kup !,np11x,np11y,np11z 
	               allocate(this % par % idwn); this % par % idwn  = idwn
	               allocate(this % par % jdwn); this % par % jdwn  = jdwn
	               allocate(this % par % kdwn); this % par % kdwn  = kdwn
                   allocate(this % par % iup);  this % par % iup   = iup
                   allocate(this % par % jup);  this % par % jup   = jup
                   allocate(this % par % kup);  this % par % kup   = kup
                   !allocate(this % np11x);  this % np11x = np11x
                   !allocate(this % np11y);  this % np11y = np11y
                   !allocate(this % np11z);  this % np11z = np11z
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
 
		   case ( 'POINT' )        
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zinj 
			  	   allocate(this % par % xinj); this % par % xinj = xinj
	               allocate(this % par % yinj); this % par % yinj = yinj
	               allocate(this % par % zinj); this % par % zinj = zinj
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'LINE' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zbot,ztop 
	               allocate(this % par % xinj); this % par % xinj = xinj  
	               allocate(this % par % yinj); this % par % yinj = yinj
	               allocate(this % par % zbot); this % par % zbot = zinj 
                   allocate(this % par % ztop); this % par % ztop = ztop  
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
		             
		   case ( 'CIRCLE' )        
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zbot,ztop,rcyr
	               allocate(this % par % xinj); this % par % xinj = xinj
	               allocate(this % par % yinj); this % par % yinj = yinj
	               allocate(this % par % zbot); this % par % zbot = zbot
                   allocate(this % par % ztop); this % par % ztop = ztop
                   allocate(this % par % rcyr); this % par % rcyr = rcyr
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
	       
		   case ( 'RADIAL' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj,yinj,zbot,ztop,rcp 
  	               allocate (this % par % xinj);this % par % xinj = xinj  
	               allocate (this % par % yinj);this % par % yinj = yinj 
	               allocate (this % par % zbot);this % par % zbot = zbot
                   allocate (this % par % ztop);this % par % ztop = ztop
                   allocate (this % par % rcp); this % par % rcp  = rcp
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
	       
		   case ( 'PLANE_RANDOM' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xdist,width,height
	               allocate (this % par % xdist) ; this % par % xdist  = xdist 
	               allocate (this % par % width) ; this % par % width  = width
	               allocate (this % par % height); this % par % height = height
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
	       
		   case ( 'PLANE'  )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xdist,width,height
	               allocate (this % par % xdist) ; this % par % xdist  = xdist 
	               allocate (this % par % width) ; this % par % width  = width
	               allocate (this % par % height); this % par % height = height
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
		   case ( 'PLANE_FLUX_WEIGHTED'  )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xdist,width,height
	               allocate (this % par % xdist) ; this % par % xdist  = xdist 
	               allocate (this % par % width) ; this % par % width  = width
	               allocate (this % par % height); this % par % height = height
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
		   
		   case ( 'LINE_BY_POINTS' )			  
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1
		           allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1  
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_2); this % par % zinj_2 = zinj_2
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
		   case ( 'LINE_FLUX_WEIGHTED' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1
		           allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1  
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_2); this % par % zinj_2 = zinj_2
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
		   case ( 'LINE_BY_POINTS_RANDOM' )			  
			  read(unit,*) pmass,zone,specie
			  read(unit,*) xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1
		           allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1  
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_2); this % par % zinj_2 = zinj_2
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'VERTICAL_LINE_FLUX_WEIGHTED' )	          
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,kup  
	               allocate(this % par % idwn);  this % par % idwn  = idwn
	               allocate(this % par % jdwn);  this % par % jdwn  = jdwn
	               allocate(this % par % kdwn);  this % par % kdwn  = kdwn
                   allocate(this % par % kup);   this % par % kup   = kup
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
		   case ( 'READ_PARTICLE_FILE' )
		       if (this%TypeInj /= 'DIRAC') then
		           stop '...reading from concentration file only accessible for DIRAC injections'
		       end if
		       ReadFromFile = .TRUE.
		       read(unit,*) file
			   inquire (file=file,exist=existeix)
			   if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop; end if
	           call open_fname_normal ( file, unit2 )
	           read(unit2,*)
			   read(unit2,*) np
		       allocate (this % par % file)
		       this % par % file = file
			   this % np     = np
			   this % pmass  =  UNEST
	           this % zone   = UNEST
	           this % specie = UNEST			   
               close(unit2)

		   case ( 'READ_CONCENTRATION_FILE' )
		       if (this%TypeInj /= 'DIRAC') then
		           stop '...reading from concentration file only accessible for DIRAC injections'
		       end if
		       ReadFromFile = .TRUE.
               read(unit,*) pmass,zone,specie
		       read(unit,*) file,const
			   inquire (file=file,exist=existeix)
			   if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop; end if
		       allocate (this % par % file)
		       allocate (this % par % const)
		       this % par % file = file
		       this % par % const = const
			   this % np     = UNEST
			   this % pmass  = pmass
	           this % zone   = zone
	           this % specie = specie

		   case ( 'VERTICAL_BLOCK_FLUX_WEIGHTED' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) idwn,jdwn,kdwn,kup  
	               allocate(this % par % idwn);  this % par % idwn  = idwn
	               allocate(this % par % jdwn);  this % par % jdwn  = jdwn
	               allocate(this % par % kdwn);  this % par % kdwn  = kdwn
                   allocate(this % par % kup);   this % par % kup   = kup
				   this % np    = np
				   this % pmass = pmass
	               this % zone   = zone
	               this % specie = specie

		   case ( 'CELLS_FILE' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) file
			  inquire (file=file,exist=existeix)
			  if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop
              end if
              allocate (this % par % file)
              this % par % file = file
			  call read_blocks_ (file,ix,iy,iz)
			       allocate(this % par % ix(size(ix))); this % par % ix = ix
			       allocate(this % par % iy(size(iy))); this % par % iy = iy
			       allocate(this % par % iz(size(iz))); this % par % iz = iz
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
		   case ( 'CELLS_FILE_FLUX_WEIGHTED' )
			  read(unit,*) pmass,zone,specie
			  read(unit,*) file
			  inquire (file=file,exist=existeix)
			  if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop
              end if
              allocate (this % par % file)
              this % par % file = file
			  call read_blocks_ (file,ix,iy,iz)
			       allocate(this % par % ix(size(ix))); this % par % ix = ix
			       allocate(this % par % iy(size(iy))); this % par % iy = iy
			       allocate(this % par % iz(size(iz))); this % par % iz = iz
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
           case ( 'CELLS_FILE_PARTICLE_NUMBER' )
              if (this%TypeInj /= 'DIRAC') then
		          stop '...reading from cells file (with particle number per cell) only accessible for DIRAC injections'
		      end if
		      ReadFromFile = .TRUE.
			  read(unit,*) zone,specie
			  read(unit,*) file
			  inquire (file=file,exist=existeix)
			  if (.not. existeix) then
			            print *, '**** file does not exist: ',file
                        stop
              end if
              allocate (this % par % file)
              this % par % file = file
			  call read_blocks_part_ (file,ix,iy,iz,np,pmass,dist,npblock,mpblock)
			       allocate(this % par % ix(size(ix))); this % par % ix = ix
			       allocate(this % par % iy(size(iy))); this % par % iy = iy
			       allocate(this % par % iz(size(iz))); this % par % iz = iz
                   allocate(this % par % npblock(size(iz))); this % par % npblock = npblock
                   allocate(this % par % mpblock(size(iz))); this % par % mpblock = mpblock
                   allocate(this % par % dist_part(3)); this % par % dist_part = dist
				   this % np     = np
				   !this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
           case ( 'HORIZONTAL_PLANE_FLUX_WEIGHTED' )
               read(unit,*) pmass,zone,specie
               read(unit,*) xinj_1,yinj_1,xinj_2,yinj_2,zinj_1
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1 
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
                   
           case ( 'WATER_TABLE_FLUX_WEIGHTED' )
               read(unit,*) pmass,zone,specie
               read(unit,*) xinj_1,yinj_1,xinj_2,yinj_2,zinj_1 !here zinj_1 is the depth bellow the water table
	  	           allocate (this % par % xinj_1); this % par % xinj_1 = xinj_1   
	               allocate (this % par % yinj_1); this % par % yinj_1 = yinj_1 
                   allocate (this % par % xinj_2); this % par % xinj_2 = xinj_2 
                   allocate (this % par % yinj_2); this % par % yinj_2 = yinj_2  
                   allocate (this % par % zinj_1); this % par % zinj_1 = zinj_1
				   this % np     = np
				   this % pmass  = pmass
	               this % zone   = zone
	               this % specie = specie
                   
	       case default
	       
		      stop '**** could not match name of injection ****'

           end select
           
           !....assign initial mass

           this%pmassINI = this%pmass

           !....read time function
           
           if (this%TypeInj == 'DIRAC') then             
              read(unit,*) this%TimeStartInj
              if(.not.ReadFromFile) read(unit,*) this%np
              this%TimeStopInj = this%TimeStartInj 
           end if

           if (this%typeInj == 'GENERAL') then
              read(unit,*) file,const
              call read_timefunction_ (file,const,1,this%timefunct)
              !read(unit,*) freq
              this % freq = 1
              do it=1,this%timefunct%nt
                if (this%timefunct%val(it) /= 0.d0) then
                   this % TimeStartInj = this%timefunct%time(it)
                   exit
                end if
              end do
              do it=this%timefunct%nt,1,-1
                if (this%timefunct%val(it) /= 0.d0) then
                   if (it==this%timefunct%nt) then
                      this % TimeStopInj = this%timefunct%time(this%timefunct%nt)
                   elseif (it<this%timefunct%nt) then 
                      this % TimeStopInj = this%timefunct%time(it+1)
                   end if
                   exit
                end if
              end do
           end if
           

	end subroutine
!   --------------------------------------------------------------------------------------
!   printer
!	-------------------------------------------------------------------------------------- 			     

    subroutine read_blocks_ (fname,ix,iy,iz)
        use gslib, only: open_fname_normal
        implicit none
        character(len=*),      intent(in)    :: fname
        integer, allocatable,  intent(inout) :: ix(:),iy(:),iz(:)
        integer                              :: nblock,ivar,i,iunit

        call open_fname_normal(fname,iunit)
        read(iunit,*) nblock
        allocate(ix(nblock),iy(nblock),iz(nblock))
        read(iunit,*) ivar
        do i=1,ivar
            read(iunit,*)
        end do
        do i=1,nblock
            read(iunit,*) ix(i),iy(i),iz(i)
        end do
        close(iunit)
    end subroutine
    
    
    subroutine read_blocks_part_ (fname,ix,iy,iz,np,pmass,dist,npblock,mpblock)
        use gslib, only: open_fname_normal,upper_case
        implicit none
        character(len=*),      intent(in)    :: fname
        integer,               intent(inout) :: np
        real*8,                intent(inout) :: pmass
        character(len=*),      intent(inout) :: dist(3)
        integer, allocatable,  intent(inout) :: ix(:),iy(:),iz(:),npblock(:)
        real*8,  allocatable,  intent(inout) :: mpblock(:)
        integer                              :: nblock,ivar,i,iunit

        call open_fname_normal(fname,iunit)
        read(iunit,*) nblock
        read(iunit,*) np
        !read(iunit,*) pmass
        do i=1,3
            read(iunit,*) dist(i)
            dist(i) = upper_case(dist(i))
        end do
        allocate(ix(nblock),iy(nblock),iz(nblock),npblock(nblock),mpblock(nblock))
        do i=1,nblock
            read(iunit,*) ix(i),iy(i),iz(i),npblock(i),mpblock(i)
        end do
        close(iunit)
    end subroutine
    


    subroutine print_source_ (this,fname)
	   use gslib, only: upper_case, open_fname
	   implicit none
	   type(source_cl),       intent(inout) :: this
	   character(len=*),      intent(in)    :: fname
	   character(len=len(trim(adjustl(this%name)))) :: name
	   real*8                               :: xinj,yinj,zinj,zbot,ztop,rcyr,rcp,xdist,width,height
	   real*8                               :: xinj_1,yinj_1,zinj_1,xinj_2,yinj_2,zinj_2
	   integer                              :: zone,specie
       integer                              :: idwn,jdwn,kdwn,iup,jup,kup
       integer                              :: unit


           call open_fname (trim(adjustl(fname)),unit)

		   write(unit,'(a26,x,i5)')  'mobile-immobile Zone....: ',this%zone
		   write(unit,'(a26,x,i5)')  'species.................: ',this%specie


           name = trim(adjustl(this % name))
		   
		   if ( trim(adjustl(name)) == 'BLOCK' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,11)                                 &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'i-up..: ',  this % par%iup,   &
								   'j-up..: ',  this % par%jup,   &
								   'k-up..: ',  this % par%kup    
               11  format (/,x,6(a8,i5,x)/)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
                   
                   
		   elseif ( trim(adjustl(name)) == 'BLOCK_FLUX_WEIGHTED' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,11)                                 &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'i-up..: ',  this % par%iup,   &
								   'j-up..: ',  this % par%jup,   &
								   'k-up..: ',  this % par%kup    

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)    
 
		   elseif ( trim(adjustl(name)) == 'POINT' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '

	               write(unit,10)   &
				                  'xinj: ',this % par%xinj,  &
								  'yinj: ',this % par%yinj,  &
								  'zinj: ',this % par%zinj
               10  format (/,x,3(a6,g15.6,x)/) 

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 
		   elseif (trim(adjustl(name)) == 'LINE' ) then	          
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
	               write(unit,12)  &
				                  'xinj: ',this % par%xinj,  &
								  'yinj: ',this % par%yinj,  &
								  'zbot: ',this % par%zbot,  &
								  'ztop: ',this % par%ztop 
               12  format (/,x,4(a6,g15.6,x)/)		             

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 
		   elseif (trim(adjustl(name)) == 'CIRCLE' ) then        
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
	               write(unit,13)  &
				                  'xinj: ',this % par%xinj,  &
								  'yinj: ',this % par%yinj,  &
								  'zbot: ',this % par%zbot,  &
								  'ztop: ',this % par%ztop,  &
								  'rcyr: ',this % par%rcyr
               13  format (/,x,5(a6,g15.6,x)/)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 	       
		   elseif ( trim(adjustl(name)) == 'RADIAL' ) then	          
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
	               write(unit,14)	  &
				                    'xinj: ',this % par%xinj,  &
									'yinj: ',this % par%yinj,  &
									'zbot: ',this % par%zbot,  &
									'ztop: ',this % par%ztop,  &
									'rcp:  ',this % par%rcp
               14  format(/,x,5(a6,g15.6,x)/)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 	       
		   elseif ( trim(adjustl(name)) == 'PLANE_RANDOM' ) then	          
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
	               write(unit,15)	 &
				                    'xdist: ',  this % par%xdist, &
									'width: ',  this % par%width, &
									'height: ', this % par%height
               15  format(/,x,3(a7,g15.6,x)/)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 
		   elseif ( trim(adjustl(name)) == 'PLANE' ) then	          
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,16)	 &
				                    'xdist: ',  this % par%xdist, &
									'width: ',  this % par%width, &
									'height: ', this % par%height
               16  format(/,x,3(a7,g15.6,x)/)
                   
                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)

		   elseif ( trim(adjustl(name)) == 'PLANE_FLUX_WEIGHTED' ) then	          
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,19)	 &
				                    'xdist: ',  this % par%xdist, &
									'width: ',  this % par%width, &
									'height: ', this % par%height
               19  format(/,x,3(a7,g15.6,x)/)
                   
                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)
                   
		   elseif ( trim(adjustl(name)) == 'LINE_BY_POINTS' ) then			  
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,17)  &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2,  &
								  'z2: ',this % par%zinj_2
               17   format(/,x,3(a7,g15.6,x),3(a7,g15.6,x)/)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 
	       elseif ( trim(adjustl(name)) == 'READ_PARTICLE_FILE' ) then
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
		           write(unit,*) 
				   write(unit,'(a40,a)') ' reading from particle file: ', trim(adjustl(this%par%file)) 
				   write(unit,*)
	       
	       elseif ( trim(adjustl(name)) == 'READ_CONCENTRATION_FILE' ) then
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
		           write(unit,*) 
				   write(unit,'(a40,a)') 'reading from concentration file: ', trim(adjustl(this%par%file)) 
				   write(unit,*)

		   elseif ( trim(adjustl(name)) == 'LINE_FLUX_WEIGHTED' ) then	          
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,17)   &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2,  &
								  'z2: ',this % par%zinj_2

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 								  
		   elseif ( trim(adjustl(name)) == 'VERTICAL_LINE_FLUX_WEIGHTED' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
	               write(unit,11)    &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'k-up..: ',  this % par%kup

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 
		   elseif ( trim(adjustl(name)) == 'VERTICAL_BLOCK_FLUX_WEIGHTED' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
	               write(unit,11)   &
				                   'i-down: ',  this % par%idwn,  &
								   'j-down: ',  this % par%jdwn,  &
								   'k-down: ',  this % par%kdwn,  &
								   'k-up..: ',  this % par%kup

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)             	
 
		   elseif ( trim(adjustl(name)) == 'CELLS_FILE_FLUX_WEIGHTED' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
		           write(unit,*) 
				   write(unit,'(a40,a)') 'reading injection cells from file: ', trim(adjustl(this%par%file)) 
				   write(unit,*)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)   
                   
		   elseif ( trim(adjustl(name)) == 'CELLS_FILE' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
		           write(unit,*) 
				   write(unit,'(a40,a)') 'reading injection cells from file: ', trim(adjustl(this%par%file)) 
				   write(unit,*)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)  
                   
		   elseif ( trim(adjustl(name)) == 'CELLS_FILE_PARTICLE_NUMBER' ) then
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
		           write(unit,*) 
				   write(unit,'(a40,a)') 'reading injection cells from file: ', trim(adjustl(this%par%file)) 
				   write(unit,*)

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)  
                   
		   elseif ( trim(adjustl(name)) == 'LINE_BY_POINTS_RANDOM' ) then			  
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a15)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,17)   &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2,  &
								  'z2: ',this % par%zinj_2

                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)       
                   
 		   elseif ( trim(adjustl(name)) == 'HORIZONTAL_PLANE_FLUX_WEIGHTED' ) then			  
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a30)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,18)   &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2
                   
                   18  format(/,x,3(a7,g15.6,x),2(a7,g15.6,x)/)
                                  
                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)
                   
                   
 		   elseif ( trim(adjustl(name)) == 'WATER_TABLE_FLUX_WEIGHTED' ) then			  
                   if (this%TypeInj == 'DIRAC') write(unit,'(a26,x,i10)')   'number of Particles.....: ',this%np
		           write(unit,'(a26,x,g16.6)') 'particle mass...........: ',this%pmass
	               write(unit,'(a26,x,a30)')   'source type.............: ',this%name
				   write(unit,'(a26,x,a15)')   'injection type..........: ',this%TypeInj
	               write(unit,'(a26,x,a15)')   'source parameters.......: '
                   write(unit,18)   &
				                  'x1: ',this % par%xinj_1,  &
								  'y1: ',this % par%yinj_1,  &
								  !'z1: ',this % par%zinj_1,  &
                                  'x2: ',this % par%xinj_2,  &
								  'y2: ',this % par%yinj_2
                   
                   20  format(/,x,3(a7,g15.6,x),2(a7,g15.6,x)/)
                                  
                   if (this%TypeInj == 'GENERAL') call print_time_function_ (fname,this%timefunct)
                   
		   else
		      stop '**** could not match name of injection ****'

           end if
 
	end subroutine
    



   end module

!********************************************************************************************
!  VECTOR OF SOURCES CLASS 
!********************************************************************************************
 module source_vect_class
 use source_class

 implicit none

 public !everything is public (inheritance of source_class)
 public :: source_vect_cl           !class
 public ::                                    & !methods
           alloc_source_vect_               , &
           delete_source_vect_            !  , &
           
 type source_vect_cl
      integer                  :: nsource = 0
      type(source_cl), pointer :: num(:)  => null()
 end type

 contains


 subroutine alloc_source_vect_ (this,n)
    implicit none
    type(source_vect_cl), intent(inout) :: this
    integer,              intent(in)    :: n
	allocate(this%num(n))
    this%nsource = n
 end subroutine

 subroutine delete_source_vect_ (this)
    use source_class
	implicit none
    type(source_vect_cl), intent(inout) :: this
	integer                             :: i,n
	if (associated(this%num)) then
	n = size(this%num)
	do i=1,n
	   call delete_source_ (this%num(i)); end do
	end if
 end subroutine

 end module