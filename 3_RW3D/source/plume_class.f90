!***********************************************************************************************
!   PLUME CLASS
!***********************************************************************************************
    module plume_class
    use array_class
    use list_class
	implicit none

	private
	public :: plume_cl     ! class
	public ::                                       & ! methods
			  print_plume_                        , &
			  moments_plume_                      , &
			  print_conc_plume_                   , &
			  print_moments_plume_                , &
			  initialize_plume_                   , &
			  add_particle_to_plume_              , &
			  print_number_of_particles_          , &
			  delete_plumeparticle_               , &
			  generate_plumeparticle_ID_          , &
			  update_particle_state_in_plume_     , &
			  get_plume_np_

	type species_cl
	      type(list_cl),  pointer, dimension(:) :: zone => null()
	end type
	
	type plume_cl
		 integer                                 :: np       = 0             !total number of particles
         real*8                                  :: mtot     = 0             !total mass over all plumes
		 real*8                                  :: time     = 0.d0          !time plume snapshot
		 integer                                 :: nspecie  = 1             !number of species
		 integer                                 :: nzone    = 0             !number of mass transfer zones
		 type(species_cl), pointer, dimension(:) :: species  => null() 
		 type(array_cl),   pointer               :: conc     => null()       !concentrations (I have never used this, better as a postprocessing)
	end type plume_cl
	

    interface print_plume_
	     module procedure print_plume_tecplot_; end interface

	interface delete_plumeparticle_ 
	     module procedure delete_plumeparticle_1,delete_plumeparticle_2; end interface     	
	
    interface print_number_of_particles_ 
         module procedure print_number_of_particles_0,print_number_of_particles_1,print_number_of_particles_2; end interface
	
    contains


   function get_plume_np_ (plume) result (np)
        use global_variables, only: nspecie,nzone
        implicit none
		type (plume_cl), intent (in) :: plume
		integer                      :: np(nspecie,0:nzone-1)
		integer                      :: ispe,izone        
   
        do ispe=1,nspecie
          do izone=0,nzone-1
            np(ispe,izone) = plume%species(ispe)%zone(izone)%np
          end do
        end do
   
   end function



!************************************************************************
!************************************************************************
   subroutine generate_plumeparticle_ID_ (ID)
         implicit none
         integer, intent(out) :: ID
         integer, save        :: IDold = 0
        ID = IDold + 1
        IDold = ID
   end subroutine


!************************************************************************
   function make_plume_conc (const,size_nx,size_ny,size_nz) result (name) !optional constructor
         use array_class
		 integer, optional, intent(in) :: size_nx,size_ny,size_nz
		 real*8,            intent(in) :: const
		 type(plume_cl)                :: name
		 integer                       :: ok
		   if (present(size_nx) .and. present(size_ny) .and. present(size_nz)) then
			  name % conc = assign_array_ (const,size_nx,size_ny,size_nz)
		   else
			  name % conc = assign_array_ (const,1,1,1)
           end if
     end function


!************************************************************************
      subroutine print_plume_tecplot_ (plume,fname)   !printer
          use gslib, only: open_fname,locate
          use code_options, only: iwcshot
          use global_variables, only: StartTimeInjection,tshot,AlwaysPrintShots,ntshot 
		  use list_class
		  character(len=*), optional, intent(in)     :: fname
		  type (plume_cl),            intent (inout) :: plume
		  integer                                    :: itnext = 1
		  type(partID_cl),  pointer                  :: part		   
          logical                                    :: exists
		  integer                                    :: iunit,ispe,izone,nt,i
		  logical, save                              :: EndPrintShot = .FALSE.
		  logical, save                              :: FirstInjection =.TRUE.      
		  
		  if (iwcshot /= 1) return
		  
		  if (EndPrintShot) return
		         
		  if (plume%np <=0) return

          !if always printing

		  if (AlwaysPrintShots) then  !if negative or start time injection always priting to file
          
		  call open_fname(fname,iunit)
		  write(iunit,2)'zone t="',plume%time,'"'
          species1: do ispe  = 1, plume%nspecie
            zones1: do izone = 0, plume%nzone-1       
               part => plume%species(ispe)%zone(izone)%head
               particles1:  do        ! print all particles one time step              
                 if (.not.associated(part)) cycle zones1        
	   	         write(iunit,3) part%xp,part%yp,part%zp,part%mp,part%rp,izone,ispe,part%id 			     
			     part => part%next
			  end do particles1
		    end do zones1
		  end do species1   	  
	      
	      return		  
		  
		  end if

          !if first injection
		  		  
		  if (FirstInjection) then  !if first injection
          
		  call open_fname(fname,iunit)
		  write(iunit,2)'zone t="',plume%time,'"'
          species2: do ispe  = 1, plume%nspecie
            zones2: do izone = 0, plume%nzone-1       
               part => plume%species(ispe)%zone(izone)%head
               particles2:  do        ! print all particles one time step              
                 if (.not.associated(part)) cycle zones2        
	   	         write(iunit,3) part%xp,part%yp,part%zp,part%mp,part%rp,izone,ispe,part%id 			     
			     part => part%next
			  end do particles2
		    end do zones2
		  end do species2   	  
	
          FirstInjection = .FALSE.
          if (plume%time > tshot(itnext)) then                  
                  do i=itnext,ntshot-1
                      if (tshot(i)<=plume%time.and.tshot(i+1)>plume%time) exit 
                  end do
                  if (i<ntshot) itnext = i+1
          end if
	      
	      return		  
		  
		  end if

          !...print when changing time interval  
    
		  if (plume%time >= tshot(itnext)) then
		  	  	         
             call open_fname(fname,iunit)
		     write(iunit,2)'zone t="',plume%time,'"'
	       2 format(a8,g15.6,a2) 		  
             species: do ispe = 1, plume%nspecie
              zones: do izone = 0, plume%nzone-1
                 part => plume%species(ispe)%zone(izone)%head
                 particles:  do        ! move all particles one time step              
                      if (.not.associated(part)) cycle zones        
	   	              write(iunit,3) part%xp,part%yp,part%zp,part%mp,part%rp,izone,ispe,part%id
	               3  format(5(x,e20.11),2(i4,x),i10)			     
			          part => part%next
			     end do particles
		     end do zones
		    end do species   
	      			     
	        itnext = itnext + 1
	        
	        if (itnext > ntshot ) EndPrintShot = .TRUE.
		
		end if
			    
	  end subroutine

!************************************************************************
	  subroutine print_conc_plume_ (this,fname,nombre) !printer
           use array_class
		   implicit none
		   type(plume_cl),             intent(in) :: this
		   character(len=*), optional, intent(in) :: fname,nombre
             if (associated(this%conc)) call print_array_ (this%conc,fname,nombre) !print in gslib format
      end subroutine


!************************************************************************
      function value_conc_plume (name,loci,locj,lock) result (value)        ! accessor member
         use array_class
         implicit none
		 integer, optional, intent(in) :: loci,locj,lock
		 type (plume_cl), intent(in)   :: name
         real*8                        :: value
         value = value_array_ (name%conc,loci,locj,lock)
	  end function 

!************************************************************************
!   Calculate Spatial Moments in Cartesian Coordinates 
!************************************************************************

      subroutine moments_plume_ (plume,ispecie,izone,xg,mom,aqmass)  ! moment calculator, cartesian coordinates
	     implicit none      
		 type (plume_cl), intent(in)   :: plume
		 integer,         intent(in)   :: izone,ispecie
		 real*8,          intent(out)  :: xg(3)
		 real*8,          intent(out)  :: mom(3,3)
		 real*8,          intent(out)  :: aqmass
		 real*8                        :: rpt
		 integer                       :: jj,np
		 type(partID_cl), pointer      :: part

         aqmass = 0.d0

         xg  = 0.d0
         mom = 0.d0

         part => plume%species(ispecie)%zone(izone)%head

         np = plume%species(ispecie)%zone(izone)%np

         do jj=1,np

               rpt = 1.d0

		       if (associated(part%rp) ) rpt = part%rp

               aqmass = aqmass + part%mp / rpt

			   xg(1) = xg(1) + part%mp * part%xp / rpt
	           xg(2) = xg(2) + part%mp * part%yp / rpt
	           xg(3) = xg(3) + part%mp * part%zp / rpt

	           mom(1,1) = mom(1,1) + part%mp * part%xp * part%xp / rpt
	           mom(2,2) = mom(2,2) + part%mp * part%yp * part%yp / rpt
	           mom(3,3) = mom(3,3) + part%mp * part%zp * part%zp / rpt
			   
			   mom(1,2) = mom(1,2) + part%mp * part%xp * part%yp / rpt
	           mom(1,3) = mom(1,3) + part%mp * part%xp * part%zp / rpt
	           mom(2,3) = mom(2,3) + part%mp * part%yp * part%zp / rpt

               part => part%next

            end do

          if (aqmass <= 0.d0) then
		       xg  = 0.d0
		       mom = 0.d0
		  else

	           xg(1) = xg(1) / aqmass
	           xg(2) = xg(2) / aqmass
	           xg(3) = xg(3) / aqmass

	           mom(1,1) = mom(1,1) / aqmass - xg(1) * xg(1)
		       mom(2,2) = mom(2,2) / aqmass - xg(2) * xg(2)
	           mom(3,3) = mom(3,3) / aqmass - xg(3) * xg(3)
	     
	           mom(1,2) = mom(1,2) / aqmass - xg(1) * xg(2)
		       mom(1,3) = mom(1,3) / aqmass - xg(1) * xg(3)
	           mom(2,3) = mom(2,3) / aqmass - xg(2) * xg(3)

		       mom(2,1) = mom(1,2)
		       mom(3,1) = mom(1,3)
		       mom(3,2) = mom(2,3)

		  end if

	  end subroutine


!************************************************************************
 subroutine print_moments_plume_ (plume,fname)     !printer: calculates and print cartesian moments
	     use gslib, only: open_fname
		 implicit none
         type (plume_cl),     intent(in)      :: plume
	     character(len=*),optional,intent(in) :: fname
	     real*8                               :: xg(3)
	     real*8                               :: mom(3,3)
	     real*8                               :: aqmass,tt
         integer                              :: iunit,i,nn,iplume,nplume,ispecie,izone,nspecie,nzone
	     logical                              :: exists
			if (plume%np == 0 ) return
			nspecie = plume%nspecie
			nzone   = plume%nzone
            call open_fname (fname,iunit)
					do ispecie=1,nspecie
				    do izone=0,nzone-1
			            call moments_plume_ (plume,ispecie,izone,xg,mom,aqmass) 
			            tt = plume % time
			            nn = plume % np
		                write(iunit,2) tt,(xg(i),i=1,3),mom(1,1),mom(2,2),mom(3,3),mom(1,2),mom(1,3),mom(2,3),nn,aqmass,'       Specie: ',ispecie
       2                format(10(g15.6,1x),i10,(x,g15.6),a15,i5)
	                end do  
					end do
	        close(iunit)
 end subroutine



!*******************************************************************************
! new
!*******************************************************************************
     subroutine initialize_plume_ (plume,time,ns,nzone)
       implicit none
       type(plume_cl), intent(inout)  :: plume
       real*8,         intent(in)     :: time 
       integer,        intent(in)     :: ns,nzone
       integer                        :: is,izone
       if (associated(plume%species)) return
           plume%nspecie  = ns
           plume%nzone    = nzone
           plume%time     = time
           plume%np       = 0
           plume%mtot     = 0.d0 
           allocate(plume%species(ns))
           do is=1,ns
             if (.not.associated(plume%species(is)%zone)) then
                  allocate(plume%species(is)%zone(0:nzone-1))
                  do izone=0,nzone-1
                    nullify (plume%species(is)%zone(izone)%head)
                    nullify (plume%species(is)%zone(izone)%tail)
                  end do
             end if
           end do
     end subroutine


!************************************************************************
!        delete a plumeparticle
!************************************************************************
    subroutine delete_plumeparticle_1 (plume,izone,ispe,part)
       implicit none
       type(plume_cl),  intent(inout) :: plume
       integer,         intent(in)    :: izone,ispe
       type(partID_cl),  pointer :: part
       type(partID_cl),  pointer :: next            
       real*8                    :: mp   
          mp=part%mp
          !save next particle
          next => part%next
          if (.not.associated(part%prev)) then  
                plume%species(ispe)%zone(izone)%head => part%next  !delete the first item in list
          else
                part%prev%next => part%next 
          end if 
          if (.not.associated(part%next)) then
                plume%species(ispe)%zone(izone)%tail => part%prev
          else
                part%next%prev => part%prev
          end if    
          !delete particle attributes memory
		  if (associated(part%id)) then
              deallocate(part%id)
              nullify(part%id)
          end if
         if (associated(part%xp)) then
              deallocate(part%xp)
              nullify(part%xp)
          end if
         if (associated(part%yp)) then
              deallocate(part%yp)
              nullify(part%yp)
          end if
         if (associated(part%zp)) then
              deallocate(part%zp)
              nullify(part%zp)
          end if
         if (associated(part%rp)) then
              deallocate(part%rp)
              nullify(part%rp)
          end if
         if (associated(part%mp)) then
              deallocate(part%mp)
              nullify(part%mp)
          end if
 		  !delete particle     
          deallocate(part)
          nullify(part)
          part => next
          plume%species(ispe)%zone(izone)%np = plume%species(ispe)%zone(izone)%np - 1
          plume%species(ispe)%zone(izone)%mtot = plume%species(ispe)%zone(izone)%mtot - mp
          plume%np = plume%np - 1
          plume%mtot = plume%mtot - mp         
     end subroutine


!************************************************************************
    subroutine delete_plumeparticle_2 (plume,react,ispe,izone)
       implicit none
       type(plume_cl), intent(inout) :: plume
       integer,        intent(in)    :: ispe,izone
       logical,        intent(in)    :: react(plume%species(ispe)%zone(izone)%np) 
       type(partID_cl), pointer      :: part
       integer                       :: np,i,id
           part => plume%species(ispe)%zone(izone)%head
           np = plume%species(ispe)%zone(izone)%np
           do i=1,np
              if (react(i)) then
                 call delete_plumeparticle_  (plume,izone,ispe,part)
              else
                 part => part%next
              end if
           end do
           
    end subroutine


!************************************************************************
    subroutine add_particle_to_plume_ (plume,id,xp,yp,zp,mp,rp,zone,species) !add a new particle to the end of the list
       implicit none
       type(plume_cl), intent(inout) :: plume
       integer,        intent(in)    :: species,zone
       real*8,         intent(in)    :: xp,yp,zp,mp,rp
       integer,        intent(in)    :: id
       type(partID_cl), pointer      :: part
	   integer                       :: err,err1,err2,err3,err4,err5,err6              
            allocate(part,stat=err)
            nullify(part%next)
            nullify(part%prev)
            allocate(part%id,stat=err1)    
			allocate(part%xp,stat=err2)
			allocate(part%yp,stat=err3)
			allocate(part%zp,stat=err4)
			allocate(part%mp,stat=err5)
			allocate(part%rp,stat=err6)
            part%id = id
			part%xp = xp
			part%yp = yp
			part%zp = zp
			part%mp = mp
			part%rp = rp
			err = err1 + err2 + err3 + err4 + err5 + err6    
 	        if (err /= 0) print *, '*** Allocation of particle in list not succesful ***'
            if(associated(plume%species(species)%zone(zone)%head)) then !plume is not empty
                plume%species(species)%zone(zone)%tail%next => part
                nullify(part%next)
                part%prev => plume%species(species)%zone(zone)%tail 
                plume%species(species)%zone(zone)%tail => part
            else                  !plume is empty
                plume%species(species)%zone(zone)%head => part
                plume%species(species)%zone(zone)%tail => part
                nullify(plume%species(species)%zone(zone)%tail%next)
            end if
            plume%species(species)%zone(zone)%np = plume%species(species)%zone(zone)%np + 1
            plume%species(species)%zone(zone)%mtot = plume%species(species)%zone(zone)%mtot + mp
            plume%np = plume%np + 1
            plume%mtot = plume%mtot + mp     
    end subroutine

!************************************************************************
    subroutine update_particle_state_in_plume_ (particle,part,plume,remove)
    use particle_class
    use global_variables, only: nspecie, nzone
    implicit none

    type(plume_cl),    intent(inout)    :: plume
    type(particle_cl), intent(inout)    :: particle
    type(partID_cl),   pointer          :: part
    logical,           intent(inout)    :: remove
    integer                             :: newzone,prezone,newspe,prespe,ip
    real*8                              :: xp,yp,zp,mp,rp

    newzone = particle%zone%num
    prezone = particle%zone%prenum
    newspe  = particle%specie%num
    prespe  = particle%specie%prenum

    if ( newspe/=prespe .OR. newzone/=prezone ) then !particle changed its state (either species or zone)

        xp = particle%position%xp(1)
        yp = particle%position%xp(2)
        zp = particle%position%xp(3)
        mp = part%mp
        rp = particle%prop%rp

        call generate_plumeparticle_ID_ (ip)
        if ( newspe<=nspecie )  call add_particle_to_plume_ (plume,ip,xp,yp,zp,mp,rp,newzone,newspe)
        particle%control%remove = .TRUE.

    end if

    remove = particle%control%remove .or. particle%control%stuck
    
    if (remove) call delete_plumeparticle_  (plume,prezone,prespe,part)
    

    end subroutine

!************************************************************************
!        function to create species and zone names to print
!************************************************************************
      function create_name (ispe,izone) result (string)
	      use gslib, only: generate_unit
	      implicit none
          integer, intent(in)  :: ispe,izone
		  character(len=10)    :: string1,string2,string
		  integer              :: unit
		  unit = generate_unit(283)
		  open(unit,status='scratch')
		  write(unit,*) ispe,izone
		  rewind(unit)
		  read(unit,*) string1,string2
		  string = '( '//trim(adjustl(string1))//' , '//trim(adjustl(string2))//' )'
		  close(unit)
      end function


!************************************************************************
    subroutine print_number_of_particles_0 (plume,dt,movez)
	  use gslib, only: open_fname,generate_unit
      use loops_particles, only: nmove
      use global_variables, only: nspecie,nzone,npInj,npWell,npPlane,npOut,npStuck
	  implicit none
	  type(plume_cl),    intent(in) :: plume
	  real*8,            intent(in) :: dt
      integer,           intent(in) :: movez
	  integer                                 :: iunit
      logical                                 :: connected,exists
      logical                                 :: FirstTime=.TRUE.
      integer, pointer                        :: np11(:) => null()
      character(len=10), pointer              :: name(:) => null()
      integer                                 :: izone,ispe,npl,i
           iunit=6
           write(iunit,'(a8,x,g15.6,x,a5,g15.6,2(3x,a7,x,i9),x,a7,g15.6)') 'TIME =',plume%time,'DT =',dt,'MOVE =',movez,'NPTOT =',plume%np,'MTOT =',plume%mtot
    end subroutine      
      
      
    subroutine print_number_of_particles_1 (plume,dt)
	  use gslib, only: open_fname,generate_unit
      use loops_particles, only: nmove
      use global_variables, only: nspecie,nzone,npInj,npWell,npPlane,npOut,npStuck
	  implicit none
	  type(plume_cl),    intent(in) :: plume
	  real*8,            intent(in) :: dt
	  integer                                 :: iunit
      logical                                 :: connected,exists
      logical                                 :: FirstTime=.TRUE.
      integer, pointer                        :: np11(:) => null()
      character(len=10), pointer              :: name(:) => null()
      integer                                 :: izone,ispe,npl,i
           iunit=6
           write(iunit,'(a8,x,g15.6,x,a5,g15.6,2(3x,a7,x,i9),x,a7,g15.6)') 'TIME =',plume%time,'DT =',dt,'MOVE =',nmove,'NPTOT =',plume%np,'MTOT =',plume%mtot
    end subroutine

    subroutine print_number_of_particles_2 (plume,fname)
	  use gslib, only: open_fname,generate_unit
      use loops_particles, only: nmove
      use global_variables, only: nspecie,nzone,npInj,npWell,npPlane,npOut,npStuck
	  implicit none
	  type(plume_cl),              intent(in) :: plume
	  character(len=*),            intent(in) :: fname
	  integer                                 :: iunit
      logical                                 :: connected,exists
      logical                                 :: FirstTime=.TRUE.
      integer, pointer                        :: np11(:) => null()
      character(len=10), pointer              :: name(:) => null()
      integer                                 :: izone,ispe,npl,i

	       inquire(file=fname,exist=exists,opened=connected)
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if  	              
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown',access='append')
           end if       
	       npl = nspecie*nzone
	       if (FirstTime) then
	          FirstTime = .FALSE.
	          if (.not.associated(np11)) allocate (np11(npl))
	          if (.not.associated(name)) allocate (name(npl))
	          i=0
	          do ispe=1,nspecie
	             do izone=0,nzone-1
	              i=i+1
	              name(i)= create_name(ispe,izone)
	             end do
	          end do
		          
	          write(iunit,*)
	          write(iunit,'(<40+16+55+6+11*npl>("-"))') 
	          write(iunit,'(95x,16x,a)') ' (SPECIES,ZONES): ' 
	          write(iunit,*)                          
	          write(iunit,'(x,a15,(x,a10),x,a15,x,6(x,a10),<npl>(a11))') '     TIME      ','MOVE','     MTOT      ',    &
                                                                         'NPTOT','NPINJ','NPWELL','NPPLANE','NPOUT','NPSTUCK',(trim(adjustl(name(i))),i=1,npl)
	          write(iunit,'(<40+55+16+6+11*npl>("-"))')                          
	       end if 
	       i=0
	       do ispe=1,nspecie
	         do izone=0,nzone-1
	           i=i+1
	           np11(i)=plume%species(ispe)%zone(izone)%np
	         end do
	       end do
	       write(iunit,1) plume%time,nmove,plume%mtot,plume%np,npInj,npWell,npPlane,npOut,npStuck,(np11(i),i=1,npl) 
	   1   format (x,g15.6,x,(x,i10),x,g15.6,6(x,i10),<npl>(x,i10) )
	   
   end subroutine

!************************************************************************ 
	end module


!****************************************************************************** 
!******************************************************************************	

