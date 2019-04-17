
     module read_input_rw3d

	 contains

!********************************************************************************************************
!    SUBROUTINE READ PARAMETERS
!********************************************************************************************************
     subroutine read_parameters  ( geo, advection, dispersion, mass_trans, reaction, plume, well, plane,  &
	                               histo, source, recirculation, StressTimes, obs_block )

     use gslib, only: generate_unit,ifcharacter,generate_name,remove_repeated_values,sortem
     use global_variables
     use velocity_times
	 use code_options
	 use heterogeneity_flags

     use gslib, only: upper_case, open_fname, open_fname_normal
	 use plane_vect_class
	 use well_vect_class
	 use source_vect_class
	 use plume_class
	 use histogram_class
	 use geometry_class
	 use advection_class
	 use dispersion_class
	 use mass_trans_class
	 use reaction_class
	 use array_class
	 use breakthru_class
	 use well_class	 
	 use to_solve !should not be here, delete later
	 use recirculation_class
     use StressTime_class
     use block_vect_class
	 
	 implicit none

	 type(source_vect_cl),       intent(out) :: source
     type(plane_vect_cl),        intent(out) :: plane
	 type(well_vect_cl),         intent(out) :: well
	 type(plume_cl),             intent(out) :: plume
	 type(histo_cl),             intent(out) :: histo(2)
	 type(geometry_cl),          intent(out) :: geo   
	 type(advection_cl),         intent(out) :: advection 
	 type(dispersion_cl),        intent(out) :: dispersion
	 type(mass_trans_cl),        intent(out) :: mass_trans
	 type(reaction_cl),          intent(out) :: reaction
	 type(recirculation_cl),     intent(out) :: recirculation
     type(StressTime_cl),        intent(out) :: StressTimes
     type(block_vect_cl),        intent(out) :: obs_block

	 character(len=200)                :: fname,file,fsour,TypeInj,string,fileTimeShots
	 character(len=1)                  :: firstchar
	 character(len=2)                  :: tipo
	 real*8                            :: const,dist
	 real*8                            :: totmass,tlen,tmult,xr0,yr0
	 integer                           :: i,ivar,iwell,nwell,iplane,nplane,flag,ninj,kinj,iunit,iblock,nblock
	 integer                           :: iunitVEL
	 integer                           :: ntstep,nx,ny,nz
	 logical                           :: exists
	 character(len=15)                 :: Kernel
	 integer                           :: ngrid,inc
	 real*8                            :: bw,tmin,tmax	 
	 integer                           :: nptot
     real*8                            :: vxav, vyav
     integer                           :: nper,iper,nstp,ispe,izone,npl
     real*8                            :: perlen,tsmult,ti,tf
     character(len=2)                  :: iss
     logical                           :: OC_plume_file
     character(len=10), pointer        :: name(:) => null()
     integer                           :: isp,isum,nt,itime,nst
     real*8                            :: v1,v2
    

     
!.....Read name file

      call getarg(1,fname)
      if(fname.eq.' ') then
	      exists=.false.
	      do while (.not.exists) 
          write (*,*) ' enter the name file: '
          read (*,'(a)') fname
	      if (trim(fname) == ' ') fname='rw3d.nam'
		  inquire (file=fname,exist=exists)
		  if (.not.exists) print *, 'file name does not exist'
		  end do      
      end if

      open(11,file=trim(fname))

      write(*,*) ' file name: ',trim(fname)
	  write(*,*)

!.....Read output files
      
	  i=0
	  do while (i < nfnam )
	  read(11,'(a1)',advance='no') firstchar
	    if (firstchar.eq.'#') then
	       read(11,*)
	    else
		   i = i + 1
	       read(11,*) files_nam(i)
		   files_nam(i)=firstchar//files_nam(i) 
        end if
	  end do

	  close(11)
      fdbg = trim(files_nam(16))
      fname_exit = trim(files_nam(17))
      
!......Reading input:

      open(10,file=trim(files_nam(1)),status='old')    ! Parameter file

      do i=1,9; read(10,*); end do

!     Read basic inputs (debug,species,phase,time)

	  read(10,*) idebug   !Debugging Degree
      
      !read species and organize arrays
      
 	  read(10,*) nspe_aq, nspe_min
 	  nspecie = nspe_aq + nspe_min
 	  allocate (namespecie(nspecie),phasespecie(nspecie))
 	  do i=1,nspe_aq
 	     phasespecie(i) = 0
 	  end do
 	  do i=nspe_aq+1,nspecie
 	     phasespecie(i) = 1
 	  end do
	  if (nspe_aq>0) then
	      read(10,*) (namespecie(i),i=1,nspe_aq)
	  else
	      read(10,*)
	  end if
	  
	  if (nspe_min>0) then
	     read(10,*) (namespecie(i),i=nspe_aq+1,nspecie)
	  else
	     read(10,*)
	  end if
	  do i=1,nspecie
	     namespecie(i)=upper_case(namespecie(i))
	  end do

	  read(10,*) tsim           !Simulation time
      !read(10,*) simul_mode     !Simulation mode: forward or backward
	  ipReStart = 1             !Restart Option (particle to restart)
	  SaveMemo = .FALSE.        !Save memory option

!     Read Geometry:

      read(10,*); read(10,*); read(10,*)
      read(10,*) nx,ny,nz;  geo%nx=nx;  geo%ny=ny;  geo%nz=nz
      read(10,*) file,const,ivar,flag;  geo%dx   = read_array_ (file,ivar,const,flag,nx,1,1)
      read(10,*) file,const,ivar,flag;  geo%dy   = read_array_ (file,ivar,const,flag,1,ny,1)
      read(10,*) file,const,ivar,flag;  geo%dz   = read_array_ (file,ivar,const,flag,1,1,nz)
	  read(10,*) file,const,ivar,flag;  call read_inactive_cells_ (geo,file,ivar,const,flag,nx,ny,nz)
	  call read_boundary_    ( geo, 10 )
	  call mesh_coordinates_ ( geo )


!     Read Time Discretization:

      read(10,*); read(10,*); read(10,*)
	  read(10,*) calcul_time_method
	  read(10,*) DtStep,courant,peclet,DaKINETIC,DaDECAY,DaMMT
	             calcul_time_method = upper_case (calcul_time_method)
     
!     Read Advection Parameters:

      read(10,*); read(10,*); read(10,*)
      read(10,*) advection%action
	  read(10,*) method_advection; method_advection = upper_case (method_advection)
      read(10,*) file,const,ivar,flag;

           dataqx%file  = file
           dataqx%const = const
           dataqx%ivar  = ivar
           dataqx%flag  = flag        
      
      ! Initialize modflow variables
      
      if (advection%action) then     
         if (flag == 2) then
            call allocate_mf2k_ (mf2k)
            call check_budget_file_ (file,mf2k)
         end if
      end if
      
      ! Initialize velocity field
   
      if (advection%action) then     
	     if (.not.associated(advection%qx%values)) advection%qx = make_array_ (0.d0,nx+1,ny,  nz )
	     if (.not.associated(advection%qy%values)) advection%qy = make_array_ (0.d0,nx,  ny+1,nz )
	     if (.not.associated(advection%qz%values)) advection%qz = make_array_ (0.d0,nx,  ny,  nz+1)
      end if

      ! Keep reading   
    
      read(10,*) file,const,ivar,flag
           dataqy%file  = file
           dataqy%const = const
           dataqy%ivar  = ivar
           dataqy%flag  = flag      
      read(10,*) file,const,ivar,flag                  
           dataqz%file  = file
           dataqz%const = const
           dataqz%ivar  = ivar
           dataqz%flag  = flag      

      read(10,*) file,const,ivar,flag;  advection%poro = read_array_ (file,ivar,const,flag,nx,ny,nz)

	  read(10,*) nper
	  
	  allocate(velotime%period(nper))
	  
      velotime%NPER = nper
      
	  do iper=1,nper
	     read(10,*) perlen,nstp,tsmult
	     velotime%period(iper)%PERLEN = perlen
	     velotime%period(iper)%NSTP   = nstp
	     velotime%period(iper)%MULT   = tsmult
      end do
      
      read(10,*) velotime%loop_per 
      
      call calculate_advection_timeshots_ (advection)

!     Read Dispersion Parameters:

      read(10,*); read(10,*); read(10,*)
      read(10,*) dispersion%action
	    if (dispersion%action) call alloc_dispersion_ (dispersion)
      read(10,*) file,const,ivar,flag
	    if (dispersion%action) dispersion%aL  = read_array_ (file,ivar,const,flag,nx,ny,nz)
      read(10,*) file,const,ivar,flag
	    if (dispersion%action) dispersion%aTH = read_array_ (file,ivar,const,flag,nx,ny,nz)
      read(10,*) file,const,ivar,flag
	    if (dispersion%action) dispersion%aTV = read_array_ (file,ivar,const,flag,nx,ny,nz)
	  read(10,*) file,const,ivar,flag 
	    if (dispersion%action) dispersion%dm = read_array_ (file,ivar,const,flag,nx,ny,nz)  
	  read(10,*) file,const,ivar,flag 
	    if (dispersion%action) dispersion%dmTH = read_array_ (file,ivar,const,flag,nx,ny,nz)  
	  read(10,*) file,const,ivar,flag 
	    if (dispersion%action) dispersion%dmTV = read_array_ (file,ivar,const,flag,nx,ny,nz)

      if (dispersion%action) then
        read(10,*) (dispersion%MultA(isp),isp=1,nspe_aq)
        read(10,*) (dispersion%MultD(isp),isp=1,nspe_aq)
        do isp=1,nspe_aq
            if (dispersion%MultA(isp) /= 1.d0) SpeciesDispersionDependent =.TRUE.
            if (dispersion%MultD(isp) /= 1.d0) SpeciesDispersionDependent =.TRUE.
        end do
      else
        read(10,*); read(10,*)
      end if

      
!     Read Mass Transfer Parameters:

      read(10,*);  read(10,*); read(10,*)
      mass_trans = read_mass_trans_ (10,geo)
      if (mass_transACTION) nzone = nzone + nzoneim

!     Read Chemical Reaction Parameters:

      read(10,*);  read(10,*); read(10,*)
      reaction = read_reaction_ (10,geo)

!     Read Pumping Wells: 

      read(10,*); read(10,*); read(10,*)        
	  read(10,*) nwell
		if(nwell>0)  call alloc_well_vect_ (well,nwell,nspecie)
	    do iwell=1,nwell
		  call read_assign_well_ (well%num(iwell),10)
		  call get_cell_well_ (well%num(iwell),geo)
        end do
      call getQwell_ (well,geo,10) !for weak-source/sink velocity interpolation scheme

!     Read Control Planes:

	  read(10,*) nplane
		if(nplane>0)  call alloc_plane_vect_ (plane,nplane,nspecie)
        do iplane=1,nplane
	      call read_assign_plane_ (plane%num(iplane),10)
        end do

!     Read Control Blocks:

	  read(10,*) nblock
		if(nblock>0)  call alloc_block_vect_ (obs_block,nblock,nspecie)
        do iblock=1,nblock
	      call read_assign_block_ (obs_block%num(iblock),10)
        end do
        
!     Read Injections:

      read(10,*); read(10,*); read(10,*)
 	  read(10,*) ninj
	    if (ninj > 0) call alloc_source_vect_ (source,ninj)

      do kinj=1,ninj
	       read(10,*) fsour,TypeInj
	       fsour   = upper_case(fsour)
	       TypeInj = upper_case(TypeInj)
		   call read_source_ (source%num(kinj),fsour,TypeInj,files_nam(1))
      end do

!     calculate StartTimeInjection

      StartTimeInjection = UNEST 
      EndTimeInjection   = UNEST    
      do kinj=1,ninj
           ti = source%num(kinj)%TimeStartInj
           tf = source%num(kinj)%TimeStopInj
           if (ti < StartTimeInjection .or. StartTimeInjection == UNEST  ) StartTimeInjection =  ti
           if (tf > EndTimeInjection   .or. EndTimeInjection   == UNEST  ) EndTimeInjection   =  tf           
      end do       


!     Read recirculation

      read(10,*); read(10,*); read(10,*)

      recirculation = read_recirculation_ (10,well)

!     Read Parameters for Snapshots:

      read(10,*); read(10,*); read(10,*)

	  read(10,*)

!     Options for moments and particle cloud at snaphots

	  read(10,*) ixmom               !Print Cartesian Spatial Moments at Snapshots	
	  !read(10,*) irmom,xr0,yr0       !Print Radial Spatial Moments at Snapshots, Xr0, Yr0 (origin coordinates)
	  irmom=0; xr0=0.d0;yr0=0.d0 
	  read(10,*) iwcshot             !Print Particle Cloud at Snapshots

      if (ixmom /= 0 .or. irmom /=0 .or. iwcshot /=0) OC_plume = .TRUE.

!     Create times for plume snaphots

      if (OC_plume) then
	          read(10,*) fileTimeShots
		      OC_plume_file = ifcharacter (fileTimeShots)
		      if (OC_plume_file) then
		         string = upper_case (fileTimeShots)
		         if (trim(adjustl(string)) =='ALWAYS') then
		           AlwaysPrintShots = .TRUE.
		         else
		           AlwaysPrintShots = .FALSE.
                   call open_fname_normal (fileTimeShots,iunit)
				   read(iunit,*)
				   read(iunit,*) ntstep
				   allocate (tshot(ntstep))
				   do i=1,ntstep
				     read(iunit,*) tshot(i)
				   end do
				 end if
		      elseif (.not.OC_plume_file) then
		         backspace(10)
	             read(10,*) tlen,ntstep,tmult 
	             allocate (tshot(ntstep))
			     call time_snapshots (tshot,tlen,ntstep,tmult)
		      end if            
      else
         !read(iunit,*)
         read(10,*)
      end if

      if(allocated(tshot)) ntshot = size(tshot)

!     Read breakthru Options:

      read(10,*)	  	         
	  read(10,*) itmom        !Print Temporal Moments of BTCs  
!	  read(10,*) ixmompl      !Print Moments of Particle Positions at Planes 
!	  read(10,*) iwcshotpl    !Print Particle Position at Planes 	
!	  read(10,*) ipldisp      !Print Apparent Dispersivity at Planes   
      ixmompl = 0;iwcshotpl=0; ipldisp=0;  
      read(10,*) iwbtc,ngrid,Kernel,bw,tmin,tmax;   if (iwbtc  == 1) call assign_histo_PDF_ (histo(1),ngrid,Kernel,bw,tmin,tmax)    !breakthru curve parameters and histogram options
	  read(10,*) iwcbtc,inc;                        if (iwcbtc == 1) call assign_histo_CDF_ (histo(1),inc)                          !cumulative breakthru curve parameters and frequency of printing
      !read(10,*) iwDbtc,ngrid,Kernel,bw,tmin,tmax;  if (iwDbtc == 1) call assign_histo_PDF_ (histo(2),ngrid,Kernel,bw,tmin,tmax)   !Derivative of breakthru curve parameters and histogram options
      iwDbtc = 0
      read(10,*) ibcbtc,conv;                       if (ibcbtc == 1) call assign_histo_CDF_ (histo(1),inc)                          !number of particles in a block and conversion into concentration

      read(10,*) iwpath,pathfreq,pathpart                  !option write particle paths and frequency of printing

      if (itmom   /= 0 .or. ixmompl /=0 .or. iwcshotpl /=0 .or.    &
	      ipldisp /= 0 .or. iwbtc   /=0 .or. iwcbtc    /=0 .or.    &
		  iwDbtc  /= 0) OC_btc = .TRUE.

!     Read Block Region Geometry:

      read(10,*)
	  read(10,*)
	  read(10,*)

	  irestime = .FALSE.
	  
	  !read(10,*,err=2,end=2) irestime

      if(irestime == .TRUE.) then
!	     print *, '  Residence Time Mode '
!	     print *, ' '
!	     read(10,*) nx,ny,nz;  geoblock%nx=nx;  geoblock%ny=ny;  geoblock%nz=nz
!         read(10,*) file,const,ivar,flag;  geoblock%dx   = read_array_ (file,ivar,const,flag,nx,1,1)
!         read(10,*) file,const,ivar,flag;  geoblock%dy   = read_array_ (file,ivar,const,flag,1,ny,1)
!         read(10,*) file,const,ivar,flag;  geoblock%dz   = read_array_ (file,ivar,const,flag,1,1,nz)
	  end if
	2 continue


!     Allocate array to save number of particles per species/zone

      allocate (np(nspecie,0:nzone-1))

      np = 0

!..... end of reading

!...... close input file

      close(10)

!.......... check input data is a reasonable problem: 

	  if (nwell==0 .and. nplane==0 ) then; iwbtc = 0; itmom = 0; endif

      if (ixmom==0 .and. irmom==0 .and. iwcshot==0 .and. iprofile==0) ntstep=0

	  !if (ntstep==0 .or. tlen==0.d0 .or. tmult==0.d0) then
	  !      ixmom   = 0
	  !      irmom   = 0
	  !      iwcshot = 0; end if

	  if (ixmom==0    .and. irmom==0   .and. iwcshot==0   .and.  &
	      itmom==0    .and. ixmompl==0 .and. iwcshotpl==0 .and.  &
		  ipldisp==0  .and. iwbtc==0   .and. iwcbtc==0    .and.  &
		  iprofile==0 .and. iwpath==0. .and. irestime==.FALSE.) then
          stop '*** stop: no reason to run program ***'; end if

!.......... create plumes and initialize some objects:

      call initialize_plume_ (plume,StartTimeInjection,nspecie,nzone) 


!.......... calculate stress times due to change in injection or flow velocity 

      StressTimes = calculate_StressTimes_ (source,advection,tsim)

      call find_loc_StressTime_ (StressTimes,plume%time)


!.......... inquire which parameters are homogeneous:

        if(reaction%action)   react_homogeneous = inquire_homogeneity_react_ (reaction)
		if(dispersion%action) disp_homogeneous  = inquire_homogeneity_disp_  (dispersion)
		!if(advection%action)  vel_homogeneous   = inquire_homogeneity_vel_   (advection)
		if(advection%action)  poro_homogeneous  = inquire_homogeneity_poro_  (advection)

!.......... inquire dimensionality of the problem

        !vxav = sum(advection%qx%values)/size(advection%qx%values)
        !vyav = sum(abs(advection%qy%values))/size(advection%qy%values)
        
        ndim = 0
        if(geo%nx > 1) then
             ndim= ndim + 1
             ActiveDim(1) = .TRUE.

             if(geo%ny > 1) then
                 ndim= ndim + 1
                 ActiveDim(2) = .TRUE.

                 if(geo%nz > 1) then
                    ndim= ndim + 1
                    ActiveDim(3) = .TRUE.

                 end if
              end if
        end if

!.......... initialize temporary files for Save-Memory-Mode:

        if (SaveMemo .and. OC_plume) then
		        if (ntstep >0) then
			    allocate (fname_plume(ntstep))
			    do i=1,ntstep
                  call generate_name (i,files_nam(4),fname_plume(i),'plume')
			      select case (ipReStart==1)
				     case (.TRUE.)
					     call open_fname_normal (fname_plume(i),iunit)
			             write(iunit,*) 
			             write(iunit,3) 'SNAPSHOT'    , &
			                            '    TIME    ', &
			                            '     XP     ', &
						   	            '     YP     ', &
							            '     ZP     ', &
							            '    MASS    ', &
							            '   RETARD   ', &
							            ' ZONE'       , &
							            'SPECIE'      , &
										'  PART'
                       3 format(a9,6(a20,x),a6,x,a8,x,a6)
	                     write(iunit,*)
				     case (.FALSE.)
					     call open_fname (fname_plume(i),iunit)
				  end select
			    end do
			    end if
	    end if
        
        if (SaveMemo .and. OC_btc) then
		        if (nwell >0) then
				allocate (fname_well(nwell))
				do iwell=1,nwell
                   call generate_name (iwell,files_nam(3),fname_well(iwell),'well')
				   select case (ipReStart==1)
				      case (.TRUE.)
					     call open_fname_normal (fname_well(iwell),iunit)
			             write(iunit,*) 
			             write(iunit,4) ' CONTROL NUM   ', &
			                            '   PARTICLE    ', &
			                            'ARRIVAL TIME   ', &
							            ' PART MASS     ', &
							            '  INJECTION    ', &
							            '   MOVE        ' 
	                     write(iunit,*)
		        	 4   format(6(a15,x))
	                     write(iunit,*)
				      case (.FALSE.)
					     call open_fname (fname_well(iwell),iunit)
				   end select
				end do
				end if

				if (nplane >0) then
                allocate (fname_plane(nplane))
				do iplane=1,nplane
				   call generate_name (iplane,files_nam(3),fname_plane(iplane),'plane')
				   select case (ipReStart==1)
				      case (.TRUE.)
					    call open_fname_normal (fname_plane(iplane),iunit )
			            write(iunit,*) 
			            write(iunit,5) ' CONTROL NUM   ', &
			                           '   PARTICLE    ', &
			                           'ARRIVAL TIME   ', &
						    	       '     DXp       ', &
							           '     DYp       ', &
							           '     DZp       ', &
							           '    PART MASS  ', &
							           '       ZONE    ', &
							           '     SPECIE    ' 
                    5   format(9(a15,x))	
			            write(iunit,*)
				      case (.FALSE.)
					    call open_fname (fname_plane(iplane),iunit )
				   end select
			 	end do
				end if
		endif


!..........initialize output files:
        
		call initialize_output_files

!....... write to debugging file:

	    iunit = generate_unit(100)
        open(iunit,file=fdbg)                                 
        
        write(iunit,*) '------------------'
        write(iunit,*) '  DEBUGGING FILE  '
        write(iunit,*) '------------------' 
        
        write(iunit,*)
        write(iunit,*) 'Problem Definition:'
        write(iunit,*)
 	    write(iunit,*)              'Number of species......................: ',nspecie
 	    write(iunit,*)
	    do i=1,nspecie
	    write(iunit,'(a16,i3,a,a)') ' Name of species ',i  ,'.....................: ',namespecie(i)
	    end do
	    write(iunit,*)
	    do i=1,nspecie
	    write(iunit,'(a17,i2,a,i2)') ' Phase of species ',i ,'.....................: ',phasespecie(i)	    
	    end do
	    write(iunit,*)
	    write(iunit,*)            'Simulation time........................: ',tsim
        write(iunit,*)        

		write(iunit,*) 'Output options: '
		write(iunit,*)
		write(iunit,*) 'Calculate spatial moments .............: ',ixmom
		!write(iunit,*) 'Calculate spatial radial moments ......: ',irmom
		write(iunit,*) 'Calculate temporal moments ............: ',itmom
		write(iunit,*) 'Write breakthrough curve values .......: ',iwbtc
		write(iunit,*) 'Write cumulative breakthrough curves...: ',iwcbtc
        !write(iunit,*) 'write Derivative breakthrough curves...: ',iwDbtc
		write(iunit,*) 'write evolution of particles ..........: ',iwcshot
		write(iunit,*) 'write particle paths ..................: ',iwpath
		!write(iunit,*) 'spatial moments at control planes......: ',ixmompl
		write(iunit,*) 'write particle position at planes......: ',iwcshotpl
		!write(iunit,*) 'write effective parameters at planes...: ',ipldisp
		write(iunit,*)

        write(iunit,*) 'Active Packages:'
        write(iunit,*)
        write(iunit,*) 'Advection..............................: ',advection%action
        write(iunit,*) 'Dispersion.............................: ',dispersion%action
        write(iunit,*) 'Mass Transfer..........................: ',mass_trans%action
        write(iunit,*) 'Sorption...............................: ',reaction%sorption%action
        write(iunit,*) 'Decay network..........................: ',reaction%decay%action
        write(iunit,*) 'Kinetic reactions......................: ',reaction%kinetic%action
        write(iunit,*)


	    write(iunit,*) 'Boundary flags: '
	    write(iunit,*)
	    write(iunit,*)  'x-direction: ',geo%ib(1,1),geo%ib(1,2)
	    write(iunit,*)  'y-direction: ',geo%ib(2,1),geo%ib(2,2)
	    write(iunit,*)  'z-direction: ',geo%ib(3,1),geo%ib(3,2)
        write(iunit,*)
        write(iunit,*) 'Geometry:'
        write(iunit,*) 
        write(iunit,*) 'Number cells x-direction...: ',geo%nx
	    write(iunit,*) 'Number cells y-direction...: ',geo%ny
	    write(iunit,*) 'Number cells z-direction...: ',geo%nz
		write(iunit,*)

	    call list_array_1D_ (geo%dx,fdbg,' cell size x-direction')
		call list_array_1D_ (geo%dy,fdbg,' cell size y-direction')
		call list_array_1D_ (geo%dz,fdbg,' cell size z-direction')

       	call open_fname (fdbg,iunit)
		
		write(iunit,*)
		write(iunit,*) 'Time discretization: '
		write(iunit,*)
		write(iunit,*) 'Method.................................: ', trim(adjustl(calcul_time_method))
        write(iunit,*) 'Grid Courant Number....................: ',courant
		write(iunit,*) 'Grid Peclet Number.....................: ',peclet
		write(iunit,*) 'Time Step..............................: ',DtStep	
		write(iunit,*) 
		write(iunit,*) 'Plume snaphots times:'
		write(iunit,*)
		if (OC_plume_file) then
		write(iunit,*) 'Reading times from file................: ', fileTimeShots
		else
        write(iunit,*) 'ntstep.................................: ', ntstep
		write(iunit,*) 'tlen...................................: ', tlen
		write(iunit,*) 'tmult..................................: ', tmult
		write(iunit,*)
		end if
		
		do i=1,size(tshot)
	         write(iunit,*) i,tshot(i)
	    end do
	    		 
	    write(iunit,*)
        write(iunit,*) 'Number of wells........................: ',nwell
		write(iunit,*)
	    do iwell=1,nwell
            call print_well_ (well%num(iwell),fdbg)
	    end do
		
	    call open_fname (fdbg,iunit)
	    write(iunit,*)
		write(iunit,*) 'Number of control planes...............: ',nplane
		write(iunit,*)
 		do iplane=1,nplane
           call print_plane_ (plane%num(iplane),fdbg)
		end do
		
	    call open_fname (fdbg,iunit)
	    write(iunit,*)
		write(iunit,*) 'Number of injections.............: ',ninj
		write(iunit,*)
        do kinj=1,ninj
		   call open_fname (fdbg,iunit)
           write(iunit,*) 'Injection number',kinj
		   write(iunit,*)
		   call print_source_ (source%num(kinj),fdbg)
        end do

        if (advection%action) call print_velocity_timeshots_ (fdbg,advection)
        call print_transp_stress_times_ (fdbg,StressTimes)
		
        call open_fname (fdbg,iunit)
        write(iunit,*)		       
		if (dispersion%action .and. idebug>0) call list_array_ (dispersion%aL,  fdbg, 'Longitudinal Dispersivity.....:')
        if (dispersion%action .and. idebug>0) call list_array_ (dispersion%aTH, fdbg, 'Trans. Disp. Horiz. Plane.....:')
        if (dispersion%action .and. idebug>0) call list_array_ (dispersion%aTV, fdbg, 'Trans. Disp. Vert. Plane......:')
		if (dispersion%action .and. idebug>0) call list_array_ (dispersion%dm,  fdbg, 'Diffusion Coefficient.........:')
        call open_fname (fdbg,iunit)
        write(iunit,*)

	  !if (advection%action.and. idebug >=2) then
   !   !if (advection%action) then
	  !     iunitVEL = generate_unit(100)
	  !     open(iunitVEL,file=trim(files_nam(14)))                 
	  !     write(iunitVEL,'(a49)') 'title="tecplot input file with Darcy velocities" '
	  !     write(iunitVEL,'(a43)') 'variables="x","y","z","qx","qy","qz","qmod" '  
   !   end if



!        if (reaction%action) call list_reaction_ (reaction, fdbg)

        if(irestime) then
          !open(iunit,file=fdbg,access='append')
!          call open_fname (fdbg,iunit)
!		  write(iunit,*) 
!		  write(iunit,*) '--------------------------------'
!		  write(iunit,*) '  GEOMETRY FOR RESIDENCE TIMES'
!		  write(iunit,*) '--------------------------------'
!		  write(iunit,*)
!          write(iunit,*) 'number blocks x-direction...: ',geoblock%nx
!	      write(iunit,*) 'number blocks y-direction...: ',geoblock%ny
!	      write(iunit,*) 'number blocks z-direction...: ',geoblock%nz
!		  write(iunit,*)
!
!	      call list_array_ (geoblock%dx,fdbg,'cell size x-direction')
!		  call list_array_ (geoblock%dy,fdbg,'cell size y-direction')
!		  call list_array_ (geoblock%dz,fdbg,'cell size z-direction')
        end if          
	          
        close(iunit)
      
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





!****************************************************************************************
!      calculates times at snapshots of particles
!      the length of a time step is calculated by multiplying the
!      length of the previous time step by "tmult"
!****************************************************************************************
      subroutine time_snapshots (tshot,tlen,ntstep,tmult)
	  implicit none
      
	  real*8,  intent(in) :: tlen,tmult
	  integer, intent(in) :: ntstep

	  real*8, intent(out) :: tshot(ntstep)

	  integer :: it
	  real*8 :: dtf

	  if (ntstep.eq.0 .or. tlen.eq.0.d0 .or. tmult.eq.0.d0) then
         return
	  else if (tmult.gt.1) then
         tshot(1)=tlen*((tmult-1.d0)/((tmult**ntstep)-1.d0))
	     dtf=tshot(1)
	     do it=2,ntstep
	        dtf=dtf*tmult
	        tshot(it)=tshot(it-1)+dtf
	     end do
	  else if (tmult.le.1.) then
	     tshot(1)=tlen/float(ntstep)
		 do it=2,ntstep
		    tshot(it)=tshot(it-1)+tshot(1)
		 end do
	  end if
	
      end subroutine

  end module


!************************************************************************************************************************************************************************************
!                 OPEN OUTPUT FILES AND INITIALIZE
!************************************************************************************************************************************************************************************
  subroutine initialize_output_files 
      use gslib, only: generate_unit,open_fname_normal
	  use code_options
	  use global_variables, only: files_nam,fname_exit
	  implicit none
	  integer :: iunit
	  logical :: existeix

 !...........open debugging file to save particles exiting the system 

        if (idebug>=0) then
	        iunit = generate_unit(100)
            open(iunit,file=fname_exit)                                 
	        write(iunit,*) 'TITLE="TECPLOT input file: Exit Particles" '
	        write(iunit,*) 'VARIABLES="Xp" "Yp" "Zp" "SPECIES" "ZONE" "PARTID" "TIME" '
		end if

!...........Output Files

      if (iwbtc == 1) then
	       iunit = generate_unit(100)
           open(iunit,file=trim(files_nam(2)))                                 
           write(iunit,*) 'TITLE="TECPLOT input file: Breakthrough Curves" '
	       write(iunit,*) 'VARIABLES="Time" "Concentrations" '		  
		   close(iunit)
      end if

      if (iwDbtc == 1) then
	       iunit = generate_unit(100)
           open(iunit,file=trim(files_nam(15)))                                 
           write(iunit,*) 'TITLE="TECPLOT input file: Derivative BTCs" '
	       write(iunit,*) 'VARIABLES="Time" "PDF" "Derivative" '		  
		   close(iunit)
      end if


      if (iwcbtc == 1) then
           iunit = generate_unit(100)
		   open(iunit,file=trim(files_nam(3))) 
           write(iunit,*) 'TITLE="TECPLOT input file: Cumulative Breakthrough Curves" '
	       write(iunit,*) 'VARIABLES="Time" "Conc" "Part"'		  		                                   
		   close(iunit)
      end if
      
      
      if (ibcbtc == 1) then
           iunit = generate_unit(100)
		   open(iunit,file=trim(files_nam(18))) 
           write(iunit,*) 'TITLE="TECPLOT input file: Cumulative Breakthrough Curves" '
	       write(iunit,*) 'VARIABLES="Time" "Np"'
		   close(iunit)
      end if

      if (iwcshot /= 0) then
           iunit = generate_unit(100)
		   open(iunit,file=trim(files_nam(4)))        ! write particle evolution 
	       write(iunit,*) 'TITLE="TECPLOT input file: Particle Evolution" '
	       write(iunit,*) 'VARIABLES= "Xp" "Yp" "Zp" "Mp" "Rp" "ZONE" "SPECIE" "ID" '
           !close(iunit)
	  end if

      if (iwpath == 1 ) then
	       iunit = generate_unit(100)
		   open(iunit,file=trim(files_nam(5)))        ! write particle paths
           write(iunit,*) 'TITLE="TECPLOT input file: Particle Paths" '
	       write(iunit,*) 'VARIABLES= "ID" "Xp" "Yp" "Zp" "Rp" "ZONE" "SPECIE" "TIME" "VELO"'
		   close(iunit)
      end if

      if (ixmom == 1) then
           iunit = generate_unit(100)
		   open(iunit,file=trim(files_nam(6)))            !  write rw3d_cart_spat_mom.out
	       write (iunit,'(/10(a10,6x),a10,(x,a15))') 			              &
               'Time','Xg','Yg','Zg','Sxx','Syy','Szz','Sxy','Sxz','Syz','Particles','Mobile Mass'
		   write (iunit,*)
		   close(iunit)
      end if

      if (iwcshotpl == 1 ) then
	       iunit = generate_unit(100)
	       open(iunit,file=trim(files_nam(8)))                 ! write rw3d_xyz_plane_location.out
	       write(iunit,*) 'TITLE="tecplot input file particle position at planes" '
	       write(iunit,'(a46)') 'VARIABLES="Xp(t)" "Yp(t)" "Zp(t)" "Tp(x)" "Mp"' 
		   close(iunit)
      end if

      if (irmom == 1) then			  	   		  
	      iunit = generate_unit(100)
		  open(iunit,file=trim(files_nam(10)))            ! write rw3d_radial_spat_mom.out
	      write(iunit,'(/10(a15,1x),a10,(x,a15))') 'Time','Rg','Tg','Zg','Srr','Stt','Szz','Srt','Srz','Stz','no. part.'
		  close(iunit)
      end if

      if (itmom == 1) then
          iunit = generate_unit(100)
		  open(iunit,file=trim(files_nam(11)))            ! write rw3d_temp_mom.out
	      write(iunit,'(/a/)') 'temporal moments'
	      write(iunit,'(10(a20,1x),/)')   '  NORMALIZED MOM_1  ' , &
		                                     '   ABSOLUTE MOM_2   ' , &
											 '      SKEWNESS      ' , &
											 '      KURTOSIS      ' , &
	                                         '   ABSOLUTE MOM_3   ' , &
											 '   ABSOLUTE MOM_4   ' , &
											 '    CENTRAL MOM_2   ' , &
											 '    CENTRAL MOM_3   ' , &
										     '    CENTRAL MOM_4   ' , &
											 '    NUMBER PARTICLES' 
		  write(iunit,*)
	      close(iunit)
	  end if

      if (ipldisp == 1) then
	     iunit = generate_unit(100)
		 open(iunit,file=trim(files_nam(12)))            ! write rw3d_var_displacement.out
         write(iunit,'(a15,x,4(a9,7x),3(a11,5x),3(a14,2x),a9)')  'Displacement','Ta','A11','A22','A33','VARxx','VARyy','VARzz','XG(t)-XG(0)','YG(t)-YG(0)','ZG(t)-ZG(0)','Particles'
         write(iunit,*)
		 close(iunit)
	  end if



  end subroutine
!************************************************************************************************************************************************************************************
