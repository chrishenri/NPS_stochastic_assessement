!***********************************************************************************************
!   BREAKTHRU CURVE CLASS
!***********************************************************************************************
    module breakthru_class
	use list_class
	implicit none

	public
	public :: breakthru_cl     ! class
	public ::                                      & ! methods
			  moments_breakthru_                 , &
			  histo_breakthru_                   , &
              print_breakthru_                   , &
			  print_moments_breakthru_           , &
			  print_pdf_breakthru_               , &
			  print_cdf_breakthru_               , &
              initialize_btc_                    , &
              delete_btcparticle_                , &
              add_particle_to_btc_

	type breakthru_cl
		 integer                  :: np = 0             !total number of particles
         type(partID_cl), pointer :: head
         type(partID_cl), pointer :: tail         
	end type breakthru_cl

    interface print_pdf_breakthru_
	     module procedure print_pdf_breakthru_Derivative; end interface
    interface print_breakthru_
	     module procedure print_breakthru_1, print_breakthru_2; end interface

	contains

!*****************************************************************
!  new subourtines for breakthru curve as a linked list
!*****************************************************************

     subroutine initialize_btc_ (btc)
       implicit none
       type(breakthru_cl), intent(inout) :: btc
           btc%np       = 0 
           nullify (btc%head)
           nullify (btc%tail)
     end subroutine

    subroutine delete_btcparticle_ (btc,part)
       use list_class
       implicit none
       type(breakthru_cl), intent(inout) :: btc
       type(partID_cl),  pointer :: part
       type(partID_cl),  pointer :: next            
          !save next particle
          next => part%next
          if (.not.associated(part%prev)) then  
                btc%head => part%next  !delete the first item in list
          else
                part%prev%next => part%next 
          end if 
          if (.not.associated(part%next)) then
                btc%tail => part%prev
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
          btc%np = btc%np - 1         
     end subroutine

      
    subroutine add_particle_to_btc_ (btc,id,tp,mp,xp,yp,zp) !add a new particle to the end of the list
       use list_class
       implicit none
       type(breakthru_cl), intent(inout) :: btc
       real*8,             intent(in)    :: tp,mp
       integer,            intent(in)    :: id
       real*8, optional,   intent(in)    :: xp,yp,zp
       type(partID_cl),     pointer      :: part
	   integer                           :: err,err1,err2,err3,err4,err5,err6             
            err=0; err1=0; err2=0
            err3=0; err4=0;err5=0;err6=0
            allocate(part,stat=err)
            nullify(part%next)
            nullify(part%prev)
            allocate(part%id,stat=err1)    
			allocate(part%tp,stat=err2)
			allocate(part%mp,stat=err3)
            part%id = id
			part%tp = tp
			part%mp = mp	
			if (present(xp)) then
                allocate(part%xp,stat=err4)    
			    allocate(part%yp,stat=err5)
			    allocate(part%zp,stat=err6)
                part%xp = xp
			    part%yp = yp
			    part%zp = zp				
			end if
			err = err1 + err2 + err3 + err4 + err5 + err6      
 	        if (err /= 0) print *, '*** Allocation of particle in btc list not succesful ***'
            if(associated(btc%head)) then !plume is not empty
                btc%tail%next => part
                nullify(part%next)
                part%prev => btc%tail 
                btc%tail => part
            else                  !plume is empty
                btc%head => part
                btc%tail => part
                nullify(btc%tail%next)
            end if
            btc%np = btc%np + 1     
    end subroutine

    subroutine print_number_of_btc_particles_ (btc,fname)
	  use gslib, only: open_fname,generate_unit
      use loops_particles, only: nmove
	  implicit none
	  type(breakthru_cl),          intent(in) :: btc
	  character(len=*), optional,  intent(in) :: fname
	  integer                                 :: iunit
      logical                                 :: connected,exists
       if(.not.present(fname)) then
           iunit=6
       else
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
	   end if
	   write(iunit,'(a8,g15.6,(a7,x,i5))') 'MOVE =',nmove,'NP =',btc%np 
   end subroutine

    subroutine print_breakthru_1 (btc,fname)
         use gslib, only: open_fname
         use list_class
		 implicit none
		 type(breakthru_cl),        intent(in) :: btc
		 character(len=*),optional, intent(in) :: fname
		 type(partID_cl), pointer              :: part 
		 integer                               :: iunit,jj
              if (btc%np <=0) return
              call open_fname(fname,iunit)             
              part => btc%head
              particles:  do                     
                 if (.not.associated(part)) exit       
	   	         write(iunit,2) part%tp,part%mp,part%ID  
			     part => part%next
			  end do particles
		  2  format(2(g15.6,x))
			 close(iunit)
	end subroutine


    subroutine print_breakthru_2 (btc,iunit)
         use gslib, only: open_fname
         use list_class
		 implicit none
		 type(breakthru_cl),        intent(in) :: btc
		 type(partID_cl), pointer              :: part
		 integer,                   intent(in) :: iunit 
		 integer                               :: jj
              if (btc%np <=0) return
              part => btc%head
              particles:  do                     
                 if (.not.associated(part)) exit       
	   	         write(iunit,2) part%tp,part%mp,part%ID  
			     part => part%next
			  end do particles
		  2  format(2(g15.6,x),i10)
	end subroutine


    subroutine print_moments_breakthru_ (name,fname)
	    use gslib, only: open_fname
		implicit none
		type(breakthru_cl)                    :: name
		character(len=*),optional, intent(in) :: fname 
		real*8                                :: alfa(4),mom(4)
	    real*8                                :: skew,kurt
        integer                               :: iunit,np
		if (name%np >=5) then
             call open_fname(fname,iunit)
		     call moments_breakthru_  (name,alfa,mom,skew,kurt)
 	  	     np = name%np
		     write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np
	    else
		     print *, '>> breakthru curve undefined or np < 5'
        end if
	 1  format (9(g20.12,1x),i20)
	end subroutine


    subroutine delete_repeated_particles_in_btc_ (btc)
 		use gslib, only: sortem
    	use list_class
 		implicit none
		type(breakthru_cl),  intent(inout) :: btc
		type(partID_cl), pointer           :: part
		real*8                             :: num(btc%np),ind(btc%np)
        integer                            :: nd,i,ipart
		real*8                             :: aa(1),numold
		logical                            :: remove(btc%np)         
          
           remove = .FALSE.
          
           nd = btc%np
    	   part => btc%head
		   do i=1,nd
		      ind(i) = dfloat(i)
		      num(i) = dfloat(part%ID)
		      part => part%next
		   end do

  		   call sortem(1,nd,num,1,ind,aa,aa,aa,aa,aa,aa)

		   numold = num(1)
		   do i=2,nd
		      ipart = dint(ind(i)+0.5)
		      if (num(i) == numold) remove(ipart)=.TRUE. 
		      numold = num(i)
		   end do


    	   part => btc%head
		   do ipart=1,nd
		      if (remove(ipart)) then
                  call delete_btcparticle_ (btc,part)
              else
		          part => part%next
		      end if
		   end do



    end subroutine

!************************************************************************
!   SUBROUTINE TO ESTIMATE CDF FROM PARTICLE TIMES AND PRINT TO FILE
!************************************************************************
    subroutine print_cdf_breakthru_ (name,histo,fname)
        use gslib, only: open_fname,sortem
		use histogram_class
		use list_class
		implicit none
		type(breakthru_cl),         intent(in)    :: name
		type(histo_cl),             intent(inout) :: histo 
		character(len=*), optional, intent(in)    :: fname
		real*8, allocatable                       :: cdf(:),rvr(:),mp(:),num(:), xp(:),yp(:),zp(:)
		integer                                   :: i,nd,iunit,inc,ipart
		real*8                                    :: cp,oldcp,numold
		real*8                                    :: aa(1) 
		type(partID_cl), pointer                  :: part

		   inc = histo%inc
           nd  = name%np
		   allocate (rvr(nd),mp(nd),num(nd),cdf(nd))
           allocate (xp(nd),yp(nd),zp(nd))
     
		   part => name%head
		   do i=1,nd
		      num(i) = dfloat(part%ID)
		      rvr(i) = part%tp
		      mp(i)  = part%mp
              
              xp(i) = part%xp
              yp(i) = part%yp
              zp(i) = part%zp
              
		      part => part%next
		   end do

!..... sort arrival times

           aa = 1.d0
           
		   !call sortem(1,nd,rvr,2,mp,num,aa,aa,aa,aa,aa)
           call sortem(1,nd,rvr,5,mp,num,xp,yp,zp,aa,aa)

!...... calculates cdf

	       cdf   = 0.0d0
           oldcp = 0.0d0
	       cp    = 0.0d0

	       do i=1,nd
	         cp      = cp + mp(i)
	         cdf(i) = (cp + oldcp) * 0.5d0
	         oldcp   = cp
	       end do

!.......write-out cdf to file

           call open_fname(fname,iunit)
           do i=1,nd,inc
                ipart = dint(num(i)+0.5)
                !write(iunit,'(2(g21.14,x),i10') rvr(i),cdf(i),ipart
                write(iunit,'(2(g21.14,x),i10,3(g21.14,x))') rvr(i),cdf(i),ipart,xp(i),yp(i),zp(i)
		   end do
           write(iunit,*)
		   close(iunit)	

		   deallocate (rvr,mp,num,cdf) 

    end subroutine



!***********************************************************************
!   SUBROUTINE TO CALCULATE TEMPORAL MOMENTS FROM PARTICLE ARRIVAL TIMES
!***********************************************************************
!
!     alfa    : vector of normalized absolute moments
!     fmom    : vector of normalized central moments
!     rvr     : vector of data values (sorted particle arrival times)
!     nd      : number of data values (particles arrived control surface)
!     fmean	  : first normalized absolute moment
!     variance: second normalized central moment
!     skew    : skewness coefficient
!     fkurt   : kurtosis coefficient
!
!******************************************************************************
!
! Skewness: characterizes the degree of asymmetry of a distribution around
! its mean. Positive skewness indicates a distribution with an asymmetric 
! tail extending toward more positive values. Negative skewness indicates a
! distribution with an asymmetric tail extending toward more negative values.
!
! Kurtosis: characterizes the relative peakedness or flatness of a distribution
! compared with the normal distribution. Positive kurtosis indicates a 
! relatively peaked distribution. Negative kurtosis indicates a relatively 
! flat distribution.  
!
!******************************************************************************     
      subroutine moments_breakthru_ (name,alfa,fmom,skew,fkurt) 
	  use list_class
	  use gslib, only: open_fname
	  use global_variables, only: fdbg 
!     -------------------------------------------------------------------------
!       SPECIFICATIONS:
!     -------------------------------------------------------------------------
	  implicit none

      type(breakthru_cl)       :: name
      type(partID_cl), pointer :: part
	  real*8, intent(out)      :: alfa(4),fmom(4)
	  real*8, intent(out)      :: skew,fkurt
	  integer                  :: nd,idat
	  integer                  :: unit,i
	  real*8                   :: fmean,vari
!     -------------------------------------------------------------------------
      nd = name%np
	  
    if (nd <= 1) return  
    
    if (nd >= 1) then 
     
!-------> CALCULATES ABSOLUTE NON-CENTRAL MOMENTS

      fmean      = 0.0d0
      vari       = 0.0d0
	  alfa       = 0.0d0
	  fmom       = 0.0d0

      part => name%head

	  do idat=1,nd

         alfa(1)= alfa(1) + part%tp
         alfa(2)= alfa(2) + part%tp*part%tp
	     alfa(3)= alfa(3) + part%tp*part%tp*part%tp
	     alfa(4)= alfa(4) + part%tp*part%tp*part%tp*part%tp

         part => part%next

      end do

      alfa(1) = alfa(1) / dfloat(nd)
	  alfa(2) = alfa(2) / dfloat(nd)
	  alfa(3) = alfa(3) / dfloat(nd)
	  alfa(4) = alfa(4) / dfloat(nd)

!-------> NORMALIZED ABSOLUTE CENTRAL MOMENTS

	  fmom(1) = alfa(1)
	  fmom(2) = alfa(2) -        alfa(1) * alfa(1)
	  fmom(3) = alfa(3) - 3.d0 * alfa(2) * alfa(1) + 			      &
                          2.d0 * alfa(1) * alfa(1) * alfa(1)
	  fmom(4) = alfa(4) - 4.d0 * alfa(3) * alfa(1) + 				  &
	                      6.d0 * alfa(2) * alfa(1) * alfa(1) -		  &
                          3.d0 * alfa(1) * alfa(1) * alfa(1) * alfa(1)

	  fmean     = fmom(1)
      vari      = fmom(2)
	  if (vari.gt.0.) then
	         skew  = fmom(3) / ((dsqrt(vari))**3.d0)
	         fkurt = fmom(4) / ((dsqrt(vari))**4.d0) - 3.0d0
	  else
	         skew  = 0.d0
	         fkurt = 0.d0
	         print *, '>> variance zero in moments of breakthru'
	  end if

    end if

    end subroutine


   subroutine print_pdf_breakthru_Derivative (name,histo,NUE,fname)
	    use global_variables, only: UNEST
		use gslib, only: open_fname
		use histogram_class
		implicit none
		type(breakthru_cl),         intent(in)     :: name
		type(histo_cl),             intent(inout)  :: histo
		integer,                    intent(in)     :: NUE
		character(len=*), optional, intent(in)     :: fname
	    real*8, allocatable                        :: pdf(:),tgrid(:),Dpdf(:)
		real*8                                     :: Deriv
	    integer                                    :: i,iunit
		   allocate (pdf(histo%ngrid),tgrid(histo%ngrid))
           call open_fname(fname,iunit) 
		   if (NUE==0) then 
		     call assign_Ndata_To_histo_PDF_ (histo,name%np)
			 call assign_NUE_To_histo_PDF_ (histo,NUE) 
			 call set_original_histo_values_ (histo)    
		     call histo_breakthru_ (name,histo%ngrid,histo%bwBTC,histo%dt,histo%tminBTC,histo%tmaxBTC,histo%Kernel,tgrid,pdf,NUE)
		     do i=1,histo%ngrid
                write(iunit,'(2(g21.14,x))') tgrid(i),pdf(i)
		     end do
           else
		     stop 'Program does not consider NUE>=1'
		   end if
           close(iunit)
           deallocate (pdf,tgrid)
	end subroutine


!***********************************************************************
!   SUBROUTINE TO CALCULATE CDF OF ARRIVAL TIMES
!***********************************************************************
    subroutine cdf_breakthru_ (name,rvr,rcdf) 
      use constants, only: pi
	  use gslib, only: sortem
	  use list_class
	  implicit none
      type(breakthru_cl), intent(in)      :: name
	  type(partID_cl), pointer            :: part
	  real*8,             intent(out)     :: rcdf(:),rvr(:)
	  real*8, allocatable                 :: mp(:)
	  integer                             :: i,j,nd
	  real*8                              :: b(1),oldcp,cp,qq

      nd = name%np
      if (nd <= 1)  return

      allocate (mp(nd))
      
      part => name%head
      
	  do i=1,nd
         rvr(i) = part%tp
         mp(i)  = part%mp
         part => part%next
      end do      
      

!..... sort arrival times

      b(1) = 1.d0

      call sortem(1,nd,rvr,1,mp,b,b,b,b,b,b)

!...... calculates cdf

	  rcdf  = 0.0d0
      oldcp = 0.0d0
	  cp    = 0.0d0

	do i=1,nd
	   cp      = cp + mp(i)
	   rcdf(i) = (cp + oldcp) * 0.5d0
	   oldcp   = cp
	end do
 
    deallocate(mp)
 
    end subroutine



!********************************************************************************************
!         SMOOTH A HISTOGRAM WITH KERNEL DENSITY FUNCTIONS
!********************************************************************************************
!
! Creates a smooth distribution with kernel distribution functions
! Either a uniform or a gaussian kernel can be specified
!
!     alfa    : vector of normalized absolute moments
!     fmom    : vector of normalized central moments
!     rvr     : vector of data values 
!     nd      : number of data values (particles arrived control surface)
!     ngrid   : number of grid values at which we evaluate the histogram
!     bw     : bandwith of kernel density function
!     flagker : kernel flag. When 1 = gaussian kernel density function
!                            When 0 = uniform kernel density function
!     flagdt  : flag for bandwith of kernel density function
!               0 = bw is specified by the user
!               1 = bw estimated as the optimal when distribution is gaussian
!     flaghist: 0 = only moment calculation  1 = moment calc. + histogram
!     fmean	: first normalized absolute moment
!     variance: second normalized central moment
!     skew    : skewness coefficient
!     fkurt   : kurtosis coefficient
!     dt      : grid interval, defines grid values at which the histogram is evaluated 
!     bw      : bandwidth of kernel density function
!
!********************************************************************************************
   subroutine histo_breakthru_ (name,ngrid,bw,dt,t1,t2,KernelType,tgrid,pdf,NUE) 
	  use gslib, only: open_fname
	  use global_variables, only: fdbg,UNEST
	  use library_histo
	  use gslib, only: sortem
	  use constants, only: pi
	  use list_class

!       SPECIFICATIONS:
!     ---------------------------------------------------------------------------------------
	  implicit none

      type(breakthru_cl), intent(in)      :: name
      integer,             intent(in)     :: ngrid
	  character(len=*),    intent(in)     :: KernelType
	  integer,             intent(in)     :: NUE   
	  real*8,              intent(inout)  :: t1,t2
	  real*8,              intent(inout)  :: bw,dt
	  real*8,              intent(inout)  :: pdf(:),tgrid(:)     
	  integer                             :: i,j,nd
	  integer                             :: unit
	  real*8                              :: tj,ti,dij,weight
	  real*8, allocatable                 :: tp(:),mp(:)
	  character(len=len(trim(adjustl(KernelType)))) :: Kernel
	  logical, save                       :: Calc_limit_1 = .FALSE.
	  logical, save                       :: Calc_limit_2 = .FALSE.
	  real*8, allocatable                 :: BO(:),G(:),W(:)
	  integer                             :: NKE,KORD,ISMO,NBO,NY
	  real*8                              :: ALPHA 
	  real*8                              :: b(1)
	  type(partID_cl), pointer            :: part
	  real*8                              :: tmin,tmax

!     ----------------------------------------------------------------------------------------
!     NOTE: bw is the bandwidth which is defined as half the width of the bin
!     ----------------------------------------------------------------------------------------
      
      Kernel = KernelType
      
	  if (NUE == 1) Kernel = 'DENEST'
	  
	  nd = name%np
	  
	  pdf   = 0.d0
      tgrid = 0.d0	 
	  
	  if (nd <= 5) then
	      call open_fname (fdbg,unit)
		  write(unit,*)
		  write(unit,*) '**** WARNING: too few values in BTCs - NOT EVALUATED'
		  write(unit,*)
		  write(*,*)    '**** WARNING: too few values in BTCs - NOT EVALUATED'
		  return
	  end if

     allocate (tp(nd),mp(nd))

     part => name%head
     
     tmin=UNEST
     tmax=UNEST
     do i=1,nd
	  tp(i) = part%tp
	  mp(i) = part%mp
      part => part%next
      if (tmin==UNEST .or. tmin>tp(i)) tmin=tp(i)
      if (tmax==UNEST .or. tmax<tp(i)) tmax=tp(i)
     end do


   select case (Kernel)

   case ('PLUGIN')

      ! calculates limits of the histogram

	  bw = optimum_bandwidth (tp,Kernel,ngrid)

      !NOTE: this bandwidth is only used to estimate the limits

      if (t1 <= 0.d0 .or. Calc_limit_1) then
	              Calc_limit_1 = .TRUE.
	              t1 = tmin - 2.0d0 * bw
				  if(t1 < 0.d0) t1=tmin - bw
				  if(t1 < 0.d0) t1=tmin
	  end if

      if (t2 <= 0.d0 .or. Calc_limit_2) then
	              Calc_limit_2 = .TRUE.
	              t2 = tmax + 2.0d0 * bw
	  endif   

      dt=(t2-t1)/dfloat(ngrid+1) 

      ! calculate the histogram
	
      do j=1,ngrid
	      tgrid(j) = t1 + (j-1) * dt
	  end do         

      b(1) =1.d0

      call sortem(1,nd,tp,0,b,b,b,b,b,b,b) 
	  
	  if(tp(1)==tp(nd)) then
	     write(*,*) ' *** BTC Not evaluated: STEP FUNCTION ***'
		 return
	  end if 

      call PLUGIN(tp,nd,tgrid,ngrid,pdf,bw,mp)

      deallocate (tp,mp)

   case ('DENEST')

      !CAUTION: When using DENEST all particles should have equal mass

	  allocate (BO(2),G(ngrid),W(0:nd)) 

      ! calculates limits of the histogram

	  bw = optimum_bandwidth (tp,Kernel,ngrid)

      !NOTE: this bandwidth is only used to estimate the limits

      if (t1 <= 0.d0 .or. Calc_limit_1) then
	              Calc_limit_1 = .TRUE.
	              t1 = tmin - 2.0d0 * bw
				  if(t1 < 0.d0) t1=tmin - bw
				  if(t1 < 0.d0) t1=tmin
	  end if

      if (t2 <= 0.d0 .or. Calc_limit_2) then
	              Calc_limit_2 = .TRUE.
	              t2 = tmax + 2.0d0 * bw
	  endif   

      dt=(t2-t1)/dfloat(ngrid+1) 

      ! calculate the histogram
	
      do j=1,ngrid
	      tgrid(j) = t1 + (j-1) * dt
	  end do         

      call sortem(1,nd,tp,0,b,b,b,b,b,b,b) 

	  if(tp(1)==tp(nd)) then
	     write(*,*) ' *** BTC Not evaluated: STEP FUNCTION ***'
		 return
	  end if 


      NKE   = 0      !TYPE OF KERNEL (1: MIN. NORMAL,2: OPTIMAL, 0: NORMAL)    
	  !NUE   = 0      !ORDER OF DERIVATIVE (0-2)
      KORD  = NUE+2  !ORDER OF KERNEL 
      ISMO  = 0      !0: OPTIMIZATION OF BANDWIDTH,ELSE NO OPTIMIZATION (SMOOTHING WITH B)       
      NBO   = 0      !TREATMENT OF BOUNDARY
      NY    = 0      !0: FIX BANDWIDTH, 1: VARIABLE BANDWIDTH (B(I)=B*G(I), I=1,...,M) 
      BO(1) = 0.     !BO(1) LEFT BOUNDARY,  BO(1)<=X(1),  (DUMMY FOR NBO=0)     
      BO(2) = 0.     !BO(2) RIGHT BOUNDARY, BO(2)>=X(N)             
      G     = 0.     !G(M): BANDW. FACTORS (FOR NY=1 AND ALPHA=0.) G(M) IS DUMMY FOR NY=0         
      ALPHA = 0.     !~=0,>=-1.,<=1.: G IS PROPORT. DENSITY**-ALPHA (COMPUTED IN THE PROGRAM)                                      
      W     = 0.     !W(0:N) WORK ARRAY                                     

	  call HOPTDE(tp,nd,NKE,NUE,KORD,NBO,NY,ISMO,BO,tgrid,ngrid,G,ALPHA,W,pdf,bw) 
	  
      deallocate (tp,mp)
	  deallocate (BO,G,W)
   
   case default  !-------> HISTOGRAM: KERNEL DENSITY ESTIMATION

      ! calculates optimal bandwith of kernel density for a gaussian distribution
      
	  if (bw <= 0.d0) bw = optimum_bandwidth (tp,Kernel,ngrid)

      if (bw <= 0.d0) then
            print *, ' Warning: Could not calculate bandwidth breakthrough curve'
			return
	  end if
      
	  ! calculates limits of the hsitogram

      if (t1 <= 0.d0 ) then
	              t1 = tmin - 2.0d0 * bw
				  if(t1 < 0.d0) t1=tmin - bw
				  if(t1 < 0.d0) t1=tmin
	  end if

      if (t2 <= 0.d0 ) then
	              t2 = tmax + 2.0d0 * bw
	  endif   

      dt=(t2-t1)/dfloat(ngrid+1) 

      ! calculate the histogram
	
      do j=1,ngrid
            
			 tgrid(j) = t1 + (j-1) * dt
	   
	     do i=1,nd

             dij=(dabs(tp(i)-tgrid(j)))/bw
             
			 !*********************************************************
			 ! Kernel density functions as defined by Bagtzoglou (1992)
			 ! Numerical methods for Part. Diff. Eq., 8, 325-240
			 !*********************************************************

		     select case (Kernel)  !select kernel density functions:
                case ('BOX') !Uniform distribution
				      if (dij <= 1.0d0) then
					         weight = 0.5d0
	                         pdf(j) = pdf(j) + weight * ( mp(i) / bw )  
					  end if               
	            case ('TRIANGLE') !Triangle distribution
				      if (dij <= 1.0d0) then
					         weight=1.d0-dij
	                         pdf(j) = pdf(j) + weight * ( mp(i) / bw )  
					  end if		             
		        case ('GAUSS') !Gaussian distribution
				    if (dij <= 1.5d0) then
					         weight= dexp(-pi*dij*dij)    
							 pdf(j) = pdf(j) + weight * ( mp(i) / bw )
					end if 
                case default   !Uniform distribution
				    if (dij <= 0.5d0) then
					         weight=1.0d0                
	                         pdf(j) = pdf(j) + weight * ( mp(i) / bw )  
					end if
		    end select

	   end do
	end do
 
  end select

 end subroutine



    function optimum_bandwidth (rvr,KernelType,ngrid) result (bw) !only used for calculating the histogram 
      use gslib, only: sortem, locate, powint
	  use constants, only: pi
      implicit none
	  integer,                          intent(in)  :: ngrid
	  real*8,                           intent(in)  :: rvr(:)
	  character(len=*),                 intent(in)  :: KernelType
	  character(len=len(trim(adjustl(KernelType)))) :: Kernel
  	  real*8, allocatable                           :: rcdf(:)
	  real*8                                        :: b(1),oldcp,cp,qq,xrange,glq,guq
	  real*8                                        :: bw  
	  integer                                       :: nd,i,j 

         nd = size(rvr)

         allocate (rcdf(nd)); rcdf  = 0.d0

         select case (Kernel)
		 
		   case ( 'BOX');  bw = (rvr(1)-rvr(nd))/dfloat(ngrid)  
		   
		   case default

              b(1) =1.d0

              call sortem(1,nd,rvr,0,b,b,b,b,b,b,b) !  SORT ARRIVAL TIMES ARRAY
              
			  !.............. CALCULATES CDF

              oldcp = 0.0d0  
	          cp    = 0.0d0
	
	          do i=1,nd
	             cp      = cp + 1.0d0 / dfloat(nd)
	             rcdf(i) = (cp + oldcp) * 0.5d0
	             oldcp   = cp
	          end do

              !................CALCULATES MEDIAN AND QUARTILES

  			  qq  = 0.25d0
              call locate(rcdf,nd,1,nd,qq,j)
              if (j>0) then; glq = powint(rcdf(j),rcdf(j+1),rvr(j),rvr(j+1),qq,1.d0)
                       else; glq = rvr(1)
              endif

              qq  = 0.75d0
              call locate(rcdf,nd,1,nd,qq,j)
              if (j<nd) then; guq = powint(rcdf(j),rcdf(j+1),rvr(j),rvr(j+1),qq,1.d0)
                        else; guq = rvr(nd)
			  endif

	          xrange=(guq-glq)/1.348979500d0
	          bw=1.06d0*xrange*(dfloat(nd)**(-0.2d0))
			  bw=bw*dsqrt(2.d0*pi)
	     
		 end select

    deallocate (rcdf)

	end function



    end module

