
!*****************************************************************************************************
!   PLANE CLASS
!*****************************************************************************************************   
    module plane_class
	use breakthru_class
	implicit none

	private
	public :: plane_cl     ! class
	public ::                               & ! methods
              assign_plane_               , & 
              print_plane_                , &
			  print_plane_dispersivity_   , &
			  print_plane_moments_        , &
			  read_assign_plane_          , &
			  print_plane_breakthru_        

	type plane_cl
		 real*8                         :: xplane=0.d0,yplane=0.d0,zplane=0.d0      !Position of a plane
		 real*8                         :: A=0.d0,B=0.d0,C=0.d0,D=0.d0              !Defines the plane as A*x+B*y+C*z+D=0
		 character(len=2)               :: tipo ='NA'                               !'xx'=parallel x direction,'yy'=parallel to y direction
		 integer                        :: out = 0                                  !True = particle is removed when crossing 
		 type(breakthru_cl), pointer    :: btc(:)                                   !Arrival time
	end type plane_cl

    interface assign_plane_
	     module procedure assign_plane_xxORyyORzz, assign_plane_general; end interface

    interface print_plane_moments_
	     module procedure print_plane_moments_,print_plane_moments_1_; end interface

    interface print_plane_dispersivity_
	     module procedure print_plane_dispersivity_,print_plane_dispersivity_1_; end interface

    interface print_plane_
	     module procedure print_plane_1_,print_plane_2_; end interface


	contains
   
	

	subroutine assign_plane_xxORyyORzz (name,a,typ,flag,curve)      !old plume constructor for planes parallel and perpendicular to x,y
	     use breakthru_class
		 use gslib, only: upper_case
		 implicit none
		 type(plane_cl),                intent(inout) :: name
		 real*8,                        intent(in)    :: a
		 character(len=*),    optional, intent(inout) :: typ
	     integer,             optional, intent(in)    :: flag
		 type(breakthru_cl),  optional, intent(in)    :: curve
            
            name%A = 0.d0
			name%B = 0.d0
			name%C = 0.d0
			name%D = 0.d0

            if(present(typ)) typ = upper_case(typ)

			if (.not.present(typ)) then
			     name%tipo = 'XX'
			else
			     name%tipo = typ; endif

			if (name % tipo == 'XX') then 
			       name % xplane = 0.d0
				   name % yplane = a
			else if (name % tipo == 'YY') then
			       name % xplane = a
				   name % yplane = 0.d0
			else if (name % tipo == 'ZZ') then
			       name % xplane = 0.d0
				   name % yplane = 0.d0
                   name % zplane = a
			end if
			if(.not.present(flag)) then
			    name%out = 0
			else
			    name%out = flag; endif
			if (.not.present(curve)) then
			    return
			else
			    name%btc = curve
			end if
	end subroutine  

	subroutine assign_plane_general (name,A,B,C,D,flag,flagchar,curve)      !plume constructor for planes defined as Ax+By+Cz+D=0
	     use breakthru_class
		 use gslib, only: upper_case
		 implicit none
		 type(plane_cl),                intent(inout) :: name
		 real*8,                        intent(in)    :: A,B,C,D
	     integer,             optional, intent(in)    :: flag
		 character(len=3),    optional, intent(inout) :: flagchar
		 type(breakthru_cl),  optional, intent(in)    :: curve	          
			       name%xplane = 0.d0
				   name%yplane = 0.d0
				   name%tipo = 'GL' 
			       name%A = A
				   name%B = B
				   name%C = C
				   name%D = D
			if(.not.present(flag)) then
			    name%out = 0  ! the default is that particles pass through the plane and are not excluded
			else
			    name%out = flag; endif
			if(present(flagchar)) flagchar = upper_case(flagchar)
			if(present(flagchar)) then  !flagchar rules and overwrites the previous flag option
                if (flagchar == 'OUT') name%out = 1
				if (flagchar == 'PAS') name%out = 0
			end if
			if (.not.present(curve)) then
			    return
			else
			    name%btc = curve
			end if
	end subroutine  

    subroutine read_assign_plane_ (name,iunit)
	     use gslib, only: ifcharacter
		 implicit none
		 type(plane_cl),  intent(inout) :: name
		 integer,         intent(in)    :: iunit
		 logical                        :: exists,oldinput
		 real*8                         :: dist,A,B,C,D
		 character(len=2)               :: tipo
		 character(len=3)               :: flagchar
		 integer                        :: flag
              inquire (unit=iunit,exist=exists)
			  if(.not.exists) then
			     print *, 'file does not exist during read_assign_plane'
				 stop
			  end if
	          read(iunit,*) dist,tipo
		      oldinput = ifcharacter (tipo)
		      backspace(iunit)
		      if (oldinput) then
	              read(10,*) dist,tipo,flag;  call assign_plane_xxORyyORzz (name,dist,tipo,flag)
		      elseif (.not.oldinput) then
		          read(10,*) A,B,C,D,flagchar; call assign_plane_general (name,A,B,C,D,flag,flagchar)
		      end if            
	end subroutine


    subroutine print_plane_1_ (name,fname)
         use gslib, only: open_fname, upper_case
		 implicit none
		 type(plane_cl), intent(inout)         :: name
		 character(len=*),optional, intent(in) :: fname
		 integer                               :: iunit          
		   call open_fname(fname,iunit)
		   name%tipo = upper_case(name%tipo)
		   if(name%tipo == 'XX') then
                write(iunit,1) 'ZONE T= "Yplane: ',name%yplane,' Type: ',name%tipo,' OUT: ',name%out, '"'
		   else if (name%tipo == 'YY') then
                write(iunit,1) 'ZONE T= "Xplane: ',name%xplane,' Type: ',name%tipo,' OUT: ',name%out, '"'
		   else if (name%tipo == 'ZZ') then
                write(iunit,1) 'ZONE T= "Zplane: ',name%zplane,' Type: ',name%tipo,' OUT: ',name%out, '"'
		   else if (name%tipo == 'GL') then
                write(iunit,2) 'ZONE T= "Ax+By+Cz+D=0 (A,B,C,D): ',name%A,name%B,name%C,name%D,' Type: ',name%tipo,' OUT: ',name%out,'"'
		   end if
		 2 format(a33,4(g15.6,x),a7,a2,x,a6,i2,a2)
         1 format((a17,g15.6,5x),a7,a2,x,a6,i2,a2)                                                       
		   close(iunit)
	end subroutine

    subroutine print_plane_2_ (name,fname,ispe)
         use gslib, only: open_fname, upper_case
		 implicit none
		 type(plane_cl), intent(inout)         :: name
		 integer,        intent(in)            :: ispe
		 character(len=*),optional, intent(in) :: fname
		 integer                               :: iunit          
		   call open_fname(fname,iunit)
		   name%tipo = upper_case(name%tipo)
		   if(name%tipo == 'XX') then
                write(iunit,1) 'ZONE T= "Yplane: ',name%yplane,' Specie: ',ispe, '"'
		   else if (name%tipo == 'YY') then
                write(iunit,1) 'ZONE T= "Xplane: ',name%xplane,' Specie: ',ispe, '"'
		   else if (name%tipo == 'ZZ') then
                write(iunit,1) 'ZONE T= "Zplane: ',name%zplane,' Specie: ',ispe, '"'
		   else if (name%tipo == 'GL') then
                write(iunit,2) 'ZONE T= "General): ',name%A,name%B,name%C,name%D,' Specie: ',ispe, '"'
		   end if
		 2 format(a19,4(f7.1,x),a9,i2,a2)
         1 format((a17,f7.1,x),a9,i2,a2)                                                       
		   close(iunit)
	end subroutine

    subroutine print_plane_moments_ (name,fname)
         use breakthru_class
         use gslib, only: open_fname, upper_case
         use global_variables, only: nspecie,UNEST
		 implicit none
		 type(plane_cl),            intent(inout) :: name
         character(len=*),optional, intent(in)    :: fname
		 real*8                                   :: mean(3),var(3)
		 real*8                                   :: alfa(4),mom(4)
	     real*8                                   :: skew,kurt,dis
         integer                                  :: iunit,np,ispe

		 do ispe=1,nspecie
		    if (name%btc(ispe)%np>5) then
            call open_fname(fname,iunit)
		    call moments_breakthru_ (name%btc(ispe),alfa,mom,skew,kurt)
		    np = name%btc(ispe)%np
		    name%tipo = upper_case(name%tipo)
            if(name%tipo == 'YY') then
		      dis = name%xplane
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','xplane= ',dis      
	        else if(name%tipo == 'XX') then
		      dis = name%yplane
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','yplane= ',dis
	        else if(name%tipo == 'ZZ') then
		      dis = name%zplane
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','zplane= ',dis
			else if(name%tipo == 'GL') then			  
			  dis = UNEST
			  write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','distan= ',dis
            end if
			end if
		  end do

	 1    format (9(g20.12,1x),i20,x,a6,a8,g15.6)
     end subroutine

    subroutine print_plane_moments_1_ (name,btc,fname)
         use breakthru_class
         use gslib, only: open_fname, upper_case
         use global_variables, only: nspecie,UNEST
		 implicit none
		 type(plane_cl),            intent(inout) :: name
		 type(breakthru_cl),        intent(in)    :: btc(nspecie)
         character(len=*),optional, intent(in)    :: fname
		 real*8                                   :: mean(3),var(3)
		 real*8                                   :: alfa(4),mom(4)
	     real*8                                   :: skew,kurt,dis
         integer                                  :: iunit,np,ispe
		   if (btc(1)%np>5) then
           call open_fname(fname,iunit)
		   do ispe=1,nspecie 
		    call moments_breakthru_ (btc(ispe),alfa,mom,skew,kurt)
		    np = btc(ispe)%np
		    name%tipo = upper_case(name%tipo)
            if(name%tipo == 'YY') then
		      dis = name%xplane
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','xplane= ',dis      
	        else if(name%tipo == 'XX') then
		      dis = name%yplane
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','yplane= ',dis
	        else if(name%tipo == 'ZZ') then
		      dis = name%zplane
		      write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','zplane= ',dis
			else if(name%tipo == 'GL') then			  
			  dis = UNEST
			  write(iunit,1) mom(1),alfa(2),skew,kurt,alfa(3),alfa(4),mom(2),mom(3),mom(4),np,'plane ','distan= ',dis
            end if
		   end do
		   end if	 
	1     format (9(g20.12,1x),i20,x,a6,a8,g15.6)
     end subroutine


   subroutine print_plane_breakthru_ (i,btc,fname)
    use gslib, only: open_fname
	use breakthru_class
	implicit none
	type(breakthru_cl),         intent(in)    :: btc
    integer,                    intent(in)    :: i
	character(len=*), optional, intent(in)    :: fname
	integer                                   :: n,iunit
	  if (btc%np<5) return
	     call open_fname(fname,iunit)
		 write(iunit,2)'zone plane="',i,'"'
	  2  format(a12,i5,a2) 
		 call print_breakthru_ (btc,fname)
   end subroutine



    subroutine print_plane_dispersivity_ (name,fname)
         use breakthru_class
         use gslib, only: open_fname, upper_case
         use global_variables, only: nspecie
		 implicit none
		 type(plane_cl),            intent(inout) :: name
         character(len=*),optional, intent(in)    :: fname
		 real*8                                   :: mean(3),var(3)
		 real*8                                   :: alfa(4),mom(4)
	     real*8                                   :: skew,kurt
		 real*8                                   :: dis,ta,a11,a22,a33
         integer                                  :: iunit,np,ispe
		   call open_fname(fname,iunit)
		   name%tipo = upper_case(name%tipo)
		 do ispe=1,nspecie
		   call calculate_variance_particle_position_from_btc_ (name%btc(ispe),mean,var)
		   call moments_breakthru_ (name%btc(ispe),alfa,mom,skew,kurt)
		   np = name%btc(ispe)%np
		   a11 = 0.d0
		   a22 = 0.d0
		   a33 = 0.d0
		   if(name%tipo == 'XX') then
		      dis = dabs(mean(2))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(1)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis               
		   else if (name%tipo == 'YY') then
              dis = dabs(mean(1))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(2)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis
		   else if (name%tipo == 'GL') then
		      dis = dabs(mean(1))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(1)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis               			   
		   end if
		   write(iunit,1) dis,ta,a11,a22,a33,var(1),var(2),var(3),mean(1),mean(2),mean(3),np
		  end do
	   1  format(11(g15.6,x),i7)
		   close(iunit)
	end subroutine


    subroutine print_plane_dispersivity_1_ (name,btc,fname)
         use breakthru_class
         use gslib, only: open_fname, upper_case
		 implicit none
		 type(plane_cl),            intent(inout) :: name
		 type(breakthru_cl),        intent(in)    :: btc
         character(len=*),optional, intent(in)    :: fname
		 real*8                                   :: mean(3),var(3)
		 real*8                                   :: alfa(4),mom(4)
	     real*8                                   :: skew,kurt
		 real*8                                   :: dis,ta,a11,a22,a33
         integer                                  :: iunit,np
		   if (btc%np<5) return
		   call open_fname(fname,iunit)
		   name%tipo = upper_case(name%tipo)
		   call calculate_variance_particle_position_from_btc_ (btc,mean,var)
		   call moments_breakthru_ (btc,alfa,mom,skew,kurt)
		   np = btc%np
		   a11 = 0.d0
		   a22 = 0.d0
		   a33 = 0.d0
		   if(name%tipo == 'XX') then
		      dis = dabs(mean(2))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(1)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis               
		   else if (name%tipo == 'YY') then
              dis = dabs(mean(1))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(2)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis
           else if (name%tipo == 'ZZ') then
              stop 'cannot calculate ZZ plane dispersivity for the moment..'
              dis = dabs(mean(3))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(2)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis
		   else if (name%tipo == 'GL') then
		      dis = dabs(mean(1))
			  ta  = mom(1)
			  if (mom(1) > 0.d0 ) a11 = 0.5d0*dis*mom(2)/(mom(1)**2)
			  if (dis    > 0.d0 ) a22 = 0.5d0*var(1)/dis
			  if (dis    > 0.d0 ) a33 = 0.5d0*var(3)/dis               			   
		   end if
		   write(iunit,1) dis,ta,a11,a22,a33,var(1),var(2),var(3),mean(1),mean(2),mean(3),np
		1  format(11(g15.6,x),i7)
		   close(iunit)
	end subroutine

!***********************************************************************
!   CALCULATE VARIANCE OF PARTICLE DISPLACEMENT
!***********************************************************************

    subroutine calculate_variance_particle_position_from_btc_ (btc,mean,var)
         use list_class
		 implicit none
		 type(breakthru_cl), intent(in)  :: btc
         real*8,             intent(out) :: mean(3),var(3)
         type(partID_cl), pointer        :: part
         integer                         :: n
         
         if (btc%np<2) return
         
         part => btc%head
         
         if (.not.associated(part)) then
          write(*,*) 'could not estimate spatial moments of a breakthru'
          return
         end if
         
         mean = 0.d0
         var  = 0.d0
         n    = 0
         
         do 
         
         n = n + 1
         
         mean(1) = mean(1)+part%xp
         mean(2) = mean(2)+part%yp
         mean(3) = mean(3)+part%zp
 
         var(1) = var(1)+part%xp*part%xp
         var(2) = var(2)+part%yp*part%yp
         var(3) = var(3)+part%zp*part%zp
        
         part => part%next         
         
         end do

         mean(1) = mean(1)/dfloat(n)
         mean(2) = mean(2)/dfloat(n)
         mean(3) = mean(3)/dfloat(n)

         var(1) = var(1)/dfloat(n)-mean(1)*mean(1)
         var(2) = var(2)/dfloat(n)-mean(2)*mean(2)
         var(3) = var(3)/dfloat(n)-mean(3)*mean(3)
          
         

	end subroutine



    end module
!******************************************************************************
!  VECTOR OF PLANES  
!******************************************************************************
 module plane_vect_class
 use plane_class

 implicit none

 public    !everything is public (inheritance of plane_class)
 public :: plane_vect_cl                        !class
 public ::                                    & !methods
           alloc_plane_vect_                , &
		   print_plane_dispersivity_vect_   , &
		   print_plane_breakthru_vect_      , &
		   nplane_                          , & 
		   print_plane_moments_vect_      

 type plane_vect_cl
      type(plane_cl),pointer :: num(:) => null()
 end type

 

 interface print_plane_moments_vect_
      module procedure print_plane_moments_vect_,print_plane_moments_vect_1_; end interface

 contains





 subroutine alloc_plane_vect_ (name,n,nspecie)
    use plane_class
	use breakthru_class
	implicit none
    type(plane_vect_cl), intent(inout) :: name
    integer,             intent(in)    :: n
    integer,             intent(in)    :: nspecie
	integer                            :: i,ispe 
	      allocate(name%num(n))
		  do i=1,n
		    allocate (name%num(i)%btc(nspecie))
			do ispe=1,nspecie
			   call initialize_btc_ (name%num(i)%btc(ispe))
			end do
		  end do	  
 end subroutine

 subroutine print_plane_dispersivity_vect_ (name,fname)
    use plane_class
	implicit none
	type(plane_vect_cl),        intent(inout) :: name
	character(len=*), optional, intent(in)    :: fname
	integer                            :: i,n
	  if (associated(name%num)) then
	  n = size(name%num)
      do i=1,n
         call print_plane_dispersivity_ (name%num(i),fname);  end do
	  end if
 end subroutine

 subroutine print_plane_moments_vect_ (name,fname)
    use gslib, only: open_fname
	use plane_class
	use breakthru_class
	use global_variables, only: nspecie,UNEST
	implicit none
	type(plane_vect_cl),        intent(inout) :: name
	character(len=*), optional, intent(in)    :: fname
	integer                                   :: i,n,iunit
	character(len=2)                          :: tipo
	real*8                                    :: dis
	integer                                   :: ispe

	  if (associated(name%num)) then
	  n = size(name%num)
	  do i=1,n
	  do ispe=1,nspecie      
	      call print_moments_breakthru_ (name%num(i)%btc(ispe),fname)
	  end do
	  end do
      call open_fname(fname,iunit)
      write(iunit,*)
	  do i=1,n
	  do ispe=1,nspecie
	      tipo = name%num(i)%tipo
		  if      (tipo == 'YY') then; dis = name%num(i)%xplane
          else if (tipo == 'XX') then; dis = name%num(i)%yplane
          else if (tipo == 'ZZ') then; dis = name%num(i)%zplane
		  else if (tipo == 'GL') then; dis = UNEST; end if
          write(iunit,*) i,tipo,dis
	  end do  
	  end do
	  end if
 end subroutine

 subroutine print_plane_moments_vect_1_ (name,i,fname)
    use gslib, only: open_fname
	use plane_class
	use breakthru_class
	use global_variables, only: nspecie,UNEST
	implicit none
	type(plane_vect_cl),        intent(inout) :: name
	integer,                    intent(in)    :: i
	character(len=*), optional, intent(in)    :: fname
	integer                                   :: n,iunit,ispe
	character(len=2)                          :: tipo
	real*8                                    :: dis
	  if (associated(name%num)) then
	  do ispe=1,nspecie
	      call print_moments_breakthru_ (name%num(i)%btc(ispe),fname)
          call open_fname(fname,iunit)
          write(iunit,*)
	      tipo = name%num(i)%tipo
		  if      (tipo == 'YY') then; dis = name%num(i)%xplane
          else if (tipo == 'XX') then; dis = name%num(i)%yplane
          else if (tipo == 'ZZ') then; dis = name%num(i)%zplane
		  else if (tipo == 'GL') then; dis = UNEST; end if
          write(iunit,*) i,tipo,dis
      end do
	  end if
 end subroutine



 subroutine print_plane_breakthru_vect_ (name,fname)
    use gslib, only: open_fname
	use breakthru_class
	use plane_class
	use global_variables, only: nspecie
	implicit none
	type(plane_vect_cl),        intent(inout) :: name
	character(len=*), optional, intent(in)    :: fname
	integer                                   :: i,n,iunit,np,ispe
	  if (.not.associated(name%num)) return
	  n = size(name%num)
	     call open_fname(fname,iunit)
	  do i=1,n
      do ispe=1,nspecie
		 np = name%num(i)%btc(ispe)%np
		 if (np > 0) write(iunit,2)'zone plane="',i,'"'
	  2  format(a12,i5,a2) 
	     close(iunit)
		 call print_breakthru_ (name%num(i)%btc(ispe),fname)
	   end do
	   end do
 end subroutine



 function nplane_ (name) result (n)
         implicit none
		 type(plane_vect_cl), intent(in) :: name
         integer                         :: n
		   if (associated(name%num)) then
		       n = size(name%num)
		   else
		       n = 0
		   end if
 end function



 end module