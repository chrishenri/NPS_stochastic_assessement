!***********************************************************************
!  LIBRARY GEOSTATISTICAL SOFTWARE
!***********************************************************************

module gslib

	interface open_fname
	     module procedure open_fname_1,open_fname_2;	end interface


contains


!***********************************************************************
!    REMOVE REPEATED VALUES FROM A VECTOR
!***********************************************************************

      subroutine remove_repeated_values (v)
          implicit none
          real*8, allocatable   :: u(:)
          real*8, allocatable, intent(inout)  :: v(:)
          integer              :: ix,nx
          logical, allocatable :: mask(:)
          integer, allocatable :: index_vector(:)
          nx = size(v)
          allocate (u(nx))
          u = v
          deallocate (v)
          ! First, find the duplicate elements
          allocate(mask(nx))
          mask = .TRUE.
          do ix = nx,2,-1
               mask(ix) = .NOT.( any(u(:ix-1)==u(ix)) )
          end  do
          ! Make an index vector
          allocate(index_vector, source=PACK([(ix, ix=1,nx) ],mask))
          ! Now copy the unique elements of a into b
          allocate(v, source=u(index_vector))
      end subroutine


!***********************************************************************
!         LOCATE VALUE WITHIN ARRAY (LOCATE)
!***********************************************************************
!
! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
! returned to indicate that x is out of range.
!
! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
      subroutine locate (xx,n,is,ie,x,j)

      implicit real*8 (a-h,o-z)

      real*8 xx(n),x

!---- modification from initial subroutine (DANI)

      if (xx(n) == x ) then
	      j = n-1
		  return
	  end if

      if (x <= 0.d0 ) then
	      j = 1
		  return
	  end if

!--------------------------------- end modification
!
! Initialize lower and upper methods:
!
      if(is.le.0) is = 1
      jl = is-1
      ju = ie
      if(xx(n).le.x) then
            j = ie
            return
      end if
!
! If we are not done then compute a midpoint:
!
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
!
! Replace the lower or upper limit with the midpoint:
!
            if((xx(ie).gt.xx(is)).eqv.			 &
               (x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
!
! Return with the array index:
!
      j = jl
      
      end subroutine locate

!***********************************************************************
!         LOCATE VALUE WITHIN ARRAY (LOCATE)
!***********************************************************************
!
! Dani: Modified to restrict "j" to be between 1 and n
!
! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
! returned to indicate that x is out of range.
!
! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
      subroutine locate2 (xx,n,is,ie,x,j)

      implicit real*8 (a-h,o-z)

      real*8 xx(n),x

!---- modification from initial subroutine (DANI)

      if (x >= xx(n)) then
	      j = n - 1
		  return
	  end if

      if (x <= 0.d0 ) then
	      j = 1
		  return
	  end if

!--------------------------------- end modification
!
! Initialize lower and upper methods:
!
      if(is.le.0) is = 1
      jl = is-1
      ju = ie
      if(xx(n).le.x) then
            j = ie
            return
      end if
!
! If we are not done then compute a midpoint:
!
 10   if(ju-jl.gt.1) then
            jm = (ju+jl)/2
!
! Replace the lower or upper limit with the midpoint:
!
            if((xx(ie).gt.xx(is)).eqv.			 &
               (x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
      endif
!
! Return with the array index:
!
      j = jl
      
      end subroutine locate2

!***********************************************************************
!         POWER INTERPOLATION	 (POWINT)
!***********************************************************************
!
! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
!                 for a value of x and a power pow.
!
!-----------------------------------------------------------------------
      real*8 function powint(xlow,xhigh,ylow,yhigh,xval,pow)

	  implicit real*8 (a-h,o-z)

      parameter (EPSLON=1.0d-20)

           if ((xhigh-xlow).lt.EPSLON) then
              powint = (yhigh+ylow)/2.0d0
           else
              powint = ylow + (yhigh-ylow)* 				   &
                    (((xval-xlow)/(xhigh-xlow))**pow)
           end if

      
      end function powint

!**********************************************************************
! 	   QUICKERSORT SUBROUTINE	(SORTEM)
!**********************************************************************
!
! This is a subroutine for sorting a real array in ascending order. This
! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
! in collected algorithms of the ACM.
!
! The method used is that of continually splitting the array into parts
! such that all elements of one part are less than all elements of the
! other, with a third part in the middle consisting of one element.  An
! element with value t is chosen arbitrarily (here we choose the middle
! element). i and j give the lower and upper limits of the segment being
! split.  After the split a value q will have been found such that 
! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
! performs operations on the two segments (i,q-1) and (q+1,j) as follows
! The smaller segment is split and the position of the larger segment is
! stored in the lt and ut arrays.  If the segment to be split contains
! two or fewer elements, it is sorted and another segment is obtained
! from the lt and ut arrays.  When no more segments remain, the array
! is completely sorted.
!
!
! INPUT PARAMETERS:
!
!   ib,ie        start and end index of the array to be sorteda
!   a            array, a portion of which has to be sorted.
!   iperm        0 no other array is permuted.
!                1 array b is permuted according to array a
!                2 arrays b,c are permuted.
!                3 arrays b,c,d are permuted.
!                4 arrays b,c,d,e are permuted.
!                5 arrays b,c,d,e,f are permuted.
!                6 arrays b,c,d,e,f,g are permuted.
!                7 arrays b,c,d,e,f,g,h are permuted.
!               >7 no other array is permuted.
!
!   b,c,d,e,f,g,h  arrays to be permuted according to array a.
!
! OUTPUT PARAMETERS:
!
!    a      = the array, a portion of which has been sorted.
!
!    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
!
! NO EXTERNAL ROUTINES REQUIRED:
!
!-----------------------------------------------------------------------
      subroutine sortem(ib,ie,a,iperm,b,c,d,e,f,g,h)

      implicit real*8 (a-h,o-z)

      dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
!
! The dimensions for lt and ut have to be at least log (base 2) n
!
      integer   lt(64),ut(64),i,j,k,m,p,q
!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
!
! If this segment has more than two elements  we split it
!
 10   if (j-i-1) 100,90,15
!
! p is the position of an arbitrary element in the segment we choose the
! middle element. Under certain circumstances it may be advantageous
! to choose p at random.
!
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
!
! Start at the beginning of the segment, search for k such that a(k)>t
!
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
!
! Such an element has now been found now search for a q such that a(q)<t
! starting at the end of the segment.
!
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
!
! a(q) has now been found. we interchange a(q) and a(k)
!
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
!
! Update q and search for another pair to interchange:
!
      q = q-1
      go to 20
 50   q = k-1
 60   continue
!
! The upwards search has now met the downwards search:
!
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
!
! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
! store the position of the largest segment in lt and ut
!
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
!
! Update m and split the new smaller segment
!
 80   m = m+1
      go to 10
!
! We arrive here if the segment has  two elements we test to see if
! the segment is properly ordered if not, we perform an interchange
!
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
!
! If lt and ut contain more segments to be sorted repeat process:
!
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      
      end subroutine sortem

!***********************************************************************
!       SUBROUTINE TO OPEN A FILE
!***********************************************************************

   subroutine open_fname_1 (fname,iunit) !the file is positioned at the end-of-file 
      implicit none
	  character(len=*),optional,  intent(in)  :: fname
	  integer,                    intent(out) :: iunit
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
   end subroutine


   subroutine open_fname_2 (fname,iunit,estatus,acceso,formato) !the file is positioned at the end-of-file 
      implicit none
	  character(len=*),           intent(in)  :: fname
	  integer,                    intent(out) :: iunit
	  character(len=*),           intent(in)  :: estatus,acceso,formato
      logical                                 :: connected,exists

	       inquire(file=fname,exist=exists,opened=connected)
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if  	              
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100) 
	           open(iunit,file=fname,status=estatus,access=acceso,form=formato)
           end if       

   end subroutine



   subroutine open_fname_binary (fname,iunit) !the file is positioned at the end-of-file 
      implicit none
	  character(len=*), optional, intent(in)  :: fname
	  integer,                    intent(out) :: iunit
      logical                                 :: connected,exists

       if(.not.present(fname)) then
           iunit=6
       else
	       inquire(file=fname,exist=exists, opened=connected)
           if(.not.exists) then 
                write(*,*) '>> file not found: ',trim(fname)
                stop
           end if  	              
	       if(connected) then
               inquire(file=fname,number=iunit)
	       else if (.not.connected) then
		       iunit = generate_unit(100)
	           open(iunit,file=fname,status='unknown',access='append', form='binary')
           end if       
	   end if
   end subroutine

   subroutine open_fname_normal (fname,iunit) !the file is positioned at the beginning-of-file 
      implicit none
	  character(len=*), optional, intent(in)  :: fname
	  integer,                    intent(out) :: iunit
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
	           open(iunit,file=fname,status='unknown',access='sequential')
           end if       
	   end if
   end subroutine

   function generate_unit (i) result (value)
       implicit none
	   integer, intent(in) :: i
	   integer             :: value
	   logical             :: connected
          value = i
		  do
		    inquire(unit=value,opened=connected)
			if (.not.connected) exit
			value = value + 1
		  end do
   end function
   
   function inquire_number_lines (fname) result (nline)
        implicit none
		character(len=*), intent(in) :: fname
		character(len=1)             :: string 
		integer                      :: nline,unit,i
		   call open_fname_normal (fname,unit)
		   i = 0
		   do 
              read(unit,*,err=10,end=10) string
		      i = i + 1
		   end do
      10   continue
	       nline = i
   end function


!*********************************************************************
!          function to change string character to upper case
!*********************************************************************
      function  upper_case (string)  result (new_string) ! like C
!     -------------------------------------------------------------
!        Convert a string or character to upper case
!         (valid for ASCII or EBCDIC processors)
!     -------------------------------------------------------------
      implicit none
      character (len = *), intent(in) :: string     ! unknown length
      character (len = len(string))   :: new_string ! same length
      character (len = 26), parameter ::               &
          UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',  &
          lower = 'abcdefghijklmnopqrstuvwxyz'
      integer :: k    ! loop counter
      integer :: loc  ! position in alphabet
      new_string = string             ! copy everything
      do k = 1, len_trim (string)     ! to change letters
         loc = index ( lower, string (k:k) )              ! locate
      if (loc /= 0 ) new_string (k:k) = UPPER(loc:loc) ! convert
      end do ! over string characters
      end function upper_case

!***********************************************************************
!          function to check if string is a character 
!***********************************************************************
      function  ifcharacter (string)  result (logic)
      implicit none
      character (len = *), intent(in) :: string     ! unknown length
	  logical                         :: logic
      character (len = 52), parameter ::               &
          char = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
      integer :: k    ! loop counter
      integer :: loc  ! position in alphabet
                do k = 1, len_trim (string)     ! to change letters
                    loc = index ( char, string (k:k) )              ! locate
                    if (loc /= 0 ) then
	                   logic = .TRUE.
	                else
	                   logic = .FALSE.; endif
                end do 
      end function ifcharacter

!***********************************************************************
!          function to check if string is a number
!***********************************************************************
      function  ifnumber (string)  result (logic)
      implicit none
      character (len = *), intent(in) :: string     ! unknown length
	  logical                         :: logic
      character (len = 11), parameter ::               &
          char = '1234567890.'
      integer :: k    ! loop counter
      integer :: loc  ! position in alphabet
                do k = 1, len_trim (string)     ! to change letters
                    loc = index ( char, string (k:k) )              ! locate
                    if (loc /= 0 ) then
	                   logic = .TRUE.
	                else
	                   logic = .FALSE.; endif
                end do 
      end function ifnumber



!***********************************************************************
!        function to count how many words in a line
!***********************************************************************
      function count_words (string_big) result (num)
        implicit none
		character (len=*), intent(in)                  :: string_big
		character (len=len(trim(adjustl(string_big)))) :: string
		integer                                        :: num,k
		   string = trim(adjustl(string_big))
		   if (len(string) == 0) then
		      num = 0
			  return
		   end if
		   num = 1
		   do k=2,len_trim(string)
		         if (string(k-1:k-1) /= ' ' .and. string(k:k) == ' ') num = num + 1; end do  
	  end function		

!***********************************************************************
!     function to locate blanks
!***********************************************************************
      function locate_blanks (string_big) result (location)
	     implicit none
		 character (len=*), intent(in)                  :: string_big
		 character (len=len(trim(adjustl(string_big)))) :: string
		 integer                                        :: num,k,inum
		 integer, allocatable                           :: location(:)
		   string = trim(adjustl(string_big))
		   num = count_words (string)
		   if (num > 0) allocate (location(1+num))
		   location(1) = 1
		   location(1+num) = len(string)
		   inum = 1
		   do k=2,len_trim(string)
		         if (string(k-1:k-1) /= ' ' .and. string(k:k) == ' ') then
				 inum=inum+1
				 location(inum)=k; endif
		   end do  
      end function

!***********************************************************************
!        function to separate string in words
!***********************************************************************
      function separate_words (string_big) result(words)
	    implicit none
		character (len=*), intent(in)                  :: string_big
		character (len=len(trim(adjustl(string_big)))) :: string
		character (len=100), allocatable               :: words(:)
		integer                                        :: nwords,iword
		integer, allocatable                           :: locblanks(:)
		   string = trim(adjustl(string_big))
		   nwords = count_words (string)
		   allocate (words(nwords))
		   words = ' '
           allocate(locblanks(1+nwords))
		   locblanks = locate_blanks(string)
		   do iword=1,nwords
		      words(iword)=string(locblanks(iword):locblanks(iword+1))
		   end do
	  end function		

!***********************************************************************
!        logical function to ask whether a word is an integer 
!***********************************************************************
      function is_this_char_an_integer (string_big) result (logic)
	     implicit none
		 character (len=*), intent(in)                  :: string_big
		 character (len=len(trim(adjustl(string_big)))) :: string
		 logical                                        :: logic
		 character (len=10)                             :: num
		 integer                                        :: loc,k
		 num = '0123456789'
         logic = .TRUE.
         loc=0
		 string = trim(adjustl(string_big))
		 do k=1,len_trim(string)
		       loc = index(num,string(k:k))
			   if (loc >= 1 .and. loc <=10) then
			       cycle
			   else
			       logic = .FALSE.
				   exit
			   end if
		 end do		     
	  end function

!************************************************************************
!        function to transform character to integer
!************************************************************************
      function char_to_integer (string) result (int)
	      implicit none
		  character(len=*), intent(in) :: string
          integer                      :: unit,int
		  unit = generate_unit(183)
		  open(unit,status='scratch')
		  write(unit,*) string
		  rewind(unit)
		  read(unit,*) int
	  end function

!************************************************************************
!        function to transform integer to character
!************************************************************************
      function integer_to_char (int) result (string)
	      implicit none
          integer, intent(in)  :: int
		  character(len=10)    :: string
		  integer              :: unit
		  unit = generate_unit(183)
		  open(unit,status='scratch')
		  write(unit,*) int
		  rewind(unit)
		  read(unit,*) string
		  close(unit)
	  end function



!************************************************************************
!        function to transform character to double real
!************************************************************************
      function char_to_dble (string) result (dble)
	      implicit none
		  character(len=*), intent(in) :: string
          integer                      :: unit
		  real*8                       :: dble
		  unit = generate_unit(183)
		  open(unit,status='scratch')
		  write(unit,*) string
		  rewind(unit)
		  read(unit,*) dble
	  end function

!***********************************************************************
!   function to generate names
!***********************************************************************
   subroutine generate_name (ir,name,long_name,string)
         implicit none
			integer,          intent(in)    :: ir
            character(len=*), intent(in)    :: name
            character(len=*), intent(inout) :: long_name
            character(len=*), optional, intent(in) :: string
            character(len=len(name))        :: type_name,short_name
            character(len=20)               :: num1
            integer                         :: iunit,n,loc
            short_name = name
            n = len_trim (name)    
            loc=index(name,'.')         
            if (loc /= 0 .and. loc < n ) then
                 type_name  = name(loc:n)
		         short_name = ' '
		         short_name = name(1:loc-1)
            end if
            iunit = generate_unit(40)
            open(iunit,status='scratch')
            write(iunit,'(i20)') ir
            rewind(iunit) 
            read(iunit,'(a20)') num1
            close(iunit) 
            long_name = ' '
			if (present(string)) then
			long_name = trim(adjustl(short_name))//'_'//string//'_' &
			               //trim(adjustl(num1))//     &
						   trim(adjustl(type_name))
			else
			long_name = trim(adjustl(short_name))//'_' &
			               //trim(adjustl(num1))//     &
						   trim(adjustl(type_name))
			end if
									  
   end subroutine


!***********************************************************************        
!  function for transforming coordinates by translation and rotation
!***********************************************************************
      function rotate_coord (x,l,m,n,x0) result (xp)
	      implicit none
		  real*8                       :: xp(3)
		  real*8,           intent(in) :: x(3),l(3),m(3),n(3) ! where l,m,n are cosenos directores
		  real*8, optional, intent(in) :: x0(3)
		  real*8                       :: dx(3)
		      dx = x
		      if (present(x0)) dx = x - x0
              xp(1) = l(1) * dx(1) + l(2) * dx(2) + l(3) * dx(3)
              xp(2) = m(1) * dx(1) + m(2) * dx(2) + m(3) * dx(3)		   
		      xp(3) = n(1) * dx(1) + n(2) * dx(2) + n(3) * dx(3)
	  end function rotate_coord

!***********************************************************************        
!  vectorial product of two vectors
!***********************************************************************
      function vectorial_product (a,b) result (v)
	      implicit none
		  real*8             :: v(3)
		  real*8, intent(in) :: a(3),b(3)
             v(1) = a(2)*b(3) - a(3)*b(2)
             v(2) = a(3)*b(1) - a(1)*b(3)
			 v(3) = a(1)*b(2) - a(2)*b(1)
	  end function vectorial_product

!***********************************************************************        
!  FUNCTION FOR GENERATING UNIFORMLY DISTRIBUTED RANDOM NUMBERS
!***********************************************************************

      function rand(i) result (value)
      implicit none
	  real*8                   :: value
      integer*4, intent(inout) :: i
      i=i*65539
      if(i < 0) i=i+2147483647+1
      value=dble(i)*0.465661287525d-9
      end function rand

!***********************************************************************        
!  FUNCTION FOR GENERATING UNIFORMLY DISTRIBUTED RANDOM NUMBERS
!***********************************************************************

      function rand2(i) result (value)
      implicit none
	  real*8                   :: value
      integer*4, intent(inout) :: i
      i=i*65539
      if(i < 0) i=i+2147483647+1
      value=dble(i)*0.465661287525d-9
      end function rand2

!***********************************************************************        
!  FUNCTION FOR GENERATING UNIFORMLY DISTRIBUTED RANDOM NUMBERS
!***********************************************************************

      function rand3(i) result (value)
      implicit none
	  real*8                   :: value
      integer*4, intent(inout) :: i
      i=i*65539
      if(i < 0) i=i+2147483647+1
      value=dble(i)*0.465661287525d-9
      end function rand3

!***********************************************************************        
!  function for generating standard normally distributed random numbers
!***********************************************************************
      function random_normal () result (value)
	    implicit none
		real   :: harvest(12)
		real*8 :: value
		integer :: seed(2)
		!call random_seed(get=seed)
		call random_number (harvest)
		value = dble(sum(harvest))-6.d0
		!print *, seed(1),seed(2),value
	  end function

!**********************************************************************
!  function to check if a real number is an integer 
!**********************************************************************
      function if_integer (a) result (logic)
         implicit none
		 real, intent(in)   :: a
		 logical            :: logic
         integer            :: i
		 real               :: residu
		    i = int(a) 
			residu = float(i) - a
			if (residu /= 0.) logic = .FALSE.
			if (residu == 0.) logic = .TRUE.                 
	  end function

!**********************************************************************
!         function to solve a*x**2 + b*x + c = 0
!**********************************************************************
      real*8 function root2 (a,b,c)
	  implicit real*8 (a-h,o-z)
  
      if (a == 0.d0 .and. b == 0.d0) stop '*** error: root2 ***'
  
      if (a == 0.d0) then
	    x1 = -c/b
	    x2 = x1
	  else if (a /= 0.d0) then 
		disc=b*b-4.d0*a*c
      	       if (disc.eq.0.d0) then
	               x1 = -0.5d0*b/a
	               x2 = x1
	           else if (disc.lt.0.d0) then
	               stop '>> error: root2 not real ***'
	           else if (disc.gt.0.d0) then
	               x1 = ( -b + dsqrt(disc) ) / ( 2.d0*a )
	               x2 = ( -b - dsqrt(disc) ) / ( 2.d0*a )
	           end if
	  end if
  
	  if (x1.lt.0.d0.and.x2.lt.0.d0) then
	      print *, '>>error: two negative solutions in root2',a,b,c
          stop
	  end if
  
	  if (x1.ge.0.d0.and.x2.ge.0.d0.and.x1.ne.x2) then
	    root2 = min(x1,x2)
	    return
	  end if
  
	  if (x1.ge.0.d0) root2 = x1
	  if (x2.ge.0.d0) root2 = x2
  
	  end function root2


!**********************************************************************
!         function to solve a*x**2 + b*x + c = 0
!**********************************************************************
      function root (a,b,c) result (value)
	  implicit real*8 (a-h,o-z)
	  real*8 :: value(2)

      if (a == 0.d0 .and. b == 0.d0) then
	      value=0.d0
		  return; endif

      if (a == 0.d0) then
	    x1 = -c/b
	    x2 = x1
	  else if (a /= 0.d0) then 
		disc=b*b-4.d0*a*c
      	       if (disc.eq.0.d0) then
	               x1 = -0.5d0*b/a
	               x2 = x1
	           else if (disc.lt.0.d0) then
	               stop '>> error: root2 not real ***'
	           else if (disc.gt.0.d0) then
	               x1 = ( -b + dsqrt(disc) ) / ( 2.d0*a )
	               x2 = ( -b - dsqrt(disc) ) / ( 2.d0*a )
	           end if
	  end if

      value(1) = x1
	  value(2) = x2

	  end function

!*********************************************************************************************
!      linear interpolation
!*********************************************************************************************
     function interpolation_linear (q,xx) result (value)
     implicit none

	  real*8, intent(in)  :: xx
	  real*8, intent(in)  :: q(2)
      real*8              :: value

      value = (1.d0-xx) * q(1) + xx * q(2)

	 end function
!*********************************************************************************************
!      trilinear interpolation
!*********************************************************************************************
      function interpolation_trilinear (q,xx,yy,zz) result (value)
	  implicit none
	  real*8, intent(in)  :: xx,yy,zz
	  real*8, intent(in)  :: q(2,2,2)
	  real*8              :: value

      value = (1.d0-xx) * (1.d0-yy) * (1.d0-zz) * q(1,1,1) +			   &
                    xx  * (1.d0-yy) * (1.d0-zz) * q(2,1,1) +			   &
                    xx  *       yy  * (1.d0-zz) * q(2,2,1) +    		   &
              (1.d0-xx) *       yy  * (1.d0-zz) * q(1,2,1) +			   &
              (1.d0-xx) * (1.d0-yy) *       zz  * q(1,1,2) +			   &
                    xx  * (1.d0-yy) *       zz  * q(2,1,2) +			   &
                    xx  *       yy  *       zz  * q(2,2,2) +			   &
              (1.d0-xx) *       yy  *       zz  * q(1,2,2)
 
      end function

!********************************************************************************************
!      mean and variance of a vector
!********************************************************************************************
      subroutine mean_var_vect (mean,var,a)
      implicit none
      real*8, intent(in) :: a(:)
	  real*8             :: value,mean,var
	  integer            :: i,n
         
		 n = size(a)
      
	     mean = 0.d0
	     var  = 0.d0
	     do i=1,n
            mean = mean + a(i)
		    var  = var  + a(i)*a(i)
	     end do
         mean = mean / dfloat(n)
         var  = var  / dfloat(n) - mean * mean

	  end subroutine

!********************************************************************************************
!      mean and variance of a vector
!********************************************************************************************
      function norm (vect) result (value)
         implicit none
		 real*8, intent(in) :: vect(:)
		 real*8             :: value
		 integer            :: n
		 n = size(vect)
		 value = dot_product ( vect , vect )
		 value = value ** (1.d0/dble(n)) 
	  end function

!######################################################################!
!                                                                      !
!  SUBROUTINE:  BesselRoot                                             !
!                                                                      !
!  THIS SUBROUTINE COMPUTES THE ROOTS OF THE ZERO-ORDER BESSEL         !
!  FUNCTION OF THE FIRST TYPE THAT ARE NEEDED FOR THE SIMULATION       !
!  OF DIFFUSION IN CYLINDRICAL GEOMETRY.  THE FIRST 20 ROOTS ARE       !
!  TAKEN FROM TABULATED VALUES AND THE REMAINING ROOTS ARE COMPUTED    !
!  FROM AN APPROXIMATION.  (SEE ABROMOWITZ AND STEGUN, 1972 TABLE 9.5  !
!  P. 409 AND FORMULA 9.5.12 P. 371 RESPECTIVELY.]                     !
!                                                                      !              
!######################################################################!
      function BesselRoot(i) result (value)	     
		 implicit none
		 integer, intent(in) :: i
         real*8, parameter   :: pi = 3.14159265359
         real*8              :: root_tab(20),term, var1, var2, var3
		 real*8              :: value

      root_tab(1)  =  2.4048255577D0
      root_tab(2)  =  5.5200781103D0
      root_tab(3)  =  8.6537279129D0
      root_tab(4)  = 11.7915344391D0
      root_tab(5)  = 14.9309177086D0
      root_tab(6)  = 18.0710639679D0
      root_tab(7)  = 21.2116366299D0
      root_tab(8)  = 24.3524715308D0
      root_tab(9)  = 27.4934791320D0
      root_tab(10) = 30.6346064684D0
      root_tab(11) = 33.7758202136D0
      root_tab(12) = 36.9170983537D0
      root_tab(13) = 40.0584257646D0
      root_tab(14) = 43.1997917132D0
      root_tab(15) = 46.3411883717D0
      root_tab(16) = 49.4826098974D0
      root_tab(17) = 52.6240518411D0
      root_tab(18) = 55.7655107550D0
      root_tab(19) = 58.9069839261D0
      root_tab(20) = 62.0484691902D0

      !================================================================!
      ! CALCULATE THE ROOTS                                            !
      !================================================================!
      select case (i<=20 .and. i>0)
                case (.TRUE.);  value = root_tab(i)   ! USE TABULATED VALUES FOR THE FIRST TWENTY ROOTS      
			    case (.FALSE.)    ! APPROXIMATE REMAINING ROOTS                            !
                   var1 = (dble(i)-0.25D0)*pi
                   var2 = 1.D0/(8.D0*var1)
                   var3 = var2*var2
                   term = - (3779.D0/15.D0) - (2.D0*6277237.D0*var3/105.D0)
                   term = (31.D0/3.D0) - (8.D0*var3*term)
                   term = 1.D0 - (4.D0*var3*term)
                   value = var1 + var2*term
			    case default; stop ' *** Coud NOT Calculate Roots of Zero-Order Bessel Function First Kind' 
        end select
      end function BesselRoot


      subroutine Delete_File (fname)
        character(len=*), intent(in) :: fname
		integer                      :: iunit
		call open_fname (fname,iunit)
        close (iunit,status='delete')
	  end subroutine


     Function erf(x) result (value)
     implicit none
	 real*8, intent(in) :: x
	 real*8             :: value,xx,y,t
	   xx = x
        y = x * x
        if (x < 0) xx = -x
        if (xx <= 0.9) then
          value = 1.12838d0 * (xx) * (((((-0.00075757d0 * y + 0.00462963d0) * y - 0.0238095d0) * y + 0.1d0) * y - 0.333333d0) * y + 1.d0)
        elseif (xx <= 3 .and. xx > 0.9) then
          t = 1.d0 / (1.d0 + 0.47047d0 * xx)
          value = 1.d0 - t * ((0.7478556d0 * t - 0.0958798d0) * t + 0.3480242d0) * dexp(-y)
        elseif (xx > 3) then
          value = 1 
        end if
        if (x < 0) value = - value
     end function


! ***************************************************************************      
! ***************************************************************************
!        Get the inverse of a square matrix
! ***************************************************************************

    ! Returns the inverse of a matrix calculated by finding the LU
    ! decomposition.  Depends on LAPACK.
      function inv(A) result(Ainv)
      !use library_lapack
      implicit none
      external DGETRF,DGETRI
      real*8, dimension(:,:), intent(in) :: A
      real*8, dimension(size(A,1),size(A,2)) :: Ainv
    
      real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: n, info
    
    ! External procedures defined in LAPACK
    ! external DGETRF
    ! external DGETRI
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A,1)
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
      call DGETRF(n, n, Ainv, n, ipiv, info)
    
      if (info /= 0) then
        stop 'Matrix is numerically singular!'
      end if
    
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
      call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
      if (info /= 0) then
        stop 'Matrix inversion failed!'
      end if
      end function inv

      
    !subroutine inverse(M,c,n)
    !!============================================================
    !! Inverse matrix
    !! Method: Based on Doolittle LU factorization for Ax=b
    !! Alex G. December 2009
    !!-----------------------------------------------------------
    !! input ...
    !! a(n,n) - array of coefficients for matrix A
    !! n      - dimension
    !! output ...
    !! c(n,n) - inverse matrix of A
    !! comments ...
    !! the original matrix a(n,n) will be destroyed 
    !! during the calculation
    !!===========================================================
    !implicit none 
    !integer, intent(in)      :: n
    !real*8,  intent(in)      :: M(n,n)
    !real*8,  intent(inout)   :: c(n,n)
    !real*8                   :: L(n,n), U(n,n), b(n), d(n), x(n), a(n,n)
    !real*8                   :: coeff
    !integer                  :: i, j, k
    !
    !a = M
    !! step 0: initialization for matrices L and U and b
    !! Fortran 90/95 aloows such operations on matrices
    !L=0.0
    !U=0.0
    !b=0.0
    !
    !! step 1: forward elimination
    !do k=1, n-1
    !   do i=k+1,n
    !      coeff=a(i,k)/a(k,k)
    !      L(i,k) = coeff
    !      do j=k+1,n
    !         a(i,j) = a(i,j)-coeff*a(k,j)
    !      end do
    !   end do
    !end do
    !
    !! Step 2: prepare L and U matrices 
    !! L matrix is a matrix of the elimination coefficient
    !! + the diagonal elements are 1.0
    !do i=1,n
    !  L(i,i) = 1.0
    !end do
    !! U matrix is the upper triangular part of A
    !do j=1,n
    !  do i=1,j
    !    U(i,j) = a(i,j)
    !  end do
    !end do
    !
    !! Step 3: compute columns of the inverse matrix C
    !do k=1,n
    !  b(k)=1.0
    !  d(1) = b(1)
    !! Step 3a: Solve Ld=b using the forward substitution
    !  do i=2,n
    !    d(i)=b(i)
    !    do j=1,i-1
    !      d(i) = d(i) - L(i,j)*d(j)
    !    end do
    !  end do
    !! Step 3b: Solve Ux=d using the back substitution
    !  x(n)=d(n)/U(n,n)
    !  do i = n-1,1,-1
    !    x(i) = d(i)
    !    do j=n,i+1,-1
    !      x(i)=x(i)-U(i,j)*x(j)
    !    end do
    !    x(i) = x(i)/u(i,i)
    !  end do
    !! Step 3c: fill the solutions x(n) into column k of C
    !  do i=1,n
    !    c(i,k) = x(i)
    !  end do
    !  b(k)=0.0
    !end do
    !end subroutine inverse
      
      
   end module

     