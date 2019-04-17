 module mf2k_class

	implicit none
	

	private
	public :: mf2k_cl
	public ::                        & ! methods
             allocate_mf2k_         ,&
             check_budget_file_     

    type mf2k_cl            
		 character(len=50), pointer  :: file => null()
		 integer, pointer            :: nx => null()
		 integer, pointer            :: ny => null()
		 integer, pointer            :: nz => null()         
	     integer, pointer            :: preCBC  => null() 
	     integer, pointer            :: unitCBC => null()
	     character(len=12), pointer  :: formCBC => null()
	end type

    
    contains
    
    subroutine allocate_mf2k_ (mf2k)
 	  implicit none
      type(mf2k_cl),    intent(inout) :: mf2k
		 allocate (mf2k%nx)
		 allocate (mf2k%ny)
		 allocate (mf2k%nz)          
	     allocate (mf2k%preCBC) 
	     allocate (mf2k%unitCBC)
	     allocate (mf2k%formCBC)
	     allocate (mf2k%file)
		 mf2k%nx = 0
		 mf2k%ny = 0
		 mf2k%nz = 0          
	     mf2k%preCBC = 0 
	     mf2k%unitCBC = 1
	     mf2k%formCBC ='UNFORMATTED'
	     mf2k%file = ' '
    end subroutine   
!*******************************************************
!-----Check for valid budget file  
!*******************************************************  
    subroutine check_budget_file_ (fname,mf2k)
	  use gslib,only: open_fname
	  implicit none
      character(len=*), intent(in)    :: fname
      type(mf2k_cl),    intent(inout) :: mf2k
      integer                         :: iprec,iunit,ncol,nrow,nlay,iufdbg

      !call open_fname (fname,iunit,'unknown','sequential','binary')
      !call open_fname (fdbg,iufdbg,'unknown','append','formatted')
      call open_fname (fname,iunit,'old','stream','unformatted') !good: from Modpath
      !call open_fname (fname,iunit,'old','sequential','unformatted')
      
      call budget_precision (iprec,iunit,ncol,nrow,nlay)
      
      if(iprec.lt.1) then
         close (iunit)
         open(unit=iunit,file=fname,status='old',form='unformatted',access='sequential',err=10)        
         call budget_precision(iprec,iunit,ncol,nrow,nlay)
      end if
      if(iprec.lt.1) then
        write(*,*) 'stopping because budget file is invalid'
        stop
      elseif(iprec.eq.1) then
        write(*,*) ' single precision budget file'
      else if(iprec.eq.2) then
        write(*,*) ' double precision budget file'
      end if
      write(*,*)
      write(*,14) nlay,nrow,ncol
      write(*,14) nlay,nrow,ncol
      write(*,*)
14    format(1x,i10,' layers',i10,' rows',i10,' columns')
     
      mf2k%preCBC = iprec
      mf2k%nx = ncol
      mf2k%ny = nrow
      mf2k%nz = nlay
      mf2k%unitCBC = iunit
      mf2k%file    = fname
      
      return
      
 10   stop 'Could not open budget file from modflow'
     
    end subroutine

!     ******************************************************************
!     Determine single or double precision file type for a MODFLOW
!     budget file:  0=unrecognized, 1=single, 2=double.
!     ******************************************************************
    subroutine budget_precision (iprec,iu,ncol,nrow,nlay)
      use gslib, only: upper_case
      implicit none
      integer, intent(in)    :: iu
      integer, intent(out)   :: iprec,ncol,nrow,nlay
      real*8                 :: deltd,pertimd,totimd,vald
      real                   :: delt,pertim,totim,val
      character(len=16)      :: text1,text2
      integer                :: kstp,kper,icode,nodes
      integer                :: nlst,n,icell,nc,nr,nl
      real, pointer          :: buff(:,:,:)
      real*8, pointer        :: buffd(:,:,:)
!
!  default is unrecognized file
     
      iprec=0
!
!  single check
!
      read(iu,err=100,end=100) kstp,kper,text1,ncol,nrow,nlay
      text1 = upper_case (text1)
      icode=0
      if(nlay.lt.0) then
        nlay=-nlay
        read(iu,err=50,end=50) icode,delt,pertim,totim
      end if
14    format(1x,i10,' layers',i10,' rows',i10,' columns')
      if(ncol.lt.1 .or. nrow.lt.1 .or. nlay.lt.1) go to 100
      if(ncol.gt.100000000 .or.nrow.gt.100000000 .or. nlay.gt.100000000) go to 100
      if(ncol*nrow.gt.100000000 .or. ncol*nlay.gt.100000000 .or. nrow*nlay.gt.100000000) go to 100
      allocate (buff(ncol,nrow,nlay))
      allocate (buffd(ncol,nrow,nlay))
      nodes=ncol*nrow*nlay
!
!  read data depending on icode.  icode 0,1, or 2 are the only allowed
!  values because the first budget terms must be from the internal
!  flow package (bcf,lpf, or huf).
!
      if(icode.eq.0 .or. icode.eq.1) then
         read(iu,err=50,end=50) buff
      else if(icode.eq.2) then
         read(iu,err=50,end=50) nlst
         if(nlst.lt.0) go to 50
         if(nlst.gt.0) then
            do n=1,nlst
               read(iu,end=50,err=50) icell,val
               if(icell.le.0 .or. icell.gt.nodes) go to 50
            end do
         end if
      else
         go to 100
      end if
!
!  read 2nd header and check for valid type.
!
      read(iu,err=50,end=50) kstp,kper,text2
      text2 = upper_case (text2)
      if(text1.eq.'         STORAGE' .and. text2.eq.'   CONSTANT HEAD') then
           iprec=1
           go to 100
      else if(text1.eq.'   CONSTANT HEAD' .and. text2.eq.'FLOW RIGHT FACE ') then
           iprec=1
           go to 100
      end if
!
!  double check
!
50    rewind(iu)
      read(iu,err=100,end=100) kstp,kper,text1,nc,nr,nl
      icode=0
      if(nl.lt.0) then
        nl=-nl
        read(iu,err=100,end=100) icode,deltd,pertimd,totimd
      end if
!
!  read data depending on icode.  icode 0,1, or 2 are the only allowed
!  values because the first budget terms must be from the internal
!  flow package (bcf,lpf, or huf).
!
      if(icode.eq.0 .or. icode.eq.1) then
         read(iu,err=100,end=100) buffd
      else if(icode.eq.2) then
         read(iu,err=100,end=100) nlst
         if(nlst.lt.0) go to 100
         if(nlst.gt.0) then
            do n=1,nlst
               read(iu,end=100,err=100) icell,vald
               if(icell.le.0 .or. icell.gt.nodes) go to 100
            end do
         end if
      else
         go to 100
      end if
!
!  read 2nd header and check for valid type.
!
      read(iu,err=100,end=100) kstp,kper,text2
      if(text1.eq.'         storage' .and.  text2.eq.'   constant head') then
           iprec=2
      else if(text1.eq.'   constant head' .and. text2.eq.'flow right face ') then
           iprec=2
      end if
!
100   rewind(iu)
      if(associated(buff))  deallocate (buff)
      if(associated(buffd)) deallocate (buffd)
      return
      end subroutine


    
    end module mf2k_class