

module list_class
   
	implicit none

	private
	public :: list_cl     ! class
	public :: partID_cl
	public ::                                             & ! methods
              initialize_list_                          , &
			  add_particle_to_list_                     , &
			  delete_particle_in_list_                  , &            
              add_move_to_plumeparticle_ 
              
    type partID_cl    
         !attributes of a particle
         integer, pointer :: id  !particle number
		 real*8,  pointer :: tp => null()
		 real*8,  pointer :: xp => null()
		 real*8,  pointer :: yp => null()
		 real*8,  pointer :: zp => null()
		 real*8,  pointer :: mp => null()
		 real*8,  pointer :: rp => null() !x,y,z=position,rp=retardation,mp=mass
		 !neighbor particles      
         type (partID_cl), pointer :: next
         type (partID_cl), pointer :: prev
    end type      
    
	type list_cl
		 integer  :: np      = 0             !total number of particles
		 real*8   :: mtot    = 0             !total mass of particles in list
         type(partID_cl), pointer :: head
         type(partID_cl), pointer :: tail         
    end type list_cl
    
    contains
    
    subroutine initialize_list_ (list) !initialize the empty list
       implicit none
       type(list_cl), pointer :: list
		  allocate(list)
		  list%np      = 0                      
          nullify (list%head,list%tail)
    end subroutine
    
    subroutine add_move_to_plumeparticle_ (part,dxp)
       implicit none
       type(partID_cl),  pointer        :: part 
	   real*8,            intent(in)    :: dxp(3)
	       part % xp = part % xp + dxp(1)
	       part % yp = part % yp + dxp(2)
	       part % zp = part % zp + dxp(3)
    end subroutine    
    
   
    subroutine add_particle_to_list_ (list,id,xp,yp,zp,mp,rp) !add a new particle to the end of the list
       implicit none
       type(list_cl), pointer    :: list
       type(partID_cl), pointer  :: part
	   integer                   :: id
	   integer                   :: err
	   real*8                    :: xp,yp,zp,mp,rp
	   integer                   :: err1,err2,err3,err4,err5,err6                
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
            if(associated(list%head)) then !plume is not empty
                list%tail%next => part
                nullify(part%next)
                part%prev => list%tail 
                list%tail => part
            else                  !plume is empty
                list%head => part
                list%tail => part
                nullify(list%tail%next)
            end if
            list%np = list%np + 1
            list%mtot = list%mtot + mp     
    end subroutine
        
    subroutine delete_particle_in_list_ (part,list)
       implicit none
       type(list_cl),    pointer :: list
       type(partID_cl),  pointer :: part
       type(partID_cl),  pointer :: next            
       real*8                    :: mp
          mp=part%mp
          !save next particle
          next => part%next
          if (.not.associated(part%prev)) then  
                list%head => part%next  !delete the first item in list
          else
                part%prev%next => part%next 
          end if 
          if (.not.associated(part%next)) then
                list%tail => part%prev
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
          list%np = list%np - 1
          list%mtot = list%mtot-mp         
     end subroutine



end module list_class