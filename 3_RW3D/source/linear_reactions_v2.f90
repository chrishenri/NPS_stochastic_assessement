module linear_reactions

implicit none
public
contains

!----------------------------------------------------------------------------------------------
!  CHANGE PARTICLE PROPERTIES FOLLOWING ACTIVATED PACKAGES (SORPTION, MASS TRANSFER, DECAY)
!            calculate transition probability matrix and perform bernouilli trial
!----------------------------------------------------------------------------------------------
    subroutine run_mass_trans_and_decay_network (particle,dt,iseed)
    use global_variables, only: nspecie, nzoneim, nspedecay
    use particle_class
	use heterogeneity_flags
	use loops_particles, only: nmove
	use to_solve
	implicit none

	type(particle_cl),  intent(inout)      :: particle
    integer,            intent(inout)      :: iseed
	real*8,             intent(inout)      :: dt
    logical                                :: heterogeneous
    real*8,  save                          :: dt_previous = -1.d0
    real*8                                 :: dtm
    integer                                :: matdim, i

    !results : eigensyst. and trans. prob matrix
    real*8, allocatable, save              :: P_mmt(:,:), eigenval_mmt(:), eigenvec_mmt(:,:), eigenvecInv_mmt(:,:)
    real*8, allocatable, save              :: P_dec(:,:), eigenval_dec(:), eigenvec_dec(:,:), eigenvecInv_dec(:,:)
    real*8, allocatable, save              :: P_cpl(:,:), eigenval_cpl(:), eigenvec_cpl(:,:), eigenvecInv_cpl(:,:)


!***********************************
    if ( (.not.mass_transACTION) .AND. (.not.decayACTION) .AND. (.not.sorptionTYPE == 'CHTM') ) return

    ! FOR CONTINOUS HISTORY TIME METHOD
    if ( sorptionACTION .AND. sorptionTYPE == 'CHTM') then
        call dtmobile_CHTM_ (dtm,particle,dt,iseed)
        dt = dtm
    end if

    ! FOR MASS TRANSFER AND/OR DECAY NETWORK
    if ( mass_transACTION .OR. decayACTION) then

      ! Check if the system to solve is heterogeneous
        if (mass_transACTION)       heterogeneous = (.not.sorption_homogeneous).OR.(.not.decay_homogeneous).OR.(.not.mass_trans_homogeneous).OR.(.not.poro_homogeneous)
        if (.not. mass_transACTION) heterogeneous = (.not.sorption_homogeneous).OR.(.not.decay_homogeneous)

      ! Save previous particle state
        particle % specie % prenum = particle % specie % num
        particle % zone   % prenum = particle % zone   % num

      ! Check conditions to (re)calculate transition matrix
        if (( .not.heterogeneous .AND. dt_previous < 0 )             .OR.  &  !homo.   conditions & first move   OR
	        ( heterogeneous      .AND. particle%control%switchcell ) .OR.  &  !hetero. conditions & switch cell  OR
	        ( dt /= dt_previous ) ) then                                      !different time step

          ! Solve decay for engaged species
            if (( decayACTION .AND. .not.mass_transACTION ) .OR. ( decayACTION .AND. mass_transACTION .AND. nspedecay/=nspecie )) then
                matdim = nspedecay
                ! Get eigensystem
                if ((.not. heterogeneous .AND. nmove == 1) .OR. (heterogeneous .AND. particle%control%switchcell)) then
                    if ( .not. allocated(eigenvec_dec) )         allocate ( eigenvec_dec(matdim,matdim) )
                    if ( .not. allocated(eigenvecInv_dec) )      allocate ( eigenvecInv_dec(matdim,matdim) )
                    if ( .not. allocated(eigenval_dec) )         allocate ( eigenval_dec(matdim) )
                    if ( .not. allocated(P_dec) )                allocate ( P_dec(matdim,matdim) )
                    call get_eigensystem_ ( particle,matdim,eigenval_dec,eigenvec_dec,eigenvecInv_dec,.TRUE.,.FALSE. )
                end if
                ! Get trans. prob. matrix
                call get_prob_matrix_ ( particle,matdim,eigenval_dec,eigenvec_dec,eigenvecInv_dec,dt,P_dec )
            end if


          ! Solve MMT for all species
            if (( .not.decayACTION .AND. mass_transACTION ) .OR. ( decayACTION .AND. mass_transACTION .AND. nspedecay/=nspecie )) then
                matdim = nspecie * (nzoneim+1)
                ! Get eigensystem
                if ((.not. heterogeneous .AND. nmove == 1) .OR. (heterogeneous .AND. particle%control%switchcell)) then
                    if ( .not. allocated(eigenvec_mmt) )         allocate ( eigenvec_mmt(matdim,matdim) )
                    if ( .not. allocated(eigenvecInv_mmt) )      allocate ( eigenvecInv_mmt(matdim,matdim) )
                    if ( .not. allocated(eigenval_mmt) )         allocate ( eigenval_mmt(matdim) )
                    if ( .not. allocated(P_mmt) )                allocate ( P_mmt(matdim,matdim) )
                    call get_eigensystem_ ( particle,matdim,eigenval_mmt,eigenvec_mmt,eigenvecInv_mmt,.FALSE.,.TRUE. )
                end if
                ! Get trans. prob. matrix
	            call get_prob_matrix_ ( particle,matdim,eigenval_mmt,eigenvec_mmt,eigenvecInv_mmt,dt,P_mmt )
            end if


          ! Solve MMT and decay network in a coupled manner
            if  (( decayACTION .AND. mass_transACTION .AND. nspedecay==nspecie )) then
                matdim = nspecie * (nzoneim+1)
                ! Get eigensystem
                if ((.not. heterogeneous .AND. nmove == 1) .OR. (heterogeneous .AND. particle%control%switchcell)) then
                    if ( .not. allocated(eigenvec_cpl) )         allocate ( eigenvec_cpl(matdim,matdim) )
                    if ( .not. allocated(eigenvecInv_cpl) )      allocate ( eigenvecInv_cpl(matdim,matdim) )
                    if ( .not. allocated(eigenval_cpl) )         allocate ( eigenval_cpl(matdim) )
                    if ( .not. allocated(P_cpl) )                allocate ( P_cpl(matdim,matdim) )
                    call get_eigensystem_ ( particle,matdim,eigenval_cpl,eigenvec_cpl,eigenvecInv_cpl,.TRUE.,.TRUE. )
                end if
                ! Get trans. prob. matrix
	            call get_prob_matrix_ ( particle,matdim,eigenval_cpl,eigenvec_cpl,eigenvecInv_cpl,dt,P_cpl )

            end if

        end if

        if (allocated(P_dec)) call new_particle_state ( particle,nspedecay,iseed,.TRUE.,.FALSE.,P_dec )
        if (allocated(P_mmt)) call new_particle_state ( particle,nspecie*(nzoneim+1),iseed,.FALSE.,.TRUE.,P_mmt )
        if (allocated(P_cpl)) call new_particle_state ( particle,nspecie*(nzoneim+1),iseed,.TRUE.,.TRUE.,P_cpl )

        !........ save previous time step
        dt_previous = dt

    end if
    end subroutine
    
    
!*****************************************************************************************
!  PERFORM BERNOUILLI TRIAL TO SELECT NEW PARTICLE STATE
!*****************************************************************************************
    subroutine new_particle_state ( particle,matdim,iseed,solve_decay,solve_mmt,P )
    use particle_class
    use global_variables,   only: nspecie, nzoneim, nspedecay, idspedecay
    use gslib,              only: rand2
    implicit none

    type(particle_cl),  intent(inout)       :: particle
    integer,            intent(inout)       :: iseed
    real*8,             intent(in)          :: P(matdim,matdim)
    integer,            intent(in)          :: matdim
    logical,            intent(in)          :: solve_decay,solve_mmt
    
    real*8                                  :: rand,sumP,test_index,ir
    integer                                 :: i,indexpartin,indexpartout,nspe,nzim,newspe
    
    if (particle%control%remove) return
    rand = rand2(iseed)

    ! Initial particle state
    if (solve_decay .AND. .not.solve_mmt)       then
        nspe = nspedecay
        nzim = 0
        do i=1,nspe
            if (particle%specie%num==idspedecay(i)) exit
        end do
        indexpartin = i

    elseif (.not.solve_decay .AND. solve_mmt)   then
        nspe = nspecie
        nzim = nzoneim
        indexpartin = (particle%specie%num-1)*(nzim+1) + particle%zone%num+1
        
    elseif (solve_decay .AND. solve_mmt)        then
        nspe = nspecie
        nzim = nzoneim
        do i=1,nspe
            if (particle%specie%num==idspedecay(i)) exit
        end do
        indexpartin = (i-1)*(nzim+1) + particle%zone%num+1
        
    end if


    sumP = 0.d0
    do i=1,matdim
        sumP = sumP + P(i,indexpartin)
    enddo
    if (rand>sumP) then       !particle out of considered chemical system
        particle % specie % num = nspecie+1
    else                      !particle stay into system => continue!

    ! Next particle state
    sumP=0.d0
    do i=1,matdim
        sumP=sumP+P(i,indexpartin)
        if (sumP>rand) exit
    end do
    indexpartout=i

    do i=1,nspe
        ir = real(i)
        test_index = real(indexpartout)/real(nzim+1)
        if ( test_index <= ir) then 
            !particle%specie%num = i
            newspe = i
            exit
        end if
    end do

    ! Final particle state
    if (solve_decay) particle%specie%num = idspedecay(newspe)
    particle%zone%num = indexpartout - (newspe-1)*(nzim+1) - 1
    
    end if      !particle spe > nspe (random number > probabilities sum)

    if (particle%specie%num == nspecie+1) particle%control%remove = .TRUE.

    end subroutine



!*****************************************************************************************
!  TRANSITION PROB. MATRIX
!*****************************************************************************************
    subroutine get_eigensystem_ (particle,matdim,eigenval,eigenvec,eigenvecInv,solve_decay,solve_mmt)
    use particle_class
    use to_solve,        only: decayTYPE
    use gslib,           only: inv
    implicit none

    type(particle_cl),   intent(inout)       :: particle
    real*8,              intent(inout)       :: eigenval(matdim), eigenvec(matdim,matdim), eigenvecInv(matdim,matdim)
    logical,             intent(in)          :: solve_decay, solve_mmt
    integer,             intent(in)          :: matdim
    real*8, dimension (matdim,matdim)        :: MAT, eigenvecL, temp
    real*8, dimension (6*matdim,6*matdim)    :: WORK
    real*8, dimension (matdim)               :: eigenvali
    integer                                  :: INFO

    ! GET THE EIGENSYSTEM OF THE PROBLEM

        !analytical solution for serial network
        if (.not.solve_mmt .AND. solve_decay .AND. decayTYPE == 'SERIAL' ) then
            call eigenserial (particle,eigenval,eigenvec,eigenvecInv)

        !using solver for a general solution
        else
            call defineMAT (particle, matdim, solve_decay, solve_mmt, MAT)
            call DGEEV( 'V', 'V', matdim, MAT, matdim, eigenval, eigenvali, &       ! Get Eigensystem using LAPACK subroutines
                        eigenvecL, matdim, eigenvec, matdim, WORK, matdim*6, INFO )
            eigenvecInv = inv(eigenvec)
            !call inverse(eigenvec,eigenvecInv,matdim)
        end if

        if (decayTYPE == 'SERIAL_MOMENTS') then 
            particle % decay % serial_mom % S       = eigenvec
            particle % decay % serial_mom % Sinv    = eigenvecInv
        end if

    end subroutine



!***************************************************************************************************
!   DEFINE REACTION MATRIX FOR MASS TRANSFER AND/OR GENERIC DECAY NETWORK
!***************************************************************************************************
    subroutine defineMAT (particle, matdim, solve_decay, solve_mmt, MAT)
    use particle_class
    use global_variables, only: nzoneim, nspecie, nspedecay, idspedecay
    use to_solve
    implicit none

    integer,              intent(in)      :: matdim
    type(particle_cl),    intent(inout)   :: particle
    real*8,               intent(inout)   :: MAT(matdim,matdim)
    logical,              intent(in)      :: solve_decay, solve_mmt

    !dummies
    real*8, allocatable                   :: k(:,:), R(:,:), y(:,:), KM(:,:,:), alpha2(:,:), beta2(:,:)
    real*8                                :: SUM
    !logical                               :: analy_sol = .FALSE.
    !logical                               :: reac_immo = .TRUE.
    integer                               :: ispe,jspe,izone,countI,countJ,I,J,m,n,kk, nspe, nzim

    !Definition of Reaction Parameters
    if (solve_decay .AND. .not.solve_mmt)       then
        allocate (k(nspedecay,1),R(nspedecay,1),y(nspedecay,nspedecay))
        allocate (alpha2(nspedecay,1), beta2(nspedecay,1))
        allocate (KM(nspedecay,nspedecay,1))
        nspe    = nspedecay !only species in decay network
        nzim    = 0         !no immobile zones
        k       = 0.d0
        y       = 0.d0
        R       = 1.d0
        alpha2  = 0.d0
        beta2   = 1.d0
        do ispe=1,nspedecay
            if (particle%zone%num == 0) then ! in mobile domain
                k(ispe,1)  = particle % decay % k(ispe)
                if (sorptionACTION .and. sorptionTYPE == 'LINEAR') R(ispe,1) = particle % sorption % linear_sorp % R(idspedecay(ispe))
            else if (mass_transACTION .AND. particle%zone%num > 0) then
                select case (mass_transTYPE)
                case ('MULTIRATE','COMPOSITE_MEDIA') 
                    k(jspe,1)  = particle % decay % kim(ispe,particle%zone%num)
                    if (sorptionACTION .and. sorptionTYPE == 'LINEAR') R(ispe,1) = particle % sorption % linear_sorp % Rim(idspedecay(ispe),particle%zone%num)
                case ('SPHERICAL_DIFFUSION','LAYERED_DIFFUSION','CYLINDRICAL_DIFFUSION','POWER_LAW','LOGNORMAL_LAW')
                    k(jspe,1)  = particle % decay % kim(ispe,1)
                    if (sorptionACTION .and. sorptionTYPE == 'LINEAR') R(ispe,1) = particle % sorption % linear_sorp % Rim(idspedecay(ispe),1)
                end select
            end if

            do jspe=1,nspedecay
                y(ispe,ispe) = -1.d0
                y(ispe,jspe)    = particle % decay % y(ispe,jspe)
                KM(ispe,jspe,1) = y(ispe,jspe) * k(jspe,1) / R(jspe,1);
            end do
        end do


    else if (solve_mmt .AND. .not.solve_decay)  then
        allocate (KM(nspecie,nspecie,nzoneim+1))
        allocate (alpha2(nspecie,nzoneim+1), beta2(nspecie,nzoneim+1))
        nspe    = nspecie !all species
        nzim    = nzoneim !all zones
        KM      = 0.d0
        alpha2  = 0.d0
        beta2   = 1.d0
        do ispe=1,nspecie; do izone=2,nzoneim+1
            alpha2(ispe,izone) = particle % zone % alpha(ispe,izone-1)
            beta2(ispe,izone)  = particle % zone % beta (ispe,izone-1)
        end do; end do


    else if (solve_mmt .AND. solve_decay)       then
        allocate (k(nspedecay,nzoneim+1),R(nspedecay,nzoneim+1),y(nspedecay,nspedecay))
        allocate (KM(nspedecay,nspedecay,nzoneim+1))
        allocate (alpha2(nspedecay,nzoneim+1), beta2(nspedecay,nzoneim+1))
        nspe    = nspedecay !all species
        nzim    = nzoneim   !all zones
        k       = 0.d0
        y       = 0.d0
        R       = 1.d0
        alpha2  = 0.d0
        beta2   = 1.d0
        do ispe=1,nspedecay
        do izone=1,nzoneim+1
            if (izone==1) then
                k(ispe,izone)  = particle % decay % k(ispe)
                if (sorptionACTION .and. sorptionTYPE == 'LINEAR') R(ispe,1) = particle % sorption % linear_sorp % R(idspedecay(ispe))
            else
                select case (mass_transTYPE)
                case ('MULTIRATE')
                    k(ispe,izone) = particle % decay % kim(ispe,izone-1)
                    !if (k(jspe,izone) > 0.d0) reac_immo = .TRUE. !check if reaction in immobile zones
                    if (sorptionACTION .and. sorptionTYPE == 'LINEAR')     R(ispe,izone) = particle % sorption % linear_sorp % Rim(idspedecay(ispe),izone-1)
                case ('SPHERICAL_DIFFUSION', 'LAYERED_DIFFUSION', 'CYLINDRICAL_DIFFUSION', 'POWER_LAW', 'LOGNORMAL_LAW')
                    k(ispe,izone) = particle % decay % kim(ispe,1)
                    !if (k(jspe,izone) > 0.d0) reac_immo = .TRUE. !check if reaction in immobile zones
                    if (sorptionACTION .and. sorptionTYPE == 'LINEAR')     R(ispe,izone) = particle % sorption % linear_sorp % Rim(idspedecay(ispe),1)
                end select
            end if
            y = 0.d0
            do jspe=1,nspedecay
                if (jspe==ispe) y(ispe,jspe) = -1.d0
                if (decayTYPE=='SERIAL') then
                    if (jspe==ispe-1) y(ispe,jspe) = particle % decay % y(ispe,jspe)
                else
                    y(ispe,jspe)    = particle % decay % y(ispe,jspe)
                end if
                KM(ispe,jspe,izone) = y(ispe,jspe) * k(jspe,1) / R(jspe,1);
            end do

            if (izone>1) then
                alpha2(ispe,izone) = particle % zone % alpha(ispe,izone-1)
                beta2(ispe,izone)  = particle % zone % beta (ispe,izone-1)
            end if
        end do; end do

    end if

    !if (solve_decay .AND. solve_mmt .AND. nspedecay == 2 .AND. nzoneim == 1 .AND. (.not.reac_immo)) analy_sol = .TRUE.

    ! General Solution: "Reaction" Matrice A^(-1)*B
    !if (.NOT.analy_sol) then
    MAT     = 0.d0
    I       = 1
    J       = 1
    countI  = 0
    countJ  = 0
    do m=1,(nspe*(nzim+1))
        countI=countI+1
        do n=1,(nspe*(nzim+1))
            countJ=countJ+1

            if (I==J) then                  !diagonal-blocks

                if (countI==countJ) then    !* diagonal of the diagonal-blocks
                    if (countI==1)  then    !** first term
                        SUM = 0;
                        do kk=2,nzim+1
                            SUM = SUM - alpha2(I,kk)*beta2(I,kk)
                        end do
                        MAT(m,n) = KM(I,I,1)+SUM;
                    else                    !** other terms
                        MAT(m,n) = KM(I,I,countJ)-alpha2(I,countJ)
                    end if
                elseif (countI==1) then
                    MAT(m,n) = alpha2(I,countJ)!-KM(I,I,countJ)
                elseif (countJ==1) then
                    MAT(m,n) = alpha2(I,countI)*beta2(I,countI)
                end if

                if (countJ==nzim+1) then !* reinitial the counter and change block in the J direction
                    countJ = 0
                    if (J+1 <= nspe) then
                        J = J+1
                    else
                        J = 1
                    end if
                end if

            else                            !for the other blocks

                if (countI==countJ)   MAT(m,n) = KM(I,J,countJ)

                if (countJ==nzim+1) then ! reinitial the counter and change block in the J direction
                    countJ = 0
                    if (J+1 <= nspe) then
                        J = J+1
                    else
                        J = 1
                    end if
                end if

            end if

        end do

        if (countI==nzim+1) then ! reinitial the counter and change block in the I direction
            countI = 0
            I      = I+1
        end if

    end do
    !end if

    end subroutine


!*****************************************************************************************
!  TRANSITION PROB. MATRIX
!*****************************************************************************************
    subroutine get_prob_matrix_ (particle,matdim,eigenval,eigenvec,eigenvecInv,dt,P)
    use particle_class
    use to_solve, only:decayTYPE
    implicit none

    type(particle_cl),   intent(inout)              :: particle
    real*8,  intent(inout)                          :: P(matdim,matdim)
    real*8,  intent(in)                             :: eigenval(matdim), eigenvec(matdim,matdim), eigenvecInv(matdim,matdim)
    real*8,  intent(in)                             :: dt
    integer, intent(in)                             :: matdim
    real*8, dimension (matdim,matdim)               :: Diag, temp
    integer                                         :: i
    
    ! GET THE TRANSITION PROB. MATRIX
    Diag = 0.d0
    do i=1,matdim
       Diag(i,i) = exp(eigenval(i)*dt)
    end do
    
    !P = matmul(eigenvec,matmul(Diag,eigenvecInv))
    call DGEMM('N','N',matdim,matdim,matdim,1.d0,Diag,matdim,eigenvecInv,matdim,0.d0,temp,matdim)   !matrices multiplication
    call DGEMM('N','N',matdim,matdim,matdim,1.d0,eigenvec,matdim,temp,matdim,0.d0,P,matdim)         !matrices multiplication
    
    if (decayTYPE == 'SERIAL_MOMENTS') particle % decay % serial_mom % P = P

    end subroutine


!***************************************************************************************************
!   GET EIGENSYSTEM ANALYTICALLY FOR SERIAL REACTION NETWORK
!***************************************************************************************************
    subroutine eigenserial(particle,D,S,Sinv)
	    use particle_class
	    use global_variables, only: nspedecay, idspedecay
        use to_solve,         only: sorptionACTION, sorptionTYPE

        type(particle_cl),    intent(inout)   :: particle
        real*8,               intent(inout)   :: S(nspedecay,nspedecay), Sinv(nspedecay,nspedecay), D(nspedecay)
        integer                               :: i,j,m
        real*8                                :: k(nspedecay), y(nspedecay), R(nspedecay) 
        real*8                                :: prodS,prodSinv

        ! reaction coefficients definition
        k = 0.d0
        y = 0.d0
        R = 1.d0
        do i=1, nspedecay
	        k(i)   = particle % decay % k(i)
	        if (i>1) y(i)   = particle % decay % y(i,i-1)
	        if (sorptionACTION .and. sorptionTYPE == 'LINEAR') R(i) = particle % sorption % linear_sorp % R(idspedecay(i))
	    enddo

        S    = 0.d0
        Sinv = 0.d0
        do i=1,nspedecay; do j=1,nspedecay
            if (i>=j) then
                ! get eigenvectors S
                prodS=1.0
                do m=j,i-1
                    prodS = prodS* ((k(m)*y(m+1))/(R(j)*k(m+1)-R(m+1)*k(j)))  ;enddo
                S(i,j)   = (R(j)**(i-j-1))*R(i)*prodS
                ! get Sinv
                prodSinv=1.0
                do m=j,i-1
                    prodSinv = prodSinv* (-(k(m)*y(m+1))/(R(m)*k(i)-R(i)*k(m)))  ;enddo
                Sinv(i,j)   = (R(i)**(i-j))*prodSinv
            end if
        end do; end do

        !get the diagonal matrix (eigenvalues of the reaction matrix)
        D = 0.d0
        do i=1,nspedecay
           D(i) = -k(i)/R(i)
        end do

    end subroutine


!*****************************************************************************************
!  CONTINUOUS HISTORY METHOD ONE SITE
!*****************************************************************************************
    subroutine dtmobile_CHTM_ (dtm,particle,dt,iseed) 
        use gslib, only: rand2
        use particle_class
        implicit none

        type(particle_cl), intent(inout)    :: particle
        integer,           intent(inout)    :: iseed
        real*8,            intent(in)       :: dt
        real*8,            intent(out)      :: dtm
        real*8                              :: Pm,Ps,ts,tm,kf,kb,dts,poro,bd

        dtm = dt

        select case (particle%zone%num == 0)

        case (.TRUE.)  ! starting from mobile state        

           do 
              Pm = rand2(iseed)
			  poro = particle % prop % poro
			  kf   = particle % sorption % CHTM % kf
			  bd   = particle % sorption % CHTM % bd
	          tm = tm -log(Pm)/kf/bd*poro
	          if ( tm+ts>dt ) then
			      dts = ts 
 			      dtm = dt-ts
                  particle % zone % num = 0
 			      return
 	          else
			      Ps = rand2(iseed)
				  kb = particle % sorption % CHTM % kb
 			      ts = ts -log(Ps)/kb  
 			      if ((tm+ts)>dt) then
 				      dtm = tm	
 				      dts = dt - tm
				      particle % zone % num = 1
 				      return
 			      end if
	          end if
	       end do

        case (.FALSE.)   !starting from sorbed state

           do ! this loop never ends
			   Ps = rand2(iseed)
			   kb = particle % sorption % CHTM % kb
	           ts = ts -log(Ps)/kb
	           if ( tm+ts>dt ) then
			      dtm = tm 
 			      dts = dt-tm
				  particle%zone%num = 1
 			      return
 	           else
                  Pm   = rand2(iseed)
			      poro = particle % prop % poro 
			      kf   = particle % sorption % CHTM % kf
			      bd   = particle % sorption % CHTM % bd
 			      tm = tm -log(Pm)/kf/poro*bd  
 			      if ((tm+ts)>dt) then
 				       dts = ts	
 				       dtm = dt - ts
				       particle % zone % num = 0
 				       return
 			      end if
	           end if
             end do
        end select

    end subroutine

!----------------------------------------------------------------------------------------------
 end module