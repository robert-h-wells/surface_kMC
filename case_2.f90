!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
! Simulates Langmuir Hinshelwood Mechanism : A(g) + B(g) -> C(g)
! Individual reactions :  A(g) + * -> A*  ,  B(g) + *-> B*
!                         A* + B* -> C*  ,  C* -> C(g) + *    
!
!     Also contains species D that just diffuses and blocks sites.  
!
! Code is updated for nearest neighbor interactions. 
!
! Contains the following
!   init_case          : sets up reaction, diffusion, and interaction energies
!   reaction_rates     : determines if nearby molecules allow reaction to occur
!   reaction_results   : gives result for specific reaction
!   ..._sohr           : same as above except it checks for atom following
!   generate_neighbors : finds neighbors for given atom
!   ..._next_neighbors : not updated (6/11/18), next neighbors finder
!   equations          : dif eq's used by LSODE
!   lateral            : calculates energy due to lateral interactions
!   rate_const         : updates rate constants due to lateral interactions
!   init_rateconst     : calculates initial rate constants
!   adsorption         : determines where adsorption takes place
!
! Author : Rob Wells              Date : 6/11/18 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!======================================================================================================!
MODULE case_2           
!======================================================================================================!
  implicit none
  save

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: nrate = 3                 ! possible rate constants of each atom
  integer, parameter :: atomType = 4              ! # of different types of atoms
  integer, parameter :: nmol = atomType + 1       ! Extra one for free site coverage, used for LSODE
  integer, parameter :: sym = 4                   ! symmetry number of surface, # of neighbors
  integer, parameter :: gas_phase = 100

  real(dp), parameter :: kBev = 8.617332478D-5
  real(dp), parameter :: kB = 1.3806488D-23
  real(dp), parameter :: pi = 3.14159265359d0

  integer, allocatable :: initial(:),init_num(:)
  integer :: natom,bound(4),nsurface,spline_points,atom_follow,atom_finish
  integer :: natom_save,nsurface_save
  real(dp), allocatable :: rate_val(:),saves(:,:,:,:)
  real(dp) :: Ea(1),Edes(atomType),Edif(atomType),ad_rate(2),des_rate(atomType)
  real(dp) :: diffusion_val(atomType),p_atm(2),area,timer(6)
  real(dp) :: lateral_rate(sym,atomType),lateral_dif_rate(sym,atomType),lateral_energy(atomType)

  type molecule
    integer :: id,pos(2)            ! pos(1) is x-axis, pos(2) is y-axis
    integer :: neighbor(4)          ! 1 right, 2 left, 3 up, 4 down 
    integer :: next_neighbor(8)
    real(dp), allocatable :: rates(:,:)
  end type molecule 

  type(molecule), allocatable :: atom(:),atom_save(:)
  
  type site
    integer :: atom_num
    integer :: atom_id
  end type site

  type(site), allocatable :: surface(:,:),surface_save(:,:)
  !$omp threadprivate(atom,surface,nsurface,atom_finish,atom_follow)

  contains
!======================================================================================================!
  SUBROUTINE init_case()  
!======================================================================================================!
    implicit none
 
    Ea = [0.20d0]                      ! Grabow, Mavrikakis : Mechanims of Methanol Synthesis
    Edes = [1.00d0,1.00d0,0.08d0,10.00d0]     ! Grabow, Mavrikakis : Mechanims of Methanol Synthesis
    Edif = [0.60d0,0.60d0,1.00d0,0.650d0]    ! Gomer_90 

    lateral_energy = [-0.10d0,-0.10d0,0.0d0,0.0d0]

  END SUBROUTINE init_case
!======================================================================================================!
  SUBROUTINE reaction_rates(n)  
!======================================================================================================!
    implicit none 
    integer, intent(in) :: n
    integer :: i,j

    do i = 1,atomType
      if (atom(n)%id == i) atom(n)%rates(3,1) = des_rate(i)
    end do 

    do j = 1,sym     ! Check if neighbor is suitable for reaction to occur 

      if (atom(n)%neighbor(j) > 0) then

        if (atom(n)%id == 1) then   
          if (atom(atom(n)%neighbor(j))%id == 2) then 
            do i = 1,sym   ! take into account neighbor's lateral interaction effects on rxn rate
              if (atom(atom(n)%neighbor(j))%rates(2,i) /= rate_val(1) .AND. atom(atom(n)%neighbor(j))%rates(2,i) /= 0.0d0) then 
                atom(n)%rates(2,j) = atom(atom(n)%neighbor(j))%rates(2,i)
              else
                atom(n)%rates(2,j) = rate_val(1)
              end if 
            end do 
          end if 
        
        else if (atom(n)%id == 2) then
          if (atom(atom(n)%neighbor(j))%id == 1) then 
            do i = 1,sym   ! take into account neighbor's lateral interaction effects on rxn rate
              if (atom(atom(n)%neighbor(j))%rates(2,i) /= rate_val(1) .AND. atom(atom(n)%neighbor(j))%rates(2,i) /= 0.0d0) then 
                atom(n)%rates(2,j) = atom(atom(n)%neighbor(j))%rates(2,i)          
              else
                atom(n)%rates(2,j) = rate_val(1)
              end if 
            end do 
          end if 
        end if
      end if
    end do

  END SUBROUTINE reaction_rates
!======================================================================================================!
  SUBROUTINE reaction_results(n,total_rate,move)
!======================================================================================================!
    implicit none
    integer, intent(inout)  :: n,move
    real(dp), intent(inout) :: total_rate
    integer                 :: i,j,x,y,choice

    choice = 0

    do i = 1,sym

      if (move == 4+i) then     ! Reaction A*+B*-> C*

        atom(n)%id = 3  ; surface(atom(n)%pos(1),atom(n)%pos(2))%atom_id = 3  ! react atom
        choice = atom(n)%neighbor(i)   !  neighbor it reacts with goes away
        total_rate = total_rate - sum(atom(choice)%rates)
 
        x = atom(choice)%pos(1) ; y = atom(choice)%pos(2) 
        surface(x,y)%atom_id = 0 ; surface(x,y)%atom_num = 0
        atom(choice)%id = 0 ; atom(choice)%pos = [gas_phase,gas_phase]

      else if (move == 8+i) then    ! Unimolecular desorption
        x = atom(n)%pos(1) ; y = atom(n)%pos(2) ; choice = n
        surface(x,y)%atom_id = 0 ; surface(x,y)%atom_num = 0
        atom(n)%id = 0 ; atom(n)%pos = [gas_phase,gas_phase] 
        total_rate = total_rate - sum(atom(choice)%rates)

        go to 323

      end if  
    end do 

    ! search for neighbors next to atom that is gone after event
    323 do j = 1,sym 
      if (atom(choice)%neighbor(j) /= 0) then
        total_rate = total_rate - sum(atom(atom(choice)%neighbor(j))%rates)
        call generate_neighbors(atom(choice)%neighbor(j),0)       
        call reaction_rates(atom(choice)%neighbor(j))
        if (all(lateral_energy /= 0.0d0)) call lateral(atom(choice)%neighbor(j))
        total_rate = total_rate + sum(atom(atom(choice)%neighbor(j))%rates)
      end if  
    end do 

    ! shift atom numbers down so that # of atoms on surface is properly counted
    do j = choice+1,nsurface
      atom(j-1) = atom(j)  
      x = atom(j-1)%pos(1) ; y = atom(j-1)%pos(2)
      if (x /= gas_phase) then 
        surface(x,y)%atom_id = atom(j-1)%id ; surface(x,y)%atom_num = j-1
      end if 
    end do 

    ! shift atoms neighbors down 
    do j = 1,nsurface
      do i = 1,sym 
        if (atom(j)%neighbor(i) >= choice) then
          atom(j)%neighbor(i) = atom(j)%neighbor(i) - 1
        end if 
      end do 
    end do 

    ! account for atom being gone
    atom(nsurface)%id = 0 ; atom(nsurface)%pos = [gas_phase,gas_phase]
    nsurface = nsurface - 1 ; 
    if (n > choice) n = n - 1

  END SUBROUTINE reaction_results
!======================================================================================================!
  SUBROUTINE reaction_results_sohr(n,total_rate,move,sums_check)
!======================================================================================================!
    implicit none
    integer, intent(inout)  :: n,move
    integer, intent(out)    :: sums_check
    real(dp), intent(inout) :: total_rate
    integer                 :: i,j,x,y,choice,choice_x,choice_y,a,b,c,d,sums

    choice = 0 ; sums_check = -1

    do i = 1,sym

      if (move == 4+i) then    ! Reaction A*+B*-> C*

        atom(n)%id = 3  ; surface(atom(n)%pos(1),atom(n)%pos(2))%atom_id = 3  ! react atom
        choice = atom(n)%neighbor(i)   !  neighbor it reacts with goes away
        total_rate = total_rate - sum(atom(choice)%rates)

        ! Determine if atom following is over and find which species it reacted with
        if (n == atom_follow) then 
          atom_finish = 1
          choice_x = atom(choice)%pos(1) ; choice_y = atom(choice)%pos(2)
        else if (choice == atom_follow) then 
          atom_finish = 1
          choice_x = atom(n)%pos(1) ; choice_y = atom(n)%pos(2)
        end if

        if (atom_finish == 1) then 
          a = choice_x+1 ; b = choice_x-1 ; c = choice_y+1 ; d = choice_y-1
          sums = 0 ; sums_check = 0

          ! Incorporate boundary conditions for nearest site
          if (a > bound(2)) a = bound(1)
          if (b < bound(1)) b = bound(2)
          if (c > bound(4)) c = bound(3)
          if (d < bound(3)) d = bound(4)

          if (surface(a,choice_y)%atom_id == 2) sums = sums + 1
          if (surface(b,choice_y)%atom_id == 2) sums = sums + 1
          if (surface(choice_x,c)%atom_id == 2) sums = sums + 1
          if (surface(choice_x,d)%atom_id == 2) sums = sums + 1

          sums_check = sums ; print*,'Sums',sums,sums_check

        end if  


        x = atom(choice)%pos(1) ; y = atom(choice)%pos(2) 
        surface(x,y)%atom_id = 0 ; surface(x,y)%atom_num = 0
        atom(choice)%id = 0 ; atom(choice)%pos = [gas_phase,gas_phase]

      else if (move == 8+i) then    ! Unimolecular desorption
        x = atom(n)%pos(1) ; y = atom(n)%pos(2) ; choice = n
        surface(x,y)%atom_id = 0 ; surface(x,y)%atom_num = 0
        atom(n)%id = 0 ; atom(n)%pos = [gas_phase,gas_phase] 
        total_rate = total_rate - sum(atom(choice)%rates)

        go to 324

      end if  
    end do 

    ! search for neighbors next to atom that is gone after event
    324 do j = 1,sym 
      if (atom(choice)%neighbor(j) /= 0) then
        total_rate = total_rate - sum(atom(atom(choice)%neighbor(j))%rates)
        call generate_neighbors(atom(choice)%neighbor(j),0)       
        call reaction_rates(atom(choice)%neighbor(j))
        if (all(lateral_energy /= 0.0d0)) call lateral(atom(choice)%neighbor(j))
        total_rate = total_rate + sum(atom(atom(choice)%neighbor(j))%rates)
      end if  
    end do 

    ! shift atom numbers down so that # of atoms on surface is properly counted
    if (atom_follow >= choice+1) then
      atom_follow = atom_follow - 1 
    end if
    do j = choice+1,nsurface
      atom(j-1) = atom(j)  
      x = atom(j-1)%pos(1) ; y = atom(j-1)%pos(2)
      if (x /= gas_phase) then 
        surface(x,y)%atom_id = atom(j-1)%id ; surface(x,y)%atom_num = j-1
      end if 
    end do 

    ! shift atoms neighbors down 
    do j = 1,nsurface
      do i = 1,sym 
        if (atom(j)%neighbor(i) >= choice) then
          atom(j)%neighbor(i) = atom(j)%neighbor(i) - 1
        end if 
      end do 
    end do 

    ! account for atom being gone
    atom(nsurface)%id = 0 ; atom(nsurface)%pos = [gas_phase,gas_phase]
    nsurface = nsurface - 1 ; 
    if (n > choice) n = n - 1

  END SUBROUTINE reaction_results_sohr
!======================================================================================================!
  SUBROUTINE generate_neighbors(n,val)  

  ! (1,1) move x+1 : (1,2) move x-1 : (1,3) move y+1 : (1,4) move y-1
!======================================================================================================! 
    implicit none
    integer, intent(in) :: n,val
    integer :: i,j,k,a,b,c,d,x,y

    ! initialize rates
    atom(n)%rates(:,:) = 0.0d0 ; atom(n)%neighbor = 0 
    atom(n)%rates(1,:) = diffusion_val(atom(n)%id)

    ! Set variables for looking at nearest sites 
    x = atom(n)%pos(1) ; y = atom(n)%pos(2)
    a = x+1  ;  b = x-1
    c = y+1  ;  d = y-1
  
    ! Incorporate boundary conditions for nearest site
    if (a > bound(2)) a = bound(1)
    if (b < bound(1)) b = bound(2)
    if (c > bound(4)) c = bound(3)
    if (d < bound(3)) d = bound(4)

    ! Check for neighbors  
    if ( surface(a,y)%atom_id > 0 ) then  
      atom(n)%rates(1,1) = 0.0d0 ; atom(n)%neighbor(1) = surface(a,y)%atom_num
    end if 
    if ( surface(b,y)%atom_id > 0 ) then 
      atom(n)%rates(1,2) = 0.0d0 ; atom(n)%neighbor(2) = surface(b,y)%atom_num
    end if 
    if ( surface(x,c)%atom_id > 0 ) then 
      atom(n)%rates(1,3) = 0.0d0 ; atom(n)%neighbor(3) = surface(x,c)%atom_num
    end if 
    if ( surface(x,d)%atom_id > 0 ) then 
      atom(n)%rates(1,4) = 0.0d0 ; atom(n)%neighbor(4) = surface(x,d)%atom_num
    end if     

  END SUBROUTINE
!======================================================================================================!
  SUBROUTINE lateral(n)                                ! modify Hamiltonian due to neighbors !
!======================================================================================================!
    implicit none 
    integer, intent(in) :: n
    integer :: i,j,k,e_add(atomType),types 

    e_add = 0 ; types = atom(n)%id

    do j = 1,sym  ! loop over neighbors to find lateral interactions
      if (atom(n)%neighbor(j) /= 0) then
        do i = 1,atomType 
          if (types == i .and. atom(atom(n)%neighbor(j))%id == i) e_add(i) = e_add(i) + 1
        end do 
      end if 
    end do

    do k = 1,atomType
      if (e_add(k) /= 0) then 
        do j = 1,sym
          if (atom(n)%rates(1,j) /= 0.0d0) then  ! diffusion
            do i = 1,4 
              if (e_add(k) == i) atom(n)%rates(1,j) = lateral_dif_rate(i,k)
            end do 
          end if 
          if (atom(n)%rates(2,j) /= 0.0d0) then  ! reaction
            do i = 1,4 
              if (e_add(k) == i) atom(n)%rates(2,j) = lateral_rate(i,k)
            end do 
          end if 
        end do 
      end if 
    end do 

  END SUBROUTINE lateral
!======================================================================================================!
  SUBROUTINE rate_const(n,prefactor,e_add)
!======================================================================================================!
    implicit none 
    integer, intent(in) :: n  ! which atom
    integer :: i,j
    real(dp), intent(in) :: prefactor,e_add(:)
    real(dp) :: T,lnrate_val(nrate,sym),e_val(nrate)

    T = 500.0d0

    ! Increase activation energies due to lateral interactions 
    ! Assuming for now that lateral interactions affect reactants and transition species the same,
    ! so the reaction rate constants will not be affected. Just diffusion
    
    do i = 1,nrate-1
      e_val(i) = Ea(i) + e_add(i)
    end do 

    ! Calculate rate constants for atom(n)
    do i = 1,1  !nrate (just doing 1 bc other rates will not be changed by lateral int)
      do j = 1,sym
        if (atom(n)%rates(i,j) /= 0.0d0) then
          lnrate_val(i,j) = log(prefactor)-e_val(i)/(kBev*T)
          atom(n)%rates(i,j) = exp(lnrate_val(i,j))
        end if 
      end do 
    end do

  END SUBROUTINE rate_const

!======================================================================================================!
  SUBROUTINE init_rateconst(T)
!======================================================================================================!
    implicit none 
    real(dp), intent(in) :: T
    integer  :: i,j
    real(dp) :: area_uc,mass

    area_uc = (3.0D0*SQRT(3.0D0)/2.0D0)*(2.85D-10)**2  ! unit cell area
    mass = 4.9859D-26                                  ! gas species mass

    ! set rate constant values
    rate_val = 0.0d0
    rate_val(1) = exp(log(10.0D0**13)-Ea(1)/(kBev*T))
    !rate_val(2) = exp(log(10.0D0**13)-Ea(1)/(kBev*T))

    diffusion_val(1) = exp(log(10.0D0**13)-Edif(1)/(kBev*T))
    diffusion_val(2) = exp(log(10.0D0**13)-Edif(2)/(kBev*T))
    diffusion_val(3) = 0.0d0
    diffusion_val(4) = exp(log(10.0D0**13)-Edif(4)/(kBev*T))
     
    ad_rate = 0.0d0
    ad_rate(1) = (p_atm(1)*101325.0d0*area_uc)/sqrt(2.0d0*pi*mass*kB*T)
    ad_rate(2) = (p_atm(2)*101325.0d0*area_uc)/sqrt(2.0d0*pi*mass*kB*T)

    des_rate(1) = exp(log(10.0D0**13)-Edes(1)/(kBev*T))
    des_rate(2) = exp(log(10.0D0**13)-Edes(2)/(kBev*T))
    des_rate(3) = exp(log(10.0D0**13)-Edes(3)/(kBev*T))
    des_rate(4) = 0.0d0

    !============================== Reaction Rates ========================================!
    do i = 1,2
      do j = 1,3
        lateral_rate(j,i) = exp(log(10.0d0**13)-(Ea(1)+dble(j)*lateral_energy(i))/(kBeV*T))
      end do 
    end do 

    !============================== Diffusion Rates =======================================!
    do i = 1,atomType 
      do j = 1,3
        lateral_dif_rate(j,i) = exp(log(10.0D0**13)-(Edif(i)+dble(j)*lateral_energy(i))/(kBev*T))
      end do 
    end do 

  END SUBROUTINE init_rateconst
!======================================================================================================!
  SUBROUTINE adsorption(surface_dim,ad_val) 
!======================================================================================================!
    implicit none 
    integer, intent(in) :: surface_dim(2)
    integer, intent(out) :: ad_val
    integer :: i,j,x,y,x2,y2,choice,a,b,c,d
    real(dp) :: z(2),w,total_rate,norm_val

    ad_val = 0 ; total_rate = 0.0d0 ; choice = 0

    do i = 1,size(ad_rate)     ! Sum all adsorption rates
      total_rate = total_rate + ad_rate(i)
    end do 

    call random_number(z)       ! randomize position on surface
    x = int( surface_dim(1)*z(1) + bound(1) )
    y = int( surface_dim(2)*z(2) + bound(3) )
    
    if (surface(x,y)%atom_id == 0) then ! attempt position
      if (nsurface + 1 > natom) print*,'Too big'

      call random_number(w)       ! choose atom to attempt adsorption
      norm_val = w*total_rate 
      do i = 1,size(ad_rate)
        if (sum(ad_rate(1:i)) >= norm_val) then
          choice = i ; exit
        end if 
      end do 

      select case (choice)

      case (1)  ! A adsorption
        ad_val = 1
        atom(nsurface+1)%pos(1) = x ; atom(nsurface+1)%pos(2) = y
        atom(nsurface+1)%id = 1 ; surface(x,y)%atom_id = atom(nsurface+1)%id
        surface(x,y)%atom_num = nsurface+1 ; nsurface = nsurface + 1

      case (2)  ! B adsorption

        ad_val = 1
        atom(nsurface+1)%pos(1) = x ; atom(nsurface+1)%pos(2) = y
        atom(nsurface+1)%id = 2 ; surface(x,y)%atom_id = atom(nsurface+1)%id
        surface(x,y)%atom_num = nsurface+1 ; nsurface = nsurface + 1

      end select 

    else if (surface(x,y)%atom_id /= 0) then 
      ad_val = 0
  
    end if 

  end subroutine adsorption
!======================================================================================================!
END MODULE case_2
!======================================================================================================!
