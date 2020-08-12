!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  Contains the tools needed to run the kmc and LSODE simulations
!
!  init_values        : allocates arrays
!  dealloc_values     : deallocates arrays
!  sohr_sim           : runs kmc simulation with atom following
!  kmc_sim            : runs kmc simulation 
!  surface_generator  : randomly sets up initial surface
!  manual_generator   : can set up specific initial surface
!  splines            : sets up splines to find coverage and neighbors
!  run                : lsode simulation
!  pulse_run          : lsode simulation for pulsed sim
!  jacnit             : dummy initialization jacobian needed for lsode 
!
!  Author : Rob Wells             Date : 11/09/18
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!======================================================================================================!
MODULE kinetic_tools
!======================================================================================================!
  use case_2

  implicit none
  contains
!======================================================================================================!
  SUBROUTINE init_values()
!======================================================================================================!
    implicit none
    integer :: i,j,n 

    n = 5000    ! array size

    ! allocate the arrays based on initalized constants
    if (.not. allocated(rate_val)) allocate(rate_val(nrate))
    if (.not. allocated(init_num)) allocate(init_num(atomType))
    if (.not. allocated(initial)) allocate(initial(nmol))

    initial = 0 ; init_num = 0 

  END SUBROUTINE init_values
!======================================================================================================!
  SUBROUTINE init_parallel_values()
!======================================================================================================!
    implicit none
    integer :: i,j

    if (.not. allocated(atom)) allocate(atom(natom))
    if (.not. allocated(surface)) allocate(surface(bound(1):bound(2),bound(3):bound(4)))
    do i = bound(1),bound(2)
      do j = bound(3),bound(4)
        surface(i,j)%atom_id = 0 ; surface(i,j)%atom_num = 0
      end do 
    end do 

    do i = 1,natom
      if (.not. allocated(atom(i)%rates)) allocate(atom(i)%rates(nrate,sym))
    end do
  END SUBROUTINE init_parallel_values
!======================================================================================================!
  SUBROUTINE dealloc_values()
!======================================================================================================!
    implicit none 

    if (allocated(atom)) deallocate(atom)
    if (allocated(rate_val)) deallocate(rate_val)
    if (allocated(initial)) deallocate(initial)
    if (allocated(surface)) deallocate(surface)
    if (allocated(init_num)) deallocate(init_num)
    if (allocated(surface_save)) deallocate(surface_save)
    if (allocated(atom_save)) deallocate(atom_save)

  END SUBROUTINE dealloc_values
!======================================================================================================!
SUBROUTINE run_sim(times,atom_select,simType,input,limiter,sim)
!======================================================================================================!
  implicit none
  integer, intent(in) :: atom_select,simType,input(:),limiter,sim
  integer :: i,j,k,t,n,z,move,values(8),surface_dim(2),ntime,counts,ad_val,sizer,atom_choice,id_val
  integer, allocatable :: seed(:),species(:,:,:)
  real(dp), intent(in) :: times(3)
  real(dp) :: randy,val,val2,time,delta_time,total_rate,part_rate,sums,timeStop,finish,start,check_rate
  real(dp) :: pick_time,follow_time
  real(dp), allocatable :: surface_species_time(:),spline_saver(:,:)

  sizer   = 20000

  !========= Random number generator ==========!
  call date_and_time(values=values)
  call random_seed(size=k)
  if (.not. allocated(seed)) allocate(seed(k))
  seed(:) = values(8)
  call random_seed(put=seed)
  !============================================!

  if (.not. allocated(species)) allocate(species(sizer,5,atomType))
  if (.not. allocated(surface_species_time)) allocate(surface_species_time(sizer))
  if (.not. allocated(spline_saver)) allocate(spline_saver(spline_points,6))

  print*,'nsurface',nsurface,sim 
  print*,natom,nsurface

  !======================= Initialize simulation constants =======================! 
  ntime = 1 ; t = 1 ; counts = 1 
  surface_dim = shape(surface)

  time = 0.0d0 ; timeStop = times(2) ; follow_time = times(3)  

  ! these atoms are not on the surface
  do i = nsurface+1,natom
    atom(i)%pos = [gas_phase,gas_phase]
    atom(i)%id = 0
  end do 

  ! Set up neighbors, determine rates, & coverages
  do i = 1,nsurface
    if (atom(i)%id == 0) then
      atom(i)%rates = 0.0d0
    else 
      call generate_neighbors(i,1)      ! find neighboring molecules 
      call reaction_rates(i)            ! determines if reactions between atoms can occur
      if (all(lateral_energy /= 0.0d0)) call lateral(i)            ! lateral interactions
    end if
  end do 

  ! initial total rate
  total_rate = 0.0d0 
  do i = 1,nsurface
    total_rate = total_rate + sum(atom(i)%rates)  
  end do
  total_rate = total_rate + (sum(ad_rate))*area

  atom_choice = 0 ; atom_finish = 0 ; atom_follow = 0

  !==================== start the simulation ====================!
  do while (time < timeStop)
    
    t = t + 1 
    
    if (1 == 0) then
      check_rate = 0.0d0 
      do i = 1,nsurface
        check_rate = check_rate + sum(atom(i)%rates)  
      end do
      check_rate = check_rate + (sum(ad_rate))*area

      if (abs(check_rate - total_rate) > check_rate*0.05d0) then 
        print*,'RATES ARENT ADDING UP',check_rate,total_rate
      end if 
    end if 
    
    !==== Solve for time at which event occurs ====!
    334  call random_number(randy)
    delta_time = -1.0D0/total_rate * LOG(randy)
    time = time + delta_time

    !====== Choose which events occurs ======!
    call random_number(randy)
    val = total_rate*randy ! first: choose which atom, weighted by total rates

    part_rate = (sum(ad_rate))*area 

    ! Adsorption happens
    if (val < part_rate) then 

      call adsorption(surface_dim,ad_val)

      if (ad_val == 1) then

        n = nsurface ; atom(n)%neighbor = 0

        call generate_neighbors(n,0)  ! determines neighbors and sets up diffusion rates
        call reaction_rates(nsurface) ! sets up desorption and reaction rates
        if (all(lateral_energy /= 0.0d0)) call lateral(n)
        total_rate = total_rate + sum(atom(n)%rates) 

        do i = 1,sym 
          if (atom(n)%neighbor(i) /= 0) then
            total_rate = total_rate - sum(atom(atom(n)%neighbor(i))%rates) ! get rid of previous contribution to total rate
            call generate_neighbors(atom(n)%neighbor(i),0)
            call reaction_rates(atom(n)%neighbor(i))
            if (all(lateral_energy /= 0.0d0)) call lateral(atom(n)%neighbor(i))
            total_rate = total_rate + sum(atom(atom(n)%neighbor(i))%rates) ! add new contribution from adsorbed atom
          end if 
        end do

      else if (ad_val == 0) then
        go to 334 

      end if

    ! Something other than diffusion happens
    else 
      do i = 1,nsurface
        part_rate = part_rate + sum(atom(i)%rates)
        if (val <= part_rate) then
          n = i 
          exit  
        end if
      end do 

      call random_number(randy)
      val2 = sum(atom(n)%rates)*randy  ! second: choose which rate from atom n is chosen
      move = 0

      do i = 1,nrate
        sums = 0.0d0
        if (i>1) sums = sum(atom(n)%rates(1:i-1,:))
        do j = 1,sym
          sums = sums +atom(n)%rates(i,j)
          if (val2 <= sums) then
            move = sym*(i-1)+j
            exit
          end if 
        end do
        if (move > 0) exit 
      end do 

      ! Update nearby neighbors if diffusion will happen and remove old spot from surface
      if (move < 5) then
        total_rate = total_rate - sum(atom(n)%rates)  ! remove rate contribution from diffusing atom
        atom(n)%rates = 0.0d0
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_id = 0
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_num = 0
        do i = 1,sym 
          if (atom(n)%neighbor(i) /= 0) then
            total_rate = total_rate - sum(atom(atom(n)%neighbor(i))%rates) ! remove rate contribution from diffusing neighbor
            call generate_neighbors(atom(n)%neighbor(i),0)
            call reaction_rates(atom(n)%neighbor(i))
            if (all(lateral_energy /= 0.0d0)) call lateral(atom(n)%neighbor(i))
            total_rate = total_rate + sum(atom(atom(n)%neighbor(i))%rates) ! add rate from neighbor with diffusing atom gone
          end if 
        end do
      end if

      ! diffuse or react according to monte carlo choice
      if (move == 1) atom(n)%pos(1) = atom(n)%pos(1) + 1
      if (move == 2) atom(n)%pos(1) = atom(n)%pos(1) - 1
      if (move == 3) atom(n)%pos(2) = atom(n)%pos(2) + 1
      if (move == 4) atom(n)%pos(2) = atom(n)%pos(2) - 1
      if (move > 4) call reaction_results(n,total_rate,move) 

      ! Implements boundary conditions
      if (atom(n)%pos(1) == bound(1)-1) atom(n)%pos(1) = bound(2)
      if (atom(n)%pos(1) == bound(2)+1) atom(n)%pos(1) = bound(1)
      if (atom(n)%pos(2) == bound(3)-1) atom(n)%pos(2) = bound(4)
      if (atom(n)%pos(2) == bound(4)+1) atom(n)%pos(2) = bound(3)

      ! Update current atom and neighbors
      if (atom(n)%pos(1) /= gas_phase) then 
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_id = atom(n)%id 
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_num = n 
        total_rate = total_rate - sum(atom(n)%rates)      
        call generate_neighbors(n,0)
        call reaction_rates(n)
        if (all(lateral_energy /= 0.0d0)) call lateral(n)
        total_rate = total_rate + sum(atom(n)%rates)
      end if 
      do i = 1,sym 
        if (atom(n)%neighbor(i) /= 0) then
          total_rate = total_rate - sum(atom(atom(n)%neighbor(i))%rates)
          call generate_neighbors(atom(n)%neighbor(i),0)
          call reaction_rates(atom(n)%neighbor(i))
          if (all(lateral_energy /= 0.0d0)) call lateral(atom(n)%neighbor(i))
          total_rate = total_rate + sum(atom(atom(n)%neighbor(i))%rates)
        end if 
      end do 

    end if ! adsorption

    ! write coverages
    if (time <= timeStop) then
      if (ntime == limiter) then
        if (counts <= sizer) then  
          call surface_species_finder(time,counts,species,surface_species_time,input(2))
          counts = counts + 1
        end if  

        ! print out surface if specified
        if (simType == 1) then
          do j = bound(4),bound(3),-1
            write(input(3),*) (surface(i,j)%atom_id,i=bound(1),bound(2))
          end do
          write(input(4),*) time
        end if
        ntime =  1    
      else 
        ntime = ntime + 1
      end if  ! limiter
    end if 

    if ( SUM(surface(bound(1):bound(2),bound(3):bound(4))%atom_id) /= SUM(atom%id) ) then
      write(*,*) 'atoms disappearing',SUM(atom%id),SUM(surface(bound(1):bound(2),bound(3):bound(4))%atom_id),time,nsurface
    end if

  end do  ! time

  print*,'Species Time: ',time-pick_time,time,pick_time,counts
  !print*,'Species Id: ',id_val
  write(input(1),*) time-pick_time,id_val,time,pick_time

  ! Create spline to average coverages and neighbors from all sims
  do i = 1,atomType-1 
    call splines(counts,species(:,:,i),surface_species_time,[times(1),times(2)],input(4+i),spline_saver)
    saves(:,:,i,sim) = spline_saver
  end do 

END SUBROUTINE run_sim  
!======================================================================================================!
SUBROUTINE sohr_sim(times,atom_select,simType,input,limiter,sim)
!======================================================================================================!
  implicit none
  integer, intent(in) :: atom_select,simType,input(:),limiter,sim
  integer :: i,j,k,t,n,z,move,values(8),surface_dim(2),ntime,counts,ad_val,sizer,atom_choice,id_val
  integer, allocatable :: seed(:),species_B(:,:),species_A(:,:),species(:,:,:)
  real(dp), intent(in) :: times(3)
  real(dp) :: randy,val,val2,time,delta_time,total_rate,part_rate,sums,timeStop,finish,start,check_rate
  real(dp) :: pick_time,follow_time
  real(dp), allocatable :: surface_species_time(:),spline_saver(:,:)

  sizer   = 20000

  !========= Random number generator ==========!
  call date_and_time(values=values)
  call random_seed(size=k)
  if (.not. allocated(seed)) allocate(seed(k))
  seed(:) = values(8)
  call random_seed(put=seed)
  !============================================!

  if (.not. allocated(species)) allocate(species(sizer,5,atomType))
  if (.not. allocated(surface_species_time)) allocate(surface_species_time(sizer))
  if (.not. allocated(spline_saver)) allocate(spline_saver(spline_points,6))

  !======================= Initialize simulation constants =======================! 
  ntime = 1 ; t = 1 ; counts = 1 
  surface_dim = shape(surface)

  time = 0.0d0 ; follow_time = times(3) ; timeStop = times(2)

  ! these atoms are not on the surface
  do i = nsurface+1,natom
    atom(i)%pos = [gas_phase,gas_phase]
    atom(i)%id = 0
  end do 

  ! Set up neighbors, determine rates, & coverages
  do i = 1,nsurface
    if (atom(i)%id == 0) then
      atom(i)%rates = 0.0d0
    else 
      call generate_neighbors(i,1)      ! find neighboring molecules 
      call reaction_rates(i)            ! determines if reactions between atoms can occur
      if (all(lateral_energy /= 0.0d0)) call lateral(i)            ! lateral interactions
    end if
  end do 

  ! initial total rate
  total_rate = 0.0d0 
  do i = 1,nsurface
    total_rate = total_rate + sum(atom(i)%rates)  
  end do
  total_rate = total_rate + (sum(ad_rate))*area

  atom_choice = 0 ; atom_finish = 0 ; atom_follow = 0

  !==================== start the simulation ====================!
  do while (atom_finish == 0)
    
    t = t + 1
    
    if (1 == 0) then
      check_rate = 0.0d0 
      do i = 1,nsurface
        check_rate = check_rate + sum(atom(i)%rates)  
      end do
      check_rate = check_rate + (sum(ad_rate))*area

      if (abs(check_rate - total_rate) > check_rate*0.05d0) then 
        print*,'RATES ARENT ADDING UP',check_rate,total_rate
      end if 
    end if 

    !====== Choose atom to follow if follow time is at begin ======!
    if (follow_time == 0.0d0 .AND. atom_choice == 0) then
      do while (atom_choice == 0)
        call random_number(randy)
        atom_follow = int(randy * nsurface) + 1
        if (atom(atom_follow)%id == atom_select ) atom_choice = atom_select
      end do
      pick_time = time 
    end if 
    
    !==== Solve for time at which event occurs ====!
    334  call random_number(randy)
    delta_time = -1.0D0/total_rate * LOG(randy)
    time = time + delta_time

    !====== Choose which events occurs ======!
    call random_number(randy)
    val = total_rate*randy ! first: choose which atom, weighted by total rates

    part_rate = (sum(ad_rate))*area 

    ! Adsorption happens
    if (val < part_rate) then 

      call adsorption(surface_dim,ad_val)

      if (ad_val == 1) then

        !==== Check to see if atom adsorbed is right atom, follow if so ====!
        if (time >= follow_time .AND. atom_choice == 0) then
          if (atom(n)%id == atom_select) then
            atom_choice = atom_select ; pick_time = time
            atom_follow = n
          end if 
        end if 

        n = nsurface ; atom(n)%neighbor = 0

        call generate_neighbors(n,0)  ! determines neighbors and sets up diffusion rates
        call reaction_rates(nsurface) ! sets up desorption and reaction rates
        if (all(lateral_energy /= 0.0d0)) call lateral(n)
        total_rate = total_rate + sum(atom(n)%rates) 

        do i = 1,sym 
          if (atom(n)%neighbor(i) /= 0) then
            total_rate = total_rate - sum(atom(atom(n)%neighbor(i))%rates) ! get rid of previous contribution to total rate
            call generate_neighbors(atom(n)%neighbor(i),0)
            call reaction_rates(atom(n)%neighbor(i))
            if (all(lateral_energy /= 0.0d0)) call lateral(atom(n)%neighbor(i))
            total_rate = total_rate + sum(atom(atom(n)%neighbor(i))%rates) ! add new contribution from adsorbed atom
          end if 
        end do

      else if (ad_val == 0) then
        go to 334 

      end if

    ! Something other than diffusion happens
    else 
      do i = 1,nsurface
        part_rate = part_rate + sum(atom(i)%rates)
        if (val <= part_rate) then
          n = i 
          exit  
        end if
      end do 

      call random_number(randy)
      val2 = sum(atom(n)%rates)*randy  ! second: choose which rate from atom n is chosen
      move = 0

      do i = 1,nrate
        sums = 0.0d0
        if (i>1) sums = sum(atom(n)%rates(1:i-1,:))
        do j = 1,sym
          sums = sums +atom(n)%rates(i,j)
          if (val2 <= sums) then
            move = sym*(i-1)+j
            exit
          end if 
        end do
        if (move > 0) exit 
      end do 

      ! Update nearby neighbors if diffusion will happen and remove old spot from surface
      if (move < 5) then
        total_rate = total_rate - sum(atom(n)%rates)  ! remove rate contribution from diffusing atom
        atom(n)%rates = 0.0d0
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_id = 0
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_num = 0
        do i = 1,sym 
          if (atom(n)%neighbor(i) /= 0) then
            total_rate = total_rate - sum(atom(atom(n)%neighbor(i))%rates) ! remove rate contribution from diffusing neighbor
            call generate_neighbors(atom(n)%neighbor(i),0)
            call reaction_rates(atom(n)%neighbor(i))
            if (all(lateral_energy /= 0.0d0)) call lateral(atom(n)%neighbor(i))
            total_rate = total_rate + sum(atom(atom(n)%neighbor(i))%rates) ! add rate from neighbor with diffusing atom gone
          end if 
        end do
      end if

      ! diffuse or react according to monte carlo choice
      if (move == 1) atom(n)%pos(1) = atom(n)%pos(1) + 1
      if (move == 2) atom(n)%pos(1) = atom(n)%pos(1) - 1
      if (move == 3) atom(n)%pos(2) = atom(n)%pos(2) + 1
      if (move == 4) atom(n)%pos(2) = atom(n)%pos(2) - 1
      if (move > 4) call reaction_results_sohr(n,total_rate,move,id_val) 

      ! Implements boundary conditions
      if (atom(n)%pos(1) == bound(1)-1) atom(n)%pos(1) = bound(2)
      if (atom(n)%pos(1) == bound(2)+1) atom(n)%pos(1) = bound(1)
      if (atom(n)%pos(2) == bound(3)-1) atom(n)%pos(2) = bound(4)
      if (atom(n)%pos(2) == bound(4)+1) atom(n)%pos(2) = bound(3)

      ! Update current atom and neighbors
      if (atom(n)%pos(1) /= gas_phase) then 
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_id = atom(n)%id 
        surface(atom(n)%pos(1),atom(n)%pos(2))%atom_num = n 
        total_rate = total_rate - sum(atom(n)%rates)      
        call generate_neighbors(n,0)
        call reaction_rates(n)
        if (all(lateral_energy /= 0.0d0)) call lateral(n)
        total_rate = total_rate + sum(atom(n)%rates)
      end if 
      do i = 1,sym 
        if (atom(n)%neighbor(i) /= 0) then
          total_rate = total_rate - sum(atom(atom(n)%neighbor(i))%rates)
          call generate_neighbors(atom(n)%neighbor(i),0)
          call reaction_rates(atom(n)%neighbor(i))
          if (all(lateral_energy /= 0.0d0)) call lateral(atom(n)%neighbor(i))
          total_rate = total_rate + sum(atom(atom(n)%neighbor(i))%rates)
        end if 
      end do 

    end if ! adsorption

    ! write coverages
    if (time <= timeStop) then
      if (ntime == limiter) then
        if (counts <= sizer) then  
          call surface_species_finder(time,counts,species,surface_species_time,input(2))
          counts = counts + 1
        end if  

        ! print out surface if specified
        if (simType == 1) then
          do j = bound(4),bound(3),-1
            write(input(3),*) (surface(i,j)%atom_id,i=bound(1),bound(2))
          end do
          write(input(4),*) time
        end if
        ntime =  1    
      else 
        ntime = ntime + 1
      end if  ! limiter
    end if 

    if ( SUM(surface(bound(1):bound(2),bound(3):bound(4))%atom_id) /= SUM(atom%id) ) then
      write(*,*) 'atoms disappearing',SUM(atom%id),SUM(surface(bound(1):bound(2),bound(3):bound(4))%atom_id),time,nsurface
    end if

  end do  ! time

  print*,'Species Time: ',time-pick_time,time,pick_time,counts
  !print*,'Species Id: ',id_val
  write(input(1),*) time-pick_time,id_val,time,pick_time

  ! Create spline to average coverages and neighbors from all sims
  do i = 1,atomType-1 
    call splines(counts,species(:,:,i),surface_species_time,[times(1),times(3)],input(4+i),spline_saver)
    saves(:,:,i,sim) = spline_saver
  end do  

END SUBROUTINE sohr_sim
!======================================================================================================!
  SUBROUTINE surface_generator(surface_dim,type)
!======================================================================================================!
    implicit none 
    integer, intent(in) :: surface_dim(2),type
    integer :: i,j,k,x,y,values(8),sums,exit_val
    integer, allocatable :: seed(:)
    real(dp) :: z(2)

    if (.not. allocated(surface_save)) allocate(surface_save(bound(1):bound(2),bound(3):bound(4)))
    if (.not. allocated(atom_save)) allocate(atom_save(natom))

    !======== Random number generator ==========!
    call date_and_time(values=values)
    call random_seed(size=k)
    if (.not. allocated(seed)) allocate(seed(k))
    seed(:) = values(8)
    call random_seed(put=seed)
    !===========================================!

    ! initialize positions of atoms on the surface
    if (type == 1) then ! start a new surface
      do j = bound(1),bound(2)
        do k = bound(3),bound(4)
          surface(j,k)%atom_id = 0 ; surface(j,k)%atom_num = 0
        end do 
      end do  
      i = 1 
    else ! use existing surface and build from there
      i = type
    end if 

    do while (i <= nsurface) 

      sums = 0 ; j = 0 ; exit_val = 0
      do while (exit_val == 0) 
        j = j+1 ; sums = sums + init_num(j) 
        if (i <= sums) then 
          atom(i)%id = j ; exit_val = 1
        else if (j > atomType) then 
          print*,'j too big' ; exit_val = 1
        end if        
      end do 
        
      call random_number(z)        ! Randomize position on surface
      x = int( surface_dim(1)*z(1) + bound(1) )
      y = int( surface_dim(2)*z(2) + bound(3) )

      if (surface(x,y)%atom_id == 0 ) then
        atom(i)%pos(1) = x ; atom(i)%pos(2) = y
        surface(x,y)%atom_id = atom(i)%id 
        surface(x,y)%atom_num = i
        i = i+1
      end if

    end do

    write(*,*) 'surface generator'

    surface_save = surface   ; atom_save = atom 
    nsurface_save = nsurface ; natom_save = natom

  END SUBROUTINE surface_generator
!======================================================================================================!
  SUBROUTINE surface_species_finder(time,counts,species,surface_species_time,fils)
!======================================================================================================!
    implicit none
    integer, intent(in)     :: fils,counts
    integer, intent(inout)  :: species(:,:,:)
    real(dp), intent(in)    :: time
    real(dp), intent(inout) :: surface_species_time(:)
    integer :: i,j,x,y,k,a,b,c,d,num_B,num_A,sums,species2(5,2) 

    num_B = 0 ; num_A = 0 ; species2 = 0

    do x = 1,bound(2)
      do y = 1,bound(4)

        if (surface(x,y)%atom_id == 1) then

          num_A = num_A + 1 ; sums = 0

          a = x+1 ; b = x-1 ; c = y+1 ; d = y-1

          ! Incorporate boundary conditions for nearest site
          if (a > bound(2)) a = bound(1)
          if (b < bound(1)) b = bound(2)
          if (c > bound(4)) c = bound(3)
          if (d < bound(3)) d = bound(4) 

          if (surface(a,y)%atom_id == 1) sums = sums + 1
          if (surface(b,y)%atom_id == 1) sums = sums + 1
          if (surface(x,c)%atom_id == 1) sums = sums + 1
          if (surface(x,d)%atom_id == 1) sums = sums + 1

          species2(sums+1,1) = species2(sums+1,1) + 1

        else if (surface(x,y)%atom_id == 2) then 

          num_B = num_B + 1 ; sums = 0

          a = x+1 ; b = x-1 ; c = y+1 ; d = y-1

          ! Incorporate boundary conditions for nearest site
          if (a > bound(2)) a = bound(1)
          if (b < bound(1)) b = bound(2)
          if (c > bound(4)) c = bound(3)
          if (d < bound(3)) d = bound(4) 

          if (surface(a,y)%atom_id == 2) sums = sums + 1
          if (surface(b,y)%atom_id == 2) sums = sums + 1
          if (surface(x,c)%atom_id == 2) sums = sums + 1
          if (surface(x,d)%atom_id == 2) sums = sums + 1

          species2(sums+1,2) = species2(sums+1,2) + 1

        end if 

      end do 
    end do    

    do i = 1,5
      species(counts,i,1) = species2(i,1)
      species(counts,i,2) = species2(i,2)
    end do
    surface_species_time(counts) = time  

    write(fils,*) time,num_A,num_B,((species2(i,j),i=1,5),j=1,2)

  END SUBROUTINE
!======================================================================================================!
  SUBROUTINE splines(counts,species,species_time,times,fil,spline_saver)
!======================================================================================================!
    implicit none
    integer, intent(in) :: counts,species(:,:),fil
    integer :: i,j,k,n,val 
    real(dp), intent(in) :: species_time(:),times(2)
    real(dp), allocatable, intent(inout) :: spline_saver(:,:)
    real(dp) :: yval,ypval,yppval,yval2,ypval2,yppval2,spline_delt
    real(dp), allocatable :: spline_x(:),spline_y(:,:),spliner(:,:),spline_time(:),spline_val(:,:)

    n = counts-1

    if (.not. allocated(spliner)) allocate(spliner(n,5))
    if (.not. allocated(spline_x)) allocate(spline_x(n))
    if (.not. allocated(spline_y)) allocate(spline_y(n,5))
    if (.not. allocated(spline_val)) allocate(spline_val(spline_points,5))
    if (.not. allocated(spline_time)) allocate(spline_time(spline_points))

    spliner = 0.0d0    ; spline_x = 0.0d0     ; spline_y = 0.0d0
    spline_val = 0.0d0 ; spline_time = 0.0d0  

    spline_saver = 0.0d0

    do i = 1,n
      do j= 1,5 
        spline_y(i,j) = real(species(i,j),dp)
      end do 
      spline_x(i) = species_time(i)
    end do  

    ! set up spline
    do j = 1,5
      call spline_cubic_set(n,spline_x,spline_y(:,j),0,spline_y(1,j),0,spline_y(n,j),spliner(:,j))
    end do 

    spline_delt = (times(2)-times(1))/dble(spline_points)
    spline_time(1) = spline_delt

    ! solve at times given 
    do i = 2,spline_points
      spline_time(i) = spline_time(i-1) + spline_delt
      do j = 1,5 
        call spline_cubic_val(n,spline_x,spline_y(:,j),spliner(:,j),spline_time(i),yval,ypval,yppval)
        if (yval < 0.0d0) yval = 0.0d0
        spline_val(i,j) = yval
      end do
    end do

    do i = 1,spline_points
      spline_saver(i,:) = [spline_time(i),(spline_val(i,j),j=1,5)] 
    end do 

    deallocate(spliner)
    deallocate(spline_x)
    deallocate(spline_y)
    deallocate(spline_val)
    deallocate(spline_time) 

  END SUBROUTINE splines
!======================================================================================================!
  SUBROUTINE splines_average(numSim,nSpecies,fil,spline_saver,atom_val)
!======================================================================================================!
    implicit none 
    integer, intent(in)   :: numSim,nSpecies,fil,atom_val
    integer :: i,j,k,l
    real(dp), intent(in)  :: spline_saver(:,:,:,:)
    real(dp), allocatable :: species(:,:),species_avg(:,:),species_total(:,:),time(:) 

    if (.not. allocated(species)) allocate(species(spline_points,nSpecies))
    if (.not. allocated(species_avg)) allocate(species_avg(spline_points,nSpecies))
    if (.not. allocated(species_total)) allocate(species_total(spline_points,nSpecies))
    if (.not. allocated(time)) allocate(time(spline_points))

    species = 0.0d0 ; species_avg = 0.0d0 ; species_total = 0.0d0

    do k = 1,nSpecies
      do l = 1,spline_points
        do j = 1,numSim  
          species_total(l,k) = species_total(l,k) + spline_saver(l,k+1,atom_val,j)
        end do 
      end do 
    end do  

    do k = 1,spline_points 
      do j = 1,nspecies
        species_avg(k,j) = species_total(k,j) / real(numSim)
      end do 
      write(fil,*) spline_saver(k,1,atom_val,1),(species_avg(k,j),j=1,nspecies)
    end do 

  END SUBROUTINE splines_average
!======================================================================================================!
END MODULE kinetic_tools
!======================================================================================================!
