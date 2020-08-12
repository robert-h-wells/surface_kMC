!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!  This program runs a kinetic monte carlo for reactions on surfaces. Variable step size algorithm.
!  Includes diffusion and lateral interactions. 
!
!  Natom is total possible atoms on surface for the simulation. Nsurface is the number of atoms on 
!  surface at the moment (needed for when reactions create extra atoms, like decomposition)    
!  
!  Author : Rob Wells           Date : 11/09/18
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

!======================================================================================================!
PROGRAM main
!======================================================================================================!
  use omp_lib
  use case_1   
  use kinetic_tools

  implicit none
  integer  :: i,j,k,ii,kmc_files(6),lsode_files(2),simType(3),initial_sum,sim,val
  integer  :: timeNum,simrun,numSim,limiter
  real(dp) :: finish,start,timeStop,times(4),T,pulse_times(6),follow_time(6),start1
  real(dp) :: pressures(3),count_rate,count,sohr_time(3)
  character(len=40) :: names1,names2,names3,names4,names5,names6

  CALL cpu_time(start1)
  990 format(a6,f9.2,a8)
  600 format(a15,4ES12.3)

  timer = 0.0d0

  ! initialize simulation times
  timeStop = 1.0d-05
  times(1:2) = [0.0d0,timeStop]      ! kmc start and end time
  times(3:4) = [1.0d-11,timeStop]    ! LSODE step size and end time 

  follow_time = [1.0d-05,2.0d-07,3.0d-07,4.0d-07,5.0d-07,6.0d-07]
  timeNum = size(follow_time)

  simrun = 1

  !================================ Kinetic Monte Carlo ================================!

  ! Initialize simulation constants
  numSim = 10               ! # of simulations
  spline_points = 600
  bound = [1,50,1,50]         ! size of surface (1-min x, 2-max x, 3-min y, 4-max y)
  area = bound(2)*bound(4)
  natom = bound(2)*bound(4)     ! # of maximum atoms on surface
  simType = [0,1,1]             ! (1) = 1 print surface, (2) = 1 create new surface
                                ! (3) = 1 normal / = 2 pulse / 3 = sohr
 
  T = 500.0d0                   ! temperature of system

  pressures = [1.0d-02,1.0d-02,1.050d-03] 
  p_atm(1) = pressures(1) ; p_atm(2) = pressures(2)

  if (.not. allocated(saves)) allocate(saves(spline_points,6,atomType,numSim))

  ! allocate arrays and set initial rate constant values
  call init_case()
  call init_values()

!$omp parallel   
  call init_parallel_values()
!$omp end parallel 

  call init_rateconst(T)

  991 format(a8,ES12.3)
  write(*,991) 'Ratea',rate_val(1)
  write(*,991) 'RateD',diffusion_val(1)

  nsurface = 0
  init_num(1) = 0 ; init_num(2) = 0
  call surface_generator(shape(surface),1)
  surface = surface_save ; atom = atom_save

  do i = 1,1 !4 !timeNum
    
    write(names1,"(A5,I0.2,A4)" ) "sohr_",i,".out"
    write(names2,"(A16,I0.2,A4)" ) "surface_species_",i,".out"
    write(names5,"(A12,I0.2,A4)" ) "spline_avgA_",i,".out"
    write(names6,"(A12,I0.2,A4)" ) "spline_avgB_",i,".out"
    open(unit=200+i,file=names1) ; open(unit=210+i,file=names2)
    open(unit=260+i,file=names5) ; open(unit=270+i,file=names6)

    if (simType(1) == 1) then 
      write(names3,"(A8,I0.2,A4)" ) "surface_",i,".out"
      write(names4,"(A13,I0.2,A4)" ) "surface_time_",i,".out"
      open(unit=220+i,file=names3) ; open(unit=230+i,file=names4)
    end if 

    kmc_files = [200+i,210+i,220+i,230+i,240+i,250+i]
    sohr_time = [times(1),times(2),follow_time(i)]

    limiter = timeStop/(1e-05)*100
    print*,'limiter',limiter

!$omp parallel do private(sim) shared(surface_save,atom_save,natom_save,nsurface_save)
! simType,kmc_files,sohr_time,

    do sim = 1,numSim
      print*,'Sim ::',sim,i,kmc_files(4:),limiter

      surface = surface_save   ; atom = atom_save 
      nsurface = nsurface_save ; natom = natom_save

      if (simType(3) == 1) call run_sim(sohr_time,1,simType(1),kmc_files,limiter,sim)  
      if (simType(3) == 3) call sohr_sim(sohr_time,1,simType(1),kmc_files,limiter,sim)

    end do ! sim
    
!$omp end parallel do

    call splines_average(numSim,5,260+i,saves,1)
    call splines_average(numSim,5,270+i,saves,2)

  end do ! time

  ! Deallocate arrays 
  call dealloc_values()

  call cpu_time(finish)
  write(*,990) 'Total Time:',finish-start1,'seconds'

END PROGRAM main
!======================================================================================================!
