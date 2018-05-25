program hydro
  use mod_primitives, only:primitives,conserved_var
  use mod_rk_step, only:step
  use mod_output, only:write_output
  use mod_setup, only:setup
  use mod_grid, only:grid

  implicit none

  real,    parameter  :: length = 1.0
  real,    parameter  :: tmax = 20.0
  integer, parameter  :: nx = 512
  integer, parameter  :: ny = 512
  integer, parameter  :: output_interval = 20

  real, allocatable :: u_con(:,:,:)
  real, allocatable :: u_prim(:,:,:)

  real :: x(1:nx), y(1:ny)
  real :: x_if(0:nx), y_if(0:ny)

  real :: time

  integer :: j, nstep

  print*,'-- Starting hydro code'

  print*,'-- Allocating memory'
  allocate(u_con(1:nx,1:ny,5))
  allocate(u_prim(1:nx,1:ny,7))

  ! Initialise grid coordinates
  call grid(x_if, x, y_if, y, length, length, nx, ny)

  ! Starting time
  time = 0.0

  ! Set up initial conditions
  print*, '-- Setting up initial conditions'
  call setup(x, y, u_prim, length, length, nx, ny)

  ! Get conserved quantities
  !$omp parallel do
  do j = 1, ny
    call conserved_var(u_prim(:,j,:), u_con(:,j,:), nx)
  enddo
  !$omp end parallel do

  ! Main time loop
  print*, '-- Integrating to t = ', tmax
  nstep = 0
  do while (time < tmax)
    call step(time, u_con, u_prim, x, y, x_if, y_if, nx, ny)

    ! Occasionally write outputs
    write(*, fmt='(a)', advance='no') '.'
    if (mod(nstep, output_interval) == 0) then
      ! if (mod(nstep, output_interval*50) == 0) then
        write(*,*) ''
        write(*, fmt='(a, f8.4, a)', advance='no') 'time = ', time, ' '
      ! endif
      ! write(*, fmt='(a)', advance='no') '.'

      ! Need to convert to primitives for output
      !$omp parallel do
      do j = 1, ny
        call primitives(u_con(:,j,:), u_prim(:,j,:), nx)
      enddo
      !$omp end parallel do
      call write_output(nstep, x, y, u_prim, nx, ny)
    endif
    nstep = nstep + 1
  enddo
  print*, ''
  print*, ''
  print*, '-- Done!'

end program hydro
