program hydro
  use mod_primitives, only:primitives,conserved_var
  use mod_rk_step, only:step
  use mod_output, only:write_output
  use mod_setup, only:setup

  implicit none

  real,    parameter  :: length = 1.0
  real,    parameter  :: tmax = 0.5
  integer, parameter  :: nx = 1000
  integer, parameter  :: output_interval = 10

  real :: u_con(1:nx,5)
  real :: u_prim(1:nx,7)

  real :: x(1:nx)
  real :: x_if(0:nx)

  real :: time

  integer :: i, nstep

  print*,'-- Starting hydro code'

  ! Starting time
  time = 0.0

  ! Calculate interface coordinates
  do i = 0, nx
    x_if(i) = real(i) * length / real(nx)
  enddo

  ! Calculate cell center coordinates (using interface coordinates)
  do i = 1, nx
    x(i) = 0.5 * (x_if(i-1) + x_if(i))
  enddo

  ! Set up initial conditions
  print*, '-- Setting up initial conditions'
  call setup(x, u_prim, nx)

  ! Get conserved quantities
  call conserved_var(u_prim, u_con, nx)

  ! Main time loop
  print*, '-- Integrating to t = ', tmax
  nstep = 0
  do while (time < tmax)
    call step(time, u_con, u_prim, x, x_if, nx)

    ! Occasionally write outputs
    if (mod(nstep, output_interval) == 0) then
      if (mod(nstep, output_interval*50) == 0) then
        write(*,*) ''
        write(*, fmt='(a, f8.4, a)', advance='no') 'time = ', time, ' '
      endif
      write(*, fmt='(a)', advance='no') '.'

      ! Need to convert to primitives for output
      call primitives(u_con, u_prim, nx)
      call write_output(nstep, x, u_prim, nx)
    endif
    nstep = nstep + 1
  enddo
  print*, ''
  print*, ''
  print*, '-- Done!'

end program hydro
