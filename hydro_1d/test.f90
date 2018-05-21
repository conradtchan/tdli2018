program test
  implicit none

  call test_boundary
  call test_eos
  call test_primitives

end program test

subroutine result(string, pass)
  character(len=*),  intent(in)   :: string
  logical,            intent(in)  :: pass

  if (pass) then
    print*, string, '  [OK]'
  else
    print*, string, '  [FAIL]'
  endif

end subroutine result

subroutine test_boundary
  use mod_boundary, only:boundary
  implicit none

  integer, parameter :: nx = 50
  real    :: q(-1:nx+2)
  real    :: valuel, valuer
  integer :: i
  logical :: pass

  print*,'-- Testing: Boundary conditions (boundary.f90)'

  ! populate q
  do i = 1, nx
    q(i) = real(i)
  enddo

  ! test periodic BCs
  pass = .true.
  call boundary(q, nx, 1)
  if (q(nx+1) /= q(1))    pass = .false.
  if (q(nx+2) /= q(2))    pass = .false.
  if (q(0)    /= q(nx))   pass = .false.
  if (q(-1)   /= q(nx-1)) pass = .false.
  call result('Periodic', pass)

  ! test reflecting BCs, positive parity
  pass = .true.
  call boundary(q, nx, 2, parity=1)
  if (q(nx+1) /= q(nx))    pass = .false.
  if (q(nx+2) /= q(nx-1))  pass = .false.
  if (q(0)    /= q(1))     pass = .false.
  if (q(-1)   /= q(2))     pass = .false.
  call result('Reflecting+', pass)

  ! test reflecting BCs, negative parity
  pass = .true.
  call boundary(q, nx, 2, parity=-1)
  if (q(nx+1) /= -q(nx))    pass = .false.
  if (q(nx+2) /= -q(nx-1))  pass = .false.
  if (q(0)    /= -q(1))     pass = .false.
  if (q(-1)   /= -q(2))     pass = .false.
  call result('Reflecting-', pass)

  ! test reflecting BCs, specific value
  pass = .true.
  valuel = 1.23456
  valuer = 9.87654
  call boundary(q, nx, 3, valuel=valuel, valuer=valuer)
  if (q(nx+1) /= valuer)  pass = .false.
  if (q(nx+2) /= valuer)  pass = .false.
  if (q(0)    /= valuel)  pass = .false.
  if (q(-1)   /= valuel)  pass = .false.
  call result('Specific', pass)

end subroutine test_boundary

subroutine test_eos
  use mod_eos, only:eos

  integer, parameter :: nn = 3
  real :: rho(nn), eps(nn)
  real :: pre(nn), cs(nn)
  real :: pre_actual(nn), cs_actual(nn)

  real :: error

  print*,'-- Testing: Equation of state (eos.f90)'

  ! Pre-calculated values
  rho(1) = 1.0
  eps(1) = 1.0
  pre_actual(1) = 0.6666666666666667
  cs_actual(1) = 1.0540925533894598

  rho(2) = 2.0
  eps(2) = 3.0
  pre_actual(2) = 4.0
  cs_actual(2) = 1.8257418583505538

  rho(3) = 4.5
  eps(3) = 6.7
  pre_actual(3) = 20.100000000000005
  cs_actual(3) = 2.7284509239574839

  call eos(rho, eps, pre, cs, nn)

  ! Compare calculated values with actual values
  error = sum(abs(pre - pre_actual))
  call result('Pressure', error < epsilon(error))
  error = sum(abs(cs - cs_actual))
  call result('Sound speed', error < epsilon(error))

end subroutine test_eos

subroutine test_primitives
  use mod_primitives, only:primitives, conserved_var
  use mod_eos, only:eos

  integer, parameter :: nn = 2

  real :: u_con(1:nn,5)
  real :: u_prim(1:nn,7)
  real :: u_prim_new(1:nn,7)

  real :: error

  print*,'-- Testing: Conservative/Primitive variable conversion (primitives.f90)'

  ! Initialise primitives, convert into conservatives, and then convert back to check
  u_prim(1,1) = 0.1
  u_prim(1,2) = 0.2
  u_prim(1,3) = 0.3
  u_prim(1,4) = 0.4
  u_prim(1,5) = 0.5

  u_prim(2,1) = 0.8
  u_prim(2,2) = 0.9
  u_prim(2,3) = 1.0
  u_prim(2,4) = 1.1
  u_prim(2,5) = 1.2

  ! pressure and cs are determined by the density and internal energy
  call eos(u_prim(1:nn,1), u_prim(1:nn,5), u_prim(1:nn,6), u_prim(1:nn,7), nn)

  call conserved_var(u_prim, u_con, nn)
  call primitives(u_con, u_prim_new, nn)

  error = sum(abs(u_prim_new - u_prim))
  call result('conversion and back', error < 1.d-15)

end subroutine test_primitives
