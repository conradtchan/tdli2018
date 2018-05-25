module mod_setup
  implicit none

contains

  subroutine setup(x, u_prim, nx)
    use mod_eos, only:eos,gamma
    real,     intent(in)  :: x(1:nx)
    real,     intent(out) :: u_prim(1:nx,7)
    integer,  intent(in)  :: nx

    real :: rho(1:nx)
    real :: v1(1:nx)
    real :: v2(1:nx)
    real :: v3(1:nx)
    real :: eps(1:nx)
    real :: pre(1:nx)

    ! Shock tube initial conditions
    rho(1:nx/2)    = 1.0
    rho(nx/2+1:nx) = 0.1
    pre(1:nx/2)    = 1.0
    pre(nx/2+1:nx) = 0.1
    v1(1:nx/2)     = 0.0
    v1(nx/2+1:nx)  = 0.0
    v2(:) = 0.0
    v3(:) = 0.0

    ! Calculate eps
    eps(1:nx) = pre(1:nx) / rho(1:nx) / (gamma - 1.0)

    u_prim(1:nx,1) = rho (1:nx)
    u_prim(1:nx,2) = v1  (1:nx)
    u_prim(1:nx,3) = v2  (1:nx)
    u_prim(1:nx,4) = v3  (1:nx)
    u_prim(1:nx,5) = eps (1:nx)

    ! Use EOS to get pressure and sound speed
    call eos(u_prim(1:nx,1), u_prim(1:nx,5), u_prim(1:nx,6), u_prim(1:nx,7), nx)

  end subroutine setup

end module mod_setup
