module mod_sweep

  implicit none

contains

  subroutine sweep(u_prim, du_con_dt, x_in, x_if_in, nx)
    use mod_boundary, only:boundary
    use mod_eos, only:eos
    use mod_reconstruction, only:recon,coeff
    use mod_riemann, only:riemann

    real,     intent(in)    :: u_prim(1:nx,7)
    real,     intent(out)   :: du_con_dt(1:nx,5)
    real,     intent(in)    :: x_in(1:nx), x_if_in(0:nx)
    integer,  intent(in)    :: nx

    real  :: rho(-1:nx+2)
    real  :: u  (-1:nx+2)
    real  :: ut (-1:nx+2)
    real  :: utt(-1:nx+2)
    real  :: eps(-1:nx+2)
    real  :: pre(-1:nx+2)
    real  :: cs (-1:nx+2)

    real  :: x   (-1:nx+2)
    real  :: x_if(-2:nx+2)

    real  :: dx(0:nx+1)
    real  :: dx_if_inv(-1:nx+1)

    real  :: rho_if(0:nx+1,2)
    real  :: u_if  (0:nx+1,2)
    real  :: ut_if (0:nx+1,2)
    real  :: utt_if(0:nx+1,2)
    real  :: eps_if(0:nx+1,2)
    real  :: pre_if(0:nx+1,2)
    real  :: cs_if (0:nx+1,2)

    real :: rhoflx  (0:nx)
    real :: momflx  (0:nx)
    real :: momtflx (0:nx)
    real :: momttflx(0:nx)
    real :: eneflx  (0:nx)

    integer :: i

    ! Copy primitive variables
    rho(1:nx) = u_prim(1:nx,1)
    ut (1:nx) = u_prim(1:nx,3)
    u  (1:nx) = u_prim(1:nx,2)
    utt(1:nx) = u_prim(1:nx,4)
    eps(1:nx) = u_prim(1:nx,5)
    pre(1:nx) = u_prim(1:nx,6)
    cs (1:nx) = u_prim(1:nx,7)

    ! Boundary conditions
    call boundary(rho(-1:nx+2), nx)
    call boundary(u  (-1:nx+2), nx)
    call boundary(ut (-1:nx+2), nx)
    call boundary(utt(-1:nx+2), nx)
    call boundary(eps(-1:nx+2), nx)
    call boundary(pre(-1:nx+2), nx)
    call boundary(cs (-1:nx+2), nx)

    ! Fill coordinate array from input
    x(1:nx) = x_in(1:nx)
    x_if(0:nx) = x_if_in(0:nx)

    ! Ghost zones, get offset using neighbouring differences
    x(0)  = x_if(0) - (x(1) - x_if(0))
    x(-1) = x_if(0) - (x(2) - x_if(0))
    x_if(-1) = x_if(0) - (x_if(1) - x_if(0))
    x_if(-2) = x_if(0) - (x_if(2) - x_if(0))

    x(nx+1) = x_if(nx) - (x(nx) - x_if(nx))
    x(nx+2) = x_if(nx) - (x(nx-1) - x_if(nx))
    x_if(nx+1) = x_if(nx) - (x_if(nx-1) - x_if(nx))
    x_if(nx+2) = x_if(nx) - (x_if(nx-2) - x_if(nx))

    ! Calculate the dx coefficients for reconstruction
    call coeff(x(-1:nx+2), x_if(-2:nx+2), dx(0:nx+1), dx_if_inv(-1:nx+1), nx)

    ! Do reconstruction on all the primitive variables
    call recon(rho(-1:nx+2),  rho_if(0:nx+1,:), dx(0:nx+1), dx_if_inv(-1:nx+1), nx)
    call recon(u(-1:nx+2),    u_if(0:nx+1,:),   dx(0:nx+1), dx_if_inv(-1:nx+1), nx)
    call recon(ut(-1:nx+2),   ut_if(0:nx+1,:),  dx(0:nx+1), dx_if_inv(-1:nx+1), nx)
    call recon(utt(-1:nx+2),  utt_if(0:nx+1,:), dx(0:nx+1), dx_if_inv(-1:nx+1), nx)
    call recon(eps(-1:nx+2),  eps_if(0:nx+1,:), dx(0:nx+1), dx_if_inv(-1:nx+1), nx)
    call recon(pre(-1:nx+2),  pre_if(0:nx+1,:), dx(0:nx+1), dx_if_inv(-1:nx+1), nx)
    call recon(cs(-1:nx+2),   cs_if(0:nx+1,:),  dx(0:nx+1), dx_if_inv(-1:nx+1), nx)

    ! Solve the Riemann problem at all the interfaces
    call riemann(rho_if(0:nx+1,:), u_if(0:nx+1,:), &
     &     ut_if(0:nx+1,:), utt_if(0:nx+1,:), &
     &     eps_if(0:nx+1,:), pre_if(0:nx+1,:), &
     &     cs_if(0:nx+1,:), &
     &     rhoflx(0:nx), momflx(0:nx), momtflx(0:nx), &
     &     momttflx(0:nx), eneflx(0:nx), nx)

    ! Time derivatives of the conserved quantities
    do i = 1, nx
      du_con_dt(i,1) = 1.0 / dx(i) * (rhoflx (i-1) - rhoflx (i))
      du_con_dt(i,2) = 1.0 / dx(i) * (momflx (i-1) - momflx (i))
      du_con_dt(i,3) = 1.0 / dx(i) * (momtflx (i-1) - momtflx (i))
      du_con_dt(i,4) = 1.0 / dx(i) * (momttflx (i-1) - momttflx (i))
      du_con_dt(i,5) = 1.0 / dx(i) * (eneflx (i-1) - eneflx (i))
    enddo

  end subroutine sweep

end module mod_sweep
