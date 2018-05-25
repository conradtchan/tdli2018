module mod_reconstruction
  implicit none

contains

  subroutine coeff(x, x_if, dx, dx_if_inv, nx)
    real,     intent(in)  :: x(-1:nx+2), x_if(-2:nx+2)
    real,     intent(out) :: dx(0:nx+1)
    real,     intent(out) :: dx_if_inv(-1:nx+1)
    integer,  intent(in)  :: nx

    integer :: i

    do i = 0, nx + 1
      dx(i) = x_if(i) - x_if(i-1)
    enddo

    do i = 01, nx + 1
      dx_if_inv = 1.0 / (x(i+1) - x(i))
    enddo

  end subroutine coeff

  subroutine recon(q, q_if, dx, dx_if_inv, nx)
    real,     intent(in)  :: q(-1:nx+2)
    real,     intent(in)  :: dx(0:nx+1)
    real,     intent(in)  :: dx_if_inv(-1:nx+1)
    real,     intent(out) :: q_if(0:nx+1,1:2)
    integer,  intent(in)  :: nx

    real    :: sl, sr ! left and right slopes
    real    :: slim   ! limited slope

    integer :: i

    do i = 0, nx + 1
      sl = (q(i)   - q(i-1)) * dx_if_inv(i-1)
      sr = (q(i+1) - q(i)  ) * dx_if_inv(i)

      slim = mc_lim(sl, sr)

      ! 2nd order reconstruction
      q_if(i,1) = q(i) - 0.5 * dx(i) * slim  ! left side of cell, right side of interface
      q_if(i,2) = q(i) + 0.5 * dx(i) * slim  ! right side of cell, left side of interface
    enddo

  end subroutine recon

  real function mc_lim(a, b)
    real, intent(in)  :: a, b

    real :: scr

    if (a*b > 0.0) then
      scr = 2.0 * min(abs(a), abs(b))
      mc_lim = sign(min(scr, 0.5*abs(a+b)), a)
    else
      mc_lim = 0.0
    endif

  end function mc_lim

end module mod_reconstruction
