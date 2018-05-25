module mod_timestep
  implicit none

  real, parameter :: cfl = 0.7

contains

  subroutine timestep(u, x_if, nn, dt)
    real,     intent(in)  :: u(1:nn,7)
    real,     intent(in)  :: x_if(0:nn)
    real,     intent(out) :: dt
    integer,  intent(in)  :: nn

    real    :: lambda_max
    integer :: i

    dt = 0.0

    ! determine the dt requirement based on all cells
    do i = 1, nn
      lambda_max = abs(sqrt(u(i,2)**2 + u(i,3)**2 + u(i,4)**2)) + u(i,7)
      dt = max(dt, lambda_max/(x_if(i)-x_if(i-1)))
    enddo

    ! up until this point, dt is inverted
    dt = cfl / dt

  end subroutine timestep

end module mod_timestep
