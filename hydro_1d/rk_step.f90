module mod_rk_step
  implicit none

contains

  subroutine step(time, u_con, u_prim, x, x_if, nx)
    use mod_primitives, only:primitives
    use mod_timestep, only:timestep
    use mod_sweep, only:sweep

    real,     intent(inout)  :: time
    real,     intent(inout)  :: u_con(1:nx,5)
    real,     intent(inout)  :: u_prim(1:nx,7)
    real,     intent(in)     :: x(1:nx)
    real,     intent(in)     :: x_if(0:nx)
    integer,  intent(in)  :: nx

    real :: u_con_old(1:nx,5)
    real :: du_con_dt(1:nx,5)

    real :: dt

    ! Get primitives from conservatives
    call primitives(u_con, u_prim, nx)

    ! Get timestep
    call timestep(u_prim, x_if, nx, dt)

    u_con_old(:,:) = u_con(:,:)

    ! 1st Runge-Kutta step
    call sweep(u_prim, du_con_dt, x, x_if, nx)
    u_con(:,:) = u_con(:,:) + dt*du_con_dt(:,:)

    ! Get updated primitives
    call primitives(u_con, u_prim, nx)

    ! 2nd Runge-Kutta step
    call sweep(u_prim, du_con_dt, x, x_if, nx)
    u_con(:,:) = 0.5 * (u_con_old(:,:) + u_con(:,:) + dt*du_con_dt(:,:))

    time = time + dt

  end subroutine step

end module mod_rk_step
