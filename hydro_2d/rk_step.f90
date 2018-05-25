module mod_rk_step
  implicit none

contains

  subroutine step(time, u_con, u_prim, x, y, x_if, y_if, nx, ny)
    use mod_primitives, only:primitives
    use mod_timestep, only:timestep

    real,     intent(inout)  :: time
    real,     intent(inout)  :: u_con(1:nx,1:ny,5)
    real,     intent(inout)  :: u_prim(1:nx,1:ny,7)
    real,     intent(in)     :: x(1:nx), y(1:ny)
    real,     intent(in)     :: x_if(0:nx), y_if(0:ny)
    integer,  intent(in)     :: nx, ny

    real :: u_con_old(1:nx,1:ny,5)
    real :: du_con_dt(1:nx,1:ny,5)

    real :: dt, dtnew

    integer :: i, j

    ! Get primitives from conservatives
    !$omp parallel do
    do j = 1, ny
      call primitives(u_con(:,j,:), u_prim(:,j,:), nx)
    enddo
    !$omp end parallel do

    ! Get timestep
    dt = 1.0e99
    do j = 1, ny
      call timestep(u_prim(:,j,:), x_if, nx, dtnew)
      dt = min(dt, dtnew)
    enddo
    do i = 1, nx
      call timestep(u_prim(i,:,:), y_if, ny, dtnew)
      dt = min(dt, dtnew)
    enddo

    u_con_old(:,:,:) = u_con(:,:,:)

    ! Reset du/dt
    du_con_dt(:,:,:) = 0.0

    ! -- 1st Runge-Kutta step
    call dudt(u_prim, du_con_dt, x, y, x_if, y_if, nx, ny)
    u_con(:,:,:) = u_con(:,:,:) + dt*du_con_dt(:,:,:)

    ! Get updated primitives
    !$omp parallel do
    do j = 1, ny
      call primitives(u_con(:,j,:), u_prim(:,j,:), nx)
    enddo
    !$omp end parallel do

    ! Reset du/dt
    du_con_dt(:,:,:) = 0.0

    ! -- 2nd Runge-Kutta step
    call dudt(u_prim, du_con_dt, x, y, x_if, y_if, nx, ny)
    u_con(:,:,:) = 0.5 * (u_con_old(:,:,:) + u_con(:,:,:) + dt*du_con_dt(:,:,:))

    ! Get updated primitives
    !$omp parallel do
    do j = 1, ny
      call primitives(u_con(:,j,:), u_prim(:,j,:), nx)
    enddo
    !$omp end parallel do

    time = time + dt

  end subroutine step

  subroutine dudt(u_prim, du_con_dt, x, y, x_if, y_if, nx, ny)
    use mod_sweep, only:sweep

    real,     intent(in)     :: u_prim(1:nx,1:ny,7)
    real,     intent(inout)  :: du_con_dt(1:nx,1:ny,5)
    real,     intent(in)     :: x(1:nx), y(1:ny)
    real,     intent(in)     :: x_if(0:nx), y_if(0:ny)
    integer,  intent(in)     :: nx, ny

    integer :: i, j

    ! x-sweep
    !$omp parallel do
    do j = 1, ny
      call sweep(u_prim(:,j,:), du_con_dt(:,j,:), x, x_if, nx, 1)
    enddo
    !$omp end parallel do

    ! y-sweep
    !$omp parallel do
    do i = 1, nx
      call sweep(u_prim(i,:,:), du_con_dt(i,:,:), y, y_if, ny, 2)
    enddo
    !$omp end parallel do

  end subroutine dudt

end module mod_rk_step
