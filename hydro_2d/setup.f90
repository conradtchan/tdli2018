module mod_setup
  implicit none

contains

  subroutine setup(x, y, u_prim, xlen, ylen, nx, ny)
    use mod_eos, only:eos,gamma
    real,     intent(in)  :: x(1:nx), y(1:ny)
    real,     intent(out) :: u_prim(1:nx,1:ny,7)
    real,     intent(in)  :: xlen, ylen
    integer,  intent(in)  :: nx, ny

    real :: rho(1:nx,1:ny)
    real :: v1(1:nx,1:ny)
    real :: v2(1:nx,1:ny)
    real :: v3(1:nx,1:ny)
    real :: eps(1:nx,1:ny)
    real :: pre(1:nx,1:ny)

    real    :: scr, dy
    integer :: i, j

    real, parameter :: pi = acos(-1.0)

    dy = ylen / real(ny)

    ! KH initial conditions
    rho(:,:) = 1.0
    pre(:,:) = 1.0
    v1(:,:) = 0.0
    v2(:,:) = 0.0
    v3(:,:) = 0.0

    do j = 1, ny
      do i = 1, nx
        if (y(j) < 0.25 * ylen + 0.05 * cos(2.0 * pi * 2 * x(i) / xlen)) then
          v1(i,j) = -0.5
          v3(i,j) = -0.01
        elseif (y(j) > 0.75 * ylen - 0.05 * cos(2.0 * pi * 2 * x(i) / xlen)) then
          v1(i,j) = -0.5
          v3(i,j) = -0.01
        else
          v1(i,j) = 0.5
          v3(i,j) = 0.01
        endif

        if (x(i) > 0.25 .and. x(i) < 0.75) then
          v3(i, j) = v3(i, j) * 2.0
        endif
        
      enddo
    enddo

    ! Calculate eps
    eps(1:nx,1:ny) = pre(1:nx,1:ny) / rho(1:nx,1:ny) / (gamma - 1.0)

    u_prim(1:nx,1:ny,1) = rho (1:nx,1:ny)
    u_prim(1:nx,1:ny,2) = v1  (1:nx,1:ny)
    u_prim(1:nx,1:ny,3) = v2  (1:nx,1:ny)
    u_prim(1:nx,1:ny,4) = v3  (1:nx,1:ny)
    u_prim(1:nx,1:ny,5) = eps (1:nx,1:ny)

    ! Use EOS to get pressure and sound speed
    do j = 1, ny
      call eos(u_prim(1:nx,j,1), u_prim(1:nx,j,5), u_prim(1:nx,j,6), u_prim(1:nx,j,7), nx)
    enddo
  end subroutine setup

end module mod_setup
