module mod_grid
  implicit none

contains

  subroutine grid(x_if, x, y_if, y, lengthx, lengthy, nx, ny)
    real,     intent(out)  :: x_if(0:nx)
    real,     intent(out)  :: x(1:nx)
    real,     intent(out)  :: y_if(0:ny)
    real,     intent(out)  :: y(1:ny)
    real,     intent(in)  :: lengthx
    real,     intent(in)  :: lengthy
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: ny

    integer :: i, j

    ! Calculate interface coordinates
    do i = 0, nx
      x_if(i) = real(i) * lengthx / real(nx)
    enddo
    do j = 0, ny
      y_if(j) = real(j) * lengthy / real(ny)
    enddo

    ! Calculate cell center coordinates (using interface coordinates)
    do i = 1, nx
      x(i) = 0.5 * (x_if(i-1) + x_if(i))
    enddo
    do j = 1, ny
      y(j) = 0.5 * (y_if(j-1) + y_if(j))
    enddo

  end subroutine grid

end module mod_grid
