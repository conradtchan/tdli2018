module mod_output
  implicit none

contains

  subroutine write_output(nstep, x, y, u_prim, nx, ny)
    integer,  intent(in) :: nstep
    real,     intent(in) :: x(1:nx), y(1:ny)
    real,     intent(in) :: u_prim(1:nx,1:ny,7)
    integer,  intent(in) :: nx, ny

    integer :: i, j
    character(len=20) :: filename

    write(filename, "('output/out',I6.6,'.dat')") nstep

    if (nstep == 0) then
      open(60, file="output/x.dat", status='new', form='formatted')
      write(60,*) x(:)
      close(60)

      open(70, file="output/y.dat", status='new', form='formatted')
      write(70,*) y(:)
      close(70)
    endif

    open(80, file=filename, status='new', form='formatted')

    do j = 1, ny
      do i = 1, nx
         write(80,*) u_prim(i,j,:)
      end do
    enddo

    close(80)

  end subroutine write_output

end module mod_output
