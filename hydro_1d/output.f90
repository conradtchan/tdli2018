module mod_output
  implicit none

contains

  subroutine write_output(nstep, x, u_prim, nx)
    integer,  intent(in) :: nstep
    real,     intent(in) :: x(1:nx)
    real,     intent(in) :: u_prim(1:nx,7)
    integer,  intent(in) :: nx

    integer :: i
    character(len=20) :: filename

    write(filename, "('output/out',I6.6,'.dat')") nstep

    open(80, file=filename,status='new', form = 'formatted')

    do i = 1, nx
       write(80,*) x(i), u_prim(i,1), u_prim(i,2), u_prim(i,3), u_prim(i,4), u_prim(i,5), u_prim(i,6), u_prim(i,7)
    end do

    close(80)

  end subroutine write_output

end module mod_output
