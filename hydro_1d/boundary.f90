module mod_boundary
  implicit none

contains

  subroutine boundary(q, nx, bndcon, parity, valuel, valuer)
    ! q is the array of solution values
    ! indices -1 and 0 are the left ghost cells
    ! indices nx+1 and nx+2 are the right ghost cells
    real,     intent(inout) :: q(-1:nx+2)
    ! size of array
    integer, intent(in)     :: nx
    ! select type of boundary condition
    integer,  intent(in)    :: bndcon
    ! for reflecting boundaries, select parity
    integer, optional, intent(in) :: parity
    ! for specific valued boundaries, select values
    real,    optional, intent(in) :: valuel, valuer

    select case(bndcon)
    case(1)  ! periodic boundary conditions
      q(nx+1) = q(1)
      q(nx+2) = q(2)
      q(0)    = q(nx)
      q(-1)   = q(nx-1)
    case(2)  ! reflecting boundary conditions
      if (parity /= 1 .and. parity /= -1) stop 'boundary(): Parity incorrect!'
      q(nx+1) = parity * q(nx)
      q(nx+2) = parity * q(nx-1)
      q(0)    = parity * q(1)
      q(-1)   = parity * q(2)
    case(3)  ! specific values
      q(nx+1) = valuer
      q(nx+2) = valuer
      q(0)    = valuel
      q(-1)   = valuel
    case default
      stop 'boundary(): Boundary condition ill defined!'
    end select

  end subroutine boundary

end module mod_boundary
