module mod_primitives
  implicit none

contains

  ! Primitives:
  !   1: density
  !   2: x-velocity
  !   3: y-velocity
  !   4: z-velocity
  !   5: internal energy
  !   6: pressure
  !   7: sound speed

  ! Conservatives:
  !   1: density
  !   2: x-momentum
  !   3: y-momentum
  !   4: z-momentum
  !   5: total energy

  subroutine primitives(u_con, u_prim, nn)
    use mod_eos, only:eos

    real,     intent(in)  :: u_con(1:nn, 5)
    real,     intent(out) :: u_prim(1:nn, 7)
    integer,  intent(in)  :: nn

    real :: pre(1:nn), cs(1:nn)

    ! Determine the primitive variables
    u_prim(1:nn,1) = u_con(1:nn,1)
    u_prim(1:nn,2) = u_con(1:nn,2)/u_con(1:nn,1)
    u_prim(1:nn,3) = u_con(1:nn,3)/u_con(1:nn,1)
    u_prim(1:nn,4) = u_con(1:nn,4)/u_con(1:nn,1)
    u_prim(1:nn,5) = u_con(1:nn,5)/u_con(1:nn,1) &
     & - 0.5*(u_prim(1:nn,2)**2 + u_prim(1:nn,3)**2 &
     & + u_prim(1:nn,4)**2)

    ! Call the equation of state to determine pressure and sound speed
    call eos(u_prim(1:nn,1), u_prim(1:nn,5),pre(1:nn),cs(1:nn), nn)

    u_prim(1:nn,6) = pre(1:nn)
    u_prim(1:nn,7) = cs(1:nn)

  end subroutine primitives

  subroutine conserved_var(u_prim, u_con, nn)
    real,     intent(in)  :: u_prim(1:nn,7)
    real,     intent(out) :: u_con(1:nn,5)
    integer,  intent(in)  :: nn

    ! Determine the conserved variables
    u_con(1:nn,1) = u_prim(1:nn,1)
    u_con(1:nn,2) = u_prim(1:nn,2)*u_prim(1:nn,1)
    u_con(1:nn,3) = u_prim(1:nn,3)*u_prim(1:nn,1)
    u_con(1:nn,4) = u_prim(1:nn,4)*u_prim(1:nn,1)
    u_con(1:nn,5) = (u_prim(1:nn,5) &
    & + 0.5*(u_prim(1:nn,2)**2 + u_prim(1:nn,3)**2 &
    & + u_prim(1:nn,4)**2))*u_prim(1:nn,1)

  end subroutine conserved_var

end module mod_primitives
