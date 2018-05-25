module mod_eos
  implicit none

  real, parameter :: gamma = 5./3.

contains

  subroutine eos(rho, eps, pre, cs, nn)
    ! Ideal gas equation of state: P = (gamma - 1) * rho * epsilon

    real, intent(in)  :: rho(1:nn)
    real, intent(in)  :: eps(1:nn)
    real, intent(out) :: pre(1:nn)
    real, intent(out) :: cs(1:nn)

    ! array size
    integer, intent(in) :: nn

    integer :: i
    
    do i = 1, nn
      pre(i) = (gamma - 1.) * rho(i) * eps(i)
      cs(i) = sqrt(gamma * (gamma - 1.) * eps(i))
    enddo
  end subroutine eos

end module mod_eos
