module mod_riemann
  implicit none

contains

  subroutine riemann(rho_if, u_if, ut_if, utt_if, eps_if, p_if, &
    & cs_if, rhoflx, momflx, momtflx, momttflx, eneflx, nx)

    real,     intent(in)  :: rho_if(0:nx+1,2)
    real,     intent(in)  :: u_if(0:nx+1,2)
    real,     intent(in)  :: ut_if(0:nx+1,2)
    real,     intent(in)  :: utt_if(0:nx+1,2)
    real,     intent(in)  :: eps_if(0:nx+1,2)
    real,     intent(in)  :: p_if(0:nx+1,2)
    real,     intent(in)  :: cs_if(0:nx+1,2)

    integer,  intent(in)  :: nx

    real,     intent(out) :: rhoflx(0:nx)
    real,     intent(out) :: momflx(0:nx)
    real,     intent(out) :: momtflx(0:nx)
    real,     intent(out) :: momttflx(0:nx)
    real,     intent(out) :: eneflx(0:nx)

    integer :: i

    real    :: lplus, lminus, ldiff, lprod
    real    :: rhol, rhor, elft, ergt, plft, prgt
    real    :: ul, ur, utl, utr, uttl, uttr, epsl, epsr
    real    :: lstar, e_hll, s1_hll

    do i = 0, nx
      rhol = rho_if (i,2)
      ul   = u_if   (i,2)
      utl  = ut_if  (i,2)
      uttl = utt_if (i,2)
      epsl = eps_if (i,2)
      plft = p_if   (i,2)

      rhor = rho_if (i+1,1)
      ur   = u_if   (i+1,1)
      utr  = ut_if  (i+1,1)
      uttr = utt_if (i+1,1)
      epsr = eps_if (i+1,1)
      prgt = p_if   (i+1,1)

      elft  = rhol * (0.5d0 * (ul**2 + utl**2 + uttl**2) + epsl)
      ergt  = rhor * (0.5d0 * (ur**2 + utr**2 + uttr**2) + epsr)

      lplus  = max(ul+cs_if(i,2),ur+cs_if(i+1,1)) ! plus characteristic
      lplus  = max(lplus, 0.0d0)
      lminus = min(ul-cs_if(i,2),ur-cs_if(i+1,1)) ! minus characteristic
      lminus = min(lminus, 0.0d0)
      ldiff = lplus - lminus
      lprod = lplus * lminus

      rhoflx(i) = &
        & (lplus * rhol * ul - lminus * rhor * ur + &
        & lprod * (rhor - rhol)) / ldiff
      momflx (i) = &
        & (lplus * (rhol * ul ** 2 + plft) - &
        & lminus * (rhor * ur ** 2 + prgt) + &
        & lprod * (rhor * ur - rhol * ul)) / &
        & ldiff
      momtflx (i) = &
        & (lplus * rhol * ul * utl - &
        & lminus * rhor * ur * utr + &
        & lprod * (rhor * utr - rhol * utl))/ &
        & ldiff
      momttflx (i) = &
        & (lplus * rhol * ul * uttl - &
        & lminus * rhor * ur * uttr + &
        & lprod * (rhor * uttr - rhol * uttl)) / &
        & ldiff
      eneflx (i) = &
        & (lplus * ((elft + plft) * ul) - &
        & lminus * ((ergt + prgt) * ur) + &
        & lprod * (ergt - elft)) / &
        & ldiff

    enddo

  end subroutine riemann

end module mod_riemann
