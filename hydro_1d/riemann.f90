module mod_riemann
  implicit none

  logical, parameter :: use_hllc = .true.

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

    real    :: lplus, lminus

    real    :: rhol, rhor, elft, ergt, plft, prgt
    real    :: ul, ur, utl, utr, uttl, uttr, epsl, epsr

    ! for HLLE
    real    :: ldiff, lprod

    ! for HLLC
    real    :: lstar, pstar
    real    :: rhostar, momstar, momtstar, momttstar, estar
    real    :: scr1, ls1, rho1

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

      elft  = rhol * (0.5 * (ul**2 + utl**2 + uttl**2) + epsl)
      ergt  = rhor * (0.5 * (ur**2 + utr**2 + uttr**2) + epsr)

      lplus  = max(ul+cs_if(i,2),ur+cs_if(i+1,1)) ! plus characteristic
      lplus  = max(lplus, 0.0)
      lminus = min(ul-cs_if(i,2),ur-cs_if(i+1,1)) ! minus characteristic
      lminus = min(lminus, 0.0)

      if (use_hllc) then
        ! HLLE
        ldiff = lplus - lminus
        lprod = lplus * lminus
        rhoflx(i) = &
          & (lplus * rhol * ul - lminus * rhor * ur + &
          & lprod * (rhor - rhol)) / ldiff
        momflx(i) = &
          & (lplus * (rhol * ul ** 2 + plft) - &
          & lminus * (rhor * ur ** 2 + prgt) + &
          & lprod * (rhor * ur - rhol * ul)) / &
          & ldiff
        momtflx(i) = &
          & (lplus * rhol * ul * utl - &
          & lminus * rhor * ur * utr + &
          & lprod * (rhor * utr - rhol * utl))/ &
          & ldiff
        momttflx(i) = &
          & (lplus * rhol * ul * uttl - &
          & lminus * rhor * ur * uttr + &
          & lprod * (rhor * uttr - rhol * uttl)) / &
          & ldiff
        eneflx(i) = &
          & (lplus * ((elft + plft) * ul) - &
          & lminus * ((ergt + prgt) * ur) + &
          & lprod * (ergt - elft)) / &
          & ldiff

      else
        ! HLLC
        lstar = (prgt - plft + rhol*ul*(lminus-ul) - rhor*ur*(lplus-ur)) &
          & / (rhol*(lminus-ul) - rhor*(lplus-ur))

        pstar = plft + rhol * (lstar - ul) * (lminus - ul)
        ! pstar = prgt + rhor * (lstar - ur) * (lplus - ur)

        if (lminus >= 0.0) then
          call flux(rhol, rhol*ul, rhol*utl, rhol*uttl, ul, utl, uttl, elft, plft, &
            & rhoflx(i), momflx(i), momtflx(i), momttflx(i), eneflx(i))
        elseif (lminus < 0.0 .and. 0.0 <= lstar) then
          ls1 = 1.0 / (lminus - lstar)
          scr1 = (lminus - ul) * ls1
          rhostar    = rhol * scr1
          momstar    = (rhol * ul * (lminus - ul) + pstar - plft) * ls1
          momtstar   = rhol * utl  * scr1
          momttstar  = rhol * uttl * scr1
          estar      = (elft * (lminus - ul) + pstar * lstar - plft * ul) * ls1
          rho1 = 1.0 / rhostar
          call flux(rhostar, momstar, momtstar, momttstar, momstar*rho1, &
            & momtstar*rho1, momttstar*rho1, estar, pstar, &
            & rhoflx(i), momflx(i), momtflx(i), momttflx(i), eneflx(i))
        elseif (lstar < 0.0 .and. 0.0 < lplus) then
          ls1 = 1.0 / (lplus - lstar)
          scr1 = (lplus - ur) * ls1
          rhostar    = rhor * scr1
          momstar    = (rhor * ur * (lplus - ur) + pstar - prgt) * ls1
          momtstar   = rhor * utr  * scr1
          momttstar  = rhor * uttr * scr1
          estar      = (ergt * (lplus - ur) + pstar * lstar - prgt * ur) * ls1
          call flux(rhostar, momstar, momtstar, momttstar, momstar*rho1, &
            & momtstar*rho1, momttstar*rho1, estar, pstar, &
            & rhoflx(i), momflx(i), momtflx(i), momttflx(i), eneflx(i))
        elseif (lplus <= 0.0) then
          call flux(rhor, rhor*ur, rhor*utr, rhor*uttr, ur, utr, uttr, ergt, prgt, &
            & rhoflx(i), momflx(i), momtflx(i), momttflx(i), eneflx(i))
        else
          stop 'Something is wrong with the HLLC solver'
        endif
      endif ! use_hllc

    enddo

  end subroutine riemann

  subroutine flux(rho, mom, momt, momtt, u, ut, utt, ene, pre, &
    & rhoflx, momflx, momtflx, momttflx, eneflx)

    real, intent(in)  :: rho, mom, momt, momtt, u, ut, utt, ene, pre
    real, intent(out) :: rhoflx, momflx, momtflx, momttflx, eneflx

    rhoflx    = rho * u
    momflx    = mom * u + pre
    momtflx   = momt * u
    momttflx  = momtt * u
    eneflx    = (ene + pre) * u

  end subroutine flux

end module mod_riemann
