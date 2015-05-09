C.....Kinetic split fluxes
      subroutine kfvs_flux(ds, priml, primr, resl, resr)
      implicit none
      include 'common.h'
      double precision ds(2), priml(nvar), primr(nvar), resl(nvar),
     &                 resr(nvar)

      integer          i
      double precision rl, ul, vl, pl, el, rr, ur, vr, pr, er,
     &                 ql2, qr2, unl, unr, vnl, vnr,
     &                 length, ct, st, Al, Bl, Ar, Br,
     &                 sl, betal, sr, betar,
     &                 Fp(4), Fm(4), Ff(4), flux(4)
      intrinsic        derf

      length= dsqrt(ds(1)**2 + ds(2)**2)
      ct    = ds(1)/length
      st    = ds(2)/length

C     Left state
      rl = priml(1)
      ul = priml(2)
      vl = priml(3)
      pl = priml(4)
      ql2= ul**2 + vl**2
      el = pl/GAMMA1 + 0.5d0*rl*ql2

C     Right state
      rr = primr(1)
      ur = primr(2)
      vr = primr(3)
      pr = primr(4)
      qr2= ur**2 + vr**2
      er = pr/GAMMA1 + 0.5d0*rr*qr2

C     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

      vnl =-ul*st + vl*ct
      vnr =-ur*st + vr*ct

c     Positive flux
      betal = 0.5d0*rl/pl
      sl    = unl*dsqrt(betal)
      Al    = 0.5d0*(1.0d0 + DERF(sl))
      Bl    = 0.5d0*dexp(-sl**2)/dsqrt(M_PI*betal)

      Fp(1) = rl*(unl*Al + Bl)
      Fp(2) = (pl + rl*unl**2)*Al + rl*unl*Bl
      Fp(3) = rl*(unl*Al + Bl)*vnl
      Fp(4) = (el + pl)*unl*Al + (el + 0.5d0*pl)*Bl

c     Negative flux
      betar = 0.5d0*rr/pr
      sr    = unr*dsqrt(betar)
      Ar    = 0.5d0*(1.0d0 - DERF(sr))
      Br    = 0.5d0*dexp(-sr**2)/dsqrt(M_PI*betar)

      Fm(1) = rr*(unr*Ar - Br)
      Fm(2) = (pr + rr*unr**2)*Ar - rr*unr*Br
      Fm(3) = rr*(unr*Ar - Br)*vnr
      Fm(4) = (er + pr)*unr*Ar - (er + 0.5d0*pr)*Br

c     Total flux
      do i=1,4
            Ff(i) = Fp(i) + Fm(i)
      enddo

      flux(1) = Ff(1)*length
      flux(2) = (ct*Ff(2) - st*Ff(3))*length
      flux(3) = (st*Ff(2) + ct*Ff(3))*length
      flux(4) = Ff(4)*length

      do i=1,4
         resl(i) = resl(i) + flux(i)
         resr(i) = resr(i) - flux(i)
      enddo

      return
      end
