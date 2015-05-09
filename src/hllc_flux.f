c HLLC flux function
c From Toro, Riemann Solvers and Numerical Methods for Fluid Dynamics,
c Chapter 10
      subroutine hllc_flux(ds, priml, primr, resl, resr)
      implicit none
      include 'common.h'
      double precision ds(2), priml(nvar), primr(nvar), resl(nvar),
     &                 resr(nvar)

      integer          i
      double precision rl, ul, vl, pl, al, El, rr, ur, vr, pr, ar, Er,
     &                 ql, qr, fl(4), fr(4), du(4),
     &                 length, ct, st, flux(nvar),
     &                 a1, a2, b1, b2, t1, fact, Plr, Q, z,
     &                 rbar, abar, Ppv, Upv, ustar, pstar, Sl, Sr, Ss,
     &                 po, prmax, prmin, Ml, Nl, Mr, Nr, gl, gr

      length= dsqrt(ds(1)**2 + ds(2)**2)
      ct    = ds(1)/length
      st    = ds(2)/length

C     Left state
      rl = priml(1)
      ul = priml(2)*ct + priml(3)*st
      vl =-priml(2)*st + priml(3)*ct
      pl = priml(4)
      al = dsqrt(GAMMA*pl/rl)
      El = pl/GAMMA1 + 0.5d0*rl*(ul**2 + vl**2)

C     Right state
      rr = primr(1)
      ur = primr(2)*ct + primr(3)*st
      vr =-primr(2)*st + primr(3)*ct
      pr = primr(4)
      ar = dsqrt(GAMMA*pr/rr)
      Er = pr/GAMMA1 + 0.5d0*rr*(ur**2 + vr**2)

C PVRS state
      rbar = 0.5d0*( rl + rr )
      abar = 0.5d0*( al + ar )
      Ppv  = 0.5d0*( pl + pr - (ur - ul)*rbar*abar )
      Upv  = 0.5d0*( ul + ur - (pr - pl)/(rbar*abar) )

      if(pl .gt. pr)then
            prmax = pl
            prmin = pr
      else
            prmax = pr
            prmin = pl
      endif

      Q = prmax/prmin

      if(Q .lt. 2.0d0 .and. prmin .lt. Ppv .and. Ppv .lt. prmax)then
C PVRS
            pstar = Ppv
            ustar = Upv
      elseif(Ppv .lt. prmin)then
C TRRS
            z   = 0.5d0*(gamma - 1.0d0)/gamma
            a1  = al + ar - 0.5d0*(gamma-1.0d0)*(ur - ul)
            a2  = al/pl**z + ar/pr**z
            pstar = (a1/a2)**(1.0d0/z)

            Plr = (pl/pr)**z
            b1 = Plr*ul/al + ur/ar + 2.0d0*(Plr-1.0d0)/(gamma-1.0d0)
            b2 = Plr/al + 1.0d0/ar
            ustar = b1/b2
      else
C TSRS
            po = dmax1(0.0d0, Ppv)

            Ml = 2.0d0/(gamma+1.0d0)/rl
            Nl = (gamma-1.0d0)*pl/(gamma+1.0d0)
            gl = dsqrt(Ml/(po + Nl))

            Mr = 2.0d0/(gamma+1.0d0)/rr
            Nr = (gamma-1.0d0)*pr/(gamma+1.0d0)
            gr = dsqrt(Mr/(po + Nr))

            a1 = gl*pl + gr*pr - (ur-ul)
            a2 = gl + gr
            pstar = a1/a2
            ustar = 0.5d0*( ul + ur + (pstar-pr)*gr - (pstar-pl)*gl )
      endif

C Wave speeds
      if(pstar .le. pl)then
          ql = 1.0d0
      else
          ql = dsqrt(1.0d0 + 0.5d0*(gamma+1.0d0)*(pstar/pl-1.0d0)/gamma)
      endif

      if(pstar .le. pr)then
          qr = 1.0d0
      else
          qr = dsqrt(1.0d0 + 0.5d0*(gamma+1.0d0)*(pstar/pr-1.0d0)/gamma)
      endif

      Sl = ul - al*ql
      Ss = ustar
      Sr = ur + ar*qr

C Computation of flux
      if (0.0d0 .le. Sl                       )then
          flux(1) = rl*ul
          flux(2) = pl + rl*ul**2
          flux(3) = rl*ul*vl
          flux(4) = (El + pl)*ul
      elseif(Sl .le. 0.0d0 .and. 0.0d0 .le. Ss)then
          fl(1) = rl*ul
          fl(2) = pl + rl*ul**2
          fl(3) = rl*ul*vl
          fl(4) = (El + pl)*ul

          fact  = rl*(Sl - ul)/(Sl - Ss)
          du(1) = fact - rl
          du(2) = fact*Ss - rl*ul
          du(3) = fact*vl - rl*vl
          t1    = El/rl + (Ss - ul)*(Ss + pl/rl/(Sl-ul))
          du(4) = fact*t1 - El

          do i=1,4
            flux(i) = fl(i) + Sl*du(i)
          enddo
      elseif(Ss .le. 0.0d0 .and. 0.0d0 .le. Sr)then
          fr(1) = rr*ur
          fr(2) = pr + rr*ur**2
          fr(3) = rr*ur*vr
          fr(4) = (Er + pr)*ur

          fact  = rr*(Sr - ur)/(Sr - Ss)
          du(1) = fact - rr
          du(2) = fact*Ss - rr*ur
          du(3) = fact*vr - rr*vr
          t1    = Er/rr + (Ss - ur)*(Ss + pr/rr/(Sr-ur))
          du(4) = fact*t1 - Er

          do i=1,4
            flux(i) = fr(i) + Sr*du(i)
          enddo
      else
          flux(1) = rr*ur
          flux(2) = pr + rr*ur**2
          flux(3) = rr*ur*vr
          flux(4) = (Er + pr)*ur
      endif

c     Total flux
      do i=1,4
         flux(i) = flux(i)*length
         resl(i) = resl(i) + flux(i)
         resr(i) = resr(i) - flux(i)
      enddo

      return
      end
