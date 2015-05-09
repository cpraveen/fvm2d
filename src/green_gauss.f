C Inviscid flow, compute only derivatives
      subroutine green_gauss_euler(coord, mcarea, prim, qx, qy)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), prim(nvar,npmax), 
     &                 qx(nvar,npmax), qy(nvar,npmax),
     &                 mcarea(npmax)

      integer          is, iv, jt, k, n1, n2, n3
      double precision dbxx(3), dbyy(3), dxt(nvar), dyt(nvar), ais, n1b6

      n1b6 = 1.0d0/6.0d0

      do is=1,np
         do iv=1,nvar
            qx(iv,is)   = 0.0d0
            qy(iv,is)   = 0.0d0
         enddo
      enddo

c     loop on global list of triangles
      do jt=1,nt
         n1          = elem(1,jt)
         n2          = elem(2,jt)
         n3          = elem(3,jt)
         call tri_gg(coord(1,n1), coord(1,n2), coord(1,n3),
     &               prim(1,n1), prim(1,n2), prim(1,n3),
     &               dbxx, dbyy, dxt, dyt, qx(1,n1), qx(1,n2), qx(1,n3),
     &               qy(1,n1), qy(1,n2), qy(1,n3))
      enddo

      do is=1,np
         ais         = n1b6/mcarea(is)
         do k=1,nvar
            qx(k,is) = qx(k,is)*ais
            qy(k,is) = qy(k,is)*ais
         enddo
      enddo

      return
      end

C Version for viscous flow, computes viscous terms also
      subroutine green_gauss_visc(coord, elarea, mcarea, prim, nut,
     &                            res, mul, mu, qx, qy)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), prim(nvar,npmax), nut(npmax),
     &                 res(nvar,npmax), mul(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax), elarea(ntmax),
     &                 mcarea(npmax), mu(npmax)

      integer          is, iv, jt, k, n1, n2, n3
      double precision dbxx(3), dbyy(3), dxt(nvar), dyt(nvar), ais, n1b6

      n1b6 = 1.0d0/6.0d0

      do is=1,np
         do iv=1,nvar
            qx(iv,is)   = 0.0d0
            qy(iv,is)   = 0.0d0
         enddo
      enddo

c     loop on global list of triangles
      do jt=1,nt
         n1 = elem(1,jt)
         n2 = elem(2,jt)
         n3 = elem(3,jt)
         call tri_gg(coord(1,n1), coord(1,n2), coord(1,n3),
     &               prim(1,n1), prim(1,n2), prim(1,n3),
     &               dbxx, dbyy, dxt, dyt, qx(1,n1), qx(1,n2), qx(1,n3),
     &               qy(1,n1), qy(1,n2), qy(1,n3))
         call tri_visc(prim(1,n1), prim(1,n2), prim(1,n3), dbxx, dbyy, 
     &                 dxt, dyt, mu(n1), mu(n2), mu(n3), 
     &                 mul(n1), mul(n2), mul(n3), 
     &                 nut(n1), nut(n2), nut(n3), 
     &                 elarea(jt), res(1,n1), res(1,n2), res(1,n3))
      enddo

      do is=1,np
         ais         = n1b6/mcarea(is)
         do k=1,nvar
            qx(k,is) = qx(k,is)*ais
            qy(k,is) = qy(k,is)*ais
         enddo
      enddo

      return
      end

c--------------------------------------------------------------------------
c     Contribution of one triangle to gradient in green-gauss formula
c--------------------------------------------------------------------------
      subroutine tri_gg(x1, x2, x3, prim1, prim2, prim3, dbxx, dbyy, 
     &                  dxt, dyt, qx1, qx2, qx3, qy1, qy2, qy3)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), x3(2), prim1(nvar), prim2(nvar),
     &                 prim3(nvar), dbxx(3), dbyy(3), dxt(nvar), 
     &                 dyt(nvar), qx1(nvar), qx2(nvar), qx3(nvar), 
     &                 qy1(nvar), qy2(nvar), qy3(nvar)

      integer          k

      dbxx(1)   = x2(2) - x3(2)
      dbxx(2)   = x3(2) - x1(2)
      dbxx(3)   = x1(2) - x2(2)
      dbyy(1)   = x3(1) - x2(1)
      dbyy(2)   = x1(1) - x3(1)
      dbyy(3)   = x2(1) - x1(1)

      do k=1,nvar
         dxt(k) = prim1(k)*dbxx(1) + prim2(k)*dbxx(2) + prim3(k)*dbxx(3)
         dyt(k) = prim1(k)*dbyy(1) + prim2(k)*dbyy(2) + prim3(k)*dbyy(3)
         qx1(k) = qx1(k) + dxt(k)
         qx2(k) = qx2(k) + dxt(k)
         qx3(k) = qx3(k) + dxt(k)
         qy1(k) = qy1(k) + dyt(k)
         qy2(k) = qy2(k) + dyt(k)
         qy3(k) = qy3(k) + dyt(k)
      enddo

      return
      end

c--------------------------------------------------------------------------
c     Contribution of one triangle to viscous flux
c--------------------------------------------------------------------------
      subroutine tri_visc(prim1, prim2, prim3, dbxx, dbyy, dxt, dyt, 
     &                    mu1, mu2, mu3, mul1, mul2, mul3, 
     &                    nut1, nut2, nut3, area, res1, res2, res3)
      implicit none
      include 'common.h'
      double precision prim1(*), prim2(*), prim3(*), dbxx(*), dbyy(*),
     &                 dxt(*), dyt(*), mu1, mu2, mu3, mul1, mul2, mul3, 
     &                 nut1, nut2, nut3, area, res1(*), res2(*), res3(*)

      double precision en(3), ext, eyt, um, vm, mut, mult, prim5t, rhot,
     &                 chi, chi3, fv1, txx, txy, tyy, gpr, efx, efy,
     &                 areat, n2b3

      n2b3 = 2.0d0/3.0d0

C     Gradient of internal energy
      en(1) = prim1(4)/prim1(1)
      en(2) = prim2(4)/prim2(1)
      en(3) = prim3(4)/prim3(1)
      ext   = en(1)*dbxx(1) + en(2)*dbxx(2) +  en(3)*dbxx(3)
      eyt   = en(1)*dbyy(1) + en(2)*dbyy(2) +  en(3)*dbyy(3)

c     Average velocity on triangle
      um  = (prim1(2) + prim2(2) + prim3(2))/3.0d0
      vm  = (prim1(3) + prim2(3) + prim3(3))/3.0d0

c     Average Reynolds number on triangle
      mut    = (mu1  + mu2  + mu3)/3.0d0
      mult   = (mul1 + mul2 + mul3)/3.0d0
      prim5t = (nut1 + nut2 + nut3)/3.0d0
      rhot   = (prim1(1) + prim2(1) + prim3(1))/3.0d0

      chi    = prim5t*rhot/mult
      chi3   = chi**3
      fv1    = chi3/(chi3 + Cv11)

c     Viscous stresses
      txx = n2b3*mut*( 2.0d0*dxt(2) - dyt(3) )
      txy =      mut*( dyt(2) + dxt(3) )
      tyy = n2b3*mut*( 2.0d0*dyt(3) - dxt(2) )

C     Viscous Energy flux
      gpr = gamma*(mult/prandtl + 
     &             rhot*prim5t*fv1/prandtl_turb)/gamma1
      efx = txx*um + txy*vm + gpr*ext
      efy = txy*um + tyy*vm + gpr*eyt

c     Add the viscous flux
      areat      = 0.25d0/area
      res1(2) = res1(2) + (txx*dbxx(1) + txy*dbyy(1))*areat
      res1(3) = res1(3) + (txy*dbxx(1) + tyy*dbyy(1))*areat
      res1(4) = res1(4) + (efx*dbxx(1) + efy*dbyy(1))*areat

      res2(2) = res2(2) + (txx*dbxx(2) + txy*dbyy(2))*areat
      res2(3) = res2(3) + (txy*dbxx(2) + tyy*dbyy(2))*areat
      res2(4) = res2(4) + (efx*dbxx(2) + efy*dbyy(2))*areat

      res3(2) = res3(2) + (txx*dbxx(3) + txy*dbyy(3))*areat
      res3(3) = res3(3) + (txy*dbxx(3) + tyy*dbyy(3))*areat
      res3(4) = res3(4) + (efx*dbxx(3) + efy*dbyy(3))*areat

      return
      end
