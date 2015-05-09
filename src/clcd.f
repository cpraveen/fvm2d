C----------------------------------------------------------------------
C.....Calculate lift and drag coefficients
C----------------------------------------------------------------------
      subroutine clcd(coord, prim, mu, qx, qy)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), prim(nvar,npmax), mu(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax)

      integer          i, ie, ip1, ip2
      double precision dx, dy, sx, sy, xf, yf, p, p1, p2, txx, txy, tyy,
     &                 txx1, txy1, tyy1, txx2, txy2, tyy2

      xf = 0.0d0
      yf = 0.0d0
      do i=nswe1,nswe2
         ie   = beindx(i)
         ip1  = edge(1,ie)
         ip2  = edge(2,ie)
         dx   = coord(1,ip2) - coord(1,ip1)
         dy   = coord(2,ip2) - coord(2,ip1)
         sx   = dy
         sy   =-dx

         p1   = prim(4,ip1)
         txx1 = 2.0d0/3.0d0*mu(ip1)*(2.0d0*qx(2,ip1) - qy(3,ip1))
         txy1 =             mu(ip1)*(      qy(2,ip1) + qx(3,ip1))
         tyy1 = 2.0d0/3.0d0*mu(ip1)*(2.0d0*qy(3,ip1) - qx(2,ip1))

         p2   = prim(4,ip2)
         txx2 = 2.0d0/3.0d0*mu(ip2)*(2.0d0*qx(2,ip2) - qy(3,ip2))
         txy2 =             mu(ip2)*(      qy(2,ip2) + qx(3,ip2))
         tyy2 = 2.0d0/3.0d0*mu(ip2)*(2.0d0*qy(3,ip2) - qx(2,ip2))

         p    = 0.5d0*(p1 + p2)
         txx  = 0.5d0*(txx1 + txx2)
         txy  = 0.5d0*(txy1 + txy2)
         tyy  = 0.5d0*(tyy1 + tyy2)

         xf = xf + p*sx - txx*sx - txy*sy
         yf = yf + p*sy - txy*sx - tyy*sy
      enddo

      cl =(-dsin(aoa)*xf + dcos(aoa)*yf)/(0.5d0*r_inf*q_inf**2)
      cd =( dcos(aoa)*xf + dsin(aoa)*yf)/(0.5d0*r_inf*q_inf**2)

      return
      end
