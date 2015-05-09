C----------------------------------------------------------------------
C.....Calculate length used for time step
C.....For each point find the minimum altitude of all triangles
C.....surrounding that point
C----------------------------------------------------------------------
      subroutine dtlength(coord, elarea, elem, drmin)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), elarea(ntmax), drmin(npmax)

      integer          ip, it, n1, n2, n3
      double precision dx1, dx2, dx3, dy1, dy2, dy3, dr1, dr2, dr3, 
     &                 h1, h2, h3


      do ip=1,np
         drmin(ip) = 1.0d15
      enddo

      do it=1,nt
         n1  = elem(1,it)
         n2  = elem(2,it)
         n3  = elem(3,it)

         dx1 = coord(1,n2) - coord(1,n3)
         dy1 = coord(2,n2) - coord(2,n3)
         dr1 = dsqrt(dx1**2 + dy1**2)
         h1  = 2.0d0*elarea(it)/dr1

         dx2 = coord(1,n3) - coord(1,n1)
         dy2 = coord(2,n3) - coord(2,n1)
         dr2 = dsqrt(dx2**2 + dy2**2)
         h2  = 2.0d0*elarea(it)/dr2

         dx3 = coord(1,n1) - coord(1,n2)
         dy3 = coord(2,n1) - coord(2,n2)
         dr3 = dsqrt(dx3**2 + dy3**2)
         h3  = 2.0d0*elarea(it)/dr3

         drmin(n1) = dmin1( h1, drmin(n1) )
         drmin(n2) = dmin1( h2, drmin(n2) )
         drmin(n3) = dmin1( h3, drmin(n3) )
      enddo

      return
      end
