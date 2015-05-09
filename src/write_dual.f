C.....Calculate face length vector. Also save dual grid into a file for
C.....visualization
      subroutine write_dual(coord, tcoord, edge, iedge, edneigh)
      implicit none
      include 'param.h'
      double precision coord(2,npmax), tcoord(2,ntmax)
      integer          edge(2,nemax), iedge(nemax), edneigh(2,nemax)

      integer          i, ii, p1, p2, n1, n2, idual
      double precision x1, y1, x2, y2, xm, ym, nx, ny, flen

      idual = 20
      open(unit=idual, file='DUAL.DAT')

      do i=1,netot
         n1 = edge(1,i)
         n2 = edge(2,i)

         p1 = edneigh(1,i)
         p2 = edneigh(2,i)

         x1 = 0.0d0
         y1 = 0.0d0
         x2 = 0.0d0
         y2 = 0.0d0

c        Mid-point of the edge
         xm = 0.5d0*(coord(1,n1) + coord(1,n2))
         ym = 0.5d0*(coord(2,n1) + coord(2,n2))

         if(p1 .ne. 0)then
            x1 = tcoord(1,p1)
            y1 = tcoord(2,p1)
         else
            print*,'flength: Fatal error at edge',i
            print*,'         p1 is zero'
            stop
         endif

         if(p2 .ne. 0)then
            x2 = tcoord(1,p2)
            y2 = tcoord(2,p2)
         else
c           Edge number i is a boundary edge
            x2 = xm
            y2 = ym
         endif

         write(idual,*)x1, y1
         write(idual,*)xm, ym
         write(idual,*)
         write(idual,*)x2, y2
         write(idual,*)xm, ym
         write(idual,*)
      enddo
      close(idual)

      return
      end
