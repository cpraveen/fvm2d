C-----------------------------------------------------------------------------
C.....Calculate element and control volume areas for median cell
C-----------------------------------------------------------------------------
      subroutine tri_coord(coord, tcoord, elem)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), tcoord(2,ntmax)

      if(cell_type .eq. median)then
         call tcoord_mc(coord, tcoord, elem)
      else
         call tcoord_bc(coord, tcoord, elem)
      endif

      return
      end

C-----------------------------------------------------------------------------
C.....Calculate element and control volume areas for median cell
C-----------------------------------------------------------------------------
      subroutine tcoord_mc(coord, tcoord, elem)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), tcoord(2,ntmax)

      integer          i, n1, n2, n3

      print*,'Finding tcoord for MEDIAN cell'

      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)
         tcoord(1,i) = (coord(1,n1) + coord(1,n2) + coord(1,n3))/3.0d0
         tcoord(2,i) = (coord(2,n1) + coord(2,n2) + coord(2,n3))/3.0d0
      enddo

      return
      end

C-----------------------------------------------------------------------------
C.....Calculate element and control volume areas for Barth cell
C-----------------------------------------------------------------------------
      subroutine tcoord_bc(coord, tcoord, elem)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), tcoord(2,ntmax)

      integer          i, n1, n2, n3

      print*,'Finding tcoord for BARTH cell'

      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)

         call circumcenter(i, elem, coord, tcoord)
      enddo

      return
      end

C-----------------------------------------------------------------------------
C.....Calculate element and control volume areas
C-----------------------------------------------------------------------------
      subroutine circumcenter(el, elem, coord, tcoord)
      implicit none
      include 'param.h'
      integer          el, elem(nvemax,ntmax)
      double precision coord(2,npmax), tcoord(2,ntmax)

      double precision dx1, dy1, dx2, dy2, dx3, dy3, l1, l2, l3, 
     &                 beta1, beta2, beta3, det, b1, b2, xc, yc, 
     &                 x1, x2, x3, y1, y2, y3
      integer          n1, n2, n3

      n1    = elem(1,el)
      n2    = elem(2,el)
      n3    = elem(3,el)

      x1    = coord(1,n1)
      y1    = coord(2,n1)
      x2    = coord(1,n2)
      y2    = coord(2,n2)
      x3    = coord(1,n3)
      y3    = coord(2,n3)

      dx1   = x2 - x3
      dy1   = y2 - y3
      l1    = dx1**2 + dy1**2

      dx2   = x3 - x1
      dy2   = y3 - y1
      l2    = dx2**2 + dy2**2

      dx3   = x1 - x2
      dy3   = y1 - y2
      l3    = dx3**2 + dy3**2

      beta1 = dmax1(0.0d0, l2 + l3 - l1)
      beta2 = dmax1(0.0d0, l3 + l1 - l2)
      beta3 = dmax1(0.0d0, l1 + l2 - l3)

C This fix is supposed to remove very small cv faces.
C I am not totally happy with this one.
c     if(beta1.lt.beta2/2.0d0 .and. beta1.lt.beta3/2.0d0) beta1=0.0d0
c     if(beta2.lt.beta3/2.0d0 .and. beta2.lt.beta1/2.0d0) beta2=0.0d0
c     if(beta3.lt.beta1/2.0d0 .and. beta3.lt.beta2/2.0d0) beta3=0.0d0

C Find circumcenter
      det   = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
      b1    = 0.5d0*( (x2 - x1)*(x2 + x1) + (y2 - y1)*(y2 + y1) )
      b2    = 0.5d0*( (x3 - x1)*(x3 + x1) + (y3 - y1)*(y3 + y1) )
      xc    = ( (y3-y1)*b1 - (y2-y1)*b2)/det
      yc    = (-(x3-x1)*b1 + (x2-x1)*b2)/det

      if(beta1 .eq. 0.0d0)then
         tcoord(1,el) = 0.5d0*(x2 + x3)
         tcoord(2,el) = 0.5d0*(y2 + y3)
      elseif(beta2 .eq. 0.0d0)then
         tcoord(1,el) = 0.5d0*(x3 + x1)
         tcoord(2,el) = 0.5d0*(y3 + y1)
      elseif(beta3 .eq. 0.0d0)then
         tcoord(1,el) = 0.5d0*(x1 + x2)
         tcoord(2,el) = 0.5d0*(y1 + y2)
      else
         tcoord(1,el) = xc
         tcoord(2,el) = yc
      endif

      return
      end
