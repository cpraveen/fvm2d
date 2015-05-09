C.....Check whether ordering of triangle is counter-clockwise
C.....Otherwise correct it
      subroutine tri_orient(elem, coord)
      implicit none
      include 'param.h'
      integer           elem(nvemax,ntmax)
      double precision  coord(2,npmax)

      integer           cw, ccw, tmp, i, p1, p2, p3
      double precision  dx1, dy1, dx2, dy2, cross

      cw = 0
      ccw= 0

      do i=1,nt
            p1    = elem(1,i)
            p2    = elem(2,i)
            p3    = elem(3,i)

            dx1   = coord(1,p2) - coord(1,p1)
            dy1   = coord(2,p2) - coord(2,p1)

            dx2   = coord(1,p3) - coord(1,p2)
            dy2   = coord(2,p3) - coord(2,p2)

            cross = dx1*dy2 - dx2*dy1

            if(cross .eq. 0.0d0)then
                  print*,'Fatal: triangle',i,' is degenerate'
                  stop
            endif

            if(cross .lt. 0.0d0)then
                  cw        = cw + 1
                  tmp       = elem(2,i)
                  elem(2,i) = elem(3,i)
                  elem(3,i) = tmp
            else
                  ccw       = ccw + 1
            endif
      enddo

      print*,'No of cw  triangles       = ',cw
      print*,'No of ccw triangles       = ',ccw

      return
      end
