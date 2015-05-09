C.....Calculate element and control volume areas for median cell
      subroutine areas(coord, tcoord, elem, elarea, cvarea, mcarea)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), elarea(ntmax), cvarea(npmax),
     &                 tcoord(2,ntmax), mcarea(npmax)

      call tri_area(coord, elem, elarea)
      call mc_area(elem, elarea, mcarea)

      if(cell_type .eq. median)then
         call cvareas_mc(cvarea, mcarea)
      else
         call cvareas_bc(coord, tcoord, elem, cvarea)
      endif

      return
      end

c---------------------------------------------------------------------------
c     Computes areas of triangles
c---------------------------------------------------------------------------
      subroutine tri_area(coord, elem, elarea)
      implicit none
      include 'param.h'
      integer          elem(3,*)
      double precision coord(2,*), elarea(*)

      integer          i, n1, n2, n3
      double precision dx1, dy1, dx2, dy2

      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)
         dx1= coord(1,n2) - coord(1,n1)
         dy1= coord(2,n2) - coord(2,n1)
         dx2= coord(1,n3) - coord(1,n1)
         dy2= coord(2,n3) - coord(2,n1)
         elarea(i) = 0.5d0*( dx1*dy2 - dx2*dy1 )
      enddo

      return
      end

c---------------------------------------------------------------------------
c.....Area of median cell
c---------------------------------------------------------------------------
      subroutine mc_area(elem, elarea, mcarea)
      implicit none
      include 'param.h'
      integer          elem(3,*)
      double precision elarea(*), mcarea(*)

      integer          i, n1, n2, n3
      double precision area3
      
      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)
         area3        = elarea(i)/3.0d0
         mcarea(n1) = mcarea(n1) + area3
         mcarea(n2) = mcarea(n2) + area3
         mcarea(n3) = mcarea(n3) + area3
      enddo

      return
      end

c---------------------------------------------------------------------------
C.....Calculate element and control volume areas for median cell
c---------------------------------------------------------------------------
      subroutine cvareas_mc(cvarea, mcarea)
      implicit none
      include 'param.h'
      double precision cvarea(*), mcarea(*)

      integer          i

      print*,'Finding control volume areas for MEDIAN cell'

      do i=1,np
         cvarea(i) = mcarea(i)
      enddo

      return
      end

c---------------------------------------------------------------------------
C.....Calculate element and control volume areas for Barth cell
c---------------------------------------------------------------------------
      subroutine cvareas_bc(coord, tcoord, elem, cvarea)
      implicit none
      include 'param.h'
      integer          elem(3,*)
      double precision coord(2,*), cvarea(*), tcoord(2,*)

      double precision x(5), y(5), area
      integer          i, j, k, n1, n2, n3

      print*,'Finding control volume areas for BARTH cell'

      do i=1,np
         cvarea(i)   = 0.0d0
      enddo

      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)

C        Control volume area
         do j=1,3
            n1 = elem(j,i)

            if(j .eq. 1)then
               n2 = elem(3,i)
            else
               n2 = elem(j-1,i)
            endif

            if(j .eq. 3)then
               n3 = elem(1,i)
            else
               n3 = elem(j+1,i)
            endif

            x(1) = coord(1,n1)
            y(1) = coord(2,n1)
            x(2) = 0.5d0*(coord(1,n1) + coord(1,n3))
            y(2) = 0.5d0*(coord(2,n1) + coord(2,n3))
            x(3) = tcoord(1,i)
            y(3) = tcoord(2,i)
            x(4) = 0.5d0*(coord(1,n1) + coord(1,n2))
            y(4) = 0.5d0*(coord(2,n1) + coord(2,n2))
            x(5) = x(1)
            y(5) = y(1)

            area = 0.0d0
            do k=1,4
               area = area + x(k)*y(k+1) - x(k+1)*y(k)
            enddo
            if(area .lt. 0.0d0)then
               print*,'Area is non-positive'
               print*,' Area         =', area
               print*,' Triangle     =',i
               stop 'areas'
            endif
            area       = 0.5d0*area
            cvarea(n1) = cvarea(n1) + area
         enddo
      enddo

      return
      end
