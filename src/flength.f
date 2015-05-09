C.....Calculate face length vector. Also save dual grid into a file for
C.....visualization
      subroutine flength(coord, tcoord, edge, edneigh, ds, dsb)
      implicit none
      include 'param.h'
      double precision coord(2,npmax), ds(2,nemax), dsb(2,npmax),
     &                 tcoord(2,ntmax)
      integer          edge(2,nemax), edneigh(2,nemax)

      integer          i, p1, p2, n1, n2
      double precision x1, y1, x2, y2, nx, ny

      do i=1,ne
         n1 = edge(1,i)
         n2 = edge(2,i)

         p1 = edneigh(1,i)
         p2 = edneigh(2,i)

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
            x2 = 0.5d0*(coord(1,n1) + coord(1,n2))
            y2 = 0.5d0*(coord(2,n1) + coord(2,n2))
         endif

         ds(1,i) = -( y2 - y1 )
         ds(2,i) =  ( x2 - x1 )

      enddo

c     Boundary edges
      do i=1,np
         dsb(1,i) = 0.0d0
         dsb(2,i) = 0.0d0
      enddo

      do i=1,ne
         n1 =  edge(1,i)
         n2 =  edge(2,i)
         if(edneigh(1,i)*edneigh(2,i) .eq. 0)then
            nx =  coord(2,n2) - coord(2,n1)
            ny =-(coord(1,n2) - coord(1,n1))

            dsb(1,n1) = dsb(1,n1) + 0.5d0*nx
            dsb(2,n1) = dsb(2,n1) + 0.5d0*ny
            dsb(1,n2) = dsb(1,n2) + 0.5d0*nx
            dsb(2,n2) = dsb(2,n2) + 0.5d0*ny
         endif
      enddo

      return
      end

C-----------------------------------------------------------------------------
C.....For barth cell, some edges may not make any contribution, eg. when
C.....both the triangles sharing this edge are right-angled. We remove such
C.....edges from the list
C-----------------------------------------------------------------------------
      subroutine reorder_edges(iedge, ds)
      implicit none
      include 'param.h'
      integer          iedge(nemax)
      double precision ds(2,nemax)

      integer          i, ecount
      double precision length

      print*,'Checking if any cell faces have zero length...'

      netot = ne

      if(cell_type .eq. barth)then

         ecount = 0
         do i=1,ne
            length = dsqrt( ds(1,i)**2 + ds(2,i)**2 )
            if(length .ne. 0.0d0)then
               ecount        = ecount + 1
               iedge(ecount) = i
            endif
         enddo

         print*,'Removed ',ne-ecount,' edges which have zero length'
         print*,'Final numbers of edges   = ',ecount
         ne    = ecount

      else

         do i=1,ne
            iedge(i) = i
         enddo

      endif

      return
      end
