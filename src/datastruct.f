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

C.....Finds elements surrounding a point
C.....Taken from Lohner
C.....esup1 stores the elements
C.....ordering is such that the elements surrounding point ipoin are stored in
C.....locations esup2(ipoin)+1 to esup2(ipoin+1)
      subroutine el_surr_point(elem, esup1, esup2)
      implicit none
      include 'param.h'
      integer esup1(mesup*npmax), esup2(npmax+1), elem(nvemax,ntmax)
      integer i, ie, inode, ipoi1, ipoin, istor

      do i=1,np+1
            esup2(i) = 0
      enddo

      do ie=1,nt
         do inode=1,3
            ipoi1        = elem(inode, ie) + 1
            esup2(ipoi1) = esup2(ipoi1) + 1
         enddo
      enddo

      do ipoin=2, np+1
         esup2(ipoin) = esup2(ipoin) + esup2(ipoin-1)
      enddo

      do ie=1, nt
         do inode=1,3
            ipoin        = elem(inode, ie)
            istor        = esup2(ipoin) + 1
            esup2(ipoin) = istor
            esup1(istor) = ie
         enddo
      enddo

      do ipoin=np+1, 2, -1
         esup2(ipoin) = esup2(ipoin-1)
      enddo

      esup2(1) = 0

      return
      end

C.....Finds points surrounding a point
C.....Taken from Lohner
C.....psup1 contains the points
C.....Neighbours of ipoin are between psup2(ipoin)+1 to psup2(ipoin+1)
      subroutine pt_surr_pt(esup1, esup2, elem, psup1, psup2)
      implicit none
      include 'param.h'
      integer esup1(mesup*npmax), esup2(npmax+1)
      integer psup1(mpsup*npmax), psup2(npmax+1), elem(nvemax,ntmax)

      integer ipoin, jpoin, inode, istor, ie, iesup, lpoin(np)

      do ipoin=1,np
         lpoin(ipoin) = 0
      enddo

      psup2(1) = 0
      istor     = 0

      do ipoin=1,np
         do iesup=esup2(ipoin)+1, esup2(ipoin+1)
            ie = esup1(iesup)
            do inode=1,3
               jpoin = elem(inode, ie)
               if(jpoin.ne.ipoin .and. lpoin(jpoin).ne.ipoin) then
                  istor = istor + 1
                  psup1(istor) = jpoin
                  lpoin(jpoin) = ipoin
               endif
            enddo
         enddo
         psup2(ipoin+1) = istor
      enddo

      return
      end

      subroutine create_edge(psup1, psup2, edge)
      implicit none
      include 'param.h'
      integer psup1(mpsup*npmax), psup2(npmax+1), edge(2,nemax)
      integer ipoin, ip, neigh

      ne = 0

      do ipoin = 1,np
         do ip=psup2(ipoin)+1, psup2(ipoin+1)
            neigh = psup1(ip)
            if(neigh.gt.ipoin) then
               ne          = ne + 1
               edge(1,ne) = ipoin
               edge(2,ne) = neigh
            endif
         enddo
      enddo

      print*, 'Number of edges           = ', ne

      return
      end

C.....For each edge find the elements to its right and left
      subroutine el_surr_edge(esup1, esup2, elem, edge, edneigh)
      implicit none
      include 'param.h'
      integer esup1(mesup*npmax), esup2(npmax+1), edneigh(2, nemax)
      integer elem(nvemax,ntmax), edge(2,nemax)
      integer i, jt, n1, n2, el, tmp

      do i=1,ne
            edneigh(1,i) = 0
            edneigh(2,i) = 0
            n1 = edge(1,i)
            n2 = edge(2,i)
            do jt=esup2(n1)+1, esup2(n1+1)
                  el = esup1(jt)
                  if( (n1.eq.elem(1,el) .and. n2.eq.elem(2,el)) .or.
     &                (n1.eq.elem(2,el) .and. n2.eq.elem(3,el)) .or.
     &                (n1.eq.elem(3,el) .and. n2.eq.elem(1,el)) )
     &                  edneigh(1,i) = el
                  if( (n2.eq.elem(1,el) .and. n1.eq.elem(2,el)) .or.
     &                (n2.eq.elem(2,el) .and. n1.eq.elem(3,el)) .or.
     &                (n2.eq.elem(3,el) .and. n1.eq.elem(1,el)) )
     &                  edneigh(2,i) = el
            enddo

            if(edneigh(1,i) .eq. 0)then
                  edneigh(1,i) = edneigh(2,i)
                  edneigh(2,i) = 0
                  tmp          = edge(1,i)
                  edge(1,i)    = edge(2,i)
                  edge(2,i)    = tmp
            endif
      enddo

      return
      end

      subroutine write_grid(coord, edge, edneigh)
      implicit none
      include 'param.h'
      double precision coord(2,npmax)
      integer          edge(2,nemax), edneigh(2,nemax)

      integer          gfile, i, n1, n2

c     Write boundary edges to bd.dat
      open(unit=10, file='BD.DAT')
      do i=1,ne
      n1 = edge(1,i)
      n2 = edge(2,i)
      if( edneigh(1,i)*edneigh(2,i) .eq. 0)then
            write(10,*)coord(1,n1), coord(2,n1)
            write(10,*)coord(1,n2), coord(2,n2)
            write(10,*)
      endif
      enddo
      close(10)

c     Write grid into file grid.dat
      gfile=15
      open(unit=gfile, file='GRID.DAT')

      do i=1,ne
            n1 = edge(1,i)
            n2 = edge(2,i)
            write(gfile,*) coord(1,n1), coord(2,n1)
            write(gfile,*) coord(1,n2), coord(2,n2)
            write(gfile,*)
      enddo

      close(gfile)

c     call system('gnuplot -noraise grid.gnu &')

      return
      end

      subroutine geom_stat(elarea, cvarea, mcarea, ds)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision elarea(*), cvarea(*), mcarea(*), ds(2,*)

      integer          i, ii
      double precision flen

c     check triangle area
      maxelarea = 0.0d0
      minelarea = 1.0d8
      do i=1,nt
         maxelarea = dmax1(maxelarea, elarea(i))
         minelarea = dmin1(minelarea, elarea(i))
      enddo

      if(minelarea .le. 0.0d0)then
         print*,'Fatal: Element area is zero/negative'
         stop
      endif

      print*,'Minimum element area     =',minelarea
      print*,'Maximum element area     =',maxelarea

c     check cv area
      maxcvarea = 0.0d0
      mincvarea = 1.0d8
      do i=1,np
         maxcvarea = dmax1(maxcvarea, cvarea(i))
         mincvarea = dmin1(mincvarea, cvarea(i))
      enddo

      if(mincvarea .le. 0.0d0)then
         print*,'Fatal: Control volume area is zero/negative'
         stop
      endif

      print*,'Minimum cv area          =',mincvarea
      print*,'Maximum cv area          =',maxcvarea

c     Find minimum/maximum face lengths
      maxflen = 0.0d0
      minflen = 1.0d8
      do ii=1,ne
         i = iedge(ii)
         flen    = dsqrt( ds(1,i)**2 + ds(2,i)**2 )
         maxflen = dmax1(maxflen, flen)   
         minflen = dmin1(minflen, flen)
      enddo
      print*,'Minimum cv length        =',minflen
      print*,'Maximum cv length        =',maxflen

      return
      end

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

      subroutine make_bdedges(ptype, edge, edneigh, beindx)
      implicit none
      include 'param.h'
      integer ptype(*), edge(2,*), edneigh(2,*), beindx(*)

      integer i, is, n1, n2, ec, nswe, nffe, swedge(nemax), 
     +        ffedge(nemax)

      nswe = 0
      nffe = 0
      do 10 i=1,netot
         is = edneigh(1,i)*edneigh(2,i)
         if(is.ne.0)goto 10
         n1 = edge(1,i)
         n2 = edge(2,i)
         if(ptype(n1).eq.solid .and. ptype(n2).eq.solid)then
            nswe = nswe + 1
            swedge(nswe) = i
         endif
         if(ptype(n1).eq.farfield .and. ptype(n2).eq.farfield)then
            nffe = nffe + 1
            ffedge(nffe) = i
         endif
10    continue

      ec   = 0

c     solid wall edges
      nswe1 = 1
      do i=1,nswe
         ec        = ec + 1
         beindx(ec) = swedge(i)
      enddo
      nswe2 = ec

c     farfield edges, characteristic bc
      nffe1 = ec + 1
      do i=1,nffe
         ec        = ec + 1
         beindx(ec) = ffedge(i)
      enddo
      nffe2 = ec

      nbe = ec
      print*,'Solid    edges =',nswe,nswe1,nswe2
      print*,'Farfield edges =',nffe,nffe1,nffe2

      return
      end
