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
