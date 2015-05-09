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
