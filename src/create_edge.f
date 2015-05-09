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
