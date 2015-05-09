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
