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
