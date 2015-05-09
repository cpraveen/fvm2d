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
