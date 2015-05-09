C.....Calculate L2 residue based on density
C.....Normalize with residue of first iteration, res1
      subroutine conres(prim, prim_old, dt)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), prim_old(nvar,npmax), dt(npmax)

      integer          i
      double precision gra, dif

      fres  = 0.0d0
      fresi = 0.0d0
      iresi = 0
      do i=1,np
         gra = prim(1,i) - prim_old(1,i)
         fres = fres + gra**2
         dif = dabs(gra)
         if(dif .gt. fresi)then
            fresi = dif
            iresi= i
         endif
      enddo
      fres = dsqrt(fres/np)

      if(iter .eq. iterlast)then
         fres1 = fres
         print*,'Residue in first iteration = ',fres1
         if(fres1 .eq. 0.0d0) stop
      endif

      fres = fres/fres1

      return
      end
