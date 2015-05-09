C----------------------------------------------------------------------
C.....Save the current solution into another array
C----------------------------------------------------------------------
      subroutine save_old(prim, prim_old)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), prim_old(nvar,npmax)

      integer          i, j

      do i=1,np
         do j=1,nvar
            prim_old(j,i) = prim(j,i)
         enddo
      enddo

      return
      end
