C------------------------------------------------------------------------------
C.....Update the solution
C.....Current implementation is single step forward Euler
C------------------------------------------------------------------------------
      subroutine update(irk, res, prim_old, prim, cvarea, dt)
      implicit none
      include 'param.h'
      include 'gdata.h'
      integer          irk
      double precision res(nvar,npmax), prim_old(nvar,npmax), 
     &                 prim(nvar,npmax), cvarea(npmax), dt(npmax)

      integer          i, j
      double precision con0_old(nvar), con1_old(nvar),
     &                 con_new(nvar,npmax)

      do i=1,np
         call prim2con(prim_old(1,i), con0_old)
         call prim2con(prim(1,i),     con1_old)
         do j=1,nvar
            con_new(j,i) = airk(irk)*con0_old(j) + 
     &                     birk(irk)*(con1_old(j) - 
     &                     (dt(i)/cvarea(i))*res(j,i))
         enddo
      enddo
                  
      if(iflow .ne. inviscid)then
         do i=1,nsp
            j = spts(i)
            con_new(2,j) = 0.0d0
            con_new(3,j) = 0.0d0
         enddo
      endif

      call con_to_prim(con_new, prim)

      return
      end
