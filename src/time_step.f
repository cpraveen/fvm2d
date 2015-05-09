C----------------------------------------------------------------------
C.....Compute local time step using cfl condition
C----------------------------------------------------------------------
      subroutine time_step(prim, mu, drmin, dt)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), mu(npmax), drmin(npmax), 
     &                 dt(npmax)

      integer          i
      double precision q, a, dtv

      if(iflow .ne. inviscid)then
c        Time-step for viscous flow
         do i=1,np
            q        = dsqrt(prim(2,i)**2 + prim(3,i)**2)
            a        = dsqrt(GAMMA*prim(4,i)/prim(1,i))
            dtv      = 2.0d0*GAMMA*mu(i)/(prim(1,i)*prandtl)
            dt(i)    = CFL*drmin(i)**2/(drmin(i)*(q + a) + dtv)
         enddo
      else
c        Time-step for inviscid flow
         do i=1,np
            q        = dsqrt(prim(2,i)**2 + prim(3,i)**2)
            a        = dsqrt(GAMMA*prim(4,i)/prim(1,i))
            dt(i)    = CFL*drmin(i)/(q + a)
         enddo
      endif

c     Find global time-step
      dtglobal = 1.0d10
      do i=1,np
         dtglobal = dmin1(dtglobal, dt(i))
      enddo

      return
      end
