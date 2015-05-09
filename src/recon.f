c------------------------------------------------------------------------------
c Limited reconstruction using van-albada limiter
c------------------------------------------------------------------------------
      subroutine recon(x1, x2, prim1, prim2, qx1, qy1, qx2, qy2, 
     &                 priml, primr)
      implicit none
      include 'param.h'
      double precision x1(2), x2(2), prim1(nvar), prim2(nvar),
     &                 qx1(nvar), qx2(nvar), qy1(nvar), qy2(nvar), 
     &                 priml(nvar), primr(nvar)

      integer          i
      double precision dx, dy, dqp, dql, si, dqr, sj, limit_albada

      dx = x2(1) - x1(1)
      dy = x2(2) - x1(2)

      if(ilimit .eq. yes)then

         do i=1,nvar
            dqp = prim2(i) - prim1(i)

            dql = ALBADA21*(dx*qx1(i) + dy*qy1(i)) + ALBADA22*dqp
            si  = limit_albada(dql, dqp)
            priml(i) = prim1(i) + 0.5d0*si

            dqr = ALBADA21*(dx*qx2(i) + dy*qy2(i)) + ALBADA22*dqp
            sj  = limit_albada( dqr, dqp )
            primr(i) = prim2(i) - 0.5d0*sj
         enddo

      else

         do i=1,nvar
            dqp = prim2(i) - prim1(i)

            dql = dx*qx1(i) + dy*qy1(i)
            si  = ALBADA11*dql + ALBADA12*dqp
            priml(i) = prim1(i) + 0.5d0*si

            dqr = dx*qx2(i) + dy*qy2(i)
            sj  = ALBADA11*dqr + ALBADA12*dqp
            primr(i) = prim2(i) - 0.5d0*sj
         enddo

      endif


      return
      end

c------------------------------------------------------------------------------
c van-albada limiter
c------------------------------------------------------------------------------
      double precision function limit_albada(ul, ur)
      implicit none
      include 'param.h'
      double precision ul, ur

      double precision top, bot

      if( ul*ur .le. 0.0d0)then
         limit_albada = 0.0d0
      else
         top = (ul**2 + EPSILON)*ur + (ur**2 + EPSILON)*ul
         bot = ul**2 + ur**2 + 2.0d0*EPSILON
         limit_albada = top/bot
      endif

      return
      end
