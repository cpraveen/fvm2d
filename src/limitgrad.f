      subroutine limitgrad(edge, coord, prim, qx, qy)
      implicit none
      include 'param.h'

      integer          edge(2,nemax)
      double precision qx(nvar,npmax), qy(nvar,npmax), coord(2,npmax),
     &                 prim(nvar,npmax)

      integer i, j, n1, n2
      double precision qmin(nvar,npmax), qmax(nvar,npmax),
     &                 phi(nvar,npmax), dx, dy, rat1, rat2, 
     &                 dqmax1, dqmax2, dqmin1, dqmin2, dq1, dq2, ppp

      do i=1,np
            do j=1,nvar
                  qmin(j,i) =  1.0d10
                  qmax(j,i) = -1.0d10
                  phi(j,i)  =  1.0d0
            enddo
      enddo

      do i=1,netot
            n1 = edge(1,i)
            n2 = edge(2,i)
            do j=1,nvar
                  qmin(j,n1) = dmin1(qmin(j,n1), prim(j,n1))
                  qmax(j,n1) = dmax1(qmax(j,n1), prim(j,n1))
                  qmin(j,n2) = dmin1(qmin(j,n2), prim(j,n2))
                  qmax(j,n2) = dmax1(qmax(j,n2), prim(j,n2))
             enddo
      enddo

      do i=1,netot
            n1 = edge(1,i)
            n2 = edge(2,i)
            dx = coord(1,n2) - coord(1,n1)
            dy = coord(2,n2) - coord(2,n1)
            do j=1,nvar
                  dq1    = qx(j,n1)*dx + qy(j,n1)*dy
                  dqmax1 = qmax(j,n1) - prim(j,n1)
                  dqmin1 = qmin(j,n1) - prim(j,n1)
                  if(dq1 .ne. 0.0d0)then
                        rat1   = dqmax1/dq1
                        rat2   = dqmin1/dq1
                        if(dq1 .gt. 0.0d0)then
                              ppp = rat1
                        else
                              ppp = rat2
                        endif
                        phi(j,n1) = dmin1(phi(j,n1), ppp)
                  endif

                  dq2    = -( qx(j,n2)*dx + qy(j,n2)*dy )
                  dqmax2 = qmax(j,n2) - prim(j,n2)
                  dqmin2 = qmin(j,n2) - prim(j,n2)
                  if(dq2 .ne. 0.0d0)then
                        rat1   = dqmax2/dq2
                        rat2   = dqmin2/dq2
                        if(dq2 .gt. 0.0d0)then
                              ppp = rat1
                        else
                              ppp = rat2
                        endif
                        phi(j,n2) = dmin1(phi(j,n2), ppp)
                  endif
            enddo

      enddo

      do i=1,np
            do j=1,nvar
                  qx(j,i) = phi(j,i)*qx(j,i) 
                  qy(j,i) = phi(j,i)*qy(j,i)
            enddo
      enddo

      return
      end
