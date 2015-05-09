      subroutine shockind(edge, coord, ds, drmin, cvarea, prim, qx, qy)
      implicit none
      include 'param.h'
      integer          edge(2,nemax)
      double precision coord(2,npmax), prim(nvar,npmax), qx(nvar,npmax),
     &                 qy(nvar,npmax), ds(2,nemax), drmin(npmax),
     &                 cvarea(npmax)

      integer          i, j, n1, n2
      double precision dx, dy, dvc, dvt, dvl, vl, dvr, vr, li, ln,
     &                 ind(npmax)

      do i=1,np
            ind(i) = 0.0d0
      enddo

      do i=1,netot
            n1 = edge(1,i)
            n2 = edge(2,i)

            dx = coord(1,n2) - coord(1,n1)
            dy = coord(2,n2) - coord(2,n1)

            dvc = prim(1,n2) - prim(1,n1)

            dvt = dx*qx(1,n1) + dy*qy(1,n1)
c           dvl = dvt
            dvl = ALBADA11*dvt + ALBADA12*dvc
            vl = prim(1,n1) + dvl

            dvt = dx*qx(1,n2) + dy*qy(1,n2)
c           dvr = dvt
            dvr = ALBADA11*dvt + ALBADA12*dvc
            vr = prim(1,n2) - dvr

            ln = dsqrt(ds(1,i)**2 + ds(2,i)**2)
            li = (vr - vl)**2*ln
            ind(n1) = ind(n1) + li
            ind(n2) = ind(n2) + li
      enddo

c     open(unit=20, file='ind.dat')
      do i=1,np
            ind(i) = ind(i)/drmin(i)/cvarea(i)**(0.75d0)
c           write(20,*) ind(i)
            if( ind(i) .gt. 1.0d0 )then
                  do j=1,nvar
                        qx(j,i) = 0.0d0
                        qy(j,i) = 0.0d0
                  enddo
            endif
      enddo
c     close(20)

      return
      end
