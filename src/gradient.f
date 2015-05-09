      subroutine gradient(coord, ds, cvarea, prim, qx, qy)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,*), ds(2,*), cvarea(*), prim(nvar,*), 
     +                 qx(nvar,*), qy(nvar,*)

      integer          i, j, ie, n1, n2
      double precision dq, qa, qa1,qa2,dx, dy, nx, ny, ai, tmpx, tmpy

      do i=1,np
         do j=1,nvar
            qx(j,i) = 0.0d0
            qy(j,i) = 0.0d0
         enddo
      enddo

c     Loop over interior non-zero edges
      do i=1,ne
         ie = iedge(i)
         n1 = edge(1,ie)
         n2 = edge(2,ie)
         do j=1,nvar
            qa       = 0.5d0*( prim(j,n1) + prim(j,n2) )
c           dq       = 0.5d0*( prim(j,n2) - prim(j,n1) )
            tmpx     = qa*ds(1,ie)
            tmpy     = qa*ds(2,ie)
            qx(j,n1) = qx(j,n1) + tmpx
            qy(j,n1) = qy(j,n1) + tmpy
            qx(j,n2) = qx(j,n2) - tmpx
            qy(j,n2) = qy(j,n2) - tmpy
         enddo
      enddo

c     Loop over all boundary edges
      do i=1,nbe
         ie = beindx(i)
         n1 = edge(1,ie)
         n2 = edge(2,ie)
         dx = coord(1,n2) - coord(1,n1)
         dy = coord(2,n2) - coord(2,n1)
         nx = 0.5d0*dy
         ny =-0.5d0*dx
         do j=1,nvar
            qa1      = ( 5.0d0*prim(j,n1) +       prim(j,n2) )/6.0d0
            qa2      = (       prim(j,n1) + 5.0d0*prim(j,n2) )/6.0d0
c           dq       = 0.5d0*( prim(j,n2) - prim(j,n1) )
            qx(j,n1) = qx(j,n1) + qa1*nx
            qy(j,n1) = qy(j,n1) + qa1*ny
            qx(j,n2) = qx(j,n2) + qa2*nx
            qy(j,n2) = qy(j,n2) + qa2*ny
         enddo
      enddo

c     Divide by area
      do i=1,np
         ai = 1.0d0/cvarea(i)
         do j=1,nvar
            qx(j,i) = ai*qx(j,i)
            qy(j,i) = ai*qy(j,i)
         enddo
         print*,i,qx(1,i),qy(1,i)
      enddo

      return
      end
