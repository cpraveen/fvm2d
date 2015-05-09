C.....Nearest distance to wall for each point, required by some
C.....turbulence models
      subroutine wall_dist(edge, ptype, coord, spts, psup1, psup2, wd)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), ptype(npmax), spts(nspmax),
     &                 psup1(mpsup*npmax), psup2(npmax+1)
      double precision coord(2,npmax), wd(npmax)

      integer          i, j, je, n1, n2, sedge(nbemax)
      double precision dx1, dy1, dx2, dy2, dx3, dy3, ds1, ds2, dl2, xi, 
     &                 xp, yp, dsperp

      print*,'Finding smallest distance to wall...'

c Find solid wall edges
      nbe = 0
      do i=1,ne
            n1 = edge(1,i)
            n2 = edge(2,i)
            if(ptype(n1) .eq. solid .and. ptype(n2) .eq. solid)then
                  nbe        = nbe + 1
                  sedge(nbe) = i
            endif
      enddo

c Find nearest distance of each point to the solid edges
      do i=1,np
            wd(i) = 1.0d10
            do j=1,nbe
                  je    = sedge(j)
                  n1    = edge(1,je)
                  n2    = edge(2,je)

                  dx1   = coord(1,n1) - coord(1,i)
                  dy1   = coord(2,n1) - coord(2,i)
                  ds1   = dsqrt(dx1**2 + dy1**2)
                  wd(i) = dmin1(wd(i), ds1)

                  dx2   = coord(1,n2) - coord(1,i)
                  dy2   = coord(2,n2) - coord(2,i)
                  ds2   = dsqrt(dx2**2 + dy2**2)
                  wd(i) = dmin1(wd(i), ds2)

                  dx3   = coord(1,n2) - coord(1,n1)
                  dy3   = coord(1,n2) - coord(1,n1)

                  dl2   = dx3**2 + dy3**2

                  xi    = -(dx1*dx3 + dy1*dy3)/dl2
                  if(xi .ge. 0.0d0 .and. xi .le. 1.0d0)then
                        xp     = coord(1,n1) + xi*dx3
                        yp     = coord(2,n1) + xi*dy3
                        dsperp = dsqrt( (coord(1,i)-xp)**2 +
     &                                  (coord(2,i)-yp)**2 )
                        wd(i)  = dmin1(wd(i), dsperp)
                  endif
            enddo
            
            if(ptype(i) .eq. solid) wd(i) = 0.0d0
      enddo

c Distance of second layer of points to wall
      do i=1,nsp
            n1     = spts(i)
            wd1(i) = 1.0d10
            do j=psup2(n1)+1, psup2(n1+1)
                  n2 = psup1(j)
                  if(ptype(n2) .ne. solid)then
                        wd1(i) = dmin1(wd1(i), wd(n2))
                  endif
            enddo
      enddo

      wdmin = 1.0d10
      wdmax = 0.0d0
      do i=1,nsp
            wdmin  = dmin1(wdmin,  wd1(i))
            wdmax  = dmax1(wdmax,  wd1(i))
      enddo

      print*,'Number of solid faces         =',nbe
      write(*,10)'Min distance of nearest point =',wdmin
      write(*,10)'Max distance of nearest point =',wdmax
10    format(a32,e10.3)

c Set distance of solid points to non-zero value, avoids division by
c zero. Note that this does not really affect anything since turbulent
c viscosity at a solid point is anyway fixed to be zero
      do i=1,np
            if(ptype(i) .eq. solid) wd(i) = 1.0d0
      enddo

      return
      end
