      subroutine vigie(coord, elem, prim, nut)
      implicit none
      include 'param.h'
      integer           elem(nvemax,ntmax)
      double precision  coord(2,npmax), prim(nvar,npmax), nut(npmax)

      integer           i, vig
      double precision  q2, mach

      vig = 55
      open(unit=vig, file='RESULT.VIG')

      write(vig,*)'points',np
      do i=1,np
            write(vig,*) coord(1,i), coord(2,i)
      enddo

      write(vig,*)'triangles',nt
      do i=1,nt
            write(vig,*) elem(1,i), elem(2,i), elem(3,i)
      enddo

      write(vig,*)'scalars  Mach'
      do i=1,np
            q2   = prim(2,i)**2 + prim(3,i)**2
            mach = dsqrt(q2*prim(1,i)/(GAMMA*prim(4,i)))
            write(vig,*) mach
      enddo

      write(vig,*)'scalars  Pressure'
      do i=1,np
            write(vig,*) prim(4,i)
      enddo

      write(vig,*)'scalars  Density'
      do i=1,np
            write(vig,*) prim(1,i)
      enddo

      write(vig,*)'scalars  Entropy'
      do i=1,np
            write(vig,*) dlog10(prim(4,i)/prim(1,i)**GAMMA/ent_inf )
      enddo

      if(iflow .eq. turbulent)then
      write(vig,*)'scalars  Viscosity'
            do i=1,np
                  write(vig,*)  nut(i)/nmax
            enddo
      endif


      write(vig,*)'vectors  vel  u  v  1e-02'
      do i=1,np
            write(vig,*) prim(2,i), prim(3,i)
      enddo
      write(vig,*)'end_block'

      close(vig)

      return
      end
