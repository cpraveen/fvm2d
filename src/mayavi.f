C.....Result in VTK format for MayaVi
      subroutine mayavi(coord, elem, prim, nut)
      implicit none
      include 'param.h'
      integer           elem(nvemax,ntmax)
      double precision  coord(2,npmax), prim(nvar,npmax), nut(npmax)

      integer           i, maya
      double precision  q2, mach

      maya = 55
      open(unit=maya, file='RESULT.VTK')

      write(maya,10)
10    format('# vtk DataFile Version 3.0')

      if(iflow .eq. inviscid)then
      write(maya,111) mach_inf, aoa_deg
111   format('Mach =', f6.3, 2x, ' AOA = ', f6.3)
      else
      write(maya,112) mach_inf, aoa_deg, Rey
112   format('Mach =', f6.3, 2x, ' AOA = ', f6.3, ' Reynolds = ', e10.4)
      endif

      write(maya,12)
12    format('ASCII')
      write(maya,13)
13    format('DATASET UNSTRUCTURED_GRID')

      write(maya,14) np
14    format('POINTS ', i8, 2x, 'float')
      do i=1,np
            write(maya,'(3e18.8)') coord(1,i), coord(2,i), 0.0
      enddo

      write(maya,15) nt, 4*nt
15    format('CELLS ', 2i10)
      do i=1,nt
            write(maya,'(i4,3i10)') 3, elem(1,i)-1, elem(2,i)-1, 
     &                              elem(3,i)-1
      enddo

      write(maya,16) nt
16    format('CELL_TYPES ', i10)
      do i=1,nt
            write(maya,'(i5)') 5
      enddo

      write(maya,17) np
17    format('POINT_DATA ', i10)
      write(maya,18) 'Mach'
18    format('SCALARS ', a10, '   float 1')
      write(maya,19)
19    format('LOOKUP_TABLE default')
      do i=1,np
            q2   = prim(2,i)**2 + prim(3,i)**2
            mach = dsqrt(q2*prim(1,i)/(GAMMA*prim(4,i)))
            write(maya,20) mach
      enddo
20    format(e18.8)

      write(maya,18) 'Pressure'
      write(maya,19)
      do i=1,np
            write(maya,20) prim(4,i)
      enddo

      write(maya,18) 'Density'
      write(maya,19)
      do i=1,np
            write(maya,20) prim(1,i)
      enddo

      write(maya,18) 'Entropy'
      write(maya,19)
      do i=1,np
            write(maya,20) dlog10(prim(4,i)/prim(1,i)**GAMMA/ent_inf)
      enddo

      if(iflow .eq. turbulent)then
      write(maya,18) 'Viscosity'
      write(maya,19)
            do i=1,np
                  write(maya,*)  nut(i)/nmax
            enddo
      endif

      close(maya)

      return
      end
