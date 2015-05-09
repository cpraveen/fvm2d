      subroutine finalize(prim, nut, coord)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision prim(nvar,npmax), nut(npmax), coord(2,npmax)
      real             et, etime, times(2)

      integer          i
      double precision u1, u2, u3, u4, u5

      call vigie(coord, elem, prim, nut)
      call mayavi(coord, elem, prim, nut)

      open(unit=20, file='SOL')
      do i=1,np
            u1 = prim(1,i)
            u2 = u1*prim(2,i)
            u3 = u1*prim(3,i)
            u4 = prim(4,i)/GAMMA1 + 0.5d0*prim(1,i)*(prim(2,i)**2 +
     &                                        prim(3,i)**2)
            u5 = nut(i)
            write(20,*)u1, u2, u3, u4, u5
      enddo
      close(20)

      open(unit=20, file='PARAM.DAT')
      write(20,9)('-',i=1,70)
9     format(70a)
      write(20,10)mach_inf,aoa_deg,Rey
10    format(' Mach no  =',f6.3,', AOA  =',f6.3, ', Rey =',e10.4)
      write(20,11)CFL,iflux,ilimit
11    format(' CFL      =',f6.3,', Flux =',i6, ', Lim =',i2)
      write(20,12) gridfile
12    format(' Gridfile = ',a30)
      write(20,9)('-',i=1,70)
      Write(20,'(" Iterations        =",i12)')iter
      write(20,'(" Cl                =",f12.6)')cl
      write(20,'(" Cd                =",f12.6)')cd
      write(20,'(" L2 residue        =",f12.6)')dlog10(fres)
      write(20,'(" Linf residue      =",e16.6)')fresi
      write(20,'(" Linf point        =",f12.6,f12.6,i8,i4)')
     &      coord(1,iresi), coord(2,iresi), iresi, ptype(iresi)
      write(20,'(" Minimum density   =",f12.6)')rmin
      write(20,'(" Maximum density   =",f12.6)')rmax
      write(20,'(" Minimum pressure  =",f12.6)')pmin
      write(20,'(" Maximum pressure  =",f12.6)')pmax
      write(20,'(" Minimum mach no   =",f12.6)')mmin
      write(20,'(" Maximum mach no   =",f12.6)')mmax
      write(20,'(" Minimum x vel     =",f12.6)')umin
      write(20,'(" Maximum x vel     =",f12.6)')umax
      write(20,'(" Minimum y vel     =",f12.6)')vmin
      write(20,'(" Maximum y vel     =",f12.6)')vmax
      write(20,'(" Minimum entropy   =",f12.6)')emin
      write(20,'(" Maximum entropy   =",f12.6)')emax
      write(20,'(" Minimum viscosity =",e16.6)')nmin
      write(20,'(" Maximum viscosity =",e16.6)')nmax

      et = etime(times)
      write(20,'(" Elapsed time      =",f12.6, " min")')et/60.0
      write(20,'(" User    time      =",f12.6, " min")')times(1)/60.0
      write(20,'(" System  time      =",f12.6, " min")')times(2)/60.0

      close(20)

      return
      end
