C.....Flux for solid boundary edges
      subroutine solid_flux(x1, x2, prim1, prim2, resl, resr)
      implicit none
      include 'common.h'
      double precision x1(2), x2(2), prim1(nvar), prim2(nvar),
     +                 resl(nvar), resr(nvar)

      double precision sx, sy, p1, p2

      sx      = 0.5d0*(x2(2) - x1(2))
      sy      =-0.5d0*(x2(1) - x1(1))
      p1      = 0.75d0*prim1(4) + 0.25d0*prim2(4)
      p2      = 0.25d0*prim1(4) + 0.75d0*prim2(4)
      resl(2) = resl(2) + p1*sx
      resl(3) = resl(3) + p1*sy
      resr(2) = resr(2) + p2*sx
      resr(3) = resr(3) + p2*sy

      return
      end
