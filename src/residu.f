C Calculates finite volume residual
      subroutine residu(coord, elarea, ds, prim, qx, qy, mu, mul,
     &                  nut, mcarea, res)
      implicit none
      include 'param.h'
      include 'gdata.h'

      double precision coord(2, npmax), elarea(ntmax), ds(2,nemax),
     &                 res(nvar,npmax), prim(nvar,npmax), 
     &                 qx(nvar,npmax), qy(nvar,npmax), mu(npmax), 
     &                 mul(npmax), nut(npmax), mcarea(npmax)

      integer          i, ip, iv, ie, n1, n2, e1, e2
      double precision priml(nvar), primr(nvar)

      do ip=1,np
         do iv=1,nvar
            res(iv,ip) = 0.0d0
         enddo
      enddo

c     Sutherland viscosity
      if(iflow .ne. inviscid)then
         call sutherland(prim, mul)
         call viscosity(prim, nut, mul, mu)
      endif

c     Find gradients for reconstruction
      if(iflow .eq. inviscid)then
         call green_gauss_euler(coord, mcarea, prim, qx, qy)
      else
         call green_gauss_visc(coord, elarea, mcarea, prim, nut,
     &                          res, mul, mu, qx, qy)
      endif

c     Inviscid flux for interior edges
      if(iflux .eq. iroe)then
         do i=1,ne
            ie = iedge(i)
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            call recon(coord(1,e1), coord(1,e2), prim(1,e1), prim(1,e2),
     &                 qx(1,e1), qy(1,e1), qx(1,e2), qy(1,e2),
     &                 priml, primr)
            call roe_flux(ds(1,ie), priml, primr, res(1,e1), res(1,e2))
         enddo
      elseif(iflux .eq. ikfvs)then
         do i=1,ne
            ie = iedge(i)
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            call recon(coord(1,e1), coord(1,e2), prim(1,e1), prim(1,e2),
     &                 qx(1,e1), qy(1,e1), qx(1,e2), qy(1,e2),
     &                 priml, primr)
            call kfvs_flux(ds(1,ie), priml, primr, res(1,e1), res(1,e2))
         enddo
      elseif(iflux .eq. ihllc)then
         do i=1,ne
            ie = iedge(i)
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            call recon(coord(1,e1), coord(1,e2), prim(1,e1), prim(1,e2),
     &                 qx(1,e1), qy(1,e1), qx(1,e2), qy(1,e2),
     &                 priml, primr)
            call hllc_flux(ds(1,ie), priml, primr, res(1,e1), res(1,e2))
         enddo
      else
         print*,'FATAL ERROR: Flux type is not defined, iflux=',iflux
         stop
      endif

C     Flux for solid boundary edges
      do i=nswe1,nswe2
         ie = beindx(i)
         n1 = edge(1,ie)
         n2 = edge(2,ie)
         call solid_flux(coord(1,n1), coord(1,n2), prim(1,n1),
     +                   prim(1,n2), res(1,n1), res(1,n2))
      enddo

C     Flux for far field points - far field bc
      do i=nffe1,nffe2
         ie = beindx(i)
         n1 = edge(1,ie)
         n2 = edge(2,ie)
         call farfield_flux(coord(1,n1), coord(1,n2), prim(1,n1),
     +                      prim(1,n2), res(1,n1), res(1,n2))
      enddo

      return
      end
