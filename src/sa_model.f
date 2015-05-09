C.....Spalart-Allmaras turbulence model
C.....Implicit treatment of destruction term
      subroutine sa_model(coord, elarea, cvarea, prim, prim_old, wd, 
     &                    nut, mul, ds, dt, qx, qy)
      implicit none
      include 'param.h'
      double precision coord(2,npmax), elarea(ntmax), cvarea(npmax),
     &                 prim(nvar,npmax), nut(npmax),
     &                 wd(npmax), mul(npmax), ds(2,nemax),
     &                 prim_old(nvar,npmax), dt(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax)

      double precision divf(npmax)


      call viscous_sa(coord, elarea, prim, nut, divf, mul)
      
      call divergence_sa(coord, ds, prim, nut, divf)

      call update_sa(divf, prim_old, prim, cvarea, dt, nut, mul, wd, 
     &               qx, qy)


      return
      end

c Viscous for SA model
      subroutine viscous_sa(coord, elarea, prim, nut, divf, mul)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), prim(nvar,npmax), nut(npmax),
     &                 divf(npmax), mul(npmax), elarea(ntmax)

      integer          is, jt, n1, n2, n3
      double precision dbxx(3), dbyy(3), nux, nuy, areat, prim5t,
     &                 divsa1, divsa2, divsa3, mult, mutb, rhot

      do is=1,np
            divf(is) = 0.0d0
      enddo

      do jt=1,nt
         areat    = 0.25d0/elarea(jt)

         n1       = elem(1,jt)
         n2       = elem(2,jt)
         n3       = elem(3,jt)

         dbxx(1)  = coord(2,n2) - coord(2,n3)
         dbxx(2)  = coord(2,n3) - coord(2,n1)
         dbxx(3)  = coord(2,n1) - coord(2,n2)
         dbyy(1)  = coord(1,n3) - coord(1,n2)
         dbyy(2)  = coord(1,n1) - coord(1,n3)
         dbyy(3)  = coord(1,n2) - coord(1,n1)

         nux      = nut(1)*dbxx(1) + 
     &              nut(2)*dbxx(2) +
     &              nut(3)*dbxx(3)
         nuy      = nut(1)*dbyy(1) +
     &              nut(2)*dbyy(2) +
     &              nut(3)*dbyy(3)

         mult     = (mul(n1)    + mul(n2)    + mul(n3))/3.0d0
         prim5t   = (nut(n1)    + nut(n2)    + nut(n3))/3.0d0
         rhot     = (prim(1,n1) + prim(1,n2) + prim(1,n3))/3.0d0
         mutb     = mult/rhot + prim5t

         divsa1   = (Cb2Sig1*mutb - Cb2Sig2*(mul(n1)/prim(1,n1) +
     &              nut(n1)))*(nux*dbxx(1) + nuy*dbyy(1))*areat

         divsa2   = (Cb2Sig1*mutb - Cb2Sig2*(mul(n2)/prim(1,n2) +
     &              nut(n2)))*(nux*dbxx(2) + nuy*dbyy(2))*areat

         divsa3   = (Cb2Sig1*mutb - Cb2Sig2*(mul(n3)/prim(1,n3) +
     &              nut(n3)))*(nux*dbxx(3) + nuy*dbyy(3))*areat

         divf(n1) = divf(n1) + divsa1
         divf(n2) = divf(n2) + divsa2
         divf(n3) = divf(n3) + divsa3
     
      enddo

      return
      end

c Convective flux for SA model
      subroutine divergence_sa(coord, ds, prim, nut, divf)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), ds(2,nemax),
     &                 prim(nvar,npmax), nut(npmax), divf(npmax)

      integer          ii, i, n1, n2
      double precision unl, unr, una, flux

      do ii=1,ne
            i    = iedge(ii)
            n1   = edge(1,i)
            n2   = edge(2,i)

            unl  = prim(2,n1)*ds(1,i) + prim(3,n1)*ds(2,i)
            unr  = prim(2,n2)*ds(1,i) + prim(3,n2)*ds(2,i)
            una  = 0.5d0*(unl + unr)

            if(una .gt. 0.0d0)then
                  flux = una*nut(n1)
            else
                  flux = una*nut(n2)
            endif

            divf(n1) = divf(n1) + flux
            divf(n2) = divf(n2) - flux
      enddo

      call outer_flux_sa(coord, prim, nut, divf)

      return
      end

C     Flux for outer boundary edges - far field bc
      subroutine outer_flux_sa(coord, prim, nut, divf)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), prim(nvar,npmax), nut(npmax), 
     +                 divf(npmax)

      integer          i, ie, n1, n2
      double precision dx, dy, sx, sy, u, v, un, flux1, flux2

c     Farfield edges
      do i=nffe1,nffe2
         ie = beindx(i)
         n1 = edge(1,ie)
         n2 = edge(2,ie)
         dx = coord(1,n2) - coord(1,n1)
         dy = coord(2,n2) - coord(2,n1)
         sx = dy
         sy =-dx

         u  = 0.5d0*(prim(2,n1) + prim(2,n2))
         v  = 0.5d0*(prim(3,n1) + prim(3,n2))
         un = u*sx + v*sy
         if(un.gt.0.0d0)then
            flux1 = un*nut(n1)
            flux2 = un*nut(n2)
         else
            flux1 = un*nut_inf
            flux2 = un*nut_inf
         endif
         divf(n1) = divf(n1) + 0.5d0*flux1
         divf(n2) = divf(n2) + 0.5d0*flux2
      enddo

      return
      end

C.....Update the turbulent viscosity
      subroutine update_sa(divf, prim_old, prim, cvarea, dt, nut,
     &                     mul, wd, qx, qy)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision divf(npmax), prim_old(nvar,npmax), 
     &                 prim(nvar,npmax), cvarea(npmax), dt(npmax),
     &                 mul(npmax), wd(npmax), qx(nvar,npmax),
     &                 qy(nvar,npmax), nut(npmax)

      integer          i, j
      double precision chi, chi3, fv1, fv2, fv3, r, Omega, S, g, fw, 
     &                 fact, n1b6, source

      n1b6 = 1.0d0/6.0d0
      do i=1,np
            chi    = nut(i)*prim(1,i)/mul(i)
            chi3   = chi**3
            fv1    = chi3/(chi3 + Cv11)
            fv2    = 1.0d0/(1.0d0 + chi/Cv2)**3
            r      = nut(i)/wd(i)/(wd(i)*kolm2)
            Omega  = qy(2,i) - qx(3,i)
            fv3    = (1.0d0 + chi*fv1)*(1.0d0 - fv2)/dmax1(chi, 0.001d0)
            S      = fv3*dabs(Omega) + r*fv2
            r      = r/S
            r      = dmin1(r,2.0d0)
            g      = r + Cw2*(r**6 - r)
            fw     = g*(Cw31/(g**6 + Cw32))**n1b6
            source = Cb1*S*nut(i)
            fact   = 1.0d0 + Cw1*fw*(nut(i)/wd(i))*(dt(i)/wd(i))
            nut(i) = (nut(i) -  (dt(i)/cvarea(i))*divf(i) +
     &                   dt(i)*source)/fact
      enddo
                  
      do i=1,nsp
            j      = spts(i)
            nut(j) = 0.0d0
      enddo

      return
      end
