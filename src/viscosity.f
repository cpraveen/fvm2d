C.....Calculate total Reynolds number
      subroutine viscosity(prim, nut, mul, mu)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), nut(npmax), mul(npmax), 
     +                 mu(npmax)

      integer          i
      double precision mut, mu_turb

      if(iflow .eq. laminar)then
         do i=1,np
            mu(i) = mul(i)
         enddo
      elseif(iflow .eq. turbulent)then
         do i=1,np
            mut   = mu_turb(mul(i), prim(1,i), nut(i))
            mu(i) = mul(i) + mut
         enddo
      endif

      return
      end

C.....Calculate Reynolds number based on Sutherland law
      subroutine sutherland(prim, mul)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), mul(npmax)

      integer          i
      double precision sutherland_viscosity

      do i=1,np
         mul(i)= sutherland_viscosity(prim(1,i), prim(4,i))
      enddo

      return
      end

C Sutherland viscosity formula
      double precision function sutherland_viscosity(density, pressure)
      implicit none
      include 'param.h'
      double precision density, pressure

      double precision temp, num, den

      temp = pressure/density/GAS_CONST
      num  = T_inf + SCONST
      den  = temp  + SCONST
      sutherland_viscosity = (temp/T_inf)**1.5d0*(num/den)/Rey

      return
      end

C Returns turbulent viscosity coefficient
      double precision function mu_turb(mu_lam, density, nu_tilde)
      implicit none
      include 'param.h'
      double precision mu_lam, density, nu_tilde

      double precision chi, chi3, fv1

      chi     = nu_tilde*density/mu_lam
      chi3    = chi**3
      fv1     = chi3/(chi3 + Cv11)
      mu_turb = density*nu_tilde*fv1

      return
      end
