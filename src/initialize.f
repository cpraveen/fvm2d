C----------------------------------------------------------------------
C.....Variables stored are primitive - density, u, v, pressure
C.....Initialize primitive variables to free stream values
C----------------------------------------------------------------------
      subroutine initialize(prim, nut, mul, mu)
      implicit none
      include 'param.h'
      double precision prim(nvar, npmax), nut(npmax), mul(npmax), 
     +                 mu(npmax)
      
      integer          i, j
      double precision u1, u2, u3, u4, u5

      q_inf   = 1.0d0
      r_inf   = 1.0d0
      p_inf   = 1.0d0/(GAMMA*mach_inf**2)
      T_inf   = p_inf/(r_inf*GAS_CONST)
      aoa     = aoa_deg*M_PI/180.0d0
      u_inf   = q_inf*dcos(aoa)
      v_inf   = q_inf*dsin(aoa)
      ent_inf = p_inf/r_inf**GAMMA
      a_inf   = dsqrt(GAMMA*p_inf/r_inf)
      H_inf   = a_inf**2/GAMMA1 + 0.5d0*q_inf

c Required by Sutherland law
      T_infd  = 300.0d0
      SCONST  = 110.4d0*T_inf/T_infd

c Primitive variables in free-stream
      prim_inf(1) = r_inf
      prim_inf(2) = u_inf
      prim_inf(3) = v_inf
      prim_inf(4) = p_inf
      nut_inf     = 0.0d0

      if(iflow.eq.inviscid) print*,'Euler computation'
      if(iflow.eq.laminar)  print*,'Laminar Navier-Stokes computation'
      if(iflow.eq.turbulent)print*,'Turbulent Navier-Stokes computation'
      print*,'Free-stream values:'
      print*,' Mach number =',mach_inf
      print*,' AOA         =',aoa_deg
      print*,' u velocity  =',u_inf
      print*,' v velocity  =',v_inf
      print*,' Pressure    =',p_inf

      if(vortex .eq. yes)then
            print*,'Using point vortex correction for far-field points'
            print*,'Vortex center = ',xref, yref
      endif

C Runge-Kutta time stepping
      NIRK    = 3
      airk(1) = 0.0d0
      airk(2) = 3.0d0/4.0d0
      airk(3) = 1.0d0/3.0d0
      birk(1) = 1.0d0
      birk(2) = 1.0d0/4.0d0
      birk(3) = 2.0d0/3.0d0

      if(istart .eq. scratch)then
            print*,'Initializing solution to free stream values'
            do j=1,np
                  do i=1,nvar
                     prim(i,j) = prim_inf(i)
                  enddo
            enddo

            if(iflow .eq. turbulent)then
                  call sutherland(prim, mul)
                  do i=1,np
                     nut(i) = 0.1d0*mul(i)/prim(1,i)
                  enddo
                  call viscosity(prim, nut, mul, mu)
                  nut_inf = nut(1)
            elseif(iflow .eq. laminar)then
                  call sutherland(prim, mul)
                  do i=1,np
                     nut(i) = 0.0d0
                  enddo
                  call viscosity(prim, nut, mul, mu)
                  nut_inf = 0.0d0
            else
                  do i=1,np
                        mul(i)    = 0.0d0
                        mu(i)     = 0.0d0
                        nut(i) = 0.0d0
                  enddo
                  nut_inf = 0.0d0
            endif

      else
            print*,'Initializing solution to old values from SOL'
            open(unit=20, file='SOL', status='old')
            do j=1,np
                  read(20,*) u1, u2, u3, u4, u5
                  prim(1,j) = u1
                  prim(2,j) = u2/u1
                  prim(3,j) = u3/u1
                  prim(4,j) = GAMMA1*( u4 - 0.5d0*(u2**2 + u3**2)/u1 )
                  nut(j)    = u5
            enddo
            close(20)
      endif

      return
      end
