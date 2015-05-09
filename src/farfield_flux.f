C.....Find flux for farfield boundary edges using Riemann invariants
C.....This is not appropriate for unsteady flows
C.....Far-field vortex correction, see Blazek
C.....Chord length is assumed to be 1, otherwise change circ
C.....(xref,yref) is center of vortex, taken to be quarter chord point
C-----------------------------------------------------------------------------
C Flux for a farfield edge: Steger-Warming flux splitting
C-----------------------------------------------------------------------------
      subroutine farfield_flux(x1, x2, prim1, prim2, res1, res2)
      implicit none
      include 'param.h'
      double precision x1(2), x2(2), prim1(nvar), prim2(nvar),
     +                 res1(nvar), res2(nvar)

      integer          i, j, k
      double precision dx, dy, dref, theta, circ, fact1, fact2, fact3, 
     +                 fact4, fact, uinf, vinf, pinf, rinf, ainf,
     +                 q2inf, flux1(nvar), flux2(nvar),
     +                 ua, va, pa, qa2, aa2, aa, ha, ra,
     +                 una, vna, ct, st, lent,
     +                 m2, t1, t2, t3, t4, t5, l1, l2, l3, l4,
     +                 lp1, lp2, lp3, lp4, lm1, lm2, lm3, lm4,
     +                 S(nvar,nvar), T(nvar,nvar), Tp(nvar,nvar),
     +                 Tm(nvar,nvar), jacp(nvar,nvar),
     +                 jacm(nvar,nvar), Uin1(nvar), Uin2(nvar),
     +                 Uout(nvar)

      if(mach_inf .lt. 1.0d0 .and. vortex .eq. yes)then
         dx    = 0.5d0*( x1(1) + x2(1) ) - xref
         dy    = 0.5d0*( x1(2) + x2(2) ) - yref
         dref  = dsqrt(dx**2 + dy**2)
         theta = datan2(dy, dx)
         circ  = 0.5d0*q_inf*Cl
         fact1 = circ*dsqrt(1.0d0 - mach_inf**2)
         fact2 = 2.0d0*M_PI*dref
         fact3 = 1.0d0 - (mach_inf*dsin(theta - aoa))**2
         fact  = fact1/(fact2*fact3)
         uinf  = u_inf + fact*dsin(theta)
         vinf  = v_inf - fact*dcos(theta)
         q2inf = uinf**2 + vinf**2
         fact4 = 1.0d0 + 0.5d0*GAMMA1*(q_inf**2 - q2inf)/a_inf**2
         rinf  = r_inf*fact4**(1.0d0/GAMMA1)
         pinf  = p_inf*(rinf/r_inf)**GAMMA
         ainf  = dsqrt(GAMMA*pinf/rinf)
      else
         uinf  = u_inf
         vinf  = v_inf
         q2inf = uinf**2 + vinf**2
         pinf  = p_inf
         rinf  = r_inf
         ainf  = a_inf
      endif

c     Conserved state for infinity
      Uout(1) = rinf
      Uout(2) = rinf*uinf
      Uout(3) = rinf*vinf
      Uout(4) = pinf/gamma1 + 0.5d0*rinf*q2inf
      

c Average state on edge
      ra = 0.5d0*(prim1(1) + prim2(1))
      ua = 0.5d0*(prim1(2) + prim2(2))
      va = 0.5d0*(prim1(3) + prim2(3))
      pa = 0.5d0*(prim1(4) + prim2(4))
      qa2= ua**2 + va**2
      aa2= GAMMA*pa/ra
      aa = dsqrt(gamma*pa/ra)
      ha = aa2/gamma1 + 0.5d0*qa2

c     Conserved state for prim1
      Uin1(1) = prim1(1)
      Uin1(2) = prim1(1)*prim1(2)
      Uin1(3) = prim1(1)*prim1(3)
      Uin1(4) = prim1(4)/gamma1 + 0.5d0*prim1(1)*(prim1(2)**2 +
     +                                            prim1(3)**2)

c     Conserved state for prim2
      Uin2(1) = prim2(1)
      Uin2(2) = prim2(1)*prim2(2)
      Uin2(3) = prim2(1)*prim2(3)
      Uin2(4) = prim2(4)/gamma1 + 0.5d0*prim2(1)*(prim2(2)**2 +
     +                                            prim2(3)**2)

      ct   =  (x2(2) - x1(2))
      st   = -(x2(1) - x1(1))
      lent = dsqrt(ct**2 + st**2)
      ct   = ct/lent
      st   = st/lent

C     Rotated velocity
      una = ua*ct + va*st
      vna =-ua*st + va*ct

C     Eigenvalues
      l1 = una
      l2 = una
      l3 = una + aa
      l4 = una - aa

C     Positive Eigenvalues
      lp1 = dmax1(l1, 0.0d0)
      lp2 = lp1
      lp3 = dmax1(l3, 0.0d0)
      lp4 = dmax1(l4, 0.0d0)

C     Negative Eigenvalues
      lm1 = l1 - lp1
      lm2 = l2 - lp2
      lm3 = l3 - lp3
      lm4 = l4 - lp4

c     Right eigenvector matrix
      t1 = 0.5d0*ra/aa
      m2 = qa2/aa2

      T(1,1) = 1.0d0
      T(2,1) = ua
      T(3,1) = va
      T(4,1) = 0.5d0*qa2

      T(1,2) = 0.0d0
      T(2,2) = ra*st
      T(3,2) = -ra*ct
      T(4,2) = -ra*vna

      T(1,3) = t1
      T(2,3) = t1*(ua + aa*ct)
      T(3,3) = t1*(va + aa*st)
      T(4,3) = t1*(ha + aa*una)

      T(1,4) = t1
      T(2,4) = t1*(ua - aa*ct)
      T(3,4) = t1*(va - aa*st)
      T(4,4) = t1*(ha - aa*una)

c     Inverse of right eigenvector matrix
      t1     = 0.5d0*gamma1*m2*aa/ra
      t2     = una/ra
      t3     = gamma1*ua/aa
      t4     = gamma1*va/aa
      t5     = gamma1/(ra*aa)

      S(1,1) = 1.0d0 - 0.5d0*gamma1*m2
      S(2,1) = vna/ra
      S(3,1) = t1 - t2
      S(4,1) = t1 + t2

      S(1,2) = t3/aa
      S(2,2) = st/ra
      S(3,2) = (ct - t3)/ra
      S(4,2) =-(ct + t3)/ra

      S(1,3) = t4/aa
      S(2,3) = -ct/ra
      S(3,3) = (st - t4)/ra
      S(4,3) =-(st + t4)/ra

      S(1,4) = -gamma1/aa2
      S(2,4) = 0.0d0
      S(3,4) = t5
      S(4,4) = t5

c     Multiply T * lambda
      do i=1,nvar
         Tp(i,1) = lp1*T(i,1)
         Tp(i,2) = lp2*T(i,2)
         Tp(i,3) = lp3*T(i,3)
         Tp(i,4) = lp4*T(i,4)
         Tm(i,1) = lm1*T(i,1)
         Tm(i,2) = lm2*T(i,2)
         Tm(i,3) = lm3*T(i,3)
         Tm(i,4) = lm4*T(i,4)
      enddo

c     Now multiply with S to get pos/neg jacobians
      do i=1,nvar
         do j=1,nvar
            jacp(i,j) = 0.0d0
            jacm(i,j) = 0.0d0
            do k=1,nvar
               jacp(i,j) = jacp(i,j) + Tp(i,k)*S(k,j)
               jacm(i,j) = jacm(i,j) + Tm(i,k)*S(k,j)
            enddo
         enddo
      enddo

c     Finally the flux jacp*Uin + jacm*Uout
      do i=1,nvar
         flux1(i) = 0.0d0
         flux2(i) = 0.0d0
         do j=1,nvar
            flux1(i) = flux1(i) + jacp(i,j)*Uin1(j) + jacm(i,j)*Uout(j)
            flux2(i) = flux2(i) + jacp(i,j)*Uin2(j) + jacm(i,j)*Uout(j)
         enddo
         flux1(i) = 0.5d0*lent*flux1(i)
         flux2(i) = 0.5d0*lent*flux2(i)
         res1(i)  = res1(i) + flux1(i)
         res2(i)  = res2(i) + flux2(i)
      enddo

      return
      end
