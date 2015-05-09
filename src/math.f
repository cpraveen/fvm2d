C----------------------------------------------------------------------
C.....Definition of some constants
C----------------------------------------------------------------------
      subroutine math
      implicit none
      include 'param.h'

      GAMMA        = 1.4d0
      GAMMA1       = GAMMA-1.0d0
      GAS_CONST    = 1.0d0
      M_PI         = 4.0d0*datan2(1.0d0, 1.0d0)

      prandtl      = 0.72d0
      prandtl_turb = 0.9d0

      ALBADA11     = 2.0d0/3.0d0
      ALBADA12     = 1.0d0 - ALBADA11

      ALBADA21     = 4.0d0/3.0d0
      ALBADA22     = 1.0d0 - ALBADA21

      xvstatus     = no

c Constants for Spalart-Allmaras model
      Cb1          = 0.1355d0
      Cb2          = 0.622d0
      sigma_sa     = 2.0d0/3.0d0
      kolm         = 0.41d0
      Cw1          = Cb1/kolm**2 + (1.0d0 + Cb2)/sigma_sa
      Cw2          = 0.3d0
      Cw3          = 2.0d0
      Cv1          = 7.1d0
      Cv2          = 5.0d0

      Cv11         = Cv1**3
      Cw31         = 1.0d0 + Cw3**6
      Cw32         = Cw3**6
      kolm2        = kolm**2
      Cb2Sig1      = (1.0d0 + Cb2)/sigma_sa
      Cb2Sig2      = Cb2/sigma_sa

      return
      end

C----------------------------------------------------------------------
C.....Rotate a vector (u,v) by angle ang
C----------------------------------------------------------------------
      subroutine rotate(u, v, ang)
      implicit none
      double precision u, v, ang

      double precision ut, vt, ct, st

      ut = u
      vt = v

      ct = dcos(ang)
      st = dsin(ang)

      u  = ut*ct + vt*st
      v  =-ut*st + vt*ct

      return
      end

C----------------------------------------------------------------------
C.....Error function, from Abromovitz and Stegun
C----------------------------------------------------------------------
      double precision function ERRF(X)
      double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

      ARG = X*X
      if(ARG .lt. 20.0d0)then
            E = dexp(-ARG)
      else
            E = 0.0d0
      endif
      VB = dabs(X)
      T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
      tmp1 = 1.061405429d0*T
      tmp2 = (tmp1 - 1.453152027d0)*T
      tmp3 = (tmp2 + 1.421413741d0)*T
      tmp4 = (tmp3 - 0.284496736d0)*T
      tmp5 = (tmp4 + 0.254829592d0)*T
      tmp6 = 1.0d0 - tmp5*E
      if(X .lt. 0.0d0)then
            ERRF = -tmp6
      else
            ERRF =  tmp6
      endif

      return
      end

      REAL*8 FUNCTION DRAND(iseed)
c     -----------------------------------------------------------------
c     Selection aleatoire d'un nombre entre 0 et 1 suivant une
c     valeur donnee iseed
c
c                  mod(iseed*7141 + 54773, 259200)
c         ran = -----------------------------------
c                            259200
c     -----------------------------------------------------------------
c
c     Parametres d'appel 
c
      INTEGER iseed
c
c     Variables locales 
c
      INTEGER ia, ic, im
      PARAMETER(ia = 7141, ic = 54773, im = 259200)
c
      iseed    = ABS(MOD(iseed*ia+ic, im))
c
      DRAND    = FLOAT(iseed)/FLOAT(im)
c
      RETURN
      END
