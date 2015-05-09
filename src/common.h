c nvar = number of variables, fixed at 5, which includes turbulent
c viscosity
      integer nvar
      parameter(nvar=4)

C     GAMMA = ratio of specific heats
C     GAMMA1= GAMMA - 1
C     GAS_CONST = gas constant, this can be set to 1.0
C     M_PI = value of pi
      double precision GAMMA, GAMMA1, GAS_CONST, M_PI
      common/const/GAMMA, GAMMA1, GAS_CONST, M_PI

      double precision Rey, prandtl, prandtl_turb, SCONST
      common/viscparam/Rey, prandtl, prandtl_turb, SCONST

      double precision Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1, 
     &                 Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2
      common/samodel/Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1,
     &               Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2
