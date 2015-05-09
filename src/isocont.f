C.....Prints data into a file for visualizing contours using gnuplot
C.....Taken from NSC2KE of Bijan Mohammadi
      SUBROUTINE isocont(ifile,F,COOR,NVAL,VAL)
      implicit none
      integer          ifile, nval
      double precision F(3),COOR(2,3),VAL(100)

      integer          IP1(3), ival, itr, k
      double precision epsi, ff1, ff2, ff3, ffma, d12, d23, val1, fk, 
     &                 fk1, fmi, fma, dif, eps, hh, x, y
C
      epsi   = 1.0d-5
      IP1(1) = 2
      IP1(2) = 3
      IP1(3) = 1
      FF1    = F(1)
      FF2    = F(2)
      FF3    = F(3)
      FFMA   = DMAX1(DABS(FF1),DABS(FF2))
      FFMA   = DMAX1(ffma,DABS(FF3))
      D12    = DABS(FF1-FF2)
      D23    = DABS(FF2-FF3)
      IF(D12+D23.LT.DMAX1(epsi,epsi*FFMA)) GOTO 1000
C  PAS DE RESTRICTION
C  ******************    
      DO 100 IVAL=1,NVAL
            VAL1 = VAL(IVAL)
            ITR  = 0
            DO 110 K=1,3
                  FK  = F(K)
                  FK1 = F(IP1(K))
                  FMI = DMIN1(FK,FK1)
                  FMA = DMAX1(FK,FK1)
                  DIF = FMA-FMI
                  IF(DIF.LT.epsi) GOTO 110
                  EPS = epsi*DIF
                  IF(VAL1.LT.FMI-EPS .OR. VAL1.GT.FMA+EPS) GOTO 110
                  HH  = DABS(FK-VAL1)/DIF
                  X   = COOR(1,K) + HH*(COOR(1,IP1(K))-COOR(1,K))
                  Y   = COOR(2,K) + HH*(COOR(2,IP1(K))-COOR(2,K))
                  IF(ITR.EQ.0) GOTO 115
                  write(ifile,*) x,y
                  write(ifile,*) 
                  GOTO 100
115               ITR = 1
                  write(ifile,*) x,y
110         CONTINUE
100   CONTINUE
1000  return      
      END
