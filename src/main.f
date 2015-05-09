C Cell-vertex finite volume code on triangles
C Performs inviscid, laminar and turbulent computations
C Spalart Allmaras turbulence model
C This is the main program for Finite Volume Solver
      program fvm
      implicit none
      include 'param.h'
      include 'gdata.h'

      double precision coord(2, npmax), elarea(ntmax), 
     &                 cvarea(npmax), ds(2,nemax), dsb(2,npmax), 
     &                 dt(npmax), res(nvar,npmax), prim(nvar,npmax), 
     &                 prim_old(nvar,npmax), qx(nvar,npmax), 
     &                 qy(nvar,npmax), drmin(npmax), mu(npmax),
     &                 nut(npmax),
     &                 mul(npmax), mcarea(npmax), wd(npmax)

      integer          ifres, irk

      call banner
      call math
      call read_input
      call read_grid(coord, elem, ptype, spts, fpts, opts, bpts)
      call prep_gnuplot
      call geometric(coord, elarea, cvarea, ds, dsb, mcarea, wd,
     &               drmin)
      call initialize(prim, nut, mul, mu)

c     call gradient(coord, ds, cvarea, prim, qx, qy)
c     stop

      iter = iterlast
      fres = MINRESIDUE + 1.0d0
      fres1= 0.0d0

      print*,'Beginning of iterations ...'

      ifres = 11
      call system("rm -f RES.DAT")
      open(unit=ifres, file='RES.DAT')
      rewind(ifres)
      do while(fres .gt. MINRESIDUE .and. iter .lt. iterlast+MAXITER)

         call save_old(prim, prim_old)
         call time_step(prim, mu, drmin, dt)

         do irk=1,NIRK
            call residu(coord, elarea, ds, prim, qx, qy, mu, mul,
     &                  nut, mcarea, res)
            call update(irk, res, prim_old, prim, cvarea, dt)
         enddo

         if(iflow .eq. turbulent)then
            call sa_model(coord, elarea, cvarea, prim, prim_old, wd, 
     &                    mul, ds, dt, qx, qy)
         endif

         call check_positivity(prim, nut, coord)

         call clcd(coord, prim, mu, qx, qy)

         call conres(prim, prim_old, dt)
         write(ifres,'(I8,E14.6,I8,3E14.6)')iter,fres,iresi,fresi,cl,cd

         iter = iter + 1

         if(mod(iter,saveinterval) .eq. 0)then
            call flush(ifres)
            call write_result(prim, nut, mu, qx, qy, coord, dsb)
            if(iter .ge. MAXITER)then
               print*,'*** Stopping execution ***'
               print*,'Maximum number of iterations reached'
               goto 100
            endif
         endif
      enddo
      close(ifres)

      call write_result(prim, nut, mu, qx, qy, coord, dsb)

      print*,'*** Stopping execution ***'
      print*,'Residue has been reduced below MINRES =',MINRESIDUE

100   call finalize(prim, nut, coord)

      stop
      end
