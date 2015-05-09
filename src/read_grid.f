C----------------------------------------------------------------------
C.....Read grid data from a file
C.....Currently supports only triangular elements
C----------------------------------------------------------------------
      subroutine read_grid(coord, elem, ptype, spts, fpts, opts, bpts)
      implicit none
      include 'param.h'
      double precision coord(2, npmax)
      integer          elem(nvemax, ntmax), ptype(npmax), spts(nspmax),
     &                 fpts(nfpmax), opts(nopmax), bpts(nbpmax)

      integer          ngrid, ip, i, j

      print*,'Reading grid from file ',gridfile

      ngrid = 10
      open(ngrid, file=gridfile, status="old")
      rewind(ngrid)
      read(ngrid,*) np, nt
      print*, 'Number of points    : ', np
      print*, 'Number of triangles : ', nt

      if(np.gt.npmax) then
            print*, 'Increase the size of npmax'
            stop
      endif

      if(nt.gt.ntmax) then
            print*, 'Increase the size of ntmax'
            stop
      endif

      do ip=1,np
            read(ngrid,*) i, coord(1,ip), coord(2,ip), ptype(ip)
      enddo

      do ip=1,nt
            read(ngrid,*) i, elem(1, ip), elem(2, ip), elem(3, ip), j
      enddo

      close(ngrid)

c     Find bounding box
      xmin = 1000000.0d0
      ymin = 1000000.0d0
      xmax =-1000000.0d0
      ymax =-1000000.0d0
      do ip=1,np
            xmin = dmin1(xmin, coord(1,ip))
            ymin = dmin1(ymin, coord(2,ip))

            xmax = dmax1(xmax, coord(1,ip))
            ymax = dmax1(ymax, coord(2,ip))
      enddo

      print*,'Bounding box:'
      print*, '   xmin = ', xmin
      print*, '   xmax = ', xmax
      print*, '   ymin = ', ymin
      print*, '   ymax = ', ymax

      nsp = 0
      nfp = 0
      nop = 0
      nbp = 0

      do ip=1,np
            if(ptype(ip) .eq. solid)then
                  nsp       = nsp + 1
                  spts(nsp) = ip
            endif
            if(ptype(ip) .eq. farfield)then
                  nfp       = nfp + 1
                  fpts(nfp) = ip
            endif
            if(ptype(ip) .eq. outflow)then
                  nop       = nop + 1
                  opts(nop) = ip
            endif
            if(ptype(ip) .ne. interior)then
                  nbp       = nbp + 1
                  bpts(nbp) = ip
            endif
      enddo

      print*,'Number of solid    points = ',nsp
      print*,'Number of farfield points = ',nfp
      print*,'Number of outflow  points = ',nop
      print*,'Number of boundary points = ',nbp

      if( nsp+nfp+nop .ne. nbp )then
            print*,'There seem to be some unrecognized point types'
            stop
      endif


      return
      end
