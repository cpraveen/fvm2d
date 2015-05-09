      subroutine write_grid(coord, edge, edneigh)
      implicit none
      include 'param.h'
      double precision coord(2,npmax)
      integer          edge(2,nemax), edneigh(2,nemax)

      integer          gfile, i, n1, n2

c     Write boundary edges to bd.dat
      open(unit=10, file='BD.DAT')
      do i=1,ne
      n1 = edge(1,i)
      n2 = edge(2,i)
      if( edneigh(1,i)*edneigh(2,i) .eq. 0)then
            write(10,*)coord(1,n1), coord(2,n1)
            write(10,*)coord(1,n2), coord(2,n2)
            write(10,*)
      endif
      enddo
      close(10)

c     Write grid into file grid.dat
      gfile=15
      open(unit=gfile, file='GRID.DAT')

      do i=1,ne
            n1 = edge(1,i)
            n2 = edge(2,i)
            write(gfile,*) coord(1,n1), coord(2,n1)
            write(gfile,*) coord(1,n2), coord(2,n2)
            write(gfile,*)
      enddo

      close(gfile)

c     call system('gnuplot -noraise grid.gnu &')

      return
      end
