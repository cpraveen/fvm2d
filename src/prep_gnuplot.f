      subroutine prep_gnuplot
      implicit none
      include 'param.h'
      integer gnu

      gnu = 10
      open(unit=gnu, file='grid.gnu')
      write(gnu,*)"set xrange[",xmin,":",xmax,"]"
      write(gnu,*)"set yrange[",ymin,":",ymax,"]"
      write(gnu,*)"set size ratio -1"
      write(gnu,*)"set nokey"
      write(gnu,*)"p \'BD.DAT\' w l"
      write(gnu,*)"pause 5"
      write(gnu,*)"p \'GRID.DAT\' w l"
      write(gnu,*)"pause 5"
      write(gnu,*)"p \'DUAL.DAT\' w l,\'BD.DAT\' w l"
      write(gnu,*)"pause 5"
      close(gnu)

      return
      end
