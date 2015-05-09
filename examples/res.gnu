       reset
       set term png
       set out 'res.png'

       set multiplot

       set size 0.5,0.5
       set origin 0.0,0.5
       set xran[-0.05:1.05]
       set nokey
#      set bmargin 1
       set title "Pressure/Friction coeffient"
       set y2tics
       set ylabel '-Cp'
       set y2label 'Cf'
#      plot 'WALL.DAT' w lp lw 2 pt 6,'WALL.DAT' u 1:3 w lp axes x1y2
       plot 'WALL.DAT' w lp lw 2 pt 6

       set size 0.5,0.5
       set origin 0.5,0.5
       set auto
       set logscale y
       set title "Residue"
       plot 'RES.DAT' u 1:2 w l
       set nologscale y

       set size 0.5,0.5
       set origin 0.07,0.0
       set xran[-0.5:1.5]
       set yran[-1.0:1.0]
       set size ratio -1
       set noxtics
       set noytics
       set nokey
       set title "Mach number"
#      set bmargin 1
       plot 'GNU.MACH' w l,'BD.DAT' w l lt 1 lw 2

       set size 0.5,0.5
       set origin 0.57,0.0
       set xran[-0.5:1.5]
       set yran[-1.0:1.0]
       set size ratio -1
       set nokey
#      set bmargin 1
       set title "Pressure"
       plot 'GNU.PRES' w l,'BD.DAT' w l lt 1 lw 2

       unset multiplot
#      set term x11
