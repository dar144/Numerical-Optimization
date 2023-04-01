#!/usr/bin/gnuplot 

set term png enhanced size 600,300 

set size ratio -1

set o "psi-1000.png"
set contour
#set cntrparam levels 40 # lub ponizsze:
set cbr [0:40]
set cntrparam levels increment -200,10,350
unset surface
set view map
unset key
sp "7_100.txt" u 1:2:3:3 w l lt -1 palette  t '' 