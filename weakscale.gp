set terminal x11
set key at 3.25,29
#set key width 1
set key horizontal
set title "Weak Scaling"
set xtic auto
set ytic auto
set xlabel "Processors"
#set xr [0.0:30.0]
set ylabel "Weak Scaling"
#set yr [-12.0:12.0]

set style line 2 lt 1 lw 2 pt 4 ps 0.75 linecolor rgb "red"
plot "weak_scaling.txt" u 1:2 w l ls 2



################################################################################

set term png
set output "wsplot.png"
replot
set term x11
